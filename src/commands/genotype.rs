use crate::{
    constants::*,
    distributions::{weighted_percentile, DepthDistribution, InsertSizeDistribution},
    extract::{count_irrs, LocusId},
    genotype::{
        get_genotype_search_params, GenotypeInterval, GenotypeProblem, GenotypeProblemDefinition,
        HeterozygousGenotypeProblemWrapper, HomozygousGenotypeProblemWrapper, LocusGenotypeResult,
        TandemRepeatGenotype,
    },
    observers::AlignmentPositionTracker,
    sa::SimulatedAnnealingOptions,
    util::load_aln_inputs,
};
use argmin::{
    core::{Executor, State},
    solver::goldensectionsearch::GoldenSectionSearch,
};
use clap::{arg, Args, Parser};
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256PlusPlus;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::fs::File;
use std::path::PathBuf;
use std::{collections::BTreeMap, path::Path};

/// Estimates the allele lengths for each tandem repeat locus
#[derive(Parser)]
pub struct GenotypeCommandArgs {
    /// Stats JSON file
    #[arg(long)]
    pub stats: Option<PathBuf>,

    /// Defs JSON file
    #[arg(long)]
    pub defs: Option<PathBuf>,

    /// Assume loci are homozygous (otherwise default is heterozygous)
    #[arg(long = "hom")]
    pub homozygous: bool,

    /// Genotype optimization options
    #[command(flatten)]
    pub gt_opts: GenotypeOptimizationArgs,

    /// Simulated annealing options
    #[command(flatten)]
    pub sa_opts: SimulatedAnnealingOptions,

    /// Number of bootstraps to perform to estimate the genotype confidence intervals
    #[arg(short = 'n', long, default_value_t = DEFAULT_N_BOOTSTRAPS)]
    pub n_bootstraps: u64,

    /// Read length
    #[arg(long, default_value_t = INITIAL_ESTIMATE_READ_LENGTH)]
    pub read_length: usize,

    /// Output alignment positions
    #[arg(short = 'o', long)]
    pub output_positions: bool,

    /// Random seed
    #[arg(short = 's', long, default_value_t = 0)]
    pub seed: u64,

    /// Enables logging
    #[arg(long, default_value_t = false)]
    pub log: bool,
}

#[derive(Clone, Debug, Args)]
pub struct GenotypeOptimizationArgs {
    /// The maximum number of iterations for golden section search of genotype
    #[arg(long = "gi", default_value_t = DEFAULT_GSS_ITERS_MAX)]
    pub gss_max_iters: u64,
}

pub fn run_genotype_command(
    args: &GenotypeCommandArgs,
    output_prefix: &Path,
) -> anyhow::Result<()> {
    // Find the paths to the input files
    let definitions_path = if let Some(definitions_path) = &args.defs {
        definitions_path
    } else {
        &output_prefix.with_extension(OUT_SUFFIX_DEFS)
    };

    let stats_path = if let Some(stats_path) = &args.stats {
        stats_path
    } else {
        &output_prefix.with_extension(OUT_SUFFIX_STATS)
    };

    // Load inputs
    let (definitions, depth_distr, insert_distr) = load_aln_inputs(
        args.defs.as_ref().unwrap_or(definitions_path),
        args.stats.as_ref().unwrap_or(stats_path),
        args.sa_opts.use_theoretical,
        args.sa_opts.depth_mean,
        args.sa_opts.insert_mean,
        args.sa_opts.insert_sd,
    )?;

    debug!(
        "Depth mean = {}, insert mean = {}, insert sd = {}",
        depth_distr.mean(),
        insert_distr.mean(),
        insert_distr.sd()
    );

    // Initialize the results map
    let mut gt_results: BTreeMap<LocusId, LocusGenotypeResult> = BTreeMap::new();

    // Process each alignment problem definition
    for aln_problem_def in definitions.iter() {
        info!(
            "[{}] Starting genotype estimation for locus",
            aln_problem_def.locus.id
        );

        // Create a new RNG for each locus so that the results are reproducible regardless
        // of position in the input file
        let mut rng = Xoshiro256PlusPlus::seed_from_u64(args.seed);

        if aln_problem_def.reads.is_empty() {
            warn!(
                "[{}] Locus does not have any reads. Skipping...",
                aln_problem_def.locus.id
            );
            continue;
        }

        // Calculate the genotype search bounds
        let motif_len = aln_problem_def.locus.motif.len() as u32;
        let read_len = args.read_length as u32;
        let depth_mean = depth_distr.mean();
        let depth_sd = depth_distr.sd();
        let num_irrs = count_irrs(&aln_problem_def.reads);
        let (genotype_init, genotype_min, genotype_max) =
            get_genotype_search_params(motif_len, read_len, depth_mean, depth_sd, num_irrs);

        let best_aln_tracker = AlignmentPositionTracker::new();

        debug!(
            "[{}] Using genotype initial estimate = {}, min = {}, max = {}",
            aln_problem_def.locus.id, genotype_init, genotype_min, genotype_max
        );

        let (genotype, genotype_conf_int) = bootstrap_genotype(
            args.n_bootstraps,
            &mut rng,
            aln_problem_def,
            &*depth_distr,
            &*insert_distr,
            &args.gt_opts,
            &args.sa_opts,
            Some(&best_aln_tracker),
            genotype_min,
            genotype_max,
            genotype_init,
            args.homozygous,
        );

        info!(
            "[{}] Genotype = {}, confidence interval = {}",
            aln_problem_def.locus.id, genotype, genotype_conf_int
        );

        // Prepare the result
        let (cost, pos) = best_aln_tracker.get_best(genotype.clone().into()).unwrap();

        let pos = if args.output_positions {
            Some(pos)
        } else {
            None
        };

        let locus_gt_result = LocusGenotypeResult {
            reference_region: format!("{}", aln_problem_def.locus),
            repeat_unit: aln_problem_def.locus.motif.clone(),
            genotype,
            genotype_conf_int,
            cost,
            alignment: pos,
        };

        // Add the locus result to the results map
        gt_results.insert(aln_problem_def.locus.id.clone(), locus_gt_result);
    }

    // Output genotypes
    info!("Writing genotype results to file");
    let gt_results_path = output_prefix.with_extension(OUT_SUFFIX_GENOTYPES);
    let file = File::create(gt_results_path)?;
    serde_json::to_writer_pretty(file, &gt_results)?;

    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub fn bootstrap_genotype(
    n_bootstraps: u64,
    rng: &mut Xoshiro256PlusPlus,
    aln_problem_def: &GenotypeProblemDefinition,
    depth_distr: &dyn DepthDistribution,
    insert_distr: &dyn InsertSizeDistribution,
    gt_opts: &GenotypeOptimizationArgs,
    sa_opts: &SimulatedAnnealingOptions,
    best_aln_tracker: Option<&AlignmentPositionTracker>,
    het_genotype_min: u32,
    het_genotype_max: u32,
    het_genotype_init: u32,
    homozygous: bool,
) -> (TandemRepeatGenotype, GenotypeInterval) {
    // Create a new RNG for each bootstrap run
    let mut bootstrap_rngs = vec![];
    for _ in 0..n_bootstraps {
        bootstrap_rngs.push(Xoshiro256PlusPlus::from_rng(&mut *rng).unwrap());
    }

    // Perform bootstrapping to build a genotype estimate distribution
    let (genotypes, costs): (Vec<TandemRepeatGenotype>, Vec<f64>) = bootstrap_rngs
        .into_par_iter()
        .map(|rng| {
            find_min_cost_genotype(
                aln_problem_def,
                depth_distr,
                insert_distr,
                gt_opts.gss_max_iters,
                sa_opts,
                best_aln_tracker,
                rng,
                het_genotype_min,
                het_genotype_max,
                het_genotype_init,
                homozygous,
            )
        })
        .unzip();

    // Calculate weights for the genotype "observations"
    let weights: Vec<f64> = costs.iter().map(|c| 1.0 / c).collect();

    // Calculate the median and confidence interval of the genotype distribution
    let hap1_gts: Vec<u32> = genotypes.iter().map(|gt| gt.0).collect();
    let hap1_gt_median = weighted_percentile(&hap1_gts, &weights, 50.0).unwrap() as u32;
    let hap1_gt_lower_bound =
        weighted_percentile(&hap1_gts, &weights, GENOTYPE_CONF_INT_LOWER_PERCENTILE).unwrap()
            as u32;
    let hap1_gt_upper_bound =
        weighted_percentile(&hap1_gts, &weights, GENOTYPE_CONF_INT_UPPER_PERCENTILE).unwrap()
            as u32;

    let hap2_gts: Vec<u32> = genotypes.iter().map(|gt| gt.1).collect();
    let hap2_gt_median = weighted_percentile(&hap2_gts, &weights, 50.0).unwrap() as u32;
    let hap2_gt_lower_bound =
        weighted_percentile(&hap2_gts, &weights, GENOTYPE_CONF_INT_LOWER_PERCENTILE).unwrap()
            as u32;
    let hap2_gt_upper_bound =
        weighted_percentile(&hap2_gts, &weights, GENOTYPE_CONF_INT_UPPER_PERCENTILE).unwrap()
            as u32;

    let genotype = TandemRepeatGenotype(hap1_gt_median, hap2_gt_median);
    let genotype_conf_int: GenotypeInterval = GenotypeInterval(
        (hap1_gt_lower_bound, hap1_gt_upper_bound),
        (hap2_gt_lower_bound, hap2_gt_upper_bound),
    );

    (genotype, genotype_conf_int)
}

#[allow(clippy::too_many_arguments)]
pub fn find_min_cost_genotype(
    aln_problem_def: &GenotypeProblemDefinition,
    depth_distr: &dyn DepthDistribution,
    insert_distr: &dyn InsertSizeDistribution,
    gss_max_iters: u64,
    sa_opts: &SimulatedAnnealingOptions,
    best_aln_tracker: Option<&AlignmentPositionTracker>,
    rng: Xoshiro256PlusPlus,
    het_genotype_min: u32,
    het_genotype_max: u32,
    het_genotype_init: u32,
    homozygous: bool,
) -> (TandemRepeatGenotype, f64) {
    let gt_problem = GenotypeProblem::new(
        aln_problem_def,
        sa_opts.base_error_rate,
        depth_distr,
        insert_distr,
        sa_opts,
        best_aln_tracker,
        rng,
    );

    if homozygous {
        let hom_genotype_init = ((het_genotype_init as f64) / 2.0).round() as u32;
        let hom_genotype_min = ((het_genotype_min as f64) / 2.0).floor() as u32;
        let hom_genotype_max = ((het_genotype_max as f64) / 2.0).ceil() as u32;

        let hom_gt_problem = HomozygousGenotypeProblemWrapper(&gt_problem);
        debug!(
            "[{}] Searching for homozygous genotype in range [{}, {}] with initial estimate {}",
            aln_problem_def.locus.id, hom_genotype_min, hom_genotype_max, hom_genotype_init
        );
        let hom_gss_solver =
            GoldenSectionSearch::new(hom_genotype_min as f64, hom_genotype_max as f64).unwrap();
        let hom_executor = Executor::new(hom_gt_problem, hom_gss_solver).configure(|state| {
            state
                .param(hom_genotype_init as f64)
                .max_iters(gss_max_iters)
        });
        let hom_result = hom_executor.run().unwrap();
        let hom_genotype = *hom_result.state().get_best_param().unwrap() as u32;
        let hom_cost = hom_result.state().get_best_cost();

        debug!(
            "[{}] Found homozygous genotype: {:?} (cost = {})",
            aln_problem_def.locus.id, hom_genotype, hom_cost
        );

        (TandemRepeatGenotype(hom_genotype, hom_genotype), hom_cost)
    } else {
        // Calculate "normal" genotype
        let motif_len = aln_problem_def.locus.motif.len() as f64;
        let reference_genotype =
            ((aln_problem_def.locus.end - aln_problem_def.locus.start) as f64 / motif_len).round();

        let het_gt_problem =
            HeterozygousGenotypeProblemWrapper(&gt_problem, reference_genotype as u32);
        debug!(
            "[{}] Searching for heteroyzgous genotype in range [{}, {}] with initial estimate {} and reference genotype {}",
            aln_problem_def.locus.id, het_genotype_min, het_genotype_max, het_genotype_init, reference_genotype
        );

        let het_gss_solver =
            GoldenSectionSearch::new(het_genotype_min as f64, het_genotype_max as f64).unwrap();

        let het_executor = Executor::new(het_gt_problem, het_gss_solver).configure(|state| {
            state
                .param(het_genotype_init as f64)
                .max_iters(gss_max_iters)
        });

        let het_result = het_executor.run().unwrap();

        let het_genotype = *het_result.state().get_best_param().unwrap() as u32;

        let het_cost = het_result.state().get_best_cost();
        debug!(
            "[{}] Found heterozygous genotype: {:?} (cost = {})",
            aln_problem_def.locus.id, het_genotype, het_cost
        );

        (
            TandemRepeatGenotype(reference_genotype as u32, het_genotype),
            het_cost,
        )
    }
}
