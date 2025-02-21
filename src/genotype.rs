use std::{
    collections::{HashMap, HashSet},
    sync::{Arc, Mutex},
};

use crate::{
    alignment::AlignmentProblem,
    catalog::TandemRepeatLocus,
    constants::{INITIAL_ESTIMATE_DEPTH_SD_BUFFER, INITIAL_ESTIMATE_MIN_DEPTH},
    distributions::{DepthDistribution, InsertSizeDistribution},
    extract::Read,
    observers::AlignmentPositionTracker,
    positions::{AlignmentPosition, RepeatAlignmentPositionSetGenerator},
    reference::TandemRepeatReference,
    sa::SimulatedAnnealingOptions,
    sequence::Sequence,
};
use anyhow::Result;
use argmin::{
    core::{CostFunction, Error, Executor, State},
    solver::simulatedannealing::{SATempFunc, SimulatedAnnealing},
};
use rand_xoshiro::{rand_core::SeedableRng, Xoshiro256PlusPlus};
use serde::{
    de::Error as _, ser::SerializeStruct, Deserialize, Deserializer, Serialize, Serializer,
};

#[derive(Clone)]
pub struct GenotypeProblem<'a> {
    aln_position_generators: Vec<RepeatAlignmentPositionSetGenerator>,
    pairs: HashSet<(usize, usize)>,
    read_lens: Vec<usize>,
    base_error_rate: f64,
    left_flank_len: usize,
    right_flank_len: usize,
    motif_len: usize,
    depth_distr: &'a dyn DepthDistribution,
    insert_distr: &'a dyn InsertSizeDistribution,
    sa_opts: SimulatedAnnealingOptions,
    best_aln_tracker: Option<&'a AlignmentPositionTracker>,
    rng: Arc<Mutex<Xoshiro256PlusPlus>>,
}

impl<'a> GenotypeProblem<'a> {
    pub fn new(
        aln_problem_def: &'a GenotypeProblemDefinition,
        base_error_rate: f64,
        depth_distr: &'a dyn DepthDistribution,
        insert_distr: &'a dyn InsertSizeDistribution,
        sa_opts: &SimulatedAnnealingOptions,
        best_aln_tracker: Option<&'a AlignmentPositionTracker>,
        rng: Xoshiro256PlusPlus,
    ) -> Self {
        let aln_position_generators = aln_problem_def.create_position_generators();
        let read_lens: Vec<usize> = aln_problem_def.reads.iter().map(|r| r.seq.len()).collect();

        let max_read_len = *read_lens.iter().max().unwrap();
        let insert_mean = insert_distr.mean() as usize + 1;
        let (left_flank_len, right_flank_len) =
            Self::get_flank_lens(&aln_position_generators, max_read_len, insert_mean);

        let motif_len = aln_problem_def.motif_len();

        let mut pairs = HashSet::new();
        for (i, j) in aln_problem_def.pairs_map.iter() {
            let mut idxs = [*i, *j];
            idxs.sort();
            pairs.insert((idxs[0], idxs[1]));
        }

        Self {
            aln_position_generators,
            pairs,
            read_lens,
            base_error_rate,
            left_flank_len,
            right_flank_len,
            motif_len,
            depth_distr,
            insert_distr,
            sa_opts: sa_opts.clone(),
            best_aln_tracker,
            rng: Arc::new(Mutex::new(rng)),
        }
    }

    fn get_flank_lens(
        aln_position_generators: &[RepeatAlignmentPositionSetGenerator],
        max_read_len: usize,
        insert_mean: usize,
    ) -> (usize, usize) {
        let left_flank_len = aln_position_generators
            .iter()
            .map(|g| g.distance_furthest_from_repeat_start())
            .max()
            .unwrap() as usize;
        let right_flank_len = aln_position_generators
            .iter()
            .map(|g| g.distance_furthest_from_repeat_end())
            .max()
            .unwrap() as usize;
        (
            left_flank_len.max(max_read_len).max(insert_mean),
            right_flank_len.max(max_read_len).max(insert_mean),
        )
    }

    pub fn build_aln_problem(&self, genotype: TandemRepeatGenotype) -> Result<AlignmentProblem> {
        let (hap1_genotype, hap2_genotype) = genotype.clone().into();

        let hap1_len =
            self.left_flank_len + hap1_genotype as usize * self.motif_len + self.right_flank_len;
        let hap2_len =
            self.left_flank_len + hap2_genotype as usize * self.motif_len + self.right_flank_len;

        let insert_mean = self.insert_distr.mean() as usize;

        let hap1_depth_range = if hap1_len <= 2 * insert_mean {
            0..hap1_len
        } else {
            insert_mean..(hap1_len - insert_mean)
        };

        let hap2_depth_range = if hap2_len <= 2 * insert_mean {
            0..hap2_len
        } else {
            insert_mean..(hap2_len - insert_mean)
        };

        let mut rng = self.rng.lock().unwrap();
        let aln_rng = Xoshiro256PlusPlus::from_rng(&mut *rng)?;

        let position_dist_sets = self
            .aln_position_generators
            .iter()
            .map(|g| {
                g.generate_position_set(
                    genotype.clone(),
                    self.left_flank_len,
                    self.right_flank_len,
                    false,
                )
            })
            .collect();

        let mut aln_problem = AlignmentProblem::new(
            &self.read_lens,
            position_dist_sets,
            self.pairs.clone(),
            self.base_error_rate,
            self.depth_distr,
            self.insert_distr,
            hap1_len,
            hap2_len,
            hap1_depth_range,
            hap2_depth_range,
            self.sa_opts.insert_sd_tolerance,
            self.sa_opts.fixed_read_threshold,
            self.sa_opts.p_unmap,
            aln_rng,
        );
        aln_problem.init();
        Ok(aln_problem)
    }

    pub fn build_solver(
        &self,
        aln_problem: &AlignmentProblem,
    ) -> Result<SimulatedAnnealing<f64, Xoshiro256PlusPlus>> {
        let problem_size = aln_problem.problem_size();
        let reanneal_best = problem_size as u64;
        let stall_best = problem_size as u64 * (problem_size as f64).ln() as u64;

        debug!(
            "Problem size: {}, reanneal_best: {}, stall_best: {}",
            problem_size, reanneal_best, stall_best
        );

        let mut rng = self.rng.lock().unwrap();
        let solver_rng = Xoshiro256PlusPlus::from_rng(&mut *rng)?;
        let solver = SimulatedAnnealing::new_with_rng(self.sa_opts.sa_init_temp, solver_rng)?
            .with_temp_func(SATempFunc::Boltzmann)
            .with_stall_accepted(self.sa_opts.sa_stall_accepted)
            .with_stall_best(stall_best)
            .with_reannealing_fixed(self.sa_opts.sa_reanneal_fixed)
            .with_reannealing_accepted(self.sa_opts.sa_reanneal_accepted)
            .with_reannealing_best(reanneal_best);
        Ok(solver)
    }
}

impl CostFunction for GenotypeProblem<'_> {
    type Param = Vec<u32>;
    type Output = f64;

    fn cost(&self, genotype: &Vec<u32>) -> Result<f64, Error> {
        let genotype: TandemRepeatGenotype = genotype.clone().into();
        let aln_problem = self.build_aln_problem(genotype.clone())?;

        let solver = self.build_solver(&aln_problem)?;

        let initial_pos = aln_problem.initial_positions();

        let executor = Executor::new(aln_problem, solver).configure(|state| {
            state
                .param(initial_pos)
                .max_iters(self.sa_opts.sa_stall_fixed)
        });

        let result = executor.run()?;
        let cost = result.state().get_best_cost();

        if let Some(best_aln_tracker) = self.best_aln_tracker {
            best_aln_tracker.update(
                genotype.clone().into(),
                cost,
                result.state().get_best_param().unwrap(),
            );
        }

        Ok(cost)
    }
}

pub struct HeterozygousGenotypeProblemWrapper<'a>(pub &'a GenotypeProblem<'a>, pub u32);
pub struct HomozygousGenotypeProblemWrapper<'a>(pub &'a GenotypeProblem<'a>);

impl CostFunction for HeterozygousGenotypeProblemWrapper<'_> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, genotype: &f64) -> Result<f64, Error> {
        let reference_genotype = self.1;
        let expanded_genotype = *genotype as u32;
        let genotype = TandemRepeatGenotype(reference_genotype, expanded_genotype);
        self.0.cost(&genotype.into())
    }
}

impl CostFunction for HomozygousGenotypeProblemWrapper<'_> {
    type Param = f64;
    type Output = f64;

    fn cost(&self, genotype: &f64) -> Result<f64, Error> {
        let expanded_genotype = *genotype as u32;
        let genotype = TandemRepeatGenotype(expanded_genotype, expanded_genotype);
        self.0.cost(&genotype.into())
    }
}

#[derive(Debug, Clone)]
pub struct GenotypeProblemDefinition {
    pub locus: TandemRepeatLocus,
    pub reference: TandemRepeatReference,
    pub reads: Vec<Read>,
    pub pairs_map: HashMap<usize, usize>,
    pub num_lowest_distances: usize,
}

impl GenotypeProblemDefinition {
    pub fn new(
        locus: TandemRepeatLocus,
        reference: TandemRepeatReference,
        num_lowest_distances: usize,
    ) -> Self {
        let reads = Vec::new();
        let mate_map = HashMap::new();
        Self {
            locus,
            reference,
            reads,
            pairs_map: mate_map,
            num_lowest_distances,
        }
    }

    pub fn add_read(&mut self, read: Read) -> Result<()> {
        let read_idx = self.reads.len();
        for i in 0..self.reads.len() {
            if read.qname == self.reads[i].qname {
                self.pairs_map.insert(read_idx, i);
                self.pairs_map.insert(i, read_idx);
            }
        }

        self.reads.push(read);
        Ok(())
    }

    pub fn create_position_generators(&self) -> Vec<RepeatAlignmentPositionSetGenerator> {
        self.reads
            .iter()
            .map(|r| r.create_position_generator(&self.reference, self.num_lowest_distances))
            .collect()
    }

    pub fn motif_len(&self) -> usize {
        self.reference.motif.len()
    }
}

#[derive(Deserialize)]
struct TempGenotypeProblemDefinition {
    id: String,
    contig: String,
    start: i64,
    end: i64,
    motif: Sequence,
    left_flank: Sequence,
    right_flank: Sequence,
    reads: Vec<Read>,
    num_lowest_distances: usize,
    mates: HashMap<usize, usize>,
}

// Custom serialization
impl Serialize for GenotypeProblemDefinition {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("OptimizationContext", 7)?; // Adjust the number of fields accordingly
        state.serialize_field("id", &self.locus.id)?;
        state.serialize_field("contig", &self.locus.contig)?;
        state.serialize_field("start", &self.locus.start)?;
        state.serialize_field("end", &self.locus.end)?;
        state.serialize_field("motif", &self.locus.motif)?;
        state.serialize_field("left_flank", &self.reference.left_flank_seq)?;
        state.serialize_field("right_flank", &self.reference.right_flank_seq)?;
        state.serialize_field("reads", &self.reads)?;
        state.serialize_field("num_lowest_distances", &self.num_lowest_distances)?;
        state.serialize_field("mates", &self.pairs_map)?;
        state.end()
    }
}

impl<'de> Deserialize<'de> for GenotypeProblemDefinition {
    fn deserialize<D>(deserializer: D) -> Result<GenotypeProblemDefinition, D::Error>
    where
        D: Deserializer<'de>,
    {
        let temp = TempGenotypeProblemDefinition::deserialize(deserializer)?;

        // Construct OptimizationContext from the temporary struct
        Ok(GenotypeProblemDefinition {
            locus: TandemRepeatLocus {
                id: temp.id,
                contig: temp.contig,
                start: temp.start,
                end: temp.end,
                motif: temp.motif.clone(),
            },
            reference: TandemRepeatReference {
                motif: temp.motif,
                left_flank_seq: temp.left_flank,
                right_flank_seq: temp.right_flank,
            },
            reads: temp.reads,
            pairs_map: temp.mates,
            num_lowest_distances: temp.num_lowest_distances,
        })
    }
}

#[derive(Clone, Debug)]
pub struct TandemRepeatGenotype(pub u32, pub u32);

// Implement Genotype into tuple
impl From<TandemRepeatGenotype> for (u32, u32) {
    fn from(g: TandemRepeatGenotype) -> Self {
        (g.0, g.1)
    }
}

// Implement Genotype into Vec<u32>
impl From<TandemRepeatGenotype> for Vec<u32> {
    fn from(g: TandemRepeatGenotype) -> Self {
        vec![g.0, g.1]
    }
}

impl From<Vec<u32>> for TandemRepeatGenotype {
    fn from(vec: Vec<u32>) -> Self {
        assert!(vec.len() == 2);
        TandemRepeatGenotype(vec[0], vec[1])
    }
}

// Implement Display for Genotype
impl std::fmt::Display for TandemRepeatGenotype {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}/{}", self.0, self.1)
    }
}

impl Serialize for TandemRepeatGenotype {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let s = format!("{}/{}", self.0, self.1);
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for TandemRepeatGenotype {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        let parts: Vec<&str> = s.split('/').collect();
        if parts.len() != 2 {
            return Err(D::Error::custom(
                "expected two integers separated by a slash",
            ));
        }
        let part0 = parts[0].parse::<u32>().map_err(D::Error::custom)?;
        let part1 = parts[1].parse::<u32>().map_err(D::Error::custom)?;
        Ok(TandemRepeatGenotype(part0, part1))
    }
}

#[derive(Clone, Debug)]
pub struct GenotypeInterval(pub (u32, u32), pub (u32, u32));

impl std::fmt::Display for GenotypeInterval {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}-{}/{}-{}", self.0 .0, self.0 .1, self.1 .0, self.1 .1)
    }
}

impl Serialize for GenotypeInterval {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let s = format!("{}-{}/{}-{}", self.0 .0, self.0 .1, self.1 .0, self.1 .1);
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for GenotypeInterval {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let s = String::deserialize(deserializer)?;
        let parts: Vec<&str> = s.split(['-', '/']).collect();
        if parts.len() != 4 {
            return Err(D::Error::custom(
                "expected four integers separated by hyphens and a slash",
            ));
        }
        let part0_0 = parts[0].parse::<u32>().map_err(D::Error::custom)?;
        let part0_1 = parts[1].parse::<u32>().map_err(D::Error::custom)?;
        let part1_0 = parts[2].parse::<u32>().map_err(D::Error::custom)?;
        let part1_1 = parts[3].parse::<u32>().map_err(D::Error::custom)?;
        Ok(GenotypeInterval((part0_0, part0_1), (part1_0, part1_1)))
    }
}

#[derive(Serialize, Clone)]
pub struct LocusGenotypeResult {
    pub reference_region: String,
    pub repeat_unit: Sequence,
    pub genotype: TandemRepeatGenotype,
    pub genotype_conf_int: GenotypeInterval,
    pub cost: f64,
    pub alignment: Option<Vec<AlignmentPosition>>,
}

/// Calculates an initial estimate of the genotype and a lower and upper bound to narrow
/// down the search space.
pub fn get_genotype_search_params(
    motif_len: u32,
    read_len: u32,
    depth_mean: f64,
    depth_sd: f64,
    num_irrs: usize,
) -> (u32, u32, u32) {
    let depth_minus_2sd =
        (depth_mean - INITIAL_ESTIMATE_DEPTH_SD_BUFFER * depth_sd).max(INITIAL_ESTIMATE_MIN_DEPTH);
    let depth_plus_2sd = depth_mean + INITIAL_ESTIMATE_DEPTH_SD_BUFFER * depth_sd;

    let genotype_init = stat_genotype_estimate(num_irrs, motif_len, read_len, depth_mean);
    let genotype_min = stat_genotype_estimate(num_irrs, motif_len, read_len, depth_plus_2sd)
        .min(genotype_init - 1);
    let genotype_max = stat_genotype_estimate(num_irrs, motif_len, read_len, depth_minus_2sd)
        .max(genotype_init + 1);

    (genotype_init, genotype_min, genotype_max)
}

// Calculates the initial estimate for the number of repeats using a
// closed-form solution:
// (((# of IRRs * read length) / mean coverage) + read length - 1) / motif length
pub fn stat_genotype_estimate(num_irrs: usize, motif_len: u32, read_len: u32, depth: f64) -> u32 {
    let initial_len_estimate =
        ((num_irrs as f64 * read_len as f64) / depth) + read_len as f64 - 1.0;
    let initial_estimate = initial_len_estimate / motif_len as f64;
    initial_estimate.round() as u32
}
