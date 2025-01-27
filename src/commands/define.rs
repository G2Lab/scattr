use crate::{
    catalog::TandemRepeatCatalog,
    constants::{DEFAULT_KMER_SIZE, OUT_SUFFIX_BAG},
    extract::{get_locus_id, Read},
    genotype::GenotypeProblemDefinition,
    reference::TandemRepeatReference,
    util::Utf8String,
};
use anyhow::{anyhow, Result};
use bio::io::fasta;
use clap::{arg, Parser};
use rayon::prelude::*;
use rust_htslib::bam::{self, Read as BAMRead};
use std::{
    collections::HashMap,
    fs::File,
    path::{Path, PathBuf},
    sync::{Arc, Mutex},
};

use log::{info, warn};

/// Generates definitions of the optimization problems needed to genotype the tandem repeat loci
#[derive(Parser)]
pub struct DefineCommandArgs {
    /// Path to TR catalog TSV file
    pub catalog: PathBuf,

    /// Path to alignment's reference genome
    pub reference: PathBuf,

    /// Path to bag of reads (BAM or CRAM)
    #[arg(long)]
    pub bag: Option<PathBuf>,

    /// Only keep the top n positions with the lowest edit distance for flanking reads
    #[arg(short = 'n', default_value_t = 1)]
    pub num_lowest_distances: usize,

    /// Flank length to consider for edit distance calculation
    /// (sane default: anything much larger than the fragment length)
    #[arg(short = 'f', long, default_value_t = 1500)]
    pub flank_len: u32,
}

pub fn run_define_command(args: &DefineCommandArgs, output_prefix: &Path) -> Result<()> {
    info!("Loading catalog");
    let catalog = TandemRepeatCatalog::from_path(&args.catalog, DEFAULT_KMER_SIZE)?;
    info!("Loaded {} loci from catalog", catalog.len());

    info!("Loading reference genome");
    let mut fasta = fasta::IndexedReader::from_file(&args.reference)?;

    info!("Loading bag of reads into memory");
    let bag_path = args
        .bag
        .clone()
        .unwrap_or(output_prefix.with_extension(OUT_SUFFIX_BAG));
    let mut reader = bam::Reader::from_path(bag_path)?;
    let records: Vec<bam::Record> = reader.records().map(|r| r.unwrap()).collect();
    info!("Loaded {} reads into bag", records.len());

    let header = reader.header().clone();
    let tid_to_name_map = header
        .target_names()
        .iter()
        .map(|s| (header.tid(s).unwrap() as i32, s.to_string()))
        .chain(std::iter::once((-1, "*".to_string())))
        .collect::<HashMap<i32, String>>();

    info!("Constructing tandem repeat references");
    let prob_def_map: HashMap<String, Arc<Mutex<GenotypeProblemDefinition>>> = {
        let mut prob_def_map = HashMap::new();
        for locus in catalog.iter() {
            let reference = TandemRepeatReference::from_fasta(
                &mut fasta,
                locus,
                args.flank_len,
                args.flank_len,
            )?;

            prob_def_map.insert(
                locus.id.clone(),
                Arc::new(Mutex::new(GenotypeProblemDefinition::new(
                    locus.clone(),
                    reference,
                    args.num_lowest_distances,
                ))),
            );
        }
        prob_def_map
    };

    info!("Processing reads");

    records.par_iter().try_for_each(|record| {
        let locus_id = get_locus_id(record)?;

        if let Some(prob_def_mutex) = prob_def_map.get(locus_id) {
            let mut prob_def = prob_def_mutex
                .lock()
                .map_err(|_| anyhow!("Could not lock problem definition mutex"))?;

            let record_contig = tid_to_name_map[&record.tid()].clone();

            let read = Read::from_record(
                record,
                &record_contig,
                &prob_def.reference,
                args.num_lowest_distances,
            )?;

            prob_def.add_read(read)?;
        } else {
            warn!(
                "Read {} is tagged with locus {}, but it's not in the catalog. Skipping...",
                record.qname().to_string(),
                locus_id
            );
        }

        anyhow::Ok(())
    })?;

    let sim_data_vec: Vec<GenotypeProblemDefinition> = prob_def_map
        .into_values()
        .map(|prob_def_mutex| {
            let prob_def = prob_def_mutex
                .lock()
                .map_err(|_| anyhow!("Could not lock problem definition mutex"))?
                .clone();
            anyhow::Ok(prob_def)
        })
        .map(|r| r.unwrap())
        .collect();

    let output_path = output_prefix.with_extension("defs.json");
    let writer = File::create(output_path)?;
    serde_json::to_writer_pretty(writer, &sim_data_vec)?;

    Ok(())
}
