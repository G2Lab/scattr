use crate::{
    catalog::TandemRepeatCatalog,
    constants::{DEFAULT_KMER_SIZE, OUT_SUFFIX_BAG, OUT_SUFFIX_BAG_UNSORTED},
    extract::{extract_bag_of_reads, RepeatPurityScoreParams},
};
use anyhow::Ok;
use bitcode::{Decode, Encode};
use clap::{arg, Parser};
use rust_htslib::bam::{self, Read};
use std::{
    io::Write,
    path::{Path, PathBuf},
    process::Command,
};

/// Extracts reads from an alignment file which are likely to be informative for genotyping
#[derive(Parser)]
pub struct ExtractCommandArgs {
    /// Path to alignment file (BAM or CRAM)
    pub alignment: PathBuf,

    /// Path to tandem repeat catalog TSV file
    pub catalog: PathBuf,

    /// K-mer size used for matching reads with candidate tandem repeats
    #[arg(short = 'k', default_value_t = DEFAULT_KMER_SIZE)]
    pub kmer_size: usize,

    /// Number of threads to use for reading the alignment file (htslib-specific)
    #[arg(short = '@', default_value_t = 1)]
    pub htslib_read_threads: usize,

    /// Size of the flanking regions to extract overlapping reads from
    /// (sane default: fragment mean + 3 * fragment sd - 2 * read length)
    #[arg(short = 'e', long, default_value_t = 300)]
    pub extension_length: u32,

    /// Flanking read minimum mapping quality
    #[arg(long = "mapq-f", default_value_t = 50)]
    pub min_flank_mapq: u8,

    /// In-repeat read maximum mapping quality
    #[arg(long = "mapq-i", default_value_t = 40)]
    pub max_irr_mapq: u8,

    /// In-repeat read minimum weighted purity score
    #[arg(short = 's', long, default_value = "0.9")]
    pub irr_score_min: f32,

    /// Repeat purity score base quality minimum
    #[arg(short = 'q', long, default_value = "20")]
    pub base_qual_min: u8,

    /// Repeat purity score match weight
    #[arg(short = 'm', long, default_value = "1.0")]
    pub match_weight: f32,

    /// Repeat purity score mismatch weight
    #[arg(short = 'x', long, default_value = "-1.0")]
    pub mismatch_weight: f32,

    /// Repeat purity score low quality mismatch weight
    #[arg(short = 'z', long, default_value = "0.5")]
    pub low_qual_mismatch_weight: f32,

    /// Sort and index the output file (requires samtools in PATH)
    #[arg(short = 'i', long)]
    pub sort_and_index: bool,
}

pub fn run_extract_command(args: &ExtractCommandArgs, output_prefix: &Path) -> anyhow::Result<()> {
    // Open TR catalog
    info!("Loading catalog");
    let catalog = TandemRepeatCatalog::from_path(&args.catalog, args.kmer_size)?;
    info!("Loaded {} loci from catalog", catalog.len());

    // Open input alignments file
    info!("Opening input alignments file");
    let mut reader: bam::IndexedReader = bam::IndexedReader::from_path(&args.alignment)?;
    reader.set_threads(args.htslib_read_threads)?;

    // Open output file
    info!("Opening output file");
    let extension = if args.sort_and_index {
        OUT_SUFFIX_BAG_UNSORTED
    } else {
        OUT_SUFFIX_BAG
    };
    let output_path = output_prefix.with_extension(extension);
    let header = bam::Header::from_template(reader.header());
    let mut writer = bam::Writer::from_path(output_path, &header, bam::Format::Bam)?;

    // Create param structs
    let purity_score_params = RepeatPurityScoreParams {
        irr_score_min: args.irr_score_min,
        base_qual_min: args.base_qual_min,
        match_weight: args.match_weight,
        mismatch_weight: args.mismatch_weight,
        low_qual_mismatch_weight: args.low_qual_mismatch_weight,
    };

    info!("Extracting bag of reads");
    let locus_reads = extract_bag_of_reads(
        &mut reader,
        &mut writer,
        catalog,
        args.extension_length,
        args.min_flank_mapq,
        args.max_irr_mapq,
        &purity_score_params,
    )?;

    drop(writer);

    // Write locus reads to file using bitcode library
    info!("Writing locus-read associations to file");
    let mut output = std::fs::File::create(output_prefix.with_extension("bag.bin"))?;
    let encoded = bitcode::encode(&locus_reads);
    output.write_all(&encoded)?;

    if args.sort_and_index {
        info!("Sorting and indexing output file");
        // Sort
        Command::new("samtools")
            .args(vec![
                "sort",
                "-o",
                output_prefix
                    .with_extension(OUT_SUFFIX_BAG)
                    .to_str()
                    .unwrap(),
                "-O",
                "bam",
                output_prefix
                    .with_extension(OUT_SUFFIX_BAG_UNSORTED)
                    .to_str()
                    .unwrap(),
                "-T",
                output_prefix
                    .with_extension("samtools_tmp")
                    .to_str()
                    .unwrap(),
                "-@",
                &args.htslib_read_threads.to_string(),
            ])
            .output()?;

        // Index
        Command::new("samtools")
            .args(vec![
                "index",
                output_prefix
                    .with_extension(OUT_SUFFIX_BAG)
                    .to_str()
                    .unwrap(),
                "-@",
                &args.htslib_read_threads.to_string(),
            ])
            .output()?;

        // Remove unsorted output file
        std::fs::remove_file(output_prefix.with_extension(OUT_SUFFIX_BAG_UNSORTED))?;
    }

    Ok(())
}
