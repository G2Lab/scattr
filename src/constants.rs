// Suffixes for output files
pub const OUT_SUFFIX_STATS: &str = "stats.json";
pub const OUT_SUFFIX_BAG_UNSORTED: &str = "bag.unsorted.bam";
pub const OUT_SUFFIX_BAG: &str = "bag.bam";
pub const OUT_SUFFIX_DEFS: &str = "defs.json";
pub const OUT_SUFFIX_GENOTYPES: &str = "genotypes.json";
pub const OUT_SUFFIX_ALIGN: &str = "align.json";
pub const OUT_SUFFIX_DEPTH_DISTR_PLOT: &str = "depth_distr.png";
pub const OUT_SUFFIX_INSERT_DISTR_PLOT: &str = "insert_distr.png";

// Stats file parameters
pub const AUTOSOME_CONTIG_NAMES: &str =  "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22";

// SAM tags for bag BAM file
pub const SAM_ID_TAG: &[u8; 2] = b"ZI";
pub const SAM_READ_TYPE_TAG: &[u8; 2] = b"ZT";
pub const SAM_READ_PAIR_TYPE_TAG: &[u8; 2] = b"ZP";

pub const DEFAULT_KMER_SIZE: usize = 8;

// Genotyping parameters
pub const INITIAL_ESTIMATE_READ_LENGTH: usize = 150;
pub const INITIAL_ESTIMATE_MIN_DEPTH: f64 = 1e-6;
pub const INITIAL_ESTIMATE_DEPTH_SD_BUFFER: f64 = 0.25;

pub const DEFAULT_N_BOOTSTRAPS: u64 = 10;
pub const GENOTYPE_CONF_INT_LOWER_PERCENTILE: f64 = 2.5;
pub const GENOTYPE_CONF_INT_UPPER_PERCENTILE: f64 = 97.5;

pub const DEFAULT_GSS_ITERS_MAX: u64 = 250;

// Alignment parameters
pub const DEFAULT_SA_BASE_ERROR_RATE: f64 = 1.0 / 1000.0;
pub const DEFAULT_SA_INIT_TEMP: f64 = 1.0;

pub const DEFAULT_SA_STALL_FIXED: u64 = u64::MAX;
pub const DEFAULT_SA_STALL_ACCEPTED: u64 = u64::MAX;
pub const DEFAULT_SA_STALL_BEST: u64 = 20_000;

pub const DEFAULT_SA_REANNEAL_FIXED: u64 = u64::MAX;
pub const DEFAULT_SA_REANNEAL_ACCEPTED: u64 = u64::MAX;
pub const DEFAULT_SA_REANNEAL_BEST: u64 = 1000;

pub const DEFAULT_FIXED_READ_THRESHOLD: usize = 4;
pub const DEFAULT_INSERT_SD_TOLERANCE: f64 = 3.0;

pub const DEFAULT_UNMAPPED_PROB: f64 = 0.01;
