use crate::sequence::Sequence;
use aho_corasick::{AhoCorasick, PatternID};
use anyhow::{Context, Ok, Result};
use bio::data_structures::interval_tree::IntervalTree;
use indexmap::IndexMap;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt;
use std::ops::Index;
use std::path::Path;

pub const KMER_SIZE: usize = 17;

/// Represents a single locus in a tandem repeat catalog
#[derive(Debug, PartialEq, Eq, Deserialize, Serialize, Clone, Hash)]
pub struct TandemRepeatLocus {
    pub id: String,
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub motif: Sequence,
}

impl fmt::Display for TandemRepeatLocus {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

/// Represents a catalog of tandem repeat loci. Loci are inserted
/// into an interval tree for each contig, allowing for efficient
/// lookup of loci overlapping a given interval.
#[derive(Debug, Clone)]
pub struct TandemRepeatCatalog {
    loci: Vec<TandemRepeatLocus>,
    /// Map from contig name to interval tree of locus indices
    trees: HashMap<String, IntervalTree<i64, usize>>,
    /// Map of motif to locus indices
    motif_map: HashMap<Vec<u8>, HashSet<usize>>,
    /// Map of locus IDs to locus indices
    id_map: HashMap<String, usize>,
}

impl TandemRepeatCatalog {
    pub fn new() -> Self {
        Self {
            loci: Vec::new(),
            trees: HashMap::new(),
            motif_map: HashMap::new(),
            id_map: HashMap::new(),
        }
    }

    pub fn contigs(&self) -> Vec<&String> {
        self.trees.keys().collect()
    }

    pub fn insert(&mut self, locus: TandemRepeatLocus) {
        let locus_index = self.loci.len();

        // Insert locus into interval tree
        self.trees
            .entry(locus.contig.clone())
            .or_default()
            .insert(locus.start..locus.end, locus_index);

        // Insert locus into motif map

        // let mut motif = locus.motif.clone();
        // let mut motif_rc = motif.revcomp();

        // for _ in 0..locus.motif.len() {
        //     self.motif_map
        //         .entry(motif.clone())
        //         .or_default()
        //         .push(locus_index);
        //     self.motif_map
        //         .entry(motif_rc.clone())
        //         .or_default()
        //         .push(locus_index);
        //     motif.rotate_left(1);
        //     motif_rc.rotate_left(1);
        // }

        let mut motif = locus.motif.as_ref().to_vec();
        let mut motif_rc = locus.motif.revcomp().as_ref().to_vec();

        for _ in 0..locus.motif.len() {
            let motif_padded: Vec<u8> = motif.iter().cycle().take(KMER_SIZE).cloned().collect();
            let motif_rc_padded: Vec<u8> =
                motif_rc.iter().cycle().take(KMER_SIZE).cloned().collect();

            self.motif_map
                .entry(motif_padded.clone())
                .or_default()
                .insert(locus_index);
            self.motif_map
                .entry(motif_rc_padded.clone())
                .or_default()
                .insert(locus_index);
            motif.rotate_left(1);
            motif_rc.rotate_left(1);
        }

        // Insert locus into ID map
        self.id_map.insert(locus.id.clone(), locus_index);

        self.loci.push(locus);
    }

    pub fn get(&self, id: &str) -> Option<&TandemRepeatLocus> {
        let index = self.id_map.get(id)?;
        Some(&self.loci[*index])
    }

    pub fn find(&self, contig: &str, start: i64, end: i64) -> Vec<&TandemRepeatLocus> {
        let tree = self.trees.get(contig);
        if let Some(tree) = tree {
            let indices = tree
                .find(start..end)
                .map(|entry| entry.data())
                .collect::<Vec<_>>();
            return indices.iter().map(|&i| &self.loci[*i]).collect();
        } else {
            Vec::new()
        }
    }

    // pub fn build_ac(&self) -> Result<AhoCorasick> {
    //     let patterns: Vec<Vec<u8>> = self.motif_map.keys().map(|s| s.as_ref().to_vec()).collect();
    //     let ac = AhoCorasick::new(patterns)
    //         .context("Failed to create AhoCorasick automaton for catalog")?;
    //     Ok(ac)
    // }

    pub fn find_motifs(&self, seq: &[u8]) -> Vec<&TandemRepeatLocus> {
        let mut result = HashSet::new();
        let mut i = 0;
        while i + KMER_SIZE <= seq.len() {
            let kmer = &seq[i..i + KMER_SIZE];
            let locus_idxs = self.motif_map.get(kmer);
            if let Some(locus_idxs) = locus_idxs {
                for &locus_idx in locus_idxs {
                    result.insert(locus_idx);
                }
            }
            i += KMER_SIZE;
        }
        result.iter().map(|&i| &self.loci[i]).collect()
    }

    pub fn iter(&self) -> impl Iterator<Item = &TandemRepeatLocus> {
        self.loci.iter()
    }

    pub fn len(&self) -> usize {
        self.loci.len()
    }

    pub fn is_empty(&self) -> bool {
        self.loci.is_empty()
    }

    pub fn from_path<P: AsRef<Path> + std::fmt::Debug + ?Sized>(path: &P) -> Result<Self> {
        let mut catalog = TandemRepeatCatalog::new();
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .context(format!("Failed to open TR catalog: {:?}", path))?;
        for r in reader.deserialize() {
            let locus: TandemRepeatLocus = r.context("Failed to parse TR catalog entry")?;
            catalog.insert(locus);
        }
        Ok(catalog)
    }

    pub fn write(&self, path: &str) -> Result<()> {
        let mut writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_path(path)
            .context(format!("Failed to open TR catalog: {}", path))?;
        for locus in self.iter() {
            writer.serialize(locus)?;
        }
        Ok(())
    }
}

impl Default for TandemRepeatCatalog {
    fn default() -> Self {
        Self::new()
    }
}

// #[cfg(test)]
// mod tests {
//     use std::collections::HashSet;

//     use aho_corasick::PatternID;

//     use super::*;

//     const READ_CATALOG_PATH: &str = "tests/read_catalog.tsv";
//     const WRITE_CATALOG_PATH: &str = "tests/write_catalog.tsv";

//     const LOCUS1: (&str, &str, i64, i64, &str) = ("chr1_100_200", "chr1", 100, 200, "CAG");
//     const LOCUS2: (&str, &str, i64, i64, &str) = ("chr1_300_400", "chr1", 300, 400, "CG");
//     const LOCUS3: (&str, &str, i64, i64, &str) = ("chrX_1000_2000", "chrX", 1000, 2000, "TAT");

//     fn make_locus(
//         (id, contig, start, end, motif): (&str, &str, i64, i64, &str),
//     ) -> TandemRepeatLocus {
//         TandemRepeatLocus {
//             id: id.to_string(),
//             contig: contig.to_string(),
//             start,
//             end,
//             motif: Sequence::from(motif),
//         }
//     }

//     fn get_test_loci() -> Vec<TandemRepeatLocus> {
//         vec![make_locus(LOCUS1), make_locus(LOCUS2), make_locus(LOCUS3)]
//     }

//     #[test]
//     fn test_catalog_insert() {
//         let mut catalog = TandemRepeatCatalog::new();
//         let locus = make_locus(LOCUS1);
//         catalog.insert(locus.clone());
//         assert_eq!(catalog.len(), 1, "Expected 1 locus in catalog");
//         assert_eq!(
//             catalog.find(&locus.contig, locus.start, locus.end)[0],
//             &locus,
//             "Expected to find locus {}",
//             locus
//         );
//         assert_eq!(
//             catalog.motif_map.len(),
//             6,
//             "Expected to find 6 motifs in motif map"
//         );
//         assert_eq!(catalog.id_map.len(), 1, "Expected to find 1 loci in id map");
//         assert_eq!(catalog.loci.len(), 1, "Expected to find 1 loci");
//         assert_eq!(catalog.trees.len(), 1, "Expected to find 1 contig in trees");

//         assert_eq!(
//             catalog
//                 .motif_map
//                 .keys()
//                 .cloned()
//                 .collect::<HashSet<Sequence>>(),
//             vec![
//                 Sequence::from("CAG"),
//                 Sequence::from("CTG"),
//                 Sequence::from("AGC"),
//                 Sequence::from("GCT"),
//                 Sequence::from("GCA"),
//                 Sequence::from("TGC")
//             ]
//             .into_iter()
//             .collect::<HashSet<_>>(),
//             "Expected to specific motifs in motif map",
//         );
//     }

//     #[test]
//     fn test_catalog_from_tsv() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         assert_eq!(catalog.len(), 3, "Expected 3 loci in catalog");
//         assert_eq!(catalog.trees.len(), 2, "Expected 2 contigs in catalog");
//     }

//     #[test]
//     fn test_catalog_to_tsv() {
//         let locus = make_locus(LOCUS1);
//         let mut catalog = TandemRepeatCatalog::new();
//         catalog.insert(locus.clone());
//         catalog.write(WRITE_CATALOG_PATH).unwrap();
//         let data = std::fs::read_to_string(WRITE_CATALOG_PATH).unwrap();
//         assert_eq!(
//             data,
//             format!(
//                 "id\tcontig\tstart\tend\tmotif\n{}\t{}\t{}\t{}\t{}\n",
//                 locus.id, locus.contig, locus.start, locus.end, locus.motif
//             ),
//             "Expected catalog to be written to file"
//         );
//     }

//     #[test]
//     fn test_catalog_find_exact() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();

//         for locus in loci {
//             let result = catalog.find(&locus.contig, locus.start, locus.end);
//             assert_eq!(
//                 result.len(),
//                 1,
//                 "Expected to find exactly 1 locus for {}",
//                 locus
//             );
//             assert_eq!(result[0], &locus, "Expected to find the locus {}", locus);
//         }
//     }

//     #[test]
//     fn test_catalog_find_left_overlap() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();

//         for locus in loci {
//             let result = catalog.find(&locus.contig, locus.start - 10, locus.end - 10);
//             assert_eq!(
//                 result.len(),
//                 1,
//                 "Expected to find exactly 1 locus for {}",
//                 locus
//             );
//             assert_eq!(result[0], &locus, "Expected to find the locus {}", locus);
//         }
//     }

//     #[test]
//     fn test_catalog_find_right_overlap() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();

//         for locus in loci {
//             let result = catalog.find(&locus.contig, locus.start + 10, locus.end + 10);
//             assert_eq!(
//                 result.len(),
//                 1,
//                 "Expected to find exactly 1 locus for {}",
//                 locus
//             );
//             assert_eq!(result[0], &locus, "Expected to find the locus {}", locus);
//         }
//     }

//     #[test]
//     fn test_catalog_find_edge_case_left() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();

//         for locus in loci {
//             let result = catalog.find(&locus.contig, locus.start - 10, locus.start);
//             assert_eq!(
//                 result.len(),
//                 0,
//                 "Expected to find no loci in {}:{}-{}",
//                 locus.contig,
//                 locus.start - 10,
//                 locus.start
//             );
//         }
//     }

//     #[test]
//     fn test_catalog_find_edge_case_right() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();

//         for locus in loci {
//             let result = catalog.find(&locus.contig, locus.end, locus.end + 10);
//             assert_eq!(
//                 result.len(),
//                 0,
//                 "Expected to find no loci in {}:{}-{}",
//                 locus.contig,
//                 locus.end,
//                 locus.end + 10
//             );
//         }
//     }

//     #[test]
//     fn test_catalog_build_ac_exact() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();
//         let ac = catalog.build_ac().unwrap();
//         for locus in loci {
//             let result: Vec<PatternID> = ac
//                 .find_iter(locus.motif.as_ref())
//                 .map(|mat| mat.pattern())
//                 .collect();
//             assert_eq!(
//                 result.len(),
//                 1,
//                 "Expected Aho-Corasick automaton to find exactly 1 match for motif {}",
//                 locus.motif
//             );
//             let p = result[0];
//             assert_eq!(
//                 catalog.find_pattern_id(&p)[0],
//                 &locus,
//                 "Expected Aho-Corasick automaton to find locus {}",
//                 locus
//             );
//         }
//     }

//     #[test]
//     fn test_catalog_build_ac_shifted() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();
//         let ac = catalog.build_ac().unwrap();
//         for locus in loci {
//             let mut motif = locus.motif.clone();
//             motif.rotate_left(1);
//             let result: Vec<PatternID> = ac
//                 .find_iter(locus.motif.as_ref())
//                 .map(|mat| mat.pattern())
//                 .collect();
//             assert_eq!(
//                 result.len(),
//                 1,
//                 "Expected Aho-Corasick automaton to find exactly 1 match for motif {}",
//                 motif
//             );
//             let p = result[0];
//             assert_eq!(
//                 catalog.find_pattern_id(&p)[0],
//                 &locus,
//                 "Expected Aho-Corasick automaton to find locus {}",
//                 locus
//             );
//         }
//     }

//     #[test]
//     fn test_catalog_iter() {
//         let catalog = TandemRepeatCatalog::from_path(READ_CATALOG_PATH).unwrap();
//         let loci = get_test_loci();
//         let result: Vec<TandemRepeatLocus> = catalog.iter().cloned().collect();
//         assert_eq!(loci, result, "Expected catalog iterator to match loci");
//     }
// }
