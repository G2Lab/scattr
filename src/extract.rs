use crate::{
    catalog::TandemRepeatCatalog,
    constants::{SAM_ID_TAG, SAM_READ_PAIR_TYPE_TAG, SAM_READ_TYPE_TAG},
    positions::{HammingDistance, RepeatAlignmentPosition, RepeatAlignmentPositionSetGenerator},
    records::{ReadKind, ReadOrder, ReadPairId, RecordManager},
    reference::TandemRepeatReference,
    sequence::Sequence,
    util::Utf8String,
};
use anyhow::anyhow;
use anyhow::{Context, Ok, Result};
use bio::alphabets::dna::revcomp;
use indicatif::{ProgressBar, ProgressStyle};
use rust_htslib::bam::FetchDefinition;
use rust_htslib::{
    bam,
    bam::{record::Aux, Read as BAMRead, Record},
};
use serde::ser::SerializeSeq;
use serde::{
    de::{self, SeqAccess, Visitor},
    Deserialize, Deserializer, Serialize, Serializer,
};
use std::fmt;
use std::hash::Hash;
use std::str::FromStr;
use std::{collections::HashMap, time::Duration};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct LocusFlankSizes {
    pub left: HashMap<LocusId, u32>,
    pub right: HashMap<LocusId, u32>,
}

#[derive(Clone)]
pub struct RepeatPurityScoreParams {
    pub irr_score_min: f32,
    pub base_qual_min: u8,
    pub match_weight: f32,
    pub mismatch_weight: f32,
    pub low_qual_mismatch_weight: f32,
}

trait RecordRange {
    fn start(&self) -> i64;
    fn stop(&self) -> i64;
}

impl RecordRange for Record {
    fn start(&self) -> i64 {
        self.pos()
    }

    fn stop(&self) -> i64 {
        // output is 0-indexed, exclusive
        self.start() + self.seq_len() as i64
    }
}

#[derive(Debug, Serialize, Deserialize, Clone, Hash)]
pub struct Read {
    pub qname: String,
    pub source_contig: String,
    pub source_pos: i64,
    pub mapq: u8,
    pub seq: Sequence,
    pub seq_len: u32,
    pub is_first_in_template: bool,
    #[serde(
        serialize_with = "serialize_vec_pos_tuple",
        deserialize_with = "deserialize_vec_pos_tuple"
    )]
    pub position_set: Vec<(RepeatAlignmentPosition, HammingDistance)>,
}

impl Read {
    pub fn from_record(
        record: &Record,
        record_contig: &str,
        reference: &TandemRepeatReference,
        num_lowest_distances: usize,
    ) -> Result<Self> {
        let qname = String::from_utf8(record.qname().to_vec())?;
        let source_contig = record_contig.to_string();
        let source_pos = record.pos();
        let mapq: u8 = record.mapq();
        let seq_len = record.seq_len() as u32;
        let is_first_in_template = record.is_first_in_template();
        let mut seq = record.seq().as_bytes();

        if record.is_reverse() {
            seq = revcomp(seq);
        }

        let generator = RepeatAlignmentPositionSetGenerator::new(
            &seq,
            reference.motif.as_ref(),
            reference.left_flank_seq.as_ref(),
            reference.right_flank_seq.as_ref(),
            num_lowest_distances,
        );
        let relative_pos = generator.relative_position_set();

        Ok(Self {
            qname,
            source_pos,
            source_contig,
            mapq,
            seq: seq.into(),
            seq_len,
            is_first_in_template,
            position_set: relative_pos,
        })
    }

    pub fn create_position_generator(
        &self,
        reference: &TandemRepeatReference,
        num_lowest_distances: usize,
    ) -> RepeatAlignmentPositionSetGenerator {
        RepeatAlignmentPositionSetGenerator::from_precomputed(
            self.seq.as_ref(),
            reference.motif.as_ref(),
            reference.left_flank_seq.as_ref(),
            reference.right_flank_seq.as_ref(),
            num_lowest_distances,
            self.position_set.clone(),
        )
    }

    fn is_irr(&self) -> bool {
        self.position_set
            .iter()
            .any(|(pos, _)| matches!(pos, RepeatAlignmentPosition::WithinRepeat(_, _)))
    }
}

pub type LocusId = String;

pub fn extract_bag_of_reads(
    reader: &mut bam::IndexedReader,
    writer: &mut bam::Writer,
    catalog: TandemRepeatCatalog,
    extension_length: u32,
    min_flank_mapq: u8,
    max_irr_mapq: u8,
    purity_score_params: &RepeatPurityScoreParams,
) -> Result<()> {
    let header = reader.header().clone();

    reader.fetch(FetchDefinition::All)?;

    let mut record_manager = RecordManager::new();

    let spinner = ProgressBar::new_spinner();
    spinner.set_style(
        ProgressStyle::with_template(
            "{spinner} [{elapsed_precise}] Processed {pos} reads. Reads per second: {per_sec} {msg}",
        )
        .unwrap(),
    );
    spinner.enable_steady_tick(Duration::from_secs(1));

    for record_result in reader.rc_records() {
        let record = record_result?;

        if !is_primary_read(&record) {
            continue;
        }

        let read_pair_id: ReadPairId = RecordManager::get_read_pair_id(&record);

        debug!("Processing read pair {}", read_pair_id);

        // Add the read to the record manager if it is part of a read pair
        // if record_manager.contains_read_pair_id(&read_pair_id) {
        //     debug!(
        //         "Adding read to record manager for read pair {}",
        //         read_pair_id
        //     );
        //     record_manager.add_read(record.as_ref());
        // }

        if record.mapq() <= max_irr_mapq {
            let candidate_irr_for_loci = catalog.find_motifs(&record.seq().as_bytes());
            for locus in candidate_irr_for_loci {
                if is_in_repeat_read(&record, &locus.motif, purity_score_params) {
                    debug!(
                        "Read {} is an in-repeat read for locus {}",
                        record.qname().to_string(),
                        locus.id
                    );
                    record_manager.add_read_to_locus_with_kind(
                        &locus.id,
                        record.as_ref(),
                        ReadKind::InRepeatRead,
                    );
                }
            }
        }

        if record.mapq() >= min_flank_mapq {
            let tid = record.tid();
            let contig = if tid >= 0 {
                Some(header.tid2name(tid as u32).to_string())
            } else {
                None
            };

            // Find all loci that the read counts as a flanking read for
            if let Some(contig) = contig {
                let extended_start = record.start() - extension_length as i64;
                let extended_end = record.stop() + extension_length as i64;
                let flank_for_loci = catalog.find(&contig, extended_start, extended_end);
                for locus in flank_for_loci {
                    debug!(
                        "Read {} is a flanking read for locus {}",
                        record.qname().to_string(),
                        locus.id
                    );
                    record_manager.add_read_to_locus_with_kind(
                        &locus.id,
                        record.as_ref(),
                        ReadKind::FlankingRead,
                    );
                }
            }
        }

        spinner.inc(1);
    }

    spinner.finish_and_clear();

    // For each locus, rescue mates that are mapped
    for locus in catalog.iter() {
        let incomplete_read_pairs_ids =
            record_manager.incomplete_read_pair_ids_for_locus(&locus.id);
        for read_pair_id in incomplete_read_pairs_ids {
            let read_pair = record_manager.get_read_pair_by_id(&read_pair_id).unwrap();

            let (read, read_order) = if read_pair.first.is_some() {
                (read_pair.first.unwrap(), ReadOrder::First)
            } else {
                (read_pair.last.unwrap(), ReadOrder::Last)
            };

            let read_kind = record_manager
                .get_read_kind(&locus.id, &read_pair_id, &read_order)
                .unwrap();

            // Only rescue mates for flanking reads (in-repeat reads mates should have been found)
            if read_kind == ReadKind::FlankingRead && !read.is_mate_unmapped() {
                let mate_tid = read.mtid();
                let mate_pos = read.mpos();
                if mate_tid >= 0 {
                    reader.fetch((mate_tid, mate_pos, mate_pos + 1))?;
                    for record_result in reader.rc_records() {
                        let record = record_result?;
                        if is_primary_read(&record) && record.qname().to_string() == read_pair.id {
                            record_manager.add_read_to_locus_with_kind(
                                &locus.id,
                                record.as_ref(),
                                ReadKind::FlankingRead,
                            );
                            break;
                        }
                    }
                } else {
                    warn!(
                        "Read pair {} has a mate with invalid tid {}",
                        read_pair.id, mate_tid
                    );
                }
            }
        }
    }

    // Rescue mates (for reads with mates that are unmapped)
    let incomplete_read_pair_ids = record_manager.incomplete_read_pair_ids();
    if !incomplete_read_pair_ids.is_empty() {
        reader.fetch(FetchDefinition::Unmapped)?;
        for record_result in reader.rc_records() {
            let record = record_result?;
            if is_primary_read(&record) {
                let read_pair_id = RecordManager::get_read_pair_id(record.as_ref());
                if incomplete_read_pair_ids.contains(&read_pair_id) {
                    let read_pair_locus_ids = record_manager.get_read_pair_locus_ids(&read_pair_id);
                    for locus_id in read_pair_locus_ids {
                        let read_pair = record_manager.get_read_pair_by_id(&read_pair_id).unwrap();

                        // Get the present read's kind
                        let read_order = if read_pair.first.is_some() {
                            ReadOrder::First
                        } else {
                            ReadOrder::Last
                        };
                        let read_kind = record_manager
                            .get_read_kind(&locus_id, &read_pair_id, &read_order)
                            .unwrap();

                        // Unmapped mates of flanking reads should be added as in repeat reads
                        if read_kind == ReadKind::FlankingRead {
                            record_manager.add_read_to_locus_with_kind(
                                &locus_id,
                                record.as_ref(),
                                ReadKind::InRepeatRead,
                            );
                        }
                    }
                }
            }
        }
    }

    // Write the read pairs to the output file
    for locus in catalog.iter() {
        let read_pairs = record_manager.get_locus_read_pairs(&locus.id);
        for read_pair in read_pairs {
            // Do not write incomplete read pairs
            if read_pair.len() != 2 {
                continue;
            }

            let first_read_kind =
                record_manager.get_read_kind(&locus.id, &read_pair.id, &ReadOrder::First);
            let last_read_kind =
                record_manager.get_read_kind(&locus.id, &read_pair.id, &ReadOrder::Last);

            match (first_read_kind, last_read_kind) {
                (Some(first_read_kind), Some(last_read_kind)) => {
                    let mut first_read = read_pair.first.as_ref().unwrap().clone();
                    let mut last_read = read_pair.last.as_ref().unwrap().clone();

                    let read_pair_kind = record_manager
                        .get_read_pair_kind(&locus.id, &read_pair.id)
                        .to_string();

                    clear_record_tags(&mut first_read)?;
                    clear_record_tags(&mut last_read)?;

                    first_read.push_aux(SAM_ID_TAG, Aux::String(&locus.id))?;
                    last_read.push_aux(SAM_ID_TAG, Aux::String(&locus.id))?;

                    first_read
                        .push_aux(SAM_READ_TYPE_TAG, Aux::String(&first_read_kind.to_string()))?;
                    last_read
                        .push_aux(SAM_READ_TYPE_TAG, Aux::String(&last_read_kind.to_string()))?;

                    first_read.push_aux(SAM_READ_PAIR_TYPE_TAG, Aux::String(&read_pair_kind))?;
                    last_read.push_aux(SAM_READ_PAIR_TYPE_TAG, Aux::String(&read_pair_kind))?;

                    writer.write(&first_read)?;
                    writer.write(&last_read)?;
                }
                _ => {
                    warn!(
                        "Read pair {} for locus {} is does not have their types defined",
                        read_pair.id, locus.id
                    );
                }
            }
        }
    }

    Ok(())
}

fn clear_record_tags(record: &mut Record) -> Result<()> {
    if record.aux(SAM_ID_TAG).is_ok() {
        record.remove_aux(SAM_ID_TAG)?;
    }
    if record.aux(SAM_READ_TYPE_TAG).is_ok() {
        record.remove_aux(SAM_READ_TYPE_TAG)?;
    }
    if record.aux(SAM_READ_PAIR_TYPE_TAG).is_ok() {
        record.remove_aux(SAM_READ_PAIR_TYPE_TAG)?;
    }
    Ok(())
}

pub fn is_primary_read(read: &Record) -> bool {
    !read.is_secondary() && !read.is_supplementary()
}

pub fn is_in_repeat_read(
    record: &Record,
    motif: &Sequence,
    params: &RepeatPurityScoreParams,
) -> bool {
    let seq = record.seq().as_bytes();
    let base_quals = record.qual();
    let mut motif = motif.clone();
    let mut motif_rc = motif.revcomp();
    for _ in 0..motif.len() {
        if seq_passes_purity_score(&seq, base_quals, motif.as_ref(), params)
            || seq_passes_purity_score(&seq, base_quals, motif_rc.as_ref(), params)
        {
            return true;
        }
        // Rotate motif
        motif.rotate_left(1);
        motif_rc.rotate_left(1);
    }
    false
}

pub fn seq_passes_purity_score(
    seq: &[u8],
    base_quals: &[u8],
    motif: &[u8],
    params: &RepeatPurityScoreParams,
) -> bool {
    let required_score = seq.len() as f32 * params.irr_score_min;

    let mut score = 0.0;

    for (i, (base, &qual)) in seq.iter().zip(base_quals).enumerate() {
        let motif_base = motif[i % motif.len()];
        if *base == motif_base {
            score += params.match_weight;
        } else if qual >= params.base_qual_min {
            score += params.mismatch_weight;
        } else {
            score += params.low_qual_mismatch_weight;
        }

        // Compute the maximum potential score remaining
        let remaining_bases = seq.len() - i - 1;
        let max_possible_score = score + remaining_bases as f32 * params.match_weight;

        // Exit early if it's impossible to reach the threshold
        if max_possible_score < required_score {
            return false;
        }

        // Exit early if the current score exceeds the threshold
        if score >= required_score {
            return true;
        }
    }

    false
}

pub fn purity_score(
    seq: &[u8],
    base_quals: &[u8],
    motif: &[u8],
    params: &RepeatPurityScoreParams,
) -> f32 {
    if seq.is_empty() {
        return 0.0;
    }

    let mut score = 0.0;

    for (i, base) in seq.iter().enumerate() {
        if *base == motif[i % motif.len()] {
            score += params.match_weight;
        } else if base_quals[i] >= params.base_qual_min {
            score += params.mismatch_weight;
        } else {
            score += params.low_qual_mismatch_weight;
        }
    }

    score / (seq.len() as f32)
}

pub fn get_locus_id(record: &Record) -> Result<&str> {
    if let Aux::String(value) = record.aux(SAM_ID_TAG)? {
        Ok(value)
    } else {
        Err(anyhow!(
            "Read with qname={} does not have ID tag",
            record.qname().to_string()
        ))
    }
}

// fn serialize_vec_pos_tuple<S>(vec: &[(i64, u32)], serializer: S) -> Result<S::Ok, S::Error>
// where
//     S: Serializer,
// {
//     let str_repr: String = vec
//         .iter()
//         .map(|&(pos, score)| format!("({}: {})", pos, score))
//         .collect::<Vec<String>>()
//         .join(", ");
//     serializer.serialize_str(&str_repr)
// }

// fn deserialize_vec_pos_tuple<'de, D>(deserializer: D) -> Result<Vec<(i64, u32)>, D::Error>
// where
//     D: Deserializer<'de>,
// {
//     struct VecTupleVisitor;

//     impl<'de> Visitor<'de> for VecTupleVisitor {
//         type Value = Vec<(i64, u32)>;

//         fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
//             formatter.write_str("a string representing a list of (i64, u32) tuples")
//         }

//         fn visit_str<E>(self, value: &str) -> std::result::Result<Vec<(i64, u32)>, E>
//         where
//             E: de::Error,
//         {
//             if value.len() == 0 {
//                 std::result::Result::Ok(Vec::new())
//             } else {
//                 value
//                     .split(", ")
//                     .map(|tuple_str| {
//                         let parts: Vec<&str> =
//                             tuple_str[1..tuple_str.len() - 1].split(": ").collect();
//                         if parts.len() != 2 {
//                             return Err(E::custom(format!("Invalid tuple format: {}", tuple_str)));
//                         }
//                         let pos = i64::from_str(parts[0]).map_err(de::Error::custom)?;
//                         let value = u32::from_str(parts[1]).map_err(de::Error::custom)?;
//                         std::result::Result::Ok((pos, value))
//                     })
//                     .collect()
//             }
//         }
//     }

//     deserializer.deserialize_str(VecTupleVisitor)
// }

pub fn count_irrs(reads: &[Read]) -> usize {
    reads.iter().filter(|r| r.is_irr()).count()
}

// Function to serialize a Vec<(RepeatAlignmentPosition, HammingDistance)>
fn serialize_vec_pos_tuple<S>(
    vec: &[(RepeatAlignmentPosition, HammingDistance)],
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut seq = serializer.serialize_seq(Some(vec.len()))?;
    for tuple in vec {
        let pos = serde_json::to_string(&tuple.0).unwrap();
        let s = format!("({}: {})", &pos[1..pos.len() - 1], tuple.1);
        seq.serialize_element(&s)?;
    }
    seq.end()
}

// Function to deserialize a Vec<(RepeatAlignmentPosition, HammingDistance)>
fn deserialize_vec_pos_tuple<'de, D>(
    deserializer: D,
) -> Result<Vec<(RepeatAlignmentPosition, HammingDistance)>, D::Error>
where
    D: Deserializer<'de>,
{
    struct VecPosTupleVisitor;

    impl<'de> Visitor<'de> for VecPosTupleVisitor {
        type Value = Vec<(RepeatAlignmentPosition, HammingDistance)>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str(
                "a list of strings representing (RepeatAlignmentPosition: HammingDistance) tuples",
            )
        }

        fn visit_seq<A>(self, mut seq: A) -> Result<Self::Value, A::Error>
        where
            A: SeqAccess<'de>,
        {
            let mut vec = Vec::new();

            while let Some(value) = seq.next_element::<String>()? {
                let value = value.trim_matches(|c| c == '(' || c == ')');
                let mut parts = value.split(": ");

                let position_str = parts
                    .next()
                    .ok_or_else(|| de::Error::custom("missing position part"))?;
                let distance_str = parts
                    .next()
                    .ok_or_else(|| de::Error::custom("missing distance part"))?;

                let repeat_pos: RepeatAlignmentPosition =
                    serde_json::from_str(&format!("\"{}\"", position_str)).map_err(|_| {
                        de::Error::custom("failed to parse RepeatAlignmentPosition")
                    })?;
                let hamming_distance: HammingDistance =
                    HammingDistance::from_str(distance_str).map_err(de::Error::custom)?;

                vec.push((repeat_pos, hamming_distance));
            }

            std::result::Result::Ok(vec)
        }
    }

    deserializer.deserialize_seq(VecPosTupleVisitor)
}
// #[cfg(test)]
// mod tests {
//     use super::*;

//     const TEST_PURITY_SCORE_PARAMS: RepeatPurityScoreParams = RepeatPurityScoreParams {
//         irr_score_min: 0.90,
//         base_qual_min: 20,
//         match_weight: 1.0,
//         mismatch_weight: -1.0,
//         low_qual_mismatch_weight: 0.5,
//     };

//     fn create_record(seq: &str, qual: &[u8]) -> Record {
//         let mut record = Record::new();
//         let seq = seq.as_bytes();
//         record.set("test".as_bytes(), None, seq, qual);
//         record
//     }

//     fn get_test_read(
//         seq_len: u32,
//         left_flank_pos: Vec<i64>,
//         right_flank_pos: Vec<i64>,
//         in_repeat_motif_pos: Vec<i64>,
//     ) -> Read {
//         Read {
//             qname: String::from("test_read"),
//             seq: Sequence::from(""), // Empty sequence
//             source_pos: 0,
//             seq_len: seq_len,
//             is_first_in_template: true,
//             kind: ReadKind::Flanking,
//             left_flank_pos: left_flank_pos.iter().map(|&p| (p, 0)).collect(),
//             right_flank_pos: right_flank_pos.iter().map(|&p| (p, 0)).collect(),
//             in_repeat_motif_pos: in_repeat_motif_pos.iter().map(|&p| (p, 0)).collect(),
//         }
//     }

//     #[test]
//     fn test_purity_score_empty() {
//         let seq = b"";
//         let base_quals = &vec![];
//         let motif = b"CAG";
//         assert_eq!(
//             purity_score(seq, base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             0.0
//         );
//     }

//     #[test]
//     fn test_purity_score_single_match() {
//         let seq = b"C";
//         let base_quals = &vec![50; 10];
//         let motif = b"C";
//         assert_eq!(
//             purity_score(seq, base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             TEST_PURITY_SCORE_PARAMS.match_weight
//         );
//     }

//     #[test]
//     fn test_purity_score_single_low_qual_match() {
//         let seq = b"C";
//         let base_quals = &vec![5; 10];
//         let motif = b"C";
//         assert_eq!(
//             purity_score(seq, base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             TEST_PURITY_SCORE_PARAMS.match_weight
//         );
//     }

//     #[test]
//     fn test_purity_score_single_high_qual_mismatch() {
//         let seq = b"A";
//         let base_quals = &vec![50; 10];
//         let motif = b"T";
//         assert_eq!(
//             purity_score(seq, base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             TEST_PURITY_SCORE_PARAMS.mismatch_weight
//         );
//     }

//     #[test]
//     fn test_purity_score_single_low_qual_mismatch() {
//         let seq = b"A";
//         let base_quals = &vec![5; 10];
//         let motif = b"T";
//         assert_eq!(
//             purity_score(seq, base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             TEST_PURITY_SCORE_PARAMS.low_qual_mismatch_weight
//         );
//     }

//     #[test]
//     fn test_purity_score_multiple() {
//         // Mismatch in first base, and low quality mismatch in 8th base
//         let seq = b"TAGCAGCTGC";
//         let mut base_quals = vec![50; 10];
//         base_quals[7] = 5;
//         let motif = b"CAG";
//         assert_eq!(
//             purity_score(seq, &base_quals, motif, &TEST_PURITY_SCORE_PARAMS),
//             (TEST_PURITY_SCORE_PARAMS.mismatch_weight
//                 + TEST_PURITY_SCORE_PARAMS.low_qual_mismatch_weight
//                 + 8.0 * TEST_PURITY_SCORE_PARAMS.match_weight)
//                 / 10.0
//         );
//     }

//     #[test]
//     fn test_seq_passes_purity_score_empty() {
//         let seq = b"";
//         let base_quals = &vec![];
//         let motif = b"CAG";
//         assert!(!seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_basic_1() {
//         let seq = b"CAGCAGCAGC";
//         let base_quals = &vec![50u8; 10];
//         let motif = b"CAG";
//         assert!(seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_basic_2() {
//         let seq = b"CAGCAGCAGC";
//         let base_quals = &vec![50u8; 10];
//         let motif = b"N";
//         assert!(!seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_low_qual_1() {
//         let seq = b"CAGCAGCAGC";
//         let base_quals = &vec![10u8; 10];
//         let motif = b"CAG";
//         assert!(seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_low_qual_2() {
//         // The score is 0.95
//         let seq = b"CAGCTGCAGC";
//         let base_quals = &vec![10u8; 10];
//         let motif = b"CAG";
//         assert!(seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_low_qual_3() {
//         // The score is 0.90
//         let seq = b"CATCATCAGC";
//         let base_quals = &vec![10u8; 20];
//         let motif = b"CAG";
//         assert!(seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_seq_passes_purity_score_low_qual_4() {
//         // The score is 0.85
//         let seq = b"CATCATCATC";
//         let base_quals = &vec![10u8; 20];
//         let motif = b"CAG";
//         assert!(!seq_passes_purity_score(
//             seq,
//             base_quals,
//             motif,
//             &TEST_PURITY_SCORE_PARAMS
//         ));
//     }

//     #[test]
//     fn test_in_repeat_read() {
//         let record = create_record("CAGCAGCAGC", &vec![10u8; 10]);
//         let motif = Sequence::from_str("CAG").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             true
//         );
//     }

//     #[test]
//     fn test_in_repeat_read_shift_1() {
//         let record = create_record("AGCAGCAGCA", &vec![10u8; 10]);
//         let motif = Sequence::from_str("CAG").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             true
//         );
//     }

//     #[test]
//     fn test_in_repeat_read_shift_2() {
//         let record = create_record("GCAGCAGCAG", &vec![10u8; 10]);
//         let motif = Sequence::from_str("CAG").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             true
//         );
//     }

//     #[test]
//     fn test_in_repeat_read_mismatch() {
//         let record = create_record("AAAAAAAAAA", &vec![10u8; 10]);
//         let motif = Sequence::from_str("G").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             false
//         );
//     }

//     #[test]
//     fn test_in_repeat_read_low_qual_match() {
//         let record = create_record("CAGCAGCAGC", &vec![10u8; 10]);
//         let motif = Sequence::from_str("CAG").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             true
//         );
//     }

//     #[test]
//     fn test_in_repeat_read_low_qual_match_rev_comp() {
//         let record = create_record("CTGCTGCTGC", &vec![10u8; 10]);
//         let motif = Sequence::from_str("CAG").unwrap();
//         assert_eq!(
//             is_in_repeat_read(&record, &motif, &TEST_PURITY_SCORE_PARAMS),
//             true
//         );
//     }

//     #[test]
//     fn test_read_num_pos_empty() {
//         let read = get_test_read(150, vec![], vec![], vec![]);
//         let flank_len = 150;
//         let motif_len = 3;
//         let num_repeats = 100;
//         assert_eq!(read.num_pos(flank_len, motif_len, num_repeats), 0);
//     }

//     #[test]
//     #[should_panic]
//     fn test_read_num_pos_empty_precondition_unmet() {
//         let read = get_test_read(150, vec![], vec![], vec![]);
//         let flank_len = 50; // shorter than seq len
//         let motif_len = 3;
//         let num_repeats = 1; // not enough repeats to make up for short flank
//         read.num_pos(flank_len, motif_len, num_repeats);
//     }

//     #[test]
//     fn test_read_num_pos_one() {
//         let read = Read {
//             qname: String::from("test_read"),
//             source_pos: 0,
//             seq_len: 5,
//             is_first_in_template: true,
//             kind: ReadKind::Flanking,
//             left_flank_pos: vec![(0, 0)],
//             right_flank_pos: vec![],
//             in_repeat_motif_pos: vec![],
//         };
//         assert_eq!(read.num_pos(10, 5, 5), 1);
//     }
// }
