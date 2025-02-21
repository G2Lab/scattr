use std::{collections::HashMap, hash::Hash};

use bio::alignment::distance::simd::levenshtein;
use bio::{alignment::distance::simd::hamming, alphabets::dna::revcomp};
use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::fmt;
use std::str::FromStr;

use crate::genotype::TandemRepeatGenotype;

const ORIENTATIONS: [Orientation; 2] = [Orientation::Forward, Orientation::Reverse];

pub type Base = u8;
pub type Position = i64;
pub type HammingDistance = u64;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Haplotype {
    Haplotype1,
    Haplotype2,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Orientation {
    Forward,
    Reverse,
}

impl Orientation {
    pub fn opposite(&self) -> Orientation {
        match self {
            Orientation::Forward => Orientation::Reverse,
            Orientation::Reverse => Orientation::Forward,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AlignmentPosition {
    Mapped(Haplotype, Orientation, Position),
    Unmapped,
}

impl AlignmentPosition {
    pub fn is_mapped(&self) -> bool {
        match self {
            AlignmentPosition::Mapped(_, _, _) => true,
            AlignmentPosition::Unmapped => false,
        }
    }

    pub fn haplotype(&self) -> Option<Haplotype> {
        match self {
            AlignmentPosition::Mapped(haplotype, _, _) => Some(*haplotype),
            AlignmentPosition::Unmapped => None,
        }
    }

    pub fn orientation(&self) -> Option<Orientation> {
        match self {
            AlignmentPosition::Mapped(_, orientation, _) => Some(*orientation),
            AlignmentPosition::Unmapped => None,
        }
    }

    pub fn pos(&self) -> Option<Position> {
        match self {
            AlignmentPosition::Mapped(_, _, position) => Some(*position),
            AlignmentPosition::Unmapped => None,
        }
    }

    pub fn is_valid_mate_position(&self, mate: &AlignmentPosition) -> bool {
        match (self, mate) {
            (AlignmentPosition::Mapped(h1, o1, _), AlignmentPosition::Mapped(h2, o2, _)) => {
                h1 == h2 && o1.opposite() == *o2
            }
            _ => false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum RepeatAlignmentPosition {
    RelativeToRepeatStart(Orientation, Position),
    RelativeToRepeatEnd(Orientation, Position),
    WithinRepeat(Orientation, Position),
    Absolute(Orientation, Position),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct RepeatAlignmentPositionSetGenerator {
    query: Vec<Base>,
    left_flank: Vec<Base>,
    right_flank: Vec<Base>,
    repeat_unit: Vec<Base>,
    relative_pos: Vec<(RepeatAlignmentPosition, HammingDistance)>,
    num_lowest_edit_scores: usize,
}

impl RepeatAlignmentPositionSetGenerator {
    pub fn new(
        query: &[Base],
        repeat_unit: &[Base],
        left_flank_seq: &[Base],
        right_flank_seq: &[Base],
        num_lowest_edit_scores: usize,
        use_levenshtein: bool,
    ) -> Self {
        let relative_pos = Self::get_relative_pos(
            left_flank_seq,
            right_flank_seq,
            repeat_unit,
            query,
            num_lowest_edit_scores,
            use_levenshtein,
        );

        Self {
            left_flank: left_flank_seq.to_vec(),
            right_flank: right_flank_seq.to_vec(),
            repeat_unit: repeat_unit.to_vec(),
            relative_pos,
            query: query.to_vec(),
            num_lowest_edit_scores,
        }
    }

    pub fn from_precomputed(
        query: &[Base],
        repeat_unit: &[Base],
        left_flank_seq: &[Base],
        right_flank_seq: &[Base],
        num_lowest_edit_scores: usize,
        relative_pos: Vec<(RepeatAlignmentPosition, HammingDistance)>,
    ) -> Self {
        Self {
            left_flank: left_flank_seq.to_vec(),
            right_flank: right_flank_seq.to_vec(),
            repeat_unit: repeat_unit.to_vec(),
            relative_pos,
            query: query.to_vec(),
            num_lowest_edit_scores,
        }
    }

    pub fn relative_position_set(&self) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        self.relative_pos.clone()
    }

    pub fn generate_position_set(
        &self,
        genotype: TandemRepeatGenotype,
        left_flank_len: usize,
        right_flank_len: usize,
        use_levenshtein: bool,
    ) -> Vec<(AlignmentPosition, HammingDistance)> {
        let (hap1_genotype, hap2_genotype) = genotype.into();
        let (hap1_repeat_start, hap2_repeat_start) =
            (left_flank_len as Position, left_flank_len as Position);
        let (hap1_repeat_end, hap2_repeat_end) = (
            hap1_repeat_start + (hap1_genotype as Position * self.repeat_unit.len() as Position),
            hap2_repeat_start + (hap2_genotype as Position * self.repeat_unit.len() as Position),
        );
        let (hap1_len, hap2_len) = (
            hap1_repeat_end + right_flank_len as Position,
            hap1_repeat_end + right_flank_len as Position,
        );

        let hap1_abs_positions = self.generate_absolute_position_set(
            hap1_repeat_start,
            hap1_repeat_end,
            left_flank_len,
            right_flank_len,
            use_levenshtein,
        );
        let hap2_abs_positions = self.generate_absolute_position_set(
            hap2_repeat_start,
            hap2_repeat_end,
            left_flank_len,
            right_flank_len,
            use_levenshtein,
        );

        let hap1_positions = Self::convert_absolute(hap1_abs_positions, Haplotype::Haplotype1);
        let hap2_positions = Self::convert_absolute(hap2_abs_positions, Haplotype::Haplotype2);

        let merged_iter = hap1_positions.into_iter().chain(hap2_positions);

        // If any absolute positions are outside the haplotype length, remove them

        merged_iter
            .filter(|(pos, _)| match pos {
                AlignmentPosition::Mapped(Haplotype::Haplotype1, _, pos) => {
                    *pos >= 0 && *pos < (hap1_len - self.query.len() as Position + 1)
                }
                AlignmentPosition::Mapped(Haplotype::Haplotype2, _, pos) => {
                    *pos >= 0 && *pos < (hap2_len - self.query.len() as Position + 1)
                }
                _ => unreachable!(),
            })
            .collect()
    }

    pub fn distance_furthest_from_repeat_start(&self) -> u32 {
        let mut max_dist = 0;
        for (pos, _) in &self.relative_pos {
            if let RepeatAlignmentPosition::RelativeToRepeatStart(_, dist) = pos {
                max_dist = max_dist.max((*dist).unsigned_abs() as u32);
            }
        }
        max_dist
    }

    pub fn distance_furthest_from_repeat_end(&self) -> u32 {
        let mut max_dist = 0;
        for (pos, _) in &self.relative_pos {
            if let RepeatAlignmentPosition::RelativeToRepeatEnd(_, dist) = pos {
                let query_end = dist + self.query.len() as Position;
                max_dist = max_dist.max((query_end).unsigned_abs() as u32);
            }
        }
        max_dist
    }

    fn convert_absolute(
        positions: Vec<(RepeatAlignmentPosition, HammingDistance)>,
        haplotype: Haplotype,
    ) -> Vec<(AlignmentPosition, HammingDistance)> {
        positions
            .into_iter()
            .map(|(pos, d)| match pos {
                RepeatAlignmentPosition::Absolute(orientation, pos) => {
                    (AlignmentPosition::Mapped(haplotype, orientation, pos), d)
                }
                _ => unreachable!(),
            })
            .collect()
    }

    fn generate_absolute_position_set(
        &self,
        repeat_start: Position,
        repeat_end: Position,
        left_flank_len: usize,
        right_flank_len: usize,
        use_levenshtein: bool,
    ) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        if ((repeat_end - repeat_start) as usize) < self.query.len() {
            self.generate_absolute_position_set_short(
                repeat_start,
                repeat_end,
                left_flank_len,
                right_flank_len,
                use_levenshtein,
            )
        } else {
            self.generate_absolute_position_set_long(repeat_start, repeat_end)
        }
    }

    fn generate_absolute_position_set_short(
        &self,
        repeat_start: Position,
        repeat_end: Position,
        left_flank_len: usize,
        right_flank_len: usize,
        use_levenshtein: bool,
    ) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        assert!((repeat_end - repeat_start) % self.repeat_unit.len() as i64 == 0);
        let num_repeats = (repeat_end - repeat_start) as usize / self.repeat_unit.len();
        let trimmed_left_flank =
            &self.left_flank[self.left_flank.len() - left_flank_len..self.left_flank.len()];
        let trimmed_right_flank = &self.right_flank[0..right_flank_len];
        let reference = Self::full_reference(
            trimmed_left_flank,
            trimmed_right_flank,
            &self.repeat_unit,
            num_repeats,
        );
        let dists = Self::sliding_distance(&reference, &self.query, use_levenshtein);

        Self::filter_by_num_lowest_edit_scores(dists, self.num_lowest_edit_scores)
    }

    fn generate_absolute_position_set_long(
        &self,
        repeat_start: Position,
        repeat_end: Position,
    ) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        assert!((repeat_end - repeat_start) % self.repeat_unit.len() as i64 == 0);
        let num_repeats = (repeat_end - repeat_start) as usize / self.repeat_unit.len();
        self.relative_pos
            .iter()
            .flat_map(|(aln_pos, d)| match aln_pos {
                RepeatAlignmentPosition::RelativeToRepeatStart(orientation, pos) => {
                    vec![(
                        RepeatAlignmentPosition::Absolute(*orientation, pos + repeat_start),
                        *d,
                    )]
                }
                RepeatAlignmentPosition::RelativeToRepeatEnd(orientation, pos) => {
                    vec![(
                        RepeatAlignmentPosition::Absolute(*orientation, pos + repeat_end),
                        *d,
                    )]
                }
                RepeatAlignmentPosition::WithinRepeat(orient, pos) => {
                    let abs_pos = self.within_repeat_absolute_positions(
                        *pos,
                        *orient,
                        repeat_start,
                        repeat_end,
                        num_repeats,
                    );
                    abs_pos.into_iter().map(|p| (p, *d)).collect()
                }
                RepeatAlignmentPosition::Absolute(_, _) => unreachable!(),
            })
            .collect()
    }

    fn within_repeat_absolute_positions(
        &self,
        pos: Position,
        orient: Orientation,
        repeat_start: Position,
        repeat_end: Position,
        num_repeats: usize,
    ) -> Vec<RepeatAlignmentPosition> {
        (0..num_repeats)
            .filter_map(|i| {
                let shift = repeat_start + (i * self.repeat_unit.len()) as Position;
                let pos = pos + shift;
                if pos + self.query.len() as Position <= repeat_end {
                    Some(RepeatAlignmentPosition::Absolute(orient, pos))
                } else {
                    None
                }
            })
            .collect()
    }

    fn full_reference(
        left_flank_seq: &[Base],
        right_flank_seq: &[Base],
        repeat_unit: &[Base],
        num_repeats: usize,
    ) -> Vec<Base> {
        left_flank_seq
            .iter()
            .chain(
                repeat_unit
                    .iter()
                    .cycle()
                    .take(repeat_unit.len() * num_repeats),
            )
            .chain(right_flank_seq.iter())
            .copied()
            .collect()
    }

    fn sliding_distance(
        reference: &[Base],
        query: &[Base],
        use_levenshtein: bool,
    ) -> HashMap<HammingDistance, Vec<RepeatAlignmentPosition>> {
        let query_orientations = Self::get_query_orientations(query);
        let mut dists = HashMap::new();
        for (orientation, query) in query_orientations {
            for i in 0..(reference.len() - query.len() + 1) {
                let ref_slice = &reference[i..i + query.len()];
                let d = if use_levenshtein {
                    levenshtein(ref_slice, &query) as u64
                } else {
                    hamming(ref_slice, &query)
                };
                dists
                    .entry(d)
                    .or_insert(vec![])
                    .push(RepeatAlignmentPosition::Absolute(
                        orientation,
                        i as Position,
                    ));
            }
        }
        dists
    }

    fn filter_by_num_lowest_edit_scores(
        dists: HashMap<HammingDistance, Vec<RepeatAlignmentPosition>>,
        num_lowest_edit_scores: usize,
    ) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        // Find the lowest edit scores
        let mut distances: Vec<HammingDistance> = dists.keys().cloned().collect();
        distances.sort();
        distances.truncate(num_lowest_edit_scores);

        // Filter the positions and return a vector
        let mut filtered = vec![];
        for d in distances {
            filtered.extend(dists.get(&d).unwrap().iter().map(|pos| (*pos, d)));
        }
        filtered
    }

    fn left_flank_reference(
        left_flank_seq: &[Base],
        repeat_unit: &[Base],
        pad: usize,
    ) -> Vec<Base> {
        left_flank_seq
            .iter()
            .chain(repeat_unit.iter().cycle().take(pad))
            .copied()
            .collect()
    }

    fn right_flank_reference(
        right_flank_seq: &[Base],
        repeat_unit: &[Base],
        pad: usize,
    ) -> Vec<Base> {
        let reversed: Vec<Base> = right_flank_seq
            .iter()
            .rev()
            .chain(repeat_unit.iter().rev().cycle().take(pad))
            .cloned()
            .collect();
        reversed.iter().rev().copied().collect()
    }

    fn in_repeat_reference(repeat_unit: &[Base], len: usize) -> Vec<Base> {
        repeat_unit.iter().cycle().take(len).copied().collect()
    }

    fn get_query_orientations(query: &[Base]) -> Vec<(Orientation, Vec<Base>)> {
        ORIENTATIONS
            .iter()
            .map(|orientation| match orientation {
                Orientation::Forward => (*orientation, query.to_vec()),
                Orientation::Reverse => (*orientation, revcomp(query)),
            })
            .collect()
    }

    fn get_relative_pos(
        left_flank_seq: &[Base],
        right_flank_seq: &[Base],
        repeat_unit: &[Base],
        query: &[Base],
        num_lowest_edit_scores: usize,
        use_levenshtein: bool,
    ) -> Vec<(RepeatAlignmentPosition, HammingDistance)> {
        // Create references to calculate hamming distances against
        let left_reference =
            Self::left_flank_reference(left_flank_seq, repeat_unit, query.len() - 1);
        let right_reference =
            Self::right_flank_reference(right_flank_seq, repeat_unit, query.len() - 1);
        let in_repeat_reference =
            Self::in_repeat_reference(repeat_unit, query.len() + repeat_unit.len() - 1);

        // Calculate sliding hamming distances between references and query
        let left_flank_pos = Self::sliding_distance_left_flank(
            &left_reference,
            query,
            left_flank_seq.len(),
            use_levenshtein,
        );
        let right_flank_pos = Self::sliding_distance_right_flank(
            &right_reference,
            query,
            query.len() - 1,
            use_levenshtein,
        );
        let in_repeat_pos = Self::sliding_distance_in_repeat(
            &in_repeat_reference,
            query,
            repeat_unit.len(),
            use_levenshtein,
        );

        // Merge into a single hashmap
        let chained_iter = left_flank_pos
            .into_iter()
            .chain(right_flank_pos)
            .chain(in_repeat_pos);
        let mut merged = HashMap::new();
        for (d, pos) in chained_iter {
            merged.entry(d).or_insert(vec![]).extend(pos);
        }

        // Filter by the lowest edit scores
        Self::filter_by_num_lowest_edit_scores(merged, num_lowest_edit_scores)
    }

    fn sliding_distance_left_flank(
        reference: &[Base],
        query: &[Base],
        repeat_start: usize,
        use_levenshtein: bool,
    ) -> HashMap<HammingDistance, Vec<RepeatAlignmentPosition>> {
        let query_orientations = Self::get_query_orientations(query);
        let mut dists = HashMap::new();
        for (orientation, target_query) in query_orientations {
            for i in 0..repeat_start {
                let ref_slice = &reference[i..i + target_query.len()];
                let d = if use_levenshtein {
                    levenshtein(ref_slice, &target_query) as u64
                } else {
                    hamming(ref_slice, &target_query)
                };
                let relative_pos = i as i64 - repeat_start as i64;
                dists.entry(d).or_insert(vec![]).push(
                    RepeatAlignmentPosition::RelativeToRepeatStart(orientation, relative_pos),
                );
                // let my_query = b"ACGCCGTGAAATGTCCTTGTCCTACAATGGCAAATTGGTATCTAATGGATGGAGGCTATTGACCTGAGTCAGTTTTTGGAGGCCCCGGACTATGACCAATCCTTCAAACACAATTGCGGGATATCTGTATCTGTATCTGTATCTGTATCT";
                // if query == my_query {
                //     println!(
                //         "Relative positions: {:?}, d: {}",
                //         RepeatAlignmentPosition::RelativeToRepeatStart(orientation, relative_pos),
                //         d
                //     );
                // }
            }
        }
        dists
    }

    fn sliding_distance_right_flank(
        reference: &[Base],
        query: &[Base],
        repeat_end: usize,
        use_levenshtein: bool,
    ) -> HashMap<HammingDistance, Vec<RepeatAlignmentPosition>> {
        let query_orientations = Self::get_query_orientations(query);
        let mut dists = HashMap::new();
        for (orientation, target_query) in query_orientations {
            // First position in the right flank where the query can start
            // (where the query overalps with the repeat by 1 base)
            let range_start = repeat_end + 1 - target_query.len();
            assert!(range_start == 0);
            let range_end = reference.len() + 1 - target_query.len();
            for i in range_start..range_end {
                let ref_slice = &reference[i..i + target_query.len()];
                let d = if use_levenshtein {
                    levenshtein(ref_slice, &target_query) as u64
                } else {
                    hamming(ref_slice, &target_query)
                };
                let rel_pos = (i as i64) - (repeat_end as i64);
                dists.entry(d).or_insert(vec![]).push(
                    RepeatAlignmentPosition::RelativeToRepeatEnd(orientation, rel_pos),
                );
            }
        }
        dists
    }

    fn sliding_distance_in_repeat(
        reference: &[Base],
        query: &[Base],
        repeat_unit_len: usize,
        use_levenshtein: bool,
    ) -> HashMap<HammingDistance, Vec<RepeatAlignmentPosition>> {
        let query_orientations = Self::get_query_orientations(query);
        let mut dists = HashMap::new();
        for (orientation, query) in query_orientations {
            for i in 0..repeat_unit_len {
                let ref_slice = &reference[i..i + query.len()];
                let d = if use_levenshtein {
                    levenshtein(ref_slice, &query) as u64
                } else {
                    hamming(ref_slice, &query)
                };
                dists
                    .entry(d)
                    .or_insert(vec![])
                    .push(RepeatAlignmentPosition::WithinRepeat(
                        orientation,
                        i as Position,
                    ));
            }
        }
        dists
    }
}

impl Serialize for RepeatAlignmentPosition {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let s = match self {
            RepeatAlignmentPosition::RelativeToRepeatStart(orientation, position) => {
                format!("RelativeToRepeatStart({:?}, {:?})", orientation, position)
            }
            RepeatAlignmentPosition::RelativeToRepeatEnd(orientation, position) => {
                format!("RelativeToRepeatEnd({:?}, {:?})", orientation, position)
            }
            RepeatAlignmentPosition::WithinRepeat(orientation, position) => {
                format!("WithinRepeat({:?}, {:?})", orientation, position)
            }
            RepeatAlignmentPosition::Absolute(orientation, position) => {
                format!("Absolute({:?}, {:?})", orientation, position)
            }
        };
        serializer.serialize_str(&s)
    }
}

impl<'de> Deserialize<'de> for RepeatAlignmentPosition {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        struct RepeatAlignmentPositionVisitor;

        impl<'de> Visitor<'de> for RepeatAlignmentPositionVisitor {
            type Value = RepeatAlignmentPosition;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("a string representing a RepeatAlignmentPosition enum variant")
            }

            fn visit_str<E>(self, v: &str) -> Result<Self::Value, E>
            where
                E: de::Error,
            {
                // Parse the input string and extract enum variant, orientation and position
                let mut parts = v.splitn(2, '(');
                let variant = parts.next().unwrap();
                let args = parts.next().ok_or_else(|| E::custom("invalid format"))?;
                let args = args.trim_end_matches(')').replace(' ', "");
                let mut args_iter = args.split(',');

                let orientation_str = args_iter
                    .next()
                    .ok_or_else(|| E::custom("missing orientation"))?;
                let position_str = args_iter
                    .next()
                    .ok_or_else(|| E::custom("missing position"))?;

                let orientation = match orientation_str {
                    "Forward" => Orientation::Forward,
                    "Reverse" => Orientation::Reverse,
                    _ => return Err(E::custom("invalid orientation")),
                };

                let position: i64 = i64::from_str(position_str).map_err(E::custom)?;

                match variant {
                    "RelativeToRepeatStart" => Ok(RepeatAlignmentPosition::RelativeToRepeatStart(
                        orientation,
                        position,
                    )),
                    "RelativeToRepeatEnd" => Ok(RepeatAlignmentPosition::RelativeToRepeatEnd(
                        orientation,
                        position,
                    )),
                    "WithinRepeat" => {
                        Ok(RepeatAlignmentPosition::WithinRepeat(orientation, position))
                    }
                    "Absolute" => Ok(RepeatAlignmentPosition::Absolute(orientation, position)),
                    _ => Err(E::custom("invalid variant")),
                }
            }
        }

        deserializer.deserialize_str(RepeatAlignmentPositionVisitor)
    }
}
