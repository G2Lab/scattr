use std::fmt;
use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
};

use rust_htslib::bam;

use crate::{extract::LocusId, util::Utf8String};

pub type ReadPairId = String;

#[derive(PartialEq, Eq, Hash, Clone)]
pub enum ReadOrder {
    First,
    Last,
    Middle,
    Unknown,
}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum ReadKind {
    FlankingRead,
    InRepeatRead,
}

impl ReadKind {
    pub fn from_string(kind: &str) -> Option<Self> {
        match kind {
            "flanking" => Some(Self::FlankingRead),
            "in_repeat" => Some(Self::InRepeatRead),
            _ => None,
        }
    }
}

impl fmt::Display for ReadKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::FlankingRead => write!(f, "flanking"),
            Self::InRepeatRead => write!(f, "in_repeat"),
        }
    }
}

#[derive(PartialEq, Eq, Hash, Clone, Debug)]
pub enum ReadPairKind {
    FlankingReadPair,
    AnchoredInRepeatReadPair,
    FullyInRepeatReadPair,
    Unknown,
}

impl ReadPairKind {
    pub fn from_string(kind: &str) -> Option<Self> {
        match kind {
            "flanking" => Some(Self::FlankingReadPair),
            "anchored" => Some(Self::AnchoredInRepeatReadPair),
            "fully_in_repeat" => Some(Self::FullyInRepeatReadPair),
            "unknown" => Some(Self::Unknown),
            _ => None,
        }
    }
}

impl fmt::Display for ReadPairKind {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::FlankingReadPair => write!(f, "flanking"),
            Self::AnchoredInRepeatReadPair => write!(f, "anchored"),
            Self::FullyInRepeatReadPair => write!(f, "fully_in_repeat"),
            Self::Unknown => write!(f, "unknown"),
        }
    }
}

#[derive(PartialEq, Eq, Clone)]
pub struct ReadPair {
    pub id: ReadPairId,
    pub first: Option<bam::Record>,
    pub last: Option<bam::Record>,
}

impl ReadPair {
    pub fn len(&self) -> usize {
        match (&self.first, &self.last) {
            (Some(_), Some(_)) => 2,
            (Some(_), None) | (None, Some(_)) => 1,
            (None, None) => 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

pub struct RecordManager {
    read_pairs: HashMap<ReadPairId, ReadPair>,
    locus_index: HashMap<LocusId, HashSet<ReadPairId>>,
    read_pair_id_index: HashMap<ReadPairId, HashSet<LocusId>>,
    locus_read_kind: HashMap<(LocusId, ReadPairId, ReadOrder), Option<ReadKind>>,
}

impl RecordManager {
    pub fn new() -> Self {
        Self {
            read_pairs: HashMap::new(),
            locus_index: HashMap::new(),
            read_pair_id_index: HashMap::new(),
            locus_read_kind: HashMap::new(),
        }
    }

    pub fn add_read(&mut self, record: &bam::Record) -> (ReadPairId, ReadOrder) {
        // Get the read pair ID and order (first or last)
        let read_pair_id = Self::get_read_pair_id(record);
        let read_order = Self::get_record_order(record);

        // Retrieve (or insert, if it doesn't exist) the read pair
        let read_pair = self
            .read_pairs
            .entry(read_pair_id.clone())
            .or_insert(ReadPair {
                id: read_pair_id.clone(),
                first: None,
                last: None,
            });

        // Update the read pair with the new read
        match read_order {
            ReadOrder::First => {
                if read_pair.first.is_none() {
                    read_pair.first = Some(record.clone());
                }
            }
            ReadOrder::Last => {
                if read_pair.last.is_none() {
                    read_pair.last = Some(record.clone());
                }
            }
            _ => {}
        }

        // // Insert the read pair ID into the read pair ID index
        // self.read_pair_id_index
        //     .insert(read_pair_id.clone(), HashSet::new());

        (read_pair.id.clone(), read_order)
    }

    pub fn get_read_pair_by_id(&self, read_pair_id: &ReadPairId) -> Option<ReadPair> {
        self.read_pairs.get(read_pair_id).cloned()
    }

    pub fn add_read_to_locus_with_kind(
        &mut self,
        locus_id: &LocusId,
        record: &bam::Record,
        kind: ReadKind,
    ) -> (ReadPairId, ReadOrder) {
        // Add the read to the record manager
        let (read_pair_id, read_order) = self.add_read(record);

        // Associate the read pair with the locus
        self.associate_read_pair_with_locus(locus_id, &read_pair_id);

        // Set the read kind
        self.set_read_kind(locus_id, &read_pair_id, &read_order, kind);

        (read_pair_id, read_order)
    }

    fn associate_read_pair_with_locus(&mut self, locus_id: &LocusId, read_pair_id: &ReadPairId) {
        self.add_to_locus_index(locus_id, read_pair_id);
        self.add_to_read_pair_id_index(read_pair_id, locus_id);
    }

    fn add_to_read_pair_id_index(&mut self, read_pair_id: &str, locus_id: &str) {
        self.read_pair_id_index
            .entry(read_pair_id.to_owned())
            .or_default()
            .insert(locus_id.to_owned());
    }

    fn add_to_locus_index(&mut self, locus_id: &str, read_pair_id: &str) {
        self.locus_index
            .entry(locus_id.to_owned())
            .or_default()
            .insert(read_pair_id.to_owned());
    }

    pub fn iter_read_pairs(&self) -> impl Iterator<Item = &ReadPair> {
        self.read_pairs.values()
    }

    pub fn get_locus_read_pairs(&self, locus_id: &LocusId) -> Vec<ReadPair> {
        self.locus_index
            .get(locus_id)
            .cloned()
            .unwrap_or_else(HashSet::new)
            .iter()
            .map(|read_pair_id| self.read_pairs.get(read_pair_id).unwrap())
            .cloned()
            .collect()
    }

    pub fn get_read_pair_locus_ids(&self, read_pair_id: &ReadPairId) -> Vec<LocusId> {
        self.read_pair_id_index
            .get(read_pair_id)
            .cloned()
            .unwrap_or_else(HashSet::new)
            .into_iter()
            .collect()
    }

    pub fn contains_read_pair_id(&self, read_pair_id: &ReadPairId) -> bool {
        self.read_pairs.contains_key(read_pair_id)
    }

    pub fn set_read_kind(
        &mut self,
        locus_id: &LocusId,
        read_pair_id: &ReadPairId,
        read_order: &ReadOrder,
        read_kind: ReadKind,
    ) {
        self.locus_read_kind.insert(
            (locus_id.clone(), read_pair_id.clone(), read_order.clone()),
            Some(read_kind),
        );
    }

    pub fn get_read_kind(
        &self,
        locus_id: &LocusId,
        read_pair_id: &ReadPairId,
        read_order: &ReadOrder,
    ) -> Option<ReadKind> {
        self.locus_read_kind
            .get(&(locus_id.clone(), read_pair_id.clone(), read_order.clone()))
            .cloned()
            .unwrap_or(None)
    }

    pub fn get_read_pair_kind(
        &self,
        locus_id: &LocusId,
        read_pair_id: &ReadPairId,
    ) -> ReadPairKind {
        let first_kind = self.get_read_kind(locus_id, read_pair_id, &ReadOrder::First);
        let last_kind = self.get_read_kind(locus_id, read_pair_id, &ReadOrder::Last);

        match (first_kind, last_kind) {
            (Some(ReadKind::FlankingRead), Some(ReadKind::FlankingRead)) => {
                ReadPairKind::FlankingReadPair
            }
            (Some(ReadKind::InRepeatRead), Some(ReadKind::InRepeatRead)) => {
                ReadPairKind::FullyInRepeatReadPair
            }
            (Some(ReadKind::FlankingRead), Some(ReadKind::InRepeatRead))
            | (Some(ReadKind::InRepeatRead), Some(ReadKind::FlankingRead)) => {
                ReadPairKind::AnchoredInRepeatReadPair
            }
            _ => ReadPairKind::Unknown,
        }
    }

    pub fn incomplete_read_pair_ids_for_locus(&self, locus_id: &LocusId) -> Vec<ReadPairId> {
        self.get_locus_read_pairs(locus_id)
            .iter()
            .filter(|read_pair| read_pair.len() != 2)
            .map(|read_pair| read_pair.id.clone())
            .collect()
    }

    pub fn incomplete_read_pair_ids(&self) -> Vec<ReadPairId> {
        self.read_pairs
            .values()
            .filter(|read_pair| read_pair.len() != 2)
            .map(|read_pair| read_pair.id.clone())
            .collect()
    }

    pub fn get_record_order(record: &bam::Record) -> ReadOrder {
        let is_first = record.is_first_in_template();
        let is_last = record.is_last_in_template();

        match (is_first, is_last) {
            (true, false) => ReadOrder::First,    // First read in template
            (false, true) => ReadOrder::Last,     // Last read in template
            (true, true) => ReadOrder::Middle,    // Both 0x40 and 0x80 are set
            (false, false) => ReadOrder::Unknown, // Neither 0x40 nor 0x80 is set
        }
    }
    pub fn get_read_pair_id(record: &bam::Record) -> ReadPairId {
        record.qname().to_string()
    }
}

impl Default for RecordManager {
    fn default() -> Self {
        Self::new()
    }
}
