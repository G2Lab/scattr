use core::f64;
use std::collections::{HashMap, HashSet};
use std::sync::{Arc, Mutex};

use anyhow::Error;
use argmin::{core::CostFunction, solver::simulatedannealing::Anneal};
use rand::distributions::{Uniform, WeightedIndex};
use rand::prelude::Distribution;
use rand::seq::SliceRandom;
use rand::Rng;
use rand_xoshiro::Xoshiro256PlusPlus;
use std::ops::Range;

use crate::distributions::{DepthDistribution, InsertSizeDistribution};
use crate::positions::{AlignmentPosition, HammingDistance, Haplotype, Position};

type InsertSize = u32;
type ReadIdx = usize;
type PairAlignmentPosition = (AlignmentPosition, AlignmentPosition);

pub const INSERT_SIZE_SD_TOLERANCE: f64 = 3.0;

pub const FIXED_READ_THRESHOLD: usize = 10000;

#[derive(Clone)]
pub struct AlignmentProblem<'a> {
    read_lens: Vec<ReadIdx>,

    position_sets: Vec<Vec<AlignmentPosition>>,

    // Map of read index -> (map of alignment position -> hamming distance)
    position_to_distance_maps: Vec<HashMap<AlignmentPosition, HammingDistance>>,

    // Pairs of read indices for free pairs (more than one possible alignment solution)
    // and fixed pairs (only one possible solution)
    unmapped_pairs: Vec<(ReadIdx, ReadIdx)>,
    fixed_pairs: Vec<(ReadIdx, ReadIdx)>,
    free_pairs: Vec<(ReadIdx, ReadIdx)>,

    free_pair_idx_distr: Option<Uniform<usize>>,

    // Map of `fixed_pairs` index -> pair position
    // fixed_pairs_position_sets: Vec<Vec<PairAlignmentPosition>>,
    fixed_pair_position: Option<Vec<PairAlignmentPosition>>,

    // Fixed contribution to the log likelihood from the fixed pairs
    ll_fixed_pair_insert_sizes: Option<f64>,
    ll_fixed_pair_hamming_distances: Option<f64>,

    // Fixed base depths vectors
    fixed_pair_depths: Option<Vec<u32>>,

    base_error_rate: f64,
    min_insert_size: u32,
    max_insert_size: u32,

    p_unmap: f64,

    depth_distr: &'a dyn DepthDistribution,
    insert_distr: &'a dyn InsertSizeDistribution,
    haplotype1_len: usize,
    haplotype2_len: usize,
    haplotype1_depth_range: Range<usize>,
    haplotype2_depth_range: Range<usize>,
    rng: Arc<Mutex<Xoshiro256PlusPlus>>,
}

impl<'a> AlignmentProblem<'a> {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        read_lens: &[usize],
        position_dist_sets: Vec<Vec<(AlignmentPosition, HammingDistance)>>,
        pairs: HashSet<(usize, usize)>,
        base_error_rate: f64,
        depth_distr: &'a dyn DepthDistribution,
        insert_distr: &'a dyn InsertSizeDistribution,
        haplotype1_len: usize,
        haplotype2_len: usize,
        haplotype1_depth_range: Range<usize>,
        haplotype2_depth_range: Range<usize>,
        insert_size_sd_tolerance: f64,
        fixed_read_threshold: usize,
        p_unmap: f64,
        rng: Xoshiro256PlusPlus,
    ) -> Self {
        let read_lens = read_lens.to_vec();

        let (position_sets, position_to_distance_maps) =
            unpack_position_dist_sets(position_dist_sets);

        let (max_insert_size, min_insert_size) =
            calculate_insert_size_bounds(insert_distr, insert_size_sd_tolerance);

        let mut pairs: Vec<(usize, usize)> = pairs
            .iter()
            .map(|&(i, j)| {
                // Let i be the read with the fewest possible positions
                if position_sets[i].len() < position_sets[j].len() {
                    (i, j)
                } else {
                    (j, i)
                }
            })
            .collect();
        pairs.sort();

        // Filter position set to guarantee that any position for a read has
        // a compatible position for its mate
        let mut filtered_position_sets =
            vec![HashSet::<AlignmentPosition>::new(); position_sets.len()];

        for &(i, j) in pairs.iter() {
            let read1_len = read_lens[i];
            let read2_len = read_lens[j];

            let read_position_set = &position_sets[i];
            let mate_position_set = &position_sets[j];

            // Iterate over the positions of i
            for read_position in read_position_set.iter() {
                let compatible_positions = filter_positions_for_mate(
                    read_position,
                    mate_position_set,
                    read1_len,
                    read2_len,
                    min_insert_size,
                    max_insert_size,
                );
                if !compatible_positions.is_empty() {
                    filtered_position_sets[i].insert(*read_position);
                    filtered_position_sets[j].extend(compatible_positions);
                }
            }
        }

        let position_sets: Vec<Vec<AlignmentPosition>> = filtered_position_sets
            .into_iter()
            .map(|set| set.into_iter().collect())
            .collect();

        // Partition pair indices into mappable and unmappable pairs indices
        let (unmapped_pair_idxs, mapped_pair_idxs): (Vec<usize>, _) =
            (0..pairs.len()).partition(|&i| {
                let (i, j) = pairs[i];
                position_sets[i].is_empty() && position_sets[j].is_empty()
            });

        // Further partition mapped pair indices into free and fixed pairs indices
        let (free_pairs_idxs, fixed_pairs_idxs): (Vec<usize>, _) =
            mapped_pair_idxs.iter().partition(|&&i| {
                let (i, j) = pairs[i];
                position_sets[i].len() > fixed_read_threshold
                    || position_sets[j].len() > fixed_read_threshold
            });

        // Collect the read indices for each of unmapped, free, and fixed pairs
        let unmapped_pairs: Vec<(ReadIdx, ReadIdx)> =
            unmapped_pair_idxs.iter().map(|i| pairs[*i]).collect();
        let free_pairs: Vec<(ReadIdx, ReadIdx)> =
            free_pairs_idxs.iter().map(|i| pairs[*i]).collect();
        let fixed_pairs: Vec<(ReadIdx, ReadIdx)> =
            fixed_pairs_idxs.iter().map(|i| pairs[*i]).collect();

        // Create a uniform distribution for sampling free pairs
        let free_pair_idx_distr = if free_pairs.is_empty() {
            // warn!("No free pairs found");
            None
        } else {
            Some(Uniform::new(0, free_pairs.len()))
        };

        debug!(
            "AlignmentProblem: haplotype lengths: {}/{}, {} reads, {} free pairs, {} fixed pairs, {} unmapped pairs, insert size tolerance: {}-{}, depth range: {}-{}/{}-{}",
            haplotype1_len,
            haplotype2_len,
            read_lens.len(),
            free_pairs.len(),
            fixed_pairs.len(),
            unmapped_pairs.len(),
            min_insert_size,
            max_insert_size,
            haplotype1_depth_range.start,
            haplotype1_depth_range.end,
            haplotype2_depth_range.start,
            haplotype2_depth_range.end
        );

        Self {
            position_sets,
            read_lens,
            position_to_distance_maps,
            unmapped_pairs,
            fixed_pairs,
            free_pairs,
            free_pair_idx_distr,
            fixed_pair_position: None,
            ll_fixed_pair_insert_sizes: None,
            ll_fixed_pair_hamming_distances: None,
            fixed_pair_depths: None,
            base_error_rate,
            min_insert_size,
            max_insert_size,
            depth_distr,
            insert_distr,
            haplotype1_len,
            haplotype2_len,
            haplotype1_depth_range,
            haplotype2_depth_range,
            p_unmap,
            rng: Arc::new(Mutex::new(rng)),
        }
    }

    pub fn init(&mut self) {
        // Sample positions for the fixed pairs
        let mut rng = self.rng.lock().unwrap();
        let fixed_pairs_position: Vec<_> = (0..self.fixed_pairs.len())
            .map(|fixed_pair_idx| self.sample_fixed_pair_position(fixed_pair_idx, &mut rng))
            .collect();
        self.fixed_pair_position = Some(fixed_pairs_position);

        // Calculate the fixed pair log likelihood contributions
        self.ll_fixed_pair_insert_sizes = Some(self.ll_fixed_pair_insert_sizes());
        self.ll_fixed_pair_hamming_distances = Some(self.ll_fixed_pair_hamming_distances());
        self.fixed_pair_depths = Some(self.get_fixed_pair_depths());

        // // Create weighted index distributions for sampling free pairs
        // let free_pairs_position_distrs = (0..self.free_pairs.len())
        //     .map(|free_pair_idx| self.free_read_pair_weighted_index(free_pair_idx))
        //     .collect();
        // self.free_pairs_position_distrs = Some(free_pairs_position_distrs);

        // // Calulcate total number of positions for free pairs
        // let total_free_positions: usize = self
        //     .free_pairs_position_sets
        //     .iter()
        //     .map(|positions| positions.len())
        //     .sum();

        // debug!("Total free pair positions: {}", total_free_positions);

        // // Find the index of free pair with the most possible positions
        // let max_free_positions_idx = self
        //     .free_pairs_position_sets
        //     .iter()
        //     .enumerate()
        //     .max_by_key(|(_, positions)| positions.len())
        //     .map(|(i, _)| i)
        //     .unwrap();

        // debug!(
        //     "Free pair with most positions: {} ({} positions)",
        //     max_free_positions_idx,
        //     self.free_pairs_position_sets[max_free_positions_idx].len()
        // );

        // debug!(
        //     "Free pair positons: {:?}",
        //     self.free_pairs_position_sets[max_free_positions_idx]
        // );
    }

    fn num_reads(&self) -> usize {
        self.read_lens.len()
    }

    fn num_free_pairs(&self) -> usize {
        self.free_pairs.len()
    }

    fn num_mapped_pairs(&self) -> usize {
        self.fixed_pairs.len() + self.free_pairs.len()
    }

    pub fn problem_size(&self) -> usize {
        self.num_free_pairs()
    }

    pub fn initial_positions(&self) -> Vec<AlignmentPosition> {
        let mut positions = vec![None; self.num_reads()];

        // Set positions for unmapped pairs
        for (read1_idx, read2_idx) in self.unmapped_pairs.iter() {
            positions[*read1_idx] = Some(AlignmentPosition::Unmapped);
            positions[*read2_idx] = Some(AlignmentPosition::Unmapped);
        }

        // Set positions for the fixed pairs
        for (pair_idx, (read1_idx, read2_idx)) in self.fixed_pairs.iter().enumerate() {
            let (read1_position, read2_position) =
                self.fixed_pair_position.as_ref().unwrap()[pair_idx];
            positions[*read1_idx] = Some(read1_position);
            positions[*read2_idx] = Some(read2_position);
        }

        // Set random positions for the free pairs
        let mut rng = self.rng.lock().unwrap();
        for (pair_idx, (read1_idx, read2_idx)) in self.free_pairs.iter().enumerate() {
            let (read1_position, read2_position) =
                self.sample_free_read_pair_position(pair_idx, &mut rng);
            positions[*read1_idx] = Some(read1_position);
            positions[*read2_idx] = Some(read2_position);
        }

        positions.into_iter().map(|pos| pos.unwrap()).collect()
    }

    fn fixed_pairs_positions(&self) -> Vec<AlignmentPosition> {
        let mut positions = vec![AlignmentPosition::Unmapped; self.num_reads()];
        for (pair_idx, (read1_idx, read2_idx)) in self.fixed_pairs.iter().enumerate() {
            let (read1_position, read2_position) =
                self.fixed_pair_position.as_ref().unwrap()[pair_idx];
            positions[*read1_idx] = read1_position;
            positions[*read2_idx] = read2_position;
        }
        positions
    }

    fn pair_log_likelihood(
        &self,
        read_idxs: (usize, usize),
        positions: PairAlignmentPosition,
    ) -> f64 {
        let insert_size = self.get_insert_size_for_pair(read_idxs, positions);
        let ll_insert_size = self.insert_distr.ln_prob(insert_size as u64);
        let mut ll_hamming_distance = 0.0;
        for &(read_idx, position) in [(read_idxs.0, positions.0), (read_idxs.1, positions.1)].iter()
        {
            let errors = *self.position_to_distance_maps[read_idx]
                .get(&position)
                .unwrap();
            let matches = self.read_lens[read_idx] as HammingDistance - errors;
            ll_hamming_distance += matches as f64 * (1.0 - self.base_error_rate).ln()
                + errors as f64 * self.base_error_rate.ln();
        }
        ll_insert_size + ll_hamming_distance
    }

    fn sample_read_pair_position(
        &self,
        read1_idx: usize,
        read2_idx: usize,
        rng: &mut Xoshiro256PlusPlus,
    ) -> PairAlignmentPosition {
        let unmap_pair = rng.gen_bool(self.p_unmap);
        if unmap_pair {
            return (AlignmentPosition::Unmapped, AlignmentPosition::Unmapped);
        }

        let (read1_len, read2_len) = (self.read_lens[read1_idx], self.read_lens[read2_idx]);

        let read1_position_set = &self.position_sets[read1_idx];
        let read2_position_set = &self.position_sets[read2_idx];

        // Sample a random position for read1
        let read1_position = read1_position_set.choose(rng).unwrap();

        // Find all positions for read2 that are compatible with the position of read1
        let compatible_positions = filter_positions_for_mate(
            read1_position,
            read2_position_set,
            read1_len,
            read2_len,
            self.min_insert_size,
            self.max_insert_size,
        );

        // Calculate the likelihood of the pair for each compatible position
        let pair_log_likelihoods = compatible_positions
            .iter()
            .map(|&read2_position| {
                let positions = (*read1_position, read2_position);
                self.pair_log_likelihood((read1_idx, read2_idx), positions)
            })
            .collect::<Vec<_>>();
        let pair_likelihoods = pair_log_likelihoods
            .iter()
            .map(|&ll| ll.exp())
            .collect::<Vec<_>>();

        let weighted_index = if pair_likelihoods.iter().sum::<f64>() == 0.0 {
            // warn!(
            //     "All position likelihoods are zero for read pair ({}, {})",
            //     read1_idx, read2_idx
            // );
            WeightedIndex::new(vec![1.0; compatible_positions.len()]).unwrap()
        } else {
            WeightedIndex::new(pair_likelihoods).unwrap()
        };
        let read2_position = compatible_positions[weighted_index.sample(rng)];

        (*read1_position, read2_position)
    }

    fn sample_fixed_pair_position(
        &self,
        fixed_pair_idx: usize,
        rng: &mut Xoshiro256PlusPlus,
    ) -> PairAlignmentPosition {
        let (read1_idx, read2_idx) = self.fixed_pairs[fixed_pair_idx];
        self.sample_read_pair_position(read1_idx, read2_idx, rng)
    }

    fn sample_free_read_pair_position(
        &self,
        free_pair_idx: usize,
        rng: &mut Xoshiro256PlusPlus,
    ) -> PairAlignmentPosition {
        let (read1_idx, read2_idx) = self.free_pairs[free_pair_idx];

        // with some probability, return unmapped
        if rng.gen_bool(0.01) {
            (AlignmentPosition::Unmapped, AlignmentPosition::Unmapped)
        } else {
            self.sample_read_pair_position(read1_idx, read2_idx, rng)
        }
    }

    fn get_insert_size_for_pair(
        &self,
        read_idxs: (usize, usize),
        positions: PairAlignmentPosition,
    ) -> InsertSize {
        match positions {
            (
                AlignmentPosition::Mapped(_, _, read1_pos),
                AlignmentPosition::Mapped(_, _, read2_pos),
            ) => {
                let read1_len = self.read_lens[read_idxs.0];
                let read2_len = self.read_lens[read_idxs.1];
                calculate_insert_size((read1_pos, read2_pos), (read1_len, read2_len))
            }
            _ => unreachable!(),
        }
    }

    fn get_insert_sizes_for_pairs(
        &self,
        pairs: &[(usize, usize)],
        positions: &[AlignmentPosition],
    ) -> Vec<InsertSize> {
        let mut insert_sizes = vec![];
        for read_idxs in pairs {
            let positions = (positions[read_idxs.0], positions[read_idxs.1]);
            if positions.0.is_mapped() || positions.1.is_mapped() {
                let insert_size = self.get_insert_size_for_pair(*read_idxs, positions);
                insert_sizes.push(insert_size);
            }
        }
        insert_sizes
    }

    pub fn get_insert_sizes(&self, positions: &[AlignmentPosition]) -> Vec<InsertSize> {
        let fixed_insert_sizes =
            self.get_insert_sizes_for_pairs(&self.fixed_pairs, &self.fixed_pairs_positions());
        let free_insert_sizes = self.get_insert_sizes_for_pairs(&self.free_pairs, positions);
        let mut insert_sizes = fixed_insert_sizes;
        insert_sizes.extend(free_insert_sizes);
        insert_sizes
    }

    fn ll_insert_sizes_for_pairs(
        &self,
        pairs: &[(usize, usize)],
        positions: &[AlignmentPosition],
    ) -> f64 {
        let insert_sizes = self.get_insert_sizes_for_pairs(pairs, positions);
        let mut ll_insert_sizes = 0.0;
        for insert_size in insert_sizes.iter() {
            ll_insert_sizes += self.insert_distr.ln_prob(*insert_size as u64);
        }
        ll_insert_sizes
    }

    fn ll_fixed_pair_insert_sizes(&self) -> f64 {
        self.ll_insert_sizes_for_pairs(&self.fixed_pairs, &self.fixed_pairs_positions())
    }

    fn ll_free_pair_insert_sizes(&self, positions: &[AlignmentPosition]) -> f64 {
        self.ll_insert_sizes_for_pairs(&self.free_pairs, positions)
    }

    fn ll_insert_sizes(&self, positions: &[AlignmentPosition]) -> f64 {
        self.ll_fixed_pair_insert_sizes.unwrap() + self.ll_free_pair_insert_sizes(positions)
    }

    fn get_read_lens_for_pairs(&self, pairs: &[(usize, usize)]) -> Vec<usize> {
        let mut read_lens = vec![];
        for (read1_idx, read2_idx) in pairs {
            read_lens.push(self.read_lens[*read1_idx]);
            read_lens.push(self.read_lens[*read2_idx]);
        }
        read_lens
    }

    fn get_hamming_distances_for_pairs(
        &self,
        pairs: &[(usize, usize)],
        positions: &[AlignmentPosition],
    ) -> Vec<HammingDistance> {
        let mut distances = vec![];
        for &(read1_idx, read2_idx) in pairs {
            let read1_position = positions[read1_idx];
            let read2_position = positions[read2_idx];

            let read1_distance = *self.position_to_distance_maps[read1_idx]
                .get(&read1_position)
                .unwrap_or(&(self.read_lens[read1_idx] as HammingDistance));
            let read2_distance = *self.position_to_distance_maps[read2_idx]
                .get(&read2_position)
                .unwrap_or(&(self.read_lens[read2_idx] as HammingDistance));

            distances.push(read1_distance);
            distances.push(read2_distance);
        }
        distances
    }

    fn ll_hamming_distances_for_pairs(
        &self,
        pairs: &[(usize, usize)],
        positions: &[AlignmentPosition],
    ) -> f64 {
        let distances = self.get_hamming_distances_for_pairs(pairs, positions);
        let read_lens = self.get_read_lens_for_pairs(pairs);
        let mut ll_hamming_distances = 0.0;
        for (distance, read_len) in distances.iter().zip(read_lens.iter()) {
            let errors = *distance;
            let matches = *read_len as HammingDistance - errors;
            ll_hamming_distances += matches as f64 * (1.0 - self.base_error_rate).ln()
                + errors as f64 * self.base_error_rate.ln()
        }
        ll_hamming_distances
    }

    fn ll_fixed_pair_hamming_distances(&self) -> f64 {
        self.ll_hamming_distances_for_pairs(&self.fixed_pairs, &self.fixed_pairs_positions())
    }

    fn ll_free_pair_hamming_distances(&self, positions: &[AlignmentPosition]) -> f64 {
        self.ll_hamming_distances_for_pairs(&self.free_pairs, positions)
    }

    fn ll_hamming_distances(&self, positions: &[AlignmentPosition]) -> f64 {
        self.ll_fixed_pair_hamming_distances.unwrap()
            + self.ll_free_pair_hamming_distances(positions)
    }

    fn calculate_depths_for_pairs(
        &self,
        pairs: &[(usize, usize)],
        positions: &[AlignmentPosition],
    ) -> Vec<u32> {
        let mut hap1_depths = vec![0; self.haplotype1_len];
        let mut hap2_depths = vec![0; self.haplotype2_len];

        for (read1_idx, read2_idx) in pairs {
            for &read_idx in [read1_idx, read2_idx] {
                if let AlignmentPosition::Mapped(hap, _, start) = positions[read_idx] {
                    let read_len = self.read_lens[read_idx];
                    let hap_depths = match hap {
                        Haplotype::Haplotype1 => &mut hap1_depths,
                        Haplotype::Haplotype2 => &mut hap2_depths,
                    };
                    let start = start as usize;
                    let end = start + read_len;
                    (start..end).for_each(|i| hap_depths[i] += 1);
                }
            }
        }

        // Trim depths to the depth range
        hap1_depths = hap1_depths[self.haplotype1_depth_range.clone()].to_vec();
        hap2_depths = hap2_depths[self.haplotype2_depth_range.clone()].to_vec();

        // Concatenate the depths
        let mut depths = hap1_depths;
        depths.extend(hap2_depths);
        depths
    }

    fn get_fixed_pair_depths(&self) -> Vec<u32> {
        self.calculate_depths_for_pairs(&self.fixed_pairs, &self.fixed_pairs_positions())
    }

    fn get_free_pair_depths(&self, positions: &[AlignmentPosition]) -> Vec<u32> {
        self.calculate_depths_for_pairs(&self.free_pairs, positions)
    }

    pub fn get_depths(&self, positions: &[AlignmentPosition]) -> Vec<u32> {
        let fixed_pair_depths = self.fixed_pair_depths.as_ref().unwrap();
        let free_pair_depths = self.get_free_pair_depths(positions);
        let depths = fixed_pair_depths
            .iter()
            .zip(free_pair_depths.iter())
            .map(|(a, b)| a + b)
            .collect();
        depths
    }

    fn ll_depths(&self, positions: &[AlignmentPosition]) -> f64 {
        let depths = self.get_depths(positions);
        let mut ll_depths = 0.0;
        for depth in depths.iter() {
            ll_depths += self.depth_distr.ln_prob(*depth as u64);
        }
        ll_depths
    }

    fn depths_len(&self) -> usize {
        self.haplotype1_depth_range.len() + self.haplotype2_depth_range.len()
    }

    fn log_likelihood(&self, positions: &[AlignmentPosition]) -> f64 {
        (self.ll_hamming_distances(positions) + self.ll_insert_sizes(positions))
            / self.num_mapped_pairs() as f64
            + self.ll_depths(positions) / self.depths_len() as f64
    }

    // fn log_likelihood(&self, positions: &[AlignmentPosition]) -> f64 {
    //     // Log likelihood of observing the read sequences give the alignment positions
    //     let mut ll_read_seqs = self.fixed_ll_seq_contrib;
    //     for &(read1_idx, read2_idx) in self.free_pairs.iter() {
    //         let read_idxs = [read1_idx, read2_idx];
    //         for &read_idx in read_idxs.iter() {
    //             let position = positions[read_idx];
    //             let read_length = self.read_lens[read_idx];
    //             let errors = *self.position_to_distance_maps[read_idx]
    //                 .get(&position)
    //                 .unwrap();
    //             ll_read_seqs += Self::log_likelihood_hamming_distance(
    //                 errors,
    //                 read_length,
    //                 self.base_error_rate,
    //             );
    //         }
    //     }
    //     if self.n_reads > 0 {
    //         ll_read_seqs /= self.n_reads as f64;
    //     }

    //     // Log likelihood of observing the insert sizes given the alignment positions
    //     let mut ll_insert_size = self.ll_fixed_pair_insert_sizes;
    //     let insert_sizes = self.get_free_pairs_insert_sizes(positions);
    //     for insert_size in insert_sizes.iter() {
    //         ll_insert_size += self.insert_distr.ln_prob(*insert_size as u64);
    //     }
    //     if !insert_sizes.is_empty() {
    //         ll_insert_size /= insert_sizes.len() as f64;
    //     }

    //     // Log likelihood of observing the depths given the alignment positions
    //     let mut ll_depth = 0.0;
    //     let depths = self.get_depths(positions);
    //     for depth in depths.iter() {
    //         ll_depth += self.depth_distr.ln_prob(*depth as u64);
    //     }
    //     if !depths.is_empty() {
    //         ll_depth /= depths.len() as f64;
    //     }

    //     ll_read_seqs + ll_insert_size + ll_depth
    // }
}

fn calculate_insert_size_bounds(
    insert_distr: &dyn InsertSizeDistribution,
    sd_tolerance: f64,
) -> (u32, u32) {
    let max_insert_size =
        (insert_distr.mean() + sd_tolerance * insert_distr.sd()).round() as InsertSize;
    let min_insert_size = (insert_distr.mean() - sd_tolerance * insert_distr.sd())
        .max(0.0)
        .round() as InsertSize;
    (max_insert_size, min_insert_size)
}

fn unpack_position_dist_sets(
    position_dist_sets: Vec<Vec<(AlignmentPosition, u64)>>,
) -> (
    Vec<Vec<AlignmentPosition>>,
    Vec<HashMap<AlignmentPosition, u64>>,
) {
    // Unpack the position distance sets
    let mut position_sets = vec![];
    let mut position_to_distance_maps = vec![];

    for position_dist_set in position_dist_sets {
        let mut position_set = vec![];
        let mut position_to_distance_map = HashMap::new();
        for (aln_pos, d) in position_dist_set {
            position_set.push(aln_pos);
            position_to_distance_map.insert(aln_pos, d);
        }
        position_sets.push(position_set);
        position_to_distance_maps.push(position_to_distance_map);
    }
    (position_sets, position_to_distance_maps)
}

fn filter_positions_for_mate(
    read_position: &AlignmentPosition,
    mate_position_set: &[AlignmentPosition],
    read_len: usize,
    mate_len: usize,
    min_insert_size: InsertSize,
    max_insert_size: InsertSize,
) -> Vec<AlignmentPosition> {
    let mut valid_mate_positions = vec![];
    if let AlignmentPosition::Mapped(_, _, read_pos) = read_position {
        for mate_position in mate_position_set {
            if let AlignmentPosition::Mapped(_, _, mate_pos) = mate_position {
                // Only allow valid mate positions with insert sizes within the tolerance
                if read_position.is_valid_mate_position(mate_position) {
                    let insert_size =
                        calculate_insert_size((*read_pos, *mate_pos), (read_len, mate_len));
                    if insert_size >= min_insert_size && insert_size <= max_insert_size {
                        valid_mate_positions.push(*mate_position);
                    }
                }
            }
        }
    }
    valid_mate_positions
}

impl CostFunction for AlignmentProblem<'_> {
    type Param = Vec<AlignmentPosition>;
    type Output = f64;

    fn cost(&self, positions: &Vec<AlignmentPosition>) -> Result<Self::Output, Error> {
        // Return the negative log likelihood as the cost
        Ok(-self.log_likelihood(positions))
    }
}

impl Anneal for AlignmentProblem<'_> {
    type Param = Vec<AlignmentPosition>;
    type Output = Vec<AlignmentPosition>;
    type Float = f64;

    fn anneal(
        &self,
        positions: &Vec<AlignmentPosition>,
        temp: f64,
    ) -> Result<Vec<AlignmentPosition>, Error> {
        let mut positions = positions.clone();
        if let Some(free_pair_idx_distr) = &self.free_pair_idx_distr {
            let mut rng = self.rng.lock().unwrap();
            // Perform modifications to a degree proportional to the current temperature `temp`.
            let num_modify = (temp * self.free_pairs.len() as f64) as usize;
            let num_modify = num_modify.clamp(1, self.free_pairs.len());
            for _ in 0..(num_modify) {
                let free_pair_idx = free_pair_idx_distr.sample(&mut *rng);
                let read_idxs = self.free_pairs[free_pair_idx];
                let new_pair_positions =
                    self.sample_free_read_pair_position(free_pair_idx, &mut rng);
                positions[read_idxs.0] = new_pair_positions.0;
                positions[read_idxs.1] = new_pair_positions.1;
            }
        }

        Ok(positions)
    }
}

fn calculate_insert_size(
    read_starts: (Position, Position),
    read_lens: (usize, usize),
) -> InsertSize {
    let (read1_start, read2_start) = read_starts;
    let (read1_len, read2_len) = read_lens;
    if read1_start <= read2_start {
        // Read1 starts before or at the same position as Read2
        let read2_end = read2_start + read2_len as Position;
        assert!(read2_end >= read1_start);
        (read2_end - read1_start) as InsertSize
    } else {
        // Read2 starts before Read1
        let read1_end = read1_start + read1_len as Position;
        assert!(read1_end >= read2_start);
        (read1_end - read2_start) as InsertSize
    }
}
