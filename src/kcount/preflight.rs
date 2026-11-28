//! Single-pass streaming preflight estimator (`grm kcount preflight`).
//!
//! Uses an `ntHash`-style O(1) rolling canonical 64-bit hash feeding a
//! HyperLogLog (`p = 14`, 16 KB) register bank plus an adaptive hash-prefix
//! subsampled exact frequency map (~2 MB bounded RAM) to estimate `F0`, `f1`,
//! `f2+`, peak coverage, and expected disk footprints across `.kset` v1,
//! Golomb-Rice sorted keys, and SPSS archives without sorting.

use crate::dna::CODE;
use crate::kcount::key::{MAX_K, MIN_K};
use crate::kcount::superkmer::StreamingSeqReader;
use crate::kset::{BUCKET_BITS, BUCKETS};
use anyhow::{Result, ensure};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::path::PathBuf;

/// 64-bit random odd seeds for A, C, G, T (canonical ntHash rolling bases).
const BASE_SEED: [u64; 4] = [
    0x3c8b_fbb3_95c6_0474,
    0x3193_c185_62a0_2b4c,
    0x2032_3ed0_8257_2324,
    0x2955_49f5_4be2_4456,
];

#[inline]
fn rc_seed(code: u8) -> u64 {
    BASE_SEED[(3 - code) as usize]
}

#[inline]
fn final_mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// O(1) per-base strand-invariant `ntHash`-style rolling canonical 64-bit hash
/// supporting any `k` in `5..=64` with zero heap allocations.
pub struct NtHash64 {
    k: usize,
    k_u32: u32,
    km1_u32: u32,
    ring: [u8; MAX_K],
    filled: usize,
    pos: usize,
    fwd: u64,
    rev: u64,
}

impl NtHash64 {
    pub fn new(k: usize) -> Self {
        assert!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
        Self {
            k,
            k_u32: k as u32,
            km1_u32: (k - 1) as u32,
            ring: [0u8; MAX_K],
            filled: 0,
            pos: 0,
            fwd: 0,
            rev: 0,
        }
    }

    #[inline]
    pub fn reset(&mut self) {
        self.filled = 0;
        self.pos = 0;
        self.fwd = 0;
        self.rev = 0;
    }

    /// Feeds one 2-bit base `c` (`0..=3`), returning `Some(canonical_hash)`
    /// once the rolling window holds `k` contiguous valid bases.
    #[inline]
    pub fn push(&mut self, c: u8) -> Option<u64> {
        let h_in = BASE_SEED[c as usize];
        let rc_in = rc_seed(c);
        if self.filled < self.k {
            self.ring[self.filled] = c;
            self.fwd = self.fwd.rotate_left(1) ^ h_in;
            self.rev ^= rc_in.rotate_left(self.filled as u32);
            self.filled += 1;
            if self.filled == self.k {
                self.pos = 0;
                return Some(final_mix64(self.fwd.min(self.rev)));
            }
            None
        } else {
            let out = self.ring[self.pos];
            self.ring[self.pos] = c;
            self.pos += 1;
            if self.pos == self.k {
                self.pos = 0;
            }
            self.fwd =
                self.fwd.rotate_left(1) ^ BASE_SEED[out as usize].rotate_left(self.k_u32) ^ h_in;
            self.rev = (self.rev ^ rc_seed(out)).rotate_right(1) ^ rc_in.rotate_left(self.km1_u32);
            Some(final_mix64(self.fwd.min(self.rev)))
        }
    }

    /// Streams canonical rolling hashes over `seq`, breaking on non-ACGT or
    /// bases below `min_qual`.
    pub fn for_each_hash(
        &mut self,
        seq: &[u8],
        qual: Option<&[u8]>,
        min_qual: u8,
        mut f: impl FnMut(u64),
    ) {
        self.reset();
        let use_qual = min_qual > 0 && qual.is_some();
        let qual_floor = 33u8.saturating_add(min_qual);
        let qslice = qual.unwrap_or(&[]);
        for (i, &b) in seq.iter().enumerate() {
            let c = CODE[b as usize];
            if c > 3 || (use_qual && qslice[i] < qual_floor) {
                self.reset();
                continue;
            }
            if let Some(h) = self.push(c) {
                f(h);
            }
        }
    }
}

/// HyperLogLog cardinality sketch with `p = 14` (`16,384` 1-byte registers).
struct HyperLogLog14 {
    regs: Vec<u8>,
}

impl HyperLogLog14 {
    const P: u32 = 14;
    const M: usize = 1 << Self::P;
    const ALPHA: f64 = 0.7213 / (1.0 + 1.079 / (Self::M as f64));

    fn new() -> Self {
        Self {
            regs: vec![0u8; Self::M],
        }
    }

    #[inline]
    fn add(&mut self, hash: u64) {
        let idx = (hash >> (64 - Self::P)) as usize;
        let rem = (hash << Self::P) | (1u64 << (Self::P - 1));
        let rho = (rem.leading_zeros() + 1) as u8;
        if rho > self.regs[idx] {
            self.regs[idx] = rho;
        }
    }

    fn estimate(&self) -> f64 {
        let mut sum = 0.0f64;
        let mut zeros = 0usize;
        for &r in &self.regs {
            if r == 0 {
                zeros += 1;
            }
            sum += 2.0f64.powi(-(r as i32));
        }
        let m = Self::M as f64;
        let raw = Self::ALPHA * m * m / sum;
        if raw <= 2.5 * m && zeros > 0 {
            m * (m / zeros as f64).ln()
        } else {
            raw
        }
    }
}

/// Adaptive hash-prefix subsampler keeping at most `MAX_SAMPLE_ENTRIES`
/// exact `(hash -> count)` entries (~2 MB bounded RAM).
struct AdaptiveFreqSketch {
    shift: u32,
    max_entries: usize,
    counts: FxHashMap<u64, u32>,
}

impl AdaptiveFreqSketch {
    const DEFAULT_MAX_ENTRIES: usize = 65_536;

    fn new() -> Self {
        Self {
            shift: 0,
            max_entries: Self::DEFAULT_MAX_ENTRIES,
            counts: FxHashMap::default(),
        }
    }

    #[inline]
    fn accepts(&self, hash: u64) -> bool {
        self.shift == 0 || (hash >> (64 - self.shift)) == 0
    }

    #[inline]
    fn add(&mut self, hash: u64) {
        if !self.accepts(hash) {
            return;
        }
        let entry = self.counts.entry(hash).or_insert(0);
        *entry = entry.saturating_add(1);
        if self.counts.len() > self.max_entries && self.shift < 48 {
            self.shift += 1;
            let s = self.shift;
            self.counts.retain(|&h, _| (h >> (64 - s)) == 0);
        }
    }

    fn summary(&self) -> (u64, u64, u64, u32, bool) {
        let scale = 1u64 << self.shift;
        let f0 = (self.counts.len() as u64).saturating_mul(scale);
        let mut f1 = 0u64;
        let mut f2_plus = 0u64;
        let mut hist = FxHashMap::<u32, u64>::default();
        for &c in self.counts.values() {
            if c == 1 {
                f1 += 1;
            } else {
                f2_plus += 1;
                *hist.entry(c).or_insert(0) += 1;
            }
        }
        let peak = hist
            .into_iter()
            .max_by_key(|&(c, n)| (n, std::cmp::Reverse(c)))
            .map(|(c, _)| c)
            .unwrap_or(if f1 > 0 { 1 } else { 0 });
        (
            f0,
            f1.saturating_mul(scale),
            f2_plus.saturating_mul(scale),
            peak,
            self.shift == 0,
        )
    }
}

/// Output report of [`preflight_inputs`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PreflightReport {
    pub k: usize,
    pub min_qual: u8,
    pub total_records: u64,
    pub total_bases: u64,
    pub valid_bases: u64,
    pub gc_fraction: f64,
    pub total_kmer_occurrences: u64,
    /// HyperLogLog cardinality estimate (`F0`).
    pub f0_hll_est: u64,
    /// Adaptive prefix-sample `F0` estimate (exact when `sample_exact == true`).
    pub f0_est: u64,
    /// Estimated singleton k-mers (`f1`).
    pub f1_est: u64,
    /// Estimated solid k-mers with abundance >= 2 (`f2+`).
    pub f2_plus_est: u64,
    pub singleton_fraction_est: f64,
    pub peak_coverage_est: u32,
    pub sample_exact: bool,
    pub subsample_shift: u32,
    /// Estimated `.kset` v1 (`GRMKSET1`) size in bytes (if `k <= 32`).
    pub est_kset_v1_bytes: Option<u64>,
    /// Estimated Golomb-Rice sorted-key table (`GRMKGR1\0`) size in bytes.
    pub est_gr_kset_bytes: u64,
    /// Estimated SPSS greedy simplitig archive (`GRMSPS1\0`) size in bytes
    /// (assuming ~0.30 B/k-mer bacterial overlap density).
    pub est_spss_archive_bytes: u64,
}

/// Single-pass streaming preflight estimator over one or more FASTA/FASTQ
/// (plain or `.gz`) files using ~2 MB of RAM without sorting.
pub fn preflight_inputs(inputs: &[PathBuf], k: usize, min_qual: u8) -> Result<PreflightReport> {
    ensure!(!inputs.is_empty(), "at least one input file is required");
    ensure!(
        (MIN_K..=MAX_K).contains(&k),
        "k must be within {MIN_K}..={MAX_K}, got {k}"
    );

    let mut nth = NtHash64::new(k);
    let mut hll = HyperLogLog14::new();
    let mut sketch = AdaptiveFreqSketch::new();

    let mut total_records = 0u64;
    let mut total_bases = 0u64;
    let mut valid_bases = 0u64;
    let mut gc_bases = 0u64;
    let mut total_occ = 0u64;
    let qual_floor = 33u8.saturating_add(min_qual);

    for path in inputs {
        let mut reader = StreamingSeqReader::open(path)?;
        reader.for_each_record(|_, seq, qual| {
            total_records += 1;
            total_bases += seq.len() as u64;
            let use_qual = min_qual > 0 && qual.is_some();
            let qslice = qual.unwrap_or(&[]);
            for (i, &b) in seq.iter().enumerate() {
                let c = CODE[b as usize];
                if c <= 3 && (!use_qual || qslice[i] >= qual_floor) {
                    valid_bases += 1;
                    if c == 1 || c == 2 {
                        gc_bases += 1;
                    }
                }
            }
            nth.for_each_hash(seq, qual, min_qual, |h| {
                total_occ += 1;
                hll.add(h);
                sketch.add(h);
            });
            Ok(())
        })?;
    }

    let f0_hll_est = hll.estimate().round() as u64;
    let (f0_sample, f1_est, f2_plus_est, peak_coverage_est, sample_exact) = sketch.summary();
    let f0_est = if sample_exact { f0_sample } else { f0_hll_est };
    let singleton_fraction_est = if f0_sample == 0 {
        0.0
    } else {
        f1_est as f64 / f0_sample as f64
    };
    let gc_fraction = if valid_bases == 0 {
        0.0
    } else {
        gc_bases as f64 / valid_bases as f64
    };

    let low_bits = (2 * k as u32 - BUCKET_BITS) as f64;
    let avg_per_bucket = (f0_est as f64 / BUCKETS as f64).max(1.0);
    let gr_bits_per_kmer = (low_bits - avg_per_bucket.log2()).max(0.0) + 1.55;
    let table_overhead = 32 + (BUCKETS as u64 + 1) * 8 + (BUCKETS as u64) * 4 + BUCKETS as u64;
    let est_gr_kset_bytes =
        table_overhead + ((f0_est as f64 * gr_bits_per_kmer) / 8.0).ceil() as u64;
    let est_kset_v1_bytes = (k <= 32).then(|| {
        let leb_bytes = ((low_bits - avg_per_bucket.log2()).max(1.0) / 7.0).ceil();
        let v1_overhead = 24 + (BUCKETS as u64 + 1) * 8 + (BUCKETS as u64) * 4;
        v1_overhead + (f0_est as f64 * leb_bytes).ceil() as u64
    });
    let est_spss_archive_bytes = 64 + ((f0_est as f64) * 0.30).ceil() as u64;

    Ok(PreflightReport {
        k,
        min_qual,
        total_records,
        total_bases,
        valid_bases,
        gc_fraction,
        total_kmer_occurrences: total_occ,
        f0_hll_est,
        f0_est,
        f1_est,
        f2_plus_est,
        singleton_fraction_est,
        peak_coverage_est,
        sample_exact,
        subsample_shift: sketch.shift,
        est_kset_v1_bytes,
        est_gr_kset_bytes,
        est_spss_archive_bytes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn nthash64_is_strictly_strand_invariant_across_k() {
        let mut rng = ChaCha8Rng::seed_from_u64(707);
        let seq: Vec<u8> = (0..1_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let rc: Vec<u8> = seq
            .iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect();
        for k in [5, 21, 31, 32, 33, 61, 64] {
            let mut n1 = NtHash64::new(k);
            let mut h_fwd = Vec::new();
            n1.for_each_hash(&seq, None, 0, |h| h_fwd.push(h));
            let mut h_rc = Vec::new();
            n1.for_each_hash(&rc, None, 0, |h| h_rc.push(h));
            h_rc.reverse();
            assert_eq!(h_fwd, h_rc, "strand mismatch at k={k}");
        }
    }
}
