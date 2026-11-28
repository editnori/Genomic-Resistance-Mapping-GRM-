//! Single-pass streaming preflight estimator (`grm kcount preflight`).
//!
//! Uses an `ntHash`-style O(1) rolling canonical 64-bit hash feeding a
//! HyperLogLog (`p = 14`, 16,384 1-byte registers = 16 KiB) plus an adaptive
//! hash-prefix subsampled exact frequency map capped at 65,536 `(u64, u32)`
//! entries (~2.1 MiB table capacity, distinct from process RSS) to estimate
//! `F0`, `f1`, `f2+`, peak coverage, GC fraction, and expected disk footprints
//! across `.kset` v1, `.grk`, and `.spss` archives without sorting.

use crate::dna::CODE;
use crate::fasta;
use crate::kcount::format::{DEFAULT_SAMPLE_NAME, InputModality};
use crate::kcount::key::{MAX_K, MIN_K};
use crate::kset::{BUCKET_BITS, BUCKETS};
use anyhow::{Result, bail, ensure};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::io::BufRead;
use std::path::{Path, PathBuf};

/// Prefix for deterministic SHA-256 derivation of the 4 odd 64-bit `ntHash`
/// base seeds (`A=0, C=1, G=2, T=3`).
#[cfg(test)]
const NTHASH_SEED_PREFIX: &str = "GRM-KCOUNT-NTHASH-SEED-V1:";

/// Four odd 64-bit base seeds for 2-bit codes `A=0, C=1, G=2, T=3`,
/// deterministically derived from `SHA-256("GRM-KCOUNT-NTHASH-SEED-V1:0")`
/// through `":3"`, taking the first 8 digest bytes as little-endian `u64` with
/// the low bit forced to 1 (`| 1`). Verified by
/// `nthash_seeds_match_sha256_derivation`.
pub const BASE_SEED: [u64; 4] = [
    0x40b5_b010_e9b6_57bd,
    0x6814_898c_c720_32c5,
    0xe578_74c7_86b0_3335,
    0x93bd_87dd_ed20_45d3,
];

/// First public-domain SplitMix64 finalizer multiplier from Steele, Lea, and
/// Flood (OOPSLA 2014, public domain).
const SPLITMIX64_MUL_1: u64 = 0xbf58_476d_1ce4_e5b9;

/// Second public-domain SplitMix64 finalizer multiplier from Steele, Lea, and
/// Flood (OOPSLA 2014, public domain).
const SPLITMIX64_MUL_2: u64 = 0x94d0_49bb_1331_11eb;

#[inline]
fn rc_seed(code: u8) -> u64 {
    BASE_SEED[(3 - code) as usize]
}

#[inline]
fn final_mix64(mut x: u64) -> u64 {
    x ^= x >> 30;
    x = x.wrapping_mul(SPLITMIX64_MUL_1);
    x ^= x >> 27;
    x = x.wrapping_mul(SPLITMIX64_MUL_2);
    x ^ (x >> 31)
}

/// Reusable streaming reader over plain or `.gz` FASTA/FASTQ preserving FASTQ
/// Phred quality bytes.
struct SeqQualReader {
    reader: Box<dyn BufRead + Send>,
    path_display: String,
    line_buf: Vec<u8>,
    seq_buf: Vec<u8>,
    plus_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    in_fasta: bool,
    saw_fastq: bool,
    eof: bool,
}

impl SeqQualReader {
    fn open(path: &Path) -> Result<Self> {
        let reader = fasta::open(path)?;
        Ok(Self {
            reader,
            path_display: path.display().to_string(),
            line_buf: Vec::with_capacity(512),
            seq_buf: Vec::with_capacity(1 << 16),
            plus_buf: Vec::with_capacity(16),
            qual_buf: Vec::with_capacity(512),
            in_fasta: false,
            saw_fastq: false,
            eof: false,
        })
    }

    fn saw_fastq(&self) -> bool {
        self.saw_fastq
    }

    fn for_each_record(
        &mut self,
        mut f: impl FnMut(&[u8], Option<&[u8]>) -> Result<()>,
    ) -> Result<()> {
        while !self.eof {
            self.line_buf.clear();
            if self.reader.read_until(b'\n', &mut self.line_buf)? == 0 {
                self.eof = true;
                break;
            }
            trim_newline(&mut self.line_buf);
            match self.line_buf.first() {
                Some(b'>') => {
                    if self.in_fasta {
                        f(&self.seq_buf, None)?;
                    }
                    self.seq_buf.clear();
                    self.in_fasta = true;
                }
                Some(b'@') if !self.in_fasta || self.seq_buf.is_empty() => {
                    self.saw_fastq = true;

                    self.seq_buf.clear();
                    let n_seq = self.reader.read_until(b'\n', &mut self.seq_buf)?;
                    ensure!(
                        n_seq > 0,
                        "{}: truncated FASTQ sequence line",
                        self.path_display
                    );
                    trim_newline(&mut self.seq_buf);

                    self.plus_buf.clear();
                    let n_plus = self.reader.read_until(b'\n', &mut self.plus_buf)?;
                    ensure!(
                        n_plus > 0 && self.plus_buf.first() == Some(&b'+'),
                        "{}: malformed FASTQ separator line",
                        self.path_display
                    );

                    self.qual_buf.clear();
                    let n_qual = self.reader.read_until(b'\n', &mut self.qual_buf)?;
                    ensure!(
                        n_qual > 0,
                        "{}: truncated FASTQ quality line",
                        self.path_display
                    );
                    trim_newline(&mut self.qual_buf);
                    ensure!(
                        self.qual_buf.len() == self.seq_buf.len(),
                        "{}: FASTQ quality length {} differs from sequence length {}",
                        self.path_display,
                        self.qual_buf.len(),
                        self.seq_buf.len()
                    );

                    f(&self.seq_buf, Some(&self.qual_buf))?;
                    self.in_fasta = false;
                }
                Some(_) if self.in_fasta => self.seq_buf.extend_from_slice(&self.line_buf),
                Some(_) => bail!("{}: sequence data before any header", self.path_display),
                None => {}
            }
        }
        if self.in_fasta {
            f(&self.seq_buf, None)?;
            self.in_fasta = false;
        }
        Ok(())
    }
}

#[inline]
fn trim_newline(line: &mut Vec<u8>) {
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
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

    /// Streams canonical rolling hashes over `seq`, breaking windows on any
    /// non-ACGT base or FASTQ base with Phred score strictly below `min_qual`.
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
            if c > 3 || (use_qual && qslice.get(i).copied().unwrap_or(0) < qual_floor) {
                self.reset();
                continue;
            }
            if let Some(h) = self.push(c) {
                f(h);
            }
        }
    }
}

/// HyperLogLog cardinality sketch with `p = 14` (`16,384` 1-byte registers = 16 KiB).
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

/// Adaptive hash-prefix subsampler keeping at most `DEFAULT_MAX_ENTRIES`
/// (`65,536`) exact `(hash -> count)` entries (~2.1 MiB table capacity).
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

/// Configuration for [`preflight_inputs_with_options`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PreflightOptions {
    pub sample: String,
    pub k: usize,
    pub min_qual: u8,
    pub modality: InputModality,
}

impl Default for PreflightOptions {
    fn default() -> Self {
        Self {
            sample: DEFAULT_SAMPLE_NAME.to_string(),
            k: 31,
            min_qual: 0,
            modality: InputModality::Unknown,
        }
    }
}

/// Output report of [`preflight_inputs`] and [`preflight_inputs_with_options`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct PreflightReport {
    /// Explicit sample label (defaults to `"sample"`; input file paths and
    /// filename stems are never copied here).
    pub sample: String,
    pub k: usize,
    pub min_qual: u8,
    pub total_records: u64,
    pub total_bases: u64,
    pub valid_bases: u64,
    pub gc_fraction: f64,
    pub total_kmer_occurrences: u64,
    /// HyperLogLog cardinality estimate (`F0`).
    pub f0_hll_est: u64,
    /// Final `F0` estimate (exact prefix-sample count when `sample_exact == true`,
    /// otherwise `f0_hll_est`).
    pub f0_est: u64,
    /// Estimated singleton k-mers (`f1`).
    pub f1_est: u64,
    /// Estimated solid k-mers with abundance >= 2 (`f2+`).
    pub f2_plus_est: u64,
    pub singleton_fraction_est: f64,
    pub peak_coverage_est: u32,
    pub sample_exact: bool,
    pub subsample_shift: u32,
    /// Modality profile (`"assembly"` or `"reads"`) used to size
    /// `est_spss_archive_bytes`.
    pub spss_archive_profile: String,
    /// Estimated `.kset` v1 (`GRMKSET1`) size in bytes (if `k <= 32`).
    pub est_kset_v1_bytes: Option<u64>,
    /// Estimated Golomb-Rice sorted-key table (`GRMKGR1\0`) size in bytes.
    pub est_gr_kset_bytes: u64,
    /// Modality-calibrated estimated SPSS greedy simplitig archive (`GRMSPS1\0`)
    /// size in bytes.
    pub est_spss_archive_bytes: u64,
}

/// Estimates `.spss` (`GRMSPS1\0`) archive size using a modality-calibrated
/// model (contig-overlap model for assemblies vs singleton error-tip model for
/// read sets).
fn estimate_spss_bytes(
    k: usize,
    sample_len: usize,
    f0_est: u64,
    total_records: u64,
    singleton_fraction_est: f64,
    profile_reads: bool,
) -> u64 {
    // 80-byte fixed header + 8-byte lengths_byte_len + sample name bytes.
    let header_overhead = 88u64 + sample_len as u64;
    if f0_est == 0 {
        return header_overhead;
    }
    if profile_reads {
        // Solid k-mers pack near 0.255 B/k-mer (contig-like). Sequencing-error
        // singletons form isolated error-tip simplitigs, adding ~(k - 1)/4 + 1.5
        // bytes across up to k k-mers per tip.
        let s_frac = singleton_fraction_est.clamp(0.0, 1.0);
        let tip_extra_bpk = ((k as f64 - 1.0) / 4.0 + 1.5) / (k as f64);
        let bpk = 0.255 + s_frac * tip_extra_bpk;
        header_overhead + ((f0_est as f64) * bpk).ceil() as u64
    } else {
        // Multi-contig assembly: ~1 simplitig per contig adds (k - 1) bases
        // (packed at 2 bits/base) plus ~2 ULEB length bytes.
        let est_simplitigs = (total_records as f64).clamp(1.0, f0_est as f64);
        let est_bases = (f0_est as f64) + (k as f64 - 1.0) * est_simplitigs;
        header_overhead + (est_bases / 4.0).ceil() as u64 + (est_simplitigs * 2.0).ceil() as u64
    }
}

/// Single-pass streaming preflight estimator over one or more FASTA/FASTQ
/// (plain or `.gz`) files with default sample name `"sample"` and automatic
/// assembly-vs-reads archive profile selection.
pub fn preflight_inputs(inputs: &[PathBuf], k: usize, min_qual: u8) -> Result<PreflightReport> {
    let opts = PreflightOptions {
        k,
        min_qual,
        ..Default::default()
    };
    preflight_inputs_with_options(inputs, &opts)
}

/// Single-pass streaming preflight estimator with explicit [`PreflightOptions`].
pub fn preflight_inputs_with_options(
    inputs: &[PathBuf],
    opts: &PreflightOptions,
) -> Result<PreflightReport> {
    let k = opts.k;
    let min_qual = opts.min_qual;
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
    let mut saw_fastq = false;
    let qual_floor = 33u8.saturating_add(min_qual);

    for path in inputs {
        let mut reader = SeqQualReader::open(path)?;
        reader.for_each_record(|seq, qual| {
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
        saw_fastq |= reader.saw_fastq();
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

    let mean_occ = if f0_est == 0 {
        0.0
    } else {
        total_occ as f64 / f0_est as f64
    };
    let mean_record_len = if total_records == 0 {
        0.0
    } else {
        total_bases as f64 / total_records as f64
    };
    let profile_reads = match opts.modality {
        InputModality::Assembly => false,
        InputModality::Reads => true,
        InputModality::Unknown => {
            saw_fastq || (mean_occ > 1.3 && f2_plus_est > 0 && mean_record_len <= 2_000.0)
        }
    };
    let spss_archive_profile = if profile_reads { "reads" } else { "assembly" }.to_string();

    let low_bits = (2 * k as u32 - BUCKET_BITS) as f64;
    let avg_per_bucket = (f0_est as f64 / BUCKETS as f64).max(1.0);
    let gr_bits_per_kmer = (low_bits - avg_per_bucket.log2()).max(0.0) + 1.55;
    let gr_overhead = 72
        + opts.sample.len() as u64
        + (BUCKETS as u64 + 1) * 8
        + (BUCKETS as u64) * 4
        + BUCKETS as u64;
    let est_gr_kset_bytes = gr_overhead + ((f0_est as f64 * gr_bits_per_kmer) / 8.0).ceil() as u64;

    let est_kset_v1_bytes = (k <= 32).then(|| {
        let leb_bytes = ((low_bits - avg_per_bucket.log2()).max(1.0) / 7.0).ceil();
        let v1_overhead =
            24 + opts.sample.len() as u64 + (BUCKETS as u64 + 1) * 8 + (BUCKETS as u64) * 4;
        v1_overhead + (f0_est as f64 * leb_bytes).ceil() as u64
    });

    let est_spss_archive_bytes = estimate_spss_bytes(
        k,
        opts.sample.len(),
        f0_est,
        total_records,
        singleton_fraction_est,
        profile_reads,
    );

    Ok(PreflightReport {
        sample: opts.sample.clone(),
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
        spss_archive_profile,
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
    use ring::digest::{SHA256, digest};

    #[test]
    fn nthash_seeds_match_sha256_derivation() {
        for (idx, &seed) in BASE_SEED.iter().enumerate() {
            let label = format!("{NTHASH_SEED_PREFIX}{idx}");
            let h = digest(&SHA256, label.as_bytes());
            let derived = u64::from_le_bytes(h.as_ref()[..8].try_into().unwrap()) | 1;
            assert_eq!(seed, derived, "ntHash seed mismatch for {label}");
            assert_eq!(seed & 1, 1, "ntHash seed {label} must be odd");
        }
        // All four base seeds must be pairwise distinct.
        for (i, &seed_i) in BASE_SEED.iter().enumerate() {
            for &seed_j in &BASE_SEED[(i + 1)..] {
                assert_ne!(seed_i, seed_j);
            }
        }
    }

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
