//! Minimizer-bucketed super-k-mer streaming counter with bounded-memory spill
//! buffers, intra-genome Rayon parallelism, exact uncapped `u32` counting,
//! single-pass histogram/digest/sketch emission, and zero-per-read-allocation
//! FASTA/FASTQ streaming.

use crate::dna::CODE;
use crate::fasta;
use crate::fnv::Fnv;
use crate::kcount::format::{
    SequenceProfile, encode_bucket_counts, encode_bucket_keys, write_gr_kset_from_buckets,
    write_ktab_from_buckets, write_profiles,
};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{BUCKET_BITS, BUCKETS};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use ring::digest::{Context as Sha256Context, SHA256};
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Maximum bases per packed 2-bit super-k-mer (fits in a single `u8` length).
const MAX_SUPERKMER_LEN: usize = 255;

/// Configuration for [`count_inputs`].
#[derive(Clone, Debug)]
pub struct CountConfig {
    pub k: usize,
    pub sample_name: String,
    pub min_count: u32,
    pub max_count: Option<u32>,
    pub min_qual: u8,
    pub sketch_size: usize,
    pub max_memory_mb: usize,
    pub partitions: usize,
    pub minimizer_len: Option<usize>,
    pub presence_only_out: bool,
    pub out_path: PathBuf,
    pub gr_kset_out: Option<PathBuf>,
    pub hist_out: Option<PathBuf>,
    pub meta_out: Option<PathBuf>,
    pub profile_out: Option<PathBuf>,
    pub tmp_dir: Option<PathBuf>,
}

impl Default for CountConfig {
    fn default() -> Self {
        Self {
            k: 31,
            sample_name: "sample".to_string(),
            min_count: 1,
            max_count: None,
            min_qual: 0,
            sketch_size: 1000,
            max_memory_mb: 512,
            partitions: 256,
            minimizer_len: None,
            presence_only_out: false,
            out_path: PathBuf::from("out.ktab"),
            gr_kset_out: None,
            hist_out: None,
            meta_out: None,
            profile_out: None,
            tmp_dir: None,
        }
    }
}

/// Exact read/contig-level QC metrics collected in the single streaming pass.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct ReadQcMetrics {
    pub total_reads: u64,
    pub total_bases: u64,
    pub valid_bases: u64,
    pub masked_bases: u64,
    pub gc_fraction: f64,
    pub usable_read_len_median: f64,
    pub usable_read_len_mean: f64,
    pub usable_read_len_p10: f64,
    pub usable_read_len_p90: f64,
    pub n50: u64,
    pub is_fastq: bool,
    pub min_qual_applied: u8,
    pub singleton_fraction: f64,
}

/// Exact k-mer abundance histogram computed from all observed k-mers prior to
/// `--min-count` filtering, plus retained counts.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct AbundanceHistogram {
    pub k: usize,
    /// Total distinct canonical k-mers observed (F0, before `min_count` filter).
    pub distinct_total: u64,
    /// Distinct canonical k-mers retained after `min_count`/`max_count` filter.
    pub distinct_retained: u64,
    /// Total canonical k-mer occurrences observed across all inputs.
    pub total_occurrences: u64,
    /// Exact singleton count (f1).
    pub f1: u64,
    /// Exact count of k-mers occurring >= 2 times (f2+).
    pub f2_plus: u64,
    /// Maximum exact k-mer count observed across the entire input.
    pub max_count: u32,
    /// Peak coverage estimate (mode abundance for counts >= 2, or 1 if only singletons).
    pub peak_coverage: u32,
    /// Exact sorted `(abundance, n_distinct_kmers)` histogram bins.
    pub bins: Vec<(u32, u64)>,
}

/// Full deterministic summary returned by [`count_inputs`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct CountSummary {
    pub sample: String,
    pub k: usize,
    pub key_width_bits: u32,
    pub min_count: u32,
    pub max_count: Option<u32>,
    pub distinct_retained: u64,
    pub max_count_retained: u32,
    pub presence_only: bool,
    pub spilled_to_disk: bool,
    pub partitions: usize,
    pub minimizer_len: usize,
    /// FNV-1a64 over sorted retained mixed keys (matches `kset::digest` when `k <= 32`).
    pub fnv64: String,
    /// SHA-256 over sorted retained mixed keys (little-endian `u64` or `u128`).
    pub keys_sha256: String,
    /// SHA-256 over sorted retained `(key_le, count_le)` pairs.
    pub counts_sha256: String,
    /// Bottom-`s` MinHash sketch of retained mixed keys as lowercase hex strings.
    pub sketch_hex: Vec<String>,
    pub histogram: AbundanceHistogram,
    pub qc: ReadQcMetrics,
}

/// Zero-per-read-allocation streaming reader over plain or `.gz` FASTA/FASTQ.
pub struct StreamingSeqReader {
    reader: Box<dyn BufRead + Send>,
    path_display: String,
    line_buf: Vec<u8>,
    header_buf: Vec<u8>,
    seq_buf: Vec<u8>,
    plus_buf: Vec<u8>,
    qual_buf: Vec<u8>,
    in_fasta: bool,
    saw_fastq: bool,
    eof: bool,
}

impl StreamingSeqReader {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = fasta::open(path)?;
        Ok(Self {
            reader,
            path_display: path.display().to_string(),
            line_buf: Vec::with_capacity(512),
            header_buf: Vec::with_capacity(256),
            seq_buf: Vec::with_capacity(1 << 16),
            plus_buf: Vec::with_capacity(16),
            qual_buf: Vec::with_capacity(512),
            in_fasta: false,
            saw_fastq: false,
            eof: false,
        })
    }

    pub fn saw_fastq(&self) -> bool {
        self.saw_fastq
    }

    /// Calls `f(name_bytes, seq_bytes, optional_qual_bytes)` for every record
    /// without allocating per read once buffers reach steady-state size.
    pub fn for_each_record(
        &mut self,
        mut f: impl FnMut(&[u8], &[u8], Option<&[u8]>) -> Result<()>,
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
                        f(&self.header_buf, &self.seq_buf, None)?;
                    }
                    extract_header_id(&self.line_buf[1..], &mut self.header_buf);
                    self.seq_buf.clear();
                    self.in_fasta = true;
                }
                Some(b'@') if !self.in_fasta || self.seq_buf.is_empty() => {
                    self.saw_fastq = true;
                    extract_header_id(&self.line_buf[1..], &mut self.header_buf);

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

                    f(&self.header_buf, &self.seq_buf, Some(&self.qual_buf))?;
                    self.in_fasta = false;
                }
                Some(_) if self.in_fasta => self.seq_buf.extend_from_slice(&self.line_buf),
                Some(_) => bail!("{}: sequence data before any header", self.path_display),
                None => {}
            }
        }
        if self.in_fasta {
            f(&self.header_buf, &self.seq_buf, None)?;
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

#[inline]
fn extract_header_id(raw: &[u8], out: &mut Vec<u8>) {
    out.clear();
    let end = raw
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(raw.len());
    out.extend_from_slice(&raw[..end]);
}

/// RAII temporary spill directory removed cleanly on drop.
struct SpillDir {
    path: PathBuf,
}

impl SpillDir {
    fn create(base: Option<&Path>) -> Result<Self> {
        let parent = base.map(PathBuf::from).unwrap_or_else(std::env::temp_dir);
        fs::create_dir_all(&parent)
            .with_context(|| format!("create spill parent {}", parent.display()))?;
        let pid = std::process::id();
        for _ in 0..1024 {
            let n = SPILL_COUNTER.fetch_add(1, Ordering::Relaxed);
            let candidate = parent.join(format!("grm-kcount-{pid}-{n}"));
            match fs::create_dir(&candidate) {
                Ok(()) => return Ok(Self { path: candidate }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => {
                    return Err(e)
                        .with_context(|| format!("create spill dir {}", candidate.display()));
                }
            }
        }
        bail!("could not allocate unique spill directory")
    }

    fn partition_path(&self, p: usize) -> PathBuf {
        self.path.join(format!("p{p:04}.skm"))
    }
}

impl Drop for SpillDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

/// Deterministic splitmix64 mixer for strand-invariant minimizer ordering.
#[inline]
fn mix_minimizer(mut x: u64) -> u64 {
    x ^= 0x5a5a_5a5a_5a5a_5a5a;
    x ^= x >> 30;
    x = x.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x ^= x >> 27;
    x = x.wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

#[inline]
pub fn default_minimizer_len(k: usize) -> usize {
    if k >= 25 {
        13
    } else if k >= 15 {
        9
    } else {
        k.min(5)
    }
}

/// Bounded-memory super-k-mer partitioner with monotonic-deque minimizers.
struct PartitionCollector {
    k: usize,
    m: usize,
    m_mask: u64,
    m_shift: u32,
    p_mask: usize,
    buffers: Vec<Vec<u8>>,
    buffered_bytes: usize,
    budget_bytes: usize,
    spill_dir: Option<SpillDir>,
    spilled_flags: Vec<bool>,
    tmp_base: Option<PathBuf>,
    run_codes: Vec<u8>,
    min_hashes: Vec<u64>,
    deque: VecDeque<(usize, u64)>,
}

impl PartitionCollector {
    fn new(
        k: usize,
        m: usize,
        partitions: usize,
        max_memory_mb: usize,
        tmp_base: Option<PathBuf>,
    ) -> Result<Self> {
        ensure!(
            partitions.is_power_of_two() && (16..=1024).contains(&partitions),
            "partitions must be a power of 2 in 16..=1024, got {partitions}"
        );
        ensure!(
            (1..=k.min(31)).contains(&m),
            "minimizer length m={m} invalid for k={k}"
        );
        let budget_bytes = (max_memory_mb.max(1) * 1024 * 1024) / 2;
        let m_bits = 2 * m as u32;
        let m_mask = if m_bits == 64 {
            u64::MAX
        } else {
            (1u64 << m_bits) - 1
        };
        Ok(Self {
            k,
            m,
            m_mask,
            m_shift: 2 * (m as u32 - 1),
            p_mask: partitions - 1,
            buffers: vec![Vec::new(); partitions],
            buffered_bytes: 0,
            budget_bytes,
            spill_dir: None,
            spilled_flags: vec![false; partitions],
            tmp_base,
            run_codes: Vec::with_capacity(4096),
            min_hashes: Vec::with_capacity(4096),
            deque: VecDeque::with_capacity(64),
        })
    }

    fn ingest_record(
        &mut self,
        seq: &[u8],
        qual: Option<&[u8]>,
        min_qual: u8,
    ) -> Result<(u64, u64, u64)> {
        let qual_floor = 33u8.saturating_add(min_qual);
        let use_qual = min_qual > 0 && qual.is_some();
        let qslice = qual.unwrap_or(&[]);

        let mut valid_bases = 0u64;
        let mut gc_bases = 0u64;
        let mut masked_bases = 0u64;

        self.run_codes.clear();
        for (idx, &b) in seq.iter().enumerate() {
            let c = CODE[b as usize];
            let low_q = use_qual && qslice[idx] < qual_floor;
            if c > 3 || low_q {
                masked_bases += 1;
                if self.run_codes.len() >= self.k {
                    self.emit_run()?;
                }
                self.run_codes.clear();
            } else {
                valid_bases += 1;
                if c == 1 || c == 2 {
                    gc_bases += 1;
                }
                self.run_codes.push(c);
            }
        }
        if self.run_codes.len() >= self.k {
            self.emit_run()?;
        }
        self.run_codes.clear();
        Ok((valid_bases, gc_bases, masked_bases))
    }

    fn emit_run(&mut self) -> Result<()> {
        let len = self.run_codes.len();
        let k = self.k;
        let m = self.m;
        debug_assert!(len >= k);

        self.min_hashes.clear();
        let mut fwd = 0u64;
        let mut rev = 0u64;
        for (i, &c_u8) in self.run_codes.iter().enumerate() {
            let c = c_u8 as u64;
            fwd = ((fwd << 2) | c) & self.m_mask;
            rev = (rev >> 2) | ((3 - c) << self.m_shift);
            if i + 1 >= m {
                self.min_hashes.push(mix_minimizer(fwd.min(rev)));
            }
        }

        let w = k - m + 1;
        self.deque.clear();
        let mut sk_start = 0usize;
        let mut cur_min_hash = 0u64;
        let mut have_sk = false;

        for i in 0..self.min_hashes.len() {
            let h = self.min_hashes[i];
            while let Some(&(_, back_h)) = self.deque.back() {
                if back_h >= h {
                    self.deque.pop_back();
                } else {
                    break;
                }
            }
            self.deque.push_back((i, h));
            if i + 1 < w {
                continue;
            }
            let pos = i + 1 - w;
            while let Some(&(front_idx, _)) = self.deque.front() {
                if front_idx < pos {
                    self.deque.pop_front();
                } else {
                    break;
                }
            }
            let win_min = self.deque.front().unwrap().1;
            if !have_sk {
                sk_start = pos;
                cur_min_hash = win_min;
                have_sk = true;
            } else if win_min != cur_min_hash || (pos - sk_start + k) > MAX_SUPERKMER_LEN {
                let end = (pos - 1) + k;
                self.push_superkmer(cur_min_hash, sk_start, end)?;
                sk_start = pos;
                cur_min_hash = win_min;
            }
        }
        if have_sk {
            self.push_superkmer(cur_min_hash, sk_start, len)?;
        }
        Ok(())
    }

    fn push_superkmer(&mut self, min_hash: u64, start: usize, end: usize) -> Result<()> {
        let sk_len = end - start;
        debug_assert!((self.k..=MAX_SUPERKMER_LEN).contains(&sk_len));
        let part = (min_hash as usize) & self.p_mask;
        let packed_bytes = sk_len.div_ceil(4);
        let buf = &mut self.buffers[part];
        buf.reserve(1 + packed_bytes);
        buf.push(sk_len as u8);
        let slice = &self.run_codes[start..end];
        let (chunks, rem) = slice.as_chunks::<4>();
        for ch in chunks {
            let b = ch[0] | (ch[1] << 2) | (ch[2] << 4) | (ch[3] << 6);
            buf.push(b);
        }
        if !rem.is_empty() {
            let mut b = 0u8;
            for (j, &code) in rem.iter().enumerate() {
                b |= code << (2 * j);
            }
            buf.push(b);
        }
        self.buffered_bytes += 1 + packed_bytes;
        if self.buffered_bytes >= self.budget_bytes {
            self.spill_all()?;
        }
        Ok(())
    }

    fn spill_all(&mut self) -> Result<()> {
        if self.spill_dir.is_none() {
            self.spill_dir = Some(SpillDir::create(self.tmp_base.as_deref())?);
        }
        let dir = self.spill_dir.as_ref().unwrap();
        for (p, buf) in self.buffers.iter_mut().enumerate() {
            if buf.is_empty() {
                continue;
            }
            let path = dir.partition_path(p);
            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .with_context(|| format!("open spill partition {}", path.display()))?;
            file.write_all(buf)?;
            buf.clear();
            self.spilled_flags[p] = true;
        }
        self.buffered_bytes = 0;
        Ok(())
    }
}

struct PartitionResult<K: KmerKey> {
    bucket_entries: Vec<Vec<(K, u32)>>,
    hist_counts: rustc_hash::FxHashMap<u32, u64>,
    distinct_total: u64,
    total_occurrences: u64,
    max_count: u32,
}

fn process_partition<K: KmerKey>(
    space: &K::Space,
    k: usize,
    min_count: u32,
    max_count: Option<u32>,
    spill_path: Option<PathBuf>,
    mem_buf: Vec<u8>,
) -> Result<PartitionResult<K>> {
    let mut keys: Vec<K> = Vec::new();
    let mut unpack = |stream: &[u8]| -> Result<()> {
        let mut pos = 0usize;
        while pos < stream.len() {
            let sk_len = stream[pos] as usize;
            pos += 1;
            let nbytes = sk_len.div_ceil(4);
            ensure!(
                sk_len >= k && pos + nbytes <= stream.len(),
                "corrupt super-k-mer partition buffer"
            );
            let packed = &stream[pos..pos + nbytes];
            pos += nbytes;
            K::for_each_mixed_in_packed_2bit(space, packed, sk_len, |mixed| keys.push(mixed));
        }
        Ok(())
    };

    if let Some(path) = spill_path {
        let mut file_bytes = Vec::new();
        let mut f = BufReader::new(
            File::open(&path).with_context(|| format!("read spill {}", path.display()))?,
        );
        f.read_to_end(&mut file_bytes)?;
        unpack(&file_bytes)?;
        let _ = fs::remove_file(&path);
    }
    unpack(&mem_buf)?;
    drop(mem_buf);

    let total_occurrences = keys.len() as u64;
    let mut scratch = Vec::new();
    K::radix_sort(&mut keys, &mut scratch);
    drop(scratch);

    let mut hist_counts = rustc_hash::FxHashMap::<u32, u64>::default();
    let mut bucket_entries: Vec<Vec<(K, u32)>> = (0..BUCKETS).map(|_| Vec::new()).collect();
    let mut distinct_total = 0u64;
    let mut max_obs = 0u32;

    let mut i = 0usize;
    while i < keys.len() {
        let key = keys[i];
        let mut j = i + 1;
        while j < keys.len() && keys[j] == key {
            j += 1;
        }
        let count = u32::try_from(j - i).unwrap_or(u32::MAX);
        distinct_total += 1;
        if count > max_obs {
            max_obs = count;
        }
        *hist_counts.entry(count).or_insert(0) += 1;

        let keep = count >= min_count && max_count.is_none_or(|mc| count <= mc);
        if keep {
            let b = key.bucket_of(k, BUCKET_BITS);
            bucket_entries[b].push((key, count));
        }
        i = j;
    }

    Ok(PartitionResult {
        bucket_entries,
        hist_counts,
        distinct_total,
        total_occurrences,
        max_count: max_obs,
    })
}

/// Runs the next-generation k-mer counting engine over one or more FASTA/FASTQ
/// (plain or `.gz`) files, dispatching to the fast `u64` path (`k <= 32`) or
/// full `u128` path (`k = 33..=64`).
pub fn count_inputs(inputs: &[PathBuf], cfg: &CountConfig) -> Result<CountSummary> {
    ensure!(!inputs.is_empty(), "at least one input file is required");
    ensure!(
        (MIN_K..=MAX_K).contains(&cfg.k),
        "k must be within {MIN_K}..={MAX_K}, got {}",
        cfg.k
    );
    ensure!(cfg.min_count >= 1, "min_count must be >= 1");
    if let Some(mc) = cfg.max_count {
        ensure!(
            mc >= cfg.min_count,
            "max_count ({mc}) must be >= min_count ({})",
            cfg.min_count
        );
    }
    if cfg.k <= MAX_K_U64 {
        count_inputs_generic::<u64>(inputs, cfg)
    } else {
        count_inputs_generic::<u128>(inputs, cfg)
    }
}

fn count_inputs_generic<K: KmerKey>(inputs: &[PathBuf], cfg: &CountConfig) -> Result<CountSummary> {
    let m = cfg
        .minimizer_len
        .unwrap_or_else(|| default_minimizer_len(cfg.k));
    let mut collector = PartitionCollector::new(
        cfg.k,
        m,
        cfg.partitions,
        cfg.max_memory_mb,
        cfg.tmp_dir.clone(),
    )?;

    let mut total_reads = 0u64;
    let mut total_bases = 0u64;
    let mut valid_bases = 0u64;
    let mut gc_bases = 0u64;
    let mut masked_bases = 0u64;
    let mut is_fastq = false;
    let mut usable_lens: Vec<u64> = Vec::new();
    let mut raw_lens: Vec<u64> = Vec::new();

    for path in inputs {
        let mut reader = StreamingSeqReader::open(path)?;
        reader.for_each_record(|_, seq, qual| {
            total_reads += 1;
            let raw_len = seq.len() as u64;
            total_bases += raw_len;
            raw_lens.push(raw_len);
            let (v, gc, masked) = collector.ingest_record(seq, qual, cfg.min_qual)?;
            valid_bases += v;
            gc_bases += gc;
            masked_bases += masked;
            usable_lens.push(v);
            Ok(())
        })?;
        is_fastq |= reader.saw_fastq();
    }

    let spilled_to_disk = collector.spill_dir.is_some();
    let spill_dir = collector.spill_dir.take();
    let spilled_flags = std::mem::take(&mut collector.spilled_flags);
    let buffers = std::mem::take(&mut collector.buffers);
    let space = K::new_space(cfg.k);

    let part_results: Vec<PartitionResult<K>> = buffers
        .into_par_iter()
        .enumerate()
        .map(|(p, mem_buf)| {
            let spill_path = if spilled_flags[p] {
                spill_dir.as_ref().map(|d| d.partition_path(p))
            } else {
                None
            };
            process_partition(
                &space,
                cfg.k,
                cfg.min_count,
                cfg.max_count,
                spill_path,
                mem_buf,
            )
        })
        .collect::<Result<Vec<_>>>()?;

    drop(spill_dir);

    let mut combined_hist = rustc_hash::FxHashMap::<u32, u64>::default();
    let mut distinct_total = 0u64;
    let mut total_occurrences = 0u64;
    let mut max_obs = 0u32;
    for pr in &part_results {
        distinct_total += pr.distinct_total;
        total_occurrences += pr.total_occurrences;
        if pr.max_count > max_obs {
            max_obs = pr.max_count;
        }
        for (&c, &n) in &pr.hist_counts {
            *combined_hist.entry(c).or_insert(0) += n;
        }
    }
    let mut hist_bins: Vec<(u32, u64)> = combined_hist.into_iter().collect();
    hist_bins.sort_unstable_by_key(|&(c, _)| c);
    let f1 = hist_bins
        .iter()
        .find(|&&(c, _)| c == 1)
        .map(|&(_, n)| n)
        .unwrap_or(0);
    let f2_plus = distinct_total.saturating_sub(f1);
    let peak_coverage = hist_bins
        .iter()
        .filter(|&&(c, _)| c >= 2)
        .max_by_key(|&&(c, n)| (n, std::cmp::Reverse(c)))
        .map(|&(c, _)| c)
        .unwrap_or(if f1 > 0 { 1 } else { 0 });

    let mut buckets_merged: Vec<Vec<(K, u32)>> = (0..BUCKETS).map(|_| Vec::new()).collect();
    for pr in part_results {
        for (b, mut chunk) in pr.bucket_entries.into_iter().enumerate() {
            if !chunk.is_empty() {
                buckets_merged[b].append(&mut chunk);
            }
        }
    }

    struct EncodedBucket {
        count: u32,
        max_count: u32,
        all_one: bool,
        keys: Vec<u8>,
        counts: Vec<u8>,
    }

    let bucket_encoded: Vec<EncodedBucket> = buckets_merged
        .par_iter_mut()
        .enumerate()
        .map(|(b, entries)| {
            entries.sort_unstable_by_key(|&(key, _)| key);
            debug_assert!(entries.windows(2).all(|w| w[0].0 < w[1].0));
            let n = entries.len() as u32;
            let mut max_c = 0u32;
            let mut all_one = true;
            let mut keys = Vec::with_capacity(entries.len());
            let mut counts = Vec::with_capacity(entries.len());
            for &(key, c) in entries.iter() {
                keys.push(key);
                counts.push(c);
                if c > max_c {
                    max_c = c;
                }
                if c != 1 {
                    all_one = false;
                }
            }
            let mut key_payload = Vec::new();
            encode_bucket_keys(cfg.k, b, &keys, &mut key_payload);
            let mut count_payload = Vec::new();
            if !counts.is_empty() {
                encode_bucket_counts(&counts, &mut count_payload);
            }
            EncodedBucket {
                count: n,
                max_count: max_c,
                all_one,
                keys: key_payload,
                counts: count_payload,
            }
        })
        .collect();

    let mut fnv = Fnv::new();
    let mut keys_sha = Sha256Context::new(&SHA256);
    let mut counts_sha = Sha256Context::new(&SHA256);
    let mut sketch_hex = Vec::with_capacity(cfg.sketch_size);
    for bucket in &buckets_merged {
        for &(key, count) in bucket {
            key.with_le_bytes(|b| {
                fnv.bytes(b);
                keys_sha.update(b);
                counts_sha.update(b);
            });
            counts_sha.update(&count.to_le_bytes());
            if sketch_hex.len() < cfg.sketch_size {
                sketch_hex.push(key.to_hex());
            }
        }
    }

    let mut bucket_counts = [0u32; BUCKETS];
    let mut key_payloads = Vec::with_capacity(BUCKETS);
    let mut count_payloads = Vec::with_capacity(BUCKETS);
    let mut distinct_retained = 0u64;
    let mut max_count_retained = 0u32;
    let mut presence_only = true;
    for (b, eb) in bucket_encoded.into_iter().enumerate() {
        bucket_counts[b] = eb.count;
        distinct_retained += eb.count as u64;
        if eb.max_count > max_count_retained {
            max_count_retained = eb.max_count;
        }
        if eb.count > 0 && !eb.all_one {
            presence_only = false;
        }
        key_payloads.push(eb.keys);
        count_payloads.push(eb.counts);
    }

    if let Some(parent) = cfg.out_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    if cfg.presence_only_out {
        write_gr_kset_from_buckets(
            &cfg.out_path,
            &cfg.sample_name,
            cfg.k,
            &bucket_counts,
            &key_payloads,
        )?;
    } else {
        let count_arg = if presence_only {
            vec![Vec::new(); BUCKETS]
        } else {
            count_payloads
        };
        write_ktab_from_buckets(
            &cfg.out_path,
            &cfg.sample_name,
            cfg.k,
            cfg.min_count,
            max_count_retained,
            presence_only,
            &bucket_counts,
            &key_payloads,
            &count_arg,
        )?;
    }

    if let Some(gr_path) = &cfg.gr_kset_out {
        if let Some(parent) = gr_path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }
        write_gr_kset_from_buckets(
            gr_path,
            &cfg.sample_name,
            cfg.k,
            &bucket_counts,
            &key_payloads,
        )?;
    }

    if let Some(prof_path) = &cfg.profile_out {
        if let Some(parent) = prof_path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }
        let profiles = build_sequence_profiles(inputs, cfg.k, cfg.min_qual, &buckets_merged)?;
        write_profiles(prof_path, cfg.k, &profiles)?;
    }

    let qc = finalize_qc(
        total_reads,
        total_bases,
        valid_bases,
        gc_bases,
        masked_bases,
        &mut usable_lens,
        &mut raw_lens,
        is_fastq,
        cfg.min_qual,
        if distinct_total == 0 {
            0.0
        } else {
            f1 as f64 / distinct_total as f64
        },
    );

    let histogram = AbundanceHistogram {
        k: cfg.k,
        distinct_total,
        distinct_retained,
        total_occurrences,
        f1,
        f2_plus,
        max_count: max_obs,
        peak_coverage,
        bins: hist_bins,
    };

    if let Some(hist_path) = &cfg.hist_out {
        write_json_atomic(hist_path, &histogram)?;
    }

    let summary = CountSummary {
        sample: cfg.sample_name.clone(),
        k: cfg.k,
        key_width_bits: K::WIDTH_BITS,
        min_count: cfg.min_count,
        max_count: cfg.max_count,
        distinct_retained,
        max_count_retained,
        presence_only,
        spilled_to_disk,
        partitions: cfg.partitions,
        minimizer_len: m,
        fnv64: fnv.hex(),
        keys_sha256: hex_bytes(keys_sha.finish().as_ref()),
        counts_sha256: hex_bytes(counts_sha.finish().as_ref()),
        sketch_hex,
        histogram,
        qc,
    };

    if let Some(meta_path) = &cfg.meta_out {
        write_json_atomic(meta_path, &summary)?;
    }

    Ok(summary)
}

fn build_sequence_profiles<K: KmerKey>(
    inputs: &[PathBuf],
    k: usize,
    min_qual: u8,
    buckets: &[Vec<(K, u32)>],
) -> Result<Vec<SequenceProfile>> {
    let mut profiles = Vec::new();
    let qual_floor = 33u8.saturating_add(min_qual);
    let space64 = (k <= MAX_K_U64).then(|| KmerSpace64::new(k));
    let space128 = (k > MAX_K_U64).then(|| KmerSpace128::new(k));
    let shift = 2 * (k - 1);

    for path in inputs {
        let mut reader = StreamingSeqReader::open(path)?;
        reader.for_each_record(|name_bytes, seq, qual| {
            let name = String::from_utf8_lossy(name_bytes).into_owned();
            let seq_len = seq.len() as u64;
            if seq.len() < k {
                profiles.push(SequenceProfile {
                    name,
                    seq_len,
                    counts: Vec::new(),
                });
                return Ok(());
            }
            let n_pos = seq.len() - k + 1;
            let mut counts = vec![0u32; n_pos];
            let use_qual = min_qual > 0 && qual.is_some();
            let qslice = qual.unwrap_or(&[]);

            if let Some(sp) = space64 {
                let mask = sp.mask();
                let (mut fwd, mut rev, mut filled) = (0u64, 0u64, 0usize);
                for (i, &b) in seq.iter().enumerate() {
                    let c = CODE[b as usize];
                    if c > 3 || (use_qual && qslice[i] < qual_floor) {
                        filled = 0;
                        continue;
                    }
                    let c = c as u64;
                    fwd = ((fwd << 2) | c) & mask;
                    rev = (rev >> 2) | ((3 - c) << shift);
                    filled += 1;
                    if filled >= k {
                        let mixed = K::from_u128(sp.mix(fwd.min(rev)) as u128);
                        let bucket = &buckets[mixed.bucket_of(k, BUCKET_BITS)];
                        if let Ok(idx) = bucket.binary_search_by_key(&mixed, |&(key, _)| key) {
                            counts[i + 1 - k] = bucket[idx].1;
                        }
                    }
                }
            } else {
                let sp = space128.unwrap();
                let mask = sp.mask();
                let (mut fwd, mut rev, mut filled) = (0u128, 0u128, 0usize);
                for (i, &b) in seq.iter().enumerate() {
                    let c = CODE[b as usize];
                    if c > 3 || (use_qual && qslice[i] < qual_floor) {
                        filled = 0;
                        continue;
                    }
                    let c = c as u128;
                    fwd = ((fwd << 2) | c) & mask;
                    rev = (rev >> 2) | ((3 - c) << shift);
                    filled += 1;
                    if filled >= k {
                        let mixed = K::from_u128(sp.mix(fwd.min(rev)));
                        let bucket = &buckets[mixed.bucket_of(k, BUCKET_BITS)];
                        if let Ok(idx) = bucket.binary_search_by_key(&mixed, |&(key, _)| key) {
                            counts[i + 1 - k] = bucket[idx].1;
                        }
                    }
                }
            }

            profiles.push(SequenceProfile {
                name,
                seq_len,
                counts,
            });
            Ok(())
        })?;
    }
    Ok(profiles)
}

#[allow(clippy::too_many_arguments)]
fn finalize_qc(
    total_reads: u64,
    total_bases: u64,
    valid_bases: u64,
    gc_bases: u64,
    masked_bases: u64,
    usable_lens: &mut [u64],
    raw_lens: &mut [u64],
    is_fastq: bool,
    min_qual_applied: u8,
    singleton_fraction: f64,
) -> ReadQcMetrics {
    usable_lens.sort_unstable();
    let usable_read_len_median = quantile_u64(usable_lens, 0.5);
    let usable_read_len_p10 = quantile_u64(usable_lens, 0.1);
    let usable_read_len_p90 = quantile_u64(usable_lens, 0.9);
    let usable_read_len_mean = if total_reads == 0 {
        0.0
    } else {
        valid_bases as f64 / total_reads as f64
    };
    let gc_fraction = if valid_bases == 0 {
        0.0
    } else {
        gc_bases as f64 / valid_bases as f64
    };

    raw_lens.sort_unstable_by(|a, b| b.cmp(a));
    let mut running = 0u64;
    let mut n50 = 0u64;
    for &len in raw_lens.iter() {
        running += len;
        if running * 2 >= total_bases {
            n50 = len;
            break;
        }
    }

    ReadQcMetrics {
        total_reads,
        total_bases,
        valid_bases,
        masked_bases,
        gc_fraction,
        usable_read_len_median,
        usable_read_len_mean,
        usable_read_len_p10,
        usable_read_len_p90,
        n50,
        is_fastq,
        min_qual_applied,
        singleton_fraction,
    }
}

fn quantile_u64(sorted: &[u64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let n = sorted.len();
    if (q - 0.5).abs() < 1e-12 {
        if n % 2 == 1 {
            sorted[n / 2] as f64
        } else {
            (sorted[n / 2 - 1] as f64 + sorted[n / 2] as f64) / 2.0
        }
    } else {
        let idx = ((n - 1) as f64 * q).round() as usize;
        sorted[idx.min(n - 1)] as f64
    }
}

pub(crate) fn hex_bytes(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

pub(crate) fn write_json_atomic<T: Serialize>(path: &Path, val: &T) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let part = crate::confine::suffixed(path, ".part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    serde_json::to_writer_pretty(&mut w, val)?;
    w.write_all(b"\n")?;
    w.flush()?;
    drop(w);
    fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))?;
    Ok(())
}
