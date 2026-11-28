//! Bounded-memory parallel minimizer-bucketed super-k-mer streaming counter.
//!
//! Architecture:
//! 1. **Stage 1 (Parallel ingest and minimizer partitioning):** Reads plain or
//!    `.gz` (including multi-member gzip) FASTA/FASTQ inputs in record-aligned
//!    chunks sized to `--max-memory-mb`, extracts canonical minimizer
//!    super-k-mers across Rayon workers with amortized reusable worker buffers,
//!    collects exact read/contig QC via length-frequency histograms, and spills
//!    packed 2-bit super-k-mer partition files (`p0000.skm`..) whenever buffered
//!    bytes reach the ingest budget.
//! 2. **Stage 2 (Bounded-concurrency partition counting):** Processes minimizer
//!    partitions with concurrency and per-worker key expansion capped to
//!    `--max-memory-mb`, sub-chunking and k-way merging any skewed partition
//!    with checked `u64` counts, enforcing the `u32::MAX` (`2^32 - 1`) count cap
//!    unless `--allow-saturation` is set, accumulating exact abundance histograms,
//!    and routing retained `(mixed_key, count)` entries into 256 super-bucket
//!    shards (each covering 4 of the 1,024 output buckets).
//! 3. **Stage 3 (Streaming single-super-bucket table finalization):** Loads one
//!    super-bucket at a time, sorts its 4 output buckets, computes FNV-1a64,
//!    key/count SHA-256 digests, and bottom-`s` MinHash sketch, and streams
//!    Golomb-Rice encoded payloads directly to verified `.ktab` (`GRMKTB1\0`)
//!    and/or `.grk` (`GRMKGR1\0`) files via atomic `<path>.part` renames.
//! 4. **Optional indexed parallel profile pass (`--profile-out`):** Loads a
//!    compact sorted-key lookup index with an 18-bit sub-bucket directory and
//!    streams per-sequence RLE + zigzag-delta `.prof.bin` (`GRMPRF1\0`) records
//!    encoded in parallel chunks in exact input order.

use crate::confine;
use crate::dna::CODE;
use crate::fasta;
use crate::fnv::Fnv;
use crate::kcount::format::{
    DEFAULT_SAMPLE_NAME, GrKsetReader, InputModality, KtabReader, MAGIC_GRKSET, MAGIC_KTAB,
    MAGIC_PROF1, MAX_COUNT, TableFormat, check_output_extension, encode_bucket_counts,
    encode_bucket_keys, hex_bytes, put_uleb,
};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{BUCKET_BITS, BUCKETS};
use anyhow::{Context, Result, bail, ensure};
use rayon::prelude::*;
use ring::digest::{Context as Sha256Context, SHA256};
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap, VecDeque};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::atomic::{AtomicI32, AtomicU64, AtomicUsize, Ordering};

static SPILL_COUNTER: AtomicU64 = AtomicU64::new(0);

/// Maximum bases per packed 2-bit super-k-mer (fits in a single `u8` length).
const MAX_SUPERKMER_LEN: usize = 255;

/// Number of super-bucket shard files used in Stage 2 -> Stage 3 streaming.
/// Each super-bucket holds `BUCKETS / SUPER_BUCKETS = 4` consecutive output
/// buckets, bounding Stage 3 working memory to 1/256 of retained k-mers while
/// avoiding high file-descriptor counts.
const SUPER_BUCKETS: usize = 256;
const BUCKETS_PER_SUPER: usize = BUCKETS / SUPER_BUCKETS;

const GRK_HEADER_LEN: usize = 72;
const GRK_DIGEST_START: usize = 40;
const GRK_DIGEST_END: usize = 72;

const KTAB_HEADER_LEN: usize = 80;
const KTAB_DIGEST_START: usize = 48;
const KTAB_DIGEST_END: usize = 80;

const PROF_HEADER_LEN: usize = 72;
const PROF_DIGEST_START: usize = 40;
const PROF_DIGEST_END: usize = 72;

const FLAG_PRESENCE_ONLY: u8 = 0x01;
const FLAG_SATURATED_COUNTS: u8 = 0x02;

/// Public-domain SplitMix64 finalizer multipliers (Steele, Lea, and Flood,
/// OOPSLA 2014, public domain) used for strand-invariant minimizer ordering.
const SPLITMIX64_MUL_1: u64 = 0xbf58_476d_1ce4_e5b9;
const SPLITMIX64_MUL_2: u64 = 0x94d0_49bb_1331_11eb;

/// User-selected input modality policy for `grm kcount count`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InputTypeChoice {
    /// Detect modality from input records (FASTQ -> `reads`, FASTA -> `assembly`).
    #[default]
    Auto,
    /// Expect FASTQ read inputs; error if FASTA records are encountered.
    Reads,
    /// Expect FASTA assembly inputs; error if FASTQ records are encountered.
    Assembly,
}

/// Configuration for [`count_inputs`].
#[derive(Clone, Debug)]
pub struct CountConfig {
    pub k: usize,
    pub sample_name: String,
    pub input_type: InputTypeChoice,
    pub min_count: u32,
    pub max_count: Option<u32>,
    pub allow_saturation: bool,
    pub min_qual: u8,
    pub sketch_size: usize,
    pub max_memory_mb: usize,
    pub partitions: usize,
    pub minimizer_len: Option<usize>,
    pub presence_only_out: bool,
    pub keep_record_ids: bool,
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
            sample_name: DEFAULT_SAMPLE_NAME.to_string(),
            input_type: InputTypeChoice::Auto,
            min_count: 1,
            max_count: None,
            allow_saturation: false,
            min_qual: 0,
            sketch_size: 1000,
            max_memory_mb: 512,
            partitions: 256,
            minimizer_len: None,
            presence_only_out: false,
            keep_record_ids: false,
            out_path: PathBuf::from("out.ktab"),
            gr_kset_out: None,
            hist_out: None,
            meta_out: None,
            profile_out: None,
            tmp_dir: None,
        }
    }
}

/// Exact read/contig-level QC metrics collected in the single ingest pass.
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
/// `--min-count`/`--max-count` filtering, plus retained count totals.
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
pub struct AbundanceHistogram {
    pub k: usize,
    /// Total distinct canonical k-mers observed (F0, before filtering).
    pub distinct_total: u64,
    /// Distinct canonical k-mers retained after `min_count`/`max_count` filter.
    pub distinct_retained: u64,
    /// Total canonical k-mer occurrences observed across all inputs.
    pub total_occurrences: u64,
    /// Exact singleton count (f1).
    pub f1: u64,
    /// Exact count of k-mers occurring >= 2 times (f2+).
    pub f2_plus: u64,
    /// Maximum k-mer count observed across the entire input (capped at `2^32 - 1`).
    pub max_count: u32,
    /// Number of observed distinct k-mers whose true occurrence count exceeded
    /// `2^32 - 1` and was clamped to `u32::MAX` (requires `--allow-saturation`).
    pub saturated_total: u64,
    /// Peak coverage estimate (mode abundance for counts >= 2, smallest count
    /// breaks ties, or 1 if only singletons, or 0 if empty).
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
    pub modality: InputModality,
    pub min_count: u32,
    pub max_count: Option<u32>,
    pub distinct_retained: u64,
    pub max_count_retained: u32,
    /// Number of retained k-mers clamped at `u32::MAX` (`2^32 - 1`), matching
    /// the `.ktab` header field.
    pub saturated_counts: u64,
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

// ---------------------------------------------------------------------------
// Async-signal-safe SIGINT / SIGTERM cancellation guard & stale spill cleanup
// ---------------------------------------------------------------------------

#[cfg(unix)]
mod signal_guard {
    use super::*;

    static INTERRUPTED_SIGNAL: AtomicI32 = AtomicI32::new(0);

    struct InstallState {
        active: usize,
        prev_int: libc::sigaction,
        prev_term: libc::sigaction,
    }

    static INSTALL_LOCK: Mutex<Option<InstallState>> = Mutex::new(None);

    extern "C" fn handle_termination_signal(sig: libc::c_int) {
        INTERRUPTED_SIGNAL.store(sig, Ordering::SeqCst);
    }

    pub(super) struct SignalGuard;

    impl SignalGuard {
        pub(super) fn install() -> Self {
            let mut guard = INSTALL_LOCK
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            if let Some(state) = guard.as_mut() {
                state.active += 1;
                return Self;
            }
            INTERRUPTED_SIGNAL.store(0, Ordering::SeqCst);
            // SAFETY: Zeroed `sigaction` structs are valid initializers before
            // filling `sa_sigaction` and calling `sigemptyset`/`sigaction`.
            let (prev_int, prev_term) = unsafe {
                let mut action: libc::sigaction = std::mem::zeroed();
                action.sa_sigaction = handle_termination_signal as *const () as usize;
                action.sa_flags = libc::SA_RESTART;
                libc::sigemptyset(&mut action.sa_mask);

                let mut old_int: libc::sigaction = std::mem::zeroed();
                let mut old_term: libc::sigaction = std::mem::zeroed();
                let _ = libc::sigaction(libc::SIGINT, &action, &mut old_int);
                let _ = libc::sigaction(libc::SIGTERM, &action, &mut old_term);
                (old_int, old_term)
            };
            *guard = Some(InstallState {
                active: 1,
                prev_int,
                prev_term,
            });
            Self
        }
    }

    impl Drop for SignalGuard {
        fn drop(&mut self) {
            let mut guard = INSTALL_LOCK
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            if let Some(state) = guard.as_mut() {
                state.active = state.active.saturating_sub(1);
                if state.active == 0 {
                    let prev_int = state.prev_int;
                    let prev_term = state.prev_term;
                    // SAFETY: Restores the previous signal dispositions saved
                    // when the guard was installed.
                    unsafe {
                        let _ = libc::sigaction(libc::SIGINT, &prev_int, std::ptr::null_mut());
                        let _ = libc::sigaction(libc::SIGTERM, &prev_term, std::ptr::null_mut());
                    }
                    *guard = None;
                }
            }
        }
    }

    #[inline]
    pub(super) fn check_cancelled() -> Result<()> {
        let sig = INTERRUPTED_SIGNAL.load(Ordering::Relaxed);
        if sig != 0 {
            bail!("kcount interrupted by signal {sig}; temporary spill files cleaned up");
        }
        Ok(())
    }

    pub(super) fn is_pid_dead(pid: u32) -> bool {
        let Ok(pid_t) = libc::pid_t::try_from(pid) else {
            return false;
        };
        if pid_t <= 0 {
            return false;
        }
        // SAFETY: `kill(pid, 0)` performs only an existence/permission check
        // and delivers no signal to the target process.
        let rc = unsafe { libc::kill(pid_t, 0) };
        if rc == 0 {
            return false;
        }
        std::io::Error::last_os_error().raw_os_error() == Some(libc::ESRCH)
    }
}

#[cfg(not(unix))]
mod signal_guard {
    use super::*;

    pub(super) struct SignalGuard;

    impl SignalGuard {
        pub(super) fn install() -> Self {
            Self
        }
    }

    #[inline]
    pub(super) fn check_cancelled() -> Result<()> {
        Ok(())
    }

    pub(super) fn is_pid_dead(_pid: u32) -> bool {
        false
    }
}

use signal_guard::{SignalGuard, check_cancelled, is_pid_dead};

/// Parses a `grm-kcount-<pid>-<nonce>` directory name, returning `Some(pid)`
/// only when both `<pid>` and `<nonce>` are strictly decimal digits and
/// `pid > 0`.
fn parse_spill_dir_pid(name: &str) -> Option<u32> {
    let rest = name.strip_prefix("grm-kcount-")?;
    let (pid_str, nonce_str) = rest.split_once('-')?;
    if pid_str.is_empty()
        || nonce_str.is_empty()
        || !pid_str.bytes().all(|b| b.is_ascii_digit())
        || !nonce_str.bytes().all(|b| b.is_ascii_digit())
    {
        return None;
    }
    let pid: u32 = pid_str.parse().ok()?;
    (pid > 0).then_some(pid)
}

/// Removes stale `grm-kcount-<old_pid>-<nonce>` spill folders in `parent` whose
/// owning `<old_pid>` is no longer alive, leaving any live process's spill
/// folders untouched.
pub fn cleanup_stale_spill_dirs(parent: &Path) -> usize {
    let Ok(entries) = fs::read_dir(parent) else {
        return 0;
    };
    let my_pid = std::process::id();
    let mut removed = 0usize;
    for entry in entries.flatten() {
        let Ok(name) = entry.file_name().into_string() else {
            continue;
        };
        let Some(pid) = parse_spill_dir_pid(&name) else {
            continue;
        };
        if pid == my_pid || !is_pid_dead(pid) {
            continue;
        }
        let path = entry.path();
        let Ok(meta) = fs::symlink_metadata(&path) else {
            continue;
        };
        if meta.file_type().is_dir()
            && !meta.file_type().is_symlink()
            && fs::remove_dir_all(&path).is_ok()
        {
            removed += 1;
        }
    }
    removed
}

/// RAII temporary spill directory named `grm-kcount-<pid>-<nonce>`, removed
/// cleanly on drop (including early error return and `SIGINT`/`SIGTERM`).
struct SpillDir {
    path: PathBuf,
}

impl SpillDir {
    fn create(base: Option<&Path>) -> Result<Self> {
        let parent = base.map(PathBuf::from).unwrap_or_else(std::env::temp_dir);
        fs::create_dir_all(&parent)
            .with_context(|| format!("create spill parent {}", parent.display()))?;
        let pid = std::process::id();
        for _ in 0..4096 {
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

    fn shard_path(&self, sb: usize) -> PathBuf {
        self.path.join(format!("s{sb:03}.shd"))
    }

    fn subrun_path(&self, p: usize, run: usize) -> PathBuf {
        self.path.join(format!("p{p:04}_r{run:04}.run"))
    }

    fn counts_tmp_path(&self) -> PathBuf {
        self.path.join("ktab_counts.tmp")
    }
}

impl Drop for SpillDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.path);
    }
}

// ---------------------------------------------------------------------------
// Exact bounded-memory length histogram for median, quantiles, and N50
// ---------------------------------------------------------------------------

const INLINE_LEN_BUCKETS: usize = 4096;

/// Exact read/contig length frequency histogram using O(distinct lengths) RAM,
/// avoiding storing a `u64` per read on multi-million-read FASTQ inputs.
#[derive(Clone, Debug)]
struct LenHistogram {
    small: Vec<u64>,
    large: BTreeMap<u64, u64>,
    total_count: u64,
}

impl LenHistogram {
    fn new() -> Self {
        Self {
            small: vec![0u64; INLINE_LEN_BUCKETS],
            large: BTreeMap::new(),
            total_count: 0,
        }
    }

    #[inline]
    fn add(&mut self, len: u64) {
        self.total_count += 1;
        if (len as usize) < INLINE_LEN_BUCKETS {
            self.small[len as usize] += 1;
        } else {
            *self.large.entry(len).or_insert(0) += 1;
        }
    }

    fn merge_from(&mut self, other: &Self) {
        self.total_count += other.total_count;
        for (i, &c) in other.small.iter().enumerate() {
            if c > 0 {
                self.small[i] += c;
            }
        }
        for (&len, &c) in &other.large {
            *self.large.entry(len).or_insert(0) += c;
        }
    }

    /// Returns the exact 0-based sorted order statistic at `target_idx`.
    fn nth(&self, target_idx: u64) -> u64 {
        debug_assert!(target_idx < self.total_count);
        let mut seen = 0u64;
        for (len, &c) in self.small.iter().enumerate() {
            if c == 0 {
                continue;
            }
            seen += c;
            if seen > target_idx {
                return len as u64;
            }
        }
        for (&len, &c) in &self.large {
            seen += c;
            if seen > target_idx {
                return len;
            }
        }
        0
    }

    /// Exact quantile matching sorting a `Vec<u64>` of lengths.
    fn quantile(&self, q: f64) -> f64 {
        if self.total_count == 0 {
            return 0.0;
        }
        let n = self.total_count;
        if (q - 0.5).abs() < 1e-12 {
            if n % 2 == 1 {
                self.nth(n / 2) as f64
            } else {
                (self.nth(n / 2 - 1) as f64 + self.nth(n / 2) as f64) / 2.0
            }
        } else {
            let idx = (((n - 1) as f64) * q).round() as u64;
            self.nth(idx.min(n - 1)) as f64
        }
    }

    /// Exact N50 matching descending-sorted sequence lengths.
    fn n50(&self, total_bases: u64) -> u64 {
        if self.total_count == 0 || total_bases == 0 {
            return 0;
        }
        let target = total_bases as u128;
        let mut running = 0u128;
        for (&len, &c) in self.large.iter().rev() {
            running += (len as u128) * (c as u128);
            if running * 2 >= target {
                return len;
            }
        }
        for (len, &c) in self.small.iter().enumerate().rev() {
            if c == 0 {
                continue;
            }
            running += (len as u128) * (c as u128);
            if running * 2 >= target {
                return len as u64;
            }
        }
        0
    }
}

// ---------------------------------------------------------------------------
// Record-aligned chunked FASTA/FASTQ reader with amortized buffer reuse
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug)]
struct ChunkRecord {
    seq_start: usize,
    seq_len: usize,
    qual_start: usize,
    has_qual: bool,
    /// Number of leading bases in `seq` that overlap the previous slice of a
    /// long FASTA contig (`0` for normal records; `k - 1` for continuation slices).
    overlap: usize,
    /// `Some(contig_id)` only when a single FASTA contig exceeded the chunk
    /// target size and was split across multiple slices.
    split_contig_id: Option<u64>,
    /// True when this slice terminates a split contig.
    split_final: bool,
}

#[derive(Default)]
struct IngestChunk {
    seq_data: Vec<u8>,
    qual_data: Vec<u8>,
    records: Vec<ChunkRecord>,
}

impl IngestChunk {
    fn clear(&mut self) {
        self.seq_data.clear();
        self.qual_data.clear();
        self.records.clear();
    }

    fn approx_bytes(&self) -> usize {
        self.seq_data.len() + self.qual_data.len()
    }
}

struct ChunkedSeqReader {
    reader: Box<dyn BufRead + Send>,
    path_display: String,
    k: usize,
    next_split_id: u64,
    line_buf: Vec<u8>,
    plus_buf: Vec<u8>,
    in_fasta: bool,
    fasta_seq: Vec<u8>,
    fasta_overlap: usize,
    fasta_split_id: Option<u64>,
    saw_fasta: bool,
    saw_fastq: bool,
    eof: bool,
}

impl ChunkedSeqReader {
    fn open(path: &Path, k: usize, next_split_id: u64) -> Result<Self> {
        let reader = fasta::open(path)?;
        Ok(Self {
            reader,
            path_display: path.display().to_string(),
            k,
            next_split_id,
            line_buf: Vec::with_capacity(512),
            plus_buf: Vec::with_capacity(64),
            in_fasta: false,
            fasta_seq: Vec::with_capacity(1 << 16),
            fasta_overlap: 0,
            fasta_split_id: None,
            saw_fasta: false,
            saw_fastq: false,
            eof: false,
        })
    }

    fn fill_chunk(&mut self, chunk: &mut IngestChunk, target_bytes: usize) -> Result<bool> {
        chunk.clear();
        let split_threshold = target_bytes.max(self.k.saturating_mul(4));

        while !self.eof {
            if chunk.approx_bytes() >= target_bytes && !self.in_fasta {
                break;
            }
            self.line_buf.clear();
            let n = read_line_trimmed(&mut self.reader, &mut self.line_buf)?;
            if n == 0 {
                self.eof = true;
                break;
            }
            match self.line_buf.first() {
                Some(b'>') => {
                    if self.in_fasta {
                        self.flush_fasta_to_chunk(chunk, true);
                        if chunk.approx_bytes() >= target_bytes {
                            self.saw_fasta = true;
                            self.in_fasta = true;
                            self.fasta_seq.clear();
                            self.fasta_overlap = 0;
                            self.fasta_split_id = None;
                            return Ok(true);
                        }
                    }
                    self.saw_fasta = true;
                    self.in_fasta = true;
                    self.fasta_seq.clear();
                    self.fasta_overlap = 0;
                    self.fasta_split_id = None;
                }
                Some(b'@') if !self.in_fasta || self.fasta_seq.is_empty() => {
                    if self.in_fasta {
                        self.flush_fasta_to_chunk(chunk, true);
                        self.in_fasta = false;
                    }
                    self.saw_fastq = true;

                    let seq_start = chunk.seq_data.len();
                    let n_seq = read_line_trimmed(&mut self.reader, &mut chunk.seq_data)?;
                    ensure!(
                        n_seq > 0,
                        "{}: truncated FASTQ sequence line",
                        self.path_display
                    );
                    let seq_len = chunk.seq_data.len() - seq_start;

                    self.plus_buf.clear();
                    let n_plus = read_line_trimmed(&mut self.reader, &mut self.plus_buf)?;
                    ensure!(
                        n_plus > 0 && self.plus_buf.first() == Some(&b'+'),
                        "{}: malformed FASTQ separator line",
                        self.path_display
                    );

                    let qual_start = chunk.qual_data.len();
                    let n_qual = read_line_trimmed(&mut self.reader, &mut chunk.qual_data)?;
                    ensure!(
                        n_qual > 0,
                        "{}: truncated FASTQ quality line",
                        self.path_display
                    );
                    let qual_len = chunk.qual_data.len() - qual_start;
                    ensure!(
                        qual_len == seq_len,
                        "{}: FASTQ quality length {qual_len} differs from sequence length {seq_len}",
                        self.path_display
                    );

                    chunk.records.push(ChunkRecord {
                        seq_start,
                        seq_len,
                        qual_start,
                        has_qual: true,
                        overlap: 0,
                        split_contig_id: None,
                        split_final: true,
                    });
                }
                Some(_) if self.in_fasta => {
                    self.fasta_seq.extend_from_slice(&self.line_buf);
                    if self.fasta_seq.len() >= split_threshold {
                        if self.fasta_split_id.is_none() {
                            self.fasta_split_id = Some(self.next_split_id);
                            self.next_split_id += 1;
                        }
                        self.flush_fasta_to_chunk(chunk, false);
                        if chunk.approx_bytes() >= target_bytes {
                            return Ok(true);
                        }
                    }
                }
                Some(_) => bail!("{}: sequence data before any header", self.path_display),
                None => {}
            }
        }

        if self.eof && self.in_fasta {
            self.flush_fasta_to_chunk(chunk, true);
            self.in_fasta = false;
        }

        Ok(!chunk.records.is_empty())
    }

    fn flush_fasta_to_chunk(&mut self, chunk: &mut IngestChunk, is_final: bool) {
        let seq_start = chunk.seq_data.len();
        chunk.seq_data.extend_from_slice(&self.fasta_seq);
        let seq_len = self.fasta_seq.len();
        chunk.records.push(ChunkRecord {
            seq_start,
            seq_len,
            qual_start: 0,
            has_qual: false,
            overlap: self.fasta_overlap,
            split_contig_id: self.fasta_split_id,
            split_final: is_final,
        });
        if is_final {
            self.fasta_seq.clear();
            self.fasta_overlap = 0;
            self.fasta_split_id = None;
        } else {
            let keep = (self.k - 1).min(self.fasta_seq.len());
            let tail_start = self.fasta_seq.len() - keep;
            self.fasta_seq.copy_within(tail_start.., 0);
            self.fasta_seq.truncate(keep);
            self.fasta_overlap = keep;
        }
    }
}

/// Appends one line (up to and including `\n`, stripping trailing `\r`/`\n`) to
/// `dst`. Returns raw bytes consumed from `reader` (`0` at EOF).
#[inline]
fn read_line_trimmed(reader: &mut dyn BufRead, dst: &mut Vec<u8>) -> Result<usize> {
    let start = dst.len();
    let n = reader.read_until(b'\n', dst)?;
    if n == 0 {
        return Ok(0);
    }
    while dst.len() > start && matches!(dst.last(), Some(b'\n' | b'\r')) {
        dst.pop();
    }
    Ok(n)
}

/// Streaming reader used by `--profile-out` to preserve exact record order and
/// optional original header IDs with amortized buffer reuse.
struct ProfileChunkRecord {
    global_idx: u64,
    name_start: usize,
    name_len: usize,
    seq_start: usize,
    seq_len: usize,
    qual_start: usize,
    has_qual: bool,
}

#[derive(Default)]
struct ProfileChunk {
    name_data: Vec<u8>,
    seq_data: Vec<u8>,
    qual_data: Vec<u8>,
    records: Vec<ProfileChunkRecord>,
}

impl ProfileChunk {
    fn clear(&mut self) {
        self.name_data.clear();
        self.seq_data.clear();
        self.qual_data.clear();
        self.records.clear();
    }

    fn approx_bytes(&self) -> usize {
        self.name_data.len() + self.seq_data.len() + self.qual_data.len()
    }
}

struct ProfileSeqReader {
    reader: Box<dyn BufRead + Send>,
    path_display: String,
    keep_record_ids: bool,
    next_global_idx: u64,
    line_buf: Vec<u8>,
    plus_buf: Vec<u8>,
    in_fasta: bool,
    fasta_name: Vec<u8>,
    fasta_seq: Vec<u8>,
    eof: bool,
}

impl ProfileSeqReader {
    fn open(path: &Path, keep_record_ids: bool, start_idx: u64) -> Result<Self> {
        let reader = fasta::open(path)?;
        Ok(Self {
            reader,
            path_display: path.display().to_string(),
            keep_record_ids,
            next_global_idx: start_idx,
            line_buf: Vec::with_capacity(512),
            plus_buf: Vec::with_capacity(64),
            in_fasta: false,
            fasta_name: Vec::with_capacity(256),
            fasta_seq: Vec::with_capacity(1 << 16),
            eof: false,
        })
    }

    fn fill_chunk(&mut self, chunk: &mut ProfileChunk, target_bytes: usize) -> Result<bool> {
        chunk.clear();
        while !self.eof {
            if chunk.approx_bytes() >= target_bytes && !self.in_fasta {
                break;
            }
            self.line_buf.clear();
            let n = read_line_trimmed(&mut self.reader, &mut self.line_buf)?;
            if n == 0 {
                self.eof = true;
                break;
            }
            match self.line_buf.first() {
                Some(b'>') => {
                    if self.in_fasta {
                        self.push_fasta(chunk);
                        if chunk.approx_bytes() >= target_bytes {
                            self.start_fasta_header();
                            return Ok(true);
                        }
                    }
                    self.start_fasta_header();
                }
                Some(b'@') if !self.in_fasta || self.fasta_seq.is_empty() => {
                    if self.in_fasta {
                        self.push_fasta(chunk);
                        self.in_fasta = false;
                    }
                    let name_start = chunk.name_data.len();
                    if self.keep_record_ids {
                        append_header_id(&self.line_buf[1..], &mut chunk.name_data);
                    }
                    let name_len = chunk.name_data.len() - name_start;

                    let seq_start = chunk.seq_data.len();
                    let n_seq = read_line_trimmed(&mut self.reader, &mut chunk.seq_data)?;
                    ensure!(
                        n_seq > 0,
                        "{}: truncated FASTQ sequence line",
                        self.path_display
                    );
                    let seq_len = chunk.seq_data.len() - seq_start;

                    self.plus_buf.clear();
                    let n_plus = read_line_trimmed(&mut self.reader, &mut self.plus_buf)?;
                    ensure!(
                        n_plus > 0 && self.plus_buf.first() == Some(&b'+'),
                        "{}: malformed FASTQ separator line",
                        self.path_display
                    );

                    let qual_start = chunk.qual_data.len();
                    let n_qual = read_line_trimmed(&mut self.reader, &mut chunk.qual_data)?;
                    ensure!(
                        n_qual > 0,
                        "{}: truncated FASTQ quality line",
                        self.path_display
                    );
                    let qual_len = chunk.qual_data.len() - qual_start;
                    ensure!(
                        qual_len == seq_len,
                        "{}: FASTQ quality length {qual_len} differs from sequence length {seq_len}",
                        self.path_display
                    );

                    let global_idx = self.next_global_idx;
                    self.next_global_idx += 1;
                    chunk.records.push(ProfileChunkRecord {
                        global_idx,
                        name_start,
                        name_len,
                        seq_start,
                        seq_len,
                        qual_start,
                        has_qual: true,
                    });
                }
                Some(_) if self.in_fasta => {
                    self.fasta_seq.extend_from_slice(&self.line_buf);
                }
                Some(_) => bail!("{}: sequence data before any header", self.path_display),
                None => {}
            }
        }
        if self.eof && self.in_fasta {
            self.push_fasta(chunk);
            self.in_fasta = false;
        }
        Ok(!chunk.records.is_empty())
    }

    fn start_fasta_header(&mut self) {
        self.fasta_name.clear();
        if self.keep_record_ids {
            append_header_id(&self.line_buf[1..], &mut self.fasta_name);
        }
        self.fasta_seq.clear();
        self.in_fasta = true;
    }

    fn push_fasta(&mut self, chunk: &mut ProfileChunk) {
        let name_start = chunk.name_data.len();
        chunk.name_data.extend_from_slice(&self.fasta_name);
        let name_len = self.fasta_name.len();
        let seq_start = chunk.seq_data.len();
        chunk.seq_data.extend_from_slice(&self.fasta_seq);
        let seq_len = self.fasta_seq.len();
        let global_idx = self.next_global_idx;
        self.next_global_idx += 1;
        chunk.records.push(ProfileChunkRecord {
            global_idx,
            name_start,
            name_len,
            seq_start,
            seq_len,
            qual_start: 0,
            has_qual: false,
        });
        self.fasta_name.clear();
        self.fasta_seq.clear();
    }
}

#[inline]
fn append_header_id(raw: &[u8], out: &mut Vec<u8>) {
    let end = raw
        .iter()
        .position(|b| b.is_ascii_whitespace())
        .unwrap_or(raw.len());
    out.extend_from_slice(&raw[..end]);
}

// ---------------------------------------------------------------------------
// Minimizer super-k-mer extraction & partition store
// ---------------------------------------------------------------------------

/// Deterministic SplitMix64 mixer (Steele, Lea, and Flood, OOPSLA 2014, public
/// domain) for strand-invariant minimizer ordering.
#[inline]
fn mix_minimizer(mut x: u64) -> u64 {
    x ^= 0x5a5a_5a5a_5a5a_5a5a;
    x ^= x >> 30;
    x = x.wrapping_mul(SPLITMIX64_MUL_1);
    x ^= x >> 27;
    x = x.wrapping_mul(SPLITMIX64_MUL_2);
    x ^ (x >> 31)
}

/// Default minimizer length `m` for a given `k` (`1 <= m <= min(k, 31)`).
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

/// Per-worker amortized reusable scratch buffers for parallel Stage 1 super-k-mer
/// extraction.
struct WorkerExtractBuf {
    k: usize,
    m: usize,
    m_mask: u64,
    m_shift: u32,
    p_mask: usize,
    run_codes: Vec<u8>,
    min_hashes: Vec<u64>,
    deque: VecDeque<(usize, u64)>,
    part_bufs: Vec<Vec<u8>>,
    total_reads: u64,
    total_bases: u64,
    valid_bases: u64,
    gc_bases: u64,
    masked_bases: u64,
    raw_lens: LenHistogram,
    usable_lens: LenHistogram,
    split_updates: Vec<(u64, u64, u64, bool)>,
}

impl WorkerExtractBuf {
    fn new(k: usize, m: usize, partitions: usize) -> Self {
        let m_bits = 2 * m as u32;
        let m_mask = (1u64 << m_bits) - 1;
        Self {
            k,
            m,
            m_mask,
            m_shift: 2 * (m as u32 - 1),
            p_mask: partitions - 1,
            run_codes: Vec::with_capacity(4096),
            min_hashes: Vec::with_capacity(4096),
            deque: VecDeque::with_capacity(64),
            part_bufs: vec![Vec::new(); partitions],
            total_reads: 0,
            total_bases: 0,
            valid_bases: 0,
            gc_bases: 0,
            masked_bases: 0,
            raw_lens: LenHistogram::new(),
            usable_lens: LenHistogram::new(),
            split_updates: Vec::new(),
        }
    }

    fn reset_batch(&mut self) {
        for buf in &mut self.part_bufs {
            buf.clear();
            if buf.capacity() > 64 * 1024 {
                *buf = Vec::new();
            }
        }
        self.split_updates.clear();
    }

    fn process_chunk(&mut self, chunk: &IngestChunk, min_qual: u8) {
        let qual_floor = 33u8.saturating_add(min_qual);
        for rec in &chunk.records {
            let seq = &chunk.seq_data[rec.seq_start..rec.seq_start + rec.seq_len];
            let qual = if rec.has_qual {
                Some(&chunk.qual_data[rec.qual_start..rec.qual_start + rec.seq_len])
            } else {
                None
            };
            let use_qual = min_qual > 0 && qual.is_some();
            let qslice = qual.unwrap_or(&[]);

            let mut slice_valid = 0u64;
            let mut slice_gc = 0u64;
            let mut slice_masked = 0u64;

            self.run_codes.clear();
            for (idx, &b) in seq.iter().enumerate() {
                let c = CODE[b as usize];
                let low_q = use_qual && qslice[idx] < qual_floor;
                let count_qc = idx >= rec.overlap;
                if c > 3 || low_q {
                    if count_qc {
                        slice_masked += 1;
                    }
                    if self.run_codes.len() >= self.k {
                        self.emit_run();
                    }
                    self.run_codes.clear();
                } else {
                    if count_qc {
                        slice_valid += 1;
                        if c == 1 || c == 2 {
                            slice_gc += 1;
                        }
                    }
                    self.run_codes.push(c);
                }
            }
            if self.run_codes.len() >= self.k {
                self.emit_run();
            }
            self.run_codes.clear();

            let non_overlap_len = (rec.seq_len - rec.overlap) as u64;
            self.total_bases += non_overlap_len;
            self.valid_bases += slice_valid;
            self.gc_bases += slice_gc;
            self.masked_bases += slice_masked;

            if let Some(cid) = rec.split_contig_id {
                self.split_updates
                    .push((cid, non_overlap_len, slice_valid, rec.split_final));
            } else {
                self.total_reads += 1;
                self.raw_lens.add(non_overlap_len);
                self.usable_lens.add(slice_valid);
            }
        }
    }

    fn emit_run(&mut self) {
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
                self.push_superkmer(cur_min_hash, sk_start, end);
                sk_start = pos;
                cur_min_hash = win_min;
            }
        }
        if have_sk {
            self.push_superkmer(cur_min_hash, sk_start, len);
        }
    }

    fn push_superkmer(&mut self, min_hash: u64, start: usize, end: usize) {
        let sk_len = end - start;
        debug_assert!((self.k..=MAX_SUPERKMER_LEN).contains(&sk_len));
        let part = (min_hash as usize) & self.p_mask;
        let packed_bytes = sk_len.div_ceil(4);
        let buf = &mut self.part_bufs[part];
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
    }
}

/// Bounded-memory partition store accumulating packed super-k-mers across Stage 1
/// and spilling to disk whenever buffered bytes reach `budget_bytes`.
struct PartitionStore {
    buffers: Vec<Vec<u8>>,
    part_total_bytes: Vec<u64>,
    buffered_bytes: usize,
    total_bytes: u64,
    budget_bytes: usize,
    spill_dir: Option<SpillDir>,
    spilled_flags: Vec<bool>,
    tmp_base: Option<PathBuf>,
}

impl PartitionStore {
    fn new(partitions: usize, budget_bytes: usize, tmp_base: Option<PathBuf>) -> Self {
        Self {
            buffers: vec![Vec::new(); partitions],
            part_total_bytes: vec![0u64; partitions],
            buffered_bytes: 0,
            total_bytes: 0,
            budget_bytes: budget_bytes.max(64 * 1024),
            spill_dir: None,
            spilled_flags: vec![false; partitions],
            tmp_base,
        }
    }

    fn absorb_worker(&mut self, worker_bufs: &mut [Vec<u8>]) -> Result<()> {
        for (p, src) in worker_bufs.iter_mut().enumerate() {
            if src.is_empty() {
                continue;
            }
            let n = src.len();
            self.buffers[p].extend_from_slice(src);
            self.part_total_bytes[p] += n as u64;
            self.buffered_bytes += n;
            self.total_bytes += n as u64;
            src.clear();
        }
        if self.buffered_bytes >= self.budget_bytes {
            self.spill_all()?;
        }
        Ok(())
    }

    fn spill_all(&mut self) -> Result<()> {
        check_cancelled()?;
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
            if buf.capacity() > 32 * 1024 {
                *buf = Vec::new();
            }
            self.spilled_flags[p] = true;
        }
        self.buffered_bytes = 0;
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Count saturation & checked u64 run-merging helpers (Finding 5)
// ---------------------------------------------------------------------------

/// Converts an exact `u64` occurrence count into a `u32` table count capped at
/// `2^32 - 1` (`MAX_COUNT`), returning `(count, saturated)`. Refuses saturation
/// with a clear error unless `allow_saturation` is true.
pub fn finalize_u32_count(raw: u64, allow_saturation: bool) -> Result<(u32, bool)> {
    if raw <= MAX_COUNT as u64 {
        Ok((raw as u32, false))
    } else if allow_saturation {
        Ok((MAX_COUNT, true))
    } else {
        bail!(
            "k-mer occurrence count ({raw}) exceeds u32::MAX ({MAX_COUNT}); re-run with --allow-saturation to clamp counts at 2^32 - 1"
        );
    }
}

/// Merges multiple strictly sorted `(key, u64_count)` runs using checked `u64`
/// addition and enforces [`finalize_u32_count`], returning `(key, u32_count)`
/// entries and the number of saturated keys.
pub fn merge_counted_runs<K: KmerKey>(
    runs: &[Vec<(K, u64)>],
    allow_saturation: bool,
) -> Result<(Vec<(K, u32)>, u64)> {
    let mut heap: BinaryHeap<Reverse<(K, usize, usize)>> = BinaryHeap::new();
    for (r_idx, run) in runs.iter().enumerate() {
        if let Some(&(k0, _)) = run.first() {
            heap.push(Reverse((k0, r_idx, 0)));
        }
    }
    let mut out = Vec::new();
    let mut saturated = 0u64;
    while let Some(Reverse((key, r_idx, pos))) = heap.pop() {
        let mut total = runs[r_idx][pos].1;
        if pos + 1 < runs[r_idx].len() {
            heap.push(Reverse((runs[r_idx][pos + 1].0, r_idx, pos + 1)));
        }
        while let Some(&Reverse((next_key, nr, np))) = heap.peek() {
            if next_key != key {
                break;
            }
            heap.pop();
            total = total
                .checked_add(runs[nr][np].1)
                .context("k-mer occurrence count overflowed u64")?;
            if np + 1 < runs[nr].len() {
                heap.push(Reverse((runs[nr][np + 1].0, nr, np + 1)));
            }
        }
        let (c32, is_sat) = finalize_u32_count(total, allow_saturation)?;
        if is_sat {
            saturated += 1;
        }
        out.push((key, c32));
    }
    Ok((out, saturated))
}

// ---------------------------------------------------------------------------
// Stage 2: Bounded-memory partition counting & 256-super-bucket shard pool
// ---------------------------------------------------------------------------

struct PartitionStats {
    hist_counts: FxHashMap<u32, u64>,
    distinct_total: u64,
    distinct_retained: u64,
    total_occurrences: u64,
    max_obs: u32,
    max_retained: u32,
    saturated_total: u64,
    saturated_retained: u64,
    all_retained_one: bool,
}

impl PartitionStats {
    fn new() -> Self {
        Self {
            hist_counts: FxHashMap::default(),
            distinct_total: 0,
            distinct_retained: 0,
            total_occurrences: 0,
            max_obs: 0,
            max_retained: 0,
            saturated_total: 0,
            saturated_retained: 0,
            all_retained_one: true,
        }
    }

    fn merge_from(&mut self, other: Self) {
        self.distinct_total += other.distinct_total;
        self.distinct_retained += other.distinct_retained;
        self.total_occurrences += other.total_occurrences;
        if other.max_obs > self.max_obs {
            self.max_obs = other.max_obs;
        }
        if other.max_retained > self.max_retained {
            self.max_retained = other.max_retained;
        }
        self.saturated_total += other.saturated_total;
        self.saturated_retained += other.saturated_retained;
        if !other.all_retained_one {
            self.all_retained_one = false;
        }
        for (c, n) in other.hist_counts {
            *self.hist_counts.entry(c).or_insert(0) += n;
        }
    }
}

/// Sink for Stage 2 retained `(K, u32)` pairs: either 256 super-bucket disk
/// shard files (bounded streaming path) or 1,024 in-memory bucket vectors (fast
/// path when the entire dataset fits comfortably below `--max-memory-mb`).
enum Stage2Sink<'a, K: KmerKey> {
    Sharded(Box<Mutex<ShardWriterPool<'a>>>),
    InMemory(Mutex<Vec<Vec<(K, u32)>>>),
}

struct ShardWriterPool<'a> {
    spill_dir: &'a SpillDir,
    buffers: Vec<Vec<u8>>,
    buffered_bytes: usize,
    budget_bytes: usize,
    written_flags: [bool; SUPER_BUCKETS],
}

impl<'a> ShardWriterPool<'a> {
    fn new(spill_dir: &'a SpillDir, budget_bytes: usize) -> Self {
        Self {
            spill_dir,
            buffers: vec![Vec::new(); SUPER_BUCKETS],
            buffered_bytes: 0,
            budget_bytes: budget_bytes.max(256 * 1024),
            written_flags: [false; SUPER_BUCKETS],
        }
    }

    fn push_entries<K: KmerKey>(&mut self, k: usize, entries: &[(K, u32)]) -> Result<()> {
        let entry_bytes = K::WIDTH_BYTES + 4;
        for &(key, count) in entries {
            let sb = key.bucket_of(k, BUCKET_BITS) / BUCKETS_PER_SUPER;
            let buf = &mut self.buffers[sb];
            key.with_le_bytes(|kb| buf.extend_from_slice(kb));
            buf.extend_from_slice(&count.to_le_bytes());
            self.buffered_bytes += entry_bytes;
        }
        if self.buffered_bytes >= self.budget_bytes {
            self.flush_all()?;
        }
        Ok(())
    }

    fn flush_all(&mut self) -> Result<()> {
        check_cancelled()?;
        for (sb, buf) in self.buffers.iter_mut().enumerate() {
            if buf.is_empty() {
                continue;
            }
            let path = self.spill_dir.shard_path(sb);
            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(&path)
                .with_context(|| format!("open output shard {}", path.display()))?;
            file.write_all(buf)?;
            buf.clear();
            if buf.capacity() > 32 * 1024 {
                *buf = Vec::new();
            }
            self.written_flags[sb] = true;
        }
        self.buffered_bytes = 0;
        Ok(())
    }
}

fn push_batch_to_sink<K: KmerKey>(
    sink: &Stage2Sink<'_, K>,
    k: usize,
    batch: &mut Vec<(K, u32)>,
) -> Result<()> {
    if batch.is_empty() {
        return Ok(());
    }
    match sink {
        Stage2Sink::Sharded(pool_mtx) => {
            let mut pool = pool_mtx
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            pool.push_entries(k, batch)?;
        }
        Stage2Sink::InMemory(buckets_mtx) => {
            let mut buckets = buckets_mtx
                .lock()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            for &(key, count) in batch.iter() {
                let b = key.bucket_of(k, BUCKET_BITS);
                buckets[b].push((key, count));
            }
        }
    }
    batch.clear();
    Ok(())
}

/// Reusable per-worker scratch buffers for Stage 2 partition counting.
struct Stage2WorkerScratch<K: KmerKey> {
    keys: Vec<K>,
    radix_scratch: Vec<K>,
    emit_batch: Vec<(K, u32)>,
}

impl<K: KmerKey> Stage2WorkerScratch<K> {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            radix_scratch: Vec::new(),
            emit_batch: Vec::with_capacity(16_384),
        }
    }
}

/// Reads packed super-k-mers sequentially from an optional spill file followed
/// by an optional in-memory buffer without loading an entire spill file at once.
struct PartitionStream {
    file: Option<BufReader<File>>,
    file_path: Option<PathBuf>,
    mem: Vec<u8>,
    mem_pos: usize,
}

impl PartitionStream {
    fn new(spill_path: Option<PathBuf>, mem: Vec<u8>) -> Result<Self> {
        let file = match &spill_path {
            Some(path) => Some(BufReader::with_capacity(
                64 * 1024,
                File::open(path).with_context(|| format!("read spill {}", path.display()))?,
            )),
            None => None,
        };
        Ok(Self {
            file,
            file_path: spill_path,
            mem,
            mem_pos: 0,
        })
    }

    /// Reads the next packed super-k-mer into `packed_buf`, returning `Ok(Some(sk_len))`
    /// or `Ok(None)` at EOF.
    fn next_superkmer(&mut self, k: usize, packed_buf: &mut [u8; 64]) -> Result<Option<usize>> {
        if let Some(ref mut f) = self.file {
            let mut len_byte = [0u8; 1];
            match f.read(&mut len_byte)? {
                0 => {
                    self.file = None;
                    if let Some(path) = self.file_path.take() {
                        let _ = fs::remove_file(&path);
                    }
                }
                1 => {
                    let sk_len = len_byte[0] as usize;
                    ensure!(sk_len >= k, "corrupt super-k-mer length {sk_len} < k={k}");
                    let nbytes = sk_len.div_ceil(4);
                    f.read_exact(&mut packed_buf[..nbytes])
                        .context("truncated super-k-mer partition file")?;
                    return Ok(Some(sk_len));
                }
                _ => unreachable!(),
            }
        }
        if self.mem_pos < self.mem.len() {
            let sk_len = self.mem[self.mem_pos] as usize;
            self.mem_pos += 1;
            let nbytes = sk_len.div_ceil(4);
            ensure!(
                sk_len >= k && self.mem_pos + nbytes <= self.mem.len(),
                "corrupt in-memory super-k-mer partition buffer"
            );
            packed_buf[..nbytes].copy_from_slice(&self.mem[self.mem_pos..self.mem_pos + nbytes]);
            self.mem_pos += nbytes;
            return Ok(Some(sk_len));
        }
        Ok(None)
    }
}

impl Drop for PartitionStream {
    fn drop(&mut self) {
        self.file = None;
        if let Some(path) = self.file_path.take() {
            let _ = fs::remove_file(&path);
        }
    }
}

/// Writes a sorted `(K, u64)` sub-run to disk in compact little-endian pairs.
fn write_subrun_file<K: KmerKey>(path: &Path, run: &[(K, u64)]) -> Result<()> {
    let mut w = BufWriter::with_capacity(
        64 * 1024,
        File::create(path).with_context(|| format!("create subrun {}", path.display()))?,
    );
    for &(key, count) in run {
        key.with_le_bytes(|kb| w.write_all(kb))?;
        w.write_all(&count.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

/// Streaming reader over one sorted `(K, u64)` sub-run (either in-memory or on disk).
enum SubRunReader<K: KmerKey> {
    Mem {
        data: Vec<(K, u64)>,
        pos: usize,
    },
    Disk {
        reader: BufReader<File>,
        path: PathBuf,
    },
}

impl<K: KmerKey> SubRunReader<K> {
    fn next_pair(&mut self) -> Result<Option<(K, u64)>> {
        match self {
            Self::Mem { data, pos } => {
                if *pos < data.len() {
                    let item = data[*pos];
                    *pos += 1;
                    Ok(Some(item))
                } else {
                    Ok(None)
                }
            }
            Self::Disk { reader, .. } => {
                let mut buf = [0u8; 24];
                let entry_len = K::WIDTH_BYTES + 8;
                let mut filled = 0usize;
                while filled < entry_len {
                    let n = reader.read(&mut buf[filled..entry_len])?;
                    if n == 0 {
                        ensure!(filled == 0, "truncated sub-run spill file");
                        return Ok(None);
                    }
                    filled += n;
                }
                let key = decode_le_key::<K>(&buf[..K::WIDTH_BYTES]);
                let count = u64::from_le_bytes(buf[K::WIDTH_BYTES..entry_len].try_into().unwrap());
                Ok(Some((key, count)))
            }
        }
    }
}

impl<K: KmerKey> Drop for SubRunReader<K> {
    fn drop(&mut self) {
        if let Self::Disk { path, .. } = self {
            let _ = fs::remove_file(path);
        }
    }
}

#[inline]
fn decode_le_key<K: KmerKey>(bytes: &[u8]) -> K {
    if K::WIDTH_BYTES == 8 {
        let v = u64::from_le_bytes(bytes[..8].try_into().unwrap());
        K::from_u128(v as u128)
    } else {
        let v = u128::from_le_bytes(bytes[..16].try_into().unwrap());
        K::from_u128(v)
    }
}

#[allow(clippy::too_many_arguments)]
fn process_partition_bounded<K: KmerKey>(
    part_idx: usize,
    space: &K::Space,
    k: usize,
    min_count: u32,
    max_count: Option<u32>,
    allow_saturation: bool,
    max_keys_per_chunk: usize,
    subrun_mem_cap_bytes: usize,
    spill_dir: Option<&SpillDir>,
    spill_path: Option<PathBuf>,
    mem_buf: Vec<u8>,
    worker: &mut Stage2WorkerScratch<K>,
    sink: &Stage2Sink<'_, K>,
) -> Result<PartitionStats> {
    check_cancelled()?;
    let mut stream = PartitionStream::new(spill_path, mem_buf)?;
    let mut stats = PartitionStats::new();
    let mut packed_buf = [0u8; 64];
    let mut sub_runs: Vec<SubRunReader<K>> = Vec::new();
    let mut mem_subrun_bytes = 0usize;
    let mut disk_subrun_seq = 0usize;
    let mut reached_eof = false;

    loop {
        check_cancelled()?;
        worker.keys.clear();
        while worker.keys.len() < max_keys_per_chunk {
            match stream.next_superkmer(k, &mut packed_buf)? {
                Some(sk_len) => {
                    let nbytes = sk_len.div_ceil(4);
                    K::for_each_mixed_in_packed_2bit(space, &packed_buf[..nbytes], sk_len, |m| {
                        worker.keys.push(m);
                    });
                }
                None => {
                    reached_eof = true;
                    break;
                }
            }
        }

        if worker.keys.is_empty() {
            break;
        }
        stats.total_occurrences += worker.keys.len() as u64;
        K::radix_sort(&mut worker.keys, &mut worker.radix_scratch);

        if reached_eof && sub_runs.is_empty() {
            let mut i = 0usize;
            while i < worker.keys.len() {
                let key = worker.keys[i];
                let mut j = i + 1;
                while j < worker.keys.len() && worker.keys[j] == key {
                    j += 1;
                }
                let raw_count = (j - i) as u64;
                record_counted_kmer(
                    key,
                    raw_count,
                    min_count,
                    max_count,
                    allow_saturation,
                    &mut stats,
                    &mut worker.emit_batch,
                )?;
                if worker.emit_batch.len() >= 16_384 {
                    push_batch_to_sink(sink, k, &mut worker.emit_batch)?;
                }
                i = j;
            }
            push_batch_to_sink(sink, k, &mut worker.emit_batch)?;
            break;
        }

        let mut collapsed: Vec<(K, u64)> = Vec::new();
        let mut i = 0usize;
        while i < worker.keys.len() {
            let key = worker.keys[i];
            let mut j = i + 1;
            while j < worker.keys.len() && worker.keys[j] == key {
                j += 1;
            }
            collapsed.push((key, (j - i) as u64));
            i = j;
        }
        let run_bytes = collapsed.len() * std::mem::size_of::<(K, u64)>();
        if let Some(dir) = spill_dir
            && mem_subrun_bytes + run_bytes > subrun_mem_cap_bytes
        {
            let run_path = dir.subrun_path(part_idx, disk_subrun_seq);
            disk_subrun_seq += 1;
            write_subrun_file(&run_path, &collapsed)?;
            drop(collapsed);
            let reader = BufReader::with_capacity(
                32 * 1024,
                File::open(&run_path)
                    .with_context(|| format!("reopen subrun {}", run_path.display()))?,
            );
            sub_runs.push(SubRunReader::Disk {
                reader,
                path: run_path,
            });
        } else {
            mem_subrun_bytes += run_bytes;
            sub_runs.push(SubRunReader::Mem {
                data: collapsed,
                pos: 0,
            });
        }

        if reached_eof {
            break;
        }
    }

    worker.keys.clear();
    worker.radix_scratch.clear();
    if worker.keys.capacity() > max_keys_per_chunk.saturating_mul(2) {
        worker.keys = Vec::new();
    }
    if worker.radix_scratch.capacity() > max_keys_per_chunk.saturating_mul(2) {
        worker.radix_scratch = Vec::new();
    }

    if !sub_runs.is_empty() {
        let mut heap: BinaryHeap<Reverse<(K, usize, u64)>> = BinaryHeap::new();
        for (r_idx, sr) in sub_runs.iter_mut().enumerate() {
            if let Some((key, count)) = sr.next_pair()? {
                heap.push(Reverse((key, r_idx, count)));
            }
        }
        while let Some(Reverse((key, r_idx, first_count))) = heap.pop() {
            let mut total = first_count;
            if let Some((nk, nc)) = sub_runs[r_idx].next_pair()? {
                heap.push(Reverse((nk, r_idx, nc)));
            }
            while let Some(&Reverse((peek_key, _, _))) = heap.peek() {
                if peek_key != key {
                    break;
                }
                let Reverse((_, nr, nc_val)) = heap.pop().unwrap();
                total = total
                    .checked_add(nc_val)
                    .context("k-mer occurrence count overflowed u64")?;
                if let Some((nk, nc)) = sub_runs[nr].next_pair()? {
                    heap.push(Reverse((nk, nr, nc)));
                }
            }
            record_counted_kmer(
                key,
                total,
                min_count,
                max_count,
                allow_saturation,
                &mut stats,
                &mut worker.emit_batch,
            )?;
            if worker.emit_batch.len() >= 16_384 {
                push_batch_to_sink(sink, k, &mut worker.emit_batch)?;
            }
        }
        push_batch_to_sink(sink, k, &mut worker.emit_batch)?;
    }

    Ok(stats)
}

#[inline]
fn record_counted_kmer<K: KmerKey>(
    key: K,
    raw_count: u64,
    min_count: u32,
    max_count: Option<u32>,
    allow_saturation: bool,
    stats: &mut PartitionStats,
    emit_batch: &mut Vec<(K, u32)>,
) -> Result<()> {
    let (count, saturated) = finalize_u32_count(raw_count, allow_saturation)?;
    stats.distinct_total += 1;
    if count > stats.max_obs {
        stats.max_obs = count;
    }
    if saturated {
        stats.saturated_total += 1;
    }
    *stats.hist_counts.entry(count).or_insert(0) += 1;

    let keep = count >= min_count && max_count.is_none_or(|mc| count <= mc);
    if keep {
        stats.distinct_retained += 1;
        if count > stats.max_retained {
            stats.max_retained = count;
        }
        if count != 1 {
            stats.all_retained_one = false;
        }
        if saturated {
            stats.saturated_retained += 1;
        }
        emit_batch.push((key, count));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Stage 3: Streaming verified .grk and .ktab writers (O(1) bucket RAM)
// ---------------------------------------------------------------------------

fn hash_file_excluding_slot(path: &Path, slot_start: usize, slot_end: usize) -> Result<[u8; 32]> {
    let mut file =
        File::open(path).with_context(|| format!("reopen {} for SHA-256", path.display()))?;
    let file_len = file.metadata()?.len();
    ensure!(
        slot_start <= slot_end && (slot_end as u64) <= file_len,
        "output file truncated before digest slot"
    );
    let mut sha = Sha256Context::new(&SHA256);
    let mut buf = [0u8; 64 * 1024];

    let mut pos = 0u64;
    let prefix_end = slot_start as u64;
    while pos < prefix_end {
        let take = ((prefix_end - pos) as usize).min(buf.len());
        file.read_exact(&mut buf[..take])?;
        sha.update(&buf[..take]);
        pos += take as u64;
    }

    file.seek(SeekFrom::Start(slot_end as u64))?;
    pos = slot_end as u64;
    while pos < file_len {
        let take = ((file_len - pos) as usize).min(buf.len());
        file.read_exact(&mut buf[..take])?;
        sha.update(&buf[..take]);
        pos += take as u64;
    }

    let digest = sha.finish();
    let mut out = [0u8; 32];
    out.copy_from_slice(digest.as_ref());
    Ok(out)
}

struct StreamingGrkWriter {
    target: PathBuf,
    part: PathBuf,
    writer: Option<BufWriter<File>>,
    name: String,
    k: usize,
    modality: InputModality,
    offsets: [u64; BUCKETS + 1],
    bucket_counts: [u32; BUCKETS],
    cur_offset: u64,
    finished: bool,
}

impl StreamingGrkWriter {
    fn create(path: &Path, name: &str, k: usize, modality: InputModality) -> Result<Self> {
        ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
        let _ = u32::try_from(name.len()).context("sample name exceeds u32 length")?;
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }
        confine::regular_or_absent(path)?;
        let part = confine::suffixed(path, ".part");
        let file = confine::create_part(&part)?;
        let mut writer = BufWriter::with_capacity(128 * 1024, file);
        let table_len = (BUCKETS + 1) * 8 + BUCKETS * 4;
        let data_start = GRK_HEADER_LEN + name.len() + table_len;
        writer.seek(SeekFrom::Start(data_start as u64))?;
        Ok(Self {
            target: path.to_path_buf(),
            part,
            writer: Some(writer),
            name: name.to_string(),
            k,
            modality,
            offsets: [0u64; BUCKETS + 1],
            bucket_counts: [0u32; BUCKETS],
            cur_offset: 0,
            finished: false,
        })
    }

    fn write_bucket(&mut self, bucket: usize, count: u32, key_payload: &[u8]) -> Result<()> {
        self.offsets[bucket] = self.cur_offset;
        self.bucket_counts[bucket] = count;
        if !key_payload.is_empty() {
            self.writer.as_mut().unwrap().write_all(key_payload)?;
            self.cur_offset = self
                .cur_offset
                .checked_add(key_payload.len() as u64)
                .context("grk payload offset overflow")?;
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.offsets[BUCKETS] = self.cur_offset;
        let n_kmers = self
            .bucket_counts
            .iter()
            .try_fold(0u64, |acc, &c| acc.checked_add(c as u64))
            .context("grk n_kmers overflow")?;
        let name_len = u32::try_from(self.name.len()).context("sample name exceeds u32")?;
        let key_width: u8 = if self.k <= MAX_K_U64 { 8 } else { 16 };

        let mut prefix = [0u8; GRK_DIGEST_START];
        prefix[0..8].copy_from_slice(MAGIC_GRKSET);
        prefix[8] = self.k as u8;
        prefix[9] = BUCKET_BITS as u8;
        prefix[10] = key_width;
        prefix[11] = self.modality.as_byte();
        prefix[16..24].copy_from_slice(&n_kmers.to_le_bytes());
        prefix[24..28].copy_from_slice(&name_len.to_le_bytes());

        let w = self.writer.as_mut().unwrap();
        w.flush()?;
        w.seek(SeekFrom::Start(0))?;
        w.write_all(&prefix)?;
        w.write_all(&[0u8; 32])?;
        w.write_all(self.name.as_bytes())?;
        for o in self.offsets {
            w.write_all(&o.to_le_bytes())?;
        }
        for c in self.bucket_counts {
            w.write_all(&c.to_le_bytes())?;
        }
        w.flush()?;

        let digest = hash_file_excluding_slot(&self.part, GRK_DIGEST_START, GRK_DIGEST_END)?;
        w.seek(SeekFrom::Start(GRK_DIGEST_START as u64))?;
        w.write_all(&digest)?;
        w.flush()?;
        drop(self.writer.take());

        confine::replace_with(&self.part, &self.target)?;
        self.finished = true;
        Ok(())
    }
}

impl Drop for StreamingGrkWriter {
    fn drop(&mut self) {
        self.writer = None;
        if !self.finished && confine::regular_or_absent(&self.part).unwrap_or(false) {
            let _ = fs::remove_file(&self.part);
        }
    }
}

struct StreamingKtabWriter {
    target: PathBuf,
    part: PathBuf,
    writer: Option<BufWriter<File>>,
    counts_tmp_path: Option<PathBuf>,
    counts_tmp_writer: Option<BufWriter<File>>,
    counts_mem_buf: Vec<u8>,
    name: String,
    k: usize,
    modality: InputModality,
    min_count: u32,
    max_count_observed: u32,
    saturated_counts: u64,
    presence_only: bool,
    key_offsets: [u64; BUCKETS + 1],
    count_offsets: [u64; BUCKETS + 1],
    bucket_counts: [u32; BUCKETS],
    cur_key_offset: u64,
    cur_count_offset: u64,
    finished: bool,
}

impl StreamingKtabWriter {
    #[allow(clippy::too_many_arguments)]
    fn create(
        path: &Path,
        name: &str,
        k: usize,
        modality: InputModality,
        min_count: u32,
        max_count_observed: u32,
        saturated_counts: u64,
        presence_only: bool,
        counts_spill_path: Option<PathBuf>,
    ) -> Result<Self> {
        ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
        let _ = u32::try_from(name.len()).context("sample name exceeds u32 length")?;
        if let Some(parent) = path.parent()
            && !parent.as_os_str().is_empty()
        {
            fs::create_dir_all(parent)?;
        }
        confine::regular_or_absent(path)?;
        let part = confine::suffixed(path, ".part");
        let file = confine::create_part(&part)?;
        let mut writer = BufWriter::with_capacity(128 * 1024, file);
        let table_len = 2 * (BUCKETS + 1) * 8 + BUCKETS * 4;
        let key_data_start = KTAB_HEADER_LEN + name.len() + table_len;
        writer.seek(SeekFrom::Start(key_data_start as u64))?;

        let (counts_tmp_path, counts_tmp_writer) = if !presence_only {
            if let Some(cp) = counts_spill_path {
                let f = File::create(&cp)
                    .with_context(|| format!("create temporary counts file {}", cp.display()))?;
                (Some(cp), Some(BufWriter::with_capacity(64 * 1024, f)))
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };

        Ok(Self {
            target: path.to_path_buf(),
            part,
            writer: Some(writer),
            counts_tmp_path,
            counts_tmp_writer,
            counts_mem_buf: Vec::new(),
            name: name.to_string(),
            k,
            modality,
            min_count,
            max_count_observed,
            saturated_counts,
            presence_only,
            key_offsets: [0u64; BUCKETS + 1],
            count_offsets: [0u64; BUCKETS + 1],
            bucket_counts: [0u32; BUCKETS],
            cur_key_offset: 0,
            cur_count_offset: 0,
            finished: false,
        })
    }

    fn write_bucket(
        &mut self,
        bucket: usize,
        count: u32,
        key_payload: &[u8],
        count_payload: &[u8],
    ) -> Result<()> {
        self.key_offsets[bucket] = self.cur_key_offset;
        self.bucket_counts[bucket] = count;
        if !key_payload.is_empty() {
            self.writer.as_mut().unwrap().write_all(key_payload)?;
            self.cur_key_offset = self
                .cur_key_offset
                .checked_add(key_payload.len() as u64)
                .context("ktab key offset overflow")?;
        }
        if !self.presence_only {
            self.count_offsets[bucket] = self.cur_count_offset;
            if !count_payload.is_empty() {
                if let Some(ref mut cw) = self.counts_tmp_writer {
                    cw.write_all(count_payload)?;
                } else {
                    self.counts_mem_buf.extend_from_slice(count_payload);
                }
                self.cur_count_offset = self
                    .cur_count_offset
                    .checked_add(count_payload.len() as u64)
                    .context("ktab count offset overflow")?;
            }
        }
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        self.key_offsets[BUCKETS] = self.cur_key_offset;
        if !self.presence_only {
            self.count_offsets[BUCKETS] = self.cur_count_offset;
        }
        let n_kmers = self
            .bucket_counts
            .iter()
            .try_fold(0u64, |acc, &c| acc.checked_add(c as u64))
            .context("ktab n_kmers overflow")?;
        ensure!(
            self.saturated_counts <= n_kmers,
            "saturated_counts {} exceeds n_kmers {n_kmers}",
            self.saturated_counts
        );

        let w = self.writer.as_mut().unwrap();
        if !self.presence_only {
            if let Some(mut cw) = self.counts_tmp_writer.take() {
                cw.flush()?;
                drop(cw);
                if let Some(cp) = self.counts_tmp_path.take() {
                    let mut cf = File::open(&cp)
                        .with_context(|| format!("reopen temporary counts {}", cp.display()))?;
                    let mut copy_buf = [0u8; 64 * 1024];
                    loop {
                        let n = cf.read(&mut copy_buf)?;
                        if n == 0 {
                            break;
                        }
                        w.write_all(&copy_buf[..n])?;
                    }
                    drop(cf);
                    let _ = fs::remove_file(&cp);
                }
            } else if !self.counts_mem_buf.is_empty() {
                w.write_all(&self.counts_mem_buf)?;
                self.counts_mem_buf.clear();
                self.counts_mem_buf.shrink_to_fit();
            }
        }

        let name_len = u32::try_from(self.name.len()).context("sample name exceeds u32")?;
        let key_width: u8 = if self.k <= MAX_K_U64 { 8 } else { 16 };
        let mut flags: u8 = 0;
        if self.presence_only {
            flags |= FLAG_PRESENCE_ONLY;
        }
        if self.saturated_counts > 0 {
            flags |= FLAG_SATURATED_COUNTS;
        }

        let mut prefix = [0u8; KTAB_DIGEST_START];
        prefix[0..8].copy_from_slice(MAGIC_KTAB);
        prefix[8] = self.k as u8;
        prefix[9] = BUCKET_BITS as u8;
        prefix[10] = key_width;
        prefix[11] = flags;
        prefix[12] = self.modality.as_byte();
        prefix[16..20].copy_from_slice(&self.min_count.to_le_bytes());
        prefix[20..24].copy_from_slice(&self.max_count_observed.to_le_bytes());
        prefix[24..32].copy_from_slice(&n_kmers.to_le_bytes());
        prefix[32..40].copy_from_slice(&self.saturated_counts.to_le_bytes());
        prefix[40..44].copy_from_slice(&name_len.to_le_bytes());

        w.flush()?;
        w.seek(SeekFrom::Start(0))?;
        w.write_all(&prefix)?;
        w.write_all(&[0u8; 32])?;
        w.write_all(self.name.as_bytes())?;
        for o in self.key_offsets {
            w.write_all(&o.to_le_bytes())?;
        }
        for o in self.count_offsets {
            w.write_all(&o.to_le_bytes())?;
        }
        for c in self.bucket_counts {
            w.write_all(&c.to_le_bytes())?;
        }
        w.flush()?;

        let digest = hash_file_excluding_slot(&self.part, KTAB_DIGEST_START, KTAB_DIGEST_END)?;
        w.seek(SeekFrom::Start(KTAB_DIGEST_START as u64))?;
        w.write_all(&digest)?;
        w.flush()?;
        drop(self.writer.take());

        confine::replace_with(&self.part, &self.target)?;
        self.finished = true;
        Ok(())
    }
}

impl Drop for StreamingKtabWriter {
    fn drop(&mut self) {
        self.writer = None;
        self.counts_tmp_writer = None;
        if let Some(cp) = self.counts_tmp_path.take() {
            let _ = fs::remove_file(&cp);
        }
        if !self.finished && confine::regular_or_absent(&self.part).unwrap_or(false) {
            let _ = fs::remove_file(&self.part);
        }
    }
}

// ---------------------------------------------------------------------------
// Fast indexed parallel streaming .prof.bin builder
// ---------------------------------------------------------------------------

/// Compact sorted-key lookup index with an 18-bit (or 10-bit when `k < 9`)
/// sub-bucket directory for fast cache-friendly lookups.
struct CompactKmerIndex<K: KmerKey> {
    k: usize,
    sub_bits: u32,
    sub_offsets: Vec<usize>,
    keys: Vec<K>,
    /// Empty when every retained key has count 1 (presence-only).
    counts: Vec<u32>,
}

impl<K: KmerKey> CompactKmerIndex<K> {
    fn load_from_table(path: &Path, is_grk: bool) -> Result<Self> {
        let mut keys: Vec<K> = Vec::new();
        let mut counts: Vec<u32> = Vec::new();
        let k: usize;

        if is_grk {
            let reader = GrKsetReader::open(path)?;
            k = reader.k;
            let n = usize::try_from(reader.n_kmers).context("n_kmers exceeds addressable RAM")?;
            keys.try_reserve_exact(n)
                .context("profile index key allocation failed")?;
            let mut buf = Vec::new();
            for b in 0..BUCKETS {
                check_cancelled()?;
                reader.read_bucket(b, &mut buf, &mut keys)?;
            }
        } else {
            let reader = KtabReader::open(path)?;
            k = reader.k;
            let n = usize::try_from(reader.n_kmers).context("n_kmers exceeds addressable RAM")?;
            keys.try_reserve_exact(n)
                .context("profile index key allocation failed")?;
            let mut kbuf = Vec::new();
            let mut cbuf = Vec::new();
            if reader.presence_only {
                for b in 0..BUCKETS {
                    check_cancelled()?;
                    reader.read_bucket_keys(b, &mut kbuf, &mut keys)?;
                }
            } else {
                counts
                    .try_reserve_exact(n)
                    .context("profile index count allocation failed")?;
                for b in 0..BUCKETS {
                    check_cancelled()?;
                    reader.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
                }
            }
        }

        let sub_bits: u32 = if 2 * (k as u32) >= 18 {
            18
        } else {
            BUCKET_BITS
        };
        let n_sub = 1usize << sub_bits;
        let mut sub_offsets = vec![0usize; n_sub + 1];
        for &key in &keys {
            let idx = key.bucket_of(k, sub_bits);
            sub_offsets[idx + 1] += 1;
        }
        for i in 1..=n_sub {
            sub_offsets[i] += sub_offsets[i - 1];
        }

        Ok(Self {
            k,
            sub_bits,
            sub_offsets,
            keys,
            counts,
        })
    }

    #[inline]
    fn get(&self, mixed: K) -> u32 {
        let idx = mixed.bucket_of(self.k, self.sub_bits);
        let lo = self.sub_offsets[idx];
        let hi = self.sub_offsets[idx + 1];
        if lo == hi {
            return 0;
        }
        let slice = &self.keys[lo..hi];
        let found = if slice.len() <= 16 {
            let mut pos = None;
            for (i, &k) in slice.iter().enumerate() {
                if k >= mixed {
                    if k == mixed {
                        pos = Some(i);
                    }
                    break;
                }
            }
            pos
        } else {
            slice.binary_search(&mixed).ok()
        };
        match found {
            Some(rel) => {
                if self.counts.is_empty() {
                    1
                } else {
                    self.counts[lo + rel]
                }
            }
            None => 0,
        }
    }
}

/// Streaming O(1)-memory RLE + zigzag-delta ULEB128 encoder matching
/// `format::encode_profile_runs` byte for byte without allocating a `Vec<u32>`.
struct RunEncoder {
    prev: i64,
    cur_val: u32,
    cur_run: u64,
    max_val: u32,
}

impl RunEncoder {
    #[inline]
    fn new() -> Self {
        Self {
            prev: 0,
            cur_val: 0,
            cur_run: 0,
            max_val: 0,
        }
    }

    #[inline]
    fn push(&mut self, val: u32, out: &mut Vec<u8>) {
        if val > self.max_val {
            self.max_val = val;
        }
        if self.cur_run == 0 {
            self.cur_val = val;
            self.cur_run = 1;
        } else if val == self.cur_val {
            self.cur_run += 1;
        } else {
            self.flush_run(out);
            self.cur_val = val;
            self.cur_run = 1;
        }
    }

    #[inline]
    fn flush_run(&mut self, out: &mut Vec<u8>) {
        if self.cur_run == 0 {
            return;
        }
        let delta = (self.cur_val as i64) - self.prev;
        let zz = ((delta << 1) ^ (delta >> 63)) as u64;
        put_uleb(out, zz);
        put_uleb(out, self.cur_run - 1);
        self.prev = self.cur_val as i64;
        self.cur_run = 0;
    }

    #[inline]
    fn finish(mut self, out: &mut Vec<u8>) -> u32 {
        self.flush_run(out);
        self.max_val
    }
}

#[allow(clippy::too_many_arguments)]
fn encode_single_profile_record<K: KmerKey>(
    rec: &ProfileChunkRecord,
    chunk: &ProfileChunk,
    index: &CompactKmerIndex<K>,
    min_qual: u8,
    keep_record_ids: bool,
    space64: Option<KmerSpace64>,
    space128: Option<KmerSpace128>,
    enc_buf: &mut Vec<u8>,
    out_buf: &mut Vec<u8>,
) -> Result<(u64, u32)> {
    let k = index.k;
    let seq = &chunk.seq_data[rec.seq_start..rec.seq_start + rec.seq_len];
    let seq_len = seq.len() as u64;
    let qual = if rec.has_qual {
        Some(&chunk.qual_data[rec.qual_start..rec.qual_start + rec.seq_len])
    } else {
        None
    };

    enc_buf.clear();
    let (n_pos, max_val) = if seq.len() < k {
        (0u64, 0u32)
    } else {
        let n_pos = (seq.len() - k + 1) as u64;
        let qual_floor = 33u8.saturating_add(min_qual);
        let use_qual = min_qual > 0 && qual.is_some();
        let qslice = qual.unwrap_or(&[]);
        let shift = 2 * (k - 1);
        let mut runner = RunEncoder::new();

        if let Some(sp) = space64 {
            let mask = sp.mask();
            let (mut fwd, mut rev, mut filled) = (0u64, 0u64, 0usize);
            for (i, &b) in seq.iter().enumerate() {
                let c = CODE[b as usize];
                if c > 3 || (use_qual && qslice[i] < qual_floor) {
                    filled = 0;
                    if i + 1 >= k {
                        runner.push(0, enc_buf);
                    }
                    continue;
                }
                let c = c as u64;
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                filled += 1;
                if i + 1 >= k {
                    let count = if filled >= k {
                        let mixed = K::from_u128(sp.mix(fwd.min(rev)) as u128);
                        index.get(mixed)
                    } else {
                        0
                    };
                    runner.push(count, enc_buf);
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
                    if i + 1 >= k {
                        runner.push(0, enc_buf);
                    }
                    continue;
                }
                let c = c as u128;
                fwd = ((fwd << 2) | c) & mask;
                rev = (rev >> 2) | ((3 - c) << shift);
                filled += 1;
                if i + 1 >= k {
                    let count = if filled >= k {
                        let mixed = K::from_u128(sp.mix(fwd.min(rev)));
                        index.get(mixed)
                    } else {
                        0
                    };
                    runner.push(count, enc_buf);
                }
            }
        }
        (n_pos, runner.finish(enc_buf))
    };

    let anon_label;
    let name_bytes: &[u8] = if keep_record_ids {
        let raw = &chunk.name_data[rec.name_start..rec.name_start + rec.name_len];
        if std::str::from_utf8(raw).is_ok() {
            raw
        } else {
            anon_label = String::from_utf8_lossy(raw).into_owned();
            anon_label.as_bytes()
        }
    } else {
        anon_label = format!("r{:06}", rec.global_idx + 1);
        anon_label.as_bytes()
    };

    let name_len = u32::try_from(name_bytes.len()).context("profile sequence name exceeds u32")?;
    out_buf.extend_from_slice(&name_len.to_le_bytes());
    out_buf.extend_from_slice(name_bytes);
    out_buf.extend_from_slice(&seq_len.to_le_bytes());
    out_buf.extend_from_slice(&n_pos.to_le_bytes());
    out_buf.extend_from_slice(&(enc_buf.len() as u64).to_le_bytes());
    out_buf.extend_from_slice(enc_buf);

    Ok((seq_len, max_val))
}

#[allow(clippy::too_many_arguments)]
fn write_profiles_indexed_streaming<K: KmerKey>(
    inputs: &[PathBuf],
    table_path: &Path,
    table_is_grk: bool,
    prof_path: &Path,
    k: usize,
    min_qual: u8,
    keep_record_ids: bool,
    budget_bytes: usize,
) -> Result<()> {
    if let Some(parent) = prof_path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    confine::regular_or_absent(prof_path)?;
    let part = confine::suffixed(prof_path, ".part");

    let index = CompactKmerIndex::<K>::load_from_table(table_path, table_is_grk)?;
    let space64 = (k <= MAX_K_U64).then(|| KmerSpace64::new(k));
    let space128 = (k > MAX_K_U64).then(|| KmerSpace128::new(k));
    let n_threads = rayon::current_num_threads().max(1);
    let chunk_target_bytes = (budget_bytes / 16).clamp(512 * 1024, 4 * 1024 * 1024);

    let file = confine::create_part(&part)?;
    let mut writer = BufWriter::with_capacity(128 * 1024, file);
    let result = (|| -> Result<()> {
        writer.write_all(&[0u8; PROF_HEADER_LEN])?;

        let mut total_sequences = 0u64;
        let mut total_bases = 0u64;
        let mut max_count_observed = 0u32;
        let mut chunk = ProfileChunk::default();

        for path in inputs {
            let mut reader = ProfileSeqReader::open(path, keep_record_ids, total_sequences)?;
            while reader.fill_chunk(&mut chunk, chunk_target_bytes)? {
                check_cancelled()?;
                let n_recs = chunk.records.len();
                let group_size = n_recs.div_ceil(n_threads).max(1);
                let chunk_ref = &chunk;

                let parts: Vec<(Vec<u8>, u64, u32)> = chunk
                    .records
                    .par_chunks(group_size)
                    .map(|slice| -> Result<(Vec<u8>, u64, u32)> {
                        let mut enc_buf = Vec::with_capacity(512);
                        let mut out_buf = Vec::with_capacity(slice.len() * 48);
                        let mut sub_bases = 0u64;
                        let mut sub_max = 0u32;
                        for rec in slice {
                            let (sl, mv) = encode_single_profile_record(
                                rec,
                                chunk_ref,
                                &index,
                                min_qual,
                                keep_record_ids,
                                space64,
                                space128,
                                &mut enc_buf,
                                &mut out_buf,
                            )?;
                            sub_bases = sub_bases
                                .checked_add(sl)
                                .context("profile total_bases overflow")?;
                            if mv > sub_max {
                                sub_max = mv;
                            }
                        }
                        Ok((out_buf, sub_bases, sub_max))
                    })
                    .collect::<Result<Vec<_>>>()?;

                for (out_buf, sub_bases, sub_max) in parts {
                    writer.write_all(&out_buf)?;
                    total_bases = total_bases
                        .checked_add(sub_bases)
                        .context("profile total_bases overflow")?;
                    if sub_max > max_count_observed {
                        max_count_observed = sub_max;
                    }
                }
                total_sequences += n_recs as u64;
            }
        }

        let mut prefix = [0u8; PROF_DIGEST_START];
        prefix[0..8].copy_from_slice(MAGIC_PROF1);
        prefix[8] = k as u8;
        prefix[12..16].copy_from_slice(&max_count_observed.to_le_bytes());
        prefix[16..24].copy_from_slice(&total_sequences.to_le_bytes());
        prefix[24..32].copy_from_slice(&total_bases.to_le_bytes());

        writer.flush()?;
        writer.seek(SeekFrom::Start(0))?;
        writer.write_all(&prefix)?;
        writer.flush()?;

        let digest = hash_file_excluding_slot(&part, PROF_DIGEST_START, PROF_DIGEST_END)?;
        writer.seek(SeekFrom::Start(PROF_DIGEST_START as u64))?;
        writer.write_all(&digest)?;
        writer.flush()?;
        Ok(())
    })();
    drop(writer);
    match result {
        Ok(()) => confine::replace_with(&part, prof_path),
        Err(e) => {
            if confine::regular_or_absent(&part).unwrap_or(false) {
                let _ = fs::remove_file(&part);
            }
            Err(e)
        }
    }
}

// ---------------------------------------------------------------------------
// Public entry point: count_inputs
// ---------------------------------------------------------------------------

/// Runs the bounded-memory parallel k-mer counting engine over one or more
/// FASTA/FASTQ (plain or `.gz`) files, dispatching to the fast `u64` path
/// (`k <= 32`) or full `u128` path (`k = 33..=64`).
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
    ensure!(
        cfg.partitions.is_power_of_two() && (16..=1024).contains(&cfg.partitions),
        "partitions must be a power of 2 in 16..=1024, got {}",
        cfg.partitions
    );
    let m = cfg
        .minimizer_len
        .unwrap_or_else(|| default_minimizer_len(cfg.k));
    ensure!(
        (1..=cfg.k.min(31)).contains(&m),
        "minimizer length m={m} invalid for k={} (must satisfy 1 <= m <= min(k, 31))",
        cfg.k
    );

    let primary_fmt = if cfg.presence_only_out {
        TableFormat::GrKset
    } else {
        TableFormat::Ktab
    };
    check_output_extension(&cfg.out_path, primary_fmt)?;
    if let Some(gr_path) = &cfg.gr_kset_out {
        check_output_extension(gr_path, TableFormat::GrKset)?;
    }

    let _sig_guard = SignalGuard::install();
    let spill_parent = cfg.tmp_dir.clone().unwrap_or_else(std::env::temp_dir);
    let _ = cleanup_stale_spill_dirs(&spill_parent);

    if cfg.k <= MAX_K_U64 {
        count_inputs_generic::<u64>(inputs, cfg, m)
    } else {
        count_inputs_generic::<u128>(inputs, cfg, m)
    }
}

fn resolve_modality(
    choice: InputTypeChoice,
    saw_fasta: bool,
    saw_fastq: bool,
) -> Result<InputModality> {
    if saw_fasta && saw_fastq {
        bail!("mixed FASTA and FASTQ inputs are not supported in a single count run");
    }
    match choice {
        InputTypeChoice::Auto => {
            if saw_fastq {
                Ok(InputModality::Reads)
            } else if saw_fasta {
                Ok(InputModality::Assembly)
            } else {
                Ok(InputModality::Unknown)
            }
        }
        InputTypeChoice::Reads => {
            ensure!(
                !saw_fasta,
                "input contains FASTA records, which conflicts with --input-type reads"
            );
            Ok(InputModality::Reads)
        }
        InputTypeChoice::Assembly => {
            ensure!(
                !saw_fastq,
                "input contains FASTQ records, which conflicts with --input-type assembly"
            );
            Ok(InputModality::Assembly)
        }
    }
}

fn count_inputs_generic<K: KmerKey>(
    inputs: &[PathBuf],
    cfg: &CountConfig,
    m: usize,
) -> Result<CountSummary> {
    let budget_bytes = cfg.max_memory_mb.max(1).saturating_mul(1024 * 1024);
    let n_threads = rayon::current_num_threads().max(1);

    // -----------------------------------------------------------------------
    // Stage 1: Parallel record-aligned ingest & minimizer partitioning
    // -----------------------------------------------------------------------
    let ingest_budget = (budget_bytes / 3).max(256 * 1024);
    let batch_raw_cap = (budget_bytes / 8).clamp(512 * 1024, 32 * 1024 * 1024);
    let w_ingest = n_threads.min((batch_raw_cap / (256 * 1024)).max(1));
    let chunk_target_bytes = (batch_raw_cap / w_ingest).clamp(256 * 1024, 4 * 1024 * 1024);

    let mut store = PartitionStore::new(cfg.partitions, ingest_budget, cfg.tmp_dir.clone());
    let mut chunks: Vec<IngestChunk> = (0..w_ingest).map(|_| IngestChunk::default()).collect();
    let mut workers: Vec<WorkerExtractBuf> = (0..w_ingest)
        .map(|_| WorkerExtractBuf::new(cfg.k, m, cfg.partitions))
        .collect();

    let mut saw_fasta = false;
    let mut saw_fastq = false;
    let mut next_split_id = 0u64;
    let mut split_acc: FxHashMap<u64, (u64, u64)> = FxHashMap::default();
    let mut split_raw_lens = LenHistogram::new();
    let mut split_usable_lens = LenHistogram::new();
    let mut split_completed_reads = 0u64;

    for path in inputs {
        let mut reader = ChunkedSeqReader::open(path, cfg.k, next_split_id)?;
        loop {
            check_cancelled()?;
            let mut filled_chunks = 0usize;
            for ch in chunks.iter_mut().take(w_ingest) {
                if reader.fill_chunk(ch, chunk_target_bytes)? {
                    filled_chunks += 1;
                } else {
                    break;
                }
            }
            if filled_chunks == 0 {
                break;
            }
            let min_qual = cfg.min_qual;
            if filled_chunks == 1 {
                workers[0].reset_batch();
                workers[0].process_chunk(&chunks[0], min_qual);
            } else {
                workers[..filled_chunks]
                    .par_iter_mut()
                    .zip(chunks[..filled_chunks].par_iter())
                    .for_each(|(w, ch)| {
                        w.reset_batch();
                        w.process_chunk(ch, min_qual);
                    });
            }
            for w in &mut workers[..filled_chunks] {
                store.absorb_worker(&mut w.part_bufs)?;
                for &(cid, raw_add, valid_add, is_final) in &w.split_updates {
                    let entry = split_acc.entry(cid).or_insert((0, 0));
                    entry.0 += raw_add;
                    entry.1 += valid_add;
                    if is_final {
                        let (tot_raw, tot_valid) = split_acc.remove(&cid).unwrap();
                        split_completed_reads += 1;
                        split_raw_lens.add(tot_raw);
                        split_usable_lens.add(tot_valid);
                    }
                }
            }
        }
        saw_fasta |= reader.saw_fasta;
        saw_fastq |= reader.saw_fastq;
        next_split_id = reader.next_split_id;
    }

    let modality = resolve_modality(cfg.input_type, saw_fasta, saw_fastq)?;

    let mut total_reads = split_completed_reads;
    let mut total_bases = 0u64;
    let mut valid_bases = 0u64;
    let mut gc_bases = 0u64;
    let mut masked_bases = 0u64;
    let mut raw_lens = split_raw_lens;
    let mut usable_lens = split_usable_lens;
    for w in &workers {
        total_reads += w.total_reads;
        total_bases += w.total_bases;
        valid_bases += w.valid_bases;
        gc_bases += w.gc_bases;
        masked_bases += w.masked_bases;
        raw_lens.merge_from(&w.raw_lens);
        usable_lens.merge_from(&w.usable_lens);
    }
    drop(chunks);
    drop(workers);

    // -----------------------------------------------------------------------
    // Decide bounded streaming path vs fast in-memory path for Stages 2 & 3
    // -----------------------------------------------------------------------
    let sort_bytes_per_key = 2 * std::mem::size_of::<K>();
    let max_part_bytes = store.part_total_bytes.iter().copied().max().unwrap_or(0);
    let est_distinct_upper = valid_bases.min(store.total_bytes.saturating_mul(4));
    let est_in_mem_bytes = est_distinct_upper
        .saturating_mul(std::mem::size_of::<(K, u32)>() as u64 + 4)
        .saturating_add(
            max_part_bytes
                .saturating_mul(4)
                .saturating_mul(sort_bytes_per_key as u64),
        );
    let use_sharded_streaming =
        store.spill_dir.is_some() || est_in_mem_bytes > (budget_bytes / 2) as u64;

    if use_sharded_streaming {
        store.spill_all()?;
        if store.spill_dir.is_none() {
            store.spill_dir = Some(SpillDir::create(cfg.tmp_dir.as_deref())?);
        }
    }

    let spilled_to_disk = store.spill_dir.is_some();
    let spill_dir = store.spill_dir.take();
    let spilled_flags = std::mem::take(&mut store.spilled_flags);
    let mem_buffers = std::mem::take(&mut store.buffers);
    drop(store);

    // -----------------------------------------------------------------------
    // Stage 2: Bounded-concurrency partition counting
    // -----------------------------------------------------------------------
    let space = K::new_space(cfg.k);
    let shard_buf_budget = (budget_bytes / 5).clamp(512 * 1024, 32 * 1024 * 1024);
    let stage2_work_budget = if use_sharded_streaming {
        (budget_bytes / 3).max(4 * 1024 * 1024)
    } else {
        (budget_bytes / 2).max(4 * 1024 * 1024)
    };
    let w_stage2 = if budget_bytes <= 64 * 1024 * 1024 {
        n_threads.min(2)
    } else {
        n_threads.min((stage2_work_budget / (16 * 1024 * 1024)).max(1))
    };
    let per_worker_budget = (stage2_work_budget / w_stage2).max(2 * 1024 * 1024);
    let max_keys_per_chunk = (per_worker_budget / sort_bytes_per_key).max(65_536);
    let subrun_mem_cap = (per_worker_budget / 2).max(1024 * 1024);

    let sink: Stage2Sink<'_, K> = if use_sharded_streaming {
        Stage2Sink::Sharded(Box::new(Mutex::new(ShardWriterPool::new(
            spill_dir.as_ref().unwrap(),
            shard_buf_budget,
        ))))
    } else {
        Stage2Sink::InMemory(Mutex::new((0..BUCKETS).map(|_| Vec::new()).collect()))
    };

    type PartitionTask = (usize, Option<PathBuf>, Vec<u8>);
    let part_tasks: Vec<PartitionTask> = mem_buffers
        .into_iter()
        .enumerate()
        .map(|(p, mem_buf)| {
            let sp = if spilled_flags[p] {
                spill_dir.as_ref().map(|d| d.partition_path(p))
            } else {
                None
            };
            (p, sp, mem_buf)
        })
        .collect();

    let combined_stats = if w_stage2 <= 1 {
        let mut worker = Stage2WorkerScratch::<K>::new();
        let mut total_stats = PartitionStats::new();
        for (p, sp, mem_buf) in part_tasks {
            let st = process_partition_bounded(
                p,
                &space,
                cfg.k,
                cfg.min_count,
                cfg.max_count,
                cfg.allow_saturation,
                max_keys_per_chunk,
                subrun_mem_cap,
                spill_dir.as_ref(),
                sp,
                mem_buf,
                &mut worker,
                &sink,
            )?;
            total_stats.merge_from(st);
        }
        total_stats
    } else {
        let task_slots: Vec<Mutex<Option<PartitionTask>>> = part_tasks
            .into_iter()
            .map(|t| Mutex::new(Some(t)))
            .collect();
        let next_idx = AtomicUsize::new(0);
        let worker_results: Vec<Result<PartitionStats>> = std::thread::scope(|scope| {
            let mut handles = Vec::with_capacity(w_stage2);
            for _ in 0..w_stage2 {
                let slots_ref = &task_slots;
                let next_ref = &next_idx;
                let sink_ref = &sink;
                let spill_ref = spill_dir.as_ref();
                handles.push(scope.spawn(move || -> Result<PartitionStats> {
                    let mut worker = Stage2WorkerScratch::<K>::new();
                    let mut local_stats = PartitionStats::new();
                    loop {
                        let idx = next_ref.fetch_add(1, Ordering::Relaxed);
                        if idx >= slots_ref.len() {
                            break;
                        }
                        let task = slots_ref[idx]
                            .lock()
                            .unwrap_or_else(|poisoned| poisoned.into_inner())
                            .take()
                            .unwrap();
                        let (p, sp, mem_buf) = task;
                        let st = process_partition_bounded(
                            p,
                            &space,
                            cfg.k,
                            cfg.min_count,
                            cfg.max_count,
                            cfg.allow_saturation,
                            max_keys_per_chunk,
                            subrun_mem_cap,
                            spill_ref,
                            sp,
                            mem_buf,
                            &mut worker,
                            sink_ref,
                        )?;
                        local_stats.merge_from(st);
                    }
                    Ok(local_stats)
                }));
            }
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .unwrap_or_else(|_| Err(anyhow::anyhow!("Stage 2 worker panicked")))
                })
                .collect()
        });
        let mut total_stats = PartitionStats::new();
        for res in worker_results {
            total_stats.merge_from(res?);
        }
        total_stats
    };

    let (shard_written_flags, mut in_mem_buckets) = match sink {
        Stage2Sink::Sharded(pool_mtx) => {
            let mut pool = pool_mtx
                .into_inner()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            pool.flush_all()?;
            (Some(pool.written_flags), None)
        }
        Stage2Sink::InMemory(buckets_mtx) => {
            let buckets = buckets_mtx
                .into_inner()
                .unwrap_or_else(|poisoned| poisoned.into_inner());
            (None, Some(buckets))
        }
    };

    let distinct_total = combined_stats.distinct_total;
    let distinct_retained = combined_stats.distinct_retained;
    let total_occurrences = combined_stats.total_occurrences;
    let max_obs = combined_stats.max_obs;
    let max_count_retained = combined_stats.max_retained;
    let saturated_total = combined_stats.saturated_total;
    let saturated_retained = combined_stats.saturated_retained;
    let presence_only = distinct_retained == 0 || combined_stats.all_retained_one;

    let mut hist_bins: Vec<(u32, u64)> = combined_stats.hist_counts.into_iter().collect();
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
        .max_by_key(|&&(c, n)| (n, Reverse(c)))
        .map(|&(c, _)| c)
        .unwrap_or(if f1 > 0 { 1 } else { 0 });

    // -----------------------------------------------------------------------
    // Stage 3: Stream final .ktab / .grk one super-bucket (4 buckets) at a time
    // -----------------------------------------------------------------------
    let counts_spill_path = if !cfg.presence_only_out && !presence_only {
        spill_dir.as_ref().map(|d| d.counts_tmp_path())
    } else {
        None
    };

    let mut primary_ktab = if !cfg.presence_only_out {
        Some(StreamingKtabWriter::create(
            &cfg.out_path,
            &cfg.sample_name,
            cfg.k,
            modality,
            cfg.min_count,
            max_count_retained,
            saturated_retained,
            presence_only,
            counts_spill_path,
        )?)
    } else {
        None
    };
    let mut primary_grk = if cfg.presence_only_out {
        Some(StreamingGrkWriter::create(
            &cfg.out_path,
            &cfg.sample_name,
            cfg.k,
            modality,
        )?)
    } else {
        None
    };
    let mut sidecar_grk = match &cfg.gr_kset_out {
        Some(p) => Some(StreamingGrkWriter::create(
            p,
            &cfg.sample_name,
            cfg.k,
            modality,
        )?),
        None => None,
    };

    let mut fnv = Fnv::new();
    let mut keys_sha = Sha256Context::new(&SHA256);
    let mut counts_sha = Sha256Context::new(&SHA256);
    let mut sketch_hex = Vec::with_capacity(cfg.sketch_size);

    let mut keys_buf: Vec<K> = Vec::new();
    let mut counts_buf: Vec<u32> = Vec::new();
    let mut key_payload_buf: Vec<u8> = Vec::new();
    let mut count_payload_buf: Vec<u8> = Vec::new();
    let mut shard_read_buf = vec![0u8; 64 * 1024];
    let entry_bytes = K::WIDTH_BYTES + 4;

    for sb in 0..SUPER_BUCKETS {
        check_cancelled()?;
        let base_bucket = sb * BUCKETS_PER_SUPER;
        let mut sub_buckets: [Vec<(K, u32)>; BUCKETS_PER_SUPER] =
            [Vec::new(), Vec::new(), Vec::new(), Vec::new()];

        if let Some(flags) = shard_written_flags {
            if flags[sb] {
                let path = spill_dir.as_ref().unwrap().shard_path(sb);
                let mut f = File::open(&path)
                    .with_context(|| format!("read output shard {}", path.display()))?;
                let mut carry = [0u8; 20];
                let mut carry_len = 0usize;
                loop {
                    let n = f.read(&mut shard_read_buf)?;
                    if n == 0 {
                        ensure!(carry_len == 0, "truncated output shard {}", path.display());
                        break;
                    }
                    let mut slice = &shard_read_buf[..n];
                    if carry_len > 0 {
                        let need = entry_bytes - carry_len;
                        let take = need.min(slice.len());
                        carry[carry_len..carry_len + take].copy_from_slice(&slice[..take]);
                        carry_len += take;
                        slice = &slice[take..];
                        if carry_len == entry_bytes {
                            let key = decode_le_key::<K>(&carry[..K::WIDTH_BYTES]);
                            let count = u32::from_le_bytes(
                                carry[K::WIDTH_BYTES..entry_bytes].try_into().unwrap(),
                            );
                            let rel = key.bucket_of(cfg.k, BUCKET_BITS) - base_bucket;
                            sub_buckets[rel].push((key, count));
                            carry_len = 0;
                        }
                    }
                    let full_len = (slice.len() / entry_bytes) * entry_bytes;
                    let mut pos = 0usize;
                    while pos < full_len {
                        let rec = &slice[pos..pos + entry_bytes];
                        let key = decode_le_key::<K>(&rec[..K::WIDTH_BYTES]);
                        let count = u32::from_le_bytes(
                            rec[K::WIDTH_BYTES..entry_bytes].try_into().unwrap(),
                        );
                        let rel = key.bucket_of(cfg.k, BUCKET_BITS) - base_bucket;
                        sub_buckets[rel].push((key, count));
                        pos += entry_bytes;
                    }
                    let rem = &slice[full_len..];
                    if !rem.is_empty() {
                        carry[..rem.len()].copy_from_slice(rem);
                        carry_len = rem.len();
                    }
                }
                drop(f);
                let _ = fs::remove_file(&path);
            }
        } else if let Some(ref mut all_buckets) = in_mem_buckets {
            for r in 0..BUCKETS_PER_SUPER {
                sub_buckets[r] = std::mem::take(&mut all_buckets[base_bucket + r]);
            }
        }

        let sb_total: usize = sub_buckets.iter().map(Vec::len).sum();
        if w_stage2 > 1 && sb_total >= 4096 {
            sub_buckets.par_iter_mut().for_each(|entries| {
                entries.sort_unstable_by_key(|&(key, _)| key);
            });
        } else {
            for entries in &mut sub_buckets {
                entries.sort_unstable_by_key(|&(key, _)| key);
            }
        }

        for (r, entries) in sub_buckets.into_iter().enumerate() {
            let b = base_bucket + r;
            debug_assert!(entries.windows(2).all(|w| w[0].0 < w[1].0));
            keys_buf.clear();
            counts_buf.clear();
            keys_buf.reserve(entries.len());
            counts_buf.reserve(entries.len());

            for (key, count) in entries {
                key.with_le_bytes(|kb| {
                    fnv.bytes(kb);
                    keys_sha.update(kb);
                    counts_sha.update(kb);
                });
                counts_sha.update(&count.to_le_bytes());
                if sketch_hex.len() < cfg.sketch_size {
                    sketch_hex.push(key.to_hex());
                }
                keys_buf.push(key);
                counts_buf.push(count);
            }

            let b_count = u32::try_from(keys_buf.len()).context("bucket count exceeds u32")?;
            key_payload_buf.clear();
            encode_bucket_keys(cfg.k, b, &keys_buf, &mut key_payload_buf);

            if let Some(ref mut kw) = primary_ktab {
                count_payload_buf.clear();
                if !presence_only && !counts_buf.is_empty() {
                    encode_bucket_counts(&counts_buf, &mut count_payload_buf);
                }
                kw.write_bucket(b, b_count, &key_payload_buf, &count_payload_buf)?;
            }
            if let Some(ref mut gw) = primary_grk {
                gw.write_bucket(b, b_count, &key_payload_buf)?;
            }
            if let Some(ref mut sg) = sidecar_grk {
                sg.write_bucket(b, b_count, &key_payload_buf)?;
            }
        }
    }

    drop(in_mem_buckets);
    drop(keys_buf);
    drop(counts_buf);
    drop(key_payload_buf);
    drop(count_payload_buf);
    drop(shard_read_buf);

    if let Some(kw) = primary_ktab {
        kw.finish()?;
    }
    if let Some(gw) = primary_grk {
        gw.finish()?;
    }
    if let Some(sg) = sidecar_grk {
        sg.finish()?;
    }
    drop(spill_dir);

    // -----------------------------------------------------------------------
    // Optional indexed parallel streaming .prof.bin pass
    // -----------------------------------------------------------------------
    if let Some(prof_path) = &cfg.profile_out {
        write_profiles_indexed_streaming::<K>(
            inputs,
            &cfg.out_path,
            cfg.presence_only_out,
            prof_path,
            cfg.k,
            cfg.min_qual,
            cfg.keep_record_ids,
            budget_bytes,
        )?;
    }

    let qc = finalize_qc(
        total_reads,
        total_bases,
        valid_bases,
        gc_bases,
        masked_bases,
        &usable_lens,
        &raw_lens,
        modality == InputModality::Reads,
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
        saturated_total,
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
        modality,
        min_count: cfg.min_count,
        max_count: cfg.max_count,
        distinct_retained,
        max_count_retained,
        saturated_counts: saturated_retained,
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

#[allow(clippy::too_many_arguments)]
fn finalize_qc(
    total_reads: u64,
    total_bases: u64,
    valid_bases: u64,
    gc_bases: u64,
    masked_bases: u64,
    usable_lens: &LenHistogram,
    raw_lens: &LenHistogram,
    is_fastq: bool,
    min_qual_applied: u8,
    singleton_fraction: f64,
) -> ReadQcMetrics {
    let usable_read_len_median = usable_lens.quantile(0.5);
    let usable_read_len_p10 = usable_lens.quantile(0.1);
    let usable_read_len_p90 = usable_lens.quantile(0.9);
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
    let n50 = raw_lens.n50(total_bases);

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

fn write_json_atomic<T: Serialize>(path: &Path, val: &T) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let mut bytes = serde_json::to_vec_pretty(val)?;
    bytes.push(b'\n');
    confine::write_atomic(path, &bytes, false)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kcount::format::{decode_profile_runs, write_gr_kset, write_ktab};

    #[test]
    fn len_histogram_matches_sorted_vec_quantiles_and_n50() {
        let samples = [
            0u64, 10, 10, 150, 150, 150, 200, 4095, 4096, 5000, 10_000, 25_000,
        ];
        let mut h = LenHistogram::new();
        for &s in &samples {
            h.add(s);
        }
        let mut sorted = samples.to_vec();
        sorted.sort_unstable();
        let total_bases: u64 = sorted.iter().sum();

        assert_eq!(h.quantile(0.5), (150.0 + 200.0) / 2.0);
        assert_eq!(
            h.quantile(0.1),
            sorted[(11.0 * 0.1f64).round() as usize] as f64
        );
        assert_eq!(
            h.quantile(0.9),
            sorted[(11.0 * 0.9f64).round() as usize] as f64
        );
        // 25_000 * 2 = 50_000 >= total_bases (48_861) -> N50 is 25_000.
        assert_eq!(total_bases, 48_861);
        assert_eq!(h.n50(total_bases), 25_000);

        // Second case where N50 lands in the second-largest entry (10_000) and
        // third case where N50 crosses from large down into small (< 4096).
        let mut h2 = LenHistogram::new();
        for &s in &[2_000u64, 5_000, 10_000, 10_000, 15_000] {
            h2.add(s);
        }
        assert_eq!(h2.n50(42_000), 10_000);

        let mut h3 = LenHistogram::new();
        for &s in &[1_000u64, 1_000, 1_000, 2_000, 2_000, 5_000] {
            h3.add(s);
        }
        assert_eq!(h3.n50(12_000), 2_000);

        let mut odd = LenHistogram::new();
        for &s in &[100u64, 200, 300] {
            odd.add(s);
        }
        assert_eq!(odd.quantile(0.5), 200.0);
        assert_eq!(odd.n50(600), 300);

        let empty = LenHistogram::new();
        assert_eq!(empty.quantile(0.5), 0.0);
        assert_eq!(empty.n50(0), 0);
    }

    #[test]
    fn count_saturation_refuses_by_default_and_clamps_when_allowed() {
        assert_eq!(
            finalize_u32_count(u32::MAX as u64, false).unwrap(),
            (u32::MAX, false)
        );
        let err = finalize_u32_count(u32::MAX as u64 + 1, false).unwrap_err();
        assert!(format!("{err:#}").contains("allow-saturation"));
        assert_eq!(
            finalize_u32_count(u32::MAX as u64 + 1, true).unwrap(),
            (u32::MAX, true)
        );

        let runs: Vec<Vec<(u64, u64)>> = vec![
            vec![(10, 3_000_000_000), (20, 5)],
            vec![(10, 2_000_000_000), (30, 7)],
        ];
        assert!(merge_counted_runs(&runs, false).is_err());
        let (merged, sat) = merge_counted_runs(&runs, true).unwrap();
        assert_eq!(sat, 1);
        assert_eq!(merged, vec![(10, u32::MAX), (20, 5), (30, 7)]);
    }

    #[test]
    fn streaming_run_encoder_matches_profile_decoder_roundtrip() {
        let cases: &[&[u32]] = &[
            &[],
            &[0],
            &[0, 0, 0, 0],
            &[1, 1, 2, 2, 2, 0, 0, 300, 65_535, 65_536, u32::MAX, 1, 1],
        ];
        for &seq in cases {
            let mut enc = Vec::new();
            let mut r = RunEncoder::new();
            for &v in seq {
                r.push(v, &mut enc);
            }
            let max_v = r.finish(&mut enc);
            assert_eq!(max_v, seq.iter().copied().max().unwrap_or(0));
            let dec = decode_profile_runs(&enc, seq.len()).unwrap();
            assert_eq!(dec, seq);
        }
    }

    #[test]
    fn streaming_writers_match_in_memory_format_writers_byte_for_byte() {
        let dir = tempfile::tempdir().unwrap();
        let k = 31usize;
        let entries: Vec<(u64, u32)> = vec![
            (10, 1),
            (25, 254),
            (100, 65_535),
            (1u64 << 52, 65_536),
            ((1u64 << 52) + 77, u32::MAX),
            ((1023u64 << 52) + 9, 3),
        ];
        let keys: Vec<u64> = entries.iter().map(|&(k, _)| k).collect();

        let ref_grk = dir.path().join("ref.grk");
        write_gr_kset(&ref_grk, "sample", k, InputModality::Assembly, &keys).unwrap();
        let stream_grk = dir.path().join("stream.grk");
        let mut gw =
            StreamingGrkWriter::create(&stream_grk, "sample", k, InputModality::Assembly).unwrap();

        let ref_ktab = dir.path().join("ref.ktab");
        write_ktab(&ref_ktab, "sample", k, InputModality::Reads, 1, 1, &entries).unwrap();
        let stream_ktab = dir.path().join("stream.ktab");
        let mut kw = StreamingKtabWriter::create(
            &stream_ktab,
            "sample",
            k,
            InputModality::Reads,
            1,
            u32::MAX,
            1,
            false,
            Some(dir.path().join("counts.tmp")),
        )
        .unwrap();

        let mut i = 0usize;
        let mut kb = Vec::new();
        let mut cb = Vec::new();
        let mut kp = Vec::new();
        let mut cp = Vec::new();
        for b in 0..BUCKETS {
            kb.clear();
            cb.clear();
            while i < entries.len() && entries[i].0.bucket_of(k, BUCKET_BITS) == b {
                kb.push(entries[i].0);
                cb.push(entries[i].1);
                i += 1;
            }
            kp.clear();
            cp.clear();
            encode_bucket_keys(k, b, &kb, &mut kp);
            encode_bucket_counts(&cb, &mut cp);
            gw.write_bucket(b, kb.len() as u32, &kp).unwrap();
            kw.write_bucket(b, kb.len() as u32, &kp, &cp).unwrap();
        }
        gw.finish().unwrap();
        kw.finish().unwrap();

        assert_eq!(fs::read(&ref_grk).unwrap(), fs::read(&stream_grk).unwrap());
        assert_eq!(
            fs::read(&ref_ktab).unwrap(),
            fs::read(&stream_ktab).unwrap()
        );
    }

    #[test]
    fn spill_dir_name_parser_validates_strict_format() {
        assert_eq!(parse_spill_dir_pid("grm-kcount-12345-0"), Some(12345));
        assert_eq!(parse_spill_dir_pid("grm-kcount-0-1"), None);
        assert_eq!(parse_spill_dir_pid("grm-kcount-abc-1"), None);
        assert_eq!(parse_spill_dir_pid("grm-kcount-123-"), None);
        assert_eq!(parse_spill_dir_pid("grm-kcount-123-4-5"), None);
        assert_eq!(parse_spill_dir_pid("other-123-0"), None);
    }
}
