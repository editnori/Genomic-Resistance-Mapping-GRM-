//! Species-level presence matrix.
//!
//! Merging every genome's kset gives each distinct k-mer a presence bitset
//! over the genomes. Many k-mers share the same bitset (every k-mer along an
//! accessory gene, every k-mer overlapping one SNP), so k-mers collapse into
//! *patterns*: one row per distinct bitset, plus the list of k-mers that
//! carry it. Learners see patterns; annotation goes back to the k-mers.
//!
//! Files in a matrix directory:
//!   matrix.json        summary and genome order (column i = genomes[i])
//!   patterns.bin       n_patterns rows of `words` little-endian u64
//!   support.bin        u32 per pattern: genomes carrying it
//!   kmer_offsets.bin   u64 per pattern + 1: slice bounds into kmers.bin
//!   kmers.bin          mixed k-mer values grouped by pattern
//!
//! `build` merges every genome's kset at once. `write_partial` and
//! `build_from_partials` reach the same matrix one block of genomes at a
//! time, so the blocks' ksets need not all exist together; see
//! docs/matrix-scaling.md.
//!
//! Every entry point refuses more than `REVIEWED_COLUMNS` genomes: past
//! that, downstream stages have not been reviewed.

use crate::bits::{self, BitMatrix, Bits, words_for};
use crate::dna::{MAX_K, MIN_K};
use crate::fnv::Fnv;
use crate::kset::{BUCKET_BITS, BUCKETS, KsetReader, radix_sort};
use crate::layout::remove_if_present;
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::os::unix::fs::FileExt;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

/// The most genomes (matrix columns) that `build`, `write_partial`,
/// `build_from_partials` and `FeatureMatrix::load` accept. `Dataset::project`,
/// the learners, annotation and reports have only been checked up to here,
/// the packed merge key's ceiling at k = 31. Raising it needs a review of
/// those stages; until then only this module's tests merge past it, through
/// private functions that skip the check (docs/matrix-scaling.md).
pub const REVIEWED_COLUMNS: usize = 4096;

fn ensure_reviewed(n: usize, max_columns: usize) -> Result<()> {
    ensure!(
        n <= max_columns,
        "{n} genomes exceed the {max_columns} matrix columns that Dataset::project, the learners and reports have been reviewed for (docs/matrix-scaling.md)"
    );
    Ok(())
}

fn ensure_unique_names<'a>(names: impl IntoIterator<Item = &'a str>) -> Result<()> {
    let mut seen = FxHashSet::default();
    for name in names {
        ensure!(seen.insert(name), "genome {name} is in two columns");
    }
    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MatrixSummary {
    pub k: usize,
    pub genomes: Vec<String>,
    pub min_support: u32,
    pub n_patterns: usize,
    pub kmers_distinct: u64,
    pub kmers_kept: u64,
    /// Present in fewer than `min_support` genomes (singletons included).
    pub kmers_rare: u64,
    /// Present in more than n - min_support genomes (near-core).
    pub kmers_core: u64,
    #[serde(default)]
    pub kset_digests: Vec<String>,
}

pub struct FeatureMatrix {
    pub summary: MatrixSummary,
    pub patterns: BitMatrix,
    pub support: Vec<u32>,
    pub kmer_offsets: Vec<u64>,
    pub kmers: Vec<u64>,
}

impl FeatureMatrix {
    pub fn n_genomes(&self) -> usize {
        self.summary.genomes.len()
    }

    pub fn kmers_of(&self, pattern: usize) -> &[u64] {
        &self.kmers[self.kmer_offsets[pattern] as usize..self.kmer_offsets[pattern + 1] as usize]
    }
}

struct BucketOut {
    rows: Vec<u64>,
    pairs: Vec<(u32, u64)>,
    distinct: u64,
    rare: u64,
    core: u64,
}

/// Which k-mers a merge keeps, and how wide their rows are.
#[derive(Clone, Copy)]
struct Keep {
    words: usize,
    min: u32,
    max: u32,
}

impl Keep {
    /// `n` must already be checked to fit u32.
    fn new(n: usize, min_support: u32) -> Keep {
        Keep {
            words: words_for(n),
            min: min_support,
            max: (n as u32).saturating_sub(min_support),
        }
    }
}

/// Collects one bucket's k-mers, visited in ascending order, into patterns
/// numbered by first appearance.
struct BucketMerge {
    keep: Keep,
    local: FxHashMap<Box<[u64]>, u32>,
    out: BucketOut,
}

impl BucketMerge {
    fn new(keep: Keep) -> Self {
        BucketMerge {
            keep,
            local: FxHashMap::default(),
            out: BucketOut {
                rows: Vec::new(),
                pairs: Vec::new(),
                distinct: 0,
                rare: 0,
                core: 0,
            },
        }
    }

    /// Counts k-mer `value` present in the genomes of `row`, keeping it
    /// when min <= support <= max.
    #[inline]
    fn add(&mut self, value: u64, row: &[u64], support: u32) -> Result<()> {
        let out = &mut self.out;
        out.distinct += 1;
        if support < self.keep.min {
            out.rare += 1;
            return Ok(());
        }
        if support > self.keep.max {
            out.core += 1;
            return Ok(());
        }
        let next = u32::try_from(self.local.len()).context("too many patterns in one bucket")?;
        let id = *self.local.entry(row.into()).or_insert_with(|| {
            out.rows.extend_from_slice(row);
            next
        });
        out.pairs.push((id, value));
        Ok(())
    }

    fn finish(self) -> BucketOut {
        self.out
    }
}

/// How a bucket's (k-mer, genome) entries are keyed for grouping. Both keys
/// sort by k-mer first, so they find the same k-mers in the same order and
/// build byte-identical matrices; only cost and the supported cohort differ.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum MergeKey {
    /// The k-mer's bits below the bucket prefix and the genome index packed
    /// into one u64 and radix sorted. Needs
    /// 2k - BUCKET_BITS + ceil(log2 n) <= 64: 4096 genomes at k = 31, 1024
    /// at k = 32.
    Packed,
    /// (k-mer bits, genome index) as a (u64, u32) tuple, comparison sorted.
    /// The same 16 bytes per entry as packed plus its radix scratch, but
    /// O(T log T) per bucket instead of a few linear passes.
    Wide,
}

/// The packed key when it holds every k-mer bit and genome index, else the
/// wide one. Genome indices and support counts are u32, which bounds both.
fn merge_key(k: usize, n: usize) -> Result<MergeKey> {
    ensure!(
        n <= u32::MAX as usize,
        "{n} genomes exceed the u32 genome index of the merge"
    );
    Ok(if 2 * k as u32 - BUCKET_BITS + genome_bits(n) <= 64 {
        MergeKey::Packed
    } else {
        MergeKey::Wide
    })
}

/// What every bucket merge of one set of ksets shares.
struct MergePlan {
    key: MergeKey,
    k: usize,
    n: usize,
    keep: Keep,
}

/// Merges ksets (one per genome, in column order) into a pattern matrix.
/// A k-mer is kept when min_support <= support <= n - min_support. More
/// than `REVIEWED_COLUMNS` ksets are refused before any is opened.
pub fn build(
    kset_paths: &[PathBuf],
    min_support: u32,
    threads_chunk: usize,
) -> Result<FeatureMatrix> {
    ensure_reviewed(kset_paths.len(), REVIEWED_COLUMNS)?;
    build_with(kset_paths, min_support, threads_chunk, None)
}

/// Opens one kset per genome, all with the same k and distinct names.
fn open_ksets(kset_paths: &[PathBuf]) -> Result<(Vec<KsetReader>, usize)> {
    ensure!(!kset_paths.is_empty(), "no ksets to merge");
    let n = kset_paths.len();
    let readers: Vec<KsetReader> = kset_paths
        .iter()
        .map(|p| {
            KsetReader::open(p).map_err(|e| {
                with_open_file_hint(e, || {
                    format!(
                        "merging {n} ksets keeps one open per genome; raise the open-file limit (ulimit -n) above {n}, or merge in column blocks (write_partial)"
                    )
                })
            })
        })
        .collect::<Result<_>>()?;
    let k = readers[0].k;
    ensure!(
        readers.iter().all(|r| r.k == k),
        "ksets mix different k values"
    );
    ensure_unique_names(readers.iter().map(|r| r.name.as_str()))?;
    Ok((readers, k))
}

/// `build` without the column cap, optionally forcing the merge key. Only
/// `build` and tests call it. Forcing `Wide` on a cohort the packed key
/// supports lets tests compare the two; forcing `Packed` beyond its ceiling
/// is an error, never a truncated key.
fn build_with(
    kset_paths: &[PathBuf],
    min_support: u32,
    threads_chunk: usize,
    force: Option<MergeKey>,
) -> Result<FeatureMatrix> {
    let (readers, k) = open_ksets(kset_paths)?;
    let n = readers.len();
    let natural = merge_key(k, n)?;
    let key = match force {
        Some(MergeKey::Packed) => {
            ensure!(
                natural == MergeKey::Packed,
                "{n} genomes at k = {k} do not fit the packed merge key"
            );
            MergeKey::Packed
        }
        Some(key) => key,
        None => natural,
    };
    if natural == MergeKey::Wide {
        eprintln!(
            "  {n} genomes at k = {k} exceed the packed merge key; merging with (k-mer, genome) tuples"
        );
    }
    let plan = MergePlan {
        key,
        k,
        n,
        keep: Keep::new(n, min_support),
    };
    let genomes = readers.iter().map(|r| r.name.clone()).collect();
    assemble(k, genomes, min_support, threads_chunk, |b| {
        merge_bucket(&plan, &readers, b)
    })
}

/// Folds every bucket's patterns, in bucket order, into one matrix whose
/// columns are `genomes`. `merge(b)` merges bucket b.
fn assemble(
    k: usize,
    genomes: Vec<String>,
    min_support: u32,
    threads_chunk: usize,
    merge: impl Fn(usize) -> Result<BucketOut> + Sync,
) -> Result<FeatureMatrix> {
    let n = genomes.len();
    let words = words_for(n);
    let mut index: FxHashMap<Box<[u64]>, u32> = FxHashMap::default();
    let mut support: Vec<u32> = Vec::new();
    let mut pairs: Vec<(u32, u64)> = Vec::new();
    let (mut distinct, mut rare, mut core) = (0u64, 0u64, 0u64);

    let chunk = threads_chunk.max(1);
    for start in (0..BUCKETS).step_by(chunk) {
        let end = (start + chunk).min(BUCKETS);
        let outs: Vec<BucketOut> = (start..end)
            .into_par_iter()
            .map(&merge)
            .collect::<Result<_>>()?;
        for out in outs {
            distinct += out.distinct;
            rare += out.rare;
            core += out.core;
            let mut local_to_global = Vec::with_capacity(out.rows.len() / words.max(1));
            for row in out.rows.chunks_exact(words) {
                // Pattern ids, and `n_patterns as u32` below, stay exact.
                ensure!(
                    support.len() < u32::MAX as usize,
                    "too many patterns for u32 ids"
                );
                let next = support.len() as u32;
                let id = *index.entry(row.into()).or_insert_with(|| {
                    support.push(bits::count(row));
                    next
                });
                local_to_global.push(id);
            }
            pairs.extend(
                out.pairs
                    .into_iter()
                    .map(|(local, kmer)| (local_to_global[local as usize], kmer)),
            );
        }
        if end.is_multiple_of(128) || end == BUCKETS {
            eprintln!(
                "  merged {end} / {BUCKETS} buckets: {distinct} distinct k-mers, {} kept, {} patterns",
                pairs.len(),
                support.len()
            );
        }
    }

    let n_patterns = support.len();
    let mut patterns = BitMatrix {
        rows: n_patterns,
        cols: n,
        words,
        data: vec![0; n_patterns * words],
    };
    for (row, &id) in index.iter() {
        patterns.data[id as usize * words..(id as usize + 1) * words].copy_from_slice(row);
    }
    drop(index);

    pairs.par_sort_unstable();
    let mut kmer_offsets = Vec::with_capacity(n_patterns + 1);
    let mut kmers = Vec::with_capacity(pairs.len());
    let mut cursor = 0;
    for id in 0..n_patterns as u32 {
        kmer_offsets.push(kmers.len() as u64);
        while cursor < pairs.len() && pairs[cursor].0 == id {
            kmers.push(pairs[cursor].1);
            cursor += 1;
        }
    }
    kmer_offsets.push(kmers.len() as u64);

    let summary = MatrixSummary {
        k,
        genomes,
        min_support,
        n_patterns,
        kmers_distinct: distinct,
        kmers_kept: kmers.len() as u64,
        kmers_rare: rare,
        kmers_core: core,
        kset_digests: vec![],
    };
    Ok(FeatureMatrix {
        summary,
        patterns,
        support,
        kmer_offsets,
        kmers,
    })
}

/// A merge keeps every input open at once, so more inputs than the
/// open-file limit fail while opening; say so, with `advice`, rather than
/// only naming the file that happened to be next.
fn with_open_file_hint(err: anyhow::Error, advice: impl FnOnce() -> String) -> anyhow::Error {
    // ENFILE (system wide) and EMFILE (per process) on Linux.
    let exhausted = err
        .chain()
        .filter_map(|c| c.downcast_ref::<std::io::Error>())
        .any(|e| matches!(e.raw_os_error(), Some(23 | 24)));
    if exhausted {
        err.context(advice())
    } else {
        err
    }
}

/// Bits needed to hold a genome index below `n`.
fn genome_bits(n: usize) -> u32 {
    (usize::BITS - (n.max(2) - 1).leading_zeros()).max(1)
}

/// Merges one bucket across all genomes.
fn merge_bucket(plan: &MergePlan, readers: &[KsetReader], bucket: usize) -> Result<BucketOut> {
    let total: usize = readers.iter().map(|r| r.bucket_len(bucket)).sum();
    let mut bytes = Vec::new();
    let mut merge = BucketMerge::new(plan.keep);
    merge_entries(
        plan,
        bucket,
        total,
        |g, out| readers[g].read_bucket(bucket, &mut bytes, out),
        |value, row, support| merge.add(value, row, support),
    )?;
    Ok(merge.finish())
}

/// Merges one bucket of `total` values, which `fill(g, out)` appends for
/// genome g, and hands each distinct k-mer, in ascending order, to
/// `sink(value, presence row, support)`. Every value in a bucket shares its
/// top BUCKET_BITS bits, so only the remaining bits are keyed, next to the
/// genome index; one sort then groups each k-mer's genomes together.
fn merge_entries(
    plan: &MergePlan,
    bucket: usize,
    total: usize,
    mut fill: impl FnMut(usize, &mut Vec<u64>) -> Result<()>,
    sink: impl FnMut(u64, &[u64], u32) -> Result<()>,
) -> Result<()> {
    let low_bits = 2 * plan.k as u32 - BUCKET_BITS;
    let low_mask = (1u64 << low_bits) - 1;
    let base = (bucket as u64) << low_bits;
    let words = plan.keep.words;
    let mut values = Vec::new();
    match plan.key {
        MergeKey::Packed => {
            let gbits = genome_bits(plan.n);
            ensure!(
                low_bits + gbits <= 64,
                "{} genomes at k = {} do not fit the packed merge key",
                plan.n,
                plan.k
            );
            let mut packed: Vec<u64> = Vec::with_capacity(total);
            for g in 0..plan.n {
                values.clear();
                fill(g, &mut values)?;
                packed.extend(values.iter().map(|&v| ((v & low_mask) << gbits) | g as u64));
            }
            radix_sort(&mut packed, &mut values);
            let gmask = (1u64 << gbits) - 1;
            group_sorted(words, base, packed.len(), sink, |i| {
                (packed[i] >> gbits, (packed[i] & gmask) as usize)
            })
        }
        MergeKey::Wide => {
            // merge_key keeps n within u32, so the genome index is exact.
            let mut entries: Vec<(u64, u32)> = Vec::with_capacity(total);
            for g in 0..plan.n {
                values.clear();
                fill(g, &mut values)?;
                entries.extend(values.iter().map(|&v| (v & low_mask, g as u32)));
            }
            entries.sort_unstable();
            group_sorted(words, base, entries.len(), sink, |i| {
                (entries[i].0, entries[i].1 as usize)
            })
        }
    }
}

/// Walks `len` entries sorted by k-mer, `entry(i)` being the i-th
/// (k-mer bits, genome index), and hands each k-mer's row to `sink`.
/// Both keys also sort a k-mer's genomes, so a genome that lists a k-mer
/// twice shows up as a repeated neighbour; it is refused, because support
/// would then count more genomes than the row holds.
fn group_sorted(
    words: usize,
    base: u64,
    len: usize,
    mut sink: impl FnMut(u64, &[u64], u32) -> Result<()>,
    entry: impl Fn(usize) -> (u64, usize),
) -> Result<()> {
    let mut row = vec![0u64; words];
    let mut i = 0;
    while i < len {
        let (key, genome) = entry(i);
        row.fill(0);
        bits::set(&mut row, genome);
        let mut support = 1u32;
        let mut last = genome;
        i += 1;
        while i < len {
            let (next, genome) = entry(i);
            if next != key {
                break;
            }
            ensure!(
                genome != last,
                "genome {genome} lists k-mer {:#x} twice",
                base | key
            );
            bits::set(&mut row, genome);
            support += 1;
            last = genome;
            i += 1;
        }
        sink(base | key, &row, support)?;
    }
    Ok(())
}

// Column-block partials: the out-of-core merge.
//
// `build` holds every kset open, and the pipeline keeps every kset on disk
// until the merge. A partial instead stores one block of columns (genomes)
// as each distinct k-mer once with its presence bitset over the block, so a
// caller can count a block of genomes, write its partial, and delete those
// ksets before counting the next block. `build_from_partials` then merges
// the blocks k-mer by k-mer with one open file per block, and produces the
// same matrix, byte for byte, that `build` produces from all the ksets.
//
// Layout (little endian):
//   magic "GRMKPRT2" | k u8 | bucket_bits u8 | reserved [u8; 6], all zero
//   | first_col u64 | width u64 | n_kmers u64 | full_rows u64
//   | width x (name_len u32 | name | kset digest u64)
//   | bucket byte offsets [(B+1) x u64] | counts [B x u32]
//   | data: per bucket, per k-mer in ascending order, LEB128 of
//     (gap from the previous k-mer, or from the bucket floor) << 1 | full,
//     then, unless full, the row's words_for(width) u64 words.
// `full` marks a k-mer present in every genome of the block, which costs no
// row at all; most of a clonal species' k-mers are expected to be full.
//
// Each partial has a receipt, `<partial>.receipt.json`, written once the
// partial is durable. It repeats the header and holds the file's SHA-256.
// `build_from_partials` accepts a partial only when its file, header and
// receipt agree and its columns (genome names and kset digests) are the
// ones the caller expects, so a partial left from other ksets, an older
// catalog, or an interrupted rewrite is refused instead of merged.

const PART_MAGIC: &[u8; 8] = b"GRMKPRT2";
const PART_HEAD: usize = 48;
/// Longest genome name a partial holds.
const MAX_NAME: usize = 4096;
/// The row slot of a k-mer present in every genome of its block.
const FULL: u32 = u32::MAX;

/// One matrix column as the caller expects it: the genome, and the
/// `kset::digest` of the kset counted for it (as recorded in counted.json).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Column {
    pub name: String,
    pub kset_digest: String,
}

/// What `write_partial` wrote. It is also saved beside the partial (see
/// `receipt_path`), and `build_from_partials` checks that the two agree.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct PartialReceipt {
    pub k: usize,
    pub first_col: usize,
    /// The block's genomes and the digests of the ksets that were read.
    pub columns: Vec<Column>,
    /// Distinct k-mers in the block.
    pub kmers: u64,
    /// Of those, present in every genome of the block.
    pub full_rows: u64,
    pub bytes: u64,
    /// SHA-256 of the partial file.
    pub sha256: String,
}

/// Where `write_partial` saves the receipt of the partial at `partial`.
pub fn receipt_path(partial: &Path) -> PathBuf {
    with_suffix(partial, ".receipt.json")
}

/// `path` with `suffix` appended to its file name, so that x.a and x.b keep
/// distinct temporary names (`with_extension` would map both to x.<ext>).
fn with_suffix(path: &Path, suffix: &str) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(suffix);
    PathBuf::from(name)
}

fn parent_dir(path: &Path) -> &Path {
    match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p,
        _ => Path::new("."),
    }
}

/// Makes renames and removals in the directory of `path` durable.
fn sync_parent(path: &Path) -> Result<()> {
    let dir = parent_dir(path);
    File::open(dir)
        .and_then(|d| d.sync_all())
        .with_context(|| format!("sync directory {}", dir.display()))
}

/// A temporary file that is removed when dropped unless kept, so that an
/// error at any point after creating it leaves nothing behind.
struct TempFile {
    path: PathBuf,
    keep: bool,
}

impl Drop for TempFile {
    fn drop(&mut self) {
        if !self.keep {
            let _ = std::fs::remove_file(&self.path);
        }
    }
}

/// Writes `bytes` to `path` through `<path>.part`: synced, renamed, and
/// the directory synced.
fn write_durably(path: &Path, bytes: &[u8]) -> Result<()> {
    let part = with_suffix(path, ".part");
    let mut tmp = TempFile {
        path: part.clone(),
        keep: false,
    };
    let mut file = File::create(&part).with_context(|| format!("create {}", part.display()))?;
    file.write_all(bytes)?;
    file.sync_all()
        .with_context(|| format!("sync {}", part.display()))?;
    drop(file);
    std::fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))?;
    tmp.keep = true;
    sync_parent(path)
}

/// Bytes an unprivileged process may still write to the filesystem that
/// holds `path` (statvfs: f_bavail x f_frsize).
#[cfg(all(target_os = "linux", target_pointer_width = "64"))]
pub fn free_bytes(path: &Path) -> Result<u64> {
    use std::os::unix::ffi::OsStrExt;
    unsafe extern "C" {
        fn statvfs(path: *const std::ffi::c_char, buf: *mut u64) -> std::ffi::c_int;
    }
    let c_path = std::ffi::CString::new(path.as_os_str().as_bytes())?;
    // On 64-bit Linux (glibc and musl) struct statvfs is 112 bytes and
    // starts f_bsize, f_frsize, f_blocks, f_bfree, f_bavail, 8 bytes each.
    let mut buf = [0u64; 32];
    // SAFETY: c_path is NUL-terminated, and buf is 8-byte aligned, larger
    // than struct statvfs, and outlives the call.
    let rc = unsafe { statvfs(c_path.as_ptr(), buf.as_mut_ptr()) };
    if rc != 0 {
        return Err(std::io::Error::last_os_error())
            .with_context(|| format!("statvfs {}", path.display()));
    }
    Ok(buf[4].saturating_mul(buf[1]))
}

#[cfg(not(all(target_os = "linux", target_pointer_width = "64")))]
pub fn free_bytes(path: &Path) -> Result<u64> {
    anyhow::bail!("free-space checks for {} need 64-bit Linux", path.display())
}

/// Fails unless the filesystem of `dir` can take `need` more bytes and
/// still keep `reserve` bytes free. A caller checks this before counting a
/// block of ksets; `write_partial` and `FeatureMatrix::save_with_reserve`
/// check their own writes.
pub fn ensure_free_space(dir: &Path, need: u64, reserve: u64) -> Result<()> {
    ensure_free_with(&free_bytes, dir, need, reserve)
}

fn ensure_free_with(
    free: &dyn Fn(&Path) -> Result<u64>,
    dir: &Path,
    need: u64,
    reserve: u64,
) -> Result<()> {
    let free = free(dir)?;
    ensure!(
        free >= need.saturating_add(reserve),
        "{} has {free} bytes free: writing {need} more would cut into the {reserve}-byte reserve",
        dir.display()
    );
    Ok(())
}

/// What one partial write may use: the columns downstream stages are
/// reviewed for, and the disk it must leave free, measured by `free`.
struct Limits<'a> {
    max_columns: usize,
    reserve: u64,
    free: &'a dyn Fn(&Path) -> Result<u64>,
}

fn put_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

fn get_varint(data: &[u8], pos: &mut usize) -> Result<u64> {
    let mut v = 0u64;
    let mut shift = 0;
    loop {
        let Some(&byte) = data.get(*pos) else {
            anyhow::bail!("truncated varint")
        };
        *pos += 1;
        ensure!(shift < 63 || byte <= 1, "varint overflow");
        v |= ((byte & 0x7f) as u64) << shift;
        if byte < 0x80 {
            return Ok(v);
        }
        shift += 7;
    }
}

/// SHA-256 of an open file, from its first byte.
fn sha256_of(file: &File) -> Result<String> {
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    let mut buf = vec![0u8; 1 << 20];
    let mut at = 0u64;
    loop {
        let got = file.read_at(&mut buf, at)?;
        if got == 0 {
            break;
        }
        ctx.update(&buf[..got]);
        at += got as u64;
    }
    Ok(ctx
        .finish()
        .as_ref()
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect())
}

/// `kset::digest` of a kset's values, read one bucket at a time.
fn kset_digest(reader: &KsetReader) -> Result<u64> {
    let (mut bytes, mut values) = (Vec::new(), Vec::new());
    let mut h = Fnv::new();
    for b in 0..BUCKETS {
        values.clear();
        reader.read_bucket(b, &mut bytes, &mut values)?;
        values.iter().for_each(|&v| h.u64(v));
    }
    Ok(h.0)
}

/// Writes columns `first_col..first_col + kset_paths.len()` of a later
/// `build_from_partials` to `out`, and its receipt to `receipt_path(out)`.
///
/// `columns` are the genomes those columns must hold and the digests their
/// ksets must have; a kset with another name or other contents is refused
/// before anything is written. Any old partial and receipt at `out` are
/// then removed, so a failed rewrite leaves no stale block. Both new files
/// are synced, with their directory, before this returns, so the caller
/// may delete the block's ksets afterwards. Writing stops, and leaves no
/// file behind, rather than leave less than `reserve_bytes` free.
///
/// `first_col` must be a multiple of 64 so blocks join on word boundaries.
/// Only the block's ksets are open, and nothing is filtered yet: support is
/// only known once every block is merged.
pub fn write_partial(
    kset_paths: &[PathBuf],
    columns: &[Column],
    first_col: usize,
    out: &Path,
    threads_chunk: usize,
    reserve_bytes: u64,
) -> Result<PartialReceipt> {
    let limits = Limits {
        max_columns: REVIEWED_COLUMNS,
        reserve: reserve_bytes,
        free: &free_bytes,
    };
    write_partial_with(kset_paths, columns, first_col, out, threads_chunk, &limits)
}

fn write_partial_with(
    kset_paths: &[PathBuf],
    columns: &[Column],
    first_col: usize,
    out: &Path,
    threads_chunk: usize,
    limits: &Limits,
) -> Result<PartialReceipt> {
    ensure!(
        first_col.is_multiple_of(64),
        "a partial must start on a 64-column boundary, not column {first_col}"
    );
    ensure!(
        columns.len() == kset_paths.len(),
        "{} ksets for {} expected columns",
        kset_paths.len(),
        columns.len()
    );
    ensure_reviewed(first_col.saturating_add(columns.len()), limits.max_columns)?;
    if let Some(c) = columns.iter().find(|c| c.name.len() > MAX_NAME) {
        anyhow::bail!(
            "a genome name of {} bytes is longer than {MAX_NAME}",
            c.name.len()
        );
    }
    let (readers, k) = open_ksets(kset_paths)?;
    let width = readers.len();
    // Each column must be the expected genome, counted from the expected
    // contents; reading every value also checks each kset's data.
    let digests: Vec<u64> = readers.par_iter().map(kset_digest).collect::<Result<_>>()?;
    for (i, ((r, c), &d)) in readers.iter().zip(columns).zip(&digests).enumerate() {
        let col = first_col + i;
        ensure!(
            r.name == c.name,
            "column {col} expects genome {}, but its kset holds {}",
            c.name,
            r.name
        );
        let found = Fnv(d).hex();
        ensure!(
            found == c.kset_digest,
            "the kset of {} (column {col}) has digest {found}, not the expected {}: it was counted from other contents",
            c.name,
            c.kset_digest
        );
    }
    let plan = MergePlan {
        key: merge_key(k, width)?,
        k,
        n: width,
        keep: Keep::new(width, 0),
    };

    // A rewrite retires the old partial first: if it then fails, no stale
    // block is left to merge.
    let receipt_at = receipt_path(out);
    for old in [receipt_at.as_path(), out] {
        remove_if_present(old).with_context(|| format!("remove {}", old.display()))?;
    }
    sync_parent(out)?;

    let dir = parent_dir(out);
    let table_len = (BUCKETS + 1) * 8 + BUCKETS * 4;
    let table_at = PART_HEAD + columns.iter().map(|c| 12 + c.name.len()).sum::<usize>();
    ensure_free_with(
        limits.free,
        dir,
        (table_at + table_len) as u64,
        limits.reserve,
    )?;
    let part = with_suffix(out, ".part");
    let mut tmp = TempFile {
        path: part.clone(),
        keep: false,
    };
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(PART_MAGIC)?;
    w.write_all(&[k as u8, BUCKET_BITS as u8, 0, 0, 0, 0, 0, 0])?;
    w.write_all(&(first_col as u64).to_le_bytes())?;
    w.write_all(&(width as u64).to_le_bytes())?;
    w.write_all(&[0u8; 16])?; // n_kmers and full_rows, written last
    for (c, d) in columns.iter().zip(&digests) {
        w.write_all(&(c.name.len() as u32).to_le_bytes())?;
        w.write_all(c.name.as_bytes())?;
        w.write_all(&d.to_le_bytes())?;
    }
    w.write_all(&vec![0u8; table_len])?;

    let mut offsets = vec![0u64; BUCKETS + 1];
    let mut counts = vec![0u32; BUCKETS];
    let (mut written, mut kmers, mut full_rows) = (0u64, 0u64, 0u64);
    let chunk = threads_chunk.max(1);
    for start in (0..BUCKETS).step_by(chunk) {
        let end = (start + chunk).min(BUCKETS);
        let encoded: Vec<(Vec<u8>, u64, u64)> = (start..end)
            .into_par_iter()
            .map(|b| encode_partial_bucket(&plan, &readers, b))
            .collect::<Result<_>>()?;
        let len: usize = encoded.iter().map(|(data, ..)| data.len()).sum();
        ensure_free_with(limits.free, dir, len as u64, limits.reserve)?;
        for (b, (data, count, full)) in (start..end).zip(encoded) {
            offsets[b] = written;
            counts[b] = u32::try_from(count).context("too many k-mers in one bucket")?;
            w.write_all(&data)?;
            written += data.len() as u64;
            kmers += count;
            full_rows += full;
        }
    }
    offsets[BUCKETS] = written;
    let mut file = w.into_inner()?;
    file.seek(SeekFrom::Start(PART_HEAD as u64 - 16))?;
    file.write_all(&kmers.to_le_bytes())?;
    file.write_all(&full_rows.to_le_bytes())?;
    file.seek(SeekFrom::Start(table_at as u64))?;
    let mut table = Vec::with_capacity(table_len);
    offsets.iter().for_each(|o| table.extend(o.to_le_bytes()));
    counts.iter().for_each(|c| table.extend(c.to_le_bytes()));
    file.write_all(&table)?;
    file.sync_all()
        .with_context(|| format!("sync {}", part.display()))?;
    // File::create is write-only; hash the synced bytes through a read handle.
    let read = File::open(&part).with_context(|| format!("read {}", part.display()))?;
    let receipt = PartialReceipt {
        k,
        first_col,
        columns: columns.to_vec(),
        kmers,
        full_rows,
        bytes: file.metadata()?.len(),
        sha256: sha256_of(&read)?,
    };
    drop(read);
    drop(file);
    std::fs::rename(&part, out).with_context(|| format!("rename to {}", out.display()))?;
    tmp.keep = true;
    sync_parent(out)?;
    // Without its receipt the partial is refused, so a failure here leaves
    // the block to be rewritten.
    write_durably(&receipt_at, &serde_json::to_vec_pretty(&receipt)?)?;
    Ok(receipt)
}

/// One bucket of a partial: its encoded bytes, k-mers, and full rows.
fn encode_partial_bucket(
    plan: &MergePlan,
    readers: &[KsetReader],
    bucket: usize,
) -> Result<(Vec<u8>, u64, u64)> {
    let total: usize = readers.iter().map(|r| r.bucket_len(bucket)).sum();
    let base = (bucket as u64) << (2 * plan.k as u32 - BUCKET_BITS);
    let (mut data, mut count, mut full) = (Vec::new(), 0u64, 0u64);
    let mut prev = base;
    let mut bytes = Vec::new();
    merge_entries(
        plan,
        bucket,
        total,
        |g, out| readers[g].read_bucket(bucket, &mut bytes, out),
        |value, row, _| {
            // From the row itself, not the entry count, so that only genomes
            // whose bit is set can make a row full.
            let is_full = bits::count(row) as usize == plan.n;
            put_varint(&mut data, ((value - prev) << 1) | is_full as u64);
            if !is_full {
                row.iter().for_each(|w| data.extend(w.to_le_bytes()));
            }
            prev = value;
            count += 1;
            full += is_full as u64;
            Ok(())
        },
    )?;
    Ok((data, count, full))
}

/// Random access to one partial's buckets.
struct PartialReader {
    file: File,
    path: PathBuf,
    k: usize,
    first_col: usize,
    width: usize,
    columns: Vec<Column>,
    n_kmers: u64,
    full_rows: u64,
    bytes: u64,
    offsets: Vec<u64>,
    counts: Vec<u32>,
    data_start: u64,
    /// The row of a full k-mer: every genome of the block.
    full_row: Vec<u64>,
}

impl PartialReader {
    /// Opens a partial whose columns end by `max_columns`, and checks its
    /// header and bucket table; `read_bucket` checks the data.
    fn open(path: &Path, max_columns: usize) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        let mut head = [0u8; PART_HEAD];
        file.read_exact(&mut head)
            .with_context(|| format!("read the header of {}", path.display()))?;
        ensure!(
            &head[..8] == PART_MAGIC,
            "{} is not a partial matrix of this version",
            path.display()
        );
        let k = head[8] as usize;
        ensure!(
            head[9] as u32 == BUCKET_BITS && (MIN_K..=MAX_K).contains(&k),
            "{}: unsupported k {k} or bucket width {}",
            path.display(),
            head[9]
        );
        ensure!(
            head[10..16].iter().all(|&b| b == 0),
            "{}: reserved header bytes are set",
            path.display()
        );
        let word = |at: usize| u64::from_le_bytes(head[at..at + 8].try_into().unwrap());
        let (first_col, width, n_kmers, full_rows) = (word(16), word(24), word(32), word(40));
        // Every column costs at least 12 bytes, and no block ends past the
        // cap, which bounds the rows a corrupt width could make.
        ensure!(
            width >= 1
                && width <= file_len / 12
                && first_col.is_multiple_of(64)
                && first_col
                    .checked_add(width)
                    .is_some_and(|end| end <= max_columns as u64),
            "{}: bad column block {first_col} + {width}",
            path.display()
        );
        ensure!(
            full_rows <= n_kmers,
            "{}: more full rows than k-mers",
            path.display()
        );
        let mut r = BufReader::new(&file);
        let mut columns = Vec::with_capacity(width as usize);
        let mut pos = PART_HEAD as u64;
        for _ in 0..width {
            let mut len = [0u8; 4];
            r.read_exact(&mut len)
                .with_context(|| format!("{} is truncated", path.display()))?;
            let len = u32::from_le_bytes(len) as u64;
            ensure!(
                len <= MAX_NAME as u64,
                "{}: a genome name of {len} bytes",
                path.display()
            );
            pos += 4 + len + 8;
            ensure!(pos <= file_len, "{} is truncated", path.display());
            let mut name = vec![0u8; len as usize];
            r.read_exact(&mut name)?;
            let mut digest = [0u8; 8];
            r.read_exact(&mut digest)?;
            columns.push(Column {
                name: String::from_utf8(name)?,
                kset_digest: Fnv(u64::from_le_bytes(digest)).hex(),
            });
        }
        let mut table = vec![0u8; (BUCKETS + 1) * 8 + BUCKETS * 4];
        let data_start = pos + table.len() as u64;
        ensure!(data_start <= file_len, "{} is truncated", path.display());
        r.read_exact(&mut table)?;
        drop(r);
        let (off_bytes, count_bytes) = table.split_at((BUCKETS + 1) * 8);
        let offsets: Vec<u64> = off_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect();
        let counts: Vec<u32> = count_bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();
        ensure!(
            offsets[0] == 0
                && offsets[BUCKETS] == file_len - data_start
                && offsets
                    .windows(2)
                    .zip(&counts)
                    .all(|(w, &c)| w[1].checked_sub(w[0]).is_some_and(|len| c as u64 <= len))
                && counts.iter().map(|&c| c as u64).sum::<u64>() == n_kmers,
            "{}: bucket table does not match its data",
            path.display()
        );
        let width = width as usize;
        let words = words_for(width);
        let mut full_row = vec![u64::MAX; words];
        if !width.is_multiple_of(64) {
            full_row[words - 1] = (1u64 << (width % 64)) - 1;
        }
        Ok(PartialReader {
            file,
            path: path.to_path_buf(),
            k,
            first_col: first_col as usize,
            width,
            columns,
            n_kmers,
            full_rows,
            bytes: file_len,
            offsets,
            counts,
            data_start,
            full_row,
        })
    }

    /// Opens a partial and checks it against its receipt.
    fn open_with_receipt(path: &Path, max_columns: usize) -> Result<Self> {
        let reader = PartialReader::open(path, max_columns)?;
        let at = receipt_path(path);
        let text = std::fs::read(&at).with_context(|| {
            format!(
                "{} has no receipt {}: it is incomplete or was not written by write_partial",
                path.display(),
                at.display()
            )
        })?;
        let receipt: PartialReceipt =
            serde_json::from_slice(&text).with_context(|| format!("read {}", at.display()))?;
        ensure!(
            receipt.bytes == reader.bytes && receipt.sha256 == sha256_of(&reader.file)?,
            "{} does not match its receipt: it changed after it was written",
            path.display()
        );
        ensure!(
            receipt.k == reader.k
                && receipt.first_col == reader.first_col
                && receipt.columns == reader.columns
                && receipt.kmers == reader.n_kmers
                && receipt.full_rows == reader.full_rows,
            "{}: its header and receipt disagree",
            path.display()
        );
        Ok(reader)
    }

    /// Appends bucket `bucket`'s k-mers to `kmers` and, for each, a slot to
    /// `slots`: FULL, or the index of its row (words_for(width) words) in
    /// `rows`. Full rows are not expanded, so `rows` takes no more memory
    /// than the bytes it was read from. Returns how many k-mers were full.
    fn read_bucket(
        &self,
        bucket: usize,
        bytes: &mut Vec<u8>,
        kmers: &mut Vec<u64>,
        slots: &mut Vec<u32>,
        rows: &mut Vec<u64>,
    ) -> Result<u64> {
        let (start, end) = (self.offsets[bucket], self.offsets[bucket + 1]);
        bytes.clear();
        bytes.resize((end - start) as usize, 0);
        self.file.read_exact_at(bytes, self.data_start + start)?;
        let low = 2 * self.k as u32 - BUCKET_BITS;
        let words = self.full_row.len();
        let bad = || {
            format!(
                "{} bucket {bucket} does not match its count and range",
                self.path.display()
            )
        };
        let mut prev = (bucket as u64) << low;
        let (mut pos, mut full) = (0, 0u64);
        for i in 0..self.counts[bucket] {
            let tag = get_varint(bytes, &mut pos).with_context(bad)?;
            let gap = tag >> 1;
            ensure!(i == 0 || gap > 0, "{}", bad());
            prev = prev.checked_add(gap).with_context(bad)?;
            ensure!(prev >> low == bucket as u64, "{}", bad());
            kmers.push(prev);
            if tag & 1 == 1 {
                slots.push(FULL);
                full += 1;
                continue;
            }
            let row = bytes.get(pos..pos + words * 8).with_context(bad)?;
            pos += words * 8;
            let at = rows.len();
            slots.push(u32::try_from(at / words).with_context(bad)?);
            rows.extend(
                row.as_chunks::<8>()
                    .0
                    .iter()
                    .map(|c| u64::from_le_bytes(*c)),
            );
            let row = &rows[at..];
            let support = bits::count(row) as usize;
            let tail_clear = row
                .iter()
                .zip(&self.full_row)
                .all(|(w, mask)| w & !mask == 0);
            ensure!(
                tail_clear && support >= 1 && support < self.width,
                "{}",
                bad()
            );
        }
        ensure!(pos == bytes.len(), "{}", bad());
        Ok(full)
    }
}

/// Merges partials, given in column order, into the matrix `build` makes
/// from all of their ksets.
///
/// `columns` are the matrix's genomes, in order, with the kset digests the
/// caller expects (counted.json). Each partial must have its receipt, match
/// it byte for byte (SHA-256), and hold exactly the expected columns of its
/// block; partials of other genomes or other kset contents are refused.
/// Blocks must tile the columns from 0 without gaps or overlaps. One file
/// is open per block, not per genome. As with `build`, the caller fills in
/// the summary's `kset_digests`.
pub fn build_from_partials(
    partial_paths: &[PathBuf],
    columns: &[Column],
    min_support: u32,
    threads_chunk: usize,
) -> Result<FeatureMatrix> {
    ensure_reviewed(columns.len(), REVIEWED_COLUMNS)?;
    build_from_partials_with(
        partial_paths,
        columns,
        min_support,
        threads_chunk,
        REVIEWED_COLUMNS,
    )
}

fn build_from_partials_with(
    partial_paths: &[PathBuf],
    columns: &[Column],
    min_support: u32,
    threads_chunk: usize,
    max_columns: usize,
) -> Result<FeatureMatrix> {
    ensure!(!partial_paths.is_empty(), "no partial matrices to merge");
    ensure_unique_names(columns.iter().map(|c| c.name.as_str()))?;
    let blocks = partial_paths.len();
    let readers: Vec<PartialReader> = partial_paths
        .iter()
        .map(|p| {
            PartialReader::open_with_receipt(p, max_columns).map_err(|e| {
                with_open_file_hint(e, || {
                    format!(
                        "merging {blocks} partials keeps one open per block; raise the open-file limit (ulimit -n) above {blocks}, or write wider blocks"
                    )
                })
            })
        })
        .collect::<Result<_>>()?;
    let k = readers[0].k;
    ensure!(
        readers.iter().all(|r| r.k == k),
        "partial matrices mix different k values"
    );
    // Every partial starts on a 64-column boundary, so contiguous blocks
    // also leave only the last one off a word boundary.
    let mut n = 0;
    for r in &readers {
        ensure!(
            r.first_col == n,
            "{} holds columns from {}, but the next column is {n}: give partials in column order, without gaps",
            r.path.display(),
            r.first_col
        );
        n += r.width;
    }
    ensure!(
        n == columns.len(),
        "the partials hold {n} columns, but {} are expected",
        columns.len()
    );
    for r in &readers {
        let want = &columns[r.first_col..r.first_col + r.width];
        if let Some(i) = (0..r.width).find(|&i| r.columns[i] != want[i]) {
            let (got, want) = (&r.columns[i], &want[i]);
            anyhow::bail!(
                "{} is stale: column {} holds {} (kset {}), but {} (kset {}) is expected; rewrite this block",
                r.path.display(),
                r.first_col + i,
                got.name,
                got.kset_digest,
                want.name,
                want.kset_digest
            );
        }
    }
    merge_key(k, n)?;
    let keep = Keep::new(n, min_support);
    let genomes = columns.iter().map(|c| c.name.clone()).collect();
    let full_seen: Vec<AtomicU64> = readers.iter().map(|_| AtomicU64::new(0)).collect();
    let m = assemble(k, genomes, min_support, threads_chunk, |b| {
        merge_partial_bucket(&readers, &full_seen, keep, b)
    })?;
    for (r, seen) in readers.iter().zip(&full_seen) {
        let seen = seen.load(Ordering::Relaxed);
        ensure!(
            seen == r.full_rows,
            "{} holds {seen} full rows, but its header says {}",
            r.path.display(),
            r.full_rows
        );
    }
    Ok(m)
}

/// Merges one bucket of every partial by k-mer, joining their rows.
fn merge_partial_bucket(
    readers: &[PartialReader],
    full_seen: &[AtomicU64],
    keep: Keep,
    bucket: usize,
) -> Result<BucketOut> {
    let mut bytes = Vec::new();
    let mut blocks = Vec::with_capacity(readers.len());
    for (r, seen) in readers.iter().zip(full_seen) {
        let (mut kmers, mut slots, mut rows) = (Vec::new(), Vec::new(), Vec::new());
        let full = r.read_bucket(bucket, &mut bytes, &mut kmers, &mut slots, &mut rows)?;
        seen.fetch_add(full, Ordering::Relaxed);
        blocks.push((kmers, slots, rows, r));
    }
    let mut cursor = vec![0usize; blocks.len()];
    let mut row = vec![0u64; keep.words];
    let mut merge = BucketMerge::new(keep);
    loop {
        let next = blocks
            .iter()
            .zip(&cursor)
            .filter_map(|((kmers, ..), &c)| kmers.get(c))
            .min()
            .copied();
        let Some(value) = next else { break };
        row.fill(0);
        for ((kmers, slots, rows, r), c) in blocks.iter().zip(cursor.iter_mut()) {
            if kmers.get(*c) == Some(&value) {
                let words = r.full_row.len();
                let src = match slots[*c] {
                    FULL => &r.full_row[..],
                    s => &rows[s as usize * words..(s as usize + 1) * words],
                };
                let at = r.first_col / 64;
                row[at..at + words].copy_from_slice(src);
                *c += 1;
            }
        }
        merge.add(value, &row, bits::count(&row))?;
    }
    Ok(merge.finish())
}

fn write_u64s(path: &Path, values: &[u64]) -> Result<()> {
    let mut w =
        BufWriter::new(File::create(path).with_context(|| format!("create {}", path.display()))?);
    for v in values {
        w.write_all(&v.to_le_bytes())?;
    }
    w.flush()?;
    Ok(())
}

fn read_u64s(path: &Path) -> Result<Vec<u64>> {
    let mut bytes = Vec::new();
    BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?)
        .read_to_end(&mut bytes)?;
    Ok(bytes
        .as_chunks::<8>()
        .0
        .iter()
        .map(|c| u64::from_le_bytes(*c))
        .collect())
}

/// The files of a saved matrix.
const MATRIX_FILES: [&str; 5] = [
    "matrix.json",
    "patterns.bin",
    "support.bin",
    "kmer_offsets.bin",
    "kmers.bin",
];

impl FeatureMatrix {
    pub fn save(&self, dir: &Path) -> Result<()> {
        std::fs::create_dir_all(dir)?;
        self.write_files(dir, &serde_json::to_vec_pretty(&self.summary)?)
    }

    /// `save`, refused before anything is written or removed when the files
    /// would leave less than `reserve_bytes` free on `dir`'s filesystem.
    /// Files of an earlier save in `dir`, which are replaced, count as free.
    pub fn save_with_reserve(&self, dir: &Path, reserve_bytes: u64) -> Result<()> {
        self.save_checked(dir, reserve_bytes, &free_bytes)
    }

    fn save_checked(
        &self,
        dir: &Path,
        reserve: u64,
        free: &dyn Fn(&Path) -> Result<u64>,
    ) -> Result<()> {
        std::fs::create_dir_all(dir)?;
        let summary = serde_json::to_vec_pretty(&self.summary)?;
        let need = summary.len() as u64
            + 8 * (self.patterns.data.len() + self.kmer_offsets.len() + self.kmers.len()) as u64
            + 4 * self.support.len() as u64;
        let replaced: u64 = MATRIX_FILES
            .iter()
            .filter_map(|f| std::fs::metadata(dir.join(f)).ok())
            .map(|m| m.len())
            .sum();
        ensure_free_with(free, dir, need.saturating_sub(replaced), reserve)?;
        self.write_files(dir, &summary)
    }

    fn write_files(&self, dir: &Path, summary_json: &[u8]) -> Result<()> {
        // A summary left from an earlier save would mark half-written
        // binaries as a complete matrix, so it goes before they are touched.
        let summary = dir.join("matrix.json");
        remove_if_present(&summary).with_context(|| format!("remove {}", summary.display()))?;
        write_u64s(&dir.join("patterns.bin"), &self.patterns.data)?;
        let support: Vec<u8> = self.support.iter().flat_map(|s| s.to_le_bytes()).collect();
        std::fs::write(dir.join("support.bin"), support)?;
        write_u64s(&dir.join("kmer_offsets.bin"), &self.kmer_offsets)?;
        write_u64s(&dir.join("kmers.bin"), &self.kmers)?;
        // Summary last: its presence marks a complete matrix.
        std::fs::write(&summary, summary_json)?;
        Ok(())
    }

    /// Loads a saved matrix. One with more than `REVIEWED_COLUMNS` genomes
    /// is refused, so that only reviewed sizes reach downstream stages.
    pub fn load(dir: &Path) -> Result<Self> {
        let summary: MatrixSummary = serde_json::from_slice(
            &std::fs::read(dir.join("matrix.json"))
                .with_context(|| format!("no matrix in {}", dir.display()))?,
        )?;
        let n = summary.genomes.len();
        ensure_reviewed(n, REVIEWED_COLUMNS)
            .with_context(|| format!("load the matrix in {}", dir.display()))?;
        let words = words_for(n);
        let data = read_u64s(&dir.join("patterns.bin"))?;
        ensure!(
            data.len() == summary.n_patterns * words,
            "patterns.bin size does not match matrix.json"
        );
        let support: Vec<u32> = std::fs::read(dir.join("support.bin"))?
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();
        ensure!(
            support.len() == summary.n_patterns,
            "support.bin size does not match matrix.json"
        );
        let kmer_offsets = read_u64s(&dir.join("kmer_offsets.bin"))?;
        let kmers = read_u64s(&dir.join("kmers.bin"))?;
        ensure!(
            kmer_offsets.len() == summary.n_patterns + 1,
            "kmer_offsets.bin is incomplete"
        );
        ensure!(
            kmer_offsets[summary.n_patterns] == kmers.len() as u64,
            "kmers.bin size does not match kmer_offsets.bin"
        );
        Ok(FeatureMatrix {
            patterns: BitMatrix {
                rows: summary.n_patterns,
                cols: n,
                words,
                data,
            },
            summary,
            support,
            kmer_offsets,
            kmers,
        })
    }
}

/// A learning problem: features restricted to the labeled genomes, with
/// patterns that became identical (or constant) on that subset collapsed.
/// This uses every labeled genome, so it is only a candidate superset: each
/// fit applies `learn::training_candidates` to its own training samples.
pub struct Dataset {
    /// Features x samples.
    pub x: BitMatrix,
    /// Bit i set = sample i is resistant (the positive class).
    pub y: Bits,
    /// Column of each sample in the species matrix.
    pub columns: Vec<usize>,
    /// Feature f covers species patterns `feature_patterns[offsets[f]..offsets[f+1]]`.
    pub feature_offsets: Vec<u32>,
    pub feature_patterns: Vec<u32>,
}

impl Dataset {
    pub fn n_samples(&self) -> usize {
        self.x.cols
    }

    pub fn n_features(&self) -> usize {
        self.x.rows
    }

    pub fn patterns_of(&self, feature: usize) -> &[u32] {
        &self.feature_patterns
            [self.feature_offsets[feature] as usize..self.feature_offsets[feature + 1] as usize]
    }

    /// Projects `matrix` onto `columns` (with labels `resistant`), keeping
    /// features present in at least `min_support` and at most
    /// n - `min_support` of those samples.
    pub fn project(
        matrix: &FeatureMatrix,
        columns: Vec<usize>,
        resistant: &[bool],
        min_support: u32,
    ) -> Dataset {
        assert_eq!(columns.len(), resistant.len());
        let n = columns.len();
        let words = words_for(n);
        let max = (n as u32).saturating_sub(min_support);
        let rows = matrix.patterns.rows;

        // Gather in parallel chunks; collapse duplicates within each chunk
        // first so the sequential merge sees far fewer rows.
        let chunk = 1 << 15;
        let parts: Vec<(Vec<u64>, Vec<Vec<u32>>)> = (0..rows.div_ceil(chunk))
            .into_par_iter()
            .map(|c| {
                let mut local: FxHashMap<Box<[u64]>, usize> = FxHashMap::default();
                let mut out_rows = Vec::new();
                let mut members: Vec<Vec<u32>> = Vec::new();
                let mut row = vec![0u64; words];
                for p in c * chunk..((c + 1) * chunk).min(rows) {
                    let src = matrix.patterns.row(p);
                    row.iter_mut().for_each(|w| *w = 0);
                    for (i, &col) in columns.iter().enumerate() {
                        if bits::get(src, col) {
                            bits::set(&mut row, i);
                        }
                    }
                    let s = bits::count(&row);
                    if s < min_support || s > max {
                        continue;
                    }
                    let next = members.len();
                    let id = *local.entry(row.as_slice().into()).or_insert_with(|| {
                        out_rows.extend_from_slice(&row);
                        members.push(Vec::new());
                        next
                    });
                    members[id].push(p as u32);
                }
                (out_rows, members)
            })
            .collect();

        let mut index: FxHashMap<Box<[u64]>, usize> = FxHashMap::default();
        let mut x = BitMatrix::with_cols(n);
        let mut members: Vec<Vec<u32>> = Vec::new();
        for (part_rows, part_members) in parts {
            for (row, m) in part_rows.chunks_exact(words).zip(part_members) {
                let next = members.len();
                let id = *index.entry(row.into()).or_insert_with(|| {
                    x.push_row(row);
                    members.push(Vec::new());
                    next
                });
                members[id].extend(m);
            }
        }
        let mut feature_offsets = Vec::with_capacity(members.len() + 1);
        let mut feature_patterns = Vec::new();
        for m in &members {
            feature_offsets.push(feature_patterns.len() as u32);
            feature_patterns.extend_from_slice(m);
        }
        feature_offsets.push(feature_patterns.len() as u32);

        let y = Bits::from_indices(
            n,
            resistant
                .iter()
                .enumerate()
                .filter(|(_, r)| **r)
                .map(|(i, _)| i),
        );
        Dataset {
            x,
            y,
            columns,
            feature_offsets,
            feature_patterns,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna::KmerSpace;
    use crate::kset;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use std::collections::{BTreeMap, HashMap, HashSet};
    use std::sync::{Mutex, MutexGuard};
    use std::time::Instant;

    fn kset_of(dir: &Path, name: &str, seq: &str, space: &KmerSpace) -> PathBuf {
        let fna = dir.join(format!("{name}.fna"));
        std::fs::write(&fna, format!(">c\n{seq}\n")).unwrap();
        let (values, _) = kset::count_file(&fna, space, 1).unwrap();
        let path = dir.join(format!("{name}.kset"));
        kset::write(&path, name, space, &values).unwrap();
        path
    }

    #[test]
    fn merge_groups_kmers_by_presence_pattern() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        // A core block shared by all, plus an accessory block in g1 and g2.
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
            kset_of(dir.path(), "g2", &format!("{accessory}NNNN{core}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        assert_eq!(m.summary.genomes, vec!["g0", "g1", "g2"]);
        // Core k-mers (support 3 > n - 1) are dropped; accessory ones form a
        // single pattern {g1, g2}.
        assert_eq!(m.n_genomes(), 3);
        assert_eq!(m.summary.n_patterns, 1);
        assert_eq!(m.support, vec![2]);
        assert!(!m.patterns.get(0, 0) && m.patterns.get(0, 1) && m.patterns.get(0, 2));
        let mut distinct = std::collections::HashSet::new();
        space.for_each_canonical(accessory.as_bytes(), |x| {
            distinct.insert(x);
        });
        assert_eq!(m.kmers_of(0).len(), distinct.len());
        assert!(m.summary.kmers_core > 0);
        for &v in m.kmers_of(0) {
            let fwd = space.decode(space.unmix(v));
            let rc = space.decode(space.reverse_complement(space.unmix(v)));
            assert!(accessory.contains(&fwd) || accessory.contains(&rc));
        }

        m.save(&dir.path().join("mx")).unwrap();
        let loaded = FeatureMatrix::load(&dir.path().join("mx")).unwrap();
        assert_eq!(loaded.patterns.data, m.patterns.data);
        assert_eq!(loaded.kmers, m.kmers);

        // On columns {g0, g1} the accessory pattern is still variable.
        let ds = Dataset::project(&loaded, vec![0, 1], &[false, true], 1);
        assert_eq!(ds.n_features(), 1);
        assert_eq!(ds.patterns_of(0), &[0]);
        // On {g1, g2} it is constant and disappears.
        let ds = Dataset::project(&loaded, vec![1, 2], &[false, true], 1);
        assert_eq!(ds.n_features(), 0);
    }

    #[test]
    fn load_rejects_truncated_files() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        assert_eq!(m.summary.n_patterns, 1);
        // An interrupted save or copy leaves a short file behind.
        for (file, cut) in [("support.bin", 4), ("kmers.bin", 8)] {
            let mx = dir.path().join(file);
            m.save(&mx).unwrap();
            FeatureMatrix::load(&mx).unwrap();
            let path = mx.join(file);
            let len = std::fs::metadata(&path).unwrap().len();
            let f = std::fs::OpenOptions::new().write(true).open(&path).unwrap();
            f.set_len(len - cut).unwrap();
            assert!(
                FeatureMatrix::load(&mx).is_err(),
                "{file} short by {cut} bytes"
            );
        }
    }

    #[test]
    fn a_failed_resave_leaves_no_summary() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(7);
        let core = "ACGTTGCAAGGCTTACCGATA";
        let accessory = "TTTGGGCCCAAATGTGTCACA";
        let paths = vec![
            kset_of(dir.path(), "g0", core, &space),
            kset_of(dir.path(), "g1", &format!("{core}NNNN{accessory}"), &space),
        ];
        let m = build(&paths, 1, 8).unwrap();
        let mx = dir.path().join("mx");
        m.save(&mx).unwrap();
        // The second save rewrites three binaries, then cannot create kmers.bin.
        std::fs::remove_file(mx.join("kmers.bin")).unwrap();
        std::fs::create_dir(mx.join("kmers.bin")).unwrap();
        assert!(m.save(&mx).is_err());
        assert!(!mx.join("matrix.json").exists());
        let err = FeatureMatrix::load(&mx).err().unwrap();
        assert!(format!("{err:#}").contains("no matrix in"), "{err:#}");
    }

    /// Synthetic cohorts: a fixed pool of k-mer values (bucket edges, the
    /// extremes of the 2k-bit space, and seeded random values) shared out by
    /// membership groups. K-mers in one group get the same presence bitset,
    /// so duplicate patterns occur both within and across buckets.
    fn pool(k: usize) -> Vec<u64> {
        let bits = 2 * k as u32;
        let low = bits - BUCKET_BITS;
        let top = u64::MAX >> (64 - bits);
        let mut values = vec![0, 1, top, top - 1];
        for b in [1u64, 511, 512, BUCKETS as u64 - 1] {
            values.extend([(b << low) - 1, b << low]);
        }
        let mut rng = ChaCha8Rng::seed_from_u64(0x6d61_7472_6978 ^ k as u64);
        values.extend((0..48).map(|_| rng.random::<u64>() & top));
        values.sort_unstable();
        values.dedup();
        values
    }

    const GROUPS: usize = 9;

    fn member(group: usize, g: usize, n: usize) -> bool {
        match group {
            0 => true,
            1 => g == n - 1,
            2 => g == 0,
            3 => g != n - 1,
            4 => g >= n / 2,
            5 => g.is_multiple_of(3),
            6 => g + 2 >= n,
            _ => (g as u64 ^ 0x9e37_79b9_7f4a_7c15).wrapping_mul(group as u64 + 1) >> 63 == 1,
        }
    }

    fn cohort(k: usize, n: usize) -> Vec<Vec<u64>> {
        let pool = pool(k);
        (0..n)
            .map(|g| {
                pool.iter()
                    .enumerate()
                    .filter(|(j, _)| member(j % GROUPS, g, n))
                    .map(|(_, &v)| v)
                    .collect()
            })
            .collect()
    }

    fn write_cohort(dir: &Path, k: usize, sets: &[Vec<u64>]) -> Vec<PathBuf> {
        let space = KmerSpace::new(k);
        sets.iter()
            .enumerate()
            .map(|(g, values)| {
                let path = dir.join(format!("g{g:05}.kset"));
                kset::write(&path, &format!("g{g:05}"), &space, values).unwrap();
                path
            })
            .collect()
    }

    /// SHA-256 over every file a saved matrix consists of.
    fn digest_of(m: &FeatureMatrix) -> String {
        let dir = tempfile::tempdir().unwrap();
        m.save(dir.path()).unwrap();
        let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
        for file in [
            "matrix.json",
            "patterns.bin",
            "support.bin",
            "kmer_offsets.bin",
            "kmers.bin",
        ] {
            ctx.update(file.as_bytes());
            ctx.update(&std::fs::read(dir.path().join(file)).unwrap());
        }
        ctx.finish()
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    }

    /// The merged matrix re-derived the plain way: each k-mer's genome set
    /// in ascending k-mer order, patterns numbered by first appearance, and
    /// each pattern's k-mers ascending.
    #[derive(Debug, PartialEq)]
    struct Expected {
        patterns: Vec<u64>,
        support: Vec<u32>,
        kmer_offsets: Vec<u64>,
        kmers: Vec<u64>,
        distinct: u64,
        rare: u64,
        core: u64,
    }

    fn oracle(sets: &[Vec<u64>], min: u32) -> Expected {
        let n = sets.len();
        let max = (n as u32).saturating_sub(min);
        let mut by_kmer: BTreeMap<u64, Vec<usize>> = BTreeMap::new();
        for (g, set) in sets.iter().enumerate() {
            for &v in set {
                by_kmer.entry(v).or_default().push(g);
            }
        }
        let mut ids: HashMap<Vec<usize>, usize> = HashMap::new();
        let mut rows: Vec<Vec<usize>> = Vec::new();
        let mut members: Vec<Vec<u64>> = Vec::new();
        let (mut distinct, mut rare, mut core) = (0, 0, 0);
        for (v, genomes) in by_kmer {
            distinct += 1;
            let s = genomes.len() as u32;
            if s < min {
                rare += 1;
            } else if s > max {
                core += 1;
            } else {
                let id = *ids.entry(genomes.clone()).or_insert_with(|| {
                    rows.push(genomes);
                    members.push(Vec::new());
                    rows.len() - 1
                });
                members[id].push(v);
            }
        }
        let mut kmer_offsets = vec![0];
        for m in &members {
            kmer_offsets.push(kmer_offsets.last().unwrap() + m.len() as u64);
        }
        Expected {
            patterns: rows
                .iter()
                .flat_map(|r| Bits::from_indices(n, r.iter().copied()).words().to_vec())
                .collect(),
            support: rows.iter().map(|r| r.len() as u32).collect(),
            kmer_offsets,
            kmers: members.concat(),
            distinct,
            rare,
            core,
        }
    }

    fn observed(m: &FeatureMatrix) -> Expected {
        Expected {
            patterns: m.patterns.data.clone(),
            support: m.support.clone(),
            kmer_offsets: m.kmer_offsets.clone(),
            kmers: m.kmers.clone(),
            distinct: m.summary.kmers_distinct,
            rare: m.summary.kmers_rare,
            core: m.summary.kmers_core,
        }
    }

    /// Checks `m` against the oracle and the summary against the inputs.
    fn assert_merges(m: &FeatureMatrix, sets: &[Vec<u64>], min: u32, what: &str) {
        let want = oracle(sets, min);
        assert_eq!(observed(m), want, "{what}");
        let n = sets.len();
        assert_eq!(m.n_genomes(), n, "{what}");
        assert_eq!(m.patterns.words, words_for(n), "{what}");
        assert_eq!(m.summary.n_patterns, want.support.len(), "{what}");
        assert_eq!(m.summary.kmers_kept, want.kmers.len() as u64, "{what}");
    }

    /// Soft open-file limit of this process, None when unlimited.
    fn open_file_limit() -> Option<u64> {
        let limits = std::fs::read_to_string("/proc/self/limits").ok()?;
        let line = limits.lines().find(|l| l.starts_with("Max open files"))?;
        line.split_whitespace().nth(3)?.parse().ok()
    }

    /// Tests that hold thousands of ksets open take this in turn, so their
    /// open files never add up past the limit they check for.
    static MANY_FILES: Mutex<()> = Mutex::new(());

    fn hold_files(n: usize) -> MutexGuard<'static, ()> {
        let guard = MANY_FILES.lock().unwrap_or_else(|e| e.into_inner());
        if let Some(limit) = open_file_limit() {
            let need = n as u64 + 256;
            assert!(
                limit >= need,
                "this test keeps {n} ksets open; the open-file limit is {limit}, run with ulimit -n {need} or more"
            );
        }
        guard
    }

    /// (k, genomes, min_support) cohorts the packed merge supports, with
    /// the digests the merge produced before the wide key existed.
    const GOLDEN: [(usize, usize, u32, &str); 5] = [
        (
            5,
            1,
            0,
            "a3ccae1edde0cb16858420135ef9508c1705e532894c09cd18ce1e08fdc3fea8",
        ),
        (
            7,
            2,
            1,
            "7bfea74d30dee9eb1c2c20b6e4b9566f98e20884979fa40d414a2d40c65f6aef",
        ),
        (
            12,
            300,
            2,
            "d35148d75fb206e3ba7ada8bb787cdfb3c80f9645ec5a8c6bb1e1a987fc2bdf9",
        ),
        (
            31,
            4096,
            1,
            "aa2dae1eacfc16501dc0fc662db3335d659c48c19a28257e3c6f30a061321114",
        ),
        (
            32,
            1024,
            1,
            "c7cc503e78cd03413fa0652968f1b82389c114f904fdb0bcfe5636cd1ece8ccc",
        ),
    ];

    #[test]
    fn supported_cohorts_keep_their_digests_under_both_keys() {
        let _files = hold_files(4096);
        for (k, n, min, digest) in GOLDEN {
            assert_eq!(merge_key(k, n).unwrap(), MergeKey::Packed);
            let dir = tempfile::tempdir().unwrap();
            let sets = cohort(k, n);
            let paths = write_cohort(dir.path(), k, &sets);
            let packed = build(&paths, min, 64).unwrap();
            assert_eq!(digest_of(&packed), digest, "k={k} n={n}");
            let wide = build_with(&paths, min, 64, Some(MergeKey::Wide)).unwrap();
            assert_eq!(digest_of(&wide), digest, "wide k={k} n={n}");
            assert_merges(&packed, &sets, min, &format!("k={k} n={n}"));
        }
    }

    #[test]
    fn merge_key_is_packed_exactly_when_the_key_fits_64_bits() {
        use MergeKey::{Packed, Wide};
        for (k, n, key) in [
            (31, 4096, Packed),
            (31, 4097, Wide),
            (31, 8193, Wide),
            (31, 1 << 20, Wide),
            (32, 1024, Packed),
            (32, 1025, Wide),
            (32, 8193, Wide),
            (30, 16384, Packed),
            (30, 16385, Wide),
            (5, 1, Packed),
        ] {
            assert_eq!(merge_key(k, n).unwrap(), key, "k={k} n={n}");
        }
        for k in MIN_K..=MAX_K {
            for n in [
                1,
                2,
                3,
                1023,
                1024,
                1025,
                4096,
                4097,
                1 << 20,
                u32::MAX as usize,
            ] {
                let mut index_bits = 1;
                while (1u64 << index_bits) < n as u64 {
                    index_bits += 1;
                }
                let fits = 2 * k + index_bits <= 64 + BUCKET_BITS as usize;
                let key = merge_key(k, n).unwrap();
                assert_eq!(key == Packed, fits, "k={k} n={n}");
            }
            assert!(merge_key(k, u32::MAX as usize + 1).is_err(), "k={k}");
        }
    }

    /// Merges a synthetic cohort through files, checks it against the
    /// oracle, and checks the forced wide key gives the same bytes.
    fn merge_at(k: usize, n: usize) {
        let _files = hold_files(n);
        let dir = tempfile::tempdir().unwrap();
        let sets = cohort(k, n);
        let paths = write_cohort(dir.path(), k, &sets);
        let what = format!("k={k} n={n}");
        // Only tests may use this private entry point past REVIEWED_COLUMNS.
        let m = build_with(&paths, 1, 64, None).unwrap();
        assert_merges(&m, &sets, 1, &what);
        let names: Vec<String> = (0..n).map(|g| format!("g{g:05}")).collect();
        assert_eq!(m.summary.genomes, names, "{what}");
        // The last genome's private k-mers need the top genome-index bit.
        let last = (0..m.summary.n_patterns)
            .find(|&p| m.support[p] == 1 && m.patterns.get(p, n - 1))
            .unwrap_or_else(|| panic!("{what}: no pattern of the last genome alone"));
        assert!(!m.kmers_of(last).is_empty());
        let wide = build_with(&paths, 1, 64, Some(MergeKey::Wide)).unwrap();
        let digest = digest_of(&m);
        assert_eq!(digest_of(&wide), digest, "{what}");
        let parts = tempfile::tempdir().unwrap();
        let (blocks, receipts) = via_partials(&paths, 1024, 1, parts.path());
        assert_eq!(digest_of(&blocks), digest, "partials {what}");
        assert_eq!(receipts.len(), n.div_ceil(1024), "{what}");
        if merge_key(k, n).unwrap() == MergeKey::Wide {
            let err = build_with(&paths, 1, 64, Some(MergeKey::Packed))
                .err()
                .unwrap();
            assert!(format!("{err:#}").contains("do not fit the packed merge key"));
        }
    }

    #[test]
    fn k31_merges_at_the_packed_ceiling_of_4096() {
        merge_at(31, 4096);
    }

    #[test]
    fn k31_merges_past_the_ceiling_at_4097() {
        merge_at(31, 4097);
    }

    #[test]
    fn k31_merges_8193_genomes() {
        merge_at(31, 8193);
    }

    #[test]
    fn k32_merges_at_the_packed_ceiling_of_1024() {
        merge_at(32, 1024);
    }

    #[test]
    fn k32_merges_past_the_ceiling_at_1025() {
        merge_at(32, 1025);
    }

    #[test]
    fn k32_merges_4097_genomes() {
        merge_at(32, 4097);
    }

    #[test]
    fn rejects_no_ksets_and_mixed_k_at_any_size() {
        let err = build(&[], 1, 8).err().unwrap();
        assert!(format!("{err:#}").contains("no ksets to merge"));
        let _files = hold_files(4097);
        let dir = tempfile::tempdir().unwrap();
        let mut paths = write_cohort(dir.path(), 31, &cohort(31, 4096));
        let odd = dir.path().join("odd.kset");
        kset::write(&odd, "odd", &KmerSpace::new(32), &pool(32)).unwrap();
        for at in [0, 1, 4096] {
            let mut mixed = paths[..at.min(paths.len())].to_vec();
            mixed.push(odd.clone());
            mixed.extend_from_slice(&paths[at.min(paths.len())..]);
            for force in [None, Some(MergeKey::Wide)] {
                let err = build_with(&mixed, 1, 8, force).err().unwrap();
                assert!(
                    format!("{err:#}").contains("ksets mix different k values"),
                    "at={at}: {err:#}"
                );
            }
        }
        paths.truncate(1);
        paths.push(odd);
        assert!(build(&paths, 1, 8).is_err());
    }

    #[test]
    fn small_and_empty_cohorts_match_the_oracle_under_both_keys() {
        let dir = tempfile::tempdir().unwrap();
        for k in [5, 7, 31, 32] {
            for n in [1, 2, 3, 64, 65] {
                let full = cohort(k, n);
                let mut gapped = full.clone();
                gapped[n / 2].clear();
                for (label, sets) in [
                    ("full", full),
                    ("gapped", gapped),
                    ("empty", vec![vec![]; n]),
                ] {
                    let sub = dir.path().join(format!("{k}-{n}-{label}"));
                    std::fs::create_dir(&sub).unwrap();
                    let paths = write_cohort(&sub, k, &sets);
                    for min in [0, 1, 2] {
                        let what = format!("k={k} n={n} {label} min={min}");
                        let packed = build(&paths, min, 8).unwrap();
                        assert_merges(&packed, &sets, min, &what);
                        let digest = digest_of(&packed);
                        let wide = build_with(&paths, min, 8, Some(MergeKey::Wide)).unwrap();
                        assert_eq!(digest_of(&wide), digest, "{what}");
                        let parts = sub.join(format!("parts-{min}"));
                        std::fs::create_dir(&parts).unwrap();
                        let (blocks, _) = via_partials(&paths, 64, min, &parts);
                        assert_eq!(digest_of(&blocks), digest, "partials {what}");
                    }
                }
            }
        }
        // An empty matrix still saves and loads.
        let paths = write_cohort(&dir.path().join("7-3-empty"), 7, &[vec![], vec![], vec![]]);
        let m = build_with(&paths, 1, 8, Some(MergeKey::Wide)).unwrap();
        assert_eq!((m.summary.n_patterns, m.kmer_offsets.clone()), (0, vec![0]));
        m.save(&dir.path().join("mx")).unwrap();
        let loaded = FeatureMatrix::load(&dir.path().join("mx")).unwrap();
        assert_eq!(loaded.summary.n_patterns, 0);
    }

    #[test]
    fn duplicate_patterns_merge_across_buckets_under_both_keys() {
        let _files = hold_files(1025);
        for (k, n) in [(31, 9), (32, 1025)] {
            let dir = tempfile::tempdir().unwrap();
            let sets = cohort(k, n);
            let paths = write_cohort(dir.path(), k, &sets);
            let low = 2 * k as u32 - BUCKET_BITS;
            for force in [None, Some(MergeKey::Wide)] {
                let m = build_with(&paths, 1, 8, force).unwrap();
                let rows: HashSet<&[u64]> = (0..m.summary.n_patterns)
                    .map(|p| m.patterns.row(p))
                    .collect();
                assert_eq!(rows.len(), m.summary.n_patterns, "k={k} n={n}: rows repeat");
                let spread = (0..m.summary.n_patterns).any(|p| {
                    let buckets: HashSet<u64> = m.kmers_of(p).iter().map(|v| v >> low).collect();
                    buckets.len() > 1
                });
                assert!(spread, "k={k} n={n}: no pattern spans buckets");
                for p in 0..m.summary.n_patterns {
                    assert!(m.kmers_of(p).windows(2).all(|w| w[0] < w[1]));
                    assert_eq!(m.support[p], bits::count(m.patterns.row(p)));
                }
                assert!(m.summary.n_patterns < m.kmers.len());
            }
        }
    }

    #[test]
    fn digest_and_genome_order_do_not_depend_on_scheduling() {
        let _files = hold_files(1025);
        for (k, n) in [(31, 300), (32, 1025)] {
            let dir = tempfile::tempdir().unwrap();
            let sets = cohort(k, n);
            let paths = write_cohort(dir.path(), k, &sets);
            let first = digest_of(&build(&paths, 1, 64).unwrap());
            for chunk in [1, 7, BUCKETS] {
                for force in [None, Some(MergeKey::Wide)] {
                    let m = build_with(&paths, 1, chunk, force).unwrap();
                    assert_eq!(digest_of(&m), first, "k={k} n={n} chunk={chunk}");
                }
            }
            // Columns follow the given order, not file names or contents.
            let reversed: Vec<PathBuf> = paths.iter().rev().cloned().collect();
            let m = build(&reversed, 1, 64).unwrap();
            assert_eq!(m.summary.genomes[0], format!("g{:05}", n - 1));
            let flipped: Vec<Vec<u64>> = sets.iter().rev().cloned().collect();
            assert_merges(&m, &flipped, 1, &format!("reversed k={k} n={n}"));
        }
    }

    /// No column cap: how tests merge past `REVIEWED_COLUMNS`. Production
    /// code only reaches the public entry points, which keep the cap.
    const ANY_COLUMNS: usize = u32::MAX as usize;

    fn uncapped() -> Limits<'static> {
        Limits {
            max_columns: ANY_COLUMNS,
            reserve: 0,
            free: &free_bytes,
        }
    }

    /// Each kset's genome and `kset::digest`, as counted.json records them.
    fn columns_of(paths: &[PathBuf]) -> Vec<Column> {
        paths
            .iter()
            .map(|p| {
                let r = KsetReader::open(p).unwrap();
                Column {
                    name: r.name.clone(),
                    kset_digest: kset::digest(&r.read_all().unwrap()),
                }
            })
            .collect()
    }

    /// Writes `paths` to `dir` as partials of `widths` columns each (the
    /// last block takes what is left), `chunk` buckets at a time.
    fn write_blocks(
        paths: &[PathBuf],
        columns: &[Column],
        widths: &[usize],
        dir: &Path,
        chunk: usize,
    ) -> (Vec<PathBuf>, Vec<PartialReceipt>) {
        let (mut parts, mut receipts, mut start) = (Vec::new(), Vec::new(), 0);
        for (j, &width) in widths.iter().enumerate() {
            let end = (start + width).min(paths.len());
            let out = dir.join(format!("b{j:04}.kpart"));
            let r = write_partial_with(
                &paths[start..end],
                &columns[start..end],
                start,
                &out,
                chunk,
                &uncapped(),
            )
            .unwrap();
            receipts.push(r);
            parts.push(out);
            start = end;
        }
        assert_eq!(start, paths.len());
        (parts, receipts)
    }

    /// Merges `paths` through partials of `width` columns written to `dir`.
    fn via_partials(
        paths: &[PathBuf],
        width: usize,
        min: u32,
        dir: &Path,
    ) -> (FeatureMatrix, Vec<PartialReceipt>) {
        let columns = columns_of(paths);
        let widths = vec![width; paths.len().div_ceil(width)];
        let (parts, receipts) = write_blocks(paths, &columns, &widths, dir, 64);
        let m = build_from_partials_with(&parts, &columns, min, 64, ANY_COLUMNS).unwrap();
        (m, receipts)
    }

    fn sha256_hex(bytes: &[u8]) -> String {
        ring::digest::digest(&ring::digest::SHA256, bytes)
            .as_ref()
            .iter()
            .map(|b| format!("{b:02x}"))
            .collect()
    }

    /// Bytes of one partial column holding a 6-byte name (g00000).
    const COLUMN_BYTES: usize = 4 + 6 + 8;

    fn read_receipt(partial: &Path) -> PartialReceipt {
        serde_json::from_slice(&std::fs::read(receipt_path(partial)).unwrap()).unwrap()
    }

    fn write_receipt(partial: &Path, receipt: &PartialReceipt) {
        std::fs::write(receipt_path(partial), serde_json::to_vec(receipt).unwrap()).unwrap();
    }

    /// No partial, receipt, or temporary file at `out`.
    fn assert_nothing_at(out: &Path) {
        for path in [
            out.to_path_buf(),
            receipt_path(out),
            with_suffix(out, ".part"),
            with_suffix(&receipt_path(out), ".part"),
        ] {
            assert!(!path.exists(), "{} exists", path.display());
        }
    }

    #[test]
    fn partial_receipts_describe_their_blocks() {
        let dir = tempfile::tempdir().unwrap();
        let n = 130;
        let sets = cohort(31, n);
        let paths = write_cohort(dir.path(), 31, &sets);
        let columns = columns_of(&paths);
        let parts = dir.path().join("parts");
        std::fs::create_dir(&parts).unwrap();
        let (m, receipts) = via_partials(&paths, 64, 1, &parts);
        assert_eq!(digest_of(&m), digest_of(&build(&paths, 1, 64).unwrap()));
        for (j, r) in receipts.iter().enumerate() {
            let cols = j * 64..((j + 1) * 64).min(n);
            let block = &sets[cols.clone()];
            let mut seen: BTreeMap<u64, usize> = BTreeMap::new();
            block
                .iter()
                .flatten()
                .for_each(|&v| *seen.entry(v).or_default() += 1);
            let want: Vec<Column> = cols
                .clone()
                .map(|g| Column {
                    name: format!("g{g:05}"),
                    kset_digest: kset::digest(&sets[g]),
                })
                .collect();
            assert_eq!((r.k, r.first_col, &r.columns), (31, j * 64, &want));
            assert_eq!(r.kmers, seen.len() as u64, "block {j}");
            let full = seen.values().filter(|&&c| c == cols.len()).count();
            assert_eq!(r.full_rows, full as u64, "block {j}");
            let part = parts.join(format!("b{j:04}.kpart"));
            let bytes = std::fs::read(&part).unwrap();
            assert_eq!(
                (r.bytes, &r.sha256),
                (bytes.len() as u64, &sha256_hex(&bytes))
            );
            // The header records the totals, and the saved receipt is this one.
            let word = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
            assert_eq!((word(32), word(40)), (r.kmers, r.full_rows));
            assert_eq!(&read_receipt(&part), r);
            assert!(!with_suffix(&part, ".part").exists());
        }

        // Partial bytes do not depend on scheduling.
        let again = dir.path().join("again.kpart");
        for chunk in [1, 7, BUCKETS] {
            let r =
                write_partial(&paths[64..128], &columns[64..128], 64, &again, chunk, 0).unwrap();
            assert_eq!(r, receipts[1], "chunk={chunk}");
        }

        // Identical genomes make every k-mer full: no rows are stored, and
        // the block is far smaller than its ksets.
        let same = dir.path().join("same");
        std::fs::create_dir(&same).unwrap();
        let clones = write_cohort(&same, 31, &vec![pool(31); 64]);
        let r = write_partial(
            &clones,
            &columns_of(&clones),
            0,
            &same.join("p.kpart"),
            64,
            0,
        )
        .unwrap();
        assert_eq!((r.kmers, r.full_rows), (pool(31).len() as u64, r.kmers));
        let header = (PART_HEAD + 64 * COLUMN_BYTES + (BUCKETS + 1) * 8 + BUCKETS * 4) as u64;
        assert!(r.bytes <= header + r.kmers * 10, "{}", r.bytes);
        let kset_bytes: u64 = clones
            .iter()
            .map(|p| std::fs::metadata(p).unwrap().len())
            .sum();
        assert!(r.bytes * 32 < kset_bytes, "{} vs {kset_bytes}", r.bytes);
    }

    #[test]
    fn partials_must_tile_the_columns_in_order_with_one_k() {
        let dir = tempfile::tempdir().unwrap();
        let paths = write_cohort(dir.path(), 31, &cohort(31, 200));
        let columns = columns_of(&paths);
        let p: Vec<PathBuf> = (0..4)
            .map(|j| {
                let out = dir.path().join(format!("p{j}.kpart"));
                let cols = j * 64..((j + 1) * 64).min(200);
                write_partial(&paths[cols.clone()], &columns[cols], j * 64, &out, 8, 0).unwrap();
                out
            })
            .collect();
        let good = build_from_partials(&p, &columns, 1, 8).unwrap();
        assert_eq!(digest_of(&good), digest_of(&build(&paths, 1, 8).unwrap()));
        let order = |idx: &[usize]| -> Vec<PathBuf> { idx.iter().map(|&i| p[i].clone()).collect() };
        for idx in [
            &[1, 0, 2, 3][..],
            &[0, 2, 3],
            &[0, 0, 1, 2, 3],
            &[1, 2, 3],
            &[0, 1, 3, 2],
        ] {
            let err = build_from_partials(&order(idx), &columns, 1, 8)
                .err()
                .unwrap();
            assert!(
                format!("{err:#}").contains("give partials in column order"),
                "{idx:?}: {err:#}"
            );
        }
        // A missing last block, or more expected columns than blocks hold.
        let err = build_from_partials(&order(&[0, 1, 2]), &columns, 1, 8)
            .err()
            .unwrap();
        assert!(
            format!("{err:#}").contains("hold 192 columns, but 200"),
            "{err:#}"
        );
        let err = build_from_partials(&p, &columns[..199], 1, 8)
            .err()
            .unwrap();
        assert!(
            format!("{err:#}").contains("hold 200 columns, but 199"),
            "{err:#}"
        );
        let err = build_from_partials(&[], &columns, 1, 8).err().unwrap();
        assert!(format!("{err:#}").contains("no partial matrices"));
        let x = dir.path().join("x.kpart");
        let err = write_partial(&paths[..10], &columns[..10], 50, &x, 8, 0)
            .err()
            .unwrap();
        assert!(format!("{err:#}").contains("64-column boundary"));
        assert_nothing_at(&x);
        let err = write_partial(&paths[..10], &columns[..9], 0, &x, 8, 0)
            .err()
            .unwrap();
        assert!(format!("{err:#}").contains("10 ksets for 9 expected columns"));
        assert_nothing_at(&x);

        let other = dir.path().join("k32");
        std::fs::create_dir(&other).unwrap();
        let k32 = write_cohort(&other, 32, &cohort(32, 8));
        let odd = dir.path().join("k32.kpart");
        write_partial(&k32, &columns_of(&k32), 192, &odd, 8, 0).unwrap();
        let mut mixed = p[..3].to_vec();
        mixed.push(odd);
        let err = build_from_partials(&mixed, &columns, 1, 8).err().unwrap();
        assert!(format!("{err:#}").contains("mix different k"), "{err:#}");
    }

    /// Opens a partial and decodes every bucket, without its receipt.
    fn read_partial(path: &Path) -> Result<()> {
        let r = PartialReader::open(path, REVIEWED_COLUMNS)?;
        let mut bufs = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
        for b in 0..BUCKETS {
            r.read_bucket(b, &mut bufs.0, &mut bufs.1, &mut bufs.2, &mut bufs.3)?;
        }
        Ok(())
    }

    #[test]
    fn partial_files_reject_corruption() {
        let dir = tempfile::tempdir().unwrap();
        let paths = write_cohort(dir.path(), 31, &cohort(31, 64));
        let columns = columns_of(&paths);
        let good = dir.path().join("good.kpart");
        let receipt = write_partial(&paths, &columns, 0, &good, 8, 0).unwrap();
        let bytes = std::fs::read(&good).unwrap();
        read_partial(&good).unwrap();
        // 64 columns with 6-byte names, then the table.
        let counts_at = PART_HEAD + 64 * COLUMN_BYTES + (BUCKETS + 1) * 8;
        let data_at = counts_at + BUCKETS * 4;
        let counts: Vec<u32> = bytes[counts_at..data_at]
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();
        let busy = counts.iter().position(|&c| c > 0).unwrap();
        let n_kmers = u64::from_le_bytes(bytes[32..40].try_into().unwrap());
        let full_rows = u64::from_le_bytes(bytes[40..48].try_into().unwrap());
        assert!(full_rows > 0 && full_rows < n_kmers);
        let patched = |at: usize, value: &[u8]| {
            let mut out = bytes.clone();
            out[at..at + value.len()].copy_from_slice(value);
            out
        };
        let recount = |delta: i64| {
            let c = (counts[busy] as i64 + delta) as u32;
            let mut out = patched(counts_at + busy * 4, &c.to_le_bytes());
            let n = (n_kmers as i64 + delta) as u64;
            out[32..40].copy_from_slice(&n.to_le_bytes());
            out
        };
        let digest_at = PART_HEAD + 4 + 6;
        let mut cases = vec![
            ("magic", patched(0, b"GRMKSET1")),
            ("version 1", patched(0, b"GRMKPRT1")),
            ("k=4", patched(8, &[4])),
            ("k=33", patched(8, &[33])),
            ("bucket bits", patched(9, &[11])),
            ("reserved", patched(12, &[1])),
            ("first_col", patched(16, &1u64.to_le_bytes())),
            (
                "first_col past the cap",
                patched(16, &4096u64.to_le_bytes()),
            ),
            ("width 0", patched(24, &0u64.to_le_bytes())),
            ("width past the cap", patched(24, &4097u64.to_le_bytes())),
            ("width huge", patched(24, &u64::MAX.to_le_bytes())),
            ("n_kmers", patched(32, &(n_kmers + 1).to_le_bytes())),
            (
                "full rows > k-mers",
                patched(40, &(n_kmers + 1).to_le_bytes()),
            ),
            ("name length", patched(PART_HEAD, &u32::MAX.to_le_bytes())),
            ("name too long", patched(PART_HEAD, &4097u32.to_le_bytes())),
            ("one more k-mer", recount(1)),
            ("one fewer k-mer", recount(-1)),
            ("padded", [bytes.as_slice(), &[0]].concat()),
        ];
        for cut in [0, 39, PART_HEAD + 5, counts_at, data_at, bytes.len() - 1] {
            cases.push(("truncated", bytes[..cut].to_vec()));
        }
        let case = dir.path().join("case.kpart");
        for (what, data) in &cases {
            std::fs::write(&case, data).unwrap();
            assert!(
                read_partial(&case).is_err(),
                "{what} ({} bytes)",
                data.len()
            );
        }
        // These pass the reader; the merge refuses them by their receipt,
        // or, with a receipt forged to match, by checking the data.
        let merge_only = [
            ("full rows - 1", patched(40, &(full_rows - 1).to_le_bytes())),
            ("a digest", patched(digest_at, &[0xff])),
        ];
        for (what, data) in &merge_only {
            std::fs::write(&case, data).unwrap();
            read_partial(&case).unwrap_or_else(|e| panic!("{what}: {e:#}"));
        }
        for (what, data) in cases.iter().chain(&merge_only) {
            std::fs::write(&case, data).unwrap();
            write_receipt(&case, &receipt);
            let merge = || build_from_partials(std::slice::from_ref(&case), &columns, 1, 8);
            assert!(merge().is_err(), "{what}");
            let forged = PartialReceipt {
                bytes: data.len() as u64,
                sha256: sha256_hex(data),
                full_rows: if *what == "full rows - 1" {
                    full_rows - 1
                } else {
                    receipt.full_rows
                },
                ..receipt.clone()
            };
            write_receipt(&case, &forged);
            assert!(merge().is_err(), "forged receipt: {what}");
        }
        // The forged receipt alone is not what refuses them.
        std::fs::copy(&good, &case).unwrap();
        write_receipt(&case, &receipt);
        build_from_partials(std::slice::from_ref(&case), &columns, 1, 8).unwrap();
    }

    #[test]
    fn receipts_bind_each_partial_to_its_file_and_columns() {
        let dir = tempfile::tempdir().unwrap();
        let sets = cohort(31, 130);
        let paths = write_cohort(dir.path(), 31, &sets);
        let columns = columns_of(&paths);
        let parts = dir.path().join("parts");
        std::fs::create_dir(&parts).unwrap();
        let (p, receipts) = write_blocks(&paths, &columns, &[64, 64, 64], &parts, 64);
        let want = digest_of(&build(&paths, 1, 64).unwrap());
        let merge = |p: &[PathBuf]| build_from_partials(p, &columns, 1, 64);
        assert_eq!(digest_of(&merge(&p).unwrap()), want);
        let refused = |p: &[PathBuf], expect: &str| {
            let err = merge(p).err().unwrap();
            assert!(format!("{err:#}").contains(expect), "{expect}: {err:#}");
        };

        // No receipt: an interrupted write, or a file from elsewhere.
        let saved = std::fs::read(receipt_path(&p[1])).unwrap();
        std::fs::remove_file(receipt_path(&p[1])).unwrap();
        refused(&p, "has no receipt");
        std::fs::write(receipt_path(&p[1]), &saved).unwrap();
        // A receipt edited, or taken from another block.
        let mut edited = receipts[1].clone();
        edited.sha256 = sha256_hex(b"other");
        write_receipt(&p[1], &edited);
        refused(&p, "does not match its receipt");
        write_receipt(&p[1], &receipts[0]);
        refused(&p, "does not match its receipt");
        // A receipt that matches the file but not its header.
        let forged = PartialReceipt {
            kmers: receipts[1].kmers + 1,
            ..receipts[1].clone()
        };
        write_receipt(&p[1], &forged);
        refused(&p, "header and receipt disagree");
        std::fs::write(receipt_path(&p[1]), &saved).unwrap();
        std::fs::write(receipt_path(&p[1]), b"{").unwrap();
        refused(&p, "receipt.json");
        std::fs::write(receipt_path(&p[1]), &saved).unwrap();
        assert_eq!(digest_of(&merge(&p).unwrap()), want);

        // The file replaced, under its old receipt, by a valid partial of
        // the same genomes counted from other contents.
        let mut changed = sets.clone();
        changed[70].pop();
        let other = dir.path().join("other");
        std::fs::create_dir(&other).unwrap();
        let changed_paths = write_cohort(&other, 31, &changed);
        let replacement = other.join("b1.kpart");
        let changed_columns = columns_of(&changed_paths);
        write_partial(
            &changed_paths[64..128],
            &changed_columns[64..128],
            64,
            &replacement,
            64,
            0,
        )
        .unwrap();
        std::fs::copy(&replacement, &p[1]).unwrap();
        refused(&p, "does not match its receipt");
        // With its own receipt too, it holds the wrong kset for column 70.
        std::fs::copy(receipt_path(&replacement), receipt_path(&p[1])).unwrap();
        refused(&p, "column 70 holds g00070");

        // A leftover temporary file is never a partial, and a rewrite
        // replaces it.
        let leftover = with_suffix(&p[1], ".part");
        std::fs::write(&leftover, [b'x'; 100]).unwrap();
        refused(
            &[p[0].clone(), leftover.clone(), p[2].clone()],
            "not a partial",
        );
        write_partial(&paths[64..128], &columns[64..128], 64, &p[1], 64, 0).unwrap();
        assert!(!leftover.exists());
        assert_eq!(digest_of(&merge(&p).unwrap()), want);

        // Two outputs that share a stem get distinct temporary names.
        let (a, b) = (parts.join("x.a"), parts.join("x.b"));
        assert_ne!(with_suffix(&a, ".part"), with_suffix(&b, ".part"));
        assert_ne!(receipt_path(&a), receipt_path(&b));
        write_partial(&paths[..64], &columns[..64], 0, &a, 64, 0).unwrap();
        write_partial(&paths[64..128], &columns[64..128], 64, &b, 64, 0).unwrap();
        assert_eq!(read_receipt(&a), receipts[0]);
        assert_eq!(read_receipt(&b), receipts[1]);
        assert_eq!(digest_of(&merge(&[a, b, p[2].clone()]).unwrap()), want);
    }

    #[test]
    fn stale_partials_are_refused_after_catalog_or_assembly_changes() {
        let dir = tempfile::tempdir().unwrap();
        let sets = cohort(31, 200);
        let paths = write_cohort(dir.path(), 31, &sets);
        let columns = columns_of(&paths);
        let parts = dir.path().join("parts");
        std::fs::create_dir(&parts).unwrap();
        let (p, _) = write_blocks(&paths, &columns, &[64; 4], &parts, 64);
        let before = digest_of(&build(&paths, 1, 64).unwrap());
        assert_eq!(
            digest_of(&build_from_partials(&p, &columns, 1, 64).unwrap()),
            before
        );

        // Assembly change: g00070 is re-assembled and recounted.
        let mut changed = sets.clone();
        changed[70].retain(|v| v % 3 != 0);
        let space = KmerSpace::new(31);
        kset::write(&paths[70], "g00070", &space, &changed[70]).unwrap();
        let new_columns = columns_of(&paths);
        assert_ne!(new_columns[70], columns[70]);
        let err = build_from_partials(&p, &new_columns, 1, 64).err().unwrap();
        let msg = format!("{err:#}");
        assert!(msg.contains("b0001.kpart is stale: column 70"), "{msg}");
        // A block cannot be written from a kset other than the expected one.
        let err = write_partial(&paths[64..128], &columns[64..128], 64, &p[1], 64, 0)
            .err()
            .unwrap();
        assert!(
            format!("{err:#}").contains("the kset of g00070 (column 70) has digest"),
            "{err:#}"
        );
        // That refusal left the old block in place, still stale.
        assert!(build_from_partials(&p, &new_columns, 1, 64).is_err());
        // Rewriting the block makes the matrix `build` makes from the new ksets.
        write_partial(&paths[64..128], &new_columns[64..128], 64, &p[1], 64, 0).unwrap();
        let after = digest_of(&build(&paths, 1, 64).unwrap());
        assert_ne!(after, before);
        assert_eq!(
            digest_of(&build_from_partials(&p, &new_columns, 1, 64).unwrap()),
            after
        );

        // Catalog change: another genome takes column 130.
        let mut recatalogued = new_columns.clone();
        recatalogued[130].name = "h00130".into();
        let err = build_from_partials(&p, &recatalogued, 1, 64).err().unwrap();
        assert!(
            format!("{err:#}").contains("column 130 holds g00130"),
            "{err:#}"
        );
        let err = write_partial(&paths[128..192], &recatalogued[128..192], 128, &p[2], 64, 0)
            .err()
            .unwrap();
        assert!(
            format!("{err:#}").contains("expects genome h00130"),
            "{err:#}"
        );
        // Two columns may not name one genome.
        let mut twice = new_columns.clone();
        twice[199].name = twice[0].name.clone();
        let err = build_from_partials(&p, &twice, 1, 64).err().unwrap();
        assert!(format!("{err:#}").contains("in two columns"), "{err:#}");

        // A rewrite that fails after creating its file (here the disk fills
        // in the middle) leaves no partial, receipt, or temporary file, so
        // a resumed merge refuses the block until it is rewritten.
        let calls = std::cell::Cell::new(0);
        let filling = |_: &Path| -> Result<u64> {
            calls.set(calls.get() + 1);
            Ok(if calls.get() <= 3 { u64::MAX } else { 0 })
        };
        let full = Limits {
            max_columns: REVIEWED_COLUMNS,
            reserve: 1,
            free: &filling,
        };
        let err = write_partial_with(
            &paths[128..192],
            &new_columns[128..192],
            128,
            &p[2],
            64,
            &full,
        )
        .err()
        .unwrap();
        assert!(format!("{err:#}").contains("1-byte reserve"), "{err:#}");
        assert_eq!(calls.get(), 4);
        assert_nothing_at(&p[2]);
        let err = build_from_partials(&p, &new_columns, 1, 64).err().unwrap();
        assert!(format!("{err:#}").contains("b0002.kpart"), "{err:#}");
        write_partial(&paths[128..192], &new_columns[128..192], 128, &p[2], 64, 0).unwrap();
        assert_eq!(
            digest_of(&build_from_partials(&p, &new_columns, 1, 64).unwrap()),
            after
        );
    }

    #[test]
    fn free_space_is_measured_and_kept_for_partials_and_saves() {
        let dir = tempfile::tempdir().unwrap();
        let free = free_bytes(dir.path()).unwrap();
        assert!(free > 0);
        if let Ok(out) = std::process::Command::new("df")
            .args(["-B1", "--output=avail"])
            .arg(dir.path())
            .output()
            && out.status.success()
        {
            let text = String::from_utf8_lossy(&out.stdout);
            let df: u64 = text.lines().nth(1).unwrap().trim().parse().unwrap();
            assert!(free.abs_diff(df) < 1 << 30, "statvfs {free}, df {df}");
        }
        assert!(free_bytes(&dir.path().join("missing")).is_err());

        let paths = write_cohort(dir.path(), 31, &cohort(31, 64));
        let columns = columns_of(&paths);
        let out = dir.path().join("p.kpart");
        let err = write_partial(&paths, &columns, 0, &out, 8, u64::MAX)
            .err()
            .unwrap();
        assert!(format!("{err:#}").contains("reserve"), "{err:#}");
        assert_nothing_at(&out);
        assert!(ensure_free_space(dir.path(), 0, u64::MAX).is_err());
        ensure_free_space(dir.path(), 1 << 20, 0).unwrap();

        let m = build(&paths, 1, 8).unwrap();
        let mx = dir.path().join("mx");
        m.save_with_reserve(&mx, 0).unwrap();
        let saved = digest_of(&FeatureMatrix::load(&mx).unwrap());
        assert_eq!(saved, digest_of(&m));
        // Refused before the earlier matrix is touched.
        let err = m.save_with_reserve(&mx, u64::MAX).err().unwrap();
        assert!(format!("{err:#}").contains("reserve"), "{err:#}");
        assert_eq!(digest_of(&FeatureMatrix::load(&mx).unwrap()), saved);
        // The files a save replaces count as free; nothing else does.
        let size: u64 = MATRIX_FILES
            .iter()
            .map(|f| std::fs::metadata(mx.join(f)).unwrap().len())
            .sum();
        let exact = |left: u64| move |_: &Path| -> Result<u64> { Ok(left) };
        m.save_checked(&mx, 0, &exact(0)).unwrap();
        let fresh = dir.path().join("fresh");
        assert!(m.save_checked(&fresh, 0, &exact(size - 1)).is_err());
        assert!(!fresh.join("matrix.json").exists());
        m.save_checked(&fresh, 0, &exact(size)).unwrap();
        assert!(
            m.save_checked(&dir.path().join("f2"), 10, &exact(size + 9))
                .is_err()
        );
    }

    #[test]
    fn merges_refuse_more_columns_than_reviewed_unless_a_test_lifts_the_cap() {
        // Refused before any file is opened: these paths do not exist.
        let fake: Vec<PathBuf> = (0..=REVIEWED_COLUMNS)
            .map(|g| PathBuf::from(format!("/nonexistent/g{g:05}.kset")))
            .collect();
        let reviewed = "4096 matrix columns that Dataset::project";
        let err = build(&fake, 1, 8).err().unwrap();
        assert!(format!("{err:#}").contains(reviewed), "{err:#}");
        let err = build(&fake[..REVIEWED_COLUMNS], 1, 8).err().unwrap();
        assert!(format!("{err:#}").contains("/nonexistent/"), "{err:#}");
        let columns: Vec<Column> = (0..=REVIEWED_COLUMNS)
            .map(|g| Column {
                name: format!("g{g:05}"),
                kset_digest: "0000000000000000".into(),
            })
            .collect();
        let out = Path::new("/nonexistent/p.kpart");
        for (first_col, width) in [(4096, 1), (4032, 65), (0, 4097)] {
            let block = first_col..first_col + width;
            let err = write_partial(&fake[..width], &columns[block], first_col, out, 8, 0)
                .err()
                .unwrap();
            assert!(
                format!("{err:#}").contains(reviewed),
                "{first_col}+{width}: {err:#}"
            );
        }
        let err = build_from_partials(&[out.to_path_buf()], &columns, 1, 8)
            .err()
            .unwrap();
        assert!(format!("{err:#}").contains(reviewed), "{err:#}");

        // A matrix past the cap, made through the test-only path, does not
        // load for downstream stages.
        let n = REVIEWED_COLUMNS + 1;
        let wide = FeatureMatrix {
            summary: MatrixSummary {
                k: 31,
                genomes: columns.iter().map(|c| c.name.clone()).collect(),
                min_support: 1,
                n_patterns: 0,
                kmers_distinct: 0,
                kmers_kept: 0,
                kmers_rare: 0,
                kmers_core: 0,
                kset_digests: vec![],
            },
            patterns: BitMatrix {
                rows: 0,
                cols: n,
                words: words_for(n),
                data: vec![],
            },
            support: vec![],
            kmer_offsets: vec![0],
            kmers: vec![],
        };
        let dir = tempfile::tempdir().unwrap();
        wide.save(dir.path()).unwrap();
        let err = FeatureMatrix::load(dir.path()).err().unwrap();
        assert!(format!("{err:#}").contains(reviewed), "{err:#}");
    }

    #[test]
    fn a_kset_that_repeats_a_kmer_is_refused_by_every_path() {
        // The review's case: g0 lists v twice and g1 lacks it. Counted
        // twice, v would reach support 2 = n, which made it a full row in a
        // partial and so set g1's bit.
        let dir = tempfile::tempdir().unwrap();
        let (k, space) = (31, KmerSpace::new(31));
        let v = 0x1234_5678u64 << 20;
        let w = v + 1000;
        let g0 = dir.path().join("g0.kset");
        kset::write(&g0, "g0", &space, &[v, v + 1]).unwrap();
        let mut bytes = std::fs::read(&g0).unwrap();
        // Bucket 0 is the only one used, so its last gap (1) ends the file.
        assert_eq!(bytes.last(), Some(&1));
        *bytes.last_mut().unwrap() = 0;
        std::fs::write(&g0, &bytes).unwrap();
        let g1 = dir.path().join("g1.kset");
        kset::write(&g1, "g1", &space, &[w]).unwrap();
        let paths = vec![g0, g1];
        let columns = vec![
            Column {
                name: "g0".into(),
                kset_digest: kset::digest(&[v, v]),
            },
            Column {
                name: "g1".into(),
                kset_digest: kset::digest(&[w]),
            },
        ];
        let repeats = |r: Result<FeatureMatrix>, what: &str| {
            let err = r.err().unwrap_or_else(|| panic!("{what} merged"));
            assert!(
                format!("{err:#}").contains("repeats a k-mer"),
                "{what}: {err:#}"
            );
        };
        for min in [0, 1] {
            repeats(build(&paths, min, 8), "packed");
            repeats(build_with(&paths, min, 8, Some(MergeKey::Wide)), "wide");
        }
        let out = dir.path().join("p.kpart");
        let err = write_partial(&paths, &columns, 0, &out, 8, 0)
            .err()
            .unwrap();
        assert!(format!("{err:#}").contains("repeats a k-mer"), "{err:#}");
        assert_nothing_at(&out);

        // Behind the reader, the merge itself refuses a repeated entry
        // under both keys, so no path can count it twice.
        let sets = [vec![v, v], vec![w]];
        for key in [MergeKey::Packed, MergeKey::Wide] {
            let err = try_merge_in_memory(key, k, &sets, 0).err().unwrap();
            assert!(
                format!("{err:#}").contains("lists k-mer"),
                "{key:?}: {err:#}"
            );
        }

        // Without the repeat, every path agrees with the oracle: v is g0's
        // alone and g1's bit stays clear.
        let clean = dir.path().join("clean");
        std::fs::create_dir(&clean).unwrap();
        let sets = vec![vec![v], vec![w]];
        let paths = write_cohort(&clean, k, &sets);
        for min in [0, 1] {
            let m = build(&paths, min, 8).unwrap();
            assert_merges(&m, &sets, min, &format!("clean min={min}"));
            let row = (0..m.summary.n_patterns)
                .find(|&p| m.kmers_of(p).contains(&v))
                .unwrap();
            assert!(m.patterns.get(row, 0) && !m.patterns.get(row, 1));
            let digest = digest_of(&m);
            let wide = build_with(&paths, min, 8, Some(MergeKey::Wide)).unwrap();
            assert_eq!(digest_of(&wide), digest);
            let parts = clean.join(format!("parts-{min}"));
            std::fs::create_dir(&parts).unwrap();
            let (blocks, receipts) = via_partials(&paths, 64, min, &parts);
            assert_eq!(digest_of(&blocks), digest, "partials min={min}");
            assert_eq!((receipts[0].kmers, receipts[0].full_rows), (2, 0));
        }
    }

    /// A cohort whose buckets 0, 1, 777 and 1023 are crowded with values
    /// each genome keeps with probability 0.6, so most of their rows are
    /// neither full nor sparse, on top of the sparse `cohort`.
    fn dense_cohort(k: usize, n: usize) -> Vec<Vec<u64>> {
        let mut sets = cohort(k, n);
        for bucket in [0, 1, 777, BUCKETS - 1] {
            for (set, extra) in sets.iter_mut().zip(bucket_sets(k, n, 1500, bucket, 600)) {
                set.extend(extra);
            }
        }
        for set in &mut sets {
            set.sort_unstable();
            set.dedup();
        }
        sets
    }

    #[test]
    fn dense_buckets_and_uneven_blocks_match_build_under_any_schedule() {
        let _files = hold_files(421);
        let (k, n) = (31, 421);
        let dir = tempfile::tempdir().unwrap();
        let sets = dense_cohort(k, n);
        let paths = write_cohort(dir.path(), k, &sets);
        let columns = columns_of(&paths);
        let built = build(&paths, 1, 64).unwrap();
        assert_merges(&built, &sets, 1, "dense");
        let digest = digest_of(&built);
        let wide = build_with(&paths, 1, 64, Some(MergeKey::Wide)).unwrap();
        assert_eq!(digest_of(&wide), digest);
        for (i, widths) in [&[64, 192, 128, 37][..], &[64; 7], &[384, 64]]
            .into_iter()
            .enumerate()
        {
            let mut first: Option<Vec<PartialReceipt>> = None;
            for chunk in [1, 7, BUCKETS] {
                let sub = dir.path().join(format!("w{i}-c{chunk}"));
                std::fs::create_dir(&sub).unwrap();
                let (parts, receipts) = write_blocks(&paths, &columns, widths, &sub, chunk);
                // Dense blocks store real rows, not only full ones.
                assert!(receipts.iter().all(|r| r.full_rows * 2 < r.kmers));
                match &first {
                    None => first = Some(receipts),
                    Some(want) => assert_eq!(&receipts, want, "{widths:?} chunk={chunk}"),
                }
                for merge_chunk in [1, 7, BUCKETS] {
                    let m = build_from_partials(&parts, &columns, 1, merge_chunk).unwrap();
                    assert_eq!(
                        digest_of(&m),
                        digest,
                        "{widths:?} chunk={chunk} merge_chunk={merge_chunk}"
                    );
                }
            }
        }
    }

    /// One bucket's values for `n` genomes without files: a pool of values
    /// inside `bucket`, each genome keeping each value with probability
    /// `keep` / 1000 (seeded).
    fn bucket_sets(k: usize, n: usize, pool: usize, bucket: usize, keep: u32) -> Vec<Vec<u64>> {
        let low = 2 * k as u32 - BUCKET_BITS;
        let mut rng = ChaCha8Rng::seed_from_u64(0x6275_636b ^ (k * n) as u64);
        let mut values: Vec<u64> = (0..pool)
            .map(|_| ((bucket as u64) << low) | (rng.random::<u64>() & ((1 << low) - 1)))
            .collect();
        values.sort_unstable();
        values.dedup();
        (0..n)
            .map(|_| {
                values
                    .iter()
                    .copied()
                    .filter(|_| rng.random_range(0..1000) < keep)
                    .collect()
            })
            .collect()
    }

    fn try_merge_in_memory(
        key: MergeKey,
        k: usize,
        sets: &[Vec<u64>],
        bucket: usize,
    ) -> Result<BucketOut> {
        let n = sets.len();
        let plan = MergePlan {
            key,
            k,
            n,
            keep: Keep::new(n, 1),
        };
        let total = sets.iter().map(Vec::len).sum();
        let mut merge = BucketMerge::new(plan.keep);
        merge_entries(
            &plan,
            bucket,
            total,
            |g, out| {
                out.extend_from_slice(&sets[g]);
                Ok(())
            },
            |value, row, support| merge.add(value, row, support),
        )?;
        Ok(merge.finish())
    }

    fn merge_in_memory(key: MergeKey, k: usize, sets: &[Vec<u64>], bucket: usize) -> BucketOut {
        try_merge_in_memory(key, k, sets, bucket).unwrap()
    }

    /// A bucket merge in the matrix's own layout, for the oracle.
    fn bucket_observed(out: &BucketOut, words: usize) -> Expected {
        let rows = out.rows.len() / words.max(1);
        let mut pairs = out.pairs.clone();
        pairs.sort_unstable();
        let mut kmer_offsets = vec![0];
        for id in 0..rows as u32 {
            kmer_offsets.push(pairs.iter().filter(|p| p.0 <= id).count() as u64);
        }
        Expected {
            patterns: out.rows.clone(),
            support: out.rows.chunks_exact(words).map(bits::count).collect(),
            kmer_offsets,
            kmers: pairs.iter().map(|p| p.1).collect(),
            distinct: out.distinct,
            rare: out.rare,
            core: out.core,
        }
    }

    #[test]
    fn bucket_merge_past_both_ceilings_matches_the_oracle() {
        // In memory, so cohorts far past the packed ceilings need no files.
        for (k, n) in [(31, 4096), (32, 1024), (31, 70_000), (32, 70_000)] {
            let bucket = 777;
            let sets = bucket_sets(k, n, 24, bucket, 500);
            let wide = merge_in_memory(MergeKey::Wide, k, &sets, bucket);
            assert_eq!(
                bucket_observed(&wide, words_for(n)),
                oracle(&sets, 1),
                "k={k} n={n}"
            );
            if merge_key(k, n).unwrap() == MergeKey::Packed {
                let packed = merge_in_memory(MergeKey::Packed, k, &sets, bucket);
                assert_eq!(
                    (packed.rows, packed.pairs),
                    (wide.rows, wide.pairs),
                    "k={k} n={n}"
                );
            }
        }
    }

    /// The benchmark's clonal-looking cohort: 4096 pooled values, nine in
    /// ten core (each genome misses one with probability dropout /
    /// 100,000), the rest in about half the genomes.
    fn clonal_pool(k: usize) -> Vec<u64> {
        let top = u64::MAX >> (64 - 2 * k);
        let mut rng = ChaCha8Rng::seed_from_u64(k as u64 ^ 0xb10c);
        let mut pool: Vec<u64> = (0..4096).map(|_| rng.random::<u64>() & top).collect();
        pool.sort_unstable();
        pool.dedup();
        pool
    }

    fn clonal_genome(pool: &[u64], k: usize, g: usize, dropout: u32) -> Vec<u64> {
        let mut rng = ChaCha8Rng::seed_from_u64(((g as u64) << 8) | k as u64);
        pool.iter()
            .enumerate()
            .filter(|(j, _)| {
                if j % 10 == 0 {
                    rng.random()
                } else {
                    rng.random_range(0..100_000) >= dropout
                }
            })
            .map(|(_, &v)| v)
            .collect()
    }

    /// What one disk-sparse run measured.
    struct BlockRun {
        matrix: FeatureMatrix,
        all_ksets: u64,
        peak_temp: u64,
        partials: u64,
        full_rows: u64,
        write_secs: f64,
        merge_secs: f64,
    }

    /// The disk-sparse sequence a caller would run, in `dir`: count a block
    /// of `width` genomes, fold it into a partial, delete the block's ksets,
    /// then the next block; finally merge the partials.
    fn run_blocks(k: usize, n: usize, width: usize, dropout: u32, dir: &Path) -> BlockRun {
        let pool = clonal_pool(k);
        let space = KmerSpace::new(k);
        let columns: Vec<Column> = (0..n)
            .map(|g| Column {
                name: format!("g{g:05}"),
                kset_digest: kset::digest(&clonal_genome(&pool, k, g, dropout)),
            })
            .collect();
        let (mut all_ksets, mut partials, mut peak_temp, mut full_rows) = (0u64, 0u64, 0u64, 0u64);
        let mut parts = Vec::new();
        let mut write_secs = 0.0;
        for (j, start) in (0..n).step_by(width).enumerate() {
            let block = dir.join(format!("block{j}"));
            std::fs::create_dir(&block).unwrap();
            let cols = start..(start + width).min(n);
            let paths: Vec<PathBuf> = cols
                .clone()
                .map(|g| {
                    let path = block.join(format!("g{g:05}.kset"));
                    let values = clonal_genome(&pool, k, g, dropout);
                    kset::write(&path, &format!("g{g:05}"), &space, &values).unwrap();
                    path
                })
                .collect();
            let block_bytes: u64 = paths
                .iter()
                .map(|p| std::fs::metadata(p).unwrap().len())
                .sum();
            all_ksets += block_bytes;
            let out = dir.join(format!("b{j:04}.kpart"));
            let t = Instant::now();
            let receipt =
                write_partial_with(&paths, &columns[cols], start, &out, 64, &uncapped()).unwrap();
            write_secs += t.elapsed().as_secs_f64();
            // The block's ksets and the finished partials coexist at most.
            peak_temp = peak_temp.max(partials + block_bytes + receipt.bytes);
            partials += receipt.bytes;
            full_rows += receipt.full_rows;
            std::fs::remove_dir_all(&block).unwrap();
            parts.push(out);
        }
        let t = Instant::now();
        let matrix = build_from_partials_with(&parts, &columns, 1, 64, ANY_COLUMNS).unwrap();
        BlockRun {
            matrix,
            all_ksets,
            peak_temp,
            partials,
            full_rows,
            write_secs,
            merge_secs: t.elapsed().as_secs_f64(),
        }
    }

    #[test]
    fn the_disk_sparse_sequence_makes_the_matrix_build_makes() {
        let _files = hold_files(300);
        for (k, n, width, dropout) in [(31, 300, 64, 200), (32, 200, 128, 0)] {
            let what = format!("k={k} n={n} width={width}");
            let dir = tempfile::tempdir().unwrap();
            let run = run_blocks(k, n, width, dropout, dir.path());
            // Only partials and receipts are left: no ksets, no .part files.
            let mut left: Vec<String> = std::fs::read_dir(dir.path())
                .unwrap()
                .map(|e| e.unwrap().file_name().into_string().unwrap())
                .collect();
            left.sort();
            let blocks = n.div_ceil(width);
            let want: Vec<String> = (0..blocks)
                .flat_map(|j| {
                    [
                        format!("b{j:04}.kpart"),
                        format!("b{j:04}.kpart.receipt.json"),
                    ]
                })
                .collect();
            assert_eq!(left, want, "{what}");
            assert!(
                run.full_rows > 0 && run.partials * 4 < run.all_ksets,
                "{what}"
            );
            // The same cohort, every kset at once, through `build`.
            let all = tempfile::tempdir().unwrap();
            let pool = clonal_pool(k);
            let sets: Vec<Vec<u64>> = (0..n)
                .map(|g| clonal_genome(&pool, k, g, dropout))
                .collect();
            let paths = write_cohort(all.path(), k, &sets);
            let built = build(&paths, 1, 64).unwrap();
            assert_merges(&built, &sets, 1, &what);
            assert_eq!(digest_of(&run.matrix), digest_of(&built), "{what}");
        }
    }

    #[test]
    fn partials_merge_under_an_open_file_limit_that_stops_build() {
        let exe = std::env::current_exe().unwrap();
        let out = std::process::Command::new("sh")
            .arg("-c")
            .arg("ulimit -Sn 96 && exec \"$0\" --exact matrix::tests::open_file_limit_probe --ignored --nocapture --test-threads=1")
            .arg(&exe)
            .env("GRM_OPEN_FILE_PROBE", "1")
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            out.status.success() && stdout.contains("1 passed"),
            "{stdout}\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        // The limited child's matrix is the one `build` makes here.
        let digest = stdout.split("PARTIAL-DIGEST ").nth(1).unwrap();
        let dir = tempfile::tempdir().unwrap();
        let paths = write_cohort(dir.path(), 31, &cohort(31, 200));
        assert_eq!(&digest[..64], digest_of(&build(&paths, 1, 8).unwrap()));
    }

    #[test]
    #[ignore = "run by partials_merge_under_an_open_file_limit_that_stops_build, under ulimit -n 96"]
    fn open_file_limit_probe() {
        if std::env::var_os("GRM_OPEN_FILE_PROBE").is_none() {
            return;
        }
        assert_eq!(open_file_limit(), Some(96));
        let dir = tempfile::tempdir().unwrap();
        let paths = write_cohort(dir.path(), 31, &cohort(31, 200));
        let err = build(&paths, 1, 8).err().unwrap();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("raise the open-file limit (ulimit -n) above 200"),
            "{msg}"
        );
        // 64-genome blocks stay under the limit.
        let parts = dir.path().join("parts");
        std::fs::create_dir(&parts).unwrap();
        let (m, _) = via_partials(&paths, 64, 1, &parts);
        println!("PARTIAL-DIGEST {}", digest_of(&m));
    }

    /// Peak resident set of this process so far, in KiB (Linux).
    fn peak_rss_kib() -> u64 {
        let status = std::fs::read_to_string("/proc/self/status").unwrap();
        let line = status.lines().find(|l| l.starts_with("VmHWM:")).unwrap();
        line.split_whitespace().nth(1).unwrap().parse().unwrap()
    }

    /// (input, k, genomes, key) per benchmark case. `files` merges a tiny
    /// cohort end to end; `memory` merges one bucket held in memory.
    const BENCH_CASES: [(&str, usize, usize, MergeKey); 18] = [
        ("files", 31, 1024, MergeKey::Packed),
        ("files", 31, 1024, MergeKey::Wide),
        ("files", 31, 4096, MergeKey::Packed),
        ("files", 31, 4096, MergeKey::Wide),
        ("files", 31, 4097, MergeKey::Wide),
        ("files", 31, 8193, MergeKey::Wide),
        ("files", 32, 1024, MergeKey::Packed),
        ("files", 32, 1024, MergeKey::Wide),
        ("files", 32, 1025, MergeKey::Wide),
        ("files", 32, 4097, MergeKey::Wide),
        ("memory", 31, 1024, MergeKey::Packed),
        ("memory", 31, 1024, MergeKey::Wide),
        ("memory", 31, 4096, MergeKey::Packed),
        ("memory", 31, 4096, MergeKey::Wide),
        ("memory", 31, 8193, MergeKey::Wide),
        ("memory", 32, 1024, MergeKey::Packed),
        ("memory", 32, 1024, MergeKey::Wide),
        ("memory", 32, 8193, MergeKey::Wide),
    ];

    /// Address-space limit of every benchmark process, in KiB.
    const BENCH_VMEM_KIB: usize = 4 << 20;

    /// (k, genomes, block width, core dropout per 100,000) per out-of-core
    /// case.
    const BLOCK_CASES: [(usize, usize, usize, u32); 5] = [
        (31, 4097, 512, 200),
        (31, 8193, 1024, 200),
        (31, 8193, 1024, 0),
        (31, 8193, 256, 200),
        (32, 4097, 1024, 200),
    ];

    /// Runs `test` in a fresh process with GRM_BENCH_CASE = `case` under an
    /// open-file limit of `nofile` and a 4 GiB address space, and returns
    /// its result row.
    fn bench_child(test: &str, case: usize, nofile: usize) -> String {
        let out = std::process::Command::new("sh")
            .arg("-c")
            .arg(format!(
                "ulimit -Sn {nofile} && ulimit -Sv {BENCH_VMEM_KIB} && exec \"$0\" --exact matrix::tests::{test} --ignored --nocapture --test-threads=1"
            ))
            .arg(std::env::current_exe().unwrap())
            .env("GRM_BENCH_CASE", case.to_string())
            .output()
            .unwrap();
        let stdout = String::from_utf8_lossy(&out.stdout);
        assert!(
            out.status.success(),
            "{test} {case}: {stdout}\n{}",
            String::from_utf8_lossy(&out.stderr)
        );
        // The harness may print the test name on the same line.
        let row = stdout.split("bench ").nth(1).unwrap();
        row.lines().next().unwrap().to_string()
    }

    #[test]
    #[ignore = "resource benchmark: cargo test --release --lib matrix::tests::bench_merge_resources -- --ignored --nocapture"]
    fn bench_merge_resources() {
        // One fresh process per case, so each peak RSS is its own. The
        // in-memory kset merge gets the cohort plus 64 open files; the
        // out-of-core merge only one block plus 64.
        println!(
            "input   k      n  key     nofile  entries  patterns  merge_s  ns/entry  setup_peak_mib  peak_mib"
        );
        for (i, (_, _, n, _)) in BENCH_CASES.iter().enumerate() {
            println!("{}", bench_child("bench_case", i, n + 64));
        }
        println!(
            "blocks  k      n  width  dropout  nofile  all_ksets_mib  peak_temp_mib  partials_mib  full_rows  write_s  merge_s  peak_mib"
        );
        for (i, (_, _, width, _)) in BLOCK_CASES.iter().enumerate() {
            println!("{}", bench_child("bench_block_case", i, width + 64));
        }
    }

    #[test]
    #[ignore = "one case of bench_merge_resources, run in its own process"]
    fn bench_block_case() {
        let Some(case) = std::env::var_os("GRM_BENCH_CASE") else {
            return;
        };
        let (k, n, width, dropout) = BLOCK_CASES[case.to_str().unwrap().parse::<usize>().unwrap()];
        let nofile = open_file_limit().unwrap();
        let dir = tempfile::tempdir().unwrap();
        let run = run_blocks(k, n, width, dropout, dir.path());
        let mib = |b: u64| b as f64 / (1 << 20) as f64;
        println!(
            "bench blocks  {k}  {n:>5}  {width:>5}  {:>6.3}%  {nofile:>6}  {:>13.1}  {:>13.1}  {:>12.2}  {:>9}  {:>7.2}  {:>7.2}  {:>8.1}",
            dropout as f64 / 1000.0,
            mib(run.all_ksets),
            mib(run.peak_temp),
            mib(run.partials),
            run.full_rows,
            run.write_secs,
            run.merge_secs,
            peak_rss_kib() as f64 / 1024.0
        );
    }

    #[test]
    #[ignore = "one case of bench_merge_resources, run in its own process"]
    fn bench_case() {
        let Some(case) = std::env::var_os("GRM_BENCH_CASE") else {
            return;
        };
        let (input, k, n, key) = BENCH_CASES[case.to_str().unwrap().parse::<usize>().unwrap()];
        let nofile = open_file_limit().unwrap();
        let (entries, patterns, setup_peak, secs) = if input == "files" {
            // 1024 pooled values, each genome keeping about half.
            let dir = tempfile::tempdir().unwrap();
            let top = u64::MAX >> (64 - 2 * k);
            let mut rng = ChaCha8Rng::seed_from_u64(k as u64);
            let mut pool: Vec<u64> = (0..1024).map(|_| rng.random::<u64>() & top).collect();
            pool.sort_unstable();
            pool.dedup();
            let sets: Vec<Vec<u64>> = (0..n)
                .map(|_| pool.iter().copied().filter(|_| rng.random()).collect())
                .collect();
            let entries: usize = sets.iter().map(Vec::len).sum();
            let paths = write_cohort(dir.path(), k, &sets);
            drop(sets);
            let setup_peak = peak_rss_kib();
            let t = Instant::now();
            let m = build_with(&paths, 1, 64, Some(key)).unwrap();
            (entries, m.summary.n_patterns, setup_peak, t.elapsed())
        } else {
            // One bucket at an assumed MTB-like density scaled down 8x:
            // 4.4M k-mers / 1024 buckets / 8 = 540 values per genome, 90%
            // kept, so most are near-core.
            let sets = bucket_sets(k, n, 600, 300, 900);
            let entries: usize = sets.iter().map(Vec::len).sum();
            let setup_peak = peak_rss_kib();
            let t = Instant::now();
            let out = merge_in_memory(key, k, &sets, 300);
            (
                entries,
                out.rows.len() / words_for(n),
                setup_peak,
                t.elapsed(),
            )
        };
        let secs = secs.as_secs_f64();
        println!(
            "bench {input:<6}  {k}  {n:>5}  {key:<6?}  {nofile:>6}  {entries:>7}  {patterns:>8}  {secs:>7.3}  {:>8.1}  {:>14.1}  {:>8.1}",
            secs * 1e9 / entries as f64,
            setup_peak as f64 / 1024.0,
            peak_rss_kib() as f64 / 1024.0
        );
    }
}
