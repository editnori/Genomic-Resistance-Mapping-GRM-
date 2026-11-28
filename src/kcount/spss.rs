//! Checked 2-bit packed Spectrum-Preserving String Set (SPSS / greedy
//! simplitig) archive format (`spss-archive`, magic `GRMSPS1\0`).
//!
//! At bacterial density near 5 million distinct 31-mers, 5.15 bytes is the
//! floor for a code that uses only the sorted values. Exploiting `(k - 1)`-base
//! overlaps between adjacent canonical k-mers collapses a presence-only k-mer
//! set into deterministic greedy simplitigs packed at 2 bits per base, while
//! leaving `.kset` v1 (`GRMKSET1`), `.grk` (`GRMKGR1\0`), and `.ktab`
//! (`GRMKTB1\0`) untouched as random-access and merge forms.
//!
//! Because `.spss` encodes presence-only k-mer sets, archiving a counted
//! `.ktab` requires explicit opt-in (`drop_counts = true` / `--drop-counts`)
//! and discards k-mer counts. Every `.spss` file carries an 80-byte fixed
//! header storing `k`, `InputModality`, `n_kmers`, `n_simplitigs`,
//! `total_bases`, the FNV-1a64 sorted-key digest, and a 32-byte SHA-256
//! header/body digest verified on full reads and `inspect --verify`.

use crate::confine;
use crate::dna::KmerSpace;
use crate::fasta;
use crate::fnv::Fnv;
use crate::kcount::format::{
    AnyTableReader, InputModality, MAGIC_GRKSET, MAGIC_KTAB, TableFormat, check_output_extension,
    get_uleb, hex_bytes, put_uleb, write_gr_kset, write_ktab,
};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{self, BUCKET_BITS, BUCKETS};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use ring::digest::{Context as Sha256Context, SHA256};
use rustc_hash::{FxBuildHasher, FxHashSet};
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;
use std::time::Instant;

pub const MAGIC_SPSS1: &[u8; 8] = b"GRMSPS1\0";

/// Fixed header length for `GRMSPS1\0` (48-byte metadata prefix + 32-byte SHA-256 digest).
const SPSS_HEADER_LEN: usize = 80;
const SPSS_DIGEST_START: usize = 48;
const SPSS_DIGEST_END: usize = 80;

/// Cap on eager `Vec::try_reserve` hints before streaming table buckets.
const MAX_EAGER_CAP: usize = 1 << 20;

/// Summary returned by [`write_spss_from_mixed`], [`archive_from_table`],
/// [`archive_from_fasta`], and [`archive_auto`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ArchiveSummary {
    pub sample: String,
    pub k: usize,
    pub modality: InputModality,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub archive_bytes: u64,
    pub bytes_per_kmer: f64,
    pub counts_dropped: bool,
    pub fnv64: String,
    pub header_sha256: String,
    pub elapsed_ms: f64,
}

/// Summary returned by [`expand_to_kset_v1`], [`expand_to_gr`], and
/// [`expand_to_ktab`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ExpandSummary {
    pub format: String,
    pub sample: String,
    pub k: usize,
    pub modality: Option<InputModality>,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub output_bytes: u64,
    pub fnv64: String,
    pub header_sha256: String,
    pub elapsed_ms: f64,
}

/// Decoded header metadata for an `spss-archive` (`GRMSPS1\0`) file.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct SpssHeader {
    pub k: usize,
    pub modality: InputModality,
    pub name: String,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub fnv64_digest: u64,
    pub digest_sha256: [u8; 32],
    pub archive_bytes: u64,
}

/// Result of full `.spss` verification (SHA-256 + simplitig decode + FNV-1a64).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SpssVerifyResult {
    pub header: SpssHeader,
    pub keys_sha256: String,
    pub sketch_hex: Option<Vec<String>>,
}

struct ParsedSpssLayout {
    header: SpssHeader,
    lengths_off: u64,
    lengths_byte_len: usize,
    packed_off: u64,
    packed_len: usize,
}

/// Allocates a zero-filled byte buffer using fallible `try_reserve_exact` so
/// damaged files can never trigger an OOM abort or capacity-overflow panic.
fn alloc_zeroed_checked(len: usize, label: &str) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    buf.try_reserve_exact(len)
        .with_context(|| format!("{label}: allocation of {len} bytes failed"))?;
    buf.resize(len, 0);
    Ok(buf)
}

/// Streams a `.spss` file in fixed 64 KiB chunks through SHA-256, covering
/// every byte except the 32-byte digest slot `[48..80]`.
fn compute_spss_sha256(file: &File, file_len: u64) -> Result<[u8; 32]> {
    ensure!(
        file_len >= SPSS_HEADER_LEN as u64,
        "SPSS archive truncated before header digest slot"
    );
    let mut sha = Sha256Context::new(&SHA256);
    let mut buf = [0u8; 64 * 1024];

    let mut pos = 0u64;
    let prefix_end = SPSS_DIGEST_START as u64;
    while pos < prefix_end {
        let take = ((prefix_end - pos) as usize).min(buf.len());
        file.read_exact_at(&mut buf[..take], pos)?;
        sha.update(&buf[..take]);
        pos += take as u64;
    }

    pos = SPSS_DIGEST_END as u64;
    while pos < file_len {
        let take = ((file_len - pos) as usize).min(buf.len());
        file.read_exact_at(&mut buf[..take], pos)?;
        sha.update(&buf[..take]);
        pos += take as u64;
    }

    let digest = sha.finish();
    let mut out = [0u8; 32];
    out.copy_from_slice(digest.as_ref());
    Ok(out)
}

/// Contiguous LSB-first 2-bit DNA base packer (4 bases per byte).
struct PackedBases {
    bytes: Vec<u8>,
    len: usize,
}

impl PackedBases {
    fn with_capacity_bases(cap_bases: usize) -> Self {
        Self {
            bytes: Vec::with_capacity(cap_bases.div_ceil(4)),
            len: 0,
        }
    }

    #[inline]
    fn push(&mut self, code: u8) {
        debug_assert!(code < 4);
        let shift = (self.len & 3) * 2;
        if shift == 0 {
            self.bytes.push(code);
        } else {
            let last = self.bytes.last_mut().unwrap();
            *last |= code << shift;
        }
        self.len += 1;
    }
}

/// Greedy ProphAsm-style simplitig builder over a sorted set of distinct raw
/// canonical k-mers. Output is deterministic because seeds and extensions
/// follow canonical ascending and A/C/G/T order.
fn build_simplitigs_generic<K: KmerKey>(
    k: usize,
    raw_sorted: &[K],
) -> (Vec<usize>, Vec<u8>, usize) {
    let space = K::new_space(k);
    let mut unvisited = FxHashSet::with_capacity_and_hasher(raw_sorted.len(), FxBuildHasher);
    for &x in raw_sorted {
        unvisited.insert(x);
    }

    let mask: u128 = if 2 * k == 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let left_shift = 2 * (k as u32 - 1);

    let mut lengths = Vec::new();
    let mut packed = PackedBases::with_capacity_bases(raw_sorted.len() + k * 64);
    let mut right_buf: Vec<u8> = Vec::with_capacity(1024);
    let mut left_buf: Vec<u8> = Vec::with_capacity(1024);

    for &seed in raw_sorted {
        if !unvisited.remove(&seed) {
            continue;
        }
        right_buf.clear();
        left_buf.clear();

        // Extend right greedily in deterministic A/C/G/T order.
        let mut cur = seed;
        loop {
            let shifted = (cur.to_u128() << 2) & mask;
            let mut extended = false;
            for b in 0u8..4 {
                let cand_fwd = K::from_u128(shifted | (b as u128));
                let cand_rc = K::reverse_complement(&space, cand_fwd);
                let cand_canon = cand_fwd.min(cand_rc);
                if unvisited.remove(&cand_canon) {
                    right_buf.push(b);
                    cur = cand_fwd;
                    extended = true;
                    break;
                }
            }
            if !extended {
                break;
            }
        }

        // Extend left greedily in deterministic A/C/G/T order.
        let mut cur = seed;
        loop {
            let shifted = cur.to_u128() >> 2;
            let mut extended = false;
            for b in 0u8..4 {
                let cand_fwd = K::from_u128(shifted | ((b as u128) << left_shift));
                let cand_rc = K::reverse_complement(&space, cand_fwd);
                let cand_canon = cand_fwd.min(cand_rc);
                if unvisited.remove(&cand_canon) {
                    left_buf.push(b);
                    cur = cand_fwd;
                    extended = true;
                    break;
                }
            }
            if !extended {
                break;
            }
        }

        for &b in left_buf.iter().rev() {
            packed.push(b);
        }
        let seed_bits = seed.to_u128();
        for i in (0..k).rev() {
            packed.push(((seed_bits >> (2 * i)) & 3) as u8);
        }
        for &b in &right_buf {
            packed.push(b);
        }
        lengths.push(left_buf.len() + k + right_buf.len());
    }

    debug_assert!(unvisited.is_empty());
    (lengths, packed.bytes, packed.len)
}

struct SpssPayload<'a> {
    name: &'a str,
    k: usize,
    modality: InputModality,
    n_kmers: u64,
    fnv64_digest: u64,
    lengths: &'a [usize],
    packed_bases: &'a [u8],
    total_bases: u64,
}

fn write_spss_file(dst: &Path, p: &SpssPayload<'_>) -> Result<(u64, [u8; 32])> {
    ensure!((MIN_K..=MAX_K).contains(&p.k), "unsupported k {}", p.k);
    let name_len = u32::try_from(p.name.len()).context("sample name exceeds u32 length")?;
    let mut len_bytes = Vec::with_capacity(p.lengths.len().saturating_mul(2));
    for &l in p.lengths {
        ensure!(l >= p.k, "simplitig length {l} shorter than k={}", p.k);
        put_uleb(&mut len_bytes, (l - p.k) as u64);
    }
    let lengths_byte_len = len_bytes.len() as u64;

    // Fixed header prefix [0..48]; bytes [48..80] hold SHA-256 of
    // [0..48] || name || lengths_byte_len || len_bytes || packed_bases.
    let mut prefix = [0u8; SPSS_DIGEST_START];
    prefix[0..8].copy_from_slice(MAGIC_SPSS1);
    prefix[8] = p.k as u8;
    prefix[9] = p.modality.as_byte();
    // [10..12] reserved = 0
    prefix[12..16].copy_from_slice(&name_len.to_le_bytes());
    prefix[16..24].copy_from_slice(&p.n_kmers.to_le_bytes());
    prefix[24..32].copy_from_slice(&(p.lengths.len() as u64).to_le_bytes());
    prefix[32..40].copy_from_slice(&p.total_bases.to_le_bytes());
    prefix[40..48].copy_from_slice(&p.fnv64_digest.to_le_bytes());

    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&prefix);
    sha.update(p.name.as_bytes());
    sha.update(&lengths_byte_len.to_le_bytes());
    sha.update(&len_bytes);
    sha.update(p.packed_bases);
    let digest = sha.finish();
    let mut digest_arr = [0u8; 32];
    digest_arr.copy_from_slice(digest.as_ref());

    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }

    let part = confine::suffixed(dst, ".part");
    confine::regular_or_absent(dst)?;
    let mut w = BufWriter::new(confine::create_part(&part)?);
    let written = (|| -> Result<()> {
        w.write_all(&prefix)?;
        w.write_all(&digest_arr)?;
        w.write_all(p.name.as_bytes())?;
        w.write_all(&lengths_byte_len.to_le_bytes())?;
        w.write_all(&len_bytes)?;
        w.write_all(p.packed_bases)?;
        w.flush()?;
        confine::replace_with(&part, dst)
    })();
    drop(w);
    if written.is_err() && confine::regular_or_absent(&part).unwrap_or(false) {
        let _ = fs::remove_file(&part);
    }
    written?;
    Ok((fs::metadata(dst)?.len(), digest_arr))
}

#[allow(clippy::too_many_arguments)]
fn finish_archive_from_raw<K: KmerKey>(
    dst: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    raw_sorted: &[K],
    fnv: Fnv,
    counts_dropped: bool,
    t0: Instant,
) -> Result<ArchiveSummary> {
    let (lengths, packed_bases, total_bases) = build_simplitigs_generic(k, raw_sorted);
    let n_kmers = raw_sorted.len() as u64;
    let n_simplitigs = lengths.len() as u64;
    let total_bases_u64 = total_bases as u64;
    let expected_bases = (k as u64 - 1)
        .checked_mul(n_simplitigs)
        .and_then(|o| n_kmers.checked_add(o))
        .context("SPSS total_bases overflow")?;
    ensure!(
        total_bases_u64 == expected_bases,
        "SPSS invariant violated: total_bases={total_bases_u64}, n_kmers={n_kmers}, n_simplitigs={n_simplitigs}"
    );

    let (archive_bytes, digest_sha256) = write_spss_file(
        dst,
        &SpssPayload {
            name,
            k,
            modality,
            n_kmers,
            fnv64_digest: fnv.0,
            lengths: &lengths,
            packed_bases: &packed_bases,
            total_bases: total_bases_u64,
        },
    )?;
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
    let bytes_per_kmer = if n_kmers == 0 {
        0.0
    } else {
        archive_bytes as f64 / n_kmers as f64
    };
    Ok(ArchiveSummary {
        sample: name.to_string(),
        k,
        modality,
        n_kmers,
        n_simplitigs,
        total_bases: total_bases_u64,
        archive_bytes,
        bytes_per_kmer,
        counts_dropped,
        fnv64: fnv.hex(),
        header_sha256: hex_bytes(&digest_sha256),
        elapsed_ms,
    })
}

/// Archives sorted distinct mixed keys directly into a verified `GRMSPS1\0`
/// SPSS archive.
pub fn write_spss_from_mixed<K: KmerKey>(
    dst: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    mixed_sorted: &[K],
) -> Result<ArchiveSummary> {
    let t0 = Instant::now();
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    ensure!(
        mixed_sorted.windows(2).all(|w| w[0] < w[1]),
        "input mixed keys must be strictly sorted and unique"
    );
    let space = K::new_space(k);
    let mut fnv = Fnv::new();
    let mut raw = Vec::new();
    raw.try_reserve(mixed_sorted.len())
        .context("SPSS raw key buffer allocation failed")?;
    for &y in mixed_sorted {
        y.with_le_bytes(|b| fnv.bytes(b));
        raw.push(K::unmix(&space, y));
    }
    let mut scratch = Vec::new();
    K::radix_sort(&mut raw, &mut scratch);
    drop(scratch);
    debug_assert!(raw.windows(2).all(|w| w[0] < w[1]));

    finish_archive_from_raw(dst, name, k, modality, &raw, fnv, false, t0)
}

/// Archives a `.kset` v1 (`GRMKSET1`), `.grk` (`GRMKGR1\0`), or `.ktab`
/// (`GRMKTB1\0`) table into a verified `GRMSPS1\0` SPSS archive.
///
/// Refuses counted `.ktab` inputs (non-presence-only or any count > 1) unless
/// `drop_counts` is `true`, because `.spss` stores presence-only k-mer sets
/// and discards counts.
pub fn archive_from_table(
    src: &Path,
    dst: &Path,
    name: &str,
    modality_override: Option<InputModality>,
    drop_counts: bool,
) -> Result<ArchiveSummary> {
    let t0 = Instant::now();
    let reader = AnyTableReader::open(src)?;
    reader.verify_integrity()?;
    let k = reader.k();
    let out_modality = modality_override.unwrap_or_else(|| reader.modality().unwrap_or_default());

    let mut counts_dropped = false;
    if let AnyTableReader::Tab(ref tab) = reader
        && (!tab.presence_only || tab.max_count_observed > 1)
    {
        ensure!(
            drop_counts,
            "{}: refusing to archive counted .ktab (max_count_observed={}) without --drop-counts; SPSS archives store presence-only k-mer sets and discard counts",
            src.display(),
            tab.max_count_observed
        );
        counts_dropped = tab.max_count_observed > 1;
    }

    let cap = usize::try_from(reader.n_kmers())
        .context("table n_kmers exceeds addressable memory")?
        .min(MAX_EAGER_CAP);

    if k <= MAX_K_U64 {
        let mut keys: Vec<u64> = Vec::new();
        keys.try_reserve(cap)
            .context("SPSS table key buffer allocation failed")?;
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        let mut counts = Vec::new();
        for b in 0..BUCKETS {
            counts.clear();
            reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if counts.iter().any(|&c| c > 1) {
                ensure!(
                    drop_counts,
                    "{}: refusing to archive counted table without --drop-counts",
                    src.display()
                );
                counts_dropped = true;
            }
        }
        let mut summary = write_spss_from_mixed(dst, name, k, out_modality, &keys)?;
        summary.counts_dropped = counts_dropped;
        summary.elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
        Ok(summary)
    } else {
        let mut keys: Vec<u128> = Vec::new();
        keys.try_reserve(cap)
            .context("SPSS table key buffer allocation failed")?;
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        let mut counts = Vec::new();
        for b in 0..BUCKETS {
            counts.clear();
            reader.read_bucket::<u128>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if counts.iter().any(|&c| c > 1) {
                ensure!(
                    drop_counts,
                    "{}: refusing to archive counted table without --drop-counts",
                    src.display()
                );
                counts_dropped = true;
            }
        }
        let mut summary = write_spss_from_mixed(dst, name, k, out_modality, &keys)?;
        summary.counts_dropped = counts_dropped;
        summary.elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
        Ok(summary)
    }
}

/// Archives one FASTA/FASTQ file (plain or `.gz`) directly into a verified
/// `GRMSPS1\0` SPSS archive at the requested `k`.
pub fn archive_from_fasta(
    src: &Path,
    dst: &Path,
    k: usize,
    name: &str,
    modality: InputModality,
) -> Result<ArchiveSummary> {
    let t0 = Instant::now();
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    if k <= MAX_K_U64 {
        let space = KmerSpace64::new(k);
        let mut raw = Vec::new();
        fasta::for_each_record(src, |_, seq| {
            space.for_each_canonical(seq, |x| raw.push(x));
        })?;
        let mut scratch = Vec::new();
        u64::radix_sort(&mut raw, &mut scratch);
        raw.dedup();
        let mut mixed: Vec<u64> = raw.iter().map(|&x| space.mix(x)).collect();
        u64::radix_sort(&mut mixed, &mut scratch);
        let mut fnv = Fnv::new();
        for &y in &mixed {
            y.with_le_bytes(|b| fnv.bytes(b));
        }
        drop(mixed);
        finish_archive_from_raw(dst, name, k, modality, &raw, fnv, false, t0)
    } else {
        let space = KmerSpace128::new(k);
        let mut raw = Vec::new();
        fasta::for_each_record(src, |_, seq| {
            space.for_each_canonical(seq, |x| raw.push(x));
        })?;
        let mut scratch = Vec::new();
        u128::radix_sort(&mut raw, &mut scratch);
        raw.dedup();
        let mut mixed: Vec<u128> = raw.iter().map(|&x| space.mix(x)).collect();
        u128::radix_sort(&mut mixed, &mut scratch);
        let mut fnv = Fnv::new();
        for &y in &mixed {
            y.with_le_bytes(|b| fnv.bytes(b));
        }
        drop(mixed);
        finish_archive_from_raw(dst, name, k, modality, &raw, fnv, false, t0)
    }
}

/// Detects whether `src` is a k-mer table (`GRMKSET1`, `GRMKGR1\0`,
/// `GRMKTB1\0`) or a sequence file (FASTA/FASTQ, plain or `.gz`) and archives
/// it with the explicit sample `name` (never copying input file names).
pub fn archive_auto(
    src: &Path,
    dst: &Path,
    k_override: Option<usize>,
    name: &str,
    modality_override: Option<InputModality>,
    drop_counts: bool,
) -> Result<ArchiveSummary> {
    let mut magic = [0u8; 8];
    let n = File::open(src)
        .with_context(|| format!("open {}", src.display()))?
        .read(&mut magic)?;
    let is_table =
        n == 8 && (&magic == b"GRMKSET1" || &magic == MAGIC_KTAB || &magic == MAGIC_GRKSET);
    if is_table {
        if let Some(req_k) = k_override {
            let reader = AnyTableReader::open(src)?;
            ensure!(
                req_k == reader.k(),
                "{}: --k {req_k} disagrees with input table k={}",
                src.display(),
                reader.k()
            );
        }
        archive_from_table(src, dst, name, modality_override, drop_counts)
    } else {
        let k = k_override.unwrap_or(31);
        let modality = modality_override.unwrap_or(InputModality::Unknown);
        archive_from_fasta(src, dst, k, name, modality)
    }
}

/// Parses and validates all header fields, arithmetic invariants, and payload
/// byte bounds against `file_len` BEFORE any dynamic allocation. Shared by
/// both `read_spss_header` (`inspect`) and `read_spss_body` (`expand`).
fn parse_spss_header_and_layout(
    path: &Path,
    file: &File,
    file_len: u64,
) -> Result<ParsedSpssLayout> {
    let min_file_len = (SPSS_HEADER_LEN as u64) + 8;
    ensure!(
        file_len >= min_file_len,
        "{} is truncated before GRMSPS1 header",
        path.display()
    );
    let mut head = [0u8; SPSS_HEADER_LEN];
    file.read_exact_at(&mut head, 0)
        .with_context(|| format!("read SPSS header of {}", path.display()))?;
    ensure!(
        &head[..8] == MAGIC_SPSS1,
        "{} is not a GRMSPS1 archive",
        path.display()
    );
    let k = head[8] as usize;
    ensure!(
        (MIN_K..=MAX_K).contains(&k),
        "{}: unsupported SPSS k {k}",
        path.display()
    );
    let modality = InputModality::from_byte(head[9])
        .with_context(|| format!("{}: invalid SPSS modality", path.display()))?;
    ensure!(
        head[10..12] == [0u8; 2],
        "{}: non-zero reserved header bytes in GRMSPS1",
        path.display()
    );

    let name_len_u32 = u32::from_le_bytes(head[12..16].try_into()?);
    let n_kmers = u64::from_le_bytes(head[16..24].try_into()?);
    let n_simplitigs = u64::from_le_bytes(head[24..32].try_into()?);
    let total_bases = u64::from_le_bytes(head[32..40].try_into()?);
    let fnv64_digest = u64::from_le_bytes(head[40..48].try_into()?);
    let mut digest_sha256 = [0u8; 32];
    digest_sha256.copy_from_slice(&head[SPSS_DIGEST_START..SPSS_DIGEST_END]);

    if n_kmers == 0 {
        ensure!(
            n_simplitigs == 0 && total_bases == 0,
            "{}: empty SPSS archive has non-zero simplitig or base counts",
            path.display()
        );
    } else {
        ensure!(
            n_simplitigs >= 1 && n_simplitigs <= n_kmers,
            "{}: invalid n_simplitigs ({n_simplitigs}) for n_kmers ({n_kmers})",
            path.display()
        );
    }

    let overhead_bases = (k as u64 - 1)
        .checked_mul(n_simplitigs)
        .with_context(|| format!("{}: SPSS simplitig overhead overflow", path.display()))?;
    let expected_bases = n_kmers
        .checked_add(overhead_bases)
        .with_context(|| format!("{}: SPSS total_bases overflow", path.display()))?;
    ensure!(
        total_bases == expected_bases,
        "{}: corrupt SPSS counts (total_bases={total_bases}, expected {expected_bases})",
        path.display()
    );

    let packed_len_u64 = total_bases.div_ceil(4);
    let len_field_off = (SPSS_HEADER_LEN as u64)
        .checked_add(name_len_u32 as u64)
        .with_context(|| format!("{}: SPSS name length overflow", path.display()))?;
    let body_start = len_field_off
        .checked_add(8)
        .with_context(|| format!("{}: SPSS body offset overflow", path.display()))?;
    ensure!(
        body_start <= file_len,
        "{}: truncated before SPSS payload (need at least {body_start} bytes, file has {file_len})",
        path.display()
    );

    let body_remaining = file_len - body_start;
    let min_body = n_simplitigs
        .checked_add(packed_len_u64)
        .with_context(|| format!("{}: SPSS minimum body length overflow", path.display()))?;
    ensure!(
        min_body <= body_remaining,
        "{}: header declares {n_simplitigs} simplitigs and {packed_len_u64} packed bytes, exceeding {body_remaining} remaining bytes",
        path.display()
    );

    let name_len = usize::try_from(name_len_u32).context("SPSS name length overflow")?;
    let mut name_bytes = alloc_zeroed_checked(name_len, "GRMSPS1 name")?;
    if !name_bytes.is_empty() {
        file.read_exact_at(&mut name_bytes, SPSS_HEADER_LEN as u64)?;
    }
    let name = String::from_utf8(name_bytes)
        .with_context(|| format!("{}: sample name is not valid UTF-8", path.display()))?;

    let mut l_len_buf = [0u8; 8];
    file.read_exact_at(&mut l_len_buf, len_field_off)?;
    let lengths_byte_len_u64 = u64::from_le_bytes(l_len_buf);
    if n_simplitigs == 0 {
        ensure!(
            lengths_byte_len_u64 == 0 && packed_len_u64 == 0,
            "{}: empty SPSS archive has non-zero payload length",
            path.display()
        );
    } else {
        ensure!(
            lengths_byte_len_u64 >= n_simplitigs,
            "{}: simplitig length stream ({lengths_byte_len_u64} bytes) too short for {n_simplitigs} simplitigs",
            path.display()
        );
    }

    let packed_off = body_start
        .checked_add(lengths_byte_len_u64)
        .with_context(|| format!("{}: SPSS packed offset overflow", path.display()))?;
    let expected_total = packed_off
        .checked_add(packed_len_u64)
        .with_context(|| format!("{}: SPSS total file size overflow", path.display()))?;
    ensure!(
        file_len == expected_total,
        "{}: file size {file_len} does not match expected {expected_total}",
        path.display()
    );

    let _ = usize::try_from(n_kmers).context("SPSS n_kmers exceeds addressable memory")?;
    let _ =
        usize::try_from(n_simplitigs).context("SPSS n_simplitigs exceeds addressable memory")?;
    let lengths_byte_len =
        usize::try_from(lengths_byte_len_u64).context("SPSS lengths_byte_len overflow")?;
    let packed_len = usize::try_from(packed_len_u64).context("SPSS packed_len overflow")?;

    Ok(ParsedSpssLayout {
        header: SpssHeader {
            k,
            modality,
            name,
            n_kmers,
            n_simplitigs,
            total_bases,
            fnv64_digest,
            digest_sha256,
            archive_bytes: file_len,
        },
        lengths_off: body_start,
        lengths_byte_len,
        packed_off,
        packed_len,
    })
}

/// Reads and structurally validates header metadata and file bounds from an
/// `spss-archive` (`GRMSPS1\0`) file without verifying the body SHA-256 or
/// expanding the packed bases.
pub fn read_spss_header(path: &Path) -> Result<SpssHeader> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let file_len = file.metadata()?.len();
    let layout = parse_spss_header_and_layout(path, &file, file_len)?;
    Ok(layout.header)
}

fn read_spss_body(path: &Path) -> Result<(SpssHeader, Vec<usize>, Vec<u8>)> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let file_len = file.metadata()?.len();
    let layout = parse_spss_header_and_layout(path, &file, file_len)?;

    let computed = compute_spss_sha256(&file, file_len)?;
    ensure!(
        computed == layout.header.digest_sha256,
        "{}: GRMSPS1 SHA-256 integrity check failed: header has {}, computed {}",
        path.display(),
        hex_bytes(&layout.header.digest_sha256),
        hex_bytes(&computed)
    );

    let mut len_bytes = alloc_zeroed_checked(layout.lengths_byte_len, "GRMSPS1 lengths")?;
    if !len_bytes.is_empty() {
        file.read_exact_at(&mut len_bytes, layout.lengths_off)?;
    }

    let n_simplitigs = usize::try_from(layout.header.n_simplitigs)
        .context("SPSS n_simplitigs exceeds addressable memory")?;
    let mut lengths = Vec::new();
    lengths
        .try_reserve(n_simplitigs)
        .context("SPSS simplitig length list allocation failed")?;

    let k = layout.header.k;
    let mut pos = 0usize;
    let mut sum_bases = 0u64;
    for _ in 0..n_simplitigs {
        let extra = get_uleb(&len_bytes, &mut pos)?;
        let extra_usize = usize::try_from(extra).context("simplitig length overflow")?;
        let l = extra_usize
            .checked_add(k)
            .context("simplitig length overflow")?;
        sum_bases = sum_bases
            .checked_add(l as u64)
            .context("SPSS decoded bases overflow")?;
        ensure!(
            sum_bases <= layout.header.total_bases,
            "{}: decoded simplitig bases ({sum_bases}) exceed header total_bases ({})",
            path.display(),
            layout.header.total_bases
        );
        lengths.push(l);
    }
    ensure!(
        pos == len_bytes.len() && sum_bases == layout.header.total_bases,
        "{}: simplitig lengths do not match total_bases",
        path.display()
    );

    let mut packed_bases = alloc_zeroed_checked(layout.packed_len, "GRMSPS1 packed bases")?;
    if !packed_bases.is_empty() {
        file.read_exact_at(&mut packed_bases, layout.packed_off)?;
    }
    let rem_bases = (layout.header.total_bases % 4) as u32;
    if rem_bases != 0 {
        let valid_mask = (1u8 << (rem_bases * 2)) - 1;
        ensure!(
            (packed_bases[layout.packed_len - 1] & !valid_mask) == 0,
            "{}: non-zero trailing padding bits in GRMSPS1 packed bases",
            path.display()
        );
    }

    Ok((layout.header, lengths, packed_bases))
}

/// Expands a verified `spss-archive` (`GRMSPS1\0`) into sorted mixed keys in
/// parallel across the 1,024 output buckets, verifying SHA-256, strict
/// uniqueness, `n_kmers`, and `fnv64`.
pub fn expand_spss_keys<K: KmerKey>(path: &Path) -> Result<(SpssHeader, Vec<K>)> {
    let (header, lengths, packed_bases) = read_spss_body(path)?;
    let k = header.k;
    let want_width = if k <= MAX_K_U64 { 8 } else { 16 };
    ensure!(
        K::WIDTH_BYTES >= want_width,
        "cannot decode SPSS k={k} into {}-byte key",
        K::WIDTH_BYTES
    );

    let space = K::new_space(k);
    let mask: u128 = if 2 * k == 128 {
        u128::MAX
    } else {
        (1u128 << (2 * k)) - 1
    };
    let shift = 2 * (k as u32 - 1);

    let mut buckets: Vec<Vec<K>> = (0..BUCKETS).map(|_| Vec::new()).collect();
    let mut base_off = 0usize;
    for &slen in &lengths {
        let mut fwd = 0u128;
        let mut rev = 0u128;
        for j in 0..slen {
            let idx = base_off + j;
            let c = ((packed_bases[idx >> 2] >> ((idx & 3) * 2)) & 3) as u128;
            fwd = ((fwd << 2) | c) & mask;
            rev = (rev >> 2) | ((3 - c) << shift);
            if j + 1 >= k {
                let canon = K::from_u128(fwd.min(rev));
                let mixed = K::mix(&space, canon);
                let b = mixed.bucket_of(k, BUCKET_BITS);
                buckets[b].push(mixed);
            }
        }
        base_off += slen;
    }

    // Sort all 1,024 buckets in parallel and verify strict uniqueness.
    buckets
        .par_iter_mut()
        .try_for_each(|bucket| -> Result<()> {
            bucket.sort_unstable();
            ensure!(
                bucket.windows(2).all(|w| w[0] < w[1]),
                "duplicate k-mer decoded from SPSS archive"
            );
            Ok(())
        })?;

    let n_kmers_usize =
        usize::try_from(header.n_kmers).context("SPSS n_kmers exceeds addressable memory")?;
    let mut all = Vec::new();
    all.try_reserve(n_kmers_usize)
        .context("SPSS expanded key buffer allocation failed")?;
    let mut fnv = Fnv::new();
    for bucket in buckets {
        for y in bucket {
            y.with_le_bytes(|b| fnv.bytes(b));
            all.push(y);
        }
    }
    ensure!(
        all.len() as u64 == header.n_kmers,
        "decoded {} k-mers, expected {}",
        all.len(),
        header.n_kmers
    );
    ensure!(
        fnv.0 == header.fnv64_digest,
        "{}: SPSS digest mismatch (expected {:016x}, got {})",
        path.display(),
        header.fnv64_digest,
        fnv.hex()
    );
    Ok((header, all))
}

/// Expands a verified `spss-archive` (`GRMSPS1\0`) with `k <= 32` back into a
/// `.kset` v1 (`GRMKSET1`) file, checking the `.kset` extension, SHA-256
/// digest, and `fnv64` digest.
pub fn expand_to_kset_v1(src: &Path, dst: &Path) -> Result<ExpandSummary> {
    let t0 = Instant::now();
    check_output_extension(dst, TableFormat::KsetV1)?;
    let header_preview = read_spss_header(src)?;
    ensure!(
        header_preview.k <= MAX_K_U64,
        "cannot expand k={} into GRMKSET1 (max k is {MAX_K_U64}); use --format grk or --format ktab",
        header_preview.k
    );
    let (header, values) = expand_spss_keys::<u64>(src)?;
    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let space = KmerSpace::new(header.k);
    kset::write(dst, &header.name, &space, &values)?;
    let output_bytes = fs::metadata(dst)?.len();
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(ExpandSummary {
        format: TableFormat::KsetV1.magic_label().to_string(),
        sample: header.name,
        k: header.k,
        modality: None,
        n_kmers: header.n_kmers,
        n_simplitigs: header.n_simplitigs,
        total_bases: header.total_bases,
        output_bytes,
        fnv64: format!("{:016x}", header.fnv64_digest),
        header_sha256: hex_bytes(&header.digest_sha256),
        elapsed_ms,
    })
}

/// Expands a verified `spss-archive` (`GRMSPS1\0`) for any `k = 5..=64` back
/// into a verified Golomb-Rice `.grk` (`GRMKGR1\0`) table, checking the `.grk`
/// extension, SHA-256 digest, and `fnv64` digest.
pub fn expand_to_gr(src: &Path, dst: &Path) -> Result<ExpandSummary> {
    let t0 = Instant::now();
    check_output_extension(dst, TableFormat::GrKset)?;
    let header_preview = read_spss_header(src)?;
    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let header = if header_preview.k <= MAX_K_U64 {
        let (h, vals) = expand_spss_keys::<u64>(src)?;
        write_gr_kset(dst, &h.name, h.k, h.modality, &vals)?;
        h
    } else {
        let (h, vals) = expand_spss_keys::<u128>(src)?;
        write_gr_kset(dst, &h.name, h.k, h.modality, &vals)?;
        h
    };
    let output_bytes = fs::metadata(dst)?.len();
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(ExpandSummary {
        format: TableFormat::GrKset.magic_label().to_string(),
        sample: header.name,
        k: header.k,
        modality: Some(header.modality),
        n_kmers: header.n_kmers,
        n_simplitigs: header.n_simplitigs,
        total_bases: header.total_bases,
        output_bytes,
        fnv64: format!("{:016x}", header.fnv64_digest),
        header_sha256: hex_bytes(&header.digest_sha256),
        elapsed_ms,
    })
}

/// Expands a verified `spss-archive` (`GRMSPS1\0`) for any `k = 5..=64` into a
/// presence-only `.ktab` (`GRMKTB1\0`) file (every k-mer count set to `1`),
/// checking the `.ktab` extension, SHA-256 digest, and `fnv64` digest.
pub fn expand_to_ktab(src: &Path, dst: &Path) -> Result<ExpandSummary> {
    let t0 = Instant::now();
    check_output_extension(dst, TableFormat::Ktab)?;
    let header_preview = read_spss_header(src)?;
    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let header = if header_preview.k <= MAX_K_U64 {
        let (h, vals) = expand_spss_keys::<u64>(src)?;
        let entries: Vec<(u64, u32)> = vals.into_iter().map(|v| (v, 1)).collect();
        write_ktab(dst, &h.name, h.k, h.modality, 1, 0, &entries)?;
        h
    } else {
        let (h, vals) = expand_spss_keys::<u128>(src)?;
        let entries: Vec<(u128, u32)> = vals.into_iter().map(|v| (v, 1)).collect();
        write_ktab(dst, &h.name, h.k, h.modality, 1, 0, &entries)?;
        h
    };
    let output_bytes = fs::metadata(dst)?.len();
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(ExpandSummary {
        format: TableFormat::Ktab.magic_label().to_string(),
        sample: header.name,
        k: header.k,
        modality: Some(header.modality),
        n_kmers: header.n_kmers,
        n_simplitigs: header.n_simplitigs,
        total_bases: header.total_bases,
        output_bytes,
        fnv64: format!("{:016x}", header.fnv64_digest),
        header_sha256: hex_bytes(&header.digest_sha256),
        elapsed_ms,
    })
}

/// Verifies an `spss-archive` (`GRMSPS1\0`) by checking its 32-byte SHA-256
/// digest, expanding every simplitig, and checking uniqueness, `n_kmers`, and
/// the `fnv64` digest.
pub fn verify_spss(path: &Path) -> Result<SpssHeader> {
    Ok(verify_spss_with_digests(path, None)?.header)
}

/// Verifies an `spss-archive` (`GRMSPS1\0`) and also returns `keys_sha256` and
/// an optional bottom-`s` MinHash sketch prefix.
pub fn verify_spss_with_digests(path: &Path, sketch: Option<usize>) -> Result<SpssVerifyResult> {
    let header = read_spss_header(path)?;
    let mut sha = Sha256Context::new(&SHA256);
    if header.k <= MAX_K_U64 {
        let (h, vals) = expand_spss_keys::<u64>(path)?;
        for &v in &vals {
            v.with_le_bytes(|b| sha.update(b));
        }
        let sketch_hex = sketch.map(|s| {
            vals.iter()
                .take(s)
                .map(|v| v.to_hex())
                .collect::<Vec<String>>()
        });
        Ok(SpssVerifyResult {
            header: h,
            keys_sha256: hex_bytes(sha.finish().as_ref()),
            sketch_hex,
        })
    } else {
        let (h, vals) = expand_spss_keys::<u128>(path)?;
        for &v in &vals {
            v.with_le_bytes(|b| sha.update(b));
        }
        let sketch_hex = sketch.map(|s| {
            vals.iter()
                .take(s)
                .map(|v| v.to_hex())
                .collect::<Vec<String>>()
        });
        Ok(SpssVerifyResult {
            header: h,
            keys_sha256: hex_bytes(sha.finish().as_ref()),
            sketch_hex,
        })
    }
}

/// Returns `Ok(true)` if `path` begins with `GRMSPS1\0`.
pub fn is_spss_file(path: &Path) -> Result<bool> {
    let mut f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut magic = [0u8; 8];
    let n = f.read(&mut magic)?;
    Ok(n == 8 && &magic == MAGIC_SPSS1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn spss_archive_and_expand_match_kset_v1_bit_for_bit() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(909);
        let mut seq: Vec<u8> = (0..120_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        // Inject a repeated segment and N break to test branching simplitigs.
        let rep = seq[10_000..11_000].to_vec();
        seq.extend_from_slice(b"NNNN");
        seq.extend_from_slice(&rep);
        let fna = dir.path().join("genome.fna");
        fs::write(&fna, [b">c1\n".as_slice(), &seq, b"\n"].concat()).unwrap();

        let space = KmerSpace::new(31);
        let (values, _) = kset::count_file(&fna, &space, 1).unwrap();
        let v1_path = dir.path().join("orig.kset");
        kset::write(&v1_path, "sample", &space, &values).unwrap();
        let orig_digest = kset::digest(&values);
        let orig_bytes = fs::read(&v1_path).unwrap();

        let spss_path = dir.path().join("sample.spss");
        let arch = archive_from_table(&v1_path, &spss_path, "sample", None, false).unwrap();
        assert_eq!(arch.fnv64, orig_digest);
        assert!(!arch.counts_dropped);
        assert_eq!(arch.header_sha256.len(), 64);
        assert!(
            arch.bytes_per_kmer < 0.30,
            "expected < 0.30 B/k-mer on assembly, got {:.4}",
            arch.bytes_per_kmer
        );

        let rebuilt_path = dir.path().join("rebuilt.kset");
        let exp = expand_to_kset_v1(&spss_path, &rebuilt_path).unwrap();
        assert_eq!(exp.fnv64, orig_digest);
        assert_eq!(exp.header_sha256, arch.header_sha256);
        let rebuilt_bytes = fs::read(&rebuilt_path).unwrap();
        assert_eq!(rebuilt_bytes, orig_bytes);
    }

    #[test]
    fn spss_empty_set_round_trips_for_u64_and_u128() {
        let dir = tempfile::tempdir().unwrap();
        for (k, modality) in [
            (31usize, InputModality::Assembly),
            (61, InputModality::Reads),
        ] {
            let spss = dir.path().join(format!("empty_{k}.spss"));
            let arch = if k <= MAX_K_U64 {
                write_spss_from_mixed::<u64>(&spss, "sample", k, modality, &[]).unwrap()
            } else {
                write_spss_from_mixed::<u128>(&spss, "sample", k, modality, &[]).unwrap()
            };
            assert_eq!(arch.n_kmers, 0);
            assert_eq!(arch.n_simplitigs, 0);
            assert_eq!(arch.total_bases, 0);
            let hdr = verify_spss(&spss).unwrap();
            assert_eq!(hdr.n_kmers, 0);
            assert_eq!(hdr.modality, modality);

            let gr = dir.path().join(format!("empty_{k}.grk"));
            let exp_gr = expand_to_gr(&spss, &gr).unwrap();
            assert_eq!(exp_gr.n_kmers, 0);
            let tab = dir.path().join(format!("empty_{k}.ktab"));
            let exp_tab = expand_to_ktab(&spss, &tab).unwrap();
            assert_eq!(exp_tab.n_kmers, 0);
            if k <= MAX_K_U64 {
                let v1 = dir.path().join(format!("empty_{k}.kset"));
                let exp_v1 = expand_to_kset_v1(&spss, &v1).unwrap();
                assert_eq!(exp_v1.n_kmers, 0);
            }
        }
    }
}
