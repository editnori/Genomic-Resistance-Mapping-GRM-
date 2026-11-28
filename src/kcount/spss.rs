//! Experimental 2-bit packed Spectrum-Preserving String Set (SPSS / greedy
//! simplitig) archive format (`spss-archive`, magic `GRMSPS1\0`).
//!
//! Sorted mixed 62-bit k-mer integers (~5M per bacterial genome) have an
//! information-theoretic lower bound of `log2(2^62 / N) + log2(e) ~= 41.2`
//! bits (~5.14 B/k-mer). Exploiting `(k - 1)`-base overlaps between adjacent
//! canonical k-mers collapses a bacterial k-mer spectrum into greedy
//! simplitigs packed at 2 bits per base (~0.25 to 0.40 B/k-mer), while
//! leaving `.kset` v1 (`GRMKSET1`) untouched as the fast merge form.

use crate::dna::KmerSpace;
use crate::fnv::Fnv;
use crate::kcount::format::{
    AnyTableReader, MAGIC_GRKSET, MAGIC_KTAB, get_uleb, put_uleb, write_ktab,
};
use crate::kcount::key::{KmerKey, KmerSpace64, KmerSpace128, MAX_K, MAX_K_U64, MIN_K};
use crate::kcount::superkmer::StreamingSeqReader;
use crate::kset::{self, BUCKET_BITS, BUCKETS};
use anyhow::{Context, Result, ensure};
use rayon::prelude::*;
use rustc_hash::{FxBuildHasher, FxHashSet};
use serde::{Deserialize, Serialize};
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::time::Instant;

pub const MAGIC_SPSS1: &[u8; 8] = b"GRMSPS1\0";

/// Summary returned by [`archive_from_table`] and [`archive_from_fasta`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ArchiveSummary {
    pub sample: String,
    pub k: usize,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub archive_bytes: u64,
    pub bytes_per_kmer: f64,
    pub fnv64: String,
    pub elapsed_ms: f64,
}

/// Summary returned by [`expand_to_kset_v1`] and [`expand_to_ktab`].
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ExpandSummary {
    pub sample: String,
    pub k: usize,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub output_bytes: u64,
    pub fnv64: String,
    pub elapsed_ms: f64,
}

/// Decoded header metadata for an `spss-archive` (`GRMSPS1\0`) file.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct SpssHeader {
    pub k: usize,
    pub name: String,
    pub n_kmers: u64,
    pub n_simplitigs: u64,
    pub total_bases: u64,
    pub fnv64_digest: u64,
    pub archive_bytes: u64,
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
/// canonical k-mers. Output is 100% deterministic because seeds and extensions
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
    n_kmers: u64,
    fnv64_digest: u64,
    lengths: &'a [usize],
    packed_bases: &'a [u8],
    total_bases: u64,
}

fn write_spss_file(dst: &Path, p: &SpssPayload<'_>) -> Result<u64> {
    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let mut len_bytes = Vec::with_capacity(p.lengths.len() * 2);
    for &l in p.lengths {
        ensure!(l >= p.k, "simplitig length {l} shorter than k={}", p.k);
        put_uleb(&mut len_bytes, (l - p.k) as u64);
    }

    let part = crate::confine::suffixed(dst, ".part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(MAGIC_SPSS1)?;
    w.write_all(&[p.k as u8, 0, 0, 0])?;
    w.write_all(&(p.name.len() as u32).to_le_bytes())?;
    w.write_all(&p.n_kmers.to_le_bytes())?;
    w.write_all(&(p.lengths.len() as u64).to_le_bytes())?;
    w.write_all(&p.total_bases.to_le_bytes())?;
    w.write_all(&p.fnv64_digest.to_le_bytes())?;
    w.write_all(p.name.as_bytes())?;
    w.write_all(&(len_bytes.len() as u64).to_le_bytes())?;
    w.write_all(&len_bytes)?;
    w.write_all(p.packed_bases)?;
    w.flush()?;
    drop(w);
    fs::rename(&part, dst).with_context(|| format!("rename to {}", dst.display()))?;
    Ok(fs::metadata(dst)?.len())
}

/// Archives sorted mixed keys directly into a `GRMSPS1\0` SPSS archive.
pub fn write_spss_from_mixed<K: KmerKey>(
    dst: &Path,
    name: &str,
    k: usize,
    mixed_sorted: &[K],
) -> Result<ArchiveSummary> {
    let t0 = Instant::now();
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    let space = K::new_space(k);
    let mut fnv = Fnv::new();
    let mut raw = Vec::with_capacity(mixed_sorted.len());
    for &y in mixed_sorted {
        y.with_le_bytes(|b| fnv.bytes(b));
        raw.push(K::unmix(&space, y));
    }
    let mut scratch = Vec::new();
    K::radix_sort(&mut raw, &mut scratch);
    drop(scratch);
    debug_assert!(raw.windows(2).all(|w| w[0] < w[1]));

    let (lengths, packed_bases, total_bases) = build_simplitigs_generic(k, &raw);
    let n_kmers = raw.len() as u64;
    let n_simplitigs = lengths.len() as u64;
    let total_bases_u64 = total_bases as u64;
    ensure!(
        total_bases_u64 == n_kmers + (k as u64 - 1) * n_simplitigs,
        "SPSS invariant violated: total_bases={total_bases_u64}, n_kmers={n_kmers}, n_simplitigs={n_simplitigs}"
    );

    let archive_bytes = write_spss_file(
        dst,
        &SpssPayload {
            name,
            k,
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
        n_kmers,
        n_simplitigs,
        total_bases: total_bases_u64,
        archive_bytes,
        bytes_per_kmer,
        fnv64: fnv.hex(),
        elapsed_ms,
    })
}

/// Archives any supported table (`.kset` v1, `.ktab`, or `.grk`) into an
/// `spss-archive` (`GRMSPS1\0`).
pub fn archive_from_table(src: &Path, dst: &Path) -> Result<ArchiveSummary> {
    let reader = AnyTableReader::open(src)?;
    let k = reader.k();
    let name = reader.name().to_string();
    if k <= MAX_K_U64 {
        let mut keys = Vec::with_capacity(reader.n_kmers() as usize);
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        let mut counts = Vec::new();
        for b in 0..BUCKETS {
            counts.clear();
            reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        }
        write_spss_from_mixed(dst, &name, k, &keys)
    } else {
        let mut keys = Vec::with_capacity(reader.n_kmers() as usize);
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        let mut counts = Vec::new();
        for b in 0..BUCKETS {
            counts.clear();
            reader.read_bucket::<u128>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        }
        write_spss_from_mixed(dst, &name, k, &keys)
    }
}

/// Archives one FASTA/FASTQ file (plain or `.gz`) directly into an
/// `spss-archive` (`GRMSPS1\0`) at the requested `k`.
pub fn archive_from_fasta(src: &Path, dst: &Path, k: usize, name: &str) -> Result<ArchiveSummary> {
    let t0 = Instant::now();
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    if k <= MAX_K_U64 {
        let space = KmerSpace64::new(k);
        let mut raw = Vec::new();
        let mut reader = StreamingSeqReader::open(src)?;
        reader.for_each_record(|_, seq, _| {
            space.for_each_canonical(seq, |x| raw.push(x));
            Ok(())
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
        let (lengths, packed_bases, total_bases) = build_simplitigs_generic(k, &raw);
        let n_kmers = raw.len() as u64;
        let n_simplitigs = lengths.len() as u64;
        let total_bases_u64 = total_bases as u64;
        let archive_bytes = write_spss_file(
            dst,
            &SpssPayload {
                name,
                k,
                n_kmers,
                fnv64_digest: fnv.0,
                lengths: &lengths,
                packed_bases: &packed_bases,
                total_bases: total_bases_u64,
            },
        )?;
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
        Ok(ArchiveSummary {
            sample: name.to_string(),
            k,
            n_kmers,
            n_simplitigs,
            total_bases: total_bases_u64,
            archive_bytes,
            bytes_per_kmer: if n_kmers == 0 {
                0.0
            } else {
                archive_bytes as f64 / n_kmers as f64
            },
            fnv64: fnv.hex(),
            elapsed_ms,
        })
    } else {
        let space = KmerSpace128::new(k);
        let mut raw = Vec::new();
        let mut reader = StreamingSeqReader::open(src)?;
        reader.for_each_record(|_, seq, _| {
            space.for_each_canonical(seq, |x| raw.push(x));
            Ok(())
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
        let (lengths, packed_bases, total_bases) = build_simplitigs_generic(k, &raw);
        let n_kmers = raw.len() as u64;
        let n_simplitigs = lengths.len() as u64;
        let total_bases_u64 = total_bases as u64;
        let archive_bytes = write_spss_file(
            dst,
            &SpssPayload {
                name,
                k,
                n_kmers,
                fnv64_digest: fnv.0,
                lengths: &lengths,
                packed_bases: &packed_bases,
                total_bases: total_bases_u64,
            },
        )?;
        let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
        Ok(ArchiveSummary {
            sample: name.to_string(),
            k,
            n_kmers,
            n_simplitigs,
            total_bases: total_bases_u64,
            archive_bytes,
            bytes_per_kmer: if n_kmers == 0 {
                0.0
            } else {
                archive_bytes as f64 / n_kmers as f64
            },
            fnv64: fnv.hex(),
            elapsed_ms,
        })
    }
}

/// Automatically detects whether `src` is a k-mer table (`GRMKSET1`,
/// `GRMKTB1\0`, `GRMKGR1\0`) or a sequence file (FASTA/FASTQ) and archives it.
pub fn archive_auto(
    src: &Path,
    dst: &Path,
    k_override: Option<usize>,
    name_override: Option<&str>,
) -> Result<ArchiveSummary> {
    let mut magic = [0u8; 8];
    let n = File::open(src)
        .with_context(|| format!("open {}", src.display()))?
        .read(&mut magic)?;
    let is_table =
        n == 8 && (&magic == b"GRMKSET1" || &magic == MAGIC_KTAB || &magic == MAGIC_GRKSET);
    if is_table {
        archive_from_table(src, dst)
    } else {
        let k = k_override.unwrap_or(31);
        let fallback_name = src.file_stem().and_then(|s| s.to_str()).unwrap_or("sample");
        let name = name_override.unwrap_or(fallback_name);
        archive_from_fasta(src, dst, k, name)
    }
}

/// Reads header metadata from an `spss-archive` (`GRMSPS1\0`) file without
/// expanding the packed bases.
pub fn read_spss_header(path: &Path) -> Result<SpssHeader> {
    let mut f = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let archive_bytes = f.metadata()?.len();
    let mut head = [0u8; 48];
    f.read_exact(&mut head)
        .with_context(|| format!("read SPSS header of {}", path.display()))?;
    ensure!(
        &head[..8] == MAGIC_SPSS1,
        "{} is not a GRMSPS1 archive",
        path.display()
    );
    let k = head[8] as usize;
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported SPSS k {k}");
    let name_len = u32::from_le_bytes(head[12..16].try_into()?) as usize;
    let n_kmers = u64::from_le_bytes(head[16..24].try_into()?);
    let n_simplitigs = u64::from_le_bytes(head[24..32].try_into()?);
    let total_bases = u64::from_le_bytes(head[32..40].try_into()?);
    let fnv64_digest = u64::from_le_bytes(head[40..48].try_into()?);
    let mut name_bytes = vec![0u8; name_len];
    f.read_exact(&mut name_bytes)?;
    let name = String::from_utf8(name_bytes)?;
    Ok(SpssHeader {
        k,
        name,
        n_kmers,
        n_simplitigs,
        total_bases,
        fnv64_digest,
        archive_bytes,
    })
}

fn read_spss_body(path: &Path) -> Result<(SpssHeader, Vec<usize>, Vec<u8>)> {
    let mut r =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let archive_bytes = fs::metadata(path)?.len();
    let mut head = [0u8; 48];
    r.read_exact(&mut head)
        .with_context(|| format!("read SPSS header of {}", path.display()))?;
    ensure!(
        &head[..8] == MAGIC_SPSS1,
        "{} is not a GRMSPS1 archive",
        path.display()
    );
    let k = head[8] as usize;
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported SPSS k {k}");
    let name_len = u32::from_le_bytes(head[12..16].try_into()?) as usize;
    let n_kmers = u64::from_le_bytes(head[16..24].try_into()?);
    let n_simplitigs = u64::from_le_bytes(head[24..32].try_into()?);
    let total_bases = u64::from_le_bytes(head[32..40].try_into()?);
    let fnv64_digest = u64::from_le_bytes(head[40..48].try_into()?);

    ensure!(
        total_bases == n_kmers + (k as u64 - 1) * n_simplitigs,
        "{}: corrupt SPSS counts",
        path.display()
    );

    let mut name_bytes = vec![0u8; name_len];
    r.read_exact(&mut name_bytes)?;
    let name = String::from_utf8(name_bytes)?;

    let mut l_len_buf = [0u8; 8];
    r.read_exact(&mut l_len_buf)?;
    let lengths_byte_len = u64::from_le_bytes(l_len_buf) as usize;
    let packed_len = (total_bases as usize).div_ceil(4);
    let expected_total = 48u64 + name_len as u64 + 8 + lengths_byte_len as u64 + packed_len as u64;
    ensure!(
        archive_bytes == expected_total,
        "{}: file size {archive_bytes} does not match expected {expected_total}",
        path.display()
    );

    let mut len_bytes = vec![0u8; lengths_byte_len];
    r.read_exact(&mut len_bytes)?;
    let mut pos = 0usize;
    let mut lengths = Vec::with_capacity(n_simplitigs as usize);
    let mut sum_bases = 0u64;
    for _ in 0..n_simplitigs {
        let extra = get_uleb(&len_bytes, &mut pos)?;
        let l = (extra as usize)
            .checked_add(k)
            .context("simplitig length overflow")?;
        sum_bases = sum_bases.checked_add(l as u64).context("bases overflow")?;
        lengths.push(l);
    }
    ensure!(
        pos == len_bytes.len() && sum_bases == total_bases,
        "{}: simplitig lengths do not match total_bases",
        path.display()
    );

    let mut packed_bases = vec![0u8; packed_len];
    r.read_exact(&mut packed_bases)?;

    let header = SpssHeader {
        k,
        name,
        n_kmers,
        n_simplitigs,
        total_bases,
        fnv64_digest,
        archive_bytes,
    };
    Ok((header, lengths, packed_bases))
}

/// Expands an `spss-archive` (`GRMSPS1\0`) into sorted mixed keys in parallel
/// across the 1,024 output buckets, verifying exact `n_kmers` and `fnv64`.
pub fn expand_spss_keys<K: KmerKey>(path: &Path) -> Result<(SpssHeader, Vec<K>)> {
    let (header, lengths, packed_bases) = read_spss_body(path)?;
    let k = header.k;
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

    // Sort all 1,024 buckets in parallel and verify uniqueness.
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

    let mut all = Vec::with_capacity(header.n_kmers as usize);
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

/// Expands an `spss-archive` (`GRMSPS1\0`) with `k <= 32` losslessly back into
/// a `.kset` v1 (`GRMKSET1`) file, verifying the `fnv64` digest.
pub fn expand_to_kset_v1(src: &Path, dst: &Path) -> Result<ExpandSummary> {
    let t0 = Instant::now();
    let header_preview = read_spss_header(src)?;
    ensure!(
        header_preview.k <= MAX_K_U64,
        "cannot expand k={} into GRMKSET1 (max k is {MAX_K_U64}); use --format ktab",
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
        sample: header.name,
        k: header.k,
        n_kmers: header.n_kmers,
        n_simplitigs: header.n_simplitigs,
        total_bases: header.total_bases,
        output_bytes,
        fnv64: format!("{:016x}", header.fnv64_digest),
        elapsed_ms,
    })
}

/// Expands an `spss-archive` (`GRMSPS1\0`) for any `k = 5..=64` into a
/// presence-only `.ktab` (`GRMKTB1\0`) file, verifying the `fnv64` digest.
pub fn expand_to_ktab(src: &Path, dst: &Path) -> Result<ExpandSummary> {
    let t0 = Instant::now();
    let header_preview = read_spss_header(src)?;
    if let Some(parent) = dst.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    let header = if header_preview.k <= MAX_K_U64 {
        let (h, vals) = expand_spss_keys::<u64>(src)?;
        let entries: Vec<(u64, u32)> = vals.into_iter().map(|v| (v, 1)).collect();
        write_ktab(dst, &h.name, h.k, 1, &entries)?;
        h
    } else {
        let (h, vals) = expand_spss_keys::<u128>(src)?;
        let entries: Vec<(u128, u32)> = vals.into_iter().map(|v| (v, 1)).collect();
        write_ktab(dst, &h.name, h.k, 1, &entries)?;
        h
    };
    let output_bytes = fs::metadata(dst)?.len();
    let elapsed_ms = t0.elapsed().as_secs_f64() * 1000.0;
    Ok(ExpandSummary {
        sample: header.name,
        k: header.k,
        n_kmers: header.n_kmers,
        n_simplitigs: header.n_simplitigs,
        total_bases: header.total_bases,
        output_bytes,
        fnv64: format!("{:016x}", header.fnv64_digest),
        elapsed_ms,
    })
}

/// Verifies an `spss-archive` (`GRMSPS1\0`) by expanding every simplitig and
/// checking the reconstructed FNV-1a64 digest against its header.
pub fn verify_spss(path: &Path) -> Result<SpssHeader> {
    let header = read_spss_header(path)?;
    if header.k <= MAX_K_U64 {
        let (h, _) = expand_spss_keys::<u64>(path)?;
        Ok(h)
    } else {
        let (h, _) = expand_spss_keys::<u128>(path)?;
        Ok(h)
    }
}

/// Returns an error if `path` does not begin with `GRMSPS1\0`.
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
        kset::write(&v1_path, "g1", &space, &values).unwrap();
        let orig_digest = kset::digest(&values);
        let orig_bytes = fs::read(&v1_path).unwrap();

        let spss_path = dir.path().join("g1.spss");
        let arch = archive_from_table(&v1_path, &spss_path).unwrap();
        assert_eq!(arch.fnv64, orig_digest);
        assert!(
            arch.bytes_per_kmer < 0.35,
            "expected < 0.35 B/k-mer, got {:.4}",
            arch.bytes_per_kmer
        );

        let rebuilt_path = dir.path().join("rebuilt.kset");
        let exp = expand_to_kset_v1(&spss_path, &rebuilt_path).unwrap();
        assert_eq!(exp.fnv64, orig_digest);
        let rebuilt_bytes = fs::read(&rebuilt_path).unwrap();
        assert_eq!(rebuilt_bytes, orig_bytes);
    }
}
