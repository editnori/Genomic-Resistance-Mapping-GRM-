//! Per-genome k-mer sets: counting, sorting, and the `.kset` file format.
//!
//! A kset holds one genome's distinct canonical k-mers, mixed (see `dna`) and
//! sorted ascending. The file splits the value range into 2^BUCKET_BITS equal
//! buckets so the species merge can read the same bucket from every genome
//! without decoding whole files.
//!
//! Layout (little endian):
//!   magic "GRMKSET1" | k u8 | bucket_bits u8 | reserved u16 | n_kmers u64
//!   | name_len u32 | name | bucket byte offsets [(B+1) x u64] | counts [B x u32]
//!   | data: per bucket, LEB128 gaps from the bucket's lower bound.

use crate::dna::{KmerSpace, MAX_K, MIN_K};
use crate::fasta::{self, SeqStats, StatsBuilder};
use crate::fnv::Fnv;
use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

pub const BUCKET_BITS: u32 = 10;
pub const BUCKETS: usize = 1 << BUCKET_BITS;
const MAGIC: &[u8; 8] = b"GRMKSET1";

/// Counts one genome file. `min_abundance` 1 keeps every k-mer (assemblies);
/// a higher value drops k-mers seen fewer times (raw reads).
pub fn count_file(
    path: &Path,
    space: &KmerSpace,
    min_abundance: u32,
) -> Result<(Vec<u64>, SeqStats)> {
    let mut values = Vec::new();
    let stats = count_into(path, space, min_abundance, &mut values, &mut Vec::new())?;
    Ok((values, stats))
}

/// `count_file` into caller-owned buffers. Reusing them across genomes
/// avoids allocating and page-faulting tens of megabytes per genome, which
/// otherwise shows up as kernel time when many genomes count in parallel.
pub fn count_into(
    path: &Path,
    space: &KmerSpace,
    min_abundance: u32,
    values: &mut Vec<u64>,
    scratch: &mut Vec<u64>,
) -> Result<SeqStats> {
    values.clear();
    let mut stats = StatsBuilder::default();
    fasta::for_each_record(path, |_, seq| {
        stats.add(seq);
        space.for_each_canonical(seq, |x| values.push(space.mix(x)));
    })?;
    radix_sort(values, scratch);
    dedup_with_min_count(values, min_abundance);
    Ok(stats.finish())
}

/// Least-significant-digit radix sort over 8-bit digits. Digits that are the
/// same for every value are skipped, so 62-bit k-mers cost at most 8 passes.
pub fn radix_sort(values: &mut Vec<u64>, scratch: &mut Vec<u64>) {
    let n = values.len();
    if n < 256 {
        values.sort_unstable();
        return;
    }
    let mut hist = vec![[0usize; 256]; 8];
    for &v in values.iter() {
        for (d, h) in hist.iter_mut().enumerate() {
            h[((v >> (8 * d)) & 0xff) as usize] += 1;
        }
    }
    scratch.clear();
    scratch.resize(n, 0);
    for (d, h) in hist.iter().enumerate() {
        if h.contains(&n) {
            continue;
        }
        let mut offsets = [0usize; 256];
        let mut running = 0;
        for (o, &c) in offsets.iter_mut().zip(h.iter()) {
            *o = running;
            running += c;
        }
        let shift = 8 * d;
        for &v in values.iter() {
            let digit = ((v >> shift) & 0xff) as usize;
            scratch[offsets[digit]] = v;
            offsets[digit] += 1;
        }
        std::mem::swap(values, scratch);
    }
}

/// Collapses a sorted vector to distinct values that occur at least `min`
/// times.
pub fn dedup_with_min_count(values: &mut Vec<u64>, min: u32) {
    if min <= 1 {
        values.dedup();
        return;
    }
    let mut write = 0;
    let mut i = 0;
    while i < values.len() {
        let v = values[i];
        let mut j = i + 1;
        while j < values.len() && values[j] == v {
            j += 1;
        }
        if (j - i) as u32 >= min {
            values[write] = v;
            write += 1;
        }
        i = j;
    }
    values.truncate(write);
}

fn bucket_of(value: u64, bits: u32) -> usize {
    (value >> (bits - BUCKET_BITS)) as usize
}

fn bucket_floor(bucket: usize, bits: u32) -> u64 {
    (bucket as u64) << (bits - BUCKET_BITS)
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
            bail!("truncated varint")
        };
        *pos += 1;
        v |= ((byte & 0x7f) as u64) << shift;
        if byte < 0x80 {
            return Ok(v);
        }
        shift += 7;
        ensure!(shift < 64, "varint overflow");
    }
}

/// Writes sorted, distinct mixed k-mers. Writes to `<path>.part` and renames,
/// so an interrupted run never leaves a truncated kset behind.
pub fn write(path: &Path, name: &str, space: &KmerSpace, values: &[u64]) -> Result<()> {
    let (offsets, counts, data) = encode(space, values);
    let part = path.with_extension("kset.part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    write_body(&mut w, name, space, values.len(), &offsets, &counts, &data)?;
    // No fsync: a kset is a re-creatable cache, and the rename below
    // already keeps an interrupted run from leaving a truncated file.
    w.into_inner()?;
    std::fs::rename(&part, path)?;
    Ok(())
}

/// `write` for the `--disk-sparse` count, with an exact-length preflight
/// before creating the temporary file. The old kset, if any, remains until
/// the rename, so `bytes` is also the peak additional file length this write
/// needs. What the preflight returns (a reserve claim, say) is held until the
/// kset is renamed into place.
///
/// The temporary is created with `confine::create_part`, never through a
/// link, and the kset only replaces a regular file or nothing.
pub fn write_checked<G>(
    path: &Path,
    name: &str,
    space: &KmerSpace,
    values: &[u64],
    preflight: impl FnOnce(&Path, u64) -> Result<G>,
) -> Result<()> {
    let (offsets, counts, data) = encode(space, values);
    let _claim = preflight(path, encoded_len(name, data.len()))?;

    let part = path.with_extension("kset.part");
    crate::confine::regular_or_absent(path)?;
    let mut w = BufWriter::new(crate::confine::create_part(&part)?);
    let written = (|| -> Result<()> {
        write_body(&mut w, name, space, values.len(), &offsets, &counts, &data)?;
        // No fsync, as in `write`.
        w.flush()?;
        crate::confine::replace_with(&part, path)
    })();
    drop(w);
    if written.is_err() && crate::confine::regular_or_absent(&part).unwrap_or(false) {
        let _ = std::fs::remove_file(&part);
    }
    written
}

/// Each bucket's offset into the varint data, its k-mer count, and the data.
fn encode(space: &KmerSpace, values: &[u64]) -> (Vec<u64>, Vec<u32>, Vec<u8>) {
    debug_assert!(values.windows(2).all(|w| w[0] < w[1]));
    let bits = space.bits();
    let mut data = Vec::with_capacity(values.len() * 6);
    let mut offsets = vec![0u64; BUCKETS + 1];
    let mut counts = vec![0u32; BUCKETS];
    let mut i = 0;
    for b in 0..BUCKETS {
        offsets[b] = data.len() as u64;
        let mut prev = bucket_floor(b, bits);
        while i < values.len() && bucket_of(values[i], bits) == b {
            put_varint(&mut data, values[i] - prev);
            prev = values[i];
            counts[b] += 1;
            i += 1;
        }
    }
    offsets[BUCKETS] = data.len() as u64;
    (offsets, counts, data)
}

/// The kset file's bytes, in order, from `encode`'s parts.
fn write_body(
    w: &mut impl Write,
    name: &str,
    space: &KmerSpace,
    n_values: usize,
    offsets: &[u64],
    counts: &[u32],
    data: &[u8],
) -> Result<()> {
    w.write_all(MAGIC)?;
    w.write_all(&[space.k as u8, BUCKET_BITS as u8, 0, 0])?;
    w.write_all(&(n_values as u64).to_le_bytes())?;
    w.write_all(&(name.len() as u32).to_le_bytes())?;
    w.write_all(name.as_bytes())?;
    for o in offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for c in counts {
        w.write_all(&c.to_le_bytes())?;
    }
    w.write_all(data)?;
    Ok(())
}

/// Bytes of a kset file whose name is `name` and whose varint data takes
/// `data_len` bytes.
fn encoded_len(name: &str, data_len: usize) -> u64 {
    (24 + name.len() + (BUCKETS + 1) * 8 + BUCKETS * 4 + data_len) as u64
}

/// The longest a kset of `n_kmers` values, named `name`, can be: each value
/// is stored as its gap from the previous one (or its bucket's floor), and a
/// gap is below 2^(2k - BUCKET_BITS) <= 2^54, so its LEB128 takes at most 8
/// bytes. Disk preflights use this, not a typical size.
pub fn max_encoded_len(name: &str, n_kmers: u64) -> u64 {
    encoded_len(name, 0).saturating_add(n_kmers.saturating_mul(8))
}

/// FNV-1a64 of a kset's sorted values. `grm count` records it in
/// counted.json, so later stages can tell two ksets with the same number of
/// k-mers apart, also after the files are deleted.
pub fn digest(values: &[u64]) -> String {
    let mut h = Fnv::new();
    for &v in values {
        h.u64(v);
    }
    h.hex()
}

fn hex_digest(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

/// SHA-256 over sorted mixed canonical k-mers encoded as little-endian `u64`
/// bytes (64-char lowercase hex).
pub fn sha256_kmers(values: &[u64]) -> String {
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    for &v in values {
        ctx.update(&v.to_le_bytes());
    }
    hex_digest(ctx.finish().as_ref())
}

fn normalize_acgtn(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(seq.len());
    for &b in seq {
        if b.is_ascii_whitespace() {
            continue;
        }
        out.push(match b {
            b'A' | b'a' => b'A',
            b'C' | b'c' => b'C',
            b'G' | b'g' => b'G',
            b'T' | b't' => b'T',
            _ => b'N',
        });
    }
    out
}

fn reverse_complement_acgtn(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Strand-canonical uppercase ACGTN contig sequence: the lexicographically
/// smaller of the normalized forward sequence and its reverse complement.
pub fn canonicalize_contig(seq: &[u8]) -> Vec<u8> {
    let fwd = normalize_acgtn(seq);
    let rc = reverse_complement_acgtn(&fwd);
    if rc < fwd { rc } else { fwd }
}

/// SHA-256 of canonicalized uppercase ACGTN contigs sorted lexicographically,
/// length-delimited so contig boundaries and copy counts stay distinct.
fn sha256_normalized_contigs(contigs: &mut [Vec<u8>]) -> String {
    contigs.sort_unstable();
    let mut ctx = ring::digest::Context::new(&ring::digest::SHA256);
    ctx.update(b"GRM-ASM-V1\0");
    ctx.update(&(contigs.len() as u64).to_le_bytes());
    for contig in contigs.iter() {
        ctx.update(&(contig.len() as u64).to_le_bytes());
        ctx.update(contig);
    }
    hex_digest(ctx.finish().as_ref())
}

/// In-memory sequence fingerprint of one FASTA assembly (plain or gzip),
/// computed without writing a `.kset` file.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct FastaFingerprint {
    pub k_digest: usize,
    pub k_sketch: usize,
    pub sketch_size: usize,
    /// Deterministic 16-char lowercase hex matching [`digest`] (and
    /// `counted.json`'s `kset_digest` when `k_digest == 31`).
    pub fnv64_k31: String,
    /// 64-char lowercase hex SHA-256 over sorted mixed canonical k-mers at
    /// `k_digest` (little-endian `u64` bytes).
    pub sha256_k31: String,
    /// 64-char lowercase hex SHA-256 over strand-canonical uppercase ACGTN
    /// contig sequences sorted lexicographically.
    pub assembly_sha256: String,
    /// Distinct canonical k-mer count at `k_digest`.
    pub n_kmers: u64,
    pub contigs: u32,
    pub bp: u64,
    /// Bottom-`sketch_size` MinHash sketch at `k_sketch` (compatible with
    /// `species-check/sketches/*.json` at `k_sketch = 21`, `sketch_size = 1000`).
    pub sketch_k21: Vec<u64>,
    /// Optional bottom-`sketch_size` MinHash sketch at `k_digest`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sketch_k31: Option<Vec<u64>>,
}

/// Computes a [`FastaFingerprint`] directly in memory from a plain or `.gz`
/// FASTA file without writing a `.kset` file.
pub fn fingerprint_fasta(
    path: &Path,
    k_digest: usize,
    k_sketch: usize,
    sketch_size: usize,
) -> Result<FastaFingerprint> {
    fingerprint_fasta_into(
        path,
        k_digest,
        k_sketch,
        sketch_size,
        &mut Vec::new(),
        &mut Vec::new(),
    )
}

/// [`fingerprint_fasta`] reusing caller-owned sort buffers across genomes.
pub fn fingerprint_fasta_into(
    path: &Path,
    k_digest: usize,
    k_sketch: usize,
    sketch_size: usize,
    values: &mut Vec<u64>,
    scratch: &mut Vec<u64>,
) -> Result<FastaFingerprint> {
    ensure!(
        (MIN_K..=MAX_K).contains(&k_digest),
        "k_digest must be within {MIN_K}..={MAX_K}, got {k_digest}"
    );
    ensure!(
        (MIN_K..=MAX_K).contains(&k_sketch),
        "k_sketch must be within {MIN_K}..={MAX_K}, got {k_sketch}"
    );
    let space_digest = KmerSpace::new(k_digest);
    let space_sketch = KmerSpace::new(k_sketch);
    let separate_sketch = k_sketch != k_digest && sketch_size > 0;
    let sketch_cap = sketch_size.saturating_mul(4).max(16);

    values.clear();
    let mut kept_sketch: Vec<u64> = if separate_sketch {
        Vec::with_capacity(sketch_cap)
    } else {
        Vec::new()
    };
    let mut cut = u64::MAX;
    let mut full_sketch = false;
    let mut stats = StatsBuilder::default();
    let mut contigs: Vec<Vec<u8>> = Vec::new();

    fasta::for_each_record(path, |_, seq| {
        stats.add(seq);
        contigs.push(canonicalize_contig(seq));
        space_digest.for_each_canonical(seq, |x| values.push(space_digest.mix(x)));
        if separate_sketch {
            space_sketch.for_each_canonical(seq, |x| {
                let h = space_sketch.mix(x);
                if !full_sketch || h < cut {
                    kept_sketch.push(h);
                    if kept_sketch.len() >= sketch_cap {
                        kept_sketch.sort_unstable();
                        kept_sketch.dedup();
                        kept_sketch.truncate(sketch_size);
                        if kept_sketch.len() == sketch_size {
                            cut = kept_sketch[sketch_size - 1];
                            full_sketch = true;
                        }
                    }
                }
            });
        }
    })?;

    radix_sort(values, scratch);
    dedup_with_min_count(values, 1);

    let sketch_k21 = if sketch_size == 0 {
        Vec::new()
    } else if separate_sketch {
        kept_sketch.sort_unstable();
        kept_sketch.dedup();
        kept_sketch.truncate(sketch_size);
        kept_sketch
    } else {
        values.iter().copied().take(sketch_size).collect()
    };
    let sketch_k31 = (sketch_size > 0).then(|| values.iter().copied().take(sketch_size).collect());

    let stats = stats.finish();
    Ok(FastaFingerprint {
        k_digest,
        k_sketch,
        sketch_size,
        fnv64_k31: digest(values),
        sha256_k31: sha256_kmers(values),
        assembly_sha256: sha256_normalized_contigs(&mut contigs),
        n_kmers: values.len() as u64,
        contigs: stats.contigs,
        bp: stats.length,
        sketch_k21,
        sketch_k31,
    })
}

/// Random access to one kset's buckets.
pub struct KsetReader {
    file: File,
    pub k: usize,
    pub name: String,
    pub n_kmers: u64,
    offsets: Vec<u64>,
    counts: Vec<u32>,
    data_start: u64,
}

impl KsetReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        let mut head = [0u8; 24];
        file.read_exact(&mut head)
            .with_context(|| format!("read the header of {}", path.display()))?;
        ensure!(&head[..8] == MAGIC, "{} is not a kset file", path.display());
        let k = head[8] as usize;
        ensure!(
            head[9] as u32 == BUCKET_BITS,
            "unsupported bucket width {}",
            head[9]
        );
        ensure!(
            k <= MAX_K && 2 * k as u32 >= BUCKET_BITS,
            "{}: unsupported k {k}",
            path.display()
        );
        let n_kmers = u64::from_le_bytes(head[12..20].try_into()?);
        let name_len = u32::from_le_bytes(head[20..24].try_into()?) as usize;
        let mut table = vec![0u8; (BUCKETS + 1) * 8 + BUCKETS * 4];
        let data_start = (24 + name_len + table.len()) as u64;
        ensure!(data_start <= file_len, "{} is truncated", path.display());
        let mut name = vec![0u8; name_len];
        file.read_exact(&mut name)?;
        file.read_exact(&mut table)?;
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
        Ok(KsetReader {
            file,
            k,
            name: String::from_utf8(name)?,
            n_kmers,
            offsets,
            counts,
            data_start,
        })
    }

    pub fn bucket_len(&self, bucket: usize) -> usize {
        self.counts[bucket] as usize
    }

    /// Appends bucket `bucket`'s values to `out`.
    pub fn read_bucket(
        &self,
        bucket: usize,
        bytes: &mut Vec<u8>,
        out: &mut Vec<u64>,
    ) -> Result<()> {
        let (start, end) = (self.offsets[bucket], self.offsets[bucket + 1]);
        bytes.clear();
        bytes.resize((end - start) as usize, 0);
        self.file.read_exact_at(bytes, self.data_start + start)?;
        let bits = 2 * self.k as u32;
        let mut prev = bucket_floor(bucket, bits);
        let mut pos = 0;
        for i in 0..self.counts[bucket] {
            let gap = get_varint(bytes, &mut pos)?;
            // Values are distinct: only the first may equal the floor.
            ensure!(
                i == 0 || gap > 0,
                "kset {} bucket {bucket} repeats a k-mer",
                self.name
            );
            prev = prev
                .checked_add(gap)
                .with_context(|| format!("kset {} bucket {bucket} overflows", self.name))?;
            out.push(prev);
        }
        ensure!(
            pos == bytes.len() && (self.counts[bucket] == 0 || bucket_of(prev, bits) == bucket),
            "kset {} bucket {bucket} does not match its count and range",
            self.name
        );
        Ok(())
    }

    pub fn read_all(&self) -> Result<Vec<u64>> {
        let mut out = Vec::with_capacity(self.n_kmers as usize);
        let mut bytes = Vec::new();
        for b in 0..BUCKETS {
            self.read_bucket(b, &mut bytes, &mut out)?;
        }
        Ok(out)
    }

    /// The `s` smallest mixed values: a bottom-s MinHash sketch.
    pub fn sketch(&self, s: usize) -> Result<Vec<u64>> {
        let mut out = Vec::with_capacity(s);
        let mut bytes = Vec::new();
        for b in 0..BUCKETS {
            if out.len() >= s {
                break;
            }
            self.read_bucket(b, &mut bytes, &mut out)?;
        }
        out.truncate(s);
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn radix_sort_matches_std_sort() {
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        for n in [0, 10, 255, 256, 10_000] {
            let mut v: Vec<u64> = (0..n)
                .map(|_| rng.random::<u64>() >> rng.random_range(0..40))
                .collect();
            let mut want = v.clone();
            want.sort_unstable();
            radix_sort(&mut v, &mut Vec::new());
            assert_eq!(v, want, "n={n}");
        }
    }

    #[test]
    fn dedup_respects_min_count() {
        let mut v = vec![1, 1, 1, 2, 3, 3, 4, 4, 4, 4];
        dedup_with_min_count(&mut v, 3);
        assert_eq!(v, vec![1, 4]);
        let mut v = vec![1, 1, 2];
        dedup_with_min_count(&mut v, 1);
        assert_eq!(v, vec![1, 2]);
    }

    #[test]
    fn write_checked_preflights_exact_bytes_before_creating_a_file() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(31);
        let values = [0, 1, 1023, u64::MAX >> 2];
        let path = dir.path().join("checked.kset");
        let seen = std::cell::Cell::new(None);
        write_checked(&path, "g", &space, &values, |_, need| {
            assert!(!path.with_extension("kset.part").exists());
            seen.set(Some(need));
            Ok(())
        })
        .unwrap();
        assert_eq!(std::fs::metadata(&path).unwrap().len(), seen.get().unwrap());

        let refused = dir.path().join("refused.kset");
        let err = write_checked(&refused, "g", &space, &values, |_, _| -> Result<()> {
            anyhow::bail!("injected full disk")
        })
        .unwrap_err();
        assert!(format!("{err:#}").contains("injected full disk"));
        assert!(!refused.exists());
        assert!(!refused.with_extension("kset.part").exists());
    }

    #[test]
    fn the_size_bound_holds_for_the_widest_gaps_at_every_k() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(3);
        for k in [crate::dna::MIN_K, 21, 31, MAX_K] {
            let space = KmerSpace::new(k);
            let mask = if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            };
            // One value per bucket, at the top of its range: the longest gaps.
            let bits = space.bits();
            let mut values: Vec<u64> = (0..BUCKETS)
                .map(|b| bucket_floor(b, bits) | (mask >> BUCKET_BITS))
                .collect();
            values.extend((0..500).map(|_| rng.random::<u64>() & mask));
            values.sort_unstable();
            values.dedup();
            let path = dir.path().join(format!("k{k}.kset"));
            write(&path, "genome", &space, &values).unwrap();
            let len = std::fs::metadata(&path).unwrap().len();
            assert!(
                len <= max_encoded_len("genome", values.len() as u64),
                "k={k}: {len} bytes"
            );
        }
    }

    #[test]
    fn writes_refuse_links_at_the_kset_and_its_temporary() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(21);
        let checked = |path: &Path| write_checked(path, "g", &space, &[1, 2], |_, _| Ok(()));
        let victim = dir.path().join("victim");
        std::fs::write(&victim, "shared store").unwrap();
        let linked = dir.path().join("linked.kset");
        std::os::unix::fs::symlink(&victim, &linked).unwrap();
        assert!(checked(&linked).is_err());
        let part_linked = dir.path().join("part.kset");
        std::os::unix::fs::symlink(&victim, part_linked.with_extension("kset.part")).unwrap();
        assert!(checked(&part_linked).is_err());
        assert_eq!(std::fs::read(&victim).unwrap(), b"shared store");
        assert!(linked.is_symlink() && !part_linked.exists());
        // A stale regular temporary from an interrupted write is replaced.
        let stale = dir.path().join("stale.kset");
        std::fs::write(stale.with_extension("kset.part"), "torn").unwrap();
        checked(&stale).unwrap();
        assert_eq!(
            KsetReader::open(&stale).unwrap().read_all().unwrap(),
            [1, 2]
        );
        assert!(!stale.with_extension("kset.part").exists());
    }

    #[test]
    fn kset_round_trips_and_sketches() {
        let dir = tempfile::tempdir().unwrap();
        let space = KmerSpace::new(31);
        let mut rng = ChaCha8Rng::seed_from_u64(5);
        let seq: Vec<u8> = (0..50_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let fna = dir.path().join("g.fna");
        std::fs::write(&fna, [b">c\n".as_slice(), &seq, b"\n"].concat()).unwrap();

        let (values, stats) = count_file(&fna, &space, 1).unwrap();
        assert_eq!(stats.length, 50_000);
        assert!(values.len() > 49_000 && values.len() <= 50_000 - 30);
        assert!(values.windows(2).all(|w| w[0] < w[1]));

        let path = dir.path().join("g.kset");
        write(&path, "g1", &space, &values).unwrap();
        let reader = KsetReader::open(&path).unwrap();
        assert_eq!(reader.name, "g1");
        assert_eq!(reader.k, 31);
        assert_eq!(reader.read_all().unwrap(), values);
        assert_eq!(reader.sketch(100).unwrap(), values[..100].to_vec());
        let total: usize = (0..BUCKETS).map(|b| reader.bucket_len(b)).sum();
        assert_eq!(total, values.len());

        // Every stored value decodes to a k-mer that occurs in the sequence.
        let text = String::from_utf8(seq).unwrap();
        for &v in values.iter().take(50) {
            let kmer = space.decode(space.unmix(v));
            let rc = space.decode(space.reverse_complement(space.unmix(v)));
            assert!(text.contains(&kmer) || text.contains(&rc));
        }
    }

    #[test]
    fn kset_round_trips_bucket_edges_at_both_k_limits() {
        // k = 5 leaves no bits below the bucket prefix; k = 32 uses all 64.
        let dir = tempfile::tempdir().unwrap();
        for k in [5, 32] {
            let space = KmerSpace::new(k);
            let low = space.bits() - BUCKET_BITS;
            let mut values = vec![0, 1, u64::MAX >> (64 - space.bits())];
            for b in [1u64, 511, 512, BUCKETS as u64 - 1] {
                values.extend([(b << low) - 1, b << low]);
            }
            values.sort_unstable();
            values.dedup();
            let path = dir.path().join(format!("k{k}.kset"));
            write(&path, "g", &space, &values).unwrap();
            let reader = KsetReader::open(&path).unwrap();
            assert_eq!(reader.k, k);
            assert_eq!(reader.read_all().unwrap(), values, "k={k}");
        }
    }

    const COUNTS: usize = 25 + (BUCKETS + 1) * 8;
    const DATA: usize = COUNTS + BUCKETS * 4;

    fn kset_bytes(dir: &Path, k: usize, values: &[u64]) -> Vec<u8> {
        let path = dir.join("valid.kset");
        write(&path, "g", &KmerSpace::new(k), values).unwrap();
        std::fs::read(path).unwrap()
    }

    fn reopen(dir: &Path, bytes: &[u8]) -> Result<KsetReader> {
        let path = dir.join("case.kset");
        std::fs::write(&path, bytes).unwrap();
        KsetReader::open(&path)
    }

    fn patched(bytes: &[u8], at: usize, value: &[u8]) -> Vec<u8> {
        let mut out = bytes.to_vec();
        out[at..at + value.len()].copy_from_slice(value);
        out
    }

    fn random_values(n: usize) -> Vec<u64> {
        let mut rng = ChaCha8Rng::seed_from_u64(17);
        let mut values: Vec<u64> = (0..n).map(|_| rng.random::<u64>() >> 2).collect();
        values.sort_unstable();
        values.dedup();
        values
    }

    #[test]
    fn open_rejects_unsupported_k() {
        let dir = tempfile::tempdir().unwrap();
        let mut bytes = kset_bytes(dir.path(), 31, &random_values(100));
        assert!(reopen(dir.path(), &bytes).is_ok());
        for k in [0, 4, 33, 255] {
            bytes[8] = k;
            assert!(reopen(dir.path(), &bytes).is_err(), "k={k}");
        }
    }

    #[test]
    fn open_rejects_truncated_or_padded_files() {
        let dir = tempfile::tempdir().unwrap();
        let bytes = kset_bytes(dir.path(), 31, &random_values(1_000));
        for cut in [0, 24, DATA - 1, DATA, bytes.len() - 1] {
            assert!(reopen(dir.path(), &bytes[..cut]).is_err(), "cut={cut}");
        }
        assert!(reopen(dir.path(), &[bytes.as_slice(), &[0]].concat()).is_err());
        let huge_name = patched(&bytes, 20, &u32::MAX.to_le_bytes());
        assert!(reopen(dir.path(), &huge_name).is_err());
    }

    #[test]
    fn a_file_shorter_than_the_header_is_named_in_the_error() {
        let dir = tempfile::tempdir().unwrap();
        let bytes = kset_bytes(dir.path(), 31, &random_values(10));
        for cut in [0, 23] {
            let err = reopen(dir.path(), &bytes[..cut]).err().unwrap();
            assert!(
                format!("{err:#}").contains("case.kset"),
                "cut={cut}: {err:#}"
            );
        }
    }

    #[test]
    fn open_rejects_inconsistent_bucket_table() {
        let dir = tempfile::tempdir().unwrap();
        let bytes = kset_bytes(dir.path(), 31, &random_values(20_000));
        let n_kmers = u64::from_le_bytes(bytes[12..20].try_into().unwrap());
        let second_end = u64::from_le_bytes(bytes[41..49].try_into().unwrap());
        let first_count = u32::from_le_bytes(bytes[COUNTS..COUNTS + 4].try_into().unwrap());
        let oversized = patched(&bytes, COUNTS, &u32::MAX.to_le_bytes());
        let rebalanced = n_kmers - first_count as u64 + u32::MAX as u64;
        let cases = [
            patched(&bytes, 12, &(n_kmers + 1).to_le_bytes()),
            patched(&bytes, 25, &1u64.to_le_bytes()),
            patched(&bytes, 33, &(second_end + 1).to_le_bytes()),
            patched(&oversized, 12, &rebalanced.to_le_bytes()),
        ];
        assert!(reopen(dir.path(), &bytes).is_ok());
        for (i, case) in cases.iter().enumerate() {
            assert!(reopen(dir.path(), case).is_err(), "case {i}");
        }
    }

    #[test]
    fn read_bucket_rejects_corrupt_payload() {
        let dir = tempfile::tempdir().unwrap();
        let (mut bytes, mut out) = (Vec::new(), Vec::new());

        let narrow = kset_bytes(dir.path(), 5, &[0, 1023]);
        let reader = reopen(dir.path(), &patched(&narrow, DATA, &[1])).unwrap();
        assert!(reader.read_bucket(0, &mut bytes, &mut out).is_err());

        let wide = kset_bytes(dir.path(), 32, &[u64::MAX]);
        let reader = reopen(dir.path(), &patched(&wide, wide.len() - 1, &[0x7f])).unwrap();
        assert!(
            reader
                .read_bucket(BUCKETS - 1, &mut bytes, &mut out)
                .is_err()
        );

        let random = kset_bytes(dir.path(), 31, &random_values(20_000));
        let first_end = u64::from_le_bytes(random[33..41].try_into().unwrap());
        let shifted = patched(&random, 33, &(first_end - 1).to_le_bytes());
        let reader = reopen(dir.path(), &shifted).unwrap();
        assert!(reader.read_bucket(1, &mut bytes, &mut out).is_err());
    }

    #[test]
    fn read_bucket_rejects_a_repeated_kmer() {
        // A gap of 0 is the bucket floor itself for the first value, and a
        // repeat after it: a merge would count the genome twice.
        let dir = tempfile::tempdir().unwrap();
        let (mut bytes, mut out) = (Vec::new(), Vec::new());
        for first in [0u64, 7] {
            let valid = kset_bytes(dir.path(), 31, &[first, first + 5]);
            let reader = reopen(dir.path(), &valid).unwrap();
            reader.read_bucket(0, &mut bytes, &mut out).unwrap();
            assert_eq!(out, [first, first + 5]);
            let repeated = patched(&valid, valid.len() - 1, &[0]);
            let reader = reopen(dir.path(), &repeated).unwrap();
            let err = reader.read_bucket(0, &mut bytes, &mut out).err().unwrap();
            assert!(format!("{err:#}").contains("repeats a k-mer"), "{err:#}");
            assert!(reader.read_all().is_err());
            out.clear();
        }
    }

    #[test]
    fn fingerprint_fasta_matches_kset_digest_and_species_sketch_for_plain_and_gz() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let c1: Vec<u8> = (0..12_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let c2: Vec<u8> = (0..8_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let fasta_text = [b">c1 desc\n".as_slice(), &c1, b"\n>c2\n", &c2, b"\n"].concat();
        let plain = dir.path().join("g.fna");
        std::fs::write(&plain, &fasta_text).unwrap();
        let gz = dir.path().join("g.fna.gz");
        let mut enc =
            flate2::write::GzEncoder::new(File::create(&gz).unwrap(), flate2::Compression::fast());
        enc.write_all(&fasta_text).unwrap();
        enc.finish().unwrap();

        let fp_plain = fingerprint_fasta(&plain, 31, 21, 1000).unwrap();
        let fp_gz = fingerprint_fasta(&gz, 31, 21, 1000).unwrap();
        assert_eq!(fp_plain, fp_gz);

        let (values31, stats31) = count_file(&plain, &KmerSpace::new(31), 1).unwrap();
        assert_eq!(fp_plain.fnv64_k31, digest(&values31));
        assert_eq!(fp_plain.sha256_k31, sha256_kmers(&values31));
        assert_eq!(fp_plain.sha256_k31.len(), 64);
        assert_eq!(fp_plain.fnv64_k31.len(), 16);
        assert_eq!(fp_plain.n_kmers, values31.len() as u64);
        assert_eq!(fp_plain.contigs, stats31.contigs);
        assert_eq!(fp_plain.bp, stats31.length);
        assert_eq!(
            fp_plain.sketch_k31.as_deref(),
            Some(&values31[..1000.min(values31.len())])
        );

        let sp_sketch =
            crate::species_check::sketch_fasta("g", &plain, &KmerSpace::new(21), 1000).unwrap();
        assert_eq!(fp_plain.sketch_k21, sp_sketch.sketch);
    }

    #[test]
    fn assembly_sha256_normalizes_layout_and_separates_from_kset_identity() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(99);
        let c1: String = (0..500)
            .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
            .collect();
        let c2: String = (0..300)
            .map(|_| b"ACGT"[rng.random_range(0..4)] as char)
            .collect();
        let rc_c2: String = reverse_complement_acgtn(c2.as_bytes())
            .into_iter()
            .map(char::from)
            .collect();

        let p1 = dir.path().join("a.fna");
        std::fs::write(&p1, format!(">h1\n{c1}\n>h2\n{c2}\n")).unwrap();
        // Reordered, different headers, lowercase, wrapped lines, reverse-complemented contig.
        let p2 = dir.path().join("b.fna");
        let wrapped_c1 = c1
            .to_lowercase()
            .as_bytes()
            .chunks(60)
            .map(|chunk| std::str::from_utf8(chunk).unwrap())
            .collect::<Vec<_>>()
            .join("\n");
        std::fs::write(&p2, format!(">other_a\n{rc_c2}\n>other_b\n{wrapped_c1}\n")).unwrap();
        // Duplicated contig: identical canonical k-mer set, different assembly content!
        let p3 = dir.path().join("c.fna");
        std::fs::write(&p3, format!(">h1\n{c1}\n>h2\n{c2}\n>h2_dup\n{c2}\n")).unwrap();

        let f1 = fingerprint_fasta(&p1, 31, 21, 200).unwrap();
        let f2 = fingerprint_fasta(&p2, 31, 21, 200).unwrap();
        let f3 = fingerprint_fasta(&p3, 31, 21, 200).unwrap();

        assert_eq!(f1.assembly_sha256, f2.assembly_sha256);
        assert_eq!(f1.fnv64_k31, f2.fnv64_k31);
        assert_eq!(f1.sha256_k31, f2.sha256_k31);
        assert_eq!(f1.sketch_k21, f2.sketch_k21);

        // Duplicated contig preserves canonical k-mer set but changes assembly content hash.
        assert_eq!(f1.fnv64_k31, f3.fnv64_k31);
        assert_eq!(f1.sha256_k31, f3.sha256_k31);
        assert_ne!(f1.assembly_sha256, f3.assembly_sha256);
        assert_eq!(f3.contigs, 3);
    }
}
