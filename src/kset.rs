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

use crate::dna::{KmerSpace, MAX_K};
use crate::fasta::{self, SeqStats, StatsBuilder};
use crate::fnv::Fnv;
use anyhow::{Context, Result, bail, ensure};
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

    let part = path.with_extension("kset.part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(MAGIC)?;
    w.write_all(&[space.k as u8, BUCKET_BITS as u8, 0, 0])?;
    w.write_all(&(values.len() as u64).to_le_bytes())?;
    w.write_all(&(name.len() as u32).to_le_bytes())?;
    w.write_all(name.as_bytes())?;
    for o in &offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for c in &counts {
        w.write_all(&c.to_le_bytes())?;
    }
    w.write_all(&data)?;
    // No fsync: a kset is a re-creatable cache, and the rename below
    // already keeps an interrupted run from leaving a truncated file.
    w.into_inner()?;
    std::fs::rename(&part, path)?;
    Ok(())
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
        for _ in 0..self.counts[bucket] {
            prev = prev
                .checked_add(get_varint(bytes, &mut pos)?)
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
}
