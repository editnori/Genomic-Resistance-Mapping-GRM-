//! Binary k-mer storage formats for `grm kcount`:
//!
//! - Counted table (`.ktab`, magic `GRMKTB1\0`): 1,024-bucket Golomb-Rice
//!   sorted mixed keys paired with exact uncapped `u32` counts (0 extra bits
//!   when presence-only, tiered 1B/3B/5B or `u16` + `0xFFFF` escape to `u32`
//!   otherwise), supporting `k = 5..=64`.
//! - Golomb-Rice sorted-key table (`.grk`, magic `GRMKGR1\0`): presence-only
//!   1,024-bucket Golomb-Rice encoding of sorted mixed keys, used as the
//!   random-access sorted-integer baseline alongside the SPSS archive format
//!   (`GRMSPS1\0` in [`crate::kcount::spss`]).
//! - FastK-style per-sequence k-mer count profile (`.prof.bin`, magic
//!   `GRMPRF1\0`): run-length and zigzag-delta encoded positional k-mer counts
//!   along each contig or read.
//! - Transparent read interop with untouched `.kset` v1 (`GRMKSET1`) via
//!   [`AnyTableReader`].

use crate::dna::KmerSpace;
use crate::kcount::key::{KmerKey, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{self, BUCKET_BITS, BUCKETS, KsetReader};
use anyhow::{Context, Result, bail, ensure};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

pub const MAGIC_GRKSET: &[u8; 8] = b"GRMKGR1\0";
pub const MAGIC_KTAB: &[u8; 8] = b"GRMKTB1\0";
pub const MAGIC_PROF1: &[u8; 8] = b"GRMPRF1\0";

const FLAG_PRESENCE_ONLY: u8 = 0x01;
const COUNT_MODE_ALL_ONE: u8 = 0;
const COUNT_MODE_TIERED: u8 = 1;
const COUNT_MODE_U16_ESC: u8 = 2;

/// Identifies which random-access binary k-mer table format a file holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TableFormat {
    KsetV1,
    GrKset,
    Ktab,
}

/// Inspects the 8-byte file magic to detect a random-access k-mer table format.
pub fn detect_format(path: &Path) -> Result<TableFormat> {
    let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let mut magic = [0u8; 8];
    file.read_exact(&mut magic)
        .with_context(|| format!("read magic of {}", path.display()))?;
    match &magic {
        b"GRMKSET1" => Ok(TableFormat::KsetV1),
        b"GRMKGR1\0" => Ok(TableFormat::GrKset),
        b"GRMKTB1\0" => Ok(TableFormat::Ktab),
        _ => bail!("{}: unrecognized k-mer table magic", path.display()),
    }
}

/// Fast LSB-first bit writer using a 128-bit accumulator.
struct BitWriter<'a> {
    out: &'a mut Vec<u8>,
    acc: u128,
    bits: u32,
}

impl<'a> BitWriter<'a> {
    fn new(out: &'a mut Vec<u8>) -> Self {
        Self {
            out,
            acc: 0,
            bits: 0,
        }
    }

    #[inline]
    fn write_u64(&mut self, v: u64, n: u32) {
        if n == 0 {
            return;
        }
        debug_assert!(n <= 64);
        let masked = if n == 64 { v } else { v & ((1u64 << n) - 1) };
        self.acc |= (masked as u128) << self.bits;
        self.bits += n;
        if self.bits >= 64 {
            self.out.extend_from_slice(&(self.acc as u64).to_le_bytes());
            self.acc >>= 64;
            self.bits -= 64;
        }
    }

    #[inline]
    fn write_u128(&mut self, v: u128, n: u32) {
        if n <= 64 {
            self.write_u64(v as u64, n);
        } else {
            self.write_u64(v as u64, 64);
            self.write_u64((v >> 64) as u64, n - 64);
        }
    }

    #[inline]
    fn write_unary(&mut self, mut q: u64) {
        while q >= 64 {
            self.write_u64(u64::MAX, 64);
            q -= 64;
        }
        if q > 0 {
            self.write_u64((1u64 << q) - 1, q as u32);
        }
        self.write_u64(0, 1);
    }

    fn finish(mut self) {
        while self.bits > 0 {
            self.out.push(self.acc as u8);
            self.acc >>= 8;
            self.bits = self.bits.saturating_sub(8);
        }
    }
}

/// Fast LSB-first bit reader using a 128-bit accumulator.
struct BitReader<'a> {
    data: &'a [u8],
    byte_pos: usize,
    acc: u128,
    bits: u32,
}

impl<'a> BitReader<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            byte_pos: 0,
            acc: 0,
            bits: 0,
        }
    }

    #[inline]
    fn refill(&mut self) {
        while self.bits <= 64 && self.byte_pos < self.data.len() {
            if self.byte_pos + 8 <= self.data.len() {
                let w = u64::from_le_bytes(
                    self.data[self.byte_pos..self.byte_pos + 8]
                        .try_into()
                        .unwrap(),
                );
                self.acc |= (w as u128) << self.bits;
                self.bits += 64;
                self.byte_pos += 8;
            } else {
                self.acc |= (self.data[self.byte_pos] as u128) << self.bits;
                self.bits += 8;
                self.byte_pos += 1;
            }
        }
    }

    #[inline]
    fn read_u64(&mut self, n: u32) -> Result<u64> {
        if n == 0 {
            return Ok(0);
        }
        if self.bits < n {
            self.refill();
            ensure!(self.bits >= n, "truncated Golomb-Rice bitstream");
        }
        let mask = if n == 64 { u64::MAX } else { (1u64 << n) - 1 };
        let v = (self.acc as u64) & mask;
        self.acc >>= n;
        self.bits -= n;
        Ok(v)
    }

    #[inline]
    fn read_u128(&mut self, n: u32) -> Result<u128> {
        if n <= 64 {
            Ok(self.read_u64(n)? as u128)
        } else {
            let lo = self.read_u64(64)? as u128;
            let hi = self.read_u64(n - 64)? as u128;
            Ok(lo | (hi << 64))
        }
    }

    #[inline]
    fn read_unary(&mut self) -> Result<u64> {
        let mut q = 0u64;
        loop {
            if self.bits == 0 {
                self.refill();
                ensure!(self.bits > 0, "truncated unary quotient");
            }
            let ones = self.acc.trailing_ones().min(self.bits);
            q = q.checked_add(ones as u64).context("unary overflow")?;
            if ones < self.bits {
                let consume = ones + 1;
                if consume == 128 {
                    self.acc = 0;
                } else {
                    self.acc >>= consume;
                }
                self.bits -= consume;
                return Ok(q);
            }
            self.acc = 0;
            self.bits = 0;
        }
    }

    fn verify_exhausted(&self) -> Result<()> {
        ensure!(
            self.byte_pos == self.data.len() && self.bits < 8 && self.acc == 0,
            "trailing non-zero bits or unconsumed bytes in bucket payload"
        );
        Ok(())
    }
}

/// Encodes one bucket's strictly sorted mixed keys using bit-optimal
/// per-bucket Golomb-Rice gap coding.
pub fn encode_bucket_keys<K: KmerKey>(k: usize, bucket: usize, keys: &[K], out: &mut Vec<u8>) {
    if keys.is_empty() {
        return;
    }
    debug_assert!(keys.windows(2).all(|w| w[0] < w[1]));
    let low_bits = 2 * k as u32 - BUCKET_BITS;
    let floor = K::bucket_floor(bucket, k, BUCKET_BITS).to_u128();
    let n = keys.len();
    let last = keys[n - 1].to_u128();
    debug_assert!(keys[0].to_u128() >= floor);
    let span = (last - floor).saturating_sub((n - 1) as u128);

    let mean_ln2 = (span as f64 / n as f64) * std::f64::consts::LN_2;
    let base_r = if mean_ln2 < 1.0 {
        0u32
    } else {
        (mean_ln2.log2().round() as u32).min(low_bits)
    };
    let r_lo = base_r.saturating_sub(1);
    let r_hi = (base_r + 1).min(low_bits);
    let mut best_r = base_r;
    let mut best_bits = u128::MAX;
    for r in r_lo..=r_hi {
        let mut q_sum = 0u128;
        let mut prev = floor;
        for (i, &key) in keys.iter().enumerate() {
            let v = key.to_u128();
            let gap = if i == 0 { v - prev } else { v - prev - 1 };
            prev = v;
            q_sum = q_sum.saturating_add(gap >> r);
        }
        let cost = (n as u128) * ((r + 1) as u128) + q_sum;
        if cost < best_bits {
            best_bits = cost;
            best_r = r;
        }
    }

    out.push(best_r as u8);
    let rem_mask = if best_r == 128 {
        u128::MAX
    } else {
        (1u128 << best_r) - 1
    };
    let mut bw = BitWriter::new(out);
    let mut prev = floor;
    for (i, &key) in keys.iter().enumerate() {
        let v = key.to_u128();
        let gap = if i == 0 { v - prev } else { v - prev - 1 };
        prev = v;
        let q = (gap >> best_r) as u64;
        bw.write_unary(q);
        bw.write_u128(gap & rem_mask, best_r);
    }
    bw.finish();
}

/// Decodes one bucket's Golomb-Rice payload into `out`, validating strict
/// monotonicity, bucket range, and payload boundaries.
pub fn decode_bucket_keys<K: KmerKey>(
    k: usize,
    bucket: usize,
    count: u32,
    bytes: &[u8],
    out: &mut Vec<K>,
) -> Result<()> {
    if count == 0 {
        ensure!(
            bytes.is_empty(),
            "empty bucket {bucket} carries trailing data"
        );
        return Ok(());
    }
    ensure!(
        !bytes.is_empty(),
        "non-empty bucket {bucket} has empty payload"
    );
    let low_bits = 2 * k as u32 - BUCKET_BITS;
    let r = bytes[0] as u32;
    ensure!(
        r <= low_bits,
        "bucket {bucket} Golomb-Rice parameter {r} exceeds low width {low_bits}"
    );
    let floor = K::bucket_floor(bucket, k, BUCKET_BITS).to_u128();
    let max_val = if low_bits == 128 {
        u128::MAX
    } else {
        floor | ((1u128 << low_bits) - 1)
    };
    let mut br = BitReader::new(&bytes[1..]);
    let mut prev = floor;
    out.reserve(count as usize);
    for i in 0..count {
        let q = br.read_unary()? as u128;
        let rem = br.read_u128(r)?;
        let shifted_q = if r == 128 {
            ensure!(q == 0, "bucket {bucket} quotient overflow");
            0
        } else {
            q.checked_shl(r)
                .filter(|&v| (v >> r) == q)
                .with_context(|| format!("bucket {bucket} gap shift overflow"))?
        };
        let gap = shifted_q | rem;
        let v = if i == 0 {
            prev.checked_add(gap)
        } else {
            prev.checked_add(1).and_then(|p| p.checked_add(gap))
        }
        .with_context(|| format!("bucket {bucket} value overflow"))?;
        ensure!(
            v <= max_val,
            "bucket {bucket} decoded key exceeds bucket upper bound"
        );
        prev = v;
        out.push(K::from_u128(v));
    }
    br.verify_exhausted()?;
    Ok(())
}

/// Encodes one bucket's exact uncapped `u32` counts with automatic per-bucket
/// mode selection (0 extra bytes for all-ones, tiered 1B/3B/5B for typical
/// counts, or 2B `u16` + `0xFFFF` escape to `u32` for deep sequencing).
pub fn encode_bucket_counts(counts: &[u32], out: &mut Vec<u8>) {
    if counts.is_empty() {
        return;
    }
    if counts.iter().all(|&c| c == 1) {
        out.push(COUNT_MODE_ALL_ONE);
        return;
    }
    let mut tiered_bytes = 0usize;
    let mut u16_bytes = 0usize;
    for &c in counts {
        tiered_bytes += if c < 254 {
            1
        } else if c < 65_535 {
            3
        } else {
            5
        };
        u16_bytes += if c < 65_535 { 2 } else { 6 };
    }
    if tiered_bytes <= u16_bytes {
        out.push(COUNT_MODE_TIERED);
        out.reserve(tiered_bytes);
        for &c in counts {
            if c < 254 {
                out.push(c as u8);
            } else if c < 65_535 {
                out.push(254);
                out.extend_from_slice(&(c as u16).to_le_bytes());
            } else {
                out.push(255);
                out.extend_from_slice(&c.to_le_bytes());
            }
        }
    } else {
        out.push(COUNT_MODE_U16_ESC);
        out.reserve(u16_bytes);
        for &c in counts {
            if c < 65_535 {
                out.extend_from_slice(&(c as u16).to_le_bytes());
            } else {
                out.extend_from_slice(&0xFFFFu16.to_le_bytes());
                out.extend_from_slice(&c.to_le_bytes());
            }
        }
    }
}

/// Decodes one bucket's exact `u32` counts into `out`.
pub fn decode_bucket_counts(
    bucket: usize,
    n: u32,
    presence_only: bool,
    bytes: &[u8],
    out: &mut Vec<u32>,
) -> Result<()> {
    if n == 0 {
        ensure!(
            bytes.is_empty(),
            "empty bucket {bucket} carries count bytes"
        );
        return Ok(());
    }
    out.reserve(n as usize);
    if presence_only {
        ensure!(
            bytes.is_empty(),
            "presence-only bucket {bucket} carries count bytes"
        );
        out.resize(out.len() + n as usize, 1);
        return Ok(());
    }
    ensure!(
        !bytes.is_empty(),
        "bucket {bucket} missing count mode header"
    );
    let mode = bytes[0];
    let mut pos = 1usize;
    match mode {
        COUNT_MODE_ALL_ONE => {
            ensure!(
                bytes.len() == 1,
                "all-ones count bucket {bucket} has trailing bytes"
            );
            out.resize(out.len() + n as usize, 1);
        }
        COUNT_MODE_TIERED => {
            for _ in 0..n {
                let Some(&tag) = bytes.get(pos) else {
                    bail!("bucket {bucket} truncated tiered count");
                };
                pos += 1;
                let c = if tag < 254 {
                    tag as u32
                } else if tag == 254 {
                    ensure!(pos + 2 <= bytes.len(), "bucket {bucket} truncated u16");
                    let v = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap()) as u32;
                    pos += 2;
                    v
                } else {
                    ensure!(pos + 4 <= bytes.len(), "bucket {bucket} truncated u32");
                    let v = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap());
                    pos += 4;
                    v
                };
                out.push(c);
            }
            ensure!(
                pos == bytes.len(),
                "bucket {bucket} trailing bytes in tiered count payload"
            );
        }
        COUNT_MODE_U16_ESC => {
            for _ in 0..n {
                ensure!(pos + 2 <= bytes.len(), "bucket {bucket} truncated u16");
                let v = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap());
                pos += 2;
                let c = if v < 0xFFFF {
                    v as u32
                } else {
                    ensure!(pos + 4 <= bytes.len(), "bucket {bucket} truncated u32 esc");
                    let full = u32::from_le_bytes(bytes[pos..pos + 4].try_into().unwrap());
                    pos += 4;
                    full
                };
                out.push(c);
            }
            ensure!(
                pos == bytes.len(),
                "bucket {bucket} trailing bytes in u16-escape count payload"
            );
        }
        other => bail!("bucket {bucket} unknown count encoding mode {other}"),
    }
    Ok(())
}

/// Writes pre-encoded 1,024 bucket payloads as a Golomb-Rice sorted-key table
/// (`GRMKGR1\0`) atomically via `<path>.part` and rename.
pub fn write_gr_kset_from_buckets(
    path: &Path,
    name: &str,
    k: usize,
    bucket_counts: &[u32; BUCKETS],
    bucket_payloads: &[Vec<u8>],
) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    ensure!(
        bucket_payloads.len() == BUCKETS,
        "expected {BUCKETS} bucket payloads"
    );
    let n_kmers: u64 = bucket_counts.iter().map(|&c| c as u64).sum();
    let key_width: u8 = if k <= MAX_K_U64 { 8 } else { 16 };
    let mut offsets = [0u64; BUCKETS + 1];
    let mut cur = 0u64;
    for (b, payload) in bucket_payloads.iter().enumerate() {
        offsets[b] = cur;
        cur += payload.len() as u64;
    }
    offsets[BUCKETS] = cur;

    let part = crate::confine::suffixed(path, ".part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(MAGIC_GRKSET)?;
    w.write_all(&[k as u8, BUCKET_BITS as u8, key_width, 0])?;
    w.write_all(&0u32.to_le_bytes())?;
    w.write_all(&n_kmers.to_le_bytes())?;
    w.write_all(&(name.len() as u32).to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?;
    w.write_all(name.as_bytes())?;
    for o in offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for &c in bucket_counts {
        w.write_all(&c.to_le_bytes())?;
    }
    for payload in bucket_payloads {
        w.write_all(payload)?;
    }
    w.flush()?;
    drop(w);
    std::fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))?;
    Ok(())
}

/// Writes sorted distinct mixed keys as a Golomb-Rice sorted-key table
/// (`GRMKGR1\0`).
pub fn write_gr_kset<K: KmerKey>(path: &Path, name: &str, k: usize, values: &[K]) -> Result<()> {
    let mut bucket_counts = [0u32; BUCKETS];
    let mut bucket_payloads = vec![Vec::new(); BUCKETS];
    let mut i = 0usize;
    for b in 0..BUCKETS {
        let start = i;
        while i < values.len() && values[i].bucket_of(k, BUCKET_BITS) == b {
            i += 1;
        }
        let slice = &values[start..i];
        bucket_counts[b] = slice.len() as u32;
        encode_bucket_keys(k, b, slice, &mut bucket_payloads[b]);
    }
    ensure!(i == values.len(), "input keys were not sorted by bucket");
    write_gr_kset_from_buckets(path, name, k, &bucket_counts, &bucket_payloads)
}

/// Writes pre-encoded 1,024 bucket key and count payloads as a `.ktab`
/// (`GRMKTB1\0`) file atomically via `<path>.part` and rename.
#[allow(clippy::too_many_arguments)]
pub fn write_ktab_from_buckets(
    path: &Path,
    name: &str,
    k: usize,
    min_count: u32,
    max_count_observed: u32,
    presence_only: bool,
    bucket_counts: &[u32; BUCKETS],
    key_payloads: &[Vec<u8>],
    count_payloads: &[Vec<u8>],
) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    ensure!(
        key_payloads.len() == BUCKETS && count_payloads.len() == BUCKETS,
        "expected {BUCKETS} bucket payloads"
    );
    let n_kmers: u64 = bucket_counts.iter().map(|&c| c as u64).sum();
    let key_width: u8 = if k <= MAX_K_U64 { 8 } else { 16 };
    let flags: u8 = if presence_only { FLAG_PRESENCE_ONLY } else { 0 };

    let mut key_offsets = [0u64; BUCKETS + 1];
    let mut cur = 0u64;
    for (b, payload) in key_payloads.iter().enumerate() {
        key_offsets[b] = cur;
        cur += payload.len() as u64;
    }
    key_offsets[BUCKETS] = cur;

    let mut count_offsets = [0u64; BUCKETS + 1];
    let mut c_cur = 0u64;
    if !presence_only {
        for (b, payload) in count_payloads.iter().enumerate() {
            count_offsets[b] = c_cur;
            c_cur += payload.len() as u64;
        }
        count_offsets[BUCKETS] = c_cur;
    }

    let part = crate::confine::suffixed(path, ".part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(MAGIC_KTAB)?;
    w.write_all(&[k as u8, BUCKET_BITS as u8, key_width, flags])?;
    w.write_all(&min_count.to_le_bytes())?;
    w.write_all(&max_count_observed.to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?;
    w.write_all(&n_kmers.to_le_bytes())?;
    w.write_all(&(name.len() as u32).to_le_bytes())?;
    w.write_all(&0u32.to_le_bytes())?;
    w.write_all(name.as_bytes())?;
    for o in key_offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for o in count_offsets {
        w.write_all(&o.to_le_bytes())?;
    }
    for &c in bucket_counts {
        w.write_all(&c.to_le_bytes())?;
    }
    for payload in key_payloads {
        w.write_all(payload)?;
    }
    if !presence_only {
        for payload in count_payloads {
            w.write_all(payload)?;
        }
    }
    w.flush()?;
    drop(w);
    std::fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))?;
    Ok(())
}

/// Writes sorted `(mixed_key, count)` pairs as `.ktab` (`GRMKTB1\0`).
pub fn write_ktab<K: KmerKey>(
    path: &Path,
    name: &str,
    k: usize,
    min_count: u32,
    entries: &[(K, u32)],
) -> Result<()> {
    let mut bucket_counts = [0u32; BUCKETS];
    let mut key_payloads = vec![Vec::new(); BUCKETS];
    let mut count_payloads = vec![Vec::new(); BUCKETS];
    let max_count_observed = entries.iter().map(|&(_, c)| c).max().unwrap_or(0);
    let presence_only = entries.iter().all(|&(_, c)| c == 1);
    let mut keys_buf = Vec::new();
    let mut counts_buf = Vec::new();
    let mut i = 0usize;
    for b in 0..BUCKETS {
        keys_buf.clear();
        counts_buf.clear();
        while i < entries.len() && entries[i].0.bucket_of(k, BUCKET_BITS) == b {
            keys_buf.push(entries[i].0);
            counts_buf.push(entries[i].1);
            i += 1;
        }
        bucket_counts[b] = keys_buf.len() as u32;
        encode_bucket_keys(k, b, &keys_buf, &mut key_payloads[b]);
        if !presence_only {
            encode_bucket_counts(&counts_buf, &mut count_payloads[b]);
        }
    }
    ensure!(
        i == entries.len(),
        "input entries were not sorted by bucket"
    );
    write_ktab_from_buckets(
        path,
        name,
        k,
        min_count,
        max_count_observed,
        presence_only,
        &bucket_counts,
        &key_payloads,
        &count_payloads,
    )
}

/// Random-access reader for Golomb-Rice sorted-key tables (`GRMKGR1\0`).
pub struct GrKsetReader {
    file: File,
    pub k: usize,
    pub key_width: u8,
    pub name: String,
    pub n_kmers: u64,
    offsets: Vec<u64>,
    counts: Vec<u32>,
    data_start: u64,
}

impl GrKsetReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        let mut head = [0u8; 32];
        file.read_exact(&mut head)
            .with_context(|| format!("read header of {}", path.display()))?;
        ensure!(
            &head[..8] == MAGIC_GRKSET,
            "{} is not a GRMKGR1 file",
            path.display()
        );
        let k = head[8] as usize;
        let bucket_bits = head[9] as u32;
        let key_width = head[10];
        ensure!(
            bucket_bits == BUCKET_BITS,
            "unsupported bucket width {bucket_bits}"
        );
        ensure!(
            (MIN_K..=MAX_K).contains(&k),
            "{}: unsupported k {k}",
            path.display()
        );
        let want_width = if k <= MAX_K_U64 { 8 } else { 16 };
        ensure!(
            key_width == want_width,
            "{}: key_width {key_width} does not match k={k}",
            path.display()
        );
        let n_kmers = u64::from_le_bytes(head[16..24].try_into()?);
        let name_len = u32::from_le_bytes(head[24..28].try_into()?) as usize;
        let table_len = (BUCKETS + 1) * 8 + BUCKETS * 4;
        let data_start = 32u64
            .checked_add(name_len as u64)
            .and_then(|v| v.checked_add(table_len as u64))
            .context("header length overflow")?;
        ensure!(data_start <= file_len, "{} is truncated", path.display());
        let mut name_bytes = vec![0u8; name_len];
        file.read_exact(&mut name_bytes)?;
        let mut table = vec![0u8; table_len];
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
                && offsets.windows(2).all(|w| w[0] <= w[1])
                && counts.iter().map(|&c| c as u64).sum::<u64>() == n_kmers,
            "{}: inconsistent GRMKGR1 bucket table",
            path.display()
        );
        Ok(Self {
            file,
            k,
            key_width,
            name: String::from_utf8(name_bytes)?,
            n_kmers,
            offsets,
            counts,
            data_start,
        })
    }

    pub fn bucket_len(&self, bucket: usize) -> usize {
        self.counts[bucket] as usize
    }

    pub fn read_bucket<K: KmerKey>(
        &self,
        bucket: usize,
        buf: &mut Vec<u8>,
        out: &mut Vec<K>,
    ) -> Result<()> {
        ensure!(
            K::WIDTH_BYTES >= self.key_width as usize,
            "cannot decode k={} into {}-byte key",
            self.k,
            K::WIDTH_BYTES
        );
        let (start, end) = (self.offsets[bucket], self.offsets[bucket + 1]);
        buf.clear();
        buf.resize((end - start) as usize, 0);
        if !buf.is_empty() {
            self.file.read_exact_at(buf, self.data_start + start)?;
        }
        decode_bucket_keys(self.k, bucket, self.counts[bucket], buf, out)
    }

    pub fn read_all<K: KmerKey>(&self) -> Result<Vec<K>> {
        let mut out = Vec::with_capacity(self.n_kmers as usize);
        let mut buf = Vec::new();
        for b in 0..BUCKETS {
            self.read_bucket(b, &mut buf, &mut out)?;
        }
        Ok(out)
    }

    /// Reads the `s` smallest mixed keys starting from bucket 0 (1-bucket
    /// MinHash sketch for typical bacterial genomes).
    pub fn sketch<K: KmerKey>(&self, s: usize) -> Result<Vec<K>> {
        let mut out = Vec::with_capacity(s);
        let mut buf = Vec::new();
        for b in 0..BUCKETS {
            if out.len() >= s {
                break;
            }
            self.read_bucket(b, &mut buf, &mut out)?;
        }
        out.truncate(s);
        Ok(out)
    }
}

/// Random-access reader for counted `.ktab` (`GRMKTB1\0`) tables.
pub struct KtabReader {
    file: File,
    pub k: usize,
    pub key_width: u8,
    pub presence_only: bool,
    pub min_count: u32,
    pub max_count_observed: u32,
    pub name: String,
    pub n_kmers: u64,
    key_offsets: Vec<u64>,
    count_offsets: Vec<u64>,
    bucket_counts: Vec<u32>,
    key_data_start: u64,
    count_data_start: u64,
}

impl KtabReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        let mut head = [0u8; 40];
        file.read_exact(&mut head)
            .with_context(|| format!("read header of {}", path.display()))?;
        ensure!(
            &head[..8] == MAGIC_KTAB,
            "{} is not a GRMKTB1 file",
            path.display()
        );
        let k = head[8] as usize;
        let bucket_bits = head[9] as u32;
        let key_width = head[10];
        let flags = head[11];
        ensure!(
            bucket_bits == BUCKET_BITS,
            "unsupported bucket width {bucket_bits}"
        );
        ensure!(
            (MIN_K..=MAX_K).contains(&k),
            "{}: unsupported k {k}",
            path.display()
        );
        let want_width = if k <= MAX_K_U64 { 8 } else { 16 };
        ensure!(
            key_width == want_width,
            "{}: key_width {key_width} does not match k={k}",
            path.display()
        );
        let presence_only = (flags & FLAG_PRESENCE_ONLY) != 0;
        let min_count = u32::from_le_bytes(head[12..16].try_into()?);
        let max_count_observed = u32::from_le_bytes(head[16..20].try_into()?);
        let n_kmers = u64::from_le_bytes(head[24..32].try_into()?);
        let name_len = u32::from_le_bytes(head[32..36].try_into()?) as usize;
        let table_len = 2 * (BUCKETS + 1) * 8 + BUCKETS * 4;
        let key_data_start = 40u64
            .checked_add(name_len as u64)
            .and_then(|v| v.checked_add(table_len as u64))
            .context("header length overflow")?;
        ensure!(
            key_data_start <= file_len,
            "{} is truncated",
            path.display()
        );
        let mut name_bytes = vec![0u8; name_len];
        file.read_exact(&mut name_bytes)?;
        let mut table = vec![0u8; table_len];
        file.read_exact(&mut table)?;
        let (koff_bytes, rest) = table.split_at((BUCKETS + 1) * 8);
        let (coff_bytes, count_bytes) = rest.split_at((BUCKETS + 1) * 8);
        let key_offsets: Vec<u64> = koff_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect();
        let count_offsets: Vec<u64> = coff_bytes
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| u64::from_le_bytes(*c))
            .collect();
        let bucket_counts: Vec<u32> = count_bytes
            .as_chunks::<4>()
            .0
            .iter()
            .map(|c| u32::from_le_bytes(*c))
            .collect();

        let key_data_len = key_offsets[BUCKETS];
        let count_data_start = key_data_start
            .checked_add(key_data_len)
            .context("key data length overflow")?;
        let count_data_len = count_offsets[BUCKETS];
        ensure!(
            key_offsets[0] == 0
                && count_offsets[0] == 0
                && key_offsets.windows(2).all(|w| w[0] <= w[1])
                && count_offsets.windows(2).all(|w| w[0] <= w[1])
                && (!presence_only || count_data_len == 0)
                && count_data_start.checked_add(count_data_len) == Some(file_len)
                && bucket_counts.iter().map(|&c| c as u64).sum::<u64>() == n_kmers,
            "{}: inconsistent GRMKTB1 bucket table",
            path.display()
        );
        Ok(Self {
            file,
            k,
            key_width,
            presence_only,
            min_count,
            max_count_observed,
            name: String::from_utf8(name_bytes)?,
            n_kmers,
            key_offsets,
            count_offsets,
            bucket_counts,
            key_data_start,
            count_data_start,
        })
    }

    pub fn bucket_len(&self, bucket: usize) -> usize {
        self.bucket_counts[bucket] as usize
    }

    pub fn read_bucket_keys<K: KmerKey>(
        &self,
        bucket: usize,
        buf: &mut Vec<u8>,
        out: &mut Vec<K>,
    ) -> Result<()> {
        ensure!(
            K::WIDTH_BYTES >= self.key_width as usize,
            "cannot decode k={} into {}-byte key",
            self.k,
            K::WIDTH_BYTES
        );
        let (start, end) = (self.key_offsets[bucket], self.key_offsets[bucket + 1]);
        buf.clear();
        buf.resize((end - start) as usize, 0);
        if !buf.is_empty() {
            self.file.read_exact_at(buf, self.key_data_start + start)?;
        }
        decode_bucket_keys(self.k, bucket, self.bucket_counts[bucket], buf, out)
    }

    pub fn read_bucket<K: KmerKey>(
        &self,
        bucket: usize,
        key_buf: &mut Vec<u8>,
        count_buf: &mut Vec<u8>,
        out_keys: &mut Vec<K>,
        out_counts: &mut Vec<u32>,
    ) -> Result<()> {
        self.read_bucket_keys(bucket, key_buf, out_keys)?;
        let (c_start, c_end) = (self.count_offsets[bucket], self.count_offsets[bucket + 1]);
        count_buf.clear();
        count_buf.resize((c_end - c_start) as usize, 0);
        if !count_buf.is_empty() {
            self.file
                .read_exact_at(count_buf, self.count_data_start + c_start)?;
        }
        decode_bucket_counts(
            bucket,
            self.bucket_counts[bucket],
            self.presence_only,
            count_buf,
            out_counts,
        )
    }

    pub fn read_all<K: KmerKey>(&self) -> Result<Vec<(K, u32)>> {
        let mut keys = Vec::with_capacity(self.n_kmers as usize);
        let mut counts = Vec::with_capacity(self.n_kmers as usize);
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        for b in 0..BUCKETS {
            self.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        }
        Ok(keys.into_iter().zip(counts).collect())
    }

    pub fn sketch<K: KmerKey>(&self, s: usize) -> Result<Vec<K>> {
        let mut out = Vec::with_capacity(s);
        let mut buf = Vec::new();
        for b in 0..BUCKETS {
            if out.len() >= s {
                break;
            }
            self.read_bucket_keys(b, &mut buf, &mut out)?;
        }
        out.truncate(s);
        Ok(out)
    }
}

/// Unified reader over `.kset` v1 (`GRMKSET1`), Golomb-Rice sorted-key tables
/// (`GRMKGR1\0`), and counted `.ktab` (`GRMKTB1\0`) files.
pub enum AnyTableReader {
    V1(KsetReader),
    Gr(GrKsetReader),
    Tab(KtabReader),
}

impl AnyTableReader {
    pub fn open(path: &Path) -> Result<Self> {
        match detect_format(path)? {
            TableFormat::KsetV1 => Ok(Self::V1(KsetReader::open(path)?)),
            TableFormat::GrKset => Ok(Self::Gr(GrKsetReader::open(path)?)),
            TableFormat::Ktab => Ok(Self::Tab(KtabReader::open(path)?)),
        }
    }

    pub fn format(&self) -> TableFormat {
        match self {
            Self::V1(_) => TableFormat::KsetV1,
            Self::Gr(_) => TableFormat::GrKset,
            Self::Tab(_) => TableFormat::Ktab,
        }
    }

    pub fn k(&self) -> usize {
        match self {
            Self::V1(r) => r.k,
            Self::Gr(r) => r.k,
            Self::Tab(r) => r.k,
        }
    }

    pub fn name(&self) -> &str {
        match self {
            Self::V1(r) => &r.name,
            Self::Gr(r) => &r.name,
            Self::Tab(r) => &r.name,
        }
    }

    pub fn n_kmers(&self) -> u64 {
        match self {
            Self::V1(r) => r.n_kmers,
            Self::Gr(r) => r.n_kmers,
            Self::Tab(r) => r.n_kmers,
        }
    }

    pub fn bucket_len(&self, bucket: usize) -> usize {
        match self {
            Self::V1(r) => r.bucket_len(bucket),
            Self::Gr(r) => r.bucket_len(bucket),
            Self::Tab(r) => r.bucket_len(bucket),
        }
    }

    /// Reads one bucket of `(key, count)` pairs into caller-owned buffers.
    /// Presence-only formats (`V1`, `Gr`) emit count `1` for every present key.
    pub fn read_bucket<K: KmerKey>(
        &self,
        bucket: usize,
        key_buf: &mut Vec<u8>,
        count_buf: &mut Vec<u8>,
        out_keys: &mut Vec<K>,
        out_counts: &mut Vec<u32>,
    ) -> Result<()> {
        match self {
            Self::V1(r) => {
                let mut tmp = Vec::with_capacity(r.bucket_len(bucket));
                r.read_bucket(bucket, key_buf, &mut tmp)?;
                out_keys.extend(tmp.iter().map(|&v| K::from_u128(v as u128)));
                out_counts.resize(out_counts.len() + tmp.len(), 1);
                Ok(())
            }
            Self::Gr(r) => {
                let start_len = out_keys.len();
                r.read_bucket(bucket, key_buf, out_keys)?;
                let added = out_keys.len() - start_len;
                out_counts.resize(out_counts.len() + added, 1);
                Ok(())
            }
            Self::Tab(r) => r.read_bucket(bucket, key_buf, count_buf, out_keys, out_counts),
        }
    }

    pub fn sketch<K: KmerKey>(&self, s: usize) -> Result<Vec<K>> {
        match self {
            Self::V1(r) => Ok(r
                .sketch(s)?
                .into_iter()
                .map(|v| K::from_u128(v as u128))
                .collect()),
            Self::Gr(r) => r.sketch(s),
            Self::Tab(r) => r.sketch(s),
        }
    }
}

/// Converts a `.kset` v1 (`GRMKSET1`) file losslessly to a Golomb-Rice
/// sorted-key table (`GRMKGR1\0`), bucket-by-bucket.
pub fn convert_v1_to_gr(src: &Path, dst: &Path) -> Result<u64> {
    let reader = KsetReader::open(src)?;
    let mut bucket_counts = [0u32; BUCKETS];
    let mut payloads = vec![Vec::new(); BUCKETS];
    let mut buf = Vec::new();
    let mut keys = Vec::new();
    for b in 0..BUCKETS {
        keys.clear();
        reader.read_bucket(b, &mut buf, &mut keys)?;
        bucket_counts[b] = keys.len() as u32;
        encode_bucket_keys(reader.k, b, &keys, &mut payloads[b]);
    }
    write_gr_kset_from_buckets(dst, &reader.name, reader.k, &bucket_counts, &payloads)?;
    Ok(reader.n_kmers)
}

/// Converts a Golomb-Rice sorted-key table (`GRMKGR1\0`) or counted `.ktab`
/// (`GRMKTB1\0`) with `k <= 32` losslessly to `.kset` v1 (`GRMKSET1`).
pub fn convert_to_v1(src: &Path, dst: &Path) -> Result<u64> {
    let reader = AnyTableReader::open(src)?;
    ensure!(
        reader.k() <= MAX_K_U64,
        "cannot convert k={} to GRMKSET1 (max k is {MAX_K_U64})",
        reader.k()
    );
    let space = KmerSpace::new(reader.k());
    let mut all = Vec::with_capacity(reader.n_kmers() as usize);
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut counts = Vec::new();
    for b in 0..BUCKETS {
        counts.clear();
        reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut all, &mut counts)?;
    }
    kset::write(dst, reader.name(), &space, &all)?;
    Ok(reader.n_kmers())
}

// ---------------------------------------------------------------------------
// FastK-style compact per-sequence k-mer count profile (.prof.bin, GRMPRF1\0)
// ---------------------------------------------------------------------------

/// One sequence's positional k-mer count profile.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SequenceProfile {
    pub name: String,
    pub seq_len: u64,
    pub counts: Vec<u32>,
}

pub(crate) fn put_uleb(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

pub(crate) fn get_uleb(data: &[u8], pos: &mut usize) -> Result<u64> {
    let mut v = 0u64;
    let mut shift = 0u32;
    loop {
        let Some(&b) = data.get(*pos) else {
            bail!("truncated varint");
        };
        *pos += 1;
        v |= ((b & 0x7f) as u64) << shift;
        if b < 0x80 {
            return Ok(v);
        }
        shift += 7;
        ensure!(shift < 64, "varint overflow");
    }
}

#[inline]
fn zigzag_encode(delta: i64) -> u64 {
    ((delta << 1) ^ (delta >> 63)) as u64
}

#[inline]
fn zigzag_decode(z: u64) -> i64 {
    ((z >> 1) as i64) ^ -((z & 1) as i64)
}

fn encode_profile_runs(counts: &[u32], out: &mut Vec<u8>) {
    out.clear();
    let mut prev: i64 = 0;
    let mut i = 0usize;
    while i < counts.len() {
        let val = counts[i];
        let mut run = 1u64;
        while i + (run as usize) < counts.len() && counts[i + (run as usize)] == val {
            run += 1;
        }
        let delta = (val as i64) - prev;
        put_uleb(out, zigzag_encode(delta));
        put_uleb(out, run - 1);
        prev = val as i64;
        i += run as usize;
    }
}

fn decode_profile_runs(bytes: &[u8], n_positions: usize) -> Result<Vec<u32>> {
    let mut out = Vec::with_capacity(n_positions);
    let mut pos = 0usize;
    let mut prev: i64 = 0;
    while pos < bytes.len() {
        let delta = zigzag_decode(get_uleb(bytes, &mut pos)?);
        let run = get_uleb(bytes, &mut pos)?
            .checked_add(1)
            .context("profile run overflow")?;
        let val = prev
            .checked_add(delta)
            .and_then(|v| u32::try_from(v).ok())
            .context("profile count out of u32 range")?;
        ensure!(
            out.len().saturating_add(run as usize) <= n_positions,
            "profile decoded more positions than header declared"
        );
        out.resize(out.len() + run as usize, val);
        prev = val as i64;
    }
    ensure!(
        out.len() == n_positions,
        "profile decoded {} positions, expected {n_positions}",
        out.len()
    );
    Ok(out)
}

/// Writes per-sequence positional k-mer count profiles into a compact
/// run-length + zigzag-delta binary file (`GRMPRF1\0`).
pub fn write_profiles(path: &Path, k: usize, profiles: &[SequenceProfile]) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    let part = crate::confine::suffixed(path, ".part");
    let mut w =
        BufWriter::new(File::create(&part).with_context(|| format!("create {}", part.display()))?);
    w.write_all(MAGIC_PROF1)?;
    w.write_all(&[k as u8, 0, 0, 0])?;
    w.write_all(&(profiles.len() as u64).to_le_bytes())?;
    let mut enc = Vec::new();
    for p in profiles {
        encode_profile_runs(&p.counts, &mut enc);
        w.write_all(&(p.name.len() as u32).to_le_bytes())?;
        w.write_all(p.name.as_bytes())?;
        w.write_all(&p.seq_len.to_le_bytes())?;
        w.write_all(&(p.counts.len() as u64).to_le_bytes())?;
        w.write_all(&(enc.len() as u64).to_le_bytes())?;
        w.write_all(&enc)?;
    }
    w.flush()?;
    drop(w);
    std::fs::rename(&part, path).with_context(|| format!("rename to {}", path.display()))?;
    Ok(())
}

/// Reads all sequence profiles from a `.prof.bin` (`GRMPRF1\0`) file.
pub fn read_profiles(path: &Path) -> Result<(usize, Vec<SequenceProfile>)> {
    let mut r =
        BufReader::new(File::open(path).with_context(|| format!("open {}", path.display()))?);
    let mut head = [0u8; 20];
    r.read_exact(&mut head)
        .with_context(|| format!("read profile header {}", path.display()))?;
    ensure!(
        &head[..8] == MAGIC_PROF1,
        "{} is not a GRMPRF1 profile file",
        path.display()
    );
    let k = head[8] as usize;
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported profile k {k}");
    let n_seq = u64::from_le_bytes(head[12..20].try_into()?) as usize;
    let mut profiles = Vec::with_capacity(n_seq);
    for _ in 0..n_seq {
        let mut len_buf = [0u8; 4];
        r.read_exact(&mut len_buf)?;
        let name_len = u32::from_le_bytes(len_buf) as usize;
        let mut name_bytes = vec![0u8; name_len];
        r.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)?;
        let mut meta = [0u8; 24];
        r.read_exact(&mut meta)?;
        let seq_len = u64::from_le_bytes(meta[0..8].try_into()?);
        let n_pos = u64::from_le_bytes(meta[8..16].try_into()?) as usize;
        let payload_len = u64::from_le_bytes(meta[16..24].try_into()?) as usize;
        let mut payload = vec![0u8; payload_len];
        r.read_exact(&mut payload)?;
        let counts = decode_profile_runs(&payload, n_pos)?;
        profiles.push(SequenceProfile {
            name,
            seq_len,
            counts,
        });
    }
    let mut extra = [0u8; 1];
    ensure!(
        r.read(&mut extra)? == 0,
        "{}: trailing bytes after last sequence profile",
        path.display()
    );
    Ok((k, profiles))
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn gr_kset_round_trips_edges_and_sketches_for_u64_and_u128() {
        let dir = tempfile::tempdir().unwrap();
        let mut rng = ChaCha8Rng::seed_from_u64(505);

        for k in [5, 21, 31, 32] {
            let mask = if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            };
            let mut vals: Vec<u64> = vec![0, 1, mask - 1, mask];
            vals.extend((0..10_000).map(|_| rng.random::<u64>() & mask));
            vals.sort_unstable();
            vals.dedup();
            let p = dir.path().join(format!("k{k}.grk"));
            write_gr_kset(&p, "sample", k, &vals).unwrap();
            let r = GrKsetReader::open(&p).unwrap();
            assert_eq!(r.k, k);
            assert_eq!(r.n_kmers as usize, vals.len());
            assert_eq!(r.read_all::<u64>().unwrap(), vals);
            let want_s = vals.len().min(64);
            assert_eq!(r.sketch::<u64>(64).unwrap(), vals[..want_s]);
        }

        for k in [33, 61, 64] {
            let mask = if k == 64 {
                u128::MAX
            } else {
                (1u128 << (2 * k)) - 1
            };
            let mut vals: Vec<u128> = vec![0, 1, mask - 1, mask];
            vals.extend((0..10_000).map(|_| {
                (((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128)) & mask
            }));
            vals.sort_unstable();
            vals.dedup();
            let p = dir.path().join(format!("k{k}.grk"));
            write_gr_kset(&p, "wide", k, &vals).unwrap();
            let r = GrKsetReader::open(&p).unwrap();
            assert_eq!(r.k, k);
            assert_eq!(r.key_width, 16);
            assert_eq!(r.read_all::<u128>().unwrap(), vals);
            assert_eq!(r.sketch::<u128>(50).unwrap(), vals[..50]);
        }
    }

    #[test]
    fn ktab_preserves_uncapped_u32_counts_and_presence_only_zero_overhead() {
        let dir = tempfile::tempdir().unwrap();
        let entries: Vec<(u64, u32)> = vec![
            (10, 1),
            (20, 2),
            (30, 3),
            (40, 253),
            (50, 254),
            (60, 255),
            (70, 256),
            (80, 300),
            (90, 40_000),
            (100, 65_534),
            (110, 65_535),
            (120, 70_000),
            (130, 1_500_000),
            (140, u32::MAX),
        ];
        let p = dir.path().join("counts.ktab");
        write_ktab(&p, "c", 31, 1, &entries).unwrap();
        let r = KtabReader::open(&p).unwrap();
        assert!(!r.presence_only);
        assert_eq!(r.max_count_observed, u32::MAX);
        assert_eq!(r.read_all::<u64>().unwrap(), entries);

        let ones: Vec<(u64, u32)> = (0..2_000u64).map(|i| (i * 17, 1)).collect();
        let p_ones = dir.path().join("ones.ktab");
        write_ktab(&p_ones, "o", 31, 1, &ones).unwrap();
        let ro = KtabReader::open(&p_ones).unwrap();
        assert!(ro.presence_only);
        assert_eq!(ro.read_all::<u64>().unwrap(), ones);
    }

    #[test]
    fn profile_round_trips_runs_and_deltas() {
        let dir = tempfile::tempdir().unwrap();
        let profs = vec![
            SequenceProfile {
                name: "contig1".to_string(),
                seq_len: 100,
                counts: vec![1, 1, 1, 1, 0, 0, 5, 5, 5, 300, 70_000, 1],
            },
            SequenceProfile {
                name: "short".to_string(),
                seq_len: 10,
                counts: vec![],
            },
        ];
        let p = dir.path().join("test.prof.bin");
        write_profiles(&p, 31, &profs).unwrap();
        let (k, got) = read_profiles(&p).unwrap();
        assert_eq!(k, 31);
        assert_eq!(got, profs);
    }
}
