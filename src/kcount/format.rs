//! Verified binary k-mer storage formats for `grm kcount`:
//!
//! - Golomb-Rice sorted-key table (`.grk`, magic `GRMKGR1\0`): presence-only
//!   1,024-bucket Golomb-Rice encoding of sorted mixed keys for `k = 5..=64`,
//!   with input modality metadata and a 32-byte SHA-256 header/body digest.
//! - Counted table (`.ktab`, magic `GRMKTB1\0`): 1,024-bucket Golomb-Rice
//!   sorted mixed keys paired with `u32` counts capped at `2^32 - 1`
//!   (`4,294,967,295`), recording `saturated_counts`, input modality, and a
//!   32-byte SHA-256 header/body digest. When all counts equal 1, `.ktab` uses
//!   presence-only mode with 8,208 bytes of extra header/table overhead over
//!   `.grk` and zero count-payload bytes.
//! - Per-sequence positional k-mer count profile (`.prof.bin`, magic
//!   `GRMPRF1\0`): run-length and zigzag-delta encoded counts along each
//!   sequence, with summary totals and a 32-byte SHA-256 header/body digest.
//! - Unified read interop with untouched `.kset` v1 (`GRMKSET1`) via
//!   [`AnyTableReader`].

use crate::confine;
use crate::dna::KmerSpace;
use crate::kcount::key::{KmerKey, MAX_K, MAX_K_U64, MIN_K};
use crate::kset::{self, BUCKET_BITS, BUCKETS, KsetReader};
use anyhow::{Context, Result, bail, ensure};
use ring::digest::{Context as Sha256Context, SHA256};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::os::unix::fs::FileExt;
use std::path::Path;

pub const MAGIC_GRKSET: &[u8; 8] = b"GRMKGR1\0";
pub const MAGIC_KTAB: &[u8; 8] = b"GRMKTB1\0";
pub const MAGIC_PROF1: &[u8; 8] = b"GRMPRF1\0";

/// Default sample name written into table headers when `--name` is omitted,
/// ensuring input file names are never copied into output headers by default.
pub const DEFAULT_SAMPLE_NAME: &str = "sample";

/// Maximum representable k-mer count in `.ktab` and `.prof.bin` (`2^32 - 1`).
pub const MAX_COUNT: u32 = u32::MAX;

const FLAG_PRESENCE_ONLY: u8 = 0x01;
const FLAG_SATURATED_COUNTS: u8 = 0x02;
const VALID_KTAB_FLAGS: u8 = FLAG_PRESENCE_ONLY | FLAG_SATURATED_COUNTS;

const COUNT_MODE_ALL_ONE: u8 = 0;
const COUNT_MODE_TIERED: u8 = 1;
const COUNT_MODE_U16_ESC: u8 = 2;

const GRK_HEADER_LEN: usize = 72;
const GRK_DIGEST_START: usize = 40;
const GRK_DIGEST_END: usize = 72;

const KTAB_HEADER_LEN: usize = 80;
const KTAB_DIGEST_START: usize = 48;
const KTAB_DIGEST_END: usize = 80;

const PROF_HEADER_LEN: usize = 72;
const PROF_DIGEST_START: usize = 40;
const PROF_DIGEST_END: usize = 72;
const PROF_MIN_RECORD_BYTES: u64 = 4 + 8 + 8 + 8;

/// Cap on eager `Vec::with_capacity` hints before reading payloads; vectors
/// grow safely via validated per-bucket/per-record reserves.
const MAX_EAGER_CAP: usize = 1 << 20;

/// Input sequencing modality stored as a 1-byte field in `.grk` and `.ktab`
/// headers so downstream analyses verify the table's recorded modality directly.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InputModality {
    #[default]
    Unknown = 0,
    Assembly = 1,
    Reads = 2,
}

impl InputModality {
    pub fn from_byte(b: u8) -> Result<Self> {
        match b {
            0 => Ok(Self::Unknown),
            1 => Ok(Self::Assembly),
            2 => Ok(Self::Reads),
            other => bail!("unrecognized input modality byte {other}"),
        }
    }

    #[inline]
    pub fn as_byte(self) -> u8 {
        self as u8
    }

    #[inline]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unknown => "unknown",
            Self::Assembly => "assembly",
            Self::Reads => "reads",
        }
    }
}

/// Identifies which random-access binary k-mer table format a file holds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TableFormat {
    KsetV1,
    GrKset,
    Ktab,
}

impl TableFormat {
    pub fn expected_extension(self) -> &'static str {
        match self {
            Self::KsetV1 => "kset",
            Self::GrKset => "grk",
            Self::Ktab => "ktab",
        }
    }

    pub fn magic_label(self) -> &'static str {
        match self {
            Self::KsetV1 => "GRMKSET1",
            Self::GrKset => "GRMKGR1",
            Self::Ktab => "GRMKTB1",
        }
    }
}

/// Ensures `path`'s file extension matches the target table format so a writer
/// never writes one format's bytes under another format's suffix.
pub fn check_output_extension(path: &Path, format: TableFormat) -> Result<()> {
    let want = format.expected_extension();
    let got = path.extension().and_then(|s| s.to_str());
    ensure!(
        got == Some(want),
        "output path {} has extension {:?}, which disagrees with {} (expected .{want})",
        path.display(),
        got,
        format.magic_label()
    );
    Ok(())
}

/// Formats bytes as lowercase hexadecimal.
pub fn hex_bytes(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut out = String::with_capacity(bytes.len() * 2);
    for &b in bytes {
        out.push(HEX[(b >> 4) as usize] as char);
        out.push(HEX[(b & 0x0f) as usize] as char);
    }
    out
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

/// Allocates a zero-filled byte buffer using fallible `try_reserve_exact` so
/// damaged files can never trigger an OOM abort or capacity-overflow panic.
fn alloc_zeroed_checked(len: usize, label: &str) -> Result<Vec<u8>> {
    let mut buf = Vec::new();
    buf.try_reserve_exact(len)
        .with_context(|| format!("{label}: allocation of {len} bytes failed"))?;
    buf.resize(len, 0);
    Ok(buf)
}

/// Streams a file in fixed 64 KiB chunks through SHA-256, covering every byte
/// of the file except the 32-byte digest slot `[slot_start..slot_end]`.
fn compute_file_sha256_excluding_slot(
    file: &File,
    file_len: u64,
    slot_start: usize,
    slot_end: usize,
) -> Result<[u8; 32]> {
    ensure!(
        slot_start <= slot_end && (slot_end as u64) <= file_len,
        "file truncated before digest slot"
    );
    let mut sha = Sha256Context::new(&SHA256);
    let mut buf = [0u8; 64 * 1024];

    let mut pos = 0u64;
    let prefix_end = slot_start as u64;
    while pos < prefix_end {
        let take = ((prefix_end - pos) as usize).min(buf.len());
        file.read_exact_at(&mut buf[..take], pos)?;
        sha.update(&buf[..take]);
        pos += take as u64;
    }

    pos = slot_end as u64;
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
        let q = gap >> best_r;
        debug_assert!(q <= u64::MAX as u128);
        bw.write_unary(q as u64);
        bw.write_u128(gap & rem_mask, best_r);
    }
    bw.finish();
}

/// Decodes one bucket's Golomb-Rice payload into `out`, validating length
/// bounds before allocating, strict monotonicity, bucket range, and payload
/// boundaries.
pub fn decode_bucket_keys<K: KmerKey>(
    k: usize,
    bucket: usize,
    count: u32,
    bytes: &[u8],
    out: &mut Vec<K>,
) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    ensure!(bucket < BUCKETS, "bucket index {bucket} out of range");
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
    // Each key takes at least (r + 1) bits (1 unary stop bit + r remainder bits).
    let min_bits = (count as u128)
        .checked_mul((r + 1) as u128)
        .context("bucket bit count overflow")?;
    let avail_bits = ((bytes.len() - 1) as u128)
        .checked_mul(8)
        .context("bucket byte count overflow")?;
    ensure!(
        min_bits <= avail_bits,
        "bucket {bucket} payload ({} bytes) cannot hold {count} keys at r={r}",
        bytes.len()
    );
    let floor = K::bucket_floor(bucket, k, BUCKET_BITS).to_u128();
    let max_val = if low_bits == 128 {
        u128::MAX
    } else {
        floor | ((1u128 << low_bits) - 1)
    };
    let mut br = BitReader::new(&bytes[1..]);
    let mut prev = floor;
    out.try_reserve(count as usize)
        .with_context(|| format!("bucket {bucket}: key buffer allocation failed"))?;
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

/// Encodes one bucket's `u32` counts (capped at `2^32 - 1`) with automatic
/// per-bucket mode selection (1-byte all-ones header, tiered 1B/3B/5B for
/// low-depth counts, or 2B `u16` + `0xFFFF` escape to `u32` for deep counts).
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

/// Decodes one bucket's `u32` counts into `out`, validating byte bounds before
/// allocating.
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
    let n_usize = usize::try_from(n).context("bucket count exceeds addressable memory")?;
    if presence_only {
        ensure!(
            bytes.is_empty(),
            "presence-only bucket {bucket} carries count bytes"
        );
        out.try_reserve(n_usize)
            .with_context(|| format!("bucket {bucket}: presence count allocation failed"))?;
        out.resize(out.len() + n_usize, 1);
        return Ok(());
    }
    ensure!(
        !bytes.is_empty(),
        "bucket {bucket} missing count mode header"
    );
    let mode = bytes[0];
    let payload_len = bytes.len() - 1;
    let mut pos = 1usize;
    match mode {
        COUNT_MODE_ALL_ONE => {
            ensure!(
                bytes.len() == 1,
                "all-ones count bucket {bucket} has trailing bytes"
            );
            out.try_reserve(n_usize)
                .with_context(|| format!("bucket {bucket}: all-ones count allocation failed"))?;
            out.resize(out.len() + n_usize, 1);
        }
        COUNT_MODE_TIERED => {
            ensure!(
                payload_len >= n_usize,
                "bucket {bucket} tiered count payload ({payload_len} bytes) too short for {n} counts"
            );
            out.try_reserve(n_usize)
                .with_context(|| format!("bucket {bucket}: tiered count allocation failed"))?;
            for _ in 0..n {
                let Some(&tag) = bytes.get(pos) else {
                    bail!("bucket {bucket} truncated tiered count");
                };
                pos += 1;
                let c = if tag < 254 {
                    tag as u32
                } else if tag == 254 {
                    ensure!(
                        pos.checked_add(2).is_some_and(|e| e <= bytes.len()),
                        "bucket {bucket} truncated u16"
                    );
                    let v = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap()) as u32;
                    pos += 2;
                    v
                } else {
                    ensure!(
                        pos.checked_add(4).is_some_and(|e| e <= bytes.len()),
                        "bucket {bucket} truncated u32"
                    );
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
            let min_bytes = (n as u64)
                .checked_mul(2)
                .context("u16 count byte length overflow")?;
            ensure!(
                (payload_len as u64) >= min_bytes,
                "bucket {bucket} u16 count payload ({payload_len} bytes) too short for {n} counts"
            );
            out.try_reserve(n_usize)
                .with_context(|| format!("bucket {bucket}: u16 count allocation failed"))?;
            for _ in 0..n {
                ensure!(
                    pos.checked_add(2).is_some_and(|e| e <= bytes.len()),
                    "bucket {bucket} truncated u16"
                );
                let v = u16::from_le_bytes(bytes[pos..pos + 2].try_into().unwrap());
                pos += 2;
                let c = if v < 0xFFFF {
                    v as u32
                } else {
                    ensure!(
                        pos.checked_add(4).is_some_and(|e| e <= bytes.len()),
                        "bucket {bucket} truncated u32 esc"
                    );
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

/// Writes pre-encoded 1,024 bucket payloads as a verified Golomb-Rice
/// sorted-key table (`GRMKGR1\0`) atomically via `<path>.part` and rename.
pub fn write_gr_kset_from_buckets(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    bucket_counts: &[u32; BUCKETS],
    bucket_payloads: &[Vec<u8>],
) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    ensure!(
        bucket_payloads.len() == BUCKETS,
        "expected {BUCKETS} bucket payloads"
    );
    let name_len = u32::try_from(name.len()).context("sample name exceeds u32 length")?;
    let n_kmers = bucket_counts
        .iter()
        .try_fold(0u64, |acc, &c| acc.checked_add(c as u64))
        .context("k-mer count overflow")?;
    let key_width: u8 = if k <= MAX_K_U64 { 8 } else { 16 };
    let mut offsets = [0u64; BUCKETS + 1];
    let mut cur = 0u64;
    for (b, payload) in bucket_payloads.iter().enumerate() {
        offsets[b] = cur;
        cur = cur
            .checked_add(payload.len() as u64)
            .context("bucket offset overflow")?;
    }
    offsets[BUCKETS] = cur;

    // Fixed header prefix [0..40]; bytes [40..72] hold the SHA-256 digest of
    // [0..40] || name || offsets || bucket_counts || bucket_payloads.
    let mut prefix = [0u8; GRK_DIGEST_START];
    prefix[0..8].copy_from_slice(MAGIC_GRKSET);
    prefix[8] = k as u8;
    prefix[9] = BUCKET_BITS as u8;
    prefix[10] = key_width;
    prefix[11] = modality.as_byte();
    // [12..16] reserved0 = 0
    prefix[16..24].copy_from_slice(&n_kmers.to_le_bytes());
    prefix[24..28].copy_from_slice(&name_len.to_le_bytes());
    // [28..40] reserved1 = 0

    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&prefix);
    sha.update(name.as_bytes());
    for o in offsets {
        sha.update(&o.to_le_bytes());
    }
    for &c in bucket_counts {
        sha.update(&c.to_le_bytes());
    }
    for payload in bucket_payloads {
        sha.update(payload);
    }
    let digest = sha.finish();

    let part = confine::suffixed(path, ".part");
    confine::regular_or_absent(path)?;
    let mut w = BufWriter::new(confine::create_part(&part)?);
    let written = (|| -> Result<()> {
        w.write_all(&prefix)?;
        w.write_all(digest.as_ref())?;
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
        confine::replace_with(&part, path)
    })();
    drop(w);
    if written.is_err() && confine::regular_or_absent(&part).unwrap_or(false) {
        let _ = std::fs::remove_file(&part);
    }
    written
}

/// Writes sorted distinct mixed keys as a verified Golomb-Rice sorted-key
/// table (`GRMKGR1\0`).
pub fn write_gr_kset<K: KmerKey>(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    values: &[K],
) -> Result<()> {
    let mut bucket_counts = [0u32; BUCKETS];
    let mut bucket_payloads = vec![Vec::new(); BUCKETS];
    let mut i = 0usize;
    for b in 0..BUCKETS {
        let start = i;
        while i < values.len() && values[i].bucket_of(k, BUCKET_BITS) == b {
            i += 1;
        }
        let slice = &values[start..i];
        bucket_counts[b] = u32::try_from(slice.len()).context("bucket count exceeds u32")?;
        encode_bucket_keys(k, b, slice, &mut bucket_payloads[b]);
    }
    ensure!(i == values.len(), "input keys were not sorted by bucket");
    write_gr_kset_from_buckets(path, name, k, modality, &bucket_counts, &bucket_payloads)
}

/// Header metadata passed to [`write_ktab_from_buckets`].
#[derive(Clone, Copy, Debug)]
pub struct KtabWriteMeta<'a> {
    pub name: &'a str,
    pub k: usize,
    pub modality: InputModality,
    pub min_count: u32,
    pub max_count_observed: u32,
    pub saturated_counts: u64,
    pub presence_only: bool,
}

/// Writes pre-encoded 1,024 bucket key and count payloads as a verified `.ktab`
/// (`GRMKTB1\0`) file atomically via `<path>.part` and rename.
pub fn write_ktab_from_buckets(
    path: &Path,
    meta: &KtabWriteMeta<'_>,
    bucket_counts: &[u32; BUCKETS],
    key_payloads: &[Vec<u8>],
    count_payloads: &[Vec<u8>],
) -> Result<()> {
    ensure!(
        (MIN_K..=MAX_K).contains(&meta.k),
        "unsupported k {}",
        meta.k
    );
    ensure!(
        key_payloads.len() == BUCKETS && count_payloads.len() == BUCKETS,
        "expected {BUCKETS} bucket payloads"
    );
    let name_len = u32::try_from(meta.name.len()).context("sample name exceeds u32 length")?;
    let n_kmers = bucket_counts
        .iter()
        .try_fold(0u64, |acc, &c| acc.checked_add(c as u64))
        .context("k-mer count overflow")?;
    ensure!(
        meta.saturated_counts <= n_kmers,
        "saturated_counts {} exceeds n_kmers {n_kmers}",
        meta.saturated_counts
    );
    if meta.presence_only {
        ensure!(
            meta.saturated_counts == 0 && meta.max_count_observed <= 1,
            "presence-only ktab cannot carry counts above 1 or saturated counts"
        );
    }
    if meta.saturated_counts > 0 {
        ensure!(
            meta.max_count_observed == MAX_COUNT,
            "saturated_counts > 0 requires max_count_observed == u32::MAX"
        );
    }
    if n_kmers == 0 {
        ensure!(
            meta.max_count_observed == 0,
            "empty ktab must have max_count_observed == 0"
        );
    }

    let key_width: u8 = if meta.k <= MAX_K_U64 { 8 } else { 16 };
    let mut flags: u8 = 0;
    if meta.presence_only {
        flags |= FLAG_PRESENCE_ONLY;
    }
    if meta.saturated_counts > 0 {
        flags |= FLAG_SATURATED_COUNTS;
    }

    let mut key_offsets = [0u64; BUCKETS + 1];
    let mut cur = 0u64;
    for (b, payload) in key_payloads.iter().enumerate() {
        key_offsets[b] = cur;
        cur = cur
            .checked_add(payload.len() as u64)
            .context("key offset overflow")?;
    }
    key_offsets[BUCKETS] = cur;

    let mut count_offsets = [0u64; BUCKETS + 1];
    let mut c_cur = 0u64;
    if !meta.presence_only {
        for (b, payload) in count_payloads.iter().enumerate() {
            count_offsets[b] = c_cur;
            c_cur = c_cur
                .checked_add(payload.len() as u64)
                .context("count offset overflow")?;
        }
        count_offsets[BUCKETS] = c_cur;
    }

    // Fixed header prefix [0..48]; bytes [48..80] hold the SHA-256 digest of
    // [0..48] || name || key_offsets || count_offsets || bucket_counts ||
    // key_payloads || count_payloads.
    let mut prefix = [0u8; KTAB_DIGEST_START];
    prefix[0..8].copy_from_slice(MAGIC_KTAB);
    prefix[8] = meta.k as u8;
    prefix[9] = BUCKET_BITS as u8;
    prefix[10] = key_width;
    prefix[11] = flags;
    prefix[12] = meta.modality.as_byte();
    // [13..16] reserved0 = 0
    prefix[16..20].copy_from_slice(&meta.min_count.to_le_bytes());
    prefix[20..24].copy_from_slice(&meta.max_count_observed.to_le_bytes());
    prefix[24..32].copy_from_slice(&n_kmers.to_le_bytes());
    prefix[32..40].copy_from_slice(&meta.saturated_counts.to_le_bytes());
    prefix[40..44].copy_from_slice(&name_len.to_le_bytes());
    // [44..48] reserved1 = 0

    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&prefix);
    sha.update(meta.name.as_bytes());
    for o in key_offsets {
        sha.update(&o.to_le_bytes());
    }
    for o in count_offsets {
        sha.update(&o.to_le_bytes());
    }
    for &c in bucket_counts {
        sha.update(&c.to_le_bytes());
    }
    for payload in key_payloads {
        sha.update(payload);
    }
    if !meta.presence_only {
        for payload in count_payloads {
            sha.update(payload);
        }
    }
    let digest = sha.finish();

    let part = confine::suffixed(path, ".part");
    confine::regular_or_absent(path)?;
    let mut w = BufWriter::new(confine::create_part(&part)?);
    let written = (|| -> Result<()> {
        w.write_all(&prefix)?;
        w.write_all(digest.as_ref())?;
        w.write_all(meta.name.as_bytes())?;
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
        if !meta.presence_only {
            for payload in count_payloads {
                w.write_all(payload)?;
            }
        }
        w.flush()?;
        confine::replace_with(&part, path)
    })();
    drop(w);
    if written.is_err() && confine::regular_or_absent(&part).unwrap_or(false) {
        let _ = std::fs::remove_file(&part);
    }
    written
}

/// Writes sorted `(mixed_key, count)` pairs as a verified `.ktab` (`GRMKTB1\0`).
pub fn write_ktab<K: KmerKey>(
    path: &Path,
    name: &str,
    k: usize,
    modality: InputModality,
    min_count: u32,
    saturated_counts: u64,
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
        bucket_counts[b] = u32::try_from(keys_buf.len()).context("bucket count exceeds u32")?;
        encode_bucket_keys(k, b, &keys_buf, &mut key_payloads[b]);
        if !presence_only {
            encode_bucket_counts(&counts_buf, &mut count_payloads[b]);
        }
    }
    ensure!(
        i == entries.len(),
        "input entries were not sorted by bucket"
    );
    let meta = KtabWriteMeta {
        name,
        k,
        modality,
        min_count,
        max_count_observed,
        saturated_counts,
        presence_only,
    };
    write_ktab_from_buckets(path, &meta, &bucket_counts, &key_payloads, &count_payloads)
}

/// Random-access reader for verified Golomb-Rice sorted-key tables (`GRMKGR1\0`).
pub struct GrKsetReader {
    file: File,
    file_len: u64,
    pub k: usize,
    pub key_width: u8,
    pub modality: InputModality,
    pub name: String,
    pub n_kmers: u64,
    pub digest_sha256: [u8; 32],
    offsets: Vec<u64>,
    counts: Vec<u32>,
    data_start: u64,
}

impl GrKsetReader {
    /// Opens a `.grk` file and validates all header fields and bucket table
    /// bounds against the file length without scanning the payload body.
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        ensure!(
            file_len >= GRK_HEADER_LEN as u64,
            "{} is truncated before GRMKGR1 header",
            path.display()
        );
        let mut head = [0u8; GRK_HEADER_LEN];
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
        let modality = InputModality::from_byte(head[11])
            .with_context(|| format!("{}: invalid modality", path.display()))?;
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
        ensure!(
            head[12..16] == [0u8; 4] && head[28..40] == [0u8; 12],
            "{}: non-zero reserved header bytes in GRMKGR1",
            path.display()
        );
        let n_kmers = u64::from_le_bytes(head[16..24].try_into()?);
        let name_len_u32 = u32::from_le_bytes(head[24..28].try_into()?);
        let mut digest_sha256 = [0u8; 32];
        digest_sha256.copy_from_slice(&head[GRK_DIGEST_START..GRK_DIGEST_END]);

        let table_len = (BUCKETS + 1) * 8 + BUCKETS * 4;
        let data_start = (GRK_HEADER_LEN as u64)
            .checked_add(name_len_u32 as u64)
            .and_then(|v| v.checked_add(table_len as u64))
            .context("header length overflow")?;
        ensure!(data_start <= file_len, "{} is truncated", path.display());
        let name_len = usize::try_from(name_len_u32).context("name length overflow")?;
        let mut name_bytes = alloc_zeroed_checked(name_len, "GRMKGR1 name")?;
        file.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)
            .with_context(|| format!("{}: sample name is not valid UTF-8", path.display()))?;
        let mut table = alloc_zeroed_checked(table_len, "GRMKGR1 bucket table")?;
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

        let payload_total = file_len - data_start;
        ensure!(
            offsets[0] == 0 && offsets[BUCKETS] == payload_total,
            "{}: inconsistent GRMKGR1 bucket offsets",
            path.display()
        );
        let mut sum_counts = 0u64;
        for b in 0..BUCKETS {
            let plen = offsets[b + 1]
                .checked_sub(offsets[b])
                .with_context(|| format!("{}: non-monotonic bucket offsets", path.display()))?;
            let c = counts[b] as u64;
            if c == 0 {
                ensure!(
                    plen == 0,
                    "{}: empty bucket {b} has non-zero payload length {plen}",
                    path.display()
                );
            } else {
                let min_plen = 1u64
                    .checked_add(c.div_ceil(8))
                    .context("bucket minimum payload overflow")?;
                ensure!(
                    plen >= min_plen,
                    "{}: bucket {b} payload ({plen} bytes) too short for {c} keys",
                    path.display()
                );
            }
            sum_counts = sum_counts
                .checked_add(c)
                .context("bucket count sum overflow")?;
        }
        ensure!(
            sum_counts == n_kmers,
            "{}: bucket counts sum ({sum_counts}) != header n_kmers ({n_kmers})",
            path.display()
        );

        Ok(Self {
            file,
            file_len,
            k,
            key_width,
            modality,
            name,
            n_kmers,
            digest_sha256,
            offsets,
            counts,
            data_start,
        })
    }

    /// Opens a `.grk` file and immediately verifies its full SHA-256 digest.
    pub fn open_verified(path: &Path) -> Result<Self> {
        let reader = Self::open(path)?;
        reader.verify_integrity()?;
        Ok(reader)
    }

    /// Recomputes the SHA-256 digest over header metadata, bucket tables, and
    /// payload bytes and ensures it matches the digest stored in the header.
    pub fn verify_integrity(&self) -> Result<[u8; 32]> {
        let got = compute_file_sha256_excluding_slot(
            &self.file,
            self.file_len,
            GRK_DIGEST_START,
            GRK_DIGEST_END,
        )?;
        ensure!(
            got == self.digest_sha256,
            "GRMKGR1 SHA-256 integrity check failed: header has {}, computed {}",
            hex_bytes(&self.digest_sha256),
            hex_bytes(&got)
        );
        Ok(got)
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
        ensure!(bucket < BUCKETS, "bucket index {bucket} out of range");
        ensure!(
            K::WIDTH_BYTES >= self.key_width as usize,
            "cannot decode k={} into {}-byte key",
            self.k,
            K::WIDTH_BYTES
        );
        let (start, end) = (self.offsets[bucket], self.offsets[bucket + 1]);
        let len = usize::try_from(end - start).context("bucket payload length overflow")?;
        buf.clear();
        buf.try_reserve(len)
            .with_context(|| format!("bucket {bucket}: payload buffer allocation failed"))?;
        buf.resize(len, 0);
        if !buf.is_empty() {
            let offset = self
                .data_start
                .checked_add(start)
                .context("bucket file offset overflow")?;
            self.file.read_exact_at(buf, offset)?;
        }
        decode_bucket_keys(self.k, bucket, self.counts[bucket], buf, out)
    }

    /// Verifies the stored SHA-256 digest and decodes every bucket in order.
    pub fn read_all<K: KmerKey>(&self) -> Result<Vec<K>> {
        self.verify_integrity()?;
        let cap = usize::try_from(self.n_kmers)
            .context("n_kmers exceeds addressable memory")?
            .min(MAX_EAGER_CAP);
        let mut out = Vec::new();
        out.try_reserve(cap).context("key list allocation failed")?;
        let mut buf = Vec::new();
        for b in 0..BUCKETS {
            self.read_bucket(b, &mut buf, &mut out)?;
        }
        Ok(out)
    }

    /// Reads the `s` smallest mixed keys starting from bucket 0.
    pub fn sketch<K: KmerKey>(&self, s: usize) -> Result<Vec<K>> {
        let cap = s.min(
            usize::try_from(self.n_kmers)
                .unwrap_or(MAX_EAGER_CAP)
                .min(MAX_EAGER_CAP),
        );
        let mut out = Vec::new();
        out.try_reserve(cap)
            .context("sketch buffer allocation failed")?;
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

/// Random-access reader for verified counted `.ktab` (`GRMKTB1\0`) tables.
pub struct KtabReader {
    file: File,
    file_len: u64,
    pub k: usize,
    pub key_width: u8,
    pub modality: InputModality,
    pub presence_only: bool,
    pub min_count: u32,
    pub max_count_observed: u32,
    pub saturated_counts: u64,
    pub name: String,
    pub n_kmers: u64,
    pub digest_sha256: [u8; 32],
    key_offsets: Vec<u64>,
    count_offsets: Vec<u64>,
    bucket_counts: Vec<u32>,
    key_data_start: u64,
    count_data_start: u64,
}

impl KtabReader {
    /// Opens a `.ktab` file and validates all header fields and bucket table
    /// bounds against the file length without scanning the payload body.
    pub fn open(path: &Path) -> Result<Self> {
        let mut file = File::open(path).with_context(|| format!("open {}", path.display()))?;
        let file_len = file.metadata()?.len();
        ensure!(
            file_len >= KTAB_HEADER_LEN as u64,
            "{} is truncated before GRMKTB1 header",
            path.display()
        );
        let mut head = [0u8; KTAB_HEADER_LEN];
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
        let modality = InputModality::from_byte(head[12])
            .with_context(|| format!("{}: invalid modality", path.display()))?;
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
        ensure!(
            (flags & !VALID_KTAB_FLAGS) == 0,
            "{}: unrecognized GRMKTB1 flags {flags:#04x}",
            path.display()
        );
        ensure!(
            head[13..16] == [0u8; 3] && head[44..48] == [0u8; 4],
            "{}: non-zero reserved header bytes in GRMKTB1",
            path.display()
        );
        let presence_only = (flags & FLAG_PRESENCE_ONLY) != 0;
        let saturated_flag = (flags & FLAG_SATURATED_COUNTS) != 0;
        let min_count = u32::from_le_bytes(head[16..20].try_into()?);
        let max_count_observed = u32::from_le_bytes(head[20..24].try_into()?);
        let n_kmers = u64::from_le_bytes(head[24..32].try_into()?);
        let saturated_counts = u64::from_le_bytes(head[32..40].try_into()?);
        let name_len_u32 = u32::from_le_bytes(head[40..44].try_into()?);
        let mut digest_sha256 = [0u8; 32];
        digest_sha256.copy_from_slice(&head[KTAB_DIGEST_START..KTAB_DIGEST_END]);

        ensure!(
            saturated_flag == (saturated_counts > 0),
            "{}: saturated_counts flag disagrees with saturated_counts={saturated_counts}",
            path.display()
        );
        ensure!(
            saturated_counts <= n_kmers,
            "{}: saturated_counts ({saturated_counts}) exceeds n_kmers ({n_kmers})",
            path.display()
        );
        if presence_only {
            ensure!(
                saturated_counts == 0 && max_count_observed <= 1,
                "{}: presence-only ktab has invalid count metadata",
                path.display()
            );
        }
        if saturated_counts > 0 {
            ensure!(
                max_count_observed == MAX_COUNT,
                "{}: saturated_counts > 0 requires max_count_observed == u32::MAX",
                path.display()
            );
        }
        if n_kmers == 0 {
            ensure!(
                max_count_observed == 0,
                "{}: empty ktab has non-zero max_count_observed",
                path.display()
            );
        }

        let table_len = 2 * (BUCKETS + 1) * 8 + BUCKETS * 4;
        let key_data_start = (KTAB_HEADER_LEN as u64)
            .checked_add(name_len_u32 as u64)
            .and_then(|v| v.checked_add(table_len as u64))
            .context("header length overflow")?;
        ensure!(
            key_data_start <= file_len,
            "{} is truncated",
            path.display()
        );
        let name_len = usize::try_from(name_len_u32).context("name length overflow")?;
        let mut name_bytes = alloc_zeroed_checked(name_len, "GRMKTB1 name")?;
        file.read_exact(&mut name_bytes)?;
        let name = String::from_utf8(name_bytes)
            .with_context(|| format!("{}: sample name is not valid UTF-8", path.display()))?;
        let mut table = alloc_zeroed_checked(table_len, "GRMKTB1 bucket table")?;
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
                && (!presence_only || count_data_len == 0)
                && count_data_start.checked_add(count_data_len) == Some(file_len),
            "{}: inconsistent GRMKTB1 section lengths",
            path.display()
        );

        let mut sum_counts = 0u64;
        for b in 0..BUCKETS {
            let klen = key_offsets[b + 1]
                .checked_sub(key_offsets[b])
                .with_context(|| format!("{}: non-monotonic key offsets", path.display()))?;
            let clen = count_offsets[b + 1]
                .checked_sub(count_offsets[b])
                .with_context(|| format!("{}: non-monotonic count offsets", path.display()))?;
            let c = bucket_counts[b] as u64;
            if c == 0 {
                ensure!(
                    klen == 0 && clen == 0,
                    "{}: empty bucket {b} has non-zero payload length",
                    path.display()
                );
            } else {
                let min_klen = 1u64
                    .checked_add(c.div_ceil(8))
                    .context("bucket minimum key payload overflow")?;
                ensure!(
                    klen >= min_klen,
                    "{}: bucket {b} key payload ({klen} bytes) too short for {c} keys",
                    path.display()
                );
                if presence_only {
                    ensure!(
                        clen == 0,
                        "{}: presence-only bucket {b} has count bytes",
                        path.display()
                    );
                } else {
                    ensure!(
                        clen >= 1,
                        "{}: non-empty counted bucket {b} has empty count payload",
                        path.display()
                    );
                }
            }
            sum_counts = sum_counts
                .checked_add(c)
                .context("bucket count sum overflow")?;
        }
        ensure!(
            sum_counts == n_kmers,
            "{}: bucket counts sum ({sum_counts}) != header n_kmers ({n_kmers})",
            path.display()
        );

        Ok(Self {
            file,
            file_len,
            k,
            key_width,
            modality,
            presence_only,
            min_count,
            max_count_observed,
            saturated_counts,
            name,
            n_kmers,
            digest_sha256,
            key_offsets,
            count_offsets,
            bucket_counts,
            key_data_start,
            count_data_start,
        })
    }

    /// Opens a `.ktab` file and immediately verifies its full SHA-256 digest.
    pub fn open_verified(path: &Path) -> Result<Self> {
        let reader = Self::open(path)?;
        reader.verify_integrity()?;
        Ok(reader)
    }

    /// Recomputes the SHA-256 digest over header metadata, bucket tables, and
    /// key/count payloads and ensures it matches the digest stored in the header.
    pub fn verify_integrity(&self) -> Result<[u8; 32]> {
        let got = compute_file_sha256_excluding_slot(
            &self.file,
            self.file_len,
            KTAB_DIGEST_START,
            KTAB_DIGEST_END,
        )?;
        ensure!(
            got == self.digest_sha256,
            "GRMKTB1 SHA-256 integrity check failed: header has {}, computed {}",
            hex_bytes(&self.digest_sha256),
            hex_bytes(&got)
        );
        Ok(got)
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
        ensure!(bucket < BUCKETS, "bucket index {bucket} out of range");
        ensure!(
            K::WIDTH_BYTES >= self.key_width as usize,
            "cannot decode k={} into {}-byte key",
            self.k,
            K::WIDTH_BYTES
        );
        let (start, end) = (self.key_offsets[bucket], self.key_offsets[bucket + 1]);
        let len = usize::try_from(end - start).context("key bucket length overflow")?;
        buf.clear();
        buf.try_reserve(len)
            .with_context(|| format!("bucket {bucket}: key payload allocation failed"))?;
        buf.resize(len, 0);
        if !buf.is_empty() {
            let offset = self
                .key_data_start
                .checked_add(start)
                .context("key bucket file offset overflow")?;
            self.file.read_exact_at(buf, offset)?;
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
        let c_len = usize::try_from(c_end - c_start).context("count bucket length overflow")?;
        count_buf.clear();
        count_buf
            .try_reserve(c_len)
            .with_context(|| format!("bucket {bucket}: count payload allocation failed"))?;
        count_buf.resize(c_len, 0);
        if !count_buf.is_empty() {
            let offset = self
                .count_data_start
                .checked_add(c_start)
                .context("count bucket file offset overflow")?;
            self.file.read_exact_at(count_buf, offset)?;
        }
        decode_bucket_counts(
            bucket,
            self.bucket_counts[bucket],
            self.presence_only,
            count_buf,
            out_counts,
        )
    }

    /// Verifies the stored SHA-256 digest and decodes every bucket in order.
    pub fn read_all<K: KmerKey>(&self) -> Result<Vec<(K, u32)>> {
        self.verify_integrity()?;
        let cap = usize::try_from(self.n_kmers)
            .context("n_kmers exceeds addressable memory")?
            .min(MAX_EAGER_CAP);
        let mut keys = Vec::new();
        keys.try_reserve(cap)
            .context("ktab key list allocation failed")?;
        let mut counts = Vec::new();
        counts
            .try_reserve(cap)
            .context("ktab count list allocation failed")?;
        let mut kbuf = Vec::new();
        let mut cbuf = Vec::new();
        for b in 0..BUCKETS {
            self.read_bucket(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
        }
        let observed_max = counts.iter().copied().max().unwrap_or(0);
        ensure!(
            observed_max == self.max_count_observed,
            "ktab header max_count_observed ({}) != decoded max ({observed_max})",
            self.max_count_observed
        );
        Ok(keys.into_iter().zip(counts).collect())
    }

    pub fn sketch<K: KmerKey>(&self, s: usize) -> Result<Vec<K>> {
        let cap = s.min(
            usize::try_from(self.n_kmers)
                .unwrap_or(MAX_EAGER_CAP)
                .min(MAX_EAGER_CAP),
        );
        let mut out = Vec::new();
        out.try_reserve(cap)
            .context("sketch buffer allocation failed")?;
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

    pub fn modality(&self) -> Option<InputModality> {
        match self {
            Self::V1(_) => None,
            Self::Gr(r) => Some(r.modality),
            Self::Tab(r) => Some(r.modality),
        }
    }

    pub fn digest_sha256(&self) -> Option<[u8; 32]> {
        match self {
            Self::V1(_) => None,
            Self::Gr(r) => Some(r.digest_sha256),
            Self::Tab(r) => Some(r.digest_sha256),
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

    /// Verifies the stored SHA-256 header digest when the format carries one.
    /// Returns `Ok(true)` if a stored digest was recomputed and matched,
    /// `Ok(false)` for `.kset` v1 (which stores no header digest), and `Err`
    /// on any integrity mismatch.
    pub fn verify_integrity(&self) -> Result<bool> {
        match self {
            Self::V1(_) => Ok(false),
            Self::Gr(r) => {
                r.verify_integrity()?;
                Ok(true)
            }
            Self::Tab(r) => {
                r.verify_integrity()?;
                Ok(true)
            }
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
                let mut tmp = Vec::new();
                tmp.try_reserve(r.bucket_len(bucket))
                    .context("kset bucket buffer allocation failed")?;
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

/// Summary returned by table format conversions.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConvertOutcome {
    pub format: TableFormat,
    pub k: usize,
    pub sample: String,
    pub modality: Option<InputModality>,
    pub n_kmers: u64,
    pub counts_dropped: bool,
}

/// Converts a `.kset` v1 (`GRMKSET1`), `.grk`, or `.ktab` file to a verified
/// Golomb-Rice sorted-key table (`GRMKGR1\0`), checking the `.grk` output
/// extension and verifying any input SHA-256 digest before converting.
pub fn convert_to_gr(
    src: &Path,
    dst: &Path,
    name: &str,
    modality: Option<InputModality>,
) -> Result<ConvertOutcome> {
    check_output_extension(dst, TableFormat::GrKset)?;
    let reader = AnyTableReader::open(src)?;
    reader.verify_integrity()?;
    let k = reader.k();
    let out_modality = modality.unwrap_or_else(|| reader.modality().unwrap_or_default());
    let mut bucket_counts = [0u32; BUCKETS];
    let mut payloads = vec![Vec::new(); BUCKETS];
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut counts = Vec::new();
    let mut counts_dropped = false;

    if k <= MAX_K_U64 {
        let mut keys: Vec<u64> = Vec::new();
        for b in 0..BUCKETS {
            keys.clear();
            counts.clear();
            reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if counts.iter().any(|&c| c > 1) {
                counts_dropped = true;
            }
            bucket_counts[b] = u32::try_from(keys.len()).context("bucket count exceeds u32")?;
            encode_bucket_keys(k, b, &keys, &mut payloads[b]);
        }
    } else {
        let mut keys: Vec<u128> = Vec::new();
        for b in 0..BUCKETS {
            keys.clear();
            counts.clear();
            reader.read_bucket::<u128>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if counts.iter().any(|&c| c > 1) {
                counts_dropped = true;
            }
            bucket_counts[b] = u32::try_from(keys.len()).context("bucket count exceeds u32")?;
            encode_bucket_keys(k, b, &keys, &mut payloads[b]);
        }
    }
    write_gr_kset_from_buckets(dst, name, k, out_modality, &bucket_counts, &payloads)?;
    Ok(ConvertOutcome {
        format: TableFormat::GrKset,
        k,
        sample: name.to_string(),
        modality: Some(out_modality),
        n_kmers: reader.n_kmers(),
        counts_dropped,
    })
}

/// Converts a `.kset` v1 (`GRMKSET1`), `.grk`, or `.ktab` file to a verified
/// `.ktab` (`GRMKTB1\0`), checking the `.ktab` output extension and verifying
/// any input SHA-256 digest before converting.
pub fn convert_to_ktab(
    src: &Path,
    dst: &Path,
    name: &str,
    modality: Option<InputModality>,
) -> Result<ConvertOutcome> {
    check_output_extension(dst, TableFormat::Ktab)?;
    let reader = AnyTableReader::open(src)?;
    reader.verify_integrity()?;
    let k = reader.k();
    let out_modality = modality.unwrap_or_else(|| reader.modality().unwrap_or_default());
    let (min_count, saturated_counts, presence_only) = match &reader {
        AnyTableReader::V1(_) | AnyTableReader::Gr(_) => (1, 0, true),
        AnyTableReader::Tab(r) => (r.min_count, r.saturated_counts, r.presence_only),
    };

    let mut bucket_counts = [0u32; BUCKETS];
    let mut key_payloads = vec![Vec::new(); BUCKETS];
    let mut count_payloads = vec![Vec::new(); BUCKETS];
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut counts = Vec::new();
    let mut max_count_observed = 0u32;

    if k <= MAX_K_U64 {
        let mut keys: Vec<u64> = Vec::new();
        for b in 0..BUCKETS {
            keys.clear();
            counts.clear();
            reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if let Some(&m) = counts.iter().max() {
                max_count_observed = max_count_observed.max(m);
            }
            bucket_counts[b] = u32::try_from(keys.len()).context("bucket count exceeds u32")?;
            encode_bucket_keys(k, b, &keys, &mut key_payloads[b]);
            if !presence_only {
                encode_bucket_counts(&counts, &mut count_payloads[b]);
            }
        }
    } else {
        let mut keys: Vec<u128> = Vec::new();
        for b in 0..BUCKETS {
            keys.clear();
            counts.clear();
            reader.read_bucket::<u128>(b, &mut kbuf, &mut cbuf, &mut keys, &mut counts)?;
            if let Some(&m) = counts.iter().max() {
                max_count_observed = max_count_observed.max(m);
            }
            bucket_counts[b] = u32::try_from(keys.len()).context("bucket count exceeds u32")?;
            encode_bucket_keys(k, b, &keys, &mut key_payloads[b]);
            if !presence_only {
                encode_bucket_counts(&counts, &mut count_payloads[b]);
            }
        }
    }

    let meta = KtabWriteMeta {
        name,
        k,
        modality: out_modality,
        min_count,
        max_count_observed,
        saturated_counts,
        presence_only,
    };
    write_ktab_from_buckets(dst, &meta, &bucket_counts, &key_payloads, &count_payloads)?;
    Ok(ConvertOutcome {
        format: TableFormat::Ktab,
        k,
        sample: name.to_string(),
        modality: Some(out_modality),
        n_kmers: reader.n_kmers(),
        counts_dropped: false,
    })
}

/// Converts a `.kset` v1 (`GRMKSET1`), `.grk`, or `.ktab` table with `k <= 32`
/// to `.kset` v1 (`GRMKSET1`), checking the `.kset` output extension and
/// verifying any input SHA-256 digest before converting.
pub fn convert_to_v1(src: &Path, dst: &Path, name: &str) -> Result<ConvertOutcome> {
    check_output_extension(dst, TableFormat::KsetV1)?;
    let reader = AnyTableReader::open(src)?;
    reader.verify_integrity()?;
    let k = reader.k();
    ensure!(
        k <= MAX_K_U64,
        "cannot convert k={k} to GRMKSET1 (max k is {MAX_K_U64})"
    );
    let space = KmerSpace::new(k);
    let cap = usize::try_from(reader.n_kmers())
        .context("n_kmers exceeds addressable memory")?
        .min(MAX_EAGER_CAP);
    let mut all = Vec::new();
    all.try_reserve(cap)
        .context("kset output buffer allocation failed")?;
    let mut kbuf = Vec::new();
    let mut cbuf = Vec::new();
    let mut counts = Vec::new();
    let mut counts_dropped = false;
    for b in 0..BUCKETS {
        counts.clear();
        reader.read_bucket::<u64>(b, &mut kbuf, &mut cbuf, &mut all, &mut counts)?;
        if counts.iter().any(|&c| c > 1) {
            counts_dropped = true;
        }
    }
    kset::write(dst, name, &space, &all)?;
    Ok(ConvertOutcome {
        format: TableFormat::KsetV1,
        k,
        sample: name.to_string(),
        modality: None,
        n_kmers: reader.n_kmers(),
        counts_dropped,
    })
}

// ---------------------------------------------------------------------------
// Compact per-sequence k-mer count profile (.prof.bin, GRMPRF1\0)
// ---------------------------------------------------------------------------

/// One sequence's positional k-mer count profile.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SequenceProfile {
    pub name: String,
    pub seq_len: u64,
    pub counts: Vec<u32>,
}

/// Header metadata for a verified `.prof.bin` (`GRMPRF1\0`) file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ProfileHeader {
    pub k: usize,
    pub max_count_observed: u32,
    pub n_sequences: u64,
    pub total_bases: u64,
    pub digest_sha256: [u8; 32],
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
        let low = (b & 0x7f) as u64;
        if shift == 63 {
            ensure!(b <= 1, "varint overflow");
        }
        v |= low << shift;
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

/// Decodes a run-length + zigzag-delta ULEB128 profile stream into `n_positions`
/// counts, checking overflow and growing fallibly so crafted lengths never panic
/// or abort.
pub fn decode_profile_runs(bytes: &[u8], n_positions: usize) -> Result<Vec<u32>> {
    if n_positions == 0 {
        ensure!(
            bytes.is_empty(),
            "empty profile sequence carries trailing payload bytes"
        );
        return Ok(Vec::new());
    }
    ensure!(
        !bytes.is_empty(),
        "non-empty profile sequence has empty payload"
    );
    let init_cap = n_positions.min(65_536);
    let mut out = Vec::new();
    out.try_reserve(init_cap)
        .context("profile count buffer allocation failed")?;
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
        let run_usize = usize::try_from(run).context("profile run exceeds addressable memory")?;
        let new_len = out
            .len()
            .checked_add(run_usize)
            .context("profile length overflow")?;
        ensure!(
            new_len <= n_positions,
            "profile decoded more positions than header declared"
        );
        out.try_reserve(run_usize)
            .context("profile count buffer expansion failed")?;
        out.resize(new_len, val);
        prev = val as i64;
    }
    ensure!(
        out.len() == n_positions,
        "profile decoded {} positions, expected {n_positions}",
        out.len()
    );
    Ok(out)
}

/// Writes per-sequence positional k-mer count profiles into a verified
/// run-length + zigzag-delta binary file (`GRMPRF1\0`).
pub fn write_profiles(path: &Path, k: usize, profiles: &[SequenceProfile]) -> Result<()> {
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported k {k}");
    let mut total_bases = 0u64;
    let mut max_count_observed = 0u32;
    let mut body = Vec::new();
    let mut enc = Vec::new();
    for p in profiles {
        let name_len = u32::try_from(p.name.len()).context("profile sequence name exceeds u32")?;
        ensure!(
            (p.counts.len() as u64) <= p.seq_len,
            "profile counts length ({}) exceeds seq_len ({})",
            p.counts.len(),
            p.seq_len
        );
        total_bases = total_bases
            .checked_add(p.seq_len)
            .context("profile total_bases overflow")?;
        if let Some(&m) = p.counts.iter().max() {
            max_count_observed = max_count_observed.max(m);
        }
        encode_profile_runs(&p.counts, &mut enc);
        body.extend_from_slice(&name_len.to_le_bytes());
        body.extend_from_slice(p.name.as_bytes());
        body.extend_from_slice(&p.seq_len.to_le_bytes());
        body.extend_from_slice(&(p.counts.len() as u64).to_le_bytes());
        body.extend_from_slice(&(enc.len() as u64).to_le_bytes());
        body.extend_from_slice(&enc);
    }

    // Fixed header prefix [0..40]; bytes [40..72] hold SHA-256([0..40] || body).
    let mut prefix = [0u8; PROF_DIGEST_START];
    prefix[0..8].copy_from_slice(MAGIC_PROF1);
    prefix[8] = k as u8;
    // [9..12] reserved0 = 0
    prefix[12..16].copy_from_slice(&max_count_observed.to_le_bytes());
    prefix[16..24].copy_from_slice(&(profiles.len() as u64).to_le_bytes());
    prefix[24..32].copy_from_slice(&total_bases.to_le_bytes());
    // [32..40] reserved1 = 0

    let mut sha = Sha256Context::new(&SHA256);
    sha.update(&prefix);
    sha.update(&body);
    let digest = sha.finish();

    let part = confine::suffixed(path, ".part");
    confine::regular_or_absent(path)?;
    let mut w = BufWriter::new(confine::create_part(&part)?);
    let written = (|| -> Result<()> {
        w.write_all(&prefix)?;
        w.write_all(digest.as_ref())?;
        w.write_all(&body)?;
        w.flush()?;
        confine::replace_with(&part, path)
    })();
    drop(w);
    if written.is_err() && confine::regular_or_absent(&part).unwrap_or(false) {
        let _ = std::fs::remove_file(&part);
    }
    written
}

fn parse_profile_header(path: &Path, file: &File, file_len: u64) -> Result<ProfileHeader> {
    ensure!(
        file_len >= PROF_HEADER_LEN as u64,
        "{} is truncated before GRMPRF1 header",
        path.display()
    );
    let mut head = [0u8; PROF_HEADER_LEN];
    file.read_exact_at(&mut head, 0)
        .with_context(|| format!("read profile header {}", path.display()))?;
    ensure!(
        &head[..8] == MAGIC_PROF1,
        "{} is not a GRMPRF1 profile file",
        path.display()
    );
    let k = head[8] as usize;
    ensure!((MIN_K..=MAX_K).contains(&k), "unsupported profile k {k}");
    ensure!(
        head[9..12] == [0u8; 3] && head[32..40] == [0u8; 8],
        "{}: non-zero reserved header bytes in GRMPRF1",
        path.display()
    );
    let max_count_observed = u32::from_le_bytes(head[12..16].try_into()?);
    let n_sequences = u64::from_le_bytes(head[16..24].try_into()?);
    let total_bases = u64::from_le_bytes(head[24..32].try_into()?);
    let mut digest_sha256 = [0u8; 32];
    digest_sha256.copy_from_slice(&head[PROF_DIGEST_START..PROF_DIGEST_END]);

    let body_bytes = file_len - PROF_HEADER_LEN as u64;
    let max_possible_seqs = body_bytes / PROF_MIN_RECORD_BYTES;
    ensure!(
        n_sequences <= max_possible_seqs,
        "{}: header declares {n_sequences} sequences but only {body_bytes} payload bytes remain",
        path.display()
    );
    if n_sequences == 0 {
        ensure!(
            body_bytes == 0 && total_bases == 0 && max_count_observed == 0,
            "{}: empty profile file has inconsistent header or trailing bytes",
            path.display()
        );
    }

    Ok(ProfileHeader {
        k,
        max_count_observed,
        n_sequences,
        total_bases,
        digest_sha256,
    })
}

/// Reads and structurally validates a `.prof.bin` (`GRMPRF1\0`) header without
/// scanning or verifying the payload body.
pub fn read_profile_header(path: &Path) -> Result<ProfileHeader> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let file_len = file.metadata()?.len();
    parse_profile_header(path, &file, file_len)
}

/// Verifies the stored 32-byte SHA-256 digest of a `.prof.bin` (`GRMPRF1\0`)
/// file and decodes every sequence profile with remaining-bytes bounds checks
/// before every allocation.
pub fn read_profiles_verified(path: &Path) -> Result<(ProfileHeader, Vec<SequenceProfile>)> {
    let file = File::open(path).with_context(|| format!("open {}", path.display()))?;
    let file_len = file.metadata()?.len();
    let hdr = parse_profile_header(path, &file, file_len)?;

    let computed =
        compute_file_sha256_excluding_slot(&file, file_len, PROF_DIGEST_START, PROF_DIGEST_END)?;
    ensure!(
        computed == hdr.digest_sha256,
        "{}: GRMPRF1 SHA-256 integrity check failed: header has {}, computed {}",
        path.display(),
        hex_bytes(&hdr.digest_sha256),
        hex_bytes(&computed)
    );

    let n_seq = usize::try_from(hdr.n_sequences).context("profile sequence count overflow")?;
    let mut profiles = Vec::new();
    profiles
        .try_reserve(n_seq.min(4096))
        .context("profile list allocation failed")?;

    let mut pos = PROF_HEADER_LEN as u64;
    let mut remaining = file_len - pos;
    let mut sum_bases = 0u64;
    let mut obs_max = 0u32;

    for idx in 0..n_seq {
        ensure!(
            remaining >= PROF_MIN_RECORD_BYTES,
            "{}: truncated before sequence record {idx}",
            path.display()
        );
        let mut len_buf = [0u8; 4];
        file.read_exact_at(&mut len_buf, pos)?;
        pos += 4;
        remaining -= 4;
        let name_len_u32 = u32::from_le_bytes(len_buf);
        let name_len_u64 = name_len_u32 as u64;
        ensure!(
            name_len_u64 <= remaining - 24,
            "{}: sequence {idx} name length {name_len_u64} exceeds remaining bytes",
            path.display()
        );
        let name_len = usize::try_from(name_len_u32).context("profile name length overflow")?;
        let mut name_bytes = alloc_zeroed_checked(name_len, "GRMPRF1 sequence name")?;
        if !name_bytes.is_empty() {
            file.read_exact_at(&mut name_bytes, pos)?;
        }
        pos += name_len_u64;
        remaining -= name_len_u64;
        let name = String::from_utf8(name_bytes).with_context(|| {
            format!("{}: sequence {idx} name is not valid UTF-8", path.display())
        })?;

        let mut meta = [0u8; 24];
        file.read_exact_at(&mut meta, pos)?;
        pos += 24;
        remaining -= 24;
        let seq_len = u64::from_le_bytes(meta[0..8].try_into()?);
        let n_pos_u64 = u64::from_le_bytes(meta[8..16].try_into()?);
        let payload_len_u64 = u64::from_le_bytes(meta[16..24].try_into()?);
        ensure!(
            payload_len_u64 <= remaining,
            "{}: sequence {idx} payload length {payload_len_u64} exceeds remaining {remaining} bytes",
            path.display()
        );
        ensure!(
            n_pos_u64 <= seq_len,
            "{}: sequence {idx} n_positions ({n_pos_u64}) exceeds seq_len ({seq_len})",
            path.display()
        );
        let n_pos = usize::try_from(n_pos_u64).context("profile n_positions overflow")?;
        let payload_len =
            usize::try_from(payload_len_u64).context("profile payload_len overflow")?;
        let mut payload = alloc_zeroed_checked(payload_len, "GRMPRF1 sequence payload")?;
        if !payload.is_empty() {
            file.read_exact_at(&mut payload, pos)?;
        }
        pos += payload_len_u64;
        remaining -= payload_len_u64;

        let counts = decode_profile_runs(&payload, n_pos)?;
        sum_bases = sum_bases
            .checked_add(seq_len)
            .context("profile total_bases overflow")?;
        if let Some(&m) = counts.iter().max() {
            obs_max = obs_max.max(m);
        }
        profiles.push(SequenceProfile {
            name,
            seq_len,
            counts,
        });
    }

    ensure!(
        remaining == 0,
        "{}: {remaining} trailing bytes after last sequence profile",
        path.display()
    );
    ensure!(
        sum_bases == hdr.total_bases,
        "{}: decoded total_bases ({sum_bases}) != header ({})",
        path.display(),
        hdr.total_bases
    );
    ensure!(
        obs_max == hdr.max_count_observed,
        "{}: decoded max_count_observed ({obs_max}) != header ({})",
        path.display(),
        hdr.max_count_observed
    );
    Ok((hdr, profiles))
}

/// Reads and verifies all sequence profiles from a `.prof.bin` (`GRMPRF1\0`) file.
pub fn read_profiles(path: &Path) -> Result<(usize, Vec<SequenceProfile>)> {
    let (hdr, profiles) = read_profiles_verified(path)?;
    Ok((hdr.k, profiles))
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
            write_gr_kset(&p, DEFAULT_SAMPLE_NAME, k, InputModality::Assembly, &vals).unwrap();
            let r = GrKsetReader::open_verified(&p).unwrap();
            assert_eq!(r.k, k);
            assert_eq!(r.modality, InputModality::Assembly);
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
            write_gr_kset(&p, "wide", k, InputModality::Reads, &vals).unwrap();
            let r = GrKsetReader::open_verified(&p).unwrap();
            assert_eq!(r.k, k);
            assert_eq!(r.key_width, 16);
            assert_eq!(r.modality, InputModality::Reads);
            assert_eq!(r.read_all::<u128>().unwrap(), vals);
            assert_eq!(r.sketch::<u128>(50).unwrap(), vals[..50]);
        }
    }

    #[test]
    fn ktab_preserves_u32_counts_saturation_and_8208_byte_presence_overhead() {
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
            (120, 65_536),
            (130, 70_000),
            (140, 1_500_000),
            (150, u32::MAX),
        ];
        let p = dir.path().join("counts.ktab");
        write_ktab(&p, "c", 31, InputModality::Reads, 1, 1, &entries).unwrap();
        let r = KtabReader::open_verified(&p).unwrap();
        assert!(!r.presence_only);
        assert_eq!(r.modality, InputModality::Reads);
        assert_eq!(r.max_count_observed, u32::MAX);
        assert_eq!(r.saturated_counts, 1);
        assert_eq!(r.read_all::<u64>().unwrap(), entries);

        let ones: Vec<(u64, u32)> = (0..2_000u64).map(|i| (i * 17, 1)).collect();
        let keys_only: Vec<u64> = ones.iter().map(|&(k, _)| k).collect();
        let p_ones = dir.path().join("ones.ktab");
        let p_gr = dir.path().join("ones.grk");
        write_ktab(&p_ones, "o", 31, InputModality::Assembly, 1, 0, &ones).unwrap();
        write_gr_kset(&p_gr, "o", 31, InputModality::Assembly, &keys_only).unwrap();
        let ro = KtabReader::open_verified(&p_ones).unwrap();
        assert!(ro.presence_only);
        assert_eq!(ro.saturated_counts, 0);
        assert_eq!(ro.read_all::<u64>().unwrap(), ones);
        let ktab_len = std::fs::metadata(&p_ones).unwrap().len();
        let grk_len = std::fs::metadata(&p_gr).unwrap().len();
        assert_eq!(ktab_len - grk_len, 8_208);
    }

    #[test]
    fn profile_round_trips_runs_and_deltas() {
        let dir = tempfile::tempdir().unwrap();
        let profs = vec![
            SequenceProfile {
                name: "contig1".to_string(),
                seq_len: 100,
                counts: vec![1, 1, 1, 1, 0, 0, 5, 5, 5, 300, 65_535, 65_536, u32::MAX, 1],
            },
            SequenceProfile {
                name: "short".to_string(),
                seq_len: 10,
                counts: vec![],
            },
        ];
        let p = dir.path().join("test.prof.bin");
        write_profiles(&p, 31, &profs).unwrap();
        let hdr = read_profile_header(&p).unwrap();
        assert_eq!(hdr.k, 31);
        assert_eq!(hdr.n_sequences, 2);
        assert_eq!(hdr.total_bases, 110);
        assert_eq!(hdr.max_count_observed, u32::MAX);
        let (k, got) = read_profiles(&p).unwrap();
        assert_eq!(k, 31);
        assert_eq!(got, profs);
    }
}
