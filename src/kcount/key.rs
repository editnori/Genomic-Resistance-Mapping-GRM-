//! Dual-width canonical k-mer representation (`k = 5..=64`), branchless O(1)
//! SWAR reverse complement, and invertible 2k-bit mixers for `u64` and `u128`.
//!
//! For `k <= 32`, [`KmerSpace64`] uses the exact same 2-bit encoding and
//! bijective mixer constants as [`crate::dna::KmerSpace`], preserving
//! bit-for-bit parity with `.kset` v1 digests and sketches. For `k = 33..=64`,
//! [`KmerSpace128`] packs `2k` bits into a `u128` with a 3-round invertible
//! 128-bit mixer so top-bit bucketing and bottom-`s` MinHash sketching remain
//! balanced and invertible across the full `2k`-bit space.

use crate::dna::{self, CODE};

/// Minimum supported k-mer length (needs at least 10 bits for 1,024 buckets).
pub const MIN_K: usize = dna::MIN_K;
/// Maximum k-mer length that fits in the fast `u64` path.
pub const MAX_K_U64: usize = dna::MAX_K;
/// Maximum supported k-mer length across both `u64` and `u128` paths.
pub const MAX_K: usize = 64;

/// First public-domain SplitMix64 finalizer multiplier from Steele, Lea, and
/// Flood (OOPSLA 2014, public domain), reused as the high 64 bits of round 0's
/// 128-bit mixer multiplier.
const SPLITMIX64_MUL_1: u64 = 0xbf58_476d_1ce4_e5b9;
/// Second public-domain SplitMix64 finalizer multiplier from Steele, Lea, and
/// Flood (OOPSLA 2014, public domain), reused as the low 64 bits of round 0's
/// 128-bit mixer multiplier.
const SPLITMIX64_MUL_2: u64 = 0x94d0_49bb_1331_11eb;

/// Odd 64-bit multiplier halves for `KmerSpace128` rounds 1 and 2,
/// deterministically derived from repository seeds
/// `SHA-256("GRM-KCOUNT-MIX128-V1:0")` through `":3"`, taking the first
/// 8 digest bytes as little-endian `u64` with the low bit forced to 1 (`| 1`).
/// Verified by `mix128_constants_match_sha256_seeds`.
#[cfg(test)]
const MIX128_SEED_PREFIX: &str = "GRM-KCOUNT-MIX128-V1:";
const MIX128_HALF_0: u64 = 0x1059_c9c3_f3a8_0405;
const MIX128_HALF_1: u64 = 0xaabe_85ea_b824_763f;
const MIX128_HALF_2: u64 = 0x2e61_9c2c_9ca4_828b;
const MIX128_HALF_3: u64 = 0x26d4_aafe_538a_10ef;

const MIX128_MUL_0: u128 = ((SPLITMIX64_MUL_1 as u128) << 64) | (SPLITMIX64_MUL_2 as u128);
const MIX128_MUL_1: u128 = ((MIX128_HALF_0 as u128) << 64) | (MIX128_HALF_1 as u128);
const MIX128_MUL_2: u128 = ((MIX128_HALF_2 as u128) << 64) | (MIX128_HALF_3 as u128);

/// Branchless O(1) SWAR reverse complement of a 2-bit packed `u64` k-mer.
#[inline]
pub fn reverse_complement_u64(mut x: u64, k: usize) -> u64 {
    debug_assert!((1..=MAX_K_U64).contains(&k));
    x = ((x & 0x3333_3333_3333_3333) << 2) | ((x >> 2) & 0x3333_3333_3333_3333);
    x = ((x & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F);
    x = (!x).swap_bytes();
    x >> (64 - 2 * k as u32)
}

/// Branchless O(1) SWAR reverse complement of a 2-bit packed `u128` k-mer.
#[inline]
pub fn reverse_complement_u128(mut x: u128, k: usize) -> u128 {
    debug_assert!((1..=MAX_K).contains(&k));
    x = ((x & 0x3333_3333_3333_3333_3333_3333_3333_3333) << 2)
        | ((x >> 2) & 0x3333_3333_3333_3333_3333_3333_3333_3333);
    x = ((x & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F) << 4)
        | ((x >> 4) & 0x0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F_0F0F);
    x = (!x).swap_bytes();
    x >> (128 - 2 * k as u32)
}

/// Fast `u64` k-mer space for `k = 5..=32`, bit-for-bit identical to
/// [`crate::dna::KmerSpace`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KmerSpace64 {
    pub k: usize,
    inner: dna::KmerSpace,
    mask: u64,
}

impl KmerSpace64 {
    pub fn new(k: usize) -> Self {
        assert!(
            (MIN_K..=MAX_K_U64).contains(&k),
            "k for u64 path must be within {MIN_K}..={MAX_K_U64}, got {k}"
        );
        let inner = dna::KmerSpace::new(k);
        let mask = inner.mask();
        Self { k, inner, mask }
    }

    #[inline]
    pub fn bits(&self) -> u32 {
        self.inner.bits()
    }

    #[inline]
    pub fn mask(&self) -> u64 {
        self.mask
    }

    #[inline]
    pub fn mix(&self, x: u64) -> u64 {
        self.inner.mix(x)
    }

    #[inline]
    pub fn unmix(&self, y: u64) -> u64 {
        self.inner.unmix(y)
    }

    #[inline]
    pub fn reverse_complement(&self, kmer: u64) -> u64 {
        reverse_complement_u64(kmer, self.k)
    }

    pub fn decode(&self, kmer: u64) -> String {
        self.inner.decode(kmer)
    }

    pub fn canonical(&self, text: &str) -> Option<u64> {
        self.inner.canonical(text)
    }

    /// Calls `f` with every canonical packed k-mer in `seq`. Any non-ACGT byte
    /// resets the rolling window cleanly.
    #[inline]
    pub fn for_each_canonical(&self, seq: &[u8], f: impl FnMut(u64)) {
        self.inner.for_each_canonical(seq, f);
    }

    /// Unpacks canonical mixed k-mers directly from a 2-bit packed super-k-mer
    /// (4 bases per byte, LSB-first).
    #[inline]
    pub fn for_each_mixed_in_packed_2bit(&self, packed: &[u8], len: usize, mut f: impl FnMut(u64)) {
        let k = self.k;
        if len < k {
            return;
        }
        assert!(
            packed.len() >= len.div_ceil(4),
            "packed 2-bit buffer too short for {len} bases"
        );
        let shift = 2 * (k - 1);
        let (mut fwd, mut rev) = (0u64, 0u64);
        for i in 0..len {
            let c = ((packed[i >> 2] >> ((i & 3) * 2)) & 3) as u64;
            fwd = ((fwd << 2) | c) & self.mask;
            rev = (rev >> 2) | ((3 - c) << shift);
            if i + 1 >= k {
                f(self.mix(fwd.min(rev)));
            }
        }
    }
}

/// Full `u128` k-mer space supporting `k = 5..=64`.
///
/// For `k <= 32`, `mix` and `unmix` match [`KmerSpace64`] bit-for-bit.
/// For `k = 33..=64`, a 3-round bijective `2k`-bit mixer provides balanced
/// avalanche across all `2k` bits while preserving exact inversion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KmerSpace128 {
    pub k: usize,
    bits: u32,
    mask: u128,
    low_space: Option<KmerSpace64>,
    shifts: [u32; 4],
    mul: [u128; 3],
    mul_inv: [u128; 3],
}

impl KmerSpace128 {
    pub fn new(k: usize) -> Self {
        assert!(
            (MIN_K..=MAX_K).contains(&k),
            "k must be within {MIN_K}..={MAX_K}, got {k}"
        );
        let bits = 2 * k as u32;
        let mask = if bits == 128 {
            u128::MAX
        } else {
            (1u128 << bits) - 1
        };
        let low_space = (k <= MAX_K_U64).then(|| KmerSpace64::new(k));
        let scale = |s: u32| (s * bits / 128).max(1);
        let mul = [
            (MIX128_MUL_0 & mask) | 1,
            (MIX128_MUL_1 & mask) | 1,
            (MIX128_MUL_2 & mask) | 1,
        ];
        let mul_inv = [
            inverse_mod_2_128(mul[0]) & mask,
            inverse_mod_2_128(mul[1]) & mask,
            inverse_mod_2_128(mul[2]) & mask,
        ];
        Self {
            k,
            bits,
            mask,
            low_space,
            shifts: [scale(59), scale(53), scale(61), scale(47)],
            mul,
            mul_inv,
        }
    }

    #[inline]
    pub fn bits(&self) -> u32 {
        self.bits
    }

    #[inline]
    pub fn mask(&self) -> u128 {
        self.mask
    }

    #[inline]
    pub fn mix(&self, mut x: u128) -> u128 {
        if let Some(low) = self.low_space {
            return low.mix(x as u64) as u128;
        }
        let [s1, s2, s3, s4] = self.shifts;
        x &= self.mask;
        x ^= x >> s1;
        x = x.wrapping_mul(self.mul[0]) & self.mask;
        x ^= x >> s2;
        x = x.wrapping_mul(self.mul[1]) & self.mask;
        x ^= x >> s3;
        x = x.wrapping_mul(self.mul[2]) & self.mask;
        x ^ (x >> s4)
    }

    #[inline]
    pub fn unmix(&self, mut y: u128) -> u128 {
        if let Some(low) = self.low_space {
            return low.unmix(y as u64) as u128;
        }
        let [s1, s2, s3, s4] = self.shifts;
        y &= self.mask;
        y = undo_xorshift_u128(y, s4, self.bits);
        y = y.wrapping_mul(self.mul_inv[2]) & self.mask;
        y = undo_xorshift_u128(y, s3, self.bits);
        y = y.wrapping_mul(self.mul_inv[1]) & self.mask;
        y = undo_xorshift_u128(y, s2, self.bits);
        y = y.wrapping_mul(self.mul_inv[0]) & self.mask;
        undo_xorshift_u128(y, s1, self.bits)
    }

    #[inline]
    pub fn reverse_complement(&self, kmer: u128) -> u128 {
        reverse_complement_u128(kmer, self.k)
    }

    pub fn decode(&self, kmer: u128) -> String {
        (0..self.k)
            .rev()
            .map(|i| b"ACGT"[((kmer >> (2 * i)) & 3) as usize] as char)
            .collect()
    }

    pub fn canonical(&self, text: &str) -> Option<u128> {
        if text.len() != self.k {
            return None;
        }
        let mut out = None;
        self.for_each_canonical(text.as_bytes(), |x| out = Some(x));
        out
    }

    /// Calls `f` with every canonical packed k-mer in `seq`. Any non-ACGT byte
    /// breaks the rolling window cleanly.
    pub fn for_each_canonical(&self, seq: &[u8], mut f: impl FnMut(u128)) {
        let k = self.k;
        let shift = 2 * (k - 1);
        let (mut fwd, mut rev, mut filled) = (0u128, 0u128, 0usize);
        for &byte in seq {
            let c = CODE[byte as usize];
            if c > 3 {
                filled = 0;
                continue;
            }
            let c = c as u128;
            fwd = ((fwd << 2) | c) & self.mask;
            rev = (rev >> 2) | ((3 - c) << shift);
            filled += 1;
            if filled >= k {
                f(fwd.min(rev));
            }
        }
    }

    /// Unpacks canonical mixed k-mers directly from a 2-bit packed super-k-mer
    /// (4 bases per byte, LSB-first).
    #[inline]
    pub fn for_each_mixed_in_packed_2bit(
        &self,
        packed: &[u8],
        len: usize,
        mut f: impl FnMut(u128),
    ) {
        let k = self.k;
        if len < k {
            return;
        }
        assert!(
            packed.len() >= len.div_ceil(4),
            "packed 2-bit buffer too short for {len} bases"
        );
        let shift = 2 * (k - 1);
        let (mut fwd, mut rev) = (0u128, 0u128);
        for i in 0..len {
            let c = ((packed[i >> 2] >> ((i & 3) * 2)) & 3) as u128;
            fwd = ((fwd << 2) | c) & self.mask;
            rev = (rev >> 2) | ((3 - c) << shift);
            if i + 1 >= k {
                f(self.mix(fwd.min(rev)));
            }
        }
    }
}

/// Generic trait unifying the `u64` (`k = 5..=32`) and `u128` (`k = 33..=64`)
/// canonical k-mer paths without runtime dispatch in inner loops.
pub trait KmerKey:
    Copy + Ord + Eq + std::hash::Hash + Default + Send + Sync + std::fmt::Debug + 'static
{
    const WIDTH_BYTES: usize;
    const WIDTH_BITS: u32;
    type Space: Copy + Send + Sync;

    fn new_space(k: usize) -> Self::Space;
    fn mix(space: &Self::Space, x: Self) -> Self;
    fn unmix(space: &Self::Space, y: Self) -> Self;
    fn reverse_complement(space: &Self::Space, x: Self) -> Self;
    fn decode(space: &Self::Space, x: Self) -> String;
    fn canonical(space: &Self::Space, text: &str) -> Option<Self>;
    fn for_each_mixed_in_packed_2bit(
        space: &Self::Space,
        packed: &[u8],
        len: usize,
        f: impl FnMut(Self),
    );
    fn bucket_of(self, k: usize, bucket_bits: u32) -> usize;
    fn bucket_floor(bucket: usize, k: usize, bucket_bits: u32) -> Self;
    fn to_u128(self) -> u128;
    fn from_u128(v: u128) -> Self;
    fn with_le_bytes<R>(self, f: impl FnOnce(&[u8]) -> R) -> R;
    fn to_hex(self) -> String;
    fn radix_sort(values: &mut Vec<Self>, scratch: &mut Vec<Self>);
}

impl KmerKey for u64 {
    const WIDTH_BYTES: usize = 8;
    const WIDTH_BITS: u32 = 64;
    type Space = KmerSpace64;

    #[inline]
    fn new_space(k: usize) -> Self::Space {
        KmerSpace64::new(k)
    }

    #[inline]
    fn mix(space: &Self::Space, x: Self) -> Self {
        space.mix(x)
    }

    #[inline]
    fn unmix(space: &Self::Space, y: Self) -> Self {
        space.unmix(y)
    }

    #[inline]
    fn reverse_complement(space: &Self::Space, x: Self) -> Self {
        space.reverse_complement(x)
    }

    fn decode(space: &Self::Space, x: Self) -> String {
        space.decode(x)
    }

    fn canonical(space: &Self::Space, text: &str) -> Option<Self> {
        space.canonical(text)
    }

    #[inline]
    fn for_each_mixed_in_packed_2bit(
        space: &Self::Space,
        packed: &[u8],
        len: usize,
        f: impl FnMut(Self),
    ) {
        space.for_each_mixed_in_packed_2bit(packed, len, f);
    }

    #[inline]
    fn bucket_of(self, k: usize, bucket_bits: u32) -> usize {
        debug_assert!(bucket_bits >= 1 && bucket_bits <= 2 * k as u32);
        (self >> (2 * k as u32 - bucket_bits)) as usize
    }

    #[inline]
    fn bucket_floor(bucket: usize, k: usize, bucket_bits: u32) -> Self {
        debug_assert!(bucket_bits >= 1 && bucket_bits <= 2 * k as u32);
        (bucket as u64) << (2 * k as u32 - bucket_bits)
    }

    #[inline]
    fn to_u128(self) -> u128 {
        self as u128
    }

    #[inline]
    fn from_u128(v: u128) -> Self {
        v as u64
    }

    #[inline]
    fn with_le_bytes<R>(self, f: impl FnOnce(&[u8]) -> R) -> R {
        f(&self.to_le_bytes())
    }

    fn to_hex(self) -> String {
        format!("{self:016x}")
    }

    fn radix_sort(values: &mut Vec<Self>, scratch: &mut Vec<Self>) {
        crate::kset::radix_sort(values, scratch);
    }
}

impl KmerKey for u128 {
    const WIDTH_BYTES: usize = 16;
    const WIDTH_BITS: u32 = 128;
    type Space = KmerSpace128;

    #[inline]
    fn new_space(k: usize) -> Self::Space {
        KmerSpace128::new(k)
    }

    #[inline]
    fn mix(space: &Self::Space, x: Self) -> Self {
        space.mix(x)
    }

    #[inline]
    fn unmix(space: &Self::Space, y: Self) -> Self {
        space.unmix(y)
    }

    #[inline]
    fn reverse_complement(space: &Self::Space, x: Self) -> Self {
        space.reverse_complement(x)
    }

    fn decode(space: &Self::Space, x: Self) -> String {
        space.decode(x)
    }

    fn canonical(space: &Self::Space, text: &str) -> Option<Self> {
        space.canonical(text)
    }

    #[inline]
    fn for_each_mixed_in_packed_2bit(
        space: &Self::Space,
        packed: &[u8],
        len: usize,
        f: impl FnMut(Self),
    ) {
        space.for_each_mixed_in_packed_2bit(packed, len, f);
    }

    #[inline]
    fn bucket_of(self, k: usize, bucket_bits: u32) -> usize {
        debug_assert!(bucket_bits >= 1 && bucket_bits <= 2 * k as u32);
        (self >> (2 * k as u32 - bucket_bits)) as usize
    }

    #[inline]
    fn bucket_floor(bucket: usize, k: usize, bucket_bits: u32) -> Self {
        debug_assert!(bucket_bits >= 1 && bucket_bits <= 2 * k as u32);
        (bucket as u128) << (2 * k as u32 - bucket_bits)
    }

    #[inline]
    fn to_u128(self) -> u128 {
        self
    }

    #[inline]
    fn from_u128(v: u128) -> Self {
        v
    }

    #[inline]
    fn with_le_bytes<R>(self, f: impl FnOnce(&[u8]) -> R) -> R {
        f(&self.to_le_bytes())
    }

    fn to_hex(self) -> String {
        format!("{self:032x}")
    }

    fn radix_sort(values: &mut Vec<Self>, scratch: &mut Vec<Self>) {
        radix_sort_u128(values, scratch);
    }
}

/// LSD 8-bit radix sort for `u128`, skipping constant byte positions.
pub fn radix_sort_u128(values: &mut Vec<u128>, scratch: &mut Vec<u128>) {
    let n = values.len();
    if n < 256 {
        values.sort_unstable();
        return;
    }
    let mut hist = vec![[0usize; 256]; 16];
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

#[inline]
fn undo_xorshift_u128(y: u128, s: u32, bits: u32) -> u128 {
    let mut x = y;
    let mut covered = s;
    while covered < bits {
        x = y ^ (x >> s);
        covered += s;
    }
    x
}

#[inline]
fn inverse_mod_2_128(a: u128) -> u128 {
    debug_assert!(a & 1 == 1);
    let mut inv = a;
    for _ in 0..7 {
        inv = inv.wrapping_mul(2u128.wrapping_sub(a.wrapping_mul(inv)));
    }
    inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;
    use ring::digest::{SHA256, digest};

    fn pack_2bit(seq: &[u8]) -> Vec<u8> {
        let mut out = vec![0u8; seq.len().div_ceil(4)];
        for (i, &b) in seq.iter().enumerate() {
            let c = CODE[b as usize];
            assert!(c <= 3, "non-ACGT byte in pack_2bit");
            out[i >> 2] |= c << ((i & 3) * 2);
        }
        out
    }

    #[test]
    fn mix128_constants_match_sha256_seeds() {
        let want_halves = [MIX128_HALF_0, MIX128_HALF_1, MIX128_HALF_2, MIX128_HALF_3];
        for (idx, &constant) in want_halves.iter().enumerate() {
            let seed = format!("{MIX128_SEED_PREFIX}{idx}");
            let h = digest(&SHA256, seed.as_bytes());
            let derived = u64::from_le_bytes(h.as_ref()[..8].try_into().unwrap()) | 1;
            assert_eq!(constant, derived, "seed {seed}");
            assert_eq!(constant & 1, 1, "seed {seed} must be odd");
        }
        assert_eq!(
            MIX128_MUL_0,
            ((SPLITMIX64_MUL_1 as u128) << 64) | (SPLITMIX64_MUL_2 as u128)
        );
        assert_eq!(
            MIX128_MUL_1,
            ((MIX128_HALF_0 as u128) << 64) | (MIX128_HALF_1 as u128)
        );
        assert_eq!(
            MIX128_MUL_2,
            ((MIX128_HALF_2 as u128) << 64) | (MIX128_HALF_3 as u128)
        );
    }

    #[test]
    fn swar_reverse_complement_matches_dna_for_all_u64_k() {
        let mut rng = ChaCha8Rng::seed_from_u64(101);
        for k in MIN_K..=MAX_K_U64 {
            let dna_space = dna::KmerSpace::new(k);
            let k64 = KmerSpace64::new(k);
            let k128 = KmerSpace128::new(k);
            let mask = k64.mask();
            assert_eq!(k64.bits(), dna_space.bits());
            assert_eq!(k128.mask(), mask as u128);
            for _ in 0..2_000 {
                let x = rng.random::<u64>() & mask;
                let rc = k64.reverse_complement(x);
                assert_eq!(rc, dna_space.reverse_complement(x));
                assert_eq!(k128.reverse_complement(x as u128), rc as u128);
                let m = k64.mix(x);
                assert_eq!(m, dna_space.mix(x));
                assert_eq!(k128.mix(x as u128), m as u128);
                assert_eq!(k64.unmix(m), x);
                assert_eq!(k128.unmix(m as u128), x as u128);
            }
        }
    }

    #[test]
    fn u128_golden_mix_unmix_vectors_for_k33_45_61_64() {
        // Golden (k, x, mix(x)) vectors pin all three 128-bit multipliers and
        // all four shift rounds against any accidental constant change.
        let vectors: &[(usize, u128, u128)] = &[
            // k = 33
            (33, 0, 0),
            (33, 1, 0x0000_0000_0000_0002_c6bd_83ad_7ea0_0ebd),
            (
                33,
                0x0000_0000_0000_0003_ffff_ffff_ffff_ffff,
                0x0000_0000_0000_0003_dca9_93c1_2f6c_9514,
            ),
            (33, 1u128 << 64, 0x0000_0000_0000_0002_a20c_7343_d591_b2ad),
            (
                33,
                (1u128 << 64) - 1,
                0x0000_0000_0000_0000_9d81_f873_9446_2c38,
            ),
            (
                33,
                0x0000_0000_0000_0003_fedc_ba98_7654_3210,
                0x0000_0000_0000_0001_49cb_45b7_f531_cbc6,
            ),
            // k = 45
            (45, 0, 0),
            (45, 1, 0x0000_0000_0250_c429_0854_2abe_8ff7_c921),
            (
                45,
                0x0000_0000_03ff_ffff_ffff_ffff_ffff_ffff,
                0x0000_0000_00b3_b096_c07c_333b_b423_f3d8,
            ),
            (45, 1u128 << 64, 0x0000_0000_0189_d23e_fba5_123c_66bd_7f76),
            (
                45,
                (1u128 << 64) - 1,
                0x0000_0000_0362_f343_c680_a589_d382_cd57,
            ),
            (
                45,
                0x0000_0000_01ab_cdef_fedc_ba98_7654_3210,
                0x0000_0000_03e5_4841_e543_016e_bc37_202b,
            ),
            // k = 61
            (61, 0, 0),
            (61, 1, 0x023e_c830_ade2_519d_0b46_8de6_0294_a4b9),
            (
                61,
                0x03ff_ffff_ffff_ffff_ffff_ffff_ffff_ffff,
                0x01d8_108a_1031_c0f9_e432_b02c_204f_b00f,
            ),
            (61, 1u128 << 64, 0x007b_080b_e46f_be7e_dc2c_95d3_b5d2_795f),
            (
                61,
                (1u128 << 64) - 1,
                0x02f0_ba31_e0a5_6ea9_1e13_356d_29bb_2061,
            ),
            (
                61,
                0x0123_4567_89ab_cdef_fedc_ba98_7654_3210,
                0x0347_a06a_fdaf_8764_194e_0aac_0d4e_f7b3,
            ),
            // k = 64
            (64, 0, 0),
            (64, 1, 0xdbc7_bf54_1232_993a_4b7d_2a48_04f5_8550),
            (64, u128::MAX, 0x5ff0_7f61_29d4_0d54_fd70_ae30_9ab7_162a),
            (64, 1u128 << 64, 0x9cbf_a7d3_433e_c957_fe92_457d_8521_eebd),
            (
                64,
                (1u128 << 64) - 1,
                0xd8af_c574_455f_ea5d_af43_82c2_0b05_80b7,
            ),
            (
                64,
                0x0123_4567_89ab_cdef_fedc_ba98_7654_3210,
                0xbe0e_2d6d_301d_2d17_99aa_e2cc_f6d5_cc1e,
            ),
        ];
        for &(k, x, want_y) in vectors {
            let space = KmerSpace128::new(k);
            assert_eq!(space.mix(x), want_y, "mix mismatch at k={k} x={x:#x}");
            assert_eq!(
                space.unmix(want_y),
                x,
                "unmix mismatch at k={k} y={want_y:#x}"
            );
        }

        // 61-mers differing only in high bits above bit 64 or only in low bits.
        let s61 = KmerSpace128::new(61);
        let prefix32 = "ACGT".repeat(8);
        let kmer_low_a = format!("{prefix32}{}", "A".repeat(29));
        let kmer_low_b = format!("{prefix32}{}C", "A".repeat(28));
        let suffix32 = "TGCA".repeat(8);
        let kmer_hi_a = format!("{}{suffix32}", "C".repeat(29));
        let kmer_hi_b = format!("{}G{suffix32}", "C".repeat(28));

        let c_la = s61.canonical(&kmer_low_a).unwrap();
        let c_lb = s61.canonical(&kmer_low_b).unwrap();
        let c_ha = s61.canonical(&kmer_hi_a).unwrap();
        let c_hb = s61.canonical(&kmer_hi_b).unwrap();
        let pack_raw = |s: &str| {
            s.bytes()
                .fold(0u128, |acc, b| (acc << 2) | (CODE[b as usize] as u128))
        };
        // Pair hi shares all 64 low raw bits and differs only above bit 64.
        assert_eq!(pack_raw(&kmer_hi_a) as u64, pack_raw(&kmer_hi_b) as u64);
        assert_ne!(pack_raw(&kmer_hi_a), pack_raw(&kmer_hi_b));

        let keys = [c_la, c_lb, c_ha, c_hb];
        for i in 0..keys.len() {
            assert_eq!(s61.unmix(s61.mix(keys[i])), keys[i]);
            for j in (i + 1)..keys.len() {
                assert_ne!(keys[i], keys[j]);
                assert_ne!(s61.mix(keys[i]), s61.mix(keys[j]));
            }
        }
    }

    #[test]
    fn u128_round_trips_511_values_at_every_k_from_33_to_64() {
        // 11 fixed boundary values plus 500 pseudorandom draws = 511 values
        // per k across all 32 values of k in 33..=64 (16,352 values total).
        let mut rng = ChaCha8Rng::seed_from_u64(3364);
        for k in (MAX_K_U64 + 1)..=MAX_K {
            let space = KmerSpace128::new(k);
            let mask = space.mask();
            assert_eq!(space.bits() as usize, 2 * k);

            let mut samples = Vec::with_capacity(511);
            samples.extend_from_slice(&[
                0u128,
                1,
                2,
                3,
                mask,
                mask - 1,
                mask >> 1,
                (mask >> 1) + 1,
                1u128 << 64,
                (1u128 << 64) - 1,
                (1u128 << 64) + 1,
            ]);
            for _ in 0..500 {
                let raw = ((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128);
                samples.push(raw & mask);
            }
            assert_eq!(samples.len(), 511);

            for x in samples {
                let y = space.mix(x);
                assert!(y <= mask, "k={k}: mix escaped mask");
                assert_eq!(space.unmix(y), x, "k={k}: unmix(mix(x)) != x");
                assert_eq!(space.mix(space.unmix(x)), x, "k={k}: mix(unmix(x)) != x");
                let rc = space.reverse_complement(x);
                assert!(rc <= mask, "k={k}: rc escaped mask");
                assert_eq!(space.reverse_complement(rc), x, "k={k}: double rc failed");
                let text = space.decode(x);
                assert_eq!(text.len(), k);
                let canon = space.canonical(&text).unwrap();
                assert_eq!(canon, x.min(rc), "k={k}: canonical mismatch");
                let rc_text = space.decode(rc);
                assert_eq!(space.canonical(&rc_text).unwrap(), canon);
            }
        }
    }

    #[test]
    fn u128_mixer_balances_top_bits_at_k61() {
        let space = KmerSpace128::new(61);
        let mut rng = ChaCha8Rng::seed_from_u64(303);
        let seq: Vec<u8> = (0..100_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let half = 1u128 << (space.bits() - 1);
        let mut low = 0usize;
        let mut total = 0usize;
        space.for_each_canonical(&seq, |x| {
            total += 1;
            low += (space.mix(x) < half) as usize;
        });
        let share = low as f64 / total as f64;
        assert!((share - 0.5).abs() < 0.01, "share was {share}");
    }

    #[test]
    fn packed_2bit_scanner_matches_canonical_mix_for_u64_and_u128() {
        let mut rng = ChaCha8Rng::seed_from_u64(707);
        for k in [5usize, 15, 31, 32] {
            let space = KmerSpace64::new(k);
            for len in [0usize, k - 1, k, k + 1, k + 19] {
                let seq: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
                let packed = pack_2bit(&seq);
                let mut want = Vec::new();
                space.for_each_canonical(&seq, |x| want.push(space.mix(x)));
                let mut got = Vec::new();
                <u64 as KmerKey>::for_each_mixed_in_packed_2bit(&space, &packed, len, |y| {
                    got.push(y);
                });
                assert_eq!(got, want, "u64 packed mismatch k={k} len={len}");
            }
        }
        for k in [33usize, 45, 61, 64] {
            let space = KmerSpace128::new(k);
            for len in [0usize, k - 1, k, k + 1, k + 19] {
                let seq: Vec<u8> = (0..len).map(|_| b"ACGT"[rng.random_range(0..4)]).collect();
                let packed = pack_2bit(&seq);
                let mut want = Vec::new();
                space.for_each_canonical(&seq, |x| want.push(space.mix(x)));
                let mut got = Vec::new();
                <u128 as KmerKey>::for_each_mixed_in_packed_2bit(&space, &packed, len, |y| {
                    got.push(y);
                });
                assert_eq!(got, want, "u128 packed mismatch k={k} len={len}");
            }
        }
    }

    #[test]
    fn kmer_key_bucketing_and_hex_cover_min_and_max_k() {
        for k in MIN_K..=MAX_K_U64 {
            let s = u64::new_space(k);
            for b in [0usize, 1, 511, 1023] {
                let floor = u64::bucket_floor(b, k, 10);
                assert_eq!(floor.bucket_of(k, 10), b, "u64 k={k} b={b}");
            }
            let top = s.mask();
            assert_eq!(top.bucket_of(k, 10), 1023, "u64 top k={k}");
            assert_eq!(top.to_hex().len(), 16);
        }
        for k in MIN_K..=MAX_K {
            let s = u128::new_space(k);
            for b in [0usize, 1, 511, 1023] {
                let floor = u128::bucket_floor(b, k, 10);
                assert_eq!(floor.bucket_of(k, 10), b, "u128 k={k} b={b}");
            }
            let top = s.mask();
            assert_eq!(top.bucket_of(k, 10), 1023, "u128 top k={k}");
            assert_eq!(top.to_hex().len(), 32);
        }
    }

    #[test]
    fn u128_radix_sort_matches_std_sort() {
        let mut rng = ChaCha8Rng::seed_from_u64(404);
        for n in [0, 17, 255, 256, 5_000] {
            let mut v: Vec<u128> = (0..n)
                .map(|_| ((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128))
                .collect();
            let mut want = v.clone();
            want.sort_unstable();
            radix_sort_u128(&mut v, &mut Vec::new());
            assert_eq!(v, want, "n={n}");
        }
    }
}
