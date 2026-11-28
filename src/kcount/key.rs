//! Dual-width canonical k-mer representation (`k = 5..=64`), branchless O(1)
//! SWAR reverse complement, and invertible 2k-bit mixers for `u64` and `u128`.
//!
//! For `k <= 32`, `KmerSpace64` uses the exact same 2-bit encoding and
//! bijective mixer constants as [`crate::dna::KmerSpace`], guaranteeing
//! bit-for-bit parity with `.kset` v1 digests and sketches. For `k = 33..=64`,
//! `KmerSpace128` packs `2k` bits into a `u128` with a 3-round invertible
//! 128-bit mixer so top-bit bucketing and bottom-`s` MinHash sketching remain
//! uniform and lossless across the full `2k`-bit space.

use crate::dna::{self, CODE};

/// Minimum supported k-mer length (needs at least 10 bits for 1,024 buckets).
pub const MIN_K: usize = dna::MIN_K;
/// Maximum k-mer length that fits in the fast `u64` path.
pub const MAX_K_U64: usize = dna::MAX_K;
/// Maximum supported k-mer length across both `u64` and `u128` paths.
pub const MAX_K: usize = 64;

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
        let bits = inner.bits();
        let mask = if bits == 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };
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

    /// Calls `f` with every canonical packed k-mer in `seq`, treating any base
    /// whose Phred+33 quality is strictly below `min_qual` like an `N` break.
    pub fn for_each_canonical_with_qual(
        &self,
        seq: &[u8],
        qual: Option<&[u8]>,
        min_qual: u8,
        mut f: impl FnMut(u64),
    ) {
        if min_qual == 0 || qual.is_none() {
            self.for_each_canonical(seq, f);
            return;
        }
        let qual = qual.unwrap();
        debug_assert_eq!(seq.len(), qual.len());
        let k = self.k;
        let shift = 2 * (k - 1);
        let qual_floor = 33u8.saturating_add(min_qual);
        let (mut fwd, mut rev, mut filled) = (0u64, 0u64, 0usize);
        for (&byte, &q) in seq.iter().zip(qual.iter()) {
            let c = CODE[byte as usize];
            if c > 3 || q < qual_floor {
                filled = 0;
                continue;
            }
            let c = c as u64;
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
    pub fn for_each_mixed_in_packed_2bit(&self, packed: &[u8], len: usize, mut f: impl FnMut(u64)) {
        let k = self.k;
        if len < k {
            return;
        }
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
/// For `k = 33..=64`, a 3-round bijective `2k`-bit mixer provides uniform
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
            (0xbf58_476d_1ce4_e5b9_94d0_49bb_1331_11eb_u128 & mask) | 1,
            (0x9e37_79b9_7f4a_7c15_f39c_c060_5ced_c835_u128 & mask) | 1,
            (0xd6e8_feb8_6659_fd93_a076_1d64_78bd_642f_u128 & mask) | 1,
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

    /// Calls `f` with every canonical packed k-mer in `seq`, treating any base
    /// whose Phred+33 quality is strictly below `min_qual` like an `N` break.
    pub fn for_each_canonical_with_qual(
        &self,
        seq: &[u8],
        qual: Option<&[u8]>,
        min_qual: u8,
        mut f: impl FnMut(u128),
    ) {
        if min_qual == 0 || qual.is_none() {
            self.for_each_canonical(seq, f);
            return;
        }
        let qual = qual.unwrap();
        debug_assert_eq!(seq.len(), qual.len());
        let k = self.k;
        let shift = 2 * (k - 1);
        let qual_floor = 33u8.saturating_add(min_qual);
        let (mut fwd, mut rev, mut filled) = (0u128, 0u128, 0usize);
        for (&byte, &q) in seq.iter().zip(qual.iter()) {
            let c = CODE[byte as usize];
            if c > 3 || q < qual_floor {
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
        (self >> (2 * k as u32 - bucket_bits)) as usize
    }

    #[inline]
    fn bucket_floor(bucket: usize, k: usize, bucket_bits: u32) -> Self {
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
        (self >> (2 * k as u32 - bucket_bits)) as usize
    }

    #[inline]
    fn bucket_floor(bucket: usize, k: usize, bucket_bits: u32) -> Self {
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

    #[test]
    fn swar_reverse_complement_matches_dna_for_all_u64_k() {
        let mut rng = ChaCha8Rng::seed_from_u64(101);
        for k in MIN_K..=MAX_K_U64 {
            let dna_space = dna::KmerSpace::new(k);
            let k64 = KmerSpace64::new(k);
            let mask = k64.mask();
            for _ in 0..2_000 {
                let x = rng.random::<u64>() & mask;
                assert_eq!(k64.reverse_complement(x), dna_space.reverse_complement(x));
                assert_eq!(k64.mix(x), dna_space.mix(x));
                assert_eq!(k64.unmix(k64.mix(x)), x);
            }
        }
    }

    #[test]
    fn u128_mixer_and_swar_rc_round_trip_across_wide_k() {
        let mut rng = ChaCha8Rng::seed_from_u64(202);
        for k in [5, 21, 31, 32, 33, 45, 61, 63, 64] {
            let space = KmerSpace128::new(k);
            let mask = space.mask();
            for _ in 0..5_000 {
                let raw = ((rng.random::<u64>() as u128) << 64) | (rng.random::<u64>() as u128);
                let x = raw & mask;
                let y = space.mix(x);
                assert!(y <= mask, "k={k}: mixed escaped mask");
                assert_eq!(space.unmix(y), x, "k={k}: unmix failed");
                let rc = space.reverse_complement(x);
                assert!(rc <= mask, "k={k}: rc escaped mask");
                assert_eq!(space.reverse_complement(rc), x, "k={k}: double rc failed");
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
