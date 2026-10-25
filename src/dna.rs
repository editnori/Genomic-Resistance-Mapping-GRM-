//! Nucleotide encoding, canonical k-mers, and the invertible k-mer mixer.
//!
//! A k-mer of length k <= 32 packs into the low 2k bits of a u64 (A=0, C=1,
//! G=2, T=3). The canonical form is the smaller of the forward k-mer and its
//! reverse complement, so both strands of an assembly yield one value.
//!
//! Stored k-mers are never the raw packed value. They pass through a bijective
//! mixer on the 2k-bit space first, for two reasons:
//! - canonical values skew toward small numbers (the min of two values), so
//!   bucketing raw values by their top bits would be badly unbalanced;
//! - a uniformly mixed, sorted k-mer set is already a bottom-s MinHash sketch:
//!   its first s values are the s smallest hashes.
//!
//! The mixer is invertible, so any stored value decodes back to its sequence.

/// Shortest k: a kset's 1,024 buckets need 10 bits of each k-mer.
pub const MIN_K: usize = 5;
pub const MAX_K: usize = 32;

const fn code_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
}

/// 2-bit code for each byte; 4 marks anything that is not A, C, G or T.
pub static CODE: [u8; 256] = code_table();

/// The k-mer length plus the mixer constants for its 2k-bit space.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct KmerSpace {
    pub k: usize,
    bits: u32,
    mask: u64,
    shifts: [u32; 3],
    mul: [u64; 2],
    mul_inv: [u64; 2],
}

impl KmerSpace {
    pub fn new(k: usize) -> Self {
        assert!(
            (MIN_K..=MAX_K).contains(&k),
            "k must be within {MIN_K}..={MAX_K}, got {k}"
        );
        let bits = 2 * k as u32;
        let mask = if bits == 64 {
            u64::MAX
        } else {
            (1u64 << bits) - 1
        };
        // splitmix64's finalizer, with shifts scaled to the narrower space.
        let scale = |s: u32| (s * bits / 64).max(1);
        let mul = [
            0xbf58_476d_1ce4_e5b9 & mask | 1,
            0x94d0_49bb_1331_11eb & mask | 1,
        ];
        KmerSpace {
            k,
            bits,
            mask,
            shifts: [scale(30), scale(27), scale(31)],
            mul,
            mul_inv: [
                inverse_mod_2_64(mul[0]) & mask,
                inverse_mod_2_64(mul[1]) & mask,
            ],
        }
    }

    /// Number of bits a k-mer occupies.
    pub fn bits(&self) -> u32 {
        self.bits
    }

    pub fn mix(&self, mut x: u64) -> u64 {
        let [s1, s2, s3] = self.shifts;
        x ^= x >> s1;
        x = x.wrapping_mul(self.mul[0]) & self.mask;
        x ^= x >> s2;
        x = x.wrapping_mul(self.mul[1]) & self.mask;
        x ^ (x >> s3)
    }

    pub fn unmix(&self, mut y: u64) -> u64 {
        let [s1, s2, s3] = self.shifts;
        y = undo_xorshift(y, s3, self.bits);
        y = y.wrapping_mul(self.mul_inv[1]) & self.mask;
        y = undo_xorshift(y, s2, self.bits);
        y = y.wrapping_mul(self.mul_inv[0]) & self.mask;
        undo_xorshift(y, s1, self.bits)
    }

    /// Calls `f` with every canonical packed k-mer in `seq`. Runs of non-ACGT
    /// bytes break the window, so no k-mer spans an N.
    pub fn for_each_canonical(&self, seq: &[u8], mut f: impl FnMut(u64)) {
        let k = self.k;
        let shift = 2 * (k - 1);
        let (mut fwd, mut rev, mut filled) = (0u64, 0u64, 0usize);
        for &byte in seq {
            let c = CODE[byte as usize];
            if c > 3 {
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

    /// Packs an ACGT string (any case) into its canonical k-mer.
    pub fn canonical(&self, text: &str) -> Option<u64> {
        if text.len() != self.k {
            return None;
        }
        let mut out = None;
        self.for_each_canonical(text.as_bytes(), |x| out = Some(x));
        out
    }

    pub fn reverse_complement(&self, kmer: u64) -> u64 {
        let mut out = 0u64;
        let mut x = kmer;
        for _ in 0..self.k {
            out = (out << 2) | (3 - (x & 3));
            x >>= 2;
        }
        out
    }

    pub fn decode(&self, kmer: u64) -> String {
        (0..self.k)
            .rev()
            .map(|i| b"ACGT"[((kmer >> (2 * i)) & 3) as usize] as char)
            .collect()
    }
}

/// Inverse of `y = x ^ (x >> s)` on a `bits`-wide value.
fn undo_xorshift(y: u64, s: u32, bits: u32) -> u64 {
    let mut x = y;
    let mut covered = s;
    while covered < bits {
        x = y ^ (x >> s);
        covered += s;
    }
    x
}

/// Multiplicative inverse of an odd number modulo 2^64 (Newton iteration;
/// each step doubles the number of correct low bits).
fn inverse_mod_2_64(a: u64) -> u64 {
    debug_assert!(a & 1 == 1);
    let mut inv = a;
    for _ in 0..6 {
        inv = inv.wrapping_mul(2u64.wrapping_sub(a.wrapping_mul(inv)));
    }
    inv
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng};
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn mixer_round_trips_and_stays_in_range() {
        let mut rng = ChaCha8Rng::seed_from_u64(7);
        for k in [5, 15, 21, 31, 32] {
            let space = KmerSpace::new(k);
            let mask = if k == 32 {
                u64::MAX
            } else {
                (1u64 << (2 * k)) - 1
            };
            for _ in 0..20_000 {
                let x = rng.random::<u64>() & mask;
                let y = space.mix(x);
                assert!(y <= mask, "k={k} mixed value escaped the space");
                assert_eq!(space.unmix(y), x, "k={k}");
            }
        }
    }

    #[test]
    fn canonical_matches_both_strands() {
        let space = KmerSpace::new(5);
        let fwd = space.canonical("AACGT").unwrap();
        let rc = space.canonical("ACGTT").unwrap();
        assert_eq!(fwd, rc);
        assert_eq!(space.decode(fwd), "AACGT");
        assert_eq!(space.decode(space.reverse_complement(fwd)), "ACGTT");
    }

    #[test]
    fn window_breaks_on_non_acgt() {
        let space = KmerSpace::new(5);
        let mut seen = Vec::new();
        space.for_each_canonical(b"ACGTANACGTA", |x| seen.push(space.decode(x)));
        // Two 5-base runs, each yields exactly one k-mer; none spans the N.
        assert_eq!(seen.len(), 2);
        let mut all = Vec::new();
        space.for_each_canonical(b"ACGTAC", |x| all.push(x));
        assert_eq!(all.len(), 2);
    }

    #[test]
    fn canonical_matches_a_naive_reference_on_both_strands() {
        // Repacks every window from scratch for every supported k, over a
        // soft-masked sequence with N runs. The reverse-complement strand must
        // yield the same values in reverse order.
        let mut rng = ChaCha8Rng::seed_from_u64(13);
        let seq: Vec<u8> = (0..2_000)
            .map(|_| b"ACGTACGTACGTACGTacgtN"[rng.random_range(0..21)])
            .collect();
        let complement = |&b: &u8| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            other => other,
        };
        let rc_seq: Vec<u8> = seq.iter().rev().map(complement).collect();
        let pack = |w: &[u8]| {
            w.iter().try_fold(0u64, |acc, &b| {
                let c = b"ACGT".iter().position(|&x| x == b.to_ascii_uppercase())?;
                Some((acc << 2) | c as u64)
            })
        };
        for k in 5..=MAX_K {
            let space = KmerSpace::new(k);
            let want: Vec<u64> = seq
                .windows(k)
                .filter_map(|w| {
                    let rc: Vec<u8> = w.iter().rev().map(complement).collect();
                    Some(pack(w)?.min(pack(&rc)?))
                })
                .collect();
            let mut got = Vec::new();
            space.for_each_canonical(&seq, |x| got.push(x));
            assert_eq!(got, want, "k={k}");
            let mut got_rc = Vec::new();
            space.for_each_canonical(&rc_seq, |x| got_rc.push(x));
            got_rc.reverse();
            assert_eq!(got_rc, want, "k={k} reverse strand");
        }
    }

    #[test]
    fn mixed_values_balance_across_top_bits() {
        // Canonical k-mers concentrate in the low half of the space; mixed
        // ones must not, or kset buckets would be lopsided.
        let space = KmerSpace::new(31);
        let mut rng = ChaCha8Rng::seed_from_u64(11);
        let seq: Vec<u8> = (0..200_000)
            .map(|_| b"ACGT"[rng.random_range(0..4)])
            .collect();
        let mut raw_low = 0usize;
        let mut mixed_low = 0usize;
        let mut total = 0usize;
        let half = 1u64 << (space.bits() - 1);
        space.for_each_canonical(&seq, |x| {
            total += 1;
            raw_low += (x < half) as usize;
            mixed_low += (space.mix(x) < half) as usize;
        });
        let raw_share = raw_low as f64 / total as f64;
        let mixed_share = mixed_low as f64 / total as f64;
        assert!(
            raw_share > 0.7,
            "raw canonical share in low half was {raw_share}"
        );
        assert!(
            (mixed_share - 0.5).abs() < 0.01,
            "mixed share was {mixed_share}"
        );
    }
}
