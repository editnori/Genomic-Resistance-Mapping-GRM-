//! Bitsets over genomes. A feature is one bitset (bit i = present in genome
//! i), and sample masks, labels, and tree nodes are bitsets over the same
//! columns, so every count a learner needs is an AND plus a popcount.

/// Number of u64 words that hold `n` bits.
pub fn words_for(n: usize) -> usize {
    n.div_ceil(64)
}

pub fn count(a: &[u64]) -> u32 {
    a.iter().map(|w| w.count_ones()).sum()
}

pub fn and_count(a: &[u64], b: &[u64]) -> u32 {
    a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum()
}

pub fn and3_count(a: &[u64], b: &[u64], c: &[u64]) -> u32 {
    a.iter()
        .zip(b)
        .zip(c)
        .map(|((x, y), z)| (x & y & z).count_ones())
        .sum()
}

pub fn get(a: &[u64], i: usize) -> bool {
    (a[i >> 6] >> (i & 63)) & 1 == 1
}

pub fn set(a: &mut [u64], i: usize) {
    a[i >> 6] |= 1 << (i & 63);
}

/// Owned bitset with a fixed length.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Bits {
    len: usize,
    words: Vec<u64>,
}

impl Bits {
    pub fn zeros(len: usize) -> Self {
        Bits {
            len,
            words: vec![0; words_for(len)],
        }
    }

    pub fn ones(len: usize) -> Self {
        let mut b = Bits {
            len,
            words: vec![u64::MAX; words_for(len)],
        };
        b.clear_tail();
        b
    }

    pub fn from_indices(len: usize, idx: impl IntoIterator<Item = usize>) -> Self {
        let mut b = Bits::zeros(len);
        for i in idx {
            b.insert(i);
        }
        b
    }

    pub fn from_words(len: usize, words: Vec<u64>) -> Self {
        assert_eq!(words.len(), words_for(len));
        let mut b = Bits { len, words };
        b.clear_tail();
        b
    }

    fn clear_tail(&mut self) {
        let rem = self.len % 64;
        if rem != 0 {
            if let Some(last) = self.words.last_mut() {
                *last &= (1u64 << rem) - 1;
            }
        }
    }

    pub fn len(&self) -> usize {
        self.len
    }

    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    pub fn words(&self) -> &[u64] {
        &self.words
    }

    pub fn insert(&mut self, i: usize) {
        debug_assert!(i < self.len);
        set(&mut self.words, i);
    }

    pub fn remove(&mut self, i: usize) {
        self.words[i >> 6] &= !(1u64 << (i & 63));
    }

    pub fn contains(&self, i: usize) -> bool {
        get(&self.words, i)
    }

    pub fn count(&self) -> u32 {
        count(&self.words)
    }

    pub fn and(&self, other: &Bits) -> Bits {
        self.zip(other, |a, b| a & b)
    }

    pub fn and_not(&self, other: &Bits) -> Bits {
        self.zip(other, |a, b| a & !b)
    }

    pub fn or(&self, other: &Bits) -> Bits {
        self.zip(other, |a, b| a | b)
    }

    pub fn not(&self) -> Bits {
        let mut b = Bits {
            len: self.len,
            words: self.words.iter().map(|w| !w).collect(),
        };
        b.clear_tail();
        b
    }

    /// `self & row` where `row` is a raw feature bitset of the same length.
    pub fn and_slice(&self, row: &[u64]) -> Bits {
        Bits {
            len: self.len,
            words: self.words.iter().zip(row).map(|(a, b)| a & b).collect(),
        }
    }

    pub fn and_not_slice(&self, row: &[u64]) -> Bits {
        Bits {
            len: self.len,
            words: self.words.iter().zip(row).map(|(a, b)| a & !b).collect(),
        }
    }

    fn zip(&self, other: &Bits, f: impl Fn(u64, u64) -> u64) -> Bits {
        assert_eq!(self.len, other.len);
        Bits {
            len: self.len,
            words: self
                .words
                .iter()
                .zip(&other.words)
                .map(|(&a, &b)| f(a, b))
                .collect(),
        }
    }

    pub fn iter_ones(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(wi, &w)| {
            let mut w = w;
            std::iter::from_fn(move || {
                if w == 0 {
                    return None;
                }
                let t = w.trailing_zeros() as usize;
                w &= w - 1;
                Some(wi * 64 + t)
            })
        })
    }
}

/// Row-major bit matrix: row r is feature r's bitset over `cols` genomes.
#[derive(Clone, Debug)]
pub struct BitMatrix {
    pub rows: usize,
    pub cols: usize,
    pub words: usize,
    pub data: Vec<u64>,
}

impl BitMatrix {
    pub fn with_cols(cols: usize) -> Self {
        BitMatrix {
            rows: 0,
            cols,
            words: words_for(cols),
            data: Vec::new(),
        }
    }

    pub fn row(&self, r: usize) -> &[u64] {
        &self.data[r * self.words..(r + 1) * self.words]
    }

    pub fn push_row(&mut self, row: &[u64]) {
        debug_assert_eq!(row.len(), self.words);
        self.data.extend_from_slice(row);
        self.rows += 1;
    }

    pub fn get(&self, r: usize, c: usize) -> bool {
        get(self.row(r), c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn set_algebra_and_counts() {
        let a = Bits::from_indices(130, [0, 5, 64, 129]);
        let b = Bits::from_indices(130, [5, 64, 100]);
        assert_eq!(a.and(&b).iter_ones().collect::<Vec<_>>(), vec![5, 64]);
        assert_eq!(a.and_not(&b).iter_ones().collect::<Vec<_>>(), vec![0, 129]);
        assert_eq!(a.or(&b).count(), 5);
        // `not` must not set bits past the length.
        assert_eq!(a.not().count(), 126);
        assert_eq!(Bits::ones(130).count(), 130);
        assert_eq!(and_count(a.words(), b.words()), 2);
        assert_eq!(
            and3_count(a.words(), b.words(), Bits::from_indices(130, [64]).words()),
            1
        );
    }
}
