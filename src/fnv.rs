//! 64-bit FNV-1a over fixed-width little-endian values, for digests that
//! later builds compare: kset contents in counted.json and report input
//! fingerprints. `DefaultHasher` or `FxHash` output may change between Rust
//! or crate releases.

pub(crate) struct Fnv(pub(crate) u64);

impl Fnv {
    pub(crate) fn new() -> Self {
        Fnv(0xcbf2_9ce4_8422_2325)
    }

    pub(crate) fn bytes(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 = (self.0 ^ u64::from(b)).wrapping_mul(0x0000_0100_0000_01b3);
        }
    }

    pub(crate) fn u64(&mut self, v: u64) {
        self.bytes(&v.to_le_bytes());
    }

    /// Length first, so adjacent strings cannot trade characters.
    pub(crate) fn str(&mut self, s: &str) {
        self.u64(s.len() as u64);
        self.bytes(s.as_bytes());
    }

    pub(crate) fn hex(&self) -> String {
        format!("{:016x}", self.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv1a_matches_published_values() {
        for (text, digest) in [
            ("", 0xcbf2_9ce4_8422_2325),
            ("a", 0xaf63_dc4c_8601_ec8c),
            ("foobar", 0x8594_4171_f739_67e8),
        ] {
            let mut h = Fnv::new();
            h.bytes(text.as_bytes());
            assert_eq!(h.0, digest, "{text:?}");
        }
        assert_eq!(Fnv(0xab).hex(), "00000000000000ab");
    }
}
