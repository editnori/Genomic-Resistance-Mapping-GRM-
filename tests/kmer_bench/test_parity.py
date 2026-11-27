"""Unit and planted-fault tests for scripts/kmer-bench/parity.py."""

import pathlib
import struct
import tempfile
import unittest

from tests.kmer_bench import common

parity = common.load_module("parity")
synthetic = common.load_module("synthetic")


class ParityReferenceTests(unittest.TestCase):
    def test_fnv1a64_matches_published_vectors_from_src_fnv(self):
        # Exact reference vectors from src/fnv.rs:41-49.
        vectors = (
            (b"", 0xCBF2_9CE4_8422_2325),
            (b"a", 0xAF63_DC4C_8601_EC8C),
            (b"foobar", 0x8594_4171_F739_67E8),
        )
        for data, want in vectors:
            h = parity.Fnv1a64()
            h.update_bytes(data)
            self.assertEqual(h.state, want, data)
        self.assertEqual(parity.Fnv1a64(0xAB).hex(), "00000000000000ab")
        # Length-prefixed str hashing distinguishes ("ab","c") from ("a","bc").
        h1, h2 = parity.Fnv1a64(), parity.Fnv1a64()
        h1.update_str("ab")
        h1.update_str("c")
        h2.update_str("a")
        h2.update_str("bc")
        self.assertNotEqual(h1.hex(), h2.hex())

    def test_canonical_encoding_and_decode_match_src_dna(self):
        # Exact vector from src/dna.rs:199-206.
        space = parity.KmerSpace(5)
        fwd = space.canonical("AACGT")
        rc = space.canonical("ACGTT")
        self.assertIsNotNone(fwd)
        self.assertEqual(fwd, rc)
        self.assertEqual(space.decode(fwd), "AACGT")
        self.assertEqual(space.decode(space.reverse_complement(fwd)), "ACGTT")
        self.assertEqual(space.decode(space.reverse_complement_swar(fwd)), "ACGTT")
        self.assertIsNone(space.canonical("AACG"))
        self.assertIsNone(space.canonical("AACGN"))

    def test_n_reset_and_lowercase_handling_match_src_dna(self):
        # Exact vectors from src/dna.rs:209-218.
        space = parity.KmerSpace(5)
        seen = [space.decode(x) for x in space.iter_canonical(b"ACGTANACGTA")]
        self.assertEqual(len(seen), 2)
        self.assertEqual(seen, ["ACGTA", "ACGTA"])
        lower = list(space.iter_canonical(b"acgtanacgta"))
        upper = list(space.iter_canonical(b"ACGTANACGTA"))
        self.assertEqual(lower, upper)
        all_two = list(space.iter_canonical(b"ACGTAC"))
        self.assertEqual(len(all_two), 2)

    def test_mix_unmix_roundtrip_and_bounds_across_all_k(self):
        rng = synthetic.DeterministicRng(7)
        for k in range(parity.MIN_K, parity.MAX_K + 1):
            space = parity.KmerSpace(k)
            mask = parity.MASK64 if k == 32 else (1 << (2 * k)) - 1
            edges = [0, 1, mask, mask - 1, mask >> 1]
            samples = edges + [rng.next_u64() & mask for _ in range(120)]
            for x in samples:
                y = space.mix(x)
                self.assertGreaterEqual(y, 0)
                self.assertLessEqual(y, mask, f"k={k} escaped mask")
                self.assertEqual(space.unmix(y), x, f"k={k} round-trip failed")

    def test_canonical_matches_naive_reference_both_strands_all_k(self):
        # Mirrors src/dna.rs:221-264 across all k in 5..=32 with soft-masked bases and N runs.
        alphabet = "ACGTACGTACGTACGTacgtN"
        rng = synthetic.DeterministicRng(13)
        seq = "".join(alphabet[rng.randrange(len(alphabet))] for _ in range(320))
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        rc_seq = seq.translate(comp)[::-1]

        def pack(window: str):
            val = 0
            for ch in window.upper():
                idx = "ACGT".find(ch)
                if idx < 0:
                    return None
                val = (val << 2) | idx
            return val

        for k in range(parity.MIN_K, parity.MAX_K + 1):
            space = parity.KmerSpace(k)
            want = []
            for i in range(len(seq) - k + 1):
                w = seq[i : i + k]
                p = pack(w)
                if p is None:
                    continue
                rc_p = pack(w.translate(comp)[::-1])
                want.append(min(p, rc_p))
            got = list(space.iter_canonical(seq))
            self.assertEqual(got, want, f"k={k} forward strand")
            got_rc = list(space.iter_canonical(rc_seq))
            got_rc.reverse()
            self.assertEqual(got_rc, want, f"k={k} reverse strand")

    def test_swar_reverse_complement_matches_scalar_all_k(self):
        rng = synthetic.DeterministicRng(99)
        for k in range(parity.MIN_K, parity.MAX_K + 1):
            space = parity.KmerSpace(k)
            mask = space.mask
            cases = [0, 1, 2, 3, mask, 0xAAAA_AAAA_AAAA_AAAA & mask, 0x5555_5555_5555_5555 & mask]
            cases.extend(rng.next_u64() & mask for _ in range(60))
            for x in cases:
                self.assertEqual(
                    space.reverse_complement_swar(x),
                    space.reverse_complement(x),
                    f"k={k} x={x:#x}",
                )

    def test_kset_v1_binary_roundtrip_edges_and_corruption(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for k in (5, 21, 31, 32):
                space = parity.KmerSpace(k)
                low = space.bits - parity.BUCKET_BITS
                top = space.mask
                values = [0, 1, top]
                for b in (1, 511, 512, parity.BUCKETS - 1):
                    values.extend([(b << low) - 1, b << low])
                values = sorted(set(v & top for v in values))
                path = root / f"k{k}.kset"
                nbytes = parity.write_kset(path, f"syn_{k}", space, values)
                self.assertEqual(nbytes, path.stat().st_size)
                loaded = parity.read_kset(path)
                self.assertEqual(loaded.k, k)
                self.assertEqual(loaded.name, f"syn_{k}")
                self.assertEqual(list(loaded.values), values)
                self.assertEqual(loaded.digest, parity.kset_digest(values))

            # Planted repeat gap (0 after first element) must be rejected (src/kset.rs:367).
            space31 = parity.KmerSpace(31)
            valid_path = root / "valid.kset"
            parity.write_kset(valid_path, "g", space31, [7, 12])
            raw = bytearray(valid_path.read_bytes())
            raw[-1] = 0
            bad_path = root / "repeat.kset"
            bad_path.write_bytes(raw)
            with self.assertRaisesRegex(ValueError, "repeats a k-mer"):
                parity.read_kset(bad_path)

            # Corrupt magic must be rejected.
            raw_magic = bytearray(valid_path.read_bytes())
            raw_magic[0:8] = b"BADMAGIC"
            bad_magic = root / "badmagic.kset"
            bad_magic.write_bytes(raw_magic)
            with self.assertRaisesRegex(ValueError, "invalid kset magic"):
                parity.read_kset(bad_magic)

    def test_external_dump_parsers_tsv_and_fasta(self):
        space = parity.KmerSpace(7)
        k1 = space.canonical("ACGTACG")
        k2 = space.canonical("TTTTGGG")
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            tsv = root / "dump.tsv"
            tsv.write_text(f"ACGTACG\t5\nCCCAAAA\t3\n", encoding="utf-8")
            fa = root / "dump.fa"
            fa.write_text(f">5\nCGTACGT\n>3\nTTTTGGG\n", encoding="utf-8")
            c_tsv = parity.parse_kmer_text_dump(tsv, k=7, min_abundance=1)
            c_fa = parity.parse_kmer_text_dump(fa, k=7, min_abundance=1)
            self.assertEqual(c_tsv, {k1: 5, k2: 3})
            self.assertEqual(c_fa, {k1: 5, k2: 3})
            self.assertEqual(parity.parse_kmer_text_dump(tsv, k=7, min_abundance=4), {k1: 5})

    def test_planted_fault_set_parity_fails_on_drop_add_or_mutate(self):
        space = parity.KmerSpace(15)
        seq = synthetic.DeterministicRng(55).random_dna(500)
        res = parity.count_sequences([("c0", seq.encode("ascii"))], k=15, min_abundance=1)
        expected = sorted(res.counts.keys())
        self.assertGreater(len(expected), 400)

        # Clean passes.
        clean = parity.verify_set_parity(expected, expected, k=15)
        self.assertTrue(clean.ok)
        self.assertEqual(clean.expected_fnv1a64, res.fnv1a64_digest)

        # Fault 1: drop 1 k-mer -> MUST fail.
        dropped = parity.verify_set_parity(expected, expected[1:], k=15)
        self.assertFalse(dropped.ok)
        self.assertEqual(dropped.missing_count, 1)
        self.assertEqual(dropped.extra_count, 0)
        self.assertNotEqual(dropped.expected_fnv1a64, dropped.observed_fnv1a64)

        # Fault 2: add 1 extra k-mer -> MUST fail.
        extra_val = 0
        while extra_val in res.counts:
            extra_val += 1
        added = parity.verify_set_parity(expected, expected + [extra_val], k=15)
        self.assertFalse(added.ok)
        self.assertEqual(added.missing_count, 0)
        self.assertEqual(added.extra_count, 1)

        # Fault 3: mutate 1 k-mer (same distinct count!) -> MUST fail.
        mutated = parity.verify_set_parity(expected, expected[1:] + [extra_val], k=15)
        self.assertFalse(mutated.ok)
        self.assertEqual(mutated.expected_distinct, mutated.observed_distinct)
        self.assertEqual(mutated.missing_count, 1)
        self.assertEqual(mutated.extra_count, 1)
        self.assertNotEqual(mutated.expected_fnv1a64, mutated.observed_fnv1a64)

    def test_planted_fault_abundance_parity_fails_on_altered_count_or_cap(self):
        seq = synthetic.DeterministicRng(77).random_dna(250)
        # Repeat a 60-bp block 3 times so counts span 1 and >=3.
        composite = seq + seq[:60] + seq[:60]
        res = parity.count_sequences([("c0", composite.encode("ascii"))], k=15, min_abundance=1)
        obs = dict(res.counts)
        self.assertTrue(parity.verify_abundance_parity(res.counts, obs, k=15, min_abundance=1).ok)

        # Fault: alter 1 abundance count by +1 -> MUST fail.
        first_km = min(obs)
        obs_bad = dict(obs)
        obs_bad[first_km] += 1
        rep = parity.verify_abundance_parity(res.counts, obs_bad, k=15, min_abundance=1)
        self.assertFalse(rep.ok)
        self.assertEqual(rep.count_mismatch, 1)
        self.assertNotEqual(rep.expected_spectrum_sha256, rep.observed_spectrum_sha256)

        # Cap semantics: count=300 vs 255 fails uncapped, passes with max_cap=255 (FastK u8).
        exp_high = {first_km: 300}
        obs_capped = {first_km: 255}
        self.assertFalse(parity.verify_abundance_parity(exp_high, obs_capped, k=15, max_cap=None).ok)
        self.assertTrue(parity.verify_abundance_parity(exp_high, obs_capped, k=15, max_cap=255).ok)

    def test_planted_fault_unitig_parity_fails_on_merge_or_split(self):
        # Mirrors src/input_build.rs:1299-1333: 1 SNP creates 4 unitigs.
        back = synthetic.DeterministicRng(21).random_dna(240)
        var = list(back)
        var[120] = "C" if var[120] == "A" else "A"
        var_s = "".join(var)
        space = parity.KmerSpace(21)
        union_kmers = set(space.iter_canonical(back)) | set(space.iter_canonical(var_s))
        unitigs = parity.compact_unitigs(union_kmers, k=21)
        self.assertEqual(len(unitigs), 4)
        self.assertTrue(parity.verify_unitig_parity(unitigs, unitigs).ok)

        # Fault 1: merge two unitigs -> MUST fail.
        merged = [unitigs[0] | unitigs[1], *unitigs[2:]]
        rep_m = parity.verify_unitig_parity(unitigs, merged)
        self.assertFalse(rep_m.ok)
        self.assertEqual(rep_m.unitigs_both, 2)
        self.assertEqual(rep_m.unitigs_expected_only, 2)
        self.assertEqual(rep_m.unitigs_observed_only, 1)

        # Fault 2: split one unitig -> MUST fail.
        big = sorted(unitigs[0])
        split = [frozenset(big[: len(big) // 2]), frozenset(big[len(big) // 2 :]), *unitigs[1:]]
        rep_s = parity.verify_unitig_parity(unitigs, split)
        self.assertFalse(rep_s.ok)
        self.assertEqual(rep_s.unitigs_expected_only, 1)

    def test_planted_fault_matrix_parity_fails_on_flipped_bit(self):
        rng = synthetic.DeterministicRng(31)
        core = rng.random_dna(120)
        acc = rng.random_dna(120)
        space = parity.KmerSpace(15)
        sets = [
            set(space.iter_canonical(core)),
            set(space.iter_canonical(core + "NNNN" + acc)),
            set(space.iter_canonical(acc + "NNNN" + core)),
            set(space.iter_canonical(core)),
        ]
        ref = parity.build_reference_matrix(sets, k=15, min_support=2)
        self.assertEqual(ref.n_patterns, 1)
        self.assertGreater(ref.kmers_kept, 50)
        self.assertTrue(parity.verify_matrix_parity(ref, ref.kmer_to_mask).ok)

        # Fault: flip 1 genome bit on 1 k-mer -> MUST fail.
        bad = dict(ref.kmer_to_mask)
        km = min(bad)
        bad[km] ^= 1
        rep = parity.verify_matrix_parity(ref, bad)
        self.assertFalse(rep.ok)
        self.assertEqual(rep.mask_mismatch, 1)


if __name__ == "__main__":
    unittest.main()
