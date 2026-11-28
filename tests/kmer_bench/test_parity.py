"""Unit, golden-grm, and planted-fault tests for scripts/kmer-bench/parity.py."""

import hashlib
import json
import pathlib
import struct
import tempfile
import unittest

from tests.kmer_bench import common

parity = common.load_module("parity")
synthetic = common.load_module("synthetic")


def _load_golden():
    return json.loads(common.GOLDEN_GRM.read_text(encoding="utf-8"))


def _build_kset_fixture_fasta(k: int) -> bytes:
    rng = synthetic.DeterministicRng(2000 + k)
    s1a = rng.random_dna(180)
    s1b = rng.random_dna(90).lower()
    s1c = rng.random_dna(70)
    s2 = rng.random_dna(160)
    return f">c1\n{s1a}NNNNN{s1b}\n@{s1c}\n>c2\n{s2}\n".encode("ascii")


def _build_6genome_fastas():
    rng = synthetic.DeterministicRng(3001)
    b6 = rng.random_dna(140)
    b5 = rng.random_dna(130)
    b4 = rng.random_dna(120)
    b3 = rng.random_dna(110)
    b2 = rng.random_dna(100)
    uniques = [rng.random_dna(90) for _ in range(6)]
    genomes = {}
    for g in range(6):
        contigs = [f">core\n{b6}"]
        if g < 5:
            contigs.append(f">s5\n{b5}")
        if g < 4:
            contigs.append(f">s4\n{b4}")
        if g < 3:
            contigs.append(f">s3\n{b3}")
        if g < 2:
            contigs.append(f">s2\n{b2}")
        contigs.append(f">u{g}\n{uniques[g]}")
        genomes[f"syn_g{g:04d}"] = ("\n".join(contigs) + "\n").encode("ascii")
    return genomes


def _build_10genome_fastas():
    rng = synthetic.DeterministicRng(3002)
    core = rng.random_dna(100)
    blocks = [rng.random_dna(80) for _ in range(5)]
    uniques = [rng.random_dna(60) for _ in range(10)]
    genomes = {}
    for g in range(10):
        contigs = [f">core\n{core}"]
        if g in (0, 1, 8, 9):
            contigs.append(f">b0\n{blocks[0]}")
        if g in (1, 2, 7, 8, 9):
            contigs.append(f">b1\n{blocks[1]}")
        if g in (0, 9):
            contigs.append(f">b2\n{blocks[2]}")
        if g in (2, 3, 4, 5, 6, 7, 8):
            contigs.append(f">b3\n{blocks[3]}")
        if g != 9:
            contigs.append(f">b4\n{blocks[4]}")
        contigs.append(f">u{g}\n{uniques[g]}")
        genomes[f"syn_g{g:04d}"] = ("\n".join(contigs) + "\n").encode("ascii")
    return genomes


def _write_synthetic_grm_matrix(matrix_dir: pathlib.Path, ref_mat, genome_names):
    """Writes a GRM binary matrix directory from a ReferenceMatrix for parser tests."""
    matrix_dir.mkdir(parents=True, exist_ok=True)
    space = parity.KmerSpace(ref_mat.k)
    n_genomes = len(genome_names)
    words = (n_genomes + 63) // 64
    pattern_to_id = {mask: idx for idx, mask in enumerate(ref_mat.pattern_masks)}
    by_pat = [[] for _ in ref_mat.pattern_masks]
    for km, mask in ref_mat.kmer_to_mask.items():
        by_pat[pattern_to_id[mask]].append(space.mix(km))
    for lst in by_pat:
        lst.sort()

    pat_bytes = bytearray()
    sup_bytes = bytearray()
    off_bytes = bytearray()
    kmer_bytes = bytearray()
    offset = 0
    for pid, mask in enumerate(ref_mat.pattern_masks):
        for wi in range(words):
            pat_bytes.extend(struct.pack("<Q", (mask >> (64 * wi)) & parity.MASK64))
        sup_bytes.extend(struct.pack("<I", mask.bit_count()))
        off_bytes.extend(struct.pack("<Q", offset))
        for mv in by_pat[pid]:
            kmer_bytes.extend(struct.pack("<Q", mv))
        offset += len(by_pat[pid])
    off_bytes.extend(struct.pack("<Q", offset))

    (matrix_dir / "patterns.bin").write_bytes(pat_bytes)
    (matrix_dir / "support.bin").write_bytes(sup_bytes)
    (matrix_dir / "kmer_offsets.bin").write_bytes(off_bytes)
    (matrix_dir / "kmers.bin").write_bytes(kmer_bytes)
    summary = {
        "k": ref_mat.k,
        "genomes": list(genome_names),
        "min_support": ref_mat.min_support,
        "n_patterns": ref_mat.n_patterns,
        "kmers_distinct": ref_mat.kmers_distinct,
        "kmers_kept": ref_mat.kmers_kept,
        "kmers_rare": ref_mat.kmers_rare,
        "kmers_core": ref_mat.kmers_core,
        "kset_digests": [],
    }
    (matrix_dir / "matrix.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


class ParityReferenceTests(unittest.TestCase):
    def test_fnv1a64_matches_published_vectors_and_little_endian_u64(self):
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
        # Explicitly verify little-endian u64 byte order (catches F09).
        val = 0x0102_0304_0506_0708
        h_u64 = parity.Fnv1a64()
        h_u64.update_u64(val)
        h_le = parity.Fnv1a64()
        h_le.update_bytes(b"\x08\x07\x06\x05\x04\x03\x02\x01")
        h_be = parity.Fnv1a64()
        h_be.update_bytes(b"\x01\x02\x03\x04\x05\x06\x07\x08")
        self.assertEqual(h_u64.hex(), h_le.hex())
        self.assertNotEqual(h_u64.hex(), h_be.hex())
        self.assertEqual(parity.kset_digest([val]), h_le.hex())
        # Length-prefixed str hashing distinguishes ("ab","c") from ("a","bc").
        h1, h2 = parity.Fnv1a64(), parity.Fnv1a64()
        h1.update_str("ab")
        h1.update_str("c")
        h2.update_str("a")
        h2.update_str("bc")
        self.assertNotEqual(h1.hex(), h2.hex())

    def test_golden_grm_mix_unmix_digest_and_kset_bytes(self):
        # Verifies parity.py against committed real grm-written values (catches F01, F02, F09, F11).
        golden = _load_golden()
        self.assertEqual(golden["format"], "grm-kmer-bench-golden-v1")
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for k_str, rows in golden["mix_vectors"].items():
                k = int(k_str)
                space = parity.KmerSpace(k)
                for row in rows:
                    x = int(row["x"])
                    want_mixed = int(row["mixed"])
                    self.assertEqual(space.mix(x), want_mixed, f"k={k} x={x}")
                    self.assertEqual(space.unmix(want_mixed), x, f"k={k} y={want_mixed}")
                    self.assertEqual(parity.kset_digest([want_mixed]), row["fnv1a64_digest"])
                    out = root / f"{row['name']}.kset"
                    nbytes = parity.write_kset(out, row["name"], space, [want_mixed])
                    raw = out.read_bytes()
                    self.assertEqual(nbytes, row["kset_bytes"])
                    self.assertEqual(raw[:24].hex(), row["header_24_hex"])
                    self.assertEqual(raw[10:12], b"\x00\x00")
                    self.assertEqual(hashlib.sha256(raw).hexdigest(), row["kset_sha256"])

    def test_golden_grm_multicontig_kset_including_at_continuation_line(self):
        # Verifies multi-contig FASTA with N run, lowercase, and '@' continuation line against grm (catches F09, F11, F15).
        golden = _load_golden()
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for k_str, fx in golden["kset_fixtures"].items():
                k = int(k_str)
                fa = root / f"{fx['name']}.fna"
                fa.write_bytes(_build_kset_fixture_fasta(k))
                res = parity.count_file(fa, k=k, min_abundance=1)
                self.assertEqual(res.distinct_kept, fx["n_kmers"], f"k={k}")
                self.assertEqual(res.fnv1a64_digest, fx["fnv1a64_digest"], f"k={k}")
                kset_path = root / f"{fx['name']}.kset"
                nbytes = parity.write_kset(kset_path, fx["name"], parity.KmerSpace(k), res.mixed_sorted)
                raw = kset_path.read_bytes()
                self.assertEqual(nbytes, fx["kset_bytes"], f"k={k}")
                self.assertEqual(raw[:24].hex(), fx["header_24_hex"], f"k={k}")
                self.assertEqual(hashlib.sha256(raw).hexdigest(), fx["kset_sha256"], f"k={k}")

    def test_fasta_at_continuation_line_not_read_as_fastq(self):
        # Directly catches F15: '@' inside a FASTA sequence must append to the record, not trigger FASTQ parsing.
        with tempfile.TemporaryDirectory() as tmp:
            fa = pathlib.Path(tmp) / "at_line.fna"
            fa.write_text(">c1\nACGTACGT\n@TTTTGGGG\nCCCCCCCC\n>c2\nAAAATTTT\n", encoding="ascii")
            recs = list(parity.iter_sequence_records(fa))
            self.assertEqual(
                recs,
                [
                    ("c1", b"ACGTACGT@TTTTGGGGCCCCCCCC"),
                    ("c2", b"AAAATTTT"),
                ],
            )
            res = parity.count_file(fa, k=5, min_abundance=1)
            space = parity.KmerSpace(5)
            # Both 'TTTTG' (after '@') and 'CCCCC' (on the line after '@' that FASTQ would skip) must be counted.
            self.assertIn(space.canonical("TTTTG"), res.counts)
            self.assertIn(space.canonical("CCCCC"), res.counts)

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

            # Writer must reject duplicate or descending k-mers (catches F27).
            space31 = parity.KmerSpace(31)
            with self.assertRaisesRegex(ValueError, "strictly ascending"):
                parity.write_kset(root / "dup.kset", "g", space31, [7, 7])
            with self.assertRaisesRegex(ValueError, "strictly ascending"):
                parity.write_kset(root / "desc.kset", "g", space31, [12, 7])

            # Planted repeat gap (0 after first element) must be rejected (F13).
            valid_path = root / "valid.kset"
            parity.write_kset(valid_path, "g", space31, [7, 12])
            raw = bytearray(valid_path.read_bytes())
            raw[-1] = 0
            bad_path = root / "repeat.kset"
            bad_path.write_bytes(raw)
            with self.assertRaisesRegex(ValueError, "repeats a k-mer"):
                parity.read_kset(bad_path)

            # Corrupt n_kmers header total must be rejected (catches F12).
            raw_tot = bytearray(valid_path.read_bytes())
            struct.pack_into("<Q", raw_tot, 12, 999)
            bad_tot = root / "badtotal.kset"
            bad_tot.write_bytes(raw_tot)
            with self.assertRaisesRegex(ValueError, "inconsistent bucket table"):
                parity.read_kset(bad_tot)

            # Corrupt varint gap spilling past bucket range must be rejected (catches F14).
            space6 = parity.KmerSpace(6)
            valid6 = root / "valid6.kset"
            # At k=6 (12 bits, 10 bucket bits), each bucket holds 4 values (0..3).
            parity.write_kset(valid6, "g", space6, [1])
            raw_rng = bytearray(valid6.read_bytes())
            raw_rng[-1] = 5  # gap=5 puts value 5 into bucket 1, violating bucket 0 range
            bad_rng = root / "badrange.kset"
            bad_rng.write_bytes(raw_rng)
            with self.assertRaisesRegex(ValueError, "does not match count and range"):
                parity.read_kset(bad_rng)

            # Corrupt magic must be rejected.
            raw_magic = bytearray(valid_path.read_bytes())
            raw_magic[0:8] = b"BADMAGIC"
            bad_magic = root / "badmagic.kset"
            bad_magic.write_bytes(raw_magic)
            with self.assertRaisesRegex(ValueError, "invalid kset magic"):
                parity.read_kset(bad_magic)

    def test_external_dump_parsers_tsv_fasta_and_sum_both_strands(self):
        # Verifies TSV and FASTA dump parsing AND summing counts across reverse-complement strands (catches F18).
        space = parity.KmerSpace(7)
        k1 = space.canonical("ACGTACG")
        k2 = space.canonical("TTTTGGG")
        # ACGTACG and its revcomp CGTACGT map to the same canonical k-mer k1.
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            tsv = root / "dump.tsv"
            tsv.write_text("ACGTACG\t5\nCCCAAAA\t3\n", encoding="utf-8")
            fa = root / "dump.fa"
            fa.write_text(">5\nCGTACGT\n>3\nTTTTGGG\n", encoding="utf-8")
            c_tsv = parity.parse_kmer_text_dump(tsv, k=7, min_abundance=1)
            c_fa = parity.parse_kmer_text_dump(fa, k=7, min_abundance=1)
            self.assertEqual(c_tsv, {k1: 5, k2: 3})
            self.assertEqual(c_fa, {k1: 5, k2: 3})
            self.assertEqual(parity.parse_kmer_text_dump(tsv, k=7, min_abundance=4), {k1: 5})

            # Both strands present in same dump: counts MUST sum (5 + 7 = 12, not 7).
            both_tsv = root / "both.tsv"
            both_tsv.write_text("ACGTACG\t5\nCGTACGT\t7\nTTTTGGG\t2\nCCCAAAA\t4\n", encoding="utf-8")
            both_fa = root / "both.fa"
            both_fa.write_text(">5\nACGTACG\n>7\nCGTACGT\n>2\nTTTTGGG\n>4\nCCCAAAA\n", encoding="utf-8")
            self.assertEqual(parity.parse_kmer_text_dump(both_tsv, k=7, min_abundance=1), {k1: 12, k2: 6})
            self.assertEqual(parity.parse_kmer_text_dump(both_fa, k=7, min_abundance=1), {k1: 12, k2: 6})

    def test_planted_fault_set_parity_fails_on_drop_add_or_mutate(self):
        seq = synthetic.DeterministicRng(55).random_dna(500)
        res = parity.count_sequences([("c0", seq.encode("ascii"))], k=15, min_abundance=1)
        expected = sorted(res.counts.keys())
        self.assertGreater(len(expected), 400)

        clean = parity.verify_set_parity(expected, expected, k=15)
        self.assertTrue(clean.ok)
        self.assertEqual(clean.expected_fnv1a64, res.fnv1a64_digest)

        dropped = parity.verify_set_parity(expected, expected[1:], k=15)
        self.assertFalse(dropped.ok)
        self.assertEqual(dropped.missing_count, 1)
        self.assertEqual(dropped.extra_count, 0)
        self.assertNotEqual(dropped.expected_fnv1a64, dropped.observed_fnv1a64)

        extra_val = 0
        while extra_val in res.counts:
            extra_val += 1
        added = parity.verify_set_parity(expected, expected + [extra_val], k=15)
        self.assertFalse(added.ok)
        self.assertEqual(added.missing_count, 0)
        self.assertEqual(added.extra_count, 1)

        mutated = parity.verify_set_parity(expected, expected[1:] + [extra_val], k=15)
        self.assertFalse(mutated.ok)
        self.assertEqual(mutated.expected_distinct, mutated.observed_distinct)
        self.assertEqual(mutated.missing_count, 1)
        self.assertEqual(mutated.extra_count, 1)
        self.assertNotEqual(mutated.expected_fnv1a64, mutated.observed_fnv1a64)

    def test_planted_fault_abundance_parity_fails_on_altered_count_or_cap(self):
        seq = synthetic.DeterministicRng(77).random_dna(250)
        composite = seq + seq[:60] + seq[:60]
        res = parity.count_sequences([("c0", composite.encode("ascii"))], k=15, min_abundance=1)
        obs = dict(res.counts)
        self.assertTrue(parity.verify_abundance_parity(res.counts, obs, k=15, min_abundance=1).ok)

        first_km = min(obs)
        obs_bad = dict(obs)
        obs_bad[first_km] += 1
        rep = parity.verify_abundance_parity(res.counts, obs_bad, k=15, min_abundance=1)
        self.assertFalse(rep.ok)
        self.assertEqual(rep.count_mismatch, 1)
        self.assertNotEqual(rep.expected_spectrum_sha256, rep.observed_spectrum_sha256)

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

        merged = [unitigs[0] | unitigs[1], *unitigs[2:]]
        rep_m = parity.verify_unitig_parity(unitigs, merged)
        self.assertFalse(rep_m.ok)
        self.assertEqual(rep_m.unitigs_both, 2)
        self.assertEqual(rep_m.unitigs_expected_only, 2)
        self.assertEqual(rep_m.unitigs_observed_only, 1)

        big = sorted(unitigs[0])
        split = [frozenset(big[: len(big) // 2]), frozenset(big[len(big) // 2 :]), *unitigs[1:]]
        rep_s = parity.verify_unitig_parity(unitigs, split)
        self.assertFalse(rep_s.ok)
        self.assertEqual(rep_s.unitigs_expected_only, 1)

    def test_golden_grm_matrix_drops_core_and_hashes_little_endian(self):
        # Verifies matrix building against committed real grm-written matrix fixtures (catches F23, F24).
        golden = _load_golden()
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            for label, builder in (
                ("cohort6_k31", _build_6genome_fastas),
                ("cohort10_k21", _build_10genome_fastas),
            ):
                fx = golden["matrix_fixtures"][label]
                k = fx["k"]
                fastas = builder()
                gids = list(fastas.keys())
                sets = []
                for gid, raw in fastas.items():
                    fp = root / f"{label}_{gid}.fna"
                    fp.write_bytes(raw)
                    sets.append(set(parity.count_file(fp, k=k, min_abundance=1).counts))

                ref = parity.build_reference_matrix(sets, k=k, min_support=fx["min_support"])
                self.assertEqual(ref.kmers_distinct, fx["kmers_distinct"], label)
                self.assertEqual(ref.kmers_kept, fx["kmers_kept"], label)
                self.assertEqual(ref.kmers_rare, fx["kmers_rare"], label)
                self.assertEqual(ref.kmers_core, fx["kmers_core"], label)
                self.assertEqual(ref.n_patterns, fx["n_patterns"], label)
                self.assertEqual(ref.digest_sha256, fx["digest_sha256"], label)

                # Explicitly assert no rare (< min_support) or core/near-core (> n - min_support) k-mers are kept (F23).
                max_sup = len(gids) - fx["min_support"]
                for mask in ref.kmer_to_mask.values():
                    self.assertGreaterEqual(mask.bit_count(), fx["min_support"])
                    self.assertLessEqual(mask.bit_count(), max_sup)
                self.assertGreater(ref.kmers_core, 0)
                self.assertGreater(ref.kmers_rare, 0)

                # Verify little-endian digest differs from big-endian digest on 10-genome (2-byte mask) matrix (F24).
                space = parity.KmerSpace(k)
                n_bytes = (len(gids) + 7) // 8
                h_be = hashlib.sha256()
                for mv in sorted(space.mix(km) for km in ref.kmer_to_mask):
                    mask = ref.kmer_to_mask[space.unmix(mv)]
                    h_be.update(struct.pack(">Q", mv))
                    h_be.update(mask.to_bytes(n_bytes, "big"))
                self.assertNotEqual(ref.digest_sha256, h_be.hexdigest())

                # Write binary GRM matrix files and verify read_grm_matrix round-trip and exact binary SHA-256.
                mat_dir = root / f"mat_{label}"
                _write_synthetic_grm_matrix(mat_dir, ref, gids)
                self.assertEqual(
                    hashlib.sha256((mat_dir / "patterns.bin").read_bytes()).hexdigest(),
                    fx["patterns_bin_sha256"],
                )
                self.assertEqual(
                    hashlib.sha256((mat_dir / "support.bin").read_bytes()).hexdigest(),
                    fx["support_bin_sha256"],
                )
                self.assertEqual(
                    hashlib.sha256((mat_dir / "kmer_offsets.bin").read_bytes()).hexdigest(),
                    fx["kmer_offsets_bin_sha256"],
                )
                self.assertEqual(
                    hashlib.sha256((mat_dir / "kmers.bin").read_bytes()).hexdigest(),
                    fx["kmers_bin_sha256"],
                )
                loaded = parity.read_grm_matrix(mat_dir)
                self.assertTrue(parity.verify_matrix_parity(ref, loaded.kmer_to_mask).ok)

    def test_golden_grm_unitigs_k31_and_snp_equality(self):
        golden = _load_golden()
        ufx = golden["unitig_fixture"]
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            fastas = _build_6genome_fastas()
            gids = list(fastas.keys())
            sets = []
            for gid, raw in fastas.items():
                fp = root / f"{gid}.fna"
                fp.write_bytes(raw)
                sets.append(set(parity.count_file(fp, k=31, min_abundance=1).counts))

            space = parity.KmerSpace(31)
            union_kmers = set().union(*sets)
            ref_u = parity.compact_unitigs(union_kmers, k=31)
            self.assertEqual(len(union_kmers), ufx["graph_kmers"])
            self.assertEqual(len(ref_u), ufx["unitigs"])

            # Synthesize a GRM exp-input-build unitigs output dir and test verify_grm_unitig_build.
            ordered_u = sorted(ref_u, key=lambda u: min(space.mix(km) for km in u))
            held_sets = [{uid for uid, u in enumerate(ordered_u) if u.issubset(gs)} for gs in sets]
            u_mat = parity.build_reference_matrix(held_sets, k=16, min_support=2)
            self.assertEqual(u_mat.n_patterns, ufx["matrix_patterns"])
            self.assertEqual(u_mat.kmers_kept, ufx["matrix_values_kept"])
            self.assertEqual(u_mat.kmers_rare, ufx["matrix_values_rare"])
            self.assertEqual(u_mat.kmers_core, ufx["matrix_values_core"])

            bdir = root / "grm_unitigs"
            _write_synthetic_grm_matrix(bdir / "matrix", u_mat, gids)
            tsv_lines = ["genome\tunitigs_held\tunitigs_in_part\tkmers_without_unitig"]
            for row in ufx["genomes_tsv_rows"]:
                tsv_lines.append("\t".join(str(x) for x in row))
            (bdir / "genomes.tsv").write_text("\n".join(tsv_lines) + "\n", encoding="utf-8")
            (bdir / "build.json").write_text(
                json.dumps({
                    "kind": "unitigs",
                    "graph_kmers": ufx["graph_kmers"],
                    "unitigs": ufx["unitigs"],
                    "pairs_held": ufx["pairs_held"],
                    "pairs_in_part": ufx["pairs_in_part"],
                    "kmers_without_unitig": ufx["kmers_without_unitig"],
                }),
                encoding="utf-8",
            )
            rep = parity.verify_grm_unitig_build(bdir, sets, gids, k=31)
            self.assertTrue(rep.ok)

            # Corrupt build.json unitig count -> verify_grm_unitig_build MUST fail.
            (bdir / "build.json").write_text(
                json.dumps({
                    "kind": "unitigs",
                    "graph_kmers": ufx["graph_kmers"],
                    "unitigs": ufx["unitigs"] + 1,
                    "pairs_held": ufx["pairs_held"],
                    "pairs_in_part": ufx["pairs_in_part"],
                    "kmers_without_unitig": 0,
                }),
                encoding="utf-8",
            )
            self.assertFalse(parity.verify_grm_unitig_build(bdir, sets, gids, k=31).ok)

        # Verify 1-SNP k=31 golden unitig equality fixture.
        sfx = golden["unitig_snp_equality"]
        rng = synthetic.DeterministicRng(sfx["seed"])
        back = rng.random_dna(sfx["length"])
        var_chars = list(back)
        pos = sfx["snp_pos"]
        var_chars[pos] = "C" if var_chars[pos] == "A" else "A"
        var_s = "".join(var_chars)
        space31 = parity.KmerSpace(31)
        union_snp = set(space31.iter_canonical(back)) | set(space31.iter_canonical(var_s))
        py_u = parity.compact_unitigs(union_snp, k=31)
        self.assertEqual(len(union_snp), sfx["kmers_grm"])
        self.assertEqual(len(py_u), sfx["unitigs_grm"])
        self.assertEqual(sorted(len(u) for u in py_u), sfx["unitig_sizes_sorted"])
        self.assertEqual(sorted(min(u) for u in py_u), sfx["unitig_min_kmers_sorted"])

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

        bad = dict(ref.kmer_to_mask)
        km = min(bad)
        bad[km] ^= 1
        rep = parity.verify_matrix_parity(ref, bad)
        self.assertFalse(rep.ok)
        self.assertEqual(rep.mask_mismatch, 1)


if __name__ == "__main__":
    unittest.main()
