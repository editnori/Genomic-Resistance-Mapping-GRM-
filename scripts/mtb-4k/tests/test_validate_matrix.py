"""validate_matrix.py against a matrix the real GRM binary builds.

usage: GRM=target/release/grm python3 -I scripts/mtb-4k/tests/test_validate_matrix.py

The pristine matrix must pass, including an independent FASTA recount that
reproduces GRM's k-mer count, FNV kset digest and sketch; each of the
review's eight corruptions must fail, as must a recount of a FASTA edited
after counting. The re-review's omissions (N3) must fail too: with every
genome recounted, against a superset matrix, and in the region every sketch
covers.
"""
import contextlib
import gzip
import hashlib
import json
import os
import random
import shutil
import struct
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
VALIDATOR = HERE.parent / "validate_matrix.py"
SLUG = "mycobacterium_tuberculosis"
FILES = ["matrix.json", "patterns.bin", "support.bin", "kmer_offsets.bin", "kmers.bin"]


def grm():
    path = os.environ.get("GRM")
    if not path or not Path(path).is_file():
        raise SystemExit("set GRM to a built grm binary (cargo build --release)")
    return path


def build_root(base: Path, n=70, seed=5):
    rng = random.Random(seed)
    root = base / "data"
    sp = root / "species" / SLUG
    sp.mkdir(parents=True)
    (root / "fasta").mkdir()
    core = "".join(rng.choice("ACGT") for _ in range(1500))
    blocks = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(6)]
    ids = [f"G{i:03d}" for i in range(n)]
    for i, g in enumerate(ids):
        seq = list(core)
        for pos in rng.sample(range(len(core)), 12):
            seq[pos] = rng.choice("ACGT")
        extra = [b for j, b in enumerate(blocks) if (i >> j) & 1]
        text = ">c1\n" + "".join(seq) + "\n>c2\n" + "NN".join(extra or ["ACGTN"]) + "\n"
        with gzip.open(root / "fasta" / f"{g}.fna.gz", "wt") as f:
            f.write(text)
    (sp / "genomes.txt").write_text("\n".join(ids) + "\n")
    count_and_matrix(root)
    return root


def count_and_matrix(root: Path, min_support=2):
    subprocess.run([grm(), "--data", str(root), "count", "--species", SLUG, "--allow-missing-catalog"], check=True, capture_output=True)
    subprocess.run([grm(), "--data", str(root), "matrix", "--species", SLUG, "--min-support", str(min_support)],
                   check=True, capture_output=True)


def build_subset(base: Path, source: Path, ids, min_support=2):
    """A root over `ids` (in the source's order), sharing its FASTA."""
    root = base / "data"
    sp = root / "species" / SLUG
    sp.mkdir(parents=True)
    (root / "fasta").mkdir()
    for g in ids:
        os.symlink(source / "fasta" / f"{g}.fna.gz", root / "fasta" / f"{g}.fna.gz")
    (sp / "genomes.txt").write_text("\n".join(ids) + "\n")
    count_and_matrix(root, min_support)
    return root


def load(m):
    s = json.loads((m / "matrix.json").read_text())
    n = len(s["genomes"])
    words = (n + 63) // 64
    pat = list(struct.unpack(f"<{s['n_patterns'] * words}Q", (m / "patterns.bin").read_bytes()))
    sup = list(struct.unpack(f"<{s['n_patterns']}I", (m / "support.bin").read_bytes()))
    off = list(struct.unpack(f"<{s['n_patterns'] + 1}Q", (m / "kmer_offsets.bin").read_bytes()))
    km = list(struct.unpack(f"<{off[-1]}Q", (m / "kmers.bin").read_bytes()))
    return s, words, pat, sup, off, km


def save(m, s, pat, sup, off, km):
    (m / "matrix.json").write_text(json.dumps(s, indent=2))
    (m / "patterns.bin").write_bytes(struct.pack(f"<{len(pat)}Q", *pat))
    (m / "support.bin").write_bytes(struct.pack(f"<{len(sup)}I", *sup))
    (m / "kmer_offsets.bin").write_bytes(struct.pack(f"<{len(off)}Q", *off))
    (m / "kmers.bin").write_bytes(struct.pack(f"<{len(km)}Q", *km))


# The review's eight corruptions (tools/adv_validator.py), applied in place.
def dup_kmer_across_patterns(s, w, pat, sup, off, km):
    x = km[off[0]]
    km[off[1]:off[2]] = sorted(km[off[1]:off[2]] + [x])
    for i in range(2, len(off)):
        off[i] += 1
    s["kmers_kept"] += 1


def dup_pattern_row(s, w, pat, sup, off, km):
    pat[w:2 * w] = pat[0:w]
    sup[1] = sup[0]


def kmer_outside_k_space(s, w, pat, sup, off, km):
    km[off[1] - 1] = (1 << 63) | 1


def swap_two_genome_columns(s, w, pat, sup, off, km):
    for r in range(len(sup)):
        if ((pat[r * w] >> 0) & 1) != ((pat[r * w] >> 1) & 1):
            pat[r * w] ^= 0b11


def empty_pattern(s, w, pat, sup, off, km):
    km[off[0]:off[2]] = sorted(km[off[0]:off[2]])
    off[1] = off[0]


def bogus_distinct_and_rare(s, w, pat, sup, off, km):
    s["kmers_distinct"] *= 10
    s["kmers_rare"] = 0
    s["kmers_core"] = 12345


def digests_swapped(s, w, pat, sup, off, km):
    s["kset_digests"][0], s["kset_digests"][1] = s["kset_digests"][1], s["kset_digests"][0]


def support_below_filter(s, w, pat, sup, off, km):
    pat[0:w] = [0] * w
    pat[0] = 1
    sup[0] = 1


MUTANTS = [dup_kmer_across_patterns, dup_pattern_row, kmer_outside_k_space, swap_two_genome_columns,
           empty_pattern, bogus_distinct_and_rare, digests_swapped, support_below_filter]


# The re-review's omissions (tools/rr_validator.py), each in the pattern
# whose smallest k-mer is largest, far above every sketch threshold.
def high_pattern(s, w, pat, sup, off, km, min_len=1):
    return max((p for p in range(len(sup)) if off[p + 1] - off[p] >= min_len), key=lambda p: km[off[p]])


def drop_pattern(s, w, pat, sup, off, km):
    p = high_pattern(s, w, pat, sup, off, km)
    n = off[p + 1] - off[p]
    del km[off[p]:off[p + 1]]
    del pat[p * w:(p + 1) * w]
    del sup[p]
    del off[p + 1]
    for i in range(p + 1, len(off)):
        off[i] -= n
    s["n_patterns"] -= 1
    s["kmers_kept"] -= n
    s["kmers_distinct"] -= n


def drop_kmer(s, w, pat, sup, off, km):
    p = high_pattern(s, w, pat, sup, off, km, min_len=2)
    del km[off[p + 1] - 1]
    for i in range(p + 1, len(off)):
        off[i] -= 1
    s["kmers_kept"] -= 1
    s["kmers_distinct"] -= 1


def kept_to_rare(s, w, pat, sup, off, km):
    p = high_pattern(s, w, pat, sup, off, km, min_len=2)
    del km[off[p + 1] - 1]
    for i in range(p + 1, len(off)):
        off[i] -= 1
    s["kmers_kept"] -= 1
    s["kmers_rare"] += 1


def fake_kmer(s, w, pat, sup, off, km):
    p = high_pattern(s, w, pat, sup, off, km)
    x = km[off[p + 1] - 1] + 1
    while x in set(km):
        x += 1
    km.insert(off[p + 1], x)
    for i in range(p + 1, len(off)):
        off[i] += 1
    s["kmers_kept"] += 1
    s["kmers_distinct"] += 1


def flip_member(s, w, pat, sup, off, km):
    p = high_pattern(s, w, pat, sup, off, km)
    row = pat[p * w:(p + 1) * w]
    c = next(c for c in range(64 * w) if (row[c // 64] >> (c % 64)) & 1)
    row[c // 64] &= ~(1 << (c % 64))
    sup[p] -= 1
    pat[p * w:(p + 1) * w] = row


def drop_lowest_pattern(s, w, pat, sup, off, km):
    """The omission inside the region every sketch covers."""
    p = min(range(len(sup)), key=lambda p: km[off[p]])
    n = off[p + 1] - off[p]
    del km[off[p]:off[p + 1]]
    del pat[p * w:(p + 1) * w]
    del sup[p]
    del off[p + 1]
    for i in range(p + 1, len(off)):
        off[i] -= n
    s["n_patterns"] -= 1
    s["kmers_kept"] -= n
    s["kmers_distinct"] -= n


OMISSIONS = [drop_pattern, drop_kmer, kept_to_rare]


class Validate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = build_root(Path(cls.tmp.name))
        cls.sp = cls.root / "species" / SLUG
        cls.pristine = {f: (cls.sp / "matrix" / f).read_bytes() for f in FILES}
        cls.n = 0
        ids = (cls.sp / "genomes.txt").read_text().split()
        base = Path(cls.tmp.name)
        picked = set(sorted(ids, key=lambda g: hashlib.sha256(g.encode()).hexdigest())[:35])
        cls.subsets = {"first 40": build_subset(base / "first40", cls.root, ids[:40]),
                       "hash-picked 35": build_subset(base / "picked35", cls.root, [g for g in ids if g in picked])}
        cls.ms3 = build_subset(base / "ms3", cls.root, ids, min_support=3)
        cls.ms3_sub = build_subset(base / "ms3-first40", cls.root, ids[:40], min_support=3)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def validate(self, *extra, root=None):
        sp = (root or self.root) / "species" / SLUG
        Validate.n += 1
        out = Path(self.tmp.name) / f"receipt-{Validate.n}.json"
        p = subprocess.run([sys.executable, "-I", str(VALIDATOR), "--matrix", str(sp / "matrix"),
                            "--counted", str(sp / "counted.json"), "--genomes", str(sp / "genomes.txt"),
                            "--out", str(out), *extra], capture_output=True, text=True)
        return p.returncode, p.stderr, out

    def restore(self):
        for f, b in self.pristine.items():
            (self.sp / "matrix" / f).write_bytes(b)

    @contextlib.contextmanager
    def mutated(self, root, mutate=None, relabel_from=None):
        """`root`'s matrix with `mutate` applied, or replaced by `relabel_from`'s
        matrix relabelled min_support 2; restored afterwards."""
        m = root / "species" / SLUG / "matrix"
        saved = {f: (m / f).read_bytes() for f in FILES}
        try:
            if relabel_from is not None:
                for f in FILES:
                    (m / f).write_bytes((relabel_from / "species" / SLUG / "matrix" / f).read_bytes())
                s = json.loads((m / "matrix.json").read_text())
                s["min_support"] = 2
                (m / "matrix.json").write_text(json.dumps(s, indent=2))
            else:
                s, w, pat, sup, off, km = load(m)
                mutate(s, w, pat, sup, off, km)
                save(m, s, pat, sup, off, km)
            yield
        finally:
            for f, b in saved.items():
                (m / f).write_bytes(b)

    def test_pristine_passes_with_an_independent_recount(self):
        self.restore()
        rc, err, out = self.validate("--fasta-dir", str(self.root / "fasta"), "--recount", "6", "--recount-ids", "G000,G069",
                                     "--require-recount")
        self.assertEqual(rc, 0, err)
        receipt = json.loads(out.read_text())
        self.assertEqual(len(receipt["checks"]["recount"]), 8)
        self.assertGreater(receipt["checks"]["sketch_content"]["kmer_bits_compared"], 0)

    def test_every_review_corruption_fails(self):
        for mutate in MUTANTS:
            with self.subTest(mutate.__name__):
                self.restore()
                m = self.sp / "matrix"
                s, w, pat, sup, off, km = load(m)
                mutate(s, w, pat, sup, off, km)
                save(m, s, pat, sup, off, km)
                rc, err, _ = self.validate()
                self.assertNotEqual(rc, 0, mutate.__name__)
                self.assertIn("INVALID", err)
        self.restore()

    def test_a_fasta_edited_after_counting_fails_the_recount(self):
        self.restore()
        fasta = self.root / "fasta" / "G003.fna.gz"
        original = fasta.read_bytes()
        text = gzip.decompress(original).decode()
        fasta.write_bytes(gzip.compress(text.replace("\n>c2", "A\n>c2", 1).encode()))
        try:
            rc, err, _ = self.validate("--fasta-dir", str(self.root / "fasta"), "--recount-ids", "G003")
        finally:
            fasta.write_bytes(original)
        self.assertNotEqual(rc, 0)
        self.assertIn("G003: recounted", err)

    def test_require_recount_needs_a_fasta_dir(self):
        self.restore()
        rc, err, _ = self.validate("--require-recount")
        self.assertNotEqual(rc, 0)
        self.assertIn("--require-recount", err)

    def test_the_review_omissions_fail_with_every_genome_recounted(self):
        self.restore()
        everyone = ("--fasta-dir", str(self.root / "fasta"), "--recount", "70", "--require-recount")
        rc, err, out = self.validate(*everyone)
        self.assertEqual(rc, 0, err)
        self.assertTrue(json.loads(out.read_text())["checks"]["recount_omissions"]["exact"])
        for mutate in OMISSIONS:
            with self.subTest(mutate.__name__), self.mutated(self.root, mutate):
                rc, err, _ = self.validate(*everyone)
                self.assertNotEqual(rc, 0, mutate.__name__)
                self.assertIn("INVALID", err)
        with self.subTest("built at min_support 3, labelled 2"), self.mutated(self.root, relabel_from=self.ms3):
            rc, err, _ = self.validate(*everyone)
            self.assertNotEqual(rc, 0)

    def test_the_exact_region_catches_omissions_without_a_recount(self):
        self.restore()
        rc, err, out = self.validate()
        self.assertEqual(rc, 0, err)
        region = json.loads(out.read_text())["checks"]["exact_region"]
        self.assertGreater(region["kept"], 0)
        with self.mutated(self.root, drop_lowest_pattern):
            rc, err, _ = self.validate()
            self.assertNotEqual(rc, 0)
            self.assertIn("missing from the matrix", err)
        with self.mutated(self.root, relabel_from=self.ms3):
            rc, err, _ = self.validate()
            self.assertNotEqual(rc, 0)
            self.assertIn("missing from the matrix", err)

    def test_a_sampled_recount_flags_a_kept_kmer_it_can_prove_is_missing(self):
        # Two recounted genomes hold the k-mer and two lack it, so it is
        # neither rare nor core at min_support 2, wherever it lies.
        self.restore()
        s, w, pat, sup, off, km = load(self.sp / "matrix")
        n = len(s["genomes"])
        p = high_pattern(s, w, pat, sup, off, km)
        cols = [c for c in range(n) if (pat[p * w + c // 64] >> (c % 64)) & 1]
        absent = [c for c in range(n) if c not in cols]
        ids = ",".join(s["genomes"][c] for c in cols[:2] + absent[:2])
        with self.mutated(self.root, drop_pattern):
            rc, err, _ = self.validate("--fasta-dir", str(self.root / "fasta"), "--recount-ids", ids)
            self.assertNotEqual(rc, 0)
            self.assertIn("leaves it out", err)

    def test_sub_cohorts_must_equal_the_superset_restricted_to_their_genomes(self):
        self.restore()
        superset = ("--superset", str(self.sp / "matrix"))
        for label, sub in self.subsets.items():
            with self.subTest(label):
                rc, err, out = self.validate(*superset, root=sub)
                self.assertEqual(rc, 0, err)
                got = json.loads(out.read_text())["checks"]["superset"]
                self.assertEqual(got["kmers_matched"], json.loads((sub / "species" / SLUG / "matrix" / "matrix.json").read_text())["kmers_kept"])
            for mutate in OMISSIONS + [fake_kmer, flip_member, drop_lowest_pattern]:
                with self.subTest(label=label, mutant=mutate.__name__), self.mutated(sub, mutate):
                    rc, err, _ = self.validate(*superset, root=sub)
                    self.assertNotEqual(rc, 0, (label, mutate.__name__))
                    self.assertIn("INVALID", err)
        with self.subTest("sub-cohort built at min_support 3, labelled 2"), \
                self.mutated(self.subsets["first 40"], relabel_from=self.ms3_sub):
            rc, err, _ = self.validate(*superset, root=self.subsets["first 40"])
            self.assertNotEqual(rc, 0)
        # An omission in the superset is caught by a sub-cohort that keeps it.
        sub = self.subsets["first 40"]
        sub_kmers = set(load(sub / "species" / SLUG / "matrix")[5])

        def drop_shared(s, w, pat, sup, off, km):
            # The last k-mer of a superset pattern that keeps others and that
            # the sub-cohort keeps too.
            p = next(p for p in range(len(sup)) if off[p + 1] - off[p] > 1 and km[off[p + 1] - 1] in sub_kmers)
            del km[off[p + 1] - 1]
            for q in range(p + 1, len(off)):
                off[q] -= 1
            s["kmers_kept"] -= 1
            s["kmers_distinct"] -= 1
        with self.subTest("superset omission"), self.mutated(self.root, drop_shared):
            rc, err, _ = self.validate(*superset, root=sub)
            self.assertNotEqual(rc, 0)
            self.assertIn("not kept in the superset", err)


if __name__ == "__main__":
    unittest.main()
