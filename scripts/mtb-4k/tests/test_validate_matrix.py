"""validate_matrix.py against a matrix the real GRM binary builds.

usage: GRM=target/release/grm python3 -I scripts/mtb-4k/tests/test_validate_matrix.py

The pristine matrix must pass, including an independent FASTA recount that
reproduces GRM's k-mer count, FNV kset digest and sketch; each of the
review's eight corruptions must fail, as must a recount of a FASTA edited
after counting.
"""
import gzip
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
    subprocess.run([grm(), "--data", str(root), "count", "--species", SLUG, "--allow-missing-catalog"], check=True, capture_output=True)
    subprocess.run([grm(), "--data", str(root), "matrix", "--species", SLUG, "--min-support", "2"], check=True, capture_output=True)
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


class Validate(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = build_root(Path(cls.tmp.name))
        cls.sp = cls.root / "species" / SLUG
        cls.pristine = {f: (cls.sp / "matrix" / f).read_bytes() for f in FILES}
        cls.n = 0

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def validate(self, *extra):
        Validate.n += 1
        out = Path(self.tmp.name) / f"receipt-{Validate.n}.json"
        p = subprocess.run([sys.executable, "-I", str(VALIDATOR), "--matrix", str(self.sp / "matrix"),
                            "--counted", str(self.sp / "counted.json"), "--genomes", str(self.sp / "genomes.txt"),
                            "--out", str(out), *extra], capture_output=True, text=True)
        return p.returncode, p.stderr, out

    def restore(self):
        for f, b in self.pristine.items():
            (self.sp / "matrix" / f).write_bytes(b)

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


if __name__ == "__main__":
    unittest.main()
