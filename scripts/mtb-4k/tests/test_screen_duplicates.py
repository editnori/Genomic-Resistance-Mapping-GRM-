"""screen_duplicates.py on synthetic cohorts.

usage: python3 -I scripts/mtb-4k/tests/test_screen_duplicates.py
"""
import json
import random
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPT = HERE.parent / "screen_duplicates.py"
sys.path.insert(0, str(HERE.parent))
import screen_duplicates as sd  # noqa: E402

SLUG = "mycobacterium_tuberculosis"
SPACE = 1 << 62


def kset(rng, core, private):
    return sorted(set(core) | {rng.randrange(SPACE) for _ in range(private)})


def write_cohort(base, name, genomes):
    sp = Path(base) / name / "species" / SLUG
    sp.mkdir(parents=True)
    (sp / "genomes.txt").write_text("".join(f"{g}\n" for g, _ in genomes))
    counted = {"k": 31, "genomes": [{"genome_id": g, "n_kmers": len(v), "kset_digest": f"d-{hash(tuple(v)) & 0xffffffff:08x}",
                                     "sketch": v[:1000]} for g, v in genomes],
               "missing_fasta": [], "failed": []}
    (sp / "counted.json").write_text(json.dumps(counted))
    return Path(base) / name


class Screen(unittest.TestCase):
    def test_finds_exact_and_near_duplicates_across_cohorts(self):
        rng = random.Random(7)
        core = [rng.randrange(SPACE) for _ in range(20000)]
        base = kset(rng, core, 3000)
        near = sorted(set(base[:-1]) | {base[0] - 1 if base[0] else base[0] + 1})  # one value differs
        prefix = [("P0", base), ("P1", kset(rng, core, 3000)), ("P2", kset(rng, core, 3000))]
        added = [("A0", list(base)), ("A1", near), ("A2", kset(rng, core, 3000)), ("A3", kset(rng, core, 3000))]
        with tempfile.TemporaryDirectory() as tmp:
            p = write_cohort(tmp, "prefix", prefix)
            a = write_cohort(tmp, "added", added)
            manifest = Path(tmp) / "fasta-manifest.tsv"
            manifest.write_text("genome_id\tsha256\nP1\tsame\nA2\tsame\nP2\tx2\nA3\tx3\n")
            records = Path(tmp) / "genomes.jsonl"
            records.write_text(json.dumps({"genome_id": "P2", "biosample_accession": "SAMN1"}) + "\n" +
                               json.dumps({"genome_id": "A3", "biosample_accession": "SAMN1"}) + "\n" +
                               json.dumps({"genome_id": "P0", "biosample_accession": None}) + "\n")
            out = Path(tmp) / "screen.json"
            subprocess.run([sys.executable, "-I", str(SCRIPT), "--prefix-root", str(p), "--added-root", str(a),
                            "--fasta-manifest", str(manifest), "--genomes-jsonl", str(records), "--out", str(out)],
                           check=True, capture_output=True)
            res = json.loads(out.read_text())
        kinds = sorted((e["kind"], sorted(m["id"] for m in e["genomes"])) for e in res["exact"])
        self.assertEqual(kinds, [("biosample_accession", ["A3", "P2"]), ("fasta_sha256", ["A2", "P1"]),
                                 ("kset_digest", ["A0", "P0"])])
        self.assertEqual(res["cross_cohort_exact"], 3)
        near_pairs = sorted(tuple(sorted((x["a"]["id"], x["b"]["id"]))) for x in res["near"])
        self.assertIn(("A0", "P0"), near_pairs)
        self.assertIn(("A1", "P0"), near_pairs)
        self.assertIn(("A0", "A1"), near_pairs)
        self.assertNotIn(("A2", "P1"), near_pairs)  # same FASTA hash, different sketches: exact list only

    def test_mash_jaccard(self):
        a = list(range(0, 2000, 2))
        self.assertEqual(sd.mash_jaccard(a, a), 1.0)
        # b swaps a's largest value for 1: of the union's 1,000 smallest
        # values (0, 1 and the evens to 1996), 999 are in both.
        b = sorted(set(a[:-1]) | {1})
        self.assertAlmostEqual(sd.mash_jaccard(a, b), 999 / 1000)


if __name__ == "__main__":
    unittest.main()
