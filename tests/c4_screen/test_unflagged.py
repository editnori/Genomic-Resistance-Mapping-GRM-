"""scripts/c4-screen/unflagged.py and unflagged.sh: the check that grm with no C4 flag writes main's bytes
(docs/screen-arms.md, The unflagged bytes).

Every folder, file and genome here is SYNTHETIC: made-up ids, contigs and JSON.
"""

import contextlib
import gzip
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.c4_screen import common

unflagged = common.load("unflagged")
CONTIGS = [("a", "ACGT" * 50), ("b", "GGCA" * 40), ("c", "TTAG" * 30)]


def run(fn, *args):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = fn(*args)
    return rc, buf.getvalue()


def tree(top, files):
    for rel, text in files.items():
        p = Path(top, rel)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)


class Compare(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.a, self.b = Path(self.tmp.name, "a"), Path(self.tmp.name, "b")
        report = {"id": "made_up__drug", "generated_at": "2026-01-01T00:00:00Z", "seconds": 1.5,
                  "splits": [{"seconds_inner": 0.25, "score": 0.5}]}
        for top in (self.a, self.b):
            tree(top, {"out/r.json": json.dumps(report), "out/x.tsv": "a\t1\n", "out/run.log": "took 1 s\n"})

    def test_equal_copies_pass(self):
        rc, text = run(unflagged.compare, self.a, self.b, [])
        self.assertEqual(rc, 0)
        self.assertIn("files: 2 in both, 0 in one only; equal bytes 2, equal but for the run clocks 0, differing 0", text)

    def test_the_run_clocks_alone_may_differ(self):
        doc = json.loads(Path(self.b, "out/r.json").read_text())
        doc["generated_at"], doc["seconds"], doc["splits"][0]["seconds_inner"] = "2026-02-02T00:00:00Z", 9.0, 3.0
        Path(self.b, "out/r.json").write_text(json.dumps(doc))
        Path(self.b, "out/run.log").write_text("took 9 s\n")
        rc, text = run(unflagged.compare, self.a, self.b, [])
        self.assertEqual(rc, 0)
        self.assertIn("equal bytes 1, equal but for the run clocks 1, differing 0", text)

    def test_any_other_change_fails_and_names_the_path_with_ids_masked(self):
        tree(self.b, {"out/x.tsv": "a\t2\n"})
        doc = json.loads(Path(self.b, "out/r.json").read_text())
        doc["splits"][0]["score"] = 0.75
        Path(self.b, "out/r.json").write_text(json.dumps(doc))
        tree(self.a, {"out/573.123.json": "{}"})
        rc, text = run(unflagged.compare, self.a, self.b, [])
        self.assertEqual(rc, 1)
        self.assertIn("differing 2", text)
        self.assertIn("in one only: out/[id].json", text)
        self.assertNotIn("573.123", text)

    def test_a_skipped_prefix_is_left_out(self):
        tree(self.b, {"root/ksets/k.bin": "x"})
        self.assertEqual(run(unflagged.compare, self.a, self.b, [])[0], 1)
        self.assertEqual(run(unflagged.compare, self.a, self.b, ["root/ksets"])[0], 0)


class Fasta(unittest.TestCase):
    def test_a_genome_is_the_same_for_its_id_and_holds_only_test_contigs_and_filler(self):
        one = unflagged.genome(CONTIGS, "1.1")
        self.assertEqual(one, unflagged.genome(CONTIGS, "1.1"))
        self.assertNotEqual(one, unflagged.genome(CONTIGS, "2.2"))
        seqs, cur = [], []
        for line in gzip.decompress(one).decode().splitlines():
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                cur = []
            else:
                cur.append(line)
        seqs.append("".join(cur))
        fillers = [s for s in seqs if s not in {c for _, c in CONTIGS}]
        self.assertEqual([len(s) for s in fillers], [unflagged.FILLER] * 2)
        self.assertTrue(all(set(s) <= set("ACGT") for s in seqs))

    def test_fasta_writes_one_file_per_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            td, ids, out = Path(tmp, "test_dna.fa"), Path(tmp, "genomes.txt"), Path(tmp, "fasta")
            td.write_text("".join(f">{n} made up\n{s}\n" for n, s in CONTIGS))
            ids.write_text("1.1\n2.2\n\n3.3\n")
            rc, text = run(unflagged.fasta, td, ids, out)
            self.assertEqual(rc, 0)
            self.assertEqual(sorted(p.name for p in out.iterdir()), ["1.1.fna.gz", "2.2.fna.gz", "3.3.fna.gz"])
            self.assertEqual(text, "synthetic genomes: 3, from 3 test_dna.fa contigs\n")


class Script(unittest.TestCase):
    def test_the_check_passes_no_c4_flag_and_runs_without_a_network(self):
        text = (common.KIT / "unflagged.sh").read_text()
        for flag in ("--exp-screen", "--exp-floor", "--screened-features"):
            self.assertNotIn(flag, text)
        calls = [line for line in text.splitlines() if '"$bin" --data' in line]
        self.assertEqual(len(calls), 4)
        self.assertTrue(all("unshare -rn" in line for line in calls))


if __name__ == "__main__":
    unittest.main()
