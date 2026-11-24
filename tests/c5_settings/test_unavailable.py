"""SYNTHETIC tests of scripts/c5-settings/score.py at a default draw that came back unavailable, and of the
scorer against the run's own, at commit 0d987ba, where no default draw is unavailable.

Every input is made up here: one species of two datasets, as test_score.write_run writes it.
"""

import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.c5_settings import common
from tests.c5_settings import test_score as ts

S, E = ts.S, ts.E
REASON = "invalid lineage outer evaluation: fold 3 validation lacks one or both classes"
FIRST = "the default's result at the same dataset, split and seed is there and available"
BESIDE = "the default's result is unavailable and no choices file stands beside it"
# The run's scorer, as commit 0d987ba holds it, and its git blob.
AT_RUN = Path(__file__).resolve().parent / "score_at_0d987ba.txt"
AT_RUN_BLOB = "6ad3de887608fb73714812573c6ed2fd390d0c96"
TABLES = ("checks.tsv", "draws.tsv", "datasets.tsv")


def rows(top):
    with open(Path(top) / "expect.tsv", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def file_of(top, arm, split, seed, dataset=common.OTHER):
    return [e["file"] for e in rows(top) if (e["dataset"], e["arm"], e["split"], e["seed"])
            == (dataset, arm, split, str(seed))][0]


def unavailable(top, name):
    """Make the result at name an unavailable one, with no outer fold, as grm learn writes it."""
    path = Path(top) / "picks" / name
    doc = json.loads(path.read_text())
    doc["unavailable_reason"], doc["folds"] = REASON, []
    path.write_text(json.dumps(doc))


def unavailable_run(top, choices_kept=False):
    """A run whose default draw on the second dataset's lineage split at seed 43 came back unavailable, with no
    choices file (or with fold 0's kept), and the twin's draw at seed 46 unavailable for the same reason."""
    args = ts.write_run(top)
    unavailable(top, file_of(top, "default", "lineage", 43))
    unavailable(top, file_of(top, "twin", "lineage", 46))
    for k in range(S.FOLDS):
        p = Path(top) / "picks" / E.choices_name(common.OTHER, "lineage", 43, k)
        if p.exists() and not (choices_kept and k == 0):
            p.unlink()
    return args


class Unavailable(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_arm_rows_beside_an_unavailable_default_are_not_scored_for_its_reason(self):
        rc, text = ts.score(unavailable_run(self.top))
        self.assertEqual(rc, 0, text)
        self.assertIn("missing draws 0", text)
        checks = ts.tsv(self.top / "score" / "checks.tsv")
        self.assertTrue(all(c["passed"] == "1" for c in checks))
        beside = [c for c in checks if c["check"] == BESIDE]
        self.assertEqual([(c["file"], c["detail"]) for c in beside],
                         [(file_of(self.top, "search", "lineage", 43), REASON)])
        draws = ts.tsv(self.top / "score" / "draws.tsv")
        at = [r for r in draws if (r["dataset"], r["split"], r["seed"]) == (common.OTHER, "lineage", "43")]
        self.assertEqual(len(at), 4 * 5)
        for r in at:
            self.assertEqual((r["scored"], r["constant"], r["why"], r["balanced_accuracy"]),
                             ("0", "0", "split unavailable", ""), r["arm"])
        twin = [r for r in draws if (r["dataset"], r["arm"], r["seed"]) == (common.OTHER, "twin", "46")]
        self.assertEqual({(r["scored"], r["why"]) for r in twin}, {("0", "split unavailable")})
        others = [r for r in draws if (r["dataset"], r["split"]) == (common.OTHER, "lineage")
                  and r["seed"] in ("42", "44") and r["arm"] in ("default", "search", "cutoff", "plain-best")]
        self.assertEqual(len(others), 2 * 4 * 5)
        self.assertTrue(all(r["scored"] == "1" and r["balanced_accuracy"] != "" for r in others))
        for r in ts.tsv(self.top / "score" / "datasets.tsv"):
            if (r["dataset"], r["split"]) == (common.OTHER, "lineage"):
                self.assertEqual((r["scored"], r["balanced_accuracy"]), ("0", ""), r["arm"])
            else:
                self.assertEqual(r["scored"], "1", (r["dataset"], r["split"], r["arm"]))

    def test_a_choices_file_beside_an_unavailable_default_fails(self):
        rc, text = ts.score(unavailable_run(self.top, choices_kept=True))
        self.assertEqual(rc, 1, text)
        failed = [c for c in ts.tsv(self.top / "score" / "checks.tsv") if c["passed"] != "1"]
        self.assertEqual([(c["file"], c["check"]) for c in failed],
                         [(file_of(self.top, "search", "lineage", 43), BESIDE)])
        self.assertFalse((self.top / "score" / "draws.tsv").exists())

    def test_an_absent_default_still_fails_the_first_choices_check(self):
        args = ts.write_run(self.top)
        (self.top / "picks" / file_of(self.top, "default", "lineage", 43)).unlink()
        rc, text = ts.score(args)
        self.assertEqual(rc, 1, text)
        failed = [c for c in ts.tsv(self.top / "score" / "checks.tsv") if c["passed"] != "1"]
        self.assertEqual([(c["file"], c["check"]) for c in failed], [(file_of(self.top, "search", "lineage", 43), FIRST)])

    def test_an_available_default_of_one_fold_still_fails_the_first_choices_check(self):
        args = ts.write_run(self.top)
        path = self.top / "picks" / file_of(self.top, "default", "lineage", 43)
        doc = json.loads(path.read_text())
        doc["folds"] = doc["folds"][:1]
        path.write_text(json.dumps(doc))
        rc, text = ts.score(args)
        self.assertEqual(rc, 1, text)
        failed = {c["check"] for c in ts.tsv(self.top / "score" / "checks.tsv") if c["passed"] != "1"}
        self.assertIn(FIRST, failed)


def run_scorer(script, args, out):
    run = subprocess.run([sys.executable, "-I", "-B", str(script)] + list(args[:3]) + [str(out)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    return run.returncode, run.stdout


class AtTheRun(unittest.TestCase):
    """The changed scorer against the run's own, on made-up species with no unavailable default draw."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)
        data = AT_RUN.read_bytes()
        self.assertEqual(hashlib.sha1(b"blob %d\0" % len(data) + data).hexdigest(), AT_RUN_BLOB)
        # The run's scorer loads C3's from beside its own folder, as the checkout's does.
        old = self.top / "old" / "scripts" / "c5-settings"
        old.mkdir(parents=True)
        (old / "score.py").write_bytes(data)
        os.symlink(common.ROOT / "scripts" / "c3-inputs", self.top / "old" / "scripts" / "c3-inputs")
        self.old = old / "score.py"

    def tearDown(self):
        self.tmp.cleanup()

    def same(self, run_dir, **kw):
        run_dir.mkdir()
        args = ts.write_run(run_dir, **kw)
        rc_old, text_old = run_scorer(self.old, args, run_dir / "old")
        rc_new, text_new = run_scorer(common.KIT / "score.py", args, run_dir / "new")
        self.assertEqual((rc_old, rc_new), (0, 0), text_old + text_new)
        self.assertEqual(text_old, text_new)
        for t in TABLES:
            with self.subTest(table=t):
                self.assertEqual((run_dir / "old" / t).read_bytes(), (run_dir / "new" / t).read_bytes())

    def test_tables_equal_the_run_scorers_byte_for_byte_on_a_whole_run(self):
        self.same(self.top / "whole")

    def test_tables_equal_the_run_scorers_byte_for_byte_with_a_missing_draw(self):
        self.same(self.top / "missing", drop={E.choices_name(common.OTHER, "lineage", 44, 1)})


if __name__ == "__main__":
    unittest.main()
