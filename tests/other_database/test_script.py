"""scripts/other-database/record.py on SYNTHETIC inputs: cells worked by hand, same bytes, stops and --check.

The inputs are built in a temporary directory by synthetic.make; none of them is data.
"""

import contextlib
import hashlib
import io
import math
import os
import pathlib
import subprocess
import sys
import tempfile
import unittest

from tests.other_database import synthetic
from tests.other_database.common import SCRIPT, load_script, read_tsv

record = load_script()


def run(tree, *args):
    """main() on a synthetic tree; returns (rc, stderr)."""
    lanes, release, data = tree
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        rc = record.main(["--lanes", lanes, "--release", release, "--data-root", data, *args])
    return rc, err.getvalue()


def files(directory):
    return {p.relative_to(directory).as_posix(): p.read_bytes() for p in sorted(pathlib.Path(directory).rglob("*")) if p.is_file()}


class SyntheticRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = cls.tmp.name
        cls.tree = synthetic.make(os.path.join(cls.root, "in"))
        cls.out = os.path.join(cls.root, "out")
        rc, err = run(cls.tree, "--out", cls.out)
        assert rc == 0, err
        cls.ds = {(r["dataset"], r["learner"]): r for r in read_tsv(os.path.join(cls.out, "datasets.tsv"))}
        cls.rd = {(r["dataset"], r["learner"]): r for r in read_tsv(os.path.join(cls.out, "readings.tsv"))}
        cls.sp = {r["species"]: r for r in read_tsv(os.path.join(cls.out, "species.tsv"))}
        cls.order = {r["species"]: r for r in read_tsv(os.path.join(cls.out, "order.tsv"))}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_scm_cells_by_hand(self):
        # SYNTHETIC: 27 of 30 R called resistant and 28 of 30 S called susceptible.
        r = self.ds[("salmonella_enterica__synthomycin", "scm")]
        self.assertEqual([r[k] for k in ("R", "S", "TP", "FN", "TN", "FP")], ["30", "30", "27", "3", "28", "2"])
        self.assertEqual(r["recall"], f"{27 / 30:.4f}")
        self.assertEqual(r["specificity"], f"{28 / 30:.4f}")
        # 5 errors at the clip 1e-4 and 55 right calls at 1 - 1e-4.
        self.assertEqual(r["log_loss"], f"{(5 * -math.log(1e-4) + 55 * -math.log(1 - 1e-4)) / 60:.4f}")
        self.assertEqual(r["log_loss"], "0.7676")
        # Clopper-Pearson of 27 of 30, as tables give it.
        self.assertEqual((r["recall_lo"], r["recall_hi"]), ("0.7347", "0.9789"))
        # The lineage split: 16 of 20 R and 18 of 20 S, 6 errors of 40.
        self.assertEqual((r["lineage_recall"], r["lineage_specificity"]), ("0.8000", "0.9000"))
        self.assertEqual(r["lineage_log_loss"], f"{(6 * -math.log(1e-4) + 34 * -math.log(1 - 1e-4)) / 40:.4f}")
        self.assertEqual((r["recall_reading"], r["log_loss_reading"], r["holds"], r["estimable"]), ("consistent", "consistent", "yes", "yes"))

    def test_cart_cells_by_hand(self):
        r = self.ds[("salmonella_enterica__synthomycin", "cart")]
        # Every genome right at 0.9 or 0.1: every loss is -ln 0.9, and the bootstrap cannot move.
        self.assertEqual((r["log_loss"], r["log_loss_lo"], r["log_loss_hi"]), (f"{-math.log(0.9):.4f}",) * 3)
        # 30 of 30: the Clopper-Pearson lower bound is 0.025 ** (1 / 30).
        self.assertEqual((r["recall_lo"], r["recall_hi"]), (f"{0.025 ** (1 / 30):.4f}", "1.0000"))
        self.assertEqual(r["lineage_log_loss"], f"{(-math.log(0.8) - math.log(0.7)) / 2:.4f}")
        self.assertEqual((r["recall_reading"], r["log_loss_reading"]), ("consistent", "better"))

    def test_labels_consolidated(self):
        # SYNTHETIC pilot genomes: R and S (conflict), intermediate only, "not defined" only, and a drug no release holds.
        r = self.ds[("escherichia_coli__synthomycin", "scm")]
        self.assertEqual([r[k] for k in ("scored", "labels_R", "labels_S", "labels_I_only", "labels_no_call", "labels_conflict")],
                         ["64", "30", "30", "1", "2", "1"])
        # The NCBI name reaches the release name through the plan's fixed table; one BioSample is lower case.
        r = self.ds[("escherichia_coli__amoxicillin_clavulanate", "scm")]
        self.assertEqual([r[k] for k in ("labels_R", "labels_S", "labels_no_call")], ["30", "30", "4"])

    def test_too_few_is_not_estimable(self):
        r = self.ds[("escherichia_coli__synthomycin", "lasso")]
        self.assertEqual((r["R"], r["S"], r["null_calls"], r["estimable"]), ("10", "10", "40", "no"))
        self.assertEqual(r["note"], "fewer than 25 R or 25 S")
        self.assertEqual(r["holds"], "")
        self.assertNotIn(("escherichia_coli__synthomycin", "lasso"), self.rd)

    def test_tie_goes_to_lower_accession(self):
        # SYNTHETIC: PRJNA9 and PRJNA10 hold 30 genomes each; as text PRJNA10 would come first.
        self.assertEqual(self.ds[("escherichia_coli__synthomycin", "scm")]["largest_bioproject"], "PRJNA9")
        self.assertEqual(self.sp["escherichia_coli"]["largest_bioproject"], "PRJNA9")
        self.assertEqual(record.largest({"PRJNA1160997": 2, "PRJNA717739": 2}), "PRJNA717739")

    def test_species_row_by_hand(self):
        s = self.sp["salmonella_enterica"]
        self.assertEqual([s[k] for k in ("C1", "C3", "C4", "P2", "D3", "scored")], ["96", "10", "20", "5", "1", "60"])
        self.assertEqual((s["largest_bioproject"], s["largest_bioproject_scored"], s["largest_bioproject_share"]), ("PRJNA100", "24", "0.400"))
        self.assertEqual((s["shared_bioproject"], s["largest_bioproject_training"], s["training_no_bioproject"]), ("24", "2", "1"))
        self.assertEqual((s["bvbrc_no_accession"], s["no_accession_share"]), ("25", "0.250"))
        # Distances 0.0006 to 0.0065 by 0.0001: 5 up to 0.001, 40 up to 0.005, 15 above.
        self.assertEqual([s[k] for k in ("dist_0.0005_0.001", "dist_0.001_0.005", "dist_above_0.005")], ["5", "40", "15"])
        self.assertEqual(self.sp["escherichia_coli"]["O1"], "")

    def test_reading_L_interval(self):
        # The release interval of the SYNTHETIC pilot dataset: vme 0.05 to 0.40, recall 0.60 to 0.95.
        r = self.rd[("escherichia_coli__synthomycin", "scm")]
        self.assertEqual((r["L_lineage_recall_lo"], r["L_lineage_recall_hi"], r["L_lineage_recall_source"]), ("0.6000", "0.9500", "release interval"))
        self.assertEqual(r["status"], "exploratory")
        r = self.rd[("salmonella_enterica__synthomycin", "scm")]
        self.assertEqual(r["L_lineage_recall_source"], "released scores")
        self.assertEqual((r["L_lineage_recall_lo"], r["L_lineage_recall_hi"]), tuple(f"{v:.4f}" for v in record.clopper_pearson(16, 20)))
        self.assertEqual((r["status"], r["errors"]), ("fixed in advance", "5"))

    def test_order_rows(self):
        p = self.order["escherichia_coli"]
        # Two SYNTHETIC requests before PLAN.md with no parameters; one after it is not counted.
        self.assertEqual((p["unrecorded_before_plan"], p["unrecorded_bytes"], p["unrecorded_first"]), ("2", "67;304", "2025-12-31T23:00:00Z"))
        self.assertEqual((p["assigned"], p["label_requests"], p["first_at_or_after_predictions"]), ("by time", "1", "yes"))
        s = self.order["salmonella_enterica"]
        self.assertEqual((s["biosamples_asked"], s["scored_biosamples"], s["asked_equals_scored"], s["unrecorded_bytes"]), ("60", "60", "yes", ""))

    def test_inputs_named_by_root(self):
        text = pathlib.Path(self.out, "inputs.sha256").read_text(encoding="utf-8")
        self.assertNotIn(self.root, text)
        self.assertTrue(all(line.split("  ")[1].split("/")[0] in ("lanes", "release", "data-root") for line in text.splitlines()))

    def test_same_bytes_on_a_second_run(self):
        again = os.path.join(self.root, "again")
        rc, err = run(self.tree, "--out", again, "--jobs", "2")
        self.assertEqual(rc, 0, err)
        self.assertEqual(files(again), files(self.out))

    def test_same_bytes_at_other_hash_seeds(self):
        lanes, release, data = self.tree
        for seed in ("0", "4242"):
            out = os.path.join(self.root, "seed" + seed)
            subprocess.run([sys.executable, "-I", "-B", str(SCRIPT), "--lanes", lanes, "--release", release, "--data-root", data, "--out", out],
                           env={"PYTHONHASHSEED": seed}, check=True, capture_output=True)
            self.assertEqual(files(out), files(self.out), seed)

    def test_check(self):
        self.assertEqual(run(self.tree, "--check", self.out), (0, ""))
        work = os.path.join(self.root, "check")
        rc, _ = run(self.tree, "--out", work)
        self.assertEqual(rc, 0)
        totals = pathlib.Path(work, "totals.tsv")
        totals.write_text(totals.read_text(encoding="utf-8").replace("\t8\t", "\t7\t", 1), encoding="utf-8")
        pathlib.Path(work, "plan", "PLAN.md").unlink()
        pathlib.Path(work, "extra.tsv").write_text("x\n", encoding="utf-8")
        rc, err = run(self.tree, "--check", work)
        self.assertEqual(rc, 1)
        self.assertIn("totals.tsv: differs", err)
        self.assertIn("plan/PLAN.md: missing", err)
        self.assertIn("extra.tsv: not written by the script", err)
        self.assertEqual(run(self.tree, "--check", os.path.join(self.root, "none"))[0], 1)

    def test_sums_whole(self):
        sums = pathlib.Path(self.out, "SHA256SUMS").read_text(encoding="utf-8").splitlines()
        named = {line.split("  ")[1]: line.split("  ")[0] for line in sums}
        on_disk = files(self.out)
        self.assertEqual(set(named), set(on_disk) - {"SHA256SUMS"})
        for name, digest in named.items():
            self.assertEqual(hashlib.sha256(on_disk[name]).hexdigest(), digest, name)


class Stops(unittest.TestCase):
    """Each SYNTHETIC fault stops the script with exit 2 and writes nothing."""

    FAULTS = {
        "pin": "differs from the lane's pin",
        "funnel": "funnel step P2 does not add up",
        "missing": "missing file",
        "d3_count": "not the genomes the filter kept",
        "stray_genome": "which is not selected",
        "plan_sha": "its sha256 differs",
        "lineage_metric": "differs from the report's",
    }

    def test_faults_stop(self):
        for fault, message in self.FAULTS.items():
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as tmp:
                tree = synthetic.make(os.path.join(tmp, "in"), fault=fault)
                out = os.path.join(tmp, "out")
                rc, err = run(tree, "--out", out)
                self.assertEqual(rc, 2)
                self.assertIn("stopped:", err)
                self.assertIn(message, err)
                self.assertFalse(os.path.exists(out))


class Pieces(unittest.TestCase):
    def test_ranks(self):
        # Sorted values 50 and 1950 of 2000; (1 - 0.95) / 2 * 2000 in floating point would give 51.
        self.assertEqual(record.ranks(2000), (49, 1949))
        self.assertEqual(record.ranks(40), (0, 38))

    def test_clopper_pearson_ends(self):
        for n in (10, 30, 64):
            lo, hi = record.clopper_pearson(n, n)
            self.assertAlmostEqual(lo, 0.025 ** (1 / n), places=9)
            self.assertEqual(hi, 1.0)
            lo, hi = record.clopper_pearson(0, n)
            self.assertEqual(lo, 0.0)
            self.assertAlmostEqual(hi, 1 - 0.025 ** (1 / n), places=9)

    def test_bootstrap_of_constant_losses(self):
        self.assertEqual(record.bootstrap([0.25] * 30), (0.25, 0.25))

    def test_rule_and_margin(self):
        # meets: holds - floor(e / 2); fails: floor(e / 2) + 1 - holds; under 3 datasets nothing.
        cases = [(5, 8, "meets", "1"), (4, 8, "fails", "1"), (8, 8, "meets", "4"), (0, 7, "fails", "4"), (2, 3, "meets", "1"), (2, 2, "not assessable", "")]
        for holds, est, result, margin in cases:
            self.assertEqual((record.rule(holds, est), record.margin(holds, est)), (result, margin), (holds, est))

    def test_loss_clip(self):
        # p is held to 1e-4 and 1 - 1e-4.
        self.assertEqual(record.loss(1, 1.0), -math.log(1 - 1e-4))
        self.assertEqual(record.loss(1, 0.0), -math.log(1e-4))
        self.assertAlmostEqual(record.loss(0, 1.0), -math.log(1e-4), places=9)
        self.assertEqual(record.loss(1, 0.5), -math.log(0.5))

    def test_parse_and_consolidate(self):
        self.assertEqual([record.parse_call(x) for x in ("Resistant", "nonsusceptible", "susceptible", "intermediate", "not defined", None)],
                         ["R", "R", "S", "I", None, None])
        self.assertEqual([record.consolidate(c) for c in ({"R", "S"}, {"R", "I"}, {"S"}, {"I"}, set())], ["conflict", "R", "S", "I only", "no call"])


if __name__ == "__main__":
    unittest.main()
