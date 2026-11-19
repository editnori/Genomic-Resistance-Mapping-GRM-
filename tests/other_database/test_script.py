"""scripts/other-database/record.py on SYNTHETIC inputs: cells worked by hand, same bytes, stops and --check.

The inputs are built in a temporary directory by synthetic.make; none of them is data.
"""

import contextlib
import hashlib
import io
import math
import os
import pathlib
import random
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
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = cls.tmp.name
        cls.tree = synthetic.make(os.path.join(cls.root, "in"))
        cls.out = os.path.join(cls.root, "out")
        rc, err = run(cls.tree, "--out", cls.out)
        assert rc == 0, err
        cls.ds = {(r["dataset"], r["learner"]): r for r in read_tsv(os.path.join(cls.out, "datasets.tsv"))}
        cls.rd = {(r["dataset"], r["learner"]): r for r in read_tsv(os.path.join(cls.out, "readings.tsv"))}
        cls.sp = {r["species"]: r for r in read_tsv(os.path.join(cls.out, "species.tsv"))}
        cls.order = {r["species"]: r for r in read_tsv(os.path.join(cls.out, "order.tsv"))}

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
        # SYNTHETIC: every genome right, R at p 0.50 to 0.79 (0.50 is called resistant) and S at p 0.01 to 0.30.
        losses = [-math.log(0.5 + 0.01 * i) for i in range(30)] + [-math.log(1 - 0.01 * j) for j in range(1, 31)]
        self.assertEqual(r["log_loss"], f"{sum(losses) / 60:.4f}")
        self.assertEqual([r[k] for k in ("TP", "FN", "TN", "FP", "specificity")], ["30", "0", "30", "0", "1.0000"])
        # 30 of 30: the Clopper-Pearson lower bound is 0.025 ** (1 / 30).
        self.assertEqual((r["recall_lo"], r["recall_hi"]), (f"{0.025 ** (1 / 30):.4f}", "1.0000"))
        # The lineage split: 19 R at 0.8, 1 R at 0.5 and 20 S at 0.3; recall 1.0 only if 0.5 counts as resistant.
        self.assertEqual(r["lineage_recall"], "1.0000")
        self.assertEqual(r["lineage_log_loss"], f"{(-19 * math.log(0.8) - math.log(0.5) - 20 * math.log(0.7)) / 40:.4f}")
        # Log loss 0.3107 with its interval 0.2671 to 0.3555 around the lineage figure 0.3017.
        self.assertLess(float(r["log_loss_lo"]), float(r["lineage_log_loss"]))
        self.assertLess(float(r["lineage_log_loss"]), float(r["log_loss_hi"]))
        self.assertEqual((r["recall_reading"], r["log_loss_reading"], r["holds"]), ("consistent", "consistent", "yes"))

    def test_labels_consolidated(self):
        # SYNTHETIC pilot genomes: R and S (conflict), intermediate only, "not defined" only, and a drug no release holds.
        r = self.ds[("escherichia_coli__synthomycin", "scm")]
        self.assertEqual([r[k] for k in ("scored", "labels_R", "labels_S", "labels_I_only", "labels_no_call", "labels_conflict")],
                         ["64", "30", "30", "1", "2", "1"])
        # The NCBI name reaches the release name through the plan's fixed table; one BioSample is lower case.
        r = self.ds[("escherichia_coli__amoxicillin_clavulanate", "scm")]
        self.assertEqual([r[k] for k in ("labels_R", "labels_S", "labels_no_call")], ["30", "30", "4"])

    def test_class_line_at_25(self):
        # SYNTHETIC S. aureus: edgemycin has exactly 25 R and 30 S, underedgemycin 24 R and 30 S.
        for ln in ("scm", "cart"):
            r = self.ds[("staphylococcus_aureus__edgemycin", ln)]
            self.assertEqual((r["R"], r["S"], r["estimable"], r["note"]), ("25", "30", "yes", ""))
            self.assertIn(("staphylococcus_aureus__edgemycin", ln), self.rd)
            r = self.ds[("staphylococcus_aureus__underedgemycin", ln)]
            self.assertEqual((r["R"], r["S"], r["estimable"], r["note"]), ("24", "30", "no", "fewer than 25 R or 25 S"))
            self.assertNotIn(("staphylococcus_aureus__underedgemycin", ln), self.rd)

    def test_each_intermediate_wording_gives_no_label(self):
        # SYNTHETIC edgemycin: genomes 25 to 29 are labelled only by the wordings of synthetic.INTERMEDIATE, in turn.
        r = self.ds[("staphylococcus_aureus__edgemycin", "scm")]
        self.assertEqual([r[k] for k in ("scored", "labels_R", "labels_S", "labels_I_only", "labels_no_call", "labels_conflict")],
                         ["60", "25", "30", "5", "0", "0"])
        for wording in synthetic.INTERMEDIATE:
            self.assertEqual(record.parse_call(wording), "I", wording)
            self.assertEqual(record.parse_call(wording.title()), "I", wording)
        r = self.ds[("staphylococcus_aureus__underedgemycin", "scm")]
        self.assertEqual((r["labels_I_only"], r["labels_no_call"]), ("0", "6"))

    def test_lineage_at_share_by_hand(self):
        clip, keep = -math.log(1e-4), -math.log(1 - 1e-4)
        # SCM: the lineage split's recall 0.8 and specificity 0.9 at the test's share 30 of 60.
        r = self.ds[("salmonella_enterica__synthomycin", "scm")]
        self.assertEqual(r["lineage_error_rate_at_share"], f"{0.5 * 0.2 + 0.5 * 0.1:.4f}")
        self.assertEqual(r["lineage_log_loss_at_share"], f"{0.15 * clip + 0.85 * keep:.4f}")
        # The same at edgemycin's share, 25 of 55.
        r = self.ds[("staphylococcus_aureus__edgemycin", "scm")]
        self.assertEqual(r["lineage_error_rate_at_share"], f"{25 / 55 * 0.2 + 30 / 55 * 0.1:.4f}")
        # CART: no lineage error; the loss of its 20 R (19 at 0.8, 1 at 0.5) and 20 S (at 0.3), weighted half and half.
        r = self.ds[("salmonella_enterica__synthomycin", "cart")]
        res = (-19 * math.log(0.8) - math.log(0.5)) / 20
        self.assertEqual(r["lineage_error_rate_at_share"], "0.0000")
        self.assertEqual(r["lineage_log_loss_at_share"], f"{0.5 * res + 0.5 * -math.log(0.7):.4f}")
        # Blank where the dataset is not estimable.
        r = self.ds[("escherichia_coli__synthomycin", "lasso")]
        self.assertEqual((r["lineage_error_rate_at_share"], r["lineage_log_loss_at_share"]), ("", ""))

    def test_label_request_in_the_plans_second_is_not_after_it(self):
        # SYNTHETIC A. baumannii: predictions written and first label request both at the second ADDENDUM2.md was fixed.
        a = self.order[synthetic.SAME_SECOND]
        self.assertEqual((a["predictions_written"], a["first_label_request"]), (synthetic.ADD2_TIME, synthetic.ADD2_TIME))
        self.assertEqual((a["first_at_or_after_predictions"], a["first_after_plan_file"]), ("yes", "no"))
        s = self.order["salmonella_enterica"]
        self.assertEqual((s["first_at_or_after_predictions"], s["first_after_plan_file"]), ("yes", "yes"))

    def test_seeds_file_on_synthetic_inputs(self):
        self.assertEqual(pathlib.Path(self.out, "seeds.tsv").read_text(encoding="utf-8"), "\t".join(record.SEED_COLS) + "\n")

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
        "log_order": "the label step is not written after the predictions",
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

    def test_bootstrap_seeds(self):
        x = [0.11, 0.23, 0.37, 0.52, 0.68, 0.85, 1.3, 0.04, 0.09, 0.16, 0.27, 2.1]
        self.assertEqual((record.SEED, record.OTHER_SEEDS), (1, (2, 3, 4, 5, 6)))
        self.assertEqual(record.bootstrap(x), record.bootstrap(x, 1))
        self.assertEqual(record.bootstrap(x), (0.2733333333333334, 0.915))
        self.assertEqual(record.bootstrap(x, 2), (0.26916666666666667, 0.93))

    def test_reading_C_draws_resistant_then_susceptible(self):
        res, sus = [0.11, 0.23, 0.37, 0.52, 0.68, 0.85, 1.3], [0.04, 0.09, 0.16, 0.27, 2.1]
        # Written in by hand from a run of the plan's code; drawing the susceptible genomes first gives 0.2577 to 0.9746.
        self.assertEqual(record.stratified(res, sus), (0.2648571428571429, 0.9787142857142858))
        rng = random.Random(1)
        vals = []
        for _ in range(2000):
            r = sum(rng.choices(res, k=len(res))) / len(res)
            s = sum(rng.choices(sus, k=len(sus))) / len(sus)
            vals.append((r + s) / 2)
        vals.sort()
        self.assertEqual(record.stratified(res, sus), (vals[49], vals[1949]))

    def test_seed_rows(self):
        def h(dataset, holds):
            return {"species": "a", "part": "eight", "dataset": dataset, "learner": "scm", "holds": holds}

        flip = ["no", "yes", "yes", "no", "yes", "yes"]
        rows = record.seed_rows([h("a__d1", flip), h("a__d2", ["yes"] * 6), h("a__d3", ["no"] * 6)])
        seeds = [f"seed_{s}" for s in record.SEED_NAMES]
        self.assertEqual(rows, [
            {"level": "dataset", "species": "a", "dataset": "a__d1", "learner": "scm", "part": "eight", **dict(zip(seeds, flip))},
            {"level": "species", "species": "a", "dataset": "", "learner": "scm", "part": "eight",
             **dict(zip(seeds, ["fails", "meets", "meets", "fails", "meets", "meets"]))},
        ])
        self.assertEqual(record.seed_rows([h("a__d2", ["yes"] * 6)]), [])

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
        self.assertEqual([record.parse_call(x) for x in ("Resistant", "nonsusceptible", "susceptible", "intermediate", "susceptible-dose dependent",
                                                         "reduced susceptibility", "not defined", None)],
                         ["R", "R", "S", "I", "I", "I", None, None])
        self.assertEqual([record.consolidate(c) for c in ({"R", "S"}, {"R", "I"}, {"S"}, {"I"}, set())], ["conflict", "R", "S", "I only", "no call"])


if __name__ == "__main__":
    unittest.main()
