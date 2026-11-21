"""scripts/kover-comparison/record.py on SYNTHETIC inputs: cells worked by hand, the same bytes, stops and --check.

The inputs are built in a temporary directory by synthetic.make; none of them is data.
"""

import contextlib
import io
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile
import unittest

from tests.kover_comparison import synthetic
from tests.kover_comparison.common import RECORD, ROOT, SCRIPT, load_script, read_tsv

record = load_script()


def run(*args):
    """main() in this process; returns (rc, stderr)."""
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        rc = record.main(list(args))
    return rc, err.getvalue()


def files(directory):
    return {p.relative_to(directory).as_posix(): p.read_bytes()
            for p in sorted(pathlib.Path(directory).rglob("*")) if p.is_file()}


class Pieces(unittest.TestCase):
    def test_measures_by_hand(self):
        m = record.measures(2, 0, 1, 1)
        self.assertEqual((m["recall"], m["specificity"], m["balanced_accuracy"], m["error_rate"]), (1.0, 0.5, 0.75, 0.25))
        self.assertAlmostEqual(m["log_loss"], 0.0001 + 9.2102 * 0.25, places=4)
        # A one-class fold has no specificity to lose.
        self.assertEqual(record.measures(3, 1, 0, 0)["specificity"], 0.0)

    def test_quantile_is_type_7(self):
        self.assertEqual(record.quantile([0.0, 1.0, 2.0, 3.0, 4.0], 0.5), 2.0)
        self.assertAlmostEqual(record.quantile([0.0, 10.0], 0.025), 0.25)
        self.assertAlmostEqual(record.quantile([0.0, 10.0, 20.0], 0.975), 19.5)

    def test_bootstrap_of_equal_values_is_that_value(self):
        self.assertEqual(record.boot_mean([0.25, 0.25, 0.25], 1), (0.25, 0.25))
        lo, hi = record.boot_mean([0.0, 1.0, 2.0, 3.0], 20261115)
        self.assertTrue(0.0 <= lo < 1.5 < hi <= 3.0)
        self.assertEqual(record.boot_mean([0.0, 1.0, 2.0, 3.0], 20261115), (lo, hi))

    def test_counts_from_shares(self):
        self.assertEqual(record.count("0.333333", 3, "x"), 1)
        self.assertEqual(record.count("0.000000", 0, "x"), 0)
        with self.assertRaises(record.Stop):
            record.count("0.400000", 3, "x")

    def test_equal_within_1e_12(self):
        self.assertEqual(record.hel([1e-13, -1e-13, 0.1, -0.1, -0.2]), (1, 2, 2))

    def test_tied_values_and_f6(self):
        self.assertEqual(record.tied_values("0.5:67;0.916667:93"), [(0.5, 67), (0.916667, 93)])
        self.assertEqual(record.f6(-1e-9), "0.000000")


class Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="kover-comparison-")
        cls.kover, cls.review = synthetic.make(os.path.join(cls.tmp, "in"))
        cls.out = os.path.join(cls.tmp, "out")
        rc, err = run("--source", cls.kover, "--review", cls.review, "--out", cls.out)
        assert rc == 0, err
        cls.per = {r["dataset"]: r for r in read_tsv(os.path.join(cls.out, "per_dataset.tsv"))}
        cls.summary = {r["measure"]: r for r in read_tsv(os.path.join(cls.out, "summary.tsv"))}
        cls.totals = {r["key"]: r["value"] for r in read_tsv(os.path.join(cls.out, "totals.tsv"))}

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_files_are_the_record_s_files(self):
        self.assertEqual(sorted(files(self.out)), sorted(files(RECORD)))

    def test_per_dataset_cells_by_hand(self):
        for ds, (ba, rec, spec, err, rules, agree) in synthetic.BY_HAND.items():
            row = self.per[ds]
            self.assertEqual(row["balanced_accuracy"], record.f6(ba), ds)
            self.assertEqual(row["recall"], record.f6(rec), ds)
            self.assertEqual(row["specificity"], record.f6(spec), ds)
            self.assertEqual(row["error_rate"], record.f6(err), ds)
            self.assertEqual(row["rules"], record.f6(rules), ds)
            self.assertEqual(row["agreeing_share"], record.f6(agree), ds)
        self.assertEqual(self.per["genus_alpha__drugone"]["log_loss"], record.f6(synthetic.ONE_ERROR_IN_4 / 5))

    def test_t_interval_by_hand(self):
        # alpha's fold differences are -0.25, 0, 0, 0, 0: mean -0.05, sd sqrt(0.0125), half 2.7764451051977987 * 0.05.
        row = self.per["genus_alpha__drugone"]
        half = 2.7764451051977987 * 0.0125 ** 0.5 / 5 ** 0.5
        self.assertEqual(row["balanced_accuracy_t_low"], record.f6(-0.05 - half))
        self.assertEqual(row["balanced_accuracy_t_high"], record.f6(-0.05 + half))
        self.assertEqual(row["balanced_accuracy_t_low"], "-0.188822")
        gamma = self.per["genus_gamma__drugthree"]
        self.assertEqual((gamma["balanced_accuracy_t_low"], gamma["balanced_accuracy_t_high"]), ("0.000000", "0.000000"))

    def test_summary_by_hand(self):
        ba = self.summary["balanced_accuracy"]
        self.assertEqual(ba["value"], record.f6(-0.1 / 3))
        self.assertEqual((ba["datasets_higher"], ba["datasets_equal"], ba["datasets_lower"]), ("0", "1", "2"))
        self.assertEqual((ba["folds_higher"], ba["folds_equal"], ba["folds_lower"]), ("0", "13", "2"))
        self.assertEqual(ba["seed"], "20261115")
        self.assertLessEqual(float(ba["low"]), float(ba["value"]))
        self.assertLessEqual(float(ba["value"]), float(ba["high"]))
        self.assertGreaterEqual(float(ba["low"]), -0.05)
        self.assertLessEqual(float(ba["high"]), 0.0)
        self.assertEqual(self.summary["rules"]["value"], "0.000000")
        self.assertEqual(self.summary["agreeing"]["value"], record.f6((0.95 + 0.95 + 1.0) / 3))
        self.assertEqual([r for r in self.summary], record.MEASURES)

    def test_tie_readings_by_hand(self):
        # alpha fold 0: tied 0.75 and 1.0 against GRM's 1.0, so equal -0.125, lowest -0.25, highest 0, over 5 folds.
        # beta has no tie rows; gamma's tied rules give one result.
        self.assertEqual(self.summary["ties_equal"]["value"], record.f6((-0.025 - 0.05 + 0.0) / 3))
        self.assertEqual(self.summary["ties_lowest"]["value"], record.f6((-0.05 - 0.05 + 0.0) / 3))
        self.assertEqual(self.summary["ties_highest"]["value"], record.f6((0.0 - 0.05 + 0.0) / 3))
        for key in ("ties_equal", "ties_lowest", "ties_highest", "error_rate"):
            self.assertEqual(self.summary[key]["reading"], "exploratory")

    def test_totals_by_hand(self):
        want = {"datasets": "3", "species": "3", "genomes": "60", "folds": "15", "kover_one_rule_folds": "10",
                "ties_more_than_one_rule": "6", "ties_different_test_results": "1", "ties_first_index_equals_run": "10",
                "datasets_identical_calls": "1", "datasets_t_interval_excludes_0": "0", "steps": "3", "steps_fits": "2",
                "steps_counts_equal": "3", "steps_pattern_equal": "1", "step1_pattern_differs": "1",
                "step1_differs_equals_a_listed_tied": "1", "kover_empty_folds": "0", "grm_rules_total": "20",
                "kover_rules_total": "20", "kover_conjunction_folds": "10"}
        for key, value in want.items():
            self.assertEqual(self.totals[key], value, key)

    def test_order_and_settings(self):
        order = read_tsv(os.path.join(self.out, "order.tsv"))
        self.assertEqual([r["event"] for r in order],
                         ["plan_hashed"] * 5 + ["first_kover_fit", "review_plan_hashed", "review_tie_summary"])
        self.assertEqual(order[0]["utc"], synthetic.PLAN_TIME)
        settings = [(r["key"], r["value"]) for r in read_tsv(os.path.join(self.out, "settings.tsv"))]
        self.assertEqual(settings, record.SETTINGS)

    def test_same_bytes_twice_and_at_other_hash_seeds(self):
        again = os.path.join(self.tmp, "again")
        self.assertEqual(run("--source", self.kover, "--review", self.review, "--out", again)[0], 0)
        self.assertEqual(files(again), files(self.out))
        for seed in ("0", "4242"):
            other = os.path.join(self.tmp, "seed" + seed)
            env = {"PYTHONHASHSEED": seed, "PATH": os.environ.get("PATH", "/usr/bin:/bin")}
            subprocess.run([sys.executable, "-I", str(SCRIPT), "--source", self.kover, "--review", self.review,
                            "--out", other], check=True, env=env, cwd=ROOT, capture_output=True)
            self.assertEqual(files(other), files(self.out), seed)

    def test_check_passes_on_its_own_output(self):
        self.assertEqual(run("--check", self.out), (0, ""))
        self.assertEqual(run("--source", self.kover, "--review", self.review, "--check", self.out), (0, ""))

    def test_check_fails_on_a_changed_cell_a_dropped_file_and_an_extra_file(self):
        for name, change in (("folds.tsv", "cell"), ("summary.tsv", "cell"), ("ties.tsv", "drop"),
                             ("README.md", "drop"), ("genomes.tsv", "extra")):
            copy = os.path.join(self.tmp, "copy-%s-%s" % (name, change))
            shutil.copytree(self.out, copy)
            path = os.path.join(copy, name)
            if change == "cell":
                text = pathlib.Path(path).read_text(encoding="utf-8")
                lines = text.split("\n")
                cells = lines[1].split("\t")
                cells[3] = cells[3] + "1"
                lines[1] = "\t".join(cells)
                pathlib.Path(path).write_text("\n".join(lines), encoding="utf-8")
            elif change == "drop":
                os.remove(path)
            else:
                pathlib.Path(path).write_text("genome\tcall\nSYNTHETIC\t1\n", encoding="utf-8")
            rc, err = run("--check", copy)
            self.assertEqual(rc, 1, (name, change))
            self.assertIn("record.py:", err)


class Stops(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="kover-comparison-stop-")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def stops(self, kover, review, needle):
        rc, err = run("--source", kover, "--review", review, "--out", os.path.join(self.tmp, "out"))
        self.assertEqual(rc, 1)
        self.assertIn("stopped", err)
        self.assertIn(needle, err)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "out")))

    def edit(self, path, old, new):
        text = pathlib.Path(path).read_text(encoding="utf-8")
        self.assertIn(old, text)
        pathlib.Path(path).write_text(text.replace(old, new, 1), encoding="utf-8")

    def test_datasets_changed_after_the_plan_hashed_them(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"), datasets_after_hash="order\tdataset\n")
        self.stops(k, r, "is not the file the plan hashed")

    def test_a_plan_without_a_setting(self):
        plan = synthetic.plan_text().replace("--max-rules 6", "--max-rules 5")
        k, r = synthetic.make(os.path.join(self.tmp, "in"), plan=plan)
        self.stops(k, r, "PLAN.md does not give")

    def test_a_settings_row_that_differs(self):
        plan = synthetic.plan_text().replace("| Selection criterion | `mean balanced error` |", "| Selection criterion | `mean error` |")
        k, r = synthetic.make(os.path.join(self.tmp, "in"), plan=plan)
        self.stops(k, r, "the row for Selection criterion differs")

    def test_a_first_fit_before_the_plan(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        self.edit(os.path.join(k, "results", "run_all.log"), synthetic.FIT_TIME, "2026-11-15T11:00:00Z")
        self.stops(k, r, "before the first Kover fit")

    def test_a_fold_measure_that_does_not_follow_from_its_counts(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        self.edit(os.path.join(k, "results", "per_fold.tsv"), "\t1.000000\t0.750000\t1.000000\t1.000000\t1.000000\t0.500000\t",
                  "\t1.000000\t0.700000\t1.000000\t1.000000\t1.000000\t0.500000\t")
        self.stops(k, r, "does not follow from the counts")

    def test_a_share_that_is_not_a_count(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        self.edit(os.path.join(k, "results", "per_fold.tsv"), "genus_beta__drugtwo\t2\t4\t2\t2\t0.750000",
                  "genus_beta__drugtwo\t2\t4\t2\t2\t0.700000")
        self.stops(k, r, "is not a count")

    def test_a_first_tied_rule_that_does_not_give_kover_s_calls(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        path = os.path.join(k, "results", "review", "fold_ties_genus_alpha__drugone.json")
        self.edit(path, '"first_index_reproduces_kover_calls": true', '"first_index_reproduces_kover_calls": false')
        self.stops(k, r, "does not give Kover's calls")

    def test_tied_counts_that_do_not_add_up(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        path = os.path.join(k, "results", "review", "fold_ties_genus_gamma__drugthree.json")
        self.edit(path, "[[0.75, 2]]", "[[0.75, 3]]")
        self.stops(k, r, "tied rule counts do not add up")

    def test_a_review_plan_changed_after_its_hash(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        self.edit(os.path.join(r, "PLAN.md"), "SYNTHETIC", "SYNTHETIC, changed")
        self.stops(k, r, "the review's PLAN.sha256 does not hold")

    def test_a_missing_input(self):
        k, r = synthetic.make(os.path.join(self.tmp, "in"))
        os.remove(os.path.join(k, "results", "final_tied.tsv"))
        self.stops(k, r, "final_tied.tsv")

    def test_a_record_whose_settings_differ(self):
        out = os.path.join(self.tmp, "rec")
        shutil.copytree(RECORD, out)
        self.edit(os.path.join(out, "settings.tsv"), "--max-rules 6", "--max-rules 5")
        rc, err = run("--check", out)
        self.assertEqual(rc, 1)
        self.assertIn("settings.tsv differs from the script's settings", err)

    def test_arguments(self):
        with contextlib.redirect_stderr(io.StringIO()):
            for args in (["--source", "x", "--check", "y"], ["--out", "y"], []):
                with self.assertRaises(SystemExit):
                    record.main(args)


if __name__ == "__main__":
    unittest.main()
