"""scripts/dense-network/record.py on SYNTHETIC inputs: cells worked by hand, the same bytes, stops and --check.

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

from tests.dense_network import synthetic
from tests.dense_network.common import RECORD, ROOT, SCRIPT, load_script, read_tsv

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


def tables(directory):
    d = pathlib.Path(directory)
    summary = {(r["split"], r["comparator"], r["measure"]): r for r in read_tsv(d / "summary.tsv")}
    totals = {r["key"]: r["value"] for r in read_tsv(d / "totals.tsv")}
    return summary, read_tsv(d / "rule.tsv"), totals


class Pieces(unittest.TestCase):
    def test_percentile_is_type_7(self):
        self.assertEqual(record.pctl([0.0, 1.0, 2.0, 3.0, 4.0], 0.5), 2.0)
        self.assertAlmostEqual(record.pctl([0.0, 10.0], 0.025), 0.25)
        self.assertAlmostEqual(record.pctl([0.0, 10.0, 20.0], 0.975), 19.5)

    def test_interval_of_two_species_by_hand(self):
        # Two species of one value each: a draw of two gives a median of 1, 2 or 3, with 1 and 3 a quarter of
        # the draws each, so the 2.5% and 97.5% percentiles of 2,000 draws are 1 and 3.
        iv = record.interval({"b": [3.0], "a": [1.0]})
        self.assertEqual((iv["median"], iv["lo"], iv["hi"]), (2.0, 1.0, 3.0))
        self.assertEqual((iv["n"], iv["species"], iv["positive"], iv["negative"]), (2, 2, 2, 0))

    def test_interval_takes_every_dataset_of_a_drawn_species(self):
        # One species of three datasets: every draw takes all three, so the interval is the median itself.
        iv = record.interval({"a": [0.1, -0.2, 0.3]})
        self.assertEqual((iv["median"], iv["lo"], iv["hi"]), (0.1, 0.1, 0.1))
        self.assertEqual((iv["positive"], iv["negative"]), (2, 1))

    def test_interval_is_the_same_whatever_order_the_species_come_in(self):
        one = record.interval({"a": [0.1, 0.2], "b": [-0.3], "c": [0.05, 0.0, 0.4]})
        two = record.interval({"c": [0.05, 0.0, 0.4], "b": [-0.3], "a": [0.1, 0.2]})
        self.assertEqual(one, two)

    def test_cells_and_f6(self):
        self.assertEqual(record.cell(0.1), "0.1")
        self.assertEqual(record.cell(None), "")
        self.assertEqual(float(record.cell(0.36276337958843685)), 0.36276337958843685)
        self.assertEqual(record.f6(-1e-9), "0.000000")
        self.assertEqual(record.f6(0.0024216), "0.002422")


class Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="dense-network-")
        cls.run_dir, cls.review = synthetic.make(os.path.join(cls.tmp, "in"))
        cls.out = os.path.join(cls.tmp, "out")
        rc, err = run("--source", cls.run_dir, "--review", cls.review, "--out", cls.out)
        assert rc == 0, err
        cls.summary, cls.rule, cls.totals = tables(cls.out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_files_are_the_record_s_files(self):
        self.assertEqual(sorted(files(self.out)), sorted(files(RECORD)))

    def test_summary_by_hand(self):
        h = synthetic.BY_HAND
        for split in record.SPLITS:
            for who, k in (("best5", "balanced_accuracy"), ("best5", "recall"), ("best5", "log_loss"),
                           ("site_best", "balanced_accuracy"), ("baseline", "balanced_accuracy")):
                row = self.summary[(split, who, k)]
                want = record.f6(h["%s_%s" % (who, k)])
                self.assertEqual((row["difference_median"], row["low"], row["high"]), (want, want, want), (split, who, k))
                self.assertEqual((row["datasets"], row["species"]), ("4", "2"))
            self.assertEqual(self.summary[(split, "best5", "balanced_accuracy")]["candidate_median"],
                             record.f6(h["network_median_balanced_accuracy"]))
            self.assertEqual(self.summary[(split, "best5", "balanced_accuracy")]["negative"], "4")
            self.assertEqual(self.summary[(split, "catalog", "balanced_accuracy")]["datasets"], "3")
        # The baseline has balanced accuracy only.
        self.assertNotIn(("lineage", "baseline", "recall"), self.summary)

    def test_rule_by_hand(self):
        h = synthetic.BY_HAND
        values = [(r["rule"], r["value"], r["holds"]) for r in self.rule]
        self.assertEqual(values, [
            ("1", record.f6(h["best5_balanced_accuracy"]), "no"), ("1", record.f6(h["best5_balanced_accuracy"]), "no"),
            ("2", record.f6(h["margin_candidate"]), "no"), ("3", record.f6(h["best5_recall"]), "no"),
            ("4", "1.000000", "yes"), ("4", "1000", "yes"), ("all", "no gain established", "no")])
        self.assertEqual(self.rule[2]["needs"], "at least %s, the comparator's" % record.f6(h["margin_best5"]))

    def test_totals_by_hand(self):
        t = self.totals
        self.assertEqual((t["datasets"], t["species"], t["datasets_largest_gap"], t["datasets_drawn"]), ("4", "2", "2", "2"))
        self.assertEqual((t["genomes_smallest"], t["genomes_largest"]), ("40", "80"))
        self.assertEqual((t["gap_largest_group_low"], t["gap_largest_group_high"]), ("0.250000", "0.300000"))
        self.assertEqual(t["catalog_datasets"], str(synthetic.BY_HAND["catalog_datasets"]))
        self.assertEqual((t["site_best_is_best5_lineage"], t["site_best_is_best5_random"]), ("0", "0"))
        for who, gap in synthetic.BY_HAND["below_best5"].items():
            self.assertEqual(t["below_best5_lineage_" + who], record.f6(gap), who)

    def test_a_network_above_the_best_of_five_is_adopted(self):
        run_dir, review = synthetic.make(os.path.join(self.tmp, "lift"), lift=True)
        out = os.path.join(self.tmp, "lift-out")
        self.assertEqual(run("--source", run_dir, "--review", review, "--out", out), (0, ""))
        rule = read_tsv(os.path.join(out, "rule.tsv"))
        self.assertEqual([r["holds"] for r in rule], ["yes"] * 7)
        self.assertEqual(rule[-1]["value"], "adopt")

    def test_order_and_settings(self):
        order = read_tsv(os.path.join(self.out, "order.tsv"))
        self.assertEqual([r["event"] for r in order], record.ORDER_EVENTS)
        self.assertEqual([r["utc"][11:] for r in order],
                         [synthetic.PLAN_TIME[11:], synthetic.FITTING, synthetic.FIRST, synthetic.LOOK, synthetic.LAST])
        self.assertEqual(order[0]["sha256"], synthetic.sha(synthetic.plan_text()))
        settings = [(r["key"], r["value"]) for r in read_tsv(os.path.join(self.out, "settings.tsv"))]
        self.assertEqual(settings, record.SETTINGS)

    def test_same_bytes_twice_and_at_other_hash_seeds(self):
        again = os.path.join(self.tmp, "again")
        self.assertEqual(run("--source", self.run_dir, "--review", self.review, "--out", again)[0], 0)
        self.assertEqual(files(again), files(self.out))
        for seed in ("0", "4242"):
            other = os.path.join(self.tmp, "seed" + seed)
            env = {"PYTHONHASHSEED": seed, "PATH": os.environ.get("PATH", "/usr/bin:/bin")}
            subprocess.run([sys.executable, "-I", str(SCRIPT), "--source", self.run_dir, "--review", self.review,
                            "--out", other], check=True, env=env, cwd=ROOT, capture_output=True)
            self.assertEqual(files(other), files(self.out), seed)

    def test_check_passes_on_its_own_output(self):
        self.assertEqual(run("--check", self.out), (0, ""))
        self.assertEqual(run("--source", self.run_dir, "--review", self.review, "--check", self.out), (0, ""))

    def test_check_fails_on_a_changed_cell_a_dropped_file_and_an_extra_file(self):
        for name, change in (("measures.tsv", "cell"), ("summary.tsv", "cell"), ("totals.tsv", "drop"),
                             ("README.md", "drop"), ("genomes.tsv", "extra")):
            copy = os.path.join(self.tmp, "copy-%s-%s" % (name, change))
            shutil.copytree(self.out, copy)
            path = os.path.join(copy, name)
            if change == "cell":
                text = pathlib.Path(path).read_text(encoding="utf-8")
                lines = text.split("\n")
                cells = lines[1].split("\t")
                cells[4] = cells[4] + "1"
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
        self.tmp = tempfile.mkdtemp(prefix="dense-network-stop-")

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def stops(self, run_dir, review, needle):
        rc, err = run("--source", run_dir, "--review", review, "--out", os.path.join(self.tmp, "out"))
        self.assertEqual(rc, 1)
        self.assertIn("stopped", err)
        self.assertIn(needle, err)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "out")))

    def make(self, **kw):
        return synthetic.make(os.path.join(self.tmp, "in"), **kw)

    def test_a_plan_whose_sha256_the_report_does_not_give(self):
        self.stops(*self.make(stated="0" * 64), "PLAN.md is not the plan whose sha256 the report gives")

    def test_a_plan_without_a_setting(self):
        plan = synthetic.plan_text().replace("at most 300 epochs", "at most 200 epochs")
        self.stops(*self.make(plan=plan), "PLAN.md does not give stopping")

    def test_a_figure_that_parts_from_the_experiment_s(self):
        def change(s, stage1, gates):
            s["paired"]["lineage"]["rf"]["recall"]["median"] += 1e-9
        self.stops(*self.make(change=change), "lineage rf recall median")

    def test_an_outcome_that_parts_from_the_experiment_s(self):
        def change(s, stage1, gates):
            s["outcome"] = "adopt"
        self.stops(*self.make(change=change), "outcome")

    def test_a_dataset_left_out(self):
        def change(s, stage1, gates):
            s["excluded"] = ["genus_beta__drugfour"]
        self.stops(*self.make(change=change), "datasets left out")

    def test_another_seed_for_the_draw(self):
        def change(s, stage1, gates):
            stage1["seed"] = 1
        self.stops(*self.make(change=change), "seed or pool differs from the plan")

    def test_gates_that_disagree(self):
        def change(s, stage1, gates):
            gates[1]["gate_b"] = False
        self.stops(*self.make(change=change), "disagree on gate B")

    def test_rows_out_of_the_plan_s_order(self):
        def change(s, stage1, gates):
            s["rows"].reverse()
        self.stops(*self.make(change=change), "rows are not the plan's datasets in order")

    def test_a_measure_out_of_range(self):
        def change(s, stage1, gates):
            s["rows"][0]["lineage"]["gbt"]["recall"] = 1.5
        self.stops(*self.make(change=change), "out of range")

    def test_a_best_of_five_that_is_not_the_best(self):
        def change(s, stage1, gates):
            row = s["rows"][2]["random"]
            row["best5"] = dict(row["cart"], learner="cart")
        self.stops(*self.make(change=change), "best5 is not the highest balanced accuracy of the five")

    def test_a_first_look_before_the_plan(self):
        review = ("The experiment ran its analysis on the outcomes of the first 2 datasets at 11:00:00Z. The analysis "
                  "script was last changed at 12:30:00Z. `ba2_candidate.py` last change 12:02:00Z;\n")
        self.stops(*self.make(review=review), "the times are not UTC times in order")

    def test_a_review_that_does_not_give_the_times(self):
        self.stops(*self.make(review="SYNTHETIC review without times\n"), "review REPORT.md: 0 matches")

    def test_a_missing_input(self):
        run_dir, review = self.make()
        os.remove(os.path.join(run_dir, "results", "gates.json"))
        self.stops(run_dir, review, "gates.json")

    def test_a_record_whose_settings_differ(self):
        out = os.path.join(self.tmp, "rec")
        shutil.copytree(RECORD, out)
        path = os.path.join(out, "settings.tsv")
        text = pathlib.Path(path).read_text(encoding="utf-8")
        self.assertIn("at most 300 epochs", text)
        pathlib.Path(path).write_text(text.replace("at most 300 epochs", "at most 200 epochs"), encoding="utf-8")
        rc, err = run("--check", out)
        self.assertEqual(rc, 1)
        self.assertIn("settings.tsv differs from the script's settings", err)

    def test_a_record_without_a_row_file(self):
        out = os.path.join(self.tmp, "rec")
        shutil.copytree(RECORD, out)
        os.remove(os.path.join(out, "measures.tsv"))
        rc, err = run("--check", out)
        self.assertEqual(rc, 1)
        self.assertIn("measures.tsv is missing", err)

    def test_arguments(self):
        with contextlib.redirect_stderr(io.StringIO()):
            for args in (["--source", "x", "--check", "y"], ["--out", "y"], ["--review", "x", "--out", "y"], []):
                with self.assertRaises(SystemExit):
                    record.main(args)


if __name__ == "__main__":
    unittest.main()
