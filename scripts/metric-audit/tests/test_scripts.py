"""nullimpact.py, transitions.py, and summarize.py on a synthetic revision-7
run: a null fold stays null in each of them instead of turning into 0 or
0.5, and an earlier 0.5 fold that is now null shows as a move to no call.

usage: python3 -I tests/test_scripts.py
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixtures as fx  # noqa: E402

AUDIT = fx.LANE
LAUNCH = "2026-01-01T00:00:00Z"


def script(name, *args):
    result = subprocess.run([sys.executable, "-I", os.path.join(AUDIT, name), *args], capture_output=True, text=True)
    if result.returncode != 0:
        raise AssertionError(result.stderr[-2000:])
    return result.stdout


def held(scheme, fold, resistant):
    return sum(1 for i in range(fx.N) if fx.FOLDS[scheme][i] == fold and fx.TRUTH[i] == resistant)


class Scripts(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.report = fx.make_report()
        self.paths = fx.audit_root(self.tmp.name, self.report, fx.manifest(self.report, "rf"),
                                   fx.make_report(revision=6, forest_variant=1),
                                   [fx.recompute_line(), fx.done_line("rf", fx.row(self.report, "random", "rf")["metrics"]["balanced_accuracy"])])
        self.run_args = ["--revision", "7", "--root", self.paths["root"]]

    def tearDown(self):
        self.tmp.cleanup()

    def out(self, name):
        return os.path.join(self.tmp.name, name)

    def test_transitions_show_each_half_fold_moving_to_no_call(self):
        script("transitions.py", self.out("t.json"), *self.run_args, "--backup", self.paths["backup"])
        t = json.load(open(self.out("t.json")))
        self.assertEqual(t["null_fold_transitions"], {"lineage/gbt R->no call": 1, "random/lasso R->no call": 1})
        flipped = t["null_fold_samples_flipped"]
        self.assertEqual((flipped["R_now_uncalled"], flipped["S_now_uncalled"]),
                         (held("lineage", 2, True) + held("random", 4, True), held("lineage", 2, False) + held("random", 4, False)))
        self.assertEqual({(r["scheme"], r["learner"], r["old_value"], r["new_value"]) for r in t["null_rows"]},
                         {("lineage", "gbt", 0.5, None), ("random", "lasso", 0.5, None)})
        changed = {(p["scheme"], p["learner"]) for p in t["pooled_changes"]}
        self.assertTrue({("lineage", "gbt"), ("random", "lasso")} <= changed, changed)

    def test_nullimpact_keeps_a_row_without_calls_null(self):
        gbt = fx.row(self.report, "random", "gbt")
        gbt["scores"] = [None] * fx.N
        gbt["fold_complexity"] = [0.0] * 5
        fx.rescore(self.report)
        fx.write_json(self.paths["report"], self.report)
        script("nullimpact.py", self.out("n.json"), *self.run_args)
        n = json.load(open(self.out("n.json")))
        rows = {(r["scheme"], r["learner"]): r for r in n["rows"]}
        self.assertEqual(set(rows), {("lineage", "gbt"), ("random", "lasso"), ("random", "gbt")})
        self.assertEqual(n["complexity_mismatch"], [])
        self.assertEqual(n["rows_without_any_call"], 1)
        empty = rows[("random", "gbt")]
        self.assertEqual((empty["n"], empty["n_without_null"], empty["vme_without_null"], empty["stored"]), (fx.N, None, None, "without_null"))
        # Scored 0.5, every isolate is a resistant call: no very major error, every susceptible a major error.
        self.assertEqual((empty["vme"], empty["me"]), (0.0, 1.0))
        lasso = rows[("random", "lasso")]
        stored = fx.row(self.report, "random", "lasso")["metrics"]
        self.assertEqual((lasso["n_without_null"], lasso["vme_without_null"]), (stored["n"], stored["vme"]))

    def test_summarize_keeps_null_fold_values_null(self):
        script("audit.py", self.out("run.json"), *self.run_args, "--backup", self.paths["backup"], "--log", self.paths["log"], "--launch", LAUNCH)
        run = json.load(open(self.out("run.json")))
        self.assertEqual(run["summary"]["issue_count"], 0, run["records"][0]["issues"])
        # A second record whose null fold carries a number and no complexity: values of both kinds must sort.
        extra = json.loads(json.dumps(run["records"][0]))
        extra["id"] = "fixture_species__other"
        extra["null_folds"][0].update(fold_ba=0.5, complexity=None)
        run["records"].append(extra)
        json.dump(run, open(self.out("run.json"), "w"))
        script("summarize.py", self.out("run.json"), self.out("tables.json"))
        tables = json.load(open(self.out("tables.json")))
        cov = tables["coverage"]["7"]
        self.assertEqual((cov["fold_ba_null"], cov["rows_with_no_call"], cov["rows_without_calls"]), (3, 4, 0))
        self.assertEqual(cov["no_call_isolates"], {"lineage/gbt": 2 * (held("lineage", 2, True) + held("lineage", 2, False)),
                                                   "random/lasso": 2 * (held("random", 4, True) + held("random", 4, False))})
        nulls = tables["null_models"]["7"]
        self.assertEqual((nulls["fold_ba_values"], nulls["complexities"]), ([None, 0.5], [None, 0.0]))
        self.assertEqual(cov["best_learner"], {"rf": 2})
        self.assertEqual(tables["vs_backup"]["learners_with_changed_scores"]["random/rf"], 2)


if __name__ == "__main__":
    unittest.main(verbosity=2)
