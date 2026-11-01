"""scripts/hyperopt/report.py on small hand-built histories.

usage: python3 -I -B tests/hyperopt/test_report.py
"""
import copy
import hashlib
import importlib.util
import io
import json
import shutil
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


report = _load("grm_hyperopt_report", ROOT / "scripts/hyperopt/report.py")


def protocol(plan="plan-A", constraints=(), names=("cart", "lasso")):
    """Task B's protocol shape, hashed as Task B hashes it."""
    p = {"schema_version": 1, "fold_plan_hash": plan,
         "objective": {"name": "auc", "direction": "maximize",
                       "constraints": [dict(c) for c in constraints]},
         "seeds": {"sampler": 3, "fit": 1}, "trials_per_family": 3,
         "families": {f: {"default": {"depth": 0}, "space": {"depth": {
             "type": "int", "low": 0, "high": 9, "log": False}}}
             for f in names}}
    p["protocol_hash"] = rehash(p)
    return p


def rehash(p):
    body = {k: v for k, v in p.items() if k != "protocol_hash"}
    return hashlib.sha256(json.dumps(
        body, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False).encode("utf-8")).hexdigest()


def trial(p, tid, scores, family="cart", status="ok", reason=None,
          config=None, complexity=1.0, cost=(None, 1.5)):
    """cost is (fit_count, seconds); fit_count None means one fit per fold."""
    fits, seconds = cost
    t = {"protocol_hash": p["protocol_hash"], "trial_id": tid,
         "family": family, "config": config or {"depth": tid},
         "complexity": complexity, "is_default": tid == 0, "seed": 1,
         "fold_plan_hash": p["fold_plan_hash"], "status": status,
         "inner_scores": [{"fold": f, "metrics": {"auc": s}}
                          for f, s in enumerate(scores)],
         "fit_count": len(scores) if fits is None else fits,
         "seconds": seconds}
    if reason is not None:
        t["reason"] = reason
    return t


def history(p, depth=3):
    return [
        trial(p, 0, [0.70, 0.72]),
        trial(p, 1, [0.80, 0.82], config={"depth": depth}, complexity=3),
        # A failed record as Task B's search writes one.
        trial(p, 2, [], status="failed",
              reason="tree | split crashed\nat node 4", complexity=None,
              cost=(0, None)),
        trial(p, 0, [0.60, 0.62], family="lasso"),
        trial(p, 1, [0.90, 0.90], family="lasso", status="unconverged",
              reason="path stopped at lambda 7"),
    ]


def run(label, plan="plan-A", depth=3, constraints=()):
    p = protocol(plan, constraints)
    return {"label": label, "protocol": p, "trials": history(p, depth)}


class Report(unittest.TestCase):
    def setUp(self):
        self.report = report.build_report(
            [run("outer0"), run("outer1", "plan-B", depth=5)])
        self.md = report.render_markdown(self.report)

    def test_every_trial_is_listed_with_status_and_reason(self):
        for r, h in zip(self.report["runs"],
                        [history(protocol()), history(protocol("plan-B"), 5)]):
            self.assertEqual(r["n_trials"], len(h))
            got = sorted((t["family"], t["trial_id"], t["status"], t["reason"])
                         for t in r["trials"])
            want = sorted((t["family"], t["trial_id"], t["status"],
                           t.get("reason")) for t in h)
            self.assertEqual(got, want)
        self.assertIn("| cart | 2 | failed | tree \\| split crashed at node 4 "
                      "| status is failed |", self.md)
        self.assertIn("| lasso | 1 | unconverged | path stopped at lambda 7 "
                      "| status is unconverged |", self.md)
        self.assertEqual(self.md.count("### All trials (5)"), 2)
        all_rows = [line for line in self.md.splitlines()
                    if line.startswith("| cart |") or line.startswith("| lasso |")]
        # Per run: 5 trial rows, 2 excluded rows, 2 selection rows, 2 budgets
        # and 1 row of unused fold scores.
        self.assertEqual(len(all_rows), 2 * (5 + 2 + 2 + 2 + 1))

    def test_unused_fold_scores_are_kept(self):
        rows = {(t["family"], t["trial_id"]): t
                for t in self.report["runs"][0]["trials"]}
        self.assertEqual(rows["lasso", 1]["fold_scores"],
                         [{"fold": 0, "score": 0.9}, {"fold": 1, "score": 0.9}])
        self.assertFalse(rows["lasso", 1]["scores_used"])
        self.assertTrue(rows["lasso", 0]["scores_used"])
        self.assertIn("### Fold scores not used (1)", self.md)
        self.assertIn("| lasso | 1 | unconverged | 0: 0.9, 1: 0.9 |", self.md)

    def test_history_digest_is_reported(self):
        r = self.report["runs"][0]
        text = "".join(json.dumps(t, sort_keys=True, separators=(",", ":"))
                       + "\n" for t in history(protocol()))
        self.assertEqual(r["history_sha256"],
                         hashlib.sha256(text.encode()).hexdigest())
        self.assertIn(f"history sha256 {r['history_sha256']}", self.md)

    def test_selection_candidates_and_budgets(self):
        r = self.report["runs"][0]
        cart = r["families"]["cart"]
        self.assertEqual(cart["selected"]["trial_id"], 1)
        self.assertEqual([c["trial_id"] for c in cart["candidate_set"]], [1])
        self.assertEqual(cart["excluded"][0]["reason"],
                         "tree | split crashed\nat node 4")
        lasso = r["families"]["lasso"]
        self.assertEqual(lasso["selected"]["trial_id"], 0)
        self.assertEqual(r["budgets"]["families"]["cart"]["failed"], 1)
        self.assertEqual({f["code"] for f in r["budgets"]["flags"]},
                         {"under_budget", "unequal_used", "unequal_ok"})
        self.assertIn("- under_budget: lasso used 2 of 3 planned trials",
                      self.md)
        selected = [t for t in r["trials"] if t["selected"]]
        self.assertEqual({(t["family"], t["trial_id"]) for t in selected},
                         {("cart", 1), ("lasso", 0)})

    def test_stability_across_runs(self):
        cart = self.report["stability"]["families"]["cart"]
        self.assertEqual(cart["selected_in"], 2)
        self.assertEqual(cart["params"]["depth"]["value_counts"],
                         [{"value": 3, "count": 1, "fraction": 0.5},
                          {"value": 5, "count": 1, "fraction": 0.5}])
        self.assertEqual(self.report["stability"]["families"]["lasso"]
                         ["default_selected"], 2)
        self.assertIn("## Stability across runs (2)", self.md)
        self.assertIn("| depth | 2 | 2 | 0.5 | 3 | 4 | 5 | 0 |", self.md)

    def test_se_is_labelled_descriptive_in_plain_text(self):
        self.assertIn("descriptive", self.report["se_label"])
        se_lines = [line for line in self.md.splitlines()
                    if line.startswith("SE note:")]
        self.assertEqual(len(se_lines), 1)
        self.assertIn("SE is descriptive", se_lines[0])
        self.assertIn("Inner folds share training data", se_lines[0])
        self.assertIn("not an independent-sample inference", se_lines[0])
        self.assertIn("SE (descriptive)", self.md)

    def test_constraints_show_as_exclusions(self):
        # auc >= 0.75 rules out cart trial 0 (mean 0.71) and lasso's only ok
        # trial (mean 0.61).
        r = run("outer0", constraints=[{"metric": "auc", "op": ">=",
                                        "bound": 0.75}])
        rep = report.build_report([r])
        md = report.render_markdown(rep)
        fams = rep["runs"][0]["families"]
        self.assertEqual(fams["cart"]["selected"]["trial_id"], 1)
        self.assertEqual(fams["lasso"]["status"], "no_feasible_trials")
        self.assertIn("Constraints: auc >= 0.75.", md)
        self.assertIn("| lasso | none (no_feasible_trials) |", md)
        # In the excluded table and, with complexity and default, in the
        # all-trials table.
        self.assertIn("| cart | 0 | ok |  | constraint auc >= 0.75 not met "
                      "(auc is 0.71) |\n", md)
        self.assertIn("| cart | 0 | ok |  | constraint auc >= 0.75 not met "
                      "(auc is 0.71) | 1 | yes |", md)
        self.assertIn("### Excluded trials (4)", md)
        self.assertIn("unequal_feasible", {f["code"] for f in
                                           rep["runs"][0]["budgets"]["flags"]})

    def test_report_is_strict_json(self):
        text = json.dumps(self.report, allow_nan=False)
        self.assertEqual(json.loads(text)["schema"], "grm-hyperopt-report/1")

    def test_refusals_name_the_run(self):
        bad = run("outer1")
        bad["trials"][1]["protocol_hash"] = "other"
        with self.assertRaises(report.HistoryError) as cm:
            report.build_report([run("outer0"), bad])
        self.assertTrue(all(p.endswith("[run outer1]")
                            for p in cm.exception.problems))
        self.assertIn("mixed_protocol_hash", cm.exception.codes)


class Cli(unittest.TestCase):
    def setUp(self):
        self.dir = ROOT / ".tmp" / "test_report_cli"
        shutil.rmtree(self.dir, ignore_errors=True)
        self.dir.mkdir(parents=True)

    def tearDown(self):
        shutil.rmtree(self.dir, ignore_errors=True)

    def write(self, name, text):
        (self.dir / name).write_text(text)
        return str(self.dir / name)

    def main(self, *runs, label="outer", out_md=None):
        """Run the CLI on (protocol text, history text) pairs; returns the
        exit code and stderr."""
        args = []
        for i, (proto, hist) in enumerate(runs):
            args += ["--run", f"{label}{i}", self.write(f"p{i}.json", proto),
                     self.write(f"h{i}.jsonl", hist)]
        err = io.StringIO()
        with redirect_stderr(err):
            rc = report.main(args + ["--out-json", str(self.dir / "r.json"),
                                     "--out-md", out_md or str(self.dir / "r.md")])
        return rc, err.getvalue()

    def refused(self, *runs, **kw):
        rc, err = self.main(*runs, **kw)
        self.assertEqual(rc, 2)
        self.assertNotIn("Traceback", err)
        self.assertFalse((self.dir / "r.json").exists())
        self.assertFalse((self.dir / "r.md").exists())
        self.assertEqual([p.name for p in self.dir.iterdir()
                          if p.name.endswith(".tmp")], [])
        return err

    def test_writes_json_and_markdown(self):
        p = protocol()
        rc, err = self.main((json.dumps(p), lines(history(p))))
        self.assertEqual((rc, err), (0, ""))
        out = json.loads((self.dir / "r.json").read_text())
        self.assertEqual(len(out["runs"][0]["trials"]), 5)
        self.assertIn("SE note:", (self.dir / "r.md").read_text())
        self.assertEqual(sorted(p.name for p in self.dir.iterdir()),
                         ["h0.jsonl", "p0.json", "r.json", "r.md"])

    def test_refuses_and_writes_nothing(self):
        p = protocol()
        h = history(p)
        h.append(h[0])
        self.assertIn("duplicate_trial_id",
                      self.refused((json.dumps(p), lines(h))))

    def test_refuses_an_edited_protocol_with_a_stale_hash(self):
        p = protocol()
        tampered = copy.deepcopy(p)
        tampered["objective"]["direction"] = "minimize"
        err = self.refused((json.dumps(tampered), lines(history(p))))
        self.assertIn("protocol_hash_mismatch", err)

    def test_refuses_a_flat_protocol(self):
        p = protocol()
        flat = {"objective": "auc", "direction": "maximize",
                "protocol_hash": p["protocol_hash"],
                "fold_plan_hash": p["fold_plan_hash"],
                "budgets": {"cart": 3, "lasso": 3}}
        self.assertIn("untrusted_protocol",
                      self.refused((json.dumps(flat), lines(history(p)))))

    def test_refuses_numbers_out_of_range_without_a_traceback(self):
        p = protocol()
        good = lines(history(p))
        huge = good.replace('"seconds": 1.5', '"seconds": 1' + "0" * 400, 1)
        self.assertIn("bad_field", self.refused((json.dumps(p), huge)))
        self.assertIn("bad_scores", self.refused(
            (json.dumps(p), good.replace("0.72", "1e999", 1))))
        self.assertIn("bad_protocol", self.refused(
            (json.dumps(p).replace('"fit": 1', '"fit": NaN'), good)))
        self.assertIn("bad_json", self.refused(
            (json.dumps(p), "[" * 100000 + "]" * 100000 + "\n")))

    def test_refuses_a_protocol_with_a_repeated_key(self):
        # Plain json.load would keep the last value, which matches the hash.
        p = protocol()
        text = json.dumps(p).replace('"trials_per_family": 3',
                                     '"trials_per_family": 9, '
                                     '"trials_per_family": 3')
        self.assertIn("bad_protocol",
                      self.refused((text, lines(history(p)))))

    def test_arithmetic_errors_are_refusals(self):
        # A backstop: validation already refuses every overflow it knows of.
        p = protocol()
        real = report.build_report
        report.build_report = lambda runs: 1 / 0
        try:
            err = self.refused((json.dumps(p), lines(history(p))))
        finally:
            report.build_report = real
        self.assertIn("refused: ZeroDivisionError", err)

    def test_refuses_the_same_history_twice(self):
        p = protocol()
        run_ = (json.dumps(p), lines(history(p)))
        err = self.refused(run_, run_)
        self.assertIn("duplicate_experiment", err)
        self.assertIn("duplicate_history", err)

    def test_refuses_names_that_are_not_utf8(self):
        # The re-check's F3 (a label byte 0xff reaches Python as \udcff)
        # and T2e (a family name holding a lone surrogate).
        p = protocol()
        err = self.refused((json.dumps(p), lines(history(p))),
                           label="lab\udcff")
        self.assertIn("bad_name", err)
        p = protocol(names=("cart", "\ud800"))
        trials = [trial(p, 0, [0.7, 0.8]), trial(p, 0, [0.6, 0.7],
                                                 family="\ud800")]
        self.assertIn("bad_name", self.refused((json.dumps(p), lines(trials))))

    def test_a_failed_write_leaves_nothing(self):
        p = protocol()
        err = self.refused((json.dumps(p), lines(history(p))),
                           out_md=str(self.dir / "missing" / "r.md"))
        self.assertIn("refused: FileNotFoundError", err)


def lines(trials):
    return "".join(json.dumps(t) + "\n" for t in trials)


if __name__ == "__main__":
    unittest.main()
