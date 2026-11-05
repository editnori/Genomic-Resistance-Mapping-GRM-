"""scripts/hyperopt/report.py on small hand-built histories.

usage: python3 -I -B tests/hyperopt/test_report.py
"""
import copy
import hashlib
import importlib.util
import io
import json
import os
import re
import shutil
import sys
import unittest
from contextlib import redirect_stderr
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


report = _load("grm_hyperopt_report", ROOT / "scripts/hyperopt/report.py")
select = report.select

# Two outer folds' plans, each with inner folds 0 and 1.
PLAN_A = {"inner_folds": 2, "outer_fold": 0}
PLAN_B = {"inner_folds": 2, "outer_fold": 1}


def protocol(plan=PLAN_A, constraints=(), names=("cart", "lasso")):
    """Task B's protocol shape, hashed as Task B hashes it. With sampler
    seed 3, cart's planned trials 1 and 2 have depth 7 and 8."""
    p = {"schema_version": 1, "fold_plan_hash": rehash(plan),
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


def planned(p, family, tid):
    return next(c for f, t, c in select.planned_trials(p)
                if (f, t) == (family, tid))


def trial(p, tid, scores, family="cart", status="ok", reason=None,
          config=None, complexity=1.0, cost=(None, 1.5)):
    """cost is (fit_count, seconds); fit_count None means one fit per fold.
    The config is the planned one unless given."""
    fits, seconds = cost
    t = {"protocol_hash": p["protocol_hash"], "trial_id": tid,
         "family": family, "config": config or planned(p, family, tid),
         "complexity": complexity, "is_default": tid == 0, "seed": 1,
         "fold_plan_hash": p["fold_plan_hash"], "status": status,
         "inner_scores": [{"fold": f, "metrics": {"auc": s}}
                          for f, s in enumerate(scores)],
         "fit_count": len(scores) if fits is None else fits,
         "seconds": seconds}
    if reason is not None:
        t["reason"] = reason
    return t


def history(p, best=1):
    """The whole plan in run order; cart trial `best` (0 or 1) scores
    0.80 and 0.82, the other 0.70 and 0.72."""
    hi, lo = [0.80, 0.82], [0.70, 0.72]
    return [
        trial(p, 0, hi if best == 0 else lo),
        trial(p, 0, [0.60, 0.62], family="lasso"),
        trial(p, 1, hi if best == 1 else lo, complexity=3),
        trial(p, 1, [0.90, 0.90], family="lasso", status="unconverged",
              reason="path stopped at lambda 7"),
        # Failed records as Task B's search writes them.
        trial(p, 2, [], status="failed",
              reason="tree | split crashed\nat node 4", complexity=None,
              cost=(0, None)),
        trial(p, 2, [], family="lasso", status="failed",
              reason="solver diverged", complexity=None, cost=(0, None)),
    ]


def run(label, plan=PLAN_A, best=1, constraints=()):
    p = protocol(plan, constraints)
    return {"label": label, "protocol": p, "trials": history(p, best),
            "fold_plan": plan}


class Report(unittest.TestCase):
    def setUp(self):
        self.report = report.build_report(
            [run("outer0"), run("outer1", PLAN_B, best=0)])
        self.md = report.render_markdown(self.report)

    def test_every_trial_is_listed_with_status_and_reason(self):
        for r, h in zip(self.report["runs"],
                        [history(protocol()), history(protocol(PLAN_B), 0)]):
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
        self.assertEqual(self.md.count("### All trials (6)"), 2)
        all_rows = [line for line in self.md.splitlines()
                    if line.startswith("| cart |") or line.startswith("| lasso |")]
        # Per run: 6 trial rows, 3 excluded rows, 2 selection rows, 2 budgets
        # and 1 row of unused fold scores.
        self.assertEqual(len(all_rows), 2 * (6 + 3 + 2 + 2 + 1))

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
        # Built from records, not read from a file: no file digest.
        self.assertIsNone(r["history_file_sha256"])
        self.assertIn(f"canonical history sha256 {r['history_sha256']}, "
                      "history file sha256 not given,", self.md)

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
                         {"unequal_ok"})
        self.assertIn("- unequal_ok: ok trials differ across families: "
                      "cart 2, lasso 1", self.md)
        selected = [t for t in r["trials"] if t["selected"]]
        self.assertEqual({(t["family"], t["trial_id"]) for t in selected},
                         {("cart", 1), ("lasso", 0)})

    def test_stability_across_runs(self):
        # outer0 selects cart trial 1 (planned depth 7), outer1 the default.
        cart = self.report["stability"]["families"]["cart"]
        self.assertEqual((cart["selected_in"], cart["default_selected"]), (2, 1))
        self.assertEqual(cart["params"]["depth"]["value_counts"],
                         [{"value": 0, "count": 1, "fraction": 0.5},
                          {"value": 7, "count": 1, "fraction": 0.5}])
        self.assertEqual(self.report["stability"]["families"]["lasso"]
                         ["default_selected"], 2)
        self.assertIn("## Stability across runs (2)", self.md)
        self.assertIn("| depth | 2 | 2 | 0.5 | 0 | 3.5 | 7 | 0 |", self.md)

    def test_scores_are_labelled_inner_tuning_never_test_performance(self):
        # Review F3: the labels change, the values do not.
        self.assertEqual(self.report["score_label"], select.SCORE_LABEL)
        self.assertIn("inner-tuning score", select.SCORE_LABEL)
        self.assertIn("None is test performance", select.SCORE_LABEL)
        self.assertEqual(self.report["outer_evaluation"], select.OUTER_EVALUATION)
        self.assertIn("untouched outer test fold", select.OUTER_EVALUATION)
        self.assertIn("remaining dependency", select.OUTER_EVALUATION)
        lines = self.md.splitlines()
        self.assertIn(f"Scores: {select.SCORE_LABEL}", lines)
        self.assertIn(f"Outer evaluation: {select.OUTER_EVALUATION}", lines)
        self.assertIn("| family | selected trial | complexity | default "
                      "| inner-tuning mean | best trial | best inner-tuning mean "
                      "| best inner-tuning SE (descriptive) | candidate set |",
                      lines)
        self.assertIn("| complexity | default | inner-tuning mean "
                      "| inner-tuning SE (descriptive) | candidate | selected |",
                      self.md)
        for bare in ("| mean |", "| best mean |", "| SE (descriptive) |"):
            self.assertNotIn(bare, self.md)
        cart = self.report["runs"][0]["families"]["cart"]
        mean, _, se = select.fold_stats([0.80, 0.82])
        self.assertEqual((cart["selected"]["inner_tuning_mean"],
                          cart["best"]["inner_tuning_mean"],
                          cart["best"]["inner_tuning_se"]), (mean, mean, se))
        self.assertIn("| cart | 1 | 3 |  | 0.81 | 1 | 0.81 | 0.01 | 1 |",
                      self.md)

    def test_json_score_keys_say_inner_tuning(self):
        # REVIEW-HPO L3: a mean, SD or SE read on its own keeps its label.
        def keys(value):
            if isinstance(value, dict):
                return set(value) | set().union(*map(keys, value.values()))
            if isinstance(value, list):
                return set().union(*map(keys, value))
            return set()

        runs = keys(self.report["runs"])
        self.assertFalse({"mean", "sd", "se"} & runs)
        self.assertTrue({"inner_tuning_mean", "inner_tuning_sd",
                         "inner_tuning_se"} <= runs)
        cart = self.report["runs"][0]["families"]["cart"]
        row = next(t for t in self.report["runs"][0]["trials"]
                   if (t["family"], t["trial_id"]) == ("cart", 1))
        mean, sd, se = select.fold_stats([0.80, 0.82])
        self.assertEqual((row["inner_tuning_mean"], row["inner_tuning_sd"],
                          row["inner_tuning_se"]), (mean, sd, se))
        self.assertEqual(cart["candidate_set"][0]["inner_tuning_mean"], mean)
        # select.py's own result keeps its keys; only the report renames.
        s = select.select_history(run("outer0")["trials"], protocol(), PLAN_A)
        self.assertEqual(s["families"]["cart"]["selected"]["mean"], mean)
        # Stability's mean is of hyperparameter values, not a score.
        self.assertIn("mean", self.report["stability"]["families"]["cart"]
                      ["params"]["depth"]["numeric"])

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
        self.assertIn("### Excluded trials (5)", md)
        self.assertIn("unequal_feasible", {f["code"] for f in
                                           rep["runs"][0]["budgets"]["flags"]})

    def test_build_report_computes_the_file_digest(self):
        # REVIEW-HPO I2: the digest is of the bytes the records are parsed
        # from, computed here; a caller's own digest is never shown.
        p = protocol()
        data = canonical_lines(history(p)).encode()
        r = report.build_report([{"label": "outer0", "protocol": p,
                                  "history_bytes": data,
                                  "fold_plan": PLAN_A}])["runs"][0]
        digest = hashlib.sha256(data).hexdigest()
        self.assertEqual((r["history_sha256"], r["history_file_sha256"]),
                         (digest, digest))
        for extra, code in (({"history_file_sha256": "0" * 64}, "untrusted_digest"),
                            ({"history_file_sha256": digest}, "untrusted_digest"),
                            ({"history_bytes": data}, "bad_run")):
            with self.subTest(code=code):
                with self.assertRaises(report.HistoryError) as cm:
                    report.build_report([dict(run("outer0"), **extra)])
                self.assertEqual(cm.exception.codes, {code})

    def test_a_search_over_a_parameter_production_lacks_is_labelled(self):
        # REVIEW-HPO I1: a protocol from before review F1 still searches the
        # lasso's lambda_min_ratio. It verifies and selects as before, and
        # both outputs say production can't set that parameter.
        p = protocol()
        p["families"]["lasso"] = {
            "default": {"lambda_min_ratio": 0.001, "n_lambdas": 30},
            "space": {"lambda_min_ratio": {"type": "float", "low": 0.0001,
                                           "high": 0.1, "log": True},
                      "n_lambdas": {"type": "int", "low": 10, "high": 100,
                                    "log": False}}}
        p["protocol_hash"] = rehash(p)
        trials = history(p)
        s = select.select_history(trials, p, PLAN_A)
        rep = report.build_report([{"label": "pre-F1", "protocol": p,
                                    "trials": trials, "fold_plan": PLAN_A}])
        r = rep["runs"][0]
        self.assertEqual(s["warnings"], r["warnings"])
        self.assertEqual([(w["code"], w["family"]) for w in r["warnings"]],
                         [("no_production_knob", "lasso")])
        self.assertIn("lasso searches lambda_min_ratio: production fixes",
                      r["warnings"][0]["message"])
        self.assertIn(f"Warnings:\n- no_production_knob: "
                      f"{r['warnings'][0]['message']}\n",
                      report.render_markdown(rep))
        self.assertEqual({f: v["selected"]["trial_id"]
                          for f, v in r["families"].items()},
                         {"cart": 1, "lasso": 0})
        # Today's protocols carry no warning.
        self.assertEqual(self.report["runs"][0]["warnings"], [])
        self.assertNotIn("Warnings:", self.md)

    def test_report_is_strict_json(self):
        text = json.dumps(self.report, allow_nan=False)
        # Schema 2 renamed the score keys (INNER_TUNING).
        self.assertEqual(json.loads(text)["schema"], "grm-hyperopt-report/2")

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
        """Run the CLI on (protocol text, history text) pairs, each with
        PLAN_A as its fold plan unless a third item gives the plan's text;
        returns the exit code and stderr."""
        args = []
        for i, (proto, hist, *plan) in enumerate(runs):
            args += ["--run", f"{label}{i}", self.write(f"p{i}.json", proto),
                     self.write(f"h{i}.jsonl", hist),
                     self.write(f"f{i}.json", plan[0] if plan
                                else json.dumps(PLAN_A))]
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
        rc, err = self.main((json.dumps(p), canonical_lines(history(p))))
        self.assertEqual((rc, err), (0, ""))
        out = json.loads((self.dir / "r.json").read_text())
        self.assertEqual(len(out["runs"][0]["trials"]), 6)
        md = (self.dir / "r.md").read_text()
        self.assertIn("SE note:", md)
        self.assertEqual(sorted(p.name for p in self.dir.iterdir()),
                         ["f0.json", "h0.jsonl", "p0.json", "r.json", "r.md"])
        # A canonical file: its sha256 is the canonical history digest.
        digest = hashlib.sha256((self.dir / "h0.jsonl").read_bytes()).hexdigest()
        run = out["runs"][0]
        self.assertEqual((run["history_sha256"], run["history_file_sha256"]),
                         (digest, digest))
        self.assertIn(f"canonical history sha256 {digest}, history file "
                      f"sha256 {digest},", md)

    def test_reports_both_digests_of_a_non_canonical_history(self):
        # Review F9: the file's sha256 and the canonical digest differ, and
        # the report names each for what it is.
        p = protocol()
        rc, err = self.main((json.dumps(p), lines(history(p))))
        self.assertEqual((rc, err), (0, ""))
        run = json.loads((self.dir / "r.json").read_text())["runs"][0]
        file_digest = hashlib.sha256(
            (self.dir / "h0.jsonl").read_bytes()).hexdigest()
        canonical = hashlib.sha256(
            canonical_lines(history(p)).encode()).hexdigest()
        self.assertNotEqual(file_digest, canonical)
        self.assertEqual((run["history_sha256"], run["history_file_sha256"]),
                         (canonical, file_digest))
        self.assertIn(f"canonical history sha256 {canonical}, history file "
                      f"sha256 {file_digest},",
                      (self.dir / "r.md").read_text())

    def test_file_digest_is_of_the_bytes(self):
        # REVIEW-HPO L1 (RV7): CRLF endings parse to the canonical records,
        # so only a digest of the bytes, not of the text, tells them apart.
        p = protocol()
        text = canonical_lines(history(p))
        rc, err = self.main((json.dumps(p), text.replace("\n", "\r\n")))
        self.assertEqual((rc, err), (0, ""))
        data = (self.dir / "h0.jsonl").read_bytes()
        self.assertIn(b"\r\n", data)
        run = json.loads((self.dir / "r.json").read_text())["runs"][0]
        self.assertEqual((run["history_sha256"], run["history_file_sha256"]),
                         (hashlib.sha256(text.encode()).hexdigest(),
                          hashlib.sha256(data).hexdigest()))

    def test_the_history_is_read_once(self):
        # REVIEW-HPO I2: the file is replaced on disk as soon as it is first
        # opened. The open stream still reads the first file, so one read
        # gives that file's records and digest; a second read would pair
        # them with the other file's digest.
        p = protocol()
        first = canonical_lines(history(p, best=1)).encode()
        path = self.dir / "h0.jsonl"
        path.write_bytes(first)
        (self.dir / "next.jsonl").write_bytes(
            canonical_lines(history(p, best=0)).encode())
        real, swapped = io.open, []

        def replacing(file, *args, **kwargs):
            stream = real(file, *args, **kwargs)
            mode = args[0] if args else kwargs.get("mode", "r")
            if os.fspath(file) == str(path) and "r" in mode and not swapped:
                swapped.append(mode)
                os.replace(self.dir / "next.jsonl", path)
            return stream

        argv = ["--run", "outer0", self.write("p0.json", json.dumps(p)),
                str(path), self.write("f0.json", json.dumps(PLAN_A)),
                "--out-json", str(self.dir / "r.json"),
                "--out-md", str(self.dir / "r.md")]
        with mock.patch("builtins.open", replacing), \
                mock.patch("io.open", replacing):
            self.assertEqual(report.main(argv), 0)
        self.assertEqual(len(swapped), 1)
        run = json.loads((self.dir / "r.json").read_text())["runs"][0]
        digest = hashlib.sha256(first).hexdigest()
        self.assertEqual((run["history_sha256"], run["history_file_sha256"]),
                         (digest, digest))
        self.assertEqual(run["families"]["cart"]["selected"]["trial_id"], 1)

    def test_refuses_histories_the_protocol_never_produced(self):
        # Review F2 and F8, through the CLI: an off-plan config with perfect
        # scores, a fold subset, a truncated history and one out of order.
        p = protocol()
        h = history(p)
        forged = copy.deepcopy(h)
        forged[2].update(config={"depth": 3},
                         inner_scores=[{"fold": f, "metrics": {"auc": 1.0}}
                                       for f in (0, 1)])
        subset = copy.deepcopy(h)
        for t in subset:
            t["inner_scores"] = t["inner_scores"][:1]
        cases = {"off_plan": forged, "off_plan_folds": subset,
                 "incomplete_history": h[:4],
                 "off_plan_order": [h[1], h[0]] + h[2:]}
        for code, trials in cases.items():
            with self.subTest(code):
                err = self.refused((json.dumps(p), lines(trials)))
                self.assertEqual(set(re.findall(r"^  ([a-z_]+):", err, re.M)),
                                 {code})

    def test_refuses_a_fold_plan_the_protocol_does_not_pin(self):
        p = protocol()
        good = lines(history(p))
        err = self.refused((json.dumps(p), good, json.dumps(PLAN_B)))
        self.assertIn(f"bad_fold_plan: the fold plan hashes to {rehash(PLAN_B)}",
                      err)
        for text in ('{"inner_folds": 2, "inner_folds": 2, "outer_fold": 0}',
                     '{"inner_folds": NaN}', "{"):
            with self.subTest(text):
                self.assertIn("bad_fold_plan",
                              self.refused((json.dumps(p), good, text)))

    def test_the_fold_plan_argument_is_required(self):
        p = protocol()
        err = io.StringIO()
        with redirect_stderr(err), self.assertRaises(SystemExit) as cm:
            report.main(["--run", "outer0", self.write("p.json", json.dumps(p)),
                         self.write("h.jsonl", lines(history(p))),
                         "--out-json", str(self.dir / "r.json"),
                         "--out-md", str(self.dir / "r.md")])
        self.assertEqual(cm.exception.code, 2)
        self.assertFalse((self.dir / "r.json").exists())

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
        trials = [trial(p, t, [0.7, 0.8], family=f)
                  for f, t, _ in select.planned_trials(p)]
        self.assertIn("bad_name", self.refused((json.dumps(p), lines(trials))))

    def test_a_failed_write_leaves_nothing(self):
        p = protocol()
        err = self.refused((json.dumps(p), lines(history(p))),
                           out_md=str(self.dir / "missing" / "r.md"))
        self.assertIn("refused: FileNotFoundError", err)


def lines(trials):
    """JSON lines with spaces after the separators: not Task B's canonical
    form."""
    return "".join(json.dumps(t) + "\n" for t in trials)


def canonical_lines(trials):
    return "".join(json.dumps(t, sort_keys=True, separators=(",", ":")) + "\n"
                   for t in trials)


if __name__ == "__main__":
    unittest.main()
