"""scripts/hyperopt/select.py on small hand-built histories.

usage: python3 -I -B tests/hyperopt/test_select.py
"""
import copy
import hashlib
import importlib.util
import json
import math
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


sel = _load("grm_hyperopt_select", ROOT / "scripts/hyperopt/select.py")

PH, FH = "proto-1", "plan-1"


def flat(direction="maximize", budgets=None, **kw):
    """The flat, test-only protocol shape."""
    return dict({"objective": "auc", "direction": direction,
                 "protocol_hash": PH, "fold_plan_hash": FH,
                 "budgets": budgets or {"cart": 4}}, **kw)


def protocol(*args, **kw):
    return sel.Protocol.from_dict(flat(*args, **kw), trusted_flat=True)


def b_hash(p, ensure_ascii=True):
    """Task B's protocol_hash, written out independently of select.py."""
    body = {k: v for k, v in p.items() if k != "protocol_hash"}
    return hashlib.sha256(json.dumps(
        body, sort_keys=True, separators=(",", ":"), ensure_ascii=ensure_ascii,
        allow_nan=False).encode("utf-8")).hexdigest()


def b_protocol(constraints=(), direction="maximize", fit=7, n=4, plan=FH,
               names=("cart",), **kw):
    """Task B's protocol shape with a true protocol_hash. Each family's
    default is depth 0; its space is depth 0..9, lr 0.0..1.0 (float) and
    kind a or b."""
    p = dict({
        "schema_version": 1, "fold_plan_hash": plan,
        "objective": {"name": "auc", "direction": direction,
                      "constraints": [dict(c) for c in constraints]},
        "seeds": {"sampler": 1, "fit": fit}, "trials_per_family": n,
        "families": {f: {"default": {"depth": 0, "kind": "a", "lr": 0.5},
                         "space": {
                             "depth": {"type": "int", "low": 0, "high": 9,
                                       "log": False},
                             "kind": {"type": "choice", "values": ["a", "b"]},
                             "lr": {"type": "float", "low": 0.0, "high": 1.0,
                                    "log": False}}}
                     for f in names}}, **kw)
    p["protocol_hash"] = b_hash(p)
    return p


def b_trial(p, tid, scores, depth=None, **kw):
    """A record that belongs to Task B protocol p."""
    fields = {"config": {"depth": tid if depth is None else depth,
                         "kind": "a", "lr": 0.5},
              "protocol_hash": p["protocol_hash"], "seed": p["seeds"]["fit"],
              "fold_plan_hash": p["fold_plan_hash"]}
    return trial(tid, scores, **dict(fields, **kw))


def trial(tid, scores, family="cart", complexity=1.0, status="ok",
          reason=None, config=None, folds=None, acc=None, **kw):
    """One trial record; scores are the objective (auc) on folds 0..k-1, and
    acc is 0.5 on every fold unless given. A trial that is not ok gets a
    reason unless one is given."""
    folds = folds if folds is not None else range(len(scores))
    acc = acc if acc is not None else [0.5] * len(scores)
    t = {
        "protocol_hash": PH, "trial_id": tid, "family": family,
        "config": config if config is not None else {"depth": tid},
        "complexity": complexity, "is_default": tid == 0, "seed": 7,
        "fold_plan_hash": FH, "status": status,
        "inner_scores": [{"fold": f, "metrics": {"auc": s, "acc": a}}
                         for f, s, a in zip(folds, scores, acc)],
        "fit_count": len(scores), "seconds": 0.5,
    }
    if status != "ok":
        t["reason"] = reason if reason is not None else f"test {status}"
    t.update(kw)
    return t


def family(result, name="cart"):
    return result["families"][name]


class Direction(unittest.TestCase):
    def test_maximize_picks_the_highest_mean(self):
        trials = [trial(0, [1, 1]), trial(1, [9, 9]), trial(2, [5, 5])]
        f = family(sel.select_history(trials, protocol("maximize")))
        self.assertEqual(f["best"]["trial_id"], 1)
        self.assertEqual(f["selected"]["trial_id"], 1)

    def test_minimize_picks_the_lowest_mean(self):
        trials = [trial(0, [1, 1]), trial(1, [9, 9]), trial(2, [5, 5])]
        f = family(sel.select_history(trials, protocol("minimize")))
        self.assertEqual(f["best"]["trial_id"], 0)
        self.assertEqual(f["selected"]["trial_id"], 0)

    def test_the_se_window_opens_toward_worse_scores(self):
        # Best mean 5, SE 1. Maximize: the window is [4, 5]; minimize, with
        # the best at 5 too, it is [5, 6].
        up = family(sel.select_history(
            [trial(0, [4, 6], complexity=3), trial(1, [4, 4], complexity=1)],
            protocol("maximize")))
        self.assertEqual(up["threshold"], 4.0)
        self.assertEqual([c["trial_id"] for c in up["candidate_set"]], [1, 0])
        down = family(sel.select_history(
            [trial(0, [4, 6], complexity=3), trial(2, [6, 6], complexity=1)],
            protocol("minimize")))
        self.assertEqual(down["threshold"], 6.0)
        self.assertEqual([c["trial_id"] for c in down["candidate_set"]], [2, 0])


class Exclusions(unittest.TestCase):
    def test_failed_and_unconverged_are_excluded_with_reasons(self):
        trials = [trial(0, [0.6, 0.6]),
                  trial(1, [0.99, 0.99], status="unconverged",
                        reason="max_iter reached at lambda 3"),
                  trial(2, [], status="failed", reason="out of memory",
                        complexity=None, fit_count=None, seconds=None),
                  trial(3, [0.9], status="failed", reason="fold 1 crashed")]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["n_trials"], 4)
        self.assertEqual(f["n_ok"], 1)
        self.assertEqual(f["selected"]["trial_id"], 0)
        self.assertEqual(f["excluded"], [
            {"trial_id": 1, "status": "unconverged",
             "reason": "max_iter reached at lambda 3",
             "excluded_because": ["status is unconverged"]},
            {"trial_id": 2, "status": "failed", "reason": "out of memory",
             "excluded_because": ["status is failed"]},
            {"trial_id": 3, "status": "failed", "reason": "fold 1 crashed",
             "excluded_because": ["status is failed"]},
        ])
        rows = {r["trial_id"]: r for r in f["trials"]}
        self.assertTrue(rows[1]["excluded"])
        self.assertIsNone(rows[1]["mean"])
        self.assertFalse(rows[1]["in_candidate_set"])

    def test_a_family_with_no_ok_trial_selects_nothing(self):
        trials = [trial(0, [], status="failed", reason="boom")]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["status"], "no_ok_trials")
        self.assertIsNone(f["selected"])
        self.assertEqual(len(f["excluded"]), 1)

    def test_a_budgeted_family_with_no_trials(self):
        trials = [trial(0, [1, 2])]
        r = sel.select_history(trials, protocol(budgets={"cart": 1, "rf": 1}))
        self.assertEqual(family(r, "rf")["status"], "no_trials")


class OneStandardError(unittest.TestCase):
    def test_exact_se_on_a_hand_built_example(self):
        # Scores 1, 3, 5, 7: mean 4, squared deviations 9+1+1+9 = 20,
        # sample variance 20/3, SD sqrt(20/3), SE sqrt(20/3)/2.
        mean, sd, se = sel.fold_stats([1, 3, 5, 7])
        self.assertEqual(mean, 4.0)
        self.assertAlmostEqual(sd, math.sqrt(20 / 3), places=15)
        self.assertAlmostEqual(se, math.sqrt(20 / 3) / 2, places=15)
        # Two folds: SE is half the gap, exactly.
        self.assertEqual(sel.fold_stats([4, 6]), (5.0, math.sqrt(2), 1.0))

    def test_threshold_uses_the_best_trials_se_and_is_inclusive(self):
        # Best: trial 1, mean 5, SE 1, so the maximize threshold is 4.
        # Trial 2 has mean exactly 4 (in); trial 3 has 3.9 (out), though its
        # own SE is large: only the best trial's SE counts.
        trials = [trial(0, [0, 0], complexity=9),
                  trial(1, [4, 6], complexity=5),
                  trial(2, [4, 4], complexity=3),
                  trial(3, [-2.1, 9.9], complexity=1)]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["best"], {"trial_id": 1, "mean": 5.0, "se": 1.0})
        self.assertEqual(f["threshold"], 4.0)
        self.assertEqual([c["trial_id"] for c in f["candidate_set"]], [2, 1])
        self.assertEqual(f["selected"]["trial_id"], 2)
        self.assertFalse(f["selected"]["is_best"])

    def test_boundary_survives_rounding(self):
        # 0.1 + 0.2 style rounding must not drop a trial on the boundary:
        # best mean 0.3, SE 0.1, candidate mean 0.2.
        trials = [trial(0, [0.2, 0.4], complexity=2),
                  trial(1, [0.2, 0.2], complexity=1)]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["selected"]["trial_id"], 1)

    def test_lower_complexity_beats_the_default(self):
        trials = [trial(0, [5, 5], complexity=2), trial(1, [4, 6], complexity=1)]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual([c["trial_id"] for c in f["candidate_set"]], [1, 0])
        self.assertEqual(f["selected"]["trial_id"], 1)

    def test_the_default_wins_a_complexity_tie(self):
        trials = [trial(0, [4.5, 4.5], complexity=2),
                  trial(1, [4, 6], complexity=2), trial(2, [5, 5], complexity=2)]
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["selected"]["trial_id"], 0)
        # The schema makes trial 0 the default, so in a history the default
        # also has the lowest id; the key itself ranks default before id.
        key = sel.preference_key
        a = {"complexity": 2, "is_default": True, "trial_id": 9}
        b = {"complexity": 2, "is_default": False, "trial_id": 1}
        c = {"complexity": 1, "is_default": False, "trial_id": 5}
        self.assertEqual(sorted([a, b, c], key=key), [c, a, b])

    def test_lowest_id_breaks_the_last_tie(self):
        trials = [trial(0, [0, 0], complexity=1), trial(5, [4, 6], complexity=2),
                  trial(3, [4, 6], complexity=2), trial(4, [5, 5], complexity=3)]
        f = family(sel.select_history(trials, protocol(budgets={"cart": 6})))
        self.assertEqual([c["trial_id"] for c in f["candidate_set"]], [3, 5, 4])
        self.assertEqual(f["selected"]["trial_id"], 3)

    def test_equal_best_means_resolve_by_preference(self):
        trials = [trial(1, [4, 6], complexity=2), trial(2, [5, 5], complexity=1)]
        trials.insert(0, trial(0, [0, 0]))
        f = family(sel.select_history(trials, protocol()))
        self.assertEqual(f["best"]["trial_id"], 2)
        self.assertEqual(f["best"]["se"], 0.0)
        self.assertEqual([c["trial_id"] for c in f["candidate_set"]], [2, 1])

    def test_se_label_is_carried(self):
        r = sel.select_history([trial(0, [1, 2])], protocol())
        for text in (r["se_label"], family(r)["se_label"]):
            self.assertIn("descriptive", text)
            self.assertIn("not an independent-sample inference", text)
        self.assertIn("descriptive", sel.select_history.__doc__ + sel.__doc__)


class Refusals(unittest.TestCase):
    def refused(self, trials, proto=None):
        """The refusal codes; a dict protocol is read as a trusted flat one
        inside the check, so its own refusal counts."""
        with self.assertRaises(sel.HistoryError) as cm:
            p = proto or protocol()
            if isinstance(p, dict) and not isinstance(p.get("objective"), dict):
                p = sel.Protocol.from_dict(p, trusted_flat=True)
            sel.select_history(trials, p)
        return cm.exception.codes

    def test_protocol_hash_mismatch(self):
        codes = self.refused([trial(0, [1, 2], protocol_hash="other")])
        self.assertEqual(codes, {"protocol_hash_mismatch"})

    def test_fold_plan_hash_mismatch(self):
        codes = self.refused([trial(0, [1, 2], fold_plan_hash="other")])
        self.assertEqual(codes, {"fold_plan_hash_mismatch"})

    def test_mixed_hashes(self):
        codes = self.refused([trial(0, [1, 2]),
                              trial(1, [1, 2], protocol_hash="p2",
                                    fold_plan_hash="f2")])
        self.assertTrue({"mixed_protocol_hash",
                         "mixed_fold_plan_hash"} <= codes)

    def test_duplicate_trial_ids(self):
        codes = self.refused([trial(0, [1, 2]), trial(1, [1, 2]),
                              trial(1, [3, 4])])
        self.assertEqual(codes, {"duplicate_trial_id"})

    def test_same_id_in_two_families_is_fine(self):
        r = sel.select_history([trial(0, [1, 2]), trial(0, [1, 2], family="rf")],
                               protocol(budgets={"cart": 1, "rf": 1}))
        self.assertEqual(set(r["families"]), {"cart", "rf"})

    def test_failed_trial_folds_must_be_within_the_plan(self):
        self.assertEqual(self.refused([trial(0, [1, 2]),
                                       trial(1, [1], status="failed",
                                             folds=[5])]), {"fold_mismatch"})

    def test_missing_fields_and_bad_types(self):
        t = trial(0, [1, 2])
        del t["seconds"]
        self.assertEqual(self.refused([t]), {"missing_field"})
        self.assertEqual(self.refused([trial(0, [1, 2], trial_id=True)]),
                         {"bad_field"})
        self.assertEqual(self.refused([trial(0, [1, 2], status="done")]),
                         {"bad_field"})
        self.assertEqual(self.refused([trial(0, [1, 2], complexity=math.nan)]),
                         {"bad_field"})
        self.assertEqual(self.refused([trial(0, [1, 2], complexity=-1)]),
                         {"bad_field"})

    def test_missing_objective_and_one_fold(self):
        t = trial(0, [1, 2])
        del t["inner_scores"][1]["metrics"]["auc"]
        self.assertEqual(self.refused([t]), {"missing_metric"})
        self.assertEqual(self.refused([trial(0, [1])]), {"too_few_folds"})

    def test_missing_constraint_metric(self):
        t = trial(0, [1, 2])
        del t["inner_scores"][0]["metrics"]["acc"]
        proto = protocol(constraints=[{"metric": "acc", "op": ">=",
                                       "bound": 0.1}])
        self.assertEqual(self.refused([t], proto), {"missing_metric"})
        # A trial-field constraint needs no inner metric.
        proto = protocol(constraints=[{"metric": "seconds", "op": "<=",
                                       "bound": 9}])
        sel.select_history([t], proto)

    def test_unsorted_and_repeated_folds(self):
        self.assertEqual(self.refused([trial(0, [1, 2], folds=[1, 0])]),
                         {"unsorted_folds"})
        self.assertEqual(self.refused([trial(0, [1, 2], folds=[0, 0])]),
                         {"duplicate_fold"})

    def test_reason_present_exactly_when_not_ok(self):
        t = trial(0, [1, 2])
        t["reason"] = None
        self.assertEqual(self.refused([t]), {"bad_field"})
        t = trial(0, [1, 2], status="failed")
        del t["reason"]
        self.assertEqual(self.refused([t]), {"bad_field"})
        self.assertEqual(self.refused([trial(0, [1, 2], status="unconverged",
                                             reason="")]), {"bad_field"})

    def test_null_costs_only_on_failed_trials(self):
        nulls = {"complexity": None, "fit_count": None, "seconds": None}
        sel.select_history([trial(0, [1, 2]),
                            trial(1, [], status="failed", **nulls)],
                           protocol())
        for status in ("ok", "unconverged"):
            for field in nulls:
                self.assertEqual(
                    self.refused([trial(0, [1, 2], status=status,
                                        **{field: None})]), {"bad_field"})

    def test_default_must_be_trial_zero(self):
        self.assertEqual(self.refused([trial(0, [1, 2]),
                                       trial(1, [1, 2], is_default=True)]),
                         {"bad_default"})

    def test_unbudgeted_family_and_bad_protocol(self):
        self.assertEqual(self.refused([trial(0, [1, 2], family="gbt")]),
                         {"unbudgeted_family"})
        self.assertEqual(self.refused([trial(0, [1, 2])],
                                      flat(direction="up")), {"bad_protocol"})

    def test_every_problem_is_reported(self):
        codes = self.refused([trial(0, [1, 2], protocol_hash="x"),
                              trial(0, [1, 2])])
        self.assertEqual(codes, {"protocol_hash_mismatch",
                                 "mixed_protocol_hash", "duplicate_trial_id"})

    def test_load_history_refuses_bad_lines(self):
        tmp = ROOT / ".tmp"
        tmp.mkdir(exist_ok=True)
        path = tmp / "test_select_history.jsonl"
        try:
            good = json.dumps(trial(0, [1, 2]))
            path.write_text(good + "\n\n{\"a\":1,\"a\":2}\n{\"a\":NaN}\n")
            with self.assertRaises(sel.HistoryError) as cm:
                sel.load_history(path)
            self.assertEqual(len(cm.exception.problems), 3)
            path.write_text(good + "\n" + good)
            with self.assertRaises(sel.HistoryError) as cm:
                sel.load_history(path)
            self.assertIn("no newline", cm.exception.problems[0])
            path.write_text(good + "\n")
            self.assertEqual(len(sel.load_history(path)), 1)
        finally:
            path.unlink(missing_ok=True)


class Constraints(unittest.TestCase):
    CONSTRAINTS = [{"metric": "acc", "op": ">=", "bound": 0.5},
                   {"metric": "seconds", "op": "<=", "bound": 2}]

    def history(self):
        # Unconstrained, trial 1 (mean 9) would be best. It fails acc >= 0.5
        # (mean acc 0.2), and trial 3 fails seconds <= 2. Trial 2 sits on
        # both bounds, which are inclusive.
        return [trial(0, [5, 5], acc=[0.4, 0.6]),
                trial(1, [9, 9], complexity=3, acc=[0.2, 0.2]),
                trial(2, [8, 8], complexity=2, seconds=2.0),
                trial(3, [8.5, 8.5], seconds=3.0)]

    def test_selection_enforces_constraints(self):
        r = sel.select_history(self.history(),
                               protocol(constraints=self.CONSTRAINTS))
        f = family(r)
        self.assertEqual((f["n_ok"], f["n_feasible"]), (4, 2))
        self.assertEqual(f["best"], {"trial_id": 2, "mean": 8.0, "se": 0.0})
        self.assertEqual(f["selected"]["trial_id"], 2)
        self.assertEqual([c["trial_id"] for c in f["candidate_set"]], [2])
        self.assertEqual(f["excluded"], [
            {"trial_id": 1, "status": "ok", "reason": None,
             "excluded_because": ["constraint acc >= 0.5 not met (acc is 0.2)"]},
            {"trial_id": 3, "status": "ok", "reason": None,
             "excluded_because": [
                 "constraint seconds <= 2 not met (seconds is 3.0)"]},
        ])
        rows = {row["trial_id"]: row for row in f["trials"]}
        self.assertEqual(rows[1]["mean"], 9.0)
        self.assertEqual(rows[1]["constraint_values"],
                         {"acc": 0.2, "seconds": 0.5})
        self.assertEqual(rows[0]["constraint_values"]["acc"], 0.5)
        self.assertFalse(rows[0]["excluded"])
        b = r["budgets"]["families"]["cart"]
        self.assertEqual((b["ok"], b["feasible"], b["infeasible"]), (4, 2, 2))
        # Without constraints, the same history selects trial 1.
        self.assertEqual(family(sel.select_history(
            self.history(), protocol()))["selected"]["trial_id"], 1)

    def test_no_feasible_trial(self):
        f = family(sel.select_history(self.history(), protocol(constraints=[
            {"metric": "auc", "op": ">=", "bound": 10}])))
        self.assertEqual(f["status"], "no_feasible_trials")
        self.assertIsNone(f["selected"])
        self.assertEqual(len(f["excluded"]), 4)

    def test_unequal_feasible_flag_only_with_constraints(self):
        trials = [trial(0, [1, 2]), trial(0, [1, 2], family="rf",
                                          acc=[0.1, 0.1])]
        budgets = {"cart": 1, "rf": 1}
        b = sel.budget_accounting(trials, protocol(
            budgets=budgets, constraints=[self.CONSTRAINTS[0]]))
        self.assertEqual([f["code"] for f in b["flags"]], ["unequal_feasible"])
        self.assertEqual(sel.budget_accounting(
            trials, protocol(budgets=budgets))["flags"], [])

    def test_bad_constraints_are_refused(self):
        for c in [{"metric": "acc", "op": "<", "bound": 1},
                  {"metric": "acc", "op": "<=", "bound": math.inf},
                  {"metric": "acc", "op": "<=", "bound": 1, "extra": 0},
                  {"metric": "", "op": "<=", "bound": 1}]:
            with self.assertRaises(sel.HistoryError) as cm:
                sel.Protocol.from_dict(flat(constraints=[c]),
                                       trusted_flat=True)
            self.assertEqual(cm.exception.codes, {"bad_protocol"})


class TaskBProtocol(unittest.TestCase):
    ACC = [{"metric": "acc", "op": ">=", "bound": 0.5}]

    def test_reads_task_b_protocol(self):
        bp = b_protocol(self.ACC, n=2, names=("cart", "rf"))
        p = sel.Protocol.from_dict(bp)
        self.assertEqual((p.objective, p.direction), ("auc", "maximize"))
        self.assertEqual(p.budgets, {"cart": 2, "rf": 2})
        self.assertEqual(p.constraints,
                         ({"metric": "acc", "op": ">=", "bound": 0.5},))
        self.assertEqual((p.protocol_hash, p.fit_seed), (bp["protocol_hash"], 7))
        r = sel.select_history([b_trial(bp, 0, [1, 2]), b_trial(bp, 1, [3, 3]),
                                b_trial(bp, 0, [1, 2], family="rf")], bp)
        self.assertEqual(family(r)["selected"]["trial_id"], 1)

    def test_refuses_bad_task_b_protocol(self):
        for kw in [{"trials_per_family": 0}, {"families": {}},
                   {"objective": {"name": "auc", "direction": "up",
                                  "constraints": []}},
                   {"seeds": {"sampler": 1}},
                   {"schema_version": 2}, {"schema_version": True},
                   {"schema_version": "1"}]:
            with self.assertRaises(sel.HistoryError) as cm:
                sel.Protocol.from_dict(b_protocol(**kw))
            self.assertEqual(cm.exception.codes, {"bad_protocol"}, kw)
        for spec in [{"type": "int", "low": 3, "high": 1},
                     {"type": "grid", "values": [1]}]:
            with self.assertRaises(sel.HistoryError) as cm:
                sel.Protocol.from_dict(b_protocol(families={"cart": {
                    "default": {"d": 1}, "space": {"d": spec}}}))
            self.assertEqual(cm.exception.problems, [
                "bad_protocol: family cart: space d must be {type: int or "
                "float, low <= high} or {type: choice, values: [...]}"])

    def test_defaults_are_fixed_or_internal_search(self):
        space = b_protocol()["families"]["cart"]["space"]
        grid = {"internal_search": True, "depth": [1, 2, 3], "kind": ["a", "b"],
                "lr": [0.25, 0.5]}
        sel.Protocol.from_dict(b_protocol(families={"cart": {
            "default": grid, "space": space}}))
        bad = [{"depth": 10, "kind": "a", "lr": 0.5},          # fixed, outside
               dict(grid, internal_search=False),              # not a flag
               {k: v for k, v in grid.items() if k != "lr"},   # lr missing
               dict(grid, lr=[]),                              # empty list
               dict(grid, lr=0.5),                             # not a list
               dict(grid, depth=[1, 1]),                       # repeated
               dict(grid, depth=[1, 10]),                      # outside
               dict(grid, depth=[1.0]),                        # int as float
               dict(grid, extra=[1])]                          # unknown key
        for default in bad:
            with self.assertRaises(sel.HistoryError) as cm:
                sel.Protocol.from_dict(b_protocol(families={"cart": {
                    "default": default, "space": space}}))
            self.assertEqual(cm.exception.codes, {"bad_protocol"}, default)
        reserved = dict(space, internal_search={"type": "choice",
                                                "values": [True]})
        with self.assertRaises(sel.HistoryError) as cm:
            sel.Protocol.from_dict(b_protocol(families={"cart": {
                "default": {"depth": 0, "kind": "a", "lr": 0.5,
                            "internal_search": True}, "space": reserved}}))
        self.assertIn("internal_search is reserved", str(cm.exception))

    def test_flat_protocol_needs_trusted_flat(self):
        with self.assertRaises(sel.HistoryError) as cm:
            sel.Protocol.from_dict(flat())
        self.assertEqual(cm.exception.codes, {"untrusted_protocol"})
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history([trial(0, [1, 2])], flat())
        self.assertEqual(cm.exception.codes, {"untrusted_protocol"})
        self.assertEqual(sel.Protocol.from_dict(flat(), trusted_flat=True)
                         .objective, "auc")

    def test_hash_is_task_bs_canonical_form(self):
        bp = b_protocol()
        self.assertEqual(sel.content_hash(bp), bp["protocol_hash"])
        # ensure_ascii matters once a value is not ASCII.
        bp["families"]["cart"]["space"]["kind"]["values"] = ["a", "bé"]
        bp["protocol_hash"] = b_hash(bp, ensure_ascii=False)
        with self.assertRaises(sel.HistoryError) as cm:
            sel.Protocol.from_dict(bp)
        self.assertEqual(cm.exception.codes, {"protocol_hash_mismatch"})
        bp["protocol_hash"] = b_hash(bp)
        sel.Protocol.from_dict(bp)

    def test_edited_direction_with_stale_hash_is_refused(self):
        # Reviewer case 1: flip the direction and keep the old hash.
        bp = b_protocol()
        trials = [b_trial(bp, 0, [1, 1]), b_trial(bp, 1, [9, 9]),
                  b_trial(bp, 2, [5, 5])]
        self.assertEqual(family(sel.select_history(trials, bp))
                         ["selected"]["trial_id"], 1)
        tampered = copy.deepcopy(bp)
        tampered["objective"]["direction"] = "minimize"
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history(trials, tampered)
        self.assertEqual(cm.exception.codes, {"protocol_hash_mismatch"})
        self.assertIn(bp["protocol_hash"], cm.exception.problems[0])

    def test_loosened_constraint_with_stale_hash_is_refused(self):
        # Reviewer case 2 (probe P5): complexity <= 0.5 leaves nothing
        # feasible; loosening the bound to 10 would select trial 1.
        bp = b_protocol([{"metric": "complexity", "op": "<=", "bound": 0.5}])
        trials = [b_trial(bp, 0, [0.7, 0.8]),
                  b_trial(bp, 1, [0.9, 0.95], complexity=3)]
        self.assertEqual(family(sel.select_history(trials, bp))["status"],
                         "no_feasible_trials")
        tampered = copy.deepcopy(bp)
        tampered["objective"]["constraints"][0]["bound"] = 10
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history(trials, tampered)
        self.assertEqual(cm.exception.codes, {"protocol_hash_mismatch"})


class TaskBRecordRules(unittest.TestCase):
    def setUp(self):
        self.p = b_protocol()

    def refused(self, *trials):
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history(list(trials), self.p)
        return cm.exception.codes

    def test_unknown_field(self):
        self.assertEqual(self.refused(b_trial(self.p, 0, [1, 2], extra=1)),
                         {"unknown_field"})
        t = b_trial(self.p, 0, [1, 2])
        t["inner_scores"][0]["note"] = "x"
        self.assertEqual(self.refused(t), {"bad_scores"})

    def test_seed_must_be_the_fit_seed(self):
        self.assertEqual(self.refused(b_trial(self.p, 0, [1, 2], seed=8)),
                         {"bad_seed"})

    def test_config_must_lie_in_the_space(self):
        good = b_trial(self.p, 0, [1, 2])
        for config in [{"depth": 10, "kind": "a", "lr": 0.5},   # above high
                       {"depth": 1.0, "kind": "a", "lr": 0.5},  # int as float
                       {"depth": 1, "kind": "a", "lr": 1},      # float as int
                       {"depth": 1, "kind": "c", "lr": 0.5},    # not a choice
                       {"depth": 1, "kind": "a"},               # missing lr
                       {"depth": 1, "kind": "a", "lr": 0.5, "x": 1}]:
            t = b_trial(self.p, 1, [1, 2], config=config)
            self.assertEqual(self.refused(good, t), {"bad_config"}, config)
        # Both bounds are inclusive.
        sel.select_history([good, b_trial(self.p, 1, [1, 2], config={
            "depth": 9, "kind": "b", "lr": 0.0})], self.p)

    def test_trial_zero_must_be_the_default_which_may_lie_outside_the_space(self):
        self.assertEqual(self.refused(b_trial(self.p, 0, [1, 2], depth=1)),
                         {"bad_config"})
        # Task B's scm and cart defaults run the learner's own search: lists
        # of in-space values plus internal_search, accepted as trial 0 only.
        default = {"internal_search": True, "depth": [1, 2, 3],
                   "kind": ["a", "b"], "lr": [0.5]}
        self.p["families"]["cart"]["default"] = default
        self.p["protocol_hash"] = b_hash(self.p)
        t0 = b_trial(self.p, 0, [1, 2], config=copy.deepcopy(default))
        self.assertEqual(family(sel.select_history([t0], self.p))
                         ["selected"]["config"], default)
        self.assertEqual(self.refused(t0, b_trial(self.p, 1, [1, 2],
                                                  config=copy.deepcopy(default))),
                         {"bad_config"})
        # A sampled trial never carries internal_search.
        self.assertEqual(self.refused(t0, b_trial(self.p, 1, [1, 2], config={
            "internal_search": True, "depth": 1, "kind": "a", "lr": 0.5})),
                         {"bad_config"})
        # Trial 0 must be the default exactly, not another grid.
        self.assertEqual(self.refused(b_trial(self.p, 0, [1, 2], config=dict(
            default, depth=[1, 2]))), {"bad_config"})


class FoldSets(unittest.TestCase):
    def test_ok_trial_on_other_folds_is_excluded_not_refused(self):
        trials = [trial(0, [1, 2, 3]), trial(1, [9, 9], folds=[0, 2]),
                  trial(2, [2, 2, 2], status="unconverged", folds=[0, 1, 2])]
        r = sel.select_history(trials, protocol())
        f = family(r)
        self.assertEqual(r["folds"], [0, 1, 2])
        self.assertEqual(f["selected"]["trial_id"], 0)
        self.assertEqual(f["excluded"][0], {
            "trial_id": 1, "status": "ok", "reason": None,
            "excluded_because": ["fold_set_mismatch: scored on folds [0, 2], "
                                 "not the history's inner folds [0, 1, 2]"]})
        row = f["trials"][1]
        self.assertEqual((row["mean"], row["se"], row["scores_used"]),
                         (None, None, False))
        self.assertEqual(row["fold_scores"], [{"fold": 0, "score": 9},
                                              {"fold": 2, "score": 9}])
        b = r["budgets"]["families"]["cart"]
        self.assertEqual((b["ok"], b["feasible"], b["infeasible"],
                          b["fold_set_mismatch"]), (2, 1, 0, 1))
        self.assertIn({"code": "fold_set_mismatch", "family": "cart",
                       "message": "cart has 1 ok trial(s) scored on folds "
                       "other than the history's [0, 1, 2], excluded"},
                      r["budgets"]["flags"])

    def test_the_fold_set_most_ok_trials_report_wins(self):
        def folds_of(*sets):
            trials = [trial(i, [5] * len(f), folds=f) for i, f in enumerate(sets)]
            return sel.select_history(
                trials, protocol(budgets={"cart": len(sets)}))["folds"]
        # Two trials on [0, 1] outvote one on the larger [0, 1, 2].
        self.assertEqual(folds_of([0, 1, 2], [0, 1], [0, 1]), [0, 1])
        # A tie goes to the larger set, then the smallest fold list.
        self.assertEqual(folds_of([0, 1], [0, 1, 2]), [0, 1, 2])
        self.assertEqual(folds_of([1, 2], [0, 1]), [0, 1])
        self.assertEqual(folds_of([0, 2], [0, 1]), [0, 1])

    def test_one_corrupt_record_does_not_empty_other_families(self):
        # The re-check's F1: one ok rf trial gains fold 9.
        names = ("cart", "gbt", "lasso", "rf", "scm")
        trials = [trial(i, [5 + i, 6 + i, 7 + i], family=f)
                  for f in names for i in range(3)]
        rf2 = next(t for t in trials if t["family"] == "rf" and t["trial_id"] == 2)
        rf2["inner_scores"].append({"fold": 9, "metrics": {"auc": 99, "acc": 0.5}})
        r = sel.select_history(trials, protocol(budgets=dict.fromkeys(names, 3)))
        self.assertEqual(r["folds"], [0, 1, 2])
        self.assertEqual({f: v["selected"]["trial_id"]
                          for f, v in r["families"].items()},
                         {"cart": 2, "gbt": 2, "lasso": 2, "rf": 1, "scm": 2})
        self.assertEqual(family(r, "rf")["excluded"][0]["excluded_because"], [
            "fold_set_mismatch: scored on folds [0, 1, 2, 9], not the "
            "history's inner folds [0, 1, 2]"])
        self.assertEqual([(f["code"], f["family"]) for f in r["budgets"]["flags"]
                          if f["code"].startswith("fold")],
                         [("fold_set_mismatch", "rf")])
        # The same extra fold on a trial that is not ok is refused (F2).
        rf2["status"], rf2["reason"] = "unconverged", "cap"
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history(trials, protocol(budgets=dict.fromkeys(names, 3)))
        self.assertEqual(cm.exception.codes, {"fold_mismatch"})

    def test_partial_scores_of_unconverged_trials_are_kept(self):
        trials = [trial(0, [1, 2]),
                  trial(1, [0.3], status="unconverged", folds=[1]),
                  trial(2, [], status="failed")]
        rows = family(sel.select_history(trials, protocol()))["trials"]
        self.assertEqual([(r["fold_scores"], r["scores_used"]) for r in rows], [
            ([{"fold": 0, "score": 1}, {"fold": 1, "score": 2}], True),
            ([{"fold": 1, "score": 0.3}], False), ([], False)])
        self.assertIsNone(rows[1]["mean"])


class OutOfRange(unittest.TestCase):
    def refused(self, *trials):
        with self.assertRaises(sel.HistoryError) as cm:
            sel.select_history(list(trials), protocol())
        return cm.exception.codes

    def test_integers_too_large_for_a_float(self):
        self.assertEqual(self.refused(trial(0, [1, 2], seconds=10 ** 400)),
                         {"bad_field"})
        self.assertEqual(self.refused(trial(0, [1, 10 ** 400])), {"bad_scores"})
        self.assertEqual(self.refused(trial(0, [1, 2], acc=[0.5, 10 ** 400])),
                         {"bad_scores"})

    def test_finite_values_that_overflow(self):
        self.assertEqual(self.refused(trial(0, [1e308, 1e308])), {"out_of_range"})
        self.assertEqual(self.refused(trial(0, [1.5e308, -1.5e308])),
                         {"out_of_range"})
        self.assertEqual(self.refused(trial(0, [1, 2], seconds=1e308),
                                      trial(1, [1, 2], seconds=1e308)),
                         {"out_of_range"})


class Budgets(unittest.TestCase):
    def test_counts_and_totals(self):
        trials = [trial(0, [1, 2]), trial(1, [1, 2], status="failed"),
                  trial(2, [1, 2], status="unconverged", seconds=2.25),
                  trial(3, [], status="failed", fit_count=None, seconds=None,
                        complexity=None),
                  trial(0, [1, 2], family="rf", fit_count=10)]
        b = sel.budget_accounting(trials, protocol(budgets={"cart": 4, "rf": 1}))
        self.assertEqual(b["families"]["cart"], {
            "planned": 4, "used": 4, "ok": 1, "failed": 2, "unconverged": 1,
            "feasible": 1, "infeasible": 0, "fold_set_mismatch": 0,
            "fit_count": 6, "seconds": 3.25, "cost_unreported": 1})
        self.assertEqual(b["families"]["rf"]["fit_count"], 10)

    def test_flags(self):
        trials = [trial(0, [1, 2]), trial(2, [1, 2]),
                  trial(0, [1, 2], family="rf"), trial(1, [1, 2], family="rf"),
                  trial(2, [1, 2], family="rf")]
        b = sel.budget_accounting(
            trials, protocol(budgets={"cart": 3, "rf": 2, "gbt": 2}))
        flags = {(f["code"], f["family"]) for f in b["flags"]}
        self.assertEqual(flags, {
            ("under_budget", "cart"), ("trial_id_gap", "cart"),
            ("over_budget", "rf"), ("no_trials", "gbt"),
            ("unequal_planned", None), ("unequal_used", None),
            ("unequal_ok", None)})

    def test_equal_budgets_raise_no_flag(self):
        trials = [trial(0, [1, 2]), trial(0, [1, 2], family="rf")]
        b = sel.budget_accounting(trials, protocol(budgets={"cart": 1, "rf": 1}))
        self.assertEqual(b["flags"], [])

    def test_failures_make_ok_counts_unequal(self):
        trials = [trial(0, [1, 2]), trial(0, [1, 2], family="rf",
                                          status="failed")]
        b = sel.budget_accounting(trials, protocol(budgets={"cart": 1, "rf": 1}))
        self.assertEqual([f["code"] for f in b["flags"]], ["unequal_ok"])


def outer(trials, i, **kw):
    """Select trials as outer fold i: a run with its own fold plan, so its
    own protocol and fold-plan hashes."""
    ph, fh = f"proto-outer{i}", f"plan-outer{i}"
    for t in trials:
        t.update(protocol_hash=ph, fold_plan_hash=fh)
    return sel.select_history(trials, protocol(protocol_hash=ph,
                                               fold_plan_hash=fh, **kw))


class Stability(unittest.TestCase):
    def run_with(self, depth_and_score, run=0):
        """Outer fold `run`: a history where trial 1 (config depth d,
        min_leaf 5, and lr unless it is None) wins clearly."""
        d, lr = depth_and_score
        config = {"depth": d, "min_leaf": 5}
        if lr is not None:
            config["lr"] = lr
        trials = [trial(0, [0, 0], config={"depth": 1, "min_leaf": 5}),
                  trial(1, [9, 9], config=config)]
        return outer(trials, run, budgets={"cart": 2})

    def test_counts_and_spread(self):
        runs = [self.run_with(x, i) for i, x in enumerate(
            [(4, "a"), (4, "a"), (8, "b"), (6, None)])]
        failed = outer([trial(0, [], status="failed"),
                        trial(1, [], status="failed")], 4, budgets={"cart": 2})
        s = sel.stability(runs + [failed], ["o0", "o1", "o2", "o3", "o4"])
        f = s["families"]["cart"]
        self.assertEqual((f["n_runs"], f["selected_in"]), (5, 4))
        self.assertEqual(f["runs_without_selection"], ["o4"])
        self.assertEqual(f["default_selected"], 0)
        self.assertEqual([(c["config"]["depth"], c["count"], c["runs"])
                          for c in f["configs"]],
                         [(4, 2, ["o0", "o1"]), (6, 1, ["o3"]), (8, 1, ["o2"])])
        self.assertEqual(f["configs"][0]["fraction"], 0.5)
        depth = f["params"]["depth"]
        self.assertEqual(depth["present"], 4)
        self.assertEqual(depth["distinct"], 3)
        self.assertEqual(depth["modal_fraction"], 0.5)
        self.assertEqual([(v["value"], v["count"]) for v in depth["value_counts"]],
                         [(4, 2), (6, 1), (8, 1)])
        self.assertEqual(depth["numeric"]["min"], 4)
        self.assertEqual(depth["numeric"]["max"], 8)
        self.assertEqual(depth["numeric"]["range"], 4)
        self.assertEqual(depth["numeric"]["median"], 5.0)
        self.assertEqual(depth["numeric"]["mean"], 5.5)
        self.assertAlmostEqual(depth["numeric"]["sd_population"],
                               math.sqrt(2.75), places=15)
        lr = f["params"]["lr"]
        self.assertEqual(lr["numeric"], {
            "n": 0, "skipped": 3, "min": None, "max": None, "range": None,
            "mean": None, "sd_population": None, "median": None})
        self.assertEqual(lr["present"], 3)
        self.assertEqual(lr["modal_fraction"], 2 / 3)
        self.assertEqual([(v["value"], v["count"], v["fraction"])
                          for v in lr["value_counts"]],
                         [("a", 2, 2 / 3), ("b", 1, 1 / 3)])
        self.assertEqual(f["params"]["min_leaf"]["distinct"], 1)

    def test_integer_and_float_values_stay_distinct(self):
        runs = [outer([trial(0, [0, 0], config={"lr": 0}),
                       trial(1, [9, 9], config={"lr": lr})], i,
                      budgets={"cart": 2}) for i, lr in enumerate((1, 1.0, 1))]
        f = sel.stability(runs)["families"]["cart"]
        self.assertEqual([(c["config"], c["count"]) for c in f["configs"]],
                         [({"lr": 1}, 2), ({"lr": 1.0}, 1)])
        values = [v["value"] for v in f["params"]["lr"]["value_counts"]]
        self.assertEqual([type(v) for v in values], [int, float])
        self.assertEqual(f["params"]["lr"]["distinct"], 2)

    def test_refuses_mixed_objectives_and_duplicate_labels(self):
        a = self.run_with((4, "a"))
        b = outer([trial(0, [1, 2])], 1, direction="minimize",
                  budgets={"cart": 2})
        with self.assertRaises(sel.HistoryError) as cm:
            sel.stability([a, b])
        self.assertEqual(cm.exception.codes, {"mixed_objective"})
        with self.assertRaises(sel.HistoryError) as cm:
            sel.stability([a, self.run_with((4, "a"), 1)], ["x", "x"])
        self.assertEqual(cm.exception.codes, {"bad_labels"})

    def test_numeric_spread_skips_internal_search_lists(self):
        grid = {"internal_search": True, "depth": [1, 2, 3], "kind": ["a", "b"],
                "lr": [0.5]}
        space = b_protocol()["families"]["cart"]["space"]

        def run(i, default_wins):
            p = b_protocol(plan=f"plan-{i}", families={"cart": {
                "default": grid, "space": space}})
            win, lose = [9, 9], [0, 0]
            return sel.select_history([
                b_trial(p, 0, win if default_wins else lose,
                        config=copy.deepcopy(grid)),
                b_trial(p, 1, lose if default_wins else win, depth=4 + i)], p)

        f = sel.stability([run(0, True), run(1, False), run(2, False)])[
            "families"]["cart"]
        self.assertEqual((f["selected_in"], f["default_selected"]), (3, 1))
        depth = f["params"]["depth"]
        self.assertEqual(depth["present"], 3)
        self.assertEqual({k: depth["numeric"][k] for k in
                          ("n", "skipped", "min", "max", "median")},
                         {"n": 2, "skipped": 1, "min": 5, "max": 6,
                          "median": 5.5})
        self.assertEqual([v["value"] for v in depth["value_counts"]],
                         [5, 6, [1, 2, 3]])
        self.assertEqual(f["params"]["lr"]["numeric"]["skipped"], 1)
        flag = f["params"]["internal_search"]
        self.assertEqual((flag["present"], flag["numeric"]["n"],
                          flag["numeric"]["skipped"]), (1, 0, 1))

    def test_refuses_the_same_history_twice(self):
        with self.assertRaises(sel.HistoryError) as cm:
            sel.stability([self.run_with((4, "a")), self.run_with((4, "a"))],
                          ["a", "b"])
        self.assertEqual(cm.exception.codes,
                         {"duplicate_experiment", "duplicate_history"})
        self.assertIn("runs a, b hold the same history",
                      cm.exception.problems[1])

    def test_refuses_one_experiment_pooled_twice(self):
        # The re-check's S8-S11: one outer fold reaches the pool twice
        # through copies that are not byte-identical.
        p = b_protocol()
        h = [b_trial(p, i, [1 + i, 2 + i]) for i in range(4)]
        rerun = [dict(t, seconds=t["seconds"] + 0.001) for t in h]
        for name, other in [("swapped lines", [h[1], h[0], h[2], h[3]]),
                            ("last trial dropped", h[:-1]),
                            ("rerun with measured seconds", rerun),
                            ("split by trial id", h[2:])]:
            first = h[:2] if name == "split by trial id" else h
            with self.assertRaises(sel.HistoryError) as cm:
                sel.stability([sel.select_history(first, p),
                               sel.select_history(other, p)], ["a", "b"])
            self.assertEqual(cm.exception.codes, {"duplicate_experiment"}, name)
            self.assertEqual(cm.exception.problems, [
                f"duplicate_experiment: runs a, b hold the same experiment "
                f"(protocol_hash {p['protocol_hash']})"], name)

    def test_refuses_runs_with_different_protocols(self):
        def run(p, seconds=0.5):
            return sel.select_history([b_trial(p, 0, [1, 2], seconds=seconds),
                                       b_trial(p, 1, [3, 4])], p)

        def codes(*runs):
            with self.assertRaises(sel.HistoryError) as cm:
                sel.stability(list(runs))
            return cm.exception.codes

        base = b_protocol(plan="plan-A")
        # Different fold plans (so different protocol hashes) pool.
        other_plan = b_protocol(plan="plan-B")
        self.assertNotEqual(base["protocol_hash"], other_plan["protocol_hash"])
        s = sel.stability([run(base), run(other_plan)])
        self.assertEqual(s["families"]["cart"]["selected_in"], 2)
        self.assertEqual(codes(run(base), run(b_protocol(plan="plan-B", n=5))),
                         {"mixed_budgets"})
        self.assertEqual(codes(run(base), run(b_protocol(
            [{"metric": "complexity", "op": "<=", "bound": 9}],
            plan="plan-B"))), {"mixed_constraints"})
        self.assertEqual(codes(run(base), run(b_protocol(plan="plan-B",
                                                         fit=8))),
                         {"mixed_protocol"})


if __name__ == "__main__":
    unittest.main()
