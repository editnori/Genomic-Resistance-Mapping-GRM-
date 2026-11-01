"""Clean and planted-leak controls for scripts/hyperopt/canaries.py.

Every canary passes its clean synthetic control and catches its planted leak,
naming the offending samples, groups, folds, arms or tasks. All data is
synthetic and seeded.

usage: python3 -I -m unittest discover -s tests/hyperopt -t .
"""
import hashlib
import enum
import importlib.util
import itertools
import json
import math
import sys
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("hyperopt_canaries",
                                              ROOT / "scripts/hyperopt/canaries.py")
c = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = c
spec.loader.exec_module(c)

K = 4
SEEDS = (0, 1, 2)


def dataset(seed=0, n_groups=24, p=30):
    """60 samples in lineages of sizes 1-4, duplicate pairs inside lineages,
    lineage-preserving outer folds, random binary labels and k-mer features."""
    rng = np.random.default_rng(seed)
    sizes = 1 + np.arange(n_groups) % 4
    lineage_of = np.repeat(np.arange(n_groups), sizes)
    n = len(lineage_of)
    lineage = [f"L{g}" for g in lineage_of]
    duplicate = [None] * n
    for g in range(n_groups):
        members = np.flatnonzero(lineage_of == g)
        if len(members) >= 2:
            duplicate[members[0]] = duplicate[members[1]] = f"D{g}"
    outer = (rng.permutation(n_groups) % K)[lineage_of]
    labels = rng.integers(0, 2, n)
    features = rng.integers(0, 2, (n, p))
    ids = [f"g{i:03d}" for i in range(n)]
    return {"n": n, "lineage": lineage, "duplicate": duplicate, "outer": outer,
            "labels": labels, "features": features, "ids": ids}


class Canary(unittest.TestCase):
    def passes(self, result):
        self.assertTrue(result.passed, f"{result.canary}: {result.offenders} {result.details}")
        json.dumps(result.details)
        return result

    def caught(self, result):
        self.assertFalse(result.passed, f"{result.canary} missed the planted leak")
        json.dumps(result.offenders)
        json.dumps(result.details)
        return result


class OuterFolds(Canary):
    def setUp(self):
        self.d = dataset()
        self.groups = {"lineage": self.d["lineage"], "duplicate": self.d["duplicate"]}

    def test_clean_lineage_folds_pass(self):
        r = self.passes(c.check_groups_within_outer_folds(self.d["outer"], self.groups, self.d["ids"]))
        self.assertEqual(r.details["groups"]["lineage"]["groups"], 24)

    def test_lineage_member_moved_to_another_fold_is_caught(self):
        outer = self.d["outer"].copy()
        member = self.d["lineage"].index("L3")  # L3 has 4 members
        outer[member] = (outer[member] + 1) % K
        r = self.caught(c.check_groups_within_outer_folds(outer, self.groups, self.d["ids"]))
        self.assertEqual([(o["kind"], o["group"]) for o in r.offenders],
                         [("lineage", "L3"), ("duplicate", "D3")])
        moved = [s for f in r.offenders[0]["folds"] for s in f["samples"] if len(f["samples"]) == 1]
        self.assertEqual(moved, [self.d["ids"][member]])

    def test_duplicate_pair_split_inside_one_lineage_is_caught(self):
        # Lineages are not checked here, so only the duplicate group can catch it.
        outer = self.d["outer"].copy()
        a = self.d["duplicate"].index("D1")
        outer[a] = (outer[a] + 1) % K
        r = self.caught(c.check_groups_within_outer_folds(outer, {"duplicate": self.d["duplicate"]}))
        self.assertEqual([o["group"] for o in r.offenders], ["D1"])

    def test_empty_and_singleton_groups_pass(self):
        groups = {"lineage": {"empty": [], "solo": [0], "pair": [1, 2]}}
        r = self.passes(c.check_groups_within_outer_folds([0, 1, 1], groups))
        self.assertEqual(r.details["groups"]["lineage"], {"groups": 3, "empty": 1})

    def test_ungrouped_samples_pass(self):
        self.passes(c.check_groups_within_outer_folds([0, 1, 2], {"duplicate": [None, None, None]}))


def inner_plans(d, split=None, test_sample_in=None):
    """Lineage-preserving inner folds and tuning holdout inside each outer training set."""
    n, outer, lineage = d["n"], d["outer"], d["lineage"]
    plans = {}
    for f in range(K):
        train_lineages = sorted({lineage[i] for i in range(n) if outer[i] != f},
                                key=lambda g: int(g[1:]))
        rank = {g: r for r, g in enumerate(train_lineages)}
        inner = [rank[lineage[i]] % 3 if outer[i] != f else None for i in range(n)]
        holdout = [("holdout" if rank[lineage[i]] % 5 == 0 else "fit") if outer[i] != f else None
                   for i in range(n)]
        plans[f] = {"inner_folds": inner, "holdout": holdout}
    if split is not None:
        # Relabel L3's first member (also in D3) in an outer fold that trains on it.
        member = lineage.index("L3")
        plan = plans[(outer[member] + 1) % K][split]
        plan[member] = ((plan[member] + 1) % 3 if split == "inner_folds"
                        else {"fit": "holdout", "holdout": "fit"}[plan[member]])
    if test_sample_in is not None:
        f = test_sample_in
        plans[f]["holdout"][int(np.flatnonzero(outer == f)[0])] = "fit"
    return plans


class InnerFolds(Canary):
    def setUp(self):
        self.d = dataset()
        self.groups = {"lineage": self.d["lineage"], "duplicate": self.d["duplicate"]}

    def check(self, plans):
        return c.check_groups_within_inner_folds(self.d["outer"], plans, self.groups, self.d["ids"])

    def test_clean_plans_pass(self):
        r = self.passes(self.check(inner_plans(self.d)))
        self.assertEqual(r.details["outer_folds_checked"], K)

    def test_lineage_split_across_inner_folds_is_caught(self):
        r = self.caught(self.check(inner_plans(self.d, split="inner_folds")))
        problems = {(o["partition"], o["problem"], o.get("group")) for o in r.offenders}
        self.assertIn(("inner_folds", "group_crosses", "L3"), problems)
        self.assertIn(("inner_folds", "group_crosses", "D3"), problems)

    def test_lineage_split_across_tuning_holdout_is_caught(self):
        r = self.caught(self.check(inner_plans(self.d, split="holdout")))
        self.assertEqual({o["partition"] for o in r.offenders}, {"holdout"})
        self.assertIn("L3", {o["group"] for o in r.offenders})

    def test_outer_test_sample_in_holdout_is_caught(self):
        r = self.caught(self.check(inner_plans(self.d, test_sample_in=2)))
        [o] = r.offenders
        self.assertEqual((o["outer_fold"], o["problem"]), (2, "test_sample_in_plan"))
        self.assertEqual(o["samples"], [self.d["ids"][int(np.flatnonzero(self.d["outer"] == 2)[0])]])

    def test_unassigned_training_sample_and_unknown_fold_are_caught(self):
        plans = inner_plans(self.d)
        train = int(np.flatnonzero(self.d["outer"] != 0)[0])
        plans[0]["inner_folds"][train] = None
        plans[9] = {}
        r = self.caught(self.check(plans))
        self.assertEqual([o["problem"] for o in r.offenders],
                         ["training_sample_unassigned", "unknown_outer_fold"])

    def test_missing_outer_fold_plans_are_caught(self):
        plans = inner_plans(self.d)
        r = self.caught(self.check({0: plans[0]}))
        self.assertEqual([(o["outer_fold"], o["problem"]) for o in r.offenders],
                         [(f, "missing_outer_fold") for f in (1, 2, 3)])
        r = self.caught(self.check({}))
        self.assertEqual([o["outer_fold"] for o in r.offenders], list(range(K)))
        plans[2] = {}
        r = self.caught(self.check(plans))
        self.assertEqual(r.offenders, [{"outer_fold": 2, "problem": "no_partitions"}])

    def test_empty_and_singleton_groups_pass(self):
        plans = {0: {"inner_folds": [None, 0, 1, 1]}, 1: {"inner_folds": [0, None, None, None]}}
        groups = {"lineage": {"none": [], "solo": [1], "pair": [2, 3]}}
        self.passes(c.check_groups_within_inner_folds([0, 1, 1, 1], plans, groups))


def make_plan_fn(outer, label_leak=False, feature_leak=False):
    """Inner folds stratified by class rank and a config that keeps the k-mer
    columns present in more than half the samples.

    The clean callback ranks and filters on training samples only. The label
    leak ranks samples within the full label array, so test labels shift
    training ranks; the feature leak runs the frequency filter on every row.
    """
    def plan(f, y, x):
        train = np.flatnonzero(outer != f)
        rows = np.arange(len(y)) if label_leak else train
        positive_rank = {int(i): r for r, i in enumerate(rows[y[rows] == 1])}
        negative_rank = {int(i): r for r, i in enumerate(rows[y[rows] == 0])}
        inner = [None] * len(y)
        for i in train:
            inner[i] = (positive_rank if y[i] == 1 else negative_rank)[int(i)] % 3
        score_rows = np.arange(len(y)) if feature_leak else train
        config = np.flatnonzero(x[score_rows].mean(axis=0) > 0.5).tolist()
        return {"inner_folds": inner, "holdout": [i % 5 == 0 for i in train.tolist()],
                "selected_config": config}
    return plan


class LabelInvariance(Canary):
    def setUp(self):
        self.d = dataset()

    def check(self, plan_fn, labels=None, **kw):
        return c.check_label_invariance(self.d["outer"], self.d["labels"] if labels is None else labels,
                                        self.d["features"], plan_fn, **kw)

    def test_training_only_plan_passes(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed)
                r = self.passes(c.check_label_invariance(d["outer"], d["labels"], d["features"],
                                                         make_plan_fn(d["outer"]), seed=seed))
                self.assertGreater(r.details["perturbed_runs"], 0)

    def test_plan_ranked_on_full_label_array_is_caught(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed)
                r = self.caught(c.check_label_invariance(d["outer"], d["labels"], d["features"],
                                                         make_plan_fn(d["outer"], label_leak=True), seed=seed))
                self.assertTrue(all(o["changed"] == ["inner_folds"] for o in r.offenders))
                self.assertEqual({o["perturbed"] for o in r.offenders}, {"labels"})

    def test_callback_that_masks_test_inputs_in_place_passes_and_leaves_caller_arrays(self):
        outer, clean = self.d["outer"], make_plan_fn(self.d["outer"])

        def masking(f, y, x):
            y[outer == f] = -1
            x[outer == f] = 0
            return clean(f, np.where(y < 0, 0, y), x)

        labels, features = self.d["labels"].copy(), self.d["features"].copy()
        self.passes(c.check_label_invariance(outer, labels, features, masking))
        self.passes(c.check_feature_invariance(outer, labels, features, masking))
        np.testing.assert_array_equal(labels, self.d["labels"])
        np.testing.assert_array_equal(features, self.d["features"])

    def test_invert_runs_once_per_fold(self):
        r = self.passes(self.check(make_plan_fn(self.d["outer"]), modes=("invert",), repeats=3))
        self.assertEqual(r.details["perturbed_runs"], K)

    def test_nondeterministic_callback_fails(self):
        rng = np.random.default_rng(1)
        r = self.caught(self.check(lambda f, y, x: {"config": int(rng.integers(1000))}))
        self.assertEqual({p["problem"] for p in r.details["problems"]}, {"callback_not_deterministic"})

    def test_single_class_labels_cannot_be_perturbed_and_fail(self):
        r = self.caught(self.check(make_plan_fn(self.d["outer"]), labels=np.zeros(self.d["n"], int)))
        self.assertEqual(len(r.details["problems"]), K)


class FeatureInvariance(Canary):
    def setUp(self):
        self.d = dataset()

    def check(self, plan_fn):
        return c.check_feature_invariance(self.d["outer"], self.d["labels"], self.d["features"], plan_fn)

    def test_training_only_plan_passes(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed)
                self.passes(c.check_feature_invariance(d["outer"], d["labels"], d["features"],
                                                       make_plan_fn(d["outer"]), seed=seed))

    def test_frequency_filter_on_all_rows_is_caught(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed)
                r = self.caught(c.check_feature_invariance(d["outer"], d["labels"], d["features"],
                                                           make_plan_fn(d["outer"], feature_leak=True),
                                                           seed=seed))
                self.assertTrue(all(o["changed"] == ["selected_config"] for o in r.offenders))

    def test_constant_features_cannot_be_perturbed_and_fail(self):
        r = self.caught(c.check_feature_invariance(self.d["outer"], self.d["labels"],
                                                   np.ones((self.d["n"], 5), int), make_plan_fn(self.d["outer"])))
        self.assertEqual(r.details["problems"],
                         [{"outer_fold": f, "problem": "no perturbation changed the test features"}
                          for f in range(K)])

    def test_missing_features_raise(self):
        with self.assertRaises(ValueError):
            c.check_feature_invariance(self.d["outer"], self.d["labels"], None, lambda f, y, x: {})
        with self.assertRaises(ValueError):
            c.check_feature_screen_training_only(self.d["outer"], self.d["labels"], None, lambda f, y, x: [])

    def test_label_leak_is_not_a_feature_leak(self):
        # Each check isolates one input: the label-leaking plan ignores test features.
        self.passes(self.check(make_plan_fn(self.d["outer"], label_leak=True)))


def make_screen(outer, leak, k=5):
    """Top-k features by |covariance with the label|; the leak screens all rows."""
    def screen(f, y, x):
        rows = np.arange(len(y)) if leak else np.flatnonzero(outer != f)
        xc = x[rows] - x[rows].mean(axis=0)
        score = np.abs(xc.T @ (y[rows] - y[rows].mean()))
        return sorted(np.argsort(-score, kind="stable")[:k].tolist())
    return screen


class FeatureScreen(Canary):
    def setUp(self):
        self.d = dataset(p=40)

    def check(self, screen):
        return c.check_feature_screen_training_only(self.d["outer"], self.d["labels"],
                                                    self.d["features"], screen)

    def test_training_only_screen_passes(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed, p=40)
                r = self.passes(c.check_feature_screen_training_only(
                    d["outer"], d["labels"], d["features"], make_screen(d["outer"], leak=False), seed=seed))
                self.assertGreater(r.details["perturbed_runs"], 0)

    def test_full_data_screen_is_caught(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                d = dataset(seed, p=40)
                r = self.caught(c.check_feature_screen_training_only(
                    d["outer"], d["labels"], d["features"], make_screen(d["outer"], leak=True), seed=seed))
                self.assertEqual({o["perturbed"] for o in r.offenders}, {"labels", "features"})
                self.assertTrue(all(o["added"] and o["removed"] for o in r.offenders))


class FoldPlanHash(Canary):
    def setUp(self):
        d = dataset()
        self.plan = {"scheme": "lineage", "assignments": dict(zip(d["ids"], d["outer"].tolist()))}

    def arms(self, *plans):
        return {f"learner{i}": {"plan": p, "plan_hash": c.fold_plan_hash(p)} for i, p in enumerate(plans)}

    def test_hash_is_sha256_of_canonical_json(self):
        text = c.canonical_json({"b": [1, 2], "a": np.int64(3)})
        self.assertEqual(text, '{"a":3,"b":[1,2]}')
        self.assertEqual(c.fold_plan_hash({"b": [1, 2], "a": 3}),
                         hashlib.sha256(text.encode("utf-8")).hexdigest())

    def test_non_ascii_text_is_escaped_before_hashing(self):
        # Pinned digest: Task B must hash the same plan to the same bytes.
        plan = {"species": "Escherichia coli \u00e9", "marker": "\U0001F9EC"}
        self.assertEqual(c.canonical_json(plan), '{"marker":"\\ud83e\\uddec","species":"Escherichia coli \\u00e9"}')
        self.assertEqual(c.fold_plan_hash(plan), "590f2a1a24e8b460af8d56b3564fc31466a2335241cfed1ee30ebd0e6bddfded")

    def test_key_order_and_numpy_types_do_not_change_the_hash(self):
        shuffled = {"assignments": {k: np.int64(v) for k, v in reversed(self.plan["assignments"].items())},
                    "scheme": "lineage"}
        self.assertEqual(c.fold_plan_hash(shuffled), c.fold_plan_hash(self.plan))

    def test_ambiguous_values_are_rejected(self):
        for bad, error in [({1: 0}, TypeError), ({"a": {1, 2}}, TypeError), ({"a": float("nan")}, ValueError)]:
            with self.assertRaises(error):
                c.canonical_json(bad)

    def test_arms_on_one_plan_pass(self):
        r = self.passes(c.check_common_fold_plan(self.arms(self.plan, self.plan, self.plan)))
        self.assertEqual(len(r.details["hashes"]), 1)

    def test_arm_on_another_plan_is_caught(self):
        other = json.loads(json.dumps(self.plan))
        other["assignments"]["g007"] = (other["assignments"]["g007"] + 1) % K
        r = self.caught(c.check_common_fold_plan(self.arms(self.plan, other, self.plan)))
        [o] = r.offenders
        self.assertEqual((o["arm"], o["problem"], o["reference_arm"]), ("learner1", "different_plan", "learner0"))

    def test_noncanonical_recorded_hash_is_caught(self):
        arms = self.arms(self.plan, self.plan)
        loose = hashlib.sha256(json.dumps(self.plan).encode()).hexdigest()
        arms["learner1"] = {"plan": self.plan, "plan_hash": loose}
        r = self.caught(c.check_common_fold_plan(arms))
        self.assertEqual([(o["arm"], o["problem"]) for o in r.offenders],
                         [("learner1", "noncanonical_hash"), ("learner1", "different_plan")])

    def test_missing_plan_and_single_arm_fail(self):
        h = c.fold_plan_hash(self.plan)
        r = self.caught(c.check_common_fold_plan({"a": {"plan_hash": h}, "b": {"plan_hash": h}}))
        self.assertEqual({o["problem"] for o in r.offenders}, {"plan_missing"})
        self.passes(c.check_common_fold_plan({"a": {"plan_hash": h}, "b": {"plan_hash": h}},
                                             require_plan=False))
        self.caught(c.check_common_fold_plan(self.arms(self.plan)))


class PooledTaskFolds(Canary):
    def setUp(self):
        d = dataset()
        self.fold = dict(zip(d["ids"], d["outer"].tolist()))
        self.lineage = dict(zip(d["ids"], d["lineage"]))
        ids = d["ids"]
        self.tasks = {"ecoli__ampicillin": {g: self.fold[g] for g in ids},
                      "ecoli__ciprofloxacin": {g: self.fold[g] for g in ids[::2]},
                      "ecoli__gentamicin": {g: self.fold[g] for g in ids[::3]}}

    def test_one_global_fold_per_genome_passes(self):
        r = self.passes(c.check_pooled_task_folds(self.tasks, self.lineage))
        self.assertEqual(r.details, {"tasks": 3, "genomes": 60})

    def test_genome_in_test_for_one_task_and_training_for_another_is_caught(self):
        self.tasks["ecoli__ciprofloxacin"]["g006"] = (self.fold["g006"] + 1) % K
        r = self.caught(c.check_pooled_task_folds(self.tasks))
        [o] = r.offenders
        self.assertEqual((o["genome"], o["problem"]), ("g006", "genome_folds_differ"))
        moved = (self.fold["g006"] + 1) % K
        leak = next(x for x in o["leaks"] if x["outer_fold"] == moved)
        self.assertEqual(leak["test_in"], ["ecoli__ciprofloxacin"])
        self.assertEqual(leak["train_in"], ["ecoli__ampicillin", "ecoli__gentamicin"])

    def test_lineage_split_across_tasks_is_caught(self):
        # Each genome keeps one fold, but two lineage mates sit in different folds
        # in different tasks: the pooled model trains on one and tests the other.
        tasks = {"a": {"x1": 0}, "b": {"x2": 1}}
        self.passes(c.check_pooled_task_folds(tasks))
        r = self.caught(c.check_pooled_task_folds(tasks, {"x1": "L", "x2": "L"}))
        self.assertEqual(r.offenders[0]["folds"],
                         [{"outer_fold": 0, "genomes": ["x1"]}, {"outer_fold": 1, "genomes": ["x2"]}])

    def test_empty_task_and_ungrouped_genomes_pass(self):
        self.passes(c.check_pooled_task_folds({"a": {}, "b": {"x": 0}}, {"x": None}))


class OneClassMetrics(Canary):
    def test_metric_values(self):
        y = [0, 0, 1, 1]
        self.assertEqual(c.auroc(y, [0.1, 0.4, 0.35, 0.8]).value, 0.75)
        self.assertEqual(c.auroc(y, [0.5, 0.5, 0.5, 0.5]).value, 0.5)
        self.assertEqual(c.auroc(y, [0.9, 0.8, 0.2, 0.1]).value, 0.0)
        self.assertEqual(c.balanced_accuracy(y, [0.1, 0.6, 0.5, 0.4]).value, 0.5)
        self.assertAlmostEqual(c.log_loss([1, 0], [0.8, 0.4]).value, -(np.log(0.8) + np.log(0.6)) / 2)
        self.assertAlmostEqual(c.brier_score([1, 0], [0.8, 0.4]).value, (0.04 + 0.16) / 2)

    def test_metric_sums_are_correctly_rounded(self):
        # The sixth re-check's inputs, where Python 3.12's compensated sum() lands
        # 1 ulp from the correctly rounded fsum value.
        brier = c.brier_score([0, 1, 0, 0], [0.75, 0.25, 8.195638656616212e-09, 6.622220886733822e-09])
        self.assertEqual(brier.value, 0.28125000000000006)
        p = [1.2664165549094176e-14, 1.2664165549094176e-14, 0.6065306596600346, 0.999999999999999]
        loss = c.log_loss([1, 1, 1, 1], p)
        # Independent of how np.log is dispatched: the result is the correctly
        # rounded sum of whatever terms this platform's np.log gives.
        terms = np.log(np.clip(np.array(p), 1e-15, 1 - 1e-15)).tolist()
        self.assertEqual(loss.value, -math.fsum(terms) / 4)
        # The pin assumes np.log rounds these four inputs correctly, as numpy
        # 1.26.4 does on x86-64 with and without AVX512. If it fails on a new
        # platform while the line above passes, check np.log, not the sum.
        self.assertEqual(loss.value, 16.125000000021682)

    def test_metric_sums_ignore_sample_order(self):
        # sum() gives 0.225 or 0.22500000000000003 depending on order here.
        y = [0, 1, 0, 0, 0]
        p = [0.75, 0.25, 1.0429197416554078e-08, 1.4973404707414789e-09, 1.1006917705286969e-10]
        for metric, expected in ((c.brier_score, 0.22500000000000003), (c.log_loss, c.log_loss(y, p).value)):
            values = {metric([y[i] for i in order], [p[i] for i in order]).value
                      for order in itertools.permutations(range(5))}
            self.assertEqual(values, {expected}, metric.__name__)

    def test_one_class_fold_is_unavailable_not_zero_or_half(self):
        m = c.fold_metrics([1, 1, 1], [0.9, 0.2, 0.6])
        for name in ("auroc", "balanced_accuracy"):
            self.assertIsNone(m[name].value)
            self.assertEqual(m[name].reason, "one class: every label is 1")
        self.assertAlmostEqual(m["log_loss"].value, -np.mean(np.log([0.9, 0.2, 0.6])))
        self.assertAlmostEqual(m["brier"].value, np.mean([0.01, 0.64, 0.16]))

    def test_empty_fold_has_no_metric(self):
        m = c.fold_metrics([], [])
        self.assertEqual({k: (v.value, v.reason) for k, v in m.items()},
                         {k: (None, "no samples") for k in m})

    def test_bad_inputs_raise(self):
        with self.assertRaises(ValueError):
            c.auroc([0, 2], [0.1, 0.2])
        with self.assertRaises(ValueError):
            c.log_loss([0, 1], [0.1, 1.2])
        with self.assertRaises(ValueError):
            c.brier_score([0, 1], [0.1])

    def folds(self):
        rng = np.random.default_rng(3)
        labels = {"two_class": rng.integers(0, 2, 20), "all_resistant": np.ones(6, int),
                  "all_susceptible": np.zeros(5, int), "empty": np.zeros(0, int)}
        return {f: {"labels": y, "metrics": c.fold_metrics(y, rng.random(len(y)))}
                for f, y in labels.items()}

    def test_honest_fold_metrics_pass(self):
        self.passes(c.check_one_class_fold_metrics(self.folds()))

    def test_one_class_fold_scored_as_half_or_zero_is_caught(self):
        folds = self.folds()
        folds["all_resistant"]["metrics"]["auroc"] = 0.5
        folds["all_susceptible"]["metrics"]["balanced_accuracy"] = {"value": 0.0, "reason": None}
        folds["two_class"]["metrics"]["auroc"] = None
        r = self.caught(c.check_one_class_fold_metrics(folds))
        self.assertEqual([(o["fold"], o["metric"], o["problem"]) for o in r.offenders],
                         [("all_resistant", "auroc", "value_on_unavailable_metric"),
                          ("all_susceptible", "balanced_accuracy", "value_on_unavailable_metric"),
                          ("two_class", "auroc", "not_finite")])

    def test_unavailable_without_reason_and_dropped_log_loss_are_caught(self):
        folds = self.folds()
        folds["all_resistant"]["metrics"]["auroc"] = None
        folds["all_resistant"]["metrics"]["log_loss"] = c.Metric(None, "one class")
        del folds["empty"]["metrics"]["brier"]
        r = self.caught(c.check_one_class_fold_metrics(folds))
        self.assertEqual([(o["fold"], o["metric"], o["problem"]) for o in r.offenders],
                         [("all_resistant", "auroc", "unavailable_without_reason"),
                          ("all_resistant", "log_loss", "not_finite"),
                          ("empty", "brier", "missing")])


def sample_permutation(labels, groups, rng, strata=None):
    """The planted null: shuffles samples and ignores groups."""
    return rng.permutation(labels)


def clustered_null(rng, per_size=4):
    """per_size groups of each size 1-5; labels and a feature both set per
    group, independently."""
    n_groups = 5 * per_size
    sizes = np.repeat([1, 2, 3, 4, 5], per_size)
    groups = np.repeat(np.arange(n_groups), sizes)
    while True:
        group_labels = rng.integers(0, 2, n_groups)
        if 0 < group_labels.sum() < n_groups:
            break
    feature = rng.normal(size=n_groups)[groups] + 0.3 * rng.normal(size=len(groups))
    return group_labels[groups], groups, None, feature


mean_difference = c.mean_difference


def dot_mean_difference(labels, feature):
    """mean_difference in numpy floating point: its last bits depend on sample order."""
    positives = labels.sum()
    return feature @ labels / positives - feature @ (1 - labels) / (len(labels) - positives)


class GroupPermutation(Canary):
    def data(self, seed):
        labels, groups, _, feature = clustered_null(np.random.default_rng(seed))
        return labels, groups, feature, np.where(groups < 10, "batch1", "batch2")

    def test_whole_groups_move_and_sizes_and_strata_hold(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, _, strata = self.data(seed)
                out = c.group_permutation(labels, groups, np.random.default_rng(seed), strata)
                for g in range(20):
                    self.assertEqual(len(set(out[groups == g])), 1)
                for s in ("batch1", "batch2"):
                    self.assertEqual(out[strata == s].sum(), labels[strata == s].sum())
                self.passes(c.check_group_level_permutation(labels, out, groups, strata))

    def test_groups_never_cross_strata(self):
        # Every group in stratum A carries 1 and every group in B carries 0, so
        # any cross-stratum swap changes a label; within-stratum swaps never do.
        groups = np.repeat(np.arange(8), 2)
        strata = np.where(groups < 4, "A", "B")
        labels = (strata == "A").astype(int)
        for seed in range(5):
            np.testing.assert_array_equal(
                c.group_permutation(labels, groups, np.random.default_rng(seed), strata), labels)
        self.assertEqual(c.check_permutation_power(labels, groups, strata).details["movable_groups"], 0)
        self.assertEqual(c.check_permutation_power(labels, groups).details["movable_groups"], 8)
        # permutation_null must report the power of the stratified null it draws.
        r = c.permutation_null(lambda y: float(y[:8].sum()), labels, groups, n_perm=19, strata=strata)
        self.assertEqual(r.power.details["movable_groups"], 0)

    def test_seeded_permutation_repeats(self):
        labels, groups, _, _ = self.data(5)
        a = c.group_permutation(labels, groups, np.random.default_rng(9))
        b = c.group_permutation(labels, groups, np.random.default_rng(9))
        np.testing.assert_array_equal(a, b)

    def test_sample_level_shuffle_is_caught(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, _, _ = self.data(seed)
                out = sample_permutation(labels, groups, np.random.default_rng(seed))
                r = self.caught(c.check_group_level_permutation(labels, out, groups))
                self.assertTrue(any(o["group_size"] > 1 for o in r.offenders))

    def test_singleton_ungrouped_and_empty_groups(self):
        labels = np.array([1, 0, 1, 1, 0])
        out = c.group_permutation(labels, [None, None, "a", "a", None], np.random.default_rng(1))
        self.assertEqual(out[2], out[3])
        self.assertEqual(sorted(out.tolist()), sorted(labels.tolist()))
        out = c.group_permutation(labels, {"empty": [], "a": [2, 3]}, np.random.default_rng(1))
        np.testing.assert_array_equal(out[[2, 3]], labels[[2, 3]])

    def test_group_spanning_two_strata_raises(self):
        with self.assertRaises(ValueError):
            c.group_permutation([0, 1], ["a", "a"], np.random.default_rng(0), strata=["x", "y"])

    def test_true_signal_gets_small_p_value(self):
        # 50 groups: 4 per size (self.data) can leave too few labellings, see below.
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, _, feature = clustered_null(np.random.default_rng(seed), per_size=10)
                feature = feature + 10.0 * labels
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                self.passes(r.power)
                self.assertLessEqual(r.p_value, 0.01)
                again = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                           n_perm=199, seed=seed)
                np.testing.assert_array_equal(r.null, again.null)

    def test_ties_with_observed_count_against_it(self):
        # A constant statistic ties every draw: p must be 1, never 1 / (n_perm + 1).
        labels, groups, _, _ = self.data(5)
        r = c.permutation_null(lambda y: 0.0, labels, groups, n_perm=99, seed=0)
        self.assertEqual((r.p_value, r.ties_with_observed), (1.0, 99))
        self.assertEqual([p["problem"] for p in self.caught(r.power).details["problems"]],
                         ["degenerate", "unreachable"])

    def test_all_distinct_group_sizes_are_underpowered(self):
        labels = np.repeat([1, 0, 1, 0, 1], [1, 2, 3, 4, 5])
        groups = np.repeat(np.arange(5), [1, 2, 3, 4, 5])
        r = c.permutation_null(lambda y: mean_difference(y, 5.0 * labels), labels, groups, n_perm=99)
        self.assertEqual(r.p_value, 1.0)
        power = self.caught(r.power)
        self.assertEqual((power.details["movable_groups"], power.details["fixed_groups"],
                          power.details["distinct_labellings"]), (0, 5, 1))
        self.assertEqual([p["problem"] for p in power.details["problems"]],
                         ["underpowered", "underpowered", "underpowered", "degenerate", "unreachable"])

    def test_identical_patterns_do_not_count_as_movable(self):
        # Four pairs that all carry label 1 can swap, but no swap changes a label.
        r = self.caught(c.check_permutation_power([1] * 8 + [0, 1], [0, 0, 1, 1, 2, 2, 3, 3, None, None]))
        self.assertEqual((r.details["movable_groups"], r.details["fixed_groups"],
                          r.details["distinct_labellings"]), (2, 4, 2))
        # clustered_null seed 0: 20 groups, but only 8 can change a label.
        labels, groups, _, _ = self.data(0)
        r = self.caught(c.check_permutation_power(labels, groups, n_perm=199))
        self.assertEqual((r.details["movable_groups"], r.details["fixed_groups"],
                          r.details["distinct_labellings"]), (8, 12, 16))

    def test_reviewer_seed_17_signal_reports_its_power(self):
        # The review's case: clustered_null seed 17 plus a 3.0 shift gave p = 0.055.
        # Its null has 576 labellings over 16 movable groups, so it can reach 0.05,
        # and a large p there is weak evidence, not a degenerate null. At alpha
        # 0.001 the same null cannot reach significance and says so.
        labels, groups, feature, _ = self.data(17)
        feature = feature + 3.0 * labels
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups, n_perm=199)
        power = self.passes(r.power)
        self.assertEqual((power.details["movable_groups"], power.details["fixed_groups"],
                          power.details["distinct_labellings"]), (16, 4, 576))
        self.assertAlmostEqual(power.details["expected_min_p"], (1 + 199 / 576) / 200)
        # 12 of 60 samples sit in fixed groups, all susceptible: 20% exactly, so no coverage note.
        self.assertEqual((power.details["fixed_samples"], power.details["fixed_share"],
                          power.details["fixed_class_counts"]), (12, 0.2, {"0": 12}))
        self.assertNotIn("coverage_note", power.details)
        # permutation_null applies its own max_fixed_share, not the default.
        tight = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups, n_perm=199,
                                   max_fixed_share=0.1)
        self.assertEqual([p["reason"] for p in self.caught(tight.power).details["problems"]],
                         ["fixed groups hold 12 of 60 samples (20%); the limit is 10%"])
        self.assertEqual(r.p_value, (1 + np.sum(r.null >= r.observed)) / 200)
        strict = self.caught(c.check_permutation_power(labels, groups, n_perm=199, alpha=0.001))
        self.assertIn("exceeds alpha 0.001", strict.details["problems"][0]["reason"])


def lineage_layout(rng, paired):
    """The re-check's layout: 20 lineages carry a 3 SD shift by label and 30
    ungrouped genomes carry none. Unique sizes 2-21 leave every lineage alone
    in its size block; paired sizes give each lineage a partner to swap with."""
    sizes = np.repeat(np.arange(2, 12), 2) if paired else np.arange(2, 22)
    group_labels = np.tile([0, 1], 10)
    labels = np.concatenate([np.repeat(group_labels, sizes), rng.integers(0, 2, 30)])
    groups = np.repeat(np.arange(20), sizes).tolist() + [None] * 30
    feature = np.concatenate([np.repeat(3.0 * group_labels + rng.normal(size=20), sizes),
                              rng.normal(size=30)])
    return labels, groups, feature + 0.3 * rng.normal(size=len(feature))


def heavy_tailed_layout(rng, n, exponent=2.3, flip=0.05):
    """The second re-check's clean GRM-like layout: lineage sizes from a discrete
    power law capped at n/8, each lineage resistant with probability 0.3, 5% of
    lineage members flipped, and a lineage-level feature independent of labels."""
    cap = max(2, n // 8)
    sizes = []
    while sum(sizes) < n:
        sizes.append(int(min(rng.zipf(exponent), cap)))
    sizes[-1] -= sum(sizes) - n
    sizes = np.array([s for s in sizes if s > 0])
    rng.shuffle(sizes)
    k = len(sizes)
    labels = np.repeat((rng.random(k) < 0.3).astype(int), sizes)
    lineage = np.repeat(np.arange(k), sizes)
    multi = sizes[lineage] > 1
    labels = np.where((rng.random(n) < flip) & multi, 1 - labels, labels)
    groups = [int(g) if m else None for g, m in zip(lineage, multi)]
    feature = np.repeat(rng.normal(size=k), sizes) + 0.3 * rng.normal(size=n)
    return labels, groups, feature


def clonal_padded_layout(rng):
    """The second re-check's 1,000-genome case: three unique-size lineages (18%)
    carry a 3 SD leak, 300 clonal discordant pairs pad the movable count, and
    220 singletons carry no signal."""
    big, lineage_labels = [90, 60, 30], [1, 1, 0]
    labels = np.concatenate([np.repeat(lineage_labels, big), np.tile([0, 1, 1, 0], 150),
                             rng.integers(0, 2, 220)])
    groups = np.repeat([0, 1, 2], big).tolist() + np.repeat(np.arange(3, 303), 2).tolist() + [None] * 220
    feature = np.concatenate([np.repeat(3.0 * np.array(lineage_labels) + rng.normal(size=3), big),
                              np.repeat(rng.normal(size=300), 2), rng.normal(size=220)])
    return labels, groups, feature


def clonal_near_tie_layout(rng, scatter):
    """The third re-check's layout: 12 singletons, one carrying a perfect leak,
    and 6 clonal discordant pairs. Only the singletons move the statistic, so
    the leak's exact p is 1/12. With scatter, the pairs sit apart in sample
    order and label swaps change the statistic's last bits."""
    labels = np.concatenate([[1], np.zeros(11, int), np.tile([0, 1], 6)])
    feature = np.concatenate([[10.0], rng.normal(size=11), np.repeat(rng.normal(size=6), 2)])
    groups = [None] * 12 + np.repeat(np.arange(6), 2).tolist()
    if scatter:
        order = rng.permutation(len(labels))
        labels, feature, groups = labels[order], feature[order], [groups[i] for i in order]
    return labels, groups, feature


class PermutationPower(Canary):
    def test_frozen_lineages_holding_the_signal_are_underpowered(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, feature = lineage_layout(np.random.default_rng(seed), paired=False)
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                power = self.caught(r.power)
                self.assertEqual((power.details["movable_samples"], power.details["fixed_samples"]),
                                 (30, 230))
                self.assertEqual([p["reason"] for p in power.details["problems"]],
                                 ["fixed groups hold 230 of 260 samples (88%); the limit is 50%"])
                # Lineages of even size are susceptible, odd size resistant.
                self.assertEqual(power.details["fixed_class_counts"], {"0": 110, "1": 120})

    def test_paired_lineages_detect_the_same_signal(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, feature = lineage_layout(np.random.default_rng(seed), paired=True)
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                self.passes(r.power)
                self.assertLessEqual(r.p_value, 0.05)

    def test_clonal_mixed_label_pairs_are_degenerate(self):
        # Duplicate genomes with labels (0, 1) or (1, 0) and one shared feature
        # value: every swap is a new labelling, yet none moves the statistic.
        labels = np.tile([0, 1, 1, 0], 6)
        groups = np.repeat(np.arange(12), 2)
        feature = np.repeat(np.random.default_rng(0).normal(size=12), 2)
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups, n_perm=199)
        power = self.caught(r.power)
        # L counts ordered patterns: (0, 1) and (1, 0) differ, so 924 labellings.
        # Movability counts compositions: every pair holds one 0 and one 1, so none moves.
        self.assertEqual((power.details["movable_groups"], power.details["distinct_labellings"],
                          power.details["fixed_class_counts"]), (0, 924, {"0": 12, "1": 12}))
        self.assertEqual((r.p_value, r.ties_with_observed), (1.0, 199))
        self.assertEqual([p["problem"] for p in power.details["problems"]],
                         ["underpowered", "underpowered", "degenerate", "unreachable"])

    def test_ties_count_equal_draws_only(self):
        # The statistic is 0 when the first three labels agree, as observed, and
        # a continuous score otherwise, so about a fifth of draws tie it and the
        # rest exceed it.
        labels = np.array([1] * 10 + [0] * 10)
        weights = np.random.default_rng(3).random(20)
        r = c.permutation_null(lambda y: 0.0 if y[0] == y[1] == y[2] else float(weights @ y),
                               labels, [None] * 20, n_perm=199)
        self.assertEqual(r.ties_with_observed, int(np.sum(r.null == r.observed)))
        self.assertGreater(np.sum(r.null > r.observed), 0)
        self.assertEqual(r.power.details["ties_with_observed"], r.ties_with_observed)
        self.passes(r.power)

    def test_heavy_tailed_clean_layout_passes(self):
        # Clean GRM-like data fixes 14-43% of samples in lone large lineages.
        # That is coverage, reported as a note, not a failure.
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, feature = heavy_tailed_layout(np.random.default_rng(1000 * 1000 + seed), 1000)
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                power = self.passes(r.power)
                share = power.details["fixed_share"]
                self.assertEqual("coverage_note" in power.details, share > 0.2)
                self.assertEqual(sum(power.details["fixed_class_counts"].values()),
                                 power.details["fixed_samples"])

    def test_fixed_share_limit_edge_and_coverage_note(self):
        # 10 movable singletons plus one lineage alone in its size block.
        at_limit = self.passes(c.check_permutation_power([0, 1] * 5 + [1] * 10, [None] * 10 + ["L"] * 10))
        self.assertEqual(at_limit.details["coverage_note"],
                         "fixed groups hold 10 of 20 samples (50%); the null tests only the other 10")
        over = self.caught(c.check_permutation_power([0, 1] * 5 + [1] * 11, [None] * 10 + ["L"] * 11))
        self.assertEqual([p["reason"] for p in over.details["problems"]],
                         ["fixed groups hold 11 of 21 samples (52%); the limit is 50%"])
        self.assertNotIn("coverage_note", over.details)

    def test_clonal_pairs_padding_a_perfect_leak_fail(self):
        # 12 clonal discordant pairs plus three singletons, one of which leaks
        # perfectly. Only the singletons can move, so p cannot go below 0.35.
        labels = np.concatenate([np.tile([0, 1, 1, 0], 6), [1, 0, 0]])
        groups = np.repeat(np.arange(12), 2).tolist() + [None] * 3
        feature = np.concatenate([np.repeat(np.random.default_rng(0).normal(size=12), 2), [10.0, 0.0, 0.0]])
        for seed in SEEDS:
            with self.subTest(seed=seed):
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                power = self.caught(r.power)
                self.assertEqual(power.details["movable_groups"], 3)
                self.assertGreater(power.details["best_attainable_p"], 0.05)
                self.assertEqual([p["problem"] for p in power.details["problems"]],
                                 ["underpowered", "underpowered", "unreachable"])

    def test_unreachable_alone_fails_a_padded_custom_null(self):
        # The same leak drawn by a custom permute: no group checks apply, and only
        # the best attainable p shows that the draws cannot reach alpha.
        labels = np.concatenate([np.tile([0, 1, 1, 0], 6), [1, 0, 0]])
        groups = np.repeat(np.arange(12), 2).tolist() + [None] * 3
        feature = np.concatenate([np.repeat(np.random.default_rng(0).normal(size=12), 2), [10.0, 0.0, 0.0]])
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups, n_perm=199,
                               permute=lambda y, g, rng, strata: c.group_permutation(y, g, rng, strata))
        power = self.caught(r.power)
        self.assertEqual([p["problem"] for p in power.details["problems"]], ["unreachable"])

    def test_clonal_pairs_cannot_hide_frozen_lineages(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                labels, groups, feature = clonal_padded_layout(np.random.default_rng(500 + seed))
                r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                                       n_perm=199, seed=seed)
                power = self.caught(r.power)
                self.assertEqual(power.details["fixed_samples"], 780)
                self.assertEqual([p["reason"] for p in power.details["problems"]],
                                 ["fixed groups hold 780 of 1000 samples (78%); the limit is 50%"])

    def test_permutation_null_applies_its_movable_floor(self):
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups, n_perm=199,
                               min_movable_groups=60)
        self.assertEqual([p["reason"] for p in self.caught(r.power).details["problems"]],
                         ["40 groups can move; the floor is 60"])

    def test_movable_floor_edge(self):
        self.passes(c.check_permutation_power([0] * 5 + [1] * 5, [None] * 10))
        r = self.caught(c.check_permutation_power([0] * 5 + [1] * 4, [None] * 9))
        self.assertEqual([p["reason"] for p in r.details["problems"]], ["9 groups can move; the floor is 10"])

    def test_permutation_null_honours_alpha(self):
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        feature = feature + 10.0 * labels
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                               n_perm=199, alpha=0.001)
        power = self.caught(r.power)
        self.assertEqual([p["reason"] for p in power.details["problems"]],
                         ["the smallest expected p-value, 0.005, exceeds alpha 0.001",
                          "the best p-value any labelling could get from these draws, 0.005, "
                          "exceeds alpha 0.001"])

    def test_custom_permute_reports_group_power_not_applicable(self):
        # clustered_null seed 0 fails the group power check (8 movable groups),
        # but a sample-level shuffle moves every sample, so that check does not apply.
        labels, groups, feature, _ = GroupPermutation.data(None, 0)
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, groups,
                               n_perm=199, permute=sample_permutation)
        power = self.passes(r.power)
        self.assertEqual(power.details["group_power"], "not applicable")
        self.assertIn("custom permute", power.details["group_power_reason"])
        self.assertNotIn("movable_groups", power.details)

    def test_constant_null_is_degenerate_even_without_ties(self):
        # Only the observed labelling scores 1, and 184,756 labellings make a
        # repeat draw unlikely: p looks tiny, but the null never varies.
        labels = np.array([0, 1] * 10)
        r = c.permutation_null(lambda y: float(np.array_equal(y, labels)), labels, [None] * 20,
                               n_perm=99)
        self.assertEqual((r.p_value, r.ties_with_observed), (0.01, 0))
        power = self.caught(r.power)
        self.assertEqual(power.details["problems"],
                         [{"problem": "degenerate", "reason": "every draw gives the same statistic"}])

    def test_half_the_draws_tying_is_degenerate(self):
        labels = np.array([0, 1] * 10)
        feature = np.arange(20.0)
        swapped = labels[::-1].copy()
        calls = []

        def alternate(y, groups, rng, strata):
            # Every other draw returns the observed labelling.
            calls.append(None)
            return y.copy() if len(calls) % 2 else swapped

        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, None, n_perm=100,
                               permute=alternate)
        self.assertEqual(r.ties_with_observed, 50)
        self.assertEqual([p["problem"] for p in self.caught(r.power).details["problems"]],
                         ["degenerate", "unreachable"])

    def near_ties_match_exact(self, statistic, offset=0.0, scale=None, label=""):
        """The re-checks' near-tie layout, 40 seeds in both sample orders: power
        fails on every seed and p stays at the exact 1/12 level. The scattered
        bound, 0.071, is the re-checks' minimum; adjacent seed 39 draws few ties
        by chance (p 0.067), which exact arithmetic gives too."""
        passed = {False: 0, True: 0}
        for seed in range(40):
            for scatter in (False, True):
                labels, groups, feature = clonal_near_tie_layout(np.random.default_rng(77 + seed), scatter)
                feature = feature + offset
                r = c.permutation_null(lambda y: statistic(y, feature), labels, groups,
                                       n_perm=999, seed=seed, scale=scale)
                with self.subTest(label, seed=seed, scatter=scatter):
                    self.assertGreaterEqual(r.p_value, 0.071 if scatter else 0.05)
                    self.assertEqual(r.power.details["best_attainable_p"], r.p_value)
                    self.assertIn("unreachable", [p["problem"] for p in r.power.details["problems"]])
                passed[scatter] += r.power.passed
        self.assertEqual(passed, {False: 0, True: 0}, label)

    def test_scattered_clonal_pairs_tie_as_adjacent_ones_do(self):
        # Exact float ties let summation-order noise above the observed value
        # lower best_p, and noise below it lower p: scattered pairs passed power
        # 18 of 40 times. An offset leaves the exact statistic unchanged but
        # grows the noise with the summed terms (25 of 40 at 1e3 in round 4);
        # fsum removes the noise at any offset.
        for offset in (0.0, 1e2, 1e3):
            self.near_ties_match_exact(mean_difference, offset, label=f"fsum, offset {offset}")

    def test_floating_point_sums_need_scale(self):
        # numpy's dot product at offset 1e3: the noise outgrows a tolerance taken
        # from the statistic, and only scale, the size of the summed terms, holds.
        passes = 0
        for seed in range(10):
            labels, groups, feature = clonal_near_tie_layout(np.random.default_rng(77 + seed), True)
            r = c.permutation_null(lambda y: dot_mean_difference(y, feature + 1e3), labels, groups,
                                   n_perm=999, seed=seed)
            passes += r.power.passed
        self.assertGreater(passes, 0)
        self.near_ties_match_exact(dot_mean_difference, 1e3, scale=24 * 1e3, label="dot with scale")

    def test_float32_statistic_uses_float32_epsilon(self):
        # A float32 statistic carries about 1e-7 relative noise; float64's
        # epsilon let scattered pairs pass power 20 of 40 times.
        def float32_difference(y, feature):
            f = feature.astype(np.float32)
            return f[y == 1].mean() - f[y == 0].mean()

        self.near_ties_match_exact(float32_difference, label="float32")
        labels, groups, feature = clonal_near_tie_layout(np.random.default_rng(77), True)
        r = c.permutation_null(lambda y: float32_difference(y, feature), labels, groups, n_perm=19)
        self.assertEqual(r.power.details["stat_rtol"], 100 * float(np.finfo(np.float32).eps))

    def test_tolerance_is_100_epsilons_of_the_statistic_dtype(self):
        self.assertEqual(c.STAT_RTOL, 100 * np.finfo(np.float64).eps)
        self.assertEqual(c.STAT_RTOL, 2.220446049250313e-14)
        self.assertEqual(c._stat_rtol(np.float32), 100 * 2.0 ** -23)
        self.assertEqual(c._stat_rtol(np.int64), 0.0)
        # Observed 1: 10 draws 30 ulps above tie it, 10 draws 200 ulps above do not.
        eps = np.finfo(np.float64).eps
        values = iter([1.0] + [1.0 + 30 * eps] * 10 + [1.0 + 200 * eps] * 10 + [0.5] * 79)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation)
        self.assertEqual((r.ties_with_observed, r.p_value), (10, 0.21))
        self.assertEqual(r.power.details["stat_rtol"], c.STAT_RTOL)
        # Integer statistics compare exactly.
        r = c.permutation_null(lambda y: int(y[:5].sum()), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation)
        self.assertEqual(r.power.details["stat_rtol"], 0.0)

    def test_scale_widens_the_tolerance_and_must_be_finite(self):
        # Observed 0.5 with draws 1e-12 above: no tie by the statistic's own
        # size, a tie at scale 1e3 (gamma 2.2e-11).
        for scale, ties in ((None, 0), (1e3, 99)):
            values = iter([0.5] + [0.5 + 1e-12] * 99)
            r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                                   permute=sample_permutation, scale=scale)
            self.assertEqual(r.ties_with_observed, ties, scale)
        # The larger of |observed| and scale counts, not their sum: at observed 1
        # and scale 1, a draw 150 ulps above is not a tie.
        eps = np.finfo(np.float64).eps
        values = iter([1.0] + [1.0 + 150 * eps] * 99)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation, scale=1.0)
        self.assertEqual(r.ties_with_observed, 0)
        for bad in (-1.0, float("nan"), float("inf"), 10 ** 400, True, np.True_, "1", 1j,
                    np.array([1.0]), Decimal("1")):
            with self.subTest(bad=repr(bad)[:20]), self.assertRaisesRegex(ValueError, "scale must be"):
                c.permutation_null(lambda y: 0.0, np.array([0, 1] * 10), [None] * 20, n_perm=9, scale=bad)

    def test_scale_is_reported_as_a_plain_float(self):
        # A float32 scale, as np.abs(f32).sum() returns, or a Fraction would
        # otherwise reach details and break json.dumps.
        for scale in (np.float32(1000), Fraction(1000), 1000, np.int64(1000), 0, None):
            with self.subTest(scale=repr(scale)):
                r = c.permutation_null(lambda y: mean_difference(y, np.arange(20.0)), np.array([0, 1] * 10),
                                       [None] * 20, n_perm=19, scale=scale)
                stat_scale = r.power.details["stat_scale"]
                self.assertIs(type(stat_scale), float)
                self.assertEqual(stat_scale, 0.0 if scale is None else float(scale))
                json.dumps(r.power.details, allow_nan=False)

    def test_scale_sets_the_tie_boundary(self):
        # Observed 0.5 at scale 1e3, so gamma = 2.2e-11: 30 draws 2e-11 below tie
        # it and 30 draws 3e-11 above do not. scale / 10 would drop the 30 ties
        # below observed and lower p from 0.61 to 0.31.
        values = iter([0.5] + [0.5 - 2e-11] * 30 + [0.5 + 3e-11] * 30 + [0.4] * 39)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation, scale=1e3)
        self.assertEqual((r.ties_with_observed, r.p_value), (30, 0.61))

    def test_scale_reaches_the_constant_null_rule(self):
        # Observed 2 with draws 0.5 + {0, 1, 2}e-12: constant at scale 1e3, not
        # by the draws' own size.
        for scale, constant in ((None, False), (1e3, True)):
            values = iter([2.0] + [0.5, 0.5 + 1e-12, 0.5 + 2e-12] * 33)
            r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                                   permute=sample_permutation, scale=scale)
            self.assertEqual({"problem": "degenerate", "reason": "every draw gives the same statistic"}
                             in r.power.details.get("problems", []), constant, scale)

    def test_statistic_types_the_tolerance_cannot_handle_are_refused(self):
        labels = np.array([0, 1] * 10)
        for value in (np.longdouble(0.5), Decimal("0.5"), Fraction(1, 2), 0.5 + 0j, np.complex128(0.5),
                      None, np.array(0.5), np.bool_(True)):
            with self.subTest(type=type(value).__name__):
                with self.assertRaisesRegex(ValueError, f"statistic returned {type(value).__name__}"):
                    c.permutation_null(lambda y: value, labels, [None] * 20, n_perm=9)
        # A refused type among the draws is caught too, and so is one in the
        # observed value alone; the observed value is checked before any draw.
        values = iter([0.5] + [0.25] * 8 + [np.longdouble(0.25)])
        with self.assertRaisesRegex(ValueError, "statistic returned longdouble"):
            c.permutation_null(lambda y: next(values), labels, [None] * 20, n_perm=9)
        for bad, name in ((np.longdouble(0.5), "longdouble"), (Decimal("0.5"), "Decimal"), (None, "NoneType")):
            calls = []

            def observed_only(y, bad=bad):
                calls.append(None)
                return bad if len(calls) == 1 else 0.25

            with self.subTest(observed=name), self.assertRaisesRegex(ValueError, f"statistic returned {name}"):
                c.permutation_null(observed_only, labels, [None] * 20, n_perm=999)
            self.assertEqual(len(calls), 1)
        # An int too large for float64 is refused by name, as observed or as a draw.
        for values in ([10 ** 400] + [1] * 9, [1] * 9 + [-10 ** 400]):
            it = iter(values)
            with self.assertRaisesRegex(ValueError, "int too large for float64"):
                c.permutation_null(lambda y: next(it), labels, [None] * 20, n_perm=9)
        # np.timedelta64 is a numpy integer subclass, but no statistic.
        with self.assertRaisesRegex(ValueError, "statistic returned timedelta64"):
            c.permutation_null(lambda y: np.timedelta64(1, "s"), labels, [None] * 20, n_perm=9)

    def test_accepted_statistic_types_set_their_own_tolerance(self):
        labels = np.array([0, 1] * 10)
        eps = {np.float16: 2.0 ** -10, np.float32: 2.0 ** -23, np.float64: 2.0 ** -52}
        cases = [(lambda y: int(y[:5].sum()), 0.0), (lambda y: bool(y[0]), 0.0),
                 (lambda y: np.int64(y[:5].sum()), 0.0), (lambda y: float(y[:5].sum()), 100 * eps[np.float64])]
        cases += [(lambda y, t=t: t(y[:5].sum()), 100 * e) for t, e in eps.items()]
        for statistic, rtol in cases:
            with self.subTest(rtol=rtol):
                r = c.permutation_null(statistic, labels, [None] * 20, n_perm=9)
                self.assertEqual(r.power.details["stat_rtol"], rtol)
        # Subclasses take their base type's dtype; a float subclass mapped to object
        # dtype got rtol 0, which ignored scale too.
        class SubFloat(float):
            pass

        class SubInt(int):
            pass

        class Code(enum.IntEnum):
            ONE = 1

        class SubFloat32(np.float32):
            pass

        for make, rtol in ((SubFloat, c.STAT_RTOL), (SubInt, 0.0), (lambda v: Code.ONE, 0.0),
                           (SubFloat32, 100 * eps[np.float32])):
            r = c.permutation_null(lambda y: make(int(y[:5].sum())), labels, [None] * 20, n_perm=9)
            self.assertEqual(r.power.details["stat_rtol"], rtol, make)
        values = iter([0.5] + [SubFloat(0.25)] * 9)
        r = c.permutation_null(lambda y: next(values), labels, [None] * 20, n_perm=9)
        self.assertEqual(r.power.details["stat_rtol"], c.STAT_RTOL)
        # The scale applies to a float subclass: draws 1e-12 above tie at scale 1e3.
        values = iter([SubFloat(0.5)] + [SubFloat(0.5 + 1e-12)] * 99)
        r = c.permutation_null(lambda y: next(values), labels, None, n_perm=99,
                               permute=sample_permutation, scale=1e3)
        self.assertEqual(r.ties_with_observed, 99)
        # Mixed types take the result type of all values, observed and draws alike.
        for first, rest in ((float, np.float32), (np.float32, float)):
            values = iter([first(0.5)] + [rest(0.25)] * 9)
            r = c.permutation_null(lambda y: next(values), labels, [None] * 20, n_perm=9)
            self.assertEqual(r.power.details["stat_rtol"], c.STAT_RTOL)

    def test_fsum_mean_difference_ignores_sample_order(self):
        # Swapping the labels of equal-feature samples anywhere in the order
        # leaves every bit; numpy's dot product does not (132 of 200 in round 4).
        rng = np.random.default_rng(0)
        changed = 0
        for _ in range(200):
            feature = rng.normal(size=40) + 1e3
            i, j = rng.choice(40, 2, replace=False)
            feature[j] = feature[i]
            labels = rng.integers(0, 2, 40)
            labels[i], labels[j] = 0, 1
            swapped = labels.copy()
            swapped[i], swapped[j] = 1, 0
            self.assertEqual(mean_difference(labels, feature), mean_difference(swapped, feature))
            changed += dot_mean_difference(labels, feature) != dot_mean_difference(swapped, feature)
        self.assertGreater(changed, 0)
        # Python 3.12's built-in sum is compensated but not exactly rounded: summed
        # in the orders (a, b, c, d) and (a, b, d, c) these values differ in the
        # last bit. Swapping the two copies of c moves c in the class order; fsum
        # gives the same bits in both classes.
        a, b, c_, d = -1.3482420665719877e+32, 1.229321383183538, 1.4949465397505963, -8.143785813409803e+31
        feature = np.array([a, b, c_, d, c_])
        first, second = np.array([1, 1, 1, 1, 0]), np.array([1, 1, 0, 1, 1])
        for labels, swapped in ((first, second), (1 - first, 1 - second)):
            self.assertEqual(mean_difference(labels, feature), mean_difference(swapped, feature))

    def test_scattered_clonal_pairs_alone_are_degenerate(self):
        # Every swap of a scattered clonal pair changes only the statistic's last
        # bits, so every draw ties the observed value. The tolerance is relative,
        # so the statistic must not be 0 in exact arithmetic, as mean_difference
        # would be here; the sum over positives is not.
        rng = np.random.default_rng(5)
        order = rng.permutation(24)
        labels = np.tile([0, 1, 1, 0], 6)[order]
        groups = np.repeat(np.arange(12), 2)[order]
        feature = np.repeat(rng.normal(size=12), 2)[order]
        r = c.permutation_null(lambda y: float(feature @ y), labels, groups, n_perm=199)
        self.assertGreater(np.sum(r.null != r.observed), 0)
        self.assertEqual((r.p_value, r.ties_with_observed), (1.0, 199))
        self.assertEqual([p["problem"] for p in self.caught(r.power).details["problems"]],
                         ["underpowered", "underpowered", "degenerate", "unreachable"])

    def test_tolerance_is_relative_to_the_observed_statistic(self):
        # Observed 1e6: 30 draws 1e-9 above (relative 1e-15, a tie), 30 draws
        # 1e-6 above and 39 draws 1e-6 below (relative 1e-12, not ties).
        values = iter([1e6] + [1e6 + 1e-9] * 30 + [1e6 + 1e-6] * 30 + [1e6 - 1e-6] * 39)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation)
        self.assertEqual((r.ties_with_observed, r.p_value), (30, 0.61))
        # The top draws, 1e6 + 1e-6, are 30 of 99 and do not tie the near-ties below them.
        self.assertEqual(r.power.details["best_attainable_p"], 0.31)
        # best_p takes its tolerance from the top draw, not from the observed value:
        # with observed 1, the 49 draws 1e-9 below a top of 1e6 still tie it.
        values = iter([1.0] + [1e6 + 1e-9] * 50 + [1e6] * 49)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=99,
                               permute=sample_permutation)
        self.assertEqual((r.p_value, r.power.details["best_attainable_p"]), (1.0, 1.0))

    def test_draws_equal_but_for_rounding_are_a_constant_null(self):
        # 0.1 + 0.2 and 0.3 differ by one ulp: the null is constant.
        labels = np.array([0, 1] * 10)
        r = c.permutation_null(lambda y: 1.0 if np.array_equal(y, labels) else (0.1 + 0.2 if y[0] else 0.3),
                               labels, [None] * 20, n_perm=99)
        self.assertGreater(r.power.details["null_variance"], 0.0)
        self.assertEqual(self.caught(r.power).details["problems"],
                         [{"problem": "degenerate", "reason": "every draw gives the same statistic"}])

    def test_best_p_exactly_at_alpha_passes(self):
        # n_perm 19 and an observed value above every draw: p = best_p = 1/20 = alpha.
        labels = np.array([1] * 10 + [0] * 10)
        feature = 10.0 * labels + np.random.default_rng(0).normal(size=20)
        r = c.permutation_null(lambda y: mean_difference(y, feature), labels, None, n_perm=19,
                               permute=sample_permutation)
        self.assertEqual((r.p_value, r.power.details["best_attainable_p"]), (0.05, 0.05))
        self.passes(r.power)

    def test_coverage_note_follows_a_custom_fixed_share_limit(self):
        labels, groups = [0, 1] * 5 + [1] * 15, [None] * 10 + ["L"] * 15
        r = self.passes(c.check_permutation_power(labels, groups, max_fixed_share=0.7))
        self.assertEqual(r.details["coverage_note"],
                         "fixed groups hold 15 of 25 samples (60%); the null tests only the other 10")
        self.assertEqual([p["cause"] for p in self.caught(c.check_permutation_power(labels, groups))
                          .details["problems"]], ["fixed_share"])

    def test_underpowered_problems_name_their_cause(self):
        # 4 movable singletons (floor), one lineage holding 10 of 14 samples
        # (share) and 6 labellings (expected p 0.175).
        r = self.caught(c.check_permutation_power([0, 1] * 2 + [1] * 10, [None] * 4 + ["L"] * 10,
                                                  n_perm=99))
        self.assertEqual([(p["problem"], p["cause"]) for p in r.details["problems"]],
                         [("underpowered", "floor"), ("underpowered", "fixed_share"),
                          ("underpowered", "expected_min_p")])

    def test_nan_observed_statistic_is_neither_leak_nor_clean(self):
        # A NaN observed value compares False with every draw, which used to give
        # the smallest p possible.
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        r = c.permutation_null(lambda y: float("nan") if np.array_equal(y, labels)
                               else mean_difference(y, feature), labels, groups, n_perm=199)
        self.assertTrue(np.isnan(r.p_value))
        power = self.caught(r.power)
        self.assertEqual((power.details["observed_finite"], power.details["non_finite_draws"],
                          power.details["best_attainable_p"]), (False, 0, None))
        self.assertEqual([p["problem"] for p in power.details["problems"]], ["undefined_statistic"])

    def test_nan_draws_are_neither_leak_nor_clean(self):
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        bad = set(np.random.default_rng(1).choice(199, 60, replace=False).tolist())
        calls = []

        def statistic(y):
            calls.append(None)
            return float("nan") if len(calls) - 2 in bad else mean_difference(y, feature)

        r = c.permutation_null(statistic, labels, groups, n_perm=199)
        self.assertTrue(np.isfinite(r.observed))
        self.assertTrue(np.isnan(r.p_value))
        self.assertIsNone(r.ties_with_observed)
        power = self.caught(r.power)
        self.assertEqual((power.details["observed_finite"], power.details["non_finite_draws"]), (True, 60))
        self.assertEqual(power.details["problems"][-1],
                         {"problem": "undefined_statistic",
                          "reason": f"the observed statistic is {r.observed} and 60 of 199 draws "
                                    f"are NaN or infinite"})

    def test_infinite_observed_statistic_is_undefined(self):
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        r = c.permutation_null(lambda y: float("inf") if np.array_equal(y, labels)
                               else mean_difference(y, feature), labels, groups, n_perm=99)
        self.assertTrue(np.isnan(r.p_value))
        self.assertEqual([p["problem"] for p in self.caught(r.power).details["problems"]],
                         ["undefined_statistic"])

    def undefined_draws(self, bad_values):
        """A real leak on clustered_null whose draws k = 0, 1, ... return bad_values[k]."""
        labels, groups, _, feature = clustered_null(np.random.default_rng(0), per_size=10)
        feature = feature + 10.0 * labels
        calls = []

        def statistic(y):
            calls.append(None)
            k = len(calls) - 2
            return bad_values[k] if 0 <= k < len(bad_values) else mean_difference(y, feature)

        return c.permutation_null(statistic, labels, groups, n_perm=199)

    def test_one_nan_draw_is_undefined(self):
        r = self.undefined_draws([float("nan")])
        self.assertTrue(np.isnan(r.p_value))
        power = self.caught(r.power)
        self.assertEqual(power.details["non_finite_draws"], 1)
        self.assertEqual([p["problem"] for p in power.details["problems"]], ["undefined_statistic"])

    def test_infinite_draws_are_undefined(self):
        for bad, count in (([float("inf")], 1), ([float("-inf")] * 60, 60)):
            with self.subTest(bad=bad[0], count=count):
                r = self.undefined_draws(bad)
                self.assertTrue(np.isnan(r.p_value))
                power = self.caught(r.power)
                self.assertEqual(power.details["non_finite_draws"], count)
                self.assertEqual([p["problem"] for p in power.details["problems"]],
                                 ["undefined_statistic"])

    def test_undefined_statistic_keeps_the_group_problems(self):
        r = c.permutation_null(lambda y: float("nan"), [0, 1] * 2 + [1] * 10, [None] * 4 + ["L"] * 10,
                               n_perm=99)
        self.assertEqual([(p["problem"], p.get("cause")) for p in self.caught(r.power).details["problems"]],
                         [("underpowered", "floor"), ("underpowered", "fixed_share"),
                          ("underpowered", "expected_min_p"), ("undefined_statistic", None)])

    def test_overflowing_null_variance_is_reported_as_none(self):
        values = iter([0.0] + [1e308, -1e308] * 50)
        r = c.permutation_null(lambda y: next(values), np.array([0, 1] * 10), None, n_perm=100,
                               permute=sample_permutation)
        self.assertIsNone(r.power.details["null_variance"])
        json.dumps(r.power.details, allow_nan=False)

    def test_constant_null_rule_scales_with_the_draws(self):
        # Observed exactly 0: the draws' own size, not observed, sets the tolerance.
        labels = np.array([0, 1] * 10)
        r = c.permutation_null(lambda y: 0.0 if np.array_equal(y, labels) else (0.1 + 0.2 if y[0] else 0.3),
                               labels, [None] * 20, n_perm=99)
        self.assertIn({"problem": "degenerate", "reason": "every draw gives the same statistic"},
                      self.caught(r.power).details["problems"])

    def test_labels_must_be_zero_or_one(self):
        # NaN != NaN, so two (0, NaN) pairs looked like different compositions.
        groups = np.repeat(np.arange(12), 2)
        good = np.tile([0, 1], 12)
        for bad in (np.nan, 2, 0.5, -1, np.inf, -np.inf):
            labels = np.tile([0.0, bad], 12)
            with self.subTest(bad=bad):
                for call in (lambda: c.check_permutation_power(labels, groups),
                             lambda: c.permutation_null(lambda y: 0.0, labels, groups, n_perm=9),
                             lambda: c.group_permutation(labels, groups, np.random.default_rng(0)),
                             lambda: c.check_group_level_permutation(labels, labels, groups),
                             lambda: c.check_group_level_permutation(labels, good, groups),
                             lambda: c.check_group_level_permutation(good, labels, groups),
                             lambda: c.check_one_class_fold_metrics({"f": {"labels": labels, "metrics": {}}})):
                    with self.assertRaisesRegex(ValueError, "labels must be 0 or 1"):
                        call()
        self.passes(c.check_permutation_power(np.array([True, False] * 10), [None] * 20))

    def test_fold_metric_labels_must_be_zero_or_one(self):
        # astype(int) turned [1, 1, NaN] into two classes, [0, 0.6] into one and
        # counted [1, 2] as two.
        metrics = {"auroc": 0.7, "balanced_accuracy": 0.6, "log_loss": 0.5, "brier": 0.2}
        for labels in ([1, 1, np.nan], [0, 0.6], [1, 2]):
            with self.subTest(labels=labels), self.assertRaisesRegex(ValueError, "labels must be 0 or 1"):
                c.check_one_class_fold_metrics({"f": {"labels": labels, "metrics": metrics}})
        self.passes(c.check_one_class_fold_metrics({"f": {"labels": [0.0, 1.0], "metrics": metrics}}))

    def test_permuted_labels_need_one_entry_per_sample(self):
        labels, groups = np.tile([0, 0, 1, 1], 6), np.repeat(np.arange(12), 2)
        for permuted in (np.append(labels, [0, 1]), labels[:-2], labels[:-1], np.append(labels, 0)):
            with self.subTest(n=len(permuted)):
                with self.assertRaisesRegex(ValueError, f"permuted has {len(permuted)} entries for 24 labels"):
                    c.check_group_level_permutation(labels, permuted, groups)

    def test_custom_permute_draws_are_checked(self):
        # A statistic that slices or pads would score a malformed draw quietly.
        labels, feature = np.array([0, 1] * 12), np.arange(26.0)
        bad_draws = {"short": lambda y: y[:-2], "long": lambda y: np.append(y, [0, 1]),
                     "doubled": lambda y: 2 * y, "nan": lambda y: np.where(np.arange(24) == 3, np.nan, y),
                     "half": lambda y: y * 0.5}
        for name, broken in bad_draws.items():
            with self.subTest(name=name):
                message = "permute returned" if name in ("short", "long") else "labels must be 0 or 1"
                with self.assertRaisesRegex(ValueError, message):
                    c.permutation_null(lambda y: float(y[:10].sum() + feature[:len(y)] @ np.nan_to_num(y)),
                                       labels, None, n_perm=9,
                                       permute=lambda y, g, rng, s, b=broken: b(rng.permutation(y)))
                with self.assertRaisesRegex(ValueError, message):
                    c.null_calibration(lambda rng: (labels, None, None, None), lambda y, data: float(y[:10].sum()),
                                       n_sim=3, permute=lambda y, g, rng, s, b=broken: b(rng.permutation(y)))

    def test_custom_permute_draws_are_checked_on_every_draw(self):
        # Permuters that are right at first and break on a later draw: a check of
        # the first draw only would let the rest reach the statistic.
        labels, feature = np.array([0, 1] * 12), 100.0 + np.array([0, 1] * 12)

        def breaks_on(k, broken, message):
            calls = []

            def permute(y, g, rng, s):
                calls.append(None)
                drawn = rng.permutation(y)
                return broken(drawn) if len(calls) == k else drawn

            with self.assertRaisesRegex(ValueError, message):
                c.permutation_null(lambda y: float(np.asarray(y)[:24] @ feature), labels, None,
                                   n_perm=99, permute=permute)
            self.assertEqual(len(calls), k)

        breaks_on(2, lambda y: 2 * y, "labels must be 0 or 1")
        breaks_on(50, lambda y: y[:-2], "permute returned 22 labels for 24 samples")
        breaks_on(99, lambda y: np.where(np.arange(24) == 0, 1 - y, y), "permute changed the class counts")

    def test_custom_permute_draws_must_be_permutations(self):
        # A resampler keeps 0/1 values and length but not the class counts; with a
        # sum statistic a perfect leak then read as clean with power passing.
        labels = np.array([0, 1] * 12)
        feature = 100.0 + labels

        def leak(y):
            return float(y @ feature)

        refused = {
            "with replacement": (lambda y, g, rng, s: rng.choice(y, size=len(y)), "changed the class counts"),
            # Off by one in each direction: one more positive, one fewer.
            "one more positive": (lambda y, g, rng, s: np.where(np.arange(24) == 0, 1, y),
                                  r"changed the class counts from \[12, 12\] to \[11, 13\]"),
            "one fewer positive": (lambda y, g, rng, s: np.where(np.arange(24) == 1, 0, y),
                                   r"changed the class counts from \[12, 12\] to \[13, 11\]"),
            "two columns": (lambda y, g, rng, s: np.column_stack([rng.permutation(y)] * 2),
                            r"permute returned shape \(24, 2\) for labels of shape \(24,\)"),
        }
        for name, (permute, message) in refused.items():
            with self.subTest(name), self.assertRaisesRegex(ValueError, message):
                c.permutation_null(leak, labels, None, n_perm=199, permute=permute)
        # Unbalanced labels pin which class holds which count: flipping 8 of 24
        # positives to 16 keeps the multiset {8, 16} of counts (MD-1).
        unbalanced = np.array([1] * 8 + [0] * 16)
        with self.subTest("flipped classes"), self.assertRaisesRegex(
                ValueError, r"changed the class counts from \[16, 8\] to \[8, 16\]"):
            c.permutation_null(lambda y: float(y @ (100.0 + unbalanced)), unbalanced, None,
                               n_perm=199, permute=lambda y, g, rng, s: 1 - rng.permutation(y))
        # The shape is checked on every draw, not only the first (MD-3): a
        # second draw shaped (24, 1) holds the right labels and counts.
        calls = []

        def second_draw_is_a_column(y, g, rng, s):
            calls.append(None)
            drawn = rng.permutation(y)
            return drawn if len(calls) == 1 else drawn.reshape(-1, 1)

        with self.subTest("later draw reshaped"), self.assertRaisesRegex(
                ValueError, r"permute returned shape \(24, 1\) for labels of shape \(24,\)"):
            c.permutation_null(leak, labels, None, n_perm=199, permute=second_draw_is_a_column)
        self.assertEqual(len(calls), 2)
        # A 0-d draw is refused loudly, as len()'s TypeError rather than a
        # ValueError (MD-5).
        with self.subTest("0-d draw"), self.assertRaisesRegex(TypeError, "unsized"):
            c.permutation_null(leak, labels, None, n_perm=199,
                               permute=lambda y, g, rng, s: np.array(1))
        kept = {
            "true permutation": lambda y, g, rng, s: rng.permutation(y),
            "stratified permutation": lambda y, g, rng, s: np.concatenate([rng.permutation(y[:12]),
                                                                           rng.permutation(y[12:])]),
            "list": lambda y, g, rng, s: rng.permutation(y).tolist(),
        }
        for name, permute in kept.items():
            with self.subTest(name):
                r = c.permutation_null(leak, labels, None, n_perm=199, permute=permute)
                self.assertEqual(r.p_value, 0.005)
                self.passes(r.power)
        r = c.permutation_null(leak, labels.astype(bool), None, n_perm=19,
                               permute=lambda y, g, rng, s: rng.permutation(y))
        self.assertEqual(r.p_value, 0.05)

    def test_strata_need_one_entry_per_sample(self):
        labels = [0, 1] * 10
        calls = []

        def custom(y, groups, rng, strata):
            calls.append(None)
            return rng.permutation(y)

        for n in (15, 19, 21, 25, 0):
            strata = ["a"] * (n // 2) + ["b"] * (n - n // 2)
            message = f"strata has {n} entries for 20 samples"
            with self.subTest(n=n):
                with self.assertRaisesRegex(ValueError, message):
                    c.check_permutation_power(labels, [None] * 20, strata=strata)
                with self.assertRaisesRegex(ValueError, message):
                    c.permutation_null(lambda y: 0.0, labels, [None] * 20, n_perm=9, strata=strata)
                # A custom permute never reaches _blocks; the length is checked first.
                with self.assertRaisesRegex(ValueError, message):
                    c.permutation_null(lambda y: 0.0, labels, [None] * 20, n_perm=9, strata=strata,
                                       permute=custom)
                with self.assertRaisesRegex(ValueError, message):
                    c.null_calibration(lambda rng: (np.array(labels), [None] * 20, strata, None),
                                       lambda y, data: 0.0, n_sim=3, permute=custom)
        self.assertEqual(calls, [])


def coarse_difference(labels, feature):
    """A valid statistic with heavy ties, like balanced accuracy on a small fold."""
    return float(np.round(mean_difference(labels, feature)))


class NullCalibration(Canary):
    def test_uniformity_decides_on_small_p_values_only(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                u = np.random.default_rng(seed).random(200)
                self.passes(c.check_pvalue_uniformity(u))
                r = self.caught(c.check_pvalue_uniformity(u ** 2))
                self.assertGreater(r.details["d_plus"], r.details["critical"])
                # p-values piled toward 1 are conservative: noted, not failed.
                r = self.passes(c.check_pvalue_uniformity(np.sqrt(u)))
                self.assertTrue(r.details["conservative"])
                self.assertGreater(r.details["d_minus"], r.details["critical"])
        self.caught(c.check_pvalue_uniformity([]))

    def test_one_sided_critical_value_decides(self):
        # D+ = 0.11 exactly: above the one-sided cutoff at alpha 0.01, n = 200
        # (0.1073), below the two-sided one (0.1151). It must fail.
        p = np.clip(np.arange(1, 201) / 200 - 0.11, 0, 1)
        r = self.caught(c.check_pvalue_uniformity(p, alpha=0.01))
        self.assertAlmostEqual(r.details["d_plus"], 0.11)

    def test_group_permutation_null_is_calibrated(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                r = self.passes(c.null_calibration(clustered_null, mean_difference, seed=seed))
                self.assertLess(r.details["share_below_0.05"], 0.1)

    def test_tied_statistic_null_is_conservative_and_passes(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                r = self.passes(c.null_calibration(clustered_null, coarse_difference, seed=seed))
                self.assertTrue(r.details["conservative"])

    def test_sample_permutation_null_on_clustered_data_is_caught(self):
        for seed in SEEDS:
            with self.subTest(seed=seed):
                r = self.caught(c.null_calibration(clustered_null, mean_difference, seed=seed,
                                                   permute=sample_permutation))
                self.assertGreater(r.details["share_below_0.05"], 0.1)

    def test_undefined_statistic_fails_calibration(self):
        r = self.caught(c.check_pvalue_uniformity([0.2, float("nan"), 0.7]))
        self.assertEqual(r.details["problems"],
                         [{"problem": "undefined_statistic", "reason": "1 of 3 p-values are NaN or infinite"}])
        r = self.caught(c.null_calibration(clustered_null, lambda y, data: float("nan"), n_sim=5))
        self.assertEqual([p["problem"] for p in r.details["problems"]], ["undefined_statistic"])

    def test_calibration_refuses_too_few_simulations(self):
        # The data equal the labels, so p-values sit at or near the floor 0.01.
        # At n_sim 2 the cutoff, 1.073, exceeds the largest possible D+, 0.99.
        def leak(rng):
            labels, groups, strata, _ = clustered_null(rng)
            return labels, groups, strata, labels.copy()

        def matches(y, data):
            return int(np.sum(y == data))

        self.assertEqual(c.min_calibration_sims(99, 0.01), 3)
        self.assertEqual(c.min_calibration_sims(19, 0.01), 3)
        self.assertEqual(c.min_calibration_sims(99, 0.001), 4)
        # The floor matters: at n_perm 19, alpha 0.003 needs n_sim > 2.90 / 0.95**2 = 3.22.
        self.assertEqual(c.min_calibration_sims(19, 0.003), 4)
        for n_sim in (-1, 0, 1, 2):
            with self.subTest(n_sim=n_sim), self.assertRaisesRegex(ValueError, "needs n_sim >= 3"):
                c.null_calibration(leak, matches, n_sim=n_sim)
        r = self.caught(c.null_calibration(leak, matches, n_sim=3))
        self.assertAlmostEqual(r.details["critical"], 0.8761, places=4)
        self.assertGreater(r.details["d_plus"], r.details["critical"])
        self.assertEqual(r.offenders[0]["problem"], "pvalues_anti_conservative")
        # Off the defaults null_calibration passes its own n_perm and alpha to the
        # minimum. At alpha 0.0025 (0.05 / 20) n_sim 3 has cutoff 0.9993 < 1, so
        # without the refusal the leak would pass with no problem attached.
        for n_perm, alpha in ((99, 0.0025), (19, 0.003)):
            with self.subTest(n_perm=n_perm, alpha=alpha):
                with self.assertRaisesRegex(ValueError, f"needs n_sim >= 4 at n_perm {n_perm} and alpha {alpha}"):
                    c.null_calibration(leak, matches, n_sim=3, n_perm=n_perm, alpha=alpha)
                r = self.caught(c.null_calibration(leak, matches, n_sim=4, n_perm=n_perm, alpha=alpha))
                self.assertEqual(r.offenders[0]["problem"], "pvalues_anti_conservative")

    def test_calibration_arguments_are_checked(self):
        # Each is refused up front, by name, not by a later ZeroDivisionError,
        # math domain error or range() TypeError.
        for kwargs, message in (({"n_perm": 0}, "n_perm >= 1"), ({"n_perm": 2.5}, "n_perm >= 1"),
                                ({"n_perm": True}, "n_perm >= 1"),
                                ({"alpha": 0.0}, "0 < alpha < 1"), ({"alpha": 1.0}, "0 < alpha < 1"),
                                ({"alpha": 1.5}, "0 < alpha < 1"),
                                ({"n_sim": 3.5}, "n_sim must be an integer"),
                                ({"n_sim": True}, "n_sim must be an integer")):
            with self.subTest(**kwargs), self.assertRaisesRegex(ValueError, message):
                c.null_calibration(clustered_null, mean_difference, **kwargs)
        # numpy integers are integers.
        r = c.null_calibration(clustered_null, mean_difference, n_sim=np.int64(5), n_perm=np.int64(19))
        self.assertEqual((r.details["n"], r.details["n_perm"]), (5, 19))

    def test_uniformity_that_could_never_fail_is_underpowered(self):
        for p in ([0.0], [0.0, 0.0]):
            r = self.caught(c.check_pvalue_uniformity(p))
            self.assertEqual(r.offenders, [])
            self.assertEqual([(q["problem"], q["cause"]) for q in r.details["problems"]],
                             [("underpowered", "too_few_pvalues")])
        r = self.caught(c.check_pvalue_uniformity([0.0, 0.0, 0.0]))
        self.assertNotIn("problems", r.details)
        # A cutoff of exactly 1 (one p-value at alpha e**-2) still cannot fail.
        r = self.caught(c.check_pvalue_uniformity([0.0], alpha=math.exp(-2)))
        self.assertEqual(r.details["critical"], 1.0)
        self.assertEqual(r.details["problems"][0]["cause"], "too_few_pvalues")

    def test_infinite_pvalue_is_undefined(self):
        r = self.caught(c.check_pvalue_uniformity([0.2, float("inf"), 0.7]))
        self.assertEqual([p["problem"] for p in r.details["problems"]], ["undefined_statistic"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
