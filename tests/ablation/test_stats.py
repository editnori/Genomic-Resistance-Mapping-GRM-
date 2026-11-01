"""Tests for scripts/ablation/stats.py, loaded by path so sibling modules are not needed.

KSchemaTests runs J's finalize() output through Task K's schema.py: this
repo's copy once K is merged, else the file named by GRM_TASK_K_SCHEMA.
"""
import importlib.util
import itertools
import json
import os
import sys
import types
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("ablation_stats", ROOT / "scripts/ablation/stats.py")
stats = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = stats
SPEC.loader.exec_module(stats)

MARGINS = dict.fromkeys(stats.ENDPOINTS, 0.02)
# Student t quantiles from published tables, independent of stats.t_quantile.
T975 = {1: 12.7062047362, 2: 4.3026527297, 4: 2.7764451052, 9: 2.2621571628, 19: 2.0930240544, 39: 2.0226909200}
T95 = {1: 6.3137515147, 2: 2.9199855804, 4: 2.1318467863, 9: 1.8331129327, 19: 1.7291328115, 39: 1.6848751217}


def sigmoid(x):
    return 1 / (1 + np.exp(-x))


def make_rows(n_clusters=6, per_cluster=4, repeats=2, folds=3, seed=0, dataset="ds"):
    """Two arms on one plan; every cluster holds both classes and sits in one fold."""
    rng = np.random.default_rng(seed)
    rows = [{"dataset": dataset, "arm": arm, "learner": "lr", "split": "lineage", "repeat": repeat,
             "sample_id": f"s{c}-{j}", "y": j % 2, "fold": (c + repeat) % folds,
             "score": float(rng.uniform(0.05, 0.95)), "fold_plan_sha256": f"plan{repeat}", "cluster": f"c{c}"}
            for arm in ("base", "treat") for repeat in range(repeats)
            for c in range(n_clusters) for j in range(per_cluster)]
    for row in rows:
        row["called"] = row["score"] >= 0.5
    return rows


def pair(rows, dataset="ds"):
    return stats.pair_arms(rows, dataset=dataset, learner="lr", split="lineage", baseline="base", treatment="treat")


def arrays(y, cluster, base, treat, null=None, dataset="sim"):
    """PairedArms built directly, for simulations that skip row validation."""
    y = np.asarray(y)
    clusters, index = np.unique(np.asarray(cluster), return_inverse=True)
    return stats.PairedArms(dataset=dataset, learner="lr", split="lineage", baseline="base", treatment="treat",
                            sample_ids=tuple(range(len(y))), y=y, cluster_index=index, clusters=tuple(clusters),
                            base=np.asarray(base, float), treat=np.asarray(treat, float),
                            null=np.full(len(y), y.mean()) if null is None else np.asarray(null, float),
                            n_repeats=1)


def clustered(rng, k, m=8, base_slope=0.3, treat_slope=0.9):
    """Labels and two arms' scores with shared cluster effects."""
    z = np.repeat(rng.normal(0, 1, k), m)
    eta = z + rng.normal(0, 1, k * m)
    y = (rng.uniform(size=k * m) < sigmoid(eta)).astype(int)
    noise = lambda: np.repeat(rng.normal(0, 0.5, k), m) + rng.normal(0, 0.5, k * m)
    return y, np.repeat(np.arange(k), m), sigmoid(base_slope * eta + noise()), sigmoid(treat_slope * eta + noise())


def balanced(k, m=4, gap=0.0, seed=0, dataset="sim"):
    """k clusters of m samples, both classes in each; the treatment moves scores towards y by gap."""
    rng = np.random.default_rng(seed)
    y = np.tile(np.arange(m) % 2, k)
    base = np.clip(0.5 + 0.2 * (y - 0.5) + rng.normal(0, 0.1, k * m), 0.01, 0.99)
    return arrays(y, np.repeat(np.arange(k), m), base, np.clip(base + gap * (y - 0.5), 0.01, 0.99), dataset=dataset)


def ks_uniform(p):
    p = np.sort(p)
    n = len(p)
    return max(np.max(np.arange(1, n + 1) / n - p), np.max(p - np.arange(n) / n))


def fake_contrast(family, dataset, p, *, endpoint="log_loss", estimate=0.02, ci=(0.01, 0.03), tost=(0.012, 0.028),
                  n_units=60, margin=0.05):
    r = dict.fromkeys(stats.RECORD_KEYS)
    defined = estimate is not None
    r.update(id=f"{dataset}/{family}:lr:lineage:base->treat:{endpoint}", family=family, endpoint=endpoint,
             estimate=estimate, ci=list(ci) if defined else None, ci_tost=list(tost) if defined else None, p=p,
             ci_method="cluster_jackknife_t" if defined else "none",
             adjust=stats.FAMILIES[family], margin=margin, unit="cluster", n_units=n_units, dataset=dataset,
             _pending=True)
    return r


def fake_pooled(family, treatment, p, *, estimate=0.02, n_units=12, margin=0.05):
    r = dict.fromkeys(stats.POOLED_KEYS)
    defined = estimate is not None
    r.update(contrast=f"{family}:lr:lineage:base->{treatment}:r2_ll", family=family, estimate=estimate,
             ci=[0.01, 0.03] if defined else None, ci_tost=[0.012, 0.028] if defined else None, p=p,
             adjust=stats.FAMILIES[family], margin=margin, unit="species", n_units=n_units, _pending=True)
    return r


def ids(records):
    return [r.get("contrast") or r["id"] for r in records]


class PairingTests(unittest.TestCase):
    def test_scores_average_repeats_and_null_is_training_prevalence(self):
        layout = {0: {"a": 0, "b": 0, "c": 1, "d": 1}, 1: {"a": 0, "c": 0, "b": 1, "d": 1}}
        y = {"a": 1, "b": 0, "c": 1, "d": 1}
        rows = [{"dataset": "ds", "arm": arm, "learner": "lr", "split": "lineage", "repeat": r, "sample_id": s,
                 "y": y[s], "fold": layout[r][s], "score": 0.2 + 0.4 * r + (0.1 if arm == "treat" else 0.0),
                 "fold_plan_sha256": f"plan{r}", "cluster": s} for arm in ("base", "treat") for r in (0, 1) for s in y]
        for row in rows:
            row["called"] = row["score"] >= 0.5
        paired = pair(rows)
        self.assertEqual(paired.sample_ids, ("a", "b", "c", "d"))
        np.testing.assert_allclose(paired.base, 0.4)
        np.testing.assert_allclose(paired.treat, 0.5)
        # Repeat 0: fold 0 trains on c, d (1.0), fold 1 on a, b (0.5). Repeat 1: fold 0 on b, d (0.5), fold 1 on a, c (1.0).
        np.testing.assert_allclose(paired.null, [0.75, 1.0, 0.5, 0.75])
        self.assertEqual(paired.n_repeats, 2)

    def test_each_mismatch_is_refused(self):
        def edit(rows, check):
            target = next(r for r in rows if r["arm"] == "treat" and r["repeat"] == 1 and r["sample_id"] == "s2-1")
            if check == "fold_plan":
                for r in rows:
                    if r["arm"] == "treat" and r["repeat"] == 1:
                        r["fold_plan_sha256"] = "other"
            elif check == "samples":
                for r in rows:
                    if r["arm"] == "treat" and r["sample_id"] == "s2-1":
                        r["sample_id"] = "s9-9"
            elif check == "y":
                target["y"] = 1 - target["y"]
            elif check == "cluster":
                target["cluster"] = "c5"
            elif check == "fold":
                target["fold"] = (target["fold"] + 1) % 3
            elif check == "duplicate":
                rows.append(dict(target))
            elif check == "called":
                target["called"] = not target["called"]
            return rows

        pair(make_rows())
        for check in ("fold_plan", "samples", "y", "cluster", "fold", "duplicate", "called"):
            with self.subTest(check=check), self.assertRaises(stats.PairingError) as caught:
                pair(edit(make_rows(), check))
            self.assertEqual(caught.exception.check, check)

    def test_strict_labels_scores_and_called(self):
        def first(rows):
            return next(r for r in rows if r["arm"] == "base" and r["repeat"] == 0 and r["sample_id"] == "s0-1")

        for field, value, check in (("y", True, "y"), ("y", 1.0, "y"), ("score", True, "score"),
                                    ("score", float("nan"), "score"), ("score", 1.5, "score"),
                                    ("called", 1, "called"), ("called", None, "called")):
            rows = make_rows()
            first(rows)[field] = value
            with self.subTest(field=field, value=value), self.assertRaises(stats.PairingError) as caught:
                pair(rows)
            self.assertEqual(caught.exception.check, check)
        rows = make_rows()
        del first(rows)["called"]
        with self.assertRaises(ValueError):
            pair(rows)

    def test_no_call_row_scores_its_training_prevalence(self):
        # Every training set in make_rows has prevalence 0.5, which is not a call on a no-call row.
        rows = make_rows()
        row = next(r for r in rows if r["arm"] == "base" and r["repeat"] == 0 and r["sample_id"] == "s0-1")
        row.update(score=0.5, called=False)
        pair(rows)
        row.update(score=0.6, called=False)
        with self.assertRaises(stats.PairingError) as caught:
            pair(rows)
        self.assertEqual(caught.exception.check, "called")
        row.update(score=0.3, called=True)
        with self.assertRaises(stats.PairingError):
            pair(rows)

    def test_thresholded_endpoints_use_called(self):
        # Both arms score alike; resistant s0-1 is a no-call (score = prevalence 0.5) in both baseline repeats.
        rows = make_rows()
        scores = {(r["repeat"], r["sample_id"]): (r["score"], r["called"]) for r in rows if r["arm"] == "base"}
        for r in rows:
            r["score"], r["called"] = scores[(r["repeat"], r["sample_id"])]
            if r["sample_id"] == "s0-1":
                r.update(score=0.5, called=r["arm"] == "treat")
            if r["sample_id"] == "s1-0" and r["arm"] == "base":  # susceptible, called in one of two repeats
                r.update(score=0.5, called=False) if r["repeat"] == 0 else r.update(score=0.7, called=True)
        paired = pair(rows)
        i, j = paired.sample_ids.index("s0-1"), paired.sample_ids.index("s1-0")
        self.assertEqual((paired.base[i], bool(paired.base_call[i]), bool(paired.treat_call[i])), (0.5, False, True))
        self.assertTrue(paired.base_call[j])  # half the repeats is a call
        n_r = int(paired.y.sum())
        n_s = len(paired.y) - n_r
        vme, me = stats.compare(paired, family="D", endpoints=["vme", "me"], margins=MARGINS, development=())
        # Only s0-1 differs among resistant samples: a false S in the baseline alone.
        self.assertAlmostEqual(vme["estimate"], 1 / n_r)
        # Only s1-0 differs among susceptible samples: a false R in the baseline, and in the treatment if called there.
        self.assertAlmostEqual(me["estimate"], (1 - int(paired.treat_call[j])) / n_s)


class EndpointTests(unittest.TestCase):
    def test_rank_metrics_match_pairwise_definitions_with_ties(self):
        y = np.array([1, 0, 1, 1, 0, 0, 1, 0])
        s = np.array([0.9, 0.9, 0.7, 0.4, 0.4, 0.2, 0.2, 0.1])
        pos, neg = s[y == 1], s[y == 0]
        pairwise = np.mean([1.0 if p > q else 0.5 if p == q else 0.0 for p in pos for q in neg])
        self.assertAlmostEqual(stats.auroc(y, s), pairwise)
        ap, before = 0.0, 0.0
        for cut in sorted(set(s), reverse=True):
            chosen = y[s >= cut]
            recall = chosen.sum() / y.sum()
            ap += (recall - before) * chosen.mean()
            before = recall
        self.assertAlmostEqual(stats.auprc(y, s), ap)

    def test_rank_deltas_are_positive_when_the_treatment_ranks_better(self):
        y = np.array([1, 0, 1, 0, 1, 0, 1, 0])
        base = np.array([0.3, 0.6, 0.5, 0.4, 0.2, 0.7, 0.8, 0.1])
        treat = np.where(y == 1, 0.9, 0.1)
        null = np.full(8, 0.5)
        for endpoint, metric in (("auroc", stats.auroc), ("auprc", stats.auprc)):
            delta = stats.endpoint_delta(endpoint, y, base, treat, null)[0]
            self.assertAlmostEqual(delta, metric(y, treat) - metric(y, base))
            self.assertGreater(delta, 0, endpoint)

    def test_threshold_deltas_match_confusion_counts(self):
        rng = np.random.default_rng(3)
        y = rng.integers(0, 2, 60)
        base, treat = rng.uniform(size=60), rng.uniform(size=60)
        null = np.full(60, y.mean())

        def counts(s):
            call = s >= 0.5
            sens, spec = call[y == 1].mean(), (~call[y == 0]).mean()
            return (sens + spec) / 2, 1 - sens, 1 - spec

        (ba_b, vme_b, me_b), (ba_t, vme_t, me_t) = counts(base), counts(treat)
        expected = {"balanced_accuracy": ba_t - ba_b, "vme": vme_b - vme_t, "me": me_b - me_t}
        for endpoint, value in expected.items():
            self.assertAlmostEqual(stats.endpoint_delta(endpoint, y, base, treat, null)[0], value)
        ll_b, ll_t, ll_0 = (stats.sample_log_loss(y, p).mean() for p in (base, treat, null))
        self.assertAlmostEqual(stats.endpoint_delta("log_loss", y, base, treat, null)[0], ll_b - ll_t)
        self.assertAlmostEqual(stats.endpoint_delta("r2_ll", y, base, treat, null)[0], (ll_b - ll_t) / ll_0)

    def test_score_of_one_half_is_a_resistant_call(self):
        y, null = np.array([1, 0]), np.full(2, 0.5)
        below = np.array([0.4999, 0.4999])
        self.assertEqual(stats.endpoint_delta("vme", y, below, np.array([0.5, 0.4999]), null)[0], 1.0)
        self.assertEqual(stats.endpoint_delta("me", y, below, np.array([0.4999, 0.5]), null)[0], -1.0)

    def test_log_loss_clips_probabilities(self):
        y = np.array([1, 0])
        self.assertAlmostEqual(stats.sample_log_loss(y, np.array([0.0, 1.0]))[0], -np.log(1e-4))
        self.assertAlmostEqual(stats.sample_log_loss(y, np.array([0.0, 1.0]))[1], -np.log(1e-4))

    def test_missing_class_is_undefined_with_reason(self):
        y = np.ones(6, dtype=int)
        s = np.linspace(0.1, 0.9, 6)
        for endpoint in ("auroc", "auprc", "balanced_accuracy", "r2_ll", "me"):
            estimate, reason = stats.endpoint_delta(endpoint, y, s, s[::-1], np.full(6, 0.5))
            self.assertIsNone(estimate, endpoint)
            self.assertIn("no susceptible (y=0) samples", reason)
        for endpoint in ("log_loss", "vme"):
            self.assertIsNotNone(stats.endpoint_delta(endpoint, y, s, s[::-1], np.full(6, 0.5))[0], endpoint)
        self.assertIn("resistant", stats.endpoint_delta("vme", 1 - y, s, s, s)[1])


class InferenceTests(unittest.TestCase):
    def test_exact_enumeration_matches_brute_force(self):
        sums = np.random.default_rng(1).normal(0.1, 1, 10)
        observed = abs(sums.sum())
        brute = np.mean([abs(np.dot(signs, sums)) >= observed - 1e-12
                         for signs in itertools.product((1, -1), repeat=10)])
        p, method = stats.sign_flip_p(sums)
        self.assertEqual(method, "exact")
        self.assertAlmostEqual(p, brute)

    def test_exact_up_to_twenty_clusters_then_monte_carlo(self):
        rng = np.random.default_rng(2)
        p, method = stats.sign_flip_p(rng.normal(0.2, 1, 20))
        self.assertEqual(method, "exact")
        self.assertAlmostEqual(p * 2 ** 20, round(p * 2 ** 20), places=6)
        p, method = stats.sign_flip_p(rng.normal(0.2, 1, 21), n_flips=20000, seed=5)
        self.assertEqual(method, "monte_carlo")
        self.assertAlmostEqual(p * 20001, round(p * 20001), places=6)

    def test_monte_carlo_p_counts_the_observed_pattern(self):
        # Only the identity and the full swap reach |observed|; 99 random draws of 2^15 patterns miss both.
        y = np.tile([1, 0], 15)
        paired = arrays(y, np.repeat(np.arange(15), 2), np.where(y == 1, 0.1, 0.9), np.where(y == 1, 0.9, 0.1))
        for endpoint in ("auroc", "auprc"):
            self.assertEqual(stats.arm_swap_p(endpoint, paired, n_flips=99, seed=1), (1 / 100, "monte_carlo"))
        self.assertEqual(stats.sign_flip_p(np.ones(30), n_flips=99, seed=1), (1 / 100, "monte_carlo"))

    def test_monte_carlo_is_seeded(self):
        sums = np.random.default_rng(4).normal(0.2, 1, 30)
        self.assertEqual(stats.sign_flip_p(sums, n_flips=999, seed=7), stats.sign_flip_p(sums, n_flips=999, seed=7))

    def test_arm_swap_equals_sign_flip_for_a_sum_endpoint(self):
        rng = np.random.default_rng(6)
        paired = arrays(*clustered(rng, 8, m=5, base_slope=0.0, treat_slope=2.0))
        for endpoint in ("log_loss", "vme", "balanced_accuracy"):
            swap = stats.arm_swap_p(endpoint, paired, n_flips=2 ** 8)
            flip = stats.randomisation_p(endpoint, paired)
            with self.subTest(endpoint=endpoint):
                self.assertEqual((swap[1], flip[1]), ("exact", "exact"))
                self.assertLess(flip[0], 1.0)  # a real difference, so swapping the arms (and their calls) matters
                self.assertAlmostEqual(swap[0], flip[0])

    def test_sign_flip_p_is_uniform_under_the_null(self):
        rng = np.random.default_rng(11)
        for k, n_flips in ((12, 20000), (30, 499)):
            p = []
            for _ in range(200):
                y, cluster, _, _ = clustered(rng, k, m=6)
                eta = np.where(y == 1, 0.8, -0.8)
                arm = lambda: sigmoid(eta + np.repeat(rng.normal(0, 0.7, k), 6) + rng.normal(0, 0.7, k * 6))
                p.append(stats.randomisation_p("log_loss", arrays(y, cluster, arm(), arm()),
                                               n_flips=n_flips, seed=len(p))[0])
            with self.subTest(clusters=k):
                self.assertLess(ks_uniform(np.array(p)), 0.115)  # KS 1% critical value at n = 200

    def test_sign_flip_has_power_against_a_real_gain(self):
        rng = np.random.default_rng(12)
        p = [stats.randomisation_p("log_loss", arrays(*clustered(rng, 30, base_slope=0.0, treat_slope=1.0)), n_flips=999, seed=i)[0] for i in range(40)]
        self.assertGreater(np.mean(np.array(p) < 0.05), 0.8)

    def test_t_quantiles_match_tables(self):
        for df in T975:
            self.assertAlmostEqual(stats.t_quantile(0.975, df), T975[df], places=8)
            self.assertAlmostEqual(stats.t_quantile(0.95, df), T95[df], places=8)

    def test_jackknife_intervals_have_the_stated_levels(self):
        for k in (5, 10, 20):
            paired = arrays(*clustered(np.random.default_rng(k), k, m=6))
            loo = []
            for i in range(k):
                keep = paired.cluster_index != i
                loo.append(stats.endpoint_delta("log_loss", paired.y[keep], paired.base[keep], paired.treat[keep],
                                                paired.null[keep])[0])
            se = np.sqrt((k - 1) / k * np.sum((np.array(loo) - np.mean(loo)) ** 2))
            r, = stats.compare(paired, family="D", endpoints=["log_loss"], margins=MARGINS, development=())
            with self.subTest(clusters=k):
                self.assertEqual(r["ci_method"], "cluster_jackknife_t")
                np.testing.assert_allclose(r["ci"], [r["estimate"] - T975[k - 1] * se, r["estimate"] + T975[k - 1] * se])
                np.testing.assert_allclose(r["ci_tost"], [r["estimate"] - T95[k - 1] * se, r["estimate"] + T95[k - 1] * se])
                self.assertTrue(r["ci"][0] < r["ci_tost"][0] < r["ci_tost"][1] < r["ci"][1])

    def test_jackknife_covers_the_population_delta(self):
        rng = np.random.default_rng(13)
        y, _, base, treat = clustered(rng, 50000)
        truth = stats.endpoint_delta("log_loss", y, base, treat, np.full(len(y), 0.5))[0]
        covered = 0
        for _ in range(100):
            r, = stats.compare(arrays(*clustered(rng, 40)), family="D", endpoints=["log_loss"], margins=MARGINS,
                               development=())
            covered += r["ci"][0] <= truth <= r["ci"][1]
        self.assertGreaterEqual(covered / 100, 0.90)

    def test_leave_one_out_recomputes_the_null_and_marks_undefined(self):
        y, base, treat = np.array([1, 0, 1, 1, 0, 1]), np.full(6, 0.6), np.full(6, 0.7)
        null = np.array([0.2, 0.3, 0.4, 0.9, 0.8, 0.7])
        # Cluster 1 holds no susceptible sample, so leaving out cluster 0 or 2 keeps one; cluster 3 holds both.
        paired = arrays(y, [0, 0, 1, 1, 2, 2], base, treat, null)
        loo = stats.leave_one_out(paired, ["r2_ll", "log_loss", "me"])
        keep = np.array([0, 1, 4, 5])
        self.assertAlmostEqual(loo["r2_ll"][1], stats.endpoint_delta("r2_ll", y[keep], base[keep], treat[keep], null[keep])[0])
        self.assertFalse(np.isnan(loo["log_loss"]).any())
        lone = arrays([1, 1, 0, 1], [0, 0, 1, 1], [0.6] * 4, [0.7] * 4)
        self.assertEqual(np.isnan(stats.leave_one_out(lone, ["me"])["me"]).tolist(), [False, True])


class MultiplicityTests(unittest.TestCase):
    def test_holm_and_bh_match_reference_values(self):
        # Hand-computed; these equal R's p.adjust(method = "holm" | "BH").
        cases = [
            ([0.01, 0.04, 0.03, 0.005], [0.03, 0.06, 0.06, 0.02], [0.02, 0.04, 0.04, 0.02]),
            ([0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205],
             [0.008, 0.056] + [0.234] * 6,
             [0.008, 0.032, 0.0672, 0.0672, 0.0672, 0.08, 0.592 / 7, 0.205]),
            ([0.5, 0.9], [1.0, 1.0], [0.9, 0.9]),
        ]
        for p, holm, bh in cases:
            np.testing.assert_allclose(stats.holm(p), holm)
            np.testing.assert_allclose(stats.bh(p), bh)

    def test_bh_never_drops_below_p(self):
        p = np.random.default_rng(14).uniform(size=500)
        self.assertTrue(np.all(np.array(stats.bh(p)) >= p))
        self.assertTrue(np.all(np.array(stats.holm(p)) >= p))

    def test_finalize_adjusts_each_family_on_its_own(self):
        records = [fake_pooled("A", "t1", 0.5), fake_pooled("B", "t1", 0.01), fake_pooled("B", "t2", 0.04),
                   fake_contrast("C", "d1", 0.01), fake_contrast("C", "d2", 0.04), fake_contrast("C", "d3", 0.03),
                   fake_contrast("D", "d1", None, endpoint="auroc")]
        out = stats.finalize(records, ids(records))
        self.assertEqual([(r["adjust"], r["p_adj"]) for r in out],
                         [("holm", 0.5), ("holm", 0.02), ("holm", 0.04), ("bh", 0.03), ("bh", 0.04), ("bh", 0.04),
                          ("none", None)])
        # The A record (p_adj 0.5) and the D record (no p) are inconclusive; there is no equivalent.
        self.assertEqual([r["decision"] for r in out], ["inconclusive", "improves", "improves", "improves", "improves",
                                                        "improves", "inconclusive"])
        self.assertEqual([tuple(r) for r in out], [stats.POOLED_KEYS] * 3 + [stats.RECORD_KEYS] * 4)

    def test_pooled_families_use_holm_not_bh(self):
        # Holm gives [0.03, 0.04, 0.04]; BH would give [0.03, 0.03, 0.03].
        records = [fake_pooled("B", f"t{i}", p) for i, p in enumerate((0.01, 0.02, 0.03))]
        out = stats.finalize(records, ids(records))
        np.testing.assert_allclose([r["p_adj"] for r in out], [0.03, 0.04, 0.04])
        self.assertEqual({r["adjust"] for r in out}, {"holm"})

    def test_finalize_decides_on_the_95pct_ci(self):
        # The 90% ci_tost excludes 0 but the 95% ci does not, so p 0.01 is not enough to improve.
        r = fake_contrast("C", "d1", 0.01, ci=(-0.001, 0.03), tost=(0.002, 0.028))
        out, = stats.finalize([r], ids([r]))
        self.assertEqual((out["p_adj"], out["decision"]), (0.01, "inconclusive"))

    def test_undefined_and_missing_hypotheses_count_at_p_one(self):
        records = [fake_contrast("C", "d1", 0.01), fake_contrast("C", "d2", 0.04),
                   fake_contrast("C", "d3", None, estimate=None)]
        planned = ids(records) + ["d4/C:lr:lineage:base->treat:log_loss"]
        out = stats.finalize(records, planned, missing={planned[-1]: "pairing refused: fold_plan"})
        np.testing.assert_allclose([r["p_adj"] for r in out[:2]], stats.bh([0.01, 0.04, 1.0, 1.0])[:2])
        # d2's raw p 0.04 would improve; its p_adj 0.08 does not.
        self.assertEqual([r["decision"] for r in out[:2]], ["improves", "inconclusive"])
        self.assertEqual((out[2]["p_adj"], out[2]["decision"]), (None, "inconclusive"))

    def test_undefined_and_missing_pooled_hypotheses_enter_holm_at_p_one(self):
        # Contract section 13: in families A and B an undefined pooled record and a missing pooled id
        # each count as one more hypothesis at p = 1, so 0.015625 is doubled, not left alone.
        for family in ("A", "B"):
            tested = fake_pooled(family, "t1", 0.015625)
            undefined = fake_pooled(family, "t2", None, estimate=None, n_units=0)
            with self.subTest(family=family, case="undefined record"):
                out = stats.finalize([tested, undefined], ids([tested, undefined]))
                self.assertEqual([(r["p_adj"], r["decision"]) for r in out],
                                 [(0.03125, "improves"), (None, "inconclusive")])
            with self.subTest(family=family, case="missing id"):
                absent = f"{family}:lr:lineage:base->t3:r2_ll"
                out, = stats.finalize([tested], ids([tested]) + [absent], missing={absent: "no datasets paired"})
                self.assertEqual((out["p_adj"], out["decision"]), (0.03125, "improves"))
            with self.subTest(family=family, case="both"):
                absent = f"{family}:lr:lineage:base->t3:r2_ll"
                out = stats.finalize([tested, undefined], ids([tested, undefined]) + [absent],
                                     missing={absent: "no datasets paired"})
                self.assertEqual(out[0]["p_adj"], 0.046875)

    def test_finalize_refusals(self):
        good = [fake_contrast("C", "d1", 0.01), fake_contrast("C", "d2", 0.04)]
        extra = "d3/C:lr:lineage:base->treat:log_loss"
        cases = {
            "incomplete family": (good, ids(good) + [extra], None),
            "unplanned record": (good, ids(good)[:1], None),
            "missing id that has a record": (good, ids(good), {ids(good)[0]: "x"}),
            "missing id without a reason": (good, ids(good) + [extra], {extra: " "}),
            "contrast in family A": ([fake_contrast("A", "d1", 0.01)], ["d1/A:lr:lineage:base->treat:log_loss"], None),
            "pooled in family C": ([fake_pooled("C", "t", 0.01)], ["C:lr:lineage:base->t:r2_ll"], None),
            "secondary endpoint in C": ([fake_contrast("C", "d1", 0.01, endpoint="auroc")],
                                        ["d1/C:lr:lineage:base->treat:auroc"], None),
            "family D with p": ([fake_contrast("D", "d1", 0.01)], ["d1/D:lr:lineage:base->treat:log_loss"], None),
            "family C without p": ([fake_contrast("C", "d1", None)], ["d1/C:lr:lineage:base->treat:log_loss"], None),
            "repeated record": (good + good[:1], ids(good), None),
        }
        for name, (records, planned, missing) in cases.items():
            with self.subTest(name), self.assertRaises(ValueError):
                stats.finalize(records, planned, missing)
        with self.assertRaises(TypeError):
            stats.finalize(good, ids(good)[0])
        with self.assertRaises(ValueError):
            stats.finalize(stats.finalize(good, ids(good)), ids(good))  # already final


class DecisionTests(unittest.TestCase):
    def test_boundaries(self):
        cases = [
            ([1e-9, 0.2], 0.049, "improves"),
            ([0.0, 0.2], 0.001, "inconclusive"),        # lower bound must exceed 0
            ([0.01, 0.2], 0.05, "inconclusive"),        # p_adj must be below 0.05
            ([-0.2, -1e-9], 0.049, "worse"),
            ([-0.2, 0.0], 0.001, "inconclusive"),       # upper bound must be below 0
            ([-0.2, -0.01], 0.05, "inconclusive"),
            ([-0.2, 0.2], 0.001, "inconclusive"),       # straddles 0
            ([-0.2, 0.2], 0.5, "inconclusive"),
            ([0.01, 0.05], None, "inconclusive"),       # family D: no p, so never improves
            ([-0.05, -0.01], None, "inconclusive"),
            (None, None, "inconclusive"),
        ]
        for ci, p_adj, expected in cases:
            with self.subTest(ci=ci, p_adj=p_adj):
                self.assertEqual(stats.decide(ci, p_adj), expected)
        self.assertEqual(stats.DECISIONS, ("improves", "worse", "inconclusive"))

    def test_interval_inside_margin_is_inconclusive_not_equivalent(self):
        # Contract section 15: a 90% interval strictly inside +/-margin at 60 or more clusters,
        # which v1.7 called equivalent, now decides inconclusive in families C and D.
        for k in (60, 120):
            r, = stats.compare(balanced(k, seed=k), family="D", endpoints=["log_loss"], margins={"log_loss": 0.5},
                               development=())
            out, = stats.finalize([r], [r["id"]])
            with self.subTest(clusters=k):
                self.assertGreaterEqual(out["n_units"], 60)
                self.assertTrue(-0.5 < out["ci_tost"][0] and out["ci_tost"][1] < 0.5)
                self.assertEqual(out["decision"], "inconclusive")
        c = fake_contrast("C", "d1", 0.5, ci=(-0.02, 0.03), tost=(-0.01, 0.02), n_units=200, margin=0.05)
        out, = stats.finalize([c], ids([c]))
        self.assertEqual((out["p_adj"], out["decision"]), (0.5, "inconclusive"))
        records = [pending(f"x{i}", 0.001 * (i % 3), margin=0.5) for i in range(40)]
        raw = stats.pool(records, species={f"x{i}": f"S{i}" for i in range(40)}, family="D", margin=0.5,
                         development=())
        self.assertTrue(-0.5 < raw["ci_tost"][0] and raw["ci_tost"][1] < 0.5)
        self.assertEqual(stats.finalize([raw], [raw["contrast"]])[0]["decision"], "inconclusive")


class RecordTests(unittest.TestCase):
    def test_compare_and_finalize_produce_documented_records(self):
        paired = pair(make_rows())
        for family, endpoints, adjust in (("C", ["log_loss"], "bh"), ("D", list(stats.ENDPOINTS), "none")):
            raw = stats.compare(paired, family=family, endpoints=endpoints, margins=MARGINS, development=set(),
                                n_flips=499, seed=3)
            for r in stats.finalize(raw, ids(raw)):
                with self.subTest(family=family, endpoint=r["endpoint"]):
                    self.assertEqual(tuple(r), stats.RECORD_KEYS)
                    self.assertEqual(r["id"], f"ds/{family}:lr:lineage:base->treat:{r['endpoint']}")
                    self.assertEqual(r["dataset"], "ds")
                    json.dumps(r, allow_nan=False)
                    self.assertEqual((r["unit"], r["n_units"], r["few_clusters"]), ("cluster", 6, True))
                    self.assertIs(type(r["n_units"]), int)
                    self.assertIs(type(r["undefined_resamples"]), int)
                    self.assertEqual((r["claim_level"], r["interval_label"]), ("predictive_increment", stats.INTERVAL_LABEL))
                    self.assertEqual((r["adjust"], r["ci_method"], r["undefined_reason"]), (adjust, "cluster_jackknife_t", None))
                    self.assertIn(r["decision"], stats.DECISIONS)
                    self.assertTrue(r["ci"][0] <= r["ci_tost"][0] <= r["ci_tost"][1] <= r["ci"][1])
                    if family == "D":
                        self.assertEqual((r["p"], r["p_adj"]), (None, None))
                    else:
                        self.assertIs(type(r["p"]), float)
                        self.assertGreaterEqual(r["p_adj"], r["p"])
        self.assertEqual(raw, stats.compare(paired, family="D", endpoints=list(stats.ENDPOINTS), margins=MARGINS,
                                            development=set(), n_flips=499, seed=3))

    def test_families_and_endpoints(self):
        paired = pair(make_rows())
        for family, endpoints in (("A", ["log_loss"]), ("B", ["log_loss"]), ("C", ["auroc"]), ("C", ["r2_ll"]),
                                  ("E", ["log_loss"])):
            with self.subTest(family=family, endpoints=endpoints), self.assertRaises(ValueError):
                stats.compare(paired, family=family, endpoints=endpoints, margins=MARGINS, development=())

    def test_ids_must_split_cleanly(self):
        paired = pair(make_rows())
        for field, value in (("dataset", "bad id"), ("dataset", "a/b"), ("learner", "l r"), ("baseline", "x/y"),
                             ("treatment", "x`y"), ("split", "`"), ("learner", "l:r"), ("split", ""),
                             ("baseline", "a->b"), ("treatment", "->")):
            changed = stats.PairedArms(**{**paired.__dict__, field: value})
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                stats.compare(changed, family="D", endpoints=["log_loss"], margins=MARGINS, development=())

    def test_margins_must_be_finite_and_positive(self):
        paired = pair(make_rows())
        for margin in (0, -0.1, float("inf"), float("nan"), True, None, "0.02"):
            with self.subTest(margin=margin), self.assertRaises(ValueError):
                stats.compare(paired, family="D", endpoints=["log_loss"], margins={"log_loss": margin}, development=())

    def test_compare_uses_exact_p_at_twenty_clusters(self):
        paired = arrays(*clustered(np.random.default_rng(8), 20, m=4))
        record, = stats.compare(paired, family="C", endpoints=["log_loss"], margins=MARGINS, development=(), n_flips=99)
        sums = np.bincount(paired.cluster_index, minlength=20,
                           weights=stats.contributions("log_loss", paired.y, paired.base, paired.treat, paired.null))
        self.assertEqual(record["p"], stats.sign_flip_p(sums, exact_max=20)[0])
        self.assertAlmostEqual(record["p"] * 2 ** 20, round(record["p"] * 2 ** 20), places=6)
        self.assertFalse(record["few_clusters"])

    def test_undefined_endpoint_record(self):
        paired = arrays([1] * 8, [0, 0, 1, 1, 2, 2, 3, 3], [0.6, 0.4, 0.3, 0.7] * 2, [0.8, 0.6, 0.2, 0.4] * 2)
        ll, au = stats.finalize(raw := stats.compare(paired, family="D", endpoints=["log_loss", "auroc"],
                                                     margins=MARGINS, development=()), ids(raw))
        self.assertIsNotNone(ll["estimate"])
        self.assertEqual(ll["undefined_resamples"], 0)
        self.assertEqual((au["estimate"], au["ci"], au["ci_tost"], au["p"], au["p_adj"]), (None,) * 5)
        self.assertEqual((au["undefined_resamples"], au["ci_method"], au["decision"]), (4, "none", "inconclusive"))
        self.assertEqual(au["undefined_reason"], "auroc undefined: no susceptible (y=0) samples")

    def test_one_cluster_is_undefined(self):
        r, = stats.compare(arrays([1, 0, 1, 0], [0] * 4, [0.6, 0.4, 0.3, 0.7], [0.8, 0.2, 0.6, 0.4]), family="C",
                           endpoints=["log_loss"], margins=MARGINS, development=())
        self.assertEqual((r["estimate"], r["ci"], r["p"], r["ci_method"]), (None, None, None, "none"))
        self.assertIn("one cluster", r["undefined_reason"])

    def test_more_than_ten_percent_undefined_resamples(self):
        # The only susceptible samples sit in cluster 0, so leaving it out undefines auroc: 1 of k replicates.
        for k, share_rule in ((5, True), (9, True), (10, False), (12, False)):
            y = np.ones(3 * k, dtype=int)
            y[:2] = 0
            paired = arrays(y, np.repeat(np.arange(k), 3), np.linspace(0.2, 0.8, 3 * k), np.linspace(0.3, 0.9, 3 * k))
            au, ll = stats.compare(paired, family="D", endpoints=["auroc", "log_loss"], margins=MARGINS, development=())
            with self.subTest(clusters=k):
                self.assertEqual(au["undefined_resamples"], 1)
                self.assertEqual((au["estimate"], au["ci"], au["ci_tost"]), (None, None, None))
                self.assertEqual(au["undefined_reason"] == stats.UNDEFINED_SHARE_REASON, share_rule)
                self.assertIn("1 of" if not share_rule else "10%", au["undefined_reason"])
                self.assertEqual((ll["undefined_resamples"], ll["undefined_reason"]), (0, None))

    def test_development_datasets_are_refused_in_confirmatory_families(self):
        paired = pair(make_rows())
        with self.assertRaises(stats.DevelopmentDatasetError):
            stats.compare(paired, family="C", endpoints=["log_loss"], margins=MARGINS, development={"ds"})
        stats.compare(paired, family="D", endpoints=["log_loss"], margins=MARGINS, development={"ds"})
        with self.assertRaises(TypeError):
            stats.compare(paired, family="C", endpoints=["log_loss"], margins=MARGINS, development="ds")


def pending(dataset, estimate, ci=None, reason=None, **changes):
    """The fields of a compare() r2_ll record that pool() and dataset_summary() read."""
    return {"dataset": dataset, "baseline": "base", "treatment": "treat", "learner": "lr", "split": "lineage",
            "endpoint": "r2_ll", "family": "D", "margin": 0.025, "estimate": estimate, "ci": ci,
            "undefined_reason": reason, **changes}


class PoolTests(unittest.TestCase):
    RECORDS = [pending("x1", 0.01, [0.005, 0.02]), pending("x2", 0.03, [0.02, 0.04]),
               pending("x3", 0.05, [0.03, 0.08]), pending("y1", 0.10, [0.06, 0.14]),
               pending("y2", None, None, "r2_ll undefined: no resistant (y=1) samples"),
               pending("z1", -0.2, [-0.3, -0.1])]
    SPECIES = {"x1": "X", "x2": "X", "x3": "X", "y1": "Y", "y2": "Y", "z1": "Z"}

    def test_pool_averages_within_species_and_lists_exclusions(self):
        raw = stats.pool(self.RECORDS[:5], species=self.SPECIES, family="A", margin=0.025, development=(),
                         missing={"w1": "pairing refused: fold_plan"})
        self.assertAlmostEqual(raw["estimate"], (0.03 + 0.10) / 2)  # the dataset mean would be 0.0475
        self.assertEqual(raw["p"], 0.5)  # two positive species means: only ++ and -- reach |observed|
        half = 0.07 / 2  # jackknife SE of a mean of two is |difference| / 2
        np.testing.assert_allclose(raw["ci"], [0.065 - T975[1] * half, 0.065 + T975[1] * half])
        np.testing.assert_allclose(raw["ci_tost"], [0.065 - T95[1] * half, 0.065 + T95[1] * half])
        pooled, = stats.finalize([raw], [raw["contrast"]])
        self.assertEqual(tuple(pooled), stats.POOLED_KEYS)
        json.dumps(pooled, allow_nan=False)
        self.assertEqual((pooled["contrast"], pooled["family"], pooled["adjust"], pooled["margin"]),
                         ("A:lr:lineage:base->treat:r2_ll", "A", "holm", 0.025))
        self.assertEqual((pooled["unit"], pooled["n_units"], pooled["n_datasets"], pooled["n_excluded"]),
                         ("species", 2, 6, 2))
        self.assertEqual(pooled["excluded"], [{"dataset": "w1", "reason": "pairing refused: fold_plan"},
                                              {"dataset": "y2", "reason": "r2_ll undefined: no resistant (y=1) samples"}])
        self.assertEqual((pooled["p_adj"], pooled["decision"]), (0.5, "inconclusive"))

    def test_one_species_gives_no_interval(self):
        raw = stats.pool(self.RECORDS[:3], species=self.SPECIES, family="B", margin=0.025, development=())
        self.assertEqual((raw["estimate"], raw["ci"], raw["ci_tost"], raw["p"], raw["n_units"], raw["n_excluded"]),
                         (None, None, None, None, 0, 3))
        self.assertTrue(all("one contributing species" in e["reason"] for e in raw["excluded"]))
        self.assertEqual(stats.finalize([raw], [raw["contrast"]])[0]["decision"], "inconclusive")

    def test_pool_with_nothing_contributing_is_inconclusive(self):
        raw = stats.pool([pending("y2", None, None, "undefined")], species=self.SPECIES, family="A", margin=0.025,
                         development=())
        pooled, = stats.finalize([raw], [raw["contrast"]])
        self.assertEqual((pooled["estimate"], pooled["ci"], pooled["ci_tost"], pooled["p"], pooled["p_adj"],
                          pooled["n_units"]), (None,) * 5 + (0,))
        self.assertEqual((pooled["n_excluded"], pooled["decision"]), (1, "inconclusive"))

    def test_dataset_summary(self):
        summary = stats.dataset_summary(self.RECORDS, 0.025)
        self.assertEqual(summary["median_dataset_delta"], 0.03)
        self.assertEqual(summary["dataset_ci_counts"],
                         {"above": 2, "below": 1, "inside": 1, "straddling": 1, "undefined": 1})

    def test_pool_refusals(self):
        for kwargs in ({"development": {"x1"}}, {"development": {"w1"}, "missing": {"w1": "absent"}}):
            with self.subTest(**kwargs), self.assertRaises(stats.DevelopmentDatasetError):
                stats.pool(self.RECORDS, species=self.SPECIES, family="A", margin=0.025, **kwargs)
        x1, y1 = self.RECORDS[0], self.RECORDS[3]
        cases = {
            "pooled family C": ([x1, y1], "C", 0.025),
            "inputs disagree on family": ([x1, dict(y1, family="C")], "A", 0.025),
            "inputs disagree on margin": ([x1, dict(y1, margin=0.03)], "A", 0.025),
            "input margin is not the pooled margin": ([x1, y1], "A", 0.03),
            "inputs disagree on learner": ([x1, dict(y1, learner="rf")], "A", 0.025),
            "endpoint is not r2_ll": ([dict(x1, endpoint="log_loss"), dict(y1, endpoint="log_loss")], "A", 0.025),
            "infinite margin": ([dict(x1, margin=float("inf")), dict(y1, margin=float("inf"))], "A", float("inf")),
            "no dataset": ([{k: v for k, v in x1.items() if k != "dataset"}], "A", 0.025),
        }
        for name, (records, family, margin) in cases.items():
            with self.subTest(name), self.assertRaises(ValueError):
                stats.pool(records, species=self.SPECIES, family=family, margin=margin, development=())


def load_k_schema():
    """Task K's schema.py, read-only: this repo's copy once K is merged, else
    the file named by GRM_TASK_K_SCHEMA. Compiled from source so no bytecode
    is written beside it."""
    for candidate in (ROOT / "scripts" / "ablation" / "schema.py", os.environ.get("GRM_TASK_K_SCHEMA")):
        if candidate and Path(candidate).is_file():
            module = types.ModuleType("grm_ablation_schema_k")
            module.__file__ = str(candidate)
            exec(compile(Path(candidate).read_text(encoding="utf-8"), str(candidate), "exec"), module.__dict__)
            return module
    return None


K_SCHEMA = load_k_schema()


def k_cases():
    """finalize() output covering each family, decision and null shape J emits."""
    r2 = [stats.compare(balanced(12, gap=0.1 + 0.05 * i, seed=i, dataset=f"d{i}"), family="D", endpoints=["r2_ll"],
                        margins=MARGINS, development=())[0] for i in range(4)]
    species = {"d0": "S0", "d1": "S1", "d2": "S2", "d3": "S2"}
    all_r = arrays([1] * 24, np.repeat(np.arange(12), 2), np.full(24, 0.6), np.full(24, 0.7), dataset="allR")
    y = np.ones(15, dtype=int)
    y[:2] = 0
    lone = arrays(y, np.repeat(np.arange(5), 3), np.linspace(0.2, 0.8, 15), np.linspace(0.3, 0.9, 15), dataset="lone")
    pooled = [stats.pool(r2, species=species, family="A", margin=0.02, development=()),
              stats.pool(r2, species=species, family="D", margin=0.02, development=()),
              stats.pool(r2[2:], species=species, family="B", margin=0.02, development=(), missing={"d9": "absent"}),
              stats.pool([stats.compare(all_r, family="D", endpoints=["r2_ll"], margins=MARGINS, development=())[0]],
                         species={"allR": "E"}, family="B", margin=0.02, development=(), missing={"d8": "absent"})]
    pooled[2]["contrast"] = pooled[2]["contrast"].replace("treat", "treat2")
    contrasts = [
        *[stats.compare(balanced(12, gap=0.3, seed=i, dataset=f"c{i}"), family="C", endpoints=["log_loss"],
                        margins=MARGINS, development=())[0] for i in range(3)],
        stats.compare(balanced(12, gap=-0.3, seed=7, dataset="cw"), family="C", endpoints=["log_loss"],
                      margins=MARGINS, development=())[0],
        stats.compare(arrays([1, 0] * 2, [0] * 4, [0.6, 0.4] * 2, [0.7, 0.3] * 2, dataset="one"), family="C",
                      endpoints=["log_loss"], margins=MARGINS, development=())[0],
        *stats.compare(balanced(12, seed=5, dataset="eq12"), family="D", endpoints=list(stats.ENDPOINTS),
                       margins=dict.fromkeys(stats.ENDPOINTS, 0.5), development=()),
        *stats.compare(balanced(60, seed=5, dataset="eq60"), family="D", endpoints=["log_loss", "auroc"],
                       margins=dict.fromkeys(stats.ENDPOINTS, 0.5), development=()),
        *stats.compare(all_r, family="D", endpoints=list(stats.ENDPOINTS), margins=MARGINS, development=()),
        *stats.compare(lone, family="D", endpoints=["auroc", "vme"], margins=MARGINS, development=()),
        *stats.compare(balanced(6, seed=6, dataset="few"), family="D", endpoints=["log_loss"],
                       margins={"log_loss": 0.5}, development=()),
        *r2,
    ]
    missing = {"c9/C:lr:lineage:base->treat:log_loss": "pairing refused: fold_plan"}
    return stats.finalize(pooled + contrasts, ids(pooled + contrasts) + list(missing), missing)


@unittest.skipIf(K_SCHEMA is None, "Task K's schema.py not found; set GRM_TASK_K_SCHEMA to its path")
class KSchemaTests(unittest.TestCase):
    def test_finalize_output_validates_with_k(self):
        records = k_cases()
        self.assertEqual({r["decision"] for r in records}, set(stats.DECISIONS))
        # Any refusal fails: a K that still expects equivalent breaks contract section 15.
        for r in records:
            kind = "pooled" if "contrast" in r else "contrast"
            with self.subTest(record=stats._name(r)):
                getattr(K_SCHEMA, f"validate_{kind}")(json.loads(json.dumps(r, allow_nan=False)))
        print(f"\nK schema {K_SCHEMA.__file__}: {len(records)} records validated", file=sys.stderr)


if __name__ == "__main__":
    unittest.main()
