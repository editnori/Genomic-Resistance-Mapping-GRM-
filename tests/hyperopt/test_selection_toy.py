"""Toy-data acceptance criteria for a learner-hyperparameter selection rule
(docs/hyperopt-selection.md), applied to select.py's current rule: the
one-standard-error window with the lower-complexity preference.

Scope: learner hyperparameters only, the per-family knobs in protocol.py.
Nothing here covers choosing the input-data configuration (cohort and genome
count, which genomes are included, k-mer length, feature families or the
retained-feature budget), which is a separate open capability.

Generic toy data only. A "configuration" is an index with a known true
held-out score; fold scores are that score plus Gaussian (or binomial) noise,
with a fold effect shared by every configuration. There are no genomes,
learners, features or real data here, and nothing is fitted. Every draw comes
from a fixed numpy seed, so results do not depend on PYTHONHASHSEED.

The Monte Carlo runs a vectorised copy of the rule (`one_se`). RuleCopy checks
that copy against select.select_history on drawn histories, flat and Task B,
so a change to select.py's rule fails that test first.

The bounds are proposals, not adopted settings. Each criterion test pins
today's value on its scenario and whether that value meets the bound, so a
miss is asserted as a miss, by its size. A future rule that changes a value
fails the pin; update the pin with the evidence the doc's A6 asks for.

usage: python3 -I -B -W error -m unittest discover -s tests/hyperopt -t . -p test_selection_toy.py
"""
import importlib.util
import inspect
import math
import sys
import unittest
from pathlib import Path

import numpy as np

from scripts.hyperopt import protocol as hp
from scripts.hyperopt import search

ROOT = Path(__file__).resolve().parents[2]


def _load(name, path):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


sel = _load("grm_hyperopt_select", ROOT / "scripts/hyperopt/select.py")

# The toy problem's shape, set to the archived synthetic run's (F4): 32
# configurations per family, 3 inner folds, and a per-fold score SD of about
# 0.045, which is the binomial SD of a balanced accuracy near 0.8 on about 84
# validation samples. RHO is the share of that variance that is a fold effect
# shared by every configuration (the archived run's median was about 0.3).
K, FOLDS, SD, RHO = 32, 3, 0.045, 0.3
N_VAL = 84            # validation samples per fold, binomial noise only
REPS = 4000           # Monte Carlo draws per scenario
SEED = 20261102

# Acceptance bounds (docs/hyperopt-selection.md, section 4). Scores are on the
# objective's own scale, here a balanced-accuracy-like score in [0, 1].
HARM_BOUND = 0.005    # A2: expected held-out loss against the fixed default
WIN_GAP = 0.08        # A3: the planted improvement over the default
WIN_SHARE = 0.5       # A3: share of WIN_GAP the rule must recover
PREF_BOUND = 0.002    # A4: cost of the preference against a uniform pick
CAL_BOUND = 0.01      # A5: bias of a reported score for the selected config

# A pinned value must reproduce within PIN_SES Monte Carlo SEs (or 1e-4, the
# pins' rounding), and every verdict must clear its bound by CLEAR_SES.
PIN_SES, CLEAR_SES = 3, 6


# ----------------------------------------------------------------- scenarios

def _capacity(default=16.0):
    """Hidden capacity of configurations 1..K-1, spread over 2..40, and of
    the default (configuration 0)."""
    return np.concatenate([[default], 2.0 + 38.0 * np.arange(K - 1) / (K - 2)])


def _plateau(h):
    """True score rising with capacity to a plateau (an underfitting side).
    A default at capacity 28 is 0.007 below the best configuration (40)."""
    return 0.88 - 0.30 * np.exp(-h / 8.0)


def _unrelated(h):
    """A complexity proxy that is a fixed shuffle of capacity: it carries no
    information about it (like trees kept under early stopping, which falls as
    the learning rate rises)."""
    order = np.array([0] + [1 + (7 * i) % (K - 1) for i in range(K - 1)])
    return h[order]


def scenario(name):
    """The toy problems. Each has the true held-out score mu of every
    configuration, the location nu its inner folds measure (nu = mu unless
    the scenario plants an inner-to-outer shift), the complexity each trial
    reports, the shared-fold share rho, and the noise model."""
    h = _capacity()
    s = {"rho": RHO, "noise": "gauss"}
    if name == "null":                      # every configuration is equal
        mu = np.full(K, 0.80)
        kappa = h
    elif name == "null_binomial":
        mu = np.full(K, 0.80)
        kappa, s["noise"] = h, "binomial"
    elif name == "default_best":            # a good default, a worse search
        mu = np.concatenate([[0.85], 0.85 - 0.003 * np.arange(1, K)])
        kappa = _unrelated(h)
    elif name == "underfit":                # simpler is worse; proxy = capacity
        h = _capacity(28.0)
        mu, kappa = _plateau(h), h
    elif name == "overfit":                 # simpler is better past an optimum
        mu = 0.85 - 0.15 * np.exp(-h / 3.0) - 0.0015 * h
        kappa = h.copy()
        kappa[0] = 30.0
        mu[0] = 0.85 - 0.15 * math.exp(-10.0) - 0.0015 * 30.0
    elif name == "unrelated_proxy":         # proxy carries no capacity signal
        h = _capacity(28.0)
        mu, kappa = _plateau(h), _unrelated(h)
    elif name == "shared_folds":            # as underfit, fold effects dominate
        h = _capacity(28.0)
        mu, kappa = _plateau(h), h
        s["rho"] = 0.8
    elif name == "clear_winner":            # one configuration beats the default
        mu = np.concatenate([[0.80], 0.80 - 0.10 * np.arange(1, K) / (K - 1)])
        mu[20] = 0.80 + WIN_GAP
        kappa = h
    else:
        raise KeyError(name)
    s.update(name=name, mu=mu, nu=mu.copy(), kappa=np.asarray(kappa, float))
    return s


# The scenarios A2 and A4 must hold on. The inner folds measure the outer
# target without bias in all of them; see the doc for the shift scenarios.
SUITE = ("null", "default_best", "underfit", "overfit", "unrelated_proxy",
         "shared_folds", "clear_winner")


def draw(s, reps=REPS, seed=SEED):
    """Inner fold scores, shape (reps, K, FOLDS)."""
    rng = np.random.default_rng([seed, SUITE.index(s["name"])
                                 if s["name"] in SUITE else 99])
    tau = SD * math.sqrt(s["rho"])
    fold = rng.normal(0.0, tau, (reps, 1, FOLDS))
    if s["noise"] == "binomial":
        p = np.clip(s["nu"][None, :, None] + fold, 0.0, 1.0)
        return rng.binomial(N_VAL, np.broadcast_to(p, (reps, K, FOLDS))) / N_VAL
    sigma = SD * math.sqrt(1.0 - s["rho"])
    return s["nu"][None, :, None] + fold + rng.normal(0.0, sigma, (reps, K, FOLDS))


# --------------------------------------------------------------------- rules

def preference_rank(kappa):
    """select.preference_key as a rank: lower complexity, then the default
    (configuration 0), then the lowest id."""
    order = sorted(range(len(kappa)), key=lambda c: (kappa[c], c != 0, c))
    rank = np.empty(len(kappa), dtype=np.int64)
    rank[order] = np.arange(len(kappa))
    return rank


def one_se(x, kappa):
    """select.py's rule, vectorised over draws. Returns the selected and the
    best-mean configuration per draw, the candidate mask and the means."""
    reps, k, n = x.shape
    m = x.mean(axis=2)
    se = x.std(axis=2, ddof=1) / math.sqrt(n)
    rank = preference_rank(kappa)
    rows = np.arange(reps)
    best = np.where(m == m.max(axis=1, keepdims=True), rank, k).argmin(axis=1)
    threshold = m[rows, best] - se[rows, best]
    slack = sel.BOUNDARY_RTOL * np.maximum(1.0, np.abs(threshold))
    cand = m >= (threshold - slack)[:, None]
    chosen = np.where(cand, rank, k).argmin(axis=1)
    return chosen, best, cand, m


def outcome(s, x=None):
    """Held-out results on one scenario, per draw: the true score of today's
    selection, of the best-mean trial and of the default; the expected true
    score of a uniform pick from today's candidate set; the reported (inner)
    mean of today's selection; and the candidate count."""
    x = draw(s) if x is None else x
    chosen, best, cand, m = one_se(x, s["kappa"])
    rows = np.arange(len(x))
    return {"one_se": s["mu"][chosen], "best_mean": s["mu"][best],
            "default": np.full(len(x), s["mu"][0]),
            "uniform": (cand * s["mu"]).sum(axis=1) / cand.sum(axis=1),
            "reported": m[rows, chosen], "n_cand": cand.sum(axis=1)}


def mean_and_mcse(v):
    v = np.asarray(v, float)
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))


# ------------------------------------------------------- select.py histories

# A toy fold plan: 12 abstract sample slots in 3 inner folds, the last three
# null because they stand for the outer test fold the selection never sees.
FOLD_PLAN = {"inner_folds": [0, 1, 2, 0, 1, 2, 0, 1, 2, None, None, None]}
FOLD_PLAN_HASH = hp.fold_plan_hash(FOLD_PLAN)

# The flat, test-only protocol (no search plan to replay), and Task B's
# protocol for the test-only generic family, both pinning FOLD_PLAN.
PROTOCOL = sel.Protocol.from_dict(
    {"objective": "score", "direction": "maximize", "protocol_hash": "toy",
     "fold_plan_hash": FOLD_PLAN_HASH, "budgets": {"toy": K}}, trusted_flat=True)
TASK_B = hp.build_protocol(FOLD_PLAN_HASH, families=("generic",),
                           trials_per_family=K, objective_name="score")
TASK_B_PLAN = search.plan_trials(TASK_B)


def history(x, kappa):
    """One toy draw (K x FOLDS) as select.py trial records, flat protocol."""
    return [{"protocol_hash": "toy", "trial_id": c, "family": "toy",
             "config": {"c": c}, "complexity": float(kappa[c]),
             "is_default": c == 0, "seed": 0, "fold_plan_hash": FOLD_PLAN_HASH,
             "status": "ok", "fit_count": FOLDS, "seconds": 0.0,
             "inner_scores": [{"fold": j, "metrics": {"score": float(x[c, j])}}
                              for j in range(FOLDS)]}
            for c in range(K)]


def task_b_history(x, kappa):
    """One toy draw as the history Task B's search writes: search.run_trial on
    every planned generic trial, in run order, with an objective that returns
    the drawn scores of configuration c = trial id."""
    ids = {hp.canonical_json(t["config"]): t["trial_id"] for t in TASK_B_PLAN}

    def objective(family, config, fold_plan, seed):
        c = ids[hp.canonical_json(config)]
        return {"status": "ok", "fit_count": FOLDS, "seconds": 0.0,
                "complexity": float(kappa[c]),
                "inner_scores": [{"fold": j, "metrics": {"score": float(x[c, j])}}
                                 for j in range(FOLDS)]}
    return [search.run_trial(TASK_B, objective, t, FOLD_PLAN) for t in TASK_B_PLAN]


def select_toy(trials, protocol=PROTOCOL):
    return sel.select_history(trials, protocol, FOLD_PLAN)["families"]


# --------------------------------------------------------------------- tests

class RuleCopy(unittest.TestCase):
    """The vectorised copy is select.py's rule on these toy histories."""

    def compare(self, x, s, make, protocol, family):
        chosen, best, cand, m = one_se(x, s["kappa"])
        for r in range(len(x)):
            f = select_toy(make(x[r], s["kappa"]), protocol)[family]
            with self.subTest(scenario=s["name"], draw=r):
                self.assertEqual(f["selected"]["trial_id"], chosen[r])
                self.assertEqual(f["best"]["trial_id"], best[r])
                self.assertEqual(sorted(c["trial_id"] for c in f["candidate_set"]),
                                 sorted(np.flatnonzero(cand[r]).tolist()))
                self.assertAlmostEqual(f["selected"]["mean"], m[r, chosen[r]],
                                       delta=1e-12)

    def test_copy_matches_select_history_on_flat_histories(self):
        for name in SUITE:
            s = scenario(name)
            self.compare(draw(s, reps=25, seed=SEED + 1), s, history, PROTOCOL, "toy")

    def test_copy_matches_select_history_on_task_b_histories(self):
        for name in SUITE:
            s = scenario(name)
            self.compare(draw(s, reps=3, seed=SEED + 4), s, task_b_history,
                         TASK_B, "generic")


class A1UntouchedOuter(unittest.TestCase):
    """A1: the rule sees only the inner history, the protocol and the fold
    plan it pins; it is deterministic. Today's rule meets A1."""

    def test_the_rule_takes_the_history_the_protocol_and_the_fold_plan(self):
        self.assertEqual(list(inspect.signature(sel.select_history).parameters),
                         ["trials", "protocol", "fold_plan"])

    def test_flat_history_same_selection_in_any_record_order(self):
        s = scenario("underfit")
        x = draw(s, reps=5, seed=SEED + 2)
        rng = np.random.default_rng(SEED + 3)
        for r in range(len(x)):
            h = history(x[r], s["kappa"])
            first = select_toy(h)["toy"]
            for _ in range(3):
                again = select_toy([h[i] for i in rng.permutation(K)])["toy"]
                self.assertEqual(again["selected"], first["selected"])
                self.assertEqual(again["candidate_set"], first["candidate_set"])

    def test_task_b_history_out_of_run_order_is_refused_not_reselected(self):
        s = scenario("underfit")
        h = task_b_history(draw(s, reps=1, seed=SEED + 5)[0], s["kappa"])
        with self.assertRaises(sel.HistoryError) as e:
            select_toy(h[::-1], TASK_B)
        self.assertIn("off_plan_order", e.exception.codes)

    def test_a_record_carrying_an_outer_score_is_refused(self):
        # The record schema allows no field beyond Task B's, so an outer
        # score cannot ride into the selection inside a trial.
        s = scenario("underfit")
        for make, protocol in ((history, PROTOCOL), (task_b_history, TASK_B)):
            h = make(draw(s, reps=1, seed=SEED + 6)[0], s["kappa"])
            h[1] = dict(h[1], outer_score=0.9)
            with self.subTest(protocol=type(protocol).__name__):
                with self.assertRaises(sel.HistoryError) as e:
                    select_toy(h, protocol)
                self.assertEqual(e.exception.codes, {"unknown_field"})

    def test_task_b_needs_the_fold_plan_it_pins(self):
        s = scenario("underfit")
        h = task_b_history(draw(s, reps=1, seed=SEED + 7)[0], s["kappa"])
        with self.assertRaises(sel.HistoryError) as e:
            sel.select_history(h, TASK_B)
        self.assertIn("missing_fold_plan", e.exception.codes)
        other = {"inner_folds": FOLD_PLAN["inner_folds"][:-1] + [0]}
        with self.assertRaises(sel.HistoryError) as e:
            sel.select_history(h, TASK_B, other)
        self.assertIn("bad_fold_plan", e.exception.codes)


class Pinned(unittest.TestCase):
    """Shared check for A2 to A5: today's value on one scenario, and whether
    it meets the criterion. `margin` is how far the value lies inside the
    bound; a negative margin is the size of the miss."""

    def pin(self, label, value, mcse, margin, today, meets, precision):
        self.assertLess(mcse, precision, f"{label}: Monte Carlo SE too large")
        self.assertAlmostEqual(
            value, today, delta=max(PIN_SES * mcse, 1e-4),
            msg=f"{label}: today's value moved from {today} to {value:.6f}")
        self.assertEqual(margin >= 0, meets,
                         f"{label}: {value:.6f}, margin to the bound {margin:+.6f}")
        self.assertGreaterEqual(abs(margin), CLEAR_SES * mcse,
                                f"{label}: verdict within {CLEAR_SES} SEs of the bound")


class A2NoHarmAgainstDefault(Pinned):
    """A2: on every suite scenario, the expected held-out loss of the selected
    configuration against the fixed default, E[mu_default - mu_selected], is
    at most HARM_BOUND. Today's rule misses on three scenarios."""

    def check(self, name, today, meets):
        o = outcome(scenario(name))
        loss, mcse = mean_and_mcse(o["default"] - o["one_se"])
        self.pin(f"A2 {name}", loss, mcse, HARM_BOUND - loss, today, meets,
                 HARM_BOUND / 4)

    def test_null(self):
        self.check("null", 0.0, True)

    def test_overfit(self):
        self.check("overfit", -0.0169, True)

    def test_clear_winner(self):
        self.check("clear_winner", -0.0704, True)

    def test_unrelated_proxy(self):
        # An uninformative proxy makes the pick a near-uniform draw from the
        # window, which here lies on the plateau.
        self.check("unrelated_proxy", -0.0006, True)

    def test_default_best_misses(self):
        """Misses by 0.0029: the default is never protected. The
        fewest-complexity candidate wins whatever the window admits, and the
        default is preferred only after complexity."""
        self.check("default_best", 0.0079, False)

    def test_underfit_misses(self):
        """Misses by 0.0019: when the true score rises with complexity, the
        fewest-complexity candidate is the worst one in the window."""
        self.check("underfit", 0.0069, False)

    def test_shared_folds_misses(self):
        """Misses by 0.0033: the window uses one trial's SE, which includes
        the fold effect every configuration shares; that effect cancels in a
        paired comparison but widens the window."""
        self.check("shared_folds", 0.0083, False)


class A3Power(Pinned):
    """A3: with one configuration WIN_GAP better than the default, the rule
    recovers at least WIN_SHARE of the gap in expectation. Today's rule
    meets A3."""

    def test_clear_winner(self):
        o = outcome(scenario("clear_winner"))
        gain, mcse = mean_and_mcse(o["one_se"] - o["default"])
        bound = WIN_SHARE * WIN_GAP
        self.pin("A3 clear_winner", gain, mcse, gain - bound, 0.0704, True,
                 bound / 10)


class A4PreferenceDirection(Pinned):
    """A4: among the candidates the rule treats as indistinguishable, its
    preference picks one no worse in expectation than a uniform pick from
    the same candidates (within PREF_BOUND), on every suite scenario. A
    preference whose ordering runs against the held-out score misses.
    Today's rule misses on three scenarios."""

    def check(self, name, today, meets):
        o = outcome(scenario(name))
        cost, mcse = mean_and_mcse(o["uniform"] - o["one_se"])
        self.pin(f"A4 {name}", cost, mcse, PREF_BOUND - cost, today, meets,
                 PREF_BOUND / 2)

    def test_null(self):
        self.check("null", 0.0, True)

    def test_overfit(self):
        self.check("overfit", -0.0002, True)

    def test_default_best(self):
        self.check("default_best", -0.0048, True)

    def test_unrelated_proxy(self):
        self.check("unrelated_proxy", 0.0007, True)

    def test_underfit_misses(self):
        """Misses by 0.0062: the true score rises with complexity, so the
        lowest-complexity candidate is the worst one in the window."""
        self.check("underfit", 0.0082, False)

    def test_shared_folds_misses(self):
        """Misses by 0.0083: as underfit, in a window widened by the shared
        fold effect."""
        self.check("shared_folds", 0.0103, False)

    def test_clear_winner_misses(self):
        """Misses by 0.0019: the better configuration is the more complex
        one, so whenever noise lifts a simpler rival into the window, the
        preference takes the rival."""
        self.check("clear_winner", 0.0039, False)


class A5Calibration(Pinned):
    """A5: a score reported for the selected configuration is within
    CAL_BOUND of its true held-out score in expectation, on every suite
    scenario. select.py reports the selected trial's inner mean, which
    misses wherever near-ties make the selection depend on noise."""

    def check(self, name, today, meets):
        o = outcome(scenario(name))
        bias, mcse = mean_and_mcse(o["reported"] - o["one_se"])
        self.pin(f"A5 {name}", bias, mcse, CAL_BOUND - abs(bias), today, meets,
                 CAL_BOUND / 4)

    def test_null_misses(self):
        """Misses by 0.0235: the reported mean is the raw inner mean of a
        trial chosen for a high inner mean (the winner's curse); with 32
        equal configurations it is biased upward by about two noise SDs."""
        self.check("null", 0.0335, False)

    def test_null_binomial_misses(self):
        """Misses by 0.0307, as test_null_misses, under binomial fold noise."""
        self.check("null_binomial", 0.0407, False)

    def test_default_best_misses(self):
        """Misses by 0.0129: near-ties at the top, so the selected trial's
        inner mean carries the winner's curse."""
        self.check("default_best", 0.0229, False)

    def test_clear_winner(self):
        # One configuration leads by about 2.6 paired SEs, so the selection
        # barely depends on noise and the reported mean is nearly unbiased.
        self.check("clear_winner", 0.0073, True)


class Harness(unittest.TestCase):
    """The toy problems are what the doc says they are."""

    def test_draws_are_fixed_by_the_seed(self):
        s = scenario("null")
        self.assertTrue(np.array_equal(draw(s, reps=3), draw(s, reps=3)))

    def test_default_is_configuration_zero_and_the_truth_is_known(self):
        for name in SUITE + ("null_binomial",):
            s = scenario(name)
            self.assertEqual(s["mu"].shape, (K,))
            self.assertEqual(s["kappa"].shape, (K,))
            self.assertTrue(np.all(s["kappa"] > 0))
        self.assertEqual(int(np.argmax(scenario("clear_winner")["mu"])), 20)
        self.assertEqual(int(np.argmax(scenario("default_best")["mu"])), 0)

    def test_task_b_plan_has_k_distinct_trials_with_the_default_first(self):
        configs = {hp.canonical_json(t["config"]) for t in TASK_B_PLAN}
        self.assertEqual(len(configs), K)
        self.assertEqual([t["trial_id"] for t in TASK_B_PLAN], list(range(K)))
        self.assertTrue(TASK_B_PLAN[0]["is_default"])

    def test_the_fold_plan_hides_the_outer_slots(self):
        self.assertEqual(hp.inner_fold_ids(FOLD_PLAN), [0, 1, 2])
        self.assertEqual(FOLD_PLAN["inner_folds"][-3:], [None] * 3)


if __name__ == "__main__":
    unittest.main()
