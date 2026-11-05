"""Toy-data acceptance criteria for a learner-hyperparameter selection rule
(docs/hyperopt-selection.md), applied to three rules:

- select.py's rule, the fixed default (the A* classes): DECISIONS F4's
  fallback, adopted on 5 Nov 2026 because the candidate below failed A6;
- the default-anchored rule, the candidate F4 evaluated (the Candidate*
  classes): keep the family default unless the best inner mean beats the
  default's by more than ANCHOR_GAMMA (0.04). It meets A1 to A5 here; it
  lost more than 0.005 to the default on the lineage split for scm in A6,
  which is not a toy test;
- the one-standard-error window with the lower-complexity preference, the
  rule select.py held before F4 (the Baseline* classes).

Scope: learner hyperparameters only, the per-family knobs in protocol.py.
Nothing here covers choosing the input-data configuration (cohort and genome
count, which genomes are included, k-mer length, feature families or the
retained-feature budget), which is a separate open capability.

Generic toy data only. A "configuration" is an index with a known true
held-out score; fold scores are that score plus Gaussian (or binomial) noise,
with a fold effect shared by every configuration. There are no genomes,
learners, features or real data here, and nothing is fitted. Every draw comes
from a fixed numpy seed, so results do not depend on PYTHONHASHSEED.

The Monte Carlo runs a vectorised copy of each rule. RuleCopy checks the
fixed-default copy (`fixed_default`) against select.select_history on drawn
histories, flat and Task B, so a change to select.py's rule fails that test
first. The candidate's and the baseline's copies (`default_anchored`,
`one_se`) are the evaluation record; select.py holds neither.

The bounds were ratified as written (DECISIONS F4, 5 Nov 2026). Each
criterion test pins the rule's value on its scenario and whether that value
meets the bound, so a miss would be asserted as a miss, by its size. A rule
that changes a value fails the pin; update the pin with the evidence the
doc's A6 asks for. A5 is met by the score a rule reports for its
selection: the selected configuration's score on an outer fold the
selection never saw, drawn here as fresh noise around its true score. Its
inner mean is reported only labelled as a selection-biased tuning score.
The fixed default has no power, so it misses A3 by design; the miss is
asserted.

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

# The candidate rule's margin, and the relative slack on the rules'
# boundaries, so a mean lying exactly on one in real arithmetic is not
# dropped by rounding (select.py's BOUNDARY_RTOL before F4).
ANCHOR_GAMMA = 0.04
BOUNDARY_RTOL = 1e-12

# A5's pins: the mean error of the outer-fold score the rule reports, per
# scenario. Each is within two Monte Carlo SEs (about 0.0007) of zero, as an
# independent outer fold must be.
OUTER_BIAS = {"null": 0.00073, "null_binomial": -0.00016,
              "default_best": -0.00009, "underfit": -0.00014,
              "overfit": 0.00128, "unrelated_proxy": 0.00045,
              "shared_folds": 0.00097, "clear_winner": -0.00078}


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


def outer_noise(s, reps=REPS, seed=SEED):
    """One untouched outer fold per draw: the noise of a score taken on
    samples no inner fold held, independent of the inner draws, with the
    per-fold SD. Shape (reps,)."""
    rng = np.random.default_rng([seed, SUITE.index(s["name"])
                                 if s["name"] in SUITE else 99, 1])
    return rng.normal(0.0, SD, reps)


# --------------------------------------------------------------------- rules

def preference_rank(kappa):
    """select.preference_key as a rank: lower complexity, then the default
    (configuration 0), then the lowest id."""
    order = sorted(range(len(kappa)), key=lambda c: (kappa[c], c != 0, c))
    rank = np.empty(len(kappa), dtype=np.int64)
    rank[order] = np.arange(len(kappa))
    return rank


def fixed_default(x):
    """select.py's rule, vectorised over draws: configuration 0 (the
    default) is selected whatever the scores, and is the only candidate.
    The best mean, ties to the lowest id (so to the default first), is
    reported for the record. Returns the selected and the best
    configuration per draw, the candidate mask and the means."""
    m = x.mean(axis=2)
    cand = np.zeros(m.shape, dtype=bool)
    cand[:, 0] = True
    return np.zeros(len(x), dtype=np.int64), m.argmax(axis=1), cand, m


def default_anchored(x):
    """The candidate rule, vectorised over draws: keep configuration 0 (the
    default) unless the best mean beats its mean by more than GAMMA. The
    best is the highest mean, ties to the lowest id, so to the default
    first. The candidates are the configurations within GAMMA of the best
    mean. Returns the selected and the best configuration per draw, the
    candidate mask and the means."""
    m = x.mean(axis=2)
    rows = np.arange(len(x))
    best = m.argmax(axis=1)
    threshold = m[rows, best] - ANCHOR_GAMMA
    slack = BOUNDARY_RTOL * np.maximum(1.0, np.abs(threshold))
    cand = m >= (threshold - slack)[:, None]
    chosen = np.where(cand[:, 0], 0, best)
    return chosen, best, cand, m


def one_se(x, kappa):
    """The replaced rule (select.py before F4), vectorised over draws, kept
    as the baseline. Returns the selected and the best-mean configuration
    per draw, the candidate mask and the means."""
    reps, k, n = x.shape
    m = x.mean(axis=2)
    se = x.std(axis=2, ddof=1) / math.sqrt(n)
    rank = preference_rank(kappa)
    rows = np.arange(reps)
    best = np.where(m == m.max(axis=1, keepdims=True), rank, k).argmin(axis=1)
    threshold = m[rows, best] - se[rows, best]
    slack = BOUNDARY_RTOL * np.maximum(1.0, np.abs(threshold))
    cand = m >= (threshold - slack)[:, None]
    chosen = np.where(cand, rank, k).argmin(axis=1)
    return chosen, best, cand, m


def outcome(s, x=None):
    """Held-out results on one scenario, per draw, for the adopted rule
    (`fixed`), the candidate (`anchored`) and the baseline (`one_se`): the
    true score of each
    selection and of the default; the expected true
    score of a uniform pick from each rule's candidate set; the score each
    rule reports (the adopted rule's outer-fold score, the baseline's inner
    mean); the adopted rule's labelled inner mean; and the candidate
    counts."""
    x = draw(s) if x is None else x
    rows = np.arange(len(x))
    mu = s["mu"]
    out = {"default": np.full(len(x), mu[0])}
    chosen, _, cand, m = fixed_default(x)
    out.update(fixed=mu[chosen],
               fixed_uniform=(cand * mu).sum(axis=1) / cand.sum(axis=1),
               fixed_reported=mu[chosen] + outer_noise(s, len(x)))
    chosen, _, cand, m = default_anchored(x)
    out.update(anchored=mu[chosen],
               anchored_uniform=(cand * mu).sum(axis=1) / cand.sum(axis=1),
               anchored_reported=mu[chosen] + outer_noise(s, len(x)),
               anchored_inner=m[rows, chosen],
               anchored_n_cand=cand.sum(axis=1))
    chosen, _, cand, m = one_se(x, s["kappa"])
    out.update(one_se=mu[chosen],
               uniform=(cand * mu).sum(axis=1) / cand.sum(axis=1),
               reported=m[rows, chosen], n_cand=cand.sum(axis=1))
    return out


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
        chosen, best, cand, m = fixed_default(x)
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
    plan it pins; it is deterministic. The fixed default meets A1."""

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

    def test_every_selection_is_the_default(self):
        s = scenario("clear_winner")
        for make, protocol, family in ((history, PROTOCOL, "toy"),
                                       (task_b_history, TASK_B, "generic")):
            h = make(draw(s, reps=1, seed=SEED + 8)[0], s["kappa"])
            r = sel.select_history(h, protocol, FOLD_PLAN)
            with self.subTest(family=family):
                f = r["families"][family]
                self.assertEqual((f["selected"]["trial_id"], f["best"]["trial_id"]),
                                 (0, 20))
                self.assertIn("fixed default (DECISIONS F4 fallback", r["rule"])


class Pinned(unittest.TestCase):
    """Shared check for A2 to A5: the rule's value on one scenario, and
    whether it meets the criterion. `margin` is how far the value lies
    inside the bound; a negative margin is the size of the miss."""

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
    """A2 for the fixed default: it is the default, so its loss against the
    default is zero on every suite scenario."""

    def test_every_scenario(self):
        for name in SUITE:
            o = outcome(scenario(name))
            loss, mcse = mean_and_mcse(o["default"] - o["fixed"])
            with self.subTest(scenario=name):
                self.pin(f"A2 {name}", loss, mcse, HARM_BOUND - loss, 0.0, True,
                         HARM_BOUND / 4)


class A3Power(Pinned):
    """A3 for the fixed default: it never moves off the default, so it
    recovers none of WIN_GAP. It misses A3 by design; that is the price of
    the fallback DECISIONS F4 set for a candidate that fails A6."""

    def test_clear_winner_misses(self):
        o = outcome(scenario("clear_winner"))
        gain, mcse = mean_and_mcse(o["fixed"] - o["default"])
        bound = WIN_SHARE * WIN_GAP
        self.pin("A3 clear_winner", gain, mcse, gain - bound, 0.0, False,
                 bound / 10)


class A4PreferenceDirection(Pinned):
    """A4 for the fixed default: its one candidate is its pick, so its cost
    against a uniform pick from the candidates is zero."""

    def test_every_scenario(self):
        for name in SUITE:
            o = outcome(scenario(name))
            cost, mcse = mean_and_mcse(o["fixed_uniform"] - o["fixed"])
            with self.subTest(scenario=name):
                self.pin(f"A4 {name}", cost, mcse, PREF_BOUND - cost, 0.0, True,
                         PREF_BOUND / 2)


class A5Calibration(Pinned):
    """A5 for the fixed default: the score it reports is the default's
    outer-fold score, unbiased up to Monte Carlo error. The outer noise is
    the candidate's (OUTER_BIAS), since both report the same kind of score
    from the same draws."""

    def test_every_scenario(self):
        for name in SUITE + ("null_binomial",):
            o = outcome(scenario(name))
            bias, mcse = mean_and_mcse(o["fixed_reported"] - o["fixed"])
            with self.subTest(scenario=name):
                self.pin(f"A5 {name}", bias, mcse, CAL_BOUND - abs(bias),
                         OUTER_BIAS[name], True, CAL_BOUND / 4)


class CandidateA2NoHarmAgainstDefault(Pinned):
    """A2: on every suite scenario, the expected held-out loss of the selected
    configuration against the fixed default, E[mu_default - mu_selected], is
    at most HARM_BOUND. The candidate meets it everywhere."""

    def check(self, name, value, meets=True):
        o = outcome(scenario(name))
        loss, mcse = mean_and_mcse(o["default"] - o["anchored"])
        self.pin(f"A2 {name}", loss, mcse, HARM_BOUND - loss, value, meets,
                 HARM_BOUND / 4)

    def test_null(self):
        self.check("null", 0.0)

    def test_default_best(self):
        # The closest call: a rival within GAMMA of a default 0.003 better
        # never displaces it, but noise sometimes lifts one past GAMMA.
        self.check("default_best", 0.00243)

    def test_underfit(self):
        self.check("underfit", -0.00130)

    def test_overfit(self):
        self.check("overfit", -0.01440)

    def test_unrelated_proxy(self):
        self.check("unrelated_proxy", -0.00137)

    def test_shared_folds(self):
        # The fold effect cancels in the comparison with the default, which
        # is paired: every configuration is scored on the same folds.
        self.check("shared_folds", -0.00044)

    def test_clear_winner(self):
        self.check("clear_winner", -0.07214)


class CandidateA3Power(Pinned):
    """A3: with one configuration WIN_GAP better than the default, the rule
    recovers at least WIN_SHARE of the gap in expectation. The
    default-anchored rule meets A3: the winner leads the default by about
    2.6 paired SEs, so it clears GAMMA in about nine draws in ten."""

    def test_clear_winner(self):
        o = outcome(scenario("clear_winner"))
        gain, mcse = mean_and_mcse(o["anchored"] - o["default"])
        bound = WIN_SHARE * WIN_GAP
        self.pin("A3 clear_winner", gain, mcse, gain - bound, 0.07214, True,
                 bound / 10)


class CandidateA4PreferenceDirection(Pinned):
    """A4: among the candidates the rule treats as indistinguishable, the
    configurations within GAMMA of the best mean, its pick (the default when
    it is one of them, else the best) is no worse in expectation than a
    uniform pick from them (within PREF_BOUND), on every suite scenario.
    The candidate meets it everywhere."""

    def check(self, name, value, meets=True):
        o = outcome(scenario(name))
        cost, mcse = mean_and_mcse(o["anchored_uniform"] - o["anchored"])
        self.pin(f"A4 {name}", cost, mcse, PREF_BOUND - cost, value, meets,
                 PREF_BOUND / 2)

    def test_null(self):
        self.check("null", 0.0)

    def test_default_best(self):
        self.check("default_best", -0.01418)

    def test_underfit(self):
        self.check("underfit", -0.00136)

    def test_overfit(self):
        # The closest call: the default (capacity 30) sits past the optimum,
        # so a window that keeps it holds better configurations on average.
        self.check("overfit", 0.00011)

    def test_unrelated_proxy(self):
        self.check("unrelated_proxy", -0.00134)

    def test_shared_folds(self):
        self.check("shared_folds", -0.00094)

    def test_clear_winner(self):
        self.check("clear_winner", -0.00511)


class CandidateA5Calibration(Pinned):
    """A5: the score reported for the selected configuration is within
    CAL_BOUND of its true held-out score in expectation, on every suite
    scenario. The default-anchored rule reports the selected configuration's
    score on an outer fold the selection never saw, so it meets A5. Its
    inner mean, biased by as much as the baseline's, is shown only labelled
    as a selection-biased tuning score (test_the_inner_mean_is_labelled)."""

    def check(self, name, value, meets=True):
        o = outcome(scenario(name))
        bias, mcse = mean_and_mcse(o["anchored_reported"] - o["anchored"])
        self.pin(f"A5 {name}", bias, mcse, CAL_BOUND - abs(bias), value, meets,
                 CAL_BOUND / 4)

    def test_null(self):
        self.check("null", OUTER_BIAS["null"])

    def test_null_binomial(self):
        self.check("null_binomial", OUTER_BIAS["null_binomial"])

    def test_default_best(self):
        self.check("default_best", OUTER_BIAS["default_best"])

    def test_underfit(self):
        self.check("underfit", OUTER_BIAS["underfit"])

    def test_overfit(self):
        self.check("overfit", OUTER_BIAS["overfit"])

    def test_unrelated_proxy(self):
        self.check("unrelated_proxy", OUTER_BIAS["unrelated_proxy"])

    def test_shared_folds(self):
        self.check("shared_folds", OUTER_BIAS["shared_folds"])

    def test_clear_winner(self):
        self.check("clear_winner", OUTER_BIAS["clear_winner"])

    def test_the_inner_mean_is_labelled(self):
        # Pinned for the record, not gated: the inner mean the rule selected
        # on overstates the selection's held-out score by up to 0.04.
        o = outcome(scenario("null"))
        bias, mcse = mean_and_mcse(o["anchored_inner"] - o["anchored"])
        self.assertAlmostEqual(bias, 0.03501, delta=max(PIN_SES * mcse, 1e-4))
        s = scenario("null")
        h = history(draw(s, reps=1, seed=SEED + 9)[0], s["kappa"])
        chosen = select_toy(h)["toy"]["selected"]
        self.assertEqual(chosen["score_label"], sel.SELECTED_SCORE_LABEL)
        self.assertIn("selection-biased tuning score", chosen["score_label"])
        self.assertIn("selection-biased tuning score", sel.SCORE_LABEL)


class BaselineA2NoHarmAgainstDefault(Pinned):
    """The replaced one-SE rule's A2, recorded: it misses on three
    scenarios."""

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


class BaselineA3Power(Pinned):
    """A3: with one configuration WIN_GAP better than the default, the rule
    recovers at least WIN_SHARE of the gap in expectation. The baseline
    meets A3."""

    def test_clear_winner(self):
        o = outcome(scenario("clear_winner"))
        gain, mcse = mean_and_mcse(o["one_se"] - o["default"])
        bound = WIN_SHARE * WIN_GAP
        self.pin("A3 clear_winner", gain, mcse, gain - bound, 0.0704, True,
                 bound / 10)


class BaselineA4PreferenceDirection(Pinned):
    """A4: among the candidates the rule treats as indistinguishable, its
    preference picks one no worse in expectation than a uniform pick from
    the same candidates (within PREF_BOUND), on every suite scenario. A
    preference whose ordering runs against the held-out score misses.
    The baseline misses on three scenarios."""

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


class BaselineA5Calibration(Pinned):
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
