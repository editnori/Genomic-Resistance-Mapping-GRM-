"""Paired statistics for feature-block ablation contrasts (GRM ablation contract v1.9).

A contrast compares a baseline arm and a treatment arm of one (dataset,
learner, split) on the same samples, fold plans, labels and lineage clusters.
Each sample's score is its mean out-of-fold probability across repeats.

Pairing rules (PAIRINGS; DECISIONS V-1 for the input comparison), named in
every contrast record's `pairing`:
    same_plan              the same samples and the same fold plan in every repeat
    k31_groups             the same, for a treatment counted at another k whose
                           fold plan was drawn on the baseline's k = 31 lineages and
                           duplicate groups, so it is the baseline's plan
    common_outer_test_set_balanced
                           a treatment with another cohort (cap or selection rule),
                           scored on a class-balanced subsample of the baseline's
                           outer test folds restricted to the genomes both cohorts
                           hold (the subsample matches each baseline fold's class
                           balance; `grm arm` draws it): its samples are a subset of
                           the baseline's, each with the baseline's fold, label and
                           cluster under the baseline's fold plan, and the baseline
                           is restricted to them. Both arms trained on other samples
                           too, so a no-call row's score, its own training
                           prevalence, is not recomputed from the paired rows.
    common_outer_test_set  the same pairing without the balanced subsample, kept so
                           records written before the rulings after lane INPUT still
                           validate; the input-comparison runner no longer uses it.
Any other pairing, such as two arms with different samples under same_plan, is
refused with PairingError.
Every estimate is oriented so that a positive value favours the treatment:

    log_loss           LL_base - LL_treat, probabilities clipped to [1e-4, 1 - 1e-4]
    r2_ll              (LL_base - LL_treat) / LL_null; the null scores each sample
                       with the resistant prevalence of its training folds
    auroc, auprc,
    balanced_accuracy  treat - base
    vme, me            base - treat; VME is false S among reference R and ME is
                       false R among reference S

Thresholded endpoints (balanced_accuracy, vme, me) use the rows' `called`, not
the score: a sample is an R call when `called` is true in at least half of its
repeats, so a no-call row is never a resistant call, whatever its prevalence
score (contract section 13).

An endpoint missing a class it needs is undefined, never 0 or 0.5. log_loss
needs no class; vme needs resistant samples, me susceptible samples, and the
rest need both.

Families (contract sections 9 and 10):
    A, B  pooled r2_ll records only, from pool(); Holm within each family
    C     per-dataset log_loss contrasts of the A and B comparisons, from compare(); BH
    D     exploratory: per-dataset contrasts on any endpoint and exploratory pools; no p
Per-dataset r2_ll records, the inputs of pool(), are family D. Development
datasets are refused in A, B and C.

Inference per dataset:
- p (family C only): within-cluster arm swap. The delta of every endpoint but
  auroc and auprc is a sum of per-sample terms whose weights depend only on y,
  so the swap is a sign flip of cluster sums: exact over all 2^k patterns at
  k <= 20 clusters, else 20,000 seeded flips with p = (1 + hits) / (1 + n).
  arm_swap_p() recomputes auroc and auprc per swap; compare() does not need it,
  since those endpoints are family D.
- ci and ci_tost: cluster jackknife-t (ci_method "cluster_jackknife_t"). Each
  cluster is left out in turn and the endpoint recomputed with both arms and
  the null; SE^2 = (k-1)/k * sum((theta_(i) - mean)^2), ci = estimate +/-
  t(0.975, k-1) SE (95%) and ci_tost = estimate +/- t(0.95, k-1) SE (90%), so
  ci contains ci_tost.
- A contrast is undefined when the endpoint misses a class, when there is one
  cluster, or when any leave-one-cluster-out estimate is undefined
  (undefined_resamples counts those). The jackknife variance needs all k
  replicates, so this is stricter than the 10% rule of contract section 9;
  above 10% the reason is that section's text. An undefined contrast has null
  estimate, ci, ci_tost, p and p_adj, ci_method "none", its reason in
  undefined_reason, and decision inconclusive.

Pooled (pool(), r2_ll only): datasets are averaged within species, then
across species. p is a sign flip of the species means; ci and ci_tost are the
species jackknife-t, which for this mean is the one-sample t interval on the
species means. With fewer than two contributing species there is no
species-level interval, so those datasets are excluded with that reason and
the record is null and inconclusive. The exact sign flip cannot give a p below
2 / 2^S with S species, so a pooled improves or worse needs at least 6 species.
The pooled p tests the null "no effect in any species" (each species' estimate
error symmetric about zero), not "the average species effect is zero" (contract
section 15). dataset_summary() gives the median
per-dataset delta and where each dataset's ci sits against the margin.

Decisions (decide(), contract section 15):
    improves      ci[0] > 0 and p_adj < 0.05
    worse         ci[1] < 0 and p_adj < 0.05
    inconclusive  otherwise
In family D (no p) the decision is always inconclusive. There is no equivalent
decision: the J re-check (grm-recheck-task-j-20261030) measured P(90% interval
inside +/-margin | true delta = margin) at 0.059 with 60 and 80 clusters, and
0.08 to 0.11 with skewed cluster effects up to 240 clusters, so contract v1.8
withdrew it. ci_tost keeps its meaning; K derives a
descriptive within-margin note from it.

finalize(records, planned_ids, missing) adjusts each family on its own and
refuses a call that lacks a planned id: every planned id is a record or is in
`missing` with a reason. Undefined records and missing ids of a family enter
its adjustment with p = 1, so dropping a hypothesis never weakens the
correction for the others. Its `unadjusted` ids (the input comparison's
sentinels) are family C but judged on their own p, outside the BH family:
p_adj = p and adjust "none". finalize_sets() takes several (records, planned,
missing) sets, such as one per species, and adjusts them as one family.

The input comparison's test (main's rulings after VERIFY-H1STATS), from
pair_units() and compare_units(): a cluster test treats every lineage of one
outer fold as independent although one pair of fitted models scores them all,
so its unit is the (repeat, outer fold) refit instead. Each arm scores R keyed
lineage repeats; per repeat J pairs the arms under the rule above (a cohort
arm's balanced subsample may differ by repeat), and each unit's value is the
mean log-loss delta over its paired test samples. With R = 3 and 5 folds the
15 unit means give a t by learner (LEARNER_TESTS, main's ruling after round
2's calibration): scm, cart and rf the paired t over the unit means (ci_method
"repeat_fold_t", SE = s / sqrt(J)); gbt and the lasso the corrected resampled
t of Nadeau and Bengio (main's ruling after the first calibration failed): SE = sqrt((1/J + rho) s^2), where s^2
is the units' sample variance and rho the mean over units of n_test / n_train
(the genomes paired in the unit over the baseline's training genomes in that
repeat and fold), df J - 1 = 14; p two-sided from Student's t, ci and ci_tost
the 95% and 90% t intervals (ci_method "repeat_fold_corrected_t", unit
"repeat_fold"). The units share genomes across repeats and folds, so the plain
1/J variance is too small (the calibration of 7 Nov 2026 measured gbt at 8.4%).
The AUROC delta is the mean over units of each unit's AUROC delta, with its 95%
interval by the same SE. Under common_outer_test_set_balanced, J also
checks the subsample is the balanced draw of the baseline's whole test fold:
every baseline fold of every repeat is scored, and each fold's class counts
are what `grm arm`'s balanced_counts gives for the baseline fold's.

Contrast record, keys in RECORD_KEYS order (contract section 1):

    id                   "<dataset>/<contrast key>"; the key is
                         "<family>:<learner>:<split>:<baseline>-><treatment>:<endpoint>"
                         and holds no whitespace, slash or backtick
    family               "C" | "D"
    baseline, treatment  arm names
    learner, split       as in the input rows
    endpoint             one of ENDPOINTS; "log_loss" in family C
    estimate             float, or None when undefined
    ci, ci_tost          [lo, hi] 95% and 90% intervals by ci_method's test, or None
    ci_method            one of CI_METHODS: "cluster_jackknife_t" (compare()), "repeat_fold_t" or
                         "repeat_fold_corrected_t" (compare_units(), by LEARNER_TESTS in the input
                         comparison), or "none" when undefined
    unit, n_units        "cluster", clusters; or "repeat_fold", (repeat, outer fold) units
    few_clusters         n_units < 20 on unit "cluster", false on "repeat_fold" (a flag for readers)
    p                    two-sided randomisation p (compare()) or t p (compare_units()) in family C;
                         None in D or when undefined
    p_adj                BH-adjusted p in family C, never below p; p itself on an unadjusted (sentinel)
                         id; None when p is None
    adjust               "bh" | "none" (family D, or an unadjusted family C id)
    margin               predeclared margin, finite and > 0; ci_tost is read against it
    decision             "improves" | "worse" | "inconclusive"
    claim_level          "predictive_increment"
    undefined_resamples  leave-one-cluster-out replicates where the endpoint is undefined
    undefined_reason     None, or why the contrast is undefined
    interval_label       INTERVAL_LABEL
    dataset              the id prefix
    pairing              the pairing rule pair_arms() checked, one of PAIRINGS
    auroc_delta, auroc_ci
                         on a log_loss record, treat - base AUROC and its 95% interval by
                         the record's ci_method test, descriptive only (no p, no decision); None
                         when undefined or on another endpoint. A log-loss improves can come
                         from calibration or overfitting; the AUROC delta shows whether
                         discrimination moved with it

Pooled record, keys in POOLED_KEYS order (contract section 2):

    contrast             the contrast key, endpoint r2_ll
    family               "A" | "B" | "D"
    unit, n_units        "species", contributing species
    n_datasets           every planned dataset, contributing or excluded
    n_excluded, excluded [{dataset, reason}] for each dataset that does not contribute
    estimate, ci, ci_tost, p, p_adj, adjust, margin, decision, interval_label
                         as above; null values and inconclusive when nothing contributes

compare() and pool() return pending records: these keys plus "_pending",
with p_adj and decision unset. finalize() returns them with exactly the keys
above. This module reads no files and makes no claims beyond the cohort it
is given.
"""
from __future__ import annotations

import math
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache

import numpy as np

EPS = 1e-4
ALPHA = 0.05
CI_LEVEL = 0.95
TOST_LEVEL = 0.90
FEW_CLUSTERS = 20
MAX_UNDEFINED_SHARE = 0.10
UNDEFINED_SHARE_REASON = "more than 10% of resamples undefined"
EXACT_MAX = 20
INTERVAL_LABEL = ("conditional on this cohort; CV intervals ignore training-set refit "
                  "variability; not external validation")
CLAIM_LEVEL = "predictive_increment"
FAMILIES = {"A": "holm", "B": "holm", "C": "bh", "D": "none"}
CONTRAST_FAMILIES = ("C", "D")
POOLED_FAMILIES = ("A", "B", "D")
CONFIRMATORY = ("A", "B", "C")
PRIMARY = "log_loss"
POOLED_ENDPOINT = "r2_ll"
JACKKNIFE = "cluster_jackknife_t"
UNIT_T = "repeat_fold_t"
UNIT_CORRECTED_T = "repeat_fold_corrected_t"
UNIT = "repeat_fold"
CI_METHODS = (JACKKNIFE, UNIT_T, UNIT_CORRECTED_T, "none")
# The input comparison's test per learner (main's ruling after round 2's calibration, 7 Nov 2026, pre-registered in
# docs/input-comparison-runner.md): the fold-mean t held size for scm, cart and rf with full power; the corrected t
# holds it for gbt and the lasso, whose fold-mean rows failed. The choice rests only on synthetic null calibration.
LEARNER_TESTS = {"scm": UNIT_T, "cart": UNIT_T, "rf": UNIT_T, "gbt": UNIT_CORRECTED_T, "lasso": UNIT_CORRECTED_T}
UNIT_INTERVAL_LABELS = {
    UNIT_T: ("conditional on this cohort; each unit is one (repeat, outer fold) refit, and the interval is the paired t "
             "over the unit means, which carries refit variability across units; not external validation"),
    UNIT_CORRECTED_T: ("conditional on this cohort; each unit is one (repeat, outer fold) refit, and the interval is "
                       "Nadeau and Bengio's corrected resampled t, which carries refit variability and the overlap of "
                       "training sets across units; not external validation")}
FIELDS = ("dataset", "arm", "learner", "split", "repeat", "sample_id", "y", "fold", "score", "called",
          "cluster", "fold_plan_sha256")
# Classes each endpoint needs: 1 resistant, 0 susceptible.
NEEDS = {"log_loss": (), "r2_ll": (0, 1), "auroc": (0, 1), "balanced_accuracy": (0, 1),
         "vme": (1,), "me": (0,), "auprc": (0, 1)}
ENDPOINTS = tuple(NEEDS)
CLASS_NAMES = {0: "susceptible (y=0)", 1: "resistant (y=1)"}
DECISIONS = ("improves", "worse", "inconclusive")
RECORD_KEYS = ("id", "family", "baseline", "treatment", "learner", "split", "endpoint", "estimate", "ci",
               "ci_tost", "ci_method", "unit", "n_units", "few_clusters", "p", "p_adj", "adjust", "margin",
               "decision", "claim_level", "undefined_resamples", "undefined_reason", "interval_label", "dataset",
               "pairing", "auroc_delta", "auroc_ci")
PAIRINGS = ("same_plan", "k31_groups", "common_outer_test_set_balanced", "common_outer_test_set")
SUBSET_PAIRINGS = ("common_outer_test_set_balanced", "common_outer_test_set")
POOLED_KEYS = ("contrast", "family", "unit", "n_units", "n_datasets", "n_excluded", "excluded", "estimate", "ci",
               "ci_tost", "p", "p_adj", "adjust", "margin", "decision", "interval_label")
DATASET_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9._+-]*")
CELLS = 1 << 22  # sign-pattern matrix entries per chunk
CALL_TOLERANCE = 1e-9  # a no-call score must equal its training prevalence to this


class PairingError(ValueError):
    """The two arms are not a paired comparison; `check` names the failed test."""

    def __init__(self, check, message):
        super().__init__(f"{check}: {message}")
        self.check = check


class DevelopmentDatasetError(ValueError):
    """A development dataset was offered to family A, B or C."""


@dataclass(frozen=True, eq=False)
class PairedArms:
    dataset: str
    learner: str
    split: str
    baseline: str
    treatment: str
    sample_ids: tuple
    y: np.ndarray              # 0/1 per sample
    cluster_index: np.ndarray  # position in clusters, per sample
    clusters: tuple
    base: np.ndarray           # mean out-of-fold probability across repeats
    treat: np.ndarray
    null: np.ndarray           # mean training-fold prevalence across repeats
    n_repeats: int
    base_call: np.ndarray = None  # R call: called in at least half the repeats; pair_arms sets it
    treat_call: np.ndarray = None
    pairing: str = "same_plan"    # the rule pair_arms checked, one of PAIRINGS

    def __post_init__(self):
        # Arrays built without rows hold only model scores, where called is score >= 0.5.
        for arm in ("base", "treat"):
            if getattr(self, f"{arm}_call") is None:
                object.__setattr__(self, f"{arm}_call", np.asarray(getattr(self, arm)) >= 0.5)


def _real(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)


def _layout(arm, table):
    if not table:
        raise PairingError("samples", f"arm {arm!r} has no rows")
    plans, samples = defaultdict(set), defaultdict(set)
    for (repeat, sample), row in table.items():
        plans[repeat].add(row["fold_plan_sha256"])
        samples[repeat].add(sample)
    for repeat, hashes in plans.items():
        if len(hashes) != 1:
            raise PairingError("fold_plan", f"arm {arm!r} repeat {repeat!r} has {len(hashes)} fold plans")
    covered = {frozenset(ids) for ids in samples.values()}
    if len(covered) != 1:
        raise PairingError("samples", f"arm {arm!r} covers different samples in different repeats")
    return {repeat: hashes.pop() for repeat, hashes in plans.items()}, covered.pop()


def _check_called(row, prevalence):
    """called = score >= 0.5, or false on a no-call row, whose score is its training prevalence.
    With prevalence None (SUBSET_PAIRINGS) a no-call row's score is not checked."""
    called, score = row["called"], float(row["score"])
    if type(called) is not bool:
        raise PairingError("called", f"sample {row['sample_id']!r} has called {called!r}; expected true or false")
    if prevalence is None:
        if called and score < 0.5:
            raise PairingError("called", f"sample {row['sample_id']!r} repeat {row['repeat']!r} is called with "
                                         f"score {score}")
        return
    if called != (score >= 0.5) and (called or abs(score - prevalence) > CALL_TOLERANCE):
        raise PairingError("called", f"sample {row['sample_id']!r} repeat {row['repeat']!r} has called={called} "
                                     f"with score {score}; a no-call row scores its training prevalence {prevalence}")


def _arm_tables(rows, dataset, learner, split, baseline, treatment, pairing):
    """{arm: {(repeat, sample_id): row}} for the two arms' rows of one (dataset, learner, split), with
    each row's fields, label and score checked; a repeated (repeat, sample_id) is refused."""
    if baseline == treatment:
        raise ValueError("baseline and treatment are the same arm")
    if pairing not in PAIRINGS:
        raise ValueError(f"unknown pairing {pairing!r}; expected one of {PAIRINGS}")
    tables = {baseline: {}, treatment: {}}
    for row in rows:
        missing = [field for field in FIELDS if field not in row]
        if missing:
            raise ValueError(f"row is missing {missing}")
        if (row["dataset"], row["learner"], row["split"]) != (dataset, learner, split) or row["arm"] not in tables:
            continue
        key = (row["repeat"], row["sample_id"])
        table = tables[row["arm"]]
        if key in table:
            raise PairingError("duplicate", f"arm {row['arm']!r} has two rows for repeat {key[0]!r}, sample {key[1]!r}")
        if type(row["y"]) is not int or row["y"] not in (0, 1):
            raise PairingError("y", f"sample {key[1]!r} has label {row['y']!r}; expected the integer 0 or 1")
        if not _real(row["score"]) or not 0.0 <= row["score"] <= 1.0:
            raise PairingError("score", f"sample {key[1]!r} has score {row['score']!r}; expected a number in [0, 1]")
        table[key] = row
    return tables


def pair_arms(rows, *, dataset, learner, split, baseline, treatment, pairing="same_plan"):
    """Validate the paired plan under `pairing` (one of PAIRINGS) and return per-sample
    arrays; raise PairingError otherwise."""
    tables = _arm_tables(rows, dataset, learner, split, baseline, treatment, pairing)
    (base_plans, base_samples), (treat_plans, treat_samples) = (_layout(arm, tables[arm]) for arm in (baseline, treatment))
    if pairing in SUBSET_PAIRINGS:
        if not treat_samples <= base_samples:
            raise PairingError("samples", f"{len(treat_samples - base_samples)} samples only in {treatment!r}; "
                                          f"{pairing} scores the treatment on the baseline's test genomes")
        tables[baseline] = {key: row for key, row in tables[baseline].items() if key[1] in treat_samples}
        base_samples = treat_samples
    if base_samples != treat_samples:
        raise PairingError("samples", f"{len(base_samples - treat_samples)} samples only in {baseline!r}, "
                                      f"{len(treat_samples - base_samples)} only in {treatment!r}")
    if base_plans != treat_plans:
        differ = sorted((r for r in base_plans.keys() | treat_plans.keys() if base_plans.get(r) != treat_plans.get(r)), key=str)
        raise PairingError("fold_plan", f"fold_plan_sha256 differs or is missing for repeats {differ}")

    repeats = sorted(base_plans, key=str)
    sample_ids = tuple(sorted(base_samples, key=str))
    labels, cluster_ids = [], []
    for sample in sample_ids:
        seen = [table[(repeat, sample)] for table in tables.values() for repeat in repeats]
        y = {row["y"] for row in seen}
        if len(y) != 1:
            raise PairingError("y", f"sample {sample!r} has labels {sorted(y)} across arms or repeats")
        cluster = {row["cluster"] for row in seen}
        if len(cluster) != 1 or None in cluster:
            raise PairingError("cluster", f"sample {sample!r} has clusters {sorted(map(str, cluster))}")
        for repeat in repeats:
            if tables[baseline][(repeat, sample)]["fold"] != tables[treatment][(repeat, sample)]["fold"]:
                raise PairingError("fold", f"sample {sample!r} repeat {repeat!r} is in different folds under one fold plan")
        labels.append(y.pop())
        cluster_ids.append(cluster.pop())

    y = np.array(labels, dtype=np.int64)
    base, treat, null = np.zeros(len(y)), np.zeros(len(y)), np.zeros(len(y))
    base_called, treat_called = np.zeros(len(y)), np.zeros(len(y))
    for repeat in repeats:
        folds = [tables[baseline][(repeat, sample)]["fold"] for sample in sample_ids]
        prevalence = np.empty(len(y))
        for fold in set(folds):
            held = np.array([f == fold for f in folds])
            if held.all():
                raise PairingError("fold_plan", f"repeat {repeat!r} has one fold, so no training samples")
            prevalence[held] = y[~held].mean()
        for table in tables.values():
            for i, sample in enumerate(sample_ids):
                _check_called(table[(repeat, sample)], None if pairing in SUBSET_PAIRINGS else prevalence[i])
        null += prevalence
        base += [float(tables[baseline][(repeat, sample)]["score"]) for sample in sample_ids]
        treat += [float(tables[treatment][(repeat, sample)]["score"]) for sample in sample_ids]
        base_called += [tables[baseline][(repeat, sample)]["called"] for sample in sample_ids]
        treat_called += [tables[treatment][(repeat, sample)]["called"] for sample in sample_ids]
    clusters = tuple(sorted(set(cluster_ids), key=str))
    position = {cluster: i for i, cluster in enumerate(clusters)}
    return PairedArms(dataset=dataset, learner=learner, split=split, baseline=baseline, treatment=treatment,
                      sample_ids=sample_ids, y=y, cluster_index=np.array([position[c] for c in cluster_ids]),
                      clusters=clusters, base=base / len(repeats), treat=treat / len(repeats),
                      null=null / len(repeats), n_repeats=len(repeats),
                      base_call=2 * base_called >= len(repeats), treat_call=2 * treat_called >= len(repeats),
                      pairing=pairing)


def sample_log_loss(y, p):
    p = np.clip(p, EPS, 1 - EPS)
    return -(y * np.log(p) + (1 - y) * np.log1p(-p))


def _midranks(x):
    order = np.argsort(x, kind="mergesort")
    ordered = x[order]
    starts = np.flatnonzero(np.r_[True, ordered[1:] != ordered[:-1]])
    ends = np.r_[starts[1:], len(x)]
    ranks = np.empty(len(x))
    ranks[order] = np.repeat((starts + ends + 1) / 2.0, ends - starts)
    return ranks


def auroc(y, s):
    """Mann-Whitney AUROC with ties counted one half; needs both classes."""
    resistant = y == 1
    n_r = resistant.sum()
    return float((_midranks(s)[resistant].sum() - n_r * (n_r + 1) / 2) / (n_r * (len(y) - n_r)))


def auprc(y, s):
    """Average precision over distinct score cuts, tied scores taken together."""
    order = np.argsort(-s, kind="mergesort")
    ordered_y, ordered_s = y[order], s[order]
    last = np.flatnonzero(np.r_[ordered_s[1:] != ordered_s[:-1], True])
    tp = np.cumsum(ordered_y)[last]
    recall = tp / ordered_y.sum()
    return float(np.sum(np.diff(recall, prepend=0.0) * tp / (last + 1)))


def undefined_reason(endpoint, y):
    n_r = int(np.sum(y))
    present = {1: n_r > 0, 0: n_r < len(y)}
    missing = [CLASS_NAMES[c] for c in NEEDS[endpoint] if not present[c]]
    return f"{endpoint} undefined: no {' or '.join(missing)} samples" if missing else None


def contributions(endpoint, y, base, treat, null, base_call=None, treat_call=None):
    """Per-sample terms that sum to the endpoint delta; None for auroc and auprc.

    The weights depend only on y, so swapping the arms inside a cluster negates
    that cluster's sum. Thresholded endpoints read the R calls (base_call,
    treat_call); without them a call is score >= 0.5, right only when no row
    is a no-call. Call only when the endpoint is defined.
    """
    resistant = y == 1
    n_r = resistant.sum()
    n_s = len(y) - n_r
    if endpoint in ("log_loss", "r2_ll"):
        terms = (sample_log_loss(y, base) - sample_log_loss(y, treat)) / len(y)
        return terms / sample_log_loss(y, null).mean() if endpoint == "r2_ll" else terms
    base_call = base >= 0.5 if base_call is None else base_call
    treat_call = treat >= 0.5 if treat_call is None else treat_call
    gain = treat_call.astype(float) - base_call.astype(float)  # resistant calls added by treatment
    if endpoint == "balanced_accuracy":
        return np.where(resistant, gain / (2 * n_r), -gain / (2 * n_s))
    if endpoint == "vme":
        return np.where(resistant, gain / n_r, 0.0)
    if endpoint == "me":
        return np.where(resistant, 0.0, -gain / n_s)
    return None


def endpoint_delta(endpoint, y, base, treat, null, base_call=None, treat_call=None):
    """Return (delta, None), or (None, reason) when a needed class is missing."""
    reason = undefined_reason(endpoint, y)
    if reason:
        return None, reason
    terms = contributions(endpoint, y, base, treat, null, base_call, treat_call)
    if terms is not None:
        return float(terms.sum()), None
    metric = auroc if endpoint == "auroc" else auprc
    return metric(y, treat) - metric(y, base), None


def sign_flip_p(sums, *, n_flips=20000, seed=0, exact_max=EXACT_MAX):
    """Two-sided p for H0: every cluster sum is symmetric about zero.

    Enumerates all 2^k sign patterns when k <= exact_max (p = hits / 2^k);
    otherwise draws n_flips seeded patterns (p = (1 + hits) / (1 + n_flips)).
    Returns (p, "exact" | "monte_carlo").
    """
    sums = np.asarray(sums, dtype=float)
    k = len(sums)
    observed = abs(sums.sum()) - 1e-9 * np.abs(sums).sum()
    rows = max(1, CELLS // k)
    hits = 0
    if k <= exact_max:
        bits = np.arange(k)
        for start in range(0, 2 ** k, rows):
            codes = np.arange(start, min(start + rows, 2 ** k))
            signs = 1.0 - 2.0 * ((codes[:, None] >> bits) & 1)
            hits += np.count_nonzero(np.abs(signs @ sums) >= observed)
        return hits / 2 ** k, "exact"
    rng = np.random.default_rng(seed)
    for start in range(0, n_flips, rows):
        signs = 1.0 - 2.0 * rng.integers(0, 2, size=(min(rows, n_flips - start), k))
        hits += np.count_nonzero(np.abs(signs @ sums) >= observed)
    return (1 + hits) / (1 + n_flips), "monte_carlo"


def arm_swap_p(endpoint, paired, *, n_flips=20000, seed=0, exact_max=EXACT_MAX):
    """Two-sided p from swapping the arms inside clusters, recomputing the endpoint.

    Exact when k <= exact_max and 2^k <= n_flips, else n_flips seeded swaps.
    Returns (p, "exact" | "monte_carlo").
    """
    k = len(paired.clusters)

    def delta(swap):
        moved = swap[paired.cluster_index]
        return endpoint_delta(endpoint, paired.y, np.where(moved, paired.treat, paired.base),
                              np.where(moved, paired.base, paired.treat), paired.null,
                              np.where(moved, paired.treat_call, paired.base_call),
                              np.where(moved, paired.base_call, paired.treat_call))[0]

    observed = delta(np.zeros(k, dtype=bool))
    if observed is None:
        raise ValueError(f"{endpoint} is undefined for {paired.dataset}")
    if k <= exact_max and 2 ** k <= n_flips:
        patterns = ((np.arange(2 ** k)[:, None] >> np.arange(k)) & 1).astype(bool)
        method = "exact"
    else:
        patterns = np.random.default_rng(seed).integers(0, 2, size=(n_flips, k)).astype(bool)
        method = "monte_carlo"
    threshold = abs(observed) * (1 - 1e-9) - 1e-15
    hits = sum(abs(delta(swap)) >= threshold for swap in patterns)
    return (hits / len(patterns), method) if method == "exact" else ((1 + hits) / (1 + n_flips), method)


def randomisation_p(endpoint, paired, *, n_flips=20000, seed=0, exact_max=EXACT_MAX):
    """Sign flip of cluster sums where the endpoint is a sum, else the arm swap."""
    terms = contributions(endpoint, paired.y, paired.base, paired.treat, paired.null, paired.base_call,
                          paired.treat_call)
    if terms is None:
        return arm_swap_p(endpoint, paired, n_flips=n_flips, seed=seed, exact_max=exact_max)
    sums = np.bincount(paired.cluster_index, weights=terms, minlength=len(paired.clusters))
    return sign_flip_p(sums, n_flips=n_flips, seed=seed, exact_max=exact_max)


def _t_central(t, df):
    """P(|T| <= t) for Student's t with integer df >= 1 (Abramowitz and Stegun 26.7.3, 26.7.4)."""
    theta = math.atan(t / math.sqrt(df))
    s, c2 = math.sin(theta), math.cos(theta) ** 2
    if df % 2:
        total = term = math.cos(theta) if df > 1 else 0.0
        for j in range(1, (df - 1) // 2):
            term *= c2 * (2 * j) / (2 * j + 1)
            total += term
        return 2 / math.pi * (theta + s * total)
    total = term = 1.0
    for j in range(1, df // 2):
        term *= c2 * (2 * j - 1) / (2 * j)
        total += term
    return s * total


@lru_cache(maxsize=None)
def t_quantile(q, df):
    """The q quantile (0.5 < q < 1) of Student's t with integer df >= 1, by bisection."""
    if not 0.5 < q < 1 or type(df) is not int or df < 1:
        raise ValueError(f"t_quantile needs 0.5 < q < 1 and an integer df >= 1, not {q!r}, {df!r}")
    target, lo, hi = 2 * q - 1, 0.0, 1.0
    while _t_central(hi, df) < target:
        lo, hi = hi, 2 * hi
    for _ in range(200):
        mid = (lo + hi) / 2
        if _t_central(mid, df) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-13 * hi:
            break
    return (lo + hi) / 2


def jackknife_t(estimate, replicates, level):
    """estimate +/- t((1 + level) / 2, k - 1) * SE_jackknife over the k leave-one-out replicates."""
    r = np.asarray(replicates, dtype=float)
    k = len(r)
    if k < 2 or not np.all(np.isfinite(r)):
        raise ValueError("a jackknife interval needs at least 2 defined leave-one-out replicates")
    se = math.sqrt((k - 1) / k * float(np.sum((r - r.mean()) ** 2)))
    half = t_quantile((1 + level) / 2, k - 1) * se
    return [float(estimate - half), float(estimate + half)]


def leave_one_out(paired, endpoints):
    """{endpoint: the k leave-one-cluster-out deltas, NaN where undefined}; both arms, calls and the null recomputed."""
    k = len(paired.clusters)
    out = {endpoint: np.full(k, np.nan) for endpoint in endpoints}
    for i in range(k):
        keep = paired.cluster_index != i
        arrays = (paired.y[keep], paired.base[keep], paired.treat[keep], paired.null[keep],
                  paired.base_call[keep], paired.treat_call[keep])
        for endpoint in endpoints:
            delta, _ = endpoint_delta(endpoint, *arrays)
            if delta is not None:
                out[endpoint][i] = delta
    return out


def _intervals(estimate, replicates):
    ci, tost = jackknife_t(estimate, replicates, CI_LEVEL), jackknife_t(estimate, replicates, TOST_LEVEL)
    if not ci[0] <= tost[0] <= tost[1] <= ci[1]:
        raise ArithmeticError(f"intervals do not nest: ci {ci}, ci_tost {tost}")
    return ci, tost


def _margin(value, what):
    if not _real(value) or not value > 0:
        raise ValueError(f"{what} needs a predeclared finite margin > 0, not {value!r}")
    return float(value)


def _check_development(family, datasets, development):
    if isinstance(development, str):
        raise TypeError("development must be a collection of dataset ids, not one string")
    refused = sorted(set(datasets) & set(development))
    if refused and family in CONFIRMATORY:
        raise DevelopmentDatasetError(f"development datasets cannot enter family {family}: {refused}")


def contrast_key(family, learner, split, baseline, treatment, endpoint):
    """The key grammar of contract section 12: every part non-empty and free of whitespace,
    ':', '/' and backticks; the arms also free of '->'."""
    parts = {"family": family, "learner": learner, "split": split, "baseline": baseline,
             "treatment": treatment, "endpoint": endpoint}
    for name, part in parts.items():
        if (not isinstance(part, str) or not part or re.search(r"[\s:/`]", part)
                or name in ("baseline", "treatment") and "->" in part):
            raise ValueError(f"contrast key part {name} {part!r} must be non-empty without whitespace, ':', '/' "
                             f"or a backtick, and an arm without '->'")
    return f"{family}:{learner}:{split}:{baseline}->{treatment}:{endpoint}"


def compare(paired, *, family, endpoints, margins, development, n_flips=20000, seed=0, exact_max=EXACT_MAX):
    """One pending contrast record per endpoint; pass them with their family's planned ids to finalize().

    Family C takes log_loss only (the per-dataset contrasts of an A or B
    comparison); family D takes any endpoint and carries no p.
    """
    if family not in CONTRAST_FAMILIES:
        raise ValueError(f"per-dataset contrasts are family C or D, not {family!r}; A and B hold pooled records")
    if not DATASET_ID.fullmatch(str(paired.dataset)):
        raise ValueError(f"dataset id {paired.dataset!r} does not match {DATASET_ID.pattern}")
    _check_development(family, [paired.dataset], development)
    endpoints = list(endpoints)
    if len(set(endpoints)) != len(endpoints):
        raise ValueError("endpoints repeat")
    for endpoint in endpoints:
        if endpoint not in NEEDS:
            raise ValueError(f"unknown endpoint {endpoint!r}; expected one of {ENDPOINTS}")
        if family == "C" and endpoint != PRIMARY:
            raise ValueError(f"family C uses {PRIMARY} only; {endpoint} is family D")
        _margin(margins.get(endpoint), endpoint)
    k = len(paired.clusters)
    loo = leave_one_out(paired, endpoints) if k >= 2 else {endpoint: np.full(k, np.nan) for endpoint in endpoints}
    auroc_delta = auroc_ci = None
    if PRIMARY in endpoints:
        delta, why = endpoint_delta("auroc", paired.y, paired.base, paired.treat, paired.null)
        replicates = leave_one_out(paired, ["auroc"])["auroc"] if k >= 2 else np.array([])
        if why is None and k >= 2 and np.all(np.isfinite(replicates)):
            auroc_delta, auroc_ci = float(delta), jackknife_t(delta, replicates, CI_LEVEL)
    records = []
    for i, endpoint in enumerate(endpoints):
        estimate, reason = endpoint_delta(endpoint, paired.y, paired.base, paired.treat, paired.null,
                                          paired.base_call, paired.treat_call)
        undefined = int(np.isnan(loo[endpoint]).sum())
        if reason is None and k < 2:
            reason = "one cluster: a jackknife interval needs at least 2"
        elif reason is None and undefined > MAX_UNDEFINED_SHARE * k:
            reason = UNDEFINED_SHARE_REASON
        elif reason is None and undefined:
            reason = f"{undefined} of {k} leave-one-cluster-out estimates undefined; the jackknife needs all"
        ci = tost = p = None
        if reason is None:
            ci, tost = _intervals(estimate, loo[endpoint])
            if family != "D":
                p = float(randomisation_p(endpoint, paired, n_flips=n_flips, seed=[seed, 1, i], exact_max=exact_max)[0])
        else:
            estimate = None
        key = contrast_key(family, paired.learner, paired.split, paired.baseline, paired.treatment, endpoint)
        records.append({
            "id": f"{paired.dataset}/{key}", "family": family, "baseline": paired.baseline,
            "treatment": paired.treatment, "learner": paired.learner, "split": paired.split, "endpoint": endpoint,
            "estimate": estimate, "ci": ci, "ci_tost": tost, "ci_method": JACKKNIFE if reason is None else "none",
            "unit": "cluster", "n_units": k, "few_clusters": k < FEW_CLUSTERS, "p": p, "p_adj": None,
            "adjust": FAMILIES[family], "margin": float(margins[endpoint]), "decision": None,
            "claim_level": CLAIM_LEVEL, "undefined_resamples": undefined, "undefined_reason": reason,
            "interval_label": INTERVAL_LABEL, "dataset": paired.dataset, "pairing": paired.pairing,
            "auroc_delta": auroc_delta if endpoint == PRIMARY else None,
            "auroc_ci": auroc_ci if endpoint == PRIMARY else None, "_pending": True,
        })
    return records


@dataclass(frozen=True, eq=False)
class PairedUnits:
    """Two arms' paired rows per repeat (pair_units): one entry per (repeat, sample) both score."""
    dataset: str
    learner: str
    split: str
    baseline: str
    treatment: str
    pairing: str
    units: tuple             # (repeat, outer fold), one per unit
    unit_index: np.ndarray   # each entry's position in units
    unit_ratio: np.ndarray   # per unit, n_test / n_train: its paired genomes over the baseline's training genomes
    y: np.ndarray
    base: np.ndarray
    treat: np.ndarray


def _round_half_up(x):
    return math.floor(x + 0.5)


def balanced_counts(have, want):
    """`grm arm`'s class-balanced subsample sizes (plan.rs balanced_counts): of `have` = [resistant,
    susceptible] common genomes in a fold whose whole baseline test fold has `want`, keep the scarcer class
    whole and match the other to want's ratio, rounding half away from zero."""
    (r, s), (wr, ws) = have, want
    if wr == 0 or ws == 0:
        return [0 if wr == 0 else r, 0 if ws == 0 else s]
    if r * ws <= s * wr:
        return [r, min(_round_half_up(r * ws / wr), s)]
    return [min(_round_half_up(s * wr / ws), r), s]


def _check_balanced(repeat, base_rows, treat_ids):
    """REVIEW-ICODE-PY L1: under common_outer_test_set_balanced the treatment is scored on the balanced draw of
    the baseline's whole test fold, in every fold: each baseline fold of the repeat holds scored treatment
    samples of both classes, in counts balanced_counts gives for some genomes the arm holds."""
    want, got = defaultdict(lambda: [0, 0]), defaultdict(lambda: [0, 0])
    for sample, row in base_rows.items():
        want[row["fold"]][1 - row["y"]] += 1
        if sample in treat_ids:
            got[row["fold"]][1 - row["y"]] += 1
    for fold in sorted(want, key=str):
        (kr, ks), (wr, ws) = got[fold], want[fold]
        if kr == 0 or ks == 0:
            raise PairingError("balance", f"repeat {repeat!r} fold {fold!r}: the treatment scores {kr} resistant and "
                                          f"{ks} susceptible genomes; a balanced subsample of the baseline's test fold "
                                          f"({wr} and {ws}) scores both classes in every fold")
        # balanced_counts keeps the scarcer class whole and rounds the other, never past what the arm holds, so
        # its output is exactly the pairs where one count is the rounded match of the other.
        if ks != _round_half_up(kr * ws / wr) and kr != _round_half_up(ks * wr / ws):
            raise PairingError("balance", f"repeat {repeat!r} fold {fold!r}: {kr} resistant and {ks} susceptible is "
                                          f"not a balanced subsample of the baseline's test fold ({wr} and {ws})")


def pair_units(rows, *, dataset, learner, split, baseline, treatment, pairing="same_plan"):
    """Pair two arms per repeat under `pairing` (one of PAIRINGS) for compare_units(); raise PairingError
    otherwise. Both arms need the same repeats with the same fold plan in each. Under same_plan and
    k31_groups the arms score the same samples, the same ones in every repeat; under SUBSET_PAIRINGS the
    treatment's samples of each repeat are a subset of the baseline's, which is restricted to them in that
    repeat (a cohort arm's subsample may differ by repeat). A sample keeps one label and one cluster
    across arms and repeats, and both arms put it in one fold per repeat."""
    tables = _arm_tables(rows, dataset, learner, split, baseline, treatment, pairing)
    plans, samples = {}, {}
    for arm, table in tables.items():
        if not table:
            raise PairingError("samples", f"arm {arm!r} has no rows")
        hashes, ids = defaultdict(set), defaultdict(set)
        for (repeat, sample), row in table.items():
            hashes[repeat].add(row["fold_plan_sha256"])
            ids[repeat].add(sample)
        for repeat, h in hashes.items():
            if len(h) != 1:
                raise PairingError("fold_plan", f"arm {arm!r} repeat {repeat!r} has {len(h)} fold plans")
        plans[arm], samples[arm] = {r: h.pop() for r, h in hashes.items()}, ids
    if plans[baseline] != plans[treatment]:
        differ = sorted((r for r in plans[baseline].keys() | plans[treatment].keys()
                         if plans[baseline].get(r) != plans[treatment].get(r)), key=str)
        raise PairingError("fold_plan", f"fold_plan_sha256 differs or is missing for repeats {differ}")
    subset = pairing in SUBSET_PAIRINGS
    for arm in (baseline,) if subset else (baseline, treatment):
        if len({frozenset(ids) for ids in samples[arm].values()}) != 1:
            raise PairingError("samples", f"arm {arm!r} covers different samples in different repeats")
    labels, clusters, n_train = {}, {}, {}
    units, unit_of, y, base, treat = {}, [], [], [], []
    for repeat in sorted(plans[baseline], key=str):
        base_ids, treat_ids = samples[baseline][repeat], samples[treatment][repeat]
        if subset and not treat_ids <= base_ids:
            raise PairingError("samples", f"repeat {repeat!r}: {len(treat_ids - base_ids)} samples only in "
                                          f"{treatment!r}; {pairing} scores the treatment on the baseline's test genomes")
        if not subset and base_ids != treat_ids:
            raise PairingError("samples", f"repeat {repeat!r}: {len(base_ids - treat_ids)} samples only in "
                                          f"{baseline!r}, {len(treat_ids - base_ids)} only in {treatment!r}")
        base_rows = {s: tables[baseline][(repeat, s)] for s in base_ids}
        # The baseline's training genomes in each fold: its whole repeat set outside the fold.
        in_fold = defaultdict(int)
        for row in base_rows.values():
            in_fold[row["fold"]] += 1
        for fold, n in in_fold.items():
            n_train[(repeat, fold)] = len(base_rows) - n
        if pairing == "common_outer_test_set_balanced":
            _check_balanced(repeat, base_rows, treat_ids)
        prevalence = {}
        if not subset:
            folds = {row["fold"] for row in base_rows.values()}
            for fold in folds:
                train = [row["y"] for row in base_rows.values() if row["fold"] != fold]
                if not train:
                    raise PairingError("fold_plan", f"repeat {repeat!r} has one fold, so no training samples")
                prevalence[fold] = sum(train) / len(train)
        for sample in sorted(treat_ids, key=str):
            b, t = base_rows[sample], tables[treatment][(repeat, sample)]
            if labels.setdefault(sample, b["y"]) != b["y"] or t["y"] != b["y"]:
                raise PairingError("y", f"sample {sample!r} has labels {sorted({labels[sample], b['y'], t['y']})} "
                                        "across arms or repeats")
            if b["cluster"] is None or t["cluster"] != b["cluster"] \
                    or clusters.setdefault(sample, b["cluster"]) != b["cluster"]:
                raise PairingError("cluster", f"sample {sample!r} has clusters "
                                              f"{sorted(map(str, {clusters.get(sample), b['cluster'], t['cluster']}))}")
            if b["fold"] != t["fold"]:
                raise PairingError("fold", f"sample {sample!r} repeat {repeat!r} is in different folds under one "
                                           "fold plan")
            for row in (b, t):
                _check_called(row, None if subset else prevalence[b["fold"]])
            unit_of.append(units.setdefault((repeat, b["fold"]), len(units)))
            y.append(b["y"])
            base.append(float(b["score"]))
            treat.append(float(t["score"]))
    order = sorted(units, key=lambda u: (str(u[0]), str(u[1])))
    position = np.array([order.index(u) for u in units])
    index = position[np.array(unit_of, dtype=np.int64)]
    n_test = np.bincount(index, minlength=len(order))
    if any(n_train[u] == 0 for u in order):
        raise PairingError("fold_plan", "a unit's fold holds every baseline genome of its repeat, so no training samples")
    return PairedUnits(dataset=dataset, learner=learner, split=split, baseline=baseline, treatment=treatment,
                       pairing=pairing, units=tuple(order), unit_index=index,
                       unit_ratio=n_test / np.array([n_train[u] for u in order], dtype=float),
                       y=np.array(y, dtype=np.int64), base=np.array(base), treat=np.array(treat))


def unit_t(values, ratio=0.0):
    """(mean, 95% interval, 90% interval, two-sided p) of the corrected resampled t over unit values (Nadeau and
    Bengio): SE = sqrt((1/J + ratio) s^2) with s^2 the sample variance, df = J - 1; ratio 0 is the plain paired
    t. With every value equal (to rounding: a spread under 1e-12 of the largest) there is no interval: all zero
    gives ([0, 0], [0, 0], p = 1), any other value None."""
    v = np.asarray(values, dtype=float)
    j = len(v)
    if j < 2:
        raise ValueError("a paired t needs at least 2 units")
    if not (math.isfinite(ratio) and ratio >= 0):
        raise ValueError(f"the correction's ratio n_test / n_train must be finite and >= 0, not {ratio!r}")
    mean = float(v.mean())
    se = math.sqrt((1 / j + ratio) * float(v.var(ddof=1)))
    if se <= 1e-12 * float(np.abs(v).max()):
        return (0.0, [0.0, 0.0], [0.0, 0.0], 1.0) if not v.any() else None
    t = abs(mean) / se
    half = {level: t_quantile((1 + level) / 2, j - 1) * se for level in (CI_LEVEL, TOST_LEVEL)}
    p = min(1.0, max(0.0, 1.0 - _t_central(t, j - 1)))
    return (mean, [mean - half[CI_LEVEL], mean + half[CI_LEVEL]], [mean - half[TOST_LEVEL], mean + half[TOST_LEVEL]],
            p)


def compare_units(paired, *, family, margins, development, test):
    """One pending log_loss record over the (repeat, fold) units (family C decides, D is descriptive), by `test`:
    UNIT_T, the paired t over the unit means, or UNIT_CORRECTED_T, Nadeau and Bengio's corrected resampled t
    (LEARNER_TESTS gives the input comparison's choice per learner); the AUROC delta's unit interval uses the same
    test. Pass the record with its family's planned ids to finalize()."""
    if test not in UNIT_INTERVAL_LABELS:
        raise ValueError(f"test {test!r} is not one of {sorted(UNIT_INTERVAL_LABELS)}")
    if family not in CONTRAST_FAMILIES:
        raise ValueError(f"per-dataset contrasts are family C or D, not {family!r}; A and B hold pooled records")
    if not DATASET_ID.fullmatch(str(paired.dataset)):
        raise ValueError(f"dataset id {paired.dataset!r} does not match {DATASET_ID.pattern}")
    _check_development(family, [paired.dataset], development)
    margin = _margin(margins.get(PRIMARY), PRIMARY)
    j = len(paired.units)
    counts = np.bincount(paired.unit_index, minlength=j)
    delta = sample_log_loss(paired.y, paired.base) - sample_log_loss(paired.y, paired.treat)
    means = np.bincount(paired.unit_index, weights=delta, minlength=j) / counts
    rho = float(np.mean(paired.unit_ratio)) if test == UNIT_CORRECTED_T else 0.0
    estimate = ci = tost = p = auroc_delta = auroc_ci = None
    reason = "one unit: a t interval needs at least 2" if j < 2 else None
    if reason is None:
        result = unit_t(means, rho)
        if result is None:
            reason = "every unit's mean delta is the same nonzero value: no t interval"
        else:
            estimate, ci, tost, p = result
    if reason is None:
        aurocs = []
        for u in range(j):
            here = paired.unit_index == u
            y = paired.y[here]
            if undefined_reason("auroc", y) is not None:
                break
            aurocs.append(auroc(y, paired.treat[here]) - auroc(y, paired.base[here]))
        else:
            result = unit_t(aurocs, rho)
            if result is not None:
                auroc_delta, auroc_ci = result[0], result[1]
            elif len(set(aurocs)) == 1:
                auroc_delta, auroc_ci = aurocs[0], [aurocs[0], aurocs[0]]
    key = contrast_key(family, paired.learner, paired.split, paired.baseline, paired.treatment, PRIMARY)
    return [{
        "id": f"{paired.dataset}/{key}", "family": family, "baseline": paired.baseline,
        "treatment": paired.treatment, "learner": paired.learner, "split": paired.split, "endpoint": PRIMARY,
        "estimate": estimate, "ci": ci, "ci_tost": tost, "ci_method": test if reason is None else "none",
        "unit": UNIT, "n_units": j, "few_clusters": False, "p": p if family != "D" else None, "p_adj": None,
        "adjust": FAMILIES[family], "margin": margin, "decision": None, "claim_level": CLAIM_LEVEL,
        "undefined_resamples": 0, "undefined_reason": reason, "interval_label": UNIT_INTERVAL_LABELS[test],
        "dataset": paired.dataset, "pairing": paired.pairing, "auroc_delta": auroc_delta, "auroc_ci": auroc_ci,
        "_pending": True,
    }]


def _shared_contrast(records):
    if not records:
        raise ValueError("no records")
    if any("dataset" not in r or "contrast" in r for r in records):
        raise ValueError("pass per-dataset contrast records from compare()")
    shared = {(r["baseline"], r["treatment"], r["learner"], r["split"], r["endpoint"]) for r in records}
    if len(shared) != 1:
        raise ValueError("records must share baseline, treatment, learner, split and endpoint")
    datasets = [r["dataset"] for r in records]
    if len(set(datasets)) != len(datasets):
        raise ValueError("a dataset appears twice")
    return shared.pop(), datasets


def dataset_summary(records, margin):
    """Median per-dataset delta and where each dataset's 95% CI sits against +/-margin."""
    _shared_contrast(records)
    counts = dict.fromkeys(("above", "below", "inside", "straddling", "undefined"), 0)
    for r in records:
        ci = r["ci"]
        bucket = ("undefined" if ci is None else "above" if ci[0] > margin else "below" if ci[1] < -margin
                  else "inside" if -margin < ci[0] and ci[1] < margin else "straddling")
        counts[bucket] += 1
    defined = [r["estimate"] for r in records if r["estimate"] is not None]
    return {"median_dataset_delta": float(np.median(defined)) if defined else None, "dataset_ci_counts": counts}


def pool(records, *, species, family, margin, development, missing=None, n_flips=20000, seed=0,
         exact_max=EXACT_MAX):
    """Pool per-dataset r2_ll records (family D, from compare()): mean within species, then across species.

    The records must agree on family, margin, endpoint, learner, split and
    arms, and their margin must be the pooled margin. `missing` maps planned
    datasets that have no record (for example a refused pairing) to a reason,
    so every planned dataset either contributes or is excluded.

    p (families A and B) is the exact species sign flip. It is valid for the
    null "no effect in any species", with each species' estimate error
    symmetric about zero; it is not a test that the average species effect is
    zero, and it rejects too often when species effects average zero but are
    skewed (contract section 15). A pooled improves (or worse) means that null
    is rejected at Holm-adjusted p < 0.05 and the mean species effect is
    positive (or negative).
    """
    if family not in POOLED_FAMILIES:
        raise ValueError(f"pooled records are family A, B or D, not {family!r}")
    margin = _margin(margin, "pooling")
    (baseline, treatment, learner, split, endpoint), datasets = _shared_contrast(records)
    if endpoint != POOLED_ENDPOINT:
        raise ValueError(f"the pooled scale is {POOLED_ENDPOINT}, not {endpoint}")
    if len({r["family"] for r in records}) != 1:
        raise ValueError("input records disagree on family")
    if {r["margin"] for r in records} != {margin}:
        raise ValueError(f"input records' margins {sorted({r['margin'] for r in records})} must all be the pooled margin {margin}")
    missing = dict(missing or {})
    if set(missing) & set(datasets):
        raise ValueError("a dataset is both recorded and missing")
    _check_development(family, datasets + list(missing), development)
    unmapped = sorted(set(datasets) - set(species))
    if unmapped:
        raise ValueError(f"no species for datasets {unmapped}")

    by_species = defaultdict(list)
    excluded = [{"dataset": dataset, "reason": str(reason)} for dataset, reason in missing.items()]
    for r in records:
        if r["estimate"] is None:
            excluded.append({"dataset": r["dataset"], "reason": r["undefined_reason"]})
        else:
            by_species[species[r["dataset"]]].append((r["dataset"], r["estimate"]))
    if len(by_species) == 1:
        for dataset, _ in by_species.popitem()[1]:
            excluded.append({"dataset": dataset, "reason": "one contributing species: no species-level interval"})
    excluded.sort(key=lambda row: row["dataset"])
    means = np.array([np.mean([value for _, value in by_species[name]]) for name in sorted(by_species)])
    s = len(means)
    estimate = ci = tost = p = None
    if s:
        estimate = float(means.mean())
        ci, tost = _intervals(estimate, (means.sum() - means) / (s - 1))
        if family != "D":
            p = float(sign_flip_p(means, n_flips=n_flips, seed=[seed, 1], exact_max=exact_max)[0])
    return {
        "contrast": contrast_key(family, learner, split, baseline, treatment, endpoint), "family": family,
        "unit": "species", "n_units": s, "n_datasets": len(datasets) + len(missing), "n_excluded": len(excluded),
        "excluded": excluded, "estimate": estimate, "ci": ci, "ci_tost": tost, "p": p, "p_adj": None,
        "adjust": FAMILIES[family], "margin": margin, "decision": None, "interval_label": INTERVAL_LABEL,
        "_pending": True,
    }


def holm(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="mergesort")
    adjusted = np.empty(len(p))
    adjusted[order] = np.maximum.accumulate((len(p) - np.arange(len(p))) * p[order]).clip(max=1.0)
    return adjusted.tolist()


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p, kind="mergesort")
    ranked = p[order] * (len(p) / np.arange(1, len(p) + 1))  # factor >= 1 keeps p_adj >= p
    adjusted = np.empty(len(p))
    adjusted[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(max=1.0)
    return adjusted.tolist()


def decide(ci, p_adj):
    """improves/worse need the 95% CI off zero and p_adj < 0.05; anything else, including
    family D (p_adj None), is inconclusive. There is no equivalent (contract section 15)."""
    if ci is not None and p_adj is not None and p_adj < ALPHA:
        if ci[0] > 0:
            return "improves"
        if ci[1] < 0:
            return "worse"
    return "inconclusive"


def _name(r):
    return r["contrast"] if "contrast" in r else r["id"]


def _planned_family(planned_id):
    family = planned_id.split("/", 1)[-1].split(":", 1)[0]
    if family not in FAMILIES:
        raise ValueError(f"planned id {planned_id!r} does not name a family")
    return family


def _check_final(r):
    defined = r["estimate"] is not None
    for key in ("estimate", "ci", "ci_tost") + (("p", "p_adj") if r["family"] != "D" else ()):
        if (r[key] is not None) != defined:
            raise ValueError(f"{_name(r)}: {key} must be {'set' if defined else 'null'}")
    if "ci_method" in r and (r["ci_method"] == "none") == defined:
        raise ValueError(f"{_name(r)}: ci_method must be {'an interval method' if defined else repr('none')}")
    if defined:
        values = [r["estimate"], *r["ci"], *r["ci_tost"]] + ([r["p"], r["p_adj"]] if r["p"] is not None else [])
        if not all(_real(v) for v in values):
            raise ValueError(f"{_name(r)}: values must be finite")
        if not r["ci"][0] <= r["ci_tost"][0] <= r["ci_tost"][1] <= r["ci"][1]:
            raise ValueError(f"{_name(r)}: intervals must nest")
        if r["p"] is not None and not r["p"] <= r["p_adj"] <= 1:
            raise ValueError(f"{_name(r)}: p_adj must lie in [p, 1]")


def finalize(records, planned_ids, missing=None, unadjusted=()):
    """Adjust p within each family, decide, and return the records with their public keys.

    `planned_ids` lists every planned contrast id ("<dataset>/<key>") and
    pooled key of the families being finalized; `missing` maps planned ids
    with no record to a reason. A call lacking a planned id, or holding an
    unplanned record, is refused, so a family is adjusted in one complete
    call. Families are adjusted separately (Holm in A and B, BH in C, none in
    D); undefined records and missing ids enter their family with p = 1.
    `unadjusted` names planned family C ids judged on their own p (the input
    comparison's sentinels, main's ruling 3 after VERIFY-H1STATS): they stay
    out of BH, and their p_adj is their p with adjust "none".
    """
    if isinstance(planned_ids, str) or isinstance(unadjusted, str):
        raise TypeError("planned_ids and unadjusted must be collections of ids, not one string")
    planned, missing, unadjusted = list(planned_ids), dict(missing or {}), set(unadjusted)
    loose = sorted(i for i in unadjusted if i not in planned or _planned_family(i) != "C")
    if loose:
        raise ValueError(f"unadjusted ids must be planned family C ids: {loose}")
    if len(set(planned)) != len(planned):
        raise ValueError("planned ids repeat")
    names = [_name(r) for r in records]
    if len(set(names)) != len(names):
        raise ValueError("contrast ids repeat")
    for r, name in zip(records, names):
        if not r.get("_pending"):
            raise ValueError(f"{name}: not a pending record from compare() or pool()")
        family, pooled = r["family"], "contrast" in r
        if family not in (POOLED_FAMILIES if pooled else CONTRAST_FAMILIES):
            raise ValueError(f"{name}: family {family!r} does not hold {'pooled' if pooled else 'per-dataset'} records")
        if family == "C" and r["endpoint"] != PRIMARY:
            raise ValueError(f"{name}: family C uses {PRIMARY} only")
        if family == "D" and r["p"] is not None:
            raise ValueError(f"{name}: family D carries no p-values")
        if family != "D" and r["estimate"] is not None and r["p"] is None:
            raise ValueError(f"{name}: family {family} needs a p-value")
    unplanned = sorted(set(names) - set(planned))
    if unplanned:
        raise ValueError(f"records not in planned_ids: {unplanned}")
    stray = sorted(set(missing) - (set(planned) - set(names)))
    if stray:
        raise ValueError(f"missing ids that are unplanned or have a record: {stray}")
    absent = sorted(set(planned) - set(names) - set(missing))
    if absent:
        raise ValueError(f"incomplete call: planned ids with no record and no reason in missing: {absent}")
    if not all(isinstance(reason, str) and reason.strip() for reason in missing.values()):
        raise ValueError("every missing id needs a non-empty reason")

    out = [{key: r[key] for key in (POOLED_KEYS if "contrast" in r else RECORD_KEYS)} for r in records]
    for family, adjust in FAMILIES.items():
        if adjust == "none":
            continue
        tested = [i for i, (r, name) in enumerate(zip(records, names))
                  if r["family"] == family and r["p"] is not None and name not in unadjusted]
        untested = (sum(r["family"] == family and r["p"] is None and name not in unadjusted
                        for r, name in zip(records, names))
                    + sum(_planned_family(pid) == family and pid not in unadjusted for pid in missing))
        if tested:
            adjusted = {"holm": holm, "bh": bh}[adjust]([records[i]["p"] for i in tested] + [1.0] * untested)
            for i, value in zip(tested, adjusted):
                out[i]["p_adj"] = float(value)
    for r, name in zip(out, names):
        if name in unadjusted:
            r["adjust"], r["p_adj"] = "none", r["p"]
    for r in out:
        r["decision"] = decide(r["ci"], r["p_adj"])
        _check_final(r)
    return out


def finalize_sets(sets, unadjusted=()):
    """finalize() over several (records, planned_ids, missing) sets as one adjustment, so BH's family C spans
    them all (main's ruling 4 after VERIFY-H1STATS: one family over every lineage contrast of a run's
    species). Returns each set's finalized records, in the sets' order."""
    records, planned, missing, owner = [], [], {}, []
    for k, (recs, ids, miss) in enumerate(sets):
        miss = dict(miss or {})
        twice = sorted(set(miss) & set(missing))
        if twice:
            raise ValueError(f"missing ids in two sets: {twice}")
        records += list(recs)
        planned += list(ids)
        missing.update(miss)
        owner += [k] * len(recs)
    out = finalize(records, planned, missing, unadjusted)
    return [[r for r, o in zip(out, owner) if o == k] for k in range(len(sets))]
