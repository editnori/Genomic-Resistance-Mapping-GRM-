"""Synthetic leakage canaries for hyperparameter search and evaluation.

Each check takes plain arrays, mappings and callbacks and returns a
CanaryResult whose offenders name the samples, groups, folds, arms or tasks
at fault. Offenders and details hold plain JSON values, so a report can embed
them as they are.

Conventions:
- Samples are positions 0..n-1. Pass sample_ids to name them in offenders.
- A per-sample fold array holds the outer test fold of each sample.
- Groups are a per-sample sequence of group ids, where None means the sample
  has no group, or a mapping of group id to member positions, which may hold
  empty groups.
- Inner partitions (inner folds, tuning holdouts) are per-sample label
  sequences with None for every sample outside the outer training set.
"""
import hashlib
import json
import math
import numbers
from collections import Counter, defaultdict
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class CanaryResult:
    canary: str
    passed: bool
    offenders: list = field(default_factory=list)
    details: dict = field(default_factory=dict)


def _result(canary, offenders, details=None, problems=()):
    """Pass only with no offenders; problems are canary-level failures."""
    details = dict(details or {})
    if problems:
        details["problems"] = list(problems)
    return CanaryResult(canary, not offenders and not problems, offenders, details)


def _plain(value):
    """Convert numpy scalars, arrays and tuples to plain JSON values."""
    if isinstance(value, np.ndarray):
        return [_plain(v) for v in value.tolist()]
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, dict):
        return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.floating):
        return float(value)
    return value


def _order(value):
    return (type(value).__name__, value)


def _names(positions, sample_ids):
    return [_plain(sample_ids[i]) if sample_ids is not None else int(i) for i in positions]


def _members(groups, n):
    """Map group id to member positions; ungrouped samples are left out."""
    if isinstance(groups, dict):
        members = {}
        for gid, idx in groups.items():
            idx = [int(i) for i in idx]
            if any(i < 0 or i >= n for i in idx):
                raise ValueError(f"group {gid!r} has a member outside 0..{n - 1}")
            members[_plain(gid)] = idx
        return members
    groups = _plain(list(groups))
    if len(groups) != n:
        raise ValueError(f"groups has {len(groups)} entries for {n} samples")
    members = defaultdict(list)
    for i, gid in enumerate(groups):
        if gid is not None:
            members[gid].append(i)
    return dict(members)


def _spread(idx, labels, sample_ids, key):
    """Members of one group split by label; more than one entry is a crossing."""
    by_label = defaultdict(list)
    for i in idx:
        by_label[labels[i]].append(i)
    return [{key: label, "samples": _names(by_label[label], sample_ids)}
            for label in sorted(by_label, key=_order)]


# 1. Groups across outer folds

def check_groups_within_outer_folds(outer_fold, groups, sample_ids=None):
    """Fail when any group has members in two outer test folds.

    groups maps a group kind ("lineage", "duplicate") to that kind's groups.
    """
    folds = _plain(list(outer_fold))
    if any(f is None for f in folds):
        raise ValueError("every sample needs an outer test fold")
    offenders, counts = [], {}
    for kind, kind_groups in groups.items():
        members = _members(kind_groups, len(folds))
        counts[kind] = {"groups": len(members),
                        "empty": sum(1 for idx in members.values() if not idx)}
        for gid in sorted(members, key=_order):
            spread = _spread(members[gid], folds, sample_ids, "outer_fold")
            if len(spread) > 1:
                offenders.append({"kind": kind, "group": gid, "folds": spread})
    return _result("groups_within_outer_folds", offenders, {"groups": counts})


# 2. Groups across inner folds and tuning holdouts

def check_groups_within_inner_folds(outer_fold, inner_plans, groups, sample_ids=None):
    """Fail when an inner partition splits a group or touches the outer test set.

    inner_plans maps an outer fold id to {partition name: per-sample labels},
    for example {"inner_folds": [...], "holdout": ["fit", None, "holdout", ...]}.
    Every outer fold needs a plan with at least one partition. Every outer
    training sample needs a label; every outer test sample, None.
    """
    folds = _plain(list(outer_fold))
    n = len(folds)
    kinds = {kind: _members(kind_groups, n) for kind, kind_groups in groups.items()}
    offenders = [{"outer_fold": f, "problem": "missing_outer_fold"}
                 for f in sorted(set(folds) - set(inner_plans), key=_order)]
    for f in sorted(inner_plans, key=_order):
        if f not in folds:
            offenders.append({"outer_fold": f, "problem": "unknown_outer_fold"})
            continue
        if not inner_plans[f]:
            offenders.append({"outer_fold": f, "problem": "no_partitions"})
        train = {i for i in range(n) if folds[i] != f}
        for name in sorted(inner_plans[f]):
            labels = _plain(list(inner_plans[f][name]))
            if len(labels) != n:
                raise ValueError(f"outer fold {f!r} partition {name!r} has "
                                 f"{len(labels)} entries for {n} samples")
            where = {"outer_fold": f, "partition": name}
            in_test = [i for i in range(n) if i not in train and labels[i] is not None]
            if in_test:
                offenders.append({**where, "problem": "test_sample_in_plan",
                                  "samples": _names(in_test, sample_ids)})
            unassigned = [i for i in sorted(train) if labels[i] is None]
            if unassigned:
                offenders.append({**where, "problem": "training_sample_unassigned",
                                  "samples": _names(unassigned, sample_ids)})
            for kind, members in kinds.items():
                for gid in sorted(members, key=_order):
                    idx = [i for i in members[gid] if i in train and labels[i] is not None]
                    spread = _spread(idx, labels, sample_ids, "label")
                    if len(spread) > 1:
                        offenders.append({**where, "problem": "group_crosses",
                                          "kind": kind, "group": gid, "labels": spread})
    return _result("groups_within_inner_folds", offenders,
                   {"outer_folds_checked": len(inner_plans)})


# 3, 4 and 6. Invariance to the outer test samples' labels and features

def _perturb_labels(labels, test, rng, mode):
    out = labels.copy()
    values = np.unique(labels)
    if mode == "invert":
        # Rotate each test label to the next observed value; binary labels swap.
        out[test] = values[(np.searchsorted(values, labels[test]) + 1) % len(values)]
    elif mode == "permute":
        out[test] = rng.permutation(labels[test])
    elif mode == "resample":
        out[test] = rng.choice(values, size=len(test))
    else:
        raise ValueError(f"unknown label perturbation {mode!r}")
    return out


def _perturb_features(features, test, rng, mode):
    out = features.copy()
    if mode == "resample":
        # Draw each test entry from its column's observed values, so binary
        # k-mer columns stay binary while column totals change.
        rows = rng.integers(0, features.shape[0], size=(len(test), features.shape[1]))
        out[test] = features[rows, np.arange(features.shape[1])]
    elif mode == "permute_rows":
        out[test] = features[rng.permutation(test)]
    else:
        raise ValueError(f"unknown feature perturbation {mode!r}")
    return out


def _diff_keys(before, after):
    if isinstance(before, dict) and isinstance(after, dict):
        return {"changed": sorted(k for k in set(before) | set(after)
                                  if canonical_json(before.get(k)) != canonical_json(after.get(k)))}
    return {"changed": ["<output>"]}


def _diff_sets(before, after):
    before, after = set(_plain(list(before))), set(_plain(list(after)))
    return {"added": sorted(after - before, key=_order),
            "removed": sorted(before - after, key=_order)}


def _invariance(canary, outer_fold, labels, features, fn, target, modes, repeats, seed, diff):
    """Call fn(outer_fold_id, labels, features) on perturbed test samples.

    fn gets fresh copies on every call, so masking its inputs in place neither
    breaks the comparison nor touches the caller's arrays.
    """
    folds = np.asarray(_plain(list(outer_fold)), dtype=object)
    labels = np.asarray(labels)
    if features is None and target == "features":
        raise ValueError(f"{canary} perturbs test features, so it needs a features array")
    features = None if features is None else np.asarray(features)

    def call(f, y, x):
        return fn(f, y.copy(), None if x is None else x.copy())

    rng = np.random.default_rng(seed)
    offenders, problems, runs = [], [], 0
    for f in sorted(set(folds.tolist()), key=_order):
        test = np.flatnonzero(folds == f)
        baseline = call(f, labels, features)
        if canonical_json(call(f, labels, features)) != canonical_json(baseline):
            problems.append({"outer_fold": f, "problem": "callback_not_deterministic"})
            continue
        perturbed_any = False
        for mode in modes:
            # invert draws nothing at random, so a repeat would rerun the same input.
            for r in range(1 if mode == "invert" else repeats):
                if target == "labels":
                    y, x = _perturb_labels(labels, test, rng, mode), features
                    perturbed = not np.array_equal(y, labels)
                else:
                    y, x = labels, _perturb_features(features, test, rng, mode)
                    perturbed = not np.array_equal(x, features)
                if not perturbed:
                    continue
                perturbed_any, runs = True, runs + 1
                after = call(f, y, x)
                if canonical_json(after) != canonical_json(baseline):
                    offenders.append({"outer_fold": f, "perturbed": target,
                                      "perturbation": mode, "repeat": r,
                                      **diff(baseline, after)})
        if not perturbed_any:
            problems.append({"outer_fold": f, "problem": f"no perturbation changed the test {target}"})
    return _result(canary, offenders, {"perturbed_runs": runs}, problems)


def check_label_invariance(outer_fold, labels, features, plan_fn,
                           modes=("invert", "permute", "resample"), repeats=2, seed=0):
    """Fail when perturbing an outer fold's test labels changes plan_fn's output.

    plan_fn(outer_fold_id, labels, features) returns the inner folds, tuning
    holdouts and selected config as JSON values; a dict names changed keys.
    plan_fn must rebuild its output from the arrays it is given. A callback
    that returns a cached or persisted plan passes without proving anything.
    """
    return _invariance("label_invariance", outer_fold, labels, features, plan_fn,
                       "labels", modes, repeats, seed, _diff_keys)


def check_feature_invariance(outer_fold, labels, features, plan_fn,
                             modes=("resample", "permute_rows"), repeats=2, seed=0):
    """Fail when perturbing an outer fold's test features changes plan_fn's output.

    plan_fn follows check_label_invariance's contract, including rebuilding
    from its inputs. features is a 2-D samples-by-features array.
    """
    return _invariance("feature_invariance", outer_fold, labels, features, plan_fn,
                       "features", modes, repeats, seed, _diff_keys)


def check_feature_screen_training_only(outer_fold, labels, features, screen_fn,
                                       repeats=2, seed=0):
    """Fail when a feature screen's selection moves with outer test labels or features.

    screen_fn(outer_fold_id, labels, features) returns the selected feature ids
    and must rerun the screen on the arrays it is given, not return a cached
    selection.
    """
    if features is None:
        raise ValueError("feature_screen_training_only perturbs test features, "
                         "so it needs a features array")
    runs = [
        _invariance("feature_screen_training_only", outer_fold, labels, features, screen_fn,
                    "labels", ("invert", "permute", "resample"), repeats, seed, _diff_sets),
        _invariance("feature_screen_training_only", outer_fold, labels, features, screen_fn,
                    "features", ("resample", "permute_rows"), repeats, seed + 1, _diff_sets),
    ]
    return _result("feature_screen_training_only",
                   runs[0].offenders + runs[1].offenders,
                   {"perturbed_runs": sum(r.details["perturbed_runs"] for r in runs)},
                   [p for r in runs for p in r.details.get("problems", [])])


# 5. Common, canonical fold-plan hashes

def canonical_json(value):
    """Sorted keys, no whitespace, ASCII, no NaN; dict keys must be strings."""
    def check(v):
        if isinstance(v, dict):
            for k, item in v.items():
                if not isinstance(k, str):
                    raise TypeError(f"canonical JSON needs string keys, got {k!r}")
                check(item)
        elif isinstance(v, list):
            for item in v:
                check(item)
        elif isinstance(v, (set, frozenset)):
            raise TypeError("canonical JSON has no unordered sets; pass a sorted list")
    value = _plain(value)
    check(value)
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False)


def fold_plan_hash(plan):
    """sha256 hex digest of the plan's canonical JSON, UTF-8 encoded."""
    return hashlib.sha256(canonical_json(plan).encode("utf-8")).hexdigest()


def check_common_fold_plan(arms, require_plan=True):
    """Fail when compared arms use different or noncanonical fold-plan hashes.

    arms maps an arm or learner name to {"plan_hash": str, "plan": object};
    the plan is optional only with require_plan=False.
    """
    offenders, problems = [], []
    if len(arms) < 2:
        problems.append({"problem": "fewer_than_two_arms", "arms": sorted(arms)})
    for arm in sorted(arms):
        entry = arms[arm]
        if "plan" in entry:
            computed = fold_plan_hash(entry["plan"])
            if computed != entry["plan_hash"]:
                offenders.append({"arm": arm, "problem": "noncanonical_hash",
                                  "recorded": entry["plan_hash"], "computed": computed})
        elif require_plan:
            offenders.append({"arm": arm, "problem": "plan_missing"})
    by_hash = defaultdict(list)
    for arm in sorted(arms):
        by_hash[arms[arm]["plan_hash"]].append(arm)
    if len(by_hash) > 1:
        # The most common hash is the reference; ties go to the first arm by name.
        reference = min(by_hash, key=lambda h: (-len(by_hash[h]), by_hash[h][0]))
        for h, names in sorted(by_hash.items()):
            if h != reference:
                offenders.extend({"arm": arm, "problem": "different_plan", "plan_hash": h,
                                  "reference_arm": by_hash[reference][0],
                                  "reference_hash": reference} for arm in names)
    return _result("common_fold_plan", offenders,
                   {"hashes": {h: names for h, names in sorted(by_hash.items())}}, problems)


# 7. One outer fold per genome in a pooled multi-task model

def check_pooled_task_folds(task_folds, groups=None):
    """Fail when a genome, or a group, sits in different outer folds across tasks.

    task_folds maps a task (species-drug pair) to {genome id: outer fold}.
    groups optionally maps genome id to a lineage or duplicate group id.
    """
    folds_of = defaultdict(dict)
    for task in sorted(task_folds, key=_order):
        for genome, fold in task_folds[task].items():
            folds_of[_plain(genome)][task] = _plain(fold)
    offenders = []
    for genome in sorted(folds_of, key=_order):
        by_task = folds_of[genome]
        if len(set(by_task.values())) > 1:
            leaks = [{"outer_fold": f,
                      "test_in": [t for t, g in by_task.items() if g == f],
                      "train_in": [t for t, g in by_task.items() if g != f]}
                     for f in sorted(set(by_task.values()), key=_order)]
            offenders.append({"genome": genome, "problem": "genome_folds_differ",
                              "folds_by_task": [{"task": t, "outer_fold": g}
                                                for t, g in by_task.items()],
                              "leaks": leaks})
    if groups is not None:
        members = defaultdict(list)
        for genome in folds_of:
            gid = groups.get(genome)
            if gid is not None:
                members[_plain(gid)].append(genome)
        for gid in sorted(members, key=_order):
            spread = defaultdict(set)
            for genome in members[gid]:
                for fold in folds_of[genome].values():
                    spread[fold].add(genome)
            if len(spread) > 1:
                offenders.append({"group": gid, "problem": "group_folds_differ",
                                  "folds": [{"outer_fold": f, "genomes": sorted(spread[f], key=_order)}
                                            for f in sorted(spread, key=_order)]})
    return _result("pooled_task_folds", offenders,
                   {"tasks": len(task_folds), "genomes": len(folds_of)})


# 8. Fold metrics that stay honest on one-class folds

@dataclass(frozen=True)
class Metric:
    value: float | None
    reason: str | None = None


def _binary(labels, probs):
    y = np.asarray(labels)
    p = np.asarray(probs, dtype=float)
    if y.shape != p.shape or y.ndim != 1:
        raise ValueError("labels and scores need the same one-dimensional shape")
    if not np.isin(y, (0, 1)).all():
        raise ValueError("labels must be 0 or 1")
    if not np.isfinite(p).all():
        raise ValueError("scores must be finite")
    return y.astype(int), p


def _one_class(y):
    if len(y) == 0:
        return "no samples"
    if y.min() == y.max():
        return f"one class: every label is {int(y[0])}"
    return None


def auroc(labels, scores):
    """Mann-Whitney AUROC with tied scores at half credit; unavailable on one class."""
    y, s = _binary(labels, scores)
    if reason := _one_class(y):
        return Metric(None, reason)
    _, inverse, counts = np.unique(s, return_inverse=True, return_counts=True)
    ranks = (np.cumsum(counts) - (counts - 1) / 2.0)[inverse]
    n_pos = int(y.sum())
    n_neg = len(y) - n_pos
    return Metric(float((ranks[y == 1].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)))


def balanced_accuracy(labels, probs, threshold=0.5):
    """Mean of sensitivity and specificity, calling score >= threshold resistant."""
    y, p = _binary(labels, probs)
    if reason := _one_class(y):
        return Metric(None, reason)
    called = p >= threshold
    return Metric(float((called[y == 1].mean() + (~called[y == 0]).mean()) / 2.0))


def log_loss(labels, probs, eps=1e-15):
    y, p = _binary(labels, probs)
    if len(y) == 0:
        return Metric(None, "no samples")
    if ((p < 0) | (p > 1)).any():
        raise ValueError("probabilities must lie in [0, 1]")
    p = np.clip(p, eps, 1 - eps)
    # fsum: the same terms in any sample order give the same bits (see STAT_TOL_EPS).
    return Metric(-math.fsum((y * np.log(p) + (1 - y) * np.log(1 - p)).tolist()) / len(y))


def brier_score(labels, probs):
    y, p = _binary(labels, probs)
    if len(y) == 0:
        return Metric(None, "no samples")
    if ((p < 0) | (p > 1)).any():
        raise ValueError("probabilities must lie in [0, 1]")
    return Metric(math.fsum(((p - y) ** 2).tolist()) / len(y))


def fold_metrics(labels, probs, threshold=0.5):
    return {"auroc": auroc(labels, probs),
            "balanced_accuracy": balanced_accuracy(labels, probs, threshold),
            "log_loss": log_loss(labels, probs),
            "brier": brier_score(labels, probs)}


def _metric_parts(value):
    if isinstance(value, Metric):
        return value.value, value.reason
    if isinstance(value, dict):
        return value.get("value"), value.get("reason")
    return value, None


def check_one_class_fold_metrics(folds):
    """Fail when a fold's reported metrics fake or drop a value.

    folds maps a fold id to {"labels": per-sample 0/1, "metrics": {name: value}},
    where a value is a Metric, a {"value", "reason"} dict, a number or None.
    On a one-class fold AUROC and balanced accuracy must be None with a reason
    and log-loss and Brier finite; on a two-class fold all four are finite.
    """
    offenders = []
    for f in sorted(folds, key=_order):
        y = _binary_labels(folds[f]["labels"]).astype(int)
        reported = folds[f]["metrics"]
        one_class = _one_class(y)
        for name in ("auroc", "balanced_accuracy", "log_loss", "brier"):
            where = {"fold": _plain(f), "metric": name}
            if name not in reported:
                offenders.append({**where, "problem": "missing"})
                continue
            value, reason = _metric_parts(reported[name])
            unavailable = one_class is not None and (len(y) == 0 or name in ("auroc", "balanced_accuracy"))
            if unavailable and value is not None:
                offenders.append({**where, "problem": "value_on_unavailable_metric",
                                  "value": _plain(value), "why_unavailable": one_class})
            elif unavailable and not reason:
                offenders.append({**where, "problem": "unavailable_without_reason"})
            elif not unavailable and (value is None or not math.isfinite(value)):
                offenders.append({**where, "problem": "not_finite", "value": _plain(value)})
    return _result("one_class_fold_metrics", offenders, {"folds": len(folds)})


# 9. Group-preserving permutation nulls

def _blocks(groups, n, strata):
    """Exchangeable blocks: groups of equal size within one stratum."""
    members = _members(groups, n)
    grouped = {i for idx in members.values() for i in idx}
    units = [idx for idx in members.values() if idx]
    units += [[i] for i in range(n) if i not in grouped]
    if len(grouped) != sum(len(idx) for idx in members.values()):
        raise ValueError("a sample belongs to two groups")
    strata = None if strata is None else _plain(list(strata))
    _check_strata_length(strata, n)
    blocks = defaultdict(list)
    for idx in units:
        stratum = None
        if strata is not None:
            if len({strata[i] for i in idx}) > 1:
                raise ValueError(f"group with samples {idx} spans two strata")
            stratum = strata[idx[0]]
        blocks[(_order(stratum), len(idx))].append(idx)
    return [np.array(block) for _, block in sorted(blocks.items())]


def _check_strata_length(strata, n):
    if strata is not None and len(strata) != n:
        raise ValueError(f"strata has {len(strata)} entries for {n} samples")


def _binary_labels(labels):
    """labels as an array, refusing anything but 0 and 1. A NaN label would make
    two identical groups look like different compositions, since NaN != NaN."""
    labels = np.asarray(labels)
    if not np.isin(labels, (0, 1)).all():
        raise ValueError("labels must be 0 or 1; NaN and other values are refused")
    return labels


def group_permutation(labels, groups, rng, strata=None):
    """Move whole label vectors between groups of the same size and stratum.

    Group sizes, the label pattern inside each group and each stratum's class
    counts stay fixed. Ungrouped samples (None) are singleton groups. labels
    must be 0 or 1.
    """
    labels = _binary_labels(labels)
    return _permute_blocks(labels, _blocks(groups, len(labels), strata), rng)


def _permute_blocks(labels, blocks, rng):
    out = labels.copy()
    for block in blocks:
        if len(block) > 1:
            # argsort of uniform keys is a uniform permutation, and faster than
            # rng.permutation on the small blocks a null draws thousands of times.
            out[block.ravel()] = labels[block[np.argsort(rng.random(len(block)))].ravel()]
    return out


def check_group_level_permutation(labels, permuted, groups, strata=None):
    """Fail when a permutation broke group structure: each block of equal-size
    groups in one stratum must keep its multiset of group label compositions."""
    labels, permuted = _binary_labels(labels), _binary_labels(permuted)
    if len(permuted) != len(labels):
        raise ValueError(f"permuted has {len(permuted)} entries for {len(labels)} labels")
    offenders = []
    for block in _blocks(groups, len(labels), strata):
        before = Counter(tuple(sorted(_plain(labels[idx].tolist()), key=_order)) for idx in block)
        after = Counter(tuple(sorted(_plain(permuted[idx].tolist()), key=_order)) for idx in block)
        if before != after:
            broken = [idx for idx in block
                      if tuple(sorted(_plain(permuted[idx].tolist()), key=_order)) not in before]
            offenders.append({"group_size": int(block.shape[1]),
                              "samples": [_plain(idx) for idx in broken] or _plain(block)})
    return _result("group_level_permutation", offenders)


# Fixed shares above this get a coverage note in details; see check_permutation_power.
COVERAGE_NOTE_SHARE = 0.2


def check_permutation_power(labels, groups, strata=None, n_perm=999, alpha=0.05,
                            min_movable_groups=10, max_fixed_share=0.5):
    """Fail as underpowered when a group-level null cannot test at alpha.

    A block holds the groups of one size in one stratum, and whole label
    vectors move between the groups of a block. Ungrouped samples count as
    singleton groups. A group is movable when its block holds a group with a
    different label composition (its sorted labels), so that a move can change
    some group's class totals. Clonal duplicates labelled (0, 1) and (1, 0)
    swap into a new labelling without moving any class total, so they count as
    fixed. L, the number of distinct reachable labellings, still counts
    ordered label vectors; it is an upper bound on the number of distinct
    statistic values.

    The null tests only the movable samples, whatever the share. Fixed samples
    add the same amount to the observed and to every null statistic, so a
    signal confined to them is invisible. details report fixed_samples,
    fixed_share and fixed_class_counts so the caller can see what is frozen,
    and a share above 20% adds a coverage_note without failing.

    The null is underpowered when any of these holds:
    - fewer than min_movable_groups groups can move. The default floor of 10
      is a judgement, not a derivation: below it the null rests on too few
      exchangeable units to read. Raise it from a power analysis if one exists.
    - fixed groups hold more than max_fixed_share of the samples. The default
      of 0.5 is also a judgement: a null that never moves most of the data is
      not a test of that data. Heavy-tailed lineage sizes leave the largest
      lineages alone in their size block, so clean GRM-like layouts commonly
      fix 20-50% of their samples.
    - the expected p-value of the most extreme labelling,
      (1 + n_perm / L) / (n_perm + 1), exceeds alpha. That covers too few
      labellings and too few permutations.

    labels must be 0 or 1, and strata, if given, needs one entry per sample.

    permutation_null adds checks on the draws it makes: degenerate,
    unreachable and undefined_statistic.
    """
    labels = _binary_labels(labels)
    return _power(labels, _blocks(groups, len(labels), strata), n_perm, alpha,
                  min_movable_groups, max_fixed_share)


def _power(labels, blocks, n_perm, alpha, min_movable_groups, max_fixed_share):
    movable = fixed = movable_samples = fixed_samples = 0
    fixed_classes = Counter()
    labellings = 1
    for block in blocks:
        patterns = Counter(tuple(_plain(labels[idx].tolist())) for idx in block)
        compositions = {tuple(sorted(pattern, key=_order)) for pattern in patterns}
        if len(compositions) > 1:
            movable += len(block)
            movable_samples += block.size
        else:
            fixed += len(block)
            fixed_samples += block.size
            fixed_classes.update(_plain(labels[block.ravel()].tolist()))
        left = len(block)
        for count in patterns.values():
            labellings *= math.comb(left, count)
            left -= count
    total = movable_samples + fixed_samples
    fixed_share = fixed_samples / total if total else 0.0
    expected_min_p = (1 + n_perm / labellings) / (n_perm + 1)
    details = {"movable_groups": movable, "fixed_groups": fixed,
               "movable_samples": movable_samples, "fixed_samples": fixed_samples,
               "fixed_share": fixed_share,
               "fixed_class_counts": {str(k): fixed_classes[k] for k in sorted(fixed_classes, key=_order)},
               "distinct_labellings": labellings if labellings < 10**15 else None,
               "distinct_labellings_log10": math.log10(labellings),
               "expected_min_p": expected_min_p, "n_perm": n_perm, "alpha": alpha,
               "min_movable_groups": min_movable_groups, "max_fixed_share": max_fixed_share}
    if COVERAGE_NOTE_SHARE < fixed_share <= max_fixed_share:
        details["coverage_note"] = (f"fixed groups hold {fixed_samples} of {total} samples "
                                    f"({fixed_share:.0%}); the null tests only the other "
                                    f"{movable_samples}")
    problems = []
    if movable < min_movable_groups:
        problems.append({"problem": "underpowered", "cause": "floor",
                         "reason": f"{movable} groups can move; the floor is {min_movable_groups}"})
    if fixed_share > max_fixed_share:
        problems.append({"problem": "underpowered", "cause": "fixed_share",
                         "reason": f"fixed groups hold {fixed_samples} of {total} samples "
                                   f"({fixed_share:.0%}); the limit is {max_fixed_share:.0%}"})
    if expected_min_p > alpha:
        problems.append({"problem": "underpowered", "cause": "expected_min_p",
                         "reason": f"the smallest expected p-value, {expected_min_p:.4g}, "
                                   f"exceeds alpha {alpha}"})
    return _result("permutation_power", [], details, problems)


# A null whose draws tie the observed statistic this often cannot rank it.
MAX_TIE_SHARE = 0.5

# Statistics are compared with a relative tolerance, as scipy.stats.permutation_test
# does: STAT_TOL_EPS machine epsilons of the statistic's own dtype, so 100 float64
# epsilons (STAT_RTOL) for a float64 statistic and 100 float32 epsilons for a float32
# one. A draw within rtol * max(|reference|, scale) of a reference statistic counts
# as equal to it; scale defaults to 0.
#
# This absorbs rounding noise of that size and no more. Rounding noise follows the
# size of the terms a statistic sums, not the size of its result: a difference of
# two means over features near 1e3 is small while its noise is about 1e3 * eps, so
# ties that hold in exact arithmetic split again. Statistics should therefore sum
# in an order-independent way (math.fsum, integer counts), as mean_difference
# does, or the caller passes scale, the magnitude of the summed terms. A float32
# statistic carries float32's epsilon, 1.2e-7, so its scale-widened tolerance can
# cover the whole null for large uncentred features (every draw then ties, which is
# conservative but detects nothing); centre float32 features or compute in float64.
STAT_TOL_EPS = 100
STAT_RTOL = STAT_TOL_EPS * np.finfo(np.float64).eps


def _stat_rtol(dtype):
    """STAT_TOL_EPS epsilons of an inexact dtype; 0 for exact types such as ints."""
    dtype = np.dtype(dtype)
    return STAT_TOL_EPS * float(np.finfo(dtype).eps) if np.issubdtype(dtype, np.inexact) else 0.0


def _gamma(reference, rtol, scale):
    """The tolerance for comparing a statistic with reference; see STAT_TOL_EPS."""
    return rtol * max(abs(reference), scale)


def _at_least(null, reference, rtol, scale):
    """Draws at least as extreme as reference, counting near-ties as ties."""
    return null >= reference - _gamma(reference, rtol, scale)


# Statistic types the tolerance can handle. Each casts to float64 without changing
# value (integers up to 2**53; larger ones can round together, which only merges
# values), so comparing in float64 with the tolerance of the returned dtype is
# sound. np.longdouble is refused: its epsilon is below one float64 ulp, so the
# cast moves values by more than its tolerance. Complex, Decimal, Fraction,
# np.timedelta64 (a numpy integer subclass) and other objects are refused too.
# Subclasses are accepted and take their base type's dtype: a Python float subclass
# is float64 and a Python int subclass (bool, IntEnum) is int64.
STAT_TYPES = (int, float, np.integer, np.float16, np.float32, np.float64)


def _value_dtype(value):
    """The dtype of one accepted statistic value, or ValueError naming its type."""
    if not isinstance(value, STAT_TYPES) or isinstance(value, np.timedelta64):
        raise ValueError(f"statistic returned {type(value).__name__}; it must return a "
                         f"Python int or float, a numpy integer, or numpy float16, "
                         f"float32 or float64 (use float() or .item())")
    if isinstance(value, np.generic):
        return value.dtype
    if isinstance(value, int):
        try:
            float(value)
        except OverflowError:
            raise ValueError(f"statistic returned an int too large for float64 "
                             f"(about 10**{value.bit_length() * 0.30103:.0f})") from None
    return np.dtype(np.float64) if isinstance(value, float) else np.dtype(np.int64)


def _stat_dtype(values):
    """The dtype that sets the tolerance; refuse a type it cannot handle."""
    return np.result_type(*{_value_dtype(value) for value in values})


def _scale(scale):
    """scale as a plain float, or 0.0 for None; refuse bools and non-numbers."""
    if scale is None:
        return 0.0
    # np.bool_ is not numbers.Real, so the second test refuses it.
    if isinstance(scale, bool) or not isinstance(scale, numbers.Real):
        raise ValueError(f"scale must be a finite number >= 0, not {type(scale).__name__}")
    try:
        scale = float(scale)
    except OverflowError:
        scale = math.inf
    if not (math.isfinite(scale) and scale >= 0):
        raise ValueError("scale must be a finite number >= 0")
    return scale


def mean_difference(labels, feature):
    """Mean feature of class 1 minus mean feature of class 0.

    Sums use math.fsum, which rounds the exact sum once. Any two labellings that
    put the same multiset of values in each class, such as swapped labels on a
    pair of clonal duplicates anywhere in the sample order, give bit-identical
    results. A numpy sum or dot product does not guarantee that.
    """
    labels = np.asarray(labels)
    feature = np.asarray(feature, dtype=float)
    positive = labels == 1
    n_pos = int(positive.sum())
    return (math.fsum(feature[positive].tolist()) / n_pos
            - math.fsum(feature[~positive].tolist()) / (len(labels) - n_pos))


def _null_power(group_power, null, observed, alpha, rtol=STAT_RTOL, scale=0.0):
    """Add the draws' own evidence to the group power check; return it with p.

    group_power is None when a custom permute replaced group_permutation; the
    group counts would then describe a null that was never drawn.
    """
    if group_power is None:
        details = {"group_power": "not applicable",
                   "group_power_reason": "a custom permute replaced group_permutation, so "
                                         "movable groups and labellings do not describe this null"}
        problems = []
    else:
        details = dict(group_power.details)
        problems = list(details.pop("problems", []))
    non_finite = int(np.sum(~np.isfinite(null)))
    if non_finite or not math.isfinite(observed):
        # A NaN compares False both ways, so it would read as the strongest leak
        # (observed) or as a draw below observed. Neither a leak nor a clean null.
        details.update({"observed_finite": math.isfinite(observed), "non_finite_draws": non_finite,
                        "ties_with_observed": None, "tie_share": None, "null_variance": None,
                        "best_attainable_p": None})
        problems.append({"problem": "undefined_statistic",
                         "reason": f"the observed statistic is {observed} and {non_finite} of "
                                   f"{len(null)} draws are NaN or infinite"})
        return _result("permutation_power", [], details, problems), math.nan
    ties = int(np.sum(np.abs(null - observed) <= _gamma(observed, rtol, scale)))
    # The observed labelling is itself reachable, so no labelling could beat the
    # top of null and observed together: this is the smallest p these draws allow.
    top = max(observed, float(null.max()))
    best_p = float((1 + np.sum(_at_least(null, top, rtol, scale))) / (len(null) + 1))
    with np.errstate(over="ignore", invalid="ignore"):
        # Finite draws near +-1e308 overflow; the range is then inf, not constant.
        variance = float(np.var(null))
        spread = float(null.max() - null.min())
    details.update({"ties_with_observed": ties, "tie_share": ties / len(null),
                    "null_variance": variance if math.isfinite(variance) else None,
                    "best_attainable_p": best_p, "stat_rtol": rtol, "stat_scale": scale})
    if ties / len(null) >= MAX_TIE_SHARE:
        problems.append({"problem": "degenerate",
                         "reason": f"{ties} of {len(null)} draws tie the observed statistic"})
    elif spread <= _gamma(max(abs(null.max()), abs(null.min())), rtol, scale):
        problems.append({"problem": "degenerate",
                         "reason": "every draw gives the same statistic"})
    if best_p > alpha:
        problems.append({"problem": "unreachable",
                         "reason": f"the best p-value any labelling could get from these draws, "
                                   f"{best_p:.4g}, exceeds alpha {alpha}"})
    p_value = float((1 + np.sum(_at_least(null, observed, rtol, scale))) / (len(null) + 1))
    return _result("permutation_power", [], details, problems), p_value


@dataclass(frozen=True)
class NullResult:
    observed: float
    null: np.ndarray
    p_value: float
    ties_with_observed: int | None
    power: CanaryResult


def permutation_null(statistic, labels, groups, n_perm=999, seed=0, strata=None,
                     permute=group_permutation, alpha=0.05, min_movable_groups=10,
                     max_fixed_share=0.5, scale=None):
    """One-sided permutation p-value, (1 + #{null >= observed}) / (n_perm + 1).

    statistic(labels) returns a number where larger means more signal. Ties
    with the observed statistic count against it, so the p-value is valid
    but conservative; ties_with_observed reports how many there were. labels
    must be 0 or 1.

    Every comparison between statistics (ties, the p-value, the best
    attainable p and the constant-null rule) treats values within
    rtol * max(|reference|, scale) of each other as equal. rtol is
    STAT_TOL_EPS epsilons of the dtype the statistic returns, so return a
    float32 result as np.float32, not cast to float; it is 0 for integers.
    The statistic must return one of STAT_TYPES (Python int or float, a numpy
    integer, numpy float16, float32 or float64); np.longdouble, complex,
    Decimal, Fraction and other objects raise ValueError. scale is a real
    number >= 0, not a bool, and details report it as a float.
    The tolerance covers rounding noise of that size only. Noise follows the
    size of the summed terms, not of the result, so the statistic must be
    order-independent: two labellings that are equal in exact arithmetic must
    give bit-identical values. Sum with math.fsum or count integers, as
    mean_difference does. A statistic that sums in floating point instead
    must pass scale, the magnitude of the terms it sums (for example max
    |feature| times the number of samples); otherwise clonal duplicates far
    apart in sample order split exact ties, and both the p-value and the
    power check can read too strong.

    Read power.passed before trusting a large p-value. Its problems name one
    of four failures:
    - underpowered: the group null fails check_permutation_power. Its cause
      key says which rule: floor, fixed_share or expected_min_p.
    - degenerate: at least half the draws tie the observed statistic, or
      every draw gives the same value. The constant-null rule also marks an
      indicator statistic, nonzero only at the observed labelling, although
      its p-value is exact; that errs conservative and needs a contrived
      statistic.
    - unreachable: the smallest p-value this test's own draws can produce,
      (1 + #{null >= max(null and observed)}) / (n_perm + 1), exceeds alpha.
      That is the p of a labelling that ties the largest observed or drawn
      statistic; a labelling strictly above every draw would get less.
      It does not prove that no labelling could reach alpha: when the draws
      miss the rarest top values, another set of draws may reach it. The
      third re-check found 15 such nulls among 105 that can reach alpha 0.05,
      at n_perm 199. The error is conservative.
    - undefined_statistic: the observed statistic or some draw is NaN or
      infinite. p_value is then NaN, so it reads neither as a leak nor as
      clean, and details give non_finite_draws and observed_finite.
    With a custom permute the group checks are reported as not applicable
    and only the draw checks run.
    """
    if n_perm < 1:
        raise ValueError("permutation_null needs n_perm >= 1")
    scale = _scale(scale)
    labels = _binary_labels(labels)
    # Checked here as well as in _blocks, because a custom permute skips _blocks.
    _check_strata_length(strata, len(labels))
    rng = np.random.default_rng(seed)
    if permute is group_permutation:
        blocks = _blocks(groups, len(labels), strata)
        group_power = _power(labels, blocks, n_perm, alpha, min_movable_groups, max_fixed_share)

        def draw():
            return _permute_blocks(labels, blocks, rng)
    else:
        group_power = None
        counts = np.bincount(labels.ravel().astype(int), minlength=2).tolist()

        def draw():
            # A custom permuter is the caller's code: every draw must be a true
            # permutation of labels, with its shape, 0/1 values and class counts.
            # A resampler such as rng.choice(y, size=len(y)) changes the counts,
            # and a sum statistic then compares the observed value with draws
            # that hold more or fewer positives.
            drawn = _binary_labels(permute(labels, groups, rng, strata))
            if len(drawn) != len(labels):
                raise ValueError(f"permute returned {len(drawn)} labels for "
                                 f"{len(labels)} samples")
            if drawn.shape != labels.shape:
                raise ValueError(f"permute returned shape {drawn.shape} for labels of "
                                 f"shape {labels.shape}")
            drawn_counts = np.bincount(drawn.ravel().astype(int), minlength=2).tolist()
            if drawn_counts != counts:
                raise ValueError(f"permute changed the class counts from {counts} to "
                                 f"{drawn_counts}; a permutation keeps them")
            return drawn
    raw_observed = statistic(labels)
    _value_dtype(raw_observed)  # fail before n_perm more statistic calls
    raw_null = [statistic(draw()) for _ in range(n_perm)]
    rtol = _stat_rtol(_stat_dtype([raw_observed] + raw_null))
    observed = float(raw_observed)
    null = np.array(raw_null, dtype=float)
    power, p_value = _null_power(group_power, null, observed, alpha, rtol, scale)
    return NullResult(observed, null, p_value, power.details["ties_with_observed"], power)


def check_pvalue_uniformity(pvalues, alpha=0.01, bins=10):
    """Fail when null p-values hold too many small values (anti-conservative).

    Leakage pushes a null's p-values toward 0, lifting their empirical CDF
    above the uniform one. D+ = max(i/n - p_(i)) measures that and decides the
    verdict against the one-sided asymptotic Kolmogorov-Smirnov critical value
    sqrt(-ln(alpha) / 2) / sqrt(n). D- = max(p_(i) - (i-1)/n) measures
    p-values piled toward 1, which a tied or discrete statistic produces
    validly; details report it as a conservativeness note, not a failure.
    A NaN p-value, from an undefined statistic, fails the check rather than
    being dropped. With so few p-values that the cutoff is at least 1, no
    sample could fail (D+ < 1), so the check fails as underpowered.
    """
    p = np.sort(np.asarray(pvalues, dtype=float))
    n = len(p)
    if n == 0:
        return _result("pvalue_uniformity", [], {}, [{"problem": "no p-values"}])
    undefined = int(np.sum(~np.isfinite(p)))
    if undefined:
        return _result("pvalue_uniformity", [], {"n": n, "undefined_pvalues": undefined},
                       [{"problem": "undefined_statistic",
                         "reason": f"{undefined} of {n} p-values are NaN or infinite"}])
    i = np.arange(1, n + 1)
    d_plus = float(max(0.0, (i / n - p).max()))
    d_minus = float(max(0.0, (p - (i - 1) / n).max()))
    critical = math.sqrt(-0.5 * math.log(alpha)) / math.sqrt(n)
    counts, _ = np.histogram(p, bins=bins, range=(0.0, 1.0))
    details = {"n": n, "d_plus": d_plus, "d_minus": d_minus, "critical": critical,
               "alpha": alpha, "conservative": d_minus > critical,
               "bin_counts": counts.tolist(), "share_below_0.05": float(np.mean(p <= 0.05))}
    offenders = [] if d_plus <= critical else [{"problem": "pvalues_anti_conservative",
                                                "d_plus": d_plus, "critical": critical}]
    problems = []
    if critical >= 1:
        problems.append({"problem": "underpowered", "cause": "too_few_pvalues",
                         "reason": f"the cutoff for {n} p-values, {critical:.4g}, is at least 1, "
                                   f"so no set of p-values could fail"})
    return _result("pvalue_uniformity", offenders, details, problems)


def min_calibration_sims(n_perm, alpha=0.01):
    """The smallest n_sim at which null_calibration can fail.

    The smallest p a null can give is 1/(n_perm + 1), so D+ is at most
    1 - 1/(n_perm + 1), reached when every p-value sits at that floor. The
    check fails when D+ exceeds sqrt(-ln(alpha) / 2) / sqrt(n_sim), so it can
    fail only when n_sim > -ln(alpha) / 2 / (1 - 1/(n_perm + 1))**2. At the
    defaults, alpha 0.01 and n_perm 99, the minimum is 3.
    """
    top = 1 - 1 / (n_perm + 1)
    return math.floor(-0.5 * math.log(alpha) / top ** 2) + 1


def null_calibration(simulate, statistic, n_sim=200, n_perm=99, seed=0, alpha=0.01,
                     permute=group_permutation):
    """Simulate true-null datasets and test that their p-values are uniform.

    simulate(rng) returns (labels, groups, strata, data) with labels drawn
    independently of data; statistic(labels, data) returns a float.

    The default n_sim=200 keeps the tests fast but is weak against a mild
    leak: a sample-level shuffle on the tests' clustered data is missed in
    about 4% of seeds at 200. Live calibration runs should use more
    simulations, for example 1000, where the D+ cutoff shrinks by sqrt(5).

    n_sim below min_calibration_sims(n_perm, alpha), 3 at the defaults,
    raises ValueError: with fewer simulations even a null whose every p-value
    sits at its floor would pass.
    """
    if isinstance(n_sim, bool) or not isinstance(n_sim, numbers.Integral):
        raise ValueError(f"n_sim must be an integer, not {type(n_sim).__name__}")
    if isinstance(n_perm, bool) or not isinstance(n_perm, numbers.Integral) or n_perm < 1:
        raise ValueError("null_calibration needs an integer n_perm >= 1")
    if not 0 < alpha < 1:
        raise ValueError("null_calibration needs 0 < alpha < 1")
    needed = min_calibration_sims(n_perm, alpha)
    if n_sim < needed:
        raise ValueError(f"null_calibration needs n_sim >= {needed} at n_perm {n_perm} and "
                         f"alpha {alpha}; with {n_sim} it could never fail")
    rng = np.random.default_rng(seed)
    pvalues = []
    for _ in range(n_sim):
        labels, groups, strata, data = simulate(rng)
        result = permutation_null(lambda y: statistic(y, data), labels, groups, n_perm,
                                  int(rng.integers(2**32)), strata, permute)
        pvalues.append(result.p_value)
    checked = check_pvalue_uniformity(pvalues, alpha)
    return _result("null_calibration", checked.offenders,
                   {**checked.details, "n_perm": n_perm}, checked.details.get("problems", ()))
