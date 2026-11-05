"""Genome-wide development folds shared by every GRM drug task.

A pooled multi-task model trains on all drug tasks at once, so a genome must
sit in the same fold for every task: global_folds gives each genome one fold.
Duplicate groups stay whole, and under the lineage scheme so do lineages,
including lineages a duplicate group links. rolling_origin_dev_folds gives
time-ordered folds that keep duplicate groups whole.

Every split function takes a required lockbox argument, a sealed manifest
from lockbox.define_lockbox or an explicit None, and refuses any genome the
lockbox claims. Pass expected_seal, the seal recorded at seal time, to refuse
a resealed or different lockbox, and releases, every earlier release used
for development since the lockbox was defined, so that duplicate closure
follows chains across releases. Remove claimed genomes first by running
lockbox.development_genomes on the full release, so that duplicate closure
sees every claimed partner. A split's own check sees only the genomes it is
given. Identifiers follow rule 6 as lockbox.py states it.

Eligibility reports every task in every fold. A task-fold whose test or
train side has fewer than min_test labelled genomes, or fewer than min_class
of either class for a 0/1-labelled task, is insufficient_evidence; it is
never dropped. Both minimums must be at least 1. Each manifest's fold_summary
gives every fold's task statuses; fit a task on a fold only if
eligible_tasks(record) names it (for list-form tasks, see its caveat). The global fold_summary also gives each
fold's largest unit and its share of all genomes, so one lineage holding a
large share is visible.
"""
import importlib
import importlib.util
import random
import sys
from collections import Counter
from pathlib import Path


def _sibling(name):
    if __package__:
        return importlib.import_module("." + name, __package__)
    key = "grm_multitask_" + name
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).with_name(name + ".py"))
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


lb = _sibling("lockbox")

KIND = "grm.multitask.splits"
VERSION = 2
SCHEMES = ("random", "lineage")
# Placeholders only: set both per protocol for each task.
DEFAULT_MIN_TEST = 10
DEFAULT_MIN_CLASS = 2
STRATIFICATION = (
    "Units (duplicate components; under the lineage scheme, lineage-and-duplicate components) go largest "
    "first, equal sizes in seeded order, each to the fold that least raises the sum over strata of the squared "
    "fold count divided by the stratum's even share, ties in seeded fold order. Strata are all genomes, each "
    "task's labelled genomes and, for 0/1-labelled tasks, each class. Best effort; balance measures the result.")
ROLLING_RULE = (
    "A duplicate component takes the latest possible collection day of its dated members. Fold i trains on "
    "components on or before cutoffs[i] and tests on those after it, up to and including cutoffs[i+1]; the "
    "last fold tests on everything after the last cutoff. Components with no dated member are left out and "
    "listed as undated.")


def _development(genomes, lockbox, expected_seal, releases):
    rows = lb.normalize_genomes(genomes)
    if lockbox is None and (expected_seal is not None or releases):
        raise lb.LockboxError("expected_seal or releases given without a lockbox")
    seal = None if lockbox is None else lb.require_development(rows, lockbox, expected_seal, releases)
    return rows, seal


def units(genomes, scheme):
    """Sorted genome-id lists that must share a fold, ordered by first id."""
    if scheme not in SCHEMES:
        raise ValueError("scheme must be one of %s, got %r" % (", ".join(SCHEMES), scheme))
    rows = lb.normalize_genomes(genomes)
    parent = {row["genome_id"]: row["genome_id"] for row in rows}

    def find(gid):
        while parent[gid] != gid:
            parent[gid] = parent[parent[gid]]
            gid = parent[gid]
        return gid

    first = {}
    for row in rows:
        keys = [("duplicate", row["duplicate_group"])] if row["duplicate_group"] else []
        if scheme == "lineage":
            keys.append(("lineage", row["lineage"]))
        for key in keys:
            if key in first:
                parent[find(row["genome_id"])] = find(first[key])
            else:
                first[key] = row["genome_id"]
    groups = {}
    for row in rows:
        groups.setdefault(find(row["genome_id"]), []).append(row["genome_id"])
    return sorted(groups.values())


def _strata(row):
    yield "genomes"
    for task, label in row["tasks"].items():
        yield "task=" + task
        if label is not None:
            yield "task=%s,class=%d" % (task, label)


def _assign(rows, k, seed, scheme):
    if type(k) is not int or k < 2:
        raise ValueError("k must be an integer of at least 2, got %r" % (k,))
    if type(seed) is not int:
        raise ValueError("seed must be an integer, got %r" % (seed,))
    groups = units(rows, scheme)
    if len(groups) < k:
        raise ValueError("%d %s units cannot fill %d folds" % (len(groups), scheme, k))
    by_id = {row["genome_id"]: row for row in rows}
    profiles = [Counter(s for gid in group for s in _strata(by_id[gid])) for group in groups]
    totals = sum(profiles, Counter())
    weight = {stratum: (k / total) ** 2 for stratum, total in totals.items()}
    rng = random.Random(seed)
    order = list(range(len(groups)))
    rng.shuffle(order)
    order.sort(key=lambda i: -len(groups[i]))
    fold_order = list(range(k))
    rng.shuffle(fold_order)
    load = [Counter() for _ in range(k)]
    assignment = {}
    for i in order:
        profile = profiles[i]
        best = min(fold_order, key=lambda f: sum(weight[s] * n * (2 * load[f][s] + n) for s, n in profile.items()))
        load[best].update(profile)
        assignment.update(dict.fromkeys(groups[i], best))
    return dict(sorted(assignment.items()))


def global_folds(genomes, k, seed, scheme, *, lockbox, expected_seal=None, releases=()):
    """{genome_id: fold} with one fold per genome, shared by every task."""
    rows, _ = _development(genomes, lockbox, expected_seal, releases)
    return _assign(rows, k, seed, scheme)


def fold_sets(assignment):
    """[{"fold", "train", "test"}] for a {genome_id: fold} assignment."""
    k = max(assignment.values()) + 1 if assignment else 0
    folds = [{"fold": f, "train": [], "test": []} for f in range(k)]
    for gid, fold in sorted(assignment.items()):
        for f in range(k):
            folds[f]["test" if f == fold else "train"].append(gid)
    return folds


def _minimum(value, name):
    if type(value) is not int or value < 1:
        raise ValueError("%s must be an integer of at least 1, got %r" % (name, value))
    return value


def _minimums(min_test, min_class, classed):
    """{task: min_test} and min_class, validated."""
    _minimum(min_class, "min_class")
    if not isinstance(min_test, dict):
        return dict.fromkeys(classed, _minimum(min_test, "min_test"))
    missing = sorted(set(classed) - set(min_test))
    if missing:
        raise ValueError("min_test has no minimum for task %s" % ", ".join(missing))
    return {task: _minimum(min_test[task], "min_test[%r]" % task) for task in classed}


def _eligibility(rows, folds, classed, min_test, min_class):
    minimums = _minimums(min_test, min_class, classed)
    by_id = {row["genome_id"]: row for row in rows}
    report = {task: {"status": "ok", "insufficient_folds": [], "folds": []} for task in classed}
    for fold in folds:
        try:
            sides = {side: lb.task_counts([by_id[gid] for gid in fold[side]], classed) for side in ("train", "test")}
        except KeyError as error:
            raise ValueError("fold %r names genome %s that is not in genomes" % (fold["fold"], error)) from None
        for task, entry in report.items():
            reasons = [reason for side in ("test", "train")
                       for reason in lb.evidence_reasons(sides[side][task], minimums[task], min_class, side)]
            entry["folds"].append({"fold": fold["fold"], "train": sides["train"][task], "test": sides["test"][task],
                                   "status": "insufficient_evidence" if reasons else "ok", "reasons": reasons})
            if reasons:
                entry["status"] = "insufficient_evidence"
                entry["insufficient_folds"].append(fold["fold"])
    return report


def eligibility_counts(genomes, folds, *, tasks=(), min_test=DEFAULT_MIN_TEST, min_class=DEFAULT_MIN_CLASS):
    """Per task and fold: train and test counts, status and reasons.

    folds is a {genome_id: fold} assignment or a list of {"fold", "train", "test"}.
    tasks names tasks to report even when no genome is labelled for them.
    """
    rows = lb.normalize_genomes(genomes)
    folds = fold_sets(folds) if isinstance(folds, dict) else folds
    return _eligibility(rows, folds, lb.task_classes(rows, tasks), min_test, min_class)


def balance(genomes, assignment):
    """Per stratum: fold counts, the even share and the largest relative deviation from it."""
    rows = lb.normalize_genomes(genomes)
    k = max(assignment.values()) + 1
    counts = {}
    for row in rows:
        for stratum in _strata(row):
            counts.setdefault(stratum, [0] * k)[assignment[row["genome_id"]]] += 1
    report = {}
    for stratum, per_fold in sorted(counts.items()):
        share = sum(per_fold) / k
        report[stratum] = {"counts": per_fold, "even_share": round(share, 6),
                           "max_relative_deviation": round(max(abs(n - share) for n in per_fold) / share, 6)}
    return report


def _fold_summary(eligibility, n_folds):
    return [{"fold": f, "tasks": {task: entry["folds"][f]["status"] for task, entry in eligibility.items()}}
            for f in range(n_folds)]


def eligible_tasks(fold_record):
    """The tasks a fold_summary record marks ok, sorted; fit only these on that fold.

    Necessary, not sufficient, for a list-form or all-None task: P holds no labels for it, so ok means
    only that each side has at least min_test genomes listing it. P cannot check that training holds
    both classes, which a fit such as Task O's requires; the caller must check that with its own labels.
    """
    tasks = fold_record.get("tasks") if isinstance(fold_record, dict) else None
    if not isinstance(tasks, dict):
        raise ValueError("expected a fold_summary record with a tasks mapping")
    return sorted(task for task, status in tasks.items() if status == "ok")


def _sealed(body):
    return dict(body, sha256=lb.sha256_json(body))


def verify_manifest(manifest, expected_seal=None):
    """The manifest's sha256, after checking it; ValueError on any edit.

    expected_seal, the lockbox seal recorded at seal time, refuses a manifest built against any
    other lockbox, including a resealed forgery, with LockboxError.
    """
    if not isinstance(manifest, dict) or manifest.get("kind") != KIND:
        raise ValueError("not a GRM split manifest")
    body = {key: value for key, value in manifest.items() if key != "sha256"}
    if manifest.get("sha256") != lb.sha256_json(body):
        raise ValueError("split manifest sha256 does not match its contents")
    if expected_seal is not None and manifest.get("lockbox_seal") != expected_seal:
        raise lb.LockboxError("split manifest was built against lockbox seal %r, not the expected seal %.12s"
                              % (manifest.get("lockbox_seal"), expected_seal))
    return manifest["sha256"]


def global_manifest(genomes, k, seed, scheme, *, lockbox, expected_seal=None, releases=(), tasks=(),
                    min_test=DEFAULT_MIN_TEST, min_class=DEFAULT_MIN_CLASS):
    """global_folds with its parameters, balance, eligibility, fold_summary and sha256."""
    rows, seal = _development(genomes, lockbox, expected_seal, releases)
    assignment = _assign(rows, k, seed, scheme)
    classed = lb.task_classes(rows, tasks)
    eligibility = _eligibility(rows, fold_sets(assignment), classed, min_test, min_class)
    groups = units(rows, scheme)
    summary = _fold_summary(eligibility, k)
    for record in summary:
        sizes = [len(group) for group in groups if assignment[group[0]] == record["fold"]]
        record.update(genomes=sum(sizes), largest_unit=max(sizes),
                      largest_unit_share=round(max(sizes) / len(rows), 6))
    return _sealed({
        "kind": KIND, "version": VERSION, "method": "global_folds", "scheme": scheme, "k": k, "seed": seed,
        "stratification": STRATIFICATION, "lockbox_seal": seal, "input_sha256": lb.sha256_json(rows),
        "n_genomes": len(rows), "n_units": len(groups), "min_test": min_test, "min_class": min_class,
        "folds": assignment, "balance": balance(rows, assignment), "eligibility": eligibility,
        "fold_summary": summary,
    })


def _rolling(rows, cutoffs):
    if isinstance(cutoffs, (str, bytes)):
        raise ValueError("cutoffs must be a list of YYYY-MM-DD dates")
    days = [lb.full_date(cutoff) for cutoff in cutoffs]
    if not days or any(a >= b for a, b in zip(days, days[1:])):
        raise ValueError("cutoffs must be a non-empty, strictly increasing list")
    by_id = {row["genome_id"]: row for row in rows}
    when, undated = {}, []
    for group in units(rows, "random"):
        dated = [lb.latest_day(by_id[gid]["collection_date"]) for gid in group if by_id[gid]["collection_date"]]
        if dated:
            when.update(dict.fromkeys(group, max(dated)))
        else:
            undated.extend(group)
    folds = []
    for i, cutoff in enumerate(days):
        until = days[i + 1] if i + 1 < len(days) else None
        folds.append({
            "fold": i, "cutoff": cutoff.isoformat(), "test_until": until and until.isoformat(),
            "train": sorted(gid for gid, day in when.items() if day <= cutoff),
            "test": sorted(gid for gid, day in when.items() if cutoff < day and (until is None or day <= until)),
        })
    return folds, sorted(undated)


def rolling_origin_dev_folds(genomes, cutoffs, *, lockbox, expected_seal=None, releases=()):
    """Time-ordered folds [{"fold", "cutoff", "test_until", "train", "test"}]; see ROLLING_RULE."""
    rows, _ = _development(genomes, lockbox, expected_seal, releases)
    return _rolling(rows, cutoffs)[0]


def rolling_origin_manifest(genomes, cutoffs, *, lockbox, expected_seal=None, releases=(), tasks=(),
                            min_test=DEFAULT_MIN_TEST, min_class=DEFAULT_MIN_CLASS):
    """rolling_origin_dev_folds with the rule, undated genomes, eligibility, fold_summary and sha256."""
    rows, seal = _development(genomes, lockbox, expected_seal, releases)
    folds, undated = _rolling(rows, cutoffs)
    eligibility = _eligibility(rows, folds, lb.task_classes(rows, tasks), min_test, min_class)
    return _sealed({
        "kind": KIND, "version": VERSION, "method": "rolling_origin_dev_folds", "rule": ROLLING_RULE,
        "cutoffs": [fold["cutoff"] for fold in folds], "lockbox_seal": seal, "input_sha256": lb.sha256_json(rows),
        "n_genomes": len(rows), "min_test": min_test, "min_class": min_class, "undated": undated, "folds": folds,
        "eligibility": eligibility, "fold_summary": _fold_summary(eligibility, len(folds)),
    })
