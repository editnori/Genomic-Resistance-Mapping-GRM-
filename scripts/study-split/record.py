#!/usr/bin/env python3
"""Build the record of the two experiments with folds by study at revision 14.

    python3 -I scripts/study-split/record.py --inputs DIR --out reference/study-split
    python3 -I scripts/study-split/record.py --inputs DIR --check reference/study-split

DIR is the reduced input set: per fit and genome the fold, the label and the
five scores; the fold files; the two design files, plans, hash lists and run
logs; and the release fields the record needs (labels, BioProjects, lineage
clusters, copy groups, the lineage split's scores). The script reads nothing
else, fits nothing and opens no network connection.

Every measure is computed again from scores and labels. Recall, specificity,
balanced accuracy and the error rate are exact fractions of the confusion
counts; log loss is the exactly rounded sum of the clipped per-genome losses
over the genome count. Differences, means over draws and medians are exact
over those values. Each species interval is exact: every multiset of the
species drawn with replacement, at the weight of its multinomial count.

The record holds no genome id, label, per-genome prediction, model,
feature or importance; a fold plan enters as its sha256 and its fold counts only.
"""
import argparse
import hashlib
import json
import math
import os
import re
import statistics
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
KINDS = ("study", "random", "matched")
DRAWS = range(5)
FOLDS = 5
CUTOFF = 0.5
CLIP = 1e-4
MEASURES = ("recall", "specificity", "balanced_accuracy", "log_loss", "error_rate")
WORSE_UP = {"recall": False, "specificity": False, "balanced_accuracy": False, "log_loss": True,
            "error_rate": True}
SECOND = {"SR": ("study", "random"), "SM": ("study", "matched"), "MR": ("matched", "random")}
FIRST = {"study": ("study-d0", "random-d0"), "lineage": ("lineage", "random-d0"),
         "study-lineage": ("study-d0", "lineage")}
THRESHOLD = {"SR": 20, "SM": 4, "MR": 20}
EXPERIMENT_LEARNERS = 3
NO_STUDY_MIN = 10
RULE_LINE = "0.10"
RULES = ("either", "both", "recall", "loss", "balanced")


class Stop(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Stop(msg)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def read_tsv(path, columns=None, name=None):
    """The rows of a TSV as dicts. With columns, the header must hold those columns, in that order, and no other."""
    name = name or path
    text = Path(path).read_text()
    lines = text.split("\n")
    need(lines and lines[-1] == "", f"{name}: no final newline")
    head = lines[0].split("\t")
    if columns is not None:
        for c in head:
            need(c in columns, f"{name}: column {c!r}, a column the script does not expect")
        for c in columns:
            need(c in head, f"{name}: no column {c!r}")
        need(head == list(columns), f"{name}: the columns stand in another order")
    out = []
    for i, line in enumerate(lines[1:-1], 2):
        cells = line.split("\t")
        need(len(cells) == len(head), f"{name}:{i}: {len(cells)} cells, header has {len(head)}")
        out.append(dict(zip(head, cells)))
    return out


# ----------------------------------------------------------------- inputs

# The reduced input set holds these files, the fit, fold and dataset tables the second design names, and nothing
# else. The three files of OPTIONAL are hashed when present and read for nothing more.
INPUT_FILES = ("checks.tsv", "SOURCES.sha256", "fits/index.tsv", "release/index.tsv", "release/release.tsv",
               "first/PLAN.md", "first/PLAN.sha256", "first/design.json", "first/eligibility-variants.tsv",
               "first/run.log", "first/study-counts.tsv", "first/times.tsv",
               "second/PLAN.md", "second/PLAN.sha256", "second/design.json", "second/plan-shapes.tsv",
               "second/run.log", "second/study-rule-variation.tsv", "second/times.tsv")
OPTIONAL = ("SOURCES.sha256", "first/study-counts.tsv", "second/plan-shapes.tsv")
COUNT_COLUMNS = [f"{l}_{k}" for l in LEARNERS for k in ("tp", "fn", "tn", "fp", "zero_folds")]
SCORE_COLUMNS = ["genome", "fold", "label"] + list(LEARNERS)
INPUT_COLUMNS = {
    "checks.tsv": ["check", "count"],
    "fits/index.tsv": ["dataset", "plan", "fold_sha256", "input_fingerprint", "revision"] + COUNT_COLUMNS,
    "release/index.tsv": ["dataset", "species", "antibiotic", "input_fingerprint", "revision", "report_sha256",
                          "lineage_status", "lineage_groups", "lineage_folds"] + COUNT_COLUMNS,
    "release/release.tsv": ["key", "value"],
    "first/eligibility-variants.tsv": ["min_studies", "min_study_share", "max_group_share", "datasets", "species",
                                       "cohort_genome_pairs"],
    "first/times.tsv": ["role", "file", "mtime"],
    "second/times.tsv": ["role", "file", "mtime"],
    "second/study-rule-variation.tsv": ["dataset", "species", "genomes", "groups", "draw", "rule",
                                        "changed_share_vs_d0"],
    "fits/": SCORE_COLUMNS,
    "release/lineage/": SCORE_COLUMNS,
    "release/genomes/": ["genome", "label", "study", "lineage", "copy_group", "random_fold"],
}


def read_input(I, rel):
    """One TSV of the input set, by its path in the set, held to the columns the script expects of it."""
    p = Path(I) / rel
    need(p.is_file(), f"inputs: {rel} missing")
    cols = INPUT_COLUMNS.get(rel) or INPUT_COLUMNS[rel.rsplit("/", 1)[0] + "/"]
    return read_tsv(p, cols, rel)


def unexpected_inputs(I, d2, rel):
    """Stop on a file of the input set the script does not expect, and on one it needs and does not find."""
    want = set(INPUT_FILES)
    for e in d2["datasets"]:
        d = e["id"]
        want.add(f"release/genomes/{d}.tsv")
        if rel[d]["lineage_status"] == "scored":
            want.add(f"release/lineage/{d}.tsv")
        for name, p in e["plans"].items():
            want.update((f"fits/{d}.{name}.tsv", f"folds/{p['file']}"))
    have = {p.relative_to(I).as_posix() for p in Path(I).rglob("*") if p.is_file()}
    for x in sorted(have - want):
        what = ("a fit file" if x.startswith("fits/") else "a fold file" if x.startswith("folds/")
                else "a dataset table" if x.startswith(("release/genomes/", "release/lineage/")) else "an input file")
        need(False, f"inputs: {x}: {what} the script does not expect")
    for x in sorted(want - have - set(OPTIONAL)):
        need(False, f"inputs: {x} missing")


def fmt(x):
    if isinstance(x, bool):
        return "yes" if x else "no"
    if isinstance(x, int):
        return str(x)
    if isinstance(x, Fraction):
        return repr(float(x))
    if isinstance(x, float):
        return repr(x)
    if x is None:
        return "-"
    return str(x)


def num(s):
    """A cell written by fmt() back to an exact value."""
    if re.fullmatch(r"-?\d+", s):
        return int(s)
    return Fraction(float(s))


def median(xs):
    xs = sorted(xs)
    need(xs, "median of nothing")
    m = len(xs)
    return Fraction(xs[(m - 1) // 2] + xs[m // 2]) / 2


def mean(xs):
    xs = list(xs)
    need(xs, "mean of nothing")
    return Fraction(sum(Fraction(x) for x in xs), len(xs))


# ----------------------------------------------------------------- columns

COLUMNS = {
    "datasets.tsv": ["dataset", "species", "antibiotic", "genomes", "resistant", "studies", "no_study",
                     "largest_study", "largest_study_share", "lineage_split"],
    "eligibility.tsv": ["dataset", "species", "genomes", "with_study", "studies", "largest_group",
                        "under_studies", "under_share", "over_group", "passes", "entered", "folds"],
    "variants.tsv": ["min_studies", "min_study_share", "max_group_share", "datasets", "species", "pairs"],
    "plans.tsv": ["dataset", "plan", "kind", "draw", "rule", "fold_sha256"]
                 + [f"size_{f}" for f in range(FOLDS)] + [f"resistant_{f}" for f in range(FOLDS)]
                 + ["study_mates", "study_mate_share", "median_study_mates", "lineage_mates",
                    "lineage_mate_share", "median_lineage_mates"],
    "fits.tsv": ["dataset", "plan", "learner", "tp", "fn", "tn", "fp", "recall", "specificity",
                 "balanced_accuracy", "log_loss", "error_rate", "null_fit"],
    "lineage.tsv": ["dataset", "learner", "tp", "fn", "tn", "fp", "recall", "specificity", "balanced_accuracy",
                    "log_loss", "error_rate", "null_fit"],
    "subsets.tsv": ["dataset", "plan", "learner", "subset", "n", "resistant", "tp", "fn", "tn", "fp", "log_loss"],
    "folds.tsv": ["dataset", "fold", "learner", "test", "test_resistant", "train_random_sum", "study_tp",
                  "study_fn", "study_tn", "study_fp", "study_log_loss", "random_tp", "random_fn", "random_tn",
                  "random_fp", "random_log_loss"],
    "distances.tsv": ["dataset", "kind", "draw_a", "draw_b", "changed"],
    "variation.tsv": ["dataset", "draw", "rule", "changed"],
    "order.tsv": ["experiment", "event", "source", "time"],
    "planfiles.tsv": ["experiment", "file", "sha256", "hashed_files", "hash_time", "listed"],
    "settings.tsv": ["key", "value"],
    "checks.tsv": ["check", "count"],
    "totals.tsv": ["key", "value"],
    # derived
    "contrasts.tsv": ["experiment", "learner", "family", "instance", "measure", "median", "low", "high",
                      "reading", "datasets", "worse", "equal", "better"],
    "first.tsv": ["learner", "contrast", "measure", "reading", "keep", "changed_without"],
    "leaveout.tsv": ["left_out", "learner", "family", "pairs", "recall_worse", "loss_worse", "cost", "answer"],
    "rules.tsv": ["learner", "family", "rule", "pairs", "cost", "answer"],
    "thresholds.tsv": ["family", "threshold", "learners"],
    "draws.tsv": ["learner", "family", "axis", "measure", "d0", "d1", "d2", "d3", "d4"],
    "parts.tsv": ["learner", "measure", "datasets", "mean_sr", "mean_sm", "mean_mr", "share", "larger", "exact"],
    "mates.tsv": ["kind", "draw", "study_share", "median_study_mates", "lineage_share", "median_lineage_mates"],
    "shapes.tsv": ["arm", "measure", "median", "low", "high"],
    "nostudy.tsv": ["learner", "family", "measure", "subset", "datasets", "mean", "median", "zero"],
    "foldcut.tsv": ["learner", "group", "folds", "datasets", "genomes", "fold_median_log_loss",
                    "pooled_log_loss", "pooled_recall", "pooled_specificity", "median_shift", "median_ratio"],
    "summary.tsv": ["key", "value"],
}

KEYS = {
    "datasets.tsv": ("dataset",), "eligibility.tsv": ("dataset",),
    "variants.tsv": ("min_studies", "min_study_share", "max_group_share"),
    "plans.tsv": ("dataset", "plan"), "fits.tsv": ("dataset", "plan", "learner"),
    "lineage.tsv": ("dataset", "learner"), "subsets.tsv": ("dataset", "plan", "learner", "subset"),
    "folds.tsv": ("dataset", "fold", "learner"), "distances.tsv": ("dataset", "kind", "draw_a", "draw_b"),
    "variation.tsv": ("dataset", "draw", "rule"),
}

ROW_FILES = ("datasets.tsv", "eligibility.tsv", "variants.tsv", "plans.tsv", "fits.tsv", "lineage.tsv",
             "subsets.tsv", "folds.tsv", "distances.tsv", "variation.tsv", "order.tsv", "planfiles.tsv",
             "settings.tsv", "checks.tsv", "totals.tsv")
DERIVED = ("contrasts.tsv", "first.tsv", "leaveout.tsv", "rules.tsv", "thresholds.tsv", "draws.tsv",
           "parts.tsv", "mates.tsv", "shapes.tsv", "nostudy.tsv", "foldcut.tsv", "summary.tsv")

# What guards each file without the archive; `--check` rebuilds every file from it. The suite pins PINNED by
# sha256, so a change to one of them is also a change to tests/study_split/test_record.py.
PINNED = ("checks.tsv", "distances.tsv", "eligibility.tsv", "folds.tsv", "inputs.sha256", "lineage.tsv", "order.tsv",
          "planfiles.tsv", "plans.tsv", "settings.tsv", "variants.tsv", "variation.tsv")
AGREED = ("datasets.tsv", "fits.tsv", "subsets.tsv", "totals.tsv")
QUOTED = ("checks.tsv", "contrasts.tsv", "datasets.tsv", "draws.tsv", "eligibility.tsv", "first.tsv", "foldcut.tsv",
          "leaveout.tsv", "mates.tsv", "nostudy.tsv", "order.tsv", "parts.tsv", "planfiles.tsv", "plans.tsv",
          "rules.tsv", "settings.tsv", "shapes.tsv", "summary.tsv", "thresholds.tsv")

DESCRIBE = {
    "datasets.tsv": "One row per dataset: species, antibiotic, genomes, resistant genomes, BioProjects, "
                    "genomes with none, the largest BioProject and its share, and whether the release scored "
                    "a lineage split",
    "eligibility.tsv": "The release's datasets against the first plan's three counts, whether each entered, and "
                       "the fold check's reason where it stopped one",
    "variants.tsv": "The first plan's eligibility variants, with the datasets, species and pairs each admits by the three "
                    "counts",
    "plans.tsv": "One row per fold plan: kind, draw, placement rule, fold file sha256, the five fold sizes "
                 "and resistant counts, and test genomes with a study or lineage-cluster mate in training",
    "fits.tsv": "One row per fit and learner: confusion counts, recall, specificity, balanced accuracy, log "
                "loss, error rate, and a null fit flag",
    "lineage.tsv": "The release's lineage split per dataset and learner, in the form of fits.tsv",
    "subsets.tsv": "Per fit and learner, the genomes with and without a BioProject apart, on datasets with at "
                   "least 10 genomes without one",
    "folds.tsv": "The first experiment per study fold and learner: test genomes, training counts, and both "
                 "arms' counts and log loss on the fold's genomes",
    "distances.tsv": "Genomes in another fold between two draws of a kind, after pairing fold labels for the "
                     "largest agreement",
    "variation.tsv": "The plan builder's table: genomes moved from draw 0 under placement rules A and B, per dataset "
                     "and draw",
    "order.tsv": "Times from the hash lists, the run logs, the file times and the second plan's text",
    "planfiles.tsv": "Each plan's file in the input set, its sha256, the files its hash list holds, and its hash time",
    "settings.tsv": "The settings the script applies, the first plan's resample count, and the two design keys "
                    "the fold builders hashed",
    "checks.tsv": "The script's checks with the count that passed each, and the count of the inputs' genome ids "
                  "found in the record",
    "totals.tsv": "Distinct genomes and BioProjects over the datasets",
    "contrasts.tsv": "Every contrast: median over datasets, exact species interval, reading, datasets worse, "
                     "equal and better",
    "first.tsv": "The first experiment's readings, with the leave-one-species-out readings that keep each",
    "leaveout.tsv": "The second experiment with each species left out, and the two largest species together",
    "rules.tsv": "The second experiment's pairs with a cost under the plan's rule and under other rules",
    "thresholds.tsv": "The learners each answer holds for, at every threshold",
    "draws.tsv": "Pairs with a cost, or reading worse, per draw",
    "parts.tsv": "The mean drop as the part matched by fold shape and the part not matched by fold shape",
    "mates.tsv": "Test genomes with a study mate or lineage-cluster mate in training, per kind and draw",
    "shapes.tsv": "Fold sizes and class balance per arm",
    "nostudy.tsv": "The three contrasts on genomes without a BioProject, and on genomes with one in all the datasets "
                   "and in the same datasets",
    "foldcut.tsv": "The first experiment per study fold, in groups by training size and class shift",
    "summary.tsv": "Single figures the page quotes",
}


def tsv(name, rows):
    cols = COLUMNS[name]
    for r in rows:
        need(len(r) == len(cols), f"{name}: row of {len(r)} cells")
        for c in r:
            need("\t" not in c and "\n" not in c and c != "", f"{name}: bad cell {c!r}")
    if name in KEYS:
        idx = [cols.index(k) for k in KEYS[name]]
        keyf = lambda r: tuple((0, int(r[i]), "") if re.fullmatch(r"\d+", r[i]) else (1, 0, r[i]) for i in idx)
        rows = sorted(rows, key=keyf)
        keys = [keyf(r) for r in rows]
        need(len(set(keys)) == len(keys), f"{name}: duplicate key")
    return ("\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows)).encode()


def parse(text, name):
    lines = text.split("\n")
    need(lines[-1] == "", f"{name}: no final newline")
    need(lines[0].split("\t") == COLUMNS[name], f"{name}: header differs")
    return [dict(zip(COLUMNS[name], line.split("\t"))) for line in lines[1:-1]]


# ----------------------------------------------------------------- exact interval

FACT = [math.factorial(i) for i in range(16)]


def interval(groups):
    """groups: one list of exact values per species. Returns (median, low, high).

    The median over a resample is the median of the datasets' values with each
    species' values repeated as often as the species is drawn. Over every
    multiset of k species with replacement, weighted by its multinomial count
    out of k**k: low is the first sorted median whose cumulative weight exceeds
    1/40, high the first that reaches 39/40.

    For speed, the cumulative count of each species at each sorted position is
    packed into one integer per species; for a multiset the sum of those
    integers holds the cumulative count of the resample at every position, and
    the rank of an order statistic is read with one addition and a bit count.
    """
    k = len(groups)
    need(k >= 1 and all(groups), "interval of an empty group")
    vals = sorted((v, s) for s, g in enumerate(groups) for v in g)
    m = len(vals)
    n = [len(g) for g in groups]
    top = k * max(n)
    w = top.bit_length() + 1
    half = 1 << (w - 1)
    rep = ((1 << (w * m)) - 1) // ((1 << w) - 1)
    flags = rep << (w - 1)
    cum = [0] * k
    packed = [0] * k
    for i, (_, s) in enumerate(vals):
        cum[s] += 1
        for t in range(k):
            if cum[t]:
                packed[t] |= cum[t] << (w * i)
    add = [(half - 1 - r) * rep for r in range(top + 1)]
    fk = FACT[k]
    acc = {}

    def rec(s, left, x, total, den):
        if s == k - 1:
            x += left * packed[s]
            total += left * n[s]
            den *= FACT[left]
            r1 = (total - 1) >> 1
            a = m - ((x + add[r1]) & flags).bit_count()
            b = a if total & 1 else m - ((x + add[r1 + 1]) & flags).bit_count()
            key = a * m + b
            acc[key] = acc.get(key, 0) + fk // den
            return
        p, ns = packed[s], n[s]
        for c in range(left + 1):
            rec(s + 1, left - c, x, total, den * FACT[c])
            x += p
            total += ns

    rec(0, k, 0, 0, 1)
    dist = sorted((Fraction(vals[key // m][0] + vals[key % m][0]) / 2, wt) for key, wt in acc.items())
    tot = k ** k
    run = 0
    lo = hi = None
    for med, wt in dist:
        run += wt
        if lo is None and 40 * run > tot:
            lo = med
        if 40 * run >= 39 * tot:
            hi = med
            break
    return median([v for v, _ in vals]), lo, hi


def reading(lo, hi, m):
    if lo > 0:
        return "worse" if WORSE_UP[m] else "better"
    if hi < 0:
        return "better" if WORSE_UP[m] else "worse"
    return "not shown"


def outcome(rr, ll, both=False):
    s = {rr, ll}
    if "better" in s:
        return "mixed" if "worse" in s else "benefit"
    if both:
        return "cost" if rr == ll == "worse" else "not shown"
    return "cost" if "worse" in s else "not shown"


def experiment_rule(per):
    if any(v in ("benefit", "mixed") for v in per):
        return "mixed"
    return "cost shown" if sum(v == "cost" for v in per) >= EXPERIMENT_LEARNERS else "cost not shown"


# ----------------------------------------------------------------- inputs

def measure_fit(y, scores):
    tp = fn = tn = fp = 0
    losses = []
    for t, s in zip(y, scores):
        call = s >= CUTOFF
        # Clip the score given to the true class, so that a wrong 0 and a wrong 1 lose the same amount.
        q = min(max(s if t else 1 - s, CLIP), 1 - CLIP)
        losses.append(-math.log(q))
        if t:
            tp += call
            fn += not call
        else:
            fp += call
            tn += not call
    nn = len(y)
    return tp, fn, tn, fp, (math.fsum(losses) / nn if nn else None)


def fit_cells(tp, fn, tn, fp, ll):
    rec = Fraction(tp, tp + fn)
    spe = Fraction(tn, tn + fp)
    return [fmt(tp), fmt(fn), fmt(tn), fmt(fp), fmt(rec), fmt(spe), fmt((rec + spe) / 2), fmt(ll),
            fmt(Fraction(fn + fp, tp + fn + tn + fp))]


def mates(fold, label):
    size = Counter(x for x in label if x is not None)
    infold = Counter((label[i], fold[i]) for i in range(len(fold)) if label[i] is not None)
    return [size[label[i]] - infold[(label[i], fold[i])] if label[i] is not None else 0 for i in range(len(fold))]


def best_agreement(a, b):
    import itertools
    t = Counter(zip(a, b))
    return max(sum(t[(f, p[f])] for f in range(FOLDS)) for p in itertools.permutations(range(FOLDS)))


def rule_of(plan, study_rule):
    kind, draw = plan.split("-d")
    if kind == "random":
        return "-"
    r = study_rule[draw]
    return "A" if r == "lane" else r


def one_dataset(task):
    inputs, e, rel, fidx, study_rule = task
    I = Path(inputs)
    d = e["id"]
    g = read_input(I, f"release/genomes/{d}.tsv")
    ids = [r["genome"] for r in g]
    y = [r["label"] == "1" for r in g]
    study = [r["study"] or None for r in g]
    lin = [r["lineage"] for r in g]
    copies = defaultdict(list)
    for i, r in enumerate(g):
        if r["copy_group"]:
            copies[r["copy_group"]].append(i)
    n = len(ids)
    res = sum(y)
    need(e["genomes"] == n and e["resistant"] == res, f"{d}: genome or resistant count differs from design")
    need(rel["revision"] == "14", f"{d}: release revision")
    studies = Counter(s for s in study if s)
    no_study = sum(1 for s in study if not s)
    big = max(studies.values())
    scored = rel["lineage_status"] == "scored"
    out = {"datasets": [[d, e["species_slug"], rel["antibiotic"], fmt(n), fmt(res), fmt(len(studies)),
                         fmt(no_study), fmt(big), fmt(Fraction(big, n)), "scored" if scored else "unavailable"]],
           "plans": [], "fits": [], "lineage": [], "subsets": [], "folds": [], "distances": [],
           "checks": Counter(), "genomes": sorted(ids), "studies": sorted(studies)}
    ck = out["checks"]
    folds, scores = {}, {}
    for name in sorted(e["plans"]):
        p = e["plans"][name]
        fx = fidx[(d, name)]
        rows = read_input(I, f"fits/{d}.{name}.tsv")
        need([r["genome"] for r in rows] == ids, f"{d} {name}: genomes or order differ from the release")
        need([r["label"] == "1" for r in rows] == y, f"{d} {name}: labels differ from the release")
        fold = [int(r["fold"]) for r in rows]
        raw = (I / "folds" / p["file"]).read_bytes()
        need(sha256_bytes(raw) == p["sha256"] == fx["fold_sha256"], f"{d} {name}: fold file sha256")
        need(fx["input_fingerprint"] == rel["input_fingerprint"] and fx["revision"] == "14",
             f"{d} {name}: fingerprint or revision")
        doc = json.loads(raw)
        need(doc["train"] == [] and doc["unused"] == [] and len(doc["test"]) == FOLDS, f"{d} {name}: fold file form")
        listed = Counter(x for f in doc["test"] for x in f)
        need(sorted(listed) == sorted(ids) and set(listed.values()) == {1}, f"{d} {name}: genomes listed once")
        pf = {x: f for f, members in enumerate(doc["test"]) for x in members}
        need([pf[x] for x in ids] == fold, f"{d} {name}: fit folds differ from the fold file")
        sizes = [fold.count(f) for f in range(FOLDS)]
        rs = [sum(1 for i in range(n) if fold[i] == f and y[i]) for f in range(FOLDS)]
        need(sizes == p["sizes"] and rs == p["resistant"], f"{d} {name}: fold counts differ from design")
        for f in range(FOLDS):
            need(0 < rs[f] < sizes[f] and 0 < res - rs[f] < n - sizes[f], f"{d} {name}: fold {f} lacks a class")
        need(all(len({fold[i] for i in m}) == 1 for m in copies.values()), f"{d} {name}: copy group split")
        if p["kind"] == "study":
            sf = defaultdict(set)
            for i in range(n):
                if study[i]:
                    sf[study[i]].add(fold[i])
            need(all(len(v) == 1 for v in sf.values()), f"{d} {name}: BioProject split")
            ck["study_plans_with_no_study_split"] += 1
        ck["fold_files_checked"] += 1
        ck["copy_groups_whole_plans"] += 1
        sm, lm = mates(fold, study), mates(fold, lin)
        ws = [i for i in range(n) if study[i]]
        smc = sum(1 for i in ws if sm[i] > 0)
        lmc = sum(1 for x in lm if x > 0)
        out["plans"].append([d, name, p["kind"], fmt(p["draw"]), rule_of(name, study_rule), p["sha256"]]
                            + [fmt(x) for x in sizes] + [fmt(x) for x in rs]
                            + [fmt(smc), fmt(Fraction(smc, n)), fmt(median([sm[i] for i in ws])),
                               fmt(lmc), fmt(Fraction(lmc, n)), fmt(median(lm))])
        folds[name] = fold
        sc = {}
        for l in LEARNERS:
            s = [float(r[l]) for r in rows]
            sc[l] = s
            if l == "scm":
                ck["scm_fits_scores_0_or_1"] += all(v in (0.0, 1.0) for v in s)
            tp, fn, tn, fp, ll = measure_fit(y, s)
            need([tp, fn, tn, fp] == [int(fx[f"{l}_{k}"]) for k in ("tp", "fn", "tn", "fp")],
                 f"{d} {name} {l}: counts differ from the evaluator's")
            ck["learner_fits_equal_evaluator_counts"] += 1
            null = l in ("gbt", "lasso") and int(fx[f"{l}_zero_folds"]) > 0
            out["fits"].append([d, name, l] + fit_cells(tp, fn, tn, fp, ll) + [fmt(null)])
            if no_study >= NO_STUDY_MIN:
                for subset, keep in (("study", True), ("no-study", False)):
                    sel = [i for i in range(n) if bool(study[i]) == keep]
                    tp, fn, tn, fp, ll = measure_fit([y[i] for i in sel], [s[i] for i in sel])
                    out["subsets"].append([d, name, l, subset, fmt(len(sel)), fmt(tp + fn), fmt(tp), fmt(fn),
                                           fmt(tn), fmt(fp), fmt(ll)])
        scores[name] = sc
    need([int(r["random_fold"]) for r in g] == folds["random-d0"], f"{d}: random draw 0 differs from the release")
    ck["random_d0_equal_release_folds"] += 1
    for i in DRAWS:
        a, b = out["plans"][[r[1] for r in out["plans"]].index(f"study-d{i}")], \
            out["plans"][[r[1] for r in out["plans"]].index(f"matched-d{i}")]
        need(a[6:16] == b[6:16], f"{d}: matched draw {i} counts differ from study draw {i}")
        ck["matched_plans_equal_study_counts"] += 1
    for kind in KINDS:
        for a in DRAWS:
            for b in DRAWS:
                if a < b:
                    ch = n - best_agreement(folds[f"{kind}-d{a}"], folds[f"{kind}-d{b}"])
                    out["distances"].append([d, kind, fmt(a), fmt(b), fmt(ch)])
    if scored:
        rows = read_input(I, f"release/lineage/{d}.tsv")
        need([r["genome"] for r in rows] == ids and [r["label"] == "1" for r in rows] == y,
             f"{d}: lineage genomes or labels")
        lf = [int(r["fold"]) for r in rows]
        groups = defaultdict(set)
        for i in range(n):
            groups[lin[i]].add(lf[i])
        need(all(len(v) == 1 for v in groups.values()), f"{d}: lineage cluster split in the lineage split")
        need(len(groups) == int(rel["lineage_groups"]), f"{d}: lineage cluster count")
        ck["lineage_splits_clusters_whole"] += 1
        for l in LEARNERS:
            tp, fn, tn, fp, ll = measure_fit(y, [float(r[l]) for r in rows])
            need([tp, fn, tn, fp] == [int(rel[f"{l}_{k}"]) for k in ("tp", "fn", "tn", "fp")],
                 f"{d} lineage {l}: counts differ from the evaluator's")
            ck["learner_fits_equal_evaluator_counts"] += 1
            null = l in ("gbt", "lasso") and int(rel[f"{l}_zero_folds"]) > 0
            out["lineage"].append([d, l] + fit_cells(tp, fn, tn, fp, ll) + [fmt(null)])
    swf, rrf = folds["study-d0"], folds["random-d0"]
    rr_size = Counter(rrf)
    for f in range(FOLDS):
        members = [i for i in range(n) if swf[i] == f]
        t_res = sum(y[i] for i in members)
        trs = sum(n - rr_size[rrf[i]] for i in members)
        for l in LEARNERS:
            tt = [y[i] for i in members]
            a = measure_fit(tt, [scores["study-d0"][l][i] for i in members])
            b = measure_fit(tt, [scores["random-d0"][l][i] for i in members])
            out["folds"].append([d, fmt(f), l, fmt(len(members)), fmt(t_res), fmt(trs)]
                                + [fmt(x) for x in a] + [fmt(x) for x in b])
    return out


def times_rows(I, first_fits, second_fits):
    """Times from each experiment's file times, the first line of its hash list, its run log and the second plan's text."""
    rows = []
    hashed = {}
    for exp in ("first", "second"):
        for r in read_input(I, f"{exp}/times.tsv"):
            rows.append([exp, r["role"], f"file time in the input set's {exp}/times.tsv", r["mtime"]])
        head = (I / exp / "PLAN.sha256").read_text().split("\n")[0]
        mt = re.search(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", head)
        need(head.startswith("#") and mt, f"{exp}: no time in the hash list's first line")
        rows.append([exp, "hash_time", f"first line of the input set's {exp}/PLAN.sha256", mt.group(1)])
        hashed[exp] = mt.group(1)
    log1 = (I / "first" / "run.log").read_text().split("\n")
    ends = [m.group(1) for m in (re.match(r"(\S+Z) \S+\.(release-random|study-whole) done", x) for x in log1) if m]
    need(len(ends) == first_fits, f"first run log: {len(ends)} fit lines, {first_fits} fold files")
    rows.append(["first", "first_fit_end", "the input set's first/run.log", ends[0]])
    rows.append(["first", "last_fit_end", "the input set's first/run.log", ends[-1]])
    log2 = (I / "second" / "run.log").read_text().split("\n")
    fits = [m.groups() for m in (re.match(r"(\S+Z) fit \S+ done start (\S+Z);", x) for x in log2) if m]
    need(len(fits) == second_fits, f"second run log: {len(fits)} fit lines, {second_fits} fold plans")
    start = [m.group(1) for m in (re.match(r"(\S+Z) fits start", x) for x in log2) if m]
    need(len(start) == 1, "second run log: fits start")
    rows.append(["second", "fits_start", "the input set's second/run.log", start[0]])
    rows.append(["second", "first_fit_end", "the input set's second/run.log", min(f[0] for f in fits)])
    rows.append(["second", "last_fit_end", "the input set's second/run.log", max(f[0] for f in fits)])
    raw = (I / "second" / "PLAN.md").read_bytes()
    plan = raw.decode()
    mt = re.findall(r"I fixed at (\d\d:\d\d:\d\dZ)", plan)
    need(len(mt) == 1, "second plan: the stated time of the rule choice")
    # The plan gives a time of day only; its date is the hash list's.
    rows.append(["second", "stated_rule_choice", f"text of the second plan, sha256 {sha256_bytes(raw)}",
                 hashed["second"][:11] + mt[0]])
    return rows


def planfile_rows(I, files):
    rows, checked = [], 0
    for exp in ("first", "second"):
        lines = [x for x in (I / exp / "PLAN.sha256").read_text().split("\n") if x and not x.startswith("#")]
        listed = {}
        for x in lines:
            h, name = x.split(None, 1)
            listed[name.strip()] = h
        plan = (I / exp / "PLAN.md").read_bytes()
        need(listed.get("PLAN.md") == sha256_bytes(plan), f"{exp}: PLAN.md differs from its hash list")
        head = (I / exp / "PLAN.sha256").read_text().split("\n")[0]
        t = re.search(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", head).group(1)
        for name, h in listed.items():
            base = Path(name).name
            local = {"results/design.json": I / exp / "design.json",
                     "results/plan-shapes.tsv": I / exp / "plan-shapes.tsv",
                     "results/study-rule-variation.tsv": I / exp / "study-rule-variation.tsv",
                     "results/counts/study-counts.tsv": I / exp / "study-counts.tsv",
                     "results/counts/eligibility-variants.tsv": I / exp / "eligibility-variants.tsv",
                     "PLAN.md": I / exp / "PLAN.md"}.get(name)
            if local is None and base in files:
                local = I / "folds" / files[base]
            if local is not None and local.is_file():
                need(sha256_bytes(local.read_bytes()) == h, f"{exp}: {name} differs from its hash list")
                checked += 1
        rows.append([exp, f"{exp}/PLAN.md", sha256_bytes(plan), fmt(len(lines)), t, "yes"])
    return rows, checked


def eligibility_rows(d1):
    rule = d1["rule"]
    lo_share = Fraction(str(rule["min_study_share"]))
    hi_group = Fraction(str(rule["max_group_share"]))
    rows = []
    for e in d1["datasets"]:
        c = e["counts"]
        gen, ws, st, lg = (int(c[k]) for k in ("genomes", "with_study", "studies", "largest_group"))
        us = st < rule["min_studies"]
        ush = Fraction(ws, gen) < lo_share
        og = Fraction(lg, gen) > hi_group
        passes = not (us or ush or og)
        rs = " | ".join(e["reasons"])
        need(us == ("studies, under" in rs) and ush == ("with a study, under" in rs)
             and og == ("genomes, over" in rs), f"{e['id']}: reasons differ from the counts")
        need(e["eligible"] == (passes and "study folds:" not in rs), f"{e['id']}: eligibility")
        folds = [x.split("study folds:", 1)[1].strip() for x in e["reasons"] if "study folds:" in x]
        rows.append([e["id"], e["species_slug"], fmt(gen), fmt(ws), fmt(st), fmt(lg), fmt(us), fmt(ush), fmt(og),
                     fmt(passes), fmt(e["eligible"]), "; ".join(folds) or "-"])
    return rows


def variant_rows(I, elig):
    rows = []
    for v in read_input(I, "first/eligibility-variants.tsv"):
        ms, sh, gr = int(v["min_studies"]), Fraction(v["min_study_share"]), Fraction(v["max_group_share"])
        keep = [r for r in elig if int(r[4]) >= ms and Fraction(int(r[3]), int(r[2])) >= sh
                and Fraction(int(r[5]), int(r[2])) <= gr]
        got = [len(keep), len({r[1] for r in keep}), sum(int(r[2]) for r in keep)]
        need(got == [int(v["datasets"]), int(v["species"]), int(v["cohort_genome_pairs"])],
             f"variant {v}: {got}")
        rows.append([v["min_studies"], v["min_study_share"], v["max_group_share"]] + [fmt(x) for x in got])
    return rows


def variation_rows(I):
    rows = []
    for r in read_input(I, "second/study-rule-variation.tsv"):
        share = Fraction(r["changed_share_vs_d0"])
        ch = round(share * int(r["genomes"]))
        need(abs(Fraction(ch, int(r["genomes"])) - share) < Fraction(1, 10**6), f"variation {r}")
        rows.append([r["dataset"], r["draw"], r["rule"], fmt(ch)])
    return rows


def build_rows(inputs, jobs):
    I = Path(inputs)
    for p in ("checks.tsv", "first/design.json", "second/design.json", "fits/index.tsv", "release/index.tsv"):
        need((I / p).is_file(), f"inputs: {p} missing")
    d1 = json.loads((I / "first" / "design.json").read_text())
    d2 = json.loads((I / "second" / "design.json").read_text())
    ids = [e["id"] for e in d2["datasets"]]
    need(sorted(ids) == sorted(e["id"] for e in d1["datasets"] if e["eligible"]),
         "the two designs hold different datasets")
    need(d1["meta_sha256_gz"] == d2["meta_sha256_gz"], "the designs read different genome tables")
    rel_rows = read_input(I, "release/index.tsv")
    rel = {r["dataset"]: r for r in rel_rows}
    for d in sorted(set(rel) - set(ids)):
        need(False, f"release/index.tsv: {d}: a dataset the script does not expect")
    need(len(rel_rows) == len(rel) == len(ids), "release/index.tsv: a dataset missing or listed twice")
    fit_rows = read_input(I, "fits/index.tsv")
    fidx = {(r["dataset"], r["plan"]): r for r in fit_rows}
    planned = {(e["id"], k) for e in d2["datasets"] for k in e["plans"]}
    for d, k in sorted(set(fidx) - planned):
        need(False, f"fits/index.tsv: {d} {k}: a fit the script does not expect")
    need(len(fit_rows) == len(fidx) == len(planned) == 15 * len(ids), "fits index: plan count")
    unexpected_inputs(I, d2, rel)
    # The line of the second plan's rule: the script reads it from its own constant, so the two must agree.
    said = re.findall(r"a median under\s+(\S+)\s+over the", (I / "second" / "PLAN.md").read_text())
    need(said == [RULE_LINE], f"second plan: the rule's line reads {said}, the script's RULE_LINE {RULE_LINE}")
    first_sha = {(e["id"], k): v["sha256"] for e in d1["datasets"] if e["eligible"]
                 for k, v in e["fold_files"].items()}
    for e in d2["datasets"]:
        need(first_sha[(e["id"], "study-whole")] == e["plans"]["study-d0"]["sha256"]
             and first_sha[(e["id"], "release-random")] == e["plans"]["random-d0"]["sha256"],
             f"{e['id']}: draw 0 is not the first experiment's fold files")
    tasks = [(str(I), e, rel[e["id"]], {k: v for k, v in fidx.items() if k[0] == e["id"]}, d2["study_rule"])
             for e in sorted(d2["datasets"], key=lambda e: e["id"])]
    if jobs > 1:
        with Pool(jobs) as pool:
            parts = pool.map(one_dataset, tasks)
    else:
        parts = [one_dataset(t) for t in tasks]
    R = {k: [] for k in ("datasets", "plans", "fits", "lineage", "subsets", "folds", "distances")}
    ck = Counter()
    genomes, studies = set(), set()
    for p in parts:
        for k in R:
            R[k] += p[k]
        ck.update(p["checks"])
        genomes.update(p["genomes"])
        studies.update(p["studies"])
    files = {e["plans"][k]["file"]: e["plans"][k]["file"] for e in d2["datasets"] for k in e["plans"]}
    second_of = {e["id"]: e for e in d2["datasets"]}
    first_names = {Path(v["file"]).name:
                   second_of[e["id"]]["plans"]["study-d0" if k == "study-whole" else "random-d0"]["file"]
                   for e in d1["datasets"] if e["eligible"] for k, v in e["fold_files"].items()}
    files.update(first_names)
    pf_rows, hashed_ok = planfile_rows(I, files)
    ck["hash_list_entries_equal_inputs"] = hashed_ok
    for r in read_input(I, "checks.tsv"):
        ck["inputs_" + r["check"]] = int(r["count"])
    elig = eligibility_rows(d1)
    release = {r["key"]: r["value"] for r in read_input(I, "release/release.tsv")}
    # The first plan's own interval: the count of species resamples it fixed, which the exact interval replaces.
    said = re.findall(r"(\d{1,3}(?:,\d{3})*) resamples\s+of\s+species", (I / "first" / "PLAN.md").read_text())
    need(len(said) == 1, f"first plan: {len(said)} counts of species resamples")
    resamples = int(said[0].replace(",", ""))
    need(release["evaluation_revision"] == "14", "release revision")
    settings = [
        ("release_revision", release["evaluation_revision"]), ("release_generated_at", release["generated_at"]),
        ("folds", fmt(d2["folds"])), ("draws", fmt(d2["draws"])), ("cutoff", fmt(CUTOFF)), ("clip", fmt(CLIP)),
        ("first_design_key", d1["key"]), ("second_design_key", d2["key"]),
        ("min_studies", fmt(d1["rule"]["min_studies"])), ("min_study_share", fmt(d1["rule"]["min_study_share"])),
        ("max_group_share", fmt(d1["rule"]["max_group_share"])),
        ("study_rule", ",".join(f"d{i}:" + rule_of(f"study-d{i}", d2["study_rule"]) for i in DRAWS)),
        ("study_rule_line", RULE_LINE), ("threshold_SR", fmt(THRESHOLD["SR"])),
        ("threshold_SM", fmt(THRESHOLD["SM"])), ("threshold_MR", fmt(THRESHOLD["MR"])),
        ("experiment_rule_learners", fmt(EXPERIMENT_LEARNERS)), ("no_study_min", fmt(NO_STUDY_MIN)),
        ("first_loss", "log_loss"), ("second_loss", "error_rate for scm, log_loss for the others"),
        ("genome_table_sha256_gz", d2["meta_sha256_gz"]),
        ("interval", "exact: every multiset of the species with replacement, low above 1/40, high at 39/40"),
        ("first_plan_resamples", fmt(resamples)),
    ]
    tables = {
        "datasets.tsv": R["datasets"], "eligibility.tsv": elig, "variants.tsv": variant_rows(I, elig),
        "plans.tsv": R["plans"], "fits.tsv": R["fits"], "lineage.tsv": R["lineage"], "subsets.tsv": R["subsets"],
        "folds.tsv": R["folds"], "distances.tsv": R["distances"], "variation.tsv": variation_rows(I),
        "order.tsv": times_rows(I, sum(len(e["fold_files"]) for e in d1["datasets"] if e["eligible"]),
                                sum(len(e["plans"]) for e in d2["datasets"])), "planfiles.tsv": pf_rows,
        "settings.tsv": [[k, v] for k, v in settings],
        "checks.tsv": [[k, fmt(v)] for k, v in sorted(ck.items())],
        "totals.tsv": [["genomes", fmt(len(genomes))], ["studies", fmt(len(studies))]],
    }
    return tables, genomes


# ----------------------------------------------------------------- derived

# scm scores only 0 or 1, so its clipped log loss is a right call's loss plus PER_ERROR times its error rate.
RIGHT_LOSS = Fraction(-math.log(1 - CLIP))
PER_ERROR = Fraction(-math.log(CLIP) - (-math.log(1 - CLIP)))


def loss_value(learner, log_loss, tp, fn, tn, fp):
    """A log loss as derive() reads it: scm's in exact form from its error rate, so that its ties are the error
    rate's, and the others' as written."""
    if learner == "scm":
        return RIGHT_LOSS + PER_ERROR * Fraction(fn + fp, tp + fn + tn + fp)
    return Fraction(float(log_loss))


def arm_loss(row, arm):
    """loss_value of one arm of a folds.tsv row."""
    counts = (int(row[f"{arm}_{k}"]) for k in ("tp", "fn", "tn", "fp"))
    return loss_value(row["learner"], row[f"{arm}_log_loss"], *counts)


def values(rows):
    V = {}
    for r in rows:
        tp, fn, tn, fp = (int(r[k]) for k in ("tp", "fn", "tn", "fp"))
        rec, spe = Fraction(tp, tp + fn), Fraction(tn, tn + fp)
        V[(r["dataset"], r.get("plan", "lineage"), r["learner"])] = {
            "recall": rec, "specificity": spe, "balanced_accuracy": (rec + spe) / 2,
            "log_loss": loss_value(r["learner"], r["log_loss"], tp, fn, tn, fp),
            "error_rate": Fraction(fn + fp, tp + fn + tn + fp),
            "null": r["null_fit"] == "yes"}
    return V


def loss_of(exp, l):
    return "error_rate" if exp == "second" and l == "scm" else "log_loss"


def measures_of(l):
    return MEASURES if l == "scm" else MEASURES[:4]


def instances(fam):
    if fam == "SM":
        return [str(i) for i in DRAWS]
    return [f"{i},{j}" for i in DRAWS for j in DRAWS]


def run_intervals(jobs_list, jobs):
    groups = [g for _, g in jobs_list]
    if jobs > 1:
        with Pool(jobs) as pool:
            res = pool.map(interval, groups, chunksize=16)
    else:
        res = [interval(g) for g in groups]
    return {k: r for (k, _), r in zip(jobs_list, res)}


def subset_value(r, m):
    """One measure of a subsets.tsv row, or None where the subset lacks the class the measure needs."""
    tp, fn, tn, fp = (int(r[k]) for k in ("tp", "fn", "tn", "fp"))
    if m == "recall":
        return Fraction(tp, tp + fn) if tp + fn else None
    if m == "specificity":
        return Fraction(tn, tn + fp) if tn + fp else None
    if m == "error_rate":
        return Fraction(fn + fp, tp + fn + tn + fp)
    return loss_value(r["learner"], r["log_loss"], tp, fn, tn, fp)


def nostudy_rows(sub, ns):
    """The three contrasts over the mean of draws, per dataset of ns: on genomes without a BioProject (no-study), on
    genomes with one in every dataset (study), and on genomes with one in the datasets the no-study row counts
    (study-same), so that the two sides of a comparison average the same datasets."""
    nos = []
    for l in LEARNERS:
        for fam, (a, b) in SECOND.items():
            for mm in ("recall", "specificity", "loss"):
                m = mm if mm != "loss" else loss_of("second", l)
                val = {}
                for subset in ("no-study", "study"):
                    for d in ns:
                        va = [subset_value(sub[(d, f"{a}-d{i}", l, subset)], m) for i in DRAWS]
                        vb = [subset_value(sub[(d, f"{b}-d{i}", l, subset)], m) for i in DRAWS]
                        val[(subset, d)] = None if None in va or None in vb else (sum(va) - sum(vb)) / 5
                same = [d for d in ns if val[("no-study", d)] is not None]
                for subset, among, of in (("no-study", ns, "no-study"), ("study", ns, "study"),
                                          ("study-same", same, "study")):
                    xs = [val[(of, d)] for d in among if val[(of, d)] is not None]
                    nos.append([l, fam, mm, subset, fmt(len(xs)), fmt(mean(xs)), fmt(median(xs)),
                                fmt(sum(1 for x in xs if x == 0))])
    return nos


def derive(t, jobs=1):
    """The derived files from the row files. t maps a file name to its rows."""
    ds = {r["dataset"]: r for r in t["datasets.tsv"]}
    datasets = sorted(ds)
    sp = {d: ds[d]["species"] for d in datasets}
    species = sorted(set(sp.values()))
    per_sp = Counter(sp.values())
    two = sorted(sorted(species, key=lambda s: (-per_sp[s], s))[:2])
    V = values(t["fits.tsv"])
    V.update(values(t["lineage.tsv"]))
    lin_ds = sorted({r["dataset"] for r in t["lineage.tsv"]})

    def diff2(d, l, fam, inst, m):
        a, b = SECOND[fam]
        if inst == "mean":
            return (sum(V[(d, f"{a}-d{i}", l)][m] for i in DRAWS) - sum(V[(d, f"{b}-d{i}", l)][m] for i in DRAWS)) / 5
        i, j = (int(inst), int(inst)) if fam == "SM" else map(int, inst.split(","))
        return V[(d, f"{a}-d{i}", l)][m] - V[(d, f"{b}-d{j}", l)][m]

    def diff1(d, l, fam, m):
        a, b = FIRST[fam]
        return V[(d, a, l)][m] - V[(d, b, l)][m]

    def grouped(fn, ds_list, drop=()):
        g = defaultdict(list)
        for d in ds_list:
            if sp[d] not in drop:
                g[sp[d]].append(fn(d))
        return [g[s] for s in sorted(g)]

    J = []
    for l in LEARNERS:
        for fam in SECOND:
            for inst in instances(fam) + ["mean"]:
                for m in measures_of(l):
                    J.append((("second", l, fam, inst, m, "all"),
                              grouped(lambda d: diff2(d, l, fam, inst, m), datasets)))
                if inst == "mean":
                    continue
                for drop in [(s,) for s in species] + [tuple(two)]:
                    scope = drop[0] if len(drop) == 1 else "two-large"
                    for m in ("recall", loss_of("second", l)):
                        J.append((("second", l, fam, inst, m, scope),
                                  grouped(lambda d: diff2(d, l, fam, inst, m), datasets, drop)))
        for fam in FIRST:
            dl = datasets if fam == "study" else lin_ds
            for m in measures_of(l):
                for drop in [()] + [(s,) for s in species]:
                    J.append((("first", l, fam, "-", m, drop[0] if drop else "all"),
                              grouped(lambda d: diff1(d, l, fam, m), dl, drop)))
    res = run_intervals(J, jobs)
    R = {}
    for key, groups in J:
        pt, lo, hi = res[key]
        m = key[4]
        vals = [v for g in groups for v in g]
        up = 1 if WORSE_UP[m] else -1
        R[key] = {"median": pt, "low": lo, "high": hi, "reading": reading(lo, hi, m), "datasets": len(vals),
                  "worse": sum(1 for v in vals if up * v > 0), "equal": sum(1 for v in vals if v == 0),
                  "better": sum(1 for v in vals if up * v < 0)}
    out = {}
    rows = []
    for key, r in R.items():
        if key[5] == "all":
            rows.append([key[0], key[1], key[2], key[3], key[4], fmt(r["median"]), fmt(r["low"]), fmt(r["high"]),
                         r["reading"], fmt(r["datasets"]), fmt(r["worse"]), fmt(r["equal"]), fmt(r["better"])])
    out["contrasts.tsv"] = rows

    def rd(l, fam, inst, m, scope="all"):
        return R[("second", l, fam, inst, m, scope)]["reading"]

    def oc(l, fam, inst, scope="all", rule="either"):
        lm = loss_of("second", l)
        if rule == "either":
            return outcome(rd(l, fam, inst, "recall", scope), rd(l, fam, inst, lm, scope))
        if rule == "both":
            return outcome(rd(l, fam, inst, "recall", scope), rd(l, fam, inst, lm, scope), both=True)
        if rule == "recall":
            return outcome(rd(l, fam, inst, "recall", scope), "not shown")
        if rule == "loss":
            return outcome("not shown", rd(l, fam, inst, lm, scope))
        return outcome(rd(l, fam, inst, "balanced_accuracy", scope), rd(l, fam, inst, lm, scope))

    def cost(l, fam, scope="all", rule="either"):
        return sum(oc(l, fam, i, scope, rule) == "cost" for i in instances(fam))

    def answer(fam, c):
        return "shown" if c >= THRESHOLD[fam] else "not shown"

    out["rules.tsv"] = [[l, fam, rule, fmt(len(instances(fam))), fmt(cost(l, fam, "all", rule)),
                         answer(fam, cost(l, fam, "all", rule))] for l in LEARNERS for fam in SECOND for rule in RULES]
    thr = []
    for fam in SECOND:
        top = len(instances(fam))
        for th in range(top // 2 + 1, top + 1):
            ls = [l for l in LEARNERS if cost(l, fam) >= th]
            thr.append([fam, fmt(th), ",".join(ls) if ls else "none"])
    out["thresholds.tsv"] = thr
    lo_rows = []
    for scope in species + ["two-large"]:
        for l in LEARNERS:
            lm = loss_of("second", l)
            for fam in SECOND:
                ins = instances(fam)
                c = cost(l, fam, scope)
                lo_rows.append([scope, l, fam, fmt(len(ins)),
                                fmt(sum(rd(l, fam, i, "recall", scope) == "worse" for i in ins)),
                                fmt(sum(rd(l, fam, i, lm, scope) == "worse" for i in ins)), fmt(c), answer(fam, c)])
        ex = sum(experiment_rule([oc(l, "SR", i, scope) for l in LEARNERS]) == "cost shown" for i in instances("SR"))
        lo_rows.append([scope, "all", "SR", fmt(len(instances("SR"))), "-", "-", fmt(ex), "-"])
    out["leaveout.tsv"] = lo_rows
    dr = []
    for l in LEARNERS:
        lm = loss_of("second", l)
        for fam in SECOND:
            axes = (("draw", 0),) if fam == "SM" else ((SECOND[fam][0], 0), (SECOND[fam][1], 1))
            for axis, pos in axes:
                for mname, test in (("cost", lambda i: oc(l, fam, i) == "cost"),
                                    ("recall", lambda i: rd(l, fam, i, "recall") == "worse"),
                                    ("specificity", lambda i: rd(l, fam, i, "specificity") == "worse"),
                                    ("loss", lambda i: rd(l, fam, i, lm) == "worse")):
                    cnt = [sum(test(i) for i in instances(fam) if int(i.split(",")[pos if fam != "SM" else 0]) == x)
                           for x in DRAWS]
                    dr.append([l, fam, axis, mname] + [fmt(x) for x in cnt])
    out["draws.tsv"] = dr
    parts = []
    for l in LEARNERS:
        for mm in ("recall", "specificity", "loss"):
            m = mm if mm != "loss" else loss_of("second", l)
            sr = [diff2(d, l, "SR", "mean", m) for d in datasets]
            sm = [diff2(d, l, "SM", "mean", m) for d in datasets]
            mr = [diff2(d, l, "MR", "mean", m) for d in datasets]
            up = 1 if WORSE_UP[m] else -1
            exact = all(a == b + c for a, b, c in zip(sr, sm, mr))
            parts.append([l, mm, fmt(len(datasets)), fmt(mean(sr)), fmt(mean(sm)), fmt(mean(mr)),
                          fmt(mean(sm) / mean(sr)), fmt(sum(1 for b, c in zip(sm, mr) if up * b > up * c)),
                          fmt(exact)])
    out["parts.tsv"] = parts
    # first experiment readings with the leave-one-species-out readings that keep each
    first = []
    for l in LEARNERS:
        for fam in FIRST:
            for m in measures_of(l):
                r0 = R[("first", l, fam, "-", m, "all")]["reading"]
                left = [s for s in species if ("first", l, fam, "-", m, s) in R]
                changed = [s for s in left if R[("first", l, fam, "-", m, s)]["reading"] != r0]
                first.append([l, fam, m, r0, fmt(len(left) - len(changed)), ",".join(changed) or "-"])
        r0 = outcome(R[("first", l, "study", "-", "recall", "all")]["reading"],
                     R[("first", l, "study", "-", "log_loss", "all")]["reading"])
        changed = [s for s in species if outcome(R[("first", l, "study", "-", "recall", s)]["reading"],
                                                 R[("first", l, "study", "-", "log_loss", s)]["reading"]) != r0]
        first.append([l, "study", "rule", r0, fmt(len(species) - len(changed)), ",".join(changed) or "-"])

    def first_rule(scope):
        return experiment_rule([outcome(R[("first", l, "study", "-", "recall", scope)]["reading"],
                                        R[("first", l, "study", "-", "log_loss", scope)]["reading"])
                                for l in LEARNERS])
    a0 = first_rule("all")
    changed = [s for s in species if first_rule(s) != a0]
    first.append(["all", "study", "experiment", a0, fmt(len(species) - len(changed)), ",".join(changed) or "-"])
    out["first.tsv"] = first
    S = {}
    S["experiment_SR_cost_shown"] = sum(experiment_rule([oc(l, "SR", i) for l in LEARNERS]) == "cost shown"
                                        for i in instances("SR"))
    # scale
    S["datasets"] = len(datasets)
    S["species"] = len(species)
    S["pairs"] = sum(int(ds[d]["genomes"]) for d in datasets)
    S["resistant_pairs"] = sum(int(ds[d]["resistant"]) for d in datasets)
    S["no_study_pairs"] = sum(int(ds[d]["no_study"]) for d in datasets)
    S["no_study_datasets_none"] = sum(1 for d in datasets if ds[d]["no_study"] == "0")
    S["cohort_median"] = median(int(ds[d]["genomes"]) for d in datasets)
    S["cohort_min"] = min(int(ds[d]["genomes"]) for d in datasets)
    S["cohort_max"] = max(int(ds[d]["genomes"]) for d in datasets)
    S["two_large"] = ",".join(two)
    S["two_large_datasets"] = sum(per_sp[s] for s in two)
    S["species_one_dataset"] = sum(1 for s in species if per_sp[s] == 1)
    S["lineage_datasets"] = len(lin_ds)
    for s in species:
        S[f"datasets_{s}"] = per_sp[s]
    S["multisets"] = math.comb(2 * len(species) - 1, len(species))
    for r in t["totals.tsv"]:
        S["distinct_" + r["key"]] = int(r["value"])
    # eligibility
    el = t["eligibility.tsv"]
    S["release_datasets"] = len(el)
    S["eligible_counts"] = sum(1 for r in el if r["passes"] == "yes")
    S["eligible_entered"] = sum(1 for r in el if r["entered"] == "yes")
    S["fail_studies"] = sum(1 for r in el if r["under_studies"] == "yes")
    S["fail_share"] = sum(1 for r in el if r["under_share"] == "yes")
    S["fail_group"] = sum(1 for r in el if r["over_group"] == "yes")
    S["left_by_folds"] = ",".join(r["dataset"] for r in el if r["passes"] == "yes" and r["entered"] == "no")
    S["variants"] = len(t["variants.tsv"])
    vd = [int(r["datasets"]) for r in t["variants.tsv"]]
    vs = [int(r["species"]) for r in t["variants.tsv"]]
    S["variants_datasets_min"], S["variants_datasets_max"] = min(vd), max(vd)
    S["variants_species_min"], S["variants_species_max"] = min(vs), max(vs)
    # rule variation
    var = defaultdict(list)
    for r in t["variation.tsv"]:
        var[(r["rule"], r["dataset"])].append(Fraction(int(r["changed"]), int(ds[r["dataset"]]["genomes"])))
    for rule in ("A", "B"):
        per = [mean(var[(rule, d)]) for d in datasets]
        S[f"rule_{rule}_median"] = median(per)
        S[f"rule_{rule}_min"] = min(per)
        S[f"rule_{rule}_max"] = max(per)
        S[f"rule_{rule}_unchanged"] = sum(1 for x in per if x == 0)
    # plans: mates, shapes
    P = {(r["dataset"], r["plan"]): r for r in t["plans.tsv"]}
    mt = []
    for kind in KINDS:
        for i in DRAWS:
            rs = [P[(d, f"{kind}-d{i}")] for d in datasets]
            mt.append([kind, fmt(i), fmt(median(num(r["study_mate_share"]) for r in rs)),
                       fmt(median(num(r["median_study_mates"]) for r in rs)),
                       fmt(median(num(r["lineage_mate_share"]) for r in rs)),
                       fmt(median(num(r["median_lineage_mates"]) for r in rs))])
        dm = [mean(Fraction(int(P[(d, f"{kind}-d{i}")]["lineage_mates"]), int(ds[d]["genomes"])) for i in DRAWS)
              for d in datasets]
        sm = [mean(Fraction(int(P[(d, f"{kind}-d{i}")]["study_mates"]), int(ds[d]["genomes"])) for i in DRAWS)
              for d in datasets]
        mt.append([kind, "mean", fmt(median(sm)), "-", fmt(median(dm)), "-"])
    out["mates.tsv"] = mt

    def lshare(d, plan):
        return Fraction(int(P[(d, plan)]["lineage_mates"]), int(ds[d]["genomes"]))
    gap = [mean(lshare(d, f"matched-d{i}") - lshare(d, f"study-d{i}") for i in DRAWS) for d in datasets]
    S["lineage_matched_minus_study_median"] = median(gap)
    S["lineage_matched_minus_study_min"] = min(gap)
    S["lineage_matched_minus_study_max"] = max(gap)
    S["lineage_matched_minus_study_positive"] = sum(1 for x in gap if x > 0)
    gap = [mean(lshare(d, f"random-d{i}") for i in DRAWS) - mean(lshare(d, f"matched-d{i}") for i in DRAWS)
           for d in datasets]
    S["lineage_random_minus_matched_median"] = median(gap)
    S["lineage_random_minus_matched_positive"] = sum(1 for x in gap if x > 0)

    def shape(r, n):
        sizes = [int(r[f"size_{f}"]) for f in range(FOLDS)]
        rs = [int(r[f"resistant_{f}"]) for f in range(FOLDS)]
        share = [Fraction(a, b) for a, b in zip(rs, sizes)]
        prev = Fraction(sum(rs), n)
        return {"largest": Fraction(max(sizes), n), "smallest": Fraction(min(sizes), n),
                "least_training": 1 - Fraction(max(sizes), n), "spread": max(share) - min(share),
                "departure": max(abs(x - prev) for x in share)}
    sh = []
    arms = [("random-d0", lambda p: p == "random-d0"), ("study-d0", lambda p: p == "study-d0")] + \
        [(k, (lambda k: lambda p: p.startswith(k + "-d"))(k)) for k in KINDS]
    for arm, sel in arms:
        recs = [shape(r, int(ds[r["dataset"]]["genomes"])) for r in t["plans.tsv"] if sel(r["plan"])]
        for m in ("largest", "smallest", "least_training", "spread", "departure"):
            xs = [x[m] for x in recs]
            sh.append([arm, m, fmt(median(xs)), fmt(min(xs)), fmt(max(xs))])
    out["shapes.tsv"] = sh
    # distances between draws
    dd = defaultdict(list)
    for r in t["distances.tsv"]:
        dd[(r["kind"], r["dataset"])].append(Fraction(int(r["changed"]), int(ds[r["dataset"]]["genomes"])))
    for kind in KINDS:
        per = [mean(dd[(kind, d)]) for d in datasets]
        S[f"distance_{kind}_median"] = median(per)
        S[f"distance_{kind}_min"] = min(per)
        S[f"distance_{kind}_max"] = max(per)
        S[f"distance_{kind}_under_tenth"] = sum(1 for x in per if x < Fraction(1, 10))
        S[f"distance_{kind}_identical_pairs"] = sum(1 for d in datasets for x in dd[(kind, d)] if x == 0)
    # genomes without a BioProject
    sub = {(r["dataset"], r["plan"], r["learner"], r["subset"]): r for r in t["subsets.tsv"]}
    ns = [d for d in datasets if int(ds[d]["no_study"]) >= NO_STUDY_MIN]
    S["no_study_datasets"] = len(ns)
    S["no_study_median_genomes"] = median(int(sub[(d, "study-d0", "scm", "no-study")]["n"]) for d in ns)
    S["no_study_median_resistant"] = median(int(sub[(d, "study-d0", "scm", "no-study")]["resistant"]) for d in ns)
    S["no_study_under_5_resistant"] = sum(1 for d in ns if int(sub[(d, "study-d0", "scm", "no-study")]["resistant"]) < 5)
    out["nostudy.tsv"] = nostudy_rows(sub, ns)
    # the first experiment per study fold
    fr = defaultdict(list)
    for r in t["folds.tsv"]:
        n = int(ds[r["dataset"]]["genomes"])
        test = int(r["test"])
        x = dict(r)
        x["ratio"] = Fraction(n - test) / Fraction(int(r["train_random_sum"]), test)
        tres = int(r["test_resistant"])
        x["shift"] = abs(Fraction(tres, test) - Fraction(int(ds[r["dataset"]]["resistant"]) - tres, n - test))
        x["d_ll"] = arm_loss(r, "study") - arm_loss(r, "random")
        fr[r["learner"]].append(x)
    base = fr["scm"]
    S["first_folds"] = len(base)
    S["first_folds_train_ge"] = sum(1 for x in base if x["ratio"] >= 1)
    S["first_pairs_train_ge"] = sum(int(x["test"]) for x in base if x["ratio"] >= 1)
    S["first_ratio_median"] = median(x["ratio"] for x in base)
    S["first_ratio_min"] = min(x["ratio"] for x in base)
    S["first_ratio_max"] = max(x["ratio"] for x in base)
    shifts = sorted(x["shift"] for x in base)
    cut, cut2 = shifts[len(shifts) // 3], shifts[2 * len(shifts) // 3]
    S["first_shift_cut"] = cut
    groups = (("train_ge", lambda x: x["ratio"] >= 1), ("train_lt", lambda x: x["ratio"] < 1),
              ("cut", lambda x: x["ratio"] >= 1 and x["shift"] < cut), ("third_low", lambda x: x["shift"] < cut),
              ("third_mid", lambda x: cut <= x["shift"] < cut2), ("third_high", lambda x: x["shift"] >= cut2))
    fc = []
    for l in LEARNERS:
        for gname, cond in groups:
            g = [x for x in fr[l] if cond(x)]
            nn = sum(int(x["test"]) for x in g)
            loss_a = sum(arm_loss(x, "study") * int(x["test"]) for x in g)
            loss_b = sum(arm_loss(x, "random") * int(x["test"]) for x in g)
            tpa, fna = sum(int(x["study_tp"]) for x in g), sum(int(x["study_fn"]) for x in g)
            tpb, fnb = sum(int(x["random_tp"]) for x in g), sum(int(x["random_fn"]) for x in g)
            tna, fpa = sum(int(x["study_tn"]) for x in g), sum(int(x["study_fp"]) for x in g)
            tnb, fpb = sum(int(x["random_tn"]) for x in g), sum(int(x["random_fp"]) for x in g)
            fc.append([l, gname, fmt(len(g)), fmt(len({x["dataset"] for x in g})), fmt(nn),
                       fmt(median(x["d_ll"] for x in g)), fmt((loss_a - loss_b) / nn),
                       fmt(Fraction(tpa, tpa + fna) - Fraction(tpb, tpb + fnb)),
                       fmt(Fraction(tna, tna + fpa) - Fraction(tnb, tnb + fpb)),
                       fmt(median(x["shift"] for x in g)), fmt(median(x["ratio"] for x in g))])
        per = defaultdict(lambda: [0] * 11)
        for x in fr[l]:
            if x["ratio"] >= 1:
                p = per[x["dataset"]]
                p[0] += arm_loss(x, "study") * int(x["test"])
                p[1] += arm_loss(x, "random") * int(x["test"])
                p[2] += int(x["test"])
                p[3] += int(x["study_tp"]); p[4] += int(x["study_fn"])
                p[5] += int(x["random_tp"]); p[6] += int(x["random_fn"])
                p[7] += int(x["study_tn"]); p[8] += int(x["study_fp"])
                p[9] += int(x["random_tn"]); p[10] += int(x["random_fp"])
        dl = [(p[0] - p[1]) / p[2] for p in per.values()]
        drc = [Fraction(p[3], p[3] + p[4]) - Fraction(p[5], p[5] + p[6]) for p in per.values() if p[3] + p[4]]
        S[f"first_ge_datasets_{l}"] = len(dl)
        S[f"first_ge_log_loss_{l}"] = median(dl)
        S[f"first_ge_worse_{l}"] = sum(1 for x in dl if x > 0)
        S[f"first_ge_recall_{l}"] = median(drc)
        dsp = [Fraction(p[7], p[7] + p[8]) - Fraction(p[9], p[9] + p[10]) for p in per.values() if p[7] + p[8]]
        S[f"first_ge_specificity_{l}"] = median(dsp)
    out["foldcut.tsv"] = fc
    # scm: the log loss summed from its scores against the exact form derive() reads
    S["scm_log_loss_per_error"] = PER_ERROR
    gaps = [abs(Fraction(float(r["log_loss"]))
                - loss_value("scm", None, *(int(r[k]) for k in ("tp", "fn", "tn", "fp"))))
            for name in ("fits.tsv", "lineage.tsv") for r in t[name] if r["learner"] == "scm"]
    S["scm_log_loss_largest_gap"] = max(gaps)
    for k in ("study", "random", "matched"):
        for l in ("gbt", "lasso"):
            S[f"null_{l}_{k}"] = ",".join(fmt(sum(V[(d, f"{k}-d{i}", l)]["null"] for d in datasets)) for i in DRAWS)
    out["summary.tsv"] = [[k, fmt(v)] for k, v in S.items()]
    return out


# ----------------------------------------------------------------- record

README = """# Study split record

This directory holds the record of the two experiments with folds by study at revision {revision}. The page that reads it is [`docs/study-split.md`](../../docs/study-split.md). `scripts/study-split/record.py` writes every file here from the reduced input set, and its `--check` mode fails on any byte it would not write.

## Files

| File | Rows | Contents |
|---|---|---|
{table}

## Rules

**Measures.** The script computes every measure from scores and labels. A score at or above {cutoff} calls resistant. Log loss clips the score given to the true class to [{clip}, 1 - {clip}], so a wrong 0 and a wrong 1 lose the same amount. Recall, specificity, balanced accuracy and the error rate are exact fractions of the confusion counts.

**Loss.** The first plan's rule reads log loss for every learner. The second plan's rule reads the error rate for scm, which scores only 0 or 1, and log loss for the others.

**scm.** A 0 or 1 score makes scm's clipped log loss a right call's loss plus {per_error} times its error rate. Every contrast, interval and count takes scm's log loss in this form, which keeps the error rate's ties. `fits.tsv` gives the log loss summed from the scores, and `summary.tsv` gives the largest gap between the two.

**Intervals.** Each species interval is exact: it takes every multiset of the species drawn with replacement, at the weight of its multinomial count. The interval holds species resampling only.

**Derived.** The script's derive() writes {derived} from the row files.

## Guards

`--check` rebuilds every file from the reduced input set and fails on any byte that differs. Without the archive, `tests/study_split` guards each file in one or more of these ways.

| Guard | Files |
|---|---|
| Rebuilt | {derived}: the suite rebuilds each from the row files, byte for byte |
| Agreed | {agreed}: the suite checks their cells against each other and the derived files |
| Quoted | {quoted}: the suite rebuilds the page's lines from them |
| Pinned | {pinned}: the suite holds each one's sha256 |
| Whole | `README.md` is what the script writes, and `SHA256SUMS` gives each file's sha256 |

## Plans

The two plans stay in the reduced input set as `first/PLAN.md` and `second/PLAN.md`. `planfiles.tsv` gives each one's sha256, the files its hash list holds and its hash time, and `--check` hashes the plans again. `settings.tsv` gives `first_design_key` and `second_design_key` as the fold builders hashed them into every keyed order; they stay as they were so that the folds rebuild.

## Commands

```sh
python3 -I scripts/study-split/record.py --inputs <dir> --out reference/study-split
python3 -I scripts/study-split/record.py --inputs <dir> --check reference/study-split
python3 -I -m unittest discover -s tests/study_split -t .
```

`<dir>` is the reduced input set, the directory that `study-split-inputs.tar.zst` unpacks to. `inputs.sha256` lists the sha256 of each of its files.
"""


def readme(files, settings):
    rows = []
    for name in sorted(DESCRIBE):
        count = files[name].count(b"\n") - 1
        rows.append(f"| `{name}` | {count:,} | {DESCRIBE[name]} |")
    names = lambda xs: ", ".join(f"`{x}`" for x in xs)
    return README.format(revision=settings["release_revision"], table="\n".join(rows), cutoff=settings["cutoff"],
                         clip=settings["clip"], per_error=f"{float(PER_ERROR):.2f}", derived=names(DERIVED),
                         agreed=names(AGREED), quoted=names(QUOTED), pinned=names(PINNED)).encode()


def input_hashes(inputs):
    I = Path(inputs)
    lines = []
    for p in sorted(x for x in I.rglob("*") if x.is_file()):
        lines.append(f"{sha256_bytes(p.read_bytes())}  {p.relative_to(I).as_posix()}\n")
    return "".join(lines).encode()


def genome_ids_in(files, genomes):
    """The genome ids of the input tables that stand as a token in any of files."""
    found = set()
    for b in files.values():
        tokens = set(re.findall(r"[A-Za-z0-9_.]+", b.decode("utf-8")))
        found |= (tokens | {x.rstrip(".") for x in tokens}) & genomes
    return found


def build(inputs, jobs=1):
    tables, genomes = build_rows(inputs, jobs)
    files = {name: tsv(name, rows) for name, rows in tables.items()}
    parsed = {name: parse(files[name].decode(), name) for name in ROW_FILES}
    for name, rows in derive(parsed, jobs).items():
        files[name] = tsv(name, rows)
    files["inputs.sha256"] = input_hashes(inputs)
    # The record holds no genome id: checks.tsv gives the count over every other file, and the build stops above 0.
    found = genome_ids_in(files, genomes)
    tables["checks.tsv"] = sorted(tables["checks.tsv"] + [["genome_ids_in_record", fmt(len(found))]])
    files["checks.tsv"] = tsv("checks.tsv", tables["checks.tsv"])
    need(not found, f"genome ids of the inputs stand in the record: {len(found)}, {min(found or [''])} first")
    settings = {r[0]: r[1] for r in tables["settings.tsv"]}
    files["README.md"] = readme(files, settings)
    total = 0
    for name, b in files.items():
        total += len(b)
    need(total < 4 * 1024 * 1024, f"record of {total} bytes, over 4 MB")
    files["SHA256SUMS"] = "".join(f"{sha256_bytes(files[n])}  {n}\n" for n in sorted(files)).encode()
    need(not genome_ids_in(files, genomes), "a genome id of the inputs stands in the README or SHA256SUMS")
    return files


def on_disk(out):
    out = Path(out)
    return {p.relative_to(out).as_posix(): p.read_bytes() for p in sorted(out.rglob("*")) if p.is_file()}


def listed_names(out):
    """The names an earlier record at out lists in its SHA256SUMS, or a stop if out is not a target to write into.

    A target is a missing or empty directory, or one that holds a SHA256SUMS this script could have written. Any
    other directory stops the script before it writes or deletes anything."""
    out = Path(out)
    if not out.exists():
        return set()
    need(out.is_dir(), f"{out}: not a directory")
    if not any(out.iterdir()):
        return set()
    sums = out / "SHA256SUMS"
    need(sums.is_file(), f"{out}: not empty and holds no SHA256SUMS, so no earlier record of this script; "
                         "nothing written")
    names = set()
    for line in sums.read_text(encoding="utf-8").split("\n"):
        if not line:
            continue
        digest, _, name = line.partition("  ")
        need(re.fullmatch(r"[0-9a-f]{64}", digest) is not None and name != "",
             f"{sums}: a line that is not a sha256 and a name; nothing written")
        parts = name.split("/")
        need(not name.startswith("/") and "" not in parts and "." not in parts and ".." not in parts,
             f"{sums}: {name}: not a name inside the record; nothing written")
        names.add(name)
    return names


def write_record(out, files, old):
    """Write files under out, then delete each name of the earlier SHA256SUMS that this build does not write."""
    out = Path(out)
    for name, b in files.items():
        p = out / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(b)
    for name in sorted(old - set(files)):
        p = out / name
        if p.is_file():
            p.unlink()
            # A directory the earlier record made and that the deletion left empty goes too.
            d = p.parent
            while d != out and not any(d.iterdir()):
                d.rmdir()
                d = d.parent


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build or check the study-split record.")
    ap.add_argument("--inputs", required=True, help="the reduced input set")
    ap.add_argument("--jobs", type=int, default=1, help="worker processes, at most 12")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out", help="write the record here")
    g.add_argument("--check", help="compare the record here with a fresh build")
    a = ap.parse_args(argv)
    try:
        need(1 <= a.jobs <= 12, "--jobs must be 1 to 12")
        need(Path(a.inputs).is_dir(), f"{a.inputs}: no such directory")
        # The target is read before the build, so that a wrong one costs nothing and loses nothing.
        old = listed_names(a.out) if a.out else set()
        files = build(a.inputs, a.jobs)
    except Stop as e:
        print(f"stop: {e}", file=sys.stderr)
        return 2
    if a.out:
        write_record(a.out, files, old)
        print(f"wrote {len(files)} files to {a.out}")
        return 0
    have = on_disk(a.check)
    bad = []
    for name in sorted(set(have) | set(files)):
        if name not in have:
            bad.append(f"{name}: missing")
        elif name not in files:
            bad.append(f"{name}: not written by the script")
        elif have[name] != files[name]:
            bad.append(f"{name}: differs")
    for line in bad:
        print(line, file=sys.stderr)
    if bad:
        return 1
    print(f"{len(files)} files equal a fresh build")
    return 0


if __name__ == "__main__":
    sys.exit(main())
