#!/usr/bin/env python3
"""Write and check reference/kover-comparison, the record of the comparison of Kover with GRM's scm learner.

On 15 Nov 2026 Kover 2.0.7 was run on the 12 smallest eligible datasets of release
rev14-274-20261112, one per species, on the release's own lineage folds, and its held-out
calls were compared with the released lineage calls of GRM's scm learner. A review then
recomputed every figure and narrowed the reading. This script writes the record of that
comparison from the reduced files the comparison kept, which hold no sequence, no rule text
and no importance value, and from the review's tie files.

Modes:

    record.py --source KOVER --review REVIEW --out DIR     write the record from the sources
    record.py --source KOVER --review REVIEW --check DIR   fail unless DIR is what the sources give
    record.py --check DIR                                  fail unless DIR's derived files, README and
                                                           SHA256SUMS are what its row files give

KOVER is the comparison's output directory (PLAN.md, PLAN.sha256, results/) and REVIEW the
review's (PLAN.sha256, results/). Neither is public, so only the last mode runs without them.
The script reads no Kover run file, fits nothing and downloads nothing. Standard library only;
run it as python3 -I.
"""

import argparse
import hashlib
import io
import json
import math
import os
import pathlib
import random
import re
import sys

EPS = 1e-4  # the program's clip of a 0 or 1 call before log loss
T975_DF4 = 2.7764451051977987  # t quantile 0.975 on 4 degrees of freedom
N_BOOT = 10000
ZERO = 1e-12  # a difference within this counts as equal, as the comparison counted it
# The seed the comparison's code gives each measure's bootstrap (BOOT_SEED + i, agreement + 99);
# the error rate takes the log loss seed, as the review's conversion did; the tie readings take 20261115.
SEEDS = {
    "balanced_accuracy": 20261115, "recall": 20261116, "specificity": 20261117, "log_loss": 20261118,
    "error_rate": 20261118, "rules": 20261119, "agreeing": 20261214,
    "ties_equal": 20261115, "ties_lowest": 20261115, "ties_highest": 20261115,
}
MEASURES = ["balanced_accuracy", "recall", "specificity", "log_loss", "error_rate", "rules", "agreeing",
            "ties_equal", "ties_lowest", "ties_highest"]
READING = {"balanced_accuracy": "fixed", "recall": "fixed", "specificity": "fixed", "log_loss": "fixed",
           "error_rate": "exploratory", "rules": "fixed", "agreeing": "fixed",
           "ties_equal": "exploratory", "ties_lowest": "exploratory", "ties_highest": "exploratory"}

COLUMNS = {
    "datasets.tsv": ["order", "dataset", "species", "genomes", "resistant", "susceptible", "lineages",
                     "lineage_folds", "genomes_listed_in_data_root", "fold_sizes"],
    "folds.tsv": ["dataset", "fold", "test_genomes", "resistant", "susceptible",
                  "grm_tp", "grm_fn", "grm_tn", "grm_fp", "grm_rules",
                  "kover_tp", "kover_fn", "kover_tn", "kover_fp", "kover_rules", "kover_kind", "kover_p", "agreeing"],
    "ties.tsv": ["dataset", "fold", "rules_tied", "rules_kept", "patterns", "first_index_equals_run",
                 "first_index_test_balanced_accuracy", "tied_test_balanced_accuracy"],
    "steps.tsv": ["dataset", "step", "kind", "p", "grm_covered", "grm_lost", "kover_covered", "kover_lost",
                  "grm_listed_tied", "counts_equal", "pattern_equal", "kover_equals_a_listed_tied"],
    "cohort.tsv": ["measure", "chosen", "eligible"],
    "order.tsv": ["utc", "event", "file", "sha256"],
    "settings.tsv": ["key", "value"],
    "mapping.tsv": ["setting", "grm", "kover", "same"],
    "per_dataset.tsv": ["dataset", "genomes", "folds", "balanced_accuracy", "balanced_accuracy_t_low",
                        "balanced_accuracy_t_high", "recall", "specificity", "log_loss", "error_rate", "rules",
                        "agreeing_share"],
    "summary.tsv": ["measure", "reading", "value", "low", "high", "seed", "datasets_higher", "datasets_equal",
                    "datasets_lower", "folds_higher", "folds_equal", "folds_lower"],
    "totals.tsv": ["key", "value"],
}
ROW_FILES = ["datasets.tsv", "folds.tsv", "ties.tsv", "steps.tsv", "cohort.tsv", "order.tsv", "settings.tsv",
             "mapping.tsv"]
DERIVED = ["per_dataset.tsv", "summary.tsv", "totals.tsv"]

# The Kover version and settings, in the plan's words where the plan gives them (checked against
# PLAN.md in the source modes).
SETTINGS = [
    ("grm_release", "rev14-274-20261112"),
    ("evaluation_revision", "14"),
    ("release_datasets", "274"),
    ("dataset_rule", "for each species, the dataset with the fewest genomes whose scm row has a lineage entry, ties by dataset id, kept at 150 genomes or fewer"),
    ("reading_rule", "a measured difference when the 95% interval over datasets lies wholly above or below 0, no difference shown when it contains 0, and no equivalence claimed either way"),
    ("grm_learner", "scm"),
    ("grm_calls", "the released lineage scm row: one 0 or 1 call per genome, from the fold that did not train on it"),
    ("folds", "the release's lineage split test_fold: fold k tests the genomes with test_fold == k and trains on the others"),
    ("kover_version", "Kover 2.0.7 (core-2.0.7, cli-2.0.1)"),
    ("kmer_tools", "GATB-core 1.2.2"),
    ("kover_dataset", "kover dataset create from-contigs --kmer-size 31 --n-cpu 4"),
    ("kover_split", "kover dataset split --id lineage_fold<k> --train-ids <train> --test-ids <test> --folds 3 --random-seed 42"),
    ("kover_learn", "kover learn scm --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 --max-rules 6 --hp-choice cv --n-cpu 1 --random-seed 42"),
    ("kover_singletons", "left out (--singleton-kmers not given)"),
    ("kover_calls", "a genome in test_correct gets its label, one in test_errors the other class"),
    ("log_loss_clip", "0.0001"),
    ("bootstrap", "10000 resamples of whole datasets, 95% percentile interval"),
    ("t_interval", "per dataset, over its 5 fold differences, t quantile 0.975 on 4 degrees of freedom"),
    ("equal_within", "1e-12"),
]
# The text each setting must match in PLAN.md (backticks removed).
SETTINGS_IN_PLAN = {
    "grm_release": "(rev14-274-20261112)",
    "release_datasets": "(268 of 274)",
    "dataset_rule": ("For each species, the one with the fewest genomes (ties: the dataset id in alphabetical order).",
                     "Keep those with at most 150 genomes."),
    "reading_rule": ("If the 95% interval of the mean balanced accuracy difference over datasets lies wholly above or below 0",
                     "If it contains 0, the report says this comparison does not show a difference",
                     "No equivalence is claimed either way."),
    "grm_learner": "The released lineage scm row",
    "grm_calls": "scores (one 0 or 1 call per genome, each from the fold that did not train on it)",
    "folds": "The report's lineage split test_fold: fold k tests the genomes with test_fold == k and trains on the others.",
    "kover_version": "Kover 2.0.7 (core-2.0.7, cli-2.0.1)",
    "kmer_tools": "GATB-core 1.2.2",
    "kover_dataset": "kover dataset create from-contigs --kmer-size 31 --n-cpu 4",
    "kover_split": "kover dataset split --id lineage_fold<k> --train-ids <train> --test-ids <test> --folds 3 --random-seed 42",
    "kover_learn": "kover learn scm --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 --max-rules 6 --hp-choice cv --n-cpu 1 --random-seed 42",
    "kover_singletons": "(--singleton-kmers not given)",
    "kover_calls": "a genome in test_correct gets its label, one in test_errors the other class",
    "log_loss_clip": "clipped to [1e-4, 1 - 1e-4]",
    "bootstrap": "10,000 bootstrap resamples of whole datasets",
}
# The plan's table of how the settings map, row by row (backticks removed). The p values row
# names the values where the plan says "the 7 values above".
MAPPING = [
    ("k-mer", "canonical 31-mers", "canonical 31-mers (DSK)", "yes"),
    ("Features offered", "training patterns with support 3 to n_train - 3, one per distinct training pattern",
     "every k-mer in at least 2 genomes of the whole dataset; presence and absence rules", "no"),
    ("Model-size limit", "6 rules", "--max-rules 6", "yes"),
    ("Kinds", "conjunction, disjunction", "both", "yes"),
    ("p values", "0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0", "the same 7 (Kover's default list adds 999999.0; not used)", "yes"),
    ("Utility", "covered - p * lost", "the same", "yes"),
    ("Rule tie", "fewer lost, lower index, presence first",
     "equal within np.isclose, then lowest training error of the rule alone, then the first in Kover's k-mer order", "no"),
    ("Inner folds", "3, whole lineages placed to balance the classes", "3, drawn at random by Kover, not grouped, not stratified", "no"),
    ("Selection criterion", "mean balanced error", "mean error rate", "no"),
    ("Selection tie", "keep earlier (conjunction, smaller p, fewer rules)",
     "fewer rules, then p nearer 1, then arrival order (--n-cpu 1 fixes it)", "no"),
    ("Fewest rules", "1", "0 (an empty model) is allowed", "no"),
    ("Seeds", "42 plus 101 times the fold", "42", "number only; the draws differ"),
]
MAPPING_REWORDED = {"p values": ("the 7 values above", None)}
# review/imbalance.txt, line by line: the 12 chosen datasets against all 268 eligible ones.
COHORT = [
    ("datasets", r"chosen12 datasets (\d+) \|", r"eligible datasets (\d+) \|"),
    ("resistant_share_median", r"chosen12 .*?resistant share median ([0-9.]+) ", r"eligible .*?resistant share median ([0-9.]+) "),
    ("resistant_majority", r"chosen12 .*?resistant-majority (\d+) ", r"eligible .*?resistant-majority (\d+) "),
    ("lineages_median", r"chosen12 .*?lineages median ([0-9.]+)$", r"eligible .*?lineages median ([0-9.]+)$"),
    ("genome_dataset_entries", r"genome-dataset entries: eligible \d+, chosen (\d+) ", r"genome-dataset entries: eligible (\d+),"),
    ("genomes_median", r"n median: eligible [0-9.]+, chosen ([0-9.]+)$", r"n median: eligible ([0-9.]+),"),
    ("species", r"species: eligible \d+, with a chosen dataset (\d+),", r"species: eligible (\d+),"),
    ("grm_lineage_scm_balanced_accuracy_median", r"balanced accuracy median: chosen ([0-9.]+),", r"balanced accuracy median: chosen [0-9.]+, eligible ([0-9.]+)$"),
]
PLAN_FILES = ["PLAN.md", "results/datasets.tsv", "scripts/select_datasets.py", "scripts/build_kmer_tools.sh",
              "scripts/smoke_synthetic.sh"]


class Stop(Exception):
    """An input or a record the script will not write from."""


# ---------------------------------------------------------------- tables

def parse(text):
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        cells = line.split("\t")
        if len(cells) != len(head):
            raise Stop("a row has %d cells where the header has %d" % (len(cells), len(head)))
        rows.append(dict(zip(head, cells)))
    return head, rows


def tsv(name, rows):
    cols = COLUMNS[name]
    out = ["\t".join(cols)]
    for r in rows:
        out.append("\t".join(str(r[c]) for c in cols))
    return "\n".join(out) + "\n"


def f6(x):
    s = "%.6f" % x
    return "0.000000" if s == "-0.000000" else s


def yes(flag):
    return "yes" if flag else "no"


def num(text):
    """'93.0' as '93'; other numbers as written."""
    return text[:-2] if text.endswith(".0") else text


# ---------------------------------------------------------------- measures

def ratio(a, b):
    return 0.0 if b == 0 else a / b


def measures(tp, fn, tn, fp):
    n = tp + fn + tn + fp
    sens, spec = ratio(tp, tp + fn), ratio(tn, tn + fp)
    errors = fn + fp
    return {
        "balanced_accuracy": (sens + spec) / 2.0,
        "recall": sens,
        "specificity": spec,
        "log_loss": (errors * -math.log(EPS) + (n - errors) * -math.log(1.0 - EPS)) / n,
        "error_rate": errors / n,
    }


def quantile(xs, q):
    # linear interpolation between order statistics (type 7), as the comparison's code
    n = len(xs)
    h = (n - 1) * q
    lo = int(math.floor(h))
    hi = min(lo + 1, n - 1)
    return xs[lo] + (h - lo) * (xs[hi] - xs[lo])


def boot_mean(values, seed):
    """The 95% percentile interval of the mean over 10,000 resamples of whole datasets, as the comparison drew them."""
    rng = random.Random(seed)
    n = len(values)
    means = []
    for _ in range(N_BOOT):
        s = 0.0
        for _ in range(n):
            s += values[rng.randrange(n)]
        means.append(s / n)
    means.sort()
    return quantile(means, 0.025), quantile(means, 0.975)


def tied_values(cell):
    """'0.5:67;0.916667:93' as [(0.5, 67), (0.916667, 93)]."""
    out = []
    for part in cell.split(";"):
        ba, count = part.split(":")
        out.append((float(ba), int(count)))
    return out


def hel(values):
    return (sum(1 for v in values if v > ZERO), sum(1 for v in values if abs(v) <= ZERO),
            sum(1 for v in values if v < -ZERO))


# ---------------------------------------------------------------- checks of the row files

def validate(tables):
    """Fail unless the row files agree with each other and with the script's settings."""
    datasets = tables["datasets.tsv"]
    folds = tables["folds.tsv"]
    order = [d["dataset"] for d in datasets]
    if [int(d["order"]) for d in datasets] != list(range(1, len(datasets) + 1)):
        raise Stop("datasets.tsv: order is not 1 to %d" % len(datasets))
    if len(set(order)) != len(order):
        raise Stop("datasets.tsv: a dataset is listed twice")
    want = [(d["dataset"], k) for d in datasets for k in range(int(d["lineage_folds"]))]
    have = [(f["dataset"], int(f["fold"])) for f in folds]
    if have != want:
        raise Stop("folds.tsv: the folds are not each dataset's lineage folds in the datasets' order")
    for d in datasets:
        rows = [f for f in folds if f["dataset"] == d["dataset"]]
        sizes = [int(f["test_genomes"]) for f in rows]
        if ",".join(str(s) for s in sizes) != d["fold_sizes"]:
            raise Stop("folds.tsv: %s fold sizes differ from datasets.tsv" % d["dataset"])
        for k, total in (("test_genomes", "genomes"), ("resistant", "resistant"), ("susceptible", "susceptible")):
            if sum(int(f[k]) for f in rows) != int(d[total]):
                raise Stop("folds.tsv: %s %s do not add up to datasets.tsv" % (d["dataset"], k))
        if int(d["resistant"]) + int(d["susceptible"]) != int(d["genomes"]):
            raise Stop("datasets.tsv: %s resistant and susceptible do not add up" % d["dataset"])
    for f in folds:
        where = "%s fold %s" % (f["dataset"], f["fold"])
        n, r, s = int(f["test_genomes"]), int(f["resistant"]), int(f["susceptible"])
        if r + s != n:
            raise Stop("folds.tsv: %s resistant and susceptible do not add up" % where)
        for side in ("grm", "kover"):
            tp, fn, tn, fp = (int(f["%s_%s" % (side, k)]) for k in ("tp", "fn", "tn", "fp"))
            if min(tp, fn, tn, fp) < 0 or tp + fn != r or tn + fp != s:
                raise Stop("folds.tsv: %s %s counts do not add up" % (where, side))
            if int(f["%s_rules" % side]) < 0:
                raise Stop("folds.tsv: %s %s rules below 0" % (where, side))
        if int(f["grm_rules"]) < 1:
            raise Stop("folds.tsv: %s GRM has no rule" % where)
        if not 0 <= int(f["agreeing"]) <= n:
            raise Stop("folds.tsv: %s agreeing out of range" % where)
        if f["kover_kind"] not in ("conjunction", "disjunction"):
            raise Stop("folds.tsv: %s kind %r" % (where, f["kover_kind"]))
    by_fold = {(f["dataset"], int(f["fold"])): f for f in folds}
    one_rule = [k for k, f in by_fold.items() if int(f["kover_rules"]) == 1]
    ties = tables["ties.tsv"]
    keys = [(t["dataset"], int(t["fold"])) for t in ties]
    if keys != [k for k in have if k in set(one_rule)]:
        raise Stop("ties.tsv: the rows are not Kover's one-rule folds in order")
    for t in ties:
        where = "%s fold %s" % (t["dataset"], t["fold"])
        f = by_fold[(t["dataset"], int(t["fold"]))]
        tied = tied_values(t["tied_test_balanced_accuracy"])
        if sum(c for _, c in tied) != int(t["rules_kept"]) or int(t["rules_kept"]) > int(t["rules_tied"]):
            raise Stop("ties.tsv: %s tied rule counts do not add up" % where)
        if t["first_index_equals_run"] != "yes":
            raise Stop("ties.tsv: %s the first tied rule does not give Kover's calls" % where)
        kover_ba = measures(*(int(f["kover_" + k]) for k in ("tp", "fn", "tn", "fp")))["balanced_accuracy"]
        if abs(float(t["first_index_test_balanced_accuracy"]) - kover_ba) > 5e-7:
            raise Stop("ties.tsv: %s first tied rule's balanced accuracy differs from Kover's fold" % where)
        if not any(abs(ba - kover_ba) <= 5e-7 for ba, _ in tied):
            raise Stop("ties.tsv: %s Kover's own result is not among the tied rules" % where)
    for st in tables["steps.tsv"]:
        if st["dataset"] not in order:
            raise Stop("steps.tsv: %s is not one of the datasets" % st["dataset"])
        equal = (st["grm_covered"], st["grm_lost"]) == (st["kover_covered"], st["kover_lost"])
        if st["counts_equal"] != yes(equal):
            raise Stop("steps.tsv: %s step %s counts_equal does not follow from the counts" % (st["dataset"], st["step"]))
    if [(r["key"], r["value"]) for r in tables["settings.tsv"]] != SETTINGS:
        raise Stop("settings.tsv differs from the script's settings")
    if [(r["setting"], r["grm"], r["kover"], r["same"]) for r in tables["mapping.tsv"]] != MAPPING:
        raise Stop("mapping.tsv differs from the plan's table")
    if [r["measure"] for r in tables["cohort.tsv"]] != [c[0] for c in COHORT]:
        raise Stop("cohort.tsv: the measures are not the review's")
    cohort = {r["measure"]: r for r in tables["cohort.tsv"]}
    if int(cohort["datasets"]["chosen"]) != len(datasets) or int(cohort["genome_dataset_entries"]["chosen"]) != sum(int(d["genomes"]) for d in datasets):
        raise Stop("cohort.tsv: the chosen datasets differ from datasets.tsv")
    utc = [r["utc"] for r in tables["order.tsv"]]
    if utc != sorted(utc) or not all(re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", u) for u in utc):
        raise Stop("order.tsv: the times are not UTC times in order")
    plan = [r for r in tables["order.tsv"] if r["event"] == "plan_hashed"]
    if [r["file"] for r in plan] != PLAN_FILES:
        raise Stop("order.tsv: the plan's hash list is not whole")
    first_fit = [r for r in tables["order.tsv"] if r["event"] == "first_kover_fit"]
    if len(first_fit) != 1 or not plan[0]["utc"] < first_fit[0]["utc"]:
        raise Stop("order.tsv: the plan was not hashed before the first Kover fit")
    # The datasets file is the one the plan hashed.
    if plan[1]["sha256"] != hashlib.sha256(tables["_datasets_bytes"]).hexdigest():
        raise Stop("datasets.tsv is not the file the plan hashed")


# ---------------------------------------------------------------- derived files

def derive(tables):
    datasets = tables["datasets.tsv"]
    folds = tables["folds.tsv"]
    ties = {(t["dataset"], int(t["fold"])): t for t in tables["ties.tsv"]}
    per_ds = {m: [] for m in MEASURES}
    fold_diffs = {m: [] for m in ("balanced_accuracy", "recall", "specificity", "log_loss", "error_rate", "rules")}
    rows_ds = []
    for d in datasets:
        ds = d["dataset"]
        rows = [f for f in folds if f["dataset"] == ds]
        diffs = {m: [] for m in fold_diffs}
        alt = {"ties_equal": [], "ties_lowest": [], "ties_highest": []}
        for f in rows:
            g = measures(*(int(f["grm_" + k]) for k in ("tp", "fn", "tn", "fp")))
            kv = measures(*(int(f["kover_" + k]) for k in ("tp", "fn", "tn", "fp")))
            g["rules"], kv["rules"] = int(f["grm_rules"]), int(f["kover_rules"])
            for m in diffs:
                diffs[m].append(kv[m] - g[m])
                fold_diffs[m].append(kv[m] - g[m])
            vals = {"ties_equal": kv["balanced_accuracy"], "ties_lowest": kv["balanced_accuracy"],
                    "ties_highest": kv["balanced_accuracy"]}
            t = ties.get((ds, int(f["fold"])))
            if t is not None and int(t["rules_kept"]) > 1:
                h = tied_values(t["tied_test_balanced_accuracy"])
                vals = {"ties_lowest": min(b for b, c in h), "ties_highest": max(b for b, c in h),
                        "ties_equal": sum(b * c for b, c in h) / float(sum(c for b, c in h))}
            for key in alt:
                alt[key].append(vals[key] - g["balanced_accuracy"])
        n = len(rows)
        row = {"dataset": ds, "genomes": d["genomes"], "folds": n}
        for m in diffs:
            mean = sum(diffs[m]) / n
            per_ds[m].append(mean)
            row[m] = f6(mean)
        ba = diffs["balanced_accuracy"]
        mean = sum(ba) / n
        sd = math.sqrt(sum((x - mean) ** 2 for x in ba) / (n - 1))
        half = T975_DF4 * sd / math.sqrt(n)
        row["balanced_accuracy_t_low"], row["balanced_accuracy_t_high"] = f6(mean - half), f6(mean + half)
        share = sum(int(f["agreeing"]) for f in rows) / int(d["genomes"])
        row["agreeing_share"] = f6(share)
        per_ds["agreeing"].append(share)
        for key in alt:
            per_ds[key].append(sum(alt[key]) / len(alt[key]))
        rows_ds.append(row)

    summary = []
    for m in MEASURES:
        vals = per_ds[m]
        lo, hi = boot_mean(vals, SEEDS[m])
        row = {"measure": m, "reading": READING[m], "value": f6(sum(vals) / len(vals)), "low": f6(lo), "high": f6(hi),
               "seed": SEEDS[m]}
        if m in fold_diffs:
            dh, de, dl = hel(vals)
            fh, fe, fl = hel(fold_diffs[m])
            row.update(datasets_higher=dh, datasets_equal=de, datasets_lower=dl, folds_higher=fh, folds_equal=fe, folds_lower=fl)
        else:
            row.update(datasets_higher="", datasets_equal="", datasets_lower="", folds_higher="", folds_equal="", folds_lower="")
        summary.append(row)

    steps = tables["steps.tsv"]
    tie_rows = tables["ties.tsv"]
    step1 = [s for s in steps if s["step"] == "1"]
    spread = [t for t in tie_rows if int(t["rules_kept"]) > 1
              and max(b for b, _ in tied_values(t["tied_test_balanced_accuracy"])) - min(b for b, _ in tied_values(t["tied_test_balanced_accuracy"])) > 1e-9]
    t_excl = [r for r in rows_ds if float(r["balanced_accuracy_t_low"]) > 0 or float(r["balanced_accuracy_t_high"]) < 0]
    sizes = [int(f["test_genomes"]) for f in folds]
    totals = [
        ("datasets", len(datasets)),
        ("species", len({d["species"] for d in datasets})),
        ("genomes", sum(int(d["genomes"]) for d in datasets)),
        ("genomes_smallest", min(int(d["genomes"]) for d in datasets)),
        ("genomes_largest", max(int(d["genomes"]) for d in datasets)),
        ("folds", len(folds)),
        ("fold_test_genomes_smallest", min(sizes)),
        ("fold_test_genomes_largest", max(sizes)),
        ("test_calls", sum(sizes)),
        ("one_class_folds", sum(1 for f in folds if int(f["resistant"]) in (0, int(f["test_genomes"])))),
        ("kover_one_rule_folds", sum(1 for f in folds if int(f["kover_rules"]) == 1)),
        ("kover_empty_folds", sum(1 for f in folds if int(f["kover_rules"]) == 0)),
        ("kover_conjunction_folds", sum(1 for f in folds if f["kover_kind"] == "conjunction")),
        ("grm_rules_total", sum(int(f["grm_rules"]) for f in folds)),
        ("kover_rules_total", sum(int(f["kover_rules"]) for f in folds)),
        ("ties_first_index_equals_run", sum(1 for t in tie_rows if t["first_index_equals_run"] == "yes")),
        ("ties_more_than_one_rule", sum(1 for t in tie_rows if int(t["rules_kept"]) > 1)),
        ("ties_different_test_results", len(spread)),
        ("datasets_identical_calls", sum(1 for r in rows_ds if r["agreeing_share"] == "1.000000")),
        ("datasets_t_interval_excludes_0", len(t_excl)),
        ("steps", len(steps)),
        ("steps_fits", len({s["dataset"] for s in steps})),
        ("steps_counts_equal", sum(1 for s in steps if s["counts_equal"] == "yes")),
        ("steps_pattern_equal", sum(1 for s in steps if s["pattern_equal"] == "yes")),
        ("step1_pattern_differs", sum(1 for s in step1 if s["pattern_equal"] == "no")),
        ("step1_differs_equals_a_listed_tied", sum(1 for s in step1 if s["pattern_equal"] == "no" and s["kover_equals_a_listed_tied"] == "yes")),
    ]
    return {
        "per_dataset.tsv": tsv("per_dataset.tsv", rows_ds),
        "summary.tsv": tsv("summary.tsv", summary),
        "totals.tsv": tsv("totals.tsv", [{"key": k, "value": v} for k, v in totals]),
    }


README_ROWS = {
    "README.md": "this list",
    "SHA256SUMS": "the sha256 of every other file of this directory",
    "cohort.tsv": "the 12 datasets against all 268 eligible datasets of the release: sizes, class balance, lineages and GRM's median lineage balanced accuracy",
    "datasets.tsv": "the 12 datasets as the plan chose them, byte for byte the file its hash list names: genomes, resistant, susceptible, lineages and test genomes per lineage fold",
    "folds.tsv": "one row per dataset and lineage fold: test genomes, the confusion counts and rule count of each side, Kover's kind and p, and the calls on which the two agree",
    "mapping.tsv": "the plan's table of how GRM's scm settings map to Kover's, with whether each is the same",
    "order.tsv": "the plan's hash list with its time, the first Kover fit, the review's plan and the review's tie summary",
    "per_dataset.tsv": "derived: per dataset, the mean over its folds of Kover minus GRM in each measure, a t interval for balanced accuracy, and the share of calls agreeing",
    "settings.tsv": "the release, the GRM learner, the Kover version and commands, the clip and the bootstrap",
    "steps.tsv": "the fits on all genomes with GRM's final kind and p held for Kover: covered and lost genomes at each greedy step, and whether the two rules' patterns over the genomes are equal",
    "summary.tsv": "derived: per measure, the mean over datasets of Kover minus GRM with its 95% bootstrap interval and seed, and the datasets and folds where Kover is higher, equal and lower",
    "ties.tsv": "exploratory: Kover's one-rule folds, the rules tied at the first step, and the test balanced accuracy each tied rule would give, as counts of rules per value",
    "totals.tsv": "derived: counts the page quotes",
}


def column_lines():
    """Each TSV of the record and its columns, in order."""
    return ["- `%s`: %s." % (name, ", ".join("`%s`" % c for c in COLUMNS[name])) for name in sorted(COLUMNS)]


def readme(files):
    lines = [
        "# Kover comparison record",
        "",
        "This directory is the record behind [`docs/kover-comparison.md`](../../docs/kover-comparison.md): Kover 2.0.7 against the released lineage calls of GRM's `scm` learner on 12 public datasets of release `rev14-274-20261112`. `scripts/kover-comparison/record.py` writes every file of it, and its `--check` mode fails on any byte it would not write.",
        "",
        "It is not `reference/kover-summary.json`. That file holds the 73 historical Kover rows of the legacy application, on other cohorts and random splits, and `grm report` reads it. This record is read by no program.",
        "",
        "## Files",
        "",
        "| File | Rows | Contents |",
        "|---|---|---|",
    ]
    for name in sorted(files, key=lambda x: (x != "README.md", x != "SHA256SUMS", x)):
        rows = "" if not name.endswith(".tsv") else str(files[name].count("\n") - 1)
        lines.append("| `%s` | %s | %s |" % (name, rows, README_ROWS[name]))
    lines += [
        "",
        "Every file is plain text, and each TSV has a header row. All differences are Kover minus GRM. A file holds counts and measures per dataset, per fold, per greedy step or per measure: no genome, no sequence, no rule and no importance value.",
        "",
        "## Columns",
        "",
    ] + column_lines() + [
        "",
        "## Rules",
        "",
        "- **Measures.** From each fold's confusion counts: balanced accuracy is the mean of recall and specificity, and log loss takes each 0 or 1 call as a probability clipped to [0.0001, 0.9999], so it is 0.0001 plus 9.2102 times the error rate.",
        "- **Over datasets.** Each dataset's value is the mean of its 5 fold differences. The figure over datasets is their mean, with a 95% percentile interval from 10,000 resamples of whole datasets at the seed `summary.tsv` gives, drawn as the comparison drew them.",
        "- **Readings.** `fixed` in `summary.tsv` is a measure the plan fixed before any Kover fit. `exploratory` is one chosen later: the error rate, and the three tie readings, which the review chose after reading the results per fold.",
        "- **Ties.** In each of Kover's one-rule folds, `ties_equal` counts every tied first rule as equally likely, and `ties_lowest` and `ties_highest` take the lowest and highest test result among them, bounds chosen by test result. GRM's side is the same in all three.",
        "- **Steps.** `pattern_equal` is an equality of hashes of the two rules' patterns over the genomes, as the comparison computed them; the record holds none of those hashes.",
        "",
        "## Commands",
        "",
        "```sh",
        "python3 -I scripts/kover-comparison/record.py --check reference/kover-comparison",
        "python3 -I scripts/kover-comparison/record.py --source \"$KOVER\" --review \"$REVIEW\" --check reference/kover-comparison",
        "```",
        "",
        "The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second rebuilds every file from the comparison's and the review's output directories, which are not public.",
        "",
    ]
    return "\n".join(lines)


def assemble(row_texts):
    """Every file of the record from the row files' texts."""
    tables = {}
    for name in ROW_FILES:
        head, rows = parse(row_texts[name])
        if head != COLUMNS[name]:
            raise Stop("%s: columns %s" % (name, head))
        tables[name] = rows
    tables["_datasets_bytes"] = row_texts["datasets.tsv"].encode("utf-8")
    validate(tables)
    files = dict(row_texts)
    files.update(derive(tables))
    files["README.md"] = ""
    files["SHA256SUMS"] = ""
    files["README.md"] = readme(files)
    del files["SHA256SUMS"]
    sums = "".join("%s  %s\n" % (hashlib.sha256(files[k].encode("utf-8")).hexdigest(), k) for k in sorted(files))
    files["SHA256SUMS"] = sums
    return files


# ---------------------------------------------------------------- the sources

def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def sha256_file(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def count(value, total, where):
    """An exact count from a 6-decimal share of total."""
    c = round(float(value) * total)
    if abs(c / total - float(value)) > 5e-7 if total else float(value) != 0.0:
        raise Stop("%s: %s of %d is not a count" % (where, value, total))
    return int(c)


def from_sources(kover, review):
    """The row files' texts from the comparison's and the review's output directories."""
    k, r = pathlib.Path(kover), pathlib.Path(review)
    res = k / "results"
    texts = {}
    texts["datasets.tsv"] = read(res / "datasets.tsv")
    _, datasets = parse(texts["datasets.tsv"])

    # The plan's hash list, and the files it names.
    lines = read(k / "PLAN.sha256").rstrip("\n").split("\n")
    m = re.search(r"recorded at (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", lines[0])
    if not m:
        raise Stop("PLAN.sha256: no time in its first line")
    plan_time = m.group(1)
    listed = [tuple(x.split("  ", 1)) for x in lines[1:]]
    if [name for _, name in listed] != PLAN_FILES:
        raise Stop("PLAN.sha256: the files differ from the plan's list")
    for digest, name in listed[:2]:
        if sha256_file(k / name) != digest:
            raise Stop("PLAN.sha256: %s is not the file the plan hashed" % name)
    order = [{"utc": plan_time, "event": "plan_hashed", "file": name, "sha256": digest} for digest, name in listed]
    m = re.match(r"run_all start (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", read(res / "run_all.log"))
    if not m:
        raise Stop("run_all.log: no start time in its first line")
    order.append({"utc": m.group(1), "event": "first_kover_fit", "file": "results/run_all.log", "sha256": ""})
    rlines = read(r / "PLAN.sha256").rstrip("\n").split("\n")
    m = re.search(r"recorded at (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", rlines[0])
    digest, name = rlines[1].split("  ", 1)
    if not m or name != "PLAN.md" or sha256_file(r / "PLAN.md") != digest:
        raise Stop("the review's PLAN.sha256 does not hold")
    order.append({"utc": m.group(1), "event": "review_plan_hashed", "file": "PLAN.md", "sha256": digest})
    summary_text = read(res / "review" / "fold_ties_summary_exploratory.txt")
    m = re.match(r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)\n", summary_text)
    if not m:
        raise Stop("fold_ties_summary_exploratory.txt: no time in its first line")
    order.append({"utc": m.group(1), "event": "review_tie_summary", "file": "results/fold_ties_summary_exploratory.txt", "sha256": ""})
    order.sort(key=lambda x: x["utc"])
    texts["order.tsv"] = tsv("order.tsv", order)

    # The plan's words.
    plan = read(k / "PLAN.md").replace("`", "")
    for key, wanted in SETTINGS_IN_PLAN.items():
        for text in (wanted,) if isinstance(wanted, str) else wanted:
            if text not in plan:
                raise Stop("PLAN.md does not give %s as %r" % (key, text))
    texts["settings.tsv"] = tsv("settings.tsv", [{"key": a, "value": b} for a, b in SETTINGS])
    for setting, grm, kv, same in MAPPING:
        row = [x for x in plan.split("\n") if x.startswith("| %s |" % setting)]
        if len(row) != 1:
            raise Stop("PLAN.md: no single row for %s" % setting)
        cells = [c.strip() for c in row[0].strip().strip("|").split("|")]
        g_plan = MAPPING_REWORDED.get(setting, (grm, None))[0]
        if cells != [setting, g_plan, kv, same]:
            raise Stop("PLAN.md: the row for %s differs" % setting)
    texts["mapping.tsv"] = tsv("mapping.tsv", [dict(zip(COLUMNS["mapping.tsv"], row)) for row in MAPPING])

    # Folds, from the comparison's per-fold table.
    _, per_fold = parse(read(res / "per_fold.tsv"))
    folds = []
    for f in per_fold:
        where = "per_fold.tsv %s fold %s" % (f["dataset"], f["fold"])
        n, rr, ss = int(f["n_test"]), int(f["resistant"]), int(f["susceptible"])
        row = {"dataset": f["dataset"], "fold": f["fold"], "test_genomes": n, "resistant": rr, "susceptible": ss}
        for side in ("grm", "kover"):
            tp = count(f[side + "_recall"], rr, where)
            tn = count(f[side + "_specificity"], ss, where)
            row.update({side + "_tp": tp, side + "_fn": rr - tp, side + "_tn": tn, side + "_fp": ss - tn})
            mm = measures(tp, rr - tp, tn, ss - tn)
            for key in ("balanced_accuracy", "log_loss"):
                if abs(mm[key] - float(f[side + "_" + key])) > 5e-7:
                    raise Stop("%s: %s %s does not follow from the counts" % (where, side, key))
        row["grm_rules"] = count(f["grm_rules"], 1, where)
        row["kover_rules"] = int(f["kover_rules"])
        row["kover_kind"] = f["kover_kind"]
        row["kover_p"] = repr(float(f["kover_p"]))
        row["agreeing"] = count(f["agreement"], n, where)
        folds.append(row)
    texts["folds.tsv"] = tsv("folds.tsv", folds)

    # Ties, from the review's per-dataset tie files.
    ties = []
    for d in datasets:
        data = json.loads(read(res / "review" / ("fold_ties_%s.json" % d["dataset"])))
        if data["dataset"] != d["dataset"]:
            raise Stop("fold_ties_%s.json names another dataset" % d["dataset"])
        for t in data["folds"]:
            if "skipped" in t:
                continue
            ties.append({
                "dataset": d["dataset"], "fold": t["fold"],
                "rules_tied": t["rules_tied_at_best_utility"], "rules_kept": t["rules_after_risk_tiebreak"],
                "patterns": t["distinct_training_patterns_among_kept"],
                "first_index_equals_run": yes(t["first_index_reproduces_kover_calls"]),
                "first_index_test_balanced_accuracy": repr(float(t["first_index_test_ba"])),
                "tied_test_balanced_accuracy": ";".join("%r:%d" % (float(b), c) for b, c in t["test_ba_of_kept_rules"]),
            })
    texts["ties.tsv"] = tsv("ties.tsv", ties)

    # Steps of the fits on all genomes, from the check with settings held equal.
    _, check = parse(read(res / "final_check.tsv"))
    _, tied = parse(read(res / "final_tied.tsv"))
    tied = {(t["dataset"], t["step"]): t for t in tied}
    steps = []
    for c in check:
        t = tied[(c["dataset"], c["step"])]
        if t["kover_equals_grm_chosen"] != c["pattern_equal"]:
            raise Stop("final_check.tsv and final_tied.tsv disagree on %s step %s" % (c["dataset"], c["step"]))
        steps.append({
            "dataset": c["dataset"], "step": c["step"], "kind": c["kind"], "p": c["p"],
            "grm_covered": c["grm_covered"], "grm_lost": c["grm_lost"], "kover_covered": c["kover_covered"],
            "kover_lost": c["kover_lost"], "grm_listed_tied": t["grm_listed_tied"],
            "counts_equal": yes(c["counts_equal"] == "True"), "pattern_equal": yes(c["pattern_equal"] == "True"),
            "kover_equals_a_listed_tied": yes(t["kover_equals_a_listed_tied"] == "True"),
        })
    texts["steps.tsv"] = tsv("steps.tsv", steps)

    # The 12 against the rest, from the review's imbalance check.
    im = read(res / "review" / "imbalance.txt").rstrip("\n").split("\n")
    cohort = []
    for key, chosen_re, eligible_re in COHORT:
        got = []
        for pattern in (chosen_re, eligible_re):
            hits = [mm.group(1) for line in im for mm in [re.search(pattern, line)] if mm]
            if len(hits) != 1:
                raise Stop("imbalance.txt: no single value for %s" % key)
            got.append(num(hits[0]))
        cohort.append({"measure": key, "chosen": got[0], "eligible": got[1]})
    texts["cohort.tsv"] = tsv("cohort.tsv", cohort)
    return texts


def from_record(directory):
    d = pathlib.Path(directory)
    texts = {}
    for name in ROW_FILES:
        path = d / name
        if not path.is_file():
            raise Stop("%s is missing" % name)
        texts[name] = read(path)
    return texts


# ---------------------------------------------------------------- main

def check(files, directory):
    d = pathlib.Path(directory)
    on_disk = {p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file()}
    bad = sorted(on_disk ^ set(files))
    for name in sorted(files):
        path = d / name
        if path.is_file() and path.read_bytes() != files[name].encode("utf-8"):
            bad.append(name)
    return sorted(set(bad))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source", help="the comparison's output directory")
    ap.add_argument("--review", help="the review's output directory")
    group = ap.add_mutually_exclusive_group(required=True)
    group.add_argument("--out", help="write the record here")
    group.add_argument("--check", help="fail unless this directory is the record")
    args = ap.parse_args(argv)
    if bool(args.source) != bool(args.review):
        ap.error("--source and --review go together")
    if args.out and not args.source:
        ap.error("--out needs --source and --review")
    try:
        texts = from_sources(args.source, args.review) if args.source else from_record(args.check)
        files = assemble(texts)
    except (Stop, OSError, KeyError, ValueError, IndexError) as exc:
        print("record.py: stopped: %s" % exc, file=sys.stderr)
        return 1
    if args.out:
        out = pathlib.Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        for name, text in files.items():
            (out / name).write_bytes(text.encode("utf-8"))
        print("record.py: wrote %d files to %s" % (len(files), out))
        return 0
    bad = check(files, args.check)
    for name in bad:
        print("record.py: differs: %s" % name, file=sys.stderr)
    if bad:
        return 1
    print("record.py: %d files as written" % len(files))
    return 0


if __name__ == "__main__":
    sys.exit(main())
