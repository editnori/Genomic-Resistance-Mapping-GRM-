#!/usr/bin/env python3
"""Write and check reference/dense-network, the record of the 15 Nov experiment that fitted one dense network.

On 15 Nov 2026 one dense feed-forward network, at settings a plan fixed before any result, was fitted on
24 public datasets of release rev14-274-20261112, on the release's own lineage and random folds and on
the forest's screened features, and its held-out calls were compared with the saved calls of GRM's five
learners. The plan's rule gave "no gain established", and a review supported that reading and no wider
one. This script writes the record of that experiment from the experiment's summary file, which holds
measures per dataset and learner and no genome, and checks every figure it derives against the figures
that file gives.

Modes:

    record.py --source RUN --review REVIEW --out DIR     write the record from the sources
    record.py --source RUN --review REVIEW --check DIR   fail unless DIR is what the sources give
    record.py --check DIR                                 fail unless DIR's derived files, README and
                                                          SHA256SUMS are what its row files give

RUN is the experiment's output directory (PLAN.md, REPORT.md, results/ba2-summary.json,
results/stage1-datasets.json, results/gates.json) and REVIEW the review's (REPORT.md). Neither is
public, so only the last mode runs without them. The script fits nothing and downloads nothing.
Standard library only; run it as python3 -I.
"""

import argparse
import hashlib
import json
import pathlib
import random
import re
import statistics
import sys

SPLITS = ["lineage", "random"]
FIVE = ["scm", "cart", "rf", "gbt", "lasso"]
COMPARATORS = FIVE + ["best5", "site_best", "mean4", "baseline", "catalog"]
LEARNERS = ["candidate"] + COMPARATORS
MEASURES = ["balanced_accuracy", "recall", "log_loss", "specificity"]
# The measures the experiment's report gives figures over datasets for; specificity's are the record's own.
REPORTED = ["balanced_accuracy", "recall", "log_loss"]
SOURCE_KEY = {"balanced_accuracy": "ba", "recall": "recall", "log_loss": "log_loss", "specificity": "specificity"}
B = 2000  # resamples of species
SEED = 20261115
MARGIN = 0.02
RECALL_FLOOR = -0.01
BUDGET_CPU_H = 8.0
BUDGET_RSS_KB = 12 * 1024 * 1024
CLOSE = 1e-12  # a derived figure differs from the experiment's own when they part by more than this

COLUMNS = {
    "datasets.tsv": ["order", "dataset", "species", "group", "genomes", "resistant", "lineage_gap", "gate_a", "gate_b"],
    "measures.tsv": ["dataset", "split", "learner", "picks", "balanced_accuracy", "recall", "log_loss", "specificity"],
    "cost.tsv": ["key", "value"],
    "order.tsv": ["utc", "event", "sha256"],
    "settings.tsv": ["key", "value"],
    "summary.tsv": ["split", "comparator", "measure", "candidate_median", "comparator_median", "difference_median",
                    "low", "high", "datasets", "species", "positive", "negative"],
    "rule.tsv": ["rule", "statistic", "value", "needs", "holds"],
    "totals.tsv": ["key", "value"],
}
ROW_FILES = ["datasets.tsv", "measures.tsv", "cost.tsv", "order.tsv", "settings.tsv"]
DERIVED = ["summary.tsv", "rule.tsv", "totals.tsv"]
COST_KEYS = ["fits", "processes", "cpu_seconds", "peak_rss_kb"]
ORDER_EVENTS = ["plan_file_time", "fitting_script_last_change", "first_network_result", "analysis_on_first_datasets", "analysis_script_last_change"]

# What the plan fixed, in its words where it gives them (checked against PLAN.md in the source modes).
SETTINGS = [
    ("grm_release", "rev14-274-20261112"),
    ("network", "dense feed-forward: 4,000 inputs, 256 ReLU units, 64 ReLU units, 1 sigmoid output"),
    ("inputs", "per outer fold, the 4,000 rows of the forest's chi-square screen over that fold's training candidates, on its training samples; each sample a 0/1 vector over those rows"),
    ("start", "He-normal, zero biases"),
    ("loss", "binary cross-entropy, unweighted, plus L2 decay 1e-4, biases left out"),
    ("optimiser", "Adam, learning rate 1e-3, beta 0.9 and 0.999, epsilon 1e-8; minibatches of 64 in a fresh shuffle each epoch"),
    ("stopping", "at most 300 epochs; log loss on the fold's boosting holdout, about 1/7 of the training samples, lineage-grouped on the lineage split; stop after 20 epochs without improvement, keep the best epoch"),
    ("repeats", "3 per fold, each with its own seed; a dataset's value is the mean over the 3"),
    ("splits", "the release's lineage split and random split, 5 outer folds each; the five learners were not fitted again"),
    ("calls", "the five's calls are the release's saved out-of-fold scores, called resistant at 0.5 or above; no learner is fitted again"),
    ("reported", "beside the rule and outside it: each of the five, the site's best (chosen on the random split), the mean of the four probability learners, the lineage baseline (lineage identity inside the release's own outer folds), the catalog rule where the release scores it, the random split, log loss and AUROC"),
    ("datasets", "the 12 lineage datasets with the largest gap, and 12 drawn with seed 20261115 from the other 113 with a gap of 0.05 or more"),
    ("gap", "on the lineage split, the balanced accuracy of calls right wherever any of the five is right, minus that of the best of the five"),
    ("comparator", "the best of the five in hindsight: per dataset, the learner with the highest lineage balanced accuracy"),
    ("measures", "balanced accuracy, recall on the resistant class, specificity, log loss with probabilities clipped to [1e-4, 1 - 1e-4], and AUROC, as the plan named them; the experiment's report gives figures over datasets for balanced accuracy, recall and log loss only; the record carries those three, and specificity worked out from the experiment's rows per dataset by the same draw"),
    ("interval", "2000 resamples of the species present, each taking every dataset of each drawn species; the median over the drawn datasets; 2.5% and 97.5% percentiles; seed 20261115"),
    ("rule_1", "median network minus comparator balanced accuracy at least 0.02, and the interval's lower end above 0"),
    ("rule_2", "median network minus baseline balanced accuracy at least the comparator's median minus baseline"),
    ("rule_3", "median network minus comparator recall at least -0.01"),
    ("rule_4", "at most 8 CPU-hours of network fitting and no process above 12 GB peak memory"),
    ("adopt", "only when all four hold; otherwise no gain established"),
]
SETTINGS_IN_PLAN = {
    "grm_release": "release",
    "network": "Layers: 4,000 inputs, 256 ReLU units, 64 ReLU units, 1 sigmoid output.",
    "inputs": ("the forest's chi-square screen of 4,000 rows over that fold's training candidates, on its training samples",
               "Each sample is a 0/1 vector over those rows."),
    "start": ("He-normal", "zero biases."),
    "loss": ("Loss: binary cross-entropy, unweighted, plus L2", "decay 1e-4", "(not the biases)."),
    "optimiser": "Optimiser: Adam, learning rate 1e-3, beta 0.9 and 0.999, epsilon 1e-8; minibatches of 64 in a fresh shuffle each epoch.",
    "stopping": ("at most 300 epochs", "about 1/7 of the training samples, lineage-grouped on the lineage split",
                 "stop after 20 epochs without improvement and keep the", "of the best epoch"),
    "repeats": ("Seeds: 3 repeats per fold.", "The candidate's value per dataset is the mean over its 3 repeats."),
    "splits": ("without refitting the five", "5 folds"),
    "calls": ("All from the release reports' saved out-of-fold scores", "called at >= 0.5. No learner is refitted."),
    "reported": ("Also reported, not in the rule: each of the five;", "chosen on the random split", "the mean of the four probability learners",
                 "lineage identity inside the release's own outer folds", "the catalog rule as the report scores it, where present.",
                 "The random split, log loss, AUROC and every other comparator are reported beside the decision and do not enter it."),
    "datasets": ("the 12 lineage datasets with the largest gap",
                 "12 drawn with random.Random(20261115).sample from the other 113 with a gap at or above 0.05"),
    "comparator": "the learner with the highest lineage balanced accuracy among scm, cart, rf, gbt, lasso",
    "measures": ("balanced accuracy, recall (sensitivity, resistant class), specificity, log loss with probabilities clipped to [1e-4, 1 - 1e-4]",
                 "scripts/ablation/stats.py), AUROC."),
    "interval": ("Interval: 2,000 resamples of the species present (with replacement), each draw taking every dataset of each drawn species; statistic the median over the drawn datasets; 2.5% and 97.5% percentiles; Python random.Random(20261115).",),
    "rule_1": "median(d) >= 0.02 and the interval's lower end > 0;",
    "rule_2": "median(candidate BA - baseline BA) >= median(comparator BA - baseline BA)",
    "rule_3": "median(candidate recall - comparator recall) >= -0.01",
    "rule_4": "at most 8 CPU-hours of candidate fitting for all of stage 1",
    "adopt": ("Adopt only when all four hold;", "Otherwise the result is \"no gain established\""),
}


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
    return "\n".join(["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]) + "\n"


def f6(x):
    s = "%.6f" % x
    return "0.000000" if s == "-0.000000" else s


def yes(flag):
    return "yes" if flag else "no"


def cell(x):
    """A measure as the record keeps it: the shortest text that reads back as the same float."""
    return "" if x is None else repr(float(x))


# ---------------------------------------------------------------- the experiment's statistics

def pctl(xs, q):
    xs = sorted(xs)
    pos = (len(xs) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (pos - lo)


def interval(by_species):
    """The median over datasets and its interval over 2,000 resamples of species, as the experiment drew them."""
    rng = random.Random(SEED)
    sp = sorted(by_species)
    meds = []
    for _ in range(B):
        draw = []
        for _ in sp:
            draw.extend(by_species[sp[rng.randrange(len(sp))]])
        meds.append(statistics.median(draw))
    allv = [v for s in sp for v in by_species[s]]
    return {"median": statistics.median(allv), "lo": pctl(meds, 0.025), "hi": pctl(meds, 0.975), "n": len(allv),
            "species": len(sp), "positive": sum(1 for v in allv if v > 0), "negative": sum(1 for v in allv if v < 0)}


# ---------------------------------------------------------------- checks of the row files

def index(tables):
    """measures.tsv as {(dataset, split, learner): row} with floats."""
    out = {}
    for r in tables["measures.tsv"]:
        key = (r["dataset"], r["split"], r["learner"])
        if key in out:
            raise Stop("measures.tsv: %s %s %s stands twice" % key)
        row = {"picks": r["picks"]}
        for m in MEASURES:
            row[m] = float(r[m]) if r[m] != "" else None
        out[key] = row
    return out


def validate(tables):
    datasets = tables["datasets.tsv"]
    if [int(d["order"]) for d in datasets] != list(range(1, len(datasets) + 1)):
        raise Stop("datasets.tsv: order is not 1 to %d" % len(datasets))
    names = [d["dataset"] for d in datasets]
    if len(set(names)) != len(names):
        raise Stop("datasets.tsv: a dataset stands twice")
    groups = [d["group"] for d in datasets]
    if groups != ["largest_gap"] * groups.count("largest_gap") + ["drawn"] * groups.count("drawn"):
        raise Stop("datasets.tsv: groups are not the largest gaps, then the drawn")
    for d in datasets:
        if not d["dataset"].startswith(d["species"] + "__"):
            raise Stop("datasets.tsv: %s is not a dataset of %s" % (d["dataset"], d["species"]))
        if not 0 < int(d["resistant"]) < int(d["genomes"]):
            raise Stop("datasets.tsv: %s resistant out of range" % d["dataset"])
        if d["gate_a"] not in ("yes", "no") or d["gate_b"] not in ("yes", "no"):
            raise Stop("datasets.tsv: %s gates are not yes or no" % d["dataset"])
    gaps = [float(d["lineage_gap"]) for d in datasets if d["group"] == "largest_gap"]
    drawn = [float(d["lineage_gap"]) for d in datasets if d["group"] == "drawn"]
    if gaps != sorted(gaps, reverse=True) or (drawn and gaps and max(drawn) > min(gaps)) or any(g < 0.05 for g in drawn):
        raise Stop("datasets.tsv: the gaps do not follow the plan's choice")
    m = index(tables)
    for d in names:
        for split in SPLITS:
            for learner in LEARNERS:
                row = m.get((d, split, learner))
                if row is None:
                    if learner == "catalog":
                        continue
                    raise Stop("measures.tsv: no row for %s %s %s" % (d, split, learner))
                for k in MEASURES:
                    v = row[k]
                    if v is None and not (learner == "baseline" and k != "balanced_accuracy"):
                        raise Stop("measures.tsv: %s %s %s has no %s" % (d, split, learner, k))
                    if v is not None and (v < 0 or (k != "log_loss" and v > 1)):
                        raise Stop("measures.tsv: %s %s %s %s out of range" % (d, split, learner, k))
                if (row["picks"] != "") != (learner in ("best5", "site_best")):
                    raise Stop("measures.tsv: %s %s %s picks" % (d, split, learner))
            best = m[(d, split, "best5")]
            if best["picks"] not in FIVE:
                raise Stop("measures.tsv: %s %s best5 picks %r" % (d, split, best["picks"]))
            if best["balanced_accuracy"] != max(m[(d, split, x)]["balanced_accuracy"] for x in FIVE):
                raise Stop("measures.tsv: %s %s best5 is not the highest balanced accuracy of the five" % (d, split))
            for who in ("best5", "site_best"):
                row = m[(d, split, who)]
                if row["picks"] not in FIVE or any(row[k] != m[(d, split, row["picks"])][k] for k in MEASURES):
                    raise Stop("measures.tsv: %s %s %s differs from the learner it picks" % (d, split, who))
        if m[(d, "lineage", "site_best")]["picks"] != m[(d, "random", "site_best")]["picks"]:
            raise Stop("measures.tsv: %s the site's best differs between splits" % d)
    extra = {k[0] for k in m} - set(names)
    if extra:
        raise Stop("measures.tsv: datasets outside datasets.tsv: %s" % sorted(extra))
    if [r["key"] for r in tables["cost.tsv"]] != COST_KEYS:
        raise Stop("cost.tsv: keys are not %s" % COST_KEYS)
    if [(r["key"], r["value"]) for r in tables["settings.tsv"]] != SETTINGS:
        raise Stop("settings.tsv differs from the script's settings")
    order = tables["order.tsv"]
    if [r["event"] for r in order] != ORDER_EVENTS:
        raise Stop("order.tsv: events are not %s" % ORDER_EVENTS)
    utc = [r["utc"] for r in order]
    if utc != sorted(utc) or not all(re.fullmatch(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ", u) for u in utc):
        raise Stop("order.tsv: the times are not UTC times in order")
    if not re.fullmatch(r"[0-9a-f]{64}", order[0]["sha256"]):
        raise Stop("order.tsv: the plan has no sha256")


# ---------------------------------------------------------------- derived files

def derive(tables):
    datasets = tables["datasets.tsv"]
    species = {d["dataset"]: d["species"] for d in datasets}
    m = index(tables)
    summary, out = [], {}
    for split in SPLITS:
        for who in COMPARATORS:
            for k in MEASURES:
                have = [d["dataset"] for d in datasets if (d["dataset"], split, who) in m
                        and m[(d["dataset"], split, who)][k] is not None]
                if not have:
                    continue
                by = {}
                for d in have:
                    by.setdefault(species[d], []).append(m[(d, split, "candidate")][k] - m[(d, split, who)][k])
                iv = interval(by)
                iv["comparator_median"] = statistics.median(m[(d, split, who)][k] for d in have)
                out[(split, who, k)] = iv
                summary.append({
                    "split": split, "comparator": who, "measure": k,
                    "candidate_median": f6(statistics.median(m[(d, split, "candidate")][k] for d in have)),
                    "comparator_median": f6(iv["comparator_median"]),
                    "difference_median": f6(iv["median"]), "low": f6(iv["lo"]), "high": f6(iv["hi"]),
                    "datasets": iv["n"], "species": iv["species"], "positive": iv["positive"], "negative": iv["negative"]})
    names = [d["dataset"] for d in datasets]
    gain = out[("lineage", "best5", "balanced_accuracy")]
    rec = out[("lineage", "best5", "recall")]
    base_c = statistics.median(m[(d, "lineage", "candidate")]["balanced_accuracy"] - m[(d, "lineage", "baseline")]["balanced_accuracy"] for d in names)
    base_b = statistics.median(m[(d, "lineage", "best5")]["balanced_accuracy"] - m[(d, "lineage", "baseline")]["balanced_accuracy"] for d in names)
    cost = {r["key"]: r["value"] for r in tables["cost.tsv"]}
    hours = float(cost["cpu_seconds"]) / 3600
    rss = int(cost["peak_rss_kb"])
    r1 = gain["median"] >= MARGIN and gain["lo"] > 0
    r2 = base_c >= base_b
    r3 = rec["median"] >= RECALL_FLOOR
    r4 = hours <= BUDGET_CPU_H and rss <= BUDGET_RSS_KB
    rule = [
        {"rule": 1, "statistic": "median network minus comparator balanced accuracy", "value": f6(gain["median"]),
         "needs": "at least %s" % MARGIN, "holds": yes(gain["median"] >= MARGIN)},
        {"rule": 1, "statistic": "lower end of its interval", "value": f6(gain["lo"]), "needs": "above 0",
         "holds": yes(gain["lo"] > 0)},
        {"rule": 2, "statistic": "median network minus baseline balanced accuracy", "value": f6(base_c),
         "needs": "at least %s, the comparator's" % f6(base_b), "holds": yes(r2)},
        {"rule": 3, "statistic": "median network minus comparator recall", "value": f6(rec["median"]),
         "needs": "at least %s" % RECALL_FLOOR, "holds": yes(r3)},
        {"rule": 4, "statistic": "CPU-hours of network fitting", "value": f6(hours), "needs": "at most %s" % BUDGET_CPU_H,
         "holds": yes(hours <= BUDGET_CPU_H)},
        {"rule": 4, "statistic": "peak memory of a network process, kB", "value": rss, "needs": "at most %d" % BUDGET_RSS_KB,
         "holds": yes(rss <= BUDGET_RSS_KB)},
        {"rule": "all", "statistic": "outcome", "value": "adopt" if (r1 and r2 and r3 and r4) else "no gain established",
         "needs": "all four", "holds": yes(r1 and r2 and r3 and r4)},
    ]
    totals = [
        ("datasets", len(datasets)),
        ("species", len(set(species.values()))),
        ("datasets_largest_gap", sum(1 for d in datasets if d["group"] == "largest_gap")),
        ("datasets_drawn", sum(1 for d in datasets if d["group"] == "drawn")),
        ("genomes_smallest", min(int(d["genomes"]) for d in datasets)),
        ("genomes_largest", max(int(d["genomes"]) for d in datasets)),
        ("gap_largest_group_low", f6(min(float(d["lineage_gap"]) for d in datasets if d["group"] == "largest_gap"))),
        ("gap_largest_group_high", f6(max(float(d["lineage_gap"]) for d in datasets if d["group"] == "largest_gap"))),
        ("gate_a", sum(1 for d in datasets if d["gate_a"] == "yes")),
        ("gate_b", sum(1 for d in datasets if d["gate_b"] == "yes")),
        ("fits", cost["fits"]),
        ("catalog_datasets", sum(1 for d in names if (d, "lineage", "catalog") in m)),
        ("site_best_is_best5_lineage", sum(1 for d in names if m[(d, "lineage", "site_best")]["balanced_accuracy"] == m[(d, "lineage", "best5")]["balanced_accuracy"])),
        ("site_best_is_best5_random", sum(1 for d in names if m[(d, "random", "site_best")]["balanced_accuracy"] == m[(d, "random", "best5")]["balanced_accuracy"])),
    ]
    for who in FIVE:
        # How far each learner, fixed in advance, sits below the hindsight best on the lineage split.
        totals.append(("below_best5_lineage_" + who, f6(statistics.median(
            m[(d, "lineage", "best5")]["balanced_accuracy"] - m[(d, "lineage", who)]["balanced_accuracy"] for d in names))))
    return {
        "summary.tsv": tsv("summary.tsv", summary),
        "rule.tsv": tsv("rule.tsv", rule),
        "totals.tsv": tsv("totals.tsv", [{"key": k, "value": v} for k, v in totals]),
    }, out, {"margin_candidate": base_c, "margin_best5": base_b, "cpu_hours": hours, "outcome": rule[-1]["value"]}


README_ROWS = {
    "README.md": "this list",
    "SHA256SUMS": "the sha256 of every other file of this directory",
    "cost.tsv": "the network's fits and processes, its CPU seconds (user plus system) and the peak memory of a process",
    "datasets.tsv": "the 24 datasets in the plan's order: species, group, genomes, resistant genomes, lineage gap and whether the export passed gates A and B",
    "measures.tsv": "one row per dataset, split and learner: balanced accuracy, recall, log loss and specificity over the dataset's genomes; best5 and site_best name the learner they pick",
    "order.tsv": "the plan's file time and sha256, the fitting script's last change, the first network result, the analysis run on the first datasets, and the analysis script's last change",
    "rule.tsv": "derived: the four parts of the plan's rule, the value each takes and whether it holds",
    "settings.tsv": "the network, its inputs and fitting, the datasets, the comparator, the measures, the interval and the rule, as the plan fixed them",
    "summary.tsv": "derived: per split, comparator and measure, the medians over datasets and the median of the network minus the comparator with its interval over species",
    "totals.tsv": "derived: counts the page quotes",
}


def column_lines():
    """Each TSV of the record and its columns, in order."""
    return ["- `%s`: %s." % (name, ", ".join("`%s`" % c for c in COLUMNS[name])) for name in sorted(COLUMNS)]


def readme(files):
    lines = [
        "# Dense network record",
        "",
        "This directory is the record behind [`docs/dense-network.md`](../../docs/dense-network.md): one dense network fitted on 24 public datasets of release `rev14-274-20261112` and compared with the saved calls of GRM's five learners. `scripts/dense-network/record.py` writes every file of it, and its `--check` mode fails on any byte it would not write. No program reads it.",
        "",
        "## Files",
        "",
        "| File | Rows | Contents |",
        "|---|---|---|",
    ]
    for name in sorted(set(files) | {"README.md"}, key=lambda x: (x != "README.md", x != "SHA256SUMS", x)):
        rows = "" if not name.endswith(".tsv") else str(files[name].count("\n") - 1)
        lines.append("| `%s` | %s | %s |" % (name, rows, README_ROWS[name]))
    lines += [
        "",
        "Every file is plain text, and each TSV has a header row. A file holds measures per dataset, split and learner, or medians and counts over datasets: no genome, no sequence, no feature and no fitted parameter.",
        "",
        "## Columns",
        "",
    ] + column_lines() + [
        "",
        "## Rules",
        "",
        "- **Network and comparator.** `candidate` is the network, the mean of its 3 repeats. `best5` is the best of the five in hindsight on that split, the plan's comparator on the lineage split.",
        "- **Other comparators.** `site_best` is the learner the release names best, chosen on the random split; `mean4` averages the scores of `cart`, `rf`, `gbt` and `lasso`; `baseline` calls by lineage identity and has balanced accuracy only; `catalog` is the catalog rule, where the release scores it.",
        "- **Over datasets.** `summary.tsv` gives the median over datasets of the network minus the comparator, with the 2.5% and 97.5% percentiles of that median over 2,000 resamples of species at seed 20261115, drawn as the experiment drew them.",
        "- **Draw.** Each resample takes every dataset of each drawn species.",
        "- **Log loss.** It is clipped at 0.0001; a positive log loss difference means the network's is higher.",
        "- **Specificity.** The experiment's report gives no specificity figure over datasets. `summary.tsv` works them out from the rows of `measures.tsv` by the same draw.",
        "- **Readings.** The plan fixed the rule on the lineage split against `best5`. Every other comparator, the random split, recall, log loss and specificity are reported beside the rule and do not enter it.",
        "- **Exploratory.** `below_best5_lineage_*` in `totals.tsv` and the counts of `site_best` were read by the review after the results, so they are exploratory.",
        "",
        "## Commands",
        "",
        "```sh",
        "python3 -I scripts/dense-network/record.py --check reference/dense-network",
        "python3 -I scripts/dense-network/record.py --source \"$RUN\" --review \"$REVIEW\" --check reference/dense-network",
        "```",
        "",
        "The first rebuilds the derived files, this README and `SHA256SUMS` from the row files. The second rebuilds every file from the experiment's and the review's output directories, which are not public, and fails if a derived figure of balanced accuracy, recall or log loss parts from the experiment's own by more than 1e-12.",
        "",
    ]
    return "\n".join(lines)


def assemble(row_texts):
    tables = {}
    for name in ROW_FILES:
        head, rows = parse(row_texts[name])
        if head != COLUMNS[name]:
            raise Stop("%s: columns %s" % (name, head))
        tables[name] = rows
    validate(tables)
    derived, stats, extra = derive(tables)
    files = dict(row_texts)
    files.update(derived)
    files["SHA256SUMS"] = ""
    files["README.md"] = readme(files)
    del files["SHA256SUMS"]
    files["SHA256SUMS"] = "".join("%s  %s\n" % (hashlib.sha256(files[k].encode("utf-8")).hexdigest(), k) for k in sorted(files))
    return files, stats, extra


# ---------------------------------------------------------------- the sources

def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def one(pattern, text, where):
    hits = re.findall(pattern, text)
    if len(hits) != 1:
        raise Stop("%s: %d matches for %r" % (where, len(hits), pattern))
    return hits[0]


def from_sources(run, review):
    a, r = pathlib.Path(run), pathlib.Path(review)
    plan_bytes = (a / "PLAN.md").read_bytes()
    plan = plan_bytes.decode("utf-8").replace("`", "")
    report = read(a / "REPORT.md")
    rev = read(r / "REPORT.md")
    digest = hashlib.sha256(plan_bytes).hexdigest()
    stated = one(r"PLAN\.md` closed all six before any candidate result; sha256 `([0-9a-f]{64})`", report, "REPORT.md")
    if stated != digest:
        raise Stop("PLAN.md is not the plan whose sha256 the report gives")
    for key, wanted in SETTINGS_IN_PLAN.items():
        for text in (wanted,) if isinstance(wanted, str) else wanted:
            if text not in plan:
                raise Stop("PLAN.md does not give %s as %r" % (key, text))
    plan_time, first = one(r"file time (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ); the first candidate result was written at (\d\d:\d\d:\d\dZ)",
                           report, "REPORT.md")
    day = plan_time[:11]
    fitting = one(r"`ba2_candidate\.py` last change (\d\d:\d\d:\d\dZ)", rev, "review REPORT.md")
    look = one(r"The \w+ ran its analysis on the outcomes of the first \d+ datasets at (\d\d:\d\d:\d\dZ)", rev, "review REPORT.md")
    last = one(r"The analysis script was last changed at (\d\d:\d\d:\d\dZ)", rev, "review REPORT.md")
    order = [{"utc": plan_time, "event": "plan_file_time", "sha256": digest},
             {"utc": day + fitting, "event": "fitting_script_last_change", "sha256": ""},
             {"utc": day + first, "event": "first_network_result", "sha256": ""},
             {"utc": day + look, "event": "analysis_on_first_datasets", "sha256": ""},
             {"utc": day + last, "event": "analysis_script_last_change", "sha256": ""}]

    summary = json.loads(read(a / "results" / "ba2-summary.json"))
    stage1 = json.loads(read(a / "results" / "stage1-datasets.json"))
    gates = {g["id"]: g for g in json.loads(read(a / "results" / "gates.json"))}
    if stage1["seed"] != SEED or stage1["pool"] != 113:
        raise Stop("stage1-datasets.json: seed or pool differs from the plan")
    if summary["excluded"]:
        raise Stop("ba2-summary.json: datasets left out: %s" % summary["excluded"])
    rows = {x["id"]: x for x in summary["rows"]}
    chosen = [(d, "largest_gap") for d in stage1["top12"]] + [(d, "drawn") for d in stage1["drawn12"]]
    if [d for d, _ in chosen] != [x["id"] for x in summary["rows"]]:
        raise Stop("ba2-summary.json: rows are not the plan's datasets in order")
    datasets, measures = [], []
    for i, (d, group) in enumerate(chosen, 1):
        x, g = rows[d], gates[d]
        if bool(g["gate_b"]) != bool(x["gate_b"]):
            raise Stop("gates.json and ba2-summary.json disagree on gate B for %s" % d)
        datasets.append({"order": i, "dataset": d, "species": x["species"], "group": group, "genomes": x["n"],
                         "resistant": x["n_resistant"], "lineage_gap": cell(stage1["gap"][d]),
                         "gate_a": yes(g["gate_a"]), "gate_b": yes(g["gate_b"])})
        for split in SPLITS:
            for learner in LEARNERS:
                v = x[split].get(learner)
                if v is None:
                    continue
                keep = MEASURES if learner != "baseline" else ["balanced_accuracy"]
                row = {"dataset": d, "split": split, "learner": learner, "picks": v.get("learner", "")}
                for k in MEASURES:
                    row[k] = cell(v[SOURCE_KEY[k]]) if k in keep else ""
                measures.append(row)
    c = summary["cost"]
    cost = [{"key": "fits", "value": c["fits"]}, {"key": "processes", "value": c["processes"]},
            {"key": "cpu_seconds", "value": cell(c["cpu_s"])}, {"key": "peak_rss_kb", "value": c["peak_rss_kb"]}]
    texts = {
        "datasets.tsv": tsv("datasets.tsv", datasets),
        "measures.tsv": tsv("measures.tsv", measures),
        "cost.tsv": tsv("cost.tsv", cost),
        "order.tsv": tsv("order.tsv", order),
        "settings.tsv": tsv("settings.tsv", [{"key": k, "value": v} for k, v in SETTINGS]),
    }
    return texts, summary


def against_experiment(stats, extra, summary):
    """Fail unless every derived figure equals the experiment's own within CLOSE."""
    for (split, who, k), iv in stats.items():
        if k not in REPORTED:
            continue
        theirs = summary["paired"][split][who]
        theirs = theirs.get(SOURCE_KEY[k]) if theirs else None
        if theirs is None:
            raise Stop("ba2-summary.json has no %s %s %s" % (split, who, k))
        for f in ("median", "lo", "hi"):
            if abs(iv[f] - theirs[f]) > CLOSE:
                raise Stop("%s %s %s %s: %r against the experiment's %r" % (split, who, k, f, iv[f], theirs[f]))
        for f in ("n", "species", "positive", "negative"):
            if iv[f] != theirs[f]:
                raise Stop("%s %s %s %s: %r against the experiment's %r" % (split, who, k, f, iv[f], theirs[f]))
        med = summary["medians"][split][who].get(SOURCE_KEY[k])
        if med is None or abs(iv["comparator_median"] - med) > CLOSE:
            raise Stop("%s %s %s median: %r against the experiment's %r" % (split, who, k, iv["comparator_median"], med))
    p = summary["primary"]
    for f in ("margin_candidate", "margin_best5", "cpu_hours"):
        if abs(extra[f] - p[f]) > CLOSE:
            raise Stop("%s: %r against the experiment's %r" % (f, extra[f], p[f]))
    if extra["outcome"] != summary["outcome"]:
        raise Stop("outcome: %r against the experiment's %r" % (extra["outcome"], summary["outcome"]))


def from_record(directory):
    d = pathlib.Path(directory)
    texts = {}
    for name in ROW_FILES:
        if not (d / name).is_file():
            raise Stop("%s is missing" % name)
        texts[name] = read(d / name)
    return texts


# ---------------------------------------------------------------- main

def check(files, directory):
    d = pathlib.Path(directory)
    on_disk = {p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file()}
    bad = sorted(on_disk ^ set(files))
    for name in sorted(files):
        if (d / name).is_file() and (d / name).read_bytes() != files[name].encode("utf-8"):
            bad.append(name)
    return sorted(set(bad))


def earlier_record(out, files):
    """Whether a directory holds nothing, or an earlier record of this script and nothing else."""
    found = list(out.iterdir())
    if not found:
        return True
    if any(not p.is_file() for p in found) or {p.name for p in found} - set(files):
        return False
    sums, readme = out / "SHA256SUMS", out / "README.md"
    if not (sums.is_file() and readme.is_file()):
        return False
    listed = {ln.split("  ", 1)[-1] for ln in sums.read_text(encoding="utf-8", errors="replace").splitlines()}
    first = readme.read_text(encoding="utf-8", errors="replace").split("\n", 1)[0]
    return "README.md" in listed and first == files["README.md"].split("\n", 1)[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--source", help="the experiment's output directory")
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
        if args.source:
            texts, summary = from_sources(args.source, args.review)
            files, stats, extra = assemble(texts)
            against_experiment(stats, extra, summary)
        else:
            files, _, _ = assemble(from_record(args.check))
    except (Stop, OSError, KeyError, ValueError, IndexError, TypeError) as exc:
        print("record.py: stopped: %s" % exc, file=sys.stderr)
        return 1
    if args.out:
        out = pathlib.Path(args.out)
        if out.exists() and not (out.is_dir() and earlier_record(out, files)):
            print("record.py: stopped: %s holds files that are not an earlier record of this script" % out,
                  file=sys.stderr)
            return 1
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
    print("record.py: %s is the record" % args.check)
    return 0


if __name__ == "__main__":
    sys.exit(main())
