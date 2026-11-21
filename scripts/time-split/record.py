#!/usr/bin/env python3
"""Build the record of the split by first-public day, rerun after the review.

    python3 -I scripts/time-split/record.py --inputs DIR --out reference/time-split
    python3 -I scripts/time-split/record.py --inputs DIR --check reference/time-split
    python3 -I scripts/time-split/record.py --check reference/time-split

With --inputs the script reads every result file of the input set, computes
each measure again from the scores and labels, holds it to the evaluator's own
counts, and writes the row files; it then derives every other file from the row
files. --out writes the record, --check rebuilds it and fails on any byte that
differs. --check alone reads the record's row files, derives the other files
again and fails on any byte that differs from them or from SHA256SUMS.

The input set DIR holds:

    design.json                                  the design: datasets, pools, fold file sha256
    results/<dataset>.first-public-<arm>.json    the 72 result files of the rerun after the review
    before/<dataset>.first-public-<arm>.json     the 72 result files of 14 Nov, read only to count what the fix changed
    studies.json                                 distinct genomes and studies behind the test sets, corrected
    plan/PLAN.md, plan/PLAN.sha256               the plan of 14 Nov and its hash list
    plan/tables.sha256                           the hash of the table layout script
    plan/scripts/analyze.py, plan/scripts/tables.py   the two scripts the hash lists name
    rerun/RERUN-PLAN.md, rerun/RERUN-PLAN.sha256 the rerun's plan and its hash
    logs/first-run.log, logs/rerun.log           the run logs of 14 Nov and of the rerun

The record holds no genome id, label or score per genome, no study accession,
no fold file and no list of test-set studies. Its rows are per dataset, per
species, per arm, per learner, per event or a single count. The script reads
nothing else, fits nothing and opens no network connection.
"""
import argparse
import hashlib
import itertools
import json
import math
import random
import re
import statistics
import sys
from fractions import Fraction
from pathlib import Path

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
ARMS = ("past", "random", "study-out")
CONTRASTS = (("total", "past", "random"), ("study", "study-out", "random"), ("time", "past", "study-out"))
# In this order, as the plan's analysis script draws its intervals.
MEASURES = ("recall", "log_loss", "specificity", "balanced_accuracy")
WORSE = {"recall": -1, "log_loss": 1, "specificity": -1, "balanced_accuracy": -1}
RULE_MEASURES = ("recall", "log_loss")
CUTOFF = 0.5
CLIP = 1e-4
RESAMPLES = 10000
RESAMPLE_SEED = 20261114
LOW_DRAW, HIGH_DRAW = 251, 9750  # the 251st and the 9,750th of the sorted draws
MIN_DATASETS, MIN_SPECIES, MIN_LEARNERS_WITH_COST = 12, 3, 3
MAX_EXACT_SPECIES = 6
RESULT_FORMAT = "grm-fold-file-result-v1"
REVISION = 14

ROW_FILES = ("fits.tsv", "datasets.tsv", "species.tsv", "review.tsv", "totals.tsv", "order.tsv", "settings.tsv",
             "checks.tsv", "inputs.sha256")
DERIVED_FILES = ("contrasts.tsv", "readings.tsv", "species-medians.tsv", "rule.tsv", "arms.tsv", "mix.tsv",
                 "summary.tsv")
FILES = ROW_FILES + DERIVED_FILES + ("README.md", "SHA256SUMS")

COLUMNS = {
    "fits.tsv": ["dataset", "species", "arm", "learner", "test", "test_resistant", "train", "train_resistant",
                 "tp", "fn", "tn", "fp", "wrong_at_clip", "scores_0_or_1", "loss_sum", "clip_loss_sum", "null_fit"],
    "datasets.tsv": ["dataset", "species", "antibiotic", "cohort", "cohort_resistant", "cut_day",
                     "share_public_by_cut", "drawn_test_studies", "studies_with_a_tested_genome", "test",
                     "test_resistant", "train", "train_resistant"]
    + [f"{a}_from_{s}" for a in ("past", "random", "study_out") for s in ("past", "later_other_studies",
                                                                          "later_test_studies")],
    "species.tsv": ["species", "datasets", "tested_pairs", "tested_resistant_pairs", "distinct_test_genomes",
                    "distinct_test_studies", "largest_study_tested_genomes"],
    "review.tsv": ["learner", "fits", "fits_changed", "scores", "scores_changed", "readings", "readings_changed",
                   "changed_readings"],
    "totals.tsv": ["name", "value"],
    "order.tsv": ["event", "time", "source"],
    "settings.tsv": ["name", "value"],
    "checks.tsv": ["check", "count"],
    "contrasts.tsv": ["dataset", "species", "learner", "contrast"] + ["d_" + m for m in MEASURES],
    "readings.tsv": ["learner", "contrast", "minuend", "subtrahend", "measure", "median", "low", "high", "reading",
                     "exact_low", "exact_high", "exact_reading", "datasets", "worse", "better", "equal",
                     "minuend_median", "subtrahend_median", "plan"],
    "species-medians.tsv": ["learner", "contrast", "measure", "species", "datasets", "median"],
    "rule.tsv": ["learner", "recall", "log_loss", "rule"],
    "arms.tsv": ["learner", "arm", "datasets", "tested_pairs", "resistant_pairs", "tp", "fn", "tn", "fp", "recall",
                 "specificity", "wrong_calls", "wrong_at_clip", "scores_0_or_1", "pooled_log_loss", "clip_share",
                 "null_fits"],
    "mix.tsv": ["arm", "datasets", "train", "from_past", "from_later_other_studies", "from_later_test_studies",
                "past_share"],
    "summary.tsv": ["name", "value"],
}
INTS = {"test", "test_resistant", "train", "train_resistant", "tp", "fn", "tn", "fp", "wrong_at_clip",
        "scores_0_or_1", "null_fit", "cohort", "cohort_resistant", "drawn_test_studies",
        "studies_with_a_tested_genome", "datasets", "tested_pairs", "tested_resistant_pairs",
        "distinct_test_genomes", "distinct_test_studies", "largest_study_tested_genomes", "fits", "fits_changed",
        "scores", "scores_changed", "readings", "readings_changed"} | {c for c in COLUMNS["datasets.tsv"] if "_from_" in c}
FLOATS = {"loss_sum", "clip_loss_sum", "share_public_by_cut"}
# Settings the script applies; settings.tsv must give each one as here.
FIXED_SETTINGS = {"evaluation_revision": str(REVISION), "cutoff": "0.5", "clip": "0.0001",
                  "loss_of_a_wrong_call_at_the_clip": repr(-math.log(CLIP)),
                  "interval_level": "0.95", "resamples": str(RESAMPLES), "resample_seed": str(RESAMPLE_SEED),
                  "low_draw": str(LOW_DRAW), "high_draw": str(HIGH_DRAW), "min_datasets": str(MIN_DATASETS),
                  "min_species": str(MIN_SPECIES), "min_learners_with_cost": str(MIN_LEARNERS_WITH_COST)}
# Shapes of what must never enter the record: genome ids, accessions of studies, samples and assemblies, sequence.
BANNED_SHAPES = (r"\b\d{3,7}\.\d{3,6}\b", r"\bPRJ[A-Z]{1,2}\d+", r"\bSAM[NED][A-Z]?\d+", r"\bGC[AF]_\d+",
                 r"\b[SE]R[RSX]\d{5,}", r"[ACGT]{20,}")


class Stop(Exception):
    pass


def need(cond, msg):
    if not cond:
        raise Stop(msg)


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def fmt(v):
    """A cell: ints as written, floats as Python's shortest round-trip form."""
    if isinstance(v, bool):
        return str(int(v))
    if isinstance(v, float):
        need(math.isfinite(v), "a cell is not a finite number")
        return repr(v)
    return str(v)


def tsv(name, rows):
    cols = COLUMNS[name]
    out = ["\t".join(cols)]
    for r in rows:
        need(set(r) == set(cols), f"{name}: a row has other cells than its columns")
        cells = [fmt(r[c]) for c in cols]
        for c in cells:
            need("\t" not in c and "\n" not in c, f"{name}: a cell holds a tab or a line break")
        out.append("\t".join(cells))
    return "\n".join(out) + "\n"


def parse_tsv(name, text):
    """A row file of the record as dicts, held to its columns; ints and floats parsed."""
    lines = text.split("\n")
    need(len(lines) >= 2 and lines[-1] == "", f"{name}: no final newline")
    head = lines[0].split("\t")
    need(head == COLUMNS[name], f"{name}: the columns are not {COLUMNS[name]}")
    rows = []
    for i, line in enumerate(lines[1:-1], 2):
        cells = line.split("\t")
        need(len(cells) == len(head), f"{name}:{i}: {len(cells)} cells, the header has {len(head)}")
        row = {}
        for c, v in zip(head, cells):
            if c in INTS:
                need(re.fullmatch(r"0|[1-9]\d*", v) is not None, f"{name}:{i}: {c} is not a count")
                row[c] = int(v)
            elif c in FLOATS:
                try:
                    row[c] = float(v)
                except ValueError:
                    raise Stop(f"{name}:{i}: {c} is not a number")
                need(fmt(row[c]) == v, f"{name}:{i}: {c} is not in its shortest form")
            else:
                row[c] = v
        rows.append(row)
    return rows


# ----------------------------------------------------------------- measures

def fit_row(dataset, species, arm, learner, truth, scores, train, train_resistant, null_fit):
    """One fit's counts and loss sums, in sample order, as the plan's analysis script sums them."""
    need(len(truth) == len(scores), f"{dataset} {arm} {learner}: {len(scores)} scores for {len(truth)} genomes")
    tp = fn = tn = fp = at_clip = zero_one = 0
    loss = clip_loss = 0.0
    for t, s in zip(truth, scores):
        need(isinstance(s, (int, float)) and not isinstance(s, bool) and 0 <= s <= 1,
             f"{dataset} {arm} {learner}: a score outside 0 to 1")
        called = s >= CUTOFF
        tp += t and called
        fn += t and not called
        tn += (not t) and not called
        fp += (not t) and called
        p = min(max(s, CLIP), 1 - CLIP)
        one = -math.log(p if t else 1 - p)
        loss -= math.log(p if t else 1 - p)
        if (t and s <= CLIP) or (not t and s >= 1 - CLIP):
            at_clip += 1
            clip_loss += one
        zero_one += s in (0.0, 1.0)
    return {"dataset": dataset, "species": species, "arm": arm, "learner": learner, "test": len(truth),
            "test_resistant": sum(truth), "train": train, "train_resistant": train_resistant, "tp": tp, "fn": fn,
            "tn": tn, "fp": fp, "wrong_at_clip": at_clip, "scores_0_or_1": zero_one, "loss_sum": loss,
            "clip_loss_sum": clip_loss, "null_fit": int(null_fit)}


def measures(row):
    recall = row["tp"] / (row["tp"] + row["fn"])
    specificity = row["tn"] / (row["tn"] + row["fp"])
    return {"recall": recall, "specificity": specificity, "balanced_accuracy": (recall + specificity) / 2,
            "log_loss": row["loss_sum"] / row["test"]}


def resampled(by_species, rng):
    """Median over datasets, and the 251st and 9,750th of 10,000 medians over species drawn with replacement."""
    species = sorted(by_species)
    point = statistics.median([v for s in species for v in by_species[s]])
    draws = []
    for _ in range(RESAMPLES):
        values = [v for s in rng.choices(species, k=len(species)) for v in by_species[s]]
        draws.append(statistics.median(values))
    draws.sort()
    return point, draws[LOW_DRAW - 1], draws[HIGH_DRAW - 1]


def exact(by_species):
    """The same two quantiles over every ordered draw of the species, each drawn as often as any other."""
    species = sorted(by_species)
    need(len(species) <= MAX_EXACT_SPECIES, f"{len(species)} species: too many ordered draws to take them all")
    meds = sorted(statistics.median([v for s in draw for v in by_species[s]])
                  for draw in itertools.product(species, repeat=len(species)))
    n = len(meds)
    return meds[math.ceil(Fraction(LOW_DRAW, RESAMPLES) * n) - 1], meds[math.ceil(Fraction(HIGH_DRAW, RESAMPLES) * n) - 1]


def read(lo, hi, measure):
    if lo > 0 or hi < 0:
        return "worse" if (lo > 0) == (WORSE[measure] > 0) else "better"
    return "not shown"


def readings(fits):
    """The 60 readings, drawn in the order of the plan's analysis script with one generator."""
    datasets = list(dict.fromkeys(r["dataset"] for r in fits))
    species_of = {r["dataset"]: r["species"] for r in fits}
    value = {(r["dataset"], r["learner"], r["arm"]): measures(r) for r in fits}
    rng = random.Random(RESAMPLE_SEED)
    out, contrast_rows, species_rows = [], [], []
    for d in datasets:
        for learner in LEARNERS:
            for name, a, b in CONTRASTS:
                row = {"dataset": d, "species": species_of[d], "learner": learner, "contrast": name}
                for m in MEASURES:
                    row["d_" + m] = value[(d, learner, a)][m] - value[(d, learner, b)][m]
                contrast_rows.append(row)
    for learner in LEARNERS:
        for name, a, b in CONTRASTS:
            for m in MEASURES:
                by_species = {}
                for d in datasets:
                    by_species.setdefault(species_of[d], []).append(
                        value[(d, learner, a)][m] - value[(d, learner, b)][m])
                values = [v for s in by_species for v in by_species[s]]
                point, lo, hi = resampled(by_species, rng)
                elo, ehi = exact(by_species)
                w = WORSE[m]
                out.append({"learner": learner, "contrast": name, "minuend": a, "subtrahend": b, "measure": m,
                            "median": point, "low": lo, "high": hi, "reading": read(lo, hi, m),
                            "exact_low": elo, "exact_high": ehi, "exact_reading": read(elo, ehi, m),
                            "datasets": len(values), "worse": sum(1 for v in values if v * w > 0),
                            "better": sum(1 for v in values if v * w < 0),
                            "equal": sum(1 for v in values if v == 0),
                            "minuend_median": statistics.median([value[(d, learner, a)][m] for d in datasets]),
                            "subtrahend_median": statistics.median([value[(d, learner, b)][m] for d in datasets]),
                            "plan": "rule" if name == "total" and m in RULE_MEASURES
                            else "fixed" if m in RULE_MEASURES else "supporting"})
                for s in sorted(by_species):
                    species_rows.append({"learner": learner, "contrast": name, "measure": m, "species": s,
                                         "datasets": len(by_species[s]), "median": statistics.median(by_species[s])})
    return out, contrast_rows, species_rows


def rule(reading_rows, n_datasets, n_species):
    per = []
    for learner in LEARNERS:
        read_ = {r["measure"]: r["reading"] for r in reading_rows
                 if r["learner"] == learner and r["contrast"] == "total" and r["measure"] in RULE_MEASURES}
        if "better" in read_.values():
            verdict = "benefit" if "worse" not in read_.values() else "mixed"
        elif "worse" in read_.values():
            verdict = "cost"
        else:
            verdict = "not shown"
        per.append({"learner": learner, "recall": read_["recall"], "log_loss": read_["log_loss"], "rule": verdict})
    cost = sum(r["rule"] == "cost" for r in per)
    odd = sum(r["rule"] in ("benefit", "mixed") for r in per)
    if n_datasets < MIN_DATASETS or n_species < MIN_SPECIES:
        answer = "descriptive only"
    elif odd:
        answer = "mixed"
    elif cost >= MIN_LEARNERS_WITH_COST:
        answer = "cost shown"
    else:
        answer = "cost not shown"
    return per, answer, cost


# ----------------------------------------------------------------- derive

def check_rows(rows):
    """What the row files must give of each other before anything is derived from them."""
    fits, datasets, species = rows["fits.tsv"], rows["datasets.tsv"], rows["species.tsv"]
    ids = [d["dataset"] for d in datasets]
    need(len(ids) == len(set(ids)) and ids, "datasets.tsv: a dataset stands twice, or none stands")
    need([(r["dataset"], r["arm"], r["learner"]) for r in fits]
         == [(d, a, l) for d in ids for a in ARMS for l in LEARNERS],
         "fits.tsv: not one row per dataset, arm and learner, in the order of datasets.tsv")
    by_id = {d["dataset"]: d for d in datasets}
    for d in datasets:
        need(re.fullmatch(re.escape(d["species"]) + r"__[a-z0-9_]+", d["dataset"]) is not None,
             f"{d['dataset']}: not species__antibiotic")
        need(0 < d["test_resistant"] < d["test"] and 0 < d["train_resistant"] < d["train"], f"{d['dataset']}: a class is empty")
        need(d["studies_with_a_tested_genome"] <= d["drawn_test_studies"], f"{d['dataset']}: a tested study not drawn")
        for a in ("past", "random", "study_out"):
            need(sum(d[f"{a}_from_{s}"] for s in ("past", "later_other_studies", "later_test_studies")) == d["train"],
                 f"{d['dataset']}: the {a} arm's sources do not add up to its training genomes")
        need(d["past_from_later_other_studies"] == d["past_from_later_test_studies"] == 0,
             f"{d['dataset']}: the past arm trains on a later genome")
        need(d["study_out_from_later_test_studies"] == 0, f"{d['dataset']}: the study-out arm trains on a test study")
    for r in fits:
        d = by_id[r["dataset"]]
        need(r["species"] == d["species"], f"fits.tsv {r['dataset']}: another species")
        need((r["test"], r["test_resistant"], r["train"], r["train_resistant"])
             == (d["test"], d["test_resistant"], d["train"], d["train_resistant"]),
             f"fits.tsv {r['dataset']} {r['arm']} {r['learner']}: counts other than datasets.tsv's")
        need(r["tp"] + r["fn"] == r["test_resistant"] and r["tn"] + r["fp"] == r["test"] - r["test_resistant"],
             f"fits.tsv {r['dataset']} {r['arm']} {r['learner']}: the confusion counts do not add up")
        need(r["wrong_at_clip"] <= r["fn"] + r["fp"] and r["scores_0_or_1"] <= r["test"],
             f"fits.tsv {r['dataset']} {r['arm']} {r['learner']}: more calls at the clip than wrong calls")
        need(0 <= r["clip_loss_sum"] <= r["loss_sum"], f"fits.tsv {r['dataset']} {r['arm']} {r['learner']}: loss sums")
        need(r["null_fit"] in (0, 1) and (r["null_fit"] == 0 or r["learner"] in ("gbt", "lasso")),
             f"fits.tsv {r['dataset']} {r['arm']} {r['learner']}: a null fit flag")
    names = [s["species"] for s in species]
    need(names == sorted(set(d["species"] for d in datasets)), "species.tsv: not one row per species of datasets.tsv, sorted")
    for s in species:
        mine = [d for d in datasets if d["species"] == s["species"]]
        need(s["datasets"] == len(mine), f"species.tsv {s['species']}: datasets")
        need(s["tested_pairs"] == sum(d["test"] for d in mine), f"species.tsv {s['species']}: tested pairs")
        need(s["tested_resistant_pairs"] == sum(d["test_resistant"] for d in mine), f"species.tsv {s['species']}: resistant pairs")
        need(s["largest_study_tested_genomes"] <= s["distinct_test_genomes"] <= s["tested_pairs"],
             f"species.tsv {s['species']}: genome counts")
        need(max(d["studies_with_a_tested_genome"] for d in mine) <= s["distinct_test_studies"],
             f"species.tsv {s['species']}: fewer studies than a dataset of it")
    need([r["learner"] for r in rows["review.tsv"]] == list(LEARNERS), "review.tsv: not one row per learner")
    for r in rows["review.tsv"]:
        need(r["fits"] == len(ids) * len(ARMS) and r["fits_changed"] <= r["fits"], f"review.tsv {r['learner']}: fits")
        need(r["scores"] == len(ARMS) * sum(d["test"] for d in datasets) and r["scores_changed"] <= r["scores"],
             f"review.tsv {r['learner']}: scores")
        need((r["fits_changed"] == 0) == (r["scores_changed"] == 0), f"review.tsv {r['learner']}: a changed fit with no changed score")
        need(r["readings"] == len(CONTRASTS) * len(MEASURES) and r["readings_changed"] <= r["readings"],
             f"review.tsv {r['learner']}: readings")
        named = [x for x in r["changed_readings"].split(";") if x]
        need(len(named) == r["readings_changed"], f"review.tsv {r['learner']}: the changed readings named are not counted")
        for x in named:
            need(re.fullmatch(r"(total|study|time) (" + "|".join(MEASURES) + ")", x) is not None,
                 f"review.tsv {r['learner']}: {x!r} is not a contrast and a measure")
    settings = {r["name"]: r["value"] for r in rows["settings.tsv"]}
    for k, v in FIXED_SETTINGS.items():
        need(settings.get(k) == v, f"settings.tsv: {k} is not {v}")
    totals = {r["name"]: r["value"] for r in rows["totals.tsv"]}
    for k in ("distinct_test_genomes", "distinct_test_studies", "sum_of_per_species_study_counts",
              "studies_tested_in_two_species", "result_files", "result_files_changed", "datasets_with_a_changed_fit"):
        need(re.fullmatch(r"0|[1-9]\d*", totals.get(k, "")) is not None, f"totals.tsv: no count {k}")
    need(int(totals["distinct_test_genomes"]) == sum(s["distinct_test_genomes"] for s in species),
         "totals.tsv: distinct test genomes are not the sum over species")
    need(int(totals["sum_of_per_species_study_counts"]) == sum(s["distinct_test_studies"] for s in species),
         "totals.tsv: the per-species study counts do not add up")
    need(int(totals["distinct_test_studies"]) + int(totals["studies_tested_in_two_species"])
         == int(totals["sum_of_per_species_study_counts"]), "totals.tsv: studies in two species do not explain the sum")
    need(int(totals["result_files"]) == len(ids) * len(ARMS), "totals.tsv: result files")
    need(int(totals["datasets_with_a_changed_fit"]) <= len(ids) and int(totals["result_files_changed"])
         <= int(totals["result_files"]), "totals.tsv: changed files or datasets")


def derive(rows):
    """Every derived file, from the row files alone."""
    check_rows(rows)
    fits, datasets = rows["fits.tsv"], rows["datasets.tsv"]
    n_datasets = len(datasets)
    n_species = len({d["species"] for d in datasets})
    reading_rows, contrast_rows, species_rows = readings(fits)
    per, answer, cost = rule(reading_rows, n_datasets, n_species)
    arm_rows = []
    for learner in LEARNERS:
        for arm in ARMS:
            rs = [r for r in fits if r["learner"] == learner and r["arm"] == arm]
            tot = {k: sum(r[k] for r in rs) for k in ("test", "test_resistant", "tp", "fn", "tn", "fp", "wrong_at_clip",
                                                      "scores_0_or_1", "null_fit")}
            loss = sum(r["loss_sum"] for r in rs)
            clip = sum(r["clip_loss_sum"] for r in rs)
            arm_rows.append({"learner": learner, "arm": arm, "datasets": len(rs), "tested_pairs": tot["test"],
                             "resistant_pairs": tot["test_resistant"], "tp": tot["tp"], "fn": tot["fn"],
                             "tn": tot["tn"], "fp": tot["fp"], "recall": tot["tp"] / (tot["tp"] + tot["fn"]),
                             "specificity": tot["tn"] / (tot["tn"] + tot["fp"]), "wrong_calls": tot["fn"] + tot["fp"],
                             "wrong_at_clip": tot["wrong_at_clip"], "scores_0_or_1": tot["scores_0_or_1"],
                             "pooled_log_loss": loss / tot["test"], "clip_share": clip / loss,
                             "null_fits": tot["null_fit"]})
    mix_rows = []
    for arm, key in (("past", "past"), ("random", "random"), ("study-out", "study_out")):
        train = sum(d["train"] for d in datasets)
        src = {s: sum(d[f"{key}_from_{s}"] for d in datasets) for s in ("past", "later_other_studies", "later_test_studies")}
        mix_rows.append({"arm": arm, "datasets": n_datasets, "train": train, "from_past": src["past"],
                         "from_later_other_studies": src["later_other_studies"],
                         "from_later_test_studies": src["later_test_studies"], "past_share": src["past"] / train})
    tr = [d["train"] for d in datasets]
    te = [d["test"] for d in datasets]
    st = [d["studies_with_a_tested_genome"] for d in datasets]
    dr = [d["drawn_test_studies"] for d in datasets]
    total = [r for r in reading_rows if r["contrast"] == "total"]

    def learners_reading(rows_, measure, word):
        return sum(1 for r in rows_ if r["measure"] == measure and r["reading"] == word)

    summary = [
        ("datasets", n_datasets), ("species", n_species), ("tested_pairs", sum(te)),
        ("resistant_pairs", sum(d["test_resistant"] for d in datasets)),
        ("susceptible_pairs", sum(te) - sum(d["test_resistant"] for d in datasets)),
        ("test_min", min(te)), ("test_max", max(te)), ("test_median", statistics.median(te)),
        ("train_min", min(tr)), ("train_max", max(tr)), ("train_median", statistics.median(tr)),
        ("studies_with_a_tested_genome_min", min(st)), ("studies_with_a_tested_genome_max", max(st)),
        ("datasets_with_one_tested_study", sum(1 for x in st if x == 1)),
        ("drawn_test_studies_min", min(dr)), ("drawn_test_studies_max", max(dr)),
        ("datasets_with_a_drawn_study_untested", sum(1 for d in datasets
                                                     if d["studies_with_a_tested_genome"] < d["drawn_test_studies"])),
        ("fits", len(fits)), ("readings", len(reading_rows)),
        ("answer", answer), ("learners_with_cost", cost), ("learners", len(LEARNERS)),
        ("total_recall_worse", learners_reading(total, "recall", "worse")),
        ("total_log_loss_worse", learners_reading(total, "log_loss", "worse")),
        ("total_specificity_worse", learners_reading(total, "specificity", "worse")),
        ("total_balanced_accuracy_worse", learners_reading(total, "balanced_accuracy", "worse")),
        ("readings_better", sum(1 for r in reading_rows if r["reading"] == "better")),
        ("readings_worse", sum(1 for r in reading_rows if r["reading"] == "worse")),
        ("exact_readings_agree", sum(1 for r in reading_rows if r["exact_reading"] == r["reading"])),
        ("ordered_species_draws", n_species ** n_species),
        ("distinct_species_draws", math.comb(2 * n_species - 1, n_species)),
        ("study_log_loss_median_above_time", sum(
            1 for l in LEARNERS
            if next(r["median"] for r in reading_rows if (r["learner"], r["contrast"], r["measure"]) == (l, "study", "log_loss"))
            > next(r["median"] for r in reading_rows if (r["learner"], r["contrast"], r["measure"]) == (l, "time", "log_loss")))),
        ("review_fits_changed", sum(r["fits_changed"] for r in rows["review.tsv"])),
        ("review_fits", sum(r["fits"] for r in rows["review.tsv"])),
        ("review_scores_changed", sum(r["scores_changed"] for r in rows["review.tsv"])),
        ("review_scores", sum(r["scores"] for r in rows["review.tsv"])),
        ("review_readings_changed", sum(r["readings_changed"] for r in rows["review.tsv"])),
        ("review_readings", sum(r["readings"] for r in rows["review.tsv"])),
        ("review_total_contrast_readings_changed", sum(
            1 for r in rows["review.tsv"] for x in r["changed_readings"].split(";") if x.startswith("total "))),
        ("record_files", len(FILES)),
    ]
    by_name = {}
    for name, value in summary:
        need(name not in by_name, f"summary.tsv: {name} twice")
        by_name[name] = value
    return {
        "contrasts.tsv": tsv("contrasts.tsv", contrast_rows),
        "readings.tsv": tsv("readings.tsv", reading_rows),
        "species-medians.tsv": tsv("species-medians.tsv", species_rows),
        "rule.tsv": tsv("rule.tsv", per),
        "arms.tsv": tsv("arms.tsv", arm_rows),
        "mix.tsv": tsv("mix.tsv", mix_rows),
        "summary.tsv": tsv("summary.tsv", [{"name": k, "value": v} for k, v in summary]),
    }


# ----------------------------------------------------------------- build from the input set

def load_json(path):
    need(path.is_file(), f"{path.name}: not in the input set")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except ValueError as e:
        raise Stop(f"{path.name}: not JSON ({e})")


def hash_list(path, name):
    """(sha256, time) of one file in a hash list: '<sha256>  <name>', and a comment 'taken at <time>'."""
    text = path.read_text(encoding="utf-8")
    times = re.findall(r"taken at (\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)", text)
    need(len(times) >= 1, f"{path.name}: no time the hashes were taken")
    found = [line.split()[0] for line in text.split("\n")
             if not line.startswith("#") and len(line.split()) == 2 and line.split()[1] == name]
    need(len(found) == 1 and re.fullmatch(r"[0-9a-f]{64}", found[0]) is not None, f"{path.name}: no one hash of {name}")
    return found[0], times[0]


def read_results(I, folder, design, ids):
    """Every result file of one folder: fit rows, and each fit's scores for the comparison of two runs."""
    fits, scores, configs = [], {}, set()
    count = {"fold_file_hash": 0, "same_test_genomes": 0, "training_counts": 0, "evaluator_counts": 0}
    for d in design:
        shared = None
        for arm in ARMS:
            path = I / folder / f"{d['id']}.first-public-{arm}.json"
            r = load_json(path)
            need(r.get("format") == RESULT_FORMAT, f"{folder}/{path.name}: not {RESULT_FORMAT}")
            need(r.get("evaluation_revision") == REVISION, f"{folder}/{path.name}: not revision {REVISION}")
            need(r.get("id") == d["id"] and r.get("scheme") == f"first-public-{arm}", f"{folder}/{path.name}: another dataset or arm")
            need(r.get("fold_file_sha256") == d["fold_files"][arm]["sha256"], f"{folder}/{path.name}: not the design's fold file")
            count["fold_file_hash"] += 1
            samples = [(s["genome_id"], s["resistant"]) for s in r["samples"]]
            for g, t in samples:
                need(isinstance(g, str) and isinstance(t, bool), f"{folder}/{path.name}: a sample without an id or a label")
                ids.add(g)
            if shared is None:
                shared = samples
            need(samples == shared, f"{folder}/{path.name}: the arms test other genomes")
            count["same_test_genomes"] += 1
            truth = [t for _, t in samples]
            need(len(truth) == d["test"]["n"] and sum(truth) == d["test"]["resistant"], f"{folder}/{path.name}: test counts")
            need(len(r["folds"]) == 1, f"{folder}/{path.name}: not one fold")
            fold = r["folds"][0]
            train, train_resistant = d["train_resistant"] + d["train_susceptible"], d["train_resistant"]
            need(fold["train"] == train and fold["train_resistant"] == train_resistant, f"{folder}/{path.name}: training counts")
            need(fold["test"] == len(truth) and fold["test_resistant"] == sum(truth), f"{folder}/{path.name}: the fold's test counts")
            count["training_counts"] += 1
            c = r["eval_config"]
            configs.add((c["seed"], c["screened_features"], c["inner_folds"]))
            learners = r["split"]["learners"]
            need([x["learner"] for x in learners] == list(LEARNERS), f"{folder}/{path.name}: not the five learners in order")
            for x in learners:
                null = x["fold_complexity"][0] == 0 and x["learner"] in ("gbt", "lasso")
                row = fit_row(d["id"], d["species_slug"], arm, x["learner"], truth, x["scores"], train, train_resistant, null)
                m, e = measures(row), x["metrics"]
                need((row["tp"], row["fn"], row["tn"], row["fp"]) == (e["tp"], e["fn"], e["tn"], e["fp"]),
                     f"{folder}/{path.name} {x['learner']}: confusion counts other than the evaluator's")
                need(abs(m["recall"] - e["sensitivity"]) < 1e-12 and abs(m["specificity"] - e["specificity"]) < 1e-12,
                     f"{folder}/{path.name} {x['learner']}: recall or specificity other than the evaluator's")
                count["evaluator_counts"] += 1
                fits.append(row)
                scores[(d["id"], arm, x["learner"])] = x["scores"]
    return fits, scores, configs, count


def build(I):
    """The whole record, as {file name: text}, from the input set."""
    I = Path(I)
    raw = load_json(I / "design.json")
    design = [d for d in raw["datasets"] if d.get("eligible") and d.get("fold_files")]
    need(design, "design.json: no eligible dataset with fold files")
    ids = set()
    fits, after, configs, count = read_results(I, "results", design, ids)
    before_fits, before, before_configs, _ = read_results(I, "before", design, ids)
    need(len(configs) == 1 and configs == before_configs, "the result files hold other settings")
    seed, screened, inner = configs.pop()

    studies = load_json(I / "studies.json")
    per_dataset = studies["per_dataset"]
    need(sorted(per_dataset) == sorted(d["id"] for d in design), "studies.json: other datasets than the design's")
    need(studies["tested_pairs"] == sum(d["test"]["n"] for d in design), "studies.json: tested pairs")
    dataset_rows = []
    for d in design:
        p = per_dataset[d["id"]]
        need(p["drawn_test_studies"] == d["test_studies"] and p["tested_genomes"] == d["test"]["n"],
             f"studies.json {d['id']}: drawn studies or tested genomes other than the design's")
        row = {"dataset": d["id"], "species": d["species_slug"], "antibiotic": d["antibiotic"], "cohort": d["cohort"],
               "cohort_resistant": d["cohort_resistant"], "cut_day": d["cut_day"],
               "share_public_by_cut": float(d["share_on_or_before_cut"]), "drawn_test_studies": d["test_studies"],
               "studies_with_a_tested_genome": p["studies_with_a_tested_genome"], "test": d["test"]["n"],
               "test_resistant": d["test"]["resistant"], "train": d["train_resistant"] + d["train_susceptible"],
               "train_resistant": d["train_resistant"]}
        need(re.fullmatch(r"\d{4}-\d\d-\d\d", d["cut_day"]) is not None, f"{d['id']}: cut day")
        for arm, key in (("past", "past"), ("random", "random"), ("study-out", "study_out")):
            for s in ("past", "later_other_studies", "later_test_studies"):
                row[f"{key}_from_{s}"] = d["train_from"][arm][s]
        dataset_rows.append(row)
    species_rows = []
    for slug in sorted(studies["species"]):
        s = studies["species"][slug]
        mine = [d for d in design if d["species_slug"] == slug]
        need(s["datasets"] == len(mine), f"studies.json {slug}: datasets")
        need(sum(s["test_genomes_per_study"]) == s["distinct_test_genomes"]
             and len(s["test_genomes_per_study"]) == s["distinct_test_studies"], f"studies.json {slug}: study sizes")
        species_rows.append({"species": slug, "datasets": len(mine), "tested_pairs": sum(d["test"]["n"] for d in mine),
                             "tested_resistant_pairs": sum(d["test"]["resistant"] for d in mine),
                             "distinct_test_genomes": s["distinct_test_genomes"],
                             "distinct_test_studies": s["distinct_test_studies"],
                             "largest_study_tested_genomes": max(s["test_genomes_per_study"])})

    # Order: the hash lists and the two run logs.
    plan_sha, plan_time = hash_list(I / "plan" / "PLAN.sha256", "PLAN.md")
    need(sha256_bytes((I / "plan" / "PLAN.md").read_bytes()) == plan_sha, "plan/PLAN.md: not the hashed plan")
    analyze_sha, _ = hash_list(I / "plan" / "PLAN.sha256", "scripts/analyze.py")
    need(sha256_bytes((I / "plan" / "scripts" / "analyze.py").read_bytes()) == analyze_sha,
         "plan/scripts/analyze.py: not the hashed analysis script")
    tables_sha, tables_time = hash_list(I / "plan" / "tables.sha256", "scripts/tables.py")
    need(sha256_bytes((I / "plan" / "scripts" / "tables.py").read_bytes()) == tables_sha,
         "plan/scripts/tables.py: not the hashed layout script")
    rerun_sha, rerun_time = hash_list(I / "rerun" / "RERUN-PLAN.sha256", "RERUN-PLAN.md")
    need(sha256_bytes((I / "rerun" / "RERUN-PLAN.md").read_bytes()) == rerun_sha, "rerun/RERUN-PLAN.md: not the hashed plan")
    stamp = r"(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)"
    first = (I / "logs" / "first-run.log").read_text(encoding="utf-8").split("\n")
    done = [re.match(stamp + r" (\S+)\.(past|random|study-out) done in \d+s$", x) for x in first]
    done = [m for m in done if m]
    need(len(done) == len(design) * len(ARMS), f"logs/first-run.log: {len(done)} results, not {len(design) * len(ARMS)}")
    need(any(re.fullmatch(stamp + " END rc=0", x) for x in first), "logs/first-run.log: no END rc=0")
    rerun = (I / "logs" / "rerun.log").read_text(encoding="utf-8").split("\n")
    runs = [re.match(stamp + " " + stamp + r" rerun-(\S+)\.(past|random|study-out) rc=(\d+) ", x) for x in rerun]
    runs = [m for m in runs if m]
    need(len(runs) == len(design) * len(ARMS), f"logs/rerun.log: {len(runs)} rerun lines, not {len(design) * len(ARMS)}")
    need(all(m.group(5) == "0" for m in runs), "logs/rerun.log: a rerun that did not exit 0")
    need(sorted((m.group(3), m.group(4)) for m in runs) == sorted((d["id"], a) for d in design for a in ARMS),
         "logs/rerun.log: other fold files than the design's")
    order = [("plan_hashed", plan_time, "the plan's hash list"),
             ("first_result", done[0].group(1), "the run log of 14 Nov"),
             ("layout_hashed", tables_time, "the layout script's hash list"),
             ("last_result", done[-1].group(1), "the run log of 14 Nov"),
             ("rerun_plan_hashed", rerun_time, "the rerun plan's hash list"),
             ("rerun_first_fit", runs[0].group(1), "the rerun's run log"),
             ("rerun_last_result", runs[-1].group(2), "the rerun's run log")]
    need(plan_time < done[0].group(1) and rerun_time < runs[0].group(1), "a plan was hashed after its first result")
    need([t for _, t, _ in order] == sorted(t for _, t, _ in order), "the events are not in time order")
    settings = [("cut_share", raw["cut_share"]), ("min_test_class", raw["min_test_class"]),
                ("min_train_class", raw["min_train_class"]), ("min_reserve", raw["min_mates"]),
                ("min_later_other_studies", raw["min_other_later"]), ("seed", seed), ("screened_features", screened),
                ("inner_folds", inner)] + list(FIXED_SETTINGS.items()) + [
        ("plan_sha256", plan_sha), ("analysis_sha256", analyze_sha), ("layout_sha256", tables_sha),
        ("rerun_plan_sha256", rerun_sha)]


    # What the fix changed: scores per fit, and the readings of the same procedure on the 14 Nov fits.
    now, then = readings(fits)[0], readings(before_fits)[0]
    review_rows, files_changed, datasets_changed = [], set(), set()
    for learner in LEARNERS:
        keys = [k for k in after if k[2] == learner]
        diffs = {k: sum(1 for x, y in zip(after[k], before[k]) if x != y) for k in keys}
        for k in keys:
            need(len(after[k]) == len(before[k]), f"{k}: the two runs score other numbers of genomes")
            if diffs[k]:
                files_changed.add(k[:2])
                datasets_changed.add(k[0])
        changed = [f"{a['contrast']} {a['measure']}" for a, b in zip(now, then)
                   if a["learner"] == learner and a["reading"] != b["reading"]]
        review_rows.append({"learner": learner, "fits": len(keys), "fits_changed": sum(1 for v in diffs.values() if v),
                            "scores": sum(len(after[k]) for k in keys), "scores_changed": sum(diffs.values()),
                            "readings": sum(1 for a in now if a["learner"] == learner),
                            "readings_changed": len(changed), "changed_readings": ";".join(changed)})
    totals = [("distinct_test_genomes", studies["distinct_test_genomes"]),
              ("distinct_test_studies", studies["distinct_test_studies"]),
              ("sum_of_per_species_study_counts", studies["sum_of_per_species_study_counts"]),
              ("studies_tested_in_two_species", len(studies["studies_tested_in_two_species"])),
              ("result_files", len(design) * len(ARMS)), ("result_files_changed", len(files_changed)),
              ("datasets_with_a_changed_fit", len(datasets_changed))]

    rows_text = {
        "fits.tsv": tsv("fits.tsv", fits),
        "datasets.tsv": tsv("datasets.tsv", dataset_rows),
        "species.tsv": tsv("species.tsv", species_rows),
        "review.tsv": tsv("review.tsv", review_rows),
        "totals.tsv": tsv("totals.tsv", [{"name": k, "value": v} for k, v in totals]),
        "order.tsv": tsv("order.tsv", [{"event": e, "time": t, "source": s} for e, t, s in order]),
        "settings.tsv": tsv("settings.tsv", [{"name": k, "value": v} for k, v in settings]),
    }
    inputs = sorted(p for p in I.rglob("*") if p.is_file())
    rows_text["inputs.sha256"] = "".join(f"{sha256_bytes(p.read_bytes())}  {p.relative_to(I).as_posix()}\n" for p in inputs)
    checks = [("result_files", count["fold_file_hash"]), ("fold_file_hash_matches_design", count["fold_file_hash"]),
              ("arms_test_the_same_genomes", count["same_test_genomes"]),
              ("training_counts_match_design", count["training_counts"]),
              ("fits_match_evaluator_counts", count["evaluator_counts"]),
              ("before_result_files", len(before_fits) // len(LEARNERS)),
              ("plans_match_their_hash", 2), ("scripts_match_their_hash", 2),
              ("rerun_runs_exit_0", len(runs))]
    full = finish(rows_text, checks, ids)
    return full


def finish(rows_text, checks, ids):
    """Derive, scan for what must not enter, and add checks.tsv, README.md and SHA256SUMS."""
    rows = {n: parse_tsv(n, rows_text[n]) for n in ROW_FILES if n not in ("checks.tsv", "inputs.sha256")}
    derived = derive(rows)
    out = dict(rows_text)
    out.update(derived)
    text = "".join(out.values())
    tokens = set(re.findall(r"[0-9A-Za-z_.]+", text)) | set(re.findall(r"[0-9.]+", text))
    found = len(ids & tokens)
    need(found == 0, f"{found} genome ids of the input set stand in the record")
    for shape in BANNED_SHAPES:
        for name, body in out.items():
            if name in ("inputs.sha256",):
                continue
            need(re.search(shape, body) is None, f"{name}: text of the shape {shape}")
    checks = list(checks) + [("genome_ids_of_the_inputs_in_the_record", found)]
    out["checks.tsv"] = tsv("checks.tsv", [{"check": k, "count": v} for k, v in checks])
    out["README.md"] = readme(out)
    out["SHA256SUMS"] = sums(out)
    return {name: out[name] for name in FILES}


def sums(out):
    return "".join(f"{sha256_bytes(out[n].encode('utf-8'))}  {n}\n" for n in sorted(out) if n != "SHA256SUMS")


# ----------------------------------------------------------------- README

DESCRIBE = {
    "fits.tsv": "Row file. One row per dataset, arm and learner",
    "datasets.tsv": "Row file. One row per dataset",
    "species.tsv": "Row file. One row per species",
    "review.tsv": "Row file. One row per learner: what the review's fix changed",
    "totals.tsv": "Row file. Counts over all datasets",
    "order.tsv": "Row file. One row per event of the plans and runs",
    "settings.tsv": "Row file. The design's minimums, the fits' settings, the settings the script applies, and the sha256 "
                    "of the two plans and two scripts",
    "checks.tsv": "The checks the build ran on the input set",
    "inputs.sha256": "The sha256 of each file of the input set",
    "contrasts.tsv": "Derived. One row per dataset, learner and contrast",
    "readings.tsv": "Derived. One row per learner, contrast and measure",
    "species-medians.tsv": "Derived. One row per learner, contrast, measure and species",
    "rule.tsv": "Derived. The plan's rule per learner",
    "arms.tsv": "Derived. One row per learner and arm, pooled over datasets",
    "mix.tsv": "Derived. One row per arm: where its training genomes come from",
    "summary.tsv": "Derived. Single figures the page quotes",
}
# Each file's columns, named by their headers; readme() fails if a header of the file is not named here.
COLUMN_TEXT = {
    "fits.tsv": "`dataset`, `species`, `arm` and `learner` name the fit. `test` and `test_resistant` count its test "
                "genomes and the resistant ones, and `train` and `train_resistant` its training genomes and the "
                "resistant ones. `tp`, `fn`, `tn` and `fp` count resistant genomes called resistant and called "
                "susceptible, and susceptible genomes called susceptible and called resistant, at a score of 0.5. "
                "`wrong_at_clip` counts wrong calls at the clip and `scores_0_or_1` scores of exactly 0 or 1. "
                "`loss_sum` is the sum of clipped log loss over the test genomes and `clip_loss_sum` the part of it "
                "from wrong calls at the clip. `null_fit` is 1 for a gbt or lasso fit of complexity 0.",
    "datasets.tsv": "`dataset`, `species` and `antibiotic` name the dataset. `cohort` and `cohort_resistant` count the "
                    "genomes `grm learn` evaluates for it and the resistant ones. `cut_day` is its cut day and "
                    "`share_public_by_cut` the share of the cohort public by that day. `drawn_test_studies` counts the "
                    "test studies drawn and `studies_with_a_tested_genome` those with a tested genome. `test`, "
                    "`test_resistant`, `train` and `train_resistant` count the test and training genomes and the "
                    "resistant ones. `past_from_past`, `past_from_later_other_studies` and "
                    "`past_from_later_test_studies` count the past arm's training genomes that are past, later genomes "
                    "of other studies and later genomes of the test studies. `random_from_past`, "
                    "`random_from_later_other_studies`, `random_from_later_test_studies`, `study_out_from_past`, "
                    "`study_out_from_later_other_studies` and `study_out_from_later_test_studies` count the same for "
                    "the random and study-out arms.",
    "species.tsv": "`species` names the species and `datasets` counts its datasets. `tested_pairs` and "
                   "`tested_resistant_pairs` count its tested genome and dataset pairs and the resistant ones, "
                   "`distinct_test_genomes` its distinct tested genomes, `distinct_test_studies` its distinct test "
                   "studies, and `largest_study_tested_genomes` the tested genomes of its largest test study.",
    "review.tsv": "`learner` names the learner. `fits` and `fits_changed` count its fits and those in which the "
                  "review's fix changed a score, and `scores` and `scores_changed` its scores and those that changed. "
                  "`readings` and `readings_changed` count its readings and those that changed, and `changed_readings` "
                  "names each changed one by contrast and measure, with `;` between them.",
    "totals.tsv": "`name` names a count over all datasets and `value` gives it.",
    "order.tsv": "`event` names the event, `time` gives its UTC time and `source` the file it is read from.",
    "settings.tsv": "`name` names a setting and `value` gives it: the design's `cut_share`, its class minimums "
                    "`min_test_class` and `min_train_class`, its pool minimums `min_reserve` and "
                    "`min_later_other_studies`, the fits' `seed`, `screened_features` and `inner_folds`, the settings "
                    "the script applies, and the sha256 of the two plans and two scripts.",
    "checks.tsv": "`check` names a check the build ran on the input set and `count` gives how many passed it. "
                  "`genome_ids_of_the_inputs_in_the_record` gives the count of the input set's genome ids found in the "
                  "record.",
    "inputs.sha256": "Each line gives the sha256 of a file of the input set and its path there.",
    "contrasts.tsv": "`dataset`, `species`, `learner` and `contrast` name the row. `d_recall`, `d_log_loss`, "
                     "`d_specificity` and `d_balanced_accuracy` give the paired difference in each measure.",
    "readings.tsv": "`learner`, `contrast` and `measure` name the reading, and `minuend` and `subtrahend` the two arms "
                    "of the contrast. `median` is the median difference over datasets, `low` and `high` its interval "
                    "from 10,000 species draws, and `reading` what that interval reads. `exact_low`, `exact_high` and "
                    "`exact_reading` give the same over every ordered species draw. `datasets` counts the datasets, and "
                    "`worse`, `better` and `equal` those whose difference is on each side of zero or at zero. "
                    "`minuend_median` and `subtrahend_median` give each arm's median, and `plan` what the plan fixed.",
    "species-medians.tsv": "`learner`, `contrast`, `measure` and `species` name the row. `datasets` counts the "
                           "species' datasets and `median` gives their median difference.",
    "rule.tsv": "`learner` names the learner. `recall` and `log_loss` give its total contrast readings, and `rule` "
                "what the plan's rule makes of them.",
    "arms.tsv": "`learner` and `arm` name the row and `datasets` counts its datasets. `tested_pairs` and "
                "`resistant_pairs` count the pooled test pairs and the resistant ones, and `tp`, `fn`, `tn` and `fp` "
                "the pooled confusion counts. `recall` and `specificity` follow from them, and `wrong_calls` is `fn` "
                "plus `fp`. `wrong_at_clip`, `scores_0_or_1` and `null_fits` add up the fits' cells. "
                "`pooled_log_loss` is the summed `loss_sum` divided by the summed `test`, and `clip_share` the share "
                "of it from wrong calls at the clip.",
    "mix.tsv": "`arm` names the arm and `datasets` counts the datasets. `train` counts its training genomes over all "
               "datasets, `from_past`, `from_later_other_studies` and `from_later_test_studies` how many are past, "
               "later genomes of other studies and later genomes of the test studies, and `past_share` the share that "
               "is past.",
    "summary.tsv": "`name` names a single figure the page quotes and `value` gives it.",
}


def readme(out):
    def count(n):
        if n == "inputs.sha256":
            return out[n].count("\n")
        return out[n].count("\n") - 1

    one_study = sum(1 for line in out["datasets.tsv"].split("\n")[1:-1]
                    if [line.split("\t")[COLUMNS["datasets.tsv"].index(c)] for c in
                        ("drawn_test_studies", "studies_with_a_tested_genome")] == ["1", "1"])
    lines = ["# Time split record", "",
             "This directory holds the record of the split by first-public day inside BV-BRC, rerun after the review. "
             "The page that reads it is [`docs/time-split.md`](../../docs/time-split.md). "
             "`scripts/time-split/record.py` writes every file here from the input set, and its `--check` mode fails "
             "on any byte it would not write.", "",
             "## Files", "", "| File | Rows | Contents |", "|---|---|---|"]
    for n in ROW_FILES + DERIVED_FILES:
        lines.append(f"| `{n}` | {count(n):,} | {DESCRIBE[n]} |")
    lines += ["", "`SHA256SUMS` gives the sha256 of every other file, and this README is what the script writes.", "",
              "## Columns", ""]
    for n in ROW_FILES + DERIVED_FILES:
        if n in COLUMNS:
            missing = [c for c in COLUMNS[n] if f"`{c}`" not in COLUMN_TEXT[n]]
            need(not missing, f"README: the columns {missing} of {n} are not named")
        lines += [f"- `{n}`: {COLUMN_TEXT[n]}"]
    lines += ["", "## Rules", "",
              "**Measures.** A score at or above 0.5 calls resistant. Recall and specificity are fractions of the "
              "confusion counts and balanced accuracy is their mean. Log loss clips the score given to the true class "
              "to [0.0001, 1 - 0.0001], and a fit's log loss is `loss_sum` divided by `test`. A wrong call at the clip "
              "is a resistant genome scored at most 0.0001 or a susceptible genome scored at least 1 - 0.0001.", "",
              "**Contrasts.** total is past minus random, study is study-out minus random, and time is past minus "
              "study-out, each a paired difference on one dataset's test genomes.", "",
              "**Intervals.** `low` and `high` are the 251st and 9,750th of 10,000 medians, each over the datasets of "
              "the species drawn with replacement, from one generator seeded 20261114 and drawn in the order of "
              "`readings.tsv`, as the plan's analysis script draws them. `exact_low` and `exact_high` take the same two "
              "quantiles over every ordered draw of the species, each drawn as often as any other.", "",
              "**Readings.** A median reads worse when its interval lies wholly on the worse side of zero, better when "
              "wholly on the other side, and not shown otherwise. Worse is lower recall, specificity or balanced "
              "accuracy, or a higher log loss. In `readings.tsv`, `plan` is `rule` for the two readings the rule reads, "
              "`fixed` for recall and log loss on the other two contrasts, and `supporting` for specificity and "
              "balanced accuracy, which the plan named as supporting figures.", "",
              "## Limits", "",
              "Each row is per dataset, species, arm, learner, contrast, measure or event, or a single count. No row "
              "is about one genome, and no cell names a genome or a study. "
              f"{one_study} datasets have one drawn and one tested study (`drawn_test_studies` and "
              "`studies_with_a_tested_genome` are 1), so their rows in `fits.tsv`, `contrasts.tsv` and `datasets.tsv` "
              "are results on one unnamed study. In `species.tsv`, `largest_study_tested_genomes` is the count of the "
              "tested genomes of one unnamed study. The result files hold no study, so the study counts rest on the input "
              "set's study count file alone. `--check` alone cannot hold a row file to the result files; that needs "
              "the input set.", "",
              "## Guards", "",
              "`--check` with the input set rebuilds every file and fails on any byte that differs. `--check` alone "
              "reads the row files, derives every derived file and this README again, and fails on any byte that "
              "differs or on a file that differs from `SHA256SUMS`. The input set is not in the repository: it holds "
              "the result files, whose scores are per genome. `--out` writes only to a new path, an empty directory or "
              "an earlier record of this script whose files agree with its own `SHA256SUMS`, and refuses a symlink.",
              "",
              "## Commands", "", "```sh",
              "python3 -I scripts/time-split/record.py --inputs <dir> --out reference/time-split",
              "python3 -I scripts/time-split/record.py --inputs <dir> --check reference/time-split",
              "python3 -I scripts/time-split/record.py --check reference/time-split",
              "python3 -I -m unittest discover -s tests/time_split -t .", "```", ""]
    return "\n".join(lines)


# ----------------------------------------------------------------- record on disk

def read_record(R):
    R = Path(R)
    need(R.is_dir(), f"{R}: no record")
    present = sorted(p.name for p in R.iterdir())
    need(present == sorted(FILES), f"{R}: files {sorted(set(present) ^ set(FILES))} differ from the record's list")
    return {n: (R / n).read_text(encoding="utf-8") for n in FILES}


def check_alone(R):
    """Derive from the record's own row files and compare every derived byte."""
    have = read_record(R)
    rows = {n: parse_tsv(n, have[n]) for n in ROW_FILES if n not in ("checks.tsv", "inputs.sha256")}
    parse_tsv("checks.tsv", have["checks.tsv"])
    need(all(re.fullmatch(r"[0-9a-f]{64}  \S+", x) for x in have["inputs.sha256"].rstrip("\n").split("\n")),
         "inputs.sha256: a line that is not a sha256 and a path")
    want = dict(have)
    want.update(derive(rows))
    want["README.md"] = readme(want)
    want["SHA256SUMS"] = sums(want)
    bad = [n for n in FILES if want[n] != have[n]]
    need(not bad, "differs from what the row files give: " + ", ".join(bad))
    return len(FILES)


def check_out(out):
    """--out takes a new path, an empty directory, or an earlier record of this script: its own files, none a
    symlink, each with the sha256 its own SHA256SUMS gives. Anything else is left as it is."""
    need(not out.is_symlink(), f"{out}: a symlink")
    if not out.exists():
        return
    need(out.is_dir(), f"{out}: not a directory")
    present = sorted(out.iterdir())
    if not present:
        return
    need(sorted(p.name for p in present) == sorted(FILES), f"{out}: holds files this script did not write")
    for p in present:
        need(not p.is_symlink() and p.is_file(), f"{out}/{p.name}: a symlink or not a plain file")
    listed = {}
    for line in (out / "SHA256SUMS").read_bytes().decode("utf-8", "replace").split("\n")[:-1]:
        m = re.fullmatch(r"([0-9a-f]{64})  (\S+)", line)
        need(m is not None and m.group(2) not in listed, f"{out}/SHA256SUMS: not one sha256 and name per line")
        listed[m.group(2)] = m.group(1)
    need(sorted(listed) == sorted(n for n in FILES if n != "SHA256SUMS"),
         f"{out}/SHA256SUMS: does not list the record's files")
    bad = [n for n, h in sorted(listed.items()) if sha256_bytes((out / n).read_bytes()) != h]
    need(not bad, f"{out}: holds files this script did not write, or changed since: " + ", ".join(bad))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--inputs")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--out")
    g.add_argument("--check")
    a = ap.parse_args(argv)
    try:
        if a.out:
            need(a.inputs is not None, "--out needs --inputs")
            out = Path(a.out)
            check_out(out)
            record = build(a.inputs)
            out.mkdir(parents=True, exist_ok=True)
            for n, text in record.items():
                (out / n).write_text(text, encoding="utf-8")
            print(f"wrote {len(record)} files to {out}")
        elif a.inputs:
            record = build(a.inputs)
            have = read_record(a.check)
            bad = [n for n in FILES if have[n] != record[n]]
            need(not bad, "differs from the input set's record: " + ", ".join(bad))
            print(f"{len(FILES)} files rebuilt from the input set, byte for byte")
        else:
            n = check_alone(a.check)
            print(f"{n} files agree with the row files and SHA256SUMS")
    except Stop as e:
        print(f"record.py: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
