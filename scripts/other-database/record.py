#!/usr/bin/env python3
"""The record of the test on data from another public database (reference/other-database).

The script reads the saved files of lanes `outside` (the pilot species) and
`outside2` (the eight other species): their scored genomes, predictions,
labels, funnel and overlap files, request logs and plan files. It also reads
the reports and interval files of the release they scored, and the data root's
drug name list. It fits nothing, runs no `grm` and downloads nothing. Every
measure is recomputed here from the predictions and labels.

    python3 -I scripts/other-database/record.py --lanes L --release R --data-root D --out DIR
    python3 -I scripts/other-database/record.py --lanes L --release R --data-root D --check DIR

L holds the lanes' output directories `outside/` and `outside2/`, R is the
unpacked release and D the public data root. With --out the script writes the
record into DIR. With --check it writes nothing: it exits 0 when DIR holds
exactly the files it would write, byte for byte, and 1 otherwise, naming each
file that differs. A missing file, a file whose sha256 differs from the
lane's own line for it, or saved files that disagree with each other stop
the script with exit 2 and a message.

The row files are species.tsv, datasets.tsv, readings.tsv, genomes.tsv,
order.tsv, plans.tsv, settings.tsv and inputs.sha256. The derived files,
learners.tsv, totals.tsv and differences.tsv, are written by derive() from
the row files' cells alone, so they can be checked from the record without
the inputs.
"""

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
import multiprocessing
import os
import random
import re
import statistics
import sys
from collections import Counter, defaultdict

SPECIES = [
    # The pilot species first, then ADDENDUM.md A1's eight species in its order.
    ("escherichia_coli", "Escherichia coli"),
    ("salmonella_enterica", "Salmonella enterica"),
    ("klebsiella_pneumoniae", "Klebsiella pneumoniae"),
    ("staphylococcus_aureus", "Staphylococcus aureus"),
    ("campylobacter_jejuni", "Campylobacter jejuni"),
    ("campylobacter_coli", "Campylobacter coli"),
    ("neisseria_gonorrhoeae", "Neisseria gonorrhoeae"),
    ("pseudomonas_aeruginosa", "Pseudomonas aeruginosa"),
    ("acinetobacter_baumannii", "Acinetobacter baumannii"),
]
PILOT = "escherichia_coli"
NAMES = dict(SPECIES)
LEARNERS = ("cart", "lasso", "scm")
MIN_CLASS = 25
CUTOFF = 0.5
EPS = 1e-4
DRAWS = 2000
SEED = 1
LEVEL = 0.95
Z = 1.96
# The plan's fixed name table (outside/scripts/survey_overlap.py, MANUAL), used after
# exact equality with a release name and the data root's own name list.
MANUAL = {
    "amoxicillin-clavulanic acid": "amoxicillin/clavulanate",
    "ampicillin-sulbactam": "ampicillin/sulbactam",
    "piperacillin-tazobactam": "piperacillin/tazobactam",
    "trimethoprim-sulfamethoxazole": "sulfamethoxazole/trimethoprim",
    "ceftazidime-avibactam": "ceftazidime/avibactam",
    "ceftolozane-tazobactam": "ceftolozane/tazobactam",
    "ticarcillin-clavulanic acid": "ticarcillin/clavulanate",
    "quinupristin-dalfopristin": "dalfopristin/quinupristin",
    "cefpodoxime-clavulanic acid": "cefpodoxime/clavulanate",
    "sulfisoxazole": "sulfisoxazole",
    "polymyxin B": "polymyxin b",
    "para-aminosalicylic acid": "para-aminosalicylic acid",
}
# The plan files as the lanes fixed them: record name, lane, file, the lane's sha256 line.
PLANS = [
    ("PLAN.md", "outside", "outside/PLAN.md", "outside/results/PLAN.sha256"),
    ("ADDENDUM.md", "outside2", "outside2/ADDENDUM.md", "outside2/ADDENDUM.sha256"),
    ("ADDENDUM2.md", "outside2", "outside2/ADDENDUM2.md", "outside2/ADDENDUM2.sha256"),
]
# Funnel steps in the plan's order, as the lanes' funnel files name them.
STEPS = [
    ("C1", "C1 NCBI AST isolates"),
    ("C2", "C2 no release drug"),
    ("C3", "C3 no GenBank assembly"),
    ("C4", "C4 in BV-BRC labelled set"),
    ("D1", "D1 training accession"),
    ("D2", "D2 training BioSample"),
    ("O1", "O1 outside-study candidate"),
    ("D4", "D4 within-batch duplicate"),
    ("P2", "P2 cap"),
    ("F1", "F1 FASTA not fetched"),
    ("Q1", "Q1 QC"),
    ("Q2_not_judged", "Q2 species check: not judged"),
    ("Q2_fail", "Q2 species check: fail"),
    ("D3", "D3 sequence duplicate"),
]
TIME = re.compile(r"\b(\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)")
SHA = re.compile(r"\b([0-9a-f]{64})\b")


class Stop(Exception):
    """Saved files that are missing or disagree: the record is not written."""


def need(cond, msg):
    if not cond:
        raise Stop(msg)


# ---------------------------------------------------------------- small helpers


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path):
    need(os.path.isfile(path), f"missing file {path}")
    with open(path, encoding="utf-8", newline="") as f:
        lines = [line for line in f.read().splitlines() if line and not line.startswith("#")]
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def comments(path):
    with open(path, encoding="utf-8") as f:
        return [line.rstrip("\n") for line in f if line.startswith("#")]


def f4(v):
    return "" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.4f}"


def f3(v):
    return "" if v is None else f"{v:.3f}"


def acc_order(text):
    """Order of BioProject keys for a tie: letters as text, then the number, lower first.

    PRJNA717739 comes before PRJNA1160997. A key that is not one accession (several
    BioProjects in one field, or a genome with none) comes after every accession."""
    m = re.fullmatch(r"([A-Za-z]+)(\d+)", text)
    return (0, m.group(1), int(m.group(2)), text) if m else (1, text, 0, text)


def largest(counter):
    """The key with the largest count, a tie going to the lower accession."""
    return min(counter, key=lambda k: (-counter[k], acc_order(k)))


def parse_call(text):
    t = (text or "").strip().lower()
    if t in ("resistant", "nonsusceptible", "non-susceptible"):
        return "R"
    if t == "susceptible":
        return "S"
    if t in ("intermediate", "susceptible-dose dependent", "reduced susceptibility"):
        return "I"
    return None


def consolidate(calls):
    if "R" in calls and "S" in calls:
        return "conflict"
    if "R" in calls:
        return "R"
    if "S" in calls:
        return "S"
    if "I" in calls:
        return "I only"
    return "no call"


def loss(y, p):
    p = min(max(p, EPS), 1 - EPS)
    return -math.log(p) if y else -math.log(1 - p)


# ---------------------------------------------------------------- intervals


def _beta_cf(a, b, x):
    """Continued fraction of the regularized incomplete beta (modified Lentz)."""
    tiny = 1e-300
    c, d = 1.0, 1.0 - (a + b) * x / (a + 1.0)
    d = 1.0 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 20000):
        m2 = 2 * m
        for aa in (m * (b - m) * x / ((a + m2 - 1.0) * (a + m2)), -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1.0))):
            d = 1.0 + aa * d
            d = 1.0 / (d if abs(d) > tiny else tiny)
            c = 1.0 + aa / c
            c = c if abs(c) > tiny else tiny
            delta = c * d
            h *= delta
        if abs(delta - 1.0) < 1e-16:
            break
    return h


def beta_cdf(a, b, x):
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x))
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_cf(a, b, x) / a
    return 1.0 - front * _beta_cf(b, a, 1.0 - x) / b


def beta_quantile(q, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if beta_cdf(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(x, n):
    """Exact 95% interval of a proportion x of n."""
    a = (1 - LEVEL) / 2
    lo = 0.0 if x == 0 else beta_quantile(a, x, n - x + 1)
    hi = 1.0 if x == n else beta_quantile(1 - a, x + 1, n - x)
    return lo, hi


def ranks(count):
    """Nearest-rank 2.5th and 97.5th percentiles, as 0-based positions in `count` sorted values."""
    return math.ceil((1 - LEVEL) / 2 * count) - 1, math.ceil((1 + LEVEL) / 2 * count) - 1


def bootstrap(losses):
    """Percentile interval of the mean loss over genome resamples, seed 1."""
    rng = random.Random(SEED)
    n = len(losses)
    means = sorted(sum(rng.choices(losses, k=n)) / n for _ in range(DRAWS))
    i, j = ranks(DRAWS)
    return means[i], means[j]


def stratified(res, sus):
    """Class-balanced log loss interval: each resample draws the resistant, then the susceptible genomes."""
    rng = random.Random(SEED)
    vals = sorted((sum(rng.choices(res, k=len(res))) / len(res) + sum(rng.choices(sus, k=len(sus))) / len(sus)) / 2 for _ in range(DRAWS))
    i, j = ranks(DRAWS)
    return vals[i], vals[j]


def clustered(groups):
    """BioProject cluster bootstrap: recall and log loss intervals over resamples of whole groups.

    groups: key -> (summed loss, genomes, resistant genomes, resistant genomes called resistant)."""
    keys = sorted(groups)
    rng = random.Random(SEED)
    recs, lls = [], []
    for _ in range(DRAWS):
        tot, n, r, tp = 0.0, 0, 0, 0
        for k in rng.choices(keys, k=len(keys)):
            g = groups[k]
            tot += g[0]
            n += g[1]
            r += g[2]
            tp += g[3]
        lls.append(tot / n)
        if r:
            recs.append(tp / r)
    recs.sort()
    lls.sort()
    i, j = ranks(len(recs))
    k, m = ranks(DRAWS)
    return recs[i], recs[j], lls[k], lls[m]


def sample_sd(xs):
    mean = sum(xs) / len(xs)
    return math.sqrt(sum((x - mean) ** 2 for x in xs) / (len(xs) - 1))


# ---------------------------------------------------------------- one dataset


def lineage_figure(report, learner):
    """The release's lineage split for one learner, from its held-out scores; None without one."""
    split = next((s for s in report["splits"] if s["scheme"] == "lineage"), None)
    entry = None if split is None else next((x for x in split["learners"] if x["learner"] == learner), None)
    if entry is None:
        return None
    need(len(entry["scores"]) == len(report["samples"]), f"{report['id']} {learner}: scores and samples differ in length")
    pairs = [(1 if s["resistant"] else 0, p) for s, p in zip(report["samples"], entry["scores"]) if p is not None]
    r = sum(y for y, _ in pairs)
    s = len(pairs) - r
    need(r and s, f"{report['id']} {learner}: the lineage split scores no genome of a class")
    tp = sum(1 for y, p in pairs if y and p >= CUTOFF)
    tn = sum(1 for y, p in pairs if not y and p < CUTOFF)
    losses = [loss(y, p) for y, p in pairs]
    recall, specificity = tp / r, tn / s
    published = entry.get("metrics", {})
    for key, mine in (("sensitivity", recall), ("specificity", specificity)):
        if published.get(key) is not None:
            need(abs(published[key] - mine) < 1e-9, f"{report['id']} {learner}: lineage {key} {mine} differs from the report's {published[key]}")
    res = [loss(y, p) for y, p in pairs if y]
    sus = [loss(y, p) for y, p in pairs if not y]
    return {
        "recall": recall,
        "specificity": specificity,
        "log_loss": sum(losses) / len(losses),
        "losses": losses,
        "tp": tp,
        "r": r,
        "balanced": (sum(res) / len(res) + sum(sus) / len(sus)) / 2,
    }


def score_dataset(task):
    """All rows of one dataset: one datasets.tsv row per learner scored, one readings.tsv row per estimable one."""
    with open(task["report"], encoding="utf-8") as f:
        report = json.load(f)
    need(report["id"] == task["dataset"], f"{task['report']} holds {report['id']}")
    intervals = {}
    if os.path.isfile(task["intervals"]):
        with open(task["intervals"], encoding="utf-8") as f:
            for split in json.load(f)["splits"]:
                for entry in split["learners"]:
                    intervals[(split.get("scheme"), entry["learner"])] = entry
    cohort = sum(1 for s in report["samples"] if s["resistant"]) / len(report["samples"])
    labels, bioproject = task["labels"], task["bioprojects"]
    genomes = sorted(labels)
    tally = Counter(labels.values())
    rows, readings = [], []
    for learner in LEARNERS:
        if learner not in task["predictions"]:
            continue
        calls = task["predictions"][learner]
        pairs, null = [], 0
        for g in genomes:
            if labels[g] not in ("R", "S"):
                continue
            call, score = calls.get(g, ("", ""))
            if call == "":
                null += 1
                continue
            need(call in ("resistant", "susceptible"), f"{task['dataset']} {g} {learner}: call {call!r}")
            c = 1 if call == "resistant" else 0
            p = float(c) if learner == "scm" else float(score)
            need((p >= CUTOFF) == bool(c), f"{task['dataset']} {g} {learner}: score {p} and call {call} disagree")
            pairs.append((g, 1 if labels[g] == "R" else 0, p, c))
        nr = sum(y for _, y, _, _ in pairs)
        ns = len(pairs) - nr
        tp = sum(1 for _, y, _, c in pairs if y and c)
        fn = nr - tp
        tn = sum(1 for _, y, _, c in pairs if not y and not c)
        fp = ns - tn
        lin = lineage_figure(report, learner)
        estimable = nr >= MIN_CLASS and ns >= MIN_CLASS and lin is not None
        if nr < MIN_CLASS or ns < MIN_CLASS:
            note = "fewer than 25 R or 25 S"
        elif lin is None:
            note = "no lineage figure in the release"
        else:
            note = ""
        groups = defaultdict(list)
        for g, y, p, c in pairs:
            groups[bioproject[g] or "none:" + g].append((y, loss(y, p), c))
        top = largest(Counter({k: len(v) for k, v in groups.items()})) if groups else None
        row = {
            "species": task["species"],
            "part": task["part"],
            "dataset": task["dataset"],
            "antibiotic": task["antibiotic"],
            "learner": learner,
            "scored": str(len(genomes)),
            "labels_R": str(tally["R"]),
            "labels_S": str(tally["S"]),
            "labels_I_only": str(tally["I only"]),
            "labels_no_call": str(tally["no call"]),
            "labels_conflict": str(tally["conflict"]),
            "null_calls": str(null),
            "R": str(nr),
            "S": str(ns),
            "TP": str(tp),
            "FN": str(fn),
            "TN": str(tn),
            "FP": str(fp),
            "estimable": "yes" if estimable else "no",
            "note": note,
            "lineage_figure": "no" if lin is None else "yes",
            "lineage_recall": "" if lin is None else f4(lin["recall"]),
            "lineage_specificity": "" if lin is None else f4(lin["specificity"]),
            "lineage_log_loss": "" if lin is None else f4(lin["log_loss"]),
            "cohort_resistant_share": f3(cohort),
            "resistant_share": f3(nr / len(pairs)) if pairs else "",
            "bioprojects": str(len(groups)),
            "largest_bioproject": "" if top is None else top,
            "largest_bioproject_share": "" if top is None else f3(len(groups[top]) / len(pairs)),
        }
        for k in ("recall", "recall_lo", "recall_hi", "specificity", "log_loss", "log_loss_lo", "log_loss_hi", "recall_reading", "log_loss_reading", "holds"):
            row[k] = ""
        rows.append(row)
        if not estimable:
            continue
        recall, specificity = tp / nr, tn / ns
        lo, hi = clopper_pearson(tp, nr)
        losses = [loss(y, p) for _, y, p, _ in pairs]
        ll = sum(losses) / len(losses)
        llo, lhi = bootstrap(losses)
        rr = "below" if hi < lin["recall"] else ("above" if lo > lin["recall"] else "consistent")
        lr = "worse" if llo > lin["log_loss"] else ("better" if lhi < lin["log_loss"] else "consistent")
        row.update(
            {
                "recall": f4(recall),
                "recall_lo": f4(lo),
                "recall_hi": f4(hi),
                "specificity": f4(specificity),
                "log_loss": f4(ll),
                "log_loss_lo": f4(llo),
                "log_loss_hi": f4(lhi),
                "recall_reading": rr,
                "log_loss_reading": lr,
                "holds": "yes" if rr != "below" and lr != "worse" else "no",
            }
        )
        # Reading L: the lineage figure with its own interval.
        iv = intervals.get(("lineage", learner))
        if iv is not None:
            l_rlo, l_rhi, source = 1 - iv["vme"]["hi"], 1 - iv["vme"]["lo"], "release interval"
        else:
            (l_rlo, l_rhi), source = clopper_pearson(lin["tp"], lin["r"]), "released scores"
        l_llo, l_lhi = bootstrap(lin["losses"])
        l_r = "below" if hi < l_rlo else ("above" if lo > l_rhi else "consistent")
        l_l = "worse" if llo > l_lhi else ("better" if lhi < l_llo else "consistent")
        # Reading C: class-balanced log loss, stratified bootstrap, the lineage figure fixed.
        res = [loss(y, p) for _, y, p, _ in pairs if y]
        sus = [loss(y, p) for _, y, p, _ in pairs if not y]
        bal = (sum(res) / len(res) + sum(sus) / len(sus)) / 2
        clo, chi = stratified(res, sus)
        c_l = "worse" if clo > lin["balanced"] else ("better" if chi < lin["balanced"] else "consistent")
        # Reading D: BioProject cluster bootstrap on this side, the lineage figure fixed.
        sums = {k: (sum(x[1] for x in v), len(v), sum(x[0] for x in v), sum(x[0] * x[2] for x in v)) for k, v in groups.items()}
        drlo, drhi, dllo, dlhi = clustered(sums)
        d_r = "below" if drhi < lin["recall"] else ("above" if drlo > lin["recall"] else "consistent")
        d_l = "worse" if dllo > lin["log_loss"] else ("better" if dlhi < lin["log_loss"] else "consistent")
        reading = {
            "species": task["species"],
            "part": task["part"],
            "dataset": task["dataset"],
            "learner": learner,
            "status": "exploratory" if task["part"] == "pilot" else "fixed in advance",
            "L_lineage_recall_lo": f4(l_rlo),
            "L_lineage_recall_hi": f4(l_rhi),
            "L_lineage_recall_source": source,
            "L_lineage_log_loss_lo": f4(l_llo),
            "L_lineage_log_loss_hi": f4(l_lhi),
            "L_recall": l_r,
            "L_log_loss": l_l,
            "L_holds": "yes" if l_r != "below" and l_l != "worse" else "no",
            "C_log_loss_value": f4(bal),
            "C_log_loss_lo": f4(clo),
            "C_log_loss_hi": f4(chi),
            "C_lineage_log_loss": f4(lin["balanced"]),
            "C_log_loss": c_l,
            "C_holds": "yes" if rr != "below" and c_l != "worse" else "no",
            "D_recall_lo": f4(drlo),
            "D_recall_hi": f4(drhi),
            "D_log_loss_lo": f4(dllo),
            "D_log_loss_hi": f4(dlhi),
            "D_recall": d_r,
            "D_log_loss": d_l,
            "D_holds": "yes" if d_r != "below" and d_l != "worse" else "no",
            "T_recall": "",
            "T_log_loss": "",
            "T_holds": "",
            "errors": str(fn + fp),
            "largest_bioproject_errors_share": "",
        }
        errors = Counter(bioproject[g] or "none:" + g for g, y, _, c in pairs if y != c)
        if fn + fp:
            reading["largest_bioproject_errors_share"] = f3(errors[top] / (fn + fp))
        if task["part"] == "pilot":
            # Reading T of the pilot's review: a two-sample 95% interval of each difference.
            se_r = math.sqrt(recall * (1 - recall) / nr + lin["recall"] * (1 - lin["recall"]) / lin["r"])
            dr = recall - lin["recall"]
            t_r = "below" if dr + Z * se_r < 0 else ("above" if dr - Z * se_r > 0 else "consistent")
            se_l = math.sqrt(sample_sd(losses) ** 2 / len(losses) + sample_sd(lin["losses"]) ** 2 / len(lin["losses"]))
            dl = ll - lin["log_loss"]
            t_l = "worse" if dl - Z * se_l > 0 else ("better" if dl + Z * se_l < 0 else "consistent")
            reading.update({"T_recall": t_r, "T_log_loss": t_l, "T_holds": "yes" if t_r != "below" and t_l != "worse" else "no"})
        readings.append(reading)
    return rows, readings


# ---------------------------------------------------------------- one species


def species_dir(lanes, slug):
    return os.path.join(lanes, "outside", "results", "pilot") if slug == PILOT else os.path.join(lanes, "outside2", "results", slug)


def report_pins(lanes, slug):
    path = os.path.join(lanes, "outside", "results", "pilot-reports.sha256") if slug == PILOT else os.path.join(species_dir(lanes, slug), "reports.sha256")
    pins = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                digest, name = line.split()
                pins[name[: -len(".json")]] = digest
    return path, pins


def funnel_counts(sdir, slug):
    removed, remaining, order = {}, {}, []
    for name in ("funnel-label-blind.tsv", "funnel-sequence.tsv"):
        for row in read_tsv(os.path.join(sdir, name)):
            step = row["step"]
            if step == "after P2":
                need(int(row["remaining"]) == remaining["P2"], f"{slug}: the sequence funnel starts at {row['remaining']}")
                continue
            code = next((c for c, label in STEPS if step.startswith(label)), None)
            need(code is not None, f"{slug}: unknown funnel step {step!r}")
            removed[code], remaining[code] = int(row["removed"]), int(row["remaining"])
            order.append(code)
    need(order[0] == "C1", f"{slug}: the funnel does not start at C1")
    for prev, code in zip(order, order[1:]):
        need(remaining[prev] - removed[code] == remaining[code], f"{slug}: funnel step {code} does not add up")
    return removed, remaining, order


def read_species(lanes, release, slug, all_drugs, rawmap, drug_of, inputs):
    sdir = species_dir(lanes, slug)
    part = "pilot" if slug == PILOT else "eight"
    files = ["selected.tsv", "filter.tsv", "predictions.tsv.gz", "labels.jsonl.gz", "funnel-label-blind.tsv", "funnel-sequence.tsv", "training-accessions.tsv"]
    for name in files:
        path = os.path.join(sdir, name)
        need(os.path.isfile(path), f"missing file {path}")
        inputs[path] = sha256_file(path)
    pins_path, pins = report_pins(lanes, slug)
    inputs[pins_path] = sha256_file(pins_path)
    selected = {r["assembly"]: r for r in read_tsv(os.path.join(sdir, "selected.tsv")) if r["assembly"]}
    filt = {r["assembly"]: r for r in read_tsv(os.path.join(sdir, "filter.tsv"))}
    training = read_tsv(os.path.join(sdir, "training-accessions.tsv"))
    removed, remaining, order = funnel_counts(sdir, slug)

    preds = defaultdict(lambda: defaultdict(dict))
    sha_of = {}
    with gzip.open(os.path.join(sdir, "predictions.tsv.gz"), "rt", encoding="utf-8", newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            need(r["genome"] in selected, f"{slug}: prediction for {r['genome']}, which is not selected")
            preds[r["dataset"]][r["learner"]][r["genome"]] = (r["call"], r["score"])
            sha_of.setdefault(r["dataset"], set()).add(r["report_sha256"])
    scored = sorted({g for d in preds.values() for calls in d.values() for g in calls})
    kept = sorted(a for a, r in filt.items() if not r["removed_at"])
    need(scored == kept, f"{slug}: the genomes scored are not the genomes the filter kept")
    need(len(scored) == remaining[order[-1]], f"{slug}: {len(scored)} genomes scored, the funnel ends at {remaining[order[-1]]}")
    need(sum(1 for r in filt.values() if r["removed_at"].startswith("D3")) == removed["D3"], f"{slug}: D3 count differs from filter.tsv")

    calls = defaultdict(set)
    with gzip.open(os.path.join(sdir, "labels.jsonl.gz"), "rt", encoding="utf-8") as f:
        for line in f:
            x = json.loads(line)
            name = x.get("antibiotic") or ""
            drug = name if name in all_drugs else rawmap.get(name.strip().lower()) or MANUAL.get(name, "")
            if not drug:
                continue
            key = (x["biosample_acc"].upper(), drug)
            c = parse_call(x.get("phenotype"))
            calls[key]  # every row of a release drug makes the pair known
            if c:
                calls[key].add(c)

    tasks = []
    for ds in sorted(preds):
        need(ds in drug_of, f"{slug}: {ds} is not a dataset of the release")
        report = os.path.join(release, "datasets", ds + ".json")
        need(os.path.isfile(report), f"missing file {report}")
        digest = sha256_file(report)
        need(pins.get(ds) == digest, f"{ds}: the release report's sha256 differs from the lane's pin")
        need(sha_of[ds] == {digest}, f"{ds}: the predictions name another report")
        inputs[report] = digest
        ivpath = os.path.join(release, "intervals", ds + ".json")
        if os.path.isfile(ivpath):
            inputs[ivpath] = sha256_file(ivpath)
        genomes = sorted({g for c in preds[ds].values() for g in c})
        tasks.append(
            {
                "species": slug,
                "part": part,
                "dataset": ds,
                "antibiotic": drug_of[ds],
                "report": report,
                "intervals": ivpath,
                "predictions": {ln: dict(c) for ln, c in preds[ds].items()},
                "labels": {g: consolidate(calls.get((selected[g]["biosample"].upper(), drug_of[ds]), set())) for g in genomes},
                "bioprojects": {g: selected[g]["bioprojects"] for g in genomes},
            }
        )

    dist = [float(filt[g]["nearest_distance"]) for g in scored]
    q = statistics.quantiles(dist, n=4)
    projects = Counter(p for g in scored for p in set(selected[g]["bioprojects"].split(";")) if p)
    top = largest(projects)
    training_projects = Counter(p for r in training for p in set(r["bioprojects"].split(";")) if p)
    row = {"species": slug, "part": part}
    for code, _ in STEPS:
        row[code] = str(removed[code]) if code in removed else ""
    row["C1"] = str(remaining["C1"])
    row.update(
        {
            "scored": str(len(scored)),
            "reports_scored": str(len(preds)),
            "dist_0.0005_0.001": str(sum(1 for x in dist if 0.0005 < x <= 0.001)),
            "dist_0.001_0.005": str(sum(1 for x in dist if 0.001 < x <= 0.005)),
            "dist_above_0.005": str(sum(1 for x in dist if x > 0.005)),
            "dist_quartiles": f"{q[0]:.6f};{statistics.median(dist):.6f};{q[2]:.6f}",
            "shared_bioproject": str(sum(1 for g in scored if selected[g]["bioproject_shared_with_training"] == "yes")),
            "outside_study_overlap": str(sum(1 for g in scored if selected[g]["outside_study_overlap"] == "yes")),
            "training_genomes": str(len(training)),
            "training_no_bioproject": str(sum(1 for r in training if not r["bioprojects"].strip())),
            "bioprojects_scored": str(len(projects)),
            "largest_bioproject": top,
            "largest_bioproject_scored": str(projects[top]),
            "largest_bioproject_share": f3(projects[top] / len(scored)),
            "largest_bioproject_training": str(training_projects.get(top, 0)),
        }
    )
    genomes = [{"species": slug, "assembly": g, "biosample": selected[g]["biosample"], "bioprojects": selected[g]["bioprojects"]} for g in scored]
    return row, tasks, genomes, {selected[g]["biosample"]: slug for g in selected}, {selected[g]["biosample"] for g in scored}


# ---------------------------------------------------------------- order of predictions and labels


def chain_times(path):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    m = re.search(r"predictions written (\S+)", text)
    need(m is not None, f"{path}: no 'predictions written' line")
    return m.group(1)


def order_rows(lanes, species_of, scored_bs, plans, inputs):
    rows = []
    requests = os.path.join(lanes, "outside2", "results", "requests.tsv")
    need(os.path.isfile(requests), f"missing file {requests}")
    inputs[requests] = sha256_file(requests)
    per = defaultdict(list)
    asked = defaultdict(set)
    with open(requests, encoding="utf-8", newline="") as f:
        for r in csv.reader(f, delimiter="\t"):
            if not r or r[0] == "time" or not r[1].startswith("labels "):
                continue
            params = json.loads(r[6])
            need(params.get("collection") == "ast", f"label request at {r[0]} is not to the AST collection")
            bs = re.findall(r'"([A-Z]+\d+)"', params["fq"])
            owners = {species_of.get(b, "?") for b in bs}
            need(len(owners) == 1 and "?" not in owners, f"label request at {r[0]} names BioSamples of {sorted(owners)}")
            sp = owners.pop()
            per[sp].append((r[0], r[4]))
            asked[sp] |= set(bs)
    fixed = {p["file"]: p["fixed"] for p in plans}
    for slug, _ in SPECIES:
        if slug == PILOT:
            log = os.path.join(lanes, "outside", "results", "logs", "chain.log")
            probes = os.path.join(lanes, "outside", "results", "probes.tsv")
            for path in (log, probes):
                need(os.path.isfile(path), f"missing file {path}")
                inputs[path] = sha256_file(path)
            probe_rows = read_tsv(probes)
            reqs = [(r["utc"], r["http"]) for r in probe_rows if r["label"].startswith("labels ")]
            assigned, n_asked, n_scored, equal, plan_file = "by time", "", "", "", "PLAN.md"
            # Requests before the plan that the log holds with no parameters (the pilot review's S5).
            bare = sorted((r["utc"], r["bytes"]) for r in probe_rows if r.get("error") == "POST-test" and r["utc"] < fixed[plan_file])
            need(bare, f"{probes}: no request logged without parameters")
            unrecorded = {"unrecorded_before_plan": str(len(bare)), "unrecorded_first": bare[0][0], "unrecorded_last": bare[-1][0],
                          "unrecorded_bytes": ";".join(b for _, b in bare)}
        else:
            log = os.path.join(lanes, "outside2", "results", "logs", slug + ".chain.log")
            need(os.path.isfile(log), f"missing file {log}")
            inputs[log] = sha256_file(log)
            reqs = per[slug]
            assigned, plan_file = "by BioSample", "ADDENDUM2.md"
            unrecorded = dict.fromkeys(("unrecorded_before_plan", "unrecorded_first", "unrecorded_last", "unrecorded_bytes"), "")
            n_asked, n_scored = str(len(asked[slug])), str(len(scored_bs[slug]))
            equal = "yes" if asked[slug] == scored_bs[slug] else "no"
        written = chain_times(log)
        times = sorted(t for t, _ in reqs)
        need(times, f"{slug}: no label request")
        rows.append(
            {
                "species": slug,
                "assigned": assigned,
                "predictions_written": written,
                "first_label_request": times[0],
                "last_label_request": times[-1],
                "label_requests": str(len(reqs)),
                "answered_200": str(sum(1 for _, s in reqs if s == "200")),
                "biosamples_asked": n_asked,
                "scored_biosamples": n_scored,
                "asked_equals_scored": equal,
                "first_at_or_after_predictions": "yes" if times[0] >= written else "no",
                "plan_file": plan_file,
                "first_after_plan_file": "yes" if times[0] > fixed[plan_file] else "no",
                **unrecorded,
            }
        )
    return rows


# ---------------------------------------------------------------- derived files


def margin(holds, estimable):
    """Datasets that must change for the plan's rule to give the other result; '' when not assessable."""
    if estimable < 3:
        return ""
    return str(holds - estimable // 2) if 2 * holds > estimable else str(estimable // 2 + 1 - holds)


def rule(holds, estimable):
    if estimable < 3:
        return "not assessable"
    return "meets" if 2 * holds > estimable else "fails"


def median_text(xs):
    if not xs:
        return ""
    m = statistics.median(xs)
    return str(int(m)) if m == int(m) else f"{m:.1f}"


def derive(datasets, readings):
    """learners.tsv, totals.tsv and differences.tsv from the cells of datasets.tsv and readings.tsv."""
    read = {(r["dataset"], r["learner"]): r for r in readings}
    learners = []
    for slug in sorted(NAMES):
        for ln in LEARNERS:
            est = [r for r in datasets if r["species"] == slug and r["learner"] == ln and r["estimable"] == "yes"]
            part = "pilot" if slug == PILOT else "eight"
            e = len(est)
            h = sum(r["holds"] == "yes" for r in est)
            spec = [(float(r["specificity"]), r["dataset"]) for r in est]
            low = min(spec) if spec else None
            row = {
                "species": slug,
                "learner": ln,
                "part": part,
                "estimable": str(e),
                "holds": str(h),
                "recall_below": str(sum(r["recall_reading"] == "below" for r in est)),
                "log_loss_worse": str(sum(r["log_loss_reading"] == "worse" for r in est)),
                "specificity_under_half": str(sum(s < 0.5 for s, _ in spec)),
                "lowest_specificity": "" if low is None else f"{low[0]:.4f}",
                "lowest_specificity_dataset": "" if low is None else low[1],
                "holds_with_specificity_under_half": ";".join(r["dataset"] for r in est if r["holds"] == "yes" and float(r["specificity"]) < 0.5),
                "rule": rule(h, e),
                "margin": margin(h, e),
            }
            for v in ("L", "C", "D", "T"):
                if v == "T" and part != "pilot":
                    row["T_holds"], row["T_rule"] = "", ""
                    continue
                hv = sum(read[(r["dataset"], ln)][v + "_holds"] == "yes" for r in est)
                row[v + "_holds"], row[v + "_rule"] = str(hv), rule(hv, e)
            shares = [(r["resistant_share"], r["cohort_resistant_share"]) for r in est]
            row.update(
                {
                    "median_bioprojects": median_text([int(r["bioprojects"]) for r in est]),
                    "share_below_cohort": str(sum(float(a) < float(b) for a, b in shares)),
                    "share_above_cohort": str(sum(float(a) > float(b) for a, b in shares)),
                    "readings_status": "exploratory" if part == "pilot" else "fixed in advance",
                }
            )
            learners.append(row)

    totals = []
    eight = [r for r in learners if r["part"] == "eight"]
    for ln in LEARNERS:
        rs = [r for r in eight if r["learner"] == ln]
        row = {"learner": ln, "species": str(len(rs))}
        for v, key in (("plan", "rule"), ("L", "L_rule"), ("C", "C_rule"), ("D", "D_rule")):
            c = Counter(r[key] for r in rs)
            row[f"{v}_meets"], row[f"{v}_fails"], row[f"{v}_not_assessable"] = str(c["meets"]), str(c["fails"]), str(c["not assessable"])
        for k in ("estimable", "holds", "recall_below", "log_loss_worse", "specificity_under_half"):
            row[k] = str(sum(int(r[k]) for r in rs))
        row["results_one_from_other"] = str(sum(r["margin"] == "1" for r in rs))
        totals.append(row)

    differences = []
    by = {(r["species"], r["learner"]): r for r in eight}
    slugs = sorted(s for s in NAMES if s != PILOT)
    pairs = [(a, b) for i, a in enumerate(LEARNERS) for b in LEARNERS[i + 1 :]]
    for slug in slugs + ["~ count over the eight species"]:
        for a, b in pairs:
            if slug.startswith("~"):
                # Species meeting the rule, counted over the species where both learners can be read.
                va, vb, signs = [], [], set()
                for key in ("rule", "L_rule", "C_rule", "D_rule"):
                    both = [sp for sp in slugs if "not assessable" not in (by[(sp, a)][key], by[(sp, b)][key])]
                    x = sum(by[(sp, a)][key] == "meets" for sp in both)
                    y = sum(by[(sp, b)][key] == "meets" for sp in both)
                    va.append(f"{x} of {len(both)}")
                    vb.append(f"{y} of {len(both)}")
                    signs.add((x > y) - (x < y))
            else:
                score = {"meets": 1, "fails": 0}
                ra = [by[(slug, a)][k] for k in ("rule", "L_rule", "C_rule", "D_rule")]
                rb = [by[(slug, b)][k] for k in ("rule", "L_rule", "C_rule", "D_rule")]
                if "not assessable" in ra + rb:
                    signs = {None}
                else:
                    signs = {(score[x] > score[y]) - (score[x] < score[y]) for x, y in zip(ra, rb)}
                va, vb = ra, rb
            differences.append(
                {
                    "species": slug.lstrip("~ ") if slug.startswith("~") else slug,
                    "learner_a": a,
                    "learner_b": b,
                    "a": ";".join(map(str, va)),
                    "b": ";".join(map(str, vb)),
                    "differs_in_all_four": "yes" if len(signs) == 1 and signs != {0} and signs != {None} else "no",
                }
            )
    return learners, totals, differences


# ---------------------------------------------------------------- writing


DATASET_COLS = [
    "species", "dataset", "learner", "part", "antibiotic", "scored", "labels_R", "labels_S", "labels_I_only", "labels_no_call",
    "labels_conflict", "null_calls", "R", "S", "TP", "FN", "TN", "FP", "estimable", "note", "recall", "recall_lo", "recall_hi",
    "specificity", "log_loss", "log_loss_lo", "log_loss_hi", "lineage_figure", "lineage_recall", "lineage_specificity",
    "lineage_log_loss", "recall_reading", "log_loss_reading", "holds", "resistant_share", "cohort_resistant_share",
    "bioprojects", "largest_bioproject", "largest_bioproject_share",
]
READING_COLS = [
    "species", "dataset", "learner", "part", "status", "L_lineage_recall_lo", "L_lineage_recall_hi", "L_lineage_recall_source",
    "L_lineage_log_loss_lo", "L_lineage_log_loss_hi", "L_recall", "L_log_loss", "L_holds", "C_log_loss_value", "C_log_loss_lo",
    "C_log_loss_hi", "C_lineage_log_loss", "C_log_loss", "C_holds", "D_recall_lo", "D_recall_hi", "D_log_loss_lo", "D_log_loss_hi",
    "D_recall", "D_log_loss", "D_holds", "T_recall", "T_log_loss", "T_holds", "errors", "largest_bioproject_errors_share",
]
SPECIES_COLS = (
    ["species", "part"]
    + [c for c, _ in STEPS]
    + ["scored", "reports_scored", "dist_0.0005_0.001", "dist_0.001_0.005", "dist_above_0.005", "dist_quartiles", "shared_bioproject",
       "outside_study_overlap", "training_genomes", "training_no_bioproject", "bvbrc_labelled", "bvbrc_no_accession", "no_accession_share",
       "bioprojects_scored", "largest_bioproject", "largest_bioproject_scored", "largest_bioproject_share", "largest_bioproject_training"]
)
LEARNER_COLS = [
    "species", "learner", "part", "estimable", "holds", "recall_below", "log_loss_worse", "specificity_under_half", "lowest_specificity",
    "lowest_specificity_dataset", "holds_with_specificity_under_half", "rule", "margin", "L_holds", "L_rule", "C_holds", "C_rule",
    "D_holds", "D_rule", "T_holds", "T_rule", "median_bioprojects", "share_below_cohort", "share_above_cohort", "readings_status",
]
TOTAL_COLS = [
    "learner", "species", "plan_meets", "plan_fails", "plan_not_assessable", "L_meets", "L_fails", "L_not_assessable", "C_meets",
    "C_fails", "C_not_assessable", "D_meets", "D_fails", "D_not_assessable", "estimable", "holds", "recall_below", "log_loss_worse",
    "specificity_under_half", "results_one_from_other",
]
DIFF_COLS = ["species", "learner_a", "learner_b", "a", "b", "differs_in_all_four"]
GENOME_COLS = ["species", "assembly", "biosample", "bioprojects"]
ORDER_COLS = [
    "species", "assigned", "predictions_written", "first_label_request", "last_label_request", "label_requests", "answered_200",
    "biosamples_asked", "scored_biosamples", "asked_equals_scored", "first_at_or_after_predictions", "plan_file", "first_after_plan_file",
    "unrecorded_before_plan", "unrecorded_first", "unrecorded_last", "unrecorded_bytes",
]
PLAN_COLS = ["file", "lane", "sha256", "fixed", "time_from"]
# The leading columns each TSV is sorted by.
KEYS = {
    "species.tsv": ("species",),
    "datasets.tsv": ("species", "dataset", "learner"),
    "readings.tsv": ("species", "dataset", "learner"),
    "learners.tsv": ("species", "learner"),
    "totals.tsv": ("learner",),
    "differences.tsv": ("species", "learner_a", "learner_b"),
    "genomes.tsv": ("species", "assembly"),
    "order.tsv": ("species",),
    "plans.tsv": ("file",),
    "settings.tsv": ("key",),
}
COLUMNS = {
    "species.tsv": SPECIES_COLS,
    "datasets.tsv": DATASET_COLS,
    "readings.tsv": READING_COLS,
    "learners.tsv": LEARNER_COLS,
    "totals.tsv": TOTAL_COLS,
    "differences.tsv": DIFF_COLS,
    "genomes.tsv": GENOME_COLS,
    "order.tsv": ORDER_COLS,
    "plans.tsv": PLAN_COLS,
    "settings.tsv": ["key", "value"],
}


def tsv(name, rows):
    cols = COLUMNS[name]
    keys = KEYS[name]
    out = io.StringIO()
    out.write("\t".join(cols) + "\n")
    for r in sorted(rows, key=lambda r: tuple(r[k] for k in keys)):
        cells = [r[c] for c in cols]
        need(all("\t" not in c and "\n" not in c for c in cells), f"{name}: a cell holds a tab or a line end")
        out.write("\t".join(cells) + "\n")
    return out.getvalue()


def parse(text):
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def settings_rows(release_id, revision):
    pairs = {
        "release": release_id,
        "evaluation_revision": str(revision),
        "learners_scored": "cart;lasso;scm",
        "learners_not_scored": "gbt;rf (the release holds their importances only)",
        "min_class": str(MIN_CLASS),
        "cutoff": str(CUTOFF),
        "clip": "0.0001",
        "level": str(LEVEL),
        "recall_interval": "Clopper-Pearson",
        "bootstrap_draws": str(DRAWS),
        "bootstrap_seed": str(SEED),
        "bootstrap_generator": "Python random.Random(seed).choices",
        "percentiles": "nearest rank: sorted values 50 and 1950 of 2000",
        "estimable": "at least 25 R and 25 S labelled genomes with a call, and a lineage figure in the release",
        "holds": "recall not below and log loss not worse than the release's lineage split",
        "pass_rule": "meets when it holds on more than half of at least 3 estimable datasets; not assessable under 3",
        "reading_L": "the lineage figure with its own interval",
        "reading_C": "class-balanced log loss, stratified bootstrap, lineage figure fixed",
        "reading_D": "BioProject cluster bootstrap of the scored side, lineage figure fixed",
        "reading_T": "pilot only: two-sample 95% interval of each difference (z 1.96)",
        "tie_rule": "a tie on the largest BioProject goes to the lower accession",
        "resample_order": "labelled genomes with a call, by assembly accession",
    }
    return [{"key": k, "value": v} for k, v in pairs.items()]


def plan_rows(lanes, inputs):
    rows, copies = [], {}
    for name, lane, rel, sums in PLANS:
        path, sums_path = os.path.join(lanes, rel), os.path.join(lanes, sums)
        for p in (path, sums_path):
            need(os.path.isfile(p), f"missing file {p}")
            inputs[p] = sha256_file(p)
        with open(sums_path, encoding="utf-8") as f:
            text = f.read()
        digest, fixed = SHA.search(text), TIME.search(text)
        need(digest is not None and fixed is not None, f"{sums_path}: no sha256 or time")
        need(digest.group(1) == inputs[path], f"{rel}: its sha256 differs from {sums}")
        with open(path, "rb") as f:
            copies["plan/" + name] = f.read()
        rows.append({"file": name, "lane": lane, "sha256": digest.group(1), "fixed": fixed.group(1), "time_from": sums})
    return rows, copies


def readme(files, rows_of, settings):
    s = {r["key"]: r["value"] for r in settings}
    table = "\n".join(f"| `{name}` | {rows_of.get(name, '')} | {DESCRIBE[name]} |" for name in sorted(files))
    return README.format(release=s["release"], table=table)


DESCRIBE = {
    "README.md": "this list",
    "species.tsv": "one row per species: each funnel step to the scored genomes, the overlap counts (each a lower bound), and the largest BioProject",
    "datasets.tsv": "one row per dataset and learner scored: labels, recall, specificity and log loss with intervals, the lineage figures, the plan's reading, the class mix",
    "readings.tsv": "one row per estimable dataset and learner: readings L, C and D (and T for the pilot), and the errors in the largest BioProject",
    "learners.tsv": "derived: one row per species and learner, under the plan's rule and each reading, with the margin",
    "totals.tsv": "derived: per learner, the counts over the eight species",
    "differences.tsv": "derived: each pair of learners, per species and over the species, under the four readings",
    "genomes.tsv": "the scored genomes as accessions only: no label and no prediction",
    "order.tsv": "per species, when predictions were written and when the label requests were sent; for the pilot, the requests before the plan with no recorded parameters",
    "plans.tsv": "the three plan files with their sha256 and the time each was fixed",
    "settings.tsv": "the release, the class line, the pass rule, and the seed and draws of every bootstrap",
    "inputs.sha256": "the sha256 of every input file the script read",
    "plan/PLAN.md": "the pilot's plan, as lane outside fixed it",
    "plan/ADDENDUM.md": "the eight species' addendum, as lane outside2 fixed it",
    "plan/ADDENDUM2.md": "the second addendum with readings L, C and D, as lane outside2 fixed it",
    "SHA256SUMS": "the sha256 of every other file of this directory",
}

README = """# The test on data from another public database: the record

This directory is the record behind [`docs/other-database-test.md`](../../docs/other-database-test.md). `scripts/other-database/record.py` writes every file of it from the saved files of lanes `outside` and `outside2`, the reports and interval files of release `{release}`, and the data root's drug name list. It fits nothing, runs no `grm` and downloads nothing. It recomputes every measure from the saved predictions and labels.

## Files

| File | Rows | Holds |
|---|---|---|
{table}

Every file is plain text. Each TSV has a header row and is sorted by its first columns. The script writes the same bytes on every run and at any hash seed.

## Rules

- **Labels.** An NCBI antibiotic name maps to a release drug by exact equality, then by the data root's `curated/drug_names.json`, then by the plan's fixed table. One label per BioSample and drug: R and S together is a conflict, then R, then S; intermediate alone gives no label.
- **Measures.** As `plan/PLAN.md` section 6 fixes them. A bootstrap draws the labelled genomes with a call in the order of their assembly accessions.
- **Readings.** L, C and D are those of `plan/ADDENDUM2.md` section B2, fixed in advance for the eight species and exploratory for the pilot. T is the pilot review's two-sample reading, exploratory, and the pilot's L, C and D are that review's readings too.
- **Ties.** A tie on the largest BioProject goes to the lower accession: letters first, then the lower number, so PRJNA717739 comes before PRJNA1160997. A field that holds several BioProjects, or a genome with none, comes after every single accession.
- **Words.** `meets` in `learners.tsv` and `totals.tsv` means the learner meets the plan's rule in the species: it holds on more than half of at least 3 estimable datasets. `margin` is how many datasets would have to change for the other result.
- **Pre-plan requests.** `unrecorded_before_plan` in `order.tsv` counts the pilot's requests before `PLAN.md` was fixed that the request log holds with no parameters.

## Commands

```sh
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --out reference/other-database
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --check reference/other-database
```

`LANES` holds the lanes' output directories `outside/` and `outside2/`, `R` is the unpacked release and `ROOT` is the public data root. `--check` exits 0 only when this directory is what the script writes from them.
"""


def build(lanes, release, root, jobs):
    """Every file of the record, as path -> bytes."""
    inputs = {}
    for p in (os.path.join(release, "index.json"), os.path.join(root, "curated", "drug_names.json")):
        need(os.path.isfile(p), f"missing file {p}")
        inputs[p] = sha256_file(p)
    with open(os.path.join(release, "index.json"), encoding="utf-8") as f:
        index = json.load(f)
    all_drugs = {d["antibiotic"] for d in index["datasets"]}
    drug_of = {d["id"]: d["antibiotic"] for d in index["datasets"]}
    with open(os.path.join(root, "curated", "drug_names.json"), encoding="utf-8") as f:
        rawmap = {n["raw"].strip().lower(): n["canonical"] for n in json.load(f)["names"] if n.get("kind") == "drug"}

    noacc_path = os.path.join(lanes, "outside2", "results", "no-accession-share.tsv")
    need(os.path.isfile(noacc_path), f"missing file {noacc_path}")
    inputs[noacc_path] = sha256_file(noacc_path)
    noacc = {r["species"]: r for r in read_tsv(noacc_path)}

    plans, copies = plan_rows(lanes, inputs)
    species_rows, tasks, genomes, species_of, scored_bs = [], [], [], {}, {}
    for slug, name in SPECIES:
        row, t, g, owners, bs = read_species(lanes, release, slug, all_drugs, rawmap, drug_of, inputs)
        need(name in noacc, f"{noacc_path}: no row for {name}")
        n = noacc[name]
        need(int(n["neither_biosample_nor_assembly"]) <= int(n["labelled_genomes"]), f"{noacc_path}: {name} does not add up")
        row.update({"bvbrc_labelled": n["labelled_genomes"], "bvbrc_no_accession": n["neither_biosample_nor_assembly"],
                    "no_accession_share": f3(int(n["neither_biosample_nor_assembly"]) / int(n["labelled_genomes"]))})
        species_rows.append(row)
        tasks += t
        genomes += g
        scored_bs[slug] = bs
        for b, sp in owners.items():
            if sp != PILOT:
                need(species_of.get(b, sp) == sp, f"BioSample {b} is selected in two species")
                species_of[b] = sp
    order = order_rows(lanes, species_of, scored_bs, plans, inputs)

    if jobs > 1:
        with multiprocessing.get_context("fork").Pool(jobs) as pool:
            results = pool.map(score_dataset, tasks, chunksize=1)
    else:
        results = [score_dataset(t) for t in tasks]
    datasets = [r for rows, _ in results for r in rows]
    readings = [r for _, rs in results for r in rs]

    settings = settings_rows(os.path.basename(os.path.normpath(release)), index["evaluation_revision"])
    learners, totals, differences = derive(datasets, readings)
    files = {
        "species.tsv": tsv("species.tsv", species_rows),
        "datasets.tsv": tsv("datasets.tsv", datasets),
        "readings.tsv": tsv("readings.tsv", readings),
        "learners.tsv": tsv("learners.tsv", learners),
        "totals.tsv": tsv("totals.tsv", totals),
        "differences.tsv": tsv("differences.tsv", differences),
        "genomes.tsv": tsv("genomes.tsv", genomes),
        "order.tsv": tsv("order.tsv", order),
        "plans.tsv": tsv("plans.tsv", plans),
        "settings.tsv": tsv("settings.tsv", settings),
    }
    roots = [(os.path.normpath(lanes), "lanes"), (os.path.normpath(release), "release"), (os.path.normpath(root), "data-root")]

    def named(path):
        path = os.path.normpath(path)
        for base, label in roots:
            if path.startswith(base + os.sep):
                return label + "/" + os.path.relpath(path, base)
        raise Stop(f"input {path} is under no named root")

    files["inputs.sha256"] = "".join(f"{d}  {p}\n" for p, d in sorted((named(p), d) for p, d in inputs.items()))
    out = {k: v.encode("utf-8") for k, v in files.items()}
    out.update(copies)
    rows_of = {k: str(len(parse(v))) for k, v in files.items() if k.endswith(".tsv")}
    rows_of["inputs.sha256"] = str(len(inputs))
    names = sorted(list(out) + ["README.md", "SHA256SUMS"])
    out["README.md"] = readme(names, rows_of, settings).encode("utf-8")
    out["SHA256SUMS"] = "".join(f"{hashlib.sha256(out[k]).hexdigest()}  {k}\n" for k in sorted(out)).encode("utf-8")
    return out


def on_disk(directory):
    found = {}
    for base, _, names in os.walk(directory):
        for n in names:
            p = os.path.join(base, n)
            with open(p, "rb") as f:
                found[os.path.relpath(p, directory).replace(os.sep, "/")] = f.read()
    return found


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--lanes", required=True)
    ap.add_argument("--release", required=True)
    ap.add_argument("--data-root", required=True)
    ap.add_argument("--jobs", type=int, default=1)
    what = ap.add_mutually_exclusive_group(required=True)
    what.add_argument("--out")
    what.add_argument("--check")
    a = ap.parse_args(argv)
    try:
        files = build(a.lanes, a.release, a.data_root, max(1, a.jobs))
    except Stop as e:
        print(f"record.py: stopped: {e}", file=sys.stderr)
        return 2
    if a.check:
        found = on_disk(a.check) if os.path.isdir(a.check) else {}
        bad = sorted(k for k in set(files) | set(found) if files.get(k) != found.get(k))
        for k in bad:
            state = "missing" if k not in found else ("not written by the script" if k not in files else "differs")
            print(f"record.py --check: {k}: {state}", file=sys.stderr)
        return 1 if bad else 0
    for k, v in files.items():
        p = os.path.join(a.out, k)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "wb") as f:
            f.write(v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
