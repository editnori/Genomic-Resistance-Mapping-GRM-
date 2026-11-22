"""The scorer of plan block C1 part B (docs/cohort-draw.md, measures and files).

  python3 -I score.py <expect.tsv> <cohort dir> <results dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> <lines.tsv> [<lines.tsv> ...]

The first form checks every result expect.py lists, then scores it. Checks,
per result: it exists; its format and revision; its fingerprint against the
release report's; its settings against the report's but the seed, which
must be the expected one; its dataset, r and draw; its cohort file's sha256
against the expected one and the file on disk; that it was not refused and
its cohort holds r of each class; its dataset counts against the report's;
that the scored genomes, the copies left out and the cohort sum to the
dataset, by class; that the new lineages lie inside the scored genomes;
that each truth scores exactly the scored genomes (all) and the new
lineages, so no truth comes from the cohort's own folds; that each estimate
scores exactly the cohort, with the balanced accuracy of its own calls; the
five learners in order.

Scores, per result, split (lineage, random, and study when the cohort has
one), truth set (all scored genomes, new lineages) and learner: the
estimate (the split's pooled balanced accuracy), the truth's balanced
accuracy, the error (their absolute difference), and two intervals on the
estimate, each with its width and whether it covers the truth, an edge
counting as covered:
- the fold range, the lowest to the highest of the five fold balanced
  accuracies, which is what a report shows;
- the levelled interval, from the split's pooled calls: recall p_R with its
  Wilson 95% interval (l_R, u_R), specificity p_S with (l_S, u_S),
  b = (p_R + p_S) / 2, lower b - sqrt((p_R - l_R)^2 + (p_S - l_S)^2) / 2,
  upper b + sqrt((u_R - p_R)^2 + (u_S - p_S)^2) / 2.
Beside them: recall, specificity, the very major and the major error rate
of the estimate and of the truth, each with its count, denominator and
Wilson 95% interval, the error rates also with their Clopper-Pearson 95%
interval, and whether the estimate's interval covers the truth's value; the
log loss of each with its denominator. A figure that cannot be computed is
empty.

Writes <out dir>/checks.tsv, <out dir>/scores.tsv (one row per result,
split, truth set and learner), <out dir>/lines.tsv (per dataset, split,
truth set, learner and line: the draws, the available draws with each
one's error, coverage and width, and over the available draws the mean
error, both coverages and both median widths) and <out dir>/count.tsv (the
count of the plan's C1 row for the datasets given). Prints the check counts
only. Exits 1 when a check fails.

The second form writes count.tsv from several species' lines.tsv files.

The count, for each line against 25 (and the twin, r 25 at seed 45, against
25), per split, truth set and learner. Draw d at a line nests inside draw d
at 25, so the lines are paired by draw. A dataset enters when at least 3 of
its draws are available at both lines; its difference is the mean over
those shared draws of the error at the line less the error at 25. A dataset
that enters on fewer, or has no line at 25, is counted as not scored.
Better, worse or same by more than MARGIN; the median of the differences
and its species interval: species drawn with replacement RESAMPLES times,
the median of their datasets' differences each time, read at the 101st and
3,900th values, each row's resamples seeded from the sha256 of
split|truth set|learner|line, so a row's interval does not move when rows
are added. Coverage at both lines is counted over the same shared results;
each line's own available draws and its coverage over all its own results
are written beside. Each row is written twice, once with the levelled
interval's coverage and once with the fold range's.

A learner meets the release row when its species interval lies wholly
below MARGIN and its coverage is within 0.05 of 25's. A line is a candidate
when at least 3 of the 5 learners meet it. One row per line decides (the
lead's ruling of 22 Nov 2026): the lineage split, the truth on all scored
genomes, the levelled interval, 3 of 5 learners. Every other row is marked
as deciding nothing: the truth on new lineages, the random and the study
split, each learner alone, the fold range's coverage and the twin. Outputs
are aggregates; no genome id is read from a result, none is written.
"""

import csv
import hashlib
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

FORMAT = "grm-cohort-draw-result-v1"
REVISION = 14
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
# The study estimate's scheme as grm learn --exp-cohort writes it (src/cutoff.rs
# STUDY_SCHEME); its rows here are named split "study".
STUDY_SCHEME = "study-whole"
SPLITS = {"lineage": "lineage", "random": "random", "study": STUDY_SCHEME}
TRUTHS = ["all", "new_lineage"]
Z = 1.959963984540054
MARGIN = 0.02
COVERAGE_BAND = Fraction(1, 20)  # 0.05, compared exactly
LEARNERS_TO_PASS = 3
MIN_DRAWS = 3
RESAMPLES = 4000
LO_INDEX, HI_INDEX = 100, 3899
BASE_LINE = "25"
TWIN_SEED = 45
# The one row that decides (the lead's ruling D2 of 22 Nov 2026).
DECIDING = {"split": "lineage", "truth_set": "all", "coverage": "levelled"}
COVERAGES = ("levelled", "fold_range")
LINE_ROW = f"{LEARNERS_TO_PASS} of {len(LEARNERS)}"


# ---------- intervals ----------

def wilson(k, n):
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def _betacf(a, b, x):
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1, a - 1
    c, d = 1.0, 1 - qab * x / qap
    d = 1 / (d if abs(d) > tiny else tiny)
    h = d
    for m in range(1, 400):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1 + aa * d
        d = 1 / (d if abs(d) > tiny else tiny)
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1 + aa * d
        d = 1 / (d if abs(d) > tiny else tiny)
        c = 1 + aa / c
        c = c if abs(c) > tiny else tiny
        de = d * c
        h *= de
        if abs(de - 1) < 1e-15:
            break
    return h


def betainc(a, b, x):
    """The regularized incomplete beta function I_x(a, b)."""
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lb = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b) + a * math.log(x) + b * math.log1p(-x)
    if x < (a + 1) / (a + b + 2):
        return math.exp(lb) * _betacf(a, b, x) / a
    return 1 - math.exp(lb) * _betacf(b, a, 1 - x) / b


def beta_quantile(q, a, b):
    lo, hi = 0.0, 1.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if betainc(a, b, mid) < q:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def clopper_pearson(k, n):
    if n == 0:
        return (None, None)
    lo = 0.0 if k == 0 else beta_quantile(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta_quantile(0.975, k + 1, n - k)
    return (lo, hi)


# ---------- figures from calls ----------

def rates(c):
    """(name, k, n) of recall, specificity, very major and major error."""
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    return [("recall", c["tp"], pos), ("specificity", c["tn"], neg), ("very_major", c["fn"], pos), ("major", c["fp"], neg)]


def balanced(c):
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    if pos == 0 or neg == 0:
        return None
    return (c["tp"] / pos + c["tn"] / neg) / 2


def levelled(c):
    """(b, lower, upper) of the levelled interval from pooled calls (the lead's ruling D1).

    b is the balanced accuracy of the calls. Each side adds the two rates'
    Wilson 95% half-widths on that side in quadrature and halves the sum,
    since b is the mean of the two rates. (None, None, None) when a class is
    empty.
    """
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    if pos == 0 or neg == 0:
        return (None, None, None)
    pr, ps = c["tp"] / pos, c["tn"] / neg
    lr, ur = wilson(c["tp"], pos)
    ls, us = wilson(c["tn"], neg)
    b = (pr + ps) / 2
    return (b, b - math.sqrt((pr - lr) ** 2 + (ps - ls) ** 2) / 2, b + math.sqrt((ur - pr) ** 2 + (us - ps) ** 2) / 2)


def covers(lo, hi, v):
    """Whether [lo, hi] covers v, an edge counting as covered; None when any is missing."""
    if lo is None or hi is None or v is None:
        return None
    return lo <= v <= hi


def fmt(x):
    if x is None:
        return ""
    if isinstance(x, bool):
        return str(int(x))
    if isinstance(x, float):
        return "" if math.isnan(x) else repr(round(x, 12))
    return str(x)


# ---------- checks ----------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check(exp, doc, report, cohort_sha):
    """[(name, ok, detail)] for one result."""
    out = []
    r = int(exp["r"])

    def add(name, ok, detail=""):
        out.append((name, bool(ok), detail))

    add("format and revision", doc.get("format") == FORMAT and doc.get("evaluation_revision") == REVISION)
    add("fingerprint is the release report's", doc.get("release_fingerprint") == report["input_fingerprint"])
    cfg, rcfg = dict(doc.get("eval_config", {})), dict(report["eval_config"])
    seed = cfg.pop("seed", None)
    rcfg.pop("seed", None)
    add("settings are the report's but the seed", cfg == rcfg)
    add("seed is the expected one", seed == int(exp["seed"]) == doc.get("seed"))
    add("dataset, r and draw", doc.get("id") == exp["dataset"] and doc.get("r") == r and doc.get("draw") == int(exp["draw"]))
    add("cohort file sha256", doc.get("cohort_file_sha256") == exp["cohort_sha256"] == cohort_sha)
    coh = doc.get("cohort") or {}
    add("not refused, r of each class", doc.get("refused") is None and coh.get("samples") == 2 * r and coh.get("resistant") == r,
        doc.get("refused") or "")
    # src/learn/eval.rs support_floor of the cohort's 2r genomes, not of the dataset's.
    add("the cohort's support floor is its own", coh.get("min_support") == max(2 * r // 200, 3))
    ds, rc = doc.get("dataset", {}), report["counts"]
    add("dataset counts are the report's", ds.get("samples") == rc["samples"] and ds.get("resistant") == rc["resistant"])
    sc, cp, nl = doc.get("scored", {}), doc.get("copies_left_out", {}), doc.get("new_lineage", {})
    add("scored, copies and cohort sum to the dataset",
        sc.get("samples", -1) + cp.get("samples", -1) + 2 * r == ds.get("samples")
        and sc.get("resistant", -1) + cp.get("resistant", -1) + r == ds.get("resistant"))
    add("new lineages inside the scored genomes",
        0 <= nl.get("samples", -1) <= sc.get("samples", -1) and 0 <= nl.get("resistant", -1) <= sc.get("resistant", -1)
        and nl.get("samples", 0) - nl.get("resistant", 0) <= sc.get("samples", 0) - sc.get("resistant", 0))
    truth = doc.get("truth") or {}
    ok = True
    if truth.get("unavailable") is None:
        rows = truth.get("learners", [])
        ok = [t["learner"] for t in rows] == LEARNERS
        for t in rows:
            a, n = t["all"], t["new_lineage"]
            ok = ok and a["n"] == sc.get("samples") and a["tp"] + a["fn"] == sc.get("resistant")
            ok = ok and n["n"] == nl.get("samples") and n["tp"] + n["fn"] == nl.get("resistant")
    add("each truth scores the scored genomes and new lineages", ok)
    ok = True
    schemes = [e["scheme"] for e in doc.get("estimates", [])]
    ok = schemes[:2] == ["random", "lineage"] and (schemes[2:] == ([STUDY_SCHEME] if exp["study"] == "1" else []))
    for e in doc.get("estimates", []):
        if e.get("unavailable") is None:
            ok = ok and [x["learner"] for x in e["learners"]] == LEARNERS
            for x in e["learners"]:
                c = x["calls"]
                ok = ok and c["n"] == 2 * r and c["tp"] + c["fn"] == r and len(x["fold_balanced_accuracy"]) == 5
    add("each estimate scores the cohort", ok)
    # The levelled interval is built around the balanced accuracy of the
    # pooled calls; grm writes the estimate from the same scores at the same
    # cutoff (binary_metrics at 0.5, cutoff::calls at CALL_CUTOFF), so the two
    # are one figure.
    ok = True
    for e in doc.get("estimates", []):
        if e.get("unavailable") is None:
            for x in e.get("learners", []):
                ba, b = x.get("balanced_accuracy"), balanced(x["calls"])
                ok = ok and isinstance(ba, (int, float)) and b is not None and abs(ba - b) <= 1e-9
    add("each estimate is the balanced accuracy of its calls", ok)
    return out


# ---------- scores ----------

def score_rows(exp, doc):
    """One row per split, truth set and learner of a result that passed its checks."""
    rows = []
    est = {e["scheme"]: e for e in doc["estimates"]}
    truth = doc.get("truth") or {}
    t_ok = truth.get("unavailable") is None
    t_rows = {t["learner"]: t for t in truth.get("learners", [])}
    for split, scheme in SPLITS.items():
        if scheme not in est:
            continue
        e = est[scheme]
        e_ok = e.get("unavailable") is None
        e_rows = {x["learner"]: x for x in e.get("learners", [])}
        for ts in TRUTHS:
            for learner in LEARNERS:
                row = {
                    "dataset": exp["dataset"], "species_slug": exp["species_slug"], "antibiotic": doc.get("antibiotic", ""),
                    "r": exp["r"], "draw": exp["draw"], "seed": exp["seed"], "split": split, "truth_set": ts,
                    "learner": learner, "estimate_unavailable": "" if e_ok else (e.get("unavailable") or ""),
                    "truth_unavailable": "" if t_ok else (truth.get("unavailable") or ""),
                }
                ex = e_rows.get(learner) if e_ok else None
                tr = t_rows.get(learner) if t_ok else None
                tc = tr[ts] if tr else None
                ba = ex["balanced_accuracy"] if ex else None
                ba = None if ba is None or (isinstance(ba, float) and math.isnan(ba)) else ba
                folds = [f for f in (ex["fold_balanced_accuracy"] if ex else []) if f is not None and not math.isnan(f)]
                lo, hi = (min(folds), max(folds)) if ex and len(folds) == 5 else (None, None)
                b, llo, lhi = levelled(ex["calls"]) if ex else (None, None, None)
                tba = balanced(tc) if tc else None
                avail = ba is not None and tba is not None
                row.update({"estimate": ba, "estimate_from_calls": b, "truth": tba,
                            "truth_null_fit": tr["null_fit"] if tr else None, "available": avail,
                            "error": abs(ba - tba) if avail else None,
                            "levelled_lo": llo, "levelled_hi": lhi,
                            "levelled_width": (lhi - llo) if llo is not None else None,
                            "levelled_covers": covers(llo, lhi, tba) if avail else None,
                            "fold_range_lo": lo, "fold_range_hi": hi,
                            "fold_range_width": (hi - lo) if lo is not None else None,
                            "fold_range_covers": covers(lo, hi, tba) if avail else None})
                for who, c in (("estimate", ex["calls"] if ex else None), ("truth", tc)):
                    if c is None:
                        for name in ("recall", "specificity", "very_major", "major"):
                            for s in ("k", "n", "value", "wilson_lo", "wilson_hi"):
                                row[f"{who}_{name}_{s}"] = None
                        for name in ("very_major", "major"):
                            row[f"{who}_{name}_cp_lo"] = row[f"{who}_{name}_cp_hi"] = None
                        row[f"{who}_log_loss"] = row[f"{who}_log_loss_n"] = None
                        continue
                    for name, k, n in rates(c):
                        wl, wh = wilson(k, n)
                        row.update({f"{who}_{name}_k": k, f"{who}_{name}_n": n, f"{who}_{name}_value": (k / n if n else None),
                                    f"{who}_{name}_wilson_lo": wl, f"{who}_{name}_wilson_hi": wh})
                        if name in ("very_major", "major"):
                            cl, ch = clopper_pearson(k, n)
                            row[f"{who}_{name}_cp_lo"], row[f"{who}_{name}_cp_hi"] = cl, ch
                    row[f"{who}_log_loss"] = c["loss_sum"] / c["n"] if c["n"] else None
                    row[f"{who}_log_loss_n"] = c["n"]
                for name in ("recall", "specificity", "very_major", "major"):
                    lo_, hi_, v = row[f"estimate_{name}_wilson_lo"], row[f"estimate_{name}_wilson_hi"], row[f"truth_{name}_value"]
                    row[f"{name}_covers"] = covers(lo_, hi_, v)
                rows.append(row)
    return rows


def line_name(r, seed):
    return f"{r}-s{TWIN_SEED}" if int(seed) == TWIN_SEED else str(r)


def lines(rows, expected):
    """Per dataset, split, truth set, learner and line: the draws, each available draw's figures, and their summary.

    `available_draws`, `errors`, `levelled_covers`, `levelled_widths`,
    `fold_range_covers` and `fold_range_widths` are comma lists in draw
    order, one entry per available draw (a fold range that cannot be taken
    leaves its entries empty). The count pairs lines on them. Every figure
    is written as fmt writes it, so the count reads the same values from
    memory and from the file.
    """
    draws = {}
    for e in expected:
        draws.setdefault((e["dataset"], line_name(e["r"], e["seed"])), set()).add(int(e["draw"]))
    cells = {}
    for row in rows:
        k = (row["dataset"], row["species_slug"], row["split"], row["truth_set"], row["learner"], line_name(row["r"], row["seed"]))
        cells.setdefault(k, []).append(row)
    out = []
    for k in sorted(cells):
        av = sorted((r for r in cells[k] if r["available"]), key=lambda r: int(r["draw"]))
        res = {"dataset": k[0], "species_slug": k[1], "split": k[2], "truth_set": k[3], "learner": k[4], "line": k[5],
               "draws": len(draws.get((k[0], k[5]), ())), "available": len(av),
               "available_draws": ",".join(str(int(r["draw"])) for r in av),
               "errors": ",".join(fmt(r["error"]) for r in av),
               "mean_error": (sum(r["error"] for r in av) / len(av)) if av else None}
        for kind in COVERAGES:
            cov = [r[f"{kind}_covers"] for r in av if r[f"{kind}_covers"] is not None]
            widths = [r[f"{kind}_width"] for r in av if r[f"{kind}_width"] is not None]
            res.update({f"{kind}_covers": ",".join(fmt(r[f"{kind}_covers"]) for r in av),
                        f"{kind}_widths": ",".join(fmt(r[f"{kind}_width"]) for r in av),
                        f"{kind}_covered": sum(cov), f"{kind}_coverage_n": len(cov),
                        f"{kind}_median_width": float(np.median(widths)) if widths else None})
        out.append(res)
    return out


def per_draw(row):
    """{draw: {"error", "levelled_covers", ...}} of one lines.tsv row, read from its comma lists."""
    def split(name):
        text = str(row.get(name) or "")
        return text.split(",") if text else []

    draws = [int(d) for d in split("available_draws")]
    cols = {name: split(name) for name in ("errors",) + tuple(f"{k}_{w}" for k in COVERAGES for w in ("covers", "widths"))}
    for name, values in cols.items():
        if len(values) != len(draws):
            raise ValueError(f"lines row {row.get('dataset')} {row.get('line')}: {name} has {len(values)} entries "
                             f"for {len(draws)} available draws")
    num = lambda v: float(v) if v != "" else None
    return {d: {"error": float(cols["errors"][i]),
                **{f"{k}_covers": (bool(int(cols[f"{k}_covers"][i])) if cols[f"{k}_covers"][i] != "" else None) for k in COVERAGES},
                **{f"{k}_width": num(cols[f"{k}_widths"][i]) for k in COVERAGES}}
            for i, d in enumerate(draws)}


def row_seed(split, truth_set, learner, line):
    """The seed of one count row's resamples: the sha256 of its key, so no other row moves it."""
    return int(hashlib.sha256(f"{split}|{truth_set}|{learner}|{line}".encode()).hexdigest(), 16)


def coverage(results, kind):
    """(covered, n, median width) of one interval kind over per-draw results."""
    cov = [x[f"{kind}_covers"] for x in results if x[f"{kind}_covers"] is not None]
    widths = [x[f"{kind}_width"] for x in results if x[f"{kind}_width"] is not None]
    return sum(cov), len(cov), (float(np.median(widths)) if widths else None)


def within_band(covered, n, base_covered, base_n):
    """Coverage within COVERAGE_BAND of 25's, compared exactly; False when either has no result."""
    if not n or not base_n:
        return False
    return abs(Fraction(covered, n) - Fraction(base_covered, base_n)) <= COVERAGE_BAND


def count(line_rows):
    """The count of the plan's C1 row: per split, truth set, learner, line against 25 and coverage, then the line rows."""
    cells = {}
    for r in line_rows:
        cells[(r["split"], r["truth_set"], r["learner"], r["line"], r["dataset"])] = r
    groups = sorted({(s, t, l) for s, t, l, _, _ in cells})
    names = sorted({k[3] for k in cells if k[3] != BASE_LINE}, key=lambda x: (int(x.split("-")[0]), x))
    learner_rows = []
    for s, t, l in groups:
        for line in names:
            pairs = []  # (dataset, species, line row, base row or None, shared draws)
            for (s2, t2, l2, ln, d), row in sorted(cells.items()):
                if (s2, t2, l2, ln) != (s, t, l, line):
                    continue
                base = cells.get((s, t, l, BASE_LINE, d))
                own = per_draw(row)
                bown = per_draw(base) if base is not None else {}
                shared = sorted(set(own) & set(bown))
                pairs.append((d, row["species_slug"], own, bown, base is not None, shared))
            entered = [p for p in pairs if len(p[5]) >= MIN_DRAWS]
            # Paired on draws (the lead's ruling D3): the mean over the shared
            # draws of the error at the line less the error at 25, rounded as
            # the files write figures so the margin reads no float noise.
            diffs = [round(sum(p[2][d]["error"] - p[3][d]["error"] for d in p[5]) / len(p[5]), 12) for p in entered]
            better = sum(1 for x in diffs if x < -MARGIN)
            worse = sum(1 for x in diffs if x > MARGIN)
            base = {"row": "learner", "split": s, "truth_set": t, "learner": l, "line": line,
                    "datasets": len(pairs), "no_base_line": sum(1 for p in pairs if not p[4]), "entered": len(entered),
                    "not_scored": len(pairs) - len(entered), "better": better, "worse": worse,
                    "same": len(entered) - better - worse,
                    "median_difference": round(float(np.median(diffs)), 12) if diffs else None,
                    "species": len({p[1] for p in entered}), "interval_lo": None, "interval_hi": None,
                    "shared_draws": sum(len(p[5]) for p in entered),
                    "line_available": sum(len(p[2]) for p in pairs),
                    "base_available": sum(len(p[3]) for p in pairs)}
            if entered:
                by_sp = {}
                for p, x in zip(entered, diffs):
                    by_sp.setdefault(p[1], []).append(x)
                sps = sorted(by_sp)
                rng = np.random.default_rng(row_seed(s, t, l, line))
                meds = []
                for _ in range(RESAMPLES):
                    pick = rng.integers(0, len(sps), len(sps))
                    # Rounded like the differences: a median of two differences
                    # whose mean is 0.02 can come out a hair below it in floating
                    # point (-0.108 and 0.148 give 0.019999999999999997), and
                    # "wholly below 0.02" must not read that noise.
                    meds.append(round(float(np.median(np.concatenate([by_sp[sps[i]] for i in pick]))), 12))
                meds.sort()
                base["interval_lo"], base["interval_hi"] = meds[LO_INDEX], meds[HI_INDEX]
            for kind in COVERAGES:
                res = dict(base, coverage=kind)
                # Over the same shared results at both lines.
                res["covered"], res["coverage_n"], res["median_width"] = coverage(
                    [p[2][d] for p in entered for d in p[5]], kind)
                res["base_covered"], res["base_coverage_n"], res["base_median_width"] = coverage(
                    [p[3][d] for p in entered for d in p[5]], kind)
                # Beside it: each line over all its own results.
                res["line_covered"], res["line_coverage_n"], _ = coverage([x for p in pairs for x in p[2].values()], kind)
                res["base_covered_own"], res["base_coverage_n_own"], _ = coverage(
                    [x for p in pairs for x in p[3].values()], kind)
                res["meets"] = (res["interval_hi"] is not None and res["interval_hi"] < MARGIN
                                and within_band(res["covered"], res["coverage_n"], res["base_covered"], res["base_coverage_n"]))
                res["decides"] = False
                learner_rows.append(res)
    out = []
    for s, t in sorted({(s, t) for s, t, _ in groups}):
        for line in names:
            for kind in COVERAGES:
                rows = [r for r in learner_rows if (r["split"], r["truth_set"], r["line"], r["coverage"]) == (s, t, line, kind)]
                meeting = [l for l in LEARNERS if any(r["learner"] == l and r["meets"] for r in rows)]
                out.extend(sorted(rows, key=lambda r: LEARNERS.index(r["learner"]) if r["learner"] in LEARNERS else len(LEARNERS)))
                out.append({"row": "line", "split": s, "truth_set": t, "learner": LINE_ROW, "line": line, "coverage": kind,
                            "learners_meeting": len(meeting), "meeting": ",".join(meeting),
                            "candidate": len(meeting) >= LEARNERS_TO_PASS,
                            "decides": ({"split": s, "truth_set": t, "coverage": kind} == DECIDING
                                        and not line.endswith(f"-s{TWIN_SEED}"))})
    return out


COUNT_COLUMNS = ["row", "split", "truth_set", "learner", "line", "coverage", "decides", "candidate", "learners_meeting",
                 "meeting", "meets", "datasets", "no_base_line", "entered", "not_scored", "better", "worse", "same",
                 "median_difference", "species", "interval_lo", "interval_hi", "shared_draws", "covered", "coverage_n",
                 "base_covered", "base_coverage_n", "median_width", "base_median_width", "line_available",
                 "line_covered", "line_coverage_n", "base_available", "base_covered_own", "base_coverage_n_own"]


def write_tsv(path, rows, columns=None):
    columns = columns or (list(rows[0].keys()) if rows else [])
    with open(path, "w", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(c)) for c in columns) + "\n")


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            print(__doc__.strip().splitlines()[3], file=sys.stderr)
            return 2
        out = Path(argv[1])
        out.mkdir(parents=True, exist_ok=True)
        rows = [r for p in argv[2:] for r in read_tsv(p)]
        write_tsv(out / "count.tsv", count(rows), COUNT_COLUMNS)
        print(f"count: {len(rows)} line rows from {len(argv) - 2} files")
        return 0
    if len(argv) != 5:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    expected = read_tsv(argv[0])
    cohorts, results, release, out = Path(argv[1]), Path(argv[2]), Path(argv[3]), Path(argv[4])
    out.mkdir(parents=True, exist_ok=True)
    reports, checks, rows = {}, [], []
    failed = 0
    for e in expected:
        path = results / e["file"]
        if e["dataset"] not in reports:
            with open(release / "datasets" / f"{e['dataset']}.json") as f:
                reports[e["dataset"]] = json.load(f)
        if not path.is_file():
            checks.append({"file": e["file"], "check": "result exists", "ok": False, "detail": ""})
            failed += 1
            continue
        with open(path) as f:
            doc = json.load(f)
        coh = cohorts / e["cohort"]
        cs = sha256_file(coh) if coh.is_file() else None
        res = [("result exists", True, "")] + check(e, doc, reports[e["dataset"]], cs)
        for name, ok, detail in res:
            checks.append({"file": e["file"], "check": name, "ok": ok, "detail": detail})
        if all(ok for _, ok, _ in res):
            rows.extend(score_rows(e, doc))
        else:
            failed += 1
    write_tsv(out / "checks.tsv", checks, ["file", "check", "ok", "detail"])
    write_tsv(out / "scores.tsv", rows)
    lr = lines(rows, expected)
    write_tsv(out / "lines.tsv", lr)
    write_tsv(out / "count.tsv", count(lr), COUNT_COLUMNS)
    print(f"checks: {len(expected) - failed} of {len(expected)} results pass every check; "
          f"{sum(1 for c in checks if not c['ok'])} failed checks of {len(checks)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
