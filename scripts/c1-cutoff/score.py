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
scores exactly the cohort; the five learners in order.

Scores, per result, split (lineage, and study when the cohort has one),
truth set (all scored genomes, new lineages) and learner: the estimate (the
split's pooled balanced accuracy) and its interval (the range of its five
fold balanced accuracies), the truth's balanced accuracy, the error (their
absolute difference) and whether the interval covers the truth; recall,
specificity, the very major and the major error rate of the estimate and
of the truth, each with its count, denominator and Wilson 95% interval, the
error rates also with their Clopper-Pearson 95% interval, and whether the
estimate's interval covers the truth's value; the log loss of each with
its denominator. A figure that cannot be computed is empty.

Writes <out dir>/checks.tsv, <out dir>/scores.tsv (one row per result,
split, truth set and learner), <out dir>/lines.tsv (per dataset, split,
truth set, learner and line: draws, draws available, mean error and
coverage) and <out dir>/count.tsv (the count of the plan's C1 row for the
datasets given). Prints the check counts only. Exits 1 when a check fails.

The second form writes count.tsv from several species' lines.tsv files.

The count: for each line against 25 (and the twin, r 25 at seed 45,
against 25), per split, truth set and learner, a dataset enters when at
least 3 of its draws are available at both; its value at a line is the
mean error over that line's available draws; better, worse or same by more
than MARGIN; the median of the differences (line less 25) and its species
interval: species drawn with replacement RESAMPLES times, the median of
their datasets' differences each time, read at the 101st and 3,900th
values. A line is a candidate when that interval lies wholly below MARGIN
and its coverage is within 0.05 of 25's. Outputs are aggregates; no genome
id is read from a result, none is written.
"""

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

FORMAT = "grm-cohort-draw-result-v1"
REVISION = 14
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
# The study estimate's scheme as grm learn --exp-cohort writes it (src/cutoff.rs
# STUDY_SCHEME); its rows here are named split "study".
STUDY_SCHEME = "study-whole"
SPLITS = {"lineage": "lineage", "study": STUDY_SCHEME}
TRUTHS = ["all", "new_lineage"]
Z = 1.959963984540054
MARGIN = 0.02
COVERAGE_BAND = 0.05
MIN_DRAWS = 3
RESAMPLES = 4000
LO_INDEX, HI_INDEX = 100, 3899
BOOT_SEED = 20261122
BASE_LINE = "25"
TWIN_SEED = 45


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
                tba = balanced(tc) if tc else None
                row.update({"estimate": ba, "interval_lo": lo, "interval_hi": hi, "truth": tba,
                            "truth_null_fit": tr["null_fit"] if tr else None})
                avail = ba is not None and tba is not None
                row["available"] = avail
                row["error"] = abs(ba - tba) if avail else None
                row["covers"] = (lo <= tba <= hi) if avail and lo is not None else None
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
                    row[f"{name}_covers"] = (lo_ <= v <= hi_) if None not in (lo_, hi_, v) else None
                rows.append(row)
    return rows


def line_name(r, seed):
    return f"{r}-s{TWIN_SEED}" if int(seed) == TWIN_SEED else str(r)


def lines(rows, expected):
    """Per dataset, split, truth set, learner and line: draws, available draws, mean error, coverage."""
    draws = {}
    for e in expected:
        draws.setdefault((e["dataset"], line_name(e["r"], e["seed"])), set()).add(int(e["draw"]))
    cells = {}
    for row in rows:
        k = (row["dataset"], row["species_slug"], row["split"], row["truth_set"], row["learner"], line_name(row["r"], row["seed"]))
        cells.setdefault(k, []).append(row)
    out = []
    for k in sorted(cells):
        rs = cells[k]
        av = [r for r in rs if r["available"]]
        cov = [r["covers"] for r in av if r["covers"] is not None]
        out.append({
            "dataset": k[0], "species_slug": k[1], "split": k[2], "truth_set": k[3], "learner": k[4], "line": k[5],
            "draws": len(draws.get((k[0], k[5]), ())), "available": len(av),
            "available_draws": ",".join(str(r["draw"]) for r in sorted(av, key=lambda r: int(r["draw"]))),
            "mean_error": (sum(r["error"] for r in av) / len(av)) if av else None,
            "covered": sum(cov), "coverage_n": len(cov),
        })
    return out


def count(line_rows):
    """The count of the plan's C1 row, per split, truth set, learner and line against 25."""
    cells = {}
    for r in line_rows:
        cells[(r["split"], r["truth_set"], r["learner"], r["line"], r["dataset"])] = r
    groups = sorted({(s, t, l) for s, t, l, _, _ in cells})
    names = sorted({k[3] for k in cells if k[3] != BASE_LINE}, key=lambda x: (int(x.split("-")[0]), x))
    rng = np.random.default_rng(BOOT_SEED)
    out = []
    for s, t, l in groups:
        for line in names:
            pairs = []
            for (s2, t2, l2, ln, d), row in cells.items():
                if (s2, t2, l2, ln) != (s, t, l, line):
                    continue
                base = cells.get((s, t, l, BASE_LINE, d))
                if base is None:
                    continue
                both = set(filter(None, str(row["available_draws"]).split(","))) & set(filter(None, str(base["available_draws"]).split(",")))
                pairs.append((d, row["species_slug"], row, base, len(both) >= MIN_DRAWS))
            entered = [p for p in pairs if p[4]]
            diffs = [float(p[2]["mean_error"]) - float(p[3]["mean_error"]) for p in entered]
            better = sum(1 for x in diffs if x < -MARGIN)
            worse = sum(1 for x in diffs if x > MARGIN)
            res = {"split": s, "truth_set": t, "learner": l, "line": line, "datasets": len(pairs), "entered": len(entered),
                   "not_scored": len(pairs) - len(entered), "better": better, "worse": worse, "same": len(entered) - better - worse,
                   "median_difference": float(np.median(diffs)) if diffs else None, "species": len({p[1] for p in entered}),
                   "interval_lo": None, "interval_hi": None}
            cov = sum(int(p[2]["covered"]) for p in entered), sum(int(p[2]["coverage_n"]) for p in entered)
            bcov = sum(int(p[3]["covered"]) for p in entered), sum(int(p[3]["coverage_n"]) for p in entered)
            res.update({"covered": cov[0], "coverage_n": cov[1], "base_covered": bcov[0], "base_coverage_n": bcov[1]})
            if entered:
                by_sp = {}
                for p, x in zip(entered, diffs):
                    by_sp.setdefault(p[1], []).append(x)
                sps = sorted(by_sp)
                meds = []
                for _ in range(RESAMPLES):
                    pick = rng.integers(0, len(sps), len(sps))
                    meds.append(float(np.median(np.concatenate([by_sp[sps[i]] for i in pick]))))
                meds.sort()
                res["interval_lo"], res["interval_hi"] = meds[LO_INDEX], meds[HI_INDEX]
            c_r = cov[0] / cov[1] if cov[1] else None
            c_b = bcov[0] / bcov[1] if bcov[1] else None
            res["candidate"] = (res["interval_hi"] is not None and res["interval_hi"] < MARGIN and c_r is not None
                                and c_b is not None and abs(c_r - c_b) <= COVERAGE_BAND)
            out.append(res)
    return out


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
        write_tsv(out / "count.tsv", count(rows))
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
    write_tsv(out / "count.tsv", count(lr))
    print(f"checks: {len(expected) - failed} of {len(expected)} results pass every check; "
          f"{sum(1 for c in checks if not c['ok'])} failed checks of {len(checks)}")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
