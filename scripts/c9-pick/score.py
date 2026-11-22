"""The scorer of plan block C9 (docs/served-pick.md, measures and files).

  python3 -I score.py <expect.tsv> <picks dir> <release dir> <out dir> [options]
  python3 -I score.py --count <out dir> <datasets.tsv> [<datasets.tsv> ...]

Options of the first form, each one of the design's open questions:
  --tie later|earlier|default   where an inner tie goes (default later);
                                default sends it to the default arm's
                                learner when that learner is tied
  --no-pick default|drop        a fold with no inner pick serves the
                                default's learner (default), or the draw
                                is left unscored
  --default per-draw|seed42     the default's learner is each draw's own
                                pick (default), or the seed-42 pick at
                                every draw

The first form checks every result expect.py lists, then scores it. Checks,
per result: it exists; its format, revision, dataset, scheme, seed and
inner folds; the five learners in order; each learner's outer calls summed
over the folds give the balanced accuracy the split's row holds; the folds'
test genomes sum to the tested genomes, by class; each inner pick covers
its outer fold's training genomes exactly and nothing else; its summed
inner calls, its inner balanced accuracy, its tied set and its pick follow
from its inner folds under the run's tie rule; at seed 42, each learner's
lineage balanced accuracy and the default's random-split figure equal the
release report's. A lineage split the run found unavailable is counted, not
failed.

Scores, per dataset, split, seed and arm (default, inner, twin): the served
calls summed over the outer folds, each fold taking the calls of the learner
it serves, and from them balanced accuracy, recall, specificity, the very
major and the major error rate, each with its count, denominator and Wilson
95% interval, and log loss with its denominator. The inner arm also counts
its folds per served learner, its ties, its folds with no pick and its
folds that serve the default's learner.

Writes <out dir>/checks.tsv, draws.tsv, datasets.tsv (per dataset, split
and arm: the mean of each measure over the draws, with the draws used) and
picks.tsv, and prints the check counts only. Exits 1 when a check fails.

The second form writes count.tsv and conditions.tsv (the plan's change
conditions 1 to 5 as the count reads them) from several species'
datasets.tsv, and prints the row count only. count.tsv holds, per
split, pair (inner against default, twin against default on the lineage
split) and measure, the datasets better, worse and the same by more than
MARGIN (lower is better for log loss and the two error rates), the median
difference, which for every measure is the arm's value less the default's,
and its species interval (RESAMPLES species resamples, read at
the 101st and 3,900th sorted values, each row seeded from the sha256 of
split|pair|measure), the one-sided sign test of better against worse, and
Holm over the lineage split's primary measure, the inner pick and the twin.
Every difference and median is rounded to 12 decimals before it is
compared.
"""

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

FORMAT = "grm-c9-pick-1"
REVISION = 14
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
RELEASE_SEED = 42
Z = 1.959963984540054
MARGIN = 0.02
RESAMPLES = 4000
LO_INDEX, HI_INDEX = 100, 3899
PRIMARY = "balanced_accuracy"
MEASURES = ["balanced_accuracy", "recall", "specificity", "vme", "me", "log_loss"]
LOWER_IS_BETTER = {"vme", "me", "log_loss"}
HOLM_ALPHA = 0.05
TIES = ("later", "earlier", "default")
NO_PICKS = ("default", "drop")
DEFAULTS = ("per-draw", "seed42")


def r12(x):
    return None if x is None else round(x, 12)


def wilson(k, n):
    if n == 0:
        return (None, None)
    p = k / n
    d = 1 + Z * Z / n
    c = (p + Z * Z / (2 * n)) / d
    h = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def balanced(c):
    if c["tp"] + c["fn"] == 0 or c["tn"] + c["fp"] == 0:
        return None
    return (c["tp"] / (c["tp"] + c["fn"]) + c["tn"] / (c["tn"] + c["fp"])) / 2


def add(a, b):
    return {k: a.get(k, 0) + b[k] for k in ("n", "tp", "fp", "tn", "fn", "loss_sum")}


ZERO = {"n": 0, "tp": 0, "fp": 0, "tn": 0, "fn": 0, "loss_sum": 0.0}


def measures(c):
    """Every measure of summed calls, with counts and denominators."""
    pos, neg = c["tp"] + c["fn"], c["tn"] + c["fp"]
    out = {"n": c["n"], "tp": c["tp"], "fp": c["fp"], "tn": c["tn"], "fn": c["fn"], "resistant": pos, "susceptible": neg,
           "balanced_accuracy": balanced(c), "log_loss": c["loss_sum"] / c["n"] if c["n"] else None}
    for name, k, n in (("recall", c["tp"], pos), ("specificity", c["tn"], neg), ("vme", c["fn"], pos), ("me", c["fp"], neg)):
        out[name] = k / n if n else None
        out[f"{name}_lo"], out[f"{name}_hi"] = wilson(k, n)
    return out


def pick(scores, tie, default_learner=None):
    """(learner, tied learners) of one fold's inner scores, rounded to 12 decimals; None when none."""
    have = {l: r12(s) for l, s in scores.items() if s is not None}
    if not have:
        return None, []
    best = max(have.values())
    tied = [l for l in LEARNERS if l in have and have[l] == best]
    if tie == "later":
        return tied[-1], tied
    if tie == "earlier":
        return tied[0], tied
    return (default_learner if default_learner in tied else tied[-1]), tied


# ---------- checks ----------

def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def check(exp, doc, release):
    """[(check name, passed, detail)] for one result against its expect row and the release report."""
    out = []

    def ok(name, cond, detail=""):
        out.append((name, bool(cond), detail))
        return bool(cond)

    ok("format and revision", doc.get("format") == FORMAT and doc.get("evaluation_revision") == REVISION)
    ok("dataset, scheme, seed and inner folds", (doc.get("dataset"), doc.get("scheme"), doc.get("seed"), doc.get("inner_folds"))
       == (exp["dataset"], exp["scheme"], int(exp["seed"]), int(exp["inner"])))
    ok("the five learners in order", [l["learner"] for l in doc.get("learners", [])] == LEARNERS)
    if doc.get("unavailable_reason") is not None:
        ok("an unavailable split has no fold", not doc.get("folds"))
        return out
    folds = doc["folds"]
    ok("folds in order", [f["fold"] for f in folds] == list(range(len(folds))) and len(folds) >= 2)
    ok("fold tests sum to the tested genomes", sum(f["test"] for f in folds) == doc["tested"]
       and sum(f["test_resistant"] for f in folds) == doc["tested_resistant"])
    good_outer = True
    for f in folds:
        good_outer &= [o["learner"] for o in f["outer"]] == LEARNERS
        good_outer &= all(o["n"] == f["test"] and o["tp"] + o["fn"] == f["test_resistant"] for o in f["outer"])
    ok("each fold's outer calls are its test genomes'", good_outer)
    if good_outer:
        for li, l in enumerate(LEARNERS):
            c = ZERO
            for f in folds:
                c = add(c, f["outer"][li])
            ok(f"{l}: the outer calls give the split's balanced accuracy",
               r12(balanced(c)) == r12(doc["learners"][li]["balanced_accuracy"]))
    tie = doc.get("tie")
    ok("tie rule", tie in ("later", "earlier"))
    for f in folds:
        inner = f["inner"]
        if int(exp["inner"]) == 0:
            ok(f"fold {f['fold']}: no inner fit", inner["status"] == "not run" and not inner["learners"])
            continue
        if inner["status"] == "no pick":
            ok(f"fold {f['fold']}: a fold with no pick says why", inner["pick"] is None and inner["reason"])
            continue
        parts = inner["inner_folds"]
        ok(f"fold {f['fold']}: the inner folds cover the training genomes exactly",
           len(parts) == int(exp["inner"]) and sum(p["test"] for p in parts) == f["train"]
           and sum(p["test_resistant"] for p in parts) == f["train_resistant"]
           and all(p["train"] + p["test"] == f["train"] for p in parts))
        scores = {}
        sums_ok = [l["learner"] for l in inner["learners"]] == LEARNERS
        for l in inner["learners"]:
            if l["error"] is not None:
                scores[l["learner"]] = None
                sums_ok &= l["balanced_accuracy"] is None
                continue
            c = ZERO
            for p, part in zip(l["folds"], parts):
                sums_ok &= p["n"] == part["test"] and p["tp"] + p["fn"] == part["test_resistant"]
                c = add(c, p)
            sums_ok &= len(l["folds"]) == len(parts)
            sums_ok &= (c["tp"], c["fp"], c["tn"], c["fn"]) == (l["tp"], l["fp"], l["tn"], l["fn"])
            sums_ok &= r12(balanced(c)) == r12(l["balanced_accuracy"])
            scores[l["learner"]] = l["balanced_accuracy"]
        ok(f"fold {f['fold']}: inner calls sum to the inner figures", sums_ok)
        p, tied = pick(scores, tie)
        ok(f"fold {f['fold']}: the pick and tied set follow from the inner figures",
           inner["status"] == "picked" and inner["pick"] == p and inner["tied"] == tied)
    if exp["split"] == "lineage":
        d = doc.get("default")
        ok("a lineage result names the default's pick", d is not None and d["learner"] in LEARNERS
           and d["learner"] in d["tied"])
        if int(exp["seed"]) == RELEASE_SEED and release is not None and d is not None:
            rel = {s["scheme"]: s for s in release["splits"]}

            def figure(scheme, learner):
                rows = [l for l in rel[scheme]["learners"] if l["learner"] == learner]
                return r12((rows[0].get("metrics") or {}).get("balanced_accuracy")) if rows else None

            ok("seed 42: each learner's lineage figure is the release report's",
               all(r12(doc["learners"][i]["balanced_accuracy"]) == figure("lineage", l) for i, l in enumerate(LEARNERS)))
            ok("seed 42: the default's random-split figure is the release report's",
               r12(d["random_balanced_accuracy"]) == figure("random", d["learner"]))
    else:
        ok("a study result names no default of its own", doc.get("default") is None)
    return out


# ---------- scores ----------

def served(doc, default_learner, tie, no_pick):
    """The inner arm's summed calls and pick counts, or None when the draw is not scored."""
    c = ZERO
    counts = {"picks_" + l: 0 for l in LEARNERS}
    counts.update({"ties": 0, "no_pick": 0, "same_as_default": 0, "folds": len(doc["folds"])})
    for f in doc["folds"]:
        inner = f["inner"]
        if inner["status"] == "picked":
            scores = {l["learner"]: l["balanced_accuracy"] for l in inner["learners"]}
            learner, tied = pick(scores, tie, default_learner)
            counts["ties"] += len(tied) > 1
        else:
            if no_pick == "drop":
                return None, counts
            learner = default_learner
            counts["no_pick"] += 1
        counts["picks_" + learner] += 1
        counts["same_as_default"] += learner == default_learner
        c = add(c, f["outer"][LEARNERS.index(learner)])
    return c, counts


def fixed(doc, learner):
    """One learner's calls summed over the outer test folds."""
    c = ZERO
    for f in doc["folds"]:
        c = add(c, f["outer"][LEARNERS.index(learner)])
    return c


def score_rows(expect, docs, names, tie, no_pick, default_rule):
    """draws.tsv's rows: per result and arm."""
    rows = []
    seed42 = {}
    for e in expect:
        if e["split"] == "lineage" and int(e["seed"]) == RELEASE_SEED:
            doc = docs[e["file"]]
            if doc.get("default"):
                seed42[e["dataset"]] = doc["default"]["learner"]
    for e in expect:
        doc = docs[e["file"]]
        sp, ab = names[e["dataset"]]
        base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "split": e["split"], "seed": int(e["seed"]),
                "role": e["role"]}
        if doc.get("unavailable_reason") is not None:
            for arm in (("default", "inner") if e["role"] != "twin" else ("twin",)):
                rows.append(dict(base, arm=arm, scored=0, why="split unavailable"))
            continue
        if e["split"] == "study" or default_rule == "seed42":
            learner = seed42.get(e["dataset"])
        else:
            learner = doc["default"]["learner"]
        if learner is None:
            for arm in (("default", "inner") if e["role"] != "twin" else ("twin",)):
                rows.append(dict(base, arm=arm, scored=0, why="no seed-42 default"))
            continue
        if e["role"] == "twin":
            rows.append(dict(base, arm="twin", scored=1, why="", learner=learner, **measures(fixed(doc, learner))))
            continue
        rows.append(dict(base, arm="default", scored=1, why="", learner=learner, **measures(fixed(doc, learner))))
        c, counts = served(doc, learner, tie, no_pick)
        if c is None:
            rows.append(dict(base, arm="inner", scored=0, why="a fold has no pick", **counts))
        else:
            rows.append(dict(base, arm="inner", scored=1, why="", learner="", **counts, **measures(c)))
    return rows


def dataset_rows(draws):
    """datasets.tsv's rows: per dataset, split and arm, the mean over its draws; unscored if any draw is."""
    keys = sorted({(r["dataset"], r["split"], r["arm"]) for r in draws})
    out = []
    for d, s, a in keys:
        rs = [r for r in draws if (r["dataset"], r["split"], r["arm"]) == (d, s, a)]
        row = {"dataset": d, "species_slug": d.split("__")[0], "species": rs[0]["species"], "antibiotic": rs[0]["antibiotic"],
               "split": s, "arm": a, "draws": len(rs), "scored_draws": sum(r["scored"] for r in rs)}
        row["scored"] = int(row["scored_draws"] == row["draws"])
        for m in MEASURES:
            vals = [r.get(m) for r in rs]
            row[m] = r12(sum(vals) / len(vals)) if row["scored"] and all(v is not None for v in vals) else None
        for k in ("n", "tp", "fp", "tn", "fn"):
            row[k] = sum(r.get(k, 0) for r in rs) if row["scored"] else None
        out.append(row)
    return out


# ---------- the count ----------

def row_seed(split, pair, measure):
    return int(hashlib.sha256(f"{split}|{pair}|{measure}".encode()).hexdigest(), 16)


def sign_p(better, worse):
    """One-sided sign test: P(at least `better` of better + worse by a fair coin)."""
    n = better + worse
    if n == 0:
        return 1.0
    return sum(math.comb(n, k) for k in range(better, n + 1)) / 2 ** n


def holm(ps):
    """Holm-adjusted p values, in the input order."""
    order = sorted(range(len(ps)), key=lambda i: ps[i])
    out, run = [None] * len(ps), 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - rank) * ps[i]))
        out[i] = run
    return out


PAIRS = {"lineage": [("inner", "default"), ("twin", "default")], "study": [("inner", "default")]}


def count(rows):
    """count.tsv's rows from several species' datasets.tsv rows."""
    cell = {(r["dataset"], r["split"], r["arm"]): r for r in rows}
    out = []
    for split, pairs in PAIRS.items():
        for arm, base in pairs:
            for m in MEASURES:
                ds = sorted({d for d, s, a in cell if s == split and a == arm})
                diffs, sps, unscored = [], [], 0
                for d in ds:
                    a, b = cell[(d, split, arm)], cell.get((d, split, base))
                    if b is None or not int(a["scored"]) or not int(b["scored"]) or a[m] in (None, "") or b[m] in (None, ""):
                        unscored += 1
                        continue
                    diffs.append(r12(float(a[m]) - float(b[m])))
                    sps.append(a["species_slug"])
                # Every difference, median and interval is the arm's value less the default's;
                # only better and worse turn round where lower is better.
                up = -1 if m in LOWER_IS_BETTER else 1
                better = sum(1 for x in diffs if up * x > MARGIN)
                worse = sum(1 for x in diffs if up * x < -MARGIN)
                row = {"split": split, "pair": f"{arm}-{base}", "measure": m, "datasets": len(ds), "entered": len(diffs),
                       "not_scored": unscored, "better": better, "worse": worse, "same": len(diffs) - better - worse,
                       "median_difference": r12(float(np.median(diffs))) if diffs else None,
                       "species": len(set(sps)), "interval_lo": None, "interval_hi": None,
                       "sign_p": sign_p(better, worse), "holm_p": None, "orientation": "lower" if m in LOWER_IS_BETTER else "higher"}
                if diffs:
                    by = {}
                    for s, x in zip(sps, diffs):
                        by.setdefault(s, []).append(x)
                    keys = sorted(by)
                    rng = np.random.default_rng(row_seed(split, row["pair"], m))
                    meds = []
                    for _ in range(RESAMPLES):
                        take = rng.integers(0, len(keys), len(keys))
                        meds.append(r12(float(np.median(np.concatenate([by[keys[i]] for i in take])))))
                    meds.sort()
                    row["interval_lo"], row["interval_hi"] = meds[LO_INDEX], meds[HI_INDEX]
                out.append(row)
    family = [r for r in out if r["split"] == "lineage" and r["measure"] == PRIMARY]
    for r, p in zip(family, holm([r["sign_p"] for r in family])):
        r["holm_p"] = p
    return out


def conditions(count_rows):
    """conditions.tsv's rows: the plan's change conditions 1 to 5, read on the count (6 is the maintainer's)."""
    at = {(r["split"], r["pair"], r["measure"]): r for r in count_rows}

    def c1(pair):
        r = at[("lineage", pair, PRIMARY)]
        return r["interval_lo"] is not None and r["interval_lo"] > 0

    def c2(pair):
        r = at[("lineage", pair, PRIMARY)]
        return r["holm_p"] is not None and r["holm_p"] < HOLM_ALPHA and r["better"] > r["worse"]

    recall = at[("lineage", "inner-default", "recall")]
    lin = at[("lineage", "inner-default", PRIMARY)]["median_difference"]
    study = at[("study", "inner-default", PRIMARY)]["median_difference"]
    sign = lambda x: None if x is None else (x > 0) - (x < 0)
    rows = [
        (1, "the lineage median difference's interval lies wholly above 0", c1("inner-default")),
        (2, "better beats worse by the one-sided sign test at p < 0.05 after Holm", c2("inner-default")),
        (3, "the twin meets neither 1 nor 2", not c1("twin-default") and not c2("twin-default")),
        (4, "the recall interval does not lie wholly below -0.02",
         recall["interval_hi"] is not None and recall["interval_hi"] >= -MARGIN),
        (5, "the study median has the lineage median's sign",
         lin is not None and study is not None and sign(lin) != 0 and sign(study) == sign(lin)),
    ]
    return [{"condition": n, "reads": text, "met": int(bool(m))} for n, text, m in rows] + [
        {"condition": 6, "reads": "the maintainer rules the change", "met": ""}]


# ---------- files ----------

def fmt(x):
    if x is None:
        return ""
    if isinstance(x, bool):
        return str(int(x))
    if isinstance(x, float):
        return "" if math.isnan(x) else repr(round(x, 12))
    return str(x)


def write_tsv(path, rows, columns=None):
    columns = columns or sorted({k for r in rows for k in r})
    with open(path, "w", newline="") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(c)) for c in columns) + "\n")


DRAW_COLUMNS = ["dataset", "species", "antibiotic", "split", "seed", "role", "arm", "scored", "why", "learner", "n", "resistant",
                "susceptible", "tp", "fp", "tn", "fn"] + [c for m in ("recall", "specificity", "vme", "me")
                                                           for c in (m, f"{m}_lo", f"{m}_hi")] + [
    "balanced_accuracy", "log_loss", "folds", "ties", "no_pick", "same_as_default"] + ["picks_" + l for l in LEARNERS]
DATASET_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "split", "arm", "draws", "scored_draws", "scored",
                   "n", "tp", "fp", "tn", "fn"] + MEASURES
COUNT_COLUMNS = ["split", "pair", "measure", "orientation", "datasets", "entered", "not_scored", "better", "worse", "same",
                 "median_difference", "species", "interval_lo", "interval_hi", "sign_p", "holm_p"]


def options(argv):
    opts = {"--tie": "later", "--no-pick": "default", "--default": "per-draw"}
    allowed = {"--tie": TIES, "--no-pick": NO_PICKS, "--default": DEFAULTS}
    rest = []
    i = 0
    while i < len(argv):
        if argv[i] in opts:
            if i + 1 >= len(argv) or argv[i + 1] not in allowed[argv[i]]:
                raise SystemExit(f"{argv[i]} takes one of {', '.join(allowed[argv[i]])}")
            opts[argv[i]] = argv[i + 1]
            i += 2
        else:
            rest.append(argv[i])
            i += 1
    return opts, rest


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            print(__doc__.strip().splitlines()[3], file=sys.stderr)
            return 2
        out = Path(argv[1])
        out.mkdir(parents=True, exist_ok=True)
        rows = [r for p in argv[2:] for r in read_tsv(p)]
        keys = [(r["dataset"], r["split"], r["arm"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise SystemExit("a dataset, split and arm is listed twice")
        counted = count(rows)
        write_tsv(out / "count.tsv", counted, COUNT_COLUMNS)
        write_tsv(out / "conditions.tsv", conditions(counted), ["condition", "reads", "met"])
        print(f"count: {len(rows)} dataset rows from {len(argv) - 2} files")
        return 0
    opts, rest = options(argv)
    if len(rest) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    expect_path, picks, release, out = Path(rest[0]), Path(rest[1]), Path(rest[2]), Path(rest[3])
    expect = read_tsv(expect_path)
    checks, docs, names = [], {}, {}
    for e in expect:
        path = picks / e["file"]
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        rel = None
        if rel_path.is_file():
            with open(rel_path) as f:
                rel = json.load(f)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
        if not path.is_file():
            checks.append((e["file"], "exists", False, ""))
            continue
        with open(path) as f:
            doc = json.load(f)
        docs[e["file"]] = doc
        checks.append((e["file"], "exists", True, ""))
        checks.append((e["file"], "release report", rel is not None, ""))
        for name, passed, detail in check(e, doc, rel):
            checks.append((e["file"], name, passed, detail))
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "checks.tsv", [{"file": f, "check": c, "passed": p, "detail": d} for f, c, p, d in checks],
              ["file", "check", "passed", "detail"])
    failed = [c for c in checks if not c[2]]
    unavailable = sum(1 for d in docs.values() if d.get("unavailable_reason") is not None)
    print(f"checks: {len(checks) - len(failed)} of {len(checks)} pass; results {len(docs)} of {len(expect)}; "
          f"lineage or study split unavailable in {unavailable}")
    if failed:
        return 1
    draws = score_rows(expect, docs, names, opts["--tie"], opts["--no-pick"], opts["--default"])
    write_tsv(out / "draws.tsv", draws, DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", dataset_rows(draws), DATASET_COLUMNS)
    picks_rows = [{k: r.get(k) for k in ["dataset", "split", "seed", "folds", "ties", "no_pick", "same_as_default"]
                   + ["picks_" + l for l in LEARNERS]} for r in draws if r["arm"] == "inner"]
    write_tsv(out / "picks.tsv", picks_rows, ["dataset", "split", "seed", "folds", "ties", "no_pick", "same_as_default"]
              + ["picks_" + l for l in LEARNERS])
    with open(out / "options.txt", "w") as f:
        for k, v in opts.items():
            f.write(f"{k} {v}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
