"""The scorer of plan block C3 (docs/input-arms.md, measures and files).

  python3 -I score.py <expect.tsv> <picks dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> <datasets.tsv> [<datasets.tsv> ...]

The first form checks every result expect.py lists, then scores it. Each
result is a grm learn --exp-pick file at no inner fold, so the checks of
scripts/c9-pick/score.py hold for it (format, revision, dataset, scheme,
seed, no inner fold, the five learners in order, each fold's calls are its
test genomes', the outer calls give the split's balanced accuracy, a study
result names no default), and at seed 42 the default's lineage figures and
served learner must be the release report's. Beside them, per result:

- an arm's folds are the default's at the same dataset, split and seed:
  every fold's training and test genomes by class, and the same split
  unavailable (the twin is the default at seeds of its own, so it has none
  to match);
- the default and the twin write no inputs file, and a twin seed is no
  draw's seed;
- an arm's inputs file has the format, revision, dataset, arm and tag of
  its run, read a table for every cohort genome, and its rows add up: the
  program names seen are those out of range, collapsed and kept, and the
  arm's rows are the program's, plus K31's for k31-program.

Scores, per dataset, arm, split, seed and learner: the learner's outer
calls summed over the folds, and from them balanced accuracy, recall,
specificity, the very major and the major error rate, each with its count,
denominator and Wilson 95% interval, and log loss. Writes checks.tsv,
draws.tsv and datasets.tsv (per dataset, split, arm and learner: the mean of
each measure over the draws), and prints the check counts only. Exits 1,
scoring nothing, when a check fails.

The second form writes count.tsv and conditions.tsv over several species'
datasets.tsv, and prints the row count only. It refuses unless every arm of
the family (unitigs, k31-counts, k31-aa, program, k31-program and the twin)
has lineage rows: the family is fixed, and the count is taken once. Per
split, arm, learner and measure, each arm against the default (k31) on the
arm's own datasets: the datasets better, worse and the same by more than
MARGIN (lower is better for log loss and the two error rates), the median
difference (always the arm's value less the default's) and its interval
from RESAMPLES species resamples (the 101st and 3,900th sorted values, each
row seeded from the sha256 of split|pair|learner|measure), the one-sided
sign test of better against worse, and Holm over the 30 lineage rows of the
primary measure, the twin's included. Every difference and median is
rounded to 12 decimals before it is compared.
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _c9():
    spec = importlib.util.spec_from_file_location("c9_score_for_c3", HERE.parent / "c9-pick" / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C9 = _c9()
LEARNERS = C9.LEARNERS
MEASURES = C9.MEASURES
LOWER_IS_BETTER = C9.LOWER_IS_BETTER
PRIMARY = C9.PRIMARY
MARGIN = C9.MARGIN
RESAMPLES = C9.RESAMPLES
LO_INDEX, HI_INDEX = C9.LO_INDEX, C9.HI_INDEX
HOLM_ALPHA = C9.HOLM_ALPHA
REVISION = C9.REVISION
RELEASE_SEED = C9.RELEASE_SEED
r12, measures, fixed, read_tsv, write_tsv = C9.r12, C9.measures, C9.fixed, C9.read_tsv, C9.write_tsv
INPUT_FORMAT = "grm-c3-inputs-1"
DEFAULT = "k31"
TWIN = "twin"
FAMILY = ["unitigs", "k31-counts", "k31-aa", "program", "k31-program", TWIN]
JOINED = {"k31-program", "k31-counts", "k31-aa"}
TOUCHED_AT_LEAST = 3


def fold_shape(doc):
    if doc.get("unavailable_reason") is not None:
        return ("unavailable",)
    return tuple((f["train"], f["train_resistant"], f["test"], f["test_resistant"]) for f in doc["folds"]) + (
        (doc["tested"], doc["tested_resistant"]),)


def check(exp, doc, base, inputs, stray_inputs, release=None, served_learner=None, draw_seeds=(), twin_seeds=()):
    """[(check name, passed, detail)] for one result.

    base is the default's result at the same dataset, split and seed (None
    when this is the default's own); inputs is the arm's inputs file (None
    when it has none); stray_inputs counts inputs files in the default's or
    the twin's folder for this dataset.
    """
    arm = exp["arm"]
    e9 = dict(exp, inner=0)
    is_default = arm == DEFAULT
    out = C9.check(e9, doc, release if is_default else None, served_learner if is_default else None)

    def ok(name, cond, detail=""):
        out.append((name, bool(cond), detail))
        return bool(cond)

    seed = int(exp["seed"])
    if is_default and seed == RELEASE_SEED and exp["split"] == "lineage":
        ok("seed 42: the default is held to the release report", release is not None and served_learner is not None)
    if arm in (DEFAULT, TWIN):
        ok("the default and the twin write no inputs file", inputs is None and stray_inputs == 0)
    if arm == TWIN:
        ok("a twin seed is the twin's and no draw's", seed in twin_seeds and seed not in draw_seeds)
    if arm not in (DEFAULT, TWIN):
        ok("the folds are the default's at the same dataset, split and seed, by class",
           base is not None and fold_shape(base) == fold_shape(doc))
    if arm not in (DEFAULT, TWIN):
        tag = "drawn" if exp["split"] == "lineage" else exp["scheme"]
        good = inputs is not None
        if good:
            ok("inputs: format, revision, dataset, arm and tag",
               (inputs.get("format"), inputs.get("evaluation_revision"), inputs.get("dataset"), inputs.get("arm"),
                inputs.get("tag")) == (INPUT_FORMAT, REVISION, exp["dataset"], arm, tag))
            ok("inputs: a table read for every cohort genome", inputs["tables_read"] == inputs["samples"] > 0)
            ok("inputs: the program names seen are those out of range, collapsed and kept",
               inputs["program_names_seen"] == inputs["program_names_out_of_range"]
               + inputs["program_names_collapsed"] + inputs["program_rows"])
            joined = arm in JOINED
            dropped = inputs["k31_rows_as_program"]
            ok("inputs: the arm's rows are its blocks', each once",
               inputs["arm_rows"] == inputs["program_rows"] + (inputs["k31_rows"] - dropped if joined else 0)
               and 0 <= dropped <= (min(inputs["program_rows"], inputs["k31_rows"]) if joined else 0))
        else:
            ok("inputs: the arm's file exists", False)
    return out


# ---------- scores ----------

def score_rows(expect, docs, names):
    """draws.tsv's rows: per result and learner."""
    rows = []
    for e in expect:
        doc = docs[e["file"]]
        sp, ab = names[e["dataset"]]
        base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "split": e["split"],
                "seed": int(e["seed"]), "role": e["role"]}
        for learner in LEARNERS:
            if doc.get("unavailable_reason") is not None:
                rows.append(dict(base, learner=learner, scored=0, why="split unavailable"))
            else:
                rows.append(dict(base, learner=learner, scored=1, why="", **measures(fixed(doc, learner))))
    return rows


def dataset_rows(draws):
    """datasets.tsv's rows: per dataset, split, arm and learner, the mean over its draws; unscored if any draw is."""
    groups = {}
    for r in draws:
        groups.setdefault((r["dataset"], r["split"], r["arm"], r["learner"]), []).append(r)
    out = []
    for (d, s, a, learner), rs in sorted(groups.items()):
        row = {"dataset": d, "species_slug": d.split("__")[0], "species": rs[0]["species"], "antibiotic": rs[0]["antibiotic"],
               "split": s, "arm": a, "learner": learner, "draws": len(rs), "scored_draws": sum(r["scored"] for r in rs)}
        row["scored"] = int(row["scored_draws"] == row["draws"])
        for m in MEASURES:
            vals = [r.get(m) for r in rs]
            row[m] = r12(sum(vals) / len(vals)) if row["scored"] and all(v is not None for v in vals) else None
        for k in ("n", "tp", "fp", "tn", "fn"):
            row[k] = sum(r.get(k, 0) for r in rs) if row["scored"] else None
        out.append(row)
    return out


# ---------- the count ----------

def row_seed(split, pair, learner, measure):
    return C9.row_seed(split, f"{pair}|{learner}", measure)


def count(rows):
    """count.tsv's rows from several species' datasets.tsv rows."""
    have = {r["arm"] for r in rows if r["split"] == "lineage"}
    missing = [a for a in FAMILY if a not in have]
    if missing:
        raise SystemExit(f"the family is not complete: no lineage result for {', '.join(missing)}; "
                         "the count is taken once, over every arm")
    cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}
    out = []
    for split, arms in (("lineage", FAMILY), ("study", [a for a in FAMILY if a != TWIN])):
        for arm in arms:
            pair = f"{arm}-{DEFAULT}"
            for learner in LEARNERS:
                for m in MEASURES:
                    ds = sorted({d for d, s, a, l in cell if s == split and a == arm and l == learner})
                    diffs, sps, unscored = [], [], 0
                    for d in ds:
                        a, b = cell[(d, split, arm, learner)], cell.get((d, split, DEFAULT, learner))
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
                    row = {"split": split, "pair": pair, "learner": learner, "measure": m, "datasets": len(ds),
                           "entered": len(diffs), "not_scored": unscored, "better": better, "worse": worse,
                           "same": len(diffs) - better - worse,
                           "median_difference": r12(float(np.median(diffs))) if diffs else None,
                           "species": len(set(sps)), "interval_lo": None, "interval_hi": None,
                           "sign_p": C9.sign_p(better, worse), "holm_p": None,
                           "orientation": "lower" if m in LOWER_IS_BETTER else "higher"}
                    if diffs and len(set(diffs)) == 1:
                        # Every resample's median is that one value.
                        row["interval_lo"] = row["interval_hi"] = diffs[0]
                    elif diffs:
                        by = {}
                        for s, x in zip(sps, diffs):
                            by.setdefault(s, []).append(x)
                        keys = sorted(by)
                        rng = np.random.default_rng(row_seed(split, pair, learner, m))
                        meds = []
                        for _ in range(RESAMPLES):
                            take = rng.integers(0, len(keys), len(keys))
                            meds.append(r12(float(np.median(np.concatenate([by[keys[i]] for i in take])))))
                        meds.sort()
                        row["interval_lo"], row["interval_hi"] = meds[LO_INDEX], meds[HI_INDEX]
                    out.append(row)
    family = [r for r in out if r["split"] == "lineage" and r["measure"] == PRIMARY]
    assert len(family) == len(FAMILY) * len(LEARNERS)
    for r, p in zip(family, C9.holm([r["sign_p"] for r in family])):
        r["holm_p"] = p
    return out


def conditions(count_rows):
    """conditions.tsv's rows: per arm, the plan's change conditions 1 to 5 as the count reads them.

    An input arm touches all five learners, so 1 and 2 need at least 3 of 5
    learners, 4 holds for all five, and 5 asks each learner that meets 1 for
    a study median of its lineage median's sign. 3 reads the twin as 1 and
    2 read an arm. 6 is the maintainer's.
    """
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in count_rows}
    sign = lambda x: None if x in (None, "") else (float(x) > 0) - (float(x) < 0)

    def c1(arm):
        rs = [at[("lineage", f"{arm}-{DEFAULT}", l, PRIMARY)] for l in LEARNERS]
        return [l for l, r in zip(LEARNERS, rs) if r["interval_lo"] not in (None, "") and float(r["interval_lo"]) > 0]

    def c2(arm):
        rs = [at[("lineage", f"{arm}-{DEFAULT}", l, PRIMARY)] for l in LEARNERS]
        return [l for l, r in zip(LEARNERS, rs) if r["holm_p"] not in (None, "") and float(r["holm_p"]) < HOLM_ALPHA
                and int(r["better"]) > int(r["worse"])]

    twin_quiet = len(c1(TWIN)) < TOUCHED_AT_LEAST and len(c2(TWIN)) < TOUCHED_AT_LEAST
    out = []
    for arm in [a for a in FAMILY if a != TWIN]:
        one, two = c1(arm), c2(arm)
        low = [l for l in LEARNERS if at[("lineage", f"{arm}-{DEFAULT}", l, "recall")]["interval_hi"] not in (None, "")
               and float(at[("lineage", f"{arm}-{DEFAULT}", l, "recall")]["interval_hi"]) < -MARGIN]
        study = at.get(("study", f"{arm}-{DEFAULT}", LEARNERS[0], PRIMARY))
        same_sign = [l for l in one
                     if sign(at[("lineage", f"{arm}-{DEFAULT}", l, PRIMARY)]["median_difference"]) not in (None, 0)
                     and sign(at[("study", f"{arm}-{DEFAULT}", l, PRIMARY)]["median_difference"])
                     == sign(at[("lineage", f"{arm}-{DEFAULT}", l, PRIMARY)]["median_difference"])] if study else []
        rows = [
            (1, "the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners",
             len(one) >= TOUCHED_AT_LEAST, len(one)),
            (2, "better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for at least 3 of 5",
             len(two) >= TOUCHED_AT_LEAST, len(two)),
            (3, "the twin meets neither 1 nor 2", twin_quiet, ""),
            (4, "no learner's recall interval lies wholly below -0.02", not low, 5 - len(low)),
            (5, "each learner meeting 1 has a study median of its lineage median's sign",
             bool(one) and len(same_sign) == len(one), len(same_sign)),
        ]
        out += [{"arm": arm, "condition": n, "reads": text, "met": int(bool(m)), "learners": k} for n, text, m, k in rows]
        out.append({"arm": arm, "condition": 6, "reads": "the maintainer rules the change", "met": "", "learners": ""})
    return out


# ---------- files ----------

DRAW_COLUMNS = ["dataset", "species", "antibiotic", "arm", "split", "seed", "role", "learner", "scored", "why", "n",
                "resistant", "susceptible", "tp", "fp", "tn", "fn"] + [c for m in ("recall", "specificity", "vme", "me")
                                                                         for c in (m, f"{m}_lo", f"{m}_hi")] + [
    "balanced_accuracy", "log_loss"]
DATASET_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "split", "arm", "learner", "draws", "scored_draws",
                   "scored", "n", "tp", "fp", "tn", "fn"] + MEASURES
COUNT_COLUMNS = ["split", "pair", "learner", "measure", "orientation", "datasets", "entered", "not_scored", "better",
                 "worse", "same", "median_difference", "species", "interval_lo", "interval_hi", "sign_p", "holm_p"]


def load(path):
    with open(path) as f:
        return json.load(f)


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            print(__doc__.strip().splitlines()[3], file=sys.stderr)
            return 2
        out = Path(argv[1])
        rows = [r for p in argv[2:] for r in read_tsv(p)]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise SystemExit("a dataset, split, arm and learner is listed twice")
        counted = count(rows)
        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "count.tsv", counted, COUNT_COLUMNS)
        write_tsv(out / "conditions.tsv", conditions(counted), ["arm", "condition", "reads", "met", "learners"])
        print(f"count: {len(rows)} dataset rows from {len(argv) - 2} files")
        return 0
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    expect_path, picks, release, out = (Path(a) for a in argv)
    expect = read_tsv(expect_path)
    draw_seeds = {int(e["seed"]) for e in expect if e["role"] == "draw"}
    twin_seeds = {int(e["seed"]) for e in expect if e["role"] == "twin"}
    served_by = {}
    if (release / "index.json").is_file():
        served_by = {x["id"]: x.get("best_learner") for x in load(release / "index.json").get("datasets", [])}
    docs, names, checks = {}, {}, []
    for e in expect:
        path = picks / e["file"]
        if path.is_file():
            docs[e["file"]] = load(path)
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        if rel_path.is_file() and e["dataset"] not in names:
            rel = load(rel_path)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in expect if e["arm"] == DEFAULT}
    for e in expect:
        if e["file"] not in docs:
            checks.append((e["file"], "exists", False, ""))
            continue
        checks.append((e["file"], "exists", True, ""))
        checks.append((e["file"], "release report", e["dataset"] in names, ""))
        rel = None
        if e["arm"] == DEFAULT and int(e["seed"]) == RELEASE_SEED and e["dataset"] in names:
            rel = load(release / "datasets" / f"{e['dataset']}.json")
        base_file = base_of.get((e["dataset"], e["split"], e["seed"]))
        base = docs.get(base_file) if e["arm"] != DEFAULT else None
        inputs = load(picks / e["inputs"]) if e["inputs"] != "-" and (picks / e["inputs"]).is_file() else None
        stray = 0
        if e["arm"] in (DEFAULT, TWIN):
            stray = sum(1 for _ in (picks / e["arm"]).glob(f"{e['dataset']}.*.inputs.json"))
        for name, passed, detail in check(e, docs[e["file"]], base, inputs, stray, rel, served_by.get(e["dataset"]),
                                          draw_seeds, twin_seeds):
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
    draws = score_rows(expect, docs, names)
    write_tsv(out / "draws.tsv", draws, DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", dataset_rows(draws), DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
