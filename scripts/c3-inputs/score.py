"""The scorer of plan block C3 (docs/input-arms.md, measures and files).

  python3 -I score.py <expect.tsv> <picks dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> [--withdrawn <arm> ...] <datasets.tsv> [<datasets.tsv> ...]

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
  its run, and its rows add up. A program arm read a table for every cohort
  genome, the program names seen are those out of range, collapsed and
  kept, and the arm's rows are the program's, plus K31's for k31-program. A
  matrix arm (k31-counts, k31-aa, unitigs) names its build, holds no
  program field, its added rows are at most its build's species patterns,
  and its rows are the added block's, plus K31's for a joined arm, each K31
  row left out as an added row's pattern counted; the unitig arm counts the
  cohort's genome and unitig pairs held whole and held in part;
- a joined arm's result (k31-program, k31-counts, k31-aa) that has folds has
  its screen count, <arm>/<dataset>.<scheme>-s<seed>.screen.json, of its
  dataset, scheme, seed, arm and first block, with one row per outer fold of
  the result: the fold's training genomes, and for the forest, boosted trees
  and lasso in that order, the forest's fitting genomes the fold's training
  genomes and the other two's fewer, the screened count the 4,000 or all
  candidates when fewer, and each block's passed count inside its
  candidates, the two summing to the screened count;
- an arm's result that grm learn refused for no candidate (a --fold-file
  fit stops when an outer fold's training genomes give no candidate; the
  run script writes <arm>/<dataset>.<scheme>-s<seed>.refused.json in place
  of the result) names its run and that reason, and the default's result at
  the same dataset, split and seed is scored.

Scores, per dataset, arm, split, seed and learner: the learner's outer
calls summed over the folds, and from them balanced accuracy, recall,
specificity, the very major and the major error rate, each with its count,
denominator and Wilson 95% interval, and log loss. A refused arm result is
scored as the constant model on the default's folds of that dataset, split
and seed (the lead's ruling 16 of 23 Nov 2026): every test genome of a fold
gets the majority class of the fold's training genomes, R on a tie, and a
score equal to their R share, clipped as log loss clips; its draws.tsv rows
say so in `constant`. Writes checks.tsv, draws.tsv and datasets.tsv (per
dataset, split, arm and learner: the mean of each measure over the draws,
and how many draws are the constant model), and prints the check counts
only. Exits 1, scoring nothing, when a check fails.

The second form writes count.tsv, count-no-constant.tsv and conditions.tsv
over several species' datasets.tsv, and prints the row count only. It
refuses unless every arm of the family (unitigs, k31-counts, k31-aa,
program, k31-program and the twin) has lineage rows or is named with
--withdrawn: the family is fixed at 30 pairs, and the count is taken once.
An arm withdrawn before it runs keeps its five pairs in the family at
p = 1 (the lead's addition to question 12); it may hold no row. Per
split, arm, learner and measure, each arm against the default (k31) on the
arm's own datasets: the datasets better, worse and the same by more than
MARGIN (lower is better for log loss and the two error rates), the median
difference (always the arm's value less the default's) and its interval
from RESAMPLES species resamples (the 101st and 3,900th sorted values, each
row seeded from the sha256 of split|pair|learner|measure), the one-sided
sign test of better against worse, and Holm over the 30 lineage rows of the
primary measure, the twin's included. Every difference and median is
rounded to 12 decimals before it is compared. Datasets with a constant-model
draw enter the count, and count.tsv gives their number per row (`constant`);
count-no-constant.tsv is the same count with them left out, with Holm over
its own 30 rows, and decides nothing.

conditions.tsv reads the plan's conditions per arm two ways. Ruling 15, the
one that decides: a learner carries the change when it meets 1 and 2 both,
and at least 3 of 5 must carry it; 3 holds when the twin meets 1 in fewer
than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart, as the
plan words it (the twin meets neither 1 nor 2); 4 holds for all five
learners; 5 holds for each learner that carries. The stricter reading,
written beside it and deciding nothing: 1 for at least 3 of 5, 2 and 5 for
all five, 3 the twin meeting neither of those, 4 for all five. 6 is the
maintainer's.
"""

import importlib.util
import json
import math
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
PROGRAM_ARMS = {"program", "k31-program"}
BUILD_OF = {"k31-counts": "counts", "k31-aa": "aa", "unitigs": "unitigs"}
FIRST_BLOCK = {"k31-program": "program", "k31-counts": "counts", "k31-aa": "aa"}
PROGRAM_FIELDS = ("tables_read", "program_names_seen", "program_names_out_of_range", "program_names_collapsed",
                  "program_rows", "k31_rows_as_program")
SCREEN_FORMAT = "grm-c3-screen-1"
SCREENED = ("rf", "gbt", "lasso")
SCREENED_FEATURES = 4000
TOUCHED_AT_LEAST = 3
REFUSED_FORMAT = "grm-c3-refused-1"
NO_CANDIDATE = "no candidate"
# Log loss clips each score to [LOSS_CLIP, 1 - LOSS_CLIP], as grm does (src/cutoff.rs).
LOSS_CLIP = 1e-4


def refused_name(file):
    """The refusal record the run script writes in place of an arm's result file."""
    return file[: -len(".pick.json")] + ".refused.json"


def constant_calls(base):
    """The constant model of ruling 16 on the default's folds `base`: one call per fold, summed."""
    c = dict(C9.ZERO)
    for f in base["folds"]:
        r, s = f["test_resistant"], f["test"] - f["test_resistant"]
        share = f["train_resistant"] / f["train"]
        p = min(max(share, LOSS_CLIP), 1 - LOSS_CLIP)
        resistant = share >= 0.5
        c = C9.add(c, {"n": r + s, "tp": r if resistant else 0, "fp": s if resistant else 0,
                       "tn": 0 if resistant else s, "fn": 0 if resistant else r,
                       "loss_sum": -(r * math.log(p) + s * math.log(1 - p))})
    return c


def fold_shape(doc):
    if doc.get("unavailable_reason") is not None:
        return ("unavailable",)
    return tuple((f["train"], f["train_resistant"], f["test"], f["test_resistant"]) for f in doc["folds"]) + (
        (doc["tested"], doc["tested_resistant"]),)


def check_refused(exp, refusal, base):
    """[(check name, passed, detail)] for an arm result grm learn refused for no candidate (ruling 16)."""
    out = []

    def ok(name, cond):
        out.append((name, bool(cond), ""))

    ok("refused: format, dataset, arm, scheme and seed of its run",
       (refusal.get("format"), refusal.get("dataset"), refusal.get("arm"), refusal.get("scheme"), refusal.get("seed"))
       == (REFUSED_FORMAT, exp["dataset"], exp["arm"], exp["scheme"], int(exp["seed"])))
    ok("refused: for no candidate", refusal.get("reason") == NO_CANDIDATE)
    ok("refused: an arm's result, not the default's or the twin's", exp["arm"] not in (DEFAULT, TWIN))
    # grm learn refuses only a --fold-file run for it; a drawn split's fold with no candidate gets null fits.
    ok("refused: a study run, the only one grm learn refuses for no candidate", exp["split"] == "study")
    ok("refused: the default's result at the same dataset, split and seed is scored",
       base is not None and base.get("unavailable_reason") is None and bool(base.get("folds"))
       and all(f["train"] > 0 for f in base["folds"]))
    return out


def screen_name(file):
    """The screen count grm learn writes beside a joined arm's result file."""
    return file[: -len(".pick.json")] + ".screen.json"


def check_screen(exp, doc, screen):
    """[(check name, passed, detail)] for a joined arm's screen count beside its result doc."""
    out = []

    def ok(name, cond):
        out.append((name, bool(cond), ""))

    if screen is None:
        ok("screen: the joined arm's screen count exists", False)
        return out
    ok("screen: format, revision, dataset, scheme, seed, arm and first block of its run",
       (screen.get("format"), screen.get("evaluation_revision"), screen.get("dataset"), screen.get("scheme"),
        screen.get("seed"), screen.get("arm"), screen.get("first_block"), screen.get("screened_features"))
       == (SCREEN_FORMAT, REVISION, exp["dataset"], exp["scheme"], int(exp["seed"]), exp["arm"],
           FIRST_BLOCK[exp["arm"]], SCREENED_FEATURES))
    folds = screen.get("folds") or []
    ok("screen: one row per outer fold of the result, its training genomes the fold's",
       [f.get("fold") for f in folds] == list(range(len(doc["folds"])))
       and all(f.get("train") == g["train"] for f, g in zip(folds, doc["folds"])))

    def inside(f):
        ls = f.get("learners") or []
        if [x.get("learner") for x in ls] != list(SCREENED):
            return False
        cf, ck = f["candidates_first"], f["candidates_k31"]
        for x in ls:
            fit = x["fitting"]
            if not (0 < fit <= f["train"]) or (x["learner"] == "rf") != (fit == f["train"]):
                return False
            if x["screened"] != min(SCREENED_FEATURES, cf + ck):
                return False
            if not (0 <= x["passed_first"] <= cf and 0 <= x["passed_k31"] <= ck
                    and x["passed_first"] + x["passed_k31"] == x["screened"]):
                return False
        return cf >= 0 and ck >= 0

    ok("screen: each learner's fitting genomes are its own, inside the fold's training genomes",
       all(all(0 < x["fitting"] <= f["train"] and (x["learner"] == "rf") == (x["fitting"] == f["train"])
               for x in f.get("learners") or []) for f in folds))
    ok("screen: the forest, boosted trees and lasso in order, each block's passed count inside its candidates",
       bool(folds) and all(inside(f) for f in folds))
    return out


def check(exp, doc, base, inputs, stray_inputs, release=None, served_learner=None, draw_seeds=(), twin_seeds=(),
          refused=False, screen=None):
    """[(check name, passed, detail)] for one result.

    base is the default's result at the same dataset, split and seed (None
    when this is the default's own); inputs is the arm's inputs file (None
    when it has none); stray_inputs counts inputs files in the default's or
    the twin's folder for this dataset. With refused, doc is the refusal
    record that stands for the arm's result.
    """
    arm = exp["arm"]
    e9 = dict(exp, inner=0)
    is_default = arm == DEFAULT
    if refused:
        out = check_refused(exp, doc, base)
    else:
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
    if arm not in (DEFAULT, TWIN) and not refused:
        ok("the folds are the default's at the same dataset, split and seed, by class",
           base is not None and fold_shape(base) == fold_shape(doc))
    if arm not in (DEFAULT, TWIN):
        tag = "drawn" if exp["split"] == "lineage" else exp["scheme"]
        good = inputs is not None
        if good:
            ok("inputs: format, revision, dataset, arm and tag",
               (inputs.get("format"), inputs.get("evaluation_revision"), inputs.get("dataset"), inputs.get("arm"),
                inputs.get("tag")) == (INPUT_FORMAT, REVISION, exp["dataset"], arm, tag))
            joined = arm in JOINED
            if arm in PROGRAM_ARMS:
                ok("inputs: a table read for every cohort genome", inputs["tables_read"] == inputs["samples"] > 0)
                ok("inputs: the program names seen are those out of range, collapsed and kept",
                   inputs["program_names_seen"] == inputs["program_names_out_of_range"]
                   + inputs["program_names_collapsed"] + inputs["program_rows"])
                dropped = inputs["k31_rows_as_program"]
                ok("inputs: the arm's rows are its blocks', each once",
                   inputs["arm_rows"] == inputs["program_rows"] + (inputs["k31_rows"] - dropped if joined else 0)
                   and 0 <= dropped <= (min(inputs["program_rows"], inputs["k31_rows"]) if joined else 0))
            else:
                ok("inputs: a matrix arm names its build and holds no program field",
                   inputs.get("added_build") == BUILD_OF[arm] and not any(k in inputs for k in PROGRAM_FIELDS)
                   and inputs["samples"] > 0)
                added = inputs.get("added_rows", -1)
                ok("inputs: the added rows are at most the build's species patterns",
                   0 <= added <= inputs.get("added_species_patterns", -1))
                dropped = inputs.get("k31_rows_as_added")
                ok("inputs: the arm's rows are its blocks', each once",
                   (dropped is not None) == joined
                   and inputs["arm_rows"] == added + (inputs["k31_rows"] - dropped if joined else 0)
                   and (not joined or 0 <= dropped <= min(added, inputs["k31_rows"])))
                pairs = (inputs.get("unitig_pairs_held"), inputs.get("unitig_pairs_in_part"))
                ok("inputs: the unitig arm, and only it, counts the pairs held whole and in part",
                   all(isinstance(x, int) and x >= 0 for x in pairs) if arm == "unitigs" else pairs == (None, None))
        else:
            ok("inputs: the arm's file exists", False)
    if arm in JOINED and not refused and doc.get("unavailable_reason") is None:
        out.extend(check_screen(exp, doc, screen))
    return out


# ---------- scores ----------

def score_rows(expect, docs, names, refused=None):
    """draws.tsv's rows: per result and learner. refused maps a refused arm result's file to the default's
    result at its dataset, split and seed, on whose folds it is scored as the constant model."""
    refused = refused or {}
    rows = []
    for e in expect:
        sp, ab = names[e["dataset"]]
        base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "split": e["split"],
                "seed": int(e["seed"]), "role": e["role"]}
        if e["file"] in refused:
            c = measures(constant_calls(refused[e["file"]]))
            for learner in LEARNERS:
                rows.append(dict(base, learner=learner, scored=1, constant=1,
                                 why="constant model: refused for no candidate", **c))
            continue
        doc = docs[e["file"]]
        for learner in LEARNERS:
            if doc.get("unavailable_reason") is not None:
                rows.append(dict(base, learner=learner, scored=0, constant=0, why="split unavailable"))
            else:
                rows.append(dict(base, learner=learner, scored=1, constant=0, why="", **measures(fixed(doc, learner))))
    return rows


def dataset_rows(draws):
    """datasets.tsv's rows: per dataset, split, arm and learner, the mean over its draws; unscored if any draw is."""
    groups = {}
    for r in draws:
        groups.setdefault((r["dataset"], r["split"], r["arm"], r["learner"]), []).append(r)
    out = []
    for (d, s, a, learner), rs in sorted(groups.items()):
        row = {"dataset": d, "species_slug": d.split("__")[0], "species": rs[0]["species"], "antibiotic": rs[0]["antibiotic"],
               "split": s, "arm": a, "learner": learner, "draws": len(rs), "scored_draws": sum(r["scored"] for r in rs),
               "constant_draws": sum(r.get("constant", 0) for r in rs)}
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


def count(rows, withdrawn=(), keep_constant=True):
    """count.tsv's rows from several species' datasets.tsv rows.

    withdrawn names arms withdrawn before they ran: each keeps its five
    lineage pairs in the family at p = 1 and may hold no row. With
    keep_constant False, a dataset with a constant-model draw is left out
    (count-no-constant.tsv, which decides nothing).
    """
    withdrawn = set(withdrawn)
    if not withdrawn <= set(FAMILY) - {TWIN}:
        raise SystemExit(f"only an arm of the family other than the twin can be withdrawn, not {sorted(withdrawn)}")
    held = sorted(withdrawn & {r["arm"] for r in rows})
    if held:
        raise SystemExit(f"a withdrawn arm has results: {', '.join(held)}")
    have = {r["arm"] for r in rows if r["split"] == "lineage"}
    missing = [a for a in FAMILY if a not in have and a not in withdrawn]
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
                    diffs, sps, unscored, constant = [], [], 0, 0
                    for d in ds:
                        a, b = cell[(d, split, arm, learner)], cell.get((d, split, DEFAULT, learner))
                        is_constant = int(a.get("constant_draws") or 0) > 0
                        if not keep_constant and is_constant:
                            unscored += 1
                            continue
                        if b is None or not int(a["scored"]) or not int(b["scored"]) or a[m] in (None, "") or b[m] in (None, ""):
                            unscored += 1
                            continue
                        diffs.append(r12(float(a[m]) - float(b[m])))
                        sps.append(a["species_slug"])
                        constant += is_constant
                    # Every difference, median and interval is the arm's value less the default's;
                    # only better and worse turn round where lower is better.
                    up = -1 if m in LOWER_IS_BETTER else 1
                    better = sum(1 for x in diffs if up * x > MARGIN)
                    worse = sum(1 for x in diffs if up * x < -MARGIN)
                    row = {"split": split, "pair": pair, "learner": learner, "measure": m,
                           "withdrawn": int(arm in withdrawn), "datasets": len(ds),
                           "entered": len(diffs), "not_scored": unscored, "constant": constant, "better": better,
                           "worse": worse, "same": len(diffs) - better - worse,
                           "median_difference": r12(float(np.median(diffs))) if diffs else None,
                           "species": len(set(sps)), "interval_lo": None, "interval_hi": None,
                           "sign_p": 1.0 if arm in withdrawn else C9.sign_p(better, worse), "holm_p": None,
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
    """conditions.tsv's rows: per arm, the plan's change conditions as the count reads them, two ways.

    Ruling 15 (it decides): a learner carries the change when it meets 1
    (its lineage median difference's interval lies wholly above 0) and 2
    (better beats worse by the one-sided sign test at p < 0.05 after Holm
    over 30) both, and at least 3 of 5 must carry it; 3 holds when the twin
    meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each
    read apart (the plan: the twin meets neither 1 nor 2, so a twin that
    meets 1 in three learners and 2 in none still voids the block); 4 holds
    for all five; 5 holds for each learner that carries. The stricter
    reading (it decides nothing): 1 for at least 3 of 5, 2 and 5 for all
    five, 3 the twin meeting neither of those, 4 for all five. 6 is the
    maintainer's.
    """
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in count_rows}
    sign = lambda x: None if x in (None, "") else (float(x) > 0) - (float(x) < 0)

    def lin(arm, learner, measure=PRIMARY):
        return at[("lineage", f"{arm}-{DEFAULT}", learner, measure)]

    def c1(arm):
        return {l for l in LEARNERS if lin(arm, l)["interval_lo"] not in (None, "") and float(lin(arm, l)["interval_lo"]) > 0}

    def c2(arm):
        return {l for l in LEARNERS if lin(arm, l)["holm_p"] not in (None, "") and float(lin(arm, l)["holm_p"]) < HOLM_ALPHA
                and int(lin(arm, l)["better"]) > int(lin(arm, l)["worse"])}

    def same_sign(arm, learners):
        out = set()
        for l in learners:
            study = at.get(("study", f"{arm}-{DEFAULT}", l, PRIMARY))
            s_lin = sign(lin(arm, l)["median_difference"])
            if study is not None and s_lin not in (None, 0) and sign(study["median_difference"]) == s_lin:
                out.add(l)
        return out

    twin_one, twin_two = c1(TWIN), c2(TWIN)
    twin_most = max(len(twin_one), len(twin_two))
    twin_strict = len(twin_one) >= TOUCHED_AT_LEAST or len(twin_two) == len(LEARNERS)
    out = []
    for arm in [a for a in FAMILY if a != TWIN]:
        one, two = c1(arm), c2(arm)
        carry = one & two
        low = [l for l in LEARNERS if lin(arm, l, "recall")["interval_hi"] not in (None, "")
               and float(lin(arm, l, "recall")["interval_hi"]) < -MARGIN]
        carried = same_sign(arm, carry)
        all_five = same_sign(arm, LEARNERS)
        rows = [
            ("ruling 15", "1 and 2", "a learner carries when its lineage median difference's interval lies wholly above 0 "
             "and better beats worse by the one-sided sign test at p < 0.05 after Holm over 30; at least 3 of 5 carry",
             len(carry) >= TOUCHED_AT_LEAST, len(carry)),
            ("ruling 15", "3", "the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read "
             "apart (learners: the larger of the twin's two counts)", twin_most < TOUCHED_AT_LEAST, twin_most),
            ("ruling 15", "4", "no learner's recall interval lies wholly below -0.02, for all five", not low, 5 - len(low)),
            ("ruling 15", "5", "each learner that carries has a study median of its lineage median's sign",
             bool(carry) and carried == carry, len(carried)),
            ("ruling 15", "6", "the maintainer rules the change", "", ""),
            ("stricter", "1", "the lineage median difference's interval lies wholly above 0, for at least 3 of 5 learners",
             len(one) >= TOUCHED_AT_LEAST, len(one)),
            ("stricter", "2", "better beats worse by the one-sided sign test at p < 0.05 after Holm over 30, for all five",
             len(two) == len(LEARNERS), len(two)),
            ("stricter", "3", "the twin meets neither 1 nor 2 as this reading reads them", not twin_strict, ""),
            ("stricter", "4", "no learner's recall interval lies wholly below -0.02, for all five", not low, 5 - len(low)),
            ("stricter", "5", "each of the five learners has a study median of its lineage median's sign",
             len(all_five) == len(LEARNERS), len(all_five)),
            ("stricter", "6", "the maintainer rules the change", "", ""),
        ]
        out += [{"arm": arm, "reading": reading, "condition": n, "reads": text,
                 "met": "" if m == "" else int(bool(m)), "learners": k} for reading, n, text, m, k in rows]
    return out


# ---------- files ----------

DRAW_COLUMNS = ["dataset", "species", "antibiotic", "arm", "split", "seed", "role", "learner", "scored", "constant", "why", "n",
                "resistant", "susceptible", "tp", "fp", "tn", "fn"] + [c for m in ("recall", "specificity", "vme", "me")
                                                                         for c in (m, f"{m}_lo", f"{m}_hi")] + [
    "balanced_accuracy", "log_loss"]
DATASET_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "split", "arm", "learner", "draws", "scored_draws",
                   "constant_draws", "scored", "n", "tp", "fp", "tn", "fn"] + MEASURES
COUNT_COLUMNS = ["split", "pair", "learner", "measure", "orientation", "withdrawn", "datasets", "entered", "not_scored",
                 "constant", "better", "worse", "same", "median_difference", "species", "interval_lo", "interval_hi",
                 "sign_p", "holm_p"]
CONDITION_COLUMNS = ["arm", "reading", "condition", "reads", "met", "learners"]


def load(path):
    with open(path) as f:
        return json.load(f)


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            print(__doc__.strip().splitlines()[3], file=sys.stderr)
            return 2
        out = Path(argv[1])
        rest, withdrawn = argv[2:], []
        while rest and rest[0] == "--withdrawn":
            if len(rest) < 2:
                raise SystemExit("--withdrawn takes an arm")
            withdrawn.append(rest[1])
            rest = rest[2:]
        if not rest:
            print(__doc__.strip().splitlines()[3], file=sys.stderr)
            return 2
        rows = [r for p in rest for r in read_tsv(p)]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise SystemExit("a dataset, split, arm and learner is listed twice")
        counted = count(rows, withdrawn)
        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "count.tsv", counted, COUNT_COLUMNS)
        write_tsv(out / "count-no-constant.tsv", count(rows, withdrawn, keep_constant=False), COUNT_COLUMNS)
        write_tsv(out / "conditions.tsv", conditions(counted), CONDITION_COLUMNS)
        print(f"count: {len(rows)} dataset rows from {len(rest)} files")
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
    docs, names, checks, refusals = {}, {}, [], {}
    for e in expect:
        path = picks / e["file"]
        if path.is_file():
            docs[e["file"]] = load(path)
        elif (picks / refused_name(e["file"])).is_file():
            refusals[e["file"]] = load(picks / refused_name(e["file"]))
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        if rel_path.is_file() and e["dataset"] not in names:
            rel = load(rel_path)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in expect if e["arm"] == DEFAULT}
    constant = {}
    for e in expect:
        is_refused = e["file"] in refusals
        if e["file"] not in docs and not is_refused:
            checks.append((e["file"], "exists", False, ""))
            continue
        checks.append((e["file"], "exists" if not is_refused else "exists, as a refusal for no candidate", True, ""))
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
        doc = refusals[e["file"]] if is_refused else docs[e["file"]]
        screen_path = picks / screen_name(e["file"])
        screen = load(screen_path) if e["arm"] in JOINED and screen_path.is_file() else None
        for name, passed, detail in check(e, doc, base, inputs, stray, rel, served_by.get(e["dataset"]),
                                          draw_seeds, twin_seeds, refused=is_refused, screen=screen):
            checks.append((e["file"], name, passed, detail))
        if is_refused:
            constant[e["file"]] = base
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "checks.tsv", [{"file": f, "check": c, "passed": p, "detail": d} for f, c, p, d in checks],
              ["file", "check", "passed", "detail"])
    failed = [c for c in checks if not c[2]]
    unavailable = sum(1 for d in docs.values() if d.get("unavailable_reason") is not None)
    print(f"checks: {len(checks) - len(failed)} of {len(checks)} pass; results {len(docs)} of {len(expect)}, "
          f"refused for no candidate {len(refusals)}; lineage or study split unavailable in {unavailable}")
    if failed:
        return 1
    draws = score_rows(expect, docs, names, constant)
    write_tsv(out / "draws.tsv", draws, DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", dataset_rows(draws), DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
