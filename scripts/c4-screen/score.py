"""The scorer of plan block C4 (docs/screen-arms.md, measures and files).

  python3 -I score.py <expect.tsv> <picks dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> [--withdrawn <arm> ...] <datasets.tsv> [<datasets.tsv> ...]

The first form checks every result expect.py lists, then scores it. The
default's and the twin's results take scripts/c3-inputs/score.py's checks,
imported: C9's checks of every result, the default at seed 42 held to the
release report, no inputs file, and the twin at seeds of its own. An arm's
result takes C9's checks, its folds must be the default's at the same
dataset, split and seed by class, and a refusal for no candidate takes C3's
check of ruling 16 (it cannot arise here: every arm screens the default's
candidates or more). Beside them, per result:

- its scheme is its split's: lineage, or the study split's study-whole;
- an arm's inputs file has the format, revision, dataset, tag, seed, arm,
  rule and budget of its run and the floor of the arm (3 for floor-3,
  max(3, n/200) at the genomes it fits for the screen arms), the result's
  floor too; its outer folds of the result's split have the result's
  training genomes; rf, gbt and lasso screen once in each outer fold of
  the run, each by the arm's rule, on its fold's candidates, rf on the
  fold's training genomes and gbt and lasso on no more, each keeping the
  budget or every candidate if fewer; the final models of a drawn run
  screen once each on the dataset's rows, rf on every genome; the adjusted
  screens name their strata, no fewer than the informative ones, and the
  chi-square screens none; each learner of each outer fold has a fit time;
- a screen arm's scm and cart, which read no screen, make the default's
  calls in every fold;
- where the dataset has fewer than 800 genomes, so max(3, n/200) is 3, the
  floor arm's result is the default's apart from its timing field.

Writes checks.tsv, draws.tsv and datasets.tsv as C3's scorer does (scores
per learner, the constant model of ruling 16 for a refusal), with each
dataset's genome count in datasets.tsv. A dataset, split and arm is scored
only with every draw expect.py lists for it present. Prints the check
counts and, when a check fails, the names of the failed checks with their
files (never a score), and exits 1, scoring nothing.

The second form writes count.tsv, count-no-constant.tsv, conditions.tsv and
count-floor-all.tsv over several species' datasets.tsv, and prints the
row count only. The lead's rulings of 23 Nov 2026 (questions 10 and 12) fix
the family at 21 pairs on the lineage primary measure: cmh-4000 and
cmh-16000 with rf, gbt and lasso, the learners a screen touches, over every
dataset; floor-3 with all five over the datasets of 800 genomes or more,
the ones its floor can move (92 of the first run's 233, fixed by cohort
size before any score); the twin with all five over every dataset, and
again over those of 800 genomes or more. C3's count runs over each set,
with its median differences, species-resample intervals, sign tests and
refusal of part of the family, and Holm runs over the 21. count.tsv's
column `set` names the datasets a row reads: `all` or `800-or-more`. A
screen arm's scm and cart are checked equal to the default's and not
counted.

conditions.tsv reads the plan's conditions per arm. A screen arm, over
every dataset: 1, 2 and 5 for each of rf, gbt and lasso, 4 for the three,
3 the twin over every dataset meeting 1 in fewer than 3 of 5 learners and 2
in fewer than 3 of 5, each read apart. The floor arm, over the datasets of
800 genomes or more: C3's two readings, imported (ruling 15 decides, the
stricter beside it), at Holm over 21, with the twin read over the same
datasets for 3. 6 is the maintainer's. count-floor-all.tsv is the floor
arm's count over every dataset, where those under 800 genomes enter as a
difference of 0; it stands beside and decides nothing.
"""

import contextlib
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _load("c3_score_for_c4", HERE.parent / "c3-inputs" / "score.py")
LISTS = _load("c4_lists_for_score", HERE / "lists.py")
EXPECT = _load("c4_expect_for_score", HERE / "expect.py")
C9 = C3.C9
LEARNERS, MEASURES, PRIMARY, MARGIN, HOLM_ALPHA = C3.LEARNERS, C3.MEASURES, C3.PRIMARY, C3.MARGIN, C3.HOLM_ALPHA
REVISION, RELEASE_SEED, DEFAULT, TWIN = C3.REVISION, C3.RELEASE_SEED, C3.DEFAULT, C3.TWIN
TOUCHED_AT_LEAST = C3.TOUCHED_AT_LEAST
read_tsv, write_tsv, load = C3.read_tsv, C3.write_tsv, C3.load
STUDY_SCHEME = EXPECT.STUDY_SCHEME
INPUT_FORMAT = "grm-c4-screen-1"
SCREEN_ARMS = ("cmh-4000", "cmh-16000")
FLOOR_ARM = "floor-3"
ARMS = SCREEN_ARMS + (FLOOR_ARM,)
assert ARMS == EXPECT.ARMS
FAMILY = list(ARMS) + [TWIN]
# The arm name an inputs file carries, its rule and its budget.
ARM_RULE = {"cmh-4000": ("cmh", "cmh", 4000), "cmh-16000": ("cmh", "cmh", 16000),
            FLOOR_ARM: ("floor", "chi-square", 4000)}
SCREENED = ["rf", "gbt", "lasso"]
TIMING = ("seconds_inner",)
DATASET_COLUMNS = C3.DATASET_COLUMNS[:4] + ["genomes"] + C3.DATASET_COLUMNS[4:]


def touched(arm):
    """The learners an arm's change reaches: the screened three for a screen arm, all five otherwise."""
    return SCREENED if arm in SCREEN_ARMS else LEARNERS


def floor_moves(genomes):
    return LISTS.default_floor(int(genomes)) != LISTS.FLOOR


def untimed(doc):
    return {k: v for k, v in doc.items() if k not in TIMING}


def same_calls(base, doc, learners):
    """Whether `learners` make the same outer calls and loss in every fold of the two results."""
    if len(base.get("folds", [])) != len(doc.get("folds", [])):
        return False
    for fb, fd in zip(base.get("folds", []), doc.get("folds", [])):
        cb = {o["learner"]: o for o in fb["outer"]}
        cd = {o["learner"]: o for o in fd["outer"]}
        if any(cb.get(l) != cd.get(l) for l in learners):
            return False
    return True


def inputs_checks(exp, doc, inputs):
    """[(check name, passed, detail)] for an arm's inputs file beside its result `doc`."""
    out = []

    def ok(name, cond):
        out.append((name, bool(cond), ""))

    if inputs is None:
        ok("inputs: the arm's file exists", False)
        return out
    arm = exp["arm"]
    name, rule, budget = ARM_RULE[arm]
    drawn = exp["split"] == "lineage"
    tag = "drawn" if drawn else exp["scheme"]
    ok("inputs: format, revision, dataset, tag, seed and arm",
       (inputs.get("format"), inputs.get("evaluation_revision"), inputs.get("dataset"), inputs.get("tag"),
        inputs.get("seed"), inputs.get("arm")) == (INPUT_FORMAT, REVISION, exp["dataset"], tag, int(exp["seed"]), name))
    ok("inputs: the arm's rule and budget", (inputs.get("rule"), inputs.get("budget")) == (rule, budget))
    samples = int(inputs.get("samples") or 0)
    floor = LISTS.FLOOR if arm == FLOOR_ARM else LISTS.default_floor(samples)
    ok("inputs: the arm's floor, and the result's", samples > 0 and inputs.get("min_support") == floor
       and doc.get("min_support") == floor)
    folds, screens, times = inputs["folds"], inputs["screens"], inputs["fit_times"]
    mine = [f for f in folds if f["split"] == exp["scheme"]]
    ok("inputs: the result's outer folds, with its training genomes",
       [f["fold"] for f in mine] == list(range(len(mine)))
       and [f["train"] for f in mine] == [f["train"] for f in doc.get("folds", [])])
    splits = {"random", "lineage", "final"} if drawn else {exp["scheme"]}
    ok("inputs: every fit is of the run's splits", all(x["split"] in splits for x in folds + screens + times))
    ok("inputs: every screen by the arm's rule", all(s["rule"] == rule for s in screens))
    ok("inputs: each screen keeps the budget, or every candidate if fewer",
       all(s["kept"] == min(budget, s["candidates"]) for s in screens))
    fold_of = {(f["split"], f["fold"]): f for f in folds}
    outer = [s for s in screens if s["split"] != "final"]
    keys = sorted((s["split"], s["fold"], s["learner"]) for s in outer)
    ok("inputs: rf, gbt and lasso screen once in each outer fold",
       keys == sorted((sp, k, l) for sp, k in fold_of for l in SCREENED))
    ok("inputs: each outer screen reads its fold's candidates, rf its training genomes, gbt and lasso no more",
       all((s["split"], s["fold"]) in fold_of and s["candidates"] == fold_of[(s["split"], s["fold"])]["candidates"]
           and (s["fitting"] == fold_of[(s["split"], s["fold"])]["train"] if s["learner"] == "rf"
                else 0 < s["fitting"] <= fold_of[(s["split"], s["fold"])]["train"]) for s in outer))
    finals = [s for s in screens if s["split"] == "final"]
    if drawn:
        ok("inputs: the final models screen once each on the dataset's rows, rf on every genome",
           sorted(s["learner"] for s in finals) == sorted(SCREENED)
           and all(s["fold"] is None and s["candidates"] == inputs["dataset_rows"]
                   and (s["fitting"] == samples if s["learner"] == "rf" else 0 < s["fitting"] <= samples)
                   for s in finals))
    else:
        ok("inputs: a fold file's run fits no final model", not finals)
    if arm in SCREEN_ARMS:
        ok("inputs: each adjusted screen names its strata, no fewer than the informative ones",
           all(s["strata"] is not None and s["informative_strata"] is not None and s["kept_at_zero"] is not None
               and 0 <= s["informative_strata"] <= s["strata"] and 0 <= s["kept_at_zero"] <= s["kept"]
               for s in screens))
    else:
        ok("inputs: the chi-square screens name no strata",
           all(s["strata"] is None and s["informative_strata"] is None and s["kept_at_zero"] is None for s in screens))
    ok("inputs: a fit time for each learner of each outer fold",
       sorted((t["split"], t["fold"], t["learner"]) for t in times)
       == sorted((sp, k, l) for sp, k in fold_of for l in LEARNERS))
    return out


def check(exp, doc, base, inputs, stray_inputs, release=None, served_learner=None, draw_seeds=(), twin_seeds=(),
          refused=False):
    """[(check name, passed, detail)] for one result, with C3's arguments: base is the default's result at the
    same dataset, split and seed (None for the default's own), inputs the arm's inputs file, stray_inputs the
    inputs files in the default's or the twin's folder for this dataset, and with refused, doc is the refusal
    record that stands for the arm's result."""
    arm = exp["arm"]
    if arm in (DEFAULT, TWIN):
        out = C3.check(exp, doc, base, inputs, stray_inputs, release, served_learner, draw_seeds, twin_seeds,
                       refused=refused)
    elif refused:
        out = C3.check_refused(exp, doc, base)
    else:
        out = C9.check(dict(exp, inner=0), doc, None, None)

    def ok(name, cond):
        out.append((name, bool(cond), ""))

    ok("the result's scheme is its split's",
       (exp["split"], exp["scheme"]) in (("lineage", "lineage"), ("study", STUDY_SCHEME))
       and doc.get("scheme") == exp["scheme"])
    if arm not in ARMS:
        return out
    if not refused:
        ok("the folds are the default's at the same dataset, split and seed, by class",
           base is not None and C3.fold_shape(base) == C3.fold_shape(doc))
        if arm in SCREEN_ARMS:
            ok("scm and cart, which read no screen, make the default's calls in every fold",
               base is not None and same_calls(base, doc, [l for l in LEARNERS if l not in SCREENED]))
        if arm == FLOOR_ARM and base is not None and not floor_moves(base["samples"]):
            ok("below 800 genomes the floor arm's result is the default's, timing apart", untimed(doc) == untimed(base))
        out += inputs_checks(exp, doc, inputs)
    return out


# ---------- scores ----------

def dataset_rows(draws, expect, genomes):
    """datasets.tsv's rows: C3's, after each draw expect lists and `draws` lacks is put in as an unscored row,
    so a dataset, split and arm is scored only with every draw present; with each dataset's genome count."""
    have = {(r["dataset"], r["split"], r["arm"], int(r["seed"]), r["learner"]) for r in draws}
    names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in draws}
    filled = list(draws)
    for e in expect:
        for learner in LEARNERS:
            if (e["dataset"], e["split"], e["arm"], int(e["seed"]), learner) not in have:
                sp, ab = names.get(e["dataset"], ("", ""))
                filled.append({"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"],
                               "split": e["split"], "seed": int(e["seed"]), "role": e["role"], "learner": learner,
                               "scored": 0, "constant": 0, "why": "missing"})
    rows = C3.dataset_rows(filled)
    for r in rows:
        r["genomes"] = genomes.get(r["dataset"])
    return rows


# ---------- the count ----------

@contextlib.contextmanager
def family(arms):
    """C3's count and conditions read their family from C3.FAMILY: C4's for the time of the call."""
    saved = C3.FAMILY
    C3.FAMILY = list(arms)
    try:
        yield
    finally:
        C3.FAMILY = saved


def arm_of(pair):
    return pair[: -len(DEFAULT) - 1]


ALL, MOVES = "all", "800-or-more"
COUNT_COLUMNS = ["set"] + C3.COUNT_COLUMNS
# The lineage primary measure's family (the lead's rulings 10 and 12): (arm, set) and its learners.
FAMILY_PARTS = [("cmh-4000", ALL), ("cmh-16000", ALL), (FLOOR_ARM, MOVES), (TWIN, ALL), (TWIN, MOVES)]
FAMILY_SIZE = 21


def moved_rows(rows):
    """The datasets.tsv rows of datasets of 800 genomes or more, the ones the floor arm can move."""
    out = []
    for r in rows:
        if r.get("genomes") in (None, ""):
            raise SystemExit(f"{r['dataset']}: no genome count, so its place in the floor arm's set is not known")
        if floor_moves(r["genomes"]):
            out.append(r)
    return out


def holm21(rows):
    fam = [r for part in FAMILY_PARTS for r in rows
           if (arm_of(r["pair"]), r["set"]) == part and r["split"] == "lineage" and r["measure"] == PRIMARY]
    assert len(fam) == FAMILY_SIZE == sum(len(touched(a)) for a, _ in FAMILY_PARTS), len(fam)
    for r in rows:
        r["holm_p"] = None
    for r, p in zip(fam, C9.holm([r["sign_p"] for r in fam])):
        r["holm_p"] = p
    return rows


def count(rows, withdrawn=(), keep_constant=True):
    """count.tsv's rows from several species' datasets.tsv rows: C3's count over every dataset for the screen
    arms and the twin, the learners each arm touches, and over the datasets of 800 genomes or more for the
    floor arm and the twin; Holm over the 21 lineage pairs of the primary measure."""
    with family(FAMILY):
        whole = C3.count(rows, withdrawn, keep_constant)
    whole = [dict(r, set=ALL) for r in whole
             if arm_of(r["pair"]) != FLOOR_ARM and r["learner"] in touched(arm_of(r["pair"]))]
    moved = moved_rows(rows)
    if not any(r["split"] == "lineage" and r["arm"] == TWIN for r in moved):
        raise SystemExit("no dataset of 800 genomes or more: the floor arm's part of the family cannot be read")
    with family([FLOOR_ARM, TWIN]):
        part = C3.count(moved, [a for a in withdrawn if a == FLOOR_ARM], keep_constant)
    return holm21(whole + [dict(r, set=MOVES) for r in part])


def floor_all_count(rows, withdrawn=()):
    """count-floor-all.tsv: the floor arm's rows over every dataset, those under 800 genomes at a difference
    of 0; decides nothing (the lead's ruling 10)."""
    if FLOOR_ARM in withdrawn:
        return []
    with family([FLOOR_ARM, TWIN]):
        out = C3.count(rows)
    out = [dict(r, set=ALL) for r in out if arm_of(r["pair"]) == FLOOR_ARM]
    for r in out:
        r["holm_p"] = None
    return out


def lineage_row(at, arm, learner, measure=PRIMARY):
    return at[("lineage", f"{arm}-{DEFAULT}", learner, measure)]


def meets_1(row):
    return row["interval_lo"] not in (None, "") and float(row["interval_lo"]) > 0


def meets_2(row):
    return row["holm_p"] not in (None, "") and float(row["holm_p"]) < HOLM_ALPHA and int(row["better"]) > int(row["worse"])


def sign(x):
    return None if x in (None, "") else (float(x) > 0) - (float(x) < 0)


def conditions(count_rows):
    """conditions.tsv's rows: per arm, the plan's change conditions as the count reads them."""
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in count_rows if r["set"] == ALL}
    twin_most = max(sum(meets_1(lineage_row(at, TWIN, l)) for l in LEARNERS),
                    sum(meets_2(lineage_row(at, TWIN, l)) for l in LEARNERS))
    out = []
    for arm in SCREEN_ARMS:
        one = [l for l in SCREENED if meets_1(lineage_row(at, arm, l))]
        two = [l for l in SCREENED if meets_2(lineage_row(at, arm, l))]
        low = [l for l in SCREENED if lineage_row(at, arm, l, "recall")["interval_hi"] not in (None, "")
               and float(lineage_row(at, arm, l, "recall")["interval_hi"]) < -MARGIN]
        five = []
        for l in SCREENED:
            study = at.get(("study", f"{arm}-{DEFAULT}", l, PRIMARY))
            s_lin = sign(lineage_row(at, arm, l)["median_difference"])
            if study is not None and s_lin not in (None, 0) and sign(study["median_difference"]) == s_lin:
                five.append(l)
        rows = [
            ("1", "the lineage median difference's interval lies wholly above 0, for each of rf, gbt and lasso",
             len(one) == len(SCREENED), len(one)),
            ("2", "better beats worse by the one-sided sign test at p < 0.05 after Holm over 21, for each of the three",
             len(two) == len(SCREENED), len(two)),
            ("3", "the twin, over every dataset, meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, "
             "each read apart (learners: the larger of the twin's two counts)", twin_most < TOUCHED_AT_LEAST,
             twin_most),
            ("4", "no recall interval of the three lies wholly below -0.02", not low, len(SCREENED) - len(low)),
            ("5", "each of the three has a study median of its lineage median's sign", len(five) == len(SCREENED),
             len(five)),
            ("6", "the maintainer rules the change", "", ""),
        ]
        out += [{"arm": arm, "reading": "the learners it touches, every dataset", "condition": n, "reads": text,
                 "met": "" if m == "" else int(bool(m)), "learners": k} for n, text, m, k in rows]
    # The floor arm and its twin over the datasets of 800 genomes or more (the lead's ruling 10).
    with family([FLOOR_ARM, TWIN]):
        floor = C3.conditions([r for r in count_rows if r["set"] == MOVES])
    for r in floor:
        r["reading"] = f"{r['reading']}, 800 genomes or more"
        # C3 writes its own family's size (30, or 10 inside family([FLOOR_ARM, TWIN])); C4's Holm runs over 21.
        r["reads"] = re.sub(r"Holm over \d+", "Holm over 21", r["reads"])
    return out + floor


# ---------- files ----------

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
        write_tsv(out / "conditions.tsv", conditions(counted), C3.CONDITION_COLUMNS)
        write_tsv(out / "count-floor-all.tsv", floor_all_count(rows, withdrawn), COUNT_COLUMNS)
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
    docs, refusals, names = {}, {}, {}
    for e in expect:
        path = picks / e["file"]
        if path.is_file():
            docs[e["file"]] = load(path)
        elif (picks / C3.refused_name(e["file"])).is_file():
            refusals[e["file"]] = load(picks / C3.refused_name(e["file"]))
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        if e["dataset"] not in names and rel_path.is_file():
            rel = load(rel_path)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in expect if e["arm"] == DEFAULT}
    checks, constant = [], {}
    for e in expect:
        is_refused = e["file"] in refusals
        if e["file"] not in docs and not is_refused:
            checks.append((e["file"], "exists", False, ""))
            continue
        checks.append((e["file"], "exists, as a refusal for no candidate" if is_refused else "exists", True, ""))
        checks.append((e["file"], "release report", e["dataset"] in names, ""))
        rel = None
        if e["arm"] == DEFAULT and int(e["seed"]) == RELEASE_SEED and e["dataset"] in names:
            rel = load(release / "datasets" / f"{e['dataset']}.json")
        base = docs.get(base_of.get((e["dataset"], e["split"], e["seed"]))) if e["arm"] != DEFAULT else None
        inputs = load(picks / e["inputs"]) if e["inputs"] != "-" and (picks / e["inputs"]).is_file() else None
        stray = 0
        if e["arm"] in (DEFAULT, TWIN):
            stray = sum(1 for _ in (picks / e["arm"]).glob(f"{e['dataset']}.*.inputs.json"))
        doc = refusals[e["file"]] if is_refused else docs[e["file"]]
        for name, passed, detail in check(e, doc, base, inputs, stray, rel, served_by.get(e["dataset"]),
                                          draw_seeds, twin_seeds, refused=is_refused):
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
        for f, c, _, _ in failed[:40]:
            print(f"failed: {f}: {c}", file=sys.stderr)
        return 1
    genomes = {e["dataset"]: docs[e["file"]]["samples"] for e in expect
               if e["arm"] == DEFAULT and e["split"] == "lineage" and int(e["seed"]) == RELEASE_SEED}
    draws = C3.score_rows(expect, docs, names, constant)
    write_tsv(out / "draws.tsv", draws, C3.DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", dataset_rows(draws, expect, genomes), DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
