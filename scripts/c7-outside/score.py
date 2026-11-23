"""The scorer of plan block C7 (docs/outside-learners.md, measures and files).

  python3 -I score.py <expect.tsv> <run dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> [--withdrawn <arm> ...] <datasets.tsv> [<datasets.tsv> ...]

The first form checks every result expect.py lists, then scores it. The
default (k31) and the twin are grm learn --exp-pick results the C3 stage of
the run wrote (<run dir>/c3/results/picks/...), checked as C3's scorer checks
them (scripts/c3-inputs/score.py), the default's seed-42 lineage figures
held to the release report's. Each library arm's result is a
grm-c7-outside-1 report (<run dir>/results/outside/<arm>/...), checked here:

- its format, revision, dataset, split, seed, arm, library and version, and
  the GRM learner it pairs with, which is the pairing's (xgboost and
  lightgbm with gbt, elasticnet with lasso);
- it holds only the listed keys, so no genome, row, feature or model field;
- its folds are the default's at the same dataset, split and seed: every
  fold's training and test genomes by class, and the same split unavailable;
- each fold's calls are its test genomes'.

Scores, per dataset, split, seed and row: the arm's calls summed over the
outer folds, and the default's and the twin's calls of the learner each arm
pairs with (gbt and lasso), with balanced accuracy, recall, specificity,
the very major and the major error rate (each with its count, denominator
and Wilson 95% interval) and log loss (scripts/c9-pick/score.py, measures).
A draw that is missing is a failed check and leaves its dataset unscored,
never a zero. A library arm's fold with no feature, or the elastic net's
null fit, was scored in the report at the fold's training share of
resistant genomes (ruling 16's constant model); its draws.tsv rows say so
in `constant`. Writes checks.tsv, draws.tsv and datasets.tsv (C3's
dataset_rows: per dataset, split, arm and learner, the mean over the draws;
the default again as k31-at-arm-seeds, over the draws at the seeds the
library arms ran, which is what each arm's pair reads, so an arm and the
default meet on the same draws; the twin's pairs read every draw)
and prints the check counts only. Exits 1, scoring nothing, when a check
fails.

The second form writes count.tsv, count-no-constant.tsv,
count-no-unconverged.tsv and conditions.tsv over several species'
datasets.tsv and prints the row count only. It refuses unless the three arms
and the twin all have lineage rows or the arm is named with --withdrawn: the
family is fixed at 5 pairs and the count is taken once. An arm withdrawn
before it runs keeps its pair at p = 1 and is named in count.tsv's
`withdrawn`; the count refuses a withdrawn arm that has results, and the twin
or the default named as withdrawn (the lead's addition, 23 Nov).
count-no-unconverged.tsv leaves out the datasets whose arm draws hold a
refit that reached saga's epoch cap, and decides nothing (question 9). Pairs, each an arm's learner against
the default's same learner on the arm's own datasets: xgboost-k31 (gbt),
lightgbm-k31 (gbt), elasticnet-k31 (lasso), and on the lineage split the
twin's twin-k31 (gbt) and twin-k31 (lasso). Per split, pair and measure:
better, worse and the same by more than MARGIN, the median difference (the
arm's value less the default's) and its interval from RESAMPLES species
resamples, the one-sided sign test, and Holm over the 5 lineage rows of the
primary measure, the twin's two included (question 10). conditions.tsv
reads, per arm, conditions 1, 2, 4 and 5 on its one pair, and condition 3 on
both of the twin's pairs: it holds when the twin meets 1 in neither pair and
2 in neither, each read apart (question 14). Beside it, in `decides` 0, each
arm's condition 3 on its own learner's twin pair alone, which decides
nothing.
"""

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _load("c3_score_for_c7", HERE.parent / "c3-inputs" / "score.py")
O = _load("c7_outside_for_score", HERE / "outside.py")
C9 = C3.C9
MEASURES, LOWER_IS_BETTER, PRIMARY = C9.MEASURES, C9.LOWER_IS_BETTER, C9.PRIMARY
MARGIN, RESAMPLES, LO_INDEX, HI_INDEX, HOLM_ALPHA = C9.MARGIN, C9.RESAMPLES, C9.LO_INDEX, C9.HI_INDEX, C9.HOLM_ALPHA
r12, measures, fixed, read_tsv, write_tsv = C9.r12, C9.measures, C9.fixed, C9.read_tsv, C9.write_tsv
DEFAULT, TWIN = C3.DEFAULT, C3.TWIN
ARMS = O.ARMS
PAIRED = sorted(set(ARMS.values()))
FAMILY = [(arm, ARMS[arm]) for arm in ("xgboost", "lightgbm", "elasticnet")] + [(TWIN, learner) for learner in PAIRED]
# The default read at the library arms' seeds only: an arm's pair meets the default on the same draws (the arms
# run at the default's three seeds, 42, 43 and 44, since the lead's ruling on question 1). The twin meets every draw.
ARM_BASE = "k31-at-arm-seeds"


# ---------- checks ----------

def check_outside(exp, doc, base):
    """(name, passed, detail) for one library arm's report against its expect row and the default's result."""
    out = []
    keys_ok = True
    try:
        O.check_report(doc)
    except O.Refuse as e:
        keys_ok = False
        out.append(("listed keys only", False, str(e)))
    if keys_ok:
        out.append(("listed keys only", True, ""))
    arm = exp["arm"]
    out.append(("format and revision", doc.get("format") == O.FORMAT and doc.get("evaluation_revision") == O.REVISION, ""))
    out.append(("dataset, split and seed", doc.get("dataset") == exp["dataset"] and doc.get("scheme") == exp["scheme"]
                and doc.get("seed") == int(exp["seed"]), ""))
    out.append(("arm and library", doc.get("arm") == arm and doc.get("library") == O.LIBRARY.get(arm)
                and doc.get("version") == O.VERSIONS.get(O.LIBRARY.get(arm)), ""))
    out.append(("pairs with its GRM learner", doc.get("pairs_with") == ARMS.get(arm), f"{doc.get('pairs_with')}"))
    same = base is not None and keys_ok and C3.fold_shape(base) == C3.fold_shape(doc)
    out.append(("folds are the default's", same, ""))
    calls_ok = keys_ok and all(f["calls"]["n"] == f["test"] and f["calls"]["tp"] + f["calls"]["fn"] == f["test_resistant"]
                               and f["calls"]["tp"] + f["calls"]["fp"] + f["calls"]["tn"] + f["calls"]["fn"] == f["test"]
                               for f in doc.get("folds", []))
    out.append(("calls are the test genomes'", calls_ok, ""))
    return out


def c3_row(e):
    """The C3 checker's expect row for a default or twin result."""
    return {"dataset": e["dataset"], "species_slug": e["species_slug"], "arm": e["arm"], "split": e["split"],
            "scheme": e["scheme"], "seed": e["seed"], "role": e["role"], "file": e["file"], "inputs": "-"}


def outside_calls(doc):
    c = dict(C9.ZERO)
    for f in doc["folds"]:
        c = C9.add(c, f["calls"])
    return c


def score_rows(expect, docs, names):
    """draws.tsv's rows: per result, the learners the pairs read. A missing result's rows are unscored."""
    rows = []
    for e in expect:
        sp, ab = names.get(e["dataset"], ("", ""))
        base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "split": e["split"],
                "seed": int(e["seed"]), "role": e["role"]}
        learners = PAIRED if e["arm"] in (DEFAULT, TWIN) else [ARMS[e["arm"]]]
        doc = docs.get(e["file"])
        for learner in learners:
            if doc is None:
                rows.append(dict(base, learner=learner, scored=0, constant=0, unconverged=0, why="draw missing"))
            elif doc.get("unavailable_reason") is not None:
                rows.append(dict(base, learner=learner, scored=0, constant=0, unconverged=0, why="split unavailable"))
            elif e["arm"] in (DEFAULT, TWIN):
                rows.append(dict(base, learner=learner, scored=1, constant=0, unconverged=0, why="",
                                 **measures(fixed(doc, learner))))
            else:
                constant = int(doc["constant_folds"] > 0)
                rows.append(dict(base, learner=learner, scored=1, constant=constant,
                                 unconverged=int(doc["unconverged_refits"] > 0),
                                 why="a fold scored at its training share" if constant else "",
                                 **measures(outside_calls(doc))))
    return rows


def dataset_rows(draws):
    """datasets.tsv's rows (C3's dataset_rows), with the default also read at the library arms' seeds as ARM_BASE."""
    arm_seeds = {(r["split"], r["seed"]) for r in draws if r["arm"] in ARMS}
    again = [dict(r, arm=ARM_BASE) for r in draws if r["arm"] == DEFAULT and (r["split"], r["seed"]) in arm_seeds]
    rows = C3.dataset_rows(draws + again)
    unconverged = {}
    for r in draws + again:
        k = (r["dataset"], r["split"], r["arm"], r["learner"])
        unconverged[k] = unconverged.get(k, 0) + int(r.get("unconverged", 0))
    for r in rows:
        r["unconverged_draws"] = unconverged[(r["dataset"], r["split"], r["arm"], r["learner"])]
    return rows


# ---------- the count ----------

def interval(diffs, sps, seed):
    """The species resample interval of the median difference."""
    if len(set(diffs)) == 1:
        return diffs[0], diffs[0]
    by = {}
    for s, x in zip(sps, diffs):
        by.setdefault(s, []).append(x)
    keys = sorted(by)
    rng = np.random.default_rng(seed)
    meds = []
    for _ in range(RESAMPLES):
        take = rng.integers(0, len(keys), len(keys))
        meds.append(r12(float(np.median(np.concatenate([by[keys[i]] for i in take])))))
    meds.sort()
    return meds[LO_INDEX], meds[HI_INDEX]


def count(rows, withdrawn=(), keep_constant=True, keep_unconverged=True):
    """count.tsv's rows from several species' datasets.tsv rows.

    withdrawn names library arms withdrawn before they ran: each keeps its pair in the family at p = 1 and holds
    no row. With keep_constant False, a dataset with a constant-model draw is left out (count-no-constant.tsv); with
    keep_unconverged False, a dataset with an unconverged refit (count-no-unconverged.tsv). Neither decides.
    """
    withdrawn = set(withdrawn)
    if not withdrawn <= set(ARMS):
        raise SystemExit(f"only a library arm ({', '.join(sorted(ARMS))}) can be withdrawn, not the twin or the "
                         f"default: {', '.join(sorted(withdrawn - set(ARMS)))}")
    held = sorted(withdrawn & {r["arm"] for r in rows})
    if held:
        raise SystemExit(f"a withdrawn arm has results: {', '.join(held)}")
    have = {(r["arm"], r["learner"]) for r in rows if r["split"] == "lineage"}
    missing = [f"{a} ({l})" for a, l in FAMILY if (a, l) not in have and a not in withdrawn]
    if missing:
        raise SystemExit(f"the family is not complete: no lineage result for {', '.join(missing)}; "
                         "the count is taken once, over every arm and the twin")
    cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}
    out = []
    for split, pairs in (("lineage", FAMILY), ("study", [p for p in FAMILY if p[0] != TWIN])):
        for arm, learner in pairs:
            pair = f"{arm}-{DEFAULT}"
            base_arm = DEFAULT if arm == TWIN else ARM_BASE
            for m in MEASURES:
                ds = sorted({d for d, s, a, l in cell if s == split and a == arm and l == learner})
                diffs, sps, unscored, constant, unconverged = [], [], 0, 0, 0
                for d in ds:
                    a, b = cell[(d, split, arm, learner)], cell.get((d, split, base_arm, learner))
                    is_constant = int(a.get("constant_draws") or 0) > 0
                    is_unconverged = int(a.get("unconverged_draws") or 0) > 0
                    if (not keep_constant and is_constant) or (not keep_unconverged and is_unconverged):
                        unscored += 1
                        continue
                    if b is None or not int(a["scored"]) or not int(b["scored"]) or a[m] in (None, "") or b[m] in (None, ""):
                        unscored += 1
                        continue
                    diffs.append(r12(float(a[m]) - float(b[m])))
                    sps.append(a["species_slug"])
                    constant += is_constant
                    unconverged += is_unconverged
                up = -1 if m in LOWER_IS_BETTER else 1
                better = sum(1 for x in diffs if up * x > MARGIN)
                worse = sum(1 for x in diffs if up * x < -MARGIN)
                row = {"split": split, "pair": pair, "learner": learner, "measure": m,
                       "withdrawn": int(arm in withdrawn), "datasets": len(ds), "entered": len(diffs),
                       "not_scored": unscored, "constant": constant, "unconverged": unconverged, "better": better,
                       "worse": worse, "same": len(diffs) - better - worse,
                       "median_difference": r12(float(np.median(diffs))) if diffs else None, "species": len(set(sps)),
                       "interval_lo": None, "interval_hi": None,
                       "sign_p": 1.0 if arm in withdrawn else C9.sign_p(better, worse), "holm_p": None,
                       "orientation": "lower" if m in LOWER_IS_BETTER else "higher"}
                if diffs:
                    row["interval_lo"], row["interval_hi"] = interval(diffs, sps, C3.row_seed(split, pair, learner, m))
                out.append(row)
    family = [r for r in out if r["split"] == "lineage" and r["measure"] == PRIMARY]
    assert len(family) == len(FAMILY)
    for r, p in zip(family, C9.holm([r["sign_p"] for r in family])):
        r["holm_p"] = p
    return out


def conditions(count_rows):
    """conditions.tsv's rows: per arm, the plan's conditions 1 to 5 as this block reads them (6 is the maintainer's)."""
    at = {(r["split"], r["pair"], r["learner"], r["measure"]): r for r in count_rows}

    def c1(arm, learner):
        r = at[("lineage", f"{arm}-{DEFAULT}", learner, PRIMARY)]
        return r["interval_lo"] is not None and r["interval_lo"] > 0

    def c2(arm, learner):
        r = at[("lineage", f"{arm}-{DEFAULT}", learner, PRIMARY)]
        return r["holm_p"] is not None and r["holm_p"] < HOLM_ALPHA and r["better"] > r["worse"]

    twin_1 = [c1(TWIN, learner) for learner in PAIRED]
    twin_2 = [c2(TWIN, learner) for learner in PAIRED]
    sign = lambda x: None if x is None else (x > 0) - (x < 0)  # noqa: E731
    out = []
    for arm in ("xgboost", "lightgbm", "elasticnet"):
        learner = ARMS[arm]
        recall = at[("lineage", f"{arm}-{DEFAULT}", learner, "recall")]
        lin = at[("lineage", f"{arm}-{DEFAULT}", learner, PRIMARY)]["median_difference"]
        study = at.get(("study", f"{arm}-{DEFAULT}", learner, PRIMARY), {}).get("median_difference")
        met = [
            (1, "the lineage median difference's interval lies wholly above 0", c1(arm, learner)),
            (2, "better beats worse by the one-sided sign test at p < 0.05 after Holm over the 5 lineage pairs",
             c2(arm, learner)),
            (3, "the twin meets 1 in neither of its pairs (gbt, lasso) and 2 in neither, each read apart",
             not any(twin_1) and not any(twin_2)),
            ("3", f"beside, deciding nothing: the twin meets 1 and 2 in neither on {learner}, this arm's learner",
             not c1(TWIN, learner) and not c2(TWIN, learner)),
            (4, "the recall interval does not lie wholly below -0.02",
             recall["interval_hi"] is not None and recall["interval_hi"] >= -MARGIN),
            (5, "the study median has the lineage median's sign",
             lin is not None and study is not None and sign(lin) != 0 and sign(study) == sign(lin)),
        ]
        for n, text, m in met:
            out.append({"arm": arm, "learner": learner, "condition": int(n), "decides": int(isinstance(n, int)),
                        "reads": text, "met": int(bool(m))})
        out.append({"arm": arm, "learner": learner, "condition": 6, "decides": 1,
                    "reads": "the maintainer rules the change", "met": ""})
    return out


def _after(columns, at, add):
    i = columns.index(at) + 1
    return columns[:i] + [add] + columns[i:]


DRAW_COLUMNS = _after(C3.DRAW_COLUMNS, "constant", "unconverged")
DATASET_COLUMNS = _after(C3.DATASET_COLUMNS, "constant_draws", "unconverged_draws")
COUNT_COLUMNS = ["split", "pair", "learner", "measure", "orientation", "withdrawn", "datasets", "entered", "not_scored",
                 "constant", "unconverged", "better", "worse", "same", "median_difference", "species", "interval_lo",
                 "interval_hi", "sign_p", "holm_p"]
CONDITION_COLUMNS = ["arm", "learner", "condition", "decides", "reads", "met"]


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            print(__doc__.strip().splitlines()[1], file=sys.stderr)
            return 2
        out = Path(argv[1])
        rest, withdrawn = argv[2:], []
        while rest and rest[0] == "--withdrawn":
            if len(rest) < 2:
                raise SystemExit("--withdrawn takes an arm")
            withdrawn.append(rest[1])
            rest = rest[2:]
        if not rest:
            print(__doc__.strip().splitlines()[1], file=sys.stderr)
            return 2
        rows = [r for p in rest for r in read_tsv(p)]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise SystemExit("a dataset, split, arm and learner is listed twice")
        counted = count(rows, withdrawn)
        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "count.tsv", counted, COUNT_COLUMNS)
        write_tsv(out / "count-no-constant.tsv", count(rows, withdrawn, keep_constant=False), COUNT_COLUMNS)
        write_tsv(out / "count-no-unconverged.tsv", count(rows, withdrawn, keep_unconverged=False), COUNT_COLUMNS)
        write_tsv(out / "conditions.tsv", conditions(counted), CONDITION_COLUMNS)
        print(f"count: {len(rows)} dataset rows from {len(rest)} files")
        return 0
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[0], file=sys.stderr)
        return 2
    expect_path, run, release, out = (Path(a) for a in argv)
    expect = read_tsv(expect_path)
    draw_seeds = {int(e["seed"]) for e in expect if e["role"] == "draw"}
    twin_seeds = {int(e["seed"]) for e in expect if e["role"] == "twin"}
    served_by = {}
    if (release / "index.json").is_file():
        served_by = {x["id"]: x.get("best_learner") for x in C3.load(release / "index.json").get("datasets", [])}
    docs, names, checks = {}, {}, []
    for e in expect:
        path = run / e["file"]
        if path.is_file():
            docs[e["file"]] = C3.load(path)
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        if rel_path.is_file() and e["dataset"] not in names:
            rel = C3.load(rel_path)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in expect if e["arm"] == DEFAULT}
    for e in expect:
        if e["file"] not in docs:
            checks.append((e["file"], "exists", False, "a missing draw leaves its dataset unscored"))
            continue
        checks.append((e["file"], "exists", True, ""))
        checks.append((e["file"], "release report", e["dataset"] in names, ""))
        doc = docs[e["file"]]
        base = docs.get(base_of.get((e["dataset"], e["split"], e["seed"]))) if e["arm"] != DEFAULT else None
        if e["arm"] in (DEFAULT, TWIN):
            rel = None
            if e["arm"] == DEFAULT and int(e["seed"]) == C3.RELEASE_SEED and e["dataset"] in names:
                rel = C3.load(release / "datasets" / f"{e['dataset']}.json")
            found = C3.check(c3_row(e), doc, None, None, 0, rel, served_by.get(e["dataset"]), draw_seeds, twin_seeds)
        else:
            found = check_outside(e, doc, base)
        for name, passed, detail in found:
            checks.append((e["file"], name, passed, detail))
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "checks.tsv", [{"file": f, "check": c, "passed": p, "detail": d} for f, c, p, d in checks],
              ["file", "check", "passed", "detail"])
    failed = [c for c in checks if not c[2]]
    unavailable = sum(1 for d in docs.values() if d.get("unavailable_reason") is not None)
    print(f"checks: {len(checks) - len(failed)} of {len(checks)} pass; results {len(docs)} of {len(expect)}; "
          f"split unavailable in {unavailable}")
    if failed:
        return 1
    draws = score_rows(expect, docs, names)
    write_tsv(out / "draws.tsv", draws, DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", dataset_rows(draws), DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
