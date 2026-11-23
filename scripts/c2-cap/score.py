"""The scorer of plan block C2 (docs/cap-arms.md, Draws, the twin and the count).

  python3 -I score.py <expect.tsv> <fits dir> <picks dir> <out dir>
  python3 -I score.py --count <out dir> [--withdrawn <arm> ...] <datasets.tsv> [<datasets.tsv> ...]

Runs on the box. The first form checks every row expect.py lists, then
scores it, and writes everything under <out dir> (score/ in run_species.sh,
never logs/). Checks, per row:

- the dataset is in the first run (scripts/c2-cap/lists.tsv: CAP, not M);
- a lineage row has its plan (the default's plans/lineage-<slot>.json, a
  wider arm's plans/lineage-<slot>.cohort.json) and its prediction rows,
  five learners, one row per scored genome and learner in each outer fold;
  a missing draw stops the scorer, never counted as zero;
- a twin row is the default's own fit (cap-2000's folder) at a slot no
  draw has; twin rows from any other folder stop the scorer;
- a wider arm's scored genomes in each outer fold are the default's test
  genomes of that fold at the same slot;
- a default lineage plan with an unavailable_reason (a draw grm arm refused
  past its replacements) makes every lineage row of the dataset, on every
  arm and the twin, unavailable: not produced, as H-1's ruling has it
  (docs/input-comparison-runner.md, A refused repeat), so its prediction
  rows are not read and it is scored 0 with why "lineage split
  unavailable", out of each arm's count; the study rows are scored as
  usual;
- a study row has its grm learn --exp-pick file (format grm-c9-pick-1, the
  dataset, scheme study-whole, seed 42, no inner fold), and a wider arm's
  test folds hold the default's tested and resistant counts, fold by fold,
  with at least the default's training genomes.

Scores, per dataset, arm, split, slot and learner: the outer calls summed
over the folds (n, tp, fp, tn, fn and the log loss of the score clipped to
[1e-4, 1 - 1e-4]), and from them the measures of scripts/c9-pick/score.py
with their Wilson intervals. Writes checks.tsv, draws.tsv and datasets.tsv
(per dataset, split, arm and learner the mean over the draws,
scripts/c3-inputs/score.py's dataset_rows) and prints the check counts
only. An arm fitted by another (expect.tsv's fit) takes that fit's rows,
named in draws.tsv's fit. A none-added row (no added genome labelled or,
on the study split, none that trains) and a same-cohort row (an arm whose
cohort is another arm's, Enterococcus faecium's cap-4000) are not scored,
so the dataset is out of that arm's count.

The second form takes the count and the conditions from an instance of
scripts/c3-inputs/score.py of its own, with the family cap-4000, cap-all
and the twin, the default cap-2000: Holm over the 15 lineage rows of the
primary measure, 3 arms by 5 learners. --withdrawn names an arm withdrawn
before any fit, as C3's scorer takes it: the arm may hold no row, its five
lineage pairs stay in the family at p = 1, count.tsv's withdrawn column is
1 on its rows, and conditions.tsv names it withdrawn and reads none of its
conditions (met and learners empty). cap-4000 is withdrawn for the run
(docs/cap-arms.md). Rows of a none-added arm never reach datasets.tsv, so
such a dataset is out of the arm's denominator. Writes count.tsv and
conditions.tsv and prints the row count only.
"""

import csv
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT, TWIN = "cap-2000", "twin"
FAMILY = ["cap-4000", "cap-all", TWIN]


def _c3():
    spec = importlib.util.spec_from_file_location("c3_score_for_c2", HERE.parent / "c3-inputs" / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.DEFAULT, module.TWIN, module.FAMILY = DEFAULT, TWIN, list(FAMILY)
    return module


C3 = _c3()
C9 = C3.C9
LEARNERS = C9.LEARNERS
MEASURES = C9.MEASURES
LOSS_CLIP = C3.LOSS_CLIP
read_tsv, write_tsv = C9.read_tsv, C9.write_tsv
DRAW_COLUMNS = ["dataset", "species", "antibiotic", "arm", "fit", "split", "slot", "role", "learner", "scored", "constant",
                "why"] + C3.DRAW_COLUMNS[C3.DRAW_COLUMNS.index("n"):]


def load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def first_run():
    """{dataset: (species, antibiotic)} of the first run."""
    return {r["dataset"]: (r["species"], r["antibiotic"]) for r in read_tsv(HERE / "lists.tsv")
            if r["first_run"] == "1" and r["in_M"] == "0"}


def loss(y, score):
    p = min(max(float(score), LOSS_CLIP), 1.0 - LOSS_CLIP)
    return -math.log(p) if y else -math.log(1.0 - p)


def predictions(path):
    """{(learner, slot): {fold: [(sample_id, y, called, score)]}} of the lineage split's rows (a row's `repeat` is
    its slot)."""
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            r = json.loads(line)
            if r["split"] != "lineage":
                continue
            out.setdefault((r["learner"], r["repeat"]), {}).setdefault(r["fold"], []).append(
                (r["sample_id"], int(r["y"]), bool(r["called"]), r["score"]))
    return out


def summed(folds):
    c = dict(C9.ZERO)
    for rows in folds.values():
        for _, y, called, score in rows:
            c["n"] += 1
            if y and called:
                c["tp"] += 1
            elif y:
                c["fn"] += 1
            elif called:
                c["fp"] += 1
            else:
                c["tn"] += 1
            c["loss_sum"] += loss(y, score)
    return c


class Lineage:
    """The lineage split's fits, read once each."""

    def __init__(self, fits):
        self.fits = Path(fits)
        self.cache = {}

    def fit(self, fit, dataset):
        key = (fit, dataset)
        if key not in self.cache:
            d = self.fits / fit / dataset
            p = d / "predictions.jsonl"
            self.cache[key] = predictions(p) if p.is_file() else None
        return self.cache[key]

    def plan(self, fit, dataset, slot):
        d = self.fits / fit / dataset / "plans"
        p = d / (f"lineage-{slot}.json" if fit == DEFAULT else f"lineage-{slot}.cohort.json")
        if not p.is_file():
            return None
        doc = load(p)
        return doc["plan"] if fit == DEFAULT else doc


def lineage_checks(e, lin, draw_slots, checks):
    """The row's {learner: {fold: rows}} or None, its checks appended."""
    where = f"{e['fit']}/{e['dataset']} lineage slot {e['slot']}"
    if e["role"] == "twin":
        ok = e["fit"] == DEFAULT and int(e["slot"]) not in draw_slots
        checks.append((where, "the twin is the default's fit at a slot no draw has", ok, ""))
        if not ok:
            return None
    plan = lin.plan(e["fit"], e["dataset"], int(e["slot"]))
    rows = lin.fit(e["fit"], e["dataset"])
    have = plan is not None and rows is not None
    checks.append((where, "the plan and the prediction rows exist", have, ""))
    if not have:
        return None
    out = {}
    for learner in LEARNERS:
        # Prediction rows carry the slot as `repeat`; the plan's `repeat` is the keyed repeat, another number
        # when grm arm replaced a refused draw (src/ablate/arm.rs).
        folds = rows.get((learner, int(e["slot"])))
        checks.append((where, f"{learner} has rows", folds is not None, ""))
        if folds is None:
            return None
        out[learner] = folds
    k = plan["k"]
    if e["fit"] == DEFAULT:
        want = {}
        for f in plan["outer"]:
            want[f] = want.get(f, 0) + 1
    else:
        want = {}
        for f in plan["outer"]:
            if f is not None:
                want[f] = want.get(f, 0) + 1
    for learner, folds in out.items():
        counts = {f: len(r) for f, r in folds.items()}
        checks.append((where, f"{learner}: one row per scored genome in each of {k} folds", counts == want, ""))
    return out


def same_scored(arm_folds, base_folds):
    a = {f: {r[0] for r in rows} for f, rows in arm_folds.items()}
    b = {f: {r[0] for r in rows} for f, rows in base_folds.items()}
    return a == b


def study_checks(e, picks, checks):
    where = f"{e['fit']}/{e['dataset']} study"
    p = Path(picks) / e["fit"] / f"{e['dataset']}.study-whole-s{C9.RELEASE_SEED}.pick.json"
    if not p.is_file():
        checks.append((where, "the pick file exists", False, ""))
        return None
    doc = load(p)
    ok = (doc.get("format") == C9.FORMAT and doc.get("evaluation_revision") == C9.REVISION
          and (doc.get("dataset"), doc.get("scheme"), doc.get("seed"), doc.get("inner_folds"))
          == (e["dataset"], "study-whole", C9.RELEASE_SEED, 0)
          and [x["learner"] for x in doc.get("learners", [])] == LEARNERS)
    checks.append((where, "the pick file is the dataset's study-whole fit at seed 42 with no inner fold", ok, ""))
    return doc if ok else None


def test_shape(doc):
    if doc.get("unavailable_reason") is not None:
        return ("unavailable",)
    return tuple((f["test"], f["test_resistant"]) for f in doc["folds"])


def score(expect, fits, picks):
    """(checks, draws.tsv rows)."""
    names = first_run()
    checks, draws = [], []
    draw_slots = {int(e["slot"]) for e in expect if e["role"] == "draw"}
    lin = Lineage(fits)
    base_lineage, base_study = {}, {}
    ok_rows = []
    refused = {}

    def lineage_refused(d):
        """True when a default lineage plan of d holds an unavailable_reason."""
        if d not in refused:
            slots = sorted({int(x["slot"]) for x in expect
                            if x["dataset"] == d and x["split"] == "lineage" and x["fit"] == DEFAULT})
            plans = [lin.plan(DEFAULT, d, s) for s in slots]
            refused[d] = any(p is not None and p.get("unavailable_reason") is not None for p in plans)
        return refused[d]

    for e in expect:
        if e["dataset"] not in names:
            checks.append((e["dataset"], "the dataset is in the first run (CAP, not M)", False, ""))
            continue
        if e["role"] in ("none-added", "same-cohort"):
            continue
        if e["split"] == "lineage" and lineage_refused(e["dataset"]):
            checks.append((f"{e['fit']}/{e['dataset']} lineage slot {e['slot']}",
                           "the lineage split is unavailable (a default plan refused past grm arm's replacements); "
                           "not produced", True, ""))
            ok_rows.append((e, None))
            continue
        if e["split"] == "lineage":
            folds = lineage_checks(e, lin, draw_slots, checks)
            if folds is None:
                continue
            if e["arm"] == DEFAULT:
                base_lineage[(e["dataset"], e["slot"])] = folds
            ok_rows.append((e, folds))
        else:
            doc = study_checks(e, picks, checks)
            if doc is None:
                continue
            if e["arm"] == DEFAULT:
                base_study[e["dataset"]] = doc
            ok_rows.append((e, doc))
    for e, got in ok_rows:
        if e["arm"] in (DEFAULT, TWIN) or got is None:
            continue
        where = f"{e['fit']}/{e['dataset']} {e['split']} slot {e['slot']}"
        if e["split"] == "lineage":
            base = base_lineage.get((e["dataset"], e["slot"]))
            ok = base is not None and all(same_scored(got[l], base[l]) for l in LEARNERS)
            checks.append((where, "the scored genomes of each outer fold are the default's", ok, ""))
        else:
            base = base_study.get(e["dataset"])
            ok = (base is not None and test_shape(got) == test_shape(base)
                  and all(a["train"] >= b["train"] for a, b in zip(got["folds"], base["folds"])))
            checks.append((where, "the test folds are the default's, with at least its training genomes", ok, ""))
    if any(not c[2] for c in checks):
        return checks, []
    for e, got in ok_rows:
        sp, ab = names[e["dataset"]]
        base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "fit": e["fit"],
                "split": e["split"], "slot": e["slot"], "role": e["role"], "constant": 0}
        for learner in LEARNERS:
            if e["split"] == "lineage" and got is None:
                draws.append(dict(base, learner=learner, scored=0, why="lineage split unavailable"))
            elif e["split"] == "lineage":
                draws.append(dict(base, learner=learner, scored=1, why="", **C9.measures(summed(got[learner]))))
            elif got.get("unavailable_reason") is not None:
                draws.append(dict(base, learner=learner, scored=0, why="split unavailable"))
            else:
                draws.append(dict(base, learner=learner, scored=1, why="", **C9.measures(C9.fixed(got, learner))))
    return checks, draws


CONDITION_COLUMNS = C3.CONDITION_COLUMNS[:1] + ["withdrawn"] + C3.CONDITION_COLUMNS[1:]


def conditions(counted, withdrawn):
    """C3's conditions, a withdrawn arm's rows named withdrawn and not read."""
    out = []
    for r in C3.conditions(counted):
        gone = r["arm"] in withdrawn
        out.append(dict(r, withdrawn=int(gone), **({"met": "", "learners": ""} if gone else {})))
    return out


def main(argv):
    if argv and argv[0] == "--count":
        if len(argv) < 3:
            raise SystemExit(__doc__)
        out = Path(argv[1])
        rest, withdrawn = argv[2:], []
        while rest and rest[0] == "--withdrawn":
            if len(rest) < 2:
                raise SystemExit("--withdrawn takes an arm")
            withdrawn.append(rest[1])
            rest = rest[2:]
        if not rest:
            raise SystemExit(__doc__)
        rows = [r for p in rest for r in read_tsv(p)]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        if len(set(keys)) != len(keys):
            raise SystemExit("a dataset, split, arm and learner is listed twice")
        counted = C3.count(rows, withdrawn)
        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "count.tsv", counted, C3.COUNT_COLUMNS)
        # C3's conditions write its family's size, 15 here: 3 arms by 5 learners.
        write_tsv(out / "conditions.tsv", conditions(counted, withdrawn), CONDITION_COLUMNS)
        print(f"count: {len(rows)} dataset rows from {len(rest)} files")
        return 0
    if len(argv) != 4:
        raise SystemExit(__doc__)
    expect_path, fits, picks, out = (Path(a) for a in argv)
    expect = read_tsv(expect_path)
    checks, draws = score(expect, fits, picks)
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "checks.tsv", [{"row": w, "check": c, "passed": p, "detail": d} for w, c, p, d in checks],
              ["row", "check", "passed", "detail"])
    failed = sum(1 for c in checks if not c[2])
    print(f"checks: {len(checks) - failed} of {len(checks)} pass; rows {len(expect)}")
    if failed:
        return 1
    write_tsv(out / "draws.tsv", draws, DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", C3.dataset_rows(draws), C3.DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
