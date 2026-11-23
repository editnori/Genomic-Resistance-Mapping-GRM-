"""The scorer of plan block C5 (docs/settings-arms.md, measures and files).

  python3 -I score.py <expect.tsv> <picks dir> <release dir> <out dir>
  python3 -I score.py --count <out dir> <datasets.tsv> [<datasets.tsv> ...]

It is block C3's scorer (scripts/c3-inputs/score.py), loaded as a module of
its own and set to C5's family: the default (`default`, the release's
settings), the search arm, the cutoff arm and the twin (the default at seeds
of its own). C3's code checks the default's and the twin's results (each a
grm learn --exp-pick result at no inner fold, held at seed 42 to the release
report), scores every draw, takes each dataset's mean over its draws and
counts; C5's own code checks the driver's choices files and reads the two
arms' calls from them.

The first form checks every result expect.py lists. The search and the
cutoff arm of a draw are read from its choices files, one per outer fold
(drive.py), as many as the default's result has folds. Per choices file: the format, revision, dataset, split, seed and
fold of its row; the sampler seed, gamma, grid and 30 random trials; the
fold's training and test genomes by class equal to the default's fold; the
five learners in order; the default's calls equal to the default's own
result on that fold (the default is never refitted with a searched
setting); each arm's calls of the fold's test genomes; the pick a trial
from 0 to 30, trial 0 unless the best inner mean beats the release's by
more than gamma, the cutoff a grid point or 0.5; scm's cutoff arm the
default's calls, since scm scores 0 or 1. A draw with a missing choices
file is not scored (its rows say `missing draw`), and its dataset's mean is
not taken: a missing draw is never read as zero. Any failed check of a file
that exists stops the scorer before it scores. Writes checks.tsv, draws.tsv
and datasets.tsv, and prints the check counts only.

The second form writes count.tsv, count-no-constant.tsv and conditions.tsv
over several species' datasets.tsv with C3's count: per split, arm and
learner, each arm against the default on the same datasets; Holm over the
15 lineage rows of the primary measure (search, cutoff and the twin by five
learners); conditions read two ways, ruling 15 deciding, with condition 3
read from the twin as corrected. It refuses unless the search, the cutoff
and the twin all have lineage rows.
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT = "default"
TWIN = "twin"
ARMS = ["search", "cutoff"]
FAMILY = ARMS + [TWIN]
CHOICE_FORMAT = "grm-c5-choices-1"
SAMPLER_SEED = 20261030
GAMMA = 0.04
GRID = [round(0.05 * k, 2) for k in range(1, 20)]
TRIALS = 30
FOLDS = 5
HALF = 0.5


def _c3():
    """C3's scorer as a module of its own, set to C5's family; C3's own runs load theirs apart."""
    spec = importlib.util.spec_from_file_location("c3_score_for_c5", HERE.parent / "c3-inputs" / "score.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.DEFAULT, module.TWIN, module.FAMILY = DEFAULT, TWIN, list(FAMILY)
    return module


C3 = _c3()
C9 = C3.C9
LEARNERS = C3.LEARNERS
CALL_KEYS = ("n", "tp", "fp", "tn", "fn", "loss_sum")
r12, read_tsv, write_tsv = C3.r12, C3.read_tsv, C3.write_tsv


def fold_file(e, k):
    return e["file"].replace("<k>", str(k))


def same_calls(a, b):
    return all(a[k] == b[k] for k in CALL_KEYS[:-1]) and abs(a["loss_sum"] - b["loss_sum"]) <= 1e-9


def check_choices(exp, files, base):
    """[(check name, passed, detail)] for a draw's choices files (present ones) against the default's result."""
    out = []

    def ok(name, cond, detail=""):
        out.append((name, bool(cond), detail))
        return bool(cond)

    ok("the default's result at the same dataset, split and seed is there and available",
       base is not None and base.get("unavailable_reason") is None and len(base.get("folds", [])) >= 2)
    for k, c in enumerate(files):
        if c is None:
            continue
        ok(f"fold {k}: format, revision, dataset, split, seed and fold",
           (c.get("format"), c.get("evaluation_revision"), c.get("dataset"), c.get("scheme"), c.get("seed"), c.get("fold"))
           == (CHOICE_FORMAT, C3.REVISION, exp["dataset"], exp["scheme"], int(exp["seed"]), k))
        ok(f"fold {k}: sampler seed, gamma, grid and 30 random trials",
           (c.get("sampler_seed"), c.get("gamma"), c.get("grid"), c.get("random_trials"))
           == (SAMPLER_SEED, GAMMA, GRID, TRIALS))
        if not ok(f"fold {k}: the five learners in order", [r["learner"] for r in c.get("learners", [])] == LEARNERS):
            continue
        if base is None or base.get("unavailable_reason") is not None or k >= len(base.get("folds", [])):
            continue
        b = base["folds"][k]
        n = c["counts"]
        ok(f"fold {k}: the training and test genomes by class are the default's fold",
           (n["train"], n["train_resistant"], n["test"], n["test_resistant"])
           == (b["train"], b["train_resistant"], b["test"], b["test_resistant"]))
        for li, r in enumerate(c["learners"]):
            l = r["learner"]
            want = {key: b["outer"][li][key] for key in CALL_KEYS}
            ok(f"fold {k} {l}: the default's calls are the default's own fit", same_calls(r["outer"]["default"], want))
            ok(f"fold {k} {l}: each arm's calls are the fold's test genomes'",
               all(r["outer"][a] is not None and r["outer"][a]["n"] == b["test"]
                   and r["outer"][a]["tp"] + r["outer"][a]["fn"] == b["test_resistant"] for a in ARMS))
            rel, best = r["release_inner_mean"], r["best_inner_mean"]
            lead = None if rel is None or best is None else r12(best - rel)
            ok(f"fold {k} {l}: the pick follows the default-anchored rule",
               0 <= r["picked_trial"] <= TRIALS
               and (r["picked_trial"] == 0) == (lead is None or lead <= GAMMA)
               and (r["picked_trial"] == 0 or r["picked_trial"] == r["plain_best_trial"]))
            ok(f"fold {k} {l}: the cutoff is a grid point or 0.5", r["cutoff"] in GRID or r["cutoff"] == HALF)
            if l == "scm":
                ok(f"fold {k} scm: the cutoff arm's calls are the default's", same_calls(r["outer"]["cutoff"], want))
    return out


def arm_doc(files, arm):
    """The arm's calls in the form of a grm learn --exp-pick result: per fold, the five learners' outer calls."""
    return {"unavailable_reason": None,
            "folds": [{"fold": k, "outer": [dict({key: r["outer"][arm][key] for key in CALL_KEYS}, learner=r["learner"])
                                           for r in c["learners"]]} for k, c in enumerate(files)]}


def missing_rows(e, names):
    sp, ab = names[e["dataset"]]
    return [{"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "split": e["split"],
             "seed": int(e["seed"]), "role": e["role"], "learner": l, "scored": 0, "constant": 0,
             "why": "missing draw"} for l in LEARNERS]


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
        counted = C3.count(rows)
        out.mkdir(parents=True, exist_ok=True)
        write_tsv(out / "count.tsv", counted, C3.COUNT_COLUMNS)
        write_tsv(out / "count-no-constant.tsv", C3.count(rows, keep_constant=False), C3.COUNT_COLUMNS)
        write_tsv(out / "conditions.tsv", C3.conditions(counted), C3.CONDITION_COLUMNS)
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
    docs, names, checks, missing = {}, {}, [], []
    for e in expect:
        rel_path = release / "datasets" / f"{e['dataset']}.json"
        if rel_path.is_file() and e["dataset"] not in names:
            rel = load(rel_path)
            names[e["dataset"]] = (rel["species"], rel["antibiotic"])
        if e["arm"] in (DEFAULT, TWIN) and (picks / e["file"]).is_file():
            docs[e["file"]] = load(picks / e["file"])
    base_of = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in expect if e["arm"] == DEFAULT}
    arm_docs, checked = {}, set()
    for e in expect:
        checks.append((e["file"], "release report", e["dataset"] in names, ""))
        if e["arm"] in (DEFAULT, TWIN):
            if e["file"] not in docs:
                missing.append(e)
                continue
            rel = None
            if e["arm"] == DEFAULT and int(e["seed"]) == C3.RELEASE_SEED and e["dataset"] in names:
                rel = load(release / "datasets" / f"{e['dataset']}.json")
            for name, passed, detail in C3.check(dict(e, inputs="-"), docs[e["file"]], None, None, 0, rel,
                                                 served_by.get(e["dataset"]), draw_seeds, twin_seeds):
                checks.append((e["file"], name, passed, detail))
            continue
        if e["role"] == "twin" or int(e["seed"]) in twin_seeds:
            checks.append((e["file"], "an arm runs at a draw's seed or the study's, never the twin's", False, ""))
            continue
        base = docs.get(base_of.get((e["dataset"], e["split"], e["seed"])))
        # One choices file per outer fold of the default's result.
        n = len(base["folds"]) if base is not None and base.get("folds") else FOLDS
        files = [load(picks / fold_file(e, k)) if (picks / fold_file(e, k)).is_file() else None for k in range(n)]
        # Both arms read the same files: they are checked once, on the first arm's row.
        if (e["file"], e["split"]) not in checked:
            checked.add((e["file"], e["split"]))
            for name, passed, detail in check_choices(e, files, base):
                checks.append((e["file"], name, passed, detail))
        if any(c is None for c in files):
            missing.append(e)
            continue
        arm_docs[e["file"] + "|" + e["arm"]] = arm_doc(files, e["arm"])
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "checks.tsv", [{"file": f, "check": c, "passed": p, "detail": d} for f, c, p, d in checks],
              ["file", "check", "passed", "detail"])
    failed = [c for c in checks if not c[2]]
    print(f"checks: {len(checks) - len(failed)} of {len(checks)} pass; results {len(expect) - len(missing)} of "
          f"{len(expect)}, missing draws {len(missing)}")
    if failed:
        return 1
    have = [e for e in expect if e not in missing]
    keyed = [dict(e, file=e["file"] + "|" + e["arm"]) if e["arm"] in ARMS else e for e in have]
    draws = C3.score_rows(keyed, {**docs, **arm_docs}, names)
    for e in missing:
        draws += missing_rows(e, names)
    for r in draws:
        r.setdefault("constant", 0)
    write_tsv(out / "draws.tsv", draws, C3.DRAW_COLUMNS)
    write_tsv(out / "datasets.tsv", C3.dataset_rows(draws), C3.DATASET_COLUMNS)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
