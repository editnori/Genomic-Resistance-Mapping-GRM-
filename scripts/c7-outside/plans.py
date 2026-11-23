"""The C7 run's check of its exports against the C3 run it reads (the lead's ruling on question 11, 23 Nov), made
after the export and before any fit.

  python3 -I plans.py <expect.tsv> <run dir> <release dir> <export dir with {seed}>

For every listed dataset and every seed the library arms run on the lineage split (expect.tsv, role draw):

- the export at that seed, <export dir> with {seed} filled in, holds a lineage plan of that seed;
- the plan's folds by class (each outer fold's training and test genomes and how many of each are resistant, and
  the genomes tested) equal the C3 default's result at the same dataset and seed (the k31 file expect.tsv lists,
  under c3/results/picks/), and a split unavailable in one is unavailable in the other;
- the plan's folds by genome equal the default's: each outer fold's genome count and the sha256 of its sorted genome
  ids, from the export's samples and lineage plan and from the grm learn report the C3 run wrote at that seed
  (c3/work/root/species/<slug>/results/ at seed 42, draws/seed-<seed>/ at any other: analysis::reports_dir). Two
  genomes of one class swapped between folds keep every count above and change two hashes. Only the pass or fail
  is printed, never an id or a hash;
- at seed 42, the export's input_fingerprint equals the release report's.

Reads nothing under a score/ folder. Prints each failed check (dataset, seed and check, no genome) and the counts,
and exits 1 when any check fails.
"""

import csv
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RELEASE_SEED = 42


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


O = _load("c7_outside_for_plans", HERE / "outside.py")


def export_shape(ex, seed):
    """The export's lineage plan at `seed` by class, as C3's fold_shape reads a result."""
    _, folds = O.export_plan(ex, "lineage", seed)
    if folds is None:
        return ("unavailable",)
    y = ex.y
    shape = tuple((int(f["train"].sum()), int((f["train"] & y).sum()), int(f["test"].sum()), int((f["test"] & y).sum()))
                  for f in folds)
    tested = sum(int(f["test"].sum()) for f in folds)
    return shape + ((tested, sum(int((f["test"] & y).sum()) for f in folds)),)


def result_shape(doc):
    """scripts/c3-inputs/score.py's fold_shape: the default's folds by class."""
    if doc.get("unavailable_reason") is not None:
        return ("unavailable",)
    return tuple((f["train"], f["train_resistant"], f["test"], f["test_resistant"]) for f in doc["folds"]) + (
        (doc["tested"], doc["tested_resistant"]),)


def report_path(run, slug, dataset, seed):
    """The grm learn report of the C3 default at `seed`: grm keeps the release seed's under results/ and any other
    seed's under draws/seed-<seed>/ (analysis::reports_dir)."""
    species = Path(run) / "c3" / "work" / "root" / "species" / slug
    return species / ("results" if seed == RELEASE_SEED else f"draws/seed-{seed}") / f"{dataset}.json"


def members(ids, folds):
    """Each outer fold's genomes as (fold, count, sha256 of its sorted ids one a line); no id is kept."""
    by = {}
    for i, f in zip(ids, folds):
        by.setdefault(int(f), []).append(i)
    return tuple((f, len(v), hashlib.sha256("\n".join(sorted(v)).encode()).hexdigest()) for f, v in sorted(by.items()))


def export_members(ex, seed):
    """The export's lineage plan at `seed` by genome."""
    plan, folds = O.export_plan(ex, "lineage", seed)
    if folds is None:
        return ("unavailable",)
    return members(ex.meta["samples"], plan["outer"])


def report_members(doc, seed):
    """A grm learn report's lineage split by genome: the fold that tested each sample (SplitResult.test_fold)."""
    if doc["eval_config"]["seed"] != seed:
        raise ValueError(f"the report is of seed {doc['eval_config']['seed']}, not {seed}")
    split = [x for x in doc["splits"] if x["scheme"] == "lineage"]
    if len(split) != 1:
        raise ValueError(f"the report holds {len(split)} lineage splits")
    if split[0].get("unavailable_reason") is not None:
        return ("unavailable",)
    ids = [x["id"] for x in doc["samples"]]
    if len(split[0]["test_fold"]) != len(ids):
        raise ValueError("the lineage split's test folds are not one a sample")
    return members(ids, split[0]["test_fold"])


def checks(expect_rows, run, release, exports):
    """[(dataset, seed, check, passed)] for every dataset and lineage seed of the library arms."""
    arm_seeds = {(r["dataset"], int(r["seed"])) for r in expect_rows
                 if r["arm"] in O.ARMS and r["split"] == "lineage" and r["role"] == "draw"}
    default = {(r["dataset"], int(r["seed"])): r["file"] for r in expect_rows
               if r["arm"] == "k31" and r["split"] == "lineage" and r["role"] == "draw"}
    slug = {r["dataset"]: r["species_slug"] for r in expect_rows}
    out = []
    for d, s in sorted(arm_seeds):
        path = Path(exports.format(seed=s)) / d
        try:
            ex = O.Export(path)
            got = export_shape(ex, s)
        except (O.Refuse, OSError, KeyError, ValueError) as e:
            out.append((d, s, f"the export at seed {s} is read: {e}", False))
            continue
        out.append((d, s, f"the export at seed {s} is read", True))
        name = default.get((d, s))
        if name is None or "score" in Path(name).parts:
            out.append((d, s, "the default's result at this seed is listed outside score/", False))
            continue
        try:
            with open(Path(run) / name, encoding="utf-8") as f:
                want = result_shape(json.load(f))
        except (OSError, KeyError, ValueError) as e:
            out.append((d, s, f"the default's result is read: {e}", False))
            continue
        out.append((d, s, "the plan's folds by class are the default's", got == want))
        report = report_path(run, slug[d], d, s)
        try:
            with open(report, encoding="utf-8") as f:
                same = export_members(ex, s) == report_members(json.load(f), s)
        except (O.Refuse, OSError, KeyError, TypeError, ValueError) as e:
            out.append((d, s, f"the default's report at seed {s} is read: {e}", False))
        else:
            out.append((d, s, "the plan's folds by genome are the default's", same))
        if s == RELEASE_SEED:
            try:
                with open(Path(release) / "datasets" / f"{d}.json", encoding="utf-8") as f:
                    released = json.load(f)["input_fingerprint"]
            except (OSError, KeyError, ValueError) as e:
                out.append((d, s, f"the release report is read: {e}", False))
                continue
            out.append((d, s, "the input fingerprint is the release report's", ex.meta["input_fingerprint"] == released))
    return out


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[3], file=sys.stderr)
        return 2
    with open(argv[0], encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    found = checks(rows, argv[1], argv[2], argv[3])
    failed = [c for c in found if not c[3]]
    for d, s, name, _ in failed:
        print(f"failed: {d} seed {s}: {name}")
    print(f"plans: {len(found) - len(failed)} of {len(found)} checks pass, "
          f"{len({(d, s) for d, s, _, _ in found})} dataset seeds")
    return 1 if failed or not found else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
