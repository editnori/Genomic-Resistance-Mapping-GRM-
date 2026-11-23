"""Shared loaders and SYNTHETIC fixtures for the input-arm scripts (scripts/c3-inputs, docs/input-arms.md).

Every result, inputs file, report and count here is made up in the test. No
genome or call comes from a real dataset, and the dataset names are made up.
The results are C9's writer's at no inner fold (tests/c9_pick/common.py).
"""

import importlib.util
import json
import pathlib

from tests.c9_pick import common as c9

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c3-inputs"
SPECIES = "synthia_coli"
DATASET = "synthia_coli__drug_a"
OTHER = "synthia_coli__drug_b"
LEARNERS = c9.LEARNERS
COLUMNS = ["dataset", "species_slug", "arm", "split", "scheme", "seed", "role", "file", "inputs"]
ALL_ARMS = ("program", "k31-program", "k31-counts", "k31-aa", "unitigs")
MATRIX_ARMS = {"k31-counts": "counts", "k31-aa": "aa", "unitigs": "unitigs"}
JOINED = ("k31-program", "k31-counts", "k31-aa")


def load(name):
    spec = importlib.util.spec_from_file_location(f"c3_inputs_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def pick(dataset=DATASET, scheme="lineage", seed=42, outer=None):
    """One SYNTHETIC result as grm learn --exp-pick --exp-pick-inner 0 writes it."""
    return c9.doc(dataset=dataset, scheme=scheme, seed=seed, inner=0, outer=outer)


def inputs(dataset=DATASET, arm="program", tag="drawn", program=4, k31=100, out_of_range=3, collapsed=1):
    """One SYNTHETIC inputs file as grm learn --exp-inputs writes it: counts only. In the
    joined arm two K31 rows equal a program column over the cohort and are left out. A matrix arm
    (k31-counts, k31-aa, unitigs) has its build's fields in place of the program's."""
    if arm in MATRIX_ARMS:
        return matrix_inputs(dataset=dataset, arm=arm, tag=tag, k31=k31)
    dropped = 2 if arm == "k31-program" else 0
    return {"format": "grm-c3-inputs-1", "evaluation_revision": 14, "dataset": dataset, "tag": tag, "arm": arm,
            "column": "symbol", "rows": "amr", "samples": 150, "min_support": 3, "tables_read": 150,
            "program_names_seen": program + out_of_range + collapsed, "program_names_out_of_range": out_of_range,
            "program_names_collapsed": collapsed, "program_rows": program, "k31_rows": k31,
            "k31_rows_as_program": dropped, "arm_rows": program + (k31 - dropped if arm == "k31-program" else 0)}


def matrix_inputs(dataset=DATASET, arm="k31-counts", tag="drawn", added=40, k31=100):
    """One SYNTHETIC inputs file of a matrix arm: counts only. In a joined arm two K31 rows
    equal an added row over the cohort and are left out."""
    joined = arm != "unitigs"
    d = {"format": "grm-c3-inputs-1", "evaluation_revision": 14, "dataset": dataset, "tag": tag, "arm": arm,
         "samples": 150, "min_support": 3, "k31_rows": k31, "arm_rows": added + (k31 - 2 if joined else 0),
         "added_build": MATRIX_ARMS[arm], "added_species_patterns": 500, "added_rows": added}
    if joined:
        d["k31_rows_as_added"] = 2
    if arm == "unitigs":
        d["unitig_pairs_held"], d["unitig_pairs_in_part"] = 9000, 120
    return d


def screen(doc, arm, first=3000, k31=2000):
    """A SYNTHETIC screen count of a joined arm's result `doc`, as grm learn writes it beside the
    result: per outer fold, each block's candidates and the three screened learners' counts."""
    folds = []
    for i, f in enumerate(doc["folds"]):
        n = f["train"]
        learners = []
        for learner, fitting in (("rf", n), ("gbt", n - round(0.15 * n)), ("lasso", n - round(0.2 * n))):
            screened = min(4000, first + k31)
            passed_first = min(first, screened * 5 // 8)
            learners.append({"learner": learner, "fitting": fitting, "screened": screened,
                             "passed_first": passed_first, "passed_k31": screened - passed_first})
        folds.append({"fold": i, "train": n, "candidates_first": first, "candidates_k31": k31, "learners": learners})
    return {"format": "grm-c3-screen-1", "evaluation_revision": 14, "dataset": doc["dataset"],
            "scheme": doc["scheme"], "seed": doc["seed"], "arm": arm,
            "first_block": {"k31-program": "program", "k31-counts": "counts", "k31-aa": "aa"}[arm],
            "screened_features": 4000, "folds": folds}


def write_run(top, rows, docs, files, releases):
    """Lay out a SYNTHETIC run: top/expect.tsv, top/picks/<arm>/<file>, the inputs files,
    top/release/datasets/<dataset>.json and top/release/index.json."""
    top = pathlib.Path(top)
    (top / "picks").mkdir(parents=True, exist_ok=True)
    (top / "release" / "datasets").mkdir(parents=True, exist_ok=True)
    with open(top / "expect.tsv", "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    for name, d in list(docs.items()) + list(files.items()):
        path = top / "picks" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            json.dump(d, f)
    for dataset, rep in releases.items():
        with open(top / "release" / "datasets" / f"{dataset}.json", "w") as f:
            json.dump(rep, f)
    index = []
    for dataset, rep in releases.items():
        rnd = [(l["metrics"]["balanced_accuracy"], i, l["learner"]) for s in rep["splits"] if s["scheme"] == "random"
               for i, l in enumerate(s["learners"])]
        index.append({"id": dataset, "best_learner": max(rnd)[2] if rnd else None})
    with open(top / "release" / "index.json", "w") as f:
        json.dump({"datasets": index}, f)
    return top / "expect.tsv", top / "picks", top / "release", top / "score"


def full_run(datasets=(DATASET, OTHER), a=(DATASET,), study=(DATASET,), seeds=(42, 43, 44), twin=(45, 46, 47),
             arms=("program", "k31-program"), change=None, change_screen=None):
    """SYNTHETIC expect rows, results, inputs files and release reports of one species' whole run.

    change(row, doc, inputs) may alter a result or its inputs file before it is written, and
    change_screen(row, screen) a joined arm's screen count, written beside its result.
    """
    expect = load("expect")
    rows = expect.rows(SPECIES, list(datasets), list(a), list(study), list(seeds), list(twin), list(arms))
    docs, files, releases = {}, {}, {}
    for r in rows:
        dataset, _, arm, split, scheme, seed, role, name, inputs_name = r
        d = pick(dataset=dataset, scheme=scheme, seed=seed)
        i = None
        if inputs_name != "-":
            i = inputs(dataset=dataset, arm=arm, tag="drawn" if split == "lineage" else scheme)
        if change:
            change(r, d, i)
        docs[name] = d
        if i is not None:
            files[inputs_name] = i
        if arm in JOINED:
            s = screen(d, arm)
            if change_screen:
                change_screen(r, s)
            files[name[: -len(".pick.json")] + ".screen.json"] = s
        if arm == "k31" and split == "lineage" and seed == 42:
            releases[dataset] = c9.release(d, dataset)
    return rows, docs, files, releases
