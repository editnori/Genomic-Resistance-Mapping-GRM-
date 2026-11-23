"""Shared loaders and SYNTHETIC fixtures for the screen-arm scripts (scripts/c4-screen, docs/screen-arms.md).

Every result, inputs file, report and count here is made up in the test. No
genome or call comes from a real dataset, and the dataset names are made up.
The results are C9's writer's at no inner fold (tests/c9_pick/common.py),
laid out as C3's tests lay them out (tests/c3_inputs/common.py).
"""

import importlib.util
import pathlib

from tests.c3_inputs import common as c3
from tests.c9_pick import common as c9

ROOT = c3.ROOT
KIT = ROOT / "scripts" / "c4-screen"
SPECIES, DATASET, OTHER = c3.SPECIES, c3.DATASET, c3.OTHER
LEARNERS = c9.LEARNERS
ARMS = ("cmh-4000", "cmh-16000", "floor-3")
pick, write_run = c3.pick, c3.write_run
# Each SYNTHETIC fold's candidates past the floor, the dataset's rows, and
# the genomes boosting and the lasso fit after their holdouts.
CANDIDATES, DATASET_ROWS = 500, 600
FITTING = {"rf": lambda t: t, "gbt": lambda t: t * 85 // 100, "lasso": lambda t: t * 80 // 100}


def load(name):
    spec = importlib.util.spec_from_file_location(f"c4_screen_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def floor_of(arm, samples):
    return 3 if arm == "floor-3" else max(3, samples // 200)


def inputs(doc, arm="cmh-4000", tag="drawn", candidates=CANDIDATES, random_folds=None):
    """One SYNTHETIC inputs file as grm learn --exp-screen or --exp-floor writes it beside the result `doc`.

    A drawn run's file holds its random split's folds (random_folds, the
    lineage result's own by default) and the lineage result's, and the final
    models; a fold file's run, its scheme's folds alone. Counts only.
    """
    rule = "chi-square" if arm == "floor-3" else "cmh"
    budget = 16000 if arm == "cmh-16000" else 4000
    samples = doc["samples"]
    splits = [(doc["scheme"], [f["train"] for f in doc["folds"]])]
    if tag == "drawn":
        splits.insert(0, ("random", random_folds or [f["train"] for f in doc["folds"]]))
    folds, screens, times = [], [], []

    def screen(split, fold, learner, fitting, rows):
        adjusted = rule == "cmh"
        kept = min(budget, rows)
        screens.append({"split": split, "fold": fold, "learner": learner, "rule": rule, "fitting": fitting,
                        "candidates": rows, "kept": kept, "kept_at_zero": 7 if adjusted else None,
                        "strata": 12 if adjusted else None, "informative_strata": 9 if adjusted else None})

    for split, trains in splits:
        for fi, train in enumerate(trains):
            folds.append({"split": split, "fold": fi, "train": train, "candidates": candidates})
            for learner in LEARNERS:
                if learner in FITTING:
                    screen(split, fi, learner, FITTING[learner](train), candidates)
                times.append({"split": split, "fold": fi, "learner": learner, "seconds": 0.5, "cpu_seconds": 1.5,
                              "peak_kb": 100_000})
    if tag == "drawn":
        for learner in ("lasso", "rf", "gbt"):
            screen("final", None, learner, FITTING[learner](samples), DATASET_ROWS)
    return {"format": "grm-c4-screen-1", "evaluation_revision": 14, "dataset": doc["dataset"], "tag": tag,
            "seed": doc["seed"], "arm": "floor" if arm == "floor-3" else "cmh", "rule": rule, "budget": budget,
            "min_support": floor_of(arm, samples), "samples": samples, "dataset_rows": DATASET_ROWS, "folds": folds,
            "screens": screens, "fit_times": times}


def full_run(datasets=(DATASET, OTHER), study=(DATASET,), seeds=(42, 43, 44), twin=(45, 46, 47), arms=ARMS,
             samples=150, change=None):
    """SYNTHETIC expect rows, results, inputs files and release reports of one species' whole run.

    samples is every dataset's genome count (a floor arm moves at 800 or
    more). change(row, doc, inputs) may alter a result or its inputs file
    before it is written.
    """
    expect = load("expect")
    rows = expect.rows(SPECIES, list(datasets), list(study), list(seeds), list(twin), list(arms))
    docs, files, releases = {}, {}, {}
    for r in rows:
        dataset, _, arm, split, scheme, seed, role, name, inputs_name = r
        d = pick(dataset=dataset, scheme=scheme, seed=seed)
        d["samples"] = samples
        d["min_support"] = floor_of(arm, samples)
        i = None
        if inputs_name != "-":
            i = inputs(d, arm=arm, tag="drawn" if split == "lineage" else scheme)
        if change:
            change(r, d, i)
        docs[name] = d
        if i is not None:
            files[inputs_name] = i
        if arm == "k31" and split == "lineage" and seed == 42:
            releases[dataset] = c9.release(d, dataset)
    return rows, docs, files, releases
