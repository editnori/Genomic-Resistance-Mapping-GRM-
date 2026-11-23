"""Shared loaders and SYNTHETIC fixtures for the settings-arm scripts (scripts/c5-settings, docs/settings-arms.md).

Every genome, label, fold, call and result here is made up in the test, and
the dataset names are made up. The driver runs against fake_grm.py, which
fits nothing.
"""

import importlib.util
import json
import os
import pathlib
import stat

from tests.c9_pick import common as c9

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c5-settings"
FAKE = pathlib.Path(__file__).resolve().parent / "fake_grm.py"
SPECIES = "synthia_coli"
DATASET = "synthia_coli__drug_a"
OTHER = "synthia_coli__drug_b"
LEARNERS = c9.LEARNERS
N_LINEAGES = 20
PER_LINEAGE = 3
FOLDS = 5


def load(name):
    spec = importlib.util.spec_from_file_location(f"c5_settings_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def fake():
    spec = importlib.util.spec_from_file_location("c5_fake_grm", FAKE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def genome(i):
    """A made-up id in the form of a BV-BRC genome id."""
    return f"{9000 + i}.{i + 1}"


def folds_doc(dataset=DATASET, scheme="lineage", seed=42):
    """A SYNTHETIC folds file as grm learn --exp-settings-folds writes it: 20 lineages of 3, half of them
    resistant, 5 outer folds of 4 lineages, and in each outer fold's training genomes 3 inner folds of whole
    lineages, each with both classes."""
    n = N_LINEAGES * PER_LINEAGE
    lineage = [i // PER_LINEAGE for i in range(n)]
    resistant = [lineage[i] % 2 == 0 for i in range(n)]
    outer = [lineage[i] % FOLDS for i in range(n)]
    inner = []
    for f in range(FOLDS):
        kept = sorted({lineage[i] for i in range(n) if outer[i] != f})
        at = {}
        for cls in (True, False):
            for k, line in enumerate(x for x in kept if (x % 2 == 0) == cls):
                at[line] = k % 3
        inner.append([at[lineage[i]] if outer[i] != f else None for i in range(n)])
    return {"format": "grm-c5-folds-1", "dataset": dataset, "scheme": scheme, "seed": seed, "folds": FOLDS,
            "inner_folds": 3, "unavailable_reason": None, "genomes": [genome(i) for i in range(n)],
            "resistant": resistant, "lineage": lineage, "outer": outer, "used": [True] * n, "inner": inner}


def world(doc, **extra):
    w = {"labels": {g: int(r) for g, r in zip(doc["genomes"], doc["resistant"])}}
    w.update(extra)
    return w


def default_pick(doc, w, seed=None):
    """The default's result for the folds file, as grm learn --exp-pick --exp-pick-inner 0 writes it,
    with each fold's calls the fake's release fit on that fold's test genomes."""
    f = fake()
    folds = []
    for k in range(doc["folds"]):
        test = [g for g, o in zip(doc["genomes"], doc["outer"]) if o == k]
        train = [i for i, o in enumerate(doc["outer"]) if o != k and doc["used"][i]]
        outer = [dict(learner=l, **f.calls(w, test, l, 0, 0.5, True)) for l in LEARNERS]
        folds.append({"fold": k, "train": len(train), "train_resistant": sum(doc["resistant"][i] for i in train),
                      "test": len(test), "test_resistant": sum(w["labels"][g] for g in test), "outer": outer,
                      "inner": {"status": "not run", "reason": None, "inner_folds": [], "learners": [], "tied": [],
                                "pick": None}})
    return {"format": "grm-c9-pick-1", "evaluation_revision": 14, "dataset": doc["dataset"], "scheme": doc["scheme"],
            "seed": doc["seed"] if seed is None else seed, "inner_folds": 0, "unavailable_reason": None,
            "folds": folds}


def setup(top, doc=None, **extra):
    """top/root (fake.json), top/grm (the fake's launcher), top/folds.json, top/pick.json; returns their paths."""
    top = pathlib.Path(top)
    doc = doc or folds_doc()
    w = world(doc, **extra)
    (top / "root").mkdir(parents=True, exist_ok=True)
    (top / "root" / "fake.json").write_text(json.dumps(w))
    grm = top / "grm"
    grm.write_text(f'#!/bin/sh\nexec python3 -I "{FAKE}" "$@"\n')
    grm.chmod(grm.stat().st_mode | stat.S_IXUSR)
    (top / "folds.json").write_text(json.dumps(doc))
    (top / "pick.json").write_text(json.dumps(default_pick(doc, w)))
    return grm, top / "root", top / "folds.json", top / "pick.json", top / "work", top / "out"


def calls_log(root):
    path = pathlib.Path(root) / "calls.jsonl"
    if not path.is_file():
        return []
    return [json.loads(line) for line in path.read_text().splitlines()]


def environ():
    return dict(os.environ)
