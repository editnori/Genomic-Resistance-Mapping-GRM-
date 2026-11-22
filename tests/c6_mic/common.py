"""Shared loaders and SYNTHETIC fixtures for the MIC fold-plan scripts (scripts/c6-mic, docs/mic-plan.md).

Every cohort, result and report here is made up in the test. No genome, MIC
or call comes from a real dataset; the ids are synthetic numbers in the
BV-BRC id shape.
"""

import importlib.util
import json
import math
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c6-mic"


def load(name):
    spec = importlib.util.spec_from_file_location(f"c6_mic_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def row(fold, log2, censor="Exact", resistant=False, released=None, ordered=None):
    """One SYNTHETIC tested genome with an MIC."""
    return {"fold": fold, "log2": log2, "censor": censor, "resistant": resistant,
            "released": released, "ordered": ordered, "median": None}


def b2(values):
    """The lower median, written here apart from score.py."""
    if not values:
        return None
    s = sorted(values)
    return s[math.ceil(len(s) / 2) - 1]


def result(rows, *, dataset="synthia_coli__drug_a", seed=42, arms=("released", "ordered"),
           r_above=(2, 2, 2, 2, 2), folds=5, unavailable=None, released_unavailable=None,
           ordered_unavailable=None, ordered_from_fold=None, fingerprint="f00d", plan=None):
    """A SYNTHETIC grm-mic-plan-result-v1 document, with the fold fields
    filled as grm learn fills them: a fold with no MIC test genome is not
    fitted, each fold's median is the lower median of the other folds'
    exact MICs, and an arm marked unavailable predicts no row.
    `ordered_from_fold` marks the ordered arm unfitted from that fold on."""
    rows = [dict(r) for r in rows]
    fold_docs = []
    median_missing = False
    for f in range(folds):
        mine = [r for r in rows if r["fold"] == f]
        train = [r for r in rows if r["fold"] != f]
        doc = {
            "train_mic": len(train),
            "train_exact": sum(1 for r in train if r["censor"] == "Exact"),
            "test_mic": len(mine),
            "test_exact": sum(1 for r in mine if r["censor"] == "Exact"),
            "test_resistant": sum(1 for r in mine if r["resistant"]),
            "breakpoint": None, "median": None, "released_trees": None, "ordered": None,
        }
        if mine:
            doc["median"] = b2([r["log2"] for r in train if r["censor"] == "Exact"])
            median_missing |= doc["median"] is None
            if r_above[f] is not None:
                doc["breakpoint"] = {"r_above": r_above[f], "agreement": 1.0, "n": len(train)}
            if "released" in arms and not released_unavailable:
                doc["released_trees"] = 7
            if "ordered" in arms and (ordered_from_fold is None or f < ordered_from_fold):
                doc["ordered"] = {"lo": -3, "hi": 5, "thresholds": 8, "fitted": 6, "constant": 1,
                                  "undecided": 1, "null_fits": 0, "pav_changed": 0}
        fold_docs.append(doc)
    if ordered_from_fold is not None and not ordered_unavailable:
        ordered_unavailable = f"fold {ordered_from_fold}'s training genomes give an empty fitting or holdout part"
    for r in rows:
        r["median"] = None if median_missing else fold_docs[r["fold"]]["median"]
        if released_unavailable or "released" not in arms:
            r["released"] = None
        if ordered_unavailable or "ordered" not in arms:
            r["ordered"] = None
    if unavailable:
        rows, fold_docs = [], []
    return {
        "format": "grm-mic-plan-result-v1",
        "evaluation_revision": 14,
        "input_fingerprint": fingerprint,
        "eval_config": {"folds": folds, "seed": 42},
        "id": dataset,
        "species_slug": dataset.split("__")[0],
        "antibiotic": dataset.split("__")[1].replace("_", " "),
        "plan": plan or {"kind": "lineage", "seed": seed},
        "arms": sorted(arms),
        "min_support": 3,
        "samples": len(rows) + 10,
        "resistant": sum(1 for r in rows if r["resistant"]),
        "n_groups": 12,
        "class_check": None,
        "n_mic": len(rows),
        "n_exact": sum(1 for r in rows if r["censor"] == "Exact"),
        "unavailable": unavailable,
        "released_unavailable": released_unavailable if "released" in arms else None,
        "ordered_unavailable": ordered_unavailable if "ordered" in arms else None,
        "median_unavailable": "fold 0's training genomes hold no exact MIC" if median_missing else None,
        "folds": fold_docs,
        "tested": rows,
        "seconds": 0.5,
    }


def expect_row(dataset="synthia_coli__drug_a", split="lineage", plan="lineage-s42", seed=42,
               arms="released,ordered"):
    sub = "lineage" if split == "lineage" else "study"
    return {"dataset": dataset, "species_slug": dataset.split("__")[0], "split": split, "plan": plan,
            "seed": str(seed), "arms": arms, "file": f"{sub}/{dataset}.{plan}.mic.json"}


def write_run(tmp, entries, fingerprint="f00d"):
    """A SYNTHETIC run under `tmp`: (expect.tsv path, results dir, release
    dir). `entries` holds (expect row, result document or None for a
    missing file)."""
    tmp = pathlib.Path(tmp)
    results, release = tmp / "results", tmp / "release"
    (release / "datasets").mkdir(parents=True, exist_ok=True)
    cols = ["dataset", "species_slug", "split", "plan", "seed", "arms", "file"]
    lines = ["\t".join(cols)]
    for exp, doc in entries:
        lines.append("\t".join(exp[c] for c in cols))
        report = release / "datasets" / f"{exp['dataset']}.json"
        report.write_text(json.dumps({"id": exp["dataset"], "input_fingerprint": fingerprint,
                                      "samples": [{"id": "SYNTHETIC"}]}))
        if doc is not None:
            path = results / exp["file"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(doc))
    expect = tmp / "expect.tsv"
    expect.write_text("\n".join(lines) + "\n")
    return expect, results, release


def read_tsv(path):
    lines = pathlib.Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return head, [dict(zip(head, line.split("\t"))) for line in lines[1:]]
