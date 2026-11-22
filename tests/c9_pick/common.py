"""Shared loaders and SYNTHETIC fixtures for the served-learner scripts (scripts/c9-pick, docs/served-pick.md).

Every result, report and count here is made up in the test. No genome or
call comes from a real dataset, and the dataset names are made up.
"""

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c9-pick"
SPECIES = "synthia_coli"
DATASET = "synthia_coli__drug_a"
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
# Each outer fold tests 20 resistant and 30 susceptible genomes and trains
# on the other 100 (40 and 60); each inner fold of those 100 tests these.
OUTER_POS, OUTER_NEG, FOLDS = 20, 30, 3
INNER_POS, INNER_NEG = (14, 13, 13), (20, 20, 20)


def load(name):
    spec = importlib.util.spec_from_file_location(f"c9_pick_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def calls(tp, tn, pos, neg, loss=0.3):
    return {"n": pos + neg, "tp": tp, "fp": neg - tn, "tn": tn, "fn": pos - tp, "loss_sum": loss * (pos + neg)}


def balanced(c):
    return (c["tp"] / (c["tp"] + c["fn"]) + c["tn"] / (c["tn"] + c["fp"])) / 2


def outer_default(f, li):
    """SYNTHETIC outer (tp, tn) of learner li on fold f."""
    return 10 + 2 * li, 20 + li + f


def inner_default(f, li):
    """SYNTHETIC inner (tp, tn) per inner fold: rf (li 2) best in every fold."""
    return [(8 + (3 if li == 2 else li % 2), 14 + (2 if li == 2 else 0)) for _ in range(3)]


def doc(dataset=DATASET, scheme="lineage", seed=42, inner=3, tie="later", outer=None, inner_calls=None,
        default=("rf", 0.8, ["rf"]), no_pick=()):
    """One SYNTHETIC result as grm learn --exp-pick writes it.

    outer(f, li) -> (tp, tn); inner_calls(f, li) -> [(tp, tn)] * 3 or None
    for a failed fit; no_pick holds the folds with no pick.
    """
    outer = outer or outer_default
    inner_calls = inner_calls or inner_default
    folds, sums = [], [dict(n=0, tp=0, fp=0, tn=0, fn=0) for _ in LEARNERS]
    for f in range(FOLDS):
        out = []
        for li, l in enumerate(LEARNERS):
            c = calls(*outer(f, li), OUTER_POS, OUTER_NEG)
            out.append(dict(learner=l, **c))
            for k in ("n", "tp", "fp", "tn", "fn"):
                sums[li][k] += c[k]
        record = {"fold": f, "train": 100, "train_resistant": 40, "test": OUTER_POS + OUTER_NEG,
                  "test_resistant": OUTER_POS, "outer": out}
        if inner == 0:
            record["inner"] = {"status": "not run", "reason": None, "inner_folds": [], "learners": [], "tied": [],
                               "pick": None}
        elif f in no_pick:
            record["inner"] = {"status": "no pick", "reason": "inner fit 0 has no candidate feature",
                               "inner_folds": [], "learners": [], "tied": [], "pick": None}
        else:
            parts = [{"train": 100 - p - q, "train_resistant": 40 - p, "test": p + q, "test_resistant": p}
                     for p, q in zip(INNER_POS, INNER_NEG)]
            rows, scores = [], []
            for li, l in enumerate(LEARNERS):
                got = inner_calls(f, li)
                if got is None:
                    rows.append({"learner": l, "folds": [], "tp": 0, "fp": 0, "tn": 0, "fn": 0,
                                 "balanced_accuracy": None, "null_fits": 0, "error": "inner fit 1: synthetic"})
                    scores.append(None)
                    continue
                cs = [calls(tp, tn, p, q) for (tp, tn), p, q in zip(got, INNER_POS, INNER_NEG)]
                s = {k: sum(c[k] for c in cs) for k in ("tp", "fp", "tn", "fn")}
                ba = balanced(s)
                rows.append({"learner": l, "folds": cs, **s, "balanced_accuracy": ba, "null_fits": 0, "error": None})
                scores.append(ba)
            have = [round(s, 12) for s in scores if s is not None]
            best = max(have)
            tied = [l for l, s in zip(LEARNERS, scores) if s is not None and round(s, 12) == best]
            record["inner"] = {"status": "picked", "reason": None, "inner_folds": parts, "learners": rows,
                               "tied": tied, "pick": tied[-1] if tie == "later" else tied[0]}
        folds.append(record)
    return {
        "format": "grm-c9-pick-1", "evaluation_revision": 14, "dataset": dataset, "scheme": scheme, "seed": seed,
        "inner_folds": inner, "tie": tie, "min_support": 3, "samples": 150, "tested": 150, "tested_resistant": 60,
        "unavailable_reason": None,
        "default": None if scheme != "lineage" else
        {"learner": default[0], "random_balanced_accuracy": default[1], "tied": default[2]},
        "learners": [{"learner": l, "balanced_accuracy": balanced(s)} for l, s in zip(LEARNERS, sums)],
        "folds": folds, "seconds_inner": 1.0,
    }


def release(d, dataset=DATASET):
    """A SYNTHETIC release report that a seed-42 lineage result `d` agrees with."""
    rnd = [{"learner": l, "metrics": {"balanced_accuracy": d["default"]["random_balanced_accuracy"]
                                      if l == d["default"]["learner"] else 0.7}} for l in LEARNERS]
    lin = [{"learner": l["learner"], "metrics": {"balanced_accuracy": l["balanced_accuracy"]}} for l in d["learners"]]
    return {"id": dataset, "species": "Synthia coli", "antibiotic": dataset.split("__")[1].replace("_", " "),
            "splits": [{"scheme": "random", "learners": rnd}, {"scheme": "lineage", "learners": lin}]}


def write_run(top, docs, expect_rows, releases):
    """Lay out a SYNTHETIC run: top/expect.tsv, top/picks/<file>, top/release/datasets/<dataset>.json and
    top/release/index.json."""
    top = pathlib.Path(top)
    (top / "picks").mkdir(parents=True, exist_ok=True)
    (top / "release" / "datasets").mkdir(parents=True, exist_ok=True)
    cols = ["dataset", "species_slug", "split", "scheme", "seed", "inner", "role", "file"]
    with open(top / "expect.tsv", "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in expect_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    for name, d in docs.items():
        with open(top / "picks" / name, "w") as f:
            json.dump(d, f)
    for dataset, rep in releases.items():
        with open(top / "release" / "datasets" / f"{dataset}.json", "w") as f:
            json.dump(rep, f)
    # The release index: the learner served is the random split's best, a tie to the later learner.
    index = []
    for dataset, rep in releases.items():
        rnd = [(l["metrics"]["balanced_accuracy"], i, l["learner"]) for s in rep["splits"] if s["scheme"] == "random"
               for i, l in enumerate(s["learners"])]
        index.append({"id": dataset, "best_learner": max(rnd)[2] if rnd else None})
    with open(top / "release" / "index.json", "w") as f:
        json.dump({"datasets": index}, f)
    return top / "expect.tsv", top / "picks", top / "release", top / "score"


def full_run(datasets=(DATASET,), study=(DATASET,), seeds=(42, 43, 44), twin=(45, 46, 47), change=None):
    """SYNTHETIC expect rows, results and release reports of one species' whole run.

    change(row, doc) may alter a result before it is written.
    """
    expect = load("expect")
    rows = expect.rows(SPECIES, list(datasets), list(study), list(seeds), list(twin))
    docs, releases = {}, {}
    for r in rows:
        dataset, _, split, scheme, seed, inner, role, name = r
        d = doc(dataset=dataset, scheme=scheme, seed=seed, inner=inner)
        if change:
            change(r, d)
        docs[name] = d
        if split == "lineage" and seed == 42:
            releases[dataset] = release(d, dataset)
    return rows, docs, releases
