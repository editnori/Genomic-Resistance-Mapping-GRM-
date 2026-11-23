"""SYNTHETIC exports, unitig builds, fold files and results for the outside-learner scripts (scripts/c7-outside,
docs/outside-learners.md).

Every matrix, label, lineage, plan, fold file, report and count here is made up in the test. No genome or call
comes from a real dataset, the dataset names are made up, and the sample names are not genome ids.
"""

import hashlib
import importlib.util
import json
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "c7-outside"
SPECIES = "synthia_coli"
DATASET = "synthia_coli__drug_a"
OTHER = "synthia_coli__drug_b"


def load(name):
    spec = importlib.util.spec_from_file_location(f"c7_outside_{name}", KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


O = load("outside")


def names(n, prefix="s"):
    return [f"{prefix}{i:04d}" for i in range(n)]


def table(n=150, rows=160, informative=6, seed=7, lineages=15, noise=0.08):
    """A SYNTHETIC 0/1 table: rows x n presence, labels from a few rows with noise, and lineage ids.
    Rows are distinct patterns, as a projected matrix's are."""
    rng = np.random.default_rng(seed)
    lineage = rng.integers(0, lineages, n)
    x = rng.random((rows, n)) < rng.uniform(0.1, 0.6, (rows, 1))
    # Rows that follow the lineages, as real ones do.
    for r in range(rows // 4):
        x[r] = np.isin(lineage, rng.choice(lineages, 3, replace=False))
    score = x[:informative].sum(axis=0) + rng.normal(0, 0.6, n)
    y = score > np.median(score)
    y ^= rng.random(n) < noise
    _, first = np.unique(x, axis=0, return_index=True)
    x = x[np.sort(first)]
    return x, y, lineage.tolist()


def packed(x):
    """Rows of booleans (k x n) as k x words little-endian 64-bit words."""
    return np.stack([O.words_of(r) for r in x]) if x.shape[0] else np.zeros((0, (x.shape[1] + 63) // 64), dtype="<u8")


def plan(scheme, y, lineage, duplicates, seed, k=5):
    """A SYNTHETIC fold plan in the export's form: outer folds, 3 inner folds and both holdouts per fold."""
    n = y.shape[0]
    rng = np.random.default_rng(seed + (1 if scheme == "lineage" else 2))
    groups = O.tuning_groups(scheme, lineage, duplicates, n)
    if scheme == "lineage":
        outer = O.grouped_folds(y, lineage, k, rng)
    else:
        outer = O.stratified_folds(y, k, rng)
        for g in duplicates:
            outer[g] = outer[g[0]]
    inner, holdouts = [], []
    for f in range(k):
        train = outer != f
        fi = O.training_folds(y, train, groups, 3, rng)
        inner.append([None if v < 0 else int(v) for v in fi])
        _, hg = O.tuning_holdout(train, y, groups, 0.15, rng)
        _, hl = O.tuning_holdout(train, y, groups, 0.2, rng)
        holdouts.append({"gbt": np.flatnonzero(hg).tolist(), "lasso": np.flatnonzero(hl).tolist()})
    p = {"format": "grm-fold-plan-v1", "scheme": scheme, "repeat": 0, "k": k, "inner_k": 3, "seed": seed,
         "repeat_seed": seed, "outer": [int(v) for v in outer], "inner": inner, "holdouts": holdouts,
         "unavailable_reason": None}
    p["sha256"] = O.canonical_sha256(p)
    return p


def write_export(top, dataset=DATASET, n=150, rows=160, seed=42, table_seed=7, duplicates=((3, 4),), flip=None,
                 change=None):
    """A SYNTHETIC `grm export` directory, top/<dataset>/, with the forest's screens computed by the export's rule.
    flip: (scheme, fold) whose test genomes' labels are turned round after the plans and screens are made, so the
    plans and screens are those of the true labels, as the export's are of its training genomes only."""
    x, y, lineage = table(n, rows, seed=table_seed)
    duplicates = [list(g) for g in duplicates]
    for g in duplicates:
        x[:, g] = x[:, [g[0]]]
    _, first = np.unique(x, axis=0, return_index=True)
    x = x[np.sort(first)]
    words = packed(x)
    splits = []
    for scheme in ("random", "lineage"):
        p = plan(scheme, y, lineage, duplicates, seed)
        outer = np.array(p["outer"])
        screens = []
        for f in range(p["k"]):
            train = outer != f
            n_cand, tops = O.screen_fold(words, train, y, 3, {"rf": train})
            screens.append({"fold": f, "n_train": int(train.sum()), "n_candidates": n_cand,
                            "screened_rows": [int(r) for r in tops["rf"]]})
        splits.append({"plan": p, "screens": screens})
    labels = y.copy()
    if flip is not None:
        scheme, f = flip
        outer = np.array(next(s["plan"]["outer"] for s in splits if s["plan"]["scheme"] == scheme))
        labels[outer == f] = ~labels[outer == f]
    d = pathlib.Path(top) / dataset
    d.mkdir(parents=True, exist_ok=True)
    raw = np.ascontiguousarray(words).tobytes()
    (d / "matrix.bin").write_bytes(raw)
    meta = {"format": "grm-dataset-export-v1", "evaluation_revision": 14, "id": dataset, "species_slug": dataset.split("__")[0],
            "antibiotic": dataset.split("__")[1], "input_fingerprint": "synthetic-" + dataset, "screened_features": 4000,
            "folds": 5, "inner_folds": 3, "seed": seed, "min_support": 3, "samples": names(n),
            "resistant": [bool(v) for v in labels], "lineage": lineage, "duplicates": duplicates,
            "rows": int(words.shape[0]), "cols": n, "words": int(words.shape[1]),
            "matrix_sha256": hashlib.sha256(raw).hexdigest(), "splits": splits}
    if change:
        change(meta)
    (d / "dataset.json").write_text(json.dumps(meta))
    return d, x, y


def write_unitigs(top, n=150, extra=10, rows=220, seed=11, table_seed=7):
    """A SYNTHETIC unitig build's matrix folder: the cohort's genomes and `extra` more, patterns over all of them,
    some of which are constant or rare on the cohort, as a species build's are."""
    x, _, _ = table(n, 160, seed=table_seed)
    rng = np.random.default_rng(seed)
    more = rng.random((rows, n)) < rng.uniform(0.02, 0.7, (rows, 1))
    cohort = np.vstack([x[:40], more, x[:5]])  # repeats collapse in the projection
    tail = rng.random((cohort.shape[0], extra)) < 0.5
    pats = np.hstack([cohort, tail])
    genomes = names(n) + names(extra, prefix="o")
    words = packed(pats)
    d = pathlib.Path(top) / "unitigs" / "matrix"
    d.mkdir(parents=True, exist_ok=True)
    (d / "patterns.bin").write_bytes(np.ascontiguousarray(words).tobytes())
    (d / "matrix.json").write_text(json.dumps({"genomes": genomes, "n_patterns": int(words.shape[0]), "k": 31}))
    return d


def write_fold_file(top, dataset=DATASET, n=150, k=3, unused=(0, 1), seed=5, table_seed=7):
    """A SYNTHETIC study-whole fold file: k test folds of whole lineages, some genomes always trained, a few unused."""
    _, y, lineage = table(n, 160, seed=table_seed)
    rng = np.random.default_rng(seed)
    ids = names(n)
    lin = sorted(set(lineage))
    rng.shuffle(lin)
    test_lin = lin[: len(lin) * 3 // 4]
    folds = [[] for _ in range(k)]
    train, gone = [], []
    for i in range(n):
        if i in unused:
            gone.append(ids[i])
        elif i in (3, 4):  # write_export's duplicate pair stays on one side
            train.append(ids[i])
        elif lineage[i] in test_lin:
            folds[test_lin.index(lineage[i]) % k].append(ids[i])
        else:
            train.append(ids[i])
    doc = {"format": "grm-fold-file-v1", "dataset": dataset, "scheme": "study-whole", "test": folds, "train": train,
           "unused": gone}
    path = pathlib.Path(top) / f"{dataset}.study-whole.json"
    path.write_text(json.dumps(doc))
    return path


# ---------------------------------------------------------------- results for the scorer

def outside_report(dataset=DATASET, scheme="lineage", seed=42, arm="xgboost", shape=None, calls=None):
    """A SYNTHETIC grm-c7-outside-1 report with the given folds by class and calls."""
    shape = shape or [(120, 60, 30, 15)] * 5
    calls = calls or [{"n": t, "tp": tr - 2, "fp": 3, "tn": t - tr - 3, "fn": 2, "loss_sum": 0.4 * t} for _, _, t, tr in shape]
    folds = []
    for f, ((n, nr, t, tr), c) in enumerate(zip(shape, calls)):
        row = {"fold": f, "train": n, "train_resistant": nr, "test": t, "test_resistant": tr, "fitting": n - 18,
               "fitting_resistant": nr - 9, "holdout": 18, "holdout_resistant": 9, "candidates": 900, "screened": 900,
               "constant": 0, "calls": c}
        if arm == "elasticnet":
            row.update(grid_points=40, grid_unconverged=0, picked_mix=0.5, picked_penalty=3, refit_converged=1, nonzero=12,
                       null_fit=0)
        else:
            row.update(rounds_run=60, rounds_kept=35)
        folds.append(row)
    return {"format": "grm-c7-outside-1", "evaluation_revision": 14, "dataset": dataset, "scheme": scheme, "seed": seed,
            "arm": arm, "pairs_with": O.ARMS[arm], "library": O.LIBRARY[arm], "version": O.VERSIONS[O.LIBRARY[arm]],
            "settings": {}, "input": "unitigs" if arm == "elasticnet" else "export", "plan_sha256": None,
            "input_fingerprint": "synthetic", "matrix_sha256": "0" * 64, "fold_file_sha256": None, "samples": 150,
            "resistant": 75, "min_support": 3, "input_rows": 900, "unavailable_reason": None,
            "tested": sum(s[2] for s in shape), "tested_resistant": sum(s[3] for s in shape), "folds": folds,
            "constant_folds": 0, "one_tree_folds": 0, "null_fit_folds": 0, "unconverged_refits": 0}
