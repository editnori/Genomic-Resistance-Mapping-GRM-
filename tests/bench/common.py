"""Shared SYNTHETIC fixture builders for tests/bench/.

All genomes, FASTA sequences, feature matrices, fold files, and reports built
here are purely synthetic (no real data, no private paths, no real accessions).
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from scripts.bench.protocol import FOLD_FILE_FORMAT, canonical_json_bytes


def pack_bitmatrix(X_rows: np.ndarray) -> tuple[bytes, int, str]:
    """Pack a (rows, cols) 0/1 uint8 array into GRM matrix.bin little-endian u64 words."""
    rows, cols = X_rows.shape
    words = (cols + 63) // 64
    out = np.zeros((rows, words), dtype="<u8")
    for r in range(rows):
        for c in range(cols):
            if X_rows[r, c]:
                out[r, c >> 6] |= np.uint64(1) << np.uint64(c & 63)
    raw = out.tobytes()
    return raw, words, hashlib.sha256(raw).hexdigest()


def canonical_plan_sha256(plan: dict[str, Any]) -> str:
    body = {k: v for k, v in plan.items() if k != "sha256"}
    return hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    ).hexdigest()


def build_synthetic_export(
    base_dir: Path,
    *,
    species_slug: str = "escherichia_coli",
    antibiotic: str = "ciprofloxacin",
    n_samples: int = 100,
    n_features: int = 32,
    seed: int = 42,
    k_folds: int = 5,
    inner_folds: int = 3,
    min_support: int = 3,
    screened_features: int = 4000,
) -> tuple[Path, dict[str, Any], np.ndarray]:
    """Build a valid SYNTHETIC `grm export` directory (dataset.json + matrix.bin).

    Constructs balanced labels across lineages so every outer fold, fitting side,
    boosting holdout, and inner fold contains both R and S classes.
    """
    drug_part = antibiotic.replace("/", "_").replace("-", "_").replace(" ", "_")
    ds_id = f"{species_slug}__{drug_part}"
    ds_dir = base_dir / ds_id
    ds_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(seed)
    # Synthetic BV-BRC-shaped IDs (used only in temp dirs to verify scan.py catches leaks)
    samples = [f"999000.{1000 + i}" for i in range(n_samples)]
    # 20 lineages of 5 samples each (or n_samples // 4), alternating labels inside each lineage
    n_lin = max(k_folds * 3, n_samples // 4)
    lineage = [i % n_lin for i in range(n_samples)]
    resistant = [(i // 2) % 2 == 0 for i in range(n_samples)]
    y_arr = np.array(resistant, dtype=bool)

    # Build synthetic binary feature matrix: row 0 = noisy signal, row 1 = strong signal, rest random
    X_rows = np.zeros((n_features, n_samples), dtype=np.uint8)
    for c in range(n_samples):
        X_rows[0, c] = int(resistant[c]) if rng.random() < 0.85 else int(not resistant[c])
        X_rows[1, c] = int(resistant[c]) if rng.random() < 0.75 else int(not resistant[c])
    for r in range(2, n_features):
        p = 0.2 + 0.6 * (r / max(n_features, 1))
        X_rows[r] = (rng.random(n_samples) < p).astype(np.uint8)

    raw_bytes, words, mat_sha = pack_bitmatrix(X_rows)
    (ds_dir / "matrix.bin").write_bytes(raw_bytes)

    def _make_split(scheme: str) -> dict[str, Any]:
        if scheme == "lineage":
            outer = [int(lin % k_folds) for lin in lineage]
        else:
            outer = [int(i % k_folds) for i in range(n_samples)]
        holdouts = []
        inner_list = []
        screens = []
        for f in range(k_folds):
            train_idx = [i for i in range(n_samples) if outer[i] != f]
            train_r = [i for i in train_idx if resistant[i]]
            train_s = [i for i in train_idx if not resistant[i]]
            # Holdout: ~15% of each class inside train_idx
            hr = max(2, len(train_r) // 6)
            hs = max(2, len(train_s) // 6)
            gbt_hold = sorted(train_r[:hr] + train_s[:hs])
            holdouts.append({"gbt": gbt_hold})

            # Inner folds strictly inside train_idx, stratified by class
            inner_arr: list[int | None] = [None] * n_samples
            for pos, idx in enumerate(train_r):
                inner_arr[idx] = pos % inner_folds
            for pos, idx in enumerate(train_s):
                inner_arr[idx] = pos % inner_folds
            inner_list.append(inner_arr)

            # Screen rows on train_idx only
            train_sub = X_rows[:, train_idx]
            sup = train_sub.sum(axis=1)
            n_tr = len(train_idx)
            cands = [r for r in range(n_features) if min_support <= sup[r] <= n_tr - min_support]
            screens.append({
                "fold": f,
                "n_train": n_tr,
                "n_candidates": len(cands),
                "screened_rows": cands[:screened_features],
            })
        plan: dict[str, Any] = {
            "scheme": scheme,
            "seed": seed,
            "repeat": 0,
            "repeat_seed": seed,
            "k": k_folds,
            "outer": outer,
            "holdouts": holdouts,
            "inner": inner_list,
            "unavailable_reason": None,
        }
        plan["sha256"] = canonical_plan_sha256(plan)
        return {"plan": plan, "screens": screens}

    fp = hashlib.sha256(f"synthetic-fp:{ds_id}:{seed}:{mat_sha}".encode("utf-8")).hexdigest()[:32]
    meta = {
        "format": "grm-dataset-export-v1",
        "evaluation_revision": 14,
        "id": ds_id,
        "species_slug": species_slug,
        "antibiotic": antibiotic,
        "input_fingerprint": fp,
        "screened_features": screened_features,
        "folds": k_folds,
        "inner_folds": inner_folds,
        "seed": seed,
        "min_support": min_support,
        "samples": samples,
        "resistant": resistant,
        "lineage": lineage,
        "duplicates": [],
        "rows": n_features,
        "cols": n_samples,
        "words": words,
        "matrix_sha256": mat_sha,
        "splits": [_make_split("random"), _make_split("lineage")],
    }
    (ds_dir / "dataset.json").write_bytes(canonical_json_bytes(meta))
    return ds_dir, meta, X_rows


def build_synthetic_fold_file(
    path: Path,
    meta: dict[str, Any],
    *,
    scheme: str = "study-whole",
    n_test_folds: int = 2,
    n_unused: int = 4,
) -> Path:
    """Build a valid SYNTHETIC `grm-fold-file-v1` partitioning all samples of `meta`."""
    samples = list(meta["samples"])
    resistant = list(meta["resistant"])
    unused = samples[:n_unused]
    rem = samples[n_unused:]
    rem_r = [g for g, y in zip(samples[n_unused:], resistant[n_unused:]) if y]
    rem_s = [g for g, y in zip(samples[n_unused:], resistant[n_unused:]) if not y]

    if n_test_folds == 1:
        # e.g., time split: 30% test, 70% train
        nr_te = max(4, len(rem_r) // 3)
        ns_te = max(4, len(rem_s) // 3)
        test = [rem_r[:nr_te] + rem_s[:ns_te]]
        train = rem_r[nr_te:] + rem_s[ns_te:]
    else:
        test = [[] for _ in range(n_test_folds)]
        train = []
        for idx, g in enumerate(rem_r):
            test[idx % n_test_folds].append(g)
        for idx, g in enumerate(rem_s):
            test[idx % n_test_folds].append(g)

    doc = {
        "format": FOLD_FILE_FORMAT,
        "dataset": meta["id"],
        "scheme": scheme,
        "test": test,
        "train": train,
        "unused": unused,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(canonical_json_bytes(doc))
    return path


def build_matching_reference_docs(meta: dict[str, Any], scheme: str = "lineage") -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    """Build synthetic C3 pick result, grm learn report, and release report matching `meta`."""
    split = next(s for s in meta["splits"] if s["plan"]["scheme"] == scheme)
    plan = split["plan"]
    outer = plan["outer"]
    y = meta["resistant"]
    folds = []
    tested = 0
    tested_r = 0
    for f in range(plan["k"]):
        tr = [i for i, o in enumerate(outer) if o != f]
        te = [i for i, o in enumerate(outer) if o == f]
        tr_r = sum(1 for i in tr if y[i])
        te_r = sum(1 for i in te if y[i])
        folds.append({
            "fold": f,
            "train": len(tr),
            "train_resistant": tr_r,
            "test": len(te),
            "test_resistant": te_r,
        })
        tested += len(te)
        tested_r += te_r
    ref_pick = {
        "unavailable_reason": None,
        "folds": folds,
        "tested": tested,
        "tested_resistant": tested_r,
    }
    ref_report = {
        "eval_config": {"seed": meta["seed"]},
        "input_fingerprint": meta["input_fingerprint"],
        "samples": [{"id": g} for g in meta["samples"]],
        "splits": [
            {"scheme": scheme, "unavailable_reason": None, "test_fold": outer},
        ],
    }
    release_report = {
        "id": meta["id"],
        "input_fingerprint": meta["input_fingerprint"],
    }
    return ref_pick, ref_report, release_report
