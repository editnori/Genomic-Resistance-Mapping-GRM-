"""Shared module loaders and synthetic fixture builders for tests/point_mutations."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS_DIR = ROOT / "scripts" / "point-mutations"
SCAN_SCRIPT = ROOT / "scripts" / "c7-outside" / "scan.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
INDEX_DOC = ROOT / "docs" / "INDEX.md"
FOUR_TRACK_DOC = ROOT / "docs" / "four-track-plan.md"
PM_PLAN_DOC = ROOT / "docs" / "point-mutation-plan.md"
SUITE = "tests/point_mutations"


def load_module(short_name: str, path: Path | None = None):
    key = f"grm_point_mutations.{short_name}"
    if key in sys.modules:
        return sys.modules[key]
    target = path if path is not None else SCRIPTS_DIR / f"{short_name}.py"
    spec = importlib.util.spec_from_file_location(key, target)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {target}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


EA = load_module("exposure_audit")
FT = load_module("features")
CMP = load_module("compare")
C7_SCAN = load_module("c7_scan", SCAN_SCRIPT)


def sample_exposure_audit(
    *,
    model_family: str = "esm2",
    sequence_overlap_status: str = "unknown",
    phenotype_annotation_status: str = "none_verified",
    eval_tuning_status: str = "frozen_preregistered",
    checked_at_utc: str = "unverified",
    claimed_ceiling: str | None = None,
) -> dict[str, Any]:
    record: dict[str, Any] = {
        "format": "grm-rep-exposure-v1",
        "model_family": model_family,
        "checkpoint_id": f"synthetic/{model_family}_checkpoint_v1",
        "checkpoint_sha256": hashlib.sha256(f"{model_family}-weights".encode("utf-8")).hexdigest(),
        "training_corpus": {
            "name": "SyntheticBacterialCorpus",
            "version_or_cutoff": "2024-01",
            "sources": ["UniRef50", "RefSeq_Prokaryotes"],
            "taxonomic_scope": "prokaryotic_proteins",
        },
        "sequence_overlap_status": sequence_overlap_status,
        "phenotype_annotation_status": phenotype_annotation_status,
        "eval_tuning_status": eval_tuning_status,
        "audit_evidence": {
            "method": "synthetic_corpus_audit",
            "checked_at_utc": checked_at_utc,
            "notes": "Synthetic unit test audit record.",
        },
    }
    if claimed_ceiling is not None:
        record["claimed_ceiling"] = claimed_ceiling
    return record


def _pack_binary_matrix(X_bin: np.ndarray) -> bytes:
    """Pack samples x features (0/1 uint8) into GRM rows x words little-endian u64 bytes."""
    n_samples, n_features = X_bin.shape
    words = (n_samples + 63) // 64
    rows_bytes = bytearray()
    for r in range(n_features):
        col = X_bin[:, r].astype(np.uint8)
        packed = np.packbits(col, bitorder="little")
        buf = np.zeros(words * 8, dtype=np.uint8)
        buf[: packed.shape[0]] = packed
        rows_bytes.extend(buf.tobytes())
    return bytes(rows_bytes)


def build_synthetic_dataset(
    base_dir: Path,
    *,
    dataset_id: str = "synthetic_species__ciprofloxacin",
    n_samples: int = 45,
    n_kmers: int = 80,
    seed: int = 42,
    n_outer_folds: int = 3,
    n_inner_folds: int = 3,
) -> dict[str, Any]:
    """Build a complete synthetic `grm-dataset-export-v1` + `grm-fold-file-v1` +
    protein variants + foundation manifest fixture in `base_dir`."""
    rng = np.random.default_rng(20261127 + seed)
    samples = [f"573.{1001 + i}" for i in range(n_samples)]
    # Construct 9 balanced lineages (5 samples each) with both classes per outer fold.
    lineages = [f"L{i // 5}" for i in range(n_samples)]
    outer = np.array([(i // 5) % n_outer_folds for i in range(n_samples)], dtype=np.int64)
    # Assign balanced resistant labels so every outer and inner fold holds both classes.
    y = np.array([(i % 2 == 0) for i in range(n_samples)], dtype=bool)
    # Flip a deterministic subset to prevent trivial exact alternation.
    flip_idx = [3, 12, 21, 30, 39]
    for idx in flip_idx:
        if idx < n_samples:
            y[idx] = not y[idx]

    # Build synthetic 31-mer matrix with predictive + noisy columns.
    X_bin = np.zeros((n_samples, n_kmers), dtype=np.uint8)
    for c in range(n_kmers):
        if c < 12:
            prob = np.where(y, 0.80, 0.15)
        elif c < 24:
            prob = np.where(y, 0.20, 0.80)
        else:
            prob = np.full(n_samples, 0.45)
        X_bin[:, c] = (rng.random(n_samples) < prob).astype(np.uint8)

    # Ensure deterministic inner folds inside each outer fold keeping lineages whole.
    inner_plans: list[list[int | None]] = []
    for f in range(n_outer_folds):
        tr_lineages = sorted({lineages[i] for i in range(n_samples) if outer[i] != f})
        lin_to_inner = {lin: idx % n_inner_folds for idx, lin in enumerate(tr_lineages)}
        inner_row: list[int | None] = [
            None if outer[i] == f else lin_to_inner[lineages[i]] for i in range(n_samples)
        ]
        inner_plans.append(inner_row)

    plan_body = {
        "scheme": "lineage",
        "seed": seed,
        "repeat": 0,
        "repeat_seed": seed,
        "k": n_outer_folds,
        "outer": [int(v) for v in outer],
        "inner": inner_plans,
        "holdouts": [{"gbt": [], "lasso": []} for _ in range(n_outer_folds)],
        "unavailable_reason": None,
    }
    plan_body["sha256"] = CMP.canonical_plan_sha256(plan_body)

    export_dir = base_dir / dataset_id
    export_dir.mkdir(parents=True, exist_ok=True)
    mbytes = _pack_binary_matrix(X_bin)
    (export_dir / "matrix.bin").write_bytes(mbytes)
    matrix_sha = hashlib.sha256(mbytes).hexdigest()

    words = (n_samples + 63) // 64
    meta = {
        "format": "grm-dataset-export-v1",
        "evaluation_revision": 14,
        "id": dataset_id,
        "input_fingerprint": hashlib.sha256(f"fp-{dataset_id}-{seed}".encode("utf-8")).hexdigest(),
        "screened_features": 4000,
        "folds": n_outer_folds,
        "inner_folds": n_inner_folds,
        "seed": seed,
        "min_support": 3,
        "samples": samples,
        "resistant": [bool(v) for v in y],
        "lineage": lineages,
        "duplicates": [],
        "rows": n_kmers,
        "cols": n_samples,
        "words": words,
        "matrix_sha256": matrix_sha,
        "splits": [
            {
                "plan": plan_body,
                "screens": [
                    {
                        "fold": f,
                        "n_train": int(np.sum(outer != f)),
                        "n_candidates": n_kmers,
                        "screened_rows": list(range(min(n_kmers, 4000))),
                    }
                    for f in range(n_outer_folds)
                ],
            }
        ],
    }
    (export_dir / "dataset.json").write_text(json.dumps(meta, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    # Build a matching study-whole grm-fold-file-v1 (group samples by contiguous study blocks).
    study_test = [
        [samples[i] for i in range(n_samples) if (i // 5) % n_outer_folds == f]
        for f in range(n_outer_folds)
    ]
    fold_doc = {
        "format": "grm-fold-file-v1",
        "dataset": dataset_id,
        "scheme": "study-whole",
        "test": study_test,
        "train": [],
        "unused": [],
    }
    fold_file = base_dir / f"{dataset_id}.study-whole.json"
    fold_file.write_text(json.dumps(fold_doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    # Build synthetic protein variants on frozen target genes + regions.
    target_genes = ["gyrA", "parC", "rpoB"]
    regions = [
        {"gene": "gyrA", "region": "QRDR", "start": 67, "end": 106},
        {"gene": "parC", "region": "QRDR", "start": 64, "end": 103},
        {"gene": "rpoB", "region": "RRDR", "start": 426, "end": 452},
    ]
    isolates_variants: list[list[dict[str, Any]]] = []
    for i in range(n_samples):
        vlist: list[dict[str, Any]] = []
        if y[i]:
            if rng.random() < 0.75:
                vlist.append({"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "missense"})
            if rng.random() < 0.45:
                vlist.append({"gene": "parC", "change": "S80I", "pos": 80, "kind": "missense"})
            if rng.random() < 0.25:
                vlist.append({"gene": "gyrA", "change": "D87N", "pos": 87, "kind": "missense"})
            if rng.random() < 0.20:
                vlist.append({"gene": "rpoB", "change": "Q432*", "pos": 432, "kind": "nonsense"})
        else:
            if rng.random() < 0.15:
                vlist.append({"gene": "gyrA", "change": "S83L", "pos": 83, "kind": "missense"})
            if rng.random() < 0.50:
                vlist.append({"gene": "gyrA", "change": "A90A", "pos": 90, "kind": "synonymous"})
            if rng.random() < 0.20:
                vlist.append({"gene": "parC", "change": "del_75_77", "pos_start": 75, "pos_end": 77, "kind": "inframe_deletion"})
        # Add an isolate-specific rare variant on some folds to exercise unseen_variant_count.
        if i % 7 == 0:
            vlist.append({"gene": "gyrA", "change": f"V{200 + i}A", "pos": 200 + i, "kind": "missense"})
        isolates_variants.append(vlist)

    variants_doc = {
        "dataset": dataset_id,
        "target_genes": target_genes,
        "regions": regions,
        "isolates": {samples[i]: isolates_variants[i] for i in range(n_samples)},
    }
    variants_file = base_dir / f"{dataset_id}.variants.json"
    variants_file.write_text(json.dumps(variants_doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    # Build a synthetic foundation-model embedding manifest (permuted sample order to test alignment!).
    n_emb = 8
    emb = rng.normal(0.0, 1.0, size=(n_samples, n_emb))
    emb[y, 0] += 1.4
    emb[y, 1] -= 0.9
    perm = rng.permutation(n_samples)
    foundation_doc = {
        "format": "grm-rep-foundation-v1",
        "dataset": dataset_id,
        "representation_kind": "embedding",
        "exposure_audit": sample_exposure_audit(
            model_family="esm2",
            sequence_overlap_status="unknown",
            phenotype_annotation_status="none_verified",
            eval_tuning_status="frozen_preregistered",
        ),
        "samples": [samples[int(p)] for p in perm],
        "feature_names": [f"esm2_dim_{d}" for d in range(n_emb)],
        "values": [[round(float(v), 6) for v in emb[int(p)]] for p in perm],
    }
    foundation_file = base_dir / f"{dataset_id}.esm2.json"
    foundation_file.write_text(json.dumps(foundation_doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    # Species JSON file shaped for scripts/c7-outside/scan.py hygiene check.
    species_file = base_dir / "synthetic_species.json"
    species_file.write_text(
        json.dumps({"species": "synthetic_species", "genomes": [{"id": s} for s in samples]}) + "\n",
        encoding="utf-8",
    )

    return {
        "dataset_id": dataset_id,
        "export_dir": export_dir,
        "export_meta": meta,
        "fold_file": fold_file,
        "fold_doc": fold_doc,
        "variants_file": variants_file,
        "isolates_variants": isolates_variants,
        "target_genes": target_genes,
        "regions": regions,
        "foundation_file": foundation_file,
        "foundation_doc": foundation_doc,
        "species_file": species_file,
        "samples": samples,
        "y": y,
    }
