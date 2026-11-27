"""Pretraining and evaluation exposure audit schema (`grm-rep-exposure-v1`).

Every pretrained foundation-model representation arm evaluated in Track B must
carry a validated `grm-rep-exposure-v1` audit record. Validation checks:

1. Model family eligibility: bacterial/prokaryotic sequence or protein models
   (`esm2`, `evo2`, `prokbert`, `bacformer`) are admitted; mammalian eukaryotic
   regulatory models (`alphagenome2`) are explicitly refused for bacterial AMR.
2. Checkpoint identity and training corpus provenance (SHA-256 pinned).
3. Three independent exposure axes:
   - `sequence_overlap_status`: whether evaluation genomes/proteins appeared in
     unsupervised pretraining corpora;
   - `phenotype_annotation_status`: whether AMR phenotypes or curated mutation
     catalogues influenced model weights;
   - `eval_tuning_status`: whether checkpoint, layer, prompt, or pooling choices
     were touched after inspecting evaluation outcomes.
4. Deterministic claim-ceiling derivation: unknown or overlapping sequence
   exposure caps claims at `retrospective_only` (never
   `certified_unseen_sequence`); phenotype/catalogue exposure caps at
   `non_independent_retrospective`; evaluation-time tuning caps at
   `contaminated_exploratory`.

CLI usage:
  python3 -I scripts/point-mutations/exposure_audit.py validate AUDIT.json
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

FORMAT = "grm-rep-exposure-v1"

ALLOWED_MODEL_FAMILIES = frozenset({"esm2", "evo2", "prokbert", "bacformer"})

EXCLUDED_MODEL_FAMILIES = {
    "alphagenome2": (
        "AlphaGenome 2 is excluded for bacterial genomes: its mammalian eukaryotic "
        "chromatin, splicing, and regulatory track domain is incompatible with "
        "prokaryotic coding point-mutation resistance mechanisms."
    ),
    "alphagenome": (
        "AlphaGenome is excluded for bacterial genomes: mammalian eukaryotic "
        "regulatory modeling is incompatible with prokaryotic AMR."
    ),
}

DISALLOWED_TAXONOMIC_SCOPES = frozenset({"eukaryotic_only", "mammalian_only", "human_only"})

SEQUENCE_OVERLAP_STATUSES = frozenset({
    "verified_disjoint",
    "partial_overlap",
    "confirmed_overlap",
    "unknown",
})

PHENOTYPE_ANNOTATION_STATUSES = frozenset({
    "none_verified",
    "indirect_catalogue",
    "supervised_phenotype",
    "unknown",
})

EVAL_TUNING_STATUSES = frozenset({
    "frozen_preregistered",
    "outer_train_only",
    "eval_touched",
    "unknown",
})

CLAIM_CEILING_ORDER = (
    "contaminated_exploratory",
    "non_independent_retrospective",
    "retrospective_only",
    "certified_unseen_sequence",
)

CLAIM_CEILING_RANK = {name: idx + 1 for idx, name in enumerate(CLAIM_CEILING_ORDER)}

REQUIRED_TOP_KEYS = frozenset({
    "format",
    "model_family",
    "checkpoint_id",
    "checkpoint_sha256",
    "training_corpus",
    "sequence_overlap_status",
    "phenotype_annotation_status",
    "eval_tuning_status",
    "audit_evidence",
})

OPTIONAL_TOP_KEYS = frozenset({"claimed_ceiling"})
ALLOWED_TOP_KEYS = REQUIRED_TOP_KEYS | OPTIONAL_TOP_KEYS

REQUIRED_CORPUS_KEYS = frozenset({
    "name",
    "version_or_cutoff",
    "sources",
    "taxonomic_scope",
})

REQUIRED_EVIDENCE_KEYS = frozenset({
    "method",
    "checked_at_utc",
    "notes",
})

SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
UTC_TS_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z$")


class AuditError(ValueError):
    """Raised when an exposure audit record is missing, malformed, or ineligible."""


def _need(cond: bool, message: str) -> None:
    if not cond:
        raise AuditError(message)


def _normalize_family_key(raw: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", raw.strip().lower())


def _require_nonempty_str(value: Any, field: str, max_len: int = 256) -> str:
    _need(isinstance(value, str), f"{field} must be a string")
    stripped = value.strip()
    _need(bool(stripped) and len(stripped) <= max_len, f"{field} must be a non-empty string (<= {max_len} chars)")
    return stripped


def derive_claim_ceiling(
    sequence_overlap_status: str,
    phenotype_annotation_status: str,
    eval_tuning_status: str,
) -> tuple[str, list[str]]:
    """Derive the maximum permissible interpretive claim and capping reasons."""
    _need(
        sequence_overlap_status in SEQUENCE_OVERLAP_STATUSES,
        f"invalid sequence_overlap_status: {sequence_overlap_status!r}",
    )
    _need(
        phenotype_annotation_status in PHENOTYPE_ANNOTATION_STATUSES,
        f"invalid phenotype_annotation_status: {phenotype_annotation_status!r}",
    )
    _need(
        eval_tuning_status in EVAL_TUNING_STATUSES,
        f"invalid eval_tuning_status: {eval_tuning_status!r}",
    )

    reasons: list[str] = []
    if eval_tuning_status == "eval_touched":
        reasons.append("eval_tuning_status is eval_touched: evaluation outcomes touched model selection")
    elif eval_tuning_status == "unknown":
        reasons.append("eval_tuning_status is unknown: evaluation-time selection provenance unverified")

    if phenotype_annotation_status == "supervised_phenotype":
        reasons.append("phenotype_annotation_status is supervised_phenotype: weights trained on AMR labels")
    elif phenotype_annotation_status == "indirect_catalogue":
        reasons.append("phenotype_annotation_status is indirect_catalogue: resistance catalogue exposure")
    elif phenotype_annotation_status == "unknown":
        reasons.append("phenotype_annotation_status is unknown: phenotype supervision exposure unverified")

    if sequence_overlap_status == "confirmed_overlap":
        reasons.append("sequence_overlap_status is confirmed_overlap: evaluation sequences overlap pretraining")
    elif sequence_overlap_status == "partial_overlap":
        reasons.append("sequence_overlap_status is partial_overlap: partial pretraining sequence overlap")
    elif sequence_overlap_status == "unknown":
        reasons.append("sequence_overlap_status is unknown: pretraining overlap unaudited")

    if eval_tuning_status in ("eval_touched", "unknown"):
        return "contaminated_exploratory", reasons
    if phenotype_annotation_status in ("indirect_catalogue", "supervised_phenotype", "unknown"):
        return "non_independent_retrospective", reasons
    if sequence_overlap_status in ("unknown", "partial_overlap", "confirmed_overlap"):
        return "retrospective_only", reasons
    return "certified_unseen_sequence", reasons


def validate_exposure_audit(record: Any) -> dict[str, Any]:
    """Strictly validate a `grm-rep-exposure-v1` audit dict and return a normalized copy."""
    _need(isinstance(record, dict), "exposure audit record must be a JSON object")

    keys = set(record.keys())
    missing = sorted(REQUIRED_TOP_KEYS - keys)
    extra = sorted(keys - ALLOWED_TOP_KEYS)
    _need(not missing, f"exposure audit missing required keys: {missing}")
    _need(not extra, f"exposure audit contains unlisted keys: {extra}")

    _need(record["format"] == FORMAT, f"exposure audit format must be {FORMAT!r}, got {record['format']!r}")

    raw_family = _require_nonempty_str(record["model_family"], "model_family", max_len=64)
    norm_family = _normalize_family_key(raw_family)
    if norm_family in EXCLUDED_MODEL_FAMILIES:
        raise AuditError(EXCLUDED_MODEL_FAMILIES[norm_family])
    _need(
        raw_family in ALLOWED_MODEL_FAMILIES,
        f"unsupported model_family {raw_family!r}; allowed: {sorted(ALLOWED_MODEL_FAMILIES)}",
    )

    checkpoint_id = _require_nonempty_str(record["checkpoint_id"], "checkpoint_id", max_len=256)
    checkpoint_sha256 = record["checkpoint_sha256"]
    _need(
        isinstance(checkpoint_sha256, str) and bool(SHA256_RE.fullmatch(checkpoint_sha256)),
        "checkpoint_sha256 must be a 64-character lowercase hexadecimal SHA-256 digest",
    )

    corpus = record["training_corpus"]
    _need(isinstance(corpus, dict), "training_corpus must be a JSON object")
    c_keys = set(corpus.keys())
    _need(c_keys == REQUIRED_CORPUS_KEYS, f"training_corpus keys must be {sorted(REQUIRED_CORPUS_KEYS)}, got {sorted(c_keys)}")
    corpus_name = _require_nonempty_str(corpus["name"], "training_corpus.name")
    corpus_ver = _require_nonempty_str(corpus["version_or_cutoff"], "training_corpus.version_or_cutoff", max_len=128)
    sources = corpus["sources"]
    _need(isinstance(sources, list) and len(sources) > 0, "training_corpus.sources must be a non-empty list")
    norm_sources = [_require_nonempty_str(s, f"training_corpus.sources[{i}]", max_len=128) for i, s in enumerate(sources)]
    tax_scope = _require_nonempty_str(corpus["taxonomic_scope"], "training_corpus.taxonomic_scope", max_len=64)
    _need(
        tax_scope.lower() not in DISALLOWED_TAXONOMIC_SCOPES,
        f"training_corpus.taxonomic_scope {tax_scope!r} is incompatible with bacterial genomes",
    )

    seq_status = record["sequence_overlap_status"]
    pheno_status = record["phenotype_annotation_status"]
    tune_status = record["eval_tuning_status"]
    _need(isinstance(seq_status, str) and seq_status in SEQUENCE_OVERLAP_STATUSES,
          f"invalid sequence_overlap_status: {seq_status!r}")
    _need(isinstance(pheno_status, str) and pheno_status in PHENOTYPE_ANNOTATION_STATUSES,
          f"invalid phenotype_annotation_status: {pheno_status!r}")
    _need(isinstance(tune_status, str) and tune_status in EVAL_TUNING_STATUSES,
          f"invalid eval_tuning_status: {tune_status!r}")

    evidence = record["audit_evidence"]
    _need(isinstance(evidence, dict), "audit_evidence must be a JSON object")
    e_keys = set(evidence.keys())
    _need(e_keys == REQUIRED_EVIDENCE_KEYS, f"audit_evidence keys must be {sorted(REQUIRED_EVIDENCE_KEYS)}, got {sorted(e_keys)}")
    ev_method = _require_nonempty_str(evidence["method"], "audit_evidence.method", max_len=256)
    ev_checked = evidence["checked_at_utc"]
    _need(
        isinstance(ev_checked, str) and (ev_checked == "unverified" or bool(UTC_TS_RE.fullmatch(ev_checked))),
        "audit_evidence.checked_at_utc must be ISO-8601 UTC (YYYY-MM-DDTHH:MM:SSZ) or 'unverified'",
    )
    if seq_status == "verified_disjoint":
        _need(ev_checked != "unverified", "verified_disjoint sequence_overlap_status requires a checked_at_utc timestamp")
    ev_notes = evidence["notes"]
    _need(isinstance(ev_notes, str) and len(ev_notes) <= 1024, "audit_evidence.notes must be a string (<= 1024 chars)")

    derived, reasons = derive_claim_ceiling(seq_status, pheno_status, tune_status)

    if "claimed_ceiling" in record:
        claimed = record["claimed_ceiling"]
        _need(
            isinstance(claimed, str) and claimed in CLAIM_CEILING_RANK,
            f"invalid claimed_ceiling: {claimed!r}",
        )
        _need(
            CLAIM_CEILING_RANK[claimed] <= CLAIM_CEILING_RANK[derived],
            f"claimed_ceiling {claimed!r} exceeds derived claim ceiling {derived!r} ({'; '.join(reasons)})",
        )
        effective_ceiling = claimed
    else:
        effective_ceiling = derived

    return {
        "format": FORMAT,
        "model_family": raw_family,
        "checkpoint_id": checkpoint_id,
        "checkpoint_sha256": checkpoint_sha256,
        "training_corpus": {
            "name": corpus_name,
            "version_or_cutoff": corpus_ver,
            "sources": norm_sources,
            "taxonomic_scope": tax_scope,
        },
        "sequence_overlap_status": seq_status,
        "phenotype_annotation_status": pheno_status,
        "eval_tuning_status": tune_status,
        "audit_evidence": {
            "method": ev_method,
            "checked_at_utc": ev_checked,
            "notes": ev_notes,
        },
        "derived_claim_ceiling": derived,
        "effective_claim_ceiling": effective_ceiling,
        "claim_ceiling_rank": CLAIM_CEILING_RANK[effective_ceiling],
        "claim_ceiling_reasons": reasons,
    }


def summarize_audit_for_report(validated: dict[str, Any]) -> dict[str, Any]:
    """Compact zero-ID audit block embedded inside `grm-rep-compare-v1` reports."""
    return {
        "format": validated["format"],
        "model_family": validated["model_family"],
        "checkpoint_id": validated["checkpoint_id"],
        "checkpoint_sha256": validated["checkpoint_sha256"],
        "training_corpus_name": validated["training_corpus"]["name"],
        "training_corpus_version": validated["training_corpus"]["version_or_cutoff"],
        "sequence_overlap_status": validated["sequence_overlap_status"],
        "phenotype_annotation_status": validated["phenotype_annotation_status"],
        "eval_tuning_status": validated["eval_tuning_status"],
        "derived_claim_ceiling": validated["derived_claim_ceiling"],
        "effective_claim_ceiling": validated["effective_claim_ceiling"],
    }


def min_claim_ceiling(ceilings: list[str]) -> str:
    """Return the strictest (lowest-rank) claim ceiling across combined arms."""
    _need(bool(ceilings), "ceilings must be non-empty")
    for c in ceilings:
        _need(c in CLAIM_CEILING_RANK, f"unknown claim ceiling: {c!r}")
    return min(ceilings, key=lambda c: CLAIM_CEILING_RANK[c])


def load_and_validate(path: str | Path) -> dict[str, Any]:
    """Load an exposure audit JSON file and validate it."""
    p = Path(path)
    try:
        raw = json.loads(p.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise AuditError(f"cannot read exposure audit {p}: {exc}") from exc
    return validate_exposure_audit(raw)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="exposure_audit.py", description="Validate grm-rep-exposure-v1 records.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    v = sub.add_parser("validate", help="Validate a JSON audit record and print normalized summary.")
    v.add_argument("path", help="Path to exposure audit JSON file.")
    args = ap.parse_args(argv)
    try:
        validated = load_and_validate(args.path)
    except AuditError as exc:
        print(f"exposure_audit.py: refused: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(validated, sort_keys=True, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
