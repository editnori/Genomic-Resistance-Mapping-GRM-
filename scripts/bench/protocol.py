"""Two-mode AMR benchmark protocol, 4-tier registry, and strict guards.

Implements:
- Adapter manifest, drug map, and 4-tier species-support registry loaders
  validated directly against src/drugs.rs.
- Canonical genome_calls.jsonl record validation.
- grm-fold-file-v1 validation (matching src/learn/foldfile.rs).
- External deduplicated cohort validation (selection_role: evaluation_only).
- Time-split provenance enforcement (historical_prospective vs
  retrospective_modern_db).
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from scripts.bench.vocab import DRUGS_RS, DrugVocabulary, load_vocabulary

BENCH_DIR = Path(__file__).resolve().parent
ADAPTERS_DIR = BENCH_DIR / "adapters"

ADAPTER_SCHEMA_VERSION = "grm-bench-adapter-v1"
RESULT_FORMAT = "grm-bench-result-v1"
FOLD_FILE_FORMAT = "grm-fold-file-v1"
EXTERNAL_COHORT_FORMAT = "grm-external-cohort-v1"
EVALUATION_REVISION = 14

MODE_PRETRAINED = "pretrained"
MODE_TRAINABLE = "trainable"
VALID_MODES = frozenset({MODE_PRETRAINED, MODE_TRAINABLE})

COMPARISON_END_TO_END = "end_to_end_pipeline"
COMPARISON_FIXED_LEARNER = "fixed_learner_representation"
VALID_COMPARISON_KINDS = frozenset({COMPARISON_END_TO_END, COMPARISON_FIXED_LEARNER})

TIER_SUPPORTED_PRIMARY = "supported_primary"
TIER_SUPPORTED_UNGRADED = "supported_ungraded"
TIER_PROXY_CLASS = "proxy_class"
TIER_UNSUPPORTED = "unsupported"
VALID_SUPPORT_TIERS = frozenset({
    TIER_SUPPORTED_PRIMARY,
    TIER_SUPPORTED_UNGRADED,
    TIER_PROXY_CLASS,
    TIER_UNSUPPORTED,
})

VALID_CALLS = frozenset({"R", "S", "I", "no_call", "failure"})
BINARY_CALLS = frozenset({"R", "S"})
VALID_MAP_KINDS = frozenset({"exact", "combination", "proxy_class"})
VALID_BOOTSTRAP_UNITS = frozenset({"species", "lineage"})

TIME_PROVENANCE_HISTORICAL = "historical_prospective"
TIME_PROVENANCE_RETRO_MODERN_DB = "retrospective_modern_db"
TIME_PROVENANCE_RETRO_SAME_DB = "retrospective_same_db"
TIME_PROVENANCE_EXTERNAL_RETRO = "external_retrospective"
VALID_TIME_PROVENANCES = frozenset({
    TIME_PROVENANCE_HISTORICAL,
    TIME_PROVENANCE_RETRO_MODERN_DB,
    TIME_PROVENANCE_RETRO_SAME_DB,
    TIME_PROVENANCE_EXTERNAL_RETRO,
})

SELECTION_ROLE_EVAL_ONLY = "evaluation_only"
SELECTION_ROLE_INTERNAL_CV = "internal_cv"
VALID_SELECTION_ROLES = frozenset({SELECTION_ROLE_EVAL_ONLY, SELECTION_ROLE_INTERNAL_CV})

PRETRAINED_TOOLS = (
    "amrfinderplus",
    "resfinder_pointfinder",
    "card_rgi",
    "staramr",
    "kleborate",
    "tbprofiler",
    "mykrobe",
)
TRAINABLE_TOOLS = (
    "kover",
    "xgboost",
    "lightgbm",
    "elasticnet",
)
ALL_TOOLS = PRETRAINED_TOOLS + TRAINABLE_TOOLS

DRAWN_SCHEMES = frozenset({"random", "lineage"})
SCHEME_RE = re.compile(r"^[a-z0-9-]{1,48}$")
SPECIES_SLUG_RE = re.compile(r"^[a-z]+_[a-z]+$")
ISO_DATE_RE = re.compile(r"^\d{4}-\d{2}-\d{2}$")

FOLD_FILE_KEYS = frozenset({"format", "dataset", "scheme", "test", "train", "unused"})
GENOME_CALL_REQUIRED_KEYS = frozenset({"genome_id", "species_slug", "reading", "tool_drug", "call"})
GENOME_CALL_ALLOWED_KEYS = GENOME_CALL_REQUIRED_KEYS | frozenset({"score", "reason"})

MANIFEST_REQUIRED_KEYS = frozenset({
    "schema_version",
    "tool",
    "display_name",
    "mode",
    "comparison_kind",
    "input_format",
    "primary_reading",
    "readings",
    "i_handling",
    "genotype_to_drug_rule",
    "software_pin",
    "database_pin",
    "bootstrap_cluster_unit",
    "species_scope",
    "supports_scores",
})

DRUG_MAP_COLUMNS = ["tool_drug", "canonical_drug", "map_kind", "notes"]
SPECIES_SUPPORT_COLUMNS = [
    "species_slug",
    "canonical_drug",
    "support_tier",
    "tool_drugs",
    "unsupported_reason",
    "notes",
]


class BenchError(Exception):
    """Raised when an adapter, input, provenance, or hygiene contract fails."""


def require(cond: bool, message: str) -> None:
    if not cond:
        raise BenchError(message)


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_json_bytes(obj: Any) -> bytes:
    """Deterministic UTF-8 JSON bytes with sorted keys and ASCII escapes."""
    text = json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False) + "\n"
    return text.encode("utf-8")


# ------------------------------------------------------------------ adapter registry


@dataclass(frozen=True)
class DrugMapEntry:
    tool_drug: str
    canonical_drug: str
    map_kind: str
    notes: str


@dataclass(frozen=True)
class SupportEntry:
    species_slug: str
    canonical_drug: str
    support_tier: str
    tool_drugs: tuple[str, ...]
    unsupported_reason: str
    notes: str


@dataclass(frozen=True)
class AdapterSpec:
    tool: str
    adapter_dir: Path
    manifest: dict[str, Any]
    drug_map: tuple[DrugMapEntry, ...]
    species_support: dict[tuple[str, str], SupportEntry]
    manifest_sha256: str
    drug_map_sha256: str
    species_support_sha256: str

    @property
    def mode(self) -> str:
        return self.manifest["mode"]

    @property
    def comparison_kind(self) -> str:
        return self.manifest["comparison_kind"]

    @property
    def primary_reading(self) -> str:
        return self.manifest["primary_reading"]

    @property
    def bootstrap_cluster_unit(self) -> str:
        return self.manifest["bootstrap_cluster_unit"]

    def support_for(self, species_slug: str, canonical_drug: str) -> SupportEntry:
        """Return pre-committed support entry; unknown pairs are strictly unsupported."""
        entry = self.species_support.get((species_slug, canonical_drug))
        if entry is not None:
            return entry
        return SupportEntry(
            species_slug=species_slug,
            canonical_drug=canonical_drug,
            support_tier=TIER_UNSUPPORTED,
            tool_drugs=(),
            unsupported_reason="not_in_precommitted_support_registry",
            notes="Pair absent from pre-committed species_support.tsv; never mapped to S.",
        )

    def tool_drugs_for_canonical(self, canonical_drug: str) -> tuple[str, ...]:
        return tuple(sorted({e.tool_drug for e in self.drug_map if e.canonical_drug == canonical_drug}))


def validate_manifest(manifest: dict[str, Any], expected_tool: str | None = None) -> None:
    require(isinstance(manifest, dict), "manifest must be a JSON object")
    missing = MANIFEST_REQUIRED_KEYS - set(manifest)
    require(not missing, f"manifest missing required keys: {sorted(missing)}")
    require(
        manifest["schema_version"] == ADAPTER_SCHEMA_VERSION,
        f"manifest schema_version {manifest['schema_version']!r} != {ADAPTER_SCHEMA_VERSION!r}",
    )
    tool = manifest["tool"]
    require(isinstance(tool, str) and bool(re.fullmatch(r"[a-z0-9_]+", tool)), f"invalid tool name: {tool!r}")
    if expected_tool is not None:
        require(tool == expected_tool, f"manifest tool {tool!r} != directory {expected_tool!r}")
    require(manifest["mode"] in VALID_MODES, f"{tool}: invalid mode {manifest['mode']!r}")
    require(
        manifest["comparison_kind"] in VALID_COMPARISON_KINDS,
        f"{tool}: invalid comparison_kind {manifest['comparison_kind']!r}",
    )
    require(
        manifest["bootstrap_cluster_unit"] in VALID_BOOTSTRAP_UNITS,
        f"{tool}: invalid bootstrap_cluster_unit {manifest['bootstrap_cluster_unit']!r}",
    )
    readings = manifest["readings"]
    require(
        isinstance(readings, list) and len(readings) >= 1 and all(isinstance(r, str) and r for r in readings),
        f"{tool}: readings must be a non-empty list of strings",
    )
    require(
        manifest["primary_reading"] in readings,
        f"{tool}: primary_reading {manifest['primary_reading']!r} not in readings",
    )
    i_hand = manifest["i_handling"]
    require(isinstance(i_hand, dict), f"{tool}: i_handling must be an object")
    require(
        isinstance(i_hand.get("tool_emits_intermediate"), bool),
        f"{tool}: i_handling.tool_emits_intermediate must be bool",
    )
    policy = i_hand.get("policy", "")
    require(
        isinstance(policy, str) and policy == "separate_count_never_map_to_s",
        f"{tool}: i_handling.policy must be 'separate_count_never_map_to_s', got {policy!r}",
    )
    require(
        isinstance(i_hand.get("description"), str) and len(i_hand["description"].strip()) >= 10,
        f"{tool}: i_handling.description must be a non-empty explanation",
    )
    rule = manifest["genotype_to_drug_rule"]
    require(
        isinstance(rule, str) and len(rule.strip()) >= 15,
        f"{tool}: genotype_to_drug_rule must be explicitly declared",
    )
    sw = manifest["software_pin"]
    require(
        isinstance(sw, dict)
        and isinstance(sw.get("package"), str)
        and bool(sw["package"])
        and isinstance(sw.get("version"), str)
        and bool(sw["version"]),
        f"{tool}: software_pin must specify package and version",
    )
    if sw.get("release_date") is not None:
        require(
            isinstance(sw["release_date"], str) and bool(ISO_DATE_RE.fullmatch(sw["release_date"])),
            f"{tool}: software_pin.release_date must be YYYY-MM-DD",
        )
    db = manifest["database_pin"]
    if manifest["mode"] == MODE_PRETRAINED:
        require(
            isinstance(db, dict)
            and isinstance(db.get("name"), str)
            and bool(db["name"])
            and isinstance(db.get("version"), str)
            and bool(db["version"]),
            f"{tool}: pretrained adapter database_pin must specify name and version",
        )
        if db.get("release_date") is not None:
            require(
                isinstance(db["release_date"], str) and bool(ISO_DATE_RE.fullmatch(db["release_date"])),
                f"{tool}: database_pin.release_date must be YYYY-MM-DD",
            )
    else:
        require(db is None or isinstance(db, dict), f"{tool}: database_pin must be null or object")
    require(isinstance(manifest["supports_scores"], bool), f"{tool}: supports_scores must be bool")


def _read_tsv(path: Path, expected_columns: list[str]) -> list[dict[str, str]]:
    require(path.is_file(), f"missing required TSV: {path}")
    raw = path.read_text(encoding="utf-8")
    require(raw.endswith("\n"), f"{path}: missing trailing newline")
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        require(header == expected_columns, f"{path}: header {header} != {expected_columns}")
        rows: list[dict[str, str]] = []
        for lineno, cells in enumerate(reader, 2):
            require(
                len(cells) == len(expected_columns),
                f"{path}:{lineno}: expected {len(expected_columns)} columns, got {len(cells)}",
            )
            rows.append(dict(zip(expected_columns, cells)))
    return rows


def load_adapter(
    tool: str,
    adapters_dir: Path = ADAPTERS_DIR,
    vocab: DrugVocabulary | None = None,
) -> AdapterSpec:
    """Load and strictly validate one adapter's manifest, drug_map, and species_support."""
    v = vocab if vocab is not None else load_vocabulary(DRUGS_RS)
    adir = Path(adapters_dir) / tool
    require(adir.is_dir(), f"adapter directory not found: {adir}")

    manifest_path = adir / "manifest.json"
    drug_map_path = adir / "drug_map.tsv"
    support_path = adir / "species_support.tsv"

    require(manifest_path.is_file(), f"{tool}: missing manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validate_manifest(manifest, expected_tool=tool)

    dm_rows = _read_tsv(drug_map_path, DRUG_MAP_COLUMNS)
    require(len(dm_rows) > 0, f"{tool}: drug_map.tsv is empty")
    drug_entries: list[DrugMapEntry] = []
    seen_dm: set[tuple[str, str]] = set()
    dm_by_canonical: dict[str, set[str]] = {}
    for row in dm_rows:
        td = row["tool_drug"].strip()
        cd = row["canonical_drug"].strip()
        mk = row["map_kind"].strip()
        require(bool(td), f"{tool}: empty tool_drug in drug_map.tsv")
        require(
            v.is_canonical(cd),
            f"{tool}: drug_map.tsv canonical_drug {cd!r} is not valid in src/drugs.rs",
        )
        require(mk in VALID_MAP_KINDS, f"{tool}: invalid map_kind {mk!r} for {td} -> {cd}")
        key = (td, cd)
        require(key not in seen_dm, f"{tool}: duplicate (tool_drug, canonical_drug) {key} in drug_map.tsv")
        seen_dm.add(key)
        dm_by_canonical.setdefault(cd, set()).add(td)
        drug_entries.append(DrugMapEntry(tool_drug=td, canonical_drug=cd, map_kind=mk, notes=row["notes"].strip()))

    sp_rows = _read_tsv(support_path, SPECIES_SUPPORT_COLUMNS)
    require(len(sp_rows) > 0, f"{tool}: species_support.tsv is empty")
    support_map: dict[tuple[str, str], SupportEntry] = {}
    for row in sp_rows:
        sp = row["species_slug"].strip()
        cd = row["canonical_drug"].strip()
        tier = row["support_tier"].strip()
        td_raw = row["tool_drugs"].strip()
        reason = row["unsupported_reason"].strip()
        notes = row["notes"].strip()
        require(bool(SPECIES_SLUG_RE.fullmatch(sp)), f"{tool}: invalid species_slug {sp!r}")
        require(
            v.is_canonical(cd),
            f"{tool}: species_support.tsv canonical_drug {cd!r} is not valid in src/drugs.rs",
        )
        require(tier in VALID_SUPPORT_TIERS, f"{tool}: invalid support_tier {tier!r} for ({sp}, {cd})")
        key = (sp, cd)
        require(key not in support_map, f"{tool}: duplicate ({sp}, {cd}) in species_support.tsv")
        if tier == TIER_UNSUPPORTED:
            require(td_raw == "-", f"{tool}: unsupported ({sp}, {cd}) must have tool_drugs='-'")
            require(
                reason != "" and reason != "-",
                f"{tool}: unsupported ({sp}, {cd}) must declare a non-empty unsupported_reason",
            )
            tds: tuple[str, ...] = ()
        else:
            require(
                reason == "-",
                f"{tool}: supported/proxy ({sp}, {cd}) must have unsupported_reason='-'",
            )
            tds = tuple(t.strip() for t in td_raw.split(",") if t.strip())
            require(len(tds) >= 1 and td_raw != "-", f"{tool}: ({sp}, {cd}) missing tool_drugs")
            mapped_tds = dm_by_canonical.get(cd, set())
            for td in tds:
                require(
                    td in mapped_tds,
                    f"{tool}: ({sp}, {cd}) references tool_drug {td!r} not mapped to {cd!r} in drug_map.tsv",
                )
        support_map[key] = SupportEntry(
            species_slug=sp,
            canonical_drug=cd,
            support_tier=tier,
            tool_drugs=tds,
            unsupported_reason=reason,
            notes=notes,
        )

    return AdapterSpec(
        tool=tool,
        adapter_dir=adir,
        manifest=manifest,
        drug_map=tuple(drug_entries),
        species_support=support_map,
        manifest_sha256=sha256_file(manifest_path),
        drug_map_sha256=sha256_file(drug_map_path),
        species_support_sha256=sha256_file(support_path),
    )


def load_all_adapters(
    adapters_dir: Path = ADAPTERS_DIR,
    vocab: DrugVocabulary | None = None,
) -> dict[str, AdapterSpec]:
    v = vocab if vocab is not None else load_vocabulary(DRUGS_RS)
    return {tool: load_adapter(tool, adapters_dir=adapters_dir, vocab=v) for tool in ALL_TOOLS}


# ------------------------------------------------------------------ genome_calls.jsonl validation


def validate_genome_call_record(record: dict[str, Any]) -> dict[str, Any]:
    """Validate one canonical genome_calls.jsonl record."""
    require(isinstance(record, dict), "genome_calls record must be a JSON object")
    keys = set(record)
    missing = GENOME_CALL_REQUIRED_KEYS - keys
    extra = keys - GENOME_CALL_ALLOWED_KEYS
    require(not missing, f"genome_calls record missing keys: {sorted(missing)}")
    require(not extra, f"genome_calls record has unlisted keys: {sorted(extra)}")
    for k in ("genome_id", "species_slug", "reading", "tool_drug"):
        require(isinstance(record[k], str) and bool(record[k].strip()), f"genome_calls.{k} must be non-empty str")
    require(bool(SPECIES_SLUG_RE.fullmatch(record["species_slug"])), f"invalid species_slug: {record['species_slug']!r}")
    call = record["call"]
    require(call in VALID_CALLS, f"invalid call {call!r}; must be one of {sorted(VALID_CALLS)}")
    if "score" in record and record["score"] is not None:
        sc = record["score"]
        require(
            isinstance(sc, (int, float)) and not isinstance(sc, bool) and math.isfinite(float(sc)) and 0.0 <= float(sc) <= 1.0,
            f"genome_calls.score must be finite float in [0, 1], got {sc!r}",
        )
    if "reason" in record and record["reason"] is not None:
        require(isinstance(record["reason"], str), "genome_calls.reason must be str or null")
    if call in ("no_call", "failure", "I"):
        require(
            isinstance(record.get("reason"), str) and bool(record["reason"].strip()),
            f"genome_calls record with call={call!r} must provide a non-empty reason",
        )
    return record


# ------------------------------------------------------------------ fold-file and external cohort


@dataclass(frozen=True)
class LoadedFoldFile:
    dataset: str
    scheme: str
    test: tuple[tuple[str, ...], ...]
    train: tuple[str, ...]
    unused: tuple[str, ...]
    sha256: str


def load_fold_file(path: Path | str, expected_dataset: str | None = None, sample_ids: list[str] | None = None) -> LoadedFoldFile:
    """Parse and strictly validate a grm-fold-file-v1 file (mirroring src/learn/foldfile.rs)."""
    p = Path(path)
    raw = p.read_bytes()
    digest = sha256_bytes(raw)
    try:
        doc = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise BenchError(f"fold file {p}: parse error: {exc}") from exc
    require(isinstance(doc, dict), f"fold file {p}: root must be an object")
    extra = set(doc) - FOLD_FILE_KEYS
    missing = FOLD_FILE_KEYS - set(doc)
    require(not extra, f"fold file {p}: unknown fields {sorted(extra)}")
    require(not missing, f"fold file {p}: missing fields {sorted(missing)}")
    require(doc["format"] == FOLD_FILE_FORMAT, f"fold file format {doc['format']!r} != {FOLD_FILE_FORMAT!r}")
    scheme = doc["scheme"]
    require(
        isinstance(scheme, str) and bool(SCHEME_RE.fullmatch(scheme)),
        f"fold file scheme {scheme!r} is not 1 to 48 lowercase letters, digits and hyphens",
    )
    require(
        scheme not in DRAWN_SCHEMES,
        f"fold file scheme {scheme!r} conflicts with drawn split name",
    )
    test_folds = doc["test"]
    require(isinstance(test_folds, list) and len(test_folds) >= 1, "fold file has no test fold")
    for idx, fold_ids in enumerate(test_folds):
        require(
            isinstance(fold_ids, list) and len(fold_ids) >= 1 and all(isinstance(g, str) and g for g in fold_ids),
            f"fold file test fold {idx} is empty or invalid",
        )
    for key in ("train", "unused"):
        require(
            isinstance(doc[key], list) and all(isinstance(g, str) and g for g in doc[key]),
            f"fold file {key} must be a list of genome ids",
        )
    dataset = doc["dataset"]
    require(isinstance(dataset, str) and bool(dataset), "fold file dataset must be non-empty")
    if expected_dataset is not None:
        require(
            dataset == expected_dataset,
            f"fold file is for dataset {dataset!r}, not {expected_dataset!r}",
        )
    all_listed: list[str] = [g for fold_ids in test_folds for g in fold_ids] + list(doc["train"]) + list(doc["unused"])
    require(len(set(all_listed)) == len(all_listed), "fold file lists a genome twice")
    if sample_ids is not None:
        require(len(set(sample_ids)) == len(sample_ids), f"{dataset}: duplicate sample id in cohort")
        cohort_set = set(sample_ids)
        listed_set = set(all_listed)
        foreign = listed_set - cohort_set
        require(not foreign, f"fold file names {len(foreign)} genome(s) outside cohort {dataset}")
        absent = cohort_set - listed_set
        require(not absent, f"fold file leaves out {len(absent)} sample(s) of cohort {dataset}")
    return LoadedFoldFile(
        dataset=dataset,
        scheme=scheme,
        test=tuple(tuple(f) for f in test_folds),
        train=tuple(doc["train"]),
        unused=tuple(doc["unused"]),
        sha256=digest,
    )


def validate_external_cohort(cohort: dict[str, Any], *, allow_selection: bool = False) -> dict[str, Any]:
    """Validate an external deduplicated cohort descriptor.

    External validation cohorts MUST declare selection_role == 'evaluation_only',
    and any attempt to use them for tuning or model/tool selection is refused.
    """
    require(isinstance(cohort, dict), "external cohort must be a JSON object")
    require(
        cohort.get("format") == EXTERNAL_COHORT_FORMAT,
        f"external cohort format {cohort.get('format')!r} != {EXTERNAL_COHORT_FORMAT!r}",
    )
    role = cohort.get("selection_role")
    require(
        role == SELECTION_ROLE_EVAL_ONLY,
        f"external validation cohort must declare selection_role={SELECTION_ROLE_EVAL_ONLY!r}, got {role!r}",
    )
    require(
        not allow_selection,
        "external validation cohorts have selection_role='evaluation_only' and cannot be used for tuning or model/tool selection",
    )
    require(isinstance(cohort.get("cohort_name"), str) and bool(cohort["cohort_name"].strip()), "missing cohort_name")
    require(isinstance(cohort.get("dataset"), str) and bool(cohort["dataset"].strip()), "missing dataset")
    sp = cohort.get("species_slug")
    require(isinstance(sp, str) and bool(SPECIES_SLUG_RE.fullmatch(sp)), f"invalid species_slug {sp!r}")
    ab = cohort.get("antibiotic")
    vocab = load_vocabulary(DRUGS_RS)
    require(isinstance(ab, str) and vocab.is_canonical(ab), f"external cohort antibiotic {ab!r} not in src/drugs.rs")
    require(cohort.get("deduplicated") is True, "external cohort must declare deduplicated=true")
    samples = cohort.get("samples")
    resistant = cohort.get("resistant")
    require(
        isinstance(samples, list) and len(samples) >= 1 and all(isinstance(g, str) and g for g in samples),
        "external cohort samples must be a non-empty list of genome ids",
    )
    require(len(set(samples)) == len(samples), "external cohort contains duplicate genome ids")
    require(
        isinstance(resistant, list) and len(resistant) == len(samples) and all(isinstance(y, bool) for y in resistant),
        "external cohort resistant labels must be one bool per sample",
    )
    if "lineage" in cohort and cohort["lineage"] is not None:
        lin = cohort["lineage"]
        require(
            isinstance(lin, list) and len(lin) == len(samples) and all(isinstance(x, (int, str)) for x in lin),
            "external cohort lineage must align 1:1 with samples",
        )
    return cohort


def assert_selection_allowed(selection_role: str, purpose: str) -> None:
    """Refuse any model/tool/threshold selection or tuning on an evaluation_only cohort."""
    require(selection_role in VALID_SELECTION_ROLES, f"unknown selection_role {selection_role!r}")
    if selection_role == SELECTION_ROLE_EVAL_ONLY:
        raise BenchError(
            f"refused: {purpose} is forbidden on selection_role={SELECTION_ROLE_EVAL_ONLY!r} cohorts"
        )


# ------------------------------------------------------------------ time-split provenance


def validate_time_provenance(
    *,
    split_scheme: str,
    declared_label: str,
    fold_cutoff_date: str | None = None,
    software_release_date: str | None = None,
    database_release_date: str | None = None,
    weights_release_date: str | None = None,
    training_labels_max_date: str | None = None,
    mode: str = MODE_PRETRAINED,
    is_bvbrc_insertion_time_split: bool = False,
) -> dict[str, Any]:
    """Enforce strict time-split provenance labeling.

    A `historical_prospective` label is accepted ONLY if:
    - `fold_cutoff_date` is a valid ISO date YYYY-MM-DD;
    - `software_release_date` is valid and <= `fold_cutoff_date`;
    - `training_labels_max_date` is valid and <= `fold_cutoff_date`;
    - For pretrained tools, `database_release_date` is valid and <= `fold_cutoff_date`
      (and `weights_release_date`, if present, is <= `fold_cutoff_date`);
    - For trainable tools, any supplied `database_release_date` or `weights_release_date`
      is <= `fold_cutoff_date`;
    - If `is_bvbrc_insertion_time_split` is True and a modern database (dated after
      `fold_cutoff_date` or with undated DB) is used, the label MUST be
      `retrospective_modern_db`.
    """
    require(declared_label in VALID_TIME_PROVENANCES, f"invalid time_provenance label {declared_label!r}")

    def _check_date(name: str, val: str | None) -> str | None:
        if val is None:
            return None
        require(isinstance(val, str) and bool(ISO_DATE_RE.fullmatch(val)), f"{name} must be YYYY-MM-DD, got {val!r}")
        return val

    cutoff = _check_date("fold_cutoff_date", fold_cutoff_date)
    sw_date = _check_date("software_release_date", software_release_date)
    db_date = _check_date("database_release_date", database_release_date)
    wt_date = _check_date("weights_release_date", weights_release_date)
    lbl_date = _check_date("training_labels_max_date", training_labels_max_date)

    modern_db_detected = False
    if cutoff is not None and mode == MODE_PRETRAINED:
        if db_date is None or db_date > cutoff:
            modern_db_detected = True
        if sw_date is None or sw_date > cutoff:
            modern_db_detected = True
        if wt_date is not None and wt_date > cutoff:
            modern_db_detected = True

    if is_bvbrc_insertion_time_split and mode == MODE_PRETRAINED and modern_db_detected:
        require(
            declared_label == TIME_PROVENANCE_RETRO_MODERN_DB,
            f"BV-BRC insertion-date time split with modern or post-cutoff database/software "
            f"(cutoff={cutoff}, sw={sw_date}, db={db_date}, weights={wt_date}) MUST be labeled "
            f"{TIME_PROVENANCE_RETRO_MODERN_DB!r}, not {declared_label!r}",
        )

    if declared_label == TIME_PROVENANCE_HISTORICAL:
        require(cutoff is not None, "historical_prospective requires fold_cutoff_date")
        require(sw_date is not None and sw_date <= cutoff, f"historical_prospective requires software_release_date ({sw_date}) <= fold_cutoff_date ({cutoff})")
        require(
            lbl_date is not None and lbl_date <= cutoff,
            f"historical_prospective requires training_labels_max_date ({lbl_date}) <= fold_cutoff_date ({cutoff})",
        )
        if mode == MODE_PRETRAINED:
            require(
                db_date is not None and db_date <= cutoff,
                f"historical_prospective pretrained run requires database_release_date ({db_date}) <= fold_cutoff_date ({cutoff})",
            )
            if wt_date is not None:
                require(
                    wt_date <= cutoff,
                    f"historical_prospective requires weights_release_date ({wt_date}) <= fold_cutoff_date ({cutoff})",
                )
        else:
            if db_date is not None:
                require(db_date <= cutoff, f"historical_prospective requires database_release_date ({db_date}) <= {cutoff}")
            if wt_date is not None:
                require(wt_date <= cutoff, f"historical_prospective requires weights_release_date ({wt_date}) <= {cutoff}")

    return {
        "split_scheme": split_scheme,
        "label": declared_label,
        "fold_cutoff_date": cutoff,
        "software_release_date": sw_date,
        "database_release_date": db_date,
        "weights_release_date": wt_date,
        "training_labels_max_date": lbl_date,
        "is_bvbrc_insertion_time_split": bool(is_bvbrc_insertion_time_split),
    }
