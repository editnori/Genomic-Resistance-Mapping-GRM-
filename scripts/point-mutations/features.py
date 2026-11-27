"""Fold-leak-free point-mutation and foundation-model feature builders (Track B).

Every encoder fits its vocabulary, summary statistics, scalers, and projections
strictly on training-fold rows and applies them unmodified to validation/test
rows. No outer test row ever enters vocabulary construction, coordinate/pocket
selection, imputation, standardization, or PCA fitting.

Builders provided:
- `ProteinVariantEncoder`: protein-aligned known-target variant presence matrix
  built strictly from training isolates and/or a frozen published target
  catalogue, with an explicit `unseen_variant_count` side-feature.
- `BurdenEncoder`: per-gene and per-predefined-region synonymous,
  nonsynonymous, truncating, and indel burden counts from frozen public region
  coordinates, with optional training-fold-only standardization.
- `FoundationScoreAdapter`: precomputed per-isolate zero-shot LLRs or
  embeddings from ESM-2, Evo 2, ProkBERT, or Bacformer aligned strictly to
  `grm-dataset-export-v1` or `grm-fold-file-v1` row manifests, gated on a
  valid `grm-rep-exposure-v1` audit record, with training-fold-only
  standardization and deterministic sign-canonical PCA projection.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np


def _load_exposure_audit():
    key = "grm_point_mutations.exposure_audit"
    if key in sys.modules:
        return sys.modules[key]
    path = Path(__file__).resolve().parent / "exposure_audit.py"
    spec = importlib.util.spec_from_file_location(key, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


EA = _load_exposure_audit()
AuditError = EA.AuditError

FOUNDATION_MANIFEST_FORMAT = "grm-rep-foundation-v1"
ALLOWED_FOUNDATION_KINDS = frozenset({"zero_shot_llr", "embedding", "combined"})
ALLOWED_VOCAB_MODES = frozenset({"train_only", "catalogue_only", "train_union_catalogue"})

BURDEN_CONSEQUENCES = ("synonymous", "nonsynonymous", "truncating", "indel")

KIND_TO_BURDEN = {
    "synonymous": "synonymous",
    "missense": "nonsynonymous",
    "nonsynonymous": "nonsynonymous",
    "substitution": "nonsynonymous",
    "nonsense": "truncating",
    "stop_gained": "truncating",
    "frameshift": "truncating",
    "start_lost": "truncating",
    "truncating": "truncating",
    "indel": "indel",
    "inframe_indel": "indel",
    "inframe_insertion": "indel",
    "inframe_deletion": "indel",
}


class FeatureError(ValueError):
    """Raised on invalid representation input, misalignment, or misuse."""


def _need(cond: bool, message: str) -> None:
    if not cond:
        raise FeatureError(message)


def _is_strict_int(v: Any) -> bool:
    return isinstance(v, int) and not isinstance(v, bool)


def _canonical_json_bytes(obj: Any) -> bytes:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")


def _parse_variant_entry(entry: Any) -> tuple[str, str, int | None, int | None, str | None]:
    """Normalize one variant dict into (gene, change, pos_start, pos_end, kind)."""
    _need(isinstance(entry, dict), "variant entry must be a dict")
    gene = entry.get("gene")
    change = entry.get("change")
    _need(isinstance(gene, str) and bool(gene.strip()), "variant entry requires non-empty string 'gene'")
    _need(isinstance(change, str) and bool(change.strip()), "variant entry requires non-empty string 'change'")
    gene = gene.strip()
    change = change.strip()

    pos = entry.get("pos")
    pos_start = entry.get("pos_start", pos)
    pos_end = entry.get("pos_end", pos_start)
    if pos_start is not None:
        _need(_is_strict_int(pos_start) and pos_start >= 1, "variant pos/pos_start must be an integer >= 1")
    if pos_end is not None:
        _need(_is_strict_int(pos_end) and pos_end >= 1, "variant pos_end must be an integer >= 1")
        if pos_start is not None:
            _need(pos_end >= pos_start, "variant pos_end must be >= pos_start")

    kind = entry.get("kind")
    if kind is not None:
        _need(isinstance(kind, str) and bool(kind.strip()), "variant kind must be a non-empty string or None")
        kind = kind.strip().lower()
    return gene, change, pos_start, pos_end, kind


def _fit_standardizer(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compute column mean and inverse population std strictly on training rows.

    Columns with zero variance on the training fold get scale = 0.0 so
    `(X - mean) * scale` maps constant columns cleanly to 0.0.
    """
    _need(X.ndim == 2 and X.shape[0] > 0, "cannot fit standardizer on empty matrix")
    mean = np.mean(X, axis=0, dtype=np.float64)
    centered = X - mean
    sd = np.sqrt(np.mean(centered * centered, axis=0, dtype=np.float64))
    scale = np.divide(1.0, sd, out=np.zeros_like(sd), where=sd > 0.0)
    return mean, scale


def _deterministic_pca_components(X_scaled: np.ndarray, n_components: int) -> np.ndarray:
    """Fit deterministic sign-canonical PCA components (k x d) on centered/scaled train rows."""
    _need(_is_strict_int(n_components) and n_components >= 1, "pca_components must be an integer >= 1")
    n_samples, n_features = X_scaled.shape
    k = min(n_components, n_samples, n_features)
    if k == 0:
        return np.zeros((0, n_features), dtype=np.float64)
    # Center explicitly on training rows (even when unstandardized).
    Xc = X_scaled - np.mean(X_scaled, axis=0, dtype=np.float64)
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    comps = np.ascontiguousarray(vt[:k, :], dtype=np.float64)
    # Canonical sign per component: largest absolute coordinate is strictly positive,
    # ties broken deterministically by lowest column index (np.argmax behavior).
    for i in range(k):
        row = comps[i]
        idx = int(np.argmax(np.abs(row)))
        if row[idx] < 0.0:
            comps[i] = -row
        # Normalize -0.0 to +0.0 for bit-for-bit determinism.
        comps[i] = np.where(comps[i] == 0.0, 0.0, comps[i])
    return comps


class ProteinVariantEncoder:
    """Protein-aligned known-target variant encoder with explicit unseen-variant side-feature.

    Vocabulary is built strictly inside the training fold (`train_only`), from a
    frozen published target catalogue (`catalogue_only`), or their union
    (`train_union_catalogue`). Test isolates are encoded strictly through the
    frozen fitted vocabulary; any distinct eligible variant outside the fitted
    vocabulary is counted in the explicit `unseen_variant_count` column.
    """

    def __init__(
        self,
        target_genes: Iterable[str] | None = None,
        frozen_catalogue_variants: Iterable[tuple[str, str]] | None = None,
        vocabulary_mode: str = "train_only",
        min_support: int = 1,
        max_support_margin: int = 0,
        include_synonymous: bool = False,
        include_unseen_side_feature: bool = True,
    ) -> None:
        _need(vocabulary_mode in ALLOWED_VOCAB_MODES, f"invalid vocabulary_mode: {vocabulary_mode!r}")
        _need(_is_strict_int(min_support) and min_support >= 1, "min_support must be an integer >= 1")
        _need(_is_strict_int(max_support_margin) and max_support_margin >= 0, "max_support_margin must be an integer >= 0")
        self.target_genes: frozenset[str] | None = (
            frozenset(g.strip() for g in target_genes) if target_genes is not None else None
        )
        if self.target_genes is not None:
            _need(all(bool(g) for g in self.target_genes), "target_genes cannot contain blank names")
        self.frozen_catalogue_variants: tuple[tuple[str, str], ...] = tuple(
            sorted({(g.strip(), c.strip()) for g, c in (frozen_catalogue_variants or ())})
        )
        for g, c in self.frozen_catalogue_variants:
            _need(bool(g) and bool(c), "frozen_catalogue_variants cannot contain blank gene/change")
            if self.target_genes is not None:
                _need(g in self.target_genes, f"catalogue variant gene {g!r} is outside target_genes")
        self.vocabulary_mode = vocabulary_mode
        self.min_support = min_support
        self.max_support_margin = max_support_margin
        self.include_synonymous = bool(include_synonymous)
        self.include_unseen_side_feature = bool(include_unseen_side_feature)

        self.fitted_: bool = False
        self.vocabulary_: tuple[tuple[str, str], ...] = ()
        self._vocab_index: dict[tuple[str, str], int] = {}
        self.feature_names_: tuple[str, ...] = ()

    def _eligible_variants(self, isolate_variants: Sequence[dict[str, Any]]) -> set[tuple[str, str]]:
        out: set[tuple[str, str]] = set()
        for entry in isolate_variants:
            gene, change, _, _, kind = _parse_variant_entry(entry)
            if self.target_genes is not None and gene not in self.target_genes:
                continue
            if not self.include_synonymous and kind == "synonymous":
                continue
            out.add((gene, change))
        return out

    def fit(self, train_isolates: Sequence[Sequence[dict[str, Any]]]) -> "ProteinVariantEncoder":
        n_train = len(train_isolates)
        _need(n_train > 0, "ProteinVariantEncoder.fit requires at least one training isolate")
        vocab_set: set[tuple[str, str]] = set()
        if self.vocabulary_mode in ("catalogue_only", "train_union_catalogue"):
            vocab_set.update(self.frozen_catalogue_variants)
        if self.vocabulary_mode in ("train_only", "train_union_catalogue"):
            counts: dict[tuple[str, str], int] = {}
            for iso in train_isolates:
                for key in self._eligible_variants(iso):
                    counts[key] = counts.get(key, 0) + 1
            hi = n_train - self.max_support_margin
            for key, cnt in counts.items():
                if cnt >= self.min_support and cnt <= hi:
                    vocab_set.add(key)

        self.vocabulary_ = tuple(sorted(vocab_set))
        self._vocab_index = {key: i for i, key in enumerate(self.vocabulary_)}
        names = [f"{gene}:{change}" for gene, change in self.vocabulary_]
        if self.include_unseen_side_feature:
            names.append("unseen_variant_count")
        self.feature_names_ = tuple(names)
        self.fitted_ = True
        return self

    def transform(self, isolates: Sequence[Sequence[dict[str, Any]]]) -> np.ndarray:
        _need(self.fitted_, "ProteinVariantEncoder must be fitted before transform")
        n = len(isolates)
        n_vocab = len(self.vocabulary_)
        n_cols = n_vocab + (1 if self.include_unseen_side_feature else 0)
        X = np.zeros((n, n_cols), dtype=np.float64)
        for row_i, iso in enumerate(isolates):
            unseen = 0
            for key in self._eligible_variants(iso):
                col = self._vocab_index.get(key)
                if col is not None:
                    X[row_i, col] = 1.0
                else:
                    unseen += 1
            if self.include_unseen_side_feature:
                X[row_i, n_vocab] = float(unseen)
        return X

    def fit_transform(self, train_isolates: Sequence[Sequence[dict[str, Any]]]) -> np.ndarray:
        return self.fit(train_isolates).transform(train_isolates)

    def state_fingerprint(self) -> str:
        _need(self.fitted_, "ProteinVariantEncoder must be fitted before fingerprinting")
        payload = {
            "vocabulary": list(self.vocabulary_),
            "feature_names": list(self.feature_names_),
        }
        return hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()


class BurdenEncoder:
    """Gene and predefined-region variant burden summarizer.

    Aggregates per-gene and per-predefined-region counts across canonical
    consequence classes (`synonymous`, `nonsynonymous`, `truncating`, `indel`).
    Target genes and active-site/region coordinates are frozen prior to fitting;
    standardization parameters are fit strictly on training rows.
    """

    def __init__(
        self,
        target_genes: Iterable[str],
        regions: Sequence[dict[str, Any]] | None = None,
        standardize: bool = True,
    ) -> None:
        genes = sorted({g.strip() for g in target_genes})
        _need(bool(genes) and all(bool(g) for g in genes), "BurdenEncoder requires non-empty target_genes")
        self.target_genes: tuple[str, ...] = tuple(genes)
        self._gene_set: frozenset[str] = frozenset(genes)

        parsed_regions: list[tuple[str, str, int, int]] = []
        seen_reg: set[tuple[str, str]] = set()
        for reg in regions or ():
            _need(isinstance(reg, dict), "each region spec must be a dict")
            _need(
                set(reg.keys()) == {"gene", "region", "start", "end"},
                f"region spec must have keys gene, region, start, end; got {sorted(reg.keys())}",
            )
            gene = reg["gene"]
            rname = reg["region"]
            start = reg["start"]
            end = reg["end"]
            _need(isinstance(gene, str) and gene.strip() in self._gene_set,
                  f"region gene {gene!r} is not in target_genes")
            _need(isinstance(rname, str) and bool(rname.strip()), "region name must be non-empty string")
            _need(_is_strict_int(start) and start >= 1, "region start must be integer >= 1")
            _need(_is_strict_int(end) and end >= start, "region end must be integer >= start")
            key = (gene.strip(), rname.strip())
            _need(key not in seen_reg, f"duplicate region specification: {key}")
            seen_reg.add(key)
            parsed_regions.append((key[0], key[1], start, end))
        self.regions: tuple[tuple[str, str, int, int], ...] = tuple(sorted(parsed_regions))
        self.standardize = bool(standardize)

        names: list[str] = []
        for g in self.target_genes:
            for cons in BURDEN_CONSEQUENCES:
                names.append(f"gene:{g}:{cons}")
        for g, rname, _, _ in self.regions:
            for cons in BURDEN_CONSEQUENCES:
                names.append(f"region:{g}:{rname}:{cons}")
        self.feature_names_: tuple[str, ...] = tuple(names)
        self._gene_index = {g: i for i, g in enumerate(self.target_genes)}
        self._regions_by_gene: dict[str, list[tuple[int, int, int]]] = {}
        for r_idx, (g, _, start, end) in enumerate(self.regions):
            self._regions_by_gene.setdefault(g, []).append((r_idx, start, end))

        self.fitted_: bool = False
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None

    def _raw_counts(self, isolates: Sequence[Sequence[dict[str, Any]]]) -> np.ndarray:
        n = len(isolates)
        n_genes = len(self.target_genes)
        n_cons = len(BURDEN_CONSEQUENCES)
        cons_idx = {c: i for i, c in enumerate(BURDEN_CONSEQUENCES)}
        X = np.zeros((n, len(self.feature_names_)), dtype=np.float64)
        region_base = n_genes * n_cons
        for row_i, iso in enumerate(isolates):
            seen_gene_variants: set[tuple[str, str, str]] = set()
            seen_reg_variants: set[tuple[int, str, str]] = set()
            for entry in iso:
                gene, change, pos_start, pos_end, kind = _parse_variant_entry(entry)
                if gene not in self._gene_set:
                    continue
                _need(kind is not None and kind in KIND_TO_BURDEN,
                      f"BurdenEncoder requires known variant kind, got {kind!r} on {gene}:{change}")
                bucket = KIND_TO_BURDEN[kind]
                ci = cons_idx[bucket]
                g_key = (gene, change, bucket)
                if g_key not in seen_gene_variants:
                    seen_gene_variants.add(g_key)
                    gi = self._gene_index[gene]
                    X[row_i, gi * n_cons + ci] += 1.0
                if pos_start is not None and gene in self._regions_by_gene:
                    v_end = pos_end if pos_end is not None else pos_start
                    for r_idx, r_start, r_end in self._regions_by_gene[gene]:
                        if pos_start <= r_end and v_end >= r_start:
                            r_key = (r_idx, change, bucket)
                            if r_key not in seen_reg_variants:
                                seen_reg_variants.add(r_key)
                                X[row_i, region_base + r_idx * n_cons + ci] += 1.0
        return X

    def fit(self, train_isolates: Sequence[Sequence[dict[str, Any]]]) -> "BurdenEncoder":
        _need(len(train_isolates) > 0, "BurdenEncoder.fit requires at least one training isolate")
        raw = self._raw_counts(train_isolates)
        if self.standardize:
            self.mean_, self.scale_ = _fit_standardizer(raw)
        else:
            self.mean_ = np.zeros(raw.shape[1], dtype=np.float64)
            self.scale_ = np.ones(raw.shape[1], dtype=np.float64)
        self.fitted_ = True
        return self

    def transform(self, isolates: Sequence[Sequence[dict[str, Any]]]) -> np.ndarray:
        _need(self.fitted_ and self.mean_ is not None and self.scale_ is not None,
              "BurdenEncoder must be fitted before transform")
        raw = self._raw_counts(isolates)
        if not self.standardize:
            return raw
        out = (raw - self.mean_) * self.scale_
        return np.where(out == 0.0, 0.0, out)

    def fit_transform(self, train_isolates: Sequence[Sequence[dict[str, Any]]]) -> np.ndarray:
        return self.fit(train_isolates).transform(train_isolates)

    def state_fingerprint(self) -> str:
        _need(self.fitted_ and self.mean_ is not None and self.scale_ is not None,
              "BurdenEncoder must be fitted before fingerprinting")
        h = hashlib.sha256()
        h.update(_canonical_json_bytes({"feature_names": list(self.feature_names_), "standardize": self.standardize}))
        h.update(np.ascontiguousarray(self.mean_, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(self.scale_, dtype="<f8").tobytes())
        return h.hexdigest()


def load_foundation_manifest(
    manifest_or_path: dict[str, Any] | str | Path,
    base_dir: str | Path | None = None,
) -> tuple[dict[str, Any], np.ndarray, dict[str, Any]]:
    """Load and validate a `grm-rep-foundation-v1` manifest and its exposure audit.

    Returns `(manifest_meta, raw_values_float64, validated_exposure_audit)`.
    """
    if isinstance(manifest_or_path, (str, Path)):
        mpath = Path(manifest_or_path)
        manifest = json.loads(mpath.read_text(encoding="utf-8"))
        root = mpath.parent if base_dir is None else Path(base_dir)
    else:
        manifest = manifest_or_path
        root = Path(base_dir) if base_dir is not None else Path(".")

    _need(isinstance(manifest, dict), "foundation manifest must be a JSON object")
    _need(
        manifest.get("format") == FOUNDATION_MANIFEST_FORMAT,
        f"foundation manifest format must be {FOUNDATION_MANIFEST_FORMAT!r}, got {manifest.get('format')!r}",
    )
    if "exposure_audit" not in manifest or manifest["exposure_audit"] is None:
        raise AuditError("foundation manifest is missing required 'exposure_audit' record")
    validated_audit = EA.validate_exposure_audit(manifest["exposure_audit"])

    rep_kind = manifest.get("representation_kind")
    _need(
        isinstance(rep_kind, str) and rep_kind in ALLOWED_FOUNDATION_KINDS,
        f"invalid representation_kind: {rep_kind!r}",
    )
    dataset = manifest.get("dataset")
    _need(isinstance(dataset, str) and bool(dataset.strip()), "foundation manifest requires non-empty 'dataset'")

    samples = manifest.get("samples")
    _need(isinstance(samples, list) and len(samples) > 0, "foundation manifest requires non-empty 'samples'")
    _need(all(isinstance(s, str) and bool(s.strip()) for s in samples), "all manifest sample ids must be non-empty strings")
    _need(len(set(samples)) == len(samples), "foundation manifest lists duplicate sample ids")

    feature_names = manifest.get("feature_names")
    _need(isinstance(feature_names, list) and len(feature_names) > 0,
          "foundation manifest requires non-empty 'feature_names'")
    _need(all(isinstance(f, str) and bool(f.strip()) for f in feature_names),
          "all manifest feature_names must be non-empty strings")
    _need(len(set(feature_names)) == len(feature_names), "foundation manifest lists duplicate feature_names")

    has_inline = "values" in manifest and manifest["values"] is not None
    has_npy = "values_npy" in manifest and manifest["values_npy"] is not None
    _need(has_inline != has_npy, "foundation manifest must specify exactly one of 'values' or 'values_npy'")

    if has_inline:
        raw_rows = manifest["values"]
        _need(isinstance(raw_rows, list) and len(raw_rows) == len(samples),
              "foundation manifest 'values' row count must match len(samples)")
        rows_float: list[list[float]] = []
        for i, row in enumerate(raw_rows):
            _need(isinstance(row, list) and len(row) == len(feature_names),
                  f"foundation manifest row {i} length does not match len(feature_names)")
            converted: list[float] = []
            for j, cell in enumerate(row):
                if cell is None:
                    converted.append(float("nan"))
                else:
                    _need(isinstance(cell, (int, float)) and not isinstance(cell, bool),
                          f"foundation manifest values[{i}][{j}] must be numeric or null")
                    converted.append(float(cell))
            rows_float.append(converted)
        arr = np.asarray(rows_float, dtype=np.float64)
    else:
        npy_rel = manifest["values_npy"]
        _need(isinstance(npy_rel, str) and bool(npy_rel.strip()), "values_npy must be a non-empty relative path")
        npy_path = (root / npy_rel).resolve()
        _need(npy_path.is_file(), f"values_npy file does not exist: {npy_path}")
        expected_sha = manifest.get("values_sha256")
        _need(
            isinstance(expected_sha, str) and bool(EA.SHA256_RE.fullmatch(expected_sha)),
            "values_npy requires a valid 64-hex values_sha256 digest",
        )
        actual_sha = hashlib.sha256(npy_path.read_bytes()).hexdigest()
        _need(actual_sha == expected_sha, f"values_npy sha256 mismatch: expected {expected_sha}, got {actual_sha}")
        loaded = np.load(npy_path, allow_pickle=False)
        _need(loaded.ndim == 2 and loaded.shape == (len(samples), len(feature_names)),
              f"values_npy shape {loaded.shape} != {(len(samples), len(feature_names))}")
        arr = np.asarray(loaded, dtype=np.float64)

    _need(not bool(np.any(np.isinf(arr))), "foundation manifest values cannot contain +/-inf")

    meta = {
        "format": FOUNDATION_MANIFEST_FORMAT,
        "dataset": dataset.strip(),
        "representation_kind": rep_kind,
        "samples": list(samples),
        "feature_names": list(feature_names),
    }
    return meta, arr, validated_audit


def align_samples_to_order(
    manifest_meta: dict[str, Any],
    values: np.ndarray,
    expected_dataset: str,
    target_samples: Sequence[str],
) -> np.ndarray:
    """Align manifest rows to `target_samples`, requiring exact 1-to-1 sample set match."""
    _need(
        manifest_meta["dataset"] == expected_dataset,
        f"foundation manifest dataset {manifest_meta['dataset']!r} != {expected_dataset!r}",
    )
    src_samples = manifest_meta["samples"]
    _need(len(target_samples) == len(set(target_samples)), "target_samples contains duplicates")
    src_index = {s: i for i, s in enumerate(src_samples)}
    missing = [s for s in target_samples if s not in src_index]
    extra = [s for s in src_samples if s not in set(target_samples)]
    _need(not missing and not extra,
          f"foundation manifest sample mismatch: missing={len(missing)}, extra={len(extra)}")
    order = np.array([src_index[s] for s in target_samples], dtype=np.int64)
    return np.ascontiguousarray(values[order], dtype=np.float64)


def fold_file_samples(fold_doc: dict[str, Any]) -> list[str]:
    """Extract canonical ordered cohort sample IDs from a `grm-fold-file-v1` document."""
    _need(isinstance(fold_doc, dict) and fold_doc.get("format") == "grm-fold-file-v1",
          "not a grm-fold-file-v1 document")
    listed: list[str] = []
    for fold_ids in fold_doc.get("test", ()):
        listed.extend(fold_ids)
    listed.extend(fold_doc.get("train", ()))
    listed.extend(fold_doc.get("unused", ()))
    _need(len(listed) == len(set(listed)), "grm-fold-file-v1 lists a sample id more than once")
    return sorted(listed)


class FoundationScoreAdapter:
    """Fold-safe adapter for precomputed foundation-model scores/embeddings.

    Ingests a `grm-rep-foundation-v1` manifest, enforces `grm-rep-exposure-v1`
    validation, aligns rows strictly to `dataset.json` or `grm-fold-file-v1`,
    and fits optional mean imputation, standardization, and deterministic PCA
    projection strictly on outer training rows.
    """

    def __init__(
        self,
        manifest_or_path: dict[str, Any] | str | Path,
        *,
        export_meta: dict[str, Any] | None = None,
        fold_doc: dict[str, Any] | None = None,
        sample_order: Sequence[str] | None = None,
        base_dir: str | Path | None = None,
        standardize: bool = True,
        pca_components: int | None = None,
        impute_missing: bool = False,
    ) -> None:
        meta, raw, audit = load_foundation_manifest(manifest_or_path, base_dir=base_dir)
        self.exposure_audit_: dict[str, Any] = audit
        self.claim_ceiling_: str = audit["effective_claim_ceiling"]
        self.manifest_meta_: dict[str, Any] = meta
        self.standardize = bool(standardize)
        self.pca_components = pca_components
        if pca_components is not None:
            _need(_is_strict_int(pca_components) and pca_components >= 1, "pca_components must be integer >= 1")
        self.impute_missing = bool(impute_missing)

        if export_meta is not None:
            _need(export_meta.get("format") == "grm-dataset-export-v1", "export_meta is not grm-dataset-export-v1")
            target_dataset = export_meta["id"]
            target_samples = list(export_meta["samples"])
            if fold_doc is not None:
                _need(fold_doc.get("format") == "grm-fold-file-v1", "fold_doc is not grm-fold-file-v1")
                _need(fold_doc["dataset"] == target_dataset, "fold_doc dataset differs from export_meta")
                _need(set(fold_file_samples(fold_doc)) == set(target_samples),
                      "fold_doc samples differ from export_meta samples")
        elif fold_doc is not None:
            target_dataset = fold_doc["dataset"]
            target_samples = list(sample_order) if sample_order is not None else fold_file_samples(fold_doc)
            _need(set(target_samples) == set(fold_file_samples(fold_doc)),
                  "sample_order does not match grm-fold-file-v1 cohort")
        else:
            target_dataset = meta["dataset"]
            target_samples = list(sample_order) if sample_order is not None else list(meta["samples"])

        aligned = align_samples_to_order(meta, raw, target_dataset, target_samples)
        if not self.impute_missing:
            _need(bool(np.all(np.isfinite(aligned))),
                  "foundation manifest contains missing/NaN values with impute_missing=False")
        self.samples_: tuple[str, ...] = tuple(target_samples)
        self._raw_matrix: np.ndarray = aligned
        self.raw_feature_names_: tuple[str, ...] = tuple(meta["feature_names"])

        self.fitted_: bool = False
        self.impute_means_: np.ndarray | None = None
        self.mean_: np.ndarray | None = None
        self.scale_: np.ndarray | None = None
        self.components_: np.ndarray | None = None
        self.pca_center_: np.ndarray | None = None
        self.feature_names_: tuple[str, ...] = ()

    @property
    def n_samples(self) -> int:
        return int(self._raw_matrix.shape[0])

    def _select_rows(self, indices: Sequence[int] | np.ndarray) -> np.ndarray:
        idx = np.asarray(indices, dtype=np.int64)
        _need(idx.ndim == 1, "indices must be a 1D index sequence")
        if idx.size > 0:
            _need(int(np.min(idx)) >= 0 and int(np.max(idx)) < self.n_samples, "sample index out of range")
        return np.ascontiguousarray(self._raw_matrix[idx], dtype=np.float64)

    def _apply_imputation(self, X: np.ndarray, means: np.ndarray) -> np.ndarray:
        if not np.any(np.isnan(X)):
            return X
        out = X.copy()
        rows, cols = np.where(np.isnan(out))
        out[rows, cols] = means[cols]
        return out

    def fit(self, train_indices: Sequence[int] | np.ndarray) -> "FoundationScoreAdapter":
        X = self._select_rows(train_indices)
        _need(X.shape[0] > 0, "FoundationScoreAdapter.fit requires at least one training sample")

        # Fit missing-value imputation strictly on training rows.
        if self.impute_missing:
            col_means = np.zeros(X.shape[1], dtype=np.float64)
            for j in range(X.shape[1]):
                col = X[:, j]
                valid = col[~np.isnan(col)]
                col_means[j] = float(np.mean(valid)) if valid.size > 0 else 0.0
            self.impute_means_ = col_means
            X = self._apply_imputation(X, col_means)
        else:
            self.impute_means_ = np.zeros(X.shape[1], dtype=np.float64)
            _need(bool(np.all(np.isfinite(X))), "non-finite values in training foundation matrix")

        if self.standardize:
            self.mean_, self.scale_ = _fit_standardizer(X)
            Xs = (X - self.mean_) * self.scale_
        else:
            self.mean_ = np.zeros(X.shape[1], dtype=np.float64)
            self.scale_ = np.ones(X.shape[1], dtype=np.float64)
            Xs = X

        if self.pca_components is not None:
            self.pca_center_ = np.mean(Xs, axis=0, dtype=np.float64)
            self.components_ = _deterministic_pca_components(Xs, self.pca_components)
            self.feature_names_ = tuple(f"pc_{i}" for i in range(self.components_.shape[0]))
        else:
            self.pca_center_ = None
            self.components_ = None
            self.feature_names_ = self.raw_feature_names_

        self.fitted_ = True
        return self

    def transform(self, indices: Sequence[int] | np.ndarray) -> np.ndarray:
        _need(
            self.fitted_ and self.mean_ is not None and self.scale_ is not None and self.impute_means_ is not None,
            "FoundationScoreAdapter must be fitted before transform",
        )
        X = self._select_rows(indices)
        if self.impute_missing:
            X = self._apply_imputation(X, self.impute_means_)
        else:
            _need(bool(np.all(np.isfinite(X))), "non-finite values in foundation matrix")
        Xs = (X - self.mean_) * self.scale_
        if self.components_ is not None and self.pca_center_ is not None:
            out = (Xs - self.pca_center_) @ self.components_.T
        else:
            out = Xs
        return np.where(out == 0.0, 0.0, out)

    def fit_transform(self, train_indices: Sequence[int] | np.ndarray) -> np.ndarray:
        return self.fit(train_indices).transform(train_indices)

    def clone_unfitted(self) -> "FoundationScoreAdapter":
        """Create a fresh unfitted adapter sharing the validated immutable raw matrix."""
        obj = object.__new__(FoundationScoreAdapter)
        obj.exposure_audit_ = self.exposure_audit_
        obj.claim_ceiling_ = self.claim_ceiling_
        obj.manifest_meta_ = self.manifest_meta_
        obj.standardize = self.standardize
        obj.pca_components = self.pca_components
        obj.impute_missing = self.impute_missing
        obj.samples_ = self.samples_
        obj._raw_matrix = self._raw_matrix
        obj.raw_feature_names_ = self.raw_feature_names_
        obj.fitted_ = False
        obj.impute_means_ = None
        obj.mean_ = None
        obj.scale_ = None
        obj.components_ = None
        obj.pca_center_ = None
        obj.feature_names_ = ()
        return obj

    def state_fingerprint(self) -> str:
        _need(
            self.fitted_ and self.mean_ is not None and self.scale_ is not None and self.impute_means_ is not None,
            "FoundationScoreAdapter must be fitted before fingerprinting",
        )
        h = hashlib.sha256()
        h.update(_canonical_json_bytes({
            "feature_names": list(self.feature_names_),
            "standardize": self.standardize,
            "pca_components": self.pca_components,
            "impute_missing": self.impute_missing,
        }))
        h.update(np.ascontiguousarray(self.impute_means_, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(self.mean_, dtype="<f8").tobytes())
        h.update(np.ascontiguousarray(self.scale_, dtype="<f8").tobytes())
        if self.pca_center_ is not None:
            h.update(np.ascontiguousarray(self.pca_center_, dtype="<f8").tobytes())
        if self.components_ is not None:
            h.update(np.ascontiguousarray(self.components_, dtype="<f8").tobytes())
        return h.hexdigest()
