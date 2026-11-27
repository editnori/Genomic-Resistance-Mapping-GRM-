"""Canonical antibiotic vocabulary parsed directly from src/drugs.rs.

Every adapter drug map and species-support registry validates its canonical
drug names against the Rust source of truth (src/drugs.rs) so Python harness
tables cannot drift from GRM's vocabulary, inhibitor ordering, or drug classes.
"""

from __future__ import annotations

import functools
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DRUGS_RS = REPO_ROOT / "src" / "drugs.rs"

_CONST_SLICE_RE = re.compile(
    r"const\s+(?P<name>[A-Z0-9_]+)\s*:\s*&\[&str\]\s*=\s*&\[(?P<body>.*?)\];",
    re.DOTALL,
)
_OVERRIDES_RE = re.compile(
    r"const\s+OVERRIDES\s*:\s*&\[\(&str,\s*&str\)\]\s*=\s*&\[(?P<body>.*?)\];",
    re.DOTALL,
)
_SINGLE_CLASS_RE = re.compile(
    r"fn\s+single_class\s*\([^)]*\)\s*->\s*Option<&'static\s+str>\s*\{(?P<body>.*?)\n\}",
    re.DOTALL,
)
_QUOTED_RE = re.compile(r'"([^"]+)"')
_TUPLE_RE = re.compile(r'\(\s*"([^"]+)"\s*,\s*"([^"]+)"\s*\)')


def _extract_slice(text: str, const_name: str) -> tuple[str, ...]:
    for m in _CONST_SLICE_RE.finditer(text):
        if m.group("name") == const_name:
            items = tuple(_QUOTED_RE.findall(m.group("body")))
            if not items:
                raise ValueError(f"src/drugs.rs: const {const_name} parsed empty")
            return items
    raise ValueError(f"src/drugs.rs: missing const {const_name}")


def _extract_overrides(text: str) -> dict[str, str]:
    m = _OVERRIDES_RE.search(text)
    if m is None:
        raise ValueError("src/drugs.rs: missing const OVERRIDES")
    pairs = _TUPLE_RE.findall(m.group("body"))
    if not pairs:
        raise ValueError("src/drugs.rs: OVERRIDES parsed empty")
    return dict(pairs)


def _extract_single_classes(text: str) -> dict[str, str]:
    m = _SINGLE_CLASS_RE.search(text)
    if m is None:
        raise ValueError("src/drugs.rs: missing fn single_class")
    body = m.group("body")
    # Match arms of the form: "a" | "b" => "class" or => { "class" }
    arm_re = re.compile(r'(?P<lhs>(?:"[^"]+"\s*\|\s*)*"[^"]+")\s*=>\s*\{?\s*"(?P<cls>[^"]+)"\s*\}?', re.DOTALL)
    out: dict[str, str] = {}
    for am in arm_re.finditer(body):
        cls = am.group("cls")
        for drug in _QUOTED_RE.findall(am.group("lhs")):
            out[drug] = cls
    if not out:
        raise ValueError("src/drugs.rs: single_class parsed empty")
    return out


class DrugVocabulary:
    """Canonical drug vocabulary loaded directly from src/drugs.rs."""

    def __init__(
        self,
        drugs: frozenset[str],
        inhibitors: frozenset[str],
        overrides: dict[str, str],
        single_classes: dict[str, str],
        source_path: Path,
    ) -> None:
        self.drugs = drugs
        self.inhibitors = inhibitors
        self.overrides = overrides
        self.single_classes = single_classes
        self.source_path = source_path

    def is_canonical(self, name: str) -> bool:
        """True iff `name` is a valid canonical GRM antibiotic or combination.

        Combinations list principal drugs alphabetically and beta-lactamase
        inhibitors alphabetically last, joined by '/'.
        """
        if not isinstance(name, str) or not name or name != name.strip():
            return False
        if not name.isascii() or name != name.lower():
            return False
        parts = name.split("/")
        if any(not p or p != p.strip() for p in parts):
            return False
        if len(set(parts)) != len(parts):
            return False
        principals = [p for p in parts if p not in self.inhibitors]
        inhibitors = [p for p in parts if p in self.inhibitors]
        if not principals:
            return False
        if any(p not in self.drugs for p in principals):
            return False
        if any(i not in self.inhibitors for i in inhibitors):
            return False
        return parts == sorted(principals) + sorted(inhibitors)

    def canonicalize(self, raw: str) -> str:
        """Normalize a known antibiotic/combination string using src/drugs.rs rules."""
        cleaned = " ".join(
            "".join(ch for ch in raw if ord(ch) < 128).strip().lower().replace("_", " ").split()
        ).replace(" / ", "/")
        if cleaned in self.overrides:
            cleaned = self.overrides[cleaned]
        parts = [p.strip() for p in cleaned.split("/") if p.strip()]
        expanded: list[str] = []
        for p in parts:
            mapped = self.overrides.get(p, p)
            expanded.extend(mapped.split("/"))
        principals = sorted({p for p in expanded if p not in self.inhibitors})
        inhibitors = sorted({p for p in expanded if p in self.inhibitors})
        canon = "/".join(principals + inhibitors)
        if not self.is_canonical(canon):
            raise ValueError(f"not a canonical GRM antibiotic in src/drugs.rs: {raw!r} -> {canon!r}")
        return canon

    def class_of(self, canonical: str) -> str:
        """Mirror src/drugs.rs::class_of for a canonical drug name."""
        if not self.is_canonical(canonical):
            raise ValueError(f"not a canonical GRM antibiotic: {canonical!r}")
        parts = canonical.split("/")
        inhibited = any(p in self.inhibitors for p in parts)
        principals = [p for p in parts if p not in self.inhibitors]
        classes = [self.single_classes[p] for p in principals if p in self.single_classes]
        if not classes:
            base = "other"
        elif len(classes) == 1:
            base = classes[0]
        elif all(c in ("trimethoprim", "sulfonamide") for c in classes):
            base = "folate pathway combination"
        elif all(c == classes[0] for c in classes):
            base = classes[0]
        else:
            base = " + ".join(classes)
        if inhibited:
            family = "cephalosporin" if base.startswith("cephalosporin") else base
            return f"{family} + inhibitor"
        return base


@functools.lru_cache(maxsize=4)
def load_vocabulary(drugs_rs: Path | str = DRUGS_RS) -> DrugVocabulary:
    path = Path(drugs_rs).resolve()
    text = path.read_text(encoding="utf-8")
    drugs = frozenset(_extract_slice(text, "GRM_VOCABULARY"))
    inhibitors = frozenset(_extract_slice(text, "INHIBITORS"))
    overrides = _extract_overrides(text)
    single_classes = _extract_single_classes(text)
    return DrugVocabulary(
        drugs=drugs,
        inhibitors=inhibitors,
        overrides=overrides,
        single_classes=single_classes,
        source_path=path,
    )


def drug_slug(canonical_drug: str) -> str:
    """Convert a canonical drug name to its GRM dataset-id drug suffix."""
    return re.sub(r"[/\-\s]+", "_", canonical_drug.strip().lower())


def dataset_id(species_slug: str, canonical_drug: str) -> str:
    """Build canonical GRM dataset id '<species_slug>__<drug_slug>'."""
    return f"{species_slug}__{drug_slug(canonical_drug)}"
