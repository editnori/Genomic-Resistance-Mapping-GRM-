#!/usr/bin/env python3
"""Standard-library metadata normalizers for public AMR portals.

Converts public antimicrobial-resistance portal metadata tables into:
1. Canonical `grm dedup` query JSONL (one deduplicated isolate per line, sorted
   deterministically by genome id).
2. Optional INSDC/ENA/NCBI synonym TSV (`kind`, `accession_a`, `accession_b`)
   for `grm dedup screen --synonym-map`, linking `SAMEA`<->`ERS`,
   `SAMN`<->`SRS`, `SAMD`<->`DRS`, run-to-sample links (`ERR`/`SRR`/`DRR`),
   `PRJEB`<->`ERP`, `PRJNA`<->`SRP`, `PRJDB`<->`DRP`, and paired GenBank/RefSeq
   assemblies (`GCA_`<->`GCF_`).

Supported portals:
- `ebi_amr` (aliases: `cabbage`, `ebi-amr`): EMBL-EBI AMR Portal / CABBAGE
- `cryptic`: CRyPTIC Mycobacterium tuberculosis compendium
- `ncbi_pd` (aliases: `narms`, `ncbi-pd`): NARMS / NCBI Pathogen Detection
- `enterobase`: EnteroBase strain and assembly exports
- `pubmlst`: PubMLST isolate metadata exports
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import pathlib
import re
import sys
from collections.abc import Iterable, Iterator, Mapping

PLACEHOLDERS = frozenset({
    "",
    "-",
    ".",
    "NA",
    "N/A",
    "NONE",
    "NULL",
    "MISSING",
    "UNKNOWN",
    "NOT APPLICABLE",
    "NOT COLLECTED",
    "NOT PROVIDED",
    "UNAVAILABLE",
})

PORTAL_ALIASES = {
    "ebi_amr": "ebi_amr",
    "ebi-amr": "ebi_amr",
    "cabbage": "ebi_amr",
    "cryptic": "cryptic",
    "ncbi_pd": "ncbi_pd",
    "ncbi-pd": "ncbi_pd",
    "narms": "ncbi_pd",
    "enterobase": "enterobase",
    "pubmlst": "pubmlst",
}

SYNONYM_HEADER = ("kind", "accession_a", "accession_b")

_SPLIT_RE = re.compile(r"[,;|\s]+")
_PMID_NUM_RE = re.compile(r"^(?:PMID:)?(\d+)$", re.IGNORECASE)


def canonical_portal(name: str) -> str:
    key = name.strip().lower().replace(" ", "_")
    if key not in PORTAL_ALIASES:
        supported = ", ".join(sorted(set(PORTAL_ALIASES.values())))
        raise ValueError(f"unsupported portal {name!r}; expected one of: {supported}")
    return PORTAL_ALIASES[key]


def genome_id_sort_key(genome_id: str) -> tuple[int, int, int, str]:
    """Matches Rust `grm::dedup::compare_genome_ids` total order."""
    if "." in genome_id:
        left, right = genome_id.split(".", 1)
        if left.isdigit() and right.isdigit():
            return (0, int(left), int(right), genome_id)
    return (1, 0, 0, genome_id)


def _clean_text(value: object) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    if text.upper() in PLACEHOLDERS:
        return None
    return text


def _split_tokens(value: object, strip_prefix: str = "") -> list[str]:
    text = _clean_text(value)
    if text is None:
        return []
    out: list[str] = []
    prefix_up = strip_prefix.upper()
    for raw in _SPLIT_RE.split(text):
        tok = raw.strip().upper()
        if prefix_up and tok.startswith(prefix_up):
            tok = tok[len(prefix_up):]
        if tok and tok not in PLACEHOLDERS:
            out.append(tok)
    return out


def _clean_pmids(value: object) -> list[str]:
    out: list[str] = []
    for tok in _split_tokens(value, strip_prefix="PMID:"):
        m = _PMID_NUM_RE.match(tok)
        if m:
            out.append(m.group(1))
        elif tok.isdigit():
            out.append(tok)
    return out


def _collapse_ws(value: object) -> str | None:
    text = _clean_text(value)
    if text is None:
        return None
    collapsed = " ".join(text.split())
    if collapsed.upper() in PLACEHOLDERS:
        return None
    return collapsed


def _normalize_species(value: object, default: str | None = None) -> str | None:
    text = _collapse_ws(value)
    if text is None:
        return default
    words = text.split()
    if len(words) >= 2 and words[0][0].isupper():
        return f"{words[0]} {words[1]}"
    return text


def _ci_lookup(row: Mapping[str, object]) -> dict[str, object]:
    return {k.strip().lower(): v for k, v in row.items() if k is not None}


def _first_col(ci_row: Mapping[str, object], names: Iterable[str]) -> str | None:
    for name in names:
        val = _clean_text(ci_row.get(name.lower()))
        if val is not None:
            return val
    return None


def _all_col_tokens(
    ci_row: Mapping[str, object],
    names: Iterable[str],
    strip_prefix: str = "",
) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for name in names:
        for tok in _split_tokens(ci_row.get(name.lower()), strip_prefix=strip_prefix):
            if tok not in seen:
                seen.add(tok)
                out.append(tok)
    return out


def _open_text(path: pathlib.Path) -> io.TextIOBase:
    with path.open("rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, mode="rt", encoding="utf-8-sig", newline="")
    return path.open(mode="r", encoding="utf-8-sig", newline="")


def read_table_rows(path: pathlib.Path, delimiter: str | None = None) -> Iterator[dict[str, object]]:
    """Reads TSV, CSV, or JSONL (plain or `.gz`) rows as dictionaries."""
    name_low = path.name.lower()
    if name_low.endswith(".gz"):
        stem_low = name_low[:-3]
    else:
        stem_low = name_low
    if stem_low.endswith(".jsonl") or stem_low.endswith(".ndjson"):
        with _open_text(path) as handle:
            for line_no, raw in enumerate(handle, 1):
                line = raw.strip()
                if not line:
                    continue
                obj = json.loads(line)
                if not isinstance(obj, dict):
                    raise ValueError(f"{path}:{line_no}: expected JSON object")
                yield obj
        return

    with _open_text(path) as handle:
        first_line = handle.readline()
        if not first_line:
            return
        if delimiter is None:
            if stem_low.endswith(".csv"):
                delim = ","
            elif stem_low.endswith(".tsv") or stem_low.endswith(".tab"):
                delim = "\t"
            else:
                delim = "\t" if "\t" in first_line else ","
        else:
            delim = delimiter
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delim)
        if reader.fieldnames is None:
            return
        for row in reader:
            yield dict(row)


def _extract_ebi_amr(ci_row: Mapping[str, object]) -> tuple[dict[str, object], list[tuple[str, str, str]]]:
    biosamples = _all_col_tokens(
        ci_row,
        ("biosample_id", "biosample_accession", "biosample"),
        strip_prefix="BIOSAMPLE:",
    )
    secondary_samples = _all_col_tokens(
        ci_row,
        ("secondary_sample_accession", "sample_accession", "ena_sample"),
    )
    runs = _all_col_tokens(
        ci_row,
        ("run_accession", "sra_accession", "ena_run", "run"),
    )
    assemblies = _all_col_tokens(
        ci_row,
        ("assembly_accession", "gca_accession", "genbank_assembly", "refseq_assembly"),
    )
    bioprojects = _all_col_tokens(
        ci_row,
        ("bioproject_accession", "bioproject", "project_accession"),
        strip_prefix="BIOPROJECT:",
    )
    secondary_studies = _all_col_tokens(
        ci_row,
        ("secondary_study_accession", "study_accession", "ena_study"),
    )
    genome_id = (
        _first_col(ci_row, ("genome_id", "isolate_id", "sample_id"))
        or (biosamples[0] if biosamples else None)
        or (secondary_samples[0] if secondary_samples else None)
        or (runs[0] if runs else None)
        or (assemblies[0] if assemblies else None)
    )
    if genome_id is None:
        raise ValueError("ebi_amr row has no usable identifier (BioSample, sample, run, or assembly)")

    species = _normalize_species(
        _first_col(ci_row, ("species", "scientific_name", "organism")),
    )
    strain = _collapse_ws(_first_col(ci_row, ("strain", "isolate", "sample_alias")))
    subject_id = _collapse_ws(_first_col(ci_row, ("subject_id", "patient_id", "host_subject_id")))
    fasta_path = _first_col(ci_row, ("fasta_path", "assembly_path"))
    pmids: list[str] = []
    for col in ("pmid", "pubmed_id", "publication"):
        pmids.extend(_clean_pmids(ci_row.get(col)))

    synonyms: list[tuple[str, str, str]] = []
    for b in biosamples:
        for s in secondary_samples:
            if b != s:
                synonyms.append(("sample", b, s))
    sample_anchors = biosamples or secondary_samples
    for r in runs:
        for anchor in sample_anchors:
            if r != anchor:
                synonyms.append(("run_to_sample", r, anchor))
    for bp in bioprojects:
        for st in secondary_studies:
            if bp != st:
                synonyms.append(("study", bp, st))

    rec = {
        "genome_id": genome_id,
        "portal": "ebi_amr",
        "species": species,
        "biosample_accessions": biosamples + [s for s in secondary_samples if s not in biosamples],
        "assembly_accessions": assemblies,
        "sra_accessions": runs,
        "bioproject_accessions": bioprojects + [s for s in secondary_studies if s not in bioprojects],
        "pmids": pmids,
        "strain": strain,
        "subject_id": subject_id,
        "fasta_path": fasta_path,
    }
    return rec, synonyms


def _extract_cryptic(ci_row: Mapping[str, object]) -> tuple[dict[str, object], list[tuple[str, str, str]]]:
    biosamples = _all_col_tokens(
        ci_row,
        ("biosample", "biosample_accession", "biosample_id"),
        strip_prefix="BIOSAMPLE:",
    )
    ena_samples = _all_col_tokens(
        ci_row,
        ("ena_sample", "ers_accession", "secondary_sample_accession", "sample_accession"),
    )
    runs = _all_col_tokens(
        ci_row,
        ("ena_run", "run_accession", "sra_accession", "err_accession"),
    )
    assemblies = _all_col_tokens(
        ci_row,
        ("assembly_accession", "gca_accession", "genbank_assembly"),
    )
    bioprojects = _all_col_tokens(
        ci_row,
        ("bioproject", "bioproject_accession"),
        strip_prefix="BIOPROJECT:",
    )
    ena_studies = _all_col_tokens(
        ci_row,
        ("ena_study", "study_accession", "secondary_study_accession"),
    )
    genome_id = (
        _first_col(ci_row, ("uniqueid", "genome_id", "sample_id"))
        or (biosamples[0] if biosamples else None)
        or (ena_samples[0] if ena_samples else None)
        or (runs[0] if runs else None)
    )
    if genome_id is None:
        raise ValueError("cryptic row has no usable identifier (UNIQUEID, BIOSAMPLE, ENA_SAMPLE, or ENA_RUN)")

    species = _normalize_species(
        _first_col(ci_row, ("species", "organism", "scientific_name")),
        default="Mycobacterium tuberculosis",
    )
    site_id = _collapse_ws(_first_col(ci_row, ("siteid", "site_id")))
    subj_raw = _collapse_ws(_first_col(ci_row, ("subjid", "subject_id", "patient_id")))
    if site_id and subj_raw and ":" not in subj_raw:
        subject_id = f"{site_id}:{subj_raw}"
    else:
        subject_id = subj_raw
    strain = _collapse_ws(_first_col(ci_row, ("strain", "isolate")))
    fasta_path = _first_col(ci_row, ("fasta_path", "assembly_path"))
    pmids: list[str] = []
    for col in ("pmid", "pubmed_id", "publication"):
        pmids.extend(_clean_pmids(ci_row.get(col)))

    synonyms: list[tuple[str, str, str]] = []
    for b in biosamples:
        for e in ena_samples:
            if b != e:
                synonyms.append(("sample", b, e))
    anchors = biosamples or ena_samples
    for r in runs:
        for a in anchors:
            if r != a:
                synonyms.append(("run_to_sample", r, a))
    for bp in bioprojects:
        for st in ena_studies:
            if bp != st:
                synonyms.append(("study", bp, st))

    rec = {
        "genome_id": genome_id,
        "portal": "cryptic",
        "species": species,
        "biosample_accessions": biosamples + [e for e in ena_samples if e not in biosamples],
        "assembly_accessions": assemblies,
        "sra_accessions": runs,
        "bioproject_accessions": bioprojects + [s for s in ena_studies if s not in bioprojects],
        "pmids": pmids,
        "strain": strain,
        "subject_id": subject_id,
        "fasta_path": fasta_path,
    }
    return rec, synonyms


def _extract_ncbi_pd(ci_row: Mapping[str, object]) -> tuple[dict[str, object], list[tuple[str, str, str]]]:
    biosamples = _all_col_tokens(
        ci_row,
        ("biosample_acc", "biosample", "biosample_accession"),
        strip_prefix="BIOSAMPLE:",
    )
    secondary_samples = _all_col_tokens(
        ci_row,
        ("sra_sample", "secondary_sample_accession", "srs_accession"),
    )
    genbank_asm = _all_col_tokens(
        ci_row,
        ("asm_acc", "assembly", "assembly_accession", "genbank_assembly_accession"),
    )
    refseq_asm = _all_col_tokens(
        ci_row,
        ("refseq_asm_acc", "refseq_assembly", "refseq_assembly_accession"),
    )
    runs = _all_col_tokens(
        ci_row,
        ("run", "sra_acc", "sra_run", "sra_accession", "run_accession"),
    )
    bioprojects = _all_col_tokens(
        ci_row,
        ("bioproject_acc", "bioproject", "bioproject_accession"),
        strip_prefix="BIOPROJECT:",
    )
    secondary_studies = _all_col_tokens(
        ci_row,
        ("sra_study", "secondary_study_accession", "study_accession"),
    )
    genome_id = (
        _first_col(ci_row, ("genome_id", "target_acc", "pdt_acc", "isolate_acc"))
        or (biosamples[0] if biosamples else None)
        or (genbank_asm[0] if genbank_asm else None)
        or (runs[0] if runs else None)
    )
    if genome_id is None:
        raise ValueError("ncbi_pd row has no usable identifier (target_acc, BioSample, Assembly, or Run)")

    species = _normalize_species(
        _first_col(ci_row, ("scientific_name", "species", "organism_group", "organism")),
    )
    strain = _collapse_ws(_first_col(ci_row, ("strain", "isolate", "isolate_identifiers")))
    subject_id = _collapse_ws(_first_col(ci_row, ("host_subject_id", "subject_id", "patient_id")))
    fasta_path = _first_col(ci_row, ("fasta_path", "assembly_path"))
    pmids: list[str] = []
    for col in ("pmid", "pubmed_id", "publication"):
        pmids.extend(_clean_pmids(ci_row.get(col)))

    synonyms: list[tuple[str, str, str]] = []
    for b in biosamples:
        for s in secondary_samples:
            if b != s:
                synonyms.append(("sample", b, s))
    anchors = biosamples or secondary_samples
    for r in runs:
        for a in anchors:
            if r != a:
                synonyms.append(("run_to_sample", r, a))
    for g in genbank_asm:
        for rf in refseq_asm:
            if g != rf:
                synonyms.append(("assembly", g, rf))
    for bp in bioprojects:
        for st in secondary_studies:
            if bp != st:
                synonyms.append(("study", bp, st))

    assemblies = genbank_asm + [a for a in refseq_asm if a not in genbank_asm]
    rec = {
        "genome_id": genome_id,
        "portal": "ncbi_pd",
        "species": species,
        "biosample_accessions": biosamples + [s for s in secondary_samples if s not in biosamples],
        "assembly_accessions": assemblies,
        "sra_accessions": runs,
        "bioproject_accessions": bioprojects + [s for s in secondary_studies if s not in bioprojects],
        "pmids": pmids,
        "strain": strain,
        "subject_id": subject_id,
        "fasta_path": fasta_path,
    }
    return rec, synonyms


def _extract_enterobase(ci_row: Mapping[str, object]) -> tuple[dict[str, object], list[tuple[str, str, str]]]:
    biosamples = _all_col_tokens(
        ci_row,
        ("sample id", "biosample", "biosample_accession", "biosample id"),
        strip_prefix="BIOSAMPLE:",
    )
    secondary_samples = _all_col_tokens(
        ci_row,
        ("secondary sample accession", "secondary_sample_accession", "ena_sample"),
    )
    runs = _all_col_tokens(
        ci_row,
        ("accession", "sra accession", "sra_accession", "run", "run_accession"),
    )
    assemblies = _all_col_tokens(
        ci_row,
        ("assembly accession", "assembly_accession", "genbank assembly", "gca"),
    )
    bioprojects = _all_col_tokens(
        ci_row,
        ("project id", "bioproject", "bioproject_accession"),
        strip_prefix="BIOPROJECT:",
    )
    secondary_studies = _all_col_tokens(
        ci_row,
        ("study accession", "study_accession", "secondary_study_accession"),
    )
    genome_id = (
        _first_col(ci_row, ("uberstrain", "assembly barcode", "barcode", "genome_id"))
        or (biosamples[0] if biosamples else None)
        or (runs[0] if runs else None)
        or _first_col(ci_row, ("name",))
    )
    if genome_id is None:
        raise ValueError("enterobase row has no usable identifier (Uberstrain, Barcode, Sample ID, or Accession)")

    species = _normalize_species(_first_col(ci_row, ("species", "organism", "scientific_name")))
    strain = _collapse_ws(_first_col(ci_row, ("name", "strain", "isolate")))
    subject_id = _collapse_ws(_first_col(ci_row, ("subject id", "subject_id", "patient id", "patient_id")))
    fasta_path = _first_col(ci_row, ("fasta_path", "assembly_path"))
    pmids: list[str] = []
    for col in ("pmid", "pubmed_id", "publication", "reference"):
        pmids.extend(_clean_pmids(ci_row.get(col)))

    synonyms: list[tuple[str, str, str]] = []
    for b in biosamples:
        for s in secondary_samples:
            if b != s:
                synonyms.append(("sample", b, s))
    anchors = biosamples or secondary_samples
    for r in runs:
        for a in anchors:
            if r != a:
                synonyms.append(("run_to_sample", r, a))
    for bp in bioprojects:
        for st in secondary_studies:
            if bp != st:
                synonyms.append(("study", bp, st))

    rec = {
        "genome_id": genome_id,
        "portal": "enterobase",
        "species": species,
        "biosample_accessions": biosamples + [s for s in secondary_samples if s not in biosamples],
        "assembly_accessions": assemblies,
        "sra_accessions": runs,
        "bioproject_accessions": bioprojects + [s for s in secondary_studies if s not in bioprojects],
        "pmids": pmids,
        "strain": strain,
        "subject_id": subject_id,
        "fasta_path": fasta_path,
    }
    return rec, synonyms


def _extract_pubmlst(ci_row: Mapping[str, object]) -> tuple[dict[str, object], list[tuple[str, str, str]]]:
    biosamples = _all_col_tokens(
        ci_row,
        ("biosample_accession", "biosample", "biosample_id"),
        strip_prefix="BIOSAMPLE:",
    )
    ena_samples = _all_col_tokens(
        ci_row,
        ("ena_sample", "secondary_sample_accession", "sample_accession"),
    )
    runs = _all_col_tokens(
        ci_row,
        ("ena_run_accession", "ena_run", "sra_accession", "run_accession"),
    )
    genbank_asm = _all_col_tokens(
        ci_row,
        ("assembly_accession", "genbank_accession", "genbank_assembly"),
    )
    refseq_asm = _all_col_tokens(
        ci_row,
        ("refseq_accession", "refseq_assembly"),
    )
    bioprojects = _all_col_tokens(
        ci_row,
        ("bioproject_accession", "bioproject"),
        strip_prefix="BIOPROJECT:",
    )
    ena_studies = _all_col_tokens(
        ci_row,
        ("ena_study", "study_accession", "secondary_study_accession"),
    )
    pub_id = _first_col(ci_row, ("pubmlst_id", "id", "genome_id"))
    isolate_name = _first_col(ci_row, ("isolate",))
    if pub_id is not None and not pub_id.upper().startswith("PUBMLST:") and pub_id.isdigit():
        genome_id = f"pubmlst:{pub_id}"
    else:
        genome_id = (
            pub_id
            or isolate_name
            or (biosamples[0] if biosamples else None)
            or (runs[0] if runs else None)
        )
    if genome_id is None:
        raise ValueError("pubmlst row has no usable identifier (id, isolate, BioSample, or run)")

    species = _normalize_species(_first_col(ci_row, ("species", "organism", "scientific_name")))
    strain = _collapse_ws(_first_col(ci_row, ("strain", "isolate", "aliases")))
    subject_id = _collapse_ws(_first_col(ci_row, ("patient_id", "subject_id", "host_id")))
    fasta_path = _first_col(ci_row, ("fasta_path", "assembly_path"))
    pmids: list[str] = []
    for col in ("pubmed_id", "pmid", "references", "publication"):
        pmids.extend(_clean_pmids(ci_row.get(col)))

    synonyms: list[tuple[str, str, str]] = []
    for b in biosamples:
        for e in ena_samples:
            if b != e:
                synonyms.append(("sample", b, e))
    anchors = biosamples or ena_samples
    for r in runs:
        for a in anchors:
            if r != a:
                synonyms.append(("run_to_sample", r, a))
    for g in genbank_asm:
        for rf in refseq_asm:
            if g != rf:
                synonyms.append(("assembly", g, rf))
    for bp in bioprojects:
        for st in ena_studies:
            if bp != st:
                synonyms.append(("study", bp, st))

    assemblies = genbank_asm + [a for a in refseq_asm if a not in genbank_asm]
    rec = {
        "genome_id": genome_id,
        "portal": "pubmlst",
        "species": species,
        "biosample_accessions": biosamples + [e for e in ena_samples if e not in biosamples],
        "assembly_accessions": assemblies,
        "sra_accessions": runs,
        "bioproject_accessions": bioprojects + [s for s in ena_studies if s not in bioprojects],
        "pmids": pmids,
        "strain": strain,
        "subject_id": subject_id,
        "fasta_path": fasta_path,
    }
    return rec, synonyms


EXTRACTORS = {
    "ebi_amr": _extract_ebi_amr,
    "cryptic": _extract_cryptic,
    "ncbi_pd": _extract_ncbi_pd,
    "enterobase": _extract_enterobase,
    "pubmlst": _extract_pubmlst,
}


def normalize_rows(
    portal: str,
    rows: Iterable[Mapping[str, object]],
) -> tuple[list[dict[str, object]], list[tuple[str, str, str]]]:
    """Normalizes raw portal rows into deduplicated `grm dedup` records and
    sorted synonym triples `(kind, accession_a, accession_b)`.

    Repeated rows for the same `genome_id` (common in per-antibiotic AST tables)
    are merged deterministically.
    """
    portal_name = canonical_portal(portal)
    extract = EXTRACTORS[portal_name]
    merged: dict[str, dict[str, object]] = {}
    synonyms_set: set[tuple[str, str, str]] = set()

    for idx, raw_row in enumerate(rows, 1):
        ci_row = _ci_lookup(raw_row)
        try:
            rec, syns = extract(ci_row)
        except ValueError as exc:
            raise ValueError(f"{portal_name} row {idx}: {exc}") from exc
        for kind, a, b in syns:
            # Canonicalize undirected pair order for sample/study/assembly; keep
            # run_to_sample oriented as (run, sample).
            if kind in ("sample", "study", "assembly") and b < a:
                a, b = b, a
            synonyms_set.add((kind, a, b))

        gid = str(rec["genome_id"])
        existing = merged.get(gid)
        if existing is None:
            merged[gid] = {
                "genome_id": gid,
                "portal": portal_name,
                "species": rec["species"],
                "biosample_accessions": set(rec["biosample_accessions"]),
                "assembly_accessions": set(rec["assembly_accessions"]),
                "sra_accessions": set(rec["sra_accessions"]),
                "bioproject_accessions": set(rec["bioproject_accessions"]),
                "pmids": set(rec["pmids"]),
                "strain": rec["strain"],
                "subject_id": rec["subject_id"],
                "fasta_path": rec["fasta_path"],
            }
        else:
            for field in (
                "biosample_accessions",
                "assembly_accessions",
                "sra_accessions",
                "bioproject_accessions",
                "pmids",
            ):
                existing[field].update(rec[field])  # type: ignore[union-attr]
            for scalar in ("species", "strain", "fasta_path"):
                if existing[scalar] is None and rec[scalar] is not None:
                    existing[scalar] = rec[scalar]
            if rec["subject_id"] is not None:
                if existing["subject_id"] is None:
                    existing["subject_id"] = rec["subject_id"]
                elif existing["subject_id"] != rec["subject_id"]:
                    raise ValueError(
                        f"{portal_name}: genome_id {gid!r} carries conflicting subject_id "
                        f"values ({existing['subject_id']!r} vs {rec['subject_id']!r})"
                    )

    records_out: list[dict[str, object]] = []
    for gid in sorted(merged, key=genome_id_sort_key):
        m = merged[gid]
        out: dict[str, object] = {
            "genome_id": gid,
            "portal": m["portal"],
        }
        if m["species"] is not None:
            out["species"] = m["species"]
        bios = sorted(m["biosample_accessions"])  # type: ignore[arg-type]
        if bios:
            out["biosample_accession"] = bios[0]
            if len(bios) > 1:
                out["biosample_aliases"] = bios[1:]
        asms = sorted(m["assembly_accessions"])  # type: ignore[arg-type]
        if asms:
            out["assembly_accession"] = ";".join(asms)
        sras = sorted(m["sra_accessions"])  # type: ignore[arg-type]
        if sras:
            out["sra_accession"] = ";".join(sras)
        bps = sorted(m["bioproject_accessions"])  # type: ignore[arg-type]
        if bps:
            out["bioproject_accession"] = bps[0]
            if len(bps) > 1:
                out["bioproject_aliases"] = bps[1:]
        pmids = sorted(m["pmids"])  # type: ignore[arg-type]
        if pmids:
            out["pmid"] = ";".join(pmids)
        if m["strain"] is not None:
            out["strain"] = m["strain"]
        if m["subject_id"] is not None:
            out["subject_id"] = m["subject_id"]
        if m["fasta_path"] is not None:
            out["fasta_path"] = m["fasta_path"]
        records_out.append(out)

    return records_out, sorted(synonyms_set)


def write_query_jsonl(path: pathlib.Path, records: Iterable[Mapping[str, object]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        for rec in records:
            fh.write(json.dumps(rec, sort_keys=True, ensure_ascii=True) + "\n")
            count += 1
    return count


def write_synonym_tsv(path: pathlib.Path, pairs: Iterable[tuple[str, str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    deduped = sorted(set(pairs))
    with path.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(SYNONYM_HEADER) + "\n")
        for kind, a, b in deduped:
            fh.write(f"{kind}\t{a}\t{b}\n")
    return len(deduped)


def merge_synonym_tsvs(paths: Iterable[pathlib.Path], out_path: pathlib.Path) -> int:
    pairs: set[tuple[str, str, str]] = set()
    for p in paths:
        with _open_text(p) as fh:
            for line_no, raw in enumerate(fh, 1):
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                cols = [c.strip() for c in line.split("\t")]
                if line_no == 1 and cols and cols[0].lower() in ("kind", "domain", "type"):
                    continue
                if len(cols) < 3:
                    raise ValueError(f"{p}:{line_no}: expected 3 tab-separated columns")
                kind, a, b = cols[0].lower(), cols[1].upper(), cols[2].upper()
                if kind in ("sample", "study", "assembly") and b < a:
                    a, b = b, a
                pairs.add((kind, a, b))
    return write_synonym_tsv(out_path, pairs)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Normalize public AMR portal metadata into grm dedup query JSONL and INSDC synonym TSV.",
    )
    parser.add_argument(
        "--portal",
        choices=sorted(PORTAL_ALIASES),
        help="Source AMR portal format (ebi_amr/cabbage, cryptic, ncbi_pd/narms, enterobase, pubmlst).",
    )
    parser.add_argument(
        "--input",
        type=pathlib.Path,
        action="append",
        default=[],
        help="Input metadata table (.tsv, .csv, .jsonl, plain or .gz). Repeatable.",
    )
    parser.add_argument(
        "--delimiter",
        default=None,
        help="Optional explicit delimiter override (default: inferred from extension or header).",
    )
    parser.add_argument(
        "--out",
        type=pathlib.Path,
        help="Output canonical query JSONL path for `grm dedup screen --query`.",
    )
    parser.add_argument(
        "--synonym-out",
        type=pathlib.Path,
        default=None,
        help="Optional output INSDC/ENA/NCBI synonym TSV path for `grm dedup screen --synonym-map`.",
    )
    parser.add_argument(
        "--merge-synonyms",
        type=pathlib.Path,
        nargs="+",
        default=None,
        help="Merge multiple synonym TSV files into `--synonym-out` without portal conversion.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.merge_synonyms is not None:
        if args.synonym_out is None:
            parser.error("--merge-synonyms requires --synonym-out")
        n = merge_synonym_tsvs(args.merge_synonyms, args.synonym_out)
        print(f"merged {n} synonym pairs -> {args.synonym_out}", file=sys.stderr)
        return 0
    if not args.portal or not args.input or args.out is None:
        parser.error("--portal, --input, and --out are required unless using --merge-synonyms")

    all_rows: list[dict[str, object]] = []
    for path in args.input:
        all_rows.extend(read_table_rows(path, delimiter=args.delimiter))
    records, synonyms = normalize_rows(args.portal, all_rows)
    n_rec = write_query_jsonl(args.out, records)
    msg = f"{canonical_portal(args.portal)}: {n_rec} isolates -> {args.out}"
    if args.synonym_out is not None:
        n_syn = write_synonym_tsv(args.synonym_out, synonyms)
        msg += f"; {n_syn} synonym pairs -> {args.synonym_out}"
    print(msg, file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
