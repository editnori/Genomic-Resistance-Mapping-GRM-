"""Pretrained AMR tool execution, parsers, genome_calls.jsonl, and split projection.

Pretrained adapters run ONCE per unique genome FASTA (not once per genome-drug
row), emit canonical genome_calls.jsonl records, and project those genome-level
calls onto:
- GRM rev14 random and lineage splits (from grm export dataset.json or compact reports),
- grm-fold-file-v1 splits (study-whole, time, or other fold-file schemes),
- external deduplicated cohorts (selection_role='evaluation_only').

Unsupported drugs, no_call, I (intermediate), and runtime failures are NEVER
mapped to Susceptible ('S').
"""

from __future__ import annotations

import csv
import io
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from scripts.bench.protocol import (
    AdapterSpec,
    BenchError,
    MODE_PRETRAINED,
    TIER_UNSUPPORTED,
    load_fold_file,
    require,
    sha256_file,
    validate_external_cohort,
    validate_genome_call_record,
)
from scripts.bench.vocab import DRUGS_RS, load_vocabulary

# Multi-drug proxy call combination precedence:
# 1. Any explicit R -> R (resistant determinant takes precedence over incomplete coverage/I).
# 2. Else any failure -> failure.
# 3. Else any no_call -> no_call.
# 4. Else any I -> I (never mapped to S).
# 5. Else all S -> S.
_CALL_PRECEDENCE = {"R": 4, "failure": 3, "no_call": 2, "I": 1, "S": 0}


def combine_proxy_calls(calls: list[tuple[str, str | None]]) -> tuple[str, str | None]:
    """Deterministically combine calls across multiple mapped tool_drugs (e.g., proxy rules).

    Each item is (call, reason). Returns (combined_call, combined_reason).
    Never maps failure, no_call, or I to S.
    """
    require(len(calls) >= 1, "combine_proxy_calls requires at least one call")
    best_call, best_reason = max(calls, key=lambda item: _CALL_PRECEDENCE[item[0]])
    if best_call in ("failure", "no_call", "I") and not best_reason:
        best_reason = f"proxy_{best_call}"
    return best_call, best_reason


# ------------------------------------------------------------------ genome_calls.jsonl I/O


def write_genome_calls_jsonl(path: Path | str, records: list[dict[str, Any]]) -> int:
    """Validate and write canonical genome_calls.jsonl records in deterministic sorted order."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    validated = [validate_genome_call_record(dict(r)) for r in records]
    validated.sort(key=lambda r: (r["species_slug"], r["genome_id"], r["reading"], r["tool_drug"]))
    # Refuse duplicate (species_slug, genome_id, reading, tool_drug) keys
    seen: set[tuple[str, str, str, str]] = set()
    lines: list[str] = []
    for r in validated:
        key = (r["species_slug"], r["genome_id"], r["reading"], r["tool_drug"])
        require(key not in seen, f"duplicate genome_calls.jsonl key: {key}")
        seen.add(key)
        lines.append(json.dumps(r, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False))
    tmp = p.with_name(p.name + ".tmp")
    tmp.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    tmp.replace(p)
    return len(validated)


def read_genome_calls_jsonl(path: Path | str) -> list[dict[str, Any]]:
    """Read and strictly validate every record in a genome_calls.jsonl file."""
    p = Path(path)
    records: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str, str]] = set()
    with p.open("r", encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.strip()
            if not line:
                continue
            try:
                rec = json.loads(line)
            except json.JSONDecodeError as exc:
                raise BenchError(f"{p}:{lineno}: invalid JSON: {exc}") from exc
            validate_genome_call_record(rec)
            key = (rec["species_slug"], rec["genome_id"], rec["reading"], rec["tool_drug"])
            require(key not in seen, f"{p}:{lineno}: duplicate call record for {key}")
            seen.add(key)
            records.append(rec)
    return records


# ------------------------------------------------------------------ unique-FASTA execution


@dataclass(frozen=True)
class GenomeTask:
    genome_id: str
    species_slug: str
    fasta_path: Path


def deduplicate_genomes(tasks: list[GenomeTask]) -> list[GenomeTask]:
    """Deduplicate genome-drug pairs down to one run per unique (species_slug, genome_id).

    Verifies that a repeated genome_id never points to conflicting species or FASTA paths.
    """
    by_id: dict[str, GenomeTask] = {}
    for t in tasks:
        if t.genome_id in by_id:
            prev = by_id[t.genome_id]
            require(
                prev.species_slug == t.species_slug and prev.fasta_path.resolve() == t.fasta_path.resolve(),
                f"conflicting FASTA/species for genome {t.genome_id}",
            )
            continue
        by_id[t.genome_id] = t
    return sorted(by_id.values(), key=lambda t: (t.species_slug, t.genome_id))


def run_pretrained_unique_genomes(
    adapter: AdapterSpec,
    tasks: list[GenomeTask],
    runner: Callable[[AdapterSpec, GenomeTask], list[dict[str, Any]]],
    out_jsonl: Path | str,
    *,
    reading: str | None = None,
) -> dict[str, Any]:
    """Run a pretrained adapter ONCE per unique genome FASTA and write genome_calls.jsonl.

    If `runner` raises an exception on a genome, records explicit `failure` records
    for all tool drugs relevant to that species (never silent S!).
    """
    require(adapter.mode == MODE_PRETRAINED, f"{adapter.tool}: expected pretrained adapter, got {adapter.mode!r}")
    use_reading = reading or adapter.primary_reading
    require(use_reading in adapter.manifest["readings"], f"{adapter.tool}: unknown reading {use_reading!r}")

    unique_tasks = deduplicate_genomes(tasks)
    all_records: list[dict[str, Any]] = []
    fasta_hashes: dict[str, str] = {}
    failed_genomes = 0

    # Tool drugs relevant for each species (for explicit failure accounting).
    species_tool_drugs: dict[str, list[str]] = {}
    for (sp, _cd), entry in adapter.species_support.items():
        if entry.support_tier != TIER_UNSUPPORTED:
            species_tool_drugs.setdefault(sp, []).extend(entry.tool_drugs)
    species_tool_drugs = {sp: sorted(set(tds)) for sp, tds in species_tool_drugs.items()}

    for task in unique_tasks:
        require(task.fasta_path.is_file(), f"missing FASTA for genome {task.genome_id}: {task.fasta_path}")
        fasta_hashes[task.genome_id] = sha256_file(task.fasta_path)
        try:
            recs = runner(adapter, task)
            for rec in recs:
                validate_genome_call_record(rec)
                require(rec["genome_id"] == task.genome_id, "runner emitted mismatched genome_id")
                require(rec["species_slug"] == task.species_slug, "runner emitted mismatched species_slug")
            all_records.extend(recs)
        except Exception as exc:  # noqa: BLE001
            failed_genomes += 1
            reason = f"runtime_failure:{type(exc).__name__}"
            tds = species_tool_drugs.get(task.species_slug) or ["__genome__"]
            for td in tds:
                all_records.append({
                    "genome_id": task.genome_id,
                    "species_slug": task.species_slug,
                    "reading": use_reading,
                    "tool_drug": td,
                    "call": "failure",
                    "reason": reason,
                })

    n_written = write_genome_calls_jsonl(out_jsonl, all_records)
    return {
        "tool": adapter.tool,
        "reading": use_reading,
        "input_tasks": len(tasks),
        "unique_genomes_run": len(unique_tasks),
        "failed_genomes": failed_genomes,
        "records_written": n_written,
        "unique_fasta_sha256_count": len(set(fasta_hashes.values())),
    }


# ------------------------------------------------------------------ 7 pretrained output parsers


def parse_amrfinderplus_tsv(
    tsv_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "tool_default",
) -> list[dict[str, Any]]:
    """Parse AMRFinderPlus TSV output into canonical genome_calls records."""
    require(
        reading in ("tool_default", "tool_default_intact", "genes_core", "genes_all_scopes"),
        f"unknown amrfinderplus reading: {reading!r}",
    )
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    hit_tokens: set[str] = set()
    for row in reader:
        etype = (row.get("Element type") or row.get("Type") or "").strip().upper()
        if etype != "AMR":
            continue
        scope = (row.get("Scope") or "").strip().lower()
        method = (row.get("Method") or "").strip().upper()
        subtype = (row.get("Element subtype") or row.get("Subtype") or "").strip().upper()

        if reading == "tool_default" and scope != "core":
            continue
        if reading == "tool_default_intact":
            if scope != "core" or "PARTIAL" in method or "INTERNAL_STOP" in method:
                continue
        if reading == "genes_core":
            if scope != "core" or subtype == "POINT" or method == "POINTX":
                continue
        if reading == "genes_all_scopes":
            if scope not in ("core", "plus") or subtype == "POINT" or method == "POINTX":
                continue

        for col in ("Subclass", "Class"):
            val = (row.get(col) or "").strip()
            if val:
                for tok in val.split("/"):
                    t = tok.strip().upper()
                    if t:
                        hit_tokens.add(t)
                hit_tokens.add(val.upper())

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        call = "R" if td.upper() in hit_tokens else "S"
        out.append({
            "genome_id": genome_id,
            "species_slug": species_slug,
            "reading": reading,
            "tool_drug": td,
            "call": call,
        })
    return out


def parse_resfinder_json(
    json_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "pheno_default",
) -> list[dict[str, Any]]:
    """Parse ResFinder/PointFinder v4 JSON output into canonical genome_calls records."""
    require(reading in ("pheno_default", "acquired_only"), f"unknown resfinder reading: {reading!r}")
    doc = json.loads(json_text)
    phenos = doc.get("phenotypes") or {}
    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        entry = phenos.get(td) or phenos.get(td.lower())
        if entry is None:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": "drug_absent_in_resfinder_output",
            })
            continue
        if reading == "acquired_only":
            # Only count acquired seq_regions, not point_mutations
            seq_regs = entry.get("amr_species_relevant", True)
            has_acq = bool(entry.get("seq_regions") or entry.get("acquired_resistant"))
            is_r = bool(seq_regs and has_acq)
        else:
            is_r = bool(entry.get("amr_resistant"))
        out.append({
            "genome_id": genome_id,
            "species_slug": species_slug,
            "reading": reading,
            "tool_drug": td,
            "call": "R" if is_r else "S",
        })
    return out


def parse_card_rgi_tsv(
    tsv_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "strict_perfect",
) -> list[dict[str, Any]]:
    """Parse CARD RGI main TSV output into canonical genome_calls records."""
    require(
        reading in ("strict_perfect", "strict_perfect_no_efflux", "perfect_only"),
        f"unknown card_rgi reading: {reading!r}",
    )
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    allowed_cutoffs = {"Perfect"} if reading == "perfect_only" else {"Perfect", "Strict"}
    hit_tokens: set[str] = set()
    for row in reader:
        cutoff = (row.get("Cut_Off") or "").strip()
        if cutoff not in allowed_cutoffs:
            continue
        mech = (row.get("Resistance Mechanism") or "").strip().lower()
        if reading == "strict_perfect_no_efflux" and mech == "antibiotic efflux":
            continue
        for col in ("Antibiotic", "Drug Class"):
            val = (row.get(col) or "").strip()
            if val:
                for tok in val.split(";"):
                    t = tok.strip().lower()
                    if t:
                        hit_tokens.add(t)

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        call = "R" if td.lower() in hit_tokens else "S"
        out.append({
            "genome_id": genome_id,
            "species_slug": species_slug,
            "reading": reading,
            "tool_drug": td,
            "call": call,
        })
    return out


def parse_staramr_summary_tsv(
    tsv_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "tool_default",
) -> list[dict[str, Any]]:
    """Parse StarAMR summary.tsv output into canonical genome_calls records."""
    require(reading == "tool_default", f"unknown staramr reading: {reading!r}")
    rows = list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))
    require(len(rows) == 1, f"staramr summary.tsv must contain 1 row for {genome_id}, got {len(rows)}")
    row = rows[0]
    qc = (row.get("Quality Module") or "Passed").strip().lower()
    pheno_raw = (row.get("Predicted Phenotype") or "").strip()
    predicted: set[str] = set()
    if pheno_raw and pheno_raw.lower() != "sensitive":
        for tok in pheno_raw.split(","):
            t = tok.strip().lower()
            if t and t != "sensitive":
                predicted.add(t)

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        if td.lower() in predicted:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "R",
            })
        elif qc == "failed":
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": "staramr_quality_module_failed",
            })
        else:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "S",
            })
    return out


def parse_kleborate_tsv(
    tsv_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "tool_default",
) -> list[dict[str, Any]]:
    """Parse Kleborate TSV output into canonical genome_calls records.

    tool_drug names correspond to Kleborate resistance column tokens (e.g.
    Bla_Carb_acquired, Bla_ESBL_acquired, Flq_mutations, Col_mutations, etc.).
    """
    require(reading in ("tool_default", "acquired_only"), f"unknown kleborate reading: {reading!r}")
    rows = list(csv.DictReader(io.StringIO(tsv_text), delimiter="\t"))
    require(len(rows) == 1, f"kleborate TSV must contain 1 row for {genome_id}, got {len(rows)}")
    row = rows[0]

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        if td not in row:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": "kleborate_column_absent",
            })
            continue
        if reading == "acquired_only" and "mutations" in td.lower():
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "S",
            })
            continue
        val = (row[td] or "").strip()
        is_r = bool(val and val != "-" and val.lower() != "none")
        out.append({
            "genome_id": genome_id,
            "species_slug": species_slug,
            "reading": reading,
            "tool_drug": td,
            "call": "R" if is_r else "S",
        })
    return out


def parse_tbprofiler_json(
    json_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "who_assoc",
) -> list[dict[str, Any]]:
    """Parse TBProfiler JSON output into canonical genome_calls records.

    Matches scripts/tbprofiler-comparison/record.py:
    - who_assoc: Assoc w R or Assoc w R - Interim
    - who_assoc_only: Assoc w R alone
    - tool_default: any reported dr_variant for the drug
    - who_assoc_small: who_assoc excluding large deletions / >=50 bp indels
    - Uncovered target loci produce no_call unless a resistant call stands.
    """
    require(
        reading in ("who_assoc", "who_assoc_only", "tool_default", "who_assoc_small"),
        f"unknown tbprofiler reading: {reading!r}",
    )
    doc = json.loads(json_text)
    r_drugs: set[str] = set()
    for var in doc.get("dr_variants") or []:
        var_type = (var.get("type") or "").lower()
        is_large = "large_deletion" in var_type or int(var.get("indel_len") or 0) >= 50
        for d_entry in var.get("drugs") or []:
            drug = (d_entry.get("drug") or "").strip().lower()
            conf = (d_entry.get("confidence") or "").strip()
            if not drug:
                continue
            if reading == "tool_default":
                r_drugs.add(drug)
            elif reading == "who_assoc":
                if conf in ("Assoc w R", "Assoc w R - Interim"):
                    r_drugs.add(drug)
            elif reading == "who_assoc_only":
                if conf == "Assoc w R":
                    r_drugs.add(drug)
            elif reading == "who_assoc_small":
                if conf in ("Assoc w R", "Assoc w R - Interim") and not is_large:
                    r_drugs.add(drug)

    uncovered_drugs: set[str] = set()
    qc = doc.get("qc") or {}
    for tqc in qc.get("target_qc") or []:
        cov = float(tqc.get("percent_depth_pass", 100.0))
        if cov < 100.0:
            for d in tqc.get("drugs") or []:
                uncovered_drugs.add(str(d).strip().lower())

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        tdl = td.lower()
        if tdl in r_drugs:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "R",
            })
        elif tdl in uncovered_drugs:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": "incomplete_target_coverage",
            })
        else:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "S",
            })
    return out


def parse_mykrobe_json(
    json_text: str,
    genome_id: str,
    species_slug: str,
    expected_tool_drugs: list[str],
    reading: str = "tool_default",
) -> list[dict[str, Any]]:
    """Parse Mykrobe JSON output into canonical genome_calls records.

    Mykrobe emits 'R', 'S', 'r' (minor/intermediate -> canonical 'I', NEVER S),
    and 'N'/'U'/'F' (-> canonical 'no_call').
    """
    require(reading in ("tool_default", "strict_r_only"), f"unknown mykrobe reading: {reading!r}")
    doc = json.loads(json_text)
    if genome_id in doc and isinstance(doc[genome_id], dict):
        sample_doc = doc[genome_id]
    elif len(doc) == 1:
        sample_doc = next(iter(doc.values()))
    else:
        raise BenchError(f"mykrobe JSON does not contain sample entry for {genome_id}")
    susc = sample_doc.get("susceptibility") or {}

    out: list[dict[str, Any]] = []
    for td in sorted(set(expected_tool_drugs)):
        entry = susc.get(td)
        if entry is None:
            out.append({
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": "drug_absent_in_mykrobe_output",
            })
            continue
        pred = str(entry.get("predict") or "").strip()
        if pred == "R":
            rec: dict[str, Any] = {
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "R",
            }
        elif pred == "S":
            rec = {
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "S",
            }
        elif pred == "r":
            # Minor-allele / intermediate call: NEVER mapped to S!
            rec = {
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "I",
                "reason": "mykrobe_minor_allele_r",
            }
        elif pred in ("N", "U", "F", ""):
            rec = {
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "no_call",
                "reason": f"mykrobe_uncalled_{pred or 'empty'}",
            }
        else:
            rec = {
                "genome_id": genome_id,
                "species_slug": species_slug,
                "reading": reading,
                "tool_drug": td,
                "call": "failure",
                "reason": f"mykrobe_unknown_token_{pred}",
            }
        out.append(rec)
    return out


PARSERS = {
    "amrfinderplus": parse_amrfinderplus_tsv,
    "resfinder_pointfinder": parse_resfinder_json,
    "card_rgi": parse_card_rgi_tsv,
    "staramr": parse_staramr_summary_tsv,
    "kleborate": parse_kleborate_tsv,
    "tbprofiler": parse_tbprofiler_json,
    "mykrobe": parse_mykrobe_json,
}


# ------------------------------------------------------------------ split projection


@dataclass(frozen=True)
class ProjectedEndpoint:
    """In-memory projected per-genome calls for one (dataset, split) before aggregation.

    Never serialized directly into grm-bench-result-v1 (contains genome_ids/lineages).
    """

    dataset: str
    species_slug: str
    canonical_drug: str
    split_scheme: str
    selection_role: str
    support_tier: str
    unsupported_reason: str
    reading: str
    genome_ids: tuple[str, ...]
    truth: tuple[bool, ...]
    lineages: tuple[int | str | None, ...]
    calls: tuple[str, ...]  # 'R' | 'S' | 'I' | 'no_call' | 'failure' | 'unsupported'
    fold_of: tuple[int, ...]  # outer test fold index (0 for single-split external cohorts)
    fold_file_sha256: str | None = None
    comparator_calls: dict[str, tuple[str, ...]] | None = None


def _index_genome_calls(
    records: list[dict[str, Any]],
    reading: str,
) -> dict[tuple[str, str, str], dict[str, Any]]:
    index: dict[tuple[str, str, str], dict[str, Any]] = {}
    for r in records:
        if r["reading"] != reading:
            continue
        index[(r["species_slug"], r["genome_id"], r["tool_drug"])] = r
    return index


def _project_calls_for_genomes(
    adapter: AdapterSpec,
    species_slug: str,
    canonical_drug: str,
    genome_ids: list[str],
    call_index: dict[tuple[str, str, str], dict[str, Any]],
    reading: str,
) -> tuple[str, str, tuple[str, ...]]:
    """Resolve calls for each genome_id using adapter's pre-committed support tier."""
    support = adapter.support_for(species_slug, canonical_drug)
    if support.support_tier == TIER_UNSUPPORTED:
        return (
            TIER_UNSUPPORTED,
            support.unsupported_reason,
            tuple("unsupported" for _ in genome_ids),
        )

    projected: list[str] = []
    for gid in genome_ids:
        sub_calls: list[tuple[str, str | None]] = []
        for td in support.tool_drugs:
            rec = call_index.get((species_slug, gid, td))
            if rec is None:
                # Check genome-level failure sentinel if present
                gfail = call_index.get((species_slug, gid, "__genome__"))
                if gfail is not None and gfail["call"] == "failure":
                    sub_calls.append(("failure", gfail.get("reason") or "genome_runtime_failure"))
                else:
                    # Missing call record for a supported drug -> failure, NEVER S!
                    sub_calls.append(("failure", f"missing_genome_call:{td}"))
            else:
                sub_calls.append((rec["call"], rec.get("reason")))
        combined, _ = combine_proxy_calls(sub_calls)
        projected.append(combined)

    return support.support_tier, "-", tuple(projected)


def project_onto_rev14_dataset(
    adapter: AdapterSpec,
    records: list[dict[str, Any]],
    dataset_doc: dict[str, Any],
    split_scheme: str,
    *,
    reading: str | None = None,
    comparator_scores: dict[str, list[float]] | None = None,
    score_cutoff: float = 0.5,
) -> ProjectedEndpoint:
    """Project pretrained genome_calls onto a GRM rev14 dataset (random or lineage split).

    Accepts either a `grm export` dataset.json or a rev14 compact report dict.
    """
    require(adapter.mode == MODE_PRETRAINED, f"{adapter.tool}: expected pretrained adapter")
    require(split_scheme in ("random", "lineage"), f"rev14 drawn split must be random or lineage, got {split_scheme!r}")
    use_reading = reading or adapter.primary_reading
    require(use_reading in adapter.manifest["readings"], f"{adapter.tool}: unknown reading {use_reading!r}")

    vocab = load_vocabulary(DRUGS_RS)
    ds_id = dataset_doc["id"]
    sp = dataset_doc["species_slug"]
    ab = dataset_doc["antibiotic"]
    require(vocab.is_canonical(ab), f"{ds_id}: antibiotic {ab!r} is not canonical in src/drugs.rs")

    if "samples" in dataset_doc:
        # grm export dataset.json
        require(
            dataset_doc.get("evaluation_revision") == 14,
            f"{ds_id}: evaluation_revision {dataset_doc.get('evaluation_revision')} != 14",
        )
        gids = list(dataset_doc["samples"])
        truth = [bool(v) for v in dataset_doc["resistant"]]
        lineages: list[int | str | None] = list(dataset_doc["lineage"])
        split_entry = next(
            (s for s in dataset_doc.get("splits", []) if s["plan"]["scheme"] == split_scheme),
            None,
        )
        require(split_entry is not None, f"{ds_id}: split {split_scheme!r} not present in dataset.json")
        require(
            not split_entry["plan"].get("unavailable_reason"),
            f"{ds_id}: split {split_scheme!r} is unavailable: {split_entry['plan'].get('unavailable_reason')}",
        )
        fold_of = [int(f) for f in split_entry["plan"]["outer"]]
    elif "sid" in dataset_doc:
        # rev14 compact report
        gids = list(dataset_doc["sid"])
        truth = [bool(v) for v in dataset_doc["y"]]
        lineages = list(dataset_doc["lineage"])
        split_obj = (dataset_doc.get("splits") or {}).get(split_scheme)
        require(split_obj is not None, f"{ds_id}: split {split_scheme!r} not present in compact report")
        fold_of = list(split_obj.get("test_fold", [0] * len(gids)))
    else:
        raise BenchError(f"{ds_id}: unrecognized rev14 dataset format")

    call_index = _index_genome_calls(records, use_reading)
    tier, unsup_reason, projected = _project_calls_for_genomes(adapter, sp, ab, gids, call_index, use_reading)

    comp_calls: dict[str, tuple[str, ...]] | None = None
    if comparator_scores is not None:
        comp_calls = {}
        for comp_name, scs in sorted(comparator_scores.items()):
            require(len(scs) == len(gids), f"{ds_id}: comparator {comp_name} length mismatch")
            comp_calls[comp_name] = tuple("R" if float(s) >= score_cutoff else "S" for s in scs)

    return ProjectedEndpoint(
        dataset=ds_id,
        species_slug=sp,
        canonical_drug=ab,
        split_scheme=split_scheme,
        selection_role="internal_cv",
        support_tier=tier,
        unsupported_reason=unsup_reason,
        reading=use_reading,
        genome_ids=tuple(gids),
        truth=tuple(truth),
        lineages=tuple(lineages),
        calls=projected,
        fold_of=tuple(fold_of),
        comparator_calls=comp_calls,
    )


def project_onto_fold_file(
    adapter: AdapterSpec,
    records: list[dict[str, Any]],
    dataset_doc: dict[str, Any],
    fold_file_path: Path | str,
    *,
    reading: str | None = None,
    comparator_calls_by_genome: dict[str, dict[str, str]] | None = None,
) -> ProjectedEndpoint:
    """Project pretrained genome_calls onto a grm-fold-file-v1 split (study-whole, time, etc.).

    Only genomes assigned to `test` folds are scored; `train` and `unused` genomes
    are excluded from test evaluation.
    """
    require(adapter.mode == MODE_PRETRAINED, f"{adapter.tool}: expected pretrained adapter")
    use_reading = reading or adapter.primary_reading
    require(use_reading in adapter.manifest["readings"], f"{adapter.tool}: unknown reading {use_reading!r}")

    vocab = load_vocabulary(DRUGS_RS)
    ds_id = dataset_doc["id"]
    sp = dataset_doc["species_slug"]
    ab = dataset_doc["antibiotic"]
    require(vocab.is_canonical(ab), f"{ds_id}: antibiotic {ab!r} is not canonical in src/drugs.rs")

    all_gids = list(dataset_doc["samples"] if "samples" in dataset_doc else dataset_doc["sid"])
    all_truth = [bool(v) for v in (dataset_doc["resistant"] if "resistant" in dataset_doc else dataset_doc["y"])]
    all_lineages = list(dataset_doc["lineage"])
    idx_of = {g: i for i, g in enumerate(all_gids)}

    ff = load_fold_file(fold_file_path, expected_dataset=ds_id, sample_ids=all_gids)

    test_gids: list[str] = []
    test_truth: list[bool] = []
    test_lineages: list[int | str | None] = []
    test_folds: list[int] = []
    for f_idx, fold_ids in enumerate(ff.test):
        for g in fold_ids:
            pos = idx_of[g]
            test_gids.append(g)
            test_truth.append(all_truth[pos])
            test_lineages.append(all_lineages[pos])
            test_folds.append(f_idx)

    call_index = _index_genome_calls(records, use_reading)
    tier, unsup_reason, projected = _project_calls_for_genomes(adapter, sp, ab, test_gids, call_index, use_reading)

    comp_calls: dict[str, tuple[str, ...]] | None = None
    if comparator_calls_by_genome is not None:
        comp_calls = {}
        for comp_name, mapping in sorted(comparator_calls_by_genome.items()):
            vals: list[str] = []
            for g in test_gids:
                require(g in mapping, f"{ds_id}: comparator {comp_name} missing test genome")
                c = mapping[g]
                require(c in ("R", "S", "I", "no_call", "failure"), f"invalid comparator call {c!r}")
                vals.append(c)
            comp_calls[comp_name] = tuple(vals)

    return ProjectedEndpoint(
        dataset=ds_id,
        species_slug=sp,
        canonical_drug=ab,
        split_scheme=ff.scheme,
        selection_role="internal_cv",
        support_tier=tier,
        unsupported_reason=unsup_reason,
        reading=use_reading,
        genome_ids=tuple(test_gids),
        truth=tuple(test_truth),
        lineages=tuple(test_lineages),
        calls=projected,
        fold_of=tuple(test_folds),
        fold_file_sha256=ff.sha256,
        comparator_calls=comp_calls,
    )


def project_onto_external_cohort(
    adapter: AdapterSpec,
    records: list[dict[str, Any]],
    cohort_doc: dict[str, Any],
    *,
    reading: str | None = None,
    comparator_calls_by_genome: dict[str, dict[str, str]] | None = None,
) -> ProjectedEndpoint:
    """Project pretrained genome_calls onto an external deduplicated validation cohort.

    Strictly enforces selection_role == 'evaluation_only'.
    """
    require(adapter.mode == MODE_PRETRAINED, f"{adapter.tool}: expected pretrained adapter")
    use_reading = reading or adapter.primary_reading
    require(use_reading in adapter.manifest["readings"], f"{adapter.tool}: unknown reading {use_reading!r}")

    validate_external_cohort(cohort_doc, allow_selection=False)
    ds_id = cohort_doc["dataset"]
    sp = cohort_doc["species_slug"]
    ab = cohort_doc["antibiotic"]
    gids = list(cohort_doc["samples"])
    truth = [bool(v) for v in cohort_doc["resistant"]]
    lineages: list[int | str | None] = (
        list(cohort_doc["lineage"]) if cohort_doc.get("lineage") is not None else [None] * len(gids)
    )

    call_index = _index_genome_calls(records, use_reading)
    tier, unsup_reason, projected = _project_calls_for_genomes(adapter, sp, ab, gids, call_index, use_reading)

    comp_calls: dict[str, tuple[str, ...]] | None = None
    if comparator_calls_by_genome is not None:
        comp_calls = {}
        for comp_name, mapping in sorted(comparator_calls_by_genome.items()):
            vals: list[str] = []
            for g in gids:
                require(g in mapping, f"{ds_id}: comparator {comp_name} missing external genome")
                c = mapping[g]
                require(c in ("R", "S", "I", "no_call", "failure"), f"invalid comparator call {c!r}")
                vals.append(c)
            comp_calls[comp_name] = tuple(vals)

    return ProjectedEndpoint(
        dataset=ds_id,
        species_slug=sp,
        canonical_drug=ab,
        split_scheme=f"external:{cohort_doc['cohort_name']}",
        selection_role="evaluation_only",
        support_tier=tier,
        unsupported_reason=unsup_reason,
        reading=use_reading,
        genome_ids=tuple(gids),
        truth=tuple(truth),
        lineages=tuple(lineages),
        calls=projected,
        fold_of=tuple(0 for _ in gids),
        comparator_calls=comp_calls,
    )
