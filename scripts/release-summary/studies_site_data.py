#!/usr/bin/env python3
"""Build web/src/data/additional-studies.json from the five committed comparison records.

Reads only aggregate TSVs under:
  - reference/time-split/
  - reference/kover-comparison/
  - reference/dense-network/
  - reference/amrfinderplus-comparison/
  - reference/tbprofiler-comparison/

Emits deterministic compact JSON (LF terminated) with zero genome IDs,
zero sequences, and zero per-isolate rows. Pass `--check FILE` to verify
an existing bundled JSON matches byte for byte.
"""
from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")


def _read_kv(path: pathlib.Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    key_col = "name" if "name" in rows[0] else "key"
    return {r[key_col]: r["value"] for r in rows}


def _read_tsv(path: pathlib.Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def build(root: pathlib.Path) -> dict:
    ref = root / "reference"

    # 1. Time split (first-public-day split)
    ts_sum = _read_kv(ref / "time-split" / "summary.tsv")
    ts_read = _read_tsv(ref / "time-split" / "readings.tsv")
    ts_by_learner: dict[str, dict] = {}
    for l in LEARNERS:
        entry: dict[str, dict] = {}
        for m in ("recall", "log_loss", "specificity", "balanced_accuracy"):
            matches = [r for r in ts_read if r["learner"] == l and r["contrast"] == "total" and r["measure"] == m]
            if len(matches) != 1:
                raise ValueError(f"expected 1 time-split row for {l}/{m}, got {len(matches)}")
            row = matches[0]
            entry[m] = {
                "median": float(row["median"]),
                "low": float(row["low"]),
                "high": float(row["high"]),
                "reading": row["reading"],
                "worse": int(row["worse"]),
                "better": int(row["better"]),
            }
        ts_by_learner[l] = entry
    time_split = {
        "datasets": int(ts_sum["datasets"]),
        "species": int(ts_sum["species"]),
        "tested_pairs": int(ts_sum["tested_pairs"]),
        "resistant_pairs": int(ts_sum["resistant_pairs"]),
        "susceptible_pairs": int(ts_sum["susceptible_pairs"]),
        "answer": ts_sum["answer"],
        "learners_with_cost": int(ts_sum["learners_with_cost"]),
        "total_log_loss_worse": int(ts_sum["total_log_loss_worse"]),
        "total_balanced_accuracy_worse": int(ts_sum["total_balanced_accuracy_worse"]),
        "total_recall_worse": int(ts_sum["total_recall_worse"]),
        "total_specificity_worse": int(ts_sum["total_specificity_worse"]),
        "by_learner": ts_by_learner,
    }

    # 2. Kover 2.0.7 vs GRM scm
    kv_tot = _read_kv(ref / "kover-comparison" / "totals.tsv")
    kv_sum = {r["measure"]: r for r in _read_tsv(ref / "kover-comparison" / "summary.tsv")}
    kover = {
        "datasets": int(kv_tot["datasets"]),
        "species": int(kv_tot["species"]),
        "genomes": int(kv_tot["genomes"]),
        "folds": int(kv_tot["folds"]),
        "kover_one_rule_folds": int(kv_tot["kover_one_rule_folds"]),
        "balanced_accuracy": {
            "value": float(kv_sum["balanced_accuracy"]["value"]),
            "low": float(kv_sum["balanced_accuracy"]["low"]),
            "high": float(kv_sum["balanced_accuracy"]["high"]),
            "higher": int(kv_sum["balanced_accuracy"]["datasets_higher"]),
            "equal": int(kv_sum["balanced_accuracy"]["datasets_equal"]),
            "lower": int(kv_sum["balanced_accuracy"]["datasets_lower"]),
        },
        "recall": {
            "value": float(kv_sum["recall"]["value"]),
            "low": float(kv_sum["recall"]["low"]),
            "high": float(kv_sum["recall"]["high"]),
        },
        "specificity": {
            "value": float(kv_sum["specificity"]["value"]),
            "low": float(kv_sum["specificity"]["low"]),
            "high": float(kv_sum["specificity"]["high"]),
        },
        "log_loss": {
            "value": float(kv_sum["log_loss"]["value"]),
            "low": float(kv_sum["log_loss"]["low"]),
            "high": float(kv_sum["log_loss"]["high"]),
        },
        "ties_equal_ba": {
            "value": float(kv_sum["ties_equal"]["value"]),
            "low": float(kv_sum["ties_equal"]["low"]),
            "high": float(kv_sum["ties_equal"]["high"]),
        },
    }

    # 3. Dense neural network
    dn_rows = _read_tsv(ref / "dense-network" / "summary.tsv")
    dn_rule = _read_tsv(ref / "dense-network" / "rule.tsv")
    dn_outcome = next(r["value"] for r in dn_rule if r["rule"] == "all")
    dn_lookup = {(r["split"], r["comparator"], r["measure"]): r for r in dn_rows}

    def _dn_cell(comp: str, measure: str) -> dict:
        r = dn_lookup[("lineage", comp, measure)]
        return {
            "diff": float(r["difference_median"]),
            "low": float(r["low"]),
            "high": float(r["high"]),
            "positive": int(r["positive"]),
            "negative": int(r["negative"]),
        }

    best5_ba = dn_lookup[("lineage", "best5", "balanced_accuracy")]
    dense_network = {
        "datasets": int(best5_ba["datasets"]),
        "species": int(best5_ba["species"]),
        "outcome": dn_outcome,
        "vs_best5": {
            "balanced_accuracy": _dn_cell("best5", "balanced_accuracy"),
            "recall": _dn_cell("best5", "recall"),
            "log_loss": _dn_cell("best5", "log_loss"),
        },
        "by_learner": {
            l: {
                "balanced_accuracy": _dn_cell(l, "balanced_accuracy"),
                "recall": _dn_cell(l, "recall"),
                "log_loss": _dn_cell(l, "log_loss"),
            }
            for l in LEARNERS
        },
    }

    # 4. AMRFinderPlus 4.2.7 program comparison
    am_rows = _read_tsv(ref / "amrfinderplus-comparison" / "summary.tsv")
    am_figs = {r["key"]: r["value"] for r in _read_tsv(ref / "amrfinderplus-comparison" / "figures.tsv")}
    am_lookup = {(r["set"], r["contrast"], r["measure"]): r for r in am_rows}
    p_vs_rule = am_lookup[("primary-new-235", "program_minus_known_gene_rule", "ba")]

    def _am_learner(l: str) -> dict:
        ba = am_lookup[("primary-new-235", f"learner_minus_program:lineage:{l}", "ba")]
        sens = am_lookup[("primary-new-235", f"learner_minus_program:lineage:{l}", "sens")]
        spec = am_lookup[("primary-new-235", f"learner_minus_program:lineage:{l}", "spec")]
        return {
            "ba_median": float(ba["median"]),
            "ba_lo": float(ba["median_species_lo"]),
            "ba_hi": float(ba["median_species_hi"]),
            "ba_ahead": int(ba["ahead"]),
            "ba_level": int(ba["level"]),
            "ba_behind": int(ba["behind"]),
            "sens_median": float(sens["median"]),
            "spec_median": float(spec["median"]),
        }

    amrfinderplus = {
        "datasets": int(p_vs_rule["datasets"]),
        "species": int(p_vs_rule["species"]),
        "program_vs_known_gene_rule": {
            "ahead": int(p_vs_rule["ahead"]),
            "level": int(p_vs_rule["level"]),
            "behind": int(p_vs_rule["behind"]),
            "median": float(p_vs_rule["median"]),
            "mean": float(p_vs_rule["mean"]),
            "ahead_by_10pts": int(am_figs["program_ahead_by_0.10_or_more"]),
            "ahead_by_10pts_quinolone": int(am_figs["of_them_quinolone"]),
        },
        "by_learner_lineage": {l: _am_learner(l) for l in LEARNERS},
    }

    # 5. TBProfiler 6.6.3 on M. tuberculosis assemblies
    tb_sum = _read_kv(ref / "tbprofiler-comparison" / "summary.tsv")
    tb_med = {(r["way"], r["set"], r["caller"], r["measure"]): r for r in _read_tsv(ref / "tbprofiler-comparison" / "medians.tsv")}
    tb_cmp = {(r["way"], r["set"], r["caller"]): r for r in _read_tsv(ref / "tbprofiler-comparison" / "comparisons.tsv")}
    tbprofiler = {
        "by_name_endpoints": int(tb_sum["by_name_endpoints"]),
        "by_name_rows": int(tb_sum["by_name_rows"]),
        "by_name_rows_with_call": int(tb_sum["by_name_rows_with_call"]),
        "by_name_rows_without_call": int(tb_sum["by_name_rows_without_call"]),
        "graded_endpoints": int(tb_sum["graded_endpoints"]),
        "graded_rows": int(tb_sum["graded_rows"]),
        "graded_rows_with_call": int(tb_sum["graded_rows_with_call"]),
        "graded_rows_without_call": int(tb_sum["graded_rows_without_call"]),
        "by_name_sensitivity_median": float(tb_med[("who_assoc", "by_name", "tool", "sensitivity")]["median"]),
        "by_name_specificity_median": float(tb_med[("who_assoc", "by_name", "tool", "specificity")]["median"]),
        "graded_sensitivity_median": float(tb_med[("who_assoc", "graded", "tool", "sensitivity")]["median"]),
        "graded_specificity_median": float(tb_med[("who_assoc", "graded", "tool", "specificity")]["median"]),
        "graded_vs_learners_lineage": {
            l: {
                "tool_higher": int(tb_cmp[("who_assoc", "graded", f"{l}:lineage")]["tool_higher"]),
                "tool_lower": int(tb_cmp[("who_assoc", "graded", f"{l}:lineage")]["tool_lower"]),
                "equal": int(tb_cmp[("who_assoc", "graded", f"{l}:lineage")]["equal"]),
                "median_learner_minus_tool": float(tb_cmp[("who_assoc", "graded", f"{l}:lineage")]["median_minus_tool"]),
            }
            for l in LEARNERS
        },
    }

    return {
        "source": "reference/{time-split,kover-comparison,dense-network,amrfinderplus-comparison,tbprofiler-comparison}",
        "evaluation_revision": 14,
        "time_split": time_split,
        "kover": kover,
        "dense_network": dense_network,
        "amrfinderplus": amrfinderplus,
        "tbprofiler": tbprofiler,
    }


def render(root: pathlib.Path) -> str:
    return json.dumps(build(root), separators=(",", ":"), ensure_ascii=False) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=pathlib.Path, default=pathlib.Path(__file__).resolve().parents[2])
    parser.add_argument("--check", type=pathlib.Path, help="Fail unless FILE matches generated bytes.")
    args = parser.parse_args(argv)
    out = render(args.root)
    if args.check is not None:
        actual = args.check.read_text(encoding="utf-8")
        if actual != out:
            print(f"{args.check} differs from generated output", file=sys.stderr)
            return 1
        return 0
    sys.stdout.write(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
