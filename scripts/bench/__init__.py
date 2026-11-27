"""Unified published AMR tool benchmark harness (Track A, rev-14-safe).

Modules:
- vocab: canonical antibiotic vocabulary parsed directly from src/drugs.rs.
- protocol: two-mode adapter contract, 4-tier support registry, time-split
  provenance, external-cohort evaluation_only guard, and schema validators.
- pretrained: deduplicated per-genome FASTA execution, tool output parsers,
  canonical genome_calls.jsonl, and projection onto rev14, grm-fold-file-v1,
  and external deduplicated cohorts.
- trainable: per-(dataset, split, outer_fold) execution consuming grm export
  or grm-fold-file-v1 with strict fold parity (scripts/c7-outside/plans.py)
  and strict inside-outer-training-fold fitting.
- score: unified grm-bench-result-v1 scoring core with exact Clopper-Pearson
  and Wilson 95% CIs, worst/best bounds under no_call + failure, pairwise
  common-support (min 25 R / 25 S) + McNemar, full-denominator coverage,
  4,000-draw SHA-256 species/lineage bootstraps, and scan.py zero-ID hygiene.
- run_lambda: Lambda execution runner with pinned run-manifest.json, offline
  mode enforcement, atomic DONE/FAILED sentinels, and separate logs/timing.tsv.
"""

__all__ = [
    "vocab",
    "protocol",
    "pretrained",
    "trainable",
    "score",
    "run_lambda",
]
