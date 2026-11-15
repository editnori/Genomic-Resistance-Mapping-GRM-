# Dataset export for a learner outside GRM

`grm export` writes what a learner outside GRM needs to be scored on the same genomes, folds and screened columns as the five learners, without refitting them. It fits nothing and changes no report.

```
grm --data <root> export --species <slug> [--only <antibiotic or id>,...] --out <dir>
```

The species must be counted and built first (`grm count`, `grm structure`, `grm matrix`), as for `grm learn`. The export reads the same inputs `grm learn` reads (`analysis::arm_inputs`) at the same settings (`EvalConfig::default()`, with `--screened-features`, 4,000 by default).

## Files, per dataset under `<dir>/<dataset id>/`

`matrix.bin` is the dataset's projected matrix (`Dataset::project`): `rows` rows of `words` 64-bit little-endian words. Bit `c % 64` of word `c / 64` of row `r` is set when feature `r` is present in sample `c`.

`dataset.json` (`format` `grm-dataset-export-v1`) holds:

| Field | Meaning |
|---|---|
| `input_fingerprint` | The fingerprint `grm learn` records as the report's `input_fingerprint` for these inputs and settings. Equal fingerprints mean the export and the report were built from the same matrix, cohort, labels, lineages, catalog calls and settings. |
| `samples`, `resistant`, `lineage`, `duplicates` | Sample ids in matrix column order, reported calls, lineage ids and k-mer-identical duplicate groups. |
| `rows`, `cols`, `words`, `matrix_sha256` | The shape of `matrix.bin` and its SHA-256. |
| `splits[].plan` | The random and lineage `FoldPlan` (`ablate::plan::fold_plan`, repeat 0 at the evaluation seed): the outer folds `evaluate` scores, each outer fold's inner folds and both tuning holdouts, and `sha256`. An unavailable plan keeps its outer folds and its `unavailable_reason`. |
| `splits[].screens` | Per outer fold of an available plan: training samples, candidate rows (`learn::training_candidates` on the fold's training samples) and the screened rows, as `matrix.bin` row indices, best first. This is the forest's screen (`learn::screen` on the training samples). |

## What it does not do

It does not reproduce boosting's or the lasso's screens, which run on the fitting side of their tuning holdouts; a learner that needs those can recompute them from the plan's holdouts with the same chi-square rule. It does not export the species matrix, k-mers or annotation.
