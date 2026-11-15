# Dataset export for a learner outside GRM

`grm export` writes what a learner outside GRM needs to be scored on the same genomes and outer folds as GRM's five learners, without refitting them: a dataset's projected matrix, its random and lineage fold plans, and the random forest's screened rows for each outer fold. Only the forest reads that screen. Boosting and the lasso screen again on the fitting side of their tuning holdouts, and the set covering machine and the tree read every candidate row. The export fits nothing and changes no report.

```
grm --data <root> export --species <slug> --out <dir> [--only <antibiotic or id>,...] [--screened-features <n>] [--force]
```

The species must be counted and built first (`grm count`, `grm structure`, `grm matrix`), as for `grm learn`. Without a matrix the export stops and names those three commands. Like `grm learn`, it also stops when `genomes.txt` adds, removes or reorders genomes after counting (`species::ensure_selection_counted`). It reads the same inputs `grm learn` reads (`analysis::arm_inputs`) at the same settings (`EvalConfig::default()`, with `--screened-features`, 4,000 by default).

`--only` takes antibiotics or dataset ids. A name that matches no dataset of the species stops the export before it writes anything, and the message lists the names that matched nothing beside the species' datasets.

## Writing

Each dataset goes in its own directory, `<dir>/<dataset id>/`. The export never writes into a dataset directory that exists. When the directory of any dataset it would export exists (an earlier export, an unfinished one, or a file at that path), it stops before writing anything for any dataset and names each such directory. With `--force` it replaces `matrix.bin` and `dataset.json` in those directories and says so for each one. Other files there stay.

Both files are written under a temporary name (`matrix.bin.tmp`, `dataset.json.tmp`) and then renamed. An earlier `dataset.json` is removed first, and the new one is renamed into place last. So a directory without `dataset.json` is an unfinished export: do not read its `matrix.bin`. Every write error names its file.

A dataset whose cohort falls short of its class minimum is not exported, just as `grm learn` writes no report for it. With `--force`, a directory such a dataset already has is left as it was, and the export says so.

## Files

`matrix.bin` is the dataset's projected matrix (`Dataset::project`): `rows` rows of `words` 64-bit little-endian words. Bit `c % 64` of word `c / 64` of row `r` is set when feature `r` is present in sample `c`.

`dataset.json` (`format` `grm-dataset-export-v1`) holds:

| Field | Meaning |
|---|---|
| `evaluation_revision` | The evaluation revision of the code that wrote the export (`analysis::EVALUATION_REVISION`). |
| `input_fingerprint` | The fingerprint `grm learn` records as the report's `input_fingerprint` for these inputs and settings. It covers the matrix, cohort, labels, lineages, catalog calls and settings, and no code. An export pairs with a report only when the fingerprints are equal and `evaluation_revision` equals the report's `evaluation.revision` too: code of another revision can project, plan or screen the same inputs differently. |
| `screened_features`, `folds`, `inner_folds`, `seed`, `min_support` | The settings, and the dataset's support floor. |
| `samples`, `resistant`, `lineage`, `duplicates` | Sample ids in matrix column order, reported calls, lineage ids and k-mer-identical duplicate groups. |
| `rows`, `cols`, `words`, `matrix_sha256` | The shape of `matrix.bin` and its SHA-256. |
| `splits[].plan` | The random and lineage `FoldPlan` (`ablate::plan::fold_plan`, repeat 0 at the evaluation seed): the outer folds `evaluate` scores, each outer fold's inner folds and both tuning holdouts, and `sha256`. An unavailable plan keeps its outer folds and its `unavailable_reason`. |
| `splits[].screens` | One entry per outer fold of an available plan, with `fold` and three more fields. `n_train` counts the fold's training samples. `n_candidates` counts its candidate rows (`learn::training_candidates` on the training samples). `screened_rows` lists the forest's screened rows (`learn::screen` on the training samples) as `matrix.bin` row indices, best first. |

### The screen, and how to recompute it

The training samples and the candidate rows are counts in the file, not lists. For outer fold `f` of a plan, they follow from the plan and the matrix:

1. The training samples are the samples `i` with `plan.outer[i] != f`. There are `n_train` of them.
2. A row of `matrix.bin` is a candidate when its support on the training samples lies in [`min_support`, `n_train` - `min_support`]. Of rows that have the same pattern on the training samples, only the lowest row is a candidate. There are `n_candidates` of them.
3. Each candidate is scored by the chi-square of its table of presence by reported call on the training samples, 0 when a margin of the table is empty. `screened_rows` holds the `screened_features` highest scores, highest first.

Among rows with an equal score, the lower matrix row comes first. The cut can fall inside a tie: when rows tie at the last place kept, the lower rows are kept and the others are left out. This is common. In the review of 15 Nov 2026, rows tied at the cut on 27 of 30 recomputed outer folds of three release datasets. A reader who orders ties another way gets another set of rows.

## Limits

- **Unavailable plans.** A plan whose `unavailable_reason` is set keeps its outer folds but has no inner folds, holdouts or screens, and its report has no `test_fold` for that split. Do not score its outer folds: GRM scored no learner on them, so a number made on them has nothing to compare with.
- **Size and memory.** `matrix.bin` takes `rows` x `words` x 8 bytes. The largest release dataset, `escherichia_coli__ampicillin` (3,078,562 rows x 1,965 genomes), writes 763.5 MB, and its export peaked at 4,962,952 kB of resident memory at 4 threads. All 274 release datasets would write 24.87 GB. The export holds the species matrix, the projected matrix and a byte copy of it in memory at once.
- **A fold plan is not a fold file.** `splits[].plan.outer` gives one fold number per sample, by position in `samples`, and the plan also carries each outer fold's inner folds and tuning holdouts. A lineage plan turned into a list of genome ids keeps its outer partition but loses its lineage-grouped tuning: inner folds and holdouts drawn again from such a list are not the plan's.
- **One fold draw.** `grm export` takes no seed flag, so its plans are the release's draw, the one `grm learn` scores without `--seed`. A report that `grm learn --seed` wrote under `draws/` was scored on other folds and has another `input_fingerprint`, so it does not pair with an export ([fold-draws.md](fold-draws.md)).

## What it does not do

It does not reproduce boosting's or the lasso's screens, which run on the fitting side of their tuning holdouts; a learner that needs those can recompute them from the plan's holdouts with the same chi-square rule. It does not write the training samples or the candidate rows as lists (see above). It does not export the species matrix, k-mers or annotation.
