# Outer folds from a file

`grm learn` draws two outer partitions for each dataset: a stratified random
split and a split that holds out whole lineages. Some questions need a
partition neither draw can make, such as "train on genomes that were public by
a given day and test on genomes deposited later". `--fold-file` takes that
partition from a file.

```
grm --data <root> learn --species <slug> --fold-file folds.json --fold-out <dir>
```

It scores the one dataset the file names and writes
`<dir>/<dataset>.<scheme>.json`. It writes no report and reads none: a
dataset's report under the data root is the same file before and after.
Without `--fold-file`, `grm learn` is unchanged.

## The file

```json
{
  "format": "grm-fold-file-v1",
  "dataset": "klebsiella_pneumoniae__meropenem",
  "scheme": "first-public-past",
  "test": [["573.1001", "573.1002"]],
  "train": ["573.2001", "573.2002", "573.2003"],
  "unused": ["573.3001"]
}
```

| Key | Meaning |
|---|---|
| `format` | `grm-fold-file-v1` |
| `dataset` | The dataset id, as in `curated/datasets.json` |
| `scheme` | The partition's name in the result: 1 to 48 lowercase letters, digits and hyphens. `random` and `lineage` are refused, because they name the splits `grm learn` draws |
| `test` | One list of genome ids per outer fold. Fold f is scored on `test[f]` |
| `train` | Genomes that train in every fold and are never scored |
| `unused` | Genomes in no fold |

Fold f trains on `train` and on every other fold's test genomes. A file with
one test fold is a plain train and test split. A file with five test folds and
empty `train` and `unused` is a five-fold cross-validation.

Every genome of the dataset's cohort must be listed exactly once. The cohort
is the set of genomes `grm learn` evaluates for the dataset: those with a
reported call that are in the species matrix. Listing the ones left out, under
`unused`, is what lets the reader tell a deliberate omission from a file made
for another cohort or cut short.

## What is supplied and what is not

Only the outer partition comes from the file. Inside each fold the evaluator
does what it does for a split it drew:

- the fold's candidate features are chosen from its training genomes alone;
- scm and cart tune on inner folds of the training genomes, boosting and the
  lasso on a holdout drawn from them;
- tuning keeps each group of k-mer-identical genomes whole, as the random
  split does (docs/duplicate-isolation.md);
- every draw is seeded from `seed` and the fold number.

So a file that repeats the random split a report drew gives that report's
random-split learner rows. Two tests hold this
(`roles_that_repeat_a_drawn_split_give_that_splits_scores`,
`a_fold_file_that_repeats_the_random_split_gives_the_reports_random_rows`).

Genomes under `unused` reach no fit. The support floor, the fewest training
genomes a candidate feature must be present in and absent from, is the
larger of 3 and `n / 200` rounded down, with `n` the genomes the file tests
or trains; a report takes `n` from its whole cohort. A group of
k-mer-identical genomes keeps only its tested and training members before
tuning reads it, so a group of unused genomes does not move tuning onto
grouped draws. The candidates, the floor and the tuning draws therefore come
from the tested and training genomes alone, and listing more genomes under
`unused` changes no score. Two tests hold this
(`genomes_a_fold_file_leaves_unused_set_neither_the_support_floor_nor_the_tuning_groups`,
`a_copy_pair_among_unused_samples_cannot_change_a_folds_fits`).

## What is refused

The command stops before any fit, and writes nothing, when:

- the file has a key it does not know, lacks one, or is not `grm-fold-file-v1`;
- the dataset is not one of the species' datasets;
- a genome is listed twice, is not in the cohort, or a cohort genome is
  missing;
- a test fold is empty, or a fold's test genomes or training genomes lack a
  class;
- a fold tests a genome and trains on a k-mer-identical copy of it. The
  refusal names both genomes. Put the copy in the same test fold or under
  `unused`;
- an inner fold or a tuning holdout would lack a class;
- the result file exists and `--force` was not given.

## The result

`<dataset>.<scheme>.json` holds:

| Key | Meaning |
|---|---|
| `format` | `grm-fold-file-result-v1` |
| `evaluation_revision` | The evaluator's revision |
| `fold_file_sha256` | SHA-256 of the fold file's bytes, as `sha256sum` prints it |
| `input_fingerprint` | The fingerprint the dataset's report records for the same inputs and settings |
| `eval_config` | The settings of the fits |
| `id`, `species_slug`, `antibiotic`, `scheme` | The dataset and the partition's name |
| `min_support` | The per-fit support floor: the larger of 3 and `n / 200` rounded down, with `n` the genomes the file tests or trains. It is the report's floor when nothing is unused |
| `counts` | Cohort genomes, resistant genomes, and how many are tested, train only, or are unused |
| `folds` | Per outer fold: training and test genomes, the resistant ones among each, and candidate features |
| `samples` | The scored genomes in order: id, reported class, outer fold |
| `split` | The learner rows, in the shape of a report's split. `test_fold` and each learner's `scores` follow `samples`, and `n_groups` counts the lineages of the tested and training genomes |
| `seconds` | Wall-clock time of the run |

The hash is of the file as it is on disk, so the same partition written with
other spacing is another hash.

## What the result is not

It is not a dataset report, and `grm report` does not read it. It holds no
lineage split, catalog baseline, MIC regression or final model, and it scores
only the genomes the file tests, so its numbers are not comparable with a
report's unless the file tests every genome.

A partition by date or by study inside one database is a retrospective test
on that database. It is not a prospective test and not external validation.

Neither `EVALUATION_REVISION` nor the report schema changes: no report gains
or loses a field, and the default path fits each fold through the same code as
before.
