# Study split record

This directory holds the record of the two experiments with folds by study at revision 14. The page that reads it is [`docs/study-split.md`](../../docs/study-split.md). `scripts/study-split/record.py` writes every file here from the reduced input set, and its `--check` mode fails on any byte it would not write.

## Files

| File | Rows | Contents |
|---|---|---|
| `checks.tsv` | 20 | The script's checks with the count that passed each, and the count of the inputs' genome ids found in the record |
| `contrasts.tsv` | 1,281 | Every contrast: median over datasets, exact species interval, reading, datasets worse, equal and better |
| `datasets.tsv` | 70 | One row per dataset: species, antibiotic, genomes, resistant genomes, BioProjects, genomes with none, the largest BioProject and its share, and whether the release scored a lineage split |
| `distances.tsv` | 2,100 | Genomes in another fold between two draws of a kind, after pairing fold labels for the largest agreement |
| `draws.tsv` | 100 | Pairs with a cost, or reading worse, per draw |
| `eligibility.tsv` | 274 | The release's datasets against the first plan's three counts, whether each entered, and the fold check's reason where it stopped one |
| `first.tsv` | 69 | The first experiment's readings, with the leave-one-species-out readings that keep each |
| `fits.tsv` | 5,250 | One row per fit and learner: confusion counts, recall, specificity, balanced accuracy, log loss, error rate, and a null fit flag |
| `foldcut.tsv` | 30 | The first experiment per study fold, in groups by training size and class shift |
| `folds.tsv` | 1,750 | The first experiment per study fold and learner: test genomes, training counts, and both arms' counts and log loss on the fold's genomes |
| `leaveout.tsv` | 160 | The second experiment with each species left out, and the two largest species together |
| `lineage.tsv` | 345 | The release's lineage split per dataset and learner, in the form of fits.tsv |
| `mates.tsv` | 18 | Test genomes with a study mate or lineage-cluster mate in training, per kind and draw |
| `nostudy.tsv` | 135 | The three contrasts on genomes without a BioProject, and on genomes with one in all the datasets and in the same datasets |
| `order.tsv` | 22 | Times from the hash lists, the run logs, the file times and the second plan's text |
| `parts.tsv` | 15 | The mean drop as the part matched by fold shape and the part not matched by fold shape |
| `planfiles.tsv` | 2 | Each plan's file in the input set, its sha256, the files its hash list holds, and its hash time |
| `plans.tsv` | 1,050 | One row per fold plan: kind, draw, placement rule, fold file sha256, the five fold sizes and resistant counts, and test genomes with a study or lineage-cluster mate in training |
| `rules.tsv` | 75 | The second experiment's pairs with a cost under the plan's rule and under other rules |
| `settings.tsv` | 23 | The settings the script applies, the first plan's resample count, and the two design keys the fold builders hashed |
| `shapes.tsv` | 25 | Fold sizes and class balance per arm |
| `subsets.tsv` | 5,850 | Per fit and learner, the genomes with and without a BioProject apart, on datasets with at least 10 genomes without one |
| `summary.tsv` | 111 | Single figures the page quotes |
| `thresholds.tsv` | 29 | The learners each answer holds for, at every threshold |
| `totals.tsv` | 2 | Distinct genomes and BioProjects over the datasets |
| `variants.tsv` | 36 | The first plan's eligibility variants, with the datasets, species and pairs each admits by the three counts |
| `variation.tsv` | 560 | The plan builder's table: genomes moved from draw 0 under placement rules A and B, per dataset and draw |

## Rules

**Measures.** The script computes every measure from scores and labels. A score at or above 0.5 calls resistant. Log loss clips the score given to the true class to [0.0001, 1 - 0.0001], so a wrong 0 and a wrong 1 lose the same amount. Recall, specificity, balanced accuracy and the error rate are exact fractions of the confusion counts.

**Loss.** The first plan's rule reads log loss for every learner. The second plan's rule reads the error rate for scm, which scores only 0 or 1, and log loss for the others.

**scm.** A 0 or 1 score makes scm's clipped log loss a right call's loss plus 9.21 times its error rate. Every contrast, interval and count takes scm's log loss in this form, which keeps the error rate's ties. `fits.tsv` gives the log loss summed from the scores, and `summary.tsv` gives the largest gap between the two.

**Intervals.** Each species interval is exact: it takes every multiset of the species drawn with replacement, at the weight of its multinomial count. The interval holds species resampling only.

**Derived.** The script's derive() writes `contrasts.tsv`, `first.tsv`, `leaveout.tsv`, `rules.tsv`, `thresholds.tsv`, `draws.tsv`, `parts.tsv`, `mates.tsv`, `shapes.tsv`, `nostudy.tsv`, `foldcut.tsv`, `summary.tsv` from the row files.

## Guards

`--check` rebuilds every file from the reduced input set and fails on any byte that differs. Without the archive, `tests/study_split` guards each file in one or more of these ways.

| Guard | Files |
|---|---|
| Rebuilt | `contrasts.tsv`, `first.tsv`, `leaveout.tsv`, `rules.tsv`, `thresholds.tsv`, `draws.tsv`, `parts.tsv`, `mates.tsv`, `shapes.tsv`, `nostudy.tsv`, `foldcut.tsv`, `summary.tsv`: the suite rebuilds each from the row files, byte for byte |
| Agreed | `datasets.tsv`, `fits.tsv`, `subsets.tsv`, `totals.tsv`: the suite checks their cells against each other and the derived files |
| Quoted | `checks.tsv`, `contrasts.tsv`, `datasets.tsv`, `draws.tsv`, `eligibility.tsv`, `first.tsv`, `foldcut.tsv`, `leaveout.tsv`, `mates.tsv`, `nostudy.tsv`, `order.tsv`, `parts.tsv`, `planfiles.tsv`, `plans.tsv`, `rules.tsv`, `settings.tsv`, `shapes.tsv`, `summary.tsv`, `thresholds.tsv`: the suite rebuilds the page's lines from them |
| Pinned | `checks.tsv`, `distances.tsv`, `eligibility.tsv`, `folds.tsv`, `inputs.sha256`, `lineage.tsv`, `order.tsv`, `planfiles.tsv`, `plans.tsv`, `settings.tsv`, `variants.tsv`, `variation.tsv`: the suite holds each one's sha256 |
| Whole | `README.md` is what the script writes, and `SHA256SUMS` gives each file's sha256 |

## Plans

The two plans stay in the reduced input set as `first/PLAN.md` and `second/PLAN.md`. `planfiles.tsv` gives each one's sha256, the files its hash list holds and its hash time, and `--check` hashes the plans again. `settings.tsv` gives `first_design_key` and `second_design_key` as the fold builders hashed them into every keyed order; they stay as they were so that the folds rebuild.

## Commands

```sh
python3 -I scripts/study-split/record.py --inputs <dir> --out reference/study-split
python3 -I scripts/study-split/record.py --inputs <dir> --check reference/study-split
python3 -I -m unittest discover -s tests/study_split -t .
```

`<dir>` is the reduced input set, the directory that `study-split-inputs.tar.zst` unpacks to. `inputs.sha256` lists the sha256 of each of its files.
