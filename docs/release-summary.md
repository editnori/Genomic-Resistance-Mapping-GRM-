# Release summary: the main numbers of a GRM release

Four scripts in `scripts/release-summary/` compute the main numbers of a GRM release from two of its public files: `index.json` and the directory `intervals/`. A fifth script puts two releases side by side. Nothing is fitted and no genome is read. The outputs for revisions 12 and 14 are committed under `reference/release-summary/`. The README there names each release asset and the sha256 of each input and output.

These are retrospective analyses with internal cross-validation on one public database. Each release has 274 datasets in 23 species. A dataset is one species and one drug.

## The inputs

`index.json` has one row per dataset. Each row gives each learner's balanced accuracy (the mean of sensitivity and specificity) on the random split and on the lineage split. It gives the learner's very major error (VME: the share of resistant genomes called susceptible) and its major error (ME: the share of susceptible genomes called resistant). The row also names the best learner and gives the known-gene rule's figures.

`intervals/<dataset>.json` gives the counts behind VME and ME (`x` errors out of `n` genomes) for each learner on each split and for the known-gene rule. Each rate has an exact 95% interval (Clopper-Pearson).

The random split assigns genomes to folds at random. The lineage split holds out whole lineages, so it is the stricter of the two. The five learners are the set covering machine (`scm`), the decision tree (`cart`), the random forest (`rf`), boosted trees (`gbt`) and the lasso (`lasso`).

The best learner of a dataset is the one of the five with the highest random-split balanced accuracy. The known-gene rule calls a genome resistant from the genes a catalog lists. It is fitted on nothing, so its one score answers for both splits.

Each script reads the release's revision from `evaluation_revision` in the index. Each understands revisions 12 and 14 and stops with a message on any other.

## The outputs

**`uncertainty.json`** comes from `python3 -I scripts/release-summary/uncertainty.py INDEX.json INTERVALS_DIR [--reps N]`. It has these parts:

- **checks.** Every intervals file has the index's revision and names its dataset, and its counts reproduce the index's VME and ME.
- **medians.** The median balanced accuracy of each learner on each split, over the datasets where the learner has a score.
- **gaps.** For each learner, the median over datasets of the random-split score minus the lineage-split score. The difference is paired: both scores come from the same dataset.
- **paired_lineage.** The median paired difference on the lineage split between the set covering machine and each other learner, with the number of datasets where the set covering machine is at least as high.
- **selection.** The best learner's lineage score against each fixed learner, and against the highest lineage score any of the five reaches on the same dataset. That highest score is an oracle: no one can know it before scoring.
- **shares** and **catalog.** The share of datasets where the best learner is ahead of the known-gene rule, and the median difference between them, on each split.
- **error_bars.** The error-limit counts described below.

Each median and share has two bootstrap 95% intervals, which are percentile intervals over resamples drawn with replacement. `ci95_datasets` resamples datasets as if they were independent. `ci95_species` resamples whole species: a resample draws 23 species with replacement and takes every dataset of each species drawn.

Datasets of one species share genomes and lineages, so the species interval is the one to quote. Each interval leaves the same number of resamples outside each end. The output records the seed (20261111) and the number of repetitions (4000 by default), so a run gives the same output for the same input.

**`summary.json`** comes from `python3 -I scripts/release-summary/summary.py INDEX.json INTERVALS_DIR`. It holds the same medians without intervals, taken over the datasets where a learner has both splits. It counts the datasets where the best learner drops by more than 5 and by more than 10 points from the random to the lineage split, and how often each learner is the best. It gives the best learner against the known-gene rule on each split, and the best learner's error-limit counts from the index's point estimates.

It also gives the medians by species, by drug class (classes with at least 5 datasets) and by cohort size (thirds of the datasets by genome count). The script does not read the intervals directory. It takes the argument so that every script has the same command line.

**`catalog.json`** comes from `python3 -I scripts/release-summary/catalog.py INDEX.json INTERVALS_DIR`. It sorts the datasets into bands by the known-gene rule's own balanced accuracy: under 0.55, 0.55 to 0.70, 0.70 to 0.85, 0.85 and up. A band holds its lower edge and not its upper one, so a rule at exactly 0.70 is in "0.70 to 0.85".

For each band it gives the datasets, how often the best learner is ahead on each split, how often the rule is ahead or within 2 points on the lineage split, and the median lineage difference (best learner minus rule). It also counts the datasets where the rule makes a constant call.

Where the rule sits under 0.55, the output says how the rule fails. It finds under 20% of the resistant genomes (sensitivity under 0.2), or it calls over 80% of the susceptible genomes resistant (specificity under 0.2), or neither.

The bands split the datasets on an outcome. A rule that scored high by chance on a small dataset is sorted into a high band, so the bands describe the data and are not a test.

**`index-counts.json`** comes from `python3 -I scripts/release-summary/index_counts.py INDEX.json INTERVALS_DIR`. For each row of each intervals file it rebuilds the two denominators and the two error counts from index fields alone, and compares them with the file's `x` and `n`. The denominators follow the release's scoring rule.

In revision 12, and on an unflagged row of revision 14, a learner row leaves out its `no_call` genomes. Those are the genomes of folds where the fitted model learned nothing. VME's `n` is the resistant genomes less the resistant no-calls, and ME's `n` is the susceptible genomes less the susceptible no-calls.

In revision 14, a learner row whose split entry carries `null_fits_scored: true` is scored on every genome ([equal-scoring.md](equal-scoring.md)). Each `n` is every genome of the class. In both cases `x` is the rate times `n`, rounded. The known-gene rule has no no-calls.

**`compare.py`** takes two outputs: `python3 -I scripts/release-summary/compare.py A/uncertainty.json B/uncertainty.json` prints the numbers of two `uncertainty.json` files side by side, one line per statistic. It recomputes nothing.

## The error limits

The error limits are VME at most 1.5% and ME at most 3%. They are an adapted research screen. They are not an applicable FDA criterion, and meeting them does not make a learner fit for clinical use.

For each learner, split and dataset, the script sorts the dataset into one of three classes by the exact 95% intervals:

- **shown to meet both:** the upper bound of VME is at most 1.5% and the upper bound of ME is at most 3%. An upper bound exactly at the limit meets it.
- **shown to miss one:** the lower bound of VME is over 1.5%, or the lower bound of ME is over 3%. A lower bound exactly at the limit does not miss it.
- **not decided:** every other dataset.

The three classes add up to the datasets with both rates. The output also counts the datasets whose point estimates meet both limits. A dataset with no very major error can be shown to meet the VME limit only if it has enough resistant genomes for the exact upper bound to fall to 1.5% or below: a dataset needs at least 245 resistant genomes. 132 datasets of each release have that many.

## The numbers for revisions 12 and 14

Every number in these tables comes from the files under `reference/release-summary/`. `tests/release_summary/test_document.py` rebuilds each table from those files and fails if the document differs. Values are rounded to 4 decimals as the files hold them. A cell "a / b" gives revision 12, then revision 14. Each interval is the species bootstrap interval.

Median balanced accuracy, from `uncertainty.json`, key `medians`:

| Learner | Split | Datasets | Revision 12 | Revision 14 |
| --- | --- | --- | --- | --- |
| Set covering machine | random | 274 / 274 | 0.9165 (0.8914 to 0.9330) | 0.9165 (0.8914 to 0.9330) |
| Set covering machine | lineage | 268 / 268 | 0.8708 (0.8418 to 0.9036) | 0.8708 (0.8418 to 0.9036) |
| Decision tree | random | 274 / 274 | 0.9092 (0.8846 to 0.9281) | 0.9092 (0.8846 to 0.9281) |
| Decision tree | lineage | 268 / 268 | 0.8684 (0.8199 to 0.8965) | 0.8684 (0.8199 to 0.8965) |
| Random forest | random | 274 / 274 | 0.8970 (0.8676 to 0.9178) | 0.8970 (0.8676 to 0.9178) |
| Random forest | lineage | 268 / 268 | 0.8124 (0.7558 to 0.8510) | 0.8124 (0.7558 to 0.8510) |
| Boosted trees | random | 274 / 274 | 0.9139 (0.8884 to 0.9377) | 0.9139 (0.8884 to 0.9377) |
| Boosted trees | lineage | 268 / 268 | 0.8459 (0.8131 to 0.8924) | 0.8409 (0.7849 to 0.8905) |
| Lasso | random | 274 / 274 | 0.9161 (0.8944 to 0.9335) | 0.9161 (0.8941 to 0.9335) |
| Lasso | lineage | 268 / 268 | 0.8567 (0.8126 to 0.9030) | 0.8565 (0.8107 to 0.9028) |
| Best learner | random | 274 / 274 | 0.9293 (0.9056 to 0.9487) | 0.9293 (0.9068 to 0.9495) |
| Best learner | lineage | 268 / 268 | 0.8790 (0.8451 to 0.9041) | 0.8768 (0.8438 to 0.9036) |

Median paired drop from the random split to the lineage split, from `uncertainty.json`, key `gaps`:

| Learner | Datasets | Revision 12 | Revision 14 |
| --- | --- | --- | --- |
| Set covering machine | 268 / 268 | 0.0242 (0.0116 to 0.0436) | 0.0242 (0.0116 to 0.0436) |
| Decision tree | 268 / 268 | 0.0293 (0.0134 to 0.0489) | 0.0293 (0.0134 to 0.0489) |
| Random forest | 268 / 268 | 0.0504 (0.0297 to 0.0782) | 0.0504 (0.0297 to 0.0782) |
| Boosted trees | 268 / 268 | 0.0463 (0.0247 to 0.0735) | 0.0519 (0.0272 to 0.0875) |
| Lasso | 268 / 268 | 0.0376 (0.0163 to 0.0624) | 0.0394 (0.0161 to 0.0658) |
| Best learner | 268 / 268 | 0.0371 (0.0177 to 0.0552) | 0.0377 (0.0193 to 0.0570) |

Median paired difference on the lineage split, from `uncertainty.json`, key `paired_lineage`:

| Set covering machine minus | Datasets | Revision 12 | Revision 14 | Set covering machine at least as high |
| --- | --- | --- | --- | --- |
| Decision tree | 268 / 268 | 0.0026 (0.0000 to 0.0085) | 0.0026 (0.0000 to 0.0085) | 180 / 180 |
| Random forest | 268 / 268 | 0.0247 (0.0144 to 0.0379) | 0.0247 (0.0144 to 0.0379) | 191 / 191 |
| Boosted trees | 268 / 268 | 0.0046 (0.0000 to 0.0124) | 0.0100 (0.0014 to 0.0191) | 162 / 174 |
| Lasso | 268 / 268 | 0.0020 (0.0000 to 0.0082) | 0.0046 (0.0000 to 0.0118) | 160 / 165 |

The cost of choosing the best learner on the random split, from `uncertainty.json`, key `selection`:

| Statistic | Datasets | Revision 12 | Revision 14 |
| --- | --- | --- | --- |
| Lineage score of the best learner | 268 / 268 | 0.8790 (0.8451 to 0.9041) | 0.8768 (0.8438 to 0.9036) |
| Highest lineage score of the five | 268 / 268 | 0.8966 (0.8591 to 0.9204) | 0.8966 (0.8591 to 0.9204) |
| Highest of the five minus the best learner | 268 / 268 | 0.0024 (0.0002 to 0.0047) | 0.0026 (0.0002 to 0.0049) |

Error-limit classes, from `uncertainty.json`, key `error_bars.rows`:

| Learner | Split | Datasets | Point within both | Shown to meet both | Shown to miss one | Not decided |
| --- | --- | --- | --- | --- | --- | --- |
| Best learner | random | 274 / 274 | 13 / 13 | 2 / 2 | 216 / 214 | 56 / 58 |
| Best learner | lineage | 268 / 268 | 6 / 6 | 0 / 0 | 239 / 239 | 29 / 29 |
| Set covering machine | random | 274 / 274 | 10 / 10 | 2 / 2 | 232 / 232 | 40 / 40 |
| Set covering machine | lineage | 268 / 268 | 5 / 5 | 0 / 0 | 240 / 240 | 28 / 28 |
| Decision tree | random | 274 / 274 | 10 / 10 | 2 / 2 | 235 / 235 | 37 / 37 |
| Decision tree | lineage | 268 / 268 | 8 / 8 | 0 / 0 | 238 / 238 | 30 / 30 |
| Random forest | random | 274 / 274 | 7 / 7 | 1 / 1 | 250 / 250 | 23 / 23 |
| Random forest | lineage | 268 / 268 | 3 / 3 | 0 / 0 | 257 / 257 | 11 / 11 |
| Boosted trees | random | 274 / 274 | 7 / 7 | 1 / 1 | 225 / 225 | 48 / 48 |
| Boosted trees | lineage | 268 / 268 | 6 / 6 | 1 / 1 | 244 / 244 | 23 / 23 |
| Lasso | random | 274 / 274 | 8 / 8 | 1 / 1 | 231 / 231 | 42 / 42 |
| Lasso | lineage | 268 / 268 | 5 / 5 | 0 / 0 | 241 / 242 | 27 / 26 |
| Known-gene rule | random | 253 / 253 | 8 / 8 | 0 / 0 | 235 / 235 | 18 / 18 |
| Known-gene rule | lineage | 249 / 249 | 8 / 8 | 0 / 0 | 231 / 231 | 18 / 18 |

The best learner against the known-gene rule, from `catalog.json`, key `bands`:

| Rule's balanced accuracy | Datasets | Best learner ahead, random split | Best learner ahead, lineage split | Rule ahead or within 2 points, lineage split | Median lineage difference |
| --- | --- | --- | --- | --- | --- |
| under 0.55 | 104 / 104 | 103 / 104 | 98 / 95 of 101 / 101 | 4 / 6 | 0.3561 / 0.3465 |
| 0.55 to 0.70 | 36 / 36 | 36 / 36 | 28 / 28 of 36 / 36 | 8 / 8 | 0.2339 / 0.2339 |
| 0.70 to 0.85 | 32 / 32 | 23 / 23 | 15 / 13 of 32 / 32 | 18 / 20 | -0.0043 / -0.0471 |
| 0.85 and up | 81 / 81 | 35 / 35 | 11 / 11 of 80 / 80 | 74 / 74 | -0.0244 / -0.0244 |
| 0.55 and up | 149 / 149 | 94 / 94 | 54 / 52 of 148 / 148 | 100 / 102 | -0.0103 / -0.0113 |
| all | 253 / 253 | 197 / 198 | 152 / 147 of 249 / 249 | 104 / 108 | 0.085 / 0.0694 |

The rule has an entry on 253 datasets. Its balanced accuracy is exactly 0.5 on 59 of them: it calls nothing resistant on 29 and everything resistant on 30. These counts are the same in both releases (`catalog.json`, key `rule`).

In both releases the counts of all 6424 error rates of the 3212 intervals rows rebuild from the index with 0 mismatches. In revision 14, 122 learner rows carry the flag, and each of them has a no-call count. (`index-counts.json`.)

Between the two releases, only rows that involve boosted trees, the lasso or the best learner change. Those are the two learners that can have a null fit, and the best learner can be one of them ([equal-scoring.md](equal-scoring.md)).

## Limits

- **The species intervals are rough.** The species bootstrap has 23 clusters, and some species hold many datasets. A percentile interval from so few clusters is approximate. Read it to one decimal of a point at most.
- **The exact intervals treat related genomes as independent.** The Clopper-Pearson interval assumes each genome is an independent draw. Genomes of one lineage or one study are not. The intervals are therefore narrower than the data support, and too many datasets are counted as shown.
- **One fold draw per dataset.** Each report scores one assignment of genomes to folds. Another draw would move each score, and the scripts cannot see by how much.
- **The best of five is chosen on the folds that score it.** The best learner is the one with the highest random-split score, and that same score is reported for it. Its random-split median is therefore biased upward. Its lineage score is not used to choose it.
- **The bands split on an outcome.** See catalog.json above.

## Changes against the earlier runs

The scripts ran from outside the repository before they moved here. Two defects were fixed when they moved. Each changed only the numbers named here.

**The bootstrap interval's lower end.** The earlier script took the lower end at the order statistic floor(0.025 (N - 1)) and the upper end at round(0.975 (N - 1)). With 4000 resamples that left 99 resamples below the interval and 100 above it. Both ends now leave round(0.025 (N - 1)) resamples outside, 100 at each end.

Only lower ends moved, each by 0.0001 to 0.001: 15 intervals of revision 12 and 14 of revision 14. No point estimate, count or upper end changed.

**Constant calls of the known-gene rule.** `rule.calls_nothing_resistant` counted every entry with sensitivity 0, and `rule.calls_everything_resistant` every entry with specificity 0. Four entries are not constant calls. Three miss every resistant genome but call some susceptible genomes resistant, and one calls every susceptible genome resistant but misses one resistant genome.

The counts now need sensitivity 0 and specificity 1, or sensitivity 1 and specificity 0. In both releases 32 became 29 and 31 became 30.

`summary.json` is byte-identical to the earlier runs. The earlier per-species sort put a median gap of exactly 0.0 after any negative gap, and it now sorts by value. No species has a negative median gap, so no byte changed.

## Commands

Run these from the repository root, with `R` an unpacked release directory (the one that holds `index.json`) and `OUT` an output directory. They need Python 3.12 and its standard library only.

```sh
python3 -I scripts/release-summary/uncertainty.py  "$R/index.json" "$R/intervals" > "$OUT/uncertainty.json"
python3 -I scripts/release-summary/summary.py      "$R/index.json" "$R/intervals" > "$OUT/summary.json"
python3 -I scripts/release-summary/catalog.py      "$R/index.json" "$R/intervals" > "$OUT/catalog.json"
python3 -I scripts/release-summary/index_counts.py "$R/index.json" "$R/intervals" > "$OUT/index-counts.json"
python3 -I scripts/release-summary/compare.py A/uncertainty.json B/uncertainty.json
```

`reference/release-summary/README.md` gives the release assets, the input sha256 values and the output sha256 values. The tests run with `python3 -I -B -W error -m unittest discover -s tests/release_summary -t .`.
