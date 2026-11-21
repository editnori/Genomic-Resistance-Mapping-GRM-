# AMRFinderPlus comparison

This page compares three callers on the same genomes and recorded labels of public BV-BRC datasets: the AMRFinderPlus program read under each report's rule, GRM's known-gene rule and GRM's five learners. This is a retrospective comparison inside one public database, BV-BRC. It is not a prospective test and not external validation. Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both read is newer than most of these genomes and may have been curated from the same public isolates, which favours both against the learners.

The program read the genomes of the 253 datasets of release `rev14-274-20261112` that have a published rule, 23,751 genomes in all. The headline set, primary-new-235, is the 235 datasets with an organism option and a lineage split whose calls were not read before the plan was hashed. Every figure on this page comes from the final scoring of 16 Nov, and the section Review gives what an independent review narrowed.

## Reading

GRM scored three callers on the same genomes and recorded labels of 235 public BV-BRC datasets with a lineage split: the AMRFinderPlus program (4.2.7, database 2026-08-07.1) read under each report's published rule, GRM's own known-gene rule, and GRM's five learners scored by cross-validation, all inside that one public database, with no test on data from outside it.

In balanced accuracy, the AMRFinderPlus program read under the report's rule was ahead of GRM's known-gene rule in 94 of 235 datasets, level in 102 and behind in 39, with a median difference of +0.0000 and a mean of +0.0642, which sits above the median because the program was ahead by 0.10 or more in 48 of the 235, 20 of them quinolone datasets (exploratory).

On the lineage split of the same 235 datasets, each GRM learner minus the program had a median balanced accuracy of +0.0048 for scm, +0.0000 for cart, -0.0068 for rf, -0.0026 for gbt and -0.0031 for lasso (each 95% species interval spans zero), with median sensitivity -0.0508, -0.0790, -0.0806, -0.0715 and -0.0631 and median specificity +0.0251, +0.0273, +0.0234, +0.0322 and +0.0266 in the same order (sensitivity and specificity exploratory).

Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both read is newer than most of these genomes and may have been curated from the same public isolates, which favours both against the learners, and a catalog row that meets a rule is a marker and an association.

## Plan and order

The plan (sha256 `5cae853c`) was hashed on 14 Nov, while the program ran and before any call of the datasets outside the first 2 was scored. Those 2 were read before the plan, so primary-new-235 leaves them out. The plan fixed the sets, the four readings of the program, the measures, the species interval and rule M.

Its headline is the program against the known-gene rule on primary-new-235 in each measure, and each learner and the served learner against the program in balanced accuracy on the lineage split. Every other summary it lists is secondary, and a figure it did not list is exploratory.

Three addenda (sha256 `d6eda978`, `c8b8500e` and `a6239776`) change where and in what order the program ran and which repeat checks stand in for lost files. Each says it changes no comparison, set, measure, interval, reading or rule of the plan. The second was hashed after part of the cohort had been scored, and it changes only the order of the genomes still to run.

| Event | Time (UTC) | File | sha256 |
|---|---|---|---|
| The plan hashed | 2026-11-14 17:36:14Z | `PLAN.md` | `5cae853c` |
| The first addendum hashed | 2026-11-15 11:48:38Z | `PLAN-ADDENDUM-20261115.md` | `d6eda978` |
| The second addendum hashed | 2026-11-15 22:16:11Z | `PLAN-ADDENDUM-20261115-2.md` | `c8b8500e` |
| The third addendum hashed | 2026-11-16 03:08:48Z | `PLAN-ADDENDUM-20261116.md` | `a6239776` |
| The final scoring started | 2026-11-16 06:18:24Z |  |  |

## Method

- **Datasets.** A dataset is one species and one antibiotic of the release with a published rule. Its genomes and labels are those the release reports for it.
- **The AMRFinderPlus program.** Version `4.2.7` with database `2026-08-07.1`, given the organism option where the species has one. A genome is called resistant when every group of the report's rule is met by at least one counted row. `tool_default`, the primary reading, counts rows of core scope; `tool_default_intact` counts no row whose method the program marks as partial or as an internal stop; `genes_core` counts gene rows of core scope only; and `genes_all_scopes` counts gene rows of core or plus scope. The record keeps each reading's counts per dataset in `callers.tsv` and the paired differences of `tool_default` alone in `differences.tsv`, and this page reads `tool_default` alone.
- **GRM's known-gene rule.** The call the release itself makes for each genome from the catalog. It is not the program's call, and the page never gives it the program's name.
- **GRM's learners.** scm, cart, rf, gbt and lasso, from the scores the release's reports hold, calling resistant at a score of 0.5 or more, on the random split, where training and test share lineages, and on the lineage split. The learner the release serves is picked on the random split.
- **Measures.** Sensitivity, specificity and balanced accuracy, their mean, from each caller's counts on each dataset. A difference is the first caller minus the second on one dataset, exact, so level means a difference of exactly 0.
- **Over datasets.** The counts ahead, level and behind, the median and quartiles, the mean and the exact two-sided sign test. The species interval draws the set's species with replacement, as many as it has, takes every dataset of each drawn species, and repeats 4,000 times; its ends are the 101st and the 3,900th sorted value. A set with fewer than 6 species gets none.
- **Generator and draw.** The plan fixed the seed string and named neither the generator nor the draw. The record's script uses Python's `random.Random` seeded with that string, `randrange` over the species sorted by name and `statistics.median`, and gives the final scoring's ends at 6 decimals.
- **Rounding.** The page prints each measure at 4 decimals, rounded from the record's cells at 6 decimals.

## The program against the known-gene rule

Each row gives the program read under the report's rule minus GRM's known-gene rule over a set's datasets. primary-new-235 is the plan's headline, and primary-237 adds the 2 datasets read before the plan.

| Set | Measure | Role | Datasets | Ahead | Level | Behind | Median (95% species interval) | Mean (95% species interval) |
|---|---|---|---|---|---|---|---|---|
| primary-new-235 | balanced accuracy | primary | 235 | 94 | 102 | 39 | +0.0000 (+0.0000 to +0.0000) | +0.0642 (+0.0254 to +0.1070) |
| primary-new-235 | sensitivity | primary | 235 | 87 | 117 | 31 | +0.0000 (+0.0000 to +0.0012) | +0.1406 (+0.0889 to +0.2286) |
| primary-new-235 | specificity | primary | 235 | 42 | 117 | 76 | +0.0000 (+0.0000 to +0.0000) | -0.0122 (-0.1267 to +0.0829) |
| primary-237 | balanced accuracy | secondary | 237 | 96 | 102 | 39 | +0.0000 (+0.0000 to +0.0000) | +0.0665 (+0.0270 to +0.1095) |
| primary-237 | sensitivity | secondary | 237 | 89 | 117 | 31 | +0.0000 (+0.0000 to +0.0014) | +0.1462 (+0.0925 to +0.2375) |
| primary-237 | specificity | secondary | 237 | 42 | 117 | 78 | +0.0000 (+0.0000 to +0.0000) | -0.0131 (-0.1245 to +0.0853) |

Of the 102 datasets level in balanced accuracy, 100 have the same call on every genome (exploratory). On the 26 quinolone datasets among the 235, the program minus the known-gene rule has a median of +0.3301 in balanced accuracy (exploratory).

## Learners against the program

Each row gives one learner minus the program on the lineage split of primary-new-235, in balanced accuracy with the counts ahead, level and behind, and in sensitivity beside specificity. The served learner is the learner the release serves, picked on the random split, and the highest of the five is an oracle no one can pick in advance.

| Learner | Role (balanced accuracy) | Ahead, level, behind | Median balanced accuracy (95% species interval) | Median sensitivity (95% species interval) | Median specificity (95% species interval) |
|---|---|---|---|---|---|
| scm | primary | 126, 5, 104 | +0.0048 (-0.0063 to +0.0466) | -0.0508 (-0.1056 to -0.0240) | +0.0251 (+0.0000 to +0.1208) |
| cart | primary | 112, 6, 117 | +0.0000 (-0.0174 to +0.0254) | -0.0790 (-0.1188 to -0.0335) | +0.0273 (+0.0000 to +0.1031) |
| rf | primary | 101, 5, 129 | -0.0068 (-0.0361 to +0.0039) | -0.0806 (-0.1501 to -0.0352) | +0.0234 (+0.0000 to +0.0886) |
| gbt | primary | 108, 5, 122 | -0.0026 (-0.0284 to +0.0218) | -0.0715 (-0.1150 to -0.0310) | +0.0322 (-0.0039 to +0.1102) |
| lasso | primary | 111, 3, 121 | -0.0031 (-0.0163 to +0.0140) | -0.0631 (-0.1239 to -0.0291) | +0.0266 (-0.0036 to +0.1216) |
| served learner | primary | 122, 5, 108 | +0.0031 (-0.0081 to +0.0411) | -0.0541 (-0.1050 to -0.0233) | +0.0322 (+0.0000 to +0.1155) |
| highest of the five | exploratory | 142, 8, 85 | +0.0208 (+0.0004 to +0.1081) | -0.0411 (-0.0702 to -0.0190) | +0.0322 (+0.0000 to +0.1205) |

Exploratory: each learner's median sensitivity is lower than the program's (-0.0806 to -0.0508), and every interval ends below zero (-0.0352 to -0.0240); each learner's median specificity is higher than the program's (+0.0234 to +0.0322), and each interval starts at zero or just below it (gbt -0.0039, lasso -0.0036).

On the lineage split, the learner the release serves, picked on the random split, was ahead of the program in 122 of 235 datasets, level in 5 and behind in 108, with a median difference of +0.0031 (species 95% -0.0081 to +0.0411); the two were within 0.02 of each other in 56 of 235 datasets. It was at most 0.02 ahead of the program, or behind it, in 132 of 235, which is within 0.02 as rule M defines it (the learner at most 0.02 ahead). Both counts at 0.02 are exploratory.

The served learner has the highest lineage balanced accuracy of the five in 100 of 235 datasets (exploratory).

The next table is the random split, secondary for the five learners and the served learner: training and test share lineages, so it is no test on unseen lineages, and the served learner is picked on it.

| Learner, random split: training and test share lineages, and the served learner is picked on this split | Role | Ahead | Level | Behind | Median balanced accuracy (95% species interval) | Mean (95% species interval) |
|---|---|---|---|---|---|---|
| scm | secondary | 162 | 6 | 67 | +0.0482 (+0.0057 to +0.1282) | +0.1231 (+0.0789 to +0.1695) |
| cart | secondary | 156 | 6 | 73 | +0.0341 (+0.0041 to +0.1272) | +0.1160 (+0.0667 to +0.1657) |
| rf | secondary | 152 | 3 | 80 | +0.0404 (+0.0004 to +0.1210) | +0.1037 (+0.0559 to +0.1529) |
| gbt | secondary | 157 | 2 | 76 | +0.0541 (+0.0048 to +0.1454) | +0.1214 (+0.0745 to +0.1695) |
| lasso | secondary | 158 | 5 | 72 | +0.0571 (+0.0070 to +0.1346) | +0.1226 (+0.0767 to +0.1730) |
| served learner | secondary | 177 | 8 | 50 | +0.0628 (+0.0127 to +0.1681) | +0.1406 (+0.0919 to +0.1886) |
| highest of the five | exploratory | 177 | 8 | 50 | +0.0628 (+0.0131 to +0.1681) | +0.1406 (+0.0925 to +0.1881) |

## Quinolone and carbapenem datasets

These sets are drug classes of the datasets. The table gives each set's counts, median and mean in balanced accuracy, with the species interval where a set has at least 6 species.

| Set | Contrast | Role | Datasets | Species | Ahead | Level | Behind | Median | Mean |
|---|---|---|---|---|---|---|---|---|---|
| rule-M-quinolone-28 | program minus known-gene rule | secondary | 28 | 14 | 24 | 2 | 2 | +0.3301 (+0.1612 to +0.4164) | +0.2763 (+0.1690 to +0.3870) |
| rule-M-quinolone-28 | served learner minus program, lineage split | secondary | 28 | 14 | 15 | 3 | 10 | +0.0067 (-0.0006 to +0.0723) | +0.0667 (+0.0174 to +0.1211) |
| rule-M-quinolone-28 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | secondary | 28 | 14 | 22 | 4 | 2 | +0.0481 (+0.0064 to +0.2297) | +0.1318 (+0.0439 to +0.2289) |
| rule-M-quinolone-new-26 | program minus known-gene rule | exploratory | 26 | 13 | 22 | 2 | 2 | +0.3301 (+0.0130 to +0.4157) | +0.2715 (+0.1580 to +0.3884) |
| rule-M-quinolone-new-26 | served learner minus program, lineage split | secondary | 26 | 13 | 13 | 3 | 10 | +0.0022 (-0.0019 to +0.0853) | +0.0691 (+0.0187 to +0.1291) |
| rule-M-quinolone-new-26 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | exploratory | 26 | 13 | 20 | 4 | 2 | +0.0477 (+0.0044 to +0.3733) | +0.1376 (+0.0413 to +0.2398) |
| quinolone-organism-option-30 | program minus known-gene rule | secondary | 30 | 14 | 26 | 2 | 2 | +0.3547 (+0.2084 to +0.4316) | +0.2900 (+0.1778 to +0.3974) |
| quinolone-organism-option-30 | served learner minus program, lineage split | secondary | 28 | 14 | 15 | 3 | 10 | +0.0067 (-0.0006 to +0.0723) | +0.0667 (+0.0186 to +0.1236) |
| quinolone-organism-option-30 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | secondary | 30 | 14 | 23 | 4 | 3 | +0.0471 (+0.0057 to +0.2124) | +0.1231 (+0.0399 to +0.2213) |
| carbapenem-20 | program minus known-gene rule | secondary | 20 | 9 | 7 | 7 | 6 | +0.0000 (-0.0714 to +0.0146) | -0.0248 (-0.0674 to +0.0139) |
| carbapenem-20 | served learner minus program, lineage split | secondary | 20 | 9 | 11 | 0 | 9 | +0.0363 (-0.1265 to +0.2200) | +0.0453 (-0.0800 to +0.1870) |
| carbapenem-20 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | secondary | 20 | 9 | 14 | 0 | 6 | +0.0693 (-0.0140 to +0.4289) | +0.1434 (+0.0304 to +0.2882) |

Rule M, a proposed screen from a review note that no one adopted, counts the quinolone datasets where the served learner is at most 0.02 ahead of the program or behind it. That was 16 of 28 where 14 were needed. Two of the 16 count by the served learner scoring below a program at 0.5000 (balanced accuracy 0.4453 and 0.3502; exploratory), and the two were within 0.02 of each other in 11 of 28. It is an adapted research screen, never an applicable FDA validation criterion.

## Datasets kept apart

Twelve datasets in two species for which the program has no organism option are reported apart: there the program reads gene rows only.

The 4 datasets with an organism option and no lineage split are read on the random split only, where training and test share lineages.

Both sets have 2 species, too few for a species interval.

| Set | Contrast | Role | Datasets | Species | Ahead | Level | Behind | Median | Mean |
|---|---|---|---|---|---|---|---|---|---|
| no-organism-option-12 | program minus known-gene rule | secondary | 12 | 2 | 0 | 10 | 2 | +0.0000 | -0.0001 |
| no-organism-option-12 | served learner minus program, lineage split | secondary | 12 | 2 | 10 | 0 | 2 | +0.2961 | +0.2296 |
| no-organism-option-12 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | secondary | 12 | 2 | 10 | 0 | 2 | +0.3247 | +0.2563 |
| organism-option-no-lineage-split-4 | program minus known-gene rule | secondary | 4 | 2 | 2 | 2 | 0 | +0.2355 | +0.2410 |
| organism-option-no-lineage-split-4 | served learner minus program, random split: training and test share lineages, and the served learner is picked on this split | secondary | 4 | 2 | 3 | 0 | 1 | +0.0301 | +0.1293 |

## Review

An independent review recomputed the figures from the final scoring's tables and the release before this record was written. It narrowed the reading in these ways:

- The counts ahead, level and behind and the median lead. The mean follows with the reason it sits above the median, and no sentence joins ahead with level.
- At 0.02 the page gives the two-sided count. The one-sided count, the learner at most 0.02 ahead of the program or behind it, appears only with that definition.
- The sentence on the catalog travels with the figures that set a learner beside the program: in the first paragraph, the section Reading and the record's README.
- Each row and figure is labelled primary, secondary or exploratory as the plan fixed it, and the 26 quinolone datasets of primary-new-235 are named by drug class alone.
- Every random-split figure carries the label that training and test share lineages and that the served learner is picked on that split.
- The two disagreement tables of the final scoring stay out of the record: their counts by the catalog's subtype and class labels, set beside a species and a drug, can read as naming how resistance forms, and some of their rows give the program's call on one genome.
- Rule M appears only in the review's words, as a proposed screen no one adopted and never an applicable FDA validation criterion.

## What the record cannot show

- The plan's sensitivity run on the newest catalog database (2026-10-29.1) was not run, so these figures hold for database 2026-08-07.1 only.
- The program's tables per genome are gone, so its figures can be checked only at the counts per dataset and per fold, and its paired lineage intervals not at all. The record holds the counts per dataset.
- Joined with the release's reports, which give each genome's label and known-gene rule call, the counts per dataset give the program's call on some genomes, for example on every genome of a dataset where the two callers agree on every genome.
- The species interval depends on a generator and a draw the plan did not name. The record holds its ends, and not the seed string the plan fixed.
- Every caller is scored against the same recorded labels of one database, BV-BRC, with one curation.
- The page gives no figure for a single dataset, and it names no gene and no catalog row.

## Files

- `reference/amrfinderplus-comparison/` holds the record: 11 files that `scripts/amrfinderplus-comparison/record.py` writes from the final scoring. Its README lists each file and its columns, and `SHA256SUMS` gives each file's sha256.
- The record holds rows per dataset and caller, per set and contrast, and per learner, and single figures. It holds no genome id, no label or call per genome and no rule.
- The final scoring's tables, the plan and the set lists are not public, so the repository holds the record alone. `--check` alone derives every derived file from the row files.
- `tests/amrfinderplus_comparison` rebuilds the section Reading and every table from the record, holds every other paragraph and list item word for word, and checks that each other number on the page stands in the record.

```sh
python3 -I scripts/amrfinderplus-comparison/record.py --check reference/amrfinderplus-comparison
python3 -I -m unittest discover -s tests/amrfinderplus_comparison -t .
```
