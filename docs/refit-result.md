# Refit result

Where a learner misses on new lineages and the AMRFinderPlus program does not, why? That is the assumption plan's diagnostic D0. Its block C8 asks whether a learner chooses the same inputs from fold to fold. This page reads both, once, over a refit of the release's rule on 238 datasets in 23 species (`docs/assumption-plan.md`, rows D0 and C8).

The rules of the reading were fixed before any D0 or C8 count of the run was read, a rehearsal on *Klebsiella aerogenes* aside, and were public before the reading program ran (`docs/assumption-plan-corrections.md`, entries "D0 causes" and "Refit reading"). The program ran once, on 22 Nov 2026. The code is in `docs/assumption-refit.md`.

D0's datasets are chosen by their result. So its counts order the plan's comparisons and test nothing. C8 on the release decides nothing alone. Nothing changes in the release. Every split here is a retrospective split inside one public database, and a split by lineage is never external validation.

Every figure is a count, a share or an overlap over datasets. No genome, sequence or feature is written, and no new marker is found, named or ranked.

## Reading

**Order.** Absent is the largest cause. It holds 107 of the 177 miss pairs (0.605). So C3, the comparison of other inputs than 31-mer presence, runs first. C4 and C5 follow in the default's order, and C7 comes after the three.

**Absent.** A miss pair is a dataset and a learner where the program beats the learner by 0.10 or more. A missed fold is Absent when no row of the 31-mer matrix agrees with the program's call at 0.95 or more on the fold's training genomes. The learner's input does not hold what the program reads as one row.

**Others.** Unused holds 33 pairs: the learner had a matching row and did not choose it. Chosen holds 28 pairs, Removed 15 pairs and Lineage 1 pair.

**Stable.** The plan's MISS set holds 204 miss pairs, and 27 of them are not read yet: they lie in 10 datasets that wait for the maintainer's label rule. With all 27 added to the next cause, Absent still leads. So the first block is stable.

**Overlap.** On the lineage split, the five folds of one dataset share few chosen features. The median overlap is 0.050 for `scm`, 0.053 for `cart`, 0.183 for `rf`, 0.075 for `gbt` and 0.066 for `lasso`. On the random split each is higher: 0.155, 0.125, 0.250, 0.123 and 0.126.

**Families.** Of the fold-chosen features that have a card, between 0.101 (`lasso`) and 0.259 (`scm`) have a card that names a family of the AMRFinderPlus catalogue, on the lineage folds.

## Data

**Refit.** The refit is `grm learn` with the release's rule at seed 42, with each fold's model saved. It ran on the 238 datasets of the release that are not marked M, one species a machine. The run's own check compared each report with the release's: 238 of 238 dataset rows and 2,790 of 2,790 learner rows read equal.

**Marks.** M marks the 36 datasets of the release where label rules (a), (b) and (d) do not give the same labels. They wait for the maintainer's rule. The 238 have the same labels under label rules (a), (b) and (d). Under rule (c), the most resistant call, a dataset holding a pair with an R and an S record is checked again.

**Program.** The program is the AMRFinderPlus program, run on the same genomes. D0 reads its per-genome calls from the tables of the comparison in `docs/amrfinderplus-comparison.md`. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalog and keeps its own name.

**MISS.** The plan's MISS set is each learner's own: the datasets where the program's balanced accuracy on the lineage split beats the learner's by 0.10 or more. It holds 75 datasets in 16 species and 204 miss pairs: 26 for `scm`, 34 for `cart`, 59 for `rf`, 48 for `gbt` and 37 for `lasso`.

**Read here.** D0 reads the 65 MISS datasets that are not M, in 16 species. They give 177 miss pairs: 23 for `scm`, 28 for `cart`, 52 for `rf`, 41 for `gbt` and 33 for `lasso`. The 177 pairs are part of the 204.

**Missed fold.** A lineage fold counts as missed when it tests both classes and the program's balanced accuracy on it beats the learner's by 0.10 or more. Each missed fold takes the first cause that holds, in the order below. A pair goes to the cause that covers most of its missed folds. A pair whose folds tie counts in each tied cause.

1. **Absent.** No feature matches the program's call on the fold's training genomes at 0.95 agreement or more.
2. **Removed.** A feature matches, but its support is under the floor, or its chi-square rank is past 4,000. The run counts the two together.
3. **Lineage.** A matching feature is offered, but the learner chose one that agrees with it at 0.95 or more and lies in at most two lineage clusters.
4. **Chosen.** The learner chose a feature that agrees with an offered match and lies in more than two lineage clusters. It still missed.
5. **Unused.** The learner had the input and did not choose it.

**Unit.** C8's unit is a dataset and a learner. Its value is the median, over the 10 pairs of the five folds, of the Jaccard overlap of the two folds' chosen features: shared features over features in either. A feature is counted by its matrix row. An overlap of 0.050 is 1 shared feature in 20.

**Chose.** For `scm`, `cart` and `lasso` a fold chose the features of its model. For `rf` and `gbt` it chose its 50 heaviest features, all that a fold model keeps.

## Causes

The 40% rule is taken of miss pairs. The largest cause among Absent, Removed, Lineage and Unused sets the order when it holds 0.40 or more. Chosen counts in the pairs and never sets the order.

| Cause | Miss pairs | Share | Largest | Block |
|---|---|---|---|---|
| Absent | 107 of 177 | 0.605 | yes | C3 |
| Removed | 15 of 177 | 0.085 | no | C4's floor arm |
| Lineage | 1 of 177 | 0.006 | no | C4's adjusted arms |
| Chosen | 28 of 177 | 0.158 | outside the four | none |
| Unused | 33 of 177 | 0.186 | no | C5 |

A tied pair counts in each tied cause, so the five causes hold 184 entries for 177 pairs, and the shares do not add to 1.

**Fixed words.** The reading program wrote these lines:

```text
D0 ordering (D0's sets are chosen by their result: it orders comparisons and tests nothing): absent holds 107 of 177 miss pairs (0.605): C3 first (sets chosen by their result; orders comparisons, tests nothing)
order of the input group: C3, then C4 and C5; then C7, after the three in every outcome (the plan's Order, item 5)
coverage: D0 reads 65 of the plan's 75 MISS datasets and 177 of its 204 miss pairs; the MISS datasets that are M wait for the maintainer's label rule, and the line is read again then
ruling 6: absent holds 107 of the plan's 204 miss pairs (0.525) and the next cause of the four 33; 107 is 0.40 or more of 204 and above 33 + 27: the first block is stable
removed, as counted, covers both the support floor (A8) and the 4,000-row screen (A9), which the copy does not split; the order keeps the plan's words, C4's floor arm (ruling 9)
```

**Each learner.** Absent is the largest cause for each of the five learners. The table counts datasets, each learner on its own miss pairs.

| Learner | Pairs | Absent | Removed | Lineage | Chosen | Unused | Tied | No fold |
|---|---|---|---|---|---|---|---|---|
| `scm` | 23 | 13 | 0 | 0 | 2 | 8 | 0 | 0 |
| `cart` | 28 | 19 | 0 | 0 | 2 | 10 | 3 | 0 |
| `rf` | 52 | 29 | 3 | 1 | 17 | 3 | 1 | 0 |
| `gbt` | 41 | 25 | 7 | 0 | 5 | 7 | 2 | 0 |
| `lasso` | 33 | 21 | 5 | 0 | 2 | 5 | 1 | 1 |
| All five | 177 | 107 | 15 | 1 | 28 | 33 | 7 | 1 |

"Tied" counts pairs whose missed folds tie between causes. "No fold" counts pairs with no missed fold and so no cause: the pair stays in the 177. No pair took the fallback fold of the entry "D0 causes".

**Gap.** A pair's gap is the program's balanced accuracy on the lineage split less the learner's. The table gives each cause's share of a learner's summed gap. A tied pair splits its gap equally. The share pooled over the five learners is not on this page, by the reading's last ruling: it is in `reference/refit-result/d0/d0-summary.tsv`.

| Learner | Absent gap | Removed gap | Lineage gap | Chosen gap | Unused gap |
|---|---|---|---|---|---|
| `scm` | 0.569 | 0.000 | 0.000 | 0.071 | 0.361 |
| `cart` | 0.638 | 0.000 | 0.000 | 0.033 | 0.330 |
| `rf` | 0.585 | 0.041 | 0.038 | 0.292 | 0.044 |
| `gbt` | 0.591 | 0.129 | 0.000 | 0.094 | 0.185 |
| `lasso` | 0.612 | 0.141 | 0.000 | 0.055 | 0.169 |

**Stability.** Ruling 6 of the reading asks whether the 27 unread miss pairs could change the first block. Absent holds 107 of the plan's 204 miss pairs (0.525), and the next cause of the four holds 33 pairs. 107 is 0.40 or more of 204 and is above 33 plus 27. So the first block is stable.

**Again.** D0 reads 65 of the plan's 75 MISS datasets and 177 of its 204 miss pairs. The other 10 datasets are M. The line is read again when the maintainer's label rule is made.

**Removed.** Removed, as counted, covers both the support floor and the 4,000-row screen. The run does not split them. The plan's words for its block stay: C4's floor arm.

**Per dataset.** Each of the 177 pairs is in `reference/refit-result/page-d0-datasets.tsv`, with species and antibiotic in full, both balanced accuracies, the gap, the missed folds and the cause.

## Overlap

| Split | Learner | Defined | First quartile | Median | Third quartile | Under two | Pattern median |
|---|---|---|---|---|---|---|---|
| lineage | `scm` | 233 of 238 | 0.000 | 0.050 | 0.333 | 0 | 0.155 |
| lineage | `cart` | 233 of 238 | 0.000 | 0.053 | 0.310 | 0 | 0.106 |
| lineage | `rf` | 233 of 238 | 0.105 | 0.183 | 0.274 | 0 | none |
| lineage | `gbt` | 233 of 238 | 0.041 | 0.075 | 0.111 | 3 | none |
| lineage | `lasso` | 233 of 238 | 0.028 | 0.066 | 0.104 | 2 | 0.083 |
| lineage | All five | 1,165 of 1,190 | 0.005 | 0.087 | 0.200 | 5 | 0.091 |
| random | `scm` | 238 of 238 | 0.000 | 0.155 | 0.333 | 0 | 0.213 |
| random | `cart` | 238 of 238 | 0.000 | 0.125 | 0.333 | 0 | 0.167 |
| random | `rf` | 238 of 238 | 0.163 | 0.250 | 0.333 | 0 | none |
| random | `gbt` | 238 of 238 | 0.070 | 0.123 | 0.163 | 0 | none |
| random | `lasso` | 238 of 238 | 0.099 | 0.126 | 0.170 | 0 | 0.146 |
| random | All five | 1,190 of 1,190 | 0.067 | 0.145 | 0.258 | 0 | 0.159 |

**Columns.** "Defined" counts the units with a value. The quartiles and the median are taken over those units, and the quartiles use the inclusive method of Python's `statistics.quantiles`. "All five" pools the units of the five learners and decides nothing.

**Under two.** A unit where one fold alone chose a feature keeps its pairs at overlap 0, so its value is 0. A pair where both folds chose nothing is left out. "Under two" counts the units where fewer than two folds chose: 3 for `gbt` and 2 for `lasso` on the lineage split. They are *Acinetobacter baumannii* with cefepime, *Escherichia coli* with amikacin and *Neisseria gonorrhoeae* with ceftriaxone for `gbt`, and *Klebsiella pneumoniae* with cefotetan and with doripenem for `lasso`.

**No lineage split.** The lineage split did not run on 5 datasets, as in the release: *Enterococcus faecium* with chloramphenicol, *Salmonella enterica* with temocillin, and *Staphylococcus aureus* with ceftaroline, with levofloxacin and with moxifloxacin. Their 25 units are undefined on the lineage split.

**Pattern.** Two matrix rows can hold the same pattern. "Pattern median" counts two features as one when their patterns agree on every sample both folds train on. It is given for `scm`, `cart` and `lasso`. The row for all five pools those three.

**No bands.** The plan gives C8 one cut point, "within 0.05", for an arm against the release. No arm has run, so no unit is counted over or under a line here.

**Families.** The table sums, over datasets, the distinct features a learner's folds chose. Only a feature of the final model has a card, so "With card" is given of the distinct features. "With family" counts the cards that name a family of the AMRFinderPlus catalogue, given of the features with a card. No family is named here.

| Folds' split | Learner | Distinct | With card | With family | Family share |
|---|---|---|---|---|---|
| lineage | `scm` | 1,929 | 936 of 1,929 | 242 of 936 | 0.259 |
| lineage | `cart` | 3,463 | 1,148 of 3,463 | 218 of 1,148 | 0.190 |
| lineage | `rf` | 36,824 | 7,903 of 36,824 | 1,159 of 7,903 | 0.147 |
| lineage | `gbt` | 37,599 | 6,178 of 37,599 | 797 of 6,178 | 0.129 |
| lineage | `lasso` | 25,351 | 4,415 of 25,351 | 447 of 4,415 | 0.101 |
| lineage | All five | 105,166 | 20,580 of 105,166 | 2,863 of 20,580 | 0.139 |
| random | `scm` | 2,291 | 1,072 of 2,291 | 255 of 1,072 | 0.238 |
| random | `cart` | 4,549 | 1,383 of 4,549 | 224 of 1,383 | 0.162 |
| random | `rf` | 32,858 | 8,224 of 32,858 | 1,164 of 8,224 | 0.142 |
| random | `gbt` | 38,810 | 8,007 of 38,810 | 886 of 8,007 | 0.111 |
| random | `lasso` | 26,545 | 5,859 of 26,545 | 487 of 5,859 | 0.083 |
| random | All five | 105,053 | 24,545 of 105,053 | 3,016 of 24,545 | 0.123 |

**Per unit.** Each dataset, learner and split is in `reference/refit-result/page-c8-datasets.tsv`: 2,380 rows, with the folds that chose, the pairs, both overlaps and the three feature counts.

## Seen after

These were seen after the reading. They change no rule and no order.

**Alone.** Absent is the one cause of 104 of the 177 pairs. Unused is the one cause of 27 pairs, Chosen of 23 pairs, Removed of 14 pairs and Lineage of 1 pair. Seven pairs tie and 1 pair has no missed fold:

| Species | Antibiotic | Learner | Causes |
|---|---|---|---|
| *Acinetobacter baumannii* | minocycline | `gbt` | Removed, Chosen, Unused |
| *Campylobacter coli* | gentamicin | `cart` | Chosen, Unused |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | `gbt` | Chosen, Unused |
| *Klebsiella pneumoniae* | doripenem | `cart` | Absent, Unused |
| *Klebsiella pneumoniae* | ofloxacin | `lasso` | no missed fold |
| *Salmonella enterica* | kanamycin | `rf` | Absent, Chosen |
| *Salmonella enterica* | sulfamethoxazole | `lasso` | Chosen, Unused |
| *Staphylococcus aureus* | tetracycline | `cart` | Absent, Unused |

**Species.** The Absent pairs are unevenly spread. *Salmonella enterica* holds 40 of the 107 and *Escherichia coli* 20. With *Staphylococcus aureus*, *Enterobacter cloacae* and *Klebsiella pneumoniae*, five species hold 96 of the 107. Seven of the 16 species have no Absent pair.

| Species | D0 datasets | D0 pairs | Absent pairs | Removed pairs | Lineage pairs | Chosen pairs | Unused pairs |
|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 2 | 8 | 0 | 1 | 0 | 2 | 7 |
| *Campylobacter coli* | 5 | 7 | 0 | 0 | 0 | 7 | 1 |
| *Campylobacter jejuni* | 2 | 2 | 0 | 0 | 0 | 2 | 0 |
| *Corynebacterium diphtheriae* | 3 | 7 | 3 | 0 | 0 | 1 | 3 |
| *Enterobacter cloacae* | 5 | 17 | 13 | 0 | 0 | 4 | 1 |
| *Enterococcus faecium* | 1 | 1 | 0 | 0 | 0 | 0 | 1 |
| *Escherichia coli* | 10 | 22 | 20 | 0 | 0 | 1 | 1 |
| *Klebsiella aerogenes* | 3 | 9 | 1 | 4 | 0 | 1 | 3 |
| *Klebsiella oxytoca* | 1 | 5 | 0 | 0 | 0 | 0 | 5 |
| *Klebsiella pneumoniae* | 9 | 24 | 9 | 7 | 0 | 2 | 6 |
| *Pseudomonas aeruginosa* | 1 | 3 | 3 | 0 | 0 | 0 | 0 |
| *Salmonella enterica* | 12 | 48 | 40 | 3 | 0 | 3 | 4 |
| *Shigella flexneri* | 2 | 4 | 4 | 0 | 0 | 0 | 0 |
| *Shigella sonnei* | 2 | 2 | 0 | 0 | 1 | 1 | 0 |
| *Staphylococcus aureus* | 5 | 15 | 14 | 0 | 0 | 1 | 1 |
| *Streptococcus pneumoniae* | 2 | 3 | 0 | 0 | 0 | 3 | 0 |

The table counts the rows of the per-dataset file. A tied pair counts in each of its causes.

**Meaning.** Absent says that no single 31-mer row follows the program's call. It does not say why. One possible reason is that the program's call can rest on any of several genes or point mutations, which one row cannot follow. C3's arms try other inputs: the program's calls themselves, unitigs, k-mer counts and amino-acid k-mers.

**Screen.** All 15 Removed pairs are on `rf`, `gbt` and `lasso`, the three learners with the 4,000-row screen. `scm` and `cart` have the floor alone and no Removed pair. This is a pattern across learners. The run does not count the screen apart from the floor.

**Lineage.** Lineage holds 1 pair of 177: `rf` on *Shigella sonnei* with colistin. By D0's rule, a learner that chose a marker of one or two lineages in place of an offered match is rare in these pairs.

**Chosen.** Chosen holds 17 of its 28 pairs on `rf`. All 7 pairs of *Campylobacter coli* are Chosen: the learner chose a feature that agrees with the program's call and still trailed the program by 0.10 or more.

**Unused.** Unused is the second cause for `scm` (8 of 23 pairs) and `cart` (10 of 28 pairs). Their folds chose the fewest distinct features on the lineage folds: 1,929 features and 3,463 features. It holds 3 of 52 pairs for `rf`.

**Random split.** Each learner's median overlap is higher on the random split than on the lineage split, and the first quartile for `scm` and `cart` is 0.000 on both. In at least a quarter of the datasets, most pairs of folds of these two learners share no chosen feature.

## Limits

**Chosen sets.** D0's datasets are those where the program beats a learner by 0.10 or more. A set chosen by its result tests nothing. The counts say where to look first.

**Catalog.** The comparison is retrospective. The catalog the program reads is newer than most of these genomes and may have been curated from the same public isolates, which favours the program against the learners.

**Coverage.** D0 reads 65 of 75 MISS datasets. The 10 that are M hold the other 27 pairs, and the line is read again when they are read.

**Match.** A match is agreement at 0.95 or more on the training genomes the program calls. The line was fixed before any count was read and was not varied. A different line would move pairs between Absent and the other causes.

**Trees.** For `rf` and `gbt` a feature used below the 50 heaviest counts as not chosen, in D0 and in C8. Their Unused and Chosen counts and their overlaps are read on that cut list.

**One seed.** The refit has one fold seed, 42, and no draws. C8's question names folds and draws, and this page answers for folds.

**Cards.** On the lineage folds, 20,580 of the 105,166 distinct fold-chosen features have a card. The family share is of those 20,580 features and says nothing of the others.

**One database.** Every genome comes from BV-BRC. The lineage split holds lineages out inside that one database.

**Alone.** C8 decides nothing alone. Its use is as the release's row, against which each arm's overlap is later counted as higher, lower or within 0.05.

## Method

The reading program stops and writes no table unless every gate passes. All five passed:

| Gate | Passed |
|---|---|
| Species | 23 of 23 species, each with `DONE` and the one code sha256 of the build |
| Datasets | 238 of 238 datasets, none twice, and the check tables `datasets.tsv` and `learners.tsv` read equal |
| Marks | the 238 datasets are exactly the label count's datasets that are not M |
| D0 sets | 65 datasets in 16 species, each species' list being the MISS set among its datasets |
| D0 pairs | 177 of 177 miss pairs, with each learner's own count |

**Run.** The program ran once, at 21:40Z on 22 Nov 2026, at commit `0f5cc5c`, and ended with status 0. The build's code sha256 is `20e08189285ece9ac29a03f989469899356d74e5876927e689f16a1a037c104c` in all 23 species.

```sh
python3 -I scripts/assumption-d0-c8/read_refit.py --out <out> \
  --code-sha256 20e08189285ece9ac29a03f989469899356d74e5876927e689f16a1a037c104c \
  --species 23 --datasets 238 --d0-datasets 65 --d0-species 16 --d0-pairs 177 \
  --d0-learner-pairs scm=23,cart=28,rf=52,gbt=41,lasso=33 \
  --marks <the label count's datasets.tsv> \
  --equal-tables datasets.tsv,learners.tsv \
  --c8-one-fold zero --c8-over-five pairs --c8-bands none --order-rest default-after \
  <the 23 species' copied directories>
```

**Files.** The reading's tables are in `reference/refit-result/`. Its `summary.txt` and `gates.txt` are the program's, with one local path replaced by the words "the label count's datasets.tsv".

**Tables.** `scripts/assumption-d0-c8/page_tables.py` prints this page's tables and its species paragraphs from those files, with shares and overlaps rounded to 3 decimals. `tests/assumption_d0_c8/test_page.py` checks each printed line against this page.

## Datasets

The table gives the refit's check per species. "In D0" counts the datasets D0 reads. "Dataset rows" and "Learner rows" count the rows of the two check tables, all equal to the release.

| Species | Datasets | In D0 | Dataset rows | Learner rows | Breakpoints | Lists | Samples | Cards |
|---|---|---|---|---|---|---|---|---|
| *Acinetobacter baumannii* | 7 | 2 | 7 | 84 | 3 | 40 | 2,969 | 0 |
| *Campylobacter coli* | 8 | 5 | 8 | 94 | 0 | 24 | 3,091 | 0 |
| *Campylobacter jejuni* | 4 | 2 | 4 | 48 | 0 | 14 | 2,030 | 0 |
| *Clostridioides difficile* | 4 | 0 | 4 | 48 | 0 | 5 | 746 | 0 |
| *Corynebacterium diphtheriae* | 3 | 3 | 3 | 36 | 0 | 15 | 708 | 0 |
| *Enterobacter cloacae* | 9 | 5 | 9 | 108 | 4 | 29 | 909 | 0 |
| *Enterococcus faecium* | 12 | 1 | 12 | 138 | 0 | 53 | 12,694 | 0 |
| *Escherichia coli* | 29 | 10 | 29 | 348 | 17 | 144 | 18,307 | 0 |
| *Klebsiella aerogenes* | 3 | 3 | 3 | 36 | 0 | 11 | 229 | 0 |
| *Klebsiella oxytoca* | 4 | 1 | 4 | 48 | 0 | 17 | 380 | 0 |
| *Klebsiella pneumoniae* | 32 | 9 | 32 | 382 | 18 | 143 | 21,108 | 2 |
| *Klebsiella quasipneumoniae* | 3 | 0 | 3 | 36 | 0 | 6 | 304 | 0 |
| *Mycobacterium tuberculosis* | 20 | 0 | 20 | 218 | 10 | 99 | 20,922 | 0 |
| *Neisseria gonorrhoeae* | 6 | 0 | 6 | 72 | 0 | 28 | 7,176 | 0 |
| *Pseudomonas aeruginosa* | 12 | 1 | 12 | 144 | 0 | 40 | 5,779 | 0 |
| *Salmonella enterica* | 25 | 12 | 25 | 291 | 0 | 118 | 23,386 | 0 |
| *Shigella flexneri* | 4 | 2 | 4 | 48 | 0 | 18 | 1,442 | 0 |
| *Shigella sonnei* | 20 | 2 | 20 | 236 | 0 | 66 | 15,295 | 0 |
| *Staphylococcus aureus* | 16 | 5 | 16 | 171 | 0 | 70 | 18,295 | 2 |
| *Staphylococcus haemolyticus* | 3 | 0 | 3 | 36 | 0 | 14 | 342 | 0 |
| *Streptococcus agalactiae* | 3 | 0 | 3 | 36 | 0 | 8 | 423 | 0 |
| *Streptococcus pneumoniae* | 8 | 2 | 8 | 96 | 0 | 22 | 14,329 | 0 |
| *Vibrio cholerae* | 3 | 0 | 3 | 36 | 0 | 7 | 239 | 0 |
| All | 238 | 65 | 238 | 2,790 | 52 | 991 | 171,103 | 4 |

**Breakpoints.** "Breakpoints" counts the datasets whose breakpoint in the rebuilt revision 14 table differs from the release report's. The run set each to the release's: 52 of 238 datasets.

**Lists.** "Lists" counts the BV-BRC CDS lists fetched again for the feature cards: 991. None is equal in bytes to an entry of the run's pinned file list. By the entry "Refit cards" that can move a card's `kmer` and `locus` and no fitted field.

**Cards.** "Cards" counts the feature cards whose `kmer` or `locus` differs from the release's: 4 cards in 3 datasets. They are *Klebsiella pneumoniae* with imipenem and with tetracycline, 1 card each, and *Staphylococcus aureus* with clindamycin, 2 cards.

**Samples.** "Samples" counts, summed over a species' datasets, the samples whose country, year, host or MLST type differs from the release report's: 171,103 in all. The check counts these fields apart, and they fail no dataset.

The 238 datasets by species, with the antibiotics in full:

*Acinetobacter baumannii* has 7 datasets: cefepime, cefotaxime, doripenem, doxycycline, minocycline, moxifloxacin, sulfamethoxazole/trimethoprim. D0 reads 2 datasets: doxycycline, minocycline.

*Campylobacter coli* has 8 datasets: azithromycin, ciprofloxacin, clindamycin, erythromycin, gentamicin, nalidixic acid, telithromycin, tetracycline. D0 reads 5 datasets: azithromycin, ciprofloxacin, erythromycin, gentamicin, nalidixic acid.

*Campylobacter jejuni* has 4 datasets: ciprofloxacin, clindamycin, nalidixic acid, tetracycline. D0 reads 2 datasets: ciprofloxacin, nalidixic acid.

*Clostridioides difficile* has 4 datasets: azithromycin, ceftriaxone, clarithromycin, moxifloxacin. D0 reads none of them.

*Corynebacterium diphtheriae* has 3 datasets: sulfamethoxazole/trimethoprim, tetracycline, trimethoprim. D0 reads 3 datasets: sulfamethoxazole/trimethoprim, tetracycline, trimethoprim.

*Enterobacter cloacae* has 9 datasets: cefepime, ciprofloxacin, gentamicin, imipenem, levofloxacin, meropenem, sulfamethoxazole/trimethoprim, tobramycin, trimethoprim. D0 reads 5 datasets: gentamicin, imipenem, meropenem, sulfamethoxazole/trimethoprim, trimethoprim.

*Enterococcus faecium* has 12 datasets: ampicillin, chloramphenicol, clindamycin, dalfopristin/quinupristin, doxycycline, erythromycin, gentamicin, kanamycin, streptomycin, teicoplanin, tetracycline, vancomycin. D0 reads 1 dataset: doxycycline.

*Escherichia coli* has 29 datasets: amikacin, ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefotaxime, cefoxitin, cefpodoxime, cefpodoxime/clavulanate, ceftazidime, ceftiofur, ceftriaxone, cefuroxime, cephalothin, chloramphenicol, colistin, doxycycline, ertapenem, erythromycin, imipenem, levofloxacin, nalidixic acid, nitrofurantoin, piperacillin/tazobactam, streptomycin, sulfamethoxazole, tetracycline, tobramycin, trimethoprim. D0 reads 10 datasets: amikacin, ampicillin/sulbactam, aztreonam, cefazolin, cefepime, ceftiofur, doxycycline, ertapenem, nitrofurantoin, trimethoprim.

*Klebsiella aerogenes* has 3 datasets: ertapenem, imipenem, meropenem. D0 reads 3 datasets: ertapenem, imipenem, meropenem.

*Klebsiella oxytoca* has 4 datasets: cefotaxime, cefoxitin, cefpodoxime, ertapenem. D0 reads 1 dataset: ertapenem.

*Klebsiella pneumoniae* has 32 datasets: ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefotetan, cefoxitin, cefpodoxime, cefpodoxime/ceftriaxone, cefpodoxime/clavulanate, ceftazidime, ceftazidime/avibactam, ceftriaxone, cefuroxime, chloramphenicol, ciprofloxacin, colistin, doripenem, doxycycline, ertapenem, fosfomycin, imipenem, levofloxacin, nitrofurantoin, ofloxacin, polymyxin b, streptomycin, sulfamethoxazole/trimethoprim, tetracycline, ticarcillin/clavulanate, tigecycline, tobramycin, trimethoprim. D0 reads 9 datasets: ampicillin/sulbactam, aztreonam, cefepime, cefotetan, doripenem, doxycycline, ofloxacin, streptomycin, ticarcillin/clavulanate.

*Klebsiella quasipneumoniae* has 3 datasets: cefotaxime, ceftazidime, sulfamethoxazole/trimethoprim. D0 reads none of them.

*Mycobacterium tuberculosis* has 20 datasets: amikacin, amoxicillin, capreomycin, ciprofloxacin, cycloserine, delamanid, ethambutol, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, nicotinamide, ofloxacin, para-aminosalicylic acid, prothionamide, pyrazinamide, rifabutin, rifampin, streptomycin. D0 reads none of them.

*Neisseria gonorrhoeae* has 6 datasets: cefixime, ceftriaxone, ciprofloxacin, erythromycin, penicillin, tetracycline. D0 reads none of them.

*Pseudomonas aeruginosa* has 12 datasets: amikacin, aztreonam, cefepime, ceftazidime, ceftazidime/avibactam, ceftolozane/tazobactam, gentamicin, imipenem, levofloxacin, meropenem, piperacillin/tazobactam, tobramycin. D0 reads 1 dataset: amikacin.

*Salmonella enterica* has 25 datasets: amoxicillin/clavulanate, ampicillin, ampicillin/sulbactam, azithromycin, cefazolin, cefotaxime, cefoxitin, ceftazidime, ceftiofur, ceftriaxone, chloramphenicol, ciprofloxacin, gentamicin, kanamycin, nalidixic acid, pefloxacin, spectinomycin, streptomycin, sulfamethoxazole, sulfamethoxazole/trimethoprim, sulfathiazole, sulfisoxazole, temocillin, tetracycline, trimethoprim. D0 reads 12 datasets: ampicillin, cefotaxime, ceftazidime, chloramphenicol, gentamicin, kanamycin, sulfamethoxazole, sulfamethoxazole/trimethoprim, sulfathiazole, sulfisoxazole, tetracycline, trimethoprim.

*Shigella flexneri* has 4 datasets: ampicillin, azithromycin, ciprofloxacin, trimethoprim. D0 reads 2 datasets: ciprofloxacin, trimethoprim.

*Shigella sonnei* has 20 datasets: amikacin, ampicillin, azithromycin, cefoperazone, cefotaxime, cefoxitin, cefpirome, ceftazidime, ceftriaxone, ciprofloxacin, colistin, ertapenem, gentamicin, nalidixic acid, streptomycin, sulfamethoxazole/trimethoprim, temocillin, tetracycline, tobramycin, trimethoprim. D0 reads 2 datasets: colistin, trimethoprim.

*Staphylococcus aureus* has 16 datasets: cefoxitin, ceftaroline, chloramphenicol, ciprofloxacin, clindamycin, daptomycin, erythromycin, fusidic acid, gentamicin, levofloxacin, moxifloxacin, mupirocin, penicillin, rifampin, sulfamethoxazole/trimethoprim, tetracycline. D0 reads 5 datasets: chloramphenicol, clindamycin, erythromycin, mupirocin, tetracycline.

*Staphylococcus haemolyticus* has 3 datasets: ciprofloxacin, fusidic acid, tetracycline. D0 reads none of them.

*Streptococcus agalactiae* has 3 datasets: clindamycin, erythromycin, tetracycline. D0 reads none of them.

*Streptococcus pneumoniae* has 8 datasets: azithromycin, cefuroxime, chloramphenicol, clindamycin, erythromycin, meropenem, sulfamethoxazole/trimethoprim, tetracycline. D0 reads 2 datasets: azithromycin, chloramphenicol.

*Vibrio cholerae* has 3 datasets: nalidixic acid, sulfamethoxazole, trimethoprim. D0 reads none of them.
