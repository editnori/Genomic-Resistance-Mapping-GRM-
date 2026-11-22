# Assumption plan: corrections

The plan in `docs/assumption-plan.md` is fixed, and its text does not change. This page lists what was found wrong in it after it was fixed, what a later design added to it, and what holds instead. Each entry gives its date. The first two were written on 22 Nov 2026, before any fit.

## Study split

- **Plan.** In Measured records, the Study split row says: "Study folds cost recall, specificity and log loss for 5 of 5 learners".
- **Found.** `docs/study-split.md` reads specificity and log loss as worse for 5 of 5 learners. It reads recall as worse for 2 of 5, the random forest and the boosted trees. For the other 3 learners recall is lower and not shown.
- **Holds.** The row is read as its source has it. No rule of the plan depends on this row.

## Fold models

- **Plan.** D0 "reads the release's fold models" and "fits nothing". C8's first arm is "the release's rule from its saved `fold_models`, with no fit", and its Code row says both inputs exist "in every report".
- **Found.** 0 of the 274 reports of the `data-rev14-274-20261112` release hold `fold_models`. A report holds them only when `grm learn` computes it again.
- **Holds.** C8's fold overlap waits for a refit of the release's rule that saves fold models. So do D0's Lineage and Unused causes. That refit is the plan's default arm, named as a refit. It changes no arm, dataset, measure or cut point.
- **Unchanged.** C8's count of chosen features whose notes carry a catalogue family reads the release's final models and needs no fit.
- **Cost.** D0's cost of about 4 xlarge-hours assumed the AMRFinderPlus program's per-genome tables were at hand. They are rebuilt for each species. One run on *Klebsiella aerogenes* took 13 min 32 s, and from it the other 15 species are put at about 18.4 xlarge-hours, $3.67.

## MIC design

- **Plan.** C6 part B names three arms: "the released censored MIC model, now on lineage and study folds; an ordered model over dilution steps; the training median". It does not say how the arms are built.
- **Design.** A design fixed the rest at 00:13:20Z on 22 Nov 2026, before the code ran on real data. Its sha256 is `d6b3c92ae30ef87dedf0b6deca36f0dd57da802698d91baa3c6c5b589e567f3c`. The plan holds none of the eight items below.
- **Timing.** One run followed at 10:32Z on 22 Nov 2026: *Campylobacter coli* with ciprofloxacin (`campylobacter_coli__ciprofloxacin`), seed 42, the released and the ordered model. No score was read from it.
- **Figures.** That run took 2.23 s and 0.78 GiB. Its input fingerprint equals the release report's.
- **Written.** This entry was written at 10:58Z on 22 Nov 2026, after that run and before the full run. No score of any arm had been read.

The eight items:

1. **Study plans.** Seeds 42 to 46 take the five study assignments of `docs/study-split.md`, unchanged. Seed 47 takes a sixth, placed by that page's rule B on the units all five keep in one fold.
2. **Cohort.** Lineage folds are drawn over all R and S genomes of a dataset. Each arm trains on the other folds' genomes with an MIC and scores the fold's own. Every arm leaves out a fold that tests no genome with an MIC.
3. **Steps.** The ordered model puts an exact MIC of v, and an MIC of at most v, at step v. It puts an MIC above v at step v + 1. Its thresholds run from the training fold's lowest step to one below its highest.
4. **Fits.** Threshold j of fold f is fitted with seed `seed + f + 1000·(j + 1)`. Its 15% holdout is stratified on the threshold's target. Its screen is the classifiers' 2×2 chi-square, top 4,000.
5. **Prior.** Each threshold's balanced fit moves back to its training share π: q = sigmoid(margin + logit π). A threshold takes q = π when its deciding genomes hold one class, or its fitting part or holdout lacks one.
6. **Gaps.** A threshold that no training genome decides takes q from the nearest lower decided threshold, else from the nearest higher. A fold with no decided threshold leaves the ordered model unavailable.
7. **Read-off.** A genome's q values are made non-increasing by pool-adjacent-violators with equal weights. The predicted step is the lowest step plus the number of thresholds whose q is at least 0.5.
8. **Median.** The training median is also written at seeds 45, 46 and 47. The plan's twin at those seeds is the released model only.

- **Holds.** Part B's release rule reads the released model's lineage EA-exact against the training median's. No item changes that reading.
- **Items.** Items 1 and 8 add rows that decide nothing. Item 2 gives the three arms the same test genomes. Items 3 to 7 define the ordered model only.
- **Seed 47.** The sixth study plan is coarser: its units merge study groups on 49 of 70 datasets.
- **Seeds.** On the study split `grm learn` refuses `--seed` beside `--fold-file`. So the fits keep the release's seeds, and only the study plan changes from seed to seed.
- **Runs.** The first run takes the datasets that are not M: 100 on the lineage split and 29 on the study split, each at 6 seeds. That is 774 runs on 12 species.
- **Cost.** The plan put part B at 15 xlarge-hours at most. From the timing run and the release's fit times, the first run is put at about 13.1 xlarge-hours, $2.63.

## MIC checks

The MIC design says each result's input fingerprint equals its release report's, and it keys the seed-47 study plan on "the unit's lowest genome id". Two species ran first, on 22 Nov 2026. *Campylobacter coli* gave 48 results on 8 datasets. *Acinetobacter baumannii* gave 30 results on 4 datasets, 1 of them also on the study split. That is 78 of the 774 results. The rehearsal met two cases where the fingerprint cannot be equal, and a review of the code named two points the design leaves open. Four changes follow.

**Key.** The seed-47 plan compares the two numbers of a genome id as numbers, never as text. The design leaves this open. The code's first draft did the same at 00:37Z on 22 Nov 2026, before any fit.

**Seeds.** The fingerprint hashes the seed, and the release holds reports at seed 42 only. A lineage result at seeds 43 to 47 must hold the report's settings but the seed, and the seed-42 result's cohort row by row. The cohort is the result's counts and each tested genome's MIC, censor and call.

**Breakpoints.** The fingerprint hashes a dataset's breakpoint, fitted on every MIC pair of the species and antibiotic. The rebuilt revision 14 table passes 79 more genomes through QC. In 33 of the 100 lineage datasets the breakpoint gains 1 to 19 MIC pairs. The R cut-off stays the same in all 33, and none of the 79 genomes is in a species' genome list.

**Binary.** The rehearsal's results count as the run's own when the run's binary has the code sha256 `20e08189285ece9ac29a03f989469899356d74e5876927e689f16a1a037c104c` and the same input files. The `grm` file's own sha256 changes with its build path.

The cohort of each of the 274 datasets equals its release report's, genome for genome and call for call.

Each of the 33 datasets is fitted again at seed 42 with the release's breakpoint. That fit must give the release's fingerprint, and it must equal the run's seed-42 result in every field but the fingerprint and the wall time. The other seeds, and the study plans of the 17 probed datasets on the study split, rest on a reading of the code: no fit reads the dataset's breakpoint, and each fold fits its own.

| Species | Antibiotics | Datasets |
|---|---|---|
| *Acinetobacter baumannii* | cefotaxime, doripenem, sulfamethoxazole/trimethoprim | 3 |
| *Escherichia coli* | amikacin, aztreonam, cefepime, cefotaxime, ceftazidime, piperacillin/tazobactam, tetracycline | 7 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefoxitin, ceftazidime, ciprofloxacin, ertapenem, imipenem, nitrofurantoin, sulfamethoxazole/trimethoprim, tetracycline, tobramycin | 13 |
| *Mycobacterium tuberculosis* | amikacin, delamanid, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, para-aminosalicylic acid, rifabutin, rifampin | 10 |

No change touches what part B measures: the models, the folds, the measures and their denominators are as designed. One thing got worse. At seeds 43 to 47 no check compares the k-mer matrix with the release's; the run builds one matrix a species, and the seed-42 check covers it. This entry was written at 13:06Z on 22 Nov 2026, after the rehearsal and the review and before any other species ran. No score of any model had been read.

## MIC list

The MIC design says `scripts/c6-mic/study_plans.py` writes the seed-47 study plan `c6mic-s47` and checks it: both classes in each fold's test and training genomes. It does not say what follows when a plan fails its check. One of the 29 study datasets fails: *Klebsiella pneumoniae* with ampicillin/sulbactam (`klebsiella_pneumoniae__ampicillin_sulbactam`), 744 genomes in 14 units. Its plan ends with "fold 4's test genomes lack a class". The other 28 plans pass.

**Change.** That dataset leaves the study list of *Klebsiella pneumoniae*, which goes from 14 datasets to 13. It stays on the lineage list with its 6 lineage results. The run holds 768 results where the design counted 774: 600 on the lineage split, as before, and 168 on the study split in place of 174.

**Loss.** The dataset's 6 study results are not run. Its plans at seeds 42 to 46 pass their checks, so 5 of the 6 could have run. The scripts take one study list a species and run every seed on it, and they stay as reviewed. Part B's release rule reads the lineage split, so no reading of that rule loses a result.

**Timing.** *Klebsiella pneumoniae* stopped at the plans step at 13:48Z on 22 Nov 2026, before any of its fits, and started again on the list of 13 at 13:59Z. This entry was written at 14:01Z on 22 Nov 2026. Nine species had finished and two more were running. No score of any model had been read.

## D0 causes

The plan's D0 puts each missed lineage fold in one of four causes, Absent, Removed, Lineage and Unused, and reads the release's fold models. The plan does not say when a fold counts as missed or how a cause is decided, and the release holds no fold models (entry "Fold models"). A review of D0's and C8's code on 22 Nov 2026 named eleven readings the plan leaves open. Each holds from 22 Nov 2026.

| Reading | Rule |
|---|---|
| Missed fold | D0 counts a lineage fold as missed when it tests both classes and the program's fold balanced accuracy beats the learner's by 0.10 or more. A miss pair with no such fold takes, of its folds that test both classes, the one with the largest gap and flags it `fallback`, but only when that gap is above 0. Otherwise the pair has no missed fold and no cause: `combine` counts it as `no_missed_fold`, and it stays in the miss pairs, the 40% denominator and the gap total. |
| Match | A feature matches the program, and a chosen feature agrees with a match, when its presence or its absence agrees at 0.95 or more on the fold's training genomes the program calls, its no-call genomes left out. |
| Chosen | A missed fold where the learner chose a feature that agrees with an offered match but lies in more than two lineage clusters goes to a fifth cause, Chosen, which counts in the 40% denominator and never sets the order. |
| Order | Each missed fold takes the first cause that holds in the order Absent, Removed, Lineage, Chosen, Unused, and the 40% is taken of miss pairs (dataset and learner), with a tie for most setting no order. |
| Ties | A dataset whose missed folds tie between causes counts in each and splits its gap equally between them. |
| Trees | For the random forest and the boosted trees a fold "chose" its 50 heaviest features, all `fold_models` keeps, so a feature used below that rank counts as not chosen in D0 and in C8. |
| Clusters | A chosen feature's lineage clusters are counted over the fold's training genomes that carry its resistance-indicating state, and a genome with no lineage counts as a cluster of its own. |
| Offered | A matching row counts as offered when the lowest row of its training pattern is a candidate, as `grm learn` sees one row per pattern. |
| Empty pairs | C8's median Jaccard leaves out a fold pair in which both folds chose nothing. |
| Pattern overlap | C8 also reports a second overlap in which two features are one when their patterns agree on every sample both folds train on. |
| Catalogue family | C8's count of features whose notes carry a catalogue family reads the final models' features that have a card, since fold features have none. Its row over all learners is named `any_summed`: it counts a dataset once and sums the five learners' features. |

**Chosen.** The lead ruled at 15:50Z on 22 Nov 2026 that the fifth cause stays. A missed fold where the learner chose a feature that agrees with an offered match and lies in more than two lineage clusters is neither Lineage nor Unused by the plan's words.

**Default.** Chosen counts in the 40% denominator and never sets the order. The largest cause is taken among Absent, Removed, Lineage and Unused. It sets the order when it alone is the largest of the four and holds 40% or more of the miss pairs. Otherwise the plan's default order C4, C3, C5 holds.

**Reading.** A pair whose missed folds tie counts in each tied cause, so Chosen can hold more miss pairs than a cause of the four that holds 40% or more. That cause still sets the order. The lead confirmed this reading at 17:05Z on 22 Nov 2026, with no D0 count read. `d0.py` already does this, and no code changed for it.

**Refit.** D0 and C8 read a refit of the release's rule at seed 42 that saves fold models, on the 238 datasets that are not M. The build's code sha256 must be `20e08189285ece9ac29a03f989469899356d74e5876927e689f16a1a037c104c`. A species with any difference from the release stops before D0 or C8 reads a report.

**Equality.** A refit report stands for the release's only when `scripts/refit/refit_check.py` finds every field of the release's report but `seconds` equal: the top-level fields, each sample's fitted fields, each split's folds, and each learner row but its `seconds`.

**Allowed.** Three new fields pass: `fold_models`, complete for the five learners on every split that ran, each sample's `rule_groups`, and `selection`. Any other new field fails. A sample's country, year, host and MLST type are counted apart and do not fail it. With the allowed fields taken out, the two reports written again in their own key order must be the same text.

**Written.** This entry was written at 16:08Z on 22 Nov 2026, while the refit, D0 and C8 were rehearsed on *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, for time and equality only. No D0 or C8 count beyond the rehearsal on *Klebsiella aerogenes* had been read. A review made the first reading and the default's tie case exact at 16:40Z the same day, and by then still no D0 or C8 count beyond that rehearsal had been read.

## C1 design

The plan's C1 part B draws cohorts of r resistant and r susceptible genomes from the BIG datasets and counts how far each cohort's estimate lies from the truth. It leaves open how a cohort is drawn, which interval's coverage counts and which row decides. A design fixed those on 22 Nov 2026, before its code ran on real data. A review of the code named twelve deviations from the plan, each holding from that day.

| Item | Deviation |
|---|---|
| Draw | Draws take r genomes per class by a keyed round robin over lineages, nested across r; the plan says only "r per class within lineage groups". |
| Interval | C1's release row reads coverage on a levelled interval from pooled recall and specificity; the five lineage folds' range is written beside it and decides nothing. |
| Scored set | The truth is scored on the release dataset's other cohort rows inside the release's 2,000-genome cap, not on every labelled genome of the species and antibiotic. |
| Copies | Genomes k-mer-identical to a cohort genome (same kset digest) are left out of the truth's scored set and counted by class; the plan does not name copies. |
| Feature value | A cohort feature's value on a scored genome is the species-matrix bit of its lowest species pattern, one choice among the patterns the cohort's projection merged. |
| Twin | The twin is r 25 at fold seed 45 on the same five cohorts, not Rules' default at seeds 45, 46 and 47. |
| Seeds | Each line uses fold seed 42 with five cohort draws in place of Rules' fold seeds 42, 43 and 44. |
| Code | The subsample arm sits on `grm learn --exp-cohort`, not on `grm arm` as the plan's Code column says. |
| First run | The first run is the 89 BIG datasets outside M in 10 species, and the 19 in M wait for the maintainer's label rule. |
| Breakpoints | The run gives the listed datasets the release's breakpoints, which the rebuilt table moves on 36 of the first run's 89 and no fit reads. |
| Pairing | A dataset's difference from 25 is its mean error difference over shared available draws; it enters on 3 or more, and coverage counts those shared results at both lines. |
| Cost | The first run is costed at 7.80 xlarge-hours central and 18.86 upper, against the plan's about 4.5. |

**Draw.** Both classes walk one lineage order, and the draw refuses a class with r or fewer genomes.

**Rulings.** The lead ruled three design items at 17:14Z on 22 Nov 2026, with no score of any cohort run read by anyone. All three sit in the scorer, `scripts/c1-cutoff/score.py`, and no fit changes. They are the rows Interval and Pairing above, and the row that decides.

**Interval.** Every row counts coverage on two intervals, each with its median width: the fold range, which is what a report shows, and the levelled interval. With recall p_R in its Wilson 95% interval (l_R, u_R) and specificity p_S in (l_S, u_S), the levelled interval runs from b − ½·√((p_R − l_R)² + (p_S − l_S)²) to b + ½·√((u_R − p_R)² + (u_S − p_S)²), where b = (p_R + p_S) / 2. An edge counts as covered.

**Level.** The release row compares coverage between lines. That comparison means something only when the interval has the same level at each line. The fold range has no nominal level and widens as r falls, so its coverage rises with its width. Each split builds its levelled interval from its own pooled counts.

**Decides.** The row that decides is the lineage split with the truth on all scored genomes. A learner meets the release row when the median over datasets of its error minus 25's has a species interval wholly below 0.02 and its levelled coverage is within 0.05 of 25's. A line is a candidate when at least 3 of the 5 learners meet it, the form Rules' condition 1 takes for a change that touches all learners.

**Others.** `count.tsv` marks that one row per line as the row that decides. Every other row is written and decides nothing: the truth on new lineages, the random and the study split, each learner alone, the fold range's coverage and the twin.

**Pairs.** Draw d at a line lies inside draw d at 25, so the count pairs the two lines on their shared draws: the error at the line less the error at 25, draw by draw. Each line's own available draws and its coverage over all its own results are written beside. A dataset with no line at 25 is counted as not scored.

**Resamples.** Each row's species resamples take their seed from the sha256 of its split, truth set, learner and line, so adding a row moves no other row's interval.

The 36 first-run datasets whose breakpoint the rebuilt table moves:

| Species | Antibiotics | Datasets |
|---|---|---|
| *Escherichia coli* | aztreonam, cefepime, cefotaxime, cefoxitin, ceftazidime, ceftriaxone, cefuroxime, ertapenem, piperacillin/tazobactam, tetracycline, tobramycin, trimethoprim | 12 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefoxitin, ceftazidime, ceftriaxone, ciprofloxacin, ertapenem, imipenem, levofloxacin, sulfamethoxazole/trimethoprim, tetracycline, tobramycin | 14 |
| *Mycobacterium tuberculosis* | amikacin, delamanid, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, para-aminosalicylic acid, rifabutin, rifampin | 10 |

The 19 BIG datasets in M, which wait for the maintainer's label rule:

| Species | Antibiotics | Datasets |
|---|---|---|
| *Acinetobacter baumannii* | amikacin, ampicillin/sulbactam, imipenem, meropenem, tetracycline, tobramycin | 6 |
| *Escherichia coli* | amoxicillin/clavulanate, ampicillin, ciprofloxacin, gentamicin, sulfamethoxazole/trimethoprim | 5 |
| *Klebsiella pneumoniae* | amikacin, gentamicin, meropenem, piperacillin/tazobactam | 4 |
| *Mycobacterium tuberculosis* | isoniazid | 1 |
| *Neisseria gonorrhoeae* | azithromycin | 1 |
| *Pseudomonas aeruginosa* | ciprofloxacin | 1 |
| *Staphylococcus aureus* | oxacillin | 1 |

**Cost.** The plan put part B at about 4.5 xlarge-hours, $0.90. A pilot on *Streptococcus pneumoniae* with erythromycin and the release's dataset sizes put the first run at 7.80 xlarge-hours central, $1.56, and 18.86 upper, $3.77. Fits make 4.83 of the 7.80, against the plan's 0.91, and the genome fetch 2.46, which the plan does not cost.

**Fetch.** The run fetches genomes 8 at a time (`--workers 8`, `grm fetch fasta`'s own default) in place of 4. On the review's 16-CPU machine, the genome step of *Escherichia coli* with aztreonam, which fetches the 2,000 release genomes and checks each, took 389.3 s at 8 against 615.5 s at 4, with no retry logged.

**Cuts.** With the pilot's genome step scaled by that ratio, the first run costs 6.90 xlarge-hours central, $1.38, and 17.51 upper, $3.50. The other cut, one species' fetch overlapping the previous species' fits, is left out.

**Written.** This entry was written at 17:56Z on 22 Nov 2026, after the review and the lead's rulings and before the first run. No score of any cohort run had been read.

## Refit cards

*Staphylococcus aureus* stopped at the refit's check at 17:42Z on 22 Nov 2026, `FAILED staphylococcus_aureus at check`. Of its 16 datasets, 15 were equal to the release in every fitted field and all 16 held complete fold models. All 171 learner rows were equal. The one field that differed was `features` of *Staphylococcus aureus* with clindamycin (`staphylococcus_aureus__clindamycin`).

**Cards.** Of that dataset's 107 feature cards, 105 are equal. The other 2 differ only in `kmer` and in `locus`: `feature`, `bits`, `n_kmers`, `n_resistant`, `n_susceptible`, `catalog_families` and `catalog_products` are equal on 107 of 107. For both cards the locus genome is the same in the two reports and its contig differs.

**Files.** For both cards the genome file fetched again differs in bytes from the release's list. Over the species, 834 of 2,000 fetched genome files differ in bytes from that list, while length, contig count and N50 equal the release's on 2,000 of 2,000. The 6 species that passed had such files too: *Klebsiella aerogenes* 27 of 93.

**Ruling.** The lead ruled at 17:50Z on 22 Nov 2026 that a refit report stands for the release's when every field but `seconds` is equal, with one more allowance: in the top-level `features` list, a card's `kmer` and `locus` may differ. Cards are compared in order, one to one, and every other card field must be equal. The count of cards whose `kmer` or `locus` differs is written per dataset.

**Reason.** Both fields are written after the fit, from genome files and CDS lists fetched again from BV-BRC. A file served today can differ from the one the release fetched. That moves where a pattern is shown and leaves the pattern, the fit and every field D0 and C8 read as they were. D0 and C8 read a card's `feature`, `bits` and `catalog_families`, and neither `kmer` nor `locus`.

**Fails.** A difference in any other card field, in the number or the order of the cards, or in a `kmer` or `locus` key outside `features` still fails the check.

**Restart.** A species stopped at the check goes on from its fitted reports on a later checkout when the built binary's sha256 is unchanged and the commits differ only under `scripts/refit/`, `docs/`, `tests/` and `.github/`. Any other changed path refuses, `src/` and `scripts/assumption-d0-c8/` among them. The check runs again and `binary.txt` gains a line `restart <commit> <time>`. A species with a D0 dataset now stops at the kit step without numpy 1.26.4.

**Written.** This entry was written at 17:59Z on 22 Nov 2026, before any species was started again. No D0 or C8 count had been read.

## C1 scorer

A second review recounted the small-cohort scorer from the rulings' words on 22 Nov 2026. Its own program gave the same value as `scripts/c1-cutoff/score.py` in 230,612 of 230,618 SYNTHETIC cells, in 2,838 of 2,838 cells of *Escherichia coli* with aztreonam, and in 28,920 of 28,920 cells of *Neisseria gonorrhoeae* with ciprofloxacin, penicillin and tetracycline. The three inputs are separate, and their cell counts do not add. The two real outputs were compared for equality only.

**Faults.** The review found three faults, fixed here: unrounded medians, an edge test on unrounded values, and tests that let planted faults pass.

**Median.** The differences were rounded to 12 decimals, and the median difference and the 4,000 resampled medians were left as computed. Differences of -0.108 and 0.148 have the median 0.02, which floating point gives as 0.019999999999999997. A species interval ending at exactly 0.02 then read as wholly below 0.02, and the learner met the release row.

**Fix.** Both medians are now rounded as the differences are. The 6 SYNTHETIC cells that differed are this fault, and with the fix 230,618 of 230,618 are equal.

**Edge.** An edge counts as covered, and the edge test compared the values as computed. Calls of 60 of 100 and 70 of 100 give a truth of 0.6499999999999999, and a fold of 13 of 20 in each class gives 0.65. A truth equal to the lowest fold then read as outside the fold range. The truth and both ends are now compared at the 12 decimals the files write, for every interval the scorer counts.

**Tests.** Of the 9 faults the review planted in the count, 3 failed a test. New SYNTHETIC tests pin the margin at -0.02, entry on 3 shared draws, the mean over draws, and the species interval's 101st and 3,900th values. `tests/c1_cutoff` goes from 44 tests to 50. With them, 12 of 12 planted faults fail a test: those 9 and 3 more that undo these fixes.

**Rulings.** The review listed seven things the scorer does that no ruling covered. The lead ruled on them at 19:01Z on 22 Nov 2026, with no score read. Each stands as the code has it.

| Item | Ruling |
|---|---|
| Random split | Its rows are written and decide nothing, as the plan's Rules show that split. |
| Coverage | In the release row, coverage at a line and at 25 runs over the shared results of the datasets that enter, those with 3 or more shared draws. The difference and the coverage then speak of the same datasets. |
| Decimals | A difference, a median and an interval's end are compared at the 12 decimals the files write. Coverage is compared with 25's in exact fractions. |
| Fold range | A result with a fold that has no balanced accuracy has no fold range and leaves the fold range's count. That coverage can run over different results at a line and at 25. It decides nothing and is written with its count of results. |
| Wilson | The Wilson interval is clipped to 0 and 1, which removes float noise only. |
| Estimate | The error uses the balanced accuracy `grm` wrote, and the levelled interval is centred on b from the calls. A result where the two differ by more than 1e-9 fails its check. |
| No dataset | A learner with no dataset entered at a line does not meet the release row. The twin's row is written and is never a candidate. |

**Cost.** The "C1 design" table's Cost row gives 7.80 and 18.86 xlarge-hours, the cost at 4 downloads at once. The run fetches 8 at once, and its cost is the Cuts paragraph's: 6.90 xlarge-hours central and 17.51 upper.

**Run.** The first run's 10 species started at `9945b97`, before these fixes, and their score step ran with the earlier scorer. Each species is scored again from its fitted results at the commit that holds these fixes. The counts are read from that second scoring alone.

**Written.** This entry was written at 19:01Z on 22 Nov 2026. No score of any cohort run had been read.

## C1 reading

The small-cohort run fitted its 10 species and scored each a second time by 19:55Z on 22 Nov 2026. A reading program, `scripts/c1-cutoff/read.py`, turns the second scoring into the result page's tables. Its writer left nine questions on how the count is read. The lead ruled on them at 20:03Z on 22 Nov 2026, with no score read.

| Question | Answer |
|---|---|
| Measures | Rules' measures run over the shared draws of the datasets that enter, with 25 beside on the same draws, as the count does. |
| Line 50 | The release row names a lower line, and 50 lies above 25. Its row is written, decides nothing and is never a candidate. |
| Twin | Rules' condition 3 is written for a gain. C1's row asks that a line be no worse than 25 by 0.02, so a twin that meets it voids nothing. |
| Twin fails | When fewer than 3 learners meet on the twin, a line that is no candidate is written "not shown to be within the margin". |
| First run | The reading is on the 89 datasets outside M, and a candidate is written "a candidate on the 89". |
| Several lines | Each lower line is read on its own against 25. The page shows each reading and names no line as the one to take. |
| Pooled | Each of Rules' measures is given pooled, as a count over its genomes, with the median over datasets beside it. |
| Differences | Each dataset's mean difference is written by name, and the reading program checks that their median is the scorer's. The differences decide nothing. |
| Conditions | Rules' six conditions are for a change to the release and are not applied to C1. Conditions 2, 4 and 5 are not computed. |
| Breakpoints | The 36 datasets with a moved breakpoint are named in the page's limits, and no count is split by them. |

**Later.** The count over the 108 datasets of BIG is taken when the maintainer's label rule lets the 19 run, with the same rules, scorer and arguments. When the twin fails the row, no line is called worse than 25 on that row.

**Arguments.** The reading program runs once, with `--measures-over shared --line-50 beside --twin-rule beside --dataset-differences yes`.

**Cells.** The entry "C1 scorer" says the twin is never a candidate. The scorer writes its `candidate` cell on the twin's line row as on any line: 1 when 3 or more learners meet. It also marks line 50's row as the row that decides.

**Kept.** The scorer's count file keeps those cells, and the scorer stays as it scored the run. The reading program carries the twin's cell as `three_or_more_meet` and puts line 50's rows among those that decide nothing. No page calls the twin or line 50 a candidate.

**Value.** In the measures table, a dataset's value is the mean over its results of each result's rate, as Rules' Draws take a dataset's value as a mean. The pooled figure is the count over genomes.

**Part A.** Part A's counts were taken on 21 Nov 2026, before any cohort fit. At a cap of 2,000 genomes, lines 5, 10, 15, 20, 25 and 50 admit 457, 373, 332, 301, 277 and 212 datasets. The six counts are separate lines over one table and do not add. The result page carries that table beside each lower line's reading.

**Review.** A review recounted each dataset's entry and Rules' measures from the rulings' words, before the reading program ran on the 10 species. Its own program gave the same value in 44,160 of 44,700 SYNTHETIC cells and in 29,175 of 29,175 cells of *Neisseria gonorrhoeae* with ciprofloxacin, penicillin and tetracycline. The two inputs are separate, and their cell counts do not add. The real output was compared for equality only.

**Faults.** The 540 SYNTHETIC cells that differed are 60 rows the program did not write, for a learner with no available result at a line. The tables also lacked the twin's fixed words and the part of the release row a learner fails. Of the review's 12 planted faults, 8 passed every test.

**Fix.** Five commits fix these. The 60 rows are written with 0 results, and 44,700 of 44,700 cells are equal.

**Part.** A learner that does not meet reads as failing the interval's upper end, the coverage or both, or as having no dataset entered, from the scorer's cells. The program stops when those parts do not give the scorer's `meets` cell, and when a species' line file is not what its score file gives.

**Tests.** `tests/c1_cutoff` goes from 84 tests to 97, and 14 of 14 planted faults fail a test: the 12 and 2 planted in input files.

**Written.** The rulings were fixed at 20:03Z on 22 Nov 2026 and this entry was written at 20:44Z. No score of any cohort run had been read.

## Refit reading

The refit of the release's rule ended on its 23 species and 238 datasets at 20:36Z on 22 Nov 2026. A reading program, `scripts/assumption-d0-c8/read_refit.py`, turns the species' D0 and C8 tables into the tables of one page. Its writer left ten questions the plan does not answer. The lead ruled on them at 17:52Z on 22 Nov 2026, before the run ended. No D0 or C8 count had been read beyond the rehearsal on *Klebsiella aerogenes* that the entry "D0 causes" names.

| Choice | Fixed |
|---|---|
| One fold | In C8, a unit where one fold alone chose a feature keeps its four pairs at overlap 0, so its median is 0. Only a unit where no fold chose is undefined. |
| Five learners | C8's row over the five learners pools the dataset and learner values. It decides nothing. |
| Bands | C8 gets no cut point beyond the plan's "within 0.05". The quartiles stand beside the median. |
| Check tables | The dataset and learner check tables gate the reading. The breakpoint and feature-list tables are counted on the species table. |
| Order's rest | The cause that sets the order moves its block to the front. The other two blocks keep the default's order. |
| Coverage | The order is set now, on 177 of the plan's 204 miss pairs, in 65 of its 75 MISS datasets. |
| Family count | C8's family count reads the fold-chosen features that have a card, and is given with both denominators. |
| Per dataset | D0 per dataset and learner goes into a file beside the page, with species and antibiotic in full. |
| Removed | Removed keeps the plan's words, "C4's floor arm". The page says it covers the support floor and the 4,000-row screen together. |
| Pooled share | The gap share pooled over the five learners is not on the page. |

**Stability.** The other 27 miss pairs lie in the 10 MISS datasets marked M. Take x as the leading cause's pairs and y as the next cause's, among the four. A named first block is stable when x over 204 is 0.40 or more and x is above y plus 27. The default order is stable when no cause of the four reaches 0.40 of 204 with 27 pairs added. Otherwise the page says the 27 pairs could change the first block.

**Again.** The line is read again when the maintainer's label rule is made.

**Units.** Beside C8's median the page gives the units where fewer than two folds chose. The first and third quartiles use the inclusive method of Python's `statistics.quantiles`.

**Cards.** Under the entry "Refit cards", a dataset whose only difference from the release is the `kmer` or `locus` of a feature card reads equal. The species table counts such cards.

**Arguments.** The reading program runs once, with `--equal-tables datasets.tsv,learners.tsv --c8-one-fold zero --c8-over-five pairs --c8-bands none --order-rest default-after`. Its gates are given 23 species, 238 datasets, and D0's 65 datasets in 16 species with 177 miss pairs: 23 for `scm`, 28 for `cart`, 52 for `rf`, 41 for `gbt` and 33 for `lasso`.

**Review.** A review recounted D0's cause counts, the order line and C8's medians from the rulings' words, before the program ran on the 23 species. Its own program gave the same value in 1,408 of 1,413 cells of a small SYNTHETIC input and in 18,449 of 18,454 cells of each of 11 SYNTHETIC inputs with the run's counts. On *Klebsiella aerogenes* with ertapenem, imipenem and meropenem it gave the same value in 545 of 550 cells. The inputs are separate, and their cell counts do not add. The real output was compared for equality only.

**Faults.** The 5 cells that differed in each input were the gap shares pooled over the five learners, written on a page table against the last ruling. The tables also lacked the stability words and the count of units where fewer than two folds chose. Three gates had no test.

**Fix.** Seven commits fix these. After them the recount agrees in 1,426 of 1,426 cells, in 18,467 of 18,467 cells of each of the 11, and in 562 of 562 cells.

**Tests.** `tests/assumption_d0_c8` goes from 41 tests to 87. Of 24 planted faults, 14 in the code and 10 in input files, each fails a test or stops at a gate.

**Written.** The rulings were fixed at 17:52Z on 22 Nov 2026 and this entry was written at 21:25Z. The run's own check had compared each species' refit with the release. No D0 or C8 count had been read beyond that rehearsal.

## C9 design

**Block.** Plan block C9 asks whether the release's served learner is the right one: the best random-split balanced accuracy among the five learners, a tie going to the later learner. The design was fixed at 19:19Z on 22 Nov 2026, before any code ran on real data, and `docs/served-pick.md` gives the method. No comparison result has been read.

**Arms.** The default serves the release's pick on every outer lineage fold, read from the same run's random split at that run's seed. The inner pick serves, per outer fold, the learner with the best balanced accuracy on inner folds of that fold's training genomes. Both arms are scored with the release's own outer fits, so the default's pooled figure is its report's lineage row.

**Partition.** The inner folds are 3 folds of the outer training genomes, each lineage whole, drawn at the fold's own seed. That is the partition the release's SCM and CART already tune on there. Each inner fit is fitted as an outer fold is, with draws seeded from the run's seed, the outer fold and the inner fold.

**Measure.** The primary measure is the pooled lineage balanced accuracy of the served learner on the outer test folds. Recall, specificity, both error rates with Wilson intervals, and log loss stand beside it. A dataset's value is the mean of seeds 42, 43 and 44. The twin is the default at seeds 45, 46 and 47.

**Study.** The 55 first-run datasets in S get one study run each, on their reviewed `study-whole` plans, at the release's seeds. `grm learn` refuses a seed beside a fold file, so the study split has no twin. Its default is the seed-42 lineage run's pick.

**Count.** The count is the plan's Rules: better, worse or the same by 0.02 at 12 decimals, the median difference and its species interval, and the one-sided sign test. Holm runs over two pairs on the lineage primary measure, the inner pick and the twin.

| Set | Datasets | Species |
|---|---:|---:|
| L | 268 of 274 | 23 |
| L not in M, the first run | 233 | 23 |
| L in M, waiting for the label rule | 35 | 9 |
| First run in S | 55 of 233 | |
| First run with an S and R genome pair | 32 of 233 | |

**Labels.** On the 233, three of the four label rules give the same labels. The 32 with an S and R pair depend on the most resistant call against the other three, and the other 201 on none of the four.

**Questions.** The design leaves eight questions to the lead, each with the answer taken and the other one. The scorer runs either answer from the same results for questions 3 to 5.

| Question | Taken | Other |
|---|---|---|
| 1. Inner folds | 3, at the release's tuning seed | A partition at a seed of its own |
| 2. Inner measure | Balanced accuracy of the summed inner calls | Mean of the inner folds' figures |
| 3. Ties | The later learner, ties listed | The earlier, or the default's learner when tied |
| 4. No pick | The default's learner on that fold, counted | The draw left unscored |
| 5. Default | Each draw's own random-split pick | The seed-42 pick at every draw |
| 6. Study draws | One run at the release's seeds, no twin | A seed under this flag only |
| 7. Holm | Two pairs, condition 1 on the served figure | Per-learner rows |
| 8. First run | L less M | All of L |

**Pilot.** One draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, 93 genomes, ran in 108 s on a 16-CPU box at 8 threads. The fit with the inner pick took 30.85 s at a peak resident size of 925,212 KB, and the run used 5.09 GB of disk. Its seed-42 reports equal the release's in all 3 datasets.

**Cost.** The inner pick took 0.19 to 0.44 of each dataset's release learning seconds. At 0.44, the first run costs 18.12 xlarge-hours central. Of those, 3.47 go to fetching genomes and building matrices again on each box, which the plan's 17 for all 268 leaves out. At the plan's 1.5 it costs 24.03.

**Written.** This entry was written at 19:54Z on 22 Nov 2026, from the design fixed at 19:19Z, the lists and the pilot. No score of any C9 run had been read.
