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

**Arms.** The default serves the release's pick on every outer lineage fold, read from the same run's random split at that run's seed. The inner pick serves, per outer fold, the learner with the best balanced accuracy on inner folds of that fold's training genomes. Both arms are scored with the release's own outer fits, so the default's pooled figure is its report's lineage row. The default's pick reads the random split over every genome, the outer test folds' genomes included, so the default's served figure can carry the pick's optimism, and the inner pick's cannot. Wherever the two are compared, a same or a worse is never read as showing that the release's pick is right (ruling 9).

**Partition.** The inner folds are 3 folds of the outer training genomes, each lineage whole, drawn at the fold's own seed. That is the partition the release's SCM and CART already tune on there. Each inner fit is fitted as an outer fold is, with its own inner folds, tuning holdouts and draws seeded from the run's seed, the outer fold and the inner fold.

**Measure.** The primary measure is the pooled lineage balanced accuracy of the served learner on the outer test folds. Recall, specificity, both error rates with Wilson intervals, and log loss stand beside it. A dataset's value is the mean of seeds 42, 43 and 44. The twin is the default at seeds 45, 46 and 47.

**Study.** The 55 first-run datasets in S get one study run each, on their reviewed `study-whole` plans, at the release's seeds. `grm learn` refuses a seed beside a fold file, so the study split has no twin. Its default is the seed-42 lineage run's pick. One study run per dataset and no study twin is a deviation from the plan's Rules, Draws, which give each arm seeds 42, 43 and 44 (question 6).

**Count.** The count is the plan's Rules: better, worse or the same by 0.02 at 12 decimals, the median difference and its species interval, and the one-sided sign test. Holm runs over two pairs on the lineage primary measure, the inner pick and the twin.

| Set | Datasets | Species |
|---|---:|---:|
| L | 268 of 274 | 23 |
| L not in M, the first run | 233 | 23 |
| L in M, waiting for the label rule | 35 | 9 |
| First run in S | 55 of 233 | |
| First run with an S and R genome pair | 32 of 233 | |

**Labels.** On the 233, three of the four label rules give the same labels. The 32 with an S and R pair depend on the most resistant call against the other three, and the other 201 on none of the four.

**Questions.** The design left eight questions to the lead, each with the answer taken and the other one. The lead ruled at 20:36Z on 22 Nov 2026, before any C9 score was read: all eight as taken, with the words of the column Ruled. For questions 3 to 5 the scorer runs every answer from the same results, and only the answer taken is read.

| Question | Taken | Other | Ruled |
|---|---|---|---|
| 1. Inner folds | 3, at the release's tuning seed | A partition at a seed of its own | As taken |
| 2. Inner measure | Balanced accuracy of the summed inner calls | Mean of the inner folds' figures | As taken |
| 3. Ties | The later learner, ties listed | The earlier, or the default's learner when tied | The later is read; the earlier and the default's are scored beside it and decide nothing |
| 4. No pick | The default's learner on that fold, counted | The draw left unscored | As taken |
| 5. Default | Each draw's own random-split pick | The seed-42 pick at every draw | Each draw's own pick is read; the seed-42 pick is written beside it and decides nothing |
| 6. Study draws | One run at the release's seeds, no twin | A seed under this flag only | As taken, named here as a deviation from Rules, Draws |
| 7. Holm | Two pairs, condition 1 on the served figure | Per-learner rows | As taken |
| 8. First run | L less M | All of L | As taken |

**Pilot.** One draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, 93 genomes, ran in 108 s on a 16-CPU box at 8 threads. The fit with the inner pick took 30.85 s at a peak resident size of 925,212 KB, and the run used 5.09 GB of disk. Its seed-42 reports equal the release's in all 3 datasets.

**Cost.** The inner pick took 0.19 to 0.44 of each dataset's release learning seconds. At 0.44, the first run costs 18.12 xlarge-hours central. Of those, 3.47 go to fetching genomes and building matrices again on each box, which the plan's 17 for all 268 leaves out. At the plan's 1.5 it costs 24.03.

**Written.** This entry was written at 19:54Z on 22 Nov 2026, from the design fixed at 19:19Z, the lists and the pilot, and the rulings of 20:36Z were added after a review of the block. No score of any C9 run had been read.

## C3 design

**Block.** Plan block C3 asks whether other inputs separate the classes better than the presence of 31-mers (K31). The design was fixed at 23:22Z on 22 Nov 2026, before any code ran on real data, and `docs/input-arms.md` gives the method. No comparison result has been read.

**Arms.** Five arms stand against K31, the release's rule refitted in the same run. Two are built behind `grm learn --exp-inputs`: the AMRFinderPlus program's calls alone, and those calls joined before the K31 rows, less a K31 row a program column repeats. K31 plus k-mer counts, unitigs and K31 plus amino-acid k-mers are designed; the lead ruled their tools at 00:47Z on 23 Nov 2026, and they wait to be built and reviewed. A program column is one element symbol of the program's Type AMR rows, kept by the cohort range K31 uses.

**Folds.** Each arm hands the release's code a new matrix. The support floor, the collapse of identical patterns and the chi-square screen to 4,000 then run on each fit's own training genomes, over both blocks of a joined arm at once. What an arm builds per dataset reads no label, and its cohort range is implied by every training subset, so no outer test genome can change a fit's candidates.

**Measure.** Per learner, the lineage balanced accuracy is primary, with recall, specificity, both error rates with Wilson intervals, and log loss beside it. A dataset's value is the mean of seeds 42, 43 and 44; the twin is K31 at seeds 45, 46 and 47. Each arm meets K31 on its own datasets and the same folds, which the scorer checks.

**Count.** Better, worse or the same by 0.02 at 12 decimals, the arm's value less the default's, the median difference with its species interval, and the one-sided sign test. Holm runs over 30 lineage pairs, 6 arm rows by 5 learners, the twin's included. The family is fixed, so the count is taken once, after every arm has run, and `score.py --count` refuses before that. By the lead's addition to question 12, the family stays 30: an arm withdrawn before it runs keeps its five pairs at p = 1, `score.py --count --withdrawn ARM` takes it, and this entry names the arm and the time. No arm has been withdrawn. By the lead's ruling 15 at 00:47Z on 23 Nov 2026, a learner carries the change when it meets conditions 1 and 2 both, and at least 3 of 5 must carry it; 4 holds for all five; 5 holds for each learner that carries. The ruling first read condition 3 on the twin the same way, the twin carrying in fewer than 3 learners. A twin that meets 1 in three learners and 2 in none would pass that reading and fail the plan's. At 02:08Z on 23 Nov 2026, before any C3 score was read, condition 3 went back to the plan's words: it holds when the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. `conditions.tsv` writes the stricter reading beside it, 2 and 5 for all five, and that one decides nothing.

**No candidate.** Only a `--fold-file` run refuses for a fold whose training genomes give no candidate feature, and it refuses the dataset's whole study run; a drawn split's fold with no candidate gets the release's fits on no feature. By the lead's ruling 16 at 00:47Z on 23 Nov 2026, an arm refused this way where the default is not is scored as the constant model on the default's folds: each fold's training majority, R on a tie, with a score equal to its R share, clipped. Those datasets enter the count, `count.tsv` gives their number, and the count without them is written beside it and decides nothing. The run script records such a refusal and goes on; before this, it stopped the species.

**Study.** The first-run datasets in S, 55 for the arms on L and 54 for the arms on A, get one study run per arm on their reviewed `study-whole` plans, at the release's seeds, with no twin. That departs from the plan's Rules, Draws, as C9's ruling 6 did (question 11).

| Set | Datasets | Species | Less M, the first run | First run in S |
|---|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 |
| A | 237 | 21 | 204 | 54 |

**Labels.** On the datasets not in M, three of the four label rules give the same labels. Of L's 233, 32 hold an S and R pair and depend on the most resistant call against the other three; of A's 204, 21 do. The rest depend on none of the four.

**Questions.** The design leaves fourteen questions to the lead, each with the answer taken and the other one. Each is one argument of the code or the run script. The lead ruled at 00:47Z on 23 Nov 2026, with no C3 score read: question 1 changed, 2 to 14 as taken, and an addition to 12.

| Question | Taken | Other |
|---|---|---|
| 1. Unitig graph | changed by the lead: the graph over every 31-mer of the species' genomes, single-genome ones included; GRM's own compaction if its unitigs equal bcalm2's or Bifrost's as a set on the pilot species, else that program | GRM's own compaction of the species matrix's 31-mers, the design's answer, not taken |
| 2. Unitig presence | every 31-mer of the unitig in the genome; the inputs file also counts the genome and unitig pairs held in part | any 31-mer; a share of them |
| 3. Unitig scope | the whole cohort, with no label | each outer training fold |
| 4. Gene caller | Prodigal 2.6.3, single mode, table 11 | Pyrodigal |
| 5. Amino-acid k | 10, spanning 30 nucleotides | 8 or 12 |
| 6. Counts | GRM's own counter, one bin, count of 2 or more, checked on SYNTHETIC sequence against a recount written apart | KMC or Jellyfish; bins 2, 3 and 5 |
| 7. Program column | element symbol | subclass; class |
| 8. Program rows | Type AMR, both scopes, every method | every Type; core scope only |
| 9. Join | one candidate pass and one screen over both blocks | a screen per block |
| 10. Collapse order | the arm's columns before K31's | K31's first |
| 11. Study draws | one run at the release's seeds, no twin | a seed under this flag only |
| 12. Holm | 30 pairs, the twin's included; an arm withdrawn before it runs keeps its 5 pairs at p = 1 | a family per arm |
| 13. First run | L and A less M | all of L and A |
| 14. Default and twin | C3's own, in the same run as its arms | shared with another block's runs |

**One run.** By the lead's ruling 17, no C3 fit starts on the two built arms alone. One run with one binary covers the default, the twin and all five arms after the three waiting arms are built and reviewed, since the default is refitted in the same run with the same binary. This replaces the design's later boxes of their own for unitigs and amino-acid k-mers.

**Pilot.** One draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, 93 genomes, ran end to end on a 16-CPU box from 00:03:48Z to 00:12:32Z on 23 Nov 2026. The program took 372.69 s for the 93 tables, each equal to its recorded sha256. The fits took 17.98 s for K31, 2.74 s for the program alone and 11.79 s for K31 plus the program, each under 940,000 KB at peak. K31's seed-42 reports equal the release's in all 3 datasets.

**Fault.** A first pilot stopped at K31 plus the program: a K31 row repeated a program column over the cohort, and the final models need each cohort pattern once. The joined arm now leaves such a K31 row out, as its training collapse already did, and counts it: 9, 9 and 7 rows in the three datasets. A synthetic test and planted fault 2d hold it.

**Cost.** Each arm's fit costs its pilot ratio to the release's learning seconds: 1.08 for K31, 0.17 for the program alone and 0.71 for K31 plus the program. The program takes 0.0178 core-hours a genome, 2.7 times the release's 0.0065, on 21,625 genomes. The built arms' part of the one run, with K31 and the twin, comes to 40.42 xlarge-hours, 24.07 of them the program; with data steps and loads, 47.12 box-hours over 23 boxes. The three waiting arms add their own, measured on their pilots.

**Written.** This entry was written at 00:16Z on 23 Nov 2026, from the design fixed at 23:22Z on 22 Nov, the lists and the second pilot. The lead's rulings of 00:47Z were carried into it from 01:20Z to 01:24Z, in a review with a recount before filing (ruling 18). No comparison result had been read.

## Mechanism design

**Blocks.** Plan blocks M1 and M2 read the AMRFinderPlus program's tables beside the release's reports and fit nothing. M1 asks whether each learner's recall differs between resistant genomes whose catalogue calls are point mutations and those whose calls are acquired genes. M2 asks whether the mix of catalogue calls among deposited genomes moves between collection years. The design was fixed at 03:22Z on 23 Nov 2026, before any script of these blocks read a real table, score, stratum or share, and `docs/mechanism-split.md` gives the method. No result has been read.

**Rule.** Both blocks reuse the AMRFinderPlus comparison's rule as it stands: each report's published rule, written for the 237 datasets of A as `scripts/mechanism/drug_rules.tsv`, and the adapter's counted rows (`tool_default`). Each run holds its own calls to the kit's adapter on every genome and its counts to the comparison's record.

**M1.** On the 204 datasets of A not in M, the lineage split is primary. A hard check recomputes each learner row's counts from its scores before any stratum is counted. A genome's stratum (point only, gene only, both, neither) comes from its rule-meeting rows, never from its label. Per learner, the median of recall on gene only less recall on point only, over datasets with 10 or more resistant genomes in both, with its species interval, the sign test and Holm over the five learners.

**M2.** Periods 2010 to 2014 against 2015 to 2019, then 2000 to 2009 against 2010 to 2019. A cell is a species with a drug class. One vote per BioProject with 5 or more genomes in a period, 3 or more such BioProjects a period, and a permutation over BioProjects, with Holm over the cells that pass. Six guards stand beside it and decide nothing.

**Not built.** M3 is not built: a feature card names a gene in 1,432 of the release's 21,443 final-model cards, in none of the cards of 6 of the 23 species with cards (the plan wrote 7; the probe's output gives 6), and by BV-BRC's annotation rather than the catalogue's name. M4, the lineage curve for Mycobacterium tuberculosis, has a script, `tb_curve.py`, with SYNTHETIC tests and no run: TBProfiler must run again on the 2,000 genomes first.

**Questions.** The plan leaves nine questions to the lead, each with the option taken and the other. A tenth came from the dry run, which stopped at M2 on a row of Subtype `POINT_DISRUPT`: the program names a disrupted core gene `<gene_family>_<change>`, which neither the catalogue's alleles nor its hierarchy hold. Taken: the gene_family before the underscore, when the catalogue lists that family with Subtype `POINT_DISRUPT`; any other symbol still stops the species. Other: leave `POINT_DISRUPT` rows out of the family mix. The plan file is unchanged.

| Question | Taken | Other | Ruled |
|---|---|---|---|
| 1. `POINT_DISRUPT` rows | gene, as the brief words it, with those genomes counted beside each stratum | point | changed: point, with those genomes counted beside each stratum and, in point only and both, the genomes whose point rows are `POINT_DISRUPT` alone; the same in M2 and the time-split link |
| 2. "Under 5" | covers zero | show 0 | changed: 0 is written 0, 1 to 4 "under 5"; withdrawn at 05:30Z by question 12 |
| 3. Species interval | none under 6 species | an interval at any count | as taken |
| 4. BioProject vote | 5 or more genomes in the period; the rest enter the share per genome only | every BioProject | as taken |
| 5. Flagged BioProject | resistant share over the datasets whose rule names the class | every dataset | as taken |
| 6. Family | the allele's gene_family, the symbol where it is one, else its first family in the hierarchy | the symbol | as taken |
| 7. Lineage clusters | joined over the species' reports, given only where they agree | per report | as taken |
| 8. Time-split determinant | a point row as one token "point", its symbol beside | each symbol | as taken |
| 9. M4 names | no variant name in any output until ruled | the catalogue's variant names | as taken |
| 10. `POINT_DISRUPT` family | the gene_family before the underscore, if the catalogue lists it with that subtype | out of the family mix | closed by question 1: a `POINT_DISRUPT` row enters no family mix, which is over Subtype `AMR` rows alone |
| 11. Log loss in M1 | recall alone, as the brief fixes M1 | the log loss of each stratum's resistant genomes beside recall | as taken |
| 12. A count of 1 to 4 written "under 5" | as ruled on question 2 | the count as it is | changed: every count as it is, 0 to 4 included, and "under 5" dropped; a recall over 10 or more resistant genomes, an M2 share only in a cell whose BioProject vote passes, M4's shares over 5 or more; a family or gene named only where 5 or more genomes carry it, the rest in one row |
| 13. The flagged guard's years | a BioProject over its genomes in every year | over the cell's periods | as taken: the flag describes the study, whose design does not change by period |

An eleventh came while the text was written. The project's wording rules put recall beside log loss, and M1 as fixed gives recall alone. The review raised two more. A count written "under 5" could be read back from the cells beside it, such as the all-genomes row less the other strata. Its recount read the flag over the cell's periods, where the script reads every year.

**Rulings.** The lead ruled on the eleven questions at 04:13Z on 23 Nov 2026, and on questions 12 and 13 at 05:30Z, before any result was read. The plan file stays unchanged; the changes are dated beside the review's report and carried into the scripts, their tests and `docs/mechanism-split.md`. The lineage guard's count of differing pairs counts each pair of genomes once, as the plan words it; the first script counted a pair once per pair of reports it differed in.

**Checks.** 121 SYNTHETIC tests in `tests/mechanism` (18 on M4). The ten planted faults and the review's 28 each fail a test; 14 of the review's passed every test before its patches. Two of the 28 no longer apply under the rulings of 05:30Z, and the review's nine faults on the rulings' rules each fail a test. The dry run on Klebsiella aerogenes (93 genomes, 3 datasets) ran end to end in 2.5 s with 30.6 MB at peak; the hard check passed on 3 of 3 datasets, 36 of 36 rows, and its scan found no genome id. Run again after the patches, the hard check passed on the same rows and the scan, which now also reads each table's counts, found no hit.

**Cost.** Where the tables are given, about 565 s of compute for the 21 species, from the dry run's steps and a timing on 2,000 SYNTHETIC genomes. Where the program must run again, about 265 CPU hours on 21,625 genomes, from the refit's own run of the program on the same 93.

**Written.** This entry was written at 03:48Z on 23 Nov 2026, after the dry run, and M4's script and question 11 were added at 04:06Z. The rulings of 04:13Z were carried in with the review's patches from 04:27Z to 04:49Z, and these checks at 05:11Z. The rulings of 05:30Z were carried in from 05:38Z to 05:47Z. No result of any block had been read.

## C3 arms

**Block.** C3's three waiting arms are built behind `grm learn --exp-inputs` under the lead's rulings of 00:47Z on 23 Nov 2026: K31 plus k-mer counts (`k31-counts`), K31 plus amino-acid k-mers (`k31-aa`) and unitigs (`unitigs`). `docs/input-arms.md` gives the method. No comparison result has been read: no pilot's `score/` folder was opened.

**Builds.** `grm exp-input-build counts|aa|unitigs` runs once a species after `grm matrix`. It reads no label, fold or seed, and writes a pattern matrix at support 2 in the release's format, with the species matrix's genomes in order, and a counts-only record. Counts: each genome's canonical 31-mers seen 2 or more times in its assembly, by the release's counter with its floor at 2. Amino-acid k-mers: Prodigal 2.6.3 in single mode at table 11, each protein held to GRM's own table 11 translation of its contig at Prodigal's coordinates, then its 10-mers over the 20 standard residues. Unitigs: the graph over every 31-mer of the species' genomes; a genome holds a unitig only when it holds every 31-mer of it, and the pairs held in part are counted. Each fit projects the build on its cohort with every label set to S. The count and amino-acid arms put their rows first, then K31 less each row whose cohort pattern equals an added row's, counted; the unitig arm takes K31's place.

**Graph.** By question 1 as the lead changed it, GRM may build the graph only if its unitigs equal an outside build's as a set. The outside build is bcalm2 2.2.3 (`-kmer-size 31 -abundance-min 1`, release binaries pinned by sha256). On *Klebsiella aerogenes*, 93 genomes, each build holds 31,083,227 31-mers and 1,370,388 unitigs, 1,370,388 in both and 0 in one only. On *Escherichia coli*, 2,000 genomes, each holds 104,518,454 31-mers and 5,493,343 unitigs, 5,493,343 in both and 0 in one only. So GRM builds it: `GRAPH` defaults to grm. `GRAPH` check runs both and stops the species on any unitig or 31-mer in one build only.

**Memory.** The first *Escherichia coli* graph build was killed by the kernel at 31.8 GB: the union of the species' 31-mers kept every bucket's read buffer. It now keeps each bucket's distinct 31-mers only. The rerun peaked at 4,268,980 KB. A step the kernel kills writes `FAILED` naming it; 32 GB is the largest box, so the restart takes a smaller `JOBS` and redoes that step alone.

**Screen count.** With a joined arm, `grm learn` writes beside each pick file how many of each block's candidates pass the 4,000 screen, per outer fold, on the forest's, boosted trees' and lasso's own fitting genomes. On the pilot each joined arm made 90 such fits (3 datasets, 2 splits, 5 folds, 3 learners), 360,000 screened rows in all. Amino-acid rows took 63,179 of them in K31 plus amino-acid k-mers, and count rows 5,172 in K31 plus counts. In K31 plus the program, program columns took 2, and in 89 of the 90 fits none passed.

**Pilot.** One draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, 93 genomes, seed 42, the default and all five arms, ran alone on a 16-CPU box from 04:13:58Z to 04:19:16Z on 23 Nov 2026, with `GRAPH` check and the program's recorded tables. The builds took 4.52 s for the counts, 36.43 s for Prodigal and the amino-acid build, and 70.09 s for bcalm2 and the unitig build. The fits took 32.51 s for K31, 3.82 s for the program alone, 28.72 s for K31 plus the program, 14.00 s for K31 plus counts, 11.37 s for K31 plus amino-acid k-mers and 16.75 s for unitigs. The peak was 2,404,604 KB, in bcalm2. K31's seed-42 reports equal the release's in all 3 datasets and 513 of 513 checks pass. This run used 1.10 to 2.05 times the CPU seconds of an earlier run of the same steps on the same box, at an earlier commit (each step over 5 s), and the cost takes this slower run.

**Second species.** *Escherichia coli*, 2,000 genomes, 10.232 Gbp, ran the three label-free builds alone, with no fit and no label read (question 18). Prodigal took 712.51 s at 12 jobs and 226,928 KB. The count build took 75.27 s at 1,645,192 KB, 189,522 patterns. The amino-acid build took 46.10 s at 1,730,124 KB, 9,795,431 proteins and 1,008,027 patterns, with 13.5 GB of disk at its peak. The unitig build, with the comparison, took 232.75 s at 4,268,980 KB, 3,351,087 patterns and 150,437 genome and unitig pairs held in part. bcalm2 alone took 689.13 s at 9,445,936 KB.

**Klebsiella pneumoniae.** 2,000 genomes, 11.126 Gbp, the largest species. Its refit peaked at 5,045,392 KB in the matrix step and 3,781,356 KB in `grm learn`. Scaled from *Escherichia coli* a Gbp, the unitig build peaks at about 4,641,973 KB and the count and amino-acid builds under 1,900,000 KB. A unitig fit holds K31's matrix and the unitig matrix, about 1 GB on disk. The expected peak is 6,173,174 KB, in that fit, of the box's 33,554,432 KB.

**Unflagged.** With no flag, at seed 42 on the pilot's 3 datasets, the pick and result files equal main's but for their seconds fields, 3 of 3 each. With `--exp-inputs program` the 6 files equal main's but for seconds, and with `k31-program` the same 6, beside 6 screen files main does not write. `check_fits.sh` holds 17 of 17 files byte for byte. `EVALUATION_REVISION` stays 14.

**Faults.** `scripts/c3-inputs/planted_faults.py` holds 31 faults, main's 16 and 15 for the three arms, and each is caught: its test passes on the clean tree and fails with the fault in, on the branch's code. In an earlier run, fault 11b, a matrix arm's floor taken over the whole cohort, was not caught, since its test's added rows were all of one density and none sat near the floor. The test now has sparse rows. Fault 17's text no longer occurred after the memory fix, and it was planted again on the new code.

**Cost.** Each arm's fit costs its ratio on the pilot to the release's learning seconds of the same 3 datasets: 1.96 for K31, 0.23 for the program alone, 1.73 for K31 plus the program, 0.84 for K31 plus counts, 0.69 for K31 plus amino-acid k-mers and 1.01 for unitigs. The builds scale with the species' Gbp at the larger of the pilot's and *Escherichia coli*'s rates. With the program run on every species of A, the one run comes to 82.50 xlarge-hours central and 200.53 upper, with lifetimes of 209.25 h over 23 boxes of 16 CPUs and 32 GB: $16.50 central, $40.11 upper and $41.85 at the lifetimes, at $0.20 an hour. With the program's tables already on the 16 refit boxes, it comes to 61.12 central and 168.39 upper, with lifetimes of 176.75 h: $12.22, $33.68 and $35.35. Per arm, central and upper: K31 12.12 and 36.37, the twin 10.01 and 30.02, K31 plus counts 5.92 and 17.40, K31 plus amino-acid k-mers 6.74 and 17.25, unitigs 7.42 and 21.09, the program alone 2.01 and 6.02, K31 plus the program 10.09 and 30.28, and the program on every species of A 24.30 and 36.45. The longest lifetime is *Klebsiella pneumoniae*'s, 52.75 h with the program and 49.50 h with its tables.

**Questions.** The build leaves six questions to the lead, each with the answer taken and the other one.

| Question | Taken | Other |
|---|---|---|
| 18. Second species for the build timings | *Escherichia coli*, 2,000 genomes, 10.232 Gbp, already on the run box, the nearest size to *Klebsiella pneumoniae* (2,000, 11.126 Gbp) | *Streptococcus pneumoniae* after a 29 GB fetch; *Campylobacter jejuni* (508) |
| 19. The unitig matrix's own support | 2, as the count and amino-acid matrices; each fit's floor then applies | 1 |
| 20. A Prodigal protein unequal to GRM's table 11 translation | stops the build | left out and counted |
| 21. The screen count's fits | each learner's outer fit | the inner settings fits and final models too |
| 22. The graph in the run | GRM's, `GRAPH` grm | `GRAPH` check on every species: bcalm2 runs too, 689.13 s and 9,445,936 KB on *Escherichia coli* |
| 23. The program's tables | the 16 species' refit boxes, resumed, each table held to its recorded sha256 | a new box per species that runs the program |

**Written.** This entry was written at 04:58Z on 23 Nov 2026, from the builds on both species, the pilot of all five arms and the fault run at the branch's code. No comparison result had been read, and no pilot's `score/` folder was opened.
