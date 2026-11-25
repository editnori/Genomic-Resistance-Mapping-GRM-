# Assumption plan: corrections

The plan in `docs/assumption-plan.md` is fixed, and its text does not change. This page lists what was found wrong in it after it was fixed, what a later design added to it, what a ruling settled where it left a rule open, and what holds instead. Each entry gives its date. The first two were written on 22 Nov 2026, before any fit.

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

**Builds.** `grm exp-input-build counts|aa|unitigs` runs once a species after `grm matrix`. It reads no label, fold or seed, and writes a pattern matrix at support 2 in the release's format, with the species matrix's genomes in order, and a counts-only record. Counts: each genome's canonical 31-mers seen 2 or more times in its assembly, by the release's counter with its floor at 2. Amino-acid k-mers: Prodigal 2.6.3 in single mode at table 11, each protein found by the sequence number in its `ID=` (Prodigal numbers records from 1 in file order) and held to GRM's own table 11 translation of that contig, read as Prodigal reads a sequence line (the bytes `A` to `z` only), at Prodigal's coordinates, then its 10-mers over the 20 standard residues. Unitigs: the graph over every 31-mer of the species' genomes; a genome holds a unitig only when it holds every 31-mer of it, and the pairs held in part are counted. Each fit projects the build on its cohort with every label set to S. The count and amino-acid arms put their rows first, then K31 less each row whose cohort pattern equals an added row's, counted; the unitig arm takes K31's place.

**Graph.** By question 1 as the lead changed it, GRM may build the graph only if its unitigs equal an outside build's as a set. The outside build is bcalm2 2.2.3 (`-kmer-size 31 -abundance-min 1`, release binaries pinned by sha256). On *Klebsiella aerogenes*, 93 genomes, each build holds 31,083,227 31-mers and 1,370,388 unitigs, 1,370,388 in both and 0 in one only. On *Escherichia coli*, 2,000 genomes, each holds 104,518,454 31-mers and 5,493,343 unitigs, 5,493,343 in both and 0 in one only. So GRM builds it: `GRAPH` defaults to grm. `GRAPH` check runs both and stops the species on any unitig or 31-mer in one build only.

**Memory.** The first *Escherichia coli* graph build was killed by the kernel at 31.8 GB: the union of the species' 31-mers kept every bucket's read buffer. It now keeps each bucket's distinct 31-mers only. The rerun peaked at 4,268,980 KB with the comparison, and GRM's build alone at 2,930,764 KB in 214.22 s. A step the kernel kills writes `FAILED` naming it, and `logs/run.log` says it was killed for memory: GNU time records such a step as the shell's exit status 137, and a fit run under `time -v` as signal 9, and the run script reads both. 32 GB is the largest box, so the restart takes a smaller `JOBS` and redoes that step alone.

**Screen count.** With a joined arm, `grm learn` writes beside each pick file how many of each block's candidates pass the 4,000 screen, per outer fold, on the forest's, boosted trees' and lasso's own fitting genomes. On the pilot each joined arm made 90 such fits (3 datasets, 2 splits, 5 folds, 3 learners), 360,000 screened rows in all. Amino-acid rows took 63,179 of them in K31 plus amino-acid k-mers, and count rows 5,172 in K31 plus counts. In K31 plus the program, program columns took 2, and in 89 of the 90 fits none passed.

**Pilot.** One draw of *Klebsiella aerogenes* with ertapenem, imipenem and meropenem, 93 genomes, seed 42, the default and all five arms, ran alone on a 16-CPU box from 04:13:58Z to 04:19:16Z on 23 Nov 2026, with `GRAPH` check and the program's recorded tables. The builds took 4.52 s for the counts, 36.43 s for Prodigal and the amino-acid build, and 70.09 s for bcalm2 and the unitig build. The fits took 32.51 s for K31, 3.82 s for the program alone, 28.72 s for K31 plus the program, 14.00 s for K31 plus counts, 11.37 s for K31 plus amino-acid k-mers and 16.75 s for unitigs. The peak was 2,404,604 KB, in bcalm2. K31's seed-42 reports equal the release's in all 3 datasets and 513 of 513 checks pass. This run used 1.10 to 2.05 times the CPU seconds of an earlier run of the same steps on the same box, at an earlier commit (each step over 5 s), and the cost takes this slower run.

**Second species.** *Escherichia coli*, 2,000 genomes, 10.232 Gbp, ran the three label-free builds alone, with no fit and no label read (question 18). Prodigal took 712.51 s at 12 jobs and 226,928 KB. The count build took 75.27 s at 1,645,192 KB, 189,522 patterns. The amino-acid build took 46.10 s at 1,730,124 KB, 9,795,431 proteins and 1,008,027 patterns, with 13.5 GB of disk at its peak. The unitig build, with the comparison, took 232.75 s at 4,268,980 KB, 3,351,087 patterns and 150,437 genome and unitig pairs held in part. bcalm2 alone took 689.13 s at 9,445,936 KB.

**Klebsiella pneumoniae.** 2,000 genomes, 11.126 Gbp, the largest species. Its refit peaked at 5,045,392 KB in the matrix step and 3,781,356 KB in `grm learn`. Scaled from *Escherichia coli* a Gbp, the unitig build peaks at about 4,641,973 KB and the count and amino-acid builds under 1,900,000 KB. A unitig fit holds K31's matrix and the unitig matrix, about 1 GB on disk. The expected peak is 6,173,174 KB, in that fit, of the box's 33,554,432 KB.

**Unflagged.** With no flag, at seed 42 on the pilot's 3 datasets, the pick and result files equal main's but for their seconds fields, 3 of 3 each. With `--exp-inputs program` the 6 files equal main's but for seconds, and with `k31-program` the same 6, beside 6 screen files main does not write. `check_fits.sh` holds 17 of 17 files byte for byte. `EVALUATION_REVISION` stays 14.

**Faults.** `scripts/c3-inputs/planted_faults.py` holds 31 faults, main's 16 and 15 for the three arms, and each is caught: its test passes on the clean tree and fails with the fault in, on the branch's code. In an earlier run, fault 11b, a matrix arm's floor taken over the whole cohort, was not caught, since its test's added rows were all of one density and none sat near the floor. The test now has sparse rows. Fault 17's text no longer occurred after the memory fix, and it was planted again on the new code. A review planted 15 faults of its own, none of these 31, and the tests caught 8. Each of the other 7 now has a test that catches it: a reverse-strand gene's start read at its left end, a unitig held with one 31-mer missing, the unitig pair counts summed without the first genome, a build's genomes held to the matrix's by number only, a learner that screened fewer than 4,000 of more candidates, a build's sha256 list of `build.json` alone, and the end's genome-id count without `logs/`.

**Cost.** Each arm's fit costs its ratio on the pilot to the release's learning seconds of the same 3 datasets: 1.96 for K31, 0.23 for the program alone, 1.73 for K31 plus the program, 0.84 for K31 plus counts, 0.69 for K31 plus amino-acid k-mers and 1.01 for unitigs. The builds scale with the species' Gbp at the larger of the pilot's and *Escherichia coli*'s rates. With the program run on every species of A, the one run comes to 82.50 xlarge-hours central and 200.53 upper, with lifetimes of 209.25 h over 23 boxes of 16 CPUs and 32 GB: $16.50 central, $40.11 upper and $41.85 at the lifetimes, at $0.20 an hour. With the program's tables already on the 16 refit boxes, it comes to 61.12 central and 168.39 upper, with lifetimes of 176.75 h: $12.22, $33.68 and $35.35. Per arm, central and upper: K31 12.12 and 36.37, the twin 10.01 and 30.02, K31 plus counts 5.92 and 17.40, K31 plus amino-acid k-mers 6.74 and 17.25, unitigs 7.42 and 21.09, the program alone 2.01 and 6.02, K31 plus the program 10.09 and 30.28, and the program on every species of A 24.30 and 36.45. The longest lifetime is *Klebsiella pneumoniae*'s, 52.75 h with the program and 49.50 h with its tables.

**Questions.** The build left six questions to the lead. The lead ruled all six at 05:05Z on 23 Nov 2026, each as taken below.

| Question | Ruled | Other |
|---|---|---|
| 18. Second species for the build timings | *Escherichia coli*, 2,000 genomes, 10.232 Gbp, already on the run box, the nearest size to *Klebsiella pneumoniae* (2,000, 11.126 Gbp) | *Streptococcus pneumoniae* after a 29 GB fetch; *Campylobacter jejuni* (508) |
| 19. The unitig matrix's own support | 2, as the count and amino-acid matrices; each fit's floor then applies | 1 |
| 20. A Prodigal protein unequal to GRM's table 11 translation | stops the build | left out and counted |
| 21. The screen count's fits | each learner's outer fit | the inner settings fits and final models too |
| 22. The graph in the run | GRM's, `GRAPH` grm | `GRAPH` check on every species: bcalm2 runs too, 689.13 s and 9,445,936 KB on *Escherichia coli* |
| 23. The program's tables | the 16 species' refit boxes, resumed, each table held to its recorded sha256 | a new box per species that runs the program |

**Written.** This entry was written at 04:58Z on 23 Nov 2026, from the builds on both species, the pilot of all five arms and the fault run at the branch's code. No comparison result had been read, and no pilot's `score/` folder was opened. A review recounted the three builds on SYNTHETIC genomes, carried in the rulings of 05:05Z, and fixed the protein read, the memory message and the scorer's log; it was carried in at 06:05Z the same day, and no comparison result had been read and no `score/` folder opened.

## Mechanism split, first run

**Run.** The first run of blocks M1 and M2 on the 21 species of A began at 06:39Z on 23 Nov 2026 from the merged scripts (`253d60e`). 17 species ended DONE. Three stopped on two defects of the scripts, each fixed with its tests on a branch of its own. No file under any species' `results/` had been opened when the defects were found, ruled on or fixed.

**Defect A, the database fetched before the program is installed.** Neisseria gonorrhoeae (2,000 genomes, no program tables on its machine, DB given as `-`) failed at step `db` at 06:40:52Z: `fetch_database.sh: amrfinder_index is not on PATH`. With DB `-`, `run_species.sh` called the kit's `fetch_database.sh`, which needs the program's `amrfinder_index`, while the program's environment was installed only inside the later `tables` step, and only when TABLES was `-`. So DB `-` could never work. The install (micromamba pinned by sha256, the kit's `spec/env-amrfinder.explicit.txt`, `amrfinder --version` 4.2.7) is now a step of its own between `amrkit` and `db`, run once whenever DB or TABLES is `-` and skipped on a restart. With TABLES and DB both given, nothing is installed, as before. The run got past it without the fix: the lead gave the five species without a database a copy of database 2026-08-07.1 taken from a refit machine, held to the kit's `spec/db-2026-08-07.1.sha256` (43 of 43 files), and Vibrio cholerae then ended DONE with 91 of 91 program tables equal to their recorded sha256.

**Defect B, a gene row named above the family level.** Streptococcus pneumoniae failed at step `history` at 06:51:33Z and Streptococcus agalactiae at 06:54:19Z: `history.py: refused: an element symbol of Subtype AMR: neither the catalogue nor the hierarchy gives it a family`. The AMRFinderPlus program names some hits at a node of the database's `ReferenceGeneHierarchy.txt` above the gene families. The lead counted them on the two machines, by database names and counts only. In Streptococcus pneumoniae, `msr` is in 7 of 2,000 genomes, 7 rows, Method `INTERNAL_STOP`, Class `MACROLIDE`; its chain is `msr > abc-f > AMR > ALL`, and 6 of its 8 child nodes are gene families. In Streptococcus agalactiae, `lsa` is in 1 of 141 genomes, 1 row, Method `BLASTX`, Class `LINCOSAMIDE/STREPTOGRAMIN`; its chain is `lsa > abc-f > AMR > ALL`, and 6 of its 6 child nodes are gene families. The other Subtype `AMR` rows resolve: 4,644 in Streptococcus pneumoniae and 242 in Streptococcus agalactiae.

**Ruling on question 14.** The lead ruled at 06:56Z on 23 Nov 2026, before any M2 result was read. A Subtype `AMR` row whose symbol is a node of the database's `ReferenceGeneHierarchy.txt` with no `gene_family` at or above it is counted in the family mix under that node's own name, as the database names it. A symbol in neither file still stops the run. The floor on names (5 or more genomes in the row's cell, the rest under "other families") holds for such a name as for any family. The step's log line states how many rows and genomes of the species were named this way. No output file changes its columns. `Database.family` takes the symbol itself as a fourth step, after the allele's family, the symbol as a family and its first family in the hierarchy, so a species that ended DONE keeps every byte of its results.

**Checks.** `tests/mechanism` holds 134 SYNTHETIC tests, 13 more than at `253d60e`. Fix A's five tests run the script's own text from the install step to the database step with stand-in tools; on the script as merged, 4 of the 5 fail. Fix B's eight tests run the three family tests again on a database with a node above the families and test the node's own name, the refusal of a symbol in neither file, the log line's count and the floor on names. Three faults were planted, each caught by its test: the install moved after the database step, a node given its first child's family, and a symbol in neither file let through. On the scripts as merged, 4 of fix B's 8 tests fail; the 4 that pass hold symbols that resolved before.

**Proof.** With both fixes, `run_species.sh` ran Streptococcus pneumoniae twice from the refit's 2,000 program tables, once with the refit's copy of database 2026-08-07.1 and once with DB `-`. The first installed nothing. The second installed AMRFinderPlus 4.2.7, then fetched the database from NCBI, 43 of 43 files equal to the kit's list. Both ended DONE, both scans passed, both history log lines counted 7 rows of 7 genomes named at a node above the gene families, as the lead counted `msr`, and the two runs' `results/SHA256SUMS` have the same sha256. No file under `results/` was opened.

**Written.** This entry was written at 07:09Z on 23 Nov 2026 and the proof added at 07:10Z. No result of any block had been read.

## C4 design

**Block.** Plan block C4 asks whether the chi-square screen or the candidate floor removes the inputs that separate the classes, or keeps lineage markers. The design was fixed at 04:03Z on 23 Nov 2026, before any code of the block existed and before anything ran on real data, and `docs/screen-arms.md` gives the method. No comparison result has been read.

**Arms.** Three arms stand against the release's rule, refitted in the same run with the same binary. `grm learn --exp-screen cmh` ranks the candidates of the forest, boosted trees and lasso by the Cochran-Mantel-Haenszel statistic over the release's lineage clusters and keeps the top 4,000 (`cmh-4000`) or, with `--screened-features 16000`, the top 16,000 (`cmh-16000`). `grm learn --exp-floor 3` fixes the support floor at 3 in place of max(3, n/200) (`floor-3`). The twin is the default at seeds 45, 46 and 47.

**Statistic.** D² / V over the informative strata of the fit's own fitting genomes, those of 2 or more genomes and both classes, in ascending cluster order; no continuity correction, two-sided, ties to the lower row, the budget filled with zeros by lower row. A stratum of one genome or one class adds nothing, so a pure lineage marker scores 0. With one stratum it ranks as the chi-square.

**Folds.** The strata read no label and each screen takes them on the genomes it fits, so no outer test genome can change a screen. The candidates are the default's. Folds come from the labels, the lineages and the seed, so each arm meets the default on the same folds, which the scorer checks. A screen arm's SCM and CART read no screen and must make the default's calls in every fold. Below 800 genomes the floor arm is the default and its results must equal the default's, timing apart.

**Count.** As C3's, the arm's value less the default's on the lineage primary measure, with Holm over 21 pairs (the lead's rulings 10 and 12): `cmh-4000` and `cmh-16000` each on the forest, boosted trees and lasso over every dataset, `floor-3` on all five learners over the datasets of 800 genomes or more, the ones its floor can move, and the twin on all five over every dataset and again over those of 800 genomes or more. A screen arm's conditions are read on the three learners it touches (1, 2 and 5 for each, 4 for the three). The floor arm's follow the lead's C3 ruling 15, over the datasets of 800 genomes or more. Condition 3 follows the plan's words as corrected for C3 at 02:08Z on 23 Nov 2026, with the twin read on each arm's datasets: the twin meets 1 in fewer than 3 of 5 learners and 2 in fewer than 3 of 5, each read apart. The floor arm's count over the whole first run, the 141 datasets it cannot move entering at 0, is written beside it and decides nothing.

**Study.** The first run's 55 datasets in S get one study run per arm on their reviewed `study-whole` plans, at the release's seeds, with no twin, as C9's ruling 6 and C3's question 11 did.

| Set | Datasets | Species | Less M, the first run | First run in S | Floor arm moves |
|---|---:|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 | 92 in 11 species, 24 in S |

**Questions.** The design leaves seventeen questions to the lead, each with the answer taken and the other one. Each is one argument of the code or the run script. The lead ruled them at 05:52Z on 23 Nov 2026, with no C4 score read: all as taken, except questions 10 and 12, changed as the table shows.

| Question | Taken | Other |
|---|---|---|
| 1. Continuity correction | none, as the default's statistic | Mantel and Haenszel's 0.5 |
| 2. Direction | two-sided | one-sided toward resistance |
| 3. Arithmetic and ties | f64 in fixed stratum order; ties to the lower row; held to exact rationals in a plain Python implementation | exact rationals in the screen |
| 4. A stratum of one genome or one class | adds 0 to D and V | pooled into one stratum |
| 5. Most strata uninformative | the budget filled with zeros by lower row, counted in the inputs file | fewer columns than the budget |
| 6. Strata | the release's lineage clusters on each fit's fitting genomes, on every split | clusters cut again per fold |
| 7. Candidates | the default's | uncollapsed rows |
| 8. Budget 16,000 | measured in the pilot | |
| 9. The floor arm's other changes | the projection's cohort range only | |
| 10. Floor arm's datasets | ruled: the 92 of 800 genomes or more, in 11 species, fixed by cohort size before any score; the whole first run, 141 at 0, beside and deciding nothing. With 141 of 233 differences fixed at 0 the median is 0 by construction, so condition 1 could never be met. The twin is read over the same 92 for the floor arm's condition 3 | taken at first: the whole first run, the 92 beside |
| 11. SCM and CART under a screen arm | refitted and checked, not counted | counted as pairs |
| 12. Holm | ruled: 21 pairs, `cmh-4000` and `cmh-16000` by 3 learners, `floor-3` by 5 over the 92, the twin by 5 over the 233 and by 5 over the 92 | taken at first: 16 pairs; or 20, with SCM and CART of the screen arms at p = 1 |
| 13. Conditions | the learners an arm touches; ruling 15 for the floor arm; condition 3 as corrected | the twin read on three learners |
| 14. Study draws | one run at the release's seeds, no twin | a seed under this flag only |
| 15. First run | L less M | all of L |
| 16. Default and twin | C4's own, in the same run as its arms | shared with C3's run |
| 17. No candidate (ruling 16) | kept by import; it cannot arise | |

**Negative control.** On a SYNTHETIC cohort of 400 genomes in 20 strata, with 1,000 lineage markers among 2,000 rows and labels permuted within strata 200 times, a permutation rejects when the screen's top 100 hold more markers than a one-sided hypergeometric test at 0.05 allows. The adjusted screen rejected in 7 of 200 against a cut of 15 fixed before it ran; the chi-square, deciding nothing, in 200 of 200.

**Pilot.** One draw of *Enterobacter cloacae* at seed 42, its 9 datasets (cefepime, ciprofloxacin, gentamicin, imipenem, levofloxacin, meropenem, sulfamethoxazole/trimethoprim, tobramycin and trimethoprim; all but trimethoprim in S), the default and the three arms on the drawn splits and the 8 study plans, ran on a 16-CPU box from 04:31:01Z on 23 Nov 2026. It stopped at its score step at 04:37:21Z, as the scorer was not yet committed, and its restart on the same binary ended `DONE` at 04:48:49Z: 68 of 68 results, 2,231 of 2,231 checks, 0 genome ids in what it keeps. Its scores were not opened. The fit steps took 54.27 s for the default with its reports, 31.67 s for `cmh-4000`, 35.59 s for `cmh-16000` and 26.36 s for `floor-3`; each fit process peaked under 1,303,000 KB. No pilot dataset reaches 800 genomes, so the floor arm's moved fits are not timed.

**Unflagged bytes.** `scripts/c4-screen/unflagged.sh` ran `grm` built at main (`6c3697b`) and at this branch, whose code, the binary less its build id, is the pilot's, with no C4 flag and no network, from 05:13:10Z to 05:17:13Z on 23 Nov 2026. On the pilot's data root, with 9 reports and 17 result files a side, the 360 files of the two copies were 317 equal in bytes and 43 equal but for the run clocks, none differing. On the root of 270 SYNTHETIC genomes made from `test_dna.fa`, the 1,141 files were 1,098 and 43, none differing.

**Faults.** 17 planted faults, each in a copy of the tree, were each caught by their test, which passes with no fault.

**Cost.** Each fit costs its pilot ratio of fit wall time to the release's learning seconds of the pilot's datasets (41.43 s): 1.3098 for the default and the twin, 0.7643 for `cmh-4000`, 0.8590 for `cmh-16000` and 0.6362 for `floor-3`, over 3 seeds on the first run's 233 datasets and one study run on its 55. The fits come to 27.10 xlarge-hours, 7.56 of them the default and 6.49 the twin; with C9's data steps and the matrix loads, 34.39 box-hours central over 23 boxes ($6.88) and 97.35 upper ($19.47), above the plan's 25. Sharing the default, the twin and the data roots with C3's run would save 18.93 of the central hours, but needs one binary for both blocks, and C3's run waits for its three unbuilt arms.

**Written.** This entry was written at 05:26Z on 23 Nov 2026, from the design fixed at 04:03Z, the lists, the pilot and the box's checks. No comparison result had been read.

## C5 design

**Block.** Plan block C5 asks whether a search over each learner's settings, made inside the training folds, beats the release's rule, which already chooses some settings inside each fit; and whether a call cutoff chosen on inner scores beats 0.5. The design was fixed at 04:08Z on 23 Nov 2026, before any code ran on real data, and `docs/settings-arms.md` gives the method. No comparison result has been read.

**Arms.** Two arms stand against the default, the release's rule refitted in the same run with the same binary, both behind `grm learn --exp-settings`. The search: per learner and outer fold, 31 trials (the release's own setting and 30 random settings from the spaces frozen in `scripts/hyperopt/protocol.py`, the same 30 everywhere, from the sampler seed 20261030), each scored by the mean balanced accuracy of 3 inner folds inside the outer fold's training genomes; the default-anchored rule at γ = 0.04 picks, a tie goes to the lowest trial id, and the pick is refitted once on the whole training set. The cutoff: the release's own fit, untouched, with its calls made at the point of the grid 0.05 to 0.95 that gives the best pooled inner balanced accuracy, a tie going to the point nearest 0.5 and then to the lower. The lasso's space is one integer, so its 30 draws hold 28 distinct settings.

**Folds.** The inner folds are the release's own tuning partition of each outer training set, grouped by lineage, so no lineage cluster is split. On the study split they keep every lineage whole but can split a BioProject, since `grm learn` holds no BioProject per genome. The outer test genomes are unused in the inner run and reach no fit; the outer run starts only after the picks and cutoffs are written. Each inner fit redoes the screen on its own fitting genomes.

**Measure and count.** As C3's: per learner, the lineage balanced accuracy is primary; a dataset's value is the mean of seeds 42, 43 and 44; the twin is the default at seeds 45, 46 and 47. Holm runs over 15 lineage pairs, the search, the cutoff and the twin by 5 learners. The conditions are read with the lead's C3 rulings 15 and 16 and condition 3 as corrected at 02:08Z on 23 Nov 2026, each of the twin's two counts read apart. The cutoff arm's scores do not move, so its log loss equals the default's and decides nothing. SCM's result is read beside the synthetic loss of 0.0126 (0.0088 to 0.0163) the plan names.

**Study.** Each of the 55 first-run datasets has a reviewed `study-whole` plan and gets one study run per arm at the release's seeds, with no twin, as C3 and C9 did (question 8).

| Set | Datasets | Species | Less M, the first run |
|---|---:|---:|---:|
| S and L | 69 | 9 | 55 in 8 species |

**Questions.** The design leaves fourteen questions to the lead, each with the answer taken and the other one, each one argument of the code or the run script. The design writes no other answer for questions 4, 7 and 14. The lead ruled all fourteen at 05:52Z on 23 Nov 2026, with no C5 score read: each as taken, with two additions. Question 3: the default-anchored rule at 0.04 is the arm read, and the plain best (the highest inner mean, ties to the lowest trial id) is refitted on the outer training genomes and scored on the outer test fold as an arm beside it, labelled "plain best, no anchor", in no Holm family and deciding nothing; where it equals the anchored pick it takes that pick's result and no refit runs. Question 13: a cut is decided on cost alone, before any score is read, and named in this entry.

| Question | Taken | Other |
|---|---|---|
| 1. The 30 settings | sampler seed 20261030, the same 30 everywhere, the release's setting as trial 0 beside them (31 trials) | a seed per dataset; 30 trials with trial 0 among them |
| 2. Inner folds | the release's tuning partition by lineage; on the study split, by lineage inside each study training fold | BioProject groups joined with their copy groups |
| 3. Rule and ties | default-anchored, γ = 0.04, ties to the lowest trial id; ruled: the plain best refitted beside it, in no family | pooled inner balanced accuracy |
| 4. What the release chooses | each sampled setting fixes every value of its space; the release's inner choices run only in trial 0 | - |
| 5. The screen | redone in each inner fit on its own fitting genomes | the outer fold's screen |
| 6. Cutoff | grid 0.05 to 0.95, pooled inner balanced accuracy, ties nearest 0.5 then lower | every midpoint between inner scores |
| 7. Log loss of the cutoff arm | the default's, deciding nothing | - |
| 8. Study draws | one run at the release's seeds, no twin | three draws |
| 9. Holm | 15 pairs, the twin's included | a family per arm |
| 10. Conditions | ruling 15 and corrected condition 3 | a per-learner reading |
| 11. Default and twin | C5's own, in the same run | shared with C3's run under one binary |
| 12. First run | S and L less M | all 69 |
| 13. Cut order past 46 xlarge-hours | the costliest learner's study search, then its seeds 43 and 44; ruled: on cost alone, before any score is read | the lead rules the rest |
| 14. Test families | the 15 lineage pairs; study split and other measures in no family | - |

**As built.** Two departures from the design, neither a change of method. The driver runs `grm learn` once per learner for the inner trials and once per learner for the outer refit, in place of one run for all five, so each learner's time and memory are measured apart. The outer settings file gives no grid, so the outer run writes calls at 0.5 and at the chosen cutoff only.

**Pilot.** *Pseudomonas aeruginosa* with piperacillin/tazobactam (210 genomes) and with aztreonam (235 genomes), seed 42, both splits, five learners, the default and both arms, on one 16-CPU, 32 GB box, 8 outer folds at once with 2 threads each: 04:43:02Z to 04:51:49Z and 04:52:59Z to 05:02:30Z on 23 Nov 2026. A fold took 274.05 s on average, under 1,890,000 KB at peak. Per fold, a learner's 31 inner trials took 39.71 s (rf) to 57.07 s (scm) and its outer refit 7.22 s to 11.47 s. 200 of 200 `grm learn` runs ended 0, the 20 choices files share one shape, and none names a genome. No score was opened.

**Run script.** `scripts/c5-settings/run_species.sh` was rehearsed on the pilot's box with *Pseudomonas aeruginosa* and piperacillin/tazobactam alone, seed 42 and twin seed 45: 05:09:33Z to 05:28:53Z on 23 Nov 2026, exit 0. The seed-42 report equalled the release's in 12 of 12 learner rows, the scorer's 323 checks passed with 7 of 7 results present, and no kept file names a genome. A second start skipped its 12 finished steps and 10 finished drives and wrote the same list of output files; a start with another seed list was refused and left `DONE` in place.

**Faults.** The 14 plants of 11 faults, each run as its own test on a copy of the tree, were all caught on 23 Nov 2026, the brief's eight among them: a setting picked on the outer test fold, an inner fold that splits a lineage, the screen fitted on test genomes, the cutoff chosen on outer scores, a tie broken by the outer score, the default refitted with a searched setting, the twin scored as an arm and a missing draw read as a zero. A first run caught 13 of 14; the test of the scorer's twin-seed check then failed for a second reason too, and now fails only on that check.

**Bytes.** With no new flag, the branch's outputs equal main's on the pilot: 6 of 6 files equal once the recorded wall-clock seconds are set to 0 (the only fields that differ between two runs of one binary), `grm predict` on `test_dna.fa` equal byte for byte, and the seed-42 report equal to the release's in 12 of 12 learner rows. `test_dna.fa` feeds only the AMRFinderPlus fixture, so the pilot's runs carry the proof for `grm learn`.

**Cost.** A fold costs its pilot wall a release-second of learning time, 8.16, times the dataset's release seconds, and never less than ten matrix loads; 1,100 folds over the 8 species, 8 at once, 6 for *Escherichia coli* and *Klebsiella pneumoniae*, whose fits hold more memory. With the default, the twin, the data steps and the fixed steps, the first run comes to 28.07 xlarge-hours central and 80.41 upper, against the plan's 46. The search is 19.61 central, the cutoff 2.24, the default 2.03, the twin 1.71 and the steps every species runs 2.49. The plain best adds one fit inside the outer run of each learner whose plain best is neither trial 0 nor the pick. Taken on every learner of every fold, it is half the outer runs' share of a fold's wall (0.1498 in the pilot): 1.64 xlarge-hours central and 4.91 upper, or 3.27 and 9.82 if each costs a whole outer run. With it, the first run is 29.71 central and 85.32 upper. The central cost is under 46, so nothing is cut. *Klebsiella pneumoniae* is 14.41 central and *Escherichia coli* 9.44. The pilot's datasets are small, so the upper is three times the central for every fit.

**Written.** This entry was written at 05:12Z on 23 Nov 2026, from the design fixed at 04:08Z, the lists and the pilot; the run script and fault lines were added at 05:48Z. No comparison result had been read.

## Mechanism split, result

**Reading.** Blocks M1 and M2 were read once on 23 Nov 2026, from 08:02:26Z to 08:02:28Z, over the second run of the 21 species of A at main `11c2996`, with the table script and its tests committed at `6372979` before any table of the run was opened (`docs/mechanism-split-result.md`). The reading met three things the method page did not foresee. None changes a figure of the M1 or M2 reading, and nothing ran again.

**A rule in two parts.** The rule of sulfamethoxazole/trimethoprim, `C:TRIMETHOPRIM & C:SULFONAMIDE`, is the only rule of the 237 with two parts. `strata.py` reads a genome's rows that meet either part (`amr_rule.rule_rows`), and the program's call needs both (`amr_rule.program_call`). So a resistant genome with a row of one part alone stands in gene only or point only, and the program does not call it. On the primary set the program calls 38,673 of the 39,814 gene-only resistant genomes and 7,959 of the 7,963 point-only ones, and the 1,141 and 4 it does not call are all on sulfamethoxazole/trimethoprim, in 8 of its 9 datasets of the primary set. No sulfamethoxazole/trimethoprim dataset holds 10 or more resistant genomes in both strata, so the M1 reading does not change. The page's strata counts include them. Whether a stratum of a two-part rule needs both parts is the maintainer's to rule.

**Two names for one class.** In Staphylococcus aureus the program's Class column names fusidic acid both `FUSIDIC ACID` and `FUSIDIC_ACID`. The cell rule, a class as the Class column names it, makes two cells: among the voting genomes `FUSIDIC ACID` has a point share of 0.000 and `FUSIDIC_ACID` a gene share of 0.000, in each period of both pairs. Neither moves. No other two of the 32 class names differ only by a space or an underscore.

**The scan over the record.** `scan.py` knows a run's tables by name. Over the record's 20 files it found 0 genome ids, 0 BioProject accessions, 0 sequences and 0 positions, and it flagged each of the record's 15 tables as a table it does not know, so it returned 1. The method page counts a table the scan does not know as a hit, and it did not foresee the scan run over the result's own tables. The 15 flags are by name alone. The floors hold on the run tables the 15 are made from (gate 5 of `page_tables.py`).

**Two sentences of the method page.** A check of the result against the run's scripts found two places where the method page's words and the script at main `11c2996` differed, and the script is what ran. The page listed six guards and the script ran seven: the one left out is the lineage clusters voting in both periods, and it moves 0 cells. The page put every resistant later test genome under the time-split link's share, and the script takes the share over those with a catalogue call: 13 of the 24 datasets hold such genomes with no call, and on 2 of the 13 the two shares differ. The link is exploratory and decides nothing. Both pages now state what ran, with the date. No figure changes.

**Written.** This entry was written at 08:18Z on 23 Nov 2026, after the reading, and the paragraph on the method page's two sentences at 09:17Z.

## C7 design

**Block.** Plan block C7 asks whether learners from outside GRM do better than GRM's own on the same folds. The design was fixed at 07:37Z on 23 Nov 2026, before any code ran on real data, and `docs/outside-learners.md` gives the method. The lead ruled its fourteen questions at 09:05Z on 23 Nov 2026 (Questions, below). No comparison result has been read.

**Arms.** Three arms stand against the release's rule as plan block C3's run refitted it. XGBoost 3.0.5 and LightGBM 4.6.0 fit the 4,000 rows `grm export` screens for each outer fold and are read against GRM's boosted trees (`gbt`); a scikit-learn 1.7.2 elastic net fits plan block C3's unitig build and is read against GRM's lasso (`lasso`). Each wheel is pinned by its sha256 on pypi.org in `scripts/c7-outside/requirements.txt`, the learner cross-check's file plus LightGBM 4.6.0. The arms are scripts. The one change under `src/` is `grm export --seed S`, 42 by default: with no flag every file is the bytes main wrote at `0d987ba`, and at seed S the plans are the ones `grm learn --seed S` scores, fold by fold, with no key added. The default and the twin are the C3 run's, read on each species' C3 machine (question 11).

**Folds.** The folds are GRM's: on the lineage split the plan `grm export --seed s` writes for the draw's seed, with its sha256 checked; on the study split the reviewed `study-whole` fold file's outer folds, with each fold's boosting holdout and 3 inner folds drawn by GRM's rule, other genomes than GRM's own study fit drew (question 2). Before any fit the run checks every export's lineage plan at seed s against the C3 default's folds at seed s by class, and the seed-42 input fingerprint against the release report's, and stops on any difference. Early stopping reads the plan's 15% boosting holdout, and the elastic net's mix and penalty are picked on the 3 inner folds. No outer test genome enters a screen, a holdout, an inner fold or a pick, and a test holds each of those apart.

**Measure.** Per arm, the lineage balanced accuracy is primary, with recall, specificity, both error rates with Wilson intervals, and log loss beside it. Each arm runs three lineage draws a dataset, at seeds 42, 43 and 44 (question 1), and is read against the default's mean at the same seeds, `k31-at-arm-seeds`. The twin is the default at seeds 45, 46 and 47.

**Count.** Better, worse or the same by 0.02 at 12 decimals, the arm's value less the default's, the median difference with its species interval, and the one-sided sign test. Holm runs over 5 lineage pairs: `xgboost-gbt`, `lightgbm-gbt`, `elasticnet-lasso`, `twin-gbt` and `twin-lasso`. An arm carries the change when its pair meets conditions 1 and 2; condition 3 reads both twin pairs for the whole block, and each arm's own learner's twin pair is written beside it, deciding nothing (question 14); 4 and 5 are read per arm. `score.py --count` refuses before all three arms and the twin have lineage rows, unless an arm is named with `--withdrawn`. The count without the datasets where an elastic net refit did not converge is written beside, `count-no-unconverged.tsv`, deciding nothing (question 9).

**Study.** The first-run datasets in S, 55, get one study run per arm at seed 42, with no twin.

| Set | Datasets | Species | Less M, the first run | First run in S |
|---|---:|---:|---:|---:|
| L | 268 | 23 | 233 | 55 |

**Questions.** The design leaves fourteen questions to the lead, each with the answer taken and the other one. The lead ruled all fourteen at 09:05Z on 23 Nov 2026, with no C7 score read: questions 4 and 11 changed, the rest as taken, five of them with a limit named on the page or a count written beside, and one addition.

| Question | Taken | Other | Ruled |
|---|---|---|---|
| 1. Seeds of the arms | a flag `grm export --seed S`, default 42, with no byte of any output moved without it; not written until ruled | seed 42 alone for good, a departure from the plan's Rules, Draws | as taken, and written: three draws, at seeds 42, 43 and 44; two Rust tests (no flag: main's bytes at `0d987ba` on the export fixture's meropenem and ciprofloxacin datasets; `--seed 43`: the plans `grm learn --seed 43` scores, fold by fold); `docs/dataset-export.md` says so; no file it writes is read by `tests/givens_export` or the Givens export |
| 2. Study plans | GRM's rule for the holdout and inner folds on the reviewed fold file's outer folds, the screen by the export's rule | `grm export --fold-file F` | as taken; the page names it: the inner folds and the holdout are other genomes than GRM's own study fit drew |
| 3. The boosters' rows | the export's forest screen on all training genomes | the chi-square screen again on the fitting genomes | as taken; the page names the limit: the boosters' stopping holdout took part in choosing the 4,000 rows and GRM's `gbt` holdout did not, so a difference reads the screen and the library together |
| 4. Settings with no GRM knob | the library's default, LightGBM's `min_data_in_leaf` 20 among them | 1 | changed: `min_data_in_leaf` 1. `min_data_in_bin` stays the library's 3: a SYNTHETIC test shows a row held by 1 or 2 of 34 fitting genomes is still split on. Both are named as deviations, and the plain-booster check holds at a floor of 1 |
| 5. XGBoost's tree method | `exact`, the learner cross-check's | `hist` | as taken |
| 6. A booster with one tree | kept, scored and counted | scored as GRM's null fit | as taken |
| 7. The elastic net's input | the chi-square top 4,000 of the fold's unitig candidates on each fit's own genomes | every candidate | as taken |
| 8. Column scale | each column over its standard deviation on the fitting genomes | raw 0 and 1 | as taken |
| 9. An unconverged elastic net refit | scored and counted | the dataset stops | as taken, and the count without the datasets with an unconverged refit written beside, `count-no-unconverged.tsv`, deciding nothing |
| 10. Holm | 5 lineage pairs | 8, the twin's five learners | as taken |
| 11. Default and twin | C7's own, in the same run with the same binary | the C3 run's, on each species' C3 box | changed: the C3 run's, a named deviation. `grm` only exports in C7; the C3 run fitted the default and the twin on the same 233 datasets and held seed 42 to the release; fitting them again costs 33.45 box-hours central and 93.11 upper. The run refuses unless every export's lineage plan at seed s equals the default's folds at seed s by class and by genome (each fold's genome count and the sha256 of its sorted ids, against the `grm learn` report the C3 run wrote at seed s), the seed-42 input fingerprint equals the release report's, and the C3 run ended `DONE` at `8bfd4ca` with its `SHA256SUMS` holding. The page lists `git diff 8bfd4ca..HEAD -- src/` with the proofs on main that a run without the flag writes the same bytes. C7 opens nothing under a C3 run's `score/` |
| 12. First stage | `scripts/c3-inputs/run_species.sh` with the unitig arm, as it is | a C7 script without C3's unitig fits | as taken, for a new box kept for a lost C3 machine; on a C3 machine nothing of stage 1 runs |
| 13. First run | L less M | all of L | as taken |
| 14. Condition 3 | both twin pairs for the whole block | each arm's own learner's twin pair | as taken, and each arm against its own learner's twin pair (`gbt` or `lasso`) written beside, deciding nothing |

**Addition.** `score.py --count --withdrawn ARM`: a withdrawn arm keeps its pair at p = 1 in the 5-pair family and is named in `count.tsv`. The count refuses a withdrawn arm that has results, and the twin or the default named as withdrawn. The elastic net runs unless the box limit says no; `run_species.sh` takes `ENET`, `run` by default or `leave-out`, which lists and fits no elastic net, and a test holds it. A run started with `leave-out` restarts on the same OUT with `run`, keeps its finished steps and booster fits, and refuses any other change to its lists. The count refuses a dataset `lists.tsv` marks M and prints how many of the 233 first-run datasets it holds. The cost with and without the elastic net is below.

**Pilot.** One lineage draw at seed 42 of *Klebsiella aerogenes* with ertapenem (83 genomes, 48 resistant), imipenem (68, 26) and meropenem (78, 26), the default and each arm, one thread a fit, on one 16-CPU box. The nine fits ran from 08:10:13Z to 08:21:55Z on 23 Nov 2026, and the elastic net as built again from 08:23:00Z to 08:34:27Z. Their reports went to a sealed folder that no one opened. CPU seconds a draw: XGBoost 0.70, 0.87 and 0.71; LightGBM 1.23, 0.85 and 1.03; the elastic net 667.37, 600.08 and 686.13, each under 320,000 KB at peak. The runner's screen equals the export's on 5 of 5 lineage folds in each dataset. A scan found 0 genome ids in the 47 files that came back and 0 in the rerun's 12. After the rulings the pilot ran again at seeds 42, 43 and 44, from 09:24:25Z to 09:38:01Z: `grm` built at the branch head in 42.38 s, the export at each seed in 3.41 to 3.55 s, the default at seeds 43 and 44 as the C3 run ran it, `plans.py` with 21 of 21 checks passing over 9 dataset seeds, the boosters at the three seeds with LightGBM's floor at 1, and the elastic net at 43 and 44. The head's export at `--seed 42` and with no flag equals `8bfd4ca`'s in 12 of 12 files. Mean CPU seconds a draw over the 9 fits of each arm: XGBoost 0.84, LightGBM 2.17 and the elastic net 644.57. The 24 fits ended rc 0, their reports sealed, and a scan found 0 genome ids in the 119 files.

**Cost.** The plan gives 15 xlarge-hours plus the elastic net after a pilot. At 2,000 genomes a SYNTHETIC elastic net fold took 603.8 s; a power law through it and the 9 pilot fits (exponent 0.4728) gives 3,019.0 s a draw, plus the projection. At three lineage draws over the 233 datasets and one study run on each of the 55, on the species' C3 machines, C7's steps take 74.25 box-hours central and 250.54 upper, $14.85 and $50.11 at xlarge's $0.20 an hour. The elastic net is 70.44 and 244.17 of that; without it the boosters, the exports, the build and the libraries take 3.81 and 6.38, $0.76 and $1.28. On a new box a species, with stage 1's 33.45 and 93.11, the run comes to 107.69 and 343.65, or 37.26 and 99.49 without the elastic net. At one draw before the rulings C7's steps were 31.61 and 119.89 on the C3 machines.

**After the review.** The lead ruled two points at 11:26Z on 23 Nov 2026, with no C7 score read. saga's stop stays at tol 1e-4 and 2,000 iterations: a tighter stop would raise the elastic net's 370.92 CPU hours, and the page names the limit, which every elastic net result carries. A missing draw stays a failed check that stops its species' score, as the scorer has it and as blocks C3, C4 and C5 have it. The design's words left only the dataset unscored; the code holds, so no dataset can drop out of a count unseen, and a species with a fit missing is run again or ruled on before any of its scores is written.

**Written.** This entry was written at 08:51Z on 23 Nov 2026, from the design fixed at 07:37Z, the lists, the pilot and the SYNTHETIC timings. The rulings of 09:05Z, the three-seed pilot and the cost at three draws were added at 09:39Z. No comparison result had been read, and no pilot's report was opened.

## C2 design

**Block.** Plan block C2 asks whether fitting more than 2,000 genomes a species helps, scored on the same held-out genomes. The design was fixed at 08:13Z on 23 Nov 2026, before any code ran on real data, and `docs/cap-arms.md` gives the method. No comparison result has been read. Nothing in `src/` changes, so the release's output bytes stay the same, and the release keeps its cap of 2,000 genomes a species.

**Arms.** Three arms stand against the release's 2,000, refitted in the same run: the top 4,000 and the top 4,096 by the release's selection rule on the rebuilt revision 14 table, and a twin of the default. The ranking does not depend on the cap, so the cohorts nest, and the release's 2,000 lie in the top 4,000 in all 8 species. The wider arms add 2,000 and 2,096 genomes in 7 species. *Enterococcus faecium* ranks 2,907, so its two wider arms are one cohort with 907 added genomes, fitted once.

**Folds.** The lineage split runs through `grm arm` with `common_outer_test_set_balanced`: every arm is scored on the default's outer test folds, the same in every arm, and an added genome only trains, never in a fold whose scored genomes share its lineage or its k-mer duplicate group. The study split is `grm learn --fold-file` on the default's `study-whole` test folds; an added genome trains in every fold or in none, and only when its BioProject is known and no tested genome has it or its k-mer digest. `check_folds.py` holds every fit to this on the run box and stops the species on a breach.

**Measure.** As C3's: per learner, the lineage balanced accuracy is primary, with recall, specificity, both error rates with Wilson intervals, and log loss beside it. A dataset's value is the mean of 3 draws, `grm arm` slots 0 to 2 at evaluation seed 42; the twin is the default at slots 3 to 5.

**Count.** C3's count on a family of its own: `cap-4000`, `cap-all` and the twin against `cap-2000`, Holm over 15 lineage pairs, 3 arms by 5 learners, taken once after every species has run. `score.py --count` refuses part of the family and a row given twice. Study rows are one fit each, read beside the count.

**Datasets.** The CAP datasets: 166 in the 8 species, 23 of them in M, so 143 in the first run, 39 of those in S. Of the 143, 115 depend on none of the four label rules and 28 on the most resistant call against the other three.

**Study.** `master-genome-meta.jsonl.gz` holds no BioProject for the genomes added in *Escherichia coli*, *Enterococcus faecium*, *Klebsiella pneumoniae* and *Streptococcus pneumoniae*. It holds one for 1,694 of the 2,000 added at 4,000 in *Mycobacterium tuberculosis*, 1,652 in *Salmonella enterica* and 1,137 in *Staphylococcus aureus*. Without question 4, a wider arm adds no study genome on 33 of the 39 datasets in S.

**Pilot.** *Enterococcus faecium* with vancomycin, ampicillin and gentamicin, one lineage slot of `cap-2000` and `cap-all`, scores moved unread, ran on a 16-CPU box from 08:23:24Z to 08:36:46Z on 23 Nov 2026, after a restart for the runner fault below. The fetch took 439.58 s for 2,907 genomes with no failure; 363 release genomes differ in bytes from the master list, none in length, contig count or N50. The `cap-all` count took 129.53 s, its matrix 105.81 s at 3,206,332 KB peak, and each fit 38.86 to 54.05 s at most 3,655,908 KB. Every fold check passed, with no balance drop and no duplicate in training. Over the 5 outer folds, the 901 added genomes labelled for vancomycin trained in 3,644 fold places and were kept out by lineage in 861; for ampicillin, 188 genomes, 753 and 187. No added genome carries a gentamicin call, so on gentamicin `cap-all`'s fit repeats the default's (question 10).

**Fault.** The first pilot fit stopped on a group file named `1000`: the runner's variable was `GROUPS`, which bash sets to the user's group ids. It is now `GROUPING`, and a runner test fails on any bash special variable.

**Checks.** `tests/c2_cap` holds 62 SYNTHETIC tests. The ten faults of `docs/cap-arms.md`, three of them planted in two places, 13 in all, are each caught by the one test named for it, and one test holds each fault's text to one occurrence in its file at this tree.

**Cost.** Scaled from the pilot: a fit costs 1.7873 s of `grm arm` wall a second of the release's learning, times its genomes over 2,000. On the species' refit boxes, which hold the release's 2,000 FASTA files, the run comes to 53.83 xlarge-hours central and 157.76 upper, $10.77 and $31.55; on new boxes, 54.50 and 158.77, $10.90 and $31.75. Per arm, central and upper hours: `cap-2000` 8.83 and 26.50, the twin 7.49 and 22.46, `cap-4000` 17.52 and 51.99, `cap-all` 18.35 and 54.42; the data steps every arm shares 0.99 and 1.41; the fetch 0.65 and 0.98 on a refit box, 1.33 and 1.99 on a new one. Nothing beyond the pilot has run.

**Questions.** The design leaves ten questions open, each with the answer taken and the other one.

| Question | Taken | Other |
|---|---|---|
| 1. Path | `grm arm` for every cohort, the default refitted, no change to `src/` | a `grm learn` flag adding training genomes to the release's drawn split, so that C3's default and twin serve |
| 2. Draws | `grm arm` slots 0 to 2, the twin 3 to 5, at evaluation seed 42 | seeds 42 to 44 and 45 to 47, which `grm arm` does not take |
| 3. Study default | fitted here | C3's study files of the default, held to their sha256 |
| 4. BioProjects of added genomes | none fetched; an added genome with no known BioProject stays unused | the BV-BRC genome records of the added genomes fetched by id |
| 5. *Enterococcus faecium* | one cohort of 2,907, fitted once, counted under both wider arms and marked | counted under `cap-all` only |
| 6. 4,000 and 4,096 | both run | one dropped; they differ by 96 genomes in 7 species |
| 7. Final models for block T1 | none built | a final fit per dataset on a passing arm's cohort |
| 8. Shell helpers | the runner's own `step` and `scrub` | a shared file after C3 finishes |
| 9. Fits at once | one dataset at a time; a pilot fit used 3.8 to 4.7 CPU seconds a wall second of 16 | four datasets at once, cutting the fits' wall about fourfold |
| 10. A wider arm with no added genome labelled | fitted; its lineage fit repeats the default's | marked none-added and not fitted, as the study split does |

**Written.** This entry was written at 08:55Z on 23 Nov 2026, from the design fixed at 08:13Z, the lists and the pilot. No comparison result had been read.

## C2 rulings

**Rulings.** The ten questions of the C2 design were ruled at 09:22Z on 23 Nov 2026, before any fit of the full run. No comparison result had been read.

| Question | Ruling |
|---|---|
| 1. Path | as taken: every arm, the default included, is fitted by `grm arm`, so the default's figures here are not the release report's figures and are never set beside them as if they were |
| 2. Draws | as taken, a deviation from the plan (below) |
| 3. Study default | as taken: fitted here |
| 4. BioProjects of added genomes | the other: the added genomes' BV-BRC genome records fetched by genome id, the genome table only, `genome_id` and `bioproject_accession` only |
| 5. *Enterococcus faecium* | the other: its one wider cohort counted under `cap-all` only |
| 6. 4,000 and 4,096 | `cap-4000` withdrawn before any fit (below) |
| 7. Final models for block T1 | as taken: none built |
| 8. Shell helpers | as taken: the runner's own `step` and `scrub` |
| 9. Fits at once | the other: up to four datasets' lineage fits at once, under a memory guard |
| 10. A wider arm with no added genome labelled | the other: none-added, not fitted and out of the count's denominator, on both splits |

**Draws.** The plan names `grm learn` seeds 42 to 44 for the draws and 45 to 47 for the twin. `grm arm` takes no seed of its own; its keyed slots are its draws, each with its own repeat seed, partitions and learner draws. So the draws are the default's slots 0 to 2 and the twin its slots 3 to 5, all at evaluation seed 42. This deviates from the plan's seeds, because the lineage split's pairing, `common_outer_test_set_balanced`, runs on `grm arm` only.

**Withdrawn.** `cap-4000` was withdrawn at 09:22Z on 23 Nov 2026, before any fit. It differs from `cap-all` by 96 genomes in 7 species and is the same cohort in *Enterococcus faecium*, so one result would be paid for twice. The family stays at 15 pairs: `score.py --count --withdrawn cap-4000` keeps its five lineage pairs at p = 1 and names it withdrawn in `count.tsv` and `conditions.tsv`. The runner builds no `cap-4000` root; the arm's code and tests stay.

**BioProjects.** `bioprojects.py` asked `https://www.bv-brc.org/api/genome/` on 23 Nov 2026, from 09:45:37Z to 09:45:53Z, for the 2,096 `cap-all` added genomes of each species with a dataset in S, and every one was answered. With a BioProject: *Escherichia coli* 985, *Klebsiella pneumoniae* 2,006, *Mycobacterium tuberculosis* 1,771, *Salmonella enterica* 1,748 and *Staphylococcus aureus* 1,137. Where the master genome rows also name one, 4,530 genomes, both name the same. Before the k-mer digest check, an added genome may train on 27 of the 39 datasets in S, 1,978 in all.

**None-added.** `cohorts.py` now writes, per dataset and wider arm, with the species and antibiotic in full, the added genomes with a label by class. Under `cap-all`, 27 of the 143 first-run datasets have none (`docs/cap-arms.md` names them). No threshold is set.

**Fits at once.** No limit in `grm arm`'s fit path depends on wall time (`docs/cap-arms.md`). On the pilot's data, 8 fits run four at once wrote the same files as each fit run alone, 56 of 56, with 80 `seconds` values set aside. The runner at this text then ran the pilot's three datasets on the pilot's binary at `AT_ONCE` 4: vancomycin's chain alone, then ampicillin's and gentamicin's at once. Gentamicin's `cap-all` was none-added and not fitted, and the other 5 fits' 20 files equal those of the pilot's fits, which ran one at a time. N is 4 for *Enterococcus faecium*, *Neisseria gonorrhoeae*, *Staphylococcus aureus* and *Streptococcus pneumoniae*; 3 for *Escherichia coli*, *Mycobacterium tuberculosis* and *Salmonella enterica*; 2 for *Klebsiella pneumoniae*.

**Cost.** On new boxes, with no `cap-4000`, no fit of a none-added dataset and lineage fits at once: 37.59 xlarge-hours central and 130.26 upper, $7.52 and $26.05, in box lifetimes of 133.25 hours ($26.65). The lineage fits take 27.37 hours central, against 60.70 one at a time. Fits are costed at the proof's alone pace, 1.49 to 2.09 times the pilot's wall. The upper is past the plan's 66 to 99.

**Checks.** The four faults this round adds, one of them planted in two places, are each caught by the one test named for it: 14 faults, 18 plantings in all.

**Written.** This entry was written at 09:59Z on 23 Nov 2026, after the rulings, from the fetch, the cohort counts and the two checks on the pilot's data. No comparison result had been read.

## C2 review

**Review.** An independent review of the C2 branch at `7d64574` recounted, with programs of its own, the cohorts by the selection rule, the study split, the BioProject answers, the none-added datasets, the scorer's count on SYNTHETIC draws and the cost. Every recount equals the branch's figures. The review read no comparison result and no score.

**Refused draws.** `grm arm` writes a lineage trace line's and prediction row's `repeat` as the slot, while the slot's plan holds the keyed repeat, which differs once a refused draw is replaced. `check_folds.py` and `score.py` read by the keyed repeat, so one replaced draw would have stopped the species. Both now read by slot. When the replacements are refused too, the scripts stopped the species, where H-1's ruling leaves the dataset's lineage runs not produced and lets the run go on. They now follow that ruling (`docs/cap-arms.md`, Refused draws), a rule the lead confirms before the run.

**Limits.** The page now gives the study split's counts per dataset behind the 1,978 added genomes that may train on 27 of the 39 datasets in S, the 2,138 tested genomes with no BioProject known, and the 19 fitted datasets that add 1 to 19 labelled genomes, 105 in all.

**Tests.** Seven of the review's 13 planted faults passed every test: genomes of one score sorted otherwise than by id as text; three breaches of the lineage fold check that no test reached (a scored genome in training, a default genome that trains otherwise than in the default, a cohort plan whose lineage count differs); a draw's counts not held to its prediction rows; and the runner's memory rule, held to its text only. Each now has a SYNTHETIC test, and with the refused-draw tests `tests/c2_cap` holds 95.

**Written.** This entry was written at 10:49Z on 23 Nov 2026, during the review. No comparison result had been read.

**After the review.** The lead ruled two points at 12:19Z on 23 Nov 2026, with no C2 score read. H-1's ruling on a refused repeat holds for C2: a dataset whose lineage draws are refused past both replacements leaves the lineage count of every arm and the twin alike, the result page names it with its reason, and its study rows are scored as before. Every arm is scored on the default's test folds, so such a dataset has no lineage test for any arm, and leaving it out favours none. The 19 fitted datasets that add 1 to 19 labelled genomes stay in `cap-all`'s count of 116: the count was fixed with no threshold before the cohorts were listed, and a threshold set now would be chosen with the lists in hand. The result page gives the count over the other 97 datasets beside it, deciding nothing.

## Label rule

**Plan.** The plan marks 36 of its 274 datasets M: "M marks a dataset holding a pair with S and I or with I and R". S is a susceptible call, I an intermediate one and R a resistant one. Its Labels paragraph says: "The 36 wait for the maintainer's rule." Step 1 of its Order says: "The M datasets wait for the maintainer's label rule."

**Pages.** Sentences like these two, that datasets of M wait for the label rule or for the maintainer's label rule, stand on 16 method and result pages, each for its own datasets in M.

**Ruled.** On 24 Nov 2026 the maintainer left the plan's open rules to the lead. The lead ruled the label rule that day, with no dataset of M fitted for any block. The release's rule stays. It is `r-or-s-wins` in `docs/input-configuration.md`: S with I is S, I with R is R, R with S has no label.

**Holds.** The 36 datasets of M are open for fits under `r-or-s-wins`. Each claim on them has two readings: one on every scored genome, and one without the genomes whose records for that antibiotic hold an intermediate call, with no refit. A claim on a dataset of M stands only if it holds both ways.

**Later records.** Block T1's Datasets row says the label rule picked applies to the later-records test too, so that test uses `r-or-s-wins`. Its claims on datasets of M take the two readings as well.

**Why.** `intermediate-aside` leaves out the borderline genomes, so the task would read easier than it is in practice. `most-resistant` changes labels outside M as well: by the plan, 33 of the 238 hold a pair with R and S. `intermediate-class` makes the task three classes. Every merged result uses `r-or-s-wins`.

**Counts with M.** The lead ruled three more points at 17:09Z on 25 Nov 2026, with no dataset of M fitted for any block.

- The plan's sentence "Each result names the rules it does not depend on" binds on M too. A result on a dataset of M names them, and says so when it depends on all four.
- A block's first-run verdict stands. No arm that missed conditions 1 to 5 on the first run comes to meet them by a later count.
- For a candidate of the entry Condition 6 the rule runs one way. The release changes to it only if conditions 1 to 5 also hold on the count over all of its block's datasets, M included, read both ways.

**Three counts.** When a block's datasets of M have run, its result page gives three counts: the first run's as it stands, the count over M alone, and the count over all its datasets. The block's rules, scorer and arguments stay unchanged.

**Not run.** No block of the plan has fitted a dataset of M. Until one has, a page's sentence that its datasets of M wait for the label rule means: the lead has made the rule, and the run has not started.

**Read again.** Three sentences say a line is read again when the label rule is made: two in `docs/refit-result.md` and one in this page's Refit reading. That reading is D0's, over the 10 MISS datasets in M. It needs their refit, which has not started.

**Written.** This entry was written at 16:49Z on 25 Nov 2026. After its independent review the lead added, at 17:12Z, the plan's words on M, the paragraphs Pages, Counts with M, Three counts and Read again, and the last sentence of Later records. No dataset of M had been fitted. The maintainer can overrule it.

## Condition 6

**Pages.** The plan's sixth condition for a release change is: "The maintainer rules the change." The result pages of blocks C2, C3, C4, C5, C7 and C9 each leave it to the maintainer.

**Ruled on 24 Nov 2026.** The lead ruled condition 6 for the comparison blocks then read, under the maintainer's word of that day.

- Blocks C2, C4, C5 and C9: no arm meets conditions 1 and 2, so there is no change to rule on, and the release keeps its rule.
- Block C3: the AMRFinderPlus program's calls alone meet conditions 1 to 5 with `cart`, `rf`, `gbt` and `lasso` (`docs/input-arms-result.md`). That arm becomes a candidate for the one release change below. The release does not change now.

**Beside those verdicts.** C2's `cap-4000` was withdrawn before any fit and has no reading. C4's `cmh-16000` meets conditions 1 and 2 on `gbt` and `lasso`, 2 of the 3 learners it touches, and the conditions ask all three. C3's verdict is read by C3's ruling 15, fixed before any C3 score was read. No arm meets the stricter reading beside it, which decides nothing.

**Ruled on 25 Nov 2026.** The lead ruled condition 6 for block C7's three arms at 14:51Z, after the one reading of block C7 (`docs/outside-learners-result.md`) and before its independent review ended.

- `xgboost`, XGBoost 3.0.5, is not a candidate. It meets conditions 1 to 4 and not 5: its study median is -0.002 on 55 study datasets, against a lineage median of 0.003. The release changes only when all six hold.
- `lightgbm`, LightGBM 4.6.0, and `elasticnet`, the scikit-learn 1.7.2 elastic net, meet conditions 1 to 5 and become candidates for the one release change. The release does not change now.

**C7's review.** It ended after the ruling. It traced the figures of the page's changed lines and of the pull request's text to the record, 1,537 with a program of its own and the other 5 by hand, and found 0 mismatch. It made 12 findings, and none changes a verdict. Commit 96860bbc answers six of them, on the page, in its tests and in the workflow's comment. The other six concern the pull request's text, a commit message or working files outside the repository, or leave the page as it stands.

**One release change.** The release changes at most once, after the block that follows C7. The plan does not hold that block: a later design adds it. Its plan is written and reviewed, and its entry stands on this page, before the lead reads its first number.

**The pick.** The plan holds no pick among arms that each met conditions 1 to 5, and no limit of one release change. The ruling of 24 Nov adds both, after blocks C2, C3, C4, C5 and C9 were read. The lead makes the pick among today's release, C3's candidate, C7's two candidates and that later block's candidates, by the plan's rules on held-out scores, and writes it down before the later-records test.

**Before the pick.** Each candidate's block first takes its count over all its datasets, M included, as the entry Label rule gives it. A candidate that does not meet conditions 1 to 5 there, read both ways, is not picked.

**Later-records test.** Block T1's Scored row stands as the plan wrote it: "Once, after every block is ruled: the release's rule and each arm that met conditions 1 to 5, on the same genomes and measure". Today those arms are C3's candidate and C7's two.

**A narrower ruling, taken back.** The ruling of 24 Nov named only the picked release and today's release for that scoring. That would have left out an arm that met conditions 1 to 5 and was not picked. The lead ruled at 16:47Z on 25 Nov 2026 that the plan's row holds, since the plan was written before any result.

**Added to T1.** The test also scores the picked release when it is not one of those arms. It scores a candidate of the later block only if that block's design entry, written before the lead reads its first number, says so and gives the candidate's rule. The test's page reports an arm it cannot score as not scored, with its reason. The later records decide no pick.

**Beside C7's candidates.** The pick weighs these, from `reference/outside-learners/`. The lineage medians are 0.003 for `lightgbm` and 0.006 for `elasticnet`, both inside the plan's band for "same" (0.02). By more than that band, `lightgbm` is better on 62 of 232 lineage datasets and worse on 14, and `elasticnet` better on 71 and worse on 19.

**Twin.** The plan reads every count against the twin, the default fitted again with only the fold seeds changed. On the same 232 datasets `twin-gbt` is better on 40 and worse on 33, and `twin-lasso` better on 28 and worse on 19.

**Costs.** On the study split `lightgbm`'s log loss is worse, by a median of 0.010. `elasticnet` kept no weight in 221 of 3,490 lineage folds, and its median fit took 5210.252 wall seconds against 4.878 for `lightgbm`. Either arm adds a library to the release. The lead rules on adding a dependency at the pick.

**Holds.** These rulings change nothing in the release and no sentence on the site. The site's sentences on the study split change only with the one release change.

**Written.** This entry was written at 16:49Z on 25 Nov 2026, after the readings of blocks C2, C3, C4, C5, C7 and C9 and with no later record read. After its independent review the lead changed, at 17:12Z, the paragraphs C7's review, One release change, The pick and Added to T1, and added Beside those verdicts, Before the pick and Twin. The maintainer can overrule each ruling in it.
