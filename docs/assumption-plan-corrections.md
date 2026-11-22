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
