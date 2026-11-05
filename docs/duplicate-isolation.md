# Duplicate isolation in the outer folds (revision 10), inner tuning (revision 11) and final-model tuning (revision 12)

From revision 10, the random split keeps each group of k-mer-identical genomes whole in one outer test fold, so no held-out classification score in a report's `splits` comes from a fit that trained on the genome's own copy. From revision 11, inner tuning keeps each group whole too (see The rule); revision 10 tuned the random split without groups. From revision 12, the final models tune on the duplicate groups as well; revision 11 tuned them by class. Only the MIC regression still draws its own folds (see Limits). Main adopted the revision-12 change on 5 Nov 2026 under Layth's delegation of that day (DECISIONS F5, UD-7). The GRM admin ruling of 2026-10-30 set the design, taken from the 2026-10-27 duplicate-leakage lane (`data/research/grm-research-dup-leakage-all23-r1`, independently reviewed twice).

## Why

Before revision 10 nothing kept copies apart. Stratified random folds could put one copy of a genome in training and the other in test. The 2026-10-27 lane reran the evaluation with copies isolated:

- *S. flexneri* ampicillin, azithromycin and ciprofloxacin lost 0.04 to 0.13 balanced accuracy. Every correct-to-wrong flip fell on a genome whose copy had been in training.
- Ciprofloxacin VME went from 0.00 to 0.30 on 33 resistant genomes.
- Elsewhere the effect was near zero: same-BioSample duplicates moved the median by 0.000, and lineage folds already hold copies together.

## The rule

- **Group.** The samples of one dataset whose genomes have the same kset digest (`matrix.json` `kset_digests`, from `counted.json`), meaning the same canonical k-mer set. Groups of one are not listed. Same-BioSample or same-assembly genomes with different k-mer sets are not grouped; GRM does not load those identifiers, and the lane measured no effect from them.
- **Move.** The random split starts from the same stratified folds as before (`stratified_folds`, seed 42). Each group moves whole into the fold holding most of its members. A tie goes to the fold of the group's lowest-index member among the tied folds. A tie broken by the lowest fold index would pile split pairs into fold 0 (the lane's first run did, and was discarded).
- **Scope.** The outer random split moves groups, from revision 10. From revision 11, inner tuning keeps them whole too: the random split's inner folds and tuning holdouts are drawn on duplicate groups (`eval::tuning_groups`), and the lineage split's on lineages, as before. Revision 10 drew the random split's inner folds and holdouts by class, without groups, so a tuning partition could put a genome and its copy on both sides. A dataset with a group can therefore score differently at revision 11, and every revision-10 report is recomputed rather than reused, because the input fingerprint does not see code changes. A grouped holdout holds out whole groups, so boosting holds out about 1/7 of training and the lasso about 1/5, not 15% and 20% of each class. The lineage outer split and the MIC regression's folds are unchanged, so the MIC regression can still split a group (see Limits). From revision 12, the final models (`eval::final_models`, fitted on every sample) tune on the duplicate groups too: scm and cart's inner folds and the boosting and lasso holdouts are drawn on the same ids the random split's inner tuning uses (`eval::final_tuning`, which is `tuning_groups("random", ..)`), whatever split the dataset's lineages allow. Lineages play no part in final tuning. Revision 11's final models drew their inner folds by class (`stratified_folds`) and their holdouts without groups, so the published hyperparameters could be tuned with a genome and its copy on both sides of a partition. A dataset with a group can therefore get different final models at revision 12 (rules, depth, lasso weights, importances), though no held-out score changes, and every revision-11 report is recomputed rather than reused. A final grouped plan that leaves an inner fold without one class fails the dataset (`invalid final-model inner evaluation`), as an invalid random split does. The 2026-10-27 lane measured duplicate-grouped tuning at about 0.013 lower balanced accuracy, for reasons unrelated to leakage. Identical genomes have identical sketches, so they always share a lineage and the lineage folds already hold each group; the evaluator refuses a group that spans lineages.
- **No duplicates, no change.** A dataset without a group gets the folds, scores and final models it got before; across each bump only the revision its report records changes. Over every sample, `training_folds` without groups is `stratified_folds`, so revision 12's final folds equal revision 11's there (`final_tuning_keeps_each_duplicate_group_whole` checks it).

The code is `learn::isolate_groups` (the move), `eval::outer_splits` (where it applies), `eval::tuning_groups` and `eval::final_tuning` (the tuning groups) and `analysis::kset_duplicates` (the groups). The tests `final_tuning_keeps_each_duplicate_group_whole` and `final_models_tune_on_the_duplicate_groups` fail if a duplicate group straddles a final inner fold or holdout, or if `evaluate`'s final models are not the ones fitted on those partitions.

## Report field

`kset_duplicates`, at the top level of a dataset report: one ascending list of sample indices per group, ordered by first member, groups disjoint. It is left out when a dataset has no group. The fold a group landed in is `splits[].test_fold` of any member.

Consumers that check exact report keys must accept it. The MTB 4k plan's validator contract (thread 234bc2b6, `report_contract.py`) must mirror it.

## Checks

- `scripts/check-results.py` checks the field's shape (ascending, disjoint, in range), that each group shares a lineage, that each group sits in one test fold of every split that ran, and that the report is from revision 10 or later. Without counted.json it cannot tell "no groups" from a missing field; that is the audit's job. It rejects an explicit empty list, which grm never writes and the audit treats as absent. It also refuses a root that mixes evaluation revisions, so a revision-12 root must hold only revision-12 reports.
- `scripts/metric-audit/audit.py`, from `DUPLICATE_ISOLATION_REVISION` (10), rebuilds the groups from `counted.json`'s kset digests, requires the report to list exactly those, and rebuilds the random folds with the same move before comparing `test_fold`. Earlier revisions must list no group.

## Scale at revision 10

Over the 277 datasets of the 23-species candidate, 58 have at least one group and 55 have a group that crossed the saved random folds. Every group is a pair.

| Species | Datasets with groups | Groups | Crossing saved random folds |
|---|---:|---:|---:|
| *S. flexneri* | 4/4 | 534 | 436 |
| *S. sonnei* | 6/20 | 1,482 | 1,223 |
| *E. coli* | 21/37 | 40 | 36 |
| *K. pneumoniae* | 22/39 | 38 | 33 |
| *P. aeruginosa* | 3/13 | 3 | 2 |
| *K. aerogenes* | 2/3 | 2 | 1 |

The other 17 species, *M. tuberculosis* included, have none. These counts describe the revision-9 cohorts. None of the genomes the species check drops is in a group, but replacements refill capped species, so recount after the revision-10 curation; the audit rebuilds each report's groups from counted.json either way.

On the revision-8 learners, isolating *S. flexneri* ciprofloxacin moved random-split VME from 0.000 to 0.24 to 0.30 across the five learners and balanced accuracy down by 0.11 to 0.16. The three *S. haemolyticus* datasets, which have no group, gave byte-identical results.

## Limits

From the independent review of `9dad7f7..2b24e50` (approved 2026-10-30; report `grm-r10-dup-isolation-review-20261030/REVIEW.md`, outside the repository):

- **A move can empty a class from a random test fold, and that fails the dataset.** Before isolation, `min_class >= folds` guaranteed every stratified fold both classes. Moving groups whole breaks that guarantee. A random split that fails its plan check is fatal for the dataset (`check_split_plan` in `evaluate` and `lasso_only`); only a failed lineage split becomes `unavailable_reason`. It does not fire on today's cohorts. The tightest is *S. flexneri* ampicillin, where moves cut one fold's susceptible genomes from 8 to 4. A smaller cohort with that duplicate rate would stop the dataset with `invalid random outer evaluation`. The test `a_move_that_empties_a_class_from_a_random_fold_fails_the_dataset` reproduces it. If it happens in a run, the run's standing rule applies: resume once, then document the dataset as a skip. Marking the random split unavailable, or refusing such a cohort at curation, would be a later revision's change.
- **The MIC regression still splits pairs.** `evaluate_mic` draws its own stratified folds, so a genome's MIC can be predicted by a fit that trained on its copy. On the revision-9 cohorts, 33 datasets have both MIC data and a group. Of the 58 pairs that carry MICs, the MIC folds split 48. That is at most 2 genomes per dataset in cohorts of 176 to 1,012, so essential agreement moves by well under a point. `mic_essential_agreement` is published.
- **Identical matrix columns are not always identical ksets.** Two genomes can differ only in k-mers dropped as rare, so the model sees them as identical. Across the six species with groups, kset identity covers 415 of the 416 identical-column groups. The exception is one *P. aeruginosa* pair with different R/S labels, which stays split in amikacin and levofloxacin (2 genomes in each).
- **Moves reshape class stratification.** In each *S. flexneri* dataset 107 to 114 genomes move, and *S. sonnei* trimethoprim moves 204. That is expected, and it is where the first limit comes from.
- **Duplicates are only as complete as k-mer identity.** Near-identical genomes (a few SNPs apart) and same-patient isolates are not grouped; the lineage split addresses relatedness.
