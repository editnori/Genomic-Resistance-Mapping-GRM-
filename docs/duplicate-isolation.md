# Duplicate isolation in the outer folds (revision 10)

From revision 10, the random split keeps each group of k-mer-identical genomes whole in one outer test fold, so no genome is scored by a fit that trained on its own copy. The GRM admin ruling of 2026-10-30 set the design, taken from the 2026-10-27 duplicate-leakage lane (`data/research/grm-research-dup-leakage-all23-r1`, independently reviewed twice).

## Why

Before revision 10 nothing kept copies apart. Stratified random folds could put one copy of a genome in training and the other in test. The 2026-10-27 lane reran the evaluation with copies isolated:

- *S. flexneri* ampicillin, azithromycin and ciprofloxacin lost 0.04 to 0.13 balanced accuracy. Every correct-to-wrong flip fell on a genome whose copy had been in training.
- Ciprofloxacin VME went from 0.00 to 0.30 on 33 resistant genomes.
- Elsewhere the effect was near zero: same-BioSample duplicates moved the median by 0.000, and lineage folds already hold copies together.

## The rule

- **Group.** The samples of one dataset whose genomes have the same kset digest (`matrix.json` `kset_digests`, from `counted.json`), meaning the same canonical k-mer set. Groups of one are not listed. Same-BioSample or same-assembly genomes with different k-mer sets are not grouped; GRM does not load those identifiers, and the lane measured no effect from them.
- **Move.** The random split starts from the same stratified folds as before (`stratified_folds`, seed 42). Each group moves whole into the fold holding most of its members. A tie goes to the fold of the group's lowest-index member among the tied folds. A tie broken by the lowest fold index would pile split pairs into fold 0 (the lane's first run did, and was discarded).
- **Scope.** Only the outer random split changes. Inner tuning folds, the lineage split, MIC cross-validation folds and the final models are unchanged, and inner tuning may still split a group: the lane found duplicate-grouped tuning cost about 0.013 balanced accuracy for reasons unrelated to leakage. Identical genomes have identical sketches, so they always share a lineage and the lineage folds already hold each group; the evaluator refuses a group that spans lineages.
- **No duplicates, no change.** A dataset without a group gets the folds, scores and report bytes it got before.

The code is `learn::isolate_groups` (the move), `eval::outer_splits` (where it applies) and `analysis::kset_duplicates` (the groups).

## Report field

`kset_duplicates`, at the top level of a dataset report: one ascending list of sample indices per group, ordered by first member, groups disjoint. It is left out when a dataset has no group. The fold a group landed in is `splits[].test_fold` of any member.

Consumers that check exact report keys must accept it. The MTB 4k plan's validator contract (thread 234bc2b6, `report_contract.py`) must mirror it.

## Checks

- `scripts/check-results.py` checks the field's shape (ascending, disjoint, in range), that each group shares a lineage, and that each group sits in one test fold of every split that ran. Without counted.json it cannot tell "no groups" from "no field". It also refuses a root that mixes evaluation revisions, so a revision-10 root must hold only revision-10 reports.
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

The other 17 species, *M. tuberculosis* included, have none. The species check changes some cohorts, and the groups are rebuilt from each run's inputs, so the revision-10 counts can differ slightly.

On the revision-8 learners, isolating *S. flexneri* ciprofloxacin moved random-split VME from 0.000 to 0.24 to 0.30 across the five learners and balanced accuracy down by 0.11 to 0.16. The three *S. haemolyticus* datasets, which have no group, gave byte-identical results.

## Limits

- Duplicates are only as complete as k-mer identity. Near-identical genomes (a few SNPs apart) and same-patient isolates are not grouped; the lineage split addresses relatedness.
- MIC cross-validation and inner tuning can still split a group.
