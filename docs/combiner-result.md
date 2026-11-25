# Combiner result

Can a combination of the five learners call the resistant genomes their mean misses? On the learner overlap's page (`docs/learner-overlap-result.md`), the second pass counts the resistant genomes the five learners' mean misses at the release's 0.5 cutoff, split by how many of the five miss each one. This page asks whether a weighted mean of the five, or a logistic combiner on their five scores, recovers some of the misses that only some of the five make, without calling more susceptible genomes resistant than it gains.

It reads the stored held-out scores of the release's lineage split and fits only the two combiners, on the genomes of the other folds. It trains no learner and scores no genome with a learner again.

## Plan

This plan follows FREEZE-COMBINE and was committed with its script and tests, before any reading.

- **Reports.** The 268 revision 14 reports the learner overlap pooled: the five learners on the lineage split and 10 or more resistant genomes. The read checks that its reports are the pooled rows of `reference/learner-overlap/datasets.tsv`, one for one.
- **Inputs.** For each genome, its label and the five learners' stored held-out scores on the lineage split, with the split's stored `test_fold`. No catalogue call, lineage, species, year or country enters a fit.
- **The incumbent.** The mean of the five with equal weights, called resistant when the five scores sum to 2.5 or more by `math.fsum`, as the second pass calls it.
- **The weighted mean.** Five weights, none negative, summing to one.
- **The logistic combiner.** An intercept and one coefficient for each of the five scores, with no interaction term.
- **Fitting.** Each combiner is fitted per report and per stored outer fold on the genomes of the other folds, by plain mean log loss with every score clipped to [0.0001, 0.9999], and applied to the held-out fold. The weighted mean is fitted by projected gradient from the equal weights; the logistic combiner by Newton's method from zero, with no penalty.
- **Cutoff.** Every form calls a genome resistant at a score of 0.5 or more.
- **Measures.** Per report and form: the counts at the cutoff, the balanced accuracy, recall, specificity and very major error rate worked from them, and the mean log loss. The very major error rate is the resistant genomes called susceptible over the resistant genomes, as the release's `vme`.
- **Decision.** An arm is kept only if the median over the reports of its balanced accuracy less the incumbent's is 0.02 or more, the species-clustered interval of that median lies wholly above 0, and the species-clustered interval of the median change in its very major error rate does not lie wholly above 0. Otherwise the reading is no difference shown.
- **Interval.** 4,000 resamples of species, the median of each resample's reports, and the 101st and 3,900th of the sorted medians. Each row's resamples are drawn from a generator seeded by the sha256 of the row's own name. A set of under 6 species gets no interval.
- **Predictions, counted and not judged.** The reports where each arm gains 0.02 or more, and the counts better, worse and the same at that margin. Of the resistant genomes the incumbent misses, how many each arm calls, split into those 1 to 4 of the five call and those none of the five calls. The weighted mean cannot call a genome none of the five calls, since its score lies within the five. Of the susceptible genomes the incumbent calls resistant, how many each arm does not, and how many susceptible genomes each arm newly calls resistant.
- **Beside.** Recall, specificity and log loss stand beside the decision and change nothing. The page gives the reports pooled and by species, and a gain on one species changes nothing.

## Reading

Not read yet.

## Limits

**Words.** A lineage split is a retrospective split inside one database, never external validation. "No difference shown" is never written as "equal". The error levels serve as an adapted research screen and are never an applicable FDA validation criterion. Coverage is the share of genomes answered, not accuracy.

**One draw.** The read uses the release's own lineage folds, the only folds whose per-genome held-out scores are stored. The learner overlap's page counted these scores, labels and folds before this plan was written.

**Rows.** The record bars a row per genome; a row per report is an aggregate and keeps its genome count. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

## Method

**Reading.** The commands, with the release's report folder as REPORTS, the reading's folder as READING and the record as OUT:

```
python3 -I scripts/mechanism/combiner.py read REPORTS --first reference/learner-overlap --second reference/learner-overlap/second-pass --out READING
python3 -I scripts/mechanism/combiner.py write READING --reports REPORTS --out OUT
python3 -I scripts/mechanism/combiner.py print OUT
```

Every key named `id` is dropped as a report is parsed. Gate 8 holds the reports' genome ids in memory only, searches every written file for each one in full, and counts tokens of a genome id's shape beside it. It prints counts only.

**Gates.** Each prints counts and stops the run if it fails:

| Gate | Stops the run when |
|---|---|
| 1. reports | the learner overlap's checks fail, or the reports read differ from the pooled rows of `reference/learner-overlap/datasets.tsv` |
| 2. cross-fitting | a fold's training genomes hold one of its held-out genomes, or the two do not add up to the report |
| 3. incumbent | the incumbent's counts, worked again, differ from `reference/learner-overlap/second-pass/vote.tsv` |
| 4. convex bound | the weighted mean calls a genome none of the five calls or misses one all five call, or a fold's weights are negative or do not sum to 1 |
| 5. convergence | a fit does not converge within its cap |
| 6. cutoff | a form leaves a genome without a call |
| 7. rows | the record holds a row per genome, or an arm's counts do not follow from the incumbent's and the transitions |
| 8. ids | a written file holds a genome id, a token of an id's shape, or an `id` or `samples` key, or a decimal cell is 10 or more in size |

**Tests.** `tests/mechanism/test_combiner.py` holds 40 SYNTHETIC tests of the script. `tests/mechanism/test_combiner_faults.py` runs `scripts/mechanism/combiner_faults.py`, which plants each fault below in the script in turn and runs the test meant to catch it. `tests/mechanism/test_combiner_page.py` holds 9 tests of this page against the script.

| Fault planted | Test that catches it |
|---|---|
| C1. Each fold's arm fitted on every genome, its own held-out genomes among them | `Crossfit.test_gate_2_counts_no_shared_genome_on_the_standard_read` |
| C2. Gate 2 blind to a fit on its own held-out genomes | `Gates.test_gate_2_stops_a_fit_on_its_own_genomes` |
| C3. A score of exactly 0.5 read as a miss | `Forms.test_a_score_of_exactly_the_cutoff_is_a_call` |
| C4. The arms called at the training prevalence (0.425 in the test's design) in place of 0.5 | `Fits.test_the_logistic_fits_each_groups_share` |
| C5. The majority of the five read in place of their mean | `Forms.test_the_incumbent_is_the_mean_not_the_majority` |
| C6. The weighted mean fitted on a class-weighted log loss | `Fits.test_the_weighted_mean_fits_plain_log_loss` |
| C7. The weights not brought back to a sum of one | `Fits.test_the_weights_stay_on_the_simplex` |
| C8. Gate 4 blind to weights that do not sum to one | `Gates.test_gate_4_stops_on_weights_off_the_simplex` |
| C9. The weighted mean free to round across the bound of the five scores | `Forms.test_the_weighted_mean_stays_inside_the_five_scores` |
| C10. A logistic fit that does not converge let through | `Gates.test_gate_5_stops_a_fit_that_does_not_converge` |
| C11. Log loss not clipped | `Forms.test_log_loss_is_clipped` |
| C12. The very major error rate worked over the susceptible genomes | `Forms.test_measures_are_worked_from_the_counts` |
| C13. Accuracy written in place of coverage | `Forms.test_coverage_is_the_share_given_a_call` |
| C14. Coverage written in place of the balanced accuracy | `Forms.test_balanced_accuracy_is_not_coverage` |
| C15. Gate 6 blind to a genome given no call | `Gates.test_gate_6_stops_on_a_genome_given_no_call` |
| C16. The interval from resamples of reports, not species | `Pooled.test_the_interval_resamples_species` |
| C17. Every row's interval drawn from one seed | `Pooled.test_the_interval_is_seeded_by_the_rows_name` |
| C18. An arm kept with a gain interval that does not lie above 0 | `Rule.test_the_gain_interval_must_lie_above_0` |
| C19. A median gain of exactly 0.02 read as under the margin | `Rule.test_the_gain_must_reach_the_margin` |
| C20. An arm kept with a higher very major error rate | `Rule.test_a_higher_vme_blocks_the_arm` |
| C21. A higher very major error rate counted as better | `Pooled.test_better_worse_and_same_at_the_margin` |
| C22. A genome one of the five calls counted with those none calls | `Counts.test_the_transitions_split_by_how_many_of_the_five_call` |
| C23. Arm counts that do not follow from the transitions let through | `Gates.test_gate_7_stops_on_arm_counts_that_do_not_add_up` |
| C24. Each genome's call kept in the reading | `Ids.test_the_reading_holds_no_row_per_genome` |
| C25. Gate 7 blind to a list as long as a report | `Gates.test_gate_7_stops_on_a_row_per_genome` |
| C26. A genome id kept as a report is parsed | `Ids.test_no_id_is_held_after_parsing` |
| C27. Gate 8 blind to a genome id with no id shape | `Gates.test_gate_8_stops_on_a_genome_id_in_the_reading` |
| C28. Gate 8 blind to an id-shaped token | `Gates.test_gate_8_stops_on_an_id_shaped_token` |
| C29. A decimal cell of 10 or more let through | `Gates.test_a_decimal_cell_stays_under_10` |
| C30. A report set that differs from the first record's let through | `Gates.test_gate_1_stops_on_a_report_set_that_differs` |
| C31. An incumbent that differs from vote.tsv let through | `Gates.test_gate_3_stops_on_an_incumbent_that_differs_from_vote_tsv` |

All 31 were caught. Each test passes as the script stands and fails with its fault planted.
