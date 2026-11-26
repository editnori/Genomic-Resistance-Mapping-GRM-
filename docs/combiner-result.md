# Combiner result

Can a combination of the five learners call the resistant genomes their mean misses? On the learner overlap's page (`docs/learner-overlap-result.md`), the second pass counts the resistant genomes the five learners' mean misses at the release's 0.5 cutoff, split by how many of the five miss each one.

This page asks two things. Does a weighted mean of the five beat that mean on the same split without raising the very major error rate? And how many of the misses that only some of the five make does it recover? A logistic combiner on their five scores is read beside it, as descriptive only.

It reads the stored held-out scores of the release's lineage split and fits only the two combiners, on the genomes of the other folds. It trains no learner and scores no genome with a learner again.

## Plan

This plan follows FREEZE-COMBINE and was committed with its script and tests, before any reading.

### Reports

- **Reports.** The 268 revision 14 reports the learner overlap pooled: the five learners on the lineage split and 10 or more resistant genomes. The read checks that its reports are the pooled rows of `reference/learner-overlap/datasets.tsv`, one for one.
- **Report count.** Gate 1 stops the read unless it reads exactly those 268 reports, and the command line gives no way to change the count. Only a SYNTHETIC test hands in its own count.
- **Folds.** A report with one fold, or a genome with no fold number, stops the read with a line that names the report by species and antibiotic.
- **Inputs.** For each genome, its label and the five learners' stored held-out scores on the lineage split, with the split's stored `test_fold`. No catalogue call, lineage, species, year or country enters a fit.

### Forms

| Form | Score |
|---|---|
| Incumbent | The mean of the five with equal weights, called resistant when the five scores sum to 2.5 or more by `math.fsum`, as the second pass calls it |
| Weighted mean | Five weights, none negative, summing to one |
| Logistic combiner | An intercept and one coefficient for each of the five scores, with no interaction term |

Every form calls a genome resistant at a score of 0.5 or more.

### Fitting

Each combiner is fitted per report and per stored outer fold on the genomes of the other folds, by plain mean log loss with every score clipped to [0.0001, 0.9999], and applied to the held-out fold. The weighted mean is fitted by projected gradient from the equal weights. The logistic combiner is fitted by Newton's method from zero, with no penalty.

This stop rule was ruled on 25 Nov 2026 (R19), after the offline figures were read.

| Fit | Ends |
|---|---|
| Weighted mean | When a full backtracking search finds no lower loss, and that point is its answer, or at its cap of 5,000 steps, where it stands |
| Logistic combiner | At the first of three: convergence by its tolerance, 100 Newton steps, or a step with no finite solution, where the last finite iterate stands |

The weighted mean starts from the equal weights and only ever lowers the loss, so its fitting loss is never above the loss at the equal weights. A fold ended at a cap or with no finite step is counted, and its report is marked for that arm. A held-out score that is not finite stops the read.

- **Logistic stop values.** The logistic fit starts with the intercept and the five coefficients at 0. It converges when no entry of the mean gradient is above 1e-10 and its last step is small.
- **Small step.** A step is small when it moved none of the six numbers by more than 1e-8 × (1 + the largest of the six in size).
- **Logistic step search.** A Newton system with a pivot under 1e-13 has no finite solution. The search halves a Newton step until the loss is no higher; once the step is under 1e-12 of its full length, it takes it whatever the loss.

### Score columns

- **One score.** A learner that gives every fitting genome of a fold one score gets a coefficient of exactly 0 and is left out of that fold's fit. This and the next rule were ruled on 26 Nov 2026 (R30), before any reading.
- **The same scores.** Learners that give every fitting genome of a fold the same scores are fitted as one column and share its coefficient in even parts. A held-out genome is scored from each of them.
- **No tolerance.** Both rules read the stored scores as they are, with no tolerance, and both can apply in one fold. Of the fits with the lowest loss, this one has the least sum of squares of the five coefficients.
- **No new number.** That fit is where a penalty on the sum of squares of the five coefficients ends as the penalty shrinks to 0. So it adds no free number, and its calls and coefficients do not depend on the learners' order.
- **Order.** How a fit ends can. On the reviewer's made-up reports, putting the learners in another order in 700 fits moved 0 of 42,000 calls and no coefficient by more than 3.3e-8, but moved 7 of 700 ends between the cap and converged.
- **Neither case.** On a fold with neither case, the fit is the one before 26 Nov, to the bit. A fold with either case is not marked; the record counts, per report, the logistic folds with a column set aside or shared.
- **First step.** If a logistic fold's first Newton step has no finite solution, that fold is not scored, and the read stops at gate 5, which names each such report by species and antibiotic. No fold is scored from coefficients left at 0.

### Counts

- **Three counts.** The count over every report and the count without the arm's flagged reports decide: the weighted mean is kept only if the rule keeps it in both. Otherwise its reading is no difference shown.
- **Few species.** A count of under 6 species has no interval, so it keeps nothing, and the weighted mean's reading is no difference shown, whatever the other count gives.
- **Second count.** Beside them, a second count without the arm's marked reports is written, and it decides nothing. The record gives each count's own outcome for the weighted mean and the arm's reading.

A report is flagged for an arm by any of three conditions, each read on the stored values with no tolerance:

| Condition | Flags the report when | Arms |
|---|---|---|
| One score (R30) | some learner gives every genome of a fold of 2 or more genomes one score | both |
| Coefficient under 0 (R37) | some fold's logistic fit gives one of the five learners a coefficient under 0; a coefficient of exactly 0 does not flag | logistic only |
| Training share (R46) | some learner gives 2 or more genomes of one fold that fold's training share of resistant genomes as their score, by `==` | both |

**Training share.** The release scores a learner that selected nothing at its fold's training share of resistant genomes. That share falls as the held-out fold's own share rises, so it carries the held-out labels, upside down. R46, ruled on 26 Nov 2026 before any reading, adds a condition for it.

**Sign.** A score that follows the training share points the wrong way for the held-out fold. With every coefficient at 0 or above, the share moves a fold's calls nearly always the wrong way: on the reviewer's made-up reports, 3 of 248 fold fits moved them the right way. A learner within a hair of the share gives no fold one score and no genome the share, so neither R30's condition nor R46's flags it. The weighted mean's weights are never under 0.

**Logistic arm.** R43, ruled on 26 Nov 2026 before any reading, reads the logistic arm as descriptive only. Its three counts, medians and intervals are written as before, and its outcome and reading say "descriptive only" in each count, whatever its figures. R37's condition, its column, the logistic arm's third count and the `flags:` line stay as counts, and decide nothing. The weighted mean is the one arm the rule can keep. It has no intercept and no weight under 0.

**Cost.** With the five learners' scores correlated at 0.8, and again at 0.95, R37's condition flagged 200 of 200 honest made-up reports, where no learner carries a fold's share. The logistic arm was then kept in 0 of 8 made-up readings whose true gain was near 0.05. With R30's condition alone, it was kept in 8 of 8.

**Break.** On the reviewer's made-up reports, a learner with real signal, whose score also falls as its training share rises, got the logistic arm kept through the fold's share alone, and R37's condition flagged 0 of its 20 reports. The same learners at each report's own share were not kept. With that learner, the weighted mean was not kept: its median gain in balanced accuracy was -0.0096 over 20 reports, on the reviewer's made-up reports.

**Weighted mean.** On the reviewer's made-up reports, a learner scored about 97 of 100 genomes of each fold at that fold's training share. The plain mean's median balanced accuracy fell by 0.0381 and 0.0320, and the weighted mean's by 0.0000. The weighted mean was kept in 10 of 10 readings of two such families, and in 0 of their 10 controls. R30's condition flagged 3 of 200 reports, as no fold had one score. R46 closes it with the training share condition above.

**Learners' code.** The release works a fold's training share by one division: resistant training genomes over training genomes, the report's genomes outside the fold. Only a boosted model with no tree, or a lasso with no weight, gives it as a score, to every genome of its fold. An empty tree or forest scores its fold 0.5. R30's condition flags both. No learner scores a part of a fold at the share. This is a reading of the code and no measurement.

**Class weights.** Four of the five, the decision tree, the forest, the boosted trees and the lasso, fit with class weights that give each class half the total weight. The rule learner gives a 0 or 1 call. Each fold's counts enter every fit of those four, and a leaf of fixed counts scores lower as the share rises. In 6 made-up families with a falling learner, two mimicking the class weights, the weighted mean was kept wrongly in 0 of 30 readings.

### Decision

- **Measures.** Per report and form: the counts at the cutoff, the balanced accuracy, recall, specificity and very major error rate worked from them, and the mean log loss.
- **Very major error rate.** The resistant genomes called susceptible over the resistant genomes, as the release's `vme`.
- **Decision.** The weighted mean is kept only if it passes the gain rule and the error rule below. Otherwise its reading is no difference shown. The logistic arm is descriptive only.
- **Gain rule.** The median over the reports of the arm's balanced accuracy less the incumbent's is 0.02 or more, and the species-clustered interval of that median lies wholly above 0.
- **Error rule.** The very major error rate is not higher on that interval. R19 confirmed the wording: the species-clustered interval of the median change in the very major error rate does not lie wholly above 0.
- **Interval.** 4,000 resamples of species, the median of each resample's reports, and the 101st and 3,900th of the sorted medians. A set of under 6 species gets no interval.
- **Seed.** Each row's resamples are drawn from Python's `random.Random`, seeded by the sha256 of the row's own name read as one whole number. The name is the arm and measure, as in `combiner|weighted|balanced_accuracy`, so an arm's three counts share one seed.
- **Beside.** Recall, specificity and log loss stand beside the decision and change nothing. The page gives the reports pooled and by species, and a gain on one species changes nothing.

### Predictions

The page counts three things for the predictions and judges none.

- **Gains.** The reports where each arm gains 0.02 or more, and the counts better, worse and the same at that margin.
- **Misses.** Of the resistant genomes the incumbent misses, how many each arm calls, split into those 1 to 4 of the five call and those none of the five calls.
- **Susceptible genomes.** Of the susceptible genomes the incumbent calls resistant, how many each arm does not, and how many susceptible genomes each arm newly calls resistant.

The weighted mean cannot call a genome none of the five calls, since its score lies within the five. The logistic combiner has an intercept and no such bound. The logistic combiner can call a genome none of the five calls only by moving what amounts to its operating point, which shows up as susceptible genomes called resistant.

**Predictions as written.** FREEZE-COMBINE wrote its first prediction before any reading:

> The weighted mean beats the equal mean by the 0.02 margin on **fewer than 1 in 5** of the reports, and not on the pooled median. About 1 in 4 that the pooled median gain clears 0.02 for either arm.

It restated its third on 25 Nov 2026, before any arm was read:

> Of the **8,067** resistant genomes the equal mean still misses, the better of the two arms recovers **under a third**, and it does not do so without giving back susceptible genomes among the 8,896 the equal mean already calls resistant.

The 8,067 of 93,786 resistant genomes are those the mean misses and 1 to 4 of the five call. The 8,896 of 121,533 susceptible genomes are those it calls resistant. Both are over 268 reports, in `reference/learner-overlap/second-pass/vote-pooled.tsv`. Each count is part of its class's total.

## Reading

Not read yet.

## Limits

**Words.** A lineage split is a retrospective split inside one database, never external validation. "No difference shown" is never written as "equal". The error levels serve as an adapted research screen and are never an applicable FDA validation criterion. Coverage is the share of genomes answered and is never accuracy.

**Per-fold inputs.** Each fold's stored scores come from base learners trained on that fold's own training genomes, so the five inputs are not one fixed model's outputs across folds. Each combiner is read as a reweighting of the release's own per-fold predictions, not as a model that could be shipped unchanged. The base learners are never refitted here. If the weighted mean is kept, shipping it would need a refit of the base learners inside one fold structure and the maintainer's word under R4.

**One draw.** The read uses the release's own lineage folds, the only folds whose per-genome held-out scores are stored. The learner overlap's page counted these scores, labels and folds before this plan was written.

**Caps.** On the reviewer's made-up reports, a logistic fit ended at its cap of 100 Newton steps had stalled at its best fit, and 0 of 4,000 calls differed from a second fit's. A fold at either cap, 5,000 steps for the weighted mean, still marks its report, and the count without marked reports decides nothing.

**Free sign.** The stored scores come from one layer of folds. On them, a combiner with a free sign and an intercept cannot be told from one that reads the fold's share, so the logistic arm is read as descriptive only. A coefficient under 0 also arises when two learners score nearly alike, with no fold's share in either.

**Incumbent.** A learner that carries the share lowers the incumbent, the mean with equal weights, so the weighted mean can gain beside it. A tree leaf can score exactly the share by chance, when the leaf's ratio of resistant to susceptible genomes is the square of the training set's. Such a report leaves the third count, and the `flags:` line counts it. A made-up learner within a hair of the share is not flagged.

**Rows.** The record bars a row per genome; a row per report is an aggregate and keeps its genome count. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

## Method

**Reading.** The commands, with the release's report directory as REPORTS, the reading's directory as READING, the sha256 the read printed as SHA256 and the record as OUT:

```
python3 -I scripts/mechanism/combiner.py read REPORTS --first reference/learner-overlap --second reference/learner-overlap/second-pass --out READING
python3 -I scripts/mechanism/combiner.py write READING --reports REPORTS --out OUT --sha256 SHA256
python3 -I scripts/mechanism/combiner.py print OUT
```

Every key named `id` is dropped as a report is parsed. Gate 8 holds the reports' genome ids in memory only, searches every written file for each one in full, and counts tokens of a genome id's shape beside it. It prints counts only.

**Gates.** Each prints counts and stops the run if it fails:

| Gate | Stops the run when |
|---|---|
| 1. reports | the learner overlap's checks fail, the reports read differ from the pooled rows of `reference/learner-overlap/datasets.tsv`, or they are not 268 |
| 2. cross-fitting | a fit is handed more or fewer genomes than its fold's training list, the list holds a held-out genome, or the two do not add up to the report |
| 3. incumbent | the incumbent's counts, worked again, differ from `reference/learner-overlap/second-pass/vote.tsv` |
| 4. convex bound | the weighted mean calls a genome none of the five calls or misses one all five call, or a fold's weights are negative or do not sum to 1 |
| 5. fit ends | a held-out score is not finite, or a logistic fold's first Newton step has no finite solution; it counts fold ends, marks, and columns set aside or shared |
| 6. cutoff | a form leaves a genome without a call |
| 7. rows | the record holds a row per genome, or an arm's counts do not follow from the incumbent's and the transitions |
| 8. ids | a written file holds a genome id, a token of an id's shape, or an `id` or `samples` key, or a decimal cell is 10 or more in size |

Before gate 1, a report with one fold or a genome with no fold number stops the read with a line of its own, `stop folds:`, naming the report. Gate 5 names in the same way each report with a fold whose first Newton step has no finite solution.

**Flags.** After the gates the read prints a line, `flags:`, with how many reports each condition flags, how many have both R30's and R37's, and how many R46's flags that R30's does not, and no score. Those counts decide nothing (R43).

**Write.** Write does not fit again. It takes the sha256 the read printed, works the sha256 of the reading's bytes before it parses them, and stops with one line, `stop sha256:`, on a difference or on a value that is not 64 hex characters. It then writes nothing. The record's `gates.txt` carries that sha256. Write does not work R30's flag or R46's out again: the sha256 covers every edit.

**Fold records.** Write works each report's sign flag, R37's condition, and its logistic mark out again from the logistic fold records it parses, and stops on a difference with a line, `stop flags:`, that names the report by species and antibiotic and gives no figure. A fold record that lacks a key write reads stops write with a line, `stop fold records:`, that names the report in the same way.

**Print.** Print reads the record's `SHA256SUMS` first. It stops with one line, `stop record:`, naming the file, when the list is missing, a file of the record is not listed, or a listed file is missing or its bytes give another sha256. Write prints the sha256 of the reading it read, so it can be matched with the read's.

**Arms.** Write stops on a reading that names any arm but the weighted mean as one the rule can keep, or that lacks R46's flag or names other flags. Write and print stop, with a line `stop rule:`, on a rule row that gives the logistic arm an outcome or a reading other than descriptive only. No option and no environment setting changes the arms the rule can keep.

**Tests.** `tests/mechanism/test_combiner.py` holds 115 SYNTHETIC tests of the script. `tests/mechanism/test_combiner_faults.py` runs `scripts/mechanism/combiner_faults.py`, which plants each fault below in the script in turn and runs the test named beside it. `tests/mechanism/test_combiner_page.py` holds 19 tests of this page against the script.

| Fault planted | Test that fails with it |
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
| C10. A logistic fold ended with no finite step counted as converged | `Gates.test_gate_5_counts_a_fold_with_no_finite_step` |
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
| C32. A held-out score that is not finite let through | `Gates.test_gate_5_stops_on_a_held_out_score_that_is_not_finite` |
| C33. A fold ended at its cap left unmarked | `Gates.test_gate_5_counts_and_marks_a_fold_at_the_cap` |
| C34. A weighted fold ended at its cap counted as converged | `Fits.test_a_fit_at_its_cap_ends_there` |
| C35. The marked reports kept in the second count | `Rule.test_the_second_count_leaves_out_the_marked_reports` |
| C36. The count without the marked reports deciding in place of the count over all | `Rule.test_the_second_count_leaves_out_the_marked_reports` |
| C37. The weighted fit started at the first learner's vertex in place of the equal weights | `Fits.test_the_weighted_fit_is_never_above_the_equal_mean` |
| C38. The record's rule for the very major error rate in other words than R19's | `Rule.test_the_record_says_the_vme_rule_in_r19s_words` |
| C39. A score column with one value over the fitting genomes fitted as a column of its own | `Fits.test_a_one_score_column_is_set_aside` |
| C40. Score columns the same on every fitting genome fitted apart | `Fits.test_columns_the_same_are_fitted_as_one` |
| C41. A shared coefficient given whole to the first column of its set | `Fits.test_a_held_out_genome_is_scored_from_both_columns` |
| C42. The same and one value read with a tolerance, not exactly | `Fits.test_the_same_and_one_score_are_read_without_tolerance` |
| C43. 1e-12 added to the diagonal of every Newton system | `Fits.test_a_fit_on_neither_case_is_the_fit_before_r30` |
| C44. A first Newton step with no finite solution ended as a later one is | `Gates.test_gate_5_stops_on_a_first_step_with_no_finite_solution` |
| C45. A fold whose first Newton step had no finite solution scored from its coefficients at 0 | `Fits.test_a_first_step_with_no_finite_solution_is_not_scored` |
| C46. Gate 5 blind to a first Newton step with no finite solution | `Gates.test_gate_5_stops_on_a_first_step_with_no_finite_solution` |
| C47. A fold with a score column set aside or shared marks its report | `Gates.test_the_reviewers_two_reports_pass_gate_5` |
| C48. A fold with a score column shared not counted with those set aside | `Gates.test_the_reviewers_two_reports_pass_gate_5` |
| C49. Gate 2 blind to a fit handed another count of genomes than its training list | `Gates.test_gate_2_stops_a_fit_handed_every_genome` |
| C50. The training count taken as the groups a fit is handed, not the genomes inside them | `Crossfit.test_each_fold_is_fitted_on_the_other_folds` |
| C51. The count without the flagged reports deciding nothing | `Rule.test_the_third_count_decides` |
| C52. A report flagged only when a learner gives the whole report one score | `Counts.test_a_report_is_flagged_by_one_fold` |
| C53. A fold of one genome flagging its report | `Counts.test_a_report_is_flagged_by_one_fold` |
| C54. The flagged reports kept in the third count | `Rule.test_the_third_count_decides` |
| C55. Write blind to a reading made under other counts of the rule or without its flags | `Gates.test_write_stops_on_a_reading_made_under_other_rules` |
| C56. An arm kept when either count that decides keeps it | `Rule.test_the_third_count_decides` |
| C57. Gate 1 blind to the count of reports it expects | `Gates.test_gate_1_holds_the_read_to_its_count` |
| C58. The command line holding the read to 267 reports, not the freeze's 268 | `Gates.test_the_command_line_holds_the_read_to_268` |
| C59. A genome with no fold number let through | `Gates.test_a_genome_with_no_fold_number_stops_the_read` |
| C60. A report with one fold let through | `Gates.test_a_report_with_one_fold_stops_the_read` |
| C61. No report flagged for a logistic fold fit that gives a learner a coefficient under 0 | `Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only` |
| C62. A coefficient of exactly 0, as a column set aside has, counted as under 0 | `Flags.test_a_column_set_aside_at_0_does_not_flag` |
| C63. A coefficient counted as under 0 only below -1e-9 | `Flags.test_a_coefficient_under_0_is_read_with_no_tolerance` |
| C64. The intercept counted among the five coefficients | `Flags.test_a_report_with_every_coefficient_at_least_0_is_not_flagged` |
| C65. The second condition read on the first fold fit only | `Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only` |
| C66. The logistic arm's third count leaving out R30's flagged reports only | `Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only` |
| C67. The weighted arm's third count also leaving out the reports flagged for a coefficient under 0 | `Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only` |
| C68. A count of 5 species given an interval | `Flags.test_a_logistic_third_count_under_6_species_gets_no_interval` |
| C69. The record counting the coefficients under 0, not the fold fits behind the flag | `Flags.test_the_record_counts_the_fold_fits_behind_the_flag` |
| C70. A report with either condition counted as one with both | `Flags.test_the_read_prints_how_many_reports_each_condition_flags` |
| C71. Write blind to a reading made without the second condition | `Flags.test_write_stops_on_a_reading_made_without_the_second_condition` |
| C72. Gate 5 naming no report whose first Newton step has no finite solution | `Flags.test_gate_5_names_each_report_with_no_finite_first_step` |
| C73. Gate 5 also naming a report whose later Newton step has no finite solution | `Flags.test_gate_5_names_each_report_with_no_finite_first_step` |
| C74. The reading's sha256 worked without its last line break | `Flags.test_read_and_write_print_the_readings_sha256` |
| C75. The logistic arm's outcome in each count taken from the rule | `Descriptive.test_the_logistic_arm_reads_descriptive_only_whatever_its_figures` |
| C76. The logistic arm's reading taken from its counts that decide | `Descriptive.test_the_logistic_arm_reads_descriptive_only_whatever_its_figures` |
| C77. The logistic arm's reports, species, medians and intervals not written | `Descriptive.test_the_logistic_arm_is_descriptive_only_in_each_count` |
| C78. Both arms named as arms the rule can keep | `Descriptive.test_one_line_names_the_arms_that_can_be_kept` |
| C79. The arms the rule can keep read from an environment setting | `Descriptive.test_no_option_or_environment_setting_changes_the_arms_that_can_be_kept` |
| C80. An option of read and write naming the arms the rule can keep | `Descriptive.test_no_option_or_environment_setting_changes_the_arms_that_can_be_kept` |
| C81. Write blind to a rule row that gives the logistic arm another outcome or reading | `Descriptive.test_write_stops_on_a_rule_row_that_gives_the_logistic_arm_another_outcome` |
| C82. Print blind to a rule.tsv that gives the logistic arm another outcome or reading | `Descriptive.test_print_stops_on_a_record_that_gives_the_logistic_arm_another_outcome` |
| C83. Write blind to a reading made under other arms that can be kept | `Descriptive.test_write_stops_on_a_reading_made_under_other_arms` |
| C84. Write not working the sign flag out again from the fold records | `Flags.test_write_works_the_sign_flag_out_again` |
| C85. Write not working the logistic mark out again from the fold records | `Flags.test_write_works_the_logistic_mark_out_again` |
| C86. Write's stop line naming the report by species only | `Flags.test_write_names_each_report_whose_flags_differ` |
| C87. Write naming only the first report whose flags differ | `Flags.test_write_names_each_report_whose_flags_differ` |
| C88. The flags line giving the sign's count where R30's belongs | `Flags.test_the_flags_line_gives_each_conditions_own_count` |
| C89. Gate 5 naming a report only with 2 or more folds whose first Newton step has no finite solution | `Flags.test_gate_5_names_a_report_with_exactly_1_fold_with_no_finite_first_step` |
| C90. Write taking a report's sign flag given as text | `Flags.test_write_stops_on_a_sign_flag_given_as_text` |
| C91. Write taking a fold fit's count given as true or false | `Flags.test_write_stops_on_a_fold_count_given_as_true_or_false` |
| C92. The training share worked as the fold's own share | `Share.test_the_training_share_is_one_division_of_two_whole_numbers` |
| C93. The training share worked as the report's share | `Share.test_the_training_share_is_one_division_of_two_whole_numbers` |
| C94. The training share rounded twice, its resistant genomes times one over its genomes | `Share.test_the_training_share_is_one_division_of_two_whole_numbers` |
| C95. A score read as the training share within a tolerance | `Share.test_two_genomes_of_a_fold_at_its_training_share_flag_the_report` |
| C96. A report flagged by 1 genome at its fold's training share | `Share.test_two_genomes_of_a_fold_at_its_training_share_flag_the_report` |
| C97. A report flagged only when a learner gives every genome of a fold its training share | `Share.test_two_genomes_of_a_fold_at_its_training_share_flag_the_report` |
| C98. The genomes at their fold's training share counted over the report, not within one fold | `Share.test_two_genomes_of_a_fold_at_its_training_share_flag_the_report` |
| C99. The share's flag left out of the weighted mean's third count | `Share.test_a_report_the_share_flags_leaves_both_arms_third_counts` |
| C100. The share's flag left out of the logistic arm's third count | `Share.test_a_report_the_share_flags_leaves_both_arms_third_counts` |
| C101. No report flagged by the share | `Share.test_a_report_the_share_flags_leaves_both_arms_third_counts` |
| C102. The flags line giving the share's whole count where the count R30's condition does not flag belongs | `Share.test_the_flags_line_gives_the_share_conditions_counts` |
| C103. Write blind to a reading made without the share condition | `Share.test_write_stops_on_a_reading_made_without_the_share_condition` |
| C104. Write not comparing the sha256 it is given with the reading's bytes | `Sha256.test_write_stops_on_an_edited_reading_given_the_reads_sha256` |
| C105. Write taking a sha256 that is not 64 hex characters | `Sha256.test_write_stops_on_a_value_not_64_hex_characters` |
| C106. Gates.txt written without the reading's sha256 | `Sha256.test_gates_txt_carries_the_readings_sha256` |
| C107. The sha256 an option write may leave out, and write unchecked without it | `Sha256.test_write_needs_the_sha256` |
| C108. Write blind to a fold record with a key missing | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C109. The fold records' stop line naming the report by species only | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C110. The logistic fold records alone checked for a key | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C111. The first fold record of each arm alone checked for a key | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C112. Write naming only the first report whose fold records lack a key | `FoldRecords.test_write_names_each_report_whose_fold_records_lack_a_key` |
| C113. Print not checking the record against SHA256SUMS | `Record.test_print_stops_on_a_file_whose_bytes_give_another_sha256` |
| C114. Print checking rule.tsv alone against SHA256SUMS | `Record.test_print_stops_on_a_file_whose_bytes_give_another_sha256` |
| C115. A record file SHA256SUMS does not list let through | `Record.test_print_stops_on_a_file_the_list_does_not_name` |
| C116. A record with no SHA256SUMS let through | `Record.test_print_stops_on_a_missing_list_or_file` |
| C117. Print's record stop line not naming the file | `Record.test_print_stops_on_a_file_whose_bytes_give_another_sha256` |
| C118. Write checking the rules before the fold records | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C119. The count of coefficients under 0 not among the fold keys write checks | `FoldRecords.test_write_stops_on_a_fold_record_with_a_key_missing` |
| C120. A rule row checked on its outcome alone | `Descriptive.test_write_stops_on_one_field_of_a_logistic_rule_row_changed` |
| C121. A rule row checked on its reading alone | `Descriptive.test_print_stops_on_one_field_of_a_logistic_rule_row_changed` |
| C122. A rule row checked in the count over all alone | `Descriptive.test_write_stops_on_one_field_of_a_logistic_rule_row_changed` |
| C123. Print checking the rule rows of the count over all alone | `Descriptive.test_print_stops_on_one_field_of_a_logistic_rule_row_changed` |

123 of 123 faults failed their test. Each test passes as the script stands and fails with its fault planted.
