# Learner overlap on the lineage split (revision 14)

Revision 14 stores, for each report's lineage split, every learner's held-out score for every genome. This page counts how many of the five learners miss each resistant genome at the release's 0.5 cutoff, and sets the genomes all five miss beside a yardstick: the count expected if the five missed independently. A declared second pass counts what the release's mean of the five misses at 0.5 and works the yardstick again with the genomes of null-fit folds left out. It fits nothing, trains nothing and scores no genome again. The lineage split is a retrospective split inside one database, BV-BRC, and never external validation.

## Reading

The five learners share a floor of misses far above the yardstick, and a larger layer of misses lies beside it, where some of the five miss a genome and the others call it. Summed over the 268 pooled reports, all five miss 6,316 of 93,786 resistant genomes, where the yardstick gives 656.85. Over the 259 reports with a yardstick above 0, the median ratio of the two is 488.8, with quartiles 20.1 to 6,778.9, and 246 of the 259 have a ratio over 1. Another 23,859 of the 93,786 are missed by 1 to 4 of the five, and on a typical report the median share missed by 1 to 4 is 0.263, against 0.057 missed by all five. By the plan's rule, a median of 10 or more, the misses stack: for the genomes all five miss, the data limits recall on this split. The 9 reports with a yardstick of 0, where one learner misses no resistant genome, are counted apart.

The partial layer is room for a combination of the five, and it has a cost on the other class: 22,372 of the 121,533 susceptible genomes are called resistant by 1 to 4 of the five. The release's mean of the five, read as a vote at 0.5, misses 14,383 of the 93,786 resistant genomes: the 6,316 all five miss and 8,067 of the 23,859 that 1 to 4 miss. It calls resistant 8,896 of the 121,533 susceptible genomes: the 3,315 all five call and 5,581 of the 22,372 that 1 to 4 call. On a typical report the mean misses a median share of 0.144 of the resistant genomes.

The null fits do not make the floor. With the genomes of null-fit folds left out, all five miss 5,823 of 88,664 resistant genomes, where the yardstick gives 376.40. The median ratio moves little, from 488.8 over 259 reports to 535.0 over 253, and over the 253 reports with both ratios it goes from 517.0 to 535.0. By the plan's rule the misses still stack.

GRM's own reader of the AMRFinderPlus catalogue calls 3,056 of the 6,316 all-five misses resistant and leaves 3,260 without a resistant call. 719 of those 3,260 are on the 19 reports that hold no catalogue rule for their antibiotic, so those 719 say nothing about the catalogue's coverage. 493 of the 6,316 sit in an outer fold where at least one of the five had a null fit, which scores every genome of the fold at its training prevalence. *Staphylococcus aureus* with clindamycin is the case where the shared misses come mostly from null fits: 229 of its 294 all-five misses sit in such folds, and its ratio is 1.5.

## Overlap

Summed over the 268 pooled reports. A resistant genome is missed when a learner's score is under 0.5, and a susceptible genome is called resistant at 0.5 or more.

| Of the five | Resistant genomes missed | Susceptible genomes called resistant |
|---|---:|---:|
| 0 | 63,611 | 95,846 |
| 1 | 9,795 | 12,151 |
| 2 | 5,720 | 4,755 |
| 3 | 4,447 | 3,249 |
| 4 | 3,897 | 2,217 |
| 5 | 6,316 | 3,315 |
| 1 to 4 | 23,859 | 22,372 |
| total | 93,786 | 121,533 |

## Yardstick

The yardstick of a report is its resistant count times the product of the five learners' miss rates. It is neither a model nor a test, and no p-value is computed. The ratio is read only where the yardstick is above 0, and each report counts once in the median and the quartiles.

| Measure | Value |
|---|---:|
| pooled reports | 268 |
| with a ratio | 259 |
| yardstick of 0, reported apart | 9 |
| missed by all five, summed | 6,316 |
| of those, a catalogue call | 3,056 |
| of those, no catalogue call | 3,260 |
| of those, on a report with no catalogue rule | 719 |
| reports with no catalogue rule | 19 |
| of those, in a fold with a null fit | 493 |
| yardstick, summed | 656.85 |
| median ratio | 488.8 |
| quartiles of the ratio | 20.1 to 6,778.9 |
| reports with a ratio over 1 | 246 |
| median share missed by all five | 0.057 |
| median share missed by 1 to 4 | 0.263 |

## Species

One row per species, over its pooled reports. A species' median weighs each of its reports once.

| Species | Reports | Resistant | Missed by all five | Yardstick | Median ratio | Share by all five | Share by 1 to 4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| *Acinetobacter baumannii* | 18 | 11,818 | 339 | 9.66 | 392.6 | 0.028 | 0.260 |
| *Campylobacter coli* | 8 | 923 | 23 | 0.09 | 105,156.4 | 0.022 | 0.391 |
| *Campylobacter jejuni* | 4 | 605 | 14 | 3.37 | 111,344.1 | 0.027 | 0.435 |
| *Clostridioides difficile* | 4 | 359 | 15 | 0.00 | 6,016.4 | 0.039 | 0.079 |
| *Corynebacterium diphtheriae* | 3 | 104 | 13 | 3.54 | 14.5 | 0.087 | 0.848 |
| *Enterobacter cloacae* | 9 | 444 | 51 | 23.74 | 5.6 | 0.132 | 0.700 |
| *Enterococcus faecium* | 13 | 6,971 | 249 | 40.22 | 2,381.4 | 0.018 | 0.231 |
| *Escherichia coli* | 36 | 11,255 | 791 | 26.11 | 536.3 | 0.071 | 0.211 |
| *Klebsiella aerogenes* | 3 | 100 | 37 | 30.43 | 1.1 | 0.692 | 0.308 |
| *Klebsiella oxytoca* | 4 | 229 | 16 | 0.69 | 143.2 | 0.088 | 0.292 |
| *Klebsiella pneumoniae* | 38 | 15,804 | 1,046 | 115.34 | 137.3 | 0.047 | 0.354 |
| *Klebsiella quasipneumoniae* | 3 | 178 | 10 | 0.00 | 29,552.7 | 0.048 | 0.048 |
| *Mycobacterium tuberculosis* | 22 | 11,683 | 1,239 | 55.54 | 112.7 | 0.121 | 0.256 |
| *Neisseria gonorrhoeae* | 7 | 2,982 | 194 | 36.18 | 101.7 | 0.031 | 0.272 |
| *Pseudomonas aeruginosa* | 13 | 2,832 | 410 | 43.57 | 40.0 | 0.137 | 0.396 |
| *Salmonella enterica* | 24 | 9,584 | 738 | 22.06 | 364.9 | 0.073 | 0.245 |
| *Shigella flexneri* | 4 | 817 | 10 | 0.04 | 16,258.5 | 0.019 | 0.087 |
| *Shigella sonnei* | 20 | 7,890 | 215 | 0.03 | 14,740.6 | 0.031 | 0.093 |
| *Staphylococcus aureus* | 14 | 6,429 | 748 | 235.46 | 1,095.6 | 0.072 | 0.348 |
| *Staphylococcus haemolyticus* | 3 | 154 | 19 | 3.72 | 22.0 | 0.154 | 0.270 |
| *Streptococcus agalactiae* | 3 | 335 | 3 | 0.00 | 29,470.0 | 0.000 | 0.067 |
| *Streptococcus pneumoniae* | 12 | 2,164 | 130 | 7.04 | 3,634.2 | 0.064 | 0.313 |
| *Vibrio cholerae* | 3 | 126 | 6 | 0.00 | 1,283,202.7 | 0.025 | 0.000 |

## Mean vote

The release ships the mean of the five learners' scores. The second pass reads it as a vote: a genome is called resistant when its five scores, added by `math.fsum`, come to 2.5 or more. Summed over the 268 pooled reports:

| Five-learner mean at 0.5 | Genomes |
|---|---:|
| resistant genomes | 93,786 |
| missed by the mean | 14,383 |
| of those, missed by all five | 6,316 |
| of those, missed by 1 to 4 | 8,067 |
| of those, missed by none | 0 |
| median share of a report's resistant genomes the mean misses | 0.144 |
| susceptible genomes | 121,533 |
| called resistant by the mean | 8,896 |
| of those, called resistant by all five | 3,315 |
| of those, called resistant by 1 to 4 | 5,581 |
| of those, called resistant by none | 0 |

A mean of five scores misses every genome all five miss and none that none of the five miss, so the rows by none are 0 and the rows by all five are those of Overlap. One row per species, over its pooled reports, with the median ratio outside null-fit folds from Null-fit cut:

| Species | Reports | Resistant | Missed by the mean | Of those, by all five | Of those, by 1 to 4 | Median ratio outside null-fit folds |
|---|---:|---:|---:|---:|---:|---:|
| *Acinetobacter baumannii* | 18 | 11,818 | 1,223 | 339 | 884 | 424.7 |
| *Campylobacter coli* | 8 | 923 | 40 | 23 | 17 | 105,156.4 |
| *Campylobacter jejuni* | 4 | 605 | 29 | 14 | 15 | 111,344.1 |
| *Clostridioides difficile* | 4 | 359 | 33 | 15 | 18 | 6,016.4 |
| *Corynebacterium diphtheriae* | 3 | 104 | 31 | 13 | 18 | 10.3 |
| *Enterobacter cloacae* | 9 | 444 | 182 | 51 | 131 | 7.2 |
| *Enterococcus faecium* | 13 | 6,971 | 712 | 249 | 463 | 2,381.4 |
| *Escherichia coli* | 36 | 11,255 | 1,624 | 791 | 833 | 536.3 |
| *Klebsiella aerogenes* | 3 | 100 | 61 | 37 | 24 | 1.3 |
| *Klebsiella oxytoca* | 4 | 229 | 65 | 16 | 49 | 143.2 |
| *Klebsiella pneumoniae* | 38 | 15,804 | 2,750 | 1,046 | 1,704 | 217.0 |
| *Klebsiella quasipneumoniae* | 3 | 178 | 14 | 10 | 4 | 29,552.7 |
| *Mycobacterium tuberculosis* | 22 | 11,683 | 2,511 | 1,239 | 1,272 | 72.9 |
| *Neisseria gonorrhoeae* | 7 | 2,982 | 579 | 194 | 385 | 101.7 |
| *Pseudomonas aeruginosa* | 13 | 2,832 | 799 | 410 | 389 | 40.0 |
| *Salmonella enterica* | 24 | 9,584 | 1,711 | 738 | 973 | 555.2 |
| *Shigella flexneri* | 4 | 817 | 41 | 10 | 31 | 16,258.5 |
| *Shigella sonnei* | 20 | 7,890 | 426 | 215 | 211 | 14,740.6 |
| *Staphylococcus aureus* | 14 | 6,429 | 1,237 | 748 | 489 | 1,095.6 |
| *Staphylococcus haemolyticus* | 3 | 154 | 47 | 19 | 28 | 22.0 |
| *Streptococcus agalactiae* | 3 | 335 | 12 | 3 | 9 | 29,470.0 |
| *Streptococcus pneumoniae* | 12 | 2,164 | 249 | 130 | 119 | 517.0 |
| *Vibrio cholerae* | 3 | 126 | 7 | 6 | 1 | 1,283,202.7 |

## Null-fit cut

A null fit leaves a learner with no model in one outer fold, and it scores every genome of that fold at its training prevalence. Where two of the five have null fits in one fold, their misses there come together whatever the genomes hold. The second pass leaves out every genome of a fold where at least one of the five had a null fit, counts each learner's misses again on the genomes left, and works the yardstick and the ratio as the first pass works them. Over the first pass's pooled reports:

| Measure | All genomes | Outside null-fit folds |
|---|---:|---:|
| resistant genomes | 93,786 | 88,664 |
| missed by all five | 6,316 | 5,823 |
| yardstick, summed | 656.85 | 376.40 |
| reports with a ratio | 259 | 253 |
| median ratio | 488.8 | 535.0 |
| quartiles of the ratio | 20.1 to 6,778.9 | 22.2 to 7,918.4 |
| reports with a ratio over 1 | 246 | 240 |
| reports with both ratios | 253 | 253 |
| median ratio over those | 517.0 | 535.0 |

The median ratio moves little, and the plan's rule gives the same reading outside the null-fit folds: the misses stack. The cut leaves out every genome of those folds, not only the shared misses, so each report's resistant count and yardstick change with it.

## Limits

The yardstick assumes the five miss independently, and five learners fitted to the same genomes and the same 31-mers are not expected to. The ratio says how far the shared misses run above that yardstick. It does not say why, and a small yardstick makes a ratio large: the species medians run from 1.1 to 1,283,202.7.

The counts are calls at the release's 0.5 cutoff, the mean vote's too. This page reads no log loss and no ranking measure, so it says nothing of how near the cutoff a missed genome's scores came. Log loss stands beside recall in the release reports, and neither is dropped for the other.

The room for a combination is an upper bound. It counts the resistant genomes at least one of the five calls, which a rule could reach only if it knew which learner to trust on each genome. The mean vote is one fixed combination, read from stored scores. No combination was fitted, and the page does not say any model would fail or succeed on the genomes all five miss.

A null fit scores all its fold's genomes alike, so where two learners have null fits in one fold, their misses there come together. The release gives null fits to the boosted trees and the lasso only (docs/equal-scoring.md), and 493 of the 6,316 all-five misses sit in such folds. Left out with every other genome of those folds, they move the median ratio from 488.8 to 535.0 (Null-fit cut).

The labels are the release's laboratory labels as BV-BRC holds them. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient. A row per report is an aggregate and keeps its genome count. The 6 reports with fewer than 2 learners on the lineage split were skipped and counted; they are named under Datasets.

## Method

The plan below was committed before the table script. The script and its tests were committed and proved on a box, and `scripts/mechanism/learner_overlap.py` then read the 274 reports and wrote `reference/learner-overlap/`. The commands, as `reference/learner-overlap/reading.txt` gives them:

```
python3 -I scripts/mechanism/learner_overlap.py plan
python3 -I scripts/mechanism/learner_overlap.py read REPORTS --out READING
python3 -I scripts/mechanism/learner_overlap.py write READING --out RECORD
python3 -I scripts/mechanism/learner_overlap.py print RECORD
```

The first read stopped at check 3 and wrote nothing, as the plan's section "Change after the first read stopped" sets down. Run again at the script's fix, the read passed every check:

```
check 1 parse: 274 files found, 274 read, 274 expected: pass
check 2 lengths: 0 reports with a score list or test_fold of another length or a bad score: pass
check 3 counts: 0 reports whose counts or labels differ from their samples; 19 used reports with no catalogue rule (catalog_resistant null on every genome): pass
check 4 learner sets: 1 sets; the release's is cart,gbt,lasso,rf,scm on 268 reports; 0 reports differ
check 5 skipped: 0 with no lineage split, 6 with fewer than 2 learners; 6 skipped and 268 used of 274 found: pass
check 6 ids: 5 files scanned, 0 id-shaped tokens or barred keys: pass
check 7 stored metrics: 0 learner rows whose counts at the cutoff differ from metrics: pass
check 8 pairs: 0 species and antibiotic pairs on two reports: pass
check 9 catalogue: 0 reports whose catalog calls differ from catalog_resistant: pass
check 10 null folds: 0 learner rows whose null-fold genomes differ from no_call
reading rule: stack
```

The plan's section "Second pass" was fixed on 25 November 2026, after a review of the first reading, because the first pass counted the learners and not the vote the release ships. Its script and tests were committed and proved on a box, and `scripts/mechanism/overlap_second.py` then read the same 274 reports once against the committed record and wrote `reference/learner-overlap/second-pass/`. The commands, as `reference/learner-overlap/second-pass/reading.txt` gives them:

```
python3 -I scripts/mechanism/overlap_second.py plan
python3 -I scripts/mechanism/overlap_second.py read REPORTS --first FIRST --out READING
python3 -I scripts/mechanism/overlap_second.py write READING --first FIRST --out RECORD
python3 -I scripts/mechanism/overlap_second.py print RECORD --first FIRST
```

Every first-pass check, run again, printed the line it printed above; check 6 runs at the write and stands here as second check 6. The second pass's own checks:

```
second check 1 first pass: 274 files found, every first-pass check passed, 268 reports used: pass
second check 2 first-pass columns: 268 reports read again, 268 rows in datasets.tsv, 0 that differ: pass
second check 3 arithmetic of the mean: 268 reports with the five, 0 where the mean misses a genome none of the five miss, calls one none call, or parts from all five: pass
second check 4 null-fit cut: 0 reports where the all-five misses outside and inside null-fit folds do not add up: pass
second check 5 the vote's edge: 0 genomes where a plain sum over 5 gives another call than math.fsum; 0 genomes whose sum lies within 1e-9 of 2.5
second check 6 ids: 5 files scanned, 0 id-shaped tokens or barred keys: pass
reading rule outside null-fit folds: stack
```

Check 6 and second check 6 scan for the shape of a genome id and for a key named `id` or `samples`. Beside them, every file written for this page was compared on the box with the reports' own genome ids, held in memory and never written, and none matched. The shape scan stays as a second net only: in a control file of 1,480 genome ids, 1,417 had the shape.

`tests/mechanism/test_learner_overlap.py` holds 38 SYNTHETIC tests of the table script, with counts worked by hand. `tests/mechanism/test_overlap_second.py` holds 22 SYNTHETIC tests of the second pass's script, worked the same way. `tests/mechanism/test_learner_overlap_page.py` holds 20 tests of this page and its two records. `scripts/mechanism/overlap_faults.py` plants each fault below in one of the two scripts in turn and runs the test that must catch it:

| Fault | Test |
|---|---|
| O1. A score of exactly 0.5 read as a miss | `Counts.test_a_score_of_exactly_the_cutoff_is_a_call` |
| O2. The random split read in place of the lineage split | `Counts.test_only_the_lineage_split_is_read` |
| O3. The catalogue reader counted among the learners | `Counts.test_the_catalogue_reader_is_not_a_learner` |
| O4. A score list of another length let through | `Stops.test_a_score_list_of_another_length_stops_the_read` |
| O5. A report with one learner kept | `Skips.test_no_lineage_split_and_one_learner_are_skipped_and_counted` |
| O6. The yardstick from the mean miss rate | `Counts.test_all_five_miss_the_same_genomes` |
| O7. The median weighted by genome count | `Pooled.test_the_median_weights_each_report_once` |
| O8. A report with a yardstick of 0 dropped from the pool | `Pooled.test_a_report_with_a_yardstick_of_0_is_reported_apart` |
| O9. An all-five miss the catalogue reader calls counted as uncalled | `Counts.test_misses_with_no_catalogue_call` |
| O10. A genome id kept as a report is parsed | `Ids.test_no_id_is_held_after_parsing` |
| O11. Check 6 blind to an id-shaped token | `Ids.test_an_id_shaped_token_in_the_reading_stops_the_write` |
| O12. Counts that differ from the samples let through | `Stops.test_counts_that_differ_from_the_samples_stop_the_read` |
| O13. A report under 10 resistant genomes pooled | `Pooled.test_a_report_under_10_resistant_stays_out_of_the_pool` |
| O14. Stored metrics that differ from the recount let through | `Stops.test_stored_metrics_that_differ_from_the_recount_stop_the_read` |
| O15. Catalogue calls that differ from catalog_resistant let through | `Stops.test_catalogue_calls_that_differ_from_catalog_resistant_stop_the_read` |
| O16. Every fold of a learner with a null fit read as null | `Counts.test_misses_in_a_null_fit_fold` |
| O17. A null catalogue call read as a call | `Counts.test_a_report_with_no_catalogue_rule_is_counted_and_marked` |
| O18. A score sum of exactly 2.5 read as a miss by the mean | `Vote.test_a_sum_of_exactly_2_5_is_a_call` |
| O19. The majority of the five read in place of their mean | `Vote.test_the_mean_is_not_the_majority` |
| O20. A plain sum over 5 read in place of math.fsum | `Vote.test_the_vote_edge_is_counted_and_the_read_goes_on` |
| O21. The null folds of one learner cut in place of all five's | `NullCut.test_null_folds_of_two_learners_are_both_cut` |
| O22. First-pass columns that differ from the record let through | `SecondStops.test_a_first_record_that_differs_stops_the_read` |
| O23. A vote that parts from the five let through | `Vote.test_a_vote_that_parts_from_the_five_stops_the_read` |
| O24. A null-fit cut that does not add up let through | `NullCut.test_a_cut_that_does_not_add_up_stops_the_read` |
| O25. The median outside null folds weighted by genome count | `SecondPooled.test_the_median_outside_null_folds_weights_each_report_once` |
| O26. A report under 10 resistant genomes pooled in the second pass | `SecondPooled.test_a_report_under_10_resistant_stays_out_of_the_pool` |
| O27. The second pass's shape scan blind to an id-shaped token | `SecondStops.test_an_id_shaped_token_in_the_reading_stops_the_write` |
| O28. A genome 3 of the five miss counted with those all five miss | `Vote.test_the_mean_misses_what_three_or_more_of_the_five_miss_here` |

All 28 were caught.

## Plan

This plan was fixed on 24 November 2026, before the table script was written and before any report was read for this question. It fits nothing, trains nothing and scores no genome again. Every figure comes from the held-out scores and labels that the 274 revision 14 release reports already store.

### Question

Each revision 14 report stores, for its lineage split, every learner's held-out score for every genome. When the learners miss a resistant genome, do they miss the same one? If the five learners miss largely the same resistant genomes, a sixth learner of any shape has little room to lift recall on them, and for those genomes the data limits recall. If they miss different genomes, a combination of the five has room. The lineage split is a retrospective split inside one database, BV-BRC. It is not external validation.

### Data

- The 274 report files of the revision 14 release's `datasets/` directory, one per species and antibiotic.
- Per report: `samples[].resistant` (the label), `samples[].catalog_resistant` (GRM's own reader of the AMRFinderPlus catalogue, the known-gene rule) and, in the split whose `scheme` is `lineage`, each learner entry's `learner`, `scores`, `metrics`, `no_call` and `fold_complexity`, and the split's `test_fold`. `counts.samples` and `counts.resistant` are checked against `samples`.
- The random split is not read. A `samples` entry also holds a genome id. The table script drops every key named `id` as a file is parsed, so no genome id is held, printed or written.
- The cutoff is the release's own. A score of 0.5 or more is a resistant call and a score under 0.5 is not. No other cutoff is tried.

### Learners

The learners are the lineage split's entries other than `catalog`, which is GRM's reader of the AMRFinderPlus catalogue and not a learner. A first look at one report's keys found `scm`, `cart`, `rf`, `gbt`, `lasso` and `catalog`. No combined or mean-vote entry stands in the split. The read prints the learner names it finds on every report. The five are the name set found on the most reports. A report with another set is listed by species and antibiotic, keeps its row and stays out of the pooled rows. `catalog` enters only through `catalog_resistant` and through check 9.

### Table `reference/learner-overlap/datasets.tsv`

One row per report used, an aggregate that keeps its genome count. Columns, in order:

- `species`, `antibiotic`: as the report names them, then `genomes`, `resistant` and `susceptible`;
- `learners`: the learner names used, comma-joined, and `learner_count`;
- `miss_0` to `miss_5`: resistant genomes scored under 0.5 by exactly k of the learners, which add to `resistant`;
- `wrong_0` to `wrong_5`: susceptible genomes scored 0.5 or more by exactly k of the learners, which add to `susceptible`;
- `miss_scm`, `miss_cart`, `miss_rf`, `miss_gbt`, `miss_lasso`: each learner's resistant genomes under 0.5;
- `miss_all5`: resistant genomes missed by all five (`miss_5`);
- `miss_all5_no_catalog_call`: of those, the ones `catalog_resistant` leaves without a resistant call;
- `miss_all5_in_null_fold`: of those, the ones in an outer fold where at least one of the five had a null fit (`fold_complexity` 0). Such a fold scores every genome at its training prevalence (docs/equal-scoring.md);
- `null_fold_pairs`: learner and fold pairs with a null fit;
- `expected_all5`: the yardstick, `resistant` times the product of the five miss rates `miss_<learner> / resistant`. It is the count the five would miss together if they missed independently. It is a yardstick. It is neither a model nor a test;
- `ratio_all5`: `miss_all5 / expected_all5`, blank when `expected_all5` is 0.

For a report with fewer than five learners, the `miss_k` and `wrong_k` columns above its learner count, `miss_all5` and the columns after it are blank. A report with more than five learners keeps its counts and names, and its `miss_k` and `wrong_k` columns are blank too (added with the table script, before the read). Every decimal is written with 12 places. `expected_all5` and `ratio_all5` are worked as exact fractions and rounded once.

### Table `reference/learner-overlap/pooled.tsv`

The same counts summed over the pooled reports: those used, carrying the five, with 10 or more resistant genomes. One row with `scope` `all`, then one row per species with `scope` `species`. Beside the sums:

- `datasets`, `datasets_with_ratio` (`expected_all5` over 0) and `datasets_expected_zero`, the reports where one learner misses no resistant genome, reported apart and never dropped silently;
- `ratio_median`, `ratio_q1`, `ratio_q3`: the median and quartiles of `ratio_all5` over the reports with a ratio, each report weighted once, none by its genome count (quartiles by Python's `statistics.quantiles` with `n=4` and `method="inclusive"`; with one report the three are its ratio);
- `datasets_ratio_over_1`: the reports whose ratio is over 1;
- `share_all5_median`: the median over the pooled reports of `miss_all5 / resistant`;
- `share_some_median`: the median over the pooled reports of `(miss_1 + miss_2 + miss_3 + miss_4) / resistant`, the resistant genomes at least one of the five calls and at least one misses.

### Checks

Each prints counts only. A check marked stop ends the read with no table written.

1. Every report file parses, and the files read equal the files found (stop). The files found must be 274 (stop).
2. For every learner entry of every report used, `len(scores) == len(samples)`, each score is a number from 0 to 1, and `len(test_fold) == len(samples)` (stop).
3. `counts.resistant` equals the labels counted in `samples`, and `counts.samples` equals its length (stop).
4. The learner name set is the same on every report used, or the reports that differ are listed by species and antibiotic and counted.
5. Reports with no lineage split, or with fewer than 2 learners, are skipped and counted, and the skipped plus the used equal the files found (stop). A report with two lineage splits stops the read.
6. No genome id stands in any written file: every file the read and the write make is scanned for a token of 2 to 7 digits, a point and 1 to 9 digits, standing alone (the shape of a BV-BRC genome id, a taxon number and a genome number), and for a key named `id` or `samples`. It must find 0 (stop; the files stay unwritten).
7. Each learner's counts at 0.5, recounted from its scores, equal its stored `metrics` `tp`, `fn`, `fp` and `tn` (stop). This proves the scores line up with `samples` and the cutoff reads as the release reads it.
8. No species and antibiotic pair stands on two reports (stop).
9. Where a `catalog` entry stands, its calls at 0.5 equal `catalog_resistant` on every genome (stop).
10. Each learner's genomes in its null-fit folds, by class, are compared with its `no_call`. The rows that differ are counted and the read goes on.

### Reading rule

The statistic is the median over the pooled reports of `ratio_all5`. The reading rule, fixed here:

- a median of 10 or more: the misses stack. The five learners miss largely the same resistant genomes, and for those genomes the data limits recall.
- a median over 1 and under 10: the five share misses more often than the yardstick, by under tenfold.
- a median of 1 or less: no stacking above the yardstick is shown.

Beside it, `share_some_median` against `share_all5_median` says whether more of a typical report's resistant genomes sit where some of the five call them, where a combination has room, or where all five miss. No p-value is computed and nothing here is a test of a null hypothesis. The page will not say that any model would fail, and it names no gene, variant or position.

### Change after the first read stopped

The first read, on 24 November 2026 at the script's first commit, stopped at check 3 and wrote nothing. Of the 274 reports, 21 hold no catalogue rule for their antibiotic: `catalog_resistant` is null on every genome and the lineage split has no `catalog` entry. Check 3 had required a true or false call, so 19 of them stopped it; the other 2 hold no learner and were already skipped. The plan read GRM's reader of the AMRFinderPlus catalogue as making a call on every genome, and on these reports it makes none. The table script was changed before the second read, which is the one this page reads:

- a null `catalog_resistant` is no resistant call, so `miss_all5_no_catalog_call` counts it;
- `datasets.tsv` gains `catalog_rule` after `learner_count`: `yes` where every call is true or false, `no` where every call is null;
- `pooled.tsv` gains `datasets_no_catalog_rule` and `miss_all5_no_catalog_rule`, the all-five misses on reports with no rule, before `expected_all5`;
- a label that is not true or false still stops the read, and so does a report that mixes null calls with true or false ones.

No other rule, column, check or reading rule changed.

### Second pass

Fixed on 25 November 2026, after a review of the first reading, before the second pass's script was written and before the reports were read again. The first pass counted the five learners one by one, and the release ships their mean, so the second pass counts what the mean misses. It also works the yardstick again with the genomes of null-fit folds left out, since a null fit is one learner with no model, not five learners agreeing. It fits nothing, and no first-pass rule, column, file or figure changes. The mean of five stored scores is arithmetic, not a fit.

- The mean vote calls a genome resistant when the sum of its five learner scores, added by Python's `math.fsum`, is 2.5 or more: the mean at the release's 0.5 cutoff. It is read only on reports that carry the five.
- `scripts/mechanism/overlap_second.py` reads the same 274 reports. It runs every first-pass check again through the first pass's own reader, then reads each used report's five score lists again.

Table `reference/learner-overlap/second-pass/vote.tsv`, one row per report the first pass used, species and antibiotic in full, with `genomes`, `resistant` and `susceptible`:

- `mean_miss`: resistant genomes the mean misses; `mean_miss_by_all5`, `mean_miss_by_1to4` and `mean_miss_by_none`: of those, the ones missed by all five learners, by 1 to 4 and by none;
- `mean_wrong`: susceptible genomes the mean calls resistant; `mean_wrong_by_all5`, `mean_wrong_by_1to4` and `mean_wrong_by_none`: of those, the ones all five, 1 to 4 or none of the five call resistant;
- `nonull_resistant`: resistant genomes outside every outer fold where at least one of the five had a null fit; `nonull_miss_scm`, `nonull_miss_cart`, `nonull_miss_rf`, `nonull_miss_gbt`, `nonull_miss_lasso` and `nonull_miss_all5`: their misses at 0.5;
- `nonull_expected_all5` and `nonull_ratio_all5`: the yardstick and the ratio worked from those columns as the first pass works them, the ratio blank where the yardstick is 0.

A report with another learner set keeps its row with these columns blank.

Table `reference/learner-overlap/second-pass/vote-pooled.tsv`: over the first pass's pooled reports, one row with `scope` `all` and one per species. It holds:

- the sums of the count columns and of `nonull_expected_all5`;
- `mean_share_median`: the median of `mean_miss / resistant`;
- `nonull_datasets_with_ratio`, `nonull_ratio_median`, `nonull_ratio_q1`, `nonull_ratio_q3` and `nonull_datasets_ratio_over_1`: as the first pass works them, each report weighted once;
- `paired_datasets`: the reports with both a first-pass ratio and a ratio outside null folds, with `paired_ratio_median` and `paired_nonull_ratio_median`, the two medians over those same reports.

Checks of the second pass, each a stop unless marked:

1. Every check of the first pass runs again and passes, and the files found are 274.
2. The reports the second pass uses are the rows of the committed `datasets.tsv`, and each report's first-pass columns, counted again, equal its row.
3. A mean of five scores cannot miss a genome none of the five miss, nor call resistant one none of the five call, and it misses every genome all five miss. So `mean_miss_by_none` and `mean_wrong_by_none` must be 0, and `mean_miss_by_all5` and `mean_wrong_by_all5` must equal `miss_5` and `wrong_5`.
4. `nonull_miss_all5` plus `miss_all5_in_null_fold` equals `miss_all5`.
5. The genomes where a plain sum, added in learner order and divided by 5, gives another call than `math.fsum` are counted, and so are the genomes whose sum lies within 1e-9 of 2.5. The read goes on.
6. Every written file is scanned for the shape of a genome id and for a key named `id` or `samples` (the first pass's check 6).

Reading: the first pass's reading rule is applied again to `nonull_ratio_median`. The page gives it beside the first pass's median and says whether the rule's reading changes. For the mean vote no rule is set: the page reports its counts.

## Datasets

One row per report used, species and antibiotic in full. A report with fewer than five learners would show its counts and names only; none was found. "No rule" marks a report with no catalogue rule for its antibiotic, and "none" a ratio where the yardstick is 0. A yardstick shown as 0.00 is above 0 and rounds to 0.00 at two places, so its ratio is read.

| Species | Antibiotic | Genomes | Resistant | Missed by all five | No catalogue call | In a null-fit fold | Yardstick | Ratio |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| *Acinetobacter baumannii* | amikacin | 1,277 | 785 | 27 | 15 | 0 | 0.07 | 360.6 |
| *Acinetobacter baumannii* | ampicillin/sulbactam | 989 | 555 | 62 | 0 | 0 | 0.15 | 424.7 |
| *Acinetobacter baumannii* | cefepime | 437 | 412 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Acinetobacter baumannii* | cefotaxime | 808 | 778 | 10 | 0 | 0 | 0.00 | 66,288.4 |
| *Acinetobacter baumannii* | ceftazidime | 1,247 | 1,087 | 33 | 0 | 0 | 0.01 | 5,413.8 |
| *Acinetobacter baumannii* | ceftriaxone | 1,070 | 1,011 | 7 | 0 | 0 | 0.00 | 6,450.2 |
| *Acinetobacter baumannii* | ciprofloxacin | 1,440 | 1,319 | 9 | 9 | 0 | 0.00 | 216,197,738.4 |
| *Acinetobacter baumannii* | doripenem | 213 | 61 | 6 | 0 | 0 | 0.86 | 7.0 |
| *Acinetobacter baumannii* | doxycycline | 167 | 107 | 9 | 1 | 0 | 0.10 | 88.6 |
| *Acinetobacter baumannii* | gentamicin | 1,342 | 1,204 | 4 | 3 | 0 | 0.00 | 7,918.4 |
| *Acinetobacter baumannii* | imipenem | 1,544 | 794 | 90 | 0 | 0 | 0.41 | 220.0 |
| *Acinetobacter baumannii* | levofloxacin | 919 | 760 | 20 | 20 | 0 | 0.00 | 82,950.9 |
| *Acinetobacter baumannii* | meropenem | 666 | 455 | 21 | 0 | 0 | 0.02 | 931.6 |
| *Acinetobacter baumannii* | minocycline | 139 | 28 | 2 | 2 | 0 | 3.49 | 0.6 |
| *Acinetobacter baumannii* | moxifloxacin | 91 | 60 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Acinetobacter baumannii* | sulfamethoxazole/trimethoprim | 1,114 | 962 | 4 | 4 | 0 | 0.01 | 592.2 |
| *Acinetobacter baumannii* | tetracycline | 962 | 759 | 15 | 5 | 0 | 0.13 | 114.9 |
| *Acinetobacter baumannii* | tobramycin | 981 | 681 | 20 | 15 | 0 | 4.42 | 4.5 |
| *Campylobacter coli* | azithromycin | 389 | 78 | 2 | 1 | 0 | 0.00 | 60,482.1 |
| *Campylobacter coli* | ciprofloxacin | 391 | 132 | 0 | 0 | 0 | 0.00 | none |
| *Campylobacter coli* | clindamycin | 389 | 94 | 11 | 10 | 0 | 0.04 | 281.4 |
| *Campylobacter coli* | erythromycin | 391 | 78 | 1 | 0 | 0 | 0.00 | 105,156.4 |
| *Campylobacter coli* | gentamicin | 391 | 65 | 4 | 0 | 0 | 0.05 | 81.2 |
| *Campylobacter coli* | nalidixic acid | 375 | 119 | 3 | 3 | 0 | 0.00 | 225,066.1 |
| *Campylobacter coli* | telithromycin | 374 | 53 | 1 | no rule | 0 | 0.00 | 131,508.0 |
| *Campylobacter coli* | tetracycline | 391 | 304 | 1 | 1 | 0 | 0.00 | 8,540,717,056.0 |
| *Campylobacter jejuni* | ciprofloxacin | 508 | 113 | 3 | 3 | 0 | 0.00 | 118,407.7 |
| *Campylobacter jejuni* | clindamycin | 508 | 25 | 6 | 5 | 3 | 3.37 | 1.8 |
| *Campylobacter jejuni* | nalidixic acid | 506 | 110 | 3 | 3 | 0 | 0.00 | 104,280.6 |
| *Campylobacter jejuni* | tetracycline | 508 | 357 | 2 | 2 | 0 | 0.00 | 203,040,595.0 |
| *Clostridioides difficile* | azithromycin | 203 | 91 | 7 | 1 | 0 | 0.00 | 3,767.5 |
| *Clostridioides difficile* | ceftriaxone | 136 | 90 | 1 | 0 | 0 | 0.00 | 8,265.3 |
| *Clostridioides difficile* | clarithromycin | 203 | 92 | 5 | 2 | 0 | 0.00 | 1,762.1 |
| *Clostridioides difficile* | moxifloxacin | 204 | 86 | 2 | 2 | 0 | 0.00 | 854,700.2 |
| *Corynebacterium diphtheriae* | sulfamethoxazole/trimethoprim | 236 | 28 | 1 | 1 | 0 | 0.00 | 1,045.3 |
| *Corynebacterium diphtheriae* | tetracycline | 243 | 46 | 4 | 4 | 1 | 0.28 | 14.5 |
| *Corynebacterium diphtheriae* | trimethoprim | 229 | 30 | 8 | 8 | 3 | 3.27 | 2.5 |
| *Enterobacter cloacae* | cefepime | 103 | 70 | 4 | 0 | 0 | 0.00 | 1,368.1 |
| *Enterobacter cloacae* | ciprofloxacin | 110 | 67 | 6 | 3 | 0 | 0.74 | 8.2 |
| *Enterobacter cloacae* | gentamicin | 104 | 27 | 5 | 0 | 1 | 1.60 | 3.1 |
| *Enterobacter cloacae* | imipenem | 86 | 58 | 8 | 0 | 0 | 4.59 | 1.7 |
| *Enterobacter cloacae* | levofloxacin | 104 | 58 | 8 | 3 | 3 | 5.52 | 1.4 |
| *Enterobacter cloacae* | meropenem | 99 | 40 | 12 | 3 | 1 | 11.19 | 1.1 |
| *Enterobacter cloacae* | sulfamethoxazole/trimethoprim | 80 | 41 | 0 | 0 | 0 | 0.00 | none |
| *Enterobacter cloacae* | tobramycin | 99 | 45 | 3 | 1 | 0 | 0.07 | 40.2 |
| *Enterobacter cloacae* | trimethoprim | 124 | 38 | 5 | 5 | 0 | 0.03 | 199.5 |
| *Enterococcus faecium* | ampicillin | 1,680 | 1,291 | 9 | 9 | 0 | 0.00 | 46,469,028.0 |
| *Enterococcus faecium* | clindamycin | 901 | 619 | 94 | 1 | 0 | 0.48 | 196.9 |
| *Enterococcus faecium* | dalfopristin/quinupristin | 1,118 | 112 | 25 | 1 | 7 | 15.10 | 1.7 |
| *Enterococcus faecium* | daptomycin | 333 | 55 | 35 | 35 | 8 | 24.45 | 1.4 |
| *Enterococcus faecium* | doxycycline | 345 | 136 | 8 | 3 | 0 | 0.00 | 2,011.2 |
| *Enterococcus faecium* | erythromycin | 1,151 | 1,063 | 18 | 0 | 0 | 0.02 | 765.7 |
| *Enterococcus faecium* | gentamicin | 1,345 | 396 | 21 | 0 | 0 | 0.01 | 2,381.4 |
| *Enterococcus faecium* | kanamycin | 944 | 612 | 11 | 0 | 0 | 0.00 | 2,901.1 |
| *Enterococcus faecium* | linezolid | 1,367 | 70 | 3 | 2 | 0 | 0.15 | 19.7 |
| *Enterococcus faecium* | streptomycin | 1,239 | 488 | 6 | 4 | 0 | 0.00 | 27,626.7 |
| *Enterococcus faecium* | teicoplanin | 1,192 | 516 | 5 | 5 | 0 | 0.00 | 4,255,237.5 |
| *Enterococcus faecium* | tetracycline | 958 | 717 | 13 | 7 | 0 | 0.00 | 221,689.1 |
| *Enterococcus faecium* | vancomycin | 1,741 | 896 | 1 | 0 | 0 | 0.00 | 825,825.1 |
| *Escherichia coli* | amikacin | 1,125 | 30 | 8 | 6 | 8 | 2.37 | 3.4 |
| *Escherichia coli* | amoxicillin/clavulanate | 1,115 | 548 | 73 | 1 | 0 | 1.42 | 51.4 |
| *Escherichia coli* | ampicillin | 1,965 | 1,596 | 47 | 2 | 0 | 0.00 | 45,886.4 |
| *Escherichia coli* | ampicillin/sulbactam | 283 | 142 | 28 | 1 | 0 | 14.17 | 2.0 |
| *Escherichia coli* | aztreonam | 699 | 225 | 15 | 0 | 0 | 0.20 | 74.7 |
| *Escherichia coli* | cefazolin | 414 | 322 | 5 | 0 | 0 | 0.05 | 108.2 |
| *Escherichia coli* | cefepime | 1,124 | 211 | 23 | 1 | 0 | 0.82 | 28.2 |
| *Escherichia coli* | cefotaxime | 1,149 | 615 | 69 | 4 | 0 | 0.06 | 1,153.8 |
| *Escherichia coli* | cefoxitin | 1,066 | 399 | 41 | 1 | 0 | 0.17 | 243.2 |
| *Escherichia coli* | cefpodoxime | 175 | 103 | 21 | 0 | 0 | 0.30 | 69.2 |
| *Escherichia coli* | cefpodoxime/clavulanate | 136 | 101 | 2 | 0 | 0 | 0.01 | 263.8 |
| *Escherichia coli* | ceftazidime | 1,466 | 601 | 33 | 3 | 0 | 0.04 | 866.7 |
| *Escherichia coli* | ceftiofur | 203 | 105 | 1 | 0 | 0 | 0.00 | 33,764.1 |
| *Escherichia coli* | ceftriaxone | 923 | 438 | 31 | 0 | 0 | 0.02 | 1,678.8 |
| *Escherichia coli* | cefuroxime | 548 | 272 | 37 | 4 | 0 | 0.12 | 305.0 |
| *Escherichia coli* | cephalothin | 472 | 278 | 20 | 3 | 0 | 0.01 | 1,866.7 |
| *Escherichia coli* | chloramphenicol | 709 | 250 | 19 | 10 | 0 | 0.03 | 577.4 |
| *Escherichia coli* | ciprofloxacin | 1,615 | 826 | 30 | 19 | 0 | 0.00 | 60,713.7 |
| *Escherichia coli* | colistin | 239 | 50 | 0 | 0 | 0 | 0.00 | none |
| *Escherichia coli* | doxycycline | 227 | 147 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Escherichia coli* | ertapenem | 820 | 205 | 17 | 17 | 0 | 0.08 | 223.9 |
| *Escherichia coli* | erythromycin | 205 | 71 | 0 | 0 | 0 | 0.00 | none |
| *Escherichia coli* | gentamicin | 1,695 | 455 | 48 | 20 | 0 | 0.10 | 488.8 |
| *Escherichia coli* | imipenem | 614 | 113 | 12 | 11 | 0 | 0.02 | 495.2 |
| *Escherichia coli* | levofloxacin | 403 | 294 | 7 | 4 | 0 | 0.00 | 23,644.7 |
| *Escherichia coli* | meropenem | 1,539 | 82 | 22 | 15 | 6 | 3.02 | 7.3 |
| *Escherichia coli* | nalidixic acid | 406 | 164 | 3 | 3 | 0 | 0.00 | 36,452.2 |
| *Escherichia coli* | nitrofurantoin | 482 | 46 | 14 | 14 | 1 | 2.77 | 5.0 |
| *Escherichia coli* | norfloxacin | 197 | 41 | 4 | 1 | 0 | 0.00 | 1,698.2 |
| *Escherichia coli* | piperacillin/tazobactam | 979 | 292 | 46 | 2 | 0 | 0.28 | 166.8 |
| *Escherichia coli* | streptomycin | 386 | 254 | 10 | 8 | 0 | 0.01 | 984.5 |
| *Escherichia coli* | sulfamethoxazole | 341 | 161 | 10 | 7 | 0 | 0.00 | 3,165.7 |
| *Escherichia coli* | sulfamethoxazole/trimethoprim | 925 | 580 | 22 | 19 | 0 | 0.01 | 2,626.3 |
| *Escherichia coli* | tetracycline | 782 | 464 | 26 | 22 | 0 | 0.00 | 13,457.6 |
| *Escherichia coli* | tobramycin | 1,070 | 343 | 27 | 19 | 0 | 0.01 | 2,046.3 |
| *Escherichia coli* | trimethoprim | 861 | 431 | 20 | 18 | 0 | 0.01 | 1,588.3 |
| *Klebsiella aerogenes* | ertapenem | 83 | 48 | 0 | 0 | 0 | 0.18 | 0.0 |
| *Klebsiella aerogenes* | imipenem | 68 | 26 | 19 | 17 | 19 | 13.23 | 1.4 |
| *Klebsiella aerogenes* | meropenem | 78 | 26 | 18 | 16 | 16 | 17.01 | 1.1 |
| *Klebsiella oxytoca* | cefotaxime | 98 | 57 | 4 | 0 | 0 | 0.03 | 126.0 |
| *Klebsiella oxytoca* | cefoxitin | 99 | 47 | 6 | 0 | 0 | 0.04 | 160.4 |
| *Klebsiella oxytoca* | cefpodoxime | 74 | 47 | 5 | 0 | 0 | 0.01 | 550.1 |
| *Klebsiella oxytoca* | ertapenem | 112 | 78 | 1 | 0 | 0 | 0.61 | 1.6 |
| *Klebsiella pneumoniae* | amikacin | 1,588 | 318 | 23 | 14 | 0 | 4.80 | 4.8 |
| *Klebsiella pneumoniae* | amoxicillin/clavulanate | 778 | 591 | 2 | 0 | 0 | 0.01 | 273.9 |
| *Klebsiella pneumoniae* | ampicillin/sulbactam | 744 | 536 | 13 | 6 | 0 | 0.01 | 893.6 |
| *Klebsiella pneumoniae* | aztreonam | 1,050 | 693 | 4 | 3 | 0 | 0.04 | 95.7 |
| *Klebsiella pneumoniae* | cefazolin | 1,141 | 763 | 21 | 17 | 0 | 0.02 | 1,214.9 |
| *Klebsiella pneumoniae* | cefepime | 1,233 | 619 | 38 | 3 | 2 | 0.46 | 81.8 |
| *Klebsiella pneumoniae* | cefotetan | 61 | 25 | 19 | 0 | 19 | 14.18 | 1.3 |
| *Klebsiella pneumoniae* | cefoxitin | 1,181 | 550 | 101 | 19 | 0 | 3.75 | 26.9 |
| *Klebsiella pneumoniae* | cefpodoxime | 221 | 180 | 5 | 0 | 0 | 0.00 | 3,029.1 |
| *Klebsiella pneumoniae* | cefpodoxime/ceftriaxone | 56 | 28 | 1 | no rule | 0 | 0.00 | 1,600.7 |
| *Klebsiella pneumoniae* | cefpodoxime/clavulanate | 207 | 182 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Klebsiella pneumoniae* | ceftazidime | 1,310 | 860 | 29 | 12 | 0 | 0.05 | 546.5 |
| *Klebsiella pneumoniae* | ceftazidime/avibactam | 259 | 79 | 3 | 0 | 0 | 0.59 | 5.1 |
| *Klebsiella pneumoniae* | ceftriaxone | 1,208 | 843 | 4 | 3 | 0 | 0.00 | 14,402.8 |
| *Klebsiella pneumoniae* | cefuroxime | 565 | 468 | 6 | 5 | 0 | 0.00 | 3,756.1 |
| *Klebsiella pneumoniae* | chloramphenicol | 409 | 302 | 15 | 0 | 0 | 0.11 | 137.3 |
| *Klebsiella pneumoniae* | ciprofloxacin | 1,706 | 1,096 | 58 | 0 | 0 | 0.01 | 4,891.2 |
| *Klebsiella pneumoniae* | colistin | 336 | 112 | 14 | 14 | 3 | 6.60 | 2.1 |
| *Klebsiella pneumoniae* | doripenem | 243 | 217 | 0 | 0 | 0 | 0.00 | none |
| *Klebsiella pneumoniae* | doxycycline | 141 | 46 | 18 | 4 | 9 | 5.32 | 3.4 |
| *Klebsiella pneumoniae* | ertapenem | 952 | 671 | 11 | 9 | 0 | 0.02 | 469.4 |
| *Klebsiella pneumoniae* | fosfomycin | 294 | 83 | 27 | 0 | 0 | 11.34 | 2.4 |
| *Klebsiella pneumoniae* | gentamicin | 1,842 | 689 | 92 | 27 | 0 | 0.18 | 498.1 |
| *Klebsiella pneumoniae* | imipenem | 1,158 | 615 | 24 | 19 | 0 | 0.03 | 806.3 |
| *Klebsiella pneumoniae* | levofloxacin | 871 | 549 | 21 | 0 | 0 | 0.00 | 40,105.0 |
| *Klebsiella pneumoniae* | meropenem | 1,765 | 736 | 74 | 25 | 32 | 0.49 | 149.6 |
| *Klebsiella pneumoniae* | nitrofurantoin | 562 | 375 | 48 | 44 | 0 | 4.15 | 11.6 |
| *Klebsiella pneumoniae* | norfloxacin | 194 | 29 | 8 | 0 | 7 | 1.30 | 6.2 |
| *Klebsiella pneumoniae* | ofloxacin | 71 | 44 | 15 | 0 | 0 | 15.62 | 1.0 |
| *Klebsiella pneumoniae* | piperacillin/tazobactam | 1,452 | 868 | 76 | 8 | 0 | 1.69 | 45.1 |
| *Klebsiella pneumoniae* | polymyxin b | 100 | 48 | 0 | 0 | 0 | 0.00 | none |
| *Klebsiella pneumoniae* | streptomycin | 80 | 54 | 0 | 0 | 0 | 0.00 | none |
| *Klebsiella pneumoniae* | sulfamethoxazole/trimethoprim | 1,546 | 1,042 | 23 | 21 | 0 | 0.00 | 46,379.8 |
| *Klebsiella pneumoniae* | tetracycline | 741 | 387 | 76 | 47 | 0 | 0.60 | 125.6 |
| *Klebsiella pneumoniae* | ticarcillin/clavulanate | 250 | 107 | 23 | 2 | 0 | 5.76 | 4.0 |
| *Klebsiella pneumoniae* | tigecycline | 924 | 257 | 95 | 30 | 0 | 38.17 | 2.5 |
| *Klebsiella pneumoniae* | tobramycin | 1,230 | 630 | 54 | 37 | 0 | 0.02 | 2,656.5 |
| *Klebsiella pneumoniae* | trimethoprim | 258 | 112 | 5 | 4 | 0 | 0.00 | 12,614.4 |
| *Klebsiella quasipneumoniae* | cefotaxime | 101 | 62 | 3 | 0 | 0 | 0.00 | 29,552.7 |
| *Klebsiella quasipneumoniae* | ceftazidime | 102 | 59 | 2 | 0 | 0 | 0.00 | 149,597.0 |
| *Klebsiella quasipneumoniae* | sulfamethoxazole/trimethoprim | 101 | 57 | 5 | 4 | 0 | 0.00 | 3,393.8 |
| *Mycobacterium tuberculosis* | amikacin | 1,813 | 548 | 67 | 67 | 0 | 0.15 | 445.9 |
| *Mycobacterium tuberculosis* | amoxicillin | 57 | 25 | 3 | 0 | 0 | 0.07 | 43.1 |
| *Mycobacterium tuberculosis* | capreomycin | 670 | 422 | 35 | no rule | 0 | 0.11 | 309.7 |
| *Mycobacterium tuberculosis* | ciprofloxacin | 101 | 73 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Mycobacterium tuberculosis* | clofazimine | 1,303 | 193 | 36 | no rule | 1 | 23.25 | 1.5 |
| *Mycobacterium tuberculosis* | cycloserine | 401 | 119 | 27 | no rule | 0 | 4.05 | 6.7 |
| *Mycobacterium tuberculosis* | delamanid | 1,248 | 448 | 33 | no rule | 0 | 4.07 | 8.1 |
| *Mycobacterium tuberculosis* | ethambutol | 1,947 | 1,338 | 164 | no rule | 0 | 0.65 | 251.5 |
| *Mycobacterium tuberculosis* | ethionamide | 1,774 | 915 | 143 | no rule | 0 | 4.55 | 31.4 |
| *Mycobacterium tuberculosis* | isoniazid | 1,954 | 1,388 | 55 | no rule | 0 | 0.00 | 45,291.7 |
| *Mycobacterium tuberculosis* | kanamycin | 1,800 | 617 | 86 | 86 | 0 | 0.23 | 372.0 |
| *Mycobacterium tuberculosis* | levofloxacin | 1,248 | 753 | 179 | 0 | 0 | 4.58 | 39.1 |
| *Mycobacterium tuberculosis* | linezolid | 1,248 | 202 | 24 | 24 | 7 | 8.11 | 3.0 |
| *Mycobacterium tuberculosis* | moxifloxacin | 1,524 | 733 | 125 | 0 | 0 | 0.19 | 645.7 |
| *Mycobacterium tuberculosis* | nicotinamide | 166 | 83 | 19 | no rule | 0 | 1.72 | 11.1 |
| *Mycobacterium tuberculosis* | ofloxacin | 643 | 355 | 33 | 0 | 0 | 0.01 | 4,431.6 |
| *Mycobacterium tuberculosis* | para-aminosalicylic acid | 1,650 | 429 | 74 | no rule | 0 | 3.34 | 22.2 |
| *Mycobacterium tuberculosis* | prothionamide | 170 | 92 | 19 | no rule | 0 | 0.10 | 182.3 |
| *Mycobacterium tuberculosis* | pyrazinamide | 492 | 381 | 6 | no rule | 0 | 0.29 | 20.7 |
| *Mycobacterium tuberculosis* | rifabutin | 1,317 | 665 | 23 | no rule | 0 | 0.00 | 20,118.5 |
| *Mycobacterium tuberculosis* | rifampin | 1,967 | 1,365 | 73 | no rule | 0 | 0.00 | 20,926.1 |
| *Mycobacterium tuberculosis* | streptomycin | 698 | 539 | 15 | 15 | 0 | 0.06 | 256.7 |
| *Neisseria gonorrhoeae* | azithromycin | 1,899 | 408 | 142 | 141 | 0 | 12.57 | 11.3 |
| *Neisseria gonorrhoeae* | cefixime | 600 | 103 | 6 | 6 | 0 | 0.06 | 101.7 |
| *Neisseria gonorrhoeae* | ceftriaxone | 1,926 | 40 | 22 | 22 | 17 | 20.35 | 1.1 |
| *Neisseria gonorrhoeae* | ciprofloxacin | 1,746 | 1,086 | 5 | 5 | 0 | 0.00 | 67,643,894.2 |
| *Neisseria gonorrhoeae* | erythromycin | 179 | 98 | 3 | 3 | 0 | 0.00 | 1,819.7 |
| *Neisseria gonorrhoeae* | penicillin | 1,354 | 888 | 9 | 9 | 0 | 3.20 | 2.8 |
| *Neisseria gonorrhoeae* | tetracycline | 1,371 | 359 | 7 | 7 | 0 | 0.00 | 1,758.2 |
| *Pseudomonas aeruginosa* | amikacin | 708 | 127 | 13 | 11 | 0 | 0.34 | 37.9 |
| *Pseudomonas aeruginosa* | aztreonam | 235 | 114 | 9 | 0 | 5 | 1.78 | 5.1 |
| *Pseudomonas aeruginosa* | cefepime | 245 | 122 | 26 | 0 | 13 | 8.78 | 3.0 |
| *Pseudomonas aeruginosa* | ceftazidime | 833 | 372 | 51 | 0 | 0 | 1.27 | 40.0 |
| *Pseudomonas aeruginosa* | ceftazidime/avibactam | 309 | 139 | 70 | 33 | 0 | 15.73 | 4.5 |
| *Pseudomonas aeruginosa* | ceftolozane/tazobactam | 259 | 132 | 14 | 7 | 7 | 7.57 | 1.8 |
| *Pseudomonas aeruginosa* | ciprofloxacin | 703 | 391 | 29 | 28 | 0 | 0.01 | 5,752.4 |
| *Pseudomonas aeruginosa* | gentamicin | 346 | 101 | 17 | 14 | 0 | 0.04 | 436.0 |
| *Pseudomonas aeruginosa* | imipenem | 447 | 235 | 28 | 28 | 0 | 0.31 | 89.6 |
| *Pseudomonas aeruginosa* | levofloxacin | 702 | 330 | 18 | 18 | 0 | 0.01 | 2,645.2 |
| *Pseudomonas aeruginosa* | meropenem | 883 | 488 | 84 | 84 | 0 | 0.87 | 96.6 |
| *Pseudomonas aeruginosa* | piperacillin/tazobactam | 210 | 87 | 24 | 0 | 0 | 6.83 | 3.5 |
| *Pseudomonas aeruginosa* | tobramycin | 645 | 194 | 27 | 25 | 0 | 0.03 | 844.4 |
| *Salmonella enterica* | amoxicillin/clavulanate | 799 | 471 | 38 | 29 | 0 | 0.00 | 19,716.3 |
| *Salmonella enterica* | ampicillin | 1,997 | 1,406 | 53 | 21 | 0 | 0.10 | 555.2 |
| *Salmonella enterica* | ampicillin/sulbactam | 183 | 126 | 7 | 3 | 0 | 0.06 | 117.4 |
| *Salmonella enterica* | azithromycin | 1,355 | 184 | 29 | 23 | 0 | 0.04 | 662.3 |
| *Salmonella enterica* | cefazolin | 185 | 119 | 5 | 4 | 0 | 1.37 | 3.7 |
| *Salmonella enterica* | cefotaxime | 941 | 319 | 38 | 27 | 3 | 2.64 | 14.4 |
| *Salmonella enterica* | cefoxitin | 1,080 | 395 | 13 | 10 | 0 | 0.00 | 238,328.2 |
| *Salmonella enterica* | ceftazidime | 783 | 30 | 13 | 6 | 7 | 6.19 | 2.1 |
| *Salmonella enterica* | ceftiofur | 689 | 461 | 10 | 4 | 0 | 0.00 | 34,814.9 |
| *Salmonella enterica* | ceftriaxone | 1,155 | 568 | 72 | 8 | 0 | 0.04 | 2,039.3 |
| *Salmonella enterica* | chloramphenicol | 1,989 | 713 | 44 | 42 | 0 | 1.25 | 35.1 |
| *Salmonella enterica* | ciprofloxacin | 1,503 | 353 | 58 | 36 | 0 | 0.06 | 910.1 |
| *Salmonella enterica* | gentamicin | 1,675 | 415 | 22 | 16 | 0 | 0.69 | 31.9 |
| *Salmonella enterica* | kanamycin | 584 | 154 | 36 | 31 | 0 | 0.50 | 71.9 |
| *Salmonella enterica* | nalidixic acid | 1,599 | 500 | 65 | 40 | 0 | 0.08 | 824.5 |
| *Salmonella enterica* | pefloxacin | 348 | 142 | 2 | no rule | 0 | 0.00 | 2,509,795.7 |
| *Salmonella enterica* | spectinomycin | 396 | 208 | 14 | 14 | 0 | 0.09 | 158.1 |
| *Salmonella enterica* | streptomycin | 1,269 | 640 | 22 | 16 | 0 | 0.04 | 570.9 |
| *Salmonella enterica* | sulfamethoxazole | 275 | 92 | 12 | 5 | 3 | 0.46 | 26.1 |
| *Salmonella enterica* | sulfamethoxazole/trimethoprim | 1,145 | 296 | 52 | 32 | 12 | 7.81 | 6.7 |
| *Salmonella enterica* | sulfathiazole | 225 | 175 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Salmonella enterica* | sulfisoxazole | 544 | 443 | 4 | 1 | 0 | 0.00 | 9,426.3 |
| *Salmonella enterica* | tetracycline | 1,869 | 1,118 | 109 | 48 | 0 | 0.62 | 174.6 |
| *Salmonella enterica* | trimethoprim | 634 | 256 | 20 | 1 | 0 | 0.02 | 1,073.2 |
| *Shigella flexneri* | ampicillin | 358 | 318 | 2 | 0 | 0 | 0.00 | 66,889.5 |
| *Shigella flexneri* | azithromycin | 361 | 186 | 6 | 6 | 0 | 0.00 | 32,467.5 |
| *Shigella flexneri* | ciprofloxacin | 361 | 33 | 2 | 0 | 0 | 0.04 | 49.4 |
| *Shigella flexneri* | trimethoprim | 362 | 280 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Shigella sonnei* | amikacin | 580 | 270 | 1 | 1 | 0 | 0.00 | 479,293.8 |
| *Shigella sonnei* | ampicillin | 1,674 | 848 | 46 | 0 | 0 | 0.01 | 5,497.5 |
| *Shigella sonnei* | azithromycin | 1,307 | 602 | 27 | 26 | 0 | 0.00 | 222,419.8 |
| *Shigella sonnei* | cefoperazone | 257 | 226 | 6 | 0 | 0 | 0.00 | 1,267,618.0 |
| *Shigella sonnei* | cefotaxime | 593 | 112 | 6 | 0 | 0 | 0.00 | 4,789.9 |
| *Shigella sonnei* | cefoxitin | 593 | 105 | 4 | 0 | 0 | 0.00 | 894.9 |
| *Shigella sonnei* | cefpirome | 333 | 97 | 3 | no rule | 0 | 0.00 | 1,721.5 |
| *Shigella sonnei* | ceftazidime | 582 | 98 | 4 | 0 | 0 | 0.00 | 3,596.0 |
| *Shigella sonnei* | ceftriaxone | 1,163 | 439 | 7 | 0 | 0 | 0.00 | 3,887,115.0 |
| *Shigella sonnei* | ciprofloxacin | 1,397 | 606 | 13 | 8 | 0 | 0.00 | 236,317.3 |
| *Shigella sonnei* | colistin | 268 | 35 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Shigella sonnei* | ertapenem | 338 | 101 | 6 | 6 | 0 | 0.00 | 1,404.6 |
| *Shigella sonnei* | gentamicin | 1,288 | 513 | 1 | 1 | 0 | 0.00 | 117,015.0 |
| *Shigella sonnei* | nalidixic acid | 854 | 544 | 8 | 3 | 0 | 0.00 | 261,731.1 |
| *Shigella sonnei* | streptomycin | 595 | 471 | 18 | 16 | 0 | 0.00 | 18,705.2 |
| *Shigella sonnei* | sulfamethoxazole/trimethoprim | 617 | 524 | 21 | 18 | 0 | 0.00 | 10,776.0 |
| *Shigella sonnei* | temocillin | 333 | 97 | 3 | no rule | 0 | 0.00 | 1,721.5 |
| *Shigella sonnei* | tetracycline | 878 | 704 | 31 | 31 | 0 | 0.00 | 86,693.9 |
| *Shigella sonnei* | tobramycin | 541 | 462 | 1 | 1 | 0 | 0.00 | 1,467.8 |
| *Shigella sonnei* | trimethoprim | 1,104 | 1,036 | 9 | 3 | 0 | 0.00 | 71,123.4 |
| *Staphylococcus aureus* | cefoxitin | 1,407 | 730 | 45 | 45 | 0 | 0.01 | 6,806.5 |
| *Staphylococcus aureus* | chloramphenicol | 1,038 | 123 | 42 | 30 | 0 | 2.06 | 20.4 |
| *Staphylococcus aureus* | ciprofloxacin | 1,624 | 625 | 45 | 45 | 0 | 0.02 | 1,901.2 |
| *Staphylococcus aureus* | clindamycin | 1,543 | 693 | 294 | 40 | 229 | 194.55 | 1.5 |
| *Staphylococcus aureus* | daptomycin | 904 | 47 | 30 | 30 | 17 | 30.94 | 1.0 |
| *Staphylococcus aureus* | erythromycin | 1,957 | 948 | 58 | 56 | 10 | 5.88 | 9.9 |
| *Staphylococcus aureus* | fusidic acid | 1,101 | 283 | 13 | 12 | 0 | 0.10 | 129.4 |
| *Staphylococcus aureus* | gentamicin | 1,913 | 367 | 16 | 12 | 0 | 0.01 | 2,517.3 |
| *Staphylococcus aureus* | mupirocin | 736 | 51 | 9 | 8 | 0 | 0.03 | 290.0 |
| *Staphylococcus aureus* | oxacillin | 1,318 | 568 | 41 | 41 | 0 | 0.00 | 19,280.7 |
| *Staphylococcus aureus* | penicillin | 1,508 | 1,252 | 35 | 35 | 0 | 0.00 | 100,543.2 |
| *Staphylococcus aureus* | rifampin | 931 | 49 | 2 | no rule | 0 | 0.00 | 10,443.5 |
| *Staphylococcus aureus* | sulfamethoxazole/trimethoprim | 1,389 | 285 | 31 | 31 | 0 | 0.01 | 3,263.5 |
| *Staphylococcus aureus* | tetracycline | 1,647 | 408 | 87 | 0 | 0 | 1.86 | 46.9 |
| *Staphylococcus haemolyticus* | ciprofloxacin | 124 | 78 | 3 | 3 | 0 | 0.00 | 31,727.2 |
| *Staphylococcus haemolyticus* | fusidic acid | 114 | 39 | 6 | 2 | 0 | 3.27 | 1.8 |
| *Staphylococcus haemolyticus* | tetracycline | 104 | 37 | 10 | 7 | 0 | 0.46 | 22.0 |
| *Streptococcus agalactiae* | clindamycin | 141 | 115 | 0 | 0 | 0 | 0.00 | 0.0 |
| *Streptococcus agalactiae* | erythromycin | 141 | 116 | 3 | 0 | 0 | 0.00 | 58,940.1 |
| *Streptococcus agalactiae* | tetracycline | 141 | 104 | 0 | 0 | 0 | 0.00 | none |
| *Streptococcus pneumoniae* | azithromycin | 1,297 | 147 | 9 | 9 | 0 | 0.02 | 517.0 |
| *Streptococcus pneumoniae* | cefepime | 1,602 | 49 | 10 | 10 | 10 | 4.43 | 2.3 |
| *Streptococcus pneumoniae* | cefotaxime | 1,884 | 43 | 5 | 5 | 0 | 0.11 | 46.0 |
| *Streptococcus pneumoniae* | ceftriaxone | 1,839 | 56 | 4 | 4 | 1 | 1.86 | 2.2 |
| *Streptococcus pneumoniae* | cefuroxime | 1,732 | 184 | 3 | 3 | 0 | 0.00 | 3,979,960.9 |
| *Streptococcus pneumoniae* | chloramphenicol | 1,859 | 134 | 24 | 24 | 0 | 0.59 | 40.5 |
| *Streptococcus pneumoniae* | clindamycin | 1,894 | 166 | 11 | 11 | 0 | 0.00 | 30,170.6 |
| *Streptococcus pneumoniae* | erythromycin | 1,991 | 460 | 22 | 22 | 0 | 0.00 | 49,922.0 |
| *Streptococcus pneumoniae* | meropenem | 1,808 | 72 | 3 | 3 | 1 | 0.03 | 88.7 |
| *Streptococcus pneumoniae* | penicillin | 1,690 | 196 | 1 | 1 | 1 | 0.00 | 6,751.3 |
| *Streptococcus pneumoniae* | sulfamethoxazole/trimethoprim | 1,777 | 343 | 10 | 10 | 0 | 0.00 | 13,097.0 |
| *Streptococcus pneumoniae* | tetracycline | 1,971 | 314 | 28 | 19 | 0 | 0.00 | 12,040.7 |
| *Vibrio cholerae* | nalidixic acid | 91 | 49 | 5 | 5 | 0 | 0.00 | 6,405.3 |
| *Vibrio cholerae* | sulfamethoxazole | 74 | 37 | 0 | 0 | 0 | 0.00 | none |
| *Vibrio cholerae* | trimethoprim | 74 | 40 | 1 | 1 | 0 | 0.00 | 2,560,000.0 |

Skipped, with fewer than 2 learners on the lineage split: *Enterococcus faecium* with chloramphenicol, *Salmonella enterica* with temocillin, *Staphylococcus aureus* with ceftaroline, *Staphylococcus aureus* with levofloxacin, *Staphylococcus aureus* with moxifloxacin and *Staphylococcus aureus* with vancomycin.
