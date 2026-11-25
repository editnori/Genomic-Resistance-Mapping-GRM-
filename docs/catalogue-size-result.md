# Catalogue size result

Is the set of ways *Mycobacterium tuberculosis* can resist an antibiotic small for some antibiotics and large for others, and does the WHO catalogue's size show it? The lineage curve's page (`docs/lineage-curve-result.md`) counts, per antibiotic, the distinct graded variants of the WHO catalogue among the resistant genomes (`catalogue_variants` in `reference/lineage-curve/page-endpoints.tsv`). An antibiotic studied for decades may have both a large catalogue and many resistant genomes on deposit, so that count may record attention rather than how many ways the antibiotic can fail. This page asks whether the count tracks the resistant share of an antibiotic's genomes or only how many resistant genomes there were to study, and says so plainly where the numbers cannot separate the two.

It fits nothing, trains nothing and scores no genome again. It reads the lineage curve's committed table beside three counts of each release report: its genomes, its resistant genomes and its susceptible genomes.

## Rule

This rule was written and committed before the release's counts were read, and it is not changed after.

- **Antibiotics.** The rows of `reference/lineage-curve/page-endpoints.tsv` whose `catalogue_variants` is 1 or more and whose `in_M` is 0. A row with `in_M` 1 waits on the maintainer's label rule: it stands in the table and is left out of every correlation. A row the tool does not read, or with no graded variant, stands in the table and is left out too. Every correlation of a set uses the same antibiotics.
- **Pairs.** The Spearman rank correlation of `catalogue_variants` against each of: the resistant share (`counts.resistant` over `counts.samples` of the antibiotic's release report), the resistant genomes (`counts.resistant`) and the genomes (`counts.samples`). Tied values take the average of their ranks, and the correlation is Pearson's over the ranks.
- **Leave one out.** For each pair, the correlation over the set with one antibiotic dropped, for each antibiotic in turn: the smallest and the largest, each with the antibiotic dropped.
- **Decision.** Let R be the correlation with the resistant genomes and S the one with the resistant share, each in absolute value and unrounded. If R is as large as S or larger, the answer is that catalogue size tracks how many resistant genomes there were to study, and target size is not separable here. If S is larger than R by 0.15 or more, the answer is that the number of graded variants tracks the resistant share, in one species, one catalogue and one database, and no causal reading follows. If S is larger than R by less than 0.15, the answer is that these numbers cannot separate the two. A difference within 0.000000001 of a bound counts as meeting it.
- **Beside.** The correlation with the genomes is reported and decides nothing. The same three pairs over the antibiotics of the set that the tool reads by their own drug (kind `primary`) are reported and decide nothing: a proxy counts another drug's graded variants.
- **Not done.** No p-value, no test, no fit, no fitted line and no model.

## Reading

Over the 16 antibiotics the rule reads, the number of graded variants has a Spearman correlation of 0.522 with the resistant share and 0.470 with the resistant genomes. The share's is the larger, by less than the 0.15 the rule asked for, so the rule's answer is that these numbers cannot separate the two. On this reading, the size of the catalogue for an antibiotic may track the resistant share or how many resistant genomes there were to study, and nothing here says which.

The rule declared three outcomes before the read. If the correlation with the resistant genomes was as strong as the share's or stronger, the answer would be that catalogue size tracks how many resistant genomes there were to study. If the share's was the stronger by 0.15 or more, the answer would be that the number of graded variants tracks the resistant share. If the share's was the stronger by less than 0.15, the answer would be that these numbers cannot separate the two. The numbers selected the third: the share's lead fell short, 0.052 against the 0.15 the second outcome needed.

With one antibiotic dropped, the share's correlation runs from 0.424 (linezolid dropped) to 0.625 (ethionamide dropped), and the resistant genomes' from 0.394 (ethionamide dropped) to 0.607 (pyrazinamide dropped). Each range is wider than the gap between the two correlations. A rank correlation over 16 antibiotics of one species is weak evidence. The correlation with the genomes is 0.068, from -0.072 (rifampin dropped) to 0.201 (pyrazinamide dropped), and decides nothing.

Beside the rule stands a second set, the 12 antibiotics the tool reads by their own drug. It was declared before the read, in the rule, and it decides nothing. The four proxy antibiotics are why the rule's set holds 16 and this one 12: the rule reads every antibiotic with a graded variant and not in set M, and this set leaves out the four that count another drug's graded variants, as ciprofloxacin and ofloxacin count those of levofloxacin and moxifloxacin, rifabutin those of rifampicin and prothionamide those of ethionamide. Over the 12, the correlation with the resistant share is 0.698, from 0.620 to 0.822 with one dropped, and with the resistant genomes 0.477, from 0.333 to 0.753, so on that set the share's lead passes 0.15. The table marks its rows "beside, decides nothing", and the rule's answer stands.

## Antibiotics

One row per report of *Mycobacterium tuberculosis*, species and antibiotic as `reference/lineage-curve/page-endpoints.tsv` names them. Genomes, resistant and susceptible are the release report's `counts.samples`, `counts.resistant` and `counts.susceptible`, and the resistant share is the resistant genomes over the genomes. The next three columns are carried from the lineage curve's page unchanged, as it defines them. Graded variants (`catalogue_variants`) counts the distinct graded variants among the resistant genomes, where a graded variant of an antibiotic is a variant the tool reports for the antibiotic's drug, passing the tool's own filter, with an entry graded "Assoc w R" or "Assoc w R - Interim". Genes named (`genes_named`) counts the genes the lineage curve's Genes table names for the antibiotic, a gene being named only where 5 or more resistant genomes carry one of its graded variants. No graded variant (`without_determinant`, `share_without`) counts the resistant genomes for which the tool reports no graded variant for the antibiotic's drug, with their share. Each count is a count: no variant and no gene is named here.

| Antibiotic | Set | Kind | Genomes | Resistant | Susceptible | Resistant share | Graded variants | Genes named | No graded variant (share) | In the correlations |
|---|---|---|---|---|---|---|---|---|---|---|
| amikacin | not M | primary | 1,813 | 548 | 1,265 | 0.302 | 4 | 2 | 116 (0.212) | yes |
| amoxicillin | not M | none | 57 | 25 | 32 | 0.439 | none | none | none | no: not read by the tool |
| capreomycin | not M | primary | 670 | 422 | 248 | 0.630 | 10 | 2 | 75 (0.178) | yes |
| ciprofloxacin | not M | proxy | 101 | 73 | 28 | 0.723 | 8 | 1 | 8 (0.110) | yes |
| clofazimine | M | primary | 1,303 | 193 | 1,110 | 0.148 | 5 | 1 | 188 (0.974) | no: set M |
| cycloserine | not M | primary | 401 | 119 | 282 | 0.297 | 0 | 0 | 119 (1.000) | no: no graded variant |
| delamanid | not M | primary | 1,248 | 448 | 800 | 0.359 | 3 | 0 | 444 (0.991) | yes |
| ethambutol | not M | primary | 1,947 | 1,338 | 609 | 0.687 | 12 | 1 | 347 (0.259) | yes |
| ethionamide | not M | primary | 1,774 | 915 | 859 | 0.516 | 136 | 2 | 259 (0.283) | yes |
| isoniazid | M | primary | 1,954 | 1,388 | 566 | 0.710 | 17 | 2 | 77 (0.055) | no: set M |
| kanamycin | not M | primary | 1,800 | 617 | 1,183 | 0.343 | 6 | 2 | 128 (0.207) | yes |
| levofloxacin | not M | primary | 1,248 | 753 | 495 | 0.603 | 16 | 2 | 311 (0.413) | yes |
| linezolid | not M | primary | 1,248 | 202 | 1,046 | 0.162 | 3 | 2 | 178 (0.881) | yes |
| moxifloxacin | not M | primary | 1,524 | 733 | 791 | 0.481 | 15 | 2 | 136 (0.186) | yes |
| nicotinamide | not M | none | 166 | 83 | 83 | 0.500 | none | none | none | no: not read by the tool |
| ofloxacin | not M | proxy | 643 | 355 | 288 | 0.552 | 12 | 2 | 34 (0.096) | yes |
| para-aminosalicylic acid | not M | primary | 1,650 | 429 | 1,221 | 0.260 | 0 | 0 | 429 (1.000) | no: no graded variant |
| prothionamide | not M | proxy | 170 | 92 | 78 | 0.541 | 14 | 2 | 26 (0.283) | yes |
| pyrazinamide | not M | primary | 492 | 381 | 111 | 0.774 | 94 | 1 | 60 (0.157) | yes |
| rifabutin | not M | proxy | 1,317 | 665 | 652 | 0.505 | 27 | 1 | 20 (0.030) | yes |
| rifampin | not M | primary | 1,967 | 1,365 | 602 | 0.694 | 40 | 1 | 72 (0.053) | yes |
| streptomycin | not M | primary | 698 | 539 | 159 | 0.772 | 15 | 3 | 44 (0.082) | yes |

Of the 22 reports, 2 are not read by the tool (amoxicillin and nicotinamide), 2 have no graded variant (cycloserine and para-aminosalicylic acid) and 2 are in set M (clofazimine and isoniazid), whose labels wait for the maintainer's label rule. All six stand in the table with their counts and are left out of the correlations, which read the other 16.

## Correlations

One row per set and pair, from `reference/catalogue-size/correlations.tsv`. Each correlation is the Spearman rank correlation of the graded variants against the named column over the set's antibiotics; one dropped gives the smallest and largest value the correlation takes with one antibiotic left out, and which one. The comparison is the rule's, over the set's correlations with the resistant genomes and the resistant share. The role marks each row as the rule's set or as beside it and deciding nothing; only the rule's set decides.

| Set | Role | Graded variants against | n | Spearman | One dropped, smallest (dropped) | One dropped, largest (dropped) | Decides | The comparison gives |
|---|---|---|---|---|---|---|---|---|
| rule | the rule's set | resistant share | 16 | 0.522 | 0.424 (linezolid) | 0.625 (ethionamide) | yes | neither by the margin |
| rule | the rule's set | resistant genomes | 16 | 0.470 | 0.394 (ethionamide) | 0.607 (pyrazinamide) | yes | neither by the margin |
| rule | the rule's set | genomes | 16 | 0.068 | -0.072 (rifampin) | 0.201 (pyrazinamide) | no | neither by the margin |
| primary only | beside, decides nothing | resistant share | 12 | 0.698 | 0.620 (linezolid) | 0.822 (ethionamide) | no | tracks share |
| primary only | beside, decides nothing | resistant genomes | 12 | 0.477 | 0.333 (linezolid) | 0.753 (pyrazinamide) | no | tracks share |
| primary only | beside, decides nothing | genomes | 12 | -0.021 | -0.194 (rifampin) | 0.212 (pyrazinamide) | no | tracks share |

The antibiotics of each set, the same for its three pairs:

| Set | n | Antibiotics |
|---|---|---|
| rule | 16 | amikacin, capreomycin, ciprofloxacin, delamanid, ethambutol, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, ofloxacin, prothionamide, pyrazinamide, rifabutin, rifampin, streptomycin |
| primary only | 12 | amikacin, capreomycin, delamanid, ethambutol, ethionamide, kanamycin, levofloxacin, linezolid, moxifloxacin, pyrazinamide, rifampin, streptomycin |

## Limits

**What the count is.** `catalogue_variants` counts the graded variants seen among the resistant genomes of these deposits, not the catalogue's entries for the drug. More resistant genomes give more chances to see a variant, so the count can grow with them whatever the number of ways the antibiotic can fail.

**Attention and target size.** The resistant genomes and the resistant share both carry how the world studied an antibiotic: which genomes were sequenced, tested and deposited. Telling attention from the number of ways an antibiotic can fail needs a source of evidence GRM does not hold.

**One species, one catalogue, one database.** The reading covers *Mycobacterium tuberculosis* alone, the WHO catalogue as TBProfiler 6.6.3 ships it (`tbdb` commit `3ed1cb99`), and BV-BRC deposits chosen for having lab results, capped at 2,000 a species: deposits, not a sample of the world's genomes. A share is no prevalence, and no causal reading follows from a correlation across antibiotics.

**What was knowable before the rule.** One of the three correlations of the rule's set could have been computed before the rule was committed: `resistant_with_result` in `reference/lineage-curve/page-endpoints.tsv`, a public table read in full before then, equals the release report's resistant count on every row the tool reads, so the correlation with the resistant genomes was knowable. The comparison the rule turns on was not, because the resistant share needed genome counts that had not been read, and the 0.15 margin and the first two outcomes were set in the brief this reading was given, before its work began; the third outcome is the gap that brief's "only if" leaves between them.

**Antibiotics that share calls.** The four proxies count another drug's graded variants, and antibiotics whose graded variants sit in the same genes are not independent of each other. A rank correlation counts each antibiotic as one.

**Care.** Nothing here says which antibiotic is a better or safer choice in care.

**Set M.** The reports of clofazimine and isoniazid wait for the maintainer's label rule. They stand in the table and are left out of the correlations.

**Rows.** The record bars a row per genome; a row per antibiotic is an aggregate and keeps its genome count. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Method

**Order.** The rule above was committed first. The script `scripts/mechanism/catalogue_size.py` and its SYNTHETIC tests were committed next and passed on one machine before the read. The read then ran once, at that commit, and did not stop. This page and `reference/catalogue-size/` were written from its output. `write` then ran four more times over the same reading, with no `read`: twice after a first review, at the work tree and at the commit, to add the role column to `correlations.tsv`, and twice after a second, to move guard 5 onto the full list of genome ids. A later review ran it again over the same reading and got the committed record byte for byte. With the role column taken out, `correlations.tsv` is byte for byte the one the first run wrote. `antibiotics.tsv` and `reading.txt` did not change, and `guards.txt` changed only in guard 5's line.

**Reading.** The commands, with the release's report folder as REPORTS, its species files as SPECIES_DIR, the reading's folder as READING and the record as OUT, copied to `reference/catalogue-size/` as written:

```
python3 -I scripts/mechanism/catalogue_size.py read REPORTS reference/lineage-curve/page-endpoints.tsv READING
python3 -I scripts/mechanism/catalogue_size.py write READING reference/lineage-curve/page-endpoints.tsv SPECIES_DIR OUT
python3 -I scripts/mechanism/catalogue_size.py print OUT
```

`read` takes three counts from each release report and nothing else; a report's `samples` list, whose entries carry the genome ids, is dropped as the file is parsed. `reading.txt` gives the first command with the sha256 of what it wrote. Guard 5 holds the genome ids of the release's species files in memory only. It searches every written file for each id in full, wherever it stands, and counts tokens of a genome id's shape beside it as a second net. It prints counts only.

**Guards.** Each printed counts only and would have stopped the run, as would a set of under 3 antibiotics or a column whose values are all tied. All six passed, in the order the script runs them:

| Guard | Printed |
|---|---|
| 1. match | 22 page rows, 22 reports of the species, 22 matched, 0 page rows with no report, 0 reports with no page row |
| 2. sums | 22 of 22 reports add up |
| 6. resistant | 20 of 20 reports the tool reads agree with the page |
| 3. sets | rule: resistant_share n=16, resistant n=16, genomes n=16; antibiotics: amikacin; capreomycin; ciprofloxacin; delamanid; ethambutol; ethionamide; kanamycin; levofloxacin; linezolid; moxifloxacin; ofloxacin; prothionamide; pyrazinamide; rifabutin; rifampin; streptomycin |
| 3. sets | primary only: resistant_share n=12, resistant n=12, genomes n=12; antibiotics: amikacin; capreomycin; delamanid; ethambutol; ethionamide; kanamycin; levofloxacin; linezolid; moxifloxacin; pyrazinamide; rifampin; streptomycin |
| 4. rank correlation | 3 of 3 hand-made cases give their known answers |
| 5. ids | 8 files; 0 hits of 23757 genome ids of 23 species files, each searched for in full; 0 tokens of a genome id's shape; 0 BioProject, sequence or position hits |

**Tests.** `tests/mechanism/test_catalogue_size.py` holds 36 SYNTHETIC tests of the script. `tests/mechanism/test_catalogue_faults.py` runs `scripts/mechanism/catalogue_faults.py`, which plants each fault below in the script in turn and runs the test meant to catch it. `tests/mechanism/test_catalogue_page.py` holds 12 tests of this page against `reference/catalogue-size/`.

| Fault planted | Test that caught it |
|---|---|
| K1. Tied values ranked by their order, not by the average of their ranks | `Spearman.test_a_case_with_ties_gives_its_known_answer` |
| K2. Pearson's correlation over the values in place of their ranks | `Spearman.test_a_case_with_no_ties_gives_its_known_answer` |
| K3. The rule compares the signed correlations | `Rule.test_the_rule_reads_absolute_values` |
| K4. The margin of 0.15 dropped | `Rule.test_a_share_stronger_by_less_than_the_margin_selects_neither` |
| K5. Equal strength read as the share's | `Rule.test_a_resistant_correlation_as_strong_as_the_share_selects_resistant_genomes` |
| K6. A row of set M in the correlations | `Write.test_set_m_is_in_the_table_and_out_of_the_correlations` |
| K7. A row with no graded variant in the correlations | `Write.test_an_antibiotic_with_no_graded_variant_is_left_out` |
| K8. A page antibiotic with no release report let through | `Read.test_a_page_antibiotic_with_no_report_stops` |
| K9. Counts that do not add up let through | `Read.test_counts_that_do_not_add_up_stop` |
| K10. A resistant count unlike the lineage curve page's let through | `Read.test_a_resistant_count_unlike_the_page_stops` |
| K11. A pair over another set of antibiotics let through | `Write.test_a_pair_over_another_set_stops` |
| K12. The hand-made cases skipped | `Write.test_a_wrong_rank_correlation_stops_before_any_correlation` |
| K13. The scan for genome ids skipped | `Write.test_a_genome_id_in_an_output_stops` |
| K14. The leave-one-out largest taken as the smallest | `Spearman.test_leave_one_out_gives_the_smallest_and_largest_with_the_antibiotic_dropped` |
| K15. A report's samples list kept when the report is parsed | `Read.test_the_samples_list_is_never_kept` |
| K16. The resistant share taken over the susceptible genomes | `Write.test_the_resistant_share_is_over_all_genomes` |
| K17. The id list searched only through tokens already shaped like an id | `Write.test_an_id_glued_to_a_letter_or_an_underscore_is_caught` |

All 17 were caught.

**Files.** The reading's tables are in `reference/catalogue-size/`. `scripts/mechanism/catalogue_size.py print` prints this page's tables from those files, rounded to 3 decimals.
