# Mechanism year result

This page reads, once, how the mix of catalogue strata among the release's genomes differs between the years the genomes were collected, and how the five learners' mean vote calls the resistant genomes of each stratum in each year. It fits nothing, trains nothing and scores no genome again: every figure comes from the held-out scores, labels and catalogue calls stored in the 274 revision 14 reports. The maintainer named this read on 23 Nov 2026 as one of three that need no fit.

The four strata of block M1 were not readable from the release reports, so the read runs on two groups taken from GRM's known-gene rule instead: the genomes it calls resistant and those it makes no call on (section Step 1).

A share over collection years describes the deposits of one public database. It is no prevalence and no trend in resistance, in the world or in any population. No line is fitted and no slope is given. A split by lineage is a retrospective split inside one database, never external validation. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Reading

The four strata of M1 were not readable from the release reports. Rebuilt from each genome's catalogue groups by GRM's known-gene rule, GRM's own reader of the AMRFinderPlus catalogue, they put no genome in point only or both, and 48,749 resistant genomes in gene only against M1's 39,814. The rebuild was not tuned. The read runs on the reader's own call instead: `reader_call`, the genomes the rule calls resistant, and `reader_no_call`, those it makes no call on. No genome of the primary set has a null call, so `reader_no_rule` is empty.

Over the 204 datasets of the primary set, the reader makes no call on 15,514 of the 62,878 resistant genomes. The five learners' mean vote on the lineage split calls 11,989 of those 15,514 resistant (recall 0.773), and 41,800 of the 47,364 resistant genomes the reader calls (0.883).

Among the 29 years that are not thin, 1985 and 1994 to 2021, the share of `reader_no_call` among a year's resistant genomes runs from 0.127 (2011, 275 of 2,163) to 0.502 (2019, 126 of 251). Within `reader_no_call`, the mean vote's recall runs from 0.317 (2019, 40 of 126) to 0.955 (2004, 257 of 269), over the 27 of those years that hold 10 or more such resistant genomes. Within `reader_call` it runs from 0.514 (1994, 18 of 35) to 0.967 (2004, 319 of 330), over all 29.

The ranges are read from the years table and say nothing about the order of the years: whether the share of `reader_no_call` grows or shrinks over the years is not read, no line is fitted and no slope is given. A share describes the deposits of one public database, not a trend in resistance in the world or in any population. 50,329 genomes carry no four-digit year and stand apart as `year_unknown`. Recall is given at the release's cutoff, beside counts; this read computes no log loss (section Limits).

`docs/mechanism-split-result.md` reads M1's strata once over the whole primary set, and its block M2 compares the AMRFinderPlus program's calls between genomes collected in 2010 to 2014 and in 2015 to 2019, one vote per BioProject. This page counts the reader's call year by year, a genome once per dataset, and tests nothing.

## Plan

The plan was fixed on 24 Nov 2026, before `scripts/mechanism/year_mix.py` read a real report. `year_mix.py plan` prints its constants.

**Question.** `docs/mechanism-split-result.md` reads the four strata of block M1 (point only, gene only, both and neither) once over the whole primary set. It does not say whether the mix changes with the collection year. This read gives, per year, each stratum's share of the year's resistant and susceptible genomes, and within each year and stratum the resistant genomes the five learners' mean vote calls. Both are measurement facts about the record GRM is scored on.

**Datasets.** The primary set is the 204 datasets of `scripts/c3-inputs/lists.tsv` with `in_A` 1 and `in_M` 0, the datasets behind the primary rows of `reference/mechanism-split/page-m1-strata.tsv`. A report in set M, outside A or not in the list is counted and not used. A primary dataset with no lineage split is counted and not used. The counts not used plus the datasets used equal the reports found.

**Step 1, the strata from the reports.** M1 built its strata from the AMRFinderPlus program's own rows. The reports hold GRM's known-gene rule instead, GRM's own reader of the AMRFinderPlus catalogue: each genome's `groups`, and per dataset `catalog.groups` with `counts_for_rule`. The rebuild is fixed here and is not changed after it is read. A genome is gene only when one of its groups counts for the dataset's rule, and neither otherwise. No field of a report marks a group as a point mutation, and the reader holds acquired genes only (`src/catalog.rs`), so the rebuild puts no genome in point only or both. The rebuild's resistant and susceptible genomes per stratum, over the primary set on the lineage split, a genome once per dataset, are set against the eight published counts. The published stratum rows must first be the same for every row of that file.

- If all eight counts are equal, the year read runs on the four strata.
- If any count differs, the rebuild is not tuned. Both sets are printed and written, and the year read runs on two groups taken straight from `samples[].catalog_resistant`: the genomes the reader calls resistant (`reader_call`) and those it does not (`reader_no_call`). A genome whose call is null is counted in a group of its own, `reader_no_rule`.

**Year.** `samples[].year` is a year when it is four ASCII digits and nothing else. Any other value, a blank or a missing year is counted as `year_unknown` and never dropped. Nothing is imputed.

**Call.** A genome is called resistant when the mean of the five learners' held-out scores on the lineage split (`scm`, `cart`, `rf`, `gbt` and `lasso`) is 0.5 or more. The cutoff, the split and the equal weights are the release's own and were not searched. No single learner is read apart.

**Tables**, in `reference/mechanism-year/`:

| File | One row per | Columns |
|---|---|---|
| `step1.tsv` | stratum of M1, and `all_strata` | `stratum`, `published_resistant`, `published_susceptible`, `rebuilt_resistant`, `rebuilt_susceptible`, `equal` |
| `years.tsv` | year and stratum | `year`, `stratum`, `genomes`, `resistant`, `susceptible`, `resistant_called_mean5`, `susceptible_called_mean5`, `recall_mean5`, `thin` |
| `totals.tsv` | stratum | `stratum`, `genomes`, `resistant`, `susceptible`, `resistant_called_mean5`, `susceptible_called_mean5`, `recall_mean5` |
| `shares.tsv` | year and stratum | `year`, `stratum`, `resistant`, `year_resistant`, `resistant_share`, `susceptible`, `year_susceptible`, `susceptible_share`, `year_genomes`, `thin` |
| `species.tsv` | species, antibiotic and stratum with 10 or more resistant genomes | `species`, `antibiotic`, `stratum`, `genomes`, `resistant`, `susceptible`, `resistant_called_mean5`, `susceptible_called_mean5`, `recall_mean5`, `year_unknown` |
| `guards.txt` | guard | the guard and its counts |

A genome counts once per dataset it appears in, so a genome tested for three antibiotics counts three times. A recall is written only over 10 or more resistant genomes, and no rate among susceptible genomes is written, only their count. A year with fewer than 100 genomes is `thin`: it keeps its rows, and the page quotes no share or recall from it. `totals.tsv` is counted from the genomes, not summed from the years. Species and antibiotics are named in full as the reports name them. No table holds a genome id, a catalogue group name or a row per genome.

**Guards**, before any table is written, counts only, and the read stops on any failure:

1. every report parses, and the reports read equal the reports found;
2. on every report with a catalogue row, the four cells of label by `catalog_resistant` equal `catalog.metrics` `tp`, `fp`, `tn` and `fn`; the reports that disagree are counted and named by species and antibiotic, and the read stops;
3. on every primary dataset, each of the five learners has a lineage-split row, its `scores` has the length of `samples`, and every score is a number from 0 to 1;
4. each stratum's genome, resistant and susceptible counts summed over the years, `year_unknown` included, equal its row in `totals.tsv`;
5. the datasets used plus those not used, by reason, equal the reports found;
6. no genome id and no catalogue group name stands in any written file, scanned on the box against every id and group name of the 274 reports; the scan prints counts only.

**What the page reads.** Descriptive figures only, fixed here. Per stratum: its lowest and highest share of a year's resistant genomes among the years that are not thin, with those years and their denominators; and the mean vote's lowest and highest recall among the years that are not thin where the stratum holds 10 or more resistant genomes. With four strata, point only and neither are read first, as the two the reader's acquired-gene catalogue cannot call by a gene of their own. With two groups, `reader_no_call` among resistant genomes is read first: the resistant genomes the reader misses. No test is run, no line is fitted, no slope is given and no cut is added after the tables are read. Recall is given at the release's cutoff; this read computes no log loss.

**Words.** A share describes the deposits of one public database. A stratum name is not a gene name; the page names no catalogue group, gene, variant or position and orders none. "No difference shown" is never written as "equal".

## Step 1

M1's published counts, the rows of `reference/mechanism-split/page-m1-strata.tsv` with set `primary`, split `lineage` and kind `stratum`, beside the rebuild from GRM's known-gene rule over the same 204 datasets, a genome once per dataset. The file's 7 rows of each stratum give the same counts.

| Stratum | Published resistant | Published susceptible | Rebuilt resistant | Rebuilt susceptible | Equal |
|---|---|---|---|---|---|
| point only | 7,963 | 6,791 | 0 | 0 | no |
| gene only | 39,814 | 11,041 | 48,749 | 20,864 | no |
| both | 9,066 | 1,479 | 0 | 0 | no |
| neither | 6,035 | 64,480 | 14,129 | 62,927 | no |
| all four | 62,878 | 83,791 | 62,878 | 83,791 | yes |

The rebuild puts no genome in point only or both; gene only holds more resistant and more susceptible genomes than M1 published, and neither more resistant and fewer susceptible, over the same 62,878 resistant and 83,791 susceptible genomes. M1 took its strata from the AMRFinderPlus program's own rows, which mark a point mutation. The reports hold GRM's known-gene rule, whose catalogue holds acquired genes only, and no field of a report marks a group as a point mutation. The rebuild was not tuned, and the four strata are not used again on this page.

Where M1 and the reports read the same thing, they agree: the 47,364 resistant and 18,678 susceptible genomes of `reader_call` equal, summed over M1's four strata, the genomes the `catalog` row of `page-m1-strata.tsv` calls resistant (124, 38,057, 8,649 and 534 resistant; 224, 9,468, 1,341 and 7,645 susceptible).

## Groups

Per group over the primary set, counted from the genomes (`totals.tsv`). Resistant called and susceptible called are the genomes the five learners' mean vote calls resistant; recall is resistant called over resistant.

| Group | Genomes | Resistant | Susceptible | Resistant called | Susceptible called | Recall |
|---|---|---|---|---|---|---|
| the reader calls resistant (`reader_call`) | 66,042 | 47,364 | 18,678 | 41,800 | 2,870 | 0.883 |
| the reader makes no call (`reader_no_call`) | 80,627 | 15,514 | 65,113 | 11,989 | 2,483 | 0.773 |
| the reader has no rule (`reader_no_rule`) | 0 | 0 | 0 | 0 | 0 |  |

## Years

Per year that is not thin, then `year_unknown`: the year's genomes, each group's count and share of the year's resistant and susceptible genomes, and the mean vote's recall within each group where it holds 10 or more resistant genomes that year. The two shares of a year add to 1, as `reader_no_rule` is empty. The 32 thin years keep their rows in `years.tsv` and `shares.tsv` and are not quoted. A genome counts once per dataset it appears in.

| Year | Genomes | Resistant | Susceptible | `reader_call` among resistant | `reader_no_call` among resistant | `reader_call` among susceptible | `reader_no_call` among susceptible | Recall, `reader_call` | Recall, `reader_no_call` |
|---|---|---|---|---|---|---|---|---|---|
| 1985 | 203 | 13 | 190 | 10 (0.769) | 3 (0.231) | 1 (0.005) | 189 (0.995) | 0.900 |  |
| 1994 | 149 | 51 | 98 | 35 (0.686) | 16 (0.314) | 0 (0.000) | 98 (1.000) | 0.514 | 0.500 |
| 1995 | 101 | 43 | 58 | 32 (0.744) | 11 (0.256) | 4 (0.069) | 54 (0.931) | 0.656 | 0.909 |
| 1996 | 257 | 145 | 112 | 101 (0.697) | 44 (0.303) | 6 (0.054) | 106 (0.946) | 0.723 | 0.864 |
| 1997 | 275 | 180 | 95 | 120 (0.667) | 60 (0.333) | 6 (0.063) | 89 (0.937) | 0.692 | 0.767 |
| 1998 | 225 | 81 | 144 | 48 (0.593) | 33 (0.407) | 7 (0.049) | 137 (0.951) | 0.562 | 0.697 |
| 1999 | 188 | 101 | 87 | 60 (0.594) | 41 (0.406) | 4 (0.046) | 83 (0.954) | 0.733 | 0.854 |
| 2000 | 199 | 118 | 81 | 72 (0.610) | 46 (0.390) | 15 (0.185) | 66 (0.815) | 0.847 | 0.848 |
| 2001 | 275 | 160 | 115 | 99 (0.619) | 61 (0.381) | 9 (0.078) | 106 (0.922) | 0.879 | 0.934 |
| 2002 | 224 | 143 | 81 | 95 (0.664) | 48 (0.336) | 5 (0.062) | 76 (0.938) | 0.863 | 0.833 |
| 2003 | 490 | 342 | 148 | 206 (0.602) | 136 (0.398) | 61 (0.412) | 87 (0.588) | 0.864 | 0.904 |
| 2004 | 881 | 599 | 282 | 330 (0.551) | 269 (0.449) | 70 (0.248) | 212 (0.752) | 0.967 | 0.955 |
| 2005 | 875 | 577 | 298 | 308 (0.534) | 269 (0.466) | 36 (0.121) | 262 (0.879) | 0.899 | 0.911 |
| 2006 | 1,475 | 933 | 542 | 679 (0.728) | 254 (0.272) | 99 (0.183) | 443 (0.817) | 0.951 | 0.858 |
| 2007 | 1,757 | 881 | 876 | 574 (0.652) | 307 (0.348) | 357 (0.408) | 519 (0.592) | 0.929 | 0.863 |
| 2008 | 3,130 | 1,505 | 1,625 | 1,197 (0.795) | 308 (0.205) | 601 (0.370) | 1,024 (0.630) | 0.882 | 0.818 |
| 2009 | 3,318 | 1,351 | 1,967 | 978 (0.724) | 373 (0.276) | 757 (0.385) | 1,210 (0.615) | 0.879 | 0.826 |
| 2010 | 3,849 | 1,637 | 2,212 | 1,199 (0.732) | 438 (0.268) | 661 (0.299) | 1,551 (0.701) | 0.892 | 0.858 |
| 2011 | 5,132 | 2,163 | 2,969 | 1,888 (0.873) | 275 (0.127) | 581 (0.196) | 2,388 (0.804) | 0.915 | 0.785 |
| 2012 | 8,853 | 4,410 | 4,443 | 3,806 (0.863) | 604 (0.137) | 989 (0.223) | 3,454 (0.777) | 0.875 | 0.637 |
| 2013 | 7,136 | 3,817 | 3,319 | 3,261 (0.854) | 556 (0.146) | 649 (0.196) | 2,670 (0.804) | 0.897 | 0.558 |
| 2014 | 10,535 | 5,696 | 4,839 | 4,873 (0.856) | 823 (0.144) | 1,420 (0.293) | 3,419 (0.707) | 0.874 | 0.685 |
| 2015 | 15,957 | 7,249 | 8,708 | 5,623 (0.776) | 1,626 (0.224) | 2,235 (0.257) | 6,473 (0.743) | 0.890 | 0.697 |
| 2016 | 12,015 | 4,522 | 7,493 | 3,545 (0.784) | 977 (0.216) | 1,598 (0.213) | 5,895 (0.787) | 0.902 | 0.811 |
| 2017 | 10,924 | 3,956 | 6,968 | 2,478 (0.626) | 1,478 (0.374) | 1,102 (0.158) | 5,866 (0.842) | 0.902 | 0.704 |
| 2018 | 4,716 | 1,954 | 2,762 | 1,498 (0.767) | 456 (0.233) | 571 (0.207) | 2,191 (0.793) | 0.930 | 0.702 |
| 2019 | 755 | 251 | 504 | 125 (0.498) | 126 (0.502) | 105 (0.208) | 399 (0.792) | 0.840 | 0.317 |
| 2020 | 1,456 | 275 | 1,181 | 202 (0.735) | 73 (0.265) | 320 (0.271) | 861 (0.729) | 0.851 | 0.849 |
| 2021 | 369 | 53 | 316 | 45 (0.849) | 8 (0.151) | 113 (0.358) | 203 (0.642) | 0.889 |  |
| year_unknown | 50,329 | 19,565 | 30,764 | 13,801 (0.705) | 5,764 (0.295) | 6,260 (0.203) | 24,504 (0.797) | 0.864 | 0.827 |

The lowest and highest of each column over the 29 years that are not thin, with the year and the denominator; a recall only where the group holds 10 or more resistant genomes that year. A tie goes to the earliest year.

| Group | Lowest share of resistant | Highest share of resistant | Lowest share of susceptible | Highest share of susceptible | Lowest recall | Highest recall |
|---|---|---|---|---|---|---|
| `reader_call` | 0.498 (2019, 125 of 251) | 0.873 (2011, 1,888 of 2,163) | 0.000 (1994, 0 of 98) | 0.412 (2003, 61 of 148) | 0.514 (1994, 18 of 35) | 0.967 (2004, 319 of 330) |
| `reader_no_call` | 0.127 (2011, 275 of 2,163) | 0.502 (2019, 126 of 251) | 0.588 (2003, 87 of 148) | 1.000 (1994, 98 of 98) | 0.317 (2019, 40 of 126) | 0.955 (2004, 257 of 269) |

## Species

`species.tsv` holds a row per species, antibiotic and group with 10 or more resistant genomes: 171 for `reader_call` and 118 for `reader_no_call`, over the 204 species and antibiotics of the primary set in 21 species. Its columns are those of `years.tsv` without the year and the thin flag, and its `year_unknown` column counts the row's genomes with no four-digit year. The page reads no figure from it: it shows which species and antibiotics stand behind the totals.

## Seen after

**One year at both ends.** The highest share of `reader_no_call` among a year's resistant genomes, 0.502 (2019, 126 of 251), and its lowest recall, 0.317 (2019, 40 of 126), fall in the same year. This was not looked into.

## Limits

**Two groups, not four.** The reader's call is not M1's strata. A resistant genome the reader makes no call on may carry a point mutation the AMRFinderPlus program reports, or nothing the catalogue holds; this read cannot tell them apart.

**Deposits.** The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species. A share is no prevalence, and a share over collection years describes deposits, not resistance in any population.

**Mixed years.** A year's genomes come from the datasets that hold that year, and the species and antibiotics deposited differ between years. A year's share mixes them, and no year table is split by species or antibiotic.

**Collection year.** The year is the one the deposit records. 50,329 genomes carry none, or none of four digits, and stand apart as `year_unknown`.

**One database.** A split by lineage is a retrospective split inside one database, never external validation.

**Recall only.** Recall is given at the release's cutoff of 0.5 and only over 10 or more resistant genomes. The plan fixed no log loss, so a recall here stands without the scores' calibration beside it. No rate among susceptible genomes is written, only their counts.

**Counting.** A genome counts once per dataset it appears in, so a genome tested for three antibiotics counts three times.

**Markers.** No catalogue group, gene, variant or position is named, ranked or ordered. The page holds aggregates over years, groups, species and antibiotics: no genome id, no BioProject accession and no row per genome. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

## Method

**Reading.** The read ran once, on 24 Nov 2026, at commit `3db9e14` with a clean tree. `read` printed counts only; its output, a file of counts, was scanned against every genome id and catalogue group name of the reports and held none. `write` checked the guards and wrote the record. `reading.txt` lists the commands with the sha256 of their inputs and outputs.

```
python3 -I scripts/mechanism/year_mix.py read --reports REPORTS --lists scripts/c3-inputs/lists.tsv --published reference/mechanism-split/page-m1-strata.tsv --out READING
python3 -I scripts/mechanism/year_mix.py write --reading READING --reports REPORTS --out OUT
python3 -I scripts/mechanism/year_mix.py print reference/mechanism-year
```

REPORTS is the folder of the release's reports, READING the read's counts and OUT the record's folder, copied to `reference/mechanism-year/` as `write` wrote it. The `print` command came after the read: it reads the record alone and prints this page's tables, rounded to 3 decimals. The read ran the script whose sha256 `reading.txt` gives; later commits added `print` and changed nothing that `read` or `write` runs.

**Guards.** `write` checked step 1 and the six guards before it wrote a table, and none stopped it. `guards.txt` holds its printout:

```
guard 1: 274 reports found, 274 read
guard 2: 253 reports with a catalogue row, their cells equal catalog.metrics; 21 with no catalogue row and no call
guard 3: 204 datasets, 1020 learner rows on the lineage split, each with one score per genome, each a number from 0 to 1
guard 4: 3 strata, each summed over 62 year rows (32 thin, year_unknown included where present) equals its total
guard 5: 204 used + 33 in_M + 31 not_in_A + 6 not_in_list + 0 no_lineage_split = 274 reports
step 1: 7 published rows the same; strata two: reader_call, reader_no_call, reader_no_rule
guard 6 (tables): 6 files, 0 with a genome id, 0 with a catalogue group name, against 23757 ids and 409 group names
```

**Tests.** `tests/mechanism/test_year_mix.py` holds 37 SYNTHETIC tests of the script, with tables worked by hand. `tests/mechanism/test_year_faults.py` runs `scripts/mechanism/year_faults.py`, which plants each fault below in the script in turn and runs the test meant to catch it. `tests/mechanism/test_year_page.py` holds 14 tests of this page against `reference/mechanism-year/`.

| Fault planted | Test that caught it |
|---|---|
| Y1. A mean of exactly 0.5 not called | `Rules.test_a_mean_of_exactly_half_is_called` |
| Y2. The five scores added one by one in order, so a mean of 0.5 can fall under it | `Rules.test_a_mean_of_exactly_half_is_called` |
| Y3. One learner's score read as the vote | `Read.test_the_call_reads_all_five_learners` |
| Y4. A year read from the first four characters | `Rules.test_a_year_not_four_digits_is_unknown` |
| Y5. Cells that differ from catalog.metrics let through | `Read.test_cells_that_disagree_with_metrics_stop_the_read` |
| Y6. A scores list of another length let through | `Read.test_a_scores_length_mismatch_stops` |
| Y7. Year counts that do not add up to the total let through | `Read.test_year_counts_that_do_not_add_up_stop_the_write` |
| Y8. A year of exactly 100 genomes flagged thin | `Read.test_a_year_under_100_genomes_is_thin` |
| Y9. The rebuild held to gene only and neither alone, so a mismatch reads as a match | `Read.test_a_rebuild_that_does_not_match_falls_back_to_two_groups` |
| Y10. A dataset of M read as primary | `Read.test_datasets_outside_the_primary_set_are_skipped_and_counted` |
| Y11. A stratum read from the label | `Read.test_the_rebuild_reads_groups_that_count_for_the_rule` |
| Y12. Every group counted for the rule | `Read.test_the_rebuild_reads_groups_that_count_for_the_rule` |
| Y13. A null call put with no call | `Read.test_a_null_call_is_its_own_group` |
| Y14. A recall written over under 10 resistant genomes | `Read.test_recall_only_over_ten_resistant` |
| Y15. The scan of the written files skipped | `Read.test_an_id_in_a_written_file_stops_the_write` |
| Y16. The species table's floor at 9 resistant genomes | `Read.test_the_species_table_keeps_strata_with_ten_or_more_resistant` |

All 16 were caught.

**Files.** The record is in `reference/mechanism-year/`, with `SHA256SUMS` over its six tables as `write` wrote them.

## Datasets

The 204 datasets of the primary set, species and antibiotic as the reports name them; each has a row in `species.tsv`. Not used: 33 reports in set M, 31 outside A and 6 not in `scripts/c3-inputs/lists.tsv`; no primary dataset lacked a lineage split.

| Species | Datasets | Antibiotics |
|---|---|---|
| *Acinetobacter baumannii* | 7 | cefepime, cefotaxime, doripenem, doxycycline, minocycline, moxifloxacin, sulfamethoxazole/trimethoprim |
| *Campylobacter coli* | 7 | azithromycin, ciprofloxacin, clindamycin, erythromycin, gentamicin, nalidixic acid, tetracycline |
| *Campylobacter jejuni* | 4 | ciprofloxacin, clindamycin, nalidixic acid, tetracycline |
| *Clostridioides difficile* | 4 | azithromycin, ceftriaxone, clarithromycin, moxifloxacin |
| *Corynebacterium diphtheriae* | 3 | sulfamethoxazole/trimethoprim, tetracycline, trimethoprim |
| *Enterobacter cloacae* | 9 | cefepime, ciprofloxacin, gentamicin, imipenem, levofloxacin, meropenem, sulfamethoxazole/trimethoprim, tobramycin, trimethoprim |
| *Enterococcus faecium* | 11 | ampicillin, clindamycin, dalfopristin/quinupristin, doxycycline, erythromycin, gentamicin, kanamycin, streptomycin, teicoplanin, tetracycline, vancomycin |
| *Escherichia coli* | 29 | amikacin, ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefotaxime, cefoxitin, cefpodoxime, cefpodoxime/clavulanate, ceftazidime, ceftiofur, ceftriaxone, cefuroxime, cephalothin, chloramphenicol, colistin, doxycycline, ertapenem, erythromycin, imipenem, levofloxacin, nalidixic acid, nitrofurantoin, piperacillin/tazobactam, streptomycin, sulfamethoxazole, tetracycline, tobramycin, trimethoprim |
| *Klebsiella aerogenes* | 3 | ertapenem, imipenem, meropenem |
| *Klebsiella oxytoca* | 4 | cefotaxime, cefoxitin, cefpodoxime, ertapenem |
| *Klebsiella pneumoniae* | 31 | ampicillin/sulbactam, aztreonam, cefazolin, cefepime, cefotetan, cefoxitin, cefpodoxime, cefpodoxime/clavulanate, ceftazidime, ceftazidime/avibactam, ceftriaxone, cefuroxime, chloramphenicol, ciprofloxacin, colistin, doripenem, doxycycline, ertapenem, fosfomycin, imipenem, levofloxacin, nitrofurantoin, ofloxacin, polymyxin b, streptomycin, sulfamethoxazole/trimethoprim, tetracycline, ticarcillin/clavulanate, tigecycline, tobramycin, trimethoprim |
| *Klebsiella quasipneumoniae* | 3 | cefotaxime, ceftazidime, sulfamethoxazole/trimethoprim |
| *Neisseria gonorrhoeae* | 6 | cefixime, ceftriaxone, ciprofloxacin, erythromycin, penicillin, tetracycline |
| *Pseudomonas aeruginosa* | 12 | amikacin, aztreonam, cefepime, ceftazidime, ceftazidime/avibactam, ceftolozane/tazobactam, gentamicin, imipenem, levofloxacin, meropenem, piperacillin/tazobactam, tobramycin |
| *Salmonella enterica* | 23 | amoxicillin/clavulanate, ampicillin, ampicillin/sulbactam, azithromycin, cefazolin, cefotaxime, cefoxitin, ceftazidime, ceftiofur, ceftriaxone, chloramphenicol, ciprofloxacin, gentamicin, kanamycin, nalidixic acid, spectinomycin, streptomycin, sulfamethoxazole, sulfamethoxazole/trimethoprim, sulfathiazole, sulfisoxazole, tetracycline, trimethoprim |
| *Shigella flexneri* | 4 | ampicillin, azithromycin, ciprofloxacin, trimethoprim |
| *Shigella sonnei* | 18 | amikacin, ampicillin, azithromycin, cefoperazone, cefotaxime, cefoxitin, ceftazidime, ceftriaxone, ciprofloxacin, colistin, ertapenem, gentamicin, nalidixic acid, streptomycin, sulfamethoxazole/trimethoprim, tetracycline, tobramycin, trimethoprim |
| *Staphylococcus aureus* | 12 | cefoxitin, chloramphenicol, ciprofloxacin, clindamycin, daptomycin, erythromycin, fusidic acid, gentamicin, mupirocin, penicillin, sulfamethoxazole/trimethoprim, tetracycline |
| *Streptococcus agalactiae* | 3 | clindamycin, erythromycin, tetracycline |
| *Streptococcus pneumoniae* | 8 | azithromycin, cefuroxime, chloramphenicol, clindamycin, erythromycin, meropenem, sulfamethoxazole/trimethoprim, tetracycline |
| *Vibrio cholerae* | 3 | nalidixic acid, sulfamethoxazole, trimethoprim |
