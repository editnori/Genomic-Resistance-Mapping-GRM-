# Mechanism year result

This page reads, once, how the mix of catalogue strata among the release's genomes changes with the year each genome was collected, and how the five learners' mean vote calls the resistant genomes of each stratum in each year. It fits nothing, trains nothing and scores no genome again: every figure comes from the held-out scores, labels and catalogue calls stored in the 274 revision 14 reports. The maintainer named this read on 23 Nov 2026 as one of three that need no fit.

A share over collection years describes the deposits of one public database. It is no prevalence and no trend in resistance, in the world or in any population. No line is fitted and no slope is given. A split by lineage is a retrospective split inside one database, never external validation. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

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
