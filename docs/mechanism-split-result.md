# Mechanism split result

This page reads blocks M1 and M2 of the mechanism split once, by the rules `docs/mechanism-split.md` fixed before any run. M1 asks whether each learner calls resistant genomes whose catalogue calls are point mutations as often as those whose calls are acquired genes. M2 asks whether the mix of catalogue calls among deposited genomes moves between collection years. Neither block fits anything: both read the AMRFinderPlus program's tables beside the release's reports.

The rules of the reading are the method page's sections "M1, recall by mechanism", "M2, the mechanism mix over collection years" and "What a table can and cannot say", with the rulings on questions 1 to 14 in `docs/assumption-plan-corrections.md` (entries "Mechanism design" and "Mechanism split, first run"). No measure, subset or cut was added. The reading ran once, with no other option.

Three callers stand apart on this page. The AMRFinderPlus program is the program itself, AMRFinderPlus 4.2.7 with database 2026-08-07.1. GRM's known-gene rule is GRM's own reader of the AMRFinderPlus catalogue, the catalog row of each report. The five learners are `scm`, `cart`, `rf`, `gbt` and `lasso`.

Every split here is a retrospective split inside one public database. A split by lineage or by study is never external validation. A share is no prevalence, and a trend describes deposits. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Reading

<!-- READING: the answer, after the reading -->

## Data

<!-- DATA: datasets, genomes and checks, after the reading -->

## M1

M1 reads each learner's recall on the resistant genomes of each stratum. A genome's stratum comes from the program's counted rows that meet the dataset's own rule, never from its label: point only (a `POINT` or `POINT_DISRUPT` row and no `AMR` row), gene only (an `AMR` row and no point row), both, and neither. A `POINT_DISRUPT` row is read with point (question 1).

The reading is, per learner on the lineage split of the primary set, recall on gene only less recall on point only, per dataset where both strata hold 10 or more resistant genomes. The median is given with a species interval from 4,000 resamples of species, read at the 101st and 3,900th values. A dataset is gene ahead or point ahead by more than 0.02. The two-sided sign test of gene ahead against point ahead runs with Holm over the five learners. A learner reads "gene only higher" or "point only higher" at a Holm p under 0.05, and "no difference shown" otherwise.

GRM's known-gene rule and the AMRFinderPlus program stand beside the learners and decide nothing. The random split stands beside, labelled "training and test share lineages", and the datasets in M stand apart. M1 gives recall alone, as the method fixes it (question 11).

<!-- M1: tables, after the reading -->

## M2

A cell is a species with a drug class as the program's Class column names it. Per period, the cell gives the share of genomes with a point row of the class, the share with a gene row (Subtype `AMR`), and among genomes with a gene row the mix of catalogue families over their Subtype `AMR` rows alone. The primary periods are 2010 to 2014 against 2015 to 2019, and the secondary periods 2000 to 2009 against 2010 to 2019. A year is four digits in `samples[].year`, and nothing is imputed.

One vote per BioProject: a BioProject votes in a period with 5 or more genomes in it, and a cell passes with 3 or more voting BioProjects in each period. A permutation of period labels over BioProjects, 1,000 times, gives the p of the point share's difference and of the family mix's total variation distance. Holm runs over the passing cells, each statistic apart, and a cell moves at a Holm p under 0.05. A cell that does not pass writes no share, mix or p. A family is named only where, in each period, 5 or more genomes of the voting BioProjects carry it or none do; the rest share one row, `other families`.

Only the BioProject vote on the primary periods over species and class cells decides. Everything else in this section and the next stands beside.

<!-- M2: tables, after the reading -->

## Guards

The guards stand beside the BioProject vote and decide nothing: the share per genome with no vote, the largest BioProject left out, flagged BioProjects left out (0.95 or more resistant labels for the class, or one family in 0.95 or more of their gene genomes, read over every year), one vote per lineage cluster where the species' reports agree on the clusters, the clusters voting in both periods, countries present in both periods only, and rows of no partial method only. Resistant genomes per dataset are a secondary reading. The time-split link is exploratory.

<!-- GUARDS: tables, after the reading -->

## Seen after

These points were seen in the tables after the reading. The rules did not ask for them. They change no reading above and decide nothing.

<!-- SEEN AFTER: after the reading -->

## Limits

**What a table can say.** A table says which catalogue calls the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. A gene family is named only as the AMRFinderPlus database names it, over 5 or more genomes.

**Deposits.** The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species. A share is no prevalence, and a trend over collection years describes deposits, not resistance in any population.

**One database.** A split by lineage or by study is a retrospective split inside one database, never external validation.

**Labels.** The datasets in M wait for the maintainer's label rule. They are marked and never join the primary reading.

**M4.** The lineage curve for Mycobacterium tuberculosis has not run. TBProfiler 6.6.3 must run again first on the 2,000 genomes.

**Markers.** No marker is found, named or ranked. The page holds aggregates over datasets, species, cells and learners only: no genome id, no BioProject accession, no sequence and no position.

## Method

<!-- METHOD: checks of the run, the reading's commands, gates, tests, after the reading -->

## Datasets

<!-- DATASETS: per species and per dataset, after the reading -->
