# Mechanism split: recall by mechanism and the mechanism mix over collection years

Plan blocks M1 and M2 read the AMRFinderPlus program's tables beside the release's reports. They fit nothing. M1 asks whether the learners call resistant genomes whose catalogue calls are point mutations as often as those whose calls are acquired genes. M2 asks whether the mix of catalogue calls among deposited genomes moves between collection years. M3 is not built, and M4 has a script that has not run on a real file. The design was fixed on 23 Nov 2026 before any script of this block read a real table, score, stratum or share. The scripts are in `scripts/mechanism/` and the tests in `tests/mechanism/`. No result of either block has been read: the one real run so far is a dry run on Klebsiella aerogenes for shape and time.

## Three callers, kept apart

| Caller | What it is | Where its call is |
|---|---|---|
| The AMRFinderPlus program | AMRFinderPlus 4.2.7, database 2026-08-07.1, `--plus` with the organism option | the program's table per genome, read by the comparison's rule (`tool_default`) |
| GRM's known-gene rule | GRM's own reader of the AMRFinderPlus catalogue, `src/catalog.rs` | each report's `catalog` row and `samples[].catalog_resistant` |
| The five learners | scm, cart, rf, gbt, lasso | each split's `learners[].scores`, resistant at 0.5 or more |

GRM's known-gene rule is never called AMRFinderPlus here.

## The comparison's rule, reused

The AMRFinderPlus comparison (`docs/amrfinderplus-comparison.md`) turned a table into a call with each report's own published rule, `catalog.rule`. This block reuses that rule and writes no new map.

- `scripts/mechanism/drug_rules.tsv` is the rule of each of the 237 datasets in A, written by `drug_rules.py` from the release's reports: species, antibiotic, the report's drug class and the rule in the catalogue's own class and subclass names. Every run holds each report's rule equal to its row.
- `amr_rule.py` holds the comparison adapter's table reader and variants and its `token_match` and `explain`, copied with their source sha256. `run_species.sh` runs the kit's own adapter unchanged, and `strata.py` stops unless its call equals the adapter's on every genome.
- Rows counted: Type AMR, Subtype `AMR`, `POINT` or `POINT_DISRUPT`, Scope core (`tool_default`). Subtype `AMR-SUSCEPTIBLE` counts for none, and an unknown subtype stops the run.
- On each dataset the program's tp, fn, tn and fp must equal its row `program:tool_default` in `reference/amrfinderplus-comparison/callers.tsv`.

## M1, recall by mechanism

**Datasets.** A, the 237 datasets of the 21 species with tables (`scripts/c3-inputs/lists.tsv`). The primary set is A not in M, 204 datasets. The 33 in M are reported apart with the label rule each depends on.

**Hard check.** Before any stratum is counted, for each dataset, split and row with metrics (the five learners and the catalog row), tp, fn, tn and fp recomputed from `scores` at 0.5 must equal the report's `metrics`, and the catalog row's calls must equal `samples[].catalog_resistant`. A dataset that fails is named and left out.

**Stratum.** Over a genome's counted rows that meet the dataset's rule: point only (a `POINT` row, no other), gene only (an `AMR` or `POINT_DISRUPT` row, no `POINT`), both, neither. Beside them: genomes whose rule-meeting rows are all partial calls, and genomes with a `POINT_DISRUPT` row. The stratum comes from the table and the rule, never from the label.

**Tables.** Per dataset (species and antibiotic in full), split and row (the five learners, GRM's known-gene rule, the program): resistant genomes and those called resistant, susceptible genomes and those called resistant, per stratum. A count under 5, zero included, is written "under 5". A recall is written only over 10 or more resistant genomes. Strata under 10 pool per species and drug class.

**Reading.** Per learner, on the lineage split of the primary set: recall on gene only less recall on point only, per dataset where both strata hold 10 or more resistant genomes. The median with a species interval (4,000 resamples of species, seeded by the row's name), gene ahead or point ahead by more than 0.02, the two-sided sign test, and Holm over the five learners. The random split is beside, labelled "training and test share lineages".

## M2, the mechanism mix over collection years

**Genomes.** Every genome of the species with a table that appears in a release report. Its year is `samples[].year`: four digits are a year, anything else is no year, and nothing is imputed. Periods: 2010 to 2014 against 2015 to 2019 (primary), 2000 to 2009 against 2010 to 2019 (secondary).

**Cell.** A species with a drug class as the program's Class column names it. Per period: the share of genomes with a `POINT` row of the class, the share with a gene row, and among genomes with a gene row the mix of catalogue families (`gene_family` of the database's `ReferenceGeneCatalog.txt`, through `ReferenceGeneHierarchy.txt` when needed).

**One vote per BioProject.** The BioProject is BV-BRC's `bioproject_accession` from the genome table the study and time splits read. A BioProject votes in a period with 5 or more genomes in it, and a cell passes with 3 or more voting BioProjects in each period. A permutation of period labels over BioProjects (1,000, seeded by the cell's name) gives the p of the point share's difference and of the family mix's total variation distance. Holm runs over the passing cells, each statistic apart; a cell moves at a Holm p under 0.05.

**Guards, beside and deciding nothing.** The share per genome with no vote; the largest BioProject left out; BioProjects with 0.95 or more resistant labels for the class, or with one family in 0.95 or more of their gene genomes, left out; one vote per lineage cluster, given only where the species' reports agree on the clusters; countries present in both periods only; rows of no partial method only. Resistant genomes per dataset are a secondary reading.

**Time-split link, exploratory.** For the 24 datasets of `docs/time-split.md`: the share of resistant later test genomes of the past arm with a catalogue call no resistant training genome carries, against each learner's past-less-random difference in recall and in log loss, by Spearman's rank correlation.

## M3 and M4

M3 would set a point share from the feature cards beside C8's families. The cards do not allow it: a card names a gene in 1,432 of the release's 21,443 final-model cards, in none of the cards of 6 of the 23 species with cards, by BV-BRC's annotation rather than the catalogue's name, and never in an rRNA gene. Nothing is built.

M4 is the lineage curve for Mycobacterium tuberculosis. `tb_curve.py` reads TBProfiler 6.6.3 run again on the 2,000 genomes with `tbdb` commit `3ed1cb99` (WHO catalogue v2), as the TBProfiler comparison ran it; that run's per-genome output was not kept. A determinant of a drug is a variant the tool reports for it, passing the tool's filter, graded Assoc w R or Assoc w R - Interim (`who_assoc`), and a report's antibiotic meets the tool's drug by the comparison's name map. Over the resistant genomes of each antibiotic in full, grouped by the tool's main lineage, by its sublineage and by the report's lineage cluster:

- the count of distinct determinants as groups are added, over 1,000 random orders seeded by the antibiotic and the grouping (median, 2.5 and 97.5 percentiles), with the order by group size, largest first, beside;
- the share of resistant genomes with no graded determinant, per group and over all;
- how many determinants are seen in 1, 2, 3 and 4 or more lineages, and the resistant genomes each count covers;
- per gene as the WHO catalogue names it, the count of its graded variants seen.

A variant's own name enters no output until the owner rules on it. A dataset in M is marked. No run has been made: the tool must run again first, and the brief gives M4 no dry run.

## What a table can and cannot say

A table says which catalogue calls the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species, so a share is no prevalence and a trend describes deposits. A split by lineage or study is a retrospective split inside one database, never external validation. No output holds a genome id, a BioProject accession, a row per genome, a sequence or a position.

## Running it

One species on a box with no agent attached, from a clean checkout:

```
setsid nohup bash scripts/mechanism/run_species.sh KIT SPECIES OUT META TABLES DB FOLDS > /dev/null 2>&1 < /dev/null &
```

KIT holds the release, `amrkit.tar` and `tool-tables.sha256`, each pinned. TABLES is a folder of the program's tables as `TABLES/SPECIES/<genome>.tsv`, each held to its recorded sha256; with `-` the program is installed from the kit and run on FASTA. The steps run once each, a restart skips finished ones, and a scan fails the run when any file under `results/` or `logs/` names a genome of the species or a BioProject accession. Over species, `combine.py --run OUT ...` reads only runs with `DONE`.

M4 runs on the machine that holds the TBProfiler run, `RUN_DIR` with its `plan.txt`, `failed.tsv` and `results/<genome>.results.json`:

```
python3 -I scripts/mechanism/tb_curve.py --species mycobacterium_tuberculosis --reports REPORTS \
  --lists scripts/c3-inputs/lists.tsv --tool-dir RUN_DIR --bed TBDB.bed --out OUT
python3 -I scripts/mechanism/scan.py --ids SPECIES.json OUT
```

It stops on a result file of another tool or database version and on a planned genome with neither a result nor a line in `failed.tsv`. Only `OUT` leaves the machine, after the scan.

## Checks

`tests/mechanism` has 82 SYNTHETIC tests, 14 of them on `tb_curve.py`. `planted_faults.py` plants ten faults in turn, from a stratum read from the label to a family named other than as the database names it, and each is caught by its test.

The dry run on Klebsiella aerogenes (93 genomes, 3 datasets) ran end to end in 2.5 s with 30.6 MB at peak. The hard check passed on 3 of 3 datasets (36 of 36 rows), and the scan found no genome id in 28 files. From it and a timing on 2,000 SYNTHETIC genomes, the 21 species need about 10 minutes of compute where the tables are given, and about 265 CPU hours where the program must run again.
