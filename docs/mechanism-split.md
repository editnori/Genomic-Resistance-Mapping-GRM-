# Mechanism split: recall by mechanism and the mechanism mix over collection years

Plan blocks M1 and M2 read the AMRFinderPlus program's tables beside the release's reports. They fit nothing. M1 asks whether the learners call resistant genomes whose catalogue calls are point mutations as often as those whose calls are acquired genes. M2 asks whether the mix of catalogue calls among deposited genomes moves between collection years. M3 is not built. M4 was read once on 24 Nov 2026, after TBProfiler 6.6.3 ran again on the 2,000 genomes, in `docs/lineage-curve-result.md`. The design was fixed on 23 Nov 2026 before any script of this block read a real table, score, stratum or share. The scripts are in `scripts/mechanism/` and the tests in `tests/mechanism/`. Blocks M1 and M2 were read once, after the second run, in `docs/mechanism-split-result.md`. Beside the dry run on Klebsiella aerogenes for shape and time, the first run on the 21 species of A began at 06:39Z on 23 Nov 2026 and stopped three species on two defects of the scripts, since fixed (`docs/assumption-plan-corrections.md`, Mechanism split, first run).

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

**Stratum.** Over a genome's counted rows that meet the dataset's rule: point only (a `POINT` or `POINT_DISRUPT` row, no `AMR` row), gene only (an `AMR` row, no `POINT` or `POINT_DISRUPT` row), both, neither. Beside each stratum: its genomes whose rule-meeting rows are all partial calls, and its genomes with a `POINT_DISRUPT` row. Beside point only and both: the genomes whose point rows are `POINT_DISRUPT` alone. The stratum comes from the table and the rule, never from the label. The lead ruled at 04:13Z on 23 Nov 2026, before any result was read, that a `POINT_DISRUPT` row is read with point (question 1); the plan took it as gene.

**Tables.** Per dataset (species and antibiotic in full), split and row (the five learners, GRM's known-gene rule, the program): resistant genomes and those called resistant, susceptible genomes and those called resistant, per stratum. A count is written as it is, 0 to 4 included (question 12, ruled at 05:30Z on 23 Nov 2026, which withdrew the ruling on question 2 that wrote 1 to 4 as "under 5"): a row per genome is barred, an aggregate keeps its genome count, and a mask that the cells beside it undo would claim a protection it does not give. A recall is written only over 10 or more resistant genomes, and no rate among susceptible genomes is written, only their counts. Strata under 10 pool per species and drug class.

**Reading.** Per learner, on the lineage split of the primary set: recall on gene only less recall on point only, per dataset where both strata hold 10 or more resistant genomes. The median with a species interval (4,000 resamples of species, seeded by the row's name), gene ahead or point ahead by more than 0.02, the two-sided sign test, and Holm over the five learners. The random split is beside, labelled "training and test share lineages".

## M2, the mechanism mix over collection years

**Genomes.** Every genome of the species with a table that appears in a release report. Its year is `samples[].year`: four digits are a year, anything else is no year, and nothing is imputed. Periods: 2010 to 2014 against 2015 to 2019 (primary), 2000 to 2009 against 2010 to 2019 (secondary).

**Cell.** A species with a drug class as the program's Class column names it. Per period: the share of genomes with a `POINT` or `POINT_DISRUPT` row of the class, with the count of genomes with a `POINT_DISRUPT` row beside it, the share with a gene row (Subtype `AMR`), and among genomes with a gene row the mix of catalogue families (`gene_family` of the database's `ReferenceGeneCatalog.txt`, through `ReferenceGeneHierarchy.txt` when needed) over their Subtype `AMR` rows alone. A `POINT_DISRUPT` row enters no family mix (questions 1 and 10, ruled at 04:13Z on 23 Nov 2026). The program names some hits at a node of `ReferenceGeneHierarchy.txt` above the gene families, such as `msr` and `lsa`, whose chains run through `abc-f` to `AMR`. Such a Subtype `AMR` row, whose symbol is a node with no `gene_family` at or above it, is counted under that node's own name, and a symbol in neither file still stops the run (question 14, ruled at 06:56Z on 23 Nov 2026). The floor on names holds for such a name as for any family, and the step's last log line counts its rows and genomes.

**One vote per BioProject.** The BioProject is BV-BRC's `bioproject_accession` from the genome table the study and time splits read. A BioProject votes in a period with 5 or more genomes in it, and a cell passes with 3 or more voting BioProjects in each period. A permutation of period labels over BioProjects (1,000, seeded by the cell's name) gives the p of the point share's difference and of the family mix's total variation distance. Holm runs over the passing cells, each statistic apart; a cell moves at a Holm p under 0.05. A voted row that does not pass writes no share, mix or p, since with one or two voting BioProjects in a period a share would be one BioProject's own, and the share per genome is written only where the cell's BioProject vote passes. A family is named only where, in each period, 5 or more genomes of the voting BioProjects carry it or none do. The rest share one row, `other families`. Each family row gives its genome count per period as it is.

**Guards, beside and deciding nothing.** The share per genome with no vote; the largest BioProject left out; BioProjects with 0.95 or more resistant labels for the class, or with one family in 0.95 or more of their gene genomes, left out; one vote per lineage cluster, given only where the species' reports agree on the clusters; the same vote with only the clusters that vote in both periods (in the run's script at main `11c2996`, named here after the reading, on 23 Nov 2026); countries present in both periods only; rows of no partial method only. The flag reads a BioProject over its genomes in every year, since it describes the study and a study's design does not change by period (question 13, ruled at 05:30Z on 23 Nov 2026). Where the reports disagree, the count written is of genome pairs, each with a cluster in both reports, that share a cluster in one and not in the other, each pair once. Resistant genomes per dataset are a secondary reading.

**Time-split link, exploratory.** For the 24 datasets of `docs/time-split.md`: of the resistant later test genomes of the past arm with a catalogue call, the share with a call no resistant training genome carries (a `POINT` or `POINT_DISRUPT` row is one call, "point"; worded after the reading, on 23 Nov 2026, as the run's script at main `11c2996` computes it), against each learner's past-less-random difference in recall and in log loss, by Spearman's rank correlation.

## M3 and M4

M3 would set a point share from the feature cards beside C8's families. The cards do not allow it: a card names a gene in 1,432 of the release's 21,443 final-model cards, in none of the cards of 6 of the 23 species with cards, by BV-BRC's annotation rather than the catalogue's name, and never in an rRNA gene. Nothing is built.

M4 is the lineage curve for Mycobacterium tuberculosis. `tb_curve.py` reads TBProfiler 6.6.3 run again on the 2,000 genomes with `tbdb` commit `3ed1cb99` (WHO catalogue v2), as the TBProfiler comparison ran it; that run's per-genome output was not kept. A determinant of a drug is a variant the tool reports for it, passing the tool's filter, graded Assoc w R or Assoc w R - Interim (`who_assoc`), and a report's antibiotic meets the tool's drug by the comparison's name map. Over the resistant genomes of each antibiotic in full, grouped by the tool's main lineage, by its sublineage and by the report's lineage cluster:

- the count of distinct determinants as groups are added, over 1,000 random orders seeded by the antibiotic and the grouping (median, 2.5 and 97.5 percentiles), with the order by group size, largest first, beside;
- the share of resistant genomes with no graded determinant, per group and over all;
- how many determinants are seen in 1, 2, 3 and 4 or more lineages, and the resistant genomes each count covers;
- per gene as the WHO catalogue names it, the count of its graded variants seen.

A gene seen in 1 to 4 resistant genomes of an antibiotic is not named: such genes share one row, `other genes`. A share is written over 5 or more resistant genomes and no curve over under 5, and every count, the summary's too, is written as it is. A variant's own name enters no output until the owner rules on it. A dataset in M is marked. The tool ran again on 24 Nov 2026 with a result for 2,000 of 2,000 genomes, and M4 was read once, with no dry run, in `docs/lineage-curve-result.md`.

## What a table can and cannot say

A table says which catalogue calls the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species, so a share is no prevalence and a trend describes deposits. A split by lineage or study is a retrospective split inside one database, never external validation. No output holds a genome id, a BioProject accession, a row per genome, a sequence or a position.

## Running it

One species on a box with no agent attached, from a clean checkout:

```
setsid nohup bash scripts/mechanism/run_species.sh KIT SPECIES OUT META TABLES DB FOLDS > /dev/null 2>&1 < /dev/null &
```

KIT holds the release, `amrkit.tar` and `tool-tables.sha256`, each pinned. TABLES is a folder of the program's tables as `TABLES/SPECIES/<genome>.tsv`, each held to its recorded sha256; with `-` the program runs on FASTA. DB is a copy of the database 2026-08-07.1; with `-` the kit's `fetch_database.sh` fetches it from NCBI. Both are held to the kit's lists. When TABLES or DB is `-`, a step of its own installs the program's environment from the kit (micromamba pinned by sha256, `spec/env-amrfinder.explicit.txt`, `amrfinder --version` 4.2.7) before the database step, since the fetch needs the program's `amrfinder_index` (the first run's defect A, 23 Nov 2026). With TABLES and DB both given, nothing is installed. The steps run once each, a restart skips finished ones, and a scan fails the run when any file under `results/` or `logs/` names a genome of the species or a BioProject accession, holds a sequence or a position (a point symbol of the program or a change as TBProfiler writes it; a bare number is not read as one), or when a table writes a rate under its floor or a name carried by under 5 genomes: a recall over under 10 resistant genomes, an M2 share or mix in a cell whose BioProject vote does not pass, a time-split or M4 share over under 5 genomes, a family or gene named where under 5 genomes carry it, a count cell that is not a count, or a table the scan does not know. A count of 0 to 4 is not a hit. Over species, `combine.py --run OUT ...` reads only runs with `DONE`.

M4 runs on the machine that holds the TBProfiler run, `RUN_DIR` with its `plan.txt`, `failed.tsv` and `results/<genome>.results.json`:

```
python3 -I scripts/mechanism/tb_curve.py --species mycobacterium_tuberculosis --reports REPORTS \
  --lists scripts/c3-inputs/lists.tsv --tool-dir RUN_DIR --bed TBDB.bed --out OUT
python3 -I scripts/mechanism/scan.py --ids SPECIES.json OUT
```

It stops on a result file of another tool or database version and on a planned genome with neither a result nor a line in `failed.tsv`. Only `OUT` leaves the machine, after the scan.

## Checks

`tests/mechanism` has 253 tests: 210 SYNTHETIC tests, 18 of them on `tb_curve.py`, 18 on the result page's table script and 1 on its planted faults, 20 on the lineage curve's table script and 1 on its planted faults, and 35 on the catalogue size reading's script (`docs/catalogue-size-result.md`) and 1 on its planted faults; 19 of the result page against its record; 12 of the lineage curve's page against its record; and 12 of the catalogue size page against its record. `planted_faults.py` plants thirteen faults in turn, each caught by its test: ten from a stratum read from the label to a family named other than as the database names it, and three from the first run, the database fetched before the program is installed, a node above the gene families given its first child's family, and a symbol in neither database file let through. A review planted 28 more, in the scripts and in `run_species.sh`, from a score of exactly 0.5 not called resistant to an M2 guard that keeps the BioProjects it flags. Of the 28, 14 passed every test at first; each now fails a test. Under the rulings of 05:30Z two of them no longer apply (a count of 1 to 4 written as a figure, and a share beside "under 5"); a ninth M2 guard fault and eight more on the rulings' rules, from a share over under 5 genomes let through the scan to a lineage pair counted once per pair of reports, each fail a test.

The dry run on Klebsiella aerogenes (93 genomes, 3 datasets) ran end to end in 2.5 s with 30.6 MB at peak. The hard check passed on 3 of 3 datasets (36 of 36 rows), and the scan found no genome id in 28 files. Run again after the review's patches, the hard check passed on the same 36 rows, and the scan, which now also reads each table's counts, found no hit in 28 files, 9 of them tables. Run again at the rulings of 05:30Z, both gave the same counts. From the first run and a timing on 2,000 SYNTHETIC genomes, the 21 species need about 10 minutes of compute where the tables are given, and about 265 CPU hours where the program must run again.
