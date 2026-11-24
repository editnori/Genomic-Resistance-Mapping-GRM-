# Lineage curve result

This page reads block M4 of the mechanism split once, by the rules `docs/mechanism-split.md` fixed before any run (section "M3 and M4"). Among the resistant genomes of each antibiotic of *Mycobacterium tuberculosis*, it counts the distinct graded variants of the WHO catalogue that appear, how that count grows as lineages are added, how many resistant genomes carry no graded variant, and how many variants are seen in one lineage only. It fits nothing: it reads TBProfiler's output beside the release's reports.

"The WHO catalogue" here is the catalogue as TBProfiler 6.6.3 ships it, `tbdb` commit `3ed1cb99`, whose graded entries carry the source label "WHO catalogue v2". A graded variant of an antibiotic is a variant the tool reports for the antibiotic's drug, passing the tool's own filter, with an entry graded "Assoc w R" or "Assoc w R - Interim". A report's antibiotic meets the tool's drug by the TBProfiler comparison's name map: rifampin is rifampicin, and ciprofloxacin, ofloxacin, rifabutin and prothionamide read another drug's calls by a shared name stem (proxy).

A curve describes the deposits of one public database. A share is no prevalence. A split by lineage is a retrospective split inside one database, never external validation. A table names no cause and no new marker. Nothing here is a change comparison: no default, no released file and no sentence on the site changes.

## Reading

(Filled from `reference/lineage-curve/` after the one reading.)

## Data

**Genomes.** TBProfiler 6.6.3 ran again on the genomes of the species' revision 14 reports, from their assemblies (FASTA), one genome per run with one thread, as the TBProfiler comparison ran it on 14 Nov 2026. That run's per-genome output was not kept, so the tool ran again.

**Reports.** A resistant genome is one the report's label calls resistant. Each report's genomes are its own; a genome can be resistant to several antibiotics and is counted once in each.

**Lineages.** Three groupings stand side by side: the tool's main lineage, its sublineage, and the report's own lineage cluster. A lineage is named as TBProfiler prints it. A genome the tool gives no lineage enters no group of that grouping and is counted in the Datasets section. A mixed call is its own group, as the tool writes it.

**Set M.** The reports of clofazimine and isoniazid are in set M: their labels wait for the maintainer's label rule. They are marked in every table and read with the release's labels, as the method page fixed.

## Curves

The curve counts the distinct graded variants among the resistant genomes as lineage groups are added one at a time. Over 1,000 random orders of the groups, seeded by the antibiotic and the grouping, each step gives the median and the 2.5 and 97.5 percentiles (the band). Beside it, the order by group size, largest first. No curve is drawn over under 5 resistant genomes.

## No graded variant

A resistant genome with no graded variant carries no variant of the antibiotic's drug that the catalogue grades. It is a count. Nobody looked for what explains these genomes. A share is written over 5 or more resistant genomes.

## Spread over lineages

Each graded variant is counted by the number of lineages it is seen in: 1, 2, 3, or 4 or more. Beside each count, the resistant genomes that carry at least one variant of that count.

## Genes

Per gene as the WHO catalogue names it, the count of its graded variants seen among the resistant genomes, and the resistant genomes with one. A gene seen in 1 to 4 resistant genomes of an antibiotic is not named: such genes share one row, `other genes`. A variant's own name enters no table.

## Seen after

These points were seen in the tables after the reading. The rules did not ask for them. They change nothing above and decide nothing.

## Limits

**What a table can say.** A table says which graded variants of the catalogue the resistant genomes carry. It does not say what causes resistance in a genome, and it names no new marker. A gene is named only as the WHO catalogue names it, in 5 or more resistant genomes of the antibiotic.

**Deposits.** The genomes are BV-BRC deposits chosen for having lab results, capped at 2,000 a species. A share is no prevalence, and a curve describes the deposits of one public database, not the variants of any population.

**One database.** A split by lineage is a retrospective split inside one database, never external validation.

**Assemblies.** The tool is built for sequencing reads, and the project holds assemblies. An assembly is one consensus sequence and cannot show a mixed population, and a target region the assembly does not cover in full gives no variant there. On these genomes the TBProfiler comparison (`docs/tbprofiler-comparison.md`) found 819 of 21,737 genome-drug rows with no call where an assembly did not cover a target region in full.

**The catalogue's version.** The grades are those of `tbdb` `3ed1cb99`. Another version of the catalogue can grade other variants. The catalogue was built from public isolates with recorded phenotypes, and some of these genomes may be among them; the overlap was not checked.

**No graded variant.** A resistant genome with no graded variant is counted, and nothing more is said of it.

**Labels.** The reports in set M wait for the maintainer's label rule. They are marked.

**Markers.** No marker is found, named or ranked. The page holds aggregates over antibiotics, lineage groups and genes only: no genome id, no BioProject accession, no variant name, no sequence and no position.

## Method

**Tool run.** On one machine, from the TBProfiler comparison's own scripts: the environment from its explicit list of 191 packages, held to their hashes, with micromamba 2.9.0-0 held to its published sha256; the annotation file the package ships copied to the place the tool reads, as that comparison did; the reference genome as a control; then `run_tbprofiler.sh` over the cohort list `cohort.py` builds from the reports.

**Checks of the run.** Before any table was read, `scripts/mechanism/curve_tables.py check` read the run's folder and printed counts only. `checks.txt` in `reference/lineage-curve/` holds its printout.

**Reading.** The reading ran once, with the command of `docs/mechanism-split.md`, then the scan of its output, then the table script, then the scan of the record. `reading.txt` lists the commands with the sha256 of each output.

**Tables.** `scripts/mechanism/curve_tables.py write` checked the reading against seven gates and wrote `reference/lineage-curve/`.

**Tests.** `tests/mechanism/test_curve_tables.py` holds SYNTHETIC tests of the table script, and `scripts/mechanism/curve_faults.py` plants each fault below in it and runs the test meant to catch it.

**Files.** The reading's tables are in `reference/lineage-curve/`. `scripts/mechanism/curve_tables.py print` prints this page's tables from those files, rounded to 3 decimals.

## Datasets
