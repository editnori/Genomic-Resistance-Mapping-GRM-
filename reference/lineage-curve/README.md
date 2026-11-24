# The lineage curve's reading of 24 November 2026: the record

This directory holds the tables behind [docs/lineage-curve-result.md](../../docs/lineage-curve-result.md): block M4 of the mechanism split, the lineage curve for Mycobacterium tuberculosis, read once over TBProfiler 6.6.3 with `tbdb` commit 3ed1cb99 run again on the 2,000 genomes of the species' 22 revision 14 reports. Every row is an aggregate over an antibiotic, a lineage group or a gene. No file holds a genome id, a BioProject accession, a variant's own name, a sequence or a position.

`scripts/mechanism/tb_curve.py` wrote the five files under `m4/`; they are copied here as it wrote them. `scripts/mechanism/curve_tables.py write` checked the reading against the seven gates of `gates.txt` and wrote `page-endpoints.tsv` from those files and `scripts/c3-inputs/lists.tsv`, species and antibiotic in full, with no new measure. `reading.txt` lists the reading's commands with the sha256 of what each wrote, and `checks.txt` the checks of the tool run, counts only. The rules of the reading are `docs/mechanism-split.md`'s, section "M3 and M4", fixed before any run.

| File | Holds |
|---|---|
| `gates.txt` | the seven gates the reading passed before a table was written |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the tool run, counts only |
| `m4/m4-curve.tsv` | per antibiotic, grouping and step: distinct graded variants as lineage groups are added, the median and the 2.5 and 97.5 percentiles over 1,000 seeded random orders, and the count in the order by group size, largest first |
| `m4/m4-coverage.tsv` | per antibiotic, grouping and group, and over all: resistant genomes, those with no graded variant, and the share where 5 or more |
| `m4/m4-convergence.tsv` | per antibiotic and grouping: graded variants seen in 1, 2, 3 and 4 or more lineages, and the resistant genomes each count covers |
| `m4/m4-genes.tsv` | per antibiotic and gene, as the WHO catalogue names it: graded variants seen and resistant genomes with one; genes in 1 to 4 resistant genomes share the row `other genes` |
| `m4/m4-summary.json` | the tool, the database, the genomes planned, with a result and failed, and per report its kind, tool drugs, set M mark and the resistant genomes with no result or no group |
| `page-endpoints.tsv` | one row per report, species and antibiotic in full, set M marked: the counts the page reads |

"The WHO catalogue" here is the catalogue as TBProfiler 6.6.3 ships it (`tbdb` 3ed1cb99). A curve describes the deposits of one public database. A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A split by lineage is a retrospective split inside one database, never external validation. A table names no cause and no new marker.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/curve_tables.py print` prints the page's tables from these files.
