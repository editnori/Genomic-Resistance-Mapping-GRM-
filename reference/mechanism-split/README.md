# The mechanism split's reading of 23 November 2026: the record

This directory holds the tables behind [docs/mechanism-split-result.md](../../docs/mechanism-split-result.md): blocks M1 (recall by mechanism) and M2 (the mechanism mix over collection years) of the mechanism split, read once over the runs of the 21 species of A, 237 datasets, 204 of them in the primary set and 33 in M. Every row is an aggregate over a dataset, a species, a cell or a learner. No file holds a genome id, a BioProject accession, a sequence or a position.

`scripts/mechanism/combine.py` wrote the five files under `combined/`; they are copied here as it wrote them. `scripts/mechanism/page_tables.py write` checked the reading against the six gates of `gates.txt` and wrote every other file from the runs' own tables, species and antibiotic in full, with no new measure. `reading.txt` lists the reading's commands with the sha256 of what each wrote, and `checks.txt` the checks of the run, counts only. The rules of the reading are `docs/mechanism-split.md`'s, fixed before any run, with the rulings in `docs/assumption-plan-corrections.md` (entries "Mechanism design" and "Mechanism split, first run").

| File | Holds |
|---|---|
| `gates.txt` | the six gates the reading passed before a table was written |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the run, counts only |
| `combined/m1-reading.tsv` | per split, set and row: recall on gene only less recall on point only over the datasets where both hold 10 or more resistant genomes; the lineage split of the primary set decides for the five learners |
| `combined/m1-neither.tsv` | per split, set and row: the median recall on the neither stratum |
| `combined/m2-reading.tsv` | per M2 row and statistic: the value, its permutation p, Holm over the passing cells and whether it moves |
| `combined/m2-summary.tsv` | per cell kind, period pair, guard and statistic: cells, passing, moving |
| `combined/timesplit-link.tsv` | per learner: Spearman's rank correlation of the unseen-determinant share with the past-less-random difference (exploratory) |
| `page-species.tsv` | per species, named in full: datasets, checks, M1 datasets with both strata at 10 or more, M2 genomes, cells, passing and moving |
| `page-datasets.tsv` | per dataset of A, species and antibiotic in full: the checks and whether it entered |
| `page-m1-datasets.tsv` | per entered dataset, split and row: resistant genomes, called and recall per stratum, and gene only less point only |
| `page-m1-strata.tsv` | per set, split, row and stratum: genome counts summed over the entered datasets, a genome once per dataset; no pooled recall |
| `page-m1-pooled.tsv` | strata under 10 resistant genomes pooled per species, drug class and set, as `strata.py` wrote them |
| `page-m2-cells.tsv` | every M2 row of the runs, species in full, with combine.py's Holm p and moves |
| `page-m2-families.tsv` | the family mix of each passing cell, the catalogue's family names over 5 or more genomes, the rest as `other families` |
| `page-m2-years.tsv` | per species, genomes per span of collection years |
| `page-m2-lineage.tsv` | per species, whether the reports agree on the lineage clusters |
| `page-timesplit.tsv` | per time-split dataset: resistant later test genomes with a determinant no resistant training genome carries |

A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A trend describes deposits. A table names no cause and no new marker. A split by lineage or study is a retrospective split inside one database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/page_tables.py print` prints the page's tables from these files.
