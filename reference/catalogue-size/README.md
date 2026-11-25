# The catalogue size reading: the record

This directory holds the tables behind [docs/catalogue-size-result.md](../../docs/catalogue-size-result.md): whether the number of graded variants of the WHO catalogue among an antibiotic's resistant genomes tracks the resistant share of its genomes or only how many resistant genomes there were to study, for the antibiotics of Mycobacterium tuberculosis. It was read once, by the rule the page's section "Rule" fixed before any count was read. Every row is an aggregate over an antibiotic or over a set of antibiotics. No file holds a genome id, a BioProject accession, a variant's own name, a gene's name, a sequence or a position.

`scripts/mechanism/catalogue_size.py read` took three counts from each of the species' release reports, `counts.samples`, `counts.resistant` and `counts.susceptible`, and nothing else. `scripts/mechanism/catalogue_size.py write` joined them to the lineage curve's `reference/lineage-curve/page-endpoints.tsv` (page `docs/lineage-curve-result.md`) and wrote the two tables. `in_M`, `kind`, `resistant_with_result`, `catalogue_variants`, `genes_named`, `without_determinant` and `share_without` are carried from that file unchanged.

| File | Holds |
|---|---|
| `antibiotics.tsv` | one row per report of the species, species and antibiotic in full: its genomes, resistant and susceptible genomes and resistant share from the release report, the lineage curve's columns carried unchanged, set M marked, and whether the antibiotic enters the correlations and, if not, why |
| `correlations.tsv` | one row per set and pair: the Spearman rank correlation of `catalogue_variants` against the resistant share, the resistant genomes or the genomes, its n, the smallest and largest value with one antibiotic dropped and which one, whether the row decides, what the rule's comparison gives over the set, and the set's antibiotics |
| `guards.txt` | the six guards the reading passed, counts only |
| `reading.txt` | the reading's commands and the sha256 of each one's output |

A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A correlation over antibiotics of one species, one catalogue and one database names no cause. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/catalogue_size.py print` prints the page's tables from these files.
