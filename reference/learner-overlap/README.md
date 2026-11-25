# The learner overlap's reading of 24 November 2026: the record

This directory holds the tables behind [docs/learner-overlap-result.md](../../docs/learner-overlap-result.md): on the lineage split of the 274 revision 14 release reports, how many of the five learners miss each resistant genome at the release's 0.5 cutoff, beside the count expected if they missed independently. Every row is an aggregate over a report, a species or the pooled reports, and keeps its genome count. No file holds a genome id. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

`scripts/mechanism/learner_overlap.py read` read the reports and `write` wrote `datasets.tsv`, `pooled.tsv` and `gates.txt`, with `SHA256SUMS` over the three. They are copied here as it wrote them. `reading.txt` lists the reading's commands with the sha256 of what each wrote; the run's directories are named REPORTS, READING and RECORD there in place of their paths. Nothing was fitted, trained or scored again.

| File | Holds |
|---|---|
| `datasets.tsv` | one row per report used, species and antibiotic in full: resistant genomes missed by exactly 0 to 5 of the five learners, susceptible genomes called resistant by exactly 0 to 5, each learner's misses, the all-five misses with no catalogue call or in a fold with a null fit, the yardstick and the ratio |
| `pooled.tsv` | the same counts summed over the pooled reports (`scope` `all`) and per species, with the median and quartiles of the ratio, each report weighted once |
| `gates.txt` | the ten checks the reading passed before a table was written, and the reading rule |
| `reading.txt` | the reading's commands and the sha256 of each one's output |

The yardstick is the count the five would miss together if they missed independently. It is neither a model nor a test, and no p-value is computed. A split by lineage is a retrospective split inside one database, never external validation. `scripts/mechanism/learner_overlap.py print` prints the page's tables from these files.
