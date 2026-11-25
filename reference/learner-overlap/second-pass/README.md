# The learner overlap's second pass of 25 November 2026: the record

This directory holds the second pass behind [docs/learner-overlap-result.md](../../../docs/learner-overlap-result.md). The first pass counted the five learners one by one, and the release ships their mean. So the second pass reads the same 274 revision 14 release reports again and counts the resistant genomes the five-learner mean misses at the release's 0.5 cutoff and the susceptible genomes it calls resistant, each split by how many of the five miss or call it. It also works the first pass's yardstick again with the genomes of null-fit folds left out. Every row is an aggregate over a report, a species or the pooled reports, and keeps its genome count. No file holds a genome id. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

`scripts/mechanism/overlap_second.py read` read the reports. It ran every first-pass check again and stopped unless each report's first-pass columns, counted again, equalled the committed `../datasets.tsv`. `write` then wrote `vote.tsv`, `vote-pooled.tsv` and `gates.txt`, with `SHA256SUMS` over the three. They are copied here as it wrote them. `reading.txt` lists the commands with the sha256 of what each wrote. The run's directories are named REPORTS, FIRST, READING and RECORD there in place of their paths. Nothing was fitted, trained or scored again: the mean of five stored scores is arithmetic.

| File | Holds |
|---|---|
| `vote.tsv` | one row per report the first pass used, species and antibiotic in full: the resistant genomes the mean misses and the susceptible genomes it calls resistant, each by all five, 1 to 4 or none of the five; the resistant genomes outside null-fit folds, each learner's misses there and the all-five misses there, with the yardstick and the ratio worked from them |
| `vote-pooled.tsv` | the same counts summed over the first pass's pooled reports (`scope` `all`) and per species, with the median share the mean misses and the median and quartiles of the ratio outside null-fit folds, each report weighted once |
| `gates.txt` | the first pass's checks run again, the second pass's six checks, and the first pass's reading rule applied to the median ratio outside null-fit folds |
| `reading.txt` | the second pass's commands and the sha256 of each one's output |

The yardstick is the count the five would miss together if they missed independently. It is neither a model nor a test, and no p-value is computed. A split by lineage is a retrospective split inside one database, never external validation. `scripts/mechanism/overlap_second.py print` prints the page's second-pass tables from these files and `../`.
