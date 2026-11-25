# The mechanism year read of 24 November 2026: the record

This directory holds the tables behind [docs/mechanism-year-result.md](../../docs/mechanism-year-result.md): the mix of catalogue strata by collection year over the 204 datasets of the primary set, read once from the 274 revision 14 reports. Every row is an aggregate over a year, a stratum, a species and antibiotic, or the whole set, and keeps its genome count. No file holds a genome id, a catalogue group name or a row per genome.

`scripts/mechanism/year_mix.py write` wrote the six tables and `SHA256SUMS`; they are copied here as it wrote them. `reading.txt` gives the commands with the sha256 of their inputs and outputs.

| File | Holds |
|---|---|
| `step1.tsv` | M1's eight published counts beside the rebuild from GRM's known-gene rule, per stratum, and the totals |
| `years.tsv` | per year and group: genomes, resistant, susceptible, those the five learners' mean vote calls resistant, recall over 10 or more resistant genomes, and the thin flag |
| `totals.tsv` | per group, counted from the genomes |
| `shares.tsv` | per year and group: the group's share of the year's resistant and susceptible genomes, with the year's denominators and the thin flag |
| `species.tsv` | per species, antibiotic and group with 10 or more resistant genomes, named in full |
| `guards.txt` | the six guards and their counts |
| `no-call-strata.tsv` | added after the reading: per M1 stratum, the resistant genomes, those GRM's known-gene rule calls, those it makes no call on and the share with no call, then the three strata with a program row and all four |

The groups are `reader_call`, `reader_no_call` and `reader_no_rule`, from each genome's call by GRM's known-gene rule, because the four strata of M1 were not readable from the reports (`step1.tsv`). A genome counts once per dataset it appears in. A year with fewer than 100 genomes is thin. A share describes the deposits of one public database: it is no prevalence and no trend. These are public bacterial genomes with public accession numbers, and the record holds no field about a patient.

`no-call-strata.tsv` was added on 25 Nov 2026, after the one reading, by `scripts/mechanism/year_mix.py strata --m1 reference/mechanism-split/page-m1-datasets.tsv --lists scripts/c3-inputs/lists.tsv --record reference/mechanism-year`. Its only new input is M1's committed public table `page-m1-datasets.tsv`: the `catalog` row on the lineage split of the primary set, summed over the 204 datasets. It opens no release report and no score and changes none of the six tables, which `SHA256SUMS` still covers as `write` wrote them. It was added because the program's no-row count and the reader's no-call count needed relating, not listing. `strata` checks that the row holds each primary dataset once and that its sums equal `step1.tsv` and `totals.tsv` before it writes.

`scripts/mechanism/year_mix.py print reference/mechanism-year` prints the page's tables from these files.
