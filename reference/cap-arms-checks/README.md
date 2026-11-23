# The cap-arms reading: two counts beside it

These two tables stand beside [reference/cap-arms/](../cap-arms/README.md), the tables behind [docs/cap-arms-result.md](../../docs/cap-arms-result.md). Neither decides anything. Every row is a count over a species, and no file holds a genome id, a BioProject, a sequence or a feature.

| File | Holds |
|---|---|
| `ties.tsv` | per species, named in full: the genomes the release's selection rule ranks on the rebuilt revision 14 table, and the genomes that share the rule's score at the 4,096th place (`tied`), those of them inside the top 4,096 (`kept`) and those past it (`left`); a species whose rule ranks 4,096 genomes or fewer has no cut. Written by `scripts/c2-cap/ties.py` |
| `added-genome-ids.tsv` | per species, named in full: the genome ids counted in the files each species' run wrote, against the release's ids and the ids of the genomes the wider cohort adds (`ids_known`); the ids the cohort lists add beyond the release (`ids_beyond_release`); and the ids found (`ids_found`). The last row holds the 8 species together. The count ran where each species ran, since the cohort lists stay there |

`SHA256SUMS` gives the sha256 of every other file here.
