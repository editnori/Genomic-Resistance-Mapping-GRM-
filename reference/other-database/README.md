# The test on data from another public database: the record

This directory is the record behind [`docs/other-database-test.md`](../../docs/other-database-test.md). `scripts/other-database/record.py` writes every file of it from the saved files of lanes `outside` and `outside2`, the reports and interval files of release `rev14-274-20261112`, and the data root's drug name list. It fits nothing, runs no `grm` and downloads nothing. It recomputes every measure from the saved predictions and labels.

## Files

| File | Rows | Holds |
|---|---|---|
| `README.md` |  | this list |
| `SHA256SUMS` |  | the sha256 of every other file of this directory |
| `datasets.tsv` | 377 | one row per dataset and learner scored: labels, recall, specificity and log loss with intervals, the lineage figures, the plan's reading, the class mix |
| `differences.tsv` | 27 | derived: each pair of learners, per species and over the species, under the four readings |
| `genomes.tsv` | 13377 | the scored genomes as accessions only: no label and no prediction |
| `inputs.sha256` | 406 | the sha256 of every input file the script read |
| `learners.tsv` | 27 | derived: one row per species and learner, under the plan's rule and each reading, with the margin |
| `order.tsv` | 9 | per species, when predictions were written and when the label requests were sent; for the pilot, the requests before the plan with no recorded parameters |
| `plan/ADDENDUM.md` |  | the eight species' addendum, as lane outside2 fixed it |
| `plan/ADDENDUM2.md` |  | the second addendum with readings L, C and D, as lane outside2 fixed it |
| `plan/PLAN.md` |  | the pilot's plan, as lane outside fixed it |
| `plans.tsv` | 3 | the three plan files with their sha256 and the time each was fixed |
| `readings.tsv` | 246 | one row per estimable dataset and learner: readings L, C and D (and T for the pilot), and the errors in the largest BioProject |
| `settings.tsv` | 22 | the release, the class line, the pass rule, and the seed and draws of every bootstrap |
| `species.tsv` | 9 | one row per species: each funnel step to the scored genomes, the overlap counts (each a lower bound), and the largest BioProject |
| `totals.tsv` | 3 | derived: per learner, the counts over the eight species |

Every file is plain text. Each TSV has a header row and is sorted by its first columns. The script writes the same bytes on every run and at any hash seed.

## Rules

- **Labels.** An NCBI antibiotic name maps to a release drug by exact equality, then by the data root's `curated/drug_names.json`, then by the plan's fixed table. One label per BioSample and drug: R and S together is a conflict, then R, then S; intermediate alone gives no label.
- **Measures.** As `plan/PLAN.md` section 6 fixes them. A bootstrap draws the labelled genomes with a call in the order of their assembly accessions.
- **Readings.** L, C and D are those of `plan/ADDENDUM2.md` section B2, fixed in advance for the eight species and exploratory for the pilot. T is the pilot review's two-sample reading, exploratory, and the pilot's L, C and D are that review's readings too.
- **Ties.** A tie on the largest BioProject goes to the lower accession: letters first, then the lower number, so PRJNA717739 comes before PRJNA1160997. A field that holds several BioProjects, or a genome with none, comes after every single accession.
- **Words.** `meets` in `learners.tsv` and `totals.tsv` means the learner meets the plan's rule in the species: it holds on more than half of at least 3 estimable datasets. `margin` is how many datasets would have to change for the other result.
- **Pre-plan requests.** `unrecorded_before_plan` in `order.tsv` counts the pilot's requests before `PLAN.md` was fixed that the request log holds with no parameters.

## Commands

```sh
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --out reference/other-database
python3 -I scripts/other-database/record.py --lanes "$LANES" --release "$R" --data-root "$ROOT" --check reference/other-database
```

`LANES` holds the lanes' output directories `outside/` and `outside2/`, `R` is the unpacked release and `ROOT` is the public data root. `--check` exits 0 only when this directory is what the script writes from them.
