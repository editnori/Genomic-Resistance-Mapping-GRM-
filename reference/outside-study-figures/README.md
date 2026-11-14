# Outside-study report figures

Each file here holds the lineage-split figures of one later release, for the datasets of the outside-study test. `docs/outside-study-results.md`, section Revision 14, reads the test's saved estimates against them. `scripts/outside-study-results/revision.py` wrote the file and reads it.

## Input

| File | Release asset | Asset sha256 | Evaluation revision | Reports read |
| --- | --- | --- | --- | --- |
| `rev14-274-20261112.tsv` | `rev14-274-20261112.tar.zst` | `c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de` | 14 | 97 reports |

The 97 reports are those of the datasets in `reference/outside-study`. The column `report_sha256` gives the hash of each report file. On 2026-11-14 each of the 97 equalled its line in the release's own checksum list, `rev14-274-20261112.sha256`.

## Command

Run from the repository root, with `R` the unpacked release directory. Python 3.12, standard library only.

```sh
python3 -I scripts/outside-study-results/revision.py reference/outside-study --figures "$R/datasets" > reference/outside-study-figures/rev14-274-20261112.tsv
```

The file was written on 2026-11-14 at the commit that added it.

## Output

The sha256 of the file, in the form `sha256sum -c` reads from this directory. `tests/outside_study_results` checks the line against the file.

```
7bb75d536159b1e303e778311c9c95ae74ee9de65eb0a5f0a164b89b08bda0ab  rev14-274-20261112.tsv
```

## Columns

One row for each dataset and learner of the runs' `metrics.tsv` files, 317 rows.

- `dataset`, `learner`: as in `metrics.tsv`. The known-gene rule is the report's `catalog` row.
- `report_sha256`, `evaluation_revision`: the report file's hash and its `evaluation.revision`.
- `lineage_n`, `lineage_tp`, `lineage_fn`, `lineage_fp`, `lineage_tn`, `lineage_sensitivity`, `lineage_specificity`: the `metrics` of the report's lineage-split row, as the kit reads them (`score.reference_row`).
- `lineage_note`: the kit's note where the report has no such row. The seven figure columns are then empty.
- `called_n`, `called_tp`, `called_fn`, `called_fp`, `called_tn`: the row's `called_metrics`, where the report has them. From revision 14 a row with a null fit keeps there its figure over the fitted folds only (`docs/equal-scoring.md`).
