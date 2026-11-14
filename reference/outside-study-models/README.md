# Outside-study final models

Each file here holds one digest of each final model in the reports of one release, for the datasets of the outside-study test. `docs/outside-study-results.md`, section Revision 14, says the runs' estimates are estimates for models the site serves now. That holds only if the reports of the served revision hold the models the runs scored, and these two files are the check. `scripts/outside-study-results/models.py` wrote them and compares them.

## Input

| File | Release asset | Asset sha256 | Evaluation revision | Reports read |
| --- | --- | --- | --- | --- |
| `rev12-274-20261105.tsv` | `grm-release-rev12-274-20261105.tar.gz` | `cefaef3957f8c9e161b32e1989cd4561c9a314430fa7eeb75b97eefebbbbf5c4` | 12 | 97 reports |
| `rev14-274-20261112.tsv` | `rev14-274-20261112.tar.zst` | `c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de` | 14 | 97 reports |

The 97 reports are those of the datasets in `reference/outside-study`. The column `report_sha256` gives the hash of each report file. `tests/outside_study_results` ties both files to their reports with hashes the repository already holds:

- Revision 12: each hash is the one the run's own record names for the report it read (`results.json.gz`, key `reports`), and the kit's pin for it in `scripts/outside-study/reports.sha256`.
- Revision 14: each hash is the report's line in `reference/outside-study-figures/rev14-274-20261112.reports.sha256`, the release's own checksum lines, and the hash in the figure file there.

## Command

Run from the repository root, with `R12` and `R14` the unpacked release directories (each holds `datasets/`). Python 3.12, standard library only.

```sh
python3 -I scripts/outside-study-results/models.py reference/outside-study --digests "$R12/datasets" > reference/outside-study-models/rev12-274-20261105.tsv
python3 -I scripts/outside-study-results/models.py reference/outside-study --digests "$R14/datasets" > reference/outside-study-models/rev14-274-20261112.tsv
```

The files were written on 2026-11-14 at the commit that added them. The script stops on a report of another dataset, of another schema than the kit's, or with an evaluation marked invalid.

To compare the two files, and to read a file again from its reports:

```sh
python3 -I scripts/outside-study-results/models.py reference/outside-study-models/rev12-274-20261105.tsv reference/outside-study-models/rev14-274-20261112.tsv
python3 -I scripts/outside-study-results/models.py reference/outside-study-models/rev14-274-20261112.tsv --reports "$R14/datasets"
```

The first exits 1 if a model's digest or a dataset's input fingerprint differs between the two files, and the suite runs it. The second exits 1 unless every row is the one the reports give. It needs the reports, so it runs where a release is unpacked and not in CI. Both files passed it on 2026-11-14.

If a digest differs, this prints each field that differs between the two reports of that dataset:

```sh
python3 -I scripts/outside-study-results/models.py --fields "$R12/datasets/<dataset>.json" "$R14/datasets/<dataset>.json"
```

## Output

The sha256 of each file, in the form `sha256sum -c` reads from this directory. `tests/outside_study_results` checks each line against its file.

```
ff50f8eb0717e5621f15de09e13402d3cd069dc05aae0a51aafde20fb8d0e21e  rev12-274-20261105.tsv
12b5bae878d7db038f5d4427fac74cf7d3f6c32095551b0d48dc1df518fb96ad  rev14-274-20261112.tsv
```

## Result

On 485 of the 485 pairs of a dataset and a learner the two files have the same `model_sha256`, and on 97 of the 97 datasets the known-gene rule has. No final model of these datasets differs between revision 12 and revision 14.

The two reports of a dataset have one `input_fingerprint` on 97 of the 97 datasets.

The runs scored 227 of these models: every set covering machine and decision tree, and the lasso where a report lists every feature it reads. For each of them, the features in both files are the ones the run wrote to `features.tsv`.

## Columns

Six rows for each dataset: the five learners (`scm`, `cart`, `rf`, `gbt`, `lasso`) and the known-gene rule (`known_gene`), 582 rows.

- `dataset`, `learner`.
- `report_sha256`, `evaluation_revision`: the report file's hash and its `evaluation.revision`.
- `input_fingerprint`: the report's own hash of what its evaluation read apart from code (`input_fingerprint`, `src/analysis.rs`): the species matrix, the settings, and each genome with its label and lineage.
- `model_sha256`: the sha256 of one JSON text with three keys. `learner` is the row's learner. `model` is everything the report holds for it under `models`. `features` is each feature id those fields name, in id order, with the `kmer`, `n_kmers` and `bits` of its entry in the report's `features` list, or with three nulls where the report does not list it. The text has its keys sorted and no spaces, and each number as Python prints it, so two reports give the same digest only if every field is equal to the last digit.
- `read_features_sha256`: the sha256 of the features the kit reads for a model (`score.learner_features`), each as `[feature, kmer, n_kmers]`, in id order. Empty for `rf`, `gbt` and `known_gene`, which have no such list.

What `model` holds for each row:

| Row | Fields of the report |
| --- | --- |
| `scm` | `models.scm`: the kind, the rules, the steps with their equivalent rules, the importances and the selected setting |
| `cart` | `models.cart`, every node, and `models.cart_depth` |
| `lasso` | `models.lasso`: the intercept, every nonzero weight and the penalty. Also `models.lasso_fit` |
| `rf` | `models.rf_importance` |
| `gbt` | `models.gbt_importance` |
| `known_gene` | `catalog.rule` and `catalog.version` |

## Limits

- A report publishes no tree of the forest or the boosted trees, only their largest importances (`docs/prediction-contract.md`). An equal digest for `rf` or `gbt` says those importances are equal. It cannot say the trees are.
- A feature the report does not list is in the digest as its id alone, with no k-mer. The id is a row of the dataset matrix. The report's `input_fingerprint` covers that matrix, so an id stands for the same row in two reports with one fingerprint. The check goes no further than the fingerprint: it reads no matrix.
- The digest does not cover a feature's annotation (`locus`, `catalog_families`, `catalog_products`), the folds or any cross-validation figure. The figure file in `reference/outside-study-figures` covers the lineage-split figures.
- The files cover the 97 datasets of the test, and not the other reports of a release.
