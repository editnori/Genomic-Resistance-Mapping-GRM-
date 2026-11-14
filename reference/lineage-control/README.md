# Within-lineage test of each rule term

Each table here holds one release's set covering machine rule terms, each with the result of a test of whether its association with the lab calls holds inside lineage clusters. The dataset inspector attaches the verdict to the caveat under the rules. `scripts/lineage-control/verdicts.py` wrote the table, and `scripts/lineage-control/site_data.py` writes the app's copy from it.

## Input

| File | Release asset | Asset sha256 | Evaluation revision | Reports read |
| --- | --- | --- | --- | --- |
| `rev14-274-20261112.tsv` | `rev14-274-20261112.tar.zst` | `c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de` | 14 | 274 reports |

The column `report_sha256` gives the hash of each report file. On 2026-11-14 each of the 274 equalled its line in the release's own checksum list, `rev14-274-20261112.sha256`. `rev14-274-20261112.reports.sha256` holds those 274 lines as the release wrote them, and `tests/lineage_control` compares each row's `report_sha256` with its line.

## Command

Run from the repository root, with `R` the unpacked release directory (the one that holds `datasets/`). Python 3.12, standard library only. The script reads reports and fits nothing.

```sh
python3 -I scripts/lineage-control/verdicts.py "$R" > reference/lineage-control/rev14-274-20261112.tsv
python3 -I scripts/lineage-control/site_data.py reference/lineage-control/rev14-274-20261112.tsv > web/src/data/lineage-control.json
```

The table was written on 2026-11-14 at the commit that added it. The second file is the release's checksum list cut down to its reports, with `L` the list:

```sh
grep -F '  ./datasets/' "$L" > reference/lineage-control/rev14-274-20261112.reports.sha256
```

## Output

The sha256 of each file, in the form `sha256sum -c` reads from this directory. `tests/lineage_control` checks each line against its file.

```
de62860a2a559d5606411d26c74bbb245d466e8224a5c1b1a0576f47ff294589  rev14-274-20261112.tsv
0aedf8780ca8e2292cfc5f227f61b87b7c6b006a6649ed4c0378514c70b72223  rev14-274-20261112.reports.sha256
```

## The test

The test was fixed on 2026-11-13, before any result was read, and the script runs it without change.

- **Context.** A term is tested on the isolates the earlier terms of its own model left undecided, the set the greedy search chose it on. In a conjunction those are the isolates every earlier term holds for; in a disjunction, the isolates no earlier term holds for.
- **Strata.** The report's lineage clusters (`samples[].lineage`). A cluster is informative when it holds both states of the term and both lab calls.
- **Statistic.** Cochran-Mantel-Haenszel over the clusters, with the continuity correction, one degree of freedom.
- **Assessable.** At least 5 informative clusters.
- **Verdict.** Benjamini-Hochberg over the assessable terms of the release at 0.05.

## Columns

One row for each rule term: 862 rule terms of 274 reports.

| Column | Meaning |
| --- | --- |
| `dataset`, `revision` | The report's id and its evaluation revision |
| `position` | The term's place in its model, from 1 |
| `feature`, `presence` | The rule: the feature's number, and 1 for presence or 0 for absence |
| `input_fingerprint` | The report's `input_fingerprint` |
| `report_sha256` | The sha256 of the report file |
| `isolates`, `context_isolates` | The report's isolates, and those the term is tested on |
| `informative_strata` | Lineage clusters that hold both states of the term and both lab calls |
| `excess` | Resistant isolates where the term holds, over what the clusters' margins give |
| `p`, `q` | The two-sided p, and its adjusted value; `q` is empty for a term that cannot be tested |
| `verdict` | See below |

| Verdict | Terms | Leading terms | Meaning |
| --- | --- | --- | --- |
| `seen` | 327 | 179 | Adjusted p under 0.05, and more resistant isolates than expected where the term holds |
| `not_shown` | 6 | 1 | Assessable, and adjusted p not under 0.05 |
| `cannot_test` | 529 | 94 | Fewer than 5 informative clusters |
| `reversed` | 0 | 0 | Adjusted p under 0.05, and fewer resistant isolates than expected |

No term of 90 datasets can be tested.

## What the site shows, and when

`web/src/data/lineage-control.json` holds, for each dataset, the report's input fingerprint and each term's feature, presence flag and verdict. The inspector shows a verdict only on a report with the same dataset id, evaluation revision and input fingerprint whose rule at that position is the same. A later release therefore shows no verdict until its own table is written here. The site has no words for `reversed` and says nothing of such a term.

## Limits

- A verdict is an association. `seen` says the association is not explained by the cut into lineage clusters. It says nothing about cause, and annotation of the k-mers does not establish mechanism.
- The model and the order of its terms were chosen on all isolates, so the test is in-sample with respect to that choice. It is not held-out evidence.
- The clusters are MinHash clusters at one threshold per dataset. They are not a phylogeny.
- `cannot_test` is not a negative result. It says these genomes hold too few clusters with both lab calls and both states of the term.
- A later term is tested only on the isolates the earlier terms left undecided, so its verdict does not describe the feature over the whole dataset.
