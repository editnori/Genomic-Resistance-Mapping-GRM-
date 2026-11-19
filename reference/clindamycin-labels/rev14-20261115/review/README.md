# The review of the clindamycin label check: its files of 15 November 2026

This folder holds the files of the independent review of the check one folder up, as far as [docs/clindamycin-labels.md](../../../../docs/clindamycin-labels.md) quotes them. The review ran from 16:29Z to 16:52Z on 15 November 2026 (UTC). It recounted the six datasets by study with its own code, scored the check's saved reports again, and fitted `staphylococcus_aureus__clindamycin` at 11 fold plans per arm. Every count, refit score and rule verdict of the check reproduced. The review narrowed five parts of the reading, and the page now carries them: the shares over samples with a PMID, the kind of each dataset, the fold plans, the 33 distinct rows, and the history.

Every file here but this README is copied byte for byte from the review's output, at the same path. `../SHA256SUMS` lists each of them under `review/`. To check them, run from the folder above:

```sh
sha256sum -c SHA256SUMS
```

## Files

| File | What it is | Written by |
|---|---|---|
| `ADDENDUM.md` | The review's plan for the fold-plan runs, exploratory, written before any of them | by hand |
| `ADDENDUM.sha256` | Its sha256, `dccb2cccb6b8313c2432d0b4b16330e88ff6a60973e44b2c39d62da06fc37b3e`, recorded 2026-11-15T16:38:12Z, before any run | `sha256sum ADDENDUM.md` |
| `results/counts/by-key.tsv` | One row per dataset and study key (a sample's sorted PMIDs, or `no-pmid`), with a total per dataset, under two study fields: `label` (the PMIDs of the records that carry the sample's label, the plan's field) and `all` | `python3 -I scripts/recount.py` |
| `results/counts/chosen.tsv` | The PMID the plan's rule picks in each dataset, under each study field and three ways of counting a sample with two PMIDs. The last two columns give the study's share of the disagreements and of the samples among the samples that carry a PMID | the same command |
| `results/foldplans/check.txt` | The `c1f36fc` build's default seed against the check's reports of both arms, the seed and sample count of each run, and the number of distinct fold plans per arm and split | `python3 -I scripts/foldplans.py` |
| `results/foldplans/spread.tsv` | Per split, row, measure and view, the value at each of the 11 seeds with its minimum, maximum and range. V1 is arm "with" on all samples, V2 the same fit on the samples outside the study, V3 arm "without" | the same command |
| `results/foldplans/tests.tsv` | For each row that passes the plan's rule in balanced accuracy or recall, the ranges of V1, V2 and V3, the change at the same seed, and the addendum's two tests | the same command |
| `results/foldplans/summary.txt` | Per split, how many of those rows meet each test | the same command |
| `results/foldplans/specificity-beside.txt` | The addendum's first test for the passing specificity rows | a command the review did not save among its scripts (file time 16:47:15Z, after `foldplans.py`) |
| `results/refit/rescore-<dataset>/table.tsv` | For the three refit datasets: each split, row, measure and part of the change (total, composition, training), with V1, V2, V3, the change and its lineage-cluster interval under bootstrap seeds 7 and 31337 | `python3 -I scripts/rescore.py --dataset D --with W.json --without WO.json --removal IDS --seed 7 --seed 31337 [--placebo IDS=REPORT ...] --out DIR` |
| `results/refit/rescore-<dataset>/rule.tsv` | The plan's rule row by row: the total change and whether it passes under each seed | the same command |
| `results/refit/rescore-campylobacter_coli/placebo.tsv` | The total and training change of each of the 10 placebo arms, per split, row and measure | the same command |
| `results/refit/placebo-check.txt` | The 480 placebo values against the check's, and for each passing row of `campylobacter_coli__clindamycin` its change against its 10 placebo changes, with its rank among the 11 | a command the review did not save among its scripts (file time 16:37:13Z) |
| `results/exploratory/october-match.txt` | Exploratory, decided after the results were seen: the search of revision 14 views of `staphylococcus_aureus__clindamycin` for the October record's two figures, over 11 pairs of runs (seeds 42 and 1 to 10) | `python3 -I scripts/october_match.py <removal ids> <with.json>=<without.json> ...` |
| `scripts/recount.py` | The recount by study from the released reports, the curated phenotypes and every raw BV-BRC record | |
| `scripts/rescore.py` | The scores of one dataset's arms, from saved reports | |
| `scripts/foldplans.py` | The fold-plan spread and tests of the addendum | |
| `scripts/october_match.py` | The October search | |

The figures of those two files can be rebuilt from files of this folder: `placebo.tsv` and `rule.tsv` of `campylobacter_coli__clindamycin`, and `spread.tsv`.

The scripts name the paths of the machine they ran on (`RELEASE` and `ROOT` in `recount.py`; `RUNS`, `LANE` and `REM` in `foldplans.py`), and they read files that are not in the repository: the release, a copy of the public data root, the check's reports and the review's own runs. To run them again, change those paths first.

## Not here

The review's own report, its other result files (the removal and root checks, the arm "with" against the release, its reruns, the interval seeds, the counts by single PMID), its other scripts, its build and its runs stay with its output on the machine that ran it. None of its files with one row per sample or genome is copied.

## Limits

- The fold-plan runs and the October search are exploratory. The fold-plan runs were fixed in `ADDENDUM.md` before they ran; the October search was decided after the results were seen.
- The fold plans were varied for `staphylococcus_aureus__clindamycin` only.
- A change when a study leaves is an association with that study's labels or genomes. It does not show which label is wrong, and it is no statement about any gene.
