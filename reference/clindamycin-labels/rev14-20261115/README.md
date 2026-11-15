# Clindamycin labels by study: the result of 15 November 2026

This folder holds the files behind [docs/clindamycin-labels.md](../../../docs/clindamycin-labels.md): the plan, its addendum, the tables the page quotes and the scripts that wrote them. The check asked whether one deposited study holds a large share of the disagreements between the reported clindamycin phenotype and GRM's known-gene rule in the six clindamycin datasets of release `rev14-274-20261112`, and whether the released figures depend on that study.

It is a check of labels and scores on public BV-BRC records. A change when a study leaves is an association with that study's labels or genomes. The check names no gene.

Every file here but this README and `SHA256SUMS` is copied unchanged from the check's own output, and `SHA256SUMS` gives the sha256 of each. To check them, run from this folder:

```sh
sha256sum -c SHA256SUMS
```

## How it was made

| | |
|---|---|
| Release | `rev14-274-20261112`, evaluation revision 14: the `samples` of its six clindamycin reports |
| Data | a copy of the public data root: `curated/phenotypes.jsonl`, `curated/drug_names.json` and the raw BV-BRC records of `bvbrc/genome_amr.jsonl` |
| Program | `grm` built at `852ebc9` with Rust 1.98.1, binary sha256 `b20b5f1415a011293a63f9fc8b90ef54b3f0935952cbcfc67aa0a36c5bcc310e` (`PLAN-ADDENDUM-20261115.md`, A1) |
| Plan | `PLAN.md`, recorded 2026-11-14T18:11:52Z, and `PLAN-ADDENDUM-20261115.md`, recorded 2026-11-15T11:36:42Z, both before any count or score was read |
| Run | 2026-11-15: counts at 11:43:22Z, fits from 11:45:10Z to 11:57:40Z, scores from 11:58:16Z to 11:58:58Z, tables at 11:59Z (UTC) |

Every `grm` call ran with no network. The scripts name the paths of the machine they ran on (`RELEASE` and `ROOT` in `count_by_study.py`, `S`, `PUB` and `OUT` in `arms.sh`, and `S` and `OUT` in `run_all.sh`), so to run them again, change those paths first.

## Files

| File | What it is | Written by |
|---|---|---|
| `PLAN.md` | The plan: datasets, the study field, the rule that picks a study, the measures and the caveat rule | by hand; hashed at 2026-11-14T18:11:52Z |
| `PLAN.sha256` | The plan's sha256, with the time it was recorded | `sha256sum PLAN.md` |
| `PLAN-ADDENDUM-20261115.md` | The addendum for the restart on another machine. It changes no rule of the plan | by hand; hashed at 2026-11-15T11:36:42Z |
| `PLAN-ADDENDUM-20261115.sha256` | The addendum's sha256, with the time it was recorded | `sha256sum PLAN-ADDENDUM-20261115.md` |
| `results/by-study.tsv` | One row per dataset and study key (the sorted set of a sample's PMIDs, or `no-pmid`), with a total per dataset; the rows of a dataset add up to it | `python3 -I scripts/count_by_study.py` |
| `results/by-pmid.tsv` | One row per dataset and single PMID; a sample with two PMIDs counts under both | the same command |
| `results/chosen.tsv` | The PMID the plan picks in each dataset, its shares of the disagreements and of the samples with exact intervals, the rates in the study and the rest, and the class counts left without it | the same command |
| `results/count-checks.txt` | The checks of `PLAN.md` section 2: each sample's curated row against its report, and what the raw study field held | the same command |
| `results/scores/<dataset>/acceptance.txt` | Arm "with" against the release, and arm "without" against arm "with", for the three refit datasets | `python3 -I scripts/score_arms.py --released … --with … --without … --removal … [--placebo …] --out results/scores/<dataset>/` |
| `results/scores/<dataset>/caveat.tsv` | The caveat rule row by row: each split, row and measure, its change, its lineage-cluster interval, whether it passes, and its placebo range | the same command |
| `results/refit-table.tsv` | For each refit dataset, split and row: recall, log loss and balanced accuracy with and without the study, the change with its interval, and the parts of the change in balanced accuracy | `python3 -I scripts/collect.py` |
| `results/caveat.tsv` | Per dataset: the existence and dependence parts of the caveat rule, and rows passing of rows tested | the same command |
| `results/logs/released-reports.sha256` | The sha256 of the three released reports that were refit | after the runs (file time 12:00:40Z) |
| `results/logs/root-checks.txt` | The comparison of each scratch data root with the public root: lines removed from `curated/phenotypes.jsonl` and every other curated file equal | after the runs (its first line gives 2026-11-15T12:00:50Z) |
| `scripts/count_by_study.py` | Step one: counts by study, with no fit | |
| `scripts/stats_util.py` | The exact interval and Fisher's exact test, standard library only | |
| `scripts/filter_phenotypes.py` | Copies `curated/phenotypes.jsonl` without the clindamycin rows of the listed genomes | |
| `scripts/placebo_ids.py` | Draws the placebo removals | |
| `scripts/arms.sh` | The scratch data roots and GRM's own commands: `grm count`, `structure`, `matrix` and `learn --only <dataset>` | |
| `scripts/run_all.sh` | Every fit, in order; it calls `arms.sh` | |
| `scripts/score_arms.py` | Scores the arms of one dataset | |
| `scripts/collect.py` | Gathers the scored datasets into `results/refit-table.tsv` and `results/caveat.tsv` | |

The three lines of `results/logs/released-reports.sha256` equal the lines of the same reports in [`reference/lineage-control/rev14-274-20261112.reports.sha256`](../../lineage-control/rev14-274-20261112.reports.sha256), the release's own checksum lines.

`PLAN.md` is here as it was hashed, so it keeps the October record's own words for that record's study. The page names each study by its PMID alone.

## Not here

The check also wrote one row per sample with its label (`results/membership/`), the genome ids each arm removed (`results/removal/`, `results/placebo/`), the scores' other files (`views.tsv`, `effects.tsv`, `placebo.tsv`, `summary.json`), the logs of the build and of the runs, an exploratory literature lookup made after the counts were read, and the check's own report. They stay with the check's output on the machine that ran it. Sample ids with their labels are public, and they are still left out of the repository.

## Limits

- The change when a study leaves is an association with that study's labels or genomes. It does not show which label is wrong. It does not separate label quality from the study's lineages, testing method or sampling. It is no statement about any gene.
- Each arm draws its own fold plan over its own samples, and the bootstrap intervals do not cover that change. Placebo arms could be drawn for `campylobacter_coli__clindamycin` only.
- The 0.05 threshold and the words "most" and "a large share" are adapted research screens fixed in the plan.
- The scripts read files that are not in the repository: the release and a copy of the public data root.
