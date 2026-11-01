# Sampling audit

This directory holds a read-only audit of how `grm curate` picks the genomes each species downloads. It ports the legacy `select_genomes` (`src/curate.rs`, unchanged from 7711459 through c90ebaf) to Python. It then reports, per task (species x antibiotic) and per cap, what the cap keeps. Nothing here changes curation, a cohort or a production file.

The ranking scores each genome by summing `1 / class size` over its reported R and S calls on the species' dataset drugs. So genomes with calls on many drugs come first, and genomes whose only labels are MIC-inferred are never picked. One drug's cohort also depends on the other drugs' labels. A 4,000 cap just takes the next 2,000 genomes of the same ranking, so a small gain from 2,000 to 4,000 says little about whether independent lineages or minority classes are exhausted. The audit makes that visible.

## Files and commands

- `legacy_select.py` is the port. `python3 -I legacy_select.py CURATED_DIR SPECIES CAP` prints the ids in the order `genomes.txt` holds them.
- `audit.py` writes the report: `python3 -I audit.py OUT.json --curated DIR --caps 2000,4000[,all] [--groups NAME=TSV ...] [--species NAME ...] [--tsv OUT.tsv]`.
- `compare_frozen.py` checks the port against a frozen run: `python3 -I compare_frozen.py OUT.json RUN CAP [--sums SHA256SUMS]`.
- `groups_from_metadata.py` writes a group TSV from one BV-BRC metadata field: `python3 -I groups_from_metadata.py OUT.tsv bvbrc/genomes.jsonl cgmlst_hc100`.
- `tests/test_sampling_audit.py` runs the fixture tests: `python3 -I scripts/sampling-audit/tests/test_sampling_audit.py`.

All four scripts need only Python 3.12 and its standard library.

## What precedes selection

Selection sees only what `curate` keeps. In order, `curate` does the following:

1. It drops records whose antibiotic name resolves to a drug class, a non-antibiotic or nothing.
2. It drops genomes without metadata, then genomes that `grm species-check` failed (before any other QC), then genomes without a two-word species.
3. Then it applies genome QC: BV-BRC quality Poor, more than 1,000 contigs, CheckM completeness below 90 or contamination above 10, no genome length, or length more than 25% off the species median. The median counts failing genomes too.
4. It drops records with neither a parsable call nor an MIC.
5. Per genome and drug, it merges records. R with S is a conflict with no reported call. R with I gives R, S with I gives S, and I alone gives I.
6. Where the inferred breakpoint applies, it labels pairs without a reported call from MIC (`label_from_mic`).
7. Datasets are the species x drug pairs with at least `--min-class` (default 25) unconflicted, reported, not MIC-inferred R calls, and as many S.

`select_genomes` then runs once for each species that has a dataset. Only that species' dataset drugs count, and only `reported_label()` rows: reported R or S, not conflicted, not MIC-inferred. A genome with no such row is never scored, whatever the cap.

## The port's input contract

The port reads the two files `curate` writes and passes to `select_genomes`.

- `phenotypes.jsonl`: one row per genome and drug after steps 1 to 6. It uses `genome_id`, `species`, `antibiotic`, `reported` (`"R"`, `"S"`, `"I"` or null), `conflict`, `label_from_mic` and `label`. Row order matters, since each genome's score adds its terms in file order.
- `datasets.json`: it uses `species` and `antibiotic` (and `id` for the report).
- A species name and a cap.

Matching the Rust exactly takes three details:

- Scores add left to right in f64. Python 3.12's `sum()` compensates rounding, so the port adds in a loop.
- The sort is score descending, then genome id ascending by UTF-8 bytes, as `f64::total_cmp` and `str::cmp` give.
- Class sizes count rows, so a duplicated row would count twice, as it would in Rust.

## Report fields

`audit.py` writes JSON with `schema` (2), `caps`, `group_sets`, `min_group_coverage`, `inputs`, a `species` list and a `tasks` list. `--tsv` adds one line per task and cap with the main figures, including each group set's `rows_without_group`, `group_coverage` and the null reasons for the ESS and both added-group counts.

Per species:

- `target_antibiotics`: the dataset drugs that score.
- `scored_genomes`: genomes with at least one reported R or S call on a dataset drug. Only these can be picked.
- `unscored_genomes`: genomes with rows on dataset drugs but only MIC-inferred, intermediate, conflicted or unlabelled ones. No cap picks them.
- `panel_breadth_scored`: the distribution of panel breadth over scored genomes. Panel breadth is the number of dataset drugs with a reported R or S call.
- `by_cap.<cap>` holds these fields:
  - `cap_binds`: whether there are more scored genomes than the cap.
  - `selected`: the number of genomes selected.
  - `cutoff_score`: the score of the last genome kept.
  - `tied_at_cutoff`: the scored genomes with exactly that score.
  - `tied_at_cutoff_selected`: how many of those were kept.
  - `tie_broken_by_genome_id`: true when genome ids, not labels, decided who is in.
  - `panel_breadth`: the breadth distribution of the selected genomes.

Per task, rows are the genome x drug rows with a `reported_label()`, the rows `grm learn` trains and scores on (`analysis::cohort`).

- `minority` is the class with fewer eligible rows. A tie goes to R and sets `minority_tied`.
- `eligible` covers all such rows of the species. `by_cap.<cap>` covers the rows of the genomes that cap selects. Each holds these fields:
  - `rows`, `R`, `S`: row and class counts.
  - `prevalence_R`: R / rows.
  - `minority_rows`: rows of the minority class.
  - `panel_breadth`: the breadth distribution over the genomes behind the rows.
  - `task_score_share_mean`: the mean share of each genome's score that comes from this drug. A low value means other drugs' labels mostly put these genomes in the cohort.
  - `non_target_rows`: the drug's other rows (MIC-inferred, intermediate, conflict, unlabelled), counted over the species for `eligible` and over the selected genomes for a cap.
- `by_cap.<cap>.minority_coverage`: selected minority rows / eligible minority rows.
- `groups.<set>` appears only with `--groups`. Ungrouped rows are counted, then left out of every other group figure. It holds these fields:
  - `rows_with_group`, `rows_without_group` and `group_coverage`: rows with and without a group id, and the share with one.
  - `groups`: the distinct groups among the grouped rows.
  - `minority_groups`: groups with at least one grouped minority row.
  - `largest_group_rows` and `largest_group_share`: the largest group's size and its share of the grouped rows. The share is null when no row is grouped.
  - `group_ess`: the inverse Simpson number of group sizes over grouped rows, (sum of n) squared / sum of n squared. It is null when `group_coverage` is under `min_group_coverage` (0.8), or when there are no rows.
  - `group_ess_null_reason`: why `group_ess` is null, for example `0 of 3572 rows grouped, under 80%`, or `no rows`. With no rows, `group_coverage` and `largest_group_share` are null too.
  - `minority_rows_without_group`, `minority_group_coverage`, `minority_group_ess` and `minority_group_ess_null_reason`: the same over minority rows only, with the same 0.8 gate on minority coverage.
- `nested` holds one entry per pair of consecutive caps (`from_cap`, `to_cap`) with these fields:
  - `genomes_added` across the species.
  - `rows_added`, `R_added`, `S_added` and `minority_rows_added` for this drug.
  - With groups, `groups.<set>.groups_added` (groups first present at the larger cap) and `minority_groups_added` (groups that first hold a minority row there). Both count grouped rows only.
  - `groups_added` is null unless both caps pass the 0.8 coverage gate, and `groups_added_null_reason` names the failing caps, for example `cap 2000: 0 of 1899 rows grouped, under 80%; cap 4000: 0 of 3572 rows grouped, under 80%`. `minority_groups_added` and `minority_groups_added_null_reason` do the same with the minority coverage gate.

The ranking does not depend on the cap, so each smaller selection is the start of the larger one. The audit asserts that.

## Caveats

Every figure is descriptive. None is a test, an interval or an estimate of model performance.

- `group_ess` is the number of equal-sized groups with the same concentration as the observed sizes. It equals an effective sample size only if rows within a group are fully dependent and groups fully independent. Neither duplicate groups nor lineage clusters guarantee that.
- `group_ess` also ignores labels, phenotype noise and structure above the grouping. It moves with the grouping's resolution (k-mer-identical duplicates versus a lineage cut), so compare it only within one grouping.
- `minority_groups` counts groups with at least one minority row. A group that holds one resistant isolate among hundreds of susceptible ones counts the same as an all-resistant group.
- "Selected" means listed in `genomes.txt`. The model's cohort is smaller: genomes with no FASTA or a failed count drop out after download, and `learn` applies `--min-class` again.
- `prevalence_R` in a selection mirrors the ranking's preference for rare classes and wide panels, not the source population. BV-BRC submissions are not a population sample to begin with.
- Panel breadth counts only dataset drugs with reported R or S calls, the quantity the score sums over. It is not the number of drugs a lab tested.

Group figures describe the grouped rows only, so read `group_coverage` before any of them. Genomes without an id are not a random subset: in the r9 metadata, whole species have no cgMLST. `groups` and `minority_groups` are always given, and they are lower bounds on the groups behind all rows, because the groups seen are a subset of the true ones.

`groups_added` is not a bound either way. An ungrouped genome at the smaller cap can belong to a group that then looks new at the larger cap. An ungrouped genome at the larger cap can hide a group that really is new. That is why it is gated on both caps.

The 0.8 gate on `group_ess` and `groups_added` is a judgement call, not a derived threshold. Below it the figure is withheld. Above it the ungrouped rows can still shift the figure.

`tie_broken_by_genome_id` flags ties in the float scores, which is what decides Rust's order. Two genomes whose exact scores are equal but whose float sums differ in the last bit are ordered by rounding, and the flag misses them.

check2k and check4k each have one such case, the same one, in A. baumannii, where no cap binds. `470.7504`, `470.7505` and `470.7494` have equal exact scores; the first two share their calls. `470.7494`'s float sum is one bit lower, so it ranks 164th, after both, instead of first by id.

At every binding cutoff in all four runs (base2k, check2k, base4k and check4k), the float and exact tie counts agree, and every flagged tie is between genomes with the same calls on the same drugs. I checked this on 2026-10-30 with exact fractions.

## Export fields still needed

The ported selection and every non-group field run fully on the existing curate outputs, `curated/phenotypes.jsonl` and `curated/datasets.json`.

Group metrics need a group id for every scored genome, not just the selected ones. `grm structure` clusters only downloaded genomes, and kset duplicate groups exist only for counted genomes. So GRM's own groups never cover the genomes a larger cap would add.

BV-BRC metadata already carries `mlst` and `cgmlst_hc5`, `cgmlst_hc20` and `cgmlst_hc100`, which `groups_from_metadata.py` turns into group files. The r9 group files came from check4k's metadata, read only:

- Path: `<r9-curate>/check4k/bvbrc/genomes.jsonl`, a symlink to `<r9-species-check-data>/bvbrc/genomes.jsonl`. `<r9-curate>` is the folder holding the four frozen r9 curate runs, `grm-r9-curate/` on the laptop, and `<r9-species-check-data>` is the r9 species-check data folder, `grm-r9-species-check-data/` beside it. Neither is in the repository.
- sha256: `d60e2b15f2cc0ed0466a1016da1761869fb42cc7e9fdd817c1ff7051e759a80e`.

Coverage depends on the denominator. Over metadata genomes (records whose BV-BRC `species` field names the species, before QC), `cgmlst_hc100` covers:

- M. tuberculosis: 65.8% of 35,014.
- E. coli: 96.1% of 12,686.
- Salmonella: 97.2% of 9,370.
- K. pneumoniae: 97.3% of 8,463.

Over the genomes the audit scores (check4k, the pool every cap draws from), the same species reach 68.9% of 21,293, 97.6% of 8,700, 98.0% of 8,169 and 98.7% of 6,091. Ten species have no cgMLST at all: C. diphtheriae, E. cloacae, K. aerogenes, N. gonorrhoeae, S. flexneri, S. sonnei, S. haemolyticus, S. agalactiae, S. pneumoniae and V. cholerae. MLST covers at least 84.9% of scored genomes in every species. These ids are BV-BRC's, unchecked by GRM, and MLST barely separates M. tuberculosis.

A full run on real cohorts still needs these exports:

1. Lineage or duplicate group ids over the whole scored pool, with method and threshold. For example, a TB-Profiler lineage, or sketch clusters and kset duplicates computed before selection.
2. A per-genome eligibility record from `curate`: metadata present, species-check verdict, resolved species, and QC verdict and reason. Today `qc.json` lists failing genomes only, and genomes without metadata are only counted.
3. Per genome and drug, the raw records dropped by name resolution or for having no measure. `funnel.json` counts them only in total.
4. The `--cap`, `--min-class` and binary commit of each curate run, written into the curated directory. Today they live in logs and directory names.
5. Per selected genome, whether it reached the matrix (FASTA present, count succeeded), so selected rows can be followed to trained rows.

## Checked against frozen runs

On 2026-10-30 the port reproduced the four frozen r9 curate runs in `<r9-curate>`, the folder named above. For base2k and check2k at cap 2,000, and base4k and check4k at cap 4,000, all 23 species matched `genomes.txt` id for id and in order. Every `genomes.txt` matched the run's `SHA256SUMS`. The inputs were each run's `curated/phenotypes.jsonl` and `datasets.json`. Those two files are not in `SHA256SUMS`, so `compare_frozen.py` records their digests. Each comparison took 8 to 26 s. The 2026-10-30 fix round re-ran all four with the same 92 of 92 result.
