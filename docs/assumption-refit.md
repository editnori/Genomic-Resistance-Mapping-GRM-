# The refit with fold models, D0 and C8

The assumption plan's diagnostic D0 asks why a learner misses on a lineage fold where the AMRFinderPlus program does not. Its block C8 asks how stable the learners' chosen inputs are across folds. Both read each outer fold's model, and none of the 274 reports of release `data-rev14-274-20261112` holds `fold_models` (`docs/assumption-plan-corrections.md`, "Fold models").

So the release's rule is fitted again with fold models saved and checked equal to the release, and D0 and C8 read that refit. No report of the release changes, `EVALUATION_REVISION` stays 14 and the report schema stays 2.

D0 is a diagnostic: its datasets are chosen by their results, and its counts order the plan's comparisons C3, C4 and C5 and test nothing. Every output is an aggregate: counts, rates and dataset names. No genome id, k-mer, sequence, position or feature list leaves the run box, and no new marker is found, named or ranked. A catalogue family is named only as the public AMRFinderPlus catalogue names it. GRM's known-gene rule, GRM's own reader of that catalogue, keeps its own name.

## The readers

`scripts/assumption-d0-c8/` holds them. They fit nothing.

**MISS.** `miss_sets.py` gives the plan's set A and its MISS set from the comparison record (`reference/amrfinderplus-comparison/`): 75 datasets in 16 species, each learner on its own set (scm 26, cart 34, rf 59, gbt 48, lasso 37).

**D0.** `d0.py dataset` reads one refit report with its fold models, the dataset's `grm export` directory and the program's per-genome calls. For each learner in MISS it finds the missed lineage folds and puts each in a cause: Absent, Removed, Lineage, Chosen or Unused, in that order. The readings the plan does not fix are in the corrections page's entry "D0 causes".

**Inputs.** D0 recomputes the floor and each learner's 4,000-row screen by `grm`'s own rules. It stops a dataset when its inputs disagree: genome ids, labels, lineages, floor, screen size, input fingerprint, the record's counts and each fold's balanced accuracy.

**Combine.** `d0.py combine` adds the datasets up per learner and writes the plan's ordering line. With `--marks` it adds the rows for the datasets marked M and the rest. With `--expect-datasets` and `--expect-pairs` it stops when a dataset table is missing or read twice. `d0.py release` reads a released report with no fold models and leaves Lineage, Chosen and Unused open.

**C8.** `c8.py` gives, per dataset and learner, the median pairwise Jaccard of the five lineage folds' and the five random folds' chosen feature sets, the distinct features chosen, and how many of those have a card that names a catalogue family. With `--final` it reads the final models of a report with no fold models.

**C8 combine.** `c8.py combine` reads each species' `c8-datasets.tsv` and writes, per learner and over the five, the median over datasets of that overlap with its quartiles and the units where fewer than two folds chose, the features summed beside it, and the units left undefined with the reason: the split did not run, no fold chose a feature, or, when asked, one fold alone chose. Three rules the plan leaves open are its arguments, with no default: a dataset where one fold alone chose (`--one-fold`), the row over the five (`--over-five`), and counts under and over two cut points (`--bands`). It refuses any cell that is not a count, an overlap, a mark, a learner or a dataset name, and never prints that cell.

**All species.** `read_refit.py` takes the species' copied directories and writes a page's tables: the species, D0 per learner and per dataset, the ordering line and C8. It stops and writes nothing unless each species holds `DONE` with the pinned code sha256 and files equal to its `SHA256SUMS`, the check tables read equal, the species and datasets number the counts given with none twice or missing, and D0's datasets are the MISS set among each species' listed datasets with the dataset and pair counts given. `d0.py combine` adds D0 up. Beside the ordering line it says, from the record's MISS set, how many of the plan's MISS datasets and miss pairs D0 reads, and whether the pairs of the MISS datasets that are M could change the order (the reading's ruling 6). Species and antibiotics are named in full, from `scripts/outside-study/frozen/datasets.tsv`. The reader and both combines use the standard library alone.

## The refit

The refit is `grm learn` at the release's seed, 42, and its settings, on the 238 datasets that are not marked M. It runs on a data root built as the MIC comparison's run builds one (`docs/mic-plan.md`, "Running one species").

That root holds the rebuilt revision 14 table with each dataset's breakpoint set to its release report's, and the catalog equal to a pinned list. The species' genomes are fetched again from BV-BRC with the release's length, contig count and N50. Count, structure and matrix then run with no network, and the build's code sha256 is recorded beside its commit.

**Equality.** `scripts/refit/refit_check.py` decides whether a refit report is the release's rule. Every top-level field of the release's report but `seconds` must be equal, with each sample's fitted fields, each split's folds, and every field but `seconds` of each learner row.

**Cards.** The top-level `features` list is compared card by card, in order. A card's `kmer` and `locus` may differ: both are written after the fit from genome files and CDS lists fetched again from BV-BRC, and a file served today can differ in bytes from the one the release fetched. Every other card field, and the number and order of the cards, must be equal. `datasets.tsv` counts the cards whose `kmer` or `locus` differs in its last column, `cards_with_other_locus`. A `kmer` or `locus` key anywhere else is compared as any field.

**Allowed.** Three new fields pass: `fold_models`, which must hold one model per fold for each of the five learners on each split that ran, each sample's `rule_groups`, and `selection`. A sample's country, year, host and MLST type are counted apart, since the refit's genome metadata is left empty and no fit or fingerprint reads them. With those fields taken out and each card's `kmer` and `locus` values set to one mark, both reports written again in their own key order must be the same text. The check writes counts and field names only.

## Running one species

`scripts/refit/run_species.sh KIT SPECIES DATASETS OUT [JOBS]` runs on a box of its own with no agent. DATASETS is the species' datasets that are not M, as a comma list. It refuses a dataset that the label count's table marks M.

Its steps up to the matrix are `scripts/c6-mic/run_species.sh`'s in the same text. The kit is pinned by sha256, OUT takes one start, and a changed build or dataset list is refused. A species with a D0 dataset is refused at the kit step, before any fetch, unless `python3 -I` imports numpy 1.26.4 on python 3.12 (`python3-numpy` by apt). Each step's log is copied to `OUT/logs` with every genome id written as `[genome id]`. Then:

| Step | Work |
|---|---|
| learn | `grm learn` over the listed datasets whose report is missing or has no `fold_models`; every listed report must then hold them |
| check | `refit_check.py`; any difference but a feature card's `kmer` and `locus` stops the species before D0 or C8 reads a report |
| C8 | `c8.py` over the species' refit reports, before the program, so a failure in the program or D0 still leaves C8's tables |
| program | for the species' MISS datasets among the listed: the program installed from the AMRFinderPlus kit and proved on the tool's own test file, then run on every genome of the species that has a table recorded for the release, each table equal to the sha256 recorded for it. The release ran the program on the genomes of its rule datasets, so a genome in none has no table: *Staphylococcus aureus* lists 2,000 genomes and 1,994 have one |
| D0 | per dataset, `grm export` on the same data root and `d0.py dataset` on the refit report, then `d0.py combine` with the species' own dataset and pair counts |

**End.** A finished species writes `DONE <species> <time> binary <sha256> code <sha256> d0 <n> datasets SHA256SUMS <sha256>`. A species with no MISS dataset among the listed skips the program and D0, and its line says `d0 skipped: no MISS dataset among the listed` in place of `d0 <n> datasets`. Any failure writes `OUT/FAILED` naming the step and exits 1. Run again, the script skips finished steps and fits only the reports that are missing or hold no `fold_models`. It stops no process but its own.

**Restart.** A restart on a later checkout goes on when the built binary's sha256 still equals `binary.txt`'s and the commits differ only under `scripts/refit/`, `docs/`, `tests/` and `.github/`. Any other changed path refuses, `src/` and `scripts/assumption-d0-c8/` among them. The restart keeps the fitted reports, runs the check again and adds a line `restart <commit> <time>` to `binary.txt`.

**Copy.** The refit reports, the data root, the exports and the raw logs hold genome ids and stay on the box. After a count of genome ids over exactly these files, they may leave it: `binary.txt`, `expect.txt`, `d0.txt`, `SHA256SUMS`, `DONE`, `FAILED`, `logs/`, `check/` and `results/`. D0's tables and the program's adapter summary are under `results/d0`, and C8's are under `results/c8`.

**Check.** `scripts/refit/check_run_species.sh KIT RELEASE WORK` checks the run script before a run, in a few minutes, with a stub in place of `grm` that writes SYNTHETIC reports. A run with no fault ends at the check. A copy of the release's reports with null fold models runs to `DONE` on *Klebsiella quasipneumoniae*, none of whose datasets is in MISS, and says D0 was skipped.

**Faults.** The same script runs a failed fit, a missing report, a report without fold models, a changed checkout or binary, a changed dataset list or kit file, and a second start. A restart after a commit under `scripts/refit/` must go on and record it, and one under `src/` or `scripts/assumption-d0-c8/` must be refused. With numpy hidden, a species with a D0 dataset must stop at the kit step and one with none must run. It lists *Acinetobacter baumannii* with amikacin, which is marked M. It puts a genome id in a retry line, which must leave none in the files that may leave the box.

**Result.** The reading of 22 Nov 2026 over the 23 species is `docs/refit-result.md`, with its tables in `reference/refit-result/`. `page_tables.py` prints that page's tables from those files, with shares and overlaps rounded to 3 decimals.

Tests: `tests/assumption_d0_c8` and `tests/refit`, on SYNTHETIC reports, exports and calls. One file differs: `tests/assumption_d0_c8/test_page.py` checks the result page against the committed tables of `reference/refit-result/`.
