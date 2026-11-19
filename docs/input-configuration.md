# Input configuration: today's values and what sets them

Layth specified a full-pipeline comparison over the input data: genome and cohort count, which genomes are included, k-mer length, feature families and the retained-feature budget. Each input now has a setting that reaches `grm`, and [input-comparison-runner.md](input-comparison-runner.md) runs the comparison main decided under Layth's 5 Nov 2026 delegation (DECISIONS V-1). It has run on synthetic data only; real cohorts are H-1's. This page gives each input's default, the flag that sets it, and what records it.

| Input | Default | Flag | Recorded in | Arms (V-1) |
|---|---|---|---|---|
| Genome and cohort count | at most 2,000 genomes per species; at least 25 per class per dataset | `grm curate --cap`, `--min-class` | qc.json, datasets.json and each report (`selection.cap`, `counts.min_class`) | cap 1,000 and 4,000; class minimum stays 25 |
| Which genomes are included | `legacy`: `curate::select_genomes`, a label-weighted ranking | `grm curate --selection legacy\|group-round-robin-v1 --selection-seed [--stratum]` | qc.json, datasets.json and each report (`selection`), and `species/<slug>/selection.json` for a group draw | `group-round-robin-v1` |
| k-mer length | 31 | `grm pipeline --k`, `grm count --k` (5 to 32) | counted.json, structure.json, matrix.json and each report (`counts.k`) | 21 and 25 |
| Feature families | the K block (k-mer patterns) only | `grm arm --blocks` (K, C, H, P) | `grm arm`'s `arm.json` (block specs) and K's run manifest | k-c, k-h, k-p, k-c-h-p (blocks K+C, K+H, K+P, K+C+H+P) |
| Retained-feature budget | 4,000 screened features per fit | `grm learn`, `grm pipeline` and `grm arm --screened-features` | each report's `eval_config` and input fingerprint; `arm.json` | 1,000 and 16,000 |

`grm pipeline` and `grm learn` at their defaults write the same report as before (the lane's byte-identity gate). Learner settings are a separate thing: `scripts/hyperopt` searches learner knobs only (`scripts/hyperopt/protocol.py`), and none of these five inputs.

## Genome and cohort count

**Default.** `grm curate` selects at most 2,000 genomes per species (`--cap`). A species × drug dataset exists only when each class has at least 25 genomes with an unconflicted reported call (`--min-class`; `QcConfig`, `src/curate.rs`).

**What reaches them.** `run_curate_with` (`src/main.rs`) refuses a class minimum below the outer fold count, puts it in `QcConfig`, which `curate::curate` applies to each dataset and records in it, and passes the cap to the selection rule once per species. `grm learn` applies each dataset's recorded class minimum again after download (`analysis::cohort`).

**What is recorded.** The class minimum is in qc.json, in each dataset of datasets.json and in each report (`counts.min_class`). The cap is now recorded with the rule: qc.json's first element carries `selection` {rule, cap, seed, always_included}, each dataset in datasets.json carries the same plus the group plan's `input_sha256`, and `grm learn` copies a dataset's `selection` into its report (`DatasetReport::selection`). A report from a dataset curated before the record existed has no `selection` field, so such a report is byte for byte what it was.

**In the pipeline.** `grm pipeline` doesn't run `curate`; it uses the `genomes.txt` the last `grm curate` wrote. A data root holds one curated cohort (`src/layout.rs`), so the comparison gives each cohort its own data root.

## Which genomes are included

**The legacy rule** (`curate::select_genomes`) scores each genome by summing 1 / (size of that call's class) over its reported R or S calls on the species' dataset drugs and keeps the top `cap`, ties broken by genome id. `--selection legacy` is the default, and its choice and `genomes.txt` bytes are unchanged (tests `legacy_selection_writes_the_same_genomes_txt_bytes` in `curate.rs` and `curate_records_the_legacy_rule_and_its_cap_without_changing_its_choice` in `main.rs`).

**The group rule** (`--selection group-round-robin-v1`) runs F's `plan_group_sampling` once per species, as its contract requires:
- eligible genomes are the ones the legacy rule ranks (`curate::labeled_genomes`: an unconflicted reported call on a species dataset drug), so both rules choose from one set;
- each genome's group is its BV-BRC MLST sequence type, the `mlst` field of its genome record, taken whole after trimming (`mlst:<value>`; a value that joins several schemes' types, such as E. coli's Achtman and Pasteur, forms one group per combination); a genome without one is its own singleton group (`singleton:<genome id>`), so none is refused. The `Selection` record names the source (`group_source` `bvbrc-genome-mlst`) and counts the species' singletons (`singletons`), per the rulings after lane INPUT;
- `--stratum` is an optional `AlwaysIncluded` JSON (genome ids, reason, `label_dependent`); each species gets the part of it that is its own, and a stratum genome eligible in no species is refused;
- `--selection-seed` seeds the draw. `genomes.txt` keeps the plan's draw order, and the whole plan, with inclusion probabilities, goes to `species/<slug>/selection.json` (removed again by a later legacy curate).

Legacy with `--selection-seed` or `--stratum`, or the group rule without a seed, is refused before anything is read. The comparison's group arm uses seed 1, a decided value.

## k-mer length

**Default.** 31 (`PIPELINE_K`, `src/main.rs`).

**Flag.** `grm pipeline --k` and `grm count --k` take 5 to 32 (`src/dna.rs`). `grm pipeline` uses its k to check the catalog, in `count_to_matrix` (catalog, stale-stage check, count) and for the recount after unreadable ksets. Stages built at another k are recounted and rebuilt at this one (`species::stale_stages`; test `the_pipeline_counts_and_builds_at_the_k_it_is_given`).

**What reads k.** The catalog gene calls (loaded in a k-mer space of that k), the stale-stage check and kset reuse, the stage folders (one k per species in a data root), structure's MinHash sketches and so the lineage labels, the kset digests behind the duplicate groups, and the report's region probe sets. A change of k can therefore change the lineages and duplicate groups the folds read; that is why the comparison's k arms take the k = 31 fold plans (`k31_groups`).

## Feature families

**Default.** The learners read the K block only: `grm learn` evaluates the dataset's projection of the species `FeatureMatrix`.

**Flag.** `grm arm --blocks` admits K, C, H and P through H's adapters (`src/ablate/arm.rs`, `block`):

| Block | What it is | Adapter, as `grm arm` builds it | State |
|---|---|---|---|
| C | curated exact calls (region probes) | `ProbeSource::from_ksets` over the species' ksets | present where the species has a callable assay |
| H | hierarchical panel calls | `PanelSource::new` over the counted catalog gene calls | present where a panel reports a target read for the drug in the species; otherwise unavailable, with the reason per panel |
| P | catalog gene-group presence | `PresenceSource::new` over the ksets | present with a catalog |
| N | gene copy counts | `CountSource` | unavailable: nothing in GRM counts gene copies |
| X | genome context | `ContextSource` | unavailable: no reviewed genome-context source |
| K | k-mer patterns | `KmerSource` over the species matrix | present |

A block whose source is unavailable makes the arm's runs `unavailable` with that reason. Building an N or X source is biological feature engineering and out of scope; the conditions a future one must meet are in `src/ablate/blocks.rs` (reviewed, sealed adapter, label-free, proxy refusal, no filled rows, `version` and `source_sha256`).

## Retained-feature budget

| Value | What it limits | Default | Flag |
|---|---|---|---|
| `screened_features` | chi-square-screened features offered to the forest, boosting and the lasso, per fit, on that fit's own fitting samples | 4,000 (`EvalConfig::default()`, `DEFAULT_SCREENED_FEATURES` in `main.rs`, held equal by a test) | `--screened-features` on `grm learn`, `grm pipeline` and `grm arm` (at least 1) |
| `K_SCREEN_CAP` | the K block's declared screen cap (`blocks.rs`) | 4,000 | none; `grm arm` applies the learners' budget only |
| matrix `min_support` | k-mers kept in the species matrix | 2 (`PIPELINE_MIN_SUPPORT`) | `grm matrix --min-support` |
| dataset `min_support` | each fit's candidate features, and `grm arm`'s encoder support floor | max(n / 200, 3) | none |
| lasso `max_nonzero` | the lasso path's stopping size | 400 | none |

Each report's `eval_config` records `screened_features`, and the input fingerprint changes with it.
