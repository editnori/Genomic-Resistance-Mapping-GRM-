# Ablation report: schema `ablation-report` 1

`scripts/ablation/schema.py` validates the records a feature-block ablation produces, and `scripts/ablation/report.py` turns a validated set into a JSON report plus a Markdown summary. Both use only the Python standard library. They make no release or clinical claim: a passing validator means the records are well formed and consistent, not that a result is true.

Both implement the GRM ablation data contract v1.10 (main orchestration ruling, 2026-10-30, sha256 `1ce28f538e5ad8a0652cedcd1642389e4fb6c96b3c9e3da1c76aa997bb98cd78`), including its §9 (v1.2), §10 (v1.3), §12 (v1.5), §13 (v1.6), §15 (v1.8) and §17 (v1.10) rulings. §11 (v1.4) and §14 (v1.7) concern L's outputs, and §16 (v1.9) concerns Task I; none changes a K check. §12 replaces §2's "pooled bootstrap" and §10's "such as" method list; §13 refines §9's undefined-resample rule and fixes what a pooled record pools; §15 withdraws `equivalent` as a decision, with §9's and §12's equivalence floors; §17 makes pooled keys planned ids, has K recompute the multiplicity adjustment, enforces the few-species floor and fixes how the report is written. Where this page or the code disagrees with the contract, the contract wins.

```sh
python3 -I -B scripts/ablation/report.py MANIFEST CONTRASTS POOLED FEATURES PLANNED PLANNED_POOLED OUT_DIR \
    [--not-produced FILE] [--predictions-root DIR]
python3 -I -B -m unittest discover -s tests/ablation -t .
```

`report.py` reads its inputs as UTF-8, and refuses a JSON object with a repeated key and input nested more than 100 levels (`report.MAX_JSON_DEPTH`), including input deep enough to hit the JSON decoder's own limit; the same reader takes the prediction files and `load_report`'s file. It renders `OUT_DIR/ablation-report.json` and `OUT_DIR/ablation-report.md` before writing either, writes both as UTF-8 to temporary names in `OUT_DIR` with the mode a plain `open()` would give (`0o666 & ~umask`), moves any earlier report files in `OUT_DIR` to backup names there, renames both new files into place and deletes the backups (§17). On any problem it prints `<path>: <message>` lines and exits 1, with `OUT_DIR` as it was: no file from the failed run, an earlier report restored byte for byte, no temporary or backup file left, and every folder this run created for `OUT_DIR` (the folder and any parents its `mkdir` made, never one that already existed, also when the path goes through `..`) removed again while empty, also when making `OUT_DIR` itself fails after its parents were made. If restoring a backup also fails, for example because renames into the final names keep failing or another process made a directory there, that backup is kept, never deleted, and the message ends `earlier report kept at <backup path>` so it can be renamed back by hand. A final path that exists but isn't a regular file (a directory, a dangling symlink or a FIFO) is refused before anything is written. This holds under any locale, including an ASCII one with UTF-8 mode off. Problems include a schema error, a NaN, an integer beyond ±2^53, unreadable or too deeply nested JSON, a missing file or a bad prediction file. The manifest is checked first, then the records, then the links between them, then the prediction files.

## General rules

- **Records are closed.** A missing or unknown key is an error, so a renamed field fails instead of being dropped.
- **Types are strict.**
  - Integers must be JSON integers within ±2^53; `true` and `1.0` are not integers.
  - Numbers must be finite and within ±2^53. NaN, infinity and huge integers raise `SchemaError`, never another exception.
  - Booleans must be `true` or `false`.
- **Formats:**
  - **ids** (dataset, run): `[A-Za-z0-9][A-Za-z0-9._+-]*`. No `/` or spaces.
  - **block ids** (§8): `[A-Za-z][A-Za-z0-9_-]*`, such as `K`, `C` or `noise`.
  - **contrast keys:** `<family>:<learner>:<split>:<baseline>-><treatment>:<endpoint>` (§12). Every part is non-empty and holds no whitespace, `/`, `:` or backtick, and the arms hold no `->`. See [Contrast keys](#contrast-keys).
  - **feature columns:** no backtick (§10).
  - **hashes:** 64 lowercase hex characters.
  - **`source_commit`:** 40 or 64 lowercase hex characters.
  - **dates:** `YYYY-MM-DD`.
- **Canonical JSON (§7).** `schema.canonical` is `json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False)`, and `schema.sha256` hashes its UTF-8 bytes. Non-string keys are refused.
- **Errors are collected.** A validator reports every problem it finds, not just the first. Cross-field and cross-record checks run only after the fields themselves pass.

## Families (§1, §9, §10)

| Family | Holds | Endpoint | Adjustment | Partition |
|---|---|---|---|---|
| A | pooled records only: the confirmatory nested-chain increments | `r2_ll` | Holm within A | confirmatory only |
| B | pooled records only: leave-one-block-out, plus C+K vs C | `r2_ll` | Holm within B | confirmatory only |
| C | per-dataset contrasts of the A and B comparisons | `log_loss` | BH across all of C | confirmatory only |
| D | exploratory, per dataset or pooled: random and source splits, other learners, secondary endpoints, sensitivity | any | none (no p-values) | any |

Per-dataset contrasts are C or D, and pooled records are A, B or D (§12). Only D appears in both, and D carries no p-values, so nothing is ever adjusted across the two kinds.

## Contrast keys

The contrast key, `<family>:<learner>:<split>:<baseline>-><treatment>:<endpoint>`, is checked against its record:

- A contrast's key must agree with the record's `family`, `learner`, `split`, `baseline`, `treatment` and `endpoint`.
- A pooled record's `family` must be the family its key implies, even when every dataset is excluded (I3).
- A planned id's key gives the family of a contrast that was never produced. So a development dataset planned in family C is refused, as is a per-dataset plan in A or B.

§12 adopts this grammar, J's format: the id is `<dataset>/<key>`, every part is non-empty, no part contains `:`, `/`, a backtick or whitespace, and `baseline` and `treatment` don't contain `->`.

## Prediction row (§4)

One row per (dataset, arm, learner, split, repeat, sample_id). Every metric and contrast is recomputed from these rows. The field names are J's.

| Field | Type |
|---|---|
| `dataset` | id |
| `arm`, `learner`, `split`, `sample_id` | non-empty string |
| `repeat`, `fold` | integer ≥ 0 |
| `y` | integer 0 or 1 |
| `score` | finite number in [0, 1] |
| `called` | boolean: `score >= 0.5`, or false for a no-call row |
| `cluster` | non-empty string or integer ≥ 0 |
| `fold_plan_sha256` | hash |

A `called` row must have `score >= 0.5`. A no-call row, whose score is the training prevalence, has `called: false` whatever its score.

`validate_predictions(rows, manifest=None)` checks a whole set:

- the key is unique;
- a sample has one label and one cluster in its dataset;
- a sample has one fold per (dataset, split, repeat), across arms and learners;
- there is one `fold_plan_sha256` per (dataset, split, repeat), identical across arms and learners;
- with a manifest, every row belongs to an `ok` run, `repeat` is less than the run's `repeats`, and the row's plan is the run's plan for that repeat.

**In the build (I4).** With `prediction_root`, `build_report` reads every `ok` run's `predictions` file (JSON Lines, a path relative to the root) and checks:

- the path stays inside the root;
- each row is a valid prediction row and belongs to an `ok` run that names that file;
- each run has rows for every repeat;
- the rows pass `validate_predictions` as one set against the manifest.

The report records each file's sha256 in `predictions_sha256`. That field is null when no root is given.

## Run manifest

| Field | Type |
|---|---|
| `schema_version` | `1` |
| `protocol_sha256`, `binary_sha256` | hash |
| `source_commit` | commit hash |
| `blocks` | list of block specs |
| `datasets` | non-empty list of datasets |
| `missingness` | list |
| `group_proxy_audit` | list |
| `runs` | non-empty list of runs |
| `compute` | object |

**Block spec (§5; H's `BlockSpec` serde JSON validates as is):**

- `id`: block id
- `kind`: one of H's kinds, `CuratedExact`, `Hierarchical`, `Presence`, `Count`, `Context` or `Anonymous`
- `version` and `source_sha256`: strings
- `priority`: integer 0–5, collapse order C=0 < H=1 < P=2 < N=3 < X=4 < K=5. It must match the kind, from CuratedExact=0 to Anonymous=5.
- `screen_cap`: integer ≥ 1 or null (K: 4000)
- `catalog_date`: date or null
- `circularity_risk`: boolean
- `unavailable_reason`: non-empty string or null

An available block (`unavailable_reason` null) needs a non-empty `version` and a 64-hex `source_sha256`. An unavailable block may have `""` for both; a non-empty `source_sha256` must still be a hash.

L's synthetic noise block (§8) has:

- `id` `noise`, `kind` `Anonymous`, `priority` 5 and `screen_cap` null;
- `version` `synthetic-v1`;
- `source_sha256` the sha256 of its canonical column JSON.

**Dataset:**

- `id`
- `species`
- `partition`: `development` or `confirmatory`
- `n_samples`: integer ≥ 0
- `group_counts`: an object from grouping name to integer ≥ 0, such as `{"lineage": 12}`

**Missingness (§6, amended by §10):** one entry per (dataset, block). In a manifest each entry carries `dataset`. The generator's per-cohort output may leave it out, and the manifest builder adds it.

- `rate_R` and `rate_S` in [0, 1];
- `rate_diff_p`: the two-sided p-value of the class missing-rate difference, in [0, 1], or null when not computed;
- `indicator_auroc`: the raw AUROC of the missing-indicator count predicting y = 1, in [0, 1] or null;
- `indicator_auroc_oriented`: `max(a, 1 - a)`, or null with the raw value;
- `flag`: boolean.

**The flag rule (§8).** `flag` must equal (`indicator_auroc_oriented` > 0.60) or (`rate_diff_p` is not null and < 0.01). Both thresholds are strict. A flag that differs either way is refused.

**Group-proxy audit:** one entry per (dataset, block, metric), with `value` (a number or null) and a boolean `flag`.

**Run:**

- `id`, `dataset` and `arm`
- `blocks`: block ids
- `learner` and `split`
- `repeats`: integer ≥ 1
- `fold_plan_sha256`: a list with one hash per repeat, or null
- `status`: `ok` or `unavailable`
- `reason`
- `metrics`: name to number or null
- `fit_count`: integer ≥ 0
- `seconds`: number ≥ 0
- `predictions`: a path, or null

An `ok` run:

- has a null `reason`;
- has a `predictions` path and a plan for every repeat;
- has at least one fit;
- uses no unavailable block.

An `unavailable` run must give its reason.

**Compute:** `planned_fit_units` (integer ≥ 0), `actual_seconds` (number ≥ 0) and `tiers_completed` (unique strings).

**Cross-checks:**

- Block, dataset and run ids are unique.
- (dataset, arm, learner, split) is unique.
- Every reference names a known dataset or block.
- An arm has the same blocks in every run.
- For each (dataset, split, repeat), every run carries the same plan, so contrasts compare paired predictions.

## Contrast record (§1)

J emits it and K validates it. It has these 24 fields, closed:

| Field | Type |
|---|---|
| `id` | `<dataset id>/<contrast key>`, such as `ecoli__amp/C:lr:lineage:base->base+C:log_loss` |
| `dataset` | the id's prefix, repeated; K checks they agree |
| `family` | `C` or `D` (§10) |
| `baseline`, `treatment` | arm names, as in the key |
| `learner`, `split`, `endpoint` | as in the key; `endpoint` is `log_loss` in family C (§9) |
| `estimate` | finite number, or null when undefined |
| `ci` | the 95% interval `[lo, hi]`, or null when undefined |
| `ci_tost` | the 90% interval `[lo, hi]`, nested in `ci`, or null when undefined |
| `ci_method` | exactly `cluster_jackknife_t` or `none` (§12) |
| `unit` | non-empty string |
| `n_units` | integer ≥ 1 |
| `few_clusters` | boolean, exactly `n_units < 20` (§13) |
| `p`, `p_adj` | number in [0, 1], or null |
| `adjust` | `bh` for C, `none` for D |
| `margin` | number > 0 |
| `decision` | `improves`, `worse` or `inconclusive`; `equivalent` is refused (§15) |
| `claim_level` | always `predictive_increment` |
| `undefined_resamples` | integer ≥ 0 |
| `undefined_reason` | non-empty string, or null when defined |
| `interval_label` | non-empty string |

**`ci_method` (§10 and §12, closes the re-check's L3).** The method is checked by equality against the closed set `schema.CI_METHODS`, which is exactly `cluster_jackknife_t` and `none`, never by a substring test. `none` is allowed only on an undefined record. A new method needs calibration evidence and a contract bump.

**Null shapes:**

- **Undefined endpoint:** `undefined_reason` is set. `estimate`, `ci`, `ci_tost`, `p` and `p_adj` are all null, and `decision` is `inconclusive`.
- **Defined endpoint:** `estimate`, `ci` and `ci_tost` are set, and the intervals nest: `ci[0] ≤ ci_tost[0] ≤ ci_tost[1] ≤ ci[1]` (§10, closes L2).
- **Family D:** `p` and `p_adj` are null.
- **Family C (defined):** `p` and `p_adj` are set, and `p_adj ≥ p` exactly, with no float slack (I5).

Equality of `p_adj` and `p` is allowed under Holm and BH, because both leave the largest p-value unchanged (§1).

**Undefined replicates (§13, refining §9).** The cluster jackknife needs every leave-one-cluster-out replicate, so there are `n_units` replicates and any undefined one makes the contrast undefined. K checks:
- a defined record has `undefined_resamples` 0;
- J's three reasons agree with the counts:
  - `one cluster: a jackknife interval needs at least 2` needs `n_units` 1;
  - `more than 10% of resamples undefined` needs `undefined_resamples > 0.1 × n_units`;
  - `<u> of <k> leave-one-cluster-out estimates undefined; the jackknife needs all` needs u = `undefined_resamples`, k = `n_units` and 0 < u ≤ 0.1 × k.

Any other reason, such as an endpoint undefined on the full data, is that endpoint's own text and is not checked against the counts.

**Multiplicity (§13, §17).** K recomputes `p_adj` for families A, B and C from the records and both planned lists: Holm over the planned pooled keys of A, and of B; BH over the planned contrast ids of C. A planned hypothesis with no record, or with an undefined record (p null), enters at p = 1. A record whose `p_adj` differs from the recomputed value by more than 1e-12 (`schema.ADJUST_TOLERANCE`) is refused, on build and on reload. `schema.holm` and `schema.bh` are J's `stats.holm` and `stats.bh` in plain Python, with the same arithmetic. The planned lists are still trusted: a hypothesis removed from both a plan and the records changes the recomputed values without trace (see the trusted-plan limit below).

**The decision rule (§1 without its `equivalent` branch, §15).** `schema.expected_decision(ci, p_adj)` matches J's `decide(ci, p_adj)`. K re-derives every decision on validation and on reload, and refuses a record whose `decision` differs:

1. `improves` when `ci[0] > 0` and `p_adj < 0.05`.
2. `worse` when `ci[1] < 0` and `p_adj < 0.05`.
3. Otherwise `inconclusive`.

Family D has no p-value, so its decision is always `inconclusive`.

**No `equivalent` (§15).** J's re-check measured P(`equivalent` | true difference = margin) above 0.05 at every cluster floor it tried, so `equivalent` is withdrawn. `schema.DECISIONS` is `improves`, `worse` and `inconclusive`, and a record that says `equivalent` is refused on validation and on reload. `ci_tost` stays required with its §1 meaning; it decides nothing and feeds only the descriptive within-margin note below. The §12 floors (`MIN_EQUIVALENCE_UNITS`) are gone from K and from J. Restoring an equivalence decision needs a contract bump backed by a calibration study.

Set `GRM_TASK_J_STATS` to J's `stats.py` and the suite also checks K's mirrored constants against J's: the decision set, the few-clusters and undefined-share constants, the method set, families, adjustments, primary endpoints and record keys, and that `expected_decision` and J's `decide` agree over a grid of intervals and adjusted p-values. It also checks that J no longer has `MIN_EQUIVALENCE_UNITS`.

## Pooled record (§2)

Pools one comparison across datasets, with species as the unit. It has these 16 fields:

| Field | Type |
|---|---|
| `contrast` | contrast key; `r2_ll` for A and B |
| `family` | `A`, `B` or `D`, the family its key implies |
| `unit` | `species` |
| `n_units` | contributing species, integer ≥ 0 |
| `n_datasets` | every planned dataset of the comparison, excluded ones included |
| `n_excluded` | integer, equal to `len(excluded)` |
| `excluded` | list of `{dataset, reason}`, one entry per dataset |
| `estimate`, `ci`, `ci_tost`, `p`, `p_adj` | as for a contrast; all null when every dataset is excluded |
| `adjust`, `margin` | `holm` for A and B, `none` for D; margin > 0 |
| `decision` | by the same rule: `improves`, `worse` or `inconclusive` (§15) |
| `interval_label` | non-empty string |

Pooled intervals use the species jackknife-t (§12, replacing §2's "pooled bootstrap"); the pooled record has no `ci_method` field, so K doesn't check the method.

**What a pooled record draws on (§13).** A pooled record links to its per-dataset inputs by learner, split, baseline, treatment and endpoint, not by family. An A or B pool on `r2_ll` therefore draws on the family D `r2_ll` contrasts of its comparison, as J's `pool()` does.

**Its denominator** is the planned datasets of those inputs. For A and B that means the confirmatory ones only, since development datasets never enter an A or B denominator.

**Confirmatory tests.** The family C `log_loss` contrasts of the same comparison are the per-dataset confirmatory tests. The summary lists them, with their decisions, beside each A and B line.

**Few species (§13, §17).** With S contributing species, the exact sign-flip test can't give a pooled p below 2/2^S. K enforces it:

- a pooled record with fewer than 2 contributing species (`schema.MIN_POOLED_SPECIES`) is undefined, with null estimate, interval and p and every dataset excluded with its reason; K refuses a defined one, in any family;
- a defined pooled record's p must be at least 2/2^S − 1e-12 (`schema.p_floor`, `schema.FLOOR_TOLERANCE`).

So a pooled `improves` or `worse` with S < 6 is refused, since its p would have to be below 0.05. Any pooled line with S < 6 still says a pooled p below 0.05 isn't reachable, and gives that smallest p.

## Feature record

| Field | Type |
|---|---|
| `dataset` | id |
| `column` | non-empty string, no backtick |
| `blocks` | non-empty list of block ids |
| `selection_freq` | number in [0, 1] |
| `sign_stable` | boolean |
| `claim_level` | see [Claim levels](#claim-levels) |
| `evidence_refs` | list of unique strings, may be empty |
| `annotations` | optional; present only after `annotate_claim` |

## Claim levels

`predictive_increment` < `association_stable` < `mechanism_candidate` < `causal_supported`.

- **Contrasts always carry `predictive_increment`** (§1, §8). They take no manual level, `evidence_refs` or `annotations`.
- **Features from the pipeline** may carry `predictive_increment` or `association_stable`.
- **`association_stable` and above** need `sign_stable`.
- **`causal_supported`** is always refused without a non-empty `evidence_refs`.
- **The manual path is `annotate_claim(record, "feature", level, evidence_refs, by, date)`.** It is the only way to set `mechanism_candidate` or `causal_supported`.
  - `evidence_refs` must be a list of strings, and the result must have at least one.
  - It appends `{by, date, level, previous_level}` to `annotations`.
  - The chain is checked:
    - the first `previous_level` is a pipeline level;
    - each later `previous_level` equals the level the annotation before set;
    - the last `level` equals `claim_level`;
    - dates never decrease and are not in the future.
  - `annotations` on a non-manual level are refused.

**Limit:** a JSON file can carry a hand-written `annotations` chain that passes these checks. The schema can't tell it from the manual path. Who may annotate is a process control, not a schema one.

## The report

`build_report(manifest, contrasts, pooled, features, planned_contrasts, planned_pooled, not_produced=None, out_dir=None, prediction_root=None)` validates every input, then runs the cross-record checks below.

- `planned_contrasts` is the protocol's full list of per-dataset contrast ids (`<dataset>/<key>`, family C or D).
- `planned_pooled` is the protocol's full list of pooled keys (a bare §12 key, family A, B or D; A and B on `r2_ll`) (§17). It is required, so a caller can't skip the check by leaving it out.
- `not_produced` optionally maps a planned id or planned pooled key with no record to its reason. J's `finalize` `missing` map is the source of these reasons and can be passed as it is.
- `prediction_root`, when given, makes the build check the prediction files (I4).

**`planned_contrasts` and `planned_pooled` are trusted input from the protocol (I2, §17).** K checks every record, denominator and adjustment against them. But a list with a hypothesis removed, together with its record and its place in each denominator, is indistinguishable from the protocol's own list, and the recomputed `p_adj` follows the shorter list. The only trace is `inputs_sha256.planned_contrasts` or `inputs_sha256.planned_pooled`. Compare those hashes with the protocol's lists to catch it.

It returns:

- `schema` and `schema_version`;
- `inputs_sha256`: the canonical hash of each input;
- `predictions_sha256`: each checked prediction file's sha256, or null;
- the `manifest`, `planned_contrasts`, `planned_pooled` and `not_produced_reasons` inputs;
- `attempts`: runs, ok and unavailable;
- `contrast_counts`: planned, produced, defined, undefined and not produced, for per-dataset contrasts;
- `pooled_counts`: planned, produced and not produced, for pooled keys (§17);
- `not_produced`: each planned id with no record, and its reason;
- `pooled_not_produced`: each planned pooled key with no record, as `{contrast, reason}`;
- `undefined_contrasts`: each undefined record and its reason;
- `unavailable_runs` and `unavailable_blocks`, each with its reason;
- `exclusions`: every pooled exclusion, with its reason;
- the `contrasts`, `pooled` and `features` records.

With `out_dir`, it also writes the JSON and Markdown files, as `report.py` does: any earlier report there is moved to backup names first and restored byte for byte if writing fails, and the backups are deleted on success. A write failure raises `SchemaError` with one `<path>: <message>` entry, leaving no new file or temporary file; a backup is left only when restoring it failed, and the message names it. `build_report` and `validate_report` check the depth of every argument's dicts, lists and tuples, dict keys included, before copying or walking it (`build_report: <argument> nested too deeply`, `report: nested too deeply`). A deep object built from another container, such as a `deque`, a `frozenset` or a `UserList`, is refused by `build_report`'s deep copy instead (`build_report: input nested too deeply`). A non-string key in an object is named in short form (`reprlib`), and its value isn't checked. Either way a deep Python object raises `SchemaError`, never `RecursionError`, from `build_report`, `validate_report` (`predictions_sha256` included) and `load_report`.

**Reload.** `validate_report(report, prediction_root=None)` and `load_report(path, prediction_root=None)` rebuild the report from its own inputs. They refuse it if any section differs, compared as canonical JSON, so `true` is not `1` and `16.0` is not `16`. An edited decision, count, interval or exclusion therefore fails on reload. With a prediction root the files are re-hashed. Without one, a recorded `predictions_sha256` is checked only for shape: one hash per `ok` run's file.

Reload reads only the JSON; it doesn't compare the Markdown file with `render_markdown(report)`, so regenerate the summary from a reloaded report rather than trusting a `.md` on its own.

**Limit:** an edit applied consistently to an input and to every derived section passes, because the report re-derives only from itself. Since §17, that no longer covers a forged `p_adj`, which the recomputation refuses. Its `inputs_sha256` then differs from the canonical hash of the original input files; compare the two to detect it.

**Cross-record checks:**

- **Every planned contrast is accounted for.**
  - A planned id's family is C or D; a development dataset can't be planned in C.
  - Each contrast id must be planned.
  - Each planned id without a record needs a reason. The reason comes from `not_produced`, or else from the exclusion of that dataset in the pooled record that draws on its comparison. With no reason, the build is refused.
  - A `not_produced` entry for an id or key that is unplanned or has a record is refused.
- **Every planned pooled key is accounted for (§17).**
  - Planned pooled keys are unique, family A, B or D, and A and B use `r2_ll`.
  - Each pooled record's key must be in `planned_pooled`.
  - Each planned pooled key without a record needs a `not_produced` reason. With none, the build is refused.
- **Contrasts:** ids are unique, each dataset is in the manifest, and both arms need an `ok` run with the contrast's learner and split.
- **Pooled records:**
  - There is one pooled record per key, and its comparison must have planned per-dataset contrasts.
  - `n_datasets` equals the planned datasets of the comparison, and every excluded dataset is one of them.
  - Every included dataset has a per-dataset record. A dataset without one must be excluded with a reason.
  - An undefined input contrast must be excluded.
  - `n_units` equals the distinct species of the included datasets, neither more nor fewer.
  - Development datasets are refused in A and B denominators.
- Features name known datasets and blocks, and (dataset, column) is unique.
- **The adjustment is recomputed** in families A to C, as described under Multiplicity.

## Summary wording rules

The Markdown summary prints, in this order:

1. the hashes, including the prediction files or "not checked";
2. every attempt, with each unavailable run and block and its reason;
3. the planned-contrast counts, each contrast not produced with its reason, and each undefined contrast with its reason; then the pooled counts, "N pooled records planned: P produced, Q not produced";
4. every exclusion;
5. the pooled records (families A, B, D), each with:
   - its family, design, adjustment, margin and denominator, as "N datasets, k contributing, e excluded" ("1 dataset" in the singular);
   - on every line with a p (A and B): "Pooled p tests 'no effect in any species'. It is not a test that the average species effect is nonzero." (§15);
   - the few-species note when S < 6;
   - for A and B, the family C `log_loss` confirmatory tests and their decisions;
   - after the records of each family, every planned pooled key of that family that wasn't produced: "- `key` (family F: design): not produced: reason.";
6. the per-dataset contrasts (families C, D);
7. for every defined record, the 95% CI and the 90% TOST interval, followed by its `interval_label`;
8. the within-margin note on qualifying per-dataset contrasts (see wording rule 2).

**Escaping.** Free text (reasons, species, arms, labels, names) is escaped: Markdown punctuation gets a backslash, including every `#`, and line breaks become spaces. Ids and keys go in code spans fenced with one more backtick than the longest backtick run inside (the re-check's L1). Contrast keys and feature columns can't hold backticks anyway (§10). Text from a record can't start a line, add a heading, open a link or forge a reading rule.

Its fixed reading rules close the summary, and the tests hold its text to them:

1. **"Not significant" never becomes "no value".** An inconclusive decision is worded by its cause:
   - when the 95% interval includes zero: "not significant: the 95% interval includes zero";
   - when it excludes zero but `p_adj` ≥ 0.05: "inconclusive: the 95% interval excludes zero, but the adjusted p-value (…) does not meet 0.05";
   - in family D: "… but exploratory family D decides no direction".

   The phrase "the interval does not settle the difference" is not used anywhere in the summary.
2. **The 90% interval inside the margin is descriptive only (§15).** A defined per-dataset contrast with `unit == "cluster"`, `n_units >= 60`, decision `inconclusive`, and `-margin < ci_tost[0]` and `ci_tost[1] < margin` (strict) gets this exact sentence after its decision: "The 90% interval lies within ±margin. This is descriptive only: GRM's calibration did not hold the false-equivalence rate at 5%, so it is not an equivalence finding." The sentence is printed verbatim; the margin's value is already on the line. It never appears on a species-pooled record, and there is no record field for it.
3. **A small significant effect is named (§9).** An `improves` or `worse` whose |estimate| is below its margin reads "…, significant, but smaller than the margin of m". It comes from §9's report-wording ruling, which §15 leaves in place; it compares the estimate with the margin only on `improves` and `worse` lines and implies no equivalence.
4. **A small leave-one-block-out drop means redundant, not unimportant.** A line whose treatment arm's blocks (from the manifest's runs) are a strict subset of the baseline's says the block is "redundant given the others, not unimportant", unless it is `worse`. This applies to any family, so it covers the family C per-dataset contrasts of a B comparison and the B pools. A B line comparing C+K with C adds a block, so it gets no such note.
5. **A block attribution is descriptive only.** The feature section opens with a fixed sentence that the tests hold word for word.
6. **Family D is exploratory.** It covers random and source splits, other learners, secondary endpoints and sensitivity analyses. Where it compares splits, the random-vs-lineage gap is a diagnostic, not causal biology.
7. **Every interval label is printed,** including on undefined records and all-excluded pools.
8. **A claim level is printed only beside an `improves` decision.** Feature lines print theirs, with any manual annotations.

The tests scan the summary, outside its reading rules, for phrasing these rules exclude. The scan covers:

- "no value", "useless", "important" and "does not settle";
- "unimportant" and "causal" without a preceding "not";
- any form of "explain";
- "mechanism";
- "no effect" and the like;
- "equivalen" in any form (CN-K-L2).

Two sentences are exempt, each once and only where it belongs (CN-K-L1): §15's pooled-p sentence on a pooled line (it quotes the sign-flip null "no effect in any species"), and §15's within-margin sentence on a per-dataset contrast line (it names equivalence). Reading rule 5 is outside the scan. Anywhere else, either sentence is flagged, so it can't appear on another kind of line or in the `Undefined:` list; the tests also hold that list to one "- `id`: reason" line per undefined contrast.

It runs over every family and decision shape, including the small-effect shape, undefined records and leave-one-block-out variants.
