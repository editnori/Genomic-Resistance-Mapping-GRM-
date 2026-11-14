# What a report keeps about its fold models and its known-gene rule

Two optional fields of a dataset report, added after the revision 14 release:

- `fold_models`: what each outer fold's model reads, so agreement between folds, and between the folds and the final models, can be counted.
- `samples[].rule_groups`: for each genome, the catalog groups whose own hit counts for the drug's rule.

Neither changes a fit, a score or a metric. The report schema stays 2 and the evaluation revision stays 14 (see Schema and revision). No report of the `data-rev14-274-20261112` release has either field; a report gets them when `grm learn` next computes it.

A feature or a catalog group named here is a marker a model or a rule reads, and an association with the reported label. It is not a cause, and annotation does not establish mechanism.

## Fold models

### Why

A report's cross-validated scores come from one model per learner and outer fold. Until now the report kept those models' scores and sizes (`splits[].learners[].fold_complexity`) and only the final models' features (`models`). Whether five folds chose the same k-mer patterns, or the pattern the final model publishes, could not be read from a report. It needed a refit.

### Shape

`fold_models` sits at the top level, after `features`:

```json
"fold_models": {
  "splits": [
    {
      "scheme": "random",
      "scm":   [{"kind": "disjunction", "p": 1.0, "rules": [{"feature": 9, "presence": true}]}],
      "cart":  [{"nodes": [{"value": 0.5, "n": 10, "n_resistant": 5, "depth": 0,
                            "split": {"feature": 7, "gain": 0.25, "absent": 1, "present": 2}}]}],
      "rf":    [[{"feature": 12, "weight": 0.75}]],
      "gbt":   [[]],
      "lasso": [{"intercept": -1.0, "coef": [[9, 1.5]], "lambda": 0.125}]
    }
  ],
  "features": [{"feature": 7, "bits": "..."}, {"feature": 9, "bits": "..."}]
}
```

`splits` has one entry per split that ran, in the order of the report's `splits`. A split with an `unavailable_reason` fitted no model and has no entry, so match entries by `scheme`. Inside an entry each learner has one model per fold, at the fold's id: entry `f` is the model that scored the samples whose `test_fold` is `f`.

Each model has the shape of the final model beside it in `models`:

| Learner | Fold model | Same shape as |
|---|---|---|
| `scm` | the kind and trade-off `p` its inner folds chose, and its rules in the order the fit added them | `models.scm` without `steps` (fold fits list no tied rules) |
| `cart` | the whole tree, grown to the depth its inner folds chose | `models.cart` |
| `rf` | the 50 heaviest features by importance, heaviest first | `models.rf_importance` |
| `gbt` | the 50 heaviest features by gain, heaviest first; an empty list for a model that kept no tree | `models.gbt_importance` |
| `lasso` | the scored refit: intercept, weights, penalty; `null` when the tuning path was empty, and no weights when the refit kept none | `models.lasso` |

The sizes agree with what the report already says: a fold's `scm` has `fold_complexity` rules, its `cart` has `fold_complexity` leaves and its `lasso` has `fold_complexity` weights. The scm's rules, the tree and the lasso refit are the whole model: read on the samples' presence patterns they give the fold's saved scores (the test `fold_models_give_each_folds_saved_scores_under_dataset_row_numbers` checks it). The forest's and boosting's lists are the top of a ranking and do not rebuild a score.

### Feature numbers

A `feature` is a row of the dataset's matrix, the numbering `models` and the `features` cards use. A fold does not fit on that matrix directly. It fits on the rows that differ on its own training samples (`learn::training_candidates`), numbered from 0, and the fold's features are renamed back to dataset rows before they are kept. The planted fault "keep the fold's own numbering" fails four tests.

One consequence matters for counting agreement. Where several rows agree on a fold's training samples, the fold sees one candidate and names the lowest row. Two folds, or a fold and the final model, can therefore name different rows for patterns that differ only on samples one of them did not train on. The test `folds_name_the_lowest_row_among_patterns_equal_on_their_training_samples` shows a fold naming row 0 and the other fold row 1 for one pattern on the first fold's training samples. So:

- equal feature numbers are the strict reading of "the same feature", and they undercount agreement;
- the patterns settle it. `fold_models.features` gives the presence pattern of every feature a fold's `scm`, `cart` or `lasso` names, ascending by feature, whether or not the feature also has a card. `bits` is encoded as on a card: base64 of little-endian 64-bit words, sample `i` at bit `i % 64` of word `i / 64`. Two features are the same for fold `f` when their bits agree on every sample whose `test_fold` is not `f`.

The forest's and boosting's features get no pattern, to bound the size. Those two learners compare by feature number alone, unless the feature is also named by another fold model or has a card.

### Size

No released report has the field, so the cost is an estimate from each report's own shape: the fold models' sizes are exact (`fold_complexity`), the bytes of one rule, node, list entry and weight are measured on the report's final models, and the pattern length follows from the sample count. The number of distinct features the folds name is not in a released report, so the patterns are a range: from "every fold names the features of the fold that names most" to "no two fold models share a feature". On a synthetic report that has the field, the estimate is 13,671 bytes and the field is 13,532.

Over the 274 reports of `data-rev14-274-20261112` (1,014,481,242 bytes; median report 678,123 bytes), counted on 14 Nov 2026:

| Part | All reports | Median report |
|---|---:|---:|
| `scm` fold models | 0.4 MB | 1.3 kB |
| `cart` fold models | 2.9 MB | 7.1 kB |
| `lasso` fold models | 2.6 MB | 8.1 kB |
| `rf` lists | 6.6 MB | 24.6 kB |
| `gbt` lists | 5.7 to 6.0 MB | 22.9 to 23.0 kB |
| patterns of the scm, cart and lasso features | 3.5 to 21.3 MB | 8.1 to 50.6 kB |
| `fold_models` in all | 21.7 to 39.9 MB | 73.5 to 116.8 kB |
| `rule_groups` | 3.7 MB | 11.0 kB |
| both fields | 25.5 to 43.6 MB | 86.4 to 129.9 kB |

That is 2.5% to 4.3% of the release's bytes, because a few reports with very many feature cards hold most of them. On a report of its own the share is larger: 16.5% to 26.6% at the median, and 26.9% to 73.7% on the report where it is largest. The two ranked lists are more than half of the fixed part; without them the fold models would cost 5.9 MB plus the patterns.

## The known-gene rule: `counts_for_rule` and `rule_groups`

The known-gene rule is GRM's own reader of the AMRFinderPlus catalog (`src/catalog.rs`). A genome has at most one hit per catalog gene group, its best-matching allele, and the hit carries that allele's class and subclass. A drug's rule names classes or subclasses (`catalog::drug_rule`), and a hit counts for the rule when its allele's class or subclass holds one of the rule's tokens (`catalog::counts_for`). Under a rule with one requirement, a genome is called resistant when one of its hits counts. Trimethoprim with a sulfonamide has two requirements and needs a counting hit for each.

### What `counts_for_rule` means

`catalog.groups` has one row per group seen in the dataset. The alleles of one group can carry different classes and subclasses, so a group has no single answer. The row's `class`, `subclass` and `counts_for_rule` describe one allele: that of the first sample, in sample order, that has the group. `n_resistant` and `n_susceptible` count every sample that has the group, whatever its allele.

So `counts_for_rule` is not a statement about every genome with the group. Where a group's alleles differ, some genomes' hits count and others do not, and the flag follows whichever sample came first. The call and its metrics never read the flag: `samples[].catalog_resistant` and `catalog.metrics` come from each genome's own hits, and they are right as published.

In the 274 reports of `data-rev14-274-20261112` (catalog 2026-08-07.1), counted on 14 Nov 2026:

- **No flag at all.** 21 reports are of a drug without a rule. They have no `catalog` block, and their 17,715 samples have a null `catalog_resistant`. The other 253 reports hold 17,072 group rows, and each row's flag agrees with the row's own class and subclass.
- **A flag that is not the genome's answer.** A call rebuilt from the group rows differs from the stored call on 3,404 of 199,799 samples, in 36 reports of four species (*E. coli*, *K. pneumoniae*, *P. aeruginosa*, *S. enterica*): 770 samples are called resistant by their own allele while every row says its group does not count, and 2,634 the other way round. Every one of them has a group whose catalog alleles do not all answer the rule alike.
- **How far that can reach.** The catalog has 1,223 groups, and the alleles of 42 differ in class or subclass. 109 reports have at least one group whose alleles answer that report's rule differently, on 38,856 samples. The 3,404 are the ones where the difference changes a call. Where a genome also has another counting group, a released report cannot show which of its hits count.

### `rule_groups`

The flag is left as it is. Its rows keep their six keys and their values, because the frozen outside-study kit rebuilds the whole `catalog.groups` table and compares it with the report's (`scripts/outside-study/score.py`, `known_gene_replay`).

The per-hit answer is a new field on each sample: `rule_groups`, the indices into `catalog_groups` of the genome's groups whose own hit counts for the rule, in the order of `groups`. It is an empty list when none counts. It is absent when the drug has no rule, where `catalog_resistant` is null, and on reports written before the field.

Under a rule with one requirement, `catalog_resistant` is true exactly when `rule_groups` is not empty. Under the two-requirement rule a listed group met one of the two, and the list does not say which.

On a synthetic dataset where one group's alleles differ (`a_groups_flag_follows_its_first_genome_and_rule_groups_each_genomes_own_hit`), the row says the group does not count while 19 of its 38 genomes' hits do; with the first genome's allele swapped the row says it counts while the other 19 do not. `rule_groups` names the same 19 genomes either way.

## Schema and revision

Neither changes.

- **Report schema, 2.** Both fields are optional additions. A schema-2 reader that ignores them reads every number as before, which is the rule revision 14 followed for `called_metrics`. Readers that check exact report keys must accept `fold_models` at the top level and `rule_groups` on a sample. `scripts/check-results.py`, `scripts/metric-audit/audit.py`, the outside-study kit and the web app's validator read named keys and ignore others.
- **Evaluation revision, 14.** The revision changes when a fit, a score or a rule changes. Here the evaluation computes what it computed and keeps more of it. A report built from the same inputs before and after differs only by the two fields and by wall-clock seconds: its other keys keep their order and values (see Checks).

One consequence follows from not bumping the revision. `grm learn` reuses a saved report whose revision and input fingerprint match, so the saved revision-14 reports stay as they are and gain the fields only when a run recomputes them (`--force`, changed inputs, or the next revision). A missing field therefore means "not recorded", never "no fold models" or "no hit counts". The index (`index.json`) holds neither field and does not change.

## Checks

- `src/learn/eval.rs`: `fold_models_give_each_folds_saved_scores_under_dataset_row_numbers` (the kept scm, tree and lasso give each fold's saved scores when read on the dataset matrix, equal a refit renamed through the fold's candidate rows, and match `fold_complexity`; the ranked lists equal the refit's importance by dataset row), `folds_name_the_lowest_row_among_patterns_equal_on_their_training_samples`, `a_split_that_did_not_run_keeps_no_fold_models`, `fold_models_are_written_in_the_final_models_shapes_under_dataset_rows`.
- `src/analysis.rs`: `fold_patterns_cover_the_features_scm_cart_and_lasso_name_and_no_others`, `reports_keep_each_available_splits_fold_models_and_their_patterns`, `a_groups_flag_follows_its_first_genome_and_rule_groups_each_genomes_own_hit`, `a_drug_without_a_rule_records_no_rule_groups`.
- `src/catalog.rs`: `a_hit_counts_by_its_own_alleles_class_and_subclass`.
- **Before and after.** The ignored test `analysis::tests::write_report_fields_fixture_root` writes a synthetic data root (150 made-up genomes, a tobramycin dataset, three made-up catalog groups) that `grm learn --species fixture` runs on without the network. Run it with a build of the earlier source and a build of this one, on two copies of the root, and compare the two reports key by key.

## Limits

- The fold lists of the forest and of boosting are a top 50 without patterns. Agreement for those two learners is by feature number, which undercounts.
- A fold model's features have no annotation unless the feature also has a card. Cards are built for the final models' features only.
- The web app does not show either field yet. Its mechanism table still shows the group row's "Counts for the rule".
- `rule_groups` says which of a genome's groups count, not which requirement of a two-requirement rule a group met, and not which allele the genome has.
- The counts above are of the revision 14 release. They describe how often a group row's flag is not a genome's own answer. They say nothing about which genes explain resistance.
