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

`splits` has one entry per split that ran, in the order of the report's `splits`. A split with an `unavailable_reason` fitted no model and has no entry, so match entries by `scheme`. Inside an entry each learner has one model per fold, at the fold's id: entry `f` is the model fitted for the samples whose `test_fold` is `f`, on the other folds' samples.

Each model has the shape of the final model beside it in `models`:

| Learner | Fold model | Same shape as |
|---|---|---|
| `scm` | the kind and trade-off `p` its inner folds chose, and its rules in the order the fit added them | `models.scm` without `steps` (fold fits list no tied rules) |
| `cart` | the whole tree, grown to the depth its inner folds chose | `models.cart` |
| `rf` | the 50 heaviest features by importance, heaviest first | `models.rf_importance` |
| `gbt` | the 50 heaviest features by gain, heaviest first; an empty list for a model that kept no tree | `models.gbt_importance` |
| `lasso` | the scored refit: intercept, weights, penalty; `null` when the tuning path was empty, and no weights when the refit kept none | `models.lasso` |

The sizes agree with what the report already says: a fold's `scm` has `fold_complexity` rules, its `cart` has `fold_complexity` leaves and its `lasso` has `fold_complexity` weights. The scm's rules, the tree and the lasso refit are the whole model: read on the samples' presence patterns they give the fold's saved scores (the test `fold_models_give_each_folds_saved_scores_under_dataset_row_numbers` checks it), except a null fit. A lasso fold model that is `null` or has no weights, like a boosting fold with an empty list, scored each test sample of its fold at the fold's training prevalence (resistant training samples over training samples, from `test_fold` and `samples[].resistant`; docs/equal-scoring.md), not at the sigmoid of its intercept. The test `a_weightless_lasso_fold_scored_its_test_samples_at_its_training_prevalence` checks it on a fold whose refit kept no weight. The forest's and boosting's lists are the top of a ranking and do not rebuild a score.

### Feature numbers

A `feature` is a row of the dataset's matrix, the numbering `models` and the `features` cards use. A fold does not fit on that matrix directly. It fits on the rows that differ on its own training samples (`learn::training_candidates`), numbered from 0, and the fold's features are renamed back to dataset rows before they are kept. A planted fault that keeps the fold's own numbering fails two tests, `fold_models_are_written_in_the_final_models_shapes_under_dataset_rows` and `fold_models_give_each_folds_saved_scores_under_dataset_row_numbers`.

One consequence matters for counting agreement. Where several rows agree on a fold's training samples, the fold sees one candidate and names the lowest row. Two folds, or a fold and the final model, can therefore name different rows for patterns that differ only on samples one of them did not train on. The test `folds_name_the_lowest_row_among_patterns_equal_on_their_training_samples` shows a fold naming row 0 and the other fold row 1 for one pattern on the first fold's training samples. So:

- equal feature numbers are the strict reading of "the same feature", and they undercount agreement;
- the patterns settle it. `fold_models.features` gives the presence pattern of every feature a fold's `scm`, `cart` or `lasso` names, ascending by feature, whether or not the feature also has a card. `bits` is encoded as on a card: base64 of little-endian 64-bit words, sample `i` at bit `i % 64` of word `i / 64`. Two features are the same for fold `f` when their bits agree on every sample whose `test_fold` is not `f`.

The forest's and boosting's features get no pattern, to bound the size. Those two learners compare by feature number alone, unless the feature is also named by another fold model or has a card.

### Size

**Measured on real reports.** The nine `enterobacter_cloacae` reports of `data-rev14-274-20261112` were rebuilt on 14 Nov 2026 from the public inputs (270 genomes), once with the source before this change and once with it. The earlier source's reports equal the released ones once wall-clock seconds are removed, so each new report is the released report plus the two fields.

| | All nine | Smallest | Median | Largest |
|---|---:|---:|---:|---:|
| `fold_models`, bytes | 499,860 | 50,707 | 55,457 | 60,568 |
| `rule_groups`, bytes | 17,075 | 1,587 | 1,915 | 2,281 |
| both fields, bytes | 516,935 | 52,311 | 57,633 | 62,849 |
| both as a share of the released report | 0.94% | 0.36% | 1.41% | 37.96% |

The fields' bytes are the same in every rebuild. A rebuilt file's own size also moves by a few tens of bytes from run to run, because wall-clock seconds are written with a varying number of digits, so the cost is given as the fields' bytes and not as a difference of two file sizes. The growth is close to constant, 52 to 63 kB, while these reports run from 164 kB to 17 MB, so the share is small on a report with many feature cards and large on a small one. Of the 499,860 bytes of `fold_models`, the forest's lists are 220,691, boosting's 139,152, the patterns 47,648 (906 features), the trees 41,452, the lasso refits 40,293 and the scm models 9,238. The two ranked lists are 72% of the field.

**Estimated over the release.** No released report has the field, so for the other reports the cost is an estimate from each report's own shape: the fold models' sizes are exact (`fold_complexity`), the bytes of one rule, node, list entry and weight are measured on the report's final models, and the pattern length follows from the sample count. Two counts are not in a released report, so each is a range. The distinct features the folds name run from "every fold names the features of the fold that names most" to "no two fold models share a feature". A fold's boosting list runs from one entry to 50: in the nine real reports the fold lists held 0 to 50 entries beside final lists of 29 to 50. The real field lay inside the range in all nine, at 1.33 to 1.65 times the low end and 0.79 to 0.90 of the high end. The range was widened after seeing them, so that fit is not a test of it. The test came after, with the range fixed first: two large `pseudomonas_aeruginosa` reports rebuilt the same way (meropenem and tobramycin, 38.9 MB and 138.3 MB). Their fields are 106,649 and 106,609 bytes, inside the ranges 53,227 to 126,140 and 54,431 to 125,944, and the two files grew by 0.31% and 0.09%.

Over the 274 reports of `data-rev14-274-20261112` (1,014,481,242 bytes; median report 678,123 bytes), counted on 14 Nov 2026:

| Part | All reports | Median report |
|---|---:|---:|
| `scm` fold models | 0.4 MB | 1.3 kB |
| `cart` fold models | 2.9 MB | 7.1 kB |
| `lasso` fold models | 2.6 MB | 8.1 kB |
| `rf` lists | 6.6 MB | 24.6 kB |
| `gbt` lists | 0.1 to 6.0 MB | 0.5 to 23.0 kB |
| patterns of the scm, cart and lasso features | 3.5 to 21.3 MB | 8.1 to 50.6 kB |
| `fold_models` in all | 16.2 to 39.9 MB | 51.6 to 116.8 kB |
| `rule_groups` | 3.7 MB | 11.0 kB |
| both fields | 19.9 to 43.6 MB | 65.5 to 129.9 kB |

That is 2.0% to 4.3% of the release's bytes, because a few reports with very many feature cards hold most of them. On a report of its own the share is larger: 12.6% to 26.6% at the median, and 23.4% to 73.7% on the report where it is largest. Without the two ranked lists the fold models would cost 5.9 MB plus the patterns.

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

The per-hit answer is a new field on each sample: `rule_groups`, the indices into `catalog_groups` of the genome's groups whose own hit counts for the rule, in the order of `groups`. It is an empty list when none counts. It is absent when the drug has no rule, where `catalog_resistant` is null, and on reports written before the field. It is the set the outside-study kit already writes by name, under the same column name, in its `known-gene/calls.tsv`.

Under a rule with one requirement, `catalog_resistant` is true exactly when `rule_groups` is not empty. Under the two-requirement rule a listed group met one of the two, and the list does not say which. So a genome can list groups and still be called susceptible: on the release's `staphylococcus_aureus__sulfamethoxazole_trimethoprim`, rebuilt with the field on 15 Nov 2026, 388 of 1,389 genomes list a group and none is called resistant. To explain a call, check each requirement against the listed groups' own hits.

On a synthetic dataset where one group's alleles differ (`a_groups_flag_follows_its_first_genome_and_rule_groups_each_genomes_own_hit`), the row says the group does not count while 19 of its 38 genomes' hits do; with the first genome's allele swapped the row says it counts while the other 19 do not. `rule_groups` names the same 19 genomes either way.

On eleven real reports rebuilt with the field on 14 Nov 2026 (nine of `enterobacter_cloacae`, and meropenem and tobramycin of `pseudomonas_aeruginosa`; 2,437 samples, 22,789 hits), every sample's `rule_groups` equals the list worked out outside the tool from that genome's own calls in `counted.json`. A hit answers the rule differently from its group's row on 54 samples: in 40 hits the row says the group counts and the genome's allele does not, and in 14 the reverse. Only 20 of those samples have a call that a reader of the rows would get wrong, the 15 and 5 the release counts above include for the two `pseudomonas_aeruginosa` reports. The other 34 are the ones a released report cannot show.

## Schema and revision

Neither changes.

- **Report schema, 2.** Both fields are optional additions. A schema-2 reader that ignores them reads every number as before, which is the rule revision 14 followed for `called_metrics`. Readers that check exact report keys must accept `fold_models` at the top level and `rule_groups` on a sample. `scripts/check-results.py`, `scripts/metric-audit/audit.py`, the outside-study kit and the web app's validator read named keys and ignore others.
- **Evaluation revision, 14.** The revision changes when a fit, a score or a rule changes. Here the evaluation computes what it computed and keeps more of it. A report built from the same inputs before and after differs only by the two fields and by wall-clock seconds: its other keys keep their order and values (see Checks).

One consequence follows from not bumping the revision. `grm learn` reuses a saved report whose revision and input fingerprint match, so the saved revision-14 reports stay as they are and gain the fields only when a run recomputes them (`--force`, changed inputs, or the next revision). A missing field therefore means "not recorded", never "no fold models" or "no hit counts". The index (`index.json`) holds neither field and does not change.

## Checks

- `src/learn/eval.rs`: `fold_models_give_each_folds_saved_scores_under_dataset_row_numbers` (the kept scm, tree and lasso give each fold's saved scores when read on the dataset matrix, a lasso with no weights at the fold's training prevalence, equal a refit renamed through the fold's candidate rows, and match `fold_complexity`; the ranked lists equal the refit's importance by dataset row), `a_weightless_lasso_fold_scored_its_test_samples_at_its_training_prevalence`, `folds_name_the_lowest_row_among_patterns_equal_on_their_training_samples`, `a_split_that_did_not_run_keeps_no_fold_models`, `fold_models_are_written_in_the_final_models_shapes_under_dataset_rows`.
- `src/analysis.rs`: `fold_patterns_cover_the_features_scm_cart_and_lasso_name_and_no_others`, `reports_keep_each_available_splits_fold_models_and_their_patterns`, `a_groups_flag_follows_its_first_genome_and_rule_groups_each_genomes_own_hit`, `a_genome_can_list_rule_groups_and_be_called_susceptible_under_two_requirements`, `a_drug_without_a_rule_records_no_rule_groups`.
- `src/catalog.rs`: `a_hit_counts_by_its_own_alleles_class_and_subclass`.
- **Before and after.** The ignored test `analysis::tests::write_report_fields_fixture_root` writes a synthetic data root (150 made-up genomes, a tobramycin dataset, three made-up catalog groups) that `grm learn --species fixture` runs on without the network. Run it with a build of the earlier source and a build of this one, on two copies of the root, and compare the two reports key by key.

- **Real reports.** Eleven released reports were rebuilt from the public inputs with the source before this change and with it (see Size and `rule_groups`). In all eleven the earlier build's report equals the released one once wall-clock seconds are removed, and the new build's report differs from it by the two fields alone. They are 11 of 274 reports and 2 of 23 species.

## Limits

- The fold lists of the forest and of boosting are a top 50 without patterns. Agreement for those two learners is by feature number, which undercounts.
- A fold model's features have no annotation unless the feature also has a card. Cards are built for the final models' features only.
- The web app does not show either field yet. Its mechanism table still shows the group row's "Counts for the rule".
- `rule_groups` says which of a genome's groups count, not which requirement of a two-requirement rule a group met, and not which allele the genome has.
- The counts above are of the revision 14 release. They describe how often a group row's flag is not a genome's own answer. They say nothing about which genes explain resistance.
