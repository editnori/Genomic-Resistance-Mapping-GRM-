# No-genome baselines (`grm baselines`)

`grm baselines` writes `report/baselines/<id>.json` for every dataset that `report/index.json` lists. Each file says what two predictors that never read a genome reach on the report's own folds and isolates:

- **Constant MIC.** For a report with an MIC model, each MIC fold's training isolates give one constant: the lower median of their exact MICs (also their mode, written beside it). The constant is scored on the same isolates as the model, the fold's test isolates with an exact MIC, by the same rule: it agrees when it is within one doubling dilution of the measured value.
- **Lineage identity.** For each split with folds, a test genome scores the resistant share of its lineage among the training genomes of the outer fold that held it out, or that fold's training prevalence when no training genome shares its lineage. The scores are called and counted as a learner row's are.

The file is arithmetic on what the report already holds. The command fits nothing, reads no network and changes no report file, so a report tree built without it is byte for byte what it was, and a second run writes the same bytes. Each file records its report's SHA-256 and evaluation revision. The command refuses a report of another evaluation revision than its own, and computes every file before it writes any, so a refusal leaves the folder as it was. It also removes any file in `report/baselines/` that the index does not list.

```bash
unshare -rn ./target/release/grm --data data baselines
```

Run it after `grm report` and before the report tree's `SOURCE.sha256` is written, as for `grm intervals` ([docs/deployment.md](deployment.md)), so a publication ships the files. Without them the dashboard works as before.

## The MIC folds

The MIC regression has its own folds. `evaluate_mic` draws them with `training_folds` over the isolates with an MIC, stratified by reported label, with the report's fold count and seed (`eval_config`). They are stratified random folds with no lineage split, and they are not the random split's folds. The report does not store them, so the command rebuilds them and checks the rebuilt folds against what the fit left in the report. A contradiction stops the command, names the dataset, and writes no file.

**Prediction check (`prediction_check`).** This is the check a dataset is confirmed by. Each fold's held-out predictions come from one boosted model, so isolates with bit-identical predictions in `mic.predicted` lie in one fold. The check groups the isolates by the bits of their prediction and requires each group of two or more to lie inside one rebuilt fold.

The rule has one exception, for folds that keep no tree. The boosted model starts every isolate at the upper median of the exact MICs in its fold's fit part (the training fold minus the 15% early-stopping holdout), and early stopping can keep no tree. Such a fold predicts its start value for every isolate, and two such folds can share that value. The rule is:

- a *constant fold* is a fold whose isolates with an MIC all have one bit-identical prediction;
- a group that lies only in constant folds is left out, whether it spans one of them or several;
- a group that spans two or more folds and touches a fold that is not constant is a contradiction;
- the folds are `confirmed` when at least one group of two or more isolates lies inside one fold that is not constant.

The file writes `groups` (the groups counted), `isolates` (their isolates) and `constant_folds` (ascending fold numbers). A dataset whose only groups lie in constant folds is unconfirmed, and the dashboard marks its constant with †.

**Threshold check (`threshold_check`).** Each rebuilt fold's training isolates give the threshold `infer_breakpoint` fits, and that must equal the threshold the report stores for the fold (`mic.fold_breakpoints`, every field, bit for bit). The file writes how many folds it `checked` and whether any was. This check confirms nothing on its own. Where every stored threshold has agreement 1, any stratified draw with the same fold sizes matches it, so a wrong seed can pass it. A dataset that stores no fold threshold is still confirmed by its predictions. The dashboard does not show this check.

**The model's own count.** The command recounts the model's agreement from `mic.predicted` under the report's rule (`evaluate_mic` rounds a prediction half away from zero) and refuses a report where that count does not give `mic.essential_agreement` exactly.

**Copies in training (`exact_with_training_copy`).** The MIC folds are random, so they can put a k-mer-identical copy of a test isolate in its training fold. The file counts the exact-MIC test isolates for which another member of their `kset_duplicates` group has an MIC and lies in a different rebuilt fold. The count is written and changes no figure.

## The lineage baseline

The lineage baseline uses the report's own outer folds for each split (`test_fold`) and the report's scoring conventions, as `eval::pool_folds` scores a learner row:

- every genome of the split is scored by the fold that held it out, and the scores are pooled;
- a score of 0.5 or more is a resistant call (`learn::binary_metrics`), so a share of exactly 0.5 calls resistant, as a null fit's prevalence does;
- the training prevalence is `learn::train_prevalence`, the share a null fit scores;
- balanced accuracy is the mean of sensitivity and specificity, and a class with no genome counts 0.

The file writes the confusion counts, sensitivity, specificity, balanced accuracy and `unseen`, the genomes whose lineage has no genome in their fold's training set, with `unseen_share` of `n`. On a lineage split every test genome's lineage is out of training by construction, so every genome takes its fold's prevalence and the baseline has nothing of the lineage to read there.

## The file

```text
schema               1
dataset              the report's id
report_sha256        SHA-256 of the report file's bytes
evaluation_revision  the report's revision
mic                  null when the report has no MIC model; otherwise:
  folds, seed          the fold count and seed the folds were rebuilt with
  prediction_check     {confirmed, groups, isolates, constant_folds}
  threshold_check      {confirmed, checked}
  exact_with_training_copy
  n, n_exact           isolates with an MIC, and with an exact MIC
  model, median, mode  {agree, n, essential_agreement} on the n_exact isolates
  fold_constants       per fold: train_exact, test_exact, median, mode,
                       mode_tied (another value was as frequent; the lower is
                       taken), upper_median_differs (the two middle values
                       differ); median and mode are null when the fold's
                       training isolates have no exact MIC
lineage              one entry per split with folds:
  scheme, folds, n, tp, fn, tn, fp, sensitivity, specificity,
  balanced_accuracy, unseen, unseen_share
```

## On the dashboard

The dashboard reads `baselines/<id>.json` as optional. It shows a baseline only when the file's dataset and revision are the report's and its `report_sha256` equals the SHA-256 of the report bytes the page read. The overview, the datasets table and the methods page read the index, which has no hash, so there a file must also give the model's essential agreement equal to the index row's. A missing file, a file pinned to another report, or a file that does not parse shows nothing and no error, and the report renders as before.

- **Dataset report, MIC.** One sentence gives the constant's agreement, the model's, and the difference on the same exact MICs; the note under it gives the mode. The heatmap's in-band sentence carries the constant's exact-MIC count beside the model's.
- **Dataset report, learners.** The random split's header gives the lineage baseline's balanced accuracy, and each learner row's random-split cell gives its difference from it in points. On the lineage split no difference is shown, and one sentence gives the unseen share.
- **Overview.** The MIC finding and its strip give the median constant beside the model's median and the count of datasets where the model is not above its constant.
- **Datasets table.** The MIC column gives the constant under the model's figure, marked c.
- **Methods.** How the constant and the lineage baseline are made, and the counts at 90% beside them.

## On the revision-14 release

These figures come from running the command on the 274 reports of `rev14-274-20261112`, and every one is a description of these cross-validated folds of public data:

- **Folds.** The prediction check confirms the MIC folds on 134 of the 134 datasets with an MIC model, with at least 3 groups on each dataset (median 44) and 7,723 in all. The threshold check passes on 124: 10 datasets store no fold threshold. The MIC baseline review of 15 November 2026 counts 7,723 of 7,724 groups inside one fold under the true draw. Its remaining group is the start value shared by the three constant folds below, which the rule here leaves out.
- **The model's start constant.** On 3 of 670 folds the model predicts only its start constant, all three in `salmonella_enterica__cefotaxime` (folds 0, 1 and 4). There the model's essential agreement on exact MICs equals the constant's, 0.960526.
- **Copies in training.** 26 exact-MIC test isolates on 14 datasets have a k-mer-identical copy in their training fold, at most 4 on one dataset. The review counts the same 26 isolates on 14 datasets.
- **Constant MIC.** On the median dataset the constant agrees on 72.2% of exact MICs and the model on 84.8%. The model is above the constant by a median of 5.4 points and is not above it on 42 of 134 datasets. Against the mode (median 67.7%) the model is not above it on 38 of 134.
- **At 90%.** Read after the scores, as a description: the model reaches 90% on 49 of the 134 datasets, the constant on 33, and both on 27. 90% is an adapted research screen, from the essential-agreement limit in FDA's 2009 guidance for growth-based susceptibility test systems, so a constant can reach it with no genome read.
- **Lineage identity.** On the random split (274 datasets) the baseline's median balanced accuracy is 0.686, with a median 8.9% of genomes unseen; the best random-split learner is above it by a median of 20.8 points and at or below it on 3 of 274. On the lineage split (268 datasets) every genome is unseen; the baseline is exactly 0.5 on 236 and between 0.216 and 0.497 on 32.

## Beside `research.json`

`research.json` (README, "Research") holds a revision-5 MIC analysis with its own no-genome baseline, the median exact MIC of each test isolate's training fold, on folds that analysis refitted for the random and lineage splits. It is packaged from saved research inputs that are not published, it pins revision 5, no release holds it, and the site answers 404 for it. The baselines file is a separate thing: it is computed from the published dataset report alone, on the report's own MIC folds and outer folds, and pinned to that report's bytes. Neither reads the other, and this command leaves `research.json`, its loader and its dashboard block as they are. Where both exist for a report, the dashboard shows the constant in the MIC section's text and the research figures in their own block headed "No-genome baseline".

## Limits

- The constant and the model are compared on random folds. A test isolate's close relatives can be in its training fold, which a genome model can read and a constant cannot. How the model compares with a constant on unseen lineages is not measured: the MIC regression has no lineage split, and measuring one needs a refit.
- A baseline says what a predictor reaches without the genome on these folds. It does not say why a model is above or below it, and it is not a clinical figure.
- The constants are fitted and scored on exact MICs only. Censored MICs enter neither the constant nor the figures shown beside it.
