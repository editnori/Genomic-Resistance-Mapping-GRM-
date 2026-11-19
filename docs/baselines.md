# No-genome baselines (`grm baselines`)

`grm baselines` writes `report/baselines/<id>.json` for every dataset that `report/index.json` lists and the command does not refuse, and one summary of the MIC constants, `report/baselines.json`, beside the index. Each file says what two predictors that never read a genome reach on the report's own folds and isolates:

- **Constant MIC.** For a report with an MIC model, each MIC fold's training isolates give one constant: the lower median of their exact MICs (also their mode, written beside it). The constant is scored on the same isolates as the model, the fold's test isolates with an exact MIC, by the same rule: it agrees when it is within one doubling dilution of the measured value.
- **Lineage identity.** For each split with folds, a test genome scores the resistant share of its lineage among the training genomes of the outer fold that held it out, or that fold's training prevalence when no training genome shares its lineage. The scores are called and counted as a learner row's are.

The file is arithmetic on what the report already holds. The command fits nothing, reads no network and changes no report file, so a report tree built without it is byte for byte what it was, and a second run writes the same bytes. Each file records its report's SHA-256 and evaluation revision.

The command refuses in two ways:

- **One dataset alone,** when its rebuilt MIC folds contradict the report (the prediction check or the threshold check below) or its saved predictions do not give the report's own figure. It prints the dataset with I and X, writes no file for it and removes the one it had, and writes the others. The run then ends with an error that gives the count of files written and of datasets refused.
- **The whole run, before any file is written,** for a `report/baselines` that is a symbolic link or not a directory (before anything is read), an index id that holds a path separator or `..`, or a report it cannot read or of another evaluation revision than its own. Every file is computed before any is written, so such a stop leaves the folder as it was.

Each file goes to a pending `<id>.json.<pid>.tmp` that is renamed into place, and a write that fails removes its pending file. A failure while writing, such as a full disk or a directory in a file's place, stops the run with the files written so far in place, each pinned to its own report. Every error names its path. After writing, the command removes the regular `*.json` and `*.tmp` files in `report/baselines/` that this run did not write, a refused dataset's old file among them; a link, a directory or any other entry stays. It holds the report directory's lock, as `grm report` does, so it waits while another writer runs.

```bash
unshare -rn ./target/release/grm --data data baselines
```

Run it after `grm report` and before the report tree's `SOURCE.sha256` is written, as for `grm intervals` ([docs/deployment.md](deployment.md)), so a publication ships the files and the summary. `grm report` touches neither: after a rebuild they describe the old build until `grm baselines` runs again, and the dashboard refuses them by their pins. The index's bytes include `generated_at`, so every run of `grm report` retires the summary.

## The MIC folds

The MIC regression has its own folds. `evaluate_mic` draws them with `training_folds` over the isolates with an MIC, stratified by reported label, with the report's fold count and seed (`eval_config`). They are stratified random folds with no lineage split, and they are not the random split's folds. The report does not store them, so the command rebuilds them and checks the rebuilt folds against what the fit left in the report. A contradiction refuses that dataset alone, as above.

**Prediction check (`prediction_check`).** This is the check a dataset is confirmed by. On the report's own folds, groups of isolates with bit-identical predictions in `mic.predicted` lie inside one fold, and a wrong rebuild splits them; correct folds can also give one value to isolates of two folds (a shared start value, or identical models), so a spanning group is evidence, not proof, of a wrong rebuild. The check groups the isolates with an MIC by the bits of their prediction and counts the groups of two or more:

- a *constant fold* is a fold whose isolates with an MIC all have one bit-identical prediction. The boosted model starts every isolate at the upper median of the exact MICs in its fold's fit part (the training fold minus the 15% early-stopping holdout), and where early stopping keeps no tree the fold predicts that start value for every isolate. A group that lies only in constant folds, inside one or across several, is left out;
- I, written `groups`, counts the other groups that lie inside one rebuilt fold, and `isolates` their isolates;
- X, written `spanning`, counts the other groups that span two or more rebuilt folds.

The rule:

- X = 0 and I ≥ 1: the file is written with `confirmed: true`;
- 0 < X ≤ I, or no group at all (I = X = 0): the file is written with `confirmed: false`, and the dashboard marks its constant with †;
- X > I: the dataset is refused alone and has no file.

The file writes `groups`, `isolates`, `spanning` and `constant_folds` (ascending fold numbers). The cut at X > I is read off one release (below), so it is a choice that separates the true draw from a wrong one there, not a proof.

**Threshold check (`threshold_check`).** Each rebuilt fold's training isolates give the threshold `infer_breakpoint` fits, and that must equal the threshold the report stores for the fold (`mic.fold_breakpoints`, every field, bit for bit). The file writes how many folds it `checked` and whether any was. This check confirms nothing on its own. Where every stored threshold has agreement 1, any stratified draw with the same fold sizes matches it, so a wrong seed can pass it. A dataset that stores no fold threshold is still confirmed by its predictions. A fold whose rebuilt threshold differs from the stored one refuses the dataset alone. The dashboard does not show this check.

**The model's own count.** The command recounts the model's agreement from `mic.predicted` under the report's rule (`evaluate_mic` rounds a prediction half away from zero) and refuses the dataset alone where that count does not give `mic.essential_agreement` exactly.

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
  prediction_check     {confirmed, groups, isolates, spanning, constant_folds}
  threshold_check      {confirmed, checked}
  exact_with_training_copy
  n, n_exact           isolates with an MIC, and with an exact MIC
  model, median, mode  {agree, n, essential_agreement} on the n_exact isolates
  fold_constants       per fold: train_exact, test_exact, median, mode,
                       mode_tied (another value was as frequent; the lower is
                       taken), upper_median_differs (the two middle values
                       differ); median and mode are log2 dilution steps
                       (the MIC is 2 to that power, in mg/L) and null when
                       the fold's training isolates have no exact MIC
lineage              one entry per split with folds:
  scheme, folds, n, tp, fn, tn, fp, sensitivity, specificity,
  balanced_accuracy, unseen, unseen_share
```

Units: counts are whole numbers of isolates or genomes; `essential_agreement` is `agree / n` and `unseen_share` is `unseen / n` (0 when `n` is 0); sensitivity, specificity and balanced accuracy are shares from 0 to 1, a class with no genome counting 0.

Form: every file, and the summary, is compact JSON with every object's keys in sorted order and no final newline. Numbers are written as the shortest decimal that reads back to the same double, and a whole float keeps its point (`1.0`). A reader in another language can check every share from its counts. To check the constants, the folds must be rebuilt with the repository's own `training_folds`, because the file holds no fold assignment.

## The summary (`report/baselines.json`)

```text
schema               1
datasets             one entry per file written with an MIC part and a
                     median constant, in index order:
  confirmed            the file's mic.prediction_check.confirmed
  constant             the file's mic.median.essential_agreement (the
                       lower median, never the mode)
  dataset              the file's dataset
  model                the file's mic.model.essential_agreement
  report_sha256        the file's report_sha256
evaluation_revision  the command's evaluation revision
index_sha256         SHA-256 of index.json's bytes as the command read them
```

The summary is what the overview, the datasets table, the methods page and the slides read: one request instead of one per dataset. Its name and place were chosen so that nothing else takes or removes it. It sits beside `index.json`, outside `baselines/`. A dataset's file is `baselines/<id>.json`, and the command refuses an id with a path separator, so no id can reach it. The command's removal reads `baselines/` only, and `grm report`'s prune reads only `datasets/` and `species/`, so neither deletes it. The deployment's tree check (`verify_tree`, [docs/deployment.md](deployment.md)) accepts it once `SOURCE.sha256` lists it and refuses the tree when it does not, as for any file.

## On the dashboard

The dashboard reads the files and the summary as optional, and checks each against what the page read:

- **The dataset report** reads its own `baselines/<id>.json`. It shows a baseline only when the file's dataset and revision are the report's and its `report_sha256` equals the SHA-256 of the report bytes the page read.
- **The overview, the datasets table, the methods page and the slides** read `baselines.json`. They take it only when its `index_sha256` equals the SHA-256 of the `index.json` bytes the page read, its revision is the index's, and every entry names a dataset the index lists with an MIC model and the same model figure. Anything else refuses the summary whole, and they show no constant.
- **Both pins need Web Crypto** (`crypto.subtle`), which a browser gives on HTTPS and on localhost. On plain HTTP from another address the page computes no digest and shows no baseline.

Where no file is published, the site shows no baseline. The MIC caption, the MIC section's note and two Methods entries carry the new wording either way, each index page requests the summary once, and each dataset report requests its own file. A file or a summary pinned to another build is refused.

- **Dataset report, MIC.** One sentence gives the constant's agreement, the model's, and the difference on the same exact MICs, and says on how many folds the model predicts only its start value where any does; the note under it gives the mode. A sentence on its own line under the heatmap's key carries the constant's share of exact MICs beside the model's.
- **Dataset report, learners.** The random split's header gives the lineage baseline's balanced accuracy, and each learner row's random-split cell gives its difference from it in points. On the lineage split no difference is shown, and one sentence gives the unseen share.
- **Overview.** The MIC finding and its strip give the median constant beside the model's median and the count of datasets where the model is not above its constant.
- **Datasets table.** The MIC column gives the constant under the model's figure, marked c.
- **Methods.** How the constant and the lineage baseline are made, and the counts at 90% beside them.

## On the revision-14 release

These figures come from running the command on the 274 reports of `rev14-274-20261112`, and every one is a description of these cross-validated folds of public data:

- **Folds.** The prediction check confirms the MIC folds on 134 of the 134 datasets with an MIC model: X = 0 on each, with I at least 3 on each dataset (median 44) and 7,723 in all. With the seed one higher, the command refuses all 134: X > I on each, I/(I+X) at most 0.385 (`neisseria_gonorrhoeae__erythromycin`, I = 5 and X = 8, median 0.113), and X at least 3. A synthetic report from the repository's own fit, whose two folds get one value from identical models, has I = 0 and X = 2: added to the release tree, it is refused alone and the 274 release files are written as before. The threshold check passes on 124: 10 datasets store no fold threshold. The MIC baseline review of 15 November 2026 counts 7,723 of 7,724 groups inside one fold under the true draw. Its remaining group is the start value shared by the three constant folds below, which the rule here leaves out.
- **The model's start constant.** On 3 of 670 folds the model predicts only its start constant, all three in `salmonella_enterica__cefotaxime` (folds 0, 1 and 4). There the model's essential agreement on exact MICs equals the constant's, 0.960526.
- **Copies in training.** 26 exact-MIC test isolates on 14 datasets have a k-mer-identical copy in their training fold, up to 4 on one dataset. The review counts the same 26 isolates on 14 datasets.
- **Constant MIC.** On the median dataset the constant agrees on 72.2% of exact MICs and the model on 84.8%. The model is above the constant by a median of 5.4 points and is not above it on 42 of 134 datasets. Against the mode (median 67.7%) the model is not above it on 38 of 134.
- **At 90%.** Read after the scores, as a description: the model reaches 90% on 49 of the 134 datasets, the constant on 33, and both on 27. 90% is an adapted research screen, from the essential-agreement limit in FDA's 2009 guidance for growth-based susceptibility test systems, so a constant can reach it with no genome read.
- **Lineage identity.** On the random split (274 datasets) the baseline's median balanced accuracy is 0.686, with a median 8.9% of genomes unseen; the best random-split learner is above it by a median of 20.8 points and at or below it on 3 of 274. On the lineage split (268 datasets) every genome is unseen; the baseline is exactly 0.5 on 236 and between 0.216 and 0.497 on 32.

## Beside `research.json`

`research.json` (README, "Research") holds a revision-5 MIC analysis with its own no-genome baseline, the median exact MIC of each test isolate's training fold, on folds that analysis refitted for the random and lineage splits. It is packaged from saved research inputs that are not published, it pins revision 5, no release holds it, and the site answers 404 for it. The baselines file is a separate thing: it is computed from the published dataset report alone, on the report's own MIC folds and outer folds, and pinned to that report's bytes. Neither reads the other, and this command leaves `research.json`, its loader and its dashboard block as they are. Where both exist for a report, the dashboard shows the constant in the MIC section's text and the research figures in their own block headed "No-genome baseline".

## Limits

- The constant and the model are compared on random folds. A test isolate's close relatives can be in its training fold, which a genome model can read and a constant cannot. How the model compares with a constant on unseen lineages is not measured: the MIC regression has no lineage split, and measuring one needs a refit.
- A baseline says what a predictor reaches without the genome on these folds. It does not say why a model is above or below it, and it is not a clinical figure.
- The constants are fitted and scored on exact MICs only. Censored MICs enter neither the constant nor the figures shown beside it.
- The fold check's cut, X > I, is read off one release. A wrong rebuild with fewer spanning groups than inside ones would be written and marked † rather than refused.
