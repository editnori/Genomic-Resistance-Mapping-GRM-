# Dense network

This page compares one dense feed-forward network with GRM's five learners, `scm`, `cart`, `rf`, `gbt` and `lasso`, on 24 public BV-BRC datasets of release `rev14-274-20261112`, on the release's own lineage and random folds. The network was fitted on 15 Nov 2026 at settings a plan fixed before its first result; the five were not fitted again, and their held-out calls are the ones the release saved. This is one experiment on public data. It establishes no gain and rules none out.

## Reading

Against the plan's comparator, the best of the five picked in hindsight per dataset, which favours the five, the median network minus comparator lineage balanced accuracy over the 24 datasets is -0.059 (-0.132 to -0.044, a 95% interval over the 12 species); the network is higher on 1 of the 24, and the rule gives "no gain established".

Against the best of the five in hindsight, the median difference in recall is -0.102 (-0.356 to 0.000) and in log loss +0.021 (-0.262 to +0.209), over the same 24 datasets.

Beside the rule, against each of the five, the median lineage balanced accuracy difference over the 24 datasets is +0.003 (-0.086 to +0.049) for `scm`, -0.012 (-0.052 to -0.001) for `cart`, -0.016 (-0.032 to +0.009) for `rf`, -0.012 (-0.035 to +0.069) for `gbt` and +0.015 (-0.007 to +0.071) for `lasso`; only `cart`'s interval lies below 0, and none above.

One network at one setting on one screen of 4,000 rows was fitted; no other setting, no other architecture and no transformer or sequence model was.

## Plan

The plan was fixed on 15 Nov 2026. Its file time is 12:19:40Z and its sha256 is `b49eaeee07a54b97828568418a810c36d28e70cc4b0b9ab5fdf1e2c513be127a`. The network's fitting script was last changed at 12:20:52Z, after the plan, and the first network result was written at 12:21:31Z. The order rests on file times alone: no file written before the first result records the plan's sha256.

The plan fixed the network, the datasets and the gates their exports had to pass, the comparator, the measures, the interval and a four-part rule on the lineage split. It reported the rest beside the rule and outside it: each of the five, the site's best (chosen on the random split), the mean of the four probability learners, the lineage baseline (lineage identity inside the release's own outer folds), the catalog rule where the release scores it, the random split, log loss and AUROC.

## Network

- **Layers.** Dense feed-forward: 4,000 inputs, 256 ReLU units, 64 ReLU units, 1 sigmoid output.
- **Inputs.** Per outer fold, the 4,000 rows of the forest's chi-square screen over that fold's training candidates, on its training samples; each sample a 0/1 vector over those rows.
- **Start.** He-normal, zero biases.
- **Loss.** Binary cross-entropy, unweighted, plus L2 decay 1e-4, biases left out.
- **Optimiser.** Adam, learning rate 1e-3, beta 0.9 and 0.999, epsilon 1e-8; minibatches of 64 in a fresh shuffle each epoch.
- **Stopping.** At most 300 epochs; log loss on the fold's boosting holdout, about 1/7 of the training samples, lineage-grouped on the lineage split; stop after 20 epochs without improvement, keep the best epoch.
- **Repeats.** 3 per fold, each with its own seed; a dataset's value is the mean over the 3.
- **Splits.** The release's lineage split and random split, 5 outer folds each; the five learners were not fitted again.

The screen ran on all of a fold's training samples, so the stopping holdout's labels took part in choosing the 4,000 rows, as the plan said in advance; no test label did. The network ran 720 fits in 1.838 CPU-hours, and the peak memory of a process was 299,972 kB.

## Datasets

The gap of a dataset is, on the lineage split, the balanced accuracy of calls right wherever any of the five is right, minus that of the best of the five. The plan took the 12 lineage datasets with the largest gap, 0.260 to 0.363, and drew 12 with seed 20261115 from the other 113 with a gap of 0.05 or more. They hold 83 to 1,957 genomes, in 12 species with 1 to 5 datasets each.

Every export passed both gates: its genomes, labels and outer folds equal the release's (gate A, 24 of 24), and so does its input fingerprint (gate B, 24 of 24).

The first 12 rows are the largest gaps, in the plan's order, and the last 12 the drawn. Balanced accuracy is on the lineage split, where the best of the five is picked per dataset.

| Dataset | Genomes | Resistant | Gap | Best of five | Its balanced accuracy | Network's |
|---|---|---|---|---|---|---|
| mycobacterium_tuberculosis__linezolid | 1,248 | 202 | 0.363 | `lasso` | 0.576 | 0.520 |
| salmonella_enterica__cefazolin | 185 | 119 | 0.357 | `scm` | 0.569 | 0.525 |
| enterobacter_cloacae__imipenem | 86 | 58 | 0.354 | `rf` | 0.559 | 0.521 |
| acinetobacter_baumannii__minocycline | 139 | 28 | 0.350 | `cart` | 0.614 | 0.562 |
| mycobacterium_tuberculosis__clofazimine | 1,303 | 193 | 0.344 | `cart` | 0.563 | 0.495 |
| mycobacterium_tuberculosis__delamanid | 1,248 | 448 | 0.343 | `gbt` | 0.575 | 0.538 |
| klebsiella_aerogenes__ertapenem | 83 | 48 | 0.308 | `rf` | 0.663 | 0.597 |
| mycobacterium_tuberculosis__pyrazinamide | 492 | 381 | 0.307 | `gbt` | 0.641 | 0.583 |
| klebsiella_oxytoca__ertapenem | 112 | 78 | 0.281 | `rf` | 0.698 | 0.689 |
| pseudomonas_aeruginosa__ceftolozane_tazobactam | 259 | 132 | 0.275 | `rf` | 0.668 | 0.641 |
| pseudomonas_aeruginosa__cefepime | 245 | 122 | 0.265 | `rf` | 0.612 | 0.649 |
| pseudomonas_aeruginosa__aztreonam | 235 | 114 | 0.260 | `rf` | 0.651 | 0.510 |
| escherichia_coli__nitrofurantoin | 482 | 46 | 0.055 | `lasso` | 0.792 | 0.668 |
| acinetobacter_baumannii__ceftazidime | 1,247 | 1,087 | 0.065 | `cart` | 0.891 | 0.740 |
| klebsiella_pneumoniae__colistin | 336 | 112 | 0.203 | `scm` | 0.734 | 0.538 |
| pseudomonas_aeruginosa__imipenem | 447 | 235 | 0.124 | `rf` | 0.708 | 0.584 |
| salmonella_enterica__sulfamethoxazole_trimethoprim | 1,145 | 296 | 0.186 | `gbt` | 0.725 | 0.715 |
| pseudomonas_aeruginosa__piperacillin_tazobactam | 210 | 87 | 0.210 | `cart` | 0.644 | 0.643 |
| neisseria_gonorrhoeae__penicillin | 1,354 | 888 | 0.156 | `cart` | 0.800 | 0.693 |
| staphylococcus_aureus__erythromycin | 1,957 | 948 | 0.069 | `scm` | 0.894 | 0.696 |
| mycobacterium_tuberculosis__ethionamide | 1,774 | 915 | 0.144 | `gbt` | 0.733 | 0.715 |
| enterococcus_faecium__dalfopristin_quinupristin | 1,118 | 112 | 0.159 | `scm` | 0.722 | 0.523 |
| klebsiella_oxytoca__cefoxitin | 99 | 47 | 0.074 | `rf` | 0.842 | 0.782 |
| enterococcus_faecium__erythromycin | 1,151 | 1,063 | 0.108 | `scm` | 0.844 | 0.634 |

## Method

The five's calls are the release's saved out-of-fold scores, called resistant at 0.5 or above; no learner is fitted again. The comparators are each of the five; the best of the five in hindsight, the plan's comparator; the site's best, the learner the release names best, chosen on the random split; the mean of four, the scores of `cart`, `rf`, `gbt` and `lasso` averaged per genome; the lineage baseline, which calls by lineage identity inside the release's own outer folds and has balanced accuracy only; and the catalog rule, on the 20 of 24 datasets where the release scores it.

The plan named balanced accuracy, recall on the resistant class, specificity, log loss with probabilities clipped to [1e-4, 1 - 1e-4], and AUROC, each over a dataset's genomes; the experiment's report gives figures over datasets for the first, second and fourth only, and those three are the measures here.

A figure over datasets is the median of the network minus the comparator, with the 2.5% and 97.5% percentiles of that median over 2,000 resamples of the 12 species, each taking every dataset of each drawn species, at seed 20261115.

The plan's rule, on the lineage split:

| Rule | Statistic | Value | Needs | Holds |
|---|---|---|---|---|
| 1 | median network minus comparator balanced accuracy | -0.059 | at least 0.02 | no |
| 1 | lower end of its interval | -0.132 | above 0 | no |
| 2 | median network minus baseline balanced accuracy | +0.115 | at least 0.203, the comparator's | no |
| 3 | median network minus comparator recall | -0.102 | at least -0.01 | no |
| 4 | CPU-hours of network fitting | 1.838 | at most 8.0 | yes |
| 4 | peak memory of a network process, kB | 299,972 | at most 12,582,912 | yes |

All four must hold to adopt the network. Here rules 1, 2 and 3 fail, so the outcome is "no gain established".

## Results

The network minus each comparator: the median over the 24 datasets (20 for the catalog rule) with its 95% interval over species. A positive log loss difference means the network's is higher.

Lineage split:

| Comparator | Balanced accuracy | Recall | Log loss | Network higher |
|---|---|---|---|---|
| `scm` | +0.003 (-0.086 to +0.049) | -0.121 (-0.248 to -0.079) | -2.650 (-3.299 to -1.545) | 12 of 24 |
| `cart` | -0.012 (-0.052 to -0.001) | -0.024 (-0.134 to +0.007) | -0.291 (-0.407 to -0.203) | 7 of 24 |
| `rf` | -0.016 (-0.032 to +0.009) | +0.018 (-0.034 to +0.072) | +0.132 (+0.057 to +0.229) | 8 of 24 |
| `gbt` | -0.012 (-0.035 to +0.069) | -0.006 (-0.070 to +0.101) | +0.056 (+0.001 to +0.193) | 10 of 24 |
| `lasso` | +0.015 (-0.007 to +0.071) | -0.002 (-0.014 to +0.105) | -0.002 (-0.027 to +0.113) | 15 of 24 |
| best of the five | -0.059 (-0.132 to -0.044) | -0.102 (-0.356 to 0.000) | +0.021 (-0.262 to +0.209) | 1 of 24 |
| site's best | -0.036 (-0.097 to +0.008) | -0.035 (-0.265 to +0.017) | +0.026 (-0.179 to +0.226) | 7 of 24 |
| mean of four | -0.007 (-0.019 to +0.014) | +0.027 (-0.042 to +0.128) | +0.127 (+0.042 to +0.241) | 10 of 24 |
| lineage baseline | +0.115 (+0.050 to +0.189) | not scored | not scored | 23 of 24 |
| catalog rule | +0.007 (-0.058 to +0.024) | -0.116 (-0.277 to +0.102) | -1.860 (-2.836 to -0.675) | 11 of 20 |

Random split:

| Comparator | Balanced accuracy | Recall | Log loss | Network higher |
|---|---|---|---|---|
| `scm` | -0.011 (-0.027 to +0.002) | +0.025 (-0.075 to +0.043) | -1.610 (-2.923 to -0.821) | 9 of 24 |
| `cart` | -0.004 (-0.031 to +0.010) | +0.004 (-0.032 to +0.040) | -0.257 (-0.332 to -0.160) | 12 of 24 |
| `rf` | -0.021 (-0.042 to -0.007) | -0.001 (-0.015 to +0.005) | +0.097 (+0.037 to +0.158) | 5 of 24 |
| `gbt` | -0.026 (-0.055 to -0.018) | -0.005 (-0.077 to +0.021) | +0.064 (+0.016 to +0.107) | 4 of 24 |
| `lasso` | -0.017 (-0.057 to -0.004) | -0.002 (-0.096 to +0.016) | +0.050 (-0.039 to +0.101) | 4 of 24 |
| best of the five | -0.043 (-0.064 to -0.026) | -0.015 (-0.093 to +0.003) | +0.037 (-0.134 to +0.100) | 1 of 24 |
| site's best | -0.043 (-0.064 to -0.026) | -0.015 (-0.093 to +0.003) | +0.037 (-0.134 to +0.100) | 1 of 24 |
| mean of four | -0.041 (-0.058 to -0.019) | -0.014 (-0.055 to +0.001) | +0.085 (+0.031 to +0.130) | 4 of 24 |
| lineage baseline | +0.066 (+0.032 to +0.120) | not scored | not scored | 19 of 24 |
| catalog rule | +0.115 (-0.015 to +0.174) | -0.028 (-0.039 to +0.370) | -2.459 (-3.175 to -0.773) | 14 of 20 |

## Review

A review on 15 Nov 2026 recomputed every figure of the experiment's answer with its own code, from the network's outputs and the release, and found each equal to the 4th decimal. It supports "no gain established" for this one network on these 24 datasets, and no wider reading. It narrowed the reading in five ways:

- **Comparator.** On the lineage split each of the five sits below the best of the five by a median of 0.024 (`rf`), 0.047 (`cart`), 0.052 (`scm`), 0.085 (`gbt`) and 0.088 (`lasso`), so none of the five could have met rule 1 either. The plan fixed this comparator before the results and gave its reason; these medians were read after them, so they are exploratory.
- **Site's best.** It has the best random-split balanced accuracy of the five on 24 of 24 datasets and the best lineage-split one on 14 of 24, so it is a comparator fixed in advance on the lineage split only.
- **Reach.** Against the comparators that do not read the lineage results, the network is shown neither higher nor lower on the lineage split, except `cart`, whose interval lies below 0. The figures do not show the network lower than GRM's learners.
- **First look.** The analysis ran on the first datasets' results at 12:25:38Z, before the rest ran, and its script last changed at 13:46:06Z; the experiment's report did not say so. The review's recompute from the plan's text gives the same figures.
- **Holdout and wording.** No test label or genome reaches a screen, so the holdout's part cannot raise test scores. The plan calls its setting a member of the dense network family; the result is about this network on this screen, and not about the family.

On the random split the intervals for `rf`, `gbt`, `lasso`, the best of the five, the site's best and the mean of four lie below 0. On the lineage split the site's best's upper end, +0.008, is near enough 0 that another order of species in the draw puts it below.

## Limits

- One network at one setting, on one screen and 24 datasets. No transformer, sequence model, other architecture or other setting was fitted, and the record says nothing about them.
- The datasets are not a sample of the release: 12 were chosen for the largest gap, and 12 drawn once from those with a gap of 0.05 or more.
- The interval resamples 12 species with 1 to 5 datasets each, so it reflects few species.
- The review did not fit the network again. It read the network's scores as written and checked the bookkeeping of all 720 fits.
- The plan named specificity and AUROC among the measures. The experiment's report gives no figure over datasets for either, and the record holds neither.
- File times show the order of the plan and the results, but they do not prove a file's content at a past time.
- Every dataset comes from one public database, so this is not external validation.

## Files

- `reference/dense-network/` holds the record: measures per dataset, split and learner, and medians and counts over datasets, with no genome, no sequence and no fitted parameter.
- `scripts/dense-network/record.py` writes it from the experiment's summary file and the review's report, checks every derived figure against the experiment's own, and rebuilds the derived files from the record alone.
- `tests/dense_network` holds this page to the record line for line.

```sh
python3 -I scripts/dense-network/record.py --check reference/dense-network
python3 -I -m unittest discover -s tests/dense_network -t .
```

The experiment's and the review's output directories are not public. With them, `--source` and `--review` rebuild every file of the record.
