# Kover comparison

This page compares Kover 2.0.7 with GRM's `scm` learner on 12 public BV-BRC datasets of release `rev14-274-20261112` (evaluation revision 14), the smallest eligible dataset of each of 12 species. Both are scored on the release's lineage folds, on the same 1,118 genomes. It is not `reference/kover-summary.json`, which holds the 73 historical Kover rows of the legacy application, on other cohorts and random splits, and which `grm report` reads. This is a comparison of two programs on 12 public datasets. It shows no equivalence and is not external validation.

## Reading

On these 12 datasets, Kover's held-out calls have lower balanced accuracy than the released lineage calls of GRM's `scm` learner, and the size of the gap turns on how Kover breaks exact ties: Kover minus GRM is -0.052 (-0.094 to -0.018, a 95% interval over the 12 datasets) as run, and -0.035 (-0.060 to -0.013) in an exploratory reading that counts each tied first rule equally in Kover's 41 one-rule folds of 60.

The gap is in specificity, -0.104 (-0.215 to -0.022), and not in recall, +0.001 (-0.047 to +0.046); log loss, +0.374 (+0.071 to +0.748), counts the same errors, since every call is 0 or 1.

Kover's balanced accuracy is lower in 8 of the 12 datasets, equal in 2 and higher in 2, and in no single dataset does its t interval exclude 0.

Which differences between the two runs make the gap is not separated, and the result is not shown to carry to larger datasets.

## Plan

The plan was fixed on 15 Nov 2026: its hash list, recorded at 12:05:29Z, gives `PLAN.md` the sha256 `e9bd5114556f1a8507d5ab1adecae5f5b6861796378597770588b9e4886b26f8`. The first Kover fit on a project genome started at 12:13:08Z, after it.

The plan chose the 12 datasets by a rule on sizes alone: of the 268 datasets of the release's 274 whose `scm` row has a lineage entry, the one with the fewest genomes in each species, ties by dataset id, kept at 150 genomes or fewer. `datasets.tsv` in the record is the file its hash list names, byte for byte.

It also fixed Kover's settings, the measures and the reading: a measured difference when the 95% interval over datasets lies wholly above or below 0, no difference shown when it contains 0, and no equivalence claimed either way.

An independent review hashed its own plan at 16:02:53Z and recomputed every figure. It chose the tie readings after reading the results per fold, so they are exploratory; its tie summary was written at 16:45:54Z.

## Datasets

Each difference is Kover minus GRM, the mean over the dataset's 5 lineage folds; the interval for balanced accuracy is a 95% t interval over those folds.

| Dataset | Genomes | Resistant | Balanced accuracy | Recall | Log loss | Specificity |
|---|---|---|---|---|---|---|
| klebsiella_pneumoniae__cefpodoxime_ceftriaxone | 56 | 28 | 0.000 (0.000 to 0.000) | 0.000 | 0.000 | 0.000 |
| mycobacterium_tuberculosis__amoxicillin | 57 | 25 | -0.038 (-0.268 to +0.192) | 0.000 | +0.377 | -0.076 |
| klebsiella_aerogenes__imipenem | 68 | 26 | -0.050 (-0.199 to +0.099) | 0.000 | +0.532 | -0.100 |
| klebsiella_oxytoca__cefpodoxime | 74 | 47 | -0.237 (-0.501 to +0.028) | +0.127 | +2.008 | -0.600 |
| vibrio_cholerae__sulfamethoxazole | 74 | 37 | 0.000 (0.000 to 0.000) | 0.000 | 0.000 | 0.000 |
| enterobacter_cloacae__sulfamethoxazole_trimethoprim | 80 | 41 | -0.082 (-0.257 to +0.093) | -0.130 | +1.105 | -0.033 |
| acinetobacter_baumannii__moxifloxacin | 91 | 60 | +0.011 (-0.020 to +0.042) | +0.022 | -0.123 | 0.000 |
| klebsiella_quasipneumoniae__cefotaxime | 101 | 62 | -0.037 (-0.126 to +0.051) | -0.025 | +0.395 | -0.050 |
| staphylococcus_haemolyticus__tetracycline | 104 | 37 | -0.076 (-0.246 to +0.093) | -0.168 | +0.452 | +0.015 |
| clostridioides_difficile__ceftriaxone | 136 | 90 | +0.003 (-0.045 to +0.051) | +0.031 | -0.040 | -0.024 |
| escherichia_coli__cefpodoxime_clavulanate | 136 | 101 | -0.113 (-0.297 to +0.071) | +0.117 | -0.003 | -0.343 |
| streptococcus_agalactiae__clindamycin | 141 | 115 | -0.0002 (-0.038 to +0.037) | +0.040 | -0.218 | -0.040 |

In 2 datasets the two make the same call on every genome.

## Method

Kover 2.0.7 (core-2.0.7, cli-2.0.1) ran with its k-mer tools built against GATB-core 1.2.2, on the decompressed public genomes, with these commands for each dataset and each lineage fold k:

```sh
kover dataset create from-contigs --kmer-size 31 --n-cpu 4
kover dataset split --id lineage_fold<k> --train-ids <train> --test-ids <test> --folds 3 --random-seed 42
kover learn scm --model-type conjunction disjunction --p 0.1 0.316 0.562 1.0 1.778 3.162 10.0 --max-rules 6 --hp-choice cv --n-cpu 1 --random-seed 42
```

K-mers seen in one genome only are left out, Kover's default. A genome in `test_correct` gets its label and one in `test_errors` the other class. GRM's side is read from the release and not fitted again: the lineage `scm` row's calls, each from the fold that did not train on that genome.

Log loss takes each call as a probability clipped at 0.0001, so it is 0.0001 plus 9.2102 times the error rate. Each dataset's figure is the mean of its fold differences, and the figure over datasets is their mean, with a 95% interval from 10,000 resamples of whole datasets; `summary.tsv` gives each measure's seed.

The plan's table of settings, with whether Kover's equals GRM's:

| Setting | GRM `scm` | Kover as run | Same |
|---|---|---|---|
| k-mer | canonical 31-mers | canonical 31-mers (DSK) | yes |
| Features offered | training patterns with support 3 to n_train - 3, one per distinct training pattern | every k-mer in at least 2 genomes of the whole dataset; presence and absence rules | no |
| Model-size limit | 6 rules | --max-rules 6 | yes |
| Kinds | conjunction, disjunction | both | yes |
| p values | 0.1, 0.316, 0.562, 1.0, 1.778, 3.162, 10.0 | the same 7 (Kover's default list adds 999999.0; not used) | yes |
| Utility | covered - p * lost | the same | yes |
| Rule tie | fewer lost, lower index, presence first | equal within np.isclose, then lowest training error of the rule alone, then the first in Kover's k-mer order | no |
| Inner folds | 3, whole lineages placed to balance the classes | 3, drawn at random by Kover, not grouped, not stratified | no |
| Selection criterion | mean balanced error | mean error rate | no |
| Selection tie | keep earlier (conjunction, smaller p, fewer rules) | fewer rules, then p nearer 1, then arrival order (--n-cpu 1 fixes it) | no |
| Fewest rules | 1 | 0 (an empty model) is allowed | no |
| Seeds | 42 plus 101 times the fold | 42 | number only; the draws differ |

The rows marked no are known sources of difference, and the plan does not separate them.

## Results

Over the 12 datasets: the mean of Kover minus GRM with its 95% bootstrap interval, and the datasets and folds where Kover is higher, equal and lower. Recall stands beside log loss.

| Measure | Kover minus GRM | Interval | Datasets | Folds | Reading |
|---|---|---|---|---|---|
| Balanced accuracy | -0.052 | -0.094 to -0.018 | 2, 2, 8 | 9, 29, 22 | planned |
| Recall | +0.001 | -0.047 to +0.046 | 5, 4, 3 | 13, 39, 8 | planned |
| Log loss | +0.374 | +0.071 to +0.748 | 6, 2, 4 | 17, 32, 11 | planned |
| Specificity | -0.104 | -0.215 to -0.022 | 1, 3, 8 | 5, 38, 17 | planned |
| Error rate | +0.041 | +0.008 to +0.081 | 6, 2, 4 | 17, 32, 11 | exploratory |
| Rules per fold | -0.317 | -0.683 to -0.033 | 2, 4, 6 | 9, 35, 16 | planned |

The two make the same call on a mean share of 0.874 of each dataset's genomes (0.799 to 0.939).

Exploratory: in 41 of Kover's 41 one-rule folds, more than one rule ties at the first step, Kover's second tie rule keeps them all, and the order of its k-mer columns decides. In 21 of them the tied rules give different test results. Only Kover's side of those folds changes below.

| Kover's side | Kover minus GRM | Interval |
|---|---|---|
| as run, the first tied rule | -0.052 | -0.094 to -0.018 |
| each tied rule equally likely | -0.035 | -0.060 to -0.013 |
| the lowest tied result, a bound chosen by test result | -0.061 | -0.104 to -0.029 |
| the highest tied result, a bound chosen by test result | -0.007 | -0.033 to +0.015 |

In 41 of the 41 folds the first tied rule gives Kover's own calls. GRM breaks its ties by feature index, an order just as arbitrary, and its fold rules are not in the release, so its side cannot be averaged the same way.

## Review

The review reproduced every figure, given the seed the code gives each measure. It narrowed the first report of the comparison in these places:

- Size. The first report gave -0.052 as the gap. Its size turns on Kover's tie order, as the table above gives, and the interval over datasets does not cover that variation.
- Cause. The first report said the gap comes from the parts of the runs that differ by design. With settings held equal, covered and lost counts agree at 19 of 19 greedy steps of the 12 fits on all genomes, and the rules' patterns at 13 of 19. The fold fits were not checked this way, so those parts are candidate causes and are not separated.
- First step. Where the patterns differ, in 3 datasets, the first step is an exact tie on both sides, and Kover's rule has the pattern of a rule GRM lists as tied in 3 of 3.
- Reach. The first report said its results do not carry to large datasets. The data show only that they are not shown to carry.
- Log loss. With calls of 0 or 1, log loss is 0.0001 plus 9.2102 times the error rate, so +0.374 is the error rate higher by 0.041 (0.008 to 0.081), the same result in other units.
- Seeds. Each measure's bootstrap uses the seed the code gives it, which `summary.tsv` lists, and no interval changes.

## Limits

- The 12 datasets are the smallest eligible one of each species: 56 to 141 genomes, a median of 85.5 against 693.5 over all 268 eligible datasets. 11 of the 23 species have none among them.
- They differ from the rest: a median resistant share of 0.563 against 0.477, 7 of 12 resistant-majority against 123 of 268, and a median of 38.5 lineages against 93.
- Test folds hold 8 to 46 genomes, so one genome moves a fold's figure in a large step, and the interval over datasets rests on 12 points.
- The interval covers the choice of datasets only. It holds no variation from Kover's tie order or from its seed.
- Three settings favor GRM under balanced accuracy on lineage folds: Kover selects by error rate, its inner folds are random and not grouped by lineage, and it allows an empty rule set, which it chose in 2 of 60 folds.
- Kover ran with one seed, 42. The spread over seeds is not measured.
- Two of the 12 datasets have a row in `reference/kover-summary.json`, on 766 and 100 genomes against 57 and 104 here. Those rows come from other cohorts and random splits and are not comparable.
- No rule is named here. A rule term is a marker and an association.

## Files

- `reference/kover-comparison/` holds the record: counts and measures per dataset, per fold, per greedy step and per measure, with no genome, no sequence and no rule.
- `scripts/kover-comparison/record.py` writes it from the comparison's reduced files and the review's tie files, and rebuilds the derived files from the record alone.
- `tests/kover_comparison` holds this page to the record, line for line.

```sh
python3 -I scripts/kover-comparison/record.py --check reference/kover-comparison
python3 -I -m unittest discover -s tests/kover_comparison -t .
```

The comparison's and the review's output directories are not public. With them, `--source` and `--review` rebuild every file of the record.
