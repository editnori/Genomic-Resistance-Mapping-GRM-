# Lane clinda: plan, fixed before any count or score is read

Lane `clinda` of the program of 14 Nov 2026. Experiment on real public data. Repository commit `c435000`. Release read: `rev14-274-20261112` (evaluation revision 14).

## 0. What this is and is not

This is a new measurement on the revision 14 cohort. It does not reproduce the lost October numbers and does not try to. The October run is history from another revision and cohort: 1,543 samples of the revision 7 candidate, 161 of 850 known-gene false-resistant calls, 128 of the 161 from a study the record calls "Wang 2021", and a set covering machine figure of 0.772809 against 0.870555 with that study left out. None of those is a target here.

The record (work program rows WP:87 and WP:252, map line 317, row S11 of the lost-results table) does not name the species, gives no PMID for the study, and does not say whether the October comparison was a refit or a rescoring. So the species and the study cannot be identified from the record on this machine. The study is defined by the rule in section 3.

Scope. This is a diagnostic and statistical check of labels and scores. A change in score when a study is left out is an association with that study's labels or genomes. It does not show which label is wrong and it is no statement about any gene. Nothing here ranks or proposes a determinant or an intervention. The known-gene rule is GRM's own reader of the AMRFinderPlus catalog. Every threshold below is an adapted research screen chosen for this lane, not a validation criterion of any regulator.

## 1. What was read before this plan was hashed

- `RULES.md`, the brief, the lost-results report and row S11, the work program digest and the repository's work program lines that mention clindamycin, map line 317.
- The release index: the keys and value types, and the ids of the datasets whose antibiotic is clindamycin. From the six clindamycin reports: key names and value types, `eval_config`, `counts.k`, `counts.min_support`, `counts.min_class`, and the `seconds` fields (for job sizing).
- The data root: file names, the key names and value types of the first record of `curated/phenotypes.jsonl`, `bvbrc/genome_amr.jsonl` and `bvbrc/genomes.jsonl`, the line counts of the six species' `genomes.txt`, the shape of `curated/drug_names.json`, and, over all 1,285,111 raw records of every drug and species together, how often `pmid` and `source` are filled (1,144,137 records have a `pmid`; 140,974 have neither; `source` is filled on 30,116 records, never without a `pmid`).
- Repository source: `src/main.rs`, `src/curate.rs` (label parsing and the PMID split), `src/analysis.rs` (cohort selection), `src/layout.rs`, `scripts/run-pipeline.sh`, `scripts/ablation/stats.py` (the log loss clip).
- One October table row in `docs/lasso-stopping.md:293` that names `staphylococcus_aureus__clindamycin` with three numbers of a lasso timing table. It is not a revision 14 count.

Not read: any revision 14 label, known-gene call, PMID value by dataset, sample count by class, metric or score of a clindamycin dataset.

## 2. Datasets

All six clindamycin datasets of the revision 14 index:

`campylobacter_coli__clindamycin`, `campylobacter_jejuni__clindamycin`, `enterococcus_faecium__clindamycin`, `staphylococcus_aureus__clindamycin`, `streptococcus_agalactiae__clindamycin`, `streptococcus_pneumoniae__clindamycin`.

A dataset's samples are the `samples` array of its released report: `id` (BV-BRC genome id), `resistant` (the reported phenotype the dataset uses) and `catalog_resistant` (the known-gene rule's call). Check against the data root: every sample must have one clindamycin row in `curated/phenotypes.jsonl` whose `reported` call equals the report's label, with `conflict` false and `label_from_mic` false. Any mismatch is reported and stops that dataset.

## 3. Study field and the rule that picks the study

**Field.** A genome's study is read from the `pmid` field of its raw BV-BRC records in `/home/user/grm-h1raw/bvbrc/genome_amr.jsonl`. The curated tables keep no study field (`curated/phenotypes.jsonl` has none, and the root holds no `curated/provenance.jsonl`), so the raw records are joined to the curated row by `genome_id`. The records used for a sample are those with its `genome_id`, an `antibiotic` name that `curated/drug_names.json` resolves to the canonical drug `clindamycin`, and a `resistant_phenotype` that parses (by the rule of `parse_call` in `src/curate.rs`) to the call that became the sample's label. These are the records that carry the reported phenotype. The `pmid` text is split on `;` and trimmed, as `curate.rs` does.

**Study key.** A sample's key is the sorted set of its distinct PMIDs joined with `+`, or `no-pmid` when its label records carry none. The table by study uses the key, so its rows add up to the dataset. The `source` field is shown beside it for information only.

**Membership.** A sample belongs to PMID p when p is in its key. A sample with two PMIDs belongs to both.

**Choice, on counts alone, before any fit.** For dataset d and PMID p, let D(d,p) be the number of samples of d that belong to p and on which the known-gene rule's call differs from the reported phenotype (false-resistant: rule resistant, reported susceptible; false-susceptible: rule susceptible, reported resistant). The chosen study of d is T(d) = the PMID with the largest D(d,p). Ties go to the larger false-resistant count, then the larger sample count, then the smaller PMID number. `no-pmid` is not a study and cannot be chosen. A dataset in which no PMID has a disagreement has no chosen study and is not refit.

**Primary pair.** The lane's headline pair is the dataset and study with the largest D(d,T(d)) over the six datasets (ties: larger share, then dataset id in alphabetical order). The other datasets are reported the same way as secondary results.

The October figure counted false-resistant calls only. The PMID with the most false-resistant calls is listed beside T(d); when it differs, that is said, and it is not refit.

## 4. Step one: counts, with no fit

For each dataset, one row per study key with: samples, reported resistant and susceptible, rule resistant and susceptible, false-resistant, false-susceptible, disagreements, disagreement rate. Totals per dataset.

For the chosen study T(d):

- share of disagreements = D(d,T) / D(d), with numerator, denominator and an exact (Clopper-Pearson) 95% interval; the same for false-resistant and for false-susceptible calls;
- share of samples = n(d,T) / n(d), for comparison;
- disagreement rate inside the study against the rate in the rest, with the odds ratio and Fisher's exact two-sided p. Genomes of one study are not independent, so the interval and the p are descriptive.

Words fixed now (adapted research screens): a share of 50% or more is "most", 25% to under 50% is "a large share", under 25% is "not a large share". The words are used only with the counts beside them and with the study's share of samples.

## 5. Step two: refit

**Build.** `grm` built from a scratch copy of the repository at `c435000` with `RUSTUP_TOOLCHAIN=1.98.1`, `cargo build --release --locked`, target directory `/home/user/scratch/clinda/target`. The binary's sha256 is recorded.

**Stages, with the repository's own commands.** For each of the six species, in a scratch data root that links the public root's files: `grm count --species S`, `grm structure --species S`, `grm matrix --species S` (defaults: k 31, minimum support 2), then `grm learn --species S --only S__clindamycin`. These are the stage functions `grm pipeline` calls; the pipeline itself is not used because it would learn every drug of the species. `RAYON_NUM_THREADS=8`. Evaluation settings are the binary's defaults, which must equal the released report's `eval_config`.

**Arm "with".** The released cohort, unchanged. Acceptance: the refit report's sample ids, labels, known-gene calls, lineage ids, `input_fingerprint` and every learner's out-of-fold scores on both splits equal the released report's. If any differs, the difference is reported, the released figures are shown beside the refit, and effects use this lane's own "with" refit.

**Arm "without".** The same data root with one file changed: `curated/phenotypes.jsonl` loses the clindamycin rows of the genomes that belong to T(d), and nothing else. The genomes stay in the species matrix and in the lineage clustering, so the k-mer patterns and lineage ids of every other sample are the same files in both arms. Then the same `grm learn` command. Proof that nothing else changed: the same binary hash, the same linked input files, a line diff of the two phenotype files (only the removed rows), equal `eval_config`, and equal label, known-gene call and lineage id for every remaining sample. What changes by construction is recorded: the sample set, the fold plans drawn over it, the per-fold feature screen and tuning, and the dataset's minimum support, which the code derives from the sample count.

**Refit set.** Every dataset with a chosen study in which both classes keep at least `min_class` (25) samples after removal. Where removal leaves a class under 25, GRM writes no report; that is recorded as "cannot be scored without the study" and no refit is made.

**Placebo removals (secondary).** To show how much scores move when the same number of genomes is removed at random, each refit dataset also gets 10 placebo arms: seed s in 1..10 removes, from the samples outside T(d), as many reported-resistant and as many reported-susceptible samples as T(d) holds (drawn with Python's `random.Random(20261114*100 + s)` over the sorted ids). If the rest holds too few of a class to match, or the removal would leave a class under 25, the dataset has no placebo and that is said. Placebo fits are launched in descending order of D(d,T), and none is launched after 21:30Z on this machine's clock; the number completed is reported.

## 6. Measures

For each dataset, each split (`random`, `lineage`), each of the five learners (`scm`, `cart`, `rf`, `gbt`, `lasso`) and the known-gene rule (`catalog` row), from the out-of-fold scores in the report, with the release's cutoff (a score of 0.5 or more is a resistant call; this is checked by recomputing the released counts from the released scores):

- recall (sensitivity, resistant is the positive class) and specificity, each with an exact Clopper-Pearson 95% interval, the release's own interval method;
- balanced accuracy, the mean of the two;
- log loss: the mean over samples of the negative log likelihood, with scores clipped to [1e-4, 1 - 1e-4] as in `scripts/ablation/stats.py`. The set covering machine and the known-gene rule give 0 or 1, so their log loss is a function of their error count and is marked as such.

No learner is singled out: all five and the rule are shown side by side, recall beside log loss and balanced accuracy.

Three views of each row:

- V1: arm "with", all samples (the released figure);
- V2: arm "with", scored only on the samples outside T(d) (no refit);
- V3: arm "without" (refit), the same samples as V2.

Effects:

- total = V3 - V1: what leaving the study out does to the figure;
- composition = V2 - V1: the part that comes from no longer scoring the study's samples;
- training = V3 - V2: the part that comes from no longer fitting on them, on identical samples.

Also shown: the "with" arm scored on T(d)'s samples alone, where both classes are present.

**Intervals for balanced accuracy, log loss and every effect.** Percentile bootstrap, 2,000 replicates, seed 20261114, resampling whole lineage clusters (`samples[].lineage`) of the full cohort with replacement; each replicate recomputes V1, V2 and V3 from the resampled samples. A replicate where a measure is undefined (a class absent) is dropped and counted. A dataset with fewer than 20 lineage clusters is flagged `few_clusters`. A second interval from resampling single samples is given beside it as a sensitivity reading. The intervals describe sampling of genomes. They do not cover the change in fold plan between arms (one seed, 42, as in the release); the placebo arms are the reading of that.

## 7. What result means the released clindamycin figures need a caveat

Decided per dataset, on the lineage-cluster bootstrap. The measures are the ones the release publishes: balanced accuracy, recall and specificity. Log loss is shown beside them and does not enter the rule.

A released clindamycin dataset needs a caveat when either holds:

1. **Existence.** Removing T(d) leaves a class under 25 samples, so the dataset cannot be scored without that study.
2. **Dependence.** For the known-gene rule or for any of the five learners, on either split, the total effect (V3 - V1) on balanced accuracy, recall or specificity is at least 0.05 in absolute value and its 95% interval excludes 0.

The caveat then names the rows that pass and the study, and says the figure is associated with that study's labels or genomes. The report gives the number of rows that pass out of the 36 tested per dataset (6 rows, 2 splits, 3 measures), so that a single passing row among 36 is read with that in mind, and it gives the placebo spread beside each passing row.

If no dataset meets either rule, the answer is "no caveat under this plan's rule", and the limits say what was not tested: datasets with no chosen study, and effects under 0.05.

The concentration counts of section 4 are reported whatever the scores do. A large share with no score change is reported as that.

## 8. After results

Anything decided after a count or score is read is labelled exploratory in the report. Foreseen exploratory step: a lookup of the chosen PMIDs in a public literature index (the PMID alone is sent) to say whether the record's "Wang 2021" could be the same study.
