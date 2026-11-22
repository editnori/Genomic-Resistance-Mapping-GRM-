<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/wordmark-dark.svg">
  <img src="assets/wordmark-light.svg" width="322" height="80" alt="GRM, Genomic Resistance Mapping. The mark is a cube seen from one corner. Each side face is a stack of three tiles in two blues that swap from row to row, and the top face is coral.">
</picture>

GRM measures how well a bacterial isolate's genome predicts its laboratory antimicrobial susceptibility result. It uses public assemblies from BV-BRC (the Bacterial and Viral Bioinformatics Resource Center) that have a reported resistant (R) or susceptible (S) call.

A Rust pipeline builds one dataset per species and antibiotic, cross-validates five learners on it, and writes a report. A React dashboard shows each report's model performance, gene and panel evidence, and MIC (minimum inhibitory concentration) predictions. The site is [genomicresistance.com](https://genomicresistance.com).

## What GRM found (Revision 14)

On public BV-BRC genomes, the best of five k-mer learners reaches a median balanced accuracy of 92.9% on random splits (274 datasets) and 87.7% with whole lineages held out (268 datasets). With lineages held out, 0 of 268 datasets are shown to meet both limits of the adapted research screen. Folds that hold whole studies out score worse than random folds on log loss and specificity for 5 of 5 learners. Every figure is retrospective. No external or clinical validation has been done.

| Question | Revision 14 finding | Details |
| --- | --- | --- |
| Coverage | 274 datasets across 23 species have a report (268 with a lineage split) | [release-summary.md](docs/release-summary.md) |
| Random vs. lineage split (best learner) | Median balanced accuracy 0.9293 random (274) vs. 0.8768 lineage (268); paired median drop 0.0377, 116 of 268 lose >5 points | [release-summary.md](docs/release-summary.md) |
| Five learners on unseen lineages | SCM 0.8708, CART 0.8684, lasso 0.8565, GBT 0.8409, RF 0.8124 (268 datasets each); SCM's paired gap to CART (0.0026) and to lasso (0.0046) does not clear zero | [release-summary.md](docs/release-summary.md) |
| Adapted research error screen (lineage, best learner) | 0 of 268 shown to meet VME ≤1.5% and ME ≤3% by 95% intervals (239 shown to miss, 29 not decided); 14 of 268 meet FDA 2009 guidance wording | [release-summary.md](docs/release-summary.md) |
| Known-gene rule vs. best learner (lineage) | Rule ahead or within 2 points on 108 of 249 (74 of 80 where rule ≥0.85); learner ahead by >2 points on 141 of 249 | [release-summary.md](docs/release-summary.md) |
| Study split inside BV-BRC | Holding BioProjects whole scores worse than random folds for 5 of 5 learners on log loss and specificity across 70 datasets in 9 species (recall worse for RF and GBT) | [study-split.md](docs/study-split.md) |
| Newly deposited BV-BRC studies | On 7,258 unseen-study genomes (rev 12 reports), lineage figures hold on 32 of 81 datasets for SCM, 23 of 81 for CART, and 7 of 25 for lasso | [outside-study-results.md](docs/outside-study-results.md) |
| Another public database (NCBI) | On 13,377 genomes (9 species), SCM meets the plan's rule in 6 of 8 species (52 of 84 datasets; plus 14 of 25 in *E. coli* pilot), CART in 5 of 8 (50 of 84; plus 12 of 25 in *E. coli* pilot), lasso in 3 of the 4 species where it can be read (17 of its 24 datasets across the eight) | [other-database-test.md](docs/other-database-test.md) |
| Input settings (4 species, H-1) | No single change shown better than default across 49 answered datasets (cap-4000 better 31 of 49; cap-1000 worse 25 of 33; budget-1000 worse 37 of 49) | [input-comparison-runner.md](docs/input-comparison-runner.md) |
| First-public-day split | Training on earlier public dates scores worse than a size-matched random draw on log loss and balanced accuracy for 5 of 5 learners across 24 datasets in 4 species (recall worse for 2 of 5) | [time-split.md](docs/time-split.md) |
| Kover 2.0.7 vs. GRM SCM | On the 12 smallest lineage-split datasets (1,118 genome-dataset entries), Kover minus GRM SCM balanced accuracy is -0.0516 (lower on 8 of 12, higher on 2, equal on 2) | [kover-comparison.md](docs/kover-comparison.md) |
| Dense neural network | Across 24 datasets in 12 species (one setting on 4,000 screened features), no gain established: median -0.0593 balanced accuracy vs. hindsight best of 5 learners (higher on 1 of 24) | [dense-network.md](docs/dense-network.md) |
| AMRFinderPlus 4.2.7 program | Across 235 datasets in 21 species, the full program leads GRM's known-gene rule on 94 of 235 (48 by >=10 points, 20 quinolones); SCM minus program median balanced accuracy is +0.0048 (95% species CI -0.0063 to +0.0466) | [amrfinderplus-comparison.md](docs/amrfinderplus-comparison.md) |
| TBProfiler 6.6.3 (*M. tuberculosis*) | On 2,000 assemblies across 14 graded TB drugs (16 named endpoints, 819 no-calls), TBProfiler scores higher than SCM on 8 of 14 (SCM minus TBProfiler median -0.0058) | [tbprofiler-comparison.md](docs/tbprofiler-comparison.md) |

The VME ≤1.5% and ME ≤3% limits are an adapted research screen. They are not an applicable FDA criterion, and meeting them does not make a learner fit for clinical use. Every number in the table above is verified against committed `reference/` and `web/src/data/` files by `tests/readme_summary`. See [docs/INDEX.md](docs/INDEX.md) for the full documentation map.

## Flow

| Step | Reads | Writes |
|---|---|---|
| Collect | BV-BRC | lab records, genome metadata, assemblies, and the AMRFinderPlus reference gene catalog |
| Curate | the lab records and genome metadata | one dataset per species and antibiotic where each class has enough genomes with a reported call |
| Describe | each genome's assembly | its 31-base DNA substrings (k-mers), its lineage, and one feature matrix per species |
| Evaluate | the matrix and the reported calls | scores for five learners and GRM's known-gene rule, on a random split and on a lineage split |
| Report | the saved scores | one report per dataset, and the index the dashboard reads |

[Stages](#stages) gives each step's command and output path.

## Use

- **Read the results.** Open [genomicresistance.com](https://genomicresistance.com), download the reports from the [`data-rev14-274-20261112` release](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev14-274-20261112), or browse [docs/INDEX.md](docs/INDEX.md).
- **Run it.** [Setup](#setup) builds the pipeline and the dashboard. [Run](#run) scores one species and serves the result on your machine.
- **Check a number.** [Checks](#checks) recomputes each learner's counts and balanced accuracy from the saved scores, apart from the code that wrote them.
- **Reuse the name.** `assets/` holds the name beside the mark, for light and dark pages. The mark alone is `web/public/favicon.svg`.

## Status

The site serves evaluation revision 14. Of its 277 curated datasets, 274 have a report and 3 have none, because one class is under the minimum of 25 genomes. The 274 reports cover 23 species. The three without a report are *Escherichia coli* amoxicillin, *Klebsiella pneumoniae* cefotaxime, and *Staphylococcus aureus* methicillin. The published report lists them in `SKIPS.json`.

Every figure is internal validation: cross-validation on public BV-BRC genomes, scored against their reported calls. No external, prospective, or clinical validation has been done. [docs/outside-study-results.md](docs/outside-study-results.md) and [docs/study-split.md](docs/study-split.md) add retrospective tests across studies inside BV-BRC, and [docs/other-database-test.md](docs/other-database-test.md) adds a test on data from another public database. None of these is external validation.

Read each learner's figures with their split and their evaluation revision. The lineage split holds out whole lineages, and a report records why when that split could not be formed. From revision 14 every learner's figures cover every isolate of the split ([docs/equal-scoring.md](docs/equal-scoring.md)).

In each of the six clindamycin datasets one deposited study holds most or a large share of the disagreements between GRM's known-gene rule and the reported phenotype, and the datasets are of four kinds: measured dependence for `staphylococcus_aureus__clindamycin` and `campylobacter_coli__clindamycin`, one study for `enterococcus_faecium__clindamycin`, class balance with dependence not measured for `campylobacter_jejuni__clindamycin` and `streptococcus_agalactiae__clindamycin`, and no dependence under the plan's rule for `streptococcus_pneumoniae__clindamycin`. [docs/clindamycin-labels.md](docs/clindamycin-labels.md) says what each released clindamycin figure rests on.

For the full revision 5 to 14 change history, retired Vast host notes, and revision-5 `research.json` packaging details, see [docs/revision-history.md](docs/revision-history.md). Revision 14 was run for the same 274 datasets on 12 November 2026, from revision 12's inputs. Of the 274 reports, 162 differ from revision 12 only in the revision number and the run time, 39 in panel rows only, 48 only in the rows of a learner with a null fit, and 25 in both. The best learner changed in 21 datasets.

## Givens

[Givens](https://github.com/editnori/givens) is a provisional shared browser of the definitions and assumptions that four projects use: GRM, SigDrug, KSD and CondEvent. Its repository is private for now. GRM's entries there cover evaluation terms only: dataset, split, evaluation revision, report figure, and how a figure reads on newly deposited genomes. They hold no sequence or gene content.

GRM writes its entries to `reference/givens-export/grm.json` with `scripts/givens-export/export.py` ([docs/givens-export.md](docs/givens-export.md)). Each count there is read from a committed summary file or quoted from a document of this repository, and an entry GRM has not counted says so. The same script writes [docs/givens-terms.md](docs/givens-terms.md), a page of those terms and counts to read. Givens keeps a copy of this file under the GRM commit it took it from, beside the file's sha256. GRM reads nothing from Givens when it runs.

## Setup

Use Linux or WSL with Rust 1.89 or newer, Node.js 22 and npm, and Bash, `flock`, and `nice` for the pipeline launcher. Python 3 is needed only for the report checks.

`.cargo/config.toml` sets `target-cpu=x86-64-v3`, so the learners' bit counting uses the CPU's popcount instruction and the binary runs on any x86-64 CPU with AVX2 (Intel Haswell, AMD Zen and newer). For a long run on the machine that builds it, `RUSTFLAGS="-C target-cpu=native" cargo build --release --locked` builds for that CPU instead; that binary may not start on another CPU. The release profile in `Cargo.toml` splits the crate into 4 codegen units. The lasso keeps its speed under any split, and [docs/lasso-layout.md](docs/lasso-layout.md) gives the reason and the command that checks it.

```bash
cargo build --release --locked
npm --prefix web ci
npm --prefix web run build
```

Download the BV-BRC antimicrobial resistance (AMR) records, the genome metadata, and the AMRFinderPlus reference gene catalog, then curate the datasets:

```bash
./target/release/grm fetch amr
./target/release/grm fetch genomes
./target/release/grm fetch catalog
./target/release/grm curate --cap 2000 --min-class 25
```

`grm count` and the pipeline stop when the catalog does not load, because counting records each genome's catalog gene calls. `--allow-missing-catalog` counts anyway, which leaves reports without the catalog baseline and gene panel detections.

Antibiotic names resolve through GRM's bundled vocabulary. To add a SigDrug checkout's reviewed vocabulary, pass `--sigdrug-vocabulary /path/to/sigdrug/data/drug/compendium/facts/vocabulary` to `curate` or set `GRM_SIGDRUG_VOCABULARY`; the report records which vocabulary was used.

Git ignores `data/` (downloads, results, and reports), `target/`, `web/dist/`, and `node_modules/`.

## Run

Run one species, then serve the dashboard:

```bash
bash scripts/run-pipeline.sh --species klebsiella_pneumoniae
./target/release/grm serve --port 8641
```

Open <http://127.0.0.1:8641>. `grm serve` listens only on 127.0.0.1 and serves the app in `web/dist` with the report in `data/report`. For frontend work, `npm --prefix web run dev` starts the Vite dev server, which maps `/data` to the same `data/report`; `GRM_REPORT_DIR=<root>/report` points it at a report built with `--data <root>`.

`scripts/run-pipeline.sh` runs `grm pipeline` one species at a time: it downloads the selected assemblies, counts, builds structure and matrix, learns, and rebuilds the dashboard report, so finished species appear while the rest run. It runs at `nice -n 10` with 64 Rayon threads; set `RAYON_NUM_THREADS` to fit the host. It raises the open-file limit to 8192 because `grm matrix` opens every genome's kset (its k-mer file) at once and can stop with "Too many open files" at a lower limit. It passes `--drop-ksets`, which deletes a species' ksets after its learning succeeds and keeps the matrix and results. It also locks `data/logs/pipeline.lock`, so a second launcher exits.

For a detached run:

```bash
mkdir -p data/logs
nohup setsid bash scripts/run-pipeline.sh \
  --species klebsiella_pneumoniae,staphylococcus_aureus \
  >data/logs/pipeline.log 2>&1 </dev/null &
tail -f data/logs/pipeline.log
```

The launcher writes its PID to `data/logs/pipeline.pid` and, on exit, its exit code to `data/logs/pipeline.exit`. A launcher killed by SIGKILL writes no exit file and leaves its PID file behind, so check that PID and the log as well.

## Deploy

A Cloudflare Worker serves genomicresistance.com from one R2 bucket, with no origin process. [docs/static-hosting.md](docs/static-hosting.md) covers the bucket, the release scripts and their checks. On 13 November 2026 the site answered HTTP 200 with the revision-14 report, 274 datasets over 23 species.

The earlier host-based `grm serve` + Cloudflare Tunnel toolkit stays in `scripts/deployment/` ([docs/deployment.md](docs/deployment.md), [docs/revision-history.md](docs/revision-history.md)).

## Stages

Every command takes `--data <directory>` (default `data`), and the output paths are relative to it:

| Stage | Command | Output |
|---|---|---|
| Download | `grm fetch amr`, `genomes`, `catalog`, `fasta` | `bvbrc/`, `catalog/`, `fasta/` |
| Curate | `grm curate` | `curated/`, `species/<slug>/genomes.txt` |
| Count | `grm count --species <slug>` | `species/<slug>/counted.json`, `species/<slug>/ksets/` |
| Structure | `grm structure --species <slug>` | `species/<slug>/structure.json` |
| Features | `grm matrix --species <slug>` | `species/<slug>/matrix/` |
| Evaluate | `grm learn --species <slug>` | `species/<slug>/results/<id>.json` |
| Evaluate on supplied folds | `grm learn --species <slug> --fold-file <file> --fold-out <dir>` | `<dir>/<dataset>.<scheme>.json` ([docs/fold-file.md](docs/fold-file.md)) |
| Report | `grm report` | `report/` |
| Research | `grm research --revision5-reports <dir> --out <file>` | the `--out` file ([docs/revision-history.md](docs/revision-history.md)) |
| Intervals | `grm intervals` | `report/intervals/<id>.json` ([docs/operations-reference.md](docs/operations-reference.md)) |
| Baselines | `grm baselines` | `report/baselines/<id>.json`, `report/baselines.json` ([docs/baselines.md](docs/baselines.md)) |
| Provenance | `grm fetch genome-dates`, then `grm provenance` | `bvbrc/genome_dates.jsonl` and its manifest, `report/provenance/<id>.json` ([docs/provenance-fields.md](docs/provenance-fields.md)) |
| Predict | `grm predict --report <file> --fasta <files>` | JSON lines on standard output, or `--out` ([docs/prediction-contract.md](docs/prediction-contract.md)) |
| Export | `grm export --species <slug> --out <dir>` | `<dir>/<id>/matrix.bin` and `dataset.json`, for a learner outside GRM ([docs/dataset-export.md](docs/dataset-export.md)) |

`--species` takes comma-separated slugs; `grm pipeline` and `grm fetch fasta` use every curated species without it. `grm learn --species <slug> --only meropenem` learns one antibiotic (`--only` also takes dataset ids), and `--force` recomputes reports even when their inputs have not changed.

To run one species by hand, raise the open-file limit first, as the launcher does:

```bash
ulimit -n 8192
./target/release/grm fetch fasta --species klebsiella_pneumoniae
./target/release/grm count --species klebsiella_pneumoniae
./target/release/grm structure --species klebsiella_pneumoniae
./target/release/grm matrix --species klebsiella_pneumoniae
./target/release/grm learn --species klebsiella_pneumoniae
./target/release/grm report
```

`grm fetch fasta` downloads only missing assemblies. `grm report` publishes saved results without recomputing them. Each build removes the report files it did not publish, so builds take turns through a lock on `data/report`, and a `grm report` started during another build waits for it to finish.

## Evaluation

A dataset is one species and one antibiotic. `grm curate` admits it when the resistant and the susceptible class each reach `--min-class` genomes, counting only reported calls and leaving out genomes reported as both R and S for that drug. Labels inferred from an MIC, using the MIC cutoff that best reproduces the reported calls, stay in the curated phenotypes, marked `label_from_mic`, and no learner trains or is scored on them. Curation records each dataset's class minimum and selects at most `--cap` genomes per species, favoring genomes with more reported calls and with calls in each drug's rarer class. Learning applies the class minimum again among the downloaded genomes, and a dataset that falls short there logs `<id>: too few labeled genomes after download, skipped` and gets no report, so a species can publish fewer reports than it has curated datasets.

Counting writes each genome's canonical k-mers, where a sequence and its reverse complement count once, to a sorted, bucketed kset. The matrix keeps k-mers present in at least two genomes and absent from at least two, and merges k-mers present in exactly the same genomes into one pattern, the yes/no feature the learners use. Feature cards link each pattern to the genomes that carry it and to its gene and product in a reference genome.

The five learners are a Set Covering Machine (a few presence or absence rules joined by AND or OR, the model Kover uses), CART (a single decision tree), a random forest, gradient-boosted trees, and lasso logistic regression. The catalog baseline finds AMRFinderPlus reference genes with GRM's own k-mer matching, in place of the AMRFinderPlus program, and applies GRM's drug rules to the hits.

Each learner is cross-validated on two splits. The random split is stratified 5-fold cross-validation, the scheme most published numbers use, with each fold keeping the dataset's R/S mix. The lineage split holds out whole lineage clusters, so a test genome has no near-identical relative in training, and the gap between the two splits shows how much a score depends on such relatives. The lineage split is unavailable when whole lineages cannot fill outer and inner folds that each hold both classes, and the report records why.

Tuning and feature screening use only each outer fold's training genomes. From evaluation revision 7, the random forest is screened too, each fold's random draws are seeded from the fold so any held-out fit can be rebuilt, and the forest's subsamples and node features come from separate random streams. A boosted model with no trees or a lasso with no weights is a null fit: it has learned nothing beyond the class balance. From evaluation revision 14 it scores each isolate of its test fold at the share of resistant genomes in that fold's training set, called resistant at 0.5 or more, so every learner's metrics cover every isolate; revisions 7 to 13 made no call for those isolates and scored only the rest ([docs/equal-scoring.md](docs/equal-scoring.md)). The headline learner has the best random-split balanced accuracy (the mean of sensitivity and specificity) among learners scored on every isolate on every split that ran, so a learner missing from a lineage split that ran cannot be the headline. Saved predictions keep the 64-bit precision the metrics used.

A report also keeps what each outer fold's model reads (`fold_models`), so agreement between folds and with the final models can be counted, and for each genome the catalog groups whose own hit counts for the drug's rule (`rule_groups`). Both are optional additions that change no score, the reports of the revision 14 release do not have them, and a report gains them when it is next computed ([docs/report-rule-fields.md](docs/report-rule-fields.md)).

MIC regression predicts each isolate's MIC and is scored separately, by agreement within one doubling dilution. The catalog baseline and the simulated molecular panels describe which target genes or wild-type regions each assembled genome carries and how often that agrees with the lab call in this cohort.

## Reuse and derived outputs

Pipeline reuse checks, input fingerprints, `--drop-ksets`, `grm learn --force`, `grm curate` reruns, and full details of `grm intervals`, `grm baselines`, and `grm provenance` live in [docs/operations-reference.md](docs/operations-reference.md).

- **Intervals (`grm intervals`):** Writes `report/intervals/<id>.json` for every dataset `report/index.json` lists: for each learner row of each split, a two-sided 95% Clopper-Pearson (exact binomial) interval for VME, over the resistant isolates the row called (`fn` of `tp + fn`), and for ME, over the susceptible ones (`fp` of `tn + fp`). See [docs/operations-reference.md](docs/operations-reference.md#intervals) for zero-denominator handling, count pins, and SciPy `beta.ppf` reference checks.
- **Baselines (`grm baselines`):** Writes `report/baselines/<id>.json` and `report/baselines.json` with the constant training-fold MIC and lineage-identity baselines ([docs/baselines.md](docs/baselines.md), [docs/operations-reference.md](docs/operations-reference.md#baselines)).
- **Provenance (`grm fetch genome-dates`, `grm provenance`):** Writes `report/provenance/<id>.json` with source-proxy dates and record metadata, pinned to report and input SHA-256s ([docs/provenance-fields.md](docs/provenance-fields.md)); not part of public site releases yet.

## Checks

`.github/workflows/ci.yml` runs the Rust, Python, deployment-toolkit and web checks on every push and pull request. [docs/ci.md](docs/ci.md) lists its jobs and what it leaves out, such as the web app's browser check.

```bash
cargo test --locked
cargo fmt --check
npm --prefix web test
npm --prefix web run typecheck
npm --prefix web run build
printf '%s\n' scripts/run-pipeline.sh scripts/deployment/*.sh scripts/deployment/tests/*.sh | xargs -n1 bash -n
bash scripts/deployment/tests/dry-run.sh
bash scripts/deployment/tests/apply.sh
```

`scripts/check-results.py` and `scripts/metric-audit/` check saved reports apart from the Rust code that wrote them. Build `foldgen` and run their tests:

```bash
cargo build --release --locked --manifest-path scripts/metric-audit/foldgen/Cargo.toml
python3 -I scripts/metric-audit/tests/test_check_results.py
python3 -I scripts/metric-audit/tests/test_audit.py
python3 -I scripts/metric-audit/tests/test_scripts.py
```

`python3 -I scripts/check-results.py data` checks the saved dataset reports and writes `data/logs/result-validation-receipt.json`; `--out <file>` writes the receipt elsewhere. See [docs/operations-reference.md](docs/operations-reference.md#deployment-tests-and-result-checking) for deployment test skip rules, `--log` vs. `--skip-log`, and `audit.py --help` inputs.
