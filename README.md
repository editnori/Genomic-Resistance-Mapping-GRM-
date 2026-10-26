# Genomic Resistance Mapping

GRM measures how well a bacterial isolate's genome predicts its laboratory antimicrobial susceptibility result. It uses public assemblies from BV-BRC (the Bacterial and Viral Bioinformatics Resource Center) that have a reported resistant (R) or susceptible (S) call. The Rust pipeline downloads and curates them, extracts each genome's 31-base DNA substrings (k-mers), groups the genomes into lineages, builds a feature matrix per species, cross-validates five learners on each species and antibiotic, and writes a report per dataset. The React dashboard shows each dataset's gene and panel evidence, MIC (minimum inhibitory concentration) predictions, and model performance.

## Setup

Use Linux or WSL with Rust 1.89 or newer, Node.js 22 and npm, and Bash, `flock`, and `nice` for the pipeline launcher. Python 3 is needed only for the report checks.

`.cargo/config.toml` sets `target-cpu=native` so the learners' bit counting uses the CPU's popcount instruction, which means the binary should be built on the machine that runs it.

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
| Report | `grm report` | `report/` |
| Research | `grm research --revision5-reports <dir> --out <file>` | the `--out` file |

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

Tuning and feature screening use only each outer fold's training genomes. From evaluation revision 7, the random forest is screened too, each fold's random draws are seeded from the fold so any held-out fit can be rebuilt, and the forest's subsamples and node features come from separate random streams. A boosted model with no trees or a lasso with no weights is a null fit: it scores only the 0.5 prior, so it makes no call for its test fold. Metrics cover only called isolates and are null when none are called. The headline learner has the best random-split balanced accuracy (the mean of sensitivity and specificity) among learners that call every isolate on every split that ran, so a learner missing from a lineage split that ran cannot be the headline. Saved predictions keep the 64-bit precision the metrics used.

MIC regression predicts each isolate's MIC and is scored separately, by agreement within one doubling dilution. The catalog baseline and the simulated molecular panels describe which target genes or wild-type regions each assembled genome carries and how often that agrees with the lab call in this cohort.

## Results

The revision-7 run covered the 13 species that revision 6 reported and recounted every genome. It produced 237 reports for those species' 239 curated datasets. Staphylococcus aureus methicillin and Escherichia coli amoxicillin have no report: among each species' 2,000 selected genomes, methicillin has 24 resistant genomes and amoxicillin has 3 labeled genomes, under the class minimum of 25. Every revision-6 dataset has a revision-7 report, and three datasets are new: Klebsiella pneumoniae doripenem, Mycobacterium tuberculosis levofloxacin, and Salmonella enterica sulfisoxazole. The result checker and metric audit under [Checks](#checks) pass on all 237 reports.

This is internal validation: cross-validation on public BV-BRC genomes scored against their reported calls. No external, prospective, or clinical validation has been done.

The original Python/Kover application and its saved results are in the [legacy archive](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/legacy-archive). `grm report` reads 73 historical Kover rows from `reference/kover-summary.json`, and the dashboard shows the one for the same species and antibiotic as historical context, because Kover used other cohorts and evaluation procedures.

## Reuse

The current source is evaluation revision 7 and writes schema 2 for dataset reports and the report index, which can hold null scores and metrics. The dashboard also reads schema 1, written by revision 6 and earlier, where every learner scored every isolate. Pulling new source recomputes nothing, so each saved report keeps the revision and cohort it was computed with until a run recomputes it.

Before reusing count, structure, and matrix, the pipeline checks the genome selection and order, k = 31, the matrix minimum support of 2, the catalog version behind the gene calls, and each genome's kset digest (a hash of its k-mers). A mismatch, including stages built by hand at another `--k` or `--min-support`, reruns counting, which reuses ksets on disk at k = 31 and recounts the rest from FASTA, then rebuilds the stages that no longer match. Direct `grm learn` refuses the same mismatches and names the command to rerun.

A learned report is reused only when its evaluation revision and input fingerprint match and no feature records a failed BV-BRC annotation download. Such a feature has no gene or product and records the error in `annotation_error`, so the next run recomputes its report. The fingerprint hashes the matrix digest, the antibiotic and its MIC breakpoint, the dataset's class minimum, the count of excluded MIC-inferred labels, the evaluation settings, the lineage threshold, the catalog version, and each evaluated genome's reported call, MIC, lineage, and gene calls in matrix order. Because it uses the kset digests recorded in the matrix, deleting ksets with `--drop-ksets` keeps results current, and editing a kset outside `grm count` goes unnoticed. Genome metadata, contig counts, and drug class change only descriptive fields and are left out, so run `grm learn --force` after changing only those. A schema change alone recomputes nothing.

Revision-6 data predates kset digests, and every revision-6 report is recomputed because the evaluation revision changed, so one pipeline run recounts, rebuilds the stages that no longer match, and relearns. Name the species with `--species`; without it, the pipeline also downloads and counts every curated species that has no counts yet. Until the run finishes, the dashboard shows revision-7 reports for finished species beside revision-6 reports for the rest.

Rerunning `grm curate` is a separate choice. Curation counts only reported calls when it admits datasets and ranks genomes, and older curated files also counted MIC-inferred labels, so a new curation can select other genomes for the pipeline to download and count. A `datasets.json` written before `min_class` was recorded reads as a class minimum of 25, the gate revision 6 applied. After raw records or curation rules change, rerun `grm curate` before the pipeline. A dataset that no longer meets its recorded class minimum loses its saved report when learning reaches it.

## Research

`research.json` is an optional dashboard file with two revision-5 analyses: GRM's MIC model beside a baseline that predicts the MIC without the genome, and selective calling, which finds the datasets where the genome alone can make the R/S call for a confident subset of isolates. `grm research` packages their saved outputs from `data/research/` without rerunning the analyses:

```bash
./target/release/grm research \
  --revision5-reports /path/to/revision5/datasets \
  --out data/report/research.json
```

It verifies the revision-5 reports against the snapshot hash list pinned in the source and records each input's hash. Each entry carries SHA-256 digests of the sample ids, labels, lineages, MICs, and held-out predictions its analysis read, and the dashboard shows its numbers only when the report on screen gives the same digests, because a later revision can keep a dataset id and change its genomes. Mismatched entries, missing entries, reports without an MIC model, and entries withheld because their MIC model read disk-diffusion zone diameters as MICs each show as their own state. Without `research.json`, the dashboard works as usual.

[docs/research/genotype-to-phenotype.md](docs/research/genotype-to-phenotype.md) covers diagnostic and surveillance uses, the published evidence, and the validation a clinical claim would still need.

## Checks

```bash
cargo test --locked
cargo fmt --check
npm --prefix web test
npm --prefix web run typecheck
npm --prefix web run build
bash -n scripts/run-pipeline.sh
```

`scripts/check-results.py` and `scripts/metric-audit/` check saved reports apart from the Rust code that wrote them. The checker recomputes each learner's confusion counts and balanced accuracy from the saved scores; the metric audit also rebuilds the folds with the `foldgen` helper, recalculates the metrics, and checks the run against its backup and log. Build `foldgen` and run their tests:

```bash
cargo build --release --locked --manifest-path scripts/metric-audit/foldgen/Cargo.toml
python3 -I scripts/metric-audit/tests/test_check_results.py
python3 -I scripts/metric-audit/tests/test_audit.py
python3 -I scripts/metric-audit/tests/test_scripts.py
```

`python3 -I scripts/check-results.py data` checks the saved dataset reports and writes `data/logs/result-validation-receipt.json`; `--out <file>` writes the receipt elsewhere, such as outside a read-only data root. Run it without `-O`, which would skip its assertions.

`--log <pipeline log>` also requires every revision-7 report to have a `done` line in that log naming its headline learner. A run that reuses a report logs no `done` line for it, so use `--log` only with the log of a run that recomputed every report, such as the first revision-7 run, and keep a copy, because the detached run replaces `data/logs/pipeline.log` each time. After a run that reused reports, use `--skip-log <pipeline log>`, which reads only the skip lines.

Both checkers take the current datasets from `data/curated/datasets.json` and list a result for a dataset no longer curated as "not in current curation" without reading it. A species is checked when it has results or both `genomes.txt` and `counted.json`, and the metric audit also checks every curated species matching a species glob other than `*`. Other curated species are listed as "not run" and do not fail. Every curated dataset of a checked species needs a result unless the skip log (`--skip-log`, or else `--log`) has its exact `too few labeled genomes after download, skipped` line with no later `done` or `recomputing report` line, so with neither option a missing result fails. A skip never exempts an existing report, and both receipts list the exclusions and skips.

`python3 -I scripts/metric-audit/audit.py --help` lists the metric audit's inputs: the evaluation revision, the report backup taken before launch, the run's pipeline log, and its UTC launch time.
