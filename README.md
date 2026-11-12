# Genomic Resistance Mapping

GRM measures how well a bacterial isolate's genome predicts its laboratory antimicrobial susceptibility result. It uses public assemblies from BV-BRC (the Bacterial and Viral Bioinformatics Resource Center) that have a reported resistant (R) or susceptible (S) call. The Rust pipeline downloads and curates them, extracts each genome's 31-base DNA substrings (k-mers), groups the genomes into lineages, builds a feature matrix per species, cross-validates five learners on each species and antibiotic, and writes a report per dataset. The React dashboard shows each dataset's gene and panel evidence, MIC (minimum inhibitory concentration) predictions, and model performance.

## Setup

Use Linux or WSL with Rust 1.89 or newer, Node.js 22 and npm, and Bash, `flock`, and `nice` for the pipeline launcher. Python 3 is needed only for the report checks.

`.cargo/config.toml` sets `target-cpu=x86-64-v3`, so the learners' bit counting uses the CPU's popcount instruction and the binary runs on any x86-64 CPU with AVX2 (Intel Haswell, AMD Zen and newer). For a long run on the machine that builds it, `RUSTFLAGS="-C target-cpu=native" cargo build --release --locked` builds for that CPU instead; that binary may not start on another CPU.

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

A Cloudflare Worker serves genomicresistance.com from one R2 bucket, with no origin process. [docs/static-hosting.md](docs/static-hosting.md) covers the bucket, the release scripts and their checks. On 11 November 2026 the site answered HTTP 200 with the revision-12 report, 274 datasets over 23 species.

Until 31 October 2026 one `grm serve` on a Vast container served the site behind a Cloudflare Tunnel. That host was destroyed that day. [docs/deployment.md](docs/deployment.md) covers that host's toolkit, which stays in the repository in case the site needs a host again: local preview from the verified backup, publishing a report copy, building and staging a release, activation, rollback, restoring after a container restart, and the tunnel connector. The publish, build, stage, activate, rollback, restore and legacy-approval scripts in `scripts/deployment/` change nothing under the deployment root without `--apply`; `connector.sh run` and `stop` act at once, and `status.sh` and `smoke.sh` write only a temporary directory. Before serving a pair, the scripts verify the build and report trees against their SHA-256 manifests and require a stage approval bound to both manifests (`status.sh` does not verify manifests). Approvals and state files are strict one-record ASCII files, written atomically and read back. Every server start is journaled, and every process it starts holds a locked spawn marker from the fork on, so an interrupted or killed start is found and stopped, and its journal is kept while the lock shows a survivor the scripts cannot see. A rollback identifies a damaged current server by its recorded start (PID, start time, token and marker), not by its executable. Records, locks and logs are opened only when they are regular files, so a link or FIFO is refused without blocking. The tunnel token is read only by cloudflared, through `--token-file`.

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
| Intervals | `grm intervals` | `report/intervals/<id>.json` |
| Provenance | `grm fetch genome-dates`, then `grm provenance` | `bvbrc/genome_dates.jsonl` and its manifest, `report/provenance/<id>.json` |
| Predict | `grm predict --report <file> --fasta <files>` | JSON lines on standard output, or `--out` ([docs/prediction-contract.md](docs/prediction-contract.md)) |

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

MIC regression predicts each isolate's MIC and is scored separately, by agreement within one doubling dilution. The catalog baseline and the simulated molecular panels describe which target genes or wild-type regions each assembled genome carries and how often that agrees with the lab call in this cohort.

## Results

The site serves the revision-12 run: 274 reports over 23 species, from 277 curated datasets. Three datasets have no report because one class has fewer than 25 genomes, the class minimum: Escherichia coli amoxicillin, Klebsiella pneumoniae cefotaxime and Staphylococcus aureus methicillin. The published report lists them in `SKIPS.json`.

Revision 14 was run for the same 274 datasets on 12 November 2026, from revision 12's inputs. Of the 274 reports, 162 differ from revision 12 only in the revision number and the run time, 39 in panel rows only, 48 only in the rows of a learner with a null fit, and 25 in both. The best learner changed in 21 datasets. The reports and their release notes are in the [`data-rev14-274-20261112` pre-release](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev14-274-20261112). An independent review recomputed the comparison with revision 12 from the files, and a second review checked the release notes. The site does not serve these reports yet.

The revision-7 run covered the 13 species that revision 6 reported and recounted every genome. It produced 237 reports for those species' 239 curated datasets. Staphylococcus aureus methicillin and Escherichia coli amoxicillin have no report: among each species' 2,000 selected genomes, methicillin has 24 resistant genomes and amoxicillin has 3 labeled genomes, under the class minimum of 25. Every revision-6 dataset has a revision-7 report, and three datasets are new: Klebsiella pneumoniae doripenem, Mycobacterium tuberculosis levofloxacin, and Salmonella enterica sulfisoxazole. The result checker and metric audit under [Checks](#checks) pass on all 237 reports.

This is internal validation: cross-validation on public BV-BRC genomes scored against their reported calls. No external, prospective, or clinical validation has been done.

The original Python/Kover application and its saved results are in the [legacy archive](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/legacy-archive). `grm report` reads 73 historical Kover rows from `reference/kover-summary.json`, relative to the working folder (`--kover` takes another path), and refuses if that file is missing or doesn't parse. `grm pipeline` has no `--kover`: it reads the file only from that path, so run it from the source folder, and it checks the file before counting. The dashboard shows the one for the same species and antibiotic as historical context, because Kover used other cohorts and evaluation procedures.

## Reuse

The current source is evaluation revision 14. The site serves revision 12. Revision 14 has reports for the same 274 datasets, in the [`data-rev14-274-20261112` pre-release](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev14-274-20261112), and the site does not serve them yet. Revision 13 was never run on its own, so no revision-13 reports exist. The detailed results above are from revision 7. Revision 9, with 275 reports over 23 species, was the last release the Vast host served, from 30 October 2026 until that host was destroyed on 31 October 2026 ([Deploy](#deploy)). Revision 8 fixes how the lasso stops: a fit that stops at a solver cap can no longer be selected, scored or saved, each report records how every lasso fit stopped, and any unconverged point makes the report invalid ([docs/lasso-stopping.md](docs/lasso-stopping.md)). Two datasets of the revision-8 run stopped on an unconverged lasso fit; revision 9 fixes the active-set solve that stalled there and changes nothing else. Revision 10 adds a species identity check that drops genomes whose sequence does not fit their labelled species ([docs/species-check.md](docs/species-check.md)), and keeps k-mer-identical duplicate genomes whole in one outer test fold, so no genome is scored by a fit that trained on its own copy ([docs/duplicate-isolation.md](docs/duplicate-isolation.md)). Revision 11 keeps those duplicates whole in the random split's inner tuning too, so no tuning partition puts a genome and its copy on both sides. Revision 12 tunes the final models on the same duplicate groups. Revision 13 reads each gene-panel target only for the drugs and organisms its sources name ([src/panels.rs](src/panels.rs)), so panel rows change and a panel outside its organisms is listed as not scored. Revision 14 scores a null fit's isolates at its fold's training prevalence, so no learner row leaves an isolate unscored, and keeps the revision-13 figure over the fitted folds in `called_metrics` ([docs/equal-scoring.md](docs/equal-scoring.md)). It still writes schema 2 for dataset reports and the report index, which can hold null scores and metrics. The dashboard also reads schema 1, written by revision 6 and earlier, where every learner scored every isolate. Pulling new source recomputes nothing, so each saved report keeps the revision and cohort it was computed with until a run recomputes it.

Before reusing count, structure, and matrix, the pipeline checks the genome selection and order, k = 31, the matrix minimum support of 2, the catalog version behind the gene calls, and each genome's kset digest (a hash of its k-mers). A mismatch, including stages built by hand at another `--k` or `--min-support`, reruns counting, which reuses ksets on disk at k = 31 and recounts the rest from FASTA, then rebuilds the stages that no longer match. Direct `grm learn` refuses the same mismatches and names the command to rerun.

A learned report is reused only when its evaluation revision and input fingerprint match, its evaluation audit is valid, and no feature records a failed BV-BRC annotation download. A report whose audit is invalid is written for review, but `grm learn` stops with an error instead of logging it as done, and the next run recomputes it. Such a feature has no gene or product and records the error in `annotation_error`, so the next run recomputes its report. The fingerprint hashes the matrix digest, the antibiotic and its MIC breakpoint, the dataset's class minimum, the count of excluded MIC-inferred labels, the evaluation settings, the lineage threshold, the catalog version, and each evaluated genome's reported call, MIC, lineage, and gene calls in matrix order. Because it uses the kset digests recorded in the matrix, deleting ksets with `--drop-ksets` keeps results current, and editing a kset outside `grm count` goes unnoticed. Genome metadata, contig counts, and drug class change only descriptive fields and are left out, so run `grm learn --force` after changing only those. A schema change alone recomputes nothing.

Revision-6 data predates kset digests, and every revision-6 report is recomputed because the evaluation revision changed, so one pipeline run recounts, rebuilds the stages that no longer match, and relearns. Name the species with `--species`; without it, the pipeline also downloads and counts every curated species that has no counts yet. Until the run finishes, the dashboard shows the new revision's reports for finished species beside revision-6 reports for the rest.

Rerunning `grm curate` is a separate choice. Curation counts only reported calls when it admits datasets and ranks genomes, and older curated files also counted MIC-inferred labels, so a new curation can select other genomes for the pipeline to download and count. A `datasets.json` written before `min_class` was recorded reads as a class minimum of 25, the gate revision 6 applied. After raw records or curation rules change, rerun `grm curate` before the pipeline. A dataset that no longer meets its recorded class minimum loses its saved report when learning reaches it.

## Research

`research.json` is an optional dashboard file with two revision-5 analyses: GRM's MIC model beside a baseline that predicts the MIC without the genome, and selective calling, which finds the datasets where the genome alone can make the R/S call for a confident subset of isolates. `grm research` packages their saved outputs from `data/research/` without rerunning the analyses:

```bash
./target/release/grm research \
  --revision5-reports /path/to/revision5/datasets \
  --out data/report/research.json
```

It verifies the revision-5 reports against the snapshot hash list pinned in the source and records each input's hash. Each entry carries SHA-256 digests of the sample ids, labels, lineages, MICs, and held-out predictions its analysis read, and the dashboard shows its numbers only when the report on screen gives the same digests, because a later revision can keep a dataset id and change its genomes. Mismatched entries, missing entries, reports without an MIC model, and entries withheld because their MIC model read disk-diffusion zone diameters as MICs each show as their own state. Without `research.json`, the dashboard works as usual.

## Intervals

`grm intervals` writes `report/intervals/<id>.json` for every dataset `report/index.json` lists: for each learner row of each split, a two-sided 95% Clopper-Pearson (exact binomial) interval for VME, over the resistant isolates the row called (`fn` of `tp + fn`), and for ME, over the susceptible ones (`fp` of `tn + fp`). It is arithmetic on the counts the report already holds; it fits nothing and changes no report file, so a report built without it is byte for byte what it was. A zero denominator gives no interval. Each file records its report's SHA-256, evaluation revision and the counts it used, and the dashboard shows an interval under its VME or ME only where those counts equal the report's, so a file left from an earlier build shows nothing. Run it after `grm report` and before the report tree's manifest is written, so the publication ships the files. Without them, the dashboard works as before. The intervals describe sampling of these isolates under internal cross-validation; they are not calibration, clinical bounds or external validation. The method and reference values are those of the revision-7 offline run (`data/research/grm-all23-saved-intervals-20261027-kimi`, SciPy `beta.ppf`), which `src/intervals.rs`'s tests reproduce to 1e-12.

## Provenance

`grm fetch genome-dates` fetches BV-BRC's `date_inserted` and `completion_date` for every genome in `species/*/genomes.txt` into `bvbrc/genome_dates.jsonl`, with `bvbrc/genome_dates.manifest.json` pinning it (the lists asked for, the genomes BV-BRC returned no record for, and the file's SHA-256). `grm curate` writes `curated/provenance.jsonl`, each genome and drug pair's record details: methods, standards, record count, MIC spread, the earliest and latest record `date_inserted`, and the sets of testing-standard years, PMIDs and sources. `grm provenance` then writes `report/provenance/<id>.json` for every dataset `report/index.json` lists, pinned to its report's SHA-256 and evaluation revision and to both inputs' SHA-256, and `grm provenance --check` fails when any pin no longer holds. None of it changes a report, a curated file `grm curate` wrote before, an input fingerprint or a learner's input; the dates are source proxies and never predictors. [docs/provenance-fields.md](docs/provenance-fields.md) lists what each field reaches. The dashboard doesn't read these files, and they aren't part of a publication yet: build a report tree bound for release without `report/provenance/`.

[docs/research/genotype-to-phenotype.md](docs/research/genotype-to-phenotype.md) covers diagnostic and surveillance uses, the published evidence, and the validation a clinical claim would still need.

[docs/outside-study-results.md](docs/outside-study-results.md) gives the result of a retrospective test on newly deposited studies from the same database: each revision-12 report's final model, applied unchanged to BV-BRC genomes added after the cohort was built. [docs/outside-study-results-by-species.md](docs/outside-study-results-by-species.md) reads it one species at a time. [docs/outside-study-test.md](docs/outside-study-test.md) is the plan, frozen before any of those genomes was fetched. The output files are in `reference/outside-study/`.

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

`bash -n` parses only its first file, so `xargs -n1` runs it once per script. The deployment tests run in temporary directories; `apply.sh` needs `cc`, starts stub servers on free loopback ports, injects signals at each step of a start and in the forked child before it execs, and takes about 15 to 20 minutes on a loaded host. Both suites count their assertions and fail when the count differs from the expected total. Their only skips depend on the environment, and each prints a `skip` line and lowers the expected total: the root-only cases in both suites, and in `apply.sh` the user and PID namespace cases where `unshare` cannot create those namespaces. A mutation harness, kept outside the repository, checked that the suites catch each reviewed deployment defect ([docs/deployment.md](docs/deployment.md#checks)).

`scripts/check-results.py` and `scripts/metric-audit/` check saved reports apart from the Rust code that wrote them. The checker recomputes each learner's confusion counts and balanced accuracy from the saved scores; the metric audit also rebuilds the folds with the `foldgen` helper, recalculates the metrics, and checks the run against its backup and log. Build `foldgen` and run their tests:

```bash
cargo build --release --locked --manifest-path scripts/metric-audit/foldgen/Cargo.toml
python3 -I scripts/metric-audit/tests/test_check_results.py
python3 -I scripts/metric-audit/tests/test_audit.py
python3 -I scripts/metric-audit/tests/test_scripts.py
```

`python3 -I scripts/check-results.py data` checks the saved dataset reports and writes `data/logs/result-validation-receipt.json`; `--out <file>` writes the receipt elsewhere, such as outside a read-only data root. Run it without `-O`, which would skip its assertions.

`--log <pipeline log>` also requires every report from revision 7 on to have a `done` line in that log naming its headline learner. A run that reuses a report logs no `done` line for it, so use `--log` only with the log of a run that recomputed every report, such as the first revision-7 run, and keep a copy, because the detached run replaces `data/logs/pipeline.log` each time. After a run that reused reports, use `--skip-log <pipeline log>`, which reads only the skip lines.

Both checkers take the current datasets from `data/curated/datasets.json` and list a result for a dataset no longer curated as "not in current curation" without reading it. A species is checked when it has results or both `genomes.txt` and `counted.json`, and the metric audit also checks every curated species matching a species glob other than `*`. Other curated species are listed as "not run" and do not fail. Every curated dataset of a checked species needs a result unless the skip log (`--skip-log`, or else `--log`) has its exact `too few labeled genomes after download, skipped` line with no later `done` or `recomputing report` line, so with neither option a missing result fails. A skip never exempts an existing report, and both receipts list the exclusions and skips.

`python3 -I scripts/metric-audit/audit.py --help` lists the metric audit's inputs: the evaluation revision, the report backup taken before launch, the run's pipeline log, and its UTC launch time.
