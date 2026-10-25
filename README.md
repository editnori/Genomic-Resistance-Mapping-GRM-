# Genomic Resistance Mapping

GRM compares genomic predictors of antimicrobial susceptibility with reported laboratory phenotypes from public BV-BRC assemblies. Rust implements data ingestion, curation, k-mer counting, population structure, feature matrices, model evaluation, and report generation. The React dashboard explores the resulting datasets, gene and panel evidence, MIC predictions, and model performance.

The current source implements **evaluation revision 7** and writes dataset/index **schema 2**. The dashboard also reads schema 1. Updating this checkout does not recompute saved results: a report's evaluation revision, cohort, and audit determine which claims its numbers support.

Revision 7 has been run over the 13 species that revision 6 reported, recounting every genome. It produced 237 reports for the 239 curated datasets of those species. Staphylococcus aureus methicillin and Escherichia coli amoxicillin have no report: among each species' 2,000 selected genomes, methicillin has 24 resistant genomes and amoxicillin has 3 labeled genomes, under the class minimum of 25. Every revision-6 dataset has a revision-7 report, and three are new: Klebsiella pneumoniae doripenem, Mycobacterium tuberculosis levofloxacin, and Salmonella enterica sulfisoxazole. The result checker and metric audit under [Checks](#checks) pass on all 237 reports. This is internal validation on public BV-BRC data; no external, prospective, or clinical validation has been done.

The original Python/Kover application and saved results are preserved in [legacy/](legacy/). Those results use different cohorts and evaluation procedures, so they are historical context.

## Setup

Use Linux or WSL, Rust with edition 2024 support, and Node.js 22 with npm. Python 3 is needed only for the independent report-audit tools. Bash, `flock`, and `nice` are required by the pipeline launcher.

The Cargo configuration enables `target-cpu=native`; build on the machine that will run the executable.

```bash
cargo build --release --locked
npm --prefix web ci
npm --prefix web run build
```

Download metadata and the reference catalog, then curate datasets:

```bash
./target/release/grm fetch amr
./target/release/grm fetch genomes
./target/release/grm fetch catalog
./target/release/grm curate --cap 2000 --min-class 25
```

Counting needs the catalog for its gene calls. `grm count` and the pipeline stop when it does not load; `--allow-missing-catalog` counts without gene calls, which leaves reports without the catalog baseline and gene panel detections.

Curation uses unconflicted reported resistant/susceptible calls to determine eligibility and select genomes. Each dataset records its class minimum. MIC-inferred labels remain distinguishable in the source data and are excluded from binary model evaluation.

Drug names resolve through GRM's bundled vocabulary. To use an available SigDrug checkout, pass `--sigdrug-vocabulary /path/to/data/drug/compendium/facts/vocabulary` to `curate`, or set `GRM_SIGDRUG_VOCABULARY`. The report records the vocabulary source.

## Run

Start with one species:

```bash
bash scripts/run-pipeline.sh --species klebsiella_pneumoniae
./target/release/grm serve --port 8641
```

Open <http://localhost:8641>. The server binds to loopback and reads `web/dist` and `data/report`. For frontend development, run `npm --prefix web run dev`; `GRM_REPORT_DIR` can select a different report directory.

The launcher runs assembly download, counting, structure, matrix construction, learning, and report generation. It defaults to 64 Rayon threads, lowers process priority, sets the open-file limit to 8192, and locks `data/logs/pipeline.lock`. Adjust `RAYON_NUM_THREADS` to the host's available capacity. A 2,000-genome matrix merge needs more than the common 1,024-file-descriptor limit.

For a detached run:

```bash
mkdir -p data/logs
nohup setsid bash scripts/run-pipeline.sh \
  --species klebsiella_pneumoniae,staphylococcus_aureus \
  >data/logs/pipeline.log 2>&1 </dev/null &
tail -f data/logs/pipeline.log
```

The launcher records its PID in `data/logs/pipeline.pid` and its final exit code in `data/logs/pipeline.exit`. Check the PID and log as well as the exit file. It deletes intermediate ksets only after learning succeeds for a species, then refreshes the dashboard report. Matrices and dataset results remain on disk.

| Stage | Command | Output under the data root |
|---|---|---|
| Download | `grm fetch amr`, `genomes`, `catalog`, `fasta` | `bvbrc/`, `catalog/`, `fasta/` |
| Curate | `grm curate` | `curated/`, species genome lists |
| Count | `grm count --species <slug>` | Ksets, `counted.json` |
| Structure | `grm structure --species <slug>` | `structure.json` |
| Features | `grm matrix --species <slug>` | Species matrix |
| Evaluate | `grm learn --species <slug>` | Species/antibiotic result JSON |
| Report | `grm report` | `report/` |
| Research views | `grm research --revision5-reports <dir> --out <file>` | Guarded research JSON |

All commands accept `--data <directory>`. `grm learn --only meropenem` selects an antibiotic; `--force` requests recomputation. `--species` accepts comma-separated slugs. The pipeline uses every curated species when that option is omitted.

## Evaluation

The counter encodes canonical k-mers, sorts them, and writes bucketed ksets. The matrix combines k-mers with identical genome-presence patterns. The pipeline counts at k=31 and builds the matrix with a minimum support of two genomes. It reuses count, structure, and matrix stages that agree with each other and with that k and minimum support; it counts and builds again stages built by hand at another `--k` or `--min-support`. Models operate on those binary patterns; feature cards connect them to observed genomes and reference annotations.

The five fitted learners are a Set Covering Machine, CART, random forest, gradient-boosted trees, and lasso logistic regression. The catalog-containment baseline uses GRM's k-mer matching and drug rules with the AMRFinderPlus reference catalog. It does not run the AMRFinderPlus program.

Random splits and lineage-grouped splits measure different forms of generalization. The lineage split keeps related isolates together and is unavailable when its outer or inner partitions cannot support two-class evaluation. Its report records the reason.

Revision 7 adds training-only feature screening to the forest, derives seeds by fold, and separates subsampling and node-feature random streams. A null boosted-tree or lasso fit makes no call for that fold. Metrics use only called isolates and are null when none are called. The headline learner must call every isolate on every split that ran; a missing learner row in an available lineage split disqualifies it. Saved predictions retain the 64-bit precision used by the metrics.

Curation admits a dataset when each class reaches its `min_class` across all curated genomes. `grm curate` then selects at most `--cap` genomes per species (2,000 by default), and learning applies the same minimum again among those genomes. A dataset below it there is skipped with the log line `<id>: too few labeled genomes after download, skipped`, so a species can publish fewer reports than it has curated datasets.

MIC regression is a separate endpoint. Reference-gene and simulated molecular-panel results describe assembled genomes and agreement with lab calls in this cohort. They do not establish assay sensitivity, a causal resistance mechanism, or clinical treatment performance.

## Reusing results

The pipeline checks genome selection, order, k, matrix minimum support, catalog provenance, and per-genome kset digests before reusing count, structure, and matrix stages. Direct `grm learn` refuses stale selections or incompatible stages and names the command to rerun.

A learned report is reused only when its evaluation revision and input fingerprint match and no feature locus records a failed BV-BRC annotation download. A locus whose reference genome's annotation could not be downloaded has no gene or product and records the error in `annotation_error`; the next run recomputes that report. The fingerprint covers the matrix provenance, dataset class minimum, excluded inferred-label count, evaluation configuration, antibiotic/MIC breakpoint, lineage threshold, catalog version, and each evaluated genome's reported call, MIC, lineage, and gene calls in matrix order.

The fingerprint uses recorded matrix digests, so deleting ksets with `--drop-ksets` does not invalidate a result. Editing a kset outside `grm count` is not detected. Genome metadata, contig counts, and drug class are not fingerprinted; use `grm learn --force` after changing only those fields. A schema change alone does not trigger recomputation.

For data predating kset digests, run the pipeline once to recount and rebuild incompatible structure/matrix stages before learning. Revision-6 data is in that state, and every revision-6 report is recomputed because the evaluation revision changed. Name the species to relearn with `--species`; without it, the pipeline also downloads and counts every curated species that has no counts yet. Where `--drop-ksets` deleted the ksets, the recount reads each genome's FASTA file again. Rerunning `grm curate` first is a separate choice: curation now admits datasets and ranks genomes by reported calls only, so it can select other genomes, which the pipeline then downloads and counts. A `datasets.json` written before `min_class` was recorded reads as a class minimum of 25, the gate revision 6 applied. Until the run finishes, the dashboard shows revision-7 reports for finished species beside revision-6 reports for the rest, because the launcher rebuilds the report after each species.

The manual sequence follows. `grm matrix` opens every genome's kset at once, so raise the open-file limit first, as the launcher does; at the common limit of 1,024, `grm matrix` on a 1,100-genome species stops with "Too many open files".

```bash
ulimit -n 8192
./target/release/grm count --species klebsiella_pneumoniae
./target/release/grm structure --species klebsiella_pneumoniae
./target/release/grm matrix --species klebsiella_pneumoniae
./target/release/grm learn --species klebsiella_pneumoniae
./target/release/grm report
```

If curation rules or raw records changed, rerun curation first and fetch the newly selected assemblies. A dataset that no longer meets its recorded class minimum loses its saved result during learning. `grm report` publishes saved results; it does not refresh their evaluations. Report builds take turns through a lock on `data/report`, so a `grm report` started during another build waits for it to finish.

## Research views

The optional `research.json` exposes the existing MIC/no-genome-baseline and selective-calling analyses. Its generator verifies a pinned revision-5 snapshot and records source hashes. The dashboard compares cohort, labels, lineages, MICs, and prediction digests before showing applicable numbers. A matching dataset name alone is insufficient. Mismatched, withheld, missing, and no-model results retain distinct states.

With the corresponding research files under `data/research/` and the pinned reports available:

```bash
./target/release/grm research \
  --revision5-reports /path/to/revision5/datasets \
  --out data/report/research.json
```

This command packages existing research outputs; it does not rerun their analyses. The ordinary dashboard works without that optional file.

See [genotype-to-phenotype research](docs/research/genotype-to-phenotype.md) for diagnostic and surveillance applications, source evidence, and remaining validation requirements.

## Checks

```bash
cargo test --locked
cargo fmt --check
npm --prefix web test
npm --prefix web run typecheck
npm --prefix web run build
bash -n scripts/run-pipeline.sh
```

Independent report checks are maintained in `scripts/check-results.py` and `scripts/metric-audit/`. Build the fold reconstruction helper and run their tests:

```bash
cargo build --release --locked --manifest-path scripts/metric-audit/foldgen/Cargo.toml
python3 -I scripts/metric-audit/tests/test_check_results.py
python3 -I scripts/metric-audit/tests/test_audit.py
python3 -I scripts/metric-audit/tests/test_scripts.py
```

Check saved dataset results with `python3 -I scripts/check-results.py data`. Add `--log data/logs/pipeline.log` to verify the logged winners only when that log comes from a run that recomputed every report, such as the first revision-7 run: the checker requires a done line for every revision-7 report, and a run that reuses a report logs no done line for it. The detached-run command above replaces the log on each run, so keep a copy of that run's log. The checker writes `data/logs/result-validation-receipt.json`. The metric audit also reconstructs folds, recalculates metrics, and checks run provenance; `python3 -I scripts/metric-audit/audit.py --help` lists its required backup, log, and launch-time inputs. Use ordinary Python without optimization, because the checker uses assertions.

Data, reports, frontend bundles, dependency folders, and build outputs are ignored by Git. Building the dashboard does not generate study results.
