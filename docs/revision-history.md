# Revision history and archived notes

This page holds the detailed revision 5 through 14 release history, the retired Vast origin host notes, and the revision-5 `research.json` packaging details that were moved out of [README.md](../README.md) to keep the top-level readme focused on the current Revision 14 release.

## Evaluation revisions (5 to 14)

The current source is evaluation revision 14, and [genomicresistance.com](https://genomicresistance.com) serves its reports from the [`data-rev14-274-20261112` release](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev14-274-20261112).

- **Revision 5 and 6 (legacy / schema 1):** Schema 1 dataset reports and report index, where every learner scored every isolate. Revision-6 data predates kset digests, so running the current pipeline on revision-6 data recounts, rebuilds mismatched stages, and relearns. No public release holds revision-5 reports.
- **Revision 7 (schema 2, October 2026):** Introduced training-only feature screening for random forest, fold-seeded random draws, separate random streams for forest subsamples and node features, and schema 2 (representing null fits as no-calls). The revision-7 run covered the 13 species that revision 6 reported and recounted every genome, producing 237 reports for 239 curated datasets (`staphylococcus_aureus__methicillin` had 24 resistant genomes and `escherichia_coli__amoxicillin` had 3 labeled genomes, under the class minimum of 25). Three datasets were new: *Klebsiella pneumoniae* doripenem, *Mycobacterium tuberculosis* levofloxacin, and *Salmonella enterica* sulfisoxazole. No public release holds the revision-7 reports; the one copy named in this repository (`data/backups/revision7-published-20261025`) was on a local machine only.
- **Revision 8 and 9 ([`data-rev9-275-20261030`](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev9-275-20261030)):** Revision 8 fixed how the lasso stops so a fit hitting a solver cap can no longer be selected, scored, or saved; each report records how every lasso fit stopped, and any unconverged point makes the report invalid ([lasso-stopping.md](lasso-stopping.md)). Two datasets of the revision-8 run stopped on an unconverged lasso fit; revision 9 fixed the active-set solve that stalled there and changed nothing else. Revision 9 (275 reports over 23 species) was the last release the Vast host served (30–31 October 2026).
- **Revision 10:** Added the species-identity check that drops genomes whose sequence does not fit their labelled species ([species-check.md](species-check.md)), and kept k-mer-identical duplicate genomes whole in one outer test fold ([duplicate-isolation.md](duplicate-isolation.md)).
- **Revision 11:** Kept k-mer-identical duplicate groups whole in the random split's inner tuning folds.
- **Revision 12 ([`data-rev12-274-20261105`](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev12-274-20261105)):** Tuned final models on the same duplicate groups. Served on the site until 13 November 2026.
- **Revision 13 (unreleased):** Restricted each gene-panel target to the drugs and organisms its sources name (`src/panels.rs`), listing out-of-scope panels as not scored. Never run on its own, so no revision-13 reports exist.
- **Revision 14 ([`data-rev14-274-20261112`](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/data-rev14-274-20261112)):** Scores a null fit's isolates at its fold's training prevalence so no learner row leaves an isolate unscored, while keeping the revision-13 figure over fitted folds in `called_metrics` ([equal-scoring.md](equal-scoring.md)). Run on 12 November 2026 from revision 12's inputs across all 274 datasets: 162 reports differ from revision 12 only in revision number and run time, 39 in panel rows only, 48 only in rows of a learner with a null fit, and 25 in both; the best learner changed in 21 datasets. Served on the site since 13 November 2026.

Pulling new source recomputes nothing; each saved report keeps the revision and cohort it was computed with until a run recomputes it.

## Retired origin host (Vast + Cloudflare Tunnel)

Until 31 October 2026 one `grm serve` on a Vast container served `genomicresistance.com` behind a Cloudflare Tunnel. That host was destroyed on 31 October 2026 and replaced by static Cloudflare Worker + R2 hosting ([static-hosting.md](static-hosting.md)).

The host toolkit under `scripts/deployment/` stays in the repository for dry-run publication verification (`ops/cloudflare/scripts/publish.sh` calls it) and in case an origin host is ever needed again. [deployment.md](deployment.md) documents its workflow: local preview, publishing, staging, activation, rollback, container-restart restoration, and the tunnel connector. The publish, build, stage, activate, rollback, restore, and legacy-approval scripts change nothing under the deployment root without `--apply`; `connector.sh run` and `stop` act immediately, and `status.sh` and `smoke.sh` write only a temporary directory. Approvals and state files are strict single-record ASCII files written atomically and read back; starts are journaled with locked spawn markers so interrupted starts are detected and stopped, regular-file checks refuse links and FIFOs without blocking, and the tunnel token is read only by `cloudflared` via `--token-file`.

## Revision-5 `research.json` packaging

`research.json` is an optional dashboard file with two revision-5 analyses: GRM's MIC model beside a no-genome baseline, and selective calling (finding datasets where the genome alone can make the R/S call for a confident subset). `grm research` packages their saved outputs from `data/research/` without rerunning the analyses:

```bash
./target/release/grm research \
  --revision5-reports /path/to/revision5/datasets \
  --out data/report/research.json
```

It verifies the revision-5 reports against the snapshot hash list pinned in `src/research.rs` and records each input's SHA-256 digest. Each entry carries digests of the sample ids, labels, lineages, MICs, and held-out predictions its analysis read, and the dashboard shows its numbers only when the report on screen matches those digests.

The saved inputs are not published, so the command cannot be run from public files. It requires, under `<data>/research/`, `RESEARCH-SNAPSHOT-SHA256SUMS` (pinned in `src/research.rs`), `grm-research-mic-first-r1/` (`receipt.md`, `out/analysis/prediction_quality.csv`, `out/pairs/pairs.tsv.gz`, `out/splits_final/<id>.json`), `grm-research-selective-calling-r1/` (`receipt.md`, `out/species_drug_summary.csv`), and the revision-5 dataset reports. No public release holds the revision-5 reports, the revision 9 and 14 report trees hold no `research.json`, and the live site answers 404 for `/data/research.json` ([static-hosting.md](static-hosting.md)). Without `research.json`, the dashboard works normally.

## Legacy Python/Kover archive

The original Python/Kover application and its saved results live in the [legacy archive](https://github.com/editnori/Genomic-Resistance-Mapping-GRM-/releases/tag/legacy-archive). `grm report` reads 73 historical Kover rows from `reference/kover-summary.json` relative to the working folder (`--kover` overrides the path) and refuses if missing or unparseable; `grm pipeline` reads only from that default path and checks it before counting. The dashboard displays the matching row as historical context because Kover used different cohorts and evaluation procedures.
