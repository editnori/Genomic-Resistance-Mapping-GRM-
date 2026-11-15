# Lasso scan at revision 14: the result of 15 November 2026

On 15 November 2026 the lasso of evaluation revision 14 was fitted again on the 18 datasets of the registered large-species scan list, [`scripts/lasso-scan/scan-list.json`](../../../scripts/lasso-scan/scan-list.json), every fit of each, and compared with the reports of the release `rev14-274-20261112`. One fit, MTB delamanid random fold 0, was also measured in depth beside the figures October's diagnosis gave for it. The section "Revision 14 scan" of [docs/lasso-stopping.md](../../../docs/lasso-stopping.md) gives what ran, what it found and its limits. This folder holds the files that section rests on.

This is a software check. The 57,573 equal scores from fits that keep a weight show that this build and the release compute the same lasso on these datasets; the other 1,045 equal scores are the training prevalence of the four folds that keep no weight, on both sides, and come from no fit. None of this says whether the lasso is well calibrated or accurate, and nothing clinical.

Every file here was copied unchanged from the output of the lane that ran the scan, keeping its paths, except the folder `review/`, which holds the files of an independent review of the scan (see [its README](review/README.md)). `SHA256SUMS` gives the sha256 of each, every line of the scan's files as it stands in that lane's checksum list; `sha256sum -c SHA256SUMS`, run in this folder, checks them.

## How it was made

| | |
|---|---|
| Code | repository `main` at commit `852ebc91c3f649d3980b1d9a6844e1daec5df2d2`, exported with `git archive`. There `src/`, `Cargo.lock`, `Cargo.toml`, `.cargo/` and `scripts/lasso-scan/` are the same as at `469a4f9`, the commit the release was built from |
| Build | Rust 1.98.1, `cargo build --release --locked --offline --bin grm --example lasso-scan`; `lasso-scan`'s sha256 is in `results/scan/meta.json` |
| Inputs | matrices built with `grm count`, `grm structure` and `grm matrix` at default settings, k = 31, from the public BV-BRC data root whose checksum list the release notes name (`scripts/build-species.sh`) |
| Release compared | `rev14-274-20261112`, its reports under `datasets/` |
| Machine | AMD Ryzen 9 9950X, 16 cores, 31 GiB, `RAYON_NUM_THREADS=16` (`results/scan/run-all.machine`) |
| Plan | `PLAN.md`, fixed at 14:11:00 UTC before any scan result was read; `results/plan.sha256` holds its sha256 and that time |
| Runs | pass 1 from 14:21:15 to 14:30:05 UTC; the second timing of MTB delamanid from 14:30:06 to 14:30:25 UTC; pass 2, with the instrumented copy, from 14:30:25 UTC; the planned one-fit measurement from 14:58:52 to 15:47:02 UTC |

## Files

| File | What it is | Written by |
|---|---|---|
| `PLAN.md` | The plan: datasets and order, settings, the scan's convergence rule, the comparison with the release, when a released figure needs a caveat, and the one-fit measurement. As hashed | the lane, before the scan |
| `results/plan.sha256` | `PLAN.md`'s sha256 and the time it was recorded | `sha256sum PLAN.md` |
| `results/scan/summary.json` | Pass 1, per dataset: input fingerprint, samples and features against the release's, each fit's audit, penalty points, IRLS steps, points with a short solve, weights and seconds, the budgets at h = 1, and whether the dataset converged by the scan's rule | `scripts/scan14.py`, run by `scripts/run-all.sh` |
| `results/scan/meta.json` | Pass 1's binary and its sha256, threads, CPU and dataset order | `scripts/scan14.py` |
| `results/scan/run-all.machine` | Time, `nproc`, CPU, load average and memory before pass 1, the load after it, and pass 2's start | `scripts/run-all.sh` |
| `results/scan/scan-time2.machine` | Time, `nproc`, CPU and load average around the second MTB delamanid run | `scripts/time-delamanid.sh` |
| `results/scan/mycobacterium_tuberculosis__delamanid.stderr.log.gz` | Pass 1's progress lines for MTB delamanid, one per solved penalty point, with its seconds; the path times are their sums. Gzipped by the lane; `zcat` reads it | the `lasso-scan` binary, under `scripts/scan14.py` |
| `results/scan/events.sha256` | sha256 of each dataset's uncompressed fit events, which are not kept here | the lane |
| `results/scan-time2/summary.json`, `meta.json`, `mycobacterium_tuberculosis__delamanid.stderr.log.gz` | The same three for the second run of MTB delamanid alone | `scripts/time-delamanid.sh` |
| `results/compare.json` | Per dataset and split: scores equal to the release's, the largest score change, calls changed, and recall, log loss, balanced accuracy, VME, ME and AUROC from the release and from the scan; whether pass 2's fits equal pass 1's | `python3 -I scripts/compare.py --scan $S/scan --scan2 $S/scan-dump --dump $S/dump --release $REL --out results/compare.json` |
| `results/tables.md` | The generated tables: per dataset, per split, and the one-fit table per penalty point | `scripts/tables.py` |
| `results/go-no-go-h1.txt` | The registered go/no-go and budgets computed at h = 1, which are exploratory | the rule's `decide()` from `scripts/lasso-scan/scan.py`, as the file's first line says |
| `results/src-dump.diff` | The only change in the instrumented copy pass 2 ran: one function added to `src/learn/eval.rs` that writes each fit's partition, path, refit and scores | `diff -ru` of the two exported trees |
| `results/dump.sha256` | sha256 of pass 2's 198 dump files, 49,691,720 bytes in all, which are not kept here | the lane |
| `results/item5/item5.json`, `item5.log` | The planned one-fit measurement: near twins, then per path point grm's objective against an independent solve by Clarabel 0.9.0 through CVXPY 1.6.0, the KKT residuals and the planned duality bound; the selection, refit and held-out calls; the second solver | `scripts/item5-run2.sh`: 15 workers of `scripts/item5-solve.py`, then `scripts/item5.py` |
| `results/item5/item5-polish.json`, `item5-polish.log` | Exploratory, chosen at 14:55 UTC after point 11's planned bound was read: a certified optimum at each path point, and grm's distance above it. Its dual point lies outside the dual box at 13 of the 30 points, by a factor of at most 1 + 3.18e-13, and certifies only after a final shrink into the box, which the review made (`review/results/item5check.json`) | `scripts/item5-polish.py` |
| `results/item5/item5-copies.json`, `item5-copies.log` | Exploratory, chosen at 15:31 UTC after the validation losses were read: those losses with the independent solution's weight on each group of exact copies moved to the group's lowest index | `scripts/item5-copies.py` |
| `results/item5/run2.machine`, `run.rc`, `saga-timebox.txt` | The one-fit run's times and load, its exit codes, and when the second solver was stopped at its 30-minute box | `scripts/item5-run2.sh`; the lane stopped the second solver |
| `scripts/` | The lane's scripts, as run: `build-species.sh` (one species' matrices), `run-all.sh` (pass 1, the second timing, pass 2), `scan14.py` (the driver), `time-delamanid.sh`, `compare.py`, `item5-solve.py`, `item5-run2.sh`, `item5.py`, `item5-polish.py`, `item5-copies.py` and `tables.py` | the lane |
| `review/` | The independent review's result files that the page quotes, the scripts that wrote them, its plan and that plan's sha256, and a README that says what each is and which sentences of `PLAN.md`, `results/tables.md` and this README the review narrowed | copied from the review's output, except its README |
| `SHA256SUMS` | sha256 of every other file here except this README | lines copied from the lane's checksum list; for `review/`, lines copied from the review's checksum list with `review/` before each path, and one computed for its README |

In the commands, `S` is the lane's scratch folder and `REL` the unpacked release. The scripts hold that machine's absolute paths. They are kept as the record of what ran, not as a kit to run as they stand. The lane's report does not record the commands that wrote `results/go-no-go-h1.txt`, `results/scan/events.sha256` and `results/dump.sha256`.

## Not kept

- The matrices of the five species, built in the lane's scratch, and their build logs.
- Pass 2's 198 dump files, 49,691,720 bytes, listed in `results/dump.sha256`. One of them, the full dump of MTB delamanid random fold 0, is what `item5.py`, `item5-polish.py` and `item5-copies.py` read; the instrumented build of `results/src-dump.diff` writes it again.
- The fit events and progress lines of the other 17 datasets, pass 2's summary and events, the 30 independent solutions of the one-fit measurement, its workers' logs, its CPU samples, and its first run, which the lane stopped at point 12.
- The lane's report and the solver environment's wheels (Clarabel 0.9.0, CVXPY 1.6.0, numpy 2.1.3 and scikit-learn 1.5.2, each checked by sha256 when installed).

All of these stayed on the machine that ran the scan. No file in this repository keeps them.
