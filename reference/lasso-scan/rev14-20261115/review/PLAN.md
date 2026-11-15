# revlasso: plan for the three-dataset replication (brief point 3.4)

Fixed before any of my runs. Its sha256 goes to `results/plan.sha256` before the first run starts. Anything decided after a result is read is labelled exploratory in `REPORT.md`.

## What is compared

My own build of `852ebc9` (`/home/user/scratch/revlasso/target`, Rust 1.98.1, the repository's `x86-64-v3`, `cargo build --release --locked --offline --bin grm --example lasso-scan`, no network) run on the lane's matrices (`/home/user/scratch/lassoscan/data`, read only), against the lane's pass 1 events (`/home/user/grm-ops/program-20261114/lassoscan/results/scan/`), the lane's pass 2 dumps (`/home/user/scratch/lassoscan/dump/`) and the release (`/home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112/datasets/`).

Two binaries of mine:

1. `lasso-scan`, unmodified, from my build.
2. An instrumented copy for scores: `852ebc9` exported again into `/home/user/scratch/revlasso/src-dump`, with the lane's patch `results/src-dump.diff` applied (one function added to `src/learn/eval.rs` that writes each fit's partition, path, refit and scores when `GRM_LASSO_DUMP` is set). I read the patch before building it: it clones values and writes them, and changes no value the fit returns. Built into `/home/user/scratch/revlasso/target-dump`.

Each binary runs once per dataset, one dataset per call, `RAYON_NUM_THREADS=8`, under `unshare -rn` (no network), outputs in `/home/user/scratch/revlasso/runs/`.

## Datasets, by a rule stated before the runs

From the registered list `scripts/lasso-scan/scan-list.json`, using only the release reports' fields:

- the dataset with the most features in the release's `counts`: `escherichia_coli__ampicillin` (3,078,562);
- the dataset with the most samples in the release's `counts`: `salmonella_enterica__ampicillin` (1,997);
- the dataset with the most lasso folds whose release `fold_complexity` is 0 (null fits), ties to registered order: `mycobacterium_tuberculosis__clofazimine` (2).

Three species: *E. coli*, *S. enterica*, *M. tuberculosis*. The first two test the largest inputs; the third tests the null-fit path.

## Measures

For each dataset:

1. The inputs event: input fingerprint, samples, resistant, features, against the release and the lane's events.
2. Every fit event with `seconds` removed (scheme, fold, audit, trace of every path and refit point, model, error), from my unmodified binary, against the lane's pass 1 events: equal as parsed JSON or not, and the first difference if any.
3. My instrumented binary's fit events (without `seconds`) against my unmodified binary's.
4. From my instrumented binary's dumps: per fold, the partition (`train`, `fit`, `val`), labels and `features` against the lane's dumps; scores against the lane's dumps, bit for bit; the pooled out-of-fold scores (null folds at training prevalence, as `pool_folds` in `src/learn/eval.rs` does) against the release's lasso `scores`, bit for bit, and calls at 0.5.

## Pass rule

The replication agrees with the lane when, for all three datasets, measures 1 to 4 are equal with no exception. Any difference is reported with the dataset, fold and field. Timings are recorded but not compared (8 threads here, 16 there).
