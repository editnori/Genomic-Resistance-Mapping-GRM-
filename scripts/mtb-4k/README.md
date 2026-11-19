# The M. tuberculosis 4,000-genome run of 6 Nov 2026

The scripts here ran the M. tuberculosis run published as release
`data-mtb4k-r12-20261106`. They reached `main` unchanged from branch
`mtb-r12/merge`, and `4k-execution-plan.md` (revision 7) is the plan that run
followed. The matrix code it needs is on `main` too, behind
`grm matrix --disk-sparse` (`docs/matrix-scaling.md`, "Pipeline wiring").

## What the run was

- **Release.** `data-mtb4k-r12-20261106`, published 2026-11-06T19:30:37Z,
  with `grm-mtb-r12-4k-results.tar.gz` (54,833,067 bytes, sha256
  `128064f5656c9cd2ffc810ae31d70936a9546960d85aa7e5c6101c101f074913`) and its
  `.sha256` file (the release's asset list).
- **Evaluation revision 12.** `revision_n.py` holds it: `REVISION_N` 12, the
  checker commit `1767b862`, and the sha256 of that commit's four checker
  files.
- **Cohorts.** Four, all M. tuberculosis (`cleaned_rerun.py`): prefix, 2,000
  genomes; added, 2,000; full, 4,000 (the cap-4000 selection); matched2k,
  2,000. The species check fails 2 genomes, `1773.41701` and `1773.40949`,
  and no cohort holds them.
- **Four stages.** Boat boxes X1 to X4 (plan, section 10). X1 prepares the
  roots, builds the four matrices with
  `grm matrix --min-support 2 --disk-sparse --block-width 128`, screens
  duplicates and builds the structures. X2 trains prefix and added, X3 trains
  full, and X4 trains matched2k, compares the cohorts and checks the result
  against revision 12's trained 2,000-genome MTB root (section 7b).
- **Source commit.** No file on the branch, and nothing in the release's asset
  list, names it. The plan's section 3 checks the frozen source clone against
  `SRC_COMMIT` in `run-mtb.conf`, and `ops/boat-kit/drive/drive.sh` writes it
  to each stage's `src-commit.txt`; neither file is in this repository. The
  plan names the source as the merge of `main` at `1767b862` with this branch
  at `bc7b742`, which is `40c0d9f`; the branch's last commit is `f46a622`.

## Running it again from `main`

`ops/boat-kit` holds the box script (`kit/mtb-box.sh`) and the driver
(`drive/drive.sh`, `KIND=mtb`). Fill these first, each in a reviewed commit
or file of its own:

1. `revision_n.py`, `VALUES`, for the revision that runs. Revision 12's
   values stay as they are. `main` is at revision 14, and its
   `scripts/check-results.py` and `scripts/metric-audit/audit.py` differ from
   the revision-12 sha256, so profiles `rN` and `rN-cleaned` refuse on `main`
   until revision 14's checker commit fills them.
2. `cleaned_rerun.py`, `VALUES`: the species check, the selections, the input
   root and its file lists, the four cohort lists and the reference root.
   Revision 12's are the 6 Nov run's own. A revision-14 rerun takes its values
   from its own first run; none can be copied from this one.
3. A run file `run-mtb.conf` and its stage table `stages-mtb.tsv`, neither in
   this repository, with at least what `drive.sh` requires: `RUNNER`,
   `NOTES`, `RUN_DIR`, `STAGES`, `KIND=mtb`, `SRC`, `SRC_COMMIT`, `KIT`,
   `KIT_REPO`, `MASTER_MANIFEST`, `MASTER_SHA256`, `LANE_PREFIX`,
   `RAM_NEED_GIB` and `DATA_ROOT`.
4. The plan's section 1 gates, read for the new revision.

Then, where `drive.sh` runs, with `SRC` a clean worktree of `main` at
`SRC_COMMIT`:

    python3 -I scripts/mtb-4k/revision_n.py status
    python3 -I scripts/mtb-4k/cleaned_rerun.py status
    python3 -I scripts/mtb-4k/revision_n.py verify-checkout "$SRC"
    RUN_CONF=run-mtb.conf bash ops/boat-kit/drive/drive.sh X1
    RUN_CONF=run-mtb.conf bash ops/boat-kit/drive/drive.sh X2   # X2 and X3 may run together
    RUN_CONF=run-mtb.conf bash ops/boat-kit/drive/drive.sh X3
    RUN_CONF=run-mtb.conf bash ops/boat-kit/drive/drive.sh X4

## Tests

`tests/` holds 20 suites (`test_*.py`) and two helpers, `fixtures.py` and
`profile_boot.py`. The suites build synthetic fixtures and read none of the
run's genomes or results. Each runs under one evaluation profile of
`report_contract.py`, picked by `MTB4K_TEST_PROFILE`: `r8` (the default),
`rN` or `rN-cleaned`. The suites that call GRM take the binary from `GRM`,
and the metric suites take `FOLDGEN`:

    cargo build --release --locked
    (cd scripts/metric-audit/foldgen && cargo build --release --locked)
    MTB4K_TEST_PROFILE=rN GRM=target/release/grm \
      FOLDGEN=scripts/metric-audit/foldgen/target/release/foldgen \
      python3 -I -B scripts/mtb-4k/tests/test_validate_matrix.py

`test_k_review.py` and `test_k_review_r3.py` also take `--receipt OUT.json`.
One test needs run data: `test_added_gate_counts_derive_from_frozen_evidence`
in `test_k_review.py` reads the frozen MTB evidence file of the run's
revision-7 candidate and skips where that file is absent.

Every profile pins the sha256 of `scripts/check-results.py` and
`scripts/metric-audit/audit.py`. On `main` neither file holds the `r8` pin
or revision 12's, so each suite that loads `report_contract.py` stops at that
check until step 1 above fills the values for the revision that runs. Under
`rN-cleaned`, three suites (`test_evaluation_profiles.py`,
`test_prepare_execution_roots.py`, `test_report_contract.py`) also fail on
the branch's last commit: their tests expect the `rN` matched2k list.

`test_validate_matrix.py` builds its matrix with `grm matrix --disk-sparse`,
as X1 did. Only a disk-sparse count records the sketches `validate_matrix.py`
checks; a stock `grm count` records none, as on `main`.
