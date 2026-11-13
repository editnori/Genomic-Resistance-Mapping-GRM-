# Continuous integration

`.github/workflows/ci.yml` runs the Rust, Python, deployment-toolkit and web checks on every push and pull request. The jobs use synthetic and unit fixtures from the tree: none reads real data, uses a secret or deploys anything, and the network serves only the source, the toolchain and packages. The workflow can only read the repository (`permissions: contents: read`), its checkouts keep no credentials, and a newer push to a pull request cancels that pull request's older run.

## Jobs

Each job runs on GitHub's hosted `ubuntu-24.04` runner. The counts are this tree's, measured on Boat with the same commands (evaluation revision 14).

| Job | Runs | Result |
|---|---|---|
| `rust-lint` | `cargo fmt --check`, then `cargo clippy --locked --all-targets -- -D warnings` | pass |
| `rust-test` | `cargo test --release --locked`, then `bash scripts/ci/rust-test-names.sh --check` | 461 passed, 0 failed, 9 ignored; the 470 names match `tests/rust-test-names.txt` |
| `rust-doc` | `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --locked` | pass |
| `metric-audit` | `cargo build --locked --release` in `scripts/metric-audit/foldgen`, then `python3 -I -B -W error <file> -v` for each test file | `test_check_results.py` 122, `test_audit.py` 111, `test_scripts.py` 3; `scripts/lasso-scan/test_scan.py` 26 with 4 skipped |
| `python`, one job per suite | `python3 -s -P -B -W error -m unittest discover -s <suite> -t . -v` under `env -i` at `PYTHONHASHSEED` 0, 1 and 4242, then once with `-I` in place of `-s -P` | `tests/hyperopt` 351, `scripts/sampling-audit/tests` 34, `tests/ablation` 280, `tests/temporal` 163, `tests/multitask` 72, `tests/input_comparison` 198 with 31 skipped, `tests/outside_study` 224, `tests/outside_study_results` 46, `tests/release_summary` 43 |
| `deployment` | `bash scripts/deployment/tests/dry-run.sh`, then `bash scripts/deployment/tests/apply.sh` | 285 passed; 564 passed |
| `web` | `npm --prefix web ci`, then `run typecheck`, `run build` and `test` (vitest), on Node 22.23.4; then `node ops/cloudflare/tests/previous-assets.mjs` and `long-key.mjs` | 43 test files passed, 1 skipped; 396 tests passed, 2 skipped; Worker tests 22/22 and 5/5 |

The `python` jobs discover `scripts/sampling-audit/tests` without `-t .` and set `GRM_TASK_J_STATS` to the tree's `scripts/ablation/stats.py`. The `tests/ablation` job also runs once without it, at `PYTHONHASHSEED=0`, so `test_schema.py` finds that file by itself. A unittest run fails if any test is skipped, except `test_scan.py`, whose 4 `RealBinary` tests need the binary and fixture root that `scripts/lasso-scan/smoke.sh` builds, and `tests/input_comparison`, whose 31 end-to-end tests need a release `grm` (`GRM_INPUT_E2E_BINARY`). That job requires exactly `OK (skipped=31)` (the matrix's `ok`), so any other skip still fails it. The end-to-end tests take about 9 minutes on a release build (Boat xlarge, 7 Nov 2026) and run on Boat. A `python` job also fails if its suite changed a tracked file. The `tests/outside_study` suite builds small git repositories in temporary folders to test the outside-study kit's frozen-file gate, so it needs `git`; its mutation check, `tests/outside_study/mutants.py`, is not a test module and runs by hand (docs/outside-study-test.md). The `tests/outside_study_results` suite scores nothing. It reads the committed output files under `reference/outside-study` and fails if the tables in `docs/outside-study-results.md` or `docs/outside-study-results-by-species.md` differ from them, if a sentence of those pages that carries a number differs from them, if a number in a page's text stands in no sentence the suite rebuilds, if a file differs from `SHA256SUMS`, or if a run's record does not name the frozen kit, population and genomes. Two of these suites also hold the files the web app bundles to the reference files. `tests/outside_study_results` fails if `web/src/data/outside-study.json` differs from what `scripts/outside-study-results/site_data.py` writes from `reference/outside-study`. `tests/input_comparison` fails if `web/src/data/input-comparison.json` differs from what `scripts/input-comparison/site_data.py` writes from the two readings under `reference/input-comparison`, or if a reading differs from the sha256 its README names. The `tests/release_summary` suite builds a small synthetic release in a temporary folder and runs the scripts of `scripts/release-summary` on it. It also reads the committed outputs under `reference/release-summary` and fails if their counts do not add up or if a file differs from the sha256 its README lists. It reads no published release, and it fails if a table of `docs/release-summary.md` differs from the committed outputs.

The deployment suites run offline, in temporary roots, with loopback ports only. They need a git checkout (`build-release.sh` archives `HEAD`), `flock`, `fuser` (psmisc), `curl`, `setsid`, `pgrep`, `python3` and `cc`; the job installs psmisc and curl and checks for each tool. `apply.sh` takes about 4 to 20 minutes, and the job's timeout is 45.

## Rust test names

`cargo test` still passes when a module's tests leave the build. Put a `#[cfg(...)]` that is never true on a module's `mod tests`, and its tests are no longer compiled, so nothing runs them and nothing fails: doing that to `src/ablate/plan.rs`'s tests drops 21 tests.

The `rust-test` job therefore compares the tests it built with `tests/rust-test-names.txt`. That file lists every test, sorted, as its crate root and name, for example `src/lib.rs ablate::encode::tests::…` or `tests/foldmask_extended.rs …`. `scripts/ci/rust-test-names.sh --check` reads the names from `cargo test --release --locked -- --list` and fails on any difference, printed as a diff: a `-` line is a test that left the build, and a `+` line is a new one.

- **Ignored tests** are listed too, each ending in ` [ignored]`, from a second listing with `--ignored`. Ignoring a test or un-ignoring one changes the list.
- **Doc tests** are listed under `doc`, without the `(line N)` that rustdoc adds, so editing lines above a doc test leaves the list alone. Two doc tests on one item give two equal lines, and the list keeps both.
- **Cargo's output:** the script stops if it reads a different number of tests than libtest counts, or a test before cargo names its target. A change in that output then fails the job instead of shortening the list.

When a change adds, removes, renames, ignores or un-ignores a test on purpose, run `bash scripts/ci/rust-test-names.sh --write` after making it, and commit the rewritten file with the change, so the review shows which tests came and went. Run it on Linux x86-64 with Rust 1.98.1 (`RUSTUP_TOOLCHAIN=1.98.1`), as the job does. With no argument, the script prints the list. The file pins names only: a test that still exists but asserts less still matches.

## Web app

The `web` job runs the scripts in `web/package.json`: `npm ci` from `web/package-lock.json`, then `typecheck` (`tsc -b`), `build` (`tsc -b && vite build`) and `test` (`vitest run`). It fails if they changed a tracked file. `src/charts/realdata.test.ts` skips, because it reads a published report from `GRM_REPORT_DIR`, which CI doesn't set; that file and its 2 tests are the skipped ones. The job then runs the Worker's two plain-Node tests, `ops/cloudflare/tests/previous-assets.mjs` and `long-key.mjs`, which import `ops/cloudflare/worker/src/index.js` against a stub bucket and need no package. `release.sh`, `routes.py` and `ranges.py` need Miniflare and `grm serve`, so they stay outside CI (docs/static-hosting.md, Checks).

`NODE_VERSION` pins Node 22.23.4, and the job checks `node --version` after the install. The README asks for Node.js 22, vitest 5 needs 22.12 or newer, and 22.23.4 is the newest 22.x release (2026-10-20). `actions/setup-node` v7.1.0 installs it, pinned at `949feb2413d6458794dcd2491c4babbbce0c15c1`, and caches npm's download cache, keyed by the lock file. The commit came from `git ls-remote` of `actions/setup-node` on 2026-11-02, the same way as the other two actions below.

## Toolchain and runner

`RUSTUP_TOOLCHAIN` pins Rust 1.98.1 for every cargo command. Each Rust job installs it with rustup's minimal profile (`rust-lint` adds rustfmt and clippy) and checks `rustc --version`. It caches cargo's registry, its git checkouts and the job's own `target/`, keyed by the toolchain, `Cargo.lock` (foldgen's own for `metric-audit`) and `.cargo/config.toml`.

`.cargo/config.toml` builds for `target-cpu=x86-64-v3`, which needs AVX2. GitHub's hosted x86-64 runners have AVX2, so the build, its build scripts and the tests run there. Each Rust job checks `/proc/cpuinfo` for `avx2` first, so a runner without it fails with that message rather than an illegal instruction.

Python is the runner's system Python 3.12 with Ubuntu 24.04's `python3-numpy`, which is numpy 1.26.4. The `python` jobs check both versions after the install, so an image that changes either fails there. `-s` and `-I` keep user site-packages out, and `env -i` with `PATH=/usr/bin:/bin` keeps out the runner's environment too.

`actions/checkout` v7.0.1 is pinned at `3d3c42e5aac5ba805825da76410c181273ba90b1` and `actions/cache` v6.1.0 at `55cc8345863c7cc4c66a329aec7e433d2d1c52a9`. Both commits came from `git ls-remote` of the action's repository on 2026-11-01. The tags are lightweight, so each names its commit directly, and GitHub's commits API returned both. To update an action, take the new tag's commit the same way and keep the version comment beside it.

## Not in CI

The web app's browser check, `web/scripts/browser-check.cjs`, is outside this workflow: it needs Chromium and a report tree. So are the input comparison's end-to-end run (`tests/input_comparison` with `GRM_INPUT_E2E_BINARY`, a release build and about 9 minutes) and the Worker's `release.sh`, `routes.py` and `ranges.py` (Miniflare and `grm serve`); they run on Boat. Checks that need real data stay out too, because ordinary CI doesn't hold that data and must not be given it:

- the feature-block registry's catalog name check, which runs `proxy_kind` over every gene group in the pinned AMRFinderPlus catalog, so that no real name trips the deny-list and drops the P block: the ignored test `ablate::blocks::tests::the_real_catalog_s_gene_groups_pass_the_deny_list`, run with `GRM_CATALOG_DIR` set to a `grm fetch catalog` folder;
- `scripts/check-results.py` and `scripts/metric-audit/audit.py` on saved reports in `data/`;
- the web suite's `src/charts/realdata.test.ts`, which reads a published report from `GRM_REPORT_DIR` and skips without it;
- the ignored test `learn::lasso::tests::mtb_delamanid_fold0_reaches_the_certified_optima`, which reads a local diagnosis folder;
- the deployment mutation check, `data/work/deployment-mutation-20261027/grm-mutate-v3.py`.

The other 7 ignored tests are opt-in runs: the full synthetic-assurance runs (`scripts/synthetic-assurance/run-full.sh`), the lasso scan smoke test (`scripts/lasso-scan/smoke.sh`), the matrix merge benchmark with its two case tests, and the open-file probe, which a default test starts in its own process.
