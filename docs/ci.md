# Continuous integration

`.github/workflows/ci.yml` runs the Rust, Python and deployment-toolkit checks on every push and pull request. The jobs use synthetic and unit fixtures from the tree: none reads real data, uses a secret or deploys anything, and the network serves only the source, the toolchain and packages. The workflow can only read the repository (`permissions: contents: read`), its checkouts keep no credentials, and a newer push to a pull request cancels that pull request's older run.

## Jobs

Each job runs on GitHub's hosted `ubuntu-24.04` runner. The counts are those at `b32456e5` (evaluation revision 11).

| Job | Runs | At `b32456e5` |
|---|---|---|
| `rust-lint` | `cargo fmt --check`, then `cargo clippy --locked --all-targets -- -D warnings` | pass |
| `rust-test` | `cargo test --release --locked` | 385 passed, 0 failed, 8 ignored |
| `rust-doc` | `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps --locked` | fails, 9 errors (see [Known failures](#known-failures-at-b32456e5)) |
| `metric-audit` | `cargo build --locked --release` in `scripts/metric-audit/foldgen`, then `python3 -I -B -W error <file> -v` for each test file | `test_check_results.py` 106, `test_audit.py` 98, `test_scripts.py` 3; `scripts/lasso-scan/test_scan.py` 26 with 4 skipped |
| `python`, one job per suite | `python3 -s -P -B -W error -m unittest discover -s <suite> -t . -v` under `env -i` at `PYTHONHASHSEED` 0, 1 and 4242, then once with `-I` in place of `-s -P` | `tests/hyperopt` 257, `scripts/sampling-audit/tests` 34, `tests/ablation` 238, `tests/temporal` 163, `tests/multitask` 72 |
| `deployment` | `bash scripts/deployment/tests/dry-run.sh`, then `bash scripts/deployment/tests/apply.sh` | 195 passed; 488 passed, flaky (see [Known failures](#known-failures-at-b32456e5)) |

The `python` jobs discover `scripts/sampling-audit/tests` without `-t .` and set `GRM_TASK_J_STATS` to the tree's `scripts/ablation/stats.py`. The `tests/ablation` job also runs once without it, at `PYTHONHASHSEED=0`, so `test_schema.py` finds that file by itself. A unittest run fails if any test is skipped, except `test_scan.py`, whose 4 `RealBinary` tests need the binary and fixture root that `scripts/lasso-scan/smoke.sh` builds. A `python` job also fails if its suite changed a tracked file.

The deployment suites run offline, in temporary roots, with loopback ports only. They need a git checkout (`build-release.sh` archives `HEAD`), `flock`, `fuser` (psmisc), `curl`, `setsid`, `pgrep`, `python3` and `cc`; the job installs psmisc and curl and checks for each tool. `apply.sh` takes about 4 to 20 minutes, and the job's timeout is 45.

## Toolchain and runner

`RUSTUP_TOOLCHAIN` pins Rust 1.98.1 for every cargo command. Each Rust job installs it with rustup's minimal profile (`rust-lint` adds rustfmt and clippy) and checks `rustc --version`. It caches cargo's registry, its git checkouts and the job's own `target/`, keyed by the toolchain, `Cargo.lock` (foldgen's own for `metric-audit`) and `.cargo/config.toml`.

`.cargo/config.toml` builds for `target-cpu=x86-64-v3`, which needs AVX2. GitHub's hosted x86-64 runners have AVX2, so the build, its build scripts and the tests run there. Each Rust job checks `/proc/cpuinfo` for `avx2` first, so a runner without it fails with that message rather than an illegal instruction.

Python is the runner's system Python 3.12 with Ubuntu 24.04's `python3-numpy`, which is numpy 1.26.4. The `python` jobs check both versions after the install, so an image that changes either fails there. `-s` and `-I` keep user site-packages out, and `env -i` with `PATH=/usr/bin:/bin` keeps out the runner's environment too.

`actions/checkout` v7.0.1 is pinned at `3d3c42e5aac5ba805825da76410c181273ba90b1` and `actions/cache` v6.1.0 at `55cc8345863c7cc4c66a329aec7e433d2d1c52a9`. Both commits came from `git ls-remote` of the action's repository on 2026-11-01. The tags are lightweight, so each names its commit directly, and GitHub's commits API returned both. To update an action, take the new tag's commit the same way and keep the version comment beside it.

## Known failures at `b32456e5`

- **`rust-doc` fails.** rustdoc reports 9 errors under `-D warnings`, in doc comments older than revision 11.
- **`apply.sh` flakes.** Right after a test kills a start with SIGKILL, the next step sometimes finds `deploy.lock` still held and refuses ("another deployment holds …/deploy.lock"), and later assertions fail after it. It shows most on small runners, and GitHub's standard runner for a public repository has 4 vCPUs. It was seen once in 4 runs on 16 vCPU hosts (477 passed, 13 failed); the other three passed 488 of 488. When only this job fails, read its first `FAIL` line before re-running it.

## Not in CI

The web app's `npm --prefix web test`, `run typecheck` and `run build` are outside this workflow. Checks that need real data stay out too, because ordinary CI doesn't hold that data and must not be given it:

- the feature-block registry's catalog name check, which runs `proxy_kind` over every gene group in the pinned AMRFinderPlus catalog, so that no real name trips the deny-list and drops the P block;
- `scripts/check-results.py` and `scripts/metric-audit/audit.py` on saved reports in `data/`;
- the ignored test `learn::lasso::tests::mtb_delamanid_fold0_reaches_the_certified_optima`, which reads a local diagnosis folder;
- the deployment mutation check, `data/work/deployment-mutation-20261027/grm-mutate-v3.py`.

The other 7 ignored tests are opt-in runs: the full synthetic-assurance runs (`scripts/synthetic-assurance/run-full.sh`), the lasso scan smoke test (`scripts/lasso-scan/smoke.sh`), the matrix merge benchmark with its two case tests, and the open-file probe, which a default test starts in its own process.
