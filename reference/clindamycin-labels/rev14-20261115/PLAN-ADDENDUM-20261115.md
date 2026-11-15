# Lane clinda: plan addendum of 15 Nov 2026, fixed before any count or score is read

`PLAN.md` (sha256 `b14ec01aaef5471aa1b8f0e3cf53fd7a3a9485db0dca8a6f5c7d9d8cffd7915b`, recorded 2026-11-14T18:11:52Z) stands as written. Its hash was checked on this machine at 2026-11-15T11:31:55Z and matches. Its method, datasets, study field, choice rule, measures and caveat rule do not change. This addendum fixes only what the restart on another machine requires. It was written before any revision 14 label, known-gene call, PMID by dataset, count by class, metric or score of a clindamycin dataset was read on 15 Nov. Nothing survives of the first run's counts: `results/` was empty.

## A1. Commit and binary

The plan names `c435000`. This machine's read-only repository is at `852ebc9`. `git diff --name-only c435000 852ebc9 -- src/ Cargo.toml Cargo.lock` lists nothing; the only test change is a Python suite (`tests/givens_export/`). So the Rust source of the plan's commit is unchanged. Build from a scratch clone at `852ebc9`: `CARGO_HOME=/home/user/scratch/clinda/cargo-home cargo fetch --locked`, then `cargo build --release --locked --offline` with `RUSTUP_TOOLCHAIN=1.98.1`, `CARGO_TARGET_DIR=/home/user/scratch/clinda/target`, `CARGO_BUILD_JOBS=8`. Done at 11:33:29Z, rc 0. Binary sha256 `b20b5f1415a011293a63f9fc8b90ef54b3f0935952cbcfc67aa0a36c5bcc310e`. Every arm uses this one binary.

## A2. Threads

`RAYON_NUM_THREADS=8` when at most one of the lanes `rev108` and `studysplit` is building or fitting, and 5 when both are, read from `uptime` and `ps` before each job. The count of threads is recorded per job. Whether a refit depends on the thread count is read from the "with" arm's acceptance test (section 5 of the plan), since the release was fitted on another machine.

## A3. No network

Every `grm` stage runs under `unshare -rn`: the process sees only a downed loopback and no DNS (tested 11:35Z). Reason: annotation reads each genome's BV-BRC features from `bvbrc/features/<id>.json` and fetches on a miss, and the program's client sends a User-Agent that names the project; this lane sends no project name to a public service. A fetch therefore fails instead of going out. If a stage fails for that reason, the failure is recorded and the dataset is not fitted by another route without saying so.

## A4. The scratch data roots

- `root-with` (arm "with", and the shared build of count, structure and matrix): `cp -rs` of the public root (real directories, one symbolic link per file), except `curated/`, which is a real copy, and `species/`, which is built here. The stages write only under `species/<slug>/`.
- `root-without-<dataset>` (arm "without", and each placebo arm): the same links, `curated/` copied with `phenotypes.jsonl` filtered (only the clindamycin rows of the removed genomes go), and `species/<slug>/` a real directory holding one link per file of `root-with/species/<slug>/` except `results/`, which is empty and real.
- Proof that nothing else changed: the sha256 of every file under `root-with/species/<slug>/` other than `results/` before and after each learn; a `diff` of the two phenotype files; `cmp` of every other `curated/` file; and `find /home/user/grm-h1raw -newer <flag>` after the runs, which must list nothing (the public root is not written).

## A5. Script correction

`scripts/count_by_study.py` prefiltered raw records with a case-sensitive `"lindamycin" in line`. It now lower-cases the line first, so a raw name in capitals is not skipped. The exact test is unchanged: the raw name must be one that `curated/drug_names.json` resolves to the drug `clindamycin`. This changes no rule of the plan.

## A6. Times and places

The placebo launch limit "none after 21:30Z on this machine's clock" means 21:30Z on 15 Nov on this machine, which ends at 23:17Z. The map row the plan cites as line 317 is line 319 of today's `PROJECT-MAP.md`; its text is the same.

## A7. Points the plan leaves to the implementation, fixed now

- A split whose report gives an `unavailable_reason` is reported as unavailable and enters no rule.
- The bootstrap draws, per replicate, as many clusters as the full cohort has, with `random.Random(20261114)` over the sorted cluster ids; the single-sample bootstrap uses `random.Random(20261115)` over the cohort order. Within a replicate, V1, V2 and V3 are computed on the same draw.
- The known-gene rule's calls do not depend on training, so its training effect is 0 by construction; it is shown and enters the rule through its total effect.
- Effects use the lane's own "with" refit when it differs from the release in any score (plan section 5); the released figure is shown beside it.
