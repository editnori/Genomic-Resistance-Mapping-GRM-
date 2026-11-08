# Boat box kit

GRM's lanes run their Boat boxes from the laptop with this kit. It has three parts:

- `runner/grm-boat.sh` runs one lane's box through the Boat CLI. Its commands include new, push, push-data, stage-local, run, extend, pull and close, which saves the box's usage and deletes it. Each step adds a line to the lane's NOTES.md.
- `drive/drive.sh` runs one stage box from start to delete through the runner, with a run file and a stage table. It keeps the budget ledger and the deadline.
- `kit/` holds the scripts the runner pushes to each box.
  - `r11-box.sh` runs a revision-N stage, and it hands an H-1 stage to `h1-box.sh`.
  - `mtb-box.sh` runs the MTB cleaned 4k run.
  - `stage-manifest.py` and `h1-manifest.py` write the manifests that pick each stage's input files.
  - `chain-lib.sh` runs a chain of steps and exits 1 when any step fails.

`ops/boat-kit/` holds KIT10's set, with KIT11's changes to the runner and drive and KIT12's k add-on boxes. The other kit files stay on the laptop, and the kit's older history stays in the laptop's kit repo at `f764602`.

## Laptop paths

The laptop runs the kit from these paths. `install/install-kit10.sh` puts this set there.

| Here | On the laptop |
|---|---|
| `runner/grm-boat.sh` | `/home/layth/review/grm-boat-runner-20261031/grm-boat.sh` |
| `drive/drive.sh` | `/home/layth/review/grm-port-r11-20261105/drive/drive.sh` |
| `kit/*` | `/home/layth/review/grm-port-r11-20261105/kit/` |

## Running from the T3 controller

KIT11 lets the runner run on the Boat T3 controller as well as the laptop. `GRM_HOST=linux` turns this mode on, and `GRM_BOAT` must name the boat binary. Staging dirs go under `GRM_STAGE_DIR`, with no `wslpath`. memcheck reads the controller's MemAvailable and its cgroup room, and each must be `GRM_MIN_AVAIL_MB` (512) or more. Without `GRM_HOST`, the runner works as before on the laptop.

`stage-local` takes a push-data set from a dir that is already on the box, such as a data snapshot's copy of H-1's data root, and only the manifest crosses. drive.sh uses it for an H-1 run box when the run file sets `DATA_SNAPSHOT_DIR`, `VOCAB_SNAPSHOT_DIR`, `VOCAB_MANIFEST` and a `grm-tmpl-*` `TEMPLATE`.

## k add-on boxes

KIT12 adds a stage kind, `addon`, for H-1's k21 and k25 arms. Each runs as a part on its own box with `run.py --reference-from`, against the main part that left the k arms out. Column 10 of its stage row names that one main part, which must be a run box of the same table. An add-on box is a run box in every step. It also waits for the main part's DONE and pulled WORK, and its H-1 flags must be the main part's.

drive.sh pins the main part's records as a one-part combine set and sends them through push-data to `/home/user/grm-combine` after the vocabulary. stage.env names the part (`H1_REFERENCE`) and the set's pin (`H1_REFERENCE_SHA256`). On the box, `h1-box.sh` checks that set against its pin before union, verify and the run, and passes `--reference-from` to each. `r11-box.sh` runs an add-on's stage as `h1-box.sh run`. The combine box then takes the add-on box as a part beside its main part.

## Installing

`install-kit10.sh` runs from KIT10's handback dir, which holds these files as `next/*.next` with its SHA256SUMS. Before it changes anything, it checks:

- the handback against SHA256SUMS
- each live file against the hash the tests ran on
- that no drive, runner or chain is running

`install/INSTALL.md` has each file's live and new hash, the steps and the rollback.

## Tests

`tests/` has the tests, and `tests/box/` has the test box's scripts:

- `make-harness.sh` packs the files into one tarball on the controller.
- `setup-mirror.sh` and `gen-port-master.sh` build the synthetic data on the box.
- `chain-tests.sh` runs every gate.
- `mutants.py` runs the mutants.
- `chain-kit11.sh` and `mutants-kit11.py` run KIT11's gates and mutants. They run in KIT10's harness with three additions: `cand/` holds KIT11's `grm-boat.sh` and `drive.sh`, `k10final/` holds KIT10's for the control runs, and the KIT11 tests (`test-runner-kit11.sh`, `test-drive-h1-kit11.sh`) sit in `tests/`.
- `chain-kit12.sh` and `mutants-kit12.py` run KIT12's gates and mutants in KIT11's harness. `cand/drive.sh` and `kit/` hold KIT12's files, `k11final/` and `k11kit/` hold KIT11's for the control runs, and the KIT12 tests (`test-h1-box-kit12.sh`, `test-drive-h1-kit12.sh`) sit in `tests/`.

The tests run from the harness's paths, which `make-harness.sh` builds. `test-restore-v7.sh.next` and `test-side-re.sh.next` read real restore sets, so round 2's chain ran neither. The other tests use synthetic files and copies from the harness.

GRM's CI doesn't read `ops/boat-kit/`.

## Drafts

`drafts/stages-h1.tsv` and `drafts/run-h1.conf` are drafts of H-1's stage table and run file. Both hold `@@NAME@@` placeholders, and `drive.sh` refuses a run file or stage row with one before it makes any box. The install leaves them out.
