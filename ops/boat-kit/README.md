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

KIT11 lets the runner run on the Boat T3 controller as well as the laptop. `GRM_HOST=linux` turns this mode on, and `GRM_BOAT` must name the boat binary. Staging dirs go under `GRM_STAGE_DIR`, with no `wslpath`. memcheck reads the controller's MemAvailable and its cgroup room, and each must be `GRM_MIN_AVAIL_MB` (512) or more. Without `GRM_HOST`, the runner works as before on the laptop. KIT12 counts each cgroup's reclaimable page cache as room: `memory.stat`'s `inactive_file` less `file_dirty` and `file_writeback`, never below 0. The work program's rule is not to treat raw `memory.current` as non-reclaimable RAM.

`stage-local` takes a push-data set from a dir that is already on the box, such as a data snapshot's copy of H-1's data root, and only the manifest crosses. drive.sh uses it for an H-1 run box when the run file sets `DATA_SNAPSHOT_DIR`, `VOCAB_SNAPSHOT_DIR`, `VOCAB_MANIFEST` and a `grm-tmpl-*` `TEMPLATE`.

A box made from H-1's data snapshot starts with its own 76 GB copy of the master root, which leaves about 202 GB free on an xlarge. The pilot and the Salmonella and TB main parts need more than that. With `DATA_SNAPSHOT_PRUNE=1`, each run box runs `h1-box.sh prune` once its FASTA set is staged. The step first checks every set against its pin, and that the metadata set was linked from that folder. It then deletes the box's copy. The staged sets are hard links, so they keep their files, and the run gets back the disk only the copy held. The snapshot itself isn't touched. A refused prune stops the stage with 81 before verify.

## k add-on boxes

KIT12 adds a stage kind, `addon`, for H-1's k21 and k25 arms. Each runs as a part on its own box with `run.py --reference-from`, against the main part that left the k arms out. Column 10 of its stage row names that one main part, which must be a run box of the same table. An add-on box is a run box in every step. It also waits for the main part's DONE and pulled WORK, and its H-1 flags must be the main part's.

drive.sh pins the main part's records as a one-part combine set and sends them through push-data to `/home/user/grm-combine` after the vocabulary. stage.env names the part (`H1_REFERENCE`) and the set's pin (`H1_REFERENCE_SHA256`). On the box, `h1-box.sh` checks that set against its pin before union, verify and the run, and passes `--reference-from` to each. `r11-box.sh` runs an add-on's stage as `h1-box.sh run`. The combine box then takes the add-on box as a part beside its main part.

## Two slow hosts in a row

drive.sh refuses a new box with 78 while the budget ledger's last two rows are slow hosts: "the last two boxes had slow hosts; stop and report". That stays the default.

KIT13 adds `SLOW_PAIR_PASS`. With `SLOW_PAIR_PASS=1`, a stage whose row is SLOW_OK (column 6) passes the stop, and NOTES.md gets a line saying so. A stage that isn't SLOW_OK still stops, and the budget cap still refuses. The ledger still records the new box's host. Unset or 0 keeps the stop, and any other value is refused with 2 before any box. A restart script sets it in the environment, and a run file that sets it wins, since drive.sh sources the run file.

H-1 needed this on 11 Nov 2026. Both S. aureus add-on boxes drew slow hosts, so its combine would have been refused, and with the other boxes closed no fast row could come. Every add-on and combine row was SLOW_OK, and H-1 had no xlarge-hour cap. The run's owner passed the stop with a copy of drive.sh that turned the stop into a note, run through a restart script. `SLOW_PAIR_PASS=1` gives that restart the kit's own drive.sh, for SLOW_OK rows only.

## Host class: what the speed test measures

`r11-box.sh host` times the source's release lib tests (`cargo test --release --lib`, built before the clock starts) and calls the host slow above `SLOW_LIB_SECS` (60 s). H-1's 28 host steps of 9 to 11 Nov 2026, from the run's drive notes:

| Class | Boxes | Lib tests |
|---|---:|---|
| fast | 22 | 5.5 to 54.1 s, median 26.4 s |
| slow | 6 | 77.6 to 936.8 s, median 168.5 s |

One CPU model spreads across the line. AMD Ryzen 9 9950X hosts took 5.5 to 93.9 s over 16 boxes, one of them slow. AMD EPYC 4585PX hosts took 27.0 to 242.3 s over 7 boxes, four of them slow.

N. gonorrhoeae's main part ran on two EPYC 4585PX hosts. The test took 54.1 s on 9 Nov (fast) and 242.3 s on 11 Nov (slow), 4.48 times as long. Over the 82 arm and dataset pairs both boxes wrote, the 11 Nov box took a median 1.29 times as long (0.84 to 6.89). The median by arm runs from 1.02 (k-h) to 2.98 (k-p), and the default arm's is 1.45. The 11 Nov box finished the stage in 1 h 24 min, between the stage table's fast estimate (1.19 h) and its slow estimate (2.04 h).

So the test over-reads. It put that host at 4.48 times its twin while the arms ran 1.29 times as long. On a stage that isn't SLOW_OK, a slow class closes the box with 80: M. tuberculosis's main part drew a host the test timed at 226.9 s, was closed and ran again on the next box. [drafts/host-class-short-arm.md](drafts/host-class-short-arm.md) is the design that classes a host by a short arm instead, with the measurement and tests it needs. It isn't built.

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

- `test-drive-h1-kit13.sh` checks the two-slow-hosts stop and `SLOW_PAIR_PASS` on the stub runner (13 checks), and `mutants-kit13.py` plants six faults in drive.sh for it to catch. Both take their paths as arguments: `test-drive-h1-kit13.sh WORK DRIVE_SH KIT_DIR` and `mutants-kit13.py DRIVE_SH KIT_DIR TEST_SH WORK OUT_JSON`.

The tests run from the harness's paths, which `make-harness.sh` builds. `test-restore-v7.sh.next` and `test-side-re.sh.next` read real restore sets, so round 2's chain ran neither. The other tests use synthetic files and copies from the harness.

GRM's CI doesn't read `ops/boat-kit/`.

## Drafts

`drafts/stages-h1.tsv` and `drafts/run-h1.conf` are drafts of H-1's stage table and run file. Both hold `@@NAME@@` placeholders, and `drive.sh` refuses a run file or stage row with one before it makes any box. The install leaves them out.

`drafts/host-class-short-arm.md` is a design, not code: a host classed by a short arm, the measurement that would set its line, and the tests it needs.
