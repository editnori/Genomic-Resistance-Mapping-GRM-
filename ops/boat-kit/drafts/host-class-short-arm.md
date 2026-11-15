# Design: class a host by a short arm (KIT13 draft, 15 Nov 2026)

Not built. This is the change, the measurement that would set its line, and the tests it needs. The kit's README, "Host class: what the speed test measures", gives the record that asks for it: the lib-test clock put N. gonorrhoeae's 11 Nov host at 4.48 times its 9 Nov twin, and the arms ran a median 1.29 times as long.

## Why it isn't built yet

- **No line can be set from what is measured.** No Boat host has run a short arm beside the lib tests. The run's per-arm seconds sit in each arm's `runs.json`, in stage folders this kit doesn't keep, and the one summary at hand is N. gonorrhoeae's two hosts, by arm. Setting a fast or slow line for a short arm now would replace one unmeasured line with another.
- **It needs a box to test.** The box side runs the release `grm` on a synthetic root. The kit's tests run on stubs, and the lane that wrote this could make no box.

`SLOW_PAIR_PASS` (README, "Two slow hosts in a row") is built. It doesn't depend on how a host is classed.

## The change

### Settings

`HOST_CLASS_BY` in the run file is `lib` (the default) or `arm`. With `lib`, nothing changes: the host step's lib-test clock and `SLOW_LIB_SECS` decide, as now. With `arm`, the run file must also set `HOST_ARM_RATIO`, a number above 1 that has no default, and drive.sh refuses `arm` without it (2, before any box). drive.sh writes `HOST_CLASS_BY` and `HOST_ARM_RATIO` to `stage.env`, and `HOST_ARM_REF_SECS` too on every box but the reference box.

### The box: `r11-box.sh host`

The build, the binary hashes and the lib-test clock stay as they are, and `receipts/host.txt` keeps its `speed_lib_tests_secs` line. With `HOST_CLASS_BY=arm`, after the lib tests:

1. `python3 -I "$S/scripts/input-comparison/synthetic_root.py" "$W/hostarm/raw"` with its default size and seed. It needs only the standard library and no network.
2. One data root for the default arm, laid out as `run.py`'s `root` lays it out (the raw folders linked, line 539) and built with the commands it gives the default: `grm --data ROOT curate --cap 2000 --min-class 25 --selection legacy` (line 540), then `count --k 31`, `structure` and `matrix` (lines 566 to 568). The groups file (STRUCTURE) is built as `run.py`'s `pool` builds it (line 820), over the default's cohort alone. These steps are not timed.
3. The timed step: one `grm arm --arm default --blocks K --screened-features 4000 --repeats 3 --pairing same_plan --groups STRUCTURE` on the root's one dataset (`salmonella_enterica__meropenem`), as `run.py` calls it (lines 954 to 965), with `THREADS` as for the stage.
4. `receipts/host.txt` gains `short_arm_secs=S`, and the step prints `HOST_ARM_SECS S`. With `HOST_ARM_REF_SECS` set, the host is slow when S is above `HOST_ARM_REF_SECS` times `HOST_ARM_RATIO`, and `HOST fast|slow` says so. A failed step stops the host step with 3, as a failed build does, and prints no `HOST` line.

### The drive

- The reference box (`REF_BOX`, the box that fixes the binary) is classed by its lib-test line, as now, since there is nothing to compare it with. Its `HOST_ARM_SECS` becomes `RUN_DIR/HOST-ARM-REF-<source commit>`, written as `BINARY-REF-<source commit>.sha` is.
- Any other box with `arm` waits for that file (79, no box), as it waits for the reference binary, and passes its seconds as `HOST_ARM_REF_SECS`.
- The host's class decides what it decides now: 80 for a slow host on a stage that isn't SLOW_OK, the lifetime set once for a SLOW_OK stage, and the ledger's `host` column, which the two-slow-hosts stop reads. NOTES.md gets both clocks on the host line.
- A box output with `arm` and no `HOST_ARM_SECS` line closes the box with 81.

## The measurement that sets `HOST_ARM_RATIO`

Before any run uses `arm`, on at least six xlarge boxes from H-1's data template, record per box: the lib-test clock, the short arm (step 3) three times, and the default arm on one fixed N. gonorrhoeae dataset from the template's data, the same dataset on every box. The ratio's line goes where the short arm's ratio to the first box separates the boxes whose real-arm ratio is above the stage table's slow allowance (`EST_SLOW / EST_FAST`, 1.19 to 2.0 in H-1's table) from those below it. If no line separates them, the short arm doesn't predict the real arm either, and `lib` stays.

## The tests it needs

On the stub runner, which gains `SR_HOST_ARM_SECS` (the `HOST_ARM_SECS` line it prints at `--log host`), in a new `test-drive-h1-kit14.sh`:

- **H1** `HOST_CLASS_BY` unset: the class is the `HOST` line, and KIT10 to KIT13's drive tests give the same verdicts.
- **H2** `arm` without `HOST_ARM_RATIO`, or with one that isn't a number above 1: 2, no box.
- **H3** `arm` on the reference box: `HOST-ARM-REF-<commit>` holds its seconds, and its class is its `HOST` line.
- **H4** `arm` on another box before the reference file exists: 79, no box.
- **H5** a box whose seconds are at or under the reference times the ratio while its `HOST` line says slow: fast, so no 80, no extend and a `fast` ledger row.
- **H6** seconds above the reference times the ratio while `HOST` says fast: slow, so 80 and closed on a stage that isn't SLOW_OK, and the lifetime set once on a SLOW_OK stage.
- **H7** `arm` and no `HOST_ARM_SECS` line: 81, closed.
- **H8** `stage.env` carries `HOST_CLASS_BY` and `HOST_ARM_RATIO`, and `HOST_ARM_REF_SECS` on every box but the reference box.
- **H9** a reference file of another source commit is not read.

On the box script, with a stub `grm` whose `arm` sleeps 2 s and whose other commands return at once, in a new `test-h1-box-kit14.sh`:

- **B1** `lib` runs no synthetic root and no `grm arm`.
- **B2** `arm` runs the root, curate, count, structure, matrix and arm in that order, and `short_arm_secs` is at least 2 and below the time of all six steps together.
- **B3** any of the six steps failing stops the host step with 3 and prints no `HOST` line.
- **B4** `receipts/host.txt` holds both `speed_lib_tests_secs` and `short_arm_secs`.

Planted faults the tests must catch: the class taken from the lib line under `arm`; the ratio test reversed; the reference box compared with itself; a reference file read across source commits; the timed span taking in the untimed steps.

One check needs a real box: the short arm with the release `grm` on a Boat host, recorded with the lib tests. That is the first box of the measurement above.
