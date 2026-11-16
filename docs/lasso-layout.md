# Lasso speed and crate layout

The lasso's speed used to depend on how the compiler split the crate into codegen units, and on code far from the lasso, with no change to the lasso's source and none to any output. This page gives the difference in the machine code, the change that removes it, and the command that checks it.

## Cause

`Wls::violators` (`src/learn/lasso.rs`) sums, for each zero column, the residual `q` over the column's rows. In a slow lasso this one loop held nearly all of the lasso's time. Where the loop kept its running sum was the register allocator's choice for the whole function the loop ended up in, and that choice changed with code elsewhere in the crate:

- fast: the sum stays in a register (`vaddsd (%r12,%rdi,8),%xmm0,%xmm0`);
- slow: every element loads the sum from the stack, adds, and stores it back (`vmovapd 0x40(%rsp),%xmm0`, `vaddsd`, `vmovapd %xmm0,0x40(%rsp)`), so each addition waits on the store before it.

Two commits with the same `lasso.rs` gave one each. Marking `violators` `#[inline(never)]` did not help: out of line, it still kept the sum on the stack with the release profile, with one codegen unit and with 16 codegen units. The likely reason is that the finished sum is held across `out.push`, which can call the allocator, and no floating-point register survives a call, so the sum must be stored around that call; the slow builds stored it on every element instead.

The additions are the same in the same order either way, so the fits are equal bit for bit.

## Fix

The gather is its own function, `column_sum`, marked `#[inline(never)]`. It makes no call but the bounds-check panic, which never returns, so nothing forces its total out of a register: with the release profile, with one codegen unit and with 16 codegen units it compiles to the same loop of seven instructions with the total in `xmm0`. `violators` calls it once per zero column. The test `column_sum_adds_in_row_order` fails if the order of the additions changes.

On `salmonella_enterica__tetracycline` at 5 threads (AMD EPYC 9455, Rust 1.98.1, November 2026), the lasso's median was 30.3 to 32.2 s in the four slow layouts and 8.3 to 8.7 s in the same layouts with `column_sum`, with every report equal outside its `seconds`.

`kkt_residual` has the same gather but writes each sum straight to `grad[k]`, with no call in between, and kept its total in a register in every build checked.

## Check

`scripts/lasso-layout/check.sh` builds `grm` from the working tree with the release profile, with one codegen unit and with 16 codegen units, and checks in each build that `column_sum` is a function of its own and that none of its instructions moves a floating-point register to or from the stack. It needs `nm` and `objdump` and exits 1 when a build fails. CI runs it on every push and pull request, in the `lasso-layout` job (docs/ci.md); the timing below is run by hand. To run it locally:

```bash
scripts/lasso-layout/check.sh
```

To time the lasso too, give it a data root where a species is counted and built:

```bash
scripts/lasso-layout/check.sh DATA_ROOT SPECIES DATASET [ROUNDS]
```

It then learns DATASET with each build in turn, ROUNDS times (default 3) at `RAYON_NUM_THREADS` (default 5), and prints the lasso's seconds of every run and each build's median; the three medians should be close. Its runs write to a new directory of links, so DATA_ROOT's reports are not touched. Timing is no CI gate, since runners are noisy: run it on a quiet machine.

No synthetic data of the repository runs a lasso long enough to time. The fixture root of `scripts/lasso-scan/smoke.sh` has 150 genomes and 5 feature rows and no FASTA, so `grm learn` does not read it, and the saved panels of `tests/data/lasso-stopping-fixtures.txt` have 60 to 508 samples and 2 to 128 features. So the timing takes a data root. The machine-code check needs none.
