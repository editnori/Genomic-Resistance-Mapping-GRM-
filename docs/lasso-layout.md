# Lasso speed and crate layout

The lasso's speed used to depend on how the compiler split the crate into codegen units, and on code far from the lasso, with no change to the lasso's source and none to any output. This page says what differed in the machine code, what holds it now, and how to check it.

## What differed

`Wls::violators` (`src/learn/lasso.rs`) sums, for each zero column, the residual `q` over the column's rows. In a slow lasso this one loop held nearly all of the lasso's time. Where the loop kept its running sum was the register allocator's choice for the whole function the loop ended up in, and that choice changed with code elsewhere in the crate:

- fast: the sum stays in a register (`vaddsd (%r12,%rdi,8),%xmm0,%xmm0`);
- slow: every element loads the sum from the stack, adds, and stores it back (`vmovapd 0x40(%rsp),%xmm0`, `vaddsd`, `vmovapd %xmm0,0x40(%rsp)`), so each addition waits on the store before it.

Two commits with the same `lasso.rs` gave one each. Marking `violators` `#[inline(never)]` did not help: out of line, it still kept the sum on the stack with the release profile, one codegen unit and 16. The likely reason is that the finished sum is held across `out.push`, which can call the allocator, and no floating-point register survives a call, so the sum must be stored around that call; the slow builds stored it on every element instead.

The additions are the same in the same order either way, so the fits are equal bit for bit.

## What holds it

The gather is its own function, `column_sum`, marked `#[inline(never)]`. It makes no call but the bounds-check panic, which never returns, so nothing forces its total out of a register: with the release profile, one codegen unit and 16 it compiles to the same seven-instruction loop with the total in `xmm0`. `violators` calls it once per zero column.

`kkt_residual` has the same gather but writes each sum straight to `grad[k]`, with no call in between, and kept its total in a register in every build checked.

## Check

`scripts/lasso-layout/check.sh` builds `grm` from the working tree with the release profile, with one codegen unit and with 16, and checks in each build that `column_sum` is a function of its own and that none of its instructions moves a floating-point register to or from the stack. It needs `nm` and `objdump` and exits 1 when a build fails. Run it after a change to `lasso.rs`, to the release profile or to the toolchain:

```bash
scripts/lasso-layout/check.sh
```

To time the lasso too, give it a data root where a species is counted and built:

```bash
scripts/lasso-layout/check.sh DATA_ROOT SPECIES DATASET [ROUNDS]
```

It then learns DATASET with each build in turn, ROUNDS times (default 3) at `RAYON_NUM_THREADS` (default 5), and prints the lasso's seconds of every run and each build's median; the three medians should be close. Its runs write to a new folder of links, so DATA_ROOT's reports are not touched. Timing is no CI gate, since runners are noisy: run it on a quiet machine.
