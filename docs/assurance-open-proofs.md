# Open proofs in the synthetic assurance

The Bend2 fold-mask work proved laws about a Bend model of the fold masks, and checked the cases the Rust run sends to Bend against that model. Those are 35,272 cases, the 33,732 ext45 cases plus the 1,540 rand cases, out of the run's 2,669,501 fold-mask cases. `docs/synthetic-assurance.md` states two equivalences as trusted, not proved: the Rust code against the model, and the Bend source against its `bendtt` kernel term. The work-program ledger also names parsing and backend, which that doc does not describe.

Each one is listed below as a candidate proof target, in the order the sources give them, and none is chosen. This is work-program row A22: "narrow fold-mask laws and concrete conformance checks passed; broader implementation equivalence is unproven" (`docs/work-program-20261027.md:149`). Line numbers are those of `b90be1ca`, and a bare line number is a line of `docs/synthetic-assurance.md`.

## What was proved

- Line 3: the tests "call the public `grm` functions directly (`grm::learn::{fold_masks, training_folds, validate_folds}`, `grm::bits`, `grm::dna::KmerSpace`) and write every input and output as evidence for an independent oracle and for the offline Bend2 checks below."
- Line 65: "Bend and Lean are not cargo dependencies and nothing in `cargo test` runs them. The full fold-mask run writes two Bend books: `foldmask/bend/ext45/` (33,732 cases in 33 shards) and `foldmask/bend/rand/` (1,500 Bend random plus 40 out-of-contract cases, 1,540 in 13 shards). Each `cases_NNNN.bend` imports `../source/MODEL.bend` and `./CHECK.bend`, which are not in this repository."
- Line 77: "run each step for `PROOF.bend` and for `bendcheck/{ext45,rand}/VECTORS.bend`". The steps are `bend FILE`, `bend FILE -o FILE.bendtt`, `bend FILE --verdict` and `bendtt FILE.bendtt` (lines 80-83).
- Line 86: "Each must exit 0 and report `ALL PROOFS CHECK`. On 2026-10-27 all 12 steps passed".
- Line 93: "The universal laws in `PROOF.bend` are about `MODEL.bend`, not about the Rust code".

Lines 67-75 pin the inputs by sha256: `MODEL.bend`, `LAWS.bend`, `PROOF.bend`, `CHECK.bend`, Bend 2.0.39, Lean 4.34.0 and the `bendtt` kernel. The ledger sums this up as "Narrow model law and concrete Rust checks passed" (`docs/work-program-20261027.md:264`).

## Candidate proof targets

1. **Rust code and `MODEL.bend`.** Line 93: "Only 35,272 of the 2,669,501 fold-mask cases reach Bend. The universal laws in `PROOF.bend` are about `MODEL.bend`, not about the Rust code; the Rust-to-model abstraction (fold IDs at or above k become `Outside`) and the Bend-to-`bendtt` translation are trusted, not proved." The target is that the Rust code the cases come from agrees with `MODEL.bend` under that abstraction for every input, beyond the 35,272 cases that reach Bend.
2. **Bend to `bendtt` (translation).** The same sentence of line 93 states the Bend-to-`bendtt` translation as trusted, not proved. The steps emit the kernel term with `bend FILE -o FILE.bendtt` and check it with `bendtt FILE.bendtt` (lines 81 and 83). The target is that the emitted kernel term states what the Bend source states.
3. **Parsing.** Only the ledger names it: "parsing/translation/backend equivalence is not proved" (`docs/work-program-20261027.md:264`). `docs/synthetic-assurance.md` has no sentence about parsing and does not say which parse is meant.
4. **Backend.** Only the same ledger line names it. The doc names the tools (lines 73-75 and 80-83) but states no backend equivalence.

## Other limits the doc states

Lines 90 to 94 also list checks that are recorded but not independently derived. The doc does not call these equivalences.

- Line 90: "The specification the Rust tests check against is written in the same test files; it is not an independent implementation. No independent oracle has yet re-derived the emitted JSONL (random, validate, training, k-mer files)."
- Line 91: "`validate_tier` records each `validate_folds` outcome without asserting it."
- Line 92: "`canonical()`, `decode()` and the window extraction of `for_each_canonical` are recorded, not checked against an independent extractor."
- Line 94: "The recorded runs used one toolchain, platform and seed set."
