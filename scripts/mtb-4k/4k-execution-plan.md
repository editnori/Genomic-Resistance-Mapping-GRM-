# MTB matched 2k / added-2k / 4k execution plan — NOT AUTHORIZED

Revision 3 (2026-10-28, plan-and-validators lane). Revision 2 replaced the
wiring lane's plan (`plans/4k-execution-plan.md`, SHA-256 `8bb780c1…6fb8`).
The re-review of the fix (`176e3cc`, rebuilt trailer-free as `79158f9`, tree
`150500fe…`) returned `WIRING_REVIEW=CHANGES_REQUIRED` with plan and script
findings N1–N8 and no open Rust finding. Revision 2.1 closed N1–N8: the
annotation copy (N1), the pilot's kill and invalidation checks (N2, N8), the
matrix validator (N3), stale commit IDs (N5), the comparator's interlocks
(N6), and the frozen-input and lane checks (N7).

Revision 3 adopts one canonical validator set, merged from the Kimi package
(K, approved at the pre-execution gate) and this directory's own (W) as the
reconciliation (`RECONCILE.md`, `MERGE_REQUIRED_WITH_FIXES`) laid out: K's
strict report contract as `report_contract.py`, K's checks around the
repo's own revision-8 checker as `run_revision8_checks.py`, W's four-cohort
design in the merged `compare_cohorts.py`, and three new scripts:
`write_run_manifest.py`, `fingerprint.py` and `diff_revision7_prefix.py`
(sections 7 and 8). The code is the merge M0 of `79158f9` with the shared
branch at `3483065`, then the N1–N8 commit, then the validator commit; none
has been reviewed, and M0 is **not** the lasso-final merge (section 7).

**No command below has been run. No biological matrix, structure, training,
report or 8k work exists. No 8k stage is planned or permitted.**

## 1. Gates

`RUN_4K=YES` only when all of these hold; otherwise `RUN_4K=NO`.

1. **Wiring re-review.** A fresh independent reviewer inspects this commit
   and records `WIRING_REVIEW=APPROVED`.
2. **Pilot.** The parent authorizes the small frozen biological pilot
   (section 4), and every pilot acceptance check passes before the prefix
   stage.
3. **Validators reviewed.** Every script in this directory is
   independently reviewed: `prepare_execution_roots.py`, `pilot_interrupt.py`,
   `validate_matrix.py`, `screen_duplicates.py`, `report_contract.py`,
   `fingerprint.py`, `write_run_manifest.py`, `run_revision8_checks.py`,
   `diff_revision7_prefix.py` and `compare_cohorts.py`. Their tests:
   `GRM=<built grm> FOLDGEN=<built foldgen> python3 -I
   scripts/mtb-4k/tests/test_<name>.py` (the K review harnesses
   `test_k_review.py` and `test_k_review_r3.py` also take `--receipt OUT`;
   `test_k_review.py`'s frozen-evidence test runs only where that evidence
   is, on Vast).
4. **Resources, immediately before every stage.** At least 45 GiB free disk
   and 64 GiB reclaim-aware available RAM. The process-group guard keeps its
   31 GiB disk and 48 GiB RAM stop floors, above the 30/40 GiB global
   reserves. Re-measure; never reuse an earlier reading. Root preparation
   checks the same 45 GiB itself: it sums the annotation bytes it will copy
   and refuses unless 45 GiB stay free afterwards (section 3).
5. **Matrix review before structure; revision-8 gate before training.**
   Structure needs a matrix-only review. Training and report generation need
   the revision-8 gate in section 7 and a separate parent go-ahead.

The 4,096-column cap is unchanged.

## 2. Frozen inputs and cohorts

- Cohort evidence: `/root/grm-research/grm-research-mtb-4k-cohort-r2`.
- Full 4k: `results/genomes-4000.txt`, 4,000 IDs, SHA-256 `80faef89571ecd2fd658204abd39fd4b55e5a86bfe8e4f379274162868f44eaf`.
- Prefix 2k: rows 1–2000, the exact byte prefix, SHA-256 `57af1c86c35fb0a4c62d3c0784a09804a919ba75941159d0c7b3cf56870aa2e3`.
- Added 2k: rows 2001–4000, SHA-256 `a980d9e014fb90972a46b3fc01e65a9654deb58935d3ddf0f54018c19e7def85`.
- FASTA: `results/fasta-manifest.tsv` (4,000 hashes, reproduced twice). FASTA are linked, never copied or modified; root preparation hashes each against the manifest, and `--verify-frozen` hashes them again after every stage.
- Frozen 2k reference matrix: `/root/grm-revision7-preflight/data/species/mycobacterium_tuberculosis/matrix` (2,000 columns, k=31, min_support=2, 315,444 patterns, 6,041,188 kept k-mers).
- Inputs linked read-only: the candidate root `/root/grm-extra-species-analysis-20261027/data` (`catalog`, `curated`, and the BV-BRC record files). Each training cohort gets lane-local copies of the cached annotations of its own genomes only (section 3).

Cohorts, each in its own data root, all built by `prepare_execution_roots.py`:

| Cohort | Genomes | Role |
|---|---|---|
| `pilot256`, `pilot256-stock` | rows 1–256 | the small biological pilot and its stock comparator (section 4); never trained, so no annotations |
| `prefix` | rows 1–2000 | the frozen 2k candidate; byte-compared with the frozen matrix |
| `added` | rows 2001–4000 | an independent cohort, not a suffix of the matrix files |
| `full` | all 4,000 | the 4k cohort |
| `matched2k` | 2,000 of the 4,000, chosen by the smallest SHA-256 of `grm-mtb-4k-matched2k-20261028:<id>`, in frozen order | the size-matched control (section 8) |

Prefix and full have 22 eligible datasets; added has 14 (amikacin,
capreomycin, delamanid, ethambutol, ethionamide, isoniazid, kanamycin,
moxifloxacin, ofloxacin, prothionamide, pyrazinamide, rifabutin, rifampin,
streptomycin). Eight are unavailable in the added cohort under the
25-per-class gate. `matched2k` eligibility is whatever the same gate gives on
its genomes; `compare_cohorts.py` derives and lists it. Amoxicillin (25
resistant in all 4k) is expected to fail it.

All training uses the same source, evaluation revision 8, `EvalConfig`,
seeds, learner set and min_support=2. Old revision-7 reports are historical
comparators only.

## 3. Root setup (future, after gates 1–4)

```bash
export WS=/root/grm-research/grm-mtb-4k-execution-src   # a reviewed checkout of the approved commit
export RUN=/root/grm-research/grm-mtb-4k-execution-r1
export SLUG=mycobacterium_tuberculosis
export P=$WS/scripts/mtb-4k
# The wiring lane's process-group guard (31 GiB disk / 48 GiB RAM floors), used read-only.
export GUARD=/root/grm-research/grm-mtb-4k-wiring-kimi-20261028/tools/resource_guard.py

test ! -e "$RUN"
python3 -c 'import os; p=os.statvfs("/root"); d=p.f_bavail*p.f_frsize; m={l.split(":",1)[0]:l.split(":",1)[1].strip() for l in open("/proc/meminfo")}; r=int(m["MemAvailable"].split()[0])*1024; print(d,r); assert d >= 45*1024**3 and r >= 64*1024**3'
python3 -I /root/grm-research/grm-research-mtb-4k-cohort-r2/scripts/verify_evidence.py --root /root/grm-research/grm-research-mtb-4k-cohort-r2
python3 -I "$P/prepare_execution_roots.py" "$RUN"

(cd "$WS" && CARGO_TARGET_DIR="$RUN/target" CARGO_BUILD_JOBS=16 RAYON_NUM_THREADS=16 cargo build --release --locked --offline)
export GRM="$RUN/target/release/grm"
sha256sum "$GRM" | tee "$RUN/logs/grm-binary.sha256"
# The committed foldgen, which rebuilds every saved test_fold (section 7).
(cd "$WS/scripts/metric-audit/foldgen" && CARGO_TARGET_DIR="$RUN/target-foldgen" cargo build --release --locked --offline)
export FOLDGEN="$RUN/target-foldgen/release/foldgen"
sha256sum "$FOLDGEN" | tee "$RUN/logs/foldgen-binary.sha256"
export SOURCE=$(git -C "$WS" rev-parse HEAD)   # the one --expected-source-commit for every check
```

Before any matrix stage, two read-only checks against the revision-7 2k
candidate must pass with these builds. Both were run by the
plan-and-validators lane on 2026-10-28 (its receipt), and are repeated here
because the execution build is new:

```bash
# fingerprint.py reproduces every revision-7 report's input_fingerprint from the frozen root
# (a view of links: the preflight root's matrix, counted.json and structure.json with the
# candidate's curated and catalog; see that lane's tools/vast_fingerprint_oracle.py).
python3 -I "$P/fingerprint.py" "$RUN/rev7-view/data" --check --out "$RUN/audit/rev7-fingerprints.json" || exit 1
# The metric audit rebuilds every saved revision-7 test_fold with this foldgen.
python3 -I "$WS/scripts/metric-audit/audit.py" "$RUN/audit/rev7-fold-rebuild.json" $SLUG --revision 7 \
  --root /root/grm-extra-species-analysis-20261027 --reshaped --foldgen "$FOLDGEN"
python3 -c 'import json,sys; r=json.load(open(sys.argv[1])); bad=[i for x in r["records"] for i in x.get("issues",[]) if "test_fold" in i["where"]]; print(len(r["records"]),len(bad)); assert len(r["records"])==22 and not bad' "$RUN/audit/rev7-fold-rebuild.json"
```

`prepare_execution_roots.py` derives every list from the frozen bytes and
checks the prefix/added hashes. It hashes every FASTA against the cohort's
FASTA manifest, then creates only new directories and FASTA links. In each
root, `catalog` and `curated` link to the candidate's read-only directories,
and `bvbrc` is a **real lane-local directory**: the two record files link to
the candidate, and `features/` holds copies of the candidate's cached
annotations for that cohort's own genomes only, each checked against
`frozen-inputs.json`. `grm learn` reads an annotation only for a genome in
the cohort it trains (`annotate::load_features`), so no other file is ever
read. The two pilot roots never train and get an empty `features/`.

**Disk (N1).** Revision 2 copied the whole cache, 1,868 files and
1,477,799,489 bytes, into all six roots: 8,866,796,934 bytes, which failed
the 45 GiB gate before the pilot. Of the 1,868 files only 191 belong to
cohort genomes. The copies are now those files only. Measured read-only on
Vast at 15:58 UTC on 2026-10-28 with the script's own selection:

| Cohort | Files | Bytes |
|---|---|---|
| prefix | 170 | 128,989,934 |
| added | 21 | 15,904,550 |
| full | 191 | 144,894,484 |
| matched2k | 101 | 76,661,462 |
| pilot256, pilot256-stock | 0 | 0 |
| **total** | 483 | **366,450,430** (0.34 GiB) |

The script refuses before it creates anything unless 45 GiB stay free after
the copy, and `preparation.json` records the bytes copied and the free space
before and after. The 45 GiB gate is still measured, not assumed: Vast had
45,345,652,736 bytes (42.2 GiB) free at that reading, below the gate before
any copy, so the pilot cannot start until space is freed.

`grm learn` fetches missing annotations from BV-BRC into the cohort's
lane-local directory only. It refuses a `bvbrc` or `features` link and stops
rather than write through it (`annotate::FeatureCacheOutsideRoot`), and it
writes reports only through real directories below the root and never
through a link (`confine`). So training never writes into the frozen
candidate root. Before any stage that writes a cohort root, and after it:

```bash
python3 -I "$P/prepare_execution_roots.py" --check-lane "$RUN/cohorts/<cohort>/data"   # any problem is a stop
python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN"                      # any change is a stop
```

`--check-lane` refuses a link anywhere but where preparation made one
(catalog, curated, the two record files and each FASTA) and any entry that is
not a real directory or regular file; it covers `grm report`, whose writes
GRM does not confine itself. `--verify-frozen` re-hashes the candidate's
catalog, curated, records and whole annotation cache, every FASTA in the
linked store, each cohort's `genomes.txt`, and checks that every FASTA link
still points where preparation put it.

Newly fetched annotations are not frozen. After training, record their
SHA-256 (`find $RUN/cohorts/*/data/bvbrc/features -newer $RUN/preparation.json -type f -exec sha256sum {} +`)
in the receipt, and treat their content as network input that changed.

## 4. The small frozen biological pilot (gate 2)

It is separate from the 2,000-genome prefix stage and uses its own two roots,
on rows 1–256 of the frozen order (SHA-256 recorded in `preparation.json`).
Width 128 gives two blocks, so a fold, a resume and the partial merge are all
exercised on real MTB data.

The kill is no longer a timed poll and a `pkill -f` pattern (N8): the pattern
also matched the guard's and `/usr/bin/time`'s command lines and killed the
guard, and nothing checked that the run was interrupted. `pilot_interrupt.py`
finds the one process whose command line starts with exactly
`$GRM --data $SP matrix`, waits until `run.json` records block 0 folded
(`folded_columns` = 128: its partial verified and its owned ksets deleted),
stops the process, checks that state, kills it and waits until it is gone.
It exits non-zero, and the pilot stops, when the run ends first, when the
point is not reached in time, or when the stopped state is not an
interrupted run.

```bash
G="python3 $GUARD"
SP=$RUN/cohorts/pilot256/data; ST=$RUN/cohorts/pilot256-stock/data
for d in "$SP" "$ST"; do python3 -I "$P/prepare_execution_roots.py" --check-lane "$d" || exit 1; done
# a) sparse, interrupted once block 0 is folded and block 1 is under way
$G --log $RUN/logs/pilot-sparse-1.log --resources $RUN/logs/pilot-sparse-1.json -- \
  /usr/bin/time -v "$GRM" --data "$SP" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128 &
python3 -I "$P/pilot_interrupt.py" --data "$SP" --grm "$GRM" --folded 128 --receipt $RUN/audit/pilot-interrupt.json || exit 1
wait   # the guard records the killed run and exits
ls $SP/species/$SLUG/ksets | wc -l          # record: block 0's owned ksets are gone, block 1's in progress
# a2) the invalidation copies, taken now, of the interrupted root and before the resume
mkdir $RUN/pilot-copies
for c in width forge link; do cp -a "$RUN/cohorts/pilot256" "$RUN/pilot-copies/$c"; done
# b) resume to completion
$G --log $RUN/logs/pilot-sparse-2.log --resources $RUN/logs/pilot-sparse-2.json -- \
  /usr/bin/time -v "$GRM" --data "$SP" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128
# c) stock comparator on the same 256 (needs ~9 GB of ksets at once; within the 45 GiB gate)
$G --log $RUN/logs/pilot-stock.log --resources $RUN/logs/pilot-stock.json -- \
  sh -c "'$GRM' --data '$ST' count --species $SLUG && '$GRM' --data '$ST' matrix --species $SLUG --min-support 2"
for f in matrix.json patterns.bin support.bin kmer_offsets.bin kmers.bin; do cmp $SP/species/$SLUG/matrix/$f $ST/species/$SLUG/matrix/$f || exit 1; done
"$GRM" --data "$SP" structure --species $SLUG; "$GRM" --data "$ST" structure --species $SLUG
cmp $SP/species/$SLUG/structure.json $ST/species/$SLUG/structure.json || exit 1
python3 -I "$P/validate_matrix.py" --matrix $SP/species/$SLUG/matrix --counted $SP/species/$SLUG/counted.json \
  --genomes $SP/species/$SLUG/genomes.txt --fasta-dir $SP/fasta --recount 16 --require-recount --out $RUN/audit/pilot-matrix-validation.json
python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN"
```

Then run the invalidation checks on the three copies from step a2, never on
the live root. Each copy is the interrupted state: block 0 folded (its owned
ksets deleted), block 1 not finished.

1. A width change (`--block-width 64` on the `width` copy) is refused while
   the width-128 run's partials exist ("unfinished disk-sparse run").
2. On the `forge` copy, block 0's partial is edited together with its
   receipt (one column loses a k-mer, the receipt's SHA-256 and byte count
   are recomputed). The resume must fail `verify_partial` with "does not
   hold that genome's kset", and leave the copy's partials in place. Revision
   2 forged the partial inside the kill window it chose, while block 0 was
   still active and its ksets present; there GRM correctly rebuilds an
   unverified partial from the ksets and finishes byte-equal to stock (N2).
   That is the expected outcome for an active block, so this check forges a
   folded one.
3. On the `link` copy, a `matrix` directory link and, separately, a
   `ksets/<id>.kset.part` link are each refused before any work, and their
   targets keep their bytes.

Acceptance, all required:

- `pilot_interrupt.py` exits 0: one process matched, the stopped state had
  `folded_columns` = 128, no saved matrix, block 0's partial and receipt,
  and the process was gone after the kill. The guard's resources record for
  step a exists and shows the killed run.
- Sparse and stock matrix files and `structure.json` are byte-identical.
- The resume completes. No owned kset is deleted before its block's verified
  partial, no `matrix-partials`, owned ksets or `pending` are left, and the
  frozen inputs are unchanged.
- `validate_matrix.py` passes, including 16 recounts.
- The three invalidation checks refuse as stated.
- Measured figures are recorded: per-genome kset bytes against the
  preflight's upper bound (12 KiB + 8 bytes per FASTA base), per-block
  partial bytes, matrix bytes, peak RSS, minimum free disk and RAM, and wall
  time. The 4k projection from them (all partials + one block's ksets +
  matrix) must fit in 15 GiB with 20% margin, or the plan stops for a new
  resource decision. It never lowers the reserve.

## 5. Matrix-only stages (future, separately authorized)

The sparse matrix runs at width 128, with 16 Rayon threads and under the
guard. **No `--allow-missing-catalog`**: the catalog must load, or the run
stops before counting.

The per-block preflight is now an upper bound, not an estimate. An uncounted
genome is bounded at 8 bytes per FASTA base plus 12 KiB, about 35 MB for a
4.4 Mbp MTB assembly, so about 4.5 GB per 128-genome block. Real 2k MTB ksets
are about 26–30 MB. Each kset and sidecar write also claims its exact length
from one reserve ledger shared by the counting threads. At a 45 GiB start,
about 6.3 GB of partials, one block's 4.5 GB bound and the about 0.4 GB
matrix fit inside the 15 GiB above the reserve. The pilot replaces these
with measured figures.

```bash
for cohort in prefix added full matched2k; do
  D=$RUN/cohorts/$cohort/data
  python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || exit 1
  python3 "$GUARD" --log $RUN/logs/$cohort-matrix.log --resources $RUN/logs/$cohort-matrix-resources.json -- \
    /usr/bin/time -v "$GRM" --data "$D" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128
  python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
    --genomes $D/species/$SLUG/genomes.txt --fasta-dir $D/fasta --recount 16 --require-recount \
    --out $RUN/audit/$cohort-matrix-validation.json || exit 1
  python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || exit 1
done
for f in matrix.json patterns.bin support.bin kmer_offsets.bin kmers.bin; do
  cmp $RUN/cohorts/prefix/data/species/$SLUG/matrix/$f /root/grm-revision7-preflight/data/species/$SLUG/matrix/$f || exit 1
done
# Containment (N3): each 2k cohort against the full cohort's matrix, after all four exist.
F=$RUN/cohorts/full/data/species/$SLUG/matrix
for cohort in prefix added matched2k; do
  D=$RUN/cohorts/$cohort/data
  python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
    --genomes $D/species/$SLUG/genomes.txt --superset "$F" --out $RUN/audit/$cohort-in-full.json || exit 1
done
```

Run the cohorts one at a time, prefix first. Any failure stops the later
ones. Re-check the resource gate before each.

`validate_matrix.py` checks binding to `counted.json` and `genomes.txt`. It
checks structure: non-empty patterns in first-appearance order, every k-mer
inside the 2k-bit space and in exactly one pattern, distinct rows, the
support filter, and that distinct equals kept plus rare plus core. It checks
every column against its genome's bottom sketch. It recounts 16 seeded
genomes plus any `--recount-ids` from FASTA, independently of GRM; each must
reproduce `counted.json`'s k-mer count, FNV kset digest and sketch, and then
its whole column.

Revision 2 called that recount "the independent content check" for added,
full and matched2k, but it never looked at the k-mers a matrix leaves out
(N3): a dropped kept pattern, a dropped kept k-mer, a kept k-mer moved to the
rare count and a matrix built at min_support 3 but labelled 2 all passed,
even with every genome recounted. Three checks now cover omissions:

- **Exact region.** Below the smallest sketch threshold every genome's
  membership is known from its sketch, so support is exact there and the kept
  set must be exactly the k-mers with support in [2, n − 2].
- **Recount bounds.** A recounted k-mer the matrix leaves out must be rare or
  core. One that at least two known genomes hold and at least two lack
  (recounted genomes, plus sketches where they reach) is neither, so it is a
  missing kept k-mer. With every genome recounted, this is exact, and the
  rare, core and distinct counts must equal the recount's.
- **Containment (`--superset`).** A k-mer kept in a sub-cohort is kept in the
  full cohort, and a k-mer kept in the full cohort whose support over a
  sub-cohort's genomes is in [2, n − 2] is kept there. So prefix, added and
  matched2k must each equal the full matrix's rows restricted to their
  columns and filtered at min_support 2, k-mer for k-mer and bit for bit,
  with the same kset digests. The prefix is byte-identical to the frozen stock
  matrix, so the full matrix is checked against stock-built content on every
  k-mer the prefix keeps; added and matched2k are checked against the full
  matrix on every k-mer either keeps. A k-mer kept in the full matrix and in
  none of the three sub-cohorts is checked only by the exact region and the
  recount bounds.

After each cohort, check the receipt: no `pending`, no `matrix-partials`, no
owned ksets, exact order and digests, and the resource minima.

### Duplicate screen (after the prefix and added matrices, before training)

```bash
python3 -I "$P/screen_duplicates.py" --prefix-root $RUN/cohorts/prefix/data --added-root $RUN/cohorts/added/data \
  --fasta-manifest /root/grm-research/grm-research-mtb-4k-cohort-r2/results/fasta-manifest.tsv \
  --genomes-jsonl $RUN/cohorts/full/data/bvbrc/genomes.jsonl --out $RUN/audit/duplicate-screen.json
```

The screen looks for exact duplicates (same FASTA SHA-256, same kset digest,
same non-empty BioSample or assembly accession) and near duplicates (bottom-
sketch Jaccard ≥ 0.998), within each cohort and across the two.
`cross_cohort_exact > 0` is a stop: the parent decides whether the cohorts
stand. Near duplicates are listed. The primary lineage split groups them,
and the comparator records the screen.

## 6. Structure (separate go-ahead after matrix review)

```bash
for cohort in prefix added full matched2k; do
  python3 "$GUARD" --log $RUN/logs/$cohort-structure.log --resources $RUN/logs/$cohort-structure-resources.json -- \
    "$GRM" --data "$RUN/cohorts/$cohort/data" structure --species $SLUG
done
```

## 7. Training and reports: the revision-8 gate

Training, reports and the comparison run only when **all** of these hold:

1. **The reviewed revision-8 lasso guards and checker are in the binary.**
   They are being finished on the shared repo's
   `recovery/working-tree-20261028` branch. When the branch was rewritten at
   12:13 UTC on 2026-10-28 to remove attribution trailers, every commit ID
   changed; the old-to-new map is at the end of
   `docs/work-program-20261027.md`. At `3483065` (fetched 15:32 UTC) the last
   code commits are `43c187b` (lasso scan build; was `6a46127`) on top of
   `855bdba` (closes the revision-8 review findings; was `5d8cd05`) on top
   of `fb7ba09` (the migration guards; was `bf01371`). The wiring commits
   are now `fd559c9` (was `dd1f6d3`) and `79158f9` (was `176e3cc`).
   `79158f9` alone does not contain the lasso commits; the merge M0 does.
   **None of these is lasso-final.** The large-species lasso scan from
   `43c187b` stopped with a no-go on MTB delamanid (fold 0 past its budget),
   and a diagnosis lane is open. Pick up the first commit on that branch that
   the lasso lane records as final, that passes its scan rule, and whose
   checker an independent review records as approved; expect it after
   `43c187b`. Do not build or adopt a parallel revision-8 checker here.
2. **The training binary comes from a reviewed merge** of that commit with
   the wiring commits, the N1–N8 commit and the validator commit, in that
   order, trailer-free. Re-run `cargo test` and the wiring adversarial checks
   on the merge, then check the validators' pins: `report_contract.PINS`
   holds the SHA-256 of `check-results.py`, `metric-audit/audit.py`,
   `summarize.py` and `result_scope.py` as reviewed at `855bdba`. If
   lasso-final moves any of them, or adds or renames a solver-status field,
   the contract stops until the pins and its exact-key checks are reviewed
   again; that stop is intended.
3. **Every cohort passes `run_revision8_checks.py`.** It is not a parallel
   checker. Before the checker runs, the run manifest and every report must
   pass `report_contract.py`, which includes each report's input fingerprint
   recomputed from the root. Then it runs the in-tree
   `scripts/check-results.py` (zero unconverged two-class lasso fits, valid
   per-fold `lasso_fits` and `models.lasso_fit`), `scripts/metric-audit/audit.py
   --revision 8` with this cohort's backup, log and launch timestamp and the
   `foldgen` built above, which rebuilds every saved test_fold, and
   `summarize.py`. Only when all three pass and K's postflight holds does it
   write the receipt that `compare_cohorts.py` requires and parses (F5): a
   placeholder, a failed run, another run manifest, other pins or an output
   changed afterwards is a stop.
4. **An invalid report is a stop.** On revision 8, `grm learn` stops a species
   when a report's `evaluation.valid` is not true. Record the dataset and its
   reasons, and stop that cohort's training for parent review. Any solver
   status field the checker does not recognize is a stop, not something to
   ignore.

```bash
for cohort in prefix added full matched2k; do
  C=$RUN/cohorts/$cohort; D=$C/data
  python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || exit 1
  mkdir "$RUN/backups/$cohort"                       # pre-launch results: none in a new root
  LAUNCH=$(date -u +%Y-%m-%dT%H:%M:%SZ); echo "$LAUNCH" > "$RUN/logs/$cohort-learn.launch"
  python3 "$GUARD" --log $RUN/logs/$cohort-learn.log --resources $RUN/logs/$cohort-learn-resources.json -- \
    /usr/bin/time -v "$GRM" --data "$D" learn --species $SLUG || exit 1
  python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || exit 1
  "$GRM" --data "$D" report || exit 1                # dashboard data in this root only
  python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || exit 1
  python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || exit 1
  python3 -I "$P/write_run_manifest.py" "$C" --cohort $cohort --source-commit "$SOURCE" --binary "$GRM" \
    --log $RUN/logs/$cohort-learn.log || exit 1
  python3 -I "$P/run_revision8_checks.py" $RUN/audit/$cohort-audit.json --root "$C" --run-manifest "$C/validation-manifest.json" \
    --expected-source-commit "$SOURCE" --binary "$GRM" --backup "$RUN/backups/$cohort" --log $RUN/logs/$cohort-learn.log \
    --launch "$LAUNCH" --foldgen "$FOLDGEN" --checker-out $RUN/audit/$cohort-check-results.json \
    --summary-out $RUN/audit/$cohort-summary.json --receipt $RUN/audit/$cohort-revision8-receipt.json || exit 1
done
# The prefix against the revision-7 2k candidate, field by field (only the lasso may differ).
python3 -I "$P/diff_revision7_prefix.py" --revision8 $RUN/cohorts/prefix/data \
  --revision7 /root/grm-extra-species-analysis-20261027/data --out $RUN/audit/prefix-vs-revision7.json || exit 1
```

`write_run_manifest.py` binds each cohort's inputs, reports and
class-gate receipts by SHA-256, recomputing every input fingerprint and
every dataset's class counts from the root first. Expected report counts
are 22 prefix, 14 added and 22 full, and whatever the class gate gives
matched2k (14 if the frozen evidence holds: amoxicillin and the other seven
sparse drugs fail it there too). MTB's unsupported probe assays must not
trigger a stock recount. Build dashboard data only in the new roots, each
after `--check-lane`. Do not publish, deploy, replace candidate reports or
touch the live site.

`diff_revision7_prefix.py` is the field-level check that
`docs/lasso-stopping.md` asks for: the prefix matrix is byte-identical to
the frozen 2k matrix, so every report's input fingerprint, samples, catalog,
panels, MIC, splits and non-lasso learner rows must equal the revision-7
candidate's; only the lasso rows and models may differ.

## 8. Matched analysis (after section 7)

```bash
python3 -I "$P/compare_cohorts.py" \
  --cohort prefix=$RUN/cohorts/prefix --cohort added=$RUN/cohorts/added \
  --cohort full=$RUN/cohorts/full --cohort matched2k=$RUN/cohorts/matched2k \
  --checker-receipt prefix=$RUN/audit/prefix-revision8-receipt.json --checker-receipt added=$RUN/audit/added-revision8-receipt.json \
  --checker-receipt full=$RUN/audit/full-revision8-receipt.json --checker-receipt matched2k=$RUN/audit/matched2k-revision8-receipt.json \
  --duplicates $RUN/audit/duplicate-screen.json --expected-source-commit "$SOURCE" --out $RUN/audit/compare-cohorts.json
```

Exit 0 means complete. Exit 1 means output was written but the comparison is
invalid (source, binary, foldgen or EvalConfig differing across cohorts, a
coverage or label mismatch, cohorts not disjoint). Exit 2 means an input or
integrity stop with no output: a manifest, hash or fingerprint binding, a
receipt that does not prove its cohort passed (F5), the matched2k list not
being `matched_selection` of the full cohort (F6), a class-gate receipt or
report count that the full cohort's labels on that cohort's genomes do not
reproduce (F3), an unexpected not-interpretable set, or an unaccepted exact
duplicate. Exit 3 means output was written and valid, but a drug was
excluded or a primary split was unavailable; stop for review before any
before/after statement.

**Decision D1, pending the user's override.** grm's random folds are a
seeded stratified shuffle over each cohort's own samples, so the same genome
lands in different random folds in the prefix and full cohorts (28–62% of
common genomes in the reconciliation's synthetic run; about 80% expected at
2k/4k). Package K treated any such mismatch as blocking, which would stop
every real comparison. The canonical comparator counts the mismatch per
row, keeps those genomes out of `matched_fold`, and marks that row's matched
random control invalid, but the comparison stays valid. Fold integrity is
not weakened: `run_revision8_checks.py` rebuilds every saved fold of every
cohort with the real `foldgen` and stops on any difference.

Prespecified rules:

- **Primary split: lineage** (grouped by `structure.json` clusters, which
  keeps near-duplicate isolates on one side of a fold). The random split is
  secondary. Lineage clusters are recomputed per cohort, so folds and groups
  differ between cohorts. Paired results align by genome ID, never by row or
  fold. A drug whose lineage split is unavailable in a cohort has no primary
  result there: it is listed under `primary_split_unavailable`, left out of
  the size-matched decomposition, and the run exits 3 (N6).
- **Primary measure:** per drug, the median of the five learners' balanced
  accuracy at 0.5 on the lineage split, only when all five have metrics.
  Per-learner rows with Clopper–Pearson 95% intervals, VME/ME, no-calls and
  per-fold denominators are recomputed from the saved held-out scores. A
  stored row that disagrees with its scores invalidates the report.
- **A drug whose report is invalid, or missing although eligible, in only
  one cohort** is dropped from every comparison that involves that cohort
  and from the size-matched decomposition. It is never imputed. It is listed
  under `excluded` with the cohort and reasons, and stays in the comparisons
  between the other cohorts. Exit 3 stops the run for review, and no
  before/after statement is written until the parent decides. Ineligibility
  under the class gate is not an exclusion; it is listed under
  `not_eligible`.
- **Not interpretable for before/after claims:** amoxicillin,
  ciprofloxacin, clofazimine, cycloserine, linezolid and nicotinamide (the
  added 2,000 hold no resistant isolate), para-aminosalicylic acid (5
  resistant) and levofloxacin (774 resistant, 1 susceptible). For these
  drugs a full-versus-prefix difference mostly reflects prevalence and
  composition, not what the added genomes teach. They are reported
  descriptively, flagged `before_after_interpretable: false`, and left out of
  any before/after statement and of multiplicity adjustment. The comparator
  derives this set from the full cohort's labels and stops if it differs
  from these eight.
- **Size-matched control:** `matched2k` is a random 2,000 of the 4,000, fixed
  by hash before any result exists. Per drug, the change from prefix to full
  is split into composition (matched2k minus prefix, at equal n) and size
  (full minus matched2k, at fixed composition). An improvement shown by
  full-versus-prefix but not by full-versus-matched2k is not attributed to
  sample size.
- **Comparisons:** full against prefix on the 22, paired on the prefix
  isolates, with fixed and broken calls at 0.5 and score deltas. Added
  against prefix on the 14, side by side as disjoint cohorts: this is
  cross-validation, not external validation or frozen-model scoring.
  Matched2k against prefix, side by side. Full against matched2k, paired.
- **Multiplicity:** exact McNemar p-values on fixed against broken calls are
  descriptive. They are Holm-adjusted within each learner and split, across
  the interpretable drugs. No confirmatory claim is made.
- **Scope:** internal cross-validation on public BV-BRC genomes. No external,
  prospective, clinical or causal claim.

Keep every command, exit code, source, binary and input hash, validator
receipt and resource measurement. A numerical result is unknown until these
steps run and pass.

## 9. Not covered

- Real power loss and filesystem-specific ordering on the Vast overlay.
  Durability was tested through fsync and rename order, not power cuts.
- Real 4k kset, partial and RAM sizes; the pilot measures them.
- Revision-8 training, checking and reporting at 4k. These wait for section 7.
- The data root and its ancestors are trusted, and bind mounts are not
  detected (N7). `grm fetch` is not confined; the plan does not run it.
