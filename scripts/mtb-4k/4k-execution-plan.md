# MTB matched 2k / added-2k / 4k execution plan — NOT AUTHORIZED

Revision 2 (2026-10-28, fix lane). It replaces the wiring lane's plan
(`plans/4k-execution-plan.md`, SHA-256 `8bb780c1…6fb8`) after the independent
review `WIRING_REVIEW=CHANGES_REQUIRED`. The code it runs is commit
`dd1f6d3` plus the fix commit that adds this directory; that pair has not
been re-reviewed.

**No command below has been run. No biological matrix, structure, training,
report or 8k work exists. No 8k stage is planned or permitted.**

## 1. Gates

`RUN_4K=YES` only when all of these hold; otherwise `RUN_4K=NO`.

1. **Wiring re-review.** A fresh independent reviewer inspects the fix commit
   on top of `dd1f6d3` and records `WIRING_REVIEW=APPROVED`.
2. **Pilot.** The parent authorizes the small frozen biological pilot
   (section 4), and every pilot acceptance check passes before the prefix
   stage.
3. **Validators reviewed.** `validate_matrix.py`, `screen_duplicates.py`,
   `compare_cohorts.py` and `prepare_execution_roots.py` in this directory
   are independently reviewed. Their tests: `python3 -I
   scripts/mtb-4k/tests/test_<name>.py` (the validator's test builds real
   matrices and needs `GRM=<built grm binary>`).
4. **Resources, immediately before every stage.** At least 45 GiB free disk
   and 64 GiB reclaim-aware available RAM. The process-group guard keeps its
   31 GiB disk and 48 GiB RAM stop floors, above the 30/40 GiB global
   reserves. Re-measure; never reuse an earlier reading.
5. **Matrix review before structure; revision-8 gate before training.**
   Structure needs a matrix-only review. Training and report generation need
   the revision-8 gate in section 7 and a separate parent go-ahead.

The 4,096-column cap is unchanged.

## 2. Frozen inputs and cohorts

- Cohort evidence: `/root/grm-research/grm-research-mtb-4k-cohort-r2`.
- Full 4k: `results/genomes-4000.txt`, 4,000 IDs, SHA-256 `80faef89571ecd2fd658204abd39fd4b55e5a86bfe8e4f379274162868f44eaf`.
- Prefix 2k: rows 1–2000, the exact byte prefix, SHA-256 `57af1c86c35fb0a4c62d3c0784a09804a919ba75941159d0c7b3cf56870aa2e3`.
- Added 2k: rows 2001–4000, SHA-256 `a980d9e014fb90972a46b3fc01e65a9654deb58935d3ddf0f54018c19e7def85`.
- FASTA: `results/fasta-manifest.tsv` (4,000 hashes, reproduced twice). FASTA are linked, never copied or modified.
- Frozen 2k reference matrix: `/root/grm-revision7-preflight/data/species/mycobacterium_tuberculosis/matrix` (2,000 columns, k=31, min_support=2, 315,444 patterns, 6,041,188 kept k-mers).
- Inputs linked read-only: the candidate root `/root/grm-extra-species-analysis-20261027/data` (`catalog`, `curated`, and the BV-BRC record files). Its annotation cache is copied, not linked (section 3).

Cohorts, each in its own data root, all built by `prepare_execution_roots.py`:

| Cohort | Genomes | Role |
|---|---|---|
| `pilot256`, `pilot256-stock` | rows 1–256 | the small biological pilot and its stock comparator (section 4) |
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
```

`prepare_execution_roots.py` derives every list from the frozen bytes and
checks the prefix/added hashes. It creates only new directories and
manifest-verified FASTA links. In each root, `catalog` and `curated` link to
the candidate's read-only directories, and `bvbrc` is a **real lane-local
directory**: the two record files link to the candidate, and `features/` is a
copy of the candidate's annotation cache, checked file by file against
`frozen-inputs.json`. `grm learn` fetches missing annotations from BV-BRC
into that lane-local directory only. It refuses a `bvbrc` or `features` link
and stops rather than write through it (`annotate::FeatureCacheOutsideRoot`).
So training never writes into the frozen candidate root. After each stage:

```bash
python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN"   # any change is a stop
```

Newly fetched annotations are not frozen. After training, record their
SHA-256 (`find $RUN/cohorts/*/data/bvbrc/features -newer $RUN/preparation.json -type f -exec sha256sum {} +`)
in the receipt, and treat their content as network input that changed.

## 4. The small frozen biological pilot (gate 2)

It is separate from the 2,000-genome prefix stage and uses its own two roots,
on rows 1–256 of the frozen order (SHA-256 recorded in `preparation.json`).
Width 128 gives two blocks, so a fold, a resume and the partial merge are all
exercised on real MTB data.

```bash
G="python3 $GUARD"
SP=$RUN/cohorts/pilot256/data; ST=$RUN/cohorts/pilot256-stock/data
# a) sparse, with a kill once block 0's partial is verified and before its ksets go
$G --log $RUN/logs/pilot-sparse-1.log --resources $RUN/logs/pilot-sparse-1.json -- \
  /usr/bin/time -v "$GRM" --data "$SP" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128 &
until ls $SP/species/$SLUG/matrix-partials/*/block-00000.partial.receipt.json >/dev/null 2>&1; do sleep 1; done
pkill -KILL -f -- "--data $SP matrix"; wait
ls $SP/species/$SLUG/ksets | wc -l          # record: block 0's ksets are still there
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

Then run the invalidation checks on copies made with `cp -a` of the
interrupted pilot root, never on the live root:

1. A width change mid-run (`--block-width 64` on a copy stopped after block 0)
   is refused while the width-128 run's partials exist.
2. A partial edited together with its receipt fails `verify_partial`, with
   "does not hold that genome's kset".
3. A `matrix` link or a `ksets/<id>.kset.part` link is refused before any
   work.

Acceptance, all required:

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
```

Run the cohorts one at a time, prefix first. Any failure stops the later
ones. Re-check the resource gate before each.

`validate_matrix.py` checks binding to `counted.json` and `genomes.txt`. It
checks structure: non-empty patterns in first-appearance order, every k-mer
inside the 2k-bit space and in exactly one pattern, distinct rows, the
support filter, and that distinct equals kept plus rare plus core. It checks
every column against its genome's bottom sketch. It also recounts 16 seeded
genomes plus any `--recount-ids` from FASTA, independently of GRM; each must
reproduce `counted.json`'s k-mer count, FNV kset digest and sketch, and then
its whole column. For the added, full and matched2k cohorts, this recount is
the independent content check; the prefix also has the byte comparison.

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
   `recovery/working-tree-20261028` branch. At this writing its head is
   `2a9f85f`. The last code commits there are `6a46127` (lasso scan build) on
   top of `5d8cd05` (closes the revision-8 review findings) on top of
   `bf01371` (the migration guards). All three descend from `33ccd2a`, so
   **neither `dd1f6d3` nor this fix commit contains them**. Pick up the first
   commit on that branch that the lasso lane records as final and whose
   checker an independent review records as approved; expect it at or after
   `6a46127`. Do not build or adopt a parallel revision-8 checker here.
2. **The training binary comes from a reviewed merge** of that commit with
   `dd1f6d3` and this fix commit. Re-run `cargo test` and the wiring
   adversarial checks on the merge.
3. **Every report passes the reviewed checker:** `scripts/check-results.py`
   on that commit validates the revision-8 lasso solver status (per-fold
   `lasso_fits` and `models.lasso_fit`). It must exit 0 with zero unconverged
   two-class lasso fits, per cohort. `scripts/metric-audit/audit.py
   --revision 8` must also exit 0, per cohort, with its own backup, log and
   launch timestamp, and `foldgen` built in the execution target. Save each
   checker's output; `compare_cohorts.py` refuses to run without one per
   cohort.
4. **An invalid report is a stop.** On revision 8, `grm learn` stops a species
   when a report's `evaluation.valid` is not true. Record the dataset and its
   reasons, and stop that cohort's training for parent review. Any solver
   status field the checker does not recognize is a stop, not something to
   ignore.

```bash
for cohort in prefix added full matched2k; do
  python3 "$GUARD" --log $RUN/logs/$cohort-learn.log --resources $RUN/logs/$cohort-learn-resources.json -- \
    /usr/bin/time -v "$GRM" --data "$RUN/cohorts/$cohort/data" learn --species $SLUG
  python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || exit 1
done
```

Expected report counts are 22 prefix, 14 added and 22 full, and whatever
the class gate gives matched2k. MTB's unsupported probe assays must not
trigger a stock recount. Build dashboard data only in the new roots. Do not
publish, deploy, replace candidate reports or touch the live site.

## 8. Matched analysis (after section 7)

```bash
python3 -I "$P/compare_cohorts.py" \
  --cohort prefix=$RUN/cohorts/prefix/data --cohort added=$RUN/cohorts/added/data \
  --cohort full=$RUN/cohorts/full/data --cohort matched2k=$RUN/cohorts/matched2k/data \
  --checker-receipt prefix=$RUN/audit/prefix-check-results.json --checker-receipt added=$RUN/audit/added-check-results.json \
  --checker-receipt full=$RUN/audit/full-check-results.json --checker-receipt matched2k=$RUN/audit/matched2k-check-results.json \
  --duplicates $RUN/audit/duplicate-screen.json --out $RUN/audit/compare-cohorts.json
```

Exit 0 means complete. Exit 2 means an input or integrity stop, with no
output. Exit 3 means output was written but a drug was excluded; stop for
review before any before/after statement.

Prespecified rules:

- **Primary split: lineage** (grouped by `structure.json` clusters, which
  keeps near-duplicate isolates on one side of a fold). The random split is
  secondary. Lineage clusters are recomputed per cohort, so folds and groups
  differ between cohorts. Paired results align by genome ID, never by row or
  fold.
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
