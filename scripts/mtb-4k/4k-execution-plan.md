# MTB matched 2k / added-2k / 4k execution plan: revision 6, cleaned cohorts

Revision 6 (2026-10-29, plan-r6 lane). Revision 5 (SHA-256 `2807de98…`) ran the frozen
cohorts on revision 8 from `c70d678`. The matrix review then found two genomes in the frozen
prefix that are not *M. tuberculosis*: `1773.41701` (most likely *M. abscessus*, `blaMAB`) and
`1773.40949` (GC 68.4%). Both are also in full and matched2k. Revision 6 runs the **cleaned
cohorts** on evaluation revision N, the first revision whose data pipeline carries the species
check. The admin's rulings of 2026-10-29 set its shape:

- **The cohorts are revision N's fresh selection.** After the species check fails the two genomes,
  `grm curate` reruns selection on the cleaned labels. The cleaned 2k is revision N's own MTB 2k
  `genomes.txt`, not the frozen prefix minus two: replacements fill the cap, and the order may
  shift near score ties. The cleaned full cohort is revision N's cap-4000 selection. **Its first
  rows must be the cleaned 2k exactly**, and a difference is a stop before any root exists.
  Added is the rows after the prefix, and matched2k is the same smallest-SHA-256 rule on the
  cleaned full list, at the prefix's size (section 2).
- **Replacements must pass, in order, the species check, curation, the FASTA fetch, and then the
  species check again**, before they enter any matrix. The plan binds the second (post-fetch)
  species check by hash. Root preparation refuses any cohort genome that check fails, has no
  verdict for, or lists as unchecked, and any selection whose `species.json` reports
  `selected_unchecked` other than 0. Borderline genomes are recorded and never excluded.
- **When it runs.** The species check's review found defects, and it moves to revision 10 with
  duplicate isolation (the orchestrator, 2026-10-29). So the revision-9 rerun of the 4k cohorts is
  uncleaned (`b2fd433`, under revision 5's plan with `--profile rN`). This plan is expected to run
  after it, on N = 10: branch (b) of section 8. The alternative, branch (a), was a single cleaned
  revision-9 run with no with/without comparison; it applies only if the check lands in the
  revision that reruns the cohorts. N stays a placeholder. Nothing here assumes 9 or 10, and
  nothing assumes a run came before.
- **Revision N's MTB 2k candidate is a reference checked after training, not a gate before it.**
  The prefix `cmp` against its matrix, the fingerprint view, the fold rebuild and the field-level
  diff against its reports (sections 7b and 9) run once that candidate exists. They are recorded
  as `PENDING REFERENCE` until then, and they are required for the final result. Nothing else is
  relaxed.

Every value that does not exist yet is an `@@NAME@@` placeholder in one of two files (section 2).
`revision_n.py` holds revision N's validator values, and `cleaned_rerun.py` holds the cleaned
cohorts' inputs. Each script refuses, naming every unfilled value it needs, until a reviewed commit
fills them. Every command that loads the report contract (`write_run_manifest.py`,
`run_revision8_checks.py`, `compare_cohorts.py`) passes `--profile rN-cleaned`, because the
default `r8` cannot load on a revision-N tree. The code is four commits on `c70d678`:
- `f6c2b87`, the revision-N validators, with R1, which removes the every-step-capped convergence
  rule;
- `e354f99`, the validator suite under `MTB4K_TEST_PROFILE`;
- `8bb1a9b`, `diff_revision7_prefix.py --profile`;
- this plan and its scripts. Section 1's gates stand, read for revision N. **No revision-6 command has been run.**

## Earlier revisions

Revision 5 (2026-10-28, plan re-check fixes lane). The focused re-check of
revision 4 (`5d2c24d`) returned `PLAN_VALIDATORS_REVIEW=APPROVED` for the
small pilot of section 4 only, with Low findings N1–N4 and Info items
N5–N13. In this revision `pilot_interrupt.py` accepts as `--launcher` only
this run's guard, checked by its command line and `--resources` before any
signal (N1). After any stop it lists every `$GRM --data $SP matrix` still
running and signals none, and step a lists them again after its `wait`
(N2). Both section-3 `python3 -c` checks run with `-I`, and
`tests/test_entry_points.py` reads every fence and every Python command
shape, `-c` included (N3, N7). A test now fails if the invalidation check's
"the refused run changed nothing" or the view's existing-path refusal is
removed (N4). A signal during the stop after a taken `--receipt` path no
longer aborts it (N5). `--verify-frozen` and `screen_duplicates.py` stop
with a line, not a traceback, on missing inputs (N13). Section 9 lists what
is deferred. No change turns a stop into a pass.

Revision 4 (2026-10-28, plan follow-ups lane). The independent review of
revision 3 (`ffa43a9`) returned `PLAN_VALIDATORS_REVIEW=APPROVED` for the small
pilot of section 4 only, with findings M1 and L1–L6. In this revision every
`python3 -I` command loads its own modules, and `tests/test_entry_points.py`
runs each one as written here (M1). The pilot's interrupt is bound to the
launch it interrupts and stops that launch whenever it stops (L1). The
comparator re-reads the outputs a checker receipt binds (L2), and any invalid
report stops it (L3). Dated statements are dated (L4). The revision-7 view and
the pilot's invalidation checks have commands (L5). The gates name the review
lines that exist and state `RUN_4K=NO` directly (L6). The foldgen build is
pinned by hash (I4). No change turns a stop into a pass.

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
branch at `3483065`, then the N1–N8 commit, then the validator commit
`ffa43a9` (reviewed as above), then revision 4's follow-up commit `5d2c24d`
(approved for the pilot only by its focused re-check), then revision 5's
commit, which needs its own re-check. M0 is **not** the lasso-final merge
(section 7).

Revision 5's statement that no command had run was true when it was written; the current run has
since executed revision 5. No 8k stage is planned or permitted.

## 1. Gates

**`RUN_4K=NO` until a lasso-final revision-8 commit exists and an
independent review records its checker as approved (section 7, item 1).** No
gate below overrides that. The gates add to it: `RUN_4K=YES` only when that
holds and all of these hold; otherwise `RUN_4K=NO`.

1. **Reviews, by the gate lines that exist.** Both are required:
   - `WIRING_REVIEW=APPROVED`, from the wiring review lineage, for the Rust in
     the binary that runs. Its last line is the re-review of `79158f9`:
     `WIRING_REVIEW=CHANGES_REQUIRED`, for plan and script findings N1–N8,
     with no open Rust finding. The plan-validators review of `ffa43a9`
     re-checked every Rust finding of both wiring reviews and found them
     closed. The training binary of section 7 is a new merge, so it needs this
     line from that lineage's review of the merge.
   - `PLAN_VALIDATORS_REVIEW=APPROVED`, from the plan-and-validators review
     lineage, for the commit that carries this plan and `scripts/mtb-4k`.
     `ffa43a9` and revision 4's follow-up commit `5d2c24d` each have it for
     the section-4 pilot only (2026-10-28). Revision 5's commit needs its
     own, from a re-check, before any stage that runs revision-5 commands.
2. **Pilot.** The parent authorizes the small frozen biological pilot
   (section 4), and every pilot acceptance check passes before the prefix
   stage.
3. **Validators reviewed.** Every script in this directory is
   independently reviewed: `prepare_execution_roots.py`, `pilot_interrupt.py`,
   `pilot_invalidation.py`, `validate_matrix.py`, `screen_duplicates.py`,
   `report_contract.py`, `fingerprint.py`, `write_run_manifest.py`,
   `run_revision8_checks.py`, `diff_revision7_prefix.py` and
   `compare_cohorts.py`. Their tests: `GRM=<built grm> FOLDGEN=<built foldgen>
   python3 -I scripts/mtb-4k/tests/test_<name>.py` (the K review harnesses
   `test_k_review.py` and `test_k_review_r3.py` also take `--receipt OUT`;
   `test_k_review.py`'s frozen-evidence test runs only where that evidence
   is, on Vast). `test_entry_points.py` runs every `python3 -I` command of
   this plan as written, from another directory, with the Python that will
   run them (M1); run it on Vast too.
4. **Resources, immediately before every stage.** At least 45 GiB free disk
   and 64 GiB reclaim-aware available RAM. The process-group guard keeps its
   31 GiB disk and 48 GiB RAM stop floors, above the 30/40 GiB global
   reserves. Re-measure; never reuse an earlier reading. Root preparation
   checks the same 45 GiB itself: it sums the annotation bytes it will copy
   and refuses unless 45 GiB stay free afterwards (section 3).
5. **Matrix review before structure; revision-8 gate before training.**
   Structure needs a matrix-only review. Training and report generation need
   the revision-8 gate in section 7 and a separate parent go-ahead. That
   gate is the `RUN_4K=NO` hold stated above; nothing here relaxes it.

The 4,096-column cap is unchanged.

## 2. Frozen inputs and cohorts

Revision 6 binds every cohort input by hash. Values that do not exist yet are placeholders,
filled by a reviewed commit and listed in section 2.1.

- **Cohort evidence:** `/root/grm-research/grm-research-mtb-4k-cohort-r2`, as in revision 5. Its
  frozen 4k list (`80faef89…`), prefix (`57af1c86…`), added (`a980d9e0…`) and FASTA manifest are
  still verified at preparation: they are the comparison basis, and the source of every FASTA
  the cleaned cohorts share with them.
- **Post-fetch species check:** `species_check.json`, at `SPECIES_CHECK_PATH` with
  `SPECIES_CHECK_SHA256`. Its MTB `fail` genomes are the exclusions. They must include `1773.41701`
  and `1773.40949` and number `EXCLUSIONS_COUNT`. Its `borderline` genomes are recorded per cohort
  in `preparation.json` and never excluded. On 2026-10-29 the species-check lane reported six:
  `1773.18404`, `1773.3676`, `1773.18177`, `1773.40722`, `1773.18213` and `1773.5309` (the last in
  added).
- **Revision N's selections:** the 2k `genomes.txt` (`RN_PREFIX_GENOMES_*`) and the cap-4000
  selection (`RN_FULL_GENOMES_*`), each with its `curated/species.json` (`RN_*_SPECIES_JSON_*`),
  whose MTB `selected_unchecked` must be 0.
- **Inputs:** `RN_INPUT_ROOT` provides catalog, curated, the BV-BRC records and cached annotations.
  Five of its files are bound by hash and checked before root preparation (review P1, re-check N1):
  `curated/datasets.json`, `curated/phenotypes.jsonl` (the labels), `catalog/version.txt`,
  `catalog/AMR_CDS.fa` and `catalog/ReferenceGeneCatalog.txt`. The last two are what `grm matrix` loads
  for gene calls; GRM identifies a catalog only by `version.txt`. `RN_INPUT_MANIFEST_SHA256` then binds
  everything a cohort root can read from the input root, as a sorted `path<TAB>sha256` manifest:
  - all of `catalog/` and `curated/`, which are linked whole, so they include `funnel.json` and
    `drug_names.json`, which `grm report` reads;
  - the two BV-BRC record files (`genomes.jsonl` feeds the duplicate screen);
  - the cohort genomes' cached annotations, which `grm learn` reads.

  The set comes from the same rules that make the roots, and a test walks the prepared roots to
  check it. Nothing the cohorts read from the input root is bound only at preparation.
- **Annotations must all be cached; `grm learn` must fetch nothing.** `load_features`
  (`src/annotate.rs:101-117`) reads `bvbrc/features/<id>.json` and fetches the annotation from BV-BRC
  whenever it cannot read or parse that file. The revision-9 run fetched 230 annotations live
  during training, an unbound and non-reproducible input. So revision N's staging must cache an
  annotation for every selected genome in `RN_INPUT_ROOT`. Preparation refuses, before creating
  anything, if any genome of any cohort has none, naming up to 20 and giving the count. The cached
  annotations are part of the input manifest.
- **Staging, not this plan, completes the species check.** Revision N's cap-4000 entrants arrive
  unchecked (`selected_unchecked` > 0). Staging must fetch them and run the species check again until
  `selected_unchecked` is 0 for both selections; the species-check lane (024208eb) produces it. This
  plan only refuses otherwise.
  These are revision N's curation, on which the fresh selection was made. Its whole set is bound in
  `frozen-inputs.json` at preparation, as revision 5 bound the candidate's.
- **FASTA:** a genome in the frozen 4k list links to its frozen-manifest FASTA. A replacement links
  to its row in the bound replacement manifest (`RN_FASTA_MANIFEST_*`), a file directly under
  `RN_FASTA_STORE`. Both are hash-checked at preparation and by every `--verify-frozen`. The
  preparation writes `$RUN/fasta-manifest.tsv` over every cohort genome for the duplicate screen.
- **Reference (pending until it exists):** revision N's trained MTB 2k candidate under
  `RN_REFERENCE_ROOT` (`data/species/mycobacterium_tuberculosis/{genomes.txt,matrix,results}`),
  read only by section 7b.

| Cohort | Genomes | Role |
|---|---|---|
| `prefix` | revision N's MTB 2k selection | the cleaned 2k; after training, checked against revision N's candidate (section 7b) |
| `added` | the cap-4000 selection's rows after the prefix | the genomes full adds; disjoint from prefix |
| `full` | revision N's cap-4000 selection, whose first rows are the prefix | the cleaned 4k |
| `matched2k` | as many as the prefix, chosen from full by the smallest SHA-256 of `grm-mtb-4k-matched2k-20261028:<id>`, in full's order | the size-matched control (section 8) |

There are no pilot roots: section 4 is not repeated. Each cohort list must equal its `CLEANED_*`
pin (hash and count). `prepare_execution_roots.py --derive-cleaned` prints those pins from the
bound inputs once the rest is filled, and writes nothing.

**Report counts: frozen sets, derived counts.** The eligible drug sets stay frozen, with a check:
22 for prefix and full, the 14 of revision 5 for added, and matched2k whatever the class gate
gives. `report_contract.py`'s `rN-cleaned` rules require the report set to equal the class gate
recomputed from the root's own phenotypes, so a cleaned cohort that moves a drug across the
25-per-class gate is a stop. The eight not-interpretable drugs of section 8 and the added cohort's
14 are prespecified analysis sets, and `compare_cohorts.py` stops if cleaning changes either. The
per-drug class counts, including the added cohort's eight class-gate receipts, are derived from
the phenotypes, not frozen. The fresh selection's added rows differ from the frozen rows by more
than the exclusions, so revision 5's frozen counts no longer apply. They are still cross-checked
three ways: `grm learn`'s gate and skip lines, the contract's recomputation from the phenotypes,
and section 8's F3 check from the full cohort's labels.

All training uses one source, evaluation revision N, `EvalConfig`, seeds, learner set and
min_support=2. Revision-8 reports, the current run's included, are not cohorts of this plan.

### 2.1 Placeholders

`python3 -I scripts/mtb-4k/revision_n.py status` and `cleaned_rerun.py status` list each value and
whether it is filled. A filling commit changes `VALUES` only; the tests hold for a wholly unfilled
or a wholly filled file. "Prepare" values are needed before section 3, and "reference"
values only for section 7b.

| Placeholder | File | Needed for | What it is |
|---|---|---|---|
| `@@REVISION_N@@` | `revision_n.py` | prepare | evaluation revision N: the first that carries the species check (expected 10; not assumed) |
| `@@RN_CHECKER_COMMIT@@` | `revision_n.py` | prepare | revision N's checker commit; `verify-checkout` requires `$WS` to descend from it |
| `@@RN_PIN_CHECK_RESULTS@@`, `@@RN_PIN_METRIC_AUDIT@@`, `@@RN_PIN_SUMMARIZE@@`, `@@RN_PIN_RESULT_SCOPE@@` | `revision_n.py` | prepare | SHA-256 of the four checker files at that commit |
| `@@RN_FIELD_CHANGES@@` | `revision_n.py` | prepare | report or status fields revision N adds, `{}` for none; a drop or rename refuses |
| `@@RN_CONVERGENCE_STATUS_RULE@@` | `revision_n.py` | prepare | `as-690815a` if revision N's checker judges convergence statuses as 690815a's does (IRLS0, R1); anything else refuses |
| `@@SPECIES_CHECK_PATH@@`, `@@SPECIES_CHECK_SHA256@@` | `cleaned_rerun.py` | prepare | the post-fetch `species_check.json` |
| `@@SPECIES_CHECK_FORMAT@@` | `cleaned_rerun.py` | prepare | `species-check-relayed-20261029` once the committed schemas match the parser; anything else refuses |
| `@@EXCLUSIONS_COUNT@@` | `cleaned_rerun.py` | prepare | how many MTB genomes it fails (at least the two) |
| `@@RN_PREFIX_GENOMES_PATH@@`, `@@RN_PREFIX_GENOMES_SHA256@@` | `cleaned_rerun.py` | prepare | revision N's MTB 2k selection |
| `@@RN_PREFIX_SPECIES_JSON_PATH@@`, `@@RN_PREFIX_SPECIES_JSON_SHA256@@` | `cleaned_rerun.py` | prepare | its `curated/species.json` |
| `@@RN_FULL_GENOMES_PATH@@`, `@@RN_FULL_GENOMES_SHA256@@` | `cleaned_rerun.py` | prepare | revision N's cap-4000 selection |
| `@@RN_FULL_SPECIES_JSON_PATH@@`, `@@RN_FULL_SPECIES_JSON_SHA256@@` | `cleaned_rerun.py` | prepare | its `curated/species.json` |
| `@@RN_INPUT_ROOT@@` | `cleaned_rerun.py` | prepare | the data root of catalog, curated, records and annotations |
| `@@RN_INPUT_DATASETS_SHA256@@`, `@@RN_INPUT_PHENOTYPES_SHA256@@`, `@@RN_INPUT_CATALOG_VERSION_SHA256@@`, `@@RN_INPUT_CATALOG_CDS_SHA256@@`, `@@RN_INPUT_CATALOG_REFERENCE_SHA256@@` | `cleaned_rerun.py` | prepare | SHA-256 of its `curated/datasets.json`, `curated/phenotypes.jsonl`, `catalog/version.txt`, `catalog/AMR_CDS.fa` and `catalog/ReferenceGeneCatalog.txt` |
| `@@RN_INPUT_MANIFEST_SHA256@@` | `cleaned_rerun.py` | prepare | SHA-256 of the manifest of every input-root file the cohort roots can read, from `--derive-cleaned` |
| `@@RN_FASTA_MANIFEST_PATH@@`, `@@RN_FASTA_MANIFEST_SHA256@@` | `cleaned_rerun.py` | prepare | the FASTA manifest of the replacements, in the cohort evidence's columns |
| `@@RN_FASTA_STORE@@` | `cleaned_rerun.py` | prepare | the directory of those FASTA |
| `@@CLEANED_PREFIX_SHA256@@`, `@@CLEANED_PREFIX_COUNT@@`, and the same for `ADDED`, `FULL`, `MATCHED2K` | `cleaned_rerun.py` | prepare | the four cleaned lists, from `--derive-cleaned` |
| `@@RN_REFERENCE_ROOT@@` | `cleaned_rerun.py` | reference | revision N's trained MTB 2k candidate root |

## 3. Root setup (future, after gates 1–4)

A fresh run root and a fresh checkout: nothing of the current run's `$RUN` is reused or written.

```bash
export WS=/root/grm-research/grm-mtb-4k-cleaned-src    # a reviewed, clean checkout of the merge of revision N's commit with this branch
export RUN=/root/grm-research/grm-mtb-4k-cleaned-r1    # new; never the current run's root
export SLUG=mycobacterium_tuberculosis
export P=$WS/scripts/mtb-4k
export PROFILE=rN-cleaned                              # every report_contract CLI passes --profile "$PROFILE"
# The wiring lane's process-group guard (31 GiB disk / 48 GiB RAM floors), used read-only.
export GUARD=/root/grm-research/grm-mtb-4k-wiring-kimi-20261028/tools/resource_guard.py

set -o pipefail                                        # a failed sha256sum stops a `| tee` line too
test ! -e "$RUN" || exit 1
python3 -I "$P/revision_n.py" status || exit 1                     # every revision-N value filled and well formed
python3 -I "$P/cleaned_rerun.py" status --group prepare || exit 1   # every cohort input filled and well formed
python3 -I "$P/revision_n.py" verify-checkout "$WS" || exit 1       # clean, descends from RN_CHECKER_COMMIT, holds the pinned checker blobs
REVISION_N=$(python3 -I "$P/revision_n.py" get REVISION_N) || exit 1
export REVISION_N
python3 -I -c 'import os; p=os.statvfs("/root"); d=p.f_bavail*p.f_frsize; m={l.split(":",1)[0]:l.split(":",1)[1].strip() for l in open("/proc/meminfo")}; r=int(m["MemAvailable"].split()[0])*1024; print(d,r); assert d >= 45*1024**3 and r >= 64*1024**3' || exit 1
python3 -I /root/grm-research/grm-research-mtb-4k-cohort-r2/scripts/verify_evidence.py --root /root/grm-research/grm-research-mtb-4k-cohort-r2 || exit 1
python3 -I "$P/prepare_execution_roots.py" "$RUN" --cohorts cleaned || exit 1   # derives the four lists and requires the CLEANED_* pins

(cd "$WS" && CARGO_TARGET_DIR="$RUN/target" CARGO_BUILD_JOBS=16 RAYON_NUM_THREADS=16 cargo build --release --locked --offline) || exit 1
export GRM="$RUN/target/release/grm"
sha256sum "$GRM" | tee "$RUN/logs/grm-binary.sha256" || exit 1
# The committed foldgen; a revision-N merge that changes its tree stops here for review (I4).
test "$(git -C "$WS" rev-parse HEAD:scripts/metric-audit/foldgen)" = dc8317faf778fd96f496cc0226515ca7b19d688d || exit 1
(cd "$WS/scripts/metric-audit/foldgen" && CARGO_TARGET_DIR="$RUN/target-foldgen" cargo build --release --locked --offline) || exit 1
export FOLDGEN="$RUN/target-foldgen/release/foldgen"
sha256sum "$FOLDGEN" | tee "$RUN/logs/foldgen-binary.sha256" || exit 1
export FOLDGEN_SHA256=$(cut -d' ' -f1 "$RUN/logs/foldgen-binary.sha256")
export SOURCE=$(git -C "$WS" rev-parse HEAD)   # the one --expected-source-commit for every check
```

`prepare_execution_roots.py --cohorts cleaned` refuses before it creates anything when:

- a value it reads is unfilled;
- an input differs from its bound SHA-256;
- the cleaned 2k is not the first rows of the cleaned 4k;
- a cohort genome has no verdict in the species check, fails it, or is on its unchecked list;
- a selection's `species.json` reports `selected_unchecked` other than 0, or a length other than
  its list's;
- the input root lacks real `catalog/` or `curated/` directories or a BV-BRC record file;
- a cohort genome has no usable cached annotation (so `grm learn` would fetch it);
- the species check does not fail both ruled genomes, or fails other than `EXCLUSIONS_COUNT` MTB
  genomes;
- the replacement FASTA manifest covers other than exactly the replacements, or a FASTA differs
  from its manifest;
- a cleaned list differs from its `CLEANED_*` pin;
- 45 GiB would not stay free after the annotation copies.

It then makes the four training roots, and no pilot roots. Each root links its FASTA, links
`catalog` and `curated` to `RN_INPUT_ROOT`, and has a lane-local `bvbrc` as in revision 5.
`preparation.json` records the exclusions, the borderline genomes per cohort, the replacements,
and how the prefix and full lists differ from the frozen rows. `frozen-inputs.json` binds the
input root, every FASTA, every cohort list and every cleaned input. `--verify-frozen` re-hashes
all of them, and `--check-lane` reads the input root and FASTA stores from that record. The
lane checks of revision 5 apply unchanged to these roots:

```bash
python3 -I "$P/prepare_execution_roots.py" --check-lane "$RUN/cohorts/<cohort>/data"   # any problem is a stop
python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN"                      # any change is a stop
```

Revision 5's two read-only checks against the revision-7 candidate move to section 7b, against
revision N's candidate. They checked that the execution build reproduces real reports'
fingerprints and folds. For this rerun the reports that matter are revision N's, which the prefix
is diffed against. The foldgen tree is the pinned `dc8317f`, the same as the one that rebuilt the
revision-7 folds, so repeating the revision-7 checks adds nothing. Revision N's `audit.py` would
also have to keep auditing revision 7, which is not this plan's to require.

**Expected time and disk for sections 3–6, from the current run's measured stages (its section 5
receipt).** Preparation hashes the 4,000 frozen FASTA, the replacements and the input root; its
time was not measured here. Section 5 took 57 min 10 s for four matrices, their validations, containment and the
duplicate screen. Section 6 took 64 s. Allow about 1 h 15 min. The peak disk growth was
6,226,591,744 B (5.80 GiB), during the full matrix. The kept full matrix was 367,757,922 B, and
the annotation copies 366,450,430 B. The cleaned cohorts drop the two genomes' roughly 9.2 M rare
k-mers each, so expect slightly less. Training (section 7) took 16 min 31 s for prefix (22
reports) and 13 min 41 s for added (14). Full had 9 of 22 reports after 17 min, so allow about
2 h 15 min for section 7 with its checks. The 45 GiB gate is re-measured before every stage, as in
revision 5.

## 4. The small frozen biological pilot (gate 2)

**Not repeated in revision 6.** The admin ruled on 2026-10-29 that the cleaned run makes fresh roots,
matrices and structures without repeating the pilot. The pilot passed on 2026-10-29 (attempt 3)
under revision 5. The review of the revision-N merge must confirm that revision N leaves matrix
building unchanged. **Rule:** no change, between revision 5's build source (`c70d678`) and the merge,
to the matrix code or anything it transitively uses, to the `matrix` and `count` wiring in
`src/main.rs`, or to `src/lib.rs`, `Cargo.toml`, `Cargo.lock` or `.cargo/config.toml`. The last sets
`target-cpu=native`, so it changes the machine code of every matrix module. Otherwise the pilot is
repeated before section 5. Binary hashes are host-specific for the same reason: the training binary
is built and its hash recorded on the host that runs it (section 3), never copied from another
host.

The current expansion of the rule is the transitive `crate::` closure of `src/matrix.rs`,
`src/kset.rs` and `src/species.rs`, generated from the code at `c70d678`:
- `src/bits.rs`, `src/bvbrc.rs`, `src/catalog.rs`, `src/cluster.rs`, `src/confine.rs`, `src/dna.rs`,
  `src/drugs.rs`, `src/fasta.rs`, `src/fnv.rs`, `src/kset.rs`, `src/layout.rs`, `src/matrix.rs`,
  `src/probes.rs`, `src/species.rs` and `src/species/`;
- plus `src/main.rs`'s wiring, `src/lib.rs`, `Cargo.toml` and `Cargo.lock`.

If the merge adds a `use` to any of these, the rule covers the new module too. A cleaned run root has no pilot cohorts.
The text and commands below are revision 5's, kept for the record.

It is separate from the 2,000-genome prefix stage and uses its own two roots,
on rows 1–256 of the frozen order (SHA-256 recorded in `preparation.json`).
Width 128 gives two blocks, so a fold, a resume and the partial merge are all
exercised on real MTB data.

The kill is no longer a timed poll and a `pkill -f` pattern (N8): the pattern
also matched the guard's and `/usr/bin/time`'s command lines and killed the
guard, and nothing checked that the run was interrupted. `pilot_interrupt.py`
takes the guard's PID (`$!`) as `--launcher` and binds only to the process
below it whose command line starts with exactly `$GRM --data $SP matrix`
(L1). A matching process outside the launch, such as a run left from an
earlier attempt, is never bound or signalled. It is a stop, because two runs
would share the root. The helper waits until `run.json` records block 0
folded (`folded_columns` = 128: its partial verified and its owned ksets
deleted), then stops the process. It scans again and requires the same single
process, then checks the stopped state and kills it. It requires the guard's
record to show a SIGKILLed run (rc 137, no guard reason) and the root to be
still in the stopped state.

It exits non-zero, and the pilot stops, when the run ends first, when the
point is not reached in time, when anything else matches, or when a check
fails. On every such stop it first stops the launch itself (freeze, kill
what runs below the guard, let the guard record that and exit), so the stop
never leaves the launched run going. A `--launcher` that is not a sibling job of the
helper is refused and nothing is signalled. So is a sibling whose command
line is not this run's guard: `python3 $GUARD` with `--resources` equal to
`--guard-resources`, and `$GRM --data $SP matrix` after its `--` (N1).

The guard starts its command in a new session, so a guard that dies before
the helper binds the run leaves `grm` running (N2). After any stop the
helper lists every `$GRM --data $SP matrix` still running and signals none.
After the helper and `wait`, step a runs `--list-runs`, which prints the same
list and exits 1 when it is not empty: any process left is a stop for the
operator, never a signal from the plan.

```bash
G="python3 $GUARD"
SP=$RUN/cohorts/pilot256/data; ST=$RUN/cohorts/pilot256-stock/data
for d in "$SP" "$ST"; do python3 -I "$P/prepare_execution_roots.py" --check-lane "$d" || exit 1; done
# a) sparse, interrupted once block 0 is folded and block 1 is under way
$G --log $RUN/logs/pilot-sparse-1.log --resources $RUN/logs/pilot-sparse-1.json -- \
  /usr/bin/time -v "$GRM" --data "$SP" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128 &
LAUNCHER=$!
python3 -I "$P/pilot_interrupt.py" --launcher "$LAUNCHER" --guard-resources $RUN/logs/pilot-sparse-1.json \
  --data "$SP" --grm "$GRM" --folded 128 --receipt $RUN/audit/pilot-interrupt.json
HELPER_RC=$?
wait "$LAUNCHER"; WAIT_RC=$?
python3 -I "$P/pilot_interrupt.py" --list-runs --data "$SP" --grm "$GRM" || exit 1   # lists, never signals (N2)
test "$HELPER_RC" -eq 0 || exit 1
test "$WAIT_RC" -eq 137 || exit 1   # the guard exits with time's status for the killed run: 128 + 9
ls $SP/species/$SLUG/ksets | wc -l          # record: block 0's owned ksets are gone, block 1's in progress
# a2) the invalidation copies, taken now, of the interrupted root and before the resume
mkdir $RUN/pilot-copies
for c in width forge link-matrix link-kset; do cp -a "$RUN/cohorts/pilot256" "$RUN/pilot-copies/$c"; done
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
# d) the invalidation checks, on the step-a2 copies only, never on the live root (L5)
$G --log $RUN/logs/pilot-invalidation.log --resources $RUN/logs/pilot-invalidation.json -- \
  python3 -I "$P/pilot_invalidation.py" --copies $RUN/pilot-copies --grm "$GRM" --width 128 \
  --receipt $RUN/audit/pilot-invalidation.json || exit 1
python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN"
```

`pilot_invalidation.py` runs the invalidation checks on the four copies from
step a2. Each copy is the interrupted state, which the helper checks first:
block 0 folded (its owned ksets deleted), block 1 not finished. Each check
runs the resume command of step b on its copy, and it must be refused and
leave the copy's species directory exactly as it was (every file's SHA-256,
every link's target), with the partials in place.

1. A width change (`--block-width 64` on the `width` copy) is refused while
   the width-128 run's partials exist ("unfinished disk-sparse run").
2. On the `forge` copy, block 0's partial is edited together with its
   receipt: one column loses a k-mer (one bit of a row that holds that column
   and another is cleared), and the receipt's SHA-256 and byte count are
   recomputed. The resume must fail `verify_partial` with "does not hold
   that genome's kset". Revision 2 forged the partial inside the kill window
   it chose, while block 0 was still active and its ksets present; there GRM
   correctly rebuilds an unverified partial from the ksets and finishes
   byte-equal to stock (N2). That is the expected outcome for an active
   block, so this check forges a folded one.
3. On the `link-matrix` copy, a `matrix` directory link to a directory beside
   the copy, and on the `link-kset` copy, a `ksets/<id>.kset.part` link (a
   block-1 genome with no kset yet) to a file beside it, are each refused
   before any work ("refusing symlinked directory"; a symlink refusal), and
   their targets keep their bytes.

Acceptance, all required:

- `pilot_interrupt.py` exits 0: one process matched, below the launcher and
  the same one at the stop, and nothing outside it. The stopped state had
  `folded_columns` = 128, no saved matrix, no `matrix.json`, and block 0's
  partial and receipt. The process was gone after the kill, the guard's
  record for step a shows the killed run (rc 137, no guard reason), `wait`
  returned 137, the root was still in the stopped state, and `--list-runs`
  found no `$GRM --data $SP matrix` still running.
- Sparse and stock matrix files and `structure.json` are byte-identical.
- The resume completes. No owned kset is deleted before its block's verified
  partial, no `matrix-partials`, owned ksets or `pending` are left, and the
  frozen inputs are unchanged.
- `validate_matrix.py` passes, including 16 recounts.
- `pilot_invalidation.py` exits 0: the three invalidation checks (four
  runs) refuse as stated and change nothing.
- Measured figures are recorded: per-genome kset bytes against the
  preflight's upper bound (12 KiB + 8 bytes per FASTA base), per-block
  partial bytes, matrix bytes, peak RSS, minimum free disk and RAM, and wall
  time. The 4k projection from them (all partials + one block's ksets +
  matrix) must fit in 15 GiB with 20% margin, or the plan stops for a new
  resource decision. It never lowers the reserve.

## 5. Matrix-only stages (future, separately authorized)

The sparse matrix runs at width 128, with 16 Rayon threads and under the guard. **No
`--allow-missing-catalog`**: the catalog must load, or the run stops before counting. The
per-block preflight bound and its arithmetic are revision 5's, and the pilot measured them.

```bash
for cohort in prefix added full matched2k; do
  D=$RUN/cohorts/$cohort/data
  python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || exit 1
  python3 "$GUARD" --log $RUN/logs/$cohort-matrix.log --resources $RUN/logs/$cohort-matrix-resources.json -- \
    /usr/bin/time -v "$GRM" --data "$D" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128 || exit 1
  python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
    --genomes $D/species/$SLUG/genomes.txt --fasta-dir $D/fasta --recount 16 --require-recount \
    --out $RUN/audit/$cohort-matrix-validation.json || exit 1
  python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || exit 1
done
# Containment (N3): each 2k cohort against the cleaned full matrix, after all four exist.
F=$RUN/cohorts/full/data/species/$SLUG/matrix
for cohort in prefix added matched2k; do
  D=$RUN/cohorts/$cohort/data
  python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
    --genomes $D/species/$SLUG/genomes.txt --superset "$F" --out $RUN/audit/$cohort-in-full.json || exit 1
done
```

Run the cohorts one at a time, prefix first. Any failure stops the later ones. Re-check the
resource gate before each.

Revision 5's byte comparison of the prefix matrix with a frozen matrix moves to section 7b. It
now compares against revision N's MTB 2k candidate, which may not be staged yet, and it never
blocks training. Until it passes, the receipt records `PENDING REFERENCE`. Without it, section 5
proves each cleaned matrix by the exact region, the 16 recounts and the recount bounds. It proves
each 2k cohort against the cleaned full matrix by containment. Revision 5's chain to stock-built
content ran through the frozen prefix, and here it runs through section 7b's `cmp`.

`validate_matrix.py` checks the binding to `counted.json` and `genomes.txt`, the structure, every
column against its genome's bottom sketch, and 16 seeded recounts from FASTA, independently of
GRM. Its three omission checks (exact region, recount bounds, containment) are revision 5's
unchanged. Their blind spot is in section 9.

After each cohort, check the receipt: no `pending`, no `matrix-partials`, no owned ksets, exact
order and digests, and the resource minima.

### Duplicate screen (after the prefix and added matrices, before training)

```bash
python3 -I "$P/screen_duplicates.py" --prefix-root $RUN/cohorts/prefix/data --added-root $RUN/cohorts/added/data \
  --fasta-manifest $RUN/fasta-manifest.tsv \
  --genomes-jsonl $RUN/cohorts/full/data/bvbrc/genomes.jsonl --out $RUN/audit/duplicate-screen.json
```

The manifest is the one preparation wrote over every cohort genome, replacements included. The
screen and its stop (`cross_cohort_exact > 0`) are revision 5's.

## 6. Structure (separate go-ahead after matrix review)

```bash
for cohort in prefix added full matched2k; do
  python3 "$GUARD" --log $RUN/logs/$cohort-structure.log --resources $RUN/logs/$cohort-structure-resources.json -- \
    "$GRM" --data "$RUN/cohorts/$cohort/data" structure --species $SLUG || exit 1
done
```

## 7. Training and reports: the revision-N gate

Training, reports and the comparison run only when **all** of these hold:

1. **Revision N's checker is reviewed and in the binary.** `$WS` is a reviewed merge of revision
   N's commit with this branch's two commits, trailer-free. `revision_n.py verify-checkout` passed
   in section 3: `HEAD` descends from `RN_CHECKER_COMMIT`, and the four checker files hold the
   pinned SHA-256 in `$WS` and at that commit. `cargo test` and the wiring adversarial checks pass
   on the merge, and so does the validator suite under `MTB4K_TEST_PROFILE=rN`
   (`tests/profile_boot.py`). For reference, revision 9's facts as relayed on 2026-10-29 are below.
   Revision N's values will differ, at least in the revision and the checker pins:
   - `b2fd433` (`EVALUATION_REVISION = 9`);
   - `check-results.py` `b7e5a8a0…`, `audit.py` `c1af2ead…`, `summarize.py` `fd5fdfab…`,
     `result_scope.py` `3bcb2db3…`;
   - no report or checker field added, renamed or redefined;
   - revision 9 judged under the revision-8 lasso rules.

   The filling commit must take its values from that commit, not from this text. A moved pin, an
   added or renamed field, or another convergence rule is a stop until `report_contract.py` is
   reviewed again; that stop is intended.
2. **Every cohort passes `run_revision8_checks.py --profile rN-cleaned`.** It is not a parallel
   checker. It runs revision 5's chain with revision N's numbers:
   - the run manifest and every report pass `report_contract.py` under the `rN-cleaned` rules,
     with each report's input fingerprint recomputed from the root;
   - the in-tree `scripts/check-results.py`;
   - `scripts/metric-audit/audit.py --revision $REVISION_N`, with this cohort's backup, log and
     launch, and the foldgen built in section 3, which rebuilds every saved test_fold;
   - `summarize.py`.

   Only when all of these pass does it write the receipt that `compare_cohorts.py` requires. The
   file and receipt keep their revision-8 names, and `evaluation_revision` in the receipt says N.
3. **An invalid report is a stop.** Record the dataset and its reasons, and stop that cohort's
   training for parent review. Any solver status field the checker does not recognize is a stop.

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
    --log $RUN/logs/$cohort-learn.log --profile "$PROFILE" || exit 1
  python3 -I "$P/run_revision8_checks.py" $RUN/audit/$cohort-audit.json --root "$C" --run-manifest "$C/validation-manifest.json" \
    --expected-source-commit "$SOURCE" --binary "$GRM" --backup "$RUN/backups/$cohort" --log $RUN/logs/$cohort-learn.log \
    --launch "$LAUNCH" --foldgen "$FOLDGEN" --expected-foldgen-sha256 "$FOLDGEN_SHA256" --checker-out $RUN/audit/$cohort-check-results.json \
    --summary-out $RUN/audit/$cohort-summary.json --receipt $RUN/audit/$cohort-revision8-receipt.json --profile "$PROFILE" || exit 1
done
```

`write_run_manifest.py` binds each cohort's inputs, reports and class-gate receipts by SHA-256,
after recomputing every input fingerprint and every dataset's class counts from the root.
Expected report counts are 22 for prefix and full and 14 for added (frozen with a check, section
2), and whatever the class gate gives for matched2k. Build dashboard data only in the new roots,
each after `--check-lane`. Do not publish, deploy, replace candidate reports or touch the live
site.

### 7b. Revision N's MTB 2k candidate (pending reference; required for the final result)

These checks need revision N's trained MTB 2k candidate (`RN_REFERENCE_ROOT`). Run them as soon as
it is staged on Vast, before or after sections 5–7 have finished; they never block training. Until
they pass, every receipt and summary of this run carries `PENDING REFERENCE`, and no final result
is stated. A failure is a stop for review, never something to skip.

```bash
set -o pipefail
REF=$(python3 -I "$P/cleaned_rerun.py" get RN_REFERENCE_ROOT) || exit 1     # unfilled: PENDING REFERENCE
INPUT_ROOT=$(python3 -I "$P/cleaned_rerun.py" get RN_INPUT_ROOT) || exit 1
R=$REF/data/species/$SLUG
sha256sum $R/genomes.txt $R/counted.json $R/structure.json $R/matrix/* $R/results/*.json | tee $RUN/audit/reference-sha256.txt || exit 1
# The reference must be revision N's candidate for this very prefix, and its matrix the prefix's, byte for byte.
cmp $RUN/cohorts/prefix/data/species/$SLUG/genomes.txt $R/genomes.txt || exit 1
for f in matrix.json patterns.bin support.bin kmer_offsets.bin kmers.bin; do
  cmp $RUN/cohorts/prefix/data/species/$SLUG/matrix/$f $R/matrix/$f || exit 1
done
# fingerprint.py reproduces each reference report's input_fingerprint from revision N's matrix and the rerun's
# inputs (RN_INPUT_ROOT's curated and catalog), in a link-only view outside both roots.
python3 -I "$P/fingerprint.py" "$RUN/rN-view/data" --make-view --matrix-from "$REF/data" --inputs-from "$INPUT_ROOT" \
  --reports-from "$REF/data" --check --expect-reports 22 --out "$RUN/audit/rN-fingerprints.json" || exit 1
# The metric audit rebuilds every saved test_fold of the reference with this foldgen.
python3 -I "$WS/scripts/metric-audit/audit.py" "$RUN/audit/rN-fold-rebuild.json" $SLUG --revision "$REVISION_N" \
  --root "$REF" --reshaped --foldgen "$FOLDGEN" || exit 1
python3 -I -c 'import json,sys; r=json.load(open(sys.argv[1])); bad=[i for x in r["records"] for i in x.get("issues",[]) if "test_fold" in i["where"]]; print(len(r["records"]),len(bad)); assert len(r["records"])==22 and not bad' "$RUN/audit/rN-fold-rebuild.json" || exit 1
# The prefix against the reference, every field (diff_reference_prefix.py).
python3 -I "$P/diff_reference_prefix.py" --prefix $RUN/cohorts/prefix/data --reference "$REF/data" \
  --out $RUN/audit/prefix-vs-reference.json || exit 1
```

**What must match, and why.** The two sides share the same genome list (checked), a byte-identical
matrix (checked), the same inputs (the fingerprint view), and the same revision, solver,
`EvalConfig` and seeds. So `diff_reference_prefix.py` requires every field to be identical, lasso
included. Only the wall-clock `seconds`, at the top level and on each learner row, is left out.
The revision-9 note of 2026-10-29 says only lasso-derived parts move between revisions:
`features`, `best_learner` and `evaluation.revision`. Within one revision none of them may differ,
so they are compared too. A degenerate lasso optimum (exact complement columns on a fold's fitting
samples) can move held-out scores by up to about 0.1 at equal objective, but only across solvers.
Within one Rust source, a lasso-only difference here is a stop for review, not noise. The
reviewer then checks whether the reference was built from another Rust tree or Rayon thread
count before calling it anything else.

## 8. Matched analysis (after section 7)

```bash
python3 -I "$P/compare_cohorts.py" --profile "$PROFILE" \
  --cohort prefix=$RUN/cohorts/prefix --cohort added=$RUN/cohorts/added \
  --cohort full=$RUN/cohorts/full --cohort matched2k=$RUN/cohorts/matched2k \
  --checker-receipt prefix=$RUN/audit/prefix-revision8-receipt.json --checker-receipt added=$RUN/audit/added-revision8-receipt.json \
  --checker-receipt full=$RUN/audit/full-revision8-receipt.json --checker-receipt matched2k=$RUN/audit/matched2k-revision8-receipt.json \
  --duplicates $RUN/audit/duplicate-screen.json --expected-source-commit "$SOURCE" --expected-foldgen-sha256 "$FOLDGEN_SHA256" \
  --out $RUN/audit/compare-cohorts.json
```

The comparison is within revision N and within this run. Prefix, added, full and matched2k come
from one selection, one binary and one checker, so no solver change sits between them. Exit codes,
stops, decision D1 and the prespecified rules are revision 5's, with two changes. First,
`matched2k` is `matched_selection` of the cleaned full cohort at the prefix's size (F6). Second,
the eight not-interpretable drugs and the added cohort's 14 eligible drugs are checked against the
same frozen sets. The counts quoted in revision 5's list (levofloxacin 774 resistant and 1
susceptible; para-aminosalicylic acid 5 resistant) are the frozen added rows'. The cleaned added
rows give their own counts, which the comparator derives and records. If cleaning moves a drug in
or out of either set, that is a stop for the parent, not a new rule.

**With/without contaminants, and what can be compared (confounding note).**

- **(b) Expected: this cleaned run on N = 10, after the uncleaned revision-9 rerun.** A
  with/without comparison would then be possible, and **cohort-level**; by the ruling below it is
  described only qualitatively, and nothing is measured. The uncleaned cohorts are
  the frozen 4k with the two genomes. The cleaned ones are a fresh selection: two genomes out,
  replacements in, and order shifts near score ties, in the prefix and in full.
  - **Solver confound.** If revision 10's solver is unchanged from 9, there is none: the same
    solver, `EvalConfig` and seeds on both sides. Check it before any comparison: `src/learn/`
    and the lasso path must be unchanged between `b2fd433` and `RN_CHECKER_COMMIT`. A
    `check-results.py` change that only admits revision 10 changes no judgement. If the solver
    did change, the confound returns for lasso only. The revision-9 note says no non-lasso bytes
    change from code, and a degenerate optimum can move lasso held-out scores by about 0.1 at
    equal objective. So bound it by reading the four non-lasso learners (scm, cart, rf, gbt) and
    the catalog row apart from lasso.
  - **Duplicate isolation is a second revision-10 change.** If it changes how duplicates are
    grouped into folds or lineages, the two sides also differ in fold construction, not only in
    cohort. That is a method change, so name it in any statement. Bound it by comparing the
    random split (seeded shuffle, untouched if isolation acts only on grouping) with the lineage
    split, once its effect is documented.
  - **What it can support.** Per drug, how the results of the uncleaned and the cleaned selections
    differ, as whole cohorts: lineage-split medians and per-learner rows, paired on the genomes
    common to both, with section 8's descriptive McNemar counts.
  - **What it cannot support.**
    - Any attribution to the two genomes alone. Replacements and order shifts change the cohort
      too, and every fold moves: grm's folds are a seeded shuffle over each cohort's own samples.
      Decision D1 measured 28–62% of common genomes changing fold, and about 80% is expected.
    - Any claim about the six borderline genomes, which both sides keep.
    - Anything external, prospective or clinical.
  - **The feature-space effect.** It is real but not a label effect. The matrix review found about
    1.58 M MTB-core k-mers kept only because the two genomes lack them: core counts were 1,613 in
    prefix against 1,584,185 in added. Removing the two drops those k-mers from the kept set, and
    with them the screening of every learner. A difference no larger than the run-to-run fold
    noise is not attributable to anything.
  - **Ruling (the GRM admin thread, 2026-10-30, after `PLAN_R6=READY`).**
    1. **No uncleaned control run on revision 10.** No quantitative with/without claim will be
       made, so there is nothing to deconfound. The difference is 2 genomes of 2,000 or 4,000, at
       bottom-sketch Jaccard 0.001, and the replacements reshuffle the cohort.
    2. **No cross-run comparator.** The with/without difference is described only qualitatively.
       The write-up says plainly that the difference is cohort-level, confounded with revision
       10's duplicate isolation and any solver change, and not measured. It gives no number and no
       direction for the effect of the two genomes. The bounding readings and the "can support"
       list above describe what a measured comparison could have shown; none is run.
    3. **What counts as the result.** The uncleaned revision-9 run is the first complete MTB 4k
       result, and the like-for-like comparison with revision 7. That comparison is its section
       7 diff (`diff_revision7_prefix.py --profile rN`, its prefix against the revision-7
       candidate); its section 8 compares within revision 9. This cleaned revision-10 run's own
       section 8, prefix against full inside the revision, becomes the canonical result once it
       exists.
- **(a) Alternative: one cleaned run that is also the cohorts' rerun.** This applies only if the
  species check lands in the revision that reruns the 4k cohorts. Then there is no with/without
  comparison: the admin dropped it as low value, since the two genomes have bottom-sketch Jaccard
  0.001 against the cohort (the matrix review). Section 8 compares prefix, added, full and matched2k
  within that revision, and nothing is confounded by a solver change.

Revision 5's text, which stands for revision 6:

Exit 0 means complete. Exit 1 means output was written but the comparison is
invalid (source, binary, foldgen or EvalConfig differing across cohorts, a
coverage or label mismatch, cohorts not disjoint). Exit 2 means an input or
integrity stop with no output. The stops are:

- a manifest, hash or fingerprint binding, or any invalid report (L3);
- a receipt that does not prove its cohort passed (F5). The receipt's own
  `valid` flag is not enough: the checker, audit and summary outputs it
  names must still hash as it says and, read again, record a pass for this
  run manifest, source, binary and foldgen (L2). Its foldgen must be the
  build recorded in section 3 (I4);
- the matched2k list not being `matched_selection` of the full cohort (F6);
- a class-gate receipt or report count that the full cohort's labels on that
  cohort's genomes do not reproduce (F3);
- an unexpected not-interpretable set, or an unaccepted exact duplicate;
- validators that cannot load (a moved pin).

Exit 3 means output was written and valid, but a primary split was
unavailable; stop for review before any before/after statement.

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
- **Any invalid report stops the comparison (L3).** Revision 3 had a
  per-drug exclusion rule: such a drug was dropped from every comparison
  involving its cohort, and the run exited 3. No reviewed step could produce
  that state. Revision-8 `grm learn` stops a species at its first invalid
  report, and `run_revision8_checks.py` gives a cohort with any invalid report
  no receipt, because the reviewed in-tree checker passes or fails a cohort
  as a whole. Making the exclusion reachable would need a per-drug checker,
  which section 7 forbids. So an invalid report, or a missing eligible one,
  in any cohort stops with exit 2 and no output, and the parent decides.
  Nothing is imputed. Ineligibility under the class gate is not an invalid
  report; it is listed under `not_eligible`.
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

- **Pending reference.** Section 7b needs revision N's trained MTB 2k candidate (`RN_REFERENCE_ROOT`),
  which may not be staged when this run trains. Until 7b passes, two links are unproven:
  - the chain from the cleaned matrices to stock-built content (revision 5 had it through the
    frozen prefix);
  - the claim that this run reproduces revision N's own candidate field for field.

  Results carry `PENDING REFERENCE` until then and are not final. A 7b failure is a stop for
  review.
- **Borderline genomes stay in.** Six MTB genomes are borderline: MTB by sequence, but GC 62.7–64.6%
  and 4.52–4.83 Mb, so they probably carry foreign contigs. They are `1773.18404`, `1773.3676`,
  `1773.18177`, `1773.40722`, `1773.18213` and `1773.5309` (in added). They are recorded per cohort
  in `preparation.json` and never excluded. A contig-level check comes in revision 10. It is out of
  scope here, and a result touched by them is not qualified beyond that record.
- **The species check's schemas are provisional.** On 2026-10-29 the relayed `species_check.json`
  (`5ba9d78e…`) and the scratch copy read on the laptop (`a11832f4…`, which adds `pairs`,
  `config.exempt_pairs`, `held_pairs` and a per-genome `inside_species`) already differed.
  `prepare_execution_roots.py` parses only these fields:
  - the MTB `species` entry (`slug`, `species`, `unchecked`, `failed`, `borderline`);
  - each genome row's `genome_id`, `species` and `verdict`;
  - the selection's `curated/species.json` list (`slug`, `selected_genomes`,
    `selected_unchecked`).

  It ignores everything else. `SPECIES_CHECK_FORMAT` is filled only when the committed code
  writes those fields so.
- **The species check is not re-derived here.** The plan trusts the bound post-fetch output's
  verdicts. Its method (GC plus a bottom-1000 MinHash at k=21 against a per-species medoid, fail
  at d > 0.13 or inside another species' 0.05 boundary) is the species-check lane's to review.
- **Fresh-selection drift.** Replacements and order shifts near score ties are accepted by
  ruling. The cleaned 2k must still be the cleaned 4k's first rows. The plan does not measure how
  either list differs from the frozen ones beyond what `preparation.json` records: who left, who
  joined, and whether the common genomes keep their order.
- Unchanged from revision 5, and still not covered:
  - real power loss and filesystem-specific ordering on the Vast overlay;
  - the data root and its ancestors are trusted, and bind mounts are not detected (N7);
  - `grm fetch` is not confined, and the plan does not run it;
  - the containment blind spot of section 5 (k-mers kept in full and in no 2k cohort are checked
    only by the exact region and the recount bounds);
  - `--check-lane` does not flag hardlinks;
  - the deferred items N6, N8, N9 (only partly closed here: `revision_n.py verify-checkout`
    requires a clean checkout of `$WS`), N11, N12 and N13 (`audit.py` still stops with a
    traceback on missing inputs).
- This plan does not run `diff_revision7_prefix.py`: its prefix is a fresh selection, not the
  frozen revision-7 prefix. Since `8bb1a9b` that script takes `--profile rN` for a revision-9 rerun
  of the uncleaned cohorts (revision 5's section 7). Its field rules are unchanged.
