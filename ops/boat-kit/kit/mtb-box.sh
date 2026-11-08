#!/usr/bin/env bash
# The MTB cleaned 4k run (plan revision 6, scripts/mtb-4k/4k-execution-plan.md at bc7b742) inside one Boat box,
# called by r11-box.sh _run when KIND=mtb. SPECIES lists this box's sections in order:
#   prepare            section 3: RUN_4K=YES, the placeholder files filled, the resource gate, the cohort evidence
#                      (verify_evidence.py from the reviewed checkout, K1), prepare_execution_roots.py --cohorts
#                      cleaned. verify-checkout and the foldgen tree pin need git, which the box's source doesn't have:
#                      drive.sh runs them on the laptop against SRC before the box (K2, main's ruling of 5 Nov)
#   restore            X2-X4: the cohort counted.json, structure.json, matrix/ and results/ of earlier boxes, and
#                      their run receipts (audit/*: revision-8 receipts, check-results, summaries, the duplicate
#                      screen; cohorts/<c>/validation-manifest.json), pushed as one push-data set and checked
#                      against its manifest (MTB_RESTORE_SUMS) first
#   matrix:<cohort>    section 5, one cohort; containment once all four exist. matrix:prefix then stops the box
#                      unless the prefix genomes.txt and its five matrix files equal revision N's g5 root's (R7's,
#                      at MTB_REFERENCE_ROOT) byte for byte (Q1, main's ruling of 5 Nov; receipt
#                      audit/prefix-matrix-vs-reference.json)
#   dupscreen          section 5's duplicate screen
#   structure          section 6, all four cohorts
#   learn:<cohort>     section 7, one cohort
#   compare            section 8
#   reference          section 7b against RN_REFERENCE_ROOT (revision N's g5 data root), with the fold rebuild (K3)
# Every section starts at the resource gate: free disk >= 45 GiB (the plan's, unchanged) and effective memory >=
# MTB_RAM_NEED_GIB (7, PF-D7, replacing the plan's 64 GiB); the guard pauses below RAM_PAUSE_GIB (4, replacing 48) and
# below DISK_PAUSE_GIB (31, unchanged). Paths: WS=/home/user/grm (the merged, reviewed checkout), RUN=W/mtb.
# MTB_SYNTH=1 is the dry run: the same grm commands on synthetic cohort roots built here from the synthetic g5 root,
# with every plan validator (prepare/check-lane/verify-frozen, validate_matrix, screen_duplicates, the contract
# CLIs, compare_cohorts, fingerprint, diff_reference_prefix) replaced by a SKIPPED-SYNTH line in logs/sections.log.
# The dry run still runs the Q1 prefix gate and the K3 fold rebuild for real, against the synthetic g5 root, and
# builds the sparse matrix when the source has --disk-sparse.
# usage: mtb-box.sh SECTION...     (W, DATA, GRM, FOLDGEN, the stage.env values from r11-box.sh)
set -u -o pipefail
K=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WS=${WS:-/home/user/grm} RUN=$W/mtb SLUG=mycobacterium_tuberculosis
P=$WS/scripts/mtb-4k PROFILE=rN-cleaned
export DISK_PAUSE_GIB=${MTB_DISK_PAUSE_GIB:-31} RAM_PAUSE_GIB=${RAM_PAUSE_GIB:-4}
gib() { awk -v g="$1" 'BEGIN { printf "%d", g * 1073741824 }'; }
slog() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$W/logs/sections.log"; }
stop() { slog "STOP $*"; echo "STOP $*" >&2; exit 1; }
skip() { slog "SKIPPED-SYNTH $*"; }
real() { [ "${MTB_SYNTH:-0}" != 1 ]; }
gate() { # the resource gate before every section, re-measured each time
  local d m
  d=$(df -B1 --output=avail /home/user | tail -1 | tr -d ' ')
  m=$(($(awk '/MemAvailable/ {print $2}' /proc/meminfo) * 1024))
  slog "gate $1: free=$d (>= $(gib 45)) mem_available=$m (>= $(gib "$MTB_RAM_NEED_GIB"))"
  [ "$d" -ge "$(gib 45)" ] && [ "$m" -ge "$(gib "$MTB_RAM_NEED_GIB")" ] || stop "resource gate before $1"
}
cohorts="prefix added full matched2k"

sec_prepare() {
  [ "${RUN_4K:-NO}" = YES ] || stop "RUN_4K is not YES (G1)"
  test ! -e "$RUN" || stop "$RUN exists"
  gate prepare
  if real; then
    python3 -I "$P/revision_n.py" status || stop "revision_n.py unfilled"
    python3 -I "$P/cleaned_rerun.py" status --group prepare || stop "cleaned_rerun.py unfilled"
    [ "$(python3 -I "$P/revision_n.py" get REVISION_N)" = "$REVISION" ] || stop "REVISION_N is not $REVISION"
    python3 -I "$P/verify_evidence.py" --root "$MTB_EVIDENCE" > "$RUN.evidence.json" || stop "cohort evidence (verify_evidence.py)"
    python3 -I "$P/prepare_execution_roots.py" "$RUN" --cohorts cleaned || stop "prepare_execution_roots"
    mv "$RUN.evidence.json" "$RUN/evidence-verification.json"
  else
    # synthetic cohorts, by the plan's rules: full = the cap-120 selection, prefix = the cap-60 selection and its
    # first rows, added = the rest, matched2k = the prefix's size by smallest sha256 of the plan's salt, in full's order
    python3 -I - "$DATA" "$RUN" "$SLUG" << 'EOF' || stop "synthetic cohort roots"
import hashlib, os, shutil, sys
d, run, slug = sys.argv[1:4]
full = open(os.path.join(d, "bindings/full-genomes.txt")).read().split()
prefix = open(os.path.join(d, "species", slug, "genomes.txt")).read().split()
assert full[:len(prefix)] == prefix, "the 2k selection is not the first rows of the 4k"
added = full[len(prefix):]
key = lambda g: hashlib.sha256(("grm-mtb-4k-matched2k-20261028:" + g).encode()).hexdigest()
pick = set(sorted(full, key=key)[:len(prefix)])
lists = {"prefix": prefix, "added": added, "full": full, "matched2k": [g for g in full if g in pick]}
for c, ids in lists.items():
    r = os.path.join(run, "cohorts", c, "data")
    os.makedirs(os.path.join(r, "species", slug))
    os.makedirs(os.path.join(r, "fasta")); os.makedirs(os.path.join(r, "bvbrc/features"))
    for x in ("catalog", "curated"):
        os.symlink(os.path.join(d, x), os.path.join(r, x))
    for x in ("genome_amr.jsonl", "genomes.jsonl"):
        os.symlink(os.path.join(d, "bvbrc", x), os.path.join(r, "bvbrc", x))
    for g in ids:
        os.symlink(os.path.join(d, "fasta", g + ".fna.gz"), os.path.join(r, "fasta", g + ".fna.gz"))
        shutil.copy2(os.path.join(d, "bvbrc/features", g + ".json"), os.path.join(r, "bvbrc/features"))
    open(os.path.join(r, "species", slug, "genomes.txt"), "w").write("\n".join(ids) + "\n")
    print(c, len(ids))
EOF
    skip "revision_n/cleaned_rerun status, verify_evidence, prepare_execution_roots (synthetic cohorts built instead)"
  fi
  mkdir -p "$RUN/logs" "$RUN/audit" "$RUN/backups"
  sha256sum "$GRM" "$FOLDGEN" > "$RUN/logs/binaries.sha256"
}
sec_restore() { # an earlier box's cohort counted/structure/matrix, checked against its receipt before use
  local src=/home/user/grm-resume
  (cd "$src" && sha256sum -c --strict --quiet "$MTB_RESTORE_SUMS") || stop "restore set differs from $MTB_RESTORE_SUMS"
  for c in $cohorts; do
    [ -d "$src/cohorts/$c" ] || continue
    cp -a "$src/cohorts/$c/data/species/$SLUG/." "$RUN/cohorts/$c/data/species/$SLUG/"
    slog "restored $c: $(ls "$RUN/cohorts/$c/data/species/$SLUG" | paste -sd' ')"
    [ -f "$src/cohorts/$c/validation-manifest.json" ] && cp -a "$src/cohorts/$c/validation-manifest.json" "$RUN/cohorts/$c/"
    # X4: compare_cohorts reads each cohort's report tree (the validation manifest hashes its index), which X2
    # and X3 built; it arrives under data/report/ and must land beside the cohort's species outputs
    if [ -d "$src/cohorts/$c/data/report" ]; then
      mkdir -p "$RUN/cohorts/$c/data/report" && cp -a "$src/cohorts/$c/data/report/." "$RUN/cohorts/$c/data/report/"
      slog "restored $c report: $(find "$RUN/cohorts/$c/data/report" -type f | wc -l) files"
    fi
    if [ -f "$RUN/cohorts/$c/validation-manifest.json" ] && [ ! -f "$RUN/cohorts/$c/data/report/index.json" ]; then
      stop "restore: cohort $c has a validation manifest but no data/report/index.json"
    fi
  done
  if [ -d "$src/audit" ]; then
    mkdir -p "$RUN/audit"; cp -an "$src/audit/." "$RUN/audit/"
    slog "restored run receipts: $(ls "$src/audit" | paste -sd' ')"
  fi
}
# synth_receipt NAME COHORT: the dry run's stand-in for a plan validator's receipt: the cohort's reports by sha256
synth_receipt() {
  python3 -I - "$RUN/audit/$1" "$2" "$RUN/cohorts/$2/data/species/$SLUG/results" << 'EOF'
import hashlib, json, os, sys
out, c, d = sys.argv[1:4]
reports = {f: hashlib.sha256(open(os.path.join(d, f), "rb").read()).hexdigest() for f in sorted(os.listdir(d))} if os.path.isdir(d) else {}
json.dump({"synthetic_stand_in": True, "cohort": c, "reports": reports}, open(out, "w"), indent=1, sort_keys=True)
EOF
}
# Q1 (main's ruling, 5 Nov): no pilot repeat; instead the prefix's sparse matrix must equal revision N's g5 matrix of
# the same genomes (R7's, built by the stock pipeline) byte for byte, checked right after matrix:prefix, before any
# other matrix, structure or learn hour. Any difference, or a file missing on either side, stops the box.
MATRIX_FILES="matrix.json patterns.bin support.bin kmer_offsets.bin kmers.bin"
prefix_vs_reference() {
  local REF=${MTB_REFERENCE_ROOT:?} out=$RUN/audit/prefix-matrix-vs-reference.json
  python3 -I - "$RUN/cohorts/prefix/data/species/$SLUG" "$REF/species/$SLUG" "$out" "$REVISION" $MATRIX_FILES << 'EOF'
import hashlib, json, os, sys
mine, ref, out, revision, names = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5:]
def sha(p):
    if not os.path.isfile(p): return None
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""): h.update(b)
    return h.hexdigest()
files = {}
for rel in ["genomes.txt"] + ["matrix/" + n for n in names]:
    a, b = sha(os.path.join(mine, rel)), sha(os.path.join(ref, rel))
    files[rel] = {"prefix": a, "reference": b, "equal": a is not None and a == b}
differ = [rel for rel, x in files.items() if not x["equal"]]
extra = {side: sorted(set(os.listdir(os.path.join(d, "matrix"))) - set(names)) if os.path.isdir(os.path.join(d, "matrix")) else None
         for side, d in (("prefix", mine), ("reference", ref))}
json.dump({"rule": "Q1: the prefix matrix equals revision %s's g5 matrix of the same genomes, byte for byte" % revision,
           "prefix": mine, "reference": ref, "files": files, "other_files_in_matrix": extra,
           "verdict": "PASS" if not differ else "STOP", "differ": differ}, open(out, "w"), indent=1, sort_keys=True)
print(" ".join(differ))
sys.exit(1 if differ else 0)
EOF
}
sec_matrix() {
  local c=$1 D=$RUN/cohorts/$1/data F diff
  gate "matrix:$c"
  if real; then python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || stop "check-lane $c"; else skip "check-lane $c"; fi
  if real || "$GRM" matrix --help 2> /dev/null | grep -q -- --disk-sparse; then
    # the dry run takes this path too when the source has --disk-sparse (the MTB merge)
    bash "$K/pg-run.sh" "$c-matrix" -- "$GRM" --data "$D" matrix --species $SLUG --min-support 2 --disk-sparse --block-width 128 || stop "matrix $c"
  else
    # r11/integration has no --disk-sparse or --block-width (the MTB branch adds them): the dry run counts, then
    # builds the stock dense matrix, under the same guard
    bash "$K/pg-run.sh" "$c-count" -- "$GRM" --data "$D" count --species $SLUG || stop "count $c"
    bash "$K/pg-run.sh" "$c-matrix" -- "$GRM" --data "$D" matrix --species $SLUG --min-support 2 || stop "matrix $c"
    skip "--disk-sparse --block-width 128 (not in r11/integration; stock count + matrix ran)"
  fi
  if [ "$c" = prefix ]; then
    if diff=$(prefix_vs_reference); then
      slog "Q1 PASS: the prefix genomes.txt and matrix ($MATRIX_FILES) equal revision $REVISION's g5 root's byte for byte (audit/prefix-matrix-vs-reference.json)"
    else
      stop "Q1: the prefix matrix is not byte-identical to revision $REVISION's g5 matrix of the same genomes (R7, $MTB_REFERENCE_ROOT); differs: ${diff:-receipt not written}; X1 stops before any other matrix or learn (audit/prefix-matrix-vs-reference.json)"
    fi
  fi
  if real; then
    python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
      --genomes $D/species/$SLUG/genomes.txt --fasta-dir $D/fasta --recount 16 --require-recount \
      --out $RUN/audit/$c-matrix-validation.json || stop "validate_matrix $c"
    python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || stop "verify-frozen"
  else skip "validate_matrix $c, verify-frozen"; fi
  if [ "$c" = matched2k ]; then
    F=$RUN/cohorts/full/data/species/$SLUG/matrix
    for x in prefix added matched2k; do
      D=$RUN/cohorts/$x/data
      if real; then
        python3 -I "$P/validate_matrix.py" --matrix $D/species/$SLUG/matrix --counted $D/species/$SLUG/counted.json \
          --genomes $D/species/$SLUG/genomes.txt --superset "$F" --out $RUN/audit/$x-in-full.json || stop "containment $x"
      else skip "containment $x"; fi
    done
  fi
}
sec_dupscreen() {
  gate dupscreen
  if real; then
    python3 -I "$P/screen_duplicates.py" --prefix-root $RUN/cohorts/prefix/data --added-root $RUN/cohorts/added/data \
      --fasta-manifest $RUN/fasta-manifest.tsv --genomes-jsonl $RUN/cohorts/full/data/bvbrc/genomes.jsonl \
      --out $RUN/audit/duplicate-screen.json || stop "duplicate screen"
  else
    python3 -I -c 'import json, sys; json.dump({"synthetic_stand_in": True, "cross_cohort_exact": 0}, open(sys.argv[1], "w"))' \
      "$RUN/audit/duplicate-screen.json"
    skip "screen_duplicates (stand-in audit/duplicate-screen.json written)"
  fi
}
sec_structure() {
  for c in $cohorts; do
    gate "structure:$c"
    bash "$K/pg-run.sh" "$c-structure" -- "$GRM" --data "$RUN/cohorts/$c/data" structure --species $SLUG || stop "structure $c"
  done
}
sec_learn() {
  local c=$1 C=$RUN/cohorts/$1 D=$RUN/cohorts/$1/data LAUNCH
  gate "learn:$c"
  if real; then python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || stop "check-lane $c"; fi
  mkdir -p "$RUN/backups/$c"
  LAUNCH=$(date -u +%Y-%m-%dT%H:%M:%SZ); echo "$LAUNCH" > "$RUN/logs/$c-learn.launch"
  RAYON_NUM_THREADS=$THREADS bash "$K/pg-run.sh" "$c-learn" -- "$GRM" --data "$D" learn --species $SLUG || stop "learn $c"
  cp "$W/logs/$c-learn.log" "$RUN/logs/$c-learn.log"
  if real; then
    python3 -I "$P/prepare_execution_roots.py" --check-annotations "$D" || stop "check-annotations $c"
    python3 -I "$P/prepare_execution_roots.py" --check-lane "$D" || stop "check-lane $c"
  else
    (cd "$D/bvbrc/features" && ls | wc -l) | grep -qx "$(wc -l < "$D/species/$SLUG/genomes.txt")" || stop "features changed in $c"
    skip "check-annotations (count only), check-lane $c"
  fi
  "$GRM" --data "$D" report > "$W/logs/$c-report.log" 2>&1 || stop "report $c"
  if real; then
    python3 -I "$P/prepare_execution_roots.py" --verify-frozen "$RUN" || stop "verify-frozen"
    python3 -I "$P/write_run_manifest.py" "$C" --cohort $c --source-commit "$SOURCE_COMMIT" --binary "$GRM" \
      --log $RUN/logs/$c-learn.log --profile "$PROFILE" || stop "write_run_manifest $c"
    python3 -I "$P/run_revision8_checks.py" $RUN/audit/$c-audit.json --root "$C" --run-manifest "$C/validation-manifest.json" \
      --expected-source-commit "$SOURCE_COMMIT" --binary "$GRM" --backup "$RUN/backups/$c" --log $RUN/logs/$c-learn.log \
      --launch "$LAUNCH" --foldgen "$FOLDGEN" --expected-foldgen-sha256 "$(sha256sum "$FOLDGEN" | cut -d' ' -f1)" \
      --checker-out $RUN/audit/$c-check-results.json --summary-out $RUN/audit/$c-summary.json \
      --receipt $RUN/audit/$c-revision8-receipt.json --profile "$PROFILE" || stop "run_revision8_checks $c"
  else
    python3 -I "$WS/scripts/check-results.py" "$D" --log "$RUN/logs/$c-learn.log" --out "$RUN/audit/$c-check-results.json" \
      > "$RUN/audit/$c-check-results.out" 2>&1 || stop "check-results $c"
    synth_receipt "$c-revision8-receipt.json" "$c"
    skip "write_run_manifest, run_revision8_checks $c (check-results.py ran; stand-in audit/$c-revision8-receipt.json written)"
  fi
  slog "learn $c: $(ls "$D/species/$SLUG/results" | wc -l) reports"
}
sec_compare() {
  gate compare
  if real; then
    local a=()
    for c in $cohorts; do a+=(--cohort "$c=$RUN/cohorts/$c" --checker-receipt "$c=$RUN/audit/$c-revision8-receipt.json"); done
    # the parent's reviewed acceptance of exact cross-cohort duplicates (plan section 8): the screen's sha256 from
    # stage.env; compare_cohorts accepts only that exact screen, and the decision record rides in audit/
    [ -n "${MTB_ACCEPT_DUPLICATES:-}" ] && { a+=(--accept-reviewed-duplicates "$MTB_ACCEPT_DUPLICATES"); slog "compare: accepting the reviewed duplicate screen $MTB_ACCEPT_DUPLICATES"; }
    python3 -I "$P/compare_cohorts.py" --profile "$PROFILE" "${a[@]}" --duplicates $RUN/audit/duplicate-screen.json \
      --expected-source-commit "$SOURCE_COMMIT" --expected-foldgen-sha256 "$(sha256sum "$FOLDGEN" | cut -d' ' -f1)" \
      --out $RUN/audit/compare-cohorts.json || stop "compare_cohorts"
  else
    # the dry run checks what compare_cohorts reads arrived: four receipts and the duplicate screen, each receipt's
    # reports equal to the cohort's reports here (prefix and added came from X2, full from X3, by restore)
    python3 -I - "$RUN" "$SLUG" << 'EOF' > "$RUN/audit/compare-inputs.json" || stop "compare inputs"
import hashlib, json, os, sys
run, slug = sys.argv[1:3]
out = {"duplicate_screen": os.path.isfile(os.path.join(run, "audit/duplicate-screen.json")), "cohorts": {}}
ok = out["duplicate_screen"]
for c in ("prefix", "added", "full", "matched2k"):
    p = os.path.join(run, "audit", c + "-revision8-receipt.json")
    d = os.path.join(run, "cohorts", c, "data/species", slug, "results")
    if not os.path.isfile(p):
        out["cohorts"][c] = "receipt missing"; ok = False; continue
    want = json.load(open(p))["reports"]
    have = {f: hashlib.sha256(open(os.path.join(d, f), "rb").read()).hexdigest() for f in sorted(os.listdir(d))}
    out["cohorts"][c] = {"reports": len(want), "match": want == have and bool(want)}
    ok = ok and want == have and bool(want)
out["ok"] = ok
print(json.dumps(out, indent=1, sort_keys=True))
sys.exit(0 if ok else 1)
EOF
    slog "compare inputs: $(python3 -I -c 'import json, sys; d = json.load(open(sys.argv[1])); print(json.dumps({"ok": d["ok"], "dupscreen": d["duplicate_screen"], "cohorts": d["cohorts"]}, sort_keys=True))' "$RUN/audit/compare-inputs.json")"
    skip "compare_cohorts (its inputs were checked instead)"
  fi
}
sec_reference() { # 7b: RN_REFERENCE_ROOT is revision N's g5 run root, pushed under /home/user/grm-ref
  local REF=${MTB_REFERENCE_ROOT:?} R f
  R=$REF/species/$SLUG
  gate reference
  sha256sum $R/genomes.txt $R/counted.json $R/structure.json $R/matrix/* $R/results/*.json > $RUN/audit/reference-sha256.txt || stop "reference files"
  cmp $RUN/cohorts/prefix/data/species/$SLUG/genomes.txt $R/genomes.txt || stop "prefix genomes.txt differs from the reference"
  for f in matrix.json patterns.bin support.bin kmer_offsets.bin kmers.bin; do
    if [ ! -e "$R/matrix/$f" ] && ! real; then slog "reference matrix has no $f (synthetic: the pipeline's dense matrix)"; continue; fi
    cmp $RUN/cohorts/prefix/data/species/$SLUG/matrix/$f $R/matrix/$f || stop "prefix matrix $f differs from the reference"
  done
  slog "reference: prefix genomes.txt and matrix equal revision $REVISION's g5 root"
  # K3 (main's ruling, 5 Nov): plan 7b's fold rebuild. audit.py --root takes a folder holding data/ and reads
  # data/species, data/report (index.json, datasets/) and data/curated/datasets.json, so the view links the pushed
  # reference's species/ and report/ and the input root's curated/ (the curation R7 learned from); every saved
  # test_fold must rebuild with this box's foldgen. The check (M1, REVIEW-MTB12) passes only with one audit record per
  # reference report (n, counted from the pushed reference), no test_fold issue and no error at all: audit.py records a
  # report it could not audit (unreadable, or a foldgen that answered nothing) as an error record and still exits 0.
  # The real run also needs n = 22, revision N's reports; the dry run's synthetic g5 root holds its own count, at its
  # own revision, and may name curated datasets its class gate skipped (allowed there only, as summary errors).
  local n rev arc synth=0
  n=$(ls "$R"/results/*.json 2> /dev/null | wc -l)
  if real; then rev=$REVISION; [ "$n" = 22 ] || stop "fold rebuild: the reference holds $n reports, not revision $REVISION's 22"; else
    synth=1
    rev=$(python3 -I -c 'import json, sys; print(json.load(open(sys.argv[1]))["evaluation"]["revision"])' "$(ls "$R"/results/*.json | head -n 1)") || stop "reference revision"
  fi
  mkdir -p "$RUN/rN-ref-root/data" && ln -s "$REF/species" "$RUN/rN-ref-root/data/species" &&
    ln -s "$REF/report" "$RUN/rN-ref-root/data/report" && ln -s "${MTB_INPUT_ROOT:?}/curated" "$RUN/rN-ref-root/data/curated" ||
    stop "fold-rebuild view"
  python3 -I "$WS/scripts/metric-audit/audit.py" "$RUN/audit/rN-fold-rebuild.json" $SLUG --revision "$rev" \
    --root "$RUN/rN-ref-root" --reshaped --foldgen "$FOLDGEN" > "$W/logs/fold-rebuild.log" 2>&1; arc=$?
  if [ $arc != 0 ] && ! real && python3 -I -c 'import json, sys; s = json.load(open(sys.argv[1]))["summary"]; sys.exit(0 if s["missing_curated_results"] and s["errors"] == s["missing_curated_results"] else 1)' "$RUN/audit/rN-fold-rebuild.json" 2> /dev/null; then
    # the synthetic g5 run made no report for a curated dataset its class gate skipped, which --reshaped can't see
    skip "audit.py's missing-curated-results exit (synthetic g5 root; the fold check below still runs)"; arc=0
  fi
  [ $arc = 0 ] || stop "fold rebuild (audit.py rc $arc; logs/fold-rebuild.log)"
  python3 -I -c 'import json,sys; r=json.load(open(sys.argv[1])); rs=r["records"]; s=r["summary"]; bad=[i for x in rs for i in x.get("issues",[]) if "test_fold" in i["where"]]; err=[x for x in rs if "error" in x]+[e for e in s["errors"] if not (sys.argv[3]=="1" and e in s["missing_curated_results"])]; print(len(rs),len(bad),len(err)); assert len(rs)==int(sys.argv[2]) and not bad and not err' \
    "$RUN/audit/rN-fold-rebuild.json" "$n" "$synth" > "$W/logs/fold-rebuild-check.log" 2>&1 ||
    stop "fold rebuild check (records, test_fold issues, errors; want $n, 0, 0): $(grep -m1 -E "^[0-9]+ [0-9]+ [0-9]+$" "$W/logs/fold-rebuild-check.log" || tail -n 1 "$W/logs/fold-rebuild-check.log")"
  slog "reference: fold rebuild of $(cut -d' ' -f1 "$W/logs/fold-rebuild-check.log") of $n reports (revision $rev), test_fold issues $(cut -d' ' -f2 "$W/logs/fold-rebuild-check.log"), errors $(cut -d' ' -f3 "$W/logs/fold-rebuild-check.log")"
  if real; then
    python3 -I "$P/fingerprint.py" "$RUN/rN-view/data" --make-view --matrix-from "$REF" --inputs-from "$MTB_INPUT_ROOT" \
      --reports-from "$REF" --check --expect-reports 22 --out "$RUN/audit/rN-fingerprints.json" || stop "fingerprint"
    python3 -I "$P/diff_reference_prefix.py" --prefix $RUN/cohorts/prefix/data --reference "$REF" \
      --out $RUN/audit/prefix-vs-reference.json || stop "diff_reference_prefix"
  else skip "fingerprint, diff_reference_prefix"; fi
}

for s in "$@"; do
  [ -e "$W/gates/sec-${s/:/-}.done" ] && { slog "$s already done"; continue; }
  slog "section $s"
  case "$s" in
    prepare) sec_prepare;; restore) sec_restore;; dupscreen) sec_dupscreen;; structure) sec_structure;;
    compare) sec_compare;; reference) sec_reference;;
    matrix:*) sec_matrix "${s#matrix:}";; learn:*) sec_learn "${s#learn:}";;
    *) stop "unknown section $s";;
  esac
  date -u +%FT%TZ > "$W/gates/sec-${s/:/-}.done"
done
