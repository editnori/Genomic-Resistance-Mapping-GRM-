#!/usr/bin/env bash
# One revision-N stage inside one Boat box (the laptop's r11-drive.sh calls each step through grm-boat.sh run).
# Ported from revision 10's driver.sh, build.sh and validate.sh (r10run; sha256s in REPORT.md): /root, ssh vast and
# KSD's HOLD are gone, every path is under /home/user, and the gates are PF-D7's.
#   host     build the pushed source here with the repo's target (x86-64-v3, .cargo/config.toml; main's correction to
#            D14, REVIEW-PORT M2), record the CPU model, avx512f and the binary
#            sha256s, check the source's EVALUATION_REVISION equals REVISION, then time the lib tests (the speed
#            check, PREFLIGHT section 5 item 2) and print HOST fast|slow.
#   verify   sha256sum -c --strict on the pushed data root against its manifest, whose sha256 must be
#            MANIFEST_SHA256; then every selected genome of the stage must have its FASTA and cached annotation, so
#            nothing is fetched from BV-BRC during the run; and the features manifest before the run is written.
#   bincheck the binary gate (D14 as main replaced it, 5 Nov), run only when this box's grm sha256 differs from the
#            run's reference box: g1a (CANARY_SPECIES) is counted and learned afresh on its own pushed inputs
#            (/home/user/grm-canary-data) and its reports must equal the reference box's (/home/user/grm-canary-ref)
#            byte for byte after the wall-clock fields only (report-bytes.py). Exit 0 on a match, 9 on a mismatch.
#   start    launch the stage detached (setsid) and return; status says how it goes
#   status   print one line: running|done|stopped ... and the last status text
#   stop     SIGTERM the stage's process groups (the TTL-near checkpoint)
# KIND=h1 (KIT10, G5): host is the same; verify and the stage are h1-box.sh's (verify, then run or combine by
# H1_ROLE; KIT12's addon role runs h1-box.sh run), whose other steps (prep, union, ckpt) the drive calls directly.
# status reports "paused" while the H-1 stage's guard holds its group stopped (gates/h1-<run|combine>.tripped ends in
# a pause; REVIEW-KIT10 L3).
# The stage, in order of SPECIES (comma list): one `scripts/run-pipeline.sh --species S` per species (PREFLIGHT
# section 2 item 2), each after the RAM gate (effective memory >= RAM_NEED_GIB, 13.5, PF-D7) and the disk gate
# (DISK_PAUSE_GIB + the species' kset estimate x 1.05 + 1 GiB), under guard.sh; then the species is validated
# (check-results.py and the metric audit on a view of it) and marked done. All species done: the features manifest
# after the run must equal the one before (no live fetch), and diff-rev.py writes the revision-9 diff of the stage's
# species against CANARY_REF as a report (classes per dataset, receipts/diff-rev9.json). It is never a gate.
# Settings come from W/stage.env (written by the laptop): STAGE SPECIES REVISION SOURCE_COMMIT MANIFEST_SHA256
# CANARY_REF CANARY_SPECIES CHANGED SKIP_SPECIES THREADS RAM_NEED_GIB RAM_PAUSE_GIB DISK_PAUSE_GIB SLOW_LIB_SECS
# RAM_WAIT_SECS AUDIT BUILD_TARGET_CPU (empty: the repo default, x86-64-v3; a value only in the dry run, to make a
# different binary).
# usage: r11-box.sh host|verify|bincheck|start|status|stop
set -u -o pipefail
W=${W:-/home/user/r11}
set -a; . "$W/stage.env"; set +a
# KIT12: an H-1 add-on box's stage is h1-box.sh run, so its guard's pause file is h1-run's
H1CMD=${H1_ROLE:-run}; [ "$H1CMD" = addon ] && H1CMD=run
KIND=${KIND:-r11}
K=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
S=${SRC:-/home/user/grm} DATA=${DATA:-/home/user/grm-data}
export W DATA
mkdir -p "$W/logs" "$W/gates" "$W/receipts"
export GRM=$S/target/release/grm FOLDGEN=$S/scripts/metric-audit/foldgen/target/release/foldgen
gib() { awk -v g="$1" 'BEGIN { printf "%d", g * 1073741824 }'; }
log() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$W/logs/driver.log"; }
status() { printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" > "$W/status.tmp" && mv "$W/status.tmp" "$W/status"; log "status: $*"; }
eff_mem() {
  local avail max anon c
  avail=$(($(awk '/MemAvailable/ {print $2}' /proc/meminfo) * 1024))
  max=$(cat /sys/fs/cgroup/memory.max 2> /dev/null || echo max)
  anon=$(awk '$1=="anon" {print $2}' /sys/fs/cgroup/memory.stat 2> /dev/null || echo 0)
  if [ "$max" != max ]; then c=$((max - anon)); [ "$c" -lt "$avail" ] && avail=$c; fi
  echo "$avail"
}
free_b() { df -B1 --output=avail "$W" | tail -1 | tr -d ' '; }
# PREFLIGHT section 3's kset estimates in GiB (measured where it says so); the disk gate keeps revision 10's form
kset_gib() {
  case $1 in
    enterococcus_faecium) echo 33;; staphylococcus_aureus) echo 31;; streptococcus_pneumoniae) echo 24;;
    pseudomonas_aeruginosa) echo 52;; shigella_sonnei) echo 53;; enterobacter_cloacae) echo 8;;
    neisseria_gonorrhoeae) echo 24;; acinetobacter_baumannii) echo 38;; mycobacterium_tuberculosis) echo 49;;
    escherichia_coli) echo 57;; klebsiella_pneumoniae) echo 62;; salmonella_enterica) echo 55;;
    *) echo 10;;
  esac
}
fail() { status "stopped stage=$STAGE $*"; echo "STOP $*" >&2; exit "${RC:-1}"; }

cmd_host() {
  local cpu avx t0 secs host
  cpu=$(grep -m1 'model name' /proc/cpuinfo | sed 's/.*: //')
  avx=$(grep -qw avx512f /proc/cpuinfo && echo yes || echo no)
  grep -q "pub const EVALUATION_REVISION: u32 = $REVISION;" "$S/src/analysis.rs" ||
    RC=2 fail "source EVALUATION_REVISION is not $REVISION: $(grep -o 'EVALUATION_REVISION: u32 = [0-9]*' "$S/src/analysis.rs")"
  # the repo's .cargo/config.toml (x86-64-v3) decides the target; no RUSTFLAGS unless the dry run names a CPU
  unset RUSTFLAGS CARGO_BUILD_RUSTFLAGS CARGO_ENCODED_RUSTFLAGS
  [ -n "${BUILD_TARGET_CPU:-}" ] && export RUSTFLAGS="-C target-cpu=$BUILD_TARGET_CPU"
  export CARGO_BUILD_JOBS=$(nproc)
  (cd "$S" && cargo build --release --locked -v) > "$W/logs/build.log" 2>&1 || RC=3 fail "cargo build failed (logs/build.log)"
  # N5: the target as the build used it. A fresh build prints grm's rustc line; if cargo found it fresh, touch the
  # crate root and build again with -v so the line is there, then take target-cpu from it
  if ! grep -q 'Running `.*rustc --crate-name grm ' "$W/logs/build.log"; then
    touch "$S/src/main.rs"; (cd "$S" && cargo build --release --locked -v) >> "$W/logs/build.log" 2>&1 || RC=3 fail "cargo build -v failed"
  fi
  BUILT_CPU=$(grep 'Running `.*rustc --crate-name grm ' "$W/logs/build.log" | tail -n 1 | grep -o 'target-cpu=[A-Za-z0-9_-]*' | tail -n 1)
  [ -n "$BUILT_CPU" ] || BUILT_CPU="target-cpu=(none in grm's rustc line: rustc's default)"
  CFG_FLAGS=$(python3 -I -c 'import sys, tomllib
c = tomllib.load(open(sys.argv[1], "rb")); f = c.get("build", {}).get("rustflags", [])
print(" ".join(f if isinstance(f, list) else f.split()))' "$S/.cargo/config.toml" 2> /dev/null)
  FEATURES=$(rustc --print cfg ${RUSTFLAGS:-$CFG_FLAGS} 2> /dev/null | sed -n 's/^target_feature="\(avx2\|avx512f\|bmi2\|popcnt\|fma\)"$/\1/p' | sort | paste -sd, -)
  (cd "$S/scripts/metric-audit/foldgen" && cargo build --release --locked) > "$W/logs/build-foldgen.log" 2>&1 ||
    RC=3 fail "foldgen build failed"
  (cd "$S" && cargo test --release --locked --lib --no-run) > "$W/logs/build-tests.log" 2>&1 || RC=3 fail "test build failed"
  t0=$(date +%s.%N)
  (cd "$S" && cargo test --release --locked --lib -q) > "$W/logs/speed-lib-tests.log" 2>&1 || RC=3 fail "lib tests failed (logs/speed-lib-tests.log)"
  secs=$(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN { printf "%.1f", b - a }')
  host=fast; awk -v s="$secs" -v m="$SLOW_LIB_SECS" 'BEGIN { exit !(s > m) }' && host=slow
  {
    echo "stage=$STAGE revision=$REVISION source_commit=$SOURCE_COMMIT"
    echo "cpu=$cpu"; echo "avx512f=$avx nproc=$(nproc) memtotal_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
    echo "rustflags_env=${RUSTFLAGS:-(none)} config_build_rustflags=[${CFG_FLAGS:-none}]"
    echo "built_with=$BUILT_CPU (from grm's rustc line in logs/build.log) target_features=$FEATURES"
    rustc --version; cargo --version
    echo "cargo_config_sha256=$(sha256sum "$S/.cargo/config.toml" | cut -d' ' -f1)"
    sha256sum "$GRM" "$FOLDGEN"
    echo "speed_lib_tests_secs=$secs slow_above=$SLOW_LIB_SECS host=$host"
  } > "$W/receipts/host.txt"
  chmod a-w "$GRM"
  cat "$W/receipts/host.txt"
  echo "HOST $host"
}

cmd_verify() {
  local m=$DATA/.grm-input-manifest.sha256 sp miss
  [ "$(sha256sum "$m" | cut -d' ' -f1)" = "$MANIFEST_SHA256" ] || RC=2 fail "manifest in the box is not $MANIFEST_SHA256"
  (cd "$DATA" && sha256sum -c --strict --quiet "$m") > "$W/logs/verify.log" 2>&1 || RC=2 fail "sha256sum -c failed in the box"
  echo "$(date -u +%FT%TZ) manifest $MANIFEST_SHA256 lines=$(wc -l < "$m") rc=0" > "$W/receipts/manifest-check.txt"
  for x in grm-resume grm-ref; do
    [ -d /home/user/$x ] || continue
    (cd /home/user/$x && sha256sum -c --strict --quiet .grm-input-manifest.sha256) || RC=2 fail "$x set check failed"
    echo "$x $(sha256sum /home/user/$x/.grm-input-manifest.sha256)" >> "$W/receipts/manifest-check.txt"
  done
  [ "$KIND" = r11 ] && [ -d /home/user/grm-resume/species ] && cp -an /home/user/grm-resume/species "$DATA/"
  sp=${SPECIES//@canary/}; sp=${sp//,,/,}
  [ "$KIND" = mtb ] && sp=mycobacterium_tuberculosis
  miss=$(python3 -I - "$DATA" "$sp" << 'EOF'
import os, sys
d, sp = sys.argv[1], [s for s in sys.argv[2].split(",") if s]
bad = []
for s in sp:
    g = os.path.join(d, "species", s, "genomes.txt")
    if not os.path.isfile(g):
        bad.append("%s: no genomes.txt" % s); continue
    for i in open(g).read().split():
        for p in ("fasta/%s.fna.gz" % i, "bvbrc/features/%s.json" % i):
            if not os.path.isfile(os.path.join(d, p)):
                bad.append(p)
print("\n".join(bad))
EOF
)
  [ -z "$miss" ] || RC=2 fail "$(wc -l <<< "$miss") inputs missing, so the run would fetch them: $(head -n 3 <<< "$miss" | tr '\n' ' ')"
  ln -sfn "$DATA" "$S/data"
  (cd "$DATA/bvbrc/features" && find . -type f -print0 | sort -z | xargs -0 -r sha256sum) > "$W/receipts/features-before.sha256"
  mkdir -p "$W/prelaunch-empty/data/species"
  echo "verified: $(wc -l < "$m") files, every selected genome has FASTA and annotation"
}

# validate_species SLUG LOG LAUNCH: check-results and the metric audit on a view holding only SLUG
validate_species() {
  local s=$1 lg=$2 launch=$3 v=$W/views/$1 o=$W/receipts/sp-$1 rc=0
  rm -rf "$v" "$o"; mkdir -p "$v/data/species" "$v/data/curated" "$o" "$W/prelaunch-empty/data/species/$s/results"
  for f in "$DATA"/curated/*; do [ "$(basename "$f")" = datasets.json ] || ln -s "$f" "$v/data/curated/"; done
  python3 -I -c 'import json, sys
rows = json.load(open(sys.argv[1])); json.dump([r for r in rows if r["species_slug"] == sys.argv[2]], open(sys.argv[3], "w"), indent=2)' \
    "$DATA/curated/datasets.json" "$s" "$v/data/curated/datasets.json"
  ln -s "$DATA/report" "$v/data/report"; ln -s "$DATA/species/$s" "$v/data/species/$s"
  nice -n 10 python3 -I "$S/scripts/check-results.py" "$v/data" --log "$lg" --out "$o/check-results.json" > "$o/check-results.out" 2>&1 || rc=1
  if [ "${AUDIT:-1}" = 1 ]; then
    nice -n 10 python3 -I "$S/scripts/metric-audit/audit.py" "$o/audit.json" "$s" --revision "$REVISION" --root "$v" \
      --backup "$W/prelaunch-empty" --log "$lg" --launch "$launch" --foldgen "$FOLDGEN" > "$o/audit.out" 2>&1 || rc=1
  fi
  echo "$rc" > "$W/gates/sp-$s.validated"
  return $rc
}

rev9_report() { # the stage's species against CANARY_REF (revision 9): classes per dataset, a report, never a gate
  local sp=$1 o=$W/receipts
  [ -d "${CANARY_REF:-/nonexistent}" ] || { log "no revision-9 tree at ${CANARY_REF:-unset}; no diff report"; return 0; }
  python3 -I "$K/diff-rev.py" "$DATA" "$CANARY_REF" "$o/diff-rev9.json" --new-revision "$REVISION" --old-revision 9 \
    --audit "$S/scripts/metric-audit/audit.py" --changed "${CHANGED:-}" ${sp//,/ } > "$o/diff-rev9.out" 2>&1
  echo $? > "$o/diff-rev9.rc"
  log "revision-9 diff (report only): $(head -n 1 "$o/diff-rev9.out" | cut -c1-300)"
}

cmd_bincheck() {
  local CD=/home/user/grm-canary-data REF=/home/user/grm-canary-ref o=$W/receipts/bincheck rc miss
  mkdir -p "$o"
  for x in "$CD" "$REF"; do
    (cd "$x" && sha256sum -c --strict --quiet .grm-input-manifest.sha256) || RC=2 fail "bincheck: $x does not match its manifest"
  done
  miss=$(for s in ${CANARY_SPECIES//,/ }; do for i in $(cat "$CD/species/$s/genomes.txt"); do
    [ -f "$CD/fasta/$i.fna.gz" ] && [ -f "$CD/bvbrc/features/$i.json" ] || echo "$s/$i"; done; done | head -n 3)
  [ -z "$miss" ] || RC=2 fail "bincheck: g1a inputs missing, so the run would fetch them: $miss"
  status "bincheck stage=$STAGE: this binary differs from the reference; g1a ($CANARY_SPECIES) reruns"
  # run-pipeline.sh's settings (12 workers, --drop-ksets, nice 10) on the canary root, under the guard, from the
  # checkout as run-pipeline.sh runs: pipeline reads reference/kover-summary.json from the working folder only and
  # refuses before counting without it (M1, REVIEW-INTEG7)
  cd "$S" || RC=2 fail "bincheck: cannot enter the checkout $S"
  RAYON_NUM_THREADS=$THREADS DATA=$CD SKIP_SPECIES= bash "$K/pg-run.sh" bincheck-g1a -- \
    nice -n 10 "$GRM" --data "$CD" pipeline --workers 12 --drop-ksets --species "$CANARY_SPECIES"
  rc=$?; [ $rc = 0 ] || RC=1 fail "bincheck: g1a pipeline exited $rc (logs/bincheck-g1a.log)"
  python3 -I "$K/report-bytes.py" "$REF" "$CD" "$CANARY_SPECIES" "$o/report-bytes.json" > "$o/report-bytes.out" 2>&1; rc=$?
  { echo "grm_sha256=$(sha256sum "$GRM" | cut -d' ' -f1)"
    echo "ref_manifest=$(sha256sum "$REF/.grm-input-manifest.sha256" | cut -d' ' -f1)"
    echo "canary_data_manifest=$(sha256sum "$CD/.grm-input-manifest.sha256" | cut -d' ' -f1)"
    cat "$o/report-bytes.out"; } > "$o/verdict.txt"
  cat "$o/verdict.txt"
  if [ $rc = 0 ]; then log "bincheck matched: $(tail -n 1 "$o/verdict.txt")"; echo "BINCHECK match"; return 0; fi
  if [ $rc = 9 ]; then status "stopped stage=$STAGE bincheck MISMATCH: $(tail -n 1 "$o/verdict.txt")"; echo "BINCHECK mismatch"; exit 9; fi
  RC=2 fail "bincheck: report-bytes.py exited $rc"
}

run_stage() {
  local sp list=() done_before="" s tag lg pg start guard rc need waited t0 launch
  status "running stage=$STAGE"
  if [ "$KIND" = h1 ]; then
    # KIT12: an add-on box runs h1-box.sh run (with --reference-from)
    case "${H1_ROLE:-run}" in run|addon|combine) ;; *) RC=2 fail "H1_ROLE is ${H1_ROLE}, not run, addon or combine";; esac
    bash "$K/h1-box.sh" "$H1CMD" >> "$W/logs/h1-stage.log" 2>&1 || RC=1 fail "H-1 ${H1_ROLE:-run} failed: $(grep '^STOP' "$W/logs/h1-stage.log" | tail -n 1)"
    date -u +%FT%TZ > "$W/done"; status "done stage=$STAGE: H-1 ${H1_ROLE:-run} of $SPEC"; return
  fi
  if [ "$KIND" = mtb ]; then
    GRM=$GRM FOLDGEN=$FOLDGEN bash "$K/mtb-box.sh" ${SPECIES//,/ } || RC=1 fail "MTB section failed: $(tail -n 1 "$W/logs/sections.log")"
    date -u +%FT%TZ > "$W/done"; status "done stage=$STAGE: sections ${SPECIES}"; return
  fi
  IFS=, read -r -a list <<< "$SPECIES"
  for s in "${list[@]}"; do
    [ "$s" = @canary ] && { log "@canary is no longer a gate (the binary is, 5 Nov); ignored"; continue; }
    done_before="$done_before,$s"
    [ -e "$W/gates/sp-$s.done" ] && { log "$s already done"; continue; }
    [[ " ${SKIP_SPECIES//,/ } " == *" $s "* ]] && { log "$s carried from an earlier box; skipped"; date -u +%FT%TZ > "$W/gates/sp-$s.done"; continue; }
    need=$(( $(gib "$DISK_PAUSE_GIB") + $(gib "$(kset_gib "$s")") * 105 / 100 + $(gib 1) ))
    waited=0 t0=$(date +%s)
    while :; do
      free=$(free_b) mem=$(eff_mem)
      [ "$free" -ge "$need" ] && [ "$mem" -ge "$(gib "$RAM_NEED_GIB")" ] && break
      [ $(( $(date +%s) - t0 )) -ge "$RAM_WAIT_SECS" ] && RC=5 fail "gates not met for $s after ${RAM_WAIT_SECS}s: free=$free need=$need eff_mem=$mem need_mem=$(gib "$RAM_NEED_GIB")"
      status "waiting stage=$STAGE species=$s free=$free need_free=$need eff_mem=$mem need_mem=$(gib "$RAM_NEED_GIB")"
      waited=1; sleep 30
    done
    log "$s gates met free=$free eff_mem=$mem (need $need, $(gib "$RAM_NEED_GIB"))$([ $waited = 1 ] && echo ' after a wait')"
    tag=$s-a1; while [ -e "$W/gates/$tag.launched" ]; do tag=$s-a$(( ${tag##*-a} + 1 )); done
    lg=$DATA/logs/pipeline-$tag.log; mkdir -p "$DATA/logs"
    launch=$(date -u +%FT%TZ)
    cd "$S"
    RAYON_NUM_THREADS=$THREADS setsid bash -c 'bash scripts/run-pipeline.sh --species "$1" > "$2" 2>&1; echo $? > "$3"' _ "$s" "$lg" "$W/gates/$tag.exit" < /dev/null > /dev/null 2>&1 &
    pg=$!
    for _ in $(seq 50); do [ "$(sed 's/^.*) //' "/proc/$pg/stat" 2> /dev/null | cut -d' ' -f3)" = "$pg" ] && break; sleep 0.1; done
    start=$(sed 's/^.*) //' "/proc/$pg/stat" | cut -d' ' -f20)
    SKIP_SPECIES=${SKIP_SPECIES:-} setsid bash "$K/guard.sh" "$pg" "$start" "$tag" < /dev/null > /dev/null 2>&1 &
    guard=$!
    printf '%s stage=%s species=%s pgid=%s starttime=%s guard=%s threads=%s grm_sha256=%s free=%s eff_mem=%s\n' "$launch" "$STAGE" "$s" "$pg" \
      "$start" "$guard" "$THREADS" "$(sha256sum "$GRM" | cut -d' ' -f1)" "$free" "$mem" > "$W/gates/$tag.launched"
    echo "$pg" >> "$W/gates/pgids"
    status "running stage=$STAGE species=$s pgid=$pg"
    wait "$pg"; rc=$(cat "$W/gates/$tag.exit" 2> /dev/null || echo 99)
    printf '%s exit=%s\n' "$(date -u +%FT%TZ)" "$rc" > "$W/gates/$tag.ended"
    [ -e "$W/gates/$tag.tripwire" ] && RC=6 fail "tripwire in $tag: $(cat "$W/gates/$tag.tripwire")"
    [ -e "$W/stop" ] && RC=7 fail "stopped by request (TTL-near checkpoint) during $s"
    [ "$rc" = 0 ] || RC=1 fail "pipeline for $s exited $rc ($lg)"
    status "validating stage=$STAGE species=$s"
    validate_species "$s" "$lg" "$launch" || RC=2 fail "validation failed for $s (receipts/sp-$s)"
    date -u +%FT%TZ > "$W/gates/sp-$s.done"
    log "$s done"
    # dry run only (unset in real runs): hold after each species so a forced TTL-near stop lands mid-stage
    for _ in $(seq "${SYNTH_PAUSE_AFTER_SPECIES:-0}"); do [ -e "$W/stop" ] && RC=7 fail "stopped by request after $s"; sleep 1; done
  done
  (cd "$DATA/bvbrc/features" && find . -type f -print0 | sort -z | xargs -0 -r sha256sum) > "$W/receipts/features-after.sha256"
  cmp -s "$W/receipts/features-before.sha256" "$W/receipts/features-after.sha256" || RC=8 fail "bvbrc/features changed during the run: a live fetch happened"
  rev9_report "${done_before#,}"
  date -u +%FT%TZ > "$W/done"
  status "done stage=$STAGE: every species validated, features unchanged"
}

case "${1:-}" in
  host) cmd_host;;
  verify) if [ "$KIND" = h1 ]; then exec bash "$K/h1-box.sh" verify; fi; cmd_verify;;
  bincheck) cmd_bincheck;;
  start)
    [ -e "$W/done" ] && { echo "already done"; exit 0; }
    rm -f "$W/stop"
    setsid nohup bash "$0" _run < /dev/null > "$W/logs/stage.out" 2>&1 &
    echo "$!" > "$W/stage.pid"; echo "started pid $(cat "$W/stage.pid")";;
  _run) run_stage;;
  status)
    st=$(cat "$W/status" 2> /dev/null || echo "none")
    # REVIEW-KIT10 L3: an H-1 stage whose guard.sh stopped its group (the last line of gates/h1-<role>.tripped is a
    # pause) reports paused with the guard's reading, so the drive stops with that reason
    gp=""; [ "$KIND" = h1 ] && gp=$(tail -n 1 "$W/gates/h1-$H1CMD.tripped" 2> /dev/null | grep '^paused ')
    if [ -e "$W/done" ]; then echo "STATE done | $st"
    elif [ -s "$W/stage.pid" ] && kill -0 "$(cat "$W/stage.pid")" 2> /dev/null; then
      if [ -n "$gp" ]; then
        why=$(awk -v d="${DISK_PAUSE_GIB:-30}" -v m="${RAM_PAUSE_GIB:-4}" '{ for (i = 1; i <= NF; i++) { split($i, a, "="); v[a[1]] = a[2] } }
          END { g = 1073741824; r = ""
                if (v["disk_free"] != "" && v["disk_free"] + 0 < d * g) r = "disk: " v["disk_free"] " B free, below DISK_PAUSE_GIB " d " GiB"
                if (v["eff_mem"] != "" && v["eff_mem"] + 0 < m * g) r = (r != "" ? r "; " : "") "memory: " v["eff_mem"] " B, below RAM_PAUSE_GIB " m " GiB"
                print (r != "" ? r : "a pause") }' <<< "$gp")
        echo "STATE paused | guard $why ($gp) | $st"
      else echo "STATE running | $st"; fi
    else echo "STATE stopped | $st"; fi;;
  stop)
    touch "$W/stop"
    for p in $(cat "$W/gates/pgids" 2> /dev/null); do kill -CONT -- "-$p" 2> /dev/null; kill -TERM -- "-$p" 2> /dev/null; done
    for _ in $(seq 60); do [ -s "$W/stage.pid" ] && kill -0 "$(cat "$W/stage.pid")" 2> /dev/null || break; sleep 2; done
    echo "stopped: $(cat "$W/status" 2> /dev/null)";;
  *) sed -n '2,/^set -u/p' "$0" | sed '$d'; exit 2;;
esac
