#!/usr/bin/env bash
# drive.sh: run one stage box of the revision-N real-data run or the MTB cleaned 4k run on Boat, through grm-boat.sh,
# from WSL. One call = one box, start to delete. Settings come from a run file (RUN_CONF, sourced; see run-r11.conf
# and run-mtb.conf) and the stage table it names. In order:
#   1. budget guard (Budget row, 5 Nov, with main's L9 ruling): the run's ledger (RUN_DIR/budget.tsv, xlarge-hours
#      from boat usage, or wall time x the type's multiplier when usage fails) at or past BUDGET_STOP (25) is a
#      report-and-continue checkpoint, logged and printed loudly; ledger + this stage's slow-host estimate past
#      BUDGET_CAP (35) is refused, and two slow hosts in a row stop the run. Exit 78.
#   2. order: the binary gate (D14 as main replaced it, 5 Nov), one reference per source commit C (M3 ruling).
#      Revision N (KIND=r11): REF_BOX (R1) is the first box; its grm sha256 becomes RUN_DIR/BINARY-REF-C.sha at its
#      host step, and its g1a (CANARY_SPECIES) reports, once back and hashed, make RUN_DIR/BINARY-REF-C. No other r11
#      box starts before that (79). MTB (KIND=mtb): the first MTB box's grm becomes BINARY-REF-C.sha (binary only;
#      MTB runs no g1a). BINARY-MISMATCH-C (or a run-wide BINARY-MISMATCH) stops every box of that source (82). A
#      stage with NEEDS=<file> ($RUN_DIR expands) waits for that file (MTB needs revision N's g5 run root). Exit 79.
#   2b. MTB source checks (K2, main's ruling of 5 Nov): the box's copy of SRC has no .git, so here, before any box,
#      SRC's HEAD must be the run file's SRC_COMMIT (the reviewed commit) in a clean worktree (L1, REVIEW-MTB12;
#      receipt src-commit.txt), SRC's own scripts/mtb-4k/revision_n.py verify-checkout SRC must pass and SRC's foldgen tree must be the pinned
#      dc8317f...; the receipts are the stage folder's verify-checkout.json and foldgen-tree.txt. Either failing is
#      exit 2 with no box.
#   3. inputs: the master manifest must hash to MASTER_SHA256; the stage's manifest is derived from it (shared files,
#      then each species' genomes.txt, FASTA and annotations, plus the canary's revision-9 report tree), written
#      under RUN_DIR/stages/<box>/ and its sha256 recorded in NOTES.md before anything is pushed (PF-D2).
#   4. plan floor and laptop check: grm-boat.sh preflight (76: wait 5 min, 3 tries; 77: stop), then new, which
#      checks both again.
#   5. box: push the source worktree and the kit, write stage.env, build and speed-check (r11-box.sh host). A slow
#      host on a stage whose slow estimate breaks its ceiling closes the box and exits 80; a slow host on a SLOW_OK
#      stage sets its remaining lifetime once so it ends at slow estimate + 30 min after creation (at most MAX_TTL,
#      6 h, PF-D6). The build uses the repo's default target (x86-64-v3, .cargo/config.toml; M2). Then the binary
#      gate: a grm sha256 equal to the reference goes straight on. An r11 box with any other pushes g1a's inputs and
#      the reference g1a reports (two push-data sets) and runs r11-box.sh bincheck, which must match byte for byte
#      after the wall-clock fields; an MTB box with any other stops. A mismatch writes BINARY-MISMATCH-C, pulls the
#      receipts, closes the box and exits 82. FORCE_BINARY_MISMATCH=1 treats an equal hash as different (dry run).
#   6. data: push-data with the stage manifest and its pin (and a resume set, if an earlier box finished species),
#      then r11-box.sh verify, then start.
#   7. watch: status every POLL s; a checkpoint copy-back every CKPT_SECS (1800, PF-D6), each pulled and checked
#      with sha256sum -c --strict on the laptop; at TTL - TTL_MARGIN a stop, a full checkpoint and close: exit 74
#      when the box finished a species (the next call resumes from them), exit 83 when it finished nothing (or is an
#      MTB box, whose sections don't resume): rerun only with TTL_OVERRIDE_S raised (6 h, PF-D6); a rerun with no
#      larger TTL is refused with 83 before any box. FORCE_TTL_NEAR_AFTER=s forces that path for a dry run.
#   Any exit after the box exists (an error, INT, TERM, HUP) runs the exit trap: it stops the stage, pulls a full
#   checkpoint when the stage had started, books usage and deletes the box. A failed `new` looks up a live box of the
#   lane's name made in this call and deletes it (M4).
#   8. end: a full checkpoint, pulled and checked, merged into RUN_DIR/stages/<box>/out with its own SHA256SUMS and a
#      second copy under COPY2 (D11) when set; boat usage into the ledger; grm-boat.sh close.
# Exit: 0 stage done; 82 binary gate mismatch (the run stops); 74 TTL checkpoint; 83 TTL reached with nothing
# finished (rerun with TTL_OVERRIDE_S); 78 budget; 79 order; 80 slow host; 81 the stage stopped
# (pipeline, gate or validation failure: see the pulled run/status); 2 usage; others from grm-boat.sh.
# usage: RUN_CONF=run-r11.conf drive.sh <box>
set -u -o pipefail
: "${RUN_CONF:?set RUN_CONF}"; . "$RUN_CONF"
BOX=${1:?box name, such as R1}
# no @@...@@ placeholder in any variable, whatever the run file named it (REVIEW-PORT L1): COPY2, CHANGED, EXTRA_ENV...
PH=$(set | grep -aE '^[A-Za-z_][A-Za-z0-9_]*=.*@@[A-Za-z0-9_-]*@@' | cut -d= -f1 | sort -u | paste -sd' ')
[ -z "$PH" ] || { echo "drive: placeholder left in: $PH" >&2; exit 2; }
# a run file with a placeholder left in it is refused before anything runs
# the evaluation revision is the source tree's; a configured REVISION must equal it (main's ruling, 5 Nov)
SRC_REV=$(sed -n 's/^pub const EVALUATION_REVISION: u32 = \([0-9][0-9]*\);$/\1/p' "${SRC:-/nonexistent}/src/analysis.rs" 2> /dev/null)
[[ "$SRC_REV" =~ ^[0-9]+$ ]] || { echo "drive: no EVALUATION_REVISION in ${SRC:-SRC unset}/src/analysis.rs" >&2; exit 2; }
[ -z "${REVISION:-}" ] && REVISION=$SRC_REV
[ "$REVISION" = "$SRC_REV" ] || { echo "drive: REVISION=$REVISION but the source says $SRC_REV" >&2; exit 2; }
[[ "${MASTER_SHA256:-}" =~ ^[0-9a-f]{64}$ ]] || { echo "drive: MASTER_SHA256 is not a sha256: ${MASTER_SHA256:-unset}" >&2; exit 2; }
for v in RUNNER NOTES RUN_DIR STAGES KIND SRC KIT KIT_REPO DATA_ROOT MASTER_MANIFEST LANE_PREFIX RAM_NEED_GIB; do
  case "${!v:-}" in ""|*@@*) echo "drive: $v is unset or a placeholder" >&2; exit 2;; esac
done
# I-g: SRC must be the top of a git worktree with a commit, and clean (a plain archive extract has no commit key)
[ "$(git -C "$SRC" rev-parse --show-toplevel 2> /dev/null)" = "$(readlink -f "$SRC")" ] &&
  [[ "$(git -C "$SRC" rev-parse HEAD 2> /dev/null)" =~ ^[0-9a-f]{40}$ ]] || { echo "drive: $SRC is not a git worktree with a commit" >&2; exit 2; }
st=$(git -C "$SRC" status --porcelain 2> /dev/null) && [ -z "$st" ] || { echo "drive: $SRC is not a clean worktree" >&2; exit 2; }
R=${RUNNER:?} NOTES=${NOTES:?} RUN_DIR=${RUN_DIR:?}
BOATX=${GRM_BOAT:-/mnt/c/Users/layth888/.ascii/bin/boat.exe}
BUDGET_CAP=${BUDGET_CAP:-35} BUDGET_STOP=${BUDGET_STOP:-25} CKPT_SECS=${CKPT_SECS:-1800} POLL=${POLL:-60}
TTL_MARGIN=${TTL_MARGIN:-900} MAX_TTL=${MAX_TTL:-21600} XL_MULT=5.56
export GRM_NOTES=$NOTES
SD=$RUN_DIR/stages/$BOX; mkdir -p "$SD" "$RUN_DIR"
# N4: one drive.sh per stage at a time; a second call refuses before it can make or delete anything
exec 8> "$SD/.drive.lock"
flock -n 8 || { echo "drive: another drive.sh is running stage $BOX (lock $SD/.drive.lock); refused" >&2; exit 84; }
BSTATE=${GRM_BOAT_STATE:-$HOME/.local/state/grm-boat}   # where grm-boat.sh records each box it makes
LEDGER=$RUN_DIR/budget.tsv
now() { date -u +%FT%TZ; }
say() { printf -- '- %s | %s | drive %s | %s\n' "$(now)" "${ID:--}" "$BOX" "$*" >> "$NOTES"; echo "drive $BOX: $*" >&2; }
die() { local rc=$1; shift; say "STOP rc=$rc: $*"; exit "$rc"; }
bt() { local c=$1; shift; timeout -k 10 180 "$BOATX" "$c" --no-update "$@" < /dev/null; }

row=$(awk -F'\t' -v b="$BOX" '$1 == b' "$STAGES")
[ -n "$row" ] || die 2 "no stage $BOX in $STAGES"
IFS=$'\t' read -r _ GATES TYPE EST_FAST EST_SLOW SLOW_OK NEEDS SPECIES <<< "$row"

# ---- 1. budget
used=$(awk -F'\t' 'NR > 1 {s += $6} END {printf "%.3f", s + 0}' "$LEDGER" 2> /dev/null || echo 0)
[ -s "$LEDGER" ] || printf 'utc\tbox\tid\tcpu\thost\txlarge_hours\tvm_seconds\tresult\n' > "$LEDGER"
if awk -v u="$used" -v s="$BUDGET_STOP" 'BEGIN { exit !(u >= s) }'; then
  # main's L9 ruling: report and continue; the 35 cap below still refuses
  msg="BUDGET CHECKPOINT: $used xlarge-hours used, at or past $BUDGET_STOP; report to main (the run continues; the cap is $BUDGET_CAP)"
  say "$msg"; printf '\n******** %s ********\n\n' "$msg" >&2
  [ -e "$RUN_DIR/BUDGET-CHECKPOINT-$BUDGET_STOP" ] || printf '%s %s\n' "$(now)" "$msg" > "$RUN_DIR/BUDGET-CHECKPOINT-$BUDGET_STOP"
fi
awk -v u="$used" -v e="$EST_SLOW" -v c="$BUDGET_CAP" 'BEGIN { exit !(u + e > c) }' && die 78 "budget: $used used + $EST_SLOW (slow-host estimate) passes the $BUDGET_CAP cap; refused"
[ "$(awk -F'\t' 'NR > 1 {print $5}' "$LEDGER" | tail -n 2 | grep -c '^slow$')" = 2 ] && die 78 "budget: the last two boxes had slow hosts; stop and report"
say "budget ok: $used xlarge-hours used, stage estimate $EST_FAST fast / $EST_SLOW slow, checkpoint $BUDGET_STOP, cap $BUDGET_CAP"

# ---- 2. order
REF_BOX=${REF_BOX-R1} CANARY_SPECIES=${CANARY_SPECIES:-}
C=$(git -C "$SRC" rev-parse HEAD)           # one reference binary per source commit (M3)
REFK=$RUN_DIR/BINARY-REF-$C MISK=$RUN_DIR/BINARY-MISMATCH-$C
[ -e "$RUN_DIR/BINARY-MISMATCH" ] && die 82 "the binary gate failed ($(head -n 1 "$RUN_DIR/BINARY-MISMATCH")); the run is stopped"
[ -e "$MISK" ] && die 82 "the binary gate failed for source $C ($(head -n 1 "$MISK")); every box of it is stopped"
if [ "$KIND" = r11 ] && [ -n "$REF_BOX" ]; then
  [ -n "$CANARY_SPECIES" ] || die 2 "REF_BOX is set but CANARY_SPECIES (g1a) is empty"
  [ "$BOX" != "$REF_BOX" ] && [ ! -e "$REFK" ] && die 79 "the reference box $REF_BOX has not fixed source $C's binary and its g1a reports yet"
fi
if [ "$KIND" = mtb ]; then
  MTB_REF_BOX=${MTB_REF_BOX:-X1}
  [ "$BOX" != "$MTB_REF_BOX" ] && [ ! -e "$REFK.sha" ] && die 79 "the MTB reference box $MTB_REF_BOX has not fixed source $C's binary yet"
fi
NEEDS=${NEEDS//\$RUN_DIR/$RUN_DIR}
[ "$NEEDS" != - ] && [ ! -e "$NEEDS" ] && die 79 "needs $NEEDS"
[ -e "$SD/DONE" ] && die 0 "already done ($(cat "$SD/DONE"))"

# ---- 2b. MTB source checks that need git (K2): on the laptop, against SRC, before any box
FOLDGEN_TREE=dc8317faf778fd96f496cc0226515ca7b19d688d   # plan section 3's pin (I4)
if [ "$KIND" = mtb ]; then
  rm -f "$SD/src-commit.txt" "$SD/verify-checkout.json" "$SD/verify-checkout.err" "$SD/foldgen-tree.txt"
  # L1 (REVIEW-MTB12): the run uses exactly the reviewed commit, from a clean worktree
  st=$(git -C "$SRC" status --porcelain --untracked-files=all 2> /dev/null); src_rc=$?
  printf '%s SRC %s HEAD %s SRC_COMMIT %s clean %s\n' "$(now)" "$SRC" "$C" "${SRC_COMMIT:-unset}" \
    "$([ $src_rc = 0 ] && [ -z "$st" ] && echo yes || echo no)" > "$SD/src-commit.txt"
  [[ "${SRC_COMMIT:-}" =~ ^[0-9a-f]{40}$ ]] || die 2 "K2: SRC_COMMIT is not set to the reviewed commit (a 40-hex sha): ${SRC_COMMIT:-unset}"
  [ "$C" = "$SRC_COMMIT" ] || die 2 "K2: $SRC is at $C, not the reviewed SRC_COMMIT $SRC_COMMIT"
  [ $src_rc = 0 ] && [ -z "$st" ] || die 2 "K2: $SRC is not a clean worktree"
  [ -f "$SRC/scripts/mtb-4k/revision_n.py" ] || die 2 "K2: $SRC has no scripts/mtb-4k/revision_n.py; it is not the MTB merge"
  (cd "$SD" && python3 -I "$SRC/scripts/mtb-4k/revision_n.py" verify-checkout "$SRC" > "$SD/verify-checkout.json" 2> "$SD/verify-checkout.err"); vrc=$?
  [ $vrc = 0 ] || die 2 "K2: revision_n.py verify-checkout refused $SRC (rc $vrc): $(python3 -I -c 'import json, sys
try: print("; ".join(json.load(open(sys.argv[1]))["problems"]))
except Exception: print(open(sys.argv[2]).read().strip()[-300:])' "$SD/verify-checkout.json" "$SD/verify-checkout.err" 2> /dev/null)"
  ft=$(git -C "$SRC" rev-parse "HEAD:scripts/metric-audit/foldgen" 2> /dev/null)
  printf '%s source %s foldgen tree %s pinned %s %s\n' "$(now)" "$C" "${ft:-none}" "$FOLDGEN_TREE" \
    "$([ "$ft" = "$FOLDGEN_TREE" ] && echo ok || echo REFUSED)" > "$SD/foldgen-tree.txt"
  [ "$ft" = "$FOLDGEN_TREE" ] || die 2 "K2: $SRC's foldgen tree is ${ft:-missing}, not the pinned $FOLDGEN_TREE"
  say "K2: verify-checkout of $SRC at $C passed (verify-checkout.json $(sha256sum "$SD/verify-checkout.json" | cut -c1-16)); foldgen tree $ft is the pin"
fi

# ---- 3. inputs
[ "$(sha256sum "$MASTER_MANIFEST" | cut -d' ' -f1)" = "$MASTER_SHA256" ] || die 2 "master manifest is not $MASTER_SHA256"
# an MTB box lists sections, not species: its inputs are MTB's selection plus EXTRA_LISTS (bindings/full-genomes.txt)
MSP=${SPECIES//@canary/}; [ "$KIND" = mtb ] && MSP=mycobacterium_tuberculosis
python3 -I "$KIT/stage-manifest.py" "$DATA_ROOT" "$MASTER_MANIFEST" "$MSP" "${CANARY_REF_PREFIX:-}" "$([ -n "${CANARY_REF_PREFIX:-}" ] && echo 1 || echo 0)" \
  "${EXTRA_LISTS:-}" > "$SD/INPUT-MANIFEST.sha256" || die 2 "stage manifest: $(tail -n 1 "$SD/INPUT-MANIFEST.sha256")"
# a side set (MTB: revision N's g5 data root as the reference; an earlier X box's cohort matrices to restore), each
# its own push-data with its own manifest of allowlisted files
# K3: the reference's report index and dataset reports go too, since plan 7b's fold rebuild (audit.py) reads
# data/report/index.json and data/report/datasets/<id>.json (an X box's restore dirs hold no top-level report/)
# X4: compare_cohorts validates every cohort's report tree against its validation manifest, so each cohort's
# data/report/ JSON (index, datasets, species, intervals) that X2 and X3 built comes back in the restore set too
SIDE_RE='^(cohorts/[a-z0-9]+/data/)?species/[a-z0-9_]+/(genomes\.txt|counted\.json|structure\.json|matrix/[^/]+|results/[^/]+\.json)$|^audit/[A-Za-z0-9][A-Za-z0-9._-]*\.(json|jsonl|tsv|txt)$|^cohorts/[a-z0-9]+/validation-manifest\.json$|^report/index\.json$|^report/datasets/[A-Za-z0-9_.-]+\.json$|^cohorts/[a-z0-9]+/data/report/(index\.json|(datasets|species|intervals)/[A-Za-z0-9_.-]+\.json)$'
side_manifest() { # dir out: species outputs, MTB cohort outputs and run receipts (audit/, validation manifests), the report index and dataset reports, never ksets
  (cd "$1" && find species cohorts audit report -type f 2> /dev/null | grep -v '/ksets/' | grep -E "$SIDE_RE" |
    LC_ALL=C sort | xargs -r sha256sum) > "$2"
  [ -s "$2" ] || die 2 "no files for a side set in $1"
}
# X4 reads X1's matrices and X2's and X3's reports and receipts: RESTORE_DIRS (space list; RESTORE_DIR still works)
# are merged into one hardlinked tree; a path two sources hold must be the same bytes
RESTORE_DIRS="${RESTORE_DIRS:-} ${RESTORE_DIR:-}"
if [ -n "${RESTORE_DIRS// /}" ]; then
  RESTORE_DIR=$SD/restore; rm -rf "$RESTORE_DIR"; mkdir -p "$RESTORE_DIR"
  for d in $RESTORE_DIRS; do
    [ -d "$d" ] || die 2 "restore source $d is not a directory"
    bad=$( (cd "$d" && find species cohorts audit -type f 2> /dev/null | grep -v '/ksets/' | grep -E "$SIDE_RE") | while read -r f; do
      if [ -e "$RESTORE_DIR/$f" ]; then cmp -s "$d/$f" "$RESTORE_DIR/$f" || { echo "disagree $f"; break; }; continue; fi
      mkdir -p "$RESTORE_DIR/$(dirname "$f")" && ln "$d/$f" "$RESTORE_DIR/$f" || { echo "link failed $f"; break; }
    done)
    [ -z "$bad" ] || die 2 "restore merge from $d: $bad"
  done
  say "restore tree from:$RESTORE_DIRS"
fi
[ -n "${REF_DIR:-}" ] && [ ! -d "$REF_DIR" ] && die 2 "reference dir $REF_DIR is not a directory"
[ -n "${REF_DIR:-}" ] && { side_manifest "$REF_DIR" "$SD/REF-MANIFEST.sha256"; say "reference set $REF_DIR pinned: $(sha256sum "$SD/REF-MANIFEST.sha256" | cut -d' ' -f1)"; }
[ -n "${RESTORE_DIR:-}" ] && { side_manifest "$RESTORE_DIR" "$SD/RESTORE-MANIFEST.sha256"; say "restore set $RESTORE_DIR pinned: $(sha256sum "$SD/RESTORE-MANIFEST.sha256" | cut -d' ' -f1)"; }
MSHA=$(sha256sum "$SD/INPUT-MANIFEST.sha256" | cut -d' ' -f1)
say "stage input manifest pinned: $SD/INPUT-MANIFEST.sha256 sha256 $MSHA, $(wc -l < "$SD/INPUT-MANIFEST.sha256") files, from master $MASTER_SHA256"
# species finished in an earlier box of this stage come back as a resume set (their outputs, hashed at pull)
SKIP="" RES=$SD/resume
if [ -d "$SD/out/data/species" ]; then
  SKIP=$(for f in "$SD"/out/run/gates/sp-*.done; do [ -e "$f" ] || continue; s=${f##*/sp-}; echo "${s%.done}"; done | paste -sd,)
fi
START_DONE=$(tr ',' '\n' <<< "$SKIP" | grep -c .)   # species already finished before this box (L5)
if [ -n "$SKIP" ]; then
  rm -rf "$RES"; mkdir -p "$RES/species"
  for s in ${SKIP//,/ }; do cp -al "$SD/out/data/species/$s" "$RES/species/" || die 2 "resume set: no outputs for $s"; done
  (cd "$RES" && find species -type f | LC_ALL=C sort | xargs sha256sum) > "$SD/RESUME-MANIFEST.sha256"
  say "resume set: species $SKIP, manifest sha256 $(sha256sum "$SD/RESUME-MANIFEST.sha256" | cut -d' ' -f1)"
fi

# ---- L5: a stage that hit TTL - margin with nothing finished reruns only with a larger TTL
TTL=$(awk -v e="$EST_FAST" -v m="$MAX_TTL" 'BEGIN { t = int((e + 0.5) * 3600); print (t > m ? m : t) }')
if [ -n "${TTL_OVERRIDE_S:-}" ]; then
  [[ "$TTL_OVERRIDE_S" =~ ^[0-9]+$ ]] && [ "$TTL_OVERRIDE_S" -le "$MAX_TTL" ] || die 2 "TTL_OVERRIDE_S must be seconds, at most $MAX_TTL"
  TTL=$TTL_OVERRIDE_S
fi
if [ -e "$SD/TTL-EXHAUSTED" ]; then
  last_ttl=$(cut -d' ' -f2 "$SD/TTL-EXHAUSTED")
  [ "$TTL" -gt "$last_ttl" ] || die 83 "the last box of $BOX reached its TTL ($last_ttl s) with nothing finished; rerun with TTL_OVERRIDE_S=$MAX_TTL (6 h, PF-D6), not the same TTL"
fi

# ---- 4. plan floor and laptop check
for t in 1 2 3 4; do
  out=$("$R" preflight 2>&1); rc=$?
  [ $rc = 0 ] && break
  [ $rc = 77 ] && die 77 "plan floor: $out"
  [ $rc = 76 ] && [ $t -lt 4 ] && { say "laptop check failed (try $t): $out; waiting 5 min"; sleep 300; continue; }
  die "$rc" "preflight: $out"
done
say "preflight: $out"

# ---- 5. box
CLOSED=0 STARTED=0 FGPID=""
# a long step runs as a waited job, so a signal reaches the exit trap now, not when the step ends
waited() { "$@" & FGPID=$!; wait "$FGPID"; local rc=$?; FGPID=""; return $rc; }
mult() { case "$TYPE" in xlarge) echo 5.56;; large) echo 2.78;; default) echo 1;; small) echo 0.5;; *) echo 5.56;; esac; }
record() { # result: boat usage into the ledger (once), then delete; the box counts as closed only on a confirmed 404
  local u vs xl est="" i gone=0
  trap '' INT TERM HUP                       # N2: nothing interrupts a delete in progress
  [ "$CLOSED" = 1 ] && return 0
  if [ "$BOOKED" != 1 ]; then
    u=$(bt usage "$ID" --json 2> /dev/null)
    vs=$(python3 -I -c 'import json, sys; v = json.load(sys.stdin).get("seconds"); assert isinstance(v, (int, float)) and v > 0; print(int(v))' <<< "$u" 2> /dev/null)
    if [ -z "$vs" ]; then   # L6: book the box's wall time at its type's rate, never 0
      vs=$(awk -v w="$(( $(date +%s) - ${T0:-$(date +%s)} ))" -v m="$(mult)" 'BEGIN { printf "%d", w * m }'); est=" (estimated from wall time: boat usage failed)"
    fi
    xl=$(awk -v s="$vs" -v m="$XL_MULT" 'BEGIN { printf "%.3f", s / m / 3600 }')
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$(now)" "$BOX" "$ID" "${CPU:-?}" "${HOST:-?}" "$xl" "$vs" "$1$est" >> "$LEDGER"
    BOOKED=1 BOOK="$xl xlarge-hours ($vs VM-s)$est"
  fi
  if GRM_BOAT_ADOPT=1 "$R" close "$ID" "$NOTES" > /dev/null 2>&1; then
    gone=1
  else
    # N3: the runner's close refuses when boat usage fails; delete directly and confirm the 404 here
    say "the runner's close refused for $ID; deleting it directly"
    bt delete "$ID" --yes > /dev/null 2>&1
    for i in $(seq 1 "${DELETE_CONFIRM_TRIES:-20}"); do
      if out=$(bt info "$ID" --json 2>&1); then sleep 6; continue; fi
      case "$out" in *[Nn]ot*[Ff]ound*|*404*|*[Nn]o\ such*) gone=1; break;; esac
      sleep 6
    done
  fi
  if [ $gone = 1 ]; then CLOSED=1; say "closed $ID: $BOOK, $1"
  else ALIVE="$ALIVE $ID"; say "BOX STILL UP: $ID could not be deleted (close and direct delete failed); delete it by hand"; fi
}
# N2, M4: whatever ends this script once `new` has been called, the box is stopped, checkpointed if it ran, booked
# and deleted; a second signal during this is ignored; the exit is 86 if a box of this call is still alive
on_exit() {
  local rc=$? ids x cp i nrc_e=""
  trap '' INT TERM HUP
  trap - EXIT
  if [ -n "$FGPID" ] && [ "$NEWPHASE" = 1 ]; then
    # N2: `new` is still running: let it finish (bounded), so the box it makes is recorded and can be adopted
    for i in $(seq 1 420); do kill -0 "$FGPID" 2> /dev/null || break; sleep 1; done
    kill -0 "$FGPID" 2> /dev/null && { pkill -TERM -P "$FGPID" 2> /dev/null; kill -TERM "$FGPID" 2> /dev/null; }
    wait "$FGPID" 2> /dev/null; nrc_e=$?; FGPID=""   # L3: new's own exit code, kept
  fi
  [ -n "$FGPID" ] && { pkill -TERM -P "$FGPID" 2> /dev/null; kill -TERM "$FGPID" 2> /dev/null; wait "$FGPID" 2> /dev/null; FGPID=""; }
  if [ -z "${ID:-}" ] && [ "$NEWCALLED" = 1 ]; then
    case "$nrc_e" in
      2|75|76|77) say "new refused (rc $nrc_e) before it made anything; nothing to look up";;   # L3
      *) adopt_orphans "$rc";;
    esac
  fi
  if [ -n "${ID:-}" ] && [ "$CLOSED" = 0 ]; then
    say "drive exiting (rc $rc) with $ID still up: stopping, pulling what it can, closing"
    if [ "$STARTED" = 1 ] && declare -F ckpt > /dev/null; then
      "$R" run "$ID" --log stop --cwd /home/user -- bash /home/user/port/r11-box.sh stop > /dev/null 2>&1
      # L2: the subshell takes signals again (on_exit ignores them), so the bound below really cuts it
      ( trap - INT TERM HUP; ckpt "abort-$(date +%s)" --all ) & cp=$!
      for i in $(seq 1 "${ABORT_CKPT_SECS:-600}"); do kill -0 "$cp" 2> /dev/null || break; sleep 1; done
      if kill -0 "$cp" 2> /dev/null; then
        kill_tree TERM "$cp"; sleep 2; kill -0 "$cp" 2> /dev/null && kill_tree KILL "$cp"
        say "abort checkpoint cut off after ${ABORT_CKPT_SECS:-600} s"
      fi
      wait "$cp" 2> /dev/null
    fi
    record "aborted rc=$rc"
  fi
  [ -n "$ALIVE" ] && { say "STOP rc=86: box(es) of this call still up:$ALIVE"; exit 86; }
  exit "$rc"
}
# L3: the box name carries this run folder's tag, so no other run folder's box can share it
RUNTAG=$(printf '%s' "$(readlink -f "$RUN_DIR")" | sha256sum | cut -c1-6)
LANE_NAME=grm-$LANE_PREFIX$(tr 'A-Z' 'a-z' <<< "$BOX")-$RUNTAG
kill_tree() { # SIG PID: the process and all its descendants, children first
  local c; for c in $(pgrep -P "$2" 2> /dev/null); do kill_tree "$1" "$c"; done; kill "-$1" "$2" 2> /dev/null; }
# L1: a create can outlive the runner (the runner's timeout runs boat.exe in its own process group): look for this
# call's box while any process for its name is alive (bounded past the runner's 300 s create timeout), then once more;
# a box adopted is deleted with its 404 confirmed by record; anything unsettled is named and the exit is 86
lookup_orphans() { # prints ids; rc 1 when boat list failed
  local lst lrc
  lst=$(bt list --all --json 2> /dev/null); lrc=$?
  { grep -o '^bx_[a-z0-9]*$' "$SD/new.out" 2> /dev/null
    awk -v n="$LANE_NAME" -v t="$T_NEW" '$2 == n && $3 >= t {print $1}' "$BSTATE/created" 2> /dev/null
    [ $lrc = 0 ] && python3 -I -c 'import json, sys
for x in json.load(sys.stdin).get("sandboxes", []):
    if x.get("name") == sys.argv[1] and (x.get("createdAt") or "") >= sys.argv[2]: print(x.get("id"))' "$LANE_NAME" "$T_NEW" <<< "$lst" 2> /dev/null
  } | sort -u
  return $lrc
}
adopt_orphans() { # RC
  local rc=$1 end=$(( T0 + ${NEW_SETTLE_SECS:-330} )) seen=" " ids x inflight last=0 lok
  while :; do
    inflight=$(pgrep -f -- "$LANE_NAME" 2> /dev/null | grep -vx "$$" | head -n 1)
    ids=$(lookup_orphans); lok=$?
    for x in $ids; do
      [[ "$seen" == *" $x "* ]] && continue; seen="$seen$x "
      ID=$x CLOSED=0 BOOKED=0; say "adopting $ID, made by this call's new before it returned (rc $rc)"; record "new interrupted (rc $rc); deleted"
    done
    [ $last = 1 ] && break
    if [ -z "$inflight" ]; then last=1; sleep 3; continue; fi            # nothing in flight: one more look, then stop
    if [ "$(date +%s)" -ge "$end" ]; then
      ALIVE="$ALIVE possible:$LANE_NAME(a create was still running at the ${NEW_SETTLE_SECS:-330} s bound)"; break
    fi
    sleep 2
  done
  [ $lok = 0 ] || ALIVE="$ALIVE possible:$LANE_NAME(boat list failed; check by hand)"
  ID=""
}
NEWCALLED=0 NEWPHASE=0 ALIVE="" BOOKED=0
trap on_exit EXIT; trap 'exit 130' INT; trap 'exit 143' TERM; trap 'exit 129' HUP   # N2: before new
T_NEW=$(date -u -d '-15 sec' +%FT%T) T0=$(date +%s)                                     # I-h: the deadline counts from here
NEWCALLED=1 NEWPHASE=1
waited "$R" new "${LANE_NAME#grm-}" --type "$TYPE" --ttl "$TTL" ${TEMPLATE:+--from-template "$TEMPLATE"} > "$SD/new.out" 2> "$SD/new.err"; nrc=$?
NEWPHASE=0
ID=$(grep -o '^bx_[a-z0-9]*$' "$SD/new.out" | tail -n 1)
if [ $nrc != 0 ] || [ -z "$ID" ]; then
  ID=""
  case $nrc in
    2|75|76|77) NEWCALLED=0   # N4: new refused before it made anything; there is nothing of ours to find
      die "$nrc" "new refused (rc $nrc): $(tail -n 1 "$SD/new.err")";;
  esac
  # M4/N4: new failed after the box may exist: only ids the runner recorded for this exact name since T_NEW
  die "$([ $nrc != 0 ] && echo $nrc || echo 86)" "new failed (rc $nrc): $(tail -n 1 "$SD/new.err")"
fi
DEADLINE=$(( T0 + TTL ))
say "box $ID type $TYPE ttl $TTL s"
bail() { local rc=$1; shift; say "STOP rc=$rc: $*"; record "stopped rc=$rc: $*"; exit "$rc"; }
# a template box holds c90ebaf's tree in /home/user/grm: clear it, so no file of another commit lingers in the build
"$R" run "$ID" --log clear-src --cwd /home/user -- rm -rf /home/user/grm > /dev/null || bail 86 "could not clear /home/user/grm"
"$R" push "$ID" "$SRC" > /dev/null || bail 86 "push source failed"
"$R" push "$ID" "$KIT_REPO" --dest /home/user/port > /dev/null || bail 86 "push kit failed"
ENV="STAGE=$BOX
KIND=$KIND
SPECIES=$SPECIES
REVISION=$REVISION
SOURCE_COMMIT=$(git -C "$SRC" rev-parse HEAD)
MANIFEST_SHA256=$MSHA
CANARY_REF=/home/user/grm-data/${CANARY_REF_PREFIX:-rev9/report}
CHANGED='${CHANGED:-}'
CANARY_SPECIES=$CANARY_SPECIES
BUILD_TARGET_CPU=${BUILD_TARGET_CPU:-}
SKIP_SPECIES=$SKIP
THREADS=${THREADS:-16}
RAM_NEED_GIB=$RAM_NEED_GIB
RAM_PAUSE_GIB=${RAM_PAUSE_GIB:-4}
DISK_PAUSE_GIB=${DISK_PAUSE_GIB:-30}
SLOW_LIB_SECS=${SLOW_LIB_SECS:-60}
RAM_WAIT_SECS=${RAM_WAIT_SECS:-1800}
AUDIT=${AUDIT:-1}
${EXTRA_ENV:-}"
"$R" run "$ID" --log stage-env --cwd /home/user -- "mkdir -p /home/user/r11 && printf '%s\n' $(printf %q "$ENV") > /home/user/r11/stage.env" > /dev/null ||
  bail 86 "could not write stage.env"
waited "$R" run "$ID" --log host --cwd /home/user -- bash /home/user/port/r11-box.sh host > "$SD/host.out" 2>&1 || bail 81 "host step failed: $(tail -n 2 "$SD/host.out" | tr '\n' ' ')"
HOST=$(sed -n 's/^HOST //p' "$SD/host.out" | tail -n 1) CPU=$(sed -n 's/^cpu=//p' "$SD/host.out" | tail -n 1)
say "host $HOST, cpu $CPU, $(grep -h -E '^(speed_lib|avx512f)' "$SD/host.out" | tr '\n' ' ')$(grep -h 'release/grm$' "$SD/host.out")"
if [ "$HOST" = slow ]; then
  if [ "$SLOW_OK" != 1 ]; then bail 80 "slow host for $BOX (slow estimate $EST_SLOW h breaks its ceiling): re-plan or recreate"; fi
  # the box ends at creation + slow estimate + 30 min; boat extend --ttl sets the remaining lifetime (not adds)
  want=$(awk -v s="$EST_SLOW" 'BEGIN { printf "%d", (s + 0.5) * 3600 }')
  [ "$want" -le "$MAX_TTL" ] || bail 80 "slow host: $EST_SLOW h + 30 min passes the $MAX_TTL s TTL ceiling"
  if [ "$want" -gt "$TTL" ]; then
    rem=$(( T0 + want - $(date +%s) ))
    bt extend "$ID" --ttl "$rem" > /dev/null || bail 86 "boat extend failed"
    DEADLINE=$((T0 + want)); say "slow host on a SLOW_OK stage: lifetime set once to $want s from creation ($rem s left)"
  fi
fi

# ---- 5b. the binary gate (one reference per source commit C)
GSHA=$(grep -h 'release/grm$' "$SD/host.out" | tail -n 1 | cut -d' ' -f1)
[[ "$GSHA" =~ ^[0-9a-f]{64}$ ]] || bail 81 "no grm sha256 in the host step's output"
echo "$GSHA $(now) $ID $CPU $C" > "$SD/binary.sha256"
REF_SHA=$(cut -d' ' -f1 "$REFK.sha" 2> /dev/null)
differs() { [ "$GSHA" != "$REF_SHA" ] || [ "${FORCE_BINARY_MISMATCH:-0}" = 1 ]; }
mismatch() { printf '%s %s %s source %s grm %s vs reference %s: %s\n' "$(now)" "$BOX" "$ID" "$C" "$GSHA" "$REF_SHA" "$1" > "$MISK"
  bail 82 "binary gate MISMATCH: $(cat "$MISK"); every box of source $C stops"; }
if [ "$KIND" = mtb ]; then
  if [ -z "$REF_SHA" ] && [ "$BOX" = "$MTB_REF_BOX" ]; then
    printf '%s %s %s %s %s\n' "$GSHA" "$BOX" "$ID" "$(now)" "$CPU" > "$REFK.sha"
    say "MTB reference binary fixed for source $C: grm $GSHA ($CPU); binary only, MTB runs no g1a"
  elif differs; then
    mismatch "an MTB box's binary differs from the MTB reference${FORCE_BINARY_MISMATCH:+ (forced)}"
  else
    say "binary gate: grm $GSHA equals source $C's MTB reference"
  fi
elif [ -n "$REF_BOX" ]; then
  if [ "$BOX" = "$REF_BOX" ] && [ ! -e "$REFK" ]; then
    if [ -n "$REF_SHA" ] && [ "$REF_SHA" != "$GSHA" ]; then
      for s in ${CANARY_SPECIES//,/ }; do [[ ",$SKIP," == *",$s,"* ]] && bail 79 "$REF_BOX resumed on another binary with part of g1a done; rerun $REF_BOX from an empty stage dir"; done
      say "reference binary replaced: the earlier $REF_BOX box made no g1a reports"
    fi
    printf '%s %s %s %s %s\n' "$GSHA" "$BOX" "$ID" "$(now)" "$CPU" > "$REFK.sha"; REF_SHA=$GSHA
    say "reference binary fixed for source $C: grm $GSHA ($CPU); its g1a reports become the reference when they are back"
  elif ! differs; then
    say "binary gate: grm $GSHA equals source $C's reference; no g1a rerun"
  else
    say "binary gate: grm $GSHA differs from the reference $REF_SHA${FORCE_BINARY_MISMATCH:+ (forced)}; g1a reruns here"
    [ -e "$REFK" ] || bail 79 "no reference g1a reports to check against"
    . "$REFK"   # ref_reports_dir ref_reports_manifest ref_reports_manifest_sha256
    python3 -I "$KIT/stage-manifest.py" "$DATA_ROOT" "$MASTER_MANIFEST" "$CANARY_SPECIES" "" 0 > "$SD/CANARY-INPUT-MANIFEST.sha256" ||
      bail 2 "g1a manifest: $(tail -n 1 "$SD/CANARY-INPUT-MANIFEST.sha256")"
    say "g1a input manifest pinned: $(sha256sum "$SD/CANARY-INPUT-MANIFEST.sha256" | cut -d' ' -f1); reference reports $ref_reports_manifest_sha256"
    "$R" push-data "$ID" "$DATA_ROOT" --manifest "$SD/CANARY-INPUT-MANIFEST.sha256" \
      --manifest-sha256 "$(sha256sum "$SD/CANARY-INPUT-MANIFEST.sha256" | cut -d' ' -f1)" --dest /home/user/grm-canary-data > /dev/null 2>&1 ||
      bail 86 "push-data of g1a's inputs failed"
    "$R" push-data "$ID" "$ref_reports_dir" --manifest "$ref_reports_manifest" --manifest-sha256 "$ref_reports_manifest_sha256" \
      --dest /home/user/grm-canary-ref > /dev/null 2>&1 || bail 86 "push-data of the reference g1a reports failed"
    waited "$R" run "$ID" --log bincheck --cwd /home/user -- bash /home/user/port/r11-box.sh bincheck > "$SD/bincheck.out" 2>&1; brc=$?
    "$R" pull "$ID" /home/user/r11/receipts/bincheck "$SD/bincheck-receipts" > /dev/null || say "bincheck receipts pull failed"
    case $brc in
      0) say "binary gate passed: $(grep 'g1a byte check' "$SD/bincheck.out" | tail -n 1)";;
      9) mismatch "$(grep 'g1a byte check' "$SD/bincheck.out" | tail -n 1)";;
      *) bail 81 "bincheck failed (rc $brc): $(tail -n 2 "$SD/bincheck.out" | tr '\n' ' ')";;
    esac
  fi
fi

# ---- 6. data
waited "$R" push-data "$ID" "$DATA_ROOT" --manifest "$SD/INPUT-MANIFEST.sha256" --manifest-sha256 "$MSHA" > "$SD/push-data.out" 2>&1 ||
  bail $? "push-data refused or failed: $(tail -n 1 "$SD/push-data.out")"
for x in REF:/home/user/grm-ref RESTORE:/home/user/grm-resume; do
  v=${x%%:*}_DIR; [ -n "${!v:-}" ] || continue
  m=$SD/${x%%:*}-MANIFEST.sha256
  "$R" push-data "$ID" "${!v}" --manifest "$m" --manifest-sha256 "$(sha256sum "$m" | cut -d' ' -f1)" --dest "${x#*:}" > /dev/null 2>&1 ||
    bail 86 "push-data of the ${x%%:*} set failed"
done
if [ -n "$SKIP" ]; then
  "$R" push-data "$ID" "$RES" --manifest "$SD/RESUME-MANIFEST.sha256" --manifest-sha256 "$(sha256sum "$SD/RESUME-MANIFEST.sha256" | cut -d' ' -f1)" \
    --dest /home/user/grm-resume > /dev/null 2>&1 || bail 86 "resume push-data failed"
fi
"$R" run "$ID" --log verify --cwd /home/user -- bash /home/user/port/r11-box.sh verify > "$SD/verify.out" 2>&1 || bail 81 "verify failed: $(tail -n 2 "$SD/verify.out" | tr '\n' ' ')"
"$R" run "$ID" --log start --cwd /home/user -- bash /home/user/port/r11-box.sh start > /dev/null || bail 86 "start failed"
say "stage started"; STARTED=1
T1=$(date +%s)   # FORCE_TTL_NEAR_AFTER counts from here

# ---- 7. watch
ckpt() { # label [--all]: build, pull, check, merge
  local l=$1 out d
  out=$("$R" run "$ID" --log "ckpt-$l" --cwd /home/user -- bash /home/user/port/ckpt.sh "$l" ${2:-} 2>&1) || { say "checkpoint $l failed on the box: $(tail -n 1 <<< "$out")"; return 1; }
  "$R" pull "$ID" "/home/user/ckpt/$l" "$SD/ckpt" > /dev/null || { say "checkpoint $l pull failed"; return 1; }
  d=$SD/ckpt/$l
  (cd "$d" && sha256sum -c --strict --quiet SHA256SUMS) || { say "checkpoint $l: sha256sum -c FAILED on the laptop"; return 1; }
  mkdir -p "$SD/out"; rsync -a "$d/" "$SD/out/" --exclude SHA256SUMS
  "$R" run "$ID" --log "ckpt-rm-$l" --cwd /home/user -- rm -rf "/home/user/ckpt/$l" > /dev/null
  say "checkpoint $l: $(grep -o 'CKPT.*' <<< "$out" | head -n 1); sha256sum -c --strict ok on the laptop"
  fix_reference
}
# the reference box's g1a reports, once all of g1a is done and back: hashed into RUN_DIR/BINARY-REF-C
fix_reference() {
  local s m=$REFK-REPORTS.sha256
  [ "$KIND" = r11 ] && [ -n "$REF_BOX" ] && [ "$BOX" = "$REF_BOX" ] && [ ! -e "$REFK" ] || return 0
  for s in ${CANARY_SPECIES//,/ }; do [ -e "$SD/out/run/gates/sp-$s.done" ] || return 0; done
  (cd "$SD/out/data" && for s in ${CANARY_SPECIES//,/ }; do find "species/$s/results" -name '*.json' -type f; done | LC_ALL=C sort | xargs -r sha256sum) > "$m"
  [ -s "$m" ] || { say "reference box: g1a done but no reports"; return 0; }
  printf 'ref_sha=%s\nref_reports_dir=%s\nref_reports_manifest=%s\nref_reports_manifest_sha256=%s\nref_source=%s\n' "$(cut -d' ' -f1 "$REFK.sha")" \
    "$SD/out/data" "$m" "$(sha256sum "$m" | cut -d' ' -f1)" "$C" > "$REFK"
  say "reference g1a reports fixed: $(wc -l < "$m") reports, manifest $(sha256sum "$m" | cut -d' ' -f1)"
}
last=$(date +%s) n=0 state=running misses=0
while :; do
  sleep "$POLL" & wait $!
  # polls log to their own notes file; a poll that gets no STATE line is a Boat hiccup, not an ended stage
  st=$("$R" --notes "$SD/poll-notes.md" run "$ID" --log status --cwd /home/user -- bash /home/user/port/r11-box.sh status 2>/dev/null | grep '^STATE' | tail -n 1)
  state=$(awk '{print $2}' <<< "$st")
  if [ -z "$state" ]; then
    misses=$((misses + 1)); [ $misses -ge 5 ] && { state=unknown; say "5 status polls in a row failed"; break; }
    continue
  fi
  misses=0
  [ "$state" = running ] || break
  near=0
  [ $(( DEADLINE - $(date +%s) )) -le "$TTL_MARGIN" ] && near=1
  [ -n "${FORCE_TTL_NEAR_AFTER:-}" ] && [ $(( $(date +%s) - T1 )) -ge "$FORCE_TTL_NEAR_AFTER" ] && near=1
  if [ $near = 1 ]; then
    say "TTL-near ($(( DEADLINE - $(date +%s) )) s left${FORCE_TTL_NEAR_AFTER:+, forced}): stopping for a checkpoint"
    "$R" run "$ID" --log stop --cwd /home/user -- bash /home/user/port/r11-box.sh stop > /dev/null
    n=$((n + 1)); ckpt "ttl-$n" --all || bail 86 "TTL checkpoint failed; the box is closed anyway"
    after=$(ls "$SD"/out/run/gates/sp-*.done 2> /dev/null | wc -l)
    if [ "$KIND" = mtb ] || [ "$after" -le "$START_DONE" ]; then
      printf '%s %s\n' "$(now)" "$TTL" > "$SD/TTL-EXHAUSTED"
      record "TTL reached with nothing finished; rerun with TTL_OVERRIDE_S=$MAX_TTL"
      say "STOP rc=83: $BOX reached TTL - ${TTL_MARGIN} s with nothing finished in this box: rerun with TTL_OVERRIDE_S=$MAX_TTL (6 h, PF-D6)"
      exit 83
    fi
    rm -f "$SD/TTL-EXHAUSTED"
    record "TTL checkpoint; resume with another drive.sh $BOX"
    exit 74
  fi
  if [ $(( $(date +%s) - last )) -ge "$CKPT_SECS" ]; then n=$((n + 1)); ckpt "c$n" || true; last=$(date +%s); fi
done
say "stage ended: $st"

# ---- 8. end
n=$((n + 1)); ckpt "final-$n" --all || bail 86 "final checkpoint failed"
(cd "$SD/out" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum) > "$SD/out/SHA256SUMS"
say "results on the laptop: $SD/out, $(wc -l < "$SD/out/SHA256SUMS") files, SHA256SUMS $(sha256sum "$SD/out/SHA256SUMS" | cut -d' ' -f1)"
if [ -n "${COPY2:-}" ]; then
  mkdir -p "$COPY2/$BOX" && rsync -a "$SD/out/" "$COPY2/$BOX/" && (cd "$COPY2/$BOX" && sha256sum -c --strict --quiet SHA256SUMS) &&
    say "second copy verified: $COPY2/$BOX" || say "second copy FAILED: $COPY2/$BOX"
fi
[ -e "$SD/out/run/receipts/diff-rev9.out" ] && say "revision-9 diff report (not a gate): $(head -n 1 "$SD/out/run/receipts/diff-rev9.out" | cut -c1-240)"
case "$state" in
  done) now > "$SD/DONE"; rm -f "$SD/TTL-EXHAUSTED"; record "done"; exit 0;;
  *) record "stopped: $(cat "$SD/out/run/status" 2> /dev/null)"; exit 81;;
esac
