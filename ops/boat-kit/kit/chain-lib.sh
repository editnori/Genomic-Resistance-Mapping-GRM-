# chain-lib.sh: sourced by a box chain script (one detached script per box, as the lane rules ask). REVIEW-DEPS L1:
# DEPS's chain always exited 0 (a failed step didn't propagate and the last command's rc won), so the runner's
# "rc=0" for it was no gate. Here every step's rc is kept and the chain's exit says whether all passed.
#   CHAIN_LOG=DIR (required): STEPS and one NAME.log per step
#   step NAME CMD...   run CMD (output to NAME.log), record "NAME rc=N secs" in STEPS; returns CMD's rc
#   must NAME CMD...   a step the rest depends on: on failure, finish at once
#   finish             record the summary line and exit: 0 when every step passed, else 1 (never the last rc); a
#                      chain that ran no step fails (REVIEW-KIT10 note: nothing ran, so nothing passed)
: "${CHAIN_LOG:?set CHAIN_LOG before sourcing chain-lib.sh}"
mkdir -p "$CHAIN_LOG"
CHAIN_FAILED="" CHAIN_N=0
step() {
  local name=$1 t0 rc; shift
  [[ "$name" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "chain: bad step name $name" >&2; CHAIN_FAILED="$CHAIN_FAILED badname"; return 2; }
  echo "$(date -u +%FT%TZ) $name start" >> "$CHAIN_LOG/STEPS"; t0=$(date +%s)
  "$@" > "$CHAIN_LOG/$name.log" 2>&1; rc=$?
  echo "$(date -u +%FT%TZ) $name rc=$rc $(( $(date +%s) - t0 ))s" >> "$CHAIN_LOG/STEPS"
  CHAIN_N=$((CHAIN_N + 1))
  [ $rc = 0 ] || CHAIN_FAILED="$CHAIN_FAILED $name"
  return $rc
}
must() { step "$@" || { echo "$(date -u +%FT%TZ) chain stopped after $1" >> "$CHAIN_LOG/STEPS"; finish; }; }
finish() {
  if [ "$CHAIN_N" = 0 ]; then
    echo "$(date -u +%FT%TZ) CHAIN FAIL: no step ran" >> "$CHAIN_LOG/STEPS"
    exit 1
  fi
  if [ -z "$CHAIN_FAILED" ]; then
    echo "$(date -u +%FT%TZ) CHAIN PASS $CHAIN_N steps" >> "$CHAIN_LOG/STEPS"; exit 0
  fi
  echo "$(date -u +%FT%TZ) CHAIN FAIL:$CHAIN_FAILED (of $CHAIN_N steps)" >> "$CHAIN_LOG/STEPS"; exit 1
}
