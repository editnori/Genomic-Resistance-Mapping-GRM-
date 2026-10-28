#!/usr/bin/env bash
# Report what serves GRM_PORT, the recorded state, any interrupted start, and
# the live health check (smoke.sh) against the live pair's approval. It does
# not verify the release or publication manifests; activate.sh and restore.sh
# do. It opens no link or special file as a record, so a FIFO cannot stall it.
# Read-only: it writes only smoke.sh's temporary directory. usage: status.sh
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
need curl python3 sha256sum fuser
show() {  # show NAME FILE FORM
    if [ ! -e "$2" ] && [ ! -L "$2" ]; then
        printf '%s: none\n' "$1"
    elif parse_fields "$2" "$3"; then
        printf '%s: %s\n' "$1" "$(for key in "${!FIELDS[@]}"; do printf '%s=%s ' "$key" "${FIELDS[$key]}"; done | tr ' ' '\n' | LC_ALL=C sort | xargs)"
    else
        printf '%s: unreadable (it %s)\n' "$1" "$PF_ERR"
    fi
}
show state/active "$D/state/active" lines
show state/previous "$D/state/previous" lines
show "serve-$PORT.pid" "$D/serve-$PORT.pid" line
show "serve-$PORT.starting" "$D/serve-$PORT.starting" line
if [ -e "$D/serve-$PORT.starting" ] && read_journal "$PORT"; then
    if alive "$J_OWNER" "$J_OWNER_START"; then
        say "pid $J_OWNER is starting $J_LABEL with $J_PUB"
    else
        orphans=$(journal_procs)
        say "an interrupted start of $J_LABEL with $J_PUB left pids: ${orphans:-none}; restore.sh --apply resolves it"
    fi
fi
pids=$(listener_pids "$PORT")
case $(wc -w <<<"$pids") in
    0) die "nothing listens on $PORT" 1 ;;
    1) ;;
    *) die "expected one process on $PORT, got: $pids" 1 ;;
esac
if service_of "$pids"; then
    label=$SVC_LABEL pub=$SVC_PUB
elif recorded_service "$pids"; then
    label=$RS_LABEL pub=$RS_PUB
    say "pid $pids is identified only by its recorded start: $SVC_ERR"
else
    die "$SVC_ERR" 1
fi
say "live pid=$pids release=$label commit=$(release_commit "$label") publication=$pub"
if (check_go "$label" "$pub") >/dev/null 2>&1; then
    check_go "$label" "$pub"
    if state_is active "$label" "$pub"; then
        say "state/active records the live pair"
    else
        say "state/active does not record the live pair; see restore.sh"
    fi
else
    say "the live pair has no valid approval; checking it against GRM_EXPECT_* (schema $EXPECT_SCHEMA, revision $EXPECT_REVISION, $EXPECT_DATASETS datasets)"
    APPROVED[$label--$pub]="$EXPECT_SCHEMA $EXPECT_REVISION $EXPECT_DATASETS"
fi
health "$label" "$pub" "$PORT"
