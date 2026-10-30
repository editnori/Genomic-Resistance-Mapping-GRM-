#!/usr/bin/env bash
# Start the recorded service when nothing serves GRM_PORT, e.g. after a Vast
# container restart. The target is state/active unless LABEL PUBLICATION are
# given (needed once for a root that activate.sh has not managed yet, or when
# state/active names a pair that no longer verifies).
# usage: restore.sh [LABEL PUBLICATION] [--apply]
# It first resolves an interrupted start on the port: if serve-PORT.starting
# names a script that is gone, --apply stops every process of that start (the
# holders of its spawn marker and its registered child), unless
# serve-PORT.pid shows it committed, and removes the journal only once the
# marker's lock shows that none is left; otherwise it keeps the journal and
# fails (common.sh resolve_start). When the target already serves and is healthy, nothing restarts; if
# state/active or serve-PORT.pid does not record it yet, --apply records it, so
# a later restore.sh without arguments finds it. Refuses when anything else
# holds the port. A failed or interrupted start stops what it started, leaving
# nothing serving and state/active as it was.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
set -- ${ARGS[@]+"${ARGS[@]}"}
need curl python3 sha256sum fuser flock setsid pgrep
case $# in
    0) verify_state active; label=$ST_LABEL; pub=$ST_PUB ;;
    2) label=$1; pub=$2; verify_target "$label" "$pub" ;;
    *) die "usage: restore.sh [LABEL PUBLICATION] [--apply]" 2 ;;
esac
note_expect_env

# Set live_pid to the healthy target serving PORT, or '' when nothing listens;
# refuse anything else.
check_port() {
    local pids
    live_pid=''
    pids=$(listener_pids "$PORT")
    case $(wc -w <<<"$pids") in
        0) return 0 ;;
        1) ;;
        *) die "expected at most one process on $PORT, got: $pids" 2 ;;
    esac
    service_of "$pids" || die "$SVC_ERR" 2
    [ "$SVC_LABEL $SVC_PUB" = "$label $pub" ] || die "port $PORT already serves $SVC_LABEL $SVC_PUB; use activate.sh to switch" 2
    health "$label" "$pub" "$PORT" >/dev/null || die "pid $pids serves $label $pub but fails its health check" 1
    live_pid=$pids
}
# state/active records the pair, and serve-PORT.pid holds exactly the record
# record_pid would write for PID now (its start time, token and marker).
recorded() {
    local want line
    state_is active "$label" "$pub" && want=$(pid_line "$label" "$pub" "$1") &&
        line=$(safe_line "$D/serve-$PORT.pid" 4096) && [ "$line" = "$want" ]
}

pending_start "$PORT"
if ((PENDING)); then
    say "serve-$PORT.starting: a start of $J_LABEL with $J_PUB by pid $J_OWNER was interrupted; its processes: ${ORPHANS:-none}"
    if ! ((APPLY)); then
        would "stop pids ${ORPHANS:-none} unless serve-$PORT.pid shows that start committed, remove the journal, then restore $label with $pub"
        would_stop_orphans
        exit 0
    fi
    lock_root
    pending_start "$PORT"
    if ((PENDING)); then resolve_start "$PORT"; fi
fi

check_port
if [ -n "$live_pid" ]; then
    if recorded "$live_pid"; then
        say "$label with $pub is already serving on $PORT (pid $live_pid), healthy and recorded in state/active"
        exit 0
    fi
    say "$label with $pub is already serving on $PORT (pid $live_pid) and healthy, but state/active or serve-$PORT.pid does not record it"
    if ! ((APPLY)); then
        would "record $label with $pub in state/active and pid $live_pid in serve-$PORT.pid"
        exit 0
    fi
    checked=$live_pid
    lock_root
    no_pending_start "$PORT"
    check_port
    [ "$live_pid" = "$checked" ] || die "the service on $PORT changed during the checks" 2
    write_state active "$label" "$pub" && record_pid "$PORT" "$label" "$pub" "$live_pid" ||
        die "cannot record $label with $pub in state/active and serve-$PORT.pid; the service was not touched"
    say "recorded $label with $pub (pid $live_pid) in state/active; nothing restarted"
    exit 0
fi

paths_ok() {
    start_paths_ok "$PORT" && plain_target "$D/serve-$PORT.pid" && plain_dir "$D/state" && plain_target "$D/state/active" ||
        die "serve-$PORT.starting, serve-$PORT.pid, logs/, logs/serve-$PORT.log, spawn/, state/ or state/active is a link, directory or special file; nothing started" 2
}
paths_ok
say "nothing listens on $PORT; target $label ($(release_commit "$label")) with $pub"
if ! ((APPLY)); then
    would "start $label with $pub on 127.0.0.1:$PORT, health-check it and record it in state/active"
    exit 0
fi

lock_root
no_pending_start "$PORT"
check_port
[ -z "$live_pid" ] || die "pid $live_pid started serving $PORT during the checks; run restore.sh again" 2
paths_ok
# A failed or interrupted restore stops what it started, found by its spawn
# marker even if START_PID was never set.
started=0
cleanup() {
    local rc=$?
    trap - EXIT
    ((started)) || exit "$rc"
    trap ':' HUP INT TERM
    stop_spawn || die "a process started for $label would not stop; stop it by hand (serve-$PORT.starting kept)" 1
    end_journal "$PORT"
    die "restored service failed its health check or the restore was interrupted (exit $rc); what it started was stopped, see $D/deploy.log" 1
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
started=1
if ! start_service "$label" "$pub"; then
    started=0
    die "cannot write serve-$PORT.starting; nothing started" 2
fi
say "started $label pid=$START_PID"
serving "$label" "$pub" "$START_PID" || exit 1
checkpoint healthy
# A full disk must not take down a service state/active already names.
write_state active "$label" "$pub" || state_is active "$label" "$pub" || exit 1
checkpoint committed
record_pid "$PORT" "$label" "$pub" "$START_PID" ||
    say "warning: could not write serve-$PORT.pid; $label with $pub serves and state/active records it; run restore.sh --apply again to record the pid"
checkpoint recorded
end_journal "$PORT" || say "warning: could not remove serve-$PORT.starting"
started=0
say "restore health pass: $label with $pub on $PORT"
