#!/usr/bin/env bash
# Replace the live grm serve on GRM_PORT with an approved release and
# publication. The running service is the rollback target, so it must itself
# be an approved, intact pair. Before it stops, it is recorded in
# state/previous and kept in state/active, so restore.sh can bring it back even
# if this script is killed outright. Once it is stopping, any error, failed
# check, INT, TERM or HUP stops the candidate and restarts it.
# usage: activate.sh LABEL PUBLICATION [--apply]
#        activate.sh --rollback LABEL PUBLICATION [--apply]   (run by rollback.sh)
# --rollback requires LABEL PUBLICATION to be state/previous. If the live pair
# then fails its checks (a damaged current release, its executable deleted or
# replaced, say), it is still replaced when it is identified by its recorded
# start (common.sh recorded_service: serve-PORT.pid's PID, start time, token
# and spawn marker) as exactly the pair state/active records; it is neither
# recorded as a rollback target nor restarted, and a failure of the target
# leaves nothing serving (exit 4). A listener identified neither way is never
# stopped.
# Exit codes: 0 switched (or dry run passed); 1 a check of either pair failed
# before any change, or the switch failed or was interrupted and the previous
# pair is serving again; 2 refused before any change (usage, the live service,
# the lock, an unresolved interrupted start, or recording state); 3 the old
# service would not stop (it keeps serving); 4 restoring the previous pair also
# failed, or was not attempted for a damaged pair: the site is down.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
rollback=0
args=()
for arg in "$@"; do
    if [ "$arg" = --rollback ]; then rollback=1; else args+=("$arg"); fi
done
parse_args ${args[@]+"${args[@]}"}
set -- ${ARGS[@]+"${ARGS[@]}"}
[ $# -eq 2 ] || die "usage: activate.sh [--rollback] LABEL PUBLICATION [--apply]" 2
label=$1
pub=$2
need curl python3 sha256sum fuser flock setsid pgrep
note_expect_env

# Identify the one listener on PORT: old_pid, old_start, old_label, old_pub,
# and by_record=1 when only its recorded start identifies it (its executable
# was deleted or replaced; --rollback only). old_mark is the spawn marker it
# holds, removed once it has stopped.
check_live() {
    local pids
    pids=$(listener_pids "$PORT")
    case $(wc -w <<<"$pids") in
        1) ;;
        0) die "nothing listens on $PORT; after a container restart use restore.sh" 2 ;;
        *) die "expected one process on $PORT, got: $pids" 2 ;;
    esac
    old_pid=$pids
    old_start=$(proc_start "$old_pid") || die "pid $old_pid exited during the checks" 2
    by_record=0
    if service_of "$old_pid"; then
        old_label=$SVC_LABEL old_pub=$SVC_PUB
    elif ((rollback)) && recorded_service "$old_pid"; then
        old_label=$RS_LABEL old_pub=$RS_PUB by_record=1 exe_reason=$SVC_ERR
    elif ((rollback)); then
        die "$SVC_ERR, and it is not identified by its recorded start: $RS_ERR; nothing was stopped" 2
    else
        die "$SVC_ERR" 2
    fi
    old_mark=$(held_mark "$old_pid" "$(env_value "$old_pid" GRM_DEPLOY_SPAWN)")
}
# A damaged live pair is replaced only when its recorded start identifies it
# as the pair state/active records.
check_damaged_identity() {
    read_state active
    [ "$ST_LABEL $ST_PUB" = "$old_label $old_pub" ] ||
        die "pid $old_pid serves $old_label with $old_pub, which fails its checks and is not the pair state/active records ($ST_LABEL with $ST_PUB); nothing was stopped" 2
    recorded_service "$old_pid" ||
        die "pid $old_pid serves $old_label with $old_pub, which fails its checks, and is not identified by its recorded start: $RS_ERR; nothing was stopped" 2
    [ "$RS_LABEL $RS_PUB" = "$old_label $old_pub" ] ||
        die "pid $old_pid serves $old_label with $old_pub, but serve-$PORT.pid records $RS_LABEL with $RS_PUB; nothing was stopped" 2
}
paths_ok() {
    start_paths_ok "$PORT" && plain_target "$D/serve-$PORT.pid" ||
        die "serve-$PORT.starting, serve-$PORT.pid, logs/, logs/serve-$PORT.log or spawn/ is a link, directory or special file; nothing was stopped" 2
}

if ((rollback)); then
    verify_state previous
    [ "$ST_LABEL $ST_PUB" = "$label $pub" ] || die "--rollback must name state/previous, $ST_LABEL with $ST_PUB" 2
else
    verify_target "$label" "$pub"
fi
no_pending_start "$PORT"
check_live
[ "$old_label $old_pub" != "$label $pub" ] || die "$label with $pub is already live" 2
# The live pair is the rollback target: before anything stops it must pass every
# check that rollback.sh and restore.sh will apply to it, its approval included.
damaged=0
if ((rollback)) && { ((by_record)) || ! reason=$(verify_target "$old_label" "$old_pub" 2>&1 >/dev/null); }; then
    if ((by_record)); then reason=$exe_reason; fi
    check_damaged_identity
    damaged=1
    say "the live pair $old_label with $old_pub fails its checks (${reason##*FAIL: }); its recorded start identifies pid $old_pid, which will be replaced but not recorded as a rollback target or restarted"
else
    verify_target "$old_label" "$old_pub"
fi
paths_ok
say "live pid=$old_pid $old_label with $old_pub; candidate $label ($(release_commit "$label")) with $pub"
if ! ((APPLY)); then
    if ((damaged)); then
        would "stop pid $old_pid, start $label with $pub on $PORT and health-check it; on failure nothing serves (the damaged pair is not restarted)"
    else
        would "record $old_label/$old_pub in state/previous and state/active, stop pid $old_pid, start $label with $pub on $PORT, health-check it, and restart $old_label on any failure"
    fi
    exit 0
fi

checked="$old_pid $old_start $old_label $old_pub $by_record"
lock_root
no_pending_start "$PORT"
check_live
[ "$old_pid $old_start $old_label $old_pub $by_record" = "$checked" ] || die "the live service changed during the checks" 2
paths_ok
if ((damaged)); then
    check_damaged_identity
else
    write_state previous "$old_label" "$old_pub" && write_state active "$old_label" "$old_pub" ||
        die "cannot record state/previous and state/active; nothing was stopped" 2
fi

# From here every exit runs recover, which restarts the old pair unless the
# switch completed (phase=done) or the old service never stopped (phase=kept).
recover() {
    local rc=$? rb_pid
    trap - EXIT
    [ "$phase" = switching ] || exit "$rc"
    # Handled, not ignored: an ignored signal would stay ignored in the restarted grm.
    trap ':' HUP INT TERM
    set +e
    say "activation of $label stopped (exit $rc); restoring $old_label with $old_pub"
    # The candidate, found by its spawn marker even if START_PID was never set.
    stop_spawn || die "a candidate process of $label would not stop; nothing restarted, serve-$PORT.starting kept for restore.sh" 4
    end_journal "$PORT"
    if alive "$old_pid" "$old_start"; then
        if owns_port "$PORT" "$old_pid" && { ((damaged)) || health "$old_label" "$old_pub" "$PORT" >/dev/null 2>&1; }; then
            die "activation interrupted before pid $old_pid stopped; $old_label with $old_pub is still serving" 1
        fi
        SIG_WHY="the live server being replaced on port $PORT, pid $old_pid with start $old_start, $old_label with $old_pub"
        ensure_stopped "$old_pid" "$old_start" || die "pid $old_pid neither serves nor stops; nothing healthy on $PORT" 4
    fi
    if ((damaged)); then
        die "rollback to $label with $pub failed; $old_label with $old_pub failed its checks and was not restarted, so nothing serves $PORT; start an intact approved pair with restore.sh LABEL PUBLICATION --apply" 4
    fi
    start_service "$old_label" "$old_pub" "$PORT" force
    rb_pid=$START_PID
    checkpoint restarted
    say "restarted $old_label pid=$rb_pid"
    if serving "$old_label" "$old_pub" "$rb_pid"; then
        forget_mark "$old_mark"
        write_state active "$old_label" "$old_pub" ||
            say "warning: could not rewrite state/active; record the serving pair with restore.sh $old_label $old_pub --apply"
        record_pid "$PORT" "$old_label" "$old_pub" "$rb_pid" || say "warning: could not write serve-$PORT.pid"
        end_journal "$PORT"
        die "activation failed; $old_label with $old_pub is serving again (pid $rb_pid)" 1
    fi
    stop_spawn || die "restarted pid $rb_pid failed and would not stop; serve-$PORT.starting kept for restore.sh" 4
    end_journal "$PORT"
    die "activation and rollback failed; nothing healthy on $PORT, see $D/deploy.log" 4
}
phase=switching
trap recover EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

say "stopping pid=$old_pid"
SIG_WHY="the live server being replaced on port $PORT, pid $old_pid with start $old_start, $old_label with $old_pub"
if ! stop_pid "$old_pid" "$old_start"; then
    if owns_port "$PORT" "$old_pid"; then
        phase=kept
        die "pid $old_pid did not stop; it is still serving, nothing started" 3
    fi
    die "pid $old_pid did not stop and no longer listens"
fi
checkpoint old-stopped
start_service "$label" "$pub" || die "cannot write serve-$PORT.starting; the candidate was not started"
say "started $label pid=$START_PID"
serving "$label" "$pub" "$START_PID" || die "live health FAIL for $label with $pub"
checkpoint healthy
write_state active "$label" "$pub" || die "cannot record $label with $pub in state/active"
checkpoint committed
record_pid "$PORT" "$label" "$pub" "$START_PID" ||
    say "warning: could not write serve-$PORT.pid; $label with $pub serves and state/active records it"
checkpoint recorded
end_journal "$PORT" || say "warning: could not remove serve-$PORT.starting"
phase=done
forget_mark "$old_mark"
if ((damaged)); then
    say "live health pass: $label with $pub; $old_label with $old_pub failed its checks and is not a rollback target"
else
    say "live health pass: $label with $pub; rollback target $old_label with $old_pub"
fi
