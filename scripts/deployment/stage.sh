#!/usr/bin/env bash
# Check a candidate on a spare loopback port, never the live one, and approve
# it for activate.sh when every check passes.
# usage: stage.sh LABEL PUBLICATION STAGE_PORT [--apply]
#   GRM_EXPECT_*         what the report must hold (see common.sh); the approval
#                        records it for activation, rollback and restore
#   GRM_PLAYWRIGHT_CORE  playwright-core module directory for the browser check
#                        (web/scripts/browser-check.cjs from src-LABEL); its
#                        Chromium needs LD_LIBRARY_PATH and node on PATH
#   GRM_SKIP_BROWSER=1   approve on the HTTP checks alone; the approval says so
# Receipts go to stage/LABEL--PUBLICATION--TIME.XXXXXX; approval to go/LABEL--PUBLICATION.
# The candidate is stopped on every exit, interruption included; if it cannot
# be stopped, the script says so and fails. A candidate left by a stage that
# was killed outright is found through serve-STAGE_PORT.starting and its spawn
# marker, and --apply stops it before staging again, or refuses, keeping the
# journal, while the marker shows a process it cannot see.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
set -- ${ARGS[@]+"${ARGS[@]}"}
[ $# -eq 3 ] || die "usage: stage.sh LABEL PUBLICATION STAGE_PORT [--apply]" 2
label=$1
pub=$2
sport=$3
[[ $sport =~ ^[0-9]+$ ]] && ((sport > 0 && sport < 65536)) || die "invalid port $sport" 2
[ "$sport" != "$PORT" ] || die "$sport is the live port; stage on a spare one" 2
need curl python3 sha256sum fuser flock setsid pgrep
verify_release "$label"
verify_publication "$pub" "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS"
APPROVED[$label--$pub]="$EXPECT_SCHEMA $EXPECT_REVISION $EXPECT_DATASETS"
R=$D/release/$label
P=$D/published/$pub
commit=$(release_commit "$label")
browser=skipped
if [ "${GRM_SKIP_BROWSER:-0}" != 1 ]; then
    : "${GRM_PLAYWRIGHT_CORE:?set GRM_PLAYWRIGHT_CORE, or GRM_SKIP_BROWSER=1 to approve without the browser check}"
    need node
    [ -f "$D/src-$label/web/scripts/browser-check.cjs" ] || die "no browser check in $D/src-$label" 2
    browser=run
fi
approval=$D/go/$label--$pub
plain_dir "$D/go" && plain_target "$approval" && start_paths_ok "$sport" ||
    die "go/, $approval, serve-$sport.starting, logs/, logs/serve-$sport.log or spawn/ is a link, directory or special file" 2
pending_start "$sport"
if ((PENDING)); then
    say "serve-$sport.starting: an interrupted stage of $J_LABEL with $J_PUB by pid $J_OWNER left: ${ORPHANS:-no process}"
    if ((APPLY)); then resolve_start "$sport"; else would "stop pids ${ORPHANS:-none} and remove the journal"; fi
fi
if ((APPLY)) || ((PENDING == 0)); then
    [ -z "$(listener_pids "$sport")" ] || die "port $sport is in use" 2
fi
say "stage $label ($commit) with $pub (schema $EXPECT_SCHEMA, revision $EXPECT_REVISION, $EXPECT_DATASETS datasets) on port $sport, browser=$browser"
if ! ((APPLY)); then
    would "serve $label with $pub on 127.0.0.1:$sport, run smoke.sh and the browser check, then write $approval"
    exit 0
fi

mkdir -p "$D/stage" "$D/go"
# A unique suffix: a stage rerun within the same second must not collide.
out=$(mktemp -d "$D/stage/$label--$pub--$(date -u +%Y%m%dT%H%M%SZ).XXXXXX") || die "cannot create receipts under $D/stage"
chmod 755 "$out"
run=${out#"$D"/}
started=0
cleanup() {
    local rc=$?
    trap - EXIT
    if ((started)); then
        trap ':' HUP INT TERM
        stop_spawn || die "a staging process of $label on port $sport would not stop; stop it by hand (serve-$sport.starting kept)" 1
        end_journal "$sport"
        ((rc)) || rc=1
        say "stage interrupted or failed (exit $rc); the staging candidate was stopped"
    fi
    exit "$rc"
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
started=1
if ! start_service "$label" "$pub" "$sport"; then
    started=0
    die "cannot write serve-$sport.starting; nothing started" 2
fi
# The candidate is its PID with its start time: once it has exited, its PID may
# belong to another process, which must get no signal.
pid=$START_PID
pstart=$SPAWN_CHILD_START
rc=0
wait_ready "$sport" || rc=1
owns_port "$sport" "$pid" || { echo "listener is not pid $pid" >"$out/listener.log"; rc=1; }
health "$label" "$pub" "$sport" >"$out/smoke.log" 2>&1 || rc=1
checkpoint healthy
if [ "$browser" = run ]; then
    node "$D/src-$label/web/scripts/browser-check.cjs" "$GRM_PLAYWRIGHT_CORE" "$out/browser" \
        "http://127.0.0.1:$sport" "$P/report" >"$out/browser.log" 2>&1 && browser=pass || { browser=fail; rc=1; }
fi
# A candidate that exits during the checks, or ignores SIGTERM (it would also
# stall activation), fails here.
if ! alive "$pid" "$pstart"; then
    echo "pid $pid exited before the checks finished" >"$out/stop.log"
    rc=1
elif ! stop_pid "$pid" "$pstart"; then
    echo "pid $pid ignored SIGTERM" >"$out/stop.log"
    rc=1
    kill_pid "$pid" "$pstart" || die "staging pid $pid on port $sport would not stop; stop it by hand" 1
fi
stop_spawn || die "a staging process of $label on port $sport would not stop; stop it by hand" 1
end_journal "$sport" || die "cannot remove serve-$sport.starting"
started=0
[ -z "$(listener_pids "$sport")" ] || die "port $sport still has a listener after the candidate stopped"
printf '%s\n' "$rc" >"$out/exit"
if ((rc)); then
    die "stage $label with $pub failed; receipts in $out"
fi
write_file "$approval" "$(printf 'stage=0 browser=%s commit=%s publication=%s release_manifest=%s report_manifest=%s schema=%s revision=%s datasets=%s receipts=%s' \
    "$browser" "$commit" "$pub" "$(release_digest "$label")" "$(report_digest "$pub")" \
    "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS" "$run")"$'\n' ||
    die "cannot write $approval; receipts in $out"
check_go "$label" "$pub"
say "stage pass: approved $label with $pub (browser=$browser); receipts in $out"
