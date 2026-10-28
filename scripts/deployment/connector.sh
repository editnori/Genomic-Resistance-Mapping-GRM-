#!/usr/bin/env bash
# Cloudflare Tunnel connector for the GRM origin. The tunnel's public hostname
# route (configured in Cloudflare, not here) points at http://127.0.0.1:8641.
# This script never reads, prints or copies the tunnel token: it checks the
# token file's ownership and mode, and cloudflared reads it via --token-file.
# usage: connector.sh check|run|stop
#   check  verify the binary against its SHA-256 and the token file's
#          permissions; changes nothing
#   run    supervise cloudflared in the foreground, restarting it after 5 s and
#          checking the binary again before every start; a second supervisor on
#          the same directory exits 0. Writes supervisor.pid, connector.pid,
#          connector.lock, connector.spawn and logs/ in the connector directory
#   stop   SIGTERM the supervisor recorded in supervisor.pid, only if that PID
#          still has the recorded start time, works in this connector
#          directory, and is the process whose descriptor holds the flock on
#          its connector.lock; then stop any cloudflared a killed supervisor
#          left (see connector.spawn below)
#   GRM_CONNECTOR_DIR      root-only directory with cloudflared,
#                          cloudflared.sha256 and logs/ (required, absolute)
#   GRM_TUNNEL_TOKEN_FILE  owned by the caller, mode 0600 or 0400, absolute
#                          (default GRM_CONNECTOR_DIR/tunnel.token)
#   GRM_TUNNEL_PROTOCOL    http2 (default; QUIC failed from the Vast host) or quic
# cloudflared.sha256 must hold exactly one line, "SHA256  cloudflared" (or
# "SHA256 *./cloudflared"), ending in a newline (LF, not CRLF), and no other byte.
#
# Every file here is opened only through lib.sh (open_fd, safe_line,
# safe_sha256, write_file), so a link, FIFO or other special file at
# connector.lock, connector.spawn, supervisor.pid, connector.pid, logs/ or in
# logs/ is refused and never blocks. The supervisor holds a flock on
# connector.spawn and forks cloudflared with it open, so cloudflared, from
# the fork on, holds it too: a cloudflared left by a supervisor that was
# killed outright still holds that lock, and the next run or stop finds it by
# the file's device and inode, stops it with its process group (cloudflared
# runs under setsid), and starts nothing while the lock shows a process it
# cannot see. A forked child execs cloudflared only while its supervisor is
# alive, and the supervisor's signal handler stops the child it last forked
# even before the PID is assigned, so a TERM at any point leaves no
# cloudflared and no supervisor waiting for one.
#
# Each run gives cloudflared a fresh random token in its environment
# (GRM_DEPLOY_SPAWN), which whatever it forks inherits. Once cloudflared has
# exited, what it left in its process group is stopped only when each member
# carries that token, holds connector.spawn or was seen in the group while
# cloudflared lived (lib.sh stop_tree); a group number can be reused.
#
# Output is best effort (say): a closed output pipe makes a write fail, not
# kill the supervisor (SIGPIPE is handled, so cloudflared still starts with
# the default disposition), so cleanup always runs.
set -euo pipefail
trap ':' PIPE
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/lib.sh"
say() { printf '%s\n' "$*" 2>/dev/null || true; }
die() {
    say "FAIL: $1" >&2
    exit "${2:-1}"
}
C=${GRM_CONNECTOR_DIR:?set GRM_CONNECTOR_DIR to the connector directory}
C=${C%/}
TOKEN=${GRM_TUNNEL_TOKEN_FILE:-$C/tunnel.token}
PROTOCOL=${GRM_TUNNEL_PROTOCOL:-http2}
LOCK_FD='' MARK_FD='' MARK_ID='' SPAWN_TOKEN=''
# run changes into the connector directory before cloudflared opens its files,
# so a relative path would be checked in one directory and used in another.
case $C in /?*) ;; *) die "GRM_CONNECTOR_DIR must be an absolute path: $C" 2 ;; esac
need() {
    local tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null || die "missing required tool $tool" 2
    done
}

private() {  # owned by us, no group or other access
    [ "$(stat -c %u "$1")" = "$(id -u)" ] || die "$1 is not owned by uid $(id -u)" 2
    case $(stat -c %a "$1") in
        [4-7]00) ;;
        *) die "$1 has mode $(stat -c %a "$1"); it must not be group- or world-accessible" 2 ;;
    esac
}

# The manifest is exactly one complete entry (lib.sh safe_line), for the
# executable that runs, whose bytes are hashed without following a link.
verify_binary() {
    local line want got
    [ -f "$C/cloudflared" ] && [ ! -L "$C/cloudflared" ] && [ -x "$C/cloudflared" ] ||
        die "$C/cloudflared is missing, not executable or a link" 2
    [ -f "$C/cloudflared.sha256" ] && [ ! -L "$C/cloudflared.sha256" ] || die "$C/cloudflared.sha256 is missing" 2
    line=$(safe_line "$C/cloudflared.sha256" 4096) && [[ $line =~ ^([0-9a-f]{64})\ [\ *](\./)?cloudflared$ ]] ||
        die "$C/cloudflared.sha256 must hold exactly one newline-terminated entry, for cloudflared, and no other byte" 2
    want=${BASH_REMATCH[1]}
    got=$(safe_sha256 "$C/cloudflared") || die "$C/cloudflared cannot be read as a regular file" 2
    [ "$got" = "$want" ] || die "cloudflared does not match its SHA-256" 2
}

check() {
    local version
    case $PROTOCOL in http2 | quic) ;; *) die "unsupported protocol $PROTOCOL" 2 ;; esac
    case $TOKEN in /*) ;; *) die "GRM_TUNNEL_TOKEN_FILE must be an absolute path: $TOKEN" 2 ;; esac
    [ -d "$C" ] && [ ! -L "$C" ] || die "$C is not a directory" 2
    private "$C"
    verify_binary
    [ -f "$TOKEN" ] && [ ! -L "$TOKEN" ] && [ -s "$TOKEN" ] || die "token file $TOKEN is missing, empty or a link" 2
    private "$TOKEN"
    private "$(dirname -- "$TOKEN")"
    version=$("$C/cloudflared" --version 2>/dev/null </dev/null | sed -n 1p) || true
    say "connector check pass: $version, protocol=$PROTOCOL, token file permissions ok"
}

# The files run writes or locks are plain, and logs/ holds only regular files.
paths_ok() {
    local f odd=''
    for f in connector.lock connector.spawn supervisor.pid connector.pid; do
        plain_target "$C/$f" || return 1
    done
    plain_dir "$C/logs" || return 1
    if [ -d "$C/logs" ]; then odd=$(find "$C/logs" -mindepth 1 -maxdepth 1 ! -type f -print -quit); fi
    [ -z "$odd" ]
}

# No supervisor holds connector.lock (it is absent, or its lock is free).
supervisor_absent() {
    local rc=0
    [ -e "$C/connector.lock" ] || [ -L "$C/connector.lock" ] || return 0
    lock_free "$C/connector.lock" || rc=$?
    case $rc in
        0) return 0 ;;
        1) return 1 ;;
        *) die "$C/connector.lock is a link, directory, FIFO or other special file" 2 ;;
    esac
}
# Stop the processes holding connector.spawn, the ones a supervisor killed
# outright left running; the caller holds connector.lock or has seen it free.
# Returns 1, so nothing starts, while any may remain.
reap_orphans() {
    local id held
    [ -e "$C/connector.spawn" ] || [ -L "$C/connector.spawn" ] || return 0
    id=$(file_id "$C/connector.spawn") || { say "FAIL: $C/connector.spawn is not a regular file" >&2; return 1; }
    lock_free "$C/connector.spawn" && return 0
    held=$(mark_holders "$id")
    say "connector processes left by a supervisor that was killed: ${held:-none this user can see}"
    mark_stop "$id" "$C/connector.spawn" || { say "FAIL: they may remain: $MARK_ERR" >&2; return 1; }
    say "stopped the connector processes left by a killed supervisor: ${held:-none}"
}

# Called on TERM, INT and HUP: stop at the next check, and stop the child
# forked last, found by $! (set by the fork itself) until child is assigned.
on_signal() {
    local p=${child:-} s=${child_start:-}
    stopping=1
    if [ -z "$p" ]; then
        p=${!:-}
        [ -n "$p" ] && [ "$(ppid_of "$p")" = "$$" ] && s=$(proc_start "$p") || return 0
    fi
    [ -n "$s" ] && alive "$p" "$s" || return 0
    if [ "$PS_PGRP" = "$p" ]; then tree_signal TERM "$p" "$p"; else tree_signal TERM "$p" ''; fi
    tree_signal CONT "$p" ''
}
on_exit() {
    local rc=$?
    trap - EXIT
    trap ':' TERM INT HUP
    set +e
    # With its process group, even if cloudflared itself has exited.
    if [ -n "${child:-}" ] && [ -n "${child_start:-}" ]; then stop_tree "$child" "$child_start" "$SPAWN_TOKEN" "$MARK_ID" group || rc=1; fi
    close_fd MARK_FD
    if [ -n "$MARK_ID" ] && ! mark_stop "$MARK_ID" "$C/connector.spawn" '' '' "$SPAWN_TOKEN"; then
        say "FAIL: connector processes may remain: $MARK_ERR" >&2
        rc=1
    fi
    rm -f -- "$C/supervisor.pid" "$C/connector.pid"
    exit "$rc"
}

run() {
    local me i
    need flock setsid pgrep python3 sha256sum
    check
    paths_ok || die "connector.lock, connector.spawn, supervisor.pid, connector.pid, logs/ or an entry of logs/ is a link, directory, FIFO or other special file" 2
    [ -d "$C/logs" ] || mkdir "$C/logs"
    cd "$C"
    umask 077
    child='' child_start='' stopping=0
    SPAWN_TOKEN=$(cat /proc/sys/kernel/random/uuid)
    trap on_signal TERM INT HUP
    open_fd LOCK_FD "$C/connector.lock" || die "$C/connector.lock cannot be opened as a regular file" 2
    flock -n "$LOCK_FD" || { say "a connector supervisor already holds $C/connector.lock"; exit 0; }
    checkpoint conn-locked
    reap_orphans || die "not starting while processes of a killed supervisor may remain" 1
    mark_lock "$C/connector.spawn" || die "cannot lock $C/connector.spawn" 1
    trap on_exit EXIT
    me=$(proc_start "$$")
    write_file "$C/supervisor.pid" "$$ $me"$'\n' || die "cannot write $C/supervisor.pid" 2
    checkpoint conn-recorded
    while ((!stopping)); do
        verify_binary
        checkpoint conn-verified
        ((!stopping)) || break
        # The child holds connector.spawn (MARK_FD) from the fork on, but not
        # connector.lock, which would let it outlive the supervisor holding it.
        (
            exec {LOCK_FD}>&- >/dev/null 2>&1 </dev/null
            child_checkpoint conn-child-forked
            alive "$$" "$me" || exit 97
            child_checkpoint conn-child-preexec
            alive "$$" "$me" || exit 97
            GRM_DEPLOY_SPAWN=$SPAWN_TOKEN exec setsid ./cloudflared tunnel --no-autoupdate --protocol "$PROTOCOL" --metrics 127.0.0.1:0 \
                --log-directory "$C/logs" --pidfile "$C/connector.pid" run --token-file "$TOKEN"
        ) &
        checkpoint conn-forked
        child=$!
        child_start=$(proc_start "$child") || child_start=''
        checkpoint conn-assigned
        while [ -n "$child_start" ] && alive "$child" "$child_start"; do
            if ((stopping)); then
                stop_tree "$child" "$child_start" "$SPAWN_TOKEN" "$MARK_ID" || die "cloudflared pid $child would not stop: $STOP_ERR" 1
                break
            fi
            wait "$child" 2>/dev/null || true
        done
        # cloudflared exited: nothing it forked into its group may outlive it.
        if [ -n "$child_start" ]; then
            stop_tree "$child" "$child_start" "$SPAWN_TOKEN" "$MARK_ID" group || die "processes of cloudflared pid $child would not stop: $STOP_ERR" 1
        fi
        child='' child_start=''
        for ((i = 0; i < 25 && !stopping; i++)); do sleep 0.2 {LOCK_FD}>&-; done
    done
}

# Whether PID (started at START) supervises this directory: it works here,
# and one of its descriptors is the open file description that holds the
# flock on this directory's connector.lock, which must be locked. A process
# that merely has connector.lock open, a command line naming this script, or
# a record naming it is not enough. Sets SUP_ERR.
supervisor_identity() {
    local pid=$1 start=$2 here=$3 lock_id fd rc
    SUP_ERR=''
    lock_id=$(file_id "$C/connector.lock") || { SUP_ERR="$C/connector.lock is missing or not a regular file"; return 1; }
    rc=0
    lock_free "$C/connector.lock" || rc=$?
    case $rc in
        0) SUP_ERR="no process holds the lock on $C/connector.lock"; return 1 ;;
        1) ;;
        *) SUP_ERR="$C/connector.lock cannot be checked"; return 1 ;;
    esac
    [ "$(readlink "/proc/$pid/cwd" 2>/dev/null)" = "$here" ] || { SUP_ERR="pid $pid works in another directory"; return 1; }
    for fd in $(fds_on "$pid" "$lock_id"); do
        if fd_holds_flock "$pid" "$fd"; then
            alive "$pid" "$start" || { SUP_ERR="pid $pid exited"; return 1; }
            return 0
        fi
    done
    SUP_ERR="pid $pid does not hold the lock on $C/connector.lock"
    return 1
}

stop() {
    local pid start line here _
    need flock pgrep python3
    [ -d "$C" ] && [ ! -L "$C" ] || die "$C is not a directory" 2
    here=$(cd -- "$C" && pwd -P)
    if [ ! -e "$C/supervisor.pid" ] && [ ! -L "$C/supervisor.pid" ]; then
        say "connector supervisor is not running"
        if supervisor_absent; then reap_orphans || die "connector processes of a killed supervisor may remain" 1; fi
        return 0
    fi
    line=$(safe_line "$C/supervisor.pid" 64) && [[ $line =~ ^([0-9]+)\ ([0-9]+)$ ]] ||
        die "$C/supervisor.pid is malformed; no signal sent" 2
    pid=${BASH_REMATCH[1]}
    start=${BASH_REMATCH[2]}
    if ! alive "$pid" "$start"; then
        say "supervisor.pid is stale: pid $pid has exited or is another process now; no signal sent"
        if supervisor_absent; then reap_orphans || die "connector processes of a killed supervisor may remain" 1; fi
        return 0
    fi
    supervisor_identity "$pid" "$start" "$here" ||
        die "pid $pid is not the supervisor of this connector directory: $SUP_ERR; no signal sent" 2
    kill -TERM "$pid"
    for _ in $(seq 1 90); do
        if ! alive "$pid" "$start"; then
            if supervisor_absent; then reap_orphans || die "the supervisor stopped, but connector processes may remain" 1; fi
            say "connector stopped; the grm service is unchanged"
            return 0
        fi
        sleep 0.5
    done
    die "supervisor $pid did not stop within 45 s; inspect $C/logs"
}

case ${1:-} in
    check) check ;;
    run) run ;;
    stop) stop ;;
    *) die "usage: connector.sh check|run|stop" 2 ;;
esac
