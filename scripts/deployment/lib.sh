# shellcheck shell=bash
# Primitives shared by common.sh (and so every deployment script), connector.sh
# and smoke.sh: safe file access, process identity, and the supervision of
# processes they start. Source it; do not run it. It sets no options and reads
# no configuration.
#
# Safe file access. No record, lock, log, pid or marker file is opened when it
# is a link, FIFO, socket, device or directory, and no open can block:
#   - reads (safe_line, safe_sha256, and parse_fields in common.sh) go through
#     GRM_PY, which opens with O_NOFOLLOW and O_NONBLOCK and checks the
#     descriptor it got;
#   - descriptors bash itself holds (locks, logs, spawn markers) come from
#     open_fd or open_append, which check the path, open it read-write (so a
#     FIFO swapped in after the check cannot block the open) and then check
#     the descriptor: it must be a regular file at that very path;
#   - whole files are written by write_file: a new temporary file, renamed
#     over the target and read back.
#
# Process identity is a PID plus its start time (proc_start), never a PID
# alone, a path or a command line.
#
# Supervised processes. A script that starts a long-lived process first takes
# a flock on a marker file and forks with that descriptor open, so the child
# holds the lock from the instant it exists, keeps it across exec, and so does
# everything it forks. mark_holders finds every process with the marker open
# (by device and inode, so a deleted or renamed marker still counts), and
# lock_free proves that none is left, including any this shell cannot see
# (another user's or a non-dumpable process). mark_stop stops them, with their
# process groups (until no member is left, even once a group's leader has
# exited), and succeeds only on that proof and when every one of them stopped.

# PATH is absent, or a regular file that is not a link.
plain_target() {
    [ ! -L "$1" ] && { [ ! -e "$1" ] || [ -f "$1" ]; }
}
# DIR is absent, or a directory that is not a link.
plain_dir() {
    [ ! -L "$1" ] && { [ ! -e "$1" ] || [ -d "$1" ]; }
}

# python3 -I -c "$GRM_PY" OP FILE [ARG]: read FILE without following a link or
# blocking on a FIFO; it must be a regular file. On failure print "ERR reason"
# and exit 1. OP is
#   sha256        print its SHA-256
#   line MAX      print its only line: at most MAX bytes of printable ASCII
#                 ending in its single newline, with no other byte after it
#   fields FORM   print "key value" lines of its key=value fields (common.sh
#                 parse_fields; FORM line or lines), at most 64 KiB
GRM_PY='
import hashlib, os, re, stat, sys

def fail(msg):
    print("ERR " + msg)
    sys.exit(1)

def open_plain(path):
    try:
        st = os.lstat(path)
    except FileNotFoundError:
        fail("is missing")
    except OSError as error:
        fail(f"cannot be read ({error.strerror})")
    if stat.S_ISLNK(st.st_mode):
        fail("is a link")
    if not stat.S_ISREG(st.st_mode):
        fail("is not a regular file")
    try:
        fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_NOCTTY)
    except OSError as error:
        fail(f"cannot be read ({error.strerror})")
    if not stat.S_ISREG(os.fstat(fd).st_mode):
        fail("is not a regular file")
    return fd

def read_capped(fd, limit, label):
    data = b""
    while len(data) <= limit:
        chunk = os.read(fd, limit + 1 - len(data))
        if not chunk:
            break
        data += chunk
    if len(data) > limit:
        fail(f"is larger than {label}")
    return data

def text(data):
    if not data:
        fail("is empty")
    for offset, byte in enumerate(data):
        if byte != 10 and not 32 <= byte <= 126:
            fail(f"holds byte 0x{byte:02x} at offset {offset}; only printable ASCII and newlines are allowed")
    if not data.endswith(b"\n"):
        fail("does not end with a newline")
    return data[:-1].decode("ascii")

op, path = sys.argv[1:3]
fd = open_plain(path)
if op == "sha256":
    digest = hashlib.sha256()
    while True:
        chunk = os.read(fd, 1 << 20)
        if not chunk:
            break
        digest.update(chunk)
    print(digest.hexdigest())
elif op == "line":
    limit = int(sys.argv[3])
    body = text(read_capped(fd, limit, f"{limit} bytes"))
    if "\n" in body:
        fail(f"has {body.count(chr(10)) + 1} lines, expected one")
    print(body)
elif op == "fields":
    body = text(read_capped(fd, 65536, "64 KiB"))
    if sys.argv[3] == "line":
        if "\n" in body:
            fail(f"has {body.count(chr(10)) + 1} lines, expected one")
        tokens = body.split(" ")
    else:
        tokens = body.split("\n")
    seen = set()
    out = []
    for token in tokens:
        if not token:
            fail("has an empty field (a blank line, or a doubled, leading or trailing space)")
        match = re.fullmatch(r"([a-z][a-z0-9_]*)=([A-Za-z0-9._:/+-]+)", token)
        if not match:
            fail(f"has a malformed field {token[:80]!r}")
        if match[1] in seen:
            fail(f"repeats field {match[1]}")
        seen.add(match[1])
        out.append(f"{match[1]} {match[2]}")
    print("\n".join(out))
else:
    fail(f"unknown operation {op}")
'

# safe_line FILE [MAX]: print FILE's only line (GRM_PY line; MAX default 4096).
safe_line() {
    local out
    out=$(python3 -I -c "$GRM_PY" line "$1" "${2:-4096}") || return 1
    printf '%s\n' "$out"
}
# safe_sha256 FILE: print the SHA-256 of a regular FILE.
safe_sha256() {
    local out
    out=$(python3 -I -c "$GRM_PY" sha256 "$1") || return 1
    printf '%s\n' "$out"
}

# open_fd VAR FILE: open FILE read-write (creating it when absent) on a new
# descriptor and store its number in VAR. Fails with nothing left open unless
# FILE was absent or a regular file that is not a link, and the descriptor is
# a regular file at that path: a link or FIFO swapped in after the check is
# caught on the descriptor, and a read-write open never blocks on a FIFO.
open_fd() {
    local __fd __want __self=$BASHPID
    plain_target "$2" || return 1
    __want=$(cd -- "$(dirname -- "$2")" 2>/dev/null && pwd -P) || return 1
    __want=$__want/${2##*/}
    { exec {__fd}<>"$2"; } 2>/dev/null || return 1
    if [ ! -L "$2" ] && [ -f "/proc/$__self/fd/$__fd" ] && [ "$(readlink "/proc/$__self/fd/$__fd")" = "$__want" ]; then
        printf -v "$1" '%s' "$__fd"
        return 0
    fi
    exec {__fd}>&-
    return 1
}
# open_append VAR FILE: open_fd, then the same file again for appending.
open_append() {
    local __rw __ap __self=$BASHPID
    open_fd __rw "$2" || return 1
    { exec {__ap}>>"/proc/$__self/fd/$__rw"; } 2>/dev/null || { exec {__rw}>&-; return 1; }
    exec {__rw}>&-
    printf -v "$1" '%s' "$__ap"
}
# close_fd VAR: close the descriptor whose number VAR holds, if any, and clear VAR.
close_fd() {
    local __n=${!1:-}
    [ -z "$__n" ] || exec {__n}>&-
    printf -v "$1" '%s' ''
}
# append_line FILE TEXT: best effort; never fails, never blocks.
append_line() {
    local fd=''
    open_append fd "$1" || return 0
    { printf '%s\n' "$2" >&"$fd"; } 2>/dev/null || true
    close_fd fd
}

# write_file FILE CONTENT: atomically replace FILE with CONTENT. FILE must be
# absent or a regular file, never a link, directory, FIFO or other special
# file, in a directory that is not a link. CONTENT goes to a new temporary file
# beside FILE, is checked there byte for byte, renamed over FILE with mv -T
# (which fails on a directory instead of moving into it), and checked byte for
# byte again in FILE through GRM_PY. Any failure returns 1; before the rename
# FILE is left as it was.
write_file() {
    local file=$1 content=$2 dir want tmp
    dir=${file%/*}
    [ -d "$dir" ] && [ ! -L "$dir" ] && plain_target "$file" || return 1
    want=$(printf '%s' "$content" | sha256sum | cut -d' ' -f1)
    tmp=$(mktemp -- "$dir/.${file##*/}.XXXXXX" 2>/dev/null) || return 1
    if printf '%s' "$content" 2>/dev/null >"$tmp" && chmod 644 "$tmp" &&
        checkpoint "write:${file##*/}" &&
        [ "$(safe_sha256 "$tmp")" = "$want" ] &&
        mv -fT -- "$tmp" "$file" 2>/dev/null &&
        checkpoint "renamed:${file##*/}" &&
        [ "$(safe_sha256 "$file")" = "$want" ]; then
        return 0
    fi
    rm -f -- "$tmp"
    return 1
}

# Test hooks (tests/apply.sh). GRM_TEST_SIGNAL_AT=NAME makes the script send
# itself GRM_TEST_SIGNAL (TERM by default; KILL and STOP are used too) at
# checkpoint NAME. GRM_TEST_CHILD_AT=NAME makes a forked child stop itself
# (SIGSTOP) at child_checkpoint NAME, before it execs. Both do nothing unless
# set.
checkpoint() {
    [ "${GRM_TEST_SIGNAL_AT:-}" = "$1" ] || return 0
    kill -s "${GRM_TEST_SIGNAL:-TERM}" "$$"
}
child_checkpoint() {
    [ "${GRM_TEST_CHILD_AT:-}" = "$1" ] || return 0
    kill -s STOP "$BASHPID"
}

# ---- Process identity --------------------------------------------------------

# proc_stat PID: set PS_STATE, PS_PPID, PS_PGRP, PS_SID and PS_START (clock
# ticks after boot, field 22 of /proc/PID/stat) and succeed, unless PID is
# gone or a zombie.
proc_stat() {
    local line
    local -a f
    { read -r line <"/proc/$1/stat"; } 2>/dev/null || return 1
    read -r -a f <<<"${line##*) }"
    [ "${f[0]:-Z}" != Z ] && [ "${f[0]}" != X ] || return 1
    PS_STATE=${f[0]} PS_PPID=${f[1]} PS_PGRP=${f[2]} PS_SID=${f[3]} PS_START=${f[19]}
}
# The start time of PID; fails when PID is gone or a zombie.
proc_start() {
    proc_stat "$1" || return 1
    printf '%s\n' "$PS_START"
}
# alive PID [START]: PID runs, and is the process that started at START.
alive() {
    proc_stat "$1" || return 1
    [ -z "${2:-}" ] || [ "$PS_START" = "$2" ]
}
# The parent PID of PID, or nothing.
ppid_of() {
    if proc_stat "$1"; then printf '%s\n' "$PS_PPID"; fi
}
# file_id FILE: "device:inode" of a regular FILE that is not a link.
file_id() {
    [ -f "$1" ] && [ ! -L "$1" ] || return 1
    stat -c %d:%i -- "$1"
}
# fd_id FD: "device:inode" of this shell's open descriptor FD.
fd_id() {
    stat -L -c %d:%i -- "/proc/$BASHPID/fd/$1"
}
# fds_on PID ID: PID's descriptors open on the file with device:inode ID.
fds_on() {
    local link out=''
    for link in "/proc/$1/fd/"*; do
        [ "$(stat -L -c %d:%i -- "$link" 2>/dev/null)" = "$2" ] && out+=" ${link##*/}"
    done
    printf '%s\n' "${out# }"
}
holds() { [ -n "$(fds_on "$1" "$2")" ]; }
# fd_holds_flock PID FD: PID's descriptor FD is the open file description
# holding an exclusive flock (the kernel lists it in /proc/PID/fdinfo/FD).
fd_holds_flock() {
    local line
    while IFS= read -r line; do
        case $line in lock:*FLOCK*WRITE*) return 0 ;; esac
    done 2>/dev/null <"/proc/$1/fdinfo/$2"
    return 1
}
# mark_holders ID: the PIDs, other than this shell, with a descriptor open on
# the file with device:inode ID, whatever its path is now. A process whose
# descriptors this user cannot read is missed; lock_free still sees it.
mark_holders() {
    local links link pid out=''
    links=$(find -L /proc/[0-9]*/fd -mindepth 1 -maxdepth 1 -inum "${1#*:}" 2>/dev/null) || true
    for link in $links; do
        pid=${link#/proc/}
        pid=${pid%%/*}
        [ "$pid" != "$$" ] && [ "$pid" != "$BASHPID" ] || continue
        [ "$(stat -L -c %d:%i -- "$link" 2>/dev/null)" = "$1" ] || continue
        case " $out " in *" $pid "*) ;; *) out+=" $pid" ;; esac
    done
    printf '%s\n' "${out# }"
}
# lock_free FILE: FILE is a regular file and no process holds a flock on it.
# Returns 1 while one does, 2 when FILE is missing or not a plain file.
lock_free() {
    local fd='' rc=0
    [ -f "$1" ] && [ ! -L "$1" ] || return 2
    open_fd fd "$1" || return 2
    flock -n "$fd" || rc=1
    close_fd fd
    return "$rc"
}

# ---- Stopping ------------------------------------------------------------------

# tree_signal SIGNAL PID GROUP: signal PID, and process group GROUP if given.
tree_signal() {
    kill -s "$1" "$2" 2>/dev/null || true
    if [ -n "$3" ]; then kill -s "$1" -- "-$3" 2>/dev/null || true; fi
}
# group_alive GROUP: a live (not zombie) process is in process group GROUP.
group_alive() {
    local p
    for p in $(pgrep -g "$1" 2>/dev/null || true); do
        proc_stat "$p" && return 0
    done
    return 1
}
# group_left GROUP START: process group GROUP, which the process that started at
# START leads or led, still has a live member. A group's ID is not reused while
# it has a member, so the members of a group whose leader has exited are still
# in that group; a live process with that PID and another start time shows that
# the group ended and the number was reused, so it counts as gone.
group_left() {
    if proc_stat "$1" && [ "$PS_START" != "$2" ]; then return 1; fi
    group_alive "$1"
}
# tree_send SIGNAL PID START GROUP: signal PID while it is the process that
# started at START, and process group GROUP, if given, while it has a member.
tree_send() {
    if alive "$2" "$3"; then kill -s "$1" "$2" 2>/dev/null || true; fi
    if [ -n "$4" ] && group_left "$4" "$3"; then kill -s "$1" -- "-$4" 2>/dev/null || true; fi
}
# stop_tree PID START [group]: stop the process that started at START, and its
# process group when it leads one (a process started through setsid), unless
# that is this shell's own group: SIGTERM (with SIGCONT, so a stopped process
# acts on it), up to 10 s, then SIGKILL, up to 10 s. The group gets each signal
# while it has a member, even after its leader has exited, so a descendant that
# ignores SIGTERM is still killed. With "group", the group that PID led is
# stopped even when PID has already exited. A PID now held by another process
# counts as stopped and gets no signal. Succeeds once all of them are gone.
stop_tree() {
    local pid=$1 start=$2 group='' sig _
    if alive "$pid" "$start"; then
        if [ "$PS_PGRP" = "$pid" ]; then group=$pid; fi
    elif [ "${3:-}" = group ] && group_left "$pid" "$start"; then
        group=$pid
    else
        return 0
    fi
    if [ -n "$group" ] && { ! proc_stat "$$" || [ "$PS_PGRP" = "$group" ]; }; then group=''; fi
    for sig in TERM KILL; do
        tree_send "$sig" "$pid" "$start" "$group"
        [ "$sig" = KILL ] || tree_send CONT "$pid" "$start" "$group"
        for _ in $(seq 1 40); do
            if ! alive "$pid" "$start" && { [ -z "$group" ] || ! group_left "$group" "$start"; }; then return 0; fi
            sleep 0.25
        done
    done
    return 1
}
# SIGTERM PID and wait up to 10 s for it to exit. Given START, a PID now held by
# another process counts as stopped and is not signalled.
stop_pid() {
    local _
    alive "$1" "${2:-}" || return 0
    kill -TERM "$1" 2>/dev/null || return 0
    for _ in $(seq 1 40); do
        alive "$1" "${2:-}" || return 0
        sleep 0.25
    done
    return 1
}
# SIGKILL PID [START]; succeeds only once the process is gone.
kill_pid() {
    local _
    alive "$1" "${2:-}" || return 0
    kill -KILL "$1" 2>/dev/null || true
    for _ in $(seq 1 20); do
        alive "$1" "${2:-}" || return 0
        sleep 0.25
    done
    return 1
}
# stop_pid, then kill_pid.
ensure_stopped() {
    stop_pid "$@" || kill_pid "$@"
}

# ---- Spawn markers -------------------------------------------------------------

# mark_lock FILE: open FILE (created if absent) on MARK_FD, set MARK_ID and take
# its flock. Returns 0 locked, 3 when another process holds the lock (MARK_FD
# stays open), 1 when FILE cannot be opened as a plain file (nothing open).
mark_lock() {
    MARK_FD='' MARK_ID=''
    open_fd MARK_FD "$1" || return 1
    MARK_ID=$(fd_id "$MARK_FD") || { close_fd MARK_FD; return 1; }
    flock -n "$MARK_FD" || return 3
}
# mark_create DIR NAME: make a new marker DIR/NAME.XXXXXX and mark_lock it; set
# MARK_PATH.
mark_create() {
    local path
    MARK_PATH='' MARK_FD='' MARK_ID=''
    plain_dir "$1" || return 1
    [ -d "$1" ] || mkdir -m 755 -- "$1" 2>/dev/null || return 1
    path=$(mktemp -- "$1/$2.XXXXXX" 2>/dev/null) || return 1
    if mark_lock "$path"; then
        MARK_PATH=$path
        return 0
    fi
    close_fd MARK_FD
    rm -f -- "$path"
    return 1
}
# mark_stop ID FILE [PID START]: stop every process holding the marker ID
# (and the registered PID, START, with the process group it leads or led, so a
# descendant that closed the marker is found even after PID exited), with their
# process groups, then require the proof that none is left: every stop
# succeeded, no holder this shell can see, the registered process and its group
# gone, and FILE, still the marker ID, with its lock free. Returns 1 with the
# reason in MARK_ERR otherwise. This shell's own descriptor on the marker
# (MARK_FD) must be closed first, or unlocked.
mark_stop() {
    local id=$1 file=$2 cpid=${3:-} cstart=${4:-} p s round pids
    MARK_ERR=''
    for round in 1 2 3; do
        pids=$(mark_holders "$id")
        if registered_left "$cpid" "$cstart"; then pids="$cpid $pids"; fi
        [ -n "${pids// /}" ] || break
        for p in $pids; do
            if [ "$p" = "$cpid" ]; then
                stop_tree "$p" "$cstart" group || MARK_ERR="pid $p or its process group would not stop"
            else
                s=$(proc_start "$p") || continue
                stop_tree "$p" "$s" || MARK_ERR="pid $p would not stop"
            fi
        done
    done
    # A process that would not stop may have closed the marker, so the proof
    # below would not see it.
    [ -z "$MARK_ERR" ] || return 1
    pids=$(mark_holders "$id")
    [ -z "$pids" ] || { MARK_ERR="pids $pids still hold the marker"; return 1; }
    if registered_left "$cpid" "$cstart"; then
        MARK_ERR="the registered pid $cpid or its process group would not stop"
        return 1
    fi
    [ "$(file_id "$file" 2>/dev/null)" = "$id" ] ||
        { MARK_ERR="${file##*/} is missing or replaced, so its lock cannot show that every process is gone"; return 1; }
    lock_free "$file" || { MARK_ERR="${file##*/} is still locked by a process this user cannot see"; return 1; }
}
# registered_left PID START: PID (if given) is still the process that started
# at START, or the process group it led has a member.
registered_left() {
    [ -n "$1" ] && [ -n "$2" ] && { alive "$1" "$2" || group_left "$1" "$2"; }
}

# cmdline_arg PID FLAG: the argument after FLAG in PID's command line.
cmdline_arg() {
    local -a argv
    local i
    { mapfile -d '' -t argv <"/proc/$1/cmdline"; } 2>/dev/null || return 0
    for ((i = 0; i + 1 < ${#argv[@]}; i++)); do
        [ "${argv[i]}" != "$2" ] || { printf '%s\n' "${argv[i + 1]}"; return 0; }
    done
}
# env_value PID NAME: the value of NAME in PID's environment, if readable.
env_value() {
    local -a env
    local entry
    { mapfile -d '' -t env <"/proc/$1/environ"; } 2>/dev/null || return 0
    for entry in "${env[@]}"; do
        [ "${entry%%=*}" != "$2" ] || { printf '%s\n' "${entry#*=}"; return 0; }
    done
}
# proc_comm PID: its command name.
proc_comm() {
    local comm=''
    { read -r comm <"/proc/$1/comm"; } 2>/dev/null || true
    printf '%s\n' "$comm"
}
