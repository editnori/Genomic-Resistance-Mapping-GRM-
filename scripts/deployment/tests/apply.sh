#!/usr/bin/env bash
# Run the deployment scripts with --apply against a throwaway root in a
# temporary directory. The releases' grm is a stub compiled from C that serves
# the app and report on 127.0.0.1 like grm serve; cloudflared, cargo and npm are
# shell stubs. Nothing here contacts the network, Cloudflare, a credential or a
# real deployment, and every process it starts is killed on exit. It needs cc,
# picks free loopback ports itself, and takes about four to six minutes. It
# interrupts starts through the scripts' GRM_TEST_SIGNAL_AT hook (common.sh).
# usage: bash scripts/deployment/tests/apply.sh
set -euo pipefail
# A non-interactive shell cannot trap a signal that was ignored when it
# started, so run as a background job (INT ignored) or under nohup (HUP
# ignored), the INT and HUP cases would test nothing. Start again with the
# default dispositions.
read -r _ sigign < <(grep '^SigIgn:' /proc/$$/status)
if ((0x$sigign & 0x3)) && [ -z "${GRM_TEST_SIGNALS_RESET:-}" ]; then
    exec env GRM_TEST_SIGNALS_RESET=1 python3 -I -c '
import os, signal, sys
for s in (signal.SIGINT, signal.SIGHUP, signal.SIGQUIT, signal.SIGTERM):
    signal.signal(s, signal.SIG_DFL)
os.execvp("bash", ["bash"] + sys.argv[1:])
' "${BASH_SOURCE[0]}" "$@"
fi
S=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
. "$S/tests/tmproot.sh"
T=$(make_test_root) || exit 2
export GRM_TEST_ROOT=$T
# Kill everything started under this test: each process inherits GRM_TEST_ROOT.
# A non-dumpable process hides its environment, so those are listed in EXTRA_PIDS.
EXTRA_PIDS=''
reap() {
    local p env
    for p in $EXTRA_PIDS; do kill -KILL "$p" 2>/dev/null || true; done
    for p in /proc/[0-9]*; do
        p=${p#/proc/}
        [ "$p" != $$ ] || continue
        env=$(tr '\0' '\n' 2>/dev/null <"/proc/$p/environ") || continue
        if grep -qxF "GRM_TEST_ROOT=$T" <<<"$env"; then kill -KILL "$p" 2>/dev/null || true; fi
    done
}
# Every assertion must run: a command outside expect that fails under set -e
# would otherwise end the suite without a FAIL line, so the EXIT trap reports
# it, and the total is checked at the end.
finished=0
on_exit() {
    local rc=$?
    if ((!finished)); then
        echo "FAIL apply suite aborted (exit $rc) after $((${pass:-0} + ${fail:-0})) assertions; a command outside expect failed"
        rc=1
    fi
    reap
    chmod -R u+w "$T"
    rm -rf "$T"
    exit "$rc"
}
trap on_exit EXIT
command -v cc >/dev/null || { echo "apply tests need cc to build the stub grm" >&2; exit 2; }
D=$T/root
# LP is the live port, SP the stage port, and XP a port where nothing listens
# (the proxy that health checks must not use).
read -r LP SP XP < <(python3 -I -c '
import socket
s = [socket.socket() for _ in range(3)]
for x in s: x.bind(("127.0.0.1", 0))
print(*[x.getsockname()[1] for x in s])')
export GRM_DEPLOY_ROOT=$D GRM_PORT=$LP GRM_SKIP_BROWSER=1
unset GRM_PLAYWRIGHT_CORE GRM_SOURCE_REPO GRM_TUNNEL_TOKEN_FILE GRM_EXPECT_SCHEMA GRM_EXPECT_REVISION GRM_EXPECT_DATASETS
pass=0
fail=0
# Every guard against a hang is "timeout -k 5 N": a script that traps TERM
# and is blocked in open() (a FIFO, say) never runs its trap, so TERM alone
# would leave it, and this suite, waiting forever.

expect() {  # expect CODE PATTERN NAME -- command...: exit code and output match
    local code=$1 pattern=$2 name=$3 got out
    shift 4
    set +e
    out=$("$@" 2>&1)
    got=$?
    set -e
    if [ "$got" = "$code" ] && grep -Eq -- "$pattern" <<<"$out"; then
        pass=$((pass + 1))
        echo "ok   $name"
    else
        fail=$((fail + 1))
        printf 'FAIL %s (exit %s, wanted %s /%s/)\n%s\n' "$name" "$got" "$code" "$pattern" "$out" | sed '2,$s/^/     /'
    fi
}
check() { local name=$1; shift 2; expect 0 '' "$name" -- "$@"; }
wait_for() {  # wait_for SECONDS command...: poll until it succeeds
    local i
    for ((i = 0; i < $1 * 10; i++)); do
        "${@:2}" 2>/dev/null && return 0
        sleep 0.1
    done
    return 1
}
listener() { fuser -n tcp "$1" 2>/dev/null | xargs || true; }
# "LABEL PUB" of the process listening on the live port.
live() {
    local pid data
    pid=$(listener "$LP")
    [ -n "$pid" ] || { echo none; return; }
    # sed reads all of its input (no q), so tr is never cut off by SIGPIPE.
    data=$(tr '\0' '\n' <"/proc/$pid/cmdline" | sed -n '/^--data$/{n;p}')
    printf '%s %s\n' "$(basename "$(dirname "$(readlink "/proc/$pid/exe")")")" "$(basename "$data")"
}
# PIDs running release LABEL's grm (one find, so a short-lived candidate is seen).
running() { find /proc -mindepth 2 -maxdepth 2 -name exe -lname "$D/release/$1/grm" 2>/dev/null | cut -d/ -f3 || true; }
has_running() { [ -n "$(running "$1")" ]; }
port_free() { [ -z "$(listener "$1")" ]; }
has_child() { pgrep -P "$1" -x sleep >/dev/null; }
field() { sed -n "s/^$2=//p" "$D/state/$1"; }
mark() { LOG0=$(wc -l <"$D/deploy.log"); }
logged() {  # the lines since mark contain PATTERN (read whole, then searched)
    local text
    text=$(sed -n "$((LOG0 + 1)),\$p" "$D/deploy.log")
    grep -q -- "$1" <<<"$text"
}
# Signal injection through the scripts' checkpoint hook (common.sh).
# at POINT SIGNAL CMD...: run CMD, which signals itself at POINT; sets rc, output in $T/at.out.
at() {
    local point=$1 sig=$2
    shift 2
    set +e
    env GRM_TEST_SIGNAL_AT="$point" GRM_TEST_SIGNAL="$sig" "$@" >"$T/at.out" 2>&1
    rc=$?
    set -e
}
stopped() { [ "$(sed 's/.*) //' "/proc/$1/stat" 2>/dev/null | cut -d' ' -f1)" = T ]; }
# paused POINT CMD...: run CMD in the background until it stops itself at POINT (bg=its PID).
paused() {
    local point=$1
    shift
    env GRM_TEST_SIGNAL_AT="$point" GRM_TEST_SIGNAL=STOP "$@" >"$T/bg.out" 2>&1 &
    bg=$!
    wait_for 30 stopped "$bg" || { echo "FAIL $* never reached $point"; fail=$((fail + 1)); }
}
resume() {  # continue the paused command; sets rc
    kill -CONT "$bg"
    set +e
    wait "$bg"
    rc=$?
    set -e
}
journal_gone() { [ ! -e "$D/serve-$1.starting" ]; }
pid_is() { grep -q "^pid=$1 " "$D/serve-$LP.pid"; }  # serve-LIVE.pid records PID
# After a failed or interrupted activation: rel1 with p237 serves under a recorded
# PID, state/active names it, and nothing of the rel2 start is left.
rel1_back() {
    test "$(live)" = "rel1 p237" && test -z "$(running rel2)" && test "$(field active label)" = rel1 &&
        journal_gone "$LP" && pid_is "$(listener "$LP")"
}
kill_live() {
    kill -KILL "$(listener "$LP")"
    wait_for 5 port_free "$LP" || true
}
listens_as() { [ "$(listener "$LP")" = "$1" ]; }
proc_alive() { [ -d "/proc/$1" ] && [ "$(sed 's/.*) //' "/proc/$1/stat" 2>/dev/null | cut -d' ' -f1)" != Z ]; }
gone() { ! proc_alive "$1"; }
# start_of PID: its start time (field 22 of /proc/PID/stat), or 0 once it is
# gone, so a process that exited early fails an assertion instead of
# aborting the suite on an unset field.
start_of() {
    local -a f
    read -r -a f <<<"$(sed 's/.*) //' "/proc/$1/stat" 2>/dev/null)" || true
    printf '%s\n' "${f[19]:-0}"
}
# The spawn marker named by serve-PORT.starting (default the live port).
jmark() { sed -n 's/.* mark=\([^ ]*\) .*/\1/p' "$D/serve-${1:-$LP}.starting" 2>/dev/null || true; }
locked() { ! flock -n "$1" true; }
# stopped_child PID: a child of PID has stopped itself (sets kid).
stopped_child() {
    local c
    for c in $(pgrep -P "$1" || true); do
        if stopped "$c"; then kid=$c; return 0; fi
    done
    return 1
}
# preexec CHILD-POINT CMD...: run CMD in the background until the child it
# forks stops itself at CHILD-POINT, before exec (bg=CMD's PID, kid=the child).
preexec() {
    local point=$1
    shift
    kid=''
    env GRM_TEST_CHILD_AT="$point" "$@" >"$T/bg.out" 2>&1 &
    bg=$!
    wait_for 30 stopped_child "$bg" || { echo "FAIL $* never paused its child at $point"; fail=$((fail + 1)); }
}
killbg() {  # SIGKILL the background command; sets rc
    kill -KILL "$bg" 2>/dev/null || true
    set +e
    wait "$bg" 2>/dev/null
    rc=$?
    set -e
}
# hidden_lock FILE: a non-dumpable process (its descriptors and environment
# unreadable to this user) takes FILE's flock (hidden=its PID). It waits for
# the lock, so the locked probe below cannot make it give up.
hidden_lock() {
    python3 -I -c '
import ctypes, fcntl, os, sys, time
ctypes.CDLL(None).prctl(4, 0)
fd = os.open(sys.argv[1], os.O_RDWR)
fcntl.flock(fd, fcntl.LOCK_EX)
time.sleep(600)
' "$1" >/dev/null 2>&1 &
    hidden=$!
    EXTRA_PIDS+=" $hidden"
    wait_for 10 locked "$1" || true
}

# Stub grm: grm --data DIR serve --app DIR --port N. GRM_STUB_FAIL, _SLOW,
# _BAD and _IGNORE_TERM, when a substring of its path, make it exit at once,
# wait 3 s before listening, serve a changed index.html, or ignore SIGTERM.
cat >"$T/stub.c" <<'C'
#include <fcntl.h>
#include <netinet/in.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <unistd.h>
static int mode(const char *name, const char *self) {
    const char *v = getenv(name);
    return v && *v && strstr(self, v);
}
int main(int argc, char **argv) {
    const char *data = 0, *app = 0;
    int port = 0, one = 1;
    for (int i = 1; i + 1 < argc; i++) {
        if (!strcmp(argv[i], "--data")) data = argv[++i];
        else if (!strcmp(argv[i], "--app")) app = argv[++i];
        else if (!strcmp(argv[i], "--port")) port = atoi(argv[++i]);
    }
    if (!data || !app || !port) return 2;
    if (mode("GRM_STUB_FAIL", argv[0])) return 1;
    if (mode("GRM_STUB_IGNORE_TERM", argv[0])) signal(SIGTERM, SIG_IGN);
    if (mode("GRM_STUB_SLOW", argv[0])) sleep(3);
    int bad = mode("GRM_STUB_BAD", argv[0]);
    signal(SIGPIPE, SIG_IGN);
    int s = socket(AF_INET, SOCK_STREAM, 0);
    setsockopt(s, SOL_SOCKET, SO_REUSEADDR, &one, sizeof one);
    struct sockaddr_in a = {0};
    a.sin_family = AF_INET;
    a.sin_port = htons(port);
    a.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if (bind(s, (struct sockaddr *)&a, sizeof a) || listen(s, 16)) { perror("bind"); return 1; }
    for (;;) {
        int c = accept(s, 0, 0);
        if (c < 0) continue;
        char req[4096], path[1024] = "", file[4096];
        ssize_t n = read(c, req, sizeof req - 1);
        req[n > 0 ? n : 0] = 0;
        sscanf(req, "GET %1023s", path);
        char *q = strchr(path, '?');
        if (q) *q = 0;
        if (!strcmp(path, "/")) snprintf(file, sizeof file, "%s/index.html", app);
        else if (!strncmp(path, "/data/", 6)) snprintf(file, sizeof file, "%s/report/%s", data, path + 6);
        else snprintf(file, sizeof file, "%s%s", app, path);
        struct stat st;
        int f = path[0] == '/' && !strstr(path, "..") ? open(file, O_RDONLY) : -1;
        if (f < 0 || fstat(f, &st) || !S_ISREG(st.st_mode)) {
            dprintf(c, "HTTP/1.0 404 Not Found\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
        } else {
            int extra = bad && !strcmp(path, "/");
            char buf[65536];
            ssize_t r;
            dprintf(c, "HTTP/1.0 200 OK\r\nContent-Length: %ld\r\nConnection: close\r\n\r\n", (long)st.st_size + extra);
            while ((r = read(f, buf, sizeof buf)) > 0 && write(c, buf, r) == r) continue;
            if (extra && write(c, "X", 1) != 1) perror("write");
        }
        if (f >= 0) close(f);
        close(c);
    }
}
C
cc -O2 -o "$T/grm-stub" "$T/stub.c"

# Fixture: a 237- and a 275-dataset publication and three releases.
make_report() {  # make_report DIR COUNT: report and DIR.sha256
    python3 -I - "$1" "$2" <<'PY'
import json, os, sys
root, n = sys.argv[1], int(sys.argv[2])
os.makedirs(os.path.join(root, "datasets"))
ids = [f"sp__d{i:03d}" for i in range(n)]
with open(os.path.join(root, "index.json"), "w") as f:
    json.dump({"schema": 2, "evaluation_revision": 7, "datasets": [{"id": i} for i in ids]}, f)
with open(os.path.join(root, "corpus.json"), "w") as f:
    json.dump({"corpus": n}, f)
for i in ids:
    with open(os.path.join(root, "datasets", i + ".json"), "w") as f:
        json.dump({"id": i}, f)
PY
    (cd "$1" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum) >"$1.sha256"
}
make_release() {
    local R=$D/release/$1
    mkdir -p "$R/dist/assets"
    cp "$T/grm-stub" "$R/grm"
    echo '<!doctype html><title>GRM</title>' >"$R/dist/index.html"
    echo "console.log('$1')" >"$R/dist/assets/index-test.js"
    echo "// browser check of $1" >"$R/browser-check.cjs"
    printf '%s\n' "$2" >"$R/COMMIT"
    (cd "$R" && find . -type f -print0 | LC_ALL=C sort -z | xargs -0 sha256sum >"$T/manifest" && mv "$T/manifest" MANIFEST.sha256)
    chmod -R a-w "$R"
}
mkdir -p "$D/release"
make_report "$T/r237" 237
make_report "$T/r275" 275
make_release rel1 1111111111111111111111111111111111111111
make_release rel2 2222222222222222222222222222222222222222
make_release rel3 3333333333333333333333333333333333333333
expect 0 'published p237' "publish 237 datasets" -- \
    bash "$S/publish.sh" "$T/r237" "$T/r237.sha256" p237 --apply
expect 0 'published p275' "publish 275 datasets with GRM_EXPECT_DATASETS=275" -- \
    env GRM_EXPECT_DATASETS=275 bash "$S/publish.sh" "$T/r275" "$T/r275.sha256" p275 --apply

# Stage both pairs; each approval records its own report and expectations.
expect 0 'stage pass' "stage rel1 with p237" -- bash "$S/stage.sh" rel1 p237 "$SP" --apply
expect 0 'stage pass' "stage rel2 with p275 (GRM_EXPECT_DATASETS=275)" -- \
    env GRM_EXPECT_DATASETS=275 bash "$S/stage.sh" rel2 p275 "$SP" --apply
check "approval binds the report manifest and 275 datasets" -- \
    grep -q "report_manifest=$(sha256sum <"$D/published/p275/SOURCE.sha256" | cut -d' ' -f1) schema=2 revision=7 datasets=275 receipts=stage/rel2--p275--" "$D/go/rel2--p275"
check "staging left no process and a free port" -- test -z "$(running rel1)$(running rel2)$(listener "$SP")"

# From here a contradictory global expectation must not matter.
export GRM_EXPECT_DATASETS=1

expect 0 'restore health pass: rel1 with p237' "restore starts the first pair" -- bash "$S/restore.sh" rel1 p237 --apply
check "state/active records 237 datasets" -- test "$(field active datasets)" = 237
live_pid=$(listener "$LP")
check "the service holds no deploy.lock descriptor" -- \
    test -z "$(find "/proc/$live_pid/fd" -lname "$D/deploy.lock" 2>/dev/null)"
check "deploy.lock is free while the service runs" -- flock -n "$D/deploy.lock" true

# Restore of a healthy service records missing state instead of exiting early.
expect 0 'healthy and recorded in state/active' "restore with recorded state changes nothing" -- bash "$S/restore.sh" --apply
rm "$D/state/active"
expect 0 'would record rel1 with p237' "restore dry run of a healthy unrecorded service" -- bash "$S/restore.sh" rel1 p237
check "the dry run wrote no state" -- test ! -e "$D/state/active"
expect 0 'recorded rel1 with p237' "restore --apply records a healthy unrecorded service" -- \
    bash "$S/restore.sh" rel1 p237 --apply
expect 0 'healthy and recorded' "a later argument-free restore finds it" -- bash "$S/restore.sh" --apply
check "the service was not restarted" -- test "$(listener "$LP")" = "$live_pid"

# The live pair must satisfy rollback's approval checks before it stops.
mv "$D/go/rel1--p237" "$T/rel1--p237"
expect 1 'no stage approval for rel1 with p237' "activate refuses when the live pair has no approval" -- \
    bash "$S/activate.sh" rel2 p275 --apply
echo "stage=0 browser=manual commit=2222222222222222222222222222222222222222" >"$D/GO-rel2"
echo "stage=0 commit=1111111111111111111111111111111111111111" >"$D/GO-rel3"
echo "stage=0 commit=1111111111111111111111111111111111111111 publication=other" >"$D/GO-rel1"
expect 1 'legacy approval' "activate refuses when the live pair has only a legacy approval" -- \
    bash "$S/activate.sh" rel2 p275 --apply
check "the live service kept running" -- test "$(listener "$LP")" = "$live_pid"
expect 1 'names another publication' "legacy migration refuses a GO file naming another publication" -- \
    env GRM_EXPECT_DATASETS=237 bash "$S/approve-legacy.sh" rel1 p237
echo "stage=0 schema1=0 browser=1042/1042 commit=1111111111111111111111111111111111111111" >"$D/GO-rel1"
expect 1 'datasets' "legacy migration checks the report against GRM_EXPECT_*" -- bash "$S/approve-legacy.sh" rel1 p237
expect 2 'only the live pair' "legacy migration refuses a pair that is not live" -- \
    env GRM_EXPECT_DATASETS=237 bash "$S/approve-legacy.sh" rel2 p237
expect 1 'names another commit' "legacy migration refuses another commit" -- \
    env GRM_EXPECT_DATASETS=237 bash "$S/approve-legacy.sh" rel3 p237
expect 0 'would write' "legacy migration dry run" -- env GRM_EXPECT_DATASETS=237 bash "$S/approve-legacy.sh" rel1 p237
check "the migration dry run wrote nothing" -- test ! -e "$D/go/rel1--p237"
expect 0 'approved rel1 with p237 from' "legacy migration of the live pair" -- \
    env GRM_EXPECT_DATASETS=237 bash "$S/approve-legacy.sh" rel1 p237 --apply
check "the migrated approval binds both manifests and the legacy file" -- \
    grep -Eq "^stage=0 browser=legacy .*release_manifest=[0-9a-f]{64} report_manifest=[0-9a-f]{64} schema=2 revision=7 datasets=237 legacy_sha256=$(sha256sum <"$D/GO-rel1" | cut -d' ' -f1)$" "$D/go/rel1--p237"
expect 0 'dry-run: would record' "activation accepts the migrated approval" -- bash "$S/activate.sh" rel2 p275
mv -f "$T/rel1--p237" "$D/go/rel1--p237"
rm "$D/GO-rel1" "$D/GO-rel2" "$D/GO-rel3"

# Failures before the old service stops change nothing: every path a start
# writes is checked, and state/previous and state/active are written and read
# back, before the stop.
live_pid=$(listener "$LP")
refused() {  # refused NAME PATTERN: activation exits 2 and the same rel1 pid keeps serving
    expect 2 "$2" "activate refuses $1; nothing stopped" -- timeout -k 5 60 bash "$S/activate.sh" rel2 p275 --apply
    check "after $1, the same rel1 pid serves and no rel2 process or start journal exists" -- \
        test "$(listener "$LP")" = "$live_pid" -a -z "$(running rel2)" -a ! -f "$D/serve-$LP.starting"
}
for name in previous active; do
    f=$D/state/$name
    if [ -e "$f" ]; then mv "$f" "$T/keep-$name"; fi
    mkdir "$f"
    refused "state/$name as a directory" 'cannot record state/previous and state/active'
    check "nothing was moved into the state/$name directory" -- test -z "$(ls -A "$f")"
    rmdir "$f"
    echo outside >"$T/outside"
    ln -s "$T/outside" "$f"
    refused "state/$name as a link" 'cannot record state/previous and state/active'
    check "the target of the state/$name link is unchanged" -- test "$(cat "$T/outside")" = outside
    rm "$f"
    ln -s "$T/nowhere" "$f"
    refused "state/$name as a dangling link" 'cannot record state/previous and state/active'
    check "the dangling state/$name link's target was not created" -- test ! -e "$T/nowhere"
    rm "$f"
    mkfifo "$f"
    refused "state/$name as a FIFO" 'cannot record state/previous and state/active'
    rm "$f"
    if [ -e "$T/keep-$name" ]; then mv "$T/keep-$name" "$f"; fi
done
if [ "$(id -u)" != 0 ]; then
    chmod 555 "$D/state"
    refused "an unwritable state/ directory" 'cannot record state/previous and state/active'
    chmod 755 "$D/state"
else
    echo "skip an unwritable state/ directory (root ignores the mode)"
fi
check "no temporary state file was left" -- test -z "$(find "$D/state" -name '.*')"

# A partial write is caught in the temporary file, before the rename.
cp "$D/state/previous" "$T/previous.before"
paused write:previous bash "$S/activate.sh" rel2 p275 --apply
tmp=$(find "$D/state" -name '.previous.*')
truncate -s 20 "$tmp"
resume
check "a truncated temporary state/previous fails activation before the stop (exit 2)" -- test "$rc" = 2
check "state/previous is unchanged and its temporary file removed" -- \
    test -z "$(find "$D/state" -name '.*')" -a "$(cat "$D/state/previous")" = "$(cat "$T/previous.before")"
check "the same rel1 pid serves" -- test "$(listener "$LP")" = "$live_pid"
# A state file that differs after its rename is caught by the read-back.
paused renamed:active bash "$S/activate.sh" rel2 p275 --apply
echo extra=1 >>"$D/state/active"
resume
check "a state/active changed after its rename fails activation before the stop (exit 2)" -- test "$rc" = 2
check "the same rel1 pid serves" -- test "$(listener "$LP")" = "$live_pid"
expect 0 'recorded rel1 with p237' "restore re-records the healthy pair over the bad state/active" -- \
    bash "$S/restore.sh" rel1 p237 --apply
# Killed between writing the temporary file and renaming it, the old record stays.
at write:active KILL bash "$S/activate.sh" rel2 p275 --apply
check "killed mid-write, state/active still names rel1 and parses" -- \
    test "$(field active label)" = rel1 -a -n "$(field active at)"
expect 0 'healthy and recorded' "after the killed write, restore finds rel1 recorded" -- bash "$S/restore.sh" --apply
find "$D/state" -name '.active.*' -delete

# The other paths a start writes.
mkdir "$D/serve-$LP.starting"
refused "serve-PORT.starting as a directory" 'is not a regular file'
rmdir "$D/serve-$LP.starting"
for kind in directory link fifo; do
    mv "$D/serve-$LP.pid" "$T/pid.keep"
    case $kind in
        directory) mkdir "$D/serve-$LP.pid" ;;
        link) ln -s "$T/pid-elsewhere" "$D/serve-$LP.pid" ;;
        fifo) mkfifo "$D/serve-$LP.pid" ;;
    esac
    refused "serve-PORT.pid as a $kind" 'link, directory or special file; nothing was stopped'
    rm -rf "$D/serve-$LP.pid"
    mv "$T/pid.keep" "$D/serve-$LP.pid"
done
check "the pid-file link's target was not created" -- test ! -e "$T/pid-elsewhere"
mkfifo "$D/logs/serve-$LP.log.fifo"
mv "$D/logs/serve-$LP.log" "$T/serve.log"
mv "$D/logs/serve-$LP.log.fifo" "$D/logs/serve-$LP.log"
refused "logs/serve-PORT.log as a FIFO" 'link, directory or special file; nothing was stopped'
rm "$D/logs/serve-$LP.log"
mv "$T/serve.log" "$D/logs/serve-$LP.log"
mv "$D/deploy.lock" "$T/deploy.lock"
mkfifo "$D/deploy.lock"
refused "deploy.lock as a FIFO" 'deploy.lock is a link, directory or special file'
rm "$D/deploy.lock"
mv "$T/deploy.lock" "$D/deploy.lock"
mv "$D/deploy.log" "$T/deploy.log"
mkfifo "$D/deploy.log"
expect 0 'healthy and recorded' "a FIFO deploy.log is skipped, not opened" -- timeout -k 5 60 bash "$S/restore.sh" --apply
rm "$D/deploy.log"
mv "$T/deploy.log" "$D/deploy.log"
expect 1 'serving again' "a candidate failing its health check restores the old pair" -- \
    env GRM_STUB_BAD=/rel2/ bash "$S/activate.sh" rel2 p275 --apply
check "rel1 with p237 serves after the failed health check" -- test "$(live)" = "rel1 p237"

mark
GRM_STUB_SLOW=/rel2/ bash "$S/activate.sh" rel2 p275 --apply >"$T/act.out" 2>&1 &
apid=$!
check "activation reaches the candidate start" -- wait_for 15 logged "started rel2"
kill -TERM "$apid" || true
set +e
wait "$apid"
rc=$?
set -e
check "SIGTERM during activation exits 1 after restoring the old pair" -- test "$rc" = 1
check "rel1 with p237 serves after the interruption" -- test "$(live)" = "rel1 p237"
check "the interrupted candidate was stopped" -- test -z "$(running rel2)"
check "the interrupted activation reports the old pair serving again" -- grep -q 'serving again' "$T/act.out"

# A closed output pipe cannot cut recovery short: the reader of activation's
# output exits once the candidate has started, the candidate then fails its
# health check, and every later write fails instead of killing the script.
set +e
GRM_STUB_SLOW=/rel2/ GRM_STUB_BAD=/rel2/ bash -c 'bash "$1/activate.sh" rel2 p275 --apply 2>&1; echo "$?" >"$2"' _ "$S" "$T/pipe.rc" |
    sed -n '/started rel2/q'
set -e
check "with its output pipe closed after the candidate started, the failed activation exits 1" -- \
    test "$(cat "$T/pipe.rc" 2>/dev/null)" = 1
check "with its output pipe closed, recovery stopped the candidate and rel1 serves again, recorded" -- rel1_back
# The same closed pipe after a switch that succeeds: the closing report fails to
# write, and say ignores that, so the activation still exits 0.
set +e
GRM_STUB_SLOW=/rel2/ bash -c 'bash "$1/activate.sh" rel2 p275 --apply 2>&1; echo "$?" >"$2"' _ "$S" "$T/pipe-ok.rc" |
    sed -n '/started rel2/q'
set -e
check "with its output pipe closed after the candidate started, a successful activation exits 0" -- \
    test "$(cat "$T/pipe-ok.rc" 2>/dev/null)" = 0
check "that activation switched to rel2 and recorded it" -- \
    eval 'test "$(live)" = "rel2 p275" -a "$(field active label)" = rel2 && journal_gone "$LP" && pid_is "$(listener "$LP")"'
expect 0 'live health pass: rel1 with p237' "rollback returns to rel1 after the switch with a closed pipe" -- bash "$S/rollback.sh" --apply

mark
GRM_STUB_FAIL=/rel2/ bash "$S/activate.sh" rel2 p275 --apply >"$T/act.out" 2>&1 &
apid=$!
check "activation reaches the candidate start" -- wait_for 15 logged "started rel2"
kill -KILL "$apid" || true
wait "$apid" 2>/dev/null || true
check "after SIGKILL nothing serves" -- test "$(live)" = none
check "state/active still names the stopped pair" -- test "$(field active label)" = rel1
expect 0 'restore health pass: rel1 with p237' "argument-free restore recovers after SIGKILL" -- \
    bash "$S/restore.sh" --apply

# Every step of a start can be interrupted. TERM (and INT, HUP) at any
# checkpoint, including between the fork and the assignment of the PID, stops
# the candidate, found by $! and its token, and restores rel1.
for point in old-stopped journaled spawned registered assigned healthy committed recorded; do
    at "$point" TERM bash "$S/activate.sh" rel2 p275 --apply
    check "activation TERM at $point exits 1 saying rel1 serves again" -- \
        grep -q 'rel1 with p237 is serving again' "$T/at.out"
    check "activation TERM at $point exit code is 1" -- test "$rc" = 1
    check "activation TERM at $point: rel1 serves, recorded; no rel2 process or journal is left" -- rel1_back
done
# Killed while recovery restarts rel1 (the candidate failed its health check):
# the restart's own journal lets restore.sh stop that orphan and start afresh.
at restarted KILL env GRM_STUB_BAD=/rel2/ bash "$S/activate.sh" rel2 p275 --apply
check "an activation killed during its rollback restart leaves a journal naming rel1 and registering its child" -- \
    test "$rc" = 137 -a -n "$(sed -n 's/.* label=\(rel1\) .* child=[0-9]* .*/\1/p' "$D/serve-$LP.starting")"
m=$(jmark)
check "the restart's marker is held by its child (or the child already exited)" -- test -e "$D/spawn/$m"
expect 0 'restore health pass: rel1 with p237' "restore stops the half-restarted rel1 and starts it afresh" -- \
    bash "$S/restore.sh" --apply
check "the half-restarted start's marker was released and removed" -- test ! -e "$D/spawn/$m"
check "after the half-restarted rel1, rel1 is back" -- rel1_back
at spawned INT bash "$S/activate.sh" rel2 p275 --apply
check "activation INT before the PID is assigned exits 1" -- test "$rc" = 1
check "after INT at spawned, rel1 is back and no rel2 is left" -- rel1_back
at assigned HUP bash "$S/activate.sh" rel2 p275 --apply
check "activation HUP after the PID is assigned exits 1" -- test "$rc" = 1
check "after HUP at assigned, rel1 is back and no rel2 is left" -- rel1_back

# KILL leaves the journal; the argument-free restore stops what the dead
# activation started (found by its token) and restarts state/active.
for point in old-stopped journaled spawned registered assigned healthy; do
    at "$point" KILL bash "$S/activate.sh" rel2 p275 --apply
    check "activation KILL at $point leaves state/active naming rel1" -- test "$rc" = 137 -a "$(field active label)" = rel1
    m=''
    case $point in
        old-stopped) ;;
        spawned)
            # Killed before it registered its child: the child never execs.
            m=$(jmark)
            check "after activation KILL at spawned the unregistered child exits without exec'ing rel2" -- \
                wait_for 10 eval '! locked "$D/spawn/$m"'
            check "no rel2 ran after activation KILL at spawned" -- test -z "$(running rel2)$(listener "$LP")"
            ;;
        healthy)
            m=$(jmark)
            check "after activation KILL at healthy an orphaned rel2 runs" -- wait_for 10 has_running rel2
            ;;
        *) m=$(jmark) ;;
    esac
    expect 0 'restore health pass: rel1 with p237' "after activation KILL at $point, restore --apply recovers rel1" -- \
        bash "$S/restore.sh" --apply
    check "after activation KILL at $point and restore, nothing of rel2 or its start is left" -- rel1_back
    if [ -n "$m" ]; then check "after activation KILL at $point, that start's marker was removed" -- test ! -e "$D/spawn/$m"; fi
done
at committed KILL bash "$S/activate.sh" rel2 p275 --apply
check "activation KILL after the state commit leaves state/active naming rel2" -- test "$(field active label)" = rel2
orphan=$(running rel2)
expect 0 'restore health pass: rel2 with p275' "restore restarts the committed rel2 from a clean start" -- \
    bash "$S/restore.sh" --apply
check "the orphan whose pid was never recorded was replaced" -- \
    test -n "$orphan" -a "$(listener "$LP")" != "$orphan" -a ! -e "/proc/$orphan/exe"
expect 0 'live health pass: rel1 with p237' "rollback returns to rel1" -- bash "$S/rollback.sh" --apply
at recorded KILL bash "$S/activate.sh" rel2 p275 --apply
committed_pid=$(listener "$LP")
expect 0 'had committed' "restore keeps a start killed after recording its pid" -- bash "$S/restore.sh" --apply
check "the committed rel2 pid was not restarted" -- test "$(listener "$LP")" = "$committed_pid" -a "$(live)" = "rel2 p275"
check "the committed start's journal is gone" -- journal_gone "$LP"
expect 0 'live health pass: rel1 with p237' "rollback returns to rel1" -- bash "$S/rollback.sh" --apply

# A child paused just before exec while its script is killed outright: it
# holds the start's marker from the fork on and is registered in the journal,
# so the resolving script stops it and removes the journal only once the
# marker's lock is free; resumed afterwards it does nothing. Resumed first, it
# sees its script gone and exits without exec'ing.
preexec_resolved() {  # preexec_resolved NAME RESOLVER-PATTERN RESOLVER... (after preexec)
    local name=$1 pattern=$2 m
    shift 2
    check "$name: the child paused before exec" -- \
        test -n "$kid" -a "$(readlink "/proc/$kid/exe")" = "$(readlink -f "$(command -v bash)")"
    wait_for 10 grep -q " child=$kid " "$D/serve-$port.starting" || true
    killbg
    m=$(jmark "$port")
    check "$name: the killed script's journal registers the paused child" -- grep -q " child=$kid " "$D/serve-$port.starting"
    check "$name: the paused child holds the marker's lock" -- locked "$D/spawn/$m"
    expect 0 "$pattern" "$name: the resolving script stops the paused child and proceeds" -- "$@"
    check "$name: the paused child is gone and its marker removed" -- \
        eval 'wait_for 10 gone "$kid" && test ! -e "$D/spawn/$m"'
    kill -CONT "$kid" 2>/dev/null || true
}
port=$LP
preexec child-preexec bash "$S/activate.sh" rel2 p275 --apply
preexec_resolved "activation killed with its child paused before exec" 'restore health pass: rel1 with p237' \
    bash "$S/restore.sh" --apply
check "after the activation's paused child was resumed, rel1 serves and no rel2 runs" -- rel1_back
preexec child-preexec bash "$S/activate.sh" rel2 p275 --apply
killbg
kill -CONT "$kid"
check "resumed after its activation was killed, the pre-exec child exits without exec'ing" -- wait_for 10 gone "$kid"
check "that child started no rel2 and nothing listens" -- test -z "$(running rel2)$(listener "$LP")"
check "the killed activation's journal is kept until it is resolved" -- test -e "$D/serve-$LP.starting"
expect 0 'restore health pass: rel1 with p237' "restore then resolves that journal and restores rel1" -- bash "$S/restore.sh" --apply
check "after the resumed child, rel1 is back" -- rel1_back
preexec child-forked bash "$S/activate.sh" rel2 p275 --apply
preexec_resolved "activation killed with its child paused at the fork" 'restore health pass: rel1 with p237' \
    bash "$S/restore.sh" --apply
check "after the child paused at the fork, rel1 is back" -- rel1_back
expect 0 'live health pass: rel2 with p275' "activate rel2 so that rollback has a target" -- bash "$S/activate.sh" rel2 p275 --apply
preexec child-preexec bash "$S/rollback.sh" --apply
preexec_resolved "rollback killed with its child paused before exec" 'restore health pass: rel2 with p275' \
    bash "$S/restore.sh" --apply
check "after the rollback's paused child, rel2 (the recorded pair) serves and no rel1 is left" -- \
    test "$(live)" = "rel2 p275" -a "$(field active label)" = rel2 -a -z "$(running rel1)" -a ! -e "$D/serve-$LP.starting"
expect 0 'live health pass: rel1 with p237' "activating rel1 then returns to it (the killed rollback recorded rel2 as previous)" -- bash "$S/activate.sh" rel1 p237 --apply
kill_live
preexec child-preexec bash "$S/restore.sh" --apply
preexec_resolved "restore killed with its child paused before exec" 'restore health pass: rel1 with p237' \
    bash "$S/restore.sh" --apply
check "after the restore's paused child, rel1 is back" -- rel1_back

# Unknown survival fails closed: a process this user cannot see holds the
# marker's lock, or the marker is gone, so restore keeps the journal, starts
# nothing and fails; once the holder is gone (or the journal is removed by
# hand) restore proceeds.
at journaled KILL bash "$S/activate.sh" rel2 p275 --apply
m=$(jmark)
hidden_lock "$D/spawn/$m"
check "a hidden process holds the interrupted start's marker" -- locked "$D/spawn/$m"
expect 1 'may remain: .*locked by a process this user cannot see; serve-.*starting is kept' \
    "restore refuses while a process it cannot see holds the marker's lock" -- bash "$S/restore.sh" --apply
check "the journal and marker are kept and nothing started" -- \
    test -e "$D/serve-$LP.starting" -a -e "$D/spawn/$m" -a -z "$(listener "$LP")$(running rel1)"
expect 2 'interrupted start of rel2' "activation still refuses that port" -- bash "$S/activate.sh" rel2 p275 --apply
kill -KILL "$hidden"
wait "$hidden" 2>/dev/null || true
expect 0 'restore health pass: rel1 with p237' "once the hidden holder is gone, restore resolves the journal" -- bash "$S/restore.sh" --apply
check "after the hidden holder, rel1 is back and the marker is removed" -- eval 'rel1_back && test ! -e "$D/spawn/$m"'
at journaled KILL bash "$S/activate.sh" rel2 p275 --apply
m=$(jmark)
rm "$D/spawn/$m"
expect 1 'is missing or replaced, so its lock cannot show' "restore refuses when the start's marker is gone" -- bash "$S/restore.sh" --apply
check "the journal of the start whose marker is gone is kept" -- test -e "$D/serve-$LP.starting" -a -z "$(listener "$LP")"
rm "$D/serve-$LP.starting"
expect 0 'restore health pass: rel1 with p237' "with the journal removed by hand, restore proceeds" -- bash "$S/restore.sh" --apply
check "after the missing marker, rel1 is back" -- rel1_back

# The state commit after the stop fails: the candidate stops and rel1 serves again.
paused healthy bash "$S/activate.sh" rel2 p275 --apply
mv "$D/state/active" "$T/active.keep"
mkdir "$D/state/active"
resume
check "a state/active commit failure after the stop exits 1" -- test "$rc" = 1
check "it reports that state/active could not be rewritten" -- grep -q 'could not rewrite state/active' "$T/bg.out"
check "after the commit failure rel1 serves and no rel2 is left" -- test "$(live)" = "rel1 p237" -a -z "$(running rel2)"
check "nothing was moved into the state/active directory" -- test -z "$(ls -A "$D/state/active")"
rmdir "$D/state/active"
mv "$T/active.keep" "$D/state/active"
expect 0 'recorded rel1 with p237|healthy and recorded' "restore records the serving rel1" -- bash "$S/restore.sh" rel1 p237 --apply

# serve-PORT.pid fails at its exact boundary, after the state commit: the
# committed pair keeps serving with a warning.
paused committed bash "$S/activate.sh" rel2 p275 --apply
rm "$D/serve-$LP.pid"
mkdir "$D/serve-$LP.pid"
resume
check "activation with a pid-file failure after the commit exits 0" -- test "$rc" = 0
check "it warns that serve-PORT.pid was not written" -- grep -q 'warning: could not write serve-.*pid' "$T/bg.out"
check "rel2 with p275 serves and state/active records it" -- test "$(live)" = "rel2 p275" -a "$(field active label)" = rel2
check "the journal of that start is gone" -- journal_gone "$LP"
rmdir "$D/serve-$LP.pid"
expect 0 'live health pass: rel1 with p237' "rollback returns to rel1" -- bash "$S/rollback.sh" --apply

# Restore from nothing: interrupted at any checkpoint it stops what it started,
# leaving the port free and state/active as it was.
kill_live
for point in journaled spawned registered assigned healthy committed recorded; do
    at "$point" TERM bash "$S/restore.sh" --apply
    check "restore TERM at $point exits 1" -- test "$rc" = 1
    check "restore TERM at $point leaves no rel1 process, listener or journal" -- \
        test -z "$(running rel1)$(listener "$LP")" -a ! -e "$D/serve-$LP.starting"
    check "restore TERM at $point leaves state/active naming rel1" -- test "$(field active label)" = rel1
done
at spawned KILL bash "$S/restore.sh" --apply
m=$(jmark)
check "after restore KILL at spawned, the unregistered child exits without exec'ing" -- wait_for 10 eval '! locked "$D/spawn/$m"'
check "after restore KILL at spawned nothing runs or listens" -- test -z "$(running rel1)$(listener "$LP")"
expect 0 'restore health pass: rel1 with p237' "after restore KILL at spawned, the next restore resolves the journal" -- \
    bash "$S/restore.sh" --apply
check "the journal and marker of restore KILL at spawned are gone" -- eval 'journal_gone "$LP" && test ! -e "$D/spawn/$m"'
kill_live
for point in healthy committed; do
    at "$point" KILL bash "$S/restore.sh" --apply
    check "after restore KILL at $point, its orphan runs" -- wait_for 10 has_running rel1
    orphan=$(running rel1)
    expect 0 'restore health pass: rel1 with p237' "after restore KILL at $point, the next restore stops the orphan and starts afresh" -- \
        bash "$S/restore.sh" --apply
    check "the orphan of restore KILL at $point is gone" -- test "$(listener "$LP")" != "$orphan" -a ! -e "/proc/$orphan/exe"
    kill_live
done
at recorded KILL bash "$S/restore.sh" --apply
committed_pid=$(listener "$LP")
expect 0 'already serving.*recorded' "a restore killed after recording its pid is kept by the next restore" -- \
    bash "$S/restore.sh" --apply
check "the recorded restore was not restarted and its journal is gone" -- \
    test "$(listener "$LP")" = "$committed_pid" -a ! -e "$D/serve-$LP.starting"

# serve-PORT.pid in restore: refused before anything starts when it is a
# directory, link or FIFO; replaced when it is a read-only file; and a
# failure at its exact boundary keeps the committed, healthy pair.
kill_live
mv "$D/serve-$LP.pid" "$T/pid.keep"
for kind in directory link fifo; do
    case $kind in
        directory) mkdir "$D/serve-$LP.pid" ;;
        link) ln -s "$T/pid-elsewhere" "$D/serve-$LP.pid" ;;
        fifo) mkfifo "$D/serve-$LP.pid" ;;
    esac
    expect 2 'nothing started' "restore refuses serve-PORT.pid as a $kind before starting" -- \
        timeout -k 5 60 bash "$S/restore.sh" --apply
    check "with serve-PORT.pid as a $kind nothing started and no journal is left" -- \
        test -z "$(running rel1)$(listener "$LP")" -a ! -e "$D/serve-$LP.starting"
    rm -rf "$D/serve-$LP.pid"
done
check "the pid-file link's target was not created" -- test ! -e "$T/pid-elsewhere"
mv "$T/pid.keep" "$D/serve-$LP.pid"
chmod 444 "$D/serve-$LP.pid"
expect 0 'restore health pass: rel1 with p237' "restore replaces a read-only serve-PORT.pid atomically" -- bash "$S/restore.sh" --apply
check "serve-PORT.pid records the restored pid" -- pid_is "$(listener "$LP")"
kill_live
paused committed bash "$S/restore.sh" --apply
rm "$D/serve-$LP.pid"
mkdir "$D/serve-$LP.pid"
resume
check "restore with a pid-file failure after the commit exits 0" -- test "$rc" = 0
check "it warns that serve-PORT.pid was not written" -- grep -q 'warning: could not write serve-.*pid' "$T/bg.out"
check "the restored rel1 serves and state/active records it" -- test "$(live)" = "rel1 p237" -a "$(field active label)" = rel1
check "the journal of that restore is gone" -- journal_gone "$LP"
rmdir "$D/serve-$LP.pid"
expect 0 'recorded rel1 with p237' "the next restore records the pid" -- bash "$S/restore.sh" --apply
check "serve-PORT.pid now records the live pid" -- pid_is "$(listener "$LP")"

# An impostor serving rel1's exact bytes takes the port while it is free; neither
# the candidate nor the restarted rel1 owns the listener, so both checks fail.
imp=$T/impostor
mkdir -p "$imp/assets"
cp "$D/release/rel1/dist/index.html" "$imp/"
cp "$D/release/rel1/dist/assets/index-test.js" "$imp/assets/"
ln -s "$D/published/p237/report" "$imp/data"
python3 -I - "$LP" "$imp" >/dev/null 2>&1 <<'PY' &
import functools, http.server, sys, time
port, root = int(sys.argv[1]), sys.argv[2]
class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args): pass
class Server(http.server.HTTPServer):
    allow_reuse_address = True
while True:
    try:
        server = Server(("127.0.0.1", port), functools.partial(Quiet, directory=root))
        break
    except OSError:
        time.sleep(0.005)
server.serve_forever()
PY
ipid=$!
expect 4 'nothing healthy' "a listener that is not the restarted pid fails the rollback" -- \
    env GRM_STUB_SLOW=/rel2/ bash "$S/activate.sh" rel2 p275 --apply
check "the impostor held the port" -- test "$(listener "$LP")" = "$ipid"
check "no rel1 or rel2 process claims success" -- test -z "$(running rel1)$(running rel2)"
kill "$ipid"
wait "$ipid" 2>/dev/null || true
wait_for 5 port_free "$LP" || true
expect 0 'restore health pass' "restore recovers once the port is free" -- bash "$S/restore.sh" --apply

# 237 -> 275 -> 237 -> 275, a restart, and back, with GRM_EXPECT_DATASETS=1 set.
expect 0 'live health pass: rel2 with p275' "activate 275 over 237" -- bash "$S/activate.sh" rel2 p275 --apply
check "state records 275 active and 237 previous" -- \
    test "$(field active datasets) $(field previous datasets)" = "275 237"
check "the new service holds no deploy.lock descriptor" -- \
    test -z "$(find "/proc/$(listener "$LP")/fd" -lname "$D/deploy.lock" 2>/dev/null)"
mkdir "$T/tmp"
expect 0 'smoke-pass .*datasets=275' "status checks the live pair against its approval" -- \
    env TMPDIR="$T/tmp" bash "$S/status.sh"
check "status and smoke.sh leave no temporary files" -- test -z "$(ls -A "$T/tmp")"
expect 0 'live health pass: rel1 with p237' "rollback 275 -> 237" -- bash "$S/rollback.sh" --apply
check "rel1 with p237 serves" -- test "$(live)" = "rel1 p237"
expect 0 'live health pass: rel2 with p275' "rollback again switches forward to 275" -- bash "$S/rollback.sh" --apply
kill -KILL "$(listener "$LP")"
wait_for 5 port_free "$LP" || true
expect 0 'restore health pass: rel2 with p275' "argument-free restore after a crash uses 275 from state" -- \
    bash "$S/restore.sh" --apply
expect 0 'live health pass: rel1 with p237' "rollback to 237 after the restore" -- bash "$S/rollback.sh" --apply

# Health requests ignore proxy settings: every proxy variable and a ~/.curlrc
# name a loopback port where nothing listens, which would fail every check,
# the rollback restart's included. The ~/.curlrc also sends connections for
# the live port to that port (connect-to), which --noproxy does not undo; only
# -q, which skips ~/.curlrc, does.
mkdir "$T/curlhome"
printf 'proxy = "http://127.0.0.1:%s"\nconnect-to = "127.0.0.1:%s:127.0.0.1:%s"\n' "$XP" "$LP" "$XP" >"$T/curlhome/.curlrc"
proxied() {
    env HOME="$T/curlhome" http_proxy="http://127.0.0.1:$XP" HTTP_PROXY="http://127.0.0.1:$XP" \
        all_proxy="http://127.0.0.1:$XP" ALL_PROXY="http://127.0.0.1:$XP" no_proxy= NO_PROXY= "$@"
}
check "a plain curl with those settings fails (the control)" -- \
    eval '! proxied curl -sf -o /dev/null --max-time 2 "http://127.0.0.1:$LP/data/index.json"'
check "a curl with --noproxy but without -q fails too: ~/.curlrc redirects it (the control for -q)" -- \
    eval '! proxied curl --noproxy "*" -sf -o /dev/null --max-time 2 "http://127.0.0.1:$LP/data/index.json"'
check "curl -q --noproxy '*' with those settings reaches the live port" -- \
    eval 'proxied curl -q --noproxy "*" -sf -o /dev/null --max-time 2 "http://127.0.0.1:$LP/data/index.json"'
expect 1 'serving again' "with proxy settings, a failed candidate's rollback restart is checked directly" -- \
    proxied env GRM_STUB_BAD=/rel2/ bash "$S/activate.sh" rel2 p275 --apply
expect 0 'live health pass: rel2 with p275' "with proxy settings, activation is checked directly" -- \
    proxied bash "$S/activate.sh" rel2 p275 --apply
expect 0 'smoke-pass .*datasets=275' "with proxy settings, status is checked directly" -- proxied bash "$S/status.sh"
expect 0 'live health pass: rel1 with p237' "with proxy settings, rollback is checked directly" -- \
    proxied bash "$S/rollback.sh" --apply

# A damaged current release does not block a rollback to an intact previous
# pair, but the damaged pair must be exactly what state/active records, is
# never a rollback target and is never restarted.
index=$D/release/rel2/dist/index.html
cp "$index" "$T/index.rel2"
damage() { chmod u+w "$D/release/rel2/dist" "$index"; echo tampered >>"$index"; }
repair() { cat "$T/index.rel2" >"$index"; }
impostor_up() {  # serve rel1's exact bytes from an unmanaged process on the live port (ipid)
    python3 -I - "$LP" "$imp" >/dev/null 2>&1 <<'PY' &
import functools, http.server, sys
port, root = int(sys.argv[1]), sys.argv[2]
class Quiet(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *args): pass
class Server(http.server.HTTPServer):
    allow_reuse_address = True
Server(("127.0.0.1", port), functools.partial(Quiet, directory=root)).serve_forever()
PY
    ipid=$!
    wait_for 10 listens_as "$ipid" || true
}
expect 0 'live health pass: rel2 with p275' "activate an intact rel2 before damaging it" -- bash "$S/activate.sh" rel2 p275 --apply
damage
damaged_pid=$(listener "$LP")
expect 1 'does not match' "activate refuses to replace a damaged live pair outside a rollback" -- \
    bash "$S/activate.sh" rel1 p237 --apply
cp "$D/state/active" "$T/active.rel2"
cp "$D/state/previous" "$D/state/active"
expect 2 'not the pair state/active records' "rollback refuses a damaged live pair that state/active does not name" -- \
    bash "$S/rollback.sh" --apply
check "the unidentified damaged pair was not stopped" -- test "$(listener "$LP")" = "$damaged_pid"
cp "$T/active.rel2" "$D/state/active"
expect 0 'fails its checks.*will be replaced' "rollback dry run explains the damaged live pair" -- bash "$S/rollback.sh"
expect 0 'on failure nothing serves' "rollback dry run says the damaged pair is not restarted" -- bash "$S/rollback.sh"
paused old-stopped bash "$S/rollback.sh" --apply
expect 2 'another deployment holds' "an explicit restore during the rollback is refused by the lock" -- \
    bash "$S/restore.sh" rel1 p237 --apply
resume
check "the rollback from a damaged rel2 exits 0" -- test "$rc" = 0
check "it reports rel2 as no rollback target" -- grep -q 'live health pass: rel1 with p237; rel2 with p275 failed its checks' "$T/bg.out"
check "rel1 with p237 serves, recorded, and nothing of rel2 is left" -- rel1_back
check "state/previous still names the intact rel1, not the damaged rel2" -- test "$(field previous label)" = rel1
expect 2 'already live' "a second rollback does not return to the damaged pair" -- bash "$S/rollback.sh" --apply

# The rollback target fails after a damaged pair stopped: nothing serves, and
# only an explicit restore of an intact pair brings the site back.
repair
expect 0 'live health pass: rel2 with p275' "activate the repaired rel2" -- bash "$S/activate.sh" rel2 p275 --apply
damage
expect 4 'was not restarted' "a failed rollback from a damaged pair leaves nothing serving (exit 4)" -- \
    env GRM_STUB_FAIL=/rel1/ bash "$S/rollback.sh" --apply
check "after the failed damaged rollback no rel1 or rel2 runs and no journal is left" -- \
    test -z "$(running rel1)$(running rel2)$(listener "$LP")" -a ! -e "$D/serve-$LP.starting"
impostor_up
expect 2 'is not grm' "rollback refuses to stop an impostor on the port" -- bash "$S/rollback.sh" --apply
expect 2 'is not grm' "an explicit restore refuses an impostor on the port" -- bash "$S/restore.sh" rel1 p237 --apply
check "the impostor still holds the port" -- test "$(listener "$LP")" = "$ipid"
kill "$ipid"
wait "$ipid" 2>/dev/null || true
wait_for 5 port_free "$LP" || true
expect 1 'does not match' "argument-free restore refuses the damaged pair state/active names" -- bash "$S/restore.sh" --apply
expect 0 'restore health pass: rel1 with p237' "an explicit restore of the intact previous pair recovers" -- \
    bash "$S/restore.sh" rel1 p237 --apply

# Killed after a damaged pair stopped: the same explicit restore recovers.
repair
expect 0 'live health pass: rel2 with p275' "activate the repaired rel2 again" -- bash "$S/activate.sh" rel2 p275 --apply
damage
at old-stopped KILL bash "$S/rollback.sh" --apply
check "a rollback killed after stopping the damaged pair leaves nothing serving" -- test "$rc" = 137 -a "$(live)" = none
expect 1 'does not match' "argument-free restore refuses the damaged pair after the kill" -- bash "$S/restore.sh" --apply
expect 0 'restore health pass: rel1 with p237' "an explicit restore of rel1 recovers after the kill" -- \
    bash "$S/restore.sh" rel1 p237 --apply
repair
chmod a-w "$index" "$D/release/rel2/dist"
expect 0 'dry-run: would record' "the repaired rel2 verifies again" -- bash "$S/activate.sh" rel2 p275

# A deleted running executable: the server is identified by its recorded
# start (serve-PORT.pid's pid, start time, token and spawn marker), never by
# the executable path, and is replaced but never restarted. A record that
# disagrees in any part, and an impostor without the marker, are refused with
# nothing stopped.
grm2=$D/release/rel2/grm
expect 0 'live health pass: rel2 with p275' "activate rel2 before deleting its running executable" -- bash "$S/activate.sh" rel2 p275 --apply
chmod u+w "$D/release/rel2"
rm "$grm2"
dpid=$(listener "$LP")
cp "$D/state/active" "$T/active.rel2"
cp "$D/serve-$LP.pid" "$T/pid.rel2"
check "the running rel2 executable is deleted" -- test "$(readlink "/proc/$dpid/exe")" = "$grm2 (deleted)"
expect 2 'deleted or replaced executable' "activate refuses a live server whose executable was deleted" -- \
    bash "$S/activate.sh" rel1 p237 --apply
expect 2 'deleted or replaced executable' "restore refuses it as the live service" -- bash "$S/restore.sh" rel1 p237 --apply
expect 0 'identified only by its recorded start' "status identifies it by its recorded start" -- bash "$S/status.sh"
forged_pid() {  # forged_pid NAME PATTERN SED: serve-PORT.pid edited by SED is refused, nothing stopped
    sed -E "$3" "$T/pid.rel2" >"$D/serve-$LP.pid"
    expect 2 "$2" "rollback refuses a deleted-executable server when serve-PORT.pid $1" -- bash "$S/rollback.sh" --apply
    check "with serve-PORT.pid that $1, the same pid serves" -- listens_as "$dpid"
}
forged_pid "records another start time" 'not identified by its recorded start: .*did not start at' 's/ start=[0-9]+/ start=1/'
forged_pid "records another pid" 'not identified by its recorded start: .*records pid' 's/^pid=[0-9]+/pid=1/'
forged_pid "records another token" 'not identified by its recorded start: .*no start token' \
    's/ token=[0-9a-f-]+/ token=00000000-0000-4000-8000-000000000000/'
forged_pid "records no marker" 'not identified by its recorded start: .*no start token' 's/ token=.*//'
forged_pid "records another marker" 'not identified by its recorded start: .*missing or replaced' 's/( mark=[0-9a-f-]+\.)[A-Za-z0-9]+/\1zzzzzz/'
forged_pid "records another pair" 'not identified by its recorded start: .*does not serve published/p237' 's/ publication=p275/ publication=p237/'
cp "$T/pid.rel2" "$D/serve-$LP.pid"
rm "$D/serve-$LP.pid"
expect 2 'not identified by its recorded start: serve-.*pid is missing' "rollback refuses a deleted-executable server without serve-PORT.pid" -- \
    bash "$S/rollback.sh" --apply
cp "$T/pid.rel2" "$D/serve-$LP.pid"
dmark=$(sed -n 's/.* mark=\([^ ]*\) .*/\1/p' "$D/serve-$LP.pid")
mv "$D/spawn/$dmark" "$T/dmark"
expect 2 'not identified by its recorded start: .*missing or replaced' "rollback refuses when the server's marker file is gone" -- \
    bash "$S/rollback.sh" --apply
mv "$T/dmark" "$D/spawn/$dmark"
check "after every refused identity, the deleted-executable server still serves" -- listens_as "$dpid"
expect 0 'recorded start identifies pid' "rollback dry run identifies the deleted-executable server by its recorded start" -- \
    bash "$S/rollback.sh"
expect 0 'live health pass: rel1 with p237; rel2 with p275 failed its checks' \
    "rollback replaces a server whose executable was deleted" -- bash "$S/rollback.sh" --apply
check "rel1 serves, recorded; nothing of rel2 is left" -- rel1_back
check "the deleted-executable server is gone and was not restarted" -- wait_for 10 gone "$dpid"
check "state/previous still names rel1, not the damaged rel2" -- test "$(field previous label)" = rel1
check "the stopped server's marker was removed" -- test ! -e "$D/spawn/$dmark"

# An impostor: grm with a deleted executable, running from the root with the
# recorded data path and a forged token naming an existing marker, recorded in
# a forged serve-PORT.pid, but not holding the marker. Rollback never stops it.
kill_live
mkdir "$T/imp-bin"
cp "$T/grm-stub" "$T/imp-bin/grm"
imark=$(basename "$(find "$D/spawn" -type f -name '*.*' -print -quit)")
(cd "$D" && GRM_DEPLOY_SPAWN=${imark%.*} exec "$T/imp-bin/grm" --data "$D/published/p275" serve \
    --app "$D/release/rel2/dist" --port "$LP" >/dev/null 2>&1 </dev/null) &
ipid=$!
wait_for 10 listens_as "$ipid" || true
rm "$T/imp-bin/grm"
printf 'pid=%s start=%s label=rel2 publication=p275 token=%s mark=%s mark_id=%s\n' "$ipid" "$(start_of "$ipid")" "${imark%.*}" "$imark" \
    "$(stat -c %d:%i "$D/spawn/$imark")" >"$D/serve-$LP.pid"
cp "$T/active.rel2" "$D/state/active"
expect 2 'not identified by its recorded start: .*does not hold spawn/' "rollback refuses an impostor that does not hold the recorded marker" -- \
    bash "$S/rollback.sh" --apply
check "the impostor was not signalled" -- listens_as "$ipid"
kill -KILL "$ipid"
wait "$ipid" 2>/dev/null || true
wait_for 5 port_free "$LP" || true
rm "$D/serve-$LP.pid"
expect 0 'restore health pass: rel1 with p237' "an explicit restore of rel1 recovers after the impostor" -- bash "$S/restore.sh" rel1 p237 --apply
cp "$T/grm-stub" "$grm2"
chmod a-w "$grm2" "$D/release/rel2"
expect 0 'dry-run: would record' "rel2 with its executable restored verifies again" -- bash "$S/activate.sh" rel2 p275

# A deleted dist/ while it serves: the release fails its checks, the record
# identifies the server, and rollback replaces it without restarting it.
expect 0 'live health pass: rel2 with p275' "activate rel2 before deleting its dist/" -- bash "$S/activate.sh" rel2 p275 --apply
dpid=$(listener "$LP")
chmod u+w "$D/release/rel2" "$D/release/rel2/dist"
mv "$D/release/rel2/dist" "$T/dist.rel2"
expect 0 'live health pass: rel1 with p237; rel2 with p275 failed its checks' \
    "rollback replaces a server whose dist/ was deleted" -- bash "$S/rollback.sh" --apply
check "after the deleted dist/, rel1 is back and the old server is gone" -- eval 'rel1_back && wait_for 10 gone "$dpid"'
mv "$T/dist.rel2" "$D/release/rel2/dist"
chmod a-w "$D/release/rel2/dist" "$D/release/rel2"
expect 0 'dry-run: would record' "rel2 with its dist/ restored verifies again" -- bash "$S/activate.sh" rel2 p275

# A candidate that ignores SIGTERM fails staging and is still killed; an
# interrupted stage stops its candidate too.
expect 1 'stage rel3 with p237 failed' "stage fails a candidate that ignores SIGTERM" -- \
    env GRM_EXPECT_DATASETS=237 GRM_STUB_IGNORE_TERM=/rel3/ bash "$S/stage.sh" rel3 p237 "$SP" --apply
check "the TERM-ignoring candidate is gone and the port is free" -- test -z "$(running rel3)$(listener "$SP")"
GRM_EXPECT_DATASETS=237 GRM_STUB_SLOW=/rel3/ bash "$S/stage.sh" rel3 p237 "$SP" --apply >"$T/stage.out" 2>&1 &
spid=$!
check "staging starts its candidate" -- wait_for 10 has_running rel3
kill -TERM "$spid" || true
set +e
wait "$spid"
rc=$?
set -e
check "an interrupted stage exits non-zero" -- test "$rc" != 0
check "an interrupted stage leaves no candidate or listener" -- test -z "$(running rel3)$(listener "$SP")"
check "an interrupted stage writes no approval" -- test ! -e "$D/go/rel3--p237"

# Staging is supervised the same way: TERM at any checkpoint stops the
# candidate, and a stage killed outright leaves a journal that the next
# stage.sh --apply on that port resolves.
stage3() { env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" "$@"; }
for point in journaled spawned registered assigned healthy; do
    at "$point" TERM env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
    check "stage TERM at $point exits non-zero" -- test "$rc" != 0
    check "stage TERM at $point leaves no candidate, listener, journal or approval" -- \
        test -z "$(running rel3)$(listener "$SP")" -a ! -e "$D/serve-$SP.starting" -a ! -e "$D/go/rel3--p237"
done
at spawned KILL env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
m=$(jmark "$SP")
check "after stage KILL at spawned the unregistered child exits without exec'ing" -- wait_for 10 eval '! locked "$D/spawn/$m"'
check "after stage KILL at spawned no candidate runs or listens" -- test -z "$(running rel3)$(listener "$SP")"
expect 0 'dry-run: would stop pids none' "stage dry run reports the journal of KILL at spawned" -- stage3
expect 0 'stage pass' "stage --apply resolves the journal of KILL at spawned and stages afresh" -- stage3 --apply
check "the journal and marker of stage KILL at spawned are gone" -- test ! -e "$D/serve-$SP.starting" -a ! -e "$D/spawn/$m"
rm "$D/go/rel3--p237"
port=$SP
preexec child-preexec env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
preexec_resolved "stage killed with its child paused before exec" 'stage pass' stage3 --apply
check "after the stage's paused child was resumed, no candidate or listener is left" -- \
    test -z "$(running rel3)$(listener "$SP")" -a ! -e "$D/serve-$SP.starting"
rm "$D/go/rel3--p237"
for point in healthy; do
    at "$point" KILL env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
    check "after stage KILL at $point its orphaned candidate runs" -- wait_for 10 has_running rel3
    orphan=$(running rel3)
    expect 0 'dry-run: would stop pids [0-9]' "stage dry run reports the orphan of KILL at $point" -- stage3
    check "the stage dry run left the orphan of KILL at $point alone" -- test "$(running rel3)" = "$orphan"
    expect 0 'stage pass' "stage --apply stops the orphan of KILL at $point and stages afresh" -- stage3 --apply
    check "the orphan of stage KILL at $point is gone and the port is free" -- \
        test ! -e "/proc/$orphan/exe" -a -z "$(running rel3)$(listener "$SP")" -a ! -e "$D/serve-$SP.starting"
    rm "$D/go/rel3--p237"
done
# Receipt directories never collide, even with a stage in the same second.
now=$(date -u +%s)
for ((i = 0; i < 4; i++)); do mkdir -p "$D/stage/rel3--p237--$(date -u -d "@$((now + i))" +%Y%m%dT%H%M%SZ)"; done
expect 0 'stage pass' "stage passes while receipt directories for this second already exist" -- stage3 --apply
rm "$D/go/rel3--p237"
mkdir "$D/go/rel3--p237"
expect 2 'link, directory or special file' "stage refuses an approval path that is a directory" -- stage3 --apply
check "nothing was moved into the approval directory" -- test -z "$(ls -A "$D/go/rel3--p237")"
rmdir "$D/go/rel3--p237"
paused healthy env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
mkdir "$D/go/rel3--p237"
resume
check "an approval path made a directory during the stage fails it" -- test "$rc" = 1 -a -z "$(ls -A "$D/go/rel3--p237")"
check "that stage leaves no candidate or listener" -- test -z "$(running rel3)$(listener "$SP")"
rmdir "$D/go/rel3--p237"

# The browser check that stage runs is the release's own copy, which its
# manifest binds to its commit; a different script in src-LABEL never runs.
# node is a stub that records the script it was given.
mkdir -p "$T/nodebin" "$T/playwright" "$D/src-rel3/web/scripts"
cat >"$T/nodebin/node" <<'SH'
#!/bin/sh
printf '%s\n' "$1" >>"$GRM_TEST_ROOT/node.runs"
exit "${GRM_STUB_NODE_RC:-0}"
SH
chmod +x "$T/nodebin/node"
echo 'process.exit(0)' >"$D/src-rel3/web/scripts/browser-check.cjs"
browser3() {  # browser3 [NAME=VALUE...]: stage rel3 with the browser check
    env -u GRM_SKIP_BROWSER PATH="$T/nodebin:$PATH" GRM_PLAYWRIGHT_CORE="$T/playwright" GRM_EXPECT_DATASETS=237 "$@" \
        bash "$S/stage.sh" rel3 p237 "$SP" --apply
}
expect 0 'stage pass: approved rel3 with p237 \(browser=pass\)' "stage runs the browser check and approves browser=pass" -- browser3
check "the check that ran is the release's manifest-bound copy, not src-rel3's" -- \
    test "$(cat "$T/node.runs")" = "$D/release/rel3/browser-check.cjs"
check "the approval records browser=pass with the SHA-256 of the check that ran" -- \
    grep -q "^stage=0 browser=pass browser_check=$(sha256sum <"$D/release/rel3/browser-check.cjs" | cut -d' ' -f1) " "$D/go/rel3--p237"
expect 0 'dry-run: would record' "activation accepts that browser=pass approval" -- bash "$S/activate.sh" rel3 p237
rm "$D/go/rel3--p237"
expect 1 'stage rel3 with p237 failed' "a failing browser check fails the stage" -- browser3 GRM_STUB_NODE_RC=1
check "the failed browser check wrote no approval and left no candidate" -- \
    test ! -e "$D/go/rel3--p237" -a -z "$(running rel3)$(listener "$SP")"

# A reused PID: in a new user and PID namespace, where this test may choose the
# next PID, stage's candidate exits at once (GRM_STUB_FAIL), stage reaps it and
# stops itself at "healthy", and a new process takes the candidate's PID.
# Resumed, stage must fail without signalling that process.
cat >"$T/pid-reuse.sh" <<'SH'
S=$1 SP=$2
env GRM_EXPECT_DATASETS=237 GRM_STUB_FAIL=/rel3/ GRM_TEST_SIGNAL_AT=healthy GRM_TEST_SIGNAL=STOP \
    bash "$S/stage.sh" rel3 p237 "$SP" --apply >"$GRM_TEST_ROOT/reuse.out" 2>&1 &
stage=$!
state() { sed 's/.*) //' "/proc/$1/stat" 2>/dev/null | cut -d' ' -f1; }
for _ in $(seq 300); do [ "$(state "$stage")" = T ] && break; sleep 0.1; done
cand=$(sed -n 's/.* child=\([0-9]*\) .*/\1/p' "$GRM_DEPLOY_ROOT/serve-$SP.starting" 2>/dev/null)
for _ in $(seq 100); do [ -n "$cand" ] && [ ! -e "/proc/$cand" ] && break; sleep 0.1; done
if [ "$(state "$stage")" != T ] || [ -z "$cand" ] || [ -e "/proc/$cand" ]; then
    echo "reuse=unset stage=$(state "$stage") cand=$cand"
    kill -KILL "$stage"
    exit 1
fi
echo "$((cand - 1))" >/proc/sys/kernel/ns_last_pid
sleep 300 &
imp=$!
if [ "$imp" != "$cand" ]; then
    echo "reuse=missed cand=$cand got=$imp"
    kill -KILL "$stage" "$imp"
    exit 1
fi
kill -CONT "$stage"
wait "$stage"
rc=$?
case $(state "$imp") in '' | Z) imp_state=gone ;; *) imp_state=alive ;; esac
kill -KILL "$imp" 2>/dev/null
echo "reuse=ok stage_rc=$rc impostor=$imp_state"
SH
# A reused process group number (the 2026-10-28 review's finding 1 and its
# probe, resolve-reuse-probe.sh), through common.sh's pending_start and
# resolve_start, which restore.sh --apply and stage.sh --apply run for an
# interrupted start. The start is journaled and its registered child, PID
# 4000, exits before it ever leads a group. An unrelated daemon then
# double-forks: its intermediate takes PID 4000, calls setsid, forks a worker
# and exits, so group 4000 has a member and no leader. resolve_start must
# signal none of it, keep the journal and fail, saying why.
cat >"$T/group-reuse.sh" <<'SH'
set -u
S=$1 ROOT=$2 N=4000
export GRM_DEPLOY_ROOT=$ROOT GRM_PORT=8999
mkdir -p "$ROOT/spawn"
source "$S/common.sh"
set +e
SPAWN_TOKEN=$(cat /proc/sys/kernel/random/uuid)
mark_create "$D/spawn" "$SPAWN_TOKEN" || { echo "no marker"; exit 3; }
SPAWN_MARK=$MARK_PATH SPAWN_MARK_ID=$MARK_ID
sleep 0.1 &
owner=$!
ostart=$(proc_start "$owner")
echo $((N - 1)) >/proc/sys/kernel/ns_last_pid
sleep 0.2 &
child=$!
cstart=$(proc_start "$child")
wait "$owner" "$child"
close_fd MARK_FD
printf 'token=%s owner=%s owner_start=%s label=rel1 publication=p1 port=8999 mark=%s mark_id=%s child=%s child_start=%s\n' \
    "$SPAWN_TOKEN" "$owner" "$ostart" "${SPAWN_MARK##*/}" "$SPAWN_MARK_ID" "$child" "$cstart" >"$D/serve-8999.starting"
echo $((N - 1)) >/proc/sys/kernel/ns_last_pid
setsid bash -c "sleep 300 </dev/null >/dev/null 2>&1 & echo \$! >$ROOT/worker.pid; exit 0"
for _ in $(seq 50); do [ -s "$ROOT/worker.pid" ] && break; sleep 0.1; done
worker=$(cat "$ROOT/worker.pid")
wstart=$(proc_start "$worker")
if [ "$child" != "$N" ] || ! proc_stat "$worker" || [ "$PS_PGRP" != "$N" ] || [ -e "/proc/$N" ]; then
    echo "setup failed: child=$child worker=$worker pgid=${PS_PGRP:-}"
    kill -KILL "$worker" 2>/dev/null
    exit 3
fi
pending_start 8999
( resolve_start 8999 ) >"$ROOT/resolve.out" 2>&1
rc=$?
if alive "$worker" "$wstart"; then w=alive; else w=KILLED; fi
kill -KILL "$worker" 2>/dev/null
if grep -q "process group $N has members (pids $worker) that cannot be tied to this start" "$ROOT/resolve.out"; then why=untied; else why=other; fi
echo "resolve_rc=$rc journal_left=$([ -e "$D/serve-8999.starting" ] && echo yes || echo no) unrelated_worker=$w reason=$why"
SH
# A process outside the start with the marker's inode open (the 2026-10-29
# review of the regression fix, finding 1), through pending_start and
# resolve_start. The start is journaled and its processes are gone, so the
# marker's lock is free. "reader": unrelated work (a backup, an indexer, an
# operator's less) has the intact marker open without its lock; resolve_start
# must signal it, resolve the journal and report no process of the start.
# "replaced": the marker is gone while its inode is still held, as when the
# marker was removed and the freed inode went to another file that unrelated
# work (the next start, say) holds locked; a reused inode cannot be forced on
# every filesystem, so here the unrelated process opens the marker and takes
# its own lock before the marker is removed. resolve_start must signal none of
# it, keep the journal and fail, saying the marker is gone.
cat >"$T/marker-holder.sh" <<'SH'
set -u
S=$1 ROOT=$2 mode=$3
export GRM_DEPLOY_ROOT=$ROOT GRM_PORT=8999
mkdir -p "$ROOT/spawn"
source "$S/common.sh"
set +e
SPAWN_TOKEN=$(cat /proc/sys/kernel/random/uuid)
mark_create "$D/spawn" "$SPAWN_TOKEN" || { echo "no marker"; exit 3; }
sleep 0.1 &
owner=$!
ostart=$(proc_start "$owner")
wait "$owner"
close_fd MARK_FD
printf 'token=%s owner=%s owner_start=%s label=rel1 publication=p1 port=8999 mark=%s mark_id=%s\n' \
    "$SPAWN_TOKEN" "$owner" "$ostart" "${MARK_PATH##*/}" "$MARK_ID" >"$D/serve-8999.starting"
if [ "$mode" = reader ]; then
    setsid bash -c 'exec sleep 300 3<"$1"' _ "$MARK_PATH" </dev/null >/dev/null 2>&1 &
else
    setsid bash -c 'exec 3<"$1"; flock -n 3 || exit 9; exec sleep 300' _ "$MARK_PATH" </dev/null >/dev/null 2>&1 &
fi
w=$!
for _ in $(seq 100); do [ "$(proc_comm "$w")" = sleep ] && holds "$w" "$MARK_ID" && break; sleep 0.05; done
wstart=$(proc_start "$w")
lk=no
for fd in $(fds_on "$w" "$MARK_ID"); do fd_holds_flock "$w" "$fd" && lk=yes; done
if ! holds "$w" "$MARK_ID" || [ "$lk" != "$([ "$mode" = reader ] && echo no || echo yes)" ]; then
    echo "setup failed: pid $w holds=$(holds "$w" "$MARK_ID" && echo yes || echo no) lock=$lk"
    kill -KILL "$w" 2>/dev/null
    exit 3
fi
[ "$mode" = reader ] || rm -f -- "$MARK_PATH"
pending_start 8999
( resolve_start 8999 ) >"$ROOT/resolve.out" 2>&1
rc=$?
if alive "$w" "$wstart"; then u=alive; else u=KILLED; fi
kill -KILL "$w" 2>/dev/null
if grep -q 'is missing or replaced, so its lock cannot show' "$ROOT/resolve.out"; then why=replaced; else why=other; fi
echo "resolve_rc=$rc journal_left=$([ -e "$D/serve-8999.starting" ] && echo yes || echo no) unrelated=$u orphans=[$ORPHANS] reason=$why"
SH
# A holder whose PID is reused while mark_stop stops an earlier holder (the
# 2026-10-29 review of 742cd4f, finding 1). Two holders of a start's lock:
# H1 ignores TERM for 6 s, H2 exits after 2 s. mark_stop lists both, then
# spends about 6 s on H1; meanwhile unrelated work Q (own session, no marker,
# no token) takes H2's PID. mark_stop must leave Q alone.
cat >"$T/pid-recycle.sh" <<'SH'
set -u
S=$1 dir=$2
source "$S/lib.sh"
mark_create "$dir" recycle || { echo "no marker"; exit 3; }
echo 199 >/proc/sys/kernel/ns_last_pid
setsid bash -c 'trap "" TERM; sleep 6; :' </dev/null >/dev/null 2>&1 &
h1=$!
echo 299 >/proc/sys/kernel/ns_last_pid
setsid sleep 2 </dev/null >/dev/null 2>&1 &
h2=$!
close_fd MARK_FD
sleep 0.3
holds "$h1" "$MARK_ID" && holds "$h2" "$MARK_ID" || { echo "setup=failed: the holders do not hold the marker"; exit 3; }
( mark_stop "$MARK_ID" "$MARK_PATH"; echo "rc=$?" >"$dir/ms.out" ) &
ms=$!
wait "$h2"
echo $((h2 - 1)) >/proc/sys/kernel/ns_last_pid
setsid sleep 300 </dev/null >/dev/null 2>&1 &
q=$!
qs=$(proc_start "$q")
wait "$ms"
if alive "$q" "$qs"; then u=alive; else u=KILLED; fi
{ kill -KILL "$q"; wait "$q"; } 2>/dev/null
setup=ok
[ "$q" = "$h2" ] || setup="failed: Q got pid $q, not $h2"
echo "setup=$setup mark_stop_$(cat "$dir/ms.out") unrelated=$u"
SH
# A process of the start that leaves the group and converts its inherited lock
# to shared (the same review, finding 2). It holds no exclusive lock, and the
# shared probe of lock_free does not conflict with it, so mark_stop must not
# succeed: it signals nothing, fails and names it.
cat >"$T/lock-shared.sh" <<'SH'
set -u
S=$1 dir=$2
source "$S/lib.sh"
mark_create "$dir" change || { echo "no marker"; exit 3; }
setsid bash -c 'flock -s "$1" || exit 9; exec sleep 300' _ "$MARK_FD" </dev/null >/dev/null 2>&1 &
c=$!
close_fd MARK_FD
for _ in $(seq 100); do [ "$(proc_comm "$c")" = sleep ] && break; sleep 0.05; done
cs=$(proc_start "$c")
lk=none
for fd in $(fds_on "$c" "$MARK_ID"); do lk=$(awk '$1 == "lock:" {print $5}' "/proc/$c/fdinfo/$fd"); lk=${lk:-none}; done
mark_stop "$MARK_ID" "$MARK_PATH"
rc=$?
if alive "$c" "$cs"; then s=alive; else s=KILLED; fi
kill -KILL "$c" 2>/dev/null
case $MARK_ERR in "pids $c hold a lock on "*" that is not the start's exclusive lock"*) why=named ;; *) why=other ;; esac
echo "child_lock=$lk mark_stop_rc=$rc child=$s reason=$why"
SH
# The marker goes during mark_stop's scan for holders (the 742cd4f review,
# F5): its inode may then be another file's, so a holder found in that scan
# must not be signalled. Here unrelated work holds the marker's inode with its
# own exclusive lock (the start's is free), and mark_holders is wrapped so that
# its first call removes the marker, between mark_stop's check of the marker
# and the end of the scan. mark_stop must leave it alone and fail.
cat >"$T/marker-scan.sh" <<'SH'
set -u
S=$1 dir=$2
source "$S/lib.sh"
mark_create "$dir" scan || { echo "no marker"; exit 3; }
close_fd MARK_FD
setsid bash -c 'exec 3<"$1"; flock -n 3 || exit 9; exec sleep 300' _ "$MARK_PATH" </dev/null >/dev/null 2>&1 &
w=$!
for _ in $(seq 100); do [ "$(proc_comm "$w")" = sleep ] && break; sleep 0.05; done
ws=$(proc_start "$w")
holds "$w" "$MARK_ID" || { echo "setup failed: pid $w does not hold the marker"; kill -KILL "$w"; exit 3; }
eval "$(declare -f mark_holders | sed '1s/mark_holders/real_mark_holders/')"
mark_holders() { rm -f -- "$MARK_PATH"; real_mark_holders "$@"; }
mark_stop "$MARK_ID" "$MARK_PATH"
rc=$?
if alive "$w" "$ws"; then u=alive; else u=KILLED; fi
{ kill -KILL "$w"; wait "$w"; } 2>/dev/null
case $MARK_ERR in *"is missing or replaced"*) why=replaced ;; *) why=other ;; esac
echo "mark_stop_rc=$rc unrelated=$u reason=$why"
SH
# A process of the start (it carries the start's token) that leaves the group
# after closing the marker ("closed") or releasing the inherited lock with
# flock -u ("unlocked"), which releases it for every process sharing it (the
# 742cd4f review, F3). Nothing ties it to the start but the token, so
# mark_stop must not succeed: it signals nothing, fails and names it.
cat >"$T/token-left.sh" <<'SH'
set -u
S=$1 dir=$2 mode=$3
source "$S/lib.sh"
tok=$(cat /proc/sys/kernel/random/uuid)
mark_create "$dir" tokleft || { echo "no marker"; exit 3; }
if [ "$mode" = closed ]; then
    env GRM_DEPLOY_SPAWN="$tok" setsid bash -c 'eval "exec $1>&-"; exec sleep 300' _ "$MARK_FD" </dev/null >/dev/null 2>&1 &
else
    env GRM_DEPLOY_SPAWN="$tok" setsid bash -c 'flock -u "$1" || exit 9; exec sleep 300' _ "$MARK_FD" </dev/null >/dev/null 2>&1 &
fi
c=$!
for _ in $(seq 100); do [ "$(proc_comm "$c")" = sleep ] && break; sleep 0.05; done
close_fd MARK_FD
cs=$(proc_start "$c")
held=no; holds "$c" "$MARK_ID" && held=yes
mark_stop "$MARK_ID" "$MARK_PATH" '' '' "$tok"
rc=$?
if alive "$c" "$cs"; then s=alive; else s=KILLED; fi
{ kill -KILL "$c"; wait "$c"; } 2>/dev/null
case $MARK_ERR in "pids $c carry this start's token"*) why=token ;; *) why=other ;; esac
echo "mode=$mode holds_marker=$held mark_stop_rc=$rc child=$s reason=$why"
SH
pidns=(unshare --user --map-root-user --pid --fork --mount-proc)
skipped=0
if "${pidns[@]}" sh -c 'echo 1 >/proc/sys/kernel/ns_last_pid' 2>/dev/null; then
    expect 0 '^reuse=ok stage_rc=1 impostor=alive$' "stage fails, and signals no process that took its exited candidate's PID" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/pid-reuse.sh" "$S" "$SP"
    check "that stage left no journal, listener or approval" -- \
        test ! -e "$D/serve-$SP.starting" -a -z "$(listener "$SP")" -a ! -e "$D/go/rel3--p237"
    check "its receipts say the candidate exited before the checks finished" -- \
        grep -rq 'exited before the checks finished' "$D/stage"
    expect 0 '^resolve_rc=1 journal_left=yes unrelated_worker=alive reason=untied$' \
        "resolving a start whose child's number now leads an unrelated group signals none of it and keeps the journal" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/group-reuse.sh" "$S" "$T/group-reuse"
    expect 0 '^resolve_rc=0 journal_left=no unrelated=alive orphans=\[\] reason=other$' \
        "resolving a start signals no process that merely has its intact marker open, and resolves the journal" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/marker-holder.sh" "$S" "$T/marker-reader" reader
    expect 0 '^resolve_rc=1 journal_left=yes unrelated=alive orphans=\[\] reason=replaced$' \
        "resolving a start whose marker is gone signals no process that holds the old inode locked, and keeps the journal" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/marker-holder.sh" "$S" "$T/marker-replaced" replaced
    expect 0 '^setup=ok mark_stop_rc=0 unrelated=alive$' \
        "mark_stop signals no process that took the PID of a holder that exited while an earlier holder was being stopped" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/pid-recycle.sh" "$S" "$T/pid-recycle"
    expect 0 '^child_lock=READ mark_stop_rc=1 child=alive reason=named$' \
        "mark_stop fails, signalling nothing, while a process of the start that left the group holds the lock as shared" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/lock-shared.sh" "$S" "$T/lock-shared"
    expect 0 '^mark_stop_rc=1 unrelated=alive reason=replaced$' \
        "mark_stop signals no holder found in a scan during which the marker went, and fails" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/marker-scan.sh" "$S" "$T/marker-scan"
    expect 0 '^mode=closed holds_marker=no mark_stop_rc=1 child=alive reason=token$' \
        "mark_stop fails, signalling nothing, while a process of the start that closed the marker and left the group carries its token" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/token-left.sh" "$S" "$T/token-closed" closed
    expect 0 '^mode=unlocked holds_marker=yes mark_stop_rc=1 child=alive reason=token$' \
        "mark_stop fails, signalling nothing, while a process of the start that released the lock and left the group carries its token" -- \
        timeout -k 5 120 "${pidns[@]}" bash "$T/token-left.sh" "$S" "$T/token-unlocked" unlocked
else
    echo "skip stage with a reused candidate PID, resolve with a reused group number or an unrelated marker holder, a reused holder PID, a shared lock, a marker gone during a scan and a token left behind (no unprivileged user and PID namespaces)"
    skipped=$((skipped + 11))
fi

# A candidate that passes every check and then exits before the stop fails
# the stage: it is killed while stage waits after its health check.
paused healthy env GRM_EXPECT_DATASETS=237 bash "$S/stage.sh" rel3 p237 "$SP" --apply
cand=$(running rel3)
check "the paused stage's candidate passed its checks and runs" -- test -n "$cand"
[ -z "$cand" ] || kill -KILL "$cand" 2>/dev/null || true
wait_for 10 eval 'test -z "$(running rel3)"' || true
resume
check "a candidate that exits after passing its checks fails the stage (exit 1)" -- test "$rc" = 1
receipts=$(sed -n 's/.*failed; receipts in //p' "$T/bg.out")
check "that stage wrote no approval, and its receipts say the candidate exited" -- \
    eval 'test ! -e "$D/go/rel3--p237" -a -n "$receipts" && grep -q "exited before the checks finished" "$receipts/stop.log"'

# Process groups, through lib.sh directly: a registered leader started under
# setsid with the start's token (as start_service starts grm) forks a
# descendant that closes the spawn marker, ignores SIGTERM and stays in the
# leader's group. mark_stop must kill it when the leader holds the marker and
# exits on SIGTERM, even without the token, because it was seen in the group
# while the leader lived ("running"), and when the leader has already exited
# and the descendant carries the token ("exited"). With the leader gone and
# neither token nor marker, nothing ties the group to the start, so mark_stop
# signals no member and fails ("untied"). A tie holds for that process: one
# that carries the token when the stop begins and then execs without it, a
# second later, is still killed ("sticky"). mark_stop must also fail when a
# stop failed even if nothing holds the marker by its final check ("failing").
cat >"$T/group-stop.sh" <<'SH'
set -u
source "$1/lib.sh"
dir=$2 mode=$3
tok=$(cat /proc/sys/kernel/random/uuid)
mark_create "$dir" grp || { echo "no marker"; exit 3; }
if [ "$mode" = failing ]; then
    sleep 300 &
    close_fd MARK_FD
    stop_tree() { kill -KILL "$1" 2>/dev/null; sleep 0.5; return 1; }
    mark_stop "$MARK_ID" "$MARK_PATH"
    echo "rc=$? err=$MARK_ERR"
    exit 0
fi
# The sticky descendant: python, until leader.go2 exists; then sleep, without the token.
cat >"$dir/drop-token.py" <<'PY'
import os, sys, time
while not os.path.exists(sys.argv[1]):
    time.sleep(0.05)
env = {k: v for k, v in os.environ.items() if k != "GRM_DEPLOY_SPAWN"}
os.execvpe("sleep", ["sleep", "300"], env)
PY
( exec env GRM_DEPLOY_SPAWN="$tok" setsid bash -c '
    case $2 in running | untied) unset GRM_DEPLOY_SPAWN ;; esac
    if [ "$2" = sticky ]; then
        (trap "" TERM; exec python3 -I "$3" "$1.go2" '"$MARK_FD"'>&-) &
    else
        (trap "" TERM; exec sleep 300 '"$MARK_FD"'>&-) &
    fi
    echo "$!" >"$1.kid"
    while [ ! -e "$1.go" ]; do sleep 0.05; done
    [ "$2" = running ] || exit 0
    wait' _ "$dir/leader" "$mode" "$dir/drop-token.py" ) &
leader=$!
for _ in $(seq 100); do [ -s "$dir/leader.kid" ] && break; sleep 0.05; done
kid=$(cat "$dir/leader.kid")
lstart=$(proc_start "$leader")
kstart=$(proc_start "$kid")
want=sleep
if [ "$mode" = sticky ]; then want=python3; fi
for _ in $(seq 100); do ! holds "$kid" "$MARK_ID" && [ "$(proc_comm "$kid")" = "$want" ] && break; sleep 0.05; done
proc_stat "$kid" && [ "$PS_PGRP" = "$leader" ] || { echo "the descendant is not in the leader's group"; exit 3; }
: >"$dir/leader.go"
if [ "$mode" != running ]; then wait "$leader"; fi
close_fd MARK_FD
if [ "$mode" = sticky ]; then (sleep 1; : >"$dir/leader.go2") & fi
mark_stop "$MARK_ID" "$MARK_PATH" "$leader" "$lstart" "$tok"
rc=$?
if alive "$kid" "$kstart"; then k=alive; else k=gone; fi
kill -KILL "$kid" 2>/dev/null
echo "rc=$rc kid=$k err=$MARK_ERR"
SH
expect 0 '^rc=0 kid=gone err=$' "mark_stop kills a TERM-ignoring descendant seen in the group of a leader that exits on SIGTERM" -- \
    timeout -k 5 90 bash "$T/group-stop.sh" "$S" "$T/groups-running" running
expect 0 '^rc=0 kid=gone err=$' "mark_stop kills a TERM-ignoring descendant carrying the token in the group of a registered leader that had exited" -- \
    timeout -k 5 90 bash "$T/group-stop.sh" "$S" "$T/groups-exited" exited
expect 0 '^rc=1 kid=alive err=process group [0-9]+ has members \(pids [0-9]+\) that cannot be tied to this start' \
    "mark_stop signals no member of an exited leader's group that has neither the token nor the marker, and fails" -- \
    timeout -k 5 90 bash "$T/group-stop.sh" "$S" "$T/groups-untied" untied
expect 0 '^rc=0 kid=gone err=$' "mark_stop still kills a TERM-ignoring descendant that carried the token and then exec'd without it" -- \
    timeout -k 5 90 bash "$T/group-stop.sh" "$S" "$T/groups-sticky" sticky
expect 0 '^rc=1 err=pid [0-9]+ would not stop$' "mark_stop fails after a failed stop even when no holder is left" -- \
    timeout -k 5 60 bash "$T/group-stop.sh" "$S" "$T/groups-failing" failing
# A lock_free probe that holds the marker's lock at the instant a mark_stop
# scans (the 2026-10-29 re-check of the marker fix, finding 1). The start's
# processes are gone, so its lock is free; a probe (connector.sh status, say,
# from a shell whose group it leads) takes it for an instant. Here the probe
# is frozen while it holds it, through a flock wrapper; mark_stop must not
# take it for a holder of the start's lock or signal it. Frozen, it holds its
# shared lock through mark_stop's steady check, so mark_stop fails closed.
cat >"$T/lock-probe.sh" <<'SH'
set -u
source "$1/lib.sh"
dir=$2
mark_create "$dir" probe || { echo "no marker"; exit 3; }
close_fd MARK_FD
( flock() { command flock "$@" && kill -s STOP "$BASHPID"; }; lock_free "$MARK_PATH" ) &
p=$!
for _ in $(seq 100); do [ "$(awk '{print $3}' "/proc/$p/stat" 2>/dev/null)" = T ] && break; sleep 0.05; done
pstart=$(proc_start "$p")
mode=none
for fd in $(fds_on "$p" "$MARK_ID"); do mode=$(awk '$1 == "lock:" {print $5}' "/proc/$p/fdinfo/$fd"); done
mark_stop "$MARK_ID" "$MARK_PATH"
rc=$?
if alive "$p" "$pstart"; then s=alive; else s=KILLED; fi
kill -KILL "$p" 2>/dev/null
echo "probe_lock=$mode mark_stop_rc=$rc probe=$s"
SH
expect 0 '^probe_lock=READ mark_stop_rc=1 probe=alive$' \
    "mark_stop signals no lock_free probe that holds the free marker lock (frozen, so it fails closed)" -- \
    timeout -k 5 90 bash "$T/lock-probe.sh" "$S" "$T/lock-probe"

# Concurrent builds under one label: one builds, the other refuses at once.
G=$T/source
mkdir -p "$G/web/scripts" "$T/bin"
echo '{}' >"$G/web/package.json"
echo '// browser check of the source commit' >"$G/web/scripts/browser-check.cjs"
git -C "$G" init -q
git -C "$G" add -A
git -C "$G" -c user.name=test -c user.email=test@example.invalid commit -qm source
cat >"$T/bin/cargo" <<'SH'
#!/bin/sh
echo "$$" >>"$GRM_TEST_ROOT/cargo-runs"
sleep 1
mkdir -p "$CARGO_TARGET_DIR/release"
printf '#!/bin/sh\n' >"$CARGO_TARGET_DIR/release/grm"
chmod +x "$CARGO_TARGET_DIR/release/grm"
SH
# npm ci stands in for a dependency's install script that rewrites the
# extracted browser check; the release must still carry the commit's.
cat >"$T/bin/npm" <<'SH'
#!/bin/sh
if [ "$1" = ci ]; then echo '// rewritten during npm ci' >scripts/browser-check.cjs; fi
if [ "$1" = run ]; then mkdir -p dist/assets && echo '<!doctype html>' >dist/index.html && echo x >dist/assets/index-a.js; fi
SH
chmod +x "$T/bin/cargo" "$T/bin/npm"
build() { env GRM_SOURCE_REPO="$G" GRM_EXTRA_PATH="$T/bin" GRM_MIN_FREE_GIB=0 bash "$S/build-release.sh" "$@"; }
build relB HEAD --apply >"$T/b1.out" 2>&1 &
b1=$!
build relB HEAD --apply >"$T/b2.out" 2>&1 &
b2=$!
set +e
wait "$b1"
r1=$?
wait "$b2"
r2=$?
set -e
check "of two concurrent builds of one label, one succeeds and one refuses" -- \
    test "$(printf '%s\n' "$r1" "$r2" | sort | xargs)" = "0 2"
check "the refused build reports the reserved label" -- grep -q 'already exists' "$T/b1.out" "$T/b2.out"
check "only one build ran cargo" -- test "$(wc -l <"$T/cargo-runs")" = 1
check "the built release verifies" -- bash -c 'source "$1/common.sh"; verify_release relB' _ "$S"
check "the release records the built commit" -- test "$(cat "$D/release/relB/COMMIT")" = "$(git -C "$G" rev-parse HEAD)"
check "npm ci rewrote the extracted browser check (the control)" -- \
    grep -q 'rewritten during npm ci' "$D/src-relB/web/scripts/browser-check.cjs"
check "the release carries the commit's browser check, listed in its manifest" -- \
    eval 'cmp -s "$G/web/scripts/browser-check.cjs" "$D/release/relB/browser-check.cjs" && grep -q " \./browser-check\.cjs$" "$D/release/relB/MANIFEST.sha256"'

# Connector supervisors: the lock stays with the supervisor, and stop signals only
# the supervisor of its own directory.
make_connector() {
    mkdir -m 700 "$1"
    printf '#!/bin/sh\ncase "$1" in --version) echo "cloudflared version test"; exit 0 ;; esac\nexec sleep 300\n' >"$1/cloudflared"
    chmod 700 "$1/cloudflared"
    (cd "$1" && sha256sum cloudflared >cloudflared.sha256)
    printf 'placeholder-not-a-credential\n' >"$1/tunnel.token"
    chmod 600 "$1/tunnel.token"
}
C1=$T/conn1
C2=$T/conn2
make_connector "$C1"
make_connector "$C2"
connector() { env GRM_CONNECTOR_DIR="$1" bash "$S/connector.sh" "${@:2}"; }
setsid bash -c 'GRM_CONNECTOR_DIR=$1 exec bash "$2/connector.sh" run' _ "$C2" "$S" >"$T/sup2.log" 2>&1 </dev/null &
check "the supervisor records its pid and start time" -- wait_for 10 test -s "$C2/supervisor.pid"
read -r sup _ <"$C2/supervisor.pid"
check "cloudflared starts" -- wait_for 10 has_child "$sup"
child=$(pgrep -P "$sup" -x sleep)
check "cloudflared does not hold connector.lock" -- \
    test -z "$(find "/proc/$child/fd" -lname "$C2/connector.lock" 2>/dev/null)"
expect 1 '' "the supervisor holds connector.lock" -- flock -n "$C2/connector.lock" true
expect 0 'already holds' "a second supervisor on one directory exits 0" -- connector "$C2" run
cp "$C2/supervisor.pid" "$C1/supervisor.pid"
expect 2 'not the supervisor of this connector directory: .*no signal sent' "stop refuses the supervisor of another directory" -- \
    connector "$C1" stop
check "the other supervisor still runs" -- kill -0 "$sup"
echo "$sup 1" >"$C1/supervisor.pid"
expect 0 'stale' "stop treats a pid with another start time as stale" -- connector "$C1" stop
check "the other supervisor still runs" -- kill -0 "$sup"
rm "$C1/supervisor.pid"
expect 0 'connector stopped' "stop signals its own supervisor" -- connector "$C2" stop
check "the supervisor and cloudflared exited" -- wait_for 5 test ! -d "/proc/$sup" -a ! -d "/proc/$child"

setsid bash -c 'GRM_CONNECTOR_DIR=$1 exec bash "$2/connector.sh" run' _ "$C2" "$S" >"$T/sup2.log" 2>&1 </dev/null &
check "the supervisor starts again" -- wait_for 10 test -s "$C2/supervisor.pid"
read -r sup _ <"$C2/supervisor.pid"
wait_for 10 has_child "$sup"
child=$(pgrep -P "$sup" -x sleep)
echo '#' >>"$C2/cloudflared"
kill -TERM "$child"
check "the supervisor checks the binary before restarting it and exits" -- wait_for 10 test ! -d "/proc/$sup"
check "the supervisor logged the changed binary" -- grep -q 'does not match its SHA-256' "$T/sup2.log"
(cd "$C2" && sha256sum cloudflared >cloudflared.sha256)

# Every fork and signal boundary of the supervisor: TERM (and INT, HUP) at each
# step exits 0 at once, leaving no supervisor record, cloudflared or lock, so
# no supervisor waits for a child it does not track.
conn_free() { ! locked "$1/connector.lock" && ! locked "$1/connector.spawn"; }
conn_run() { timeout -k 5 30 env GRM_CONNECTOR_DIR="$C2" bash "$S/connector.sh" run; }
for point in conn-locked conn-recorded conn-verified conn-forked conn-assigned; do
    at "$point" TERM timeout -k 5 30 env GRM_CONNECTOR_DIR="$C2" bash "$S/connector.sh" run
    check "connector TERM at $point: the supervisor exits 0 without waiting" -- test "$rc" = 0
    check "connector TERM at $point leaves no supervisor, cloudflared or lock" -- \
        eval 'conn_free "$C2" && test ! -e "$C2/supervisor.pid"'
done
for sig in INT HUP; do
    at conn-forked "$sig" timeout -k 5 30 env GRM_CONNECTOR_DIR="$C2" bash "$S/connector.sh" run
    check "connector $sig between the fork and the PID assignment exits 0 and leaves nothing" -- \
        eval 'test "$rc" = 0 && conn_free "$C2" && test ! -e "$C2/supervisor.pid"'
done

# A child paused before exec while its supervisor is killed outright never
# execs cloudflared: resumed, it sees the supervisor gone; or the next stop
# finds it by connector.spawn and stops it.
conn_bg() { preexec conn-child-preexec env GRM_CONNECTOR_DIR="$C2" bash "$S/connector.sh" run; }
conn_bg
check "the connector child paused before exec" -- test -n "$kid"
killbg
kill -CONT "$kid"
check "resumed after its supervisor was killed, the connector child exits without exec'ing" -- wait_for 10 gone "$kid"
expect 0 'stale' "connector stop treats the killed supervisor's record as stale" -- connector "$C2" stop
check "no cloudflared or lock is left" -- conn_free "$C2"
conn_bg
killbg
expect 0 'stopped the connector processes left by a killed supervisor' \
    "connector stop stops the paused child of a killed supervisor" -- connector "$C2" stop
check "the paused connector child is gone" -- wait_for 10 gone "$kid"
kill -CONT "$kid" 2>/dev/null || true
check "after the paused connector child, no cloudflared or lock is left" -- conn_free "$C2"

# A supervisor killed outright leaves cloudflared holding connector.spawn: the
# next stop, or the next run, stops it first.
sup_up() {
    rm -f "$C2/supervisor.pid"
    setsid bash -c 'GRM_CONNECTOR_DIR=$1 exec bash "$2/connector.sh" run' _ "$C2" "$S" >"$T/sup2.log" 2>&1 </dev/null &
    wait_for 10 test -s "$C2/supervisor.pid" || true
    read -r sup _ <"$C2/supervisor.pid"
    wait_for 10 has_child "$sup" || true
    child=$(pgrep -P "$sup" -x sleep || true)
}
sup_up
kill -KILL "$sup"
check "a killed supervisor's cloudflared keeps running, holding connector.spawn" -- \
    eval 'proc_alive "$child" && locked "$C2/connector.spawn"'
expect 0 'stopped the connector processes left by a killed supervisor' "connector stop stops a killed supervisor's cloudflared" -- \
    connector "$C2" stop
check "the orphaned cloudflared is gone" -- wait_for 10 gone "$child"
sup_up
old=$child
kill -KILL "$sup"
sup_up
check "the next run stopped the orphaned cloudflared first and started its own" -- \
    eval 'gone "$old" && test -n "$child" && proc_alive "$child"'
check "the next run reported the orphan it stopped" -- grep -q 'stopped the connector processes left by a killed supervisor' "$T/sup2.log"
expect 0 'connector stopped' "connector stop stops that supervisor" -- connector "$C2" stop
check "it left no cloudflared or lock" -- wait_for 10 conn_free "$C2"
# The same with the next run's output pipe closed at once: its report of the
# orphan fails to write instead of killing it before it stops the orphan.
sup_up
old=$child
kill -KILL "$sup"
rm -f "$C2/supervisor.pid"
setsid bash -c 'GRM_CONNECTOR_DIR=$1 bash "$2/connector.sh" run 2>&1 | true' _ "$C2" "$S" </dev/null &
wait_for 10 test -s "$C2/supervisor.pid" || true
sup=$(sed -n 's/ .*//p' "$C2/supervisor.pid" 2>/dev/null || true)
check "with its output pipe closed, the next run stopped the orphaned cloudflared" -- wait_for 10 gone "$old"
check "with its output pipe closed, the next run started its own cloudflared" -- eval 'test -n "$sup" && wait_for 10 has_child "$sup"'
expect 0 'connector stopped' "connector stop stops the supervisor whose output pipe is closed" -- connector "$C2" stop
check "that supervisor left no cloudflared or lock" -- wait_for 10 conn_free "$C2"

# A holder of connector.spawn this user cannot see: run starts nothing.
hidden_lock "$C2/connector.spawn"
expect 1 'may remain' "connector run starts nothing while a hidden process holds connector.spawn" -- conn_run
check "no supervisor record or cloudflared was started" -- eval 'test ! -e "$C2/supervisor.pid" && ! locked "$C2/connector.lock"'
kill -KILL "$hidden"
wait "$hidden" 2>/dev/null || true

# Links, FIFOs and special files at every connector path are refused before
# anything is opened, so nothing blocks.
for f in connector.lock connector.spawn supervisor.pid connector.pid logs/cloudflared.log; do
    rm -f "$C2/$f"
    mkfifo "$C2/$f"
    expect 2 'link, directory, FIFO or other special file' "connector run refuses $f as a FIFO without blocking" -- conn_run
    rm "$C2/$f"
done
ln -s "$T/outside" "$C2/supervisor.pid"
expect 2 'link, directory, FIFO or other special file' "connector run refuses supervisor.pid as a link" -- conn_run
rm "$C2/supervisor.pid"
mv "$C2/logs" "$T/conn-logs"
ln -s "$T/conn-logs" "$C2/logs"
expect 2 'link, directory, FIFO or other special file' "connector run refuses logs/ as a link" -- conn_run
rm "$C2/logs"
mv "$T/conn-logs" "$C2/logs"
check "none of the refused runs started cloudflared" -- conn_free "$C2"

# Identity is the open file description holding connector.lock's flock, with
# the recorded start time and working directory; not a path, a command line
# or an open descriptor 9. A forged holder waits for the lock (flock -w), and
# the test sees it hold the lock in its fdinfo instead of probing with flock,
# which could take the lock at the forger's instant and make it give up.
fd_locked() { grep -q '^lock:.*FLOCK.*WRITE' "/proc/$1/fdinfo/$2" 2>/dev/null; }
(cd "$C2" && exec 9<>connector.lock && exec bash -c 'sleep 300; :' "$S/connector.sh" run) >/dev/null 2>&1 </dev/null &
forged=$!
wait_for 5 test -e "/proc/$forged/fd/9" || true
fstart=$(start_of "$forged")
printf '%s %s\n' "$forged" "$fstart" >"$C2/supervisor.pid"
expect 2 'no process holds the lock on .*connector.lock; no signal sent' \
    "connector stop refuses a shell in the directory with fd 9 on connector.lock while the lock is free" -- connector "$C2" stop
check "the forged supervisor got no signal" -- proc_alive "$forged"
(cd "$T" && exec 9<>"$C2/connector.lock" && flock -w 10 9 && exec bash -c 'sleep 300; :' "$S/connector.sh" run) >/dev/null 2>&1 </dev/null &
forged3=$!
wait_for 10 fd_locked "$forged3" 9 || true
printf '%s %s\n' "$forged3" "$(start_of "$forged3")" >"$C2/supervisor.pid"
expect 2 'works in another directory; no signal sent' \
    "connector stop refuses the holder of connector.lock when it works in another directory" -- connector "$C2" stop
check "that lock holder got no signal" -- proc_alive "$forged3"
pkill -KILL -P "$forged3" || true
kill -KILL "$forged3" 2>/dev/null || true
wait "$forged3" 2>/dev/null || true
wait_for 5 eval '! locked "$C2/connector.lock"' || true
sup_up
cp "$C2/supervisor.pid" "$T/supervisor.real"
printf '%s %s\n' "$forged" "$fstart" >"$C2/supervisor.pid"
expect 2 'does not hold the lock on .*connector.lock; no signal sent' \
    "connector stop refuses that shell while the real supervisor holds the lock" -- connector "$C2" stop
(cd "$C2" && exec 9<>"$C1/connector.lock" && flock -w 10 9 && exec bash -c 'sleep 300; :' "$S/connector.sh" run) >/dev/null 2>&1 </dev/null &
forged2=$!
wait_for 10 fd_locked "$forged2" 9 || true
printf '%s %s\n' "$forged2" "$(start_of "$forged2")" >"$C2/supervisor.pid"
expect 2 'does not hold the lock on .*connector.lock; no signal sent' \
    "connector stop refuses a process in the directory holding another directory's lock" -- connector "$C2" stop
check "the real and forged supervisors got no signal" -- eval 'proc_alive "$sup" && proc_alive "$forged" && proc_alive "$forged2"'
read -r rsup _ <"$T/supervisor.real"
printf '%s 1\n' "$rsup" >"$C2/supervisor.pid"
expect 0 'stale' "connector stop treats the real supervisor's PID with another start time as stale" -- connector "$C2" stop
check "the real supervisor got no signal" -- proc_alive "$sup"
cp "$T/supervisor.real" "$C2/supervisor.pid"
expect 0 'connector stopped' "connector stop signals the real lock holder" -- connector "$C2" stop
check "the real supervisor and its cloudflared exited" -- eval 'wait_for 10 gone "$sup" && gone "$child" && conn_free "$C2"'
pkill -KILL -P "$forged" || true
pkill -KILL -P "$forged2" || true
kill -KILL "$forged" "$forged2" 2>/dev/null || true

# A cloudflared that leaves a process in its group. With kid.mode in the
# connector directory, the stub forks a descendant that closes every inherited
# descriptor (connector.spawn included) and stays in its process group, waits
# until that descendant runs sleep (so its environment and signal handling are
# final), writes kid.pid and stub.pid, and exits. "keep" keeps the run's token;
# "ignore" keeps it and ignores SIGTERM; "untied" runs the stub again without
# the token first, so the descendant never carries it.
C3=$T/conn3
make_connector "$C3"
cat >"$C3/cloudflared" <<'SH'
#!/bin/sh
case "$1" in --version) echo "cloudflared version test"; exit 0 ;; esac
mode=$(cat kid.mode 2>/dev/null || true)
[ -n "$mode" ] || exec sleep 300
if [ "$mode" = untied ] && [ -n "${GRM_DEPLOY_SPAWN:-}" ]; then exec env -u GRM_DEPLOY_SPAWN "$0" "$@"; fi
rm -f kid.mode
python3 -I -c '
import os, signal, sys
if sys.argv[1] == "ignore":
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
os.closerange(3, 65536)
os.execvp("sleep", ["sleep", "300"])
' "$mode" </dev/null >/dev/null 2>&1 &
kid=$!
for _ in $(seq 100); do [ "$(cat "/proc/$kid/comm" 2>/dev/null)" = sleep ] && break; sleep 0.05; done
echo "$kid" >kid.pid.tmp && mv kid.pid.tmp kid.pid
echo "$$" >stub.pid
exit 0
SH
(cd "$C3" && sha256sum cloudflared >cloudflared.sha256)
sup3_up() {  # start a supervisor on C3 in its own session (sup=its PID), its exit code to sup.rc
    rm -f "$C3/supervisor.pid" "$C3/sup.rc" "$C3/stub.pid" "$C3/kid.pid"
    setsid bash -c 'GRM_CONNECTOR_DIR=$1 bash "$2/connector.sh" run; echo "$?" >"$1/sup.rc"' _ "$C3" "$S" >"$T/sup3.log" 2>&1 </dev/null &
    wait_for 10 test -s "$C3/supervisor.pid" || true
    sup=$(sed -n 's/ .*//p' "$C3/supervisor.pid" 2>/dev/null || true)
}
stub_done() {  # wait for the stub to write its descendant's PID and exit (kid=that PID, stub_ok=1)
    kid='' stub_ok=0
    wait_for 15 test -s "$C3/stub.pid" || return 0
    kid=$(cat "$C3/kid.pid")
    if wait_for 5 gone "$(cat "$C3/stub.pid")"; then stub_ok=1; fi
}
conn3_run() { timeout -k 5 30 env GRM_CONNECTOR_DIR="$C3" bash "$S/connector.sh" run; }
echo keep >"$C3/kid.mode"
sup3_up
stub_done
check "a cloudflared that leaves a descendant in its group starts and exits" -- test "$stub_ok" = 1
check "the supervisor then stops that descendant, which carries the run's token" -- eval 'test -n "$kid" && wait_for 15 gone "$kid"'
expect 0 'a connector supervisor holds' "connector status reports the running supervisor" -- connector "$C3" status
# Ended by a signal it does not handle (USR1, which still runs its EXIT trap)
# while it waits for a descendant that ignores SIGTERM, the supervisor stops
# that descendant on its way out, though cloudflared itself has exited.
# The mode is set only once the supervisor has restarted a plain cloudflared,
# so the descendant comes from the start that follows this TERM.
rm -f "$C3/stub.pid" "$C3/kid.pid"
wait_for 15 eval 'pgrep -P "$sup" -f "^sleep 300$" >/dev/null' || true
echo ignore >"$C3/kid.mode"
pkill -TERM -P "$sup" -f '^sleep 300$' || true
stub_done
check "the restarted cloudflared leaves a descendant that ignores SIGTERM and exits" -- test "$stub_ok" = 1
sleep 1
kill -USR1 "$sup" 2>/dev/null || true
check "a supervisor ended by USR1 stops that descendant on its way out" -- \
    eval 'test -n "$kid" && wait_for 30 gone "$kid" && wait_for 10 gone "$sup"'
[ -z "$kid" ] || kill -KILL "$kid" 2>/dev/null || true
# A descendant without the run's token, in the group of a cloudflared that has
# exited, cannot be told apart from an unrelated group that took the number:
# the supervisor signals nothing, records the failed stop and exits 3. Nothing
# restarts it; status reports it, and run refuses until it is resolved by hand.
echo untied >"$C3/kid.mode"
sup3_up
stub_done
check "a cloudflared that leaves a descendant without the run's token starts and exits" -- test "$stub_ok" = 1
check "the supervisor exits 3 when it cannot tie that descendant to its run" -- \
    eval 'wait_for 15 test -s "$C3/sup.rc" && test "$(cat "$C3/sup.rc")" = 3'
check "it did not signal that descendant" -- eval 'test -n "$kid" && proc_alive "$kid"'
check "its log names the process and says the tunnel is down with nothing to restart it" -- \
    eval 'grep -q "(pids $kid) that cannot be tied to this start" "$T/sup3.log" && grep -q "nothing restarts it, so the tunnel is down" "$T/sup3.log"'
expect 3 'records a failed stop: .*cannot be tied to this start' "connector status reports the recorded failed stop (exit 3)" -- connector "$C3" status
expect 3 'not starting: .*records a failed stop' "connector run refuses to start while connector.failed exists" -- conn3_run
check "that refused run started no supervisor or cloudflared" -- eval 'conn_free "$C3" && test ! -e "$C3/supervisor.pid"'
[ -z "$kid" ] || kill -KILL "$kid" 2>/dev/null || true
wait_for 5 gone "$kid" || true
rm -f "$C3/connector.failed"
expect 1 'no connector supervisor holds' "after the manual recovery, status reports no supervisor and no failure" -- connector "$C3" status
# A cloudflared that exits before the supervisor has read its start time: the
# supervisor stops itself right after the fork (conn-forked) until the stub has
# exited. What the stub left in its group carries the run's token and must
# still be stopped. The supervisor stops itself again at its next fork, so TERM
# and CONT end it there.
rm -f "$C3/supervisor.pid" "$C3/stub.pid" "$C3/kid.pid"
echo keep >"$C3/kid.mode"
paused conn-forked env GRM_CONNECTOR_DIR="$C3" bash "$S/connector.sh" run
stub_done
check "with its supervisor paused right after the fork, cloudflared leaves a descendant and exits" -- test "$stub_ok" = 1
kill -CONT "$bg" 2>/dev/null || true
check "the supervisor still stops that descendant, though it never read cloudflared's start time" -- \
    eval 'test -n "$kid" && wait_for 15 gone "$kid"'
wait_for 15 stopped "$bg" || true
kill -TERM "$bg" 2>/dev/null || true
kill -CONT "$bg" 2>/dev/null || true
check "TERM then ends that supervisor and its next cloudflared" -- eval 'wait_for 15 gone "$bg" && wait_for 10 conn_free "$C3"'
[ -z "$kid" ] || kill -KILL "$kid" 2>/dev/null || true
wait "$bg" 2>/dev/null || true

EXPECTED_ASSERTIONS=$((488 - skipped))  # less the user and PID namespace cases where those are unavailable
if [ "$(id -u)" = 0 ]; then EXPECTED_ASSERTIONS=$((EXPECTED_ASSERTIONS - 2)); fi  # the unwritable state/ case is skipped
if [ "$((pass + fail))" != "$EXPECTED_ASSERTIONS" ]; then
    fail=$((fail + 1))
    echo "FAIL apply suite ran $((pass + fail - 1)) assertions, expected $EXPECTED_ASSERTIONS"
fi
finished=1
echo "apply tests: $pass passed, $fail failed"
if ((fail)); then exit 1; fi
