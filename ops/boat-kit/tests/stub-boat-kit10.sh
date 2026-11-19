#!/usr/bin/env bash
# KIT10's stand-in for boat.exe (GRM_BOAT=this), from KIT9's stub-boat-kit9.sh: info answers a ready grm-kit10-stub
# with archiveAfter from STUB_ROOT/archive, exec runs its bash -c script locally under STUB_ROOT and answers in boat's
# JSON (and --detach / --status for run), scp copies between the cwd and STUB_ROOT, extend sets STUB_ROOT/archive,
# limits and list answer for preflight. Nothing remote. Knobs:
#   STUB_SCP_RATE=B/s       copy at that rate (chunks once a second) instead of at once
#   STUB_SCP_SLOW_TRIES=N   only the first N scp calls are slow, delayed or hung (STUB_ROOT/scp.count)
#   STUB_SCP_HANG=1         write nothing and sleep (a silent scp)
#   STUB_SCP_DELAY=s        wait s seconds, writing nothing, before reading the source (a slow start-up)
#   STUB_STAT_FAIL=1        an exec that reads the copy's size (stat -c %s) fails without JSON
#   STUB_ARCHIVE_BACK=1     extend moves archiveAfter back an hour (a misbehaving API)
#   STUB_ARCHIVE_SAME=1     extend returns ok and leaves archiveAfter as it was (REVIEW-KIT10 L2)
#   STUB_ARCHIVE_SHORT=s    extend returns ok and sets archiveAfter s seconds short of the asked end (L2)
#   STUB_SCP_TAMPER=1       KIT11: scp adds a byte to what it copied (a file changed on the way)
c=$1; shift; [ "${1:-}" = --no-update ] && shift
R=${STUB_ROOT:?}
echo "$c $*" >> "$R/boat-calls.log"
case "$c" in
  info)
    a=$(cat "$R/archive" 2> /dev/null || echo null); [ "$a" = null ] || a="\"$a\""
    echo "{\"sandbox\": {\"state\": \"ready\", \"name\": \"grm-kit10-stub\", \"snapshots\": false, \"archiveAfter\": $a}}";;
  extend)
    t=""; while [ $# -gt 0 ]; do [ "$1" = --ttl ] && t=$2; shift; done
    if [ -n "${STUB_ARCHIVE_SAME:-}" ]; then :
    elif [ -n "${STUB_ARCHIVE_BACK:-}" ]; then date -u -d '-1 hour' +%FT%T.000Z > "$R/archive"
    else date -u -d "+$(( t - ${STUB_ARCHIVE_SHORT:-0} )) sec" +%FT%T.000Z > "$R/archive"; fi
    echo '{"ok": true}';;
  limits) echo '{"creditBalanceHours": 1000.5, "planIncluded": {"vmHoursRemaining": 900}}';;
  list) echo '{"sandboxes": []}';;
  scp)
    s=$1 d=$2; s=${s/bx_stub1:/$R}; d=${d/bx_stub1:/$R}; mkdir -p "$(dirname "$d")"
    n=$(( $(cat "$R/scp.count" 2> /dev/null || echo 0) + 1 )); echo $n > "$R/scp.count"
    slow=1; [ -n "${STUB_SCP_SLOW_TRIES:-}" ] && [ "$n" -gt "$STUB_SCP_SLOW_TRIES" ] && slow=0
    if [ $slow = 1 ] && [ -n "${STUB_SCP_HANG:-}" ]; then : > "$d"; exec sleep 100000; fi
    [ $slow = 1 ] && [ -n "${STUB_SCP_DELAY:-}" ] && sleep "$STUB_SCP_DELAY"
    if [ $slow = 1 ] && [ -n "${STUB_SCP_RATE:-}" ]; then
      exec python3 -I -c 'import sys, time
s, d, r = sys.argv[1], sys.argv[2], int(sys.argv[3])
with open(s, "rb") as f, open(d, "wb") as o:
    while True:
        b = f.read(r)
        if not b: break
        o.write(b); o.flush(); time.sleep(1)' "$s" "$d" "$STUB_SCP_RATE"
    fi
    cp "$s" "$d" || exit 1; [ -n "${STUB_SCP_TAMPER:-}" ] && echo x >> "$d"; exit 0;;
  exec)
    det=0; st=""
    while [ $# -gt 0 ] && [ "$1" != -- ]; do case "$1" in --detach) det=1;; --status) st=$2; shift;; esac; shift; done
    [ $# -gt 0 ] && shift
    if [ -n "$st" ]; then
      if kill -0 "$st" 2> /dev/null; then echo '{"running": true}'; else echo '{"running": false, "exitCode": null}'; fi; exit 0
    fi
    if [ $det = 1 ]; then   # bash /home/user/grm-logs/NAME.sh: the launcher, already written under STUB_ROOT
      f=${2/\/home\/user/$R/home/user}; setsid bash "$f" < /dev/null > /dev/null 2>&1 & echo "{\"pid\": $!}"; exit 0
    fi
    script=$(printf '%s' "$3" | sed 's/^echo \([A-Za-z0-9+/=]*\) | base64 -d | bash$/\1/' | base64 -d \
      | sed "s#/home/user#$R/home/user#g; s#/tmp/grm-#$R/tmp/grm-#g")
    if [ -n "${STUB_STAT_FAIL:-}" ] && grep -q 'stat -c %s' <<< "$script"; then echo "stub: exec failed (vsock)" >&2; exit 1; fi
    out=$(mkdir -p "$R/tmp" "$R/home/user"; bash -c "$script" 2>&1); rc=$?
    python3 -I -c 'import json, sys; print(json.dumps({"stdout": sys.argv[1], "stderr": "", "exitCode": int(sys.argv[2])}))' "$out" "$rc";;
  *) echo "stub: $c not handled" >&2; exit 1;;
esac
