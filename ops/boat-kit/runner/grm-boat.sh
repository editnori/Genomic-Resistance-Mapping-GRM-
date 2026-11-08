#!/usr/bin/env bash
# grm-boat.sh: run one GRM lane's Rust gates in its own Boat sandbox from WSL, with the Boat CLI problems the GRM lanes
# hit handled once (GRM-BOAT-LANES-20261031.md). Every step appends one line to the caller's NOTES.md:
# "- UTC | sandbox id | step | result | detail (sha256s)".
#   grm-boat.sh [--notes PATH] new <lane> [--ttl s] [--type T] [--from-template NAME] [--template]
#       laptop memory check (decision 433) and plan floor, cap check against GRM_CAP, then create grm-<lane>, xlarge
#       unless --type (--no-env --no-snapshots --no-monitoring); prints the id.
#       --from-template starts it from a saved grm-tmpl-* snapshot; --template (lane tmpl-*) leaves out --no-snapshots.
#   grm-boat.sh [--notes PATH] toolchain <id>
#       Rust 1.98.1 (minimal) as the default under CARGO_HOME=/usr/local/cargo; fails unless rustc prints 1.98.1
#   grm-boat.sh [--notes PATH] push <id> <worktree> [--with-tmp] [--code-only] [--dest DIR]
#       git ls-files -co --exclude-standard minus target/ (and .tmp/ unless --with-tmp), through the data guard, tarred
#       in a Windows staging dir, boat scp, sha256 on both sides, extracted into DIR (default /home/user/grm).
#       --code-only: committed files only (clean worktree), minus every *.gz and every data/ path (template fill).
#   grm-boat.sh [--notes PATH] push-data <id> <dir> --manifest FILE --manifest-sha256 HEX [--dest DIR]
#       public GRM stage data (PF-D2). FILE is copied once and must hash to HEX; each line "<sha256>  <path>" must name
#       an allowlisted class (BV-BRC FASTA and metadata, annotation cache, curated, catalog, matrices, report trees,
#       run receipts: text under audit/, receipts/ or tables/ and cohorts/<c>/validation-manifest.json, at most 64 MiB)
#       with no denied segment, in <dir> (whose own path is denied-name free too). The tar is built from the list with
#       no symlink followed, then read back: exactly those regular members, each hashing to its line, are what is sent.
#       Refusal: 97, nothing sent. In the box: tar sha256, extract into DIR (default /home/user/grm-data, new or
#       empty), sha256sum -c --strict and the count; any failure there removes DIR and the tar (86). The manifest
#       stays as .grm-input-manifest.sha256. push's data guard is unchanged.
#   grm-boat.sh [--notes PATH] run <id> [--log NAME] [--cwd DIR] [--argv] [--] <command...>
#       one argument is a legacy bash script unless --argv is set; multiple arguments are executed directly with their
#       boundaries preserved. The command runs detached in DIR with the toolchain env and ulimit -n 65536; output goes to
#       /home/user/grm-logs/NAME.log (rc in NAME.rc, /usr/bin/time -v in NAME.time); polls --status until it exits,
#       prints the log tail and returns the command's exit code. DIR must exist (push creates the default); a cd
#       failure logs the directory and bootstrap action, returns 86 and does not run the command.
#   grm-boat.sh [--notes PATH] pull <id> <remote-path> <local-dir>
#       tar on the sandbox, boat scp back, sha256 on both sides, extract into local-dir
#   grm-boat.sh [--notes PATH] snapshot <id> <grm-tmpl-name>      save a grm-tmpl-* sandbox as a named template
#   grm-boat.sh close <id> <notes.md>
#       boat usage appended to the notes, boat delete --yes, then boat info must answer 404
#   grm-boat.sh preflight      the laptop memory check, the plan floor and the live grm- count; creates nothing
# NOTES.md comes from --notes or GRM_NOTES (close takes it as its argument). Only ids this script created are used
# (state in ~/.local/state/grm-boat/created; GRM_BOAT_ADOPT=1 accepts a grm-* sandbox made by hand).
# Exit codes: run returns the command's own; 97 the data guard or push-data refused (nothing sent); 86 a Boat step or working-directory
# change failed;
# 75 the cap stayed full for 30 min; 76 the laptop memory check failed; 77 the plan is under its floor; 2 usage or
# refused id.
set -u -o pipefail
BOAT=${GRM_BOAT:-/mnt/c/Users/layth888/.ascii/bin/boat.exe}
LANES_MD=${GRM_LANES_MD:-/home/layth/review/GRM-BOAT-LANES-20261031.md}   # holds the live GRM_CAP= line
STAGE_WIN=${GRM_STAGE_WIN:-'C:\Users\layth888\.t3\scratch\2026-10-28-test-3a358d9a\boat-grm'}
STATE=${GRM_BOAT_STATE:-$HOME/.local/state/grm-boat}
LOCK=${GRM_BOAT_LOCK:-/var/tmp/grm-boat-create.lock}   # one cap check and create at a time on this laptop
BT=${GRM_BT_TIMEOUT:-180}       # seconds for one boat call; a WSL interop call that stalls is killed after this
RDIR=/home/user/grm RLOG=/home/user/grm-logs
TC=1.98.1 TYPE=xlarge
ROUTE_RC=86 GUARD_RC=97 CAP_RC=75 MEM_RC=76 PLAN_RC=77
PWSH=${GRM_POWERSHELL:-/mnt/c/Windows/System32/WindowsPowerShell/v1.0/powershell.exe}
MIN_FREE_GB=${GRM_MIN_FREE_GB:-4} MAX_COMMIT_GB=${GRM_MAX_COMMIT_GB:-95}   # decision 433's laptop check
MIN_PLAN_VMH=${GRM_MIN_PLAN_VMH:-150}   # user rule, 23:4xZ Nov 1: under this many plan VM-hours, start no new box;
# since Layth's prepaid pack (6 Nov) the floor applies to plan plus pack (boat limits: creditBalanceHours)
NOTES=${GRM_NOTES:-}
# every remote step starts with this: the image's cargo home, the pinned toolchain's PATH and the open-file limit
RENV='export CARGO_HOME=/usr/local/cargo PATH=/usr/local/cargo/bin:$PATH
[ -d /usr/local/rustup ] && export RUSTUP_HOME=/usr/local/rustup
ulimit -n 65536'

now() { date -u +%FT%TZ; }
usage() { sed -n '2,/^set -u/p' "$0" | sed '$d' | sed 's/^# \{0,1\}//' >&2; exit 2; }
note() { [ -z "$NOTES" ] || printf -- '- %s | %s | %s | %s | %s\n' "$(now)" "$1" "$2" "$3" "${4:-}" >> "$NOTES"; }
# fail RC ID STEP MESSAGE: record the failure in the notes, say it on stderr, exit RC
fail() { local rc=$1 id=$2 step=$3; shift 3; note "$id" "$step" "FAIL rc=$rc" "$*"; echo "grm-boat: $step: $*" >&2; exit "$rc"; }
# one boat call through WSL interop: no stdin, no self-update, killed (rc 124) when it stalls
bt() { local c=$1; shift; timeout -k 10 "$BT" "$BOAT" "$c" --no-update "$@" < /dev/null; }
# bt_retry: the same, three tries, for calls that change nothing
bt_retry() { local t out; for t in 1 2 3; do out=$(bt "$@" 2>&1) && { printf '%s\n' "$out"; return 0; }; sleep 5; done; printf '%s\n' "$out"; return 1; }
oneline() { printf %s "$1" | tr -s '\r\n\t' '   ' | cut -c1-"${2:-240}"; }
sha() { sha256sum "$1" | cut -d' ' -f1; }

# binfo ID -> "state name snapshots archiveAfter"; rc 3 when Boat no longer knows the sandbox (404), and on any
# failure the error text instead (also left in INFO_ERR)
binfo() {
  local out rc
  out=$(bt info "$1" --json 2>&1); rc=$?
  if [ $rc != 0 ]; then
    INFO_ERR=$(oneline "$out" 160); printf '%s\n' "$INFO_ERR"
    case "$out" in *[Nn]ot*[Ff]ound*|*404*|*[Nn]o\ such*) return 3;; esac; return 1
  fi
  printf '%s' "$out" | python3 -I -c 'import json, sys
d = json.load(sys.stdin); d = d.get("sandbox", d)
print(*["-" if d.get(k) is None else str(d.get(k)).lower() for k in ("state", "name", "snapshots", "archiveAfter")])'
}
# list_live: live sandboxes (up, pending, stopping) as "id name snapshots createdAt" lines; rc 1 when Boat can't be read
list_live() {
  bt_retry list --filter rpt --json | python3 -I -c 'import json, sys
d = json.load(sys.stdin)
for s in d.get("sandboxes", []):
    print(s.get("id"), s.get("name") or "-", str(s.get("snapshots")).lower(), s.get("createdAt") or "-")
p = d.get("pageInfo") or {}
if p.get("hasMore") or p.get("hasNextPage"): print("MORE - - -")'
}
cap_now() { sed -n 's/^GRM_CAP=\([0-9][0-9]*\).*/\1/p' "$LANES_MD" | head -n 1; }
# cap_wait: GRM_CAP is re-read from the lanes file at every check; wait in 60 s steps for up to 30 min
cap_wait() {
  local i cap n ls
  for i in $(seq 0 30); do
    cap=$(cap_now); [ -n "$cap" ] || fail 2 - cap "no GRM_CAP= line in $LANES_MD"
    n=""; if ls=$(list_live); then
      grep -q '^MORE ' <<< "$ls" && fail $ROUTE_RC - cap "boat list has more than one page; count by hand"
      n=$(awk '$2 ~ /^grm-/' <<< "$ls" | wc -l)
    fi
    if [ -n "$n" ] && [ "$n" -lt "$cap" ]; then CAPMSG="cap: $n live grm- < GRM_CAP $cap"; return 0; fi
    [ "$i" = 30 ] && break
    echo "grm-boat: cap check: ${n:-unreadable} live grm- sandboxes, GRM_CAP $cap; waiting 60 s ($((i + 1))/30)" >&2
    sleep 60
  done
  fail $CAP_RC - cap "${n:-unreadable} live grm- sandboxes against GRM_CAP $cap for 30 min; none freed"
}
# host_guard: before every new sandbox, decision 433's laptop check (Windows free physical memory and commit charge)
# and the user rule's floor on VM-hours left (boat limits --all: creditBalanceHours, the plan's remainder plus the
# prepaid pack; planIncluded.vmHoursRemaining if the account reports no credit balance); sets HOSTMSG
host_guard() {
  local m free commit lim vmh
  m=$(timeout -k 5 60 "$PWSH" -NoProfile -NonInteractive -Command '$o = Get-CimInstance Win32_OperatingSystem; "{0:F2} {1:F2}" -f ($o.FreePhysicalMemory / 1MB), (($o.TotalVirtualMemorySize - $o.FreeVirtualMemory) / 1MB)' < /dev/null 2> /dev/null | tr -d '\r')
  read -r free commit <<< "$m"
  [[ "$free" =~ ^[0-9]+(\.[0-9]+)?$ && "$commit" =~ ^[0-9]+(\.[0-9]+)?$ ]] \
    || fail $MEM_RC - memcheck "could not read Windows memory ($(oneline "$m" 80))"
  awk -v f="$free" -v c="$commit" -v mf="$MIN_FREE_GB" -v mc="$MAX_COMMIT_GB" 'BEGIN { exit !(f >= mf && c <= mc) }' \
    || fail $MEM_RC - memcheck "Windows free $free GB (needs $MIN_FREE_GB or more), commit $commit GB (needs $MAX_COMMIT_GB or less)"
  lim=$(bt_retry limits --all --json) || fail $ROUTE_RC - plan "boat limits failed: $(oneline "$lim")"
  vmh=$(python3 -c 'import json, sys
l = json.load(sys.stdin)
c = l.get("creditBalanceHours")
print(c if isinstance(c, (int, float)) and not isinstance(c, bool) else l["planIncluded"]["vmHoursRemaining"])' <<< "$lim" 2> /dev/null)
  [[ "$vmh" =~ ^[0-9]+(\.[0-9]+)?$ ]] || fail $ROUTE_RC - plan "no creditBalanceHours or planIncluded.vmHoursRemaining in boat limits"
  awk -v v="$vmh" -v m="$MIN_PLAN_VMH" 'BEGIN { exit !(v >= m) }' \
    || fail $PLAN_RC - plan "the plan and pack have $vmh VM-hours left, under $MIN_PLAN_VMH: start no new box and tell the user"
  HOSTMSG="host: free $free GB, commit $commit GB; plan $vmh VM-h left"
}
# own ID: refuse any sandbox this script did not create or whose name is not grm-*; sets ID NAME LANE SNAP
own() {
  local info st
  ID=$1
  [[ "$ID" =~ ^bx_[a-z0-9]+$ ]] || fail 2 "$ID" check "not a sandbox id"
  grep -q "^$ID " "$STATE/created" 2> /dev/null || [ "${GRM_BOAT_ADOPT:-}" = 1 ] ||
    fail 2 "$ID" check "this script did not create $ID (GRM_BOAT_ADOPT=1 to use a grm- sandbox you made by hand)"
  info=$(binfo "$ID") || fail $ROUTE_RC "$ID" check "boat info failed: $info"
  read -r st NAME SNAP _ <<< "$info"
  case "$NAME" in grm-*) ;; *) fail 2 "$ID" check "$ID is named $NAME, not grm-*: another lane's sandbox";; esac
  LANE=${NAME#grm-}
}
# rexec ID SECS SCRIPT: run a bash script on the sandbox and wait. The script crosses WSL interop and the Boat API as
# base64 inside one argument, so no quoting is lost on the way. Off a TTY boat exec answers in JSON (stdout, stderr,
# exitCode): print stdout then stderr and return exitCode. SECS stays at 600 or less: above that the output streams.
rexec() {
  local b64 out rc; b64=$(printf '%s\n' "$3" | base64 -w0)
  out=$(BT=$(($2 + 60)) bt exec "$1" --json --timeout "$2" -- bash -c "echo $b64 | base64 -d | bash" 2>&1); rc=$?
  printf '%s' "$out" | python3 -I -c 'import json, sys
t = sys.stdin.read()
try:
    d = json.loads(t)
except ValueError:
    sys.stdout.write(t); sys.exit(int(sys.argv[1]) or 1)
sys.stdout.write((d.get("stdout") or "") + (d.get("stderr") or ""))
if d.get("stdoutTruncated") or d.get("stderrTruncated"): print("[boat exec output truncated]")
if d.get("timedOut"): print("[boat exec timed out]")
e = d.get("exitCode"); sys.exit(e if isinstance(e, int) else 1)' "$rc"
}
stage_dir() { STAGE=$(wslpath -u "$STAGE_WIN")/stage-$LANE; mkdir -p "$STAGE"; }

# ---------------------------------------------------------------- new
cmd_new() {
  local lane="" ttl=7200 from="" tmpl=0 name out id t0 info st snap i args
  while [ $# -gt 0 ]; do case "$1" in
    --ttl) ttl=${2:-}; shift 2 || usage;;
    --from-template) from=${2:-}; shift 2 || usage;;
    --type) TYPE=${2:-}; shift 2 || usage;;
    --template) tmpl=1; shift;;
    -*) usage;;
    *) [ -z "$lane" ] || usage; lane=$1; shift;;
  esac; done
  [[ "$lane" =~ ^[a-z0-9][a-z0-9-]{0,100}$ ]] || fail 2 - new "lane must be lowercase letters, digits and dashes"
  [[ "$ttl" =~ ^[0-9]+$ ]] || fail 2 - new "--ttl takes seconds"
  [[ "$TYPE" =~ ^[a-z0-9-]+$ ]] || fail 2 - new "--type takes a Boat sandbox type, such as default or xlarge"
  case "$from" in ""|grm-tmpl-*) ;; *) fail 2 - new "--from-template takes a grm-tmpl-* snapshot name";; esac
  if [ $tmpl = 1 ]; then
    case "$lane" in tmpl-*) ;; *) fail 2 - new "--template needs a lane named tmpl-<what>";; esac
    [ -z "$from" ] || fail 2 - new "--template and --from-template don't mix"
  fi
  name=grm-$lane
  mkdir -p "$STATE"
  exec 9> "$LOCK"; flock 9
  host_guard
  cap_wait
  if out=$(list_live) && awk -v n="$name" '$2 == n' <<< "$out" | grep -q .; then fail 2 - new "a live sandbox is already named $name"; fi
  args=(--no-env --no-snapshots --no-monitoring --type "$TYPE" --name "$name" --ttl "$ttl" --json)
  [ $tmpl = 1 ] && args=(--no-env --no-monitoring --type "$TYPE" --name "$name" --ttl "$ttl" --json)
  [ -n "$from" ] && args=(--from "$from" "${args[@]}")
  t0=$(date -u -d '-10 sec' +%FT%T.000Z)
  out=$(BT=300 bt new "${args[@]}" 2>&1)
  id=$(printf '%s\n' "$out" | grep -o '"id": *"bx_[a-z0-9]*"' | head -n 1 | grep -o 'bx_[a-z0-9]*')
  # a create that timed out can still have made the sandbox (bheavy.sh saw one hang 600 s): adopt it by its name
  [ -n "$id" ] || id=$(list_live | awk -v n="$name" -v t="$t0" '$2 == n && $4 >= t {print $1; exit}')
  [ -n "$id" ] || fail $ROUTE_RC - "new $name" "boat new gave no id: $(oneline "$out")"
  echo "$id $name $(now)" >> "$STATE/created"
  flock -u 9; exec 9>&-
  st=-; for i in $(seq 1 50); do
    if info=$(binfo "$id"); then read -r st _ _ _ <<< "$info"; case "$st" in ready|idle|running) break;; esac; fi
    sleep 6
  done
  case "$st" in ready|idle|running) ;; *) fail $ROUTE_RC "$id" "new $name" "state $st after 5 min; close it";; esac
  snap=$(list_live | awk -v i="$id" '$1 == i {print $3}')
  if [ $tmpl = 0 ] && [ "$snap" != false ]; then fail $ROUTE_RC "$id" "new $name" "boat list shows snapshots=${snap:-?}, not false; close it"; fi
  note "$id" "new $name" ok "type=$TYPE ttl=$ttl state=$st snapshots=$snap${from:+ from=$from}; $CAPMSG; $HOSTMSG"
  echo "$id"
}

# ---------------------------------------------------------------- preflight
cmd_preflight() {
  local n
  [ $# = 0 ] || usage
  host_guard
  n=$(list_live | awk '$2 ~ /^grm-/' | wc -l)
  echo "grm-boat: preflight ok; $HOSTMSG; $n live grm- sandboxes, GRM_CAP $(cap_now)"
}

# ---------------------------------------------------------------- toolchain
cmd_toolchain() {
  [ $# = 1 ] || usage
  own "$1"
  local out rc v
  out=$(rexec "$ID" 600 "set -e
$RENV
rustup toolchain install $TC --profile minimal --no-self-update 2>&1 | tail -n 3
rustup default $TC 2>&1 | tail -n 1
echo \"nofile \$(ulimit -n)\"
echo \"RUSTC \$(rustc --version)\"
echo \"CARGO \$(cargo --version)\"" 2>&1); rc=$?
  printf '%s\n' "$out"
  v=$(sed -n 's/^RUSTC //p' <<< "$out" | tr -d '\r')
  [ $rc = 0 ] || fail $ROUTE_RC "$ID" toolchain "rc $rc: $(oneline "$out")"
  case "$v" in "rustc $TC "*) ;; *) fail $ROUTE_RC "$ID" toolchain "rustc --version printed '${v:-nothing}', not $TC";; esac
  note "$ID" toolchain ok "$v; $(sed -n 's/^nofile /ulimit -n /p' <<< "$out" | tr -d '\r')"
}

# ---------------------------------------------------------------- push
# guard LIST TRACKED: print every path the account rules keep out of a sandbox
guard() {
  python3 -I - "$1" "$2" << 'EOF'
import re, sys
paths = open(sys.argv[1], "rb").read().decode().split("\0")
tracked = set(open(sys.argv[2], "rb").read().decode().split("\0"))
seq = re.compile(r"\.(fasta|fna|fa|faa|ffn|frn|fas)(\.(gz|bz2|xz|zst))?$", re.I)
key = re.compile(r"^(\.env(\..*)?|.*\.(pem|key|p12|pfx|jks|keystore|ppk|asc|gpg)|id_(rsa|dsa|ecdsa|ed25519)(\.pub)?"
                 r"|\.netrc|\.npmrc|\.pypirc|\.git-credentials|credentials(\.json)?|.*secret.*\.(json|toml|ya?ml|txt))$", re.I)
for p in filter(None, paths):
    base = p.rsplit("/", 1)[-1]
    if p.startswith("data/"):
        print("data/ tree:", p)
    elif "/data/" in "/" + p and p not in tracked:
        print("untracked data/ file:", p)
    elif seq.search(base):
        print("sequence file:", p)
    elif key.match(base):
        print("env or key-like file:", p)
EOF
}
cmd_push() {
  local wt="" tmp=0 code=0 dest=$RDIR list tracked bad n lsha rsha out rc ts tarw rtar
  [ $# -ge 2 ] || usage
  own "$1"; shift
  while [ $# -gt 0 ]; do case "$1" in
    --with-tmp) tmp=1; shift;;
    --code-only) code=1; shift;;
    --dest) dest=${2:-}; shift 2 || usage;;
    -*) usage;;
    *) [ -z "$wt" ] || usage; wt=$1; shift;;
  esac; done
  [[ "$dest" =~ ^/(home/user|tmp)(/[A-Za-z0-9._-]+)+$ ]] || fail 2 "$ID" push "--dest must be a folder under /home/user or /tmp"
  case "$wt" in //*|\\\\*|*wsl.localhost*|*wsl\$*) fail 2 "$ID" push "give the worktree as a WSL path, not \\\\wsl.localhost";; esac
  wt=$(git -C "$wt" rev-parse --show-toplevel 2> /dev/null) || fail 2 "$ID" push "$wt is not a git worktree"
  stage_dir; ts=$(date -u +%Y%m%dT%H%M%SZ)
  list=$STAGE/files-$ts.lst tracked=$STAGE/tracked-$ts.lst
  trap 'rm -rf "$STAGE"' EXIT
  git -C "$wt" ls-files -z > "$tracked"
  if [ $code = 1 ]; then
    [ -z "$(git -C "$wt" status --porcelain --untracked-files=no)" ] || fail 2 "$ID" push "--code-only needs a clean worktree"
    git -C "$wt" ls-tree -r -z --name-only HEAD | python3 -I -c 'import sys
ps = sys.stdin.buffer.read().split(b"\0")
keep = [p for p in ps if p and not p.endswith(b".gz") and b"/data/" not in b"/" + p and not p.startswith(b".tmp/")]
sys.stdout.buffer.write(b"\0".join(keep) + b"\0")' > "$list.all"
  else
    git -C "$wt" ls-files -co --exclude-standard -z | python3 -I -c 'import sys
tmp = sys.argv[1] == "1"
ps = sys.stdin.buffer.read().split(b"\0")
keep = [p for p in ps if p and b"/target/" not in b"/" + p and (tmp or b"/.tmp/" not in b"/" + p)]
sys.stdout.buffer.write(b"\0".join(keep) + b"\0")' "$tmp" > "$list.all"
  fi
  # tracked files deleted in the worktree drop out; tar would fail on them
  (cd "$wt" && while IFS= read -r -d '' p; do [ -e "$p" ] || [ -L "$p" ] && printf '%s\0' "$p"; done < "$list.all") > "$list"
  n=$(tr -cd '\0' < "$list" | wc -c)
  [ "$n" -gt 0 ] || fail 2 "$ID" push "empty file set"
  bad=$(guard "$list" "$tracked")
  if [ -n "$bad" ]; then
    printf '%s\n' "$bad" | head -n 40 >&2
    fail $GUARD_RC "$ID" push "data guard refused $(wc -l <<< "$bad") paths, first: $(head -n 1 <<< "$bad"); nothing sent"
  fi
  tar -C "$wt" --null -T "$list" -cf "$STAGE/push-$ts.tar" || fail $ROUTE_RC "$ID" push "local tar failed"
  lsha=$(sha "$STAGE/push-$ts.tar")
  # boat scp can't read \\wsl.localhost paths, and C:\ reads as a host "C": copy from the staging dir by bare name
  tarw=push-$ts.tar rtar=/tmp/grm-push-$ts.tar
  out=$(cd "$STAGE" && BT=1200 bt scp "$tarw" "$ID:$rtar" 2>&1) || fail $ROUTE_RC "$ID" push "boat scp failed: $(oneline "$out")"
  out=$(rexec "$ID" 600 "set -e
rsha=\$(sha256sum $rtar | cut -d' ' -f1); echo \"RSHA \$rsha\"
[ \"\$rsha\" = $lsha ] || { rm -f $rtar; exit 3; }
mkdir -p $dest && tar -xf $rtar -C $dest && rm -f $rtar
echo \"FILES \$(find $dest -path $dest/target -prune -o -type f -print | wc -l)\"" 2>&1); rc=$?
  rsha=$(sed -n 's/^RSHA //p' <<< "$out" | tr -d '\r')
  [ $rc = 0 ] || fail $ROUTE_RC "$ID" push "remote check or extract failed (rc $rc), tar sha256 local $lsha remote ${rsha:-?}: $(oneline "$out")"
  rm -rf "$STAGE"; trap - EXIT
  note "$ID" "push $(basename "$wt")@$(git -C "$wt" rev-parse --short HEAD)$([ $code = 1 ] && echo ' code-only')$([ $tmp = 1 ] && echo ' with .tmp')" ok \
    "$n files -> $dest; tar sha256 local $lsha remote $rsha; staging removed"
  echo "pushed $n files, tar sha256 $lsha (both sides)"
}

# ---------------------------------------------------------------- push-data
# pd_tool DIR MANIFEST_COPY PIN TAR: PF-D2's checks, made on what is sent (REVIEW-PORT M1, L3, L4, I4). The
# manifest copy is read once; its bytes must hash to PIN, and those bytes are the list. Every line must be
# "<sha256>  <relative path>" for a path in an allowlisted public GRM data class with no denied segment, and DIR's own
# absolute path may hold no denied segment either ("review" excepted there: GRM's roots sit under it). Each listed
# file is opened component by component with O_NOFOLLOW (so no symlink anywhere on the way), must be a regular file,
# and is copied into TAR from that one open descriptor, as a plain regular member under its listed name; the
# manifest copy goes in as INPUT-MANIFEST.sha256. Then TAR is closed and read back: its members must be exactly the
# list plus the manifest, every one a regular file, each member's sha256 (hashed from the tar) equal to its line, and
# the manifest member equal to PIN; that read-back also hashes every byte of TAR. Prints refusals, one per line, and
# "OK <files> <bytes> <tar sha256>" as the last line only when nothing was refused; any other exit or output refuses
# (fail closed).
pd_tool() {
  python3 -I - "$@" << 'EOF'
import hashlib, io, os, re, stat, sys, tarfile
_, d, man, pin, out = sys.argv
problems = []
def done():
    for x in problems: print(x)
    sys.exit(1)
line_re = re.compile(r"([0-9a-f]{64})  ([^\0\n\r\\]+)")
seg = r"[a-z0-9][a-z0-9._-]*"
gid = r"[0-9]+\.[0-9]+"
classes = [
    r"fasta/" + gid + r"\.fna(\.gz)?",                              # BV-BRC FASTA
    r"bvbrc/(genome_amr|genomes)\.jsonl",                           # BV-BRC metadata
    r"(bvbrc/)?features/" + gid + r"\.json",                        # annotation cache
    r"curated/[a-z0-9_]+\.jsonl?", r"species-check/species_check\.json",
    r"bindings/[a-z0-9][a-z0-9._-]*\.(txt|json|tsv)", r"species/[a-z0-9_]+/genomes\.txt",   # curated
    r"catalog/(AMR_CDS\.fa|ReferenceGeneCatalog\.txt|version\.txt)",                        # catalog
    r"species/[a-z0-9_]+/(counted|structure)\.json", r"species/[a-z0-9_]+/matrix/[a-z0-9_.-]+",  # matrices
    r"species/[a-z0-9_]+/results/[A-Za-z0-9_.-]+\.json", r"report/([A-Za-z0-9_.-]+/){0,2}[A-Za-z0-9_.-]+\.(json|md)",  # reports
]
# GRM's own run outputs a later stage reads (main's ruling, 5 Nov): text only, at most 64 MiB each
receipts = [r"(audit|receipts|tables)/[A-Za-z0-9][A-Za-z0-9._-]*\.(json|jsonl|tsv|txt|sha256)",
            r"cohorts/[a-z0-9]+/validation-manifest\.json"]
pre = r"(%s/){0,3}" % seg
allow = re.compile(pre + "(%s)" % "|".join("(?:%s)" % c for c in classes + receipts))
data_only = re.compile(pre + "(%s)" % "|".join("(?:%s)" % c for c in classes))
receipt = re.compile(pre + "(%s)" % "|".join("(?:%s)" % c for c in receipts))
RECEIPT_MAX = 64 * 1024 * 1024
# the account rules' names, wherever they appear in a segment (I4), plus env and key files
deny = re.compile(r"sealed|held-?out|ks30|ks50|vumc|aua2027|gold-review|clinical|private|secret|credential|ksd|sigdrug"
                  r"|^review$|^layth_review$|^\.env|\.(pem|key|p12|pfx|jks|ppk|asc|gpg)$|^id_(rsa|dsa|ecdsa|ed25519)", re.I)
root = os.path.realpath(d)
for x in root.split("/"):
    if x and x.lower() != "review" and deny.search(x):
        problems.append("the data folder's own path has a refused segment %r: %s" % (x, root))
raw = open(man, "rb").read()
if hashlib.sha256(raw).hexdigest() != pin:
    problems.append("manifest sha256 %s is not the pinned %s" % (hashlib.sha256(raw).hexdigest(), pin)); done()
try:
    text = raw.decode("utf-8")
except UnicodeDecodeError:
    problems.append("manifest is not UTF-8"); done()
lines = text.split("\n")
if lines and lines[-1] == "": lines.pop()
if not lines: problems.append("manifest lists no file")
want, order = {}, []
for n, line in enumerate(lines, 1):
    m = line_re.fullmatch(line)
    if not m:
        problems.append("line %d is not '<sha256>  <path>' (binary-mode '*', CR or backslash refused): %r" % (n, line[:120])); continue
    h, p = m.groups(); parts = p.split("/")
    if p.startswith("/") or any(x in ("", ".", "..") for x in parts):
        problems.append("line %d: not a plain relative path: %s" % (n, p)); continue
    if p in want:
        problems.append("line %d: listed twice: %s" % (n, p)); continue
    bad = [x for x in parts if deny.search(x)]
    if bad:
        problems.append("line %d: refused segment %r: %s" % (n, bad[0], p)); continue
    if not allow.fullmatch(p) or p == "INPUT-MANIFEST.sha256":
        problems.append("line %d: not in a public GRM data class: %s" % (n, p)); continue
    want[p] = h; order.append(p)
if problems: done()
# build: each file read once, from a descriptor opened without following any symlink
def open_nofollow(rel):
    fd = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    try:
        parts = rel.split("/")
        for x in parts[:-1]:
            nfd = os.open(x, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=fd); os.close(fd); fd = nfd
        return os.open(parts[-1], os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=fd)
    finally:
        os.close(fd)
with tarfile.open(out, "w", format=tarfile.GNU_FORMAT) as t:
    for p in order:
        try:
            fd = open_nofollow(p)
        except OSError as e:
            problems.append("%s: cannot open without following a link (%s)" % (p, e.strerror)); continue
        with os.fdopen(fd, "rb") as f:
            st = os.fstat(f.fileno())
            if not stat.S_ISREG(st.st_mode):
                problems.append("%s: not a regular file" % p); continue
            if receipt.fullmatch(p) and not data_only.fullmatch(p) and st.st_size > RECEIPT_MAX:
                problems.append("%s: a run receipt over 64 MiB" % p); continue
            ti = tarfile.TarInfo(p); ti.size = st.st_size; ti.mode = 0o644; ti.mtime = 0; ti.type = tarfile.REGTYPE
            t.addfile(ti, f)
    ti = tarfile.TarInfo("INPUT-MANIFEST.sha256"); ti.size = len(raw); ti.mode = 0o644; ti.mtime = 0
    t.addfile(ti, io.BytesIO(raw))
if problems: done()
# verify what will be sent: the tar itself, read back once, front to back (stream mode, no seek), and every byte of
# that read hashed: the whole-tar sha256 printed below is of exactly the bytes whose members were checked (N1)
class Hashing:
    def __init__(self, f): self.f, self.h, self.n = f, hashlib.sha256(), 0
    def read(self, k=-1):
        b = self.f.read(k); self.h.update(b); self.n += len(b); return b
seen, nbytes = {}, 0
tfd = os.open(out, os.O_RDONLY | os.O_NOFOLLOW)
with os.fdopen(tfd, "rb") as raw_tar:
  size_now = os.fstat(raw_tar.fileno()).st_size
  hr = Hashing(raw_tar)
  with tarfile.open(fileobj=hr, mode="r|") as t:
    for mem in t:
      if not mem.isreg():
          problems.append("tar member %s is not a regular file" % mem.name); continue
      if mem.name in seen:
          problems.append("tar member %s twice" % mem.name); continue
      f = t.extractfile(mem); h = hashlib.sha256(); size = 0; head = b""
      for chunk in iter(lambda: f.read(1 << 20), b""):
          h.update(chunk); size += len(chunk)
          if receipt.fullmatch(mem.name) and not data_only.fullmatch(mem.name): head += chunk
      seen[mem.name] = h.hexdigest(); nbytes += size
      if mem.name != "INPUT-MANIFEST.sha256" and receipt.fullmatch(mem.name) and not data_only.fullmatch(mem.name):
          try: head.decode("utf-8")
          except UnicodeDecodeError: problems.append("%s: a run receipt that is not UTF-8 text" % mem.name)
  while hr.read(1 << 20): pass          # the end-of-archive blocks and padding: hashed too
  if hr.n != size_now: problems.append("the read-back covered %d of the tar's %d bytes" % (hr.n, size_now))
tar_sha = hr.h.hexdigest()
if set(seen) != set(want) | {"INPUT-MANIFEST.sha256"}:
    problems.append("tar members differ from the manifest: extra %s, missing %s" % (
        sorted(set(seen) - set(want) - {"INPUT-MANIFEST.sha256"})[:3], sorted(set(want) - set(seen))[:3]))
for p, h in want.items():
    if p in seen and seen[p] != h: problems.append("%s: the bytes in the tar hash %s, not the listed %s" % (p, seen[p], h))
if seen.get("INPUT-MANIFEST.sha256") != pin: problems.append("the manifest in the tar is not the pinned one")
if problems: done()
print("OK %d %d %s" % (len(want), nbytes, tar_sha))
EOF
}
# PD_BOX: push-data's box side (N1), run with d, t, vsha, want and n set before it
PD_BOX='set -u
made=0
bye() { rm -f "$t"; [ $made = 1 ] && rm -rf "$d"; echo "REFUSED $1"; exit $2; }
rsha=$(sha256sum "$t" | cut -d" " -f1); echo "RSHA $rsha"
[ "$rsha" = "$vsha" ] || bye "tar sha256 is not the verified one" 3
if [ -e "$d" ] && [ -n "$(ls -A "$d")" ]; then rm -f "$t"; echo "REFUSED dest not empty"; exit 4; fi
mkdir -p "$d" || bye "mkdir failed" 5
made=1
tar --no-same-owner --no-overwrite-dir -xf "$t" -C "$d" || bye "extract failed" 5
rm -f "$t"
[ -f "$d/INPUT-MANIFEST.sha256" ] && [ ! -L "$d/INPUT-MANIFEST.sha256" ] || bye "no manifest in the tar" 5
mv "$d/INPUT-MANIFEST.sha256" "$d/.grm-input-manifest.sha256" || bye "manifest rename failed" 5
m=$(sha256sum "$d/.grm-input-manifest.sha256" | cut -d" " -f1); echo "MSHA $m"
[ "$m" = "$want" ] || bye "manifest is not the pinned one" 6
[ -z "$(find "$d" -mindepth 1 ! -type f ! -type d | head -n 1)" ] || bye "a link or special file came out of the tar" 6
c=$(wc -l < "$d/.grm-input-manifest.sha256")
[ "$c" = "$n" ] || bye "manifest lists $c files, expected $n" 6
got=$(cd "$d" && find . -mindepth 1 -type f -printf "%P\n" | LC_ALL=C sort)
exp=$( { cut -c67- "$d/.grm-input-manifest.sha256"; echo .grm-input-manifest.sha256; } | LC_ALL=C sort)
[ "$got" = "$exp" ] || bye "the extracted files are not exactly the manifest" 6
gotd=$(cd "$d" && find . -mindepth 1 -type d -printf "%P\n" | LC_ALL=C sort)
expd=$(cut -c67- "$d/.grm-input-manifest.sha256" | awk -F/ "{p=\"\"; for (i = 1; i < NF; i++) { p = (i == 1 ? \$1 : p \"/\" \$i); print p } }" | LC_ALL=C sort -u)
[ "$gotd" = "$expd" ] || bye "the extracted folders are not exactly the manifest" 6
(cd "$d" && sha256sum -c --strict --quiet .grm-input-manifest.sha256) || bye "sha256sum -c failed in the box" 6
echo "CHECKED $c"'
cmd_push_data() {
  local dir="" man="" want="" dest=/home/user/grm-data out rc n bytes lsha vsha rsha ts tarw rtar avail last priv
  [ $# -ge 2 ] || usage
  own "$1"; shift
  while [ $# -gt 0 ]; do case "$1" in
    --manifest) man=${2:-}; shift 2 || usage;;
    --manifest-sha256) want=${2:-}; shift 2 || usage;;
    --dest) dest=${2:-}; shift 2 || usage;;
    -*) usage;;
    *) [ -z "$dir" ] || usage; dir=$1; shift;;
  esac; done
  [ -n "$dir" ] && [ -n "$man" ] || usage
  [[ "$want" =~ ^[0-9a-f]{64}$ ]] || fail 2 "$ID" push-data "--manifest-sha256 takes the manifest's 64-hex sha256, written in NOTES.md before the push"
  [[ "$dest" =~ ^/(home/user|tmp)(/[A-Za-z0-9._-]+)+$ ]] && [[ "/$dest/" != */../* ]] && [[ "/$dest/" != */./* ]] ||
    fail 2 "$ID" push-data "--dest must be a folder under /home/user or /tmp, without . or .. segments"
  case "$dir$man" in *wsl.localhost*|*wsl\$*|//*|\\\\*) fail 2 "$ID" push-data "give WSL paths, not \\\\wsl.localhost";; esac
  [ -d "$dir" ] && [ -f "$man" ] || fail 2 "$ID" push-data "$dir is not a directory or $man is not a file"
  dir=$(readlink -f "$dir")
  case "${dir,,}" in */ksd*|*sigdrug*) fail $GUARD_RC "$ID" push-data "$dir belongs to another project; nothing sent";; esac
  stage_dir; ts=$(date -u +%Y%m%dT%H%M%SZ)
  # N1: the tar is built and verified in a private folder only this user can enter; the Windows staging folder only
  # ever holds a copy made from one hashed stream
  priv=$(umask 077; mktemp -d "${TMPDIR:-/tmp}/grm-push-data.XXXXXXXX") || fail $ROUTE_RC "$ID" push-data "mktemp -d failed"
  chmod 700 "$priv"
  trap 'rm -rf "$STAGE" "$priv"' EXIT
  # the manifest is copied once; only the copy is read from here on
  cp -- "$man" "$priv/manifest.sha256" || fail $ROUTE_RC "$ID" push-data "could not copy the manifest"
  out=$(pd_tool "$dir" "$priv/manifest.sha256" "$want" "$priv/data.tar" 2>&1); rc=$?
  last=$(tail -n 1 <<< "$out")
  if [ $rc != 0 ] || [[ ! "$last" =~ ^OK\ [0-9]+\ [0-9]+\ [0-9a-f]{64}$ ]]; then
    printf '%s\n' "$out" | head -n 40 >&2
    fail $GUARD_RC "$ID" push-data "refused (checker rc $rc, $(grep -vc '^OK ' <<< "$out") problems), first: $(oneline "$(head -n 1 <<< "$out")" 200); nothing sent"
  fi
  read -r _ n bytes vsha <<< "$last"
  # free space is measured now that the tar's size is known (I-f)
  avail=$(df -B1 --output=avail "$STAGE" | tail -n 1 | tr -d ' ')
  [ "$avail" -gt $(( $(stat -c %s "$priv/data.tar") + bytes / 10 + 4 * 1024 ** 3 )) ] ||
    fail $ROUTE_RC "$ID" push-data "staging has $avail B free for a $(stat -c %s "$priv/data.tar") B tar; nothing sent"
  # one read of the verified tar feeds both the staged copy and its sha256; it must be the verified one
  lsha=$(tee "$STAGE/data-$ts.tar" < "$priv/data.tar" | sha256sum | cut -d' ' -f1)
  [ "$lsha" = "$vsha" ] || fail $GUARD_RC "$ID" push-data "the tar changed after it was verified (copied $lsha, verified $vsha); nothing sent"
  rm -f "$priv/data.tar"
  tarw=data-$ts.tar rtar=/tmp/grm-data-$ts.tar
  # at least 1 MB/s through WSL interop and the Boat API, on top of the usual 1200 s
  out=$(cd "$STAGE" && BT=$((1200 + bytes / 1000000)) bt scp "$tarw" "$ID:$rtar" 2>&1) || {
    rexec "$ID" 60 "rm -f $rtar" > /dev/null 2>&1; fail $ROUTE_RC "$ID" push-data "boat scp failed: $(oneline "$out")"; }
  # in the box: the tar must hash to the verified sha256 before anything is extracted; then the manifest must be the
  # pinned one, every extracted file and folder must be exactly the manifest's (its own name excluded at the top
  # only), no link or special file, and every file re-hashed; any failure removes DEST and the tar
  out=$(rexec "$ID" 600 "d=$(printf %q "$dest") t=$rtar vsha=$vsha want=$want n=$n
$PD_BOX" 2>&1); rc=$?
  rsha=$(sed -n 's/^RSHA //p' <<< "$out" | tr -d '\r')
  [ $rc = 0 ] || fail $ROUTE_RC "$ID" push-data "refused in the box (rc $rc, destination and tar removed): $(oneline "$(grep -E '^(REFUSED|RSHA)' <<< "$out" | tr '\n' ' ')")"
  [ "$(sed -n 's/^CHECKED //p' <<< "$out" | tr -d '\r')" = "$n" ] || fail $ROUTE_RC "$ID" push-data "the box did not confirm $n checked files: $(oneline "$out")"
  rm -rf "$STAGE" "$priv"; trap - EXIT
  note "$ID" "push-data $(basename "$dir")" ok "$n files, $bytes B -> $dest; tar built in a private folder and verified member by member against manifest $want; verified tar sha256 $vsha = staged stream $lsha = box $rsha; exact file set and sha256sum -c --strict in the box; staging removed"
  echo "pushed $n data files ($bytes B) to $dest, manifest $want, tar sha256 $lsha (both sides)"
}

# ---------------------------------------------------------------- run
cmd_run() {
  local log="" cwd=$RDIR cmd out pid rc="" st sx="" t0 fails=0 secs="" tm polls=0
  local force_argv=0 mode=shell quoted
  local -a argv
  [ $# -ge 2 ] || usage
  own "$1"; shift
  while [ $# -gt 0 ]; do case "$1" in
    --log) log=${2:-}; shift 2 || usage;;
    --cwd) cwd=${2:-}; shift 2 || usage;;
    --argv) force_argv=1; shift;;
    --) shift; break;;
    *) break;;
  esac; done
  argv=("$@")
  [ ${#argv[@]} -gt 0 ] && [ -n "${argv[0]}" ] || usage
  [ $force_argv = 0 ] && [ ${#argv[@]} = 1 ] || mode=argv
  cmd="$*"
  printf -v quoted '%q ' "${argv[@]}"
  log=${log:-run-$(date -u +%Y%m%dT%H%M%SZ)}
  [[ "$log" =~ ^[A-Za-z0-9._-]+$ ]] || fail 2 "$ID" run "--log takes a plain file name"
  # the launcher writes NAME.rc ("rc seconds") when the command ends, so the exit code survives a lost poll
  out=$(rexec "$ID" 60 "mkdir -p $RLOG && rm -f $RLOG/$log.rc && cat > $RLOG/$log.sh << 'GRMEOF'
$RENV
cwd=$(printf %q "$cwd")
if ! cd "\$cwd" 2>$RLOG/$log.log; then
  printf 'grm-boat: run: cannot enter working directory %s; push the worktree to this directory or pass --cwd DIR for an existing directory\n' "\$cwd" >>$RLOG/$log.log
  echo '86 0' > $RLOG/$log.rc
  exit 86
fi
run_mode=$mode
argv=($quoted)
if [ \"\$run_mode\" = argv ]; then set -- \"\${argv[@]}\"; else set -- bash -c \"\${argv[0]}\"; fi
t0=\$(date +%s)
if [ -x /usr/bin/time ]; then /usr/bin/time -v -o $RLOG/$log.time -- \"\$@\" > $RLOG/$log.log 2>&1
else ( exec -- \"\$@\" ) > $RLOG/$log.log 2>&1; fi
rc=\$?
echo \"\$rc \$(( \$(date +%s) - t0 ))\" > $RLOG/$log.rc
exit \$rc
GRMEOF
echo ok" 2>&1) || fail $ROUTE_RC "$ID" "run $log" "could not write the launcher: $(oneline "$out")"
  t0=$(date +%s)
  out=$(bt exec "$ID" --detach --json -- bash "$RLOG/$log.sh" 2>&1) || fail $ROUTE_RC "$ID" "run $log" "boat exec --detach failed: $(oneline "$out")"
  pid=$(printf '%s' "$out" | python3 -I -c 'import json, sys, re
t = sys.stdin.read()
try:
    d = json.loads(t); print(d.get("pid") or d.get("processId") or d.get("id") or "")
except Exception:
    m = re.search(r"\d+", t); print(m.group(0) if m else "")' 2> /dev/null)
  [ -n "$pid" ] || fail $ROUTE_RC "$ID" "run $log" "no pid from boat exec --detach: $(oneline "$out")"
  echo "grm-boat: $log running as pid $pid on $ID: $cmd" >&2
  while :; do
    # 5 s steps for the first 30 s, so a short step doesn't wait out a 15 s poll
    polls=$((polls + 1)); if [ $polls -le 6 ]; then sleep 5; else sleep 15; fi
    if out=$(bt exec "$ID" --status "$pid" --json 2>&1); then
      fails=0
      # {"running": false, "status": "exited", "exitCode": N, ...}; "success" only says the launch worked
      read -r st sx <<< "$(printf '%s' "$out" | python3 -I -c 'import json, sys
d = json.load(sys.stdin)
r = d.get("running")
if r is None: r = str(d.get("status") or "").lower() in ("running", "started", "starting")
print("running" if r else "ended", d.get("exitCode"))' 2> /dev/null || echo unknown)"
      [ "$st" = ended ] && break
      if [ "$st" = unknown ] && rexec "$ID" 30 "test -s $RLOG/$log.rc" > /dev/null 2>&1; then break; fi
    else
      fails=$((fails + 1))
      if [ $fails -ge 8 ]; then
        binfo "$ID" > /dev/null; [ $? = 3 ] && fail $ROUTE_RC "$ID" "run $log" "the sandbox is gone (404) while pid $pid ran"
        fail $ROUTE_RC "$ID" "run $log" "8 status polls in a row failed: $(oneline "$out")"
      fi
    fi
  done
  out=$(rexec "$ID" 60 "for t in 1 2 3 4 5; do [ -s $RLOG/$log.rc ] && break; sleep 2; done
echo \"RC \$(cat $RLOG/$log.rc 2>/dev/null)\"
echo \"LOGSHA \$(sha256sum $RLOG/$log.log | cut -d' ' -f1)\"
echo \"PEAK \$(sed -n 's/.*Maximum resident set size (kbytes): //p' $RLOG/$log.time 2>/dev/null)\"
echo ----; tail -n 40 $RLOG/$log.log" 2>&1)
  printf '%s\n' "$out" | sed -n '/^----/,$p' | sed 1d
  read -r rc secs <<< "$(sed -n 's/^RC //p' <<< "$out" | tr -d '\r')"
  tm=$(sed -n 's/^PEAK //p' <<< "$out" | tr -d '\r'); [ -n "$tm" ] && tm="; peak RSS $((tm / 1024)) MiB"
  # no rc file (killed before the launcher could write it): fall back on the exit code --status reported
  [ -n "${rc:-}" ] || { [[ "$sx" =~ ^[0-9]+$ ]] && rc=$sx secs="?"; }
  [ -n "${rc:-}" ] || fail $ROUTE_RC "$ID" "run $log" "pid $pid ended but left no rc file and --status gave no exit code"
  note "$ID" "run $log: $(oneline "$cmd" 120)" "rc=$rc" "${secs}s on the sandbox ($(( $(date +%s) - t0 ))s wall)$tm; $RLOG/$log.log sha256 $(sed -n 's/^LOGSHA //p' <<< "$out" | tr -d '\r')"
  echo "grm-boat: $log ended rc=$rc after ${secs}s" >&2
  return "$rc"
}

# ---------------------------------------------------------------- pull
cmd_pull() {
  [ $# = 3 ] || usage
  own "$1"
  local rp=$2 ld=$3 ts rtar out rc rsha lsha n
  stage_dir; ts=$(date -u +%Y%m%dT%H%M%SZ); rtar=/tmp/grm-pull-$ts.tar
  trap 'rm -rf "$STAGE"' EXIT
  out=$(rexec "$ID" 600 "set -e; p=$(printf %q "$rp"); [ -e \"\$p\" ]
tar -C \"\$(dirname \"\$p\")\" -cf $rtar \"\$(basename \"\$p\")\"
echo \"RSHA \$(sha256sum $rtar | cut -d' ' -f1)\"" 2>&1) || fail $ROUTE_RC "$ID" "pull $rp" "remote tar failed: $(oneline "$out")"
  rsha=$(sed -n 's/^RSHA //p' <<< "$out" | tr -d '\r')
  out=$(cd "$STAGE" && BT=1200 bt scp "$ID:$rtar" "pull-$ts.tar" 2>&1) || fail $ROUTE_RC "$ID" "pull $rp" "boat scp failed: $(oneline "$out")"
  lsha=$(sha "$STAGE/pull-$ts.tar")
  [ "$lsha" = "$rsha" ] || fail $ROUTE_RC "$ID" "pull $rp" "tar sha256 local $lsha remote $rsha differ"
  mkdir -p "$ld" && tar -xf "$STAGE/pull-$ts.tar" -C "$ld" --no-same-owner || fail $ROUTE_RC "$ID" "pull $rp" "local extract failed"
  n=$(tar -tf "$STAGE/pull-$ts.tar" | grep -vc '/$')
  rexec "$ID" 60 "rm -f $rtar" > /dev/null 2>&1
  rm -rf "$STAGE"; trap - EXIT
  note "$ID" "pull $rp" ok "$n files -> $ld; tar sha256 local $lsha remote $rsha"
  echo "pulled $n files into $ld, tar sha256 $lsha (both sides)"
}

# ---------------------------------------------------------------- snapshot
cmd_snapshot() {
  [ $# = 2 ] || usage
  own "$1"
  local nm=$2 out
  case "$NAME" in grm-tmpl-*) ;; *) fail 2 "$ID" snapshot "only grm-tmpl-* sandboxes are saved as templates";; esac
  case "$nm" in grm-tmpl-*) ;; *) fail 2 "$ID" snapshot "a template's name starts with grm-tmpl-";; esac
  out=$(BT=1800 bt snapshot "$ID" "$nm" --json 2>&1) || fail $ROUTE_RC "$ID" "snapshot $nm" "$(oneline "$out")"
  printf '%s\n' "$out"
  note "$ID" "snapshot $nm" ok "$(oneline "$out" 300)"
}

# ---------------------------------------------------------------- close
cmd_close() {
  [ $# = 2 ] || usage
  NOTES=$(readlink -f "$2"); [ -n "$NOTES" ] || usage
  own "$1"
  local u us rc i gone=0 out
  u=$(bt_retry usage "$ID" --json) || fail $ROUTE_RC "$ID" close "boat usage failed: $(oneline "$u")"
  printf '\nboat usage %s (%s), %s:\n```\n%s\n```\n' "$ID" "$NAME" "$(now)" "$u" >> "$NOTES"
  # "seconds" are billed VM-seconds: wall time times billingMultiplier (5.56 for xlarge)
  us=$(printf '%s' "$u" | python3 -I -c 'import json, sys
d = json.load(sys.stdin); m = d.get("billingMultiplier") or 1; s = d.get("seconds") or 0; w = round(s / m)
print("%s %d s wall (%dm%02ds), %d VM-s billed at x%.2f, $%.4f" % (d.get("sandboxType"), w, w // 60, w % 60, s, m, d.get("dollars") or 0))' 2> /dev/null) ||
    us=$(oneline "$u" 160)
  printf '%s\n' "$us"
  out=$(bt delete "$ID" --yes 2>&1) || fail $ROUTE_RC "$ID" close "boat delete --yes failed: $(oneline "$out")"
  for i in $(seq 1 20); do binfo "$ID" > /dev/null; rc=$?; [ $rc = 3 ] && { gone=1; break; }; sleep 6; done
  [ $gone = 1 ] || fail $ROUTE_RC "$ID" close "boat delete ran, but boat info still answers (rc $rc)"
  out=$(list_live) && awk -v i="$ID" '$1 == i' <<< "$out" | grep -q . && fail $ROUTE_RC "$ID" close "boat list still shows $ID"
  out=$(grep -o 'not found (404)' <<< "$INFO_ERR" || echo "$INFO_ERR")
  note "$ID" "close $NAME" ok "usage: $us; deleted --yes; boat info: $out; boat list shows no $NAME"
  echo "deleted $ID ($NAME); boat info: $out"
}

# ---------------------------------------------------------------- main
if [ "${1:-}" = --notes ]; then NOTES=${2:-}; shift 2 || usage; fi
cmd=${1:-}; shift || usage
if [ "$cmd" != close ] && [ "$cmd" != preflight ]; then
  [ -n "$NOTES" ] || { echo "grm-boat: name a NOTES.md with --notes PATH or GRM_NOTES" >&2; exit 2; }
  mkdir -p "$(dirname "$NOTES")" && NOTES=$(readlink -f "$NOTES")
fi
case "$cmd" in
  new|toolchain|push|run|pull|snapshot|close|preflight) "cmd_$cmd" "$@";;
  push-data) cmd_push_data "$@";;
  *) usage;;
esac
