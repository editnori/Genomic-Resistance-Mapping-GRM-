# shellcheck shell=bash
# Shared checks for the GRM deployment scripts. Source it; do not run it.
#
# Configuration comes only from the environment, so no script names a host, a
# home directory, a run date or a credential:
#   GRM_DEPLOY_ROOT      absolute deployment root (required)
#   GRM_PORT             live loopback port (default 8641)
#   GRM_EXPECT_SCHEMA    report index schema (default 2)
#   GRM_EXPECT_REVISION  report evaluation revision (default 7)
#   GRM_EXPECT_DATASETS  dataset report count (default 237)
# GRM_EXPECT_* is what publish.sh, stage.sh and approve-legacy.sh require of a
# report. The approval they lead to records those values, and activate.sh,
# rollback.sh, restore.sh and status.sh check each pair against its own
# approval, so a 275-dataset candidate and a 237-dataset rollback target need
# no shared setting.
#
# Layout under the root:
#   published/NAME/{report/,SOURCE.sha256}       immutable report publications
#   release/LABEL/{grm,dist/,COMMIT,MANIFEST.sha256}   immutable builds
#   src-LABEL/, src-LABEL.commit                 extracted source of a build
#   stage/LABEL--NAME--TIME/                     candidate check receipts
#   go/LABEL--NAME                               approval of one build and one report copy
#   state/active, state/previous                 what activate.sh switched
#   serve-PORT.starting                          journal of a start in progress
#   serve-PORT.pid                               the committed server on PORT
#   spawn/TOKEN.XXXXXX                           marker held by the processes of one start
#   deploy.log, deploy.lock, logs/
#
# The scripts that deploy change the root or the service only with --apply.
# Without it they run every read-only check, print what they would do, and
# write nothing under the root.
#
# Test hooks: GRM_TEST_SIGNAL_AT and GRM_TEST_CHILD_AT (lib.sh checkpoint and
# child_checkpoint). They do nothing unless set.

set -euo pipefail

D=${GRM_DEPLOY_ROOT:?set GRM_DEPLOY_ROOT to the absolute deployment root}
D=${D%/}
case $D in
    /?*) ;;
    *) echo "GRM_DEPLOY_ROOT must be an absolute path: $D" >&2; exit 2 ;;
esac
# /proc reports a process's cwd and executable by physical path.
if [ -d "$D" ]; then D=$(cd -- "$D" && pwd -P); fi
PORT=${GRM_PORT:-8641}
EXPECT_SCHEMA=${GRM_EXPECT_SCHEMA:-2}
EXPECT_REVISION=${GRM_EXPECT_REVISION:-7}
EXPECT_DATASETS=${GRM_EXPECT_DATASETS:-237}
for value in "$PORT" "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS"; do
    [[ $value =~ ^[0-9]+$ ]] || { echo "not a number: $value" >&2; exit 2; }
done
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/lib.sh"
APPLY=0
LOCKED=0
LOCK_FD=''
ARGS=()
declare -A APPROVED=() FIELDS=()
APPROVAL_KEYS="stage browser commit publication release_manifest report_manifest schema revision datasets"
STATE_KEYS="label publication commit release_manifest report_manifest schema revision datasets"
PF_ERR=''
START_PID=''
SPAWN_TOKEN='' SPAWN_JOURNAL=0 SPAWN_MARK='' SPAWN_MARK_ID='' SPAWN_CHILD='' SPAWN_CHILD_START=''
MARK_FD=''

# Split --apply from the positional arguments: parse_args "$@"; set -- "${ARGS[@]}"
parse_args() {
    local arg
    for arg in "$@"; do
        case $arg in
            --apply) APPLY=1 ;;
            --*) echo "unknown option $arg" >&2; exit 2 ;;
            *) ARGS+=("$arg") ;;
        esac
    done
}

# The log copy is best effort: a full disk must not stop a rollback, and a
# deploy.log that is a link, FIFO or directory is never opened (append_line).
say() {
    local line
    line="$(date -u +%FT%TZ) $*"
    printf '%s\n' "$line"
    if ((APPLY)) && [ -d "$D" ]; then append_line "$D/deploy.log" "$line"; fi
}
die() {
    say "FAIL: $1" >&2
    exit "${2:-1}"
}
would() { printf 'dry-run: would %s\n' "$*"; }

need() {
    local tool
    for tool in "$@"; do
        command -v "$tool" >/dev/null || die "missing required tool $tool" 2
    done
}

# Labels and publication names become path components, so keep them plain.
valid_name() { [[ $1 =~ ^[A-Za-z0-9][A-Za-z0-9._-]*$ ]] && [[ $1 != *--* ]]; }
check_name() {
    valid_name "$1" || die "invalid name '$1'" 2
}

# DIR holds only directories and regular files, every file is listed once in
# MANIFEST with a matching SHA-256, and nothing listed is missing. Links,
# devices, FIFOs and sockets are refused. Only a release manifest, stored as
# DIR/MANIFEST.sha256, is itself left out of the listing. Sets VT_LISTED.
verify_tree() {
    local dir=$1 manifest=$2 actual odd
    local skip=()
    [ -d "$dir" ] && [ ! -L "$dir" ] || die "missing directory $dir"
    [ -f "$manifest" ] && [ ! -L "$manifest" ] && [ -s "$manifest" ] || die "missing manifest $manifest"
    odd=$(cd "$dir" && find . ! -type f ! -type d -print -quit)
    [ -z "$odd" ] || die "$dir holds a link or special file: ${odd#./}"
    if [ "$manifest" = "$dir/MANIFEST.sha256" ]; then skip=(! -path ./MANIFEST.sha256); fi
    VT_LISTED=$(sed -E 's/^[0-9a-f]{64} [ *]//; s#^\./##' "$manifest" | LC_ALL=C sort)
    actual=$(cd "$dir" && find . -type f "${skip[@]}" | sed 's#^\./##' | LC_ALL=C sort)
    # Compared before hashing, so sha256sum reads nothing outside DIR.
    [ "$VT_LISTED" = "$actual" ] || die "$dir holds files its manifest does not list, or lacks listed files"
    (cd "$dir" && sha256sum --quiet --strict -c -- "$manifest") || die "$dir does not match $manifest"
}

# The report's index has schema $2, revision $3 and $4 datasets, and lists
# exactly the dataset reports present.
check_report_index() {
    python3 -I - "$@" <<'PY' || die "report $1 does not match schema $2, revision $3, $4 datasets"
import json, os, sys
root = sys.argv[1]
schema, revision, count = map(int, sys.argv[2:])
with open(os.path.join(root, "index.json")) as f:
    index = json.load(f)
ids = [d["id"] for d in index["datasets"]]
files = sorted(f[:-5] for f in os.listdir(os.path.join(root, "datasets")) if f.endswith(".json"))
problems = []
if index.get("schema") != schema:
    problems.append(f"index schema {index.get('schema')!r}, expected {schema}")
if index.get("evaluation_revision") != revision:
    problems.append(f"evaluation_revision {index.get('evaluation_revision')!r}, expected {revision}")
if len(ids) != count or len(set(ids)) != count:
    problems.append(f"index lists {len(ids)} datasets ({len(set(ids))} distinct), expected {count}")
if sorted(ids) != files:
    problems.append(f"index datasets differ from the {len(files)} files in datasets/")
for problem in problems:
    print(problem, file=sys.stderr)
sys.exit(1 if problems else 0)
PY
}

publication_tree() {
    check_name "$1"
    local P=$D/published/$1
    [ -d "$P" ] && [ ! -L "$P" ] || die "missing directory $P"
    verify_tree "$P/report" "$P/SOURCE.sha256"
}

# verify_publication NAME SCHEMA REVISION DATASETS
verify_publication() {
    publication_tree "$1"
    check_report_index "$D/published/$1/report" "$2" "$3" "$4"
}

verify_release() {
    check_name "$1"
    local R=$D/release/$1 file
    [ -d "$R" ] && [ ! -L "$R" ] || die "release $1 is missing"
    verify_tree "$R" "$R/MANIFEST.sha256"
    for file in grm COMMIT dist/index.html; do
        grep -qxF "$file" <<<"$VT_LISTED" || die "release $1 is incomplete: its manifest does not list $file"
    done
    [ -x "$R/grm" ] || die "release $1 grm is not executable"
    release_commit "$1" >/dev/null
}

# These read a release or publication that may not have been verified yet
# (status.sh), so a link or special file is refused rather than opened.
release_commit() {
    local commit
    commit=$(safe_line "$D/release/$1/COMMIT" 100) || commit=''
    [[ $commit =~ ^[0-9a-f]{40}$ ]] || die "release $1 has no full commit id"
    printf '%s\n' "$commit"
}
release_digest() { safe_sha256 "$D/release/$1/MANIFEST.sha256"; }
report_digest() { safe_sha256 "$D/published/$1/SOURCE.sha256"; }

# parse_fields FILE line|lines: parse FILE's key=value fields into FIELDS, or
# return 1 with the reason in PF_ERR. The whole file is read once, without
# following a link or blocking on a FIFO (lib.sh GRM_PY), and must be a regular
# file of at most 64 KiB holding only printable ASCII and newlines, ending with
# a newline: a NUL, CR, tab, other control byte or non-ASCII byte anywhere
# refuses it before any field is parsed. "line" is one line of fields
# separated by single spaces; "lines" is one field per line with no blank
# line. Every field must match key=value with a plain value, and no key may
# repeat.
parse_fields() {
    local out key value
    FIELDS=()
    PF_ERR=''
    out=$(python3 -I -c "$GRM_PY" fields "$1" "$2") || { PF_ERR=${out#ERR }; return 1; }
    while read -r key value; do FIELDS[$key]=$value; done <<<"$out"
}

# stage.sh approves one build of a release with one copy of a publication in
# go/LABEL--NAME, a single line of unique fields:
#   stage=0 browser=pass|skipped commit=SHA publication=NAME
#   release_manifest=SHA256 report_manifest=SHA256 schema=N revision=N
#   datasets=N [receipts=PATH]
# approve-legacy.sh writes browser=legacy with legacy_sha256=SHA256 instead of
# receipts. Every field is compared exactly, and the report's expected schema,
# revision and dataset count are kept in APPROVED[LABEL--NAME]. A hand-written
# GO-LABEL alone approves nothing.
check_go() {
    local label=$1 pub=$2 file=$D/go/$1--$2 key
    if [ ! -e "$file" ] && [ ! -L "$file" ]; then
        [ ! -e "$D/GO-$label" ] || die "no stage approval for $label with $pub: $D/GO-$label is a legacy approval that names no publication or report copy; stage the pair, or run approve-legacy.sh $label $pub while it serves"
        die "no stage approval for $label with $pub; run stage.sh"
    fi
    parse_fields "$file" line || die "$file is not a one-line approval of unique key=value fields: it $PF_ERR"
    for key in "${!FIELDS[@]}"; do
        case " $APPROVAL_KEYS receipts legacy_sha256 " in
            *" $key "*) ;;
            *) die "$file has unknown field $key" ;;
        esac
    done
    for key in $APPROVAL_KEYS; do
        [ -n "${FIELDS[$key]:-}" ] || die "$file lacks $key"
    done
    [ "${FIELDS[stage]}" = 0 ] || die "$file does not approve: stage=${FIELDS[stage]}"
    case ${FIELDS[browser]} in
        pass | skipped) ;;
        legacy) [[ ${FIELDS[legacy_sha256]:-} =~ ^[0-9a-f]{64}$ ]] || die "$file is a legacy migration without legacy_sha256" ;;
        *) die "$file has browser=${FIELDS[browser]}" ;;
    esac
    [ "${FIELDS[commit]}" = "$(release_commit "$label")" ] || die "$file does not approve commit $(release_commit "$label")"
    [ "${FIELDS[publication]}" = "$pub" ] || die "$file approves another publication"
    [ "${FIELDS[release_manifest]}" = "$(release_digest "$label")" ] || die "$file approves another build of $label"
    [ "${FIELDS[report_manifest]}" = "$(report_digest "$pub")" ] || die "$file approves another copy of $pub"
    for key in schema revision datasets; do
        [[ ${FIELDS[$key]} =~ ^[0-9]+$ ]] || die "$file has $key=${FIELDS[$key]}"
    done
    APPROVED[$label--$pub]="${FIELDS[schema]} ${FIELDS[revision]} ${FIELDS[datasets]}"
}

# An intact, approved pair whose report matches what its approval recorded.
verify_target() {
    local s r n
    verify_release "$1"
    publication_tree "$2"
    check_go "$1" "$2"
    read -r s r n <<<"${APPROVED[$1--$2]}"
    check_report_index "$D/published/$2/report" "$s" "$r" "$n"
}

# state/NAME records a pair with what its approval bound: one key=value per
# line for each of STATE_KEYS, then at=TIME.
expected_record() {
    local s r n
    read -r s r n <<<"${APPROVED[$1--$2]:?no approval loaded for $1 with $2}"
    printf 'label=%s publication=%s commit=%s release_manifest=%s report_manifest=%s schema=%s revision=%s datasets=%s ' \
        "$1" "$2" "$(release_commit "$1")" "$(release_digest "$1")" "$(report_digest "$2")" "$s" "$r" "$n"
}
fields_record() {
    local key
    for key in $STATE_KEYS; do printf '%s=%s ' "$key" "${FIELDS[$key]:-}"; done
}

# Read state/NAME into ST_LABEL, ST_PUB and ST_RECORD.
read_state() {
    local file=$D/state/$1 key
    [ -e "$file" ] || [ -L "$file" ] || die "no state/$1 recorded"
    parse_fields "$file" lines || die "state/$1 is malformed: it $PF_ERR"
    for key in "${!FIELDS[@]}"; do
        case " $STATE_KEYS at " in
            *" $key "*) ;;
            *) die "state/$1 has unknown field $key" ;;
        esac
    done
    for key in $STATE_KEYS at; do
        [ -n "${FIELDS[$key]:-}" ] || die "state/$1 lacks $key"
    done
    ST_LABEL=${FIELDS[label]}
    ST_PUB=${FIELDS[publication]}
    check_name "$ST_LABEL"
    check_name "$ST_PUB"
    ST_RECORD=$(fields_record)
}

# Read state/NAME, verify its pair, and require the record to match that
# pair's build, report copy and approval.
verify_state() {
    read_state "$1"
    verify_target "$ST_LABEL" "$ST_PUB"
    [ "$ST_RECORD" = "$(expected_record "$ST_LABEL" "$ST_PUB")" ] ||
        die "state/$1 does not match the current build, report copy or approval of $ST_LABEL with $ST_PUB"
}

# Whether state/NAME already records LABEL with PUBLICATION (loaded approval).
state_is() {
    parse_fields "$D/state/$1" lines || return 1
    [ -n "${FIELDS[at]:-}" ] && [ "${#FIELDS[@]}" = 9 ] &&
        [ "$(fields_record)" = "$(expected_record "$2" "$3")" ]
}

# Record LABEL with NAME in state/FILE and read it back: the result must be a
# regular, readable file holding exactly the expected fields.
write_state() {
    local s r n text
    read -r s r n <<<"${APPROVED[$2--$3]:?no approval loaded for $2 with $3}"
    text=$(printf 'label=%s\npublication=%s\ncommit=%s\nrelease_manifest=%s\nreport_manifest=%s\nschema=%s\nrevision=%s\ndatasets=%s\nat=%s' \
        "$2" "$3" "$(release_commit "$2")" "$(release_digest "$2")" "$(report_digest "$3")" "$s" "$r" "$n" "$(date -u +%FT%TZ)") ||
        return 1
    [ -e "$D/state" ] || [ -L "$D/state" ] || mkdir "$D/state" 2>/dev/null || true
    write_file "$D/state/$1" "$text"$'\n' && state_is "$1" "$2" "$3"
}

note_expect_env() {
    if [ -n "${GRM_EXPECT_SCHEMA:-}${GRM_EXPECT_REVISION:-}${GRM_EXPECT_DATASETS:-}" ]; then
        say "note: GRM_EXPECT_* is not used here; each pair is checked against the schema, revision and dataset count its approval records"
    fi
}

# Take the root lock on LOCK_FD, once per script. Every process a script forks
# for the long term closes LOCK_FD, so none holds the lock after the script.
lock_root() {
    ((LOCKED)) && return 0
    open_fd LOCK_FD "$D/deploy.lock" || die "$D/deploy.lock is a link, directory or special file" 2
    flock -n "$LOCK_FD" || die "another deployment holds $D/deploy.lock" 2
    LOCKED=1
}

# PIDs with a TCP socket on local port $1 (all belong to the listener).
listener_pids() {
    fuser -n tcp "$1" 2>/dev/null | xargs || true
}
owns_port() { [ "$(listener_pids "$1")" = "$2" ]; }

# service_of PID: set SVC_LABEL and SVC_PUB when PID is a grm serve of a
# release of this root, running its executable from the root and serving one
# of its publications; else return 1 with the reason in SVC_ERR. A server
# whose executable was deleted or replaced fails here; rollback identifies it
# by its recorded start instead (recorded_service).
service_of() {
    local pid=$1 exe cwd data
    SVC_LABEL='' SVC_PUB='' SVC_ERR=''
    [ "$(proc_comm "$pid")" = grm ] || { SVC_ERR="pid $pid on the port is not grm"; return 1; }
    exe=$(readlink "/proc/$pid/exe" 2>/dev/null) || exe=''
    cwd=$(readlink "/proc/$pid/cwd" 2>/dev/null) || cwd=''
    [ "$cwd" = "$D" ] || { SVC_ERR="pid $pid runs from ${cwd:-<unknown>}, not $D; it is not a managed release"; return 1; }
    case $exe in
        *" (deleted)") SVC_ERR="pid $pid runs a deleted or replaced executable, $exe"; return 1 ;;
        "$D"/release/*/grm) ;;
        *) SVC_ERR="pid $pid runs unmanaged executable ${exe:-<unknown>}"; return 1 ;;
    esac
    SVC_LABEL=${exe#"$D"/release/}
    SVC_LABEL=${SVC_LABEL%/grm}
    data=$(cmdline_arg "$pid" --data)
    case $data in
        "$D"/published/*) ;;
        *) SVC_ERR="pid $pid serves unmanaged data ${data:-<none>}"; return 1 ;;
    esac
    SVC_PUB=${data#"$D"/published/}
    SVC_PUB=${SVC_PUB%/}
    valid_name "$SVC_LABEL" && valid_name "$SVC_PUB" || { SVC_ERR="pid $pid serves an invalid release or publication name"; return 1; }
}

# Supervised starts (start_service). Before anything is forked, a fresh
# marker spawn/TOKEN.XXXXXX is created and flock'ed on MARK_FD, and
# serve-PORT.starting journals the random token, the marker's name and
# device:inode, the starting script's PID and start time, and the pair. The
# child is forked with the marker open, so from the fork on it, the grm it
# execs and anything that forks hold the marker's lock. The journal is then
# rewritten with the child's PID and start time; the child waits for that
# registration (spawn_gate) and for its parent to be alive before it execs,
# so a start whose script dies before registering never execs. grm runs in
# its own session (setsid), so its process group is stopped with it. If the
# script dies outright, the journal remains and names a dead owner; restore.sh
# (live port) and stage.sh (its port) then stop every holder of the marker and
# the registered child (resolve_start). The journal is removed only when the
# start committed (serve-PORT.pid records its token and a live process holds
# the marker) or when the marker's lock is free: a survivor this script cannot
# see keeps the journal, and the scripts refuse that port until it is
# resolved by hand.

# The paths a start on PORT writes are plain (see plain_target).
start_paths_ok() {
    plain_target "$D/serve-$1.starting" && plain_dir "$D/logs" && plain_target "$D/logs/serve-$1.log" && plain_dir "$D/spawn"
}

journal_text() {  # journal_text PORT LABEL NAME OWNER_START
    printf 'token=%s owner=%s owner_start=%s label=%s publication=%s port=%s mark=%s mark_id=%s' \
        "$SPAWN_TOKEN" "$$" "$4" "$2" "$3" "$1" "${SPAWN_MARK##*/}" "$SPAWN_MARK_ID"
}

# read_journal PORT: parse serve-PORT.starting into J_TOKEN, J_OWNER,
# J_OWNER_START, J_LABEL, J_PUB, J_MARK, J_MARK_ID and (once registered)
# J_CHILD and J_CHILD_START, or return 1 with PF_ERR.
read_journal() {
    local key
    J_CHILD='' J_CHILD_START=''
    parse_fields "$D/serve-$1.starting" line || return 1
    for key in token owner owner_start label publication port mark mark_id; do
        [ -n "${FIELDS[$key]:-}" ] || { PF_ERR="lacks $key"; return 1; }
    done
    if [ -n "${FIELDS[child]:-}${FIELDS[child_start]:-}" ]; then
        [[ ${FIELDS[child]:-} =~ ^[0-9]+$ && ${FIELDS[child_start]:-} =~ ^[0-9]+$ ]] ||
            { PF_ERR="has a malformed child registration"; return 1; }
        [ "${#FIELDS[@]}" = 10 ] || { PF_ERR="has unknown fields"; return 1; }
    else
        [ "${#FIELDS[@]}" = 8 ] || { PF_ERR="has unknown fields"; return 1; }
    fi
    [[ ${FIELDS[token]} =~ ^[0-9a-f-]{36}$ && ${FIELDS[owner]} =~ ^[0-9]+$ && ${FIELDS[owner_start]} =~ ^[0-9]+$ ]] &&
        [ "${FIELDS[port]}" = "$1" ] || { PF_ERR="has a malformed token, owner or port"; return 1; }
    [[ ${FIELDS[mark]} =~ ^[0-9a-f-]{36}\.[A-Za-z0-9]{6,}$ && ${FIELDS[mark]} == "${FIELDS[token]}".* && ${FIELDS[mark_id]} =~ ^[0-9]+:[0-9]+$ ]] ||
        { PF_ERR="has a malformed spawn marker"; return 1; }
    J_TOKEN=${FIELDS[token]}
    J_OWNER=${FIELDS[owner]}
    J_OWNER_START=${FIELDS[owner_start]}
    J_LABEL=${FIELDS[label]}
    J_PUB=${FIELDS[publication]}
    J_MARK=${FIELDS[mark]}
    J_MARK_ID=${FIELDS[mark_id]}
    J_CHILD=${FIELDS[child]:-}
    J_CHILD_START=${FIELDS[child_start]:-}
}

# The processes of the journaled start: the marker's holders and the
# registered child while it runs.
journal_procs() {
    local out
    out=$(mark_holders "$J_MARK_ID")
    if [ -n "$J_CHILD" ] && alive "$J_CHILD" "$J_CHILD_START" && [[ " $out " != *" $J_CHILD "* ]]; then
        out="$J_CHILD${out:+ $out}"
    fi
    printf '%s\n' "$out"
}

# pending_start PORT: set PENDING=1, J_* and ORPHANS (the processes of that
# start) when serve-PORT.starting records a start whose script is gone, else
# PENDING=0. Refuses a malformed journal and a start still in progress.
pending_start() {
    local f=$D/serve-$1.starting
    PENDING=0
    ORPHANS=''
    [ -e "$f" ] || [ -L "$f" ] || return 0
    read_journal "$1" || die "$f $PF_ERR; inspect it and remove it by hand" 2
    if alive "$J_OWNER" "$J_OWNER_START"; then
        die "pid $J_OWNER is still starting $J_LABEL with $J_PUB on port $1 ($f)" 2
    fi
    PENDING=1
    ORPHANS=$(journal_procs)
}

# Refuse to act on PORT while an interrupted start there is unresolved.
no_pending_start() {
    pending_start "$1"
    ((PENDING == 0)) || die "an interrupted start of $J_LABEL with $J_PUB on port $1 left serve-$1.starting (pids: ${ORPHANS:-none}); run ${2:-restore.sh --apply} first" 2
}

# resolve_start PORT, after pending_start set PENDING=1: if serve-PORT.pid
# records the journal's token for a live PID holding the start's marker, the
# start had committed and only the journal goes. Otherwise every process of
# the start is stopped, and the journal (with the marker) is removed only once
# the marker's lock is free; else it is kept and the script fails.
resolve_start() {
    if pid_record "$1" && [ "$PID_TOKEN" = "$J_TOKEN" ] && alive "$PID_PID" "$PID_START" && holds "$PID_PID" "$J_MARK_ID"; then
        rm -f -- "$D/serve-$1.starting" || die "cannot remove $D/serve-$1.starting" 1
        say "the interrupted start of $J_LABEL with $J_PUB on port $1 had committed (pid $PID_PID); removed its journal"
        return 0
    fi
    say "stopping pids ${ORPHANS:-none} left by an interrupted start of $J_LABEL with $J_PUB on port $1"
    mark_stop "$J_MARK_ID" "$D/spawn/$J_MARK" "$J_CHILD" "$J_CHILD_START" ||
        die "processes of the interrupted start of $J_LABEL with $J_PUB on port $1 may remain: $MARK_ERR; serve-$1.starting is kept: find and stop them, then remove it and spawn/$J_MARK by hand" 1
    rm -f -- "$D/spawn/$J_MARK" "$D/serve-$1.starting" || die "cannot remove $D/serve-$1.starting" 1
    say "resolved the interrupted start of $J_LABEL with $J_PUB on port $1 (stopped: ${ORPHANS:-none})"
}

# spawn_gate PORT OWNER_START, in the forked child: wait (up to 30 s) until
# serve-PORT.starting registers this process as the child of this start,
# while the script that forked it lives. Fails, so the child exits without
# exec'ing, if the script dies or the journal goes or names another start.
spawn_gate() {
    local self=$BASHPID me i
    me=$(proc_start "$self") || return 1
    for ((i = 0; i < 300; i++)); do
        alive "$$" "$2" || return 1
        read_journal "$1" 2>/dev/null && [ "$J_TOKEN" = "$SPAWN_TOKEN" ] || return 1
        if [ -n "$J_CHILD" ]; then
            [ "$J_CHILD $J_CHILD_START" = "$self $me" ]
            return
        fi
        sleep 0.1
    done
    return 1
}

# start_service LABEL NAME [PORT [force]]: start LABEL serving NAME on PORT
# (default the live port), detached, from the root, as a supervised start, and
# set START_PID. Returns 1, having left nothing running, when the marker or
# the journal cannot be written, unless "force" (a rollback restart) starts it
# anyway: without a journal its child does not wait for a registration, and
# without a marker only this shell's $! finds it. setsid does not fork here
# (the child leads no process group), so the PID is grm's. The server must not
# inherit the deployment lock (LOCK_FD), or it would hold the lock for its
# lifetime.
start_service() {
    local R=$D/release/$1 P=$D/published/$2 port=${3:-$PORT} me logfd='' cstart
    START_PID='' SPAWN_JOURNAL=0 SPAWN_MARK='' SPAWN_MARK_ID='' SPAWN_CHILD='' SPAWN_CHILD_START=''
    SPAWN_TOKEN=$(cat /proc/sys/kernel/random/uuid) || return 1
    me=$(proc_start "$$") || return 1
    if mark_create "$D/spawn" "$SPAWN_TOKEN"; then
        SPAWN_MARK=$MARK_PATH SPAWN_MARK_ID=$MARK_ID
        if write_file "$D/serve-$port.starting" "$(journal_text "$port" "$1" "$2" "$me")"$'\n'; then SPAWN_JOURNAL=1; fi
    fi
    if ((!SPAWN_JOURNAL)); then
        if [ "${4:-}" != force ]; then
            close_fd MARK_FD
            [ -z "$SPAWN_MARK" ] || rm -f -- "$SPAWN_MARK"
            SPAWN_MARK='' SPAWN_MARK_ID=''
            return 1
        fi
        if [ -n "$SPAWN_MARK" ]; then
            say "warning: could not write serve-$port.starting; starting $1 without a journal"
        else
            say "warning: could not create a spawn marker; starting $1 without a journal, found only as this script's child"
        fi
    fi
    checkpoint journaled
    if plain_dir "$D/logs" && { [ -d "$D/logs" ] || mkdir "$D/logs" 2>/dev/null; }; then
        open_append logfd "$D/logs/serve-$port.log" || logfd=''
    fi
    (
        [ -z "$LOCK_FD" ] || exec {LOCK_FD}>&-
        if [ -n "$logfd" ]; then exec >&"$logfd" 2>&1 {logfd}>&-; else exec >/dev/null 2>&1; fi
        exec </dev/null
        child_checkpoint child-forked
        if ((SPAWN_JOURNAL)); then spawn_gate "$port" "$me" || exit 97; fi
        child_checkpoint child-preexec
        alive "$$" "$me" || exit 97
        cd "$D" || exit 97
        GRM_DEPLOY_SPAWN=$SPAWN_TOKEN exec setsid nohup "$R/grm" --data "$P" serve --app "$R/dist" --port "$port"
    ) &
    checkpoint spawned
    SPAWN_CHILD=$!
    cstart=$(proc_start "$SPAWN_CHILD") || cstart=0
    SPAWN_CHILD_START=$cstart
    close_fd MARK_FD
    close_fd logfd
    if ((SPAWN_JOURNAL)) && ! write_file "$D/serve-$port.starting" \
        "$(journal_text "$port" "$1" "$2" "$me") child=$SPAWN_CHILD child_start=$cstart"$'\n'; then
        stop_spawn && end_journal "$port" ||
            die "cannot register the child of the start of $1 on port $port, and it would not stop; serve-$port.starting kept" 1
        say "cannot register the child of the start of $1 on port $port in serve-$port.starting; it was stopped"
        return 1
    fi
    checkpoint registered
    START_PID=$SPAWN_CHILD
    checkpoint assigned
}

# Stop every process of the last start_service, found by its marker (so even
# before START_PID or SPAWN_CHILD is set) and its registered child, and remove
# the marker once its lock shows none is left. Returns 1, keeping the marker,
# if any may remain. A forced start without a marker is found only as this
# shell's background child ($!).
stop_spawn() {
    local child=${!:-} start
    close_fd MARK_FD
    if [ -n "$SPAWN_MARK_ID" ]; then
        mark_stop "$SPAWN_MARK_ID" "$SPAWN_MARK" "$SPAWN_CHILD" "$SPAWN_CHILD_START" || {
            say "processes of start $SPAWN_TOKEN may remain: $MARK_ERR"
            return 1
        }
        rm -f -- "$SPAWN_MARK"
        return 0
    fi
    [ -n "$child" ] && [ "$(ppid_of "$child")" = "$$" ] && start=$(proc_start "$child") || return 0
    stop_tree "$child" "$start"
}

# Remove serve-PORT.starting if it journals this script's last start.
end_journal() {
    read_journal "$1" 2>/dev/null && [ "$J_TOKEN" = "$SPAWN_TOKEN" ] && rm -f -- "$D/serve-$1.starting"
}

# held_mark PID TOKEN: "NAME DEVICE:INODE" of the marker spawn/TOKEN.* that PID
# holds open, if the file is still that marker.
held_mark() {
    local link target id
    [ -n "$2" ] || return 0
    for link in "/proc/$1/fd/"*; do
        target=$(readlink "$link" 2>/dev/null) || continue
        case $target in
            "$D/spawn/$2".*) ;;
            *) continue ;;
        esac
        id=$(stat -L -c %d:%i -- "$link" 2>/dev/null) || continue
        [ "$(file_id "$target" 2>/dev/null)" = "$id" ] || continue
        printf '%s %s\n' "${target##*/}" "$id"
        return 0
    done
}
# forget_mark "NAME DEVICE:INODE": remove that marker once no process holds it.
forget_mark() {
    local name=${1% *} id=${1#* }
    [ -n "$1" ] && [[ $name =~ ^[0-9a-f-]{36}\.[A-Za-z0-9]+$ ]] || return 0
    [ "$(file_id "$D/spawn/$name" 2>/dev/null)" = "$id" ] && [ -z "$(mark_holders "$id")" ] &&
        lock_free "$D/spawn/$name" && rm -f -- "$D/spawn/$name" || true
}

# serve-PORT.pid: "pid=N start=N label=L publication=P [token=T mark=M
# mark_id=D:I]", written after the pair it names is committed in
# state/active. The token is read from the server's environment and the
# marker from its open descriptors, so the record says what the process
# itself holds. It is a record, not a lock: restore.sh checks that a healthy
# service is recorded, resolve_start that a start committed, and rollback
# identifies a damaged server by it (recorded_service).
pid_line() {  # pid_line LABEL NAME PID: the record for PID, or fail if PID is gone
    local start token mark
    start=$(proc_start "$3") || return 1
    token=$(env_value "$3" GRM_DEPLOY_SPAWN)
    [[ $token =~ ^[0-9a-f-]{36}$ ]] || token=''
    mark=$(held_mark "$3" "$token")
    [ -n "$mark" ] || token=''
    printf 'pid=%s start=%s label=%s publication=%s%s' "$3" "$start" "$1" "$2" \
        "${token:+ token=$token mark=${mark% *} mark_id=${mark#* }}"
}
record_pid() {  # record_pid PORT LABEL NAME PID
    local line
    line=$(pid_line "$2" "$3" "$4") || return 1
    write_file "$D/serve-$1.pid" "$line"$'\n'
}
# pid_record PORT: parse serve-PORT.pid into PID_PID, PID_START, PID_LABEL,
# PID_PUB, PID_TOKEN, PID_MARK and PID_MARK_ID.
pid_record() {
    PID_PID='' PID_START='' PID_LABEL='' PID_PUB='' PID_TOKEN='' PID_MARK='' PID_MARK_ID=''
    parse_fields "$D/serve-$1.pid" line || return 1
    [[ ${FIELDS[pid]:-} =~ ^[0-9]+$ && ${FIELDS[start]:-} =~ ^[0-9]+$ ]] || return 1
    PID_PID=${FIELDS[pid]}
    PID_START=${FIELDS[start]}
    PID_LABEL=${FIELDS[label]:-}
    PID_PUB=${FIELDS[publication]:-}
    PID_TOKEN=${FIELDS[token]:-}
    PID_MARK=${FIELDS[mark]:-}
    PID_MARK_ID=${FIELDS[mark_id]:-}
}

# recorded_service PID: set RS_LABEL and RS_PUB when PID is identified, by
# its recorded start rather than by its executable, as the server a supervised
# start of this root launched: serve-PORT.pid names PID with its start time,
# the token in PID's environment and the marker spawn/NAME, which still
# exists and which PID holds open; PID is grm, runs from the root and serves
# the recorded publication. Else return 1 with the reason in RS_ERR. A
# process that lacks any of these is never identified this way.
recorded_service() {
    local pid=$1
    RS_LABEL='' RS_PUB='' RS_ERR=''
    pid_record "$PORT" || { RS_ERR="serve-$PORT.pid is missing or unreadable"; return 1; }
    [ "$PID_PID" = "$pid" ] || { RS_ERR="serve-$PORT.pid records pid $PID_PID"; return 1; }
    alive "$pid" "$PID_START" || { RS_ERR="pid $pid did not start at the time serve-$PORT.pid records"; return 1; }
    [[ $PID_TOKEN =~ ^[0-9a-f-]{36}$ && $PID_MARK == "$PID_TOKEN".* && $PID_MARK =~ ^[0-9a-f-]{36}\.[A-Za-z0-9]+$ && $PID_MARK_ID =~ ^[0-9]+:[0-9]+$ ]] ||
        { RS_ERR="serve-$PORT.pid records no start token and spawn marker"; return 1; }
    [ "$(env_value "$pid" GRM_DEPLOY_SPAWN)" = "$PID_TOKEN" ] || { RS_ERR="pid $pid lacks the recorded start token"; return 1; }
    [ "$(file_id "$D/spawn/$PID_MARK" 2>/dev/null)" = "$PID_MARK_ID" ] || { RS_ERR="spawn/$PID_MARK is missing or replaced"; return 1; }
    holds "$pid" "$PID_MARK_ID" || { RS_ERR="pid $pid does not hold spawn/$PID_MARK"; return 1; }
    [ "$(proc_comm "$pid")" = grm ] || { RS_ERR="pid $pid is not grm"; return 1; }
    [ "$(readlink "/proc/$pid/cwd" 2>/dev/null)" = "$D" ] || { RS_ERR="pid $pid does not run from $D"; return 1; }
    valid_name "$PID_LABEL" && valid_name "$PID_PUB" || { RS_ERR="serve-$PORT.pid records an invalid pair"; return 1; }
    [ "$(cmdline_arg "$pid" --data)" = "$D/published/$PID_PUB" ] || { RS_ERR="pid $pid does not serve published/$PID_PUB"; return 1; }
    RS_LABEL=$PID_LABEL RS_PUB=$PID_PUB
}

wait_ready() {
    local _
    for _ in $(seq 1 40); do
        curl -sf -o /dev/null --max-time 2 "http://127.0.0.1:$1/data/index.json" && return 0
        sleep 0.25
    done
    return 1
}

# Byte and index checks of LABEL serving NAME on a port, against the schema,
# revision and dataset count in APPROVED[LABEL--NAME].
health() {
    local s r n
    read -r s r n <<<"${APPROVED[$1--$2]:?no approval loaded for $1 with $2}"
    GRM_EXPECT_SCHEMA=$s GRM_EXPECT_REVISION=$r GRM_EXPECT_DATASETS=$n \
        bash "$HERE/smoke.sh" "$3" "$D/published/$2/report" "$D/release/$1/dist"
}

# PID, and no other process, listens on the live port and serves LABEL with
# NAME. The health output goes to deploy.log when it is a plain file.
serving() {
    local fd='' rc=0
    open_append fd "$D/deploy.log" || fd=''
    if [ -n "$fd" ]; then
        wait_ready "$PORT" && owns_port "$PORT" "$3" && health "$1" "$2" "$PORT" >&"$fd" 2>&1 || rc=1
    else
        wait_ready "$PORT" && owns_port "$PORT" "$3" && health "$1" "$2" "$PORT" >/dev/null 2>&1 || rc=1
    fi
    close_fd fd
    return "$rc"
}
