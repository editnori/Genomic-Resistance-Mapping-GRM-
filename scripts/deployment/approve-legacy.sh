#!/usr/bin/env bash
# Migrate a hand-written GO-LABEL approval, which names no publication or report
# copy, to go/LABEL--PUBLICATION for the pair serving GRM_PORT now. Only the
# live pair qualifies: it must be this release and publication, intact, with a
# report matching GRM_EXPECT_*, and pass smoke.sh. GO-LABEL must be one line of
# unique key=value fields in printable ASCII (no NUL or other control byte
# anywhere), say stage=0, name the release's commit, and name no other
# publication. Every other pair needs stage.sh.
# usage: approve-legacy.sh LABEL PUBLICATION [--apply]
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 2 ] || die "usage: approve-legacy.sh LABEL PUBLICATION [--apply]" 2
label=$1
pub=$2
need curl python3 sha256sum fuser
check_name "$label"
check_name "$pub"
legacy=$D/GO-$label
file=$D/go/$label--$pub
[ -f "$legacy" ] && [ ! -L "$legacy" ] || die "no legacy approval $legacy" 2
[ ! -e "$file" ] && [ ! -L "$file" ] || die "$file already exists" 2
plain_dir "$D/go" || die "$D/go is not a directory" 2
verify_release "$label"
verify_publication "$pub" "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS"
commit=$(release_commit "$label")

# The whole file is checked as one line of unique key=value fields before any
# field is read, so a NUL, control byte or second line cannot hide a field.
parse_fields "$legacy" line || die "$legacy is not one line of unique key=value fields: it $PF_ERR"
[ "${FIELDS[stage]:-}" = 0 ] || die "$legacy says stage=${FIELDS[stage]:-<none>}; it must say stage=0 and commit=$commit"
[ "${FIELDS[commit]:-}" = "$commit" ] || die "$legacy names another commit: commit=${FIELDS[commit]:-<none>}"
[ -z "${FIELDS[publication]+x}" ] || [ "${FIELDS[publication]}" = "$pub" ] ||
    die "$legacy names another publication: publication=${FIELDS[publication]}"

pids=$(listener_pids "$PORT")
[ "$(wc -w <<<"$pids")" = 1 ] || die "expected one process on $PORT, got: ${pids:-none}; only the live pair can be migrated" 2
service_of "$pids" || die "$SVC_ERR" 2
[ "$SVC_LABEL $SVC_PUB" = "$label $pub" ] || die "port $PORT serves $SVC_LABEL $SVC_PUB; only the live pair can be migrated, stage $label with $pub instead" 2
APPROVED[$label--$pub]="$EXPECT_SCHEMA $EXPECT_REVISION $EXPECT_DATASETS"
health "$label" "$pub" "$PORT" >/dev/null || die "pid $pids serves $label with $pub but fails its health check" 1
say "$legacy approves $label ($commit); pid $pids serves it with $pub and passes the health check"
if ! ((APPLY)); then
    would "write $file binding this build of $label, this copy of $pub, schema $EXPECT_SCHEMA, revision $EXPECT_REVISION and $EXPECT_DATASETS datasets"
    exit 0
fi
mkdir -p "$D/go"
write_file "$file" "$(printf 'stage=0 browser=legacy commit=%s publication=%s release_manifest=%s report_manifest=%s schema=%s revision=%s datasets=%s legacy_sha256=%s' \
    "$commit" "$pub" "$(release_digest "$label")" "$(report_digest "$pub")" \
    "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS" "$(safe_sha256 "$legacy")")"$'\n' ||
    die "cannot write $file"
check_go "$label" "$pub"
say "approved $label with $pub from $legacy"
