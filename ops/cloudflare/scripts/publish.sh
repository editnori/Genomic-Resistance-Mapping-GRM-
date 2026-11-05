#!/usr/bin/env bash
# Upload a verified report tree to the bucket as a new, sealed publication.
# usage: publish.sh SOURCE_REPORT_DIR SOURCE_MANIFEST NAME [--apply]
# The source must pass scripts/deployment/publish.sh's checks (the tree
# against SOURCE_MANIFEST, its index against GRM_EXPECT_*), run as a dry run
# against an empty root. GRM_EXPECT_REVISION and GRM_EXPECT_DATASETS are both
# required; GRM_EXPECT_SCHEMA defaults to 2. pub/NAME must not exist yet: the
# claim goes first, then the files, then every object is read back and
# hashed, and SOURCE.sha256 goes last as the seal.
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 3 ] || die "usage: publish.sh SOURCE_REPORT_DIR SOURCE_MANIFEST NAME [--apply]" 2
[ -n "${GRM_EXPECT_REVISION:-}" ] && [ -n "${GRM_EXPECT_DATASETS:-}" ] ||
    die "set both GRM_EXPECT_REVISION and GRM_EXPECT_DATASETS for this publication" 2
schema=${GRM_EXPECT_SCHEMA:-2}
src=$(cd -- "$1" && pwd)
manifest=$(cd -- "$(dirname -- "$2")" && pwd)/$(basename -- "$2")
name=$3
check_name "$name"

empty=$(mktemp -d)
GRM_DEPLOY_ROOT=$empty GRM_EXPECT_SCHEMA=$schema bash "$TOOLKIT/publish.sh" "$src" "$manifest" "$name" >/dev/null ||
    { rm -rf "$empty"; die "$src fails scripts/deployment/publish.sh's checks"; }
rm -rf "$empty"
files=$(wc -l <"$manifest")
digest=$(sha256sum <"$manifest" | cut -d' ' -f1)
say "source $src matches $manifest ($files files, $digest); index schema $schema, revision $GRM_EXPECT_REVISION, $GRM_EXPECT_DATASETS datasets"
P=pub/$name
if exists "$P/CLAIM" || exists "$P/SOURCE.sha256"; then die "$BUCKET/$P already exists; publish under a new name" 2; fi
if ! ((APPLY)); then
    would "claim $BUCKET/$P, upload $files files to $P/report/, read each back, then seal $P/SOURCE.sha256"
    exit 0
fi

claim=$(mktemp)
printf '{"publication":"%s","report_manifest":"%s","schema":%d,"revision":%d,"datasets":%d,"files":%d,"time":"%s"}\n' \
    "$name" "$digest" "$schema" "$GRM_EXPECT_REVISION" "$GRM_EXPECT_DATASETS" "$files" "$(date -u +%FT%TZ)" >"$claim"
put_checked "$P/CLAIM" "$claim"
rm -f "$claim"
say "claimed $BUCKET/$P"
upload_tree "$src" "$manifest" "$P/report/"
say "uploaded $files files to $BUCKET/$P/report/"
verify_objects "$P/report/" "$manifest"
put_checked "$P/SOURCE.sha256" "$manifest"
say "published $name files=$files report_manifest=$digest"
