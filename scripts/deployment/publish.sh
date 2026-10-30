#!/usr/bin/env bash
# Copy a verified report tree into a new, read-only publication.
# usage: publish.sh SOURCE_REPORT_DIR SOURCE_MANIFEST NAME [--apply]
# SOURCE_MANIFEST lists every report file as "SHA256  ./path" (SOURCE.sha256 in
# a verified backup such as data/backups/revision7-published-20261025). The
# source must match it before anything is copied, its index must match
# GRM_EXPECT_*, and published/NAME must not exist yet.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 3 ] || die "usage: publish.sh SOURCE_REPORT_DIR SOURCE_MANIFEST NAME [--apply]" 2
need python3 sha256sum
src=$(cd -- "$1" && pwd)
manifest=$(cd -- "$(dirname -- "$2")" && pwd)/$(basename -- "$2")
name=$3
check_name "$name"
P=$D/published/$name

verify_tree "$src" "$manifest"
check_report_index "$src" "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS"
[ ! -e "$P" ] && [ ! -e "$P.partial" ] || die "$P already exists; publish under a new name" 2
say "source $src matches $manifest ($(wc -l <"$manifest") files)"
if ! ((APPLY)); then
    would "copy $src to $P/report, verify it and make it read-only"
    exit 0
fi

mkdir -p "$D/published"
mkdir "$P.partial"
cp -a "$src" "$P.partial/report"
cp "$manifest" "$P.partial/SOURCE.sha256"
verify_tree "$P.partial/report" "$P.partial/SOURCE.sha256"
chmod -R a-w "$P.partial/report" "$P.partial/SOURCE.sha256"
mv -T "$P.partial" "$P"
verify_publication "$name" "$EXPECT_SCHEMA" "$EXPECT_REVISION" "$EXPECT_DATASETS"
say "published $name files=$(wc -l <"$P/SOURCE.sha256")"
