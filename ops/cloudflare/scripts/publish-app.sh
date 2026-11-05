#!/usr/bin/env bash
# Upload a built web app (web/dist) to the bucket as a new, sealed release.
# usage: publish-app.sh DIST_DIR LABEL COMMIT [--apply]
# DIST_DIR holds only directories and regular files and has index.html;
# COMMIT is the full commit id it was built from. app/LABEL must not exist
# yet: the claim goes first, then dist/ and COMMIT, then every object is read
# back and hashed, and MANIFEST.sha256 goes last as the seal.
source "$(dirname -- "${BASH_SOURCE[0]}")/site-lib.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 3 ] || die "usage: publish-app.sh DIST_DIR LABEL COMMIT [--apply]" 2
dist=$(cd -- "$1" && pwd)
label=$2 commit=$3
check_name "$label"
[[ $commit =~ ^[0-9a-f]{40}$ ]] || die "COMMIT must be a full commit id" 2
[ -f "$dist/index.html" ] && [ ! -L "$dist/index.html" ] || die "$dist has no index.html" 2
odd=$(cd "$dist" && find . ! -type f ! -type d -print -quit)
[ -z "$odd" ] || die "$dist holds a link or special file: ${odd#./}" 2

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
(cd "$dist" && find . -type f | LC_ALL=C sort | xargs -d '\n' sha256sum) >"$work/MANIFEST.sha256"
printf '%s\n' "$commit" >"$work/COMMIT"
files=$(wc -l <"$work/MANIFEST.sha256")
digest=$(sha256sum <"$work/MANIFEST.sha256" | cut -d' ' -f1)
say "app $dist: $files files, manifest $digest, commit $commit"
A=app/$label
if exists "$A/CLAIM" || exists "$A/MANIFEST.sha256"; then die "$BUCKET/$A already exists; publish under a new label" 2; fi
if ! ((APPLY)); then
    would "claim $BUCKET/$A, upload $files files to $A/dist/ and COMMIT, read each back, then seal $A/MANIFEST.sha256"
    exit 0
fi

printf '{"app":"%s","commit":"%s","app_manifest":"%s","files":%d,"time":"%s"}\n' \
    "$label" "$commit" "$digest" "$files" "$(date -u +%FT%TZ)" >"$work/CLAIM"
put_checked "$A/CLAIM" "$work/CLAIM"
say "claimed $BUCKET/$A"
upload_tree "$dist" "$work/MANIFEST.sha256" "$A/dist/"
put_checked "$A/COMMIT" "$work/COMMIT"
say "uploaded $files files to $BUCKET/$A/dist/"
verify_objects "$A/dist/" "$work/MANIFEST.sha256"
put_checked "$A/MANIFEST.sha256" "$work/MANIFEST.sha256"
say "published app $label files=$files app_manifest=$digest commit=$commit"
