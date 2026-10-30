#!/usr/bin/env bash
# Build grm and the dashboard from one commit and install them as a new,
# immutable release. Build on the host that serves it: .cargo/config.toml sets
# target-cpu=native.
# usage: build-release.sh LABEL COMMIT [--apply]
#   GRM_SOURCE_REPO   git checkout to archive COMMIT from (default: this checkout)
#   GRM_EXTRA_PATH    prepended to PATH, e.g. for cargo and node outside PATH
#   GRM_BUILD_JOBS    cargo jobs (default 8)
#   GRM_MIN_FREE_GIB  free space required on the root's filesystem (default 40)
# Writes only src-LABEL, src-LABEL.commit, target-LABEL (removed after the
# build) and release/LABEL, none of which may exist beforehand. Creating
# src-LABEL reserves the label before anything is extracted, so of two builds
# under one label only the first proceeds. The release holds grm, dist/, COMMIT
# and, when the commit has one, web/scripts/browser-check.cjs as
# browser-check.cjs, which stage.sh runs. That file is taken from the commit
# itself, not from src-LABEL, where npm ci and npm run build have run.
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/common.sh"
parse_args "$@"
set -- "${ARGS[@]}"
[ $# -eq 2 ] || die "usage: build-release.sh LABEL COMMIT [--apply]" 2
if [ -n "${GRM_EXTRA_PATH:-}" ]; then export PATH=$GRM_EXTRA_PATH:$PATH; fi
label=$1
check_name "$label"
repo=${GRM_SOURCE_REPO:-$(cd -- "$HERE/../.." && pwd)}
min_gib=${GRM_MIN_FREE_GIB:-40}
need git tar sha256sum df
commit=$(git -C "$repo" rev-parse --verify --quiet "$2^{commit}") || die "$2 is not a commit in $repo" 2
S=$D/src-$label
R=$D/release/$label
T=$D/target-$label
for path in "$S" "$S.commit" "$R" "$R.partial" "$T"; do
    [ ! -e "$path" ] || die "$path already exists; build under a new label" 2
done
[ -d "$D" ] || die "deployment root $D does not exist" 2
avail=$(df -B1 --output=avail "$D" | tail -1)
[ "$avail" -ge $((min_gib * 1024 * 1024 * 1024)) ] || die "below $min_gib GiB free under $D" 2
say "build $label from $commit ($repo)"
if ! ((APPLY)); then
    need cargo npm
    would "extract $commit to $S, build it in $T, install release $R"
    exit 0
fi

need cargo npm nice
mkdir -p "$D/release"
mkdir "$S" 2>/dev/null || die "$S already exists; another build may hold label $label" 2
for path in "$S.commit" "$R" "$R.partial" "$T"; do
    [ ! -e "$path" ] || { rmdir "$S"; die "$path already exists; build under a new label" 2; }
done
git -C "$repo" archive "$commit" | tar -x -C "$S"
printf '%s\n' "$commit" >"$S.commit"
# The browser check is meant to verify the build independently, but npm ci and
# npm run build run in src-LABEL/web, where a dependency's install script could
# rewrite web/scripts/browser-check.cjs. So its blob in the commit (a regular
# file, not a link) is noted before anything runs, and the release gets that
# blob, checked by its ID.
check_mode='' check_blob=''
read -r check_mode _ check_blob _ < <(git -C "$repo" ls-tree --full-tree "$commit" web/scripts/browser-check.cjs) || true
case $check_mode in 100644 | 100755) ;; *) check_blob='' ;; esac
export CARGO_TARGET_DIR=$T CARGO_BUILD_JOBS=${GRM_BUILD_JOBS:-8} CARGO_INCREMENTAL=0
(cd "$S" && nice -n 15 cargo build --release --locked)
(cd "$S/web" && nice -n 15 npm ci --no-audit --no-fund && nice -n 15 npm run build)
mkdir "$R.partial"
cp "$T/release/grm" "$R.partial/grm"
cp -a "$S/web/dist" "$R.partial/dist"
cp "$S.commit" "$R.partial/COMMIT"
# stage.sh runs the browser check from here, where the manifest binds it to COMMIT.
if [ -n "$check_blob" ]; then
    git -C "$repo" cat-file blob "$check_blob" >"$R.partial/browser-check.cjs"
    [ "$(git -C "$repo" hash-object --no-filters -- "$R.partial/browser-check.cjs")" = "$check_blob" ] ||
        die "the browser check copied into $R.partial is not blob $check_blob of $commit"
fi
(cd "$R.partial" && find . -type f ! -path ./MANIFEST.sha256 -print0 | LC_ALL=C sort -z | xargs -0 sha256sum >MANIFEST.sha256)
chmod -R a-w "$R.partial"
mv -T "$R.partial" "$R"
rm -rf "$T" "$S/web/node_modules" "$S/web/dist"
verify_release "$label"
say "release $label commit=$commit grm=$(sha256sum "$R/grm" | cut -c1-16)"
