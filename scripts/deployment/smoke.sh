#!/usr/bin/env bash
# Check that a running grm serve returns the exact bytes of an app build and a
# report tree, and that the served index has the expected schema, evaluation
# revision and dataset count. It only sends GET requests to 127.0.0.1 and
# writes only a temporary directory (mktemp), removed on exit.
# usage: smoke.sh PORT REPORT_DIR DIST_DIR
#   GRM_EXPECT_SCHEMA (default 2), GRM_EXPECT_REVISION (7), GRM_EXPECT_DATASETS (237)
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
port=$1
report=$2
dist=$3
schema=${GRM_EXPECT_SCHEMA:-2}
revision=${GRM_EXPECT_REVISION:-7}
count=${GRM_EXPECT_DATASETS:-237}
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

check() {
    local url=$1 file=$2 got want
    # status.sh runs this on trees it has not verified; never open a link or FIFO.
    want=$(safe_sha256 "$file") || { echo "FAIL $file is not a regular file"; exit 1; }
    # Never through a proxy: -q skips ~/.curlrc, --noproxy '*' the proxy variables.
    curl -q --noproxy '*' -sS --fail --max-time 20 "http://127.0.0.1:$port$url" -o "$tmp/body" || {
        echo "FAIL $url unreachable"
        exit 1
    }
    got=$(sha256sum <"$tmp/body" | cut -d' ' -f1)
    if [ "$got" != "$want" ]; then
        echo "FAIL $url got=$got want=$want"
        exit 1
    fi
    echo "ok $url $got"
}

check / "$dist/index.html"
check /data/corpus.json "$report/corpus.json"
shopt -s nullglob
datasets=("$report/datasets"/*.json)
((${#datasets[@]})) || { echo "FAIL no dataset reports in $report/datasets"; exit 1; }
check "/data/datasets/${datasets[0]##*/}" "${datasets[0]}"
assets=("$dist"/assets/index-*.js)
((${#assets[@]})) || { echo "FAIL no index JavaScript in $dist/assets"; exit 1; }
check "/assets/${assets[0]##*/}" "${assets[0]}"
check /data/index.json "$report/index.json"
python3 -I - "$tmp/body" "$schema" "$revision" "$count" <<'PY' || { echo "FAIL served index.json"; exit 1; }
import json, sys
index = json.load(open(sys.argv[1]))
schema, revision, count = map(int, sys.argv[2:])
got = (index.get("schema"), index.get("evaluation_revision"), len(index.get("datasets", [])))
if got != (schema, revision, count):
    sys.exit(f"schema/revision/datasets {got}, expected {(schema, revision, count)}")
PY
echo "smoke-pass port=$port schema=$schema revision=$revision datasets=$count"
