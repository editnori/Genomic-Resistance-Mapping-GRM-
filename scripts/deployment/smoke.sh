#!/usr/bin/env bash
# Check that a running grm serve returns the exact bytes of an app build and a
# report tree, and that the served index has the expected schema, evaluation
# revision and dataset count. From revision 8 on, every dataset the served
# index lists must also record that revision, in its index row and in its
# served report. It only sends GET requests to 127.0.0.1 and writes only a
# temporary directory (mktemp), removed on exit.
# usage: smoke.sh PORT REPORT_DIR DIST_DIR
#   GRM_EXPECT_SCHEMA (default 2), GRM_EXPECT_REVISION (7), GRM_EXPECT_DATASETS
#   (237, revision 7's count; required for any other revision); each a whole
#   number of at most 18 ASCII digits. Without GRM_EXPECT_DATASETS it says that
#   the count it checks is that default.
set -euo pipefail
source "$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/lib.sh"
port=$1
report=$2
dist=$3
schema=${GRM_EXPECT_SCHEMA:-2}
revision=${GRM_EXPECT_REVISION:-7}
count=${GRM_EXPECT_DATASETS:-237}
[ "$revision" = 7 ] || [ -n "${GRM_EXPECT_DATASETS:-}" ] || { echo "FAIL GRM_EXPECT_REVISION=$revision needs GRM_EXPECT_DATASETS, that revision's report count: the default, 237, is revision 7's"; exit 2; }
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# Digits only, as common.sh requires of GRM_EXPECT_*. The index check reads
# these with Python's int(), which also takes "10\r" (a CRLF env file), "1_0"
# and full-width digits, but bash's -ge below does not, so such a revision
# skipped the per-dataset check and passed. The class lists the digits, so no
# locale widens it. At most 18 of them: from 2^63 on, -ge errors, which counts
# as false, so such a revision skipped the check the same way.
for pair in "GRM_EXPECT_SCHEMA=$schema" "GRM_EXPECT_REVISION=$revision" "GRM_EXPECT_DATASETS=$count"; do
    [[ ${pair#*=} =~ ^[0123456789]{1,18}$ ]] ||
        { echo "FAIL ${pair%%=*} must be a whole number in ASCII digits, at most 18 of them, not $(printf %q "${pair#*=}")"; exit 2; }
done
# The count's source is check_expect's (common.sh): GRM_EXPECT_DATASETS, or
# revision 7's 237. Say so when it is the default.
[ -n "${GRM_EXPECT_DATASETS:-}" ] || echo "note: GRM_EXPECT_DATASETS is unset, so the count checked is the default, 237, revision 7's"

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
# From revision 8 on, as check_report_index (common.sh) requires of a
# publication: each dataset the served index lists records the revision as its
# evaluation.revision, an integer, in its index row and in the report served
# for it. grm report stamps the index with its own revision, so the index alone
# can pass reports of an earlier revision.
if [ "$revision" -ge 8 ]; then
    cp "$tmp/body" "$tmp/index.json"
    python3 -I - "$tmp/index.json" >"$tmp/ids" <<'PY' || { echo "FAIL served index.json dataset ids"; exit 1; }
import json, re, sys
for d in json.load(open(sys.argv[1]))["datasets"]:
    i = d.get("id") if isinstance(d, dict) else None
    if not isinstance(i, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", i):
        sys.exit(f"dataset id {i!r} is not a plain name")
    print(i)
PY
    mkdir "$tmp/datasets"
    while read -r id; do
        curl -q --noproxy '*' -sS --fail --max-time 20 "http://127.0.0.1:$port/data/datasets/$id.json" -o "$tmp/datasets/$id.json" || {
            echo "FAIL /data/datasets/$id.json unreachable"
            exit 1
        }
    done <"$tmp/ids"
    python3 -I - "$tmp/index.json" "$tmp/datasets" "$revision" <<'PY' || { echo "FAIL served dataset revisions"; exit 1; }
import json, os, sys
index = json.load(open(sys.argv[1]))
root, revision = sys.argv[2], int(sys.argv[3])
stale = {}
audits = [("with index", d["id"], d) for d in index["datasets"]]
for d in index["datasets"]:
    try:
        with open(os.path.join(root, d["id"] + ".json")) as f:
            audits.append(("whose served report has", d["id"], json.load(f)))
    except (OSError, ValueError):
        stale.setdefault("datasets whose served report is not readable JSON", []).append(d["id"])
for where, name, record in audits:
    audit = record.get("evaluation") if isinstance(record, dict) else None
    found = audit.get("revision") if isinstance(audit, dict) else None
    if type(found) is not int or found != revision:
        got = "missing" if found is None else repr(found)
        stale.setdefault(f"datasets {where} evaluation.revision {got}, expected {revision}", []).append(name)
for what, names in sorted(stale.items()):
    print(f"{what}: {', '.join(sorted(names))}", file=sys.stderr)
sys.exit(1 if stale else 0)
PY
fi
echo "smoke-pass port=$port schema=$schema revision=$revision datasets=$count"
