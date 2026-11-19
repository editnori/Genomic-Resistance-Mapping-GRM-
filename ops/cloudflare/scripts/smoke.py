"""Checks a served GRM site against the manifests of the pair it should serve.

usage: python3 -I smoke.py BASE_URL APP_MANIFEST REPORT_MANIFEST SCHEMA REVISION DATASETS [--sweep]

GETs /, the app's index-*.js, /data/index.json, /data/corpus.json, the first
dataset report, /data/SKIPS.json when the publication has one, a client
route (/species, which must come back as index.html) and /data/research.json
(which must be 404, as the app expects). Each body must hash to its manifest
entry, and the served index must have the schema, revision and dataset count
given. --sweep also GETs every report file the manifest lists. Prints one
line per failure and a summary; exits 1 on any failure. Sends no credential.
"""
import hashlib, json, re, sys, urllib.error, urllib.request


def manifest(path):
    out = {}
    for line in open(path, encoding="ascii"):
        m = re.fullmatch(r"([0-9a-f]{64}) [ *](?:\./)?(.+)\n", line)
        out[m.group(2)] = m.group(1)
    return out


def fetch(url):
    request = urllib.request.Request(url, headers={"User-Agent": "grm-site-smoke/1", "Cache-Control": "no-cache"})
    try:
        with urllib.request.urlopen(request, timeout=120) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()


def main():
    args = [a for a in sys.argv[1:] if a != "--sweep"]
    sweep = "--sweep" in sys.argv[1:]
    base, app_path, report_path = args[0].rstrip("/"), args[1], args[2]
    schema, revision, datasets = map(int, args[3:6])
    app, report = manifest(app_path), manifest(report_path)
    failures, passed = [], 0

    def expect(path, want_status, want_sha=None):
        nonlocal passed
        status, body = fetch(base + path)
        sha = hashlib.sha256(body).hexdigest()
        if status != want_status or (want_sha and sha != want_sha):
            failures.append(f"{path}: status {status} sha {sha[:12]}, expected {want_status} {want_sha and want_sha[:12]}")
        else:
            passed += 1
        return body

    expect("/", 200, app["index.html"])
    expect("/species", 200, app["index.html"])
    scripts = sorted(p for p in app if re.fullmatch(r"assets/index-[^/]+\.js", p))
    if not scripts:
        failures.append("app manifest lists no assets/index-*.js")
    for p in scripts:
        expect("/" + p, 200, app[p])
    index = json.loads(expect("/data/index.json", 200, report["index.json"]))
    ids = [d["id"] for d in index.get("datasets", [])]
    got = (index.get("schema"), index.get("evaluation_revision"), len(ids), len(set(ids)))
    if got != (schema, revision, datasets, datasets):
        failures.append(f"index schema/revision/datasets/distinct {got}, expected {(schema, revision, datasets, datasets)}")
    expect("/data/corpus.json", 200, report["corpus.json"])
    if ids:
        expect(f"/data/datasets/{ids[0]}.json", 200, report[f"datasets/{ids[0]}.json"])
    if "SKIPS.json" in report:
        expect("/data/SKIPS.json", 200, report["SKIPS.json"])
    expect("/data/research.json", 404)
    species = sum(1 for p in report if p.startswith("species/"))
    if sweep:
        for p in sorted(report):
            expect("/data/" + p, 200, report[p])
    for f in failures:
        print("smoke-fail", f)
    print(f"smoke {'pass' if not failures else 'FAIL'} url={base} schema={got[0]} revision={got[1]} datasets={got[2]} "
          f"species={species} skips_sha256={report.get('SKIPS.json', '-')} requests_ok={passed} failures={len(failures)}")
    sys.exit(1 if failures else 0)


main()
