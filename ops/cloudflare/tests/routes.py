"""Compares two servers of the same app and report tree route by route.

usage: python3 -I routes.py GRM_URL WORKER_URL DIST_DIR REPORT_DIR OUT.json

GRM_URL is `grm serve`, WORKER_URL the Worker (wrangler dev), both serving
DIST_DIR and REPORT_DIR. Every built file, every report file, the app's
client routes, missing paths, traversal attempts, queries, HEAD, other
methods and Range forms go to both, with the path sent byte for byte. For
each, the status, the headers below and the body's SHA-256 must agree. The
Worker's ETag, and headers neither server controls (Date, Server,
Connection), are not compared. Paths the Workers runtime rewrites before the
Worker sees them (`..`, `.` and `%2e%2e` segments resolved, a `#fragment`
dropped) are listed apart, as known differences, with what each server
answered. So are answers that differ only in a Content-Length the Worker
sends and grm serve leaves out (tiny_http chunks larger bodies, so a HEAD of
one carries no length). Exits 1 on any other difference.
"""
import hashlib, http.client, json, os, sys, urllib.parse

HEADERS = ["content-type", "cache-control", "accept-ranges", "content-range", "content-length", "allow"]
# The runtime parses the request URL (WHATWG), which resolves dot segments
# and drops a fragment, so these never reach the Worker as sent. grm serve
# answers them as literal paths. Browsers make the same rewrite before sending.
NORMALISED = {"/data/../x", "/data/%2e%2e/x", "/data/%2E%2E/index.json", "/data/datasets/../index.json",
              "/data/./index.json", "/assets/../index.html", "/data/index.json#frag"}


def request(base, method, path, headers=()):
    u = urllib.parse.urlsplit(base)
    c = http.client.HTTPConnection(u.hostname, u.port, timeout=120)
    c.putrequest(method, path, skip_accept_encoding=True)
    for name, value in headers:
        c.putheader(name, value)
    c.endheaders()
    r = c.getresponse()
    body = r.read()
    got = {h: r.getheader(h) for h in HEADERS}
    if method != "HEAD" and got["content-length"] is None:
        got["content-length"] = str(len(body))
    c.close()
    return {"status": r.status, **got, "body_len": len(body), "sha256": hashlib.sha256(body).hexdigest()}


def files(root):
    out = []
    for d, _, names in os.walk(root):
        for n in names:
            out.append(os.path.relpath(os.path.join(d, n), root).replace(os.sep, "/"))
    return sorted(out)


def main():
    grm, worker, dist, report, out = sys.argv[1:6]
    index = json.load(open(os.path.join(report, "index.json")))
    ids = [d["id"] for d in index["datasets"]]
    species = sorted(p[len("species/"):-5] for p in files(report) if p.startswith("species/"))
    cases = [("GET", "/" + p, ()) for p in files(dist)]
    cases += [("GET", "/data/" + p, ()) for p in files(report)]
    cases += [("GET", p, ()) for p in [
        "/", "/index.html", "/datasets", "/species", "/methods", "/tests", "/overview", "/present",
        f"/datasets?dataset={ids[0]}", "/species/", "/a/b/c/d", "/data", "/data?x=1",
    ]]
    cases += [("GET", "/species/" + s, ()) for s in species]
    cases += [("GET", p, ()) for p in [
        "/data/research.json", "/data/missing.json", "/data/", "/data/datasets/", "/data/datasets",
        "/data/datasets/missing.json", "/assets/missing.js", "/assets/", "/assets",
        "/data/index.json?x=1", "//etc/passwd", "/data//index.json",
        "/data/datasets//" + ids[0] + ".json", "/%2e%2e/etc/passwd", "/data/..%2findex.json",
        "/data/%2e%2e%2findex.json", "/data/index%2ejson", "/data/index.json/", "/favicon.ico",
    ] + sorted(NORMALISED)]
    for p in ["/", "/data/index.json", f"/data/datasets/{ids[0]}.json", "/species", "/data/research.json", "/assets/missing.js"]:
        cases.append(("HEAD", p, ()))
    for m in ["POST", "PUT", "DELETE", "OPTIONS", "PATCH"]:
        cases += [(m, "/", ()), (m, "/data/index.json", ())]
    ranges = ["bytes=0-0", "bytes=2-5", "bytes=-3", "bytes=7-", "bytes=0-", "bytes=99999999-", "bytes=5-2",
              "bytes=0-1,4-5", "items=0-1", "bytes=abc", "bytes=-0", "bytes=-", "bytes=+2-+5", "bytes= 2-5",
              "bytes=0-99999999999", "bytes=-99999999999", "bytes=18446744073709551616-"]
    for path in ["/data/index.json", "/", "/species", "/assets/missing.js", "/data/research.json"]:
        cases += [("GET", path, (("Range", r),)) for r in ranges]
        cases.append(("HEAD", path, (("Range", "bytes=2-5"),)))
    cases.append(("GET", "/data/index.json", (("Range", "bytes=0-1"), ("Range", "bytes=2-3"))))
    results, same, differ, known, length_only = [], 0, [], [], []
    for method, path, headers in cases:
        a, b = request(grm, method, path, headers), request(worker, method, path, headers)
        row = {"method": method, "path": path, "headers": list(headers), "grm": a, "worker": b}
        results.append(row)
        if a == b:
            same += 1
        elif path in NORMALISED:
            known.append(row)
        elif a["content-length"] is None and {k for k in a if a[k] != b[k]} == {"content-length"}:
            length_only.append(row)
        else:
            differ.append(row)
    json.dump({"cases": len(cases), "same": same, "differ": differ, "known": known, "length_only": length_only,
               "results": results},
              open(out, "w"), indent=1)
    for row in differ:
        print("DIFF", row["method"], row["path"], row["headers"], {k: (row["grm"][k], row["worker"][k])
              for k in row["grm"] if row["grm"][k] != row["worker"][k]})
    for row in known:
        print("KNOWN", row["method"], row["path"], "grm", row["grm"]["status"], "worker", row["worker"]["status"],
              "worker body = index.html" if row["worker"]["sha256"] == hashlib.sha256(open(os.path.join(dist, "index.html"), "rb").read()).hexdigest() else "")
    for row in length_only:
        print("LENGTH", row["method"], row["path"], row["headers"], "worker Content-Length", row["worker"]["content-length"])
    print(f"routes cases={len(cases)} same={same} differ={len(differ)} known={len(known)} length_only={len(length_only)}")
    sys.exit(1 if differ else 0)


main()
