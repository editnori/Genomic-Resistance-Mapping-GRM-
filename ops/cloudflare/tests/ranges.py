"""Range requests on one large file, checked against the file on disk.

usage: python3 -I ranges.py BASE_URL URL_PATH LOCAL_FILE OUT.json

A full GET (streamed and hashed), HEAD, the first and last byte, suffixes,
open ranges, ranges across the whole file from a fixed seed, ranges past the
end (clamped, as grm serve clamps them) and unsatisfiable ones (416). Every
206 body must equal those bytes of LOCAL_FILE, with the Content-Range grm
serve sends. A Content-Length, when one is sent, must equal the length; grm
serve sends none on larger bodies (tiny_http chunks them), so the summary
counts the answers that carried one. Exits 1 on any failure.
"""
import hashlib, http.client, json, os, random, sys, time, urllib.parse


def get(base, method, path, rng=None):
    u = urllib.parse.urlsplit(base)
    c = http.client.HTTPConnection(u.hostname, u.port, timeout=600)
    c.putrequest(method, path, skip_accept_encoding=True)
    if rng:
        c.putheader("Range", rng)
    c.endheaders()
    r = c.getresponse()
    h, n = hashlib.sha256(), 0
    body = b""
    while chunk := r.read(1 << 20):
        h.update(chunk)
        n += len(chunk)
        if n <= (1 << 22):
            body += chunk
    head = {k: r.getheader(k) for k in ("content-range", "content-length", "content-type", "accept-ranges")}
    c.close()
    return r.status, head, n, h.hexdigest(), body


def main():
    base, path, local, out = sys.argv[1:5]
    size = os.path.getsize(local)
    data = open(local, "rb")
    whole = hashlib.sha256()
    while chunk := data.read(1 << 20):
        whole.update(chunk)

    def want(start, length):
        data.seek(start)
        return data.read(length)

    rows, failures = [], []

    def record(row, ok, head, **extra):
        rows.append(row | {"ok": ok, "length_header": head["content-length"] is not None})
        if not ok:
            failures.append(rows[-1] | extra | {"head": head})

    def case(rng, status, start=None, length=None):
        t = time.time()
        s, head, n, sha, body = get(base, "GET", path, rng)
        ok = s == status
        if status == 206:
            ok = ok and n == length and body == want(start, length)
            ok = ok and head["content-range"] == f"bytes {start}-{start + length - 1}/{size}"
            ok = ok and head["content-length"] in (None, str(length))
        if status == 416:
            ok = ok and head["content-range"] == f"bytes */{size}"
        record({"range": rng, "status": s, "bytes": n, "seconds": round(time.time() - t, 3)}, ok, head, want_status=status)

    t = time.time()
    s, head, n, sha, _ = get(base, "GET", path)
    ok = s == 200 and n == size and sha == whole.hexdigest() and head["content-length"] in (None, str(size))
    record({"range": None, "status": s, "bytes": n, "sha256": sha, "seconds": round(time.time() - t, 3)}, ok, head)
    s, head, n, _, _ = get(base, "HEAD", path)
    ok = s == 200 and n == 0 and head["content-length"] in (None, str(size)) and head["accept-ranges"] == "bytes"
    record({"range": "HEAD", "status": s, "bytes": n}, ok, head)
    case("bytes=0-0", 206, 0, 1)
    case(f"bytes={size - 1}-{size - 1}", 206, size - 1, 1)
    case("bytes=-1", 206, size - 1, 1)
    case("bytes=-1048576", 206, size - 1048576, 1048576)
    case(f"bytes={size - 4096}-", 206, size - 4096, 4096)
    case(f"bytes={size - 10}-{size + 1000}", 206, size - 10, 10)
    case(f"bytes={size}-", 416)
    case(f"bytes={size + 1}-{size + 5}", 416)
    case("bytes=10-5", 416)
    case("bytes=0-1,5-6", 416)
    rnd = random.Random(20261105)
    for _ in range(64):
        start = rnd.randrange(size)
        length = min(rnd.choice([1, 17, 4096, 65536, 1 << 20, 3 << 20]), size - start)
        case(f"bytes={start}-{start + length - 1}", 206, start, length)
    lengths = sum(1 for r in rows if r["length_header"])
    json.dump({"size": size, "sha256": whole.hexdigest(), "cases": len(rows), "with_content_length": lengths,
               "failures": failures, "rows": rows}, open(out, "w"), indent=1)
    print(f"ranges url={base}{path} size={size} sha256={whole.hexdigest()} cases={len(rows)} "
          f"failures={len(failures)} with_content_length={lengths}")
    sys.exit(1 if failures else 0)


main()
