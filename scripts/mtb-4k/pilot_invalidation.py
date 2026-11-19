#!/usr/bin/env python3
"""The pilot's three invalidation checks, on copies of the interrupted root (stdlib, Linux).

usage: python3 -I pilot_invalidation.py --copies DIR --grm GRM --width 128 --receipt OUT.json
           [--slug mycobacterium_tuberculosis] [--timeout 1800] [--grm-arg=ARG ...]

DIR holds four `cp -a` copies of the interrupted pilot cohort (each with data/ inside), taken
after pilot_interrupt.py and before the resume: width, forge, link-matrix and link-kset. Each must
still be the interrupted state (pilot_interrupt.interrupted_problems at --width). Then, each on
its own copy and never on the live root, the resume
`GRM --data COPY/data matrix --species SLUG --min-support 2 --disk-sparse --block-width W` must
be refused (--grm-arg adds an argument the interrupted run was started with to every resume; the
plan passes none, tests pass --allow-missing-catalog):

1. width: at --block-width W-64 (W+64 when W is 64; GRM's widths are multiples of 64), with
   "unfinished disk-sparse run".
2. forge: block 0's folded partial loses one k-mer from one column (a row that holds that column
   and at least one other is cleared at that column's bit), and its receipt's SHA-256 and byte
   count are recomputed, so the receipt agrees with the forged file. The resume must fail
   verify_partial with "does not hold that genome's kset".
3. link-matrix: species/SLUG/matrix is a link to a directory beside the copy that holds one
   sentinel file; link-kset: ksets/<id>.kset.part, for a block-1 genome with no kset yet, is a
   link to a sentinel file beside the copy. Each resume must be refused before any work, and
   each target must keep its bytes (the directory holds only its sentinel).

For every check, the copy's species directory (each file's SHA-256 and size, each link's target)
must be the same after the refused run as before it, so the partials stay in place. It exits 0
only when all four are refused as stated, and 1 otherwise; the receipt records each run.
"""
import argparse
import hashlib
import importlib
import json
import os
import struct
import subprocess
import sys
import time
from pathlib import Path

LANE = Path(__file__).resolve().parent
if sys.path[:1] != [str(LANE)]:  # M1: python3 -I leaves this directory off sys.path
    sys.path.insert(0, str(LANE))
pi = importlib.import_module("pilot_interrupt")
if Path(pi.__file__).resolve().parent != LANE:
    sys.exit(f"pilot_invalidation: STOP: pilot_interrupt loaded from {pi.__file__}, not from {LANE}")

COPIES = ("width", "forge", "link-matrix", "link-kset")
MAGIC = b"GRMKPRT2"
HEAD = 48
SENTINEL = b"pilot_invalidation sentinel: must keep these bytes\n"


def tree(species: Path):
    """Every entry below species, links not followed: {relative path: (kind, detail)}."""
    out = {}
    for top, dirs, files in os.walk(species):
        for name in dirs + files:
            path = Path(top) / name
            rel = str(path.relative_to(species))
            st = path.lstat()
            if path.is_symlink():
                out[rel] = ("link", os.readlink(path))
            elif path.is_dir():
                out[rel] = ("dir", None)
            elif path.is_file():
                out[rel] = ("file", (st.st_size, hashlib.sha256(path.read_bytes()).hexdigest()))
            else:
                out[rel] = ("other", st.st_mode)
    return out


def changes(before, after):
    return sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))


def varint(data, pos):
    value, shift = 0, 0
    while True:
        byte = data[pos]
        pos += 1
        value |= (byte & 0x7F) << shift
        if byte < 0x80:
            return value, pos
        shift += 7


def drop_one_kmer(data: bytearray):
    """Clear one column's bit in one non-full row that holds at least one other column; return
    (column, byte offset of the changed word). The partial format is src/matrix.rs write_partial."""
    if data[:8] != MAGIC:
        raise ValueError("not a GRMKPRT2 partial")
    bucket_bits = data[9]
    width = struct.unpack_from("<Q", data, 24)[0]
    pos = HEAD
    for _ in range(width):
        pos += 4 + struct.unpack_from("<I", data, pos)[0] + 8
    buckets = 1 << bucket_bits
    offsets = struct.unpack_from(f"<{buckets + 1}Q", data, pos)
    counts = struct.unpack_from(f"<{buckets}I", data, pos + (buckets + 1) * 8)
    start = pos + (buckets + 1) * 8 + buckets * 4
    words = (width + 63) // 64
    for bucket in range(buckets):
        p = start + offsets[bucket]
        for _ in range(counts[bucket]):
            entry, p = varint(data, p)
            if entry & 1:
                continue  # a full row: every column holds it, no row words follow
            row = [struct.unpack_from("<Q", data, p + 8 * w)[0] for w in range(words)]
            if sum(bin(w).count("1") for w in row) >= 2:
                for w, value in enumerate(row):
                    if value:
                        bit = value & -value
                        struct.pack_into("<Q", data, p + 8 * w, value & ~bit)
                        return w * 64 + bit.bit_length() - 1, p + 8 * w
            p += 8 * words
    raise ValueError("no row holds two columns; nothing to forge")


def resume(grm, data_root, slug, width, timeout, extra=()):
    cmd = [grm, "--data", str(data_root), "matrix", "--species", slug, "--min-support", "2", "--disk-sparse",
           "--block-width", str(width), *extra]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return {"cmd": cmd, "rc": p.returncode, "stderr": p.stderr[-2000:]}
    except subprocess.TimeoutExpired as exc:  # subprocess.run has killed it
        return {"cmd": cmd, "rc": None, "stderr": f"timed out after {exc.timeout} s"}


def run_dir(species: Path):
    found, _ = pi.run_state(species)
    if found is None:
        raise ValueError(f"{species} holds no single sparse run directory")
    return found


def check(copy: Path, args, mutate, width, want):
    """Mutate the copy, run the resume, and require a refusal with `want` that changes nothing."""
    species = copy / "data" / "species" / args.slug
    row = {"copy": str(copy)}
    state = pi.snapshot(species, *pi.run_state(species))
    problems = pi.interrupted_problems(state, args.width)
    if problems:
        row["error"] = "the copy is not the interrupted state: " + "; ".join(problems)
        return row
    row["mutation"] = mutate(copy, species) if mutate else None
    before = tree(species)
    row["run"] = resume(args.grm, copy / "data", args.slug, width, args.timeout, args.grm_arg)
    row["changed"] = changes(before, tree(species))[:20]
    reasons = []
    if row["run"]["rc"] in (0, None):
        reasons.append(f"the resume was not refused (rc {row['run']['rc']})")
    if want not in row["run"]["stderr"]:
        reasons.append(f"the refusal does not say {want!r}")
    if row["changed"]:
        reasons.append(f"the refused run changed the copy: {row['changed']}")
    if reasons:
        row["error"] = "; ".join(reasons)
    return row


def forge(copy, species):
    folder = run_dir(species)
    part = folder / "block-00000.partial"
    receipt_path = folder / "block-00000.partial.receipt.json"
    data = bytearray(part.read_bytes())
    column, offset = drop_one_kmer(data)
    part.write_bytes(bytes(data))
    receipt = json.loads(receipt_path.read_text())
    receipt["sha256"], receipt["bytes"] = hashlib.sha256(bytes(data)).hexdigest(), len(data)
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    return {"partial": str(part), "column": column, "word_offset": offset, "sha256": receipt["sha256"]}


def link_matrix(copy, species):
    target = copy / "outside-matrix"
    target.mkdir()
    (target / "sentinel").write_bytes(SENTINEL)
    os.symlink(target, species / "matrix", target_is_directory=True)
    return {"link": str(species / "matrix"), "target": str(target)}


def link_kset(copy, species, width):
    genomes = (species / "genomes.txt").read_text().split()
    ksets = species / "ksets"
    free = [g for g in genomes[width:2 * width] if not any(os.path.lexists(ksets / f"{g}{s}") for s in (".kset", ".kset.part", ".json", ".json.part"))]
    if not free:
        raise ValueError("every block-1 genome already has a kset or sidecar; no free name to link")
    target = copy / "outside-kset"
    target.write_bytes(SENTINEL)
    link = ksets / f"{free[-1]}.kset.part"
    os.symlink(target, link)
    return {"link": str(link), "target": str(target)}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--copies", required=True, type=Path)
    ap.add_argument("--grm", required=True)
    ap.add_argument("--width", required=True, type=int)
    ap.add_argument("--receipt", required=True, type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    ap.add_argument("--timeout", type=float, default=1800)
    ap.add_argument("--grm-arg", action="append", default=[])
    args = ap.parse_args(argv)
    with args.receipt.open("x") as out:
        receipt = {"copies": str(args.copies), "grm": args.grm, "width": args.width, "started": time.time(), "checks": {}}
        try:
            missing = [c for c in COPIES if not (args.copies / c / "data").is_dir()]
            if missing or args.width < 64 or args.width % 64:
                raise ValueError(f"need copies {list(COPIES)} with data/ under {args.copies}; missing {missing}")
            plan = [
                ("width", None, args.width - 64 if args.width > 64 else args.width + 64, "unfinished disk-sparse run"),
                ("forge", forge, args.width, "does not hold that genome's kset"),
                ("link-matrix", link_matrix, args.width, "refusing symlinked directory"),
                ("link-kset", lambda c, s: link_kset(c, s, args.width), args.width, "symlink"),
            ]
            for name, mutate, width, want in plan:
                receipt["checks"][name] = row = check(args.copies / name, args, mutate, width, want)
                if name.startswith("link") and "error" not in row:
                    target = Path(row["mutation"]["target"])
                    kept = (target / "sentinel").read_bytes() == SENTINEL and sorted(os.listdir(target)) == ["sentinel"] \
                        if target.is_dir() else target.read_bytes() == SENTINEL
                    row["target_kept"] = kept
                    if not kept:
                        row["error"] = f"the link target {target} changed"
            errors = {name: row["error"] for name, row in receipt["checks"].items() if "error" in row}
            receipt["error"] = errors or None
        except (OSError, ValueError, KeyError, struct.error) as exc:
            receipt["error"] = f"{type(exc).__name__}: {exc}"
        receipt["refused_as_stated"] = receipt["error"] is None
        receipt["finished"] = time.time()
        json.dump(receipt, out, indent=1)
        out.write("\n")
    if receipt["error"]:
        print(f"pilot_invalidation: STOP: {receipt['error']}", file=sys.stderr)
        return 1
    print(json.dumps({name: {"rc": row["run"]["rc"], "changed": row["changed"]} for name, row in receipt["checks"].items()}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
