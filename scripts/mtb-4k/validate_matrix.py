#!/usr/bin/env python3
"""Independent check of a saved GRM species matrix (stdlib only).

usage: python3 -I validate_matrix.py --matrix DIR --counted counted.json
           --genomes genomes.txt --out receipt.json [--k 31] [--min-support 2]
           [--fasta-dir DIR --recount N [--recount-seed S] [--recount-ids a,b]]
           [--require-recount] [--superset MATRIX_DIR]

It rereads the five matrix files without GRM code and checks:

1. Binding to counted.json and genomes.txt: genome order, k, min_support,
   kset digests, no pending or failed genomes.
2. Structure, as `matrix::assemble` writes it: sizes; every pattern has at
   least one k-mer, its k-mers strictly ascending; patterns in order of
   their smallest k-mer (first appearance); every k-mer inside the 2k-bit
   space and in exactly one pattern; pattern rows distinct; support equal
   to the row's popcount and inside [min_support, n - min_support]; tail
   bits clear; kmers_distinct = kmers_kept + kmers_rare + kmers_core.
3. Content of every column against counted.json's bottom sketches: a
   genome's sketch is its kset's smallest SKETCH mixed values, so for every
   kept k-mer at or below the sketch's largest value, the genome's bit must
   be set exactly when the k-mer is in the sketch (the whole kset when the
   sketch is shorter than SKETCH).
4. Omissions in the region every sketch covers: at or below the smallest
   sketch threshold, every genome's membership is known from its sketch, so
   each k-mer's support there is exact. The matrix must keep exactly the
   k-mers whose support is in [min_support, n - min_support]: a kept k-mer
   left out, or moved to the rare or core count, or a matrix built at
   another min_support, fails here when the region holds such a k-mer.
5. Optionally (--fasta-dir with --recount or --recount-ids), an independent
   recount: each chosen genome's FASTA is counted here at k with GRM's
   canonical k-mers and mixer (dna.rs), and must give counted.json's k-mer
   count, FNV-1a kset digest and sketch; then its column must hold exactly
   the kept k-mers it contains. Every recounted k-mer the matrix does not
   keep must be rare or core: one that at least min_support known genomes
   hold and at least min_support known genomes lack (recounted genomes,
   plus sketches where they reach) cannot be either, so it is a missing
   kept k-mer. When every genome is recounted, support is exact for every
   k-mer: the kept set, kmers_rare, kmers_core and kmers_distinct must all
   equal the recount's. --require-recount makes the recount mandatory (the
   plan's added, full and matched2k cohorts).
6. Optionally (--superset), containment in a matrix built over a superset
   of these genomes at the same k and min_support (the plan's full cohort):
   every genome here must be in it with the same kset digest. A k-mer kept
   here is kept there, and a k-mer kept there whose support restricted to
   these genomes is in [min_support, n - min_support] is kept here; so this
   matrix must equal the superset's rows restricted to these columns and
   filtered at min_support, k-mer for k-mer and bit for bit. A k-mer kept in
   either matrix and left out of the other fails, wherever it lies.

Any failure exits non-zero before the receipt is written.
"""

import argparse
import bisect
import gzip
import hashlib
import json
import random
import sys
from array import array
from collections import Counter
from pathlib import Path

FILES = ["matrix.json", "patterns.bin", "support.bin", "kmer_offsets.bin", "kmers.bin"]
SKETCH = 1000
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
U64 = (1 << 64) - 1


class Invalid(Exception):
    pass


def check(cond, message):
    if not cond:
        raise Invalid(message)


def integers(path: Path, code: str, width: int) -> array:
    values = array(code)
    check(values.itemsize == width, f"array code {code} is not {width} bytes here")
    data = path.read_bytes()
    check(len(data) % width == 0, f"{path}: partial integer")
    values.frombytes(data)
    if sys.byteorder != "little":
        values.byteswap()
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class KmerSpace:
    """dna.rs KmerSpace: canonical packed k-mers and the bijective mixer."""

    def __init__(self, k: int):
        check(5 <= k <= 32, f"k = {k} is outside 5..=32")
        self.k = k
        bits = 2 * k
        self.mask = U64 if bits == 64 else (1 << bits) - 1

        def scale(s):
            return max(s * bits // 64, 1)

        self.shifts = (scale(30), scale(27), scale(31))
        self.mul = ((0xBF58476D1CE4E5B9 & self.mask) | 1, (0x94D049BB133111EB & self.mask) | 1)

    def mix(self, x: int) -> int:
        s1, s2, s3 = self.shifts
        m0, m1 = self.mul
        mask = self.mask
        x ^= x >> s1
        x = (x * m0) & mask
        x ^= x >> s2
        x = (x * m1) & mask
        return x ^ (x >> s3)

    def kset(self, seqs) -> list:
        """Sorted distinct mixed canonical k-mers of the records `seqs`."""
        k, mask, mix = self.k, self.mask, self.mix
        shift = 2 * (k - 1)
        table = bytearray([4] * 256)  # dna.rs CODE: 4 for anything but ACGT
        for code, base in enumerate(b"ACGT"):
            table[base] = table[base + 32] = code
        table = bytes(table)
        out = set()
        add = out.add
        for seq in seqs:
            codes = seq.translate(table)
            fwd = rev = 0
            filled = 0
            for c in codes:
                if c > 3:
                    filled = 0
                    continue
                fwd = ((fwd << 2) | c) & mask
                rev = (rev >> 2) | ((3 - c) << shift)
                filled += 1
                if filled >= k:
                    add(mix(fwd if fwd < rev else rev))
        return sorted(out)


def fasta_records(path: Path):
    """fasta.rs for_each_record for FASTA: records may span lines; bytes
    other than ACGT (either case) break k-mer windows."""
    raw = path.read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    records, current = [], None
    for line in raw.split(b"\n"):
        line = line.rstrip(b"\r")
        if line.startswith(b">"):
            if current is not None:
                records.append(b"".join(current))
            current = []
        elif line:
            check(current is not None, f"{path}: sequence data before any header")
            current.append(line)
    if current is not None:
        records.append(b"".join(current))
    return records


def fnv_digest(values) -> str:
    h = FNV_OFFSET
    for byte in array("Q", values).tobytes() if sys.byteorder == "little" else b"".join(
        v.to_bytes(8, "little") for v in values
    ):
        h = ((h ^ byte) * FNV_PRIME) & U64
    return f"{h:016x}"


def fasta_path(fasta_dir: Path, genome_id: str) -> Path:
    # bvbrc::fasta_path; fasta_records reads it gzipped or plain.
    path = fasta_dir / f"{genome_id}.fna.gz"
    check(path.exists(), f"no FASTA for {genome_id} in {fasta_dir}")
    return path


def load_matrix(directory: Path):
    """(summary, words per row, patterns, support, offsets, kmers) of a saved matrix."""
    summary = json.loads((directory / "matrix.json").read_text())
    words = (len(summary["genomes"]) + 63) // 64
    return (summary, words, integers(directory / "patterns.bin", "Q", 8), integers(directory / "support.bin", "I", 4),
            integers(directory / "kmer_offsets.bin", "Q", 8), integers(directory / "kmers.bin", "Q", 8))


def row_int(patterns, words, p) -> int:
    return int.from_bytes(patterns[p * words:(p + 1) * words].tobytes(), "little")


def restriction(columns, n_super):
    """A function from a superset row (as an int) to the row over `columns`
    (superset column indices, in this matrix's order), as an int."""
    m = len(columns)
    start = columns[0] if columns else 0
    if columns == list(range(start, start + m)):
        mask = (1 << m) - 1
        return lambda row: (row >> start) & mask
    # Byte tables: for superset byte j, the selected bits it holds, already
    # shifted to their place in the restricted row.
    nbytes = (n_super + 7) // 8
    place = {c: i for i, c in enumerate(columns)}
    tables = []
    for j in range(nbytes):
        bits = [(b, place[8 * j + b]) for b in range(8) if 8 * j + b in place]
        if not bits:
            tables.append(None)
            continue
        table = [0] * 256
        for v in range(256):
            x = 0
            for b, at in bits:
                if (v >> b) & 1:
                    x |= 1 << at
            table[v] = x
        tables.append((j, table))
    tables = [t for t in tables if t is not None]

    def restrict(row):
        data = row.to_bytes(nbytes, "little")
        return sum(table[data[j]] for j, table in tables)
    return restrict


def check_superset(summary, words, patterns, offsets, kmers, superset_dir: Path, min_support: int, k: int):
    """Section 6 of the module docstring. Returns the receipt entry."""
    sup, sup_words, sup_patterns, sup_support, sup_offsets, sup_kmers = load_matrix(superset_dir)
    check(sup["k"] == k, f"superset k {sup['k']} differs from {k}")
    check(sup["min_support"] == min_support, f"superset min_support {sup['min_support']} differs from {min_support}")
    index = {g: i for i, g in enumerate(sup["genomes"])}
    columns = []
    for g, digest in zip(summary["genomes"], summary["kset_digests"]):
        check(g in index, f"{g} is not in the superset matrix")
        check(sup["kset_digests"][index[g]] == digest, f"{g}: kset digest differs from the superset's")
        columns.append(index[g])
    n = len(columns)
    lo, hi = min_support, n - min_support
    restrict = restriction(columns, len(sup["genomes"]))
    pattern_of = {}
    for p in range(len(offsets) - 1):
        for x in kmers[offsets[p]:offsets[p + 1]]:
            pattern_of[x] = p
    rows = {}
    matched = 0
    for p in range(len(sup_offsets) - 1):
        r = restrict(row_int(sup_patterns, sup_words, p))
        c = r.bit_count()
        kept = lo <= c <= hi
        for x in sup_kmers[sup_offsets[p]:sup_offsets[p + 1]]:
            q = pattern_of.get(x)
            if not kept:
                check(q is None, f"k-mer {x} is kept here, but its support over these genomes in the superset is {c}, outside the filter")
                continue
            check(q is not None, f"k-mer {x} has support {c} over these genomes in the superset but is missing here")
            if q not in rows:
                rows[q] = row_int(patterns, words, q)
            check(rows[q] == r, f"k-mer {x}: pattern {q} differs from the superset's row restricted to these genomes")
            matched += 1
    check(matched == len(kmers), f"{len(kmers) - matched} k-mers kept here are not kept in the superset")
    return {"superset": str(superset_dir), "superset_genomes": len(sup["genomes"]), "kmers_matched": matched,
            "superset_sha256": {name: sha256(superset_dir / name) for name in FILES}}


def validate(args) -> dict:
    summary = json.loads((args.matrix / "matrix.json").read_text())
    counted = json.loads(args.counted.read_text())
    expected = [line.strip() for line in args.genomes.read_text().splitlines() if line.strip()]
    n = len(expected)
    k = args.k
    check(0 < n <= 4096, f"{n} genomes is outside 1..=4096")
    check(len(set(expected)) == n, "genomes.txt repeats a genome")
    check(summary["genomes"] == expected, "matrix genomes differ from genomes.txt (order included)")
    check(summary["k"] == counted["k"] == k, f"k differs: matrix {summary['k']}, counted {counted['k']}, expected {k}")
    check(summary["min_support"] == args.min_support, f"min_support {summary['min_support']} != {args.min_support}")
    check(not counted.get("pending", []), "counted.json has pending genomes")
    check(not counted["failed"], "counted.json has failed genomes")
    genomes = counted["genomes"]
    check([g["genome_id"] for g in genomes] == expected, "counted.json genomes differ from genomes.txt")
    check(summary["kset_digests"] == [g["kset_digest"] for g in genomes], "kset digests differ from counted.json")

    _, words, patterns, support, offsets, kmers = load_matrix(args.matrix)
    n_patterns = summary["n_patterns"]
    check(len(patterns) == n_patterns * words, "patterns.bin size does not match n_patterns")
    check(len(support) == n_patterns, "support.bin size does not match n_patterns")
    check(len(offsets) == n_patterns + 1, "kmer_offsets.bin size does not match n_patterns")
    check(offsets[0] == 0 and offsets[-1] == len(kmers) == summary["kmers_kept"], "offsets or kmers_kept do not match kmers.bin")
    check(
        summary["kmers_distinct"] == summary["kmers_kept"] + summary["kmers_rare"] + summary["kmers_core"],
        "kmers_distinct is not kmers_kept + kmers_rare + kmers_core",
    )

    space = KmerSpace(k)
    limit = space.mask
    tail_mask = 0 if n % 64 == 0 else U64 & ~((1 << (n % 64)) - 1)
    lo, hi = args.min_support, n - args.min_support

    # Per-genome sketch thresholds: kept k-mers at or below them are checked.
    sketches = []
    for g in genomes:
        s = g.get("sketch")
        check(isinstance(s, list), f"counted.json has no sketch for {g['genome_id']}")
        check(len(s) == min(g["n_kmers"], SKETCH), f"{g['genome_id']}: sketch length {len(s)} for {g['n_kmers']} k-mers")
        check(all(0 <= v <= limit for v in s), f"{g['genome_id']}: sketch value outside the 2k-bit space")
        check(all(a < b for a, b in zip(s, s[1:])), f"{g['genome_id']}: sketch not strictly ascending")
        sketches.append((s[-1] if len(s) == SKETCH else limit, set(s)))
    threshold = max(t for t, _ in sketches)
    exact_limit = min(t for t, _ in sketches)

    seen_rows = {}
    low = []  # (k-mer, pattern) for k-mers at or below every sketch threshold
    previous_first = -1
    for p in range(n_patterns):
        row = patterns[p * words:(p + 1) * words]
        if tail_mask:
            check(row[-1] & tail_mask == 0, f"tail bits set in pattern {p}")
        observed = sum(w.bit_count() for w in row)
        check(observed == support[p], f"support mismatch in pattern {p}")
        check(lo <= observed <= hi, f"support outside the filter in pattern {p}: {observed}")
        key = row.tobytes()
        other = seen_rows.setdefault(hashlib.blake2b(key, digest_size=16).digest(), p)
        check(other == p or patterns[other * words:(other + 1) * words].tobytes() != key, f"patterns {other} and {p} have the same row")
        start, end = offsets[p], offsets[p + 1]
        check(start < end, f"pattern {p} has no k-mers")
        segment = kmers[start:end]
        check(segment[0] > previous_first, f"pattern {p} is out of first-appearance order")
        previous_first = segment[0]
        check(all(a < b for a, b in zip(segment, segment[1:])), f"unordered or repeated k-mers in pattern {p}")
        check(segment[-1] <= limit, f"pattern {p} holds a k-mer outside the 2k-bit space")
        for x in segment:
            if x > threshold:
                break
            low.append((x, p))
    ordered = sorted(kmers)
    check(all(a < b for a, b in zip(ordered, ordered[1:])), "a k-mer is in more than one pattern")
    del ordered

    def bit(p, col):
        return (patterns[p * words + col // 64] >> (col % 64)) & 1

    low.sort()
    compared = 0
    for col, (t, members) in enumerate(sketches):
        for x, p in low:
            if x > t:
                break
            compared += 1
            check(
                bit(p, col) == (x in members),
                f"column {col} ({expected[col]}): k-mer {x} is {'set' if bit(p, col) else 'clear'} in pattern {p}, but the genome's sketch says otherwise",
            )

    # 4. The region every sketch covers: exact support, so an exact kept set.
    sketch_count = Counter(x for g in genomes for x in g["sketch"])
    region = Counter({x: c for x, c in sketch_count.items() if x <= exact_limit})
    want_kept = {x for x, c in region.items() if lo <= c <= hi}
    have_kept = {x for x, _ in low if x <= exact_limit}
    missing = sorted(want_kept - have_kept)
    if missing:
        raise Invalid(f"k-mer {missing[0]} has support {region[missing[0]]} (exact: every sketch covers it) but is missing from the matrix")
    check(have_kept <= want_kept, "a kept k-mer below every sketch threshold is in no sketch at the support the matrix gives it")
    exact_region = {"limit": exact_limit, "kmers": len(region), "kept": len(want_kept),
                    "rare": sum(c < lo for c in region.values()), "core": sum(c > hi for c in region.values())}

    recount = []
    ids = list(args.recount_ids)
    if args.recount:
        pool = [g for g in expected if g not in ids]
        ids += random.Random(args.recount_seed).sample(pool, min(args.recount, len(pool)))
    if args.require_recount:
        check(args.fasta_dir is not None and ids, "--require-recount needs --fasta-dir and genomes to recount")
    omissions = None
    if ids:
        check(args.fasta_dir is not None, "--recount needs --fasta-dir")
        check(len(set(ids)) == len(ids), "a genome is named twice for the recount")
        index = {g: i for i, g in enumerate(expected)}
        kept_set = set(kmers)
        held = Counter()  # recounted genomes holding each k-mer the matrix does not keep
        thresholds = []  # sketch thresholds of the recounted genomes
        absent_sets = []
        for gid in ids:
            check(gid in index, f"recount genome {gid} is not in the matrix")
            col = index[gid]
            record = genomes[col]
            values = space.kset(fasta_records(fasta_path(args.fasta_dir, gid)))
            check(len(values) == record["n_kmers"], f"{gid}: recounted {len(values)} k-mers, counted.json has {record['n_kmers']}")
            digest = fnv_digest(values)
            check(digest == record["kset_digest"], f"{gid}: recounted kset digest {digest}, counted.json has {record['kset_digest']}")
            check(values[:SKETCH] == record["sketch"], f"{gid}: recounted sketch differs from counted.json")
            members = set(values)
            present = 0
            for p in range(n_patterns):
                b = bit(p, col)
                for x in kmers[offsets[p]:offsets[p + 1]]:
                    check((x in members) == b, f"{gid}: k-mer {x} of pattern {p} disagrees with the recount")
                present += b * (offsets[p + 1] - offsets[p])
            absent = members - kept_set
            held.update(absent)
            absent_sets.append(absent)
            thresholds.append(sketches[col][0])
            recount.append({"genome": gid, "column": col, "n_kmers": len(values), "kset_digest": digest,
                            "kept_kmers_present": present, "not_kept": len(absent)})
        # Known holders and non-holders of each recounted k-mer the matrix
        # leaves out: the recounted genomes, plus every other genome whose
        # sketch reaches it (x <= its threshold).
        all_thresholds = sorted(t for t, _ in sketches)
        recounted_thresholds = sorted(thresholds)
        everyone = len(ids) == n
        for x, h in held.items():
            a = len(ids) - h
            if x <= threshold and not everyone:
                # Genomes outside the recount whose sketch reaches x: a
                # recounted genome's sketch holds x exactly when it does.
                reach = len(all_thresholds) - bisect.bisect_left(all_thresholds, x)
                reach_recounted = len(recounted_thresholds) - bisect.bisect_left(recounted_thresholds, x)
                holders = sketch_count.get(x, 0) - sum(1 for absent, t in zip(absent_sets, thresholds) if t >= x and x in absent)
                h += holders
                a += (reach - reach_recounted) - holders
            check(h < lo or a < args.min_support,
                  f"k-mer {x} is held by at least {h} genomes and absent from at least {a}, so it is kept at "
                  f"min_support {args.min_support}, but the matrix leaves it out")
        omissions = {"recounted_kmers_not_kept": len(held), "exact": everyone}
        if everyone:
            # Every genome recounted: support is exact for every k-mer, so the
            # summary's rare, core and distinct counts are too.
            rare = sum(1 for c in held.values() if c < lo)
            core = sum(1 for c in held.values() if c > hi)
            check(rare == summary["kmers_rare"], f"kmers_rare {summary['kmers_rare']} differs from the full recount's {rare}")
            check(core == summary["kmers_core"], f"kmers_core {summary['kmers_core']} differs from the full recount's {core}")
            check(len(held) + len(kmers) == summary["kmers_distinct"],
                  f"kmers_distinct {summary['kmers_distinct']} differs from the full recount's {len(held) + len(kmers)}")

    superset = None
    if args.superset:
        superset = check_superset(summary, words, patterns, offsets, kmers, args.superset, args.min_support, k)

    return {
        "matrix": str(args.matrix),
        "genomes": n,
        "k": k,
        "min_support": args.min_support,
        "n_patterns": n_patterns,
        "kmers_kept": len(kmers),
        "checks": {
            "binding": "genome order, k, min_support, kset digests, no pending or failed genomes",
            "structure": "non-empty ascending patterns in first-appearance order, k-mers in the 2k-bit space and in one pattern, distinct rows, support = popcount within the filter, tail bits clear, distinct = kept + rare + core",
            "sketch_content": {"columns": n, "kmer_bits_compared": compared},
            "exact_region": exact_region,
            "recount": recount,
            "recount_omissions": omissions,
            "superset": superset,
        },
        "sha256": {name: sha256(args.matrix / name) for name in FILES},
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--counted", required=True, type=Path)
    parser.add_argument("--genomes", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--min-support", type=int, default=2)
    parser.add_argument("--fasta-dir", type=Path)
    parser.add_argument("--recount", type=int, default=0)
    parser.add_argument("--recount-seed", type=int, default=20261028)
    parser.add_argument("--recount-ids", type=lambda s: [x for x in s.split(",") if x], default=[])
    parser.add_argument("--require-recount", action="store_true")
    parser.add_argument("--superset", type=Path, help="matrix directory built over a superset of these genomes")
    args = parser.parse_args()
    try:
        receipt = validate(args)
    except (Invalid, KeyError, ValueError, OSError) as err:
        print(f"validate_matrix: INVALID: {err}", file=sys.stderr)
        return 1
    with args.out.open("xb") as handle:
        handle.write((json.dumps(receipt, indent=2) + "\n").encode())
    print(json.dumps(receipt, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
