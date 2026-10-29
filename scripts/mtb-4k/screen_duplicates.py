#!/usr/bin/env python3
"""Duplicate-isolate screen across the prefix and added MTB cohorts (stdlib only).

usage: python3 -I screen_duplicates.py --prefix-root ROOT --added-root ROOT
           [--fasta-manifest fasta-manifest.tsv] [--genomes-jsonl genomes.jsonl]
           [--near-jaccard 0.998] --out OUT.json

Run after both cohorts' matrices exist (counted.json holds each genome's
kset digest and bottom-1000 sketch) and before any training. It reports,
within each cohort and across the two:

- exact duplicates: the same genome ID, the same FASTA SHA-256 (manifest),
  the same kset digest (identical k-mer content), or the same non-empty
  BioSample or assembly accession (BV-BRC genome records);
- near duplicates: bottom-sketch Jaccard at or above --near-jaccard
  (0.998, about 30 SNP-equivalents per million bases at k = 31). Candidate
  pairs come from ten global value ranges of the pooled sketches: two
  sketches whose union bottom-1000 differs in fewer than ten values agree
  exactly in at least one range, so every pair above about 0.99 is a
  candidate; each candidate's Jaccard is then computed exactly.

`cross_cohort_exact` > 0 is a stop for a parent decision before training
(the plan's gate); near duplicates are reported for the lineage split,
which groups them, and for the sensitivity analysis.
"""

import argparse
import csv
import hashlib
import json
import sys
from itertools import combinations
from pathlib import Path

SLUG = "mycobacterium_tuberculosis"
SKETCH = 1000
RANGES = 10


def load(root: Path, cohort: str):
    sp = root / "species" / SLUG
    genomes = [x.strip() for x in (sp / "genomes.txt").read_text().splitlines() if x.strip()]
    counted = json.loads((sp / "counted.json").read_text())
    if counted.get("pending"):
        raise SystemExit(f"{root}: counted.json has pending genomes; finish the matrix first")
    by_id = {g["genome_id"]: g for g in counted["genomes"]}
    if list(by_id) != genomes:
        raise SystemExit(f"{root}: counted.json does not list genomes.txt in order")
    return [{"cohort": cohort, "id": g, "digest": by_id[g]["kset_digest"], "sketch": by_id[g]["sketch"],
             "n_kmers": by_id[g]["n_kmers"]} for g in genomes]


def mash_jaccard(a, b, s=SKETCH):
    """Jaccard estimate from two sorted bottom-s sketches: of the s smallest
    values in their union, the share both hold."""
    i = j = seen = shared = 0
    while seen < s and (i < len(a) or j < len(b)):
        if j >= len(b) or (i < len(a) and a[i] < b[j]):
            i += 1
        elif i >= len(a) or b[j] < a[i]:
            j += 1
        else:
            shared += 1
            i += 1
            j += 1
        seen += 1
    return shared / seen if seen else 1.0


def exact_groups(genomes, key):
    groups = {}
    for g in genomes:
        k = key(g)
        if k:
            groups.setdefault(k, []).append(g)
    return [v for v in groups.values() if len(v) > 1]


def screen(genomes, fasta_sha, records, near):
    exact = []
    kinds = {
        "genome_id": lambda g: g["id"],
        "fasta_sha256": lambda g: fasta_sha.get(g["id"]),
        "kset_digest": lambda g: g["digest"],
        "biosample_accession": lambda g: (records.get(g["id"]) or {}).get("biosample_accession"),
        "assembly_accession": lambda g: (records.get(g["id"]) or {}).get("assembly_accession"),
    }
    for kind, key in kinds.items():
        for group in exact_groups(genomes, key):
            exact.append({"kind": kind, "value": key(group[0]),
                          "genomes": [{"cohort": g["cohort"], "id": g["id"]} for g in group]})
    pooled = sorted(v for g in genomes for v in g["sketch"])
    bounds = [pooled[len(pooled) * r // RANGES] for r in range(1, RANGES)] if pooled else []
    buckets = {}
    for idx, g in enumerate(genomes):
        parts = [[] for _ in range(RANGES)]
        r = 0
        for v in g["sketch"]:
            while r < len(bounds) and v >= bounds[r]:
                r += 1
            parts[r].append(v)
        for r, part in enumerate(parts):
            if part:
                buckets.setdefault((r, tuple(part)), []).append(idx)
    candidates = set()
    for members in buckets.values():
        if len(members) > 1:
            candidates.update(combinations(sorted(members), 2))
    pairs = []
    for i, j in sorted(candidates):
        jac = mash_jaccard(genomes[i]["sketch"], genomes[j]["sketch"])
        if jac >= near:
            pairs.append({"a": {"cohort": genomes[i]["cohort"], "id": genomes[i]["id"]},
                          "b": {"cohort": genomes[j]["cohort"], "id": genomes[j]["id"]}, "jaccard": jac})
    cross = lambda members: len({m["cohort"] for m in members}) > 1
    return {
        "exact": exact,
        "near": pairs,
        "candidate_pairs_compared": len(candidates),
        "cross_cohort_exact": sum(1 for e in exact if cross(e["genomes"])),
        "within_cohort_exact": sum(1 for e in exact if not cross(e["genomes"])),
        "cross_cohort_near": sum(1 for p in pairs if p["a"]["cohort"] != p["b"]["cohort"]),
        "within_cohort_near": sum(1 for p in pairs if p["a"]["cohort"] == p["b"]["cohort"]),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--prefix-root", required=True, type=Path)
    ap.add_argument("--added-root", required=True, type=Path)
    ap.add_argument("--fasta-manifest", type=Path)
    ap.add_argument("--genomes-jsonl", type=Path)
    ap.add_argument("--near-jaccard", type=float, default=0.998)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, KeyError) as exc:  # a missing or unreadable input is a stop, not a traceback
        print(f"screen_duplicates: STOP: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


def run(args):
    genomes = load(args.prefix_root, "prefix") + load(args.added_root, "added")
    ids = {g["id"] for g in genomes}
    fasta_sha = {}
    if args.fasta_manifest:
        with args.fasta_manifest.open(newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                if row["genome_id"] in ids:
                    fasta_sha[row["genome_id"]] = row["sha256"]
    records = {}
    if args.genomes_jsonl:
        with args.genomes_jsonl.open() as f:
            for line in f:
                if '"genome_id"' not in line:
                    continue
                r = json.loads(line)
                if r.get("genome_id") in ids:
                    records[r["genome_id"]] = r
    out = screen(genomes, fasta_sha, records, args.near_jaccard)
    out["inputs"] = {
        "prefix_genomes": sum(1 for g in genomes if g["cohort"] == "prefix"),
        "added_genomes": sum(1 for g in genomes if g["cohort"] == "added"),
        "fasta_manifest": str(args.fasta_manifest) if args.fasta_manifest else None,
        "fasta_sha256_known": len(fasta_sha),
        "genome_records_known": len(records),
        "near_jaccard": args.near_jaccard,
    }
    with args.out.open("x") as f:
        json.dump(out, f, indent=1)
        f.write("\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("exact", "near")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
