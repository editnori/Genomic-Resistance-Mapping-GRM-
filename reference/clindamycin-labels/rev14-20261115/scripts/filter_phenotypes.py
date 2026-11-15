#!/usr/bin/env python3
"""Arm "without" (PLAN.md section 5): copy curated/phenotypes.jsonl, dropping only the clindamycin rows of the
listed genomes. Every other line is written byte for byte as read.

Run: python3 -I scripts/filter_phenotypes.py <in phenotypes.jsonl> <removal ids> <out phenotypes.jsonl>
Prints the number of lines read, dropped and written, and the number of listed ids that had no row.
"""
import json
import sys

DRUG = "clindamycin"


def main():
    src, ids_path, dst = sys.argv[1:4]
    ids = set(open(ids_path).read().split())
    seen = set()
    n_in = n_drop = n_out = 0
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        for line in fi:
            n_in += 1
            if DRUG.encode() in line:
                r = json.loads(line)
                if r["antibiotic"] == DRUG and r["genome_id"] in ids:
                    n_drop += 1
                    seen.add(r["genome_id"])
                    continue
            fo.write(line)
            n_out += 1
    print(f"read {n_in} dropped {n_drop} written {n_out} listed {len(ids)} listed-without-row {len(ids - seen)}")
    if ids - seen or n_drop != len(ids):
        sys.exit("every listed genome must have exactly one clindamycin row")


if __name__ == "__main__":
    main()
