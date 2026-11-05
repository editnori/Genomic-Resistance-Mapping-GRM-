#!/usr/bin/env python3
"""Writes an audit.py group TSV from one field of BV-BRC genome metadata
(bvbrc/genomes.jsonl as `grm fetch` saves it), such as mlst, cgmlst_hc5,
cgmlst_hc20 or cgmlst_hc100. Genomes with the field empty or absent are left
out, so the audit counts them in rows_without_group. Each group id is
prefixed with the field name. Reads METADATA only and writes only OUT.tsv.

usage: python3 -I groups_from_metadata.py OUT.tsv METADATA FIELD
"""
import json
import sys


def groups(path, field):
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            g = json.loads(line)
            value = g.get(field)
            if value in (None, "") or isinstance(value, (list, dict)):
                continue
            out[g["genome_id"]] = f"{field}:{value}"
    return out


def main(argv):
    if len(argv) != 4:
        sys.exit(__doc__.rsplit("usage:", 1)[1].strip())
    out, path, field = argv[1:]
    rows = groups(path, field)
    with open(out, "w", encoding="utf-8") as f:
        f.write("genome_id\tgroup_id\n")
        for gid in sorted(rows):
            f.write(f"{gid}\t{rows[gid]}\n")
    print(f"{len(rows)} genomes with {field}")


if __name__ == "__main__":
    main(sys.argv)
