"""scan.py: fail when a file names a genome of the species or a BioProject
accession.

    python3 -I scan.py --ids SPECIES.json PATH [PATH ...]

SPECIES.json is the release's species file (its `genomes[].id`). Every file
under each PATH is read. A token of digits, a dot and digits that is a genome
id of the list, or any text of the form PRJ<letters><digits>, is a hit. The
output gives counts per file only, never the hit itself. Exit 1 on any hit.

Python 3.12, standard library only; run with -I.
"""

import argparse
import json
import os
import re
import sys

TOKEN = re.compile(rb"(?<![\w.])\d+\.\d+(?![\w])")
BIOPROJECT = re.compile(rb"PRJ[A-Z]{1,3}\d+")


def load_ids(path):
    with open(path, encoding="utf-8") as f:
        return {g["id"].encode() for g in json.load(f)["genomes"]}


def scan_file(path, ids):
    with open(path, "rb") as f:
        data = f.read()
    genomes = sum(1 for m in TOKEN.finditer(data) if m.group(0) in ids)
    projects = len(BIOPROJECT.findall(data))
    return genomes, projects


def files_under(paths):
    for p in paths:
        if os.path.isfile(p):
            yield p
            continue
        for root, _, names in os.walk(p):
            for n in sorted(names):
                yield os.path.join(root, n)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids", required=True)
    ap.add_argument("paths", nargs="+")
    args = ap.parse_args(argv)
    ids = load_ids(args.ids)
    files = hits = 0
    for path in files_under(args.paths):
        files += 1
        g, b = scan_file(path, ids)
        if g or b:
            hits += 1
            print(f"scan: {path}: {g} genome ids, {b} BioProject accessions")
    print(f"scan: {files} files, {hits} with a hit, against {len(ids)} genome ids")
    return 1 if hits else 0


if __name__ == "__main__":
    sys.exit(main())
