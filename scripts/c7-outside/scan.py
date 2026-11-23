"""The genome-id scan of a C7 run (docs/outside-learners.md, The run), as C3's run ends with.

  python3 -I scan.py <species json> <dir or file> [...]

<species json> is the release's species file (species/<species>.json), whose
"genomes" list holds each genome's "id". Every file under each dir is read
as text and each token that looks like a genome id (digits, a point, digits,
not inside a longer word or number) is looked up in that list. Prints the
count of hits and the count of files read, never a token, and exits 1 when
any id is found.
"""

import json
import os
import re
import sys

TOKEN = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")


def ids_of(species_json):
    with open(species_json, encoding="utf-8") as f:
        return {g["id"] for g in json.load(f)["genomes"]}


def files(paths):
    for top in paths:
        if os.path.isfile(top):
            yield top
            continue
        for d, _, names in os.walk(top):
            for name in sorted(names):
                yield os.path.join(d, name)


def scan(ids, paths):
    """(hits, files read)."""
    hit = n = 0
    for path in files(paths):
        with open(path, errors="replace") as f:
            hit += sum(1 for m in TOKEN.finditer(f.read()) if m.group(0) in ids)
        n += 1
    return hit, n


def main(argv):
    if len(argv) < 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    ids = ids_of(argv[0])
    if not ids:
        print("the species file lists no genome", file=sys.stderr)
        return 2
    hit, n = scan(ids, argv[1:])
    print(f"genome ids: {hit} in {n} files")
    return int(hit != 0)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
