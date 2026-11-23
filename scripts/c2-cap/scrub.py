"""Genome ids stay on the box (docs/cap-arms.md, Files).

  python3 -I scrub.py <release species json> <cohorts dir> < raw log > log
  python3 -I scrub.py --scan <release species json> <cohorts dir> <dir> [<dir> ...]

The ids are the release's list of the species and every cohort list
(<cohorts dir>/*.txt, when the folder exists). The first form copies a log
with each of those ids written as [genome id], as scripts/c3-inputs/run_species.sh
scrubs its logs. The second counts the ids named in every file under the
folders, prints the count and exits 1 when it is not 0.
"""

import json
import os
import re
import sys
from pathlib import Path

TOKEN = re.compile(rb"(?<![\w.])\d+\.\d+(?!\d)")


def ids(release_json, cohorts):
    out = set()
    if os.path.isfile(release_json):
        with open(release_json, encoding="utf-8") as f:
            out |= {g["id"] for g in json.load(f)["genomes"]}
    if os.path.isdir(cohorts):
        for p in sorted(Path(cohorts).glob("*.txt")):
            out |= {line.strip() for line in p.read_text(encoding="utf-8").splitlines() if line.strip()}
    return {i.encode() for i in out}


def scrub(stream_in, stream_out, known):
    for line in stream_in:
        stream_out.write(TOKEN.sub(lambda m: b"[genome id]" if m.group(0) in known else m.group(0), line))


def scan(dirs, known):
    hit = 0
    for top in dirs:
        for d, _, files in os.walk(top):
            for name in files:
                with open(os.path.join(d, name), "rb") as f:
                    hit += sum(1 for m in TOKEN.finditer(f.read()) if m.group(0) in known)
    return hit


def main(argv):
    if len(argv) >= 4 and argv[1] == "--scan":
        hit = scan(argv[4:], ids(argv[2], argv[3]))
        print(f"genome ids in {', '.join(argv[4:])}: {hit}")
        return 1 if hit else 0
    if len(argv) != 3:
        raise SystemExit(__doc__)
    scrub(sys.stdin.buffer, sys.stdout.buffer, ids(argv[1], argv[2]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
