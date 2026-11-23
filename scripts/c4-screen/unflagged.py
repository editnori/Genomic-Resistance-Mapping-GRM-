"""The two tools of unflagged.sh, the check that grm with no C4 flag writes main's bytes.

  python3 -I unflagged.py fasta <test_dna.fa> <genomes.txt> <fasta dir>
  python3 -I unflagged.py compare <dir a> <dir b> [<skip prefix> ...]

fasta: one SYNTHETIC genome per id of <genomes.txt>, gzipped as
<fasta dir>/<id>.fna.gz: the contigs of the AMRFinderPlus program's
bundled test_dna.fa, each kept or left out by a draw seeded by the id, then
two contigs of 20,000 bases of DNA drawn from the same seed. No genome's
own sequence is read. Prints counts only.

compare: every file under <dir a> against the same path under <dir b>, logs
(*.log) left out: equal bytes, or JSON equal once the run clocks
(generated_at, seconds, seconds_inner) are set to null. Prints counts only,
and for a file that differs or stands in one folder only its path with every
token of a genome id's shape written [id]. Exits 1 when a file differs or
stands in one folder only.
"""

import gzip
import hashlib
import io
import json
import os
import random
import re
import sys
from pathlib import Path

CLOCKS = {"generated_at", "seconds", "seconds_inner"}
ID = re.compile(r"\d+\.\d+")
FILLER = 20_000


def read_fasta(path):
    out, name, seq = [], None, []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(seq)))
            name, seq = line[1:].split()[0], []
        elif line.strip():
            seq.append(line.strip())
    if name is not None:
        out.append((name, "".join(seq)))
    return out


def genome(contigs, gid):
    """The gzipped FASTA bytes of one SYNTHETIC genome, the same for the same id."""
    rng = random.Random(int(hashlib.sha256(f"c4-unflagged|{gid}".encode()).hexdigest(), 16))
    parts = [s for _, s in contigs if rng.random() < 0.5]
    parts += ["".join(rng.choice("ACGT") for _ in range(FILLER)) for _ in range(2)]
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:
        for i, s in enumerate(parts):
            gz.write(f">contig{i + 1}\n".encode())
            for j in range(0, len(s), 80):
                gz.write((s[j:j + 80] + "\n").encode())
    return buf.getvalue()


def fasta(test_dna, genomes, out_dir):
    contigs = read_fasta(test_dna)
    ids = [x.strip() for x in Path(genomes).read_text().splitlines() if x.strip()]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for gid in ids:
        (out / f"{gid}.fna.gz").write_bytes(genome(contigs, gid))
    print(f"synthetic genomes: {len(ids)}, from {len(contigs)} test_dna.fa contigs")
    return 0


def unclock(v):
    if isinstance(v, dict):
        return {k: (None if k in CLOCKS else unclock(x)) for k, x in v.items()}
    if isinstance(v, list):
        return [unclock(x) for x in v]
    return v


def files(top, skip):
    out = set()
    for dirpath, _, names in os.walk(top):
        for n in names:
            rel = os.path.relpath(os.path.join(dirpath, n), top)
            # Logs carry times and stay out of the comparison.
            if not rel.endswith(".log") and not any(rel.startswith(s) for s in skip):
                out.add(rel)
    return out


def compare(a, b, skip):
    fa, fb = files(a, skip), files(b, skip)
    same_bytes = same_unclocked = 0
    differ = []
    for rel in sorted(fa & fb):
        x, y = Path(a, rel).read_bytes(), Path(b, rel).read_bytes()
        if x == y:
            same_bytes += 1
            continue
        if rel.endswith(".json"):
            try:
                if unclock(json.loads(x)) == unclock(json.loads(y)):
                    same_unclocked += 1
                    continue
            except ValueError:
                pass
        differ.append(rel)
    only = sorted(fa ^ fb)
    print(f"files: {len(fa & fb)} in both, {len(only)} in one only; equal bytes {same_bytes}, "
          f"equal but for the run clocks {same_unclocked}, differing {len(differ)}")
    for rel in differ[:20]:
        print("differs: " + ID.sub("[id]", rel))
    for rel in only[:20]:
        print("in one only: " + ID.sub("[id]", rel))
    return 0 if not differ and not only else 1


def main(argv):
    if argv[:1] == ["fasta"] and len(argv) == 4:
        return fasta(*argv[1:])
    if argv[:1] == ["compare"] and len(argv) >= 3:
        return compare(argv[1], argv[2], argv[3:])
    print(__doc__.strip().splitlines()[2], file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
