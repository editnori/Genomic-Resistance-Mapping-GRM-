"""Take score.py's three tables out of one species' score.tar by name, and nothing else (docs/cap-arms-result.md, Method).

  python3 -I take_tables.py <score.tar> <species folder>

c2-drive.sh brings each species' score/ back as score.tar, made by
`tar -cf score-<species>.tar score` and never unpacked. page_tables.py
reads score/checks.tsv, score/draws.tsv and score/datasets.tsv as files.
This step writes those three, and only those, to <species folder>/score/,
from the members named score/<table> or ./score/<table>.

It walks the archive's 512-byte headers itself and seeks past every
member's bytes without reading them, so score/score.log, and any other
member, is never read and never written to disk. It refuses:
- a compressed archive (gzip, bzip2, xz or zstd), which it would have to
  unpack whole to reach a member;
- a header it cannot parse, or an archive that ends before its end block;
- a wanted table given twice, missing, or not a plain file;
- a header of another kind (a link, a long name, a pax header), since its
  names could not be held to the three;
- a <species folder>/score/ that already holds any of the three, or a
  score.log.

Prints one line per table with its bytes and sha256, and exits 1 on a
refusal, writing nothing.
"""

import hashlib
import sys
import tarfile
from pathlib import Path

TABLES = ("checks.tsv", "draws.tsv", "datasets.tsv")
BLOCK = 512
COMPRESSED = {b"\x1f\x8b": "gzip", b"BZh": "bzip2", b"\xfd7zXZ\x00": "xz", b"\x28\xb5\x2f\xfd": "zstd"}
PLAIN = (tarfile.REGTYPE, tarfile.AREGTYPE)


class Refuse(Exception):
    pass


def table_of(name):
    """The table a member's name gives, or None: score/<table> or ./score/<table> only."""
    for prefix in ("score/", "./score/"):
        if name.startswith(prefix) and name[len(prefix):] in TABLES:
            return name[len(prefix):]
    return None


def take(f):
    """{table: bytes} from an open, seekable archive. Reads every header and the three tables' bytes only."""
    start = f.read(6)
    for magic, kind in COMPRESSED.items():
        if start.startswith(magic):
            raise Refuse(f"the archive is {kind}-compressed and would have to be unpacked whole")
    f.seek(0)
    found = {}
    while True:
        head = f.read(BLOCK)
        if len(head) < BLOCK:
            raise Refuse("the archive ends before its end block")
        if head == bytes(BLOCK):
            break
        try:
            info = tarfile.TarInfo.frombuf(head, "utf-8", "surrogateescape")
        except tarfile.HeaderError as e:
            raise Refuse(f"a header cannot be parsed: {e}") from None
        if info.type not in PLAIN + (tarfile.DIRTYPE,):
            raise Refuse(f"a member of type {info.type!r}: only plain files and folders are taken")
        data = f.tell()
        table = table_of(info.name)
        if table is not None:
            if table in found:
                raise Refuse(f"score/{table} is in the archive twice")
            if info.type not in PLAIN:
                raise Refuse(f"score/{table} is not a plain file")
            found[table] = f.read(info.size)
            if len(found[table]) != info.size:
                raise Refuse(f"score/{table} is cut short")
        f.seek(data + -(-info.size // BLOCK) * BLOCK)
    missing = [t for t in TABLES if t not in found]
    if missing:
        raise Refuse(f"the archive holds no {', '.join('score/' + t for t in missing)}")
    return found


def write(found, folder):
    score = Path(folder) / "score"
    there = [n for n in TABLES + ("score.log",) if (score / n).exists()]
    if there:
        raise Refuse(f"{score} already holds {', '.join(there)}")
    score.mkdir(parents=True, exist_ok=True)
    lines = []
    for t in TABLES:
        with open(score / t, "xb") as out:
            out.write(found[t])
        lines.append(f"score/{t}\t{len(found[t])} bytes\t{hashlib.sha256(found[t]).hexdigest()}")
    return lines


def main(argv):
    if len(argv) != 2:
        sys.stderr.write(__doc__)
        return 2
    try:
        with open(argv[0], "rb") as f:
            found = take(f)
        lines = write(found, argv[1])
    except Refuse as e:
        print(f"refused: {e}", file=sys.stderr)
        return 1
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
