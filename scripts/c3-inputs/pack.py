"""What leaves the run's machine after plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I pack.py check <OUT> <species json>
  python3 -I pack.py tar <OUT> <species json> <tar file>
  python3 -I pack.py mask <species json> < <log> > <log, each id as [genome id]>

The pack is every file under OUT but work/, flags/, run.lock, DONE, FAILED
and pack.tar: results/
(the pick files, counts per fold; the builds' records; m/floor.tsv), check/,
score/ (both readings' tables), logs/ (each written with every genome id as
[genome id]), binary.txt, expect.tsv, expect-first.tsv, lists.txt,
m-checked.txt and SHA256SUMS. check refuses it, naming the path and why
with no genome id, when it would carry:

- a path under work/, or a link, which could reach into it;
- a score file (grm learn --exp-scores, format grm-c3-scores-1, or a name
  ending .scores.json), a study fit's fold result (grm-fold-file-result-v1),
  a report (a lineage split's per-genome scores, or a list of samples), or
  a list file of genomes (work/m/list);
- any file that names a genome id of the species: <species json> is the
  release's species file, whose genomes' ids are the ones looked for, as
  run_species.sh's end looks for them;
- a path (a file or folder name) that names one of those ids;
- a compressed file, or a file that is not UTF-8 text or holds a NUL byte
  (UTF-16, say), which check cannot read.

An id is the form a genome id takes, digits, a point and digits, not inside
a longer number and not after a letter; after "_" or "." it is still read,
and so is a part of a longer run of points and digits. Each match is held
to the species file's ids, so a figure or a version that names no genome
passes. A refusal never prints an id: in a path it reads [genome id].

The pack holds both readings, as the count machine needs them. It goes
unopened from the run's machine to the one machine that takes the count
over all the species (ruling R27), and to no other; it is no page file.
Only the count tables and the page tables leave that machine.

tar runs check, then writes <tar file> of the listed files and nothing
else, in name order, with no owner, group or time, and prints its sha256.
The tar file is written outside OUT or as OUT/pack.tar.

mask copies its input to its output with each genome id of the species,
read as check reads one, written as [genome id]; any other byte is kept.
With no species file yet (before the release is unpacked) it copies the
input as it is. run_m.sh's scrub runs it on every log.
"""

import hashlib
import io
import json
import os
import re
import sys
import tarfile
from pathlib import Path

LEFT = {"work", "flags"}
# DONE and FAILED are written after the pack, and pack.tar is the pack itself.
LEFT_FILES = {"run.lock", "pack.tar", "DONE", "FAILED"}
# Each start not after a letter or digit: an id after "_" or "." is read, and matches may overlap.
GENOME_ID = re.compile(r"(?<![0-9A-Za-z])(?=(\d+\.\d+(?!\d)))")
COMPRESSED = {b"\x1f\x8b": "gzip", b"BZh": "bzip2", b"\xfd7zXZ\x00": "xz", b"\x28\xb5\x2f\xfd": "zstd",
              b"PK\x03\x04": "zip", b"PK\x05\x06": "zip", b"\x04\x22\x4d\x18": "lz4", b"7z\xbc\xaf\x27\x1c": "7z"}
COMPRESSED_NAMES = (".gz", ".tgz", ".bz2", ".xz", ".zst", ".zip", ".lz4", ".7z", ".Z")
GENOME_FORMATS = {"grm-c3-scores-1": "a score file", "grm-fold-file-result-v1": "a study fit's fold result"}


class Stop(Exception):
    pass


def species_ids(path):
    with open(path) as f:
        return {g["id"] for g in json.load(f)["genomes"]}


def listing(out):
    """Every file the pack would carry, as paths relative to OUT, in name order; links are kept to be refused."""
    out = Path(out)
    files = []
    for d, dirs, names in os.walk(out):
        rel = Path(d).relative_to(out)
        if rel == Path("."):
            dirs[:] = [x for x in dirs if x not in LEFT]
            names = [n for n in names if n not in LEFT_FILES]
        for x in list(dirs):
            if (Path(d) / x).is_symlink():
                files.append(str(rel / x))
                dirs.remove(x)
        files += [str(rel / n) for n in names]
    return sorted(os.path.normpath(f) for f in files)


def found(text, ids):
    """The genome ids of the species that the text names, each place it names one."""
    return [m.group(1) for m in GENOME_ID.finditer(text) if m.group(1) in ids]


def masked(text, ids):
    """The text with each place check reads a genome id of the species written as [genome id], and only there."""
    spans = [(m.start(), m.start() + len(m.group(1))) for m in GENOME_ID.finditer(text) if m.group(1) in ids]
    out, at = [], 0
    for start, end in spans:
        if start < at:
            # Overlaps the id before it: [genome id] already covers its start.
            start = at
        out.append(text[at:start] + ("[genome id]" if end > start else ""))
        at = max(at, end)
    return "".join(out) + text[at:]


def compressed(rel, data):
    """The kind of compression of a file, by its first bytes or its name, or None."""
    for magic, kind in COMPRESSED.items():
        if data.startswith(magic):
            return kind
    return "by its name" if rel.endswith(COMPRESSED_NAMES) else None


def genome_doc(text):
    """Why a JSON file holds genomes one by one, or None."""
    try:
        doc = json.loads(text)
    except ValueError:
        return None
    if not isinstance(doc, dict):
        return None
    if doc.get("format") in GENOME_FORMATS:
        return GENOME_FORMATS[doc["format"]]
    if isinstance(doc.get("genomes"), list) and doc["genomes"] and isinstance(doc["genomes"][0], dict) \
            and "scores" in doc["genomes"][0]:
        return "a score file"
    if isinstance(doc.get("samples"), list) and doc["samples"]:
        return "a report or a result that lists its samples"
    return None


def check(out, ids, files=None):
    """The files of the pack. Raises Stop naming the first path refused and why."""
    out = Path(out)
    files = listing(out) if files is None else files
    for rel in files:
        parts = Path(rel).parts
        name = masked(rel, ids)
        if not parts or parts[0] in LEFT or ".." in parts:
            raise Stop(f"{name}: under work/ or flags/, which stay on the machine")
        path = out / rel
        if path.is_symlink():
            raise Stop(f"{name}: a link, which could reach into work/")
        if rel.endswith(".scores.json"):
            raise Stop(f"{name}: a score file")
        if any(parts[i:i + 2] == ("m", "list") for i in range(len(parts) - 1)):
            raise Stop(f"{name}: a list file of genomes")
        in_path = len(found(rel, ids))
        if in_path:
            raise Stop(f"{name}: its path names {in_path} genome ids of the species")
        data = path.read_bytes()
        kind = compressed(rel, data)
        if kind:
            raise Stop(f"{name}: a compressed file ({kind}), which check cannot read")
        try:
            text = data.decode("utf-8")
        except UnicodeDecodeError:
            raise Stop(f"{name}: not UTF-8 text, which check cannot read") from None
        if "\x00" in text:
            raise Stop(f"{name}: holds a NUL byte (UTF-16, say), which check cannot read as text")
        why = genome_doc(text) if rel.endswith(".json") else None
        if why:
            raise Stop(f"{name}: {why}")
        hits = len(found(text, ids))
        if hits:
            raise Stop(f"{name}: names {hits} genome ids of the species")
    return files


def write_tar(out, files, tar_path):
    out = Path(out)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for rel in files:
            data = (out / rel).read_bytes()
            info = tarfile.TarInfo(rel)
            info.size, info.mtime, info.mode = len(data), 0, 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            tar.addfile(info, io.BytesIO(data))
    Path(tar_path).write_bytes(buf.getvalue())
    return hashlib.sha256(buf.getvalue()).hexdigest()


def mask_stream(ids, src, dst):
    """Each line of src to dst with each genome id of the species as [genome id]; other bytes as they are."""
    for line in src:
        text = line.decode("utf-8", errors="surrogateescape")
        dst.write(masked(text, ids).encode("utf-8", errors="surrogateescape"))


def main(argv):
    if len(argv) == 2 and argv[0] == "mask":
        ids = species_ids(argv[1]) if os.path.isfile(argv[1]) else set()
        mask_stream(ids, sys.stdin.buffer, sys.stdout.buffer)
        return 0
    if len(argv) not in (3, 4) or argv[0] not in ("check", "tar") or (argv[0] == "tar") != (len(argv) == 4):
        print("\n".join(__doc__.strip().splitlines()[2:5]), file=sys.stderr)
        return 2
    out, ids = Path(argv[1]), species_ids(argv[2])
    try:
        files = check(out, ids)
    except Stop as e:
        print(f"pack: refused: {e}")
        return 1
    if argv[0] == "check":
        print(f"pack: {len(files)} files, none under work/, no score file, report, fold result or list, "
              "no compressed or non-UTF-8 file, no genome id of the species in a path or a file")
        return 0
    tar_path = Path(argv[3]).resolve()
    if out.resolve() in tar_path.parents and tar_path != out.resolve() / "pack.tar":
        print("pack: the tar file is written outside OUT or as OUT/pack.tar")
        return 2
    print(f"pack: {len(files)} files in {tar_path.name}, sha256 {write_tar(out, files, tar_path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
