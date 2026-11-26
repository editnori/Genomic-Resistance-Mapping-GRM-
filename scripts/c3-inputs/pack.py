"""What leaves the run's machine after plan block C3's run on its datasets of M (docs/input-arms.md, The run on the datasets of M).

  python3 -I pack.py check <OUT> <species json>
  python3 -I pack.py tar <OUT> <species json> <tar file>

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
  run_species.sh's end looks for them.

tar runs check, then writes <tar file> of the listed files and nothing
else, in name order, with no owner, group or time, and prints its sha256.
The tar file is written outside OUT or as OUT/pack.tar.
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
GENOME_ID = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
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
        if not parts or parts[0] in LEFT or ".." in parts:
            raise Stop(f"{rel}: under work/ or flags/, which stay on the machine")
        path = out / rel
        if path.is_symlink():
            raise Stop(f"{rel}: a link, which could reach into work/")
        if rel.endswith(".scores.json"):
            raise Stop(f"{rel}: a score file")
        if any(parts[i:i + 2] == ("m", "list") for i in range(len(parts) - 1)):
            raise Stop(f"{rel}: a list file of genomes")
        text = path.read_bytes().decode("utf-8", errors="replace")
        why = genome_doc(text) if rel.endswith(".json") else None
        if why:
            raise Stop(f"{rel}: {why}")
        hits = sum(1 for m in GENOME_ID.finditer(text) if m.group(0) in ids)
        if hits:
            raise Stop(f"{rel}: names {hits} genome ids of the species")
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


def main(argv):
    if len(argv) not in (3, 4) or argv[0] not in ("check", "tar") or (argv[0] == "tar") != (len(argv) == 4):
        print("\n".join(__doc__.strip().splitlines()[2:4]), file=sys.stderr)
        return 2
    out, ids = Path(argv[1]), species_ids(argv[2])
    try:
        files = check(out, ids)
    except Stop as e:
        print(f"pack: refused: {e}")
        return 1
    if argv[0] == "check":
        print(f"pack: {len(files)} files, none under work/, no score file, report, fold result or list, "
              "no genome id of the species")
        return 0
    tar_path = Path(argv[3]).resolve()
    if out.resolve() in tar_path.parents and tar_path != out.resolve() / "pack.tar":
        print("pack: the tar file is written outside OUT or as OUT/pack.tar")
        return 2
    print(f"pack: {len(files)} files in {tar_path.name}, sha256 {write_tar(out, files, tar_path)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
