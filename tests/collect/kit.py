"""What the tests of scripts/collect share: the tool, a temporary folder, and one small synthetic project.

Everything here is synthetic. No test reads a real record, a release or a genome.
"""

import hashlib
import io
import os
import pathlib
import shutil
import subprocess
import sys
import tarfile
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
TOOL = ROOT / "scripts" / "collect" / "collect.py"
BARRED = ("sealed", "held-out", "ks30", "ks50", "vumc", "clinical", "private", "ksd", "sigdrug")
EMPTY_SHA256 = "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"

# The synthetic project. The sizes in the comments are counted by hand and the tests use them as written.
FILES = {
    "reports/alpha.json": b'{"dataset": "alpha"}\n',  # 21 bytes
    "reports/beta.json": b'{"dataset": "beta"}\n',  # 20
    "notes/ledger.md": b"# Ledger\n\nrow 1\n",  # 16
    "notes/messages/m1.txt": b"a message\n",  # 10
    "notes/alpha-copy.json": b'{"dataset": "alpha"}\n',  # 21, the content of reports/alpha.json again
    "notes/empty.txt": b"",  # 0
    "inputs/genome1.fna": b">g1\nACGT\n",  # 9
    "build/target/debug/app": b"\x7fELF" + b"\0" * 96,  # 100
    "build/target/debug/deps/lib.rlib": b"r" * 50,  # 50
    "build/__pycache__/m.pyc": b"pyc",  # 3
    "env/lib/site.py": b"import sys\n",  # 11
    "scratch/run.tmp": b"tmp",  # 3
}
RULES = ("# synthetic rules\n"
         "public-result\tproj/reports/*\n"
         "public-input\tproj/inputs/*\n"
         "regenerable\tproj/env/\n"
         "regenerable\t*.tmp\n"
         "internal-record\t*\n")
# The list of the synthetic project: (path, bytes, hashed content or None, class), in the list's order.
ROWS = [
    ("proj/build/__pycache__/", 3, None, "regenerable"),
    ("proj/build/target/", 150, None, "regenerable"),  # 100 + 50, over two levels
    ("proj/env/", 11, None, "regenerable"),
    ("proj/inputs/genome1.fna", 9, FILES["inputs/genome1.fna"], "public-input"),
    ("proj/notes/alpha-copy.json", 21, FILES["notes/alpha-copy.json"], "internal-record"),
    ("proj/notes/empty.txt", 0, b"", "internal-record"),
    ("proj/notes/ledger.md", 16, FILES["notes/ledger.md"], "internal-record"),
    ("proj/notes/messages/m1.txt", 10, FILES["notes/messages/m1.txt"], "internal-record"),
    ("proj/reports/alpha.json", 21, FILES["reports/alpha.json"], "public-result"),
    ("proj/reports/beta.json", 20, FILES["reports/beta.json"], "public-result"),
    ("proj/scratch/run.tmp", 3, None, "regenerable"),
]
INTERNAL = ["proj/notes/alpha-copy.json", "proj/notes/empty.txt", "proj/notes/ledger.md", "proj/notes/messages/m1.txt"]
PUBLIC = ["proj/reports/alpha.json", "proj/reports/beta.json"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def noise(seed, n):
    """n bytes that do not compress, the same on every run."""
    return hashlib.shake_256(str(seed).encode()).digest(n)


def list_text(rows=ROWS):
    lines = ["path\tbytes\tsha256\tclass"]
    lines += [f"{p}\t{n}\t{'-' if data is None else sha(data)}\t{cls}" for p, n, data, cls in rows]
    return "\n".join(lines) + "\n"


def tmp(case):
    """A new empty folder, removed when the test ends. A random name could hold a barred word, so it is redrawn."""
    while True:
        d = os.path.realpath(tempfile.mkdtemp(prefix="collect-"))
        if not any(w in d.lower() for w in BARRED):
            break
        os.rmdir(d)
    case.addCleanup(restore_and_remove, d)
    return pathlib.Path(d)


def restore_and_remove(d):
    for base, dirs, files in os.walk(d):
        for name in dirs + files:
            p = os.path.join(base, name)
            if not os.path.islink(p):
                os.chmod(p, 0o700)
    shutil.rmtree(d)


def write(base, files):
    for rel, data in files.items():
        p = base / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(data)


def project(case):
    """(folder, the project's folder, the rules file) with the synthetic project written."""
    d = tmp(case)
    write(d / "proj", FILES)
    (d / "proj" / "empty").mkdir()
    (d / "rules.tsv").write_text(RULES, encoding="utf-8")
    return d, d / "proj", d / "rules.tsv"


def run(*args, cwd=None):
    """(exit status, stdout, stderr) of the tool."""
    p = subprocess.run([sys.executable, "-I", "-B", "-W", "error", str(TOOL), *map(str, args)],
                       capture_output=True, text=True, cwd=cwd, stdin=subprocess.DEVNULL)
    return p.returncode, p.stdout, p.stderr


def listed(case):
    """The synthetic project with its list written: (folder, project, list)."""
    d, proj, rules = project(case)
    rc, out, err = run("list", "--rules", rules, "--root", f"proj={proj}", "--out", d / "list.tsv")
    case.assertEqual((rc, err), (0, ""), out)
    return d, proj, d / "list.tsv"


def published(d, proj):
    """A sha256sum-format file naming the two reports, as a release's checksum file would."""
    text = "".join(f"{sha(FILES[f'reports/{n}'])}  {n}\n" for n in ("alpha.json", "beta.json"))
    (d / "release.sha256").write_text(text, encoding="utf-8")
    return d / "release.sha256"


def packed(case, cls="internal-record", compress="gz", *more):
    """The synthetic project packed: (folder, project, the parts.sha256 file)."""
    d, proj, lst = listed(case)
    key = ("--published", published(d, proj)) if cls == "public-result" else ()
    rc, out, err = run("pack", "--list", lst, "--class", cls, "--name", "demo", "--out", d / "out",
                       "--root", f"proj={proj}", "--compress", compress, *key, *more)
    case.assertEqual((rc, err), (0, ""), out)
    return d, proj, d / "out" / f"demo.{cls}.parts.sha256"


def parts_of(sums):
    return [sums.parent / line[66:] for line in sums.read_text(encoding="utf-8").splitlines()]


def members(part):
    """{name: content} of the regular files of one part; any other member is named with None."""
    out = {}
    with tarfile.open(part) as tar:
        for m in tar:
            out[m.name] = tar.extractfile(m).read() if m.isreg() else None
    return out


def unpack(sums, into):
    """Unpacks every part into one new folder and returns the NAME.CLASS folder."""
    into.mkdir()
    for part in parts_of(sums):
        with tarfile.open(part) as tar:
            tar.extractall(into, filter="data")
    return into / sums.name[:-len(".parts.sha256")]


def rewrite(part, change):
    """Rewrites one part through change(name, data) -> [(name, data)] and returns its new sha256."""
    old = members(part)
    mode = "w:gz" if part.name.endswith(".gz") else "w:xz" if part.name.endswith(".xz") else "w"
    with tarfile.open(part, mode, format=tarfile.PAX_FORMAT) as tar:
        for name, data in old.items():
            for new_name, new_data in change(name, data):
                info = tarfile.TarInfo(new_name)
                info.size = len(new_data)
                tar.addfile(info, io.BytesIO(new_data))
    return sha(part.read_bytes())


def resum(sums):
    """Rewrites the parts.sha256 file to the parts as they now are, as someone hiding a change would."""
    sums.write_text("".join(f"{sha(p.read_bytes())}  {p.name}\n" for p in parts_of(sums)), encoding="utf-8")
