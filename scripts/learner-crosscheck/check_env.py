"""check_env.py

Checks that the Python it runs under has exactly the packages the three
requirement files beside it pin, at the pinned versions, and prints them. It
exits 1 if a pinned package is missing or at another version. It needs no
library.

It also checks the requirement files themselves: every requirement is
`name==version` with one sha256 hash.

Then it prints what the environment built or brought that no requirement file
pins: its pip, the pyscm-ml wheel built from the pinned archive into
ENV/built (README.md), and the compiler named in that wheel's extension
module. It exits 1 if that wheel is missing, or if the installed pyscm-ml is
not that wheel's files.
"""

import hashlib
import pathlib
import re
import struct
import sys
import zipfile
from importlib import metadata

HERE = pathlib.Path(__file__).resolve().parent
FILES = ("requirements-build.txt", "requirements.txt", "requirements-scm.txt")
BUILT = "pyscm-ml"  # the one pinned package built here, from its source archive
REQUIREMENT = re.compile(r"^([A-Za-z0-9_.-]+)==([A-Za-z0-9_.]+) \\\n    --hash=sha256:([0-9a-f]{64})$", re.M)


def pins(directory=HERE):
    """{file: [(name, version, sha256)]}, or a ValueError if a line is no pinned requirement."""
    out = {}
    for name in FILES:
        text = (pathlib.Path(directory) / name).read_text(encoding="utf-8")
        body = "\n".join(line for line in text.splitlines() if not line.startswith("#"))
        found = REQUIREMENT.findall(body)
        if REQUIREMENT.sub("", body).strip() or not found:
            raise ValueError(f"{name}: a line is no `name==version` with one sha256 hash")
        out[name] = found
    return out


def comments(data):
    """The strings of a 64-bit little-endian ELF file's .comment section, where compilers write their name."""
    if data[:6] != b"\x7fELF\x02\x01":
        raise ValueError("not a 64-bit little-endian ELF file")
    (offset,) = struct.unpack_from("<Q", data, 0x28)
    size, count, names_at = struct.unpack_from("<HHH", data, 0x3A)

    def section(index):
        name, _, _, _, start, length = struct.unpack_from("<IIQQQQ", data, offset + index * size)
        return name, data[start : start + length]

    names = section(names_at)[1]
    for index in range(count):
        name, body = section(index)
        if names[name : names.index(b"\0", name)] == b".comment":
            return sorted({text.decode("ascii") for text in body.split(b"\0") if text})
    return []


def built(prefix, distribution):
    """(rows, problems): the wheel of BUILT in PREFIX/built, its sha256, and the compiler of its extension modules.

    `distribution` is the installed package (importlib.metadata). Each file of the wheel outside its
    .dist-info folder must be the installed file of the same path, byte for byte.
    """
    name = BUILT.replace("-", "_")
    wheels = sorted((pathlib.Path(prefix) / "built").glob(f"{name}-{distribution.version}-*.whl"))
    if len(wheels) != 1:
        return [], [f"{BUILT}: {len(wheels)} wheels of {distribution.version} in {pathlib.Path(prefix) / 'built'}; README.md says how to build one"]
    wheel, problems, compilers = wheels[0], [], set()
    with zipfile.ZipFile(wheel) as archive:
        members = [m for m in archive.namelist() if not m.endswith("/") and ".dist-info/" not in m]
        for member in members:
            body = archive.read(member)
            installed = pathlib.Path(distribution.locate_file(member))
            if not installed.is_file() or installed.read_bytes() != body:
                problems.append(f"{BUILT}: the installed {member} is not the one in {wheel.name}")
            if member.endswith(".so"):
                compilers.update(comments(body))
    if not compilers:
        problems.append(f"{BUILT}: no extension module of {wheel.name} names its compiler")
    digest = hashlib.sha256(wheel.read_bytes()).hexdigest()
    rows = [
        f"{BUILT} wheel\tbuilt from the pinned archive\t{wheel.name} sha256 {digest}",
        f"{BUILT} compiler\tbuilt from the pinned archive\t{'; '.join(sorted(compilers))}",
    ]
    return rows, problems


def main():
    wrong = []
    for name, entries in pins().items():
        for package, version, _ in entries:
            try:
                have = metadata.version(package)
            except metadata.PackageNotFoundError:
                have = None
            print(f"{package}\t{version}\t{have}")
            if have != version:
                wrong.append(f"{package}: {name} pins {version}, installed is {have}")
    try:
        pip = metadata.version("pip")
    except metadata.PackageNotFoundError:
        pip = None
    print(f"pip\tnot pinned\t{pip}")
    try:
        rows, problems = built(sys.prefix, metadata.distribution(BUILT))
    except metadata.PackageNotFoundError:
        rows, problems = [], []  # reported above as a pinned package that is missing
    for row in rows:
        print(row)
    wrong += problems
    for line in wrong:
        print(line, file=sys.stderr)
    sys.exit(1 if wrong else 0)


if __name__ == "__main__":
    main()
