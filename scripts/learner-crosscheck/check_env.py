"""check_env.py

Checks that the Python it runs under has exactly the packages the three
requirement files beside it pin, at the pinned versions, and prints them. It
exits 1 if a pinned package is missing or at another version. It needs no
library.

It also checks the requirement files themselves: every requirement is
`name==version` with one sha256 hash.
"""

import pathlib
import re
import sys
from importlib import metadata

HERE = pathlib.Path(__file__).resolve().parent
FILES = ("requirements-build.txt", "requirements.txt", "requirements-scm.txt")
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
    for line in wrong:
        print(line, file=sys.stderr)
    sys.exit(1 if wrong else 0)


if __name__ == "__main__":
    main()
