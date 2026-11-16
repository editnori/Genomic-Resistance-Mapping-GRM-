"""Paths and readers shared by the lasso layout suite. It reads files of the tree and builds nothing."""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[2]
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
PAGE = ROOT / "docs" / "lasso-layout.md"
SCRIPT = ROOT / "scripts" / "lasso-layout" / "check.sh"
LASSO = ROOT / "src" / "learn" / "lasso.rs"
CARGO = ROOT / "Cargo.toml"
SUITE = "tests/lasso_layout"
JOB = "lasso-layout"
COMMAND = "bash scripts/lasso-layout/check.sh"


def text(path):
    return pathlib.Path(path).read_text(encoding="utf-8")


def jobs():
    """Each job of the workflow as (key, text), in file order."""
    body = text(WORKFLOW).split("\njobs:\n", 1)[1]
    found = []
    for chunk in re.split(r"\n  (?=[a-z][a-z0-9_-]*:\n)", "\n" + body):
        m = re.match(r"([a-z][a-z0-9_-]*):\n", chunk)
        if m:
            found.append((m.group(1), chunk))
    return found


def job(key):
    matches = [chunk for name, chunk in jobs() if name == key]
    return matches[0] if len(matches) == 1 else None


def steps(job_text):
    """A job's steps, each as its lines without comments and blank lines."""
    found = []
    inside = False
    for line in job_text.splitlines():
        if line == "    steps:":
            inside = True
            continue
        if not inside or not line.strip() or line.lstrip().startswith("#"):
            continue
        if line.startswith("      - "):
            found.append([line])
        elif found and line.startswith("        "):
            found[-1].append(line)
        else:
            inside = False
    return found


def section(page_text, header):
    """The lines of a page under a header, up to the next header of any level."""
    part = page_text.split(f"\n{header}\n", 1)[1]
    return re.split(r"\n#+ ", part, maxsplit=1)[0]
