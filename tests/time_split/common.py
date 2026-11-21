"""Paths and loaders shared by the time-split suite."""

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RECORD = ROOT / "reference" / "time-split"
PAGE = ROOT / "docs" / "time-split.md"
SCRIPT = ROOT / "scripts" / "time-split" / "record.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
SUITE = "tests/time_split"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
ARMS = ("past", "random", "study-out")


def load_script():
    spec = importlib.util.spec_from_file_location("time_split_record", SCRIPT)
    module = sys.modules.get(spec.name)
    if module is None:
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    return module


def read_tsv(path):
    """A TSV of the record as a list of dicts, all cells as text."""
    text = pathlib.Path(path).read_text(encoding="utf-8")
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


def pairs(path):
    """A name and value TSV as a dict."""
    return {r["name"]: r["value"] for r in read_tsv(path)}
