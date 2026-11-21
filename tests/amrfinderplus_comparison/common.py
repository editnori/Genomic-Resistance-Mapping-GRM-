"""Paths and loaders shared by the AMRFinderPlus comparison suite."""

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RECORD = ROOT / "reference" / "amrfinderplus-comparison"
PAGE = ROOT / "docs" / "amrfinderplus-comparison.md"
SCRIPT = ROOT / "scripts" / "amrfinderplus-comparison" / "record.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
SUITE = "tests/amrfinderplus_comparison"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")


def load_script():
    spec = importlib.util.spec_from_file_location("amrfinderplus_comparison_record", SCRIPT)
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
    """A key and value TSV as a dict."""
    return {r["key"]: r["value"] for r in read_tsv(path)}
