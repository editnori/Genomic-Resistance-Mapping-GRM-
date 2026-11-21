"""Paths and loaders shared by the TBProfiler comparison suite."""

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RECORD = ROOT / "reference" / "tbprofiler-comparison"
PAGE = ROOT / "docs" / "tbprofiler-comparison.md"
SCRIPT = ROOT / "scripts" / "tbprofiler-comparison" / "record.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
SUITE = "tests/tbprofiler_comparison"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
SPLITS = ("lineage", "random")
WAY = "who_assoc"


def load_script():
    spec = importlib.util.spec_from_file_location("tbprofiler_comparison_record", SCRIPT)
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
