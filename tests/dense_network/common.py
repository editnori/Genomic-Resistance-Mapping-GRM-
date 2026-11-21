"""Paths and loaders shared by the dense-network suite."""

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
RECORD = ROOT / "reference" / "dense-network"
PAGE = ROOT / "docs" / "dense-network.md"
SCRIPT = ROOT / "scripts" / "dense-network" / "record.py"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
SUITE = "tests/dense_network"


def load_script():
    spec = importlib.util.spec_from_file_location("dense_network_record", SCRIPT)
    module = sys.modules.get(spec.name)
    if module is None:
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
    return module


def read_tsv(path, header_only=False):
    """A TSV of the record as a list of dicts (or its header alone)."""
    text = pathlib.Path(path).read_text(encoding="utf-8")
    lines = text.rstrip("\n").split("\n")
    head = lines[0].split("\t")
    if header_only:
        return head
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]
