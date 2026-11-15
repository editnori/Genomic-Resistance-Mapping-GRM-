"""Paths and loaders shared by the clindamycin label suite."""

import importlib.util
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "clindamycin-labels" / "rev14-20261115"
BUNDLED = ROOT / "web" / "src" / "data" / "clindamycin-labels.json"
PAGE = ROOT / "docs" / "clindamycin-labels.md"
RELEASE_SUMS = ROOT / "reference" / "lineage-control" / "rev14-274-20261112.reports.sha256"


def load_site_data():
    spec = importlib.util.spec_from_file_location("clindamycin_labels_site_data", ROOT / "scripts" / "clindamycin-labels" / "site_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def bundled():
    return json.loads(BUNDLED.read_text(encoding="utf-8"))
