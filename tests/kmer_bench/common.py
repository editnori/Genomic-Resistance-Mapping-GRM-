"""Shared loader and repository paths for tests/kmer_bench."""

import importlib.util
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SUITE = "tests/kmer_bench"
SCRIPT_DIR = ROOT / "scripts" / "kmer-bench"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
INDEX_DOC = ROOT / "docs" / "INDEX.md"
GOLDEN_GRM = ROOT / SUITE / "golden_grm.json"


def load_module(name: str):
    """Loads one module from scripts/kmer-bench/ cleanly under -I and -s -P."""
    mod_name = f"_grm_kmer_bench_{name}"
    if mod_name in sys.modules:
        return sys.modules[mod_name]
    path = SCRIPT_DIR / f"{name}.py"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod
