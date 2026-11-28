"""Wiring and hygiene tests for tests/kmer_bench and scripts/kmer-bench."""

import ast
import json
import re
import unittest

from tests.kmer_bench import common

runners = common.load_module("runners")


def suite_count() -> int:
    loader = unittest.TestLoader()
    return loader.discover(start_dir=str(common.ROOT / common.SUITE), top_level_dir=str(common.ROOT)).countTestCases()


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        text = common.WORKFLOW.read_text(encoding="utf-8")
        job = text.split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertRegex(
            job,
            rf"\n          # SYNTHETIC [^\n]+\n(          # [^\n]+\n)*          - suite: "
            rf"{re.escape(common.SUITE)}\n            top: \.\n(?!            ok:)",
        )
        self.assertEqual(job.count(f"suite: {common.SUITE}\n"), 1)

    def test_python_row_gives_the_count(self):
        rows = [line for line in common.CI_DOC.read_text(encoding="utf-8").split("\n") if line.startswith("| `python`")]
        self.assertEqual(len(rows), 1)
        m = re.search(rf"`{re.escape(common.SUITE)}` (\d+)(,| \|)", rows[0])
        self.assertIsNotNone(m, "the python row does not give the suite's count")
        self.assertEqual(int(m.group(1)), suite_count())

    def test_ci_doc_names_the_suite_once_in_its_text(self):
        text = common.CI_DOC.read_text(encoding="utf-8")
        self.assertEqual(text.count(f"The `{common.SUITE}` suite"), 1)

    def test_golden_grm_fixture_exists_and_passes_hygiene(self):
        self.assertTrue(common.GOLDEN_GRM.is_file())
        data = json.loads(common.GOLDEN_GRM.read_text(encoding="utf-8"))
        self.assertEqual(data.get("format"), "grm-kmer-bench-golden-v1")
        self.assertEqual(runners.scan_json_for_leaks(data), [])

    def test_scripts_and_tests_use_standard_library_only(self):
        allowed = {
            "argparse", "ast", "collections", "contextlib", "dataclasses", "gzip", "hashlib",
            "importlib", "io", "json", "os", "pathlib", "re", "shutil", "statistics",
            "struct", "subprocess", "sys", "tempfile", "threading", "time",
            "typing", "unittest", "tests",
        }
        py_files = sorted(common.SCRIPT_DIR.glob("*.py")) + sorted((common.ROOT / common.SUITE).glob("*.py"))
        self.assertGreaterEqual(len(py_files), 8)
        for path in py_files:
            with self.subTest(file=path.relative_to(common.ROOT).as_posix()):
                tree = ast.parse(path.read_text(encoding="utf-8"))
                imports = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
                imports |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
                self.assertEqual(imports - allowed, set())

    def test_no_em_dashes_private_paths_or_skips(self):
        checked = [common.GOLDEN_GRM] + sorted(common.SCRIPT_DIR.glob("*.py")) + sorted((common.ROOT / common.SUITE).glob("*.py"))
        for path in checked:
            text = path.read_text(encoding="utf-8")
            rel = path.relative_to(common.ROOT).as_posix()
            self.assertNotIn("\u2014", text, f"em dash in {rel}")
            if path.suffix == ".py" and path.parent == common.ROOT / common.SUITE:
                self.assertNotIn("skip" + "Test", text, rel)
                self.assertNotIn("unittest." + "skip", text, rel)


if __name__ == "__main__":
    unittest.main()
