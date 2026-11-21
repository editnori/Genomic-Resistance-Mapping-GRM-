"""The suite's wiring: the workflow's python job runs it, and the python row of docs/ci.md gives its count."""

import ast
import pathlib
import re
import unittest

from tests.amrfinderplus_comparison import common
from tests.amrfinderplus_comparison.common import CI_DOC, ROOT, SUITE, WORKFLOW


def suite_count():
    loader = unittest.TestLoader()
    return loader.discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        job = text.split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertRegex(job, rf"\n          # SYNTHETIC [^\n]+\n          - suite: {re.escape(SUITE)}\n"
                              rf"            top: \.\n")
        self.assertEqual(job.count(f"suite: {SUITE}\n"), 1)
        self.assertLess(job.index("suite: tests/dense_network\n"), job.index(f"suite: {SUITE}\n"))

    def test_python_row_gives_the_count(self):
        rows = [line for line in CI_DOC.read_text(encoding="utf-8").split("\n") if line.startswith("| `python`")]
        self.assertEqual(len(rows), 1)
        m = re.search(rf"`{re.escape(SUITE)}` (\d+)(,| \|)", rows[0])
        self.assertIsNotNone(m, "the python row does not give the suite's count")
        self.assertEqual(int(m.group(1)), suite_count())

    def test_ci_doc_names_the_suite_once_in_its_text(self):
        text = CI_DOC.read_text(encoding="utf-8")
        self.assertEqual(text.count(f"The `{SUITE}` suite"), 1)

    def test_paths_inside_the_repository(self):
        for name in ("RECORD", "PAGE", "SCRIPT", "WORKFLOW", "CI_DOC"):
            path = getattr(common, name)
            self.assertTrue(path.is_relative_to(ROOT), name)
            self.assertTrue(path.exists(), name)
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")):
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(re.search(r"[\"'](/(home|Users|mnt|tmp|root|opt|srv)/|[A-Za-z]:\\)", source), path.name)
            self.assertIsNone(re.search(r"[\"']\.\./", source), path.name)

    def test_script_uses_the_standard_library_only(self):
        tree = ast.parse(common.SCRIPT.read_text(encoding="utf-8"))
        imports = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
        imports |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
        self.assertEqual(imports - {"argparse", "fractions", "gzip", "hashlib", "json", "math", "pathlib", "random",
                                    "re", "statistics", "sys"}, set())


if __name__ == "__main__":
    unittest.main()
