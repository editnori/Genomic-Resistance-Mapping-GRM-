"""The suite's wiring: the workflow's python job runs it, the python row of docs/ci.md gives its count, and the scripts use the standard library alone."""

import ast
import re
import sys
import unittest

from tests.mechanism.common import KIT, ROOT

SUITE = "tests/mechanism"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
OWN = {"amr_rule", "common", "strata", "history", "combine", "scan", "drug_rules", "planted_faults", "tests"}


def suite_count():
    return unittest.TestLoader().discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        job = WORKFLOW.read_text(encoding="utf-8").split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertRegex(job, rf"\n          # SYNTHETIC [^\n]+\n          - suite: {re.escape(SUITE)}\n            top: \.\n")
        self.assertEqual(job.count(f"suite: {SUITE}\n"), 1)

    def test_python_row_gives_the_count(self):
        rows = [line for line in CI_DOC.read_text(encoding="utf-8").split("\n") if line.startswith("| `python`")]
        self.assertEqual(len(rows), 1)
        m = re.search(rf"`{re.escape(SUITE)}` (\d+)(,| \|)", rows[0])
        self.assertIsNotNone(m, "the python row does not give the suite's count")
        self.assertEqual(int(m.group(1)), suite_count())

    def test_ci_doc_names_the_suite_once_in_its_text(self):
        self.assertEqual(CI_DOC.read_text(encoding="utf-8").count(f"The `{SUITE}` suite"), 1)

    def test_the_scripts_use_the_standard_library_alone(self):
        paths = sorted(KIT.glob("*.py")) + sorted((ROOT / SUITE).glob("*.py"))
        for path in paths:
            tree = ast.parse(path.read_text(encoding="utf-8"))
            for node in ast.walk(tree):
                names = []
                if isinstance(node, ast.Import):
                    names = [a.name.split(".")[0] for a in node.names]
                elif isinstance(node, ast.ImportFrom) and node.level == 0:
                    names = [node.module.split(".")[0]]
                for name in names:
                    self.assertTrue(name in sys.stdlib_module_names or name in OWN, f"{path.name}: {name}")


if __name__ == "__main__":
    unittest.main()
