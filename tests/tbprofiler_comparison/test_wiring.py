"""The suite's wiring: the workflow's python job runs it, and the python row of docs/ci.md gives its count."""

import pathlib
import re
import unittest

from tests.tbprofiler_comparison import common
from tests.tbprofiler_comparison.common import CI_DOC, ROOT, SUITE, WORKFLOW


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
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")) + [common.SCRIPT]:
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(re.search(r"[\"'](/(home|Users|mnt|tmp|root|opt|srv)/|[A-Za-z]:\\)", source), path.name)
            self.assertIsNone(re.search(r"[\"']\.\./", source), path.name)

    def test_the_script_uses_the_standard_library_alone(self):
        source = common.SCRIPT.read_text(encoding="utf-8")
        imported = {a or b for a, b in re.findall(r"^(?:import ([a-zA-Z_]+)\s*$|from ([a-zA-Z_]+)[\w.]* import )", source,
                                                  re.M)}
        self.assertGreaterEqual(len(imported), 8)
        self.assertLessEqual(imported, {"argparse", "collections", "gzip", "hashlib", "json", "math", "random", "re",
                                        "sys", "pathlib"})
        self.assertNotIn("urllib", source)
        self.assertNotIn("socket", source)
        self.assertNotIn("subprocess", source)


if __name__ == "__main__":
    unittest.main()
