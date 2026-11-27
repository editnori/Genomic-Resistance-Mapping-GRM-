"""The suite's wiring: the workflow's python job runs it, and the python row of docs/ci.md gives its count."""

import ast
import pathlib
import re
import unittest

from tests.matched_accuracy import common, stage2_common
from tests.matched_accuracy.common import CI_DOC, ROOT, SUITE, WORKFLOW

SCRIPTS = [common.SCRIPT, common.FAULTS, stage2_common.SCRIPT, stage2_common.FAULTS]
# numpy runs stage 2's resampling, and the tests that compute its interval apart from the script.
NUMPY = {"stage2.py", "test_stage2_measure.py", "test_stage2_plain.py", "test_stage2_made_up.py"}


def suite_count():
    loader = unittest.TestLoader()
    return loader.discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        job = text.split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertRegex(job, rf"\n          # SYNTHETIC [^\n]+\n(          # [^\n]+\n)*          - suite: "
                              rf"{re.escape(SUITE)}\n            top: \.\n(?!            ok:)")
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

    def test_ci_doc_says_the_suite_reads_no_release_report_and_the_committed_sums_list(self):
        text = CI_DOC.read_text(encoding="utf-8")
        line = next(line for line in text.split("\n") if line.startswith(f"The `{SUITE}` suite"))
        self.assertIn("reads no score and no release report", line)
        self.assertIn("reads the committed sums list of the release's reports", line)
        self.assertNotIn("no release file", line)

    def test_paths_inside_the_repository(self):
        for name in ("SCRIPT", "FAULTS", "PAGE", "INDEX", "WORKFLOW", "CI_DOC", "REFERENCE"):
            path = getattr(common, name)
            self.assertTrue(path.is_relative_to(ROOT), name)
            self.assertTrue(path.exists(), name)
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")) + SCRIPTS:
            source = path.read_text(encoding="utf-8")
            self.assertIsNone(re.search(r"[\"'](/(home|Users|mnt|tmp|root|opt|srv)/|[A-Za-z]:\\)", source), path.name)
            self.assertIsNone(re.search(r"[\"']\.\./", source), path.name)

    def test_the_suite_and_its_scripts_use_the_standard_library_alone_and_numpy_where_named(self):
        allowed = {"argparse", "ast", "collections", "contextlib", "copy", "csv", "datetime", "fractions", "functools",
                   "gzip", "hashlib", "importlib", "io", "json", "math", "operator", "os", "pathlib", "random", "re",
                   "shutil", "stat", "subprocess", "sys", "tempfile", "unittest", "tests"}
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")) + SCRIPTS:
            source = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                extra = {"numpy"} if path.name in NUMPY else set()
                self.assertLessEqual(set(re.findall(r"^(?:import|from) ([a-zA-Z_]+)", source, re.M)), allowed | extra)
                self.assertNotIn("url" + "lib", source)
                self.assertNotIn("sock" + "et", source)

    def test_stage_2_uses_numpy_for_the_resampling_and_its_version_alone(self):
        tree = ast.parse(stage2_common.SCRIPT.read_text(encoding="utf-8"))
        where = set()
        for f in tree.body:
            if isinstance(f, ast.FunctionDef):
                if any(isinstance(n, ast.Name) and n.id == "numpy" for n in ast.walk(f)):
                    where.add(f.name)
        self.assertEqual(where, {"check_numpy", "interval", "numpy_parts"})
        top = [n for f in tree.body if not isinstance(f, ast.FunctionDef) for n in ast.walk(f)
               if isinstance(n, ast.Name) and n.id == "numpy"]
        self.assertEqual(top, [])

    def test_no_test_of_the_suite_is_skipped(self):
        # The job requires unittest's last line to be OK, with no skip.
        for path in sorted(pathlib.Path(ROOT / SUITE).glob("*.py")):
            source = path.read_text(encoding="utf-8")
            self.assertNotIn("skip" + "Test", source, path.name)
            self.assertNotIn("unittest." + "skip", source, path.name)


if __name__ == "__main__":
    unittest.main()
