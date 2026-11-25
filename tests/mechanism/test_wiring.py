"""The suite's wiring: the workflow's python job runs it and its comment and docs/ci.md give the faults the suite plants, each count of a runner's faults in docs/ci.md and docs/mechanism-split.md is its own list's, the python row of docs/ci.md and the Checks of docs/mechanism-split.md give the suite's own test counts, and the scripts use the standard library alone."""

import ast
import re
import sys
import unittest

from tests.mechanism.common import KIT, ROOT, load

SUITE = "tests/mechanism"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
SPLIT_DOC = ROOT / "docs" / "mechanism-split.md"
OWN = {"amr_rule", "common", "strata", "history", "combine", "scan", "drug_rules", "planted_faults", "tests",
       "learner_overlap"}


UNITS = ("zero one two three four five six seven eight nine ten eleven twelve thirteen fourteen fifteen sixteen "
         "seventeen eighteen nineteen").split()
TENS = "zero ten twenty thirty forty fifty sixty seventy eighty ninety".split()
COUNT = r"([A-Za-z]+(?:-[a-z]+)?|\d+)"


def suite_count():
    return unittest.TestLoader().discover(start_dir=str(ROOT / SUITE), top_level_dir=str(ROOT)).countTestCases()


def in_words(n):
    """A count under 100 as the docs write it in words: 17 is seventeen and 49 is forty-nine."""
    if n < 20:
        return UNITS[n]
    tens, unit = divmod(n, 10)
    return TENS[tens] + (f"-{UNITS[unit]}" if unit else "")


def as_number(written):
    """A count written in digits or in words, a sentence's first word capitalised, or None."""
    words = {in_words(n): n for n in range(100)}
    return int(written) if written.isdigit() else words.get(written.lower())


def runners():
    return {path.stem: len(load(path.stem).FAULTS) for path in sorted(KIT.glob("*_faults.py"))}


class Wiring(unittest.TestCase):
    def test_in_the_python_job(self):
        job = WORKFLOW.read_text(encoding="utf-8").split("\n  python:\n", 1)[1].split("\n    steps:", 1)[0]
        self.assertRegex(job, rf"\n          # SYNTHETIC [^\n]+\n          - suite: {re.escape(SUITE)}\n            top: \.\n")
        self.assertEqual(job.count(f"suite: {SUITE}\n"), 1)
        comment = re.search(rf"\n          (# SYNTHETIC [^\n]+)\n          - suite: {re.escape(SUITE)}\n", job).group(1)
        self.assertTrue(runners())
        planted, how_many = sum(runners().values()), in_words(len(runners()))
        self.assertIn(f" plants the {planted} faults of its {how_many} runners (scripts/mechanism/*_faults.py) ", comment)
        self.assertEqual(CI_DOC.read_text(encoding="utf-8").count(
            f"how many faults its {how_many} runners plant, {planted} in all."), 1)

    def test_each_count_of_a_runners_faults_is_its_own_lists(self):
        # in docs/ci.md, the clause of the suite's text that ends at each runner's path; in docs/mechanism-split.md,
        # the sentence that gives planted_faults.py's faults and the two parts they come from
        text = CI_DOC.read_text(encoding="utf-8")
        start = text.index(f"The `{SUITE}` suite")
        suite = text[start:text.index("The `tests/", start + 1)]
        for stem, n in runners().items():
            written = re.findall(rf"\b{COUNT} faults [^()]*\(`scripts/mechanism/{stem}\.py`\)", suite)
            self.assertEqual([as_number(w) for w in written], [n], f"docs/ci.md, {stem}.py")
        parts = re.findall(rf"`planted_faults\.py` plants {COUNT} faults in turn, each caught by its test: "
                           rf"{COUNT} from [^;]*?, and {COUNT} from the first run", SPLIT_DOC.read_text(encoding="utf-8"))
        self.assertEqual(len(parts), 1, "docs/mechanism-split.md, planted_faults.py")
        total, first, then = (as_number(w) for w in parts[0])
        self.assertEqual(total, runners()["planted_faults"])
        self.assertEqual((first or 0) + (then or 0), total)

    def test_the_checks_give_the_suites_own_test_counts(self):
        lines = [line for line in SPLIT_DOC.read_text(encoding="utf-8").split("\n")
                 if line.startswith(f"`{SUITE}` has ")]
        self.assertEqual(len(lines), 1)
        loader = unittest.TestLoader()
        n = {name: loader.loadTestsFromName(f"tests.mechanism.test_{name}").countTestCases() for name in (
            "tb_curve", "page_tables", "page_faults", "curve_tables", "curve_faults", "learner_overlap", "overlap_second",
            "overlap_faults", "catalogue_size", "catalogue_faults", "page", "lineage_page", "learner_overlap_page",
            "catalogue_page")}
        total = suite_count()
        synthetic = total - n["page"] - n["lineage_page"] - n["learner_overlap_page"] - n["catalogue_page"]
        parts = (f"`{SUITE}` has {total} tests: {synthetic} SYNTHETIC tests, {n['tb_curve']} of them on `tb_curve.py`, ",
                 f", {n['page_tables']} on the result page's table script and {n['page_faults']} on its planted faults, ",
                 f", {n['curve_tables']} on the lineage curve's table script and {n['curve_faults']} on its planted "
                 "faults, ",
                 f", {n['learner_overlap']} on the learner overlap's table script (`docs/learner-overlap-result.md`), "
                 f"{n['overlap_second']} on its second pass and {n['overlap_faults']} on their planted faults, ",
                 f", and {n['catalogue_size']} on the catalogue size reading's script "
                 f"(`docs/catalogue-size-result.md`) and {n['catalogue_faults']} on its planted faults; ",
                 f"; {n['page']} of the result page against its record; ",
                 f"; {n['lineage_page']} of the lineage curve's page against its record; ",
                 f"; {n['learner_overlap_page']} of the learner overlap's page against its records; ",
                 f"; and {n['catalogue_page']} of the catalogue size page against its record.")
        for part in parts:
            self.assertIn(part, lines[0])

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
