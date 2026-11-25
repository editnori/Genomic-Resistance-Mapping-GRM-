"""docs/combiner-result.md against its script, scripts/mechanism/combiner.py, before any reading.

This test file reads the page and the script's own constants, gates and fault list. It reads no record, fits
nothing and reads no genome.
"""

import re
import unittest

from tests.mechanism import common

cb = common.load("combiner")
cm = common.load("common")
faults = common.load("combiner_faults")

PAGE = common.ROOT / "docs" / "combiner-result.md"
SCRIPT = common.KIT / "combiner.py"
RECORD = common.ROOT / "reference" / "combiner"
SECTIONS = ["Plan", "Reading", "Limits", "Method"]
FIXED = ("A lineage split is a retrospective split inside one database, never external validation.",
         "\"No difference shown\" is never written as \"equal\".",
         "The error levels serve as an adapted research screen and are never an applicable FDA validation criterion.",
         "Coverage is the share of genomes answered, not accuracy.")
COMMANDS = ("python3 -I scripts/mechanism/combiner.py read REPORTS --first reference/learner-overlap "
            "--second reference/learner-overlap/second-pass --out READING",
            "python3 -I scripts/mechanism/combiner.py write READING --reports REPORTS --out OUT",
            "python3 -I scripts/mechanism/combiner.py print OUT")


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def ordinal(n):
    return f"{n:,}" + ({1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th") if n % 100 not in (11, 12, 13) else "th")


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.sections = sections(self.text)

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_the_fixed_words_stand_in_the_limits(self):
        for words in FIXED:
            self.assertIn(words, self.sections["Limits"])

    def test_no_difference_shown_is_never_written_as_equal(self):
        rest = self.text.replace("\"equal\"", "")
        self.assertEqual(re.findall(r"\bequal\b(?! weights)", rest, flags=re.I), [])
        self.assertIn("no difference shown", self.sections["Plan"])

    def test_the_plan_gives_the_scripts_constants(self):
        plan = self.sections["Plan"]
        lo_value, hi_value = cm.LO_INDEX + 1, cm.HI_INDEX + 1
        for words in (f"a score of {cb.CUTOFF} or more", f"sum to {cb.VOTE_SUM} or more by `math.fsum`",
                      f"clipped to [{cb.CLIP:.4f}, {1 - cb.CLIP:.4f}]", f"is {cb.MARGIN} or more",
                      f"{cm.RESAMPLES:,} resamples of species", f"the {ordinal(lo_value)} and {ordinal(hi_value)}",
                      f"under {cm.MIN_SPECIES} species", "five learners on the lineage split and 10 or more resistant"):
            self.assertIn(words, plan)
        self.assertEqual(cb.FORMS, ("mean", "weighted", "logistic"))
        self.assertEqual(lo_value, 101)

    def test_the_reading_waits_for_its_record(self):
        if not RECORD.exists():
            self.assertEqual(self.sections["Reading"].strip(), "Not read yet.")

    def test_the_commands_are_the_scripts(self):
        method = self.sections["Method"]
        for command in COMMANDS:
            self.assertIn(command + "\n", method)

    def test_the_gate_table_names_each_gate_of_the_script(self):
        gates = sorted(set(re.findall(r"\"gate (\d) ([a-z -]+):", SCRIPT.read_text(encoding="utf-8"))))
        self.assertEqual([n for n, _ in gates], [str(i) for i in range(1, 9)])      # one name for each gate
        rows = re.findall(r"^\| (\d)\. ([a-z -]+) \|", self.sections["Method"], flags=re.M)
        self.assertEqual(rows, gates)

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        rows = re.findall(r"^\| (C\d+)\. (.+) \| `([A-Za-z]+\.test_\w+)` \|$", method, flags=re.M)
        self.assertEqual(rows, [(fault, words[0].upper() + words[1:], test.rsplit(".", 2)[-2] + "." + test.rsplit(".", 1)[-1])
                                for fault, words, test, _ in faults.FAULTS])
        self.assertIn(f"All {len(faults.FAULTS)} were caught.", method)

    def test_the_test_counts_are_the_suites(self):
        loader = unittest.TestLoader()
        script = loader.loadTestsFromName("tests.mechanism.test_combiner").countTestCases()
        page = loader.loadTestsFromName("tests.mechanism.test_combiner_page").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_combiner.py` holds {script} SYNTHETIC tests of the script.", method)
        self.assertIn(f"`tests/mechanism/test_combiner_page.py` holds {page} tests of this page against the script.",
                      method)


if __name__ == "__main__":
    unittest.main()
