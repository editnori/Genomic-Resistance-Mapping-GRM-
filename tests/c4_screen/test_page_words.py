"""docs/screen-arms-result.md: its sections, the sentences it must carry whatever the count says, and its tables of
gates and faults as the scripts hold them. Reads the page and the scripts; reads no score.
"""

import re
import unittest

from tests.c4_screen import common

PAGE = common.ROOT / "docs" / "screen-arms-result.md"
page = common.load("page_tables")
faults = common.load("page_faults")
SECTIONS = ["Reading", "Data", "Count", "Per arm", "Beside it", "Seen after", "Limits", "Method", "Datasets"]
FIXED = [
    "Every split here is a retrospective split inside one public database. A split by lineage or by study is never "
    "external validation.",
    '"No difference shown" is the rule\'s reading',
    "It is no statement that two screens are equal.",
    "Recall stands beside log loss wherever log loss is given.",
    "The release keeps its screen, its budget and its floor whatever the count says, until the maintainer rules.",
    "The Cochran-Mantel-Haenszel screen adjusts the ranking of the candidates for lineage strata and names no marker.",
    "The datasets of M wait for the label rule.",
]


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")

    def test_the_sections_stand_in_order(self):
        self.assertEqual(re.findall(r"^## (.+)$", self.text, flags=re.M), SECTIONS)

    def test_the_fixed_sentences_stand(self):
        for sentence in FIXED:
            self.assertEqual(self.text.count(sentence), 1, sentence)

    def test_no_long_dash(self):
        self.assertNotIn("—", self.text)

    def test_the_gate_table_is_the_scripts(self):
        rows = re.findall(r"^\| (\d+)\. (\w+) \| (.+) \|$", self.text, flags=re.M)
        self.assertEqual([(int(n), name, words) for n, name, words in rows], list(page.GATES))

    def test_the_fault_table_is_the_runners(self):
        rows = re.findall(r"^\| ([PR]\d+) \| (.+?) \| (.+?) \| `(.+?)` \|$", self.text, flags=re.M)
        want = []
        for n, what, _, _, _, gate, test in faults.FAULTS:
            label = {"check": "check script", "table": "a table"}.get(gate, gate)
            want.append((n, what, label, ".".join(test.rsplit(".", 2)[-2:])))
        self.assertEqual(rows, want)


if __name__ == "__main__":
    unittest.main()
