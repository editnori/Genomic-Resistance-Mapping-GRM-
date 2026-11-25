"""docs/outside-learners-result.md: the page's form and the sentences it must hold whatever the count says.

Reads the page and scripts/c7-outside/lists.tsv; reads no score.
"""

import csv
import re
import unittest

from tests.c7_outside import common

PAGE = common.ROOT / "docs" / "outside-learners-result.md"
SECTIONS = ["Reading", "Data", "Count", "Per arm", "Beside it", "Seen after", "Limits", "Method", "Datasets"]
FIXED = [
    "Every split here is a retrospective split inside one public database. A split by lineage or by study is never "
    "external validation.",
    '"No difference shown" is the rule\'s reading and no statement that two learners are equal.',
    "Recall stands beside log loss wherever log loss is given.",
    "The release keeps its five learners and their settings whatever the count says, until the maintainer rules.",
    "Each booster runs GRM's settings where GRM has them, on the rows GRM's screen kept, so the pair reads the library "
    "and leaves the screen as it is.",
    "The elastic net stops at tol 1e-4 or 2,000 epochs, and every elastic net result carries that limit as "
    "`docs/outside-learners.md` words it: saga's tol of 1e-4 stops on the change in the weights and bounds no "
    "distance to the optimum",
    "No feature, unitig or k-mer is named or ranked.",
    "The 35 datasets of M wait for the label rule.",
    "Condition 6 is the maintainer's.",
]
# The one sentence a withdrawn elastic net would have changed. No arm was withdrawn.
ARMS_SENTENCE = "**Arms.** All three arms ran: `xgboost`, `lightgbm` and `elasticnet`."


def text():
    return PAGE.read_text(encoding="utf-8")


class Page(unittest.TestCase):
    def test_the_sections_stand_in_order(self):
        heads = re.findall(r"^## (.+)$", text(), flags=re.M)
        self.assertEqual(heads, SECTIONS)

    def test_each_fixed_sentence_is_there(self):
        page = text()
        for sentence in FIXED:
            with self.subTest(sentence=sentence[:40]):
                self.assertIn(sentence, page)

    def test_no_em_dash_and_no_external_validation_claim(self):
        page = text()
        self.assertNotIn("—", page)
        for m in re.finditer(r"external validation", page):
            before = page[max(0, m.start() - 12):m.start()]
            self.assertTrue(before.endswith("never ") or before.endswith("not "), before)

    def test_every_first_run_species_and_antibiotic_is_named_in_full(self):
        page = text().split("\n## Datasets\n", 1)[1]
        with open(common.KIT / "lists.tsv", newline="", encoding="utf-8") as handle:
            run = [r for r in csv.DictReader(handle, delimiter="\t") if r["first_run"] == "1"]
        self.assertEqual(len(run), 233)
        self.assertEqual(len({r["species"] for r in run}), 23)
        for r in run:
            with self.subTest(dataset=r["dataset"]):
                row = [x for x in page.splitlines() if x.startswith(f"| *{r['species']}* |")]
                self.assertEqual(len(row), 1)
                self.assertIn(f"{r['antibiotic']} ({int(r['samples']):,})", row[0])

    def test_the_arms_sentence_stands_once_with_no_arm_withdrawn(self):
        page = text()
        self.assertEqual(page.count(ARMS_SENTENCE), 1)
        self.assertEqual(page.count("All three arms ran"), 1)
        self.assertIn(ARMS_SENTENCE + " No arm was withdrawn.\n", page)

    def test_no_part_2_mark_is_left(self):
        self.assertNotIn("<!--", text())

    def test_the_gates_table_is_the_scripts(self):
        page = common.load("page_tables")
        rows = re.findall(r"^\| (\d+)\. ([a-zA-Z]+) \|", text(), flags=re.M)
        self.assertEqual([(int(n), name) for n, name in rows], [(n, name) for n, name, _ in page.GATES])


if __name__ == "__main__":
    unittest.main()
