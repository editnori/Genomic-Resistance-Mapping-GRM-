"""docs/settings-arms-result.md: the page's form and the sentences it must hold whatever the count says.

Reads the page and scripts/c5-settings/lists.tsv; reads no score.
"""

import csv
import re
import unittest

from tests.c5_settings import common

PAGE = common.ROOT / "docs" / "settings-arms-result.md"
SECTIONS = ["Reading", "Data", "Count", "Per arm", "Beside it", "Seen after", "Limits", "Method", "Datasets"]
FIXED = [
    "Every split here is a retrospective split inside one public database. A split by lineage or by study is never "
    "external validation.",
    '"No difference shown" is the rule\'s reading and no statement that two rules are equal.',
    "Recall stands beside log loss wherever log loss is given.",
    "The cutoff arm moves no score, so its log loss equals the default's by construction.",
    "The release keeps its settings and its 0.5 cutoff whatever the count says, until the maintainer rules.",
    "The search covers 55 datasets in 8 species, those in both S and L and not in M.",
    "The 14 datasets of M wait for the label rule.",
    "*Klebsiella pneumoniae* with polymyxin B leaves the lineage count, as it left block C3's",
    "Its default draw on the lineage split at seed 43 came back unavailable, with no outer fold, and the result gives "
    "the reason as `invalid lineage outer evaluation: fold <k> validation lacks one or both classes`.",
    "The twin's draw at seed 46 is unavailable for the same reason.",
    "`scripts/c5-settings/score.py` changed after the run and before any score was opened.",
    "A column that reads 0 in every row is not printed and stays in its file only",
    "A time can equal a genome id as text, so the id check counts apart",
]


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

    def test_the_refusal_columns_are_not_said_to_read_0_by_construction(self):
        page = text()
        for m in re.finditer(r"by construction", page):
            line = page[page.rfind("\n", 0, m.start()) + 1:page.find("\n", m.end())]
            self.assertIn("cutoff arm", line)

    def test_no_em_dash_and_no_external_validation_claim(self):
        page = text()
        self.assertNotIn("—", page)
        for m in re.finditer(r"external validation", page):
            before = page[max(0, m.start() - 12):m.start()]
            self.assertTrue(before.endswith("never ") or before.endswith("not "), before)

    def test_every_first_run_species_and_antibiotic_is_named_in_full(self):
        page = text()
        with open(common.KIT / "lists.tsv", newline="", encoding="utf-8") as handle:
            run = [r for r in csv.DictReader(handle, delimiter="\t") if r["first_run"] == "1" and r["in_M"] == "0"]
        self.assertEqual(len(run), 55)
        for r in run:
            with self.subTest(dataset=r["dataset"]):
                self.assertIn(f"*{r['species']}*", page)
                row = [x for x in page.splitlines() if x.startswith(f"| *{r['species']}* |")]
                self.assertEqual(len(row), 1)
                self.assertIn(f"{r['antibiotic']} ({int(r['samples']):,})", row[0])

    def test_the_maintainer_rules_and_no_person_is_named(self):
        page = text()
        self.assertIn("until the maintainer rules", page)
        self.assertIn("Condition 6 is the maintainer's.", page)

    def test_the_plain_best_and_the_cutoff_log_loss_decide_nothing(self):
        page = text()
        self.assertIn("It stands in no Holm family and no condition, and decides nothing.", page)
        self.assertIn("Its log loss rows count every dataset as the same and decide nothing.", page)
        self.assertIn("The stricter reading stands beside ruling 15 and decides nothing.", page)


if __name__ == "__main__":
    unittest.main()
