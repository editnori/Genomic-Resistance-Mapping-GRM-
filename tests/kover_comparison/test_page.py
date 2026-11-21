"""docs/kover-comparison.md: the page is what lines.page() builds from the record, every number in it is in the
record, and its form: the opening, one Reading section of at most 4 sentences, the sections in order."""

import json
import math
import re
import unittest

from tests.kover_comparison import lines as builder
from tests.kover_comparison.common import LEGACY, PAGE, RECORD

HEADERS = ["# Kover comparison", "## Reading", "## Plan", "## Datasets", "## Method", "## Results", "## Review",
           "## Limits", "## Files"]
TOKEN = re.compile(r"(?<![\w.])[-+]?\d[\d,]*(?:\.\d+)?")
# The Reading as committed, word for word: a rewording of the page and builder must change this too.
READING = [
    ("On these 12 datasets, under the plan's settings, Kover minus GRM's `scm` in balanced accuracy is "
     '-0.052 (-0.094 to -0.018, a 95% interval over the 12 datasets) as run and -0.035 (-0.060 to -0.013) '
     "when each tied first rule counts equally, an exploratory reading of Kover's 41 one-rule folds of 60."),
    ("Kover's balanced accuracy gap is in specificity, -0.104 (-0.215 to -0.022), not recall, +0.001 "
     '(-0.047 to +0.046); log loss, +0.374 (+0.071 to +0.748), counts the same errors because every call '
     'is 0 or 1.'),
    ("Kover's balanced accuracy is lower in 8 of the 12 datasets, equal in 2 and higher in 2."),
    ('Which differences between the runs make the gap is not separated, and the result, from the smallest '
     'eligible dataset of each of 12 species, is not shown to carry to larger datasets.'),
]


def record_numbers():
    """Every number written in a file of the record, and the counts of the historical file the page cites."""
    found = set()
    for path in sorted(RECORD.iterdir()):
        text = path.read_text(encoding="utf-8")
        # Numbers as words see them, and every run of digits, so that 16 in 2026-11-15T16:02:53Z counts.
        for t in TOKEN.findall(text) + re.findall(r"\d+(?:\.\d+)?", text):
            found.add(t.lstrip("+").replace(",", ""))
    legacy = json.loads(LEGACY.read_text(encoding="utf-8"))
    found.add(str(len(legacy)))
    found.update(str(int(float(r["ds_n_examples"]))) for r in legacy)
    return found


def computed(rec):
    """The numbers the page computes from cells of the record, with the cells they come from."""
    clip = float(rec.settings["log_loss_clip"])
    return {"%.4f" % math.log((1 - clip) / clip),
            str(int(rec.cohort["species"]["eligible"]) - int(rec.cohort["species"]["chosen"]))}


def found_in(token, numbers):
    """Whether a number of the page is a number of the record, as written or rounded to the page's decimals."""
    t = token.lstrip("+").replace(",", "")
    if t in numbers:
        return True
    if "." not in t:
        return False
    places = len(t.split(".")[1])
    for x in numbers:
        try:
            v = float(x)
        except ValueError:
            continue
        s = "%.*f" % (places, v)
        if s == "-0." + "0" * places:
            s = s[1:]
        if s == t:
            return True
    return False


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.lines = cls.text.split("\n")
        cls.rec = builder.Record()

    def test_page_is_built_from_the_record(self):
        built = builder.page(self.rec).split("\n")
        for i, (have, want) in enumerate(zip(self.lines, built), 1):
            self.assertEqual(have, want, "line %d" % i)
        self.assertEqual(len(self.lines), len(built))

    def test_every_number_is_in_the_record(self):
        numbers = record_numbers() | computed(self.rec)
        prose = re.sub(r"`[^`]*`", " ", re.sub(r"```.*?```", " ", self.text, flags=re.S))
        tokens = TOKEN.findall(prose)
        self.assertGreater(len(tokens), 250)
        for token in tokens:
            self.assertTrue(found_in(token, numbers), token)

    def test_a_number_not_in_the_record_is_caught(self):
        numbers = record_numbers()
        self.assertFalse(found_in("-0.999", numbers))
        self.assertFalse(found_in("0.5555", numbers))
        self.assertTrue(found_in("-0.104", numbers))
        self.assertTrue(found_in("1,118", numbers))

    def test_headers_in_order(self):
        self.assertEqual([x for x in self.lines if x.startswith("#")], HEADERS)

    def test_opening(self):
        first = self.text.split("\n\n")[1]
        self.assertNotIn("\n", first)
        self.assertTrue(first.endswith(builder.LIMIT.format(d=12)))
        self.assertIn("This is a comparison of two programs on 12 public datasets. It shows no equivalence and is not "
                      "external validation.", first)
        self.assertIn("It is not `reference/kover-summary.json`", first)
        self.assertIn("12 public BV-BRC datasets", first)

    def test_reading_is_at_most_4_sentences(self):
        sentences = builder.reading(self.text)
        self.assertLessEqual(len(sentences), 4)
        self.assertEqual(self.text.count("\n## Reading\n"), 1)
        for s in sentences:
            # One sentence per line: no stop followed by a capital inside a line.
            self.assertIsNone(re.search(r"[.!?] [A-Z]", s), s)
            self.assertTrue(s.endswith("."), s)

    def test_reading_gives_the_review_s_figures(self):
        r = builder.reading(self.text)
        self.assertIn("under the plan's settings", r[0])
        self.assertIn("-0.052 (-0.094 to -0.018, a 95% interval over the 12 datasets) as run and", r[0])
        self.assertIn("-0.035 (-0.060 to -0.013) when each tied first rule counts equally, an exploratory reading", r[0])
        self.assertIn("specificity, -0.104 (-0.215 to -0.022), not recall, +0.001 (-0.047 to +0.046)", r[1])
        self.assertIn("lower in 8 of the 12 datasets", r[2])
        # The plan kept the per-dataset t interval out of the reading.
        self.assertNotIn("t interval", " ".join(r))
        # Recall stands beside log loss.
        self.assertIn("log loss", r[1])

    def test_reading_word_for_word(self):
        self.assertEqual(builder.reading(self.text), READING)

    def test_reading_numbers_carry_their_denominators(self):
        r = builder.reading(self.text)
        self.assertIn("41 one-rule folds of 60", r[0])
        self.assertIn("over the 12 datasets", r[0])
        self.assertIn("of the 12 datasets", r[2])

    def test_datasets_named(self):
        table = [x for x in self.lines if x.startswith("| ") and "__" in x]
        names = [x.split(" | ")[0][2:] for x in table]
        self.assertEqual(names, [d["dataset"] for d in self.rec.datasets])
        self.assertEqual(len(names), 12)

    def test_version_and_settings(self):
        self.assertIn("Kover 2.0.7 (core-2.0.7, cli-2.0.1)", self.text)
        for key in ("kover_dataset", "kover_split", "kover_learn"):
            self.assertIn("\n" + self.rec.settings[key] + "\n", self.text)

    def test_tables_are_whole(self):
        for line in self.lines:
            if line.startswith("|"):
                self.assertTrue(line.endswith("|"), line)
        for block in re.findall(r"(?:^\|.*\n)+", self.text, flags=re.M):
            rows = block.rstrip("\n").split("\n")
            widths = {r.count("|") for r in rows}
            self.assertEqual(len(widths), 1, rows[0])
            self.assertRegex(rows[1], r"^\|(---\|)+$")

    def test_what_it_cannot_show(self):
        limits = self.text.split("\n## Limits\n", 1)[1].split("\n## ", 1)[0]
        self.assertIn("The interval covers the choice of datasets only.", limits)
        self.assertIn("The spread over seeds is not measured.", limits)
        self.assertIn("not comparable", limits)
        self.assertIn("not shown to carry to larger datasets", self.text)
        # The first review's seed finding as it found it: the intervals move, and no reading changes.
        self.assertIn("Under the plan's seed the intervals for recall, specificity and log loss move in the third "
                      "decimal, and no reading changes.", self.text)
        self.assertNotIn("no interval changes", self.text)


if __name__ == "__main__":
    unittest.main()
