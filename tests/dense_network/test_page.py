"""docs/dense-network.md: the page is what lines.page() builds from the record, every number in it is in the
record, and its form: the opening, one Reading section of at most 4 sentences, the sections in order."""

import re
import unittest

from tests.dense_network import lines as builder
from tests.dense_network.common import PAGE, RECORD

HEADERS = ["# Dense network", "## Reading", "## Plan", "## Network", "## Datasets", "## Method", "## Results",
           "## Review", "## Limits", "## Files"]
TOKEN = re.compile(r"(?<![\w.])[-+]?\d[\d,]*(?:\.\d+)?")


def record_numbers():
    """Every number written in a file of the record."""
    found = set()
    for path in sorted(RECORD.iterdir()):
        text = path.read_text(encoding="utf-8")
        # Numbers as words see them, and every run of digits, so that 19 in 2026-11-15T12:19:40Z counts.
        for t in TOKEN.findall(text) + re.findall(r"\d+(?:\.\d+)?", text):
            found.add(t.lstrip("+").replace(",", ""))
    return found


def computed(rec):
    """The numbers the page computes from cells of the record."""
    low, high = re.search(r"([\d.]+)% and ([\d.]+)% percentiles", rec.settings["interval"]).groups()
    return {"%g" % (float(high) - float(low))} | {str(x) for x in rec.per_species()}


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
        self.assertFalse(found_in("0.4321", numbers))
        self.assertTrue(found_in("-0.059", numbers))
        self.assertTrue(found_in("1,957", numbers))

    def test_headers_in_order(self):
        self.assertEqual([x for x in self.lines if x.startswith("#")], HEADERS)

    def test_opening(self):
        first = self.text.split("\n\n")[1]
        self.assertNotIn("\n", first)
        self.assertTrue(first.endswith(" " + builder.LIMIT))
        self.assertEqual(builder.LIMIT, "This is one experiment on public data. It establishes no gain and rules none out.")
        self.assertTrue(first.startswith("This page compares one dense feed-forward network with GRM's five learners, "
                                         "`scm`, `cart`, `rf`, `gbt` and `lasso`, on 24 public BV-BRC datasets"))
        self.assertIn("on the release's own lineage and random folds", first)

    def test_reading_is_at_most_4_sentences(self):
        sentences = builder.reading(self.text)
        self.assertLessEqual(len(sentences), 4)
        self.assertEqual(self.text.count("\n## Reading\n"), 1)
        for s in sentences:
            # One sentence per line: no stop followed by a capital inside a line.
            self.assertIsNone(re.search(r"[.!?] [A-Z]", s), s)
            self.assertTrue(s.endswith("."), s)

    def test_reading_gives_the_figures(self):
        r = builder.reading(self.text)
        self.assertIn("over the 24 datasets is -0.059 (-0.132 to -0.044, a 95% interval over the 12 species)", r[0])
        self.assertIn("the network is higher on 1 of the 24, and the plan's rule gives \"no gain established\".", r[0])
        self.assertIn("in recall is -0.102 (-0.356 to 0.000) and in log loss +0.021 (-0.262 to +0.209)", r[1])
        self.assertIn("-0.012 (-0.052 to -0.001) for `cart`", r[2])
        self.assertIn("so only `cart`'s interval lies below 0 and none lies above it.", r[2])
        self.assertIn("no transformer or sequence model was fitted.", r[3])

    def test_reading_numbers_carry_their_denominators_and_each_learner_its_effect(self):
        r = builder.reading(self.text)
        self.assertIn("over the 24 datasets", r[0])
        self.assertIn("over the 12 species", r[0])
        self.assertIn("over the same 24 datasets", r[1])
        self.assertIn("over the 24 datasets", r[2])
        for who in builder.FIVE:
            self.assertEqual(len(re.findall(r"\([-+0-9.]+ to [-+0-9.]+\) for `%s`" % who, r[2])), 1, who)
        # Recall stands beside log loss.
        self.assertRegex(r[1], r"recall is .+ and in log loss ")

    def test_datasets_named(self):
        table = [x for x in self.lines if x.startswith("| ") and "__" in x]
        names = [x.split(" | ")[0][2:] for x in table]
        self.assertEqual(names, [d["dataset"] for d in self.rec.datasets])
        self.assertEqual(len(names), 24)

    def test_the_network_as_the_plan_fixed_it(self):
        net = self.text.split("\n## Network\n", 1)[1].split("\n## ", 1)[0]
        for key in ("network", "inputs", "start", "loss", "optimiser", "stopping", "repeats", "splits"):
            value = self.rec.settings[key]
            self.assertIn(value[0].upper() + value[1:] + ".\n", net, key)
        self.assertIn("4,000 inputs, 256 ReLU units, 64 ReLU units, 1 sigmoid output", net)

    def test_plan_and_its_order(self):
        plan = self.text.split("\n## Plan\n", 1)[1].split("\n## ", 1)[0]
        self.assertIn("Its file time is 12:19:40Z and its sha256 is `%s`." % self.rec.order["plan_file_time"]["sha256"], plan)
        self.assertIn("the first network result was written at 12:21:31Z", plan)
        self.assertIn("The order rests on file times alone", plan)

    def test_tables_are_whole(self):
        for line in self.lines:
            if line.startswith("|"):
                self.assertTrue(line.endswith("|"), line)
        blocks = re.findall(r"(?:^\|.*\n)+", self.text, flags=re.M)
        self.assertEqual(len(blocks), 4)
        for block in blocks:
            rows = block.rstrip("\n").split("\n")
            widths = {r.count("|") for r in rows}
            self.assertEqual(len(widths), 1, rows[0])
            self.assertRegex(rows[1], r"^\|(---\|)+$")

    def test_results_give_every_comparator_on_both_splits(self):
        results = self.text.split("\n## Results\n", 1)[1].split("\n## ", 1)[0]
        for split in ("Lineage split:", "Random split:"):
            block = results.split(split, 1)[1].split("\n\n", 2)[1]
            names = [x.split(" | ")[0][2:] for x in block.strip("\n").split("\n")[2:]]
            self.assertEqual(names, [builder.table_name(w) for w in builder.ROWS], split)

    def test_what_the_review_narrowed(self):
        review = self.text.split("\n## Review\n", 1)[1].split("\n## ", 1)[0]
        self.assertIn("It supports \"no gain established\" for this one network on these 24 datasets, and no wider reading.",
                      review)
        for word in ("Comparator", "Site's best", "Reach", "First look", "Holdout and wording"):
            self.assertIn("- **%s.**" % word, review)
        self.assertIn("goes beyond the figures", review)
        self.assertIn("so they are exploratory", review)

    def test_what_it_cannot_show(self):
        limits = self.text.split("\n## Limits\n", 1)[1].split("\n## ", 1)[0]
        self.assertIn("No transformer, sequence model, other architecture or other setting was fitted", limits)
        self.assertIn("The datasets are not a sample of the release", limits)
        self.assertIn("so this is not external validation", limits)
        self.assertIn("The review did not fit the network again.", limits)


if __name__ == "__main__":
    unittest.main()
