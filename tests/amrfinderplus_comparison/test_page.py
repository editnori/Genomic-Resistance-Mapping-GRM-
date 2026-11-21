"""docs/amrfinderplus-comparison.md against the committed record: the limit and catalog sentences, the section
Reading, the paragraphs that give the record's figures, every table, every other paragraph and list item word for
word, and every other number of the page."""

import re
import unittest

from tests.amrfinderplus_comparison import lines
from tests.amrfinderplus_comparison.common import PAGE, RECORD, read_tsv
from tests.amrfinderplus_comparison.held import HELD

TITLE = "# AMRFinderPlus comparison"
LIMIT = ("This is a retrospective comparison inside one public database, BV-BRC. "
         "It is not a prospective test and not external validation.")
SECTIONS = ["Reading", "Plan and order", "Method", "The program against the known-gene rule",
            "Learners against the program", "Quinolone and carbapenem datasets", "Datasets kept apart", "Review",
            "What the record cannot show", "Files"]
TABLE_SECTION = {"order": "Plan and order", "program": "The program against the known-gene rule",
                 "learners": "Learners against the program", "random": "Learners against the program",
                 "sets": "Quinolone and carbapenem datasets", "apart": "Datasets kept apart"}
# The review's label for every random-split figure, held here apart from lines.py so that a builder cannot drop it.
RANDOM_LABEL = "random split: training and test share lineages, and the served learner is picked on this split"
BUILT_SECTION = [("program_paragraph", "The program against the known-gene rule"),
                 ("learner_paragraphs", "Learners against the program"),
                 ("rule_m_paragraph", "Quinolone and carbapenem datasets"),
                 ("apart_paragraphs", "Datasets kept apart")]


def outside_code(text):
    out, fenced = [], False
    for line in text.split("\n"):
        if line.startswith("```"):
            fenced = not fenced
            continue
        if not fenced:
            out.append(line)
    return out


def blocks(text):
    """The page's paragraphs and list items outside code, tables and headers, one string each, in page order."""
    out = []
    for para in "\n".join(outside_code(text)).split("\n\n"):
        para = para.strip("\n")
        if not para or para.startswith(("#", "|")):
            continue
        if para.startswith("- "):
            out += ["- " + x for x in para[2:].split("\n- ")]
        else:
            out.append(para)
    return out


def section(text, title):
    return text.split(f"\n## {title}\n", 1)[1].split("\n## ", 1)[0]


def page_tables(text):
    tables, cur = [], []
    for line in outside_code(text) + [""]:
        if line.startswith("|"):
            cur.append(line)
        elif cur:
            tables.append(cur)
            cur = []
    return tables


def sentences(text):
    return [s for b in blocks(text) for s in re.split(r"(?<=[a-z0-9)][.])\s+(?=[A-Z0-9])", b)]


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.rec = lines.Record()
        cls.summary = {(r["set"], r["contrast"], r["measure"]): r for r in read_tsv(RECORD / "summary.tsv")}

    def test_first_paragraph(self):
        paragraphs = [p for p in self.text.split("\n\n") if p.strip()]
        self.assertEqual(paragraphs[0], TITLE)
        first = paragraphs[1]
        self.assertIn(LIMIT, first)
        self.assertIn(lines.CATALOG, first)
        self.assertEqual(self.text.count(LIMIT), 1)

    def test_reading_carries_the_catalog_sentence(self):
        self.assertIn(lines.CATALOG, section(self.text, "Reading"))
        self.assertIn(lines.MARKER, section(self.text, "Reading"))

    def test_sections_in_order(self):
        body = outside_code(self.text)
        self.assertEqual([x[3:] for x in body if x.startswith("## ")], SECTIONS)
        self.assertEqual([x for x in body if x.startswith("# ")], [TITLE])
        self.assertFalse([x for x in body if x.startswith("###")])

    def test_reading_is_rebuilt_from_the_record(self):
        body = section(self.text, "Reading").strip()
        self.assertEqual(body.split("\n\n"), lines.reading(self.rec))

    def test_reading_has_five_sentences(self):
        body = section(self.text, "Reading").strip()
        found = [x for p in body.split("\n\n") for x in re.split(r"(?<=[a-z)][.])\s+(?=[A-Z])", p)]
        self.assertEqual(len(found), 5)
        for x in found:
            self.assertTrue(x.endswith("."), x)

    def test_reading_figures_stand_in_the_summary(self):
        """Read back from the page's words, apart from the builder in lines.py."""
        body = section(self.text, "Reading")
        m = re.search(r"was ahead of GRM's known-gene rule in (\d+) of (\d+) datasets, level in (\d+) and behind in "
                      r"(\d+), with a median difference of ([+-]\d\.\d{4}) and a mean of ([+-]\d\.\d{4})", body)
        self.assertIsNotNone(m)
        r = self.summary[(lines.NEW, lines.PR, "ba")]
        self.assertEqual(m.groups()[:4], (r["ahead"], r["datasets"], r["level"], r["behind"]))
        self.assertEqual(m.groups()[4:], ("%+.4f" % float(r["median"]), "%+.4f" % float(r["mean"])))
        self.assertLess(body.index(" ahead of GRM's"), body.index("median difference"))
        self.assertLess(body.index("median difference"), body.index("a mean of"))
        third = [p for p in body.split("\n\n") if p.startswith("On the lineage split")][0]
        named = re.findall(r"([+-]\d\.\d{4}) for (scm|cart|rf|gbt|lasso)", third)
        self.assertEqual([l for _, l in named], list(lines.LEARNERS))
        for value, l in named:
            self.assertEqual(value, "%+.4f" % float(self.summary[(lines.NEW, lines.lin(l), "ba")]["median"]), l)
        for word, m_ in (("sensitivity", "sens"), ("specificity", "spec")):
            found = re.search(r"median %s ((?:[+-]\d\.\d{4}(?:, | and )?){5})" % word, third).group(1)
            got = re.findall(r"[+-]\d\.\d{4}", found)
            want = ["%+.4f" % float(self.summary[(lines.NEW, lines.lin(l), m_)]["median"]) for l in lines.LEARNERS]
            self.assertEqual(got, want, word)

    def test_built_paragraphs_stand_in_their_sections(self):
        for name, title in BUILT_SECTION:
            built = getattr(lines, name)(self.rec)
            found = blocks(section(self.text, title))
            for para in built:
                self.assertEqual(found.count(para), 1, para[:60])

    def test_every_paragraph_and_list_item_is_built_or_held(self):
        known = set(lines.built(self.rec)) | set(HELD)
        self.assertEqual([b for b in blocks(self.text) if b not in known], [])

    def test_held_lines_stand_once_and_in_order(self):
        self.assertEqual(len(HELD), len(set(HELD)))
        found = blocks(self.text)
        for line in HELD:
            self.assertEqual(found.count(line), 1, line)
        self.assertEqual([b for b in found if b in set(HELD)], HELD)
        self.assertFalse(set(lines.built(self.rec)) & set(HELD), "a held line is also rebuilt")

    def test_every_table_is_built_from_the_record(self):
        built = lines.tables(self.rec)
        found = page_tables(self.text)
        self.assertEqual(len(found), len(built))
        for name, table in built.items():
            self.assertEqual(sum(1 for t in found if t == table), 1, name)
            self.assertIn("\n" + "\n".join(table) + "\n", section(self.text, TABLE_SECTION[name]), name)

    def test_every_number_is_in_the_record(self):
        found, _ = lines.values(RECORD)
        missing = [(x, line[:50]) for line in outside_code(self.text) for x in lines.numbers(line) if x not in found]
        self.assertEqual(missing, [])

    def test_hashes_are_the_records(self):
        _, hashes = lines.values(RECORD)
        found = re.findall(r"`([0-9a-f]{6,64})`", self.text)
        self.assertGreaterEqual(len(found), 4)
        for h in found:
            self.assertTrue(len(h) >= 8 and any(x.startswith(h) for x in hashes), h)

    def test_one_sided_count_only_with_its_definition(self):
        r = self.summary[(lines.NEW, lines.lin("served_best"), "ba")]
        one, two, n = r["learner_at_most_0.02_ahead"], r["within_0.02_either_way"], r["datasets"]
        with_one = [s for s in sentences(self.text) if f"{one} of {n}" in s]
        self.assertTrue(with_one)
        for s in with_one:
            self.assertIn("at most 0.02 ahead of the program, or behind it", s)
            self.assertIn("as rule M defines it", s)
        self.assertIn(f"the two were within 0.02 of each other in {two} of {n} datasets", self.text)
        for s in sentences(self.text):
            if re.search(r"within 0\.02", s) and "either way" not in s and "of each other" not in s:
                self.assertRegex(s, r"at most 0\.02 ahead|as rule M defines it", s)

    def test_served_learner_and_random_split_named_as_the_review_names_them(self):
        self.assertIn("the learner the release serves, picked on the random split", self.text)
        self.assertEqual(lines.RANDOM_LABEL, RANDOM_LABEL)
        labelled = 0
        for line in outside_code(self.text):
            if line.startswith("|") and "random split" in line:
                self.assertIn(RANDOM_LABEL, line)
                labelled += 1
        self.assertGreaterEqual(labelled, 7)
        self.assertIn("training and test share lineages, so it is no test on unseen lineages, and the served learner "
                      "is picked on it", self.text)

    def test_what_the_record_cannot_show(self):
        body = section(self.text, "What the record cannot show")
        st = self.rec.settings
        self.assertIn(f"The plan's sensitivity run on the newest catalog database ({st['newest_database_not_run']}) "
                      f"was not run, so these figures hold for database {st['database_version']} only.", body)
        self.assertIn("The program's tables per genome are gone, so its figures can be checked only at the counts per "
                      "dataset and per fold.", body)
        self.assertIn("The species interval depends on a generator and a draw the plan did not name.", body)
        method = section(self.text, "Method")
        self.assertIn("named neither the generator nor the draw", method)
        self.assertIn("`random.Random`", method)
        self.assertIn("`randrange` over the species sorted by name", method)

    def test_order_follows_the_record(self):
        self.assertEqual([r["event"] for r in self.rec.order], list(lines.ORDER_LABELS))
        times = [r["utc"] for r in self.rec.order]
        self.assertEqual(times, sorted(times))

    def test_file_count(self):
        count = len(list(RECORD.iterdir()))
        self.assertIn(f"holds the record: {count} files that", self.text)

    def test_page_names_its_commands(self):
        self.assertIn("python3 -I scripts/amrfinderplus-comparison/record.py --check "
                      "reference/amrfinderplus-comparison\n", self.text)
        self.assertIn("python3 -I -m unittest discover -s tests/amrfinderplus_comparison -t .\n", self.text)


if __name__ == "__main__":
    unittest.main()
