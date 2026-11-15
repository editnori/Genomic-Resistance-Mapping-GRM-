"""docs/other-database-test.md: every line that carries a number is rebuilt from the record, once; the number-free limits and definitions are held word for word; and the page's form."""

import re
import unittest

from tests.other_database import lines as builder
from tests.other_database.common import PAGE, RECORD, ROOT, load_script, read_tsv

record = load_script()


def blocks(text):
    """The page as (kind, lines): header, table, list, code or paragraph."""
    out, cur, kind, fence = [], [], None, False
    for line in text.split("\n"):
        if line.startswith("```"):
            if fence:
                out.append(("code", cur + [line]))
                cur, kind, fence = [], None, False
            else:
                if cur:
                    out.append((kind, cur))
                cur, kind, fence = [line], "code", True
            continue
        if fence:
            cur.append(line)
            continue
        if not line.strip():
            if cur:
                out.append((kind, cur))
            cur, kind = [], None
            continue
        k = "header" if line.startswith("#") else "table" if line.startswith("|") else "list" if line.startswith("- ") else "paragraph"
        if kind not in (None, k) or k == "header":
            if cur:
                out.append((kind, cur))
            cur = []
        cur.append(line)
        kind = k
    if cur:
        out.append((kind, cur))
    return out


def words(text):
    return len(text.split())


def cells(row):
    return [c.strip() for c in row.strip().strip("|").split("|")]


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.lines = cls.text.split("\n")
        cls.blocks = blocks(cls.text)
        cls.built = builder.lines()
        cls.held = builder.held()

    def test_every_numbered_line_is_built(self):
        built = set(self.built)
        for line in self.lines:
            if builder.carries_number(line) and not re.fullmatch(r"\|(---\|)+", line):
                self.assertIn(line, built)

    def test_every_built_line_is_on_the_page_once(self):
        numbered = [x for x in self.built if builder.carries_number(x)]
        self.assertEqual(len(numbered), len(set(numbered)), "the builder gives a line twice")
        for line in numbered:
            self.assertEqual(self.lines.count(line), 1, line)
        for line in self.built:
            self.assertIn(line, self.lines)

    def test_built_lines_in_page_order(self):
        # Each built line after the one before it; a table's header row may stand on the page more than once.
        at = -1
        for x in self.built:
            self.assertIn(x, self.lines[at + 1:], x)
            at = self.lines.index(x, at + 1)

    def test_number_free_lines_held_word_for_word(self):
        # The opening's denial of external validation, the Plan opener, the rule's limit, the reading of the kept sentences,
        # the definitions of readings L, C and D and the number-free limits: each on the page once, as lines.py words it.
        self.assertEqual(len(self.held), len(set(self.held)))
        for line in self.held:
            self.assertEqual(self.lines.count(line), 1, line)
            self.assertFalse(builder.carries_number(line), line)
        self.assertTrue(any("It is not external validation." in x and x.startswith("It is " + builder.NAME) for x in self.held))
        self.assertTrue(any(x.startswith("Reading L keeps") and "D resamples whole BioProjects" in x for x in self.held))

    def test_limits_whole(self):
        # The Limits list: its count and each item, in order, and nothing more.
        record_ = builder.Record()
        at = self.lines.index("## Limits")
        block = next(ls for kind, ls in blocks("\n".join(self.lines[at + 1:])) if kind == "list")
        self.assertEqual(block, builder.limits(record_))
        self.assertEqual(len(block), 10)

    def test_headers(self):
        heads = [b[1][0] for b in self.blocks if b[0] == "header"]
        self.assertLessEqual(len(heads), 12)
        self.assertEqual(sum(h.startswith("# ") for h in heads), 1)
        self.assertTrue(self.lines[0].startswith("# "))
        for h in heads:
            title = h.lstrip("#").strip()
            self.assertIn(words(title), (1, 2), h)
            self.assertNotRegex(title, r"^(What|How|Why)\b")
            self.assertFalse(builder.carries_number(title), h)

    def test_lists(self):
        items = [x for b in self.blocks if b[0] == "list" for x in b[1]]
        self.assertLessEqual(len(items), 30)
        for x in items:
            self.assertLessEqual(words(x[2:]), 45, x)

    def test_paragraphs(self):
        for kind, ls in self.blocks:
            if kind == "paragraph":
                self.assertEqual(len(ls), 1, "a paragraph is one line")
                self.assertLessEqual(words(ls[0]), 80, ls[0])

    def test_tables(self):
        for kind, ls in self.blocks:
            if kind != "table":
                continue
            head = cells(ls[0])
            self.assertRegex(ls[1], r"^\|(---\|)+$")
            self.assertEqual(ls[1].count("|") - 1, len(head))
            for h in head:
                self.assertIn(words(h), (1, 2), h)
            for row in ls[2:]:
                self.assertEqual(len(cells(row)), len(head), row)
                for c in cells(row):
                    self.assertLessEqual(words(c), 30, c)

    def test_recall_log_loss_and_specificity_side_by_side(self):
        seen = 0
        for kind, ls in self.blocks:
            head = cells(ls[0]) if kind == "table" else []
            if "Recall" in head:
                seen += 1
                at = sorted(head.index(x) for x in ("Recall", "Log loss", "Specificity"))
                self.assertEqual(at, list(range(at[0], at[0] + 3)), ls[0])
        self.assertEqual(seen, 4)

    def test_no_line_numbers(self):
        self.assertNotRegex(self.text, r"\.(py|rs|md|yml|tsv):\d")
        self.assertNotRegex(self.text, r"\bline \d")

    def test_number_free_claims(self):
        """Sentences that carry no number but state a fact of the record, rebuilt from it."""
        species = {r["species"]: r for r in read_tsv(RECORD / "species.tsv")}
        # The funnel steps that removed no isolate in any species, Q2's two parts as one step.
        zero = []
        for code, _ in record.STEPS:
            step = "Q2" if code.startswith("Q2") else code
            if all(r[code] == "0" for r in species.values()) and step not in zero:
                zero.append(step)
            elif any(r[code] not in ("0", "") for r in species.values()) and step in zero:
                zero.remove(step)
        self.assertIn(f"steps {builder.joined(zero)} removed none in any species.", self.text)
        # The species where most scored genomes share a BioProject with training genomes.
        most = [sp for sp in [builder.PILOT] + builder.EIGHT if 2 * int(species[sp]["shared_bioproject"]) > int(species[sp]["scored"])]
        self.assertIn(f"- In {builder.joined(builder.short(sp) for sp in most)} most scored genomes share a BioProject with training genomes.", self.lines)
        # The status of the readings: the pilot's chosen after its results, the eight's fixed by ADDENDUM2.md.
        status = {part: {r["status"] for r in read_tsv(RECORD / "readings.tsv") if r["part"] == part} for part in ("pilot", "eight")}
        self.assertEqual(status, {"pilot": {"exploratory"}, "eight": {"fixed in advance"}})
        self.assertIn("The pilot's readings other than the plan's rule are exploratory: the pilot's review chose them after its results were seen.", self.text)
        self.assertIn("Each reading of the pilot other than the plan's rule is exploratory.", self.text)
        self.assertIn("`ADDENDUM2.md` fixed readings L, C and D before any of their labels was read.", self.text)
        # "SCM, CART and the lasso are scored; RF and GBT are not"
        settings = {r["key"]: r["value"] for r in read_tsv(RECORD / "settings.tsv")}
        self.assertEqual(settings["learners_scored"], "cart;lasso;scm")
        self.assertTrue(settings["learners_not_scored"].startswith("gbt;rf"))
        learners = {r["learner"] for r in read_tsv(RECORD / "datasets.tsv")}
        self.assertEqual(learners, {"scm", "cart", "lasso"})

    def test_files_named_exist(self):
        for path in ("reference/other-database", "scripts/other-database/record.py", "tests/other_database"):
            self.assertIn(f"`{path}", self.text)
            self.assertTrue((ROOT / path).exists(), path)

    def test_commands_match_the_script(self):
        code = "\n".join(x for kind, ls in self.blocks if kind == "code" for x in ls)
        for flag in re.findall(r"(--[a-z-]+)", code):
            self.assertIn(f'"{flag}"', (ROOT / "scripts" / "other-database" / "record.py").read_text(encoding="utf-8"), flag)
        self.assertIn("discover -s tests/other_database -t .", code)


if __name__ == "__main__":
    unittest.main()
