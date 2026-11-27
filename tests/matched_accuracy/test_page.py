"""docs/matched-accuracy.md and the words of stage 1's files: "Not run yet.", bare headers, no result, stage 2 not
authorised, the words the page never uses, its row in docs/INDEX.md and the workflow's comment."""

import re
import unittest

from tests.matched_accuracy import common
from tests.matched_accuracy.test_reference import AMENDMENT_SHA256, PROTOCOL_SHA256

# Words the page never uses, and words sent back by the reader it is written for. Written so that this file
# holds none of them.
BANNED = (r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|provenance|receipts?|manifest|"
          r"fold[e]rs?|safe|help|explain|becau[s]e|cau[s]e|due to|clinic|patient|sc[o]pe|threshold|mechanis|confer|"
          r"proves?|pro[o]f|guarantee|h[e]ld|spot[- ]checked|la[n]es?|census|landed|lands|cleanly|traps?|snags?|"
          r"breakthroughs?|synced|frames?|keepers?|guarded|own[e]rs?)\b")
ROLE = r"(?i)\b(mach[i]nes?|b[o]x(es)?|contr[o]ller|threa[d]s?|work[e]rs?|assist[a]nt|mod[e]ls?|vend[o]rs?|agen[t]s?)\b"
MAKER = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|g[p]t-?\d|gem[i]ni)\b"
MARKER = r"(?i)\bk-?mers?\b|\b(genes?|alleles?|mutations?|variants?|markers?|determinants?|snps?)\b"
# The numbers the page may give: the revision, the checks' figures, the steps and the stages. None is a result.
NUMBERS = {"14", "1", "2", "3", "4", "5", "6", "7", "19:10:00", "2026", "1,285,111", "252,148"}


def prose(text):
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    return re.sub(r"`[^`]*`", "CODE", text).replace("*", "")


def sentences(text):
    out = []
    for line in prose(text).split("\n"):
        line = line.strip()
        if not line or line.startswith(("|", "#", "Record: ")):
            continue
        line = re.sub(r"^(- |\d+\. )", "", line)
        out += [s for s in re.split(r"(?<=[^A-Z][.!?])\s+", line) if s]
    return out


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = common.PAGE.read_text(encoding="utf-8")
        cls.lines = cls.page.rstrip("\n").split("\n")
        cls.body = "\n".join(cls.lines[:-1])

    def test_not_run_yet(self):
        self.assertEqual(self.page.count("Not run yet."), 1)
        self.assertTrue(self.lines[2].startswith("Not run yet."))

    def test_under_500_words(self):
        self.assertLess(len(self.page.split()), 500)

    def test_the_record_line_gives_the_protocols_and_the_amendments_sha256(self):
        self.assertTrue(self.lines[-1].startswith("Record: "))
        for sha in (PROTOCOL_SHA256, AMENDMENT_SHA256):
            self.assertIn(f"sha256 {sha}", self.lines[-1])
            self.assertEqual(self.page.count(sha), 1)

    def test_paths_and_hashes_only_in_the_record_line(self):
        for line in self.lines[:-1]:
            self.assertIsNone(re.search(r"\.(py|md|tsv|json|jsonl|txt|gz)\b|\b[0-9a-f]{16,}\b|/", line), line)

    def test_headers_are_bare_nouns(self):
        headers = re.findall(r"^#+ (.+)$", self.page, re.M)
        self.assertEqual(headers[0], "Matched accuracy")
        self.assertLessEqual(len(headers), 12)
        for h in headers:
            self.assertIsNone(re.search(r"[:?!.,;()0-9]", h), h)
            self.assertIsNone(re.match(r"(?i)(the|a|an|each|every|how|why|what) ", h), h)
            self.assertLessEqual(len(h.split()), 4, h)

    def test_lengths(self):
        items = 0
        for block in re.split(r"\n\s*\n", self.body):
            for item in re.split(r"\n(?=- |\d+\. )", block):
                if re.match(r"(- |\d+\. )", item):
                    items += 1
                    self.assertLess(len(item.split()), 45, item)
                elif not item.startswith(("#", "|")):
                    self.assertLess(len(item.split()), 85, item)
        self.assertLessEqual(items, 30)

    def test_no_result(self):
        found = set(re.findall(r"(?<![A-Za-z\d])(?:\d[\d,:]*\d|\d)(?![A-Za-z\d])", self.body))
        self.assertLessEqual(found, NUMBERS)
        self.assertIsNone(re.search(r"\d\.\d|%", self.body))
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(was|were|found|came out|turned out|shows?|showed)\b", s), s)

    def test_stage_2_is_not_authorised(self):
        said = [s for s in sentences(self.page) if re.search(r"(?i)\bstage 2\b", s)]
        self.assertEqual(len(said), 2)
        for s in said:
            self.assertIn("is not authorised", s)

    def test_retrospective_inside_one_database_and_external_validation_only_denied(self):
        found = 0
        for s in sentences(self.page):
            for m in re.finditer(r"(?i)prospective|external validation", s):
                if s[m.start() - 2:m.start()].lower() == "re":
                    continue
                found += 1
                self.assertRegex(s[:m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertEqual(found, 2)
        self.assertIn("retrospective read inside one database", self.page)

    def test_words_it_never_uses(self):
        for pattern in (BANNED, ROLE, MAKER, MARKER):
            self.assertIsNone(re.search(pattern, self.page, re.I), pattern)
        self.assertNotIn(chr(0x2014), self.page)
        self.assertIsNone(re.search(r"(?i)\brather than\b", self.page))

    def test_in_the_index(self):
        index = common.INDEX.read_text(encoding="utf-8")
        self.assertEqual(index.count("[matched-accuracy.md](matched-accuracy.md)"), 1)
        row = next(line for line in index.split("\n") if line.startswith("| [matched-accuracy.md]"))
        self.assertIn("stage 2 is not authorised. Not run yet.", row)

    def test_the_workflows_comment_names_the_page(self):
        text = common.WORKFLOW.read_text(encoding="utf-8")
        block = text.split("\n          - suite: tests/matched_accuracy\n", 1)[0].rsplit("\n          - suite:", 1)[1]
        comments = [line.strip() for line in block.split("\n") if line.strip().startswith("#")]
        self.assertTrue(comments[0].startswith("# SYNTHETIC"), comments)
        self.assertIn("docs/matched-accuracy.md", comments[0])


class FileWords(unittest.TestCase):
    """The words of every file of stage 1: no maker, no long dash, and no real id, accession or sequence."""

    @classmethod
    def setUpClass(cls):
        paths = [common.SCRIPT, common.FAULTS, common.PAGE] + sorted((common.ROOT / common.SUITE).glob("*.py"))
        cls.texts = {p.relative_to(common.ROOT).as_posix(): p.read_text(encoding="utf-8") for p in paths}

    def test_no_maker_and_no_long_dash(self):
        for name, text in self.texts.items():
            with self.subTest(file=name):
                self.assertIsNone(re.search(MAKER, text))
                self.assertNotIn(chr(0x2014), text)
                self.assertIsNone(re.search(r"(?i)\brather than\b", text))

    def test_no_id_accession_or_sequence(self):
        for name, text in self.texts.items():
            with self.subTest(file=name):
                self.assertIsNone(re.search(r"(?i)(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_|[SED]R[RSXP])\d+", text))
                self.assertIsNone(re.search(r"\b\d{3,7}\.\d{3,6}\b", text))
                self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text))

    def test_no_name_of_a_marker_in_the_script_and_its_runner(self):
        for name in ("scripts/matched-accuracy/stage1.py", "scripts/matched-accuracy/planted_faults.py"):
            self.assertIsNone(re.search(MARKER, self.texts[name]), name)


if __name__ == "__main__":
    unittest.main()
