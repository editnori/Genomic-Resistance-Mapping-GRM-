"""docs/time-split-stagger.md and the words of the staggered split's scripts: the frozen plan's sha256, "Not run
yet.", bare headers, no result, the words the page never uses, its row in docs/INDEX.md and the workflow's comment."""

import re
import unittest

from tests.time_split import stagger_fixture as fx
from tests.time_split.common import WORKFLOW
from tests.time_split.test_words import BANNED, LABEL, MAKER, NAMES, ROLE, prose, sentences

sd = fx.design_script()
# Words sent back by the reader the page is written for, and words likely to lose them.
# Written so that this file holds none of them.
SENT_BACK = (r"\b(held|spot[- ]checked|la[n]es?|census|receipts?|provenance|landed|cleanly|traps?|snags?|breakthroughs?|"
             r"synced|frames?|keepers?|guarded|folders?|own[e]rs?)\b")
SCRIPTS = (fx.DESIGN, fx.ANALYZE, fx.FAULTS)


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = fx.PAGE.read_text(encoding="utf-8")
        cls.lines = cls.page.rstrip("\n").split("\n")

    def test_not_run_yet(self):
        self.assertEqual(self.page.count("Not run yet."), 1)
        self.assertTrue(self.lines[2].startswith("Not run yet."))

    def test_the_record_line_quotes_the_plans_sha256(self):
        self.assertTrue(self.lines[-1].startswith("Record: "))
        self.assertIn(f"sha256 {sd.PLAN_SHA256}", self.lines[-1])
        self.assertEqual(self.page.count(sd.PLAN_SHA256), 1)

    def test_paths_and_hashes_only_in_the_record_line(self):
        for line in self.lines[:-1]:
            self.assertIsNone(re.search(r"\.(py|md|tsv|json|txt)\b|\b[0-9a-f]{16,}\b|/", line), line)

    def test_headers_are_bare_nouns(self):
        headers = re.findall(r"^#+ (.+)$", self.page, re.M)
        self.assertEqual(headers[0], "Staggered split by time")
        for h in headers:
            self.assertIsNone(re.search(r"[:?!.,;()0-9]", h), h)
            self.assertIsNone(re.match(r"(?i)(the|a|an|each|every|how|why|what) ", h), h)
            self.assertLessEqual(len(h.split()), 4, h)

    def test_no_result(self):
        self.assertIsNone(re.search(r"[+-]\d*\.\d+|\[\s*[+-]?\d", self.page))
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(was|were|found|came out|turned out)\b", s), s)

    def test_lengths(self):
        for block in re.split(r"\n\s*\n", self.page):
            for item in re.split(r"\n(?=- |\d+\. )", block):
                if item.startswith(("- ", "1. ", "2. ", "3. ")):
                    self.assertLess(len(item.split()), 45, item)
                elif not item.startswith(("#", "|", "Record: ")):
                    self.assertLess(len(item.split()), 85, item)

    def test_words_it_never_uses(self):
        text = prose(self.page)
        for pattern in (BANNED, SENT_BACK, ROLE, NAMES, MAKER):
            self.assertIsNone(re.search(pattern, self.page, re.I), pattern)
        self.assertIsNone(re.search(LABEL, text))
        self.assertIsNone(re.search(r"\bk-?mers?\b|\b(genes?|alleles?|mutations?|variants?|markers?|determinants?|"
                                    r"snps?)\b", self.page, re.I))

    def test_not_shown_is_never_equal(self):
        self.assertIn("not shown", self.page)
        self.assertIsNone(re.search(r"(?i)\b(equal\w*|the same score|no difference)\b", self.page))

    def test_prospective_and_external_validation_only_denied(self):
        found = 0
        for s in sentences(self.page):
            for m in re.finditer(r"(?i)prospective|external validation", s):
                found += 1
                self.assertRegex(s[:m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertEqual(found, 2)
        self.assertIn("retrospective split inside one database", self.page)

    def test_in_the_index(self):
        index = (fx.ROOT / "docs" / "INDEX.md").read_text(encoding="utf-8")
        self.assertEqual(index.count("[time-split-stagger.md](time-split-stagger.md)"), 1)

    def test_the_workflows_comment_names_the_page(self):
        text = WORKFLOW.read_text(encoding="utf-8")
        comment = text.split("\n          - suite: tests/time_split\n", 1)[0].rsplit("\n", 1)[1]
        self.assertTrue(comment.strip().startswith("# SYNTHETIC"), comment)
        self.assertIn("docs/time-split-stagger.md", comment)


class ScriptWords(unittest.TestCase):
    def test_marks_and_names(self):
        for path in SCRIPTS + (fx.PAGE,):
            text = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                self.assertNotIn(chr(0x2014), text)
                self.assertIsNone(re.search(r"\brather than\b", text, re.I))
                for pattern in (ROLE, NAMES, MAKER):
                    self.assertIsNone(re.search(pattern, text), pattern)
                self.assertIsNone(re.search(LABEL, text))

    def test_no_id_accession_or_sequence(self):
        for path in SCRIPTS + (fx.PAGE,):
            text = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                self.assertIsNone(re.search(r"\bPRJ[A-Z]{1,2}\d+|\bSAM[NED][A-Z]?\d+|\bGC[AF]_\d+|\b[SE]R[RSX]\d{5,}",
                                            text))
                self.assertIsNone(re.search(r"\b\d{3,7}\.\d{3,6}\b", text))
                self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text))

    def test_the_scripts_use_the_standard_library_alone(self):
        allowed = {"contextlib", "datetime", "hashlib", "importlib", "io", "json", "math", "os", "random", "re",
                   "statistics", "subprocess", "sys", "pathlib"}
        for path in SCRIPTS:
            source = path.read_text(encoding="utf-8")
            with self.subTest(file=path.name):
                self.assertLessEqual(set(re.findall(r"^(?:import|from) ([a-zA-Z_]+)", source, re.M)), allowed)
                self.assertNotIn("urllib", source)
                self.assertNotIn("socket", source)


if __name__ == "__main__":
    unittest.main()
