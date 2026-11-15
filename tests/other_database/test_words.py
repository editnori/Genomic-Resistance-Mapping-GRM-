"""The words of the page and of the record's README: what the test is called, and words it never uses."""

import re
import unittest

from tests.other_database.common import PAGE, RECORD

BANNED = r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|carrier|provenance|receipt|folder|safe)"
LEARNER = r"\b(SCM|CART|lasso|Lasso)\b"


def prose(text):
    """Text outside code blocks and inline code, with emphasis marks dropped."""
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"`[^`]*`", "CODE", text)
    return text.replace("*", "")


def sentences(text):
    """Sentences of the prose lines and list items; a stop after one capital letter (E. coli) does not end one."""
    out = []
    for line in prose(text).split("\n"):
        line = line.strip()
        if not line or line.startswith(("|", "#")):
            continue
        line = line[2:] if line.startswith("- ") else line
        out += [s for s in re.split(r"(?<=[^A-Z][.!?])\s+", line) if s]
    return out


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.texts = {"page": PAGE.read_text(encoding="utf-8"), "README": (RECORD / "README.md").read_text(encoding="utf-8")}

    def test_banned_words(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BANNED, text, re.I), name)

    def test_external_validation_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"external validation", s, re.I):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 1)

    def test_called_a_test_on_data_from_another_public_database(self):
        self.assertIn("a test on data from another public database", self.texts["page"])
        self.assertIn("an adapted research screen", self.texts["page"])

    def test_carries_over_never(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\bcarr(y|ies|ied) over\b", text, re.I), name)

    def test_marks(self):
        for name, text in self.texts.items():
            self.assertNotIn("—", text, name)
            self.assertIsNone(re.search(r"\brather than\b", text, re.I), name)

    def test_no_sentence_ends_in_not_x(self):
        for name, text in self.texts.items():
            for s in sentences(text):
                self.assertIsNone(re.search(r",\s*not\s+[^,]*$", s.rstrip(".:;")), s)

    def test_pages_give_list_show_or_hold(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\b(page|file|table|test|record|review|suite)s? (says?|speaks?|thinks?|tells?|claims?|argues?|believes?)\b", text, re.I), name)

    def test_no_learner_ranked_above_another(self):
        for s in sentences(self.texts["page"]):
            if len(set(m.lower() for m in re.findall(LEARNER, s))) >= 2:
                self.assertIsNone(re.search(r"\b(better|best|worse than|outperform\w*|superior|beats?|stronger|weaker|wins?)\b", s, re.I), s)

    def test_names_no_sequence(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text), name)
            self.assertIsNone(re.search(r"\bk-?mers?\b", text, re.I), name)


if __name__ == "__main__":
    unittest.main()
