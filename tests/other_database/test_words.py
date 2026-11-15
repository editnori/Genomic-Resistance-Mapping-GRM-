"""The words of the page and of the record's README: what the test is called, and words it never uses."""

import re
import unittest

from tests.other_database.common import PAGE, RECORD

BANNED = r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|carrier|provenance|receipt|folder|safe|scope|notice|under review|transfer)"
# "The reports work on genomes from another database", and the like.
WORKS = r"\bwork(s|ed)?\s+(\w+\s+){0,3}?(data|genomes?)\b"
# Labels of the reviews' findings and notes, which a public reader cannot look up.
LABEL = r"\b(S|N|SF)\d{1,2}\b"
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
            self.assertIsNone(re.search(WORKS, text, re.I), name)

    def test_no_lane_and_no_review_label(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\blanes?\b", prose(text), re.I), name)
            # Funnel steps (C2, D3, F1, O1, P2, Q1, Q2) are named in the page; review labels never are.
            self.assertIsNone(re.search(LABEL, prose(text)), name)

    def test_owner_only_in_the_one_quotation(self):
        self.assertNotIn("owner", self.texts["page"].lower())
        readme = self.texts["README"]
        self.assertEqual(readme.lower().count("owner"), 1)
        self.assertEqual(readme.count('"the owner\'s note"'), 1)

    def test_holds_only_for_the_rule(self):
        # What a BioProject or a dataset has, it contains; "holds" is the rule's word. The kept sentence of the review keeps its own.
        kept = "most scored genomes come from a BioProject that also holds training genomes. It is not external validation."
        for name, text in self.texts.items():
            body = prose(text).replace(kept, "")
            self.assertIsNone(re.search(r"\b(BioProjects?|datasets?|record|which|that|it|they)( also)? holds? (a |\d|training|few|no |every)", body), name)
            self.assertIsNone(re.search(r"\bhold (a median|\d+ BioProjects)", body), name)

    def test_external_validation_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"external validation", s, re.I):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 1)

    def test_called_by_its_full_name(self):
        # ADDENDUM2.md B4's words, in the opening and in the first limit.
        self.assertGreaterEqual(self.texts["page"].count("a test on data from another public database that shares isolates with the training database"), 2)

    def test_calls_are_an_adapted_research_screen(self):
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

    def test_no_file_is_the_subject_of_a_verb_of_speech(self):
        speech = r"(says?|speaks?|thinks?|tells?|claims?|argues?|believes?|states?|explains?|mentions?|notes?)"
        for name, text in self.texts.items():
            self.assertIsNone(re.search(rf"\b(page|file|table|test|record|review|suite|README|addendum|plan)s? {speech}\b", text, re.I), name)
            self.assertIsNone(re.search(rf"\.(md|tsv|py)`? {speech}\b", text), name)

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
