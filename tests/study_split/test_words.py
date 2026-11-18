"""The words of the page, the record's README, the script and this suite: what the test is called, and words they never use."""

import re
import unittest

from tests.study_split.common import PAGE, RECORD, ROOT, SCRIPT, SUITE

BANNED = (r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|carrier|provenance|receipt|folder|"
          r"safe|help|explain|because|caus|due to|no tuning|prospective|clinic|patient|scope|transfer)")
# Labels of the reviews' findings and notes, which a public reader cannot look up.
LABEL = r"\b(S|N|SF|F)\d{1,2}\b"
# The names of the work's parts and of the machines it ran on, written so that this file holds none of them.
NAMES = r"(?i)\b(study(split|draws|doc)|rev[ -]?study(draws)?|box[-]node|boa[t]|fix1\d\d|rev1\d\d)\b"
VENDOR = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t)\b"
LEARNER = r"\b(scm|cart|rf|gbt|lasso)\b"
FILE_VERBS = ("gives", "give", "lists", "list", "shows", "show", "holds", "hold", "is", "are")


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


def code_files():
    return {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8")
            for p in [SCRIPT] + sorted((ROOT / SUITE).glob("*.py"))}


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.texts = {"page": PAGE.read_text(encoding="utf-8"),
                     "README": (RECORD / "README.md").read_text(encoding="utf-8")}
        cls.page = cls.texts["page"]

    def test_banned_words(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BANNED, text, re.I), name)

    def test_no_name_of_the_work_or_label(self):
        for name, text in list(self.texts.items()) + list(code_files().items()):
            self.assertIsNone(re.search(NAMES, text), name)
            self.assertIsNone(re.search(VENDOR, text), name)
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\blanes?\b", prose(text), re.I), name)
            self.assertIsNone(re.search(LABEL, prose(text)), name)
        for name, text in code_files().items():
            self.assertIsNone(re.search(LABEL, text), name)

    def test_external_validation_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"external validation", s, re.I):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 1)
        self.assertIn("a retrospective split by study inside one database", self.page)

    def test_the_owner_s_phrases(self):
        self.assertIn("the part not matched by fold shape", self.page.lower())
        for m in re.finditer(r"in line with", self.page):
            self.assertEqual(self.page[m.start() - 3:m.start()], "is ", self.page[m.start() - 40:m.end()])
        self.assertIn("Tuning keeps copy groups whole and keeps no study whole in any arm; the release's lineage split, "
                      "the first experiment's third arm, tunes on lineage groups.", self.page)
        self.assertIn("The design does not separate study from lineage", self.page)

    def test_no_drop_given_to_study_alone(self):
        # The design does not separate study from lineage: no sentence gives the drop to study membership alone.
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(study|BioProject)( membership)? alone\b", s), s)
            self.assertIsNone(re.search(r"(?i)\b(attribut\w*|owing to|driven by|accounts? for)\b", s), s)

    def test_no_time_split_and_no_stand_in_for_a_lineage_split(self):
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(time split|temporal|by date|deposit\w*)\b", s), s)
            if re.search(r"(?i)\blineage\b", s):
                self.assertIsNone(re.search(r"(?i)\bstands? in for\b", s), s)

    def test_marks(self):
        for name, text in self.texts.items():
            self.assertNotIn("—", text, name)
            self.assertIsNone(re.search(r"\brather than\b", text, re.I), name)

    def test_no_sentence_ends_in_not_x(self):
        for name, text in self.texts.items():
            for s in sentences(text):
                self.assertIsNone(re.search(r",\s*not\s+[^,]*$", s.rstrip(".:;")), s)
            self.assertIsNone(re.search(r"\bIt is [^.]*\. It is not\b", prose(text)), name)

    def test_files_give_list_show_or_hold(self):
        verbs = 0
        for name, text in self.texts.items():
            for m in re.finditer(r"(`[^`]+\.(tsv|md|sha256|py)`|`SHA256SUMS`|\]\([\w.-]+\.md\)|\bREADME)"
                                 r"( in the record)? ([a-z]+)\b", text):
                if m.group(4) in ("and", "in", "of", "from", "beside", "that", "with", "to", "for", "on", "or"):
                    continue
                verbs += 1
                # A script also writes; a file only gives, lists, shows or holds.
                allowed = FILE_VERBS + (("writes",) if m.group(2) == "py" else ())
                self.assertIn(m.group(4), allowed, m.group(0))
        self.assertGreaterEqual(verbs, 6)

    def test_no_file_is_the_subject_of_a_verb_of_speech(self):
        speech = r"(says?|speaks?|thinks?|tells?|claims?|argues?|believes?|states?|explains?|mentions?|notes?|records?)"
        for name, text in self.texts.items():
            self.assertIsNone(re.search(rf"\b(page|file|table|test|record|review|suite|README|plan|script)s? "
                                        rf"{speech}\b", text, re.I), name)
            self.assertIsNone(re.search(rf"\.(md|tsv|py)`? {speech}\b", text), name)

    def test_no_learner_ranked_above_another(self):
        # A ranking word fails beside any learner name, one or more.
        rank = (r"\b(best|worst|worse than|better than|outperform\w*|superior|inferior|beats?|stronger|weaker|wins?|"
                r"rank\w*|top|leads?|ahead|first of|last of)\b")
        named = 0
        for s in sentences(self.page):
            if re.search(LEARNER, s):
                named += 1
                self.assertIsNone(re.search(rank, s, re.I), s)
        self.assertGreaterEqual(named, 20)
        self.assertIsNotNone(re.search(rank, "rf ranks first of the five learners.", re.I))

    def test_names_no_sequence_or_determinant(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text), name)
            self.assertIsNone(re.search(r"\bk-?mers?\b", text, re.I), name)
            self.assertIsNone(re.search(r"\b(genes?|alleles?|mutations?|determinants?|mechanisms?)\b", text, re.I),
                              name)


if __name__ == "__main__":
    unittest.main()
