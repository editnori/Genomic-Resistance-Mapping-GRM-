"""The words of the page, the record, the script and this suite: what the test is called, and words they never use."""

import re
import unittest

from tests.time_split.common import CI_DOC, PAGE, RECORD, ROOT, SCRIPT, SUITE

BANNED = (r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|carrier|carries|provenance|receipt|"
          r"manifest|folder|safe|help|explain|because|caus|due to|no tuning|clinic|patient|scope|transfer|threshold|"
          r"mechanis|confer|proves?|proof|guarantee)")
# Words for the people, places and tools of the work, written so that this file holds none of them.
ROLE = r"(?i)\b(la[n]es?|own[e]rs?|mach[i]nes?|b[o]x(es)?|contr[o]ller|threa[d]s?|work[e]rs?|assist[a]nt|mod[e]ls?|vend[o]rs?|agen[t]s?)\b"
# Names of the work's parts and of the hosts it ran on.
NAMES = r"(?i)\b(times[p]lit|rect[i]me|fix1[0-9]{2}|rev1[0-9]{2}|own[e]r2|bo[a]t|bx_[a-z0-9]+|grm-[o]ps|sandb[o]x)\b"
MAKER = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|g[p]t-?\d|gem[i]ni)\b"
# Labels of the review's findings and notes, which a public reader cannot look up.
LABEL = r"\b(B|S|N|SF|F)\d{1,2}\b"
LEARNER = r"\b(scm|cart|rf|gbt|lasso)\b"
FILE_VERBS = ("gives", "give", "lists", "list", "shows", "show", "holds", "hold", "is", "are")


def prose(text):
    """Text outside code blocks and inline code, with emphasis marks dropped."""
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    text = re.sub(r"`[^`]*`", "CODE", text)
    return text.replace("*", "")


def sentences(text):
    """Sentences of the prose lines and list items; a stop after one capital letter does not end one."""
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


def record_files():
    return {"reference/time-split/" + p.name: p.read_text(encoding="utf-8") for p in sorted(RECORD.iterdir())}


def ci_sentence():
    text = CI_DOC.read_text(encoding="utf-8")
    start = text.index(f"The `{SUITE}` suite")
    end = text.find(" The `tests/", start + 1)
    line_end = text.index("\n", start)
    return text[start:line_end] if end < 0 or end > line_end else text[start:end]


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = PAGE.read_text(encoding="utf-8")
        cls.texts = {"page": cls.page, "README": (RECORD / "README.md").read_text(encoding="utf-8"),
                     "docs/ci.md": ci_sentence()}

    def test_banned_words(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BANNED, text, re.I), name)

    def test_no_name_of_the_work_its_hosts_or_a_vendor(self):
        everything = list(self.texts.items()) + list(code_files().items()) + list(record_files().items())
        for name, text in everything:
            self.assertIsNone(re.search(NAMES, text), name)
            self.assertIsNone(re.search(MAKER, text), name)
            self.assertIsNone(re.search(ROLE, text), name)
        for name, text in self.texts.items():
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
                self.assertIsNone(re.search(r"(?i)clinical validation", s), s)
        self.assertEqual(found, 1)

    def test_prospective_only_in_the_limit_sentence(self):
        found = [s for s in sentences(self.page) if re.search(r"(?i)prospective", s)]
        # The sentence split keeps "BV-BRC." and the next sentence together, so the two limit sentences come as one.
        self.assertEqual(found, ["This is a split by first-public day inside one database, BV-BRC. It is not a prospective "
                                 "test and not external validation."])

    def test_the_outside_study_wording(self):
        self.assertEqual(self.page.count("a retrospective test on newly deposited studies from the same database"), 1)
        for s in sentences(self.page):
            if re.search(r"(?i)\bdeposit", s):
                self.assertIn("a retrospective test on newly deposited studies from the same database", s)

    def test_settings_are_selected_inside_each_fit(self):
        self.assertIn("Each fit selects candidate features and tunes its settings inside its own training genomes.",
                      self.page)
        self.assertIsNone(re.search(r"(?i)\b(untuned|without tuning|not tuned)\b", self.page))

    def test_marks(self):
        for name, text in list(self.texts.items()) + list(code_files().items()) + list(record_files().items()):
            self.assertNotIn(chr(0x2014), text, name)
            self.assertIsNone(re.search(r"\brather than\b", text, re.I), name)

    def test_files_give_list_show_or_hold(self):
        verbs = 0
        for name, text in self.texts.items():
            for m in re.finditer(r"(`[^`]+\.(tsv|md|sha256|py)`|`SHA256SUMS`|\]\([\w.-]+\.md\)|\bREADME)"
                                 r" ([a-z]+)\b", text):
                if m.group(3) in ("and", "in", "of", "from", "beside", "that", "with", "to", "for", "on", "or", "again"):
                    continue
                verbs += 1
                allowed = FILE_VERBS + (("writes",) if m.group(2) == "py" else ())
                self.assertIn(m.group(3), allowed, m.group(0))
        self.assertGreaterEqual(verbs, 3)

    def test_no_file_is_the_subject_of_a_verb_of_speech(self):
        speech = r"(says?|speaks?|thinks?|tells?|claims?|argues?|believes?|states?|explains?|mentions?)"
        for name, text in self.texts.items():
            self.assertIsNone(re.search(rf"\b(page|file|table|test|record|review|suite|README|plan|script|list)s? "
                                        rf"{speech}\b", text, re.I), name)

    def test_no_learner_ranked_above_another(self):
        rank = (r"\b(best|worst|outperform\w*|superior|inferior|beats?|stronger|weaker|wins?|rank\w*|top|leads?|ahead|"
                r"first of|last of)\b")
        named = 0
        for s in sentences(self.page):
            if re.search(LEARNER, s):
                named += 1
                self.assertIsNone(re.search(rank, s, re.I), s)
        self.assertGreaterEqual(named, 8)

    def test_recall_beside_log_loss(self):
        found = 0
        for s in sentences(self.page):
            if re.search(r"(?i)\blog[ _]loss\b", s):
                found += 1
                self.assertRegex(s, r"(?i)\brecall\b", s)
        self.assertGreaterEqual(found, 12)

    def test_names_no_sequence_or_determinant(self):
        for name, text in list(self.texts.items()) + list(record_files().items()):
            self.assertIsNone(re.search(r"\b[ACGT]{8,}\b", text), name)
            self.assertIsNone(re.search(r"\bk-?mers?\b|\b\d+-mers?\b", text, re.I), name)
            self.assertIsNone(re.search(r"\b(genes?|alleles?|mutations?|variants?|markers?|determinants?|snps?)\b",
                                        text, re.I), name)
            self.assertIsNone(re.search(r"\b(bla|aac|aph|ant|tet|sul|dfr|gyr|par|mcr|qnr|erm|van|mec|omp|acr)[A-Z(-]",
                                        text), name)
            self.assertIsNone(re.search(r"\b[A-Z][0-9]{1,4}[A-Z*]\b", prose(text) if name in self.texts else text), name)

    def test_names_no_study_genome_or_sample(self):
        for name, text in list(self.texts.items()) + list(record_files().items()):
            self.assertIsNone(re.search(r"\bPRJ[A-Z]{1,2}\d+|\bSAM[NED][A-Z]?\d+|\bGC[AF]_\d+|\b[SE]R[RSX]\d{5,}", text),
                              name)
            self.assertIsNone(re.search(r"\b\d{3,7}\.\d{3,6}\b", text), name)

    def test_exploratory_is_named(self):
        self.assertGreaterEqual(len(re.findall(r"(?i)\bexploratory\b", self.page)), 4)


if __name__ == "__main__":
    unittest.main()
