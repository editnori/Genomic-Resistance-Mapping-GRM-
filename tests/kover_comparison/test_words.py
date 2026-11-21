"""The words of the page, the record and the code: what the comparison is called, and what they never hold."""

import re
import unittest

from tests.kover_comparison.common import PAGE, RECORD, ROOT, SCRIPT, SUITE

BANNED = (r"\b(validated|generali[sz]\w*|robust\w*|reliab\w*|real-world|real world|carriers?|provenance|receipts?|"
          r"manifests?|folders?|no tuning|prospective\w*|clinic\w*|patients?|mechanisms?|determinants?|genes?|alleles?|"
          r"mutations?|variants?|resistance arises|confers?|equivalent to|as good as|on par)\b")
# The names of the work's parts and of the hosts it ran on, written so that this file holds none of them.
NAMES = (r"(?i)\b(reck[o]v|revk[o]ver|k[o]ver2|b[o]x(es)?|c[o]ntroller|thr[e]ads?|w[o]rkers?|[o]wner\d?|l[a]nes?|"
         r"sub-?[a]gents?|boa[t]|bx_\w+)\b")
VENDOR = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|ass[i]stant|l[l]m|m[o]del name)\b"
GENOME_ID = r"(?<![\d.])\d{2,7}\.\d{3,6}(?![\d.])"
ACCESSION = r"\b(GC[AF]_\d{9}(\.\d+)?|SAM[ND]\w?\d{6,}|SAMEA\d{6,}|[SE]RR\d{6,}|PRJ[NED][AB]\d+)\b"
SEQUENCE = r"[ACGT]{20,}"


def prose(text):
    """Text outside code blocks and inline code."""
    text = re.sub(r"```.*?```", " ", text, flags=re.S)
    return re.sub(r"`[^`]*`", "CODE", text)


def sentences(text):
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
    return {p.name: p.read_text(encoding="utf-8") for p in sorted(RECORD.iterdir())}


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = PAGE.read_text(encoding="utf-8")
        cls.texts = {"page": cls.page, "README": (RECORD / "README.md").read_text(encoding="utf-8")}
        cls.all = dict(cls.texts, **record_files(), **code_files())

    def test_banned_words(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BANNED, text, re.I), name)
        self.assertIsNotNone(re.search(BANNED, "A gene confers it.", re.I))

    def test_no_name_of_the_work_or_vendor(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(NAMES, text), name)
            self.assertIsNone(re.search(VENDOR, text), name)

    def test_machine_only_in_set_covering_machine(self):
        for name, text in self.all.items():
            for m in re.finditer(r"(?i)\bmachines?\b", text):
                self.assertEqual(text[m.start() - 13:m.end()].lower(), "set covering machine", name)

    def test_external_validation_and_equivalence_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"(?i)(external|clinical) validation|equivalence", s):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 2)

    def test_marks(self):
        for name, text in self.all.items():
            self.assertNotIn("\u2014", text, name)
            self.assertNotIn("\u2013", text, name)

    def test_no_sequence_genome_or_accession(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(SEQUENCE, text), name)
            self.assertIsNone(re.search(ACCESSION, text), name)
        for name, text in dict(self.texts, **record_files()).items():
            self.assertIsNone(re.search(GENOME_ID, text), name)
        self.assertIsNotNone(re.search(GENOME_ID, "a genome 1280.12345 here"))
        self.assertIsNone(re.search(GENOME_ID, "-0.104287 and 0.0001"))

    def test_no_rule_or_weight_content(self):
        for name, text in dict(self.texts, **record_files()).items():
            self.assertIsNone(re.search(r"(?i)(kmer_rule|rule_terms?\b|presence\(|absence\()", text), name)

    def test_marker_and_association_wording(self):
        self.assertIn("A rule term is a marker and an association.", self.page)
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(causes resistance|drives resistance|explains resistance)\b", s), s)

    def test_planned_and_exploratory_readings_named(self):
        self.assertIn("in an exploratory reading", self.page)
        self.assertIn("so they are exploratory", self.page)
        self.assertIn("| Error rate |", self.page)
        row = [x for x in self.page.split("\n") if x.startswith("| Error rate |")][0]
        self.assertTrue(row.endswith("| exploratory |"))

    def test_recall_beside_log_loss(self):
        rows = [x.split(" | ")[0][2:] for x in self.page.split("\n") if x.startswith("| ") and " to " in x
                and x.split(" | ")[0][2:] in ("Balanced accuracy", "Recall", "Log loss", "Specificity")]
        self.assertEqual(rows.index("Log loss") - rows.index("Recall"), 1)


if __name__ == "__main__":
    unittest.main()
