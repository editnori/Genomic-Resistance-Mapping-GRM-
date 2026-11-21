"""The words of the page, the record and the code: the three callers kept apart, the review's wording, and what they
never hold."""

import re
import unittest

from tests.amrfinderplus_comparison.common import PAGE, RECORD, ROOT, SCRIPT, SUITE

BANNED = (r"\b(validated|generali[sz]\w*|robust\w*|reliab\w*|real-world|real world|carriers?|provenance|receipts?|"
          r"manifests?|folders?|no tuning|clinic\w*|patients?|mechanisms?|determinants?|alleles?|mutations?|"
          r"mutational|variants?|resistance arises|confers?|equivalent to|as good as|on par|in general|"
          r"outperform\w*|superior|inferior|matches or exceeds|matched or exceeded|pass(es|ed)?)\b")
# The names of the work's parts and of the hosts it ran on, written so that this file holds none of them.
NAMES = (r"(?i)\b(amrf[u]ll|reca[m]r|reva[m]r|b[o]x(es)?|c[o]ntroller|thr[e]ads?|w[o]rkers?|[o]wner\d?|l[a]nes?|"
         r"sub-?[a]gents?|boa[t]|bx_\w+|scr[a]tch)\b")
VENDOR = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|ass[i]stant|l[l]m|m[o]del name)\b"
GENOME_ID = r"(?<![\d.])\d{2,7}\.\d{3,6}(?![\d.])"
ACCESSION = r"\b(GC[AF]_\d{9}(\.\d+)?|SAM[ND]\w?\d{6,}|SAMEA\d{6,}|[SE]RR\d{6,}|PRJ[NED][AB]\d+)\b"
SEQUENCE = r"(?i)[acgt]{20,}"
# The program's name stands only before "program" or in the page's and the record's titles.
PROGRAM_NAME = r"AMRFinderPlus(?! (program|comparison)\b)"
RULE_AS_PROGRAM = (r"(?i)(AMRFinderPlus('s)?\W+(known-gene )?rule\b|known-gene rule\W{1,3}(\(|is |of )?"
                   r"(the )?AMRFinderPlus|AMRFinderPlus's known-gene)")
# A catalog row's class and subclass labels and the program's method codes; the page names a method only in words.
RULE_CONTENT = r"(\bC:[A-Z]|\bS:[A-Z]|(?i:subclass)|\b(POINT_DISRUPT|INTERNAL_STOP|PARTIAL\w*|EXACT[NPX]|BLAST[NPX])\b)"
CAUSE = r"(?i)\b(due to|caus\w*|because of|explain\w*|reason\w* (in|for) the biology|drives?|driven)\b"


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
        self.assertIsNotNone(re.search(BANNED, "A marker confers it.", re.I))
        self.assertIsNotNone(re.search(BANNED, "The rule passes.", re.I))
        self.assertIsNotNone(re.search(BANNED, "It matches or exceeds the rule.", re.I))

    def test_no_196_lead(self):
        for name, text in self.texts.items():
            self.assertNotIn("196", text, name)
            self.assertIsNone(re.search(r"(?i)\bmatch\w* or exceed", text), name)

    def test_no_name_of_the_work_or_vendor(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(NAMES, text), name)
            self.assertIsNone(re.search(VENDOR, text), name)
        self.assertIsNotNone(re.search(NAMES, "the run b" + "ox"))

    def test_known_gene_rule_never_takes_the_programs_name(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(PROGRAM_NAME, prose(text)), name)
            self.assertIsNone(re.search(RULE_AS_PROGRAM, prose(text)), name)
        self.assertGreaterEqual(self.page.count("GRM's known-gene rule"), 4)
        self.assertIsNotNone(re.search(RULE_AS_PROGRAM, "GRM's known-gene rule (AMRFinderPlus) read"))
        self.assertIsNotNone(re.search(RULE_AS_PROGRAM, "the AMRFinderPlus rule"))
        self.assertIsNotNone(re.search(PROGRAM_NAME, "the AMRFinderPlus known-gene rule"))

    def test_external_validation_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"(?i)(external|clinical) validation|equivalence", s):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 2)

    def test_rules_g_and_m_in_the_reviews_words(self):
        self.assertIsNone(re.search(r"\b[Rr]ule G\b", self.page))
        body = [s for s in sentences(self.page) if s.startswith("Rule M")]
        self.assertTrue(body)
        self.assertIn("Rule M, a proposed screen from a review note that no one adopted,", self.page)
        self.assertIn("It is an adapted research screen, never an applicable FDA validation criterion.", self.page)
        for s in sentences(self.page):
            if "FDA" in s:
                self.assertRegex(s, r"never an applicable FDA validation criterion", s)

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
        # Built at run time, so that this file holds no run of the four letters, in either case.
        self.assertIsNotNone(re.search(SEQUENCE, "kind " + "ac" + "gt" * 10))
        self.assertIsNotNone(re.search(SEQUENCE, "kind " + "AC" + "GT" * 10))
        self.assertIsNone(re.search(SEQUENCE, "kind " + "ac" + "gt" * 8))

    def test_no_rule_content_and_no_disagreement_table(self):
        for name, text in dict(self.texts, **record_files()).items():
            self.assertIsNone(re.search(RULE_CONTENT, text), name)
        self.assertFalse([p.name for p in RECORD.iterdir() if "disagree" in p.name.lower()])
        head = (RECORD / "datasets.tsv").read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
        self.assertNotIn("rule", head)
        self.assertIsNotNone(re.search(RULE_CONTENT, "a row marked PARTIAL_CONTIG_END"))
        self.assertIsNone(re.search(RULE_CONTENT, "the program marks a row as partial"))

    def test_marker_and_association_wording_and_no_cause(self):
        self.assertIn("A catalog row that meets a rule is a marker and an association.", self.page)
        for s in sentences(self.page):
            self.assertIsNone(re.search(CAUSE, s), s)
            # The one "because" a quinolone sentence holds is the mean's reason, a count, as the review asks.
            if "quinolone" in s:
                self.assertNotRegex(s, r"(?i)\bbecause\b(?! the program is ahead by 0\.10 or more in)", s)


if __name__ == "__main__":
    unittest.main()
