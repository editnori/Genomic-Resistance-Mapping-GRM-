"""The words of the page, the record and the code: what the experiment shows, and what they never hold."""

import re
import unittest

from tests.dense_network.common import PAGE, RECORD, ROOT, SCRIPT, SUITE

BANNED = (r"\b(validated|generali[sz]\w*|robust\w*|reliab\w*|real-world|real world|carriers?|provenance|receipts?|"
          r"manifests?|folders?|no tuning|prospective\w*|clinic\w*|patients?|mechanisms?|determinants?|genes?|alleles?|"
          r"mutations?|variants?|markers?|resistance arises|confers?|equivalent to|as good as|on par|weights?|"
          r"importance|k-?mers?|outperform\w*|superior|inferior|state of the art|deep learning)\b")
# The names of the work's parts and of the hosts it ran on, written so that this file holds none of them.
NAMES = (r"(?i)\b(reck[o]v|revarc[h]|fixarc[h]|arc[h]|b[o]x(es)?|c[o]ntroller|thr[e]ads?|w[o]rkers?|[o]wner\d?|"
         r"l[a]nes?|sub-?[a]gents?|boa[t]|bx_\w+|smok[e])\b")
VENDOR = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|ass[i]stant|l[l]m|m[o]del name)\b"
GENOME_ID = r"(?<![\d.])\d{2,7}\.\d{3,6}(?![\d.])"
ACCESSION = r"\b(GC[AF]_\d{9}(\.\d+)?|SAM[ND]\w?\d{6,}|SAMEA\d{6,}|[SE]RR\d{6,}|PRJ[NED][AB]\d+)\b"
SEQUENCE = r"(?i)[acgt]{20,}"
# What the experiment cannot carry: a claim about networks in general, or a gain ruled out.
WIDER = (r"(?i)\b(dense networks|neural networks|networks (do|does|did|can|cannot)|cannot help|does not help|"
         r"no network|rules? out (a|any|the) gain|shows? no gain|worse than GRM|adds? nothing|no use|of no value|useless)\b")
RECHECK = r"(?i)(re-?check|23 of 27|earlier claims)"


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
        cls.records = record_files()
        cls.all = dict(cls.texts, **cls.records, **code_files())

    def test_banned_words(self):
        for name, text in dict(self.texts, **self.records).items():
            self.assertIsNone(re.search(BANNED, prose(text), re.I), name)
        self.assertIsNotNone(re.search(BANNED, "A gene confers it.", re.I))
        self.assertIsNotNone(re.search(BANNED, "the weights of the network", re.I))

    def test_no_name_of_the_work_or_vendor(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(NAMES, text), name)
            self.assertIsNone(re.search(VENDOR, text), name)

    def test_no_machine(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(r"(?i)\bmachines?\b", text), name)

    def test_external_validation_and_equivalence_only_denied(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                for m in re.finditer(r"(?i)(external|clinical) validation|equivalence", s):
                    found += 1
                    self.assertRegex(s[: m.start()], r"(?i)\b(no|not|never)\b", s)
        self.assertGreaterEqual(found, 1)

    def test_marks(self):
        for name, text in self.all.items():
            self.assertNotIn("\u2014", text, name)
            self.assertNotIn("\u2013", text, name)

    def test_no_sequence_genome_or_accession(self):
        for name, text in self.all.items():
            self.assertIsNone(re.search(SEQUENCE, text), name)
            self.assertIsNone(re.search(ACCESSION, text), name)
        for name, text in dict(self.texts, **self.records).items():
            self.assertIsNone(re.search(GENOME_ID, text), name)
        self.assertIsNotNone(re.search(GENOME_ID, "a genome 1280.12345 here"))
        self.assertIsNone(re.search(GENOME_ID, "-0.059308 and 0.0001"))
        # Built at run time, so that this file holds no run of the four letters.
        self.assertIsNotNone(re.search(SEQUENCE, "event " + "ac" + "gt" * 10))

    def test_no_reading_wider_than_the_experiment(self):
        for name, text in self.texts.items():
            for s in sentences(text):
                self.assertIsNone(re.search(WIDER, s), s)
        self.assertIsNotNone(re.search(WIDER, "Dense networks cannot help here."))
        self.assertIsNotNone(re.search(WIDER, "The network adds nothing to GRM's learners."))
        self.assertIn("no transformer or sequence model was.", self.page)
        self.assertIn("and not about the family", self.page)

    def test_no_recheck_of_earlier_claims(self):
        for name, text in dict(self.texts, **self.records).items():
            self.assertIsNone(re.search(RECHECK, text), name)

    def test_association_wording(self):
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(causes?|drives?|explains?|leads? to) resistance\b", s), s)

    def test_planned_and_exploratory_readings_named(self):
        self.assertIn("so they are exploratory", self.page)
        self.assertIn("Beside the rule, against each of the five", self.page)
        self.assertIn("`below_best5_lineage_*`", self.texts["README"])
        self.assertIn("so they are exploratory", self.texts["README"])

    def test_recall_beside_log_loss(self):
        heads = [x for x in self.page.split("\n") if x.startswith("| Comparator |")]
        self.assertEqual(len(heads), 2)
        for head in heads:
            cols = [c.strip() for c in head.strip("|").split("|")]
            self.assertEqual(cols.index("Log loss") - cols.index("Recall"), 1)


if __name__ == "__main__":
    unittest.main()
