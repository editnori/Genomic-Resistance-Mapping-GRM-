"""The words of the page, the record, the script and this suite: what the comparison is called, the review's wording,
and words and shapes they never use."""

import re
import unittest

from tests.tbprofiler_comparison.common import CI_DOC, PAGE, RECORD, ROOT, SCRIPT, SUITE, load_script

record = load_script()
BANNED = (r"\b(validated|generalis|generaliz|robust|reliab|real-world|real world|carrier|carries|provenance|receipt|"
          r"manifest|folder|safe|help|explain|caus|due to|clinic|patient|scope|transfer|threshold|confer|proves?|proof|"
          r"guarantee|pass(es|ed|ing)?|reference method|october|rented|work program|test isolate|seconds|minutes?|"
          r"hours?|core-hours?|cpu|wall[- ]?clock|run ?time|isolate'?s name|sample name)\b")
# A run of bases, in either case. The test checks the pattern on runs it builds, so this file holds none.
SEQUENCE = r"(?i)[acgt]{20,}"
# Words for the people, places and tools of the work, written so that this file holds none of them.
ROLE = r"(?i)\b(la[n]es?|own[e]rs?|mach[i]nes?|b[o]x(es)?|contr[o]ller|threa[d]s?|work[e]rs?|assist[a]nt|mod[e]ls?|vend[o]rs?|agen[t]s?)\b"
# Names of the work's parts and of the hosts it ran on.
NAMES = r"(?i)\b(rect[b]|tbwh[o]|revtbwh[o]|own[e]r2|bo[a]t|bx_[a-z0-9]+|grm-[o]ps|sandb[o]x|h3[7]r[va])\b"
MAKER = r"(?i)\b(cl[a]ude|anthr[o]pic|co-author[e]d-by|generated[ ]with|open[a]i|chatg[p]t|g[p]t-?\d|gem[i]ni)\b"
# Labels of the review's findings and notes, which a public reader cannot look up.
LABEL = r"\b(B|S|N|SF|F)\d{1,2}\b"
LEARNER = r"\b(scm|cart|rf|gbt|lasso)\b"
FILE_VERBS = ("gives", "give", "lists", "list", "shows", "show", "holds", "hold", "is", "are")
# Words that would give a biological reason for a call.
BIOLOGY = (r"\b(efflux|pumps?|promoters?|operons?|plasmids?|enzymes?|proteins?|binding|compensat\w*|heteroresist\w*|"
           r"intrinsic\w*|express\w*|frameshifts?|epista\w*|target sites?|loss of function|silent|synonymous|"
           r"lineage-specific|fitness|mechanis\w*)\b")
TB_GENE = (r"\b(kat|rpo|inh|emb|pnc|gyr|eth|rpl|ddn|fbi|fgd|atp|mmp|tly|gid|thy|fol|rib|ald|eis|ahp|fab|ubi|pep|msh|"
           r"whi|Rv)[A-Z0-9]|\brr[sl]\b")
SCREEN = ("the 0.02 margin is an adapted research screen the plan proposed and no one has adopted, and it is never an "
          "applicable FDA validation criterion")
CATALOG = "GRM's own reader of the AMRFinderPlus catalog"
NO_MECHANISM = "a call is a marker read against a phenotype and establishes no mechanism"
BECAUSE = "because its database grades no entry for the other 2 drugs and it calls neither resistant"
SETTINGS = "production selects its settings inside each fit"


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


def tables(text):
    """The header line of each table outside code."""
    rows = prose(text).split("\n")
    return [a for a, b in zip(rows, rows[1:]) if a.startswith("| ") and b.startswith("|---")]


def code_files():
    return {p.relative_to(ROOT).as_posix(): p.read_text(encoding="utf-8")
            for p in [SCRIPT] + sorted((ROOT / SUITE).glob("*.py"))}


def record_files():
    return {"reference/tbprofiler-comparison/" + p.name: p.read_text(encoding="utf-8") for p in sorted(RECORD.iterdir())}


def ci_sentence():
    text = CI_DOC.read_text(encoding="utf-8")
    start = text.index(f"The `{SUITE}` suite")
    end = text.find(" The `tests/", start + 1)
    line_end = text.index("\n", start)
    return text[start:line_end] if end < 0 or end > line_end else text[start:end]


def without(text, *allowed):
    for words in allowed:
        text = text.replace(words, " ")
    return text


class Words(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.page = PAGE.read_text(encoding="utf-8")
        cls.texts = {"page": cls.page, "README": (RECORD / "README.md").read_text(encoding="utf-8"),
                     "docs/ci.md": ci_sentence()}

    def everything(self):
        return list(self.texts.items()) + list(code_files().items()) + list(record_files().items())

    def test_banned_words(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BANNED, text, re.I), name)

    def test_no_name_of_the_work_its_hosts_or_a_vendor(self):
        for name, text in self.everything():
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
        found = [s for s in sentences(self.page) if re.search(r"(?i)\bprospective", s)]
        self.assertEqual(found, ["It is not a prospective test and not external validation."])

    def test_the_screen_keeps_the_reviews_words(self):
        self.assertEqual(self.page.count(SCREEN), 1)
        for name, text in self.texts.items():
            self.assertNotRegex(without(text, SCREEN), r"\bFDA\b", name)
            for s in sentences(text):
                if re.search(r"(?i)(research screen|0\.02 margin)", s):
                    self.assertIn(SCREEN, s.replace("CODE", ""))

    def test_the_catalog_row_is_named_one_way(self):
        for name, text in self.texts.items():
            self.assertNotIn("AMRFinderPlus", without(text, CATALOG), name)
        self.assertEqual(self.page.count(CATALOG), 1)

    def test_settings_are_selected_inside_each_fit(self):
        self.assertEqual(self.page.count(SETTINGS), 1)
        self.assertIsNone(re.search(r"(?i)\b(untuned|without tuning|not tuned|no tuning)\b", self.page))

    def test_no_biological_reason(self):
        for name, text in self.texts.items():
            self.assertIsNone(re.search(BIOLOGY, without(text, NO_MECHANISM), re.I), name)
        self.assertEqual(self.page.count(NO_MECHANISM), 1)
        # The one reason the page gives for a zero sensitivity is the database's grade table.
        for s in sentences(self.page):
            if re.search(r"(?i)\b(because|since|owing|reason|why)\b", s):
                self.assertIn(BECAUSE, s)
        self.assertEqual(self.page.count("because"), 1)

    def test_marks(self):
        for name, text in self.everything():
            self.assertNotIn(chr(0x2014), text, name)
            self.assertIsNone(re.search(r"\brather[ ]than\b", text, re.I), name)

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
        self.assertGreaterEqual(named, 1)

    def test_no_learners_pooled(self):
        # Each effect is per learner: no figure for the five learners taken together.
        for s in sentences(self.page):
            self.assertIsNone(re.search(r"(?i)\b(the learners'|all learners|the five learners') (median|mean|figure|"
                                        r"sensitivity|specificity|balanced accuracy)", s), s)
            self.assertIsNone(re.search(r"(?i)\b(pooled|averaged over the learners|on average)\b", s), s)

    def test_sensitivity_beside_specificity(self):
        found = 0
        for name, text in self.texts.items():
            for s in sentences(text):
                if re.search(r"(?i)\bsensitivit", s):
                    found += 1
                    self.assertRegex(s, r"(?i)\bspecificit", s)
            for head in tables(text):
                if re.search(r"(?i)sensitivit", head):
                    self.assertRegex(head, r"(?i)specificit", head)
        self.assertGreaterEqual(found, 6)

    def test_the_base_run_pattern_checks_itself(self):
        for run in ("acgt" * 6 + "a", "ACGT" * 6 + "A", "aCgT" * 5, "t" * 20):
            self.assertRegex(run, SEQUENCE, run)
        for run in ("acgt" * 4 + "acg", "acgtn" * 5, "ACGU" * 6):
            self.assertNotRegex(run, SEQUENCE, run)
        self.assertIn(SEQUENCE, record.BANNED_SHAPES)

    def test_names_no_sequence_or_determinant(self):
        for name, text in self.everything():
            self.assertIsNone(re.search(SEQUENCE, text), name)
        for name, text in list(self.texts.items()) + list(record_files().items()):
            self.assertIsNone(re.search(r"\bk-?mers?\b|\b\d+-mers?\b", text, re.I), name)
            self.assertIsNone(re.search(r"\b(genes?|alleles?|mutations?|variants?|determinants?|snps?|indels?)\b",
                                        text, re.I), name)
            self.assertIsNone(re.search(r"\b(bla|aac|aph|ant|tet|sul|dfr|gyr|par|mcr|qnr|erm|van|mec|omp|acr)[A-Z(-]",
                                        text), name)
            self.assertIsNone(re.search(TB_GENE, text), name)
            self.assertIsNone(re.search(r"\b[pcn]\.[A-Z*-]?-?\d|\d[ACGT]>[ACGT]", text), name)
            self.assertIsNone(re.search(r"\b[A-Z][0-9]{1,4}[A-Z*]\b", prose(text) if name in self.texts else text), name)

    def test_names_no_study_genome_or_sample(self):
        for name, text in list(self.texts.items()) + list(record_files().items()):
            self.assertIsNone(re.search(r"\bPRJ[A-Z]{1,2}\d+|\bSAM[NED][A-Z]?\d+|\bGC[AF]_\d+|\b[SE]R[RSX]\d{5,}", text),
                              name)
            self.assertIsNone(re.search(r"\b\d{3,7}\.\d{3,6}\b", text), name)
            for shape in record.BANNED_SHAPES:
                self.assertIsNone(re.search(shape, text), (name, shape))

    def test_exploratory_is_named(self):
        # The kanamycin counts and the no-call bounds stand in the section Exploratory alone, after its label.
        self.assertGreaterEqual(len(re.findall(r"\*\*Exploratory\.\*\*", self.page)), 2)
        exploratory = self.page.split("\n## Exploratory\n", 1)[1].split("\n## ", 1)[0]
        for words in ("set aside", "no call counted wrong", "kanamycin report"):
            self.assertIn(words, exploratory)
            self.assertNotIn(words, self.page.replace(exploratory, ""))


if __name__ == "__main__":
    unittest.main()
