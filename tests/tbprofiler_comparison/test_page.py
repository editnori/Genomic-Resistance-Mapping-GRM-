"""docs/tbprofiler-comparison.md against the committed record: the limit sentence, the section Reading, the
paragraphs and tables built from the record, the review's sentences the page must keep, every other paragraph and
list item word for word, and every other number of the page."""

import re
import unittest

from tests.tbprofiler_comparison import lines
from tests.tbprofiler_comparison.common import LEARNERS, PAGE, RECORD, WAY, pairs, read_tsv
from tests.tbprofiler_comparison.held import HELD
from tests.tbprofiler_comparison.lines import f3, n

TITLE = "# TBProfiler comparison"
LIMIT = ("This is a retrospective comparison on assemblies of public BV-BRC genomes. "
         "It is not a prospective test and not external validation.")
SECTIONS = ["Reading", "Plan and order", "Method", "Reports", "Medians", "Learners beside the tool", "Ways of reading",
            "Disagreements by grade", "Exploratory", "Review", "Limits and the record"]
TABLE_SECTION = {"order": "Plan and order", "grades": "Method", "primary": "Reports", "proxy": "Reports",
                 "medians": "Medians", "comparison": "Learners beside the tool", "gap": "Learners beside the tool",
                 "ways": "Ways of reading", "disagreement": "Disagreements by grade", "bounds": "Exploratory"}
BUILT_SECTION = ["Reports", "Learners beside the tool", "Exploratory"]
# The review's sentence on the catalogue's overlap, beside every figure that sets the tool beside a learner.
OVERLAP = ("The catalogue behind TBProfiler was built from public isolates with recorded phenotypes, and some of these "
           "2,000 genomes may be among them, which would favour the tool; the overlap was not checked, and each learner "
           "was scored only on genomes held out of its fit.")
OVERLAP_SECTIONS = ["Medians", "Learners beside the tool", "Ways of reading", "Exploratory"]
# The review's other sentences, word for word.
KEPT = {
    "diagnostic": ("This is a diagnostic comparison of a published catalogue-based caller: a call is a marker read "
                   "against a phenotype and establishes no mechanism."),
    "input limit": ("TBProfiler is built for sequencing reads and was run here on assemblies, which cannot show a mixed "
                    "population, and it gave no call in 819 of 21,737 genome-drug rows where an assembly did not cover a "
                    "target region in full; these figures say nothing about the tool on reads."),
    "lowest layer": ("The tool's raw output and the assemblies it read no longer exist; the per-genome call file, kept "
                     "outside this repository, is the lowest layer these counts can be rebuilt from."),
    "interval": ("Each of the 2,000 resamples draws 881 group numbers with Python's "
                 "`random.Random(20261114).randrange(881)`, one call per draw, with the groups numbered in the text "
                 "order of their smallest genome id. The interval's ends are the 51st and the 1,950th of the 2,000 "
                 "sorted values."),
    "screen": ("On the lineage split and the same genomes, TBProfiler's balanced accuracy was ahead of or within 0.02 of "
               "the best of the five learners on 12 of 16 endpoints and 12 of the 14 graded ones; the 0.02 margin is an "
               "adapted research screen the plan proposed and no one has adopted, and it is never an applicable FDA "
               "validation criterion."),
    "two drugs": ("The shipped database grades no entry for cycloserine or para-aminosalicylic acid (`grades.tsv`), so "
                  "the tool cannot call either resistant; this was written into the plan's addendum, hashed at "
                  "16:51:41Z, before the evaluator first read a recorded phenotype after 17:28:42Z."),
}


def outside_code(text):
    out, fenced = [], False
    for line in text.split("\n"):
        if line.startswith("```"):
            fenced = not fenced
            continue
        if not fenced:
            out.append(line)
    return out


def blocks(text):
    """The page's paragraphs and list items outside code, tables and headers, one string each, in page order."""
    out = []
    for para in "\n".join(outside_code(text)).split("\n\n"):
        para = para.strip("\n")
        if not para or para.startswith(("#", "|")):
            continue
        if para.startswith("- "):
            out += ["- " + x for x in para[2:].split("\n- ")]
        else:
            out.append(para)
    return out


def section(text, title):
    return text.split(f"\n## {title}\n", 1)[1].split("\n## ", 1)[0]


def page_tables(text):
    tables, cur = [], []
    for line in outside_code(text) + [""]:
        if line.startswith("|"):
            cur.append(line)
        elif cur:
            tables.append(cur)
            cur = []
    return tables


def sentences(paragraph):
    return re.split(r"(?<=\.)\s+(?=[A-Z])", paragraph)


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.rec = lines.Record()

    def test_first_paragraph(self):
        paragraphs = [p for p in self.text.split("\n\n") if p.strip()]
        self.assertEqual(paragraphs[0], TITLE)
        first = paragraphs[1]
        self.assertTrue(first.endswith(LIMIT), "the limit sentence does not close the first paragraph")
        self.assertEqual(self.text.count(LIMIT), 1)
        self.assertTrue(first.startswith("This page counts the calls of TBProfiler "))

    def test_sections_in_order(self):
        body = outside_code(self.text)
        self.assertEqual([x[3:] for x in body if x.startswith("## ")], SECTIONS)
        self.assertEqual([x for x in body if x.startswith("# ")], [TITLE])
        self.assertFalse([x for x in body if x.startswith("###")])

    # ------------------------------------------------------------- Reading

    def test_reading_is_rebuilt_from_the_record(self):
        body = section(self.text, "Reading").strip()
        self.assertEqual(body.split("\n\n"), lines.reading(self.rec))

    def test_reading_has_four_sentences(self):
        body = section(self.text, "Reading").strip()
        found = [x for p in body.split("\n\n") for x in sentences(p)]
        self.assertEqual(len(found), 4)
        for x in found:
            self.assertTrue(x.endswith("."), x)

    def test_reading_figures_are_the_records(self):
        # Read straight from the record's files, apart from lines.py, so a figure changed with its builder fails.
        body = section(self.text, "Reading").strip().split("\n\n")
        s, t = pairs(RECORD / "summary.tsv"), pairs(RECORD / "totals.tsv")
        med = {(r["way"], r["set"], r["caller"], r["measure"]): r for r in read_tsv(RECORD / "medians.tsv")}
        three = lambda key, m: tuple(f3(med[(WAY, key, "tool", m)][x]) for x in ("median", "min", "max"))
        self.assertIn(f"assemblies of the {n(t['distinct_genomes'])} genomes in its {t['reports_of_the_species']} "
                      f"revision 14", body[0])
        self.assertIn(f"on the {s['by_name_endpoints']} reports whose drug it names", body[0])
        self.assertIn(f"made a call on {n(s['by_name_rows_with_call'])} of {n(s['by_name_rows'])} genome-drug rows and "
                      f"none on the other {n(s['by_name_rows_without_call'])}", body[1])
        self.assertIn("median sensitivity %s (%s to %s)" % three("by_name", "sensitivity"), body[1])
        self.assertIn("beside median specificity %s (%s to %s)" % three("by_name", "specificity"), body[1])
        self.assertIn("they were %s (%s to %s)" % three("graded", "sensitivity") +
                      " and %s (%s to %s)" % three("graded", "specificity"), body[1])
        self.assertIn(f"with a call on {n(s['graded_rows_with_call'])} of {n(s['graded_rows'])} rows", body[1])
        self.assertIn(f"on the {s['graded_endpoints']} within these {s['by_name_endpoints']}", body[1])
        self.assertIn(f"whose drug its database grades (none for the other {s['ungraded_in_by_name']})", body[1])
        cmp_ = {(r["set"], r["caller"]): r for r in read_tsv(RECORD / "comparisons.tsv") if r["way"] == WAY}
        higher = [cmp_[("graded", f"{l}:lineage")]["tool_higher"] for l in LEARNERS]
        self.assertIn(f"higher than scm's on {higher[0]} of the {s['graded_endpoints']} M. tuberculosis reports whose "
                      f"drug its database grades,", body[2])
        self.assertIn(", ".join(f"{l}'s on {h}" for l, h in zip(LEARNERS[1:-1], higher[1:-1])) +
                      f" and {LEARNERS[-1]}'s on {higher[-1]}, and lower on the rest;", body[2])
        self.assertIn(f"left {n(s['graded_rows_without_call'])} of {n(s['graded_rows'])} rows of the "
                      f"{s['graded_endpoints']} graded reports without a call", body[3])
        # Read alone, the sentence that sets the tool beside the learners keeps the catalogue's possible overlap.
        self.assertTrue(body[2].endswith("; its catalogue may include some of these genomes, which would favour it."),
                        body[2])
        self.assertIn("; its catalogue may include some of these genomes, which would favour it, and each learner saw "
                      "only genomes held out of its fit.", body[3])

    def test_reading_states_how_the_counts_relate(self):
        body = section(self.text, "Reading")
        s = self.rec.summary
        self.assertIn("within these", body)
        self.assertEqual(int(s["by_name_rows_with_call"]) + int(s["by_name_rows_without_call"]), int(s["by_name_rows"]))
        self.assertIn(f"of {n(s['by_name_rows'])} genome-drug rows and none on the other", body)

    # ------------------------------------------------------------- built and held lines

    def test_built_paragraphs_stand_in_their_sections(self):
        for para, title in zip(lines.built(self.rec), BUILT_SECTION):
            self.assertEqual(blocks(section(self.text, title)).count(para), 1, para[:60])
            self.assertEqual(blocks(self.text).count(para), 1, para[:60])

    def test_every_paragraph_and_list_item_is_built_or_held(self):
        known = set(lines.reading(self.rec)) | set(lines.built(self.rec)) | set(HELD) | {OVERLAP}
        self.assertEqual([b for b in blocks(self.text) if b not in known], [])

    def test_held_lines_stand_once_and_in_order(self):
        self.assertEqual(len(HELD), len(set(HELD)))
        found = blocks(self.text)
        for line in HELD:
            self.assertEqual(found.count(line), 1, line)
        self.assertEqual([b for b in found if b in set(HELD)], HELD)
        built = set(lines.reading(self.rec)) | set(lines.built(self.rec)) | {OVERLAP}
        self.assertFalse(built & set(HELD), "a held line is also built")

    def test_overlap_sentence_beside_every_comparison_with_a_learner(self):
        self.assertEqual(blocks(self.text).count(OVERLAP), len(OVERLAP_SECTIONS))
        for title in OVERLAP_SECTIONS:
            self.assertEqual(blocks(section(self.text, title)).count(OVERLAP), 1, title)
        # Every table and paragraph that sets the tool beside a learner stands in a section with the sentence.
        for name in ("medians", "comparison", "gap", "ways"):
            self.assertIn(TABLE_SECTION[name], OVERLAP_SECTIONS, name)
        self.assertIn("best lineage learner minus the tool", lines.set_aside(self.rec))
        self.assertIn(BUILT_SECTION[2], OVERLAP_SECTIONS)
        reading = section(self.text, "Reading").strip().split("\n\n")
        for k in (2, 3):
            self.assertIn("its catalogue may include some of these genomes, which would favour it", reading[k])

    def test_the_reviews_sentences_are_kept(self):
        for name, words in KEPT.items():
            self.assertEqual(self.text.count(words), 1, name)
        first_two = [p for p in self.text.split("\n\n") if p.strip()][1:3]
        self.assertIn(KEPT["input limit"], first_two[1])
        self.assertIn(KEPT["diagnostic"], first_two[1])
        self.assertIn(KEPT["interval"], section(self.text, "Method"))
        self.assertIn(KEPT["screen"], section(self.text, "Learners beside the tool"))
        self.assertIn(KEPT["two drugs"], section(self.text, "Plan and order"))
        self.assertIn(KEPT["lowest layer"], section(self.text, "Limits and the record"))

    def test_review_figures_are_the_records(self):
        s, st = self.rec.summary, self.rec.settings
        screen = {r["set"]: r for r in read_tsv(RECORD / "screen.tsv") if r["way"] == WAY}
        self.assertEqual((screen["by_name"]["within_best"], screen["by_name"]["endpoints"]), ("12", "16"))
        self.assertEqual((screen["graded"]["within_best"], screen["graded"]["endpoints"]), ("12", "14"))
        self.assertEqual((s["by_name_rows_without_call"], s["by_name_rows"]), ("819", "21737"))
        self.assertEqual((st["resamples"], self.rec.totals["lineage_groups"], st["seed"], st["low_rank"], st["high_rank"]),
                         ("2000", "881", "20261114", "51", "1950"))
        order = {r["event"]: r["time"] for r in self.rec.order}
        self.assertEqual((order["addendum_hashed"], order["last_step_before_scoring"]),
                         ("2026-11-14T16:51:41Z", "2026-11-14T17:28:42Z"))
        ungraded = sorted(r["tool_drug"] for r in self.rec.grades if r["counted_grade_entries"] == "0")
        self.assertEqual(ungraded, ["cycloserine", "para-aminosalicylic_acid"])
        # Every interval had a value in all of its resamples, as the Method says.
        self.assertEqual(st["min_resamples_with_value"], "40")
        self.assertEqual({r["resamples_with_value"] for r in read_tsv(RECORD / "intervals.tsv")}, {st["resamples"]})

    def test_every_table_is_built_from_the_record(self):
        built = lines.tables(self.rec)
        found = page_tables(self.text)
        self.assertEqual(len(found), len(built))
        for name, table in built.items():
            self.assertEqual(sum(1 for t in found if t == table), 1, name)
            self.assertIn("\n" + "\n".join(table) + "\n", section(self.text, TABLE_SECTION[name]), name)

    def test_proxies_stand_apart(self):
        proxies = [r["antibiotic"] for r in self.rec.kind("proxy")]
        self.assertEqual(len(proxies), 4)
        rows = {name: [x.split(" | ")[0][2:] for x in table[2:]] for name, table in lines.tables(self.rec).items()}
        self.assertEqual(rows["proxy"], proxies)
        for name in ("primary", "gap", "disagreement"):
            self.assertFalse(set(rows[name]) & set(proxies), name)
        for r in read_tsv(RECORD / "medians.tsv"):
            self.assertEqual(r["endpoints"], self.rec.summary[f"{r['set']}_endpoints"], (r["way"], r["set"]))

    def test_every_number_is_in_the_record(self):
        found, _ = lines.values(RECORD)
        missing = [(x, line[:50]) for line in outside_code(self.text) for x in lines.numbers(line) if x not in found]
        self.assertEqual(missing, [])

    def test_hashes_are_the_records(self):
        st = self.rec.settings
        order = section(self.text, "Plan and order")
        self.assertIn(f"The plan (sha256 `{st['plan_sha256'][:8]}`) was hashed at", order)
        self.assertIn(f"its addendum (sha256 `{st['addendum_sha256'][:8]}`)", order)
        _, hashes = lines.values(RECORD)
        for h in re.findall(r"`([0-9a-f]{6,64})`", self.text):
            self.assertTrue(len(h) >= 8 and (any(x.startswith(h) for x in hashes) or h == st["database_commit"]), h)

    def test_order_follows_the_record(self):
        order = [r["event"] for r in self.rec.order]
        self.assertEqual(order, list(lines.ORDER_LABELS))
        times = [r["time"] for r in self.rec.order]
        self.assertEqual(times, sorted(times))

    def test_exploratory_figures_are_labelled(self):
        body = section(self.text, "Exploratory")
        self.assertTrue(blocks(body)[0].startswith("**Exploratory.**"))
        self.assertTrue(lines.set_aside(self.rec).startswith("**Exploratory.**"))
        self.assertIn("\n" + "\n".join(lines.bounds_table(self.rec)) + "\n", body)
        self.assertIn(lines.set_aside(self.rec), blocks(body))
        self.assertNotIn(self.rec.set_aside["dataset"] + " report", self.text.replace(body, ""))
        self.assertNotIn("set aside", self.text.replace(body, ""))


if __name__ == "__main__":
    unittest.main()
