"""docs/study-split.md: every line that carries a number is rebuilt from the record, once and in order; the number-free lines it holds are word for word; and the page's form."""

import re
import unittest
from fractions import Fraction

from tests.study_split import lines as builder
from tests.study_split.common import PAGE, RECORD, ROOT, SCRIPT, load_script

record = load_script()
HEADERS = ["# Study split", "## Plans", "## Folds", "## Rule", "## First experiment", "## Second experiment",
           "## Choices", "## Limits", "## Files"]


def blocks(text):
    """The page as (kind, lines): header, table, list, code or paragraph."""
    out, cur, kind, fence = [], [], None, False
    for line in text.split("\n"):
        if line.startswith("```"):
            if fence:
                out.append(("code", cur + [line]))
                cur, kind, fence = [], None, False
            else:
                if cur:
                    out.append((kind, cur))
                cur, kind, fence = [line], "code", True
            continue
        if fence:
            cur.append(line)
            continue
        if not line.strip():
            if cur:
                out.append((kind, cur))
            cur, kind = [], None
            continue
        k = "header" if line.startswith("#") else "table" if line.startswith("|") else \
            "list" if line.startswith("- ") else "paragraph"
        if kind not in (None, k) or k == "header":
            if cur:
                out.append((kind, cur))
            cur = []
        cur.append(line)
        kind = k
    if cur:
        out.append((kind, cur))
    return out


def words(text):
    return len(text.split())


def cells(row):
    return [c.strip() for c in row.strip().strip("|").split("|")]


class Page(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = PAGE.read_text(encoding="utf-8")
        cls.lines = cls.text.split("\n")
        cls.blocks = blocks(cls.text)
        cls.rec = builder.Record()
        cls.built = builder.lines(cls.rec)
        cls.held = builder.held(cls.rec)

    def test_every_numbered_line_is_built(self):
        built = set(self.built)
        for line in self.lines:
            if builder.carries_number(line) and not re.fullmatch(r"\|(---\|)+", line):
                self.assertIn(line, built)

    def test_every_built_line_is_on_the_page_once(self):
        numbered = [x for x in self.built if builder.carries_number(x)]
        self.assertEqual(len(numbered), len(set(numbered)), "the builder gives a line twice")
        for line in numbered:
            self.assertEqual(self.lines.count(line), 1, line)
        for line in self.built:
            self.assertIn(line, self.lines)

    def test_built_lines_in_page_order(self):
        at = -1
        for x in self.built:
            self.assertIn(x, self.lines[at + 1:], x)
            at = self.lines.index(x, at + 1)

    def test_number_free_lines_held_word_for_word(self):
        self.assertEqual(len(self.held), len(set(self.held)))
        for line in self.held:
            self.assertEqual(self.lines.count(line), 1, line)
            self.assertFalse(builder.carries_number(line), line)
        opening = self.lines[:self.lines.index("## Plans")]
        self.assertIn(dict(builder.groups(self.rec))["open_denial"][0], opening)

    def test_every_line_is_built_or_held(self):
        """Each paragraph, list item and table line outside the code block is built from the record or held."""
        known = set(self.built) | set(self.held)
        for kind, ls in self.blocks:
            if kind in ("header", "code"):
                continue
            for line in ls:
                if not re.fullmatch(r"\|(---\|)+", line):
                    self.assertIn(line, known)

    def test_files_count_the_held_lines(self):
        item = self.lines[self.lines.index("## Files"):]
        said = [int(m) for x in item for m in re.findall(r"holds its other (\d+) paragraphs and list items", x)]
        self.assertEqual(said, [len(self.held)])
        self.assertEqual(len(self.held), builder.held_count(self.rec))

    def test_number_words_count_as_numbers(self):
        for text in ("five study folds", "Two placement rules", "the first of twenty", "thirty", "one dataset"):
            self.assertTrue(builder.carries_number(text), text)
        for text in ("the first experiment", "someone", "anyone", "a third arm", "sha256"):
            self.assertFalse(builder.carries_number(text), text)

    def test_limits_whole(self):
        at = self.lines.index("## Limits")
        block = next(ls for kind, ls in blocks("\n".join(self.lines[at + 1:])) if kind == "list")
        self.assertEqual(block, builder.limits(self.rec))
        self.assertEqual(len(block), 14)

    def test_headers(self):
        heads = [b[1][0] for b in self.blocks if b[0] == "header"]
        self.assertEqual(heads, HEADERS)
        self.assertLessEqual(len(heads), 12)
        for h in heads:
            title = h.lstrip("#").strip()
            self.assertIn(words(title), (1, 2), h)
            self.assertNotRegex(title, r"^(What|How|Why)\b")

    def test_bold_lead_ins(self):
        leads = re.findall(r"^(?:- )?\*\*([^*]+)\*\*", self.text, re.M)
        self.assertTrue(leads)
        for lead in leads:
            self.assertTrue(lead.endswith("."), lead)
            self.assertIn(words(lead), (1, 2), lead)
            self.assertNotRegex(lead, r"^(What|How|Why)\b")

    def test_lists(self):
        items = [x for b in self.blocks if b[0] == "list" for x in b[1]]
        self.assertLessEqual(len(items), 30)
        for x in items:
            self.assertLessEqual(words(x[2:]), 45, x)

    def test_paragraphs(self):
        for kind, ls in self.blocks:
            if kind == "paragraph":
                self.assertEqual(len(ls), 1, "a paragraph is one line")
                self.assertLessEqual(words(ls[0]), 80, ls[0])

    def test_tables(self):
        for kind, ls in self.blocks:
            if kind != "table":
                continue
            head = cells(ls[0])
            self.assertRegex(ls[1], r"^\|(---\|)+$")
            self.assertEqual(ls[1].count("|") - 1, len(head))
            for h in head:
                self.assertIn(words(h), (1, 2), h)
            for row in ls[2:]:
                self.assertEqual(len(cells(row)), len(head), row)
                for c in cells(row):
                    self.assertLessEqual(words(c), 30, c)

    def test_recall_specificity_and_loss_side_by_side_per_learner(self):
        seen = 0
        for kind, ls in self.blocks:
            head = cells(ls[0]) if kind == "table" else []
            if "Recall" in head:
                seen += 1
                at = head.index("Recall")
                self.assertEqual(head[at + 1], "Specificity", ls[0])
                self.assertIn(head[at + 2], ("Loss", "Log loss"), ls[0])
                self.assertEqual(head[0], "Learner", ls[0])
                self.assertEqual(list(dict.fromkeys(cells(r)[0] for r in ls[2:])), list(builder.LEARNERS))
        self.assertEqual(seen, 7)

    def test_no_line_numbers(self):
        self.assertNotRegex(self.text, r"\.(py|rs|md|yml|tsv|json):\d")
        self.assertNotRegex(self.text, r"\b[Ll]ines? \d")

    def test_no_figure_for_a_single_dataset(self):
        # The one dataset named is the one the fold check stopped, with its genome count.
        for r in self.rec.datasets:
            self.assertNotIn(r["dataset"], self.text)
        self.assertEqual(re.findall(r"\b[a-z]+_[a-z]+__[a-z_]+\b", self.text), [self.rec.s("left_by_folds")])

    def test_scm_loss_said_where_scm_first_appears(self):
        first = next(ls[0] for kind, ls in self.blocks if kind == "paragraph" and re.search(r"\bscm\b", ls[0]))
        self.assertIn("this page reads scm's loss as its error rate", first)

    def test_number_free_claims(self):
        """Sentences that carry no number but state a fact of the record, rebuilt from it."""
        rec = self.rec
        # "Matched plans equal the study draws on every measure, as they copy them."
        for m in ("largest", "smallest", "least_training", "spread", "departure"):
            self.assertEqual([rec.shapes[("matched", m)][k] for k in ("median", "low", "high")],
                             [rec.shapes[("study", m)][k] for k in ("median", "low", "high")], m)
        self.assertIn("Matched plans equal the study draws on every measure, as they copy them.", self.text)
        # "Study folds differ from random folds in size and in class balance."
        for m in ("largest", "spread", "departure"):
            self.assertGreater(Fraction(rec.shapes[("study-d0", m)]["median"]),
                               Fraction(rec.shapes[("random-d0", m)]["median"]), m)
        self.assertIn("Study folds differ from random folds in size and in class balance.", self.text)
        # "... and four of them use a more uneven rule than the first."
        self.assertEqual(rec.settings["study_rule"], "d0:A,d1:B,d2:B,d3:B,d4:B")
        for m in ("largest", "spread", "departure"):
            self.assertGreater(Fraction(rec.shapes[("study", m)]["median"]),
                               Fraction(rec.shapes[("study-d0", m)]["median"]), m)
        # "scm's figures are its error rate ...; both give the same readings." The counts of datasets
        # agree too: the script reads scm's log loss in exact form from its error rate.
        pairs = 0
        for (exp, l, fam, inst, m), r in rec.contrasts.items():
            if l == "scm" and m == "log_loss":
                e = rec.contrasts[(exp, l, fam, inst, "error_rate")]
                self.assertEqual([r[k] for k in ("reading", "worse", "equal", "better")],
                                 [e[k] for k in ("reading", "worse", "equal", "better")], (exp, fam, inst))
                pairs += 1
        self.assertEqual(pairs, len(record.FIRST) + sum(len(record.instances(f)) + 1 for f in record.SECOND))
        self.assertIn("- scm's figures are its error rate, apart from the log loss in the first experiment's two tables, "
                      "which the first plan fixed; both give the same readings.", self.lines)
        # "No test genome in a study plan has a study mate in training."
        self.assertEqual({r["study_mates"] for r in rec.plans if r["kind"] == "study"}, {"0"})
        self.assertIn("No test genome in a study plan has a study mate in training.", self.text)
        # "Study plans also keep more lineage-cluster mates out of training."
        share = {k: Fraction(rec.mates[(k, "mean")]["lineage_share"]) for k in ("study", "matched", "random")}
        self.assertLess(share["study"], share["matched"])
        self.assertLess(share["matched"], share["random"])
        # "Each experiment hashed its plan before its first fit." and the order of the variant count.
        o = rec.order
        for exp in ("first", "second"):
            self.assertLess(o[(exp, "hash_time")], o[(exp, "first_fit_end")])
        self.assertLess(o[("first", "fold_files_last")], o[("first", "variant_table")])
        self.assertLess(o[("first", "variant_table")], o[("first", "hash_time")])
        self.assertLess(o[("second", "plan_builder")], o[("second", "stated_rule_choice")])
        self.assertIn("Each experiment hashed its plan before its first fit.", self.text)
        # Five learners, five folds, five draws; scm scores only 0 or 1.
        self.assertEqual((len(record.LEARNERS), rec.settings["folds"], rec.settings["draws"]), (5, "5", "5"))
        self.assertEqual(rec.checks["scm_fits_scores_0_or_1"], len(rec.plans))
        self.assertIn("scm scores only 0 or 1", self.text)
        # The interval rests on species resampling only: the script resamples nothing else.
        self.assertIn("The interval holds species resampling only.", self.text)

    def test_files_named_exist(self):
        for path in ("reference/study-split", "scripts/study-split/record.py", "tests/study_split"):
            self.assertIn(f"`{path}", self.text)
            self.assertTrue((ROOT / path).exists(), path)
        for name in re.findall(r"`([\w./-]+\.(?:tsv|md|sha256))`", self.text):
            self.assertTrue((RECORD / name).exists(), name)
        self.assertIn("SHA256SUMS", self.text)
        self.assertTrue((ROOT / "docs" / "fold-file.md").exists())

    def test_commands_match_the_script(self):
        code = "\n".join(x for kind, ls in self.blocks if kind == "code" for x in ls)
        source = SCRIPT.read_text(encoding="utf-8")
        for flag in re.findall(r"(--[a-z-]+)", code):
            self.assertIn(f'"{flag}"', source, flag)
        self.assertIn("python3 -I scripts/study-split/record.py --inputs <dir> --check reference/study-split", code)
        self.assertIn("python3 -I -m unittest discover -s tests/study_split -t .", code)


if __name__ == "__main__":
    unittest.main()
