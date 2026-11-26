"""reference/input-arms-m/ against the pages under docs/, in the record's three states.

The record stands in one of three states, by the tables it holds (m_page.py's state). m_page.py writes two pages
whole from it, and this test compares each with its bytes:

- none: README.md and first-run-times.tsv alone; neither written page exists; the result page says the datasets of M
  wait, and input-arms.md that the run has not happened;
- floor: the floor's two tables too; docs/input-arms-m-left-out.md is the bytes `m_page.py left` gives, and
  docs/input-arms-m-count.md does not exist; input-arms.md names the left-out page;
- count: all 7 tables; both written pages are the bytes `m_page.py left` and `m_page.py page` give; the result page
  names the count page; and every file equals its SHA256SUMS line, which SHA256SUMS lists once each.

Every other page under docs/ is read line by line as cells: the line stripped, one outer pipe on each side dropped
when there, the rest split on pipes and each cell stripped. No line of 2 or more cells may be the header of a table of
this run, in any state; after the floor, nor a row of the left-out or species table; after the count, nor a row of any
table print gives, nor of any count in the page's forms on any split and measure, which hold the written pages' rows.
No check reads Markdown.

The committed tree is checked in the state its record gives. The other states are proven on copies of the tree's
pages but the two written pages, with a record that starts as the committed one before the run and SYNTHETIC tables
written by m_reading.py's own forms (tests/c3_inputs/m_tables.py), run through the copy's m_page.py as the lead will
run it.
"""

import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_tables

m_page = m_tables.m_page
RECORD = Path("reference") / "input-arms-m"
ARMS = "docs/input-arms.md"
RESULT = "docs/input-arms-result.md"
LEFT, COUNT = m_page.LEFT_PAGE, m_page.COUNT_PAGE
# The written pages each state may hold, which the record gives, and the command that gives each.
WRITTEN = {"none": (), "floor": (LEFT,), "count": (LEFT, COUNT)}
GIVES = {LEFT: "left", COUNT: "page"}
RULE = re.compile(r":?-+:?")
BUILT = {}


class Off(Exception):
    """The tree does not stand as the state its record gives."""


def cells(line):
    """A line as cells: stripped, one outer pipe on each side dropped when there, split on pipes, each cell stripped."""
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return tuple(cell.strip() for cell in line.split("|"))


def table_cells(lines):
    """The cells of each line of these tables but their rule rows."""
    return {cells(line) for line in lines if line.startswith("|") and not all(RULE.fullmatch(c) for c in cells(line))}


def forbidden(record, m_set, state):
    """(the headers' cells, the rows' cells) that no page but the written ones may hold in this state."""
    heads = {tuple(h) for h in m_page.headers()}
    if state == "none":
        return heads, set()
    printed = m_page.tables(record, m_set)
    floor_rows = table_cells(printed["left-out"] + printed["species"]) - heads
    if state == "floor":
        return heads, floor_rows
    count = m_page.read_tables(record, ["m-count.tsv"])["m-count.tsv"]
    lines = [line for name, t in printed.items() for line in t]
    for c in m_page.COUNTS:
        lines += m_page.print_count(count, c) + m_page.print_marked(count, c)
        for split, measure in sorted({(r["split"], r["measure"]) for _, r in count.rows}):
            lines += m_page.print_other(count, c, split, measure)
    count_rows = table_cells(lines) - heads
    return heads, floor_rows | count_rows


def check_tree(root, lists=None):
    """The state of the tree at <root>, after every check of that state; Off names the first check that fails."""
    root = Path(root)
    record = root / RECORD
    m_set = m_page.m_datasets(lists or root / "scripts" / "c3-inputs" / "lists.tsv")
    state = m_page.state(record)
    for name in (LEFT, COUNT):
        if (root / name).exists() and name not in WRITTEN[state]:
            raise Off(f"{name} exists, and the record in the state {state} cannot give it")
    pages = {p.relative_to(root).as_posix(): p.read_text(encoding="utf-8")
             for p in sorted((root / "docs").glob("*.md")) if p.relative_to(root).as_posix() not in GIVES}
    arms, result = pages[ARMS], pages[RESULT]
    if state == "none":
        have = sorted(p.name for p in record.iterdir())
        if have != ["README.md", "first-run-times.tsv"]:
            raise Off(f"the record holds {', '.join(have)}, not README.md and first-run-times.tsv alone")
        if f"The {len(m_set)} datasets of M wait" not in result:
            raise Off(f"{RESULT} does not say the {len(m_set)} datasets of M wait")
        if "The run has not happened" not in arms:
            raise Off(f"{ARMS} does not say the run has not happened")
    try:
        given = {}
        for name in WRITTEN[state]:
            given[name] = m_page.left_page(record, m_set) if name == LEFT else m_page.count_page(record, m_set)
        heads, rows = forbidden(record, m_set, state)
    except m_page.Stop as e:
        raise Off(f"m_page: stopped: {e}") from None
    for name, want in given.items():
        path = root / name
        if not path.is_file():
            raise Off(f"{name} does not exist, and the record in the state {state} gives it")
        have = path.read_bytes()
        if have != want:
            diff = m_page.first_difference(want, have) or (1, "differs")
            raise Off(f"{name} {m_page.difference_words(*diff, GIVES[name])}")
    if state != "none" and f"`{LEFT}`" not in arms:
        raise Off(f"{ARMS} does not name `{LEFT}`")
    if state == "count" and f"`{COUNT}`" not in result:
        raise Off(f"{RESULT} does not name `{COUNT}`")
    for name, text in pages.items():
        for n, line in enumerate(text.split("\n"), 1):
            got = cells(line)
            if len(got) < 2:
                continue
            if got in heads:
                raise Off(f"{name} line {n} holds the cells of a header of this run's tables")
            if got in rows:
                raise Off(f"{name} line {n} holds the cells of a row of this run's tables")
    if state != "none":
        sums = [line.partition("  ") for line in (record / "SHA256SUMS").read_text(encoding="utf-8").splitlines()]
        listed = {name: digest for digest, _, name in sums}
        files = sorted(p.name for p in record.iterdir() if p.name not in ("SHA256SUMS", "README.md"))
        if sorted(listed) != files or len(sums) != len(listed):
            raise Off("SHA256SUMS does not list exactly the record's files but itself and README.md, once each")
        for name in files:
            if hashlib.sha256((record / name).read_bytes()).hexdigest() != listed[name]:
                raise Off(f"{name} differs from its SHA256SUMS line")
    return state


def lay_out(root, source):
    """A copy at <root>: <source>'s pages but the two written ones, its record's README.md and first-run-times.tsv
    alone, its m_page.py and pack.py, and a SYNTHETIC lists.tsv, whose path it returns."""
    (root / "docs").mkdir(parents=True)
    for p in (Path(source) / "docs").glob("*.md"):
        if f"docs/{p.name}" not in GIVES:
            shutil.copy(p, root / "docs" / p.name)
    record = root / RECORD
    record.mkdir(parents=True)
    for name in ("README.md", "first-run-times.tsv"):
        shutil.copy(Path(source) / RECORD / name, record / name)
    kit = root / "scripts" / "c3-inputs"
    kit.mkdir(parents=True)
    for name in ("m_page.py", "pack.py"):
        shutil.copy(Path(source) / "scripts" / "c3-inputs" / name, kit / name)
    return m_tables.lists_file(kit / "lists.tsv")


def setUpModule():
    BUILT["tmp"] = tempfile.TemporaryDirectory()
    top = Path(BUILT["tmp"].name)
    BUILT["seven"] = m_tables.seven(top / "default")
    BUILT["all in"] = m_tables.seven(top / "all-in", m_tables.ALL_IN)


def tearDownModule():
    BUILT.pop("tmp").cleanup()


def text(lines):
    return ("\n".join(lines) + "\n").encode("utf-8")


class Copy(unittest.TestCase):
    """Each test on its own copy: the tree's pages but the written ones, the record as it stands before the run, the
    tree's m_page.py and pack.py, and a SYNTHETIC lists.tsv."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "tree"
        self.record = self.root / RECORD
        self.lists = lay_out(self.root, common.ROOT)

    def again(self, name):
        """A copy made from this copy as setUp makes one from the tree, which it then stands for."""
        source = self.root
        self.root = Path(self.tmp.name) / name
        self.record = self.root / RECORD
        self.lists = lay_out(self.root, source)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *argv):
        """(exit, the bytes on standard output) of the copy's m_page.py, run as the lead runs it."""
        done = subprocess.run([sys.executable, "-I", "-B", str(self.root / "scripts" / "c3-inputs" / "m_page.py"),
                               *map(str, argv)], cwd=self.root, capture_output=True, timeout=300)
        return done.returncode, done.stdout

    def page(self, name):
        return self.root / name

    def add(self, name, lines):
        """These lines put at the end of a page; returns the line number of the first."""
        path = self.page(name)
        kept = path.read_text(encoding="utf-8")
        if not kept.endswith("\n"):
            kept += "\n"
        path.write_text(kept + "\n" + "\n".join(lines) + "\n", encoding="utf-8")
        return len(kept.split("\n")) + 1

    def written(self, name):
        rc, out = self.cli(GIVES[name], self.record)
        self.assertEqual(rc, 0, out)
        self.page(name).write_bytes(out)
        return out

    def to_floor(self, which="seven"):
        """floor through the copy's CLI, then the left-out page written from what left gives. Returns its bytes."""
        self.assertEqual(self.cli("floor", BUILT[which], self.record),
                         (0, b"floor: 2 tables copied; SHA256SUMS lists 3 files\n"))
        return self.written(LEFT)

    def to_count(self, which="seven"):
        """The floor, then record through the copy's CLI and the count page written from what page gives. Returns its
        bytes."""
        self.to_floor(which)
        self.assertEqual(self.cli("record", BUILT[which], self.record),
                         (0, b"record: 5 tables copied; SHA256SUMS lists 8 files\n"))
        return self.written(COUNT)

    def fresh(self):
        """A new copy inside one test, for a rule tried several ways."""
        self.tearDown()
        self.setUp()

    def off(self, pattern):
        with self.assertRaisesRegex(Off, pattern):
            check_tree(self.root, self.lists)


def forms(line):
    """One table line as a copy may give it: as it is, respaced, without its outer pipes or with one alone, indented,
    and with tabs about its pipes."""
    inner = line[2:-2]
    return [line, "|" + "|".join(inner.split(" | ")) + "|", "|   " + "   |   ".join(inner.split(" | ")) + "   |",
            inner, "| " + inner, inner + " |", "  " + line + "  ", "|\t" + "\t|\t".join(inner.split(" | ")) + "\t|"]


class Cells(unittest.TestCase):
    def test_a_line_reads_as_its_cells_with_the_outer_pipes_optional_and_each_cell_stripped(self):
        for line in forms("| a b | c | `d` |"):
            self.assertEqual(cells(line), ("a b", "c", "`d`"), line)
        self.assertEqual(cells("| a | b"), ("a", "b"))
        self.assertEqual(cells("a | | b"), ("a", "", "b"))
        self.assertEqual(cells("||a|"), ("", "a"))
        self.assertEqual(cells("Words alone."), ("Words alone.",))
        self.assertEqual(cells("|"), ("",))
        self.assertEqual(table_cells(["| a | b |", "|---|---|", "| :--- | --: |", "| c | d |", "words"]),
                         {("a", "b"), ("c", "d")})


class Committed(unittest.TestCase):
    def test_the_committed_record_and_pages_stand_as_the_records_state_gives(self):
        self.assertIn(check_tree(common.ROOT), ("none", "floor", "count"))

    def test_no_header_of_this_run_is_a_header_on_the_committed_pages(self):
        # The header rows of every table on the committed pages but the written ones, by their cells, against the
        # headers of this run's 11 tables.
        heads = {tuple(h) for h in m_page.headers()}
        self.assertEqual(len(heads), 11)
        seen = 0
        for p in sorted((common.ROOT / "docs").glob("*.md")):
            if f"docs/{p.name}" in GIVES:
                continue
            lines = p.read_text(encoding="utf-8").split("\n")
            for n, line in enumerate(lines[1:], 2):
                if line.startswith("|") and all(RULE.fullmatch(c) for c in cells(line)):
                    seen += 1
                    self.assertNotIn(cells(lines[n - 2]), heads, f"{p.name} line {n - 1}")
        self.assertGreater(seen, 0)


class BeforeTheRun(Copy):
    """The state none on made-up pages, so the test stands after the tree's pages change."""

    def setUp(self):
        super().setUp()
        for p in (self.root / "docs").glob("*.md"):
            p.unlink()
        self.page(ARMS).write_text("# SYNTHETIC\n\n**Not run.** The run has not happened.\n", encoding="utf-8")
        self.page(RESULT).write_text("# SYNTHETIC\n\n**Labels.** The 6 datasets of M wait for a rule.\n",
                                     encoding="utf-8")

    def test_none_passes_with_the_two_sentences_and_no_table_of_m(self):
        self.assertEqual(check_tree(self.root, self.lists), "none")

    def test_none_refuses_pages_that_no_longer_say_m_waits(self):
        self.page(RESULT).write_text("# SYNTHETIC\n\n**Labels.** The 6 datasets of M ran.\n", encoding="utf-8")
        self.off(f"^{RESULT} does not say the 6 datasets of M wait$")

    def test_none_refuses_a_page_without_the_run_not_happened(self):
        self.page(ARMS).write_text("# SYNTHETIC\n\n**Not run.** The run is done.\n", encoding="utf-8")
        self.off(f"^{ARMS} does not say the run has not happened$")

    def test_none_refuses_a_record_holding_more_than_the_readme_and_times(self):
        (self.record / "SHA256SUMS").write_text("", encoding="utf-8")
        self.off("^the record holds README.md, SHA256SUMS, first-run-times.tsv, not README.md and "
                 "first-run-times.tsv alone$")

    def test_none_refuses_either_written_page(self):
        for name in (LEFT, COUNT):
            with self.subTest(name=name):
                self.page(name).write_text("# SYNTHETIC\n", encoding="utf-8")
                self.off(f"^{name} exists, and the record in the state none cannot give it$")
                self.page(name).unlink()

    def test_none_refuses_a_header_of_this_run_on_a_page_in_any_form(self):
        for h in m_page.headers():
            for line in forms("| " + " | ".join(h) + " |"):
                with self.subTest(header=h[:2], line=line[:12]):
                    kept = self.page(RESULT).read_text(encoding="utf-8")
                    n = self.add(RESULT, [line, "|---|---|"])
                    self.off(f"^{RESULT} line {n} holds the cells of a header of this run's tables$")
                    self.page(RESULT).write_text(kept, encoding="utf-8")
        self.assertEqual(check_tree(self.root, self.lists), "none")


class FloorTaken(Copy):
    def test_floor_passes_with_the_left_out_page_written(self):
        left = self.to_floor()
        self.assertEqual(m_page.state(self.record), "floor")
        self.assertEqual(check_tree(self.root, self.lists), "floor")
        self.assertEqual(left, m_page.left_page(self.record, m_tables.M_SET))
        self.assertEqual(self.cli("onpage", self.record, self.page(LEFT)),
                         (0, f"onpage: {self.page(LEFT)} is the bytes `m_page.py left` gives\n".encode("utf-8")))

    def test_floor_with_no_dataset_left_out_passes_with_its_one_line(self):
        left = self.to_floor("all in")
        self.assertTrue(left.endswith(b"\n\n0 of the 6 datasets of M are left out of the second reading.\n"))
        self.assertEqual(check_tree(self.root, self.lists), "floor")

    def test_floor_refuses_a_left_out_page_that_does_not_exist(self):
        self.to_floor()
        self.page(LEFT).unlink()
        self.off(f"^{LEFT} does not exist, and the record in the state floor gives it$")

    def test_floor_refuses_a_page_with_a_row_changed(self):
        # A row changed, respaced, ended early, run on, or two rows in the other order: the line it names is the first
        # that differs.
        left = self.to_floor().decode("utf-8").split("\n")[:-1]
        gives = "the bytes `m_page.py left` gives"
        for lines, words in (
                ([*left[:10], left[10].replace("| 30 R, 40 S |", "| 31 R, 40 S |")], f"line 11 differs from {gives}"),
                ([*left[:8], left[8] + " ", *left[9:]], f"line 9 differs from {gives}"),
                ([*left[:8], left[8].replace(" | ", "|"), *left[9:]], f"line 9 differs from {gives}"),
                (left[:10], f"ends before line 11 of {gives}"),
                (left + left[8:9], f"runs on at line 12, past the end of {gives}"),
                ([*left[:8], left[9], left[8], left[10]], f"line 9 differs from {gives}")):
            with self.subTest(words=words, lines=len(lines)):
                self.page(LEFT).write_bytes(text(lines))
                self.off(f"^{re.escape(LEFT)} {re.escape(words)}$")
        self.page(LEFT).write_bytes(text(left))
        self.assertEqual(check_tree(self.root, self.lists), "floor")

    def test_floor_refuses_the_count_page(self):
        left = self.to_floor()
        self.page(COUNT).write_bytes(left)
        self.off(f"^{COUNT} exists, and the record in the state floor cannot give it$")

    def test_floor_refuses_a_design_page_that_does_not_name_the_left_out_page(self):
        self.to_floor()
        path = self.page(ARMS)
        path.write_text(path.read_text(encoding="utf-8").replace(f"`{LEFT}`", "the left-out page"), encoding="utf-8")
        self.off(f"^{ARMS} does not name `{LEFT}`$")

    def test_floor_reads_each_line_of_another_page_as_cells(self):
        # A row of the left-out table in each form a copy may give it, on each page but the written ones.
        left = self.to_floor().decode("utf-8").split("\n")
        for name in (ARMS, RESULT, "docs/ci.md"):
            for line in forms(left[8]):
                with self.subTest(name=name, line=line):
                    path = self.page(name)
                    kept = path.read_text(encoding="utf-8")
                    n = self.add(name, [line])
                    self.off(f"^{name} line {n} holds the cells of a row of this run's tables$")
                    path.write_text(kept, encoding="utf-8")
        self.assertEqual(check_tree(self.root, self.lists), "floor")

    def test_floor_refuses_a_count_table_on_any_page(self):
        # Each header of this run, the left-out table's among them, and each row of the left-out and species tables.
        self.to_floor()
        printed = m_page.tables(self.record, m_tables.M_SET)
        heads = ["| " + " | ".join(h) + " |" for h in m_page.headers()]
        for line, kind in [(h, "header") for h in heads] + [(r, "row") for r in printed["left-out"][2:]
                                                              + printed["species"][2:]]:
            for name in (ARMS, RESULT):
                with self.subTest(line=line[:30], name=name):
                    path = self.page(name)
                    kept = path.read_text(encoding="utf-8")
                    n = self.add(name, [line])
                    self.off(f"^{name} line {n} holds the cells of a {kind} of this run's tables$")
                    path.write_text(kept, encoding="utf-8")

    def test_floor_refuses_a_record_that_breaks_a_gate(self):
        self.to_floor()
        (self.record / "m-left-out.tsv").write_bytes((self.record / "m-left-out.tsv").read_bytes() + b"\n")
        self.off("^m_page: stopped: m-left-out.tsv differs from its SHA256SUMS line$")


class CountTaken(Copy):
    def test_count_passes_with_both_written_pages(self):
        page = self.to_count()
        self.assertEqual(m_page.state(self.record), "count")
        self.assertEqual(check_tree(self.root, self.lists), "count")
        self.assertEqual(page, m_page.count_page(self.record, m_tables.M_SET))
        self.assertEqual(self.page(LEFT).read_bytes(), m_page.left_page(self.record, m_tables.M_SET))

    def test_a_copy_made_from_a_tree_in_each_state_stands_in_each_state(self):
        # A tree taken to the floor or the count, as the lead leaves it after each step; a copy made from it as setUp
        # makes one from the committed tree takes its own record, from other tables, to the floor or the count.
        for first in ("floor", "count"):
            for then in ("floor", "count"):
                with self.subTest(tree=first, copy=then):
                    self.fresh()
                    (self.to_floor if first == "floor" else self.to_count)()
                    self.assertEqual(check_tree(self.root), first)
                    self.again("copy")
                    (self.to_floor if then == "floor" else self.to_count)("all in")
                    self.assertEqual(check_tree(self.root, self.lists), then)

    def test_count_refuses_a_written_page_that_does_not_exist(self):
        for name in (LEFT, COUNT):
            with self.subTest(name=name):
                self.fresh()
                self.to_count()
                self.page(name).unlink()
                self.off(f"^{name} does not exist, and the record in the state count gives it$")

    def test_count_refuses_a_count_page_unlike_the_bytes_page_gives(self):
        page = self.to_count().decode("utf-8").split("\n")[:-1]
        gives = "the bytes `m_page.py page` gives"
        n = next(i for i, line in enumerate(page) if line.startswith("| first | unitigs | `scm` |"))
        changed = page[n].replace("| first |", "| second |", 1)
        for lines, words in (([*page[:n], changed, *page[n + 1:]], f"line {n + 1} differs from {gives}"),
                             ([*page[:n], page[n] + " ", *page[n + 1:]], f"line {n + 1} differs from {gives}"),
                             (page[:-2], f"ends before line {len(page) - 1} of {gives}"),
                             (page + page[n:n + 1], f"runs on at line {len(page) + 1}, past the end of {gives}"),
                             (page + page[5:], f"runs on at line {len(page) + 1}, past the end of {gives}")):
            with self.subTest(words=words):
                self.page(COUNT).write_bytes(text(lines))
                self.off(f"^{re.escape(COUNT)} {re.escape(words)}$")
        # The 8 tables in another order: the first line that differs is the first run-in line moved.
        blocks = "\n".join(page).split("\n\n")
        self.assertEqual(len(blocks), 20)
        self.page(COUNT).write_bytes(("\n\n".join(blocks[:3] + blocks[5:7] + blocks[3:5] + blocks[7:]) + "\n").encode())
        self.off(f"^{re.escape(COUNT)} line 7 differs from {re.escape(gives)}$")

    def test_count_refuses_a_result_page_that_does_not_name_the_count_page(self):
        self.to_count()
        path = self.page(RESULT)
        path.write_text(path.read_text(encoding="utf-8").replace(f"`{COUNT}`", "the count page"), encoding="utf-8")
        self.off(f"^{RESULT} does not name `{COUNT}`$")

    def test_count_refuses_a_table_that_stays_in_the_record_on_any_page(self):
        # The header of each table print gives, on either page.
        self.to_count()
        printed = m_page.tables(self.record, m_tables.M_SET)
        for name, lines in printed.items():
            for page in (RESULT, ARMS):
                with self.subTest(name=name, page=page):
                    path = self.page(page)
                    kept = path.read_text(encoding="utf-8")
                    header = next(i for i, line in enumerate(lines) if line.startswith("|"))
                    n = self.add(page, lines[header:])
                    self.off(f"^{page} line {n} holds the cells of a header of this run's tables$")
                    path.write_text(kept, encoding="utf-8")

    def test_count_refuses_a_row_of_a_table_that_stays_in_the_record_on_any_page(self):
        # A row of each table print gives with no header, respaced, on either page; and the page's own form of a count
        # on another split or measure, or of the first run's count, which the count page does not give.
        self.to_count()
        printed = m_page.tables(self.record, m_tables.M_SET)
        for name, lines in printed.items():
            for page in (RESULT, ARMS):
                with self.subTest(name=name, page=page):
                    path = self.page(page)
                    kept = path.read_text(encoding="utf-8")
                    n = self.add(page, [forms(lines[-1])[1]])
                    self.off(f"^{page} line {n} holds the cells of a row of this run's tables$")
                    path.write_text(kept, encoding="utf-8")
        count = m_page.read_tables(self.record, ["m-count.tsv"])["m-count.tsv"]
        for lines in (m_page.print_other(count, "first run", "study", "balanced_accuracy"),
                      m_page.print_other(count, "M", "lineage", "specificity"),
                      m_page.print_other(count, "M", "study", "log_loss"),
                      m_page.print_count(count, "first run"), m_page.print_marked(count, "M")):
            with self.subTest(row=lines[2][:40]):
                path = self.page(RESULT)
                kept = path.read_text(encoding="utf-8")
                n = self.add(RESULT, ["**More.** Another count.", "", lines[-1]])
                self.off(f"^{RESULT} line {n + 2} holds the cells of a row of this run's tables$")
                path.write_text(kept, encoding="utf-8")

    def test_count_refuses_a_row_of_either_written_page_on_any_other_page(self):
        # Every table line of the two written pages is among the cells no other page may hold.
        page = self.to_count().decode("utf-8").split("\n")
        left = self.page(LEFT).read_text(encoding="utf-8").split("\n")
        heads, rows = forbidden(self.record, m_tables.M_SET, "count")
        lines = [line for line in page + left if line.startswith("|")]
        self.assertEqual(len(lines), 363 + 8 * 2 + 3 + 2)
        self.assertLessEqual(table_cells(lines), heads | rows)
        for line in (page[-4], left[-2]):
            path = self.page(RESULT)
            kept = path.read_text(encoding="utf-8")
            n = self.add(RESULT, [line])
            self.off(f"^{RESULT} line {n} holds the cells of a row of this run's tables$")
            path.write_text(kept, encoding="utf-8")

    def test_count_refuses_a_left_out_table_changed_on_the_design_page(self):
        self.to_count()
        path = self.page(LEFT)
        path.write_bytes(path.read_bytes().replace(b"| 30 R, 40 S |", b"| 31 R, 40 S |"))
        self.off(f"^{re.escape(LEFT)} line 11 differs from the bytes `m_page.py left` gives$")

    def test_count_refuses_a_file_unlike_its_sha256sums_line(self):
        self.to_count()
        (self.record / "m-both.tsv").write_bytes((self.record / "m-both.tsv").read_bytes() + b"\n")
        self.off("^m_page: stopped: m-both.tsv differs from its SHA256SUMS line$")

    def test_count_checks_the_sums_list_itself_as_well(self):
        # The record's own file check is turned off here, so the test's list check is the one that stops: a file the
        # list lacks, a file twice and a name of no file.
        for change in (lambda lines: [l for l in lines if not l.endswith("  m-size.tsv")],
                       lambda lines: lines + [lines[0]],
                       lambda lines: lines + ["0" * 64 + "  m-other.tsv"]):
            with self.subTest(change=change):
                self.fresh()
                self.to_count()
                path = self.record / "SHA256SUMS"
                path.write_text("".join(l + "\n" for l in change(path.read_text().splitlines())))
                gate, m_page.gate_files = m_page.gate_files, lambda record: None
                try:
                    self.off("^SHA256SUMS does not list exactly the record's files but itself and README.md, once "
                             "each$")
                finally:
                    m_page.gate_files = gate

    def test_count_checks_each_sum_itself_as_well(self):
        # The record's own file check is turned off here, so the test's sum check is the one that stops; the times file
        # is one m_page.py sums and reads no further.
        for upto in ("floor", "count"):
            with self.subTest(upto=upto):
                self.fresh()
                (self.to_floor if upto == "floor" else self.to_count)()
                path = self.record / "first-run-times.tsv"
                path.write_bytes(path.read_bytes() + b"\n")
                gate, m_page.gate_files = m_page.gate_files, lambda record: None
                try:
                    self.off("^first-run-times.tsv differs from its SHA256SUMS line$")
                finally:
                    m_page.gate_files = gate


if __name__ == "__main__":
    unittest.main()
