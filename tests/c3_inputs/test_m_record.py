"""reference/input-arms-m/ against docs/input-arms.md and docs/input-arms-result.md, in the record's three states.

The record stands in one of three states, by the tables it holds (m_page.py's state):

- none: README.md and first-run-times.tsv alone; the result page says the datasets of M wait, and input-arms.md that
  the run has not happened;
- floor: the floor's two tables too; input-arms.md holds the left-out table as `m_page.py print` gives it and names
  `reference/input-arms-m/m-left-out.tsv`, and no count table stands on any page;
- count: all 7 tables; input-arms.md still holds the left-out table, the result page holds the blocks
  `m_page.py page` gives, each whole, as a block of its own and in its order, with each block's first line once; no
  page holds the header or a row of a table that stays in the record, in print's form or the page's; and every file
  equals its SHA256SUMS line, which SHA256SUMS lists once each.

Before the count only the header lines tell a count table: no record exists then to compare a row with.

The committed tree is checked in the state its record gives. The other states are proven on copies of the tree's
pages, with a record that starts as the committed one before the run and SYNTHETIC tables written by m_reading.py's
own forms (tests/c3_inputs/m_tables.py), run through the copy's m_page.py as the lead will run it. Every line the
committed record already gives, by print and, after the count, by page, is taken off the copied pages first, so the
copies stand whatever the tree's state. A committed file unlike its SHA256SUMS line fails the test of the committed
tree, which names it; the copies are still made.
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
LEFT_PATH = "`reference/input-arms-m/m-left-out.tsv`"
BUILT = {}


class Off(Exception):
    """The tree does not stand as the state its record gives."""


def check_tree(root, lists=None):
    """The state of the tree at <root>, after every check of that state; Off names the first check that fails."""
    root = Path(root)
    record = root / RECORD
    m_set = m_page.m_datasets(lists or root / "scripts" / "c3-inputs" / "lists.tsv")
    pages = {p.relative_to(root).as_posix(): p.read_text(encoding="utf-8")
             for p in sorted((root / "docs").glob("*.md"))}
    arms, result = pages[ARMS], pages[RESULT]
    state = m_page.state(record)
    if state in ("none", "floor"):
        heads = m_page.count_headers()
        for name, text in pages.items():
            if set(text.split("\n")) & set(heads):
                raise Off(f"{name} holds a count table before the count")
    if state == "none":
        have = sorted(p.name for p in record.iterdir())
        if have != ["README.md", "first-run-times.tsv"]:
            raise Off(f"the record holds {', '.join(have)}, not README.md and first-run-times.tsv alone")
        if f"The {len(m_set)} datasets of M wait" not in result:
            raise Off(f"{RESULT} does not say the {len(m_set)} datasets of M wait")
        if "The run has not happened" not in arms:
            raise Off(f"{ARMS} does not say the run has not happened")
        return state
    try:
        m_page.tables(record, m_set)
    except m_page.Stop as e:
        raise Off(f"m_page: stopped: {e}") from None
    ok, line = m_page.onpage(record, root / ARMS, m_set)
    if not ok:
        raise Off(line)
    if LEFT_PATH not in arms:
        raise Off(f"{ARMS} does not name {LEFT_PATH}")
    if state == "floor":
        return state
    heads = m_page.record_headers()
    for name, text in pages.items():
        if set(text.split("\n")) & set(heads):
            raise Off(f"{name} holds a table that stays in the record")
    blocks = m_page.page(record, m_set)
    rows = record_rows(record, m_set, blocks)
    for name, text in pages.items():
        if set(text.split("\n")) & rows:
            raise Off(f"{name} holds a row of a table that stays in the record")
    last = -1
    for name, lines in blocks.items():
        k, at, why = m_page.find_block(result, lines)
        if k < len(lines):
            raise Off(f"{RESULT} lacks the block {name} as page gives it")
        if why:
            raise Off(f"{RESULT} holds the block {name} with no blank line {why} it")
        if result.split("\n").count(lines[0]) != 1:
            raise Off(f"{RESULT} holds the first line of the block {name} more than once")
        if at < last:
            raise Off(f"{RESULT} holds the block {name} before the one page gives ahead of it")
        last = at
    sums = [line.partition("  ") for line in (record / "SHA256SUMS").read_text(encoding="utf-8").splitlines()]
    listed = {name: digest for digest, _, name in sums}
    files = sorted(p.name for p in record.iterdir() if p.name not in ("SHA256SUMS", "README.md"))
    if sorted(listed) != files or len(sums) != len(listed):
        raise Off("SHA256SUMS does not list exactly the record's files but itself and README.md, once each")
    for name in files:
        if hashlib.sha256((record / name).read_bytes()).hexdigest() != listed[name]:
            raise Off(f"{name} differs from its SHA256SUMS line")
    return state


def record_rows(record, m_set, blocks):
    """The row lines of the five tables that stay in the record, as print gives them, and each count's rows on every
    split and measure in the page's own forms, less the lines of page's blocks."""
    printed = m_page.tables(record, m_set)
    rows = {line for name in ("count", "beside", "conditions", "size", "datasets") for line in printed[name][4:]}
    count = m_page.read_tables(record, ["m-count.tsv"])["m-count.tsv"]
    for c in m_page.COUNTS:
        rows |= set(m_page.print_count(count, c)[2:])
        for split, measure in {(r["split"], r["measure"]) for _, r in count.rows}:
            rows |= set(m_page.print_other(count, c, split, measure)[2:])
    return rows - {line for lines in blocks.values() for line in lines}


def real_lines(root, lists=None):
    """Every line but blank ones that the record at <root> gives: print's tables and, after the count, page's
    blocks. A record that breaks the file gate is read with that gate off, and one that breaks another gives none;
    the test of the committed tree names either."""
    record = Path(root) / RECORD
    if m_page.state(record) == "none":
        return set()
    m_set = m_page.m_datasets(lists or Path(root) / "scripts" / "c3-inputs" / "lists.tsv")

    def given():
        out = list(m_page.tables(record, m_set).values())
        if m_page.state(record) == "count":
            out += list(m_page.page(record, m_set).values())
        return {line for lines in out for line in lines if line}

    try:
        return given()
    except m_page.Stop:
        pass
    gate, m_page.gate_files = m_page.gate_files, lambda record: None
    try:
        return given()
    except m_page.Stop:
        return set()
    finally:
        m_page.gate_files = gate


def lay_out(root, source, real):
    """A copy at <root>: <source>'s pages less the lines in <real>, its record's README.md and first-run-times.tsv
    alone, its m_page.py and pack.py, and a SYNTHETIC lists.tsv, whose path it returns."""
    (root / "docs").mkdir(parents=True)
    for p in (Path(source) / "docs").glob("*.md"):
        lines = p.read_text(encoding="utf-8").split("\n")
        (root / "docs" / p.name).write_text("\n".join(line for line in lines if line not in real), encoding="utf-8")
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
    BUILT["real"] = real_lines(common.ROOT, common.KIT / "lists.tsv")


def tearDownModule():
    BUILT.pop("tmp").cleanup()


def printed_tables(out):
    """{name: lines} from the output of `m_page.py print`."""
    parts = re.split(r"^\[([a-z-]+)\]\n", out, flags=re.M)
    assert parts[0] == "", parts[0]
    return {name: body[:-2].split("\n") for name, body in zip(parts[1::2], parts[2::2])}


class Copy(unittest.TestCase):
    """Each test on its own copy: the tree's pages less every line the committed record gives by print and page, the
    record as it stands before the run, the tree's m_page.py and pack.py, and a SYNTHETIC lists.tsv."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name) / "tree"
        self.record = self.root / RECORD
        self.lists = lay_out(self.root, common.ROOT, BUILT["real"])

    def again(self, name):
        """A copy made from this copy as setUp makes one from the tree, which it then stands for."""
        source, real = self.root, real_lines(self.root)
        self.root = Path(self.tmp.name) / name
        self.record = self.root / RECORD
        self.lists = lay_out(self.root, source, real)

    def tearDown(self):
        self.tmp.cleanup()

    def cli(self, *argv):
        done = subprocess.run([sys.executable, "-I", "-B", str(self.root / "scripts" / "c3-inputs" / "m_page.py"),
                               *map(str, argv)], cwd=self.root, capture_output=True, text=True, timeout=300)
        return done.returncode, done.stdout

    def page(self, name):
        return self.root / name

    def change(self, name, old, new):
        path = self.page(name)
        text = path.read_text(encoding="utf-8")
        self.assertEqual(text.count(old), 1, old)
        path.write_text(text.replace(old, new), encoding="utf-8")

    def to_floor(self, which="seven"):
        """floor through the copy's CLI, then the left-out table put on input-arms.md after the paragraph that names
        its file. Returns print's tables."""
        self.assertEqual(self.cli("floor", BUILT[which], self.record),
                         (0, "floor: 2 tables copied; SHA256SUMS lists 3 files\n"))
        rc, out = self.cli("print", self.record)
        self.assertEqual(rc, 0)
        tables = printed_tables(out)
        path = self.page(ARMS)
        parts = path.read_text(encoding="utf-8").split("\n\n")
        i = next(k for k, part in enumerate(parts) if LEFT_PATH in part)
        path.write_text("\n\n".join(parts[:i + 1] + ["\n".join(tables["left-out"])] + parts[i + 1:]),
                        encoding="utf-8")
        return tables

    def to_count(self, which="seven"):
        """The floor, then record through the copy's CLI and what its page command prints put at the result page's
        end. Returns page's blocks."""
        self.to_floor(which)
        self.assertEqual(self.cli("record", BUILT[which], self.record),
                         (0, "record: 5 tables copied; SHA256SUMS lists 8 files\n"))
        rc, out = self.cli("page", self.record)
        self.assertEqual(rc, 0)
        blocks = m_page.page(self.record, m_tables.M_SET)
        self.assertEqual(out, "\n\n".join("\n".join(lines) for lines in blocks.values()) + "\n")
        path = self.page(RESULT)
        path.write_text(path.read_text(encoding="utf-8").rstrip("\n") + "\n\n## M, SYNTHETIC\n\n" + out,
                        encoding="utf-8")
        return blocks

    def printed(self):
        rc, out = self.cli("print", self.record)
        self.assertEqual(rc, 0)
        return printed_tables(out)

    def fresh(self):
        """A new copy inside one test, for a rule tried several ways."""
        self.tearDown()
        self.setUp()

    def off(self, pattern):
        with self.assertRaisesRegex(Off, pattern):
            check_tree(self.root, self.lists)


class Committed(unittest.TestCase):
    def test_the_committed_record_and_pages_stand_as_the_records_state_gives(self):
        self.assertIn(check_tree(common.ROOT), ("none", "floor", "count"))


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
        self.change(RESULT, "The 6 datasets of M wait", "The 6 datasets of M ran")
        self.off(f"^{RESULT} does not say the 6 datasets of M wait$")

    def test_none_refuses_a_page_without_the_run_not_happened(self):
        self.change(ARMS, "The run has not happened.", "The run is done.")
        self.off(f"^{ARMS} does not say the run has not happened$")

    def test_none_refuses_a_record_holding_more_than_the_readme_and_times(self):
        (self.record / "SHA256SUMS").write_text("", encoding="utf-8")
        self.off("^the record holds README.md, SHA256SUMS, first-run-times.tsv, not README.md and "
                 "first-run-times.tsv alone$")

    def test_none_refuses_a_count_table_on_a_page(self):
        self.change(RESULT, "**Labels.**", "| Reading | Arm | Learner | Entered | Better | Worse | Same |\n\n**Labels.**")
        self.off(f"^{RESULT} holds a count table before the count$")


class FloorTaken(Copy):
    def test_floor_passes_with_the_left_out_table_on_the_page(self):
        tables = self.to_floor()
        self.assertEqual(m_page.state(self.record), "floor")
        self.assertEqual(check_tree(self.root, self.lists), "floor")
        self.assertEqual(list(tables), ["left-out", "species", "count"])
        self.assertEqual(tables["count"], [m_page.NOT_TAKEN])
        self.assertEqual(self.cli("onpage", self.record, self.page(ARMS)),
                         (0, f"onpage: {self.page(ARMS)} holds the left-out table as print gives it\n"))

    def test_floor_with_no_dataset_left_out_passes_with_its_one_line(self):
        self.assertEqual(self.to_floor("all in")["left-out"],
                         ["0 of the 6 datasets of M are left out of the second reading."])
        self.assertEqual(check_tree(self.root, self.lists), "floor")

    def test_floor_refuses_a_page_with_the_path_alone(self):
        self.assertEqual(self.cli("floor", BUILT["seven"], self.record)[0], 0)
        self.off(f"^onpage: {re.escape(str(self.page(ARMS)))} lacks the left-out table's header row$")
        self.assertEqual(self.cli("onpage", self.record, self.page(ARMS)),
                         (1, f"onpage: {self.page(ARMS)} lacks the left-out table's header row\n"))

    def test_floor_refuses_a_page_with_a_row_changed(self):
        self.to_floor()
        self.change(ARMS, "| left out, under the floor | 30 R, 40 S |", "| left out, under the floor | 31 R, 40 S |")
        self.off("lacks the left-out table's row of Synth d with drug 10$")

    def test_floor_refuses_a_page_that_does_not_name_the_file(self):
        self.to_floor()
        self.change(ARMS, f"It goes on the page as {LEFT_PATH} before", "It goes on the page before")
        self.off(f"^{ARMS} does not name {re.escape(LEFT_PATH)}$")

    def test_floor_refuses_a_count_table_on_any_page(self):
        tables = self.to_floor()
        self.change(ARMS, tables["left-out"][-1], tables["left-out"][-1] + "\n\n" + m_page.count_headers()[0])
        self.off(f"^{ARMS} holds a count table before the count$")

    def test_floor_refuses_a_record_that_breaks_a_gate(self):
        self.to_floor()
        (self.record / "m-left-out.tsv").write_bytes((self.record / "m-left-out.tsv").read_bytes() + b"\n")
        self.off("^m_page: stopped: m-left-out.tsv differs from its SHA256SUMS line$")


class CountTaken(Copy):
    def test_count_passes_with_the_blocks_of_page_on_the_result_page(self):
        blocks = self.to_count()
        self.assertEqual(m_page.state(self.record), "count")
        self.assertEqual(check_tree(self.root, self.lists), "count")
        self.assertEqual(tuple(blocks), ("reads", "both", "species", "count-all", "recall-all", "study-all", "count-m",
                                         "files"))
        self.assertEqual(len(self.printed()), 8)

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

    def test_count_refuses_a_result_page_without_one_block(self):
        for name in m_page.ON_PAGE:
            with self.subTest(name=name):
                self.fresh()
                blocks = self.to_count()
                self.change(RESULT, "\n".join(blocks[name]), "")
                self.off(f"^{RESULT} lacks the block {name} as page gives it$")

    def test_count_refuses_a_block_that_runs_on_past_its_last_row_or_has_no_blank_line_before_it(self):
        # The first run's rows in the page's form under the block on all datasets; the block on M alone's rows again
        # under it; a line of words under the files line; the reads line run into the block after it.
        for name, change in (
                ("count-all", lambda b, f: (b["count-all"][-1], b["count-all"][-1] + "\n" + "\n".join(f[2:]))),
                ("count-m", lambda b, f: (b["count-m"][-1], b["count-m"][-1] + "\n" + "\n".join(b["count-m"][4:]))),
                ("files", lambda b, f: (b["files"][0], b["files"][0] + "\nMore words.")),
                ("species", lambda b, f: (b["both"][-1] + "\n\n" + b["species"][0],
                                          b["both"][-1] + "\n\nWords.\n" + b["species"][0]))):
            with self.subTest(name=name):
                self.fresh()
                blocks = self.to_count()
                first = m_page.print_count(m_page.read_tables(self.record, ["m-count.tsv"])["m-count.tsv"], "first run")
                self.change(RESULT, *change(blocks, first))
                side = "before" if name == "species" else "after"
                if name == "count-all":
                    # The first run's rows are rows of the record's count table in the page's form.
                    self.off(f"^{RESULT} holds a row of a table that stays in the record$")
                    self.change(RESULT, "\n".join(first[2:]), "")
                    self.change(RESULT, blocks["count-all"][-1] + "\n\n", blocks["count-all"][-1] + "\nSo it reads.\n\n")
                self.off(f"^{RESULT} holds the block {name} with no blank line {side} it$")

    def test_count_refuses_the_blocks_twice(self):
        blocks = self.to_count()
        path = self.page(RESULT)
        path.write_text(path.read_text(encoding="utf-8") + "\n" + "\n\n".join("\n".join(b) for b in blocks.values())
                        + "\n", encoding="utf-8")
        self.off(f"^{RESULT} holds the first line of the block reads more than once$")

    def test_count_refuses_a_row_of_a_table_that_stays_in_the_record_on_any_page(self):
        # Each table's rows with no header, on either page; and the first run's count in the page's own form, as a
        # block of its own under the header the page's count blocks share.
        for name in ("count", "beside", "conditions", "size", "datasets"):
            for page in (RESULT, ARMS):
                with self.subTest(name=name, page=page):
                    self.fresh()
                    self.to_count()
                    path = self.page(page)
                    path.write_text(path.read_text(encoding="utf-8") + "\n" + "\n".join(self.printed()[name][4:])
                                    + "\n", encoding="utf-8")
                    self.off(f"^{page} holds a row of a table that stays in the record$")
        self.fresh()
        self.to_count()
        count = m_page.read_tables(self.record, ["m-count.tsv"])["m-count.tsv"]
        for c, split, measure in (("first run", "lineage", "balanced_accuracy"), ("M", "lineage", "recall"),
                                  ("first run", "study", "balanced_accuracy")):
            with self.subTest(count=c, split=split, measure=measure):
                lines = (m_page.print_count(count, c) if (split, measure) == ("lineage", "balanced_accuracy")
                         else m_page.print_other(count, c, split, measure))
                path = self.page(RESULT)
                kept = path.read_text(encoding="utf-8")
                path.write_text(kept + "\n**More.** Another count.\n\n" + "\n".join(lines) + "\n", encoding="utf-8")
                self.off(f"^{RESULT} holds a row of a table that stays in the record$")
                path.write_text(kept, encoding="utf-8")

    def test_count_refuses_a_result_page_with_a_row_changed(self):
        blocks = self.to_count()
        row = blocks["count-m"][-1]
        cells = row.split(" | ")
        cells[6] = "0.999" if cells[6] != "0.999" else "0.998"
        self.change(RESULT, row, " | ".join(cells))
        self.off(f"^{RESULT} lacks the block count-m as page gives it$")

    def test_count_refuses_the_blocks_in_another_order(self):
        blocks = self.to_count()
        first, second = ("\n".join(blocks[name]) for name in ("count-all", "recall-all"))
        self.change(RESULT, first + "\n\n" + second, second + "\n\n" + first)
        self.off(f"^{RESULT} holds the block recall-all before the one page gives ahead of it$")

    def test_count_refuses_a_table_that_stays_in_the_record_on_any_page(self):
        for name in ("count", "beside", "conditions", "size", "datasets"):
            for page in (RESULT, ARMS):
                with self.subTest(name=name, page=page):
                    self.fresh()
                    self.to_count()
                    path = self.page(page)
                    path.write_text(path.read_text(encoding="utf-8") + "\n" + "\n".join(self.printed()[name]) + "\n",
                                    encoding="utf-8")
                    self.off(f"^{page} holds a table that stays in the record$")

    def test_count_refuses_a_left_out_table_changed_on_the_design_page(self):
        self.to_count()
        self.change(ARMS, "| left out, under the floor | 30 R, 40 S |", "| left out, under the floor | 31 R, 40 S |")
        self.off("lacks the left-out table's row of Synth d with drug 10$")

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

    def test_a_record_file_unlike_its_line_fails_the_committed_tree_test_and_not_the_setup(self):
        # A tree whose record breaks a gate still gives the lines the copies take off: those of the file gate turned
        # off, or none. The test of that tree names the file.
        # Each change leaves what print and page give as it was: a species' slug, which no table prints; a study-split
        # figure, which the datasets table does not print; a line more in the times file, which m_page.py only sums.
        def study(text):
            lines = text.split("\n")
            n = next(i for i, l in enumerate(lines) if l.split("\t")[4:5] == ["study"])
            cells = lines[n].split("\t")
            at = lines[0].split("\t").index("arm_value")
            cells[at] = "0.5" if cells[at] != "0.5" else "0.25"
            lines[n] = "\t".join(cells)
            return "\n".join(lines)

        for upto, name, change in (("floor", "m-floor.tsv", lambda t: t.replace("\tsynth_c\t", "\tsynth_x\t", 1)),
                                   ("floor", "first-run-times.tsv", lambda t: t + "\n"),
                                   ("count", "m-datasets.tsv", study),
                                   ("count", "first-run-times.tsv", lambda t: t + "\n")):
            with self.subTest(upto=upto, name=name):
                self.fresh()
                (self.to_floor if upto == "floor" else self.to_count)()
                want = real_lines(self.root)
                self.assertIn(m_page.BOTH_LINE if upto == "count" else m_page.table(m_page.SPECIES_HEADER, [])[0],
                              want)
                m_tables.edit(self.record, name, change, sums=False)
                self.assertEqual(real_lines(self.root), want)
                with self.assertRaisesRegex(Off, f"^m_page: stopped: {name} differs from its SHA256SUMS line$"):
                    check_tree(self.root)
        (self.record / "m-size.tsv").write_text("x\n")
        self.assertEqual(real_lines(self.root), set())

    def test_count_checks_each_sum_itself_as_well(self):
        # The record's own file check is turned off here, so the test's sum check is the one that stops; the times file
        # is one m_page.py sums and reads no further.
        self.to_count()
        path = self.record / "first-run-times.tsv"
        path.write_bytes(path.read_bytes() + b"\n")
        gate, m_page.gate_files = m_page.gate_files, lambda record: None
        try:
            self.off("^first-run-times.tsv differs from its SHA256SUMS line$")
        finally:
            m_page.gate_files = gate


if __name__ == "__main__":
    unittest.main()
