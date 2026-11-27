"""Shared paths, the script and SYNTHETIC releases for tests/matched_accuracy.

Every release here is made up by tests/matched_accuracy/made_up.py: made-up species, antibiotics, genome ids
(gz-<species>-<n>, not the shape of a real genome id), BioProject accessions (pz-<species>-<n>), lab records,
label table, genome tables and scores. A test that needs the shape of a real genome id says so.
"""

import contextlib
import copy
import csv
import io
import os
import pathlib
import tempfile

from tests.matched_accuracy import made_up

ROOT = pathlib.Path(__file__).resolve().parents[2]
SUITE = "tests/matched_accuracy"
SCRIPT = ROOT / "scripts" / "matched-accuracy" / "stage1.py"
FAULTS = ROOT / "scripts" / "matched-accuracy" / "planted_faults.py"
PAGE = ROOT / "docs" / "matched-accuracy.md"
INDEX = ROOT / "docs" / "INDEX.md"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
REFERENCE = ROOT / "reference" / "matched-accuracy"
REPORT_SUMS = ROOT / "reference" / "lineage-control" / "rev14-274-20261112.reports.sha256"
stage1 = made_up.load_script()
STOP = "stage1: stopped: "


def release():
    """A fresh copy of the made-up release, to change in one test."""
    return copy.deepcopy(made_up.made_up())


def report_of(rel, dataset):
    return next(r for r in rel["reports"] if r["id"] == dataset)


def row_of(rel, dataset, genome):
    return next(s for s in report_of(rel, dataset)["samples"] if s["id"] == genome)


def records_of(rel, genome, raw):
    return [r for r in rel["records"] if r["genome_id"] == genome and r["antibiotic"].strip() == raw]


def report_bytes(dataset, change):
    """An `after` that rewrites one written report's bytes and its line of the made-up sums list."""
    def after(top):
        path = top / "reports" / f"{dataset}.json"
        path.write_bytes(change(path.read_bytes()))
        lines = []
        for p in sorted((top / "reports").glob("*.json")):
            lines.append(f"{made_up.sha(p.read_bytes())}  ./datasets/{p.name}\n")
        sums = "".join(lines).encode("utf-8")
        (top / "reports.sha256").write_bytes(sums)
        return {"REPORT_SUMS_SHA256": made_up.sha(sums)}
    return after


def outcome(fn, *args, **kw):
    """("ok", value) or ("stopped", reason): a stop becomes a value a test can compare, never an error."""
    try:
        return ("ok", fn(*args, **kw))
    except stage1.Stop as e:
        return ("stopped", str(e))


class Run:
    """One run of stage1.main on a release written in a new temporary directory.

    `after(top)` changes the written files and may return pins to set. A Python error inside main becomes
    code "error" and a line naming it, so a test meets it as an assertion, never as an error of its own."""

    def __init__(self, rel=None, hc20=True, after=None, pins=None, extra=None, umask=None):
        self._tmp = tempfile.TemporaryDirectory()
        self.top = pathlib.Path(self._tmp.name) / "in"
        self.out = pathlib.Path(self._tmp.name) / "out"
        argv, self.pins = made_up.write_release(self.top, rel or release(), hc20)
        if after:
            self.pins.update(after(self.top) or {})
        self.pins.update(pins or {})
        printed = io.StringIO()
        old_mask = None if umask is None else os.umask(umask)
        with contextlib.redirect_stdout(printed):
            try:
                self.code = made_up.run(stage1, argv + ["--out", str(self.out)] + (extra or []), self.pins)
            except Exception as e:  # noqa: BLE001 - a Python error is a result the test compares
                self.code = "error"
                print(f"python error: {type(e).__name__}")
            finally:
                if old_mask is not None:
                    os.umask(old_mask)
        self.lines = printed.getvalue().splitlines()

    def close(self):
        self._tmp.cleanup()

    def wrote(self):
        """True when the run wrote anything under its output directory."""
        return self.out.exists() and any(self.out.iterdir())

    def export(self, name):
        return (self.out / "export" / name).read_bytes()

    def table(self, name):
        text = self.export(name).decode("utf-8")
        if name == "inputs.tsv":
            text = text.rsplit("scores opened: no\n", 1)[0]
        return list(csv.DictReader(io.StringIO(text), delimiter="\t"))

    def counts(self, level, name, reading):
        return next(r for r in self.table("counts.tsv")
                    if (r["level"], r["name"], r["reading"]) == (level, name, reading))

    def typing(self, level, name, reading, step):
        return next(r for r in self.table("typing.tsv")
                    if (r["level"], r["name"], r["reading"], r["step"]) == (level, name, reading, step))

    def row_list(self, reading):
        return (self.out / "local" / f"reading-{reading}.tsv").read_bytes()


class Case:
    """A mixin: run() closes each run when the test ends."""

    def run_stage1(self, rel=None, **kw):
        r = Run(rel, **kw)
        self.addCleanup(r.close)
        return r

    def assert_stops(self, run, words):
        """Exit 1, one printed line naming the check, and nothing written."""
        self.assertEqual(run.code, 1, run.lines)
        self.assertEqual(len(run.lines), 1, run.lines)
        self.assertTrue(run.lines[0].startswith(STOP), run.lines)
        self.assertIn(words, run.lines[0])
        self.assertFalse(run.wrote())

    def assert_passes(self, run):
        self.assertEqual(run.code, 0, run.lines)
