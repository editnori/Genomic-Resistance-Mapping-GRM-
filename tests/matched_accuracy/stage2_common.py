"""Shared runs of stage 2 for tests/matched_accuracy: a SYNTHETIC release of made_up2.py, stage 1 on it, then stage 2.

Every release here is made up (made_up2.py): species, antibiotics, genome ids (gz-<species>-<n>, not the shape of a
real genome id), BioProjects, lab records, tables, catalogue calls and scores. The pins of both scripts are set to
the made-up files' by the module constants, a route the command line does not reach (test_stage2_reference.py).
"""

import contextlib
import csv
import hashlib
import io
import json
import pathlib
import tempfile

from tests.matched_accuracy import made_up2

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "matched-accuracy" / "stage2.py"
FAULTS = ROOT / "scripts" / "matched-accuracy" / "stage2_faults.py"
STAGE1_COPY = ROOT / "reference" / "matched-accuracy" / "stage1-rev14-274-20261112"
stage1, stage2 = made_up2.load_scripts()
STOP = "stage2: stopped: "
EXPORT = ["cells.tsv", "inputs.tsv", "interval.tsv", "measure.tsv", "null_fits.tsv", "typing.tsv"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


class Run:
    """One run of stage 1, then stage2.main, on a release written in a new temporary directory.

    `before(release)` changes the release before stage 1 runs; `after(run)` changes the written files after it,
    before the pins of stage 2 are read from them; `pins` (a dict, or a function of the pins read) and `inner`
    then override stage 2's own pins and those of the stage 1 module it imports; `argv` maps the argument list to
    another. A Python error inside stage 2
    becomes code "error" and a line naming it, so a test meets it as an assertion."""

    def __init__(self, plan=None, before=None, after=None, pins=None, inner=None, argv=None):
        self._tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(self._tmp.name)
        self.top, self.stage1_out, self.out = base / "in", base / "stage1", base / "stage2"
        self.release = made_up2.build(made_up2.designed() if plan is None else plan)
        if before:
            before(self.release)
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            self.stage1_code, _ = made_up2.stage1_run(stage1, self.top, self.release, self.stage1_out)
        self.stage1_lines = printed.getvalue().splitlines()
        if self.stage1_code != 0:
            raise AssertionError(f"stage 1 did not pass on the made-up release: {self.stage1_lines}")
        if after:
            after(self)
        own, inner_pins = made_up2.stage2_pins(self.top, self.stage1_out)
        own.update(pins(own) if callable(pins) else pins or {})
        inner_pins.update(inner or {})
        self.argv = made_up2.stage2_argv(self.top, self.stage1_out, self.out)
        if argv:
            self.argv = argv(self.argv)
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            try:
                self.code = made_up2.stage2_run(stage2, self.argv, own, inner_pins)
            except Exception as e:  # noqa: BLE001 - a Python error is a result the test compares
                self.code = "error"
                print(f"python error: {type(e).__name__}: {e}")
        self.lines = printed.getvalue().splitlines()

    def close(self):
        self._tmp.cleanup()

    # The files of the run.
    def row_list(self, reading):
        return self.stage1_out / "local" / f"reading-{reading}.tsv"

    def report(self, dataset):
        return self.top / "reports" / f"{dataset}.json"

    def counts(self):
        return self.stage1_out / "export" / "counts.tsv"

    def typing_table(self):
        return self.stage1_out / "export" / "typing.tsv"

    def rewrite_table(self, path, change):
        """Rewrite one of stage 1's tables: change(head, fields) changes a data line's fields in place."""
        lines = path.read_text(encoding="utf-8").split("\n")
        head = lines[0].split("\t")
        for i, line in enumerate(lines[1:-1], 1):
            f = line.split("\t")
            change(head, f)
            lines[i] = "\t".join(f)
        path.write_text("\n".join(lines), encoding="utf-8")

    def rewrite_list(self, reading, change):
        """Rewrite a row list: change(lines) gets and returns its lines, each with its line break."""
        path = self.row_list(reading)
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        path.write_text("".join(change(lines)), encoding="utf-8")

    def rewrite_report(self, dataset, change, sums=True):
        """Rewrite one report: change(doc) changes it in place. With sums, its line of the made-up sums list too."""
        path = self.report(dataset)
        doc = json.loads(path.read_bytes())
        change(doc)
        path.write_bytes(json.dumps(doc, sort_keys=True).encode("utf-8"))
        if sums:
            self.rewrite_sums()

    def rewrite_sums(self):
        """Write the made-up sums list again from the reports' files as they now are."""
        lines = []
        for p in sorted((self.top / "reports").glob("*.json")):
            lines.append(f"{sha(p.read_bytes())}  ./datasets/{p.name}\n")
        (self.top / "reports.sha256").write_text("".join(lines), encoding="utf-8")

    def wrote_export(self):
        return (self.out / "export").exists()

    def wrote(self):
        return self.out.exists() and any(self.out.iterdir())

    def export(self, name):
        return (self.out / "export" / name).read_bytes()

    def table(self, name):
        text = self.export(name).decode("utf-8")
        if name == "inputs.tsv":
            text = text.rsplit(stage2.OPENED + "\n", 1)[0]
        return list(csv.DictReader(io.StringIO(text), delimiter="\t"))

    def measure(self, reading, level, name):
        return next(r for r in self.table("measure.tsv")
                    if (r["reading"], r["level"], r["name"]) == (reading, level, name))


class Case:
    """A mixin: run() closes each run when the test ends."""

    def run_stage2(self, plan=None, **kw):
        r = Run(plan, **kw)
        self.addCleanup(r.close)
        return r

    def assert_stops(self, run, words):
        """Before the line: exit 1, one printed line naming the check, and nothing written."""
        self.assertEqual(run.code, 1, run.lines)
        self.assertEqual(len(run.lines), 1, run.lines)
        self.assertTrue(run.lines[0].startswith(STOP), run.lines)
        self.assertIn(words, run.lines[0])
        self.assertFalse(run.wrote(), run.lines)

    def assert_stops_after_the_line(self, run, words):
        """After the line: exit 1, the line and one line naming the check, and no export."""
        self.assertEqual(run.code, 1, run.lines)
        self.assertEqual(run.lines[0], stage2.OPENED, run.lines)
        self.assertEqual(len(run.lines), 2, run.lines)
        self.assertTrue(run.lines[1].startswith(STOP), run.lines)
        self.assertIn(words, run.lines[1])
        self.assertFalse(run.wrote_export(), run.lines)

    def assert_passes(self, run):
        self.assertEqual(run.code, 0, run.lines)
