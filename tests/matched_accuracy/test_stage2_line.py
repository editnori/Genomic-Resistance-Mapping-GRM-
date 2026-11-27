"""The line "scores opened: yes" of scripts/matched-accuracy/stage2.py: every check runs before it, no score can be
reached before it, and the second loader returns the named values of listed rows alone.

Every report here is SYNTHETIC (made_up2.py): its scores are made-up numbers, and its roc, metrics, conformal,
no_call, called_metrics, lasso_fits and top-level catalogue carry made-up marker texts."""

import ast
import contextlib
import io
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.matched_accuracy import made_up2, stage2_common as c
from tests.matched_accuracy.stage2_common import stage2


def reachable(obj, seen=None):
    """Every object reachable from obj through dict keys and values, list, tuple and set items and attributes."""
    seen = set() if seen is None else seen
    if id(obj) in seen:
        return
    seen.add(id(obj))
    yield obj
    items = []
    if isinstance(obj, dict):
        items = list(obj.keys()) + list(obj.values())
    elif isinstance(obj, (list, tuple, set, frozenset)):
        items = list(obj)
    elif hasattr(obj, "__dict__") and not isinstance(obj, type):
        items = list(vars(obj).values())
    for x in items:
        yield from reachable(x, seen)


def tree():
    return ast.parse(c.SCRIPT.read_text(encoding="utf-8"))


def functions():
    return {f.name: f for f in tree().body if isinstance(f, ast.FunctionDef)}


def called(node):
    """The names of the plain calls under node, in source order."""
    calls = [n for n in ast.walk(node) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)]
    return [n.func.id for n in sorted(calls, key=lambda n: (n.lineno, n.col_offset))]


class Order(c.Case, unittest.TestCase):
    """A run records each check, the line and each read of scores, in the order they happen."""

    def test_every_check_runs_before_the_line_and_every_score_after_it(self):
        events = []

        def recorder(name, real):
            def wrapped(*a, **kw):
                events.append(name)
                return real(*a, **kw)
            return wrapped

        patches = [mock.patch.object(stage2, n, recorder(n, getattr(stage2, n)))
                   for n in stage2.CHECKS + ("count_bad", "open_scores", "read_scores", "pick")]
        with contextlib.ExitStack() as stack:
            for p in patches:
                stack.enter_context(p)
            run = self.run_stage2()
        self.assert_passes(run)
        self.assertEqual(run.lines[0], stage2.OPENED)
        line = events.index("open_scores")
        self.assertEqual(events.count("open_scores"), 1)
        for n in stage2.CHECKS + ("count_bad",):
            self.assertIn(n, events, n)
            self.assertLess(max(i for i, e in enumerate(events) if e == n), line, n)
        # count_bad reads each of the 10 listed reports once, before the line.
        self.assertEqual(events.count("count_bad"), 10)
        reads = [i for i, e in enumerate(events) if e in ("read_scores", "pick")]
        self.assertEqual(events.count("read_scores"), 10)
        self.assertGreater(min(reads), line)

    def test_after_the_line_only_the_bytes_the_scores_the_export_scan_and_a_write_can_stop_the_run(self):
        # A profile of one whole run: every function of either script called after the line that holds a stop.
        stops = {}
        for script in (c.SCRIPT, made_up2.STAGE1):
            for f in ast.walk(ast.parse(script.read_text(encoding="utf-8"))):
                if isinstance(f, ast.FunctionDef) and any(
                        isinstance(n, ast.Raise) and isinstance(n.exc, ast.Call) and getattr(n.exc.func, "id", "")
                        == "Stop" for n in ast.walk(f)):
                    stops[(script.name, f.name)] = f
        calls = []

        def profile(frame, event, arg):
            if event == "call" and Path(frame.f_code.co_filename).name in ("stage1.py", "stage2.py"):
                calls.append((Path(frame.f_code.co_filename).name, frame.f_code.co_name))
        old = sys.getprofile()
        sys.setprofile(profile)
        try:
            run = self.run_stage2()
        finally:
            sys.setprofile(old)
        self.assert_passes(run)
        line = calls.index(("stage2.py", "open_scores"))
        before = {x for x in calls[:line] if x in stops}
        after = {x for x in calls[line + 1:] if x in stops}
        self.assertLessEqual({("stage2.py", n) for n in stage2.CHECKS + ("count_bad",)} & set(stops), before)
        # The re-check of a report's bytes, the score stops, the export's scan and a write that fails.
        self.assertEqual(after, {("stage2.py", "read_scores"), ("stage2.py", "component"),
                                 ("stage2.py", "scan_export"), ("stage1.py", "scan"), ("stage2.py", "write_export")})

    def test_only_check_calls_the_checks_and_run_calls_check_before_the_line_and_measure_after_it(self):
        fs = functions()
        for n in stage2.CHECKS:
            callers = sorted(name for name, f in fs.items() if n in called(f) and name != n)
            self.assertEqual(callers, ["check"], n)
        self.assertEqual(sorted(name for name, f in fs.items() if "count_bad" in called(f)), ["check_values"])
        body = called(fs["run"])
        self.assertEqual(body[:3], ["check", "open_scores", "measure"])
        for n in ("read_scores", "pick"):
            callers = sorted(name for name, f in fs.items() if n in called(f) and name != n)
            self.assertEqual(callers, {"read_scores": ["measure"], "pick": ["read_scores"]}[n], n)
        self.assertEqual(sorted(name for name, f in fs.items() if "measure" in called(f)), ["run"])
        self.assertEqual(sorted(name for name, f in fs.items() if "open_scores" in called(f)), ["run"])
        # check() reaches no score: it calls no function of the second loader or of the measure.
        self.assertEqual(set(called(fs["check"])) & {"read_scores", "pick", "open_scores", "measure", "row_mean",
                                                     "component"}, set())

    def test_what_check_returns_holds_no_score_and_no_number_but_integers(self):
        run = self.run_stage2(made_up2.randomised(5))
        self.assert_passes(run)
        scores = set()
        for doc in run.release["reports"]:
            for sp in doc["splits"]:
                for e in sp["learners"]:
                    scores |= {x for x in e["scores"] if isinstance(x, float)}
        self.assertGreater(len(scores), 100)
        args = stage2.parse(run.argv[:-1] + [run.argv[-1] + "-again"])
        own, inner = made_up2.stage2_pins(run.top, run.stage1_out)
        old_own, old_inner = made_up2.set_pins(stage2, own), made_up2.set_pins(stage2.stage1, inner)
        try:
            checked = stage2.check(args)
        finally:
            made_up2.set_pins(stage2, old_own)
            made_up2.set_pins(stage2.stage1, old_inner)
        floats = [x for x in reachable(checked) if isinstance(x, float)]
        self.assertEqual(floats, [])


class LineWrite(c.Case, unittest.TestCase):
    def test_a_line_file_that_cannot_be_written_stops_with_the_stop_line_alone(self):
        real = os.open

        def refused(path, *a, **kw):
            if Path(path).name == stage2.OPENED_FILE:
                raise PermissionError(13, "made-up refusal")
            return real(path, *a, **kw)
        with mock.patch.object(stage2.os, "open", side_effect=refused):
            run = self.run_stage2()
        self.assertEqual(run.code, 1, run.lines)
        self.assertEqual(run.lines, [c.STOP + "the line could not be written (PermissionError)"])
        self.assertFalse((run.out / stage2.OPENED_FILE).exists())
        self.assertFalse(run.wrote_export())


class Gate(unittest.TestCase):
    """The second loader takes the line's token alone, and its file must hold the line."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = Path(self._tmp.name)
        doc = made_up2.build(made_up2.designed())["reports"][0]
        self.path = self.dir / "report.json"
        self.data = json.dumps(doc).encode("utf-8")
        self.path.write_bytes(self.data)
        self.listed = [(doc["samples"][0]["id"], 0)]

    def outcome(self, opened):
        try:
            return ("ok", stage2.read_scores(self.path, c.sha(self.data), self.listed, opened))
        except stage2.Stop as e:
            return ("stopped", str(e))

    def test_no_score_without_the_line(self):
        self.assertEqual(self.outcome(None), ("stopped", "a score was asked for before the line"))
        self.assertEqual(self.outcome(stage2.Opened(self.dir / stage2.OPENED_FILE)),
                         ("stopped", "a score was asked for before the line"))
        (self.dir / stage2.OPENED_FILE).write_text("scores opened: no\n")
        self.assertEqual(self.outcome(stage2.Opened(self.dir / stage2.OPENED_FILE)),
                         ("stopped", "a score was asked for before the line"))

    def test_the_line_is_written_then_printed_then_scores_open(self):
        printed = io.StringIO()
        with contextlib.redirect_stdout(printed):
            opened = stage2.open_scores(self.dir / "out")
        self.assertEqual(printed.getvalue(), stage2.OPENED + "\n")
        self.assertEqual((self.dir / "out" / stage2.OPENED_FILE).read_text(), stage2.OPENED + "\n")
        self.assertEqual(self.outcome(opened)[0], "ok")

    def test_a_report_changed_after_its_check_stops(self):
        with contextlib.redirect_stdout(io.StringIO()):
            opened = stage2.open_scores(self.dir / "out")
        self.path.write_bytes(self.data + b" ")
        self.assertEqual(self.outcome(opened), ("stopped", "report.json changed after it was checked"))


class SecondLoader(unittest.TestCase):
    """What read_scores returns: listed rows' five scores and calls, five fold_complexity lists, the catalogue as
    present or None, and nothing else."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        d = Path(self._tmp.name)
        plan = made_up2.randomised(11, species=1)
        self.doc = made_up2.build(plan)["reports"][0]
        # Every score of the report a number of its own: learner j, row i.
        for sp in self.doc["splits"]:
            for j, e in enumerate(sp["learners"]):
                e["scores"] = [0.5 + j / 100 + i / 100000 + (0.25 if sp["scheme"] == "random" else 0)
                               for i in range(len(e["scores"]))]
        data = json.dumps(self.doc).encode("utf-8")
        (d / "r.json").write_bytes(data)
        n = len(self.doc["samples"])
        self.listed = [(self.doc["samples"][i]["id"], i) for i in range(0, n, 2)]
        with contextlib.redirect_stdout(io.StringIO()):
            opened = stage2.open_scores(d / "out")
        self.got = stage2.read_scores(d / "r.json", c.sha(data), self.listed, opened)

    def test_the_named_values_alone(self):
        self.assertEqual(set(self.got), {"rows", "fold_complexity", "catalog"})
        self.assertEqual(set(self.got["rows"]), {g for g, _ in self.listed})
        self.assertEqual(set(self.got["fold_complexity"]), set(stage2.FIVE))
        self.assertEqual(self.got["catalog"], "present")
        lineage = next(s for s in self.doc["splits"] if s["scheme"] == "lineage")
        by = {e["learner"]: e for e in lineage["learners"]}
        for g, i in self.listed:
            five, call = self.got["rows"][g]
            self.assertEqual(five, tuple(by[n]["scores"][i] for n in stage2.FIVE))
            self.assertEqual(call, self.doc["samples"][i]["catalog_resistant"])

    def test_no_score_of_an_unlisted_row_no_catalogue_score_and_no_value_of_the_named_keys_is_reachable(self):
        listed = {i for _, i in self.listed}
        lineage = next(s for s in self.doc["splits"] if s["scheme"] == "lineage")
        barred = set()
        for sp in self.doc["splits"]:
            for e in sp["learners"]:
                barred |= {x for i, x in enumerate(e["scores"])
                           if not (sp is lineage and e["learner"] in stage2.FIVE and i in listed)}
        self.assertTrue(any(e["learner"] == "catalog" for e in lineage["learners"]))
        texts = set(made_up2.SENTINELS.values())
        values = list(reachable(self.got))
        self.assertEqual([x for x in values if isinstance(x, float) and x in barred], [])
        self.assertEqual([x for x in values if isinstance(x, str) and x in texts], [])
        self.assertNotIn(made_up2.BALANCED, [x for x in values if isinstance(x, float)])
        for key in ("roc", "metrics", "conformal", "no_call", "called_metrics", "lasso_fits",
                    "fold_balanced_accuracy", "groups", "rule", "version"):
            self.assertNotIn(key, [x for x in values if isinstance(x, str)], key)

    def test_a_report_without_a_drug_rule_gives_none(self):
        doc = json.loads(json.dumps(self.doc))
        doc["catalog"] = None
        with tempfile.TemporaryDirectory() as d:
            data = json.dumps(doc).encode("utf-8")
            (Path(d) / "r.json").write_bytes(data)
            with contextlib.redirect_stdout(io.StringIO()):
                opened = stage2.open_scores(Path(d) / "out")
            got = stage2.read_scores(Path(d) / "r.json", c.sha(data), self.listed, opened)
        self.assertIsNone(got["catalog"])


if __name__ == "__main__":
    unittest.main()
