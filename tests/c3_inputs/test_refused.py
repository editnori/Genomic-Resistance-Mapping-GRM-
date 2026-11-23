"""scripts/c3-inputs/refused.py: the record of an arm's study fit refused for no candidate (ruling 16).

Every log, dataset name and count is SYNTHETIC. The refusal's sentence is filled in from grm's own format
string in src/learn/eval.rs, so a change to grm's words fails here before a run.
"""

import contextlib
import io
import json
import re
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common

refused = common.load("refused")
score = common.load("score")


def grm_sentence(scheme="study-whole", fold=2, train=87, floor=3):
    """grm's refusal sentence for a fold with no candidate, from its format string, filled in."""
    src = (common.ROOT / "src" / "learn" / "eval.rs").read_text()
    found = re.findall(r'"(\{scheme\} outer fold \{f\}\'s \{\} training samples give no candidate feature[^"]*)"', src)
    if len(found) != 1:
        raise AssertionError(f"grm's sentence for a fold with no candidate: {len(found)} found in eval.rs")
    text = found[0]
    for key, value in (("{scheme}", scheme), ("{f}", fold), ("{min_support}", floor), ("{}", train)):
        text = text.replace(key, str(value))
    return text


def log(sentence):
    return f"SYNTHETIC log line\nError: {common.DATASET}\n\nCaused by:\n    {sentence}\n"


class Record(unittest.TestCase):
    def test_grms_own_sentence_gives_the_record_with_its_counts(self):
        rec = refused.record(log(grm_sentence()), common.DATASET, "program")
        self.assertEqual(rec, {"format": "grm-c3-refused-1", "dataset": common.DATASET, "arm": "program",
                               "scheme": "study-whole", "seed": 42, "reason": "no candidate", "fold": 2,
                               "training_genomes": 87, "support_floor": 3, "detail": grm_sentence()})

    def test_the_default_and_the_twin_get_no_record(self):
        for arm in ("k31", "twin"):
            with self.subTest(arm=arm):
                self.assertIsNone(refused.record(log(grm_sentence()), common.DATASET, arm))

    def test_another_refusal_or_scheme_gets_no_record(self):
        for text in ("study-whole outer fold 0's training samples lack one or both classes",
                     grm_sentence(scheme="lineage"), grm_sentence(floor=3).replace("at least 3 of", "at least 4 of"), ""):
            with self.subTest(text=text[:40]):
                self.assertIsNone(refused.record(log(text), common.DATASET, "program"))

    def test_the_record_written_is_the_one_score_reads(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fit.log"
            path.write_text(log(grm_sentence()))
            out = Path(tmp) / f"{common.DATASET}.study-whole-s42.refused.json"
            with contextlib.redirect_stdout(io.StringIO()):
                rc = refused.main([str(path), common.DATASET, "program", str(out)])
            self.assertEqual(rc, 0)
            self.assertEqual(sorted(p.name for p in Path(tmp).iterdir()), ["fit.log", out.name])
            rec = json.loads(out.read_text())
        exp = {"dataset": common.DATASET, "arm": "program", "split": "study", "scheme": "study-whole", "seed": "42"}
        base = common.pick(scheme="study-whole")
        checks = score.check_refused(exp, rec, base)
        self.assertEqual([c for c in checks if not c[1]], [])
        self.assertEqual(len(checks), 5)
        # Counts and grm's sentence only: no path.
        self.assertFalse([v for v in rec.values() if "/" in str(v)])

    def test_nothing_is_written_for_another_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fit.log"
            path.write_text(log("study-whole outer fold 1's test samples lack one or both classes"))
            out = Path(tmp) / "out.refused.json"
            rc = refused.main([str(path), common.DATASET, "program", str(out)])
            self.assertEqual(rc, 1)
            self.assertEqual([p.name for p in Path(tmp).iterdir()], ["fit.log"])


if __name__ == "__main__":
    unittest.main()
