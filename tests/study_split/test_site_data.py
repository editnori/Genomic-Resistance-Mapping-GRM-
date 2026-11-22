"""The study-split numbers the web app bundles (web/src/data/study-split.json)."""

import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

from tests.study_split.common import RECORD, ROOT

BUNDLED = ROOT / "web" / "src" / "data" / "study-split.json"


def load_site_data():
    spec = importlib.util.spec_from_file_location("study_split_site_data", ROOT / "scripts" / "study-split" / "site_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_site_data()


class SiteData(unittest.TestCase):
    def test_bundled_file_is_what_the_record_gives(self):
        self.assertEqual(BUNDLED.read_text(encoding="utf-8"), S.build(RECORD))
        self.assertEqual(S.main(["site_data.py", str(RECORD), "--check", str(BUNDLED)]), 0)

    def test_check_fails_when_a_record_file_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            copy = pathlib.Path(tmp) / "record"
            shutil.copytree(RECORD, copy)
            p = copy / "summary.tsv"
            p.write_text(p.read_text(encoding="utf-8").replace("datasets\t70", "datasets\t71"), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(copy), "--check", str(BUNDLED)]), 1)

    def test_bundled_counts_match_headline_record_values(self):
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        self.assertEqual((data["datasets"], data["species"], data["pairs"]), (70, 9, 57914))
        self.assertEqual(sum(x["datasets"] for x in data["species_counts"]), 70)
        for learner in S.LEARNERS:
            self.assertEqual(data["first"][learner]["rule"], "cost")
            self.assertEqual(data["second"][learner]["SR"]["cost"], 25)
            self.assertEqual(data["second"][learner]["SM"]["cost"], 5)
        raw = BUNDLED.read_text(encoding="utf-8")
        self.assertNotIn("SAMN", raw)
        self.assertNotIn("GCA_", raw)
