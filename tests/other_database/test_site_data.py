"""The other-database numbers the web app bundles (web/src/data/other-database.json)."""

import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

from tests.other_database.common import RECORD, ROOT

BUNDLED = ROOT / "web" / "src" / "data" / "other-database.json"


def load_site_data():
    spec = importlib.util.spec_from_file_location("other_database_site_data", ROOT / "scripts" / "other-database" / "site_data.py")
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
            p = copy / "totals.tsv"
            p.write_text(p.read_text(encoding="utf-8").replace("\t52\t", "\t53\t"), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(copy), "--check", str(BUNDLED)]), 1)

    def test_bundled_counts_match_headline_record_values_and_hold_no_genomes(self):
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        self.assertEqual(sum(s["scored"] for s in data["species"]), 13377)
        self.assertEqual(len(data["species"]), 9)
        self.assertEqual((data["totals"]["scm"]["plan_meets"], data["totals"]["scm"]["holds"], data["totals"]["scm"]["estimable"]), (6, 52, 84))
        self.assertEqual((data["totals"]["cart"]["plan_meets"], data["totals"]["cart"]["holds"], data["totals"]["cart"]["estimable"]), (5, 50, 84))
        self.assertEqual((data["totals"]["lasso"]["plan_meets"], data["totals"]["lasso"]["holds"], data["totals"]["lasso"]["estimable"]), (3, 17, 24))
        for l in S.LEARNERS:
            eight = [r for r in data["by_species"] if r["part"] == "eight" and r["learner"] == l]
            self.assertEqual(sum(r["holds"] for r in eight), data["totals"][l]["holds"])
            self.assertEqual(sum(r["estimable"] for r in eight), data["totals"][l]["estimable"])
        raw = BUNDLED.read_text(encoding="utf-8")
        self.assertNotIn("SAMN", raw)
        self.assertNotIn("GCA_", raw)
