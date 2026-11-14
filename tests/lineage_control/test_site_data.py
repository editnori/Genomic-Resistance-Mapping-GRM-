"""The within-lineage verdicts the web app bundles (web/src/data/lineage-control.json).

The file must be exactly what scripts/lineage-control/site_data.py writes
from the table under reference/lineage-control, so the site attaches the
table's verdicts and no others.
"""

import importlib.util
import json
import pathlib
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
TABLE = ROOT / "reference" / "lineage-control" / "rev14-274-20261112.tsv"
BUNDLED = ROOT / "web" / "src" / "data" / "lineage-control.json"


def load_site_data():
    spec = importlib.util.spec_from_file_location("lineage_control_site_data", ROOT / "scripts" / "lineage-control" / "site_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_site_data()


class BundledVerdicts(unittest.TestCase):
    def changed(self, tmp, change):
        """A copy of the table with one cell of its first data row changed."""
        lines = TABLE.read_text(encoding="utf-8").split("\n")
        head = lines[0].split("\t")
        cells = lines[1].split("\t")
        column, value = change(dict(zip(head, cells)))
        cells[head.index(column)] = value
        lines[1] = "\t".join(cells)
        path = pathlib.Path(tmp) / TABLE.name
        path.write_text("\n".join(lines), encoding="utf-8")
        return path

    def test_bundled_file_is_what_the_table_gives(self):
        self.assertEqual(BUNDLED.read_text(encoding="utf-8"), S.build(TABLE))
        self.assertEqual(S.main(["site_data.py", str(TABLE), "--check", str(BUNDLED)]), 0)

    def test_check_fails_when_a_verdict_a_rule_or_a_fingerprint_changes(self):
        changes = [
            lambda row: ("verdict", "cannot_test" if row["verdict"] == "seen" else "seen"),
            lambda row: ("feature", str(int(row["feature"]) + 1)),
            lambda row: ("presence", "0" if row["presence"] == "1" else "1"),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for change in changes:
                self.assertEqual(S.main(["site_data.py", str(self.changed(tmp, change)), "--check", str(BUNDLED)]), 1)

    def test_a_table_that_is_not_one_release_in_order_is_refused(self):
        refused = [
            lambda row: ("revision", "15"),  # two revisions in one table
            lambda row: ("position", "2"),  # a term out of order
            lambda row: ("verdict", "proved"),  # no such verdict
            lambda row: ("input_fingerprint", "0" * 16),  # two fingerprints for one dataset
        ]
        with tempfile.TemporaryDirectory() as tmp:
            for change in refused:
                with self.assertRaises(SystemExit):
                    S.build(self.changed(tmp, change))

    def test_the_bundle_names_its_release_and_keeps_each_term_in_the_models_order(self):
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        self.assertEqual((data["release"], data["evaluation_revision"], data["min_strata"], data["alpha"], data["terms"]), ("rev14-274-20261112", 14, 5, 0.05, 862))
        self.assertEqual(data["counts"], {"seen": 327, "not_shown": 6, "cannot_test": 529, "reversed": 0})
        self.assertEqual(len(data["datasets"]), 274)
        first = data["datasets"]["acinetobacter_baumannii__amikacin"]
        self.assertEqual(first["fingerprint"], "dd48cc69d891d65d")
        self.assertEqual(first["terms"][:2], [[181099, True, "seen"], [47267, True, "cannot_test"]])


if __name__ == "__main__":
    unittest.main()
