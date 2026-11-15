"""The within-lineage verdicts the web app bundles (web/src/data/lineage-control.json).

The file must be exactly what scripts/lineage-control/site_data.py writes
from the table under reference/lineage-control, so the site attaches the
table's verdicts and no others. The counts the bundle gives are compared with
the table's rows counted here; the one table whose expected bundle is typed
out is a synthetic one of three rows.
"""

import collections
import csv
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
        with TABLE.open(encoding="utf-8", newline="") as f:
            table = list(csv.DictReader(f, delimiter="\t"))
        level = (S.verdicts.MIN_STRATA, S.verdicts.ALPHA)
        self.assertEqual((data["release"], data["evaluation_revision"], data["min_strata"], data["alpha"], data["terms"]), (TABLE.stem, int(table[0]["revision"]), *level, len(table)))
        # The counts are the table's rows counted, every verdict of the script named, and they add up to the terms.
        counted = collections.Counter(r["verdict"] for r in table)
        self.assertEqual(data["counts"], {v: counted[v] for v in S.verdicts.VERDICTS})
        self.assertEqual(sum(data["counts"].values()), data["terms"])
        by_dataset = collections.defaultdict(list)
        for r in table:
            by_dataset[r["dataset"]].append(r)
        self.assertEqual(list(data["datasets"]), sorted(by_dataset))
        for dataset, terms in by_dataset.items():
            entry = data["datasets"][dataset]
            self.assertEqual(entry["fingerprint"], terms[0]["input_fingerprint"], dataset)
            self.assertEqual(entry["terms"], [[int(r["feature"]), r["presence"] == "1", r["verdict"]] for r in terms], dataset)

    def test_a_synthetic_table_of_three_rows_gives_the_bundle_typed_out_here(self):
        # Synthetic: two made-up datasets, three terms. The expected text is typed by hand, so the script's
        # layout is held by something other than the script.
        head = "\t".join(S.verdicts.COLUMNS)
        made_up = [
            "b__drug\t14\t1\t7\t0\tffff000011112222\tsha-b\t40\t40\t2\t0.5\t0.4\t\tcannot_test",
            "a__drug\t14\t1\t12\t1\t0123456789abcdef\tsha-a\t60\t60\t6\t9\t0.001\t0.002\tseen",
            "a__drug\t14\t2\t3\t0\t0123456789abcdef\tsha-a\t60\t25\t5\t1\t0.3\t0.3\tnot_shown",
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "made-up-release.tsv"
            path.write_text("\n".join([head, *made_up]) + "\n", encoding="utf-8")
            self.assertEqual(
                S.build(path),
                '{"release":"made-up-release","evaluation_revision":14,"min_strata":5,"alpha":0.05,"terms":3,'
                '"counts":{"seen":1,"not_shown":1,"cannot_test":1,"reversed":0},'
                '"datasets":{"a__drug":{"fingerprint":"0123456789abcdef","terms":[[12,true,"seen"],[3,false,"not_shown"]]},'
                '"b__drug":{"fingerprint":"ffff000011112222","terms":[[7,false,"cannot_test"]]}}}\n',
            )


if __name__ == "__main__":
    unittest.main()
