"""The outside-study rows the web app bundles (web/src/data/outside-study.json).

The file must be exactly what scripts/outside-study-results/site_data.py
writes from reference/outside-study, so the site and the documents read the
same rows by the same rules.
"""

import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
RESULTS = ROOT / "reference" / "outside-study"
BUNDLED = ROOT / "web" / "src" / "data" / "outside-study.json"


def load_site_data():
    spec = importlib.util.spec_from_file_location("outside_study_site_data", ROOT / "scripts" / "outside-study-results" / "site_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_site_data()


class BundledRows(unittest.TestCase):
    def test_bundled_file_is_what_the_files_give(self):
        self.assertEqual(BUNDLED.read_text(encoding="utf-8"), S.build(RESULTS))
        self.assertEqual(S.main(["site_data.py", str(RESULTS), "--check", str(BUNDLED)]), 0)

    def test_check_fails_when_a_reading_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp) / "outside-study"
            shutil.copytree(RESULTS, root)
            metrics = root / "primary" / "escherichia_coli" / "metrics.tsv"
            lines = metrics.read_text(encoding="utf-8").split("\n")
            head = lines[0].split("\t")
            at = head.index("reading_specificity")
            row = next(i for i, line in enumerate(lines[1:], 1) if line.split("\t")[head.index("learner")] == "scm")
            cells = lines[row].split("\t")
            cells[at] = "consistent" if cells[at] == "below" else "below"
            lines[row] = "\t".join(cells)
            metrics.write_text("\n".join(lines), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(root), "--check", str(BUNDLED)]), 1)

    def test_rows_keep_the_populations_apart(self):
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        self.assertEqual([p["id"] for p in data["populations"]], ["primary", "same-study"])
        for p in data["populations"]:
            ids = {d["id"] for d in p["datasets"]}
            self.assertTrue(all(r["dataset"] in ids for r in p["rows"]))
            self.assertTrue(all(("figures" in r) == (r["learner"] == "scm") for r in p["rows"]))


if __name__ == "__main__":
    unittest.main()
