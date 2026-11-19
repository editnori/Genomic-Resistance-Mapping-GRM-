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
            self.assertTrue(all(("report_at_edge" in r) == (r["learner"] == "known_gene") for r in p["rows"]))

    def test_bundled_rows_give_the_documents_limits(self):
        # docs/outside-study-results.md: 12 of the 39 known-gene holds have a report figure at 0% or 100%;
        # the same-study population is four studies, each at least 99.9% of a dataset's resistant genomes.
        data = {p["id"]: p for p in json.loads(BUNDLED.read_text(encoding="utf-8"))["populations"]}
        rule = [r for r in data["primary"]["rows"] if r["learner"] == "known_gene"]
        held = [r for r in rule if {r["sensitivity"], r["specificity"]} <= {"consistent", "above"}]
        self.assertEqual((len(held), sum(r["report_at_edge"] for r in held)), (39, 12))
        same = [d for d in data["same-study"]["datasets"] if d["estimates"]]
        self.assertEqual(sum(d["top_study_resistant_share"] >= 0.999 for d in same), 49)
        self.assertEqual(len({d["top_study_resistant"] for d in same}), 4)


if __name__ == "__main__":
    unittest.main()
