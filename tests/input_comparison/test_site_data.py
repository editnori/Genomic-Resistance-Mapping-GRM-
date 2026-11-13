"""The input comparison's numbers the web app bundles (web/src/data/input-comparison.json).

The file must be exactly what scripts/input-comparison/site_data.py writes from the two reference readings, each
reference reading must be the one its README names, the second reading must give the runner document's "By dataset"
table, and the reading by learner must give the numbers the document's "On H-1" lines name.
"""

import hashlib
import importlib.util
import json
import pathlib
import re
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
READING = ROOT / "reference" / "input-comparison" / "second-reading.json"
LEARNERS = ROOT / "reference" / "input-comparison" / "learner-reading.json"
README = ROOT / "reference" / "input-comparison" / "README.md"
BUNDLED = ROOT / "web" / "src" / "data" / "input-comparison.json"
RUNNER = ROOT / "docs" / "input-comparison-runner.md"


def load_site_data():
    spec = importlib.util.spec_from_file_location("input_comparison_site_data", ROOT / "scripts" / "input-comparison" / "site_data.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


S = load_site_data()


class BundledReading(unittest.TestCase):
    def test_bundled_file_is_what_the_reading_gives(self):
        self.assertEqual(BUNDLED.read_text(encoding="utf-8"), S.build(READING, LEARNERS))
        self.assertEqual(S.main(["site_data.py", str(READING), str(LEARNERS), "--check", str(BUNDLED)]), 0)

    def test_check_fails_when_the_reading_changes(self):
        reading = json.loads(READING.read_text(encoding="utf-8"))
        reading["panel"]["by_dataset"]["cap-1000"]["worse"] += 1
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "second-reading.json"
            changed.write_text(json.dumps(reading), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(changed), str(LEARNERS), "--check", str(BUNDLED)]), 1)

    def test_firm_decisions_must_add_up(self):
        reading = json.loads(READING.read_text(encoding="utf-8"))
        reading["panel"]["firm"]["worse"] += 1
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "second-reading.json"
            changed.write_text(json.dumps(reading), encoding="utf-8")
            with self.assertRaises(SystemExit):
                S.build(changed, LEARNERS)

    def test_check_fails_when_the_reading_by_learner_changes(self):
        reading = json.loads(LEARNERS.read_text(encoding="utf-8"))
        reading["tables"]["lineage"]["cap-4000"]["scm"]["better"] += 1
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "learner-reading.json"
            changed.write_text(json.dumps(reading), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(READING), str(changed), "--check", str(BUNDLED)]), 1)

    def test_the_references_are_the_files_the_readme_names(self):
        text = README.read_text(encoding="utf-8")
        for name, path in (("second-reading.json", READING), ("learner-reading.json", LEARNERS)):
            named = re.search(rf"`{re.escape(name)}` with sha256\s+`([0-9a-f]{{64}})`", text)
            self.assertIsNotNone(named, name)
            self.assertEqual(hashlib.sha256(path.read_bytes()).hexdigest(), named.group(1), name)

    def test_bundled_learner_part_gives_the_documents_on_h1_numbers(self):
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))
        arms = {a["id"]: a["learners"] for a in data["arms"]}
        scm = arms["cap-4000"]["scm"]
        self.assertEqual((round(scm["share"]["units"], 3), round(scm["share"]["genomes"], 3)), (0.739, 0.694))
        self.assertEqual((scm["better"], scm["worse"], scm["unchanged"], scm["datasets"], round(scm["p"], 4)), (32, 16, 1, 49, 0.0293))
        self.assertEqual([arms["cap-1000"][l]["worse"] for l in S.LEARNERS], [21, 24, 28, 25, 19])
        self.assertEqual({arms["cap-1000"][l]["datasets"] for l in S.LEARNERS}, {33})
        named = data["learner_reading"]["named"][0]
        self.assertEqual(named["dataset"], "mycobacterium_tuberculosis__clofazimine")
        self.assertEqual((named["default"]["tp"], named["default"]["fn"], named["arm_calls"]["tp"], named["arm_calls"]["fn"]), (381, 198, 105, 474))
        self.assertEqual(round(named["estimate"], 4), 2.7935)
        # The calls are summed over the repeats: 579 calls on 193 resistant genomes, three times each.
        self.assertEqual((named["repeats"], (named["default"]["tp"] + named["default"]["fn"]) // named["repeats"]), (3, 193))
        self.assertIn("decides nothing", data["learner_reading"]["exploratory"])

    def test_bundled_counts_give_the_documents_by_dataset_table(self):
        doc = RUNNER.read_text(encoding="utf-8").splitlines()
        data = json.loads(BUNDLED.read_text(encoding="utf-8"))

        def cell(c):
            return "none made" if c["datasets"] == 0 else f"{c['better']}, {c['worse']} of {c['datasets']}"

        for a in data["arms"]:
            row = f"| {a['id']} | " + " | ".join(cell(c) for c in a["species"]) + f" | {cell(a['all'])} | {a['all']['p']:.4f} | {a['all']['median']:.4f} |"
            self.assertIn(row, doc)


if __name__ == "__main__":
    unittest.main()
