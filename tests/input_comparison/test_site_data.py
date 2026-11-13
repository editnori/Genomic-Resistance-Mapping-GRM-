"""The input comparison's numbers the web app bundles (web/src/data/input-comparison.json).

The file must be exactly what scripts/input-comparison/site_data.py writes from the reference reading, and the
reference reading must be the one its README names and give the runner document's "By dataset" table.
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
        self.assertEqual(BUNDLED.read_text(encoding="utf-8"), S.build(READING))
        self.assertEqual(S.main(["site_data.py", str(READING), "--check", str(BUNDLED)]), 0)

    def test_check_fails_when_the_reading_changes(self):
        reading = json.loads(READING.read_text(encoding="utf-8"))
        reading["panel"]["by_dataset"]["cap-1000"]["worse"] += 1
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "second-reading.json"
            changed.write_text(json.dumps(reading), encoding="utf-8")
            self.assertEqual(S.main(["site_data.py", str(changed), "--check", str(BUNDLED)]), 1)

    def test_firm_decisions_must_add_up(self):
        reading = json.loads(READING.read_text(encoding="utf-8"))
        reading["panel"]["firm"]["worse"] += 1
        with tempfile.TemporaryDirectory() as tmp:
            changed = pathlib.Path(tmp) / "second-reading.json"
            changed.write_text(json.dumps(reading), encoding="utf-8")
            with self.assertRaises(SystemExit):
                S.build(changed)

    def test_the_reference_is_the_file_its_readme_names(self):
        named = re.search(r"`second-reading.json` with sha256\s+`([0-9a-f]{64})`", README.read_text(encoding="utf-8"))
        self.assertIsNotNone(named)
        self.assertEqual(hashlib.sha256(READING.read_bytes()).hexdigest(), named.group(1))

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
