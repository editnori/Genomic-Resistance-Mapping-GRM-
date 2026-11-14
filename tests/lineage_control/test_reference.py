"""The committed within-lineage verdict table (reference/lineage-control) and its README.

The suite reads no release and tests no term again. It checks that the table
holds one release's terms under the script's columns, that each verdict
follows from the row's own figures, that each report hash is the release's
own checksum line, and that the README gives the files' hashes and counts.
"""

import collections
import csv
import hashlib
import importlib.util
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "lineage-control"
RELEASE = "rev14-274-20261112"
TABLE = REFERENCE / f"{RELEASE}.tsv"
REPORTS = REFERENCE / f"{RELEASE}.reports.sha256"


def load_verdicts():
    spec = importlib.util.spec_from_file_location("lineage_control_verdicts", ROOT / "scripts" / "lineage-control" / "verdicts.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = load_verdicts()


def rows():
    with TABLE.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


class Reference(unittest.TestCase):
    def test_the_table_holds_one_release_under_the_script_columns(self):
        with TABLE.open(encoding="utf-8") as f:
            self.assertEqual(f.readline().rstrip("\n").split("\t"), list(V.COLUMNS))
        table = rows()
        self.assertEqual((len(table), len({r["dataset"] for r in table}), {r["revision"] for r in table}), (862, 274, {"14"}))
        # A dataset's terms are numbered from 1 in the model's order, under one fingerprint and one report.
        by_dataset = collections.defaultdict(list)
        for r in table:
            by_dataset[r["dataset"]].append(r)
        for dataset, terms in by_dataset.items():
            self.assertEqual([int(r["position"]) for r in terms], list(range(1, len(terms) + 1)), dataset)
            self.assertEqual(len({(r["input_fingerprint"], r["report_sha256"], r["isolates"]) for r in terms}), 1, dataset)
            self.assertEqual(int(terms[0]["context_isolates"]), int(terms[0]["isolates"]), dataset)
        self.assertEqual(list(by_dataset), sorted(by_dataset))

    def test_the_counts_are_the_ones_the_site_and_the_readme_give(self):
        table = rows()
        self.assertEqual(collections.Counter(r["verdict"] for r in table), {"seen": 327, "not_shown": 6, "cannot_test": 529})
        leading = collections.Counter(r["verdict"] for r in table if r["position"] == "1")
        self.assertEqual(leading, {"seen": 179, "not_shown": 1, "cannot_test": 94})
        verdicts = collections.defaultdict(set)
        for r in table:
            verdicts[r["dataset"]].add(r["verdict"])
        self.assertEqual(sum(v == {"cannot_test"} for v in verdicts.values()), 90)

    def test_each_verdict_follows_from_its_row(self):
        for r in rows():
            where = f"{r['dataset']} term {r['position']}"
            strata = int(r["informative_strata"])
            if strata < V.MIN_STRATA or r["p"] == "":
                self.assertEqual((r["verdict"], r["q"]), ("cannot_test", ""), where)
                continue
            q, excess = float(r["q"]), float(r["excess"])
            self.assertGreaterEqual(q, float(r["p"]) * (1 - 1e-5), where)
            want = "not_shown" if q >= V.ALPHA else "seen" if excess > 0 else "reversed"
            self.assertEqual(r["verdict"], want, where)

    def test_each_report_hash_is_the_release_checksum_line(self):
        lines = {}
        for line in REPORTS.read_text(encoding="utf-8").splitlines():
            digest, name = line.split("  ")
            self.assertRegex(name, r"^\./datasets/[^/]+\.json$")
            lines[name[len("./datasets/") : -len(".json")]] = digest
        self.assertEqual(len(lines), 274)
        table = rows()
        self.assertEqual({r["dataset"] for r in table}, set(lines))
        for r in table:
            self.assertEqual(r["report_sha256"], lines[r["dataset"]], r["dataset"])

    def test_the_readme_gives_each_file_hash(self):
        readme = (REFERENCE / "README.md").read_text(encoding="utf-8")
        listed = dict(reversed(m.groups()) for m in re.finditer(r"^([0-9a-f]{64})  (\S+)$", readme, re.M))
        names = sorted(p.name for p in REFERENCE.iterdir() if p.name != "README.md")
        self.assertEqual(sorted(listed), names)
        for name in names:
            self.assertEqual(hashlib.sha256((REFERENCE / name).read_bytes()).hexdigest(), listed[name], name)

    def test_the_readme_gives_the_counts(self):
        readme = (REFERENCE / "README.md").read_text(encoding="utf-8")
        for words in ("862 rule terms of 274 reports", "| `seen` | 327 | 179 |", "| `not_shown` | 6 | 1 |", "| `cannot_test` | 529 | 94 |", "| `reversed` | 0 | 0 |", "90 datasets"):
            self.assertIn(words, readme)


if __name__ == "__main__":
    unittest.main()
