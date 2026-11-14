"""The committed within-lineage verdict table (reference/lineage-control) and its README.

The suite reads no release and tests no term again. It checks that the table
holds one release's terms under the script's columns, that each verdict
follows from the row's own figures, that each report hash is the release's
own checksum line, and that the README gives the files' hashes and the
table's counts. No count is typed here: each is counted from the table, and
the table's bytes are held by the hash the README gives.

`scripts/lineage-control/verdicts.py RELEASE_DIR --check TABLE` compares the
table with what the script prints for the unpacked release (docs/ci.md, "Not
in CI"); tests/lineage_control/test_verdicts.py tests that mode on a
synthetic release.
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


def checksum_lines():
    """The release's checksum line of each report: dataset id to sha256."""
    lines = {}
    for line in REPORTS.read_text(encoding="utf-8").splitlines():
        digest, name = line.split("  ")
        if not (name.startswith("./datasets/") and name.endswith(".json") and "/" not in name[len("./datasets/") :]):
            raise AssertionError(f"{REPORTS.name}: {name!r} is not a report of the release")
        lines[name[len("./datasets/") : -len(".json")]] = digest
    return lines


def count_lines(table):
    """The README's lines that give a count, written from the table's rows."""
    revisions = sorted({r["revision"] for r in table})
    reports = len({r["dataset"] for r in table})
    terms = collections.Counter(r["verdict"] for r in table)
    leading = collections.Counter(r["verdict"] for r in table if r["position"] == "1")
    verdicts = collections.defaultdict(set)
    for r in table:
        verdicts[r["dataset"]].add(r["verdict"])
    untestable = sum(v == {"cannot_test"} for v in verdicts.values())
    tested = len(table) - terms["cannot_test"]
    return {
        "input": f"| `{TABLE.name}` | `{RELEASE}.tar.zst` | `{{asset}}` | {', '.join(revisions)} | {reports} reports |",
        "checked": f"On 2026-11-14 each of the {reports} equalled its line in the release's own checksum list, `{RELEASE}.sha256`.",
        "kept": f"`{REPORTS.name}` holds those {reports} lines as the release wrote them",
        "rows": f"One row for each rule term: {len(table)} rule terms of {reports} reports.",
        "verdicts": [f"| `{v}` | {terms[v]} | {leading[v]} |" for v in V.VERDICTS],
        "tested": f"Of the {tested} terms that could be tested, {terms['seen']} are `seen`.",
        "untestable": f"No term of {untestable} datasets can be tested.",
    }


class Reference(unittest.TestCase):
    def test_the_table_holds_one_release_under_the_script_columns(self):
        with TABLE.open(encoding="utf-8") as f:
            self.assertEqual(f.readline().rstrip("\n").split("\t"), list(V.COLUMNS))
        table = rows()
        self.assertGreater(len(table), 0)
        self.assertEqual(len({r["revision"] for r in table}), 1)
        # One row a term and no empty report: as many datasets as the release's checksum list has reports.
        self.assertEqual(len({r["dataset"] for r in table}), len(checksum_lines()))
        self.assertEqual({r["verdict"] for r in table} - set(V.VERDICTS), set())
        # A dataset's terms are numbered from 1 in the model's order, under one fingerprint and one report.
        by_dataset = collections.defaultdict(list)
        for r in table:
            by_dataset[r["dataset"]].append(r)
        for dataset, terms in by_dataset.items():
            self.assertEqual([int(r["position"]) for r in terms], list(range(1, len(terms) + 1)), dataset)
            self.assertEqual(len({(r["input_fingerprint"], r["report_sha256"], r["isolates"]) for r in terms}), 1, dataset)
            self.assertEqual(int(terms[0]["context_isolates"]), int(terms[0]["isolates"]), dataset)
        self.assertEqual(list(by_dataset), sorted(by_dataset))

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
        lines = checksum_lines()
        self.assertEqual(len(lines), len(REPORTS.read_text(encoding="utf-8").splitlines()))
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

    def test_the_readme_gives_the_counts_the_table_gives(self):
        # Every line of the README that holds a count is written again here from the table's rows and compared
        # with the README's own, so a README left behind by a new table fails, and so does a count typed wrong.
        readme = (REFERENCE / "README.md").read_text(encoding="utf-8")
        want = count_lines(rows())
        asset = re.search(r"^\| `" + re.escape(TABLE.name) + r"` \| `[^`]+` \| `([0-9a-f]{64})` \|", readme, re.M)
        self.assertIsNotNone(asset)
        self.assertIn(want["input"].format(asset=asset.group(1)), readme.splitlines())
        self.assertIn(want["checked"], readme)
        self.assertIn(want["kept"], readme)
        self.assertEqual(re.findall(r"^One row for each rule term: .*$", readme, re.M), [want["rows"]])
        self.assertEqual(re.findall(r"^\| `\w+` \| \d+ \| \d+ \|", readme, re.M), want["verdicts"])
        self.assertEqual(re.findall(r"^Of the \d+ terms that could be tested, .*$", readme, re.M), [want["tested"]])
        self.assertEqual(re.findall(r"^No term of \d+ datasets can be tested\.$", readme, re.M), [want["untestable"]])

    def test_the_readme_says_the_fourth_verdict_is_an_addition_and_what_the_site_does_with_it(self):
        readme = (REFERENCE / "README.md").read_text(encoding="utf-8")
        self.assertIn("`reversed`, is an addition to those rules", readme)
        self.assertIn("The inspector has no words for `reversed`", readme)
        # The page's own list of worded verdicts is the three fixed ones.
        app = (ROOT / "web" / "src" / "data" / "derive" / "lineage.ts").read_text(encoding="utf-8")
        self.assertIn('const WORDED: readonly string[] = ["seen", "not_shown", "cannot_test"];', app)
        self.assertEqual(V.VERDICTS, ("seen", "not_shown", "cannot_test", "reversed"))


if __name__ == "__main__":
    unittest.main()
