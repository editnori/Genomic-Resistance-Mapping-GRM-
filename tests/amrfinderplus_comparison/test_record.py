"""The committed record under reference/amrfinderplus-comparison: its files, its sums, its roles and the review's
headline figures."""

import hashlib
import io
import re
import unittest
from contextlib import redirect_stderr, redirect_stdout

from tests.amrfinderplus_comparison.common import LEARNERS, PAGE, RECORD, ROOT, SCRIPT, SUITE, load_script, pairs, read_tsv

record = load_script()
SUITE_DIR = ROOT / SUITE

# The sha256 of SHA256SUMS: a change to any file of the record changes it.
SHA256SUMS_SHA256 = "a0cdfccfd7979b3e6367932c62a85a98aec940bc95fc02fd0c566dce8a63a7db"
FILES = {"README.md", "SHA256SUMS", "callers.tsv", "datasets.tsv", "differences.tsv", "figures.tsv", "intervals.tsv",
         "order.tsv", "settings.tsv", "sets.tsv", "summary.tsv"}
BUDGET = 3 * 1024 * 1024
NEW = "primary-new-235"
PR = "program_minus_known_gene_rule"
# The review's figures, which bind: (set, contrast, measure) -> cells.
REVIEW = {
    (NEW, PR, "ba"): {"datasets": "235", "ahead": "94", "level": "102", "behind": "39", "median": "0.000000",
                      "mean": "0.064236"},
    (NEW, "learner_minus_program:lineage:served_best", "ba"): {
        "ahead": "122", "level": "5", "behind": "108", "median": "0.003133", "median_species_lo": "-0.008135",
        "median_species_hi": "0.041124", "within_0.02_either_way": "56", "learner_at_most_0.02_ahead": "132"},
    ("rule-M-quinolone-28", "learner_minus_program:lineage:served_best", "ba"): {
        "datasets": "28", "learner_at_most_0.02_ahead": "16", "within_0.02_either_way": "11"},
    ("rule-M-quinolone-new-26", PR, "ba"): {"datasets": "26", "median": "0.330055"},
}
MEDIANS = {"ba": ["0.004751", "0.000000", "-0.006751", "-0.002604", "-0.003077"],
           "sens": ["-0.050847", "-0.078981", "-0.080645", "-0.071496", "-0.063131"],
           "spec": ["0.025111", "0.027344", "0.023417", "0.032154", "0.026634"]}


class Record(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.summary = {(r["set"], r["contrast"], r["measure"]): r for r in read_tsv(RECORD / "summary.tsv")}

    def test_files_and_sums(self):
        self.assertEqual({p.name for p in RECORD.iterdir()}, FILES)
        sums = (RECORD / "SHA256SUMS").read_bytes()
        self.assertEqual(hashlib.sha256(sums).hexdigest(), SHA256SUMS_SHA256)
        listed = dict(reversed(line.split("  ", 1)) for line in sums.decode("utf-8").splitlines())
        self.assertEqual(set(listed), FILES - {"SHA256SUMS"})
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)

    def test_check_rebuilds_every_derived_file(self):
        err, out = io.StringIO(), io.StringIO()
        with redirect_stderr(err), redirect_stdout(out):
            code = record.main(["--check", str(RECORD)])
        self.assertEqual(code, 0, err.getvalue())

    def test_shares_between_0_and_1_and_differences_between_minus_1_and_1(self):
        # Held apart from record.RANGES: every column that names a measure, an end or a share.
        share = re.compile(r"(sensitivity|specificity|balanced_accuracy|_exact_|_lineage_lo|_lineage_hi|_p$|^share_)")
        diff = re.compile(r"^(first_minus_second|lineage_lo|lineage_hi|median|q1|q3|mean|median_species_\w+|"
                          r"mean_species_\w+)$")
        checked = 0
        for name in ("callers.tsv", "differences.tsv", "intervals.tsv", "summary.tsv"):
            for r in read_tsv(RECORD / name):
                for col, cell in r.items():
                    lo = 0 if share.search(col) else -1 if diff.match(col) else None
                    if lo is not None and cell != "":
                        checked += 1
                        self.assertTrue(lo <= float(cell) <= 1, (name, col, cell))
        self.assertGreater(checked, 50000)

    def test_text_files_end_with_one_newline(self):
        for p in [PAGE, SCRIPT] + sorted(RECORD.iterdir()) + sorted(SUITE_DIR.glob("*.py")):
            text = p.read_text(encoding="utf-8")
            if text:
                self.assertTrue(text.endswith("\n") and not text.endswith("\n\n"), p.name)
                self.assertNotIn("\r", text, p.name)

    def test_size(self):
        total = 0
        for p in RECORD.iterdir():
            self.assertLess(p.stat().st_size, BUDGET, p.name)
            total += p.stat().st_size
        self.assertLess(total, BUDGET)

    def test_rows_and_renames(self):
        diffs = read_tsv(RECORD / "differences.tsv")
        self.assertEqual({r["reading"] for r in diffs}, {"tool_default"})
        callers = read_tsv(RECORD / "callers.tsv")
        readings = {r["caller"].split(":", 1)[1] for r in callers if r["caller"].startswith("program:")}
        self.assertEqual(readings, {"tool_default", "tool_default_intact", "genes_core", "genes_all_scopes"})
        for r in self.summary.values():
            empty = r["contrast"] == PR
            for col in ("learner_at_most_0.02_ahead", "share_learner_at_most_0.02_ahead", "within_0.02_either_way"):
                self.assertEqual(r[col] == "", empty, (r["set"], r["contrast"], col))

    def test_roles_are_the_plans(self):
        roles = {k: r["role"] for k, r in self.summary.items()}
        self.assertEqual(set(roles.values()), {"primary", "secondary", "exploratory"})
        for (name, con, m), r in roles.items():
            self.assertEqual(r, record.role(name, con, m), (name, con, m))
        primary = {k for k, v in roles.items() if v == "primary"}
        want = {(NEW, PR, m) for m in ("ba", "sens", "spec")}
        want |= {(NEW, "learner_minus_program:lineage:" + p, "ba") for p in LEARNERS + ("served_best",)}
        self.assertEqual(primary, want)
        for r in read_tsv(RECORD / "figures.tsv"):
            self.assertEqual(r["role"], "exploratory", r["key"])
        for r in self.summary.values():
            if ":random:" in r["contrast"]:
                self.assertNotEqual(r["role"], "primary")

    def test_the_reviews_figures(self):
        for key, cells in REVIEW.items():
            for col, value in cells.items():
                self.assertEqual(self.summary[key][col], value, (key, col))
        for m, values in MEDIANS.items():
            got = [self.summary[(NEW, "learner_minus_program:lineage:" + l, m)]["median"] for l in LEARNERS]
            self.assertEqual(got, values, m)
        figures = {r["key"]: r["value"] for r in read_tsv(RECORD / "figures.tsv")}
        self.assertEqual(figures["program_ahead_by_0.10_or_more"], "48")
        self.assertEqual(figures["of_them_quinolone"], "20")
        self.assertEqual(figures["level_with_the_same_call_on_every_genome"], "100")

    def test_sets(self):
        sets = {r["set"]: r for r in read_tsv(RECORD / "sets.tsv")}
        self.assertEqual([r["set"] for r in read_tsv(RECORD / "sets.tsv")], record.SET_NAMES)
        for name in record.SET_NAMES:
            self.assertEqual(sets[name]["datasets"], name.rsplit("-", 1)[1], name)
        self.assertEqual(sets["no-organism-option-12"]["species"], "2")
        datasets = read_tsv(RECORD / "datasets.tsv")
        self.assertEqual(len(datasets), int(pairs(RECORD / "settings.tsv")["datasets_scored"]))
        self.assertEqual(set(datasets[0]), set(record.COLUMNS["datasets.tsv"]))

    def test_settings_name_the_generator_and_the_draw(self):
        st = pairs(RECORD / "settings.tsv")
        self.assertEqual(st["generator"], "random.Random(seed string)")
        self.assertEqual(st["draw"], "randrange over the species of the set sorted by name")
        self.assertEqual(st["median"], "statistics.median")
        self.assertEqual((st["program_version"], st["database_version"]), ("4.2.7", "2026-08-07.1"))
        self.assertEqual(st["newest_database_not_run"], "2026-10-29.1")

    def test_readme_carries_the_limit_and_the_catalog_sentence(self):
        readme = (RECORD / "README.md").read_text(encoding="utf-8")
        self.assertIn("This is a retrospective comparison inside one public database, BV-BRC. It is not a prospective "
                      "test and not external validation.", readme)
        self.assertIn("Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both "
                      "read is newer than most of these genomes and may have been curated from the same public "
                      "isolates, which favours both against the learners.", readme)


if __name__ == "__main__":
    unittest.main()
