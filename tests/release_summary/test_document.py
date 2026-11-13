"""docs/release-summary.md and the reference outputs it quotes.

Each table of the document is rebuilt here from the files under reference/release-summary, and the document must hold
it line for line. The sentences that carry a count are checked one by one."""

import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "release-summary"
DOCUMENT = ROOT / "docs" / "release-summary.md"
A, B = "rev12-274-20261105", "rev14-274-20261112"
NAMES = {"scm": "Set covering machine", "cart": "Decision tree", "rf": "Random forest", "gbt": "Boosted trees",
         "lasso": "Lasso", "best": "Best learner", "catalog": "Known-gene rule"}


def load(release, name):
    return json.loads((REFERENCE / release / name).read_text(encoding="utf-8"))


def cell(r):
    lo, hi = r["ci95_species"]
    return f"{r['value']:.4f} ({lo:.4f} to {hi:.4f})"


def tables():
    """{name: the table's lines}, each built from the reference outputs of both releases."""
    ua, ub = load(A, "uncertainty.json"), load(B, "uncertainty.json")
    ca, cb = load(A, "catalog.json"), load(B, "catalog.json")
    out = {}
    rows = ["| Learner | Split | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- | --- |"]
    for l in ["scm", "cart", "rf", "gbt", "lasso", "best"]:
        for s in ("random", "lineage"):
            a, b = ua["medians"][f"{l}:{s}"], ub["medians"][f"{l}:{s}"]
            rows.append(f"| {NAMES[l]} | {s} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["medians"] = rows
    rows = ["| Learner | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- |"]
    for l in ["scm", "cart", "rf", "gbt", "lasso", "best"]:
        a, b = ua["gaps"][f"{l}:gap"], ub["gaps"][f"{l}:gap"]
        rows.append(f"| {NAMES[l]} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["gaps"] = rows
    rows = ["| Set covering machine minus | Datasets | Revision 12 | Revision 14 | Set covering machine at least as high |",
            "| --- | --- | --- | --- | --- |"]
    for l in ["cart", "rf", "gbt", "lasso"]:
        a, b = ua["paired_lineage"][f"scm-{l}:lineage"], ub["paired_lineage"][f"scm-{l}:lineage"]
        rows.append(f"| {NAMES[l]} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} | "
                    f"{a['first_at_least_second']} / {b['first_at_least_second']} |")
    out["paired"] = rows
    rows = ["| Statistic | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- |"]
    for k, label in (("median_lineage_of_best", "Lineage score of the best learner"),
                     ("median_lineage_of_oracle", "Highest lineage score of the five"),
                     ("median_oracle_minus_best", "Highest of the five minus the best learner")):
        a, b = ua["selection"][k], ub["selection"][k]
        rows.append(f"| {label} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["selection"] = rows
    rows = ["| Learner | Split | Datasets | Point within both | Shown to meet both | Shown to miss one | Not decided |",
            "| --- | --- | --- | --- | --- | --- | --- |"]
    for who in ["best", "scm", "cart", "rf", "gbt", "lasso", "catalog"]:
        for s in ("random", "lineage"):
            a, b = ua["error_bars"]["rows"][f"{who}:{s}"], ub["error_bars"]["rows"][f"{who}:{s}"]
            f = [f"{a[k]} / {b[k]}" for k in ("n", "point_meets_both", "shown_to_meet_both", "shown_to_miss_one", "not_decided")]
            rows.append(f"| {NAMES[who]} | {s} | " + " | ".join(f) + " |")
    out["error_limits"] = rows
    rows = ["| Rule's balanced accuracy | Datasets | Best learner ahead, random split | Best learner ahead, lineage split "
            "| Rule ahead or within 2 points, lineage split | Median lineage difference |",
            "| --- | --- | --- | --- | --- | --- |"]
    for k in ca["bands"]:
        a, b = ca["bands"][k], cb["bands"][k]
        f = lambda x: f"{a[x]} / {b[x]}"
        rows.append(f"| {k} | {f('datasets')} | {f('learner_ahead_random')} | {f('learner_ahead_lineage')} of "
                    f"{f('with_a_lineage_score')} | {f('rule_ahead_or_within_2pt_lineage')} | "
                    f"{a['median_lineage_difference']} / {b['median_lineage_difference']} |")
    out["bands"] = rows
    return out


class Document(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = DOCUMENT.read_text(encoding="utf-8")

    def test_each_table_is_the_reference_files(self):
        for name, rows in tables().items():
            with self.subTest(table=name):
                self.assertIn("\n".join(rows) + "\n", self.text)

    def test_counted_sentences(self):
        ea, eb = load(A, "uncertainty.json")["error_bars"], load(B, "uncertainty.json")["error_bars"]
        self.assertEqual(ea["resistant_needed_for_zero_misses_to_show_vme"], eb["resistant_needed_for_zero_misses_to_show_vme"])
        self.assertEqual(ea["datasets_with_that_many_resistant"], eb["datasets_with_that_many_resistant"])
        ra, rb = load(A, "catalog.json")["rule"], load(B, "catalog.json")["rule"]
        self.assertEqual(ra, rb)
        self.assertEqual(ra["calls_nothing_resistant"] + ra["calls_everything_resistant"], ra["balanced_accuracy_exactly_half"])
        ia, ib = load(A, "index-counts.json"), load(B, "index-counts.json")
        sa = load(A, "uncertainty.json")
        sentences = [
            f"Each release has {sa['datasets']} datasets in {sa['species']} species.",
            f"a dataset needs at least {ea['resistant_needed_for_zero_misses_to_show_vme']} resistant genomes. "
            f"{ea['datasets_with_that_many_resistant']} datasets of each release have that many.",
            f"The rule has an entry on {ra['with_an_entry']} datasets. Its balanced accuracy is exactly 0.5 on "
            f"{ra['balanced_accuracy_exactly_half']} of them: it calls nothing resistant on {ra['calls_nothing_resistant']} "
            f"and everything resistant on {ra['calls_everything_resistant']}.",
            f"In both releases the counts of all {ia['rates_compared']} error rates of the {ia['rows']} intervals rows "
            f"rebuild from the index with {ia['mismatches']} mismatches.",
            f"In revision 14, {ib['flagged_rows']} learner rows carry the flag, and each of them has a no-call count.",
        ]
        self.assertEqual((ia["mismatches"], ib["mismatches"]), (0, 0))
        self.assertEqual(ib["flagged_rows"], ib["flagged_with_no_call"])
        for s in sentences:
            with self.subTest(sentence=s[:60]):
                self.assertIn(s, " ".join(self.text.split()))


if __name__ == "__main__":
    unittest.main()
