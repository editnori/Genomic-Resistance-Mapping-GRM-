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
    ga, gb = load(A, "guidance.json"), load(B, "guidance.json")
    rows = ["| Split | Learner | Datasets | Meet the guidance's wording | Shown to meet the release's screen |",
            "| --- | --- | --- | --- | --- |"]
    for s in ("random", "lineage"):
        for who in ["best", "scm", "cart", "rf", "gbt", "lasso", "catalog"]:
            a, b = ga["rows"][f"{who}:{s}"], gb["rows"][f"{who}:{s}"]
            rows.append(f"| {s} | {NAMES[who]} | {a['n']} / {b['n']} | {a['meets_wording']} / {b['meets_wording']} | "
                        f"{a['screen']['shown_to_meet_both']} / {b['screen']['shown_to_meet_both']} |")
    out["guidance"] = rows
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

    def test_guidance_sentences(self):
        ga, gb = load(A, "guidance.json"), load(B, "guidance.json")
        ra, rb = ga["rows"], gb["rows"]
        learner_rows = [f"{l}:{s}" for l in ("best", "scm", "cart", "rf", "gbt", "lasso") for s in ("random", "lineage")]
        both = lambda f: f(ra) if f(ra) == f(rb) else self.fail(f"the releases differ: {f(ra)} and {f(rb)}")
        meets = both(lambda r: (r["best:random"]["meets_wording"], r["best:lineage"]["meets_wording"]))
        shown = both(lambda r: (r["best:random"]["screen"]["shown_to_meet_both"], r["best:lineage"]["screen"]["shown_to_meet_both"]))
        self.assertTrue(all(ra[k]["meets_wording"] == rb[k]["meets_wording"] for k in ra))
        self.assertEqual(ga["best_meets_wording"], gb["best_meets_wording"])
        # revision 12 leaves out no-call genomes on boosted trees and lasso rows only; revision 14 on none
        fewer = {k: r["fewer_than_every_genome"] for k, r in ra.items() if r["fewer_than_every_genome"]}
        self.assertEqual(sorted(fewer), ["gbt:lineage", "gbt:random", "lasso:lineage", "lasso:random"])
        self.assertEqual(sum(fewer.values()), load(A, "index-counts.json")["unflagged_with_no_call"])
        self.assertEqual(sum(r["fewer_than_every_genome"] for r in rb.values()), 0)
        flagged = load(B, "index-counts.json")["flagged_rows"]
        table_8 = both(lambda r: (r["best:random"]["vme_n_under_table_8"], r["best:lineage"]["vme_n_under_table_8"],
                                  r["catalog:random"]["vme_n_under_table_8"], r["catalog:lineage"]["vme_n_under_table_8"]))
        first = ga["criteria"]["table_8_first_resistant"]
        self.assertEqual(first, gb["criteria"]["table_8_first_resistant"])
        tb = lambda r, k: r[k]["n"] - r[k]["without_tuberculosis"]["n"]
        tb_learners = both(lambda r: {tb(r, k) for k in learner_rows})
        tb_rule = both(lambda r: {tb(r, "catalog:random"), tb(r, "catalog:lineage")})
        self.assertEqual((len(tb_learners), len(tb_rule)), (1, 1))
        for r in (ra, rb):
            self.assertTrue(all(x["meets_wording"] == x["without_tuberculosis"]["meets_wording"] for x in r.values()))
        sentences = [
            f"Its first row is {first} resistant organisms with 0 discrepancies and an interval of 0.00 to 7.40%.",
            f"By the guidance's wording, the best learner meets on {meets[0]} datasets of the random split and {meets[1]} of "
            f"the lineage split, in both releases. By the release's screen, it is shown to meet on {shown[0]} and {shown[1]}.",
            f"Of the rows read here, {sum(fewer.values())} count fewer genomes than their class: boosted trees on "
            f"{fewer['gbt:random']} datasets of the random split and {fewer['gbt:lineage']} of the lineage split, and the "
            f"lasso on {fewer['lasso:random']} and {fewer['lasso:lineage']}. None of them is the best learner's row:",
            f"In revision 14, the {flagged} learner rows that carry the flag are scored on every genome, and every row read "
            f"here counts every genome of its class.",
            f"The best learner counts fewer than {first} resistant genomes on {table_8[0]} datasets of the random split and "
            f"{table_8[1]} of the lineage split in both releases, and the known-gene rule on {table_8[2]} and {table_8[3]}.",
            f"There are {tb_learners.pop()} of them on every learner row of both splits, and {tb_rule.pop()} on each split of "
            f"the known-gene rule. None of them meets the wording on any row of either release,",
        ]
        text = " ".join(self.text.split())
        for s in sentences:
            with self.subTest(sentence=s[:60]):
                self.assertIn(s, text)


if __name__ == "__main__":
    unittest.main()
