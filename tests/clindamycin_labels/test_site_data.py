"""The clindamycin kinds the web app bundles (web/src/data/clindamycin-labels.json).

The file must be exactly what scripts/clindamycin-labels/site_data.py writes
from reference/clindamycin-labels/rev14-20261115 and the release's report
checksums beside it, so the site shows the reference files' counts and no
others. The writer is also run on copies of those files with one cell
changed: --check must then fail, and files that no longer agree with each
other must stop the writer. The kind rule is tested here on cases made up
for it; that it gives the page's kind for the six datasets is test_page.py's
check.
"""

import pathlib
import re
import shutil
import tempfile
import unittest

from . import common

S = common.load_site_data()


def append_copy_of(folder, rel, where):
    """Append to a TSV of the copied folder a second copy of the one row that `where` picks."""
    path = pathlib.Path(folder) / rel
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    head = lines[0].split("\t")
    picked = [line for line in lines[1:] if where(dict(zip(head, line.split("\t"))))]
    if len(picked) != 1:
        raise AssertionError(f"{rel}: {len(picked)} rows picked")
    path.write_text("\n".join(lines + picked) + "\n", encoding="utf-8")


def edit(folder, rel, dataset, column, change, where=lambda row: True):
    """Change one column of the rows of `dataset` in a TSV of the copied folder."""
    path = pathlib.Path(folder) / rel
    lines = path.read_text(encoding="utf-8").split("\n")
    head = lines[0].split("\t")
    hits = 0
    for n, line in enumerate(lines[1:], start=1):
        cells = line.split("\t")
        if len(cells) != len(head):
            continue
        row = dict(zip(head, cells))
        if row["dataset"] == dataset and where(row):
            cells[head.index(column)] = change(row[column])
            lines[n] = "\t".join(cells)
            hits += 1
    if hits == 0:
        raise AssertionError(f"{rel}: no row of {dataset} to change")
    path.write_text("\n".join(lines), encoding="utf-8")


class Copy:
    """A copy of the reference folder under its own name, in a temporary
    folder laid out as the repository is, with the release's checksum lines
    beside it."""

    def __enter__(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = pathlib.Path(self.tmp.name) / common.REFERENCE.relative_to(common.ROOT)
        shutil.copytree(common.REFERENCE, self.path)
        sums = pathlib.Path(self.tmp.name) / common.RELEASE_SUMS.relative_to(common.ROOT)
        sums.parent.mkdir(parents=True)
        shutil.copyfile(common.RELEASE_SUMS, sums)
        self.release_sums = sums
        return self.path

    def __exit__(self, *exc):
        self.tmp.cleanup()


SA = "staphylococcus_aureus__clindamycin"
CC = "campylobacter_coli__clindamycin"
CJ = "campylobacter_jejuni__clindamycin"
SP = "streptococcus_pneumoniae__clindamycin"
SG = "streptococcus_agalactiae__clindamycin"
EF = "enterococcus_faecium__clindamycin"


def row_is(split, row, measure):
    return lambda r: (r["split"], r["row"], r["measure"]) == (split, row, measure)


def check(folder):
    """--check on a folder: 1 when the bundled file differs from what the
    folder gives, or when the writer refuses the folder."""
    try:
        return S.main(["site_data.py", str(folder), "--check", str(common.BUNDLED)])
    except SystemExit as stop:
        if not isinstance(stop.code, str):
            raise
        return 1


class BundledKinds(unittest.TestCase):
    def test_bundled_file_is_what_the_reference_folder_gives(self):
        self.assertEqual(common.BUNDLED.read_text(encoding="utf-8"), S.build(common.REFERENCE))
        self.assertEqual(S.main(["site_data.py", str(common.REFERENCE), "--check", str(common.BUNDLED)]), 0)

    def test_a_second_run_writes_the_same_bytes(self):
        self.assertEqual(S.build(common.REFERENCE), S.build(common.REFERENCE))

    def test_check_fails_when_a_count_of_the_folder_changes(self):
        changes = [
            ("results/chosen.tsv", SA, "disagreements_in_study", lambda v: str(int(v) - 1)),
            ("results/chosen.tsv", CC, "samples_in_study", lambda v: str(int(v) + 1)),
            ("results/chosen.tsv", CJ, "chosen_pmid", lambda v: "1" + v),
            ("results/by-study.tsv", SP, "samples", lambda v: str(int(v) + 1)),
        ]
        for rel, dataset, column, change in changes:
            with self.subTest(column=column), Copy() as copy:
                if column == "chosen_pmid":
                    edit(copy, "results/caveat.tsv", dataset, column, change)
                    edit(copy, rel, dataset, column, change)
                elif rel == "results/by-study.tsv":
                    edit(copy, rel, dataset, column, change, where=lambda row: row["study"] == "no-pmid")
                else:
                    edit(copy, rel, dataset, column, change)
                self.assertEqual(check(copy), 1)

    def test_files_that_do_not_agree_stop_the_writer(self):
        refused = [
            # A passing row more in the caveat table than in the score table.
            ("results/caveat.tsv", SA, "rows_passing", lambda v: str(int(v) + 1), lambda row: True, "its dependence part"),
            # The caveat table calls a scored dataset unscored.
            ("results/caveat.tsv", SP, "existence", lambda v: "yes", lambda row: True, "existence part"),
            # Scored without its study, yet a class is left under the minimum.
            ("results/chosen.tsv", CJ, "refit", lambda v: "yes", lambda row: True, "refit 'yes' with 8 resistant"),
            # A verdict the rule does not give.
            ("results/caveat.tsv", SP, "caveat", lambda v: "yes", lambda row: True, "is not the rule's"),
            # The known-gene rule's rows differ between the two splits.
            (f"results/scores/{SA}/caveat.tsv", SA, "total", lambda v: "0.0715", lambda row: row["row"] == "catalog" and row["split"] == "lineage" and row["measure"] == "recall", "differs between the splits"),
            # One passing change turned the other way: the changes go both ways.
            (f"results/scores/{CC}/caveat.tsv", CC, "total", lambda v: v.lstrip("-"), lambda row: row["split"] == "random" and row["row"] == "rf" and row["measure"] == "recall", "go both ways"),
            # Not a count.
            ("results/chosen.tsv", SA, "reported_R_left", lambda v: "325.0", lambda row: True, "is not a count"),
        ]
        for rel, dataset, column, change, where, message in refused:
            with self.subTest(rel=rel, column=column), Copy() as copy:
                edit(copy, rel, dataset, column, change, where)
                with self.assertRaisesRegex(SystemExit, re.escape(message)):
                    S.build(copy)

    def test_a_scored_dataset_where_only_the_known_gene_rules_rows_pass_stops_the_writer(self):
        # SP's known-gene rule recall made to pass on both splits, and no learner row.
        rel = f"results/scores/{SP}/caveat.tsv"
        with Copy() as copy:
            rule_recall = lambda row: row["row"] == "catalog" and row["measure"] == "recall"
            for column, value in (("total", "0.0600"), ("cluster_low", "0.0100"), ("cluster_high", "0.1100"), ("passes", "yes")):
                edit(copy, rel, SP, column, lambda v, value=value: value, rule_recall)
            with self.assertRaisesRegex(SystemExit, re.escape("only the known-gene rule's rows pass")):
                S.build(copy)

    def test_a_known_gene_rule_row_repeated_within_one_split_stops_the_writer(self):
        rel = f"results/scores/{SA}/caveat.tsv"
        with Copy() as copy:
            append_copy_of(copy, rel, row_is("lineage", "catalog", "recall"))
            with self.assertRaisesRegex(SystemExit, re.escape("row lineage catalog recall appears twice")):
                S.build(copy)

    def test_a_passes_cell_the_plans_rule_does_not_give_stops_the_writer(self):
        cases = [
            # A change whose interval holds 0, marked as passing.
            (f"results/scores/{SP}/caveat.tsv", SP, row_is("lineage", "rf", "recall"), "yes"),
            # A change of 0.05 or more with an interval that excludes 0, marked as not passing.
            (f"results/scores/{CC}/caveat.tsv", CC, row_is("lineage", "rf", "recall"), "no"),
        ]
        for rel, dataset, where, value in cases:
            with self.subTest(dataset=dataset), Copy() as copy:
                edit(copy, rel, dataset, "passes", lambda v: value, where)
                with self.assertRaisesRegex(SystemExit, re.escape("the plan's rule gives the other")):
                    S.build(copy)

    def test_the_plans_rule_for_one_row(self):
        self.assertTrue(S.passes_rule(-0.05, -0.1, -0.01, 0.05))
        self.assertTrue(S.passes_rule(0.06, 0.001, 0.2, 0.05))
        self.assertFalse(S.passes_rule(0.0499, 0.01, 0.1, 0.05))
        self.assertFalse(S.passes_rule(0.1, 0.0, 0.2, 0.05))
        self.assertFalse(S.passes_rule(-0.1, -0.2, 0.0, 0.05))

    def test_by_study_parts_that_do_not_add_up_stop_the_writer(self):
        cases = [
            # One more disagreement among the samples with no PMID than the TOTAL row holds.
            (CC, lambda row: row["study"] == "no-pmid", "disagreements", "do not add up to its TOTAL"),
            # One more sample under another PMID.
            (SA, lambda row: row["study"] == "26806258", "samples", "do not add up to its TOTAL"),
        ]
        for dataset, where, column, message in cases:
            with self.subTest(dataset=dataset, column=column), Copy() as copy:
                edit(copy, "results/by-study.tsv", dataset, column, lambda v: str(int(v) + 1), where)
                with self.assertRaisesRegex(SystemExit, re.escape(message)):
                    S.build(copy)

    def test_a_report_sha256_that_is_not_the_releases_line_stops_the_writer(self):
        with Copy() as copy:
            path = copy / "results" / "logs" / "released-reports.sha256"
            text = path.read_text(encoding="utf-8")
            digest = text.split("  ")[0]
            path.write_text(text.replace(digest, digest[:-1] + ("0" if digest[-1] != "0" else "1"), 1), encoding="utf-8")
            with self.assertRaisesRegex(SystemExit, re.escape("is not its line in ../../lineage-control/rev14-274-20261112.reports.sha256")):
                S.build(copy)

    def test_release_checksum_lines_that_lack_a_dataset_stop_the_writer(self):
        with Copy() as copy:
            sums = copy.parent.parent / "lineage-control" / "rev14-274-20261112.reports.sha256"
            lines = sums.read_text(encoding="utf-8").splitlines()
            sums.write_text("\n".join(line for line in lines if EF not in line) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(SystemExit, re.escape(f"has no line for {EF}")):
                S.build(copy)

    def test_a_plan_that_names_no_release_or_threshold_stops_the_writer(self):
        for old in ("Release read: `rev14-274-20261112` (evaluation revision 14).", "is at least 0.05 in absolute value and its 95% interval excludes 0"):
            with self.subTest(old=old), Copy() as copy:
                plan = copy / "PLAN.md"
                text = plan.read_text(encoding="utf-8")
                self.assertIn(old, text)
                plan.write_text(text.replace(old, "(removed)"), encoding="utf-8")
                with self.assertRaises(SystemExit):
                    S.build(copy)

    def test_the_index_line_must_name_the_plans_revision(self):
        with Copy() as copy:
            path = copy / "results" / "count-checks.txt"
            path.write_text(path.read_text(encoding="utf-8").replace("evaluation_revision=14", "evaluation_revision=15"), encoding="utf-8")
            with self.assertRaises(SystemExit):
                S.build(copy)


class KindRule(unittest.TestCase):
    """The rule of site_data.py's docstring, on made-up cases."""

    def test_a_dataset_not_scored_is_one_study_when_nothing_is_left(self):
        self.assertEqual(S.kind_of(False, 0, 0, None), "one_study")

    def test_a_dataset_not_scored_with_samples_left_is_class_balance(self):
        for left in [(8, 469), (38, 18), (0, 3), (24, 0)]:
            self.assertEqual(S.kind_of(False, *left, None), "class_balance", left)

    def test_a_scored_dataset_is_measured_dependence_when_a_test_passes(self):
        self.assertEqual(S.kind_of(True, 60, 248, 1), "measured_dependence")
        self.assertEqual(S.kind_of(True, 60, 248, 33), "measured_dependence")

    def test_a_scored_dataset_with_no_passing_test_is_no_dependence(self):
        self.assertEqual(S.kind_of(True, 41, 392, 0), "no_dependence")

    def test_the_writer_names_exactly_the_four_kinds(self):
        self.assertEqual(S.KINDS, ("measured_dependence", "one_study", "class_balance", "no_dependence"))
        kinds = {d["kind"] for d in common.bundled()["datasets"].values()}
        self.assertEqual(kinds, set(S.KINDS))

    def test_the_no_pmid_sentence_prints_for_the_two_campylobacter_datasets_and_no_other(self):
        # The rule of site_data.py's docstring, on the bundled counts; the web app prints the sentence by the same rule.
        holds = sorted(d for d, entry in common.bundled()["datasets"].items() if S.no_pmid_contrast(entry))
        self.assertEqual(holds, [CC, CJ])

    def test_the_no_pmid_rule_needs_each_of_its_three_parts(self):
        coli = dict(common.bundled()["datasets"][CC])
        self.assertTrue(S.no_pmid_contrast(coli))
        # 1. Most samples carry no PMID.
        self.assertFalse(S.no_pmid_contrast({**coli, "no_pmid_samples": coli["samples"] // 2}))
        # 2. The study disagrees more often than the samples with no PMID.
        self.assertFalse(S.no_pmid_contrast({**coli, "no_pmid_disagreements": coli["no_pmid_samples"]}))
        # 3. Among samples with a PMID, the study's share of the disagreements is close to or under its share of them.
        self.assertFalse(S.no_pmid_contrast({**coli, "other_pmid_disagreements": 0}))
        jejuni = dict(common.bundled()["datasets"][CJ])
        self.assertTrue(S.no_pmid_contrast(jejuni))
        self.assertFalse(S.no_pmid_contrast({**jejuni, "study_disagreements": jejuni["study_disagreements"] + 1, "other_pmid_disagreements": 0}))

    def test_direction_needs_every_passing_change_on_one_side(self):
        self.assertEqual(S.direction_of([0.1, 0.2], "x"), "rises")
        self.assertEqual(S.direction_of([-0.1, -0.2], "x"), "falls")
        self.assertIsNone(S.direction_of([], "x"))
        with self.assertRaises(SystemExit):
            S.direction_of([0.1, -0.1], "x")


if __name__ == "__main__":
    unittest.main()
