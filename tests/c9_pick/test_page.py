"""docs/served-pick-result.md against the committed tables of reference/served-pick/.

This is the one test file of the suite that reads committed result tables:
the aggregates the table program wrote from the one reading, by answer,
split, species, antibiotic and learner. It fits nothing and reads no genome.
"""

import csv
import hashlib
import re
import unittest

from tests.c9_pick import common

page_tables = common.load("page_tables")
planted = common.load("planted_faults")

RECORD = common.ROOT / "reference" / "served-pick"
PAGE = common.ROOT / "docs" / "served-pick-result.md"
L_TSV = common.KIT / "l.tsv"
SECTIONS = ["Reading", "Data", "Count", "Beside it", "Seen after", "Limits", "Method", "Datasets"]
TABLES = {"Reading": ["conditions"], "Count": ["twin", "count"], "Beside it": ["beside", "pooled", "picks", "pairs"],
          "Datasets": ["species", "datasets"]}
COMPARING = ["Reading", "Count", "Beside it", "Seen after", "Limits", "Datasets"]
d3, of, interval, pval = page_tables.d3, page_tables.of, page_tables.interval, page_tables.pval


def rows(name, path=None):
    with open(path or RECORD / name, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_row(split, pair, measure="balanced_accuracy", answer="read"):
    found = [r for r in rows("page-count.tsv")
             if (r["answer"], r["split"], r["pair"], r["measure"]) == (answer, split, pair, measure)]
    assert len(found) == 1, (answer, split, pair, measure)
    return found[0]


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def blocks(text):
    return [b for b in re.split(r"\n\s*\n", text) if b.strip()]


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text()
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_every_printed_table_stands_whole_in_its_section(self):
        tables = page_tables.tables(RECORD)
        self.assertEqual(sorted(tables), sorted(t for names in TABLES.values() for t in names))
        for section, names in TABLES.items():
            for name in names:
                printed = tables[name]
                self.assertGreater(len(printed), 2, name)
                self.assertIn("\n".join(printed) + "\n", self.sections[section], name)
                for line in printed[:1] + printed[2:]:
                    self.assertEqual(self.lines.count(line), 1, f"{name}: {line}")

    def test_the_reading_states_the_count_condition_by_condition(self):
        reading = self.sections["Reading"]
        conds = {r["condition"]: r["read"] for r in rows("page-conditions.tsv")}
        met = [c for c in "12345" if conds[c] == "1"]
        unmet = [c for c in "12345" if conds[c] == "0"]
        self.assertEqual((met, unmet), (["3", "4"], ["1", "2", "5"]))
        self.assertIn("Of the plan's six conditions, 3 and 4 are met and 1, 2 and 5 are not.", reading)
        inner, twin = count_row("lineage", "inner-default"), count_row("lineage", "twin-default")
        self.assertEqual((inner["meets_1"], inner["meets_2"], twin["meets_1"], twin["meets_2"]), ("0",) * 4)
        self.assertIn(f"On the lineage split the inner pick's median difference in balanced accuracy is "
                      f"{d3(inner['median_difference'])}, with a species interval of "
                      f"{interval(inner['interval_lo'], inner['interval_hi'])}.", reading)
        self.assertIn(f"The inner pick is better on {inner['better']} of {inner['entered']} datasets, worse on "
                      f"{inner['worse']} and the same on {inner['same']}. The one-sided sign test gives p "
                      f"{pval(inner['sign_p'])}, and {pval(inner['holm_p'])} after Holm.", reading)
        self.assertIn(f"is better on {twin['better']} of {twin['entered']}, worse on {twin['worse']} and the same on "
                      f"{twin['same']}, with a median of {d3(twin['median_difference'])}, an interval of "
                      f"{interval(twin['interval_lo'], twin['interval_hi'])} and a Holm p of {pval(twin['holm_p'])}.",
                      reading)
        recall = count_row("lineage", "inner-default", "recall")
        self.assertIn(f"The inner pick's recall interval is {interval(recall['interval_lo'], recall['interval_hi'])}.",
                      reading)
        study = count_row("study", "inner-default")
        self.assertNotEqual(float(study["median_difference"]) > 0, float(inner["median_difference"]) > 0)
        self.assertIn(f"The study split's median is {d3(study['median_difference'])} on {study['entered']} datasets, "
                      f"against the lineage split's {d3(inner['median_difference'])}. The signs differ.", reading)

    def test_the_folds_the_inner_pick_serves_stand_in_the_reading(self):
        picks = {(r["split"]): r for r in rows("page-picks.tsv") if r["answer"] == "read" and r["arm"] == "inner"}
        lin, st = picks["lineage"], picks["study"]
        self.assertIn(f"The inner pick served another learner than the default's on {of(lin['differs'], lin['folds'])} "
                      f"lineage folds and on {of(st['differs'], st['folds'])} study folds. On "
                      f"{page_tables.count_words(lin['no_pick'])} lineage folds it had no pick",
                      self.sections["Reading"])

    def test_the_datasets_are_the_first_runs_and_named_in_full(self):
        lrows = rows(None, L_TSV)
        run = [r for r in lrows if r["first_run"] == "1"]
        in_m = [r for r in lrows if r["in_M"] == "1"]
        study = [r for r in run if r["in_S"] == "1"]
        species = {r["species"] for r in run}
        self.assertEqual((len(lrows), len(in_m), len(run), len(species), len(study)), (268, 35, 233, 23, 55))
        self.assertIn(f"The first run: {len(run)} datasets in {len(species)} species, the {len(lrows)} of L, the "
                      f"reports whose lineage split ran, less the {len(in_m)} in M,", self.sections["Data"])
        self.assertIn(f"Of the {len(run)}, {len(study)} in {len({r['species'] for r in study})} species are in S",
                      self.sections["Data"])
        named = {(r["species"], r["antibiotic"], r["split"]) for r in rows("page-datasets.tsv")}
        want = {(r["species"], r["antibiotic"], "lineage") for r in run}
        want |= {(r["species"], r["antibiotic"], "study") for r in study}
        self.assertEqual(named, want)
        datasets = self.sections["Datasets"]
        for sp, ab, split in want:
            self.assertIn(f"| *{sp}* | {ab} | {split} |", datasets)

    def test_ruling_9_stands_where_the_arms_are_compared(self):
        for name in COMPARING:
            body = self.sections[name]
            self.assertRegex(body, r"[Rr]uling 9", name)
            self.assertIn("optimism", body, name)
        self.assertIn("never read as showing that the release's pick is right", self.sections["Limits"])

    def test_recall_stands_beside_log_loss(self):
        for block in blocks(self.text):
            if "log loss" in block:
                self.assertIn("recall", block, block[:80])

    def test_a_split_is_never_called_external_validation(self):
        for sentence in re.split(r"(?<=\.)\s", self.text):
            if "external validation" in sentence:
                self.assertRegex(sentence, r"\b(never|not)\b", sentence)

    def test_the_seen_after_figures_follow_from_the_dataset_table(self):
        seen = self.sections["Seen after"]
        lineage = [r for r in rows("page-datasets.tsv") if r["split"] == "lineage"]
        figures = {}
        for key, col in (("twin", "twin_difference"), ("inner", "difference")):
            scored = [r for r in lineage if r[col] != ""]
            self.assertEqual(len(scored), 232)
            low = min(scored, key=lambda r: float(r[col]))
            high = max(scored, key=lambda r: float(r[col]))
            big = sum(1 for r in scored if abs(float(r[col])) > 0.05)
            figures[key] = (big, len(scored), low, high, col)
        big, n, low, high, col = figures["twin"]
        self.assertIn(f"The twin moves {big} of {n} lineage datasets by more than 0.05 in balanced accuracy, from "
                      f"{d3(low[col])} (*{low['species']}* with {low['antibiotic']}) to {d3(high[col])} "
                      f"(*{high['species']}* with {high['antibiotic']}).", seen)
        big, n, low, high, col = figures["inner"]
        self.assertIn(f"The inner pick moves {big} of {n} by more than 0.05, from {d3(low[col])} (*{low['species']}* "
                      f"with {low['antibiotic']}) to {d3(high[col])} (*{high['species']}* with {high['antibiotic']}).",
                      seen)
        species = rows("page-species.tsv")
        low = min(species, key=lambda r: float(r["lineage_inner_median"]))
        high = max(species, key=lambda r: float(r["lineage_inner_median"]))
        self.assertIn(f"run from {d3(low['lineage_inner_median'])} (*{low['species']}*, {low['datasets']} datasets) to "
                      f"{d3(high['lineage_inner_median'])} (*{high['species']}*, {high['datasets']} datasets)", seen)
        spec = count_row("lineage", "inner-default", "specificity")
        self.assertIn(f"specificity is better on {spec['better']} of {spec['entered']} datasets and worse on "
                      f"{spec['worse']}, with a median of {d3(spec['median_difference'])} and an interval of "
                      f"{interval(spec['interval_lo'], spec['interval_hi'])}.", seen)

    def test_the_checks_and_the_reading_stand_as_their_files_give(self):
        checks = (RECORD / "checks.txt").read_text()
        method, limits = self.sections["Method"], self.sections["Limits"]
        for in_file, on_page, where in (
            ("CHECKS 6 of 6", "All six passed.", method),
            ("1683 of 1683 lines hold", "1,683 of 1,683 lines", method),
            ("1453 of 1453 results found, 0 more", "1,453 of 1,453 results found and none more", method),
            ("0 genome ids in 2671 files, against 23757 release ids of 23 species",
             "0 in 2,671 files, against the 23,757 genome ids of the release's 23 species lists", method),
            ("IDCOUNT ids=0 files=488 known=23757", "23,757 genome ids over the reading and those tables found 0 in "
             "488 files", method),
            ("species per commit: 0fd05e9 22 (the rest); dd17f98 1 (klebsiella_pneumoniae)",
             "22 species ran on the binary built from main at `0fd05e9`. *Klebsiella pneumoniae*", limits),
            ("32 of 32 first-run seed-42 results equal the second run's with seconds_inner left out",
             "32 seed-42 results equal the second run's in 32 of 32, with the timing field `seconds_inner` left out",
             limits),
            ("'invalid lineage outer evaluation: fold 4 validation lacks one or both classes'",
             "`invalid lineage outer evaluation: fold 4 validation lacks one or both classes`", self.sections["Data"]),
        ):
            self.assertIn(in_file, checks)
            self.assertIn(on_page, where)
        reading = (RECORD / "reading.txt").read_text().splitlines()
        self.assertEqual(sum(1 for line in reading if line == "  printed: count: 809 dataset rows from 23 files"), 4)
        self.assertEqual(3 * 233 + 2 * 55, 809)
        self.assertIn("Each count printed `count: 809 dataset rows from 23 files`", method)
        started = reading[0].split(", ")[1].split(".")[0]
        self.assertEqual((started, reading[-1]), ("2026-11-23T06:43:51Z", "ended 2026-11-23T06:44:30Z"))
        self.assertIn("The count was read once, from 06:43:51Z to 06:44:30Z on 23 Nov 2026", self.text)
        gates = (RECORD / "gates.txt").read_text().splitlines()
        self.assertEqual(len(gates), len(page_tables.GATES))
        for n, line in enumerate(gates, 1):
            self.assertTrue(line.startswith(f"{n} "))
            name, words = line.split(" ", 1)[1].split(": ", 1)
            self.assertEqual(self.lines.count(f"| {n}. {name} | {words} |"), 1, line)
        self.assertIn(f"All {len(gates)} passed", method)

    def test_the_fault_table_names_each_planted_faults_test(self):
        mine = [f for f in planted.FAULTS if f[2] == "scripts/c9-pick/page_tables.py"]
        self.assertEqual([f[0] for f in mine], ["7", "8", "8b", "9", "10", "11", "12"])
        method = self.sections["Method"]
        for fault, words, _, _, _, _, test in mine:
            self.assertIn(f"| {fault}. {words[0].upper()}{words[1:]} | `{test.split('test_page_tables.')[1]}` |",
                          method)
        self.assertIn(f"All {len(mine)} were caught.", method)

    def test_the_page_gives_the_suites_own_test_counts(self):
        loader = unittest.TestLoader()
        here = loader.loadTestsFromName("tests.c9_pick.test_page").countTestCases()
        tables = loader.loadTestsFromName("tests.c9_pick.test_page_tables").countTestCases()
        self.assertIn(f"`tests/c9_pick/test_page_tables.py` holds {tables} SYNTHETIC tests", self.sections["Method"])
        self.assertIn(f"`tests/c9_pick/test_page.py` holds {here} tests of this page", self.sections["Method"])

    def test_no_mark_of_the_skeleton_is_left(self):
        self.assertNotIn("[AFTER THE", self.text)
        self.assertNotRegex(self.text, r"\{[a-z_0-9]+\}")
        self.assertNotIn("\u2014", self.text)

    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = sorted(
            str(p.relative_to(RECORD)) for p in RECORD.rglob("*") if p.is_file() and p.name != "SHA256SUMS"
        )
        self.assertEqual(sorted(listed), found)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)


if __name__ == "__main__":
    unittest.main()
