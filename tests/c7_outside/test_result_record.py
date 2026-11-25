"""docs/outside-learners-result.md against the committed tables of reference/outside-learners/.

This is the one test file of the suite that reads committed result tables: the aggregates page_tables.py wrote from
the one reading, by pair, split, measure, species, antibiotic and arm. It fits nothing and reads no genome, report or
fold. Each table on the page must be the one `page_tables.py print` prints from the record, and each figure in the
page's sentences is rebuilt here from the record's files.
"""

import ast
import csv
import hashlib
import re
import unittest

from tests.c7_outside import common

page_tables = common.load("page_tables")
page_faults = common.load("page_faults")

RECORD = common.ROOT / "reference" / "outside-learners"
PAGE = common.ROOT / "docs" / "outside-learners-result.md"
HERE = common.ROOT / "tests" / "c7_outside" / "test_result_record.py"
# Where each printed table stands. `datasets` is printed too, and stays in page-datasets.tsv: the section Datasets
# gives the first run's datasets per species instead.
TABLES = {"Reading": ["conditions", "withdrawn"], "Data": ["checks"], "Count": ["count"],
          "Per arm": ["measures", "folds", "libraries", "binaries", "times"],
          "Beside it": ["no-constant", "no-unconverged"], "Seen after": ["species"]}
ARMS = (("xgboost", "xgboost-gbt"), ("lightgbm", "lightgbm-gbt"), ("elasticnet", "elasticnet-lasso"))
TWINS = ("twin-gbt", "twin-lasso")
d3, of, interval, pval, words = (page_tables.d3, page_tables.of, page_tables.interval, page_tables.pval,
                                 page_tables.count_words)


def rows(name):
    with open(RECORD / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_row(split, pair, measure="balanced_accuracy", file="count.tsv"):
    found = [r for r in rows("page-count.tsv")
             if (r["file"], r["split"], r["pair"], r["measure"]) == (file, split, pair, measure)]
    assert len(found) == 1, (file, split, pair, measure)
    return found[0]


def read(pair, condition):
    found = [r for r in rows("page-conditions.tsv")
             if (r["pair"], r["condition"], r["decides"]) == (pair, condition, "1")]
    assert len(found) == 1, (pair, condition)
    return {"1": "met", "0": "not met"}[found[0]["met"]]


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def blocks(text):
    return [b for b in re.split(r"\n\s*\n", text) if b.strip()]


class Record(unittest.TestCase):
    def test_the_record_holds_the_files_its_checksums_give(self):
        listed = {}
        for line in (RECORD / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            listed[name] = digest
        found = sorted(p.relative_to(RECORD).as_posix() for p in RECORD.rglob("*")
                       if p.is_file() and p.name != "SHA256SUMS")
        self.assertEqual(sorted(listed), found)
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((RECORD / name).read_bytes()).hexdigest(), digest, name)

    def test_the_record_passed_every_gate_and_run_check(self):
        gates = (RECORD / "gates.txt").read_text().splitlines()
        self.assertEqual([line.split(":")[0] for line in gates], [f"{n} {name}" for n, name, _ in page_tables.GATES])
        checks = (RECORD / "checks.txt").read_text().splitlines()
        self.assertEqual(sum(line.startswith("PASS ") for line in checks), 6)
        self.assertFalse([line for line in checks if line.startswith("FAIL")])
        self.assertEqual(checks[-1], "CHECKS 6 of 6")


class ResultPage(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.lines = self.text.split("\n")
        self.sections = sections(self.text)
        self.tables = page_tables.tables(RECORD)

    def test_every_printed_table_stands_whole_in_its_section(self):
        self.assertEqual(sorted(self.tables), sorted([t for names in TABLES.values() for t in names] + ["datasets"]))
        for section, names in TABLES.items():
            for name in names:
                self.assertIn("\n".join(self.tables[name]) + "\n", self.sections[section], name)

    def test_each_printed_line_stands_as_often_as_the_tables_print_it(self):
        on_page = [name for names in TABLES.values() for name in names]
        printed = [line for name in on_page for line in self.tables[name] if not line.startswith("|---")]
        for line in set(printed):
            want = sum(line in self.tables[name] for name in on_page)
            got = self.lines.count(line) if line.startswith("|") else self.text.count(line)
            self.assertEqual(got, want, line)

    def test_the_reading_gives_each_arm_condition_by_condition(self):
        reading = self.sections["Reading"]
        for arm, pair in ARMS:
            c, recall, study = count_row("lineage", pair), count_row("lineage", pair, "recall"), count_row("study", pair)
            agree = (float(study["median_difference"]) > 0) == (float(c["median_difference"]) > 0)
            self.assertEqual(agree, read(pair, "5") == "met", pair)
            self.assertEqual(c["meets_1"] == "1", read(pair, "1") == "met", pair)
            self.assertEqual(c["meets_2"] == "1", read(pair, "2") == "met", pair)
            self.assertIn(
                f"**`{arm}`.** `{pair}` is better on {c['better']} of {c['entered']} lineage datasets, worse on "
                f"{c['worse']} and the same on {c['same']}. Its median difference in balanced accuracy is "
                f"{d3(c['median_difference'])}, with a species interval of "
                f"{interval(c['interval_lo'], c['interval_hi'])} (condition 1 {read(pair, '1')}). The one-sided sign "
                f"test gives p {pval(c['sign_p'])}, and {pval(c['holm_p'])} after Holm (condition 2 "
                f"{read(pair, '2')}). Its recall interval is {interval(recall['interval_lo'], recall['interval_hi'])} "
                f"(condition 4 {read(pair, '4')}). Its study median is {d3(study['median_difference'])} on "
                f"{study['entered']} of {study['datasets']} study datasets, against the lineage median's "
                f"{d3(c['median_difference'])}: the signs {'agree' if agree else 'differ'} (condition 5 "
                f"{read(pair, '5')}).", reading)

    def test_the_reading_gives_the_twin_and_the_result(self):
        reading = self.sections["Reading"]
        for pair in TWINS:
            c = count_row("lineage", pair)
            self.assertEqual((c["meets_1"], c["meets_2"]), ("0", "0"), pair)
            self.assertIn(f"`{pair}` is better on {c['better']} of {c['entered']}, worse on {c['worse']} and the same "
                          f"on {c['same']}, with a median of {d3(c['median_difference'])}, an interval of "
                          f"{interval(c['interval_lo'], c['interval_hi'])} and a Holm p of {pval(c['holm_p'])}.",
                          reading)
        self.assertEqual({read(pair, "3") for _, pair in ARMS}, {"met"})
        self.assertIn("The twin meets condition 1 in neither pair and 2 in neither, so condition 3 is met for every "
                      "arm.", reading)
        whole = [arm for arm, pair in ARMS if all(read(pair, k) == "met" for k in "12345")]
        short = {arm: [k for k in "12345" if read(pair, k) != "met"] for arm, pair in ARMS if arm not in whole}
        self.assertEqual((whole, short), (["lightgbm", "elasticnet"], {"xgboost": ["5"]}))
        self.assertIn("**Result.** `lightgbm` and `elasticnet` meet conditions 1 to 5. `xgboost` meets 1 to 4 and not "
                      "5: its study median has the other sign. Condition 6 is the maintainer's for all three",
                      reading)

    def test_the_clock_is_the_commands_tables_first_start_and_last_end(self):
        steps = [line.split(" | ") for line in self.sections["Method"].splitlines()
                 if re.match(r"\| [a-z]+ \| `python3 -I ", line)]
        self.assertEqual([s[0] for s in steps], ["| check", "| gather", "| count", "| write", "| print"])
        self.assertIn(f"The run was read once, from {steps[0][2]} to {steps[-1][3]} on 25 Nov 2026, with no other "
                      f"option, subset or measure", self.text)
        self.assertEqual(self.text.count("The run was read once"), 1)

    def test_the_commands_table_gives_reading_txts_sha256(self):
        method = self.sections["Method"]
        reading = (RECORD / "reading.txt").read_text().splitlines()
        sums = [line for line in reading if re.match(r"[0-9a-f]{64}  ", line)]
        counted = [line for line in sums if line.split("  ", 1)[1].startswith("count/")]
        self.assertEqual(len(counted), 4)
        self.assertIn(f"`reading.txt` gives the sha256 of the {words(len(sums) - len(counted))} files `gather` "
                      f"copied and of the four files the count wrote:", method)
        for line in counted:
            self.assertIn(f"| `{line.split('  ', 1)[1]}` | `{line[:64]}` |", method)
        self.assertIn(next(line for line in reading if line.startswith("score.py --count: ")).split(": ", 1)[1],
                      method)

    def test_the_unscored_dataset_is_named_in_full(self):
        unscored = [r for r in rows("page-datasets.tsv") if r["difference"] == ""]
        self.assertEqual(len(unscored), 5)
        named = {(r["species"], r["antibiotic"], r["genomes"], r["split"]) for r in unscored}
        self.assertEqual(len(named), 1)
        species, antibiotic, genomes, split = named.pop()
        self.assertEqual(split, "lineage")
        self.assertEqual({(r["draws"], r["arm_draws_scored"], r["default_draws_scored"]) for r in unscored},
                         {("3", "2", "2")})
        self.assertIn(f"One lineage dataset of the 233 is not scored: *{species}* with {antibiotic} ({genomes} "
                      f"genomes). Of its 3 lineage draws, 2 are scored", self.sections["Data"])
        for _, pair in ARMS:
            self.assertEqual(count_row("study", pair)["entered"], "55")
        for pair in [p for _, p in ARMS] + list(TWINS):
            self.assertEqual((count_row("lineage", pair)["entered"], count_row("lineage", pair)["datasets"]),
                             ("232", "233"), pair)

    def test_the_beside_sentences_give_the_count_files(self):
        beside = self.sections["Beside it"]
        whole = count_row("lineage", "elasticnet-lasso")
        for name, column, words_ in (("count-no-constant.tsv", "constant", "hold a draw scored at its training share"),
                                     ("count-no-unconverged.tsv", "unconverged",
                                      "hold an elastic net refit that reached 2,000 epochs")):
            c = count_row("lineage", "elasticnet-lasso", file=name)
            self.assertEqual(int(c["entered"]), int(whole["entered"]) - int(whole[column]), name)
            self.assertIn(f"{whole[column]} {words_}", beside)
            self.assertIn(f"enters {c['entered']} of {c['datasets']}: better on {c['better']}, worse on {c['worse']} "
                          f"and the same on {c['same']}, with a median of {d3(c['median_difference'])} and an "
                          f"interval of {interval(c['interval_lo'], c['interval_hi'])}.", beside)
            for pair in ("xgboost-gbt", "lightgbm-gbt") + TWINS:
                a, b = count_row("lineage", pair), count_row("lineage", pair, file=name)
                self.assertEqual([a[k] for k in ("entered", "better", "worse", "same")],
                                 [b[k] for k in ("entered", "better", "worse", "same")], (name, pair))

    def test_the_seen_after_figures_follow_from_the_tables(self):
        seen = self.sections["Seen after"]
        times = {r["arm"]: r for r in rows("page-times.tsv")}
        self.assertIn(f"Over {words(times['elasticnet']['fits'])} fits an arm, the elastic net's median fit took "
                      f"{d3(times['elasticnet']['wall_s_median'])} wall seconds, against "
                      f"{d3(times['xgboost']['wall_s_median'])} for XGBoost and "
                      f"{d3(times['lightgbm']['wall_s_median'])} for LightGBM.", seen)
        loss, recall = count_row("lineage", "elasticnet-lasso", "log_loss"), count_row("lineage", "elasticnet-lasso",
                                                                                          "recall")
        self.assertIn(f"the elastic net's log loss is better on {loss['better']} of {loss['entered']} datasets and "
                      f"worse on {loss['worse']}, with a median of {d3(loss['median_difference'])} and an interval of "
                      f"{interval(loss['interval_lo'], loss['interval_hi'])}, while its recall is better on "
                      f"{recall['better']} and worse on {recall['worse']}, with a median of "
                      f"{d3(recall['median_difference'])}.", seen)
        xl, ll = count_row("study", "xgboost-gbt", "log_loss"), count_row("study", "lightgbm-gbt", "log_loss")
        xr, lr = count_row("study", "xgboost-gbt", "recall"), count_row("study", "lightgbm-gbt", "recall")
        self.assertTrue(float(xl["median_difference"]) > 0 and float(ll["median_difference"]) > 0)
        self.assertEqual(d3(xr["median_difference"]), d3(lr["median_difference"]))
        self.assertIn(f"XGBoost's median is {d3(xl['median_difference'])} "
                      f"({interval(xl['interval_lo'], xl['interval_hi'])}) and LightGBM's {d3(ll['median_difference'])} "
                      f"({interval(ll['interval_lo'], ll['interval_hi'])}), on {xl['entered']} of {xl['datasets']}, "
                      f"with a recall median of {d3(xr['median_difference'])} beside each.", seen)
        moved = {p: int(count_row("lineage", p)["better"]) + int(count_row("lineage", p)["worse"])
                 for p in [p for _, p in ARMS] + list(TWINS)}
        self.assertIn(f"the twin moves {moved['twin-gbt']} of 232 lineage datasets by more than 0.02 on `gbt` and "
                      f"{moved['twin-lasso']} on `lasso`. The arms move {moved['xgboost-gbt']} (`xgboost-gbt`), "
                      f"{moved['lightgbm-gbt']} (`lightgbm-gbt`) and {moved['elasticnet-lasso']} "
                      f"(`elasticnet-lasso`) of 232.", seen)

    def test_the_run_checks_table_is_checks_txt(self):
        method = self.sections["Method"]
        checks = (RECORD / "checks.txt").read_text().splitlines()
        for line in checks:
            if line.startswith("PASS "):
                name, printed = line[5:].split(": ", 1)
                self.assertEqual(self.lines.count(f"| {name} | {printed} |"), 1, name)
        self.assertIn("`CHECKS 6 of 6`: all six passed.", method)
        values = dict(line.strip().split(": ", 1) for line in checks if line.startswith("  species per "))
        self.assertEqual(values, {"species per export commit": "da70f4b 23", "species per grm sha256": "d3b53997a02c 23",
                                  "species per code sha256": "d3de72fcc10f 23",
                                  "species per export grm sha256": "637b2da2b59b 23",
                                  "species per elastic net pass": "run 23"})
        self.assertIn("The `DONE` lines name export commit `da70f4b`, grm sha256 `d3b53997a02c`, code sha256 "
                      "`d3de72fcc10f` and export grm sha256 `637b2da2b59b` in 23 of 23 species each, and the elastic "
                      "net run in 23 of 23.", method)

    def test_the_fault_table_names_each_page_fault_and_its_test(self):
        method = self.sections["Method"]
        for n, gate, what, _, _, _, test in page_faults.FAULTS:
            self.assertEqual(self.lines.count(f"| {n} | {gate} | {what} | `{test.split('tests.c7_outside.')[1]}` | "
                                              f"caught |"), 1, n)
        self.assertIn(f"All {len(page_faults.FAULTS)} were caught, before any score was opened.", method)

    def test_recall_stands_beside_log_loss_in_every_block(self):
        given = [b for b in blocks(self.text) if "log loss" in b and re.search(r"\d\.\d{3}", b)]
        self.assertGreaterEqual(len(given), 2)
        for block in given:
            self.assertIn("recall", block, block[:80])

    def test_the_method_gives_this_files_own_test_count(self):
        tree = ast.parse(HERE.read_text())
        n = sum(1 for c in tree.body if isinstance(c, ast.ClassDef)
                for f in c.body if isinstance(f, ast.FunctionDef) and f.name.startswith("test_"))
        self.assertIn(f"`tests/c7_outside/test_result_record.py` holds {n} tests of this page against "
                      f"`reference/outside-learners/`", self.sections["Method"])


if __name__ == "__main__":
    unittest.main()
