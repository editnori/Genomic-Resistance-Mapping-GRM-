"""The committed result of the learner cross-check (reference/learner-crosscheck): its shape and its counts.

The suite fits nothing and needs no reference library. It reads the two tables with the kit's own table.py, the panel
file with the kit's own decoder, and requires that the tables agree with each other, with the panel file, with
totals.json and detail.json, and with the counts measured on 14 November 2026, which are written out below. A new run
that gives other counts fails here until those are changed on purpose."""

import collections
import hashlib
import re
import unittest

from . import kit

table, panels = kit.table, kit.panels

# The 16 panels of tests/data/lasso-stopping-fixtures.txt, in file order: samples, features, training samples.
PANELS = {
    "campylobacter_jejuni__ciprofloxacin_lineage_fold0": (508, 128, 397),
    "campylobacter_jejuni__ciprofloxacin_random_fold0": (508, 128, 406),
    "constant_balanced": (80, 4, 60),
    "constant_imbalanced": (80, 4, 60),
    "enterobacter_cloacae__tobramycin_lineage_fold0": (99, 128, 69),
    "enterobacter_cloacae__tobramycin_random_fold0": (99, 128, 79),
    "imbalanced_5pct": (200, 16, 150),
    "marker_dup_constant": (120, 12, 90),
    "one_class_train": (60, 4, 40),
    "planted_logistic": (300, 40, 225),
    "staphylococcus_haemolyticus__fusidic_acid_lineage_fold0": (114, 128, 90),
    "staphylococcus_haemolyticus__fusidic_acid_random_fold0": (114, 128, 91),
    "streptococcus_agalactiae__clindamycin_lineage_fold0": (141, 128, 101),
    "streptococcus_agalactiae__clindamycin_random_fold0": (141, 128, 112),
    "tied_counts": (100, 3, 86),
    "xor_balanced": (80, 2, 80),
}

# Measured on 14 November 2026 at c435000: (learner, check) -> units compared, matched, with a shown cause, without
# one; then panels matched, with a shown cause, without one, not run.
MEASURED = {
    ("scm", "fits"): (224, 194, 30, 0, 12, 4, 0, 0),
    ("cart", "nodes"): (400, 395, 5, 0, 12, 4, 0, 0),
    ("rf", "trees"): (3200, 2774, 426, 0, 3, 13, 0, 0),
    ("gbt", "tree by tree"): (5200, 4766, 434, 0, 7, 9, 0, 0),
    ("gbt", "end to end"): (15, 11, 4, 0, 11, 4, 0, 1),
}

# The causes PLAN.md allows a difference, as compare.py writes them, and the rows each has in the committed table.
ROUNDING = "float64 rounding broke an exact tie: GRM's own gains for the two features differ"
CAUSES = {
    ("scm", "fits", "the library's stump at a feature's largest value, which rules out every remaining sample or none, has the higher utility"): 30,
    ("cart", "nodes", ROUNDING): 5,
    ("rf", "trees", ROUNDING): 507,
    ("gbt", "tree by tree", "equal split"): 64,
    ("gbt", "tree by tree", "float32 rounding"): 2,
    ("gbt", "tree by tree", "leaf set to 0"): 368,
    ("gbt", "end to end", "a tree differs (equal split)"): 1,
    ("gbt", "end to end", "a tree differs (float32 rounding)"): 1,
    ("gbt", "end to end", "no tree kept"): 2,
}

# The record added after the plan on 15 November 2026 (end-to-end-trees.tsv; PLAN-note-2026-11-15.md): each tree that
# both end-to-end fits keep and that is not identical, as rows per panel and cause. "" is a tree with no cause shown.
FOLLOWS = "an earlier tree moved the margins"
TREES = {
    ("campylobacter_jejuni__ciprofloxacin_lineage_fold0", "equal split"): 5,
    ("campylobacter_jejuni__ciprofloxacin_random_fold0", "equal split"): 5,
    ("enterobacter_cloacae__tobramycin_random_fold0", "equal split"): 2,
    ("planted_logistic", ""): 1,
    ("planted_logistic", FOLLOWS): 9,
    ("staphylococcus_haemolyticus__fusidic_acid_random_fold0", "float32 rounding"): 1,
    ("staphylococcus_haemolyticus__fusidic_acid_random_fold0", FOLLOWS): 2,
    ("tied_counts", "float32 rounding"): 11,
    ("tied_counts", ""): 2,
    ("tied_counts", FOLLOWS): 364,
}

# What a test of the run record asks for when the result no longer describes the code.
RERUN = "rerun scripts/learner-crosscheck/run.sh ENV OUT and replace the committed result with its files (scripts/learner-crosscheck/README.md)"

UNIT = {
    ("scm", "fits"): r"(conjunction|disjunction) p=(0\.1|0\.316|0\.562|1|1\.778|3\.162|10)",
    ("cart", "nodes"): r"node (0|[1-9]\d*)",
    ("rf", "trees"): r"tree (\d|[1-9]\d|1\d\d)",
    ("gbt", "tree by tree"): r"round (\d|[1-9]\d|[1-3]\d\d)",
    ("gbt", "end to end"): r"fit",
}


class Shape(unittest.TestCase):
    def setUp(self):
        self.summary, self.differences = kit.summary(), kit.differences()

    def test_the_tables_have_no_problem(self):
        fixture = [panel["name"] for panel in panels.read(kit.FIXTURES)]
        self.assertEqual(table.problems(self.summary, self.differences, fixture), [])

    def test_one_row_per_panel_and_check(self):
        self.assertEqual(len(self.summary), 16 * 5)
        self.assertEqual([row["panel"] for row in self.summary[::5]], list(PANELS))
        for start in range(0, 80, 5):
            self.assertEqual([(row["learner"], row["check"], row["units"]) for row in self.summary[start : start + 5]], list(table.CHECKS))

    def test_the_panels_are_the_fixture_file_s(self):
        decoded = {panel["name"]: (panel["n"], len(panel["rows"]), panels.count(panel["fit"])) for panel in panels.read(kit.FIXTURES)}
        self.assertEqual(decoded, PANELS)
        self.assertEqual(list(decoded), list(PANELS))
        for row in self.summary:
            self.assertEqual((row["samples"], row["features"], row["fit_samples"]), PANELS[row["panel"]], row["panel"])

    def test_a_training_mask_lies_within_the_samples_and_a_row_within_the_panel(self):
        for panel in panels.read(kit.FIXTURES):
            everyone = (1 << panel["n"]) - 1
            for bits in [panel["y"], panel["fit"], *panel["rows"]]:
                self.assertEqual(bits & ~everyone, 0, panel["name"])

    def test_every_difference_names_a_unit_of_its_check(self):
        for row in self.differences:
            key = (row["learner"], row["check"])
            self.assertRegex(row["unit"], f"^{UNIT[key]}$", row)
            self.assertTrue(row["what"] and row["evidence"], row)

    def test_a_unit_has_one_row_where_the_check_writes_one(self):
        seen = collections.Counter((row["panel"], row["learner"], row["check"], row["unit"]) for row in self.differences)
        for (panel, learner, check, unit), rows in seen.items():
            # a forest tree has one row for each of its nodes that differs
            if (learner, check) != ("rf", "trees"):
                self.assertEqual(rows, 1, (panel, learner, check, unit))

    def test_the_rows_are_in_panel_order(self):
        order = [row["panel"] for row in self.differences]
        self.assertEqual(sorted(order, key=list(PANELS).index), order)


class Counts(unittest.TestCase):
    def setUp(self):
        self.summary, self.differences, self.detail = kit.summary(), kit.differences(), kit.load("detail.json")

    def rows(self, panel, learner, check):
        return [row for row in self.differences if (row["panel"], row["learner"], row["check"]) == (panel, learner, check)]

    def test_totals_json_is_what_the_summary_adds_up_to(self):
        self.assertEqual(kit.load("totals.json"), table.totals(self.summary))

    def test_the_totals_are_the_measured_ones(self):
        got = {}
        for row in table.totals(self.summary):
            by = row["panels"]
            got[(row["learner"], row["check"])] = (
                row["compared"], row["units_matched"], row["units_explained"], row["units_unexplained"],
                by["matched"], by["explained"], by["unexplained"], by["not applicable"],
            )  # fmt: skip
        self.assertEqual(got, MEASURED)

    def test_no_difference_lacks_a_cause(self):
        self.assertEqual([row for row in self.differences if not row["cause"]], [])
        self.assertEqual([row for row in self.summary if row["unexplained"] or row["outcome"] == "unexplained"], [])

    def test_the_causes_are_the_plan_s_and_their_rows_the_measured_ones(self):
        self.assertEqual(collections.Counter((row["learner"], row["check"], row["cause"]) for row in self.differences), CAUSES)
        self.assertEqual(len(self.differences), sum(CAUSES.values()))

    def test_the_units_with_a_difference_are_the_units_the_rows_name(self):
        named = collections.defaultdict(set)
        for row in self.differences:
            named[(row["panel"], row["learner"], row["check"])].add(row["unit"])
        for row in self.summary:
            key = (row["panel"], row["learner"], row["check"])
            with self.subTest(key=key):
                self.assertEqual(row["explained"] + row["unexplained"], len(named[key]))
                self.assertEqual(row["outcome"] == "matched", row["matched"] == row["compared"] and row["outcome"] != "not applicable")

    def test_units_compared(self):
        for row in self.summary:
            key, detail = (row["learner"], row["check"]), self.detail[row["panel"]]
            with self.subTest(panel=row["panel"], check=key):
                if key == ("scm", "fits"):
                    self.assertEqual(row["compared"], 14)
                elif key == ("cart", "nodes"):
                    self.assertEqual(row["compared"], detail["cart nodes"]["nodes_searched"])
                elif key == ("rf", "trees"):
                    self.assertEqual(row["compared"], 200)
                elif key == ("gbt", "tree by tree"):
                    self.assertEqual(row["compared"], detail["gbt tree by tree"]["trees_kept"])
                    self.assertLessEqual(row["compared"], detail["gbt tree by tree"]["rounds_grown"])
                    self.assertEqual(detail["gbt tree by tree"]["rounds_grown"], 400)
                else:
                    self.assertEqual(row["compared"], 0 if row["outcome"] == "not applicable" else 1)

    def test_the_set_covering_machine_s_detail(self):
        for name in PANELS:
            scm, rows = self.detail[name]["scm fits"], self.rows(name, "scm", "fits")
            with self.subTest(panel=name):
                # the independent implementation gave every fit; a fit differs only by a dominated step
                self.assertEqual(scm["independent_same"], 14)
                self.assertEqual(scm["steps_other"], 0)
                self.assertEqual(scm["steps_certified"] + scm["steps_dominated"], scm["steps_read_against_library"])
                self.assertLessEqual(scm["steps_certified_same_effect_only"], scm["steps_certified"])
                self.assertEqual(scm["library_free_fits"] + scm["fits_library_refuses"], 14)
                self.assertEqual(scm["steps_read_against_library"], 0 if scm["fits_library_refuses"] else scm["rules"])
                self.assertLessEqual(len(rows), scm["steps_dominated"])
                self.assertEqual(bool(rows), bool(scm["steps_dominated"]))
        # the library refuses the one panel whose training samples have one class
        self.assertEqual([name for name in PANELS if self.detail[name]["scm fits"]["fits_library_refuses"]], ["one_class_train"])
        self.assertEqual({row["what"] for row in self.differences if row["learner"] == "scm"}, {"steps [1] are not the library's best rule"})

    def test_the_decision_tree_s_detail(self):
        for name in PANELS:
            cart, rows = self.detail[name]["cart nodes"], self.rows(name, "cart", "nodes")
            with self.subTest(panel=name):
                agreed = cart["nodes_both_leaf"] + cart["nodes_same_feature"] + cart["nodes_other_member_of_a_tie"]
                self.assertEqual(agreed + len(rows), cart["nodes_searched"])
                self.assertLessEqual(cart["nodes_searched"], cart["nodes"])
                self.assertEqual((cart["stop_rule_failures"], cart["node_record_failures"], cart["free_tree_parts_elsewhere"]), (0, 0, 0))
                self.assertLessEqual(cart["largest_score_gap"], 1e-12)
                self.assertEqual(cart["free_tree_same"], cart["free_tree_parts_at_exact_ties"] == 0)

    def test_the_forest_s_detail(self):
        for name in PANELS:
            forest, rows = self.detail[name]["rf trees"], self.rows(name, "rf", "trees")
            row = next(r for r in self.summary if (r["panel"], r["learner"]) == (name, "rf"))
            with self.subTest(panel=name):
                agreed = forest["nodes_both_leaf"] + forest["nodes_same_feature"] + forest["nodes_other_member_of_a_tie"]
                self.assertEqual(agreed, forest["nodes_agree"])
                # one row per node that differs; the trees that hold them are the trees with a shown cause
                self.assertEqual(forest["nodes_searched"] - forest["nodes_agree"], len(rows))
                self.assertEqual(len({r["unit"] for r in rows}), row["explained"])
                self.assertLessEqual(forest["largest_score_gap"], 1e-12)
                self.assertEqual(forest["subsample_size"], int(PANELS[name][2] * 0.632 + 0.5))

    def test_the_boosted_trees_detail(self):
        for name in PANELS:
            by_tree, end = self.detail[name]["gbt tree by tree"], self.detail[name]["gbt end to end"]
            rows = self.rows(name, "gbt", "tree by tree")
            with self.subTest(panel=name):
                self.assertEqual(collections.Counter(row["cause"] for row in rows), by_tree.get("causes", {}))
                self.assertLessEqual(by_tree["largest_margin_gap"], 1e-9)
                fit = self.rows(name, "gbt", "end to end")
                if "reason" in end:
                    self.assertEqual(fit, [])
                    continue
                self.assertEqual(end["fitting_samples"] + end["validation_samples"], PANELS[name][2])
                same_count = end["grm_trees"] == end["xgboost_trees"]
                # matched end to end is the plan's rule: the same number of trees and probabilities within 1e-4
                self.assertEqual(fit == [], same_count and end["largest_probability_gap"] <= 1e-4)
        self.assertEqual([name for name in PANELS if "reason" in self.detail[name]["gbt end to end"]], ["one_class_train"])

    def test_what_a_match_holds_is_the_measured(self):
        # the README's table of what a match holds although the two sides differ, as counted on 15 November 2026
        self.assertEqual(kit.detail_sum("cart nodes", "nodes_other_member_of_a_tie"), 55)
        self.assertEqual((kit.detail_sum("rf trees", "nodes_other_member_of_a_tie"), kit.detail_sum("rf trees", "nodes_agree")), (6871, 108575))
        matched = [row["panel"] for row in self.summary if row["check"] == "end to end" and row["outcome"] == "matched"]
        record = [row for row in kit.trees() if row["panel"] in matched]
        explored = kit.load("exploratory.json")["end_to_end_trees"]
        self.assertEqual((len(matched), len({row["panel"] for row in record}), len(record)), (11, 4, 20))
        self.assertEqual(sum(explored[panel]["trees_both_keep"] for panel in matched), 612)

    def test_a_panel_with_no_tree_kept_compares_no_tree(self):
        empty = [row["panel"] for row in self.summary if row["check"] == "tree by tree" and row["compared"] == 0]
        self.assertEqual(empty, ["constant_balanced", "constant_imbalanced", "xor_balanced"])


class Record(unittest.TestCase):
    def test_versions_are_the_pinned_ones(self):
        versions = kit.load("versions.json")
        pins = kit.module("check_env").pins()
        pinned = {name: version for file in ("requirements.txt", "requirements-scm.txt") for name, version, _ in pins[file]}
        self.assertEqual(versions["libraries"], pinned)
        self.assertTrue(versions["python"].startswith("3.12."))
        self.assertIs(versions["xgboost_build"]["USE_CUDA"], False)

    def test_the_environment_held_the_pinned_versions(self):
        pins = kit.module("check_env").pins()
        want = [f"{name}\t{version}\t{version}" for file in pins.values() for name, version, _ in file]
        lines = (kit.REFERENCE / "environment.tsv").read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[: len(want)], want)
        # then the environment's pip, which no file pins, and the pyscm-ml wheel built from its pinned archive:
        # the wheel's sha256 and the compiler its extension module names
        version = dict((name, version) for file in pins.values() for name, version, _ in file)["pyscm-ml"]
        self.assertEqual(len(lines), len(want) + 3)
        self.assertRegex(lines[-3], r"^pip\tnot pinned\t\d+(\.\d+)+$")
        wheel = rf"^pyscm-ml wheel\tbuilt from the pinned archive\tpyscm_ml-{re.escape(version)}-cp312-cp312-linux_x86_64\.whl sha256 [0-9a-f]{{64}}$"
        self.assertRegex(lines[-2], wheel)
        self.assertRegex(lines[-1], r"^pyscm-ml compiler\tbuilt from the pinned archive\tGCC: \S.*$")

    def test_the_panel_file_is_the_one_that_was_read(self):
        versions = kit.load("versions.json")
        self.assertEqual(versions["fixtures_sha256"], hashlib.sha256(kit.FIXTURES.read_bytes()).hexdigest())
        self.assertEqual(list(versions["grm_fits_sha256"]), sorted(PANELS))
        for digest in versions["grm_fits_sha256"].values():
            self.assertRegex(digest, "^[0-9a-f]{64}$")

    def test_the_run_record(self):
        # run.json is the run of 14 November 2026; rerun.json the run of 15 November that wrote the record added after the plan
        versions = kit.load("versions.json")
        for name, day in (("run.json", "2026-11-14"), ("rerun.json", "2026-11-15")):
            with self.subTest(name=name):
                run = kit.load(name)
                self.assertEqual(list(run), ["head", "src_tree", "src_and_lockfile_as_committed", "rustc", "harness_identical_at_1_and_8_threads", "started", "ended"])
                self.assertRegex(run["head"], "^[0-9a-f]{40}$")
                self.assertRegex(run["src_tree"], "^[0-9a-f]{40}$")
                self.assertIs(run["src_and_lockfile_as_committed"], True)
                self.assertIs(run["harness_identical_at_1_and_8_threads"], True)
                self.assertTrue(run["rustc"].startswith("rustc 1.98.1 "), run["rustc"])
                stamp = rf"^{day}T\d\d:\d\d:\d\dZ$"
                self.assertRegex(run["started"], stamp)
                self.assertRegex(run["ended"], stamp)
                self.assertLess(run["started"], run["ended"])
        self.assertEqual(versions["evaluation_revision"], 14)

    def test_the_result_is_of_the_src_tree_under_test(self):
        # run.json names the tree of src/ that run.sh built. A change to a learner gives another tree, and the
        # committed tables then describe code that is no longer here.
        tree = kit.git("rev-parse", "HEAD:src")
        for name in ("run.json", "rerun.json"):
            run = kit.load(name)
            self.assertEqual(
                run["src_tree"], tree,
                f"{name} was written for the src/ tree {run['src_tree']}, and src/ at the head is the tree {tree}: {RERUN}",
            )  # fmt: skip

    def test_src_has_no_change_that_no_commit_holds(self):
        changed = kit.git("status", "--porcelain", "--untracked-files=all", "--", "src")
        self.assertEqual(changed, "", f"src/ has changes outside a commit, so no run.json can name its tree: commit them, then {RERUN}")

    def test_the_revision_measured_is_not_ahead_of_the_code(self):
        source = (kit.ROOT / "src" / "analysis.rs").read_text(encoding="utf-8")
        found = re.findall(r"^pub const EVALUATION_REVISION: u32 = (\d+);$", source, re.M)
        self.assertEqual(len(found), 1)
        self.assertLessEqual(kit.load("versions.json")["evaluation_revision"], int(found[0]))

    def test_the_clarifications_name_the_plan_they_clarify(self):
        plan = hashlib.sha256((kit.REFERENCE / "PLAN.md").read_bytes()).hexdigest()
        self.assertIn(f"`PLAN.md` (sha256 `{plan}`) is unchanged", (kit.REFERENCE / "PLAN-clarifications.md").read_text(encoding="utf-8"))

    def test_the_exploratory_counts_are_labelled_and_add_up(self):
        e = kit.load("exploratory.json")
        self.assertTrue(e["label"].startswith("exploratory: written after the first comparison was read"))
        for learner, causes in (("cart", ("cart", "nodes", ROUNDING)), ("rf", ("rf", "trees", ROUNDING))):
            ties = e["rounding_ties"][learner]
            self.assertEqual(ties["nodes"], CAUSES[causes])
            self.assertEqual(ties["features_present_where_the_other_is_absent"] + ties["other_exact_ties"], ties["nodes"])
        splits = e["equal_splits_total"]
        self.assertEqual(splits["same_way_round"] + splits["other_way_round"], CAUSES[("gbt", "tree by tree", "equal split")])
        self.assertEqual(e["end_to_end_trees_total"]["trees_both_keep"], kit.detail_sum("gbt end to end", "grm_trees"))
        self.assertEqual(sorted(e["end_to_end_trees"]), sorted(set(PANELS) - {"one_class_train"}))

    def test_the_note_beside_the_plan_names_the_plan_and_the_record(self):
        note = (kit.REFERENCE / "PLAN-note-2026-11-15.md").read_text(encoding="utf-8")
        for name in ("PLAN.md", "PLAN-clarifications.md"):
            digest = hashlib.sha256((kit.REFERENCE / name).read_bytes()).hexdigest()
            self.assertIn(f"`{name}` (sha256 `{digest}`)", note)
        for name in ("end-to-end-trees.tsv", "rerun.json", "detail.json", "environment.tsv"):
            self.assertIn(f"`{name}`", note)

    def test_the_lasso_suite_s_lines(self):
        lines = (kit.REFERENCE / "lasso-suite.txt").read_text(encoding="utf-8").splitlines()
        self.assertEqual(lines[0], "running 29 tests")
        tests = [line for line in lines if line.startswith("test learn::lasso::")]
        self.assertEqual(collections.Counter(line.split(" ... ", 1)[1].split(",")[0] for line in tests), {"ok": 28, "ignored": 1})
        self.assertEqual([line.split(" ... ")[0] for line in tests if " ... ignored" in line], ["test learn::lasso::tests::mtb_delamanid_fold0_reaches_the_certified_optima"])
        self.assertTrue(lines[-1].startswith("test result: ok. 28 passed; 0 failed; 1 ignored; 0 measured; "), lines[-1])


class TreeRecord(unittest.TestCase):
    """end-to-end-trees.tsv, added after the plan: a cause for each end-to-end tree that differs."""

    def setUp(self):
        self.rows, self.detail = kit.trees(), kit.load("detail.json")
        self.constants = kit.constants("check_gbt")

    def test_the_trees_are_the_measured_ones(self):
        self.assertEqual(collections.Counter((row["panel"], row["cause"]) for row in self.rows), TREES)
        self.assertEqual(self.constants["FOLLOWS"], FOLLOWS)

    def test_each_tree_both_keep_that_is_not_identical_has_one_row(self):
        explored = kit.load("exploratory.json")["end_to_end_trees"]
        count = collections.Counter(row["panel"] for row in self.rows)
        self.assertLessEqual(set(count), set(explored))
        for panel, trees in explored.items():
            self.assertEqual(count[panel], trees["trees_both_keep"] - trees["identical"], panel)
        names = list(PANELS)
        self.assertEqual([row["panel"] for row in self.rows], sorted((row["panel"] for row in self.rows), key=names.index))
        for panel in count:
            rounds = [int(row["unit"].split()[1]) for row in self.rows if row["panel"] == panel]
            self.assertEqual(rounds, sorted(set(rounds)), panel)
            self.assertLess(rounds[-1], min(self.detail[panel]["gbt end to end"]["grm_trees"], self.detail[panel]["gbt end to end"]["xgboost_trees"]))
        for row in self.rows:
            self.assertRegex(row["unit"], "^" + UNIT[("gbt", "tree by tree")] + "$")

    def test_a_tree_follows_an_earlier_one_exactly_where_the_margins_had_moved(self):
        tolerance = self.constants["MARGIN_TOLERANCE"]
        for row in self.rows:
            with self.subTest(panel=row["panel"], unit=row["unit"]):
                self.assertEqual(row["cause"] == FOLLOWS, row["fitting_margin_gap"] > tolerance)
                split = r"^at root[.a-z]* GRM splits on feature \d+ and XGBoost on feature \d+$"
                if row["cause"] == FOLLOWS:
                    self.assertEqual(row["evidence"], f"before this round the two margins on the fitting samples differ by up to {row['fitting_margin_gap']!r}")
                    self.assertRegex(row["what"], split + r"|^[1-9]\d* leaf values? differs?, the first at root[.a-z]*$")
                else:
                    self.assertRegex(row["what"], split)
                if row["unit"] == "round 0":
                    self.assertEqual(row["fitting_margin_gap"], 0.0)

    def test_the_first_tree_that_differs_has_the_cause_the_plan_gave_the_panel(self):
        for row in kit.differences():
            if row["check"] == "end to end" and row["cause"].startswith("a tree differs ("):
                first = next(tree for tree in self.rows if tree["panel"] == row["panel"])
                self.assertEqual(f"a tree differs ({first['cause']})", row["cause"])
                self.assertEqual(int(first["unit"].split()[1]), self.detail[row["panel"]]["gbt end to end"]["first_tree_not_identical"])
                self.assertTrue(row["what"].endswith(f"; the first tree that is not identical is {first['unit']}: {first['what']}"), row["what"])
                self.assertEqual(first["evidence"], row["evidence"])

    def test_the_probabilities_on_the_fitting_samples_are_no_further_apart_than_on_all(self):
        for panel, checks in self.detail.items():
            end = checks["gbt end to end"]
            if "largest_probability_gap" in end:
                self.assertLessEqual(end["largest_probability_gap_on_the_fitting_samples"], end["largest_probability_gap"], panel)


if __name__ == "__main__":
    unittest.main()
