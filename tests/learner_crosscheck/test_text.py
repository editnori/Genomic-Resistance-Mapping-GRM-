"""What is written about the learner cross-check, held to the committed result.

Three texts state the result: the Methods view of the site (web/src/views/MethodsView.tsx, the part "Learner
cross-check"), the README of reference/learner-crosscheck and the README of the kit. The whole Methods part is held
here, sentence by sentence and in order: each sentence that carries a number is rebuilt from the result files, and the
others are written out, so no sentence can be added, dropped or changed unseen. Each table of the result's README is
rebuilt line for line, and its sha256 lines are checked against the files."""

import hashlib
import html
import re
import sys
import unittest

from . import kit

METHODS = kit.ROOT / "web" / "src" / "views" / "MethodsView.tsx"
SCOPE = "says nothing about clinical performance"
# The commit that first holds the kit; the run of record (run.json) used a working tree on c435000, which holds none of it.
KIT_COMMIT = "8a74f167ed60bbebf81eeb525626bc130ef79455"
MONTHS = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")


def day(stamp):
    year, month, number = stamp[:10].split("-")
    return f"{int(number)} {MONTHS[int(month) - 1]} {year}"


def methods_part():
    """The text of the Methods part, its tags taken out and its white space joined."""
    source = METHODS.read_text(encoding="utf-8")
    found = re.findall(r'<Part title="Learner cross-check">(.*?)</Part>', source, re.S)
    if len(found) != 1:
        raise AssertionError(f"{len(found)} parts named Learner cross-check")
    return " ".join(html.unescape(re.sub(r"<[^<>]*>", " ", found[0])).split())


def tree_differs():
    """The two panels whose end-to-end difference is a tree that differs, in panel order."""
    return [row["panel"] for row in kit.differences() if row["check"] == "end to end" and row["cause"].startswith("a tree differs (")]


def sentences():
    """Every sentence of the Methods part, in order: those with a number from the result files, the others as written."""
    totals = {(row["learner"], row["check"]): row for row in kit.load("totals.json")}
    versions, run = kit.load("versions.json"), kit.load("run.json")
    libraries = versions["libraries"]
    scm, cart, forest = totals[("scm", "fits")], totals[("cart", "nodes")], totals[("rf", "trees")]
    by_tree, end = totals[("gbt", "tree by tree")], totals[("gbt", "end to end")]
    causes = {}
    for row in kit.differences():
        causes[(row["check"], row["cause"])] = causes.get((row["check"], row["cause"]), 0) + 1
    gaps = [panel["gbt end to end"].get("largest_probability_gap", 0) for panel in kit.load("detail.json").values()]
    lasso = (kit.REFERENCE / "lasso-suite.txt").read_text(encoding="utf-8").splitlines()[-1]
    passed = re.fullmatch(r"test result: ok\. (\d+) passed; 0 failed; 1 ignored; .*", lasso).group(1)
    no_tree = causes[("end to end", "no tree kept")]
    # the record added on 15 November 2026: each end-to-end tree that differs, with a cause of its own
    record, constants, detail = kit.trees(), kit.constants("check_gbt"), kit.load("detail.json")
    tolerance, close, follows = constants["MARGIN_TOLERANCE"], constants["PROBABILITY_TOLERANCE"], constants["FOLLOWS"]
    one, other = tree_differs()
    first = [row for row in record if row["panel"] == one]
    second = [row for row in record if row["panel"] == other]
    same = [row for row in second if row["cause"] != follows]
    return [
        "GRM runs its own code for each learner.",
        f"On {day(run['started'])} a check compared each learner with a standard version of the same method on {len(versions['grm_fits_sha256'])} small test panels, "
        f"at revision {versions['evaluation_revision']} of the evaluation (commit {run['head'][:7]}).",
        f"It was a software check on small test panels, and it {SCOPE}.",
        "Its plan was fixed before any comparison was read.",
        "The kit and the result table are in the repository, under scripts/learner-crosscheck and reference/learner-crosscheck.",
        "Under the plan's rules, no difference was left without a shown cause.",
        f"The set covering machine gave the same model as an independent implementation that follows GRM's tie order and stop rule, on all {scm['compared']} fits.",
        f"pyscm-ml {libraries['pyscm-ml']} was asked for its best rules at {kit.detail_sum('scm fits', 'steps_read_against_library')} of the steps GRM took, "
        f"and GRM's rule was among them at {kit.detail_sum('scm fits', 'steps_certified')}.",
        f"At the other {kit.detail_sum('scm fits', 'steps_dominated')} the library preferred a rule GRM does not have, one that rejects every sample or none.",
        f"The decision tree agreed with scikit-learn {libraries['scikit-learn']} on {cart['units_matched']} of {cart['compared']} nodes: "
        "both made a leaf, or both split on a feature of the largest gain, GRM on the lowest-numbered one.",
        f"Of the {forest['compared']:,} forest trees, {forest['units_matched']:,} agreed on every node.",
        f"On the other {cart['units_explained']} nodes, and on at least one node of each of the other {forest['units_explained']} trees, "
        "several features had exactly equal gain and GRM did not split on the lowest-numbered one.",
        "Its code gives equal gains to the lowest-numbered feature, but rounding in its own arithmetic made these equal gains differ in the last digits.",
        f"Set against XGBoost {libraries['xgboost']} one tree at a time, {by_tree['units_matched']:,} of {by_tree['compared']:,} boosted trees were identical.",
        f"Of the other {by_tree['units_explained']}, {causes[('tree by tree', 'leaf set to 0')]} came from the one panel with a single class in training, "
        "where XGBoost sets a leaf below its minimum weight to zero and GRM does not.",
        f"Another {causes[('tree by tree', 'equal split')]} divided the samples into the same two groups with another feature, "
        f"and {causes[('tree by tree', 'float32 rounding')]} followed from XGBoost's single-precision arithmetic.",
        f"Run end to end on the {end['compared']} panels with both classes in training, "
        f"the two gave the same number of trees and probabilities within 0.0001 of each other on {end['units_matched']}.",
        f"On {no_tree} panels GRM kept no tree and XGBoost kept one, and the probabilities were equal.",
        "GRM keeps that difference on purpose: a model with no tree gives no call.",
        f"On the other {end['units_explained'] - no_tree} some trees differed, and a record added after the plan gives each of those trees a cause of its own.",
        f"On one, {len(first)} trees each split a node on another feature that divides the fitting samples into the same two groups, and before each of them the two fits "
        f"still gave the fitting samples the same scores within {tolerance:.5f}, so none followed from an earlier tree.",
        f"Its probabilities differed by up to {detail[one]['gbt end to end']['largest_probability_gap']:.3f}, the largest difference in the check, "
        f"and only outside the fitting samples: on those they agreed within {close:.4f}.",
        f"On the other, {len(second)} trees differed.",
        f"Of these, {len(same)} grew from scores that agreed within {tolerance:.5f}: {sum(row['cause'] == 'float32 rounding' for row in same)} differed through XGBoost's "
        f"single-precision arithmetic, and {sum(not row['cause'] for row in same)} have no cause shown.",
        f"The other {len(second) - len(same)} followed an earlier tree that had moved the scores.",
        "The L1 logistic regression was not compared again.",
        f"Its own tests ran at the same commit and {passed} passed; they require each saved case to converge on the best known fit.",
        "One more of its tests needs files that are not in the repository and did not run.",
        # history: the check of October 2026 has no file behind it, so these two are held as written
        "An earlier check, in October 2026 under revision 7, reported that the set covering machine, the decision tree and the random forest "
        "matched their reference on all 16 cases.",
        "It also showed that the L1 solver of that time stopped early on some cases, and revision 8 replaced that solver.",
        "The files of that check are lost, so that statement cannot be checked, and it stays here only as history.",
    ]


class Methods(unittest.TestCase):
    def test_the_part_is_these_sentences_and_no_other(self):
        text = methods_part()
        for sentence in sentences():
            self.assertTrue(text.startswith(sentence), f"where the part has {text[:160]!r}..., this test holds the sentence {sentence!r}")
            text = text.removeprefix(sentence)
            self.assertRegex(text, "^( |$)", f"the part goes on {text[:80]!r} after {sentence!r}")
            text = text.removeprefix(" ")
        self.assertEqual(text, "", "the part has text after the last sentence this test holds")

    def test_what_the_sentences_take_for_granted(self):
        differences, detail = kit.differences(), kit.load("detail.json")
        totals = {(row["learner"], row["check"]): row for row in kit.load("totals.json")}
        # "No difference was left without a shown cause."
        self.assertEqual(sum(row["units_unexplained"] + row["panels"]["unexplained"] for row in totals.values()), 0)
        self.assertEqual([row for row in differences if not row["cause"]], [])
        # "the same model as an independent implementation that follows GRM's tie order and stop rule, on all ... fits"
        self.assertEqual(kit.detail_sum("scm fits", "independent_same"), totals[("scm", "fits")]["compared"])
        # "at the other ..." steps: certified and dominated are all the steps read
        read = kit.detail_sum("scm fits", "steps_read_against_library")
        self.assertEqual(kit.detail_sum("scm fits", "steps_certified") + kit.detail_sum("scm fits", "steps_dominated"), read)
        # "the one panel with a single class in training"
        self.assertEqual({row["panel"] for row in differences if row["cause"] == "leaf set to 0"}, {"one_class_train"})
        # the three causes of a tree that differs are all the causes there are
        by_tree = {row["cause"] for row in differences if row["check"] == "tree by tree"}
        self.assertEqual(by_tree, {"leaf set to 0", "equal split", "float32 rounding"})
        # "the 15 panels with both classes in training": the one not run is the one-class panel
        self.assertEqual(totals[("gbt", "end to end")]["panels"]["not applicable"], 1)
        self.assertIn("reason", detail["one_class_train"]["gbt end to end"])
        # "GRM kept no tree and XGBoost kept one, and the probabilities were equal"
        for row in differences:
            if row["cause"] == "no tree kept":
                end = detail[row["panel"]]["gbt end to end"]
                self.assertEqual((end["grm_trees"], end["xgboost_trees"], end["largest_probability_gap"]), (0, 1, 0))
        # the end-to-end causes are those two and "a tree differs", on two panels
        end = {row["cause"] for row in differences if row["check"] == "end to end"}
        self.assertEqual(end, {"no tree kept", "a tree differs (equal split)", "a tree differs (float32 rounding)"})
        one, other = tree_differs()
        constants, record = kit.constants("check_gbt"), kit.trees()
        tolerance, follows = constants["MARGIN_TOLERANCE"], constants["FOLLOWS"]
        # "On one, ... trees each split a node on another feature that divides the fitting samples into the same two groups, and
        # before each of them the two fits still gave the fitting samples the same scores ..., so none followed from an earlier tree."
        first = [row for row in record if row["panel"] == one]
        self.assertEqual({row["cause"] for row in first}, {"equal split"})
        self.assertLessEqual(max(row["fitting_margin_gap"] for row in first), tolerance)
        # "Its probabilities differed by up to ..., the largest difference in the check, and only outside the fitting samples"
        gaps = [panel["gbt end to end"].get("largest_probability_gap", 0) for panel in detail.values()]
        self.assertEqual(detail[one]["gbt end to end"]["largest_probability_gap"], max(gaps))
        self.assertLessEqual(detail[one]["gbt end to end"]["largest_probability_gap_on_the_fitting_samples"], constants["PROBABILITY_TOLERANCE"])
        # "On the other, ... trees differed. Of these, ... grew from scores that agreed within ...: ... single-precision ..., and ...
        # have no cause shown. The other ... followed an earlier tree that had moved the scores."
        second = [row for row in record if row["panel"] == other]
        self.assertEqual({row["cause"] for row in second}, {"float32 rounding", "", follows})
        for row in second:
            self.assertEqual(row["cause"] == follows, row["fitting_margin_gap"] > tolerance, row)
        # the first of them is the tree that differences.tsv gives the panel's cause from
        self.assertEqual(second[0]["cause"], "float32 rounding")
        self.assertEqual(int(second[0]["unit"].split()[1]), detail[other]["gbt end to end"]["first_tree_not_identical"])
        # "a model with no tree gives no call" is the code's own statement
        gbt = (kit.ROOT / "src" / "learn" / "gbt.rs").read_text(encoding="utf-8")
        self.assertIn("//! A model with no tree is a no-call, by design", gbt)

    def test_the_part_says_what_the_check_is_and_is_not(self):
        text = methods_part()
        self.assertIn(f"It was a software check on small test panels, and it {SCOPE}.", text)
        self.assertIn("Under the plan's rules, no difference was left without a shown cause.", text)
        self.assertIn("under scripts/learner-crosscheck and reference/learner-crosscheck.", text)
        self.assertIn("The files of that check are lost, so that statement cannot be checked, and it stays here only as history.", text)
        self.assertTrue((kit.ROOT / "scripts" / "learner-crosscheck").is_dir())

    def test_no_text_calls_the_panels_what_the_panel_file_does_not(self):
        self.assertNotIn("synthetic", kit.FIXTURES.read_text(encoding="ascii").lower())
        for text in (methods_part(), (kit.REFERENCE / "README.md").read_text(encoding="utf-8"), (kit.KIT / "README.md").read_text(encoding="utf-8")):
            self.assertNotIn("synthetic", text.lower())


class ResultReadme(unittest.TestCase):
    def setUp(self):
        self.text = (kit.REFERENCE / "README.md").read_text(encoding="utf-8")

    def test_each_table_is_the_result_s(self):
        for make in kit.TABLES:
            with self.subTest(table=make.__name__):
                self.assertIn("\n\n" + "\n".join(make()) + "\n\n", self.text)

    def test_no_other_table_row_carries_a_count(self):
        rebuilt = {line for make in kit.TABLES for line in make()}
        other = [line for line in self.text.splitlines() if line.startswith("|") and line not in rebuilt]
        # the table under "How it was made": its header, its rule and seven rows
        self.assertEqual([line.split("|")[1].strip() for line in other], ["", "---", "GRM", "Kit", "Compiler", "Panels", "References", "Command", "Run again"])

    def test_the_sha256_lines(self):
        block = re.findall(r"\n```\n(.*?)\n```\n", self.text, re.S)
        self.assertEqual(len(block), 1)
        listed = dict(reversed(line.split("  ")) for line in block[0].splitlines())
        files = sorted(path.name for path in kit.REFERENCE.iterdir() if path.name != "README.md")
        self.assertEqual(sorted(listed), files)
        self.assertEqual(list(listed), files)
        for name in files:
            self.assertEqual(listed[name], hashlib.sha256((kit.REFERENCE / name).read_bytes()).hexdigest(), name)

    def test_how_it_was_made(self):
        run, rerun, versions = kit.load("run.json"), kit.load("rerun.json"), kit.load("versions.json")
        libraries = versions["libraries"]
        for expected in (
            f"| GRM | `main` at commit `{run['head']}`, evaluation revision {versions['evaluation_revision']}, `src/` tree `{run['src_tree']}`, unchanged by the kit |",
            f"| Kit | first committed at `{KIT_COMMIT}`, on top of `{run['head'][:7]}`, which holds none of its files.",
            f"| Compiler | {run['rustc'].removeprefix('rustc ').join(['rustc ', ''])}, release build |",
            f"| Panels | the {len(versions['grm_fits_sha256'])} of `tests/data/lasso-stopping-fixtures.txt`, sha256 `{versions['fixtures_sha256']}` |",
            f"Python {versions['python']} with scikit-learn {libraries['scikit-learn']}, XGBoost {libraries['xgboost']} (the CPU wheel) and pyscm-ml {libraries['pyscm-ml']}, "
            f"on numpy {libraries['numpy']} and scipy {libraries['scipy']}",
            f"from {run['started'][11:19]} to {run['ended'][11:19]} UTC on {day(run['started'])} (`run.json`)",
            f"| Run again | `scripts/learner-crosscheck/run.sh ENV OUT` at commit `{rerun['head']}`, from {rerun['started'][11:19]} to {rerun['ended'][11:19]} UTC "
            f"on {day(rerun['started'])}, in another environment built from the same pinned files (`rerun.json`) |",
            f"# Learner cross-check: the result of {day(run['started'])}",
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected, self.text)

    def test_it_says_what_the_check_is_and_is_not(self):
        self.assertIn(f"This is a software check on small test panels. It {SCOPE}", self.text)
        self.assertIn("**Under the plan's rules, no difference was left without a shown cause.**", self.text)
        self.assertIn("They are labelled exploratory and no outcome depends on them.", self.text)
        self.assertIn("Exploratory, written after the first comparison was read:", self.text)

    def test_the_lasso_lines(self):
        lasso = (kit.REFERENCE / "lasso-suite.txt").read_text(encoding="utf-8").splitlines()[-1]
        counts = re.fullmatch(r"test result: ok\. (\d+) passed; (\d+) failed; (\d+) ignored; .*", lasso).groups()
        self.assertIn("{} passed, {} failed, {} ignored (`lasso-suite.txt`)".format(*counts), self.text)


class KitReadme(unittest.TestCase):
    def setUp(self):
        self.text = (kit.KIT / "README.md").read_text(encoding="utf-8")

    def test_it_names_every_file_of_the_kit(self):
        files = sorted(path.name for path in kit.KIT.iterdir() if path.is_file() and path.name != "README.md")
        self.assertEqual(
            files,
            ["check_env.py", "check_fits.sh", "check_gbt.py", "check_scm.py", "check_trees.py", "compare.py", "exact.py", "exploratory.py", "main.rs", "panels.py",
             "requirements-build.txt", "requirements-scm.txt", "requirements.txt", "run.sh", "scm_reference.py", "table.py"],
        )  # fmt: skip
        for name in files:
            self.assertIn(name, self.text)

    def test_the_install_commands_check_every_hash(self):
        commands = [line for line in self.text.splitlines() if line.startswith("ENV/bin/pip ")]
        self.assertEqual(len(commands), 4)
        # the three requirement files, each with its hashes; pyscm-ml is built into ENV/built, where check_env.py reads it
        for command, name in zip(commands, kit.module("check_env").FILES):
            self.assertIn("--no-deps --require-hashes", command)
            self.assertTrue(command.endswith(f"-r scripts/learner-crosscheck/{name}"), command)
        self.assertTrue(commands[2].startswith("ENV/bin/pip wheel --no-deps --require-hashes --no-build-isolation --no-binary pyscm-ml -w ENV/built "))
        # then that wheel, and nothing from an index
        self.assertEqual(commands[3], "ENV/bin/pip install --no-deps --no-index ENV/built/pyscm_ml-1.1.2-*.whl")

    def test_it_names_the_versions_that_ran(self):
        libraries = kit.load("versions.json")["libraries"]
        for expected in (f"pyscm-ml {libraries['pyscm-ml']}", f"scikit-learn {libraries['scikit-learn']}", f"XGBoost {libraries['xgboost']}"):
            self.assertIn(expected, self.text)

    def test_it_says_what_the_check_is_and_is_not(self):
        self.assertIn(f"It is a software check on small test panels. It {SCOPE}", self.text)
        self.assertIn("It reads no genome data.", self.text)

    def test_the_example_is_in_the_cargo_file_and_the_suite_in_ci(self):
        cargo = (kit.ROOT / "Cargo.toml").read_text(encoding="utf-8")
        self.assertIn('[[example]]\nname = "learner-crosscheck"\npath = "scripts/learner-crosscheck/main.rs"\n', cargo)
        workflow = (kit.ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        self.assertTrue("          - suite: tests/learner_crosscheck\n            top: .\n" in workflow, "the python matrix of ci.yml lacks the suite")
        self.assertTrue("`tests/learner_crosscheck`" in (kit.ROOT / "docs" / "ci.md").read_text(encoding="utf-8"), "docs/ci.md does not name the suite")

    def test_ci_runs_check_fits_in_the_rust_test_job(self):
        workflow = (kit.ROOT / ".github" / "workflows" / "ci.yml").read_text(encoding="utf-8")
        jobs = re.split(r"\n  (?=[a-z][a-z-]*:\n)", workflow.split("\njobs:\n", 1)[1])
        job = [text for text in jobs if text.startswith("rust-test:\n")]
        self.assertEqual(len(job), 1, "ci.yml has no job rust-test")
        # the job's steps, each its lines without comments and blank lines; the step is these two lines and no other
        steps = []
        for line in job[0].splitlines()[1:]:
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            if line.startswith("      - "):
                steps.append([line])
            elif steps and line.startswith("        "):
                steps[-1].append(line)
        step = ["      - name: Check the learner cross-check's fits", "        run: bash scripts/learner-crosscheck/check_fits.sh"]
        self.assertIn(step, steps, "the rust-test job of ci.yml does not run check_fits.sh")
        self.assertNotRegex(job[0], r"\n    (if|continue-on-error):")
        self.assertIn("`bash scripts/learner-crosscheck/check_fits.sh`", (kit.ROOT / "docs" / "ci.md").read_text(encoding="utf-8"))

    def test_docs_ci_gives_the_suite_s_count(self):
        from . import test_kit, test_table

        loader = unittest.TestLoader()
        tests = sum(loader.loadTestsFromModule(module).countTestCases() for module in (test_kit, test_table, sys.modules[__name__]))
        self.assertIn(f"`tests/learner_crosscheck` {tests},", (kit.ROOT / "docs" / "ci.md").read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
