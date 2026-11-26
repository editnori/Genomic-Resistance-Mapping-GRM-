"""docs/combiner-result.md against its script, scripts/mechanism/combiner.py, before any reading.

This test file reads the page and the script's own constants, gates and fault list. It reads no record, fits
nothing and reads no genome.
"""

import hashlib
import math
import re
import unittest

from tests.mechanism import common

cb = common.load("combiner")
cm = common.load("common")
faults = common.load("combiner_faults")
lo = common.load("learner_overlap")

PAGE = common.ROOT / "docs" / "combiner-result.md"
SCRIPT = common.KIT / "combiner.py"
RECORD = common.ROOT / "reference" / "combiner"
VOTE_POOLED = common.ROOT / "reference" / "learner-overlap" / "second-pass" / "vote-pooled.tsv"
# FREEZE-COMBINE's first prediction and its third as restated on 25 Nov 2026, word for word, before any reading
FIRST_PREDICTION = ("The weighted mean beats the equal mean by the 0.02 margin on **fewer than 1 in 5** of the reports, "
                    "and not on the pooled median. About 1 in 4 that the pooled median gain clears 0.02 for either arm.")
THIRD_PREDICTION = ("Of the **8,067** resistant genomes the equal mean still misses, the better of the two arms "
                    "recovers **under a third**, and it does not do so without giving back susceptible genomes among "
                    "the 8,896 the equal mean already calls resistant.")
SECTIONS = ["Plan", "Reading", "Limits", "Method"]
FIXED = ("A lineage split is a retrospective split inside one database, never external validation.",
         "\"No difference shown\" is never written as \"equal\".",
         "The error levels serve as an adapted research screen and are never an applicable FDA validation criterion.",
         "Coverage is the share of genomes answered and is never accuracy.")
COMMANDS = ("python3 -I scripts/mechanism/combiner.py read REPORTS --first reference/learner-overlap "
            "--second reference/learner-overlap/second-pass --out READING",
            "python3 -I scripts/mechanism/combiner.py write READING --reports REPORTS --out OUT --sha256 SHA256",
            "python3 -I scripts/mechanism/combiner.py print OUT")


def sections(text):
    parts = re.split(r"^## (.+)$", text, flags=re.M)
    return dict(zip(parts[1::2], parts[2::2]))


def power(x):
    """A power of ten as the page writes it, 1e-10 for 1e-10."""
    e = round(math.log10(x))
    assert float(f"1e{e}") == x, x
    return f"1e{e}"


def ordinal(n):
    return f"{n:,}" + ({1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th") if n % 100 not in (11, 12, 13) else "th")


class Page(unittest.TestCase):
    def setUp(self):
        self.text = PAGE.read_text(encoding="utf-8")
        self.sections = sections(self.text)

    def test_the_sections_stand_in_order(self):
        self.assertEqual(list(self.sections), SECTIONS)

    def test_the_fixed_words_stand_in_the_limits(self):
        for words in FIXED:
            self.assertIn(words, self.sections["Limits"])
        # FREEZE-COMBINE's caveat, which the page must state in its own words
        for words in ("not one fixed model's outputs across folds", "a reweighting of the release's own per-fold "
                      "predictions, not as a model that could be shipped unchanged", "would need a refit of the base "
                      "learners inside one fold structure and the maintainer's word under R4"):
            self.assertIn(words, self.sections["Limits"])

    def test_no_difference_shown_is_never_written_as_equal(self):
        # the freeze's predictions are quoted word for word, and only there does "the equal mean" stand
        rest = "\n".join(line for line in self.text.replace("\"equal\"", "").split("\n") if not line.startswith("> "))
        self.assertEqual(re.findall(r"\bequal\b(?! weights)", rest, flags=re.I), [])
        self.assertIn("no difference shown", self.sections["Plan"])

    def test_the_plan_gives_the_scripts_constants(self):
        plan = self.sections["Plan"]
        lo_value, hi_value = cm.LO_INDEX + 1, cm.HI_INDEX + 1
        for words in (f"a score of {cb.CUTOFF} or more", f"sum to {cb.VOTE_SUM} or more by `math.fsum`",
                      f"clipped to [{cb.CLIP:.4f}, {1 - cb.CLIP:.4f}]", f"is {cb.MARGIN} or more",
                      f"{cm.RESAMPLES:,} resamples of species", f"the {ordinal(lo_value)} and {ordinal(hi_value)}",
                      f"under {cm.MIN_SPECIES} species", "five learners on the lineage split and 10 or more resistant"):
            self.assertIn(words, plan)
        self.assertEqual(cb.FORMS, ("mean", "weighted", "logistic"))
        self.assertEqual(lo_value, 101)

    def test_the_stop_rule_is_dated_and_its_caps_are_a_limit(self):
        plan, limits = self.sections["Plan"], self.sections["Limits"]
        self.assertIn("This stop rule was ruled on 25 Nov 2026 (R19), after the offline figures were read.", plan)
        for words in (f"its cap of {cb.W_CAP:,} steps", f"{cb.L_CAP} Newton steps", "a step with no finite solution, "
                      "where the last finite iterate stands", "A held-out score that is not finite stops the read.",
                      "a second count without the arm's marked reports is written, and it decides nothing"):
            self.assertIn(words, plan)
        # R37: a fit at its cap had stalled at its best fit, in place of "depends on the cap"
        self.assertIn(f"a logistic fit ended at its cap of {cb.L_CAP} Newton steps had stalled at its best fit, and 0 "
                      "of 4,000 calls differed from a second fit's.", limits)
        self.assertIn(f"A fold at either cap, {cb.W_CAP:,} steps for the weighted mean, still marks its report, and the "
                      "count without marked reports decides nothing.", limits)
        self.assertNotIn("depends on the cap", self.text)
        self.assertIn(f"R19 confirmed the wording: {cb.VME_RULE}.", plan)

    def test_the_plan_states_r30s_rules(self):
        plan = self.sections["Plan"]
        for words in ("gets a coefficient of exactly 0 and is left out of that fold's fit",
                      "are fitted as one column and share its coefficient in even parts",
                      "Both rules read the stored scores as they are, with no tolerance, and both can apply in one fold",
                      "the least sum of squares of the five coefficients",
                      "its calls and coefficients do not depend on the learners' order",
                      "moved 7 of 700 ends between the cap and converged",
                      "On a fold with neither case, the fit is the one before 26 Nov, to the bit",
                      "A fold with either case is not marked",
                      "that fold is not scored, and the read stops at gate 5",
                      "No fold is scored from coefficients left at 0",
                      "the weighted mean is kept only if the rule keeps it in both",
                      f"gives every genome of a fold of {cb.FLAG_FOLD} or more genomes one score",
                      "The record gives each count's own outcome for the weighted mean and the arm's reading.",
                      f"Gate 1 stops the read unless it reads exactly those {cb.FROZEN_REPORTS} reports",
                      "A report with one fold, or a genome with no fold number, stops the read"):
            self.assertIn(words, plan)
        self.assertEqual(cb.COUNTS_OF_RULE, ("all", "without_marked", "without_flagged"))
        self.assertEqual(cb.plan()["counts"], {"all": "decides", "without_marked": "decides nothing",
                                               "without_flagged": "decides"})

    def test_the_plan_states_r37s_rules(self):
        plan, limits = self.sections["Plan"], self.sections["Limits"]
        rows = re.findall(r"^\| (One score|Coefficient under 0|Training share) \((R\d+)\) \| (.+) \| "
                          r"(both|logistic only) \|$", plan, flags=re.M)
        self.assertEqual([(name, ruling, arms) for name, ruling, _, arms in rows],
                         [("One score", "R30", "both"), ("Coefficient under 0", "R37", "logistic only"),
                          ("Training share", "R46", "both")])
        self.assertEqual(cb.FLAGS, {"weighted": ("one_score_fold", "training_share_fold"),
                                    "logistic": ("one_score_fold", "coefficient_under_0", "training_share_fold")})
        self.assertEqual(cb.plan()["flags"], {a: list(cb.FLAGS[a]) for a in cb.ARMS})
        self.assertIn("a coefficient under 0; a coefficient of exactly 0 does not flag", rows[1][2])
        for words in ("the count without the arm's flagged reports decide",
                      f"A count of under {cm.MIN_SPECIES} species has no interval, so it keeps nothing, and the "
                      "weighted mean's reading is no difference shown, whatever the other count gives.",
                      "each read on the stored values with no tolerance",
                      "The weighted mean's weights are never under 0."):
            self.assertIn(words, plan)
        self.assertIn("A coefficient under 0 also arises when two learners score nearly alike, with no fold's share in "
                      "either.", limits)

    def test_the_plan_states_r43s_rules(self):
        # R43: the logistic arm is descriptive only, and why, with the review's figures from made-up reports; the
        # weighted mean is the one arm the rule can keep; the learners' code; the share's sign; and the incumbent
        plan, limits = self.sections["Plan"], self.sections["Limits"]
        self.assertEqual((cb.CAN_BE_KEPT, cb.DESCRIPTIVE), (("weighted",), "descriptive only"))
        self.assertEqual((cb.plan()["can_be_kept"], cb.plan()["descriptive_only"]), (["weighted"], ["logistic"]))
        for words in ("reads the logistic arm as descriptive only. Its three counts, medians and intervals are written "
                      "as before, and its outcome and reading say \"descriptive only\" in each count, whatever its "
                      "figures.",
                      "R37's condition, its column, the logistic arm's third count and the `flags:` line stay as "
                      "counts, and decide nothing.",
                      "correlated at 0.8, and again at 0.95, R37's condition flagged 200 of 200 honest made-up reports",
                      "kept in 0 of 8 made-up readings whose true gain was near 0.05",
                      "got the logistic arm kept through the fold's share alone, and R37's condition flagged 0 of its "
                      "20 reports",
                      "The weighted mean is the one arm the rule can keep. It has no intercept and no weight under 0.",
                      "its median gain in balanced accuracy was -0.0096 over 20 reports, on the reviewer's made-up "
                      "reports.",
                      "Four of the five, the decision tree, the forest, the boosted trees and the lasso, fit with class "
                      "weights that give each class half the total weight. The rule learner gives a 0 or 1 call.",
                      "An empty tree or forest scores its fold 0.5. R30's condition flags both. No learner scores a "
                      "part of a fold at the share. This is a reading of the code and no measurement.",
                      "With every coefficient at 0 or above, the share moves a fold's calls nearly always the wrong "
                      "way: on the reviewer's made-up reports, 3 of 248 fold fits moved them the right way.",
                      "The logistic arm is descriptive only."):
            self.assertIn(words, plan)
        self.assertIn("A learner that carries the share lowers the incumbent, the mean with equal weights, so the "
                      "weighted mean can gain beside it.", limits)
        # "either arm" stands only in the freeze's prediction, quoted as written, and no line gives the logistic arm
        # a rule of its own
        rest = [line for line in self.text.split("\n") if not line.startswith("> ")]
        self.assertEqual([line for line in rest if "either arm" in line or "an arm is kept" in line.lower()], [])
        self.assertNotIn("the wrong way or nowhere", self.text)

    def test_write_works_the_flags_out_again(self):
        method = self.sections["Method"]
        for words in ("Write works each report's sign flag, R37's condition, and its logistic mark out again from the "
                      "logistic fold records it parses, and stops on a difference with a line, `stop flags:`, that "
                      "names the report by species and antibiotic and gives no figure.",
                      "Write and print stop, with a line `stop rule:`, on a rule row that gives the logistic arm an "
                      "outcome or a reading other than descriptive only.",
                      "a line, `flags:`, with how many reports each condition flags, how many have both R30's and "
                      "R37's, and how many R46's flags that R30's does not, and no score. Those counts decide nothing "
                      "(R43)."):
            self.assertIn(words, method)
        self.assertNotIn("It trusts the read's R30 flag and weighted mark", self.text)
        script = SCRIPT.read_text(encoding="utf-8")
        for printed in ('wrote {READING}, sha256 {written}"', 'print(f"write: read {READING}, sha256 {digest}")',
                        'f"flags: {len(rows)} reports; ', 'f"stop flags: {d.get(\'species\')} with {d.get(\'antibiotic\')}: ',
                        'f"stop rule: the {arm} arm has an outcome or a reading other than {DESCRIPTIVE} in {where}: STOP"'):
            self.assertEqual(script.count(printed), 1, printed)

    def test_the_plan_states_r46s_rules(self):
        # R46: the gap with the reviewer's figures, the rule that closes it, the learners' code, its limit, and write
        # and print as the script stands
        plan, limits, method = self.sections["Plan"], self.sections["Limits"], self.sections["Method"]
        self.assertEqual(cb.SHARE_GENOMES, 2)
        self.assertEqual(cb.plan()["share_genomes"], cb.SHARE_GENOMES)
        self.assertIn("one division", cb.plan()["training_share_fold"])
        self.assertIn(f"some learner gives {cb.SHARE_GENOMES} or more genomes of one fold that fold's training share "
                      "of resistant genomes as their score, by `==`", plan)
        for words in ("a learner scored about 97 of 100 genomes of each fold at that fold's training share. The plain "
                      "mean's median balanced accuracy fell by 0.0381 and 0.0320, and the weighted mean's by 0.0000. "
                      "The weighted mean was kept in 10 of 10 readings of two such families, and in 0 of their 10 "
                      "controls. R30's condition flagged 3 of 200 reports, as no fold had one score. R46 closes it "
                      "with the training share condition above.",
                      "The release works a fold's training share by one division: resistant training genomes over "
                      "training genomes, the report's genomes outside the fold. Only a boosted model with no tree, or "
                      "a lasso with no weight, gives it as a score, to every genome of its fold.",
                      "Each fold's counts enter every fit of those four, and a leaf of fixed counts scores lower as the "
                      "share rises. In 6 made-up families with a falling learner, two mimicking the class weights, the "
                      "weighted mean was kept wrongly in 0 of 30 readings.",
                      "gives no fold one score and no genome the share, so neither R30's condition nor R46's flags it"):
            self.assertIn(words, plan)
        self.assertIn("A tree leaf can score exactly the share by chance, when the leaf's ratio of resistant to "
                      "susceptible genomes is the square of the training set's. Such a report leaves the third count, "
                      "and the `flags:` line counts it. A made-up learner within a hair of the share is not flagged.",
                      limits)
        for words in ("It takes the sha256 the read printed, works the sha256 of the reading's bytes before it parses "
                      "them, and stops with one line, `stop sha256:`, on a difference or on a value that is not 64 hex "
                      "characters. It then writes nothing. The record's `gates.txt` carries that sha256. Write does not "
                      "work R30's flag or R46's out again: the sha256 covers every edit.",
                      "A fold record that lacks a key write reads stops write with a line, `stop fold records:`, that "
                      "names the report in the same way.",
                      "Print reads the record's `SHA256SUMS` first. It stops with one line, `stop record:`, naming the "
                      "file, when the list is missing, a file of the record is not listed, or a listed file is missing "
                      "or its bytes give another sha256.",
                      "or that lacks R46's flag or names other flags."):
            self.assertIn(words, method)
        for gone in ("The training share enters a stored score in one place", "No learner is built to fall",
                     "sets no rule for it", "R30's condition alone misses it"):
            self.assertNotIn(gone, self.text)
        script = SCRIPT.read_text(encoding="utf-8")
        for printed in ('"stop sha256: the value given is not 64 hex characters: STOP"',
                        'f"stop sha256: {READING} is not the reading the read wrote: its bytes give another sha256: STOP"',
                        'f"stop fold records: {d.get(\'species\')} with {d.get(\'antibiotic\')}: {what}: STOP"',
                        'f"stop record: {SUMS} is missing: STOP"', 'f"stop record: {name} is not listed in {SUMS}: STOP"',
                        'f"stop record: {name} gives another sha256 than {SUMS} lists: STOP"'):
            self.assertEqual(script.count(printed), 1, printed)

    def test_the_page_keeps_its_size(self):
        # R37: 30 list items or fewer and 12 headers or fewer
        lines = self.text.split("\n")
        self.assertLessEqual(sum(1 for line in lines if re.match(r"^\s*(- |\d+\. )", line)), 30)
        self.assertLessEqual(sum(1 for line in lines if line.startswith("#")), 12)

    def test_the_plan_names_the_logistic_stop_values(self):
        plan = self.sections["Plan"]
        script = SCRIPT.read_text(encoding="utf-8")
        whatever = re.findall(r"or t_ < (1e-\d+):", script)
        self.assertEqual(len(whatever), 1)
        for words in ("starts with the intercept and the five coefficients at 0",
                      f"no entry of the mean gradient is above {power(cb.L_TOL)}",
                      f"by more than {power(cb.L_STEP)} × (1 + the largest of the six in size)",
                      f"a pivot under {power(cb.PIVOT)} has no finite solution",
                      f"once the step is under {whatever[0]} of its full length, it takes it whatever the loss",
                      f"{cb.L_CAP} Newton steps"):
            self.assertIn(words, plan)
        self.assertEqual(cb.plan()["logistic_fit"]["start"], "zero")

    def test_the_plan_names_the_intervals_seed(self):
        plan = self.sections["Plan"]
        script = SCRIPT.read_text(encoding="utf-8")
        self.assertEqual(script.count('interval(dict(by), f"combiner|{arm}|{name}")'), 1)
        name = "combiner|weighted|balanced_accuracy"
        self.assertEqual(cm.seed_of(name), int(hashlib.sha256(name.encode("utf-8")).hexdigest(), 16))
        self.assertIn("seeded by the sha256 of the row's own name read as one whole number", plan)
        self.assertIn(f"`{name}`", plan)

    def test_the_predictions_are_quoted_as_written(self):
        # each figure the record holds is checked against it: the 8,067 resistant genomes the mean misses that 1 to
        # 4 of the five call, and the 8,896 susceptible genomes it calls resistant, over the 268 reports
        plan = self.sections["Plan"]
        self.assertIn("> " + FIRST_PREDICTION + "\n", plan)
        self.assertIn("> " + THIRD_PREDICTION + "\n", plan)
        pooled = [r for r in lo.read_tsv(str(VOTE_POOLED)) if r["scope"] == "all"]
        self.assertEqual(len(pooled), 1)
        r = pooled[0]
        self.assertEqual((f"{int(r['mean_miss_by_1to4']):,}", f"{int(r['mean_wrong']):,}", int(r["datasets"])),
                         ("8,067", "8,896", cb.FROZEN_REPORTS))
        self.assertEqual(f"{cb.MARGIN}", "0.02")
        self.assertIn(f"The 8,067 of {int(r['resistant']):,} resistant genomes are those the mean misses and 1 to 4 of "
                      f"the five call. The 8,896 of {int(r['susceptible']):,} susceptible genomes are those it calls "
                      f"resistant. Both are over {r['datasets']} reports", plan)

    def test_the_reading_waits_for_its_record(self):
        if not RECORD.exists():
            self.assertEqual(self.sections["Reading"].strip(), "Not read yet.")

    def test_the_commands_are_the_scripts(self):
        method = self.sections["Method"]
        for command in COMMANDS:
            self.assertIn(command + "\n", method)

    def test_the_gate_table_names_each_gate_of_the_script(self):
        gates = sorted(set(re.findall(r"\"gate (\d) ([a-z -]+):", SCRIPT.read_text(encoding="utf-8"))))
        self.assertEqual([n for n, _ in gates], [str(i) for i in range(1, 9)])      # one name for each gate
        rows = re.findall(r"^\| (\d)\. ([a-z -]+) \|", self.sections["Method"], flags=re.M)
        self.assertEqual(rows, gates)

    def test_the_fault_table_names_each_planted_faults_test(self):
        method = self.sections["Method"]
        rows = re.findall(r"^\| (C\d+)\. (.+) \| `([A-Za-z][A-Za-z0-9]*\.test_\w+)` \|$", method, flags=re.M)
        self.assertEqual(rows, [(fault, words[0].upper() + words[1:], test.rsplit(".", 2)[-2] + "." + test.rsplit(".", 1)[-1])
                                for fault, words, test, _ in faults.FAULTS])
        self.assertIn(f"{len(faults.FAULTS)} of {len(faults.FAULTS)} faults failed their test.", method)

    def test_the_test_counts_are_the_suites(self):
        loader = unittest.TestLoader()
        script = loader.loadTestsFromName("tests.mechanism.test_combiner").countTestCases()
        page = loader.loadTestsFromName("tests.mechanism.test_combiner_page").countTestCases()
        method = self.sections["Method"]
        self.assertIn(f"`tests/mechanism/test_combiner.py` holds {script} SYNTHETIC tests of the script.", method)
        self.assertIn(f"`tests/mechanism/test_combiner_page.py` holds {page} tests of this page against the script.",
                      method)


if __name__ == "__main__":
    unittest.main()
