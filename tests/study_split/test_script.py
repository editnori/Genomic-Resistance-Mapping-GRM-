"""scripts/study-split/record.py on SYNTHETIC inputs: cells worked by hand, the same bytes, stops and --check.

The inputs are built in a temporary directory by synthetic.make; none of them is data.
"""

import contextlib
import io
import math
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction

from tests.study_split import agree, synthetic
from tests.study_split.common import ROOT, SCRIPT, load_script, read_tsv

record = load_script()
D1, D2, D3, D4 = synthetic.ENTERED


def run(inputs, *args):
    """main() on a synthetic input set; returns (rc, stderr)."""
    err = io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(io.StringIO()):
        rc = record.main(["--inputs", inputs, *args])
    return rc, err.getvalue()


def files(directory):
    return {p.relative_to(directory).as_posix(): p.read_bytes()
            for p in sorted(pathlib.Path(directory).rglob("*")) if p.is_file()}


def rows(out, name, *key):
    return {tuple(r[k] for k in key): r for r in read_tsv(os.path.join(out, name))}


class Pieces(unittest.TestCase):
    def test_measure_fit_by_hand(self):
        # Two resistant genomes scored 0.9 and 0.4, two susceptible scored 0.5 and 0.0: a score of 0.5 calls
        # resistant, and the 0.0 of a susceptible genome is clipped to lose -log(1 - 0.0001).
        tp, fn, tn, fp, ll = record.measure_fit([True, True, False, False], [0.9, 0.4, 0.5, 0.0])
        self.assertEqual((tp, fn, tn, fp), (1, 1, 1, 1))
        self.assertAlmostEqual(ll, -(math.log(0.9) + math.log(0.4) + math.log(0.5) + math.log(0.9999)) / 4, places=15)
        # A wrong 1 and a wrong 0 lose the same amount.
        self.assertEqual(record.measure_fit([True], [0.0])[4], record.measure_fit([False], [1.0])[4])
        self.assertAlmostEqual(record.measure_fit([True], [0.0])[4], -math.log(0.0001), places=12)

    def test_scm_log_loss_in_exact_form(self):
        # Two 0 or 1 fits of 40 genomes with 3 errors each, in a different order: the summed losses may part by a
        # rounding unit, and the exact form gives one value within 1e-12 of both.
        a = record.measure_fit([True] * 20 + [False] * 20, [0.0] * 3 + [1.0] * 17 + [0.0] * 20)
        b = record.measure_fit([False] * 20 + [True] * 20, [1.0] * 2 + [0.0] * 18 + [1.0] * 19 + [0.0])
        self.assertEqual(a[:4], (17, 3, 20, 0))
        self.assertEqual(b[:4], (19, 1, 18, 2))
        va, vb = (record.loss_value("scm", x[4], *x[:4]) for x in (a, b))
        self.assertEqual(va, vb)
        self.assertEqual(va, record.RIGHT_LOSS + record.PER_ERROR * Fraction(3, 40))
        for x in (a, b):
            self.assertLess(abs(float(va) - x[4]), 1e-12)
        # Other learners keep the log loss as written.
        self.assertEqual(record.loss_value("rf", a[4], *a[:4]), Fraction(a[4]))
        self.assertEqual(record.arm_loss({"learner": "scm", "study_tp": "17", "study_fn": "3", "study_tn": "20",
                                          "study_fp": "0", "study_log_loss": "9"}, "study"), va)

    def test_fit_cells_by_hand(self):
        cells = record.fit_cells(3, 1, 6, 2, 0.25)
        self.assertEqual(cells, ["3", "1", "6", "2", "0.75", "0.75", "0.75", "0.25", "0.25"])

    def test_interval_by_hand(self):
        # Three species with one value each: of the 27 ordered draws, 7 give median 1, 13 give 2 and 7 give 3.
        self.assertEqual(record.interval([[1], [2], [3]]), (2, 1, 3))
        # A species with two datasets enters a draw with both: AA gives [0, 0, 4, 4], median 2, in 1 draw of 4.
        self.assertEqual(record.interval([[0, 4], [1]]), (1, 1, 2))
        self.assertEqual(record.interval([[Fraction(5, 2)]]), (Fraction(5, 2), Fraction(5, 2), Fraction(5, 2)))

    def test_interval_against_every_ordered_draw(self):
        r = random.Random("study split interval")
        for _ in range(40):
            k = r.randint(1, 5)
            groups = [[Fraction(r.randint(-6, 6), r.randint(1, 3)) for _ in range(r.randint(1, 3))] for _ in range(k)]
            with self.subTest(groups=groups):
                pt, lo, hi = record.interval(groups)
                self.assertEqual((lo, hi), synthetic.brute_interval(groups))
                self.assertEqual(pt, record.median(v for g in groups for v in g))

    def test_reading(self):
        self.assertEqual(record.reading(Fraction(1, 10), Fraction(2, 10), "log_loss"), "worse")
        self.assertEqual(record.reading(Fraction(1, 10), Fraction(2, 10), "recall"), "better")
        self.assertEqual(record.reading(Fraction(-2, 10), Fraction(-1, 10), "recall"), "worse")
        self.assertEqual(record.reading(Fraction(-2, 10), Fraction(-1, 10), "error_rate"), "better")
        # A bound at exactly 0 reads not shown.
        self.assertEqual(record.reading(Fraction(0), Fraction(1, 10), "log_loss"), "not shown")
        self.assertEqual(record.reading(Fraction(-1, 10), Fraction(0), "recall"), "not shown")

    def test_outcome_and_experiment_rule(self):
        self.assertEqual(record.outcome("worse", "not shown"), "cost")
        self.assertEqual(record.outcome("not shown", "worse"), "cost")
        self.assertEqual(record.outcome("worse", "better"), "mixed")
        self.assertEqual(record.outcome("not shown", "better"), "benefit")
        self.assertEqual(record.outcome("not shown", "not shown"), "not shown")
        self.assertEqual(record.outcome("worse", "not shown", both=True), "not shown")
        self.assertEqual(record.outcome("worse", "worse", both=True), "cost")
        self.assertEqual(record.experiment_rule(["cost"] * 3 + ["not shown"] * 2), "cost shown")
        self.assertEqual(record.experiment_rule(["cost"] * 2 + ["not shown"] * 3), "cost not shown")
        self.assertEqual(record.experiment_rule(["cost"] * 4 + ["mixed"]), "mixed")

    def test_mates_and_agreement(self):
        # Study a holds genomes 0 to 2 and b genome 3; genome 4 has none and so no mate.
        self.assertEqual(record.mates([0, 0, 1, 1, 1], ["a", "a", "a", "b", None]), [1, 1, 2, 0, 0])
        a = [0, 0, 1, 1, 2, 2, 3, 3, 4, 4]
        self.assertEqual(record.best_agreement(a, [(x + 2) % 5 for x in a]), 10)
        self.assertEqual(record.best_agreement(a, [(x + 2) % 5 for x in a[:-1]] + [0]), 9)

    def test_fmt_median_and_num(self):
        self.assertEqual([record.fmt(x) for x in (True, False, None, 7, Fraction(1, 8), 0.1)],
                         ["yes", "no", "-", "7", "0.125", "0.1"])
        self.assertEqual(record.median([4, 1, 3, 2]), Fraction(5, 2))
        self.assertEqual(record.num("0.125"), Fraction(1, 8))
        self.assertEqual(record.num("-12"), -12)

    def test_nostudy_on_the_same_datasets_by_hand(self):
        # Two datasets. In a, the genomes without a BioProject hold no resistant genome, so their recall stands in b
        # alone. With a BioProject, recall is 6/10 in study plans and 8/10 in matched plans in a, 5/10 against 9/10
        # in b. Study minus matched: -1/5 in a and -2/5 in b. Over both datasets the mean is -3/10; over the same
        # datasets as the genomes without a BioProject, b alone, it is -2/5.
        def r(learner, tp, fn, tn, fp):
            return {"learner": learner, "tp": str(tp), "fn": str(fn), "tn": str(tn), "fp": str(fp), "log_loss": "0.5"}
        sub = {}
        for d, study, matched, none in (("a", (6, 4), (8, 2), (0, 0)), ("b", (5, 5), (9, 1), (1, 1))):
            for kind in record.KINDS:
                for i in record.DRAWS:
                    for l in record.LEARNERS:
                        tp, fn = study if kind == "study" else matched
                        sub[(d, f"{kind}-d{i}", l, "study")] = r(l, tp, fn, 10, 0)
                        sub[(d, f"{kind}-d{i}", l, "no-study")] = r(l, *none, 3, 0)
        got = {(x[0], x[1], x[2], x[3]): x[4:6] for x in record.nostudy_rows(sub, ["a", "b"])}
        self.assertEqual(got[("rf", "SM", "recall", "no-study")], ["1", "0.0"])
        self.assertEqual(got[("rf", "SM", "recall", "study")], ["2", repr(-0.3)])
        self.assertEqual(got[("rf", "SM", "recall", "study-same")], ["1", repr(-0.4)])
        # Specificity stands in both datasets, so both sides count 2.
        self.assertEqual(got[("rf", "SM", "specificity", "study-same")], ["2", "0.0"])

    def test_rule_of(self):
        rule = {"0": "lane", "1": "B"}
        self.assertEqual([record.rule_of(p, rule) for p in ("random-d1", "study-d0", "matched-d1")], ["-", "A", "B"])


class SyntheticRecord(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.addClassCleanup(cls.tmp.cleanup)
        cls.root = cls.tmp.name
        cls.inputs = synthetic.make(os.path.join(cls.root, "in"))
        cls.out = os.path.join(cls.root, "out")
        rc, err = run(cls.inputs, "--out", cls.out)
        assert rc == 0, err
        cls.fits = rows(cls.out, "fits.tsv", "dataset", "plan", "learner")
        cls.summary = {r["key"]: r["value"] for r in read_tsv(os.path.join(cls.out, "summary.tsv"))}

    def test_fits_by_hand(self):
        # SYNTHETIC D1: 19 resistant and 21 susceptible genomes; PRJSYN01, its first BioProject, holds 2 and 3.
        f = self.fits
        self.assertEqual([f[(D1, "random-d0", "scm")][k] for k in ("tp", "fn", "tn", "fp", "error_rate")],
                         ["18", "1", "21", "0", "0.025"])
        self.assertEqual([f[(D1, "study-d3", "scm")][k] for k in ("tp", "fn", "tn", "fp", "recall", "error_rate")],
                         ["17", "2", "18", "3", repr(17 / 19), "0.125"])
        self.assertAlmostEqual(float(f[(D1, "random-d0", "rf")]["log_loss"]), -math.log(0.8), places=14)
        self.assertAlmostEqual(float(f[(D1, "study-d0", "rf")]["log_loss"]),
                               (35 * -math.log(0.8) + 5 * -math.log(0.2)) / 40, places=14)
        self.assertAlmostEqual(float(f[(D1, "matched-d2", "scm")]["log_loss"]),
                               (-math.log(1e-4) + 39 * -math.log(1 - 1e-4)) / 40, places=14)
        lin = rows(self.out, "lineage.tsv", "dataset", "learner")
        self.assertEqual([lin[(D1, "cart")][k] for k in ("tp", "fn", "tn", "fp")], ["17", "2", "21", "0"])
        self.assertNotIn((D4, "cart"), lin)

    def test_contrast_and_interval_by_hand(self):
        # scm's error rate, study minus random: 1/10 in D1, D3 and D4 (5 errors against 1 of 40) and 1/15 in D2
        # (5 against 1 of 60). The median is 1/10; only the draw of alphabacter three times gives the median
        # (1/15 + 1/10) / 2 = 1/12, 1 of 27, which is above 1/40, so the interval is [1/12, 1/10].
        c = rows(self.out, "contrasts.tsv", "experiment", "learner", "family", "instance", "measure")
        for inst in ("0,0", "3,1", "mean"):
            r = c[("second", "scm", "SR", inst, "error_rate")]
            self.assertEqual([r[k] for k in ("median", "low", "high", "reading", "worse")],
                             ["0.1", repr(1 / 12), "0.1", "worse", "4"])
        r = c[("second", "scm", "MR", "4,2", "error_rate")]
        self.assertEqual([r[k] for k in ("median", "low", "high", "reading", "equal")], ["0.0", "0.0", "0.0",
                                                                                         "not shown", "4"])
        # Recall, study minus random: -1/R in each dataset, with R 19, 27, 20 and 20.
        r = c[("second", "scm", "SR", "2,2", "recall")]
        self.assertEqual([r[k] for k in ("median", "low", "high", "reading")],
                         ["-0.05", "-0.05", repr(float(-(Fraction(1, 19) + Fraction(1, 27)) / 2)), "worse"])
        rules = rows(self.out, "rules.tsv", "learner", "family", "rule")
        self.assertEqual((rules[("scm", "SR", "either")]["cost"], rules[("scm", "SR", "either")]["answer"]),
                         ("25", "shown"))
        self.assertEqual((rules[("scm", "MR", "either")]["cost"], rules[("scm", "MR", "either")]["answer"]),
                         ("0", "not shown"))

    def test_no_study_and_parts_by_hand(self):
        # D2 alone holds 10 genomes without a BioProject (4 resistant); no plan gets one of them wrong. With one,
        # scm's error rate is 5 of 50 in study plans against 1 of 50 in matched plans: 0.08.
        ns = rows(self.out, "nostudy.tsv", "learner", "family", "measure", "subset")
        self.assertEqual([ns[("scm", "SM", "loss", s)]["mean"] for s in ("no-study", "study")], ["0.0", "0.08"])
        self.assertEqual(ns[("scm", "SM", "loss", "no-study")]["datasets"], "1")
        sub = rows(self.out, "subsets.tsv", "dataset", "plan", "learner", "subset")
        self.assertEqual([sub[(D2, "study-d0", "scm", "no-study")][k] for k in ("n", "resistant", "fn", "fp")],
                         ["10", "4", "0", "0"])
        parts = rows(self.out, "parts.tsv", "learner", "measure")
        self.assertEqual([parts[("scm", "loss")][k] for k in ("mean_mr", "share", "larger", "exact")],
                         ["0.0", "1.0", "4", "yes"])

    def test_counts_by_hand(self):
        s = self.summary
        self.assertEqual([s[k] for k in ("datasets", "species", "pairs", "no_study_pairs", "multisets",
                                         "release_datasets", "eligible_counts", "fail_studies", "fail_share",
                                         "fail_group", "distinct_genomes", "distinct_studies", "two_large")],
                         ["4", "3", "180", "10", "10", "8", "5", "1", "1", "2", "160", "30",
                          "alphabacter_primus,betabacter_secundus"])
        self.assertEqual(s["left_by_folds"], "deltabacter_quartus__onemycin")
        # Rule A leaves D1 unchanged and moves i of the others' genomes in draw i; rule B moves 2i + 4.
        self.assertEqual((s["rule_A_unchanged"], s["rule_B_max"]), ("1", repr(float(Fraction(9, 40)))))
        el = rows(self.out, "eligibility.tsv", "dataset")
        self.assertEqual([el[("betabacter_secundus__twomycin",)][k] for k in ("under_share", "over_group", "passes")],
                         ["yes", "yes", "no"])
        self.assertEqual(el[("deltabacter_quartus__onemycin",)]["folds"], "fold 2 test genomes lack a class")
        var = rows(self.out, "variants.tsv", "min_studies", "min_study_share", "max_group_share")
        self.assertEqual([var[("5", "0.8", "0.5")][k] for k in ("datasets", "species", "pairs")], ["5", "4", "210"])
        ck = {r["check"]: r["count"] for r in read_tsv(os.path.join(self.out, "checks.tsv"))}
        self.assertEqual([ck[k] for k in ("fold_files_checked", "learner_fits_equal_evaluator_counts",
                                          "matched_plans_equal_study_counts", "hash_list_entries_equal_inputs",
                                          "inputs_synthetic_inputs", "genome_ids_in_record")],
                         ["60", "315", "20", "4", "1", "0"])
        ds = rows(self.out, "datasets.tsv", "dataset")
        self.assertEqual([ds[(D2,)][k] for k in ("studies", "no_study", "largest_study_share", "lineage_split")],
                         ["10", "10", repr(1 / 12), "scored"])
        self.assertEqual(ds[(D4,)]["lineage_split"], "unavailable")
        plans = rows(self.out, "plans.tsv", "dataset", "plan")
        self.assertEqual([plans[(D1, f"study-d{i}")]["study_mates"] for i in range(5)], ["0"] * 5)
        self.assertEqual([plans[(D1, f"study-d{i}")]["rule"] for i in range(5)], ["A", "B", "B", "B", "B"])

    def test_order_by_hand(self):
        o = {(r["experiment"], r["event"]): r["time"] for r in read_tsv(os.path.join(self.out, "order.tsv"))}
        self.assertEqual(o[("second", "stated_rule_choice")], f"{synthetic.DAY}T{synthetic.FIXED}")
        self.assertEqual(o[("first", "hash_time")], f"{synthetic.DAY}T09:40:00Z")
        self.assertEqual((o[("first", "first_fit_end")], o[("first", "last_fit_end")]),
                         (f"{synthetic.DAY}T10:00:00Z", f"{synthetic.DAY}T10:07:00Z"))
        self.assertEqual((o[("second", "fits_start")], o[("second", "last_fit_end")]),
                         (f"{synthetic.DAY}T11:40:00Z", f"{synthetic.DAY}T11:41:00Z"))
        src = {(r["experiment"], r["event"]): r["source"] for r in read_tsv(os.path.join(self.out, "order.tsv"))}
        plan = pathlib.Path(self.inputs, "second", "PLAN.md").read_bytes()
        self.assertEqual(src[("second", "stated_rule_choice")],
                         "text of the second plan, sha256 " + synthetic.sha(plan))
        self.assertEqual(src[("first", "fold_builder")], "file time in the input set's first/times.tsv")
        self.assertEqual(src[("second", "hash_time")], "first line of the input set's second/PLAN.sha256")
        pf = rows(self.out, "planfiles.tsv", "experiment")
        self.assertEqual([pf[(e,)]["file"] for e in ("first", "second")], ["first/PLAN.md", "second/PLAN.md"])
        st = {r["key"]: r["value"] for r in read_tsv(os.path.join(self.out, "settings.tsv"))}
        self.assertEqual(st["first_plan_resamples"], "10000")

    def test_cells_agree(self):
        self.assertEqual(agree.problems(self.out), [])

    def test_holds_no_genome_id(self):
        text = "".join(v.decode("utf-8") for k, v in files(self.out).items() if k != "inputs.sha256")
        self.assertNotIn("9001.0", text)
        self.assertNotIn("PRJSYN", text)

    def test_inputs_named_relative(self):
        text = pathlib.Path(self.out, "inputs.sha256").read_text(encoding="utf-8")
        self.assertNotIn(self.root, text)
        self.assertEqual(len(text.splitlines()), len(files(self.inputs)))

    def test_same_bytes_on_a_second_run(self):
        again = os.path.join(self.root, "again")
        rc, err = run(self.inputs, "--out", again, "--jobs", "2")
        self.assertEqual(rc, 0, err)
        self.assertEqual(files(again), files(self.out))

    def test_same_bytes_at_other_hash_seeds(self):
        # Without -I, which would drop PYTHONHASHSEED; -s and -P keep site-packages and the cwd out.
        seen = set()
        for seed in ("0", "4242"):
            env = {"PYTHONHASHSEED": seed, "PATH": "/usr/bin:/bin"}
            out = os.path.join(self.root, "seed" + seed)
            subprocess.run([sys.executable, "-s", "-P", "-B", str(SCRIPT), "--inputs", self.inputs, "--out", out],
                           env=env, check=True, capture_output=True)
            self.assertEqual(files(out), files(self.out), seed)
            seen.add(subprocess.run([sys.executable, "-s", "-P", "-c", "print(hash('study split'))"], env=env,
                                    check=True, capture_output=True, text=True).stdout)
        self.assertEqual(len(seen), 2, "the two seeds gave one hash")

    def test_check(self):
        self.assertEqual(run(self.inputs, "--check", self.out), (0, ""))
        work = os.path.join(self.root, "check")
        self.assertEqual(run(self.inputs, "--out", work)[0], 0)
        path = pathlib.Path(work, "totals.tsv")
        path.write_text(path.read_text(encoding="utf-8").replace("\t160\n", "\t161\n"), encoding="utf-8")
        pathlib.Path(work, "variation.tsv").unlink()
        pathlib.Path(work, "extra.tsv").write_text("x\n", encoding="utf-8")
        rc, err = run(self.inputs, "--check", work)
        self.assertEqual(rc, 1)
        self.assertIn("totals.tsv: differs", err)
        self.assertIn("variation.tsv: missing", err)
        self.assertIn("extra.tsv: not written by the script", err)
        self.assertEqual(run(self.inputs, "--check", os.path.join(self.root, "none"))[0], 1)

    def test_out_removes_stale_files(self):
        # Names the earlier SHA256SUMS lists and this build does not write go, with a directory left empty.
        work = os.path.join(self.root, "stale")
        self.assertEqual(run(self.inputs, "--out", work)[0], 0)
        sums = pathlib.Path(work, "SHA256SUMS")
        for name in ("old.tsv", "old/old-PLAN.md"):
            pathlib.Path(work, name).parent.mkdir(exist_ok=True)
            pathlib.Path(work, name).write_text("x\n", encoding="utf-8")
            sums.write_text(sums.read_text(encoding="utf-8") + "0" * 64 + f"  {name}\n", encoding="utf-8")
        self.assertEqual(run(self.inputs, "--out", work)[0], 0)
        self.assertEqual(files(work), files(self.out))
        self.assertFalse(pathlib.Path(work, "old").exists())

    def test_out_leaves_a_file_no_sha256sums_lists(self):
        work = os.path.join(self.root, "unlisted")
        self.assertEqual(run(self.inputs, "--out", work)[0], 0)
        pathlib.Path(work, "notes.txt").write_text("kept\n", encoding="utf-8")
        pathlib.Path(work, "sub").mkdir()
        pathlib.Path(work, "sub", "draft.tsv").write_text("kept too\n", encoding="utf-8")
        self.assertEqual(run(self.inputs, "--out", work)[0], 0)
        have = files(work)
        self.assertEqual((have.pop("notes.txt"), have.pop("sub/draft.tsv")), (b"kept\n", b"kept too\n"))
        self.assertEqual(have, files(self.out))

    def test_out_refuses_a_target_without_sha256sums(self):
        work = pathlib.Path(self.root, "other")
        (work / "deep").mkdir(parents=True)
        (work / "keep.txt").write_text("keep\n", encoding="utf-8")
        (work / "deep" / "README.md").write_text("another record\n", encoding="utf-8")
        before = files(work)
        rc, err = run(self.inputs, "--out", str(work))
        self.assertEqual(rc, 2)
        self.assertIn("holds no SHA256SUMS", err)
        self.assertEqual(files(work), before)
        # A SHA256SUMS that names a file outside the record stops it too, before anything is written.
        (work / "SHA256SUMS").write_text("0" * 64 + "  ../outside.tsv\n", encoding="utf-8")
        rc, err = run(self.inputs, "--out", str(work))
        self.assertEqual(rc, 2)
        self.assertIn("not a name inside the record", err)
        self.assertEqual(set(files(work)), set(before) | {"SHA256SUMS"})

    def test_out_refuses_another_record(self):
        # Another record holds a SHA256SUMS of its own; the script must not take it for an earlier record of its own.
        source = ROOT / "reference" / "other-database"
        self.assertTrue((source / "SHA256SUMS").is_file())
        copy = pathlib.Path(self.root, "other-database")
        shutil.copytree(source, copy)
        before = files(copy)
        rc, err = run(self.inputs, "--out", str(copy))
        self.assertEqual(rc, 2)
        self.assertIn("lists no README.md titled", err)
        self.assertEqual(files(copy), before)
        shutil.rmtree(copy)
        # An earlier record of this script whose README.md lost its title is refused the same way.
        work = pathlib.Path(self.root, "retitled")
        self.assertEqual(run(self.inputs, "--out", str(work))[0], 0)
        readme = work / "README.md"
        readme.write_text("# Another record\n" + readme.read_text(encoding="utf-8").split("\n", 1)[1], encoding="utf-8")
        (work / "stale.tsv").write_text("x\n", encoding="utf-8")
        sums = work / "SHA256SUMS"
        sums.write_text(sums.read_text(encoding="utf-8") + "0" * 64 + "  stale.tsv\n", encoding="utf-8")
        before = files(work)
        rc, err = run(self.inputs, "--out", str(work))
        self.assertEqual(rc, 2)
        self.assertIn("lists no README.md titled", err)
        self.assertEqual(files(work), before)

    def test_out_one_level_too_high_leaves_the_other_records(self):
        # The repository's reference/ holds this record beside the others, and no SHA256SUMS of its own.
        source = ROOT / "reference"
        self.assertFalse((source / "SHA256SUMS").exists())
        others = [p.name for p in source.iterdir() if p.name != "study-split"]
        self.assertGreaterEqual(len(others), 13)
        copy = pathlib.Path(self.root, "reference")
        shutil.copytree(source, copy)
        before = files(copy)
        rc, err = run(self.inputs, "--out", str(copy))
        self.assertEqual(rc, 2)
        self.assertIn("holds no SHA256SUMS", err)
        self.assertEqual(files(copy), before)
        shutil.rmtree(copy)


class Stops(unittest.TestCase):
    """Each SYNTHETIC fault stops the script with exit 2 and writes nothing."""

    def test_faults_stop(self):
        for fault, message in synthetic.FAULTS.items():
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as tmp:
                inputs = synthetic.make(os.path.join(tmp, "in"), fault=fault)
                out = os.path.join(tmp, "out")
                rc, err = run(inputs, "--out", out)
                self.assertEqual(rc, 2)
                self.assertIn("stop:", err)
                self.assertIn(message, err)
                self.assertFalse(os.path.exists(out))

    def unexpected(self, fault):
        """--out and --check on inputs with one thing the script does not expect: exit 2, the thing named."""
        with tempfile.TemporaryDirectory() as tmp:
            clean = synthetic.make(os.path.join(tmp, "clean"))
            record_dir = os.path.join(tmp, "record")
            self.assertEqual(run(clean, "--out", record_dir)[0], 0)
            before = files(record_dir)
            inputs = synthetic.make(os.path.join(tmp, "in"), fault=fault)
            for mode, target in (("--out", os.path.join(tmp, "out")), ("--check", record_dir)):
                rc, err = run(inputs, mode, target)
                self.assertEqual(rc, 2, mode)
                self.assertIn(synthetic.FAULTS[fault], err, mode)
            self.assertFalse(os.path.exists(os.path.join(tmp, "out")))
            self.assertEqual(files(record_dir), before)

    def test_an_unexpected_input_file_stops(self):
        self.unexpected("extra_input")

    def test_an_unexpected_fit_file_stops(self):
        self.unexpected("extra_fit")

    def test_an_unexpected_column_stops(self):
        self.unexpected("extra_column")

    def test_an_unexpected_dataset_table_stops(self):
        self.unexpected("extra_dataset")

    def test_jobs_and_missing_inputs_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(run(tmp, "--jobs", "13", "--out", os.path.join(tmp, "o"))[0], 2)
            rc, err = run(tmp, "--out", os.path.join(tmp, "o"))
            self.assertEqual(rc, 2)
            self.assertIn("checks.tsv missing", err)
            self.assertEqual(run(os.path.join(tmp, "none"), "--out", os.path.join(tmp, "o"))[0], 2)


if __name__ == "__main__":
    unittest.main()
