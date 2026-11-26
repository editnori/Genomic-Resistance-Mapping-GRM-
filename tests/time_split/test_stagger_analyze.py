"""scripts/time-split/stagger_analyze.py on SYNTHETIC designs and result files: each stop, the mean over a dataset's
pairs, the median over datasets, the interval's two ends and its seed, each reading and verdict, the second
reading's three cases, the list files, stale same size against recent(j), the scan of what it writes, the
predictions and their outcome, and whether a verdict stands, as stands.tsv and predictions.tsv both give it."""

import json
import math
import random
import shutil
import statistics
import tempfile
import unittest
import unittest.mock
from pathlib import Path

from tests.time_split import stagger_fixture as fx

an = fx.analyze_script()
sd = fx.design_script()


def samples_of(tag, r, s):
    """r resistant then s susceptible made-up test genomes."""
    return [(f"syn-{tag}-r{i:02d}", True) for i in range(r)] + [(f"syn-{tag}-s{i:02d}", False) for i in range(s)]


def write_case(base, datasets, clock="public-day", minimum=10, lists=None, shares=None):
    """A hand design and its result files. Each dataset: id, species, tests {window: samples}, scores
    {(window, origin): one score list for every learner, or {learner: scores}}, pairs [(older, window)]; and, for
    stale same size, recent {(window, origin of 2 or more): scores} and same_size_pairs [(older, window)], with
    recent(1) in the file of all(1), as the design writes it. shares: {clock: share set}, 40, 60 and 80 where absent.
    Returns (design dir, result dir, list dir or None)."""
    base = Path(base)
    design_dir, results = base / "design", base / "results"
    design_dir.mkdir(parents=True)
    results.mkdir()
    design = {"plan_sha256": sd.PLAN_SHA256, "key": sd.KEY, "min_train_class": 25,
              "clocks": {c: {"shares": list((shares or {}).get(c, (40, 60, 80))), "test_minimum": minimum,
                             "descriptive_only": False} for c in sd.CLOCKS}, "datasets": []}
    for d in datasets:
        row = {"id": d["id"], "species_slug": d["species"], "antibiotic": "drug", "clock": clock, "built": True,
               "why_not": "", "windows": [], "pairs": {c: [] for c in sd.CONTRASTS}}
        for j, samples in sorted(d["tests"].items()):
            w = {"window": j, "test": len(samples), "test_resistant": sum(t for _, t in samples), "copies": 0,
                 "fits": []}
            for k in range(1, j + 1):
                if (j, k) not in d["scores"]:
                    continue
                scheme = f"stagger-{clock}-w{j}-o{k}-all"
                fit = {"origin": k, "arm": "all", "train": 60 + k, "train_resistant": 30, "scheme": scheme,
                       "file": f"{d['id']}.{scheme}.json", "sha256": fx.sha(d["id"] + scheme)}
                scores = d["scores"][(j, k)]
                by_learner = scores if isinstance(scores, dict) else {l: scores for l in fx.LEARNERS}
                (results / fit["file"]).write_text(json.dumps(fx.result(d["id"], scheme, fit["sha256"], samples,
                                                                        by_learner, fit["train"], 30)))
                w["fits"].append(fit)
            recent = {k: s for (jj, k), s in d.get("recent", {}).items() if jj == j}
            if recent:
                shared = next(f for f in w["fits"] if f["origin"] == 1)
                w["fits"].append({**shared, "arm": "recent", "same_file_as": "all"})
            for k, scores in sorted(recent.items()):
                scheme = f"stagger-{clock}-w{j}-o{k}-recent"
                fit = {"origin": k, "arm": "recent", "train": 50 + k, "train_resistant": 25, "scheme": scheme,
                       "file": f"{d['id']}.{scheme}.json", "sha256": fx.sha(d["id"] + scheme)}
                (results / fit["file"]).write_text(json.dumps(fx.result(d["id"], scheme, fit["sha256"], samples,
                                                                        {l: scores for l in fx.LEARNERS},
                                                                        fit["train"], 25)))
                w["fits"].append(fit)
            row["windows"].append(w)
        row["pairs"]["stale-all"] = [{"older": k, "fresh": j, "window": j, "horizon": j - k} for k, j in d["pairs"]]
        row["pairs"]["stale-same-size"] = [{"older": k, "fresh": j, "window": j, "horizon": j - k}
                                           for k, j in d.get("same_size_pairs", [])]
        design["datasets"].append(row)
    (design_dir / "design.json").write_text(json.dumps(design))
    list_dir = None
    if lists is not None:
        list_dir = base / "lists"
        list_dir.mkdir()
        for name, gids in lists.items():
            (list_dir / f"{name}.txt").write_text("".join(g + "\n" for g in gids))
    return design_dir, results, list_dir


def scores_for(samples, missed_r=(), called_s=()):
    """0.9 for a resistant genome and 0.1 for a susceptible one, but 0.1 for the R ids in missed_r and 0.9 for the
    S ids in called_s."""
    return [(0.1 if g in missed_r else 0.9) if t else (0.9 if g in called_s else 0.1) for g, t in samples]


def table(path):
    lines = Path(path).read_text().splitlines()
    head = lines[0].split("\t")
    return [dict(zip(head, line.split("\t"))) for line in lines[1:]]


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp(prefix="stagger-an-"))

    def tearDown(self):
        shutil.rmtree(self.tmp)

    def one_pair(self, lists=None):
        """One dataset, one pair on window 2: 20 R and 20 S; the older fit misses 5 R, the fresh fit none."""
        samples = samples_of("p", 20, 20)
        missed = {g for g, t in samples if t and g.endswith(("r00", "r01", "r02", "r03", "r04"))}
        return write_case(self.tmp, [{"id": "hand_species__drug", "species": "hand_species", "tests": {2: samples},
                                      "scores": {(2, 1): scores_for(samples, missed), (2, 2): scores_for(samples)},
                                      "pairs": [(1, 2)]}], lists=lists), samples, missed

    def run_it(self, design_dir, results, list_dir=None, out="out"):
        out = self.tmp / out
        argv = [design_dir, results, out] + (["--listed", list_dir] if list_dir else [])
        code, text, err = fx.run(an, *argv)
        return code, text, err, out


class Stops(Case):
    def stops(self, change, message):
        (design_dir, results, _), samples, _ = self.one_pair()
        change(results, samples)
        code, text, err, out = self.run_it(design_dir, results)
        self.assertEqual(code, 1)
        self.assertEqual(text, "")
        self.assertIn(message, err)
        self.assertFalse(out.exists())

    def edit(self, results, origin, fn):
        path = results / f"hand_species__drug.stagger-public-day-w2-o{origin}-all.json"
        r = json.loads(path.read_text())
        r = fn(r) or r
        path.write_text(json.dumps(r))

    def test_the_case_runs_without_a_change(self):
        (design_dir, results, _), _, _ = self.one_pair()
        code, _, err, out = self.run_it(design_dir, results)
        self.assertEqual(code, 0, err)
        self.assertTrue((out / "verdicts.tsv").exists())

    def test_a_result_that_does_not_name_its_fold_files_sha256(self):
        def change(results, _):
            self.edit(results, 1, lambda r: r.update(fold_file_sha256=fx.sha("another fold file")))
        self.stops(change, "does not name its fold file's sha256")

    def test_two_fits_that_scored_other_genomes(self):
        def change(results, samples):
            other = [("syn-p-other", True)] + samples[1:]
            self.edit(results, 1, lambda r: r.update(samples=[{"genome_id": g, "resistant": t, "fold": 0}
                                                              for g, t in other]))
        self.stops(change, "the two fits scored other genomes")

    def test_the_same_genomes_with_other_labels(self):
        def change(results, samples):
            # The first R and the first S swap labels: the counts stay, the labels do not.
            swapped = [(g, not t) if i in (0, 20) else (g, t) for i, (g, t) in enumerate(samples)]
            path = results / "hand_species__drug.stagger-public-day-w2-o1-all.json"
            old = json.loads(path.read_text())
            scores = {x["learner"]: x["scores"] for x in old["split"]["learners"]}
            path.write_text(json.dumps(fx.result(old["id"], old["scheme"], old["fold_file_sha256"], swapped, scores,
                                                 old["folds"][0]["train"], old["folds"][0]["train_resistant"])))
        self.stops(change, "the same genomes, other labels")

    def test_training_counts_other_than_the_designs(self):
        def change(results, _):
            self.edit(results, 2, lambda r: r["folds"][0].update(train=r["folds"][0]["train"] - 1))
        self.stops(change, "training counts other than the design's")

    def test_recall_worked_out_again_other_than_the_evaluators(self):
        def change(results, _):
            def fn(r):
                m = r["split"]["learners"][2]["metrics"]
                m["sensitivity"] = m["sensitivity"] - 0.01
            self.edit(results, 1, fn)
        self.stops(change, "rf: recall worked out again is not the evaluator's")

    def test_a_missing_result_file(self):
        self.stops(lambda results, _: (results / "hand_species__drug.stagger-public-day-w2-o2-all.json").unlink(),
                   "no result file")

    def test_a_barred_result_path(self):
        (design_dir, results, _), _, _ = self.one_pair()
        code, _, err, _ = self.run_it(design_dir, self.tmp / "held-out-results")
        self.assertEqual((code, err), (1, "stagger_analyze: stopped: the result dir path holds 'held-out'; nothing "
                                          "under it is opened\n"))


class Values(Case):
    def test_a_dataset_with_two_pairs_at_a_horizon_gives_their_mean(self):
        w2, w3 = samples_of("m", 2, 2), samples_of("n", 2, 2)
        # Window 2: the older fit finds 1 R of 2, the fresh fit 2 of 2, a difference of -0.5. Window 3: no difference.
        write_case(self.tmp, [{"id": "mean_species__drug", "species": "mean_species", "tests": {2: w2, 3: w3},
                               "scores": {(2, 1): scores_for(w2, {"syn-m-r01"}), (2, 2): scores_for(w2),
                                          (3, 2): scores_for(w3), (3, 3): scores_for(w3)},
                               "pairs": [(1, 2), (2, 3)]}])
        code, _, err, out = self.run_it(self.tmp / "design", self.tmp / "results")
        self.assertEqual(code, 0, err)
        rows = [r for r in table(out / "contrasts.tsv") if r["horizon"] == "1"]
        self.assertEqual(len(rows), 5)
        self.assertEqual({(r["pairs"], r["test"], r["d_recall"], r["d_specificity"]) for r in rows},
                         {("2", "8", "-0.250000", "0.000000")})
        older = [r for r in table(out / "arms.tsv") if r["horizon"] == "1" and r["arm"] == "older"]
        self.assertEqual({r["recall"] for r in older}, {"0.750000"})

    def test_measures_at_the_cutoff_and_the_clip(self):
        samples = [("syn-q-1", True), ("syn-q-2", True), ("syn-q-3", False), ("syn-q-4", False)]
        got = an.measure(samples, [0.5, 0.0, 0.49, 1.0])
        self.assertEqual((got["tp"], got["fn"], got["tn"], got["fp"]), (1, 1, 1, 1))
        # Log loss: 0.0 and 1.0 are clipped to 1e-4 and 1 - 1e-4.
        want = -(math.log(0.5) + math.log(1e-4) + math.log(0.51) + math.log(1e-4)) / 4
        self.assertAlmostEqual(got["log_loss"], want, places=12)
        kept = an.measure(samples, [0.5, 0.0, 0.49, 1.0], keep={0, 2})
        self.assertEqual((kept["recall"], kept["specificity"]), (1.0, 1.0))


class Scripted:
    """A generator whose draws are written out: 250 of species a, 1 of b, 9,498 of c, 1 of d and 250 of e."""

    def __init__(self):
        self.plan = ["a"] * 250 + ["b"] + ["c"] * 9498 + ["d"] + ["e"] * 250

    def choices(self, population, k):
        return [self.plan.pop()] * k


class Interval(unittest.TestCase):
    def test_the_251st_and_the_9750th_of_the_sorted_draws(self):
        by_species = {"a": [0.0], "b": [1.0], "c": [2.0], "d": [3.0], "e": [4.0]}
        rng = Scripted()
        self.assertEqual(an.interval(by_species, rng), (2.0, 1.0, 3.0))
        self.assertEqual(rng.plan, [])

    def test_the_seed_and_the_order_of_draws(self):
        rows = []
        for i in range(14):
            d = int(fx.sha(f"value {i}")[:6], 16) / 16 ** 6 - 0.5
            rows.append({"species": ("sp_a", "sp_b", "sp_c")[i % 3], "dataset": f"ds_{i}",
                         "learners": {l: {"d": {m: d + n / 100 for n, m in enumerate(an.MEASURES)}}
                                      for l in fx.LEARNERS}})
        design = {"clocks": {c: {"shares": [40, 60, 80], "descriptive_only": False} for c in sd.CLOCKS}}
        made = []

        class Recording(random.Random):
            def __init__(self, seed):
                made.append(seed)
                super().__init__(seed)

        with unittest.mock.patch.object(an.random, "Random", Recording):
            reading_rows = an.analyze(design, {("public-day", "first", "stale-all", 1): rows}, ("first",))[0]
        # One generator a clock and reading, each seeded 20261125.
        self.assertEqual(made, [20261125, 20261125])
        # Worked out here on its own: one generator, recall then log loss for the first learner.
        rng = random.Random(20261125)
        for row, (n, m) in zip(reading_rows[:2], enumerate(an.MEASURES)):
            by = {}
            for r in rows:
                by.setdefault(r["species"], []).append(r["learners"]["scm"]["d"][m])
            names = sorted(by)
            draws = sorted(statistics.median([v for s in rng.choices(names, k=3) for v in by[s]])
                           for _ in range(10000))
            self.assertEqual((row["learner"], row["measure"]), ("scm", m))
            self.assertEqual((row["low"], row["high"]), (f"{draws[250]:.6f}", f"{draws[9749]:.6f}"))

    def test_the_point_is_the_median_over_datasets(self):
        # Three datasets in two species: the median of 0.0, 0.1 and 0.5 is 0.1; their mean is 0.2, and the median
        # of the two species' medians is 0.275.
        point = an.interval({"a": [0.0, 0.1], "b": [0.5]}, random.Random(an.SEED))[0]
        self.assertEqual(point, 0.1)

    def test_the_point_of_an_even_count_of_datasets_is_the_mean_of_the_two_middle_values(self):
        # Four datasets in two species, 0.0, 0.25, 0.75 and 1.0: the median is 0.5. The lower of the two middle
        # values is 0.25 and the higher 0.75.
        point = an.interval({"a": [0.0, 0.75], "b": [0.25, 1.0]}, random.Random(an.SEED))[0]
        self.assertEqual(point, 0.5)


class DescriptiveClock(unittest.TestCase):
    def test_a_clock_the_rule_marks_descriptive_gives_descriptive_only(self):
        rows = [{"species": ("sp_a", "sp_b", "sp_c")[i % 3], "dataset": f"ds_{i}",
                 "learners": {l: {"d": {m: 0.2 * an.WORSE[m] for m in an.MEASURES}} for l in fx.LEARNERS}}
                for i in range(14)]
        values = {("public-day", "first", "stale-all", 1): rows}
        for descriptive, want in ((False, "cost shown"), (True, "descriptive only")):
            design = {"clocks": {c: {"shares": [50, 70], "descriptive_only": descriptive} for c in sd.CLOCKS}}
            verdicts = an.analyze(design, values, ("first",))[3]
            self.assertEqual(verdicts[("first", "public-day", "stale-all", 1)], want)


class Readings(unittest.TestCase):
    def test_each_reading_of_an_interval(self):
        cases = [((-0.3, -0.1, "recall"), "worse"), ((0.1, 0.3, "recall"), "better"),
                 ((-0.1, 0.1, "recall"), "not shown"), ((0.0, 0.2, "recall"), "not shown"),
                 ((0.1, 0.3, "log_loss"), "worse"), ((-0.3, -0.1, "log_loss"), "better"),
                 ((-0.2, 0.0, "log_loss"), "not shown"), ((-0.3, -0.1, "specificity"), "worse")]
        for args, want in cases:
            self.assertEqual(an.read(*args), want, args)

    def test_each_rule_for_a_learner(self):
        cases = [(("worse", "not shown"), "cost"), (("not shown", "worse"), "cost"), (("worse", "worse"), "cost"),
                 (("better", "not shown"), "benefit"), (("not shown", "better"), "benefit"),
                 (("better", "worse"), "mixed"), (("not shown", "not shown"), "not shown")]
        for args, want in cases:
            self.assertEqual(an.learner_rule(*args), want, args)

    def rules(self, costs, better=0):
        rules, reads = {}, {}
        for i, l in enumerate(fx.LEARNERS):
            if i < better:
                reads[l] = ("better", "not shown")
            elif i < better + costs:
                reads[l] = ("worse", "not shown")
            else:
                reads[l] = ("not shown", "not shown")
            rules[l] = an.learner_rule(*reads[l])
        return rules, reads

    def test_each_verdict(self):
        self.assertEqual(an.verdict(12, 3, *self.rules(3), False), "cost shown")
        self.assertEqual(an.verdict(12, 3, *self.rules(2), False), "not shown")
        self.assertEqual(an.verdict(12, 3, *self.rules(4, better=1), False), "mixed")
        self.assertEqual(an.verdict(12, 3, *self.rules(0), False), "not shown")
        self.assertEqual(an.verdict(11, 3, *self.rules(5), False), "descriptive only")
        self.assertEqual(an.verdict(30, 2, *self.rules(5), False), "descriptive only")
        self.assertEqual(an.verdict(30, 5, *self.rules(5), True), "descriptive only")

    def test_not_shown_is_never_written_as_equal(self):
        for args in [(-0.1, 0.1, m) for m in an.MEASURES]:
            self.assertNotEqual(an.read(*args), "equal")
        source = fx.ANALYZE.read_text()
        self.assertNotIn('"equal"', source)


class SecondCase(unittest.TestCase):
    samples = samples_of("c", 20, 20)

    def listed(self, r, s):
        return set([g for g, t in self.samples if t][:r] + [g for g, t in self.samples if not t][:s])

    def test_no_genome_listed(self):
        keep, why = an.second_case(self.samples, {"syn-c-elsewhere"}, 10)
        self.assertEqual((len(keep), why), (40, ""))

    def test_one_to_four_of_a_class(self):
        for r, s in ((1, 0), (4, 5), (5, 4), (4, 4), (6, 0)):
            self.assertEqual(an.second_case(self.samples, self.listed(r, s), 10), (None, an.SHORT), (r, s))

    def test_five_or_more_of_each_class(self):
        keep, why = an.second_case(self.samples, self.listed(5, 5), 10)
        self.assertEqual((len(keep), why), (30, ""))
        self.assertFalse({i for i, (g, _) in enumerate(self.samples) if g in self.listed(5, 5)} & keep)

    def test_five_or_more_but_a_class_falls_under_the_minimum(self):
        self.assertEqual(an.second_case(self.samples, self.listed(5, 5), 16), (None, an.UNDER))
        self.assertEqual(an.second_case(self.samples, self.listed(5, 5), 15)[1], "")

    def test_the_resistant_class_alone_falls_under_the_minimum_once_out(self):
        # 20 R and 30 S, 5 of each listed: 15 R and 25 S are left, under a minimum of 16 in R alone.
        samples = samples_of("u", 20, 30)
        listed = set([g for g, t in samples if t][:5] + [g for g, t in samples if not t][:5])
        self.assertEqual(an.second_case(samples, listed, 16), (None, an.UNDER))
        keep, why = an.second_case(samples, listed, 15)
        self.assertEqual((len(keep), why), (40, ""))

    def test_the_susceptible_class_alone_falls_under_the_minimum_once_out(self):
        # 30 R and 20 S, 5 of each listed: 25 R and 15 S are left, under a minimum of 16 in S alone.
        samples = samples_of("v", 30, 20)
        listed = set([g for g, t in samples if t][:5] + [g for g, t in samples if not t][:5])
        self.assertEqual(an.second_case(samples, listed, 16), (None, an.UNDER))
        keep, why = an.second_case(samples, listed, 15)
        self.assertEqual((len(keep), why), (40, ""))


class SecondReading(Case):
    def test_the_listed_genomes_leave_both_fits(self):
        samples = samples_of("p", 20, 20)
        missed = [g for g, t in samples if t][:5]
        listed = missed + [g for g, t in samples if not t][:5]
        (design_dir, results, list_dir), _, _ = self.one_pair(lists={"hand_species__drug": listed})
        code, text, err, out = self.run_it(design_dir, results, list_dir)
        self.assertEqual(code, 0, err)
        rows = {(r["reading"], r["learner"]): r for r in table(out / "contrasts.tsv")}
        # First reading: older recall 15 of 20, fresh 20 of 20. Second: 15 of 15 and 15 of 15.
        self.assertEqual((rows["first", "rf"]["d_recall"], rows["first", "rf"]["test"]), ("-0.250000", "40"))
        self.assertEqual((rows["second", "rf"]["d_recall"], rows["second", "rf"]["test"]), ("0.000000", "30"))
        self.assertIn("lists: 1 list files read; they name 1 datasets of the design\n", text)

    def test_a_short_list_leaves_the_pair_out_with_no_figure(self):
        # two_species__drug has no list file: it lists no genome, so its second reading is its first.
        a, b = samples_of("a", 12, 12), samples_of("b", 12, 12)
        datasets = [{"id": f"{n}_species__drug", "species": f"{n}_species", "tests": {2: s},
                     "scores": {(2, 1): scores_for(s, {s[0][0]}), (2, 2): scores_for(s)}, "pairs": [(1, 2)]}
                    for n, s in (("one", a), ("two", b))]
        lists = {"one_species__drug": [a[0][0], a[1][0], a[12][0]]}
        design_dir, results, list_dir = write_case(self.tmp, datasets, lists=lists)
        code, _, err, out = self.run_it(design_dir, results, list_dir)
        self.assertEqual(code, 0, err)
        lacking = table(out / "lacking.tsv")
        self.assertEqual(lacking, [{"clock": "public-day", "contrast": "stale-all", "horizon": "1",
                                    "dataset": "one_species__drug", "species": "one_species", "reason": an.SHORT}])
        for name in ("arms.tsv", "contrasts.tsv"):
            rows = table(out / name)
            second = [r for r in rows if r["reading"] == "second"]
            self.assertTrue(second)
            self.assertEqual({r["dataset"] for r in second}, {"two_species__drug"}, name)
            first = [{**r, "reading": "second"} for r in rows
                     if r["reading"] == "first" and r["dataset"] == "two_species__drug"]
            self.assertEqual(second, first, name)
        second = [r for r in table(out / "readings.tsv") if r["reading"] == "second"]
        self.assertEqual({r["datasets"] for r in second}, {"1"})
        self.assertEqual({r["datasets"] for r in table(out / "readings.tsv") if r["reading"] == "first"}, {"2"})


class Lists(Case):
    a, b = samples_of("a", 12, 12), samples_of("b", 12, 12)

    def two_datasets(self, lists):
        datasets = [{"id": f"{n}_species__drug", "species": f"{n}_species", "tests": {2: s},
                     "scores": {(2, 1): scores_for(s), (2, 2): scores_for(s)}, "pairs": [(1, 2)]}
                    for n, s in (("one", self.a), ("two", self.b))]
        return write_case(self.tmp, datasets, lists=lists)

    def stopped(self, design_dir, results, list_dir, message, out="out"):
        code, text, err, out = self.run_it(design_dir, results, list_dir, out=out)
        self.assertEqual((code, text), (1, ""))
        self.assertEqual(err, f"stagger_analyze: stopped: {message}\n")
        self.assertFalse(out.exists())
        return err

    def test_a_list_file_that_names_no_dataset_of_the_design_stops(self):
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": [self.a[0][0]],
                                                           "three_species__drug": [self.a[1][0]],
                                                           "two_species__drug": [self.b[0][0]]})
        self.stopped(design_dir, results, list_dir,
                     "1 of 3 list files name no dataset of the design: three_species__drug.txt")

    def test_a_list_file_named_in_upper_case_stops(self):
        # Names are matched letter for letter: the upper-case name names no dataset of the design.
        design_dir, results, list_dir = self.two_datasets({"ONE_SPECIES__DRUG": [self.a[0][0]],
                                                           "two_species__drug": [self.b[0][0]]})
        self.stopped(design_dir, results, list_dir,
                     "1 of 2 list files name no dataset of the design: ONE_SPECIES__DRUG.txt")

    def test_a_list_file_with_no_id_stops(self):
        # No line, then blank lines alone: with nothing listed, the second reading would be the first again (R41).
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": [self.a[0][0]],
                                                           "two_species__drug": []})
        for i, text in enumerate(("", "\n  \n\t\n")):
            with self.subTest(text=repr(text)):
                (Path(list_dir) / "two_species__drug.txt").write_text(text, encoding="utf-8")
                self.stopped(design_dir, results, list_dir,
                             "1 of 2 list files name no genome id: two_species__drug.txt", out=f"out{i}")

    def test_a_list_file_whose_ids_are_in_no_test_set_of_its_dataset_stops(self):
        # A test genome of the other dataset, and an id of no test set at all: no fit of two_species__drug scored
        # either. The stop names the file and no id.
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": [self.a[0][0]],
                                                           "two_species__drug": [self.a[1][0], "syn-nowhere-0001"]})
        err = self.stopped(design_dir, results, list_dir,
                           "1 of 2 list files name no test genome of their dataset: two_species__drug.txt")
        self.assertNotIn("syn-", err)

    def test_a_list_file_for_a_dataset_the_design_leaves_out_stops(self):
        # three_species__drug is in the design with no pair at any horizon. one_species__drug is left out on the
        # collection-year clock and kept on the public-day clock, so its list file passes.
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": [self.a[0][0]],
                                                           "three_species__drug": [self.a[1][0]]})
        path = Path(design_dir) / "design.json"
        design = json.loads(path.read_text())
        for dataset, clock in (("three_species__drug", "public-day"), ("one_species__drug", "collection-year")):
            design["datasets"].append({"id": dataset, "species_slug": dataset.split("__")[0], "antibiotic": "drug",
                                       "clock": clock, "built": True, "why_not": "no pair at any horizon",
                                       "windows": [], "pairs": {c: [] for c in sd.CONTRASTS}})
        path.write_text(json.dumps(design))
        self.stopped(design_dir, results, list_dir,
                     "1 of 2 list files are for a dataset the design leaves out: three_species__drug.txt")

    def test_a_byte_order_mark_is_dropped(self):
        # 20 R and 20 S, 5 of each listed, the first listed id after the mark. With the mark dropped, that id
        # matches its genome: 5 R and 5 S leave the second reading, which scores 15 and 15.
        samples = samples_of("p", 20, 20)
        listed = [g for g, t in samples if t][:5] + [g for g, t in samples if not t][:5]
        (design_dir, results, list_dir), _, _ = self.one_pair(lists={"hand_species__drug": listed})
        path = Path(list_dir) / "hand_species__drug.txt"
        self.assertEqual(path.read_text(encoding="utf-8").split("\n")[0], samples[0][0])
        path.write_text("﻿" + path.read_text(encoding="utf-8"), encoding="utf-8")
        code, _, err, out = self.run_it(design_dir, results, list_dir)
        self.assertEqual(code, 0, err)
        self.assertEqual(table(out / "lacking.tsv"), [])
        rows = {(r["reading"], r["test"], r["test_resistant"]) for r in table(out / "contrasts.tsv")}
        self.assertEqual(rows, {("first", "40", "20"), ("second", "30", "15")})

    def test_a_directory_named_as_a_list_file_stops(self):
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": [self.a[0][0]]})
        (Path(list_dir) / "two_species__drug.txt").mkdir()
        self.stopped(design_dir, results, list_dir, "two_species__drug.txt in the list dir is not a file")

    def test_a_list_dir_with_no_list_file_stops(self):
        # No list file, and then a file of another kind alone: a second reading with nothing listed is the first
        # reading again, so both would read as made.
        design_dir, results, list_dir = self.two_datasets({})
        for name in (None, "one_species__drug.tsv"):
            with self.subTest(file=name):
                if name:
                    (Path(list_dir) / name).write_text(self.a[0][0] + "\n", encoding="utf-8")
                self.assertEqual([p.name for p in Path(list_dir).glob("*.txt")], [])
                code, text, err, out = self.run_it(design_dir, results, list_dir)
                self.assertEqual((code, text), (1, ""))
                self.assertEqual(err, "stagger_analyze: stopped: the list dir holds no list file, so no second "
                                      "reading can be made from it\n")
                self.assertFalse(out.exists())

    def test_it_prints_the_files_and_the_datasets_they_name(self):
        # Two list files naming 3 genomes and 1: 2 files, 2 datasets, and no count of genomes.
        listed = [self.a[0][0], self.a[1][0], self.a[12][0]]
        design_dir, results, list_dir = self.two_datasets({"one_species__drug": listed,
                                                           "two_species__drug": [self.b[0][0]]})
        code, text, err, _ = self.run_it(design_dir, results, list_dir)
        self.assertEqual(code, 0, err)
        self.assertEqual([x for x in text.splitlines() if x.startswith("lists")],
                         ["lists: 2 list files read; they name 2 datasets of the design"])
        self.assertNotIn("genome", text)


class SameSize(Case):
    def test_stale_same_size_scores_against_recent_j(self):
        # Window 2, 20 R and 20 S. all(1), which is recent(1), misses 5 R; all(2) misses none; recent(2) misses
        # the same 5 R. Stale all reads 15 of 20 against 20 of 20, stale same size 15 of 20 against 15 of 20.
        samples = samples_of("z", 20, 20)
        missed = {g for g, t in samples if t and g.endswith(("r00", "r01", "r02", "r03", "r04"))}
        write_case(self.tmp, [{"id": "size_species__drug", "species": "size_species", "tests": {2: samples},
                               "scores": {(2, 1): scores_for(samples, missed), (2, 2): scores_for(samples)},
                               "recent": {(2, 2): scores_for(samples, missed)},
                               "pairs": [(1, 2)], "same_size_pairs": [(1, 2)]}])
        code, _, err, out = self.run_it(self.tmp / "design", self.tmp / "results")
        self.assertEqual(code, 0, err)
        rows = {(r["contrast"], r["learner"]): r for r in table(out / "contrasts.tsv")}
        self.assertEqual({rows["stale-same-size", l]["d_recall"] for l in fx.LEARNERS}, {"0.000000"})
        self.assertEqual({rows["stale-all", l]["d_recall"] for l in fx.LEARNERS}, {"-0.250000"})
        fresh = {r["recall"] for r in table(out / "arms.tsv") if r["contrast"] == "stale-same-size"
                 and r["arm"] == "fresh"}
        self.assertEqual(fresh, {"0.750000"})


class Scan(unittest.TestCase):
    def test_a_genome_id_stops_the_writing(self):
        with self.assertRaises(an.Stop):
            an.scan({"arms.tsv": "x\tsyn-z-0001\n"}, {"syn-z-0001"})
        an.scan({"arms.tsv": "x\tsyn-z-00010\n"}, {"syn-z-0001"})

    def test_the_shape_of_an_id_an_accession_or_a_sequence(self):
        # Each made up of zeros, or of one letter, so that no real id, accession or sequence can equal it.
        zeros = "0" * 12
        for text in ("000." + "0" * 4, "PRJNA" + zeros, "SAMN" + zeros, "GCA_" + zeros, "SRR" + zeros, "A" * 24):
            with self.assertRaises(an.Stop, msg=text):
                an.scan({"rule.tsv": f"a\t{text}\n"}, set())
        an.scan({"rule.tsv": "a\t-0.250000\t0.750000\n"}, set())


def reading_world(verdicts, reads):
    """A reading's verdicts, learner rules and readings, keyed as prediction_met reads them."""
    rules = {k: an.learner_rule(*v) for k, v in reads.items()}
    return verdicts, rules, reads


# The horizons each clock has, as build() gives them: P3 is read at every one of the collection-year clock.
TWO = {"public-day": [1], "collection-year": [1, 2]}
THREE = {"public-day": [1], "collection-year": [1, 2, 3]}


class Predictions(unittest.TestCase):
    def reads(self, clock, contrast, values):
        return {(clock, contrast, 1, l): v for l, v in zip(fx.LEARNERS, values)}

    def p3(self, changes=(), drop=()):
        """Collection-year verdicts at horizons 1 to 3 on both contrasts, each not shown, but for changes
        ((contrast, horizon, verdict), ...) and with the (contrast, horizon) of drop given no dataset."""
        v = {("collection-year", c, k): "not shown" for c in sd.CONTRASTS for k in THREE["collection-year"]}
        for c, k, x in changes:
            v[("collection-year", c, k)] = x
        for c, k in drop:
            del v[("collection-year", c, k)]
        return an.prediction_met("P3", v, {}, {}, THREE)

    def test_p1(self):
        ns, w = ("not shown", "not shown"), ("worse", "not shown")
        v = {("public-day", "stale-same-size", 1): "not shown"}
        self.assertTrue(an.prediction_met("P1", *reading_world(v, self.reads("public-day", "stale-same-size",
                                                                             [ns, ns, ns, ns, w])), TWO))
        self.assertFalse(an.prediction_met("P1", *reading_world(v, self.reads("public-day", "stale-same-size",
                                                                              [ns, ns, ns, w, w])), TWO))
        d = {("public-day", "stale-same-size", 1): "descriptive only"}
        self.assertFalse(an.prediction_met("P1", *reading_world(d, self.reads("public-day", "stale-same-size",
                                                                              [ns] * 5)), TWO))

    def test_p2(self):
        ll, rc, ns = ("not shown", "worse"), ("worse", "not shown"), ("not shown", "not shown")
        v = {("public-day", "stale-all", 1): "cost shown"}
        self.assertTrue(an.prediction_met("P2", *reading_world(v, self.reads("public-day", "stale-all",
                                                                             [ll, ll, ll, ns, ns])), TWO))
        # Three costs, one through recall alone: two through log loss.
        self.assertFalse(an.prediction_met("P2", *reading_world(v, self.reads("public-day", "stale-all",
                                                                              [ll, ll, rc, ns, ns])), TWO))

    def test_p3_not_shown_on_both_contrasts_at_every_horizon(self):
        self.assertTrue(self.p3())
        self.assertEqual(an.read_at("P3", THREE), [(c, k) for c in sd.CONTRASTS for k in (1, 2, 3)])

    def test_p3_needs_a_verdict_at_horizon_1_on_both_contrasts(self):
        for c in sd.CONTRASTS:
            self.assertFalse(self.p3(changes=[(c, 1, "descriptive only")]), c)
            self.assertFalse(self.p3(drop=[(c, 1)]), c)

    def test_p3_a_cost_at_a_later_horizon_is_not_met(self):
        for c in sd.CONTRASTS:
            for k in (2, 3):
                self.assertFalse(self.p3(changes=[(c, k, "cost shown")]), (c, k))

    def test_p3_a_mixed_verdict_at_a_later_horizon_is_not_met(self):
        for c in sd.CONTRASTS:
            for k in (1, 2, 3):
                self.assertFalse(self.p3(changes=[(c, k, "mixed")]), (c, k))

    def test_p3_a_descriptive_later_horizon_gives_no_verdict(self):
        # A later horizon with no verdict, descriptive only or with no dataset, leaves P3 met on the rest.
        self.assertTrue(self.p3(changes=[("stale-same-size", 2, "descriptive only")]))
        self.assertTrue(self.p3(changes=[("stale-all", 3, "descriptive only")], drop=[("stale-same-size", 3)]))
        self.assertFalse(self.p3(changes=[("stale-same-size", 2, "descriptive only"), ("stale-all", 3, "mixed")]))

    def test_the_table_names_each_horizon_with_no_verdict(self):
        v = {("first", "collection-year", c, 1): "not shown" for c in sd.CONTRASTS}
        v[("first", "collection-year", "stale-all", 2)] = "descriptive only"
        v[("first", "collection-year", "stale-same-size", 2)] = "not shown"
        v[("first", "collection-year", "stale-all", 3)] = "not shown"
        rows = {r["prediction"]: r for r in an.predictions(v, {}, {}, ("first",), THREE)}
        self.assertEqual(rows["P3"]["no_verdict"], "first stale-all horizon 2; first stale-same-size horizon 3")
        self.assertEqual(rows["P3"]["first"], "met")
        self.assertEqual(rows["P1"]["no_verdict"], "first stale-same-size horizon 1")
        self.assertEqual(rows["P2"]["no_verdict"], "first stale-all horizon 1")
        v = {(r, "collection-year", c, k): "not shown" for r in an.READINGS for c in sd.CONTRASTS for k in (1, 2)}
        rows = {r["prediction"]: r for r in an.predictions(v, {}, {}, an.READINGS, TWO)}
        self.assertEqual(rows["P3"]["no_verdict"], "none")
        self.assertEqual(an.PREDICTION_COLUMNS[-1], "no_verdict")

    def test_no_verdict_names_the_second_readings_horizons_too(self):
        # Both readings made. The second has no verdict at horizon 2: no dataset on stale all, descriptive only on
        # stale same size. P1 and P2 have no dataset in either reading.
        v = {(r, "collection-year", c, k): "not shown" for r in an.READINGS for c in sd.CONTRASTS for k in (1, 2)}
        del v[("second", "collection-year", "stale-all", 2)]
        v[("second", "collection-year", "stale-same-size", 2)] = "descriptive only"
        rows = {r["prediction"]: r for r in an.predictions(v, {}, {}, an.READINGS, TWO)}
        self.assertEqual((rows["P3"]["first"], rows["P3"]["second"], rows["P3"]["outcome"]), ("met", "met", "met"))
        self.assertEqual(rows["P3"]["no_verdict"], "second stale-all horizon 2; second stale-same-size horizon 2")
        self.assertEqual(rows["P1"]["no_verdict"], "first stale-same-size horizon 1; second stale-same-size horizon 1")
        self.assertEqual(rows["P2"]["no_verdict"], "first stale-all horizon 1; second stale-all horizon 1")

    def outcomes(self, first, second=None):
        """P3's row with the first reading met or not, and the second met, not met or not made (None)."""
        v = {}
        for reading, met in (("first", first), ("second", second)):
            if met is not None:
                for c in sd.CONTRASTS:
                    for k in (1, 2):
                        v[(reading, "collection-year", c, k)] = "not shown" if met else "cost shown"
        made = an.READINGS if second is not None else an.READINGS[:1]
        row = {r["prediction"]: r for r in an.predictions(v, {}, {}, made, TWO)}["P3"]
        return row["first"], row["second"], row["outcome"]

    def test_met_only_when_both_readings_are_made_and_both_meet_it(self):
        self.assertEqual(self.outcomes(True, True), ("met", "met", "met"))
        self.assertNotEqual(self.outcomes(True)[2], "met")
        self.assertEqual(an.outcome({"first": "met", "second": "met"}), "met")
        self.assertNotEqual(an.outcome({"first": "met"}), "met")

    def test_not_met_when_a_reading_made_fails_it(self):
        self.assertEqual(self.outcomes(True, False), ("met", "not met", "not met"))
        self.assertEqual(self.outcomes(False, True), ("not met", "met", "not met"))
        self.assertEqual(self.outcomes(False, False), ("not met", "not met", "not met"))
        self.assertEqual(self.outcomes(False), ("not met", "not read", "not met"))

    def test_otherwise_not_read(self):
        self.assertEqual(self.outcomes(True), ("met", "not read", "not read"))
        self.assertEqual(an.outcome({"first": "met"}), "not read")


def horizon_datasets(windows, cost_at=None):
    """12 datasets in 3 species on the collection-year clock, with test sets of 20 R and 20 S on windows 2 to
    windows and pairs at every horizon on both contrasts. A fit scores each resistant genome 0.9 and each
    susceptible one 0.1, but an older fit at horizon cost_at scores each resistant genome 0.6: log loss worse and
    recall the same for all 5 learners, a cost shown there."""
    datasets = []
    for i in range(12):
        tests = {j: samples_of(f"h{i:02d}w{j}", 20, 20) for j in range(2, windows + 1)}
        scores = {(j, k): [0.6 if t else 0.1 for _, t in s] if j - k == cost_at else scores_for(s)
                  for j, s in tests.items() for k in range(1, j + 1)}
        pairs = [(k, j) for j in tests for k in range(1, j)]
        datasets.append({"id": f"sp{i % 3}_species__drug_{i:02d}", "species": f"sp{i % 3}_species", "tests": tests,
                         "scores": scores, "recent": {(j, k): x for (j, k), x in scores.items() if k >= 2},
                         "pairs": pairs, "same_size_pairs": pairs})
    return datasets


class Horizons(Case):
    """Prediction 3 through build(): read at every horizon the collection-year clock's share set gives (R38)."""

    def p3(self, shares, windows, cost_at=None, name="case"):
        write_case(self.tmp / name, horizon_datasets(windows, cost_at), clock="collection-year", shares=shares)
        code, text, err = fx.run(an, self.tmp / name / "design", self.tmp / name / "results", self.tmp / name / "out")
        self.assertEqual(code, 0, err)
        p3 = {r["prediction"]: r for r in table(self.tmp / name / "out" / "predictions.tsv")}["P3"]
        verdicts = {(r["contrast"], r["horizon"]): r["verdict"] for r in table(self.tmp / name / "out" / "verdicts.tsv")
                    if r["clock"] == "collection-year"}
        return p3["first"], p3["no_verdict"], verdicts

    def test_p3_a_cost_at_the_last_horizon_through_build_is_not_met(self):
        # Both clocks with shares 40, 60 and 80: horizons 1 and 2. Not shown at horizon 1, a cost at horizon 2.
        first, no_verdict, verdicts = self.p3({c: (40, 60, 80) for c in sd.CLOCKS}, 3, cost_at=2)
        self.assertEqual(verdicts, {("stale-all", "1"): "not shown", ("stale-all", "2"): "cost shown",
                                    ("stale-same-size", "1"): "not shown", ("stale-same-size", "2"): "cost shown"})
        self.assertEqual((first, no_verdict), ("not met", "none"))

    def test_p3_reads_the_collection_year_clocks_own_horizons(self):
        # The public-day clock has 1 horizon and the collection-year clock 2, with a cost at horizon 2: not met.
        # Then the other way: the collection-year clock's 1 horizon, not shown, is all P3 reads, and it is met.
        first, no_verdict, _ = self.p3({"public-day": (50, 70), "collection-year": (40, 60, 80)}, 3, cost_at=2,
                                       name="fewer-public")
        self.assertEqual((first, no_verdict), ("not met", "none"))
        first, no_verdict, verdicts = self.p3({"public-day": (40, 60, 80), "collection-year": (50, 70)}, 2,
                                              name="fewer-collection")
        self.assertEqual(verdicts, {("stale-all", "1"): "not shown", ("stale-same-size", "1"): "not shown"})
        self.assertEqual((first, no_verdict), ("met", "none"))


class Stands(unittest.TestCase):
    design = {"clocks": {c: {"shares": [50, 70]} for c in sd.CLOCKS}}

    def test_a_verdict_stands_only_when_both_readings_give_it(self):
        v = {("first", "public-day", "stale-all", 1): "cost shown", ("second", "public-day", "stale-all", 1): "cost shown",
             ("first", "public-day", "stale-same-size", 1): "cost shown",
             ("second", "public-day", "stale-same-size", 1): "not shown",
             ("first", "collection-year", "stale-all", 1): "not shown"}
        rows = {(r["clock"], r["contrast"]): r["stands"] for r in an.stands(self.design, v, an.READINGS)}
        self.assertEqual(rows, {("public-day", "stale-all"): "cost shown",
                                ("public-day", "stale-same-size"): "no verdict stands",
                                ("collection-year", "stale-all"): "no verdict stands"})

    def test_without_lists_the_second_reading_is_not_made(self):
        v = {("first", "public-day", "stale-all", 1): "not shown"}
        rows = an.stands(self.design, v, an.READINGS[:1])
        self.assertEqual([(r["second"], r["stands"]) for r in rows], [("not read", "second reading not made")])


class Agreement(unittest.TestCase):
    """12 datasets in 3 species of 20 R and 20 S whose older fit scores each resistant genome 0.6 and the fresh fit
    0.9: log loss worse for every learner, recall not shown, so stale all at horizon 1 shows a cost and P2's first
    reading meets it. Run without lists, and with a list file a dataset that names 5 R and 5 S of its test
    genomes."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-agree-"))
        datasets, lists = [], {}
        for i in range(12):
            s = samples_of(f"g{i:02d}", 20, 20)
            older = [0.6 if t else 0.1 for _, t in s]
            datasets.append({"id": f"sp{i % 3}_species__drug_{i:02d}", "species": f"sp{i % 3}_species",
                             "tests": {2: s}, "scores": {(2, 1): older, (2, 2): scores_for(s)}, "pairs": [(1, 2)]})
            lists[datasets[-1]["id"]] = [g for g, t in s if t][:5] + [g for g, t in s if not t][:5]
        design_dir, results, list_dir = write_case(cls.tmp, datasets, lists=lists)
        cls.got = {}
        for name, extra in (("unlisted", []), ("listed", ["--listed", list_dir])):
            code, text, err = fx.run(an, design_dir, results, cls.tmp / name, *extra)
            assert code == 0, err
            cls.got[name] = ({r["prediction"]: r for r in table(cls.tmp / name / "predictions.tsv")},
                             {(r["clock"], r["contrast"], r["horizon"]): r for r in
                              table(cls.tmp / name / "stands.tsv")}, text,
                             {(r["reading"], r["test"], r["test_resistant"])
                              for r in table(cls.tmp / name / "contrasts.tsv")})

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_without_lists_the_prediction_is_not_read_and_no_verdict_stands(self):
        predictions, stands, text, _ = self.got["unlisted"]
        p2, row = predictions["P2"], stands[("public-day", "stale-all", "1")]
        self.assertEqual((row["first"], row["second"], row["stands"]),
                         ("cost shown", "not read", "second reading not made"))
        self.assertEqual((p2["first"], p2["second"], p2["outcome"]), ("met", "not read", "not read"))
        self.assertIn("P2 (" + p2["rule"] + "): not read\n", text)
        # P1 reads stale same size, which has no dataset: no verdict, so the first reading does not meet it.
        p1 = predictions["P1"]
        self.assertEqual((p1["first"], p1["outcome"], p1["no_verdict"]),
                         ("not met", "not met", "first stale-same-size horizon 1"))

    def test_with_both_readings_the_prediction_and_the_verdict_agree(self):
        predictions, stands, text, tested = self.got["listed"]
        p2, row = predictions["P2"], stands[("public-day", "stale-all", "1")]
        self.assertEqual((row["first"], row["second"], row["stands"]), ("cost shown", "cost shown", "cost shown"))
        self.assertEqual((p2["first"], p2["second"], p2["outcome"]), ("met", "met", "met"))
        # The second reading is not the first again: each dataset's 10 listed test genomes are out of it.
        self.assertEqual(tested, {("first", "40", "20"), ("second", "30", "15")})
        self.assertIn("lists: 12 list files read; they name 12 datasets of the design\n", text)


class MadeUpWorld(unittest.TestCase):
    """folds, made-up results and the analysis, with lists, on the world of 3 species and 5 drugs."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="stagger-world-"))
        root, fetch, built, labels = fx.write_world(cls.tmp / "w")
        design_dir = cls.tmp / "design"
        code, _, err = fx.run(sd, "folds", root, fetch, built, design_dir)
        assert code == 0, err
        design = fx.write_results(design_dir, cls.tmp / "results", labels)
        lists = cls.tmp / "lists"
        lists.mkdir()
        cls.short = "alpha_species__drug_one"
        tests = set()
        for d in design["datasets"]:
            if d["id"] == cls.short and not d["why_not"]:
                for w in d["windows"]:
                    body = json.loads((design_dir / "fold-files" / w["fits"][0]["file"]).read_text())
                    tests.update(body["test"][0][:2])
        (lists / f"{cls.short}.txt").write_text("".join(g + "\n" for g in sorted(tests)))
        cls.out = cls.tmp / "out"
        cls.code, cls.text, cls.err = fx.run(an, design_dir, cls.tmp / "results", cls.out, "--listed", lists)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp)

    def test_it_ran_and_wrote_each_table(self):
        self.assertEqual(self.code, 0, self.err)
        names = sorted(p.name for p in self.out.iterdir())
        self.assertEqual(names, ["arms.tsv", "contrasts.tsv", "lacking.tsv", "predictions.tsv", "readings.tsv",
                                 "rule.tsv", "stands.tsv", "verdicts.tsv"])

    def test_readings_count_the_datasets_with_a_difference_of_zero(self):
        head = (self.out / "readings.tsv").read_text().splitlines()[0].split("\t")
        self.assertEqual(head[-3:], ["datasets_worse", "datasets_better", "datasets_zero"])
        self.assertNotIn("equal", " ".join(head))
        for r in table(self.out / "readings.tsv"):
            self.assertEqual(int(r["datasets_worse"]) + int(r["datasets_better"]) + int(r["datasets_zero"]),
                             int(r["datasets"]))

    def test_no_genome_id_and_no_study_in_what_it_writes(self):
        for p in self.out.iterdir():
            self.assertNotIn("syn-", p.read_text(), p.name)

    def test_the_short_list_dataset_is_lacking_in_the_second_reading(self):
        lacking = table(self.out / "lacking.tsv")
        self.assertTrue(lacking)
        self.assertEqual({(r["dataset"], r["reason"]) for r in lacking}, {(self.short, an.SHORT)})
        second = {r["dataset"] for r in table(self.out / "contrasts.tsv") if r["reading"] == "second"}
        self.assertNotIn(self.short, second)

    def test_three_predictions_each_met_or_not_met(self):
        rows = table(self.out / "predictions.tsv")
        self.assertEqual([r["prediction"] for r in rows], ["P1", "P2", "P3"])
        for r in rows:
            self.assertIn(r["outcome"], ("met", "not met"))

    def test_a_row_per_clock_reading_contrast_and_horizon_in_the_verdicts(self):
        rows = table(self.out / "verdicts.tsv")
        keys = [(r["clock"], r["reading"], r["contrast"], r["horizon"]) for r in rows]
        self.assertEqual(len(keys), len(set(keys)))
        self.assertEqual({r["reading"] for r in rows}, {"first", "second"})


if __name__ == "__main__":
    unittest.main()
