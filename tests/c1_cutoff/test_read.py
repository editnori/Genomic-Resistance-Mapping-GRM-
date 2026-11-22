"""scripts/c1-cutoff/read.py: the gates and the page's tables of the small-cohort run (docs/cohort-draw.md).

Every species, dataset, result and score here is SYNTHETIC: made-up species
directories in the form run_species.sh and the second scoring write, scored
by this tree's score.py from SYNTHETIC results (tests/c1_cutoff/common.py).
Three species, synthia_alpha (two datasets, one with a study plan),
synthia_beta and synthia_gamma (one each); lines 10, 25 and 50, three draws,
and the twin, so 12 results a dataset, 4 datasets and 48 results. The
estimate is the same at every draw and for every learner, so the figures
are worked by hand below.
"""

import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
import shutil
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from tests.c1_cutoff import common

score = common.load("score")
read = common.load("read")

SCORER = str(common.KIT / "score.py")
SCORER_SHA = hashlib.sha256(Path(SCORER).read_bytes()).hexdigest()
CODE = "c0de" * 16
BINARY = "b1" * 32
LINES, DRAWS = (10, 25, 50), 3
SPECIES = {"synthia_alpha": ["synthia_alpha__drug_a", "synthia_alpha__drug_b"],
           "synthia_beta": ["synthia_beta__drug_c"], "synthia_gamma": ["synthia_gamma__drug_d"]}
STUDY = {"synthia_alpha__drug_b"}
NAMES = {"synthia_alpha__drug_a": "drug a/drug e", "synthia_alpha__drug_b": "drug b",
         "synthia_beta__drug_c": "drug c", "synthia_gamma__drug_d": "drug d"}
# The estimate's calls at each r: (tp, tn) of r resistant and r susceptible.
EST = {10: (10, 10), 25: (23, 22), 50: (45, 40)}
TMP = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def seal_second(second):
    with open(second / "SHA256SUMS", "w") as f:
        for name in read.SECOND_FILES:
            f.write(f"{sha(second / name)}  {name}\n")


def seal_first(first):
    names = sorted(str(p.relative_to(first)) for p in first.rglob("*") if p.is_file()
                   and p.name not in ("SHA256SUMS", "DONE", "FAILED") and not str(p.relative_to(first)).startswith("logs"))
    with open(first / "SHA256SUMS", "w") as f:
        for name in names:
            f.write(f"{sha(first / name)}  {name}\n")
    sp = first.name.split(".")[0]
    (first / "DONE").write_text(f"DONE {sp} 2026-11-22T00:00:00Z binary {BINARY} code {CODE} "
                                f"SHA256SUMS {sha(first / 'SHA256SUMS')}\n")


def build(root):
    """The SYNTHETIC copies: root/back/<species> and root/back2/<species>, and root/names.tsv."""
    root.mkdir(parents=True)
    with open(root / "names.tsv", "w") as f:
        f.write("dataset_id\tspecies_slug\tantibiotic\n")
        for d, a in NAMES.items():
            f.write(f"{d}\t{d.split('__')[0]}\t{a}\n")
    for sp, datasets in SPECIES.items():
        work, first, second = root / "work" / sp, root / "back" / sp, root / "back2" / sp
        for p in (work / "cohorts", first / "results" / "cohorts", first / "score", first / "logs",
                  second / "logs"):
            p.mkdir(parents=True)
        common.write_release(work / "release", [common.report(common.samples(100, 100), dataset=d) for d in datasets])
        rows = []
        for d in datasets:
            for r in LINES:
                for j in range(DRAWS):
                    for seed in ((42, 45) if r == 25 else (42,)):
                        cohort = json.dumps({"SYNTHETIC": True, "dataset": d, "r": r, "draw": j}).encode()
                        e = common.expect_row(r=r, draw=j, seed=seed, dataset=d, study="1" if d in STUDY else "0",
                                              cohort_sha=hashlib.sha256(cohort).hexdigest())
                        (work / "cohorts" / e["cohort"]).write_bytes(cohort)
                        tp, tn = EST[r]
                        doc = common.result(r=r, draw=j, seed=seed, dataset=d, samples_=200, resistant=100,
                                            study=d in STUDY, cohort_sha=e["cohort_sha256"], est_tp=tp, est_tn=tn,
                                            ba=(tp / r + tn / r) / 2)
                        doc["antibiotic"] = NAMES[d]
                        (first / "results" / "cohorts" / e["file"]).write_text(json.dumps(doc))
                        rows.append(e)
        with open(first / "expect.tsv", "w") as f:
            f.write("\t".join(read.EXPECT_COLUMNS) + "\n")
            for e in rows:
                f.write("\t".join(e[c] for c in read.EXPECT_COLUMNS) + "\n")
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = score.main([str(first / "expect.tsv"), str(work / "cohorts"), str(first / "results" / "cohorts"),
                             str(work / "release"), str(second / "score")])
        assert rc == 0, out.getvalue()
        (second / "logs" / "score.log").write_text(out.getvalue())
        seal_second(second)
        (first / "binary.txt").write_text(f"commit {'9' * 40}\ngrm {BINARY}\ncode {CODE}\n")
        (first / "lists.txt").write_text(f"species\t{sp}\ndatasets\t{','.join(datasets)}\nlines\t10,25,50\n"
                                         f"draws\t{DRAWS}\ntwin\t1\n")
        (first / "results" / "draws.tsv").write_text("SYNTHETIC\n")
        (first / "score" / "count.tsv").write_text("SYNTHETIC first scoring, never opened\n")
        (first / "logs" / "run.log").write_text("SYNTHETIC\n")
        seal_first(first)


def rescore(root, sp):
    """Score species sp of root again into its second copy and seal both copies."""
    first, second = root / "back" / sp, root / "back2" / sp
    shutil.rmtree(second / "score")
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = score.main([str(first / "expect.tsv"), str(root / "work" / sp / "cohorts"), str(first / "results" / "cohorts"),
                         str(root / "work" / sp / "release"), str(second / "score")])
    assert rc == 0, out.getvalue()
    (second / "logs" / "score.log").write_text(out.getvalue())
    seal_second(second)
    seal_first(first)


def setUpModule():
    """The template, and the two runs the table tests read: the default, and the other answer to each question.

    The other run's copy has no lineage split at draw 2 of r 10 of synthia_beta
    drug c, so that dataset has 2 shared draws at line 10 and is not scored there.
    """
    global TMP, DEFAULT, OTHER
    TMP = tempfile.TemporaryDirectory()
    try:
        base = Path(TMP.name)
        build(base / "template")
        shutil.copytree(base / "template", base / "default")
        DEFAULT = (base / "default",) + run(args(base / "default"))
        shutil.copytree(base / "template", base / "other")
        f = base / "other" / "back" / "synthia_beta" / "results" / "cohorts" / "synthia_beta__drug_c.r10.d2.s42.cohort.json"
        doc = json.loads(f.read_text())
        doc["estimates"][1] = {"scheme": "lineage", "n_groups": 2, "unavailable": "too few lineages", "learners": []}
        f.write_text(json.dumps(doc))
        rescore(base / "other", "synthia_beta")
        OTHER = (base / "other",) + run(args(base / "other", **{
            "--measures-over": "shared", "--line-50": "beside", "--twin-rule": "void", "--dataset-differences": "yes"}))
    except BaseException:
        TMP.cleanup()
        raise


def tearDownModule():
    TMP.cleanup()


def args(root, **over):
    a = {"--out": str(root / "out"), "--code-sha256": CODE, "--scorer": SCORER, "--scorer-sha256": SCORER_SHA,
         "--species": "3", "--datasets": "4", "--results": "48", "--per-dataset": "12",
         "--names": str(root / "names.tsv"), "--measures-over": "available", "--line-50": "release-row",
         "--twin-rule": "beside", "--dataset-differences": "no"}
    pairs = over.pop("pairs", [(root / "back" / sp, root / "back2" / sp) for sp in SPECIES])
    a.update(over)
    out = [x for k, v in a.items() for x in (k, v)]
    for f, s in pairs:
        out += ["--pair", str(f), str(s)]
    return out


def run(argv):
    err, out = io.StringIO(), io.StringIO()
    with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
        rc = read.main(argv)
    return rc, out.getvalue(), err.getvalue()


def tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


class Copy(unittest.TestCase):
    """Each test works on its own copy of the SYNTHETIC template."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.root = Path(self._tmp.name) / "c"
        shutil.copytree(Path(TMP.name) / "template", self.root)

    def tearDown(self):
        self._tmp.cleanup()

    def stops(self, argv, text):
        rc, out, err = run(argv)
        self.assertEqual(rc, 2, err)
        self.assertTrue(err.startswith("stop: "), err)
        self.assertIn(text, err)
        self.assertFalse((self.root / "out").exists() and any((self.root / "out").iterdir()))
        return err


class Reading(unittest.TestCase):
    """The default run: --measures-over available --line-50 release-row --twin-rule beside --dataset-differences no."""

    def out(self, name):
        return tsv(DEFAULT[0] / "out" / name)

    def test_a_clean_run_writes_every_table_and_exits_0(self):
        root, rc, out, err = DEFAULT
        self.assertEqual(rc, 0, err)
        for name in read.OUTPUTS:
            self.assertTrue((root / "out" / name).exists(), name)
        self.assertIn("7 counts: 3 species, 4 datasets, 48 results, none twice", out)
        self.assertIn("species: 3: Synthia alpha, Synthia beta, Synthia gamma", out)
        self.assertEqual((root / "out" / "summary.txt").read_text(), out)
        self.assertNotIn("9000.", "".join(p.read_text() for p in (root / "out").rglob("*") if p.is_file()))

    def test_the_deciding_rows_are_the_scorers_own_cells(self):
        count, dec = self.out("count/count.tsv"), self.out("page-deciding.tsv")
        want = [c for c in count if (c["split"], c["truth_set"], c["coverage"]) == ("lineage", "all", "levelled")
                and c["line"] in ("10", "50")]
        self.assertEqual(len(dec), len(want))
        self.assertEqual(len(dec), 12)  # two lines, five learners and the line row each
        for c in want:
            d = next(x for x in dec if (x["line"], x["learner"]) == (c["line"], c["learner"]))
            for k in read.DECIDING_COLUMNS:
                if k in c:
                    self.assertEqual(d[k], c[k], (c["line"], c["learner"], k))

    def test_the_lines_read_as_worked_by_hand(self):
        # Truth on all scored genomes, lineage split: r 10 scores 178 genomes,
        # 89 resistant, called 84 and 79 right; r 25 148, 74, 69 and 64; r 50
        # 98, 49, 44 and 39. Estimates: 1.0 at 10, 0.90 at 25, 0.85 at 50.
        truth = {r: (Fraction(n - 5, n) + Fraction(n - 10, n)) / 2 for r, n in ((10, 89), (25, 74), (50, 49))}
        est = {r: (Fraction(tp, r) + Fraction(tn, r)) / 2 for r, (tp, tn) in EST.items()}
        err = {r: abs(est[r] - truth[r]) for r in LINES}
        dec = self.out("page-deciding.tsv")
        for line, diff, cls in (("10", err[10] - err[25], "worse"), ("50", err[50] - err[25], "same")):
            for l in common.LEARNERS:
                d = next(x for x in dec if (x["line"], x["learner"]) == (line, l))
                self.assertEqual((d["datasets"], d["entered"], d["not_scored"]), ("4", "4", "0"))
                self.assertEqual(d[cls], "4")
                self.assertAlmostEqual(float(d["median_difference"]), float(diff), places=9)
                self.assertEqual((d["interval_lo"], d["interval_hi"]), (d["median_difference"],) * 2)
                self.assertEqual(d["species"], "3")
                self.assertEqual((d["covered"], d["coverage_n"], d["base_covered"], d["base_coverage_n"]),
                                 ("12", "12", "12", "12"))
                self.assertEqual((d["coverage_at_line"], d["coverage_minus_25"]), ("1.000000", "0.000000"))
                self.assertEqual(d["reading"], "does not meet: the interval's upper end" if line == "10" else "meets")
            row = next(x for x in dec if (x["line"], x["row"]) == (line, "line"))
            self.assertEqual(row["learners_meeting"], "0" if line == "10" else "5")
            # The twin meets on 5 of 5 here, so a line that is no candidate is "not a candidate" alone (ruling 3).
            self.assertEqual(row["reading"], read.NOT_CANDIDATE if line == "10"
                             else read.CANDIDATE + "; the twin: 5 of 5 learners meet")
        self.assertIn(f"line 50: 5 of 5 learners meet (scm,cart,rf,gbt,lasso): {read.CANDIDATE}",
                      (DEFAULT[0] / "out" / "summary.txt").read_text())

    def test_every_other_row_is_written_and_marked_as_deciding_nothing(self):
        count, dec, twin, other = (self.out(n) for n in ("count/count.tsv", "page-deciding.tsv", "page-twin.tsv",
                                                         "page-other.tsv"))
        self.assertEqual(len(dec) + len(twin) + len(other), len(count))
        self.assertTrue(all(r["decides"].startswith("no: ") and len(r["decides"]) > 4 for r in other))
        why = lambda **k: {r["decides"] for r in other if all(r[a] == b for a, b in k.items())}
        self.assertEqual(why(split="random", truth_set="all", coverage="levelled", line="10"), {"no: the random split"})
        self.assertEqual(why(split="study", truth_set="all", coverage="levelled", line="10"), {"no: the study split"})
        self.assertEqual(why(split="lineage", truth_set="new_lineage", coverage="levelled", line="10"),
                         {"no: truth on new lineages"})
        self.assertEqual(why(split="lineage", truth_set="all", coverage="fold_range", line="10"),
                         {"no: the fold range's coverage"})
        self.assertEqual(why(split="random", truth_set="all", coverage="levelled", line="25-s45"),
                         {"no: the random split, the twin"})
        # The study split counts only the dataset with a study plan.
        st = [r for r in other if r["split"] == "study" and r["row"] == "learner"]
        self.assertTrue(st and all(r["datasets"] == "1" for r in st))
        self.assertTrue(all(r["labels"] == read.LABELS for r in other + dec + twin))

    def test_the_twin_stands_apart_and_is_never_a_candidate(self):
        twin = self.out("page-twin.tsv")
        self.assertEqual(len(twin), 6)
        line = next(r for r in twin if r["row"] == "line")
        # Seed 45 gives the same calls here, so all five learners meet and it is still no candidate.
        # The scorer's cell counts 3 or more meeting; the table names it so, and no twin cell says candidate.
        self.assertEqual((line["line"], line["learners_meeting"], line["three_or_more_meet"]), ("25-s45", "5", "1"))
        self.assertNotIn("candidate", twin[0])
        self.assertNotIn("candidate", self.out("page-other.tsv")[0])
        self.assertEqual(line["reading"], "the count with no change; never a candidate; " + read.TWIN_HOLDS)
        self.assertTrue(all(r["same"] == "4" for r in twin if r["row"] == "learner"))

    def test_species_and_datasets_are_named_in_full(self):
        sp = self.out("page-species.tsv")
        self.assertEqual([(r["species"], r["datasets"], r["results"], r["study_datasets"], r["antibiotics"]) for r in sp],
                         [("Synthia alpha", "2", "24", "1", "drug a/drug e, drug b"), ("Synthia beta", "1", "12", "0", "drug c"),
                          ("Synthia gamma", "1", "12", "0", "drug d")])
        ds = self.out("page-datasets.tsv")
        self.assertEqual([(r["species"], r["antibiotic"], r["results"], r["study_plan"]) for r in ds],
                         [("Synthia alpha", "drug a/drug e", "12", "no"), ("Synthia alpha", "drug b", "12", "yes"),
                          ("Synthia beta", "drug c", "12", "no"), ("Synthia gamma", "drug d", "12", "no")])

    def test_the_entry_lists_each_dataset_with_no_difference_unless_asked(self):
        ent = self.out("page-entry.tsv")
        self.assertEqual(len(ent), 4 * 3 * 5)  # datasets, lines 10, 50 and the twin, learners
        self.assertTrue(all(r["entered"] == "yes" and r["shared_draws"] == "3" and r["difference"] == "" for r in ent))

    def test_rules_measures_carry_their_denominators(self):
        m = self.out("page-measures.tsv")
        get = lambda **k: next(r for r in m if all(r[a] == b for a, b in k.items()))
        # The truth's recall at 25, lineage split, all scored genomes: 69 of 74 in each of
        # 4 datasets x 3 draws.
        r = get(split="lineage", truth_set="all", learner="rf", line="25", who="truth", measure="recall")
        self.assertEqual((r["pooled_count"], r["pooled_denominator"], r["results"], r["datasets"]),
                         ("828", "888", "12", "4"))
        self.assertEqual(r["pooled_value"], "%.6f" % (69 / 74))
        self.assertEqual(r["median_over_datasets"], "%.6f" % (69 / 74))
        # The estimate's major error rate at 10: 0 of 10 a result.
        r = get(split="lineage", truth_set="all", learner="rf", line="10", who="estimate", measure="major")
        self.assertEqual((r["pooled_count"], r["pooled_denominator"], r["pooled_value"]), ("0", "120", "0.000000"))
        # Log loss: the truth's 12.5 over 178 genomes at r 10.
        r = get(split="lineage", truth_set="all", learner="rf", line="10", who="truth", measure="log_loss")
        self.assertEqual((r["pooled_count"], r["pooled_denominator"], r["pooled_value"]), ("", "2136", "%.6f" % (12.5 / 178)))
        # Whether the estimate's Wilson interval holds the truth's recall: at 10 the estimate is 10 of 10.
        r = get(split="lineage", truth_set="all", learner="rf", line="10", who="truth", measure="recall")
        lo, hi = score.wilson(10, 10)
        inside = lo <= 84 / 89 <= hi
        self.assertEqual((r["truth_inside_estimate_wilson"], r["of_results"]), ("12" if inside else "0", "12"))
        self.assertEqual({r["side"] for r in m}, {"line"})
        self.assertTrue(all(r["decides"] == "no: beside the error" and r["labels"] == read.LABELS for r in m))


class OtherAnswers(unittest.TestCase):
    """The run with the other answer to each question and a dataset with 2 shared draws at line 10."""

    def out(self, name):
        return tsv(OTHER[0] / "out" / name)

    def test_it_runs(self):
        self.assertEqual(OTHER[1], 0, OTHER[3])
        self.assertIn("rules chosen: --measures-over shared, --line-50 beside, --twin-rule void, "
                      "--dataset-differences yes", OTHER[2])

    def test_twin_rule_void_makes_the_block_decide_nothing_when_the_twin_meets(self):
        rows = [r for r in self.out("page-deciding.tsv") if r["row"] == "line"]
        self.assertTrue(rows)
        for r in rows:
            self.assertTrue(r["reading"].startswith("decides nothing: the twin meets"), r["reading"])
        self.assertIn("under --twin-rule void the block decides nothing", OTHER[2])

    def test_line_50_beside_moves_it_to_the_rows_that_decide_nothing(self):
        self.assertEqual({r["line"] for r in self.out("page-deciding.tsv")}, {"10"})
        fifty = [r for r in self.out("page-other.tsv") if r["line"] == "50"
                 and (r["split"], r["truth_set"], r["coverage"]) == ("lineage", "all", "levelled")]
        self.assertEqual(len(fifty), 6)
        self.assertEqual({r["decides"] for r in fifty}, {"no: line 50 is not a lower line"})

    def test_a_dataset_with_two_shared_draws_is_listed_as_not_scored(self):
        ent = self.out("page-entry.tsv")
        beta = [r for r in ent if r["dataset"] == "synthia_beta__drug_c" and r["line"] == "10"]
        self.assertEqual({(r["entered"], r["line_draws"], r["shared_draws"], r["difference"]) for r in beta},
                         {("no: 2 shared draws, under 3", "0,1", "2", "")})
        dec = self.out("page-deciding.tsv")
        self.assertEqual({(r["entered"], r["not_scored"]) for r in dec if r["line"] == "10" and r["row"] == "learner"},
                         {("3", "1")})

    def test_the_differences_listed_have_the_scorers_median(self):
        ent, dec = self.out("page-entry.tsv"), self.out("page-deciding.tsv")
        for l in common.LEARNERS:
            md = next(r for r in dec if (r["line"], r["learner"]) == ("10", l))["median_difference"]
            got = {r["difference"] for r in ent if r["line"] == "10" and r["learner"] == l and r["entered"] == "yes"}
            self.assertEqual(got, {md})

    def test_measures_over_shared_draws_pair_each_line_with_25(self):
        m = self.out("page-measures.tsv")
        self.assertEqual({r["side"] for r in m}, {"line", "25"})
        self.assertFalse([r for r in m if r["line"] == "25"])
        key = lambda r: (r["split"], r["truth_set"], r["learner"], r["line"], r["side"], r["who"], r["measure"])
        r = next(x for x in m if key(x) == ("lineage", "all", "rf", "10", "25", "truth", "recall"))
        # 3 datasets enter at line 10, each with 3 shared draws: 9 results at 25.
        self.assertEqual((r["pooled_count"], r["pooled_denominator"], r["results"]), ("621", "666", "9"))
        r = next(x for x in m if key(x) == ("lineage", "all", "rf", "50", "25", "truth", "recall"))
        self.assertEqual((r["pooled_count"], r["pooled_denominator"], r["results"]), ("828", "888", "12"))


def count_row(**cells):
    """One SYNTHETIC row of the scorer's count.tsv as read.py reads it: every column text, empty unless given."""
    row = {c: "" for c in score.COUNT_COLUMNS}
    row.update({"split": "lineage", "truth_set": "all", "coverage": "levelled"})
    row.update({k: str(v) for k, v in cells.items()})
    return row


def line_rows(line, meeting, entered=4, interval_hi="0.01", covered=(9, 10), base=(9, 10)):
    """SYNTHETIC deciding-view rows of one line: five learner rows, the first `meeting` of which meet, and its line row."""
    rows = []
    for i, l in enumerate(common.LEARNERS):
        meets = i < meeting
        rows.append(count_row(row="learner", learner=l, line=line, meets=int(meets), datasets=4, entered=entered,
                              not_scored=4 - entered, interval_hi=interval_hi if meets else "0.05",
                              covered=covered[0], coverage_n=covered[1], base_covered=base[0],
                              base_coverage_n=base[1]))
    rows.append(count_row(row="line", learner=score.LINE_ROW, line=line, learners_meeting=meeting,
                          meeting=",".join(common.LEARNERS[:meeting]), candidate=int(meeting >= 3),
                          decides=int(not line.endswith("-s45"))))
    return rows


def score_row(dataset, r, draw, learner="scm", seed=42, split="lineage", truth_set="all", available=True, **cells):
    """One SYNTHETIC scores.tsv row as read.py reads it: 10 resistant and 10 susceptible genomes each side, all right."""
    row = {"dataset": dataset, "species_slug": dataset.split("__")[0], "r": str(r), "draw": str(draw),
           "seed": str(seed), "split": split, "truth_set": truth_set, "learner": learner,
           "available": "1" if available else "0"}
    for who in ("estimate", "truth"):
        for name in read.RATES:
            right = name in ("recall", "specificity")
            row[f"{who}_{name}_k"], row[f"{who}_{name}_n"] = ("10" if right else "0"), "10"
        row[f"{who}_log_loss"], row[f"{who}_log_loss_n"] = "0.1", "20"
    for name in read.RATES:
        row[f"{name}_covers"] = "1"
    row.update({k: str(v) for k, v in cells.items()})
    return row


def lines_row(dataset, line, draws, learner="scm", split="lineage", truth_set="all", error="0.01"):
    """One SYNTHETIC lines.tsv row: its available draws, each with the same error."""
    return {"dataset": dataset, "species_slug": dataset.split("__")[0], "split": split, "truth_set": truth_set,
            "learner": learner, "line": line, "available_draws": ",".join(map(str, draws)),
            "errors": ",".join([error] * len(draws))}


class Measures(unittest.TestCase):
    """read.py's entry and Rules' measures on SYNTHETIC rows of the scorer's files."""

    def test_a_learner_with_no_result_available_at_a_line_keeps_its_rows(self):
        d = "synthia_alpha__drug_a"
        scores = [score_row(d, r, j) for r in (10, 25) for j in range(3)]
        scores += [score_row(d, 10, j, learner="gbt", available=False) for j in range(3)]
        scores += [score_row(d, 25, j, learner="gbt") for j in range(3)]
        lrows = [lines_row(d, "10", range(3)), lines_row(d, "25", range(3)), lines_row(d, "10", [], learner="gbt"),
                 lines_row(d, "25", range(3), learner="gbt")]
        ent = read.entry(lrows, ["10", "25"])
        for over, sides in (("shared", ("line", "25")), ("available", ("line",))):
            meas = read.measures(scores, ent, over, ["10", "25"])
            gbt = [m for m in meas if m["learner"] == "gbt" and m["line"] == "10"]
            self.assertEqual(sorted({m["side"] for m in gbt}), sorted(sides), over)
            self.assertEqual(len(gbt), len(sides) * 10, over)
            self.assertTrue(all(m["results"] == 0 and m["pooled_value"] == "" for m in gbt), over)
            scm = [m for m in meas if m["learner"] == "scm" and m["line"] == "10" and m["side"] == "line"]
            self.assertTrue(scm and all(m["results"] == 3 for m in scm), over)


    def test_the_truth_inside_the_estimates_wilson_interval_counts_only_results_it_lies_inside(self):
        d = "synthia_alpha__drug_a"
        scores = [score_row(d, 10, 0, recall_covers="1"), score_row(d, 10, 1, recall_covers="0"),
                  score_row(d, 10, 2, recall_covers="")] + [score_row(d, 25, j) for j in range(3)]
        ent = read.entry([lines_row(d, "10", range(3)), lines_row(d, "25", range(3))], ["10", "25"])
        m = next(m for m in read.measures(scores, ent, "shared", ["10", "25"])
                 if (m["line"], m["side"], m["who"], m["measure"]) == ("10", "line", "truth", "recall"))
        self.assertEqual((m["truth_inside_estimate_wilson"], m["of_results"], m["results"]), (1, 2, 3))

    def test_the_median_over_datasets_is_a_median_and_unavailable_results_stay_out(self):
        # Recall 1 of 10, 2 of 10 and 9 of 10 in three datasets: median 0.2, mean 0.4. A fourth result,
        # not available, with 10 of 10 on another 10 genomes, enters nothing.
        ds = ["synthia_alpha__drug_a", "synthia_beta__drug_c", "synthia_gamma__drug_d"]
        scores = [score_row(d, 10, 0, truth_recall_k=k) for d, k in zip(ds, (1, 2, 9))]
        scores.append(score_row(ds[0], 10, 1, available=False))
        for over in ("available", "shared"):
            lrows = [lines_row(d, ln, range(1) if ln == "10" else ()) for d in ds for ln in ("10", "25")]
            meas = read.measures(scores, read.entry(lrows, ["10", "25"]), over, ["10", "25"])
            if over == "available":
                m = next(m for m in meas if (m["line"], m["who"], m["measure"]) == ("10", "truth", "recall"))
                self.assertEqual((m["pooled_count"], m["pooled_denominator"], m["results"], m["datasets"]),
                                 (12, 30, 3, 3))
                self.assertEqual((m["pooled_value"], m["median_over_datasets"]), ("0.400000", "0.200000"))
            else:  # no dataset has 3 shared draws, so nothing enters the shared measures
                m = next(m for m in meas if (m["line"], m["side"], m["who"], m["measure"]) == ("10", "line", "truth", "recall"))
                self.assertEqual((m["results"], m["pooled_denominator"]), (0, 0))

    def test_a_dataset_with_no_line_at_25_is_listed_as_such_and_checked_against_the_count(self):
        d, e = "synthia_alpha__drug_a", "synthia_beta__drug_c"
        lrows = [lines_row(d, "10", range(3)), lines_row(d, "25", range(3)), lines_row(e, "10", range(3))]
        ent = read.entry(lrows, ["10", "25"])
        self.assertEqual(ent[("lineage", "all", "scm", "10", e)]["why"], "no line at 25")
        self.assertEqual(ent[("lineage", "all", "scm", "10", d)]["why"], "")
        row = count_row(row="learner", learner="scm", line="10", datasets=2, entered=1, not_scored=1,
                        no_base_line=1, median_difference="0.0")
        read.check_entry(ent, [row], True)
        with self.assertRaises(read.Stop):
            read.check_entry(ent, [dict(row, no_base_line="0")], False)
        with self.assertRaises(read.Stop):  # the median the scorer wrote must be the listed differences'
            read.check_entry(ent, [dict(row, median_difference="0.01")], True)


class Rulings(unittest.TestCase):
    """The rulings of 20:03Z against tables() and summary(), on SYNTHETIC count rows."""

    def tables(self, rows, **over):
        a = dict(line_50="beside", twin_rule="beside", dataset_differences="yes",
                 names=str(Path(TMP.name) / "template" / "names.tsv"))
        a.update(over)
        ns = argparse.Namespace(**a)
        out = read.tables(ns, [], rows, {}, ["10", "15", "20", "25", "25-s45", "50"])
        text = read.summary(argparse.Namespace(measures_over="shared", **a), [], [], out[0], out[7], out[6])
        return out, text

    def test_ruling_3_a_twin_that_fails_makes_a_line_not_shown_to_be_within_the_margin(self):
        (dec, twin, *_), text = self.tables(line_rows("10", 1) + line_rows("15", 4) + line_rows("25-s45", 2)
                                            + line_rows("50", 5))
        reading = {r["line"]: r["reading"] for r in dec if r["row"] == "line"}
        self.assertEqual(reading["10"], read.NOT_CANDIDATE + ": " + read.NOT_SHOWN)
        self.assertEqual(reading["15"], read.CANDIDATE + "; the twin: 2 of 5 learners meet")
        t = next(r for r in twin if r["row"] == "line")
        self.assertEqual(t["reading"], "the count with no change; never a candidate; " + read.TWIN_FAILS)
        self.assertIn("twin 25-s45: 2 of 5 learners meet; the count with no change, never a candidate; "
                      + read.TWIN_FAILS, text)
        self.assertNotIn("worse", text)

    def test_ruling_3_a_twin_that_meets_holds_the_row_and_is_never_a_candidate(self):
        (dec, twin, *_), text = self.tables(line_rows("10", 1) + line_rows("15", 4) + line_rows("25-s45", 3)
                                            + line_rows("50", 5))
        reading = {r["line"]: r["reading"] for r in dec if r["row"] == "line"}
        self.assertEqual(reading["10"], read.NOT_CANDIDATE)
        self.assertEqual(reading["15"], read.CANDIDATE + "; the twin: 3 of 5 learners meet")
        t = next(r for r in twin if r["row"] == "line")
        self.assertEqual((t["learners_meeting"], t["three_or_more_meet"]), ("3", "1"))
        self.assertNotIn("candidate", t)
        self.assertEqual(t["reading"], "the count with no change; never a candidate; " + read.TWIN_HOLDS)
        twin_lines = [x for x in text.splitlines() if x.startswith("twin ")]
        self.assertEqual(twin_lines, ["twin 25-s45: 3 of 5 learners meet; the count with no change, never a "
                                      "candidate; " + read.TWIN_HOLDS])

    def test_ruling_10_a_learner_that_does_not_meet_names_the_part_that_failed(self):
        rows = line_rows("10", 1)
        rows[1].update(interval_hi="0.02")  # cart: the upper end at 0.02 is not wholly below it
        rows[2].update(interval_hi="0.019999999999", covered="7", coverage_n="10")  # rf: 0.7 against 0.9
        rows[3].update(interval_hi="0.019999999999", covered="2", coverage_n="5", base_covered="11",
                       base_coverage_n="20")  # gbt: 0.4 against 0.55 is 0.15 apart
        rows[4].update(entered="0", not_scored="4", interval_hi="")  # lasso: no dataset entered
        (dec, *_), _ = self.tables(rows + line_rows("25-s45", 3))
        got = {r["learner"]: r["reading"] for r in dec if r["row"] == "learner"}
        self.assertEqual(got, {"scm": "meets", "cart": "does not meet: the interval's upper end",
                               "rf": "does not meet: the coverage", "gbt": "does not meet: the coverage",
                               "lasso": "does not meet: no dataset entered"})
        both = line_rows("15", 0)
        both[0].update(covered="1", coverage_n="10")
        (dec, *_), _ = self.tables(both + line_rows("25-s45", 3))
        self.assertEqual(dec[0]["reading"], "does not meet: the interval's upper end and the coverage")
        # A coverage exactly 0.05 from 25's is within it, as the scorer compares it.
        edge = line_rows("20", 1)
        edge[0].update(covered="17", coverage_n="20", base_covered="18", base_coverage_n="20")
        (dec, *_), _ = self.tables(edge + line_rows("25-s45", 3))
        self.assertEqual(dec[0]["reading"], "meets")

    def test_ruling_10_a_meets_cell_that_its_cells_do_not_give_stops(self):
        rows = line_rows("10", 1)
        rows[0].update(interval_hi="0.03")
        with self.assertRaises(read.Stop):
            self.tables(rows + line_rows("25-s45", 3))
        rows = line_rows("10", 0)
        rows[0].update(interval_hi="0.01")
        with self.assertRaises(read.Stop):
            self.tables(rows + line_rows("25-s45", 3))

    def test_the_coverage_fractions_are_the_shared_results_at_the_line_less_25s(self):
        rows = line_rows("10", 0, covered=(9, 10), base=(7, 10))
        rows[0].update(base_covered_own="5", base_coverage_n_own="10")
        (dec, *_), _ = self.tables(rows + line_rows("25-s45", 3))
        self.assertEqual((dec[0]["coverage_at_line"], dec[0]["coverage_at_25"], dec[0]["coverage_minus_25"]),
                         ("0.900000", "0.700000", "0.200000"))

    def test_the_twin_is_read_on_its_levelled_row_alone(self):
        rows = line_rows("10", 1) + line_rows("25-s45", 1)
        rows += [dict(r, coverage="fold_range") for r in line_rows("25-s45", 4)]
        for rule in ("void", "beside"):
            (dec, twin, other, *rest), text = self.tables(rows, twin_rule=rule)
            self.assertFalse(rest[3], rule)  # void: the twin does not meet on its deciding row
            self.assertEqual(next(r for r in dec if r["row"] == "line")["reading"],
                             read.NOT_CANDIDATE + ": " + read.NOT_SHOWN, rule)
            self.assertIn(read.TWIN_FAILS, next(r for r in twin if r["row"] == "line")["reading"])

    def test_ruling_2_line_50_beside_is_no_candidate_in_any_cell_or_summary_line(self):
        (dec, twin, other, *_), text = self.tables(line_rows("10", 1) + line_rows("25-s45", 3) + line_rows("50", 5))
        self.assertFalse([r for r in dec + twin if r["line"] == "50"])
        fifty = [r for r in other if r["line"] == "50"]
        self.assertEqual(len(fifty), 6)
        self.assertTrue(all(r["decides"] == "no: line 50 is not a lower line" and "candidate" not in r for r in fifty))
        self.assertFalse([x for x in text.splitlines() if "50" in x and "candidate" in x
                          and "never a candidate" not in x])


class Gates(Copy):
    def test_a_scorer_with_another_sha256_stops(self):
        self.stops(args(self.root, **{"--scorer-sha256": "0" * 64}), "its sha256 is not the given")

    def test_a_failed_species_stops(self):
        (self.root / "back" / "synthia_beta" / "FAILED").write_text("FAILED synthia_beta at fit-s42\n")
        self.stops(args(self.root), "synthia_beta: FAILED is present")

    def test_no_done_stops(self):
        (self.root / "back" / "synthia_beta" / "DONE").unlink()
        self.stops(args(self.root), "synthia_beta: no DONE")

    def test_a_wrong_code_sha256_stops(self):
        self.stops(args(self.root, **{"--code-sha256": "f" * 64}), f"DONE has code sha256 {CODE}, the given one is")

    def test_binary_txt_that_names_another_code_stops(self):
        first = self.root / "back" / "synthia_gamma"
        (first / "binary.txt").write_text(f"commit {'9' * 40}\ngrm {BINARY}\ncode {'d' * 64}\n")
        seal_first(first)
        self.stops(args(self.root), "binary.txt does not name DONE's binary and code")

    def test_a_first_copy_changed_after_its_sha256sums_stops(self):
        with open(self.root / "back" / "synthia_alpha" / "expect.tsv", "a") as f:
            f.write("")
        (self.root / "back" / "synthia_alpha" / "lists.txt").write_text("species\tsynthia_alpha\n")
        self.stops(args(self.root), "lists.txt is missing, not listed or differs from its sha256 in SHA256SUMS")

    def test_a_species_twice_stops(self):
        shutil.copytree(self.root / "back" / "synthia_beta", self.root / "back" / "synthia_beta.again")
        pairs = [(self.root / "back" / sp, self.root / "back2" / sp) for sp in SPECIES]
        pairs.append((self.root / "back" / "synthia_beta.again", self.root / "back2" / "synthia_beta"))
        self.stops(args(self.root, pairs=pairs, **{"--species": "4"}), "species synthia_beta twice")

    def test_a_missing_species_stops(self):
        pairs = [(self.root / "back" / sp, self.root / "back2" / sp) for sp in ("synthia_alpha", "synthia_gamma")]
        self.stops(args(self.root, pairs=pairs), "2 species, want 3")

    def test_the_dataset_and_result_counts_are_arguments(self):
        self.stops(args(self.root, **{"--datasets": "5"}), "4 datasets, want 5")
        self.stops(args(self.root, **{"--results": "47"}), "48 results, want 47")
        self.stops(args(self.root, **{"--per-dataset": "30"}), "gives 12 results a dataset, want 30")

    def test_a_second_scoring_that_differs_from_its_sha256sums_stops(self):
        with open(self.root / "back2" / "synthia_alpha" / "score" / "scores.tsv", "a") as f:
            f.write("\n")
        self.stops(args(self.root), "score/scores.tsv is missing or differs from its sha256 in SHA256SUMS")

    def test_a_file_in_score_that_sha256sums_does_not_list_stops(self):
        (self.root / "back2" / "synthia_alpha" / "score" / "extra.tsv").write_text("SYNTHETIC\n")
        self.stops(args(self.root), "score/ holds ['extra.tsv'], which SHA256SUMS does not list")

    def test_a_failed_check_line_stops(self):
        second = self.root / "back2" / "synthia_gamma"
        (second / "logs" / "score.log").write_text("checks: 11 of 12 results pass every check; 1 failed checks of 180\n")
        seal_second(second)
        self.stops(args(self.root), "11 of 12 results pass every check with 1 failed checks; want 12 of 12 and 0")

    def test_a_check_row_that_is_not_ok_stops(self):
        second = self.root / "back2" / "synthia_gamma"
        p = second / "score" / "checks.tsv"
        lines = p.read_text().splitlines(keepends=True)
        cells = lines[3].split("\t")
        cells[2] = "0"
        lines[3] = "\t".join(cells)
        p.write_text("".join(lines))
        seal_second(second)
        self.stops(args(self.root), "fails the check")

    def test_a_dataset_short_of_its_results_stops(self):
        first = self.root / "back" / "synthia_gamma"
        p = first / "expect.tsv"
        lines = p.read_text().splitlines(keepends=True)
        p.write_text("".join(lines[:-1]))
        seal_first(first)
        self.stops(args(self.root), "synthia_gamma: synthia_gamma__drug_d holds 11 results of the 12")

    def test_a_result_missing_from_the_second_scoring_stops(self):
        # expect.tsv lists 12 results; the second scoring lacks one result's 20 score rows.
        second = self.root / "back2" / "synthia_gamma"
        p = second / "score" / "scores.tsv"
        lines = p.read_text().splitlines(keepends=True)
        p.write_text("".join(lines[:-20]))
        seal_second(second)
        self.stops(args(self.root), "score/scores.tsv: 20 rows of expect.tsv's results missing, 0 rows not in them")

    def test_an_antibiotic_other_than_the_names_files_stops(self):
        names = (self.root / "names.tsv").read_text().replace("drug c", "drug x")
        (self.root / "names.tsv").write_text(names)
        self.stops(args(self.root), "synthia_beta__drug_c has antibiotic ['drug c'], the names file 'drug x'")

    def test_a_dataset_the_names_file_lacks_stops(self):
        names = [x for x in (self.root / "names.tsv").read_text().splitlines(keepends=True) if "drug_d" not in x]
        (self.root / "names.tsv").write_text("".join(names))
        self.stops(args(self.root), "the names file does not name synthia_gamma__drug_d under synthia_gamma")

    def test_a_count_that_is_not_the_scorers_count_of_its_lines_stops(self):
        # The first species, so the gate stops after one count.
        second = self.root / "back2" / "synthia_alpha"
        p = second / "score" / "count.tsv"
        lines = p.read_text().splitlines(keepends=True)
        cells = lines[1].split("\t")
        cells[-1] = "SYNTHETIC\n"
        lines[1] = "\t".join(cells)
        p.write_text("".join(lines))
        seal_second(second)
        self.stops(args(self.root), "count.tsv is not what the scorer's --count writes from its lines.tsv")

    def test_a_lines_row_missing_stops(self):
        second = self.root / "back2" / "synthia_beta"
        p = second / "score" / "lines.tsv"
        lines = p.read_text().splitlines(keepends=True)
        p.write_text("".join(lines[:-1]))
        seal_second(second)
        self.stops(args(self.root), "score/lines.tsv: 1 rows missing, 0 rows not expected")

    def test_an_output_that_exists_stops(self):
        (self.root / "out").mkdir()
        (self.root / "out" / "gates.txt").write_text("SYNTHETIC\n")
        rc, _, err = run(args(self.root))
        self.assertEqual(rc, 2)
        self.assertIn("gates.txt exists", err)
        self.assertEqual(os.listdir(self.root / "out"), ["gates.txt"])


if __name__ == "__main__":
    unittest.main()
