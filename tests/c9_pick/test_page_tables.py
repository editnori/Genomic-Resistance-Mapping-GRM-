"""scripts/c9-pick/page_tables.py: the gates and the tables of docs/served-pick-result.md.

Every reading, list row, species, antibiotic and dataset name here is
SYNTHETIC, made up in the test: three made-up species, six made-up
datasets, one of them in M. The scorer's own functions write the reading's
files, so they have the form a real reading has. One reading is built and
written once for the suite; each refusal test changes a copy of it and
names the gate that must stop it.
"""

import contextlib
import csv
import hashlib
import io
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.c9_pick import common

score = common.load("score")
page = common.load("page_tables")

L_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
             "in_M", "in_S", "first_run", "pairs_S+R", "labels_depend_on"]
# dataset, species, antibiotic, in M, in S, first run
L = [
    ("synthia_coli__drug_a", "Synthia coli", "drug a", 0, 1, 1),
    ("synthia_coli__drug_b", "Synthia coli", "drug b", 0, 0, 1),
    ("synthia_coli__drug_m", "Synthia coli", "drug m", 1, 1, 0),
    ("fakeria_maxima__drug_c", "Fakeria maxima", "drug c", 0, 1, 1),
    ("fakeria_maxima__drug_d", "Fakeria maxima", "drug d", 0, 0, 1),
    ("mockella_minor__drug_e", "Mockella minor", "drug e", 0, 0, 1),
]
NAMES = {d: (sp, ab) for d, sp, ab, *_ in L}
# Per dataset: the default's learner per draw, and the resistant genomes called resistant, of 40, per draw for the
# default, the inner pick and the twin; each draw also calls 45 of 60 susceptible genomes susceptible. A step of
# 4 in 40 moves balanced accuracy by 0.05. None is a draw whose split did not run. "study" is (default, inner).
SPEC = {
    "synthia_coli__drug_a": {"learners": ["rf"] * 3, "default": [28] * 3, "inner": [32] * 3, "twin": [28, 29, 28],
                             "study": (30, 26)},
    "synthia_coli__drug_b": {"learners": ["rf", "rf", "gbt"], "default": [30] * 3, "inner": [30] * 3, "twin": [30] * 3},
    "synthia_coli__drug_m": {"learners": ["rf"] * 3, "default": [30] * 3, "inner": [30] * 3, "twin": [30] * 3,
                             "study": (30, 30)},
    "fakeria_maxima__drug_c": {"learners": ["gbt"] * 3, "default": [25] * 3, "inner": [21] * 3, "twin": [26] * 3,
                               "study": (25, 29)},
    "fakeria_maxima__drug_d": {"learners": ["lasso"] * 3, "default": [33] * 3, "inner": [37] * 3, "twin": [33] * 3},
    "mockella_minor__drug_e": {"learners": ["cart", None, "cart"], "default": [30, None, 30], "inner": [34, None, 34],
                               "twin": [30] * 3},
}
# What each answer beside changes: the earlier tie turns drug a's gain into a loss, the default's tie removes drug
# d's gain, and the seed-42 default serves rf at drug b's third draw.
ANSWER_SPEC = {
    "read": {},
    "tie-earlier": {"synthia_coli__drug_a": {"inner": [24] * 3}},
    "tie-default": {"fakeria_maxima__drug_d": {"inner": [33] * 3}},
    "default-seed42": {"synthia_coli__drug_b": {"learners": ["rf"] * 3, "default": [30, 30, 31]}},
}


def spec(answer, d, change=None):
    s = dict(SPEC[d])
    s.update(ANSWER_SPEC[answer].get(d, {}))
    s.update((change or {}).get((answer, d), {}))
    return s


def picks(answer, d, i, learner):
    """SYNTHETIC pick counts of one inner draw over its 5 outer folds."""
    other = "scm" if learner != "scm" else "cart"
    if d == "mockella_minor__drug_e" and i == 2:
        served, same, no_pick = {learner: 4, other: 1}, 3, 1
    elif answer == "tie-earlier":
        served, same, no_pick = {learner: 2, other: 3}, 2, 0
    else:
        served, same, no_pick = {learner: 3, other: 2}, 3, 0
    out = {"folds": 5, "ties": 1, "no_pick": no_pick, "same_as_default": same}
    out.update({"picks_" + l: served.get(l, 0) for l in common.LEARNERS})
    return out


def row(d, split, seed, role, arm, tp, learner="", pick=None):
    sp, ab = NAMES[d]
    base = {"dataset": d, "species": sp, "antibiotic": ab, "split": split, "seed": seed, "role": role, "arm": arm}
    if tp is None:
        return dict(base, scored=0, why="split unavailable")
    c = {"n": 100, "tp": tp, "fp": 15, "tn": 45, "fn": 40 - tp, "loss_sum": 40.0 - tp / 2}
    r = dict(base, scored=1, why="", learner=learner, **score.measures(c))
    if arm == "inner":
        r.update(pick)
    return r


def draws_of(answer, d, change=None):
    """One dataset's draws.tsv rows, in the order score.py writes them."""
    s = spec(answer, d, change)
    out = []
    for i, seed in enumerate((42, 43, 44)):
        learner, dflt, inner = s["learners"][i], s["default"][i], s["inner"][i]
        if dflt is None:
            out += [row(d, "lineage", seed, "draw", arm, None) for arm in ("default", "inner")]
            continue
        out.append(row(d, "lineage", seed, "draw", "default", dflt, learner=learner))
        out.append(row(d, "lineage", seed, "draw", "inner", inner, pick=picks(answer, d, i, learner)))
    for i, seed in enumerate((45, 46, 47)):
        out.append(row(d, "lineage", seed, "twin", "twin", s["twin"][i], learner=s["learners"][0]))
    if "study" in s:
        learner = s["learners"][0]
        out.append(row(d, "study", 42, "study", "default", s["study"][0], learner=learner))
        out.append(row(d, "study", 42, "study", "inner", s["study"][1], pick=picks(answer, d, 0, learner)))
    return out


def write_species(folder, answer, datasets, change=None):
    folder = Path(folder)
    shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True)
    draws = [r for d in sorted(datasets) for r in draws_of(answer, d, change)]
    score.write_tsv(folder / "draws.tsv", draws, score.DRAW_COLUMNS)
    score.write_tsv(folder / "datasets.tsv", score.dataset_rows(draws), score.DATASET_COLUMNS)
    score.write_tsv(folder / "picks.tsv", [{k: r.get(k) for k in page.PICK_COLUMNS} for r in draws if r["arm"] == "inner"],
                    page.PICK_COLUMNS)
    (folder / "options.txt").write_text(page.options_text(page.ANSWERS[answer]))


def write_l(path):
    with open(path, "w") as f:
        f.write("\t".join(L_COLUMNS) + "\n")
        for d, sp, ab, m, s, first in L:
            cells = [d, d.split("__")[0], sp, ab, 150, 60, 90, 30, m, s, first, 0, "none"]
            f.write("\t".join(str(c) for c in cells) + "\n")


def build(top):
    """A SYNTHETIC reading of the first run: top/reading and top/l.tsv."""
    top = Path(top)
    write_l(top / "l.tsv")
    reading = top / "reading"
    for answer in page.ANSWERS:
        paths = []
        for slug in ("fakeria_maxima", "mockella_minor", "synthia_coli"):
            datasets = [d for d, *_, first in L if d.startswith(slug) and first]
            write_species(reading / answer / "score" / slug, answer, datasets)
            paths.append(str(reading / answer / "score" / slug / "datasets.tsv"))
        with contextlib.redirect_stdout(io.StringIO()):
            assert score.main(["--count", str(reading / answer / "count"), *paths]) == 0
    (reading / "reading.txt").write_text("SYNTHETIC reading: no command ran\n")
    return reading, top / "l.tsv"


def table(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


TOP = OUT = None


def setUpModule():
    global TOP, OUT
    TOP = Path(tempfile.mkdtemp())
    build(TOP)
    OUT = TOP / "out"
    page.write(TOP / "reading", TOP / "l.tsv", OUT)


def tearDownModule():
    shutil.rmtree(TOP, ignore_errors=True)


def at(rows, **keys):
    got = [r for r in rows if all(r[k] == v for k, v in keys.items())]
    assert len(got) == 1, (keys, len(got))
    return got[0]


class Gates(unittest.TestCase):
    def copy(self):
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp, True)
        shutil.copytree(TOP / "reading", tmp / "reading")
        shutil.copy(TOP / "l.tsv", tmp / "l.tsv")
        return tmp

    def assert_stops(self, tmp, gate):
        with self.assertRaises(page.Stop) as e:
            page.write(tmp / "reading", tmp / "l.tsv", tmp / "out")
        self.assertTrue(str(e.exception).startswith(f"gate {gate} "), str(e.exception))
        self.assertFalse((tmp / "out").exists(), "nothing is written when a gate stops")

    def test_a_whole_reading_passes_every_gate_and_writes_its_files(self):
        files = sorted(p.relative_to(OUT).as_posix() for p in OUT.rglob("*") if p.is_file())
        want = ["README.md", "SHA256SUMS", "gates.txt", "page-conditions.tsv", "page-count.tsv", "page-datasets.tsv",
                "page-pick-pairs.tsv", "page-picks.tsv", "page-pooled.tsv", "page-species.tsv", "reading.txt"]
        want += [f"count/{a}/{n}" for a in page.ANSWERS for n in ("conditions.tsv", "count.tsv")]
        self.assertEqual(files, sorted(want))
        sums = (OUT / "SHA256SUMS").read_text().splitlines()
        self.assertEqual(len(sums), len(want) - 1)
        for line in sums:
            digest, name = line.split("  ")
            self.assertEqual(hashlib.sha256((OUT / name).read_bytes()).hexdigest(), digest, name)
        self.assertEqual(len((OUT / "gates.txt").read_text().splitlines()), 8)
        self.assertIn("first run's 5 datasets in 3 species, 2 of them with a study split", (OUT / "README.md").read_text())

    def test_main_prints_its_gate_count_only(self):
        tmp = self.copy()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            rc = page.main(["write", str(tmp / "reading"), str(tmp / "l.tsv"), str(tmp / "out")])
        self.assertEqual((rc, buf.getvalue()), (0, "gates: 8 of 8 pass; wrote 19 files\n"))
        with contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(page.main(["write", str(tmp / "reading"), str(tmp / "l.tsv"), str(tmp / "out")]), 1)

    def test_a_species_missing_from_the_reading_is_refused(self):
        # Planted fault 9: a species dropped.
        tmp = self.copy()
        for a in page.ANSWERS:
            shutil.rmtree(tmp / "reading" / a / "score" / "mockella_minor")
        self.assert_stops(tmp, 2)

    def test_an_answer_scored_under_other_options_is_refused(self):
        # Planted fault 8b: the read folder holds the earlier tie's scoring and count whole.
        tmp = self.copy()
        shutil.rmtree(tmp / "reading" / "read")
        shutil.copytree(tmp / "reading" / "tie-earlier", tmp / "reading" / "read")
        self.assert_stops(tmp, 3)

    def test_a_dataset_of_m_is_refused(self):
        # Planted fault 10: a dataset of M counted, in every answer alike.
        tmp = self.copy()
        for a in page.ANSWERS:
            write_species(tmp / "reading" / a / "score" / "synthia_coli", a,
                          ["synthia_coli__drug_a", "synthia_coli__drug_b", "synthia_coli__drug_m"])
        self.assert_stops(tmp, 4)

    def test_a_missing_draw_is_refused(self):
        tmp = self.copy()
        path = tmp / "reading" / "read" / "score" / "fakeria_maxima" / "draws.tsv"
        lines = path.read_text().splitlines(keepends=True)
        path.write_text("".join(l for l in lines if "\t46\ttwin\t" not in l or "drug_d" not in l))
        self.assert_stops(tmp, 5)

    def test_a_dataset_figure_other_than_its_draws_mean_is_refused(self):
        tmp = self.copy()
        path = tmp / "reading" / "read" / "score" / "fakeria_maxima" / "datasets.tsv"
        rows = table(path)
        rows[0]["balanced_accuracy"] = repr(float(rows[0]["balanced_accuracy"]) + 0.001)
        score.write_tsv(path, rows, score.DATASET_COLUMNS)
        self.assert_stops(tmp, 5)

    def test_a_tie_answer_whose_default_differs_from_the_read_ones_is_refused(self):
        tmp = self.copy()
        change = {("tie-default", "synthia_coli__drug_b"): {"default": [30, 31, 30]}}
        write_species(tmp / "reading" / "tie-default" / "score" / "synthia_coli", "tie-default",
                      ["synthia_coli__drug_a", "synthia_coli__drug_b"], change)
        self.assert_stops(tmp, 6)

    def test_a_count_other_than_the_scorers_is_refused(self):
        tmp = self.copy()
        path = tmp / "reading" / "default-seed42" / "count" / "count.tsv"
        rows = table(path)
        rows[0]["better"] = str(int(rows[0]["better"]) + 1)
        score.write_tsv(path, rows, score.COUNT_COLUMNS)
        self.assert_stops(tmp, 7)

    def test_an_out_dir_that_holds_a_file_is_refused(self):
        tmp = self.copy()
        (tmp / "out").mkdir()
        (tmp / "out" / "old.txt").write_text("x")
        with self.assertRaises(page.Stop):
            page.write(tmp / "reading", tmp / "l.tsv", tmp / "out")


class Tables(unittest.TestCase):
    def test_the_conditions_and_the_count_are_the_scorers(self):
        conds = table(OUT / "page-conditions.tsv")
        for a in page.ANSWERS:
            got = [r[a] for r in conds]
            self.assertEqual(got, [r["met"] for r in table(OUT / "count" / a / "conditions.tsv")], a)
        count = table(OUT / "page-count.tsv")
        for a in page.ANSWERS:
            mine = [{k: v for k, v in r.items() if k not in ("answer", "meets_1", "meets_2")} for r in count
                    if r["answer"] == a]
            self.assertEqual(mine, table(OUT / "count" / a / "count.tsv"), a)
        primary = at(count, answer="read", split="lineage", pair="inner-default", measure="balanced_accuracy")
        # Lineage, inner less default: drug a +0.05, drug b 0, drug c -0.05, drug d +0.05; drug e unscored.
        self.assertEqual((primary["better"], primary["worse"], primary["same"], primary["not_scored"]), ("2", "1", "1", "1"))

    def test_the_inner_figure_is_the_inner_picks_and_the_twin_stands_apart(self):
        # Planted fault 7: a twin row read as the arm.
        r = at(table(OUT / "page-datasets.tsv"), dataset="synthia_coli__drug_a", split="lineage")
        self.assertEqual((r["default"], r["inner"], r["difference"]), ("0.725", "0.775", "0.05"))
        self.assertAlmostEqual(float(r["twin"]), (0.725 * 2 + 0.7375) / 3, places=12)
        self.assertAlmostEqual(float(r["twin_difference"]), 0.0125 / 3, places=12)

    def test_the_read_answer_is_the_reading_and_the_beside_answers_stand_apart(self):
        # Planted fault 8: a beside answer read as the read one.
        r = at(table(OUT / "page-datasets.tsv"), dataset="synthia_coli__drug_a", split="lineage")
        self.assertEqual(r["inner"], "0.775")
        sp = at(table(OUT / "page-species.tsv"), species_slug="synthia_coli")
        self.assertEqual((sp["lineage_inner_better"], sp["lineage_inner_worse"]), ("1", "0"))
        earlier = at(table(OUT / "page-count.tsv"), answer="tie-earlier", split="lineage", pair="inner-default",
                     measure="balanced_accuracy")
        self.assertEqual((earlier["better"], earlier["worse"]), ("1", "2"))
        picks = table(OUT / "page-picks.tsv")
        self.assertEqual(at(picks, answer="read", split="lineage", arm="inner")["same_as_default"], "42")
        self.assertEqual(at(picks, answer="tie-earlier", split="lineage", arm="inner")["same_as_default"], "29")

    def test_a_study_row_stays_on_the_study_split(self):
        # Planted fault 11: drug a gains 0.05 on the lineage split and loses 0.05 on the study split.
        sp = at(table(OUT / "page-species.tsv"), species_slug="synthia_coli")
        self.assertEqual([sp[f"lineage_inner_{k}"] for k in ("entered", "better", "worse", "same", "median")],
                         ["2", "1", "0", "1", "0.025"])
        self.assertEqual([sp[f"study_inner_{k}"] for k in ("entered", "better", "worse", "same", "median")],
                         ["1", "0", "1", "0", "-0.05"])
        fk = at(table(OUT / "page-species.tsv"), species_slug="fakeria_maxima")
        self.assertEqual([fk[f"lineage_inner_{k}"] for k in ("better", "worse", "same")], ["1", "1", "0"])
        self.assertEqual([fk[f"study_inner_{k}"] for k in ("better", "worse", "same")], ["1", "0", "0"])
        study = at(table(OUT / "page-datasets.tsv"), dataset="synthia_coli__drug_a", split="study")
        self.assertEqual((study["draws"], study["difference"], study["twin"]), ("1", "-0.05", ""))

    def test_a_draw_not_scored_leaves_its_dataset_unscored_not_zero(self):
        # Planted fault 12: a missing draw read as zero.
        r = at(table(OUT / "page-datasets.tsv"), dataset="mockella_minor__drug_e", split="lineage")
        self.assertEqual((r["default"], r["inner"], r["difference"], r["twin_difference"]), ("", "", "", ""))
        self.assertEqual((r["draws"], r["inner_draws_scored"], r["default_learners"]), ("3", "2", "cart,cart"))
        self.assertEqual((r["folds"], r["folds_no_pick"], r["folds_differing"]), ("10", "1", "3"))
        self.assertEqual(r["twin"], "0.75")
        mm = at(table(OUT / "page-species.tsv"), species_slug="mockella_minor")
        self.assertEqual((mm["lineage_inner_entered"], mm["lineage_inner_not_scored"], mm["lineage_inner_median"]),
                         ("0", "1", ""))

    def test_every_species_and_antibiotic_is_named_in_full(self):
        sp = table(OUT / "page-species.tsv")
        self.assertEqual([r["species"] for r in sp], ["Fakeria maxima", "Mockella minor", "Synthia coli"])
        self.assertEqual([r["datasets"] for r in sp], ["2", "1", "2"])
        ds = table(OUT / "page-datasets.tsv")
        self.assertEqual(sorted({(r["species"], r["antibiotic"]) for r in ds}),
                         sorted((sp_, ab) for d, sp_, ab, m, s, first in L if first))
        self.assertEqual(len(ds), 5 + 2)

    def test_the_pick_tables_count_folds_per_learner(self):
        picks = table(OUT / "page-picks.tsv")
        inner = at(picks, answer="read", split="lineage", arm="inner")
        # 14 scored draws of 5 folds: drug e's second draw did not run.
        self.assertEqual((inner["draws"], inner["folds"], inner["no_pick"], inner["tied"]), ("14", "70", "1", "14"))
        self.assertEqual(int(inner["differs"]), 70 - 42 - 1)
        self.assertEqual(sum(int(inner[l]) for l in common.LEARNERS), 70)
        dflt = at(picks, answer="read", split="lineage", arm="default")
        self.assertEqual([dflt[l] for l in common.LEARNERS], ["0", "10", "25", "20", "15"])
        seed42 = at(picks, answer="default-seed42", split="lineage", arm="default")
        self.assertEqual((seed42["rf"], seed42["gbt"]), ("30", "15"))
        pairs = table(OUT / "page-pick-pairs.tsv")
        rf = at(pairs, split="lineage", default_learner="rf")
        self.assertEqual((rf["draws"], rf["folds"], rf["rf"], rf["scm"]), ("5", "25", "15", "10"))
        self.assertEqual(at(pairs, split="study", default_learner="gbt")["folds"], "5")

    def test_pooled_calls_carry_their_denominators(self):
        pooled = table(OUT / "page-pooled.tsv")
        inner = at(pooled, split="lineage", arm="inner")
        # Drugs a to d scored at 3 draws of 100 genomes; drug e is not scored.
        self.assertEqual((inner["datasets"], inner["scored"], inner["n"], inner["resistant"]), ("5", "4", "1200", "480"))
        self.assertEqual(inner["tp"], str(3 * (32 + 30 + 21 + 37)))
        self.assertAlmostEqual(float(inner["recall"]), 3 * (32 + 30 + 21 + 37) / 480, places=12)
        self.assertLess(float(inner["recall_lo"]), float(inner["recall"]))
        twin = at(pooled, split="lineage", arm="twin")
        self.assertEqual((twin["scored"], twin["n"]), ("5", "1500"))
        self.assertEqual({r["split"] + "/" + r["arm"] for r in pooled},
                         {"lineage/default", "lineage/inner", "lineage/twin", "study/default", "study/inner"})


class Printing(unittest.TestCase):
    def test_print_gives_each_table_with_names_in_full_and_no_genome(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.assertEqual(page.main(["print", str(OUT)]), 0)
        text = buf.getvalue()
        for name in ("conditions", "twin", "count", "beside", "pooled", "picks", "pairs", "species", "datasets"):
            self.assertIn(f"[{name}]\n| ", text)
        self.assertIn("| *Synthia coli* | drug a | lineage | `rf`, `rf`, `rf` | 0.725 | 0.775 | 0.050 | 0.004 | 3 of 3 |",
                      text)
        self.assertIn("| *Mockella minor* | drug e | lineage | `cart`, `cart` | none | none | none | none | 2 of 3 | 1 |",
                      text)
        self.assertIn("| *Fakeria maxima* | 2 | 1 | 1 / 1 / 0 |", text)
        self.assertIn("| inner pick less default | ", text)
        self.assertNotRegex(text, r"\d+\.\d+\.\d+")

    def test_each_page_table_prints_from_the_out_dir_alone(self):
        got = page.tables(OUT)
        conditions = got["conditions"]
        self.assertEqual(len(conditions), 2 + 6)
        self.assertTrue(conditions[-1].endswith("| the maintainer's | the maintainer's | the maintainer's | the maintainer's |"))
        self.assertEqual(len(got["count"]), 2 + 3 * len(score.MEASURES))
        self.assertEqual(len(got["datasets"]), 2 + 7)


if __name__ == "__main__":
    unittest.main()
