"""scripts/mechanism/year_mix.py: the mix of catalogue strata by collection year (docs/mechanism-year-result.md).

Every report, genome id, group name, year, score and count here is made up
in the test; the species and antibiotic names are shapes only. The answers
were worked by hand from the rows below.
"""

import contextlib
import csv
import io
import json
import os
import tempfile
import unittest

from tests.mechanism import common as c

ym = c.load("year_mix")

SPECIES = "Synthia coli"
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
GROUPS = ["synA family", "synB", "synC (synC-1 family)"]
COUNTS_FOR_RULE = [True, False, True]
LISTS_HEADER = ["dataset", "species_slug", "species", "antibiotic", "in_M", "in_A"]
PUBLISHED_ROWS = ("scm", "cart", "rf", "gbt", "lasso", "catalog", "program")
# point_only, gene_only, both, neither: (resistant, susceptible)
MISMATCH = {"point_only": (5, 4), "gene_only": (3, 1), "both": (1, 0), "neither": (3, 3)}
MATCH = {"point_only": (0, 0), "gene_only": (3, 1), "both": (0, 0), "neither": (3, 3)}

# (genome id, resistant, reader's call, groups, year, score of all five or the five scores)
DRUG_A = [
    ("9999.1001", True, True, [0], "2015", 0.9),
    ("9999.1002", True, False, [1], "2015", 0.5),
    ("9999.1003", True, False, [], "2016", 0.2),
    ("9999.1004", False, False, [], "2016", 0.1),
    ("9999.1005", False, False, [2], "2015", 0.6),
    ("9999.1006", False, False, [], "", 0.4),
    ("9999.1007", True, True, [0, 1], "2015-2016", 0.7),
]
DRUG_B = [
    ("9999.1001", True, True, [0], "2015", 0.8),
    ("9999.1008", False, False, [], "2016", 0.3),
    ("9999.1009", True, False, [], "2016", [0.0, 0.6, 0.6, 0.6, 0.7]),
]
OTHER = [("9999.2001", True, True, [0], "2015", 0.9), ("9999.2002", False, False, [], "2016", 0.1)]
IDS = sorted({r[0] for rows in (DRUG_A, DRUG_B, OTHER) for r in rows})


def slug(antibiotic):
    return "synthia_coli__" + antibiotic.replace(" ", "_")


def make_report(antibiotic, rows, catalog=True, lineage=True, species=SPECIES):
    samples = [{"id": g, "resistant": r, "catalog_resistant": call, "groups": list(gr), "year": y, "lineage": 1,
                "country": "Nowhere", "host": "none", "mlst": "1", "mic": None, "from_mic": False}
               for g, r, call, gr, y, _ in rows]
    per = [s if isinstance(s, list) else [s] * 5 for *_, s in rows]
    learners = [{"learner": name, "scores": [p[j] for p in per], "metrics": {}} for j, name in enumerate(LEARNERS)]
    learners.append({"learner": "catalog", "scores": [1.0 if r[2] else 0.0 for r in rows], "metrics": {}})
    splits = [{"scheme": "random", "learners": [dict(lr, scores=list(reversed(lr["scores"]))) for lr in learners]}]
    if lineage:
        splits.append({"scheme": "lineage", "learners": learners})
    rep = {"id": slug(antibiotic), "species": species, "antibiotic": antibiotic, "samples": samples,
           "catalog_groups": list(GROUPS), "splits": splits, "catalog": None}
    if catalog:
        tp = sum(1 for r in rows if r[1] and r[2])
        fn = sum(1 for r in rows if r[1] and not r[2])
        fp = sum(1 for r in rows if not r[1] and r[2])
        tn = sum(1 for r in rows if not r[1] and not r[2])
        rep["catalog"] = {"version": "synthetic", "rule": [["C:SYN"]],
                          "groups": [{"group": g, "counts_for_rule": f, "class": "SYN", "subclass": "SYN",
                                      "n_resistant": 0, "n_susceptible": 0} for g, f in zip(GROUPS, COUNTS_FOR_RULE)],
                          "metrics": {"tp": tp, "fp": fp, "tn": tn, "fn": fn, "n": len(rows)}}
    return rep


def world():
    """Six reports: two primary, one in M, one outside A, one not in the list,
    one primary with no lineage split."""
    reports = {
        slug("drug a"): make_report("drug a", DRUG_A),
        slug("drug b"): make_report("drug b", DRUG_B),
        slug("drug m"): make_report("drug m", OTHER),
        slug("drug n"): make_report("drug n", OTHER),
        slug("drug x"): make_report("drug x", OTHER),
        slug("drug c"): make_report("drug c", OTHER, lineage=False),
    }
    lists = [(slug(a), "1", "0") for a in ("drug a", "drug b", "drug c")]
    lists += [(slug("drug m"), "1", "1"), (slug("drug n"), "0", "0")]
    return reports, lists


def thick():
    """One primary dataset: 100 genomes in 2018 (10 resistant the reader calls,
    7 of them at a mean of 0.9) and 99 in 2019 (9 resistant it misses, all at 0.9)."""
    rows = [(f"9999.3{i:03d}", True, True, [0], "2018", 0.9 if i < 7 else 0.1) for i in range(10)]
    rows += [(f"9999.4{i:03d}", False, False, [], "2018", 0.1) for i in range(90)]
    rows += [(f"9999.5{i:03d}", True, False, [], "2019", 0.9) for i in range(9)]
    rows += [(f"9999.6{i:03d}", False, False, [], "2019", 0.1) for i in range(90)]
    return {slug("drug t"): make_report("drug t", rows)}, [(slug("drug t"), "1", "0")]


def build(tmp, reports, lists, published=MISMATCH, differ=False):
    rdir = os.path.join(tmp, "datasets")
    os.makedirs(rdir)
    for stem, rep in reports.items():
        with open(os.path.join(rdir, stem + ".json"), "w", encoding="utf-8") as f:
            f.write(rep if isinstance(rep, str) else json.dumps(rep))
    lpath = os.path.join(tmp, "lists.tsv")
    with open(lpath, "w", encoding="utf-8") as f:
        f.write("\t".join(LISTS_HEADER) + "\n")
        for stem, a, m in lists:
            antibiotic = stem.split("__")[1].replace("_", " ")
            f.write("\t".join([stem, "synthia_coli", SPECIES, antibiotic, m, a]) + "\n")
    ppath = os.path.join(tmp, "page-m1-strata.tsv")
    with open(ppath, "w", encoding="utf-8") as f:
        f.write("set\tsplit\trow\tstratum\tkind\tresistant\tsusceptible\n")
        for row in PUBLISHED_ROWS:
            for s in ym.FOUR:
                r, n = published[s]
                if differ and row == "lasso" and s == "neither":
                    r += 1
                f.write(f"primary\tlineage\t{row}\t{s}\tstratum\t{r}\t{n}\n")
                f.write(f"M\tlineage\t{row}\t{s}\tstratum\t99\t99\n")
    return rdir, lpath, ppath


def run(*argv):
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        rc = ym.main(list(argv))
    return rc, out.getvalue()


def tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return [list(r) for r in csv.reader(f, delimiter="\t")]


class Rules(unittest.TestCase):
    def test_a_year_not_four_digits_is_unknown(self):
        self.assertEqual(ym.year_key("2015"), "2015")
        self.assertEqual(ym.year_key("1999"), "1999")
        for text in ("2015-2016", "", None, "15", " 2015", "2015 ", "2015\n", "２０１５", "2015.0",
                     "20150", 2015):
            with self.subTest(text=text):
                self.assertEqual(ym.year_key(text), "year_unknown")

    def test_a_mean_of_exactly_half_is_called(self):
        self.assertTrue(ym.mean_called([0.5] * 5))
        self.assertTrue(ym.mean_called([1.0, 1.0, 0.5, 0.0, 0.0]))
        self.assertTrue(ym.mean_called([0.0, 0.6, 0.6, 0.6, 0.7]))
        self.assertTrue(ym.mean_called([0.7, 0.6, 0.6, 0.6, 0.0]))
        self.assertFalse(ym.mean_called([0.5, 0.5, 0.5, 0.5, 0.4999]))

    def test_the_call_is_the_mean_of_the_five(self):
        self.assertFalse(ym.mean_called([0.9, 0.0, 0.0, 0.0, 0.0]))
        self.assertTrue(ym.mean_called([0.0, 0.6, 0.6, 0.6, 0.7]))

    def test_the_plan_prints_the_columns_the_writer_writes(self):
        plan = ym.plan()
        self.assertEqual(plan["learners"], list(LEARNERS))
        self.assertEqual((plan["split"], plan["cutoff"], plan["thin_under_genomes"]), ("lineage", 0.5, 100))
        with tempfile.TemporaryDirectory() as tmp:
            rdir, lpath, ppath = build(tmp, *world())
            reading, out = os.path.join(tmp, "reading.json"), os.path.join(tmp, "out")
            run("read", "--reports", rdir, "--lists", lpath, "--published", ppath, "--out", reading)
            run("write", "--reading", reading, "--reports", rdir, "--out", out)
            for name, header in plan["columns"].items():
                self.assertEqual(tsv(os.path.join(out, name))[0], header, name)


class Read(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.reading = os.path.join(self.tmp.name, "reading.json")
        self.out = os.path.join(self.tmp.name, "out")

    def go(self, reports, lists, published=MISMATCH, differ=False):
        self.rdir, lpath, ppath = build(self.tmp.name, reports, lists, published, differ)
        return run("read", "--reports", self.rdir, "--lists", lpath, "--published", ppath, "--out", self.reading)

    def write(self):
        return run("write", "--reading", self.reading, "--reports", self.rdir, "--out", self.out)

    def table(self, name):
        return tsv(os.path.join(self.out, name))

    def loaded(self):
        with open(self.reading, encoding="utf-8") as f:
            return json.load(f)

    def test_a_rebuild_that_does_not_match_falls_back_to_two_groups(self):
        rc, text = self.go(*world())
        self.assertEqual(rc, 0, text)
        self.assertEqual(self.write()[0], 0)
        self.assertEqual(self.table("step1.tsv"), [
            ym.STEP1_HEADER,
            ["point_only", "5", "4", "0", "0", "no"],
            ["gene_only", "3", "1", "3", "1", "yes"],
            ["both", "1", "0", "0", "0", "no"],
            ["neither", "3", "3", "3", "3", "yes"],
            ["all_strata", "12", "8", "6", "4", "no"],
        ])
        self.assertIn("step 1: point_only published 5 resistant 4 susceptible; rebuilt 0 resistant 0 susceptible;"
                      " equal no", text)
        self.assertIn("the year read runs on the reader's call and no call", text)
        self.assertEqual(self.loaded()["mode"], "two")
        self.assertEqual([r[1] for r in self.table("years.tsv")[1:4]], ["reader_call", "reader_no_call",
                                                                         "reader_no_rule"])

    def test_the_two_group_tables_stand_as_worked_by_hand(self):
        self.assertEqual(self.go(*world())[0], 0)
        self.assertEqual(self.write()[0], 0)
        self.assertEqual(self.table("years.tsv")[1:], [
            ["2015", "reader_call", "2", "2", "0", "2", "0", "", "yes"],
            ["2015", "reader_no_call", "2", "1", "1", "1", "1", "", "yes"],
            ["2015", "reader_no_rule", "0", "0", "0", "0", "0", "", "yes"],
            ["2016", "reader_call", "0", "0", "0", "0", "0", "", "yes"],
            ["2016", "reader_no_call", "4", "2", "2", "1", "0", "", "yes"],
            ["2016", "reader_no_rule", "0", "0", "0", "0", "0", "", "yes"],
            ["year_unknown", "reader_call", "1", "1", "0", "1", "0", "", "yes"],
            ["year_unknown", "reader_no_call", "1", "0", "1", "0", "0", "", "yes"],
            ["year_unknown", "reader_no_rule", "0", "0", "0", "0", "0", "", "yes"],
        ])
        self.assertEqual(self.table("totals.tsv")[1:], [
            ["reader_call", "3", "3", "0", "3", "0", ""],
            ["reader_no_call", "7", "3", "4", "2", "1", ""],
            ["reader_no_rule", "0", "0", "0", "0", "0", ""],
        ])
        self.assertEqual(self.table("shares.tsv")[1:4], [
            ["2015", "reader_call", "2", "3", "0.666666666667", "0", "1", "0.000000000000", "4", "yes"],
            ["2015", "reader_no_call", "1", "3", "0.333333333333", "1", "1", "1.000000000000", "4", "yes"],
            ["2015", "reader_no_rule", "0", "3", "0.000000000000", "0", "1", "0.000000000000", "4", "yes"],
        ])
        self.assertEqual(self.table("species.tsv"), [ym.SPECIES_HEADER])

    def test_a_rebuild_that_matches_runs_on_four_strata(self):
        rc, text = self.go(*world(), published=MATCH)
        self.assertEqual(rc, 0, text)
        self.assertNotIn("reader's call and no call", text)
        self.assertEqual(self.write()[0], 0)
        self.assertEqual([r[-1] for r in self.table("step1.tsv")[1:]], ["yes"] * 5)
        self.assertEqual(self.table("years.tsv")[1:5], [
            ["2015", "point_only", "0", "0", "0", "0", "0", "", "yes"],
            ["2015", "gene_only", "3", "2", "1", "2", "1", "", "yes"],
            ["2015", "both", "0", "0", "0", "0", "0", "", "yes"],
            ["2015", "neither", "1", "1", "0", "1", "0", "", "yes"],
        ])
        self.assertEqual(self.table("totals.tsv")[1:], [
            ["point_only", "0", "0", "0", "0", "0", ""],
            ["gene_only", "4", "3", "1", "3", "1", ""],
            ["both", "0", "0", "0", "0", "0", ""],
            ["neither", "6", "3", "3", "2", "0", ""],
        ])

    def test_the_rebuild_reads_groups_that_count_for_the_rule(self):
        # 9999.1002 is resistant with a group that does not count: neither.
        # 9999.1005 is susceptible with a group that counts: gene only.
        self.assertEqual(self.go(*world())[0], 0)
        step = {r["stratum"]: r for r in self.loaded()["step1"]}
        self.assertEqual((step["gene_only"]["rebuilt_resistant"], step["gene_only"]["rebuilt_susceptible"]), (3, 1))
        self.assertEqual((step["neither"]["rebuilt_resistant"], step["neither"]["rebuilt_susceptible"]), (3, 3))

    def test_the_call_reads_all_five_learners(self):
        # 9999.1009: scm 0.0 and the other four 0.6, 0.6, 0.6, 0.7, a mean of 0.5: called.
        self.assertEqual(self.go(*world())[0], 0)
        cell = [x for x in self.loaded()["cells"] if (x["antibiotic"], x["year"]) == ("drug b", "2016")]
        self.assertEqual([(x["stratum"], x["resistant"], x["resistant_called_mean5"]) for x in cell],
                         [("reader_no_call", 1, 1)])

    def test_a_genome_counts_once_per_dataset_it_appears_in(self):
        self.assertEqual(self.go(*world())[0], 0)
        totals = self.loaded()["totals"]
        self.assertEqual(sum(t["genomes"] for t in totals.values()), len(DRUG_A) + len(DRUG_B))

    def test_datasets_outside_the_primary_set_are_skipped_and_counted(self):
        rc, text = self.go(*world())
        self.assertEqual(rc, 0, text)
        reading = self.loaded()
        self.assertEqual(reading["used"], 2)
        self.assertEqual(reading["not_used"], {"in_M": 1, "not_in_A": 1, "not_in_list": 1, "no_lineage_split": 1})
        self.assertIn("guard 5: 2 used + 1 in_M + 1 not_in_A + 1 not_in_list + 1 no_lineage_split = 6 reports", text)
        self.assertEqual({x["antibiotic"] for x in reading["cells"]}, {"drug a", "drug b"})

    def test_cells_that_disagree_with_metrics_stop_the_read(self):
        reports, lists = world()
        reports[slug("drug b")]["catalog"]["metrics"]["tp"] += 1
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 2: 1 reports' cells differ from catalog.metrics: Synthia coli with drug b", text)
        self.assertFalse(os.path.exists(self.reading))

    def test_a_report_outside_the_primary_set_is_held_to_its_metrics_too(self):
        reports, lists = world()
        reports[slug("drug m")]["catalog"]["metrics"]["fn"] += 1
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("guard 2: 1 reports' cells differ from catalog.metrics: Synthia coli with drug m", text)

    def test_a_scores_length_mismatch_stops(self):
        reports, lists = world()
        reports[slug("drug a")]["splits"][1]["learners"][4]["scores"].pop()
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 3: Synthia coli with drug a, lasso: 6 scores for 7 samples", text)
        self.assertFalse(os.path.exists(self.reading))

    def test_a_score_that_is_not_a_number_from_0_to_1_stops(self):
        for bad in (1.5, -0.1, None, True, "0.4"):
            with self.subTest(bad=bad):
                reports, lists = world()
                reports[slug("drug b")]["splits"][1]["learners"][2]["scores"][0] = bad
                with tempfile.TemporaryDirectory() as tmp:
                    rdir, lpath, ppath = build(tmp, reports, lists)
                    rc, text = run("read", "--reports", rdir, "--lists", lpath, "--published", ppath, "--out",
                                   os.path.join(tmp, "r.json"))
                self.assertEqual(rc, 2)
                self.assertIn("guard 3: Synthia coli with drug b, rf: 1 scores are not a number from 0 to 1", text)

    def test_a_missing_learner_stops(self):
        reports, lists = world()
        del reports[slug("drug a")]["splits"][1]["learners"][1]
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("guard 3: Synthia coli with drug a has no cart row on the lineage split", text)

    def test_a_report_that_does_not_parse_stops(self):
        reports, lists = world()
        reports[slug("drug n")] = "{\"samples\": ["
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 1: 1 of 6 reports do not parse: synthia_coli__drug_n", text)

    def test_a_primary_dataset_with_no_report_stops(self):
        reports, lists = world()
        del reports[slug("drug c")]
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 2)
        self.assertIn("guard 5: 1 primary datasets of the list have no report: synthia_coli__drug_c", text)

    def test_a_null_call_is_its_own_group(self):
        reports, lists = world()
        rows = [(g, r, None, gr, y, s) for g, r, _, gr, y, s in DRUG_B]
        reports[slug("drug b")] = make_report("drug b", rows, catalog=False)
        rc, text = self.go(reports, lists)
        self.assertEqual(rc, 0, text)
        self.assertIn("1 with no catalogue row and no call", text)
        self.assertEqual(self.loaded()["totals"]["reader_no_rule"], {
            "genomes": 3, "resistant": 2, "susceptible": 1, "resistant_called_mean5": 2, "susceptible_called_mean5": 0})

    def test_published_rows_that_differ_between_learners_stop(self):
        rc, text = self.go(*world(), differ=True)
        self.assertEqual(rc, 2)
        self.assertIn("stop: step 1: the published stratum rows differ between the file's 7 rows", text)

    def test_year_counts_that_do_not_add_up_stop_the_write(self):
        self.assertEqual(self.go(*world())[0], 0)
        reading = self.loaded()
        reading["totals"]["reader_no_call"]["genomes"] += 1
        with open(self.reading, "w", encoding="utf-8") as f:
            json.dump(reading, f)
        rc, text = self.write()
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 4: reader_no_call genomes sums to 7 over the years, its total is 8", text)
        self.assertFalse(os.path.exists(os.path.join(self.out, "years.tsv")))

    def test_a_cell_dropped_from_a_year_stops_the_write(self):
        self.assertEqual(self.go(*world())[0], 0)
        reading = self.loaded()
        reading["cells"] = [x for x in reading["cells"] if x["year"] != "year_unknown"]
        with open(self.reading, "w", encoding="utf-8") as f:
            json.dump(reading, f)
        rc, text = self.write()
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 4: reader_call genomes sums to 2 over the years, its total is 3", text)

    def test_strata_that_do_not_follow_from_step1_stop_the_write(self):
        self.assertEqual(self.go(*world())[0], 0)
        reading = self.loaded()
        reading["mode"] = "four"
        with open(self.reading, "w", encoding="utf-8") as f:
            json.dump(reading, f)
        rc, text = self.write()
        self.assertEqual(rc, 2)
        self.assertIn("stop: write: the reading's strata do not follow from step 1", text)

    def test_a_year_under_100_genomes_is_thin(self):
        self.assertEqual(self.go(*thick())[0], 0)
        self.assertEqual(self.write()[0], 0)
        thin = {(r[0], r[1]): (r[-2], r[-1]) for r in self.table("shares.tsv")[1:]}
        self.assertEqual(thin[("2018", "reader_call")], ("100", "no"))
        self.assertEqual(thin[("2019", "reader_call")], ("99", "yes"))
        self.assertEqual({r[0]: r[-1] for r in self.table("years.tsv")[1:]}, {"2018": "no", "2019": "yes"})

    def test_recall_only_over_ten_resistant(self):
        self.assertEqual(self.go(*thick())[0], 0)
        self.assertEqual(self.write()[0], 0)
        years = {(r[0], r[1]): r for r in self.table("years.tsv")[1:]}
        self.assertEqual(years[("2018", "reader_call")][2:8], ["10", "10", "0", "7", "0", "0.700000000000"])
        self.assertEqual(years[("2019", "reader_no_call")][2:8], ["99", "9", "90", "9", "0", ""])

    def test_the_species_table_keeps_strata_with_ten_or_more_resistant(self):
        self.assertEqual(self.go(*thick())[0], 0)
        self.assertEqual(self.write()[0], 0)
        self.assertEqual(self.table("species.tsv")[1:], [
            [SPECIES, "drug t", "reader_call", "10", "10", "0", "7", "0", "0.700000000000", "0"]])

    def test_shares_of_a_year_add_to_one_with_their_denominators(self):
        self.assertEqual(self.go(*thick())[0], 0)
        self.assertEqual(self.write()[0], 0)
        rows = self.table("shares.tsv")[1:]
        for year in ("2018", "2019"):
            mine = [r for r in rows if r[0] == year]
            self.assertEqual(sum(int(r[2]) for r in mine), int(mine[0][3]))
            self.assertEqual(sum(int(r[5]) for r in mine), int(mine[0][6]))
            self.assertAlmostEqual(sum(float(r[4]) for r in mine), 1.0, places=9)
            self.assertAlmostEqual(sum(float(r[7]) for r in mine), 1.0, places=9)

    def test_no_id_or_group_name_in_any_written_file(self):
        rc, text = self.go(*world())
        self.assertIn("guard 6 (reading): 1 files, 0 with a genome id, 0 with a catalogue group name, against 11 ids"
                      " and 3 group names", text)
        rc, text = self.write()
        self.assertEqual(rc, 0, text)
        self.assertIn("guard 6 (tables): 6 files, 0 with a genome id, 0 with a catalogue group name", text)
        for name in os.listdir(self.out):
            with open(os.path.join(self.out, name), encoding="utf-8") as f:
                body = f.read()
            for word in IDS + GROUPS + ["samples"]:
                self.assertNotIn(word, body, name)
        with open(self.reading, encoding="utf-8") as f:
            body = f.read()
        for word in IDS + GROUPS + ["\"samples\"", "\"id\""]:
            self.assertNotIn(word, body)

    def plant(self, word):
        """The thick dataset's antibiotic replaced by `word` in the reading:
        species.tsv writes it in its one row."""
        self.assertEqual(self.go(*thick())[0], 0)
        reading = self.loaded()
        for x in reading["cells"]:
            x["antibiotic"] = word
        with open(self.reading, "w", encoding="utf-8") as f:
            json.dump(reading, f)
        return self.write()

    def test_an_id_in_a_written_file_stops_the_write(self):
        rc, text = self.plant("9999.3001")
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 6 (tables): 6 files, 1 with a genome id, 0 with a catalogue group name", text)
        self.assertNotIn("9999.3001", text)
        self.assertEqual(os.listdir(self.out), [])

    def test_a_group_name_in_a_written_file_stops_the_write(self):
        rc, text = self.plant(GROUPS[2])
        self.assertEqual(rc, 2)
        self.assertIn("stop: guard 6 (tables): 6 files, 0 with a genome id, 1 with a catalogue group name", text)
        self.assertNotIn(GROUPS[2], text)
        self.assertEqual(os.listdir(self.out), [])

    def test_the_scan_reads_whole_names_not_parts_of_words(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.tsv")
            with open(path, "w", encoding="utf-8") as f:
                f.write("catalogue\t0.99991001\tsynBx\t19999.1001\n")
            self.assertEqual(ym.scan_hits([path], {"9999.1001"}, {"cat", "synB"}), (0, 0))
            with open(path, "w", encoding="utf-8") as f:
                f.write("a cat\t9999.1001\n")
            self.assertEqual(ym.scan_hits([path], {"9999.1001"}, {"cat"}), (1, 1))

    def test_the_messages_print_no_id_or_group_name(self):
        rc, text = self.go(*world())
        rc2, text2 = self.write()
        self.assertEqual((rc, rc2), (0, 0))
        for word in IDS + GROUPS:
            self.assertNotIn(word, text + text2)

    def test_the_same_reports_give_the_same_bytes(self):
        self.assertEqual(self.go(*world())[0], 0)
        self.assertEqual(self.write()[0], 0)
        again = os.path.join(self.tmp.name, "again")
        self.assertEqual(run("write", "--reading", self.reading, "--reports", self.rdir, "--out", again)[0], 0)
        for name in sorted(os.listdir(self.out)):
            with open(os.path.join(self.out, name), "rb") as a, open(os.path.join(again, name), "rb") as b:
                self.assertEqual(a.read(), b.read(), name)
        self.assertEqual(sorted(os.listdir(self.out)), sorted(list(ym.TABLES) + ["SHA256SUMS"]))

    def test_write_refuses_a_folder_that_is_not_empty(self):
        self.assertEqual(self.go(*world())[0], 0)
        os.makedirs(self.out)
        with open(os.path.join(self.out, "x"), "w", encoding="utf-8") as f:
            f.write("x\n")
        rc, text = self.write()
        self.assertEqual(rc, 2)
        self.assertIn("is not empty", text)


if __name__ == "__main__":
    unittest.main()
