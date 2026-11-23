"""scripts/mechanism/history.py: block M2, the mechanism mix over collection years (docs/mechanism-split.md).

Every report, table, genome id, BioProject, year and count is SYNTHETIC.
"""

import contextlib
import hashlib
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.mechanism import common as c

history = c.load("history")
amr = c.load("amr_rule")
scan = c.load("scan")
OTHER = "synthia_coli__drug_b"


def project(prefix, n, year, kind, bp, country="Aland", lineage=0, kinds=None):
    """n genomes of one SYNTHETIC BioProject; odd members resistant."""
    out = []
    for i in range(1, n + 1):
        out.append({"id": f"{prefix}.{i}", "kind": (kinds[(i - 1) % len(kinds)] if kinds else kind),
                    "resistant": i % 2 == 1, "year": year, "country": country, "lineage": lineage, "bp": bp})
    return out


def scenario():
    """Earlier (2010 to 2014): three BioProjects of gene rows. Later (2015
    to 2019): three of point and gene rows. Beside them: years that are no
    year, thin years, and genomes with no BioProject."""
    g = (project("8101", 7, "2014", "gene", "PRJNA900001", lineage=1)
         + project("8102", 5, "2012", None, "PRJNA900002", lineage=2, kinds=["gene", "gene_b"])
         + project("8103", 5, "2010", "gene", "PRJNA900003", lineage=3)
         + project("8104", 5, "2015", "both", "PRJNA900004", lineage=4)
         + project("8105", 5, "2017", "both", "PRJNA900005", lineage=5)
         + project("8106", 5, "2019", "both", "PRJNA900006", country="Borduria", lineage=6)
         + project("8107", 5, "2016-03", "gene", "PRJNA900007", lineage=7)
         + project("8108", 6, None, "gene", None, lineage=8)
         + project("8109", 3, "2020", "both", "PRJNA900004", lineage=4)
         + project("8110", 2, "1998", "gene", "PRJNA900001", lineage=1))
    return g


def layout(top, genomes, other=None, folds=None):
    top = Path(top)
    c.write_species(top, {c.DATASET: (genomes, {}), OTHER: (other or genomes, {})})
    c.write_db(top / "db")
    sha = c.write_meta(top / "meta.jsonl.gz", {g["id"]: g["bp"] for g in genomes})
    return sha


def run(top, sha, folds=None, design_sha=None):
    top = Path(top)
    argv = ["--species", c.SPECIES, "--reports", str(top / "reports"), "--lists", str(top / "lists.tsv"),
            "--drug-rules", str(top / "drug_rules.tsv"), "--tables", str(top / "tool" / c.SPECIES),
            "--db", str(top / "db"), "--meta", str(top / "meta.jsonl.gz"), "--out", str(top / "out")]
    if folds:
        argv += ["--folds", str(folds), "--design", str(folds.parent / "design.sha256")]
    with mock.patch.object(history, "META_SHA256", sha), \
            mock.patch.object(history, "DESIGN_SHA256", design_sha or history.DESIGN_SHA256), \
            contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
        rc = history.main(argv)
    return rc, err.getvalue(), top / "out"


def row(out, guard, periods="primary", cell="species_class", cls="QUINOLONE", dataset=""):
    for r in c.read_tsv(out / "m2-cells.tsv"):
        if (r["guard"], r["periods"], r["cell"], r["class"], r["dataset"]) == (guard, periods, cell, cls, dataset):
            return r
    raise AssertionError(f"no row {guard} {periods} {cell} {cls} {dataset}")


class Years(unittest.TestCase):
    def test_a_year_is_four_digits_and_nothing_else(self):
        self.assertEqual(history.year_of("2015"), 2015)
        self.assertEqual(history.year_of(" 2015 "), 2015)
        for text in (None, "", "2015-06", "c. 2014", "15", "20150", "2014/2015", 2015):
            self.assertIsNone(history.year_of(text), text)

    def test_period_ends_are_inclusive_and_exact(self):
        p = history.period_of
        self.assertEqual([p(y, "primary") for y in (2009, 2010, 2014, 2015, 2019, 2020)],
                         [None, "earlier", "earlier", "later", "later", None])
        self.assertEqual([p(y, "secondary") for y in (1999, 2000, 2009, 2010, 2019, 2020)],
                         [None, "earlier", "earlier", "later", "later", None])
        self.assertIsNone(p(None, "primary"))


class Families(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.db = amr.Database(c.write_db(Path(self.tmp.name) / "db"))

    def test_family_is_the_databases_own_name(self):
        self.assertEqual(self.db.family("qnrS1"), "qnrS")
        self.assertEqual(self.db.family("blaKPC-2"), "blaKPC")
        self.assertEqual(self.db.family("tetS"), "tetS")
        self.assertEqual(self.db.family("blaX-9"), "blaX")  # through the hierarchy
        with self.assertRaises(amr.Refusal):
            self.db.family("madeUp-1")

    def test_a_disrupted_core_gene_enters_no_family_mix(self):
        # the lead's rulings on questions 1 and 10: a POINT_DISRUPT row is read with point
        for subtype in ("POINT_DISRUPT", "POINT"):
            with self.assertRaises(amr.Refusal) as cm:
                self.db.family("synR_Q11Ter42", subtype)
            self.assertNotIn("Q11", str(cm.exception))
        with self.assertRaises(amr.Refusal):
            self.db.family("synR_Q11Ter42", "AMR")
        rows = [{"symbol": "synR_Q11Ter42", "scope": "core", "type": "AMR", "subtype": "POINT_DISRUPT",
                 "class": "QUINOLONE", "subclass": "QUINOLONE", "method": "BLASTX"},
                {"symbol": "qnrS1", "scope": "core", "type": "AMR", "subtype": "AMR",
                 "class": "QUINOLONE", "subclass": "QUINOLONE", "method": "EXACTX"}]
        self.assertEqual(history.genome_cells(rows, self.db)["QUINOLONE"], (True, True, frozenset({"qnrS"}), True))
        self.assertEqual(history.genome_cells(rows[:1], self.db)["QUINOLONE"], (True, False, frozenset(), True))

    def test_a_rule_names_a_class_by_c_or_by_a_subclass_of_it(self):
        self.assertTrue(self.db.rule_names_class([["C:QUINOLONE"]], "QUINOLONE"))
        self.assertTrue(self.db.rule_names_class([["S:CARBAPENEM"]], "BETA-LACTAM"))
        self.assertFalse(self.db.rule_names_class([["S:CARBAPENEM"]], "QUINOLONE"))


class Vote(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)

    def go(self, genomes=None, other=None):
        genomes = genomes or scenario()
        sha = layout(self.top, genomes, other)
        return run(self.top, sha)

    def test_one_vote_per_bioproject(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        r = row(out, "bioproject")
        self.assertEqual((r["passes"], r["earlier_genomes"], r["later_genomes"]), ("1", "17", "15"))
        self.assertEqual((r["earlier_point"], r["later_point"], r["d_point"]),
                         ("0.000000000000", "1.000000000000", "1.000000000000"))
        # mean of the three earlier BioProjects' mixes: qnrS 1, 0.6, 1 and qnrB 0, 0.4, 0
        self.assertEqual(r["tvd"], "0.133333333333")

    def test_a_period_boundary_is_exact(self):
        rc, err, out = self.go()
        r = row(out, "bioproject")
        # 2014 is earlier and 2015 later; 2020 and 1998 enter no period
        self.assertEqual((r["earlier_genomes"], r["later_genomes"]), ("17", "15"))
        s = row(out, "bioproject", periods="secondary")
        self.assertEqual((s["earlier_genomes"], s["later_genomes"], s["passes"]), ("0", "32", "0"))

    def test_no_year_is_imputed(self):
        rc, err, out = self.go()
        years = {r["years"]: r["genomes"] for r in c.read_tsv(out / "m2-years.tsv")}
        self.assertEqual(years["no year"], "11")
        self.assertEqual(years["2020 on"], "under 5")
        self.assertEqual(years["before 2000"], "under 5")
        self.assertEqual(years["2010 to 2014"], "17")
        self.assertEqual(row(out, "bioproject")["with_year_share"], f"{37 / 48:.12f}")

    def test_permutation_moves_bioprojects_not_genomes(self):
        rc, err, out = self.go()
        p = float(row(out, "bioproject")["p_point"])
        # 2 of the 20 ways to label 3 of 6 BioProjects later give |d| = 1
        self.assertGreater(p, 0.05)
        self.assertLess(p, 0.2)

    def test_permutation_is_seeded_by_the_cell_name(self):
        rc, err, out = self.go()
        first = (out / "m2-cells.tsv").read_bytes()
        rc, err, out = run(self.top, hashlib.sha256((self.top / "meta.jsonl.gz").read_bytes()).hexdigest())
        self.assertEqual((out / "m2-cells.tsv").read_bytes(), first)

    def test_family_mix_names_the_databases_families(self):
        rc, err, out = self.go()
        fams = {r["family"]: r for r in c.read_tsv(out / "m2-families.tsv") if r["periods"] == "primary"}
        self.assertEqual(set(fams), {"qnrB", "qnrS"})
        self.assertEqual(fams["qnrS"]["earlier_mix"], f"{(1 + 0.6 + 1) / 3:.12f}")
        self.assertEqual(fams["qnrS"]["later_mix"], "1.000000000000")

    def test_largest_bioproject_left_out(self):
        rc, err, out = self.go()
        r = row(out, "largest_out")
        self.assertEqual((r["passes"], r["earlier_genomes"]), ("0", "10"))

    def test_countries_in_both_periods_only(self):
        rc, err, out = self.go()
        r = row(out, "countries_both")
        self.assertEqual((r["passes"], r["later_genomes"]), ("0", "10"))

    def test_lineage_clusters_vote_when_reports_agree(self):
        rc, err, out = self.go()
        self.assertEqual(row(out, "lineage")["passes"], "1")
        lin = c.read_tsv(out / "m2-lineage.tsv")[0]
        self.assertEqual((lin["given"], lin["differing_genome_pairs"]), ("1", "0"))

    def test_lineage_guard_not_given_when_reports_differ(self):
        g = scenario()
        other = [dict(x) for x in g]
        other[0]["lineage"] = 2  # one genome of the first BioProject sits with the second in one report
        rc, err, out = self.go(g, other)
        self.assertEqual(rc, 0, err)
        lin = c.read_tsv(out / "m2-lineage.tsv")[0]
        # in one report it shares a cluster with 8 genomes, in the other with 5
        self.assertEqual((lin["given"], lin["differing_genome_pairs"]), ("0", str(8 + 5)))
        self.assertIn("not given", row(out, "lineage")["note"])

    def test_a_year_that_differs_between_reports_stops_the_species(self):
        g = scenario()
        other = [dict(x) for x in g]
        other[3]["year"] = "2011"
        rc, err, out = self.go(g, other)
        self.assertEqual(rc, 1)
        self.assertIn("year differs", err)

    def test_the_bvbrc_table_is_held_to_its_sha256(self):
        layout(self.top, scenario())
        rc, err, out = run(self.top, "0" * 64)
        self.assertEqual(rc, 1)
        self.assertIn("sha256", err)

    def test_resistant_genomes_per_dataset_are_secondary(self):
        rc, err, out = self.go()
        r = row(out, "bioproject", cell="dataset_resistant", cls="", dataset=c.DATASET)
        self.assertEqual(r["genomes"], "28")

    def test_no_genome_id_or_bioproject_in_any_output(self):
        rc, err, out = self.go()
        ids = self.top / "species.json"
        ids.write_text(json.dumps({"genomes": [{"id": g["id"]} for g in scenario()]}))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(scan.main(["--ids", str(ids), str(out)]), 0)


class Guards(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.db = amr.Database(c.write_db(Path(self.tmp.name) / "db"))

    def test_flagged_by_resistance_or_by_one_family(self):
        bp = {f"1.{i}": "ABC"[i // 5] for i in range(15)}
        genomes = {g: {"bioproject": b, "year": 2012, "country": None} for g, b in bp.items()}
        fam = {g: (frozenset({"qnrS"}) if b == "B" else frozenset({["qnrS", "qnrB"][int(g[2:]) % 2]}))
               for g, b in bp.items()}
        cells = {g: {"QUINOLONE": (False, True, fam[g])} for g in bp}
        labels = {g: (1 if b == "A" else int(g[2:]) % 2) for g, b in bp.items()}
        # A: every labelled pair resistant. B: every genome with a gene row carries qnrS. C: neither.
        got = history.flagged_bioprojects("QUINOLONE", sorted(bp), genomes, cells, [([["C:QUINOLONE"]], labels)],
                                          self.db)
        self.assertEqual(got, {"A", "B"})
        # pairs count only from datasets whose rule names the class
        got = history.flagged_bioprojects("QUINOLONE", sorted(bp), genomes, cells, [([["C:TETRACYCLINE"]], labels)],
                                          self.db)
        self.assertEqual(got, {"B"})

    def test_intact_rows_drop_partial_calls(self):
        tmp = Path(self.tmp.name)
        c.write_table(tmp / "t.tsv", "partial")
        rows = amr.read_table(tmp / "t.tsv")
        self.assertEqual(history.genome_cells(rows, self.db)["QUINOLONE"], (False, True, frozenset({"qnrB"}), False))
        self.assertEqual(history.genome_cells(rows, self.db, "tool_default_intact"), {})

    def test_plus_scope_and_another_class_are_not_counted(self):
        tmp = Path(self.tmp.name)
        c.write_table(tmp / "t.tsv", "plus_only")
        self.assertEqual(history.genome_cells(amr.read_table(tmp / "t.tsv"), self.db), {})
        c.write_table(tmp / "u.tsv", "neither")
        self.assertEqual(set(history.genome_cells(amr.read_table(tmp / "u.tsv"), self.db)), {"TETRACYCLINE"})


class TimeSplit(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)

    def go(self, train_kind="gene", test_kind="both"):
        g = scenario()
        for x in g:
            if x["year"] in ("2010", "2012", "2014"):
                x["kind"] = train_kind if x["kind"] == "gene" else x["kind"]
            if x["year"] in ("2015", "2017", "2019"):
                x["kind"] = test_kind if x["kind"] == "both" else x["kind"]
        sha = layout(self.top, g)
        folds = self.top / "ts" / "fold-files"
        folds.mkdir(parents=True)
        train = [x["id"] for x in g if x["year"] in ("2010", "2012", "2014")]
        test = [x["id"] for x in g if x["year"] in ("2015", "2017", "2019")]
        doc = {"format": "grm-fold-file-v1", "dataset": c.DATASET, "scheme": "first-public-past", "train": train,
               "test": [test], "unused": []}
        (folds / f"{c.DATASET}.past.json").write_text(json.dumps(doc))
        h = hashlib.sha256((folds / f"{c.DATASET}.past.json").read_bytes()).hexdigest()
        design = self.top / "ts" / "design.sha256"
        design.write_text(f"# made up\n{h}  results/fold-files/{c.DATASET}.past.json\n")
        return run(self.top, sha, folds, hashlib.sha256(design.read_bytes()).hexdigest())

    def test_unseen_determinant_share(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        rows = {r["reading"]: r for r in c.read_tsv(out / "timesplit.tsv")}
        r = rows["point_as_one"]
        self.assertEqual((r["resistant_later_test"], r["with_determinant"], r["with_unseen_determinant"]), ("9", "9", "9"))
        self.assertEqual(r["share_unseen"], "1.000000000000")

    def test_a_point_row_seen_in_training_is_not_unseen(self):
        rc, err, out = self.go(train_kind="both")
        self.assertEqual(rc, 0, err)
        r = {r["reading"]: r for r in c.read_tsv(out / "timesplit.tsv")}["point_as_one"]
        self.assertEqual((r["with_unseen_determinant"], r["share_unseen"]), ("0", "0.000000000000"))

    def test_a_point_disrupt_row_seen_in_training_is_read_as_point(self):
        # the lead's ruling on question 1: a POINT_DISRUPT row is point in the time-split link too
        rc, err, out = self.go(train_kind="disrupt", test_kind="point")
        self.assertEqual(rc, 0, err)
        rows = {r["reading"]: r for r in c.read_tsv(out / "timesplit.tsv")}
        self.assertEqual((rows["point_as_one"]["with_unseen_determinant"], rows["point_as_one"]["share_unseen"]),
                         ("0", "0.000000000000"))
        self.assertEqual(rows["point_by_symbol"]["with_unseen_determinant"], "9")

    def test_a_fold_file_not_in_the_design_stops(self):
        g = scenario()
        sha = layout(self.top, g)
        folds = self.top / "ts" / "fold-files"
        folds.mkdir(parents=True)
        (folds / f"{c.DATASET}.past.json").write_text("{}")
        design = self.top / "ts" / "design.sha256"
        design.write_text("# made up\n")
        rc, err, out = run(self.top, sha, folds, hashlib.sha256(design.read_bytes()).hexdigest())
        self.assertEqual(rc, 1)
        self.assertIn("not in design.sha256", err)


if __name__ == "__main__":
    unittest.main()
