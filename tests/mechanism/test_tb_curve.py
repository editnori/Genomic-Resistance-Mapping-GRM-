"""scripts/mechanism/tb_curve.py, block M4 (docs/mechanism-split.md).

Every genome id, result file, gene, change and lineage here is SYNTHETIC.
"""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

tb = c.load("tb_curve")
TB = "synthia_tuberculosis"
DRUGS = ("rifampicin", "isoniazid", "levofloxacin", "moxifloxacin", "ethionamide")


def variant(gene, change, drug, grade="Assoc w R", filt="pass", kind="drug_resistance"):
    return {"gene_name": gene, "change": change, "filter": filt, "freq": 1.0, "type": "missense_variant",
            "ref": "A", "alt": "G", "sv": False,
            "drugs": [{"type": kind, "drug": drug, "confidence": grade}]}


def result(gid, main, sub, variants, tool="6.6.3", db="3ed1cb99"):
    return {"id": gid, "pipeline": {"software_version": tool, "db_version": {"commit": db}},
            "main_lineage": main, "sub_lineage": sub, "dr_variants": variants,
            "qc": {"target_qc": [], "num_sequences": 1, "num_bases": 1}}


def write_run(top, genomes, failed=()):
    """genomes: list of (id, main, sub, variants)."""
    run = Path(top) / "tool"
    (run / "results").mkdir(parents=True)
    (run / "plan.txt").write_text("".join(g[0] + "\n" for g in genomes) + "".join(f + "\n" for f in failed))
    if failed:
        (run / "failed.tsv").write_text("".join(f"{f}\ttimeout\n" for f in failed))
    for gid, main, sub, vs in genomes:
        (run / "results" / f"{gid}.results.json").write_text(json.dumps(result(gid, main, sub, vs)))
    bed = Path(top) / "tbdb.bed"
    bed.write_text("".join(f"Chromosome\t1\t2\tsyn{i}\tsynG{i}\t{d}\n" for i, d in enumerate(DRUGS)))
    return run, bed


def write_reports(top, antibiotics, samples):
    """samples: list of (id, resistant, cluster)."""
    rdir = Path(top) / "reports"
    rdir.mkdir()
    lists = []
    for ab, in_m in antibiotics:
        ds = f"{TB}__{ab}"
        rep = {"id": ds, "species_slug": TB, "species": "Synthia tuberculosis", "antibiotic": ab,
               "samples": [{"id": gid, "resistant": res, "lineage": cl} for gid, res, cl in samples]}
        (rdir / f"{ds}.json").write_text(json.dumps(rep))
        lists.append({"dataset": ds, "species_slug": TB, "species": "Synthia tuberculosis", "antibiotic": ab,
                      "in_M": "1" if in_m else "0", "in_A": "0", "labels_depend_on": ""})
    lp = Path(top) / "lists.tsv"
    c.write_tsv(lp, ["dataset", "species_slug", "species", "antibiotic", "in_M", "in_A", "labels_depend_on"], lists)
    return rdir, lp


def run_tb(top, run, bed, rdir, lp):
    out = Path(top) / "out"
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
        rc = tb.main(["--species", TB, "--reports", str(rdir), "--lists", str(lp), "--tool-dir", str(run),
                      "--bed", str(bed), "--out", str(out)])
    return rc, out, err.getvalue()


def lineage_set(n_per, mains, dets_of):
    """n_per genomes in each main lineage; dets_of(lineage index, i) -> variants."""
    genomes, samples = [], []
    for li, main in enumerate(mains):
        for i in range(n_per):
            gid = f"9{li}{i:03d}.{li + 1}"
            genomes.append((gid, main, f"{main}.{1 + i % 2}", dets_of(li, i)))
            samples.append((gid, True, li))
    return genomes, samples


class DrugMap(unittest.TestCase):
    def test_names_synonyms_and_proxies(self):
        self.assertEqual(tb.drug_map("rifampin", DRUGS), ("primary", ("rifampicin",)))
        self.assertEqual(tb.drug_map("Isoniazid", DRUGS), ("primary", ("isoniazid",)))
        self.assertEqual(tb.drug_map("ciprofloxacin", DRUGS), ("proxy", ("levofloxacin", "moxifloxacin")))
        self.assertEqual(tb.drug_map("prothionamide", DRUGS), ("proxy", ("ethionamide",)))
        self.assertEqual(tb.drug_map("amoxicillin", DRUGS), ("none", ()))

    def test_a_proxy_the_tool_lacks_refuses(self):
        with self.assertRaises(tb.Refusal):
            tb.drug_map("rifabutin", ("isoniazid",))


class Determinants(unittest.TestCase):
    def test_only_passing_counted_resistance_grades_count(self):
        vs = [variant("synA", "p.Ser1Leu", "rifampicin"),
              variant("synA", "p.Asp2Val", "rifampicin", grade="Assoc w R - Interim"),
              variant("synA", "p.His3Tyr", "rifampicin", grade="Uncertain significance"),
              variant("synA", "p.Leu4Pro", "rifampicin", filt="hard_fail"),
              variant("synB", "p.Ser5Thr", "rifampicin", kind="other_annotation")]
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.json"
            p.write_text(json.dumps(result("1.1", "lineage4", "lineage4.1", vs)))
            main, sub, per = tb.genome_determinants(str(p), "1.1")
        self.assertEqual((main, sub), ("lineage4", "lineage4.1"))
        self.assertEqual(per, {"rifampicin": {("synA", "p.Ser1Leu"), ("synA", "p.Asp2Val")}})

    def test_another_tool_database_or_genome_refuses(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.json"
            for kw, gid in (({"tool": "6.6.2"}, "1.1"), ({"db": "0000000"}, "1.1"), ({}, "2.2")):
                p.write_text(json.dumps(result("1.1", "lineage4", "lineage4.1", [], **kw)))
                with self.assertRaises(tb.Refusal) as cm:
                    tb.genome_determinants(str(p), gid)
                self.assertNotIn("1.1", str(cm.exception))


class Curve(unittest.TestCase):
    def groups(self):
        return {"L1": {"dets": {("g", "a"), ("g", "b")}, "genomes": 30, "without": 2},
                "L2": {"dets": {("g", "b"), ("h", "c")}, "genomes": 20, "without": 0},
                "L3": {"dets": {("h", "d")}, "genomes": 6, "without": 6}}

    def test_largest_first_and_the_random_band(self):
        rows = tb.curve_rows("S", "rifampin", False, "main_lineage", self.groups())
        self.assertEqual([r["determinants_largest_first"] for r in rows], [2, 3, 4])
        self.assertEqual([r["step"] for r in rows], [1, 2, 3])
        last = rows[-1]
        self.assertEqual((last["determinants_median"], last["determinants_lo"], last["determinants_hi"]), (4, 4, 4))
        first = rows[0]
        self.assertEqual((first["determinants_lo"], first["determinants_hi"]), (1, 2))
        for r in rows:
            self.assertLessEqual(r["determinants_lo"], r["determinants_median"])
            self.assertLessEqual(r["determinants_median"], r["determinants_hi"])

    def test_orders_are_seeded_by_antibiotic_and_grouping(self):
        a = tb.curve_rows("S", "rifampin", False, "main_lineage", self.groups())
        b = tb.curve_rows("S", "rifampin", False, "main_lineage", self.groups())
        self.assertEqual(a, b)

    def test_convergence_buckets_and_cover(self):
        per = []
        for li in range(5):
            dets = {("g", "everywhere")}
            if li < 2:
                dets.add(("g", "two"))
            if li == 0:
                dets.add(("h", "one"))
            per += [(f"L{li}", dets)] * 2
        per.append(("L9", set()))
        rows = {r["lineages_seen"]: r for r in tb.convergence_rows("S", "rifampin", False, "main_lineage", per)}
        self.assertEqual({k: rows[k]["determinants"] for k in tb.BUCKETS}, {"1": 1, "2": 1, "3": 0, "4+": 1})
        self.assertEqual(rows["4+"]["resistant_covered"], "10")
        self.assertEqual(rows["2"]["resistant_covered"], "under 5")
        self.assertEqual(rows["4+"]["share_covered"], f"{10 / 11:.12f}")


class Genes(unittest.TestCase):
    def test_a_gene_seen_in_under_5_genomes_is_not_named(self):
        dets = {"synR": {("synR", "a"), ("synR", "b")}, "synK": {("synK", "c")}, "synG": {("synG", "d")}}
        sets = [{"synR"}] * 6 + [{"synR", "synK"}] * 2 + [{"synG"}]
        rows = {r["gene"]: r for r in tb.gene_rows("Synthia s0", "rifampin", False, dets, sets)}
        self.assertEqual(set(rows), {"synR", tb.FOLDED_GENES})
        self.assertEqual((rows["synR"]["catalogue_variants"], rows["synR"]["resistant_with_gene"]), (2, "8"))
        folded = rows[tb.FOLDED_GENES]
        self.assertEqual((folded["catalogue_variants"], folded["resistant_with_gene"]), (2, "under 5"))

    def test_a_gene_in_5_genomes_is_named(self):
        rows = tb.gene_rows("Synthia s0", "rifampin", False, {"synK": {("synK", "c")}}, [{"synK"}] * 5)
        self.assertEqual([(r["gene"], r["resistant_with_gene"]) for r in rows], [("synK", "5")])


class SmallCurve(unittest.TestCase):
    def test_no_curve_over_under_5_resistant_genomes(self):
        groups = {"lineage1": {"dets": {("synR", "a")}, "genomes": 3, "without": 0},
                  "lineage2": {"dets": set(), "genomes": 1, "without": 1}}
        self.assertEqual(tb.curve_rows("Synthia s0", "rifampin", False, "main_lineage", groups), [])
        groups["lineage2"]["genomes"] = 2
        rows = tb.curve_rows("Synthia s0", "rifampin", False, "main_lineage", groups)
        self.assertEqual([r["resistant"] for r in rows], ["5", "5"])


class EndToEnd(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def build(self, failed=(), in_m=False):
        def dets(li, i):
            vs = []
            if i % 3:
                vs.append(variant("synR", f"p.Ser{450 + li}Leu", "rifampicin"))
            if li == 0 and i % 2:
                vs.append(variant("synK", "p.Ser315Thr", "isoniazid"))
            vs.append(variant("synG", "p.Asp94Gly", "levofloxacin", grade="Uncertain significance"))
            return vs
        genomes, samples = lineage_set(8, ("lineage1", "lineage2", "lineage4"), dets)
        samples += [("8000.1", False, 7), ("8001.1", False, 7)]
        genomes += [("8000.1", "lineage4", "lineage4.9", []), ("8001.1", "lineage4", "lineage4.9", [])]
        run, bed = write_run(self.top, genomes, failed)
        rdir, lp = write_reports(self.top, [("rifampin", in_m), ("amoxicillin", False), ("ciprofloxacin", False)],
                                 samples + [(f, True, 3) for f in failed])
        return genomes, run, bed, rdir, lp

    def test_runs_and_names_no_variant_or_genome(self):
        genomes, run, bed, rdir, lp = self.build()
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        names = sorted(p.name for p in out.iterdir())
        self.assertEqual(names, ["m4-convergence.tsv", "m4-coverage.tsv", "m4-curve.tsv", "m4-genes.tsv",
                                 "m4-summary.json"])
        text = "".join(p.read_text() for p in out.iterdir())
        for gid, *_ in genomes:
            self.assertNotIn(gid, text)
        self.assertNotIn("p.", text)
        for word in ("Ser", "Asp94", "450", "315"):
            self.assertNotIn(word, text)
        summary = json.loads((out / "m4-summary.json").read_text())
        self.assertEqual([e["antibiotic"] for e in summary["endpoints"]], ["ciprofloxacin", "rifampin"])
        self.assertEqual(summary["not_scored_reports"][0]["antibiotic"], "amoxicillin")

    def test_gene_rows_count_catalogue_variants_and_genomes(self):
        *_, run, bed, rdir, lp = self.build()
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        genes = c.read_tsv(out / "m4-genes.tsv")
        rif = {r["gene"]: r for r in genes if r["antibiotic"] == "rifampin"}
        self.assertEqual(set(rif), {"synR"})
        self.assertEqual(rif["synR"]["catalogue_variants"], "3")
        self.assertEqual(rif["synR"]["resistant_with_gene"], "15")
        self.assertFalse([r for r in genes if r["antibiotic"] == "ciprofloxacin"])

    def test_no_determinant_share_per_group_and_overall(self):
        *_, run, bed, rdir, lp = self.build()
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        cov = [r for r in c.read_tsv(out / "m4-coverage.tsv") if r["antibiotic"] == "rifampin"]
        overall = [r for r in cov if r["grouping"] == "all"]
        self.assertEqual(len(overall), 1)
        self.assertEqual((overall[0]["resistant"], overall[0]["without_determinant"]), ("24", "9"))
        self.assertEqual(overall[0]["share_without"], f"{9 / 24:.12f}")
        main = {r["group"]: r for r in cov if r["grouping"] == "main_lineage"}
        self.assertEqual(set(main), {"lineage1", "lineage2", "lineage4"})
        self.assertEqual(main["lineage1"]["without_determinant"], "under 5")
        cipro = [r for r in c.read_tsv(out / "m4-coverage.tsv")
                 if r["antibiotic"] == "ciprofloxacin" and r["grouping"] == "all"]
        self.assertEqual(cipro[0]["without_determinant"], "24")

    def test_curve_reaches_every_distinct_determinant(self):
        *_, run, bed, rdir, lp = self.build()
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        curve = [r for r in c.read_tsv(out / "m4-curve.tsv")
                 if r["antibiotic"] == "rifampin" and r["grouping"] == "main_lineage"]
        self.assertEqual([r["step"] for r in curve], ["1", "2", "3"])
        self.assertEqual(curve[-1]["determinants_largest_first"], "3")
        self.assertEqual(curve[-1]["determinants_lo"], "3")
        self.assertEqual({r["grouping"] for r in c.read_tsv(out / "m4-curve.tsv")},
                         {"main_lineage", "sub_lineage", "cluster"})

    def test_failed_genomes_are_counted_not_read(self):
        *_, run, bed, rdir, lp = self.build(failed=("7000.1", "7001.1"))
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        summary = json.loads((out / "m4-summary.json").read_text())
        self.assertEqual(summary["genomes_failed"], "under 5")
        rif = [e for e in summary["endpoints"] if e["antibiotic"] == "rifampin"][0]
        self.assertEqual((rif["resistant_with_result"], rif["resistant_no_result"]), ("24", "under 5"))

    def test_a_genome_neither_run_nor_failed_refuses(self):
        *_, run, bed, rdir, lp = self.build()
        (run / "plan.txt").write_text((run / "plan.txt").read_text() + "7002.1\n")
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 1)
        self.assertIn("refused", err)
        self.assertNotIn("7002.1", err)

    def test_a_dataset_in_m_is_marked(self):
        *_, run, bed, rdir, lp = self.build(in_m=True)
        rc, out, err = run_tb(self.top, run, bed, rdir, lp)
        self.assertEqual(rc, 0, err)
        for name in ("m4-curve.tsv", "m4-coverage.tsv", "m4-convergence.tsv", "m4-genes.tsv"):
            marks = {(r["antibiotic"], r["in_M"]) for r in c.read_tsv(out / name)}
            self.assertIn(("rifampin", "1"), marks, name)
            self.assertNotIn(("ciprofloxacin", "1"), marks, name)


if __name__ == "__main__":
    unittest.main()
