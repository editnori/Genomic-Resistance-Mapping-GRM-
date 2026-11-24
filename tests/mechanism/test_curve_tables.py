"""scripts/mechanism/curve_tables.py, the table script of block M4 (docs/lineage-curve-result.md).

Every genome id, result file, report, gene, change, lineage and count here is
SYNTHETIC: tb_curve.py runs on made-up TBProfiler result files and the table
script reads what it wrote.
"""

import contextlib
import io
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

tb = c.load("tb_curve")
ct = c.load("curve_tables")
scan = c.load("scan")
SLUG = "synthia_tuberculosis"
NAME = "Synthia tuberculosis"
DRUGS = ("rifampicin", "isoniazid", "levofloxacin", "moxifloxacin")
GENES = ("synR", "synK", "synG")
OTHER_DB = "0badc0de"


def variant(gene, change, drug, grade="Assoc w R"):
    return {"gene_name": gene, "change": change, "filter": "pass", "freq": 1.0, "type": "missense_variant",
            "ref": "A", "alt": "G", "sv": False, "drugs": [{"type": "drug_resistance", "drug": drug, "confidence": grade}]}


def result(gid, main, sub, variants, tool="6.6.3", db="3ed1cb99"):
    return {"id": gid, "pipeline": {"software_version": tool, "db_version": {"commit": db}},
            "main_lineage": main, "sub_lineage": sub, "dr_variants": variants}


def genomes():
    """(id, main, sub, variants, resistant, cluster): 3 main lineages of 8 resistant genomes, 2 susceptible."""
    out = []
    for li, main in enumerate(("lineage1", "lineage2", "lineage4")):
        for i in range(8):
            vs = []
            if i % 3:
                vs.append(variant("synR", f"p.Ser{450 + li}Leu", "rifampicin"))
            if li == 0 and i < 3:
                vs.append(variant("synK", "p.Ser315Thr", "rifampicin"))
            out.append((f"9{li}{i:03d}.{li + 1}", main, f"{main}.{1 + i % 2}", vs, True, li))
    out += [("8000.1", "lineage4", "lineage4.9", [], False, 7), ("8001.1", "lineage4", "lineage4.9", [], False, 7)]
    return out


class Reading:
    """A SYNTHETIC reading folder: the tool run, the reports, lists.tsv, tb_curve.py's output and checks.txt."""

    def __init__(self, top, in_m_run=False, in_m_lists=None, bad_version=False, bad_db=False, extra_planned=False):
        self.top = Path(top)
        run = self.top / "run"
        (run / "results").mkdir(parents=True)
        gs = genomes()
        (run / "plan.txt").write_text("".join(g[0] + "\n" for g in gs))
        (run / "failed.tsv").write_text("")
        for gid, main, sub, vs, _, _ in gs:
            (run / "results" / f"{gid}.results.json").write_text(json.dumps(result(gid, main, sub, vs)))
        (run / "run.start").write_text("2026-01-01T00:00:00Z\n")
        (run / "run.end").write_text("2026-01-01T00:30:00Z\n")
        self.bed = self.top / "tbdb.bed"
        self.bed.write_text("".join(f"Chromosome\t1\t2\tsyn{i}\t{GENES[i % 3]}\t{d}\n" for i, d in enumerate(DRUGS)))
        rdir = self.top / "reports"
        rdir.mkdir()
        self.lists = self.top / "lists.tsv"
        self.write_lists(in_m_run)
        for ab in ("rifampin", "amoxicillin", "ciprofloxacin"):
            rep = {"id": f"{SLUG}__{ab}", "species_slug": SLUG, "species": NAME, "antibiotic": ab,
                   "samples": [{"id": g[0], "resistant": g[4], "lineage": g[5]} for g in gs]}
            (rdir / f"{SLUG}__{ab}.json").write_text(json.dumps(rep))
        self.reading = self.top / "reading"
        self.reading.mkdir()
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
            rc = tb.main(["--species", SLUG, "--reports", str(rdir), "--lists", str(self.lists), "--tool-dir", str(run),
                          "--bed", str(self.bed), "--out", str(self.reading / "m4")])
        assert rc == 0, err.getvalue()
        if bad_version:
            gid = gs[0][0]
            (run / "results" / f"{gid}.results.json").write_text(json.dumps(result(gid, "lineage1", "lineage1.1", [],
                                                                                     tool="6.6.2")))
        if bad_db:
            gid = gs[1][0]
            (run / "results" / f"{gid}.results.json").write_text(json.dumps(result(gid, "lineage1", "lineage1.2", [],
                                                                                     db=OTHER_DB)))
        if extra_planned:
            (run / "plan.txt").write_text((run / "plan.txt").read_text() + "7002.1\n")
        self.run = run
        self.checks_ok = ct.write_checks(ct.run_counts(run), self.reading / "checks.txt")
        (self.reading / "reading.txt").write_text("python3 -I scripts/mechanism/tb_curve.py ... (SYNTHETIC)\n")
        if in_m_lists is not None:
            self.write_lists(in_m_lists)

    def write_lists(self, in_m):
        rows = [{"dataset": f"{SLUG}__{ab}", "species_slug": SLUG, "species": NAME, "antibiotic": ab,
                 "in_M": "1" if (in_m and ab == "rifampin") else "0", "in_A": "0",
                 "labels_depend_on": "the maintainer's rule" if (in_m and ab == "rifampin") else "none of the four"}
                for ab in ("rifampin", "amoxicillin", "ciprofloxacin")]
        c.write_tsv(self.lists, ["dataset", "species_slug", "species", "antibiotic", "in_M", "in_A", "labels_depend_on"],
                    rows)

    def write(self):
        out = self.top / "out"
        return ct.write(self.reading, self.lists, self.bed, out), out

    def edit(self, name, fn):
        path = self.reading / "m4" / name
        rows = c.read_tsv(path)
        header = path.read_text().split("\n", 1)[0].split("\t")
        c.write_tsv(path, header, fn(rows))


def stopped(case, reading, gate):
    with case.assertRaises(ct.Stop) as cm:
        reading.write()
    case.assertIn(f"gate {gate}", str(cm.exception))
    return str(cm.exception)


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def reading(self, **kw):
        sub = self.top / f"r{len(list(self.top.iterdir()))}"
        sub.mkdir()
        return Reading(sub, **kw)


class Check(Base):
    def test_a_clean_run_passes_and_names_no_genome(self):
        r = self.reading()
        self.assertTrue(r.checks_ok)
        text = (r.reading / "checks.txt").read_text()
        for line in ("planned: 26", "results: 26", "failed_lines: 0", "neither: 0", "versions_named: 26",
                     "versions_other: 0", "failed_by_reason: none", "wall_seconds: 1800"):
            self.assertIn(line + "\n", text)
        self.assertNotIn("FAIL", text)
        for g in genomes():
            self.assertNotIn(g[0], text)

    def test_failures_are_counted_by_reason(self):
        r = self.reading()
        ids = [g[0] for g in genomes()[:3]]
        for gid in ids:
            (r.run / "results" / f"{gid}.results.json").unlink()
        (r.run / "failed.tsv").write_text(f"{ids[0]}\t124\n{ids[1]}\tno_result\n{ids[2]}\t124\n")
        counts = ct.run_counts(r.run)
        self.assertEqual((counts["results"], counts["failed_lines"], counts["accounted"]), (23, 3, 26))
        self.assertEqual(counts["failed_by_reason"], "124=2,no_result=1")
        self.assertTrue(all(ok for _, ok in ct.checks_pass(counts)))

    def test_a_result_of_another_tool_version_fails_the_check(self):
        r = self.reading(bad_version=True)
        self.assertFalse(r.checks_ok)
        counts = ct.run_counts(r.run)
        self.assertEqual((counts["versions_named"], counts["versions_other"]), (25, 1))

    def test_a_result_of_another_database_commit_fails_the_check(self):
        r = self.reading(bad_db=True)
        self.assertFalse(r.checks_ok)
        counts = ct.run_counts(r.run)
        self.assertEqual((counts["versions_named"], counts["versions_other"]), (25, 1))

    def test_a_planned_genome_with_neither_fails_the_check(self):
        r = self.reading(extra_planned=True)
        self.assertFalse(r.checks_ok)
        self.assertEqual(ct.run_counts(r.run)["neither"], 1)


class Gates(Base):
    def test_a_clean_reading_writes_the_record(self):
        r = self.reading()
        files, out = r.write()
        self.assertEqual(files, ["README.md", "checks.txt", "gates.txt", "m4/m4-convergence.tsv", "m4/m4-coverage.tsv",
                                 "m4/m4-curve.tsv", "m4/m4-genes.tsv", "m4/m4-summary.json", "page-endpoints.tsv",
                                 "reading.txt"])
        for name in ct.M4_FILES:
            self.assertEqual((out / "m4" / name).read_bytes(), (r.reading / "m4" / name).read_bytes())
        for line in (out / "SHA256SUMS").read_text().splitlines():
            digest, name = line.split("  ", 1)
            self.assertEqual(ct.sha256_of(out / name), digest, name)
        self.assertEqual(len((out / "gates.txt").read_text().splitlines()), len(ct.GATES))
        self.assertTrue((out / "README.md").read_text().startswith("# The lineage curve's reading of 1 January 2026: "))
        ids = self.top / "species.json"
        ids.write_text(json.dumps({"genomes": [{"id": g[0]} for g in genomes()]}))
        with contextlib.redirect_stdout(io.StringIO()) as printed:
            self.assertEqual(scan.main(["--ids", str(ids), str(out)]), 0, printed.getvalue())
        self.assertIn("5 tables", printed.getvalue())

    def test_the_endpoints_hold_every_report_named_in_full(self):
        r = self.reading()
        _, out = r.write()
        rows = {x["antibiotic"]: x for x in c.read_tsv(out / "page-endpoints.tsv")}
        self.assertEqual(sorted(rows), ["amoxicillin", "ciprofloxacin", "rifampin"])
        self.assertEqual({x["species"] for x in rows.values()}, {NAME})
        rif = rows["rifampin"]
        self.assertEqual((rif["kind"], rif["tool_drugs"], rif["resistant_with_result"]), ("primary", "rifampicin", "24"))
        # synR: 3 changes over 3 lineages, in 15 genomes; synK: 1 change in 3 genomes, folded
        self.assertEqual((rif["catalogue_variants"], rif["genes_named"]), ("4", "1"))
        self.assertEqual((rif["without_determinant"], rif["share_without"]), ("8", f"{8 / 24:.12f}"))
        self.assertEqual((rif["groups_main"], rif["one_lineage_main"]), ("3", "4"))
        self.assertEqual(rows["amoxicillin"]["kind"], "none")
        self.assertEqual(rows["amoxicillin"]["resistant_with_result"], "")

    def test_a_dataset_in_m_is_marked_in_the_page_tables(self):
        r = self.reading(in_m_run=True)
        _, out = r.write()
        rows = {x["antibiotic"]: x for x in c.read_tsv(out / "page-endpoints.tsv")}
        self.assertEqual((rows["rifampin"]["in_M"], rows["ciprofloxacin"]["in_M"]), ("1", "0"))
        self.assertEqual(rows["rifampin"]["labels_depend_on"], "the maintainer's rule")
        printed = "\n".join(ct.tables(out)["reading"])
        self.assertIn("| rifampin | primary | M |", printed)
        self.assertIn("| ciprofloxacin | proxy | not M |", printed)

    def test_a_dataset_of_m_left_unmarked_is_refused(self):
        r = self.reading(in_m_run=False, in_m_lists=True)
        self.assertIn("rifampin", stopped(self, r, 4))

    def test_a_share_under_5_resistant_genomes_is_refused(self):
        r = self.reading()

        def small(rows):
            for x in rows:
                if x["grouping"] == "sub_lineage" and x["antibiotic"] == "rifampin" and x["group"] == "lineage1.1":
                    x["resistant"], x["without_determinant"] = "4", "1"
                    x["share_without"] = f"{1 / 4:.12f}"
            return rows
        r.edit("m4-coverage.tsv", small)
        stopped(self, r, 5)

    def test_a_curve_over_under_5_is_refused(self):
        r = self.reading()

        def small(rows):
            rows[0]["resistant"] = "4"
            return rows
        r.edit("m4-curve.tsv", small)
        stopped(self, r, 5)

    def test_a_convergence_share_under_5_resistant_genomes_with_a_group_is_refused(self):
        r = self.reading()
        path = r.reading / "m4" / "m4-summary.json"
        d = json.loads(path.read_text())
        for e in d["endpoints"]:
            if e["antibiotic"] == "rifampin":
                e["resistant_with_no_group"]["main_lineage"] = "20"
        path.write_text(json.dumps(d))
        convergence = c.read_tsv(r.reading / "m4" / "m4-convergence.tsv")
        shares = [x for x in convergence
                  if x["antibiotic"] == "rifampin" and x["grouping"] == "main_lineage" and x["share_covered"]]
        self.assertEqual(len(shares), 4)
        self.assertIn(": 4 rows write", stopped(self, r, 5))

    def test_a_variants_own_name_in_an_output_is_refused(self):
        for case in ("gene", "summary"):
            with self.subTest(case=case):
                r = self.reading()
                if case == "gene":
                    def named(rows):
                        for x in rows:
                            if x["gene"] == "synR":
                                x["gene"] = "synR_S450L"
                        return rows
                    r.edit("m4-genes.tsv", named)
                else:
                    path = r.reading / "m4" / "m4-summary.json"
                    d = json.loads(path.read_text())
                    d["note"] = "synR p.Ser450Leu"
                    path.write_text(json.dumps(d))
                stopped(self, r, 6)

    def test_a_gene_in_1_to_4_resistant_genomes_named_is_refused(self):
        r = self.reading()

        def named(rows):
            folded = [x for x in rows if x["gene"] == "other genes"]
            rows = [x for x in rows if x["gene"] != "other genes"]
            for x in folded:
                rows.append(dict(x, gene="synK"))
            return rows
        r.edit("m4-genes.tsv", named)
        stopped(self, r, 7)

    def test_a_result_of_another_tool_version_is_refused(self):
        r = self.reading(bad_version=True)
        stopped(self, r, 2)

    def test_a_result_of_another_database_commit_is_refused(self):
        r = self.reading(bad_db=True)
        stopped(self, r, 2)

    def test_a_summary_of_another_database_commit_is_refused(self):
        r = self.reading()
        path = r.reading / "m4" / "m4-summary.json"
        d = json.loads(path.read_text())
        self.assertEqual((d["tool"], d["database_commit"]), ("6.6.3", "3ed1cb99"))
        d["database_commit"] = OTHER_DB
        path.write_text(json.dumps(d))
        self.assertIn("m4-summary.json", stopped(self, r, 2))

    def test_a_planned_genome_with_neither_a_result_nor_a_failed_line_is_refused(self):
        r = self.reading(extra_planned=True)
        stopped(self, r, 3)

    def test_a_reading_with_another_file_under_m4_is_refused(self):
        r = self.reading()
        shutil.copyfile(r.reading / "m4" / "m4-genes.tsv", r.reading / "m4" / "m4-genes-old.tsv")
        stopped(self, r, 1)


class Printing(Base):
    def test_every_table_prints_from_the_record_alone(self):
        r = self.reading()
        _, out = r.write()
        t = ct.tables(out)
        self.assertEqual(sorted(t), ["curve_ends", "curves_main", "datasets", "genes", "reading", "spread", "without"])
        text = "\n".join(line for lines in t.values() for line in lines)
        self.assertIn("| rifampin | *synR* | 3 | 15 |", text)
        self.assertIn("other genes", text)
        self.assertNotIn("p.Ser", text)
        main = [x for x in t["curves_main"] if x.startswith("| rifampin |")]
        self.assertEqual([x.split(" | ")[1] for x in main], ["1 of 3", "2 of 3", "3 of 3"])
        self.assertTrue(main[-1].endswith("| 4 (4 to 4) | 4 |"), main[-1])
        reading = [x for x in t["reading"] if x.startswith("| rifampin |")][0]
        self.assertIn("| 8 (0.333) |", reading)


if __name__ == "__main__":
    unittest.main()
