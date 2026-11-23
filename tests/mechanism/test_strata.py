"""scripts/mechanism/strata.py: block M1, recall by mechanism (docs/mechanism-split.md).

Every report, table, genome id and count is SYNTHETIC.
"""

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

strata = c.load("strata")
scan = c.load("scan")


def run(top, out="out"):
    top = Path(top)
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as err:
        rc = strata.main(["--species", c.SPECIES, "--reports", str(top / "reports"), "--lists", str(top / "lists.tsv"),
                          "--drug-rules", str(top / "drug_rules.tsv"), "--tables", str(top / "tool" / c.SPECIES),
                          "--calls", str(top / "calls.tsv"), "--callers", str(top / "callers.tsv"),
                          "--out", str(top / out)])
    return rc, err.getvalue(), top / out


def cells(out):
    return {(r["dataset"], r["split"], r["row"], r["stratum"]): r for r in c.read_tsv(out / "strata.tsv")}


def mixed():
    """One dataset with every kind of genome."""
    g = (c.genomes("9001", 10, "point", True) + c.genomes("9002", 12, "gene", True)
         + c.genomes("9003", 6, "both", True) + c.genomes("9004", 9, "neither", True)
         + c.genomes("9005", 6, "partial", True) + c.genomes("9006", 5, "disrupt", True)
         + c.genomes("9007", 6, "gene", False) + c.genomes("9008", 20, "neither", False))
    return g


class Strata(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)

    def go(self, datasets=None, **kw):
        datasets = datasets or {c.DATASET: (mixed(), {})}
        c.write_species(self.top, datasets, **kw)
        return run(self.top)

    def test_point_rows_make_point_only(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        got = cells(out)
        # 10 POINT genomes and 5 POINT_DISRUPT genomes: both read with point (the lead's ruling on question 1)
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_only")]["resistant"], "15")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "both")]["resistant"], "6")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "gene_only")]["resistant"], "18")

    def test_stratum_comes_from_the_table_not_the_label(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        got = cells(out)
        # gene 12 + partial 6 resistant (POINT_DISRUPT is read with point); the 6 susceptible genomes with a gene row are gene only too
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "gene_only")]["resistant"], "18")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "gene_only")]["susceptible"], "6")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "neither")]["susceptible"], "20")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "neither")]["resistant"], "9")

    def test_partial_calls_and_point_disrupt_are_counted_apart(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        got = cells(out)
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "partial_only:gene_only")]["resistant"], "6")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "partial_only:point_only")]["resistant"], "0")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "with_point_disrupt:point_only")]["resistant"], "5")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "with_point_disrupt:gene_only")]["resistant"], "0")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_disrupt_alone:point_only")]["resistant"], "5")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_disrupt_alone:both")]["resistant"], "0")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "all")]["resistant"], "48")

    def test_point_disrupt_is_read_with_point_and_counted_beside(self):
        # a genome with a POINT row and a POINT_DISRUPT row: point only, with POINT_DISRUPT, not alone
        c.ROWS["point_and_disrupt"] = c.ROWS["point"] + c.ROWS["disrupt"]
        self.addCleanup(c.ROWS.pop, "point_and_disrupt")
        g = (c.genomes("9101", 10, "point_and_disrupt", True) + c.genomes("9102", 10, "disrupt", True)
             + c.genomes("9103", 10, "gene", True) + c.genomes("9104", 5, "neither", False))
        rc, err, out = self.go({c.DATASET: (g, {})})
        self.assertEqual(rc, 0, err)
        got = cells(out)
        key = (c.DATASET, "random", "gbt")
        self.assertEqual(got[key + ("point_only",)]["resistant"], "20")
        self.assertEqual(got[key + ("gene_only",)]["resistant"], "10")
        self.assertEqual(got[key + ("with_point_disrupt:point_only",)]["resistant"], "20")
        self.assertEqual(got[key + ("point_disrupt_alone:point_only",)]["resistant"], "10")

    def test_recall_only_over_ten_resistant_and_small_counts_hidden(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        got = cells(out)
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_only")]["recall"], "1.000000000000")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "neither")]["recall"], "")  # 9 resistant
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "both")]["recall"], "")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_only")]["susceptible"], "0")
        pooled = {(r["split"], r["row"], r["stratum"]): r for r in c.read_tsv(out / "pooled.tsv")}
        self.assertEqual(pooled[("lineage", "scm", "neither")]["resistant_pairs"], "9")
        self.assertEqual(pooled[("lineage", "scm", "neither")]["recall"], "")

    def test_scores_are_read_from_their_own_split_and_learner(self):
        def scores(split, learner, i, g):
            hit = split == "lineage" and learner == "cart" and g["kind"] in ("point", "disrupt")
            return 0.9 if hit else 0.1
        rc, err, out = self.go({c.DATASET: (mixed(), {"scores": scores})})
        self.assertEqual(rc, 0, err)
        got = cells(out)
        self.assertEqual(got[(c.DATASET, "lineage", "cart", "point_only")]["resistant_called"], "15")
        self.assertEqual(got[(c.DATASET, "lineage", "scm", "point_only")]["resistant_called"], "0")
        self.assertEqual(got[(c.DATASET, "random", "cart", "point_only")]["resistant_called"], "0")
        self.assertEqual(got[(c.DATASET, "lineage", "cart", "point_only")]["recall"], "1.000000000000")
        self.assertEqual(got[(c.DATASET, "random", "cart", "point_only")]["recall"], "0.000000000000")

    def test_program_and_catalog_rows(self):
        rc, err, out = self.go()
        self.assertEqual(rc, 0, err)
        got = cells(out)
        self.assertEqual(got[(c.DATASET, "lineage", "program", "neither")]["resistant_called"], "0")
        self.assertEqual(got[(c.DATASET, "lineage", "program", "gene_only")]["susceptible_called"], "6")
        self.assertEqual(got[(c.DATASET, "lineage", "catalog", "point_only")]["resistant_called"], "15")

    def test_hard_check_fails_when_metrics_differ(self):
        c.write_species(self.top, {c.DATASET: (mixed(), {})})
        path = self.top / "reports" / f"{c.DATASET}.json"
        rep = json.loads(path.read_text())
        rep["splits"][1]["learners"][2]["metrics"]["tp"] += 1
        path.write_text(json.dumps(rep))
        rc, err, out = run(self.top)
        self.assertEqual(rc, 0, err)
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["hard_check"], d["entered"], d["learner_rows_equal"], d["learner_rows"]), ("0", "0", "11", "12"))
        self.assertIn("counts from scores differ from metrics", d["why"])
        self.assertEqual(c.read_tsv(out / "strata.tsv"), [])

    def test_catalog_calls_must_equal_catalog_resistant(self):
        c.write_species(self.top, {c.DATASET: (mixed(), {})})
        path = self.top / "reports" / f"{c.DATASET}.json"
        rep = json.loads(path.read_text())
        rep["samples"][0]["catalog_resistant"] = not rep["samples"][0]["catalog_resistant"]
        path.write_text(json.dumps(rep))
        rc, err, out = run(self.top)
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["hard_check"], d["entered"]), ("0", "0"))

    def test_program_call_must_equal_the_adapters(self):
        def flip(d, g, call):
            return "susceptible" if g["id"] == "9001.1" else call
        rc, err, out = self.go(adapter_change=flip)
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["program_check"], d["entered"]), ("0", "0"))

    def test_counts_must_equal_the_record(self):
        def bump(rec):
            rec["fn"] += 1
        rc, err, out = self.go(record_change=bump)
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["record_check"], d["entered"]), ("0", "0"))

    def test_rule_must_equal_the_drug_map(self):
        rc, err, out = self.go(rules_change=lambda d, r: "C:TETRACYCLINE")
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["rule_check"], d["entered"]), ("0", "0"))

    def test_a_genome_the_program_failed_on_is_left_out(self):
        g = mixed() + [{"id": "9009.1", "kind": None, "resistant": True}]
        rc, err, out = self.go({c.DATASET: (g, {})})
        self.assertEqual(rc, 0, err)
        d = c.read_tsv(out / "datasets.tsv")[0]
        self.assertEqual((d["genomes"], d["with_table"], d["entered"]), ("75", "74", "1"))
        self.assertEqual(cells(out)[(c.DATASET, "lineage", "scm", "all")]["resistant"], "48")

    def test_a_dataset_in_m_is_set_apart_with_its_label_rule(self):
        other = "synthia_coli__drug_b"
        lists = {other: c.lists_row(other, in_m=True, depends="intermediate read as resistant")}
        rc, err, out = self.go({c.DATASET: (mixed(), {}), other: (mixed(), {})}, lists=lists)
        self.assertEqual(rc, 0, err)
        rows = c.read_tsv(out / "strata.tsv")
        self.assertEqual({r["set"] for r in rows if r["dataset"] == other}, {"M"})
        self.assertEqual({r["labels_depend_on"] for r in rows if r["dataset"] == other}, {"intermediate read as resistant"})
        self.assertEqual({r["set"] for r in rows if r["dataset"] == c.DATASET}, {"primary"})

    def test_no_genome_id_in_any_output(self):
        rc, err, out = self.go()
        ids = self.top / "species.json"
        ids.write_text(json.dumps({"genomes": [{"id": g["id"]} for g in mixed()]}))
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(scan.main(["--ids", str(ids), str(out)]), 0)

    def test_a_missing_table_refuses_without_naming_the_genome(self):
        c.write_species(self.top, {c.DATASET: (mixed(), {})})
        (self.top / "tool" / c.SPECIES / "9001.3.tsv").unlink()
        rc, err, out = run(self.top)
        self.assertEqual(rc, 1)
        self.assertIn("has no tool output", err)
        self.assertNotIn("9001.3", err)


if __name__ == "__main__":
    unittest.main()
