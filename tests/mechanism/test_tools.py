"""scripts/mechanism/scan.py, drug_rules.py, amr_rule.py and run_species.sh (docs/mechanism-split.md).

Every genome id, accession, report and table is SYNTHETIC.
"""

import contextlib
import hashlib
import io
import json
import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from tests.mechanism import common as c

scan = c.load("scan")
rules = c.load("drug_rules")
amr = c.load("amr_rule")
RUNNER = c.KIT / "run_species.sh"


def quiet(fn, *a):
    with contextlib.redirect_stdout(io.StringIO()) as out, contextlib.redirect_stderr(io.StringIO()) as err:
        rc = fn(*a)
    return rc, out.getvalue() + err.getvalue()


class Scan(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)
        (self.top / "ids.json").write_text(json.dumps({"genomes": [{"id": "9999.11"}, {"id": "9999.12"}]}))
        (self.top / "res").mkdir()

    def scan(self, text, name="t.txt"):
        (self.top / "res" / name).write_text(text)
        return quiet(scan.main, ["--ids", str(self.top / "ids.json"), str(self.top / "res")])

    def test_a_genome_id_fails_and_is_not_printed(self):
        rc, out = self.scan("a\t9999.11\n")
        self.assertEqual(rc, 1)
        self.assertNotIn("9999.11", out)
        self.assertIn("1 genome ids", out)

    def test_a_bioproject_accession_fails(self):
        rc, out = self.scan("a\tPRJNA123456\n")
        self.assertEqual(rc, 1)
        self.assertNotIn("PRJNA123456", out)

    def test_every_bioproject_prefix_fails(self):
        for acc in ("PRJEB12345", "PRJDB678", "PRJNA1"):
            rc, out = self.scan(f"a\t{acc}\n")
            self.assertEqual(rc, 1, acc)
            self.assertNotIn(acc, out)

    def test_figures_and_other_ids_pass(self):
        rc, out = self.scan("0.500000000000\t9999.1\t19999.11\t9999.110\tunder 5\n")
        self.assertEqual(rc, 0, out)


strata = c.load("strata")
history = c.load("history")
M4_GENES = ["species", "antibiotic", "in_M", "gene", "catalogue_variants", "resistant_with_gene"]


class ScanTables(unittest.TestCase):
    """Rows planted in each table the scan reads, each of which could stand
    for 1 to 4 genomes or for one BioProject."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)
        (self.top / "ids.json").write_text(json.dumps({"genomes": [{"id": "9999.11"}]}))

    def scan(self, name, header, rows):
        res = self.top / "res"
        shutil.rmtree(res, ignore_errors=True)
        res.mkdir()
        lines = ["\t".join(header)] + ["\t".join(r.get(h, "") for h in header) for r in rows]
        (res / name).write_text("\n".join(lines) + "\n")
        return quiet(scan.main, ["--ids", str(self.top / "ids.json"), str(res)])

    def assertHit(self, name, header, row, rule):
        rc, out = self.scan(name, header, [row])
        self.assertEqual(rc, 1, out)
        self.assertIn(rule, out)

    def test_a_count_of_1_to_4_shown_fails(self):
        ok = {"resistant": "0", "resistant_called": "under 5", "susceptible": "12", "susceptible_called": "5"}
        self.assertEqual(self.scan("strata.tsv", strata.STRATA_HEADER, [ok])[0], 0)
        for bad in ("1", "4", "-1", "05", "4.0", "x"):
            self.assertHit("strata.tsv", strata.STRATA_HEADER, dict(ok, susceptible_called=bad),
                           "susceptible_called: not 0, under 5 or 5 and more")

    def test_a_recall_over_under_10_fails(self):
        row = {"resistant": "9", "resistant_called": "9", "susceptible": "0", "susceptible_called": "0",
               "recall": "1.000000000000"}
        self.assertHit("strata.tsv", strata.STRATA_HEADER, row, "recall: written over under 10")
        self.assertEqual(self.scan("strata.tsv", strata.STRATA_HEADER, [dict(row, resistant="10")])[0], 0)

    def test_a_pooled_row_over_one_dataset_with_one_genome_fails_unless_masked(self):
        row = {"datasets": "1", "resistant_pairs": "1", "resistant_called": "1", "susceptible_pairs": "0",
               "susceptible_called": "0", "recall": ""}
        self.assertHit("pooled.tsv", strata.POOLED_HEADER, row, "resistant_pairs: not 0, under 5 or 5 and more")
        masked = dict(row, resistant_pairs="under 5", resistant_called="under 5")
        self.assertHit("pooled.tsv", strata.POOLED_HEADER, dict(masked, recall="1.000000000000"),
                       "recall: written over under 10")
        self.assertEqual(self.scan("pooled.tsv", strata.POOLED_HEADER, [masked])[0], 0)

    def m2_row(self, **kw):
        row = {k: "0" for k in ("genomes", "earlier_genomes", "later_genomes", "earlier_units", "later_units",
                                "earlier_point_disrupt", "later_point_disrupt")}
        row.update(genomes="40", earlier_genomes="20", later_genomes="20", guard="pooled", passes="")
        row.update(kw)
        return row

    def test_a_share_beside_under_5_genomes_fails(self):
        self.assertEqual(self.scan("m2-cells.tsv", history.CELLS_HEADER,
                                   [self.m2_row(earlier_point="0.5", d_point="0.1")])[0], 0)
        self.assertHit("m2-cells.tsv", history.CELLS_HEADER,
                       self.m2_row(earlier_genomes="under 5", earlier_point="1.000000000000"),
                       "earlier_point: written beside a count under 5")

    def test_a_voted_row_that_does_not_pass_with_a_share_fails(self):
        voted = self.m2_row(guard="bioproject", passes="0", earlier_units="under 5", later_units="under 5")
        self.assertEqual(self.scan("m2-cells.tsv", history.CELLS_HEADER, [voted])[0], 0)
        self.assertHit("m2-cells.tsv", history.CELLS_HEADER, dict(voted, later_gene="0.250000000000"),
                       "a voted row that does not pass writes a share")
        self.assertEqual(self.scan("m2-cells.tsv", history.CELLS_HEADER,
                                   [dict(voted, passes="1", later_gene="0.25")])[0], 0)

    def test_a_family_seen_in_under_5_genomes_named_fails(self):
        row = {"family": "qnrB", "earlier_genomes": "under 5", "later_genomes": "0", "earlier_mix": "",
               "later_mix": "0.000000000000"}
        self.assertHit("m2-families.tsv", history.FAMILIES_HEADER, row, "a family under 5 genomes is named")
        fold = dict(row, family=history.FOLDED_FAMILIES)
        self.assertEqual(self.scan("m2-families.tsv", history.FAMILIES_HEADER, [fold])[0], 0)
        self.assertHit("m2-families.tsv", history.FAMILIES_HEADER, dict(fold, earlier_mix="0.1"),
                       "earlier_mix: written beside a count under 5")

    def test_a_gene_seen_in_under_5_genomes_named_fails(self):
        row = {"gene": "synR", "catalogue_variants": "1", "resistant_with_gene": "under 5"}
        self.assertHit("m4-genes.tsv", M4_GENES, row, "a gene under 5 genomes is named")
        fold = dict(row, gene="genes seen in under 5 genomes")
        self.assertEqual(self.scan("m4-genes.tsv", M4_GENES, [fold])[0], 0)

    def test_a_table_of_another_name_or_without_a_column_fails(self):
        self.assertHit("per-bioproject.tsv", ["a"], {"a": "1"}, "a table the scan does not know")
        self.assertHit("strata.tsv", [h for h in strata.STRATA_HEADER if h != "recall"], {"resistant": "5"},
                       "a column the scan reads is missing")
        self.assertEqual(self.scan("timing.tsv", ["step", "wall_seconds", "peak_kb"], [{"step": "1"}])[0], 0)

    def test_the_output_names_the_rule_not_the_value(self):
        rc, out = self.scan("strata.tsv", strata.STRATA_HEADER, [{"resistant": "3"}])
        self.assertEqual(rc, 1)
        self.assertNotIn("\t3", out)
        self.assertNotRegex(out, r"\b3\b")


class Rules(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.top = Path(self.tmp.name)
        rule = [["C:QUINOLONE", "S:FLUOROQUINOLONE"], ["C:AMINOGLYCOSIDE"]]
        c.write_species(self.top, {c.DATASET: (c.genomes("9001", 3, "gene", True), {"rule": rule})})

    def test_the_drug_map_round_trips_and_is_checked(self):
        args = ["--reports", str(self.top / "reports"), "--lists", str(self.top / "lists.tsv")]
        rc, out = quiet(rules.main, args + ["--out", str(self.top / "r.tsv")])
        self.assertEqual(rc, 0, out)
        row = c.read_tsv(self.top / "r.tsv")[0]
        self.assertEqual(row["rule"], "C:QUINOLONE | S:FLUOROQUINOLONE & C:AMINOGLYCOSIDE")
        self.assertEqual(amr.parse_rule(row["rule"]), [["C:QUINOLONE", "S:FLUOROQUINOLONE"], ["C:AMINOGLYCOSIDE"]])
        self.assertEqual(quiet(rules.main, args + ["--check", str(self.top / "r.tsv")])[0], 0)
        text = (self.top / "r.tsv").read_text().replace("C:AMINOGLYCOSIDE", "C:TETRACYCLINE")
        (self.top / "r.tsv").write_text(text)
        rc, out = quiet(rules.main, args + ["--check", str(self.top / "r.tsv")])
        self.assertEqual(rc, 1)
        self.assertIn(c.DATASET, out)

    def test_a_bad_token_is_refused(self):
        with self.assertRaises(amr.Refusal):
            amr.parse_rule("X:QUINOLONE")


class Rule(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.dir = Path(self.tmp.name)

    def rows(self, kind):
        c.write_table(self.dir / f"{kind}.tsv", kind)
        return amr.read_table(self.dir / f"{kind}.tsv")

    def test_every_group_must_be_met(self):
        both = [["C:QUINOLONE"], ["C:TETRACYCLINE"]]
        self.assertFalse(amr.program_call(self.rows("gene"), both))
        self.assertTrue(amr.program_call(self.rows("gene") + self.rows("neither"), both))

    def test_tool_default_counts_point_and_disrupt_in_core_only(self):
        self.assertTrue(amr.program_call(self.rows("point"), c.RULE))
        self.assertTrue(amr.program_call(self.rows("disrupt"), c.RULE))
        self.assertFalse(amr.program_call(self.rows("plus_only"), c.RULE))
        self.assertFalse(amr.program_call(self.rows("point"), c.RULE, "genes_core"))
        self.assertFalse(amr.program_call(self.rows("partial"), c.RULE, "tool_default_intact"))

    def test_an_unknown_subtype_or_a_bad_table_is_refused(self):
        p = self.dir / "bad.tsv"
        p.write_text("\t".join(c.TABLE_HEADER) + "\n" + "\t".join(
            {"Type": "AMR", "Subtype": "NEWTHING", "Scope": "core"}.get(h, "x") for h in c.TABLE_HEADER) + "\n")
        with self.assertRaises(amr.Refusal):
            amr.read_table(p)
        p.write_text("")
        with self.assertRaises(amr.Refusal):
            amr.read_table(p)

    def test_susceptible_alleles_count_for_none(self):
        p = self.dir / "s.tsv"
        p.write_text("\t".join(c.TABLE_HEADER) + "\n" + "\t".join(
            {"Type": "AMR", "Subtype": "AMR-SUSCEPTIBLE", "Scope": "core", "Class": "QUINOLONE"}.get(h, "x")
            for h in c.TABLE_HEADER) + "\n")
        self.assertEqual(amr.read_table(p), [])


class Runner(unittest.TestCase):
    text = RUNNER.read_text()

    def test_kit_files_are_pinned(self):
        for sha in ("c67155f598dbb1ebe1a4ab908e793124ebc058430a33c1c67d7e2058910256de",
                    "da07f277a7d271f1ae006919d2c2df4185737526348302f665a9bcb3a083d0e6",
                    "9291b83f4c7d01a0384f35446e6de8f986dff5f659439ee3d673bccdd09427ba",
                    c.load("history").META_SHA256):
            self.assertIn(sha, self.text)

    def test_steps_in_order_and_the_scan_before_done(self):
        steps = re.findall(r"^step (\w+) ", self.text, re.M)
        self.assertEqual(steps, ["release", "amrkit", "db", "tables", "adapter", "rules", "strata", "history"])
        self.assertLess(self.text.index("scan.py"), self.text.index('> "$OUT/DONE"'))
        self.assertIn("trap 'fail at' ERR", self.text)
        self.assertIn("flock -n 9", self.text)

    def scan_step(self, scan_rc):
        """Run the script's scan step alone, with a stand-in scan.py that
        writes its arguments and exits `scan_rc`."""
        block = self.text[self.text.index("STEP=scan\n"):self.text.index('touch "$OUT/flags/scan.done"\n')]
        block += 'touch "$OUT/flags/scan.done"\n'
        tmp = Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, tmp)
        (tmp / "here").mkdir()
        (tmp / "here" / "scan.py").write_text(
            "import sys\nopen(sys.argv[-1] + '/../args', 'w').write(' '.join(sys.argv[1:]))\n"
            f"sys.exit({scan_rc})\n")
        out = tmp / "out"
        (out / "flags").mkdir(parents=True)
        (out / "results").mkdir()
        (out / "logs").mkdir()
        script = (f"set -Eeuo pipefail\nsay() {{ :; }}\nHERE={tmp / 'here'}\nR={tmp}\nSP=synthia\nOUT={out}\n" + block)
        r = subprocess.run(["bash", "-c", script], capture_output=True, text=True, timeout=60)
        return r.returncode, (out / "flags" / "scan.done").exists(), (out / "args").read_text().split()

    def test_the_scan_reads_results_and_logs(self):
        rc, done, args = self.scan_step(0)
        self.assertEqual((rc, done), (0, True))
        self.assertTrue(args[-2].endswith("/out/results") and args[-1].endswith("/out/logs"), args)

    def test_a_failed_scan_stops_the_run_before_its_flag(self):
        rc, done, _ = self.scan_step(1)
        self.assertNotEqual(rc, 0)
        self.assertFalse(done)

    def test_the_port_names_the_kit_files_it_copies(self):
        self.assertIn("amr_rule.ADAPTER_SHA256", self.text)
        self.assertIn("amr_rule.GRM_EVAL_SHA256", self.text)
        self.assertEqual(len(amr.ADAPTER_SHA256), 64)
        self.assertEqual(hashlib.sha256(b"").hexdigest() != amr.ADAPTER_SHA256, True)


if __name__ == "__main__":
    unittest.main()
