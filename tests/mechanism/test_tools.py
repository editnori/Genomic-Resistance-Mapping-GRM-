"""scripts/mechanism/scan.py, drug_rules.py, amr_rule.py and run_species.sh (docs/mechanism-split.md).

Every genome id, accession, report and table is SYNTHETIC.
"""

import contextlib
import hashlib
import io
import json
import re
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

    def scan(self, text):
        (self.top / "res" / "t.tsv").write_text(text)
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

    def test_figures_and_other_ids_pass(self):
        rc, out = self.scan("0.500000000000\t9999.1\t19999.11\t9999.110\tunder 5\n")
        self.assertEqual(rc, 0, out)


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

    def test_the_port_names_the_kit_files_it_copies(self):
        self.assertIn("amr_rule.ADAPTER_SHA256", self.text)
        self.assertIn("amr_rule.GRM_EVAL_SHA256", self.text)
        self.assertEqual(len(amr.ADAPTER_SHA256), 64)
        self.assertEqual(hashlib.sha256(b"").hexdigest() != amr.ADAPTER_SHA256, True)


if __name__ == "__main__":
    unittest.main()
