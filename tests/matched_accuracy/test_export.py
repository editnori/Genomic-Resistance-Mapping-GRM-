"""The export and the row lists of scripts/matched-accuracy/stage1.py, on the SYNTHETIC release of made_up.py."""

import hashlib
import re
import stat
import unittest

from tests.matched_accuracy import common

s1 = common.stage1
FILES = ["checks.tsv", "counts.tsv", "inputs.tsv", "shares.tsv", "typing.tsv"]


class Export(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.made = common.Run()
        cls.rel = common.release()

    @classmethod
    def tearDownClass(cls):
        cls.made.close()

    def ran(self):
        self.assertEqual(self.made.code, 0, self.made.lines)

    def test_the_export_holds_five_files_and_the_local_directory_two(self):
        self.ran()
        self.assertEqual(sorted(p.name for p in (self.made.out / "export").iterdir()), FILES)
        self.assertEqual(sorted(p.name for p in (self.made.out / "local").iterdir()), ["reading-A.tsv", "reading-B.tsv"])
        self.assertEqual(sorted(p.name for p in self.made.out.iterdir()), ["export", "local"])

    def test_no_export_file_holds_a_decimal_point_a_genome_id_or_an_accession(self):
        self.ran()
        ids = {s["id"] for r in self.rel["reports"] for s in r["samples"]} | {"gz-zz-01"}
        accessions = {p for g in self.rel["genomes"] for p in s1.bioproject_parts(g["bioproject_accession"])}
        self.assertEqual(len(ids), 91)
        self.assertEqual(len(accessions), 10)
        for name in FILES:
            text = self.made.export(name).decode("utf-8")
            with self.subTest(file=name):
                self.assertIsNone(re.search(r"\d\.\d", text))
                self.assertIsNone(re.search(r"\b(gz|pz)-", text))
                self.assertEqual([x for x in ids | accessions if x in text], [])

    def test_the_row_lists_are_mode_0600_in_a_directory_of_mode_0700(self):
        self.ran()
        local = self.made.out / "local"
        self.assertEqual(stat.S_IMODE(local.stat().st_mode), 0o700)
        for name in ("reading-A.tsv", "reading-B.tsv"):
            self.assertEqual(stat.S_IMODE((local / name).stat().st_mode), 0o600, name)

    def test_the_row_list_is_sorted_utf8_lf_and_tab_separated(self):
        self.ran()
        for reading in ("A", "B"):
            data = self.made.row_list(reading)
            with self.subTest(reading=reading):
                self.assertNotIn(b"\r", data)
                self.assertTrue(data.endswith(b"\n"))
                lines = data.decode("utf-8").split("\n")[:-1]
                self.assertEqual(lines, sorted(lines))
                self.assertTrue(all(len(line.split("\t")) == 8 for line in lines))
                self.assertEqual(len(lines), len(set(lines)))

    def test_the_inputs_table_gives_each_sha256_and_ends_with_scores_opened_no(self):
        self.ran()
        data = self.made.export("inputs.tsv")
        self.assertTrue(data.endswith(b"\nscores opened: no\n"))
        inputs = {r["input"]: r for r in self.made.table("inputs.tsv")}
        for reading in ("A", "B"):
            self.assertEqual(inputs[f"row list, reading {reading}"]["sha256"],
                             hashlib.sha256(self.made.row_list(reading)).hexdigest())
        top = self.made.top
        for name, path in (("lab records", "genome_amr.jsonl"), ("drug names", "drug_names.json"),
                           ("label table", "phenotypes.jsonl"), ("genome table, second copy", "master-genome-meta.jsonl.gz"),
                           ("genome table, 5 Nov", "genomes.jsonl")):
            self.assertEqual(inputs[name]["sha256"], hashlib.sha256((top / path).read_bytes()).hexdigest(), name)
        self.assertEqual(inputs["reports"]["sha256"], self.made.pins["REPORT_SUMS_SHA256"])

    def test_the_dropped_figure_is_named_as_recorded_never_as_measured(self):
        self.ran()
        note = {r["input"]: r["note"] for r in self.made.table("inputs.tsv")}["lab records"]
        self.assertIn("252148 records dropped at the cut is a recorded figure of ruling R49, not measured here", note)
        for name in FILES:
            text = self.made.export(name).decode("utf-8")
            for m in re.finditer(r"measured", text):
                self.assertEqual(text[m.start() - 4:m.start()], "not ", name)

    def test_nothing_under_local_is_copied_into_the_export(self):
        self.ran()
        export = b"".join(self.made.export(n) for n in FILES)
        for reading in ("A", "B"):
            for line in self.made.row_list(reading).split(b"\n")[:-1]:
                self.assertNotIn(line, export)
                self.assertNotIn(line.split(b"\t")[1], export)

    def test_the_export_has_no_row_per_genome_pair_group_or_bioproject(self):
        self.ran()
        for name in ("counts.tsv", "typing.tsv"):
            levels = {r["level"] for r in self.made.table(name)}
            self.assertEqual(levels, {"dataset", "species", "overall"}, name)
        self.assertEqual(len(self.made.table("counts.tsv")), (12 + 4 + 1) * 2)
        self.assertEqual(len(self.made.table("typing.tsv")), (11 + 4 + 1) * 2 * 2)
        self.assertEqual([r["reading"] for r in self.made.table("shares.tsv")], ["A", "B"])

    def test_both_readings_are_always_computed_and_written(self):
        self.ran()
        for name in ("counts.tsv", "typing.tsv", "shares.tsv"):
            self.assertEqual({r["reading"] for r in self.made.table(name)}, {"A", "B"}, name)


class WithoutHc20(common.Case, unittest.TestCase):
    def test_without_the_5_nov_table_the_typing_gives_the_mlst_step_alone(self):
        run = self.run_stage1(hc20=False)
        self.assertEqual(run.code, 0, run.lines)
        self.assertEqual({r["step"] for r in run.table("typing.tsv")}, {"mlst"})
        row = {r["input"]: r for r in run.table("inputs.tsv")}["genome table, 5 Nov"]
        self.assertEqual((row["sha256"], row["count"], row["note"]), ("-", "0", "not given: the HC20 read is dropped"))


if __name__ == "__main__":
    unittest.main()
