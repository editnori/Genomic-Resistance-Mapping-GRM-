"""The export of scripts/matched-accuracy/stage2.py on the SYNTHETIC designed release of made_up2.py: its files,
integers and names alone, both readings, no row per genome, pair, group or BioProject, and the scan before any
file is written. Every id and BioProject here is made up (gz-<species>-<n>, and made_up2.bioproject)."""

import os
import re
import unittest
from pathlib import Path
from unittest import mock

from tests.matched_accuracy import made_up2, stage2_common as c
from tests.matched_accuracy.stage2_common import stage1, stage2

NUMBERS = {
    "measure.tsv": stage2.MEASURE[3:],
    "interval.tsv": ("lower_ppm", "upper_ppm", "species", "draws"),
    "cells.tsv": stage2.CELL_COLUMNS[3:],
    "typing.tsv": stage2.TYPING_COLUMNS[5:],
    "null_fits.tsv": stage2.NULL_COLUMNS[3:],
    "inputs.tsv": ("count",),
}
WORDS = {
    "reading": set(stage2.READINGS), "level": set(stage2.LEVELS), "step": {"mlst", "mlst_hc20"},
    "match": {"mlst", "mlst_hc20", "primary"}, "read": {"primary", "mlst", "mlst_hc20", "null_fits"},
    "resistant_member_call": {"resistant", "susceptible"} | set(stage2.APART),
    "susceptible_member_call": {"resistant", "susceptible"} | set(stage2.APART),
    "chance": {"above", "below", "unresolved", ""},
}


class Export(c.Case, unittest.TestCase):
    def setUp(self):
        self.made = self.run_stage2()
        self.assert_passes(self.made)

    def test_the_output_holds_the_line_and_the_six_files_and_nothing_under_local(self):
        self.assertEqual(sorted(p.name for p in self.made.out.iterdir()), ["export", stage2.OPENED_FILE])
        self.assertEqual(sorted(p.name for p in (self.made.out / "export").iterdir()), c.EXPORT)
        self.assertEqual((self.made.out / stage2.OPENED_FILE).read_text(), stage2.OPENED + "\n")
        self.assertFalse((self.made.out / "local").exists())

    def test_no_export_file_holds_a_decimal_point(self):
        for name in c.EXPORT:
            self.assertNotIn(b".", self.made.export(name), name)

    def test_every_cell_is_an_integer_a_name_of_the_release_or_a_fixed_word(self):
        names = {"all"} | {r["name"] for r in self.made.table("measure.tsv")}
        self.assertEqual({n for n in names if n != "all"},
                         {f"synthia_{n}" for n in ("kappa", "lambda", "mu", "nu", "xi", "omicron", "rho", "tau")} |
                         {d["id"] for d in self.made.release["reports"]})
        for file in c.EXPORT:
            for r in self.made.table(file):
                for column, value in r.items():
                    with self.subTest(file=file, column=column):
                        if column in NUMBERS[file]:
                            self.assertRegex(value, r"\A(?:-?[0-9]+)?\Z")
                        elif column in WORDS:
                            self.assertIn(value, WORDS[column])
                        elif column == "name":
                            self.assertIn(value, names)
                        else:
                            self.assertEqual(file, "inputs.tsv")
                            self.assertIn(column, ("input", "sha256", "note"))

    def test_no_genome_id_bioproject_or_shape_in_any_file(self):
        ids = {s["id"] for d in self.made.release["reports"] for s in d["samples"]}
        projects = {g["bioproject_accession"] for g in self.made.release["genomes"]}
        self.assertGreater(len(ids), 80)
        for name in c.EXPORT:
            data = self.made.export(name)
            self.assertEqual(stage1.hits(data, ids, projects), 0, name)
            self.assertNotIn(b"gz-", data, name)
            self.assertIsNone(re.search(rb"PRJ", data, re.IGNORECASE), name)

    def test_no_row_per_genome_pair_group_or_bioproject(self):
        # One row a level and name: overall, 8 species and 10 datasets, for each reading.
        per = 2 * (1 + 8 + 10)
        self.assertEqual(len(self.made.table("measure.tsv")), per)
        self.assertEqual(len(self.made.table("null_fits.tsv")), per)
        self.assertEqual(len(self.made.table("typing.tsv")), 5 * per)
        self.assertEqual(len(self.made.table("cells.tsv")), 2 * 6)
        self.assertEqual(len(self.made.table("interval.tsv")), 2 * 7)
        for name in c.EXPORT:
            head = self.made.export(name).decode().split("\n", 1)[0].split("\t")
            self.assertEqual([h for h in head if re.search(r"genome|bioproject|group\b|pair\b|gene|feature", h)], [],
                             name)

    def test_both_readings_are_written(self):
        for name in [n for n in c.EXPORT if n != "inputs.tsv"]:
            self.assertEqual({r["reading"] for r in self.made.table(name)}, {"A", "B"}, name)

    def test_inputs_give_numpy_as_three_integers_and_end_with_the_line(self):
        rows = {r["input"]: r for r in self.made.table("inputs.tsv")}
        self.assertEqual([rows[f"numpy version, {p}"]["count"] for p in ("major", "minor", "micro")], ["1", "26", "4"])
        self.assertTrue(self.made.export("inputs.tsv").endswith(b"\n" + stage2.OPENED.encode() + b"\n"))
        self.assertEqual(rows["row list, reading A"]["sha256"], c.sha(self.made.row_list("A").read_bytes()))
        self.assertEqual(rows["count table of stage 1"]["sha256"], c.sha(self.made.counts().read_bytes()))
        self.assertEqual(rows["typing table of stage 1"]["sha256"], c.sha(self.made.typing_table().read_bytes()))
        self.assertEqual(rows["row list, reading A"]["count"], "82")


class Scan(c.Case, unittest.TestCase):
    def test_the_scan_gets_every_read_id_and_bioproject_and_runs_its_positive_control(self):
        calls = []
        real = stage1.scan

        def spy(files, ids, accessions, *rest, **kw):
            calls.append((sorted(files), set(ids), set(accessions), rest, kw))
            return real(files, ids, accessions, *rest, **kw)
        with mock.patch.object(stage2.stage1, "scan", side_effect=spy):
            made = self.run_stage2()
        self.assert_passes(made)
        self.assertEqual(len(calls), 1)
        files, ids, accessions, rest, kw = calls[0]
        self.assertEqual((files, rest, kw), (c.EXPORT, (), {}))
        every = {s["id"] for d in made.release["reports"] for s in d["samples"]}
        self.assertLessEqual(every, ids)
        self.assertEqual(accessions, {g["bioproject_accession"] for g in made.release["genomes"]
                                      if g["genome_id"] in {line.split("\t")[1] for line in
                                                            made.row_list("A").read_text().splitlines()}})

    def test_a_hit_names_the_file_and_its_count_and_writes_no_export(self):
        real = stage2.export_files

        def leaky(checked, result):
            out = real(checked, result)
            out["measure.tsv"] += b"gz-ka-01\n"
            return out
        with mock.patch.object(stage2, "export_files", side_effect=leaky):
            made = self.run_stage2()
        self.assertEqual(made.code, 1, made.lines)
        self.assertEqual(made.lines, [stage2.OPENED, "measure.tsv 1"])
        self.assertFalse(made.wrote_export())

    def test_a_decimal_point_in_an_export_file_stops_after_the_line(self):
        with mock.patch.object(stage2, "ppm", side_effect=lambda v: round(stage2.fractions.Fraction(v) * 10 ** 6) / 1):
            made = self.run_stage2()
        self.assert_stops_after_the_line(made, "measure.tsv holds a decimal point")

    def test_a_full_stop_in_the_cells_table_alone_stops_after_the_line(self):
        with mock.patch.object(stage2, "CALLS", {True: "resistant", False: "susceptible" + "."}):
            made = self.run_stage2()
        self.assert_stops_after_the_line(made, "cells.tsv holds a decimal point")


def left_out_and_unlisted(release):
    """A `before`: the designed release gains a report of a ninth species that holds none of the five learners, so
    stage 1 leaves it out, and a genome-table row of a made-up genome that no report holds."""
    rows = [{"id": made_up2.made_up.gid("ch", n), "resistant": n % 2 == 0, "fold": n % 5, "lineage": 30 + n,
             "mlst": "ST-1", "from_mic": False} for n in range(1, 8)]
    release["reports"].append(made_up2.made_up.report("synthia_chi__drug_a", "Synthia chi", "synthia_chi", "drug_a",
                                                        rows, None, learners=()))
    release["hc20"].append({"genome_id": made_up2.made_up.gid("hz", 1), "genome_name": "made-up", "mlst": "made-up",
                            "cgmlst_hc20": "hc-1"})


class Inputs(c.Case, unittest.TestCase):
    """inputs.tsv on the designed release with a report left out and a genome-table row no report holds: each row
    gives its own file's sha256 and count, and the scan reads the ids of the report left out and of the table."""

    def test_each_input_row_gives_its_own_files_sha256_and_count(self):
        made = self.run_stage2(before=left_out_and_unlisted)
        self.assert_passes(made)
        rows = {r["input"]: r for r in made.table("inputs.tsv")}
        reports = sorted((made.top / "reports").glob("*.json"))
        self.assertEqual(len(reports), 11)
        self.assertEqual((rows["reports"]["sha256"], rows["reports"]["count"]),
                         (c.sha((made.top / "reports.sha256").read_bytes()), "11"))
        self.assertTrue(rows["reports"]["note"].endswith("listed in a row list: 10"), rows["reports"]["note"])
        lines = {}
        for reading in stage2.READINGS:
            row, data = rows[f"row list, reading {reading}"], made.row_list(reading).read_bytes()
            lines[reading] = str(len(data.decode("utf-8").splitlines()))
            self.assertEqual((row["sha256"], row["count"]), (c.sha(data), lines[reading]), reading)
        self.assertNotEqual(rows["row list, reading A"]["sha256"], rows["row list, reading B"]["sha256"])
        self.assertEqual(lines["A"], "82")
        table = (made.top / "genomes.jsonl").read_bytes()
        held = len(table.decode("utf-8").splitlines())
        self.assertEqual(held, len(made.release["hc20"]))
        self.assertEqual((rows["genome table, 5 Nov"]["sha256"], rows["genome table, 5 Nov"]["count"]),
                         (c.sha(table), str(held)))
        self.assertTrue(rows["genome table, 5 Nov"]["note"].endswith("listed genomes it lacks: 1"),
                        rows["genome table, 5 Nov"]["note"])

    def test_an_id_of_the_report_left_out_or_of_the_genome_table_alone_is_a_hit(self):
        for planted in (made_up2.made_up.gid("ch", 1), made_up2.made_up.gid("hz", 1)):
            with self.subTest(planted=planted):
                real = stage2.export_files

                def leaky(checked, result, planted=planted):
                    out = real(checked, result)
                    out["measure.tsv"] += planted.encode("utf-8") + b"\n"
                    return out
                with mock.patch.object(stage2, "export_files", side_effect=leaky):
                    made = self.run_stage2(before=left_out_and_unlisted)
                listed = {line.split("\t")[1] for reading in stage2.READINGS
                          for line in made.row_list(reading).read_text(encoding="utf-8").splitlines()}
                self.assertNotIn(planted, listed)
                self.assertEqual(made.code, 1, made.lines)
                self.assertEqual(made.lines, [stage2.OPENED, "measure.tsv 1"])
                self.assertFalse(made.wrote_export())


class FailedWrite(c.Case, unittest.TestCase):
    def test_an_export_whose_third_file_cannot_be_written_leaves_no_export(self):
        real, tried = os.open, []

        def full(path, *a, **kw):
            p = Path(path)
            if (p.parent.name, p.parent.parent.name) == ("export", "stage2"):
                tried.append(p.name)
                if len(tried) == 3:
                    raise OSError("made-up full disk")
            return real(path, *a, **kw)
        with mock.patch.object(stage2.os, "open", side_effect=full):
            made = self.run_stage2()
        self.assertEqual(tried, ["measure.tsv", "interval.tsv", "cells.tsv"])
        self.assertEqual(made.code, 1, made.lines)
        self.assertEqual(made.lines, [stage2.OPENED, c.STOP + "the export could not be written (OSError); of the 3 "
                                                              "paths it made, 0 could not be removed"])
        self.assertFalse(made.wrote_export())
        self.assertEqual(sorted(p.name for p in made.out.iterdir()), [stage2.OPENED_FILE])


if __name__ == "__main__":
    unittest.main()
