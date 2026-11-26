"""The run on the SYNTHETIC release of tests/matched_accuracy/made_up.py, each count worked out by hand here.

The plan of 10 rows (made_up.py) gives groups at positions {1,2} {3,4} {5,6} {7,8,9} and {10}: each pair of
positions is one R and one S row; {7,8,9} is R, S, R (2 pairs); {10} is one S row, which leaves as reason 6.
So the plan alone gives 4 groups, 9 rows and 1 + 1 + 1 + 2 = 5 pairs. Each dataset changes it as noted below.
"""

import hashlib
import unittest

from tests.matched_accuracy import common, made_up

COLUMNS = ("out_copy_group", "out_bioproject", "out_bioproject_not_in_table", "out_no_lab_record",
           "out_blank_method_or_standard", "out_assay_not_single", "out_intermediate", "out_one_lab_outcome",
           "rows", "genomes", "groups", "pairs")
# (dataset, reading): the 12 columns above, worked out by hand.
HAND = {
    # alpha drug_a: positions 7 and 9 are a copy group (reason 1, 2 rows); 8 and 10 are left alone (reason 6).
    ("synthia_alpha__drug_a", "A"): (2, 0, 0, 0, 0, 0, 0, 2, 6, 6, 3, 3),
    ("synthia_alpha__drug_a", "B"): (2, 0, 0, 0, 0, 0, 0, 2, 6, 6, 3, 3),
    # alpha drug_b: positions 3 and 8 add an I record. A is the plan. B takes out 3 and 8, so {4}, {7,9} (R, R)
    # and {10} leave (4 rows); {1,2} and {5,6} stay: 2 groups, 4 rows, 2 pairs.
    ("synthia_alpha__drug_b", "A"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    ("synthia_alpha__drug_b", "B"): (0, 0, 0, 0, 0, 0, 2, 4, 4, 4, 2, 2),
    # alpha drug_c: the plan.
    ("synthia_alpha__drug_c", "A"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    ("synthia_alpha__drug_c", "B"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    # beta drug_a: genome 1 has a blank BioProject and genome 5 two (reason 2). {2}, {6} and {10} leave;
    # {3,4} and {7,8,9} stay: 2 groups, 5 rows, 1 + 2 = 3 pairs.
    ("synthia_beta__drug_a", "A"): (0, 2, 0, 0, 0, 0, 0, 3, 5, 5, 2, 3),
    ("synthia_beta__drug_a", "B"): (0, 2, 0, 0, 0, 0, 0, 3, 5, 5, 2, 3),
    # beta drug_b: genome 11 (position 1) is not in the second copy (reason 2, counted apart). {2} and {10} leave;
    # {3,4} {5,6} {7,8,9} stay: 3 groups, 7 rows, 4 pairs.
    ("synthia_beta__drug_b", "A"): (0, 1, 1, 0, 0, 0, 0, 2, 7, 7, 3, 4),
    ("synthia_beta__drug_b", "B"): (0, 1, 1, 0, 0, 0, 0, 2, 7, 7, 3, 4),
    # beta drug_c: as beta drug_b; position 9 has from_mic true and stays eligible.
    ("synthia_beta__drug_c", "A"): (0, 1, 1, 0, 0, 0, 0, 2, 7, 7, 3, 4),
    ("synthia_beta__drug_c", "B"): (0, 1, 1, 0, 0, 0, 0, 2, 7, 7, 3, 4),
    # gamma drug_a: position 3 has two methods (reason 5). Position 5's "clsi" and "CLSI " are one tuple.
    # Positions 7 and 8 have a blank year, 9 has 2018: {7,8} is a group, {9} leaves. {4} and {10} leave too.
    ("synthia_gamma__drug_a", "A"): (0, 0, 0, 0, 0, 1, 0, 3, 6, 6, 3, 3),
    ("synthia_gamma__drug_a", "B"): (0, 0, 0, 0, 0, 1, 0, 3, 6, 6, 3, 3),
    # gamma drug_b: position 1 has no record (reason 3), 3 a blank method and 6 a blank standard (reason 4).
    # Position 2's two equal records are one tuple. {2} {4} {5} {10} leave; {7,8,9} stays: 1 group, 3 rows, 2 pairs.
    ("synthia_gamma__drug_b", "A"): (0, 0, 0, 1, 2, 0, 0, 4, 3, 3, 1, 2),
    ("synthia_gamma__drug_b", "B"): (0, 0, 0, 1, 2, 0, 0, 4, 3, 3, 1, 2),
    # gamma drug_c: position 1 adds an I record. B takes it out, and {2} leaves: 3 groups, 7 rows, 4 pairs.
    ("synthia_gamma__drug_c", "A"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    ("synthia_gamma__drug_c", "B"): (0, 0, 0, 0, 0, 0, 1, 2, 7, 7, 3, 4),
    # delta drug_a: genomes 1 to 6 are pz-de-1 and 7 to 10 pz-de-2; no group crosses that line, so the plan.
    ("synthia_delta__drug_a", "A"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    ("synthia_delta__drug_a", "B"): (0, 0, 0, 0, 0, 0, 0, 1, 9, 9, 4, 5),
    # delta drug_b: 12 rows; genome 21 (S) joins {7,8,9} and genome 22 (R) joins {10}:
    # {1,2} {3,4} {5,6} give 3 pairs, {7,8,9,21} 2 x 2 = 4 and {10,22} 1: 5 groups, 12 rows, 8 pairs.
    ("synthia_delta__drug_b", "A"): (0, 0, 0, 0, 0, 0, 0, 0, 12, 12, 5, 8),
    ("synthia_delta__drug_b", "B"): (0, 0, 0, 0, 0, 0, 0, 0, 12, 12, 5, 8),
}
# Species and overall: the sums of the rows above; genomes are distinct genomes, so alpha's three datasets on
# genomes 1 to 10 give 9 (A: 1 to 9) and beta's give 5 + 7 = 12.
HAND_UP = {
    ("species", "synthia_alpha", "A"): (2, 0, 0, 0, 0, 0, 0, 4, 24, 9, 11, 13),
    ("species", "synthia_alpha", "B"): (2, 0, 0, 0, 0, 0, 2, 7, 19, 9, 9, 10),
    ("species", "synthia_beta", "A"): (0, 4, 2, 0, 0, 0, 0, 7, 19, 12, 8, 11),
    ("species", "synthia_beta", "B"): (0, 4, 2, 0, 0, 0, 0, 7, 19, 12, 8, 11),
    ("species", "synthia_gamma", "A"): (0, 0, 0, 1, 2, 1, 0, 8, 18, 18, 8, 10),
    ("species", "synthia_gamma", "B"): (0, 0, 0, 1, 2, 1, 1, 9, 16, 16, 7, 9),
    ("species", "synthia_delta", "A"): (0, 0, 0, 0, 0, 0, 0, 1, 21, 21, 9, 13),
    ("species", "synthia_delta", "B"): (0, 0, 0, 0, 0, 0, 0, 1, 21, 21, 9, 13),
    # 120 genome rows; 8 of them in the report left out; 112 - 10 (reasons 1 to 5) = 102; A 102 - 20 = 82;
    # B 102 - 3 - 24 = 75.
    ("overall", "all", "A"): (2, 4, 2, 1, 2, 1, 0, 20, 82, 60, 36, 47),
    ("overall", "all", "B"): (2, 4, 2, 1, 2, 1, 3, 24, 75, 58, 33, 43),
}
# Typing, delta (the other species have one MLST type and one HC20 group, so their typing equals their counts).
# delta drug_a: MLST 1, 2 ST-1; 3 ST-3; 4 ST-4; 5 ST-5; 6 none; 7, 8 ST-7; 9 ST-9; 10 ST-10. HC20 1, 2 hc-1;
# 3 to 6 hc-3; 7 hc-7; 8 hc-8; 9 hc-9; 10 not in the 5 Nov table.
#   mlst: 6 lacks a value; {1,2} and {7,8} stay (4 rows, 2 pairs); primary on 1, 2, 7, 8: 2 groups, 2 pairs.
#   mlst_hc20: 6 and 10 lack a value; {1,2} alone stays; 7 and 8 differ in HC20.
# delta drug_b: MLST ST-2, but position 9 and genome 21 are ST-6; HC20 hc-2, but position 4 is blank.
#   mlst: {1,2} {3,4} {5,6} {7,8} {9,21} {10,22}: 6 groups, 12 rows, 6 pairs; primary 5 groups, 8 pairs.
#   mlst_hc20: 4 lacks a value, so {3} leaves: 5 groups, 10 rows, 5 pairs; primary on those 10 rows: {1,2} {5,6}
#   {7,8,9,21} {10,22}, 4 groups, 1 + 1 + 4 + 1 = 7 pairs; MLST on them 5 groups, 5 pairs.
TYPING_COLUMNS = ("rows_lacking_value", "rows", "genomes", "groups", "pairs", "primary_groups", "primary_pairs",
                  "mlst_groups", "mlst_pairs", "rows_lost_from_primary")
HAND_TYPING = {
    ("dataset", "synthia_delta__drug_a", "mlst"): (1, 4, 4, 2, 2, 2, 2, "", "", 5),
    ("dataset", "synthia_delta__drug_a", "mlst_hc20"): (2, 2, 2, 1, 1, 1, 1, 1, 1, 7),
    ("dataset", "synthia_delta__drug_b", "mlst"): (0, 12, 12, 6, 6, 5, 8, "", "", 0),
    ("dataset", "synthia_delta__drug_b", "mlst_hc20"): (1, 10, 10, 5, 5, 4, 7, 5, 5, 2),
    ("species", "synthia_delta", "mlst"): (1, 16, 16, 8, 8, 7, 10, "", "", 5),
    ("species", "synthia_delta", "mlst_hc20"): (3, 12, 12, 6, 6, 5, 8, 6, 6, 9),
    # Overall: alpha, beta and gamma add their counts (A: 61 rows, 39 genomes, 27 groups, 34 pairs;
    # B: 54 rows, 37 genomes, 24 groups, 30 pairs) to delta's.
    ("overall", "all", "mlst", "A"): (1, 77, 55, 35, 42, 34, 44, "", "", 5),
    ("overall", "all", "mlst_hc20", "A"): (3, 73, 51, 33, 40, 32, 42, 33, 40, 9),
    ("overall", "all", "mlst", "B"): (1, 70, 53, 32, 38, 31, 40, "", "", 5),
    ("overall", "all", "mlst_hc20", "B"): (3, 66, 49, 30, 36, 29, 38, 30, 36, 9),
}
# The eligible positions of each dataset by hand, for the row lists. Positions 11 and 12 of delta drug_b are
# genomes 21 and 22.
ELIGIBLE = {
    ("synthia_alpha__drug_a", "A"): (1, 2, 3, 4, 5, 6), ("synthia_alpha__drug_a", "B"): (1, 2, 3, 4, 5, 6),
    ("synthia_alpha__drug_b", "A"): tuple(range(1, 10)), ("synthia_alpha__drug_b", "B"): (1, 2, 5, 6),
    ("synthia_alpha__drug_c", "A"): tuple(range(1, 10)), ("synthia_alpha__drug_c", "B"): tuple(range(1, 10)),
    ("synthia_beta__drug_a", "A"): (3, 4, 7, 8, 9), ("synthia_beta__drug_a", "B"): (3, 4, 7, 8, 9),
    ("synthia_beta__drug_b", "A"): tuple(range(3, 10)), ("synthia_beta__drug_b", "B"): tuple(range(3, 10)),
    ("synthia_beta__drug_c", "A"): tuple(range(3, 10)), ("synthia_beta__drug_c", "B"): tuple(range(3, 10)),
    ("synthia_gamma__drug_a", "A"): (1, 2, 5, 6, 7, 8), ("synthia_gamma__drug_a", "B"): (1, 2, 5, 6, 7, 8),
    ("synthia_gamma__drug_b", "A"): (7, 8, 9), ("synthia_gamma__drug_b", "B"): (7, 8, 9),
    ("synthia_gamma__drug_c", "A"): tuple(range(1, 10)), ("synthia_gamma__drug_c", "B"): tuple(range(3, 10)),
    ("synthia_delta__drug_a", "A"): tuple(range(1, 10)), ("synthia_delta__drug_a", "B"): tuple(range(1, 10)),
    ("synthia_delta__drug_b", "A"): tuple(range(1, 13)), ("synthia_delta__drug_b", "B"): tuple(range(1, 13)),
}
FIRST = {"alpha": 1, "beta_drug_a": 1, "beta": 11, "gamma_drug_a": 1, "gamma_drug_b": 11, "gamma_drug_c": 21,
         "delta_drug_a": 1, "delta_drug_b": 11}


def expected_line(dataset, p):
    """The row list line of position p of a dataset, from the plan and each dataset's note."""
    species, drug = dataset[len("synthia_"):].split("__")
    first = FIRST.get(f"{species}_{drug}", FIRST.get(species, 1))
    tag = species[:2]
    fold, lineage = ((3, 13), (4, 14))[p - 11] if p > 10 else made_up.PLAN[p - 1][:2]
    genome = made_up.gid(tag, first + p - 1)
    n = first + p - 1
    bioproject = {"alpha": "pz-al-1", "beta": "pz-be-1" if n <= 10 else "pz-be-2",
                  "gamma": f"pz-ga-{(n - 1) // 10 + 1}", "delta": "pz-de-1" if n <= 6 else "pz-de-2"}[species]
    year = "" if (dataset, p) in (("synthia_gamma__drug_a", 7), ("synthia_gamma__drug_a", 8)) else "2018"
    return "\t".join([dataset, genome, str(fold), str(lineage), bioproject, "Broth dilution", "CLSI", year]) + "\n"


class MadeUp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.made = common.Run()

    @classmethod
    def tearDownClass(cls):
        cls.made.close()

    def ran(self):
        self.assertEqual(self.made.code, 0, self.made.lines)

    def test_the_release_has_12_reports_in_4_species_with_each_case_the_proof_names(self):
        rel = common.release()
        self.assertEqual(len(rel["reports"]), 12)
        self.assertEqual(len({r["species"] for r in rel["reports"]}), 4)
        self.assertEqual([r["id"] for r in rel["reports"] if not r["splits"][0]["learners"]], ["synthia_delta__drug_c"])
        self.assertEqual([r["id"] for r in rel["reports"] if r.get("kset_duplicates")], ["synthia_alpha__drug_a"])
        self.assertEqual(len([x for x in rel["records"] if x["resistant_phenotype"] == "Intermediate"]), 3)
        self.assertEqual([g["genome_id"] for g in rel["genomes"] if g["bioproject_accession"] == ""], ["gz-be-01"])
        self.assertEqual(len(common.records_of(rel, "gz-ga-03", "Drug A")), 2)

    def test_the_run_prints_its_counts_in_four_lines(self):
        self.ran()
        sums = {n: hashlib.sha256(self.made.row_list(n)).hexdigest() for n in ("A", "B")}
        self.assertEqual(self.made.lines, [
            "stage 1: 11 of 12 reports enter, 1 left out; 127 lab records, each inserted before 2026-11-05T19:10:00",
            f"reading A: 82 eligible rows, 60 genomes, 36 groups, 47 pairs, 11 datasets, 4 species; row list sha256 "
            f"{sums['A']}",
            f"reading B: 75 eligible rows, 58 genomes, 33 groups, 43 pairs, 11 datasets, 4 species; row list sha256 "
            f"{sums['B']}",
            "scores opened: no"])

    def test_the_counts_of_each_dataset_equal_the_hand_count(self):
        self.ran()
        for (dataset, reading), want in HAND.items():
            with self.subTest(dataset=dataset, reading=reading):
                got = self.made.counts("dataset", dataset, reading)
                self.assertEqual(tuple(int(got[c]) for c in COLUMNS), want)
                rows = 12 if dataset == "synthia_delta__drug_b" else 10
                self.assertEqual((got["reports"], got["reports_left_out"], got["genome_rows"], got["datasets"],
                                  got["species"]), ("1", "0", str(rows), "1", "1"))
                self.assertEqual(got["from_mic_rows"], "1" if dataset == "synthia_beta__drug_c" else "0")

    def test_the_counts_of_each_species_and_overall_equal_the_hand_count(self):
        self.ran()
        for (level, name, reading), want in HAND_UP.items():
            with self.subTest(name=name, reading=reading):
                got = self.made.counts(level, name, reading)
                self.assertEqual(tuple(int(got[c]) for c in COLUMNS), want)
        a = self.made.counts("overall", "all", "A")
        self.assertEqual((a["reports"], a["reports_left_out"], a["genome_rows"], a["from_mic_rows"], a["datasets"],
                          a["species"]), ("12", "1", "120", "1", "11", "4"))

    def test_a_report_without_learners_is_counted_and_left_out(self):
        self.ran()
        for reading in ("A", "B"):
            got = self.made.counts("dataset", "synthia_delta__drug_c", reading)
            self.assertEqual({c: got[c] for c in ("reports", "reports_left_out", "genome_rows", "rows", "groups")},
                             {"reports": "1", "reports_left_out": "1", "genome_rows": "8", "rows": "0", "groups": "0"})
        delta = self.made.counts("species", "synthia_delta", "A")
        self.assertEqual((delta["reports"], delta["reports_left_out"], delta["genome_rows"], delta["datasets"]),
                         ("3", "1", "30", "2"))

    def test_the_shares_equal_the_hand_count(self):
        # A: 8 BioProjects (pz-al-1, pz-be-1, pz-be-2, pz-ga-1 to 3, pz-de-1, pz-de-2); pz-al-1 has the most rows,
        # 6 + 9 + 9 = 24 of 82, and all of alpha's weight, 1/4. pz-de-2 has 3 + 12 = 15 rows and 1/32 + 1/8 = 5/32.
        # B: pz-al-1 has 6 + 4 + 9 = 19 of 75 rows and still 1/4.
        self.ran()
        self.assertEqual(self.made.table("shares.tsv"), [
            dict(zip(common.stage1.SHARE_COLUMNS, ("A", "8", "82", "24", "1", "4", "1", "4"))),
            dict(zip(common.stage1.SHARE_COLUMNS, ("B", "8", "75", "19", "1", "4", "1", "4")))])

    def test_the_typing_counts_equal_the_hand_count(self):
        self.ran()
        for k, want in HAND_TYPING.items():
            level, name, step = k[:3]
            for reading in (k[3:] or ("A", "B")):
                with self.subTest(name=name, step=step, reading=reading):
                    got = self.made.typing(level, name, reading, step)
                    self.assertEqual(tuple(got[c] for c in TYPING_COLUMNS), tuple(str(x) for x in want))

    def test_the_typing_of_one_type_equals_the_counts(self):
        self.ran()
        for (dataset, reading), want in HAND.items():
            if "delta" in dataset:
                continue
            for step in ("mlst", "mlst_hc20"):
                with self.subTest(dataset=dataset, reading=reading, step=step):
                    got = self.made.typing("dataset", dataset, reading, step)
                    self.assertEqual((got["rows_lacking_value"], got["rows"], got["groups"], got["pairs"],
                                      got["primary_groups"], got["primary_pairs"], got["rows_lost_from_primary"]),
                                     ("0", str(want[8]), str(want[10]), str(want[11]), str(want[10]), str(want[11]),
                                      "0"))

    def test_the_row_lists_hold_the_hand_counted_rows(self):
        self.ran()
        for reading in ("A", "B"):
            lines = sorted(expected_line(d, p) for (d, r), ps in ELIGIBLE.items() if r == reading for p in ps)
            with self.subTest(reading=reading):
                self.assertEqual(len(lines), {"A": 82, "B": 75}[reading])
                self.assertEqual(self.made.row_list(reading).decode("utf-8"), "".join(lines))

    def test_the_inputs_and_checks_equal_the_hand_count(self):
        # Lab records: one a genome row of the 12 reports (120), plus alpha drug_b's 2 I records, gamma drug_a's
        # 2 more, gamma drug_b's 1 less and 1 more, gamma drug_c's I record, and one record of a drug class and one
        # of a genome in no report: 120 + 2 + 2 + 0 + 1 + 2 = 127.
        # Genome table, second copy: 10 + 19 + 30 + 30 + 1 = 90 rows (beta's genome 11 is not in it).
        self.ran()
        inputs = {r["input"]: (r["count"], r["note"]) for r in self.made.table("inputs.tsv")}
        self.assertEqual({k: v[0] for k, v in inputs.items()},
                         {"reports": "12", "lab records": "127", "drug names": "4", "label table": "111",
                          "genome table, second copy": "90", "genome table, 5 Nov": "90", "row list, reading A": "82",
                          "row list, reading B": "75"})
        self.assertTrue(inputs["genome table, second copy"][1].endswith("genomes of the reports that enter not in it: 1"))
        self.assertTrue(inputs["genome table, 5 Nov"][1].endswith("genomes of the reports that enter not in it: 1"))
        checks = dict((r["check"], r["result"]) for r in self.made.table("checks.tsv"))
        self.assertEqual(checks["reports"], "11 of 12 reports hold the five learners under the lineage split; "
                                            "left out, holding none of them: 1")
        self.assertEqual(checks["lab calls"], "rows with from_mic false, each equal to the label table: 111; rows with "
                                              "from_mic true, counted apart: 1")
        # 90 genomes in the reports and gz-zz-01 in the records: 91 ids; pz-be-9 and pz-zz-1 with the 8 above: 10.
        self.assertTrue(checks["scan"].startswith("91 genome ids and 10 BioProject accessions read by the run"))
        self.assertEqual(checks["scores opened"], "no")


if __name__ == "__main__":
    unittest.main()
