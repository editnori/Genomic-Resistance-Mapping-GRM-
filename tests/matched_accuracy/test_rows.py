"""The rows of scripts/matched-accuracy/stage1.py: fields, lab records, BioProjects, reasons, groups, readings
and weights, each on a few SYNTHETIC rows written here (gz-<n> genome ids and pz-<n> BioProjects, made up).

The genome table tests need BioProjects of an accession's shape, PRJ, two letters and digits, since ruling R53
reads no other shape as a BioProject: acc() builds them from pieces, with the letters ZZ, which no registry uses."""

import collections
import fractions
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.matched_accuracy import common, made_up
from tests.matched_accuracy.common import outcome

s1 = common.stage1
F = fractions.Fraction
RESOLVE = {"Drug A": ("drug", "drug_a"), "Drug B": ("drug", "drug_b"), "Drug Class Q": ("class", None)}
ASSAY = ("Broth dilution", "CLSI", "2018")


def acc(n):
    """A made-up BioProject of an accession's shape (see the module's note)."""
    return "PRJ" + "ZZ" + f"9{n:03d}"


def write_lines(d, name, rows):
    path = Path(d) / name
    path.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")
    return path


def records(rows, wanted):
    """read_records on rows written to a new file, with the pinned line count set to theirs."""
    with tempfile.TemporaryDirectory() as d:
        path = write_lines(d, "records.jsonl", rows)
        with mock.patch.object(s1, "RECORDS_LINES", len(rows)):
            return outcome(s1.read_records, path, RESOLVE, wanted)


def genome_rows(rows, wanted):
    with tempfile.TemporaryDirectory() as d:
        return outcome(s1.read_genomes, write_lines(d, "genomes.jsonl", rows), wanted)


def row(dataset="d1", species="sp1", genome="gz-1", fold=0, lineage=1, bioproject="pz-1", resistant=True,
        intermediate=False, mlst="ST-1", assay=ASSAY):
    method, standard, year = assay
    return {"dataset": dataset, "species": species, "genome": genome, "fold": fold, "lineage": lineage,
            "bioproject": bioproject, "method": method, "standard": standard, "year": year, "resistant": resistant,
            "mlst": mlst, "intermediate": intermediate}


class Fields(unittest.TestCase):
    def test_value_text_trims_and_joins_as_grm_reads_a_field(self):
        self.assertEqual(s1.value_text(" Drug A "), "Drug A")
        self.assertIsNone(s1.value_text("   "))
        self.assertIsNone(s1.value_text(None))
        self.assertEqual(s1.value_text(["x", " ", " y "]), "x;y")
        self.assertEqual(s1.value_text(2018), "2018")

    def test_values_gives_each_text_of_a_list_apart(self):
        self.assertEqual(s1.values(["x", " y ", "", None]), {"x", "y"})
        self.assertEqual(s1.values(" x "), {"x"})
        self.assertEqual(s1.values(None), set())

    def test_parse_call_reads_each_word_as_grm_curate_does(self):
        for text, call in ((" Resistant ", "R"), ("Nonsusceptible", "R"), ("non-susceptible", "R"),
                           ("SUSCEPTIBLE", "S"), ("Intermediate", "I"), ("Susceptible-dose dependent", "I"),
                           ("Reduced susceptibility", "I"), ("", None), ("Not defined", None), (None, None)):
            self.assertEqual(s1.parse_call(text), call, text)

    def test_a_bioproject_value_splits_at_spaces_commas_and_semicolons(self):
        self.assertEqual(s1.bioproject_parts("pz-1"), {"pz-1"})
        self.assertEqual(s1.bioproject_parts("pz-1, pz-9"), {"pz-1", "pz-9"})
        self.assertEqual(s1.bioproject_parts(["pz-1;pz-2", " pz-3 "]), {"pz-1", "pz-2", "pz-3"})
        self.assertEqual(s1.bioproject_parts("  "), set())

    def test_an_insertion_time_is_its_first_19_characters(self):
        self.assertEqual(s1.stamp_heads("2026-11-05T19:09:59.999Z"), ["2026-11-05T19:09:59"])
        self.assertIsNone(s1.stamp_heads("2026-11-05"))
        self.assertIsNone(s1.stamp_heads(""))
        self.assertLess(max(s1.stamp_heads("2026-11-05T19:09:59Z")), s1.SEAL)
        self.assertGreaterEqual(max(s1.stamp_heads("2026-11-05T19:10:00+05:00")), s1.SEAL)

    def test_a_list_of_insertion_times_gives_each_text_apart(self):
        self.assertEqual(s1.stamp_heads(["2025-06-22T08:30:00Z", " 2026-11-05T19:10:00Z "]),
                         ["2025-06-22T08:30:00", "2026-11-05T19:10:00"])
        for v in ([], ["2025-06-22T08:30:00Z", ""], ["2025-06-22T08:30:00Z", None], ["2025-06-22T08:30:00Z", "x"],
                  ["2025-06-22T08:30:00Z", ["2025-06-22"]]):
            with self.subTest(v=v):
                self.assertIsNone(s1.stamp_heads(v))

    def test_an_impossible_date_or_time_is_no_insertion_time(self):
        for v in ("2025-13-01T00:00:00Z", "0000-00-00T00:00:00Z", "2025-02-30T08:30:00Z", "2025-06-22T24:00:00Z",
                  "2025-06-22t08:30:00Z", "\uff12\uff10\uff12\uff15-06-22T08:30:00Z"):
            with self.subTest(v=v):
                self.assertIsNone(s1.stamp_heads(v))


class Records(unittest.TestCase):
    WANT = {("gz-1", "drug_a"), ("gz-2", "drug_a"), ("gz-3", "drug_a"), ("gz-1", "drug_b")}

    def test_a_record_reaches_its_antibiotic_through_the_drug_names_file(self):
        got = records([made_up.record("gz-1", "Drug A"), made_up.record("gz-2", " Drug A "),
                       made_up.record("gz-1", "Drug B", call="Susceptible"),
                       made_up.record("gz-3", "Drug Class Q"), made_up.record("gz-9", "Drug A")], self.WANT)
        self.assertEqual(got[0], "ok")
        pairs, ids, n = got[1]
        self.assertEqual(set(pairs), {("gz-1", "drug_a"), ("gz-2", "drug_a"), ("gz-1", "drug_b")})
        self.assertEqual(ids, {"gz-1", "gz-2", "gz-3", "gz-9"})
        self.assertEqual(n, 5)

    def test_a_blank_method_or_standard_marks_the_pair(self):
        got = records([made_up.record("gz-1", "Drug A", standard=None), made_up.record("gz-2", "Drug A", method=" "),
                       made_up.record("gz-3", "Drug A")], self.WANT)
        pairs = got[1][0]
        self.assertEqual({k[0]: p["blank"] for k, p in pairs.items()}, {"gz-1": True, "gz-2": True, "gz-3": False})

    def test_a_standard_is_trimmed_and_uppercased_as_grm_curate_does(self):
        got = records([made_up.record("gz-1", "Drug A", standard="clsi"),
                       made_up.record("gz-1", "Drug A", method=" Broth dilution ", standard=" CLSI ")], self.WANT)
        self.assertEqual(got[1][0][("gz-1", "drug_a")]["tuples"], {("Broth dilution", "CLSI", "2018")})

    def test_two_years_of_one_method_and_standard_are_two_tuples(self):
        got = records([made_up.record("gz-1", "Drug A", year="2018"), made_up.record("gz-1", "Drug A", year="2020"),
                       made_up.record("gz-2", "Drug A", year=" ")], self.WANT)
        pairs = got[1][0]
        self.assertEqual(pairs[("gz-1", "drug_a")]["tuples"],
                         {("Broth dilution", "CLSI", "2018"), ("Broth dilution", "CLSI", "2020")})
        # A blank year is the level "unknown version".
        self.assertEqual(pairs[("gz-2", "drug_a")]["tuples"], {("Broth dilution", "CLSI", None)})

    def test_a_list_of_two_methods_in_one_record_gives_two_tuples(self):
        got = records([made_up.record("gz-1", "Drug A", method=["Broth dilution", " Disk diffusion "])], self.WANT)
        self.assertEqual(got[1][0][("gz-1", "drug_a")]["tuples"],
                         {("Broth dilution", "CLSI", "2018"), ("Disk diffusion", "CLSI", "2018")})

    def test_an_intermediate_call_by_curates_words_marks_the_pair(self):
        got = records([made_up.record("gz-1", "Drug A", call="Susceptible-dose dependent"),
                       made_up.record("gz-2", "Drug A", call=" intermediate "),
                       made_up.record("gz-3", "Drug A", call="Resistant")], self.WANT)
        pairs = got[1][0]
        self.assertEqual({k[0]: p["intermediate"] for k, p in pairs.items()}, {"gz-1": True, "gz-2": True, "gz-3": False})


class Labels(unittest.TestCase):
    def test_a_conflict_or_a_label_from_an_mic_in_the_table_is_no_call(self):
        rows = [{"genome_id": "gz-1", "antibiotic": "drug_a", "reported": "R", "conflict": False, "label": "R",
                 "label_from_mic": False},
                {"genome_id": "gz-2", "antibiotic": "drug_a", "reported": "S", "conflict": True, "label": None,
                 "label_from_mic": False},
                {"genome_id": "gz-3", "antibiotic": "drug_a", "reported": None, "conflict": False, "label": "R",
                 "label_from_mic": True},
                {"genome_id": "gz-4", "antibiotic": "drug_b", "reported": "S", "conflict": False, "label": "S",
                 "label_from_mic": False}]
        with tempfile.TemporaryDirectory() as d:
            got = outcome(s1.read_labels, write_lines(d, "labels.jsonl", rows),
                          {("gz-1", "drug_a"), ("gz-2", "drug_a"), ("gz-3", "drug_a")})
        self.assertEqual(got, ("ok", ({("gz-1", "drug_a"): "R", ("gz-2", "drug_a"): None, ("gz-3", "drug_a"): None},
                                      {"gz-1", "gz-2", "gz-3", "gz-4"})))

    def test_a_row_with_from_mic_true_is_counted_apart_and_not_checked(self):
        entering = [{"id": "d1", "antibiotic": "drug_a",
                     "samples": [{"id": "gz-1", "resistant": True, "from_mic": True},
                                 {"id": "gz-2", "resistant": False, "from_mic": False}]}]
        self.assertEqual(outcome(s1.check_labels, entering, {("gz-2", "drug_a"): "S"}), ("ok", (1, 1)))


class Genomes(unittest.TestCase):
    @staticmethod
    def read(rows, wanted):
        return genome_rows([{"genome_id": g, "bioproject_accession": b, "date_inserted": made_up.STAMP}
                            for g, b in rows], wanted)

    def test_a_genome_gets_one_bioproject_or_none(self):
        rows = [("gz-1", acc(1)), ("gz-2", f"{acc(1)}, {acc(9)}"), ("gz-3", ""), ("gz-4", acc(1)), ("gz-4", acc(2)),
                ("gz-5", acc(5)), ("gz-5", f" {acc(5)} "), ("gz-8", acc(8))]
        got = self.read(rows, {"gz-1", "gz-2", "gz-3", "gz-4", "gz-5", "gz-6"})
        self.assertEqual(got[0], "ok")
        project, why, ids, accessions, n = got[1]
        self.assertEqual(project, {"gz-1": acc(1), "gz-2": None, "gz-3": None, "gz-4": None, "gz-5": acc(5),
                                   "gz-6": None})
        self.assertEqual({g: w for g, w in why.items() if w}, {"gz-6": "not_in_table"})
        self.assertEqual(ids, {"gz-1", "gz-2", "gz-3", "gz-4", "gz-5", "gz-8"})
        self.assertEqual(accessions, {acc(1), acc(9), acc(2), acc(5), acc(8)})
        self.assertEqual(n, 8)

    def test_a_blank_row_beside_a_row_naming_a_bioproject_is_skipped(self):
        # As the study split's reader skips it (scripts/c2-cap/check_folds.py, bioprojects); two BioProjects
        # still give none, and a genome with blank rows alone is in the table but has none.
        rows = [("gz-1", ""), ("gz-1", acc(1)), ("gz-2", None), ("gz-2", acc(2)), ("gz-3", [" "]), ("gz-3", acc(3)),
                ("gz-4", acc(4)), ("gz-4", " "), ("gz-4", acc(5)), ("gz-5", ""), ("gz-5", None)]
        project, why, ids, accessions, n = self.read(rows, {"gz-1", "gz-2", "gz-3", "gz-4", "gz-5"})[1]
        self.assertEqual(project, {"gz-1": acc(1), "gz-2": acc(2), "gz-3": acc(3), "gz-4": None, "gz-5": None})
        self.assertEqual(why, {"gz-1": None, "gz-2": None, "gz-3": None, "gz-4": None, "gz-5": None})
        self.assertEqual(n, 11)

    def test_a_value_of_one_part_without_an_accessions_shape_reads_as_no_bioproject(self):
        # Ruling R53: a part is an accession when the whole part is PRJ, two letters and digits, without case.
        rows = [("gz-1", "-"), ("gz-2", "n/a"), ("gz-3", acc(3).lower()), ("gz-4", "PRJ" + "Z" + "9004"),
                ("gz-5", "PRJ" + "ZZZ" + "9005"), ("gz-6", "x" + acc(6)), ("gz-7", acc(7) + "x"),
                ("gz-8", f"-, {acc(8)}"), ("gz-9", "PRJ" + "ZZ"), ("gz-10", "PRJ" + "ZZ" + "\uff19\uff10"),
                ("gz-11", "-"), ("gz-11", "-")]
        wanted = {f"gz-{i}" for i in range(1, 12)}
        project, why, ids, accessions, n = self.read(rows, wanted)[1]
        self.assertEqual({g for g, b in project.items() if b}, {"gz-3"})
        self.assertEqual(project["gz-3"], acc(3).lower())
        self.assertEqual({g for g, w in why.items() if w == "no_accession"},
                         {"gz-1", "gz-2", "gz-4", "gz-5", "gz-6", "gz-7", "gz-9", "gz-10", "gz-11"})
        self.assertIsNone(why["gz-8"])
        # The scan looks up accession-shaped parts alone, those of a value of two parts too.
        self.assertEqual(accessions, {acc(3).lower(), acc(8)})


class HC20(unittest.TestCase):
    def test_a_genome_gets_one_hc20_group_or_none(self):
        # gz-1's two rows differ, so none; gz-2's two rows are one value once trimmed; gz-3's is blank; gz-4 has none.
        rows = [{"genome_id": "gz-1", "cgmlst_hc20": "hc-1"}, {"genome_id": "gz-1", "cgmlst_hc20": "hc-2"},
                {"genome_id": "gz-2", "cgmlst_hc20": "hc-2"}, {"genome_id": "gz-2", "cgmlst_hc20": " hc-2 "},
                {"genome_id": "gz-3", "cgmlst_hc20": " "}]
        with tempfile.TemporaryDirectory() as d:
            got = outcome(s1.read_hc20, write_lines(d, "hc20.jsonl", rows), {"gz-1", "gz-2", "gz-3", "gz-4"})
        self.assertEqual(got, ("ok", ({"gz-1": None, "gz-2": "hc-2", "gz-3": None, "gz-4": None},
                                      {"gz-1", "gz-2", "gz-3"}, 5, 1)))


class Reasons(unittest.TestCase):
    def report(self, n, copies=None):
        return {"id": "d1", "species_slug": "sp1", "antibiotic": "drug_a", "kset_duplicates": copies,
                "samples": [{"id": f"gz-{i}", "resistant": i % 2 == 0, "from_mic": False, "lineage": 1,
                             "mlst": "ST-1"} for i in range(n)],
                "lineage_split": {"test_fold": [0] * n}}

    @staticmethod
    def pair(*tuples, blank=False, intermediate=False):
        return {"tuples": set(tuples), "blank": blank, "intermediate": intermediate}

    def test_a_row_leaves_under_the_first_reason_that_applies(self):
        other = ("Disk diffusion", "CLSI", "2018")
        rep = self.report(7, copies=[[0, 1]])
        pairs = {("gz-1", "drug_a"): self.pair(ASSAY),
                 ("gz-4", "drug_a"): self.pair(ASSAY, other, blank=True),
                 ("gz-5", "drug_a"): self.pair(ASSAY, other),
                 ("gz-6", "drug_a"): self.pair(ASSAY)}
        project = {"gz-0": None, "gz-1": "pz-1", "gz-2": None, "gz-3": "pz-1", "gz-4": "pz-1", "gz-5": "pz-1",
                   "gz-6": "pz-1"}
        why = {g: "not_in_table" if g == "gz-0" else None for g in project}
        rows, out = s1.base_rows([rep], pairs, project, why)
        # gz-0: a copy group, no BioProject and no record; gz-1: a copy group; gz-2: no BioProject and no record;
        # gz-3: no record; gz-4: a blank and two tuples; gz-5: two tuples; gz-6 stays.
        self.assertEqual(out["d1"], collections.Counter({"copy_group": 2, "bioproject": 1, "no_lab_record": 1,
                                                         "blank_method_or_standard": 1, "assay_not_single": 1}))
        self.assertEqual([r["genome"] for r in rows], ["gz-6"])

    def test_a_pair_with_two_tuples_leaves_under_reason_5(self):
        rep = self.report(1)
        pairs = {("gz-0", "drug_a"): self.pair(ASSAY, ASSAY[:2] + ("2020",))}
        rows, out = s1.base_rows([rep], pairs, {"gz-0": "pz-1"}, {"gz-0": None})
        self.assertEqual((rows, dict(out["d1"])), ([], {"assay_not_single": 1}))

    def test_a_genome_the_second_copy_lacks_is_counted_apart_under_reason_2(self):
        rep = self.report(2)
        rows, out = s1.base_rows([rep], {}, {"gz-0": None, "gz-1": None}, {"gz-0": "not_in_table", "gz-1": None})
        self.assertEqual(dict(out["d1"]), {"bioproject": 2, "bioproject_not_in_table": 1})

    def test_a_genome_whose_bioproject_has_no_accessions_shape_is_counted_apart_under_reason_2(self):
        rep = self.report(3)
        why = {"gz-0": "no_accession", "gz-1": "not_in_table", "gz-2": None}
        rows, out = s1.base_rows([rep], {}, {"gz-0": None, "gz-1": None, "gz-2": None}, why)
        self.assertEqual(dict(out["d1"]), {"bioproject": 3, "bioproject_no_accession": 1,
                                           "bioproject_not_in_table": 1})

    def test_a_kept_row_carries_its_fields_and_a_blank_mlst_as_none(self):
        rep = self.report(1)
        rep["samples"][0]["mlst"] = " "
        rows, _ = s1.base_rows([rep], {("gz-0", "drug_a"): self.pair(ASSAY, intermediate=True)}, {"gz-0": "pz-1"},
                               {"gz-0": None})
        self.assertEqual(rows, [row(dataset="d1", species="sp1", genome="gz-0", fold=0, lineage=1, bioproject="pz-1",
                                    resistant=True, intermediate=True, mlst=None)])


class Groups(unittest.TestCase):
    def test_a_group_is_keyed_by_its_seven_fields(self):
        r = row(genome="gz-1", resistant=True)
        self.assertEqual(len(s1.groups_of([r, row(genome="gz-2", resistant=False)])[0]), 1)
        for field, value in (("dataset", "d2"), ("fold", 1), ("lineage", 2), ("bioproject", "pz-2"),
                             ("method", "Disk diffusion"), ("standard", "EUCAST"), ("year", None)):
            with self.subTest(field=field):
                s = dict(row(genome="gz-2", resistant=False), **{field: value})
                both, rest = s1.groups_of([r, s])
                self.assertEqual((both, len(rest)), ({}, 2))

    def test_a_group_with_one_lab_outcome_leaves(self):
        both, rest = s1.groups_of([row(genome="gz-1"), row(genome="gz-2")])
        self.assertEqual((both, [r["genome"] for r in rest]), ({}, ["gz-1", "gz-2"]))

    def test_the_pairs_of_a_group_are_its_resistant_rows_times_its_susceptible_rows(self):
        for r, s, pairs in ((1, 1, 1), (2, 2, 4), (1, 3, 3), (3, 2, 6)):
            group = [row(resistant=True)] * r + [row(resistant=False)] * s
            self.assertEqual(s1.pairs_in(group), pairs, (r, s))


class Readings(unittest.TestCase):
    def test_reading_b_takes_out_every_pair_with_an_intermediate_record_before_grouping(self):
        base = [row(genome="gz-1", resistant=True), row(genome="gz-2", resistant=False),
                row(genome="gz-3", resistant=False, lineage=2), row(genome="gz-4", resistant=True, lineage=2,
                                                                     intermediate=True),
                row(genome="gz-5", resistant=False, lineage=2, intermediate=True)]
        a_groups, a_out, a_rows = s1.reading(base, "A")
        self.assertEqual(sorted(len(g) for g in a_groups.values()), [2, 3])
        self.assertEqual(dict(a_out), {})
        b_groups, b_out, b_rows = s1.reading(base, "B")
        self.assertEqual([[r["genome"] for r in g] for g in b_groups.values()], [["gz-1", "gz-2"]])
        self.assertEqual({d: dict(c) for d, c in b_out.items()}, {"d1": {"intermediate": 2, "one_lab_outcome": 1}})
        self.assertEqual([r["genome"] for r in b_rows], ["gz-1", "gz-2", "gz-3"])


class Weights(unittest.TestCase):
    @staticmethod
    def groups(spec):
        """{key: rows} from (dataset, species, BioProject, rows) tuples, one group each, in fold order."""
        out = {}
        for fold, (dataset, species, bp, n) in enumerate(spec):
            rows = [row(dataset=dataset, species=species, fold=fold, bioproject=bp, genome=f"gz-{fold}-{i}",
                        resistant=i % 2 == 0) for i in range(n)]
            out[s1.key(rows[0])] = rows
        return out

    def test_groups_datasets_and_species_weigh_equally_and_exactly(self):
        g = self.groups([("d1", "x", "pz-1", 2), ("d1", "x", "pz-1", 2), ("d2", "x", "pz-2", 2),
                         ("d3", "y", "pz-3", 2)])
        w = s1.weights(g)
        self.assertEqual(sorted(w.values()), [F(1, 8), F(1, 8), F(1, 4), F(1, 2)])
        self.assertTrue(all(isinstance(v, F) for v in w.values()))
        self.assertEqual(sum(w.values()), 1)

    def test_the_largest_bioprojects_share_is_two_pairs_of_integers(self):
        # pz-1 is the largest by rows (6 of 10), with weight 1/4; pz-3 has the most weight, 1/2.
        g = self.groups([("d1", "x", "pz-1", 6), ("d1", "x", "pz-2", 2), ("d3", "y", "pz-3", 2)])
        self.assertEqual(s1.shares(g), {"bioprojects": 3, "eligible_rows": 10, "largest_bioproject_rows": 6,
                                        "its_weight_numerator": 1, "its_weight_denominator": 4,
                                        "largest_weight_numerator": 1, "largest_weight_denominator": 2})

    def test_among_bioprojects_tied_on_rows_the_share_gives_the_largest_weight(self):
        # pz-1 and pz-2 tie on 4 rows each; pz-1 weighs 1/2 (d1's one group), pz-2 1/4 (one of d2's two groups).
        got = s1.shares(self.groups([("d1", "x", "pz-1", 4), ("d2", "x", "pz-2", 4), ("d2", "x", "pz-4", 2)]))
        self.assertEqual((got["largest_bioproject_rows"], got["its_weight_numerator"], got["its_weight_denominator"]),
                         (4, 1, 2))

    def test_weights_that_do_not_add_up_to_1_stop(self):
        g = self.groups([("d1", "x", "pz-1", 2)])
        with mock.patch.object(s1, "weights", return_value={k: F(1, 2) for k in g}):
            self.assertEqual(outcome(s1.shares, g), ("stopped", "the analysis weights do not add up to 1"))


class Typing(unittest.TestCase):
    def test_a_narrow_match_is_counted_beside_the_primary_match_on_its_surviving_rows(self):
        rows = [row(genome="gz-1", resistant=True, mlst="ST-1"), row(genome="gz-2", resistant=False, mlst="ST-1"),
                row(genome="gz-3", resistant=True, mlst="ST-2"), row(genome="gz-4", resistant=False, mlst=None)]
        hc20 = {"gz-1": "hc-1", "gz-2": "hc-2", "gz-3": "hc-1", "gz-4": "hc-1"}
        names = {"dataset": {"d1"}, "species": {"sp1"}, "overall": {"all"}}
        groups = s1.groups_of(rows)[0]
        self.assertEqual(sum(len(g) for g in groups.values()), 4)
        got = s1.typing(rows, groups, hc20, names)
        mlst = got["mlst"]["dataset"]["d1"]
        self.assertEqual({k: mlst[k] for k in ("lacking_value_of_rows_entering_groups", "rows", "groups", "pairs",
                                               "primary_groups", "primary_pairs", "rows_lost_from_primary")},
                         {"lacking_value_of_rows_entering_groups": 1, "rows": 2, "groups": 1, "pairs": 1,
                          "primary_groups": 1, "primary_pairs": 1, "rows_lost_from_primary": 2})
        hc = got["mlst_hc20"]["overall"]["all"]
        self.assertEqual((hc["lacking_value_of_rows_entering_groups"], hc["rows"], hc["groups"], hc["mlst_groups"],
                          hc["rows_lost_from_primary"]), (1, 0, 0, 0, 4))

    def test_a_row_lacking_a_value_is_counted_among_the_rows_entering_groups_and_among_the_eligible_rows(self):
        # gz-1, gz-2 and gz-4 form a group with both outcomes; gz-3 sits alone and leaves by the one-outcome rule.
        # gz-3 and gz-4 lack an MLST type: 2 of the 4 rows that enter the groups, 1 of the 3 eligible rows.
        rows = [row(genome="gz-1", resistant=True, mlst="ST-1"), row(genome="gz-2", resistant=False, mlst="ST-1"),
                row(genome="gz-3", resistant=True, lineage=2, mlst=None),
                row(genome="gz-4", resistant=False, mlst=None)]
        names = {"dataset": {"d1"}, "species": {"sp1"}, "overall": {"all"}}
        got = s1.typing(rows, s1.groups_of(rows)[0], None, names)
        for level, name in (("dataset", "d1"), ("species", "sp1"), ("overall", "all")):
            c = got["mlst"][level][name]
            with self.subTest(level=level):
                self.assertEqual((c["lacking_value_of_rows_entering_groups"], c["lacking_value_of_eligible_rows"],
                                  c["rows"], c["rows_lost_from_primary"]), (2, 1, 2, 1))

    def test_without_hc20_only_the_mlst_step_is_counted(self):
        names = {"dataset": {"d1"}, "species": {"sp1"}, "overall": {"all"}}
        got = s1.typing([row()], {}, None, names)
        self.assertEqual(list(got), ["mlst"])


if __name__ == "__main__":
    unittest.main()
