"""scripts/c3-inputs/m_list.py: the list of genomes with an intermediate call and the floor, from the label table.

Every data root here is SYNTHETIC: a made-up table, records, drug names and
matrix, with ids of the form a genome id takes. No record or call is real.
"""

import contextlib
import hashlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from tests.c3_inputs import common, m_common

m_list = common.load("m_list")
D = m_common.M_SET
SP = m_common.SPECIES
DRUG = "drug m"


def root(top, r_only=30, s_only=30, change_table=None, drop_record=False, lists_change=None):
    """A data root with r_only R genomes, s_only S genomes and seven of each other kind:
    S with I (listed, S), I with R (listed, R), nonsusceptible with intermediate (listed, R), I alone (no label),
    R with S (a conflict), one outside the matrix and one labelled from an MIC."""
    top = Path(top)
    gids = m_common.ids(r_only + s_only + 7)
    kinds = ["R"] * r_only + ["S"] * s_only + ["S+I", "I+R", "N+I", "I", "R+S", "out", "mic"]
    words = {"R": ["Resistant"], "S": ["Susceptible"], "I": ["Intermediate"], "N": ["Nonsusceptible"]}
    records, table = [], []
    for g, k in zip(gids, kinds):
        calls = {"S+I": ["S", "I"], "I+R": ["I", "R"], "N+I": ["N", "I"], "R+S": ["R", "S"], "out": ["R"],
                 "mic": ["R"]}.get(k, [k])
        for c in calls:
            for w in words[c]:
                records.append({"genome_id": g, "antibiotic": " Drug M ", "resistant_phenotype": w})
        reported, conflict = {"S+I": ("S", False), "I+R": ("R", False), "N+I": ("R", False), "I": ("I", False),
                              "R+S": (None, True), "out": ("R", False), "mic": ("R", False)}.get(k, (k, False))
        table.append({"genome_id": g, "species": "Synthia coli", "antibiotic": DRUG, "reported": reported,
                      "conflict": conflict, "label_from_mic": k == "mic"})
    # A record of a drug class, which names no antibiotic, and one of another drug.
    records.append({"genome_id": gids[0], "antibiotic": "Aminoglycosides", "resistant_phenotype": "Susceptible"})
    records.append({"genome_id": gids[1], "antibiotic": "Drug Z", "resistant_phenotype": "Intermediate"})
    table.append({"genome_id": gids[1], "species": "Synthia coli", "antibiotic": "drug z", "reported": "I",
                  "conflict": False, "label_from_mic": False})
    if drop_record:
        records = records[1:]
    if change_table:
        change_table(table, gids)
    counts = {}
    for r in records:
        raw = r["antibiotic"].strip()
        counts[raw] = counts.get(raw, 0) + 1
    if drop_record:
        counts["Drug M"] += 1
    names = [{"raw": "Drug M", "kind": "drug", "canonical": DRUG, "records": counts["Drug M"], "class": None},
             {"raw": "Aminoglycosides", "kind": "class", "canonical": None, "records": 1, "class": "aminoglycoside"},
             {"raw": "Drug Z", "kind": "drug", "canonical": "drug z", "records": 1, "class": None}]
    (top / "curated").mkdir(parents=True)
    (top / "curated" / "drug_names.json").write_text(json.dumps({"vocabulary": "synthetic", "names": names}))
    (top / "curated" / "datasets.json").write_text(json.dumps([{"id": D, "species": "Synthia coli", "species_slug": SP,
                                                                 "antibiotic": DRUG}]))
    with open(top / "curated" / "phenotypes.jsonl", "w") as f:
        for t in table:
            f.write(json.dumps(t) + "\n")
    with open(top / "records.jsonl", "w") as f:
        for r in records:
            f.write(json.dumps(r) + "\n")
    matrix = [g for g, k in zip(gids, kinds) if k != "out"]
    (top / "species" / SP / "matrix").mkdir(parents=True)
    (top / "species" / SP / "matrix" / "matrix.json").write_text(json.dumps({"genomes": matrix}))
    lrow = {"dataset": D, "species_slug": SP, "species": "Synthia coli", "antibiotic": DRUG, "in_M": "1",
            "samples": str(r_only + s_only + 3), "resistant": str(r_only + 2), "susceptible": str(s_only + 1)}
    if lists_change:
        lists_change(lrow)
    return top, top / "records.jsonl", [lrow], gids


class Calls(unittest.TestCase):
    def test_the_label_rule_r_or_s_wins(self):
        self.assertEqual(m_list.consolidate({"S", "I"}), ("S", False))
        self.assertEqual(m_list.consolidate({"I", "R"}), ("R", False))
        self.assertEqual(m_list.consolidate({"I"}), ("I", False))
        self.assertEqual(m_list.consolidate({"R", "S"}), (None, True))
        self.assertEqual(m_list.consolidate({"R", "S", "I"}), (None, True))

    def test_the_call_words_are_parse_calls(self):
        for word, call in (("Resistant", "R"), (" non-susceptible ", "R"), ("Nonsusceptible", "R"),
                           ("Susceptible", "S"), ("Intermediate", "I"), ("Susceptible-dose dependent", "I"),
                           ("Reduced susceptibility", "I"), ("Not defined", None), (None, None)):
            self.assertEqual(m_list.parse_call(word), call, word)

    def test_a_record_field_is_read_as_grm_reads_it(self):
        self.assertEqual(m_list.value_text(["Resistant", " ", "Intermediate"]), "Resistant;Intermediate")
        self.assertEqual(m_list.value_text("  "), None)
        self.assertEqual(m_list.value_text(573), "573")


class List(unittest.TestCase):
    def test_the_list_holds_s_with_i_and_i_with_r_and_not_i_alone_or_r_with_s(self):
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, gids = root(d)
            rows, lists, lines = m_list.build(top, records, lrows, SP, [D])
            n = 60
            self.assertEqual(lists[D], [gids[n], gids[n + 1], gids[n + 2]])
            r = rows[0]
            self.assertEqual((r["cohort"], r["resistant"], r["susceptible"]), (63, 32, 31))
            self.assertEqual((r["listed"], r["listed_resistant"], r["listed_susceptible"]), (3, 2, 1))
            self.assertEqual((r["second_resistant"], r["second_susceptible"], r["enters"]), (30, 30, 1))
            self.assertTrue(lines[0].startswith("records: every pair"))
            self.assertFalse(any(g in "".join(lines) for g in gids))

    def test_records_that_consolidate_to_another_call_stop(self):
        def change(table, gids):
            table[60]["reported"] = "R"
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, _ = root(d, change_table=change)
            with self.assertRaisesRegex(m_list.Stop, "do not consolidate"):
                m_list.build(top, records, lrows, SP, [D])

    def test_records_whose_name_counts_differ_from_drug_names_stop(self):
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, _ = root(d, drop_record=True)
            with self.assertRaisesRegex(m_list.Stop, "name counts differ"):
                m_list.build(top, records, lrows, SP, [D])

    def test_a_cohort_unlike_lists_tsv_stops(self):
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, _ = root(d, lists_change=lambda r: r.update(resistant=str(int(r["resistant"]) + 1)))
            with self.assertRaisesRegex(m_list.Stop, "lists.tsv gives"):
                m_list.build(top, records, lrows, SP, [D])

    def test_a_dataset_outside_m_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, _ = root(d, lists_change=lambda r: r.update(in_M="0"))
            with self.assertRaisesRegex(m_list.Stop, "not a dataset of M"):
                m_list.build(top, records, lrows, SP, [D])


class Floor(unittest.TestCase):
    def test_twenty_five_of_each_class_enter_and_twenty_four_do_not(self):
        for r_only, s_only, enters in ((25, 25, 1), (24, 30, 0), (30, 24, 0)):
            with tempfile.TemporaryDirectory() as d:
                top, records, lrows, _ = root(d, r_only=r_only, s_only=s_only)
                row = m_list.build(top, records, lrows, SP, [D])[0][0]
                self.assertEqual((row["second_resistant"], row["second_susceptible"]), (r_only, s_only))
                self.assertEqual(row["enters"], enters, (r_only, s_only))

    def test_the_main_form_writes_the_list_and_a_floor_table_of_counts(self):
        with tempfile.TemporaryDirectory() as d:
            top, records, lrows, gids = root(d)
            lists_path = Path(d) / "lists.tsv"
            m_common.write_tsv(lists_path, list(lrows[0]), lrows)
            with contextlib.redirect_stdout(io.StringIO()) as out:
                rc = m_list.main([str(top), str(records), str(lists_path), SP, D, str(Path(d) / "list"),
                                  str(Path(d) / "m" / "floor.tsv")])
            self.assertEqual(rc, 0)
            self.assertEqual((Path(d) / "list" / f"{D}.txt").read_text().split(), gids[60:63])
            floor = (Path(d) / "m" / "floor.tsv").read_text()
            self.assertEqual(floor.splitlines()[0].split("\t"), m_list.FLOOR_COLUMNS)
            self.assertFalse(any(g in floor or g in out.getvalue() for g in gids))
            self.assertEqual(out.getvalue().splitlines()[0], m_list.records_line(records))


class RecordsLine(unittest.TestCase):
    """The records file's sha256 is printed and flushed before any check, so a stop still shows it."""

    def lay_out(self, d, bad):
        top, records, lrows, gids = root(d)
        with open(records, "a") as f:
            f.write(bad)
        lists_path = Path(d) / "lists.tsv"
        m_common.write_tsv(lists_path, list(lrows[0]), lrows)
        argv = [str(top), str(records), str(lists_path), SP, D, str(Path(d) / "list"), str(Path(d) / "m" / "floor.tsv")]
        want = f"the calls are read from records.jsonl, sha256 {hashlib.sha256(records.read_bytes()).hexdigest()}"
        return argv, want, gids

    def test_a_stop_on_a_bad_record_still_prints_the_records_sha256_first(self):
        with tempfile.TemporaryDirectory() as d:
            argv, want, gids = self.lay_out(d, json.dumps({"genome_id": "90007.99999", "antibiotic": " "}) + "\n")
            with contextlib.redirect_stdout(io.StringIO()) as out:
                self.assertEqual(m_list.main(argv), 1)
            lines = out.getvalue().splitlines()
            self.assertTrue(lines[0].endswith(want), lines)
            self.assertEqual(lines[1:], ["m_list: stopped: a record has no antibiotic name"])
            self.assertFalse((Path(d) / "m" / "floor.tsv").exists())

    def test_a_record_that_is_not_json_ends_the_run_after_the_line_is_out(self):
        with tempfile.TemporaryDirectory() as d:
            argv, want, _ = self.lay_out(d, "{not json\n")
            done = subprocess.run([sys.executable, "-I", "-B", str(common.KIT / "m_list.py")] + argv,
                                  capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertTrue(done.stdout.splitlines()[0].endswith(want), done.stdout)
            self.assertIn("JSONDecodeError", done.stderr)


if __name__ == "__main__":
    unittest.main()
