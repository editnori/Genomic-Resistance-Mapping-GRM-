"""scripts/mechanism/catalogue_size.py, the catalogue size reading (docs/catalogue-size-result.md).

Every report, genome id, antibiotic, species file and count here is SYNTHETIC:
the reports, the lineage curve's endpoints and the species files are built in
the test. Only test_the_page_states_the_rule_the_script_applies reads a
committed file, the page's section "Rule".
"""

import contextlib
import hashlib
import io
import json
import math
import re
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from tests.mechanism import common as c

cs = c.load("catalogue_size")
ENDPOINT_COLUMNS = c.load("curve_tables").ENDPOINT_COLUMNS
PAGE = c.ROOT / "docs" / "catalogue-size-result.md"
SLUG = "synthia_tuberculosis"
NAME = "Synthia tuberculosis"
# antibiotic, in_M, kind, genomes, resistant, graded variants when they follow the resistant genomes, and when they
# follow the resistant share
DRUGS = [
    ("drug a", "0", "primary", 100, 10, 1, 3),
    ("drug b", "0", "primary", 400, 20, 2, 2),
    ("drug c", "0", "primary", 100, 30, 3, 6),
    ("drug d", "0", "primary", 1000, 40, 4, 1),
    ("drug e", "0", "primary", 200, 50, 5, 5),
    ("drug f", "0", "proxy", 500, 60, 6, 4),
    ("drug m", "1", "primary", 300, 70, 7, 7),
    ("drug z", "0", "primary", 150, 15, 0, 0),
    ("drug n", "0", "none", 120, 12, None, None),
]
IN_RULE = ["drug a", "drug b", "drug c", "drug d", "drug e", "drug f"]
# the share's ranks over IN_RULE are 3, 2, 6, 1, 5, 4: squared rank differences against 1 to 6 sum to 26
SHARE_RHO = 1 - 6 * 26 / (6 * 35)


def dataset(antibiotic):
    return f"{SLUG}__{antibiotic.replace(' ', '_')}"


def gid(k, i):
    return f"9{k}{i:02d}.1"


class World:
    """A SYNTHETIC species: page-endpoints.tsv, the release reports, the species files, and the reading folders."""

    def __init__(self, top, by="resistant", drop_report=None, extra_report=False, bad_sum=None, bad_page=None,
                 name_of=None):
        self.top = Path(top)
        self.reports, self.species = self.top / "reports", self.top / "species"
        self.reports.mkdir()
        self.species.mkdir()
        self.endpoints = self.top / "page-endpoints.tsv"
        self.reading, self.out = self.top / "reading", self.top / "out"
        self.ids = []
        rows = []
        for k, (ab, in_m, kind, genomes, resistant, by_r, by_s) in enumerate(DRUGS):
            ab = (name_of or {}).get(ab, ab)
            cv = by_r if by == "resistant" else by_s
            read = kind != "none"
            rwr = resistant + (1 if bad_page == ab else 0)
            rows.append({"dataset": dataset(ab), "species": NAME, "antibiotic": ab, "in_M": in_m,
                         "labels_depend_on": "the maintainer's rule" if in_m == "1" else "none of the four",
                         "kind": kind, "tool_drugs": ab if read else "", "resistant_with_result": rwr if read else "",
                         "resistant_no_result": 0 if read else "", "catalogue_variants": cv if read else "",
                         "genes_named": (1 if cv else 0) if read else "",
                         "without_determinant": resistant // 5 if read else "",
                         "share_without": cs.f12((resistant // 5) / rwr) if read else "",
                         "groups_main": 3 if read else ""})
            samples = [{"id": gid(k, i), "label": "R" if i < 3 else "S"} for i in range(5)]
            self.ids += [s["id"] for s in samples]
            susceptible = genomes - resistant + (1 if bad_sum == ab else 0)
            if drop_report != ab:
                self.report(dataset(ab), genomes, resistant, susceptible, samples)
        if extra_report:
            self.report(f"{SLUG}__drug_x", 50, 5, 45, [])
        self.report("othera_coli__drug_a", 80, 8, 72, [{"id": "7777.1"}])
        cs.write_tsv(self.endpoints, rows, ENDPOINT_COLUMNS)
        (self.species / f"{SLUG}.json").write_text(json.dumps({"genomes": [{"id": i} for i in self.ids]}))
        (self.species / "othera_coli.json").write_text(json.dumps({"genomes": [{"id": "7777.1"}]}))

    def report(self, name, genomes, resistant, susceptible, samples):
        body = {"species": NAME, "counts": {"samples": genomes, "resistant": resistant, "susceptible": susceptible,
                                            "k": 31, "features": 9},
                "samples": samples, "splits": [{"samples": [s["id"] for s in samples]}]}
        (self.reports / f"{name}.json").write_text(json.dumps(body))

    def read(self):
        lines = cs.read(self.reports, self.endpoints, self.reading)
        text = ["python3 -I scripts/mechanism/catalogue_size.py read REPORTS ENDPOINTS READING"]
        for name in ("release-counts.tsv", "read.txt"):
            text.append(f"{hashlib.sha256((self.reading / name).read_bytes()).hexdigest()}  {name}")
        (self.reading / "reading.txt").write_text("\n".join(text) + "\n")
        return lines

    def write(self):
        self.read()
        return cs.write(self.reading, self.endpoints, self.species, self.out)

    def table(self):
        return {r["antibiotic"]: r for r in cs.load(self.out / "antibiotics.tsv")}

    def correlations(self):
        return cs.load(self.out / "correlations.tsv")


class Case(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.top = Path(self._tmp.name)

    def tearDown(self):
        self._tmp.cleanup()

    def world(self, **kw):
        return World(self.top, **kw)


class Spearman(unittest.TestCase):
    def test_a_case_with_no_ties_gives_its_known_answer(self):
        self.assertAlmostEqual(cs.spearman([1, 2, 3, 4, 100], [1, 30, 20, 500, 400]), 0.8, places=12)

    def test_a_reversed_case_gives_minus_one(self):
        self.assertAlmostEqual(cs.spearman([1, 2, 3, 4, 5], [50, 40, 30, 20, 10]), -1.0, places=12)

    def test_a_case_with_ties_gives_its_known_answer(self):
        self.assertAlmostEqual(cs.spearman([1, 2, 2, 4], [1, 2, 3, 4]), math.sqrt(0.9), places=12)

    def test_tied_values_take_the_average_of_their_ranks(self):
        self.assertEqual(cs.ranks([10, 20, 20, 30, 20]), [1.0, 3.0, 3.0, 5.0, 3.0])
        self.assertEqual(cs.ranks([5, 5]), [1.5, 1.5])

    def test_a_column_of_tied_values_stops(self):
        with self.assertRaisesRegex(cs.Stop, "all tied"):
            cs.spearman([1, 2, 3], [4, 4, 4])

    def test_leave_one_out_gives_the_smallest_and_largest_with_the_antibiotic_dropped(self):
        # over all four 1 - 6 * 2 / (4 * 15) = 0.8; dropping a or b gives 0.5, dropping c or d gives 1
        x, y = [1, 2, 3, 4], [1, 2, 4, 3]
        self.assertAlmostEqual(cs.spearman(x, y), 0.8, places=12)
        self.assertEqual(cs.leave_one_out(x, y, ["a", "b", "c", "d"]), (0.5, "a", 1.0, "c"))

    def test_the_proof_passes_on_the_three_cases(self):
        self.assertEqual(len(cs.CASES), 3)
        self.assertTrue(any(len(set(x)) < len(x) for _, x, _, _ in cs.CASES))
        self.assertEqual(cs.prove_spearman(), "4 rank correlation: 3 of 3 hand-made cases give their known answers")


class Rule(Case):
    def test_a_resistant_correlation_as_strong_as_the_share_selects_resistant_genomes(self):
        self.assertEqual(cs.branch(0.6, 0.6), "tracks resistant genomes")
        self.assertEqual(cs.branch(-0.6, 0.6), "tracks resistant genomes")
        self.assertEqual(cs.branch(0.9, 0.2), "tracks resistant genomes")

    def test_a_share_stronger_by_the_margin_selects_share(self):
        self.assertEqual(cs.branch(0.4, 0.6), "tracks share")
        self.assertLess(0.6 - 0.45, 0.15)
        self.assertEqual(cs.branch(0.45, 0.6), "tracks share")

    def test_a_share_stronger_by_less_than_the_margin_selects_neither(self):
        self.assertEqual(cs.branch(0.5, 0.6), "neither by the margin")
        self.assertEqual(cs.branch(0.46, 0.6), "neither by the margin")

    def test_the_rule_reads_absolute_values(self):
        self.assertEqual(cs.branch(0.3, -0.6), "tracks share")
        self.assertEqual(cs.branch(-0.7, 0.6), "tracks resistant genomes")

    def test_the_page_states_the_rule_the_script_applies(self):
        text = PAGE.read_text(encoding="utf-8")
        section = re.split(r"^## ", text, flags=re.M)
        rule = [s for s in section if s.startswith("Rule\n")]
        self.assertEqual(len(rule), 1)
        self.assertEqual(rule[0][len("Rule\n"):].strip(), cs.RULE)
        self.assertIn(f"by {cs.MARGIN} or more", cs.RULE)
        self.assertIn(f"by less than {cs.MARGIN}", cs.RULE)
        self.assertEqual(float(re.search(r"within (0\.0+1) of a bound", cs.RULE).group(1)), cs.TOL)
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            self.assertEqual(cs.main(["plan"]), 0)
        self.assertEqual(out.getvalue(), cs.RULE + "\n")

    def test_a_synthetic_reading_selects_resistant_genomes(self):
        w = self.world(by="resistant")
        lines, _ = w.write()
        rule = {r["x"]: r for r in w.correlations() if r["set"] == "rule"}
        self.assertEqual(rule["resistant"]["rho"], "1.000000000000")
        self.assertEqual(rule["resistant_share"]["rho"], cs.f12(SHARE_RHO))
        self.assertEqual({r["comparison_gives"] for r in rule.values()}, {"tracks resistant genomes"})
        self.assertEqual({x: r["decides"] for x, r in rule.items()},
                         {"resistant_share": "yes", "resistant": "yes", "genomes": "no"})

    def test_a_synthetic_reading_selects_share(self):
        w = self.world(by="share")
        w.write()
        rule = {r["x"]: r for r in w.correlations() if r["set"] == "rule"}
        self.assertEqual(rule["resistant_share"]["rho"], "1.000000000000")
        self.assertEqual(rule["resistant"]["rho"], cs.f12(SHARE_RHO))
        self.assertEqual({r["comparison_gives"] for r in rule.values()}, {"tracks share"})


class Read(Case):
    def test_a_page_antibiotic_with_no_report_stops(self):
        w = self.world(drop_report="drug c")
        with self.assertRaisesRegex(cs.Stop, r"guard 1 match: 1 page rows with no report \(synthia_tuberculosis__drug_c\)"):
            w.read()
        self.assertFalse(w.reading.exists())

    def test_a_report_with_no_page_row_stops(self):
        w = self.world(extra_report=True)
        with self.assertRaisesRegex(cs.Stop, r"guard 1 match: 0 page rows with no report \(none\), 1 reports"):
            w.read()

    def test_a_report_of_another_species_is_not_matched(self):
        w = self.world()
        lines = w.read()
        self.assertEqual(lines[0], "1 match: 9 page rows, 9 reports of the species, 9 matched, 0 page rows with no "
                                   "report, 0 reports with no page row")

    def test_counts_that_do_not_add_up_stop(self):
        w = self.world(bad_sum="drug e")
        with self.assertRaisesRegex(cs.Stop, r"guard 2 sums: 1 reports do not add up \(synthia_tuberculosis__drug_e\)"):
            w.read()

    def test_a_resistant_count_unlike_the_page_stops(self):
        w = self.world(bad_page="drug b")
        with self.assertRaisesRegex(cs.Stop, r"guard 6 resistant: 1 reports differ from the page"):
            w.read()

    def test_the_samples_list_is_never_kept(self):
        w = self.world()
        path = w.reports / f"{dataset('drug a')}.json"
        with open(path, encoding="utf-8") as handle:
            parsed = json.load(handle, object_pairs_hook=cs.drop_samples)
        self.assertNotIn("samples", parsed)
        self.assertEqual(parsed["counts"]["samples"], 100)
        self.assertNotIn(gid(0, 0), json.dumps(parsed["counts"]))
        self.assertEqual(cs.counts_of(path), {"samples": 100, "resistant": 10, "susceptible": 90})
        w.read()
        for f in w.reading.iterdir():
            text = f.read_text(encoding="utf-8")
            self.assertFalse([i for i in w.ids if i in text], f.name)

    def test_the_read_writes_the_three_counts_per_report(self):
        w = self.world()
        lines = w.read()
        rows = {r["dataset"]: r for r in cs.load(w.reading / "release-counts.tsv")}
        self.assertEqual(list(cs.load(w.reading / "release-counts.tsv")[0]), cs.RELEASE_COLUMNS)
        self.assertEqual(len(rows), 9)
        self.assertEqual(rows[dataset("drug d")], {"dataset": dataset("drug d"), "genomes": "1000", "resistant": "40",
                                                   "susceptible": "960"})
        self.assertEqual(lines[1:], ["2 sums: 9 of 9 reports add up",
                                     "6 resistant: 8 of 8 reports the tool reads agree with the page"])

    def test_a_reading_folder_not_empty_stops(self):
        w = self.world()
        w.reading.mkdir()
        (w.reading / "left.txt").write_text("x")
        with self.assertRaisesRegex(cs.Stop, "not empty"):
            w.read()


class Write(Case):
    def test_set_m_is_in_the_table_and_out_of_the_correlations(self):
        w = self.world()
        w.write()
        m = w.table()["drug m"]
        self.assertEqual((m["in_M"], m["genomes"], m["resistant"], m["susceptible"], m["catalogue_variants"]),
                         ("1", "300", "70", "230", "7"))
        self.assertEqual((m["in_correlation"], m["left_out_because"]), ("0", "set M"))
        for r in w.correlations():
            self.assertNotIn("drug m", r["antibiotics"].split("; "))

    def test_an_antibiotic_with_no_graded_variant_is_left_out(self):
        w = self.world()
        w.write()
        t = w.table()
        self.assertEqual((t["drug z"]["in_correlation"], t["drug z"]["left_out_because"]), ("0", "no graded variant"))
        self.assertEqual((t["drug n"]["in_correlation"], t["drug n"]["left_out_because"]),
                         ("0", "not read by the tool"))
        self.assertEqual((t["drug n"]["genomes"], t["drug n"]["catalogue_variants"]), ("120", ""))
        self.assertEqual(sorted(a for a, r in t.items() if r["in_correlation"] == "1"), IN_RULE)

    def test_the_resistant_share_is_over_all_genomes(self):
        w = self.world()
        w.write()
        t = w.table()
        self.assertEqual(t["drug a"]["resistant_share"], "0.100000000000")
        self.assertEqual(t["drug d"]["resistant_share"], "0.040000000000")

    def test_the_columns_of_the_page_are_carried_unchanged(self):
        w = self.world()
        w.write()
        page = {r["antibiotic"]: r for r in cs.load(w.endpoints)}
        t = w.table()
        self.assertEqual(list(cs.load(w.out / "antibiotics.tsv")[0]), cs.TABLE_COLUMNS)
        self.assertEqual(list(t), [d[0] for d in DRUGS])
        for ab, r in t.items():
            self.assertEqual(r["species"], NAME)
            for key in cs.CARRIED:
                self.assertEqual(r[key], page[ab][key], (ab, key))

    def test_every_correlation_of_a_set_uses_the_same_antibiotics(self):
        w = self.world()
        lines, _ = w.write()
        rows = w.correlations()
        self.assertEqual(list(rows[0]), cs.CORRELATION_COLUMNS)
        self.assertEqual([(r["set"], r["x"]) for r in rows], [(s, x) for s in cs.SETS for x in cs.PAIRS])
        rule = [r for r in rows if r["set"] == "rule"]
        primary = [r for r in rows if r["set"] == "primary only"]
        self.assertEqual({(r["n"], r["antibiotics"]) for r in rule}, {("6", "; ".join(IN_RULE))})
        self.assertEqual({(r["n"], r["antibiotics"]) for r in primary}, {("5", "; ".join(IN_RULE[:5]))})
        self.assertEqual({r["decides"] for r in primary}, {"no"})
        self.assertEqual({r["role"] for r in rule}, {"the rule's set"})
        self.assertEqual({r["role"] for r in primary}, {"beside, decides nothing"})
        self.assertIn("3 sets: rule: resistant_share n=6, resistant n=6, genomes n=6; antibiotics: drug a; drug b; "
                      "drug c; drug d; drug e; drug f", lines)

    def test_a_pair_over_another_set_stops(self):
        w = self.world()

        def wider(rows, name, x, _orig=cs.rows_for):
            return [r for r in rows if r["catalogue_variants"] not in ("", "0")] if x == "genomes" else \
                _orig(rows, name, x)

        with mock.patch.object(cs, "rows_for", wider), self.assertRaisesRegex(cs.Stop, "guard 3 sets"):
            w.write()
        self.assertFalse(w.out.exists())

    def test_a_wrong_rank_correlation_stops_before_any_correlation(self):
        w = self.world()
        with mock.patch.object(cs, "spearman", cs.pearson), self.assertRaisesRegex(cs.Stop, "guard 4"):
            w.write()
        self.assertFalse(w.out.exists())

    def test_a_genome_id_in_an_output_stops(self):
        w = self.world(name_of={"drug b": f"drug {gid(1, 2)}"})
        with self.assertRaisesRegex(cs.Stop, r"guard 5 ids: \d+ hits"):
            w.write()

    def test_the_scan_reads_every_species_file(self):
        w = self.world()
        lines, _ = w.write()
        # OUT's five files before SHA256SUMS, and READING's three
        self.assertEqual(lines[-1], f"5 ids: 8 files, 0 hits, against {len(w.ids) + 1} genome ids of 2 species files")

    def test_a_reading_txt_whose_sha256_differs_stops(self):
        w = self.world()
        w.read()
        with open(w.reading / "read.txt", "a", encoding="utf-8") as handle:
            handle.write("x\n")
        with self.assertRaisesRegex(cs.Stop, "read.txt does not have the sha256"):
            cs.write(w.reading, w.endpoints, w.species, w.out)

    def test_an_out_folder_not_empty_stops(self):
        w = self.world()
        w.out.mkdir()
        (w.out / "left.tsv").write_text("x")
        with self.assertRaisesRegex(cs.Stop, "not empty"):
            w.write()

    def test_the_record_lists_every_file_in_sha256sums(self):
        w = self.world()
        lines, _ = w.write()
        listed = dict(line.split("  ", 1)[::-1] for line in (w.out / "SHA256SUMS").read_text().splitlines())
        self.assertEqual(sorted(listed), ["README.md", "antibiotics.tsv", "correlations.tsv", "guards.txt",
                                          "reading.txt"])
        for name, digest in listed.items():
            self.assertEqual(hashlib.sha256((w.out / name).read_bytes()).hexdigest(), digest)
        self.assertEqual((w.out / "guards.txt").read_text().splitlines(), lines)
        self.assertEqual([line.split(" ", 1)[0] for line in lines], ["1", "2", "6", "3", "3", "4", "5"])

    def test_print_gives_the_page_tables(self):
        w = self.world()
        w.write()
        t = cs.tables(w.out)
        self.assertEqual(list(t), ["antibiotics", "correlations", "sets", "guards"])
        self.assertEqual(len(t["antibiotics"]), 2 + len(DRUGS))
        self.assertIn("| drug m | M | primary | 300 | 70 | 230 | 0.233 | 7 | 1 | 14 (0.200) | no: set M |",
                      t["antibiotics"])
        self.assertIn("| drug n | not M | none | 120 | 12 | 108 | 0.100 | none | none | none | no: not read by the tool |",
                      t["antibiotics"])
        self.assertIn("| rule | the rule's set | resistant genomes | 6 | 1.000 | 1.000 (drug a) | 1.000 (drug a) | yes | "
                      "tracks resistant genomes |", t["correlations"])
        self.assertEqual({line.split(" | ")[1] for line in t["correlations"][2:] if line.startswith("| primary only |")},
                         {"beside, decides nothing"})
        self.assertEqual(len(t["correlations"]), 2 + 6)
        self.assertEqual(t["sets"][2], "| rule | 6 | drug a, drug b, drug c, drug d, drug e, drug f |")


if __name__ == "__main__":
    unittest.main()
