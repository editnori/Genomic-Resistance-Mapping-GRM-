"""scripts/mechanism/overlap_second.py, the learner overlap's second pass (docs/learner-overlap-result.md).

Every report, genome id, label, score, species and antibiotic here is SYNTHETIC and made up in the test, with the
fixtures of tests/mechanism/test_learner_overlap.py. Each count below is worked by hand. With the fixture's scores
of 0.1 and 0.9, a genome k of the five miss has the score sum 4.5 - 0.8k, so the mean misses it when k is 3 or more.
"""

import contextlib
import hashlib
import io
import json
import os
import tempfile
import unittest
from unittest import mock

from tests.mechanism import common as c
from tests.mechanism.test_learner_overlap import AUREUS, COLI, FIVE, Report, lay_out, standard

lo = c.load("learner_overlap")
os2 = c.load("overlap_second")


def run(module, *argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = module.main([str(x) for x in argv])
    return rc, buf.getvalue()


def edit_tsv(path, change):
    rows = lo.read_tsv(path)
    with open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
    rows = change(rows)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(head) + "\n" + "".join("\t".join(r[h] for h in head) + "\n" for r in rows))


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = self.tmp.name
        self.first = os.path.join(self.top, "first")
        self.reading = os.path.join(self.top, "reading2")
        self.record = os.path.join(self.top, "record2")

    def tearDown(self):
        self.tmp.cleanup()

    def first_pass(self, reports):
        d = lay_out(self.top, reports)
        rc, out = run(lo, "read", d, "--out", os.path.join(self.top, "reading1"), "--expect", len(reports))
        self.assertEqual(rc, 0, out)
        rc, out = run(lo, "write", os.path.join(self.top, "reading1"), "--out", self.first)
        self.assertEqual(rc, 0, out)
        return d

    def read(self, reports=None, before=None):
        reports = standard() if reports is None else reports
        d = self.first_pass(reports)
        if before:
            before(d)
        return run(os2, "read", d, "--first", self.first, "--out", self.reading, "--expect", len(reports))

    def read_write(self, reports=None):
        rc, out = self.read(reports)
        self.assertEqual(rc, 0, out)
        rc, out2 = run(os2, "write", self.reading, "--first", self.first, "--out", self.record)
        self.assertEqual(rc, 0, out2)
        rows = {r["antibiotic"]: r for r in lo.read_tsv(os.path.join(self.record, "vote.tsv"))}
        pooled = lo.read_tsv(os.path.join(self.record, "vote-pooled.tsv"))
        return out + out2, rows, pooled

    def stops(self, needle, reports=None, before=None):
        rc, out = self.read(reports, before)
        self.assertEqual(rc, 1, out)
        self.assertIn(needle, out)
        self.assertFalse(os.path.exists(self.reading))
        return out


def mean_cols(row):
    return [row[f"mean_{w}"] for w in ("miss", "miss_by_all5", "miss_by_1to4", "miss_by_none",
                                       "wrong", "wrong_by_all5", "wrong_by_1to4", "wrong_by_none")]


class Vote(Case):
    def test_the_mean_misses_what_three_or_more_of_the_five_miss_here(self):
        _, rows, _ = self.read_write()
        # drug a: genomes 0 to 2 missed by all five, genome 3 by three, genomes 4 and 5 by one; susceptible genome 12
        # called by all five, genome 13 by one
        self.assertEqual(mean_cols(rows["drug a"]), ["4", "3", "1", "0", "1", "1", "0", "0"])
        self.assertEqual(mean_cols(rows["drug b"]), ["0"] * 8)
        self.assertEqual(mean_cols(rows["drug c"]), ["3", "0", "3", "0", "0", "0", "0", "0"])
        self.assertEqual(mean_cols(rows["drug d"]), ["1", "1", "0", "0", "0", "0", "0", "0"])
        self.assertEqual(mean_cols(rows["drug h"]), ["1", "1", "0", "0", "0", "0", "0", "0"])

    def test_a_sum_of_exactly_2_5_is_a_call(self):
        r = Report(COLI, "drug a", 10, 5)
        for n in FIVE:
            r.scores[n][0] = 0.5
        _, rows, _ = self.read_write([r])
        self.assertEqual(mean_cols(rows["drug a"]), ["0"] * 8)

    def test_the_mean_is_not_the_majority(self):
        r = Report(COLI, "drug a", 10, 5)
        for n, x, y in zip(FIVE, (0.0, 0.6, 0.6, 0.6, 0.6), (1.0, 0.4, 0.4, 0.4, 0.4)):
            r.scores[n][0], r.scores[n][10] = x, y
        _, rows, _ = self.read_write([r])
        # resistant genome 0: four of five call it, the mean is 0.48; susceptible genome 10: one calls it, mean 0.52
        self.assertEqual(mean_cols(rows["drug a"]), ["1", "0", "1", "0", "1", "0", "1", "0"])

    def test_the_vote_edge_is_counted_and_the_read_goes_on(self):
        r = Report(COLI, "drug a", 10, 5)
        for n, x in zip(FIVE, (0.1, 0.1, 0.95, 0.95, 0.4)):
            r.scores[n][0] = x
        out, rows, _ = self.read_write([r])
        # added in learner order the sum is 2.4999999999999996; math.fsum gives 2.5, a call
        self.assertIn("second check 5 the vote's edge: 1 genomes where a plain sum over 5 gives another call than "
                      "math.fsum; 1 genomes whose sum lies within 1e-9 of 2.5", out)
        self.assertEqual(mean_cols(rows["drug a"]), ["0"] * 8)

    def test_a_report_without_the_five_keeps_a_blank_row(self):
        _, rows, _ = self.read_write()
        g = rows["drug g"]
        self.assertEqual((g["species"], g["genomes"], g["resistant"], g["susceptible"]), (AUREUS, "15", "10", "5"))
        self.assertEqual({g[col] for col in os2.VOTE_COLUMNS[5:]}, {""})
        self.assertEqual(sorted(rows), ["drug a", "drug b", "drug c", "drug d", "drug g", "drug h"])

    def test_a_vote_that_parts_from_the_five_stops_the_read(self):
        with mock.patch.object(os2, "vote", lambda scores: True):
            self.stops("second check 3 arithmetic of the mean: 5 reports with the five, 5 where the mean misses")


class NullCut(Case):
    def test_the_null_fit_cut_is_worked_by_hand(self):
        _, rows, _ = self.read_write()
        a = rows["drug a"]
        # gbt has a null fit in fold 1, which holds resistant genomes 1, 6 and 11; nine stay
        self.assertEqual([a["nonull_resistant"]] + [a[f"nonull_miss_{n}"] for n in FIVE] + [a["nonull_miss_all5"]],
                         ["9", "5", "3", "3", "2", "2", "2"])
        # 9 x 5/9 x 3/9 x 3/9 x 2/9 x 2/9 = 20/729; 2 / (20/729) = 72.9
        self.assertEqual((a["nonull_expected_all5"], a["nonull_ratio_all5"]), ("0.027434842250", "72.900000000000"))
        b = rows["drug b"]
        self.assertEqual((b["nonull_resistant"], b["nonull_expected_all5"], b["nonull_ratio_all5"]),
                         ("10", "0.000100000000", "0.000000000000"))
        self.assertEqual((rows["drug c"]["nonull_expected_all5"], rows["drug c"]["nonull_ratio_all5"]),
                         ("0.000000000000", ""))

    def test_null_folds_of_two_learners_are_both_cut(self):
        r = Report(COLI, "drug n", 10, 5, null={"gbt": (1,), "lasso": (2,)}, misses={n: (0, 1, 2, 3) for n in FIVE})
        _, rows, _ = self.read_write([r])
        n = rows["drug n"]
        # folds 1 and 2 hold resistant genomes 1, 6, 2 and 7; of the misses 0 to 3, genomes 0 and 3 stay
        self.assertEqual([n["nonull_resistant"]] + [n[f"nonull_miss_{x}"] for x in FIVE] + [n["nonull_miss_all5"]],
                         ["6", "2", "2", "2", "2", "2", "2"])
        # 6 x (2/6)^5 = 2/81; 2 / (2/81) = 81
        self.assertEqual((n["nonull_expected_all5"], n["nonull_ratio_all5"]), ("0.024691358025", "81.000000000000"))

    def test_a_cut_that_does_not_add_up_stops_the_read(self):
        with mock.patch.object(os2, "null_union", lambda learners: set()):
            self.stops("second check 4 null-fit cut: 1 reports where the all-five misses outside and inside null-fit "
                       "folds do not add up: STOP")


class SecondStops(Case):
    def test_a_first_record_that_differs_stops_the_read(self):
        def change(d):
            edit_tsv(os.path.join(self.first, "datasets.tsv"),
                     lambda rows: [dict(r, miss_1="4") if r["antibiotic"] == "drug b" else r for r in rows])
        self.stops("1 that differ: STOP", before=change)

    def test_a_first_record_without_a_report_stops_the_read(self):
        def change(d):
            edit_tsv(os.path.join(self.first, "datasets.tsv"),
                     lambda rows: [r for r in rows if r["antibiotic"] != "drug h"])
        self.stops("1 reports stand in one of the read and the record only", before=change)

    def test_a_failed_first_pass_check_stops_the_second_read(self):
        def change(d):
            name = "synthia_coli__drug_b.json"
            with open(os.path.join(d, name), encoding="utf-8") as f:
                doc = json.load(f)
            doc["counts"]["resistant"] = 11
            with open(os.path.join(d, name), "w", encoding="utf-8") as f:
                json.dump(doc, f)
        self.stops("check 3 counts: 1 reports whose counts or labels differ from their samples", before=change)

    def test_an_id_shaped_token_in_the_reading_stops_the_write(self):
        rc, out = self.read([Report(COLI, "drug a", 10, 5)])
        self.assertEqual(rc, 0, out)
        with open(os.path.join(self.reading, "note.txt"), "w", encoding="utf-8") as f:
            f.write("a genome 990001.11 here\n")
        rc, out = run(os2, "write", self.reading, "--first", self.first, "--out", self.record)
        self.assertEqual(rc, 1, out)
        self.assertIn("second check 6 ids: 1 id-shaped tokens or barred keys: STOP", out)
        self.assertFalse(os.path.exists(self.record))
        self.assertEqual([n for n in os.listdir(self.top) if n.startswith(".overlap2-")], [])


class SecondPooled(Case):
    def test_the_pooled_row_is_worked_by_hand(self):
        _, _, pooled = self.read_write()
        p = pooled[0]
        # drugs a, b, c and h: drug d has 4 resistant genomes and drug g four learners
        self.assertEqual((p["scope"], p["datasets"], p["resistant"], p["susceptible"]), ("all", "4", "42", "48"))
        self.assertEqual(mean_cols(p), ["8", "4", "4", "0", "1", "1", "0", "0"])
        self.assertEqual((p["nonull_resistant"], p["nonull_miss_all5"], p["nonull_expected_all5"]),
                         ("39", "3", "0.027634842250"))
        # mean shares 4/12, 0, 3/10, 1/10: the median is 0.2
        self.assertEqual(p["mean_share_median"], "0.200000000000")
        # ratios outside null folds 0, 72.9 and 10000; drug c has none
        self.assertEqual((p["nonull_datasets_with_ratio"], p["nonull_ratio_q1"], p["nonull_ratio_q3"],
                          p["nonull_datasets_ratio_over_1"]), ("3", "36.450000000000", "5036.450000000000", "2"))
        self.assertEqual((p["paired_datasets"], p["paired_ratio_median"], p["paired_nonull_ratio_median"]),
                         ("3", "72.000000000000", "72.900000000000"))
        self.assertEqual([r["scope"] + ":" + r["species"] for r in pooled], ["all:", f"species:{COLI}"])

    def test_the_median_outside_null_folds_weights_each_report_once(self):
        _, _, pooled = self.read_write()
        self.assertEqual(pooled[0]["nonull_ratio_median"], "72.900000000000")

    def test_a_report_under_10_resistant_stays_out_of_the_pool(self):
        _, rows, pooled = self.read_write()
        self.assertEqual(rows["drug d"]["mean_miss"], "1")
        self.assertEqual((pooled[0]["datasets"], pooled[0]["mean_miss"]), ("4", "8"))

    def test_the_paired_medians_take_only_reports_with_both_ratios(self):
        # drug p: all five miss genome 0 in fold 0, where gbt has a null fit, and gbt misses nothing else, so the
        # report keeps its first-pass ratio (625) and has none outside null folds
        p = Report(COLI, "drug p", 10, 5, null={"gbt": (0,)},
                   misses=dict({n: (0, 1) for n in FIVE if n != "gbt"}, gbt=(0,)))
        _, rows, pooled = self.read_write(standard() + [p])
        self.assertEqual((rows["drug p"]["nonull_expected_all5"], rows["drug p"]["nonull_ratio_all5"]),
                         ("0.000000000000", ""))
        first = lo.read_tsv(os.path.join(self.first, "pooled.tsv"))[0]
        self.assertEqual((first["datasets"], first["datasets_with_ratio"]), ("5", "4"))
        self.assertEqual((pooled[0]["nonull_datasets_with_ratio"], pooled[0]["paired_datasets"],
                          pooled[0]["paired_ratio_median"]), ("3", "3", "72.000000000000"))


class SecondRecord(Case):
    def test_the_record_holds_its_files_checksums_and_gates(self):
        out, _, _ = self.read_write()
        self.assertEqual(sorted(os.listdir(self.record)), sorted(list(os2.RECORD_FILES) + ["SHA256SUMS"]))
        with open(os.path.join(self.record, "SHA256SUMS"), encoding="utf-8") as f:
            for line in f.read().splitlines():
                digest, name = line.split("  ")
                with open(os.path.join(self.record, name), "rb") as g:
                    self.assertEqual(hashlib.sha256(g.read()).hexdigest(), digest, name)
        with open(os.path.join(self.record, "gates.txt"), encoding="utf-8") as f:
            gates = f.read().splitlines()
        self.assertEqual([g.split(":")[0] for g in gates if g.startswith("second check")],
                         [f"second check {n} " + w for n, w in
                          zip(range(1, 7), ("first pass", "first-pass columns", "arithmetic of the mean",
                                            "null-fit cut", "the vote's edge", "ids"))])
        self.assertTrue(any(g.startswith("check 1 parse: 8 files found") for g in gates))
        self.assertFalse(any("STOP" in g for g in gates))
        self.assertEqual(gates[-1], "reading rule outside null-fit folds: stack")

    def test_no_genome_id_stands_in_any_written_file(self):
        reports = standard()
        self.read_write(reports)
        written = lo.files_under(self.reading) + lo.files_under(self.record)
        self.assertEqual(len(written), 5)
        for path in written:
            with open(path, encoding="utf-8") as f:
                text = f.read()
            self.assertNotIn('"id"', text)
            for r in reports:
                for g in r.ids:
                    self.assertNotIn(g, text, os.path.basename(path))

    def test_a_second_run_writes_the_same_bytes(self):
        self.read_write()
        again = os.path.join(self.top, "again")
        rc, out = run(os2, "write", self.reading, "--first", self.first, "--out", again)
        self.assertEqual(rc, 0, out)
        for name in os2.RECORD_FILES + ("SHA256SUMS",):
            with open(os.path.join(self.record, name), "rb") as f, open(os.path.join(again, name), "rb") as g:
                self.assertEqual(f.read(), g.read(), name)

    def test_print_gives_the_page_tables(self):
        self.read_write()
        rc, out = run(os2, "print", self.record, "--first", self.first)
        self.assertEqual(rc, 0, out)
        for line in ("## mean vote", "| missed by the mean | 8 |", "| of those, missed by all five | 4 |",
                     "| of those, missed by none | 0 |", "| called resistant by the mean | 1 |",
                     "## null-fit cut", "| resistant genomes | 42 | 39 |", "| missed by all five | 4 | 3 |",
                     "| median ratio | 72.0 | 72.9 |", "| median ratio over those | 72.0 | 72.9 |",
                     "## mean vote by species", f"| *{COLI}* | 4 | 42 | 8 | 4 | 4 | 72.9 |"):
            self.assertIn(line + "\n", out)

    def test_the_plan_names_the_vote_and_the_cut(self):
        rc, out = run(os2, "plan")
        self.assertEqual(rc, 0, out)
        plan = json.loads(out)
        self.assertEqual((plan["vote"], plan["vote_sum"], plan["cutoff"], plan["scheme"]),
                         ("math.fsum(five scores) >= 2.5", 2.5, 0.5, "lineage"))
        self.assertEqual(plan["vote_columns"][-2:], ["nonull_expected_all5", "nonull_ratio_all5"])


if __name__ == "__main__":
    unittest.main()
