"""scripts/mechanism/learner_overlap.py, the table script of docs/learner-overlap-result.md.

Every report, genome id, label, score, species and antibiotic here is
SYNTHETIC and made up in the test. Each count below is worked by hand from
the fixture's design.
"""

import contextlib
import io
import json
import os
import tempfile
import unittest

from tests.mechanism import common as c

lo = c.load("learner_overlap")
FIVE = ("scm", "cart", "rf", "gbt", "lasso")
COLI, AUREUS = "Synthia coli", "Synthia aureus"


def confusion(labels, scores):
    out = {"tp": 0, "fn": 0, "fp": 0, "tn": 0}
    for t, x in zip(labels, scores):
        out[("tp" if x >= 0.5 else "fn") if t else ("fp" if x >= 0.5 else "tn")] += 1
    return out


def scores_for(labels, misses=(), wrongs=(), low=0.1, high=0.9):
    """A learner's scores: a resistant genome in `misses` and a susceptible one outside `wrongs` score low."""
    return [(low if i in misses else high) if t else (high if i in wrongs else low) for i, t in enumerate(labels)]


class Report:
    """A SYNTHETIC revision 14 report, built field by field so a test can break one."""

    def __init__(self, species, antibiotic, resistant, susceptible, misses=None, wrongs=None, no_catalog=(),
                 learners=FIVE, null=None, tag=1):
        self.labels = [True] * resistant + [False] * susceptible
        n = len(self.labels)
        self.catalog = [t and i not in no_catalog for i, t in enumerate(self.labels)]
        self.folds = [i % 5 for i in range(n)]
        misses, wrongs, null = misses or {}, wrongs or {}, null or {}
        self.scores = {name: scores_for(self.labels, misses.get(name, ()), wrongs.get(name, ())) for name in learners}
        self.random = {name: scores_for(self.labels) for name in learners}
        self.complexity = {name: [0.0 if f in null.get(name, ()) else 3.0 for f in range(5)] for name in learners}
        self.no_call = {}
        for name in learners:
            nc = {"resistant": 0, "susceptible": 0}
            for t, f in zip(self.labels, self.folds):
                if f in null.get(name, ()):
                    nc["resistant" if t else "susceptible"] += 1
            self.no_call[name] = nc
        self.species, self.antibiotic, self.tag = species, antibiotic, tag
        self.schemes = ("random", "lineage")
        self.counts = {"samples": n, "resistant": resistant, "susceptible": susceptible}
        self.metric_change = {}
        self.catalog_scores = None
        self.ids = [f"99{tag:04d}.{i + 11}" for i in range(n)]

    def entries(self, scheme):
        out = []
        for name, sc in (self.random if scheme == "random" else self.scores).items():
            m = confusion(self.labels, sc)
            m.update(self.metric_change.get(name, {}))
            out.append({"learner": name, "scores": sc, "metrics": m, "no_call": self.no_call[name],
                        "fold_complexity": self.complexity[name]})
        cat = self.catalog_scores or [1.0 if x else 0.0 for x in self.catalog]
        out.append({"learner": "catalog", "scores": cat, "metrics": confusion(self.labels, cat), "fold_complexity": []})
        return out

    def doc(self):
        return {"id": f"{self.species}__{self.antibiotic}", "schema": 2, "species": self.species,
                "antibiotic": self.antibiotic, "counts": self.counts,
                "samples": [{"id": g, "resistant": t, "catalog_resistant": k, "lineage": i % 3, "year": None,
                             "country": None, "groups": []}
                            for i, (g, t, k) in enumerate(zip(self.ids, self.labels, self.catalog))],
                "splits": [{"scheme": s, "folds": 5, "n_groups": 9, "test_fold": self.folds,
                            "learners": self.entries(s)} for s in self.schemes]}


def standard():
    """Seven reports, each worked by hand in the tests."""
    a = Report(COLI, "drug a", 12, 8, tag=1, no_catalog=(0,), null={"gbt": (1,)},
               misses={"scm": (0, 1, 2, 3, 4, 5), "cart": (0, 1, 2, 3), "rf": (0, 1, 2, 3), "gbt": (0, 1, 2),
                       "lasso": (0, 1, 2)},
               wrongs={"scm": (12, 13), "cart": (12,), "rf": (12,), "gbt": (12,), "lasso": (12,)})
    b = Report(COLI, "drug b", 10, 30, tag=2, misses={n: (j,) for j, n in enumerate(FIVE)})
    cc = Report(COLI, "drug c", 10, 5, tag=3, misses={n: (0, 1, 2) for n in FIVE[1:]})
    d = Report(COLI, "drug d", 4, 4, tag=4, misses={n: (0,) for n in FIVE})
    h = Report(COLI, "drug h", 10, 5, tag=5, misses={n: (0,) for n in FIVE})
    e = Report(AUREUS, "drug e", 10, 5, tag=6)
    e.schemes = ("random",)
    f = Report(AUREUS, "drug f", 10, 5, tag=7, learners=("scm",))
    g = Report(AUREUS, "drug g", 10, 5, tag=8, learners=FIVE[:4], misses={n: (0,) for n in FIVE[:4]})
    return [a, b, cc, d, h, e, f, g]


def run(*argv):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        rc = lo.main([str(x) for x in argv])
    return rc, buf.getvalue()


def lay_out(top, reports):
    d = os.path.join(top, "reports")
    os.makedirs(d)
    for r in reports:
        with open(os.path.join(d, f"{r.species}__{r.antibiotic}.json".replace(" ", "_").lower()), "w",
                  encoding="utf-8") as f:
            json.dump(r.doc(), f)
    return d


class Case(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.top = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def read(self, reports, expect=None):
        d = lay_out(self.top, reports)
        return run("read", d, "--out", os.path.join(self.top, "reading"), "--expect",
                   len(reports) if expect is None else expect)

    def read_write(self, reports=None):
        reports = standard() if reports is None else reports
        rc, out = self.read(reports)
        self.assertEqual(rc, 0, out)
        rc, out2 = run("write", os.path.join(self.top, "reading"), "--out", os.path.join(self.top, "record"))
        self.assertEqual(rc, 0, out2)
        rows = {r["antibiotic"]: r for r in lo.read_tsv(os.path.join(self.top, "record", "datasets.tsv"))}
        pooled = lo.read_tsv(os.path.join(self.top, "record", "pooled.tsv"))
        return out + out2, rows, pooled

    def stops(self, reports, needle, expect=None):
        rc, out = self.read(reports, expect)
        self.assertEqual(rc, 1, out)
        self.assertIn(needle, out)
        self.assertFalse(os.path.exists(os.path.join(self.top, "reading")))
        return out


class Counts(Case):
    def test_all_five_miss_the_same_genomes(self):
        _, rows, _ = self.read_write()
        a = rows["drug a"]
        self.assertEqual([a[f"miss_{k}"] for k in range(6)], ["6", "2", "0", "1", "0", "3"])
        self.assertEqual([a[f"wrong_{k}"] for k in range(6)], ["6", "1", "0", "0", "0", "1"])
        self.assertEqual([a[f"miss_{n}"] for n in FIVE], ["6", "4", "4", "3", "3"])
        self.assertEqual((a["miss_all5"], a["expected_all5"], a["ratio_all5"]),
                         ("3", "0.041666666667", "72.000000000000"))

    def test_disjoint_misses_stack_on_no_genome(self):
        _, rows, _ = self.read_write()
        b = rows["drug b"]
        self.assertEqual([b[f"miss_{k}"] for k in range(6)], ["5", "5", "0", "0", "0", "0"])
        self.assertEqual((b["miss_all5"], b["expected_all5"], b["ratio_all5"]), ("0", "0.000100000000", "0.000000000000"))

    def test_a_score_of_exactly_the_cutoff_is_a_call(self):
        r = Report(COLI, "drug a", 10, 5, misses={n: (0,) for n in FIVE})
        r.scores["scm"][0] = 0.5
        r.scores["cart"][1] = 0.49999
        r.scores["rf"][12] = 0.5
        _, rows, _ = self.read_write([r])
        row = rows["drug a"]
        self.assertEqual((row["miss_4"], row["miss_1"], row["miss_scm"], row["miss_cart"]), ("1", "1", "0", "2"))
        self.assertEqual((row["wrong_1"], row["wrong_0"]), ("1", "4"))

    def test_the_k_counts_add_to_their_class(self):
        _, rows, _ = self.read_write()
        for row in rows.values():
            n = int(row["learner_count"])
            self.assertEqual(sum(int(row[f"miss_{k}"]) for k in range(n + 1)), int(row["resistant"]), row["antibiotic"])
            self.assertEqual(sum(int(row[f"wrong_{k}"]) for k in range(n + 1)), int(row["susceptible"]))
            for k in range(n + 1, 6):
                self.assertEqual((row[f"miss_{k}"], row[f"wrong_{k}"]), ("", ""))

    def test_only_the_lineage_split_is_read(self):
        # every learner calls every genome right on the random split; the lineage split's misses must stand
        _, rows, _ = self.read_write()
        self.assertEqual(rows["drug h"]["miss_5"], "1")
        self.assertEqual(rows["drug c"]["miss_4"], "3")

    def test_the_catalogue_reader_is_not_a_learner(self):
        out, rows, _ = self.read_write()
        self.assertEqual(rows["drug a"]["learners"], ",".join(FIVE))
        self.assertEqual(rows["drug a"]["learner_count"], "5")
        self.assertIn("the release's is cart,gbt,lasso,rf,scm on 5 reports", out)

    def test_misses_with_no_catalogue_call(self):
        # drug a: genome 0 is missed by all five and the catalogue reader does not call it
        _, rows, _ = self.read_write()
        self.assertEqual(rows["drug a"]["miss_all5_no_catalog_call"], "1")
        self.assertEqual(rows["drug h"]["miss_all5_no_catalog_call"], "0")

    def test_misses_in_a_null_fit_fold(self):
        # drug a: gbt's fold 1 is null; of the all-five misses 0, 1 and 2 only genome 1 is in fold 1
        out, rows, _ = self.read_write()
        self.assertEqual((rows["drug a"]["miss_all5_in_null_fold"], rows["drug a"]["null_fold_pairs"]), ("1", "1"))
        self.assertIn("check 10 null folds: 0 learner rows", out)

    def test_a_null_fold_that_differs_from_no_call_is_counted_and_the_read_goes_on(self):
        r = Report(COLI, "drug a", 10, 5)
        r.no_call["lasso"] = {"resistant": 2, "susceptible": 0}
        out, _, _ = self.read_write([r])
        self.assertIn("check 10 null folds: 1 learner rows", out)
        self.assertIn("null folds differ: Synthia coli with drug a: lasso", out)


class Skips(Case):
    def test_no_lineage_split_and_one_learner_are_skipped_and_counted(self):
        out, rows, _ = self.read_write()
        self.assertIn("check 5 skipped: 1 with no lineage split, 1 with fewer than 2 learners; 2 skipped and 6 used "
                      "of 8 found: pass", out)
        self.assertNotIn("drug e", rows)
        self.assertNotIn("drug f", rows)

    def test_a_report_with_another_learner_set_is_listed_and_kept_out_of_the_pool(self):
        out, rows, pooled = self.read_write()
        self.assertIn("1 reports differ", out)
        self.assertIn("differs: Synthia aureus with drug g: scm,cart,rf,gbt", out)
        g = rows["drug g"]
        self.assertEqual((g["learner_count"], g["miss_4"], g["miss_5"], g["miss_all5"], g["miss_lasso"]),
                         ("4", "1", "", "", ""))
        self.assertNotIn(AUREUS, [p["species"] for p in pooled])


class Stops(Case):
    def test_a_score_list_of_another_length_stops_the_read(self):
        r = Report(COLI, "drug a", 10, 5)
        r.scores["rf"] = r.scores["rf"][:-1]
        self.stops([r], "check 2 lengths: 1 reports")

    def test_a_score_outside_0_to_1_stops_the_read(self):
        r = Report(COLI, "drug a", 10, 5)
        r.scores["gbt"][3] = 1.5
        self.stops([r], "check 2 lengths: 1 reports")

    def test_counts_that_differ_from_the_samples_stop_the_read(self):
        r = Report(COLI, "drug a", 10, 5)
        r.counts["resistant"] = 11
        self.stops([r], "check 3 counts: 1 reports")

    def test_stored_metrics_that_differ_from_the_recount_stop_the_read(self):
        r = Report(COLI, "drug a", 10, 5, misses={n: (0,) for n in FIVE})
        r.metric_change["cart"] = {"fn": 2, "tp": 8}
        self.stops([r], "check 7 stored metrics: 1 learner rows")

    def test_catalogue_calls_that_differ_from_catalog_resistant_stop_the_read(self):
        r = Report(COLI, "drug a", 10, 5)
        r.catalog_scores = [0.0] + [1.0 if x else 0.0 for x in r.catalog[1:]]
        self.stops([r], "check 9 catalogue: 1 reports")

    def test_a_file_that_does_not_parse_stops_the_read(self):
        d = lay_out(self.top, [Report(COLI, "drug a", 10, 5)])
        with open(os.path.join(d, "broken.json"), "w", encoding="utf-8") as f:
            f.write("{\"species\": ")
        rc, out = run("read", d, "--out", os.path.join(self.top, "reading"), "--expect", 2)
        self.assertEqual(rc, 1)
        self.assertIn("check 1 parse: 2 files found, 1 read, 2 expected: STOP", out)

    def test_another_file_count_than_expected_stops_the_read(self):
        self.stops([Report(COLI, "drug a", 10, 5)], "check 1 parse: 1 files found, 1 read, 274 expected: STOP",
                   expect=274)

    def test_two_lineage_splits_stop_the_read(self):
        r = Report(COLI, "drug a", 10, 5)
        r.schemes = ("lineage", "lineage")
        self.stops([r], "check 5 skipped")

    def test_a_pair_on_two_reports_stops_the_read(self):
        a, b = Report(COLI, "drug a", 10, 5, tag=1), Report(COLI, "drug a", 10, 5, tag=2)
        d = os.path.join(self.top, "reports")
        os.makedirs(d)
        for name, r in (("one.json", a), ("two.json", b)):
            with open(os.path.join(d, name), "w", encoding="utf-8") as f:
                json.dump(r.doc(), f)
        rc, out = run("read", d, "--out", os.path.join(self.top, "reading"), "--expect", 2)
        self.assertEqual(rc, 1)
        self.assertIn("check 8 pairs: 1 species and antibiotic pairs on two reports: STOP", out)


class Pooled(Case):
    def test_the_pooled_row_is_worked_by_hand(self):
        # pooled: drugs a, b, c and h (10 or more resistant, all five); d has 4 resistant, g four learners
        _, _, pooled = self.read_write()
        p = pooled[0]
        self.assertEqual((p["scope"], p["datasets"], p["resistant"], p["miss_all5"]), ("all", "4", "42", "4"))
        self.assertEqual((p["datasets_with_ratio"], p["datasets_expected_zero"], p["datasets_ratio_over_1"]),
                         ("3", "1", "2"))
        # ratios 72 (a), 0 (b) and 10,000 (h); c's yardstick is 0, as scm misses none
        self.assertEqual((p["ratio_median"], p["ratio_q1"], p["ratio_q3"]),
                         ("72.000000000000", "36.000000000000", "5036.000000000000"))
        self.assertEqual(p["expected_all5"], "0.041866666667")
        # shares by all five 0.25, 0, 0, 0.1; by 1 to 4: 0.25, 0.5, 0.3, 0
        self.assertEqual((p["share_all5_median"], p["share_some_median"]), ("0.050000000000", "0.275000000000"))

    def test_the_median_weights_each_report_once(self):
        # drug b holds 40 genomes and a ratio of 0; weighted by genomes the median would be 0
        _, _, pooled = self.read_write()
        self.assertEqual(pooled[0]["ratio_median"], "72.000000000000")

    def test_a_report_with_a_yardstick_of_0_is_reported_apart(self):
        _, rows, pooled = self.read_write()
        self.assertEqual((rows["drug c"]["expected_all5"], rows["drug c"]["ratio_all5"]), ("0.000000000000", ""))
        self.assertEqual((pooled[0]["datasets"], pooled[0]["datasets_expected_zero"]), ("4", "1"))

    def test_a_report_under_10_resistant_stays_out_of_the_pool(self):
        _, rows, pooled = self.read_write()
        self.assertEqual(rows["drug d"]["ratio_all5"], "256.000000000000")
        self.assertEqual(pooled[0]["datasets"], "4")
        self.assertEqual(pooled[0]["resistant"], str(12 + 10 + 10 + 10))

    def test_one_row_per_species_beside_the_whole(self):
        _, _, pooled = self.read_write()
        self.assertEqual([(p["scope"], p["species"]) for p in pooled], [("all", ""), ("species", COLI)])
        self.assertEqual(pooled[1]["datasets"], pooled[0]["datasets"])

    def test_the_reading_rule_bands(self):
        self.assertEqual([lo.reading_rule(x) for x in (72.0, 10.0, 9.99, 1.01, 1.0, 0.0)],
                         ["stack", "stack", "above", "above", "none shown", "none shown"])


class Ids(Case):
    def test_no_id_is_held_after_parsing(self):
        d = lay_out(self.top, [Report(COLI, "drug a", 10, 5)])
        doc = lo.load_report(os.path.join(d, os.listdir(d)[0]))

        def keys(v):
            if isinstance(v, dict):
                for k, x in v.items():
                    yield k
                    yield from keys(x)
            elif isinstance(v, list):
                for x in v:
                    yield from keys(x)
        self.assertNotIn("id", set(keys(doc)))
        self.assertEqual(len(doc["samples"]), 15)

    def test_no_genome_id_stands_in_any_written_file(self):
        reports = standard()
        out, _, _ = self.read_write(reports)
        self.assertIn("check 6 ids: 5 files scanned, 0 id-shaped tokens or barred keys: pass", out)
        written = lo.files_under(os.path.join(self.top, "reading")) + lo.files_under(os.path.join(self.top, "record"))
        self.assertEqual(len(written), 5)
        for path in written:
            with open(path, encoding="utf-8") as f:
                text = f.read()
            for r in reports:
                for g in r.ids:
                    self.assertNotIn(g, text, os.path.basename(path))

    def test_an_id_shaped_token_in_the_reading_stops_the_write(self):
        rc, out = self.read([Report(COLI, "drug a", 10, 5)])
        self.assertEqual(rc, 0, out)
        with open(os.path.join(self.top, "reading", "note.txt"), "w", encoding="utf-8") as f:
            f.write("a genome 990001.11 here\n")
        rc, out = run("write", os.path.join(self.top, "reading"), "--out", os.path.join(self.top, "record"))
        self.assertEqual(rc, 1, out)
        self.assertIn("check 6 ids: 1 id-shaped tokens or barred keys: STOP", out)
        self.assertFalse(os.path.exists(os.path.join(self.top, "record")))
        self.assertEqual([n for n in os.listdir(self.top) if n.startswith(".overlap-")], [])

    def test_the_id_shape(self):
        hits = lo.ID_SHAPE.findall("1280.12345 470.1 72.000000000000 0.5 0.041666666667 a1773.1 1773.10019x")
        self.assertEqual(hits, ["1280.12345", "470.1"])


class Record(Case):
    def test_the_same_reading_writes_the_same_bytes(self):
        self.read_write()
        rc, out = run("write", os.path.join(self.top, "reading"), "--out", os.path.join(self.top, "again"))
        self.assertEqual(rc, 0, out)
        for name in lo.RECORD_FILES + ("SHA256SUMS",):
            with open(os.path.join(self.top, "record", name), "rb") as f1, \
                    open(os.path.join(self.top, "again", name), "rb") as f2:
                self.assertEqual(f1.read(), f2.read(), name)

    def test_checksums_and_columns(self):
        self.read_write()
        rec = os.path.join(self.top, "record")
        self.assertEqual(sorted(os.listdir(rec)), sorted(lo.RECORD_FILES + ("SHA256SUMS",)))
        with open(os.path.join(rec, "SHA256SUMS"), encoding="utf-8") as f:
            for line in f.read().splitlines():
                digest, name = line.split("  ")
                self.assertEqual(lo.sha256(os.path.join(rec, name)), digest)
        with open(os.path.join(rec, "datasets.tsv"), encoding="utf-8") as f:
            self.assertEqual(f.readline().rstrip("\n").split("\t"), lo.DATASET_COLUMNS)
        with open(os.path.join(rec, "pooled.tsv"), encoding="utf-8") as f:
            self.assertEqual(f.readline().rstrip("\n").split("\t"), lo.POOLED_COLUMNS)

    def test_write_refuses_an_out_that_holds_files(self):
        self.read_write()
        rc, out = run("write", os.path.join(self.top, "reading"), "--out", os.path.join(self.top, "record"))
        self.assertEqual(rc, 1)
        self.assertIn("exists and is not an empty directory", out)

    def test_gates_hold_every_check_and_the_rule(self):
        self.read_write()
        with open(os.path.join(self.top, "record", "gates.txt"), encoding="utf-8") as f:
            lines = f.read().splitlines()
        heads = [line.split(" ", 2)[1] for line in lines if line.startswith("check ")]
        self.assertEqual(heads, [str(n) for n in range(1, 11)])
        self.assertEqual(lines[-1], "reading rule: stack")

    def test_the_page_tables_name_species_and_antibiotic_in_full(self):
        self.read_write()
        t = lo.tables(os.path.join(self.top, "record"))
        self.assertEqual(sorted(t), ["datasets", "overlap", "species", "yardstick"])
        self.assertIn("| *Synthia coli* | drug a | 20 | 12 | 3 | 1 | 1 | 0.04 | 72.0 |", t["datasets"])
        self.assertIn("| *Synthia aureus* | drug g | 15 | 10 | | | | | |", t["datasets"])
        self.assertIn("| *Synthia coli* | drug c | 15 | 10 | 0 | 0 | 0 | 0.00 | none |", t["datasets"])
        self.assertIn("| median ratio | 72.0 |", t["yardstick"])
        self.assertIn("| 5 | 4 | 1 |", t["overlap"])

    def test_plan_prints_the_fixed_rules(self):
        rc, out = run("plan")
        self.assertEqual(rc, 0)
        p = json.loads(out)
        self.assertEqual((p["cutoff"], p["scheme"], p["separate"], p["five"], p["min_resistant"], p["stack"]),
                         (0.5, "lineage", ["catalog"], list(FIVE), 10, 10))
        self.assertEqual(p["datasets_columns"], lo.DATASET_COLUMNS)


if __name__ == "__main__":
    unittest.main()
