"""Tests of scripts/refit/refit_check.py, the gate that a refit report is the
release's rule, on SYNTHETIC reports built here.

Run as CI runs the python suites:
    python3 -s -P -B -W error -m unittest discover -s tests/refit -t . -v
"""

import copy
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CHECK = os.path.join(ROOT, "scripts", "refit", "refit_check.py")
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
FOLDS = 5
N = 10
D = "synthetic_species__drug_a"


def release_report():
    """A SYNTHETIC report in the release's shape: the keys refit_check reads."""
    samples = [OrderedDict([("id", f"g{i}"), ("resistant", i % 2 == 0), ("lineage", i % 3), ("from_mic", False),
                            ("mic", None), ("groups", [i % 4]), ("catalog_resistant", i % 2 == 0),
                            ("country", "SYNTHETIC"), ("year", 2000 + i), ("host", None), ("mlst", None)])
               for i in range(N)]

    def learner(name, k):
        return OrderedDict([("learner", name), ("seconds", 0.5 + k),
                            ("metrics", OrderedDict([("balanced_accuracy", 0.8 - 0.01 * k), ("recall", 0.75)])),
                            ("fold_balanced_accuracy", [0.7 + 0.01 * f for f in range(FOLDS)]),
                            ("fold_complexity", [3] * FOLDS), ("scores", [0.1 * (i % 10) for i in range(N)]),
                            ("roc", [[0.0, 0.0], [1.0, 1.0]]), ("conformal", None), ("no_call", 0)])

    splits = [OrderedDict([("scheme", scheme), ("folds", FOLDS), ("n_groups", 3),
                           ("test_fold", [i % FOLDS for i in range(N)]),
                           ("learners", [learner(n, k) for k, n in enumerate(LEARNERS + ("catalog",))])])
              for scheme in ("random", "lineage")]
    return OrderedDict([("schema", 2), ("evaluation", OrderedDict([("revision", 14)])),
                        ("input_fingerprint", "SYNTHETIC-fingerprint"), ("id", D),
                        ("species_slug", "synthetic_species"), ("antibiotic", "drug a"),
                        ("counts", OrderedDict([("genomes", N), ("lineages", 3), ("min_support", 3)])),
                        ("seconds", 12.5), ("samples", samples), ("splits", splits),
                        ("models", OrderedDict([("scm", None)])), ("features", [card(i) for i in range(4)])])


def card(i):
    """A SYNTHETIC feature card in the release's key order."""
    locus = OrderedDict([("genome_id", f"g{i}"), ("contig", f"SYNTHETIC-contig-{i}"), ("start", 100 * i),
                         ("end", 100 * i + 31), ("kmers_placed", 1), ("gene", None), ("product", "SYNTHETIC"),
                         ("patric_id", f"SYNTHETIC-peg-{i}"), ("gene_strand", "+"), ("codons", None),
                         ("nearest", None), ("context", "SYNTHETIC")])
    return OrderedDict([("feature", 10 + i), ("n_kmers", 1 + i), ("kmer", f"SYNTHETICKMER{i}"),
                        ("catalog_families", ["SYNTHETIC-family"] if i == 0 else []),
                        ("catalog_products", []), ("locus", locus), ("n_resistant", 4 + i),
                        ("n_susceptible", 6 - i), ("bits", "SYNTHETIC-bits-%d" % i)])


def refit_of(release):
    """The release report as a refit writes it: other seconds, and fold
    models, rule groups and a selection added."""
    r = copy.deepcopy(release)
    r["seconds"] = 99.0
    for s in r["splits"]:
        for lr in s["learners"]:
            lr["seconds"] = lr["seconds"] + 7
    for s in r["samples"]:
        s["rule_groups"] = []
    r["selection"] = OrderedDict([("rule", "SYNTHETIC")])
    r["fold_models"] = OrderedDict([("features", []), ("splits", [
        OrderedDict([("scheme", s["scheme"]), *((n, [None] * FOLDS) for n in LEARNERS)]) for s in r["splits"]])])
    return r


class RefitCheck(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.refit = os.path.join(self.tmp.name, "refit")
        self.release = os.path.join(self.tmp.name, "release")
        self.out = os.path.join(self.tmp.name, "check")
        os.makedirs(self.refit)
        os.makedirs(os.path.join(self.release, "datasets"))
        self.base = release_report()
        self.write(self.release + "/datasets", self.base)

    def write(self, directory, report, name=D):
        with open(os.path.join(directory, name + ".json"), "w", encoding="utf-8") as f:
            json.dump(report, f)

    def run_check(self, refit=None, listed=D):
        if refit is not None:
            self.write(self.refit, refit)
        p = subprocess.run([sys.executable, "-I", "-B", CHECK, self.refit, self.release, listed, self.out],
                           capture_output=True, text=True)
        with open(os.path.join(self.out, "datasets.tsv"), encoding="utf-8", newline="") as f:
            ds = list(csv.DictReader(f, delimiter="\t"))
        with open(os.path.join(self.out, "learners.tsv"), encoding="utf-8", newline="") as f:
            ln = list(csv.DictReader(f, delimiter="\t"))
        return p, ds, ln

    def test_a_refit_with_only_the_allowed_differences_is_equal(self):
        p, ds, ln = self.run_check(refit_of(self.base))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual((ds[0]["fitted_fields"], ds[0]["fold_models"], ds[0]["text_equal_less_allowed"]),
                         ("equal", "complete", "yes"))
        self.assertEqual(ds[0]["samples_with_rule_groups"], str(N))
        self.assertEqual(ds[0]["cards_with_other_locus"], "0")
        # 2 splits x (5 learners and the catalog row).
        self.assertEqual(len(ln), 12)
        self.assertTrue(all(r["verdict"] == "equal" for r in ln))

    def test_a_score_raised_by_a_hundredth_names_the_learner(self):
        r = refit_of(self.base)
        lineage = next(s for s in r["splits"] if s["scheme"] == "lineage")
        rf = next(x for x in lineage["learners"] if x["learner"] == "rf")
        rf["scores"][3] += 0.01
        p, ds, ln = self.run_check(r)
        self.assertEqual(p.returncode, 1)
        bad = [x for x in ln if x["verdict"] != "equal"]
        self.assertEqual([(x["split"], x["learner"], x["fields_that_differ"]) for x in bad], [("lineage", "rf", "1")])
        self.assertAlmostEqual(float(bad[0]["largest_abs_difference"]), 0.01, places=9)
        self.assertIn("differs: " + D + " lineage rf, 1 fields", p.stdout)

    def test_fold_models_short_of_a_fold_for_any_learner_fail(self):
        for name in LEARNERS:
            r = refit_of(self.base)
            r["fold_models"]["splits"][1][name] = [None] * (FOLDS - 1)
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fold_models"]), (1, "incomplete"), name)
        r = refit_of(self.base)
        del r["fold_models"]["splits"][0]
        p, ds, _ = self.run_check(r)
        self.assertEqual((p.returncode, ds[0]["fold_models"]), (1, "incomplete"))

    def test_no_fold_models_fails(self):
        r = refit_of(self.base)
        del r["fold_models"]
        p, ds, _ = self.run_check(r)
        self.assertEqual((p.returncode, ds[0]["fold_models"]), (1, "no"))

    def test_another_new_key_fails(self):
        r = refit_of(self.base)
        r["SYNTHETIC_extra"] = 1
        p, ds, _ = self.run_check(r)
        self.assertEqual((p.returncode, ds[0]["unexpected_new_keys"]), (1, "SYNTHETIC_extra"))

    def test_a_missing_report_fails(self):
        p, ds, _ = self.run_check(None)
        self.assertEqual((p.returncode, ds[0]["fitted_fields"], ds[0]["cards_with_other_locus"]), (1, "no report", ""))

    def test_cards_whose_locus_and_kmer_differ_pass_and_are_counted(self):
        # The refetched genome file moved where two patterns are shown.
        r = refit_of(self.base)
        r["features"][1]["kmer"] = "SYNTHETICKMER-other"
        r["features"][1]["locus"]["contig"] = "SYNTHETIC-contig-other"
        r["features"][1]["locus"]["start"] = 7
        r["features"][3]["locus"]["patric_id"] = "SYNTHETIC-peg-other"
        p, ds, ln = self.run_check(r)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual((ds[0]["fitted_fields"], ds[0]["fields_that_differ"], ds[0]["text_equal_less_allowed"],
                          ds[0]["cards_with_other_locus"]), ("equal", "-", "yes", "2"))
        self.assertIn("datasets with a card whose kmer or locus differs 1, such cards 2", p.stdout)
        self.assertTrue(all(x["verdict"] == "equal" for x in ln))

    def test_a_card_that_differs_in_another_field_fails(self):
        for field, value in (("bits", "SYNTHETIC-bits-other"), ("catalog_families", []), ("n_resistant", 99),
                             ("feature", 77), ("catalog_products", ["SYNTHETIC"])):
            r = refit_of(self.base)
            r["features"][0][field] = value
            r["features"][2]["locus"]["end"] = 1
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fitted_fields"], ds[0]["fields_that_differ"],
                              ds[0]["text_equal_less_allowed"], ds[0]["cards_with_other_locus"]),
                             (1, "differs", "features", "no", "1"), field)
            self.assertIn("report fields differ: " + D + ": features", p.stdout, field)

    def test_a_missing_or_an_added_card_fails(self):
        for change in ("missing", "added"):
            r = refit_of(self.base)
            if change == "missing":
                del r["features"][2]
            else:
                r["features"].append(card(9))
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fields_that_differ"], ds[0]["text_equal_less_allowed"]),
                             (1, "features", "no"), change)

    def test_two_cards_swapped_fail(self):
        r = refit_of(self.base)
        r["features"][0], r["features"][1] = r["features"][1], r["features"][0]
        p, ds, _ = self.run_check(r)
        self.assertEqual((p.returncode, ds[0]["fields_that_differ"], ds[0]["text_equal_less_allowed"]),
                         (1, "features", "no"))

    def test_a_card_without_its_kmer_or_locus_key_fails(self):
        for key in ("kmer", "locus"):
            r = refit_of(self.base)
            del r["features"][1][key]
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fields_that_differ"], ds[0]["text_equal_less_allowed"]),
                             (1, "features", "no"), key)

    def test_a_locus_or_kmer_key_outside_features_that_differs_fails(self):
        base = release_report()
        base["models"]["scm"] = OrderedDict([("kmer", "SYNTHETICKMER0"), ("locus", card(0)["locus"])])
        self.write(self.release + "/datasets", base)
        for key in ("kmer", "locus"):
            r = refit_of(base)
            if key == "kmer":
                r["models"]["scm"]["kmer"] = "SYNTHETICKMER-other"
            else:
                r["models"]["scm"]["locus"]["contig"] = "SYNTHETIC-contig-other"
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fields_that_differ"], ds[0]["text_equal_less_allowed"],
                              ds[0]["cards_with_other_locus"]), (1, "models", "no", "0"), key)

    def test_descriptive_fields_are_counted_apart_and_pass(self):
        r = refit_of(self.base)
        for s in r["samples"][:4]:
            s["country"] = None
        r["samples"][5]["mlst"] = "SYNTHETIC"
        p, ds, _ = self.run_check(r)
        self.assertEqual(p.returncode, 0, p.stdout)
        self.assertEqual((ds[0]["fitted_fields"], ds[0]["samples_with_other_descriptive_fields"]), ("equal", "5"))

    def test_a_fitted_sample_field_or_report_field_fails(self):
        for change in ("lineage", "fingerprint", "test_fold"):
            r = refit_of(self.base)
            if change == "lineage":
                r["samples"][2]["lineage"] = 9
                want = "samples"
            elif change == "fingerprint":
                r["input_fingerprint"] = "SYNTHETIC-other"
                want = "input_fingerprint"
            else:
                r["splits"][1]["test_fold"][0] = 4
                want = "splits"
            p, ds, _ = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fields_that_differ"]), (1, want), change)

    def test_key_order_counts_in_the_text_check(self):
        r = refit_of(self.base)
        moved = OrderedDict([("id", r["id"])] + [(k, v) for k, v in r.items() if k != "id"])
        p, ds, _ = self.run_check(moved)
        self.assertEqual((p.returncode, ds[0]["fitted_fields"], ds[0]["text_equal_less_allowed"]), (1, "equal", "no"))

    def test_a_dropped_null_field_or_a_new_sample_field_fails_by_the_text(self):
        # The field checks read a missing key as null and compare a sample on
        # its named fields only: the text check is what fails these two.
        base = release_report()
        base["unsupported_assays"] = None
        self.write(self.release + "/datasets", base)
        for change in ("null_dropped", "sample_field"):
            r = refit_of(base)
            if change == "null_dropped":
                del r["unsupported_assays"]
            else:
                r["samples"][3]["SYNTHETIC_study"] = "SYNTHETIC"
            p, ds, ln = self.run_check(r)
            self.assertEqual((p.returncode, ds[0]["fitted_fields"], ds[0]["unexpected_new_keys"],
                              ds[0]["text_equal_less_allowed"]), (1, "equal", "-", "no"), change)
            self.assertTrue(all(x["verdict"] == "equal" for x in ln), change)

    def test_fold_models_in_another_split_order_fail(self):
        r = refit_of(self.base)
        r["fold_models"]["splits"].reverse()
        p, ds, _ = self.run_check(r)
        self.assertEqual((p.returncode, ds[0]["fold_models"]), (1, "incomplete"))

    def test_the_tables_hold_no_sample_id_or_value(self):
        r = refit_of(self.base)
        r["samples"][1]["lineage"] = 7
        r["features"][1]["kmer"] = "SYNTHETICKMER-other"
        r["features"][2]["bits"] = "SYNTHETIC-bits-other"
        p, _, _ = self.run_check(r)
        texts = [("stdout", p.stdout)]
        for name in ("datasets.tsv", "learners.tsv"):
            with open(os.path.join(self.out, name), encoding="utf-8") as f:
                texts.append((name, f.read()))
        for name, text in texts:
            for value in ("g1", "SYNTHETIC-fingerprint", "SYNTHETICKMER", "SYNTHETIC-contig", "SYNTHETIC-bits"):
                self.assertNotIn(value, text, name)


if __name__ == "__main__":
    unittest.main()
