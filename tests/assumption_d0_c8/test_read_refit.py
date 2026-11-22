"""SYNTHETIC tests of the refit's reader (scripts/assumption-d0-c8/read_refit.py)
and of `c8.py combine`, on made-up species directories (synth_refit.py).

Run as CI runs the python suites:
    python3 -s -P -B -W error -m unittest discover -s tests/assumption_d0_c8 -t . -v
"""

import contextlib
import csv
import io
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
SCRIPTS = os.path.join(ROOT, "scripts", "assumption-d0-c8")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import c8  # noqa: E402
import foldread as fr  # noqa: E402
import read_refit as rr  # noqa: E402
import synth_refit as S  # noqa: E402


def table(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def edit(path, old, new, count=1):
    with open(path, encoding="utf-8") as f:
        text = f.read()
    if old not in text:
        raise AssertionError(f"{old!r} not in {path}")
    with open(path, "w", encoding="utf-8") as f:
        f.write(text.replace(old, new, count))


class Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.paths = S.build(self.tmp.name)
        self.out = os.path.join(self.tmp.name, "out")

    def species(self, slug):
        return os.path.join(self.tmp.name, "back", slug)

    def run_reader(self, **over):
        err, out = io.StringIO(), io.StringIO()
        with contextlib.redirect_stderr(err), contextlib.redirect_stdout(out):
            code = rr.main(S.args(self.paths, self.out, **over))
        return code, out.getvalue(), err.getvalue()

    def stops(self, pattern, **over):
        code, _, err = self.run_reader(**over)
        self.assertEqual(code, 2, err)
        self.assertRegex(err, pattern)
        # A stopped reading writes nothing.
        self.assertFalse(os.path.exists(self.out) and os.listdir(self.out), "a stopped reading wrote files")
        return err

    def page(self, name):
        return table(os.path.join(self.out, name))


class Clean(Base):
    def test_reads_clean_and_names_in_full(self):
        code, out, err = self.run_reader()
        self.assertEqual(code, 0, err)
        self.assertIn("no single cause holds 0.40 of 6 miss pairs: C4, C3, C5", out)
        self.assertIn("C4, C3, C5 (the plan's default)", out)
        sp = {r["species"]: r for r in self.page("page-species.tsv")}
        self.assertEqual(sorted(sp), ["Synthia alpha", "Synthia beta", "Synthia gamma"])
        self.assertEqual(sp["Synthia alpha"]["antibiotics"], "drugamycin; fakecillin; mockazole/trimethoprim")
        self.assertEqual(sp["Synthia gamma"]["d0_antibiotics"], "none")
        self.assertEqual((sp["Synthia alpha"]["breakpoints_set_to_the_release"],
                          sp["Synthia alpha"]["feature_lists_equal_to_master"]), ("1", "0"))
        c8rows = self.page("page-c8-datasets.tsv")
        self.assertEqual(len(c8rows), 6 * 5 * 2)
        self.assertIn("nalidixic acid", {r["antibiotic"] for r in c8rows})
        d0rows = self.page("page-d0-datasets.tsv")
        self.assertEqual(len(d0rows), 6)
        self.assertEqual({(r["species"], r["antibiotic"]) for r in d0rows},
                         {("Synthia alpha", "drugamycin"), ("Synthia alpha", "fakecillin"),
                          ("Synthia beta", "drugamycin")})
        for name in ("page-d0-learners.tsv", "page-d0-datasets.tsv"):
            with open(os.path.join(self.out, name), encoding="utf-8") as f:
                head = f.readline()
            self.assertIn("amrfinderplus_program", head)
            self.assertNotIn("known", head)

    def test_d0_learner_table_reuses_d0_combine(self):
        self.assertEqual(self.run_reader()[0], 0)
        by = {r["learner"]: r for r in self.page("page-d0-learners.tsv")}
        g = "share_of_gap_to_the_amrfinderplus_program_"
        self.assertEqual((by["scm"]["miss_pairs"], by["scm"]["datasets_absent"], by["scm"][g + "absent"]),
                         ("2", "2", "1.000000"))
        self.assertEqual((by["cart"]["miss_pairs"], by["cart"]["pairs_with_no_missed_fold"], by["cart"][g + "absent"]),
                         ("1", "1", "0.000000"))
        rf = by["rf"]
        self.assertEqual((rf["datasets_removed"], rf["datasets_unused"], rf["datasets_lineage"], rf["ties"]),
                         ("1", "1", "1", "1"))
        self.assertEqual((rf[g + "removed"], rf[g + "unused"], rf[g + "lineage"]), ("0.250000", "0.250000", "0.500000"))
        self.assertEqual((by["gbt"]["miss_pairs"], by["gbt"][g + "absent"]), ("0", "NA"))
        every = by["all five (dataset and learner pairs)"]
        self.assertEqual([every[f"datasets_{c}"] for c in ("absent", "removed", "lineage", "chosen", "unused")],
                         ["2", "1", "1", "0", "2"])
        self.assertEqual((every["miss_pairs"], every[g + "absent"]), ("6", "0.318182"))
        order = {r["cause"]: r for r in self.page("page-d0-order.tsv")}
        self.assertEqual((order["absent"]["largest_of_the_four"], order["absent"]["sets_the_order"]), ("yes", "no"))
        self.assertEqual(order["chosen"]["sets_the_order"], "never (entry D0 causes)")
        # d0.py combine's own tables are kept beside the page's.
        self.assertTrue(os.path.isfile(os.path.join(self.out, "d0", "d0-summary.tsv")))

    def test_order_line_when_a_cause_sets_it(self):
        alpha = self.species("synthia_alpha")
        p = os.path.join(alpha, "results", "d0", "d0-datasets.tsv")
        rows = table(p)
        for r in rows:
            if r["dataset"] == "synthia_alpha__drugamycin" and r["learner"] == "lasso":
                r["cause"], r["folds_unused"], r["folds_absent"] = "absent", "0", "1"
        with open(p, "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(list(rows[0]))
            w.writerows([list(r.values()) for r in rows])
        S.seal(alpha)
        code, out, err = self.run_reader()
        self.assertEqual(code, 0, err)
        self.assertIn("absent holds 3 of 6 miss pairs (0.500): C3 first", out)
        self.assertIn("order of the input group: C3 first", out)
        shutil.rmtree(self.out)
        code, out, err = self.run_reader(order_rest="default-after")
        self.assertEqual(code, 0, err)
        self.assertIn("order of the input group: C3, then C4 and C5", out)


class SpeciesGates(Base):
    def test_a_missing_species_stops(self):
        self.stops(r"2 species, want 3", dirs=self.paths["dirs"][:2])

    def test_a_doubled_species_stops(self):
        twin = shutil.copytree(self.species("synthia_beta"), os.path.join(self.tmp.name, "back", "beta-again"))
        self.stops(r"species synthia_beta twice", dirs=self.paths["dirs"] + [twin], species="4")

    def test_no_done_or_a_failed_species_stops(self):
        os.remove(os.path.join(self.species("synthia_gamma"), "DONE"))
        self.stops(r"synthia_gamma: no DONE")
        S.write_species(os.path.join(self.tmp.name, "back"), "synthia_gamma")
        with open(os.path.join(self.species("synthia_gamma"), "FAILED"), "w", encoding="utf-8") as f:
            f.write("FAILED synthia_gamma at learn 2026-11-22T00:00:00Z\n")
        self.stops(r"FAILED is present \(FAILED synthia_gamma at learn")

    def test_a_wrong_code_sha256_stops(self):
        S.write_species(os.path.join(self.tmp.name, "back"), "synthia_beta", code="d1" * 32)
        self.stops(r"synthia_beta: DONE has code sha256 d1d1")

    def test_binary_txt_must_name_dones_code(self):
        edit(os.path.join(self.species("synthia_beta"), "binary.txt"), "code " + S.CODE, "code " + "e" * 64)
        S.seal(self.species("synthia_beta"))
        self.stops(r"binary.txt does not name DONE's binary and code")

    def test_a_file_changed_after_sha256sums_stops(self):
        edit(os.path.join(self.species("synthia_alpha"), "results", "c8", "c8-datasets.tsv"), "0.600000", "0.700000")
        self.stops(r"results/c8/c8-datasets.tsv is missing or differs from its sha256")

    def test_check_tables_must_read_equal(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "check", "datasets.tsv"), "synthia_alpha__fakecillin\tequal\t-",
             "synthia_alpha__fakecillin\tdiffers\tsamples")
        S.seal(a)
        self.stops(r"check/datasets.tsv: synthia_alpha__fakecillin does not read equal")
        S.write_species(os.path.join(self.tmp.name, "back"), "synthia_alpha")
        edit(os.path.join(a, "check", "learners.tsv"), "synthia_alpha__drugamycin\tlineage\trf\tequal\t0\t0",
             "synthia_alpha__drugamycin\tlineage\trf\tdiffers\t1\t0.002")
        S.seal(a)
        self.stops(r"check/learners.tsv: synthia_alpha__drugamycin lineage rf does not read equal")

    def with_card_column(self, slug, counts, name="cards_with_other_locus"):
        """Rewrite the species' check/datasets.tsv as the later check writes it, and seal."""
        a = self.species(slug)
        path = os.path.join(a, "check", "datasets.tsv")
        with open(path, encoding="utf-8") as f:
            lines = f.read().splitlines()
        self.assertEqual(len(lines) - 1, len(counts))
        with open(path, "w", encoding="utf-8") as f:
            f.write(lines[0] + "\t" + name + "\n")
            f.write("".join(f"{line}\t{n}\n" for line, n in zip(lines[1:], counts)))
        S.seal(a)

    def test_check_datasets_table_with_and_without_the_card_column(self):
        # One species from the later check (a card column), the others from the earlier one (none).
        self.with_card_column("synthia_alpha", ["0", "2", "5"])
        code, _, err = self.run_reader()
        self.assertEqual(code, 0, err)
        sp = {r["species"]: r for r in self.page("page-species.tsv")}
        self.assertEqual((sp["Synthia alpha"]["datasets_with_other_locus_cards"],
                          sp["Synthia alpha"]["cards_with_other_locus"]), ("2", "7"))
        self.assertEqual((sp["Synthia beta"]["datasets_with_other_locus_cards"],
                          sp["Synthia beta"]["cards_with_other_locus"]), ("0", "0"))

    def test_a_card_count_that_is_no_count_or_another_new_column_stops(self):
        self.with_card_column("synthia_alpha", ["0", "some", "5"])
        self.stops(r"check/datasets.tsv: synthia_alpha__\w+: a count is not a count")
        S.write_species(os.path.join(self.tmp.name, "back"), "synthia_alpha")
        self.with_card_column("synthia_alpha", ["0", "0", "0"], name="cards_with_other_bits")
        self.stops(r"check/datasets.tsv: its columns are not the run's")

    def test_a_check_row_missing_stops(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "check", "learners.tsv"), "synthia_alpha__drugamycin\trandom\tgbt\tequal\t0\t0\n", "")
        S.seal(a)
        self.stops(r"synthia_alpha__drugamycin random lacks one of the five learners")

    def test_release_specs_and_feature_lists_are_counted_and_gated_only_when_named(self):
        self.assertEqual(self.run_reader()[0], 0)
        shutil.rmtree(self.out)
        self.stops(r"1 breakpoints differ", equal_tables="datasets.tsv,learners.tsv,release-specs.tsv")
        self.stops(r"0 of 11 feature lists equal the master list",
                   equal_tables="datasets.tsv,learners.tsv,feature-lists.tsv")
        self.stops(r"--equal-tables must hold", equal_tables="datasets.tsv")


class DatasetGates(Base):
    def test_the_dataset_count_is_an_argument(self):
        self.stops(r"6 datasets, want 7", datasets="7")

    def test_marks_name_a_missing_dataset(self):
        S.write_marks(self.paths["marks"], extra=["synthia_gamma__zzzmycin"])
        self.stops(r"1 missing \(the first synthia_gamma__zzzmycin\)")

    def test_a_dataset_missing_from_c8_stops(self):
        b = self.species("synthia_beta")
        p = os.path.join(b, "results", "c8", "c8-datasets.tsv")
        with open(p, encoding="utf-8") as f:
            lines = [x for x in f if "nalidixic" not in x]
        with open(p, "w", encoding="utf-8") as f:
            f.writelines(lines)
        S.seal(b)
        self.stops(r"the C8 table does not hold the listed datasets exactly")

    def test_a_dataset_needs_its_name(self):
        edit(self.paths["names"], "synthia_beta__nalidixic_acid\tsynthia_beta\tnalidixic acid\n", "")
        self.stops(r"synthia_beta__nalidixic_acid: no species and antibiotic name")


class D0Gates(Base):
    def test_d0_counts_are_arguments(self):
        self.stops(r"D0: 3 datasets in 2 species, want 4 in 2", d0_datasets="4")
        self.stops(r"D0: 3 datasets in 2 species, want 3 in 3", d0_species="3")
        self.stops(r"d0.py combine: gate failed: 6 miss pairs, expected 7", d0_pairs="7")
        self.stops(r"D0: rf has 2 miss pairs, want 3", d0_learner_pairs="scm=2,cart=1,rf=3,gbt=0,lasso=1")
        self.assertEqual(self.run_reader(d0_learner_pairs="scm=2,cart=1,rf=2,gbt=0,lasso=1")[0], 0)

    def test_d0_list_must_be_the_miss_set_among_the_listed(self):
        b = self.species("synthia_beta")
        with open(os.path.join(b, "d0.txt"), "w", encoding="utf-8") as f:
            f.write("")
        shutil.rmtree(os.path.join(b, "results", "d0"))
        S.seal(b)
        self.stops(r"synthia_beta: d0.txt lists 0 datasets, the MISS set among its listed holds 1")

    def test_a_skipped_d0_with_d0_tables_stops(self):
        g = self.species("synthia_gamma")
        shutil.copytree(os.path.join(self.species("synthia_beta"), "results", "d0"), os.path.join(g, "results", "d0"))
        S.seal(g)
        self.stops(r"D0 was skipped and results/d0 holds tables")

    def test_a_d0_row_outside_the_records_miss_set_stops(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "results", "d0", "d0-datasets.tsv"),
             "synthia_alpha__fakecillin\tsynthia_alpha\tgbt\tfold models\t0.800000\t0.800000\t0.000000\tno",
             "synthia_alpha__fakecillin\tsynthia_alpha\tgbt\tfold models\t0.800000\t0.800000\t0.000000\tyes")
        S.seal(a)
        self.stops(r"synthia_alpha__fakecillin gbt: in the miss set by D0's table and not by the record")

    def test_a_cause_read_without_fold_models_stops(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "results", "d0", "d0-datasets.tsv"), "\tfold models\t", "\trelease report, no fold models\t")
        S.seal(a)
        self.stops(r"not read from fold models")

    def test_a_d0_dataset_twice_stops(self):
        a = self.species("synthia_alpha")
        p = os.path.join(a, "results", "d0", "d0-datasets.tsv")
        with open(p, encoding="utf-8") as f:
            lines = f.readlines()
        with open(p, "w", encoding="utf-8") as f:
            f.writelines(lines + lines[1:2])
        S.seal(a)
        self.stops(r"synthia_alpha__drugamycin scm twice")


class C8Gates(Base):
    def test_feature_text_is_refused_and_not_printed(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "results", "c8", "c8-datasets.tsv"), "0.600000", "ACGTTGCAACGTTGCA")
        S.seal(a)
        err = self.stops(r"not an overlap; a C8 table holds counts")
        self.assertNotIn("ACGT", err)
        S.write_species(os.path.join(self.tmp.name, "back"), "synthia_alpha")
        edit(os.path.join(a, "results", "c8", "c8-datasets.tsv"), "\tscm\t", "\tblaKPC\t")
        S.seal(a)
        err = self.stops(r"column learner: not a learner")
        self.assertNotIn("blaKPC", err)

    def test_pairs_that_folds_cannot_give_stop(self):
        a = self.species("synthia_alpha")
        edit(os.path.join(a, "results", "c8", "c8-datasets.tsv"), "\t5\t10\t0.600000", "\t5\t5\t0.600000")
        S.seal(a)
        self.stops(r"5 pairs cannot come from 5 folds")

    def test_a_split_that_did_not_run_holds_no_value(self):
        b = self.species("synthia_beta")
        edit(os.path.join(b, "results", "c8", "c8-datasets.tsv"), "\t0\t0\tNA\tNA\tNA\tNA\tNA\n",
             "\t0\t0\tNA\tNA\t0\t0\t0\n")
        S.seal(b)
        self.stops(r"random: a split that did not run holds a value")

    def test_a_c8_row_twice_stops(self):
        a = self.species("synthia_alpha")
        p = os.path.join(a, "results", "c8", "c8-datasets.tsv")
        with open(p, encoding="utf-8") as f:
            lines = f.readlines()
        with open(p, "w", encoding="utf-8") as f:
            f.writelines(lines[:2] + lines[1:])
        S.seal(a)
        self.stops(r"synthia_alpha__\w+ scm twice")

    def test_fewer_than_two_folds_that_chose(self):
        code, _, err = self.run_reader(c8_one_fold="zero")
        self.assertEqual(code, 0, err)
        s = {(r["learner"], r["scheme"]): r for r in self.page("page-c8-summary.tsv")}
        scm = s[("scm", "lineage")]
        self.assertEqual((scm["units"], scm["defined"], scm["undefined_one_fold_chose"], scm["median"]),
                         ("6", "6", "0", "0.350000"))
        self.assertEqual((scm["q1"], scm["q3"]), ("0.225000", "0.550000"))
        cart = s[("cart", "lineage")]
        self.assertEqual((cart["defined"], cart["undefined_no_fold_chose"]), ("5", "1"))
        rows = {(r["dataset"], r["learner"], r["scheme"]): r for r in self.page("page-c8-datasets.tsv")}
        one = rows[("synthia_gamma__fakecillin", "scm", "lineage")]
        self.assertEqual((one["folds_choosing"], one["median_jaccard"], one["defined"]), ("1", "0.000000", "yes"))
        shutil.rmtree(self.out)
        code, _, err = self.run_reader(c8_one_fold="undefined")
        self.assertEqual(code, 0, err)
        s = {(r["learner"], r["scheme"]): r for r in self.page("page-c8-summary.tsv")}
        scm = s[("scm", "lineage")]
        self.assertEqual((scm["units"], scm["defined"], scm["undefined_one_fold_chose"], scm["median"]),
                         ("6", "5", "1", "0.400000"))
        rows = {(r["dataset"], r["learner"], r["scheme"]): r for r in self.page("page-c8-datasets.tsv")}
        one = rows[("synthia_gamma__fakecillin", "scm", "lineage")]
        self.assertEqual((one["median_jaccard"], one["defined"]), ("NA", "one fold chose a feature"))

    def test_split_not_run_counts_apart_and_random_decides_nothing(self):
        self.assertEqual(self.run_reader()[0], 0)
        s = {(r["learner"], r["scheme"]): r for r in self.page("page-c8-summary.tsv")}
        self.assertEqual((s[("rf", "random")]["undefined_split_not_run"], s[("rf", "random")]["defined"]), ("1", "5"))
        self.assertIn("decides nothing", s[("rf", "random")]["reads"])
        self.assertIn("decides nothing alone", s[("rf", "lineage")]["reads"])
        scm = s[("scm", "lineage")]
        self.assertEqual((scm["distinct_summed"], scm["with_card_summed"],
                          scm["with_amrfinderplus_catalog_family_summed"], scm["catalog_family_share_of_carded"]),
                         ("21", "10", "5", "0.500000"))
        self.assertEqual((s[("rf", "lineage")]["pattern_median"], s[("lasso", "lineage")]["pattern_median"]),
                         ("NA", "0.500000"))

    def test_over_five_and_bands_are_arguments(self):
        self.assertEqual(self.run_reader(c8_bands="0.2,0.8")[0], 0)
        s = {(r["learner"], r["scheme"]): r for r in self.page("page-c8-summary.tsv")}
        five = s[(c8.OVER_FIVE, "lineage")]
        self.assertEqual((five["unit"], five["units"], five["defined"], five["undefined_no_fold_chose"]),
                         ("dataset and learner", "30", "29", "1"))
        scm = s[("scm", "lineage")]
        self.assertEqual((scm["under_band_low"], scm["over_band_high"]), ("1", "0"))
        shutil.rmtree(self.out)
        self.assertEqual(self.run_reader(c8_over_five="dataset-median")[0], 0)
        s = {(r["learner"], r["scheme"]): r for r in self.page("page-c8-summary.tsv")}
        five = s[(c8.OVER_FIVE, "lineage")]
        self.assertEqual((five["unit"], five["units"], five["defined"], five["undefined_no_learner_defined"],
                          five["median"]), ("dataset", "6", "6", "0", "0.500000"))
        self.assertEqual(s[(c8.OVER_FIVE, "random")]["undefined_no_learner_defined"], "1")
        shutil.rmtree(self.out)
        self.stops(r"--bands must be none or LOW,HIGH", c8_bands="0.2")

    def test_a_dataset_marked_m_in_c8_stops(self):
        b = self.species("synthia_beta")
        edit(os.path.join(b, "results", "c8", "c8-datasets.tsv"), "synthia_beta\tno\tscm", "synthia_beta\tyes\tscm")
        S.seal(b)
        self.stops(r"the C8 table holds a dataset marked M")


class Outputs(Base):
    def test_an_existing_output_stops_and_open_rules_have_no_default(self):
        self.assertEqual(self.run_reader()[0], 0)
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            self.assertEqual(rr.main(S.args(self.paths, self.out)), 2)
        self.assertIn("exists", err.getvalue())
        for opt in ("c8_one_fold", "c8_over_five", "c8_bands", "order_rest"):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                rr.main(S.args(self.paths, os.path.join(self.tmp.name, "o2"), **{opt: None}))

    def test_outputs_hold_no_genome_or_feature_text(self):
        self.assertEqual(self.run_reader()[0], 0)
        for dirpath, _, files in os.walk(self.out):
            for n in files:
                with open(os.path.join(dirpath, n), encoding="utf-8") as f:
                    text = f.read()
                self.assertNotRegex(text, r"[ACGT]{8}|bla[A-Z]", n)

    def test_the_reader_runs_without_numpy(self):
        code = ("import sys; sys.modules['numpy'] = None; sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[2]);"
                "import read_refit, synth_refit as S; import json;"
                "p = json.loads(sys.argv[3]);"
                "sys.exit(read_refit.main(S.args(p, sys.argv[4])))")
        import json
        r = subprocess.run([sys.executable, "-I", "-c", code, SCRIPTS, HERE, json.dumps(self.paths), self.out],
                           capture_output=True, text=True, check=False)
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("no single cause holds 0.40", r.stdout)


class C8Combine(unittest.TestCase):
    def test_choosing_inverts_the_kept_pairs(self):
        self.assertEqual([c8.choosing(5, p) for p in (0, 4, 7, 9, 10)], [[0], [1], [2], [3], [4, 5]])
        self.assertEqual(c8.choosing(5, 5), [])
        self.assertEqual(c8.choosing(2, 1), [1, 2])

    def test_two_folds_cannot_tell_one_fold_that_chose(self):
        s = {"folds": 2, "pairs": 1, "median_jaccard": 0.0, "median_jaccard_pattern": None, "distinct": 1,
             "with_card": 0, "with_family": 0}
        with self.assertRaisesRegex(c8.CombineError, "cannot tell one fold"):
            c8.define(s, "undefined")

    def test_species_name_in_full(self):
        self.assertEqual(c8.species_name("klebsiella_aerogenes"), "Klebsiella aerogenes")

    def test_combine_cli_writes_its_two_tables(self):
        with tempfile.TemporaryDirectory() as tmp:
            paths = S.build(tmp)
            dirs = [os.path.join(d, "results", "c8") for d in paths["dirs"]]
            out = os.path.join(tmp, "c8out")
            with contextlib.redirect_stdout(io.StringIO()):
                code = c8.main(["combine", "--out", out, "--one-fold", "zero", "--over-five", "pairs", "--bands",
                                "none", "--names", paths["names"], "--expect-datasets", "6"] + dirs)
            self.assertEqual(code, 0)
            self.assertEqual(len(table(os.path.join(out, "c8-combined-datasets.tsv"))), 60)
            self.assertEqual(len(table(os.path.join(out, "c8-combined-summary.tsv"))), 2 * (len(fr.LEARNERS) + 1))
            with contextlib.redirect_stderr(io.StringIO()) as err:
                code = c8.main(["combine", "--out", out, "--one-fold", "zero", "--over-five", "pairs", "--bands",
                                "none", "--names", paths["names"]] + dirs + dirs[:1])
            self.assertEqual(code, 2)
            self.assertIn("species synthia_alpha again", err.getvalue())


if __name__ == "__main__":
    unittest.main()
