"""The committed release summary outputs (reference/release-summary) and their README.

The suite recomputes nothing from a release. It checks that each output parses, names its release's revision and 274
datasets, that its counts partition as they must, and that the README lists the sha256 of each output. guidance.json must
read the same rows as uncertainty.json, sort them into the same screen classes, and count the no-call rows that
index-counts.json finds."""

import hashlib
import json
import pathlib
import re
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "release-summary"
RELEASES = {"rev12-274-20261105": 12, "rev14-274-20261112": 14}
OUTPUTS = ["uncertainty.json", "summary.json", "catalog.json", "index-counts.json", "guidance.json"]
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]


def load(release, name):
    return json.loads((REFERENCE / release / name).read_text(encoding="utf-8"))


class Reference(unittest.TestCase):
    def test_every_output_is_there_and_parses(self):
        for release in RELEASES:
            self.assertEqual(sorted(p.name for p in (REFERENCE / release).iterdir()), sorted(OUTPUTS), release)
            for name in OUTPUTS:
                self.assertIsInstance(load(release, name), dict, f"{release}/{name}")

    def test_revision_and_datasets(self):
        for release, revision in RELEASES.items():
            with self.subTest(release=release):
                for name in ("uncertainty.json", "summary.json", "index-counts.json", "guidance.json"):
                    doc = load(release, name)
                    self.assertEqual((doc["revision"], doc["datasets"]), (revision, 274), name)
                # catalog.json names no revision; it reads the same index as summary.json beside it
                catalog, summary = load(release, "catalog.json"), load(release, "summary.json")
                self.assertEqual(catalog["datasets"], 274)
                self.assertEqual(catalog["rule"]["with_an_entry"], summary["catalog"]["with_rule"])
                unc = load(release, "uncertainty.json")
                self.assertEqual((unc["species"], summary["species"]), (23, 23))
                self.assertEqual((unc["seed"], unc["reps"]), (20261111, 4000))
                self.assertEqual(unc["checks"]["intervals_files"], 274)
                self.assertEqual(unc["checks"]["problem_count"], 0, unc["checks"]["problems"])

    def test_error_limit_rows_partition(self):
        for release in RELEASES:
            bars = load(release, "uncertainty.json")["error_bars"]
            self.assertEqual((bars["vme_bar"], bars["me_bar"]), (0.015, 0.03))
            self.assertEqual(len(bars["rows"]), 2 * (len(LEARNERS) + 2))
            for key, row in bars["rows"].items():
                with self.subTest(release=release, row=key):
                    self.assertEqual(row["shown_to_meet_both"] + row["shown_to_miss_one"] + row["not_decided"], row["n"])
                    self.assertLessEqual(row["point_meets_both"], min(row["point_meets_vme"], row["point_meets_me"]))
                    self.assertLessEqual(max(row["shown_to_miss_vme"], row["shown_to_miss_me"]), row["shown_to_miss_one"])
                    self.assertLessEqual(row["shown_to_miss_one"], row["shown_to_miss_vme"] + row["shown_to_miss_me"])
                    # an upper bound within both limits means a point within both
                    self.assertLessEqual(row["shown_to_meet_both"], row["point_meets_both"])

    def test_intervals_are_ordered(self):
        for release in RELEASES:
            unc = load(release, "uncertainty.json")
            for block in ("medians", "gaps", "paired_lineage", "shares", "catalog"):
                for key, r in unc[block].items():
                    for name in ("ci95_datasets", "ci95_species"):
                        with self.subTest(release=release, stat=key, interval=name):
                            lo, hi = r[name]
                            self.assertLessEqual(lo, hi)

    def test_known_gene_bands_partition(self):
        for release in RELEASES:
            c = load(release, "catalog.json")
            b = {k: v["datasets"] for k, v in c["bands"].items()}
            with self.subTest(release=release):
                self.assertEqual(b["under 0.55"] + b["0.55 to 0.70"] + b["0.70 to 0.85"] + b["0.85 and up"],
                                 c["rule"]["with_an_entry"])
                self.assertEqual(b["all"], c["rule"]["with_an_entry"])
                self.assertEqual(b["0.55 to 0.70"] + b["0.70 to 0.85"] + b["0.85 and up"], b["0.55 and up"])
                # a constant call has balanced accuracy exactly 0.5
                self.assertLessEqual(c["rule"]["calls_nothing_resistant"] + c["rule"]["calls_everything_resistant"],
                                     c["rule"]["balanced_accuracy_exactly_half"])
                self.assertEqual(sum(c["at_chance_by"]["failure"].values()), b["under 0.55"])
                for field in ("species", "drug_class"):
                    self.assertEqual(sum(r["under_0.55"] for r in c["at_chance_by"][field]), b["under 0.55"])
                for name, band in c["bands"].items():
                    self.assertLessEqual(band["learner_ahead_lineage"], band["with_a_lineage_score"], name)
                    self.assertLessEqual(band["with_a_lineage_score"], band["datasets"], name)

    def test_summary_counts_partition(self):
        for release in RELEASES:
            s = load(release, "summary.json")
            with self.subTest(release=release):
                self.assertEqual(sum(s["best_learner_counts"].values()), s["datasets"])
                self.assertEqual(s["best"]["n_both"] + len(s["no_lineage"]), s["datasets"])
                self.assertEqual(sum(r["datasets"] for r in s["per_species"]), s["datasets"])
                self.assertEqual(len(s["per_species"]), s["species"])
                for split, r in s["catalog_vs_best"].items():
                    self.assertEqual(r["learner_ahead_over_2pt"] + r["catalog_within_2pt_or_ahead"], r["n"], split)
                    self.assertLessEqual(r["learner_ahead"] + r["catalog_ahead"], r["n"], split)
                for split, r in s["error_bars_best_learner"].items():
                    self.assertLessEqual(r["both"], min(r["vme_le_1.5"], r["me_le_3"]), split)
                sizes = {k: v for k, v in s["lineage_by_size"].items() if k != "cuts"}
                self.assertEqual(sum(n for _, n in sizes.values()), s["best"]["n_both"])

    def test_index_counts(self):
        for release, revision in RELEASES.items():
            c = load(release, "index-counts.json")
            with self.subTest(release=release):
                self.assertEqual(c["mismatches"], 0, c["examples"])
                self.assertEqual(c["files"], 274)
                self.assertEqual(c["rows"], c["learner_rows"] + c["rule_rows"])
                self.assertEqual(c["rates_compared"] + c["rates_null_in_file"], 2 * c["rows"])
                self.assertEqual(c["rows"], load(release, "uncertainty.json")["checks"]["rows"])
                if revision == 12:
                    self.assertEqual((c["flagged_rows"], c["flag_keys_seen"]), (0, {}))
                else:
                    self.assertEqual(c["flagged_rows"], 122)
                    self.assertEqual(c["flagged_with_no_call"], c["flagged_rows"])
                    self.assertEqual(c["unflagged_with_no_call"], 0)
                    self.assertEqual(c["flag_keys_seen"], {"split:true": 122})

    def test_guidance_criteria(self):
        for release in RELEASES:
            c = load(release, "guidance.json")["criteria"]
            with self.subTest(release=release):
                w = c["wording"]
                self.assertEqual((w["vme_upper_below"], w["vme_lower_below"], w["me_rate_below"], w["strict"]),
                                 (0.075, 0.015, 0.03, True))
                self.assertEqual((c["screen"]["vme_upper_at_most"], c["screen"]["me_upper_at_most"]), (0.015, 0.03))
                self.assertEqual(c["without_tuberculosis"], "mycobacterium_tuberculosis")
                # the guidance's Table 8 starts at 48 resistant organisms, with 0 discrepancies and (0.00, 7.40)
                self.assertEqual(c["table_8_first_resistant"], 48)

    def test_guidance_parts_add_up(self):
        for release in RELEASES:
            g = load(release, "guidance.json")
            self.assertEqual(len(g["rows"]), 2 * (len(LEARNERS) + 2))
            for key, r in g["rows"].items():
                with self.subTest(release=release, row=key):
                    self.assertLessEqual(r["meets_wording"], min(r["vme_upper"], r["vme_lower"], r["me"], r["vme"]))
                    self.assertLessEqual(r["vme"], min(r["vme_upper"], r["vme_lower"]))
                    self.assertGreaterEqual(r["vme"], r["vme_upper"] + r["vme_lower"] - r["n"])
                    self.assertLessEqual(max(r["vme_upper"], r["vme_lower"], r["me"]), r["n"])
                    self.assertEqual(sum(r["screen"].values()), r["n"])
                    t = r["without_tuberculosis"]
                    self.assertLessEqual(t["meets_wording"], min(t["n"], r["meets_wording"]))
                    self.assertLessEqual(r["meets_wording"] - t["meets_wording"], r["n"] - t["n"])
                    # with fewer than 48 resistant genomes no exact upper bound falls under 7.5%
                    self.assertLessEqual(r["vme_upper"], r["n"] - r["vme_n_under_table_8"])
                    self.assertLessEqual(r["fewer_than_every_genome"], r["n"])
                    # the screen asks more: an upper bound at or under 1.5% and 3% is under 7.5% and 1.5% as well
                    self.assertLessEqual(r["screen"]["shown_to_meet_both"], r["meets_wording"])

    def test_guidance_list_is_the_best_learners_count(self):
        for release in RELEASES:
            g = load(release, "guidance.json")
            for split, listed in g["best_meets_wording"].items():
                with self.subTest(release=release, split=split):
                    self.assertEqual(len(listed), g["rows"][f"best:{split}"]["meets_wording"])
                    self.assertEqual(len({e["dataset"] for e in listed}), len(listed))
                    for e in listed:
                        self.assertIn(e["learner"], LEARNERS)
                        self.assertLess(e["vme"]["hi"], 0.075, e["dataset"])
                        self.assertLess(e["vme"]["lo"], 0.015, e["dataset"])
                        self.assertLessEqual(e["vme"]["x"], e["vme"]["n"])
                        self.assertEqual(e["me"]["rate"], e["me"]["x"] / e["me"]["n"])
                        self.assertLess(e["me"]["rate"], 0.03, e["dataset"])

    def test_guidance_reads_the_rows_uncertainty_reads(self):
        for release in RELEASES:
            g, u = load(release, "guidance.json"), load(release, "uncertainty.json")
            self.assertEqual(list(g["rows"]), list(u["error_bars"]["rows"]))
            for key, r in g["rows"].items():
                b = u["error_bars"]["rows"][key]
                with self.subTest(release=release, row=key):
                    self.assertEqual(r["n"], b["n"])
                    self.assertEqual(r["screen"], {k: b[k] for k in ("shown_to_meet_both", "shown_to_miss_one", "not_decided")})

    def test_guidance_counts_the_no_call_rows_index_counts_finds(self):
        for release, revision in RELEASES.items():
            g, c = load(release, "guidance.json"), load(release, "index-counts.json")
            with self.subTest(release=release):
                # a learner row counts fewer than every genome exactly when it leaves out no-call genomes
                fewer = sum(g["rows"][f"{l}:{s}"]["fewer_than_every_genome"] for l in LEARNERS for s in ("random", "lineage"))
                self.assertEqual(fewer, c["unflagged_with_no_call"])
                self.assertEqual(g["rows"]["catalog:random"]["fewer_than_every_genome"], 0)
                self.assertEqual(g["rows"]["catalog:lineage"]["fewer_than_every_genome"], 0)
                if revision == 14:
                    self.assertEqual(sum(r["fewer_than_every_genome"] for r in g["rows"].values()), 0)

    def test_readme_lists_the_sha256_of_each_output(self):
        readme = (REFERENCE / "README.md").read_text(encoding="utf-8")
        listed = dict((path, digest) for digest, path in re.findall(r"^([0-9a-f]{64})  (\S+)$", readme, re.M))
        expected = sorted(f"{release}/{name}" for release in RELEASES for name in OUTPUTS)
        self.assertEqual(sorted(listed), expected)
        for path, digest in listed.items():
            self.assertEqual(hashlib.sha256((REFERENCE / path).read_bytes()).hexdigest(), digest, path)


if __name__ == "__main__":
    unittest.main()
