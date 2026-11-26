"""The release summary scripts (scripts/release-summary) on a small synthetic release.

The release is built here, in a temporary folder: 3 species (alpha with 5 datasets, beta with 5, gamma with 2), 12
datasets, both splits, the five learners, a known-gene rule on 8 datasets, one learner row with genomes in null-fit
folds, and labels inferred from a MIC on 4 datasets. Its index also holds a funnel block and 4 species rows, one of
them a species whose one curated dataset has no report. Each expected value below is worked out by hand from the
tables in this file; the comment beside it shows how.
"""

import json
import pathlib
import subprocess
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts" / "release-summary"
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]

SPECIES = ["alpha"] * 5 + ["beta"] * 5 + ["gamma"] * 2
IDS = [f"{s}__drug{i}" for i, s in enumerate(SPECIES)]
# Balanced accuracy of the set covering machine on each dataset, random and lineage split.
SCM_RANDOM = [0.90, 0.80, 0.70, 0.60, 0.95, 0.85, 0.75, 0.65, 0.92, 0.82, 0.72, 0.62]
SCM_LINEAGE = [0.80, 0.80, 0.66, 0.58, 0.75, 0.75, 0.75, 0.58, 0.86, 0.78, 0.71, 0.59]
# random minus lineage: 0.10 0.00 0.04 0.02 0.20 0.10 0.00 0.07 0.06 0.04 0.01 0.03
NO_LASSO = 9  # beta__drug9 has no lasso entry
N_RESISTANT = [200] * 11 + [245]  # gamma__drug11 has 245 resistant genomes
N_SUSCEPTIBLE = 200
# labels_from_mic_excluded of each dataset's index entry: above zero on drug1, drug4, drug9 and drug11, 3 + 7 + 1 + 12 = 23
MIC_EXCLUDED = [0, 3, 0, 0, 7, 0, 0, 0, 0, 1, 0, 12]
ABSENT = object()  # in a mic list: the index entry has no labels_from_mic_excluded

# The index's funnel block, with the keys src/curate.rs:474-506 writes. 120 genomes with a record = 100 passing + 5
# without metadata + 15 failing (7 + 2 + 6); 600 pairs = 400 + 120 + 50 + 20 + 10. The last two keys are not read.
FUNNEL = {
    "records": 1000, "records_class_label": 30, "records_not_antibiotic": 20, "records_unresolved_name": 10,
    "records_no_measure": 40, "genomes_with_records": 120, "genomes_without_metadata": 5,
    "genomes_failed_qc": {"BV-BRC quality: Poor": 7, "no species": 2, "species check: sequence is another species": 6},
    "genomes_passing_qc": 100, "pairs": 600, "pairs_conflict": 10, "pairs_intermediate_only": 20,
    "pairs_reported": 400, "pairs_labeled_from_mic": 120, "pairs_unlabeled": 50, "pairs_intermediate_mixed": 0,
    "pairs_call_against_mic": 3,
}
# The index's species rows: slug, datasets_curated, datasets_reported, genomes_labeled, genomes_analyzed. One alpha
# dataset falls short on the genomes kept, and delta's one curated dataset has no report and no genome matrix.
SPECIES_ROWS = [("alpha", 6, 5, 900, 400), ("beta", 5, 5, 700, 400), ("gamma", 2, 2, 300, 250), ("delta", 1, 0, 40, 0)]

# The known-gene rule: balanced accuracy, sensitivity, specificity. The values sit on the band edges.
CATALOG = {
    0: (0.5, 0.0, 1.0),  # calls nothing resistant
    1: (0.5, 1.0, 0.0),  # calls everything resistant
    2: (0.49, 0.0, 0.98),  # misses every resistant genome, but is not a constant call
    3: (0.5499, 0.3, 0.7998),  # just under 0.55
    4: (0.55, 0.5, 0.6),  # exactly 0.55: band "0.55 to 0.70"
    5: (0.70, 0.6, 0.8),  # exactly 0.70: band "0.70 to 0.85"
    6: (0.85, 0.8, 0.9),  # exactly 0.85: band "0.85 and up"
    7: (0.6999, 0.6, 0.7998),  # just under 0.70
}

# Intervals rows of the set covering machine on the random split, by dataset: (x, lo, hi) for VME, then for ME.
# Each n is the class size, 200 here. The other rows use DEFAULT_ROW.
DEFAULT_ROW = ((2, 0.001, 0.04), (4, 0.005, 0.05))
SCM_RANDOM_ROWS = {
    0: ((0, 0.0, 0.015), (0, 0.0, 0.03)),  # both upper bounds exactly at the limits: shown to meet
    1: ((6, 0.0150001, 0.06), (4, 0.001, 0.05)),  # VME lower bound just over 1.5%: shown to miss
    2: ((2, 0.001, 0.04), (10, 0.0300001, 0.09)),  # ME lower bound just over 3%: shown to miss
    3: ((3, 0.003, 0.0150001), (0, 0.0, 0.018)),  # VME point exactly 1.5%, upper bound just over: not decided
    4: ((4, 0.015, 0.05), (7, 0.03, 0.07)),  # lower bounds exactly at the limits: not decided
}


def learner_values(i):
    """{learner: (random, lineage)} on dataset i."""
    r, g = SCM_RANDOM[i], SCM_LINEAGE[i]
    v = {
        "scm": (r, g),
        "cart": (r, round(g - 0.05, 4)),
        "rf": (round(r + 0.01, 4), round(g + 0.01, 4)),
        "gbt": (0.9 if SPECIES[i] == "gamma" else 0.5, 0.5),
        "lasso": (r, g),
    }
    if i == NO_LASSO:
        del v["lasso"]
    return v


def build(folder, revision, flag=None, scored_on=None, mic=MIC_EXCLUDED, funnel=FUNNEL):
    """Write index.json and intervals/ for the synthetic release into folder.

    mic gives each dataset's labels_from_mic_excluded, or ABSENT to leave the field out of its index entry. funnel is
    the index's funnel block, or None to leave the block out.

    One learner row has genomes in null-fit folds: lasso, lineage split, gamma__drug10, with 4 resistant and 6
    susceptible. flag says whether its index entry carries null_fits_scored: true; by default it does from revision
    14. scored_on chooses the denominators its intervals row uses: "all" for every genome, "called" for the genomes
    outside null-fit folds; by default "all" when the row is flagged."""
    folder = pathlib.Path(folder)
    (folder / "intervals").mkdir(parents=True)
    if flag is None:
        flag = revision >= 14
    if scored_on is None:
        scored_on = "all" if flag else "called"
    datasets = []
    for i, ident in enumerate(IDS):
        n_res, n_sus = N_RESISTANT[i], N_SUSCEPTIBLE
        values = learner_values(i)
        splits = {"random": [], "lineage": []}
        entries = []
        for learner, (ba_r, ba_l) in values.items():
            entry = {"learner": learner}
            for scheme, ba in (("random", ba_r), ("lineage", ba_l)):
                (vx, vlo, vhi), (mx, mlo, mhi) = DEFAULT_ROW
                if learner == "scm" and scheme == "random" and i in SCM_RANDOM_ROWS:
                    (vx, vlo, vhi), (mx, mlo, mhi) = SCM_RANDOM_ROWS[i]
                no_call = {"resistant": 0, "susceptible": 0}
                vn, mn = n_res, n_sus
                split_entry = {}
                if learner == "lasso" and scheme == "lineage" and i == 10:
                    no_call = {"resistant": 4, "susceptible": 6}
                    if flag:
                        split_entry["null_fits_scored"] = True
                    if scored_on == "called":
                        vn, mn = n_res - 4, n_sus - 6
                if learner == "cart" and scheme == "lineage" and i == 11:
                    # a row with no-calls and no flag keeps the subtraction in every revision
                    no_call = {"resistant": 2, "susceptible": 0}
                    vn = n_res - 2
                split_entry = {"balanced_accuracy": ba, "vme": vx / vn, "me": mx / mn, "no_call": no_call, **split_entry}
                entry[scheme] = split_entry
                splits[scheme].append({"learner": learner, "vme": {"x": vx, "n": vn, "lo": vlo, "hi": vhi},
                                       "me": {"x": mx, "n": mn, "lo": mlo, "hi": mhi}})
            entries.append(entry)
        catalog = None
        if i in CATALOG:
            ba, sens, spec = CATALOG[i]
            vx, mx = round((1 - sens) * n_res), round((1 - spec) * n_sus)
            catalog = {"balanced_accuracy": ba, "sensitivity": sens, "specificity": spec, "vme": vx / n_res,
                       "me": mx / n_sus}
            for scheme in splits:
                splits[scheme].append({"learner": "catalog", "vme": {"x": vx, "n": n_res, "lo": 0.0, "hi": 1.0},
                                       "me": {"x": mx, "n": n_sus, "lo": 0.0, "hi": 1.0}})
        best = max(values, key=lambda l: values[l][0])
        datasets.append({
            "id": ident, "species": SPECIES[i].title(), "species_slug": SPECIES[i], "drug_class": f"class{i % 2}",
            "n": n_res + n_sus, "n_resistant": n_res, "learners": entries, "best_learner": best,
            "best_random": values[best][0], "best_lineage": values[best][1], "catalog": catalog,
            "mic_essential_agreement": None,
        })
        if mic[i] is not ABSENT:
            datasets[-1]["labels_from_mic_excluded"] = mic[i]
        intervals = {"dataset": ident, "evaluation_revision": revision, "method": "clopper-pearson",
                     "splits": [{"scheme": s, "learners": rows} for s, rows in splits.items()]}
        (folder / "intervals" / f"{ident}.json").write_text(json.dumps(intervals))
    species = [{"species": slug.title(), "slug": slug, "datasets_curated": curated, "datasets_reported": reported,
                "genomes_labeled": labeled, "genomes_analyzed": analyzed, "lineages": None, "lineage_threshold": None}
               for slug, curated, reported, labeled, analyzed in SPECIES_ROWS]
    index = {"evaluation_revision": revision, "generated_at": "2026-01-01T00:00:00Z", "species": species,
             "datasets": datasets}
    if funnel is not None:
        index["funnel"] = funnel
    (folder / "index.json").write_text(json.dumps(index))
    return folder


def run(script, *args, check=True):
    p = subprocess.run([sys.executable, "-I", str(SCRIPTS / script), *map(str, args)], capture_output=True, text=True)
    if check and p.returncode != 0:
        raise AssertionError(f"{script} exited {p.returncode}: {p.stderr}")
    return p


def run_json(script, folder, *extra):
    folder = pathlib.Path(folder)
    return json.loads(run(script, folder / "index.json", folder / "intervals", *extra).stdout)


class Synthetic(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        base = pathlib.Path(cls.tmp.name)
        cls.rev12 = build(base / "rev12", 12)
        cls.rev14 = build(base / "rev14", 14)
        cls.unc = run_json("uncertainty.py", cls.rev14)
        cls.summary = run_json("summary.py", cls.rev14)
        cls.catalog = run_json("catalog.py", cls.rev14)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_the_release_passes_its_own_checks(self):
        c = self.unc["checks"]
        self.assertEqual(c["problem_count"], 0, c["problems"])
        self.assertEqual(c["intervals_files"], 12)
        # 12 datasets x 2 splits x 5 learners, less lasso on one dataset (2), plus the rule on 8 datasets x 2 splits
        self.assertEqual(c["rows"], 120 - 2 + 16)
        self.assertEqual((self.unc["datasets"], self.unc["species"], self.unc["revision"]), (12, 3, 14))

    def test_medians_per_learner_and_split(self):
        m = self.unc["medians"]
        # SCM_RANDOM sorted: .60 .62 .65 .70 .72 .75 | .80 .82 .85 .90 .92 .95; median (.75 + .80) / 2
        self.assertEqual((m["scm:random"]["value"], m["scm:random"]["n"]), (0.775, 12))
        # SCM_LINEAGE sorted: .58 .58 .59 .66 .71 .75 | .75 .75 .78 .80 .80 .86; median .75
        self.assertEqual(m["scm:lineage"]["value"], 0.75)
        self.assertEqual(m["cart:random"]["value"], 0.775)
        self.assertEqual(m["cart:lineage"]["value"], 0.70)  # every lineage value 0.05 lower
        self.assertEqual(m["rf:random"]["value"], 0.785)
        self.assertEqual(m["rf:lineage"]["value"], 0.76)
        self.assertEqual(m["gbt:random"]["value"], 0.5)  # 10 of 12 at 0.5
        # lasso leaves out beta__drug9 (.82 random, .78 lineage): 11 values, the 6th of the sorted list
        self.assertEqual((m["lasso:random"]["value"], m["lasso:random"]["n"]), (0.75, 11))
        self.assertEqual((m["lasso:lineage"]["value"], m["lasso:lineage"]["n"]), (0.75, 11))
        s = self.summary["per_learner"]["scm"]
        self.assertEqual((s["n_both"], s["median_random"], s["median_lineage"]), (12, 0.775, 0.75))
        self.assertEqual(self.summary["best_learner_counts"], {"scm": 0, "cart": 0, "rf": 10, "gbt": 2, "lasso": 0})

    def test_paired_drop_from_random_to_lineage(self):
        # gaps sorted: 0 0 .01 .02 .03 .04 | .04 .06 .07 .10 .10 .20; median .04
        self.assertEqual(self.unc["gaps"]["scm:gap"]["value"], 0.04)
        self.assertEqual(self.unc["gaps"]["cart:gap"]["value"], 0.09)
        s = self.summary["per_learner"]["scm"]
        self.assertEqual(s["median_gap"], 0.04)
        self.assertEqual(s["share_gap_over_5pt"], 0.417)  # .06 .07 .10 .10 .20: 5 of 12
        # best learner: rf on alpha and beta (the scm gaps), gbt on gamma (0.9 - 0.5 twice)
        # sorted: 0 0 .02 .04 .04 .06 | .07 .10 .10 .20 .40 .40; median (.06 + .07) / 2
        b = self.summary["best"]
        self.assertEqual((b["n_both"], b["median_gap"], b["drop_over_5pt"], b["lineage_at_least_random"]), (12, 0.065, 7, 2))

    def test_paired_difference_between_learners(self):
        p = self.unc["paired_lineage"]
        self.assertEqual((p["scm-cart:lineage"]["value"], p["scm-cart:lineage"]["first_at_least_second"]), (0.05, 12))
        self.assertEqual((p["scm-rf:lineage"]["value"], p["scm-rf:lineage"]["first_at_least_second"]), (-0.01, 0))
        self.assertEqual(p["scm-gbt:lineage"]["value"], 0.25)  # SCM_LINEAGE - 0.5, median .75 - .5
        lasso = p["scm-lasso:lineage"]
        self.assertEqual((lasso["value"], lasso["n"], lasso["first_at_least_second"]), (0.0, 11, 11))

    def test_species_bootstrap_resamples_whole_species(self):
        r = self.unc["medians"]["gbt:random"]
        # gbt is 0.9 on gamma's 2 datasets and 0.5 on the other 10. Resampling 3 species: only a draw of gamma three
        # times (1/27, 3.7% of draws) gives a median of 0.9, so the upper 2.5% end is 0.9 and the lower end 0.5.
        self.assertEqual(r["ci95_species"], [0.5, 0.9])
        # Resampling 12 datasets: a median above 0.5 needs 6 or more gamma draws at p = 1/6 (under 1% of draws).
        self.assertEqual(r["ci95_datasets"], [0.5, 0.5])

    def test_same_seed_same_interval(self):
        self.assertEqual((self.unc["seed"], self.unc["reps"]), (20261111, 4000))
        again = run_json("uncertainty.py", self.rev14)
        self.assertEqual(again, self.unc)
        short = run_json("uncertainty.py", self.rev14, "--reps", "200")
        self.assertEqual((short["seed"], short["reps"]), (20261111, 200))
        self.assertEqual(short, run_json("uncertainty.py", self.rev14, "--reps", "200"))

    def test_error_limit_classes_at_their_edges(self):
        row = self.unc["error_bars"]["rows"]["scm:random"]
        self.assertEqual(row["n"], 12)
        self.assertEqual(row["shown_to_meet_both"], 1)  # drug0: upper bounds exactly 1.5% and 3%
        self.assertEqual(row["shown_to_miss_one"], 2)  # drug1 (VME lower bound just over), drug2 (ME)
        self.assertEqual((row["shown_to_miss_vme"], row["shown_to_miss_me"]), (1, 1))
        self.assertEqual(row["not_decided"], 9)  # drug3, drug4 and the 7 default rows
        self.assertEqual(row["shown_to_meet_both"] + row["shown_to_miss_one"] + row["not_decided"], row["n"])
        # points: drug1 VME 6/200, drug4 VME 4/200 and ME 7/200, drug2 ME 10/200 miss; drug3 VME 3/200 = 1.5% meets
        self.assertEqual((row["point_meets_vme"], row["point_meets_me"], row["point_meets_both"]), (10, 10, 9))
        self.assertEqual(row["median_resistant_scored"], 200)
        # 1 - 0.025 ** (1 / n) <= 0.015 first holds at n = 245; gamma__drug11 alone has 245 resistant genomes
        self.assertEqual(self.unc["error_bars"]["resistant_needed_for_zero_misses_to_show_vme"], 245)
        self.assertEqual(self.unc["error_bars"]["datasets_with_that_many_resistant"], 1)

    def test_known_gene_bands_at_their_edges(self):
        c = self.catalog
        self.assertEqual(c["rule"], {"with_an_entry": 8, "calls_nothing_resistant": 1, "calls_everything_resistant": 1,
                                     "balanced_accuracy_exactly_half": 2})
        bands = {k: v["datasets"] for k, v in c["bands"].items()}
        # under 0.55: 0.5, 0.5, 0.49, 0.5499; 0.55 to 0.70: 0.55, 0.6999; 0.70 to 0.85: 0.70; 0.85 and up: 0.85
        self.assertEqual(bands, {"under 0.55": 4, "0.55 to 0.70": 2, "0.70 to 0.85": 1, "0.85 and up": 1,
                                 "0.55 and up": 4, "all": 8})
        self.assertEqual(c["at_chance_by"]["failure"], {"misses_resistance": 2, "calls_all_resistant": 1, "neither": 1})
        self.assertEqual(self.summary["catalog"]["with_rule"], 8)

    def test_denominators_revision_12(self):
        out = run_json("index_counts.py", self.rev12)
        self.assertEqual((out["revision"], out["rule"]), (12, "no_call genomes left out"))
        self.assertEqual((out["rows"], out["learner_rows"], out["rule_rows"]), (134, 118, 16))
        self.assertEqual((out["rates_compared"], out["mismatches"]), (268, 0), out["examples"])
        self.assertEqual((out["flagged_rows"], out["unflagged_with_no_call"]), (0, 2))

    def test_denominators_revision_14(self):
        out = run_json("index_counts.py", self.rev14)
        self.assertEqual((out["revision"], out["rule"]), (14, "flagged rows scored on every genome"))
        self.assertEqual((out["rates_compared"], out["mismatches"]), (268, 0), out["examples"])
        self.assertEqual((out["flagged_rows"], out["flagged_with_no_call"], out["unflagged_with_no_call"]), (1, 1, 1))
        self.assertEqual(out["flag_keys_seen"], {"split:true": 1})

    def test_each_rule_finds_the_other_rules_denominators_wrong(self):
        with tempfile.TemporaryDirectory() as tmp:
            # revision 14, flagged row, but intervals that leave out the 4 + 6 no-call genomes (196 and 194)
            wrong14 = run_json("index_counts.py", build(pathlib.Path(tmp) / "a", 14, scored_on="called"))
            self.assertEqual(wrong14["mismatches"], 2)
            self.assertIn("gamma__drug10 lineage lasso vme: index 2/200, file 2/196", wrong14["examples"])
            self.assertIn("gamma__drug10 lineage lasso me: index 4/200, file 4/194", wrong14["examples"])
            # revision 12 index with the flag and intervals over every genome: revision 12 ignores the flag
            wrong12 = run_json("index_counts.py", build(pathlib.Path(tmp) / "b", 12, flag=True))
            self.assertEqual((wrong12["mismatches"], wrong12["flagged_rows"]), (2, 0))
            self.assertEqual(wrong12["flag_keys_seen"], {"split:true": 1})
            self.assertIn("gamma__drug10 lineage lasso vme: index 2/196, file 2/200", wrong12["examples"])
            self.assertIn("gamma__drug10 lineage lasso me: index 4/194, file 4/200", wrong12["examples"])

    def test_labels_inferred_from_a_mic_are_counted(self):
        for folder in (self.rev12, self.rev14):
            out = run_json("index_counts.py", folder)
            with self.subTest(revision=out["revision"]):
                # every entry holds the field; drug1, drug4, drug9 and drug11 are above zero: 3 + 7 + 1 + 12
                self.assertEqual(out["labels_from_mic_excluded"], {"datasets_with_field": 12, "datasets_without_field": 0,
                                                                   "datasets_above_zero": 4, "labels": 23})
                self.assertEqual(list(out)[-1], "labels_from_mic_excluded")
        with tempfile.TemporaryDirectory() as tmp:
            zero = run_json("index_counts.py", build(pathlib.Path(tmp) / "zero", 14, mic=[0] * 12))
        # zero on every entry is a count of zero, and it is not the same as no field
        self.assertEqual(zero["labels_from_mic_excluded"], {"datasets_with_field": 12, "datasets_without_field": 0,
                                                            "datasets_above_zero": 0, "labels": 0})

    def test_a_missing_field_is_never_counted_as_zero(self):
        with tempfile.TemporaryDirectory() as tmp:
            none = run_json("index_counts.py", build(pathlib.Path(tmp) / "none", 12, mic=[ABSENT] * 12))
            # drug1 and drug3 have no field: 10 hold it, and drug4, drug9 and drug11 are above zero, 7 + 1 + 12
            some = [ABSENT if i in (1, 3) else v for i, v in enumerate(MIC_EXCLUDED)]
            part = run_json("index_counts.py", build(pathlib.Path(tmp) / "part", 14, mic=some))
        self.assertEqual(none["labels_from_mic_excluded"], {"datasets_with_field": 0, "datasets_without_field": 12,
                                                            "datasets_above_zero": None, "labels": None})
        self.assertEqual(part["labels_from_mic_excluded"], {"datasets_with_field": 10, "datasets_without_field": 2,
                                                            "datasets_above_zero": 3, "labels": 20})
        # the field changes no other count of the file
        self.assertEqual({k: v for k, v in none.items() if k != "labels_from_mic_excluded"},
                         {k: v for k, v in run_json("index_counts.py", self.rev12).items() if k != "labels_from_mic_excluded"})

    def test_a_mic_value_that_is_not_a_count_is_refused(self):
        for value in (None, -1, True, 1.5, "3"):
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                folder = build(pathlib.Path(tmp) / "bad", 14, mic=[value if i == 4 else v for i, v in enumerate(MIC_EXCLUDED)])
                p = run("index_counts.py", folder / "index.json", folder / "intervals", check=False)
                self.assertEqual((p.returncode, p.stdout), (1, ""))
                self.assertIn(f"alpha__drug4 has labels_from_mic_excluded {value!r}, which is not a count", p.stderr)

    def test_unknown_revision_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = build(pathlib.Path(tmp) / "rev13", 13)
            for script in ("uncertainty.py", "summary.py", "catalog.py", "index_counts.py", "path.py"):
                p = run(script, folder / "index.json", folder / "intervals", check=False)
                self.assertEqual(p.returncode, 1, script)
                self.assertEqual(p.stdout, "", script)
                self.assertIn("evaluation revision 13", p.stderr, script)
                self.assertIn("understands revisions 12, 14", p.stderr, script)
            doc = pathlib.Path(tmp) / "u13.json"
            doc.write_text(json.dumps({**self.unc, "revision": 13}))
            p = run("compare.py", doc, doc, check=False)
            self.assertEqual((p.returncode, p.stdout), (1, ""))
            self.assertIn("understands revisions 12, 14", p.stderr)

    def test_compare_prints_both_revisions(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = pathlib.Path(tmp) / "a.json"
            a.write_text(json.dumps(run_json("uncertainty.py", self.rev12, "--reps", "50")))
            b = pathlib.Path(tmp) / "b.json"
            b.write_text(json.dumps(self.unc))
            lines = run("compare.py", a, b).stdout.splitlines()
        self.assertEqual(lines[0], "revisions 12 and 14; datasets 12 and 12; reps 50 and 4000")

    def path_on(self, change):
        """path.py on the revision 14 index after change(index), with the release's own intervals directory."""
        index = json.loads((self.rev14 / "index.json").read_text())
        change(index)
        with tempfile.TemporaryDirectory() as tmp:
            path = pathlib.Path(tmp) / "index.json"
            path.write_text(json.dumps(index))
            return run("path.py", path, self.rev14 / "intervals", check=False)

    def test_path_sums(self):
        expected = {
            "datasets": 12, "species_rows": 4,
            "funnel": {
                "records": 1000,
                "records_dropped_by_name": {"drug_class": 30, "not_an_antibiotic": 20, "name_not_resolved": 10,
                                            "total": 60},  # 30 + 20 + 10
                "records_no_call_no_mic": 40,
                "genomes_with_a_record": 120,
                "genomes_without_metadata": 5,
                "genomes_failing_a_check": {"reasons": FUNNEL["genomes_failed_qc"], "total": 15},  # 7 + 2 + 6
                "genomes_passing": 100,
                "pairs": 600,
                "pair_kinds": {"reported_resistant_or_susceptible": 400, "label_from_mic": 120, "no_label": 50,
                               "intermediate_only": 20, "conflict": 10},
            },
            # curated 6 + 5 + 2 + 1, reported 5 + 5 + 2 + 0, labeled 900 + 700 + 300 + 40, analyzed 400 + 400 + 250 + 0;
            # alpha and beta sit at the highest, 400
            "species": {"datasets_curated": 14, "datasets_with_report": 12, "genomes_with_any_label": 1940,
                        "genomes_analyzed": 1050, "genomes_analyzed_highest": 400, "species_at_highest": 2},
            # n: 11 datasets of 200 + 200 and gamma__drug11 of 245 + 200; n_resistant: 11 x 200 + 245
            "dataset_rows": {"n": 4845, "n_resistant": 2445},
            "checks": {"genomes": {"with_a_record": 120, "passing_plus_without_metadata_plus_failing": 120},
                       "pairs": {"pairs": 600, "sum_of_five_kinds": 600}},
        }
        for folder, revision in ((self.rev12, 12), (self.rev14, 14)):
            with self.subTest(revision=revision):
                out = run_json("path.py", folder)
                self.assertEqual(out, {"revision": revision, **expected})
                self.assertEqual(list(out), ["revision", "datasets", "species_rows", "funnel", "species", "dataset_rows",
                                             "checks"])

    def test_path_check_fails_on_a_wrong_count(self):
        genomes = "the genomes check fails: with_a_record is {} and passing_plus_without_metadata_plus_failing is {}"
        pairs = "the pairs check fails: pairs is {} and sum_of_five_kinds is {}"
        plants = [
            ({"genomes_passing_qc": 101}, genomes.format(120, 121)),
            ({"genomes_without_metadata": 4}, genomes.format(120, 119)),
            ({"genomes_with_records": 119}, genomes.format(119, 120)),
            ({"genomes_failed_qc": {**FUNNEL["genomes_failed_qc"], "no species": 3}}, genomes.format(120, 121)),
            ({"pairs_conflict": 11}, pairs.format(600, 601)),
            ({"pairs_labeled_from_mic": 119}, pairs.format(600, 599)),
            ({"pairs": 599}, pairs.format(599, 600)),
            # a pair the label rule sets aside for an intermediate call is none of the five kinds
            ({"pairs": 605, "pairs_intermediate_mixed": 5}, pairs.format(605, 600)),
        ]
        for plant, message in plants:
            with self.subTest(plant=plant):
                p = self.path_on(lambda index: index["funnel"].update(plant))
                self.assertEqual((p.returncode, p.stdout), (1, ""))
                self.assertIn(message, p.stderr)

    def test_path_without_funnel_is_null(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = run_json("path.py", build(pathlib.Path(tmp) / "nofunnel", 14, funnel=None))
        # no block is null, never zeros, and the checks need the block
        self.assertIsNone(out["funnel"])
        self.assertIsNone(out["checks"])
        self.assertEqual(list(out), ["revision", "datasets", "species_rows", "funnel", "species", "dataset_rows", "checks"])
        # the block changes no other count of the file
        full = run_json("path.py", self.rev14)
        self.assertEqual({k: v for k, v in out.items() if k not in ("funnel", "checks")},
                         {k: v for k, v in full.items() if k not in ("funnel", "checks")})

    def test_path_refuses_a_value_that_is_not_a_count(self):
        for value in (None, -1, True, 1.5, "3"):
            cases = [
                (lambda index: index["funnel"].update(records=value), f"funnel has records {value!r}"),
                (lambda index: index["funnel"]["genomes_failed_qc"].update({"no species": value}),
                 f"funnel genomes_failed_qc has no species {value!r}"),
                (lambda index: index["species"][2].update(genomes_analyzed=value),
                 f"species row gamma has genomes_analyzed {value!r}"),
                (lambda index: index["datasets"][4].update(n=value), f"dataset alpha__drug4 has n {value!r}"),
            ]
            for change, message in cases:
                with self.subTest(value=value, message=message[:30]):
                    p = self.path_on(change)
                    self.assertEqual((p.returncode, p.stdout), (1, ""))
                    self.assertIn(f"{message}, which is not a count", p.stderr)
        # a key the funnel block lacks is no count either, and the reasons must be a block of counts
        p = self.path_on(lambda index: index["funnel"].pop("pairs_unlabeled"))
        self.assertEqual((p.returncode, p.stdout), (1, ""))
        self.assertIn("funnel has pairs_unlabeled None, which is not a count", p.stderr)
        p = self.path_on(lambda index: index["funnel"].update(genomes_failed_qc=[7, 2, 6]))
        self.assertEqual((p.returncode, p.stdout), (1, ""))
        self.assertIn("funnel has genomes_failed_qc [7, 2, 6], which is not a block of counts", p.stderr)


if __name__ == "__main__":
    unittest.main()
