"""docs/release-summary.md and the reference outputs it quotes.

Each table of the document is rebuilt here from the files under reference/release-summary, and the document must hold
it line for line. The sentences that carry a count are checked one by one. The section "Curation path" quotes three
committed pages by line; each quote must stand on the line the section names."""

import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ROOT / "reference" / "release-summary"
DOCUMENT = ROOT / "docs" / "release-summary.md"
A, B = "rev12-274-20261105", "rev14-274-20261112"
NAMES = {"scm": "Set covering machine", "cart": "Decision tree", "rf": "Random forest", "gbt": "Boosted trees",
         "lasso": "Lasso", "best": "Best learner", "catalog": "Known-gene rule"}
CAP = 2000  # the cap's default, src/main.rs:39-40
MIN_CLASS = 25  # the rule of 25 resistant and 25 susceptible genomes, src/curate.rs:811 and README.md:54
SPECIES_CHECK = "species check: sequence is another species"  # src/curate.rs:543
# The counts of committed pages that the section "Curation path" sets beside the index, and the lines they stand on.
PAIRS_PAGE, POOL_PAGE, CAP_LEFT_OUT, CHECK_FAILED, CHECK_JUDGED = 1064823, 80967, 57206, 96, 35231


def load(release, name):
    return json.loads((REFERENCE / release / name).read_text(encoding="utf-8"))


def cell(r):
    lo, hi = r["ci95_species"]
    return f"{r['value']:.4f} ({lo:.4f} to {hi:.4f})"


def n(x):
    return f"{x:,}"


def path_rows():
    """The table of the section "Curation path", from revision 14's path.json."""
    p = load(B, "path.json")
    f, k, s = p["funnel"], p["funnel"]["pair_kinds"], p["species"]
    return [
        "| Step | Count | Left out |",
        "| --- | --- | --- |",
        f"| Lab records | {n(f['records'])} | |",
        f"| Genomes those records cover | {n(f['genomes_with_a_record'])} | |",
        f"| Genomes that pass quality checks | {n(f['genomes_passing'])} | {n(f['genomes_failing_a_check']['total'])} "
        f"genomes fail a check and {n(f['genomes_without_metadata'])} have no metadata |",
        f"| Genome-and-drug pairs on those genomes | {n(f['pairs'])} | {n(f['records_dropped_by_name']['total'])} records "
        f"on any genome name a drug class, no antibiotic or an unresolved name; {n(f['records_no_call_no_mic'])} records "
        f"on passing genomes have no call and no MIC |",
        f"| Pairs with a reported resistant or susceptible call | {n(k['reported_resistant_or_susceptible'])} | "
        f"{n(k['label_from_mic'])} pairs have a label worked out from a MIC, {n(k['no_label'])} have no label, "
        f"{n(k['intermediate_only'])} are intermediate only and {n(k['conflict'])} conflict |",
        f"| Datasets curated | {n(s['datasets_curated'])} | species and drug combinations with under {MIN_CLASS} resistant "
        f"or under {MIN_CLASS} susceptible genomes; the index does not count them |",
        f"| Genomes analyzed | {n(s['genomes_analyzed'])} | genomes past the cap of {n(CAP)} a species; the index does not "
        f"count them |",
        f"| Datasets with a report | {n(s['datasets_with_report'])} | {n(s['datasets_curated'] - s['datasets_with_report'])} "
        f"of {n(s['datasets_curated'])} datasets have one class under {MIN_CLASS} among the genomes kept |",
    ]


def tables():
    """{name: the table's lines}, each built from the reference outputs of both releases."""
    ua, ub = load(A, "uncertainty.json"), load(B, "uncertainty.json")
    ca, cb = load(A, "catalog.json"), load(B, "catalog.json")
    out = {}
    rows = ["| Learner | Split | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- | --- |"]
    for l in ["scm", "cart", "rf", "gbt", "lasso", "best"]:
        for s in ("random", "lineage"):
            a, b = ua["medians"][f"{l}:{s}"], ub["medians"][f"{l}:{s}"]
            rows.append(f"| {NAMES[l]} | {s} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["medians"] = rows
    rows = ["| Learner | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- |"]
    for l in ["scm", "cart", "rf", "gbt", "lasso", "best"]:
        a, b = ua["gaps"][f"{l}:gap"], ub["gaps"][f"{l}:gap"]
        rows.append(f"| {NAMES[l]} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["gaps"] = rows
    rows = ["| Set covering machine minus | Datasets | Revision 12 | Revision 14 | Set covering machine at least as high |",
            "| --- | --- | --- | --- | --- |"]
    for l in ["cart", "rf", "gbt", "lasso"]:
        a, b = ua["paired_lineage"][f"scm-{l}:lineage"], ub["paired_lineage"][f"scm-{l}:lineage"]
        rows.append(f"| {NAMES[l]} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} | "
                    f"{a['first_at_least_second']} / {b['first_at_least_second']} |")
    out["paired"] = rows
    rows = ["| Statistic | Datasets | Revision 12 | Revision 14 |", "| --- | --- | --- | --- |"]
    for k, label in (("median_lineage_of_best", "Lineage score of the best learner"),
                     ("median_lineage_of_oracle", "Highest lineage score of the five"),
                     ("median_oracle_minus_best", "Highest of the five minus the best learner")):
        a, b = ua["selection"][k], ub["selection"][k]
        rows.append(f"| {label} | {a['n']} / {b['n']} | {cell(a)} | {cell(b)} |")
    out["selection"] = rows
    rows = ["| Learner | Split | Datasets | Point within both | Shown to meet both | Shown to miss one | Not decided |",
            "| --- | --- | --- | --- | --- | --- | --- |"]
    for who in ["best", "scm", "cart", "rf", "gbt", "lasso", "catalog"]:
        for s in ("random", "lineage"):
            a, b = ua["error_bars"]["rows"][f"{who}:{s}"], ub["error_bars"]["rows"][f"{who}:{s}"]
            f = [f"{a[k]} / {b[k]}" for k in ("n", "point_meets_both", "shown_to_meet_both", "shown_to_miss_one", "not_decided")]
            rows.append(f"| {NAMES[who]} | {s} | " + " | ".join(f) + " |")
    out["error_limits"] = rows
    rows = ["| Rule's balanced accuracy | Datasets | Best learner ahead, random split | Best learner ahead, lineage split "
            "| Rule ahead or within 2 points, lineage split | Median lineage difference |",
            "| --- | --- | --- | --- | --- | --- |"]
    for k in ca["bands"]:
        a, b = ca["bands"][k], cb["bands"][k]
        f = lambda x: f"{a[x]} / {b[x]}"
        rows.append(f"| {k} | {f('datasets')} | {f('learner_ahead_random')} | {f('learner_ahead_lineage')} of "
                    f"{f('with_a_lineage_score')} | {f('rule_ahead_or_within_2pt_lineage')} | "
                    f"{a['median_lineage_difference']} / {b['median_lineage_difference']} |")
    out["bands"] = rows
    ga, gb = load(A, "guidance.json"), load(B, "guidance.json")
    rows = ["| Split | Learner | Datasets | Meet the guidance's wording | Shown to meet the release's screen |",
            "| --- | --- | --- | --- | --- |"]
    for s in ("random", "lineage"):
        for who in ["best", "scm", "cart", "rf", "gbt", "lasso", "catalog"]:
            a, b = ga["rows"][f"{who}:{s}"], gb["rows"][f"{who}:{s}"]
            rows.append(f"| {s} | {NAMES[who]} | {a['n']} / {b['n']} | {a['meets_wording']} / {b['meets_wording']} | "
                        f"{a['screen']['shown_to_meet_both']} / {b['screen']['shown_to_meet_both']} |")
    out["guidance"] = rows
    out["path"] = path_rows()
    return out


class Document(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.text = DOCUMENT.read_text(encoding="utf-8")

    def test_each_table_is_the_reference_files(self):
        for name, rows in tables().items():
            with self.subTest(table=name):
                self.assertIn("\n".join(rows) + "\n", self.text)

    def test_counted_sentences(self):
        ea, eb = load(A, "uncertainty.json")["error_bars"], load(B, "uncertainty.json")["error_bars"]
        self.assertEqual(ea["resistant_needed_for_zero_misses_to_show_vme"], eb["resistant_needed_for_zero_misses_to_show_vme"])
        self.assertEqual(ea["datasets_with_that_many_resistant"], eb["datasets_with_that_many_resistant"])
        ra, rb = load(A, "catalog.json")["rule"], load(B, "catalog.json")["rule"]
        self.assertEqual(ra, rb)
        self.assertEqual(ra["calls_nothing_resistant"] + ra["calls_everything_resistant"], ra["balanced_accuracy_exactly_half"])
        ia, ib = load(A, "index-counts.json"), load(B, "index-counts.json")
        sa = load(A, "uncertainty.json")
        sentences = [
            f"Each release has {sa['datasets']} datasets in {sa['species']} species.",
            f"a dataset needs at least {ea['resistant_needed_for_zero_misses_to_show_vme']} resistant genomes. "
            f"{ea['datasets_with_that_many_resistant']} datasets of each release have that many.",
            f"The rule has an entry on {ra['with_an_entry']} datasets. Its balanced accuracy is exactly 0.5 on "
            f"{ra['balanced_accuracy_exactly_half']} of them: it calls nothing resistant on {ra['calls_nothing_resistant']} "
            f"and everything resistant on {ra['calls_everything_resistant']}.",
            f"In both releases the counts of all {ia['rates_compared']} error rates of the {ia['rows']} intervals rows "
            f"rebuild from the index with {ia['mismatches']} mismatches.",
            f"In revision 14, {ib['flagged_rows']} learner rows carry the flag, and each of them has a no-call count.",
        ]
        ma, mb = ia["labels_from_mic_excluded"], ib["labels_from_mic_excluded"]
        self.assertEqual(ma, mb)
        sentences.append(
            f"In both releases the index records `labels_from_mic_excluded` on {ma['datasets_with_field']} of its "
            f"{ia['datasets']} datasets. On {ma['datasets_above_zero']} of them GRM inferred at least one label from a MIC "
            f"and left it out of training and scoring, {ma['labels']} labels in all.")
        self.assertEqual((ia["mismatches"], ib["mismatches"]), (0, 0))
        self.assertEqual(ib["flagged_rows"], ib["flagged_with_no_call"])
        for s in sentences:
            with self.subTest(sentence=s[:60]):
                self.assertIn(s, " ".join(self.text.split()))

    def test_guidance_sentences(self):
        ga, gb = load(A, "guidance.json"), load(B, "guidance.json")
        ra, rb = ga["rows"], gb["rows"]
        learner_rows = [f"{l}:{s}" for l in ("best", "scm", "cart", "rf", "gbt", "lasso") for s in ("random", "lineage")]
        both = lambda f: f(ra) if f(ra) == f(rb) else self.fail(f"the releases differ: {f(ra)} and {f(rb)}")
        meets = both(lambda r: (r["best:random"]["meets_wording"], r["best:lineage"]["meets_wording"]))
        shown = both(lambda r: (r["best:random"]["screen"]["shown_to_meet_both"], r["best:lineage"]["screen"]["shown_to_meet_both"]))
        self.assertTrue(all(ra[k]["meets_wording"] == rb[k]["meets_wording"] for k in ra))
        self.assertEqual(ga["best_meets_wording"], gb["best_meets_wording"])
        # revision 12 leaves out no-call genomes on boosted trees and lasso rows only; revision 14 on none
        fewer = {k: r["fewer_than_every_genome"] for k, r in ra.items() if r["fewer_than_every_genome"]}
        self.assertEqual(sorted(fewer), ["gbt:lineage", "gbt:random", "lasso:lineage", "lasso:random"])
        self.assertEqual(sum(fewer.values()), load(A, "index-counts.json")["unflagged_with_no_call"])
        self.assertEqual(sum(r["fewer_than_every_genome"] for r in rb.values()), 0)
        flagged = load(B, "index-counts.json")["flagged_rows"]
        table_8 = both(lambda r: (r["best:random"]["vme_n_under_table_8"], r["best:lineage"]["vme_n_under_table_8"],
                                  r["catalog:random"]["vme_n_under_table_8"], r["catalog:lineage"]["vme_n_under_table_8"]))
        first = ga["criteria"]["table_8_first_resistant"]
        self.assertEqual(first, gb["criteria"]["table_8_first_resistant"])
        tb = lambda r, k: r[k]["n"] - r[k]["without_tuberculosis"]["n"]
        tb_learners = both(lambda r: {tb(r, k) for k in learner_rows})
        tb_rule = both(lambda r: {tb(r, "catalog:random"), tb(r, "catalog:lineage")})
        self.assertEqual((len(tb_learners), len(tb_rule)), (1, 1))
        for r in (ra, rb):
            self.assertTrue(all(x["meets_wording"] == x["without_tuberculosis"]["meets_wording"] for x in r.values()))
        sentences = [
            f"Its first row is {first} resistant organisms with 0 discrepancies and an interval of 0.00 to 7.40%.",
            f"By the guidance's wording, the best learner meets on {meets[0]} datasets of the random split and {meets[1]} of "
            f"the lineage split, in both releases. By the release's screen, it is shown to meet on {shown[0]} and {shown[1]}.",
            f"Of the rows read here, {sum(fewer.values())} count fewer genomes than their class: boosted trees on "
            f"{fewer['gbt:random']} datasets of the random split and {fewer['gbt:lineage']} of the lineage split, and the "
            f"lasso on {fewer['lasso:random']} and {fewer['lasso:lineage']}. None of them is the best learner's row:",
            f"In revision 14, the {flagged} learner rows that carry the flag are scored on every genome, and every row read "
            f"here counts every genome of its class.",
            f"The best learner counts fewer than {first} resistant genomes on {table_8[0]} datasets of the random split and "
            f"{table_8[1]} of the lineage split in both releases, and the known-gene rule on {table_8[2]} and {table_8[3]}.",
            f"There are {tb_learners.pop()} of them on every learner row of both splits, and {tb_rule.pop()} on each split of "
            f"the known-gene rule. None of them meets the wording on any row of either release,",
        ]
        text = " ".join(self.text.split())
        for s in sentences:
            with self.subTest(sentence=s[:60]):
                self.assertIn(s, text)

    def test_path_sentences(self):
        p = load(B, "path.json")
        f, k, s = p["funnel"], p["funnel"]["pair_kinds"], p["species"]
        d, failing = f["records_dropped_by_name"], f["genomes_failing_a_check"]
        # no species sits above the cap, and some sit at it
        self.assertEqual(s["genomes_analyzed_highest"], CAP)
        # the reasons as the index writes them, most genomes first
        reasons = sorted(failing["reasons"].items(), key=lambda r: (-r[1], r[0]))
        listed = ", ".join(f"`{why}` {n(c)}" for why, c in reasons[:-1]) + f" and `{reasons[-1][0]}` {n(reasons[-1][1])}"
        sentences = [
            f"every antimicrobial resistance record of the download, {n(f['records'])} (`src/curate.rs:576-579`).",
            f"each genome a record names, counted once, {n(f['genomes_with_a_record'])} (`src/curate.rs:598-604`).",
            f"a genome with no metadata leaves first, {n(f['genomes_without_metadata'])} of them (`src/curate.rs:606-609`).",
            f"The index writes the reasons of the {n(failing['total'])} genomes that fail a check as {listed}. "
            f"{n(f['genomes_passing'])} genomes pass (`src/curate.rs:636`).",
            f"{n(d['total'])} records leave there: {n(d['drug_class'])} name a drug class, {n(d['not_an_antibiotic'])} name "
            f"something that is no antibiotic and {n(d['name_not_resolved'])} carry a name that is not resolved.",
            f"Records with no call and no MIC are counted on passing genomes only, {n(f['records_no_call_no_mic'])} of them "
            f"(`src/curate.rs:673-677`). The other records of passing genomes join into one pair for each genome and drug, "
            f"{n(f['pairs'])} pairs (`src/curate.rs:679-680`).",
            f"{n(k['reported_resistant_or_susceptible'])} have a reported resistant or susceptible call, "
            f"{n(k['label_from_mic'])} have a label worked out from a MIC, {n(k['no_label'])} have no label, "
            f"{n(k['intermediate_only'])} are intermediate only and {n(k['conflict'])} conflict. The five add up to the "
            f"{n(f['pairs'])} pairs.",
            f"at least {MIN_CLASS} of its genomes have an unconflicted reported resistant call and {MIN_CLASS} a susceptible "
            f"one, {n(s['datasets_curated'])} datasets (`src/curate.rs:802-813`). This rule of {MIN_CLASS} and {MIN_CLASS} "
            f"runs before the cap.",
            f"the cap keeps at most {n(CAP)} genomes a species (`src/main.rs:39-40`).",
            f"The genomes analyzed are those of each species' genome matrix, {n(s['genomes_analyzed'])} in all "
            f"(`src/report.rs:783-787`, `src/report.rs:995`). The highest count of any species is "
            f"{n(s['genomes_analyzed_highest'])}, and {s['species_at_highest']} of the {p['species_rows']} species sit at it.",
            f"{n(s['datasets_with_report'])} of the {n(s['datasets_curated'])} datasets have a report, and the other "
            f"{n(s['datasets_curated'] - s['datasets_with_report'])} have one class under {MIN_CLASS} genomes (`README.md:54`).",
            f"The species rows also count {n(s['genomes_with_any_label'])} genomes with a label of any kind on any drug, "
            f"intermediate calls and labels worked out from a MIC included (`src/report.rs:984-995`). That count is no step "
            f"of the path.",
        ]
        text = " ".join(self.text.split())
        for sentence in sentences:
            with self.subTest(sentence=sentence[:60]):
                self.assertIn(sentence, text)

    def test_path_differences_cite_their_lines(self):
        p = load(B, "path.json")
        f, s = p["funnel"], p["species"]
        failed = f["genomes_failing_a_check"]["reasons"][SPECIES_CHECK]
        cited = [
            ("docs/assumption-plan.md", 40, f"{n(PAIRS_PAGE)} pairs"),
            ("docs/assumption-plan.md", 7, f"of {n(POOL_PAGE)} genomes with a scorable label"),
            ("docs/assumption-plan.md", 7, f"The release's cap of {n(CAP)} leaves out {n(CAP_LEFT_OUT)}."),
            ("docs/species-check.md", 100, f"it judges {n(CHECK_JUDGED)} genomes"),
            ("docs/species-check.md", 100, f"{CHECK_FAILED} fail"),
            ("docs/species-check.md", 65, f"the {n(POOL_PAGE - CAP_LEFT_OUT)} genomes of the frozen 23-species candidate"),
            ("docs/species-check.md", 114, f"{n(s['genomes_analyzed'])} in all, and the file judges every one"),
            ("docs/species-check.md", 114, "None fails."),
            ("README.md", 54, f"Of its {s['datasets_curated']} curated datasets, {s['datasets_with_report']} have a report and "
                              f"{s['datasets_curated'] - s['datasets_with_report']} have none, because one class is under "
                              f"the minimum of {MIN_CLASS} genomes."),
        ]
        for path, number, quote in cited:
            with self.subTest(line=f"{path}:{number}", quote=quote[:40]):
                self.assertIn(quote, (ROOT / path).read_text(encoding="utf-8").split("\n")[number - 1])
        # each difference is still there; the section asserts neither side
        self.assertNotEqual(PAIRS_PAGE, f["pairs"])
        self.assertNotEqual(POOL_PAGE - CAP_LEFT_OUT, s["genomes_analyzed"])
        self.assertNotEqual(CHECK_FAILED, failed)
        sentences = [
            f"A committed page counts {n(POOL_PAGE)} genomes with a scorable label and says the release's cap of {n(CAP)} "
            f"leaves out {n(CAP_LEFT_OUT)} (`docs/assumption-plan.md:7`). The index has no key for either count,",
            f"`docs/assumption-plan.md:40` counts {n(PAIRS_PAGE)} pairs; the index counts {n(f['pairs'])}.",
            f"`docs/assumption-plan.md:7` gives {n(POOL_PAGE)} genomes less the {n(CAP_LEFT_OUT)} the cap leaves out, which "
            f"is {n(POOL_PAGE - CAP_LEFT_OUT)}; the index counts {n(s['genomes_analyzed'])} genomes analyzed.",
            f"`docs/species-check.md:100` counts {CHECK_FAILED} failures among {n(CHECK_JUDGED)} genomes the species check "
            f"judged; the index counts {failed} genomes that fail the species check.",
            f"One lead, and no more: `docs/species-check.md:65` gives {n(POOL_PAGE - CAP_LEFT_OUT)} as the genome count of "
            f"the frozen 23-species candidate of revision 9, the same figure as {n(POOL_PAGE)} less {n(CAP_LEFT_OUT)}.",
            f"The check judged every genome in the reports, {n(s['genomes_analyzed'])}, and none of them fails "
            f"(`docs/species-check.md:114`).",
        ]
        text = " ".join(self.text.split())
        for sentence in sentences:
            with self.subTest(sentence=sentence[:60]):
                self.assertIn(sentence, text)


if __name__ == "__main__":
    unittest.main()
