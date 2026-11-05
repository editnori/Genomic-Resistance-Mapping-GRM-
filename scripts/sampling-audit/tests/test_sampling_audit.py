"""The legacy-selection port and the sampling audit on hand-built fixtures.

Each expected value below is derived by hand from src/curate.rs
(`select_genomes`, `Phenotype::reported_label`) or from the definitions in
README.md, not from running the code under test.

usage: python3 -I scripts/sampling-audit/tests/test_sampling_audit.py
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
import audit  # noqa: E402
import legacy_select as ls  # noqa: E402

SP = "A b"


def ph(gid, drug, call, *, from_mic=False, conflict=False, species=SP):
    """A phenotype row as curate writes it: MIC-inferred rows have no
    reported call, conflicted rows have neither a reported call nor a label."""
    return {
        "genome_id": gid,
        "species": species,
        "antibiotic": drug,
        "reported": None if from_mic or conflict else call,
        "conflict": conflict,
        "label_from_mic": from_mic,
        "label": None if conflict else call,
    }


def ds(drug, species=SP):
    return {"id": f"{ls.slug(species)}__{ls.slug(drug)}", "species": species, "antibiotic": drug}


class LegacyMatch(unittest.TestCase):
    def test_rust_unit_test_rare_class_first_and_mic_rows_never_picked(self):
        # curate.rs selection_prefers_rare_reported_class_carriers, row for row.
        rows = [ph(f"s{i}", "x", "S") for i in range(9)]
        rows.append(ph("z-r", "x", "R"))
        rows += [ph(f"m{i}", "x", "R", from_mic=True) for i in range(10)]
        self.assertEqual(ls.select_genomes(rows, [ds("x")], SP, 2), ["z-r", "s0"])
        everything = ls.select_genomes(rows, [ds("x")], SP, 100)
        self.assertEqual(len(everything), 10)
        self.assertFalse(any(g.startswith("m") for g in everything))

    def test_many_panel_genome_outranks_rare_single_call(self):
        # Drug x: R={a}, S={b,c,d}; y: R={b,c}, S={d,e}; z: R={b}, S={c,d,e}.
        # a = 1/1 = 1; b = 1/3 + 1/2 + 1/1 = 1.8333; c = 1/3 + 1/2 + 1/3 = 1.1667;
        # d = 1/3 + 1/2 + 1/3 = 1.1667; e = 1/2 + 1/3 = 0.8333.
        rows = [ph("a", "x", "R"), ph("b", "x", "S"), ph("c", "x", "S"), ph("d", "x", "S"),
                ph("b", "y", "R"), ph("c", "y", "R"), ph("d", "y", "S"), ph("e", "y", "S"),
                ph("b", "z", "R"), ph("c", "z", "S"), ph("d", "z", "S"), ph("e", "z", "S")]
        sets = [ds("x"), ds("y"), ds("z")]
        # c and d tie exactly (same addends in the same order); c sorts first.
        self.assertEqual(ls.select_genomes(rows, sets, SP, 10), ["b", "c", "d", "a", "e"])
        self.assertEqual(ls.select_genomes(rows, sets, SP, 2), ["b", "c"])

    def test_only_dataset_drugs_and_the_named_species_count(self):
        # Drug w has no dataset, so b's w call neither scores nor grows a class.
        rows = [ph("a", "x", "R"), ph("b", "x", "S"), ph("c", "x", "S"),
                ph("b", "w", "R"), ph("q", "x", "R", species="C d")]
        # a = 1/1 (q is another species); b = c = 1/2.
        self.assertEqual(ls.select_genomes(rows, [ds("x"), ds("w", "C d")], SP, 10), ["a", "b", "c"])

    def test_intermediate_conflicted_and_mic_rows_do_not_score(self):
        rows = [ph("a", "x", "R"), ph("b", "x", "S"),
                ph("c", "x", "I"), ph("d", "x", None, conflict=True), ph("e", "x", "S", from_mic=True)]
        self.assertEqual(ls.select_genomes(rows, [ds("x")], SP, 10), ["a", "b"])
        # A reported call with label_from_mic set (not written by curate) is
        # still skipped, as reported_label() requires !label_from_mic.
        odd = dict(ph("f", "x", "R"), label_from_mic=True)
        self.assertEqual(ls.select_genomes(rows + [odd], [ds("x")], SP, 10), ["a", "b"])

    def test_mic_only_genome_is_not_picked_even_when_the_cap_is_free(self):
        # m has MIC-inferred labels on both drugs and no reported call.
        rows = [ph("a", "x", "R"), ph("b", "x", "S"), ph("a", "y", "S"), ph("b", "y", "R"),
                ph("m", "x", "R", from_mic=True), ph("m", "y", "R", from_mic=True)]
        self.assertEqual(ls.select_genomes(rows, [ds("x"), ds("y")], SP, 1000), ["a", "b"])

    def test_ties_break_on_genome_id_bytes(self):
        # Equal scores: ids sort by UTF-8 bytes, so "B" < "a" < "a.10" < "a.9" < "é".
        ids = ["é", "a.9", "a", "B", "a.10"]
        rows = [ph(g, "x", "S") for g in ids] + [ph("r", "x", "R")]
        self.assertEqual(ls.select_genomes(rows, [ds("x")], SP, 10), ["r", "B", "a", "a.10", "a.9", "é"])

    def test_scores_sum_left_to_right_as_rust_does(self):
        # g has ten calls in classes of ten: 1/10 added ten times in f64 is
        # 0.9999999999999999, below h's single 1/1. Python 3.12's sum()
        # compensates to exactly 1.0, which would tie and put g first by id.
        drugs = [f"d{i}" for i in range(10)]
        rows = []
        for d in drugs:
            rows.append(ph("g", d, "S"))
            rows += [ph(f"s{d}-{j}", d, "S") for j in range(9)]
            rows += [ph(f"r{d}-{j}", d, "R") for j in range(10)]
        rows += [ph("h", "solo", "R")] + [ph(f"t{j}", "solo", "S") for j in range(10)]
        sets = [ds(d) for d in drugs] + [ds("solo")]
        acc = 0.0
        for _ in range(10):
            acc += 0.1
        self.assertLess(acc, 1.0)
        if sys.version_info >= (3, 12):
            self.assertEqual(sum([0.1] * 10), 1.0)
        self.assertEqual(ls.select_genomes(rows, sets, SP, 2), ["h", "g"])

    def test_file_order_of_addends_is_kept(self):
        # Class sizes p=4, q=3, r=3. A's rows come p, q, r and B's r, q, p:
        # (1/4 + 1/3) + 1/3 < (1/3 + 1/3) + 1/4 in f64, so Rust ranks B
        # first; a port that sorted addends or summed exactly would tie them
        # and put A first by id.
        self.assertLess((0.25 + 1 / 3) + 1 / 3, (1 / 3 + 1 / 3) + 0.25)
        rows = [ph("A", "p", "S"), ph("A", "q", "S"), ph("A", "r", "S"),
                ph("B", "r", "S"), ph("B", "q", "S"), ph("B", "p", "S"),
                ph("p1", "p", "S"), ph("p2", "p", "S"), ph("q1", "q", "S"), ph("r1", "r", "S")]
        sets = [ds("p"), ds("q"), ds("r")]
        self.assertEqual(ls.select_genomes(rows, sets, SP, 2), ["B", "A"])

    def test_cli_prints_genomes_txt_order(self):
        rows = [ph("a", "x", "R"), ph("b", "x", "S"), ph("c", "x", "S")]
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "phenotypes.jsonl"), "w") as f:
                f.write("".join(json.dumps(r) + "\n" for r in rows))
            with open(os.path.join(d, "datasets.json"), "w") as f:
                json.dump([ds("x")], f)
            out = subprocess.run([sys.executable, "-I", os.path.join(LANE, "legacy_select.py"), d, SP, "2"],
                                 check=True, capture_output=True, text=True).stdout
        self.assertEqual(out, "a\nb\n")

    def test_slug_matches_curate(self):
        self.assertEqual(ls.slug("Mycobacterium tuberculosis"), "mycobacterium_tuberculosis")
        self.assertEqual(ls.slug("ampicillin/sulbactam"), "ampicillin_sulbactam")
        self.assertEqual(ls.slug("  A--b  "), "a_b")


def fixture():
    """Two drugs. x: R={a,b}, S={c,d,e,f}; y: R={c}, S={a,e,g}; MIC-inferred
    x row for m and y row for f, intermediate x row for g.
    Class sizes: (x,R)=2, (x,S)=4, (y,R)=1, (y,S)=3.
    Scores: c = 1/4 + 1 = 1.25; a = 1/2 + 1/3 = 0.8333; b = 0.5;
    e = 1/4 + 1/3 = 0.5833; g = 1/3; d = f = 0.25.
    Ranking: c, a, e, b, g, d, f (d before f by id); m is unscored.
    """
    rows = [
        ph("a", "x", "R"), ph("b", "x", "R"), ph("c", "x", "S"), ph("d", "x", "S"),
        ph("e", "x", "S"), ph("f", "x", "S"), ph("m", "x", "R", from_mic=True), ph("g", "x", "I"),
        ph("c", "y", "R"), ph("a", "y", "S"), ph("e", "y", "S"), ph("g", "y", "S"),
        ph("f", "y", "R", from_mic=True),
    ]
    return rows, [ds("x"), ds("y")]


GROUPS = {"a": "G1", "b": "G1", "c": "G2", "d": "G2", "e": "G3", "g": "G3"}  # f ungrouped


def added(groups, minority_groups):
    """A nested-cap group delta where both caps pass the coverage gate."""
    return {"groups_added": groups, "groups_added_null_reason": None,
            "minority_groups_added": minority_groups, "minority_groups_added_null_reason": None}


class Report(unittest.TestCase):
    def setUp(self):
        rows, sets = fixture()
        self.ranking = ls.ranking(rows, sets, SP)
        self.rep = audit.build_report(rows, sets, [2, 4, audit.ALL], {"dup": GROUPS})
        self.tasks = {t["antibiotic"]: t for t in self.rep["tasks"]}
        self.sp = self.rep["species"][0]

    def test_ranking_is_the_hand_ranking(self):
        self.assertEqual([g for _, g in self.ranking], ["c", "a", "e", "b", "g", "d", "f"])

    def test_eligible_counts_prevalence_and_minority(self):
        x = self.tasks["x"]
        self.assertEqual(x["minority"], "R")
        self.assertFalse(x["minority_tied"])
        e = x["eligible"]
        self.assertEqual((e["rows"], e["R"], e["S"], e["minority_rows"]), (6, 2, 4, 2))
        self.assertAlmostEqual(e["prevalence_R"], 2 / 6)
        self.assertEqual(e["non_target_rows"], {"mic_inferred": 1, "intermediate": 1, "conflict": 0, "unlabelled": 0})
        y = self.tasks["y"]
        self.assertEqual((y["eligible"]["rows"], y["eligible"]["R"], y["minority"]), (4, 1, "R"))

    def test_selected_rows_per_cap(self):
        # Cap 2 = {c, a}: x rows c(S), a(R); y rows c(R), a(S).
        # Cap 4 = {c, a, e, b}: x rows a, b, c, e; y rows c, a, e.
        x2, x4 = self.tasks["x"]["by_cap"]["2"], self.tasks["x"]["by_cap"]["4"]
        self.assertEqual((x2["rows"], x2["R"], x2["S"]), (2, 1, 1))
        self.assertEqual((x4["rows"], x4["R"], x4["S"]), (4, 2, 2))
        self.assertAlmostEqual(x2["minority_coverage"], 0.5)
        self.assertAlmostEqual(x4["minority_coverage"], 1.0)
        y4 = self.tasks["y"]["by_cap"]["4"]
        self.assertEqual((y4["rows"], y4["R"], y4["S"]), (3, 1, 2))
        # At "all", f is in: its MIC-inferred y row is in the cohort but not a task row.
        self.assertEqual(self.tasks["y"]["by_cap"]["all"]["non_target_rows"]["mic_inferred"], 1)
        self.assertEqual(self.tasks["y"]["by_cap"]["4"]["non_target_rows"]["mic_inferred"], 0)

    def test_panel_breadth(self):
        # Breadth (dataset drugs with a reported R/S): a=2, b=1, c=2, d=1, e=2, f=1, g=1.
        e = self.tasks["x"]["eligible"]["panel_breadth"]
        self.assertEqual(e["histogram"], {"1": 3, "2": 3})  # x rows a, b, c, d, e, f
        self.assertEqual((e["min"], e["median"], e["max"]), (1, 1.5, 2))
        s2 = self.sp["by_cap"]["2"]["panel_breadth"]
        self.assertEqual(s2["histogram"], {"2": 2})
        self.assertEqual(self.sp["panel_breadth_scored"]["histogram"], {"1": 4, "2": 3})
        self.assertEqual((self.sp["scored_genomes"], self.sp["unscored_genomes"]), (7, 1))

    def test_task_score_share(self):
        # Cap 2, task x: c's x share = (1/4)/1.25 = 0.2; a's = (1/2)/(5/6) = 0.6.
        self.assertAlmostEqual(self.tasks["x"]["by_cap"]["2"]["task_score_share_mean"], 0.4)

    def test_cutoff_and_id_tie_break(self):
        self.assertEqual(self.sp["by_cap"]["2"]["tie_broken_by_genome_id"], False)
        self.assertFalse(self.sp["by_cap"]["all"]["cap_binds"])
        rows, sets = fixture()
        rep = audit.build_report(rows, sets, [6])
        c = rep["species"][0]["by_cap"]["6"]
        # Cap 6 ends on d (0.25), tied with f, which is left out by id.
        self.assertEqual((c["cutoff_score"], c["tied_at_cutoff"], c["tied_at_cutoff_selected"]), (0.25, 2, 1))
        self.assertTrue(c["tie_broken_by_genome_id"])

    def test_nested_cap_deltas(self):
        # 2 -> 4 adds genomes e, b: x rows e(S), b(R); y row e(S).
        n = {(d["from_cap"], d["to_cap"]): d for d in self.tasks["x"]["nested"]}
        d = n[(2, 4)]
        self.assertEqual((d["genomes_added"], d["rows_added"], d["R_added"], d["S_added"], d["minority_rows_added"]),
                         (2, 2, 1, 1, 1))
        # 4 -> all adds g, d, f: x rows d(S), f(S) (g's x row is I).
        d = n[(4, audit.ALL)]
        self.assertEqual((d["genomes_added"], d["rows_added"], d["minority_rows_added"]), (3, 2, 0))
        y = {(d["from_cap"], d["to_cap"]): d for d in self.tasks["y"]["nested"]}
        self.assertEqual((y[(2, 4)]["rows_added"], y[(2, 4)]["minority_rows_added"]), (1, 0))

    def test_group_metrics(self):
        # Task x eligible rows a,b (G1) c,d (G2) e (G3) f (ungrouped): 5 of 6
        # grouped (0.83 >= 0.8), so figures are over sizes 2,2,1 and f is apart.
        g = self.tasks["x"]["eligible"]["groups"]["dup"]
        self.assertEqual((g["groups"], g["minority_groups"], g["largest_group_rows"], g["rows_without_group"],
                          g["rows_with_group"]), (3, 1, 2, 1, 5))
        self.assertAlmostEqual(g["group_coverage"], 5 / 6)
        self.assertAlmostEqual(g["largest_group_share"], 2 / 5)
        self.assertAlmostEqual(g["group_ess"], 25 / 9)
        self.assertIsNone(g["group_ess_null_reason"])
        self.assertAlmostEqual(g["minority_group_ess"], 1.0)  # both R rows (a, b) are in G1
        g2 = self.tasks["x"]["by_cap"]["2"]["groups"]["dup"]  # c (G2), a (G1)
        self.assertEqual((g2["groups"], g2["minority_groups"]), (2, 1))
        self.assertAlmostEqual(g2["group_ess"], 2.0)

    def test_group_deltas(self):
        # Every cap here passes the coverage gate (x at "all" has 5 of 6 rows grouped).
        n = {(d["from_cap"], d["to_cap"]): d for d in self.tasks["x"]["nested"]}
        # 2 -> 4 adds e (G3) and b (G1, already present, already minority).
        self.assertEqual(n[(2, 4)]["groups"]["dup"], added(1, 0))
        # 4 -> all adds d (G2, present) and f (ungrouped, so no group).
        self.assertEqual(n[(4, audit.ALL)]["groups"]["dup"], added(0, 0))
        # y: 2 -> 4 adds e (G3): a new group, not minority (y minority is R, only c).
        y = {(d["from_cap"], d["to_cap"]): d for d in self.tasks["y"]["nested"]}
        self.assertEqual(y[(2, 4)]["groups"]["dup"], added(1, 0))

    def test_minority_groups_added_counts_new_minority_groups(self):
        # x: R={a (G1), b (G2)}, S={c, d, e}; scores a = b = 1/2, c = d = e = 1/3.
        rows = [ph("a", "x", "R"), ph("b", "x", "R"), ph("c", "x", "S"), ph("d", "x", "S"), ph("e", "x", "S")]
        rep = audit.build_report(rows, [ds("x")], [1, 2], {"lin": {"a": "G1", "b": "G2", "c": "G1"}})
        d = rep["tasks"][0]["nested"][0]
        self.assertEqual((d["minority_rows_added"], d["groups"]["lin"]), (1, added(1, 1)))

    def test_groups_added_needs_both_caps_grouped(self):
        # The re-check's case: u (R) scores 1, v, w, z (S) score 1/3 each, so
        # cap 1 = {u} and cap 2 = {u, v}. If u is truly in v's group, nothing
        # is new; an ungrouped u would make G1 look new.
        rows = [ph("u", "x", "R"), ph("v", "x", "S"), ph("w", "x", "S"), ph("z", "x", "S")]

        def delta(groups):
            return audit.build_report(rows, [ds("x")], [1, 2], {"lin": groups})["tasks"][0]["nested"][0]["groups"]["lin"]

        self.assertEqual(delta({"u": "G1", "v": "G1"}), added(0, 0))
        small_ungrouped = delta({"v": "G1"})
        self.assertIsNone(small_ungrouped["groups_added"])
        self.assertEqual(small_ungrouped["groups_added_null_reason"],
                         "cap 1: 0 of 1 rows grouped, under 80%; cap 2: 1 of 2 rows grouped, under 80%")
        self.assertIsNone(small_ungrouped["minority_groups_added"])
        self.assertEqual(small_ungrouped["minority_groups_added_null_reason"],
                         "cap 1: 0 of 1 rows grouped, under 80%; cap 2: 0 of 1 rows grouped, under 80%")
        # Only the smaller cap failing is enough: u (R) alone at cap 1, then
        # five grouped S genomes make cap 6 pass at 5 of 6.
        wide = [ph("u", "x", "R")] + [ph(i, "x", "S") for i in "vwxyz"]
        d = audit.build_report(wide, [ds("x")], [1, 6], {"lin": {i: "G1" for i in "vwxyz"}})["tasks"][0]["nested"][0]
        self.assertEqual((d["groups"]["lin"]["groups_added"], d["groups"]["lin"]["groups_added_null_reason"]),
                         (None, "cap 1: 0 of 1 rows grouped, under 80%"))
        # An ungrouped genome at the larger cap could hide a new group.
        large_ungrouped = delta({"u": "G1"})
        self.assertIsNone(large_ungrouped["groups_added"])
        self.assertEqual(large_ungrouped["groups_added_null_reason"], "cap 2: 1 of 2 rows grouped, under 80%")
        self.assertEqual((large_ungrouped["minority_groups_added"], large_ungrouped["minority_groups_added_null_reason"]),
                         (0, None))

    def test_minority_gate_fires_when_overall_coverage_passes(self):
        # R={a, b}, S={c..j}, every genome but b grouped, each in its own group:
        # 9 of 10 rows grouped, but only 1 of 2 minority rows.
        ids = "abcdefghij"
        rows = [ph(i, "x", "R" if i in "ab" else "S") for i in ids]
        groups = {i: f"G{i}" for i in ids if i != "b"}
        g = audit.build_report(rows, [ds("x")], [audit.ALL], {"lin": groups})["tasks"][0]["eligible"]["groups"]["lin"]
        self.assertAlmostEqual(g["group_coverage"], 0.9)
        self.assertAlmostEqual(g["group_ess"], 9.0)
        self.assertIsNone(g["group_ess_null_reason"])
        self.assertAlmostEqual(g["minority_group_coverage"], 0.5)
        self.assertIsNone(g["minority_group_ess"])
        self.assertEqual(g["minority_group_ess_null_reason"], "1 of 2 rows grouped, under 80%")

    def test_zero_rows(self):
        g, members, minority_members = audit.group_summary([], "R", {"a": "G1"})
        self.assertEqual((g["group_coverage"], g["largest_group_share"], g["group_ess"], g["minority_group_ess"]),
                         (None, None, None, None))
        self.assertEqual((g["group_ess_null_reason"], g["minority_group_ess_null_reason"]), ("no rows", "no rows"))
        self.assertEqual((members, minority_members), (set(), set()))
        # Cap 0 selects nothing, so the step from it has no groups_added either.
        rows = [ph("a", "x", "R"), ph("b", "x", "S")]
        rep = audit.build_report(rows, [ds("x")], [0, audit.ALL], {"lin": {"a": "G1", "b": "G2"}})
        t = rep["tasks"][0]
        self.assertEqual(t["by_cap"]["0"]["groups"]["lin"]["group_ess_null_reason"], "no rows")
        d = t["nested"][0]["groups"]["lin"]
        self.assertEqual((d["groups_added"], d["groups_added_null_reason"]), (None, "cap 0: no rows"))

    def coverage_report(self, groups):
        # x: R={a, b}, S={c, d, e}; minority R. "all" selects every row.
        rows = [ph("a", "x", "R"), ph("b", "x", "R"), ph("c", "x", "S"), ph("d", "x", "S"), ph("e", "x", "S")]
        return audit.build_report(rows, [ds("x")], [audit.ALL], {"lin": groups})

    def coverage_case(self, groups):
        return self.coverage_report(groups)["tasks"][0]["eligible"]["groups"]["lin"]

    def test_all_rows_ungrouped(self):
        g = self.coverage_case({})
        self.assertEqual((g["rows_with_group"], g["rows_without_group"], g["group_coverage"], g["groups"],
                          g["minority_groups"], g["largest_group_rows"], g["largest_group_share"]),
                         (0, 5, 0.0, 0, 0, 0, None))
        self.assertIsNone(g["group_ess"])
        self.assertEqual(g["group_ess_null_reason"], "0 of 5 rows grouped, under 80%")
        self.assertIsNone(g["minority_group_ess"])
        self.assertEqual(g["minority_rows_without_group"], 2)
        # The TSV says so too: no ESS, and the reason beside it.
        head, line = audit.tsv_lines(self.coverage_report({}))
        t = dict(zip(head.split("\t"), line.split("\t")))
        self.assertEqual((t["lin_rows_without_group"], t["lin_group_coverage"], t["lin_groups"], t["lin_group_ess"],
                          t["lin_group_ess_null_reason"]), ("5", "0", "0", "", "0 of 5 rows grouped, under 80%"))

    def test_partly_ungrouped_withholds_ess(self):
        # a, b, c grouped (3 of 5 = 0.6 < 0.8); both minority rows grouped.
        g = self.coverage_case({"a": "G1", "b": "G1", "c": "G2"})
        self.assertEqual((g["rows_without_group"], g["groups"], g["minority_groups"], g["largest_group_rows"]),
                         (2, 2, 1, 2))
        self.assertAlmostEqual(g["group_coverage"], 0.6)
        self.assertAlmostEqual(g["largest_group_share"], 2 / 3)
        self.assertIsNone(g["group_ess"])
        self.assertEqual(g["group_ess_null_reason"], "3 of 5 rows grouped, under 80%")
        self.assertEqual((g["minority_group_coverage"], g["minority_group_ess"]), (1.0, 1.0))
        self.assertIsNone(g["minority_group_ess_null_reason"])

    def test_coverage_at_threshold_reports_ess(self):
        # 4 of 5 grouped is exactly 0.8: sizes 1, 2, 1 -> 16 / 6.
        g = self.coverage_case({"a": "G1", "b": "G2", "c": "G2", "d": "G3"})
        self.assertAlmostEqual(g["group_ess"], 16 / 6)
        self.assertAlmostEqual(g["largest_group_share"], 2 / 4)

    def test_fully_grouped(self):
        # Sizes G1=1 (a), G2=2 (b, c), G3=2 (d, e): ESS 25 / 9; minority a, b in G1, G2.
        g = self.coverage_case({"a": "G1", "b": "G2", "c": "G2", "d": "G3", "e": "G3"})
        self.assertEqual((g["rows_without_group"], g["group_coverage"], g["groups"], g["minority_groups"]),
                         (0, 1.0, 3, 2))
        self.assertAlmostEqual(g["group_ess"], 25 / 9)
        self.assertAlmostEqual(g["largest_group_share"], 2 / 5)
        self.assertAlmostEqual(g["minority_group_ess"], 2.0)

    def test_without_groups(self):
        rows, sets = fixture()
        rep = audit.build_report(rows, sets, [2, 4])
        t = rep["tasks"][0]
        self.assertNotIn("groups", t["eligible"])
        self.assertNotIn("groups", t["by_cap"]["2"])
        self.assertNotIn("groups", t["nested"][0])
        self.assertEqual(rep["group_sets"], [])
        self.assertEqual(len(audit.tsv_lines(rep)), 1 + 2 * 2)

    def test_mic_only_rows_in_report(self):
        rows, sets = fixture()
        rows = rows + [ph("n", "x", "S", from_mic=True), ph("n", "y", "S", from_mic=True)]
        rep = audit.build_report(rows, sets, [audit.ALL])
        sp = rep["species"][0]
        self.assertEqual((sp["scored_genomes"], sp["unscored_genomes"]), (7, 2))  # m and n
        x = rep["tasks"][0]
        self.assertEqual(x["eligible"]["rows"], 6)
        self.assertEqual(x["eligible"]["non_target_rows"]["mic_inferred"], 2)
        # Even with no cap, m and n are not selected, so their rows are not in the cohort.
        self.assertEqual(x["by_cap"]["all"]["non_target_rows"]["mic_inferred"], 0)

    def test_tied_minority(self):
        rows = [ph("a", "x", "R"), ph("b", "x", "S")]
        t = audit.build_report(rows, [ds("x")], [1])["tasks"][0]
        self.assertTrue(t["minority_tied"])
        self.assertEqual(t["minority"], "R")

    def test_caps_parse_and_order(self):
        self.assertEqual(audit.parse_caps("4000,all,2000"), [2000, 4000, audit.ALL])
        with self.assertRaises(ValueError):
            audit.parse_caps("2000,2000")

    def test_cli_writes_json_and_tsv(self):
        rows, sets = fixture()
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "phenotypes.jsonl"), "w") as f:
                f.write("".join(json.dumps(r) + "\n" for r in rows))
            with open(os.path.join(d, "datasets.json"), "w") as f:
                json.dump(sets, f)
            with open(os.path.join(d, "dup.tsv"), "w") as f:
                f.write("genome_id\tgroup_id\n" + "".join(f"{k}\t{v}\n" for k, v in GROUPS.items()))
            out, tsv = os.path.join(d, "out.json"), os.path.join(d, "out.tsv")
            subprocess.run([sys.executable, "-I", os.path.join(LANE, "audit.py"), out, "--curated", d,
                            "--caps", "2,4", "--groups", f"dup={d}/dup.tsv", "--tsv", tsv], check=True)
            with open(out) as f:
                rep = json.load(f)
            with open(tsv) as f:
                lines = f.read().splitlines()
        self.assertEqual(rep["tasks"][0]["eligible"]["groups"]["dup"]["groups"], 3)
        self.assertEqual(len(lines), 1 + 2 * 2)
        head = lines[0].split("\t")
        for col in ("dup_rows_without_group", "dup_group_coverage", "dup_group_ess", "dup_group_ess_null_reason"):
            self.assertIn(col, head)
        # Task x at cap 4 (a, b, c, e): all grouped. Task y at cap 2 (c, a): all grouped.
        x4 = dict(zip(head, lines[2].split("\t")))
        self.assertEqual((x4["antibiotic"], x4["cap"], x4["dup_rows_without_group"], x4["dup_group_coverage"],
                          x4["dup_group_ess_null_reason"]), ("x", "4", "0", "1", ""))

    def test_groups_from_metadata_round_trips_through_read_groups(self):
        meta = [{"genome_id": "a", "mlst": "MLST.x.1"}, {"genome_id": "b", "mlst": ""},
                {"genome_id": "c"}, {"genome_id": "d", "mlst": "MLST.x.1"}]
        with tempfile.TemporaryDirectory() as d:
            src, out = os.path.join(d, "genomes.jsonl"), os.path.join(d, "mlst.tsv")
            with open(src, "w") as f:
                f.write("".join(json.dumps(m) + "\n" for m in meta))
            subprocess.run([sys.executable, "-I", os.path.join(LANE, "groups_from_metadata.py"), out, src, "mlst"],
                           check=True, capture_output=True)
            self.assertEqual(audit.read_groups(out), {"a": "mlst:MLST.x.1", "d": "mlst:MLST.x.1"})

    def test_groups_file_rejects_two_groups_for_one_genome(self):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "g.tsv")
            with open(path, "w") as f:
                f.write("a\tG1\na\tG2\n")
            with self.assertRaises(ValueError):
                audit.read_groups(path)


if __name__ == "__main__":
    unittest.main()
