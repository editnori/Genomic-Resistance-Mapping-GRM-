"""Plant each fault of stage 2 of the matched accuracy protocol and run the test meant to fail with it.

  python3 -I stage2_faults.py <out dir> [fault ...]

The runner is stage 1's, scripts/matched-accuracy/planted_faults.py, given this table: it copies this checkout's
scripts/, tests/, reference/ and docs/ to <out dir>/tree, plants one fault at a time in
scripts/matched-accuracy/stage2.py there, runs the one test named for it as CI runs the suite, and puts the text
back. Each test must pass with no fault planted. A fault counts only when its test ran alone and failed by an
assertion. Writes <out dir>/faults.tsv and <out dir>/logs/, prints one line per fault and a last line with the
count, and exits 1 when a fault does not fail its test by an assertion. Every test it runs is SYNTHETIC or reads
committed files; it reads no data and no score.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
S = "scripts/matched-accuracy/stage2.py"
T = "tests.matched_accuracy."
MEASURE = T + "test_stage2_measure."
MADE_UP = T + "test_stage2_made_up."
STOPS = T + "test_stage2_stops."
LINE = T + "test_stage2_line."
EXPORT = T + "test_stage2_export."

# fault, rule, what it plants, file, text, planted text, test
FAULTS = [
    ("S01", "row mean", "the five scores added by sum, compensated, in place of left to right", S,
     "    return functools.reduce(operator.add, [component(x) for x in five]) / 5\n",
     "    return sum([component(x) for x in five]) / 5\n",
     MEASURE + "RowMean.test_the_five_scores_are_added_left_to_right_then_divided_by_5"),
    ("S02", "row mean", "the five scores added by math.fsum in place of left to right", S,
     "    return functools.reduce(operator.add, [component(x) for x in five]) / 5\n",
     "    return math.fsum([component(x) for x in five]) / 5\n",
     MEASURE + "RowMean.test_a_row_of_the_five_values_loses_to_a_row_of_five_halves"),
    ("S03", "tie", "an exact tie counted 1", S,
     "    return 2 if r > s else 1 if r == s else 0\n",
     "    return 2 if r > s else 2 if r == s else 0\n",
     MADE_UP + "Designed.test_the_one_tie_counts_one_half"),
    ("S04", "tie", "an exact tie counted 0", S,
     "    return 2 if r > s else 1 if r == s else 0\n",
     "    return 2 if r > s else 0 if r == s else 0\n",
     MEASURE + "Pairs.test_a_higher_resistant_mean_counts_1_a_tie_one_half_and_lower_0"),
    ("S05", "weights", "pairs pooled across the groups of a dataset", S,
     "        dataset_value[d] = mean([fractions.Fraction(h, 2 * p) for h, p, _ in gs])\n",
     "        dataset_value[d] = fractions.Fraction(sum(h for h, _, _ in gs), 2 * sum(p for _, p, _ in gs))\n",
     MEASURE + "Weights.test_pairs_are_not_pooled_across_groups"),
    ("S06", "weights", "groups pooled across the datasets of a species", S,
     "    species_value = {sp: mean([dataset_value[d] for d in ds]) for sp, ds in datasets_of.items()}\n",
     "    species_value = {sp: mean([fractions.Fraction(h, 2 * p) for d in ds for h, p, _ in by_dataset[d]])\n"
     "                     for sp, ds in datasets_of.items()}\n",
     MEASURE + "Weights.test_groups_are_not_pooled_across_datasets"),
    ("S07", "weights", "datasets pooled across species", S,
     "    overall = mean([species_value[sp] for sp in species_value]) if species_value else None\n",
     "    overall = mean([dataset_value[d] for d in dataset_value]) if species_value else None\n",
     MEASURE + "Weights.test_datasets_are_not_pooled_across_species"),
    ("S08", "interval", "a seed other than 42", S,
     "SEED = 42\n",
     "SEED = 7\n",
     MEASURE + "Interval.test_the_bounds_are_those_of_pcg64_seed_42_and_the_linear_percentiles"),
    ("S09", "interval", "the nearest percentile in place of the linear one", S,
     "    lower, upper = numpy.percentile(draws, PERCENTILES)\n",
     "    lower, upper = numpy.percentile(draws, PERCENTILES, method=\"nearest\")\n",
     MEASURE + "Interval.test_the_bounds_are_those_of_pcg64_seed_42_and_the_linear_percentiles"),
    ("S10", "interval", "species drawn in the order they are met, not sorted by name", S,
     "    names = sorted(lv[\"species\"])\n",
     "    names = list(lv[\"species\"])\n",
     MEASURE + "Interval.test_species_are_drawn_in_sorted_order_of_their_names"),
    ("S11", "interval", "an interval from 5 species", S,
     "if len(names) >= MIN_SPECIES else None",
     "if len(names) >= MIN_SPECIES - 1 else None",
     MADE_UP + "FiveSpecies.test_5_species_give_the_estimate_and_no_interval"),
    ("S12", "the line", "a score read before the line", S,
     "    hc20, hc20_ids, hc20_rows, hc20_missing = read_hc20(args.hc20[0], wanted)\n",
     "    hc20, hc20_ids, hc20_rows, hc20_missing = read_hc20(args.hc20[0], wanted)\n"
     "    pick(json.loads(Path(files[0][1]).read_bytes()), [])\n",
     LINE + "Order.test_every_check_runs_before_the_line_and_every_score_after_it"),
    ("S13", "the line", "the line printed before the checks", S,
     "    checked = check(args)\n    opened = open_scores(args.out)\n",
     "    opened = open_scores(args.out)\n    checked = check(args)\n",
     LINE + "Order.test_only_check_calls_the_checks_and_run_calls_check_before_the_line_and_measure_after_it"),
    ("S14", "ten pairs", "the 10-pair rule taken off", S,
     "    if s is None or s[\"value\"] is None or (level != \"overall\" and s[\"pairs\"] < MIN_PAIRS):\n",
     "    if s is None or s[\"value\"] is None:\n",
     MADE_UP + "Designed.test_a_dataset_under_10_pairs_gives_counts_and_no_estimate_and_enters_every_mean"),
    ("S15", "export", "an estimate written with a decimal point", S,
     "    return round(fractions.Fraction(value) * PPM)\n",
     "    return round(fractions.Fraction(value) * PPM) / 1\n",
     EXPORT + "Export.test_no_export_file_holds_a_decimal_point"),
    ("S16", "export", "the decimal point check taken off", S,
     "        if b\".\" in data:\n",
     "        if b\".\" in data and False:\n",
     EXPORT + "Scan.test_a_decimal_point_in_an_export_file_stops_after_the_line"),
    ("S17", "typing", "the typing support check taken off", S,
     "        check_typing(rows[name], groups[name], hc20, typing_table, name)\n",
     "",
     STOPS + "Typing.test_an_hc20_group_other_than_the_one_stage_1_read_stops"),
    ("S18", "typing", "the typing support checked for the MLST read alone", S,
     "    for step, _, _ in TYPING:\n        for level in LEVELS:\n",
     "    for step, _, _ in TYPING[:1]:\n        for level in LEVELS:\n",
     STOPS + "Typing.test_typing_support_that_differs_from_stage_1s_typing_table_at_any_level_stops"),
    ("S19", "counts", "the counts checked at the dataset level alone", S,
     "    zero = dict.fromkeys(SUPPORT, 0)\n    for level in LEVELS:\n",
     "    zero = dict.fromkeys(SUPPORT, 0)\n    for level in LEVELS[2:]:\n",
     STOPS + "Groups.test_counts_that_differ_from_stage_1s_count_table_at_any_level_stop"),
    ("S20", "groups", "a group with one lab outcome passes", S,
     "    if rest:\n",
     "    if rest and False:\n",
     STOPS + "Groups.test_a_group_with_one_lab_outcome_stops"),
    ("S21", "row list", "lines out of stage 1's sort pass", S,
     "        if i and line < lines[i - 1]:\n",
     "        if i and line < lines[i - 1] and False:\n",
     STOPS + "Lines.test_lines_out_of_stage_1s_sort_stop"),
    ("S22", "rows", "a listed row whose lineage cluster differs passes", S,
     "        if s[\"lineage\"] != r[\"lineage\"]:\n",
     "        if s[\"lineage\"] != r[\"lineage\"] and False:\n",
     STOPS + "Rows.test_a_listed_row_whose_lineage_cluster_differs_from_its_reports_stops"),
    ("S23", "learners", "a report with the five in another order passes", S,
     "        if [n for n in names if n in FIVE] != list(FIVE):\n",
     "        if sorted(n for n in names if n in FIVE) != sorted(FIVE):\n",
     STOPS + "Learners.test_a_listed_report_with_the_five_in_another_order_stops"),
    ("S24", "null fits", "the null-fit read looks at scm, cart and rf alone", S,
     "            if not any(complexity[g[0][\"dataset\"]][n][k[1]] == 0 for n in FIVE)}\n",
     "            if not any(complexity[g[0][\"dataset\"]][n][k[1]] == 0 for n in FIVE[:3])}\n",
     MADE_UP + "Designed.test_the_null_fit_read_leaves_out_xis_fold_4_group"),
    ("S25", "chance", "the word beside an interval taken from the rounded bounds", S,
     "    return \"above\" if lower > 0.5 else \"below\" if upper < 0.5 else \"unresolved\"\n",
     "    return \"above\" if ppm(lower) > PPM // 2 else \"below\" if ppm(upper) < PPM // 2 else \"unresolved\"\n",
     MEASURE + "Chance.test_a_lower_bound_just_above_one_half_that_rounds_to_500000_gives_above"),
]


def _runner():
    spec = importlib.util.spec_from_file_location("matched_accuracy_faults_of_stage2", HERE / "planted_faults.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(argv):
    runner = _runner()
    runner.FAULTS = FAULTS
    runner.__doc__ = __doc__
    return runner.main(argv)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
