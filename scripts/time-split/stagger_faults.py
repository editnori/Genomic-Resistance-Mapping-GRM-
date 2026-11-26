#!/usr/bin/env python3
"""Plant each fault of the staggered split's rules and run the test that fails with it.

  python3 -I stagger_faults.py <out dir> [fault ...]
  python3 -I stagger_faults.py --list

Copies this checkout's last commit to <out dir>/tree (git archive) and, one fault at a time (every fault, or those
named), changes one piece of text there, runs the one test named for it, and puts the text back. A fault is caught
when its test fails; each test must also pass with no fault planted, and each run must run exactly one test. The
out dir must be new or empty: nothing is deleted. An out dir whose path holds `latetest`, `sealed` or `held-out`,
as given or where it leads, stops it before it makes any directory, with the stop line of the other two scripts.
Writes <out dir>/faults.tsv (fault, rule, file, test, result) and <out dir>/logs/, prints one line per fault and
`caught X of Y`, and exits 1 when a fault is not caught. Every test it runs is SYNTHETIC; it reads no data.
"""

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
_spec = importlib.util.spec_from_file_location("stagger_design_for_faults", HERE / "stagger_design.py")
design_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(design_mod)

DESIGN = "scripts/time-split/stagger_design.py"
ANALYZE = "scripts/time-split/stagger_analyze.py"
RUNNER = "scripts/time-split/stagger_faults.py"
TD = "tests.time_split.test_stagger_design."
TA = "tests.time_split.test_stagger_analyze."
TF = "tests.time_split.test_stagger_faults."

# (fault, rule it breaks, file, text, planted text, test)
FAULTS = [
    ("cut-floor", "R33 point 1: the cut point at place ceil(s x n)", DESIGN,
     "return [order[-(-n * s // 100) - 1] for s in shares]",
     "return [order[n * s // 100 - 1] for s in shares]",
     TD + "CutPoints.test_seven_dates_round_the_place_up"),
    ("window-edge", "item 2: window j runs after c(j) up to c(j+1)", DESIGN,
     "        if date <= cuts[j]:\n            return j",
     "        if date < cuts[j]:\n            return j",
     TD + "Windows.test_after_one_cut_up_to_the_next_the_last_to_the_end"),
    ("study-known-late", "item 3: a study is known when a genome carries it by c(j)", DESIGN,
     "date[g] <= cuts[j - 1]}",
     "date[g] < cuts[j - 1]}",
     TD + "HandWorld.test_a_study_known_by_an_origin"),
    ("no-study-tested", "item 3: a genome with no study is tested nowhere", DESIGN,
     "== j and study.get(g) and study[g] not in known]",
     "== j and study.get(g) not in known]",
     TD + "HandWorld.test_a_genome_with_no_study_is_tested_nowhere"),
    ("undated-trains", "item 1: an undated genome has no training set", DESIGN,
     'arms[("all", k)] = [g for g in dated if date[g] <= cuts[k - 1] and g not in copies]',
     'arms[("all", k)] = [g for g in cohort if (date[g] or cuts[0]) <= cuts[k - 1] and g not in copies]',
     TD + "HandWorld.test_an_undated_genome_trains_nowhere_and_is_tested_nowhere"),
    ("undated-knows-study", "item 3: an undated genome makes no study known", DESIGN,
     "date.get(g) is not None and date[g] <= cuts[j - 1]}",
     "(date.get(g) is None or date[g] <= cuts[j - 1])}",
     TD + "HandWorld.test_an_undated_genome_makes_no_study_known"),
    ("year-not-four-digits", "item 1: the collection year is a four-digit year", DESIGN,
     "return int(text) if text is not None and YEAR.fullmatch(text) else None",
     "return int(text[:4]) if text is not None and YEAR.match(text) else None",
     TD + "Clocks.test_the_collection_year_is_four_digits"),
    ("one-clock-dates", "item 1: each clock reads its own dates", DESIGN,
     '"dates": {"public-day": day, "collection-year": year}',
     '"dates": {"public-day": day, "collection-year": day}',
     TD + "Clocks.test_each_clock_its_own_dates_and_hash"),
    ("tie-hash-no-clock", "item 6: the keyed hash to the letter", DESIGN,
     'hashlib.sha256(("stagger-20261125|" + clock + "|" + dataset + "|" + genome).encode())',
     'hashlib.sha256(("stagger-20261125|" + dataset + "|" + genome).encode())',
     TD + "Ties.test_the_plans_keyed_hash_to_the_letter"),
    ("tie-larger-hash", "item 6: a tie goes to the smaller keyed hash", DESIGN,
     'members = sorted((g for g in pool if labels[g] == c), key=key)',
     'members = sorted((g for g in pool if labels[g] == c), key=key, reverse=True)',
     TD + "Ties.test_a_tie_goes_to_the_smaller_hash_of_its_clock"),
    ("recent-earliest", "item 4: recent(k) takes the latest", DESIGN,
     "members.sort(key=lambda g: date[g], reverse=True)",
     "members.sort(key=lambda g: date[g])",
     TD + "HandWorld.test_recent_k_takes_the_latest"),
    ("recent-size", "R33 point 4: recent(k) holds as many of each class as all(1)", DESIGN,
     'size = dict(zip("RS", counts(arms[("all", 1)], labels)))',
     'size = dict(zip("RS", counts(arms[("all", j)], labels)))',
     TD + "HandWorld.test_recent_k_holds_as_many_of_each_class_as_all_1"),
    ("copies-train", "item 7 and R33 point 3: copies of a window's test genomes train nowhere on it", DESIGN,
     "copies = {g for g in cohort if g not in tested and digests[g] in seen}",
     "copies = set()",
     TD + "HandWorld.test_copies_per_window"),
    ("same-genomes-pair", "a pair needs two training sets that differ", DESIGN,
     "            if set(older) == set(fresh):\n                continue\n",
     "",
     TD + "Pairs.test_two_training_sets_of_the_same_genomes_make_no_pair"),
    ("test-minimum-one-class", "R33 point 5: the test minimum of each class", DESIGN,
     'if min(counts(w["test"], labels)) < test_minimum:',
     'if counts(w["test"], labels)[0] < test_minimum:',
     TD + "Pairs.test_a_test_set_under_the_minimum_in_either_class"),
    ("train-minimum-24", "R33 point 5: 25 of each class in a training set", DESIGN,
     "if min(counts(older, labels)) < MIN_TRAIN_CLASS:",
     "if min(counts(older, labels)) < MIN_TRAIN_CLASS - 1:",
     TD + "Pairs.test_a_training_set_under_25_in_either_class"),
    ("one-arm-both-contrasts", "item 7: stale same size compares recent(k) with recent(j)", DESIGN,
     'ARM_OF = {"stale-all": "all", "stale-same-size": "recent"}',
     'ARM_OF = {"stale-all": "all", "stale-same-size": "all"}',
     TD + "Pairs.test_each_contrast_reads_its_own_arms"),
    ("rule-step-1-or", "variant rule step 1: 12 datasets in 3 species", DESIGN,
     'if e["datasets"] >= MIN_DATASETS and e["species"] >= MIN_SPECIES]',
     'if e["datasets"] >= MIN_DATASETS or e["species"] >= MIN_SPECIES]',
     TD + "VariantRule.test_step_1_needs_12_datasets_and_3_species"),
    ("rule-step-2-fewest", "variant rule step 2: the most origins", DESIGN,
     'most = max(len(e["shares"]) for e in kept)',
     'most = min(len(e["shares"]) for e in kept)',
     TD + "VariantRule.test_step_2_the_most_origins"),
    ("rule-step-3-lowest", "variant rule step 3: the highest test minimum", DESIGN,
     'highest = max(e["minimum"] for e in kept)',
     'highest = min(e["minimum"] for e in kept)',
     TD + "VariantRule.test_step_3_the_highest_test_minimum"),
    ("rule-step-4-fewest", "variant rule step 4: the most datasets", DESIGN,
     'most = max(e["datasets"] for e in kept)',
     'most = min(e["datasets"] for e in kept)',
     TD + "VariantRule.test_step_4_the_most_datasets"),
    ("rule-tie-last", "variant rule: a tie goes to the first in the plan's order", DESIGN,
     "taken = kept[0]",
     "taken = kept[-1]",
     TD + "VariantRule.test_a_tie_goes_to_the_first_in_the_plans_order"),
    ("rule-fallback-fewest", "variant rule: none at 12 in 3, the most datasets", DESIGN,
     'most = max(e["datasets"] for e in entries)',
     'most = min(e["datasets"] for e in entries)',
     TD + "VariantRule.test_none_reaches_12_in_3_the_most_datasets"),
    ("rule-fallback-not-descriptive", "variant rule: none at 12 in 3, descriptive only", DESIGN,
     "descriptive = not kept",
     "descriptive = False",
     TD + "VariantRule.test_none_reaches_12_in_3_the_clock_is_descriptive_only"),
    ("rule-reads-same-size", "R33 point 6: the rule counts stale all at horizon 1", DESIGN,
     'CHOICE_CONTRAST, CHOICE_HORIZON = "stale-all", 1',
     'CHOICE_CONTRAST, CHOICE_HORIZON = "stale-same-size", 1',
     TD + "VariantRule.test_the_rule_reads_stale_all_at_horizon_1"),
    ("check60-share", "the 60% check cuts at 60%", DESIGN,
     "CHECK_SHARE = 60",
     "CHECK_SHARE = 50",
     TD + "Check60.test_the_cut_days_of_the_table"),
    ("check60-passes-all", "the 60% check fails on a cut day other than the record's", DESIGN,
     'row["same_cut_day"] = int(cut == ref["cut_day"] and row["undated"] == 0)',
     'row["same_cut_day"] = int(row["undated"] == 0)',
     TD + "Check60.test_a_cut_day_other_than_the_tables"),
    ("barred-word-dropped", "no path holding latetest, sealed or held-out is opened", DESIGN,
     'BARRED = ("latetest", "sealed", "held-out")',
     'BARRED = ("latetest", "held-out")',
     TD + "BarredPaths.test_each_word_and_role"),
    ("barred-link-followed", "a path that leads to a barred place is not opened", DESIGN,
     "texts.append(str(p.resolve()).lower())",
     "texts.append(str(p).lower())",
     TD + "BarredPaths.test_a_link_that_leads_into_a_sealed_place"),
    ("recent-1-own-file", "recent(1) is all(1): one file", DESIGN,
     'if k == 1 and arm == "recent":',
     "if False:",
     TD + "Folds.test_recent_1_and_all_1_share_one_file"),
    ("sha256-other-bytes", "each fold file named with the sha256 of its bytes", DESIGN,
     '"sha256": hashlib.sha256(data).hexdigest()}',
     '"sha256": hashlib.sha256(json.dumps(body).encode()).hexdigest()}',
     TD + "Folds.test_each_file_named_with_the_sha256_of_its_bytes"),
    ("pair-tests-other-genomes", "both fits of a pair test the same genomes", DESIGN,
     'body = fold_body(spec["id"], scheme, cohort, w["test"], w["arms"][base])',
     'body = fold_body(spec["id"], scheme, cohort, w["test"][:len(w["test"]) - (base[0] == "recent")], '
     'w["arms"][base])',
     TD + "Folds.test_both_fits_of_a_pair_test_the_same_genomes"),
    ("tested-also-unused", "every cohort genome once in a fold file", DESIGN,
     '"unused": [g for g in cohort if g not in tested and g not in trained]',
     '"unused": [g for g in cohort if g not in trained]',
     TD + "Folds.test_every_cohort_genome_once_in_each_fold_file"),
    ("design-names-genomes", "the design holds no genome id", DESIGN,
     '"copies": len(w["copies"]),',
     '"copies": sorted(w["copies"]),',
     TD + "Folds.test_the_design_holds_no_genome_id"),
    ("unbuilt-folds", "no fold file for a species not built", DESIGN,
     'elif not sp["built"]:',
     "elif False:",
     TD + "Folds.test_a_species_not_built_gets_no_fold_file"),
    ("out-dir-overwritten", "an out dir that holds a file is not written", DESIGN,
     "need(not out.exists() or (out.is_dir() and not any(out.iterdir())),",
     "need(True,",
     TD + "Folds.test_an_out_dir_that_holds_a_file_is_refused"),
    ("check60-count-not-24", "R38: the 60% check stops unless its table has 24 rows", DESIGN,
     "need(len(reference) == datasets,",
     "need(True,",
     TD + "Check60.test_a_table_of_23_datasets_stops"),
    ("check60-count-of-the-test", "R38: the 24 is one line of the script, and main gives it", DESIGN,
     "CHECK_DATASETS = 24",
     "CHECK_DATASETS = 2",
     TD + "Check60.test_the_count_is_24_whatever_the_table"),
    ("check60-dataset-twice", "R38: each row of the 60% check's table a different dataset", DESIGN,
     "need(not twice,",
     "need(True,",
     TD + "Check60.test_one_dataset_twice_and_one_left_out_stops"),
    ("check60-dataset-not-in-metadata", "R38: each row of the 60% check's table a dataset of the metadata", DESIGN,
     "need(not lacking,",
     "need(True,",
     TD + "Check60.test_a_dataset_the_metadata_lacks_stops"),
    ("check60-undated-ignored", "the 14 Nov rule: every cohort genome dated (the reviewer's fault)", DESIGN,
     'int(cut == ref["cut_day"] and row["undated"] == 0)',
     'int(cut == ref["cut_day"])',
     TD + "Check60.test_an_undated_cohort_genome_fails_its_dataset"),
    ("count-pairs-as-datasets", "a dataset enters at a horizon when it has a pair there (the reviewer's fault)", DESIGN,
     '"datasets": len(hits),',
     '"datasets": sum(n for _, _, n in hits),',
     TD + "Counts.test_a_dataset_with_two_pairs_at_a_horizon_counts_once"),
    ("species-as-datasets", "12 datasets in 3 species: species apart from datasets (the reviewer's fault)", DESIGN,
     '"species": len({s["species_slug"] for s, _, _ in hits}),',
     '"species": len(hits),',
     TD + "Counts.test_two_datasets_of_one_species_count_one_species"),
    ("result-sha-unchecked", "item 5: each result names its fold file's sha256", ANALYZE,
     'r.get("fold_file_sha256") == fit["sha256"]',
     'r.get("fold_file_sha256") is not None',
     TA + "Stops.test_a_result_that_does_not_name_its_fold_files_sha256"),
    ("pair-genomes-unchecked", "item 5: both fits of a pair scored the same genomes", ANALYZE,
     "need([g for g, _ in older[0]] == [g for g, _ in fresh[0]],",
     "need(len(older[0]) == len(fresh[0]),",
     TA + "Stops.test_two_fits_that_scored_other_genomes"),
    ("pair-labels-unchecked", "item 5: with the same labels", ANALYZE,
     'need(older[0] == fresh[0], f"',
     'need(True, f"',
     TA + "Stops.test_the_same_genomes_with_other_labels"),
    ("train-count-unchecked", "item 5: the training counts equal the design's", ANALYZE,
     'need(fold["train"] == fit["train"] and fold["train_resistant"] == fit["train_resistant"],',
     'need(fold["train_resistant"] == fit["train_resistant"],',
     TA + "Stops.test_training_counts_other_than_the_designs"),
    ("recall-loose", "item 5: recall worked out again equals the evaluator's", ANALYZE,
     'abs(m["recall"] - e["sensitivity"]) < 1e-12',
     'abs(m["recall"] - e["sensitivity"]) < 0.1',
     TA + "Stops.test_recall_worked_out_again_other_than_the_evaluators"),
    ("missing-result-unchecked", "item 5: a fit with no result file stops the reading", ANALYZE,
     'need(path.exists(), f"{name}: no result file")',
     "pass",
     TA + "Stops.test_a_missing_result_file"),
    ("pair-mean-first", "a dataset with several pairs at a horizon gives their mean", ANALYZE,
     '"d": {m: mean([o[m] - f[m] for o, f in zip(older, fresh)]) for m in MEASURES}}',
     '"d": {m: older[0][m] - fresh[0][m] for m in MEASURES}}',
     TA + "Values.test_a_dataset_with_two_pairs_at_a_horizon_gives_their_mean"),
    ("interval-ends", "the 251st and 9,750th of 10,000 sorted draws", ANALYZE,
     "return point, draws[LOW_DRAW - 1], draws[HIGH_DRAW - 1]",
     "return point, draws[LOW_DRAW], draws[HIGH_DRAW]",
     TA + "Interval.test_the_251st_and_the_9750th_of_the_sorted_draws"),
    ("seed", "seed 20261125", ANALYZE,
     "RESAMPLES, SEED = 10000, 20261125",
     "RESAMPLES, SEED = 10000, 20261126",
     TA + "Interval.test_the_seed_and_the_order_of_draws"),
    ("generator-per-interval", "one generator a clock and reading, drawn in order", ANALYZE,
     "point, lo, hi = interval(by_species, rng)",
     "point, lo, hi = interval(by_species, random.Random(SEED))",
     TA + "Interval.test_the_seed_and_the_order_of_draws"),
    ("descriptive-clock-ignored", "a clock the rule marks descriptive only", ANALYZE,
     'descriptive_clock = design["clocks"][clock]["descriptive_only"]',
     "descriptive_clock = False",
     TA + "DescriptiveClock.test_a_clock_the_rule_marks_descriptive_gives_descriptive_only"),
    ("read-sign", "worse, better or not shown, by each measure's side", ANALYZE,
     'return "worse" if (lo > 0) == (WORSE[m] > 0) else "better"',
     'return "worse" if lo > 0 else "better"',
     TA + "Readings.test_each_reading_of_an_interval"),
    ("learner-rule-recall-only", "a learner's cost through recall or log loss", ANALYZE,
     'return "cost" if "worse" in got else "not shown"',
     'return "cost" if got[0] == "worse" else "not shown"',
     TA + "Readings.test_each_rule_for_a_learner"),
    ("verdict-two-costs", "cost shown needs at least 3 of 5 learners", ANALYZE,
     ">= MIN_LEARNERS_WITH_COST:",
     ">= MIN_LEARNERS_WITH_COST - 1:",
     TA + "Readings.test_each_verdict"),
    ("verdict-species-dropped", "descriptive only under 3 species", ANALYZE,
     "if descriptive_clock or datasets < MIN_DATASETS or species < MIN_SPECIES:",
     "if descriptive_clock or datasets < MIN_DATASETS:",
     TA + "Readings.test_each_verdict"),
    ("verdict-better-ignored", "any learner better makes the contrast mixed", ANALYZE,
     '        return "mixed"\n    if sum(',
     '        return "cost shown"\n    if sum(',
     TA + "Readings.test_each_verdict"),
    ("no-list-dropped", "R33 point 9: a pair whose test set lists no genome stays", ANALYZE,
     'if not on_list:\n        return set(range(len(samples))), ""',
     "if not on_list:\n        return None, SHORT",
     TA + "SecondCase.test_no_genome_listed"),
    ("list-minimum-4", "R33 point 9: at least 5 of each class listed", ANALYZE,
     "LIST_MINIMUM = 5",
     "LIST_MINIMUM = 4",
     TA + "SecondCase.test_one_to_four_of_a_class"),
    ("left-under-minimum", "R33 point 9: each class keeps the test minimum once listed genomes are out", ANALYZE,
     "if min(r, s) < test_minimum:",
     "if min(r, s) < 0:",
     TA + "SecondCase.test_five_or_more_but_a_class_falls_under_the_minimum"),
    ("listed-out-of-fresh-only", "R33 point 8: listed genomes out for both fits", ANALYZE,
     "older = [measure(o[0], o[1][l], keep) for o, _, keep in kept]",
     "older = [measure(o[0], o[1][l]) for o, _, keep in kept]",
     TA + "SecondReading.test_the_listed_genomes_leave_both_fits"),
    ("short-list-counted", "R33 point 9: a short-list pair gives no figure", ANALYZE,
     "                        if keep is None:\n"
     "                            reasons.append(why)\n"
     "                        else:\n"
     "                            kept.append((older, fresh, keep))",
     "                        if keep is None and why != SHORT:\n"
     "                            reasons.append(why)\n"
     "                        else:\n"
     "                            kept.append((older, fresh, keep if keep is not None else set(range(len(fresh[0])))))",
     TA + "SecondReading.test_a_short_list_leaves_the_pair_out_with_no_figure"),
    ("short-list-figure-written", "R33 point 9: the lacking table gives a reason and no figure", ANALYZE,
     '"reason": "; ".join(r for r in (SHORT, UNDER) if r in reasons)})',
     '"reason": "; ".join(reasons) + f" ({len(listed)} listed)"})',
     TA + "SecondReading.test_a_short_list_leaves_the_pair_out_with_no_figure"),
    ("id-scan-off", "item 7: no genome id leaves the run", ANALYZE,
     'need(not tokens & ids, f"{name}: a genome id stands in it; nothing is written")',
     'need(True, f"{name}: a genome id stands in it; nothing is written")',
     TA + "Scan.test_a_genome_id_stops_the_writing"),
    ("shape-scan-short", "item 7: no accession or sequence leaves the run", ANALYZE,
     "for shape in SHAPES:",
     "for shape in SHAPES[3:]:",
     TA + "Scan.test_the_shape_of_an_id_an_accession_or_a_sequence"),
    ("p1-three", "prediction 1: not shown for at least 4 of 5", ANALYZE,
     "for l in LEARNERS) >= 4",
     "for l in LEARNERS) >= 3",
     TA + "Predictions.test_p1"),
    ("p1-descriptive-met", "a prediction is not met on a descriptive-only contrast", ANALYZE,
     'if v is None or v == "descriptive only":',
     "if v is None:",
     TA + "Predictions.test_p1"),
    ("p2-any-cost", "prediction 2: a cost through log loss", ANALYZE,
     'rules[(clock, contrast, h, l)] == "cost" and reads[(clock, contrast, h, l)][1] == "worse"',
     'rules[(clock, contrast, h, l)] == "cost"',
     TA + "Predictions.test_p2"),
    ("p3-descriptive-met", "R38: prediction 3 needs a verdict at horizon 1 on both contrasts", ANALYZE,
     '        if any(got[(c, 1)] == "descriptive only" for c in CONTRASTS):\n            return False\n',
     "",
     TA + "Predictions.test_p3_needs_a_verdict_at_horizon_1_on_both_contrasts"),
    ("p3-horizon-1-only", "R38: prediction 3 read at every horizon the clock has", ANALYZE,
     "return [(c, k) for c in CONTRASTS for k in horizons[clock]]",
     "return [(c, k) for c in CONTRASTS for k in (1,)]",
     TA + "Predictions.test_p3_a_cost_at_a_later_horizon_is_not_met"),
    ("p3-mixed-later-met", "R38: a mixed verdict at any horizon, prediction 3 not met", ANALYZE,
     'return all(v in ("not shown", "descriptive only") for v in got.values())',
     'return all(v != "cost shown" for v in got.values())',
     TA + "Predictions.test_p3_a_mixed_verdict_at_a_later_horizon_is_not_met"),
    ("p3-later-descriptive-fails", "R38: a horizon that is descriptive only gives no verdict", ANALYZE,
     'return all(v in ("not shown", "descriptive only") for v in got.values())',
     'return all(v == "not shown" for v in got.values())',
     TA + "Predictions.test_p3_a_descriptive_later_horizon_gives_no_verdict"),
    ("p3-no-verdict-unnamed", "R38: the predictions table names each horizon with no verdict", ANALYZE,
     'unread += [f"{reading} {c} horizon {k}" for c, k in no_verdict(pid, v, horizons)]',
     'unread += [f"{reading} {c} horizon {k}" for c, k in no_verdict(pid, v, horizons) if k == 1]',
     TA + "Predictions.test_the_table_names_each_horizon_with_no_verdict"),
    ("prediction-first-only", "R38: the outcome reads both readings, as stands.tsv does", ANALYZE,
     '"outcome": outcome(got),',
     '"outcome": got["first"],',
     TA + "Agreement.test_without_lists_the_prediction_is_not_read_and_no_verdict_stands"),
    ("prediction-one-reading-met", "R38: met only when both readings are made and both meet it", ANALYZE,
     'if all(got.get(r) == "met" for r in READINGS):',
     'if all(x == "met" for x in got.values()):',
     TA + "Predictions.test_met_only_when_both_readings_are_made_and_both_meet_it"),
    ("prediction-failed-unread", "R38: not met when a reading that was made fails it", ANALYZE,
     '    if "not met" in got.values():\n        return "not met"\n',
     "",
     TA + "Predictions.test_not_met_when_a_reading_made_fails_it"),
    ("prediction-unread-not-met", "R38: otherwise not read", ANALYZE,
     '    return "not read"\n',
     '    return "not met"\n',
     TA + "Predictions.test_otherwise_not_read"),
    ("list-names-no-dataset", "R38: each list file names a dataset of the design", ANALYZE,
     "need(not strays,",
     "need(True,",
     TA + "Lists.test_a_list_file_that_names_no_dataset_of_the_design_stops"),
    ("lists-genome-count", "R38: it prints the list files and the datasets they name, and no genome count", ANALYZE,
     "they name {len(set(lists) & named)} datasets",
     "they name {sum(map(len, lists.values()))} datasets",
     TA + "Lists.test_it_prints_the_files_and_the_datasets_they_name"),
    ("point-mean", "the median over datasets (the reviewer's fault)", ANALYZE,
     "    point = statistics.median([v for s in species for v in by_species[s]])",
     "    point = statistics.mean([v for s in species for v in by_species[s]])",
     TA + "Interval.test_the_point_is_the_median_over_datasets"),
    ("second-left-from-totals", "R33 point 9: the test minimum once listed genomes are out (the reviewer's fault)",
     ANALYZE,
     "    r = sum(1 for _, t in samples if t) - lr",
     "    r = sum(1 for _, t in samples if t)",
     TA + "SecondCase.test_the_resistant_class_alone_falls_under_the_minimum_once_out"),
    ("same-size-fresh-all", "stale same size: recent(k) against recent(j) (the reviewer's fault)", ANALYZE,
     'loaded[fits[(j, j, arm)]["file"]]',
     'loaded[fits[(j, j, "all")]["file"]]',
     TA + "SameSize.test_stale_same_size_scores_against_recent_j"),
    # The runner's own text is written in two pieces here, so that it stands once in this file.
    ("faults-out-unguarded", "R38: a barred out dir stops the runner before it makes any directory", RUNNER,
     'design_mod.guard(argv[0], "out' ' dir")',
     "Path(argv[0])",
     TF + "OutDir.test_a_barred_out_dir_stops_before_any_directory_is_made"),
    ("stands-first-only", "R33 point 10: a verdict stands only when both readings give it", ANALYZE,
     '"stands": first if first == second else (',
     '"stands": first if True else (',
     TA + "Stands.test_a_verdict_stands_only_when_both_readings_give_it"),
]


def run_test(tree, name, log):
    """As CI runs the suite: discover from the top, so `tests` imports under -P. (passed, ran exactly one test)."""
    cmd = ["python3", "-s", "-P", "-B", "-W", "error", "-m", "unittest", "discover", "-s", "tests/time_split", "-t",
           ".", "-k", name]
    env = {"PATH": "/usr/bin:/bin", "HOME": os.environ.get("HOME", "/"), "PYTHONHASHSEED": "0"}
    with open(log, "w") as f:
        p = subprocess.run(cmd, cwd=tree, env=env, stdout=f, stderr=subprocess.STDOUT)
    text = Path(log).read_text(errors="replace")
    return p.returncode == 0, "\nRan 1 test in " in text


def main(argv):
    if argv == ["--list"]:
        for n, rule, path, _, _, test in FAULTS:
            print(f"{n}\t{rule}\t{path}\t{test}")
        print(f"faults: {len(FAULTS)}")
        return 0
    if not argv:
        print(__doc__.strip().split("\n\n")[1], file=sys.stderr)
        return 2
    chosen = argv[1:]
    unknown = sorted(set(chosen) - {f[0] for f in FAULTS})
    if unknown:
        print(f"no fault {', '.join(unknown)}", file=sys.stderr)
        return 2
    try:
        design_mod.guard(argv[0], "out dir")
    except design_mod.Stop as e:
        print(f"stagger_faults: stopped: {e}", file=sys.stderr)
        return 1
    out = Path(argv[0]).resolve()
    if out.exists() and any(out.iterdir()):
        print("the out dir holds a file; nothing is written", file=sys.stderr)
        return 2
    (out / "logs").mkdir(parents=True, exist_ok=True)
    tree = out / "tree"
    tree.mkdir()
    archive = subprocess.run(["git", "-C", str(REPO), "archive", "HEAD"], check=True, capture_output=True).stdout
    subprocess.run(["tar", "-x", "-C", str(tree)], input=archive, check=True)
    rows = []
    for n, rule, path, old, new, test in FAULTS:
        if chosen and n not in chosen:
            continue
        target = tree / path
        text = target.read_text(encoding="utf-8")
        if text.count(old) != 1:
            print(f"fault {n}: the text to change occurs {text.count(old)} times in {path}", file=sys.stderr)
            return 1
        passed, ran = run_test(tree, test, out / "logs" / f"{n}-clean.log")
        if not (passed and ran):
            result = "its test fails with no fault" if ran else "not one test ran"
        else:
            target.write_text(text.replace(old, new), encoding="utf-8")
            try:
                passed, ran = run_test(tree, test, out / "logs" / f"{n}.log")
            finally:
                target.write_text(text, encoding="utf-8")
            result = "caught" if ran and not passed else ("not one test ran" if not ran else "not caught")
        rows.append((n, rule, path, test, result))
        print(f"{n}: {result} by {test}", flush=True)
    with open(out / "faults.tsv", "w", encoding="utf-8") as f:
        f.write("fault\trule\tfile\ttest\tresult\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    caught = sum(r[4] == "caught" for r in rows)
    print(f"caught {caught} of {len(rows)}")
    return 0 if caught == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
