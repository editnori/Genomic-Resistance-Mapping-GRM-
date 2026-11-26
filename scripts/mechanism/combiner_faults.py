"""combiner_faults.py: the faults of the combiner reading's script (combiner.py, C1 to C91), each planted in turn, and
the SYNTHETIC test that must catch it. C39 to C60 are the rules of ruling R30, C61 to C74 those of ruling R37, and
C75 to C91 those of ruling R43, the last four the third review's N11, N12, N15 and N16.

    python3 -I scripts/mechanism/combiner_faults.py [--out faults.tsv]

Run from the repository root. For each fault the named test runs once as the script stands, where it must pass,
and once with the fault planted, where it must fail. Each fault changes one piece of the script's logic. Exit 1
unless every fault is caught. The faults are planted in memory only.

Python 3.12, standard library only; run with -I.
"""

import argparse
import fractions
import hashlib
import json
import math
import os
import re
import sys
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for p in (HERE, ROOT):
    if p not in sys.path:
        sys.path.insert(0, p)

from tests.mechanism import common as fixtures  # noqa: E402

planted = fixtures.load("planted_faults")
lo = fixtures.load("learner_overlap")
common = fixtures.load("common")
cb = fixtures.load("combiner")

T = "tests.mechanism.test_combiner."


def every_genome(folds, f):
    # fault C1: each fold's arm fitted on every genome, its own held-out genomes among them
    return list(range(len(folds)))


def majority(scores):
    # fault C5: three or more of the five calls in place of their mean
    return sum(1 for x in scores if x >= cb.CUTOFF) >= 3


def balanced_loss_grad(w, groups):
    # fault C6: the weighted mean fitted on a log loss with each class weighted to one half
    totals = {True: sum(k for _, t, k in groups if t), False: sum(k for _, t, k in groups if not t)}
    total, grad = 0.0, [0.0] * len(w)
    for s, t, k in groups:
        p = math.fsum(a * b for a, b in zip(w, s))
        c = cb.clipped(p)
        share = k / (2 * totals[t])
        total += share * (-math.log(c) if t else -math.log(1.0 - c))
        if cb.CLIP < p < 1.0 - cb.CLIP:
            d = share * (-1.0 / p if t else 1.0 / (1.0 - p))
            for j in range(len(w)):
                grad[j] += d * s[j]
    return total, grad


def clipped_only(v):
    # fault C7: the weights kept non-negative but not brought back to a sum of one
    return [max(x, 0.0) for x in v]


def plain_weighted(w, s):
    # fault C9: the weighted mean as a plain sum, free to round across the bound of the five scores
    return math.fsum(a * b for a, b in zip(w, s))


def always_converged(groups, _orig=cb.fit_logistic):
    # fault C10: a logistic fold ended with no finite step counted as converged
    beta, steps, end = _orig(groups)
    return beta, steps, "converged" if end == "no_finite_step" else end


def weighted_cap_converged(groups, _orig=cb.fit_weighted):
    # fault C34: a weighted fold ended at its cap counted as converged
    w, steps, end = _orig(groups)
    return w, steps, "converged" if end == "cap" else end


def second_count_decides(rows, arm):
    # fault C36: the count without the marked reports decides, and the count over all stands beside it
    return [("all", cb.unmarked(rows, arm), True), ("without_marked", rows, False),
            ("without_flagged", cb.unflagged(rows, arm), True)]


def vme_of_susceptible(c, name, _orig=cb.measure):
    # fault C12: the very major error rate worked as fp / (fp + tn)
    return 1.0 - _orig(c, "specificity") if name == "vme" else _orig(c, name)


def coverage_as_accuracy(c, genomes):
    # fault C13: the share of genomes called right written in place of the share given a call
    return fractions.Fraction(c["tp"] + c["tn"], genomes)


def accuracy_as_coverage(c, name, _orig=cb.measure):
    # fault C14: the share of genomes given a call written in place of the balanced accuracy
    if name == "balanced_accuracy":
        return float(fractions.Fraction(int(c["called"]), int(c["tp"]) + int(c["fn"]) + int(c["fp"]) + int(c["tn"])))
    return _orig(c, name)


def report_interval(by_species, name):
    # fault C16: the interval from resamples of reports in place of species
    flat = [v for k in sorted(by_species) for v in by_species[k]]
    return common.species_interval({f"report {i}": [v] for i, v in enumerate(flat)}, name)


def one_seed(by_species, name):
    # fault C17: every row's interval drawn from one seed, not the row's own name
    return common.species_interval(by_species, "combiner")


def k1_with_all5(labels, mean_calls, arm_calls, called_by, _orig=cb.transitions):
    # fault C22: a genome one of the five calls counted with those none calls
    return _orig(labels, mean_calls, arm_calls, [0 if k <= 1 else k for k in called_by])


def per_genome(labels, S, folds, _orig=cb.report_counts):
    # fault C24: each genome's call kept in the reading
    row, gate = _orig(labels, S, folds)
    return dict(row, calls=[int(t) for t in labels]), gate


def no_aside(groups, _orig=cb.column_sets):
    # fault C39: a score column with one value over the fitting genomes fitted as a column of its own
    aside, sets = _orig(groups)
    return [], sorted(sets + [[j] for j in aside])


def no_sets(groups, _orig=cb.column_sets):
    # fault C40: score columns the same on every fitting genome fitted apart
    aside, sets = _orig(groups)
    return aside, [[j] for st in sets for j in st]


def first_takes_all(coef, sets):
    # fault C41: the coefficient of a set of columns the same given whole to its first column
    beta = [coef[0]] + [0.0] * len(cb.FIVE)
    for c, st in zip(coef[1:], sets):
        beta[1 + st[0]] = c
    return beta


def near_sets(groups):
    # fault C42: "the same" and "one value" read with a tolerance, not with ==
    def close(a, b):
        return math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-12)

    cols = range(len(cb.FIVE))
    aside = [j for j in cols if all(close(s[j], groups[0][0][j]) for s, _, _ in groups)]
    sets = []
    for j in cols:
        if j in aside or any(j in st for st in sets):
            continue
        sets.append([j] + [m for m in cols if m > j and m not in aside and all(close(s[m], s[j]) for s, _, _ in groups)])
    return aside, sets


def ridged(H, g, _orig=cb.solve):
    # fault C43: 1e-12 added to the diagonal of every Newton system, so no fit is the one before R30
    return _orig([[v + (1e-12 if i == j else 0.0) for j, v in enumerate(row)] for i, row in enumerate(H)], g)


def marked_by_columns(recs):
    # fault C47: a fold with a score column set aside or shared marks its report
    return any(r["end"] != "converged" or cb.aside_or_shared(r) for r in recs)


def old_fold_ok(rec, genomes):
    # fault C49: gate 2 reads the fold's training list only, not what its fit was handed
    return rec["shared"] == 0 and rec["listed"] + rec["tested"] == genomes


def third_count_beside(rows, arm):
    # fault C51: the count without the flagged reports stands beside the count over all and decides nothing
    return [("all", rows, True), ("without_marked", cb.unmarked(rows, arm), False),
            ("without_flagged", cb.unflagged(rows, arm), False)]


def whole_report_flag(S, folds):
    # fault C52: a report flagged only when a learner gives every genome of the report one score
    return any(all(s[j] == S[0][j] for s in S) for j in range(len(cb.FIVE)))


def old_rules(reading):
    # fault C55: write checks the cutoff, forms, five learners and clip only, not the counts of the rule or the flags
    return ((reading.get("cutoff"), reading.get("forms"), reading.get("five"), reading.get("clip"))
            == (cb.CUTOFF, list(cb.FORMS), list(cb.FIVE), cb.CLIP))


def under_0_or_at_0(beta):
    # fault C62: a coefficient of exactly 0, as a column set aside has, counted as under 0
    return sum(1 for b in beta[1:] if b <= 0.0)


def under_0_with_tolerance(beta):
    # fault C63: a coefficient counted as under 0 only below -1e-9
    return sum(1 for b in beta[1:] if b < -1e-9)


def intercept_counted(beta):
    # fault C64: the intercept counted among the five coefficients
    return sum(1 for b in beta if b < 0.0)


def first_fold_only(recs):
    # fault C65: the second condition read on the first fold fit only
    return cb.under_0(recs[0])


def coefficients_not_folds(recs):
    # fault C69: the record counts the coefficients under 0, not the fold fits that give one
    return sum(r.get("coefficients_under_0", 0) for r in recs)


def either_as_both(rows):
    # fault C70: a report with either condition counted as one with both
    return (sum(1 for r in rows if r["one_score_fold"]), sum(1 for r in rows if r["coefficient_under_0"]),
            sum(1 for r in rows if r["one_score_fold"] or r["coefficient_under_0"]))


def later_steps_named(gate):
    # fault C73: gate 5 also names a report whose logistic fold fit met no finite solution at a later step
    return gate["no_first_step"] > 0 or gate["ends"]["logistic"]["no_finite_step"] > 0


def outcome_from_the_rule(r):
    # fault C75: the logistic arm's outcome in each count taken from the rule, as before R43
    return dict(r, decides="no", reading=cb.DESCRIPTIVE)


def reading_from_the_counts(r):
    # fault C76: the logistic arm's reading taken from its counts that decide, as before R43
    return dict(r, decides="no", outcome=cb.DESCRIPTIVE)


def figures_dropped(r):
    # fault C77: the logistic arm's reports, species, medians and intervals not written
    return dict(r, decides="no", outcome=cb.DESCRIPTIVE, reading=cb.DESCRIPTIVE, datasets=0, species=0,
                median_gain=None, gain_lo=None, gain_hi=None, vme_median=None, vme_lo=None, vme_hi=None,
                gain_at_least_margin=0)


def kept_by_the_environment(arm):
    # fault C79: the arms the rule can keep read from an environment setting, the weighted mean when it is not set
    return arm in os.environ.get("COMBINER_CAN_BE_KEPT", "weighted").split(",")


def parser_with_arms(_orig=cb.parser):
    # fault C80: read and write take an option that names the arms the rule can keep
    ap = _orig()
    commands = next(a for a in ap._actions if isinstance(a, argparse._SubParsersAction)).choices
    for name in ("read", "write"):
        commands[name].add_argument("--can-be-kept", default=",".join(cb.CAN_BE_KEPT))
    return ap


def write_unchecked(*args, _orig=cb.write_record):
    # fault C81: write blind to a rule row that gives the logistic arm another outcome or reading
    with mock.patch.object(cb, "not_descriptive", lambda rules: []):
        return _orig(*args)


def tables_unchecked(record, _orig=cb.tables):
    # fault C82: print blind to a rule.tsv that gives the logistic arm another outcome or reading
    with mock.patch.object(cb, "not_descriptive", lambda rules: []):
        return _orig(record)


def mark_only(d, _orig=cb.worked_out):
    # fault C84: write works the logistic mark out again, and not the sign flag
    return [what for what in _orig(d) if what != "its sign flag"]


def sign_only(d, _orig=cb.worked_out):
    # fault C85: write works the sign flag out again, and not the logistic mark
    return [what for what in _orig(d) if what != "its logistic mark"]


def species_only(d, what):
    # fault C86: the stop line names the report by species only
    verb = "differ" if len(what) > 1 else "differs"
    return f"stop flags: {d.get('species')}: {' and '.join(what)} {verb} from its logistic fold records: STOP"


def first_report_only(reading, _orig=cb.flags_differ):
    # fault C87: write names only the first report whose flags differ
    return _orig(reading)[:1]


def sign_count_twice(rows, _orig=cb.flag_counts):
    # fault C88, the review's N11: the flags line gives the sign's count where R30's belongs
    one_score, under, both = _orig(rows)
    return under, under, both


def text_flag_taken(reading):
    # fault C90, the review's N15: write takes a report's sign flag given as text
    return (reading.get("flags") == {a: list(cb.FLAGS[a]) for a in cb.ARMS}
            and all(isinstance(d.get("coefficient_under_0"), (bool, str))
                    and all(type(r.get("coefficients_under_0")) is int for r in d.get("logistic_folds") or [])
                    for d in reading.get("datasets") or []))


def bool_count_taken(reading):
    # fault C91, the review's N16: write takes a fold fit's count given as true or false
    return (reading.get("flags") == {a: list(cb.FLAGS[a]) for a in cb.ARMS}
            and all(isinstance(d.get("coefficient_under_0"), bool)
                    and all(isinstance(r.get("coefficients_under_0"), int) for r in d.get("logistic_folds") or [])
                    for d in reading.get("datasets") or []))


def plain_load(path):
    # fault C26: a report parsed with every genome id kept
    with open(path, encoding="utf-8") as f:
        return json.load(f)


FAULTS = [
    ("C1", "each fold's arm fitted on every genome, its own held-out genomes among them",
     T + "Crossfit.test_gate_2_counts_no_shared_genome_on_the_standard_read",
     lambda: mock.patch.object(cb, "training_genomes", every_genome)),
    ("C2", "gate 2 blind to a fit on its own held-out genomes", T + "Gates.test_gate_2_stops_a_fit_on_its_own_genomes",
     lambda: mock.patch.object(cb, "fold_ok", lambda rec, genomes: True)),
    ("C3", "a score of exactly 0.5 read as a miss", T + "Forms.test_a_score_of_exactly_the_cutoff_is_a_call",
     lambda: mock.patch.object(cb, "is_call", lambda score: score > cb.CUTOFF)),
    ("C4", "the arms called at the training prevalence (0.425 in the test's design) in place of 0.5",
     T + "Fits.test_the_logistic_fits_each_groups_share",
     lambda: mock.patch.object(cb, "CUTOFF", 0.425)),
    ("C5", "the majority of the five read in place of their mean",
     T + "Forms.test_the_incumbent_is_the_mean_not_the_majority",
     lambda: mock.patch.object(cb, "mean_call", majority)),
    ("C6", "the weighted mean fitted on a class-weighted log loss", T + "Fits.test_the_weighted_mean_fits_plain_log_loss",
     lambda: mock.patch.object(cb, "weighted_loss_grad", balanced_loss_grad)),
    ("C7", "the weights not brought back to a sum of one", T + "Fits.test_the_weights_stay_on_the_simplex",
     lambda: mock.patch.object(cb, "project_simplex", clipped_only)),
    ("C8", "gate 4 blind to weights that do not sum to one", T + "Gates.test_gate_4_stops_on_weights_off_the_simplex",
     lambda: mock.patch.object(cb, "SUM_TOL", math.inf)),
    ("C9", "the weighted mean free to round across the bound of the five scores",
     T + "Forms.test_the_weighted_mean_stays_inside_the_five_scores",
     lambda: mock.patch.object(cb, "weighted_score", plain_weighted)),
    ("C10", "a logistic fold ended with no finite step counted as converged",
     T + "Gates.test_gate_5_counts_a_fold_with_no_finite_step",
     lambda: mock.patch.object(cb, "fit_logistic", always_converged)),
    ("C11", "log loss not clipped", T + "Forms.test_log_loss_is_clipped",
     lambda: mock.patch.object(cb, "clipped", lambda p: p)),
    ("C12", "the very major error rate worked over the susceptible genomes",
     T + "Forms.test_measures_are_worked_from_the_counts",
     lambda: mock.patch.object(cb, "measure", vme_of_susceptible)),
    ("C13", "accuracy written in place of coverage", T + "Forms.test_coverage_is_the_share_given_a_call",
     lambda: mock.patch.object(cb, "coverage", coverage_as_accuracy)),
    ("C14", "coverage written in place of the balanced accuracy", T + "Forms.test_balanced_accuracy_is_not_coverage",
     lambda: mock.patch.object(cb, "measure", accuracy_as_coverage)),
    ("C15", "gate 6 blind to a genome given no call", T + "Gates.test_gate_6_stops_on_a_genome_given_no_call",
     lambda: mock.patch.object(cb, "coverage", lambda c, genomes: fractions.Fraction(1))),
    ("C16", "the interval from resamples of reports, not species", T + "Pooled.test_the_interval_resamples_species",
     lambda: mock.patch.object(cb, "interval", report_interval)),
    ("C17", "every row's interval drawn from one seed", T + "Pooled.test_the_interval_is_seeded_by_the_rows_name",
     lambda: mock.patch.object(cb, "interval", one_seed)),
    ("C18", "an arm kept with a gain interval that does not lie above 0",
     T + "Rule.test_the_gain_interval_must_lie_above_0",
     lambda: mock.patch.object(cb, "keep", lambda m, g, v: m >= cb.MARGIN and not (v is not None and v > 0))),
    ("C19", "a median gain of exactly 0.02 read as under the margin", T + "Rule.test_the_gain_must_reach_the_margin",
     lambda: mock.patch.object(cb, "keep", lambda m, g, v: m > cb.MARGIN and g is not None and g > 0
                               and not (v is not None and v > 0))),
    ("C20", "an arm kept with a higher very major error rate", T + "Rule.test_a_higher_vme_blocks_the_arm",
     lambda: mock.patch.object(cb, "keep", lambda m, g, v: m >= cb.MARGIN and g is not None and g > 0)),
    ("C21", "a higher very major error rate counted as better", T + "Pooled.test_better_worse_and_same_at_the_margin",
     lambda: mock.patch.object(cb, "LOWER_BETTER", ("log_loss",))),
    ("C22", "a genome one of the five calls counted with those none calls",
     T + "Counts.test_the_transitions_split_by_how_many_of_the_five_call",
     lambda: mock.patch.object(cb, "transitions", k1_with_all5)),
    ("C23", "arm counts that do not follow from the transitions let through",
     T + "Gates.test_gate_7_stops_on_arm_counts_that_do_not_add_up",
     lambda: mock.patch.object(cb, "arithmetic_holds", lambda r: True)),
    ("C24", "each genome's call kept in the reading", T + "Ids.test_the_reading_holds_no_row_per_genome",
     lambda: mock.patch.object(cb, "report_counts", per_genome)),
    ("C25", "gate 7 blind to a list as long as a report", T + "Gates.test_gate_7_stops_on_a_row_per_genome",
     lambda: mock.patch.object(cb, "MAX_LIST", 10 ** 9)),
    ("C26", "a genome id kept as a report is parsed", T + "Ids.test_no_id_is_held_after_parsing",
     lambda: mock.patch.object(cb, "load_report", plain_load)),
    ("C27", "gate 8 blind to a genome id with no id shape", T + "Gates.test_gate_8_stops_on_a_genome_id_in_the_reading",
     lambda: mock.patch.object(cb, "ids_of", lambda directory: set())),
    ("C28", "gate 8 blind to an id-shaped token", T + "Gates.test_gate_8_stops_on_an_id_shaped_token",
     lambda: mock.patch.object(lo, "ID_SHAPE", re.compile(r"(?!x)x"))),
    ("C29", "a decimal cell of 10 or more let through", T + "Gates.test_a_decimal_cell_stays_under_10",
     lambda: mock.patch.object(cb, "DECIMAL_LIMIT", 10 ** 9)),
    ("C30", "a report set that differs from the first record's let through",
     T + "Gates.test_gate_1_stops_on_a_report_set_that_differs",
     lambda: mock.patch.object(cb, "same_reports", lambda wanted, read_keys: True)),
    ("C31", "an incumbent that differs from vote.tsv let through",
     T + "Gates.test_gate_3_stops_on_an_incumbent_that_differs_from_vote_tsv",
     lambda: mock.patch.object(cb, "incumbent_agrees", lambda row, vote: True)),
    ("C32", "a held-out score that is not finite let through",
     T + "Gates.test_gate_5_stops_on_a_held_out_score_that_is_not_finite",
     lambda: mock.patch.object(cb, "finite", lambda *values: True)),
    ("C33", "a fold ended at its cap left unmarked", T + "Gates.test_gate_5_counts_and_marks_a_fold_at_the_cap",
     lambda: mock.patch.object(cb, "marked", lambda recs: False)),
    ("C34", "a weighted fold ended at its cap counted as converged", T + "Fits.test_a_fit_at_its_cap_ends_there",
     lambda: mock.patch.object(cb, "fit_weighted", weighted_cap_converged)),
    ("C35", "the marked reports kept in the second count",
     T + "Rule.test_the_second_count_leaves_out_the_marked_reports",
     lambda: mock.patch.object(cb, "unmarked", lambda rows, arm: rows)),
    ("C36", "the count without the marked reports deciding in place of the count over all",
     T + "Rule.test_the_second_count_leaves_out_the_marked_reports",
     lambda: mock.patch.object(cb, "rule_counts", second_count_decides)),
    ("C37", "the weighted fit started at the first learner's vertex in place of the equal weights",
     T + "Fits.test_the_weighted_fit_is_never_above_the_equal_mean",
     lambda: mock.patch.object(cb, "start_weights", lambda: [1.0, 0.0, 0.0, 0.0, 0.0])),
    ("C38", "the record's rule for the very major error rate in other words than R19's",
     T + "Rule.test_the_record_says_the_vme_rule_in_r19s_words",
     lambda: mock.patch.object(cb, "VME_RULE", "the very major error rate is not higher")),
    ("C39", "a score column with one value over the fitting genomes fitted as a column of its own",
     T + "Fits.test_a_one_score_column_is_set_aside",
     lambda: mock.patch.object(cb, "column_sets", no_aside)),
    ("C40", "score columns the same on every fitting genome fitted apart",
     T + "Fits.test_columns_the_same_are_fitted_as_one",
     lambda: mock.patch.object(cb, "column_sets", no_sets)),
    ("C41", "a shared coefficient given whole to the first column of its set",
     T + "Fits.test_a_held_out_genome_is_scored_from_both_columns",
     lambda: mock.patch.object(cb, "spread", first_takes_all)),
    ("C42", "the same and one value read with a tolerance, not exactly",
     T + "Fits.test_the_same_and_one_score_are_read_without_tolerance",
     lambda: mock.patch.object(cb, "column_sets", near_sets)),
    ("C43", "1e-12 added to the diagonal of every Newton system",
     T + "Fits.test_a_fit_on_neither_case_is_the_fit_before_r30",
     lambda: mock.patch.object(cb, "solve", ridged)),
    ("C44", "a first Newton step with no finite solution ended as a later one is",
     T + "Gates.test_gate_5_stops_on_a_first_step_with_no_finite_solution",
     lambda: mock.patch.object(cb, "failed_step", lambda step: "no_finite_step")),
    ("C45", "a fold whose first Newton step had no finite solution scored from its coefficients at 0",
     T + "Fits.test_a_first_step_with_no_finite_solution_is_not_scored",
     lambda: mock.patch.object(cb, "scored", lambda end: True)),
    ("C46", "gate 5 blind to a first Newton step with no finite solution",
     T + "Gates.test_gate_5_stops_on_a_first_step_with_no_finite_solution",
     lambda: mock.patch.object(cb, "ends_pass", lambda not_finite, no_first_step, fits: not_finite == 0 and fits > 0)),
    ("C47", "a fold with a score column set aside or shared marks its report",
     T + "Gates.test_the_reviewers_two_reports_pass_gate_5",
     lambda: mock.patch.object(cb, "marked", marked_by_columns)),
    ("C48", "a fold with a score column shared not counted with those set aside",
     T + "Gates.test_the_reviewers_two_reports_pass_gate_5",
     lambda: mock.patch.object(cb, "aside_or_shared", lambda rec: rec.get("columns_aside", 0) > 0)),
    ("C49", "gate 2 blind to a fit handed another count of genomes than its training list",
     T + "Gates.test_gate_2_stops_a_fit_handed_every_genome",
     lambda: mock.patch.object(cb, "fold_ok", old_fold_ok)),
    ("C50", "the training count taken as the groups a fit is handed, not the genomes inside them",
     T + "Crossfit.test_each_fold_is_fitted_on_the_other_folds",
     lambda: mock.patch.object(cb, "handed", len)),
    ("C51", "the count without the flagged reports deciding nothing",
     T + "Rule.test_the_third_count_decides",
     lambda: mock.patch.object(cb, "rule_counts", third_count_beside)),
    ("C52", "a report flagged only when a learner gives the whole report one score",
     T + "Counts.test_a_report_is_flagged_by_one_fold",
     lambda: mock.patch.object(cb, "one_score_fold", whole_report_flag)),
    ("C53", "a fold of one genome flagging its report",
     T + "Counts.test_a_report_is_flagged_by_one_fold",
     lambda: mock.patch.object(cb, "FLAG_FOLD", 1)),
    ("C54", "the flagged reports kept in the third count",
     T + "Rule.test_the_third_count_decides",
     lambda: mock.patch.object(cb, "unflagged", lambda rows, arm: rows)),
    ("C55", "write blind to a reading made under other counts of the rule or without its flags",
     T + "Gates.test_write_stops_on_a_reading_made_under_other_rules",
     lambda: mock.patch.object(cb, "same_rules", old_rules)),
    ("C56", "an arm kept when either count that decides keeps it",
     T + "Rule.test_the_third_count_decides",
     lambda: mock.patch.object(cb, "arm_reading", lambda outcomes: "kept" if "kept" in outcomes
                               else "no difference shown")),
    ("C57", "gate 1 blind to the count of reports it expects",
     T + "Gates.test_gate_1_holds_the_read_to_its_count",
     lambda: mock.patch.object(cb, "held_to", lambda read_count, expected: True)),
    ("C58", "the command line holding the read to 267 reports, not the freeze's 268",
     T + "Gates.test_the_command_line_holds_the_read_to_268",
     lambda: mock.patch.object(cb, "FROZEN_REPORTS", 267)),
    ("C59", "a genome with no fold number let through",
     T + "Gates.test_a_genome_with_no_fold_number_stops_the_read",
     lambda: mock.patch.object(cb, "fold_number", lambda g: True)),
    ("C60", "a report with one fold let through",
     T + "Gates.test_a_report_with_one_fold_stops_the_read",
     lambda: mock.patch.object(cb, "enough_folds", lambda folds: True)),
    ("C61", "no report flagged for a logistic fold fit that gives a learner a coefficient under 0",
     T + "Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only",
     lambda: mock.patch.object(cb, "coefficient_under_0", lambda recs: False)),
    ("C62", "a coefficient of exactly 0, as a column set aside has, counted as under 0",
     T + "Flags.test_a_column_set_aside_at_0_does_not_flag",
     lambda: mock.patch.object(cb, "negatives", under_0_or_at_0)),
    ("C63", "a coefficient counted as under 0 only below -1e-9",
     T + "Flags.test_a_coefficient_under_0_is_read_with_no_tolerance",
     lambda: mock.patch.object(cb, "negatives", under_0_with_tolerance)),
    ("C64", "the intercept counted among the five coefficients",
     T + "Flags.test_a_report_with_every_coefficient_at_least_0_is_not_flagged",
     lambda: mock.patch.object(cb, "negatives", intercept_counted)),
    ("C65", "the second condition read on the first fold fit only",
     T + "Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only",
     lambda: mock.patch.object(cb, "coefficient_under_0", first_fold_only)),
    ("C66", "the logistic arm's third count leaving out R30's flagged reports only",
     T + "Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only",
     lambda: mock.patch.object(cb, "FLAGS", {"weighted": ("one_score_fold",), "logistic": ("one_score_fold",)})),
    ("C67", "the weighted arm's third count also leaving out the reports flagged for a coefficient under 0",
     T + "Flags.test_one_fold_fit_under_0_flags_the_logistic_arm_only",
     lambda: mock.patch.object(cb, "FLAGS", {a: ("one_score_fold", "coefficient_under_0") for a in cb.ARMS})),
    ("C68", "a count of 5 species given an interval",
     T + "Flags.test_a_logistic_third_count_under_6_species_gets_no_interval",
     lambda: mock.patch.object(common, "MIN_SPECIES", 5)),
    ("C69", "the record counting the coefficients under 0, not the fold fits behind the flag",
     T + "Flags.test_the_record_counts_the_fold_fits_behind_the_flag",
     lambda: mock.patch.object(cb, "under_0_folds", coefficients_not_folds)),
    ("C70", "a report with either condition counted as one with both",
     T + "Flags.test_the_read_prints_how_many_reports_each_condition_flags",
     lambda: mock.patch.object(cb, "flag_counts", either_as_both)),
    ("C71", "write blind to a reading made without the second condition",
     T + "Flags.test_write_stops_on_a_reading_made_without_the_second_condition",
     lambda: mock.patch.object(cb, "second_condition_written", lambda reading: True)),
    ("C72", "gate 5 naming no report whose first Newton step has no finite solution",
     T + "Flags.test_gate_5_names_each_report_with_no_finite_first_step",
     lambda: mock.patch.object(cb, "first_step_named", lambda gate: False)),
    ("C73", "gate 5 also naming a report whose later Newton step has no finite solution",
     T + "Flags.test_gate_5_names_each_report_with_no_finite_first_step",
     lambda: mock.patch.object(cb, "first_step_named", later_steps_named)),
    ("C74", "the reading's sha256 worked without its last line break",
     T + "Flags.test_read_and_write_print_the_readings_sha256",
     lambda: mock.patch.object(cb, "reading_sha256", lambda raw: hashlib.sha256(raw.rstrip(b"\n")).hexdigest())),
    ("C75", "the logistic arm's outcome in each count taken from the rule",
     T + "Descriptive.test_the_logistic_arm_reads_descriptive_only_whatever_its_figures",
     lambda: mock.patch.object(cb, "descriptive", outcome_from_the_rule)),
    ("C76", "the logistic arm's reading taken from its counts that decide",
     T + "Descriptive.test_the_logistic_arm_reads_descriptive_only_whatever_its_figures",
     lambda: mock.patch.object(cb, "descriptive", reading_from_the_counts)),
    ("C77", "the logistic arm's reports, species, medians and intervals not written",
     T + "Descriptive.test_the_logistic_arm_is_descriptive_only_in_each_count",
     lambda: mock.patch.object(cb, "descriptive", figures_dropped)),
    ("C78", "both arms named as arms the rule can keep",
     T + "Descriptive.test_one_line_names_the_arms_that_can_be_kept",
     lambda: mock.patch.object(cb, "CAN_BE_KEPT", cb.ARMS)),
    ("C79", "the arms the rule can keep read from an environment setting",
     T + "Descriptive.test_no_option_or_environment_setting_changes_the_arms_that_can_be_kept",
     lambda: mock.patch.object(cb, "can_be_kept", kept_by_the_environment)),
    ("C80", "an option of read and write naming the arms the rule can keep",
     T + "Descriptive.test_no_option_or_environment_setting_changes_the_arms_that_can_be_kept",
     lambda: mock.patch.object(cb, "parser", parser_with_arms)),
    ("C81", "write blind to a rule row that gives the logistic arm another outcome or reading",
     T + "Descriptive.test_write_stops_on_a_rule_row_that_gives_the_logistic_arm_another_outcome",
     lambda: mock.patch.object(cb, "write_record", write_unchecked)),
    ("C82", "print blind to a rule.tsv that gives the logistic arm another outcome or reading",
     T + "Descriptive.test_print_stops_on_a_record_that_gives_the_logistic_arm_another_outcome",
     lambda: mock.patch.object(cb, "tables", tables_unchecked)),
    ("C83", "write blind to a reading made under other arms that can be kept",
     T + "Descriptive.test_write_stops_on_a_reading_made_under_other_arms",
     lambda: mock.patch.object(cb, "kept_arms_written", lambda reading: True)),
    ("C84", "write not working the sign flag out again from the fold records",
     T + "Flags.test_write_works_the_sign_flag_out_again",
     lambda: mock.patch.object(cb, "worked_out", mark_only)),
    ("C85", "write not working the logistic mark out again from the fold records",
     T + "Flags.test_write_works_the_logistic_mark_out_again",
     lambda: mock.patch.object(cb, "worked_out", sign_only)),
    ("C86", "write's stop line naming the report by species only",
     T + "Flags.test_write_names_each_report_whose_flags_differ",
     lambda: mock.patch.object(cb, "flags_stop", species_only)),
    ("C87", "write naming only the first report whose flags differ",
     T + "Flags.test_write_names_each_report_whose_flags_differ",
     lambda: mock.patch.object(cb, "flags_differ", first_report_only)),
    ("C88", "the flags line giving the sign's count where R30's belongs",
     T + "Flags.test_the_flags_line_gives_each_conditions_own_count",
     lambda: mock.patch.object(cb, "flag_counts", sign_count_twice)),
    ("C89", "gate 5 naming a report only with 2 or more folds whose first Newton step has no finite solution",
     T + "Flags.test_gate_5_names_a_report_with_exactly_1_fold_with_no_finite_first_step",
     lambda: mock.patch.object(cb, "first_step_named", lambda gate: gate["no_first_step"] > 1)),
    ("C90", "write taking a report's sign flag given as text",
     T + "Flags.test_write_stops_on_a_sign_flag_given_as_text",
     lambda: mock.patch.object(cb, "second_condition_written", text_flag_taken)),
    ("C91", "write taking a fold fit's count given as true or false",
     T + "Flags.test_write_stops_on_a_fold_count_given_as_true_or_false",
     lambda: mock.patch.object(cb, "second_condition_written", bool_count_taken)),
]


def run_all():
    rows = []
    for n, name, test_id, plant in FAULTS:
        clean, with_fault = planted.outcome(test_id), planted.outcome(test_id, plant)
        rows.append({"fault": n, "planted": name, "test": test_id.rsplit(".", 2)[-2] + "." + test_id.rsplit(".", 1)[-1],
                     "as_written": clean, "with_fault": with_fault,
                     "caught": "yes" if clean == "pass" and with_fault == "fail" else "no"})
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rows = run_all()
    header = ["fault", "planted", "test", "as_written", "with_fault", "caught"]
    lines = ["\t".join(header)] + ["\t".join(str(r[h]) for h in header) for r in rows]
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    caught = sum(r["caught"] == "yes" for r in rows)
    print(f"planted faults: {caught} of {len(rows)} caught")
    return 0 if caught == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
