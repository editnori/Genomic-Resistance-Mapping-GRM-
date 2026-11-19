"""Every line of docs/study-split.md that carries a number, and the number-free lines it holds, rebuilt from reference/study-split.

groups(rec) gives the page's built blocks in page order, as (name, lines):
a paragraph is one line, a table or a list several. lines(rec) is every
line of them. The page test holds each line of the page that carries a
number to lines(rec), and each line of lines(rec) that carries a number to
the page, once and in this order. The words around the figures are typed
here; every figure comes from a cell of the record.

held(rec) gives the number-free lines the page test also holds word for
word, once: the opening's denial and lineage sentence, the tuning sentence,
the plans' exploratory paragraph, the division of the drop, the reading of
the part not matched by fold shape, and every number-free item of the Limits
list. limits(rec) gives the Limits list itself, whole and in order.
Rewording one of these lines means rewording it here too.
"""

import re
from fractions import Fraction

from . import common

LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
GRADED = LEARNERS[1:]
FAMILIES = ("SR", "SM", "MR")
FAMILY = {"SR": "study-random", "SM": "study-matched", "MR": "matched-random"}
DRAWS = range(5)
NUMBER = re.compile(r"(?<![A-Za-z0-9_])\d")
WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve",
         "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen", "twenty")
TENS = ("thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety")
NUMBER_WORD = re.compile(r"(?i)\b(" + "|".join(WORDS + TENS) + r")\b")


def carries_number(line):
    """A digit that starts a number, or a number word from zero to twenty or a ten."""
    return NUMBER.search(line) is not None or NUMBER_WORD.search(line) is not None


def word(x):
    """A count from zero to twenty as its word: 5 -> five."""
    return WORDS[int(x)]


def num(x):
    return Fraction(x) if "/" in str(x) else Fraction(float(x)) if "." in str(x) or "e" in str(x) else Fraction(int(x))


def signed(x):
    """A difference with its sign and four decimals: +0.2841, -0.0181."""
    return f"{float(x):+.4f}"


def bound(x):
    """An interval bound: four decimals, and 0 for an exact zero."""
    x = float(x)
    return "0" if x == 0 else f"{x:.4f}"


def share(x):
    """A share with three decimals, and 0 for an exact zero."""
    x = float(x)
    return "0" if x == 0 else f"{x:.3f}"


def plain(x):
    """A count or median with no trailing zeros: 100, 98.5, 95.25."""
    return f"{float(x):g}"


def n(x):
    return f"{int(x):,}"


def pts(x):
    """The size of a difference, four decimals, no sign."""
    return f"{abs(float(x)):.4f}"


def clock(stamp):
    """2026-11-15T17:23:24.89Z -> 17:23:24.89Z; 2026-11-15T11:38:56.00Z -> 11:38:56Z."""
    t = stamp.split("T")[1]
    return t[:-4] + "Z" if t.endswith(".00Z") else t


def joined(items):
    items = list(items)
    if not items:
        return "none"
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def short(slug):
    genus, species = slug.split("_", 1)
    return f"*{genus[0].upper()}. {species.replace('_', ' ')}*"


def by_count(counts, total, unit):
    """[(learner, count)] -> '25 of 25 pairs for rf and gbt, 19 for cart and 6 for scm'."""
    seen = []
    for _, c in sorted(counts, key=lambda lc: -lc[1]):
        if c not in seen:
            seen.append(c)
    parts = []
    for i, c in enumerate(seen):
        who = joined(ln for ln, x in counts if x == c)
        parts.append(f"{c} of {total} {unit} for {who}" if i == 0 else f"{c} for {who}")
    return joined(parts)


class Record:
    def __init__(self, directory=None):
        d = directory or common.RECORD
        self.dir = d
        rd = common.read_tsv
        self.summary = {r["key"]: r["value"] for r in rd(d / "summary.tsv")}
        self.settings = {r["key"]: r["value"] for r in rd(d / "settings.tsv")}
        self.checks = {r["check"]: int(r["count"]) for r in rd(d / "checks.tsv")}
        self.order = {(r["experiment"], r["event"]): r["time"] for r in rd(d / "order.tsv")}
        self.planfiles = {r["experiment"]: r for r in rd(d / "planfiles.tsv")}
        self.contrasts = {(r["experiment"], r["learner"], r["family"], r["instance"], r["measure"]): r
                          for r in rd(d / "contrasts.tsv")}
        self.first = {(r["learner"], r["contrast"], r["measure"]): r for r in rd(d / "first.tsv")}
        self.rules = {(r["learner"], r["family"], r["rule"]): r for r in rd(d / "rules.tsv")}
        self.thresholds = rd(d / "thresholds.tsv")
        self.leaveout = {(r["left_out"], r["learner"], r["family"]): r for r in rd(d / "leaveout.tsv")}
        self.draws = {(r["learner"], r["family"], r["axis"], r["measure"]): [int(r[f"d{i}"]) for i in DRAWS]
                      for r in rd(d / "draws.tsv")}
        self.parts = {(r["learner"], r["measure"]): r for r in rd(d / "parts.tsv")}
        self.mates = {(r["kind"], r["draw"]): r for r in rd(d / "mates.tsv")}
        self.shapes = {(r["arm"], r["measure"]): r for r in rd(d / "shapes.tsv")}
        self.nostudy = {(r["learner"], r["family"], r["measure"], r["subset"]): r for r in rd(d / "nostudy.tsv")}
        self.foldcut = {(r["learner"], r["group"]): r for r in rd(d / "foldcut.tsv")}
        self.datasets = rd(d / "datasets.tsv")
        self.eligibility = rd(d / "eligibility.tsv")
        self.plans = rd(d / "plans.tsv")
        self.files = len((d / "SHA256SUMS").read_text(encoding="utf-8").strip().split("\n")) + 1
        # The page writes the pairs of each family as "of 25" and "of 5": the record's count of them.
        total = {f: {r["pairs"] for (_, fam, _), r in self.rules.items() if fam == f} for f in FAMILIES}
        assert total == {"SR": {"25"}, "SM": {"5"}, "MR": {"25"}}, total

    def s(self, key):
        return self.summary[key]

    def c2(self, l, fam, inst, m):
        return self.contrasts[("second", l, fam, inst, m)]

    def c1(self, l, fam, m):
        return self.contrasts[("first", l, fam, "-", m)]

    def loss(self, l):
        return "error_rate" if l == "scm" else "log_loss"

    def pairs(self, fam):
        return [str(i) for i in DRAWS] if fam == "SM" else [f"{i},{j}" for i in DRAWS for j in DRAWS]

    def worse(self, l, fam, m):
        return sum(self.c2(l, fam, i, m)["reading"] == "worse" for i in self.pairs(fam))

    def species(self):
        out = {}
        for r in self.datasets:
            out[r["species"]] = out.get(r["species"], 0) + 1
        return out


def cell(r):
    return f"{signed(r['median'])} [{bound(r['low'])}, {bound(r['high'])}] {r['reading']}"


def row(*cells):
    return "| " + " | ".join(cells) + " |"


def head(*cells):
    return [row(*cells), "|" + "---|" * len(cells)]


# ----------------------------------------------------------------- opening and plans

def opening(rec):
    s, st = rec.s, rec.settings
    yield "open_release", [
        f"At evaluation revision {st['release_revision']}, {s('datasets')} of the release's "
        f"{s('release_datasets')} datasets, in {s('species')} species, pass the first plan's rule for entering a "
        "split by study (the Folds section gives it). They hold "
        f"{n(s('pairs'))} genome and dataset pairs, {n(s('resistant_pairs'))} of them resistant, in cohorts of "
        f"{n(s('cohort_min'))} to {n(s('cohort_max'))} genomes (median {plain(s('cohort_median'))})."]
    yield "open_study", [
        "Two experiments scored each pair with folds that keep each BioProject whole and with random folds on the "
        "same genomes. A BioProject stands in for a study, and an approximate one: the files cannot show whether "
        f"two BioProjects hold one study, and {n(s('no_study_pairs'))} pairs have no BioProject."]
    yield "open_arms", [
        f"The first experiment compares one assignment of BioProjects to {word(st['folds'])} study folds with the "
        f"release's own random folds and its lineage split. The second draws {word(st['draws'])} study assignments "
        f"and {word(st['draws'])} random draws per dataset and adds a third arm: matched random folds, which copy "
        "each study plan's fold sizes and resistant counts."]
    yield "open_learners", [
        f"{word(len(LEARNERS)).capitalize()} learners fit each fold plan: scm (set covering machine), cart (one "
        "decision tree), rf (random forest), gbt (boosted trees) and lasso. scm scores only 0 or 1, so its log loss "
        f"is {float(s('scm_log_loss_per_error')):.2f} times its error rate plus a constant, and outside the first "
        "experiment's two tables this page reads scm's loss as its error rate."]
    yield "open_denial", [
        "Each experiment is a retrospective split by study inside one database, BV-BRC, and never external "
        "validation. A test genome's study mates are the other genomes of its BioProject in the dataset, and its "
        "lineage-cluster mates those of its cluster in the release's lineage split. The design does not separate "
        "study from lineage: study folds keep every study mate of a test genome out of training, and about half of "
        "its lineage-cluster mates as well."]


def plans(rec):
    o = rec.order
    dates = {v.split("T")[0] for v in o.values()}
    assert len(dates) == 1, dates
    yield "plans_order", [
        "Each experiment hashed its plan before its first fit. `order.tsv` in the record lists every time on this "
        f"page with its source, all on {dates.pop()} (UTC). The Rule section defines the contrasts, the species "
        "interval and a cost."]
    yield "plans_table", head("Plan", "Hashed", "First result", "Fixed") + [
        row("First", clock(o[("first", "hash_time")]), clock(o[("first", "first_fit_end")]),
            "The eligibility rule, one study assignment, recall and log loss for every learner, the species "
            "interval and the rule for a cost"),
        row("Second", clock(o[("second", "hash_time")]), clock(o[("second", "first_fit_end")]),
            "Five draws of each arm, the three contrasts and their thresholds, scm's error rate, and each answer "
            "with one species left out")]
    yield "plans_runs", [
        f"The second experiment's fits started at {clock(o[('second', 'fits_start')])}. The last result came at "
        f"{clock(o[('first', 'last_fit_end')])} in the first experiment and at "
        f"{clock(o[('second', 'last_fit_end')])} in the second."]
    yield "plans_exploratory", [
        "**Exploratory.** Every reading the plans did not fix is exploratory: the first experiment's fold cuts, the "
        "lineage mates, the distances between assignments, the draw split, the two largest species left out "
        "together, the other outcome rules and thresholds, and the shares of the drop. The first plan gave its "
        "lineage contrasts no rule, and the second lists the genomes without a BioProject as a description with no "
        "rule."]


# ----------------------------------------------------------------- folds

def folds(rec):
    s, st, ck, o = rec.s, rec.settings, rec.checks, rec.order
    yield "folds_groups", [
        "The fold builder places groups of genomes. A copy group is a set of genomes the release marks as identical "
        "copies. A group joins a BioProject with every copy group that holds one of its genomes, and so with any "
        "other BioProject such a copy group reaches; a genome with no BioProject is a group of its own unless a "
        "copy group joins it."]
    yield "folds_rule", [
        f"A dataset enters when it holds at least {st['min_studies']} BioProjects, at least "
        f"{st['min_study_share']} of its genomes carry one, and no group holds over {st['max_group_share']} of its "
        f"genomes. Of the release's {s('release_datasets')} datasets, {s('eligible_counts')} pass; "
        f"{s('fail_studies')} fail on BioProjects, {s('fail_share')} on the share and {s('fail_group')} on the "
        "largest group, some on more than one."]
    yield "folds_variants", [
        f"**Variant count.** The first plan gives {s('variants')} variants of its eligibility rule, counted, in the "
        "plan's words, before it fixed its numbers. The file times place the count at "
        f"{clock(o[('first', 'variant_table')])}, after the fold builder ({clock(o[('first', 'fold_builder')])}), "
        f"which already holds {st['min_studies']}, {st['min_study_share']} and {st['max_group_share']}, and after "
        f"the fold files ({clock(o[('first', 'fold_files_first')])} to {clock(o[('first', 'fold_files_last')])}), "
        "so the numbers and the folds came first."]
    yield "folds_variants_use", [
        f"The variants admit {s('variants_datasets_min')} to {s('variants_datasets_max')} datasets in "
        f"{s('variants_species_min')} to {s('variants_species_max')} species. The count reads no score: the script "
        f"rebuilds all {s('variants')} from each dataset's three eligibility counts, and all of it comes before the "
        "hash."]
    left = [r for r in rec.eligibility if r["passes"] == "yes" and r["entered"] == "no"]
    assert len(left) == 1 and left[0]["dataset"] == s("left_by_folds"), left
    sp = rec.species()
    counts = sorted(sp.items(), key=lambda kv: (-kv[1], kv[0]))
    many = [f"{c} of {short(k)}" for k, c in counts if c > 1]
    ones = [short(k) for k, c in counts if c == 1]
    yield "folds_left", [
        f"One dataset passes and did not enter, {left[0]['dataset']} ({left[0]['genomes']} genomes): its fold check "
        f"gives \"{left[0]['folds']}\"."]
    yield "folds_species", [
        f"That leaves {s('datasets')} datasets in {s('species')} species: " + ", ".join(many)
        + f", and 1 each of {joined(ones)}."]
    yield "folds_arms", [
        f"The second experiment holds {n(len(rec.plans))} fold plans, {word(st['draws'])} of each kind per dataset. "
        "Random draw 0 is the release's own random folds, and each matched plan copies its study draw's "
        f"{word(st['folds'])} fold sizes and resistant counts; the script checks the first in all "
        f"{ck['random_d0_equal_release_folds']} datasets and the second in all "
        f"{ck['matched_plans_equal_study_counts']} matched plans."]
    rules = dict(x.split(":") for x in st["study_rule"].split(","))
    later = sorted(int(k[1:]) for k, v in rules.items() if v == "B")
    assert rules["d0"] == "A" and later == [1, 2, 3, 4], rules
    yield "folds_placement", [
        f"Two placement rules assign the groups of a study plan to its {word(st['folds'])} folds. Rule A, the first "
        "experiment's placement, takes groups largest first, with a keyed order only among groups of equal size; "
        "rule B takes every group in a keyed random order. Both then put each group in the fold that least raises "
        "the class imbalance summed over folds. Rule A places study draw 0, the first experiment's assignment, and "
        "rule B draws 1 to 4."]
    yield "folds_movement", [
        "The movement figure is the share of a cohort that a rule moves from draw 0, from the plan builder's table "
        f"of draws 1 to 4: a median {share(s('rule_A_median'))} under rule A ({s('rule_A_unchanged')} datasets not "
        f"at all) and {share(s('rule_B_median'))} under rule B ({share(s('rule_B_min'))} to "
        f"{share(s('rule_B_max'))}). The table's file time is {clock(o[('second', 'rule_variation')])}."]
    yield "folds_rule_time", [
        f"**Rule time.** The second plan gives {clock(o[('second', 'stated_rule_choice')])} as the time it fixed, "
        f"before measuring rule A, that rule A moving a median under {st['study_rule_line']} of a cohort from draw 0 "
        "is too little. No file time or log line records that time; only the plan's text gives it. The plan "
        f"builder's file time, {clock(o[('second', 'plan_builder')])}, comes before it, and the second experiment's "
        "review found rule B already in the builder."]
    yield "folds_rule_score", [
        "Neither placement rule nor the movement figure reads a score, but the first experiment's results, complete "
        f"at {clock(o[('first', 'last_fit_end')])}, were on disk when the rule was chosen."]
    yield "folds_checks", [
        f"`checks.tsv` lists the script's checks: every plan keeps copy groups whole "
        f"({n(ck['copy_groups_whole_plans'])} plans), every study plan keeps BioProjects whole "
        f"({ck['study_plans_with_no_study_split']}), and every fold holds both classes in its test and training "
        f"genomes ({n(ck['fold_files_checked'])} fold files). Its rows that begin `inputs_` give the checks run when "
        "the reduced input set was written."]
    yield "folds_shape_lead", [
        "Study folds differ from random folds in size and in class balance. The table gives medians over fold plans "
        "with their range. Largest is the largest test fold's share of the cohort and Training the smallest "
        "training share; Spread is the range of the resistant share over a plan's "
        f"{word(st['folds'])} test folds, and Departure the largest gap between a test fold's resistant share and "
        "the cohort's. Matched plans equal the study draws on every measure, as they copy them."]
    sh = rec.shapes

    def sc(arm, m):
        r = sh[(arm, m)]
        return f"{share(r['median'])} ({share(r['low'])} to {share(r['high'])})"
    names = (("largest", "Largest"), ("least_training", "Training"), ("spread", "Spread"), ("departure", "Departure"))
    yield "folds_shape_table", head("Measure", "Release random", "Random draws", "First study", "Study draws") + [
        row(label, sc("random-d0", m), sc("random", m), sc("study-d0", m), sc("study", m)) for m, label in names]
    ms = rec.mates
    rand = [ms[("random", str(i))] for i in DRAWS]
    mat = [ms[("matched", str(i))] for i in DRAWS]
    stu = [ms[("study", str(i))] for i in DRAWS]
    assert all(num(r["study_share"]) == 0 for r in stu + [ms[("study", "mean")]])

    def rng(rs, key):
        xs = [float(r[key]) for r in rs]
        return f"{plain(min(xs))} to {plain(max(xs))}"
    yield "folds_mates", [
        "No test genome in a study plan has a study mate in training. In random and matched plans "
        f"{share(ms[('random', 'mean')]['study_share'])} and {share(ms[('matched', 'mean')]['study_share'])} of "
        "test genomes have one (median over datasets of the mean over draws), with a median "
        f"{rng(rand, 'median_study_mates')} and {rng(mat, 'median_study_mates')} study mates each over the draws."]
    yield "folds_lineage", [
        "Exploratory: study plans also keep more lineage-cluster mates out of training. The clusters are the ones the "
        f"release's lineage split holds whole; the script checks this in {ck['lineage_splits_clusters_whole']} "
        "datasets. The median count of "
        f"a test genome's lineage mates in training is {rng(stu, 'median_lineage_mates')} over the study draws, "
        f"against {rng(mat, 'median_lineage_mates')} in matched and {rng(rand, 'median_lineage_mates')} in random "
        "plans."]
    yield "folds_lineage_share", [
        f"The share of test genomes with any lineage mate in training is {share(ms[('study', 'mean')]['lineage_share'])} "
        f"in study plans, against {share(ms[('matched', 'mean')]['lineage_share'])} in matched and "
        f"{share(ms[('random', 'mean')]['lineage_share'])} in random plans. Matched minus study is positive in "
        f"{s('lineage_matched_minus_study_positive')} of {s('datasets')} datasets (median "
        f"{float(s('lineage_matched_minus_study_median')):+.3f}), and random minus matched in "
        f"{s('lineage_random_minus_matched_positive')} (median {float(s('lineage_random_minus_matched_median')):+.3f})."]
    assert s("distance_study_identical_pairs") == "1", s("distance_study_identical_pairs")
    yield "folds_distance", [
        "Exploratory: the five study assignments share most of their placement. Between two study draws a median "
        f"{share(s('distance_study_median'))} of a cohort sits in another fold ({share(s('distance_study_min'))} to "
        f"{share(s('distance_study_max'))}, after pairing fold labels for the largest agreement); the share is under "
        f"0.10 in {s('distance_study_under_tenth')} of {s('datasets')} datasets, and in "
        f"{s('distance_study_identical_pairs')} dataset two study draws are identical. Between two random draws the "
        f"median is {share(s('distance_random_median'))}, and between two matched draws "
        f"{share(s('distance_matched_median'))}."]
    yield "folds_pairs", [
        f"So the second experiment's {int(st['draws']) ** 2} study-random draw pairs, each one study draw beside one "
        f"random draw, cross {word(st['draws'])} close study assignments with {word(st['draws'])} random draws, all "
        f"on the same {s('datasets')} datasets in {s('species')} species."]


# ----------------------------------------------------------------- rule

def rule(rec):
    s, st = rec.s, rec.settings
    a, b, c, d = map(int, re.search(r"low above (\d+)/(\d+), high at (\d+)/(\d+)", st["interval"]).groups())
    level = (Fraction(c, d) - Fraction(a, b)) * 100
    assert level.denominator == 1, level
    level = level.numerator
    yield "rule_interval", [
        f"Each contrast is a paired difference per dataset on the same genomes, read as its median over the "
        f"{s('datasets')} datasets with a {level}% species interval. The script computes each interval exactly, over "
        f"all {n(s('multisets'))} multisets of the {s('species')} species drawn with replacement, each at the weight "
        f"of its multinomial count. The first plan fixed {n(st['first_plan_resamples'])} resamples, and the script "
        "uses the exact count, as the first experiment's review did."]
    yield "rule_reading", [
        "A median reads worse when its interval lies wholly on the worse side of zero, better when it lies wholly on "
        "the other side, and not shown otherwise; a bound at exactly 0 reads not shown. Worse is lower recall, "
        "specificity or balanced accuracy, or a higher loss."]
    th = {f: st[f"threshold_{f}"] for f in FAMILIES}
    pairs = {f: len(rec.pairs(f)) for f in FAMILIES}
    yield "rule_cost", [
        "For one learner a contrast gives a cost when recall or the loss reads worse and neither reads better. The "
        f"first plan's answer for the experiment is cost shown when at least {st['experiment_rule_learners']} of the "
        f"{len(LEARNERS)} learners show a cost and none reads better. The second plan asks, per learner, for a cost "
        f"in at least {th['SR']} of {pairs['SR']} study-random draw pairs, {th['SM']} of {pairs['SM']} "
        f"study-matched draw pairs and {th['MR']} of {pairs['MR']} matched-random draw pairs."]
    yield "rule_limits", [
        "The interval holds species resampling only. It holds nothing of the study assignment or the random seed, "
        f"and nothing of other species or datasets. With {s('two_large_datasets')} of the {s('datasets')} datasets "
        f"in two species and {s('species_one_dataset')} species with one dataset each, a reading beyond the interval "
        f"shows only that the median keeps its sign over resampled sets of these {s('species')} species."]


# ----------------------------------------------------------------- first experiment

def first(rec):
    s = rec.s
    readings = {l: {m: rec.c1(l, "study", m)["reading"] for m in ("recall", "specificity", "log_loss")}
                for l in LEARNERS}
    assert all(readings[l]["log_loss"] == "worse" and readings[l]["specificity"] == "worse" for l in LEARNERS)
    rw = [l for l in LEARNERS if readings[l]["recall"] == "worse"]
    rn = [l for l in LEARNERS if readings[l]["recall"] == "not shown"]
    assert set(rw) | set(rn) == set(LEARNERS)
    costs = [l for l in LEARNERS if rec.first[(l, "study", "rule")]["reading"] == "cost"]
    exp = rec.first[("all", "study", "experiment")]
    assert exp["reading"] == "cost shown" and exp["keep"] == s("species"), exp
    yield "first_answer", [
        f"On the {n(s('pairs'))} pairs, {word(rec.settings['folds'])} folds that keep each BioProject whole, under "
        "one assignment, score lower "
        f"than the release's random folds for all {len(LEARNERS)} learners. Log loss and specificity read worse for "
        f"all {len(LEARNERS)}; recall reads worse for {joined(rw)} and is lower, but not shown, for {joined(rn)}. The "
        f"plan's answer is cost shown, with a cost for {len(costs)} of {len(LEARNERS)} learners, and it stands with "
        f"any one of the {s('species')} species left out."]
    yield "first_table", head("Learner", "Recall", "Specificity", "Log loss") + [
        row(l, cell(rec.c1(l, "study", "recall")), cell(rec.c1(l, "study", "specificity")),
            cell(rec.c1(l, "study", "log_loss"))) for l in LEARNERS]
    er = rec.c1("scm", "study", "error_rate")
    moved = [(l, m) for l in LEARNERS for m in ("recall", "specificity", "log_loss")
             if rec.first[(l, "study", m)]["keep"] != s("species")]
    assert len(moved) == 1, moved
    l, m = moved[0]
    without = [short(x) for x in rec.first[(l, "study", m)]["changed_without"].split(",")]
    yield "first_scm", [
        f"scm's log loss rise is a rise of {signed(er['median'])} [{bound(er['low'])}, {bound(er['high'])}] in its "
        f"error rate. With one species left out, {l}'s {m.replace('_', ' ')} reading changes without "
        f"{joined(without).replace(' and ', ' or ')}, and no other reading of the table changes."]
    fc = rec.foldcut
    cut = fc[("scm", "cut")]
    yield "first_cut", [
        "Study folds also differ from random folds in size and in class balance between training and test. The drop "
        f"remains where both are small: on the {s('first_folds_train_ge')} of {s('first_folds')} study folds that "
        f"train on at least as many genomes as the random folds ({n(s('first_pairs_train_ge'))} pairs), and on the "
        f"{cut['folds']} of those whose resistant share lies within {share(s('first_shift_cut'))} of their training "
        f"genomes' ({cut['datasets']} datasets, {n(cut['genomes'])} pairs)."]
    v = {(l, g, m): num(fc[(l, g)][f"pooled_{m}"]) for l in LEARNERS for g in ("train_ge", "cut", "third_low",
                                                                              "third_high")
         for m in ("recall", "specificity", "log_loss")}
    assert all(v[(l, "cut", "specificity")] > v[(l, "train_ge", "specificity")]
               and 0 < v[(l, "cut", "log_loss")] < v[(l, "train_ge", "log_loss")]
               and v[(l, "cut", "recall")] < v[(l, "train_ge", "recall")] < 0 for l in LEARNERS)
    yield "first_cut_compare", [
        f"On the {cut['folds']} folds the specificity drop and the log loss rise are smaller than on the "
        f"{s('first_folds_train_ge')}, and the recall drop is larger, for every learner."]
    assert all(v[(l, "third_high", m)] < v[(l, "third_low", m)] < 0 for l in LEARNERS
               for m in ("recall", "specificity"))
    assert all(v[(l, "third_high", "log_loss")] > v[(l, "third_low", "log_loss")] > 0 for l in LEARNERS)
    yield "first_shift", [
        "A larger class shift is in line with a larger drop in recall and specificity and a larger log loss rise, for "
        "every learner. The lowest and highest thirds of study folds by shift train on similar counts (median ratio "
        f"to the random folds {float(fc[('scm', 'third_low')]['median_ratio']):.3f} and "
        f"{float(fc[('scm', 'third_high')]['median_ratio']):.3f}). The first experiment's review chose these cuts "
        "after reading the results; they carry no rule and no interval and rest on one assignment."]
    yield "first_cut_key", [
        "Each row pools one group of study folds and gives the study folds' figure minus the random folds' on the "
        "same genomes."]
    groups = (("train_ge", "with larger training"), ("cut", "with both"), ("third_low", "in the lowest third"),
              ("third_high", "in the highest third"))
    for l in LEARNERS:
        for g, _ in groups:
            assert fc[(l, g)]["folds"] == fc[("scm", g)]["folds"]
    yield "first_cut_table", head("Learner", "Folds", "Recall", "Specificity", "Log loss") + [
        row(l, f"{fc[(l, g)]['folds']} {label}", signed(fc[(l, g)]["pooled_recall"]),
            signed(fc[(l, g)]["pooled_specificity"]), signed(fc[(l, g)]["pooled_log_loss"]))
        for l in LEARNERS for g, label in groups]
    # scm scores only 0 or 1, so its pooled log loss difference is its error rate difference times the loss per error.
    per = num(s("scm_log_loss_per_error"))
    rises = [signed(num(fc[("scm", g)]["pooled_log_loss"]) / per) for g, _ in groups]
    yield "first_cut_scm", [
        f"scm's {word(len(rises))} log loss rises are rises of {joined(rises)} in its error rate."]
    lin = int(s("lineage_datasets"))
    nl = [l for l in LEARNERS if rec.first[(l, "study-lineage", "log_loss")]["reading"] == "not shown"]
    keep_all = s("species")
    rb = [l for l in LEARNERS if rec.first[(l, "study-lineage", "recall")]["reading"] == "better"
          and rec.first[(l, "study-lineage", "recall")]["keep"] == keep_all]
    sw = [l for l in LEARNERS if rec.first[(l, "study-lineage", "specificity")]["reading"] == "worse"
          and rec.first[(l, "study-lineage", "specificity")]["keep"] == keep_all]
    turn = []
    for l in LEARNERS:
        ms = [m.replace("_", " ") for m in ("recall", "log_loss", "specificity")
              if rec.first[(l, "study-lineage", m)]["keep"] != keep_all]
        if ms:
            turn.append(f"{l}'s {joined(ms)}")
    yield "first_lineage", [
        f"The first plan gave its lineage contrasts no rule. On the {lin} datasets with a lineage split, study minus "
        f"lineage is not shown on log loss for {len(nl)} of {len(LEARNERS)} learners."]
    yield "first_lineage_more", [
        f"Study folds keep more recall than the lineage split for {joined(rb)} and lose more specificity for "
        f"{joined(sw)}, with every species left out. {joined(turn)} readings each change with one species left out."]
    yield "first_lineage_key", [
        "Each cell gives the median of study minus lineage, its reading and the species leave-outs, of "
        f"{s('species')}, that keep the reading."]

    def lc(l, m):
        r = rec.c1(l, "study-lineage", m)
        return f"{signed(r['median'])} {r['reading']}, {rec.first[(l, 'study-lineage', m)]['keep']} of {keep_all}"
    yield "first_lineage_table", head("Learner", "Recall", "Specificity", "Log loss") + [
        row(l, lc(l, "recall"), lc(l, "specificity"), lc(l, "log_loss")) for l in LEARNERS]


# ----------------------------------------------------------------- second experiment

def second(rec):
    s, st = rec.s, rec.settings
    loss_pairs = {l: rec.worse(l, "SR", rec.loss(l)) for l in LEARNERS}
    assert set(loss_pairs.values()) == {25}, loss_pairs
    med = {l: [float(rec.c2(l, "SR", i, rec.loss(l))["median"]) for i in rec.pairs("SR")] for l in LEARNERS}
    sr = len(rec.pairs("SR"))
    yield "second_s1", [
        f"Study folds score lower than random folds on the same genomes for all {len(LEARNERS)} learners, in each of "
        f"the {sr} draw pairs. The loss reads worse in {loss_pairs['scm']} of {sr} draw pairs for every learner. "
        f"Draw pair medians of scm's error rate rise by {signed(min(med['scm']))} to {signed(max(med['scm']))}, and "
        "of log loss by " + joined(f"{signed(min(med[l]))} to {signed(max(med[l]))} for {l}" for l in GRADED) + "."]
    rmed = {l: [float(rec.c2(l, "SR", i, "recall")["median"]) for i in rec.pairs("SR")] for l in LEARNERS}
    assert all(max(x) < 0 for x in rmed.values())
    yield "second_s2", [
        "Recall is lower in every draw pair as a median, by "
        + joined(f"{pts(max(rmed[l]))} to {pts(min(rmed[l]))} for {l}" for l in LEARNERS) + "."]
    yield "second_s3", [
        "Recall reads worse beyond the species interval in "
        + by_count([(l, rec.worse(l, "SR", "recall")) for l in LEARNERS], sr, "draw pairs")
        + ". Specificity reads worse in "
        + by_count([(l, rec.worse(l, "SR", "specificity")) for l in LEARNERS], sr, "draw pairs") + "."]
    total = {f: len(rec.pairs(f)) for f in FAMILIES}
    yield "second_table", head("Learner", "Contrast", "Recall", "Specificity", "Loss", "Draw pairs") + [
        row(l, FAMILY[f], cell(rec.c2(l, f, "mean", "recall")), cell(rec.c2(l, f, "mean", "specificity")),
            cell(rec.c2(l, f, "mean", rec.loss(l))),
            f"recall {rec.worse(l, f, 'recall')} of {total[f]}, loss {rec.worse(l, f, rec.loss(l))} of {total[f]}")
        for l in LEARNERS for f in FAMILIES]
    ex = s("experiment_SR_cost_shown")
    yield "second_answers_lead", [
        f"Under the plan's rule the answers are as below. The first experiment's rule gives cost shown in {ex} of "
        f"{sr} study-random draw pairs."]

    def ans(l, f):
        r = rec.rules[(l, f, "either")]
        return f"{r['cost']} of {r['pairs']}, {r['answer']}"
    yield "second_answers", head("Learner", "Study-random", "Study-matched", "Matched-random") + [
        row(l, ans(l, "SR"), ans(l, "SM"), ans(l, "MR")) for l in LEARNERS]
    mr = {l: {m: rec.c2(l, "MR", "mean", m) for m in ("recall", "specificity", rec.loss(l))} for l in LEARNERS}
    assert all(mr[l]["log_loss"]["reading"] == "worse" for l in GRADED)
    assert mr["scm"]["error_rate"]["reading"] == "not shown"
    rw = [l for l in LEARNERS if mr[l]["recall"]["reading"] == "worse"]
    rn = [l for l in LEARNERS if mr[l]["recall"]["reading"] == "not shown"]
    spw = [l for l in LEARNERS if mr[l]["specificity"]["reading"] == "worse"]
    yield "second_matched", [
        "Matched folds score a little lower than random folds. Log loss rises by "
        + joined(f"{signed(mr[l]['log_loss']['median'])} for {l}" for l in GRADED)
        + ". Recall falls by " + joined(f"{pts(mr[l]['recall']['median'])} for {l}" for l in rw)
        + f", with no difference shown for {joined(rn)}. Specificity reads worse for {joined(spw)}. scm's error rate "
        "shows no difference beyond the species interval."]
    sm = {l: rec.c2(l, "SM", "mean", rec.loss(l)) for l in LEARNERS}
    assert all(r["reading"] == "worse" for r in sm.values())
    sm_spec = [l for l in LEARNERS if rec.c2(l, "SM", "mean", "specificity")["reading"] == "worse"]
    assert sm_spec == list(LEARNERS)
    yield "second_study", [
        "Study folds score lower than matched folds for every learner on the loss: log loss "
        + joined(f"{signed(sm[l]['median'])} for {l}" for l in GRADED)
        + f", and scm's error rate {signed(sm['scm']['median'])}. Recall reads worse in "
        + by_count([(l, rec.worse(l, "SM", "recall")) for l in LEARNERS], total["SM"], "draw pairs")
        + f". Specificity, over the mean of draws, reads worse for all {len(LEARNERS)} learners."]
    pt = rec.parts
    assert all(r["exact"] == "yes" for r in pt.values())
    rl = [float(pt[(l, m)]["share"]) for l in LEARNERS for m in ("recall", "loss")]
    sp = [float(pt[(l, "specificity")]["share"]) for l in LEARNERS]
    yield "second_parts", [
        "Exploratory: per dataset, study minus random equals study minus matched plus matched minus random, exactly. "
        f"On average over datasets the part not matched by fold shape is {share(min(rl))} to {share(max(rl))} of the "
        f"mean drop, per learner and on recall and the loss; on specificity it is {share(min(sp))} to "
        f"{share(max(sp))}. The table gives each share and the datasets, of {s('datasets')}, where that part is the "
        "larger."]

    def pc(l, m):
        return f"{share(pt[(l, m)]['share'])}, {pt[(l, m)]['larger']} datasets"
    yield "second_parts_table", head("Learner", "Recall", "Specificity", "Loss") + [
        row(l, pc(l, "recall"), pc(l, "specificity"), pc(l, "loss")) for l in LEARNERS]
    cost = {l: int(rec.rules[(l, "MR", "either")]["cost"]) for l in LEARNERS}
    shown = [l for l in LEARNERS if rec.rules[(l, "MR", "either")]["answer"] == "shown"]
    assert shown == ["rf", "gbt"], shown
    g0 = rec.draws[("gbt", "MR", "matched", "cost")]
    rate = Fraction(int(st["threshold_MR"]), 25)
    assert Fraction(g0[0], 5) < rate <= Fraction(sum(g0[1:]), 20), g0
    assert rec.leaveout[("two-large", "gbt", "MR")]["answer"] == "not shown"
    rest = sorted((l for l in LEARNERS if l not in shown), key=lambda l: -cost[l])
    mrn = len(rec.pairs("MR"))
    yield "second_shape", [
        f"The part matched by fold shape, the smaller one, shows a cost for rf in {cost['rf']} of {mrn} "
        f"matched-random draw pairs. gbt shows one in {cost['gbt']} of {mrn}, which holds only with the four more "
        "uneven study assignments and not without the two largest species; both are exploratory readings. "
        + joined(f"{l}'s {cost[l]}" for l in rest) + f" stay under {st['threshold_MR']}."]


# ----------------------------------------------------------------- choices

def runs(xs):
    """[13, 14, 15] -> '13 to 15'."""
    return f"{xs[0]} to {xs[-1]}" if len(xs) > 1 else f"{xs[0]}"


def choices(rec):
    s, st = rec.s, rec.settings
    yield "choices_open", [
        "Every reading under Choices is exploratory apart from the one-species leave-outs, which the second plan "
        "fixed."]
    th = {}
    for r in rec.thresholds:
        th.setdefault(r["family"], []).append((int(r["threshold"]), r["learners"]))
    for f in ("SR", "SM"):
        assert all(ls == ",".join(LEARNERS) for _, ls in th[f]), th[f]
    bands = []
    for t, ls in th["MR"]:
        if bands and bands[-1][1] == ls:
            bands[-1][0].append(t)
        else:
            bands.append(([t], ls))
    mr_text = joined(f"for {joined(ls.split(',')) if len(ls.split(',')) > 1 else ls + ' alone'} at {runs(ts)}"
                     for ts, ls in bands)
    cost = {l: int(rec.rules[(l, "MR", "either")]["cost"]) for l in LEARNERS}
    line = int(st["threshold_MR"])
    near = sorted(LEARNERS, key=lambda l: (abs(cost[l] - line), LEARNERS.index(l)))[:2]
    near = [l for l in LEARNERS if l in near]
    total = {f: len(rec.pairs(f)) for f in FAMILIES}
    yield "choices_thresholds", [
        "**Thresholds.** The study-random answer holds for every learner at any threshold from "
        f"{runs([t for t, _ in th['SR']])} of {total['SR']} draw pairs, and the study-matched answer at "
        f"{runs([t for t, _ in th['SM']])} of {total['SM']}."]
    yield "choices_thresholds_mr", [
        f"The matched-random answer holds {mr_text}. The plan fixed {line} of {total['MR']} before any fit, and "
        + joined(f"{l}'s {cost[l]}" for l in near) + " are the counts nearest to it."]
    yield "choices_draws_lead", [
        "**Draw split.** Study draw 0 is the first experiment's assignment, and matched draw 0 copies it; draws 1 to "
        "4 use rule B, which gives more uneven folds. Each cell gives the draw pairs reading worse with draw 0 of "
        "the study or matched arm, then with draws 1 to 4."]
    axis = {"SR": "study", "MR": "matched"}
    d0 = total["SR"] // len(DRAWS)

    def dc(l, f, m):
        x = rec.draws[(l, f, axis[f], m)]
        return f"{x[0]} of {d0}, {sum(x[1:])} of {total['SR'] - d0}"
    yield "choices_draws_table", head("Learner", "Contrast", "Recall", "Specificity", "Loss", "Cost") + [
        row(l, FAMILY[f], dc(l, f, "recall"), dc(l, f, "specificity"), dc(l, f, "loss"), dc(l, f, "cost"))
        for l in LEARNERS for f in ("SR", "MR")]

    def span(l, first):
        xs = [float(rec.c2(l, "SR", i, rec.loss(l))["median"]) for i in rec.pairs("SR") if (i[0] == "0") == first]
        return min(xs), max(xs)
    parts = []
    for l in LEARNERS:
        a, b = span(l, True), span(l, False)
        assert a[0] <= b[1] and b[0] <= a[1], (l, a, b)
        parts.append((l, f"{signed(a[0])} to {signed(a[1])} against {signed(b[0])} to {signed(b[1])}"))
    yield "choices_draws_loss", [
        "With rule B draws the study-random loss medians per draw pair overlap those of draw 0. With draw 0 against "
        f"draws 1 to 4, scm's error rate is {parts[0][1]}."]
    yield "choices_draws_logloss", [
        "Log loss with draw 0 against draws 1 to 4 is " + joined(f"{t} for {l}" for l, t in parts[1:]) + "."]

    def moved(f, m):
        more, less = [], []
        for l in LEARNERS:
            x = rec.draws[(l, f, axis[f], m)]
            a, b = Fraction(x[0], d0), Fraction(sum(x[1:]), total["SR"] - d0)
            (more if b > a else less if b < a else []).append(l)
        return more, less
    cm, cl = moved("MR", "cost")
    rm, rl = moved("SR", "recall")
    yield "choices_draws_more", [
        f"With rule B draws the matched-random cost shows more often for {joined(cm)} and less often for "
        f"{joined(cl)}, and study-random recall reads worse more often for {joined(rm)} and less often for "
        f"{joined(rl)}."]
    two = s("two_large").split(",")
    lo = rec.leaveout
    rest = int(s("datasets")) - int(s("two_large_datasets"))
    sr = [int(lo[("two-large", l, "SR")]["cost"]) for l in LEARNERS]
    smc = {int(lo[("two-large", l, "SM")]["cost"]) for l in LEARNERS}
    base = {(l, f): rec.rules[(l, f, "either")]["answer"] for l in LEARNERS for f in FAMILIES}
    changed = [(l, f) for l in LEARNERS for f in FAMILIES if lo[("two-large", l, f)]["answer"] != base[(l, f)]]
    assert changed == [("gbt", "MR")] and len(smc) == 1, changed
    yield "choices_two", [
        f"**Two species.** {joined(short(x) for x in two)} hold {s('two_large_datasets')} of the {s('datasets')} "
        f"datasets. With both left out ({rest} datasets in {int(s('species')) - 2} species), gbt's matched-random "
        f"cost falls to {lo[('two-large', 'gbt', 'MR')]['cost']} of {total['MR']} and is not shown, and rf's to "
        f"{lo[('two-large', 'rf', 'MR')]['cost']}, still {lo[('two-large', 'rf', 'MR')]['answer']}. Every other "
        f"answer stands, study-random at {min(sr)} to {max(sr)} of {total['SR']} and study-matched at {smc.pop()} "
        f"of {total['SM']}."]
    one = [(k[0], k[1], r) for k, r in lo.items() if k[0] != "two-large" and k[1] != "all"
           and r["answer"] != base[(k[1], k[2])]]
    yield "choices_one", [
        f"With one species left out, {len(one)} answers change: "
        + joined(f"{l}'s matched-random answer without {short(sp)} ({r['cost']} of {r['pairs']}, {r['answer']})"
                 for sp, l, r in sorted(one, key=lambda x: LEARNERS.index(x[1]))) + "."]
    ru = rec.rules
    both = {f: {l: int(ru[(l, f, "both")]["cost"]) for l in LEARNERS} for f in ("SR", "SM")}
    alone = {f: {l: int(ru[(l, f, "recall")]["cost"]) for l in LEARNERS} for f in ("SR", "SM")}
    assert both == alone
    held_sm = [l for l in LEARNERS if ru[(l, "SM", "both")]["answer"] == "shown"]
    held_sr = [l for l in LEARNERS if ru[(l, "SR", "both")]["answer"] == "shown"]

    def others(f, held_):
        rest_ = [l for l in LEARNERS if l not in held_]
        groups = []
        for l in sorted(rest_, key=lambda l: -both[f][l]):
            if groups and groups[-1][1] == both[f][l]:
                groups[-1][0].append(l)
            else:
                groups.append(([l], both[f][l]))
        return ", ".join(f"{joined(ls)} {c}" for ls, c in groups)
    lasso_bal = ru[("lasso", "MR", "balanced")]
    yield "choices_both", [
        "**Both measures.** Under a stricter rule, a cost only when recall and the loss both read worse, "
        f"study-matched holds for {joined(held_sm)} alone ({both['SM'][held_sm[0]]} of {total['SM']}; "
        f"{others('SM', held_sm)}) and study-random for {joined(held_sr)} ({both['SR'][held_sr[0]]} of "
        f"{total['SR']}; {others('SR', held_sr)}). Recall alone gives the same counts for both contrasts."]
    yield "choices_balanced", [
        "With balanced accuracy in place of recall, lasso's matched-random count is "
        f"{lasso_bal['cost']} of {lasso_bal['pairs']}, {lasso_bal['answer']}, against {cost['lasso']} under the "
        "plan's rule."]
    ns = rec.nostudy
    for l in LEARNERS:
        for m in ("recall", "specificity", "loss"):
            a = abs(float(ns[(l, "SM", m, "no-study")]["mean"]))
            assert all(a < abs(float(ns[(l, "SM", m, side)]["mean"])) for side in ("study", "study-same")), (l, m)
            assert ns[(l, "SM", m, "no-study")]["datasets"] == ns[(l, "SM", m, "study-same")]["datasets"], (l, m)
    yield "choices_nostudy", [
        "**No BioProject.** Exploratory, with no rule in the plan. Genomes without a BioProject have no study mate to "
        f"lose in any arm. On the {s('no_study_datasets')} datasets with at least {st['no_study_min']} of them "
        f"(median {plain(s('no_study_median_genomes'))} genomes, {plain(s('no_study_median_resistant'))} resistant; "
        f"{s('no_study_under_5_resistant')} with under 5 resistant), study minus matched is closer to zero on those "
        "genomes than on genomes with a BioProject, for every learner and measure."]
    yield "choices_nostudy_set", [
        "These genomes are a different and smaller set, and the comparison does not separate study from lineage."]
    dn = {m: ns[("scm", "SM", m, "no-study")]["datasets"] for m in ("recall", "specificity", "loss")}
    assert dn["loss"] == s("no_study_datasets")
    yield "choices_nostudy_key", [
        "Each cell gives the mean over datasets of study minus matched, over the mean of draws, on genomes without a "
        f"BioProject against genomes with one, in the same datasets: {dn['recall']} for recall and "
        f"{dn['specificity']} for specificity, where both classes stand among the genomes without a BioProject, and "
        f"{dn['loss']} for the loss."]

    def nc(l, m):
        return (f"{signed(ns[(l, 'SM', m, 'no-study')]['mean'])} against "
                f"{signed(ns[(l, 'SM', m, 'study-same')]['mean'])}")
    yield "choices_nostudy_table", head("Learner", "Recall", "Specificity", "Loss") + [
        row(l, nc(l, "recall"), nc(l, "specificity"), nc(l, "loss")) for l in LEARNERS]


# ----------------------------------------------------------------- limits and files

def limits(rec):
    """The Limits list, whole and in page order."""
    s, st = rec.s, rec.settings
    d = word(st["draws"])
    return [
        "- Both experiments use one database, BV-BRC, with one curation of its labels and BioProjects.",
        "- A BioProject is an approximate study, and no file gives the laboratory, place or year behind one.",
        "- The matched arm copies fold sizes and class counts only, and study folds also keep about half of the "
        "lineage-cluster mates out of training.",
        f"- The {d} study assignments share most of their placement, and {word(int(st['draws']) - 1)} of them use "
        "a more uneven rule than the first.",
        f"- The species interval rests on {s('species')} species, with {s('two_large_datasets')} of the "
        f"{s('datasets')} datasets in two of them and {s('species_one_dataset')} species with one dataset each.",
        f"- Recall and specificity rest on one cutoff: a score of {st['cutoff']} or more calls resistant.",
        "- scm's figures are its error rate, apart from the log loss in the first experiment's two tables, which the "
        "first plan fixed; both give the same readings.",
        "- Study folds differ from random folds in fold size and in class balance between training and test, and "
        "the first experiment has no matched arm to hold both fixed.",
        "- The page gives no study-split figure for a single dataset.",
        "- The reviews chose the first experiment's fold cuts, the lineage mates, the draw split and the two largest "
        "species left out together after reading the results.",
        f"- The first experiment rests on one study assignment and one random seed, and the second on {d} close "
        f"assignments and {d} random draws.",
        "- Nothing on disk shows what the first experiment opened before its plan, so its statement that it read "
        "no score of the release reports rests on its word.",
        "- No file time or log line records when the second plan fixed its line for rule A.",
        "- The files cannot show whether genomes without a BioProject have study mates.",
    ]


FILES_DIR = ("`<dir>` is the reduced input set, the directory that `study-split-inputs.tar.zst` unpacks to; the "
             "repository does not hold that archive. The record's `inputs.sha256` lists the sha256 of each of its "
             "files.")


def files_items(rec):
    """The Files list without its last item, which counts the held lines."""
    pf = rec.planfiles
    return [
        f"- `reference/study-split/` holds the record: {rec.files} files that `scripts/study-split/record.py` writes "
        "from the reduced inputs. Its README lists each file and what guards it, and `SHA256SUMS` gives each file's "
        "sha256.",
        "- The record holds no genome identifier, label, prediction per genome, model, feature or importance, and "
        f"`checks.tsv` gives {rec.checks['genome_ids_in_record']} as the count of the inputs' genome identifiers "
        "found in its files. Each fold plan enters by its sha256 and fold counts alone.",
        "- The two plans stay in the reduced inputs. `planfiles.tsv` gives each plan's sha256 and the "
        f"{n(pf['first']['hashed_files'])} and {n(pf['second']['hashed_files'])} files the two hash lists hold, "
        "and `--check` hashes the plans again.",
    ]


def held_count(rec):
    """The page's paragraphs and list items the suite holds word for word: HELD, and number-free list items."""
    return len(HELD) + sum(not carries_number(x) for x in limits(rec) + files_items(rec))


def files(rec):
    return files_items(rec) + [
        "- `tests/study_split` rebuilds from the record every paragraph, list item and table of this page that "
        f"carries a number, in figures or words, and holds its other {held_count(rec)} paragraphs and list items "
        "word for word; [fold-file.md](fold-file.md) gives the fold file format.",
    ]


def groups(rec=None):
    rec = rec or Record()
    out = []
    for part in (opening, plans, folds, rule, first, second, choices):
        out.extend(part(rec))
    out.append(("limits", limits(rec)))
    out.append(("files", files(rec)))
    return out


def lines(rec=None):
    return [x for _, ls in groups(rec) for x in ls]


# The page's number-free paragraphs, each held word for word. A deliberate rewording of one is a change here too.
HELD = [
    ("folds_tuning", "Every arm selects settings inside each training set. Tuning keeps copy groups whole and keeps "
                     "no study whole in any arm; the release's lineage split, the first experiment's third arm, "
                     "tunes on lineage groups."),
    ("first_divide", "The design does not measure how the drop divides between study mates, lineage mates and the "
                     "class shift."),
    ("second_key", "The table gives each contrast as the mean over draws: its median over datasets with the species "
                   "interval and reading, beside the draw pairs in which recall and the loss read worse. Loss is the "
                   "error rate for scm and log loss for the others."),
    ("second_line", "The part not matched by fold shape comes with test genomes losing genomes of their own study, "
                    "and about half of the genomes of their own lineage cluster, from training. The design does not "
                    "separate study from lineage."),
    ("files_dir", FILES_DIR),
]


def held(rec=None):
    """The number-free lines the page test holds word for word."""
    rec = rec or Record()
    out = [x for _, x in HELD]
    out += [x for x in limits(rec) + files(rec) if not carries_number(x)]
    return out


def held_groups():
    return {k: [x] for k, x in HELD}
