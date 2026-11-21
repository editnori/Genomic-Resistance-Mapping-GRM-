"""docs/kover-comparison.md rebuilt from reference/kover-comparison.

page() gives the whole page. Every figure in it comes from a cell of the record, or, for the
historical rows, from reference/kover-summary.json, which the page sets apart from the record.
The words around the figures are typed here, so rewording the page means rewording it here too.
The page test holds the page to page() line for line.
"""

import json
import math
import re

from tests.kover_comparison import common

NUMBER = re.compile(r"(?<![A-Za-z0-9_])\d")
LIMIT = "This is a comparison of two programs on {d} public datasets. It shows no equivalence and is not external validation."


def carries_number(line):
    return NUMBER.search(line) is not None


def f3(x):
    s = "%.3f" % float(x)
    return "0.000" if s == "-0.000" else s


def s3(x):
    """A difference, signed, at 3 decimals; at 4 when it is not 0 and rounds to 0.000."""
    if 0 < abs(float(x)) < 0.0005:
        return "%+.4f" % float(x)
    v = f3(x)
    return v if v.startswith("-") or v == "0.000" else "+" + v


def n(x):
    return f"{int(x):,}"


def day(utc):
    y, m, d = utc[:10].split("-")
    month = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"][int(m) - 1]
    return f"{int(d)} {month} {y}"


def clock(utc):
    return utc[11:]


def cells(row):
    return "| " + " | ".join(row) + " |"


class Record:
    def __init__(self, directory=None):
        d = directory or common.RECORD
        self.datasets = common.read_tsv(d / "datasets.tsv")
        self.per_dataset = {r["dataset"]: r for r in common.read_tsv(d / "per_dataset.tsv")}
        self.summary = {r["measure"]: r for r in common.read_tsv(d / "summary.tsv")}
        self.totals = {r["key"]: r["value"] for r in common.read_tsv(d / "totals.tsv")}
        self.cohort = {r["measure"]: r for r in common.read_tsv(d / "cohort.tsv")}
        self.order = common.read_tsv(d / "order.tsv")
        self.settings = {r["key"]: r["value"] for r in common.read_tsv(d / "settings.tsv")}
        self.mapping = common.read_tsv(d / "mapping.tsv")
        self.legacy = json.loads(common.LEGACY.read_text(encoding="utf-8"))

    def event(self, name):
        rows = [r for r in self.order if r["event"] == name]
        return rows[0]

    def figure(self, measure):
        r = self.summary[measure]
        return s3(r["value"]), s3(r["low"]), s3(r["high"])

    def historical(self):
        """The datasets of the record that have a row in reference/kover-summary.json, with both sizes."""
        out = []
        for d in self.datasets:
            species = d["species"]
            drug = d["dataset"].split("__", 1)[1].replace("_", " ")
            rows = [r for r in self.legacy if r["species"] == species and r["antibiotic"].lower() == drug]
            if rows:
                out.append((d["dataset"], int(rows[0]["ds_n_examples"]), int(d["genomes"])))
        return out


def direction(r):
    if float(r["high"]) < 0:
        return "lower"
    if float(r["low"]) > 0:
        return "higher"
    return "not shown to differ in"


def page(rec=None):
    rec = rec or Record()
    t, s, c = rec.totals, rec.settings, rec.cohort
    D, F = t["datasets"], t["folds"]
    version = re.match(r"Kover (\S+)", s["kover_version"]).group(1)
    full_version = s["kover_version"][len("Kover "):]
    ba, ba_lo, ba_hi = rec.figure("balanced_accuracy")
    te, te_lo, te_hi = rec.figure("ties_equal")
    sp, sp_lo, sp_hi = rec.figure("specificity")
    rc, rc_lo, rc_hi = rec.figure("recall")
    ll, ll_lo, ll_hi = rec.figure("log_loss")
    bar = rec.summary["balanced_accuracy"]
    clip = s["log_loss_clip"]
    mult = "%.4f" % math.log((1 - float(clip)) / float(clip))
    resamples = re.match(r"(\d+) resamples", s["bootstrap"]).group(1)
    level = re.search(r"(\d+)% percentile", s["bootstrap"]).group(1)
    seed = re.search(r"--random-seed (\d+)$", s["kover_learn"]).group(1)
    plan = rec.event("plan_hashed")
    fit = rec.event("first_kover_fit")
    review = rec.event("review_plan_hashed")
    ties = rec.event("review_tie_summary")
    cap = re.search(r"kept at (\d+) genomes", s["dataset_rule"]).group(1)
    hist = rec.historical()
    words = {2: "Two", 3: "Three"}

    out = []
    add = out.append

    add("# Kover comparison")
    add("")
    add(f"This page compares Kover {version} with GRM's `scm` learner on {D} public BV-BRC datasets of release "
        f"`{s['grm_release']}` (evaluation revision {s['evaluation_revision']}), the smallest eligible dataset of each "
        f"of {t['species']} species. Both are scored on the release's lineage folds, on the same {n(t['genomes'])} genomes. "
        f"It is not `reference/kover-summary.json`, which holds the {len(rec.legacy)} historical Kover rows of the legacy "
        f"application, on other cohorts and random splits, and which `grm report` reads. " + LIMIT.format(d=D))
    add("")

    add("## Reading")
    add("")
    add(f"On these {D} datasets, Kover's held-out calls have {direction(bar)} balanced accuracy than the released lineage "
        f"calls of GRM's `scm` learner, and the size of the gap turns on how Kover breaks exact ties: Kover minus GRM is "
        f"{ba} ({ba_lo} to {ba_hi}, a {level}% interval over the {D} datasets) as run, and {te} ({te_lo} to {te_hi}) in "
        f"an exploratory reading that counts each tied first rule equally in Kover's {t['kover_one_rule_folds']} "
        f"one-rule folds of {F}.")
    add("")
    add(f"The gap is in specificity, {sp} ({sp_lo} to {sp_hi}), and not in recall, {rc} ({rc_lo} to {rc_hi}); "
        f"log loss, {ll} ({ll_lo} to {ll_hi}), counts the same errors, since every call is 0 or 1.")
    add("")
    excl = int(t["datasets_t_interval_excludes_0"])
    tail = "in no single dataset does its t interval exclude 0" if excl == 0 else f"in {excl} datasets its t interval excludes 0"
    add(f"Kover's balanced accuracy is lower in {bar['datasets_lower']} of the {D} datasets, equal in "
        f"{bar['datasets_equal']} and higher in {bar['datasets_higher']}, and {tail}.")
    add("")
    add("Which differences between the two runs make the gap is not separated, and the result is not shown to carry to larger datasets.")
    add("")

    add("## Plan")
    add("")
    add(f"The plan was fixed on {day(plan['utc'])}: its hash list, recorded at {clock(plan['utc'])}, gives `PLAN.md` the "
        f"sha256 `{plan['sha256']}`. The first Kover fit on a project genome started at {clock(fit['utc'])}, after it.")
    add("")
    add(f"The plan chose the {D} datasets by a rule on sizes alone: of the {c['datasets']['eligible']} datasets of the "
        f"release's {s['release_datasets']} whose `scm` row has a lineage entry, the one with the fewest genomes in each "
        f"species, ties by dataset id, kept at {cap} genomes or fewer. `datasets.tsv` in the record is the file its hash "
        f"list names, byte for byte.")
    add("")
    add(f"It also fixed Kover's settings, the measures and the reading: {s['reading_rule']}.")
    add("")
    add(f"An independent review hashed its own plan at {clock(review['utc'])} and recomputed every figure. It chose the "
        f"tie readings after reading the results per fold, so they are exploratory; its tie summary was written at "
        f"{clock(ties['utc'])}.")
    add("")

    add("## Datasets")
    add("")
    add(f"Each difference is Kover minus GRM, the mean over the dataset's {rec.per_dataset[rec.datasets[0]['dataset']]['folds']} "
        f"lineage folds; the interval for balanced accuracy is a {level}% t interval over those folds.")
    add("")
    add(cells(["Dataset", "Genomes", "Resistant", "Balanced accuracy", "Recall", "Log loss", "Specificity"]))
    add("|---|---|---|---|---|---|---|")
    for d in rec.datasets:
        p = rec.per_dataset[d["dataset"]]
        add(cells([d["dataset"], n(d["genomes"]), n(d["resistant"]),
                   f"{s3(p['balanced_accuracy'])} ({s3(p['balanced_accuracy_t_low'])} to {s3(p['balanced_accuracy_t_high'])})",
                   s3(p["recall"]), s3(p["log_loss"]), s3(p["specificity"])]))
    add("")
    add(f"In {t['datasets_identical_calls']} datasets the two make the same call on every genome.")
    add("")

    add("## Method")
    add("")
    add(f"Kover {full_version} ran with its k-mer tools built against {s['kmer_tools']}, on the decompressed public "
        f"genomes, with these commands for each dataset and each lineage fold k:")
    add("")
    add("```sh")
    add(s["kover_dataset"])
    add(s["kover_split"])
    add(s["kover_learn"])
    add("```")
    add("")
    add("K-mers seen in one genome only are left out, Kover's default. A genome in `test_correct` gets its label and "
        "one in `test_errors` the other class. GRM's side is read from the release and not fitted again: the lineage "
        "`scm` row's calls, each from the fold that did not train on that genome.")
    add("")
    add(f"Log loss takes each call as a probability clipped at {clip}, so it is {clip} plus {mult} times the error rate. "
        f"Each dataset's figure is the mean of its fold differences, and the figure over datasets is their mean, with a "
        f"{level}% interval from {n(resamples)} resamples of whole datasets; `summary.tsv` gives each measure's seed.")
    add("")
    add("The plan's table of settings, with whether Kover's equals GRM's:")
    add("")
    add(cells(["Setting", "GRM `scm`", "Kover as run", "Same"]))
    add("|---|---|---|---|")
    for m in rec.mapping:
        add(cells([m["setting"], m["grm"], m["kover"], m["same"]]))
    add("")
    add("The rows marked no are known sources of difference, and the plan does not separate them.")
    add("")

    add("## Results")
    add("")
    add(f"Over the {D} datasets: the mean of Kover minus GRM with its {level}% bootstrap interval, and the datasets and "
        f"folds where Kover is higher, equal and lower. Recall stands beside log loss.")
    add("")
    add(cells(["Measure", "Kover minus GRM", "Interval", "Datasets", "Folds", "Reading"]))
    add("|---|---|---|---|---|---|")
    for key, label in (("balanced_accuracy", "Balanced accuracy"), ("recall", "Recall"), ("log_loss", "Log loss"),
                       ("specificity", "Specificity"), ("error_rate", "Error rate"), ("rules", "Rules per fold")):
        r = rec.summary[key]
        add(cells([label, s3(r["value"]), f"{s3(r['low'])} to {s3(r['high'])}",
                   f"{r['datasets_higher']}, {r['datasets_equal']}, {r['datasets_lower']}",
                   f"{r['folds_higher']}, {r['folds_equal']}, {r['folds_lower']}",
                   "planned" if r["reading"] == "fixed" else "exploratory"]))
    add("")
    ag = rec.summary["agreeing"]
    add(f"The two make the same call on a mean share of {f3(ag['value'])} of each dataset's genomes "
        f"({f3(ag['low'])} to {f3(ag['high'])}).")
    add("")
    add(f"Exploratory: in {t['ties_more_than_one_rule']} of Kover's {t['kover_one_rule_folds']} one-rule folds, more "
        f"than one rule ties at the first step, Kover's second tie rule keeps them all, and the order of its k-mer "
        f"columns decides. In {t['ties_different_test_results']} of them the tied rules give "
        f"different test results. Only Kover's side of those folds changes below.")
    add("")
    add(cells(["Kover's side", "Kover minus GRM", "Interval"]))
    add("|---|---|---|")
    for key, label in (("balanced_accuracy", "as run, the first tied rule"), ("ties_equal", "each tied rule equally likely"),
                       ("ties_lowest", "the lowest tied result, a bound chosen by test result"),
                       ("ties_highest", "the highest tied result, a bound chosen by test result")):
        v, lo, hi = rec.figure(key)
        add(cells([label, v, f"{lo} to {hi}"]))
    add("")
    add(f"In {t['ties_first_index_equals_run']} of the {t['kover_one_rule_folds']} folds the first tied rule gives "
        f"Kover's own calls. GRM breaks its ties by feature index, an order just as arbitrary, and its fold rules are "
        f"not in the release, so its side cannot be averaged the same way.")
    add("")

    add("## Review")
    add("")
    add("The review reproduced every figure, given the seed the code gives each measure. It narrowed the first report of the comparison in these places:")
    add("")
    add(f"- Size. The first report gave {ba} as the gap. Its size turns on Kover's tie order, as the table above gives, "
        f"and the interval over datasets does not cover that variation.")
    add(f"- Cause. The first report said the gap comes from the parts of the runs that differ by design. With settings "
        f"held equal, covered and lost counts agree at {t['steps_counts_equal']} of {t['steps']} greedy steps of the "
        f"{t['steps_fits']} fits on all genomes, and the rules' patterns at {t['steps_pattern_equal']} of {t['steps']}. "
        f"The fold fits were not checked this way, so those parts are candidate causes and are not separated.")
    add(f"- First step. Where the patterns differ, in {t['step1_pattern_differs']} datasets, the first step is an exact "
        f"tie on both sides, and Kover's rule has the pattern of a rule GRM lists as tied in "
        f"{t['step1_differs_equals_a_listed_tied']} of {t['step1_pattern_differs']}.")
    add("- Reach. The first report said its results do not carry to large datasets. The data show only that they are not shown to carry.")
    add(f"- Log loss. With calls of 0 or 1, log loss is {clip} plus {mult} times the error rate, so {ll} is the error "
        f"rate higher by {f3(rec.summary['error_rate']['value'])} ({f3(rec.summary['error_rate']['low'])} to "
        f"{f3(rec.summary['error_rate']['high'])}), the same result in other units.")
    add("- Seeds. Each measure's bootstrap uses the seed the code gives it, which `summary.tsv` lists, and no interval changes.")
    add("")

    add("## Limits")
    add("")
    add(f"- The {D} datasets are the smallest eligible one of each species: {t['genomes_smallest']} to "
        f"{t['genomes_largest']} genomes, a median of {c['genomes_median']['chosen']} against "
        f"{c['genomes_median']['eligible']} over all {c['datasets']['eligible']} eligible datasets. "
        f"{int(c['species']['eligible']) - int(c['species']['chosen'])} of the {c['species']['eligible']} species have none among them.")
    add(f"- They differ from the rest: a median resistant share of {c['resistant_share_median']['chosen']} against "
        f"{c['resistant_share_median']['eligible']}, {c['resistant_majority']['chosen']} of {D} resistant-majority against "
        f"{c['resistant_majority']['eligible']} of {c['datasets']['eligible']}, and a median of "
        f"{c['lineages_median']['chosen']} lineages against {c['lineages_median']['eligible']}.")
    add(f"- Test folds hold {t['fold_test_genomes_smallest']} to {t['fold_test_genomes_largest']} genomes, so one genome "
        f"moves a fold's figure in a large step, and the interval over datasets rests on {D} points.")
    add("- The interval covers the choice of datasets only. It holds no variation from Kover's tie order or from its seed.")
    add(f"- Three settings favor GRM under balanced accuracy on lineage folds: Kover selects by error rate, its inner folds "
        f"are random and not grouped by lineage, and it allows an empty rule set, which it chose in "
        f"{t['kover_empty_folds']} of {F} folds.")
    add(f"- Kover ran with one seed, {seed}. The spread over seeds is not measured.")
    add(f"- {words[len(hist)]} of the {D} datasets have a row in `reference/kover-summary.json`, on "
        f"{' and '.join(str(h) for _, h, _ in hist)} genomes against {' and '.join(str(g) for _, _, g in hist)} here. "
        f"Those rows come from other cohorts and random splits and are not comparable.")
    add("- No rule is named here. A rule term is a marker and an association.")
    add("")

    add("## Files")
    add("")
    add("- `reference/kover-comparison/` holds the record: counts and measures per dataset, per fold, per greedy step and per measure, with no genome, no sequence and no rule.")
    add("- `scripts/kover-comparison/record.py` writes it from the comparison's reduced files and the review's tie files, and rebuilds the derived files from the record alone.")
    add("- `tests/kover_comparison` holds this page to the record, line for line.")
    add("")
    add("```sh")
    add("python3 -I scripts/kover-comparison/record.py --check reference/kover-comparison")
    add("python3 -I -m unittest discover -s tests/kover_comparison -t .")
    add("```")
    add("")
    add("The comparison's and the review's output directories are not public. With them, `--source` and `--review` rebuild every file of the record.")
    return "\n".join(out) + "\n"


def reading(text):
    """The sentences of the page's Reading section."""
    body = text.split("\n## Reading\n", 1)[1].split("\n## ", 1)[0]
    return [line for line in body.split("\n") if line.strip()]


if __name__ == "__main__":
    import sys
    sys.stdout.write(page())
