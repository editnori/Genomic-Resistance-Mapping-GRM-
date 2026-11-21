"""docs/dense-network.md rebuilt from reference/dense-network.

page() gives the whole page. Every figure in it is a cell of the record, rounded, or a count of its
rows. The words around the figures are typed here, so rewording the page means rewording it here too.
The page test holds the page to page() line for line.
"""

import re

from tests.dense_network import common

LIMIT = "This is one experiment on public data. It establishes no gain and rules none out."
FIVE = ["scm", "cart", "rf", "gbt", "lasso"]
FIXED = FIVE + ["site_best", "mean4"]  # comparators that do not read the lineage results
ROWS = FIVE + ["best5", "site_best", "mean4", "baseline", "catalog"]
TABLE_NAME = {"best5": "best of the five", "site_best": "site's best", "mean4": "mean of four",
              "baseline": "lineage baseline", "catalog": "catalog rule"}
TEXT_NAME = {"best5": "the best of the five", "site_best": "the site's best", "mean4": "the mean of four",
             "baseline": "the lineage baseline", "catalog": "the catalog rule"}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def f3(x):
    s = "%.3f" % float(x)
    return "0.000" if s == "-0.000" else s


def s3(x):
    """A difference at 3 decimals, signed; zero unsigned."""
    v = f3(x)
    return v if v.startswith("-") or v == "0.000" else "+" + v


def n(x):
    return "{:,}".format(int(x))


def day(utc):
    y, m, d = utc[:10].split("-")
    return "%d %s %s" % (int(d), MONTHS[int(m) - 1], y)


def clock(utc):
    return utc[11:]


def cells(row):
    return "| " + " | ".join(row) + " |"


def listing(words):
    words = list(words)
    return words[0] if len(words) == 1 else ", ".join(words[:-1]) + " and " + words[-1]


def code(who):
    return "`%s`" % who


def text_name(who):
    return TEXT_NAME.get(who, code(who))


def table_name(who):
    return TABLE_NAME.get(who, code(who))


def sentence_start(text):
    return text[0].upper() + text[1:]


class Record:
    def __init__(self, directory=None):
        d = directory or common.RECORD
        self.datasets = common.read_tsv(d / "datasets.tsv")
        self.measures = {(r["dataset"], r["split"], r["learner"]): r for r in common.read_tsv(d / "measures.tsv")}
        self.summary = {(r["split"], r["comparator"], r["measure"]): r for r in common.read_tsv(d / "summary.tsv")}
        self.rule = common.read_tsv(d / "rule.tsv")
        self.totals = {r["key"]: r["value"] for r in common.read_tsv(d / "totals.tsv")}
        self.settings = {r["key"]: r["value"] for r in common.read_tsv(d / "settings.tsv")}
        self.order = {r["event"]: r for r in common.read_tsv(d / "order.tsv")}
        self.cost = {r["key"]: r["value"] for r in common.read_tsv(d / "cost.tsv")}

    def diff(self, split, who, measure):
        r = self.summary[(split, who, measure)]
        return "%s (%s to %s)" % (s3(r["difference_median"]), s3(r["low"]), s3(r["high"]))

    def below(self, split, among):
        """Comparators whose interval for the network minus them in balanced accuracy lies wholly below 0."""
        return [w for w in among if float(self.summary[(split, w, "balanced_accuracy")]["high"]) < 0]

    def above(self, split, among):
        return [w for w in among if float(self.summary[(split, w, "balanced_accuracy")]["low"]) > 0]

    def per_species(self):
        counts = {}
        for d in self.datasets:
            counts[d["species"]] = counts.get(d["species"], 0) + 1
        return min(counts.values()), max(counts.values())


def page(rec=None):
    rec = rec or Record()
    t, s = rec.totals, rec.settings
    D, S = t["datasets"], t["species"]
    best = rec.summary[("lineage", "best5", "balanced_accuracy")]
    plan, fitting = rec.order["plan_file_time"], rec.order["fitting_script_last_change"]
    first, look = rec.order["first_network_result"], rec.order["analysis_on_first_datasets"]
    last = rec.order["analysis_script_last_change"]
    outcome = '"%s"' % rec.rule[-1]["value"]
    resamples = re.match(r"(\d+) resamples", s["interval"]).group(1)
    seed = re.search(r"seed (\d+)$", s["interval"]).group(1)
    rows = re.search(r"the ([\d,]+) rows", s["inputs"]).group(1)
    pool = re.search(r"the other (\d+)", s["datasets"]).group(1)
    floor = re.search(r"a gap of ([\d.]+) or more", s["datasets"]).group(1)
    lo_sp, hi_sp = rec.per_species()
    out = []
    add = out.append

    add("# Dense network")
    add("")
    add("This page compares one dense feed-forward network with GRM's five learners, `scm`, `cart`, `rf`, `gbt` and "
        "`lasso`, on %s public BV-BRC datasets of release `%s`, on the release's own lineage and random folds. The "
        "network was fitted on %s at settings a plan fixed before its first result; the five were not fitted again, "
        "and their held-out calls are the ones the release saved. %s" % (D, s["grm_release"], day(plan["utc"]), LIMIT))
    add("")

    add("## Reading")
    add("")
    add("Against the plan's comparator, the best of the five picked in hindsight per dataset, which favours the five, "
        "the median network minus comparator lineage balanced accuracy over the %s datasets is %s (%s to %s, a 95%% "
        "interval over the %s species); the network is higher on %s of the %s, and the rule gives %s."
        % (D, s3(best["difference_median"]), s3(best["low"]), s3(best["high"]), S, best["positive"], D, outcome))
    add("")
    add("Against the best of the five in hindsight, the median difference in recall is %s and in log loss %s, over the "
        "same %s datasets." % (rec.diff("lineage", "best5", "recall"), rec.diff("lineage", "best5", "log_loss"), D))
    add("")
    per = listing("%s for %s" % (rec.diff("lineage", w, "balanced_accuracy"), code(w)) for w in FIVE)
    below = rec.below("lineage", FIVE)
    above = rec.above("lineage", FIVE)
    tail = ("only %s's interval lies below 0" % listing(code(w) for w in below)) if below else "no interval lies below 0"
    tail += (", and none above" if not above else ", and %s's above" % listing(code(w) for w in above))
    add("Beside the rule, against each of the five, the median lineage balanced accuracy difference over the %s "
        "datasets is %s; %s." % (D, per, tail))
    add("")
    add("One network at one setting on one screen of %s rows was fitted; no other setting, no other architecture and "
        "no transformer or sequence model was." % rows)
    add("")

    add("## Plan")
    add("")
    add("The plan was fixed on %s. Its file time is %s and its sha256 is `%s`. The network's fitting script was last "
        "changed at %s, after the plan, and the first network result was written at %s. The order rests on file times "
        "alone: no file written before the first result records the plan's sha256."
        % (day(plan["utc"]), clock(plan["utc"]), plan["sha256"], clock(fitting["utc"]), clock(first["utc"])))
    add("")
    add("The plan fixed the network, the datasets and the gates their exports had to pass, the comparator, the "
        "measures, the interval and a four-part rule on the lineage split. It reported the rest %s." % s["reported"])
    add("")

    add("## Network")
    add("")
    for key, label in (("network", "Layers"), ("inputs", "Inputs"), ("start", "Start"), ("loss", "Loss"),
                       ("optimiser", "Optimiser"), ("stopping", "Stopping"), ("repeats", "Repeats"),
                       ("splits", "Splits")):
        add("- **%s.** %s." % (label, sentence_start(s[key])))
    add("")
    add("The screen ran on all of a fold's training samples, so the stopping holdout's labels took part in choosing "
        "the %s rows, as the plan said in advance; no test label did. The network ran %s fits in %s CPU-hours, and "
        "the peak memory of a process was %s kB."
        % (rows, n(rec.cost["fits"]), f3(float(rec.cost["cpu_seconds"]) / 3600), n(rec.cost["peak_rss_kb"])))
    add("")

    add("## Datasets")
    add("")
    add("The gap of a dataset is, %s. The plan took the %s lineage datasets with the largest gap, %s to %s, and drew "
        "%s with seed %s from the other %s with a gap of %s or more. They hold %s to %s genomes, in %s species with %s "
        "to %s datasets each."
        % (s["gap"], t["datasets_largest_gap"], f3(t["gap_largest_group_low"]), f3(t["gap_largest_group_high"]),
           t["datasets_drawn"], seed, pool, floor, n(t["genomes_smallest"]), n(t["genomes_largest"]), S, lo_sp, hi_sp))
    add("")
    add("Every export passed both gates: its genomes, labels and outer folds equal the release's (gate A, %s of %s), "
        "and so does its input fingerprint (gate B, %s of %s)." % (t["gate_a"], D, t["gate_b"], D))
    add("")
    add("The first %s rows are the largest gaps, in the plan's order, and the last %s the drawn. Balanced accuracy is "
        "on the lineage split, where the best of the five is picked per dataset."
        % (t["datasets_largest_gap"], t["datasets_drawn"]))
    add("")
    add(cells(["Dataset", "Genomes", "Resistant", "Gap", "Best of five", "Its balanced accuracy", "Network's"]))
    add("|---|---|---|---|---|---|---|")
    for d in rec.datasets:
        c = rec.measures[(d["dataset"], "lineage", "candidate")]
        b = rec.measures[(d["dataset"], "lineage", "best5")]
        add(cells([d["dataset"], n(d["genomes"]), n(d["resistant"]), f3(d["lineage_gap"]), code(b["picks"]),
                   f3(b["balanced_accuracy"]), f3(c["balanced_accuracy"])]))
    add("")

    add("## Method")
    add("")
    add("%s. The comparators are each of the five; the best of the five in hindsight, the plan's comparator; the "
        "site's best, the learner the release names best, chosen on the random split; the mean of four, the scores of "
        "`cart`, `rf`, `gbt` and `lasso` averaged per genome; the lineage baseline, which calls by lineage identity "
        "inside the release's own outer folds and has balanced accuracy only; and the catalog rule, on the %s of %s "
        "datasets where the release scores it." % (sentence_start(s["calls"]), t["catalog_datasets"], D))
    add("")
    add("The plan named balanced accuracy, recall on the resistant class, specificity, log loss with probabilities "
        "clipped to [1e-4, 1 - 1e-4], and AUROC, each over a dataset's genomes; the experiment's report gives figures "
        "over datasets for the first, second and fourth only. The record works out the same figures for specificity "
        "from the experiment's rows per dataset, by the same draw, so these four are the measures here.")
    add("")
    add("A figure over datasets is the median of the network minus the comparator, with the 2.5%% and 97.5%% percentiles of that "
        "median over %s resamples of the %s species, each taking every dataset of each drawn species, at seed %s."
        % (n(resamples), S, seed))
    add("")
    add("The plan's rule, on the lineage split:")
    add("")
    add(cells(["Rule", "Statistic", "Value", "Needs", "Holds"]))
    add("|---|---|---|---|---|")
    for r in rec.rule[:-1]:
        if r["value"].isdigit():
            value = n(r["value"])
        elif r["rule"] == "4":
            value = f3(r["value"])
        else:
            value = s3(r["value"])
        needs = re.sub(r"-?\d+\.\d{4,}", lambda m: f3(m.group(0)), r["needs"])
        needs = re.sub(r"\b\d{5,}\b", lambda m: n(m.group(0)), needs)
        add(cells([r["rule"], r["statistic"], value, needs, r["holds"]]))
    add("")
    failing = sorted({r["rule"] for r in rec.rule[:-1] if r["holds"] == "no"})
    add("All four must hold to adopt the network. Here rules %s fail, so the outcome is %s."
        % (listing(failing), outcome))
    add("")

    add("## Results")
    add("")
    add("The network minus each comparator: the median over the %s datasets (%s for the catalog rule) with its 95%% "
        "interval over species. A positive log loss difference means the network's is higher." % (D, t["catalog_datasets"]))
    add("")
    for split, word in (("lineage", "Lineage split"), ("random", "Random split")):
        add(word + ":")
        add("")
        add(cells(["Comparator", "Balanced accuracy", "Recall", "Log loss", "Specificity", "Network higher"]))
        add("|---|---|---|---|---|---|")
        for who in ROWS:
            r = rec.summary[(split, who, "balanced_accuracy")]
            figures = [rec.diff(split, who, k) if (split, who, k) in rec.summary else "not scored"
                       for k in ("balanced_accuracy", "recall", "log_loss", "specificity")]
            add(cells([table_name(who)] + figures + ["%s of %s" % (r["positive"], r["datasets"])]))
        add("")
    spec = rec.summary[("lineage", "best5", "specificity")]
    add("Exploratory: against the best of the five on the lineage split, the network's lower balanced accuracy sits "
        "in recall, %s, and not in specificity, %s, where the network is higher on %s of the %s datasets."
        % (rec.diff("lineage", "best5", "recall"), rec.diff("lineage", "best5", "specificity"), spec["positive"],
           spec["datasets"]))
    add("")

    add("## Review")
    add("")
    add("A review on %s recomputed every figure of the experiment's answer with its own code, from the network's "
        "outputs and the release, and found each equal to the 4th decimal. It supports %s for this one network on "
        "these %s datasets, and no wider reading. It narrowed the reading in five ways:" % (day(plan["utc"]), outcome, D))
    add("")
    gaps = sorted(FIVE, key=lambda w: float(t["below_best5_lineage_" + w]))
    add("- **Comparator.** On the lineage split each of the five sits below the best of the five by a median of %s, "
        "so none of the five could have met rule 1 either. The plan fixed this comparator before the results and gave "
        "its reason; these medians were read after them, so they are exploratory."
        % listing("%s (%s)" % (f3(t["below_best5_lineage_" + w]), code(w)) for w in gaps))
    add("- **Site's best.** It has the best random-split balanced accuracy of the five on %s of %s datasets and the "
        "best lineage-split one on %s of %s, so it is a comparator fixed in advance on the lineage split only."
        % (t["site_best_is_best5_random"], D, t["site_best_is_best5_lineage"], D))
    fixed_below = rec.below("lineage", FIXED)
    random_below = rec.below("random", [w for w in ROWS if w not in ("baseline", "catalog")])
    site_high = s3(rec.summary[("lineage", "site_best", "balanced_accuracy")]["high"])
    add("- **Reach.** Against the comparators that do not read the lineage results, the network is shown neither "
        "higher nor lower on the lineage split, except %s, whose interval lies below 0. The figures do not show the "
        "network lower than GRM's learners." % listing(text_name(w) for w in fixed_below))
    add("- **First look.** The analysis ran on the first datasets' results at %s, before the rest ran, and its "
        "script last changed at %s; the experiment's report did not say so. The review's recompute from the plan's "
        "text gives the same figures." % (clock(look["utc"]), clock(last["utc"])))
    add("- **Holdout and wording.** No test label or genome reaches a screen, so the holdout's part cannot raise "
        "test scores. The plan calls its setting a member of the dense network family; the result is about this "
        "network on this screen, and not about the family.")
    add("")
    add("On the random split the intervals for %s lie below 0. On the lineage split the site's best's upper end, %s, "
        "is near enough 0 that another order of species in the draw puts it below."
        % (listing(text_name(w) for w in random_below), site_high))
    add("")

    add("## Limits")
    add("")
    add("- One network at one setting, on one screen and %s datasets. No transformer, sequence model, other "
        "architecture or other setting was fitted, and the record says nothing about them." % D)
    add("- The datasets are not a sample of the release: %s were chosen for the largest gap, and %s drawn once from "
        "those with a gap of %s or more." % (t["datasets_largest_gap"], t["datasets_drawn"], floor))
    add("- The interval resamples %s species with %s to %s datasets each, so it reflects few species." % (S, lo_sp, hi_sp))
    add("- The review did not fit the network again. It read the network's scores as written and checked the "
        "bookkeeping of all %s fits." % n(rec.cost["fits"]))
    add("- The plan named specificity and AUROC among the measures. The experiment's report gives no figure over "
        "datasets for either. The record's specificity figures are its own, and it holds no AUROC.")
    add("- File times show the order of the plan and the results, but they do not prove a file's content at a past "
        "time.")
    add("- Every dataset comes from one public database, so this is not external validation.")
    add("")

    add("## Files")
    add("")
    add("- `reference/dense-network/` holds the record: measures per dataset, split and learner, and medians and "
        "counts over datasets, with no genome, no sequence and no fitted parameter.")
    add("- `scripts/dense-network/record.py` writes it from the experiment's summary file and the review's report, "
        "checks every derived figure against the experiment's own, and rebuilds the derived files from the record "
        "alone.")
    add("- `tests/dense_network` holds this page to the record line for line.")
    add("")
    add("```sh")
    add("python3 -I scripts/dense-network/record.py --check reference/dense-network")
    add("python3 -I -m unittest discover -s tests/dense_network -t .")
    add("```")
    add("")
    add("The experiment's and the review's output directories are not public. With them, `--source` and `--review` "
        "rebuild every file of the record.")
    return "\n".join(out) + "\n"


def reading(text):
    body = text.split("\n## Reading\n", 1)[1].split("\n## ", 1)[0]
    return [line for line in body.split("\n") if line.strip()]


if __name__ == "__main__":
    import sys
    sys.stdout.write(page())
