"""The lines of docs/amrfinderplus-comparison.md that the suite rebuilds from the record, and every spelling of the
record's numbers."""

import re

from tests.amrfinderplus_comparison.common import LEARNERS, RECORD, pairs, read_tsv

PR = "program_minus_known_gene_rule"
NEW = "primary-new-235"
PICKS = LEARNERS + ("served_best", "highest_of_five")
PICK_WORDS = {"served_best": "served learner", "highest_of_five": "highest of the five"}
MEASURE_WORDS = {"ba": "balanced accuracy", "sens": "sensitivity", "spec": "specificity"}
RANDOM_LABEL = "random split: training and test share lineages, and the served learner is picked on this split"
CATALOG = ("Neither the AMRFinderPlus program nor GRM's known-gene rule reads a label, but the catalog both read is "
           "newer than most of these genomes and may have been curated from the same public isolates, which favours "
           "both against the learners.")
MARKER = "a catalog row that meets a rule is a marker and an association."
# The fourth Reading sentence: the catalog sentence, then the marker clause.
FOURTH = CATALOG[:-1] + ", and " + MARKER
ORDER_LABELS = {
    "plan hashed": "The plan hashed",
    "first addendum hashed": "The first addendum hashed",
    "second addendum hashed": "The second addendum hashed",
    "third addendum hashed": "The third addendum hashed",
    "final scoring started": "The final scoring started",
}
SET_TABLES = {
    "sets": ("rule-M-quinolone-28", "rule-M-quinolone-new-26", "quinolone-organism-option-30", "carbapenem-20"),
    "apart": ("no-organism-option-12", "organism-option-no-lineage-split-4"),
}


class Record:
    def __init__(self, root=RECORD):
        self.summary = {(r["set"], r["contrast"], r["measure"]): r for r in read_tsv(root / "summary.tsv")}
        self.sets = {r["set"]: r for r in read_tsv(root / "sets.tsv")}
        self.figures = {r["key"]: r["value"] for r in read_tsv(root / "figures.tsv")}
        self.roles = {r["key"]: r["role"] for r in read_tsv(root / "figures.tsv")}
        self.order = read_tsv(root / "order.tsv")
        self.settings = pairs(root / "settings.tsv")


def lin(pick):
    return "learner_minus_program:lineage:" + pick


def rnd(pick):
    return "learner_minus_program:random:" + pick


def sg(x):
    """A difference as the page writes it: signed, 4 decimals."""
    return "%+.4f" % float(x)


def n(x):
    """An integer as the page writes it."""
    return f"{int(x):,}"


def joined(items):
    items = list(items)
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def interval(r, key):
    lo, hi = r[key + "_species_lo"], r[key + "_species_hi"]
    if lo == "":
        return sg(r[key])
    return "%s (%s to %s)" % (sg(r[key]), sg(lo), sg(hi))


def spans_zero(r):
    return float(r["median_species_lo"]) <= 0 <= float(r["median_species_hi"])


# ----------------------------------------------------------------- built paragraphs

def reading(rec):
    """The four sentences of the section Reading, one paragraph each."""
    st, f = rec.settings, rec.figures
    pr = rec.summary[(NEW, PR, "ba")]
    if not float(pr["mean"]) > float(pr["median"]):
        raise ValueError("the mean is not above the median; the second sentence needs other words")
    if rec.roles["program_ahead_by_0.10_or_more"] != "exploratory" or rec.roles["of_them_quinolone"] != "exploratory":
        raise ValueError("the figures behind the mean are no longer exploratory")
    first = (f"GRM scored three callers on the same genomes and recorded labels of {pr['datasets']} public BV-BRC "
             f"datasets with a lineage split: the AMRFinderPlus program ({st['program_version']}, database "
             f"{st['database_version']}) read under each report's published rule, GRM's own known-gene rule, and GRM's "
             f"five learners scored by cross-validation, all inside that one public database, with no test on data "
             f"from outside it.")
    second = (f"In balanced accuracy, the AMRFinderPlus program read under the report's rule was ahead of GRM's "
              f"known-gene rule in {pr['ahead']} of {pr['datasets']} datasets, level in {pr['level']} and behind in "
              f"{pr['behind']}, with a median difference of {sg(pr['median'])} and a mean of {sg(pr['mean'])}, which "
              f"sits above the median because the program was ahead by {st['large']} or more in "
              f"{f['program_ahead_by_0.10_or_more']} of the {pr['datasets']}, {f['of_them_quinolone']} of them "
              f"quinolone datasets (exploratory).")
    rows = {m: [rec.summary[(NEW, lin(l), m)] for l in LEARNERS] for m in ("ba", "sens", "spec")}
    if not all(spans_zero(r) for r in rows["ba"]):
        raise ValueError("a learner's species interval no longer spans zero")
    if [r["role"] for r in rows["ba"]] != ["primary"] * 5 or {r["role"] for r in rows["sens"] + rows["spec"]} != {"exploratory"}:
        raise ValueError("the roles of the third sentence changed")
    ba = joined(f"{sg(r['median'])} for {l}" for l, r in zip(LEARNERS, rows["ba"]))
    third = (f"On the lineage split of the same {pr['datasets']} datasets, each GRM learner minus the program had a "
             f"median balanced accuracy of {ba} (each {float(st['interval_level']) * 100:.0f}% species interval spans "
             f"zero), with median sensitivity {joined(sg(r['median']) for r in rows['sens'])} and median specificity "
             f"{joined(sg(r['median']) for r in rows['spec'])} in the same order (sensitivity and specificity "
             f"exploratory).")
    return [first, second, third, FOURTH]


def program_paragraph(rec):
    pr = rec.summary[(NEW, PR, "ba")]
    q = rec.summary[("rule-M-quinolone-new-26", PR, "ba")]
    if q["role"] != "exploratory" or rec.roles["level_with_the_same_call_on_every_genome"] != "exploratory":
        raise ValueError("the roles of the program paragraph changed")
    return [f"Of the {pr['level']} datasets level in balanced accuracy, "
            f"{rec.figures['level_with_the_same_call_on_every_genome']} have the same call on every genome "
            f"(exploratory). On the {q['datasets']} quinolone datasets among the {pr['datasets']}, the program minus "
            f"the known-gene rule has a median of {sg(q['median'])} in balanced accuracy (exploratory)."]


def learner_paragraphs(rec):
    st = rec.settings
    sens = [rec.summary[(NEW, lin(l), "sens")] for l in LEARNERS]
    spec = [rec.summary[(NEW, lin(l), "spec")] for l in LEARNERS]
    if not all(float(r["median"]) < 0 and float(r["median_species_hi"]) < 0 for r in sens):
        raise ValueError("a learner's sensitivity is no longer lower with its interval below zero")
    if not all(float(r["median"]) > 0 and float(r["median_species_lo"]) <= 0 for r in spec):
        raise ValueError("a learner's specificity interval no longer starts at or below zero")
    below = [(l, r) for l, r in zip(LEARNERS, spec) if float(r["median_species_lo"]) < 0]
    smed = sorted(float(r["median"]) for r in sens)
    shi = sorted(float(r["median_species_hi"]) for r in sens)
    pmed = sorted(float(r["median"]) for r in spec)
    first = (f"Exploratory: each learner's median sensitivity is lower than the program's ({sg(smed[0])} to "
             f"{sg(smed[-1])}), and every interval ends below zero ({sg(shi[0])} to {sg(shi[-1])}); each learner's "
             f"median specificity is higher than the program's ({sg(pmed[0])} to {sg(pmed[-1])}), and each interval "
             f"starts at zero or just below it ("
             + ", ".join(f"{l} {sg(r['median_species_lo'])}" for l, r in below) + ").")
    r = rec.summary[(NEW, lin("served_best"), "ba")]
    if r["role"] != "primary":
        raise ValueError("the served learner's lineage row is no longer primary")
    w = st["within"]
    second = (f"On the lineage split, the learner the release serves, picked on the random split, was ahead of the "
              f"program in {r['ahead']} of {r['datasets']} datasets, level in {r['level']} and behind in {r['behind']}, "
              f"with a median difference of {sg(r['median'])} (species {float(st['interval_level']) * 100:.0f}% "
              f"{sg(r['median_species_lo'])} to {sg(r['median_species_hi'])}); the two were within {w} of each other in "
              f"{r['within_0.02_either_way']} of {r['datasets']} datasets. It was at most {w} ahead of the program, or "
              f"behind it, in {r['learner_at_most_0.02_ahead']} of {r['datasets']}, which is within {w} as rule M "
              f"defines it (the learner at most {w} ahead). Both counts at {w} are exploratory.")
    third = (f"The served learner has the highest lineage balanced accuracy of the five in "
             f"{rec.figures['served_learner_highest_of_five_on_the_lineage_split']} of {r['datasets']} datasets "
             f"(exploratory).")
    return [first, second, third]


def rule_m_paragraph(rec):
    st, f = rec.settings, rec.figures
    r = rec.summary[("rule-M-quinolone-28", lin("served_best"), "ba")]
    w = st["within"]
    if f["served_learner_behind_by_more_than_0.02_where_the_program_is_at_0.5"] != "2":
        raise ValueError("rule M's edge count is not two; the paragraph needs other words")
    edge = [f["served_learner_balanced_accuracy_%d_of_them" % i] for i in (1, 2)]
    return [f"Rule M, a proposed screen from a review note that no one adopted, counts the quinolone datasets where the "
            f"served learner is at most {w} ahead of the program or behind it. That was "
            f"{r['learner_at_most_0.02_ahead']} of {r['datasets']} where {st['rule_m_needs']} were needed. Two of the "
            f"{r['learner_at_most_0.02_ahead']} count by the served learner scoring below a program at "
            f"{float(st['score_cutoff']):.4f} (balanced accuracy {float(edge[0]):.4f} and {float(edge[1]):.4f}; "
            f"exploratory), and the two were within {w} of each other in {r['within_0.02_either_way']} of "
            f"{r['datasets']}. It is an "
            f"adapted research screen, never an applicable FDA validation criterion."]


def apart_paragraphs(rec):
    s12, s4 = rec.sets["no-organism-option-12"], rec.sets["organism-option-no-lineage-split-4"]
    if (s12["datasets"], s12["species"]) != ("12", "2"):
        raise ValueError("the set without an organism option is not 12 datasets in 2 species")
    return ["Twelve datasets in two species for which the program has no organism option are reported apart: there "
            "the program reads gene rows only.",
            f"The {s4['datasets']} datasets with an organism option and no lineage split are read on the random split "
            f"only, where training and test share lineages."]


def built(rec):
    return (reading(rec) + program_paragraph(rec) + learner_paragraphs(rec) + rule_m_paragraph(rec)
            + apart_paragraphs(rec))


# ----------------------------------------------------------------- tables

def order_table(rec):
    out = ["| Event | Time (UTC) | File | sha256 |", "|---|---|---|---|"]
    for r in rec.order:
        out.append("| %s | %s | %s | %s |" % (ORDER_LABELS[r["event"]], r["utc"].replace("T", " "),
                                               "`%s`" % r["file"] if r["file"] else "",
                                               "`%s`" % r["sha256"][:8] if r["sha256"] else ""))
    return out


def program_table(rec):
    out = ["| Set | Measure | Role | Datasets | Ahead | Level | Behind | Median (95% species interval) "
           "| Mean (95% species interval) |", "|---|---|---|---|---|---|---|---|---|"]
    for name in (NEW, "primary-237"):
        for m in ("ba", "sens", "spec"):
            r = rec.summary[(name, PR, m)]
            out.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                name, MEASURE_WORDS[m], r["role"], r["datasets"], r["ahead"], r["level"], r["behind"],
                interval(r, "median"), interval(r, "mean")))
    return out


def learner_table(rec):
    out = ["| Learner | Role (balanced accuracy) | Ahead, level, behind | Median balanced accuracy (95% species interval) "
           "| Median sensitivity (95% species interval) | Median specificity (95% species interval) |",
           "|---|---|---|---|---|---|"]
    for p in PICKS:
        ba, se, sp = (rec.summary[(NEW, lin(p), m)] for m in ("ba", "sens", "spec"))
        out.append("| %s | %s | %s, %s, %s | %s | %s | %s |" % (
            PICK_WORDS.get(p, p), ba["role"], ba["ahead"], ba["level"], ba["behind"], interval(ba, "median"),
            interval(se, "median"), interval(sp, "median")))
    return out


def random_table(rec):
    out = ["| Learner, %s | Role | Ahead | Level | Behind | Median balanced accuracy (95%% species interval) "
           "| Mean (95%% species interval) |" % RANDOM_LABEL, "|---|---|---|---|---|---|---|"]
    for p in PICKS:
        r = rec.summary[(NEW, rnd(p), "ba")]
        out.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            PICK_WORDS.get(p, p), r["role"], r["ahead"], r["level"], r["behind"], interval(r, "median"),
            interval(r, "mean")))
    return out


def contrast_words(con):
    if con == PR:
        return "program minus known-gene rule"
    if con == lin("served_best"):
        return "served learner minus program, lineage split"
    if con == rnd("served_best"):
        return "served learner minus program, " + RANDOM_LABEL
    raise ValueError(con)


def set_table(rec, names):
    out = ["| Set | Contrast | Role | Datasets | Species | Ahead | Level | Behind | Median | Mean |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for name in names:
        for con in (PR, lin("served_best"), rnd("served_best")):
            r = rec.summary.get((name, con, "ba"))
            if r is None:
                continue
            out.append("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                name, contrast_words(con), r["role"], r["datasets"], r["species"], r["ahead"], r["level"],
                r["behind"], interval(r, "median"), interval(r, "mean")))
    return out


def tables(rec):
    return {"order": order_table(rec), "program": program_table(rec), "learners": learner_table(rec),
            "random": random_table(rec), "sets": set_table(rec, SET_TABLES["sets"]),
            "apart": set_table(rec, SET_TABLES["apart"])}


# ----------------------------------------------------------------- numbers

INT = re.compile(r"0|[1-9]\d*")
FLOAT = re.compile(r"-?\d+\.\d+(e-?\d+)?|-?\d+e-?\d+")
TIME = re.compile(r"(\d{4})-(\d\d)-(\d\d)T(\d\d:\d\d:\d\d)Z")
VERSION = re.compile(r"\d+\.\d+\.\d+|\d{4}-\d\d-\d\d\.\d+")


def spellings(text):
    out = {text}
    if INT.fullmatch(text):
        v = int(text)
        out |= {str(v), f"{v:,}"}
    elif FLOAT.fullmatch(text):
        v = float(text)
        for f in ("%.4f", "%+.4f", "%.2f", "%.1f"):
            out.add(f % v)
        if round(100 * v, 9) == int(round(100 * v, 9)):
            out.add(str(int(round(100 * v, 9))))
        if v == int(v):
            out |= {str(int(v)), f"{int(v):,}"}
    m = TIME.fullmatch(text)
    if m:
        y, mo, d, t = m.groups()
        out |= {y, str(int(d)), f"{y}-{mo}-{d}", t + "Z", f"{y}-{mo}-{d} {t}Z"}
    return out


def values(root=RECORD):
    """Every spelling of every cell of the record's tables, and the record's sha256 values."""
    out, hashes = set(), set()
    for path in sorted(root.glob("*.tsv")):
        for row in read_tsv(path):
            for c in row.values():
                out |= spellings(c)
                if re.fullmatch(r"[0-9a-f]{64}", c):
                    hashes.add(c)
    return out, hashes


NUMBER = re.compile(r"(?<![\w.])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?")
CLOCK = re.compile(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\dZ|\d\d:\d\d:\d\dZ|\d{4}-\d\d-\d\d")


def numbers(line):
    """The numbers a line gives: versions, times, then figures; inline code is left out."""
    line = re.sub(r"`[^`]*`", " ", line)
    found = VERSION.findall(line)
    line = VERSION.sub(" ", line)
    found += CLOCK.findall(line)
    line = CLOCK.sub(" ", line)
    found += NUMBER.findall(line)
    return found
