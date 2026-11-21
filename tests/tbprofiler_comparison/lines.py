"""The lines of docs/tbprofiler-comparison.md that the suite rebuilds from the record, and every spelling of the
record's numbers."""

import re

from tests.tbprofiler_comparison.common import LEARNERS, RECORD, SPLITS, WAY, pairs, read_tsv

SETS = ("by_name", "graded")
WAYS = ("who_assoc", "who_assoc_only", "tool_default", "who_assoc_small")
CALLERS = ("tool",) + tuple(f"{l}:{s}" for s in SPLITS for l in LEARNERS)
ORDER_LABELS = {
    "plan_hashed": "The plan hashed",
    "addendum_hashed": "The plan's addendum hashed",
    "hashes_checked": "Both hashes checked by the scoring run",
    "last_step_before_scoring": "Last step before the evaluator read a phenotype",
    "scoring_done": "Scoring done",
}
BOUND_WORDS = (("sensitivity_lower", "Sensitivity, no call counted wrong"),
               ("sensitivity_upper", "Sensitivity, no call counted right"),
               ("specificity_lower", "Specificity, no call counted wrong"),
               ("specificity_upper", "Specificity, no call counted right"))


class Record:
    def __init__(self, root=RECORD):
        self.reports = read_tsv(root / "reports.tsv")
        self.report = {r["dataset"]: r for r in self.reports}
        self.callers = {(r["way"], r["dataset"], r["caller"]): r for r in read_tsv(root / "callers.tsv")}
        self.measures = {(r["way"], r["dataset"], r["caller"]): r for r in read_tsv(root / "measures.tsv")}
        self.endpoints = {(r["way"], r["dataset"]): r for r in read_tsv(root / "endpoints.tsv")}
        self.medians = {(r["way"], r["set"], r["caller"], r["measure"]): r for r in read_tsv(root / "medians.tsv")}
        self.comparisons = {(r["way"], r["set"], r["caller"]): r for r in read_tsv(root / "comparisons.tsv")}
        self.screen = {(r["way"], r["set"]): r for r in read_tsv(root / "screen.tsv")}
        self.bounds = {(r["way"], r["set"], r["bound"]): r for r in read_tsv(root / "bounds.tsv")}
        self.intervals = {(r["kind"], r["set"], r["dataset"], r["caller"], r["measure"]): r
                          for r in read_tsv(root / "intervals.tsv")}
        self.disagreements = {r["dataset"]: r for r in read_tsv(root / "disagreements.tsv")}
        self.grades = read_tsv(root / "grades.tsv")
        self.order = read_tsv(root / "order.tsv")
        self.summary = pairs(root / "summary.tsv")
        self.totals = pairs(root / "totals.tsv")
        self.settings = pairs(root / "settings.tsv")
        self.set_aside = pairs(root / "set-aside.tsv")

    def kind(self, kind):
        return [r for r in self.reports if r["kind"] == kind]


def n(x):
    """An integer as the page writes it."""
    return f"{int(x):,}"


def f3(x):
    return "%.3f" % float(x)


def joined(items):
    items = list(items)
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def drug(name):
    return name.replace("_", " ")


def caller_name(c):
    if c == "tool":
        return "TBProfiler"
    learner, split = c.split(":")
    return f"{learner}, {split} split"


def set_name(rec, s):
    count = rec.summary[f"{s}_endpoints"]
    return f"the plan's {count}" if s == "by_name" else f"the addendum's {count}"


def with_interval(value, row):
    """A value at 3 decimals with its lineage interval, or the value alone when the interval has no ends."""
    if row is None or row["low"] == "-":
        return f3(value)
    return f"{f3(value)} [{f3(row['low'])}, {f3(row['high'])}]"


def spread(row, interval=True):
    """A median, its range and, for the primary way, its lineage interval."""
    out = f"{f3(row['median'])} ({f3(row['min'])} to {f3(row['max'])})"
    if interval and row.get("low", "-") != "-":
        out += f" [{f3(row['low'])}, {f3(row['high'])}]"
    return out


# ----------------------------------------------------------------- Reading

def reading(rec):
    """The four sentences of the section Reading, one paragraph each."""
    s, t, st = rec.summary, rec.totals, rec.settings
    a, b = s["by_name_endpoints"], s["graded_endpoints"]
    if t["primary_reports"] != a:
        raise ValueError("a primary report is not scored; the reports whose drug the tool names are not the set")
    if int(s["graded_rows"]) > int(s["by_name_rows"]) or int(b) > int(a):
        raise ValueError("the graded set is not within the set by name")
    for key in ("by_name", "graded"):
        if int(s[f"{key}_rows_with_call"]) + int(s[f"{key}_rows_without_call"]) != int(s[f"{key}_rows"]):
            raise ValueError("rows with a call plus rows without do not make all rows")
    ungraded = [d for d in rec.report if rec.endpoints.get((WAY, d), {}).get("by_name") == "yes"
                and rec.endpoints[(WAY, d)]["graded"] == "no"]
    if len(ungraded) != int(s["ungraded_in_by_name"]):
        raise ValueError("summary.tsv and endpoints.tsv count other ungraded reports")
    for d in ungraded:
        c = rec.callers[(WAY, d, "tool")]
        if c["tp"] != "0" or c["fp"] != "0":
            raise ValueError("the tool calls an ungraded drug resistant; the second sentence's reason fails")
    m = lambda key, measure: rec.medians[(WAY, key, "tool", measure)]
    ns, np_, gs, gp = (m("by_name", "sensitivity"), m("by_name", "specificity"), m("graded", "sensitivity"),
                       m("graded", "specificity"))
    first = (f"GRM ran TBProfiler {st['tool_version']}, a public catalogue-based caller, with the database it ships on "
             f"assemblies of the {n(t['distinct_genomes'])} genomes in its {t['reports_of_the_species']} revision 14 "
             f"M. tuberculosis reports and counted its calls against the phenotypes those reports record on the {a} "
             f"reports whose drug it names, beside the calls each of the five learners made on the same genomes when "
             f"those genomes sat outside its fit.")
    second = (f"On the {a} M. tuberculosis reports whose drug it names, TBProfiler made a call on "
              f"{n(s['by_name_rows_with_call'])} of {n(s['by_name_rows'])} genome-drug rows and none on the other "
              f"{n(s['by_name_rows_without_call'])}, with median sensitivity {f3(ns['median'])} (range {f3(ns['min'])} "
              f"to {f3(ns['max'])}) beside median specificity {f3(np_['median'])} ({f3(np_['min'])} to "
              f"{f3(np_['max'])}), and on the {b} within these {a} whose drug its database grades, with a call on "
              f"{n(s['graded_rows_with_call'])} of {n(s['graded_rows'])} rows, the two medians were {f3(gs['median'])} "
              f"({f3(gs['min'])} to {f3(gs['max'])}) and {f3(gp['median'])} ({f3(gp['min'])} to {f3(gp['max'])}), "
              f"because its database grades no entry for the other {s['ungraded_in_by_name']} drugs and it calls "
              f"neither resistant.")
    rows = [rec.comparisons[(WAY, "graded", f"{l}:lineage")] for l in LEARNERS]
    for r in rows:
        if r["equal"] != "0" or int(r["tool_higher"]) + int(r["tool_lower"]) != int(b):
            raise ValueError("a learner ties the tool on a graded report; 'lower on the rest' no longer holds")
    pieces = [f"{l}'s on {r['tool_higher']}" for l, r in zip(LEARNERS, rows)]
    third = (f"On the lineage split and the same genomes, TBProfiler's balanced accuracy was higher than "
             f"{LEARNERS[0]}'s on {rows[0]['tool_higher']} of the {b} graded M. tuberculosis reports, "
             + ", ".join(pieces[1:-1]) + f" and {pieces[-1]}, and lower on the rest.")
    fourth = (f"These TBProfiler figures are retrospective counts from assemblies, which cannot show a mixed population "
              f"and left {n(s['graded_rows_without_call'])} of the {n(s['graded_rows'])} rows of the {b} graded reports "
              f"without a call, so they say nothing about the tool on sequencing reads, and its catalogue drew on public "
              f"isolates with recorded phenotypes that may include some of these genomes, which would favour it, while "
              f"each learner saw only genomes held out of its fit.")
    return [first, second, third, fourth]


# ----------------------------------------------------------------- built paragraphs

def not_named(rec):
    rows = rec.kind("not named")
    if len(rows) != int(rec.summary["not_named_endpoints"]):
        raise ValueError("reports.tsv and summary.tsv count other reports the tool does not name")
    named = joined(f"{r['antibiotic']} ({n(r['genomes'])} genomes)" for r in rows)
    return (f"The tool names no drug for {len(rows)} reports, {named}. They are counted here and enter no table of "
            f"calls and no median.")


def catalog(rec):
    s = rec.summary
    if int(s["catalog_calls_agree"]) + int(s["catalog_tool_only_right"]) + int(s["catalog_only_right"]) != int(
            s["catalog_rows"]):
        raise ValueError("the catalog row's counts do not sum to its rows")
    return (f"GRM's own reader of the AMRFinderPlus catalog has a row on {s['catalog_endpoints']} of the "
            f"{s['by_name_endpoints']} reports. On them, the tool's balanced accuracy minus the catalog row's has median "
            f"{f3(s['tool_minus_catalog_median'])}. Of the {n(s['catalog_rows'])} genome-drug rows with a call there, "
            f"the two calls agree on {n(s['catalog_calls_agree'])}, only the tool is right on "
            f"{n(s['catalog_tool_only_right'])} and only the catalog row on {n(s['catalog_only_right'])}, which sum to "
            f"{n(s['catalog_rows'])}.")


def set_aside(rec):
    a, s = rec.set_aside, rec.summary
    if a["entries_set_aside"] != "1":
        raise ValueError("the rule set aside other than one entry; the paragraph speaks of one")
    tp, fp, tn, fn = (int(a[f"{k}_as_run"]) for k in ("tp", "fp", "tn", "fn"))
    tp2, fp2, tn2, fn2 = (int(a[f"{k}_counted"]) for k in ("tp", "fp", "tn", "fn"))
    if (tp, fn) != (tp2, fn2):
        raise ValueError("the set-aside entry moves a resistant genome; the paragraph says it does not")
    return (f"**Exploratory.** The tool's rule set aside {a['entries_set_aside']} entry, noted in "
            f"{a['result_files_with_a_note']} result files, on the {a['dataset']} report. As run, the tool called "
            f"{n(tp)} of {n(tp + fn)} resistant genomes resistant and {n(tn)} of {n(tn + fp)} susceptible genomes "
            f"susceptible (sensitivity {'%.4f' % float(s['set_aside_sensitivity_as_run'])}, specificity "
            f"{'%.4f' % float(s['set_aside_specificity_as_run'])}). With the entry counted it calls {n(tn2)} of "
            f"{n(tn2 + fp2)} susceptible genomes susceptible (sensitivity "
            f"{'%.4f' % float(s['set_aside_sensitivity_counted'])}, specificity "
            f"{'%.4f' % float(s['set_aside_specificity_counted'])}), and the best lineage learner minus the tool moves "
            f"from {'%+.4f' % float(s['set_aside_gap_as_run'])} to {'%+.4f' % float(s['set_aside_gap_counted'])}.")


def built(rec):
    return [not_named(rec), catalog(rec), set_aside(rec)]


# ----------------------------------------------------------------- tables

def order_table(rec):
    out = ["| Event | Time (UTC) | Source |", "|---|---|---|"]
    for r in rec.order:
        out.append("| %s | %s | %s |" % (ORDER_LABELS[r["event"]], r["time"].replace("T", " "), r["source"]))
    return out


def grades_table(rec):
    out = ["| Tool drug | Assoc w R | Assoc w R - Interim | Counted entries |", "|---|---|---|---|"]
    for r in rec.grades:
        out.append("| %s | %s | %s | %s |" % (drug(r["tool_drug"]), n(r["assoc_w_r"]), n(r["assoc_w_r_interim"]),
                                             n(r["counted_grade_entries"])))
    return out


def counts_cells(rec, d):
    c = rec.callers[(WAY, d, "tool")]
    m = rec.measures[(WAY, d, "tool")]
    iv = lambda measure: rec.intervals.get(("endpoint", "-", d, "tool", measure))
    return [f"{n(c['with_call'])} of {n(c['genomes'])}", f"{n(c['tp'])} of {n(int(c['tp']) + int(c['fn']))}",
            with_interval(m["sensitivity"], iv("sensitivity")), f"{n(c['tn'])} of {n(int(c['tn']) + int(c['fp']))}",
            with_interval(m["specificity"], iv("specificity")),
            with_interval(m["balanced_accuracy"], iv("balanced_accuracy"))]


def primary_table(rec):
    out = ["| Report | Graded | Rows with a call | Resistant called resistant | Sensitivity | Susceptible called "
           "susceptible | Specificity | Balanced accuracy |", "|---|---|---|---|---|---|---|---|"]
    for r in rec.kind("primary"):
        out.append("| %s | %s | %s |" % (r["antibiotic"], r["graded"], " | ".join(counts_cells(rec, r["dataset"]))))
    return out


def proxy_table(rec):
    out = ["| Report | Read from | Rows with a call | Resistant called resistant | Sensitivity | Susceptible called "
           "susceptible | Specificity | Balanced accuracy |", "|---|---|---|---|---|---|---|---|"]
    for r in rec.kind("proxy"):
        if rec.endpoints[(WAY, r["dataset"])]["by_name"] != "no":
            raise ValueError("a proxy report stands in a set")
        out.append("| %s | %s | %s |" % (r["antibiotic"], joined(drug(x) for x in r["tool_drugs"].split(",")),
                                         " | ".join(counts_cells(rec, r["dataset"]))))
    return out


def medians_table(rec):
    out = ["| Caller | Set | Sensitivity | Specificity | Balanced accuracy |", "|---|---|---|---|---|"]
    for c in CALLERS:
        for s in SETS:
            out.append("| %s | %s | %s |" % (caller_name(c), set_name(rec, s), " | ".join(
                spread(rec.medians[(WAY, s, c, m)]) for m in ("sensitivity", "specificity", "balanced_accuracy"))))
    return out


def comparison_table(rec):
    out = ["| Learner | Set | Tool higher | Tool lower | Tool ahead or within 0.02 | Learner minus tool, median | "
           "Right only for the learner | Right only for the tool |", "|---|---|---|---|---|---|---|---|"]
    for c in CALLERS[1:]:
        for s in SETS:
            r = rec.comparisons[(WAY, s, c)]
            diff = "%+.3f [%+.3f, %+.3f]" % (float(r["median_minus_tool"]), float(r["low"]), float(r["high"]))
            out.append("| %s | %s | %s | %s | %s | %s | %s | %s |" % (
                caller_name(c), set_name(rec, s), r["tool_higher"], r["tool_lower"], r["tool_ahead_or_within"], diff,
                n(r["caller_only_right"]), n(r["tool_only_right"])))
    return out


def gap_table(rec):
    out = ["| Report | Graded | TBProfiler | Best lineage learner | Its balanced accuracy | Gap | "
           "Tool ahead or within 0.02 |",
           "|---|---|---|---|---|---|---|"]
    rows = [rec.endpoints[(WAY, r["dataset"])] for r in rec.reports
            if rec.endpoints.get((WAY, r["dataset"]), {}).get("by_name") == "yes"]
    for e in sorted(rows, key=lambda e: (float(e["gap"]), e["dataset"])):
        out.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            rec.report[e["dataset"]]["antibiotic"], e["graded"], f3(e["tool_balanced_accuracy"]),
            e["best_lineage_learner"], f3(e["best_lineage_balanced_accuracy"]), "%+.4f" % float(e["gap"]),
            e["within_margin"]))
    return out


def ways_table(rec):
    out = ["| Way of reading | Set | TBProfiler sensitivity | TBProfiler specificity | Ahead of or within 0.02 of the "
           "best lineage learner |", "|---|---|---|---|---|"]
    for w in WAYS:
        for s in SETS:
            sc = rec.screen[(w, s)]
            out.append("| %s | %s | %s | %s | %s of %s |" % (
                w, set_name(rec, s), spread(rec.medians[(w, s, "tool", "sensitivity")], False),
                spread(rec.medians[(w, s, "tool", "specificity")], False), sc["within_best"], sc["endpoints"]))
    return out


def disagreement_table(rec):
    out = ["| Report | Resistant call on Assoc w R, label susceptible | Resistant call on Interim alone, label "
           "susceptible | Susceptible call with an entry of another grade, label resistant | Susceptible call with "
           "nothing reported, label resistant | No call, label resistant | No call, label susceptible |",
           "|---|---|---|---|---|---|---|"]
    for r in rec.kind("primary"):
        d = rec.disagreements[r["dataset"]]
        out.append("| %s | %s |" % (r["antibiotic"], " | ".join(n(d[k]) for k in (
            "tool_r_label_s_assoc", "tool_r_label_s_interim", "tool_s_label_r_other", "tool_s_label_r_none",
            "no_call_label_r", "no_call_label_s"))))
    return out


def bounds_table(rec):
    out = ["| Set | " + " | ".join(w for _, w in BOUND_WORDS) + " |", "|---|---|---|---|---|"]
    for s in SETS:
        out.append("| %s | %s |" % (set_name(rec, s), " | ".join(
            spread(rec.bounds[(WAY, s, b)], False) for b, _ in BOUND_WORDS)))
    return out


def tables(rec):
    return {"order": order_table(rec), "grades": grades_table(rec), "primary": primary_table(rec),
            "proxy": proxy_table(rec), "medians": medians_table(rec), "comparison": comparison_table(rec),
            "gap": gap_table(rec), "ways": ways_table(rec), "disagreement": disagreement_table(rec),
            "bounds": bounds_table(rec)}


# ----------------------------------------------------------------- numbers

INT = re.compile(r"0|[1-9]\d*")
FLOAT = re.compile(r"-?\d+\.\d+(e-?\d+)?|-?\d+e-?\d+")
TIME = re.compile(r"(\d{4})-(\d\d)-(\d\d)T(\d\d:\d\d:\d\d)Z")


def spellings(text):
    out = {text}
    if INT.fullmatch(text):
        v = int(text)
        out |= {str(v), f"{v:,}"}
    elif FLOAT.fullmatch(text):
        v = float(text)
        for f in ("%.4f", "%+.4f", "%.3f", "%+.3f", "%.2f"):
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
    """Every spelling of every cell of the record's tables, the resistant and susceptible genomes with a call of
    each row of counts (tp plus fn, tn plus fp), and the record's sha256 values."""
    out, hashes = set(), set()
    for path in sorted(root.glob("*.tsv")):
        for row in read_tsv(path):
            for c in row.values():
                out |= spellings(c)
                if re.fullmatch(r"[0-9a-f]{64}", c):
                    hashes.add(c)
            if {"tp", "fp", "tn", "fn"} <= set(row):
                out |= spellings(str(int(row["tp"]) + int(row["fn"]))) | spellings(str(int(row["tn"]) + int(row["fp"])))
    return out, hashes


VERSION = re.compile(r"(?<![\w.])\d+\.\d+\.\d+(?![\w.])")
NUMBER = re.compile(r"(?<![\w.])[-+]?(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d+)?")
CLOCK = re.compile(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\dZ|\d\d:\d\d:\d\dZ|\d{4}-\d\d-\d\d")


def numbers(line):
    """The numbers a line of prose or a table row gives: times and versions first, then figures; inline code is
    left out."""
    line = re.sub(r"`[^`]*`", " ", line)
    found = CLOCK.findall(line)
    line = CLOCK.sub(" ", line)
    found += VERSION.findall(line)
    line = VERSION.sub(" ", line)
    found += NUMBER.findall(line)
    return found
