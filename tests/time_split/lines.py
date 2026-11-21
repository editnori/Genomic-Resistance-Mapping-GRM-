"""The lines of docs/time-split.md that the suite rebuilds from the record, and every spelling of the record's numbers."""

import re

from tests.time_split.common import ARMS, LEARNERS, RECORD, pairs, read_tsv

MEASURES = ("recall", "log_loss", "balanced_accuracy", "specificity")
ORDER_LABELS = {
    "plan_hashed": "The plan of 14 Nov and its analysis script hashed",
    "first_result": "First result of the run of 14 Nov",
    "layout_hashed": "The table layout script hashed",
    "last_result": "Last result of the run of 14 Nov",
    "rerun_plan_hashed": "The rerun's plan hashed",
    "rerun_first_fit": "First fit of the rerun started",
    "rerun_last_result": "Last result of the rerun",
}


class Record:
    def __init__(self, root=RECORD):
        self.readings = {(r["learner"], r["contrast"], r["measure"]): r for r in read_tsv(root / "readings.tsv")}
        self.arms = {(r["learner"], r["arm"]): r for r in read_tsv(root / "arms.tsv")}
        self.mix = read_tsv(root / "mix.tsv")
        self.species = read_tsv(root / "species.tsv")
        self.order = read_tsv(root / "order.tsv")
        self.review = read_tsv(root / "review.tsv")
        self.rule = {r["learner"]: r for r in read_tsv(root / "rule.tsv")}
        self.summary = pairs(root / "summary.tsv")
        self.totals = pairs(root / "totals.tsv")
        self.settings = pairs(root / "settings.tsv")


def n(x):
    """An integer as the page writes it."""
    return f"{int(x):,}"


def joined(items):
    items = list(items)
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


def cell(r):
    return "%+.4f [%+.4f, %+.4f] %s" % (float(r["median"]), float(r["low"]), float(r["high"]), r["reading"])


def species_name(slug):
    return "*" + slug.replace("_", " ").capitalize() + "*"


def reading(rec):
    """The four sentences of the section Reading, one paragraph each."""
    s, L = rec.summary, len(LEARNERS)
    worse = lambda m: [l for l in LEARNERS if rec.readings[(l, "total", m)]["reading"] == "worse"]
    recall = worse("recall")
    others = [l for l in LEARNERS if l not in recall]
    better = int(s["readings_better"])
    first = (f"In {s['datasets']} datasets of {s['species']} species, scored on later genomes of studies absent from "
             f"training, training only on genomes public by each dataset's cut day read worse than a random draw of the "
             f"same size on log loss for {len(worse('log_loss'))} of {L} learners and on balanced accuracy for "
             f"{len(worse('balanced_accuracy'))} of {L}, and on recall for {len(recall)} of {L} ({joined(recall)}; "
             f"not shown for {joined(others)}).")
    if better == 0:
        first = first[:-1] + "; no learner read better on any measure or contrast."
    pieces = []
    for l in LEARNERS:
        rc, ll = rec.readings[(l, "total", "recall")], rec.readings[(l, "total", "log_loss")]
        pieces.append("%+.4f and %+.4f for %s" % (float(rc["median"]), float(ll["median"]), l))
    second = (f"Median differences over the {s['datasets']} datasets, past minus random, in recall and log loss were "
              + ", ".join(pieces[:-1]) + ", and " + pieces[-1] + ".")
    past, rand = rec.arms[("scm", "past")], rec.arms[("scm", "random")]
    third = (f"scm scores only 0 or 1, so its log loss is a count of wrong calls at the clip: {n(past['wrong_calls'])} "
             f"of {n(past['tested_pairs'])} tested pairs in the past arm against {n(rand['wrong_calls'])} of "
             f"{n(rand['tested_pairs'])} in the random arm.")
    fourth = (f"The plan's rule, fixed before any fit, reads {s['answer']}, with a cost for {s['learners_with_cost']} of "
              f"{s['learners']} learners; the interval resamples only the {s['species']} species, so it is coarse.")
    return [first, second, third, fourth]


def total_table(rec):
    out = ["| Learner | Recall | Log loss | Balanced accuracy | Specificity | Rule |", "|---|---|---|---|---|---|"]
    for l in LEARNERS:
        out.append("| %s | %s | %s |" % (l, " | ".join(cell(rec.readings[(l, "total", m)]) for m in MEASURES),
                                         rec.rule[l]["rule"]))
    return out


def other_table(rec):
    out = ["| Learner | Contrast | Recall | Log loss | Balanced accuracy | Specificity |", "|---|---|---|---|---|---|"]
    for l in LEARNERS:
        for c in ("study", "time"):
            out.append("| %s | %s | %s |" % (l, c, " | ".join(cell(rec.readings[(l, c, m)]) for m in MEASURES)))
    return out


def arms_table(rec):
    out = ["| Learner | Arm | Resistant called resistant | Recall | Susceptible called susceptible | Specificity "
           "| Wrong calls | Log loss from wrong calls at the clip |", "|---|---|---|---|---|---|---|---|"]
    for l in LEARNERS:
        for a in ARMS:
            r = rec.arms[(l, a)]
            susceptible = int(r["tested_pairs"]) - int(r["resistant_pairs"])
            out.append("| %s | %s | %s of %s | %.4f | %s of %s | %.4f | %s | %.1f%% |" % (
                l, a, n(r["tp"]), n(r["resistant_pairs"]), float(r["recall"]), n(r["tn"]), n(susceptible),
                float(r["specificity"]), n(r["wrong_calls"]), 100 * float(r["clip_share"])))
    return out


def mix_table(rec):
    out = ["| Arm | Training genomes | Past | Later, other studies | Later, test studies | Share past |",
           "|---|---|---|---|---|---|"]
    for r in rec.mix:
        out.append("| %s | %s | %s | %s | %s | %.2f%% |" % (r["arm"], n(r["train"]), n(r["from_past"]),
                                                          n(r["from_later_other_studies"]),
                                                          n(r["from_later_test_studies"]), 100 * float(r["past_share"])))
    return out


def species_table(rec):
    out = ["| Species | Datasets | Tested pairs | Resistant pairs | Distinct tested genomes | Test studies "
           "| Genomes of the largest test study |", "|---|---|---|---|---|---|---|"]
    for r in rec.species:
        out.append("| %s | %s | %s | %s | %s | %s | %s |" % (
            species_name(r["species"]), r["datasets"], n(r["tested_pairs"]), n(r["tested_resistant_pairs"]),
            n(r["distinct_test_genomes"]), r["distinct_test_studies"], n(r["largest_study_tested_genomes"])))
    return out


def order_table(rec):
    out = ["| Event | Time (UTC) | Source |", "|---|---|---|"]
    for r in rec.order:
        out.append("| %s | %s | %s |" % (ORDER_LABELS[r["event"]], r["time"].replace("T", " "), r["source"]))
    return out


def review_table(rec):
    out = ["| Learner | Fits changed | Scores changed | Readings changed |", "|---|---|---|---|"]
    for r in rec.review:
        out.append("| %s | %s of %s | %s of %s | %s of %s |" % (
            r["learner"], r["fits_changed"], r["fits"], n(r["scores_changed"]), n(r["scores"]),
            r["readings_changed"], r["readings"]))
    return out


def tables(rec):
    return {"order": order_table(rec), "total": total_table(rec), "other": other_table(rec),
            "arms": arms_table(rec), "mix": mix_table(rec), "species": species_table(rec),
            "review": review_table(rec)}


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
        for f in ("%.4f", "%+.4f", "%.3f", "%.2f", "%.1f"):
            out.add(f % v)
        for f in ("%.2f", "%.1f"):
            out.add(f % (100 * v))
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


def numbers(line):
    """The numbers a line of prose or a table row gives: times first, then figures; inline code is left out."""
    line = re.sub(r"`[^`]*`", " ", line)
    found = []
    for m in re.finditer(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\dZ|\d\d:\d\d:\d\dZ|\d{4}-\d\d-\d\d", line):
        found.append(m.group(0))
    line = re.sub(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\dZ|\d\d:\d\d:\d\dZ|\d{4}-\d\d-\d\d", " ", line)
    found += NUMBER.findall(line)
    return found


def carries_number(line):
    return bool(numbers(line))
