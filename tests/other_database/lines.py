"""Every line of docs/other-database-test.md that carries a number, and the number-free lines it holds, rebuilt from reference/other-database.

lines(record) gives the lines in page order. The page test holds each line of
the page that carries a number to this list, and each line of this list to
the page, once. Words around the figures are typed here; every figure comes
from a cell of the record.

held(record) gives the lines with no number that the page test also holds
word for word, once and in order: the opening's second paragraph, the Plan
section's first paragraph, the specificity paragraph of the Rule section, the
paragraph that reads the kept sentences with three limits, the definitions
of readings L, C and D, and every item of the Limits list. limits(record)
gives the Limits list itself, which the page must give whole and in this
order. Rewording one of these lines means rewording it here too.
"""

import math
import re
import statistics

from . import common

EIGHT = [
    "salmonella_enterica", "klebsiella_pneumoniae", "staphylococcus_aureus", "campylobacter_jejuni", "campylobacter_coli",
    "neisseria_gonorrhoeae", "pseudomonas_aeruginosa", "acinetobacter_baumannii",
]
PILOT = "escherichia_coli"
LEARNERS = [("scm", "SCM"), ("cart", "CART"), ("lasso", "Lasso")]
# The datasets the review of the eight species lists as clear failures, in its order.
CLEAR = [
    ("salmonella_enterica__ciprofloxacin", "scm"), ("salmonella_enterica__ciprofloxacin", "cart"),
    ("klebsiella_pneumoniae__ceftazidime_avibactam", "scm"), ("klebsiella_pneumoniae__ceftazidime_avibactam", "cart"),
    ("neisseria_gonorrhoeae__tetracycline", "scm"), ("klebsiella_pneumoniae__cefpodoxime", "scm"),
    ("salmonella_enterica__nalidixic_acid", "scm"), ("acinetobacter_baumannii__ampicillin_sulbactam", "scm"),
    ("acinetobacter_baumannii__ampicillin_sulbactam", "cart"),
]
NAME = "a test on data from another public database that shares isolates with the training database"
NUMBER = re.compile(r"(?<![A-Za-z0-9_])\d")


def carries_number(line):
    return NUMBER.search(line) is not None


def short(slug):
    genus, species = slug.split("_", 1)
    return f"*{genus[0].upper()}. {species}*"


def drug(dataset_or_antibiotic):
    return dataset_or_antibiotic


def n(x):
    return f"{int(x):,}"


def of(x, y):
    return f"{int(x)} of {int(y)}"


def joined(items):
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + " and " + items[-1]


def joined_comma(items):
    """'a, b, and c': for lists whose items hold 'and' themselves."""
    items = list(items)
    return items[0] if len(items) == 1 else ", ".join(items[:-1]) + ", and " + items[-1]


class Record:
    def __init__(self, directory=None):
        d = directory or common.RECORD
        self.species = {r["species"]: r for r in common.read_tsv(d / "species.tsv")}
        self.datasets = {(r["dataset"], r["learner"]): r for r in common.read_tsv(d / "datasets.tsv")}
        self.readings = {(r["dataset"], r["learner"]): r for r in common.read_tsv(d / "readings.tsv")}
        self.learners = {(r["species"], r["learner"]): r for r in common.read_tsv(d / "learners.tsv")}
        self.totals = {r["learner"]: r for r in common.read_tsv(d / "totals.tsv")}
        self.differences = common.read_tsv(d / "differences.tsv")
        self.order = {r["species"]: r for r in common.read_tsv(d / "order.tsv")}
        self.plans = {r["file"]: r for r in common.read_tsv(d / "plans.tsv")}
        self.settings = {r["key"]: r["value"] for r in common.read_tsv(d / "settings.tsv")}
        self.seeds = common.read_tsv(d / "seeds.tsv")
        self.genome_columns = common.read_tsv(d / "species.tsv", header_only=True)
        self.seed_columns = common.read_tsv(d / "seeds.tsv", header_only=True)

    def est(self, slug, ln):
        return sorted((r for r in self.datasets.values() if r["species"] == slug and r["learner"] == ln and r["estimable"] == "yes"), key=lambda r: r["dataset"])

    def antibiotic(self, dataset):
        return next(r["antibiotic"] for (d, _), r in self.datasets.items() if d == dataset)


def time_of(stamp):
    """2026-11-15T14:30:42Z -> 14:30:42Z."""
    return stamp.split("T")[1]


def band_edges(columns):
    """The distance band edges, from the record's column names: dist_0.0005_0.001 -> ('0.0005', '0.001')."""
    col = next(c for c in columns if c.startswith("dist_") and c.count("_") == 2 and not c.startswith("dist_above"))
    _, lo, hi = col.split("_")
    return lo, hi


def most_shared(rec):
    """The species where most scored genomes are in a BioProject that also contains training genomes."""
    return [sp for sp in [PILOT] + EIGHT if 2 * int(rec.species[sp]["shared_bioproject"]) > int(rec.species[sp]["scored"])]


def limits(rec):
    """The Limits list, whole and in page order."""
    s, t = rec.settings, rec.totals
    eight_low = [r for sp in EIGHT for ln, _ in LEARNERS for r in rec.est(sp, ln) if r["holds"] == "yes" and float(r["specificity"]) < 0.5]
    pilot_low = [r for r in rec.est(PILOT, "scm") + rec.est(PILOT, "cart") + rec.est(PILOT, "lasso")
                 if r["holds"] == "yes" and float(r["specificity"]) < 0.5]
    pilot_above = [r for r in pilot_low if float(r["lineage_specificity"]) > 0.5]
    assert pilot_above == pilot_low, "every pilot result that holds with specificity under 0.5 has a lineage split above it"
    shares = [float(rec.species[sp]["no_accession_share"]) for sp in [PILOT] + EIGHT]
    sa = sorted({r["dataset"] for r in rec.datasets.values() if r["species"] == "staphylococcus_aureus" and r["lineage_figure"] == "no"})
    return [
        f"- It is {NAME}. It is not external validation, and it is not the outside-study test.",
        f"- The rule does not read specificity: in the eight species {len(eight_low)} dataset results hold with specificity under 0.5, "
        f"and in the pilot {len(pilot_above)} more hold with specificity under 0.5 where the lineage split is above it.",
        f"- Every overlap count is a lower bound: {min(shares):.3f} to {max(shares):.3f} of labelled BV-BRC genomes have no accession, and "
        "isolates BV-BRC took in after the data root's copy are not excluded.",
        f"- In {joined(short(sp) for sp in most_shared(rec))} most scored genomes share a BioProject with training genomes.",
        "- The class mix differs from the training cohorts, and log loss depends on it.",
        "- The calls are an adapted research screen.",
        f"- RF and GBT are not scored, and the lasso cannot be read in {t['lasso']['plan_not_assessable']} species.",
        f"- {len(sa)} {short('staphylococcus_aureus')} reports have no lineage figure in the release, so no learner is read on them.",
        "- Reading D shows little where a dataset contains few BioProjects.",
        "- Intervals over genomes treat genomes as independent, and genomes cluster by BioProject.",
    ]


def held(rec=None):
    """The number-free lines the page test holds word for word."""
    rec = rec or Record()
    out = []
    build(rec, out)
    return out


def lines(rec=None):
    return build(rec or Record(), [])


def build(rec, kept):
    """The page's built lines in order; number-free lines held word for word are also put in `kept`."""
    s = rec.settings
    out = []

    def hold(line):
        assert not carries_number(line), line
        out.append(line)
        kept.append(line)

    eight_scored = sum(int(rec.species[sp]["scored"]) for sp in EIGHT)
    pilot_scored = int(rec.species[PILOT]["scored"])
    scored_sets = sorted(s["learners_scored"].split(";"))

    # ---- the test in two paragraphs
    out.append(
        f"The released reports of `{s['release']}` (evaluation revision {s['evaluation_revision']}) were scored on "
        f"{n(pilot_scored + eight_scored)} NCBI Pathogen Detection genomes with antibiogram records: {n(pilot_scored)} {short(PILOT)} genomes "
        f"in a pilot and {n(eight_scored)} in {len(EIGHT)} more species. SCM, CART and the lasso are scored; RF and GBT are not, as the release "
        "contains their importances only."
    )
    assert scored_sets == ["cart", "lasso", "scm"]
    hold(
        f"It is {NAME}. The genomes are outside BV-BRC's labelled set by BioSample and assembly accession, yet in "
        f"{joined(short(sp) for sp in most_shared(rec))} most scored genomes come from a BioProject that also contains training genomes. "
        "It is not external validation."
    )

    # ---- plan and order
    hold(
        "The test ran in two parts, `outside` for the *E. coli* pilot and `outside2` for eight more species, and a review recomputed each "
        "part from its saved files. Neither the parts' output directories nor the reviews' reports are public; the record contains every "
        "figure this page takes from them, and `--check` needs those directories."
    )
    p, a1, a2 = rec.plans["PLAN.md"], rec.plans["ADDENDUM.md"], rec.plans["ADDENDUM2.md"]
    op = rec.order[PILOT]
    assert op["plan_file"] == "PLAN.md" and op["predictions_written"] > p["fixed"] and op["first_after_plan_file"] == "yes"
    assert p["fixed"] < a1["fixed"] < a2["fixed"]
    out.append(
        f"The pilot's plan, `PLAN.md`, was fixed at {p['fixed']}, before the pilot ran. The eight species' addenda were fixed after it: "
        f"`ADDENDUM.md` at {a1['fixed']} and `ADDENDUM2.md`, which adds readings L, C and D, at {a2['fixed']}. Each time is the one written "
        "in the sha256 file beside the plan file, and `plans.tsv` in the record gives each file's sha256."
    )
    o8 = [rec.order[sp] for sp in EIGHT]
    assert all(r["asked_equals_scored"] == "yes" and r["first_at_or_after_predictions"] == "yes" and r["first_after_plan_file"] == "yes" for r in o8)
    out.append(
        f"The review of the eight species assigned each of the {sum(int(r['label_requests']) for r in o8)} label requests "
        f"({sum(int(r['answered_200']) for r in o8)} answered with status 200) to a species by the BioSamples it names. "
        f"In all {len(o8)} species the BioSamples asked equal the scored BioSamples, and the first label request came at or after the time "
        "the species' predictions were written, and after `ADDENDUM2.md` was fixed."
    )
    assert op["first_at_or_after_predictions"] == "yes"
    out.append(
        f"The pilot's {op['label_requests']} label requests ran from {time_of(op['first_label_request'])} to {time_of(op['last_label_request'])}, "
        f"at or after its predictions were written ({time_of(op['predictions_written'])}). Before `PLAN.md` was fixed, "
        f"{op['unrecorded_before_plan']} requests of the pilot went to NCBI Pathogen Detection's search service between "
        f"{time_of(op['unrecorded_first'])} and {time_of(op['unrecorded_last'])} with no parameters saved, and answered "
        f"{joined(op['unrecorded_bytes'].split(';'))} bytes. Whether any carried a phenotype cannot be shown."
    )
    every = [rec.order[sp] for sp in [PILOT] + EIGHT]
    same = sum(r["first_label_request"] == r["predictions_written"] for r in every)
    out.append(
        f"In {same} of the {len(every)} species the first label request came in the same second as the predictions were written, so the clock "
        "alone cannot order the two. Each run's log shows the label step starting after the prediction step ended."
    )

    # ---- genomes
    lo, hi = band_edges(rec.genome_columns)
    out.append(
        "C2 removes isolates with no release drug, C3 those with no GenBank assembly, C4 those in BV-BRC's labelled set by BioSample "
        "or assembly accession, O1 the outside-study candidates, P2 those past the per-species cap, Q1 those failing QC, and D3 those "
        f"at {lo} or less from a training genome on the nearest-training distance."
    )
    out.append("| Species | Isolates | C2 | C3 | C4 | O1 | P2 | Q1 | D3 | Genomes |")
    for sp in [PILOT] + EIGHT:
        r = rec.species[sp]
        cells = [n(r["C1"])] + [n(r[c]) if r[c] != "" else "" for c in ("C2", "C3", "C4", "O1", "P2", "Q1", "D3")] + [n(r["scored"])]
        out.append(f"| {short(sp)} | " + " | ".join(cells) + " |")
    out.append(
        "Share is the largest BioProject's part of the scored genomes and Training the training genomes it contains. Overlap counts scored "
        f"genomes in a BioProject that also contains a training genome, Neighbours those within {hi} of a training genome, and Accession gap the "
        "share of BV-BRC's labelled genomes with neither a BioSample nor an assembly accession. Each overlap count is a lower bound."
    )
    out.append("| Species | Genomes | BioProject | Share | Training | Overlap | Neighbours | Accession gap |")
    for sp in [PILOT] + EIGHT:
        r = rec.species[sp]
        out.append(
            f"| {short(sp)} | {n(r['scored'])} | {r['largest_bioproject']} | {r['largest_bioproject_share']} | {r['largest_bioproject_training']} | "
            f"{n(r['shared_bioproject'])} | {n(r['dist_' + lo + '_' + hi])} | {r['no_accession_share']} |"
        )
    se, cc, cj = rec.species["salmonella_enterica"], rec.species["campylobacter_coli"], rec.species["campylobacter_jejuni"]
    assert cc["largest_bioproject"] == cj["largest_bioproject"]
    out.append(
        f"In {short('salmonella_enterica')}, {n(se['largest_bioproject_scored'])} of {n(se['scored'])} scored genomes are in "
        f"{se['largest_bioproject']}, which also contains {se['largest_bioproject_training']} of the {n(se['training_genomes'])} training genomes. "
        f"In {short('campylobacter_coli')}, {n(cc['largest_bioproject_scored'])} of {n(cc['scored'])} are in {cc['largest_bioproject']}, "
        f"which contains {cc['largest_bioproject_training']} of its {n(cc['training_genomes'])} training genomes. So in these two species the scored "
        "genomes are different isolates, most of them from BioProjects that also contain training genomes."
    )
    out.append(
        f"In {short('campylobacter_jejuni')}, {n(cj['largest_bioproject_scored'])} of {n(cj['scored'])} scored genomes are in the same "
        f"{cj['largest_bioproject']}, which contains {cj['largest_bioproject_training']} of its {n(cj['training_genomes'])} training genomes. So the "
        f"difference between its overlap count and that of {short('campylobacter_coli')} turns on a single training genome."
    )

    # ---- rule
    lr = rec.learners
    out.append(
        f"A dataset is estimable for a learner when it has at least {s['min_class']} resistant and {s['min_class']} susceptible labelled "
        "genomes with a call, and the release has a lineage figure for it. It holds when recall is not below and log loss is not worse "
        f"than the release's lineage split: recall by a {round(float(s['level']) * 100)}% Clopper-Pearson interval, log loss by a "
        f"{round(float(s['level']) * 100)}% bootstrap of {n(s['bootstrap_draws'])} genome resamples with seed {s['bootstrap_seed']}."
    )
    sa = "staphylococcus_aureus"
    no_lineage = sorted({r["dataset"] for r in rec.datasets.values() if r["species"] == sa and r["lineage_figure"] == "no"})
    both = [r for r in rec.datasets.values() if r["species"] == sa and r["learner"] == "lasso"
            and int(r["R"]) >= int(s["min_class"]) and int(r["S"]) >= int(s["min_class"])]
    held_sa = sum(r["holds"] == "yes" for r in both)
    unread = sum(r["estimable"] == "no" for r in both)
    assert lr[(sa, "lasso")]["rule"] == "not assessable" and unread > 0
    # Without the condition the unread datasets would count, holding or not; the result is the same either way.
    from_rule = {rule_of(held_sa, len(both)), rule_of(held_sa + unread, len(both))}
    assert from_rule == {"meets"}, from_rule
    out.append(
        f"The plan's rule asks for {s['min_class']} of each class; the condition that the release has a lineage figure was added after the "
        f"{short(sa)} labels were read, when {len(no_lineage)} of its reports were found to have none, and it moves the lasso in "
        f"{short(sa)} from meeting the rule to not assessable."
    )
    m = re.fullmatch(r"meets when it holds on more than half of at least (\d+) estimable datasets; not assessable under (\d+)", s["pass_rule"])
    out.append(
        f"A learner meets the plan's rule in a species when it holds on more than half of at least {m.group(1)} estimable datasets; "
        f"with fewer it is not assessable. The resistant and susceptible calls are an adapted research screen: SCM's p is its call, "
        f"1 or 0, and CART and the lasso call a genome resistant at p {s['cutoff']} or more. Log loss clips p at {s['clip']}, "
        f"so SCM's log loss is its error rate times {-math.log(float(s['clip'])):.2f}."
    )
    assert "specificity" not in s["holds"] and "recall" in s["holds"] and "log loss" in s["holds"]
    hold(
        "The rule reads recall and log loss only. Specificity is not part of it, so a dataset can hold while its learner calls every "
        "susceptible genome resistant."
    )

    # ---- results: the review's sentences, kept whole
    hold(
        "Read them with three limits given in full below: the rule does not read specificity, intervals over genomes treat genomes as "
        "independent although genomes cluster by BioProject, and log loss depends on the class mix."
    )
    out.append(
        "This is a test on data from another public database: NCBI Pathogen Detection genomes outside BV-BRC's labelled set by BioSample "
        f"and assembly accession, scored with the released reports of `{s['release']}`. In {short('salmonella_enterica')} "
        f"({se['largest_bioproject_share']}) and {short('campylobacter_coli')} ({cc['largest_bioproject_share']}), most scored genomes come "
        "from a BioProject that also holds training genomes. It is not external validation."
    )
    t = rec.totals

    def rule_count(ln):
        return t[ln]["plan_meets"], str(int(t[ln]["plan_meets"]) + int(t[ln]["plan_fails"]))

    def fails(ln):
        return [sp for sp in EIGHT if lr[(sp, ln)]["rule"] == "fails"]

    kp = "klebsiella_pneumoniae"
    poly = ("klebsiella_pneumoniae__polymyxin_b", "scm")
    assert rec.datasets[poly]["FP"] == rec.datasets[poly]["S"] and rec.datasets[poly]["holds"] == "yes"
    assert lr[(kp, "scm")]["margin"] == "1"
    out.append(
        f"SCM meets the plan's rule in {of(*rule_count('scm'))} species. It fails in "
        + " and ".join(f"{short(sp)} ({of(lr[(sp, 'scm')]['holds'], lr[(sp, 'scm')]['estimable'])}{' datasets hold' if i == 0 else ''})" for i, sp in enumerate(fails("scm")))
        + f". Over its {t['scm']['estimable']} estimable datasets, recall is below the lineage figure in {t['scm']['recall_below']} and "
        f"log loss is worse in {t['scm']['log_loss_worse']}. In {short(kp)} the result rests on one dataset "
        f"({of(lr[(kp, 'scm')]['holds'], lr[(kp, 'scm')]['estimable'])} hold), polymyxin B, where SCM calls every susceptible genome resistant."
    )
    out.append(
        f"CART meets the rule in {of(*rule_count('cart'))} species. It fails in "
        + joined(f"{short(sp)} ({of(lr[(sp, 'cart')]['holds'], lr[(sp, 'cart')]['estimable'])})" for sp in fails("cart"))
        + f". Over {t['cart']['estimable']} datasets, recall is below in {t['cart']['recall_below']} and log loss is worse in {t['cart']['log_loss_worse']}."
    )
    lasso_meets = [sp for sp in EIGHT if lr[(sp, "lasso")]["rule"] == "meets"]
    lasso_read = [sp for sp in EIGHT if lr[(sp, "lasso")]["rule"] != "not assessable"]
    poly_lasso = rec.datasets[("klebsiella_pneumoniae__polymyxin_b", "lasso")]
    assert poly_lasso["holds"] == "yes"

    def lasso_cell(sp):
        r = lr[(sp, "lasso")]
        extra = f" (one of them polymyxin B, specificity {poly_lasso['specificity']})" if sp == kp else ""
        return f"{short(sp)} {of(r['holds'], r['estimable'])}{extra}"

    out.append(
        f"The lasso meets the rule in {len(lasso_meets)} of the {len(lasso_read)} species where it can be read: "
        + joined(lasso_cell(sp) for sp in lasso_meets)
        + ". It fails in "
        + joined(f"{short(sp)} ({of(lr[(sp, 'lasso')]['holds'], lr[(sp, 'lasso')]['estimable'])})" for sp in fails("lasso"))
        + f". It is not assessable in {t['lasso']['plan_not_assessable']} species. Over {t['lasso']['estimable']} datasets, recall is below in "
        f"{t['lasso']['recall_below']} and log loss is worse in {t['lasso']['log_loss_worse']}. RF and GBT are not scored."
    )
    for ln, name in LEARNERS:
        tt = t[ln]
        read = int(tt["plan_meets"]) + int(tt["plan_fails"])
        where = f"{read} species" if read == len(EIGHT) else f"the {read} species where it can be read"
        out.append(
            f"**{name}.** Over its {tt['estimable']} estimable datasets in the eight species it holds on {tt['holds']}, recall is below the lineage "
            f"figure in {tt['recall_below']}, log loss is worse in {tt['log_loss_worse']}, and specificity is under 0.5 in "
            f"{tt['specificity_under_half']}. It meets the plan's rule in {tt['plan_meets']} of {where}."
        )
        out.append("| Species | Learner | Datasets | Hold count | Recall | Log loss | Specificity | Rule | Margin | Readings |")
        for sp in EIGHT:
            r = lr[(sp, ln)]
            low = f"{r['lowest_specificity']} ({rec.antibiotic(r['lowest_specificity_dataset'])})" if r["lowest_specificity"] else ""
            readings = ", ".join(f"{v} {r[v + '_rule']}" for v in ("L", "C", "D"))
            out.append(
                f"| {short(sp)} | {name} | {r['estimable']} | {r['holds']} | {r['recall_below']} | {r['log_loss_worse']} | {low} | {r['rule']} | "
                f"{r['margin']} | {readings} |"
            )
    pc, pl, ps = (rec.datasets[("klebsiella_pneumoniae__polymyxin_b", x)] for x in ("cart", "lasso", "scm"))
    assert {pc["holds"], pl["holds"], ps["holds"]} == {"yes"}
    out.append(
        f"In {short(kp)}, polymyxin B holds for all three learners with specificity {ps['specificity']} for SCM, {pc['specificity']} for CART and "
        f"{pl['specificity']} for the lasso, against {ps['lineage_specificity']}, {pc['lineage_specificity']} and {pl['lineage_specificity']} in the "
        f"lineage split. SCM calls all {ps['S']} susceptible genomes resistant, and the dataset holds because its log loss, "
        f"{ps['log_loss']}, is under the lineage figure, {ps['lineage_log_loss']}. With SCM on {of(lr[(kp, 'scm')]['holds'], lr[(kp, 'scm')]['estimable'])} "
        f"and the lasso on {of(lr[(kp, 'lasso')]['holds'], lr[(kp, 'lasso')]['estimable'])}, polymyxin B decides both results."
    )
    pa = "pseudomonas_aeruginosa"
    ca, cs = rec.datasets[(f"{pa}__cefepime", "cart")], rec.datasets[(f"{pa}__cefepime", "scm")]
    assert ca["specificity"] == cs["specificity"] and ca["holds"] == cs["holds"] == "yes"
    assert int(lr[(pa, "scm")]["margin"]) > 1 and int(lr[(pa, "cart")]["margin"]) > 1
    out.append(
        f"{short(pa)} cefepime also holds with specificity {cs['specificity']} for SCM and CART, against {cs['lineage_specificity']} and "
        f"{ca['lineage_specificity']} in the lineage split, and decides neither result (margins {lr[(pa, 'scm')]['margin']} and {lr[(pa, 'cart')]['margin']})."
    )
    close = [(ln, name, sp) for ln, name in LEARNERS for sp in EIGHT if lr[(sp, ln)]["margin"] == "1"]
    assessable = sum(lr[(sp, ln)]["rule"] != "not assessable" for ln, _ in LEARNERS for sp in EIGHT)
    by = {}
    for ln, name, sp in close:
        r = lr[(sp, ln)]
        by.setdefault(name, []).append(f"{short(sp)} ({of(r['holds'], r['estimable'])}{', fails' if r['rule'] == 'fails' else ''})")
    named = {"SCM": "SCM", "CART": "CART", "Lasso": "the lasso"}
    out.append(
        f"{len(close)} of the {assessable} results under the rule are one dataset from the other result: "
        + joined_comma(f"{named[name]} in {joined(cells)}" for name, cells in by.items())
        + "."
    )
    seed_cols = [c for c in rec.seed_columns if c.startswith("seed_")]
    first, others = seed_cols[0], seed_cols[1:]
    assert first == f"seed_{s['bootstrap_seed']}"
    moved = [r for r in rec.seeds if r["level"] == "dataset"]
    assert len(moved) == 1, moved
    mv = moved[0]
    mr = rec.datasets[(mv["dataset"], mv["learner"])]
    assert mv[first] == mr["holds"] == "no"
    out.append(
        f"At seeds {others[0].split('_')[1]} to {others[-1].split('_')[1]} in place of the plan's seed {s['bootstrap_seed']}, one "
        f"{dict(LEARNERS)[mv['learner']]} dataset, {short(mv['species'])} {mr['antibiotic']}, holds at "
        f"{sum(mv[c] == 'yes' for c in others)} of the {len(others)}, and "
        + ("no species result changes." if not any(r["level"] == "species" for r in rec.seeds) else "a species result changes.")
    )

    # ---- clear failures
    out.append("| Dataset | Learner | Recall | Specificity | Log loss |")
    for ds, ln in CLEAR:
        r = rec.datasets[(ds, ln)]
        out.append(
            f"| {short(r['species'])} {r['antibiotic']} | {dict(LEARNERS)[ln]} | {r['recall']} ({r['lineage_recall']}) | "
            f"{r['specificity']} ({r['lineage_specificity']}) | {r['log_loss']} ({r['lineage_log_loss']}) |"
        )

    # ---- readings
    assert "own interval" in s["reading_L"] and "class-balanced" in s["reading_C"] and "stratified" in s["reading_C"]
    assert s["reading_D"].startswith("BioProject cluster bootstrap")
    hold(
        "Reading L keeps the lineage figure with its own interval, C uses a class-balanced log loss with a stratified bootstrap, and D "
        "resamples whole BioProjects. Each keeps the rest of the rule."
    )

    def counts(ln):
        return ", ".join(t[ln][f"{v}_meets"] for v in ("plan", "L", "C", "D")[:-1]) + " and " + t[ln]["D_meets"]

    out.append(
        f"Under the plan's rule and readings L, C and D, SCM meets the rule in {counts('scm')} of {len(EIGHT)} species, CART in {counts('cart')} "
        f"of {len(EIGHT)}, and the lasso in {counts('lasso')} of the {len(lasso_read)} species where it can be read."
    )
    differ = sum(r["differs_in_all_four"] == "yes" for r in rec.differences)
    per_species = sum(1 for r in rec.differences if not r["species"].startswith("count"))
    out.append(
        f"No difference between two learners holds in all four readings. Of the {len(rec.differences)} comparisons, {per_species} of a pair "
        f"in one species and {len(rec.differences) - per_species} of a pair counted over the species, {differ} put the same learner ahead "
        "under the plan's rule and all three readings."
    )
    meds = [f"{lr[(sp, 'scm')]['median_bioprojects']} in {short(sp)}" for sp in EIGHT]
    cc_scm = rec.est("campylobacter_coli", "scm")
    floor = sum(float(rec.readings[(r["dataset"], "scm")]["D_log_loss_lo"]) <= float(s["clip"]) for r in cc_scm)
    out.append(
        "Reading D cannot show a difference where a dataset contains few BioProjects. The median count of BioProjects over SCM's estimable "
        f"datasets is {joined(meds)}."
    )
    ccl = lr[("campylobacter_coli", "scm")]
    out.append(
        f"In {short('campylobacter_coli')}, D turns SCM's {of(ccl['holds'], ccl['estimable'])} into {of(ccl['D_holds'], ccl['estimable'])}, and in "
        f"{floor} of its {len(cc_scm)} datasets the cluster interval of log loss reaches {s['clip']}, the value of a resample with no error."
    )

    # ---- class mix
    out.append("| Species | SCM | CART | Lasso |")
    for sp in EIGHT:
        cells = []
        for ln, _ in LEARNERS:
            r = lr[(sp, ln)]
            cells.append(f"{r['share_below_cohort']} below, {r['share_above_cohort']} above" if int(r["estimable"]) else "")
        out.append(f"| {short(sp)} | " + " | ".join(cells) + " |")
    sel = lr[("salmonella_enterica", "scm")]
    out.append(
        f"SCM's log loss counts its errors, which favours a species with few resistant genomes. In {short('salmonella_enterica')} the "
        f"resistant share is below the cohort's in {of(sel['share_below_cohort'], sel['estimable'])} SCM datasets; under reading C, which "
        f"balances the classes, SCM still holds on {of(sel['C_holds'], sel['estimable'])}."
    )
    ng = "neisseria_gonorrhoeae"
    tet, cip, azi = (rec.datasets[(f"{ng}__{x}", "scm")] for x in ("tetracycline", "ciprofloxacin", "azithromycin"))
    assert lr[(ng, "scm")]["rule"] == "fails"
    assert tet["FP"] == tet["S"] and int(cip["FP"]) > int(cip["FN"]) and int(azi["FN"]) > int(azi["FP"])
    # The mix does not explain tetracycline and ciprofloxacin: their log loss is above the lineage split's at their own share.
    assert all(float(r["log_loss"]) > float(r["lineage_log_loss_at_share"]) for r in (tet, cip))
    out.append(
        f"In {short(ng)} SCM fails. In tetracycline and ciprofloxacin its errors are false resistant calls, which the mix does not explain. "
        f"In tetracycline SCM calls all {tet['S']} susceptible genomes resistant: specificity {tet['specificity']}, against "
        f"{tet['lineage_specificity']} in the lineage split, and log loss {tet['log_loss']}. In ciprofloxacin {cip['FP']} of its "
        f"{int(cip['FP']) + int(cip['FN'])} errors are false resistant calls."
    )
    out.append(
        f"In azithromycin {azi['FN']} of its {int(azi['FN']) + int(azi['FP'])} errors are false susceptible calls, with recall "
        f"{azi['recall']} against {azi['lineage_recall']} in the lineage split. The resistant share is "
        + joined(f"{r['resistant_share']} against {r['cohort_resistant_share']} in {r['antibiotic']}" for r in (azi, cip, tet))
        + "."
    )
    # The mix explains azithromycin: at the lineage split's own rates and this share, log loss is above the one measured.
    azi_c = rec.readings[(azi["dataset"], "scm")]
    assert azi["log_loss_reading"] == "worse" and float(azi["lineage_log_loss_at_share"]) > float(azi["log_loss"]) and azi_c["C_holds"] == "yes"
    out.append(
        f"In azithromycin the mix does explain the worse log loss. At the lineage split's recall, {azi['lineage_recall']}, and specificity, "
        f"{azi['lineage_specificity']}, a resistant share of {azi['resistant_share']} gives an error rate of "
        f"{float(azi['lineage_error_rate_at_share']):.3f} and a log loss of {float(azi['lineage_log_loss_at_share']):.2f}, above the "
        f"{azi['log_loss']} measured. Under reading C, which balances the classes, SCM holds on azithromycin."
    )

    # ---- pilot
    ps_ = lr[(PILOT, "scm")]
    pc_ = lr[(PILOT, "cart")]
    pl_ = lr[(PILOT, "lasso")]
    scm_est = rec.est(PILOT, "scm")
    line = int(ps_["estimable"]) // 2 + 1  # the fewest datasets that hold for the rule to be met
    assert ps_["estimable"] == pc_["estimable"] and int(ps_["holds"]) == line + 1 and int(pc_["holds"]) == line - 1
    out.append(
        f"On {n(pilot_scored)} {short(PILOT)} genomes from NCBI Pathogen Detection whose BioSample and assembly accession match no genome of "
        f"BV-BRC's labelled set, SCM recall was not below its lineage-split figure on {of(int(ps_['estimable']) - int(ps_['recall_below']), ps_['estimable'])} "
        f"estimable datasets and CART recall on {of(int(pc_['estimable']) - int(pc_['recall_below']), pc_['estimable'])}; under the plan's fixed "
        f"rule SCM holds on {ps_['holds']} and CART on {pc_['holds']}, one dataset either side of the line."
    )
    share = statistics.median(float(r["resistant_share"]) for r in scm_est)
    cohort = statistics.median(float(r["cohort_resistant_share"]) for r in scm_est)
    out.append(
        "This margin changes with how the lineage figure's own uncertainty and the clustering of genomes by BioProject are handled, so the "
        "pilot does not show that the two learners differ. The pilot has fewer resistant genomes than the training cohorts (median share "
        f"{share:.2f} against {cohort:.2f}), which lowers log loss, and SCM's log loss is an error rate times {-math.log(float(s['clip'])):.2f}. "
        f"The lasso was scored on {pl_['estimable']} estimable datasets, all sparse fits."
    )
    out.append("| Reading | SCM | CART | Lasso |")
    rows = [("Plan's rule", "holds"), ("L: lineage interval", "L_holds"), ("T: two-sample interval", "T_holds"),
            ("C: balanced log loss", "C_holds"), ("D: BioProject resamples", "D_holds")]
    for label, key in rows:
        out.append(f"| {label} | " + " | ".join(of(lr[(PILOT, ln)][key], lr[(PILOT, ln)]["estimable"]) for ln, _ in LEARNERS) + " |")
    bps = [int(r["bioprojects"]) for r in scm_est]
    out.append(
        f"SCM's estimable pilot datasets contain a median of {lr[(PILOT, 'scm')]['median_bioprojects']} BioProjects, at least {min(bps)}, and the "
        f"pilot's resistant share is below the cohort's in {ps_['share_below_cohort']} of the {ps_['estimable']}."
    )

    # ---- limits
    for item in limits(rec):
        if carries_number(item):
            out.append(item)
        else:
            hold(item)
    return out


def rule_of(holds, estimable):
    """The plan's pass rule, as record.py's rule() gives it."""
    if estimable < 3:
        return "not assessable"
    return "meets" if 2 * holds > estimable else "fails"
