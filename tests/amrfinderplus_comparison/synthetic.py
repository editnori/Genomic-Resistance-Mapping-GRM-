"""SYNTHETIC inputs for the record script: a final scoring, its set lists, a plan with three addenda and the run
logs, laid out as the script reads them. Every name, count and score is made up here; no real input is read."""

import gzip
import hashlib
import json
import pathlib
import random
from fractions import Fraction

from tests.amrfinderplus_comparison.common import LEARNERS, load_script

SEED_STRING = "20991231:synthetic:<set>:<variant>:<contrast>:<measure>"
PLAN_TIME = "2026-11-14T10:00:00Z"
TIMES = ["2026-11-14T10:00:00Z", "2026-11-15T09:00:00Z", "2026-11-15T20:00:00Z", "2026-11-16T01:00:00Z"]
FINAL_START = "2026-11-16T05:00:00Z"
READINGS = ["tool_default", "tool_default_intact", "genes_core", "genes_all_scopes"]
# Species: name, organism option. Seven with an option, so primary-new-235 gets a species interval.
SPECIES = [("alpha_one", True), ("beta_two", True), ("gamma_three", True), ("delta_four", True),
           ("epsilon_five", True), ("zeta_six", True), ("eta_seven", True), ("theta_eight", False)]
DRUGS = [("ciprofloxacin", "quinolone"), ("meropenem", "carbapenem"), ("amikacin", "aminoglycoside"),
         ("tetracycline", "tetracycline")]


def plan_text():
    return "\n".join([
        "# A made-up plan",
        "",
        "Release `rev14-274-20261112`, with the newest database 2026-10-29.1 named for a later run.",
        "Learners are read at 0.5 with `>=`; a set with fewer than 6 species gets no interval.",
        "4,000 resamples, seed string `%s`, and the per-dataset seed string `20991231:<dataset id>`." % SEED_STRING,
        "within 0.02 = best learner minus program <= 0.02 in balanced accuracy; at least 14 of the 28.",
        "",
    ])


def datasets():
    """The made-up datasets: (dataset, species, drug_class, organism option, lineage split, read before)."""
    out = []
    for i, (sp, org) in enumerate(SPECIES):
        for j, (drug, cls) in enumerate(DRUGS):
            if not org and j > 1:
                continue
            lin = not (org and i == 6 and j in (0, 2))
            before = org and lin and i == 1 and j == 0
            out.append(("%s__%s" % (sp, drug), sp, cls, org, lin, before))
    return sorted(out)


def f6(x):
    s = "%.6f" % x
    return "0.000000" if s == "-0.000000" else s


def caller_row(ds, caller, rng, program_half=False):
    resistant, susceptible = rng.randint(8, 40), rng.randint(8, 40)
    if program_half:
        tp, tn = resistant, 0
    else:
        tp, tn = rng.randint(0, resistant), rng.randint(0, susceptible)
    sens, spec = Fraction(tp, resistant), Fraction(tn, susceptible)
    row = {"dataset_id": ds, "caller": caller, "genomes": resistant + susceptible, "scored": resistant + susceptible,
           "no_call_resistant": 0, "no_call_susceptible": 0, "resistant": resistant, "susceptible": susceptible,
           "tp": tp, "fn": resistant - tp, "tn": tn, "fp": susceptible - tn, "sensitivity": f6(sens),
           "specificity": f6(spec), "balanced_accuracy": f6((sens + spec) / 2)}
    for k in ("sens_exact", "sens_lineage", "spec_exact", "spec_lineage", "ba_lineage"):
        row[k + "_lo"], row[k + "_hi"] = "0.100000", "0.900000"
    return row


def write_tsv(path, cols, rows, gz=False):
    text = "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if gz:
        with gzip.open(path, "wt", encoding="utf-8", newline="") as fh:
            fh.write(text)
    else:
        path.write_text(text, encoding="utf-8")


CALLER_COLS = ["dataset_id", "caller", "genomes", "scored", "no_call_resistant", "no_call_susceptible", "resistant",
               "susceptible", "tp", "fn", "tn", "fp", "sensitivity", "sens_exact_lo", "sens_exact_hi", "sens_lineage_lo",
               "sens_lineage_hi", "specificity", "spec_exact_lo", "spec_exact_hi", "spec_lineage_lo", "spec_lineage_hi",
               "balanced_accuracy", "ba_lineage_lo", "ba_lineage_hi"]
DIFF_COLS = ["dataset_id", "variant", "contrast", "measure", "first", "second", "first_minus_second", "lineage_lo",
             "lineage_hi", "first_only_right", "second_only_right", "mcnemar_p", "calls_agree", "scored"]
COHORT_COLS = ["set", "planned", "variant", "contrast", "measure", "units", "not_on_every_genome", "first_ahead",
               "second_ahead", "equal", "within_0.02", "median", "q1", "q3", "mean", "min", "max", "sign_test_p",
               "species", "median_species_lo", "median_species_hi", "mean_species_lo", "mean_species_hi",
               "share_within_0.02", "share_within_species_lo", "share_within_species_hi", "lineage_interval_above_zero",
               "lineage_interval_below_zero", "lineage_interval_spans_zero"]


def build(top):
    """Write a made-up comparison directory under top; return the final scoring's directory."""
    record = load_script()
    rng = random.Random(20991231)
    top = pathlib.Path(top)
    final = top / "results" / "final"
    rows = datasets()
    callers, diffs, head, by_dataset = [], [], [], []
    for ds, sp, cls, org, lin, before in rows:
        served = LEARNERS[rng.randrange(5)]
        head.append({"dataset_id": ds, "species_slug": sp, "drug_class": cls,
                     "primary_set": "yes" if org and lin else "no", "read_before_this_plan": "yes" if before else "no",
                     "served_best": served})
        by_dataset.append({"variant": "tool_default", "dataset_id": ds, "point_screened": "yes" if org else "no"})
        by_dataset.append({"variant": "genes_core", "dataset_id": ds, "point_screened": "yes" if org else "no"})
        half = cls == "quinolone" and sp in ("alpha_one", "gamma_three")
        names = ["program:" + r for r in READINGS] + ["grm_known_gene_rule"]
        names += ["learner:%s:%s" % (l, s) for s in (["random", "lineage"] if lin else ["random"]) for l in LEARNERS]
        base = caller_row(ds, "program:tool_default", rng, program_half=half)
        rr, ss = base["resistant"], base["susceptible"]
        here = {}
        for name in sorted(names):
            if half and name.startswith("program:"):
                tp, tn = rr, 0
            elif name.startswith("program:") or (name == "grm_known_gene_rule" and rng.random() < 0.4):
                tp, tn = base["tp"], base["tn"]
            else:
                tp, tn = rng.randint(0, rr), rng.randint(0, ss)
            sens, spec = Fraction(tp, rr), Fraction(tn, ss)
            here[name] = dict(base, caller=name, tp=tp, fn=rr - tp, tn=tn, fp=ss - tn, sensitivity=f6(sens),
                              specificity=f6(spec), balanced_accuracy=f6((sens + spec) / 2))
        callers += [here[k] for k in sorted(here)]
        rates = {k: record.rates(v) for k, v in here.items()}

        def diff(variant, con, m, first, second):
            d = rates[first][m] - rates[second][m]
            scored = here[first]["scored"]
            a = rng.randint(0, scored // 3)
            b = rng.randint(0, scored // 3)
            lo = float(d) - 0.05 - rng.random() * 0.1
            hi = float(d) + 0.05 + rng.random() * 0.1
            diffs.append({"dataset_id": ds, "variant": variant, "contrast": con, "measure": m, "first": first,
                          "second": second, "first_minus_second": f6(d), "lineage_lo": f6(lo), "lineage_hi": f6(hi),
                          "first_only_right": a, "second_only_right": b, "mcnemar_p": "0.5", "calls_agree": scored - a - b,
                          "scored": scored})

        for variant in ("genes_core", "tool_default"):
            for m in ("ba", "sens", "spec"):
                diff(variant, "program_minus_known_gene_rule", m, "program:" + variant, "grm_known_gene_rule")
        for s in (["random", "lineage"] if lin else ["random"]):
            for l in LEARNERS:
                diff("tool_default", "learner_minus_program:" + s, "ba", "learner:%s:%s" % (l, s),
                     "program:tool_default")
    write_tsv(final / "by-dataset-callers.tsv.gz", CALLER_COLS, callers, gz=True)
    write_tsv(final / "differences.tsv.gz", DIFF_COLS, diffs, gz=True)
    write_tsv(final / "per-dataset-headline.tsv", ["dataset_id", "species_slug", "drug_class", "primary_set",
                                                   "read_before_this_plan", "served_best"], head)
    write_tsv(final / "by-dataset.tsv", ["variant", "dataset_id", "point_screened"], by_dataset)
    (final / "run-manifest.json").write_text(json.dumps({"tool": "AMRFinderPlus", "software_version": "4.2.7",
                                                         "database_version": "2026-08-07.1"}), encoding="utf-8")
    (final / "adapter-summary.json").write_text(json.dumps({"datasets_scored": 253, "genomes_read": 23751}),
                                                encoding="utf-8")

    # The plan's set lists, from the same definitions the script holds them to.
    tables = {
        "datasets.tsv": [{"dataset": h["dataset_id"], "species": h["species_slug"], "drug_class": h["drug_class"],
                          "organism_option": b["point_screened"],
                          "lineage_split": "yes" if r[4] else "no",
                          "read_before_the_plan": h["read_before_this_plan"], "served_learner": h["served_best"]}
                         for h, b, r in zip(head, by_dataset[::2], rows)],
        "callers.tsv": [dict(r, dataset=r["dataset_id"]) for r in callers],
        "differences.tsv": [dict(r, dataset=r["dataset_id"], reading=r["variant"]) for r in diffs
                            if r["variant"] == "tool_default"],
    }
    cohort = record.Cohort(tables)
    lists = top / "results" / "plan-lists"
    lists.mkdir(parents=True, exist_ok=True)
    for name in record.PLAN_LISTS:
        (lists / (name + ".txt")).write_text("".join(x + "\n" for x in cohort.sets[name]), encoding="utf-8")

    # The final scoring's own summary, as its scoring script would give it.
    intervals = {}
    for name, con, m in record.targets():
        seed = SEED_STRING.replace("<set>", name).replace("<variant>", "tool_default").replace("<contrast>", con)
        ends = record.redraw(cohort.by_species(name, con, m), seed.replace("<measure>", m))
        keys = ["median_species_lo", "median_species_hi", "mean_species_lo", "mean_species_hi", "share_species_lo",
                "share_species_hi"]
        intervals[(name, con, m)] = {k: "" if ends is None or (con == record.PR and k.startswith("share")) else f6(e)
                                     for k, e in zip(keys, ends or [0] * 6)}
    summary = record.summarize(cohort, intervals)
    out = []
    for r in summary:
        iv = intervals[(r["set"], r["contrast"], r["measure"])]
        out.append({"set": r["set"], "planned": r["datasets"], "variant": "tool_default", "contrast": r["contrast"],
                    "measure": r["measure"], "units": r["datasets"], "not_on_every_genome": 0,
                    "first_ahead": r["ahead"], "second_ahead": r["behind"], "equal": r["level"],
                    "within_0.02": r["learner_at_most_0.02_ahead"] if r["contrast"] != record.PR else 1,
                    "median": r["median"], "q1": r["q1"], "q3": r["q3"], "mean": r["mean"], "min": "-1.000000",
                    "max": "1.000000", "sign_test_p": r["sign_test_p"], "species": r["species"],
                    "median_species_lo": iv["median_species_lo"], "median_species_hi": iv["median_species_hi"],
                    "mean_species_lo": iv["mean_species_lo"], "mean_species_hi": iv["mean_species_hi"],
                    "share_within_0.02": r["share_learner_at_most_0.02_ahead"] if r["contrast"] != record.PR else "0.5",
                    "share_within_species_lo": iv["share_species_lo"] or ("0.1" if iv["median_species_lo"] else ""),
                    "share_within_species_hi": iv["share_species_hi"] or ("0.9" if iv["median_species_lo"] else ""),
                    "lineage_interval_above_zero": r["own_interval_above_zero"],
                    "lineage_interval_below_zero": r["own_interval_below_zero"],
                    "lineage_interval_spans_zero": r["own_interval_spans_zero"]})
    write_tsv(final / "cohort-summary.tsv.gz", COHORT_COLS, out, gz=True)

    # The plan, its addenda, their hash files and the logs.
    logs = top / "results" / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    texts = [plan_text()] + ["# A made-up addendum %d\n\nIt changes no comparison.\n" % i for i in (1, 2, 3)]
    for (event, name, sums, when), text, utc in zip(record.ORDER, texts, TIMES):
        (top / name).write_text(text, encoding="utf-8")
        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        (top / sums).write_text("%s  %s\n%s  scripts/made_up.py\n" % (digest, name, "0" * 64), encoding="utf-8")
        (top / when).write_text(utc + "\n", encoding="utf-8")
    (top / record.FINAL_LOG).write_text("%s scope final: %d datasets\nmore lines\n" % (FINAL_START, len(rows)),
                                        encoding="utf-8")
    return final
