#!/usr/bin/env python3
"""Step one of PLAN.md: counts by study, with no fit.

Reads the six released clindamycin reports (samples: id, resistant, catalog_resistant, lineage), the curated
phenotype table and the raw BV-BRC AMR records of the public data root, and writes under results/:

  by-study.tsv            one row per dataset and study key (the rows of a dataset add up to it)
  by-pmid.tsv             one row per dataset and single PMID (a sample with two PMIDs counts under both)
  chosen.tsv              the chosen study T(d) per dataset, its shares and intervals, the primary pair
  membership/<d>.tsv      one row per sample: label, rule call, study key, source field, lineage, in_chosen
  removal/<d>.ids         the genome ids that belong to T(d), one per line, sorted
  count-checks.txt        the checks of PLAN.md section 2 and what the raw field held

Run: python3 -I scripts/count_by_study.py
"""
import collections
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from stats_util import clopper_pearson, fisher_two_sided, odds_ratio  # noqa: E402

RELEASE = "/home/user/grm-ops/owner2/work/rev14-full/rev14-274-20261112"
ROOT = "/home/user/grm-h1raw"
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")
DRUG = "clindamycin"
DATASETS = [
    "campylobacter_coli__clindamycin", "campylobacter_jejuni__clindamycin", "enterococcus_faecium__clindamycin",
    "staphylococcus_aureus__clindamycin", "streptococcus_agalactiae__clindamycin",
    "streptococcus_pneumoniae__clindamycin",
]


def parse_call(text):
    """src/curate.rs parse_call."""
    if text is None:
        return None
    t = text.strip().lower()
    if t in ("resistant", "nonsusceptible", "non-susceptible"):
        return "R"
    if t == "susceptible":
        return "S"
    if t in ("intermediate", "susceptible-dose dependent", "reduced susceptibility"):
        return "I"
    return None


def parts(value):
    """src/curate.rs parts: a value fetch joined from a list (a;b) counts as each part."""
    if not value:
        return []
    return [p.strip() for p in str(value).split(";") if p.strip()]


def pmid_order(p):
    return (0, int(p), p) if p.isdigit() else (1, 0, p)


def fmt(x, nd=4):
    if x is None or x != x:
        return "NA"
    if x in (float("inf"), float("-inf")):
        return "inf"
    return f"{x:.{nd}f}"


def main():
    os.makedirs(os.path.join(OUT, "membership"), exist_ok=True)
    os.makedirs(os.path.join(OUT, "removal"), exist_ok=True)
    checks = []

    # The index must list exactly these six clindamycin datasets.
    index = json.load(open(os.path.join(RELEASE, "index.json")))
    in_index = sorted(d["id"] for d in index["datasets"] if d["antibiotic"] == DRUG)
    assert in_index == sorted(DATASETS), in_index
    checks.append(f"index: evaluation_revision={index['evaluation_revision']} datasets={len(index['datasets'])} "
                  f"clindamycin datasets={len(in_index)}")

    reports = {d: json.load(open(os.path.join(RELEASE, "datasets", d + ".json"))) for d in DATASETS}
    wanted = {}  # genome id -> list of datasets
    for d, rep in reports.items():
        ids = [s["id"] for s in rep["samples"]]
        assert len(ids) == len(set(ids)), d
        for g in ids:
            wanted.setdefault(g, []).append(d)

    # Raw names that the curated name table resolves to the drug.
    names = json.load(open(os.path.join(ROOT, "curated", "drug_names.json")))["names"]
    raw_names = sorted(n["raw"] for n in names if n["canonical"] == DRUG and n["kind"] == "drug")
    other_kind = sorted((n["raw"], n["kind"], n["canonical"]) for n in names
                        if DRUG in n["raw"].lower() and not (n["canonical"] == DRUG and n["kind"] == "drug"))
    checks.append(f"drug_names.json: raw names resolved to drug '{DRUG}': {raw_names}")
    checks.append(f"drug_names.json: other raw names that hold '{DRUG}' (not used): {other_kind}")
    raw_set = set(raw_names)

    # Curated rows of the drug for the wanted genomes.
    curated = collections.defaultdict(list)
    with open(os.path.join(ROOT, "curated", "phenotypes.jsonl")) as fh:
        for line in fh:
            if DRUG not in line:
                continue
            r = json.loads(line)
            if r["antibiotic"] == DRUG and r["genome_id"] in wanted:
                curated[r["genome_id"]].append(r)

    # Raw records of the drug for the wanted genomes.
    raw = collections.defaultdict(list)
    n_raw = 0
    with open(os.path.join(ROOT, "bvbrc", "genome_amr.jsonl")) as fh:
        for line in fh:
            n_raw += 1
            # Prefilter only; the exact test is the raw name set below. Lower-cased so that a raw name in
            # capitals is not skipped (fix of 15 Nov: the 14 Nov text tested "lindamycin" case-sensitively).
            if "lindamycin" not in line.lower():
                continue
            r = json.loads(line)
            if r["antibiotic"] in raw_set and r["genome_id"] in wanted:
                raw[r["genome_id"]].append(r)
    checks.append(f"genome_amr.jsonl: {n_raw} records read; clindamycin records of the six datasets' samples: "
                  f"{sum(len(v) for v in raw.values())} over {len(raw)} genomes")

    study_rows, pmid_rows, chosen_rows = [], [], []
    for d in DATASETS:
        rep = reports[d]
        samples = rep["samples"]
        n = len(samples)
        mism = 0
        member = []
        source_differs = 0
        records_used = 0
        no_label_record = 0
        for s in samples:
            g = s["id"]
            label = "R" if s["resistant"] else "S"
            rule = "R" if s["catalog_resistant"] else "S"
            rows = curated.get(g, [])
            ok = (len(rows) == 1 and rows[0]["reported"] == label and rows[0]["conflict"] is False
                  and rows[0]["label_from_mic"] is False)
            if not ok:
                mism += 1
            label_records = [r for r in raw.get(g, []) if parse_call(r.get("resistant_phenotype")) == label]
            if not label_records:
                no_label_record += 1
            records_used += len(label_records)
            pmids = sorted({p for r in label_records for p in parts(r.get("pmid"))}, key=pmid_order)
            sources = sorted({p for r in label_records for p in parts(r.get("source"))})
            if sources and set(sources) != set(pmids):
                source_differs += 1
            key = "+".join(pmids) if pmids else "no-pmid"
            member.append({"id": g, "label": label, "rule": rule, "key": key, "pmids": pmids,
                           "sources": sources, "lineage": s["lineage"], "n_label_records": len(label_records)})
        checks.append(f"{d}: samples={n}; curated-row mismatches={mism}; samples with no raw record carrying the "
                      f"label={no_label_record}; label records used={records_used}; samples whose source field "
                      f"differs from pmid={source_differs}")
        if mism:
            checks.append(f"{d}: STOPPED by PLAN.md section 2 (curated rows do not match the report)")
            continue

        def tally(ms):
            c = collections.Counter()
            for m in ms:
                c["n"] += 1
                c["rep_" + m["label"]] += 1
                c["rule_" + m["rule"]] += 1
                if m["rule"] == "R" and m["label"] == "S":
                    c["fr"] += 1
                if m["rule"] == "S" and m["label"] == "R":
                    c["fs"] += 1
            c["dis"] = c["fr"] + c["fs"]
            return c

        total = tally(member)
        by_key = collections.defaultdict(list)
        by_pmid = collections.defaultdict(list)
        for m in member:
            by_key[m["key"]].append(m)
            for p in m["pmids"]:
                by_pmid[p].append(m)

        def row(name, c):
            return [d, name, c["n"], c["rep_R"], c["rep_S"], c["rule_R"], c["rule_S"], c["fr"], c["fs"], c["dis"],
                    fmt(c["dis"] / c["n"] if c["n"] else None),
                    fmt(c["dis"] / total["dis"] if total["dis"] else None),
                    fmt(c["n"] / total["n"])]

        keys = sorted(by_key, key=lambda k: (-tally(by_key[k])["dis"], -len(by_key[k]), k))
        for k in keys:
            study_rows.append(row(k, tally(by_key[k])))
        study_rows.append(row("TOTAL", total))

        pm = {p: tally(ms) for p, ms in by_pmid.items()}
        order = sorted(pm, key=lambda p: (-pm[p]["dis"], -pm[p]["fr"], -pm[p]["n"], pmid_order(p)))
        for p in order:
            pmid_rows.append(row(p, pm[p]))

        # The choice, PLAN.md section 3.
        eligible = [p for p in order if pm[p]["dis"] > 0]
        if not eligible:
            chosen_rows.append({"dataset": d, "chosen": None, "total": total})
            checks.append(f"{d}: no PMID has a disagreement; no chosen study")
            continue
        T = eligible[0]
        most_fr = sorted(pm, key=lambda p: (-pm[p]["fr"], -pm[p]["dis"], -pm[p]["n"], pmid_order(p)))[0]
        c = pm[T]
        rest = collections.Counter({k: total[k] - c[k] for k in total})
        a, b = c["dis"], c["n"] - c["dis"]
        cc, dd = rest["dis"], rest["n"] - rest["dis"]
        chosen_rows.append({
            "dataset": d, "chosen": T, "total": total, "study": c, "rest": rest,
            "most_fr": most_fr, "most_fr_count": pm[most_fr]["fr"],
            "share_dis": (c["dis"], total["dis"], clopper_pearson(c["dis"], total["dis"])),
            "share_fr": (c["fr"], total["fr"], clopper_pearson(c["fr"], total["fr"]) if total["fr"] else None),
            "share_fs": (c["fs"], total["fs"], clopper_pearson(c["fs"], total["fs"]) if total["fs"] else None),
            "share_n": (c["n"], total["n"]),
            "rate_in": (a, a + b), "rate_out": (cc, cc + dd),
            "or": odds_ratio(a, b, cc, dd), "fisher_p": fisher_two_sided(a, b, cc, dd),
            "n_pmids": len(pm), "runner_up": (eligible[1], pm[eligible[1]]["dis"]) if len(eligible) > 1 else None,
            "left_R": total["rep_R"] - c["rep_R"], "left_S": total["rep_S"] - c["rep_S"],
            "min_class": rep["counts"]["min_class"],
        })

        with open(os.path.join(OUT, "membership", d + ".tsv"), "w") as fh:
            fh.write("genome_id\treported\tknown_gene_rule\tstudy_key\tsource_field\tlineage\tn_label_records\t"
                     "in_chosen_study\n")
            for m in member:
                fh.write("\t".join([m["id"], m["label"], m["rule"], m["key"], "+".join(m["sources"]) or "none",
                                    str(m["lineage"]), str(m["n_label_records"]),
                                    "1" if T in m["pmids"] else "0"]) + "\n")
        with open(os.path.join(OUT, "removal", d + ".ids"), "w") as fh:
            for g in sorted(m["id"] for m in member if T in m["pmids"]):
                fh.write(g + "\n")

    head = ["dataset", "study", "samples", "reported_R", "reported_S", "rule_R", "rule_S", "false_resistant",
            "false_susceptible", "disagreements", "disagreement_rate", "share_of_dataset_disagreements",
            "share_of_dataset_samples"]
    with open(os.path.join(OUT, "by-study.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in study_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    with open(os.path.join(OUT, "by-pmid.tsv"), "w") as fh:
        fh.write("\t".join(head) + "\n")
        for r in pmid_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # The primary pair: largest D(d, T(d)); ties: larger share, then dataset id.
    with_choice = [r for r in chosen_rows if r["chosen"]]
    primary = sorted(with_choice, key=lambda r: (-r["share_dis"][0], -(r["share_dis"][0] / r["share_dis"][1]),
                                                 r["dataset"]))[0]["dataset"] if with_choice else None

    def word(share):
        return "most" if share >= 0.5 else ("a large share" if share >= 0.25 else "not a large share")

    with open(os.path.join(OUT, "chosen.tsv"), "w") as fh:
        fh.write("\t".join([
            "dataset", "primary_pair", "chosen_pmid", "pmids_in_dataset", "samples_in_study", "samples_in_dataset",
            "share_of_samples", "disagreements_in_study", "disagreements_in_dataset", "share_of_disagreements",
            "share_ci95_low", "share_ci95_high", "plan_word", "false_resistant_in_study",
            "false_resistant_in_dataset", "fr_share", "fr_ci95_low", "fr_ci95_high", "false_susceptible_in_study",
            "false_susceptible_in_dataset", "fs_share", "fs_ci95_low", "fs_ci95_high", "rate_in_study",
            "rate_in_rest", "odds_ratio", "fisher_two_sided_p", "runner_up_pmid", "runner_up_disagreements",
            "pmid_with_most_false_resistant", "its_false_resistant", "reported_R_left", "reported_S_left",
            "min_class", "refit"]) + "\n")
        for r in chosen_rows:
            d = r["dataset"]
            if not r["chosen"]:
                fh.write("\t".join([d, "0", "none"] + ["NA"] * 32) + "\n")
                continue
            k, n_, ci = r["share_dis"]

            def tri(t):
                if t[2] is None:
                    return [str(t[0]), str(t[1]), "NA", "NA", "NA"]
                return [str(t[0]), str(t[1]), fmt(t[0] / t[1]), fmt(t[2][0]), fmt(t[2][1])]
            refit = "yes" if min(r["left_R"], r["left_S"]) >= r["min_class"] else "no: a class falls under min_class"
            fh.write("\t".join([
                d, "1" if d == primary else "0", r["chosen"], str(r["n_pmids"]), str(r["share_n"][0]),
                str(r["share_n"][1]), fmt(r["share_n"][0] / r["share_n"][1]), str(k), str(n_), fmt(k / n_),
                fmt(ci[0]), fmt(ci[1]), word(k / n_)] + tri(r["share_fr"]) + tri(r["share_fs"]) + [
                f"{r['rate_in'][0]}/{r['rate_in'][1]}={fmt(r['rate_in'][0] / r['rate_in'][1])}",
                f"{r['rate_out'][0]}/{r['rate_out'][1]}="
                + (fmt(r['rate_out'][0] / r['rate_out'][1]) if r['rate_out'][1] else "NA"),
                fmt(r["or"], 3), f"{r['fisher_p']:.3g}",
                r["runner_up"][0] if r["runner_up"] else "none",
                str(r["runner_up"][1]) if r["runner_up"] else "NA",
                r["most_fr"], str(r["most_fr_count"]), str(r["left_R"]), str(r["left_S"]), str(r["min_class"]),
                refit]) + "\n")
    checks.append(f"primary pair: {primary}")
    with open(os.path.join(OUT, "count-checks.txt"), "w") as fh:
        fh.write("\n".join(checks) + "\n")
    print("\n".join(checks))


if __name__ == "__main__":
    main()
