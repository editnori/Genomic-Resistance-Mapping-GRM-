"""tb_curve.py: block M4 of the mechanism split, the lineage curve for
Mycobacterium tuberculosis.

    python3 -I tb_curve.py --species SLUG --reports DIR --lists lists.tsv \
        --tool-dir RUN_DIR --bed TBDB.bed --out DIR

Reads TBProfiler 6.6.3 output (database `tbdb` commit 3ed1cb99, WHO catalogue
v2) run again on the species' genomes, beside the release reports for the
labels and the lineage clusters. A determinant of a drug is a variant the
tool reports for that drug, passing its own filter, with a resistance entry
graded "Assoc w R" or "Assoc w R - Interim" (the `who_assoc` reading of the
tbwho lane's adapter). The report antibiotic is matched to the tool's drug by
name, with rifampin = rifampicin and the four proxy stems, as that adapter
matches it (`tb_adapter.py` of the TBProfiler comparison, sha256
ff8054be2e43b830b48706a64879de1fd0bc967b7e5e19655fa5be2318c20988; its
`drug_map`, `norm`, grades and filter are copied here unchanged).

Per antibiotic in full, over the resistant genomes of the release label:
  - the accumulation curve of distinct determinants as lineage groups are
    added, in 1,000 seeded random orders (a median and a 2.5-to-97.5 band)
    and in the order by group size, largest first;
  - the share of resistant genomes with no graded determinant, per group and
    over all;
  - the convergence, how many determinants are seen in 1, 2, 3 and 4 or more
    lineages, and the share of resistant genomes they cover;
  - per gene the WHO catalogue names, a count of its distinct graded variants.

Groups: the tool's main lineage (strict), its sublineage (lenient), and the
report's own lineage cluster beside. A genome the tool gives no lineage
enters no group of that grouping and is counted in the summary; a mixed call
is its own group, as the tool writes it.

What leaves the box: counts and shares with species, antibiotic and gene
names only. A variant's own name goes into no output (lane question 9). A
cell under 5 genomes is written "under 5". No genome id, no row per genome,
no sequence, no position.

Python 3.12, standard library only; run with `python3 -I`.
"""

import argparse
import collections
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import FLOOR, r12, read_tsv, seed_of, shown, write_tsv  # noqa: E402

WANT_TOOL = "6.6.3"
WANT_DB = "3ed1cb99"
ASSOC, INTERIM = "Assoc w R", "Assoc w R - Interim"
COUNTED = (ASSOC, INTERIM)
SYNONYMS = {"rifampin": "rifampicin"}
PROXIES = {  # report antibiotic -> tool drugs, by a shared name stem
    "ciprofloxacin": ("levofloxacin", "moxifloxacin"),
    "ofloxacin": ("levofloxacin", "moxifloxacin"),
    "rifabutin": ("rifampicin",),
    "prothionamide": ("ethionamide",),
}
ADAPTER_SHA256 = "ff8054be2e43b830b48706a64879de1fd0bc967b7e5e19655fa5be2318c20988"
ORDERS = 1000
GROUPINGS = ("main_lineage", "sub_lineage", "cluster")
CONVERGENCE_GROUPINGS = ("main_lineage", "sub_lineage")
BUCKETS = ("1", "2", "3", "4+")
FOLDED_GENES = "genes seen in under 5 genomes"


class Refusal(Exception):
    """An input breaks the contract. Nothing is written after one."""


def refuse(message):
    raise Refusal(message)


def norm(name):
    return name.lower().replace("-", " ").replace("_", " ").strip()


def drug_map(antibiotic, tool_drugs):
    """(kind, tool drugs) for a report antibiotic, from names alone."""
    by_norm = {norm(d): d for d in tool_drugs}
    key = norm(SYNONYMS.get(norm(antibiotic), antibiotic))
    if key in by_norm:
        return "primary", (by_norm[key],)
    if norm(antibiotic) in PROXIES:
        drugs = PROXIES[norm(antibiotic)]
        missing = [d for d in drugs if d not in tool_drugs]
        if missing:
            refuse(f"proxy for {antibiotic} names {missing}, which the tool does not have")
        return "proxy", drugs
    return "none", ()


def read_bed(path):
    """The tool's drug names, from the sixth column of tbdb.bed."""
    drugs = set()
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, start=1):
            row = line.rstrip("\n").split("\t")
            if len(row) < 6:
                refuse(f"{path}:{n}: fewer than 6 columns")
            for drug in row[5].split(","):
                if drug:
                    drugs.add(drug)
    if not drugs:
        refuse(f"{path}: names no drug")
    return sorted(drugs)


def genome_determinants(path, gid):
    """(main lineage, sublineage, {tool drug: {(gene, change)}}) for one result
    file. The change is read for a determinant's identity and never written."""
    with open(path, encoding="utf-8") as f:
        d = json.load(f)
    if d["id"] != gid:
        refuse(f"{os.path.basename(path)}: holds the result of another genome")
    sw = d["pipeline"]["software_version"]
    db = d["pipeline"]["db_version"]["commit"]
    if sw != WANT_TOOL or db != WANT_DB:
        refuse(f"{os.path.basename(path)}: tool {sw}, database {db}; wanted {WANT_TOOL} and {WANT_DB}")
    per_drug = collections.defaultdict(set)
    for v in d["dr_variants"]:
        if v["filter"] != "pass":
            continue
        key = (v["gene_name"], v["change"])
        for a in v["drugs"]:
            if a.get("type") != "drug_resistance":
                continue
            if a["confidence"] in COUNTED:
                per_drug[a["drug"]].add(key)
    main = d.get("main_lineage") or None
    sub = d.get("sub_lineage") or None
    return main, sub, dict(per_drug)


def percentile(values, q):
    """Nearest-rank percentile of a sorted list, q in 0..1."""
    i = int(q * (len(values) - 1) + 0.5)
    return values[min(len(values) - 1, max(0, i))]


def curve_rows(species, antibiotic, in_m, grouping, groups):
    """`groups` maps a group key to the set of determinants its resistant
    genomes carry. The accumulation curve over 1,000 random orders and the
    largest-first order."""
    keys = sorted(groups, key=lambda k: (-groups[k]["genomes"], k))
    n = len(keys)
    rng = random.Random(seed_of(f"{antibiotic}|{grouping}"))
    steps = [[] for _ in range(n)]
    order = list(keys)
    for _ in range(ORDERS):
        rng.shuffle(order)
        seen = set()
        for step, k in enumerate(order):
            seen |= groups[k]["dets"]
            steps[step].append(len(seen))
    largest = []
    seen = set()
    for k in keys:
        seen |= groups[k]["dets"]
        largest.append(len(seen))
    resistant = sum(g["genomes"] for g in groups.values())
    if resistant < FLOOR:
        # a curve over 1 to 4 genomes would stand for those genomes: none is written
        return []
    rows = []
    for step in range(n):
        vals = sorted(steps[step])
        mid = vals[len(vals) // 2] if len(vals) % 2 else (vals[len(vals) // 2 - 1] + vals[len(vals) // 2]) / 2
        rows.append({"species": species, "antibiotic": antibiotic, "in_M": "1" if in_m else "0",
                     "grouping": grouping, "groups": n, "resistant": shown(resistant), "step": step + 1,
                     "determinants_median": r12(mid), "determinants_lo": percentile(vals, 0.025),
                     "determinants_hi": percentile(vals, 0.975), "determinants_largest_first": largest[step]})
    return rows


def coverage_rows(species, antibiotic, in_m, grouping, groups):
    rows = []
    for k in sorted(groups):
        g = groups[k]
        res = g["genomes"]
        without = g["without"]
        rows.append({"species": species, "antibiotic": antibiotic, "in_M": "1" if in_m else "0",
                     "grouping": grouping, "group": k, "resistant": shown(res),
                     "without_determinant": shown(without),
                     "share_without": f"{r12(without / res):.12f}" if res >= FLOOR else ""})
    return rows


def convergence_rows(species, antibiotic, in_m, grouping, per_genome):
    """`per_genome` is a list of (group key, determinant set) for the resistant
    genomes with a group under this grouping. Determinants by the number of
    lineages they are seen in, and the resistant genomes each bucket covers."""
    lineages_of = collections.defaultdict(set)
    for key, dets in per_genome:
        for det in dets:
            lineages_of[det].add(key)
    bucket_of = {}
    for det, lins in lineages_of.items():
        c = len(lins)
        bucket_of[det] = "4+" if c >= 4 else str(c)
    covered = {b: 0 for b in BUCKETS}
    determinants = {b: 0 for b in BUCKETS}
    for b in BUCKETS:
        determinants[b] = sum(1 for det, bb in bucket_of.items() if bb == b)
    for key, dets in per_genome:
        for b in BUCKETS:
            if any(bucket_of.get(det) == b for det in dets):
                covered[b] += 1
    total = len(per_genome)
    rows = []
    for b in BUCKETS:
        rows.append({"species": species, "antibiotic": antibiotic, "in_M": "1" if in_m else "0",
                     "grouping": grouping, "lineages_seen": b, "determinants": determinants[b],
                     "resistant_covered": shown(covered[b]),
                     "share_covered": f"{r12(covered[b] / total):.12f}" if total >= FLOOR else ""})
    return rows


def gene_rows(species, antibiotic, in_m, dets_by_gene, gene_sets):
    """Per gene, its graded variants seen and the resistant genomes with one.
    `gene_sets` holds each resistant genome's genes. A gene seen in 1 to 4
    genomes is not named: it joins one row under FOLDED_GENES."""
    base = {"species": species, "antibiotic": antibiotic, "in_M": "1" if in_m else "0"}
    genomes_by_gene = collections.Counter(g for genes in gene_sets for g in genes)
    named = [g for g in sorted(dets_by_gene) if genomes_by_gene[g] >= FLOOR]
    folded = {g for g in dets_by_gene if genomes_by_gene[g] < FLOOR}
    rows = [dict(base, gene=g, catalogue_variants=len(dets_by_gene[g]),
                 resistant_with_gene=shown(genomes_by_gene[g])) for g in named]
    if folded:
        rows.append(dict(base, gene=FOLDED_GENES, catalogue_variants=sum(len(dets_by_gene[g]) for g in folded),
                         resistant_with_gene=shown(sum(1 for genes in gene_sets if genes & folded))))
    return rows


def load_reports(directory, species):
    reports = []
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".json"):
            continue
        with open(os.path.join(directory, name), encoding="utf-8") as f:
            rep = json.load(f)
        if rep.get("species_slug") != species:
            continue
        reports.append(rep)
    if not reports:
        refuse(f"{species}: no report in {directory}")
    return reports


def read_plan(tool_dir):
    with open(os.path.join(tool_dir, "plan.txt"), encoding="utf-8") as f:
        plan = [line.strip() for line in f if line.strip()]
    failed = set()
    fpath = os.path.join(tool_dir, "failed.tsv")
    if os.path.isfile(fpath):
        with open(fpath, encoding="utf-8") as f:
            failed = {line.rstrip("\n").split("\t")[0] for line in f if line.strip()}
    return plan, failed


def run(args):
    tool_drugs = read_bed(args.bed)
    reports = load_reports(args.reports, args.species)
    species = reports[0].get("species", args.species)
    in_m = {}
    for row in read_tsv(args.lists):
        if row["species_slug"] == args.species:
            in_m[row["dataset"]] = row.get("in_M") == "1"

    plan, failed = read_plan(args.tool_dir)
    evidence, lineage = {}, {}
    for gid in plan:
        rpath = os.path.join(args.tool_dir, "results", gid + ".results.json")
        if os.path.isfile(rpath):
            if gid in failed:
                refuse("a genome is listed as failed and has a result file")
            main, sub, per_drug = genome_determinants(rpath, gid)
            evidence[gid] = per_drug
            lineage[gid] = {"main_lineage": main, "sub_lineage": sub}
        elif gid not in failed:
            refuse("a genome has no result file and is not listed as failed")

    curve, coverage, convergence, genes = [], [], [], []
    summary = {"tool": WANT_TOOL, "database_commit": WANT_DB, "orders": ORDERS, "species": args.species,
               "tool_drugs": tool_drugs, "genomes_planned": shown(len(plan)),
               "genomes_with_result": shown(len(evidence)), "genomes_failed": shown(len(failed)), "endpoints": [],
               "not_scored_reports": []}
    for rep in reports:
        antibiotic = rep["antibiotic"]
        kind, drugs = drug_map(antibiotic, tool_drugs)
        if kind == "none":
            summary["not_scored_reports"].append({"antibiotic": antibiotic, "reason": "the tool has no drug of this name"})
            continue
        is_m = in_m.get(rep["id"], False)
        # Resistant genomes of the release label with a result, their determinants and lineages.
        members, no_result = [], 0
        for s in rep["samples"]:
            if not s.get("resistant"):
                continue
            gid = s["id"]
            if gid in failed or gid not in evidence:
                no_result += 1
                continue
            dets = set()
            for drug in drugs:
                dets |= evidence[gid].get(drug, set())
            members.append({"main_lineage": lineage[gid]["main_lineage"], "sub_lineage": lineage[gid]["sub_lineage"],
                            "cluster": None if s.get("lineage") is None else f"cluster_{s['lineage']}", "dets": dets})

        dets_by_gene = collections.defaultdict(set)
        gene_sets = []
        for m in members:
            for gene, change in m["dets"]:
                dets_by_gene[gene].add((gene, change))
            gene_sets.append({gene for gene, _ in m["dets"]})
        genes += gene_rows(species, antibiotic, is_m, dets_by_gene, gene_sets)

        for grouping in GROUPINGS:
            groups = collections.defaultdict(lambda: {"dets": set(), "genomes": 0, "without": 0})
            for m in members:
                key = m[grouping]
                if key is None:
                    continue
                g = groups[key]
                g["dets"] |= m["dets"]
                g["genomes"] += 1
                if not m["dets"]:
                    g["without"] += 1
            if groups:
                curve += curve_rows(species, antibiotic, is_m, grouping, groups)
                coverage += coverage_rows(species, antibiotic, is_m, grouping, groups)
        # The overall no-determinant share, over every resistant genome with a result.
        total = len(members)
        without = sum(1 for m in members if not m["dets"])
        coverage.append({"species": species, "antibiotic": antibiotic, "in_M": "1" if is_m else "0",
                         "grouping": "all", "group": "all", "resistant": shown(total),
                         "without_determinant": shown(without),
                         "share_without": f"{r12(without / total):.12f}" if total >= FLOOR else ""})
        for grouping in CONVERGENCE_GROUPINGS:
            per_genome = [(m[grouping], m["dets"]) for m in members if m[grouping] is not None]
            if per_genome:
                convergence += convergence_rows(species, antibiotic, is_m, grouping, per_genome)
        summary["endpoints"].append({"antibiotic": antibiotic, "kind": kind, "tool_drugs": list(drugs),
                                     "in_M": is_m, "resistant_with_result": shown(total),
                                     "resistant_no_result": shown(no_result),
                                     "resistant_with_no_group": {gr: shown(sum(1 for m in members if m[gr] is None))
                                                                 for gr in GROUPINGS}})

    os.makedirs(args.out, exist_ok=True)
    write_tsv(os.path.join(args.out, "m4-curve.tsv"),
              ["species", "antibiotic", "in_M", "grouping", "groups", "resistant", "step",
               "determinants_median", "determinants_lo", "determinants_hi", "determinants_largest_first"], curve)
    write_tsv(os.path.join(args.out, "m4-coverage.tsv"),
              ["species", "antibiotic", "in_M", "grouping", "group", "resistant", "without_determinant",
               "share_without"], coverage)
    write_tsv(os.path.join(args.out, "m4-convergence.tsv"),
              ["species", "antibiotic", "in_M", "grouping", "lineages_seen", "determinants", "resistant_covered",
               "share_covered"], convergence)
    write_tsv(os.path.join(args.out, "m4-genes.tsv"),
              ["species", "antibiotic", "in_M", "gene", "catalogue_variants", "resistant_with_gene"], genes)
    with open(os.path.join(args.out, "m4-summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, sort_keys=True)
        f.write("\n")
    print(f"{len(summary['endpoints'])} endpoints, {len(evidence)} genomes with a result; "
          f"{len(curve)} curve rows, {len(genes)} gene rows")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("--species", "--reports", "--lists", "--tool-dir", "--bed", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args(argv)
    try:
        return run(args)
    except Refusal as e:
        print(f"tb_curve.py: refused: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
