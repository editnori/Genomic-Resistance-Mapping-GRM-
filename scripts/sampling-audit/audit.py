#!/usr/bin/env python3
"""Read-only audit of the legacy genome selection, per task (species x
antibiotic) and per cap.

For every dataset in datasets.json this reports the eligible rows (genomes
with an unconflicted, reported R or S call for the drug, the only rows a
model trains and scores on), the rows a cap selects through the ported
`select_genomes`, their class counts and prevalence, each genome's panel
breadth, and, when group ids are given, how many duplicate or lineage groups
those rows span. Consecutive caps are nested (the ranking does not depend on
the cap), so it also reports what each larger cap adds. Per species it
reports where each cap cuts the ranking and whether genome-id tie breaks
decide who is in. Every figure is descriptive; README.md lists the caveats.

It reads phenotypes.jsonl and datasets.json as `grm curate` writes them and
writes only the output paths it is given.

usage: python3 -I audit.py OUT.json --curated DIR [--caps 2000,4000]
           [--species NAME ...] [--groups NAME=TSV ...] [--tsv OUT.tsv]
--curated is a curate output directory holding phenotypes.jsonl and
datasets.json (or give --phenotypes and --datasets). --caps are per-species
caps; "all" means every scored genome. --groups adds one group set, a TSV
with a genome_id and a group_id column (a header line starting with
genome_id is skipped); repeat it for duplicate and lineage groups.
"""
import argparse
import json
import os
import statistics
import sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import legacy_select as ls  # noqa: E402

SCHEMA = 2
ALL = "all"
# Group ESS is withheld when fewer of the rows than this have a group id.
MIN_GROUP_COVERAGE = Fraction(4, 5)


def parse_caps(text):
    caps = []
    for part in text.split(","):
        part = part.strip()
        caps.append(ALL if part == ALL else int(part))
    if any(c != ALL and c < 0 for c in caps):
        raise ValueError("caps must be non-negative")
    if len(set(caps)) != len(caps):
        raise ValueError("caps repeat")
    # "all" is the largest cap.
    return sorted(caps, key=lambda c: (c == ALL, 0 if c == ALL else c))


def read_groups(path):
    """{genome_id: group_id} from a two-column TSV."""
    groups = {}
    with open(path, encoding="utf-8") as f:
        for n, line in enumerate(f, 1):
            line = line.rstrip("\r\n")
            if not line or (n == 1 and line.split("\t", 1)[0] == "genome_id"):
                continue
            fields = line.split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError(f"{path}:{n}: expected genome_id<TAB>group_id")
            gid, group = fields[0], fields[1]
            if groups.get(gid, group) != group:
                raise ValueError(f"{path}:{n}: {gid} is in two groups")
            groups[gid] = group
    return groups


def distribution(values):
    """Histogram and summary of a list of integers; None when empty."""
    if not values:
        return None
    hist = {}
    for v in values:
        hist[v] = hist.get(v, 0) + 1
    return {
        "n": len(values),
        "min": min(values),
        "median": statistics.median(values),
        "mean": sum(values) / len(values),
        "max": max(values),
        "histogram": {str(k): hist[k] for k in sorted(hist)},
    }


def inverse_simpson(sizes):
    """(sum n)^2 / sum n^2: the number of equal groups with the same concentration."""
    total = sum(sizes)
    if total == 0:
        return None
    return total * total / sum(n * n for n in sizes)


def coverage_gate(grouped, total):
    """(share of rows with a group id, reason the ESS is withheld or None)."""
    if total == 0:
        return None, "no rows"
    if Fraction(grouped, total) < MIN_GROUP_COVERAGE:
        return grouped / total, f"{grouped} of {total} rows grouped, under {float(MIN_GROUP_COVERAGE):.0%}"
    return grouped / total, None


def group_summary(rows, minority, groups):
    """Group metrics over task rows [(genome_id, call)]. A row whose genome has
    no group id is counted apart and left out of every group figure, so "no
    group data" never reads as "all singletons"."""
    sizes, minority_sizes = {}, {}
    missing = minority_missing = n_minority = 0
    for gid, call in rows:
        n_minority += call == minority
        group = groups.get(gid)
        if group is None:
            missing += 1
            minority_missing += call == minority
            continue
        sizes[group] = sizes.get(group, 0) + 1
        if call == minority:
            minority_sizes[group] = minority_sizes.get(group, 0) + 1
    grouped = len(rows) - missing
    coverage, reason = coverage_gate(grouped, len(rows))
    minority_coverage, minority_reason = coverage_gate(n_minority - minority_missing, n_minority)
    largest = max(sizes.values(), default=0)
    return {
        "rows_with_group": grouped,
        "rows_without_group": missing,
        "group_coverage": coverage,
        "groups": len(sizes),
        "minority_groups": len(minority_sizes),
        "largest_group_rows": largest,
        "largest_group_share": largest / grouped if grouped else None,
        "group_ess": None if reason else inverse_simpson(list(sizes.values())),
        "group_ess_null_reason": reason,
        "minority_rows_without_group": minority_missing,
        "minority_group_coverage": minority_coverage,
        "minority_group_ess": None if minority_reason else inverse_simpson(list(minority_sizes.values())),
        "minority_group_ess_null_reason": minority_reason,
    }, set(sizes), set(minority_sizes)


def pair_reason(small, small_reason, large, large_reason):
    """Why a figure comparing two caps is withheld, or None when both pass."""
    reasons = [f"cap {cap}: {r}" for cap, r in ((small, small_reason), (large, large_reason)) if r]
    return "; ".join(reasons) or None


def rows_summary(rows, minority, breadth, share, group_sets):
    """Counts, prevalence, panel breadth, score share and group metrics of
    task rows [(genome_id, call)]. Returns (summary, {set: (groups, minority_groups)})."""
    n_r = sum(1 for _, c in rows if c == ls.R)
    n_s = len(rows) - n_r
    out = {
        "rows": len(rows),
        "R": n_r,
        "S": n_s,
        "prevalence_R": n_r / len(rows) if rows else None,
        "minority_rows": n_r if minority == ls.R else n_s,
        "panel_breadth": distribution([breadth[gid] for gid, _ in rows]),
        "task_score_share_mean": (sum(share[gid] for gid, _ in rows) / len(rows)) if rows else None,
    }
    members = {}
    if group_sets:
        out["groups"] = {}
        for name, groups in group_sets.items():
            out["groups"][name], g, mg = group_summary(rows, minority, groups)
            members[name] = (g, mg)
    return out, members


def cap_slice(ranking, cap):
    return ranking if cap == ALL else ranking[:cap]


def cutoff(ranking, cap):
    """Where a cap cuts the ranking, and whether ids decide a boundary tie."""
    if cap == ALL or len(ranking) <= cap:
        return {"cap_binds": False, "selected": len(ranking)}
    if cap == 0:
        return {"cap_binds": True, "selected": 0}
    edge = ranking[cap - 1][0]
    tied_total = sum(1 for s, _ in ranking if s == edge)
    tied_selected = sum(1 for s, _ in ranking[:cap] if s == edge)
    return {
        "cap_binds": True,
        "selected": cap,
        "cutoff_score": edge,
        "tied_at_cutoff": tied_total,
        "tied_at_cutoff_selected": tied_selected,
        "tie_broken_by_genome_id": tied_total > tied_selected,
    }


def species_report(species, phenotypes, datasets, caps, group_sets):
    drugs = ls.targets(datasets, species)
    class_sizes, labels = ls.labels_by_genome(phenotypes, datasets, species)
    ranking = ls.ranking(phenotypes, datasets, species)
    scores = {gid: s for s, gid in ranking}
    breadth = {gid: len({a for a, _ in lst}) for gid, lst in labels.items()}

    # Genomes with rows on a target drug but no reported R/S label on any:
    # MIC-inferred, intermediate, conflicted or unlabelled only.
    on_targets = {p["genome_id"] for p in phenotypes if p["species"] == species and p["antibiotic"] in drugs}
    unscored = on_targets - set(labels)

    selections = {cap: [gid for _, gid in cap_slice(ranking, cap)] for cap in caps}
    for small, large in zip(caps, caps[1:]):
        assert selections[large][: len(selections[small])] == selections[small], "selections are not nested"

    sp = {
        "species": species,
        "slug": ls.slug(species),
        "target_antibiotics": sorted(drugs),
        "scored_genomes": len(ranking),
        "unscored_genomes": len(unscored),
        "panel_breadth_scored": distribution(list(breadth.values())),
        "by_cap": {},
    }
    for cap in caps:
        sel = selections[cap]
        sp["by_cap"][str(cap)] = {
            **cutoff(ranking, cap),
            "panel_breadth": distribution([breadth[g] for g in sel]),
        }

    tasks = []
    by_drug = {}
    other = {}
    for p in phenotypes:
        if p["species"] != species or p["antibiotic"] not in drugs:
            continue
        call = ls.reported_label(p)
        if call is not None:
            by_drug.setdefault(p["antibiotic"], []).append((p["genome_id"], call))
        else:
            kind = "mic_inferred" if p["label_from_mic"] else "conflict" if p["conflict"] else (
                "intermediate" if p["reported"] == "I" else "unlabelled")
            other.setdefault(p["antibiotic"], []).append((p["genome_id"], kind))

    ids = {d["antibiotic"]: d["id"] for d in datasets if d["species"] == species}
    for drug in sorted(drugs):
        rows = by_drug.get(drug, [])
        n_r = sum(1 for _, c in rows if c == ls.R)
        n_s = len(rows) - n_r
        minority = ls.R if n_r <= n_s else ls.S
        share = {gid: (1.0 / class_sizes[(drug, c)]) / scores[gid] for gid, c in rows}
        excluded = other.get(drug, [])

        def excluded_counts(chosen):
            counts = {"mic_inferred": 0, "intermediate": 0, "conflict": 0, "unlabelled": 0}
            for gid, kind in excluded:
                if chosen is None or gid in chosen:
                    counts[kind] += 1
            return counts

        eligible, _ = rows_summary(rows, minority, breadth, share, group_sets)
        task = {
            "species": species,
            "antibiotic": drug,
            "dataset_id": ids[drug],
            "minority": minority,
            "minority_tied": n_r == n_s,
            "eligible": {**eligible, "non_target_rows": excluded_counts(None)},
            "by_cap": {},
            "nested": [],
        }
        members = {}
        sel_rows = {}
        for cap in caps:
            chosen = set(selections[cap])
            sel_rows[cap] = [(g, c) for g, c in rows if g in chosen]
            summary, members[cap] = rows_summary(sel_rows[cap], minority, breadth, share, group_sets)
            summary["minority_coverage"] = (
                summary["minority_rows"] / eligible["minority_rows"] if eligible["minority_rows"] else None)
            summary["non_target_rows"] = excluded_counts(chosen)
            task["by_cap"][str(cap)] = summary
        for small, large in zip(caps, caps[1:]):
            before = set(sel_rows[small])
            added = [r for r in sel_rows[large] if r not in before]
            delta = {
                "from_cap": small,
                "to_cap": large,
                "genomes_added": len(selections[large]) - len(selections[small]),
                "rows_added": len(added),
                "R_added": sum(1 for _, c in added if c == ls.R),
                "S_added": sum(1 for _, c in added if c == ls.S),
                "minority_rows_added": sum(1 for _, c in added if c == minority),
            }
            if group_sets:
                delta["groups"] = {}
                for name in group_sets:
                    g_small, mg_small = members[small][name]
                    g_large, mg_large = members[large][name]
                    # An ungrouped genome at either cap can make a known group
                    # look new or hide a new one, so both caps must pass the
                    # coverage gate (the ESS null reasons are that gate).
                    at_small = task["by_cap"][str(small)]["groups"][name]
                    at_large = task["by_cap"][str(large)]["groups"][name]
                    reason = pair_reason(small, at_small["group_ess_null_reason"],
                                         large, at_large["group_ess_null_reason"])
                    minority_reason = pair_reason(small, at_small["minority_group_ess_null_reason"],
                                                  large, at_large["minority_group_ess_null_reason"])
                    delta["groups"][name] = {
                        "groups_added": None if reason else len(g_large - g_small),
                        "groups_added_null_reason": reason,
                        "minority_groups_added": None if minority_reason else len(mg_large - mg_small),
                        "minority_groups_added_null_reason": minority_reason,
                    }
            task["nested"].append(delta)
        tasks.append(task)
    return sp, tasks


def build_report(phenotypes, datasets, caps, group_sets=None, species=None):
    group_sets = group_sets or {}
    names = sorted({d["species"] for d in datasets})
    if species:
        missing = set(species) - set(names)
        if missing:
            raise ValueError(f"no dataset for species: {sorted(missing)}")
        names = [n for n in names if n in set(species)]
    report = {
        "schema": SCHEMA,
        "selection": "legacy select_genomes (src/curate.rs at c90ebaf)",
        "caps": caps,
        "group_sets": sorted(group_sets),
        "min_group_coverage": float(MIN_GROUP_COVERAGE),
        "species": [],
        "tasks": [],
    }
    for name in names:
        sp, tasks = species_report(name, phenotypes, datasets, caps, group_sets)
        report["species"].append(sp)
        report["tasks"].extend(tasks)
    return report


def fmt(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def tsv_lines(report):
    """One line per task and cap, the main figures of the JSON report."""
    sets = report["group_sets"]
    head = ["species", "antibiotic", "cap", "minority", "eligible_rows", "eligible_R", "eligible_S",
            "eligible_prevalence_R", "selected_rows", "selected_R", "selected_S", "selected_prevalence_R",
            "minority_coverage", "breadth_median", "task_score_share_mean", "rows_added", "minority_rows_added"]
    for s in sets:
        head += [f"{s}_rows_without_group", f"{s}_group_coverage", f"{s}_groups", f"{s}_minority_groups",
                 f"{s}_largest_share", f"{s}_group_ess", f"{s}_group_ess_null_reason",
                 f"{s}_groups_added", f"{s}_groups_added_null_reason",
                 f"{s}_minority_groups_added", f"{s}_minority_groups_added_null_reason"]
    lines = ["\t".join(head)]
    for t in report["tasks"]:
        e = t["eligible"]
        deltas = {str(d["to_cap"]): d for d in t["nested"]}
        for cap, c in t["by_cap"].items():
            d = deltas.get(cap)
            row = [t["species"], t["antibiotic"], cap, t["minority"], e["rows"], e["R"], e["S"],
                   e["prevalence_R"], c["rows"], c["R"], c["S"], c["prevalence_R"], c["minority_coverage"],
                   c["panel_breadth"]["median"] if c["panel_breadth"] else None, c["task_score_share_mean"],
                   d["rows_added"] if d else None, d["minority_rows_added"] if d else None]
            for s in sets:
                g = c["groups"][s]
                a = d["groups"][s] if d else {}
                row += [g["rows_without_group"], g["group_coverage"], g["groups"], g["minority_groups"],
                        g["largest_group_share"], g["group_ess"], g["group_ess_null_reason"],
                        a.get("groups_added"), a.get("groups_added_null_reason"),
                        a.get("minority_groups_added"), a.get("minority_groups_added_null_reason")]
            lines.append("\t".join(fmt(v) for v in row))
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("out")
    ap.add_argument("--curated")
    ap.add_argument("--phenotypes")
    ap.add_argument("--datasets")
    ap.add_argument("--caps", default="2000,4000")
    ap.add_argument("--species", action="append")
    ap.add_argument("--groups", action="append", default=[], metavar="NAME=TSV")
    ap.add_argument("--tsv")
    a = ap.parse_args(argv)
    phen = a.phenotypes or (a.curated and os.path.join(a.curated, "phenotypes.jsonl"))
    dsets = a.datasets or (a.curated and os.path.join(a.curated, "datasets.json"))
    if not phen or not dsets:
        ap.error("give --curated, or --phenotypes and --datasets")
    group_sets = {}
    for spec in a.groups:
        name, sep, path = spec.partition("=")
        if not sep or not name or name in group_sets:
            ap.error(f"--groups {spec!r}: expected a new NAME=TSV")
        group_sets[name] = read_groups(path)
    datasets = ls.read_datasets(dsets)
    wanted = set(a.species) if a.species else None
    phenotypes = ls.read_phenotypes(phen, wanted)
    report = build_report(phenotypes, datasets, parse_caps(a.caps), group_sets, a.species)
    report["inputs"] = {"phenotypes": os.path.abspath(phen), "datasets": os.path.abspath(dsets),
                        "groups": {n: os.path.abspath(s.partition("=")[2]) for n, s in zip(group_sets, a.groups)}}
    with open(a.out, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=1)
        f.write("\n")
    if a.tsv:
        with open(a.tsv, "w", encoding="utf-8") as f:
            f.write("\n".join(tsv_lines(report)) + "\n")


if __name__ == "__main__":
    main()
