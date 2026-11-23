"""history.py: block M2 of the mechanism split, the mechanism mix over
collection years, for one species.

    python3 -I history.py --species SLUG --reports DIR --lists lists.tsv \
        --drug-rules drug_rules.tsv --tables DIR --db DIR --meta FILE \
        --out DIR [--folds DIR --design design.sha256]

Genomes: every genome of the species with a table that appears in the
samples of any of the species' release reports in DIR. Year: samples[].year,
a text of four digits, else no year; nothing is imputed. BioProject:
`bioproject_accession` of the BV-BRC genome table (--meta), held to its
sha256. Rows: the AMRFinderPlus program's counted rows of `tool_default`.
Family: the database's own `gene_family` (amr_rule.Database.family).

Writes aggregate tables only (PLAN.md, M2): per cell, period pair and
guard, the shares, the differences and the permutation p values, with no
genome id and no BioProject accession; the family mix of the primary guard;
the year counts; the lineage agreement; and, when fold files are given, the
time-split link for the datasets of the species that have them.

Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import gzip
import hashlib
import json
import os
import random
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import amr_rule as ar  # noqa: E402
from common import r12, read_tsv, seed_of, shown, write_tsv  # noqa: E402

META_SHA256 = "84ae57a777cb234aa0a39102ab7f9ca9448e52d1474654179f26165644c25f64"
DESIGN_SHA256 = "cc4635ed1e9741fbd58996d3d1731fd15211363ac3d20c836ef372718308878e"
# name -> ((first, last) earlier, (first, last) later), both ends inclusive
PERIODS = {
    "primary": ((2010, 2014), (2015, 2019)),
    "secondary": ((2000, 2009), (2010, 2019)),
}
MIN_UNIT_GENOMES = 5   # a unit votes in a period with 5 or more genomes
MIN_UNITS = 3          # a cell passes with 3 or more voting units per period
PERMUTATIONS = 1000
FLAG_SHARE = 0.95
GUARDS = ("bioproject", "pooled", "largest_out", "flagged_out", "lineage", "lineage_both", "countries_both",
          "intact")
YEAR_RE = re.compile(r"[0-9]{4}")

CELLS_HEADER = ["species_slug", "species", "cell", "class", "dataset", "antibiotic", "periods", "guard",
                "genomes", "with_year_share", "earlier_genomes", "later_genomes", "earlier_units", "later_units",
                "passes", "earlier_point", "later_point", "d_point", "p_point", "earlier_point_disrupt",
                "later_point_disrupt", "earlier_gene", "later_gene", "d_gene", "tvd", "p_tvd", "note"]
FAMILIES_HEADER = ["species_slug", "species", "class", "periods", "family", "earlier_genomes", "later_genomes",
                   "earlier_mix", "later_mix"]
# A family seen in 1 to 4 genomes of a period's voting BioProjects is not
# named: it joins one row under this label, with its genome counts as they are.
FOLDED_FAMILIES = "other families"
YEARS_HEADER = ["species_slug", "species", "years", "genomes"]
LINEAGE_HEADER = ["species_slug", "species", "reports", "report_pairs", "differing_genome_pairs", "given"]
TIMESPLIT_HEADER = ["dataset", "species_slug", "species", "antibiotic", "reading", "resistant_later_test",
                    "with_determinant", "with_unseen_determinant", "no_determinant", "share_unseen"]


# ---------------------------------------------------------------- genomes

def year_of(text):
    """A year from samples[].year: four ASCII digits and nothing else, not
    even a blank, else None (the plan: anything else is no year)."""
    if isinstance(text, str) and YEAR_RE.fullmatch(text):
        return int(text)
    return None


def period_of(year, periods):
    """'earlier', 'later' or None for a year under a period pair."""
    if year is None:
        return None
    (a0, a1), (b0, b1) = PERIODS[periods]
    if a0 <= year <= a1:
        return "earlier"
    if b0 <= year <= b1:
        return "later"
    return None


def year_bucket(year):
    if year is None:
        return "no year"
    if year < 2000:
        return "before 2000"
    if year >= 2020:
        return "2020 on"
    if year < 2010:
        return "2000 to 2009"
    if year < 2015:
        return "2010 to 2014"
    return "2015 to 2019"


def bioproject_unit(genome):
    """The voting unit of the BioProject vote: the genome's BioProject."""
    return genome["bioproject"]


def read_meta(path, wanted):
    """genome id -> bioproject_accession for the genomes in `wanted`."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    if h.hexdigest() != META_SHA256:
        ar.refuse("the BV-BRC genome table differs from its recorded sha256")
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            row = json.loads(line)
            gid = str(row.get("genome_id", ""))
            if gid in wanted:
                bp = (row.get("bioproject_accession") or "").strip()
                out[gid] = bp or None
    return out


def load_reports(directory, species):
    reports = []
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".json"):
            continue
        with open(os.path.join(directory, name), encoding="utf-8") as f:
            rep = json.load(f)
        if rep.get("species_slug") != species:
            ar.refuse(f"{name}: a report of another species is in the species' folder")
        reports.append(rep)
    if not reports:
        ar.refuse(f"{species}: no report")
    return reports


def genome_table(reports, tables_dir, meta):
    """id -> {year, country, bioproject}, over the samples of every report,
    for genomes with a table; and the count of genomes with none."""
    seen, no_table = {}, set()
    for rep in reports:
        for s in rep["samples"]:
            gid = s["id"]
            y = year_of(s.get("year"))
            if gid in seen:
                if seen[gid]["year"] != y:
                    ar.refuse(f"{rep['species_slug']}: a genome's year differs between two reports")
                continue
            if gid in no_table:
                continue
            if not os.path.isfile(os.path.join(tables_dir, gid + ".tsv")):
                no_table.add(gid)
                continue
            seen[gid] = {"year": y, "country": s.get("country") or None, "bioproject": None}
    for gid, g in seen.items():
        g["bioproject"] = meta.get(gid)
    return seen, len(no_table)


def lineage_components(reports, genomes):
    """Join lineage clusters over the reports. Returns (component per genome
    or None where the relation differs, report pairs, differing pairs)."""
    parent = {}

    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    by_report = []
    for i, rep in enumerate(reports):
        m = {s["id"]: s.get("lineage") for s in rep["samples"] if s["id"] in genomes and s.get("lineage") is not None}
        by_report.append(m)
        first = {}
        for gid, c in m.items():
            if c in first:
                parent[find(gid)] = find(first[c])
            else:
                first[c] = gid
    pairs, differ = 0, 0
    for i in range(len(by_report)):
        for j in range(i + 1, len(by_report)):
            a, b = by_report[i], by_report[j]
            shared = [g for g in a if g in b]
            if not shared:
                continue
            pairs += 1
            na, nb, nab = collections.Counter(), collections.Counter(), collections.Counter()
            for g in shared:
                na[a[g]] += 1
                nb[b[g]] += 1
                nab[(a[g], b[g])] += 1
            same_a = sum(n * (n - 1) // 2 for n in na.values())
            same_b = sum(n * (n - 1) // 2 for n in nb.values())
            same_ab = sum(n * (n - 1) // 2 for n in nab.values())
            differ += same_a + same_b - 2 * same_ab
    if differ:
        return None, pairs, differ
    return {g: find(g) for g in genomes if g in parent}, pairs, 0


# ---------------------------------------------------------------- calls

def genome_cells(rows, db, variant="tool_default", rule=None):
    """class -> (point, gene, families, POINT_DISRUPT) for one genome. A row
    of Subtype POINT or POINT_DISRUPT is point and enters no family mix; a
    row of Subtype AMR is gene, with its family (the lead's rulings on
    questions 1 and 10). With `rule`, only the rows that meet the rule count
    and the class is the single key ''."""
    out = {}
    use = ar.rule_rows(rows, rule, variant) if rule is not None else ar.counted(rows, variant)
    for r in use:
        keys = [""] if rule is not None else ar.classes_of(r)
        for c in keys:
            p, g, fams, d = out.get(c, (False, False, frozenset(), False))
            if ar.is_point(r):
                p = True
                d = d or ar.is_disrupt(r)
            else:
                g = True
                fams = fams | {db.family(r["symbol"], r["subtype"])}
            out[c] = (p, g, fams, d)
    return out


def mix_of(families):
    """A genome's mix vector over its families of a cell: weight 1 in all."""
    return {f: 1.0 / len(families) for f in families} if families else None


# ---------------------------------------------------------------- the vote

def units_of(gids, genomes, periods, unit_fn, cell_of):
    """Voting units of a cell: [(period, unit key, [cell values])], each
    with 5 or more genomes with a year in that period; and the genomes that
    enter only the pooled share."""
    groups = collections.defaultdict(list)
    for gid in gids:
        per = period_of(genomes[gid]["year"], periods)
        if per is None:
            continue
        groups[(per, unit_fn(genomes[gid]) if unit_fn else None)].append(cell_of(gid))
    units = []
    for (per, key), vals in sorted(groups.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
        if key is not None and len(vals) >= MIN_UNIT_GENOMES:
            units.append((per, key, vals))
    return units


def unit_summary(vals):
    """(point share, gene share, mix or None) of one unit's genomes."""
    n = len(vals)
    point = sum(v[0] for v in vals) / n
    gene = sum(v[1] for v in vals) / n
    mixes = [mix_of(v[2]) for v in vals if v[1] and v[2]]
    if not mixes:
        return point, gene, None
    m = collections.Counter()
    for x in mixes:
        for f, w in x.items():
            m[f] += w / len(mixes)
    return point, gene, dict(m)


def mean_mix(mixes):
    mixes = [m for m in mixes if m is not None]
    if not mixes:
        return None
    out = collections.Counter()
    for m in mixes:
        for f, w in m.items():
            out[f] += w / len(mixes)
    return dict(out)


def tvd(a, b):
    if a is None or b is None:
        return None
    keys = set(a) | set(b)
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in keys)


def statistics_of(summaries, labels):
    """(earlier point, later point, d point, earlier gene, later gene,
    d gene, tvd, earlier mix, later mix) for unit summaries under labels."""
    e = [s for s, lab in zip(summaries, labels) if lab == "earlier"]
    l_ = [s for s, lab in zip(summaries, labels) if lab == "later"]
    if not e or not l_:
        return None
    ep, lp = sum(s[0] for s in e) / len(e), sum(s[0] for s in l_) / len(l_)
    eg, lg = sum(s[1] for s in e) / len(e), sum(s[1] for s in l_) / len(l_)
    em, lm = mean_mix([s[2] for s in e]), mean_mix([s[2] for s in l_])
    return ep, lp, lp - ep, eg, lg, lg - eg, tvd(em, lm), em, lm


def permutation_p(units, observed, seed_name, statistic):
    """p from PERMUTATIONS shuffles of period labels over the voting units,
    so a unit's genomes move together and the number of later units is
    fixed. 'point' is two-sided on the difference, 'tvd' one-sided."""
    if observed is None:
        return None
    summaries = [unit_summary(u[2]) for u in units]
    rng = random.Random(seed_of(seed_name))
    lab = [u[0] for u in units]
    hits = 0
    obs = r12(abs(observed) if statistic == "point" else observed)
    for _ in range(PERMUTATIONS):
        rng.shuffle(lab)
        st = statistics_of(summaries, lab)
        v = st[2] if statistic == "point" else st[6]
        if v is None:
            continue
        if r12(abs(v) if statistic == "point" else v) >= obs:
            hits += 1
    return r12((1 + hits) / (PERMUTATIONS + 1))


def fmt(x):
    return "" if x is None else f"{r12(x):.12f}"


def entry_passes(units):
    """The plan's entry rule: 3 or more voting units in each period."""
    ne = sum(1 for u in units if u[0] == "earlier")
    nl = sum(1 for u in units if u[0] == "later")
    return ne >= MIN_UNITS and nl >= MIN_UNITS


def vote_row(base, guard, units, seed_stub, note=""):
    """One output row of a voted guard from its units."""
    ne = sum(1 for u in units if u[0] == "earlier")
    nl = sum(1 for u in units if u[0] == "later")
    passes = entry_passes(units)
    row = dict(base, guard=guard, earlier_units=shown(ne), later_units=shown(nl), passes=int(passes),
               earlier_genomes=shown(sum(len(u[2]) for u in units if u[0] == "earlier")),
               later_genomes=shown(sum(len(u[2]) for u in units if u[0] == "later")),
               earlier_point_disrupt=shown(sum(v[3] for u in units if u[0] == "earlier" for v in u[2])),
               later_point_disrupt=shown(sum(v[3] for u in units if u[0] == "later" for v in u[2])), note=note)
    if not passes:
        # under 3 voting BioProjects in a period, a share would be one or two
        # BioProjects' own: no share, mix or p is written
        row["note"] = (note + "; " if note else "") + "under 3 voting BioProjects in a period, no share written"
        return row, None
    summaries = [unit_summary(u[2]) for u in units]
    labels = [u[0] for u in units]
    ep, lp, dp, eg, lg, dg, tv, em, lm = statistics_of(summaries, labels)
    row.update(earlier_point=fmt(ep), later_point=fmt(lp), d_point=fmt(dp), earlier_gene=fmt(eg),
               later_gene=fmt(lg), d_gene=fmt(dg), tvd=fmt(tv))
    row["p_point"] = fmt(permutation_p(units, dp, seed_stub + f"|{guard}|point", "point"))
    row["p_tvd"] = fmt(permutation_p(units, tv, seed_stub + f"|{guard}|tvd", "tvd"))
    counts = tuple(collections.Counter(f for u in units if u[0] == period for v in u[2] for f in v[2])
                   for period in ("earlier", "later"))
    holders = tuple([v for u in units if u[0] == period for v in u[2]] for period in ("earlier", "later"))
    return row, (em, lm, counts, holders)


def family_rows(base, em, lm, counts, holders):
    """m2-families rows of a passing cell: a family by name where its genome
    count in each period is 0 or 5 or more, the rest folded into one row.
    `holders` are the period's genome cells; each genome is in one unit."""
    named, folded = [], set()
    for fam in sorted(set(em) | set(lm)):
        n = (counts[0].get(fam, 0), counts[1].get(fam, 0))
        if any(0 < k < MIN_UNIT_GENOMES for k in n):
            folded.add(fam)
        else:
            named.append(fam)
    out = [dict(base, family=fam, earlier_genomes=shown(counts[0].get(fam, 0)),
                later_genomes=shown(counts[1].get(fam, 0)), earlier_mix=fmt(em.get(fam, 0.0)),
                later_mix=fmt(lm.get(fam, 0.0))) for fam in named]
    if folded:
        row = dict(base, family=FOLDED_FAMILIES)
        for period, mix, held in (("earlier", em, holders[0]), ("later", lm, holders[1])):
            k = sum(1 for v in held if v[2] & folded)
            row[f"{period}_genomes"] = shown(k)
            row[f"{period}_mix"] = fmt(sum(mix.get(f, 0.0) for f in folded))
        out.append(row)
    return out


def pooled_row(base, gids, genomes, periods, cell_of, vote_passes):
    """Guard 1: the pooled share per genome, no vote and no p. Its shares are
    written only where the cell's BioProject vote passes the plan's entry
    rule (the lead's ruling on question 12)."""
    per = {"earlier": [], "later": []}
    for gid in gids:
        p = period_of(genomes[gid]["year"], periods)
        if p:
            per[p].append(cell_of(gid))
    row = dict(base, guard="pooled", earlier_genomes=shown(len(per["earlier"])),
               later_genomes=shown(len(per["later"])), passes="", note="per genome; no vote, no p",
               earlier_point_disrupt=shown(sum(v[3] for v in per["earlier"])),
               later_point_disrupt=shown(sum(v[3] for v in per["later"])))
    if not vote_passes:
        row["note"] += "; the cell's BioProject vote does not pass, no share written"
        return row
    st = statistics_of([unit_summary(per["earlier"]), unit_summary(per["later"])], ["earlier", "later"])
    ep, lp, dp, eg, lg, dg, tv, _, _ = st
    row.update(earlier_point=fmt(ep), later_point=fmt(lp), d_point=fmt(dp), earlier_gene=fmt(eg),
               later_gene=fmt(lg), d_gene=fmt(dg), tvd=fmt(tv))
    return row


# ---------------------------------------------------------------- guards

def flagged_bioprojects(cls, gids, genomes, cells, labelled, db):
    """Guard 3: BioProjects with 0.95 or more resistant over their labelled
    genome and dataset pairs whose rule names the class, or with 0.95 or more
    of their genomes with a gene row of the cell carrying one family."""
    by_bp = collections.defaultdict(list)
    for gid in gids:
        bp = genomes[gid]["bioproject"]
        if bp:
            by_bp[bp].append(gid)
    flagged = set()
    for bp, members in by_bp.items():
        res = tot = 0
        for rule, labels in labelled:
            if not db.rule_names_class(rule, cls):
                continue
            for gid in members:
                if gid in labels:
                    tot += 1
                    res += labels[gid]
        if tot and res / tot >= FLAG_SHARE:
            flagged.add(bp)
            continue
        with_gene = [cells[gid][cls][2] for gid in members if cls in cells[gid] and cells[gid][cls][1]]
        if with_gene:
            fam = collections.Counter(f for fams in with_gene for f in fams)
            if max(fam.values()) / len(with_gene) >= FLAG_SHARE:
                flagged.add(bp)
    return flagged


def countries_in_both(gids, genomes, periods):
    have = {"earlier": set(), "later": set()}
    for gid in gids:
        p = period_of(genomes[gid]["year"], periods)
        c = genomes[gid]["country"]
        if p and c:
            have[p].add(c)
    return have["earlier"] & have["later"]


def cell_rows(base, gids, genomes, cells_main, cells_intact, cls, periods, lineage, lineage_note, labelled, db,
              seed_stub, guards=GUARDS):
    """Every guard's row for one cell and one period pair."""
    def cell_of_from(cells):
        def f(gid):
            return cells[gid].get(cls, (False, False, frozenset(), False))
        return f

    cell_of = cell_of_from(cells_main)
    out, mixes = [], None
    units = units_of(gids, genomes, periods, bioproject_unit, cell_of)
    if "bioproject" in guards:
        row, mixes = vote_row(base, "bioproject", units, seed_stub)
        out.append(row)
    if "pooled" in guards:
        out.append(pooled_row(base, gids, genomes, periods, cell_of, entry_passes(units)))
    if "largest_out" in guards:
        count = collections.Counter(genomes[g]["bioproject"] for g in gids
                                    if genomes[g]["bioproject"] and period_of(genomes[g]["year"], periods))
        if count:
            top = max(count.values())
            largest = sorted(bp for bp, n in count.items() if n == top)[0]
            keep = [g for g in gids if genomes[g]["bioproject"] != largest]
            out.append(vote_row(base, "largest_out", units_of(keep, genomes, periods, bioproject_unit, cell_of),
                                seed_stub)[0])
        else:
            out.append(dict(base, guard="largest_out", note="no BioProject in the periods"))
    if "flagged_out" in guards:
        flagged = flagged_bioprojects(cls, gids, genomes, cells_main, labelled, db)
        keep = [g for g in gids if genomes[g]["bioproject"] not in flagged]
        out.append(vote_row(base, "flagged_out", units_of(keep, genomes, periods, bioproject_unit, cell_of),
                            seed_stub, note=f"BioProjects flagged: {shown(len(flagged))}")[0])
    if "lineage" in guards:
        if lineage is None:
            out.append(dict(base, guard="lineage", note=lineage_note))
            out.append(dict(base, guard="lineage_both", note=lineage_note))
        else:
            with_lin = [g for g in gids if g in lineage]
            units = units_of_ids(with_lin, genomes, periods, lineage, cell_of)
            out.append(vote_row(base, "lineage", units, seed_stub, note=lineage_note)[0])
            both = {u[1] for u in units if u[0] == "earlier"} & {u[1] for u in units if u[0] == "later"}
            out.append(vote_row(base, "lineage_both", [u for u in units if u[1] in both], seed_stub,
                                note="clusters voting in both periods")[0])
    if "countries_both" in guards:
        countries = countries_in_both(gids, genomes, periods)
        keep = [g for g in gids if genomes[g]["country"] in countries]
        out.append(vote_row(base, "countries_both", units_of(keep, genomes, periods, bioproject_unit, cell_of),
                            seed_stub, note=f"countries in both periods: {shown(len(countries))}")[0])
    if "intact" in guards and cells_intact is not None:
        units = units_of(gids, genomes, periods, bioproject_unit, cell_of_from(cells_intact))
        out.append(vote_row(base, "intact", units, seed_stub, note="rows of no partial method")[0])
    return out, mixes


def units_of_ids(gids, genomes, periods, unit_by_id, cell_of):
    """units_of with the unit read from a map by genome id (lineage clusters)."""
    groups = collections.defaultdict(list)
    for gid in gids:
        per = period_of(genomes[gid]["year"], periods)
        if per is None:
            continue
        groups[(per, unit_by_id[gid])].append(cell_of(gid))
    return [(per, key, vals) for (per, key), vals in sorted(groups.items(), key=lambda kv: (kv[0][0], str(kv[0][1])))
            if len(vals) >= MIN_UNIT_GENOMES]


# ---------------------------------------------------------------- time split

def timesplit_rows(species, species_name, reports, rules, tables, db, folds_dir):
    out = []
    for rep in reports:
        d = rep["id"]
        path = os.path.join(folds_dir, d + ".past.json")
        if not os.path.isfile(path) or d not in rules:
            continue
        with open(path, encoding="utf-8") as f:
            fold = json.load(f)
        if fold.get("format") != "grm-fold-file-v1" or fold.get("dataset") != d:
            ar.refuse(f"{d}: the past arm's fold file is not grm-fold-file-v1 for this dataset")
        rule = ar.parse_rule(rules[d]["rule"])
        label = {s["id"]: s["resistant"] for s in rep["samples"]}
        train = [g for g in fold["train"] if label.get(g) and g in tables]
        test = [g for fold_ in fold["test"] for g in fold_ if label.get(g) and g in tables]

        for reading in ("point_as_one", "point_by_symbol"):
            def dets(gid):
                out_ = set()
                for r in ar.rule_rows(tables[gid], rule):
                    if ar.is_point(r):
                        out_.add("point" if reading == "point_as_one" else "point:" + r["symbol"])
                    else:
                        out_.add(db.family(r["symbol"], r["subtype"]))
                return out_

            seen = set()
            for g in train:
                seen |= dets(g)
            with_det = unseen = none = 0
            for g in test:
                ds = dets(g)
                if not ds:
                    none += 1
                    continue
                with_det += 1
                unseen += bool(ds - seen)
            out.append({"dataset": d, "species_slug": species, "species": species_name,
                        "antibiotic": rep.get("antibiotic", ""), "reading": reading,
                        "resistant_later_test": shown(len(test)), "with_determinant": shown(with_det),
                        "with_unseen_determinant": shown(unseen), "no_determinant": shown(none),
                        "share_unseen": fmt(unseen / with_det) if with_det >= MIN_UNIT_GENOMES else ""})
    return out


def sha256_of(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def check_design(folds_dir, design):
    if sha256_of(design) != DESIGN_SHA256:
        ar.refuse("design.sha256 differs from its recorded sha256")
    want = {}
    with open(design, encoding="utf-8") as f:
        lines = f.read().splitlines()
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        digest, name = line.split(None, 1)
        name = name.strip()
        if name.startswith("results/fold-files/"):
            want[name[len("results/fold-files/"):]] = digest
    for name in os.listdir(folds_dir):
        if name not in want:
            ar.refuse(f"fold file {name} is not in design.sha256")
        if sha256_of(os.path.join(folds_dir, name)) != want[name]:
            ar.refuse(f"fold file {name} differs from design.sha256")


# ---------------------------------------------------------------- run

def run(args):
    reports = load_reports(args.reports, args.species)
    species_name = reports[0].get("species", "")
    ids = {s["id"] for rep in reports for s in rep["samples"]}
    meta = read_meta(args.meta, ids)
    genomes, no_table = genome_table(reports, args.tables, meta)
    db = ar.Database(args.db)
    tables = {gid: ar.read_table(os.path.join(args.tables, gid + ".tsv")) for gid in genomes}
    cells_main = {gid: genome_cells(rows, db) for gid, rows in tables.items()}
    cells_intact = {gid: genome_cells(rows, db, "tool_default_intact") for gid, rows in tables.items()}
    lineage, pairs, differ = lineage_components(reports, genomes)
    lineage_note = (f"same-cluster relation differs in {differ} genome pairs; not given" if lineage is None
                    else f"clusters joined over {len(reports)} reports")
    lists = {r["dataset"]: r for r in read_tsv(args.lists) if r["species_slug"] == args.species}
    rules = {r["dataset"]: r for r in read_tsv(args.drug_rules) if r["species_slug"] == args.species}
    primary = [rep for rep in reports if rep["id"] in lists and lists[rep["id"]]["in_A"] == "1"
               and lists[rep["id"]]["in_M"] != "1"]
    labelled = [(rep["catalog"]["rule"], {s["id"]: int(s["resistant"]) for s in rep["samples"]}) for rep in primary
                if rep.get("catalog")]
    gids = sorted(genomes)
    with_year = sum(1 for g in gids if genomes[g]["year"] is not None)
    classes = sorted({c for v in cells_main.values() for c in v})
    rows, fam_rows = [], []
    for cls in classes:
        for periods in PERIODS:
            base = {"species_slug": args.species, "species": species_name, "cell": "species_class", "class": cls,
                    "periods": periods, "genomes": shown(len(gids)),
                    "with_year_share": fmt(with_year / len(gids)) if gids else ""}
            stub = f"m2|{args.species}|{cls}|{periods}"
            out, mixes = cell_rows(base, gids, genomes, cells_main, cells_intact, cls, periods, lineage,
                                   lineage_note, labelled, db, stub)
            rows.extend(out)
            if mixes and mixes[0] is not None and mixes[1] is not None:
                fam_rows.extend(family_rows({"species_slug": args.species, "species": species_name, "class": cls,
                                             "periods": periods}, *mixes))
    # secondary: resistant genomes only, per dataset of A not in M, with the
    # dataset's own rule-meeting rows
    for rep in primary:
        d = rep["id"]
        res = sorted(s["id"] for s in rep["samples"] if s["resistant"] and s["id"] in genomes)
        cells_d = {g: genome_cells(tables[g], db, rule=rep["catalog"]["rule"]) for g in res}
        for periods in PERIODS:
            base = {"species_slug": args.species, "species": species_name, "cell": "dataset_resistant",
                    "class": "", "dataset": d, "antibiotic": rep.get("antibiotic", ""), "periods": periods,
                    "genomes": shown(len(res)),
                    "with_year_share": fmt(sum(genomes[g]["year"] is not None for g in res) / len(res)) if res else ""}
            out, _ = cell_rows(base, res, genomes, cells_d, None, "", periods, None, "", [], db,
                               f"m2|{args.species}|{d}|{periods}", guards=("bioproject", "pooled"))
            rows.extend(out)
    years = collections.Counter(year_bucket(genomes[g]["year"]) for g in gids)
    year_rows = [{"species_slug": args.species, "species": species_name, "years": b, "genomes": shown(years[b])}
                 for b in ("no year", "before 2000", "2000 to 2009", "2010 to 2014", "2015 to 2019", "2020 on")]
    year_rows.append({"species_slug": args.species, "species": species_name, "years": "no table",
                      "genomes": shown(no_table)})
    lin_rows = [{"species_slug": args.species, "species": species_name, "reports": len(reports),
                 "report_pairs": pairs, "differing_genome_pairs": differ, "given": int(lineage is not None)}]
    ts_rows = []
    if args.folds:
        check_design(args.folds, args.design)
        ts_rows = timesplit_rows(args.species, species_name, reports, rules, tables, db, args.folds)
    os.makedirs(args.out, exist_ok=True)
    write_tsv(os.path.join(args.out, "m2-cells.tsv"), CELLS_HEADER, rows)
    write_tsv(os.path.join(args.out, "m2-families.tsv"), FAMILIES_HEADER, fam_rows)
    write_tsv(os.path.join(args.out, "m2-years.tsv"), YEARS_HEADER, year_rows)
    write_tsv(os.path.join(args.out, "m2-lineage.tsv"), LINEAGE_HEADER, lin_rows)
    write_tsv(os.path.join(args.out, "timesplit.tsv"), TIMESPLIT_HEADER, ts_rows)
    passing = sum(1 for r in rows if r.get("guard") == "bioproject" and r.get("cell") == "species_class"
                  and r.get("periods") == "primary" and r.get("passes") == 1)
    cells = sum(1 for r in rows if r.get("guard") == "bioproject" and r.get("cell") == "species_class"
                and r.get("periods") == "primary")
    print(f"history: {args.species}: {len(gids)} genomes with a table, {len(classes)} classes, {cells} cells, "
          f"{passing} pass the BioProject vote on the primary periods; {len(rows)} rows; "
          f"lineage given {int(lineage is not None)}; time-split rows {len(ts_rows)}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser()
    for a in ("--species", "--reports", "--lists", "--drug-rules", "--tables", "--db", "--meta", "--out"):
        ap.add_argument(a, required=True)
    ap.add_argument("--folds")
    ap.add_argument("--design")
    args = ap.parse_args(argv)
    if bool(args.folds) != bool(args.design):
        print("history.py: --folds and --design go together", file=sys.stderr)
        return 2
    try:
        return run(args)
    except ar.Refusal as e:
        print(f"history.py: refused: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
