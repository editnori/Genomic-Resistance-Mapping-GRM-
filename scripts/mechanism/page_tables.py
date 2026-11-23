"""page_tables.py: the tables of docs/mechanism-split-result.md, blocks M1 and
M2 of the mechanism split read once, as files and as the page prints them.

    python3 -I page_tables.py write READING LISTS CONTRASTS COMMIT OUT
    python3 -I page_tables.py print OUT

READING is the reading's directory:

    reading.txt     the reading's commands, each with the sha256 of what it wrote
    checks.txt      the checks of the run, counts only
    combined/       what combine.py --run ... wrote: m1-reading.tsv,
                    m1-neither.tsv, m2-reading.tsv, m2-summary.tsv,
                    timesplit-link.tsv
    runs/<slug>/    one species' copy of run_species.sh's OUT: DONE and
                    results/ (m1/ and m2/ and SHA256SUMS)

LISTS is scripts/c3-inputs/lists.tsv, CONTRASTS the time split's
reference/time-split/contrasts.tsv, COMMIT the commit every DONE must name.

combine.py's own five tables are the reading; this script does not compute
them again in another way. `write` checks them and copies them under
combined/ as they are. It writes beside them the page's other tables, each a
join or a concatenation of the runs' own tables, species and antibiotic in
full, with no new measure: page-species.tsv, page-datasets.tsv,
page-m1-datasets.tsv, page-m1-strata.tsv, page-m1-pooled.tsv,
page-m2-cells.tsv, page-m2-families.tsv, page-m2-years.tsv,
page-m2-lineage.tsv and page-timesplit.tsv. Then gates.txt, README.md and
SHA256SUMS. It passes every gate of GATES before it writes a file and stops
at the first that fails, naming it. OUT must be absent or empty.

`print` prints the page's tables from OUT alone, figures rounded to 3
decimals. Every row is an aggregate over a dataset, a species, a cell or a
learner: no file holds a genome id, a BioProject accession, a sequence or a
position.

Python 3.12, standard library only; run with -I.
"""

import contextlib
import csv
import hashlib
import io
import os
import sys
import tempfile
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import combine  # noqa: E402
import scan  # noqa: E402

LEARNERS = combine.LEARNERS
ROWS = combine.ROWS
COMBINED = ("m1-reading.tsv", "m1-neither.tsv", "m2-reading.tsv", "m2-summary.tsv", "timesplit-link.tsv")
RESULTS = ("m1/datasets.tsv", "m1/nocall.tsv", "m1/pooled.tsv", "m1/strata.tsv", "m1/summary.json",
           "m2/m2-cells.tsv", "m2/m2-families.tsv", "m2/m2-lineage.tsv", "m2/m2-years.tsv", "m2/timesplit.tsv")
STRATA = ("point_only", "gene_only", "both", "neither")
BESIDE = ("partial_only:point_only", "partial_only:gene_only", "partial_only:both",
          "with_point_disrupt:point_only", "with_point_disrupt:gene_only", "with_point_disrupt:both",
          "point_disrupt_alone:point_only", "point_disrupt_alone:both", "all")
RECALL_FLOOR = 10
# The one M2 reading that decides: the BioProject vote on the primary periods over species and class cells.
DECIDING = ("species_class", "primary", "bioproject")
GUARDS = ("bioproject", "pooled", "largest_out", "flagged_out", "lineage", "lineage_both", "countries_both", "intact")
GATES = [
    (1, "layout", "the reading holds reading.txt, checks.txt, combine.py's five tables under combined/ and runs/"),
    (2, "species", "runs/ holds exactly the species of A in lists.tsv, none dropped and none more"),
    (3, "runs", "each run's DONE names its species, the reading's commit, 10 results and the sha256 of its "
                "SHA256SUMS; SHA256SUMS lists the ten tables and each file equals its line; no run is a first-run copy"),
    (4, "datasets", "each run's datasets.tsv holds exactly its species' datasets of A, species and antibiotic as "
                    "lists.tsv names them, set M exactly where lists.tsv has in_M 1; strata.tsv and pooled.tsv hold "
                    "only entered datasets, each with the set datasets.tsv gives it"),
    (5, "floors", "every run table passes scan.py's table checks: a recall only over 10 or more resistant genomes, an "
                  "M2 share or mix only in a cell whose BioProject vote passes, a family named only over 5 or more "
                  "genomes"),
    (6, "combined", "combined/ is byte for byte what combine.py writes from runs/, species in sorted order, and the "
                    "time split's contrasts"),
]
LEARNER_WORDS = {"scm": "set covering machine (`scm`)", "cart": "decision tree (`cart`)",
                 "rf": "random forest (`rf`)", "gbt": "boosted trees (`gbt`)", "lasso": "lasso (`lasso`)",
                 "catalog": "GRM's known-gene rule", "program": "the AMRFinderPlus program"}
STRATUM_WORDS = {"point_only": "point only", "gene_only": "gene only", "both": "both", "neither": "neither",
                 "all": "all genomes"}
GUARD_WORDS = {"bioproject": "one vote per BioProject (decides)", "pooled": "share per genome, no vote",
               "largest_out": "largest BioProject left out", "flagged_out": "flagged BioProjects left out",
               "lineage": "one vote per lineage cluster", "lineage_both": "lineage clusters voting in both periods",
               "countries_both": "countries in both periods only", "intact": "rows of no partial method only"}
PERIOD_WORDS = {"primary": "2010 to 2014 against 2015 to 2019", "secondary": "2000 to 2009 against 2010 to 2019"}


class Stop(Exception):
    pass


# ---------- reading files ----------

def load(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))


def sha256_of(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def a_rows(lrows):
    return [r for r in lrows if r["in_A"] == "1"]


def a_species(lrows):
    """The species slugs of A, sorted."""
    return sorted({r["species_slug"] for r in a_rows(lrows)})


def in_m(lrow):
    return lrow["in_M"] == "1"


def species_name(lrows, slug):
    """A species named in full, as lists.tsv names it."""
    names = {r["species"] for r in lrows if r["species_slug"] == slug}
    if len(names) != 1:
        raise Stop(f"lists.tsv names {slug} {len(names)} ways")
    return names.pop()


# ---------- gates ----------

def gate_layout(reading):
    for need in ("reading.txt", "checks.txt"):
        if not (reading / need).is_file():
            raise Stop(f"gate 1 layout: no {need}")
    for name in COMBINED:
        if not (reading / "combined" / name).is_file():
            raise Stop(f"gate 1 layout: combined/{name} is missing")
    extra = sorted(p.name for p in (reading / "combined").iterdir() if p.name not in COMBINED)
    if extra:
        raise Stop(f"gate 1 layout: combined/ holds files combine.py does not write: {', '.join(extra)}")
    if not (reading / "runs").is_dir():
        raise Stop("gate 1 layout: no runs/")


def gate_species(reading, lrows):
    """The species of runs/, which must be exactly those of A."""
    got = sorted(p.name for p in (reading / "runs").iterdir() if p.is_dir())
    want = a_species(lrows)
    if got != want:
        missing, more = sorted(set(want) - set(got)), sorted(set(got) - set(want))
        raise Stop(f"gate 2 species: runs/ holds {len(got)} species of {len(want)}; missing {missing}, more {more}")
    return got


def parse_done(text):
    """DONE <species> <time> commit <sha> results <n> SHA256SUMS <sha256>."""
    p = text.split()
    if len(p) != 9 or p[0] != "DONE" or p[3] != "commit" or p[5] != "results" or p[7] != "SHA256SUMS":
        return None
    return {"species": p[1], "time": p[2], "commit": p[4], "results": p[6], "sums": p[8]}


def gate_runs(reading, species, commit):
    for slug in species:
        run = reading / "runs" / slug
        if not (run / "DONE").is_file():
            raise Stop(f"gate 3 runs: {slug} has no DONE")
        done = parse_done((run / "DONE").read_text(encoding="utf-8"))
        if done is None or done["species"] != slug:
            raise Stop(f"gate 3 runs: {slug}: DONE is not a DONE line of this species")
        if done["commit"] != commit:
            raise Stop(f"gate 3 runs: {slug}: DONE names commit {done['commit'][:12]}, not the reading's "
                       f"{commit[:12]}; a run of another commit is not read")
        sums = run / "results" / "SHA256SUMS"
        if not sums.is_file() or sha256_of(sums) != done["sums"]:
            raise Stop(f"gate 3 runs: {slug}: results/SHA256SUMS differs from the sha256 DONE names")
        lines = sums.read_text(encoding="utf-8").splitlines()
        listed = {}
        for line in lines:
            digest, name = line.split("  ", 1)
            listed[name[2:] if name.startswith("./") else name] = digest
        if sorted(listed) != sorted(RESULTS) or done["results"] != str(len(RESULTS)):
            raise Stop(f"gate 3 runs: {slug}: SHA256SUMS does not list the run's ten tables")
        here = sorted(p.relative_to(run / "results").as_posix() for p in (run / "results").rglob("*")
                      if p.is_file() and p.name != "SHA256SUMS")
        if here != sorted(RESULTS):
            raise Stop(f"gate 3 runs: {slug}: results/ holds other files than its SHA256SUMS lists")
        for name, digest in listed.items():
            if sha256_of(run / "results" / name) != digest:
                raise Stop(f"gate 3 runs: {slug}: results/{name} differs from its SHA256SUMS line")


def read_runs(reading, species):
    """{slug: {table name: rows}} of the runs' tables."""
    out = {}
    for slug in species:
        res = reading / "runs" / slug / "results"
        out[slug] = {name: load(res / name) for name in RESULTS if name.endswith(".tsv")}
    return out


def gate_datasets(runs, lrows):
    for slug, t in runs.items():
        want = {r["dataset"]: r for r in a_rows(lrows) if r["species_slug"] == slug}
        got = t["m1/datasets.tsv"]
        if sorted(r["dataset"] for r in got) != sorted(want) or len(got) != len(want):
            raise Stop(f"gate 4 datasets: {slug}: datasets.tsv holds other datasets than A's")
        sets = {}
        for r in got:
            lr = want[r["dataset"]]
            if (r["species"], r["antibiotic"], r["species_slug"]) != (lr["species"], lr["antibiotic"], slug):
                raise Stop(f"gate 4 datasets: {r['dataset']} is named otherwise than in lists.tsv")
            if r["set"] != ("M" if in_m(lr) else "primary"):
                raise Stop(f"gate 4 datasets: {r['dataset']} is in set {r['set']}, lists.tsv has in_M {lr['in_M']}")
            if in_m(lr) and r["labels_depend_on"] != lr["labels_depend_on"]:
                raise Stop(f"gate 4 datasets: {r['dataset']}: the label rule it waits for is not lists.tsv's")
            sets[r["dataset"]] = (r["set"], r["entered"] == "1")
        for name in ("m1/strata.tsv", "m1/pooled.tsv"):
            for r in t[name]:
                if name == "m1/strata.tsv":
                    s = sets.get(r["dataset"])
                    if s is None or not s[1]:
                        raise Stop(f"gate 4 datasets: {slug}: strata.tsv holds a dataset not entered")
                    if r["set"] != s[0]:
                        raise Stop(f"gate 4 datasets: {r['dataset']}: strata.tsv puts it in set {r['set']}, "
                                   f"datasets.tsv in {s[0]}")
                elif r["set"] not in ("primary", "M"):
                    raise Stop(f"gate 4 datasets: {slug}: pooled.tsv has a set that is neither primary nor M")
        pooled_sets = {(r["drug_class"], r["set"]) for r in t["m1/pooled.tsv"]}
        entered_sets = {(r["drug_class"], r["set"]) for r in got if r["entered"] == "1"}
        if not pooled_sets <= entered_sets:
            raise Stop(f"gate 4 datasets: {slug}: pooled.tsv pools a drug class and set with no entered dataset")


def floor_hits(path):
    """scan.py's table checks on one table: rule -> rows."""
    return scan.check_table(str(path))


def gate_floors(reading, species):
    for slug in species:
        for name in RESULTS:
            if not name.endswith(".tsv"):
                continue
            hits = floor_hits(reading / "runs" / slug / "results" / name)
            if hits:
                raise Stop(f"gate 5 floors: {slug}/{name}: {', '.join(sorted(hits))}")


def gate_combined(reading, species, contrasts):
    with tempfile.TemporaryDirectory() as tmp:
        argv = ["--contrasts", str(contrasts), "--out", tmp]
        for slug in species:
            argv += ["--run", str(reading / "runs" / slug)]
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            try:
                rc = combine.main(argv)
            except SystemExit as e:
                rc = e.code
        if rc != 0:
            raise Stop("gate 6 combined: combine.py does not run over runs/")
        for name in COMBINED:
            if (Path(tmp) / name).read_bytes() != (reading / "combined" / name).read_bytes():
                raise Stop(f"gate 6 combined: combined/{name} is not what combine.py writes from runs/")


def check_reading(reading, lrows, contrasts, commit):
    """Every gate of GATES, in order; the runs' tables. Raises Stop naming the gate."""
    reading = Path(reading)
    gate_layout(reading)
    species = gate_species(reading, lrows)
    gate_runs(reading, species, commit)
    runs = read_runs(reading, species)
    gate_datasets(runs, lrows)
    gate_floors(reading, species)
    gate_combined(reading, species, contrasts)
    return species, runs


# ---------- the page's tables ----------

def entered_strata(runs):
    """{(dataset, split, row, stratum): strata row} over every run."""
    return {(r["dataset"], r["split"], r["row"], r["stratum"]): r
            for t in runs.values() for r in t["m1/strata.tsv"]}


def num(x):
    return None if x in ("", None) else float(x)


def r12(x):
    return None if x is None else round(x + 0.0, 12)


def fmt(x):
    if x is None:
        return ""
    if isinstance(x, float):
        return f"{r12(x):.12f}"
    return str(x)


def page_datasets(runs, lrows):
    """Every dataset of A with its checks, species and antibiotic in full."""
    out = []
    for slug in sorted(runs):
        for r in sorted(runs[slug]["m1/datasets.tsv"], key=lambda r: r["dataset"]):
            row = dict(r)
            row["species"] = species_name(lrows, slug)
            out.append(row)
    return out


DATASETS_COLUMNS = ["dataset", "species", "antibiotic", "drug_class", "set", "labels_depend_on", "genomes",
                    "with_table", "rule_check", "hard_check", "learner_rows_equal", "learner_rows", "program_check",
                    "record_check", "entered", "why"]


def page_m1_datasets(runs, lrows):
    """Per entered dataset, split and row: the four strata's resistant genomes and recall, and recall on gene
    only less recall on point only where both hold 10 or more resistant genomes (combine.py's paired value)."""
    out = []
    at = entered_strata(runs)
    for slug in sorted(runs):
        for d in sorted({r["dataset"] for r in runs[slug]["m1/strata.tsv"]}):
            for split in ("lineage", "random"):
                for row in ROWS:
                    cells = {s: at.get((d, split, row, s)) for s in STRATA}
                    if any(c is None for c in cells.values()):
                        continue
                    first = cells["point_only"]
                    o = {"dataset": d, "species": species_name(lrows, slug), "antibiotic": first["antibiotic"],
                         "drug_class": first["drug_class"], "set": first["set"],
                         "labels_depend_on": first["labels_depend_on"], "split": split, "row": row}
                    for s in STRATA:
                        o[f"{s}_resistant"] = cells[s]["resistant"]
                        o[f"{s}_called"] = cells[s]["resistant_called"]
                        o[f"{s}_recall"] = cells[s]["recall"]
                    g, p = num(cells["gene_only"]["recall"]), num(cells["point_only"]["recall"])
                    o["gene_less_point"] = r12(g - p) if g is not None and p is not None else None
                    out.append(o)
    return out


M1_DATASETS_COLUMNS = (["dataset", "species", "antibiotic", "drug_class", "set", "labels_depend_on", "split", "row"]
                       + [f"{s}_{c}" for s in STRATA for c in ("resistant", "called", "recall")]
                       + ["gene_less_point"])


def page_m1_strata(runs):
    """Per set, split, row and stratum (and the counts beside): the entered datasets and their genome counts summed
    over datasets, a genome once per dataset. Counts only: no recall is pooled over datasets here."""
    sums = {}
    for slug in sorted(runs):
        for r in runs[slug]["m1/strata.tsv"]:
            key = (r["set"], r["split"], r["row"], r["stratum"])
            s = sums.setdefault(key, {"datasets": set(), "with_resistant": 0, "ten_or_more": 0, "resistant": 0,
                                      "resistant_called": 0, "susceptible": 0, "susceptible_called": 0,
                                      "kind": r["kind"]})
            s["datasets"].add(r["dataset"])
            n = int(r["resistant"])
            s["with_resistant"] += n > 0
            s["ten_or_more"] += n >= RECALL_FLOOR
            for k in ("resistant", "resistant_called", "susceptible", "susceptible_called"):
                s[k] += int(r[k])
    order = {s: i for i, s in enumerate(STRATA + BESIDE)}
    out = []
    for key in sorted(sums, key=lambda k: ({"primary": 0, "M": 1}[k[0]], {"lineage": 0, "random": 1}[k[1]],
                                           ROWS.index(k[2]), order.get(k[3], 99))):
        s = sums[key]
        out.append({"set": key[0], "split": key[1], "row": key[2], "stratum": key[3], "kind": s["kind"],
                    "datasets": len(s["datasets"]), "with_resistant": s["with_resistant"],
                    "ten_or_more": s["ten_or_more"], "resistant": s["resistant"],
                    "resistant_called": s["resistant_called"], "susceptible": s["susceptible"],
                    "susceptible_called": s["susceptible_called"]})
    return out


M1_STRATA_COLUMNS = ["set", "split", "row", "stratum", "kind", "datasets", "with_resistant", "ten_or_more",
                     "resistant", "resistant_called", "susceptible", "susceptible_called"]


def page_m1_pooled(runs, lrows):
    out = []
    for slug in sorted(runs):
        for r in runs[slug]["m1/pooled.tsv"]:
            out.append(dict(r, species=species_name(lrows, slug)))
    return out


M1_POOLED_COLUMNS = ["species", "drug_class", "set", "split", "row", "stratum", "datasets", "resistant_pairs",
                     "resistant_called", "susceptible_pairs", "susceptible_called", "recall"]


def reading_key(r):
    """A row of combine.py's m2-reading.tsv by cell, guard and statistic."""
    return (r["species_slug"], r["cell"], r["class"], r["dataset"], r["periods"], r["guard"], r["statistic"])


def m2_key(r):
    """A run's cell row as it is looked up in the reading, its guard included: a row of one guard never takes
    another guard's figures."""
    return (r["species_slug"], r["cell"], r["class"], r["dataset"], r["periods"], r["guard"])


def page_m2_cells(runs, combined_m2, lrows):
    """Every M2 row of the runs, with combine.py's Holm p and moves joined by cell and guard."""
    joined = {}
    for r in combined_m2:
        if reading_key(r) in joined:
            raise Stop("combined m2-reading.tsv holds a cell, guard and statistic twice")
        joined[reading_key(r)] = r
    out = []
    for slug in sorted(runs):
        for r in runs[slug]["m2/m2-cells.tsv"]:
            o = dict(r, species=species_name(lrows, slug))
            for stat in ("point", "tvd"):
                j = joined.get(m2_key(r) + (stat,))
                if j is None:
                    raise Stop(f"combined m2-reading.tsv has no {stat} row of a {slug} cell")
                o[f"holm_p_{stat}"] = j.get("holm_p", "")
                o[f"moves_{stat}"] = j.get("moves", "")
                o["decides"] = j["decides"]
            out.append(o)
    return out


M2_CELLS_COLUMNS = ["species", "cell", "class", "dataset", "antibiotic", "periods", "guard", "decides", "genomes",
                    "with_year_share", "earlier_genomes", "later_genomes", "earlier_units", "later_units", "passes",
                    "earlier_point", "later_point", "d_point", "p_point", "holm_p_point", "moves_point",
                    "earlier_point_disrupt", "later_point_disrupt", "earlier_gene", "later_gene", "d_gene", "tvd",
                    "p_tvd", "holm_p_tvd", "moves_tvd", "note"]


def concat(runs, lrows, name):
    out = []
    for slug in sorted(runs):
        for r in runs[slug][name]:
            out.append(dict(r, species=species_name(lrows, slug)))
    return out


def deciding(rows):
    return [r for r in rows if (r["cell"], r["periods"], r["guard"]) == DECIDING]


def page_species(runs, lrows, cells):
    out = []
    for slug in sorted(runs):
        t = runs[slug]
        ds = t["m1/datasets.tsv"]
        mine = [r for r in cells if r["species_slug"] == slug]
        dec = deciding(mine)
        sec = [r for r in mine if (r["cell"], r["periods"], r["guard"]) == ("species_class", "secondary",
                                                                            "bioproject")]
        lineage = {}
        for r in t["m1/strata.tsv"]:
            if r["split"] == "lineage" and r["row"] == "program" and r["kind"] == "stratum":
                lineage.setdefault(r["dataset"], {})[r["stratum"]] = int(r["resistant"])
        both_ten = [d for d, s in lineage.items()
                    if s.get("gene_only", 0) >= RECALL_FLOOR and s.get("point_only", 0) >= RECALL_FLOOR]
        sets = {r["dataset"]: r["set"] for r in ds}
        years = {r["years"]: int(r["genomes"]) for r in t["m2/m2-years.tsv"]}
        genomes = sum(v for k, v in years.items() if k != "no table")
        out.append({
            "species": species_name(lrows, slug), "species_slug": slug,
            "datasets": len(ds), "in_m": sum(r["set"] == "M" for r in ds),
            "entered": sum(r["entered"] == "1" for r in ds),
            "hard_check_pass": sum(r["hard_check"] == "1" for r in ds),
            "learner_rows_equal": sum(int(r["learner_rows_equal"]) for r in ds),
            "learner_rows": sum(int(r["learner_rows"]) for r in ds),
            "m1_primary_both_ten": sum(sets[d] == "primary" for d in both_ten),
            "m1_m_both_ten": sum(sets[d] == "M" for d in both_ten),
            "m2_genomes": genomes, "m2_no_table": years.get("no table", 0),
            "m2_genomes_2010_2014": years.get("2010 to 2014", 0), "m2_genomes_2015_2019": years.get("2015 to 2019", 0),
            "m2_classes": len(dec), "m2_pass_primary": sum(r["passes"] == "1" for r in dec),
            "m2_pass_secondary": sum(r["passes"] == "1" for r in sec),
            "m2_moves_point": sum(r["moves_point"] == "1" for r in dec),
            "m2_moves_tvd": sum(r["moves_tvd"] == "1" for r in dec)})
    return out


SPECIES_COLUMNS = ["species", "species_slug", "datasets", "in_m", "entered", "hard_check_pass", "learner_rows_equal",
                   "learner_rows", "m1_primary_both_ten", "m1_m_both_ten", "m2_genomes", "m2_no_table",
                   "m2_genomes_2010_2014", "m2_genomes_2015_2019", "m2_classes", "m2_pass_primary",
                   "m2_pass_secondary", "m2_moves_point", "m2_moves_tvd"]
FAMILIES_COLUMNS = ["species", "class", "periods", "family", "earlier_genomes", "later_genomes", "earlier_mix",
                    "later_mix"]
YEARS_COLUMNS = ["species", "years", "genomes"]
LINEAGE_COLUMNS = ["species", "reports", "report_pairs", "differing_genome_pairs", "given"]
TIMESPLIT_COLUMNS = ["dataset", "species", "antibiotic", "reading", "resistant_later_test", "with_determinant",
                     "with_unseen_determinant", "no_determinant", "share_unseen"]


# ---------- writing ----------

README = """# The mechanism split's reading of 23 November 2026: the record

This directory holds the tables behind [docs/mechanism-split-result.md](../../docs/mechanism-split-result.md): blocks M1 (recall by mechanism) and M2 (the mechanism mix over collection years) of the mechanism split, read once over the runs of the {species} species of A, {datasets} datasets, {primary} of them in the primary set and {m} in M. Every row is an aggregate over a dataset, a species, a cell or a learner. No file holds a genome id, a BioProject accession, a sequence or a position.

`scripts/mechanism/combine.py` wrote the five files under `combined/`; they are copied here as it wrote them. `scripts/mechanism/page_tables.py write` checked the reading against the six gates of `gates.txt` and wrote every other file from the runs' own tables, species and antibiotic in full, with no new measure. `reading.txt` lists the reading's commands with the sha256 of what each wrote, and `checks.txt` the checks of the run, counts only. The rules of the reading are `docs/mechanism-split.md`'s, fixed before any run, with the rulings in `docs/assumption-plan-corrections.md` (entries "Mechanism design" and "Mechanism split, first run").

| File | Holds |
|---|---|
| `gates.txt` | the six gates the reading passed before a table was written |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `checks.txt` | the checks of the run, counts only |
| `combined/m1-reading.tsv` | per split, set and row: recall on gene only less recall on point only over the datasets where both hold 10 or more resistant genomes; the lineage split of the primary set decides for the five learners |
| `combined/m1-neither.tsv` | per split, set and row: the median recall on the neither stratum |
| `combined/m2-reading.tsv` | per M2 row and statistic: the value, its permutation p, Holm over the passing cells and whether it moves |
| `combined/m2-summary.tsv` | per cell kind, period pair, guard and statistic: cells, passing, moving |
| `combined/timesplit-link.tsv` | per learner: Spearman's rank correlation of the unseen-determinant share with the past-less-random difference (exploratory) |
| `page-species.tsv` | per species, named in full: datasets, checks, M1 datasets with both strata at 10 or more, M2 genomes, cells, passing and moving |
| `page-datasets.tsv` | per dataset of A, species and antibiotic in full: the checks and whether it entered |
| `page-m1-datasets.tsv` | per entered dataset, split and row: resistant genomes, called and recall per stratum, and gene only less point only |
| `page-m1-strata.tsv` | per set, split, row and stratum: genome counts summed over the entered datasets, a genome once per dataset; no pooled recall |
| `page-m1-pooled.tsv` | strata under 10 resistant genomes pooled per species, drug class and set, as `strata.py` wrote them |
| `page-m2-cells.tsv` | every M2 row of the runs, species in full, with combine.py's Holm p and moves |
| `page-m2-families.tsv` | the family mix of each passing cell, the catalogue's family names over 5 or more genomes, the rest as `other families` |
| `page-m2-years.tsv` | per species, genomes per span of collection years |
| `page-m2-lineage.tsv` | per species, whether the reports agree on the lineage clusters |
| `page-timesplit.tsv` | per time-split dataset: resistant later test genomes with a determinant no resistant training genome carries |

A share is no prevalence: the genomes are BV-BRC deposits chosen for having lab results. A trend describes deposits. A table names no cause and no new marker. A split by lineage or study is a retrospective split inside one database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/mechanism/page_tables.py print` prints the page's tables from these files.
"""


def write_tsv(path, rows, columns):
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            cells = [fmt(r.get(c)) for c in columns]
            if any("\t" in c or "\n" in c for c in cells):
                raise Stop(f"{path}: a cell holds a tab or a newline")
            handle.write("\t".join(cells) + "\n")


def write(reading, l_path, contrasts, commit, out):
    reading, out = Path(reading), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    species, runs = check_reading(reading, lrows, contrasts, commit)
    combined_m2 = load(reading / "combined" / "m2-reading.tsv")
    cells = page_m2_cells(runs, combined_m2, lrows)
    out.mkdir(parents=True, exist_ok=True)
    (out / "combined").mkdir()
    for name in COMBINED:
        (out / "combined" / name).write_bytes((reading / "combined" / name).read_bytes())
    for name in ("reading.txt", "checks.txt"):
        (out / name).write_bytes((reading / name).read_bytes())
    write_tsv(out / "page-species.tsv", page_species(runs, lrows, cells), SPECIES_COLUMNS)
    write_tsv(out / "page-datasets.tsv", page_datasets(runs, lrows), DATASETS_COLUMNS)
    write_tsv(out / "page-m1-datasets.tsv", page_m1_datasets(runs, lrows), M1_DATASETS_COLUMNS)
    write_tsv(out / "page-m1-strata.tsv", page_m1_strata(runs), M1_STRATA_COLUMNS)
    write_tsv(out / "page-m1-pooled.tsv", page_m1_pooled(runs, lrows), M1_POOLED_COLUMNS)
    write_tsv(out / "page-m2-cells.tsv", cells, M2_CELLS_COLUMNS)
    write_tsv(out / "page-m2-families.tsv", concat(runs, lrows, "m2/m2-families.tsv"), FAMILIES_COLUMNS)
    write_tsv(out / "page-m2-years.tsv", concat(runs, lrows, "m2/m2-years.tsv"), YEARS_COLUMNS)
    write_tsv(out / "page-m2-lineage.tsv", concat(runs, lrows, "m2/m2-lineage.tsv"), LINEAGE_COLUMNS)
    write_tsv(out / "page-timesplit.tsv", concat(runs, lrows, "m2/timesplit.tsv"), TIMESPLIT_COLUMNS)
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    a = a_rows(lrows)
    (out / "README.md").write_text(README.format(species=len(species), datasets=len(a),
                                                 primary=sum(not in_m(r) for r in a), m=sum(in_m(r) for r in a)),
                                   encoding="utf-8")
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w", encoding="utf-8") as handle:
        for name in files:
            handle.write(f"{sha256_of(out / name)}  {name}\n")
    return files


# ---------- printing ----------

def d3(x):
    if x in ("", None):
        return "none"
    value = Decimal(str(x)).quantize(Decimal("0.001"), rounding=ROUND_HALF_EVEN)
    if value == 0:
        value = Decimal("0.000")
    return str(value)


def pval(x):
    if x in ("", None):
        return "none"
    return "< 0.001" if float(x) < 0.001 else d3(x)


def count_words(x):
    return f"{int(x):,}"


def of(part, whole):
    return f"{count_words(part)} of {count_words(whole)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def print_m1_reading(rows):
    """The five learners on the lineage split of the primary set, then the two catalogue callers beside."""
    out = []
    for r in rows:
        if (r["split"], r["set"]) != ("lineage", "primary"):
            continue
        out.append([LEARNER_WORDS[r["row"]], r["datasets"], r["species"], d3(r["median"]), interval(r["lo"], r["hi"]),
                    r["gene_ahead"], r["point_ahead"], r["same"], pval(r["sign_p"]),
                    pval(r["holm_p"]) if r["decides"] == "1" else "beside",
                    r["reading"] if r["decides"] == "1" else "decides nothing"])
    return table(["Row", "Datasets", "Species", "Median gene less point", "Species interval", "Gene ahead",
                  "Point ahead", "Same", "Sign p", "Holm p", "Reading"], out)


def print_m1_beside(rows):
    """Every other split and set: beside, deciding nothing."""
    out = []
    for r in rows:
        if (r["split"], r["set"]) == ("lineage", "primary"):
            continue
        out.append([r["split"] + (f" ({r['label']})" if r["label"] else ""), r["set"], LEARNER_WORDS[r["row"]],
                    r["datasets"], r["species"], d3(r["median"]), interval(r["lo"], r["hi"]), r["gene_ahead"],
                    r["point_ahead"], r["same"], pval(r["sign_p"])])
    return table(["Split", "Set", "Row", "Datasets", "Species", "Median gene less point", "Species interval",
                  "Gene ahead", "Point ahead", "Same", "Sign p"], out)


def print_m1_neither(rows):
    out = [[r["split"], r["set"], LEARNER_WORDS[r["row"]], r["datasets"], d3(r["median_recall"])]
           for r in rows if r["split"] == "lineage"]
    return table(["Split", "Set", "Row", "Datasets", "Median recall on neither"], out)


def print_m1_strata(rows):
    """Resistant genomes per stratum on the lineage split, summed over the entered datasets (row-free counts, read
    off the program's row)."""
    out = []
    for r in rows:
        if r["split"] != "lineage" or r["row"] != "program":
            continue
        out.append([r["set"], STRATUM_WORDS.get(r["stratum"], r["stratum"].replace("_", " ")),
                    r["datasets"], r["with_resistant"], r["ten_or_more"], count_words(r["resistant"]),
                    count_words(r["susceptible"])])
    return table(["Set", "Stratum", "Datasets", "With a resistant genome", "With 10 or more resistant",
                  "Resistant genomes", "Susceptible genomes"], out)


def print_m2_summary(rows):
    out = []
    for r in rows:
        if r["cell"] != "species_class":
            continue
        out.append([PERIOD_WORDS[r["periods"]], GUARD_WORDS[r["guard"]], r["statistic"], r["cells"], r["pass"],
                    r["move"], "yes" if r["decides"] == "1" else "no"])
    return table(["Periods", "Guard", "Statistic", "Cells", "Pass", "Move", "Decides"], out)


def print_m2_secondary(rows):
    out = [[PERIOD_WORDS[r["periods"]], GUARD_WORDS[r["guard"]], r["statistic"], r["cells"], r["pass"], r["move"]]
           for r in rows if r["cell"] == "dataset_resistant"]
    return table(["Periods", "Guard", "Statistic", "Dataset cells", "Pass", "Move"], out)


def print_m2_cells(rows):
    """The deciding cells that pass the BioProject vote."""
    out = []
    for r in rows:
        if (r["cell"], r["periods"], r["guard"]) != DECIDING or r["passes"] != "1":
            continue
        out.append([f"*{r['species']}*", r["class"], f"{r['earlier_units']} / {r['later_units']}",
                    f"{count_words(r['earlier_genomes'])} / {count_words(r['later_genomes'])}",
                    f"{d3(r['earlier_point'])} to {d3(r['later_point'])}", pval(r["holm_p_point"]),
                    f"{d3(r['earlier_gene'])} to {d3(r['later_gene'])}", d3(r["tvd"]), pval(r["holm_p_tvd"]),
                    "point" * (r["moves_point"] == "1") + (" and " if r["moves_point"] == r["moves_tvd"] == "1" else "")
                    + "mix" * (r["moves_tvd"] == "1") or "no"])
    return table(["Species", "Class", "BioProjects voting", "Genomes", "Point share", "Holm p",
                  "Gene share", "Mix distance", "Holm p", "Moves"], out)


def print_species(rows):
    out = []
    for r in rows:
        out.append([f"*{r['species']}*", r["datasets"], r["in_m"], r["entered"],
                    of(r["learner_rows_equal"], r["learner_rows"]), r["m1_primary_both_ten"],
                    count_words(r["m2_genomes"]), r["m2_classes"], r["m2_pass_primary"], r["m2_pass_secondary"],
                    r["m2_moves_point"], r["m2_moves_tvd"]])
    return table(["Species", "Datasets", "In M", "Entered", "Rows equal", "M1 primary, both 10 or more",
                  "M2 genomes", "Classes", "Pass, primary", "Pass, secondary", "Point moves", "Mix moves"], out)


def print_datasets(rows, m1):
    lineage = {(r["dataset"], r["row"]): r for r in m1 if r["split"] == "lineage"}
    out = []
    for r in rows:
        prog = lineage.get((r["dataset"], "program"))
        diffs = [d3(lineage[(r["dataset"], x)]["gene_less_point"]) if (r["dataset"], x) in lineage else "none"
                 for x in LEARNERS]
        out.append([f"*{r['species']}*", r["antibiotic"], r["set"], "yes" if r["entered"] == "1" else "no",
                    prog["gene_only_resistant"] if prog else "none", prog["point_only_resistant"] if prog else "none",
                    prog["both_resistant"] if prog else "none", prog["neither_resistant"] if prog else "none"]
                   + diffs)
    return table(["Species", "Antibiotic", "Set", "Entered", "Gene only", "Point only", "Both", "Neither"]
                 + [f"`{x}`" for x in LEARNERS], out)


def tables(folder):
    folder = Path(folder)
    m1 = load(folder / "combined" / "m1-reading.tsv")
    summary = load(folder / "combined" / "m2-summary.tsv")
    return {
        "m1_reading": print_m1_reading(m1),
        "m1_beside": print_m1_beside(m1),
        "m1_neither": print_m1_neither(load(folder / "combined" / "m1-neither.tsv")),
        "m1_strata": print_m1_strata(load(folder / "page-m1-strata.tsv")),
        "m2_summary": print_m2_summary(summary),
        "m2_secondary": print_m2_secondary(summary),
        "m2_cells": print_m2_cells(load(folder / "page-m2-cells.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "datasets": print_datasets(load(folder / "page-datasets.tsv"), load(folder / "page-m1-datasets.tsv")),
    }


def main(argv):
    if len(argv) == 6 and argv[0] == "write":
        try:
            files = write(*argv[1:])
        except Stop as e:
            print(f"stopped: {e}", file=sys.stderr)
            return 1
        print(f"gates: {len(GATES)} of {len(GATES)} pass; wrote {len(files) + 1} files")
        return 0
    if len(argv) == 2 and argv[0] == "print":
        for name, lines in tables(argv[1]).items():
            print(f"[{name}]")
            print("\n".join(lines))
            print()
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
