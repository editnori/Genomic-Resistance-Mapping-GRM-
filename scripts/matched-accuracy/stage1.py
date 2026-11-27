"""Stage 1 of the matched accuracy protocol (docs/matched-accuracy.md): count the eligible rows and freeze them.

  python3 -I stage1.py --reports DIR --records FILE --drug-names FILE SHA256 --labels FILE SHA256
                       --genomes FILE --hc20 FILE SHA256 --out DIR
  python3 -I stage1.py ... --hc20 - --out DIR

Stage 1 alone. It reads no score: the report loader keeps a named list of
keys and, of each learner's list of scores, its length alone. Stage 2 is
not written here.

Inputs. DIR holds the release's reports, `<dataset>.json`, each equal to
its line of the release's own sums list, committed as
reference/lineage-control/rev14-274-20261112.reports.sha256. --records is
the lab records table cut at the seal (bvbrc/genome_amr.jsonl), --genomes
the second copy of BV-BRC's genome table (gzip), which gives the
BioProject. --drug-names and --labels are curated/drug_names.json and the
pinned label table curated/phenotypes.jsonl, each with its sha256. --hc20
is the 5 Nov genome table with its sha256, or `-` for no HC20 read. Every
sha256 is checked before any other read.

Steps.
1. Check each sha256, then each report: the lineage split's 5 folds,
   test_fold and every scores list as long as the genome rows, each
   lineage cluster and each copy group in one fold. A report whose lineage
   split holds none of scm, cart, rf, gbt and lasso is counted and left out.
2. Check the lab records: 1,285,111 lines, each inserted before
   2026-11-05T19:10:00 (the first 19 characters of each text of its
   insertion time, a real date and time, compared as text). The 252,148
   records dropped at the cut are a recorded figure, never measured here.
3. Check each row's lab call against the label table, for rows with
   from_mic false. Rows with from_mic true are counted apart.
4. Take each row out under the first reason that applies: 1 a copy group,
   2 no single BioProject (a value of one part without an accession's shape,
   PRJ, two letters and digits, reads as none and is counted apart), 3 no lab
   record for the pair, 4 a record with a blank method or standard, 5 more
   than one tuple of method, standard and year. Group the rest by dataset, fold, lineage cluster, BioProject,
   method, standard and year; a group with one lab outcome leaves (6).
5. Reading A is those rows. Reading B first takes out every pair with an
   intermediate record, then forms the groups again.
6. Write one row list a reading under <out>/local/, mode 0600, and counts
   alone under <out>/export/, after a scan for every genome id and
   accession-shaped BioProject part the run read, with a positive control.

Every check that fails prints one line and exits 1, writing nothing; a scan
hit prints one line a file.
Python 3.12, standard library only; run with -I.
"""

import argparse
import collections
import datetime
import fractions
import gzip
import hashlib
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SEAL = "2026-11-05T19:10:00"
STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d", re.ASCII)
SHA = re.compile(r"[0-9a-f]{64}")
RECORDS_SHA256 = "d231fa768a31fc4ce4c8a6bb7024d05d632621846070688adb892fbfe340fd78"
RECORDS_LINES = 1285111
DROPPED_RECORDED = 252148
GENOMES_SHA256 = "84ae57a777cb234aa0a39102ab7f9ca9448e52d1474654179f26165644c25f64"
REPORT_SUMS = ROOT / "reference" / "lineage-control" / "rev14-274-20261112.reports.sha256"
REPORT_SUMS_SHA256 = "0aedf8780ca8e2292cfc5f227f61b87b7c6b006a6649ed4c0378514c70b72223"

FIVE = ("scm", "cart", "rf", "gbt", "lasso")
LEARNER_NAMES = FIVE + ("catalog",)
FOLDS = 5
# The loader's named keys. Every other key of a report is dropped before it returns.
REPORT_KEYS = ("id", "species", "species_slug", "antibiotic", "kset_duplicates")
ROW_KEYS = ("id", "resistant", "from_mic", "lineage", "mlst")
SPLIT_KEYS = ("folds", "test_fold", "unavailable_reason")

R_WORDS = {"resistant", "nonsusceptible", "non-susceptible"}
S_WORDS = {"susceptible"}
I_WORDS = {"intermediate", "susceptible-dose dependent", "reduced susceptibility"}

REASONS = ("copy_group", "bioproject", "no_lab_record", "blank_method_or_standard", "assay_not_single")
COUNT_COLUMNS = ("level", "name", "reading", "reports", "reports_left_out", "genome_rows", "from_mic_rows",
                 "out_copy_group", "out_bioproject", "out_bioproject_not_in_table", "out_bioproject_no_accession",
                 "out_no_lab_record", "out_blank_method_or_standard", "out_assay_not_single", "out_intermediate",
                 "out_one_lab_outcome", "rows", "genomes", "groups", "pairs", "datasets", "species")
OUT_COLUMNS = tuple(c for c in COUNT_COLUMNS if c.startswith("out_"))
TYPING_COLUMNS = ("level", "name", "reading", "step", "lacking_value_of_rows_entering_groups",
                  "lacking_value_of_eligible_rows", "rows", "genomes", "groups", "pairs", "datasets", "species",
                  "primary_groups", "primary_pairs", "mlst_groups", "mlst_pairs", "rows_lost_from_primary")
SHARE_COLUMNS = ("reading", "bioprojects", "eligible_rows", "largest_bioproject_rows",
                 "its_weight_numerator", "its_weight_denominator",
                 "largest_weight_numerator", "largest_weight_denominator")
ROW_LIST = {"A": "reading-A.tsv", "B": "reading-B.tsv"}
ID_SHAPED = re.compile(r"\d+\.\d+")
ACCESSION_SHAPED = re.compile(r"\b(?:PRJ[A-Z]{1,2}|SAM[A-Z]{1,2}|GC[AF]_)\d+")
# A BioProject part is an accession when the whole part is PRJ, two letters and digits, matched without case.
BIOPROJECT = re.compile(r"PRJ[A-Z]{2}[0-9]+", re.IGNORECASE | re.ASCII)


class Stop(Exception):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def value_text(v):
    """A record field as grm reads it (src/bvbrc.rs, value_text): trimmed text, a list joined with ';', None if blank."""
    if v is None:
        return None
    if isinstance(v, str):
        return v.strip() or None
    if isinstance(v, list):
        parts = [p for p in (value_text(x) for x in v) if p is not None]
        return ";".join(parts) if parts else None
    if isinstance(v, bool):
        return "true" if v else "false"
    return repr(v) if isinstance(v, float) else str(v)


def values(v):
    """The distinct texts of a field: as value_text, but a list gives each of its texts apart."""
    if isinstance(v, list):
        out = set()
        for x in v:
            out |= values(x)
        return out
    t = value_text(v)
    return {t} if t is not None else set()


def parse_call(text):
    """src/curate.rs, parse_call: R, S, I or None, on the trimmed lowercase text."""
    if text is None:
        return None
    t = text.strip().lower()
    return "R" if t in R_WORDS else "S" if t in S_WORDS else "I" if t in I_WORDS else None


def stamp_heads(v):
    """The first 19 characters of each text of an insertion time, a list giving each of its texts apart, or None
    when the value is blank or a text of it is not a real date and time to the second in ASCII digits."""
    if isinstance(v, list):
        heads = [stamp_heads(x) for x in v]
        return None if not heads or None in heads else [h for hs in heads for h in hs]
    t = value_text(v)
    head = t[:19] if t else ""
    if not STAMP.fullmatch(head):
        return None
    try:
        datetime.datetime.strptime(head, "%Y-%m-%dT%H:%M:%S")
    except ValueError:
        return None
    return [head]


def jsonl(path, what):
    """Each line's JSON object, gzip or plain. A line that is not a JSON object, or a gzip file that does not read
    to its end, stops the run."""
    with open(path, "rb") as f:
        packed = f.read(2) == b"\x1f\x8b"
    opener = gzip.open if packed else open
    with opener(path, "rb") as f:
        n = 0
        while True:
            try:
                line = f.readline()
            except (OSError, EOFError):
                raise Stop(f"{what} is not a whole gzip file") from None
            if not line:
                return
            n += 1
            try:
                row = json.loads(line)
            except ValueError:
                raise Stop(f"line {n} of {what} is not JSON") from None
            if not isinstance(row, dict):
                raise Stop(f"line {n} of {what} is not a JSON object")
            yield row


# 1. The sha256 of every input, before any other read.

def report_files(directory):
    """[(dataset, path)] of the reports, each checked against its line of the release's sums list."""
    if sha256(REPORT_SUMS) != REPORT_SUMS_SHA256:
        raise Stop("the release's sums list of reports differs from its sha256")
    want = {}
    for line in REPORT_SUMS.read_text(encoding="utf-8").splitlines():
        digest, _, name = line.partition("  ")
        if not SHA.fullmatch(digest) or not name.startswith("./datasets/") or not name.endswith(".json"):
            raise Stop("a line of the release's sums list of reports is not a sha256 and a report name")
        want[name[len("./datasets/"):]] = digest
    directory = Path(directory)
    if not directory.is_dir():
        raise Stop("the reports directory is missing")
    have = sorted(p.name for p in directory.iterdir() if p.name.endswith(".json"))
    if sorted(want) != have:
        raise Stop(f"the reports directory holds {len(have)} reports; the release's sums list names {len(want)} "
                   f"and {len(set(have) ^ set(want))} differ")
    out = []
    for name in have:
        if sha256(directory / name) != want[name]:
            raise Stop(f"{name} differs from its line of the release's sums list")
        out.append((name[:-len(".json")], directory / name))
    return out


def check_sums(args):
    """{input: sha256}, each checked against its pin or the sha256 given beside it."""
    pins = [("lab records", args.records, RECORDS_SHA256), ("genome table, second copy", args.genomes, GENOMES_SHA256),
            ("drug names", args.drug_names[0], args.drug_names[1]), ("label table", args.labels[0], args.labels[1])]
    if args.hc20 != ["-"]:
        pins.append(("genome table, 5 Nov", args.hc20[0], args.hc20[1]))
    out = {}
    for what, path, want in pins:
        if not Path(path).is_file():
            raise Stop(f"the {what} file is missing")
        got = sha256(path)
        if got != want:
            raise Stop(f"the {what} file differs from its sha256")
        out[what] = got
    return out


# 2. The report loader, and the checks of each report.

def load_report(path):
    """The report's named keys alone. Of each learner of the lineage split: its name and its count of scores."""
    try:
        with open(path, encoding="utf-8") as f:
            doc = json.load(f)
    except ValueError:
        raise Stop(f"{Path(path).name} is not a JSON report") from None
    try:
        return project(doc)
    finally:
        del doc


def project(doc):
    if not isinstance(doc, dict) or not isinstance(doc.get("samples"), list) or not isinstance(doc.get("splits"), list):
        raise Stop("a report has no list of samples or of splits")
    lineage = [s for s in doc["splits"] if isinstance(s, dict) and s.get("scheme") == "lineage"]
    if len(lineage) != 1:
        raise Stop(f"{doc.get('id')}: the report holds {len(lineage)} lineage splits")
    split = lineage[0]
    learners = split.get("learners")
    if not isinstance(learners, list) or not all(isinstance(e, dict) for e in learners):
        raise Stop(f"{doc.get('id')}: the lineage split's learners are not a list")
    out = {k: doc.get(k) for k in REPORT_KEYS}
    out["kset_duplicates"] = [list(g) if isinstance(g, list) else g for g in doc["kset_duplicates"]] \
        if isinstance(doc.get("kset_duplicates"), list) else doc.get("kset_duplicates")
    out["samples"] = [{k: s.get(k) for k in ROW_KEYS} if isinstance(s, dict) else None for s in doc["samples"]]
    kept = {k: split.get(k) for k in SPLIT_KEYS}
    if isinstance(kept["test_fold"], list):
        kept["test_fold"] = list(kept["test_fold"])
    kept["learners"] = [{"learner": e.get("learner"),
                         "scores_length": len(e["scores"]) if isinstance(e.get("scores"), list) else None}
                        for e in learners]
    out["lineage_split"] = kept
    return out


def is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def enters(rep):
    """True when the lineage split holds the five learners, False when it holds none of them."""
    rid = rep["id"]
    names = [e["learner"] for e in rep["lineage_split"]["learners"]]
    if len(names) != len(set(names)) or not set(names) <= set(LEARNER_NAMES):
        raise Stop(f"{rid}: the lineage split's learners are not each named once from {', '.join(LEARNER_NAMES)}")
    five = set(names) & set(FIVE)
    if five and five != set(FIVE):
        raise Stop(f"{rid}: the lineage split holds {len(five)} of the five learners")
    return bool(five)


def check_report(rep):
    """The checks that stop the run, on a report that enters."""
    rid, rows, split = rep["id"], rep["samples"], rep["lineage_split"]
    n = len(rows)
    for s in rows:
        if s is None or not isinstance(s["id"], str) or not s["id"].strip() or not isinstance(s["resistant"], bool) \
                or not isinstance(s["from_mic"], bool) or not is_int(s["lineage"]) \
                or not (s["mlst"] is None or isinstance(s["mlst"], str)):
            raise Stop(f"{rid}: a genome row lacks a text id, a true or false call and from_mic, or an integer lineage")
    if len({s["id"] for s in rows}) != n:
        raise Stop(f"{rid}: a genome id is repeated in the report")
    fold = split["test_fold"]
    if not isinstance(fold, list) or len(fold) != n:
        raise Stop(f"{rid}: test_fold is not a list as long as the report's {n} genome rows")
    for e in split["learners"]:
        if e["scores_length"] != n:
            raise Stop(f"{rid}: the scores of {e['learner']} are not a list as long as the report's {n} genome rows")
    if split["folds"] != FOLDS or not all(is_int(f) for f in fold) or set(fold) != set(range(FOLDS)):
        raise Stop(f"{rid}: the lineage split has other than {FOLDS} folds")
    seen = {}
    for s, f in zip(rows, fold):
        if seen.setdefault(s["lineage"], f) != f:
            raise Stop(f"{rid}: a lineage cluster sits in two folds")
    groups = rep["kset_duplicates"] or []
    if not isinstance(groups, list) or not all(isinstance(g, list) and g and all(is_int(i) and 0 <= i < n for i in g)
                                               for g in groups):
        raise Stop(f"{rid}: kset_duplicates is not a list of lists of sample indices")
    members = [i for g in groups for i in g]
    if len(members) != len(set(members)):
        raise Stop(f"{rid}: a genome sits in two copy groups")
    for g in groups:
        if len({fold[i] for i in g}) != 1:
            raise Stop(f"{rid}: a copy group sits in two folds")


def load_reports(files):
    """(reports that enter, reports left out), each sorted by dataset."""
    entering, left_out, ids = [], [], set()
    for name, path in files:
        rep = load_report(path)
        for k in ("id", "species", "species_slug", "antibiotic"):
            if not isinstance(rep[k], str) or not rep[k].strip():
                raise Stop(f"{name}: the report's {k} is not text")
        if rep["id"] in ids:
            raise Stop(f"{rep['id']}: two reports give this dataset")
        ids.add(rep["id"])
        if enters(rep):
            check_report(rep)
            entering.append(rep)
        else:
            if not all(s is not None and isinstance(s["id"], str) and s["id"].strip() for s in rep["samples"]):
                raise Stop(f"{rep['id']}: a genome row of a report left out lacks a text id")
            left_out.append(rep)
    return entering, left_out


# 3. The lab records, the label table and the genome tables.

def drug_names(path):
    """raw name -> (kind, canonical), as m_list.py reads curated/drug_names.json."""
    try:
        with open(path, encoding="utf-8") as f:
            doc = json.load(f)
    except ValueError:
        raise Stop("the drug names file is not JSON") from None
    names = doc.get("names") if isinstance(doc, dict) else None
    if not isinstance(names, list) or not all(
            isinstance(n, dict) and isinstance(n.get("raw"), str) and isinstance(n.get("kind"), str)
            and (n.get("canonical") is None or isinstance(n.get("canonical"), str)) for n in names):
        raise Stop("the drug names file is not a list of names, each with a raw text, a kind and a canonical name")
    return {n["raw"]: (n["kind"], n["canonical"]) for n in names}


def read_records(path, resolve, wanted):
    """({(genome, antibiotic): record facts} for the wanted pairs, every genome id read, line count).

    A record's name reaches its antibiotic through drug_names.json, as in m_list.py. Facts of a pair: its tuples
    of method, uppercased standard and year (None for a blank year), whether a record has a blank method or
    standard, and whether a record's call is intermediate. Every line must be inserted before the seal."""
    pairs, ids, n = {}, set(), 0
    for rec in jsonl(path, "the lab records"):
        n += 1
        heads = stamp_heads(rec.get("date_inserted"))
        if heads is None:
            raise Stop("a lab record's date_inserted is blank or not a time to the second")
        if max(heads) >= SEAL:
            raise Stop(f"a lab record was inserted at or after {SEAL}")
        raw = value_text(rec.get("antibiotic"))
        if raw is None:
            raise Stop("a lab record has no antibiotic name")
        kind, canonical = resolve.get(raw, (None, None))
        if kind is None:
            raise Stop("a lab record's antibiotic name is not in the drug names file")
        genome = value_text(rec.get("genome_id"))
        if genome is not None:
            ids.add(genome)
        if kind != "drug" or (genome, canonical) not in wanted:
            continue
        p = pairs.setdefault((genome, canonical), {"tuples": set(), "blank": False, "intermediate": False})
        methods = values(rec.get("laboratory_typing_method"))
        standards = {s.upper() for s in values(rec.get("testing_standard"))}
        years = values(rec.get("testing_standard_year")) or {None}
        if not methods or not standards:
            p["blank"] = True
        p["tuples"] |= {(m, s, y) for m in methods for s in standards for y in years}
        if parse_call(value_text(rec.get("resistant_phenotype"))) == "I":
            p["intermediate"] = True
    if n != RECORDS_LINES:
        raise Stop(f"the lab records hold {n} lines, not {RECORDS_LINES}")
    return pairs, ids, n


def read_labels(path, wanted):
    """({(genome, antibiotic): the table's call or None} for the wanted pairs, every genome id read)."""
    table, ids = {}, set()
    for p in jsonl(path, "the label table"):
        genome = value_text(p.get("genome_id"))
        if genome is not None:
            ids.add(genome)
        antibiotic = p.get("antibiotic")
        if antibiotic is not None and not isinstance(antibiotic, str):
            raise Stop("a row of the label table names an antibiotic that is not text")
        key = (genome, antibiotic)
        if key not in wanted:
            continue
        if key in table:
            raise Stop("the label table holds a genome and antibiotic pair twice")
        table[key] = None if p.get("conflict") or p.get("label_from_mic") else p.get("reported")
    return table, ids


def check_labels(entering, table):
    """(rows that equal the table, rows labelled from an MIC). A row with from_mic false must equal the table."""
    equal = mic = 0
    for rep in entering:
        for s in rep["samples"]:
            if s["from_mic"]:
                mic += 1
                continue
            if table.get((s["id"], rep["antibiotic"])) != ("R" if s["resistant"] else "S"):
                raise Stop(f"{rep['id']}: a row's lab call differs from the label table's call for its pair")
            equal += 1
    return equal, mic


def bioproject_parts(v):
    return {p for t in values(v) for p in re.split(r"[\s,;]+", t) if p}


def is_accession(part):
    return BIOPROJECT.fullmatch(part) is not None


def read_genomes(path, wanted):
    """({genome: BioProject or None}, {genome: why it has none}, ids read, accession-shaped parts read, rows).

    A row whose BioProject is blank is skipped, as the study split's reader skips it
    (scripts/c2-cap/check_folds.py, bioprojects). A genome gets None when the copy lacks it ("not_in_table"),
    when its one part has no accession's shape ("no_accession"), or when it has no row with a value, a value of
    more than one part, or two rows that differ. Every row must be inserted before the seal."""
    found, ids, accessions, n = {}, set(), set(), 0
    for row in jsonl(path, "the genome table, second copy"):
        n += 1
        genome = value_text(row.get("genome_id"))
        if genome is None:
            raise Stop("a row of the genome table, second copy, has no genome_id")
        heads = stamp_heads(row.get("date_inserted"))
        if heads is None or max(heads) >= SEAL:
            raise Stop(f"a row of the genome table, second copy, is not inserted before {SEAL}")
        parts = bioproject_parts(row.get("bioproject_accession"))
        ids.add(genome)
        accessions |= {p for p in parts if is_accession(p)}
        if genome in wanted and parts:
            found.setdefault(genome, []).append(frozenset(parts))
    project, why = {}, {}
    for g in wanted:
        rows = found.get(g)
        one = rows is not None and len(set(rows)) == 1 and len(rows[0]) == 1
        part = next(iter(rows[0])) if one else None
        project[g] = part if part is not None and is_accession(part) else None
        why[g] = "not_in_table" if g not in ids else "no_accession" if part is not None and project[g] is None \
            else None
    return project, why, ids, accessions, n


def read_hc20(path, wanted):
    """({genome: HC20 group or None}, ids read, rows, genomes the table lacks)."""
    found, ids, n = {}, set(), 0
    for row in jsonl(path, "the genome table, 5 Nov"):
        n += 1
        genome = value_text(row.get("genome_id"))
        if genome is None:
            raise Stop("a row of the genome table, 5 Nov, has no genome_id")
        ids.add(genome)
        if genome in wanted:
            found.setdefault(genome, []).append(frozenset(values(row.get("cgmlst_hc20"))))
    out = {}
    for g in wanted:
        rows = found.get(g)
        out[g] = next(iter(rows[0])) if rows and len(set(rows)) == 1 and len(rows[0]) == 1 else None
    return out, ids, n, sum(1 for g in wanted if g not in found)


# 4. The rows, the readings and the groups.

def base_rows(entering, pairs, project, why):
    """(rows that pass reasons 1 to 5, {dataset: Counter of reasons}). A row leaves under the first that applies."""
    rows, out = [], {}
    for rep in entering:
        c = collections.Counter()
        copies = {i for g in rep["kset_duplicates"] or [] for i in g}
        fold = rep["lineage_split"]["test_fold"]
        for i, s in enumerate(rep["samples"]):
            pair = pairs.get((s["id"], rep["antibiotic"]))
            if i in copies:
                c["copy_group"] += 1
            elif project[s["id"]] is None:
                c["bioproject"] += 1
                if why[s["id"]]:
                    c["bioproject_" + why[s["id"]]] += 1
            elif pair is None:
                c["no_lab_record"] += 1
            elif pair["blank"]:
                c["blank_method_or_standard"] += 1
            elif len(pair["tuples"]) != 1:
                c["assay_not_single"] += 1
            else:
                method, standard, year = next(iter(pair["tuples"]))
                rows.append({"dataset": rep["id"], "species": rep["species_slug"], "genome": s["id"],
                             "fold": fold[i], "lineage": s["lineage"], "bioproject": project[s["id"]],
                             "method": method, "standard": standard, "year": year, "resistant": s["resistant"],
                             "mlst": (s["mlst"] or "").strip() or None, "intermediate": pair["intermediate"]})
        out[rep["id"]] = c
    return rows, out


def key(r, extra=()):
    return (r["dataset"], r["fold"], r["lineage"], r["bioproject"], r["method"], r["standard"], r["year"]) + \
        tuple(r[k] for k in extra)


def groups_of(rows, extra=()):
    """{group key: rows} of the groups that hold both lab outcomes, and the rows of the other groups."""
    by = collections.defaultdict(list)
    for r in rows:
        by[key(r, extra)].append(r)
    both = {k: g for k, g in by.items() if len({r["resistant"] for r in g}) == 2}
    rest = [r for k, g in by.items() if k not in both for r in g]
    return both, rest


def pairs_in(group):
    r = sum(1 for x in group if x["resistant"])
    return r * (len(group) - r)


def reading(base, name):
    """One reading: (eligible groups, {dataset: Counter of the reading's own reasons}, rows that enter grouping)."""
    out = collections.defaultdict(collections.Counter)
    rows = base
    if name == "B":
        rows = []
        for r in base:
            if r["intermediate"]:
                out[r["dataset"]]["intermediate"] += 1
            else:
                rows.append(r)
    both, rest = groups_of(rows)
    for r in rest:
        out[r["dataset"]]["one_lab_outcome"] += 1
    return both, out, rows


def support(groups):
    """{level: {name: counts}} of rows, genomes, groups, pairs, datasets and species over eligible groups."""
    acc = {"dataset": collections.defaultdict(lambda: [0, set(), 0, 0, set(), set()]),
           "species": collections.defaultdict(lambda: [0, set(), 0, 0, set(), set()]),
           "overall": collections.defaultdict(lambda: [0, set(), 0, 0, set(), set()])}
    for g in groups.values():
        d, sp = g[0]["dataset"], g[0]["species"]
        for level, name in (("dataset", d), ("species", sp), ("overall", "all")):
            a = acc[level][name]
            a[0] += len(g)
            a[1] |= {r["genome"] for r in g}
            a[2] += 1
            a[3] += pairs_in(g)
            a[4].add(d)
            a[5].add(sp)
    return {level: {n: {"rows": a[0], "genomes": len(a[1]), "groups": a[2], "pairs": a[3], "datasets": len(a[4]),
                        "species": len(a[5])} for n, a in names.items()} for level, names in acc.items()}


def weights(groups):
    """{group key: weight}: groups equal inside a dataset, datasets equal inside a species, species equal overall."""
    by_dataset = collections.Counter(k[0] for k in groups)
    species_of = {g[0]["dataset"]: g[0]["species"] for g in groups.values()}
    datasets_in = collections.Counter(species_of.values())
    n_species = len(datasets_in)
    return {k: fractions.Fraction(1, by_dataset[k[0]] * datasets_in[species_of[k[0]]] * n_species) for k in groups}


def shares(groups):
    """The largest BioProject's share of eligible rows and of analysis weight, as integers alone."""
    w = weights(groups)
    rows, weight = collections.Counter(), collections.defaultdict(fractions.Fraction)
    for k, g in groups.items():
        rows[k[3]] += len(g)
        weight[k[3]] += w[k]
    if sum(w.values()) != 1:
        raise Stop("the analysis weights do not add up to 1")
    top = max(rows.values())
    its = max(weight[b] for b in rows if rows[b] == top)
    most = max(weight.values())
    return {"bioprojects": len(rows), "eligible_rows": sum(rows.values()), "largest_bioproject_rows": top,
            "its_weight_numerator": its.numerator, "its_weight_denominator": its.denominator,
            "largest_weight_numerator": most.numerator, "largest_weight_denominator": most.denominator}


def typing(rows, groups, hc20, names):
    """{step: {level: {name: counts}}}: each narrower match beside the wider matches on its surviving rows.

    rows: the reading's rows that enter grouping. groups: the reading's eligible groups. names: {level: names}
    of the datasets and species that enter, and "all". Rows lacking a value are counted twice: among the rows
    that enter the groups, and among the eligible rows of the reading."""
    eligible_rows = {level: {n: c["rows"] for n, c in ns.items()} for level, ns in support(groups).items()}
    eligible = {(r["dataset"], r["genome"]) for g in groups.values() for r in g}
    steps = [("mlst", ("mlst",))] + ([("mlst_hc20", ("mlst", "hc20"))] if hc20 is not None else [])
    out = {}
    for step, extra in steps:
        typed = []
        lacking = collections.defaultdict(collections.Counter)
        lacking_eligible = collections.defaultdict(collections.Counter)
        for r in rows:
            r2 = dict(r, hc20=hc20.get(r["genome"]) if hc20 is not None else None)
            if any(r2[k] is None for k in extra):
                for level, name in (("dataset", r["dataset"]), ("species", r["species"]), ("overall", "all")):
                    lacking[level][name] += 1
                    if (r["dataset"], r["genome"]) in eligible:
                        lacking_eligible[level][name] += 1
            else:
                typed.append(r2)
        narrow, _ = groups_of(typed, extra)
        surviving = [r for g in narrow.values() for r in g]
        counts = support(narrow)
        primary = support(groups_of(surviving)[0])
        mlst = support(groups_of(surviving, ("mlst",))[0]) if step == "mlst_hc20" else None
        table = {}
        for level, ns in names.items():
            table[level] = {}
            for n in ns:
                c = counts.get(level, {}).get(n, {"rows": 0, "genomes": 0, "groups": 0, "pairs": 0, "datasets": 0,
                                                   "species": 0})
                p = primary.get(level, {}).get(n, {"groups": 0, "pairs": 0})
                m = mlst.get(level, {}).get(n, {"groups": 0, "pairs": 0}) if mlst is not None else None
                table[level][n] = dict(c, lacking_value_of_rows_entering_groups=lacking[level][n],
                                       lacking_value_of_eligible_rows=lacking_eligible[level][n],
                                       primary_groups=p["groups"],
                                       primary_pairs=p["pairs"], mlst_groups="" if m is None else m["groups"],
                                       mlst_pairs="" if m is None else m["pairs"],
                                       rows_lost_from_primary=eligible_rows[level].get(n, 0) - c["rows"])
        out[step] = table
    return out


# 5. The row lists, the export and the scan.

def row_list(groups):
    """The frozen row list's bytes: one tab-separated line an eligible row, sorted, UTF-8, LF."""
    lines = []
    for g in groups.values():
        for r in g:
            fields = [r["dataset"], r["genome"], str(r["fold"]), str(r["lineage"]), r["bioproject"], r["method"],
                      r["standard"], r["year"] or ""]
            if any(re.search(r"[\t\r\n]", f) for f in fields):
                raise Stop(f"{r['dataset']}: a field of an eligible row holds a tab or a line break")
            lines.append("\t".join(fields) + "\n")
    return "".join(sorted(lines)).encode("utf-8")


def tsv(columns, rows):
    out = ["\t".join(columns)]
    for r in rows:
        out.append("\t".join(str(r.get(c, "")) for c in columns))
    return ("\n".join(out) + "\n").encode("utf-8")


def count_rows(entering, left_out, base_reasons, readings):
    """The rows of counts.tsv: per dataset, per species and overall, for each reading."""
    reports = {rep["id"]: rep for rep in entering + left_out}
    zero = {c: 0 for c in COUNT_COLUMNS[3:]}
    lines = []
    for name in ("A", "B"):
        groups, own = readings[name]["groups"], readings[name]["reasons"]
        sup = support(groups)
        per = {"dataset": {}, "species": {}, "overall": {}}
        for d, rep in reports.items():
            sp = rep["species_slug"]
            c = dict(zero)
            c["reports"] = 1
            c["genome_rows"] = len(rep["samples"])
            if d in base_reasons:
                c["from_mic_rows"] = sum(1 for s in rep["samples"] if s["from_mic"])
                b, o = base_reasons[d], own.get(d, collections.Counter())
                for r in REASONS:
                    c["out_" + r] = b[r]
                c["out_bioproject_not_in_table"] = b["bioproject_not_in_table"]
                c["out_bioproject_no_accession"] = b["bioproject_no_accession"]
                c["out_intermediate"] = o["intermediate"]
                c["out_one_lab_outcome"] = o["one_lab_outcome"]
            else:
                c["reports_left_out"] = 1
            per["dataset"][d] = (c, sp)
        for level in ("species", "overall"):
            for d, (c, sp) in per["dataset"].items():
                n = sp if level == "species" else "all"
                t = per[level].setdefault(n, (dict(zero), None))[0]
                for k in ("reports", "reports_left_out", "genome_rows", "from_mic_rows") + OUT_COLUMNS:
                    t[k] += c[k]
        for level in ("dataset", "species", "overall"):
            for n in sorted(per[level]):
                c = per[level][n][0]
                c.update(sup.get(level, {}).get(n, {"rows": 0, "genomes": 0, "groups": 0, "pairs": 0, "datasets": 0,
                                                    "species": 0}))
                lines.append(dict(c, level=level, name=n, reading=name))
    order = {"dataset": 0, "species": 1, "overall": 2}
    return sorted(lines, key=lambda r: (order[r["level"]], r["name"], r["reading"]))


def typing_rows(readings):
    lines = []
    order = {"dataset": 0, "species": 1, "overall": 2}
    for name in ("A", "B"):
        for step, table in readings[name]["typing"].items():
            for level, names in table.items():
                for n, c in names.items():
                    lines.append(dict(c, level=level, name=n, reading=name, step=step))
    return sorted(lines, key=lambda r: (order[r["level"]], r["name"], r["reading"], r["step"]))


def export_files(facts, readings, lists):
    """{file name: bytes} of the export: counts alone, with no genome id, accession or decimal point."""
    inputs = [
        {"input": "reports", "sha256": facts["sums"]["reports"], "count": facts["reports"],
         "note": "each report equals its line of the release's sums list; the sha256 is that list's"},
        {"input": "lab records", "sha256": facts["sums"]["lab records"], "count": facts["record_lines"],
         "note": f"lines, each inserted before {SEAL}; {DROPPED_RECORDED} records dropped at the cut is a recorded "
                 f"figure of ruling R49, not measured here"},
        {"input": "drug names", "sha256": facts["sums"]["drug names"], "count": facts["drug_names"],
         "note": "names"},
        {"input": "label table", "sha256": facts["sums"]["label table"], "count": facts["label_pairs"],
         "note": "pairs of the release's rows read from it"},
        {"input": "genome table, second copy", "sha256": facts["sums"]["genome table, second copy"],
         "count": facts["genome_rows"],
         "note": f"rows, each inserted before {SEAL}; a row with a blank BioProject is skipped; distinct genome "
                 f"ids it lacks: {facts['not_in_copy_every']} of the {facts['genomes_every']} of every report, "
                 f"{facts['not_in_copy']} of the {facts['genomes_entering']} of the reports that enter; genomes of "
                 f"the reports that enter whose one BioProject part has no accession's shape: "
                 f"{facts['no_accession']}"},
    ]
    if facts["hc20_rows"] is None:
        inputs.append({"input": "genome table, 5 Nov", "sha256": "-", "count": 0,
                       "note": "not given: the HC20 read is dropped and typing.tsv holds no mlst_hc20 step"})
    else:
        inputs.append({"input": "genome table, 5 Nov", "sha256": facts["sums"]["genome table, 5 Nov"],
                       "count": facts["hc20_rows"],
                       "note": f"rows; the table has no insertion time, so no row of it is checked against "
                               f"{SEAL}; it gives the mlst_hc20 step of typing.tsv; genomes of the reports that "
                               f"enter not in it: {facts['not_in_hc20']}"})
    for name in ("A", "B"):
        inputs.append({"input": f"row list, reading {name}", "sha256": hashlib.sha256(lists[name]).hexdigest(),
                       "count": lists[name].count(b"\n"), "note": "eligible rows; the list stays local"})
    checks = [
        {"check": "lab records", "result": f"{facts['record_lines']} lines, each inserted before {SEAL}"},
        {"check": "reports", "result": f"{facts['entering']} of {facts['reports']} reports hold the five learners "
                                       f"under the lineage split; left out, holding none of them: {facts['left_out']}"},
        {"check": "folds", "result": f"each of the {facts['entering']} has {FOLDS} folds, and test_fold and every "
                                     f"scores list as long as its genome rows"},
        {"check": "lineage clusters", "result": "each in one fold"},
        {"check": "copy groups", "result": "each in one fold"},
        {"check": "lab calls", "result": f"rows with from_mic false, each equal to the label table: "
                                         f"{facts['labels_equal']}; rows with from_mic true, counted apart: "
                                         f"{facts['from_mic']}"},
        {"check": "scan", "result": f"{facts['scan_ids']} genome ids and {facts['scan_accessions']} "
                                    f"accession-shaped BioProject parts read by the run, looked up in every export "
                                    f"file after a positive control, before any file was written"},
        {"check": "scores opened", "result": "no"},
    ]
    out = {
        "counts.tsv": tsv(COUNT_COLUMNS, readings["counts"]),
        "typing.tsv": tsv(TYPING_COLUMNS, readings["typing_rows"]),
        "shares.tsv": tsv(SHARE_COLUMNS, [dict(readings[n]["shares"], reading=n) for n in ("A", "B")]),
        "checks.tsv": tsv(("check", "result"), checks),
        "inputs.tsv": tsv(("input", "sha256", "count", "note"), inputs) + b"scores opened: no\n",
    }
    return out


def hits(data, ids, accessions):
    """The count of genome ids, accessions, id-shaped and accession-shaped tokens in one file's bytes."""
    text = data.decode("utf-8")
    n = sum(text.count(x) for x in ids if x in text) + sum(text.count(x) for x in accessions if x in text)
    return n + len(ID_SHAPED.findall(text)) + len(ACCESSION_SHAPED.findall(text))


def scan(files, ids, accessions):
    """[(file name, hits)] of the export files that hold any. A positive control runs first."""
    first = sorted(files)[0]
    planted = files[first] + min(ids).encode("utf-8")
    if hits(planted, ids, accessions) == 0:
        raise Stop("the scan's positive control found no genome id in a file that holds one")
    if accessions and hits(files[first] + min(accessions).encode("utf-8"), ids, accessions) == 0:
        raise Stop("the scan's positive control found no accession in a file that holds one")
    return [(name, n) for name in sorted(files) if (n := hits(files[name], ids, accessions))]


def write(out, lists, files):
    """Write the row lists and the export. On an OSError, remove every directory and file it made, then stop."""
    made = []
    try:
        missing = []
        p = out
        while not p.exists():
            missing.append(p)
            p = p.parent
        for d in reversed(missing):
            d.mkdir()
            made.append(d)
        local = out / "local"
        local.mkdir(mode=0o700)
        made.append(local)
        for name in ("A", "B"):
            path = local / ROW_LIST[name]
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            made.append(path)
            with os.fdopen(fd, "wb") as f:
                f.write(lists[name])
            os.chmod(path, 0o600)
        export = out / "export"
        export.mkdir()
        made.append(export)
        for name, data in files.items():
            fd = os.open(export / name, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666)
            made.append(export / name)
            with os.fdopen(fd, "wb") as f:
                f.write(data)
    except OSError as e:
        left = 0
        for m in reversed(made):
            try:
                m.rmdir() if m.is_dir() else m.unlink()
            except OSError:
                left += 1
        raise Stop(f"the output could not be written ({type(e).__name__}); of the {len(made)} paths it made, "
                   f"{left} could not be removed") from None


def build(args):
    """(row lists, export files, ids and accessions read, printout). Raises Stop on any check."""
    out = Path(args.out)
    if out.exists() and not out.is_dir():
        raise Stop("the output path is not a directory")
    if (out / "local").exists() or (out / "export").exists():
        raise Stop("the output directory already holds local or export")
    if args.hc20 != ["-"] and (len(args.hc20) != 2 or not SHA.fullmatch(args.hc20[1])):
        raise Stop("--hc20 takes a file and its sha256, or -")
    for what, given in (("drug names", args.drug_names[1]), ("label table", args.labels[1])):
        if not SHA.fullmatch(given):
            raise Stop(f"the {what} sha256 is not 64 lowercase hex digits")
    files = report_files(args.reports)
    sums = check_sums(args)
    sums["reports"] = REPORT_SUMS_SHA256

    entering, left_out = load_reports(files)
    if not entering:
        raise Stop("no report holds the five learners")
    wanted_genomes = {s["id"] for rep in entering for s in rep["samples"]}
    wanted_pairs = {(s["id"], rep["antibiotic"]) for rep in entering for s in rep["samples"]}
    resolve = drug_names(args.drug_names[0])
    pairs, record_ids, record_lines = read_records(args.records, resolve, wanted_pairs)
    table, label_ids = read_labels(args.labels[0], {(s["id"], rep["antibiotic"]) for rep in entering
                                                    for s in rep["samples"] if not s["from_mic"]})
    labels_equal, from_mic = check_labels(entering, table)
    project, why, genome_ids, accessions, genome_rows = read_genomes(args.genomes, wanted_genomes)
    every_genome = {s["id"] for rep in entering + left_out for s in rep["samples"]}
    hc20 = hc20_ids = hc20_rows = not_in_hc20 = None
    if args.hc20 != ["-"]:
        hc20, hc20_ids, hc20_rows, not_in_hc20 = read_hc20(args.hc20[0], wanted_genomes)

    base, base_reasons = base_rows(entering, pairs, project, why)
    names = {"dataset": {rep["id"] for rep in entering}, "species": {rep["species_slug"] for rep in entering},
             "overall": {"all"}}
    readings = {}
    for name in ("A", "B"):
        groups, own, rows = reading(base, name)
        if not groups:
            raise Stop(f"no group qualifies in reading {name}")
        readings[name] = {"groups": groups, "reasons": own, "shares": shares(groups),
                          "typing": typing(rows, groups, hc20, names)}
    readings["counts"] = count_rows(entering, left_out, base_reasons, readings)
    readings["typing_rows"] = typing_rows(readings)
    lists = {name: row_list(readings[name]["groups"]) for name in ("A", "B")}

    ids = every_genome | record_ids | label_ids | genome_ids
    if hc20_ids:
        ids |= hc20_ids
    ids.discard(None)
    facts = {"sums": sums, "reports": len(files), "entering": len(entering), "left_out": len(left_out),
             "record_lines": record_lines, "drug_names": len(resolve), "label_pairs": len(table),
             "genome_rows": genome_rows, "genomes_every": len(every_genome), "genomes_entering": len(wanted_genomes),
             "not_in_copy_every": len(every_genome - genome_ids), "not_in_copy": len(wanted_genomes - genome_ids),
             "no_accession": sum(1 for w in why.values() if w == "no_accession"), "hc20_rows": hc20_rows,
             "not_in_hc20": not_in_hc20, "labels_equal": labels_equal, "from_mic": from_mic,
             "scan_ids": len(ids), "scan_accessions": len(accessions)}
    export = export_files(facts, readings, lists)
    said = [f"stage 1: {len(entering)} of {len(files)} reports enter, {len(left_out)} left out; "
            f"{record_lines} lab records, each inserted before {SEAL}"]
    for name in ("A", "B"):
        c = next(r for r in readings["counts"] if r["level"] == "overall" and r["reading"] == name)
        said.append(f"reading {name}: {c['rows']} eligible rows, {c['genomes']} genomes, {c['groups']} groups, "
                    f"{c['pairs']} pairs, {c['datasets']} datasets, {c['species']} species; row list sha256 "
                    f"{hashlib.sha256(lists[name]).hexdigest()}")
    said.append("scores opened: no")
    return lists, export, ids, accessions, said


def parse(argv):
    p = argparse.ArgumentParser(prog="stage1.py", description="Stage 1 of the matched accuracy protocol.")
    p.add_argument("--reports", required=True)
    p.add_argument("--records", required=True)
    p.add_argument("--drug-names", required=True, nargs=2, metavar=("FILE", "SHA256"))
    p.add_argument("--labels", required=True, nargs=2, metavar=("FILE", "SHA256"))
    p.add_argument("--genomes", required=True)
    p.add_argument("--hc20", required=True, nargs="+", metavar="FILE SHA256 | -")
    p.add_argument("--out", required=True)
    return p.parse_args(argv)


def main(argv):
    args = parse(argv)
    try:
        lists, export, ids, accessions, said = build(args)
        found = scan(export, ids, accessions)
        if not found:
            write(Path(args.out), lists, export)
    except Stop as e:
        print(f"stage1: stopped: {e}")
        return 1
    if found:
        for name, n in found:
            print(f"{name} {n}")
        return 1
    print("\n".join(said))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
