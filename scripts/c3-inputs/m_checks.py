"""The three checks of plan block C3's run on its datasets of M, before any dataset of M is read.

  python3 -I m_checks.py species <lists.tsv>
  python3 -I m_checks.py build <page-binaries.tsv> <species slug> <commit> <code sha256> [<cargo home>]
  python3 -I m_checks.py results <reference dir> <species slug> <datasets.tsv>
  python3 -I m_checks.py counts <expect.tsv> <picks dir> <scores dir> <reports dir> <study dir>

species prints the check species by the entry's rule: among the species
that hold a dataset of M and a first-run dataset on lists.tsv, the one with
the fewest first-run datasets, the first by name on a tie.

build (Same build) holds the code sha256 of the first-run commit's own
build to the species' row of reference/input-arms/page-binaries.tsv: the
row's commit must be <commit> and its code sha256 <code sha256>. The code
sha256 moves with cargo's home, whose registry paths are built into the
binary, so a difference names the home the build used (<cargo home>, as
binary-first.txt records it).

results (Same results) holds score.py's datasets.tsv over the check
species' first-run datasets, fitted by the changed binary, to that
species' rows of reference/input-arms/page-datasets.tsv, itself held to the
directory's SHA256SUMS: page_tables.py's own page_datasets builds the rows,
and every cell must be equal at 12 decimals, with no row missing or extra.

counts (Same counts) holds, for every result expect.tsv lists, each fold's
calls taken again from the genomes' scores to the counts the fit wrote in
its pick file: n, tp, fp, tn and fn equal, the log loss sum within 1e-9.
The scores are the lineage fits' score files (grm learn --exp-scores,
<scores dir>/<arm>/<dataset>.lineage-s<seed>.scores.json), the default's
and the twin's kept reports (<reports dir>/<arm>/<dataset>.lineage-s<seed>.json),
which must also hold the score file's scores, and the study fits'
results (<study dir>/<arm>/<dataset>.study-whole.json). A result refused
for no candidate has no scores and is passed over; its default is checked.

Each form prints what it checked, counts and sha256 only, and names what
differs, with no genome id. It exits 1 on any difference.
"""

import csv
import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIRST_RUN_COMMIT = "8bfd4ca8a64678c69431a6bc22b869af6d5cb81b"
SCORES_FORMAT = "grm-c3-scores-1"
FOLD_FORMAT = "grm-fold-file-result-v1"
LEARNERS = ["scm", "cart", "rf", "gbt", "lasso"]
CALL_CUTOFF = 0.5
# Log loss clips each score to [LOSS_CLIP, 1 - LOSS_CLIP], as grm does (src/cutoff.rs).
LOSS_CLIP = 1e-4
LOSS_TOLERANCE = 1e-9
DEFAULT, TWIN = "k31", "twin"


class Stop(Exception):
    pass


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def load(path):
    with open(path) as f:
        return json.load(f)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# ---------- the check species ----------

def check_species(lists_rows):
    """The species with the fewest first-run datasets among those with a dataset of M and a first-run dataset."""
    m = {r["species_slug"] for r in lists_rows if r["in_M"] == "1"}
    first = {}
    for r in lists_rows:
        if r["first_run"] == "1":
            first[r["species_slug"]] = first.get(r["species_slug"], 0) + 1
    both = sorted((n, s) for s, n in first.items() if s in m)
    if not both:
        raise Stop("no species holds a dataset of M and a first-run dataset")
    return both[0][1], both[0][0]


# ---------- Same build ----------

def same_build(binaries_rows, slug, commit, code, home=None):
    rows = [r for r in binaries_rows if r["species_slug"] == slug]
    if len(rows) != 1:
        raise Stop(f"page-binaries.tsv holds {len(rows)} rows for {slug}")
    row = rows[0]
    if row["commit"] != commit:
        raise Stop(f"Same build: the record's commit for {slug} is {row['commit']}, the build's {commit}")
    built = f"built {code}" + (f" with cargo's home {home}" if home else " with cargo's home not named")
    if row["code"] != code:
        raise Stop(f"Same build: the code sha256 differs for {slug}: recorded {row['code']}, {built}")
    return f"Same build: {slug} at {commit}, code sha256 {code}, equal to page-binaries.tsv" \
        + (f", cargo's home {home}" if home else "")


# ---------- Same results ----------

def held_sums(reference):
    """The reference directory's SHA256SUMS, checked for page-datasets.tsv."""
    reference = Path(reference)
    sums = {}
    for line in (reference / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        sums[name.lstrip("*")] = digest
    for name in ("page-datasets.tsv",):
        if sums.get(name) != sha256(reference / name):
            raise Stop(f"{name} differs from SHA256SUMS")
    return sums


class _Data:
    """The one field of page_tables.py's reading that page_datasets reads."""

    def __init__(self, slug, rows):
        self.species = [slug]
        self.score = {slug: {"datasets.tsv": rows}}


def same_results(reference, slug, dataset_rows):
    page = _load("c3_page_tables_for_m", "page_tables.py")
    held_sums(reference)
    recorded = [r for r in read_tsv(Path(reference) / "page-datasets.tsv") if r["dataset"].split("__")[0] == slug]
    if not recorded:
        raise Stop(f"page-datasets.tsv holds no row of {slug}")
    columns = page.DATASET_TABLE_COLUMNS
    rebuilt = [{c: page.C9.fmt(r.get(c)) for c in columns} for r in page.page_datasets(_Data(slug, dataset_rows))]
    key = lambda r: (r["dataset"], r["split"], r["pair"], r["learner"])
    have, want = {key(r): r for r in rebuilt}, {key(r): r for r in recorded}
    if len(have) != len(rebuilt):
        raise Stop("Same results: a dataset, split, pair and learner is given twice")
    missing, extra = sorted(set(want) - set(have)), sorted(set(have) - set(want))
    if missing or extra:
        raise Stop(f"Same results: {len(missing)} recorded rows not given back, {len(extra)} rows not in the record; "
                   f"first {(missing or extra)[0]}")
    for k in sorted(want):
        for c in columns:
            a, b = have[k][c], want[k][c]
            if a == b:
                continue
            try:
                equal = page.r12(float(a)) == page.r12(float(b))
            except ValueError:
                equal = False
            if not equal:
                raise Stop(f"Same results: {k} {c}: the record gives {b}, the changed binary {a}")
    return f"Same results: {len(want)} of {len(want)} rows of {slug} in page-datasets.tsv given back at 12 decimals"


# ---------- Same counts ----------

def calls(truth, scores):
    """src/cutoff.rs, calls: counts at CALL_CUTOFF and the clipped log loss sum, in the order given."""
    c = {"n": len(truth), "tp": 0, "fp": 0, "tn": 0, "fn": 0, "loss_sum": 0.0}
    for t, s in zip(truth, scores):
        p = min(max(s, LOSS_CLIP), 1 - LOSS_CLIP)
        called = s >= CALL_CUTOFF
        c[("tp" if called else "fn") if t else ("fp" if called else "tn")] += 1
        c["loss_sum"] -= math.log(p) if t else math.log(1 - p)
    return c


def from_scores(doc):
    """[(genome id, resistant, fold, [five scores])] from a score file."""
    if doc.get("format") != SCORES_FORMAT or doc.get("learners") != LEARNERS:
        raise Stop("a score file is not grm-c3-scores-1 with the five learners in order")
    return [(g["genome_id"], g["resistant"], g["fold"], list(g["scores"])) for g in doc["genomes"]]


def _five(learners):
    by = {r["learner"]: r for r in learners}
    if not all(l in by for l in LEARNERS):
        raise Stop("a split does not hold the five learners")
    return [by[l]["scores"] for l in LEARNERS]


def from_report(doc):
    """[(genome id, resistant, fold, [five scores])] from a report's lineage split, in sample order."""
    split = next((s for s in doc["splits"] if s["scheme"] == "lineage"), None)
    if split is None:
        raise Stop("a report has no lineage split")
    if split.get("unavailable_reason") is not None:
        return []
    scores = _five(split["learners"])
    return [(s["id"], s["resistant"], split["test_fold"][i], [x[i] for x in scores])
            for i, s in enumerate(doc["samples"])]


def from_study(doc):
    """[(genome id, resistant, fold, [five scores])] from a study fit's fold result, by tested position."""
    if doc.get("format") != FOLD_FORMAT:
        raise Stop("a study result is not a fold result")
    scores = _five(doc["split"]["learners"])
    return [(s["genome_id"], s["resistant"], s["fold"], [x[i] for x in scores]) for i, s in enumerate(doc["samples"])]


def fold_calls(genomes, folds, left_out=frozenset()):
    """Per fold, per learner: the calls on that fold's genomes not in left_out, in their order; and per fold
    the test genomes and resistant ones."""
    out = []
    for f in range(folds):
        members = [g for g in genomes if g[2] == f and g[0] not in left_out]
        truth = [g[1] for g in members]
        out.append(([calls(truth, [g[3][j] for g in members]) for j in range(len(LEARNERS))],
                    len(members), sum(truth)))
    return out


def compare(doc, genomes, where):
    """The differences between a pick file's per-fold counts and those taken again from the genomes' scores."""
    if doc.get("unavailable_reason") is not None:
        return [] if not genomes else [f"{where}: the pick file's split is unavailable, the scores hold genomes"]
    folds = doc["folds"]
    if any(not (0 <= g[2] < len(folds)) for g in genomes):
        return [f"{where}: a genome's fold is not one of the pick file's {len(folds)}"]
    diffs = []
    if (len(genomes), sum(g[1] for g in genomes)) != (doc["tested"], doc["tested_resistant"]):
        diffs.append(f"{where}: {len(genomes)} genomes scored, the pick file tests {doc['tested']}")
    for f, (got, test, test_r) in enumerate(fold_calls(genomes, len(folds))):
        rec = folds[f]
        if (test, test_r) != (rec["test"], rec["test_resistant"]):
            diffs.append(f"{where}: fold {f} tests {test} genomes, {test_r} R; the pick file {rec['test']}, "
                         f"{rec['test_resistant']} R")
        for j, (c, w) in enumerate(zip(got, rec["outer"])):
            if w.get("learner") != LEARNERS[j]:
                diffs.append(f"{where}: fold {f} names learner {w.get('learner')} at {LEARNERS[j]}'s place")
                continue
            for k in ("n", "tp", "fp", "tn", "fn"):
                if c[k] != w[k]:
                    diffs.append(f"{where}: fold {f}, {LEARNERS[j]}, {k}: {c[k]} from the scores, {w[k]} written")
            if abs(c["loss_sum"] - w["loss_sum"]) > LOSS_TOLERANCE * max(1.0, abs(w["loss_sum"])):
                diffs.append(f"{where}: fold {f}, {LEARNERS[j]}, loss_sum: {c['loss_sum']!r} from the scores, "
                             f"{w['loss_sum']!r} written")
    return diffs


def score_path(scores, e):
    return Path(scores) / e["arm"] / f"{e['dataset']}.lineage-s{e['seed']}.scores.json"


def report_path(reports, e):
    return Path(reports) / e["arm"] / f"{e['dataset']}.lineage-s{e['seed']}.json"


def study_path(study, e):
    return Path(study) / e["arm"] / f"{e['dataset']}.{e['scheme']}.json"


def genomes_of(e, scores, study):
    """The genomes' scores behind one expect row's result, from its score file or study result."""
    if e["split"] == "lineage":
        doc = load(score_path(scores, e))
        want = (e["dataset"], "lineage", int(e["seed"]), DEFAULT if e["arm"] == TWIN else e["arm"])
        if (doc["dataset"], doc["scheme"], doc["seed"], doc["arm"]) != want:
            raise Stop(f"{score_path(scores, e).name}: not the score file of {e['arm']}'s fit")
        return from_scores(doc)
    doc = load(study_path(study, e))
    if (doc["id"], doc["scheme"]) != (e["dataset"], e["scheme"]):
        raise Stop(f"{study_path(study, e).name}: not the study result of {e['dataset']}")
    return from_study(doc)


def same_counts(expect, picks, scores, reports, study):
    """(fits checked, reports checked, refused passed over). Raises Stop naming the first differences."""
    picks = Path(picks)
    fits = kept = refused = 0
    diffs = []
    for e in expect:
        path = picks / e["file"]
        if not path.is_file():
            if (picks / (e["file"][: -len(".pick.json")] + ".refused.json")).is_file() and e["arm"] != DEFAULT:
                refused += 1
                continue
            raise Stop(f"{e['file']}: no pick file")
        doc = load(path)
        try:
            genomes = genomes_of(e, scores, study)
        except FileNotFoundError:
            raise Stop(f"{e['file']}: its genomes' scores were not kept")
        diffs += compare(doc, genomes, e["file"])
        fits += 1
        if e["split"] == "lineage" and e["arm"] in (DEFAULT, TWIN):
            rpath = report_path(reports, e)
            if not rpath.is_file():
                raise Stop(f"{e['file']}: the fit's report was not kept")
            from_rep = from_report(load(rpath))
            if from_rep != genomes:
                diffs.append(f"{rpath.name}: the kept report's genomes, folds or scores differ from the score file's")
            diffs += compare(doc, from_rep, rpath.name)
            kept += 1
    if diffs:
        raise Stop(f"Same counts: {len(diffs)} differences; " + "; ".join(diffs[:5]))
    return fits, kept, refused


def main(argv):
    if not argv or argv[0] not in ("species", "build", "results", "counts"):
        print("\n".join(__doc__.strip().splitlines()[2:6]), file=sys.stderr)
        return 2
    form, args = argv[0], argv[1:]
    try:
        if form == "species" and len(args) == 1:
            slug, n = check_species(read_tsv(args[0]))
            print(slug)
            print(f"check species: {slug}, {n} first-run datasets", file=sys.stderr)
        elif form == "build" and len(args) in (4, 5):
            print(same_build(read_tsv(args[0]), *args[1:]))
        elif form == "results" and len(args) == 3:
            print(same_results(args[0], args[1], read_tsv(args[2])))
        elif form == "counts" and len(args) == 5:
            fits, kept, refused = same_counts(read_tsv(args[0]), *args[1:])
            print(f"Same counts: {fits} fits and {kept} kept reports give back every fold's counts; "
                  f"{refused} results refused for no candidate hold no scores")
        else:
            print("\n".join(__doc__.strip().splitlines()[2:6]), file=sys.stderr)
            return 2
    except Stop as e:
        print(f"m_checks {form}: stopped: {e}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
