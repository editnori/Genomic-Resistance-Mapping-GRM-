"""The tables of docs/cap-arms-result.md: plan block C2's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <count dir> <runs dir> <lists.tsv> <checks.txt> <out dir>
  python3 -I page_tables.py print <out dir>

<count dir> holds what `score.py --count <count dir> --withdrawn cap-4000
<each species' score/datasets.tsv>` wrote: count.tsv and conditions.tsv,
and any other table it writes beside them. <runs dir> holds one folder per
species, named by its slug, in run_species.sh's OUT layout: score/
(score.py's checks.tsv, draws.tsv and datasets.tsv; score/score.log is
never read), results/ (cohorts.py's added-labelled.tsv and cohorts.json;
and, where the run wrote them, fetch.json, folds-lineage-<dataset>.tsv,
folds-study.tsv and study-files-cap-all.tsv), expect.tsv, binary.txt,
logs/timing.tsv and logs/folds-<dataset>.log. <lists.tsv> is
scripts/c2-cap/lists.tsv: the first run is its rows with first_run 1 and
in_M 0. <checks.txt> is run_checks.py's printout over the run's copies,
made before any score was read: gate 1 holds that it says each of
run_checks.py's checks passed, and the page is not written otherwise.

`write` passes every gate of GATES before it writes a file, and stops at
the first that fails, naming it. It then writes <out dir> (absent or
empty): the scorer's count files under count/, copied as they are,
checks.txt, copied, the page's tables as page-*.tsv, gates.txt, README.md
and SHA256SUMS. Where
score.py's own tables are the page's files (count.tsv, conditions.tsv),
the script checks them against score.py's own count of the datasets.tsv
files and copies them; every other table is summed here from score.py's
tables, cohorts.py's counts and the run's counts-only files. Every row is
an aggregate over a dataset, a species, a split, an arm or a learner: no
file holds a genome, a BioProject, a sequence or a feature.

The count that decides is count.tsv, the count of every dataset fitted.
Beside it, deciding nothing: cap-all's count without the datasets that add
fewer than 20 labelled genomes under cap-all (the lead's ruling of 12:19Z
on 23 Nov 2026), made with score.py's own counting code and Holm over its
own rows; the study split's rows; the stricter reading of the conditions;
condition 5, which reads the study split (docs/cap-arms.md, Limits: no
condition reads it); and every table per species and per dataset.
cap-4000 was withdrawn before any fit: its rows are named withdrawn and
none of its conditions is read.

`print` prints the page's tables from <out dir> alone, figures rounded to
3 decimals.
"""

import csv
import hashlib
import importlib.util
import json
import re
import statistics
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name):
    spec = importlib.util.spec_from_file_location(f"c2_cap_{name}_for_page", HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score = _load("score")
RUN_CHECKS = tuple(_load("run_checks").CHECKS)
C3 = score.C3
LEARNERS = score.LEARNERS
MEASURES = score.MEASURES
PRIMARY = C3.PRIMARY
MARGIN = C3.MARGIN
HOLM_ALPHA = C3.HOLM_ALPHA
DEFAULT, TWIN, ALL, MID = "cap-2000", "twin", "cap-all", "cap-4000"
WITHDRAWN = (MID,)
WIDER = (MID, ALL)
FEW = 20
CAP_GATE = 4096
NOT_SCORED = ("none-added", "same-cohort")
UNAVAILABLE = "lineage split unavailable"
SCORE_FILES = ("checks.tsv", "draws.tsv", "datasets.tsv")
COUNT_FILES = ("count.tsv", "conditions.tsv")
SCORED_SAME = "the scored genomes of each outer fold are the default's"
STUDY_SAME = "the test folds are the default's, with at least its training genomes"
GATES = [
    (1, "inputs", "the run's checks printout holds one pass line for each of run_checks.py's checks, in its order, "
                  "and no line starting FAIL; the count directory holds count.tsv and conditions.tsv; each species' "
                  "folder holds score/checks.tsv, score/draws.tsv, score/datasets.tsv, results/added-labelled.tsv, "
                  "results/cohorts.json, expect.tsv and binary.txt"),
    (2, "species", "the runs directory holds every species of lists.tsv's first run once, and no other"),
    (3, "first run", "each species' expect.tsv lists exactly its first-run datasets, and no table holds a dataset of M "
                     "or of another species; species and antibiotic equal lists.tsv's"),
    (4, "cohorts", "each cohorts.json is its species', every cohort holds at most 4,096 genomes, and each wider "
                   "cohort holds the default's whole list"),
    (5, "binary", "each binary.txt names a commit and a code sha256; the species that share them are named together"),
    (6, "none-added", "a none-added or same-cohort dataset is in no arm's rows and no denominator, on either split"),
    (7, "withdrawn", "cap-4000 holds no row; its five lineage pairs stay in the family at p = 1, named withdrawn, and "
                     "none of its conditions is read"),
    (8, "twin", "the twin's rows are the default's fit at slots no draw has, and no condition reads the twin as an arm"),
    (9, "study", "a study row stays on the study split, with the study role and no slot; no lineage row has the "
                 "study role"),
    (10, "unavailable", "a dataset whose lineage split is unavailable is unscored alike on every arm's lineage rows, "
                        "and its study rows are scored as usual"),
    (11, "missing draw", "every row expect.tsv lists for scoring has its five learners' draws; an unscored draw or "
                         "dataset holds no figure; datasets.tsv is the mean of draws.tsv"),
    (12, "scored genomes", "score.py's checks all passed, and each scored cap-all row has its check that its scored "
                           "genomes are the default's"),
    (13, "few added", "the datasets that add fewer than 20 labelled genomes under cap-all are in count.tsv, which is "
                      "not the count without them"),
    (14, "count", "count.tsv's lineage denominators hold lineage rows only, and count.tsv and conditions.tsv are "
                  "what score.py --count --withdrawn cap-4000 writes from the species' datasets.tsv files"),
]
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
PAIR_WORDS = {f"{MID}-{DEFAULT}": "cap-4000 less cap-2000", f"{ALL}-{DEFAULT}": "cap-all less cap-2000",
              f"{TWIN}-{DEFAULT}": "twin less cap-2000"}


class Stop(Exception):
    pass


# ---------- reading files ----------

def load(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_json(path):
    with open(path, encoding="utf-8") as handle:
        return json.load(handle)


def num(x):
    return None if x in ("", None) else float(x)


def whole(x):
    return None if x in ("", None) else int(x)


def r12(x):
    return None if x is None else round(x, 12)


def same_cell(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 2e-12
    except (TypeError, ValueError):
        return False


def same_rows(got, want, columns):
    return len(got) == len(want) and all(same_cell(g.get(c, ""), w.get(c, "")) for g, w in zip(got, want) for c in columns)


def first_run(lrows):
    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]


def binary_of(text):
    """{commit, grm, code} of a binary.txt; a restart line is kept apart."""
    out = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("commit", "grm", "code") and parts[0] not in out:
            out[parts[0]] = parts[1]
    return out


def read_species(folder):
    """One species' files as the page reads them; optional files are None when absent."""
    res = folder / "results"
    run = {name: load(folder / "score" / name) for name in SCORE_FILES}
    run["added"] = load(res / "added-labelled.tsv")
    run["cohorts"] = load_json(res / "cohorts.json")
    run["expect"] = load(folder / "expect.tsv")
    run["binary"] = binary_of((folder / "binary.txt").read_text(encoding="utf-8"))
    run["fetch"] = load_json(res / "fetch.json") if (res / "fetch.json").is_file() else None
    run["folds"] = [r for p in sorted(res.glob("folds-lineage-*.tsv")) for r in load(p)]
    run["study_folds"] = load(res / "folds-study.tsv") if (res / "folds-study.tsv").is_file() else []
    run["study_files"] = (load(res / f"study-files-{ALL}.tsv") if (res / f"study-files-{ALL}.tsv").is_file() else [])
    timing = folder / "logs" / "timing.tsv"
    run["timing"] = load(timing) if timing.is_file() else []
    run["refused"] = {}
    for p in sorted((folder / "logs").glob("folds-*.log")) if (folder / "logs").is_dir() else []:
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
            m = re.match(r"^(\S+): lineage split unavailable: (\d+) default slot plans refused", line)
            if m:
                run["refused"][m.group(1)] = int(m.group(2))
    return run


# ---------- rules the gates and tables share ----------

def none_added(run):
    """{(arm, dataset)} that cohorts.py marks with no added genome labelled."""
    return {(r["arm"], r["dataset"]) for r in run["added"] if r["none_added"] == "1"}


def not_scored(run):
    """{(dataset, split, arm)} that expect.tsv lists as not scored, and every wider arm's rows cohorts.py marks
    none-added, on both splits."""
    out = {(e["dataset"], e["split"], e["arm"]) for e in run["expect"] if e["role"] in NOT_SCORED}
    for arm, d in none_added(run):
        for split in ("lineage", "study"):
            out.add((d, split, arm))
    return out


def few_added(runs, names):
    """[row] of the first-run datasets fitted under cap-all that add fewer than FEW labelled genomes, by rule from
    added-labelled.tsv, named in full."""
    out = []
    for slug in sorted(runs):
        for r in runs[slug]["added"]:
            if r["arm"] != ALL or r["dataset"] not in names or r["none_added"] == "1":
                continue
            if int(r["added_labelled"]) < FEW:
                sp, ab = names[r["dataset"]]
                out.append({"dataset": r["dataset"], "species": sp, "antibiotic": ab,
                            "added_labelled": int(r["added_labelled"]), "added_resistant": int(r["added_resistant"]),
                            "added_susceptible": int(r["added_susceptible"])})
    return out


def unavailable(run):
    """{dataset} whose lineage rows score.py scored 0 for an unavailable lineage split."""
    return {r["dataset"] for r in run["draws.tsv"] if r["split"] == "lineage" and r["why"] == UNAVAILABLE}


def all_datasets(runs):
    return [r for slug in sorted(runs) for r in runs[slug]["datasets.tsv"]]


def without_few(rows, few):
    """The datasets.tsv rows less cap-all's rows of the few-added datasets: the count beside."""
    gone = {r["dataset"] for r in few}
    return [r for r in rows if not (r["arm"] == ALL and r["dataset"] in gone)]


def counted(rows):
    """score.py --count's tables of these rows, as its files hold them."""
    try:
        rows_ = C3.count(rows, WITHDRAWN)
    except SystemExit as e:
        raise Stop(f"gate 14 count: score.py's count refuses these rows: {e}") from None
    count_rows = [{c: score.C9.fmt(r.get(c)) for c in C3.COUNT_COLUMNS} for r in rows_]
    cond_rows = [{c: score.C9.fmt(r.get(c)) for c in score.CONDITION_COLUMNS}
                 for r in score.conditions(rows_, WITHDRAWN)]
    return count_rows, cond_rows


def typed_draw(r):
    """A draws.tsv row as C3.dataset_rows takes it."""
    out = dict(r)
    out["scored"] = int(r["scored"])
    out["constant"] = whole(r.get("constant")) or 0
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = whole(r.get(k)) or 0
    return out


# ---------- gates ----------

def check_run_checks(checks_path):
    """Gate 1's first part: run_checks.py's printout says each of its checks passed, one line each. Raises Stop."""
    lines = Path(checks_path).read_text(encoding="utf-8").splitlines()
    passed = [line[len("pass "):].split(":", 1)[0] for line in lines if line.startswith("pass ")]
    if len(lines) != len(RUN_CHECKS) or passed != list(RUN_CHECKS) or any(line.startswith("FAIL") for line in lines):
        raise Stop(f"gate 1 inputs: the run's checks printout does not say that each of run_checks.py's "
                   f"{len(RUN_CHECKS)} checks passed")


def check_reading(count_dir, runs_dir, lrows):
    """Every gate of GATES, in order. Returns (runs, count rows, condition rows, beside count files); raises Stop
    naming the gate."""
    count_dir, runs_dir = Path(count_dir), Path(runs_dir)
    run_rows = first_run(lrows)
    names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in run_rows}
    in_lists = {r["dataset"]: (r["species"], r["antibiotic"]) for r in lrows}
    species = sorted({r["species_slug"] for r in run_rows})
    by_species = {s: sorted(r["dataset"] for r in run_rows if r["species_slug"] == s) for s in species}
    species_name = {r["species_slug"]: r["species"] for r in run_rows}

    # 1 inputs
    for name in COUNT_FILES:
        if not (count_dir / name).is_file():
            raise Stop(f"gate 1 inputs: the count directory has no {name}")
    if not runs_dir.is_dir():
        raise Stop("gate 1 inputs: no runs directory")
    for folder in sorted(p for p in runs_dir.iterdir() if p.is_dir()):
        for need in [f"score/{n}" for n in SCORE_FILES] + ["results/added-labelled.tsv", "results/cohorts.json",
                                                            "expect.tsv", "binary.txt"]:
            if not (folder / need).is_file():
                raise Stop(f"gate 1 inputs: {folder.name}/{need} is missing")

    # 2 species
    found = sorted(p.name for p in runs_dir.iterdir() if p.is_dir())
    if found != species:
        missing, more = sorted(set(species) - set(found)), sorted(set(found) - set(species))
        raise Stop(f"gate 2 species: {len(found)} species folders for the {len(species)} of the first run; "
                   f"missing {missing}, more {more}")
    runs = {slug: read_species(runs_dir / slug) for slug in found}
    rows = all_datasets(runs)
    seen = {}
    for r in rows:
        seen.setdefault(r["dataset"], set()).add(r["species_slug"])
    if any(len(v) > 1 for v in seen.values()):
        raise Stop("gate 2 species: a dataset's rows come from two species")

    # 3 first run
    for slug, run in runs.items():
        listed = sorted({e["dataset"] for e in run["expect"]})
        if listed != by_species[slug]:
            raise Stop(f"gate 3 first run: {slug}'s expect.tsv lists other datasets than its first run; "
                       f"more {sorted(set(listed) - set(by_species[slug]))}, "
                       f"less {sorted(set(by_species[slug]) - set(listed))}")
        for table in ("datasets.tsv", "draws.tsv"):
            for r in run[table]:
                if r["dataset"] not in by_species[slug]:
                    raise Stop(f"gate 3 first run: {slug}/{table} holds {r['dataset']}, not of its first run "
                               "(in M, or another species')")
                if (r["species"], r["antibiotic"]) != in_lists.get(r["dataset"]):
                    raise Stop(f"gate 3 first run: {slug}/{table}: {r['dataset']} is named otherwise than in lists.tsv")

    # 4 cohorts
    for slug, run in runs.items():
        doc = run["cohorts"]
        if doc.get("slug") != slug or doc.get("species") != species_name[slug]:
            raise Stop(f"gate 4 cohorts: {slug}/results/cohorts.json is another species'")
        arms = doc.get("arms", {})
        if DEFAULT not in arms or ALL not in arms:
            raise Stop(f"gate 4 cohorts: {slug}/results/cohorts.json has no {DEFAULT} or {ALL}")
        for arm, a in arms.items():
            if int(a["genomes"]) > CAP_GATE:
                raise Stop(f"gate 4 cohorts: {slug}: {arm} holds {int(a['genomes']):,} genomes, past {CAP_GATE:,}")
            if arm != DEFAULT and int(a["holds_default"]) != int(arms[DEFAULT]["genomes"]):
                raise Stop(f"gate 4 cohorts: {slug}: {arm} drops a genome of the default's list")

    # 5 binary
    for slug, run in runs.items():
        if not run["binary"].get("commit") or not run["binary"].get("code"):
            raise Stop(f"gate 5 binary: {slug}/binary.txt names no commit or no code sha256")

    # 6 none-added
    for slug, run in runs.items():
        gone = not_scored(run)
        for table in ("datasets.tsv", "draws.tsv"):
            for r in run[table]:
                if (r["dataset"], r["split"], r["arm"]) in gone:
                    raise Stop(f"gate 6 none-added: {slug}/{table} holds {r['arm']} on {r['dataset']} "
                               f"({r['split']}), which is none-added or same-cohort there")

    # 7 withdrawn
    got_count, got_cond = load(count_dir / "count.tsv"), load(count_dir / "conditions.tsv")
    for slug, run in runs.items():
        for table in ("expect", "draws.tsv", "datasets.tsv"):
            if any(r["arm"] == MID for r in run[table]):
                raise Stop(f"gate 7 withdrawn: {slug}'s {table} holds a {MID} row")
    pair = f"{MID}-{DEFAULT}"
    mine = [r for r in got_count if r["pair"] == pair]
    family = [r for r in mine if r["split"] == "lineage" and r["measure"] == PRIMARY]
    if (len(family) != len(LEARNERS) or any(r["withdrawn"] != "1" or r["entered"] != "0" or num(r["sign_p"]) != 1.0
                                             or r["holm_p"] in ("", None) for r in family)
            or any(r["withdrawn"] != "1" for r in mine)):
        raise Stop(f"gate 7 withdrawn: count.tsv does not hold {MID}'s five lineage pairs at p = 1, named withdrawn")
    if any(r["withdrawn"] == "1" for r in got_count if r["pair"] != pair):
        raise Stop("gate 7 withdrawn: count.tsv names another arm withdrawn")
    cond_mid = [r for r in got_cond if r["arm"] == MID]
    if not cond_mid or any(r["withdrawn"] != "1" or r["met"] != "" or r["learners"] != "" for r in cond_mid):
        raise Stop(f"gate 7 withdrawn: conditions.tsv reads a condition of {MID}")

    # 8 twin
    for slug, run in runs.items():
        draw_slots = {e["slot"] for e in run["expect"] if e["role"] == "draw"}
        for table in ("expect", "draws.tsv"):
            for r in run[table]:
                if r["arm"] == TWIN and (r["fit"] != DEFAULT or r["role"] != "twin" or r["slot"] in draw_slots
                                         or r["split"] != "lineage"):
                    raise Stop(f"gate 8 twin: {slug}/{table}: a twin row of {r['dataset']} is not the default's fit "
                               "at a slot of its own")
    if any(r["arm"] == TWIN for r in got_cond):
        raise Stop("gate 8 twin: conditions.tsv reads the twin as an arm")

    # 9 study
    for slug, run in runs.items():
        for r in run["draws.tsv"]:
            if (r["split"] == "study") != (r["role"] == "study") or (r["split"] == "study") != (r["slot"] == "-"):
                raise Stop(f"gate 9 study: {slug}/draws.tsv: {r['dataset']} {r['arm']} has a study row off the "
                           "study split or a lineage row with a study role")

    # 10 unavailable
    for slug, run in runs.items():
        gone = unavailable(run)
        for d in gone:
            lin = [r for r in run["draws.tsv"] if r["dataset"] == d and r["split"] == "lineage"]
            if any(r["why"] != UNAVAILABLE or r["scored"] != "0" for r in lin):
                raise Stop(f"gate 10 unavailable: {slug}: {d}'s lineage split is unavailable on one arm and scored "
                           "on another")
            if any(r["scored"] != "0" for r in run["datasets.tsv"] if r["dataset"] == d and r["split"] == "lineage"):
                raise Stop(f"gate 10 unavailable: {slug}: {d}'s lineage rows are scored in datasets.tsv")
            study = [r for r in run["draws.tsv"] if r["dataset"] == d and r["split"] == "study"]
            if any(r["why"] == UNAVAILABLE for r in study):
                raise Stop(f"gate 10 unavailable: {slug}: {d}'s study rows are read as lineage-unavailable")

    # 11 missing draw
    for slug, run in runs.items():
        gone = not_scored(run)
        have = {}
        for r in run["draws.tsv"]:
            have.setdefault((r["dataset"], r["split"], r["arm"], r["slot"]), []).append(r["learner"])
        for e in run["expect"]:
            if (e["dataset"], e["split"], e["arm"]) in gone or e["role"] in NOT_SCORED:
                continue
            if sorted(have.get((e["dataset"], e["split"], e["arm"], e["slot"]), [])) != sorted(LEARNERS):
                raise Stop(f"gate 11 missing draw: {slug}: {e['arm']} {e['dataset']} {e['split']} slot {e['slot']} "
                           "has not its five learners' draws")
        for table in ("draws.tsv", "datasets.tsv"):
            for r in run[table]:
                if r["scored"] != "1" and any(r.get(m) not in ("", None) for m in MEASURES):
                    raise Stop(f"gate 11 missing draw: {slug}/{table}: an unscored row of {r['dataset']} holds a figure")
        again = {(r["dataset"], r["split"], r["arm"], r["learner"]): r
                 for r in C3.dataset_rows([typed_draw(r) for r in run["draws.tsv"]])}
        got = run["datasets.tsv"]
        if len(got) != len(again):
            raise Stop(f"gate 11 missing draw: {slug}/datasets.tsv holds {len(got)} rows, not the {len(again)} of its draws")
        for r in got:
            g = again.get((r["dataset"], r["split"], r["arm"], r["learner"]))
            if g is None or any(not same_cell(r[c], score.C9.fmt(g.get(c))) for c in C3.DATASET_COLUMNS):
                raise Stop(f"gate 11 missing draw: {slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']} is not "
                           "the mean of its draws")

    # 12 scored genomes
    for slug, run in runs.items():
        failed = [c for c in run["checks.tsv"] if c["passed"] != "1"]
        if failed:
            raise Stop(f"gate 12 scored genomes: {slug}/checks.tsv holds {len(failed)} failed checks")
        gone = unavailable(run)
        passed = {(c["row"], c["check"]) for c in run["checks.tsv"]}
        for e in run["expect"]:
            if e["arm"] in (DEFAULT, TWIN) or e["role"] in NOT_SCORED:
                continue
            if e["split"] == "lineage" and e["dataset"] not in gone:
                where, check = f"{e['fit']}/{e['dataset']} lineage slot {e['slot']}", SCORED_SAME
            elif e["split"] == "study":
                where, check = f"{e['fit']}/{e['dataset']} study slot {e['slot']}", STUDY_SAME
            else:
                continue
            if (where, check) not in passed:
                raise Stop(f"gate 12 scored genomes: {slug}: {where} has no passed check that its scored genomes "
                           "are the default's")

    # 13 few added
    few = few_added(runs, names)
    fitted = {(r["dataset"], r["learner"]) for r in rows if r["split"] == "lineage" and r["arm"] == ALL}
    for r in few:
        if any((r["dataset"], l) not in fitted for l in LEARNERS):
            raise Stop(f"gate 13 few added: {r['species']} {r['antibiotic']} adds {r['added_labelled']} labelled "
                       "genomes and is not in cap-all's rows")
    if few:
        beside, _ = counted(without_few(rows, few))
        if same_rows(got_count, beside, C3.COUNT_COLUMNS):
            raise Stop("gate 13 few added: count.tsv is the count without the datasets that add fewer than 20 "
                       "labelled genomes")

    # 14 count
    for r in got_count:
        if r["split"] != "lineage":
            continue
        arm = r["pair"].rsplit(f"-{DEFAULT}", 1)[0]
        lineage = {x["dataset"] for x in rows if x["split"] == "lineage" and x["arm"] == arm
                   and x["learner"] == r["learner"]}
        if int(r["datasets"]) != len(lineage):
            raise Stop(f"gate 14 count: count.tsv's lineage {r['pair']} {r['learner']} {r['measure']} counts "
                       f"{r['datasets']} datasets, not its {len(lineage)} lineage datasets")
    want_count, want_cond = counted(rows)
    if not same_rows(got_count, want_count, C3.COUNT_COLUMNS):
        raise Stop("gate 14 count: count.tsv is not score.py's count of the species' datasets.tsv files")
    if not same_rows(got_cond, want_cond, score.CONDITION_COLUMNS):
        raise Stop("gate 14 count: conditions.tsv is not score.py's conditions of its count")
    others = sorted(p.name for p in count_dir.iterdir() if p.is_file() and p.name not in COUNT_FILES)
    return runs, got_count, got_cond, others


# ---------- the page's tables ----------

def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.py's conditions read them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    return lo is not None and lo > 0, holm is not None and holm < HOLM_ALPHA and int(row["better"]) > int(row["worse"])


def decides(r):
    """Why a conditions.tsv row decides or not."""
    if r["withdrawn"] == "1":
        return "withdrawn before any fit: not read"
    if r["reading"] != "ruling 15":
        return "decides nothing: the stricter reading"
    if r["condition"] == "5":
        return "decides nothing: it reads the study split"
    if r["condition"] == "6":
        return "the maintainer's"
    return "decides"


def page_conditions(conds):
    return [dict(r, decides=decides(r)) for r in conds]


CONDITION_PAGE_COLUMNS = score.CONDITION_COLUMNS + ["decides"]


DECIDING = ("1 and 2", "3", "4")


def page_reading(cond_rows):
    """cap-all's deciding conditions as conditions.tsv reads them under ruling 15, and whether all of them hold.
    Condition 5 reads the study split, which decides nothing in C2; 6 is the maintainer's."""
    at = {r["condition"]: r for r in cond_rows if r["arm"] == ALL and r["reading"] == "ruling 15"}
    out = [{"arm": ALL, "condition": c, "reads": at[c]["reads"], "met": at[c]["met"], "learners": at[c]["learners"]}
           for c in DECIDING]
    out.append({"arm": ALL, "condition": "all deciding", "reads": "conditions " + ", ".join(DECIDING) + " all met",
                "met": int(all(r["met"] == "1" for r in out)), "learners": ""})
    return out


READING_COLUMNS = ["arm", "condition", "reads", "met", "learners"]


def page_count(count_rows, beside=False):
    out = []
    for r in count_rows:
        row = dict(r, meets_1="", meets_2="")
        if r["split"] == "lineage" and r["measure"] == PRIMARY and r["withdrawn"] != "1":
            c1, c2 = meets(r)
            row["meets_1"], row["meets_2"] = int(c1), int(c2)
        if r["withdrawn"] == "1":
            row["reading"] = "withdrawn before any fit"
        elif beside:
            row["reading"] = "beside: decides nothing"
        elif r["split"] == "study":
            row["reading"] = "the study split: decides nothing"
        else:
            row["reading"] = "the count"
        out.append(row)
    return out


COUNT_PAGE_COLUMNS = C3.COUNT_COLUMNS + ["meets_1", "meets_2", "reading"]


def pairs(rows, arm, split="lineage", measure=PRIMARY):
    """[(dataset, learner, difference)] of one arm against the default, as C3.count takes them, and the pairs not
    scored."""
    cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}
    out, unscored = [], 0
    for (d, s, a, l), r in sorted(cell.items()):
        if s != split or a != arm:
            continue
        b = cell.get((d, s, DEFAULT, l))
        if b is None or r["scored"] != "1" or b["scored"] != "1" or r[measure] in ("", None) or b[measure] in ("", None):
            unscored += 1
            continue
        out.append((d, l, r12(float(r[measure]) - float(b[measure]))))
    return out, unscored


def tally(values, measure=PRIMARY):
    up = -1 if measure in C3.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def cohort_sizes(run):
    """{arm: genomes} as fitted: cohorts.json's, less the added genomes fetch_check took out of every wider cohort
    when cohorts.json was copied before it."""
    arms = run["cohorts"]["arms"]
    left = (run["fetch"] or {}).get("left_wider_cohorts", {}) or {}
    out = {}
    for arm, a in arms.items():
        n = int(a["genomes"])
        if "fetch_failed" not in a and arm != DEFAULT:
            n -= int(left.get(arm, 0))
        out[arm] = n
    return out


def short_because(arm, size, failed):
    """Why cap-all's cohort as fitted is short of 4,096: the rule's own count (cohorts.json's genomes, with any
    fetch_failed it already took out put back) short of it, failed added genomes, or both, each named."""
    if size >= CAP_GATE:
        return ""
    why = []
    if int(arm["genomes"]) + int(arm.get("fetch_failed", 0)) < CAP_GATE:
        why.append("the rule ranks fewer genomes")
    if failed:
        why.append(f"{failed} added genomes failed the fetch")
    return " and ".join(why)


def page_species(runs, lrows):
    run_rows = first_run(lrows)
    out = []
    for slug in sorted(runs):
        run = runs[slug]
        mine = [r for r in run_rows if r["species_slug"] == slug]
        arms = run["cohorts"]["arms"]
        size = cohort_sizes(run)
        failed = int((run["fetch"] or {}).get("failed_added", 0) or 0)
        row = {"species": mine[0]["species"], "species_slug": slug, "datasets": len(mine),
               "study_datasets": sum(1 for r in mine if r["in_S"] == "1"), "ranked": run["cohorts"].get("ranked"),
               "genomes_cap_2000": size[DEFAULT], "genomes_cap_4000_withdrawn": size.get(MID),
               "genomes_cap_all": size[ALL], "added_cap_all": int(arms[ALL]["added"]) - (int(arms[ALL]["genomes"]) - size[ALL]),
               "fetch_failed_added": failed,
               "short_of_4096": CAP_GATE - size[ALL] if size[ALL] < CAP_GATE else 0,
               "short_because": short_because(arms[ALL], size[ALL], failed),
               "none_added": sum(1 for arm, d in none_added(run) if arm == ALL and d in {m["dataset"] for m in mine})}
        rows = run["datasets.tsv"]
        for arm in (ALL, TWIN):
            got, unscored = pairs(rows, arm)
            b, w, s = tally([x for _, _, x in got])
            key = "cap_all" if arm == ALL else "twin"
            row.update({f"{key}_pairs": len(got), f"{key}_not_scored": unscored, f"{key}_better": b,
                        f"{key}_worse": w, f"{key}_same": s})
        row["reading"] = "decides nothing"
        out.append(row)
    return out


SPECIES_COLUMNS = ["species", "species_slug", "datasets", "study_datasets", "ranked", "genomes_cap_2000",
                   "genomes_cap_4000_withdrawn", "genomes_cap_all", "added_cap_all", "fetch_failed_added",
                   "short_of_4096", "short_because", "none_added"] + [
    f"{k}_{c}" for k in ("cap_all", "twin") for c in ("pairs", "not_scored", "better", "worse", "same")] + ["reading"]


def fit_genomes(run, dataset, arm):
    """The dataset's genomes in the arm's lineage fit, from the fold check's first slot and fold: the scored and
    the trained, and for a wider arm the added genomes kept out."""
    rows = sorted((r for r in run["folds"] if r["dataset"] == dataset and r["arm"] == arm),
                  key=lambda r: (int(r["slot"]), int(r["fold"])))
    if not rows:
        return None
    r = rows[0]
    return sum(int(r[k]) for k in ("scored", "trained", "kept_out_lineage", "kept_out_duplicate"))


def role_of(run, dataset, split, arm, gone):
    if split == "lineage" and dataset in gone:
        return UNAVAILABLE
    roles = {e["role"] for e in run["expect"] if e["dataset"] == dataset and e["split"] == split and e["arm"] == arm}
    if (arm, dataset) in none_added(run):
        return "none-added"
    return ",".join(sorted(roles)) or "not listed"


def page_datasets(runs, lrows):
    meta = {r["dataset"]: r for r in lrows}
    out = []
    for slug in sorted(runs):
        run = runs[slug]
        gone = unavailable(run)
        cells = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in run["datasets.tsv"]}
        added = {r["dataset"]: r for r in run["added"] if r["arm"] == ALL}
        train = {r["dataset"]: r for r in run["study_folds"] if r["arm"] == ALL}
        for d in sorted({e["dataset"] for e in run["expect"]}):
            m, a = meta[d], added.get(d, {})
            for split in ("lineage", "study"):
                if not any(e["dataset"] == d and e["split"] == split for e in run["expect"]):
                    continue
                for learner in LEARNERS:
                    base = cells.get((d, split, DEFAULT, learner), {})
                    wide = cells.get((d, split, ALL, learner), {})
                    row = {"dataset": d, "species": m["species"], "antibiotic": m["antibiotic"], "split": split,
                           "learner": learner, "role": role_of(run, d, split, ALL, gone),
                           "draws": base.get("draws", ""), "default_draws_scored": base.get("scored_draws", ""),
                           "cap_all_draws_scored": wide.get("scored_draws", ""),
                           "default": num(base.get(PRIMARY)), "cap_all": num(wide.get(PRIMARY)),
                           "default_recall": num(base.get("recall")), "cap_all_recall": num(wide.get("recall")),
                           "genomes": int(m["samples"]),
                           "fit_genomes_cap_2000": fit_genomes(run, d, DEFAULT) if split == "lineage" else None,
                           "fit_genomes_cap_all": fit_genomes(run, d, ALL) if split == "lineage" else None,
                           "added_labelled": whole(a.get("added_labelled")),
                           "added_resistant": whole(a.get("added_resistant")),
                           "added_susceptible": whole(a.get("added_susceptible")),
                           "few_added": int(a.get("none_added") == "0" and int(a.get("added_labelled", 0)) < FEW),
                           "study_added_train": whole(train.get(d, {}).get("train_added")) if split == "study" else None,
                           "reading": "decides nothing"}
                    for k, x, y in (("difference", "cap_all", "default"), ("recall_difference", "cap_all_recall",
                                                                           "default_recall")):
                        row[k] = r12(row[x] - row[y]) if row[x] is not None and row[y] is not None else None
                    if split == "lineage":
                        twin = cells.get((d, split, TWIN, learner), {})
                        row["twin"] = num(twin.get(PRIMARY))
                        row["twin_difference"] = (r12(row["twin"] - row["default"])
                                                  if row["twin"] is not None and row["default"] is not None else None)
                    out.append(row)
    return out


DATASET_COLUMNS = ["dataset", "species", "antibiotic", "split", "learner", "role", "draws", "default_draws_scored",
                   "cap_all_draws_scored", "default", "cap_all", "difference", "twin", "twin_difference",
                   "default_recall", "cap_all_recall", "recall_difference", "genomes", "fit_genomes_cap_2000",
                   "fit_genomes_cap_all", "added_labelled", "added_resistant", "added_susceptible", "few_added",
                   "study_added_train", "reading"]


def page_not_fitted(runs, lrows):
    """The datasets an arm does not score, named in full with the reason."""
    meta = {r["dataset"]: r for r in lrows}
    out = []
    for slug in sorted(runs):
        run = runs[slug]
        listed = {e["dataset"] for e in run["expect"]}
        study = {e["dataset"] for e in run["expect"] if e["split"] == "study"}
        for arm, d in sorted(none_added(run)):
            if d in listed and arm == ALL:
                out.append({"dataset": d, "species": meta[d]["species"], "antibiotic": meta[d]["antibiotic"],
                            "arm": arm, "split": "lineage and study" if d in study else "lineage",
                            "reason": "none-added: no added genome labelled", "refused_slots": ""})
        for e in run["expect"]:
            if e["role"] in NOT_SCORED and (e["arm"], e["dataset"]) not in none_added(run):
                why = ("none-added: no added genome trains on the study split" if e["role"] == "none-added"
                       else "same-cohort: its cohort is " + e["fit"] + "'s")
                key = (e["dataset"], e["arm"], e["split"], why)
                if key not in {(r["dataset"], r["arm"], r["split"], r["reason"]) for r in out}:
                    out.append({"dataset": e["dataset"], "species": meta[e["dataset"]]["species"],
                                "antibiotic": meta[e["dataset"]]["antibiotic"], "arm": e["arm"], "split": e["split"],
                                "reason": why, "refused_slots": ""})
        for d in sorted(unavailable(run)):
            out.append({"dataset": d, "species": meta[d]["species"], "antibiotic": meta[d]["antibiotic"],
                        "arm": "every arm and the twin", "split": "lineage",
                        "reason": "lineage split unavailable: the default's lineage draws refused past grm arm's two "
                                  "replacements", "refused_slots": run["refused"].get(d, "")})
    return out


NOT_FITTED_COLUMNS = ["dataset", "species", "antibiotic", "arm", "split", "reason", "refused_slots"]


def page_few(few):
    return [dict(r, reading="in the deciding count; the count without them decides nothing") for r in few]


FEW_COLUMNS = ["dataset", "species", "antibiotic", "added_labelled", "added_resistant", "added_susceptible", "reading"]


def page_study(runs, lrows):
    """Per study dataset: the added genomes that may train under cap-all, and the role; deciding nothing."""
    meta = {r["dataset"]: r for r in lrows}
    out = []
    for slug in sorted(runs):
        run = runs[slug]
        files = {r["dataset"]: r for r in run["study_files"] if r["arm"] == ALL}
        folds = {r["dataset"]: r for r in run["study_folds"] if r["arm"] == ALL}
        for e in run["expect"]:
            if e["split"] != "study" or e["arm"] != ALL:
                continue
            f, k = files.get(e["dataset"], {}), folds.get(e["dataset"], {})
            out.append({"dataset": e["dataset"], "species": meta[e["dataset"]]["species"],
                        "antibiotic": meta[e["dataset"]]["antibiotic"], "role": e["role"],
                        "added_labelled": whole(f.get("added")), "may_train": whole(f.get("train_added")),
                        "trained": whole(k.get("train_added")), "unused_bioproject": whole(f.get("unused_bioproject")),
                        "unused_unknown": whole(f.get("unused_unknown")),
                        "unused_duplicate": whole(f.get("unused_duplicate")), "tested": whole(k.get("tested")),
                        "reading": "the study split: decides nothing"})
    return out


STUDY_COLUMNS = ["dataset", "species", "antibiotic", "role", "added_labelled", "may_train", "trained",
                 "unused_bioproject", "unused_unknown", "unused_duplicate", "tested", "reading"]


def page_checks(runs):
    """score.py's checks per name: passed over made, summed over the species."""
    out = {}
    for slug in sorted(runs):
        for c in runs[slug]["checks.tsv"]:
            r = out.setdefault(c["check"], {"check": c["check"], "made": 0, "passed": 0, "species": set()})
            r["made"] += 1
            r["passed"] += c["passed"] == "1"
            r["species"].add(slug)
    return [dict(r, species=len(r["species"])) for _, r in sorted(out.items())]


CHECK_COLUMNS = ["check", "passed", "made", "species"]


def page_binaries(runs):
    """Per commit and code sha256, its species named; and under each, every binary sha256 with its species (the
    lead's answer 3: two binaries under one code sha256 are named, not failed)."""
    groups, grms = {}, {}
    for slug in sorted(runs):
        b = runs[slug]["binary"]
        groups.setdefault((b["commit"], b["code"]), []).append(runs[slug]["cohorts"]["species"])
        grms.setdefault((b["commit"], b["code"]), {}).setdefault(b.get("grm", ""), []).append(
            runs[slug]["cohorts"]["species"])
    return [{"commit": c, "code": k, "species": len(v), "names": ", ".join(v),
             "binaries": len(grms[(c, k)]),
             "binary_names": "; ".join(f"{g}: {', '.join(s)}" for g, s in sorted(grms[(c, k)].items()))}
            for (c, k), v in sorted(groups.items())]


BINARY_COLUMNS = ["commit", "code", "species", "names", "binaries", "binary_names"]


FIT = re.compile(r"^(fit|study)-(cap-2000|cap-4000|cap-all)-(.+)$")


def page_times(runs):
    """Per arm and split: the fits' wall and CPU seconds and peak memory, as sums and medians. One fit serves the
    five learners, so a time per learner is not recorded."""
    cells = {}
    for slug in sorted(runs):
        for r in runs[slug]["timing"]:
            m = FIT.match(r["step"])
            if not m or r["wall_s"] in ("-", ""):
                continue
            split = "lineage" if m.group(1) == "fit" else "study"
            c = cells.setdefault((m.group(2), split), {"wall": [], "cpu": [], "peak": []})
            c["wall"].append(float(r["wall_s"]))
            c["cpu"].append(float(r["cpu_s"]))
            c["peak"].append(int(r["peak_kb"]))
    out = []
    for (arm, split), c in sorted(cells.items()):
        out.append({"arm": arm, "split": split, "fits": len(c["wall"]), "wall_s_sum": round(sum(c["wall"]), 2),
                    "wall_s_median": round(statistics.median(c["wall"]), 2), "cpu_s_sum": round(sum(c["cpu"]), 2),
                    "cpu_s_median": round(statistics.median(c["cpu"]), 2), "peak_kb_median": statistics.median(c["peak"]),
                    "peak_kb_max": max(c["peak"])})
    return out


TIME_COLUMNS = ["arm", "split", "fits", "wall_s_sum", "wall_s_median", "cpu_s_sum", "cpu_s_median", "peak_kb_median",
                "peak_kb_max"]


# ---------- writing ----------

README = """# The cap-arms reading: the record

This directory holds the tables behind [docs/cap-arms-result.md](../../docs/cap-arms-result.md): plan block C2, `cap-all` (the top 4,096 genomes a species by the release's selection rule) against `cap-2000` (the release's 2,000), scored on the default's held-out genomes, beside a twin, with `cap-4000` withdrawn before any fit. The first run holds {datasets} datasets in {species} species, {study} of them with a study split. Every row is an aggregate over a dataset, a species, a split, an arm or a learner. No file holds a genome, a BioProject, a sequence or a feature.

`scripts/c2-cap/score.py --count` wrote the files under `count/`, copied here as they are; `scripts/c2-cap/page_tables.py write` checked them against the scorer's own count of the species' `datasets.tsv` files and wrote every other file here from the run's `score/` and `results/` tables, after the {gates} gates of `gates.txt` passed. The rules of the reading were fixed before any C2 score was read (`docs/assumption-plan-corrections.md`, entries "C2 design", "C2 rulings" and "C2 review").

| File | Holds |
|---|---|
| `gates.txt` | the gates the program passed before it wrote a table |
| `checks.txt` | the printout of `scripts/c2-cap/run_checks.py`, the checks of the run's copies made before any score was read, counts only; gate 1 holds that each check passed |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the conditions as the scorer reads them{others} |
| `page-reading.tsv` | `cap-all`'s deciding conditions under ruling 15, 1 and 2, 3 and 4, and whether all are met; 6 is the maintainer's |
| `page-conditions.tsv` | the conditions per arm and reading, with which of them decide: ruling 15's 1 and 2, 3 and 4 for `cap-all`; condition 5 reads the study split and decides nothing; the stricter reading decides nothing; `cap-4000`'s rows are withdrawn and not read; 6 is the maintainer's |
| `page-count.tsv` | every row of `count.tsv`, with conditions 1 and 2 read on the lineage primary rows; study rows and `cap-4000`'s rows marked |
| `page-count-without-few.tsv` | beside it, deciding nothing: the same count with `cap-all`'s rows of the datasets that add fewer than 20 labelled genomes left out, Holm over its own rows |
| `page-few-added.tsv` | the datasets that add fewer than 20 labelled genomes under `cap-all`, named in full, with the added genomes by class; they stay in the count |
| `page-species.tsv` | per species, named in full: the cohorts' sizes, the added genomes, and the dataset and learner pairs better, worse and the same; deciding nothing |
| `page-datasets.tsv` | per dataset, split and learner, species and antibiotic in full: `cap-2000`'s and `cap-all`'s balanced accuracy and recall, the differences, the twin's, the draws scored, the genomes and the added genomes by class, and the role |
| `page-not-fitted.tsv` | the datasets an arm does not score, named in full with the reason |
| `page-study.tsv` | per study dataset: the added genomes that may train under `cap-all` and those that trained; deciding nothing |
| `page-checks.tsv` | the scorer's checks per name, passed over made |
| `page-binaries.tsv` | each binary by commit and code sha256, with its species named, and under it each binary sha256 with its species |
| `page-times.tsv` | the fits' wall and CPU seconds and peak memory per arm and split, as sums and medians |

A dataset's figure is the mean of its draws: `grm arm` slots 0 to 2 for each arm and 3 to 5 for the twin, at evaluation seed 42; the study split has one fit per dataset and arm. An added genome only trains and is never scored, so every arm is read on the default's held-out genomes. A split by lineage or study is a retrospective split inside one database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c2-cap/page_tables.py print` prints the page's tables from these files.
"""


def write_tsv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(score.C9.fmt(r.get(c)) for c in columns) + "\n")


def write(count_dir, runs_dir, l_path, checks_path, out):
    count_dir, out = Path(count_dir), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    check_run_checks(checks_path)
    runs, count_rows, cond_rows, others = check_reading(count_dir, runs_dir, lrows)
    names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in first_run(lrows)}
    few = few_added(runs, names)
    rows = all_datasets(runs)
    beside, _ = counted(without_few(rows, few))
    out.mkdir(parents=True, exist_ok=True)
    (out / "count").mkdir()
    for name in COUNT_FILES + tuple(others):
        (out / "count" / name).write_bytes((count_dir / name).read_bytes())
    (out / "checks.txt").write_bytes(Path(checks_path).read_bytes())
    write_tsv(out / "page-reading.tsv", page_reading(cond_rows), READING_COLUMNS)
    write_tsv(out / "page-conditions.tsv", page_conditions(cond_rows), CONDITION_PAGE_COLUMNS)
    write_tsv(out / "page-count.tsv", page_count(count_rows), COUNT_PAGE_COLUMNS)
    write_tsv(out / "page-count-without-few.tsv", page_count(beside, beside=True), COUNT_PAGE_COLUMNS)
    write_tsv(out / "page-few-added.tsv", page_few(few), FEW_COLUMNS)
    write_tsv(out / "page-species.tsv", page_species(runs, lrows), SPECIES_COLUMNS)
    write_tsv(out / "page-datasets.tsv", page_datasets(runs, lrows), DATASET_COLUMNS)
    write_tsv(out / "page-not-fitted.tsv", page_not_fitted(runs, lrows), NOT_FITTED_COLUMNS)
    write_tsv(out / "page-study.tsv", page_study(runs, lrows), STUDY_COLUMNS)
    write_tsv(out / "page-checks.tsv", page_checks(runs), CHECK_COLUMNS)
    write_tsv(out / "page-binaries.tsv", page_binaries(runs), BINARY_COLUMNS)
    write_tsv(out / "page-times.tsv", page_times(runs), TIME_COLUMNS)
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    run_rows = first_run(lrows)
    beside_files = "".join(f"; `count/{n}`, written beside by the scorer, decides nothing" for n in others)
    (out / "README.md").write_text(README.format(datasets=len(run_rows), species=len(runs), gates=len(GATES),
                                                 study=sum(1 for r in run_rows if r["in_S"] == "1"),
                                                 others=beside_files), encoding="utf-8")
    files = sorted(p.relative_to(out).as_posix() for p in out.rglob("*") if p.is_file())
    with open(out / "SHA256SUMS", "w", encoding="utf-8") as handle:
        for name in files:
            handle.write(f"{hashlib.sha256((out / name).read_bytes()).hexdigest()}  {name}\n")
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
    p = float(x)
    return "< 0.001" if p < 0.001 else d3(p)


def count_words(x):
    return "none" if x in ("", None) else f"{int(float(x)):,}"


def of(part, whole_):
    return f"{count_words(part)} of {count_words(whole_)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def yes(x):
    return {"1": "met", "0": "not met", "": "not read"}[str(x)]


def table(header, rows):
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def print_conditions(rows):
    return table(["Arm", "Reading", "Condition", "Reads", "Met", "Learners", "Decides"],
                 [[f"`{r['arm']}`", r["reading"], r["condition"], r["reads"],
                   "the maintainer's" if r["condition"] == "6" else yes(r["met"]),
                   r["learners"] or "none", r["decides"]] for r in rows])


def print_reading(rows):
    return table(["Arm", "Condition", "Reads", "Met", "Learners"],
                 [[f"`{r['arm']}`", r["condition"], r["reads"], yes(r["met"]), r["learners"] or "none"] for r in rows])


def print_family(rows):
    out = []
    for r in rows:
        if r["split"] != "lineage" or r["measure"] != PRIMARY:
            continue
        if r["withdrawn"] == "1":
            out.append([PAIR_WORDS[r["pair"]], f"`{r['learner']}`", "withdrawn", "withdrawn", "none", "none", "none",
                        "none", "none", pval(r["sign_p"]), pval(r["holm_p"])])
            continue
        out.append([PAIR_WORDS[r["pair"]], f"`{r['learner']}`", yes(r["meets_1"]), yes(r["meets_2"]),
                    of(r["better"], r["entered"]), of(r["worse"], r["entered"]), of(r["same"], r["entered"]),
                    d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"]),
                    pval(r["holm_p"])])
    return table(["Pair", "Learner", "Condition 1", "Condition 2", "Better", "Worse", "Same", "Median", "Interval",
                  "Sign p", "Holm p"], out)


def print_count(rows, measures=MEASURES):
    out = []
    for r in rows:
        if r["withdrawn"] == "1" or r["measure"] not in measures:
            continue
        out.append([r["split"], PAIR_WORDS[r["pair"]], f"`{r['learner']}`", MEASURE_WORDS[r["measure"]], r["orientation"],
                    of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"], d3(r["median_difference"]),
                    interval(r["interval_lo"], r["interval_hi"]), r["species"], r["reading"]])
    return table(["Split", "Pair", "Learner", "Measure", "Better is", "Entered", "Better", "Worse", "Same", "Median",
                  "Interval", "Species", "Reading"], out)


def print_few(rows):
    return table(["Species", "Antibiotic", "Added labelled", "Resistant", "Susceptible"],
                 [[f"*{r['species']}*", r["antibiotic"], r["added_labelled"], r["added_resistant"],
                   r["added_susceptible"]] for r in rows])


def print_species(rows):
    out = []
    for r in rows:
        out.append([f"*{r['species']}*", r["datasets"], r["study_datasets"], count_words(r["genomes_cap_2000"]),
                    count_words(r["genomes_cap_all"]), count_words(r["added_cap_all"]),
                    r["short_because"] or "none", r["none_added"],
                    f"{r['cap_all_better']} / {r['cap_all_worse']} / {r['cap_all_same']} of {r['cap_all_pairs']}",
                    f"{r['twin_better']} / {r['twin_worse']} / {r['twin_same']} of {r['twin_pairs']}"])
    return table(["Species", "Datasets", "In S", "`cap-2000` genomes", "`cap-all` genomes", "Added", "Short of 4,096",
                  "None-added", "`cap-all` better / worse / same", "Twin better / worse / same"], out)


def print_datasets(rows):
    by = {}
    for r in rows:
        by.setdefault((r["species"], r["antibiotic"], r["split"]), []).append(r)
    out = []
    for (sp, ab, split), rs in by.items():
        first = rs[0]
        added = ("none" if first["added_labelled"] in ("", None) else
                 f"{count_words(first['added_labelled'])} ({first['added_resistant']} R, {first['added_susceptible']} S)")
        out.append([f"*{sp}*", ab, split, first["role"], count_words(first["genomes"]), added]
                   + [d3(r["difference"]) for r in rs])
    return table(["Species", "Antibiotic", "Split", "`cap-all` role", "Genomes", "Added labelled"]
                 + [f"`{l}` difference" for l in LEARNERS], out)


def print_not_fitted(rows):
    return table(["Species", "Antibiotic", "Arm", "Split", "Reason"],
                 [[f"*{r['species']}*", r["antibiotic"], r["arm"], r["split"], r["reason"]] for r in rows])


def print_checks(rows):
    return table(["Check", "Passed", "Species"], [[r["check"], of(r["passed"], r["made"]), r["species"]] for r in rows])


def print_binaries(rows):
    return table(["Commit", "Code sha256", "Species", "Binary sha256s"],
                 [[f"`{r['commit'][:7]}`", f"`{r['code'][:12]}`", r["names"],
                   "; ".join(f"`{x.split(': ', 1)[0][:12]}`: {x.split(': ', 1)[1]}" if ": " in x else x
                             for x in r.get("binary_names", "").split("; ") if x) or "none"] for r in rows])


def print_times(rows):
    return table(["Arm", "Split", "Fits", "Wall s, sum", "Wall s, median", "CPU s, sum", "CPU s, median",
                  "Peak KB, median", "Peak KB, largest"],
                 [[f"`{r['arm']}`", r["split"], r["fits"], r["wall_s_sum"], r["wall_s_median"], r["cpu_s_sum"],
                   r["cpu_s_median"], count_words(r["peak_kb_median"]), count_words(r["peak_kb_max"])] for r in rows])


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    return {
        "reading": print_reading(load(folder / "page-reading.tsv")),
        "conditions": print_conditions(load(folder / "page-conditions.tsv")),
        "family": print_family(count),
        "count": print_count(count),
        "without few": print_count(load(folder / "page-count-without-few.tsv"), (PRIMARY, "recall")),
        "few added": print_few(load(folder / "page-few-added.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "not fitted": print_not_fitted(load(folder / "page-not-fitted.tsv")),
        "datasets": print_datasets([r for r in load(folder / "page-datasets.tsv")]),
        "checks": print_checks(load(folder / "page-checks.tsv")),
        "binaries": print_binaries(load(folder / "page-binaries.tsv")),
        "times": print_times(load(folder / "page-times.tsv")),
    }


def main(argv):
    if len(argv) == 6 and argv[0] == "write":
        try:
            files = write(argv[1], argv[2], argv[3], argv[4], argv[5])
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
