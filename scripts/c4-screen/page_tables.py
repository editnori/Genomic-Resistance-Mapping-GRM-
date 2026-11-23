"""The tables of docs/screen-arms-result.md: plan block C4's one reading, as files and as the page prints them.

  python3 -I page_tables.py write <reading dir> <lists.tsv> <checks.txt> <out dir>
  python3 -I page_tables.py print <out dir>

<lists.tsv> is scripts/c4-screen/lists.tsv, the list of L: the first run
is its rows with first_run 1, which are the rows not in M. <checks.txt> is
the printout of scripts/c4-screen/run_checks.py over the run's copies, made
before any score was read. The reading directory holds:

  reading.txt            the reading's commands, each with the sha256 of what it wrote
  count/                 score.py --count's count.tsv, count-no-constant.tsv,
                         conditions.tsv and count-floor-all.tsv, and
                         withdrawn.txt: the arms named with --withdrawn, one a
                         line, empty when none
  score/<species slug>/  one species' run as run_checks.py gather copies it:
                         score.py's checks.tsv, draws.tsv and datasets.tsv,
                         fit_times.py's fits.tsv (the run's logs/fits.tsv) and
                         the run's binary.txt; never score.log

`write` passes every gate of GATES before it writes a file, and stops at
the first that fails, naming it. It then writes <out dir> (absent or
empty): score.py's four count tables and withdrawn.txt under count/, copied
once the count gate has held them to score.py's own count of the
datasets.tsv files; reading.txt and checks.txt, copied; the page's tables
as page-*.tsv; gates.txt, README.md and SHA256SUMS. Where score.py's tables
already are the page's (the count, the count beside it and the
conditions), this script checks and copies them, and its page-*.tsv files
add only marks: which rows decide, which stand in the Holm family, and
conditions 1 and 2 read on each family row. Every row is an aggregate over
a dataset, a species, an arm or a learner: no file holds a genome, a
sequence or a feature.

`print` prints the page's tables from <out dir> alone, figures rounded to
3 decimals.
"""

import csv
import hashlib
import importlib.util
import re
import statistics
import sys
from decimal import ROUND_HALF_EVEN, Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


score = _load("c4_screen_score_for_page", HERE / "score.py")
fit_times = _load("c4_screen_fit_times_for_page", HERE / "fit_times.py")

LEARNERS = list(score.LEARNERS)
MEASURES = list(score.MEASURES)
PRIMARY = score.PRIMARY
MARGIN = score.MARGIN
DEFAULT, TWIN = score.DEFAULT, score.TWIN
ARMS, SCREEN_ARMS, FLOOR_ARM = tuple(score.ARMS), tuple(score.SCREEN_ARMS), score.FLOOR_ARM
SCREENED = tuple(score.SCREENED)
ALL, MOVES = score.ALL, score.MOVES
SEEDS, TWIN_SEEDS, STUDY_SEED = (42, 43, 44), (45, 46, 47), 42
COUNT_FILES = ("count.tsv", "count-no-constant.tsv", "conditions.tsv", "count-floor-all.tsv")
SPECIES_FILES = ("checks.tsv", "draws.tsv", "datasets.tsv", "fits.tsv", "binary.txt")
REFUSED_WHY = "constant model: refused for no candidate"
# The count's parts per split, (arm, set), as score.py --count writes them.
PARTS = {"lineage": list(score.FAMILY_PARTS), "study": [("cmh-4000", ALL), ("cmh-16000", ALL), (FLOOR_ARM, MOVES)]}
# The deciding readings of conditions.tsv; "stricter" stands beside and decides nothing.
SCREEN_READING = "the learners it touches, every dataset"
FLOOR_READING = "ruling 15, 800 genomes or more"
GATES = [
    (1, "files", "the reading holds reading.txt, count/ with score.py --count's four tables and withdrawn.txt, and "
                 "score/; no score.log is in it; the run's checks.txt says every check passed"),
    (2, "species", "score/ holds each species of lists.tsv's first run once and no other, each with its five files"),
    (3, "datasets", "each datasets.tsv holds exactly its species' first-run datasets, none in M: the default, each arm "
                    "not withdrawn and the twin on the lineage split, and the default and each arm not withdrawn on "
                    "the study split for those in S, each for the five learners; species, antibiotic and genome count "
                    "equal lists.tsv's, and the floor arm moves the datasets of 800 genomes or more"),
    (4, "study", "every study row is a study run at seed 42, and no study row stands on the lineage split"),
    (5, "twin", "the twin stands on the lineage split at seeds 45, 46 and 47 alone, counted as the pair twin-k31 over "
                "every dataset and over those of 800 genomes or more, and conditions.tsv reads it as no arm"),
    (6, "missing", "a draw not scored has no figure, and a dataset with a draw not scored is not scored and has no "
                   "figure: a missing draw is never read as zero"),
    (7, "constant", "each constant-model draw is an arm's study run refused for no candidate, scored and marked, and "
                    "count-no-constant.tsv counts none"),
    (8, "draws", "each draws.tsv holds each dataset's draws at seeds 42, 43 and 44, the twin's at 45, 46 and 47 and "
                 "the study run at seed 42, each for the five learners, and nothing else; datasets.tsv is their mean"),
    (9, "checks", "every check of each checks.tsv passed, and the figures agree: a screen arm's scm and cart, and "
                  "the floor arm below 800 genomes, equal the default's in datasets.tsv"),
    (10, "screen", "count.tsv and count-no-constant.tsv count a screen arm on rf, gbt and lasso only"),
    (11, "floor", "count.tsv counts the floor arm over the datasets of 800 genomes or more only, the ones lists.tsv "
                  "says it moves; count-floor-all.tsv holds the floor arm over every dataset with no Holm p, and "
                  "decides nothing"),
    (12, "withdrawn", "withdrawn.txt names arms of the block only, and the Holm family is the 21 lineage pairs "
                      "score.py holds, a withdrawn arm's at p = 1 before Holm"),
    (13, "binary", "each binary.txt names a commit, a binary and a code sha256; the page names each binary with its "
                   "species unless every species ran one code sha256"),
    (14, "conditions", "conditions 1, 2, 3 and 4 of each screen arm, and 1 and 2, 3 and 4 of the floor arm's ruling 15, "
                       "as read here off count.tsv, equal conditions.tsv's"),
    (15, "count", "count.tsv, count-no-constant.tsv, conditions.tsv and count-floor-all.tsv are what score.py --count "
                  "writes from the datasets.tsv files with withdrawn.txt's arms"),
    (16, "fits", "each fits.tsv is fit_times.py's table of the arms not withdrawn and names neither the default nor "
                 "the twin"),
]
GATE_NAMES = {n: name for n, name, _ in GATES}
MEASURE_WORDS = {"balanced_accuracy": "balanced accuracy", "recall": "recall", "specificity": "specificity",
                 "vme": "very major error rate", "me": "major error rate", "log_loss": "log loss"}
SET_WORDS = {ALL: "every dataset", MOVES: "800 genomes or more"}


class Stop(Exception):
    pass


def stop(n, text):
    raise Stop(f"gate {n} {GATE_NAMES[n]}: {text}")


# ---------- reading files ----------

def load(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def header(path):
    with open(path, encoding="utf-8") as handle:
        return handle.readline().rstrip("\n").split("\t")


def num(x):
    return None if x in ("", None) else float(x)


def r12(x):
    return None if x is None else round(x, 12)


def fmt(x):
    return score.C9.fmt(x)


def first_run(lrows):
    """The first run: L's rows with first_run 1, none of them in M."""
    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]


def touched(arm):
    """The learners a pair is counted on: rf, gbt and lasso for a screen arm, all five for the floor arm and the twin."""
    return list(SCREENED) if arm in SCREEN_ARMS else list(LEARNERS)


def read_withdrawn(path):
    return [w.strip() for w in Path(path).read_text(encoding="utf-8").splitlines() if w.strip()]


def expected(rows, slug, arms):
    """expect.py's rows for one species' first-run rows, as dicts: the run's results, draw by draw."""
    datasets = [r["dataset"] for r in rows]
    study = [r["dataset"] for r in rows if r["in_S"] == "1"]
    out = score.EXPECT.rows(slug, datasets, study, list(SEEDS), list(TWIN_SEEDS), list(arms))
    return [dict(zip(score.EXPECT.COLUMNS, map(str, r))) for r in out]


def typed_draw(r):
    """A draws.tsv row as score.dataset_rows takes it."""
    out = dict(r)
    out["seed"] = int(r["seed"])
    out["scored"] = int(r["scored"])
    out["constant"] = int(r["constant"] or 0)
    for m in MEASURES:
        out[m] = num(r.get(m))
    for k in ("n", "tp", "fp", "tn", "fn"):
        out[k] = int(r[k]) if r.get(k) not in ("", None) else 0
    return out


def same_cell(a, b):
    if a == b:
        return True
    try:
        return abs(float(a) - float(b)) <= 2e-12
    except ValueError:
        return False


def fmt_rows(rows, columns):
    return [{c: fmt(r.get(c)) for c in columns} for r in rows]


# ---------- gates ----------

def check_reading(reading, lrows, checks_path):
    """Every gate of GATES, in order. Returns what the tables read; raises Stop naming the first gate that fails."""
    reading = Path(reading)
    run = first_run(lrows)
    lists_of = {r["dataset"]: r for r in lrows}

    # 1 files
    if not (reading / "reading.txt").is_file():
        stop(1, "no reading.txt")
    for name in COUNT_FILES + ("withdrawn.txt",):
        if not (reading / "count" / name).is_file():
            stop(1, f"count/{name} is missing")
    if not (reading / "score").is_dir():
        stop(1, "score/ is missing")
    logs = sorted(p.relative_to(reading).as_posix() for p in reading.rglob("score.log"))
    if logs:
        stop(1, f"score.log stays out of the record, and the reading holds {', '.join(logs)}")
    lines = Path(checks_path).read_text(encoding="utf-8").splitlines()
    tail = re.fullmatch(r"CHECKS (\d+) of (\d+)", lines[-1]) if lines else None
    if not tail or tail.group(1) != tail.group(2) or any(line.startswith("FAIL") for line in lines):
        stop(1, "the run's checks.txt does not say that every check passed")

    # 2 species
    species = sorted({r["species_slug"] for r in run})
    found = sorted(p.name for p in (reading / "score").iterdir() if p.is_dir())
    strays = sorted(p.name for p in (reading / "score").iterdir() if not p.is_dir())
    if found != species or strays:
        stop(2, f"score/ holds {len(found)} species of the first run's {len(species)}; missing "
                f"{sorted(set(species) - set(found))}, more {sorted(set(found) - set(species))}, files {strays}")
    for slug in species:
        for name in SPECIES_FILES:
            if not (reading / "score" / slug / name).is_file():
                stop(2, f"score/{slug}/{name} is missing")
    files = {slug: {name: load(reading / "score" / slug / name) for name in SPECIES_FILES if name.endswith(".tsv")}
             for slug in species}
    counts = {name: load(reading / "count" / name) for name in COUNT_FILES}
    withdrawn = read_withdrawn(reading / "count" / "withdrawn.txt")
    arms_run = [a for a in ARMS if a not in withdrawn]
    by_species = {slug: [r for r in run if r["species_slug"] == slug] for slug in species}
    expect = {slug: expected(by_species[slug], slug, arms_run) for slug in species}

    # 3 datasets
    for slug in species:
        rows = files[slug]["datasets.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], r["learner"]) for r in rows]
        want = {(e["dataset"], e["split"], e["arm"], l) for e in expect[slug] for l in LEARNERS}
        if len(set(keys)) != len(keys) or set(keys) != want:
            more = sorted({k[0] for k in set(keys) - want})
            less = sorted({k[0] for k in want - set(keys)})
            stop(3, f"score/{slug}/datasets.tsv holds other rows than the first run's; more {more}, less {less}")
        for r in rows:
            mine = lists_of[r["dataset"]]
            if (r["species_slug"], r["species"], r["antibiotic"]) != (slug, mine["species"], mine["antibiotic"]):
                stop(3, f"{r['dataset']} is named otherwise than in lists.tsv")
            if r["genomes"] != mine["samples"]:
                stop(3, f"{r['dataset']} has {r['genomes']} genomes in datasets.tsv and {mine['samples']} in lists.tsv")
    for r in run:
        if r["floor_moves"] != str(int(score.floor_moves(r["samples"]))):
            stop(3, f"lists.tsv says the floor arm {'moves' if r['floor_moves'] == '1' else 'leaves'} {r['dataset']} "
                    f"at {r['samples']} genomes")

    # 4 study
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            if (r["split"] == "study") != (r["role"] == "study"):
                stop(4, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']} seed {r['seed']}: a {r['role']} row on "
                        f"the {r['split']} split")
            if r["role"] == "study" and int(r["seed"]) != STUDY_SEED:
                stop(4, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']}: a study run at seed {r['seed']}")

    # 5 twin
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            is_twin = r["arm"] == TWIN
            if is_twin != (r["role"] == "twin") or (is_twin and (r["split"] != "lineage"
                                                                 or int(r["seed"]) not in TWIN_SEEDS)):
                stop(5, f"score/{slug}/draws.tsv: {r['dataset']} {r['arm']} as {r['role']} on the {r['split']} split "
                        f"at seed {r['seed']}")
    named = sorted({r["arm"] for r in counts["conditions.tsv"]} - set(ARMS))
    if named:
        stop(5, f"conditions.tsv reads {', '.join(named)} as an arm")
    for name in ("count.tsv", "count-no-constant.tsv"):
        for r in counts[name]:
            if score.arm_of(r["pair"]) == TWIN and (r["split"] != "lineage" or r["set"] not in (ALL, MOVES)):
                stop(5, f"{name} counts the twin on the {r['split']} split over {r['set']}")

    # 6 missing
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            if r["scored"] != "1" and any(r[m] != "" for m in MEASURES):
                stop(6, f"score/{slug}/draws.tsv: {r['dataset']} {r['split']} {r['arm']} seed {r['seed']} "
                        f"{r['learner']}: a draw not scored has a figure")
        for r in files[slug]["datasets.tsv"]:
            whole_draws = r["scored_draws"] == r["draws"]
            if r["scored"] != str(int(whole_draws)) or (r["scored"] != "1" and any(r[k] != "" for k in MEASURES)):
                stop(6, f"score/{slug}/datasets.tsv: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: "
                        f"{r['scored_draws']} of {r['draws']} draws scored, read as scored {r['scored']}")

    # 7 constant
    for slug in species:
        for r in files[slug]["draws.tsv"]:
            marked = r["constant"] == "1"
            if marked != (r["why"] == REFUSED_WHY) or (marked and (r["split"] != "study" or r["arm"] not in ARMS
                                                                   or r["scored"] != "1")):
                stop(7, f"score/{slug}/draws.tsv: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: a constant "
                        f"model marked {r['constant']} for '{r['why']}'")
    for r in counts["count-no-constant.tsv"]:
        if r["constant"] not in ("", "0"):
            stop(7, f"count-no-constant.tsv counts {r['constant']} constant-model datasets for {r['pair']}")

    # 8 draws
    for slug in species:
        draws = files[slug]["draws.tsv"]
        keys = [(r["dataset"], r["split"], r["arm"], int(r["seed"]), r["role"], r["learner"]) for r in draws]
        want = {(e["dataset"], e["split"], e["arm"], int(e["seed"]), e["role"], l) for e in expect[slug]
                for l in LEARNERS}
        if len(set(keys)) != len(keys) or set(keys) != want:
            stop(8, f"score/{slug}/draws.tsv holds other draws than the first run's: {len(set(keys) - want)} more, "
                    f"{len(want - set(keys))} less, {len(keys) - len(set(keys))} twice")
        genomes = {r["dataset"]: int(r["samples"]) for r in by_species[slug]}
        again = score.dataset_rows([typed_draw(r) for r in draws], expect[slug], genomes)
        again_at = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in again}
        for r in files[slug]["datasets.tsv"]:
            g = again_at[(r["dataset"], r["split"], r["arm"], r["learner"])]
            for c in score.DATASET_COLUMNS:
                if not same_cell(r[c], fmt(g.get(c))):
                    stop(8, f"score/{slug}: {r['dataset']} {r['split']} {r['arm']} {r['learner']}: {c} is not the "
                            "mean of its draws")

    # 9 checks
    run_checks = []
    for slug in species:
        rows = files[slug]["checks.tsv"]
        bad = [r for r in rows if r["passed"] != "1"]
        if not rows or bad:
            stop(9, f"score/{slug}/checks.tsv: {len(bad)} of {len(rows)} checks did not pass")
        run_checks += rows
    agree = figure_checks(files, lists_of, arms_run)
    for name, passed, made in agree:
        if passed != made:
            stop(9, f"{name}: {passed} of {made}")

    # 10 screen
    for name in ("count.tsv", "count-no-constant.tsv"):
        for r in counts[name]:
            if score.arm_of(r["pair"]) in SCREEN_ARMS and r["learner"] not in SCREENED:
                stop(10, f"{name} counts {r['pair']} on {r['learner']}")

    # 11 floor
    floor_on = FLOOR_ARM not in withdrawn
    moved = {"lineage": sum(r["floor_moves"] == "1" for r in run),
             "study": sum(r["floor_moves"] == "1" and r["in_S"] == "1" for r in run)}
    every = {"lineage": len(run), "study": sum(r["in_S"] == "1" for r in run)}
    for name in ("count.tsv", "count-no-constant.tsv"):
        for r in counts[name]:
            if score.arm_of(r["pair"]) != FLOOR_ARM:
                continue
            want = moved[r["split"]] if floor_on else 0
            if r["set"] != MOVES or int(r["datasets"]) != want:
                stop(11, f"{name} counts the floor arm on the {r['split']} split over {r['set']}, {r['datasets']} "
                         f"datasets, not the {want} of 800 genomes or more")
    for r in counts["count-floor-all.tsv"]:
        if (score.arm_of(r["pair"]) != FLOOR_ARM or r["set"] != ALL or r["holm_p"] != ""
                or int(r["datasets"]) != every[r["split"]]):
            stop(11, f"count-floor-all.tsv holds {r['pair']} over {r['set']}, {r['datasets']} datasets, Holm p "
                     f"'{r['holm_p']}'")
    if floor_on and not counts["count-floor-all.tsv"]:
        stop(11, "count-floor-all.tsv is empty")

    # 12 withdrawn
    if any(w not in ARMS for w in withdrawn) or len(set(withdrawn)) != len(withdrawn):
        stop(12, f"withdrawn.txt names {withdrawn}, not arms of the block, each once")
    family = [r for r in counts["count.tsv"] if r["holm_p"] != ""]
    want = sorted((f"{a}-{DEFAULT}", s, l) for a, s in score.FAMILY_PARTS for l in touched(a))
    got = sorted((r["pair"], r["set"], r["learner"]) for r in family)
    if got != want or any(r["split"] != "lineage" or r["measure"] != PRIMARY for r in family):
        stop(12, f"the Holm family holds {len(family)} rows, not the {score.FAMILY_SIZE} lineage pairs of score.py")
    for r in family:
        gone = score.arm_of(r["pair"]) in withdrawn
        if (r["withdrawn"] == "1") != gone or (gone and (r["datasets"] != "0" or num(r["sign_p"]) != 1.0)):
            stop(12, f"{r['pair']} {r['learner']}: withdrawn {r['withdrawn']}, {r['datasets']} datasets, sign p "
                     f"{r['sign_p']}")

    # 13 binary
    bins = []
    for slug in species:
        fields, restarts = {}, 0
        for line in (reading / "score" / slug / "binary.txt").read_text(encoding="utf-8").splitlines():
            key, _, value = line.partition(" ")
            if key == "restart":
                restarts += 1
            elif key in ("commit", "grm", "code"):
                if key in fields:
                    stop(13, f"score/{slug}/binary.txt names its {key} twice")
                fields[key] = value.strip()
        if not (re.fullmatch(r"[0-9a-f]{40}", fields.get("commit", ""))
                and re.fullmatch(r"[0-9a-f]{64}", fields.get("grm", ""))
                and re.fullmatch(r"[0-9a-f]{64}", fields.get("code", ""))):
            stop(13, f"score/{slug}/binary.txt does not name a commit, a binary and a code sha256")
        bins.append({"species": by_species[slug][0]["species"], "species_slug": slug, **fields,
                     "restarts": restarts})

    # 14 conditions
    read = read_conditions(counts["count.tsv"])
    for r in counts["conditions.tsv"]:
        key = (r["arm"], r["reading"], r["condition"])
        if key in read and r["met"] != str(int(read[key])):
            stop(14, f"{r['arm']}, {r['reading']}, condition {r['condition']}: conditions.tsv says {r['met']}, read "
                     f"here {int(read[key])}")
    if len([r for r in counts["conditions.tsv"] if (r["arm"], r["reading"], r["condition"]) in read]) != len(read):
        stop(14, "conditions.tsv lacks a condition read here")

    # 15 count
    rows = [r for slug in species for r in files[slug]["datasets.tsv"]]
    try:
        counted = score.count(rows, withdrawn)
        beside = score.count(rows, withdrawn, keep_constant=False)
        floor_all = score.floor_all_count(rows, withdrawn)
        conds = score.conditions(counted)
    except SystemExit as e:
        stop(15, f"score.py refuses the count: {e}")
    want = {"count.tsv": fmt_rows(counted, score.COUNT_COLUMNS),
            "count-no-constant.tsv": fmt_rows(beside, score.COUNT_COLUMNS),
            "conditions.tsv": fmt_rows(conds, score.C3.CONDITION_COLUMNS),
            "count-floor-all.tsv": fmt_rows(floor_all, score.COUNT_COLUMNS)}
    for name in COUNT_FILES:
        if counts[name] != want[name]:
            stop(15, f"count/{name} is not what score.py --count writes from the datasets.tsv files")

    # 16 fits
    fits = {}
    for slug in species:
        path = reading / "score" / slug / "fits.tsv"
        rows = files[slug]["fits.tsv"]
        if header(path) != fit_times.COLUMNS:
            stop(16, f"score/{slug}/fits.tsv is not fit_times.py's table")
        arms = sorted({r["arm"] for r in rows})
        if arms != sorted(arms_run) or any(r["learner"] not in LEARNERS for r in rows):
            stop(16, f"score/{slug}/fits.tsv holds the arms {arms}, not the arms run, {sorted(arms_run)}")
        for r in rows:
            if num(r["seconds"]) is None or num(r["cpu_seconds"]) is None or num(r["peak_kb"]) is None:
                stop(16, f"score/{slug}/fits.tsv: a fit with no time")
        fits[slug] = rows

    return {"run": run, "species": species, "by_species": by_species, "files": files, "counts": counts,
            "withdrawn": withdrawn, "arms_run": arms_run, "run_checks": run_checks, "agree": agree,
            "binaries": bins, "fits": fits}


def figure_checks(files, lists_of, arms_run):
    """[(check, passed, made)]: a screen arm's scm and cart, and the floor arm below 800 genomes, equal the
    default's in datasets.tsv, every column but the arm's name. A study run refused for no candidate is scored as
    the constant model on every learner (ruling 16), and score.py makes neither check of it: nor does this."""
    screen = [a for a in SCREEN_ARMS if a in arms_run]
    floor = FLOOR_ARM in arms_run
    s_pass = s_made = f_pass = f_made = 0
    for slug in sorted(files):
        cell = {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in files[slug]["datasets.tsv"]}

        def same(a, b):
            return a is not None and b is not None and all(a[c] == b[c] for c in a if c != "arm")

        def refused(arm, d, split, learner):
            mine = cell.get((d, split, arm, learner))
            return mine is not None and int(mine["constant_draws"] or 0) > 0

        for (d, split, arm, learner), base in sorted(cell.items()):
            if arm != DEFAULT:
                continue
            for a in screen:
                if learner in SCREENED or refused(a, d, split, learner):
                    continue
                s_made += 1
                s_pass += same(cell.get((d, split, a, learner)), base)
            if floor and lists_of[d]["floor_moves"] == "0" and not refused(FLOOR_ARM, d, split, learner):
                f_made += 1
                f_pass += same(cell.get((d, split, FLOOR_ARM, learner)), base)
    return [("a screen arm's scm and cart figures equal the default's", s_pass, s_made),
            ("below 800 genomes the floor arm's figures equal the default's", f_pass, f_made)]


def meets(row):
    """(condition 1, condition 2) of one lineage primary row of count.tsv, as score.py reads them."""
    lo, holm = num(row["interval_lo"]), num(row["holm_p"])
    c1 = lo is not None and lo > 0
    c2 = holm is not None and holm < score.HOLM_ALPHA and int(row["better"]) > int(row["worse"])
    return c1, c2


def read_conditions(count_rows):
    """{(arm, reading, condition): met}: the deciding conditions 1 to 4, read here off count.tsv's rows."""
    at = {(r["split"], r["pair"], r["set"], r["learner"], r["measure"]): r for r in count_rows}

    def row(arm, learner, set_, measure=PRIMARY):
        key = ("lineage", f"{arm}-{DEFAULT}", set_, learner, measure)
        if key not in at:
            stop(14, f"count.tsv has no lineage row for {arm} with {learner} on {measure} over {set_}")
        return at[key]

    def low(r):
        hi = num(r["interval_hi"])
        return hi is not None and hi < -MARGIN

    out = {}
    for arm, set_ in [(a, ALL) for a in SCREEN_ARMS] + [(FLOOR_ARM, MOVES)]:
        learners = touched(arm)
        one = {l for l in learners if meets(row(arm, l, set_))[0]}
        two = {l for l in learners if meets(row(arm, l, set_))[1]}
        twin_most = max(sum(meets(row(TWIN, l, set_))[0] for l in LEARNERS),
                        sum(meets(row(TWIN, l, set_))[1] for l in LEARNERS))
        recall = not any(low(row(arm, l, set_, "recall")) for l in learners)
        if arm in SCREEN_ARMS:
            out.update({(arm, SCREEN_READING, "1"): len(one) == len(learners),
                        (arm, SCREEN_READING, "2"): len(two) == len(learners),
                        (arm, SCREEN_READING, "3"): twin_most < score.TOUCHED_AT_LEAST,
                        (arm, SCREEN_READING, "4"): recall})
        else:
            out.update({(arm, FLOOR_READING, "1 and 2"): len(one & two) >= score.TOUCHED_AT_LEAST,
                        (arm, FLOOR_READING, "3"): twin_most < score.TOUCHED_AT_LEAST,
                        (arm, FLOOR_READING, "4"): recall})
    return out


# ---------- the page's tables ----------

def page_conditions(conds):
    """conditions.tsv's rows, each marked: the deciding readings 1, the stricter reading 0."""
    return [{"arm": r["arm"], "reading": r["reading"], "decides": int(not r["reading"].startswith("stricter")),
             "condition": r["condition"], "reads": r["reads"], "met": r["met"], "learners": r["learners"]}
            for r in conds]


CONDITION_COLUMNS = ["arm", "reading", "decides", "condition", "reads", "met", "learners"]


def page_count(counts):
    """count.tsv's rows, the count that decides, each marked: in the Holm family or not, and conditions 1 and 2
    read on each family row."""
    out = []
    for r in counts["count.tsv"]:
        row = {"decides": 1, **r, "family": int(r["holm_p"] != ""), "meets_1": "", "meets_2": ""}
        if row["family"]:
            c1, c2 = meets(r)
            row["meets_1"], row["meets_2"] = int(c1), int(c2)
        out.append(row)
    return out


COUNT_COLUMNS = ["decides"] + score.COUNT_COLUMNS + ["family", "meets_1", "meets_2"]


def page_beside(counts):
    """count-no-constant.tsv's and count-floor-all.tsv's rows, each beside the count and deciding nothing."""
    out = []
    for beside, name in (("no-constant", "count-no-constant.tsv"), ("floor-all", "count-floor-all.tsv")):
        out += [{"beside": beside, "decides": 0, **r} for r in counts[name]]
    return out


BESIDE_COLUMNS = ["beside", "decides"] + score.COUNT_COLUMNS


def cells(rows):
    return {(r["dataset"], r["split"], r["arm"], r["learner"]): r for r in rows}


def diffs(rows, split, arm, learner, datasets, measure=PRIMARY):
    """[(dataset, difference)] of one split, pair and learner over `datasets`, the arm's value less the default's at
    12 decimals as score.py takes them; the datasets not scored; the datasets with a constant-model draw."""
    cell = cells(rows)
    out, unscored, constant = [], 0, 0
    for d in sorted(datasets):
        a, b = cell.get((d, split, arm, learner)), cell.get((d, split, DEFAULT, learner))
        if a is None:
            continue
        va, vb = num(a[measure]), (num(b[measure]) if b else None)
        if b is None or a["scored"] != "1" or b["scored"] != "1" or va is None or vb is None:
            unscored += 1
            continue
        out.append((d, r12(va - vb)))
        constant += int(a["constant_draws"] or 0) > 0
    return out, unscored, constant


def tally(values, measure=PRIMARY):
    up = -1 if measure in score.C3.LOWER_IS_BETTER else 1
    better = sum(1 for x in values if up * x > MARGIN)
    worse = sum(1 for x in values if up * x < -MARGIN)
    return better, worse, len(values) - better - worse


def median(values):
    return r12(statistics.median(values)) if values else None


def set_datasets(rows, set_):
    """The species' first-run datasets of a set: every one, or those of 800 genomes or more."""
    if set_ == ALL:
        return [r["dataset"] for r in rows]
    return [r["dataset"] for r in rows if r["floor_moves"] == "1"]


def page_species(found):
    """Per species, split, pair, set and learner, the primary measure's count over the species' datasets of the
    set; it decides nothing."""
    out = []
    for slug in found["species"]:
        rows = found["files"][slug]["datasets.tsv"]
        mine = found["by_species"][slug]
        for split, parts in PARTS.items():
            for arm, set_ in parts:
                if arm in found["withdrawn"]:
                    continue
                pool = set_datasets([r for r in mine if split == "lineage" or r["in_S"] == "1"], set_)
                for learner in touched(arm):
                    got, unscored, constant = diffs(rows, split, arm, learner, pool)
                    values = [x for _, x in got]
                    b, w, s = tally(values)
                    out.append({"species": mine[0]["species"], "species_slug": slug, "split": split,
                                "pair": f"{arm}-{DEFAULT}", "set": set_, "learner": learner, "decides": 0,
                                "datasets": len(pool), "entered": len(values), "not_scored": unscored,
                                "constant": constant, "better": b, "worse": w, "same": s, "median": median(values)})
    return out


SPECIES_COLUMNS = ["species", "species_slug", "split", "pair", "set", "learner", "decides", "datasets", "entered",
                   "not_scored", "constant", "better", "worse", "same", "median"]


def page_datasets(found):
    """Per dataset, split, pair and learner the pair is counted on: the default's value and the arm's or the
    twin's on the primary measure and recall, the differences, the draws scored, the genome count, whether the floor
    arm moves the dataset, whether this cell enters the deciding count, and a study run refused."""
    out = []
    for slug in found["species"]:
        cell = cells(found["files"][slug]["datasets.tsv"])
        for lrow in found["by_species"][slug]:
            d = lrow["dataset"]
            for split in ("lineage", "study"):
                if (d, split, DEFAULT, LEARNERS[0]) not in cell:
                    continue
                for arm in [a for a in ARMS if a not in found["withdrawn"]] + [TWIN]:
                    if (d, split, arm, LEARNERS[0]) not in cell:
                        continue
                    for learner in touched(arm):
                        base, mine = cell[(d, split, DEFAULT, learner)], cell[(d, split, arm, learner)]
                        row = {"dataset": d, "species": lrow["species"], "antibiotic": lrow["antibiotic"],
                               "genomes": lrow["samples"], "floor_moves": lrow["floor_moves"], "split": split,
                               "pair": f"{arm}-{DEFAULT}", "learner": learner,
                               "counted": int(arm != FLOOR_ARM or lrow["floor_moves"] == "1"),
                               "draws": mine["draws"], "default_draws_scored": base["scored_draws"],
                               "draws_scored": mine["scored_draws"],
                               "refused": int(int(mine["constant_draws"] or 0) > 0)}
                        for m, names in ((PRIMARY, ("default_value", "value", "difference")),
                                         ("recall", ("default_recall", "recall_value", "recall_difference"))):
                            va = num(mine[m]) if mine["scored"] == "1" else None
                            vb = num(base[m]) if base["scored"] == "1" else None
                            row[names[0]], row[names[1]] = vb, va
                            row[names[2]] = r12(va - vb) if va is not None and vb is not None else None
                        out.append(row)
    return out


DATASET_COLUMNS = ["dataset", "species", "antibiotic", "genomes", "floor_moves", "split", "pair", "learner", "counted",
                   "draws", "default_draws_scored", "draws_scored", "refused", "default_value", "value", "difference",
                   "default_recall", "recall_value", "recall_difference"]


def page_checks(found):
    """score.py's checks over every species, as checks passed of checks made per check, and the page's own
    agreement of the figures."""
    by = {}
    for r in found["run_checks"]:
        n = by.setdefault(r["check"], [0, 0])
        n[0] += r["passed"] == "1"
        n[1] += 1
    out = [{"source": "score.py checks.tsv", "check": c, "passed": p, "made": m} for c, (p, m) in sorted(by.items())]
    out += [{"source": "datasets.tsv", "check": c, "passed": p, "made": m} for c, p, m in found["agree"]]
    return out


CHECK_COLUMNS = ["source", "check", "passed", "made"]


def page_fit_times(found):
    """Per arm and learner, from fit_times.py's tables: the fits, and the sums and medians of their wall and CPU
    seconds and of the process's peak after each, in KB."""
    groups = {}
    for slug in found["species"]:
        for r in found["fits"][slug]:
            groups.setdefault((r["arm"], r["learner"]), []).append(r)
    out = []
    for arm in found["arms_run"]:
        for learner in LEARNERS:
            rs = groups.get((arm, learner), [])
            if not rs:
                continue
            sec = [float(r["seconds"]) for r in rs]
            cpu = [float(r["cpu_seconds"]) for r in rs]
            peak = [float(r["peak_kb"]) for r in rs]
            out.append({"arm": arm, "learner": learner, "fits": len(rs), "seconds_sum": r12(sum(sec)),
                        "seconds_median": median(sec), "cpu_seconds_sum": r12(sum(cpu)),
                        "cpu_seconds_median": median(cpu), "peak_kb_median": median(peak)})
    return out


FIT_COLUMNS = ["arm", "learner", "fits", "seconds_sum", "seconds_median", "cpu_seconds_sum", "cpu_seconds_median",
               "peak_kb_median"]
BINARY_COLUMNS = ["species", "species_slug", "commit", "grm", "code", "restarts"]


def binary_codes(bins):
    """The code sha256 values the species ran on."""
    return sorted({b["code"] for b in bins})


# ---------- writing ----------

README = """# The screen and floor reading: the record

This directory holds the tables behind [docs/screen-arms-result.md](../../docs/screen-arms-result.md): plan block C4, the Cochran-Mantel-Haenszel screen over lineage strata at budgets 4,000 and 16,000 (`cmh-4000`, `cmh-16000`) and the candidate floor fixed at 3 (`floor-3`), each against the release's chi-square screen at 4,000 with the floor max(3, n/200), beside a twin, read once on the first run's {datasets} datasets in {species} species, {study} of them with a study split. Every row is an aggregate over a dataset, a species, an arm or a learner. No file holds a genome, a sequence or a feature.

`scripts/c4-screen/score.py --count` wrote the four tables under `count/`, and `withdrawn.txt` beside them names the arms it took as withdrawn, one a line. `checks.txt` is the printout of `scripts/c4-screen/run_checks.py`, the checks of the run's copies made before any score was read. `scripts/c4-screen/page_tables.py write` checked the reading against the {gates} gates of `gates.txt`, and wrote every other file here once they passed. `reading.txt` lists the reading's commands, each with the sha256 of what it wrote. The rules of the reading were fixed at 05:52Z on 23 Nov 2026, before any C4 score was read (`docs/assumption-plan-corrections.md`, entry "C4 design").

| File | Holds |
|---|---|
| `checks.txt` | the checks of the run's copies before the reading, counts only |
| `gates.txt` | the gates the table script passed before it wrote a table |
| `reading.txt` | the reading's commands and the sha256 of each one's output |
| `count/count.tsv`, `count/conditions.tsv` | the count that decides, and the plan's conditions per arm as `score.py` read them |
| `count/count-no-constant.tsv`, `count/count-floor-all.tsv` | the counts beside it, each deciding nothing |
| `page-conditions.tsv` | the conditions per arm, `decides` 1 for the reading that decides and 0 for the stricter reading beside it; condition 6 is the maintainer's |
| `page-count.tsv` | every row of the count that decides, with the 21 rows of the Holm family marked and conditions 1 and 2 read on each |
| `page-beside.tsv` | every row of the two counts beside it, each marked as deciding nothing |
| `page-species.tsv` | per species, named in full, split, pair, set and learner: better, worse and the same on the primary measure; it decides nothing |
| `page-datasets.tsv` | per dataset, species and antibiotic in full, split, pair and learner: the default's value and the arm's on balanced accuracy and recall, the differences, the draws scored, the genome count, whether the floor arm moves the dataset and whether the cell is counted, and a study run refused |
| `page-checks.tsv` | `score.py`'s checks of every result as checks passed of checks made, and the table script's own check that a screen arm's `scm` and `cart`, and the floor arm below 800 genomes, equal the default |
| `page-fit-times.tsv` | per arm and learner, the outer fits' wall and CPU seconds and peak memory, as sums and medians |
| `page-binaries.tsv` | each species' commit, binary and code sha256, and its restarts |

{binaries}

A dataset's figure is the mean of its draws: seeds 42, 43 and 44 on the lineage split, the twin's 45, 46 and 47, and one study run at the release's seeds. A screen arm is counted on `rf`, `gbt` and `lasso`, the learners its screen reaches; its `scm` and `cart` read no screen and are checked equal to the default's, not counted. The floor arm is counted on all five learners over the datasets of 800 genomes or more, the ones its floor can move. A split by lineage or study is a retrospective split inside one public database, never external validation.

`SHA256SUMS` gives the sha256 of every other file here. `scripts/c4-screen/page_tables.py print` prints the page's tables from these files.
"""


def binaries_sentence(bins):
    codes = binary_codes(bins)
    commits = sorted({b["commit"] for b in bins})
    if len(codes) == 1:
        return (f"All {len(bins)} species ran one binary by its code sha256, `{codes[0]}`, built from "
                f"{len(commits)} commit{'s' if len(commits) != 1 else ''}: "
                + ", ".join(f"`{c[:7]}`" for c in commits) + ".")
    return (f"The {len(bins)} species ran {len(codes)} binaries by their code sha256; `page-binaries.tsv` names each "
            "binary with its species.")


def write_tsv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(fmt(r.get(c)) for c in columns) + "\n")


def write(reading, l_path, checks_path, out):
    reading, out = Path(reading), Path(out)
    lrows = load(l_path)
    if out.exists() and any(out.iterdir()):
        raise Stop(f"{out} is not empty")
    found = check_reading(reading, lrows, checks_path)
    (out / "count").mkdir(parents=True, exist_ok=True)
    for name in COUNT_FILES + ("withdrawn.txt",):
        (out / "count" / name).write_bytes((reading / "count" / name).read_bytes())
    (out / "reading.txt").write_bytes((reading / "reading.txt").read_bytes())
    (out / "checks.txt").write_bytes(Path(checks_path).read_bytes())
    write_tsv(out / "page-conditions.tsv", page_conditions(found["counts"]["conditions.tsv"]), CONDITION_COLUMNS)
    write_tsv(out / "page-count.tsv", page_count(found["counts"]), COUNT_COLUMNS)
    write_tsv(out / "page-beside.tsv", page_beside(found["counts"]), BESIDE_COLUMNS)
    write_tsv(out / "page-species.tsv", page_species(found), SPECIES_COLUMNS)
    write_tsv(out / "page-datasets.tsv", page_datasets(found), DATASET_COLUMNS)
    write_tsv(out / "page-checks.tsv", page_checks(found), CHECK_COLUMNS)
    write_tsv(out / "page-fit-times.tsv", page_fit_times(found), FIT_COLUMNS)
    write_tsv(out / "page-binaries.tsv", found["binaries"], BINARY_COLUMNS)
    with open(out / "gates.txt", "w", encoding="utf-8") as handle:
        for n, name, words in GATES:
            handle.write(f"{n} {name}: {words}\n")
    run = found["run"]
    (out / "README.md").write_text(README.format(
        datasets=len(run), species=len(found["species"]), study=sum(1 for r in run if r["in_S"] == "1"),
        gates=len(GATES), binaries=binaries_sentence(found["binaries"])), encoding="utf-8")
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
    return "< 0.001" if float(x) < 0.001 else d3(x)


def count_words(x):
    return f"{int(x):,}"


def of(part, whole_):
    return f"{count_words(part)} of {count_words(whole_)}"


def interval(lo, hi):
    return "none" if lo in ("", None) else f"{d3(lo)} to {d3(hi)}"


def yes(x):
    return {"1": "met", "0": "not met", "": "the maintainer's"}[str(x)]


def table(head, rows):
    lines = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    lines += ["| " + " | ".join(str(c) for c in row) + " |" for row in rows]
    return lines


def pair_words(pair):
    arm = score.arm_of(pair)
    return "twin less default" if arm == TWIN else f"`{arm}` less default"


def print_conditions(rows):
    return table(["Arm", "Reading", "Decides", "Condition", "Reads", "Met", "Learners"],
                 [[f"`{r['arm']}`", r["reading"], "yes" if r["decides"] == "1" else "no, beside", r["condition"],
                   r["reads"], yes(r["met"]), r["learners"] or "none"] for r in rows])


def print_family(rows):
    out = []
    for r in rows:
        if r["family"] != "1":
            continue
        out.append([pair_words(r["pair"]), SET_WORDS[r["set"]], f"`{r['learner']}`",
                    "yes" if r["withdrawn"] == "1" else "no", of(r["better"], r["entered"]), r["worse"], r["same"],
                    d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"]),
                    pval(r["holm_p"]), yes(r["meets_1"]), yes(r["meets_2"])])
    return table(["Pair", "Datasets", "Learner", "Withdrawn", "Better", "Worse", "Same", "Median", "Interval", "Sign p",
                  "Holm p", "Condition 1", "Condition 2"], out)


def print_count(rows, measures=(PRIMARY, "recall")):
    out = []
    for r in rows:
        if r["measure"] not in measures:
            continue
        out.append([r["split"], pair_words(r["pair"]), SET_WORDS[r["set"]], f"`{r['learner']}`",
                    MEASURE_WORDS[r["measure"]], r["orientation"], of(r["entered"], r["datasets"]), r["better"],
                    r["worse"], r["same"], r["constant"], d3(r["median_difference"]),
                    interval(r["interval_lo"], r["interval_hi"])])
    return table(["Split", "Pair", "Datasets", "Learner", "Measure", "Better is", "Entered", "Better", "Worse", "Same",
                  "Constant model", "Median", "Interval"], out)


def print_measures(rows):
    """Per split, pair, set and learner: the median difference of each measure beside balanced accuracy, recall
    beside log loss."""
    at = {(r["split"], r["pair"], r["set"], r["learner"], r["measure"]): r for r in rows}
    keys = []
    for r in rows:
        k = (r["split"], r["pair"], r["set"], r["learner"])
        if k not in keys:
            keys.append(k)
    out = []
    for k in keys:
        cells_ = [k[0], pair_words(k[1]), SET_WORDS[k[2]], f"`{k[3]}`"]
        for m in MEASURES:
            r = at[k + (m,)]
            cells_.append(f"{d3(r['median_difference'])} ({r['better']} / {r['worse']})")
        out.append(cells_)
    return table(["Split", "Pair", "Datasets", "Learner"] + [MEASURE_WORDS[m].capitalize() for m in MEASURES], out)


def print_beside(rows):
    out = []
    for r in rows:
        if r["measure"] not in (PRIMARY, "recall"):
            continue
        out.append([r["beside"], r["split"], pair_words(r["pair"]), SET_WORDS[r["set"]], f"`{r['learner']}`",
                    MEASURE_WORDS[r["measure"]], of(r["entered"], r["datasets"]), r["better"], r["worse"], r["same"],
                    d3(r["median_difference"]), interval(r["interval_lo"], r["interval_hi"]), pval(r["sign_p"])])
    return table(["Beside", "Split", "Pair", "Datasets", "Learner", "Measure", "Entered", "Better", "Worse", "Same",
                  "Median", "Interval", "Sign p"], out)


def print_checks(rows):
    return table(["Source", "Check", "Passed"], [[r["source"], r["check"], of(r["passed"], r["made"])] for r in rows])


def print_species(rows):
    """Per species, split, pair and set: better, worse and the same summed over the pair's learners, as dataset and
    learner pairs of the denominator."""
    groups = {}
    for r in rows:
        groups.setdefault((r["species"], r["split"], r["pair"], r["set"]), []).append(r)
    out = []
    for (sp, split, pair, set_), rs in groups.items():
        pairs = sum(int(r["datasets"]) for r in rs)
        b, w, s = (sum(int(r[k]) for r in rs) for k in ("better", "worse", "same"))
        out.append([f"*{sp}*", split, pair_words(pair), SET_WORDS[set_], len(rs),
                    f"{b} / {w} / {s} of {pairs}", sum(int(r["not_scored"]) for r in rs)])
    return table(["Species", "Split", "Pair", "Datasets", "Learners", "Better / worse / same, dataset and learner pairs",
                  "Not scored"], out)


def print_datasets(rows):
    groups = {}
    for r in rows:
        groups.setdefault((r["species"], r["antibiotic"], r["genomes"], r["floor_moves"], r["split"]), {}).setdefault(
            r["pair"], []).append(r)
    out = []
    for (sp, ab, genomes, moves, split), pairs in groups.items():
        cells_ = [f"*{sp}*", ab, count_words(genomes), "yes" if moves == "1" else "no", split]
        for arm in ARMS + (TWIN,):
            rs = pairs.get(f"{arm}-{DEFAULT}")
            cells_.append(" / ".join(d3(r["difference"]) for r in rs) if rs else "")
        refused = [score.arm_of(p) for p, rs in pairs.items() if any(r["refused"] == "1" for r in rs)]
        cells_.append(", ".join(f"`{a}`" for a in refused) or "none")
        out.append(cells_)
    return table(["Species", "Antibiotic", "Genomes", "Floor arm moves it", "Split",
                  "`cmh-4000`: `rf` / `gbt` / `lasso`", "`cmh-16000`: `rf` / `gbt` / `lasso`",
                  "`floor-3`: five learners", "Twin: five learners", "Study run refused"], out)


def print_fit_times(rows):
    return table(["Arm", "Learner", "Fits", "Seconds, sum", "Seconds, median", "CPU seconds, sum", "CPU seconds, median",
                  "Peak KB, median"],
                 [[f"`{r['arm']}`", f"`{r['learner']}`", count_words(r["fits"]), d3(r["seconds_sum"]),
                   d3(r["seconds_median"]), d3(r["cpu_seconds_sum"]), d3(r["cpu_seconds_median"]),
                   count_words(round(float(r["peak_kb_median"])))] for r in rows])


def print_binaries(rows):
    codes = binary_codes(rows)
    if len(codes) < 2:
        return [f"All {len(rows)} species ran one binary by its code sha256, `{codes[0][:16]}`."]
    return table(["Species", "Commit", "Code sha256", "Restarts"],
                 [[f"*{r['species']}*", f"`{r['commit'][:7]}`", f"`{r['code'][:16]}`", r["restarts"]] for r in rows])


def tables(folder):
    folder = Path(folder)
    count = load(folder / "page-count.tsv")
    return {
        "conditions": print_conditions(load(folder / "page-conditions.tsv")),
        "family": print_family(count),
        "count": print_count(count),
        "measures": print_measures(count),
        "beside": print_beside(load(folder / "page-beside.tsv")),
        "checks": print_checks(load(folder / "page-checks.tsv")),
        "fit-times": print_fit_times(load(folder / "page-fit-times.tsv")),
        "binaries": print_binaries(load(folder / "page-binaries.tsv")),
        "species": print_species(load(folder / "page-species.tsv")),
        "datasets": print_datasets(load(folder / "page-datasets.tsv")),
    }


def main(argv):
    if len(argv) == 5 and argv[0] == "write":
        try:
            files = write(argv[1], argv[2], argv[3], argv[4])
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
