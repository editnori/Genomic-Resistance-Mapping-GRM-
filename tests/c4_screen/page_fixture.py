"""SYNTHETIC readings of plan block C4 for the page's table script (scripts/c4-screen/page_tables.py) and the run's
check script (scripts/c4-screen/run_checks.py).

Every species, antibiotic, dataset, genome count, figure, fit time, binary and genome id here is made up in the
test. The species are named "Synthia alpha", "Synthia beta" and "Synthia gamma" and the antibiotics "synthetic
drug 1" and on, so that no figure here can be taken for a real one. The reading's tables are what score.py's own
functions write from SYNTHETIC draws: score.dataset_rows, and score.py --count run on the datasets.tsv files.

Each figure is a multiple of 1/64, so a mean of equal draws and a difference are exact, and every count row's
differences are one value: the count draws no species resample.
"""

import contextlib
import hashlib
import io
import json
import shutil
from pathlib import Path

from tests.c4_screen import common

score = common.load("score")
expect = common.load("expect")
fit_times = common.load("fit_times")
C3 = score.C3

LEARNERS = list(score.LEARNERS)
MEASURES = list(score.MEASURES)
SEEDS, TWIN_SEEDS = [42, 43, 44], [45, 46, 47]
LIST_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
                "in_M", "in_S", "first_run", "default_floor", "floor_moves", "arms"]
SPECIES = [("synthia_alpha", "Synthia alpha"), ("synthia_beta", "Synthia beta"), ("synthia_gamma", "Synthia gamma")]
# Per species: four first-run datasets, two of 800 genomes or more, two of them in S.
GENOMES = [900, 150, 1200, 300]
IN_S = [1, 1, 0, 0]
M_DATASET = "synthia_alpha__synthetic_drug_m"
# The dataset whose lineage split is unavailable in every draw, and the study run refused for no candidate.
UNSCORED = "synthia_beta__synthetic_drug_2"
REFUSED = ("synthia_gamma__synthetic_drug_1", "cmh-16000")
# Each arm's offset from the default, in 64ths, per learner: cmh-4000 better on the learners its screen reaches.
OFFSET = {"cmh-4000": {"rf": 2, "gbt": 2, "lasso": 2}}
COMMIT = "a" * 40
GRM = "b" * 64
CODE = "c" * 64
REFUSED_WHY = "constant model: refused for no candidate"


def list_rows():
    """lists.tsv's rows: the first run's twelve datasets and one dataset of M."""
    rows = []
    for k, (slug, name) in enumerate(SPECIES):
        for i, g in enumerate(GENOMES):
            rows.append({"dataset": f"{slug}__synthetic_drug_{i + 1}", "species_slug": slug, "species": name,
                         "antibiotic": f"synthetic drug {i + 1}", "samples": g, "resistant": g // 3,
                         "susceptible": g - g // 3, "lineages": 12, "in_M": 0, "in_S": IN_S[i], "first_run": 1,
                         "default_floor": max(3, g // 200), "floor_moves": int(max(3, g // 200) != 3),
                         "arms": "k31,cmh-4000,cmh-16000,floor-3,twin"})
    rows.append({"dataset": M_DATASET, "species_slug": "synthia_alpha", "species": "Synthia alpha",
                 "antibiotic": "synthetic drug m", "samples": 500, "resistant": 100, "susceptible": 400,
                 "lineages": 9, "in_M": 1, "in_S": 0, "first_run": 0, "default_floor": 3, "floor_moves": 0,
                 "arms": "k31"})
    return rows


def write_lists(path, rows=None):
    rows = list_rows() if rows is None else rows
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(LIST_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in LIST_COLUMNS) + "\n")


def run_of(rows, slug):
    return [r for r in rows if r["species_slug"] == slug and r["first_run"] == 1]


def expect_rows(rows, slug, arms):
    run = run_of(rows, slug)
    out = expect.rows(slug, [r["dataset"] for r in run], [r["dataset"] for r in run if r["in_S"]], SEEDS, TWIN_SEEDS,
                      list(arms))
    return [dict(zip(expect.COLUMNS, r)) for r in out]


def draw_rows(rows, slug, arms, cart_off=False):
    """draws.tsv's rows for one species: every figure the default's base for its dataset plus the arm's offset."""
    lrow = {r["dataset"]: r for r in rows}
    index = {r["dataset"]: i for i, r in enumerate(rows)}
    out = []
    for e in expect_rows(rows, slug, arms):
        d = e["dataset"]
        base = {"dataset": d, "species": lrow[d]["species"], "antibiotic": lrow[d]["antibiotic"], "arm": e["arm"],
                "split": e["split"], "seed": int(e["seed"]), "role": e["role"]}
        for learner in LEARNERS:
            if d == UNSCORED and e["split"] == "lineage":
                out.append(dict(base, learner=learner, scored=0, constant=0, why="split unavailable"))
                continue
            off = OFFSET.get(e["arm"], {}).get(learner, 0)
            if cart_off and e["arm"] == "cmh-4000" and learner == "cart":
                off = 1
            value = (24 + index[d] + off) / 64
            refused = (d, e["arm"]) == REFUSED and e["split"] == "study"
            row = dict(base, learner=learner, scored=1, constant=int(refused), why=REFUSED_WHY if refused else "",
                       n=60, resistant=20, susceptible=40, tp=15, fp=5, tn=35, fn=5)
            for m in ("recall", "specificity", "vme", "me"):
                row.update({m: value, f"{m}_lo": value - 1 / 64, f"{m}_hi": value + 1 / 64})
            row.update(balanced_accuracy=value, log_loss=value)
            out.append(row)
    return out


def fit_rows(rows, slug, arms):
    out = []
    for r in run_of(rows, slug):
        for arm in arms:
            for learner in LEARNERS:
                out.append({"arm": arm, "rule": "cmh" if arm.startswith("cmh") else "chi-square",
                            "budget": 16000 if arm == "cmh-16000" else 4000, "min_support": 3, "tag": "drawn",
                            "seed": 42, "dataset": r["dataset"], "split": "lineage", "fold": 0, "learner": learner,
                            "seconds": 1.5 if arm == "cmh-4000" else 0.5, "cpu_seconds": 3.0, "peak_kb": 2048})
    return out


def check_rows(slug, arms):
    names = ["exists", "release report", "scm and cart, which read no screen, make the default's calls in every fold",
             "below 800 genomes the floor arm's result is the default's, timing apart"]
    return [{"file": f"{arm}/{slug}__synthetic_drug_1.lineage-s42.pick.json", "check": c, "passed": 1, "detail": ""}
            for arm in arms for c in names]


def build(root, withdrawn=(), codes=None, cart_off=False):
    """A SYNTHETIC reading under root/reading, lists.tsv at root/lists.tsv and the run's checks.txt at
    root/checks.txt. codes maps a species slug to its code sha256 where it is not CODE."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    rows = list_rows()
    write_lists(root / "lists.tsv", rows)
    (root / "checks.txt").write_text("".join(f"PASS SYNTHETIC check {i}\n" for i in range(1, 12)) + "CHECKS 11 of 11\n",
                                     encoding="utf-8")
    reading = root / "reading"
    arms = [a for a in score.ARMS if a not in withdrawn]
    paths = []
    for slug, _ in SPECIES:
        out = reading / "score" / slug
        out.mkdir(parents=True)
        draws = draw_rows(rows, slug, arms, cart_off)
        score.write_tsv(out / "draws.tsv", draws, C3.DRAW_COLUMNS)
        genomes = {r["dataset"]: r["samples"] for r in run_of(rows, slug)}
        score.write_tsv(out / "datasets.tsv", score.dataset_rows(draws, expect_rows(rows, slug, arms), genomes),
                        score.DATASET_COLUMNS)
        score.write_tsv(out / "checks.tsv", check_rows(slug, arms), ["file", "check", "passed", "detail"])
        score.write_tsv(out / "fits.tsv", fit_rows(rows, slug, arms), fit_times.COLUMNS)
        code = (codes or {}).get(slug, CODE)
        (out / "binary.txt").write_text(f"commit {COMMIT}\ngrm {GRM}\ncode {code}\n", encoding="utf-8")
        paths.append(str(out / "datasets.tsv"))
    (reading / "count").mkdir()
    argv = ["--count", str(reading / "count")]
    for a in withdrawn:
        argv += ["--withdrawn", a]
    with contextlib.redirect_stdout(io.StringIO()):
        assert score.main(argv + paths) == 0
    (reading / "count" / "withdrawn.txt").write_text("".join(f"{a}\n" for a in withdrawn), encoding="utf-8")
    (reading / "reading.txt").write_text("SYNTHETIC reading: score.py --count over three made-up species\n",
                                         encoding="utf-8")
    return root


def copy(src, dst):
    shutil.copytree(src, dst)
    return Path(dst)


def read(path):
    return C3.read_tsv(path)


def rewrite(path, rows, columns=None):
    """Write rows (dicts of strings) back to a tsv in their own column order."""
    columns = columns or list(rows[0])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(columns) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in columns) + "\n")


def header(path):
    with open(path, encoding="utf-8") as f:
        return f.readline().rstrip("\n").split("\t")


# ---------- the run's copies, for run_checks.py ----------

GENOME_IDS = ["990001.11", "990001.12", "990002.21"]


def build_copies(root, codes=None, commits=None, refused=True):
    """SYNTHETIC copies of the run under root/copies/<species slug>/, as run_species.sh lays out its OUT, a
    SYNTHETIC release under root/release with made-up genome ids, and lists.tsv at root/lists.tsv."""
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    rows = list_rows()
    write_lists(root / "lists.tsv", rows)
    (root / "release" / "species").mkdir(parents=True)
    for slug, _ in SPECIES:
        doc = {"genomes": [{"id": g} for g in GENOME_IDS]}
        (root / "release" / "species" / f"{slug}.json").write_text(json.dumps(doc), encoding="utf-8")
        out = root / "copies" / slug
        exp = expect_rows(rows, slug, score.ARMS)
        picks = out / "results" / "picks"
        n_picks = 0
        for e in exp:
            (picks / e["arm"]).mkdir(parents=True, exist_ok=True)
            is_refused = refused and (e["dataset"], e["arm"]) == REFUSED and e["split"] == "study"
            if is_refused:
                (picks / C3.refused_name(e["file"])).write_text('{"reason": "no candidate"}\n', encoding="utf-8")
            else:
                (picks / e["file"]).write_text('{"format": "synthetic"}\n', encoding="utf-8")
                n_picks += 1
                if e["inputs"] != "-":
                    (picks / e["inputs"]).write_text('{"fit_times": []}\n', encoding="utf-8")
        for sub in ("check", "score", "logs"):
            (out / sub).mkdir()
        (out / "check" / "refit.tsv").write_text("dataset\tsame\nsynthetic\t1\n", encoding="utf-8")
        draws = draw_rows(rows, slug, score.ARMS)
        if not refused:
            for r in draws:
                if r.get("constant"):
                    r.update(constant=0, why="")
        score.write_tsv(out / "score" / "draws.tsv", draws, C3.DRAW_COLUMNS)
        score.write_tsv(out / "score" / "datasets.tsv", [], score.DATASET_COLUMNS)
        score.write_tsv(out / "score" / "checks.tsv", check_rows(slug, score.ARMS), ["file", "check", "passed",
                                                                                     "detail"])
        (out / "score" / "score.log").write_text("checks: SYNTHETIC\n", encoding="utf-8")
        score.write_tsv(out / "logs" / "fits.tsv", fit_rows(rows, slug, score.ARMS), fit_times.COLUMNS)
        (out / "logs" / "timing.tsv").write_text("step\tseconds\nsynthetic\t1\n", encoding="utf-8")
        (out / "logs" / "run.log").write_text("SYNTHETIC run log\n", encoding="utf-8")
        code = (codes or {}).get(slug, CODE)
        commit = (commits or {}).get(slug, COMMIT)
        (out / "binary.txt").write_text(f"commit {commit}\ngrm {GRM}\ncode {code}\n", encoding="utf-8")
        run = run_of(rows, slug)
        (out / "lists.txt").write_text(
            f"species\t{slug}\ndatasets\t{','.join(r['dataset'] for r in run)}\nseeds\t42,43,44\ntwin\t45,46,47\n"
            f"arms\tcmh-4000,cmh-16000,floor-3\n", encoding="utf-8")
        with open(out / "expect.tsv", "w", encoding="utf-8") as f:
            f.write("\t".join(expect.COLUMNS) + "\n")
            for e in exp:
                f.write("\t".join(str(e[c]) for c in expect.COLUMNS) + "\n")
        sums(out)
        done(out, slug, n_picks)
    return root


def listed(out):
    """The files run_species.sh's SHA256SUMS covers."""
    out = Path(out)
    names = [p.relative_to(out).as_posix() for top in ("results", "check", "score") for p in (out / top).rglob("*")
             if p.is_file() and p.suffix in (".json", ".tsv", ".txt")]
    return sorted(names + ["binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv", "logs/fits.tsv"])


def sums(out):
    out = Path(out)
    lines = [f"{hashlib.sha256((out / n).read_bytes()).hexdigest()}  {n}\n" for n in listed(out)]
    (out / "SHA256SUMS").write_text("".join(lines), encoding="utf-8")


def done(out, slug, n_picks, time="2026-11-23T13:00:00Z"):
    out = Path(out)
    lines = dict(line.split(" ", 1) for line in (out / "binary.txt").read_text(encoding="utf-8").splitlines())
    sha = hashlib.sha256((out / "SHA256SUMS").read_bytes()).hexdigest()
    (out / "DONE").write_text(f"DONE {slug} {time} binary {lines['grm']} code {lines['code']} results {n_picks} "
                              f"SHA256SUMS {sha}\n", encoding="utf-8")
