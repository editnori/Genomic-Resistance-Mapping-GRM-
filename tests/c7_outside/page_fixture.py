"""A SYNTHETIC C7 run for page_tables.py and run_checks.py: copies as the run's come back, a release's species lists,
the reading run_checks.py gather writes from them, and the count score.py writes from that.

Every species, dataset, antibiotic, genome id, call, fold count, time and sha256 here is made up in the test: no
genome, call or score comes from a real dataset or a real run.
"""

import hashlib
import io
import json
import shutil
import tarfile
from pathlib import Path

from tests.c7_outside import common

score = common.load("score")
checks = common.load("run_checks")
page = common.load("page_tables")
expect = common.load("expect")
O, C9 = score.O, score.C9

LIST_COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
                "in_M", "in_S", "first_run", "labels_depend_on", "arms", "splits"]
# species slug, species, [(antibiotic, genomes, in S, in M)]
SPECIES = [
    ("synthia_coli", "Synthia coli", [("drug a", 150, 1, 0), ("drug b", 120, 0, 0), ("drug m", 90, 0, 1)]),
    ("fictus_aureus", "Fictus aureus", [("drug a", 140, 1, 0), ("drug c", 110, 0, 0)]),
    ("madeup_pneumoniae", "Madeup pneumoniae", [("drug b", 100, 0, 0)]),
]
SEEDS, TWIN_SEEDS, ARM_SEEDS = [42, 43, 44], [45, 46, 47], [42, 43, 44]
# The study split is unavailable for this dataset: its default rows and its arms' study reports hold no score.
UNAVAILABLE = "fictus_aureus__drug_a"
# Calls per fold of five, the same in every draw of a learner, so each pair's differences are one value.
DEFAULT_CALLS = {"gbt": {"n": 30, "tp": 12, "fp": 4, "tn": 11, "fn": 3, "loss_sum": 12.0},
                 "lasso": {"n": 30, "tp": 11, "fp": 3, "tn": 12, "fn": 4, "loss_sum": 13.0}}
ARM_CALLS = {"xgboost": {"n": 30, "tp": 14, "fp": 2, "tn": 13, "fn": 1, "loss_sum": 9.0},
             "lightgbm": {"n": 30, "tp": 12, "fp": 4, "tn": 11, "fn": 3, "loss_sum": 12.0},
             "elasticnet": {"n": 30, "tp": 9, "fp": 5, "tn": 10, "fn": 6, "loss_sum": 16.0}}
# One booster fold that kept one tree, one constant fold, and one elastic net refit that reached its epoch limit.
ONE_TREE = ("xgboost", "synthia_coli__drug_a", "lineage", 43)
CONSTANT = ("lightgbm", "fictus_aureus__drug_c", "lineage", 44)
UNCONVERGED = ("elasticnet", "madeup_pneumoniae__drug_b", "lineage", 42)
LIBRARIES = {"python": "3.12.3", "numpy": "1.26.4", "scipy": "1.14.1", "sklearn": "1.7.2", "xgboost": "3.0.5",
             "lightgbm": "4.6.0"}
# Made-up sha256 values, hashed from words so that no run of hex letters reads as nucleotides.
GRM, CODE, EXPORT = (hashlib.sha256(f"SYNTHETIC {w}".encode()).hexdigest() for w in ("grm", "code", "export"))
COMMIT = hashlib.sha1(b"SYNTHETIC commit").hexdigest()
DONE_AT = "2026-11-23T10:00:00Z"
WITHDRAWN = {"elasticnet": ("SYNTHETIC: the elastic net's cost", "2026-11-23T09:00Z")}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def slug_of(name):
    return name.replace(" ", "_")


def lists_rows():
    out = []
    for slug, name, drugs in SPECIES:
        for drug, n, in_s, in_m in drugs:
            out.append({"dataset": f"{slug}__{slug_of(drug)}", "species_slug": slug, "species": name,
                        "antibiotic": drug, "samples": n, "resistant": n // 2, "susceptible": n - n // 2, "lineages": 9,
                        "in_M": in_m, "in_S": in_s, "first_run": 1 - in_m, "labels_depend_on": "",
                        "arms": "k31,twin,xgboost,lightgbm,elasticnet",
                        "splits": "lineage,study" if in_s and not in_m else "lineage"})
    return out


def genome_ids(slug):
    base = {"synthia_coli": 90001, "fictus_aureus": 90002, "madeup_pneumoniae": 90003}[slug]
    return [f"{base}.{i}" for i in range(1, 6)]


def write_tsv(path, rows, columns):
    with open(path, "w", newline="", encoding="utf-8") as handle:
        handle.write("\t".join(columns) + "\n")
        for r in rows:
            handle.write("\t".join(str(r[c]) for c in columns) + "\n")


def times(calls, k):
    return {key: v * k for key, v in calls.items()}


def default_rows(e, names):
    sp, ab = names[e["dataset"]]
    base = {"dataset": e["dataset"], "species": sp, "antibiotic": ab, "arm": e["arm"], "split": e["split"],
            "seed": int(e["seed"]), "role": e["role"]}
    out = []
    for learner in score.PAIRED:
        if e["dataset"] == UNAVAILABLE and e["split"] == "study":
            out.append(dict(base, learner=learner, scored=0, constant=0, unconverged=0, why="split unavailable"))
        else:
            out.append(dict(base, learner=learner, scored=1, constant=0, unconverged=0, why="",
                            **C9.measures(times(DEFAULT_CALLS[learner], 5))))
    return out


def report(e):
    arm, d, scheme, seed = e["arm"], e["dataset"], e["scheme"], int(e["seed"])
    doc = common.outside_report(d, scheme, seed, arm, calls=[dict(ARM_CALLS[arm]) for _ in range(5)])
    split = e["split"]
    if d == UNAVAILABLE and split == "study":
        doc.update(unavailable_reason="SYNTHETIC: one study", folds=[], tested=0, tested_resistant=0, input_rows=None)
        return doc
    if (arm, d, split, seed) == ONE_TREE:
        doc["folds"][1]["rounds_kept"] = 1
    if (arm, d, split, seed) == CONSTANT:
        doc["folds"][2]["constant"] = 1
    if (arm, d, split, seed) == UNCONVERGED:
        doc["folds"][3]["refit_converged"] = 0
    doc.update(page.fold_sums(doc))
    return doc


def species_copy(folder, slug, lrows, arms, elasticnet):
    """One species' copy as the driver brings it back: top files, logs/, results/outside/ and score.tar."""
    mine = [r for r in lrows if r["species_slug"] == slug and r["first_run"] == 1]
    ds = [r["dataset"] for r in mine]
    st = [r["dataset"] for r in mine if r["in_S"] == 1]
    names = {r["dataset"]: (r["species"], r["antibiotic"]) for r in lrows}
    rows = checks.expect_rows(expect.rows(slug, ds, st, SEEDS, TWIN_SEEDS, ARM_SEEDS, arms))
    copy = folder / slug
    (copy / "logs").mkdir(parents=True)
    (copy / "results" / "outside").mkdir(parents=True)
    (copy / "lists.txt").write_text(f"species\t{slug}\ndatasets\t{','.join(ds)}\nseeds\t42,43,44\ntwin\t45,46,47\n"
                                    f"arm_seeds\t42,43,44\narms\t{','.join(arms)}\n", encoding="utf-8")
    (copy / "libraries.txt").write_text("".join(f"{k} {v}\n" for k, v in LIBRARIES.items()), encoding="utf-8")
    write_tsv(copy / "expect.tsv", rows, expect.COLUMNS)
    (copy / "export-binary.txt").write_text(f"commit {COMMIT}\ngrm {EXPORT}\n", encoding="utf-8")
    docs, draws, timing, check_rows = {}, [], [], []
    for e in rows:
        check_rows.append({"file": e["file"], "check": "exists", "passed": 1, "detail": ""})
        if e["arm"] in score.ARMS:
            doc = report(e)
            docs[e["file"]] = doc
            path = copy / e["file"]
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(O.dump(doc), encoding="utf-8")
            timing.append({"fit": f"{e['arm']}-{e['dataset']}-{e['scheme']}-s{e['seed']}",
                           "wall_s": {"xgboost": 2.5, "lightgbm": 1.5, "elasticnet": 40.0}[e["arm"]],
                           "cpu_s": {"xgboost": 2.0, "lightgbm": 1.25, "elasticnet": 39.0}[e["arm"]], "peak_kb": 50000})
            check_rows.append({"file": e["file"], "check": "listed keys only", "passed": 1, "detail": ""})
        else:
            draws += default_rows(e, names)
    arm_rows = [e for e in rows if e["arm"] in score.ARMS]
    draws += score.score_rows(arm_rows, docs, names)
    write_tsv(copy / "logs" / "timing.tsv", timing, ["fit", "wall_s", "cpu_s", "peak_kb"])
    seeds = {(e["dataset"], e["seed"]) for e in arm_rows if e["split"] == "lineage"}
    at42 = sum(1 for _, s in seeds if s == "42")
    (copy / "logs" / "plans.log").write_text(
        f"plans: {3 * len(seeds) + at42} of {3 * len(seeds) + at42} checks pass, {len(seeds)} dataset seeds\n",
        encoding="utf-8")
    (copy / "logs" / "run.log").write_text("SYNTHETIC run log: fits end\n", encoding="utf-8")
    tables = {}
    for name, items, columns in (("checks.tsv", check_rows, ["file", "check", "passed", "detail"]),
                                 ("draws.tsv", draws, score.DRAW_COLUMNS),
                                 ("datasets.tsv", score.dataset_rows(draws), score.DATASET_COLUMNS)):
        buffer = io.StringIO()
        buffer.write("\t".join(columns) + "\n")
        for r in items:
            buffer.write("\t".join(C9.fmt(r.get(c)) for c in columns) + "\n")
        tables[f"score/{name}"] = buffer.getvalue().encode("utf-8")
    tables["score/score.log"] = b"SYNTHETIC score log\n"
    pack(copy / "score.tar", tables)
    listed = sorted([p.relative_to(copy).as_posix() for p in (copy / "results").rglob("*.json")]
                    + [n for n in tables if n.endswith(".tsv")]
                    + ["libraries.txt", "expect.tsv", "lists.txt", "export-binary.txt", "logs/timing.tsv"])
    sums = "".join(f"{sha(tables[n] if n in tables else (copy / n).read_bytes())}  {n}\n" for n in listed)
    (copy / "SHA256SUMS").write_text(sums, encoding="utf-8")
    (copy / "DONE").write_text(
        f"DONE {slug} {DONE_AT} grm {GRM} code {CODE} export-grm {EXPORT} elasticnet {elasticnet} results "
        f"{len(docs)} SHA256SUMS {sha(sums.encode('utf-8'))}\n", encoding="utf-8")
    (copy / "back.tar").write_bytes(b"")
    return copy


def pack(path, files):
    with tarfile.open(path, "w") as archive:
        for name, data in sorted(files.items()):
            info = tarfile.TarInfo(name)
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))


def resum(copy):
    """Write a copy's SHA256SUMS and DONE again after a test changed a file SHA256SUMS lists."""
    tables = checks.score_files(copy)
    lines = []
    for line in (copy / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        name = line[66:]
        data = tables[name] if name.startswith("score/") else (copy / name).read_bytes()
        lines.append(f"{sha(data)}  {name}\n")
    sums = "".join(lines)
    (copy / "SHA256SUMS").write_text(sums, encoding="utf-8")
    done = (copy / "DONE").read_text(encoding="utf-8").split()
    done[-1] = sha(sums.encode("utf-8"))
    (copy / "DONE").write_text(" ".join(done) + "\n", encoding="utf-8")


def build(folder, withdrawn=False):
    """folder/lists.tsv, folder/release/species/, folder/back/ (one copy per species), folder/reading/ (gathered and
    counted). With withdrawn, the elastic net is withdrawn and every copy is of the first pass alone."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    lrows = lists_rows()
    write_tsv(folder / "lists.tsv", lrows, LIST_COLUMNS)
    arms = ("xgboost", "lightgbm") if withdrawn else ("xgboost", "lightgbm", "elasticnet")
    for slug, _, _ in SPECIES:
        (folder / "release" / "species").mkdir(parents=True, exist_ok=True)
        (folder / "release" / "species" / f"{slug}.json").write_text(
            json.dumps({"genomes": [{"id": g} for g in genome_ids(slug)]}), encoding="utf-8")
        species_copy(folder / "back", slug, lrows, arms, "leave-out" if withdrawn else "run")
    w = WITHDRAWN if withdrawn else {}
    rep = checks.gather(folder / "back", folder / "lists.tsv", folder / "release", folder / "reading", w)
    if rep.failed:
        raise AssertionError("the SYNTHETIC copies fail run_checks.py:\n" + "\n".join(rep.lines))
    count(folder / "reading", sorted(w))
    return folder


def count(reading, withdrawn):
    """score.py --count's four files from the reading's datasets.tsv files, as part 2 runs it."""
    reading = Path(reading)
    paths = sorted((reading / "score").glob("*/datasets.tsv"))
    rows = [r for p in paths for r in C9.read_tsv(p)]
    out = reading / "count"
    out.mkdir()
    counted = score.count(rows, withdrawn)
    C9.write_tsv(out / "count.tsv", counted, score.COUNT_COLUMNS)
    C9.write_tsv(out / "count-no-constant.tsv", score.count(rows, withdrawn, keep_constant=False), score.COUNT_COLUMNS)
    C9.write_tsv(out / "count-no-unconverged.tsv", score.count(rows, withdrawn, keep_unconverged=False),
                 score.COUNT_COLUMNS)
    C9.write_tsv(out / "conditions.tsv", score.conditions(counted), score.CONDITION_COLUMNS)
    with open(reading / "reading.txt", "a", encoding="utf-8") as handle:
        handle.write("score.py --count: SYNTHETIC\n")
        for p in sorted(out.iterdir()):
            handle.write(f"{sha(p.read_bytes())}  count/{p.name}\n")


def copy_of(src, dst):
    shutil.copytree(src, dst)
    return Path(dst)
