"""SYNTHETIC releases for stage 2 of the matched accuracy protocol, and one run of stage 1 then stage 2 on each.

  python3 -I tests/matched_accuracy/made_up2.py <top> [--species N] [--seed S]

writes a made-up release under <top>/in, runs scripts/matched-accuracy/stage1.py on it into <top>/stage1, then
scripts/matched-accuracy/stage2.py into <top>/stage2. Each script's pins are set to the made-up files' by a route
the command line does not reach: the module's own constants, set in this process (a real run checks the real
ones). Without --seed it writes the designed release of 8 species worked by hand in
tests/matched_accuracy/test_stage2_made_up.py; with --species 5, its first 5 species.

Everything here is made up: the species, the antibiotics, the genome ids (gz-<species>-<n>, which is not the
shape of a real genome id), the BioProject accessions (an accession's shape, which stage 1's rule needs; made_up.py
says why), the lab records, the label table, the genome tables, the catalogue calls and every score.
"""

import hashlib
import importlib.util
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
STAGE1 = ROOT / "scripts" / "matched-accuracy" / "stage1.py"
STAGE2 = ROOT / "scripts" / "matched-accuracy" / "stage2.py"


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


made_up = _load("matched_accuracy_made_up_of_stage2", HERE / "made_up.py")
FIVE = made_up.FIVE
TAGS = ("ka", "la", "mu", "nu", "xi", "om", "rh", "ta", "ph", "ch", "ps", "et")
NAMES = ("kappa", "lambda", "mu", "nu", "xi", "omicron", "rho", "tau", "phi", "chi", "psi", "eta")
DRUGS = ("drug_a", "drug_b", "drug_c")
# Values no output of the second loader may carry: each of roc, metrics, conformal, no_call, called_metrics,
# lasso_fits and fold_balanced_accuracy, and the top-level catalogue's own values.
SENTINELS = {"roc": "roc-sentinel", "metrics": "metrics-sentinel", "conformal": "conformal-sentinel",
             "no_call": "no-call-sentinel", "called_metrics": "called-metrics-sentinel",
             "lasso_fits": "lasso-fits-sentinel", "catalog": "catalog-sentinel"}
BALANCED = 0.3141592653


def bioproject(species, n):
    """A made-up BioProject with an accession's shape, built from pieces: the letters ZZ, which no registry uses,
    then 7, the species' number and the project's number."""
    return "PRJ" + "ZZ" + f"7{species:02d}{n:02d}"


def row(call, score, mlst="ST-1", hc20="hc-1", catalog=None, intermediate=False):
    """One planned row: its lab call (True resistant), its five scores (one number gives all five), MLST, HC20
    (None: the 5 Nov table lacks the genome), catalogue call ("rule" takes the lab call) and an I record."""
    five = tuple(score) if isinstance(score, (tuple, list)) else (score,) * 5
    return {"resistant": call, "five": five, "mlst": mlst, "hc20": hc20,
            "catalog": call if catalog == "rule" else catalog, "intermediate": intermediate}


def R(score, **kw):
    return row(True, score, **dict({"catalog": "rule"}, **kw))


def S(score, **kw):
    return row(False, score, **dict({"catalog": "rule"}, **kw))


def group(fold, lineage, rows):
    return {"fold": fold, "lineage": lineage, "rows": rows}


def dataset(drug, groups, project=1, rule=True, null_folds=None):
    """A dataset: its groups, its BioProject number, whether the report holds a drug rule, and {learner: folds}
    whose fold_complexity is 0, a null fit."""
    return {"drug": drug, "groups": groups, "project": project, "rule": rule, "null_folds": null_folds or {}}


def designed():
    """The 8 species worked by hand in test_stage2_made_up.py. A group's value is given beside it."""
    return [
        # kappa: one dataset; group 2 holds the one exact tie of the release, R 0.6 against S 0.6, and its S row
        # at 0.3 carries an I record. Two rows of group 1 are of ST-2, so the MLST match splits it in two groups
        # of value 1.
        [dataset("drug_a", [
            group(0, 10, [R(0.9, mlst="ST-2"), R(0.7), S(0.8, mlst="ST-2"), S(0.2), S(0.1)]),  # 5/6, 6 pairs
            group(1, 11, [R(0.6), R(0.7), S(0.6), S(0.3, intermediate=True)]),          # 7/8, 4 pairs
        ])],
        # lambda: two datasets, the first with three MLST types and the second in BioProject 2.
        [dataset("drug_a", [
            group(0, 10, [R(0.9), R(0.8), R(0.7, mlst="ST-2"), S(0.5), S(0.4, mlst="ST-2"),
                          S(0.3, mlst="ST-3"), S(0.2, mlst="ST-3")]),                    # 1, 12 pairs
        ]), dataset("drug_b", [
            group(2, 12, [R(0.4), R(0.3), S(0.5), S(0.2), S(0.1)]),                      # 2/3, 6 pairs
            group(3, 13, [R(0.2), R(0.9), S(0.3), S(0.1)]),                              # 3/4, 4 pairs
        ], project=2)],
        # mu: one dataset of 6 pairs, under 10.
        [dataset("drug_a", [
            group(0, 10, [R(0.7), R(0.4), S(0.5), S(0.1)]),                              # 3/4, 4 pairs
            group(1, 11, [R(0.2), S(0.6), S(0.3)]),                                      # 0, 2 pairs
        ])],
        # nu: a report with no drug rule.
        [dataset("drug_a", [
            group(0, 10, [R(0.8), R(0.6), R(0.5), S(0.7), S(0.4), S(0.3), S(0.2)]),      # 5/6, 12 pairs
        ], rule=False)],
        # xi: gbt has a null fit in fold 4.
        [dataset("drug_a", [
            group(0, 10, [R(0.9), R(0.5), S(0.6), S(0.4), S(0.3)]),                      # 5/6, 6 pairs
            group(4, 14, [R(0.3), R(0.2), S(0.4), S(0.1)]),                              # 1/2, 4 pairs
        ], null_folds={"gbt": (4,)})],
        # omicron: one row lacks an MLST type; one genome the 5 Nov table lacks; the R row at 0.1 is the one row
        # of ST-2 in its group, so the MLST match leaves it out.
        [dataset("drug_a", [
            group(0, 10, [R(0.9), R(0.8), S(0.7, mlst=None), S(0.6), S(0.5, hc20=None)]),  # 1, 6 pairs
            group(1, 11, [R(0.1, mlst="ST-2"), R(0.5), S(0.4), S(0.3)]),                 # 1/2, 4 pairs
        ])],
        # rho: two datasets; one S row with an unevaluable catalogue call, one R row called susceptible.
        [dataset("drug_a", [
            group(0, 10, [R(0.6), R(0.4, catalog=False), S(0.5), S(0.3), S(0.2)]),       # 5/6, 6 pairs
            group(2, 12, [R(0.7), R(0.8), S(0.9, catalog=True), S(0.1)]),                # 1/2, 4 pairs
        ]), dataset("drug_b", [
            group(3, 13, [R(0.5), R(0.5), R(0.5), S(0.4), S(0.4), S(0.4), S(0.6, catalog=None)]),  # 3/4, 12 pairs
        ])],
        # tau: five scores that differ by learner, each row's mean as a number gives.
        [dataset("drug_a", [
            group(0, 10, [R((0.3, 0.3, 0.3, 0.3, 0.3)), R(0.2), S(0.8), S(0.6), S(0.4)]),  # 0, 6 pairs
            group(1, 11, [R((0.9, 0.8, 1.0, 0.9, 0.9)), R(0.95), S(0.1), S(0.2)]),       # 1, 4 pairs
        ])],
    ]


def randomised(seed, species=7):
    """A release of `species` species with groups, scores, ties, typing values, catalogue calls and null fits
    drawn from random.Random(seed)."""
    rng = random.Random(seed)
    out = []
    for _ in range(species):
        datasets = []
        for d in range(rng.randint(1, 3)):
            groups = []
            for g in range(rng.randint(1, 4)):
                rows = []
                for call in [True] * rng.randint(1, 4) + [False] * rng.randint(1, 4):
                    five = tuple(round(rng.random(), rng.choice((2, 6, 17))) for _ in FIVE)
                    if rows and rng.random() < 0.15:
                        five = rng.choice(rows)["five"]   # an exact tie with an earlier row of the group
                    rows.append(row(call, five, mlst=rng.choice(("ST-1", "ST-1", "ST-2", None)),
                                    hc20=rng.choice(("hc-1", "hc-1", "hc-2", None)),
                                    catalog=rng.choice((True, False, True, None)) if d != 2 else None,
                                    intermediate=rng.random() < 0.05))
                groups.append(group(g % 5, 10 + g, rows))
            nulls = {"gbt": (rng.randrange(5),)} if rng.random() < 0.4 else {}
            if rng.random() < 0.3:
                nulls["lasso"] = (rng.randrange(5),)
            datasets.append(dataset(DRUGS[d], groups, project=rng.randint(1, 2), rule=d != 2, null_folds=nulls))
        out.append(datasets)
    return out


def learner(name, scores, complexity):
    e = {"conformal": [SENTINELS["conformal"]], "fold_balanced_accuracy": [BALANCED] * 5,
         "fold_complexity": complexity, "learner": name, "metrics": {"note": SENTINELS["metrics"]},
         "no_call": {"resistant": 0, "susceptible": 0, "note": SENTINELS["no_call"]},
         "roc": [{"note": SENTINELS["roc"]}], "scores": scores, "seconds": 1.25}
    if name in ("gbt", "lasso") and 0 in complexity:
        e["called_metrics"] = {"note": SENTINELS["called_metrics"]}
    if name == "lasso":
        e["lasso_fits"] = [{"note": SENTINELS["lasso_fits"]}]
    return e


def build(plan, filler_score=0.123456789):
    """The release as made_up.write_release takes it. Each dataset fills every fold its groups leave empty with
    one S row in a lineage of its own, which stage 1 leaves out under the one-outcome rule; its scores are made-up
    numbers that no listed row carries."""
    reports, records, labels, genomes, hc20 = [], [], [], [], []
    for k, datasets in enumerate(plan):
        tag, name = TAGS[k], NAMES[k]
        species, slug = f"Synthia {name}", f"synthia_{name}"
        n = 0
        for ds in datasets:
            rows = [dict(r, fold=g["fold"], lineage=g["lineage"]) for g in ds["groups"] for r in g["rows"]]
            for f in sorted(set(range(5)) - {r["fold"] for r in rows}):
                rows.append(dict(row(False, filler_score + f / 1000, catalog="rule" if ds["rule"] else None),
                                 fold=f, lineage=90 + f))
            for r in rows:
                n += 1
                r["id"] = made_up.gid(tag, n)
                if not ds["rule"]:
                    r["catalog"] = None
            dataset_id = f"{slug}__{ds['drug']}"
            base = [{"id": r["id"], "resistant": r["resistant"], "fold": r["fold"], "lineage": r["lineage"],
                     "mlst": r["mlst"], "from_mic": False} for r in rows]
            doc = made_up.report(dataset_id, species, slug, ds["drug"], base, None,
                                 FIVE + (("catalog",) if ds["rule"] else ()))
            lineage = next(s for s in doc["splits"] if s["scheme"] == "lineage")
            entries = []
            for j, lname in enumerate(FIVE):
                complexity = [3.0 if lname != "rf" else 50.0] * 5
                for f in ds["null_folds"].get(lname, ()):
                    complexity[f] = 0.0
                entries.append(learner(lname, [r["five"][j] for r in rows], complexity))
            if ds["rule"]:
                entries.append({"conformal": [], "fold_balanced_accuracy": [], "fold_complexity": [],
                                "learner": "catalog", "metrics": {"note": SENTINELS["metrics"]}, "no_call": None,
                                "roc": [], "scores": [0.987654321 if r["catalog"] else 0.012345678 for r in rows],
                                "seconds": 0.5})
            lineage["learners"] = entries
            for s, r in zip(doc["samples"], rows):
                s["catalog_resistant"] = r["catalog"]
            doc["catalog"] = {"groups": [SENTINELS["catalog"]], "metrics": {"note": SENTINELS["catalog"]},
                              "rule": SENTINELS["catalog"], "version": SENTINELS["catalog"]} if ds["rule"] else None
            reports.append(doc)
            raw = made_up.RAW[ds["drug"]]
            for r in rows:
                call = "Resistant" if r["resistant"] else "Susceptible"
                records.append(made_up.record(r["id"], raw, call=call))
                if r["intermediate"]:
                    records.append(made_up.record(r["id"], raw, call="Intermediate"))
                labels.append({"genome_id": r["id"], "species": species, "antibiotic": ds["drug"],
                               "reported": call[0], "conflict": False, "label": call[0], "label_from_mic": False})
                genomes.append({"genome_id": r["id"], "bioproject_accession": bioproject(k + 1, ds["project"]),
                                "mlst": "made-up", "date_inserted": made_up.STAMP})
                if r["hc20"] is not None:
                    hc20.append({"genome_id": r["id"], "genome_name": "made-up", "mlst": "made-up",
                                 "cgmlst_hc20": r["hc20"]})
    return {"reports": reports, "records": records, "labels": labels, "genomes": genomes, "hc20": hc20,
            "drug_names": made_up.DRUG_NAMES}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def set_pins(module, pins):
    """Set module constants; returns their old values to put back."""
    old = {k: getattr(module, k) for k in pins}
    for k, v in pins.items():
        setattr(module, k, v)
    return old


def stage1_run(stage1, top, release, out):
    """Write the release under <top> and run stage 1 into <out> with its pins set to the made-up files'."""
    argv, pins = made_up.write_release(top, release, hc20=True)
    code = made_up.run(stage1, argv + ["--out", str(out)], pins)
    return code, pins


def stage2_argv(top, stage1_out, out):
    top, stage1_out = Path(top), Path(stage1_out)
    return ["--reports", str(top / "reports"), "--rows-a", str(stage1_out / "local" / "reading-A.tsv"),
            "--rows-b", str(stage1_out / "local" / "reading-B.tsv"),
            "--hc20", str(top / "genomes.jsonl"), "--out", str(out)]


def stage2_pins(top, stage1_out):
    """The pins of stage 2 for a made-up run: the row lists' sha256, the count and typing tables' paths and sha256,
    the genome table's sha256 and stage 1's sums list of reports. Two dicts: stage 2's own constants, and those of
    the stage 1 module stage 2 imports."""
    top, stage1_out = Path(top), Path(stage1_out)
    export = stage1_out / "export"
    own = {"ROW_LIST_SHA256": {"A": sha(stage1_out / "local" / "reading-A.tsv"),
                               "B": sha(stage1_out / "local" / "reading-B.tsv")},
           "COUNTS": export / "counts.tsv", "COUNTS_SHA256": sha(export / "counts.tsv"),
           "TYPING_TABLE": export / "typing.tsv", "TYPING_TABLE_SHA256": sha(export / "typing.tsv"),
           "HC20_SHA256": sha(top / "genomes.jsonl")}
    inner = {"REPORT_SUMS": top / "reports.sha256", "REPORT_SUMS_SHA256": sha(top / "reports.sha256")}
    return own, inner


def stage2_run(stage2, argv, own, inner):
    """Run stage2.main with its pins set; put them back after."""
    old_own, old_inner = set_pins(stage2, own), set_pins(stage2.stage1, inner)
    try:
        return stage2.main(argv)
    finally:
        set_pins(stage2, old_own)
        set_pins(stage2.stage1, old_inner)


def load_scripts():
    return _load("matched_accuracy_stage1_for_made_up2", STAGE1), _load("matched_accuracy_stage2", STAGE2)


def main(argv):
    species, seed = 8, None
    args = list(argv)
    if not args:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    top = Path(args.pop(0))
    while args:
        flag = args.pop(0)
        if flag == "--species" and args:
            species = int(args.pop(0))
        elif flag == "--seed" and args:
            seed = int(args.pop(0))
        else:
            print(__doc__.strip().splitlines()[2], file=sys.stderr)
            return 2
    plan = designed()[:species] if seed is None else randomised(seed, species)
    stage1, stage2 = load_scripts()
    code, _ = stage1_run(stage1, top / "in", build(plan), top / "stage1")
    if code != 0:
        return code
    own, inner = stage2_pins(top / "in", top / "stage1")
    for k, v in sorted(own.items()):
        print(f"made-up pin: {k} = {v if not isinstance(v, Path) else v.name}", flush=True)
    return stage2_run(stage2, stage2_argv(top / "in", top / "stage1", top / "stage2"), own, inner)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
