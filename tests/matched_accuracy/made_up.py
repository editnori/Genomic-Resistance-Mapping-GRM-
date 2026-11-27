"""A SYNTHETIC release for stage 1 of the matched accuracy protocol, and the writer of every input it reads.

  python3 -I tests/matched_accuracy/made_up.py <top> <out> [--no-hc20]

writes the made-up release under <top> and runs scripts/matched-accuracy/stage1.py on it, with the script's
pinned sha256 values and line count set to those of the made-up files (a real run checks the real ones).

Everything here is made up: the species, the antibiotics, the genome ids (gz-<species>-<n>, which is not the
shape of a real genome id), the BioProject accessions, the lab records, the label table, the genome tables and
every score. No row comes from a real file. A BioProject here has an accession's shape, PRJ, two letters and
digits, because ruling R53 reads a value of any other shape as no BioProject; the test needs the shape, so
bioproject() builds it from pieces, with the letters ZZ, which no registry uses, and the digits 9, the
species' number and the project's number.

The made-up release holds 12 reports in 4 species. Most datasets follow one plan of 10 rows:

  position  1  2  3  4  5  6  7  8  9  10
  fold      0  0  1  1  2  2  3  3  3  4
  lineage  10 10 11 11 12 12 13 13 13 14
  call      R  S  R  S  R  S  R  S  R  S

Each dataset changes the plan where its note says. The counts each dataset gives are worked by hand in
tests/matched_accuracy/test_made_up.py.
"""

import gzip
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "matched-accuracy" / "stage1.py"
FIVE = ("scm", "cart", "rf", "gbt", "lasso")
PLAN = [(0, 10, True), (0, 10, False), (1, 11, True), (1, 11, False), (2, 12, True), (2, 12, False),
        (3, 13, True), (3, 13, False), (3, 13, True), (4, 14, False)]
STAMP = "2025-06-22T08:30:00Z"
ASSAY = ("Broth dilution", "CLSI", "2018")
# Values no row of the loader's output may carry: each score and each value of roc, metrics, conformal and no_call.
SENTINELS = ("roc-sentinel", "metrics-sentinel", "conformal-sentinel", "no-call-sentinel")
SCORE = 0.7071067


def load_script():
    spec = importlib.util.spec_from_file_location("matched_accuracy_stage1", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def gid(tag, n):
    return f"gz-{tag}-{n:02d}"


SPECIES_NUMBER = {"al": 1, "be": 2, "ga": 3, "de": 4, "zz": 8}


def bioproject(tag, n):
    """A made-up BioProject with an accession's shape (see the module's note), such as BioProject 1 of alpha."""
    return "PRJ" + "ZZ" + f"9{SPECIES_NUMBER[tag]}{n:02d}"


def plan_rows(tag, first, mlst="ST-1", change=None):
    """The plan's 10 rows for genomes first..first+9; `change` maps a position (1 to 10) to the fields it sets."""
    change = change or {}
    rows = []
    for p, (fold, lineage, call) in enumerate(PLAN, 1):
        r = {"id": gid(tag, first + p - 1), "resistant": call, "fold": fold, "lineage": lineage, "mlst": mlst,
             "from_mic": False}
        r.update(change.get(p, {}))
        rows.append(r)
    return rows


def learner(name, n, k):
    return {"conformal": {"alpha": 0.1, "note": SENTINELS[2]}, "fold_balanced_accuracy": [0.5] * 5,
            "fold_complexity": [3] * 5, "learner": name,
            "metrics": {"balanced_accuracy": 0.5, "note": SENTINELS[1]},
            "no_call": {"share": 0.25, "note": SENTINELS[3]},
            "roc": {"auc": 0.5, "points": [[0.0, 0.0], [1.0, 1.0]], "note": SENTINELS[0]},
            "scores": [round(SCORE + k / 1000 + i / 100000, 7) for i in range(n)], "seconds": 1.25}


def report(dataset, species, slug, antibiotic, rows, copies=None, learners=FIVE + ("catalog",)):
    """One report in the release's shape (SHAPES.md): every top-level key, every row key, two splits."""
    n = len(rows)
    names = list(learners)
    lineage = {"scheme": "lineage", "folds": 5, "n_groups": len({r["lineage"] for r in rows}),
               "learners": [learner(x, n, k) for k, x in enumerate(names)]}
    if names:
        lineage["test_fold"] = [r["fold"] for r in rows]
    else:
        lineage["unavailable_reason"] = "fewer lineage clusters than folds"
    rand = {"scheme": "random", "folds": 5, "n_groups": n, "test_fold": [i % 5 for i in range(n)],
            "learners": [learner(x, n, 10 + k) for k, x in enumerate(FIVE + ("catalog",))]}
    doc = {
        "antibiotic": antibiotic, "breakpoint": None, "catalog": {"rules": []}, "catalog_groups": [],
        "counts": {"features": 4, "k": 31, "kmers": 100, "labels_from_mic": sum(r["from_mic"] for r in rows),
                   "labels_from_mic_excluded": 0, "lineages": len({r["lineage"] for r in rows}), "min_class": 4,
                   "min_support": 3, "patterns": 4, "resistant": sum(r["resistant"] for r in rows), "samples": n,
                   "susceptible": sum(not r["resistant"] for r in rows)},
        "drug_class": "made-up class", "eval_config": {"seed": 42}, "evaluation": {"revision": 14},
        "features": [], "id": dataset, "input_fingerprint": "made-up", "lineage_threshold": 0.02, "mic": None,
        "models": [], "panels": [],
        "samples": [{"catalog_resistant": None, "country": "Nowhere", "from_mic": r["from_mic"], "groups": [],
                     "host": "made-up host", "id": r["id"], "lineage": r["lineage"], "mic": None, "mlst": r["mlst"],
                     "resistant": r["resistant"], "year": "2019"} for r in rows],
        "schema": 1, "seconds": 2.5, "species": species, "species_slug": slug, "splits": [lineage, rand],
        "unsupported_assays": [],
    }
    if copies:
        doc["kset_duplicates"] = copies
    return doc


def record(genome, drug, call="Resistant", method=ASSAY[0], standard=ASSAY[1], year=ASSAY[2], stamp=STAMP,
           **more):
    r = {"genome_id": genome, "genome_name": "made-up", "taxon_id": "1", "antibiotic": drug,
         "resistant_phenotype": call, "measurement": None, "measurement_sign": None, "measurement_value": None,
         "measurement_unit": None, "laboratory_typing_method": method, "laboratory_typing_platform": None,
         "vendor": None, "testing_standard": standard, "testing_standard_year": year, "source": None, "pmid": None,
         "date_inserted": stamp}
    r.update(more)
    return r


DRUG_NAMES = {"names": [
    {"raw": "Drug A", "kind": "drug", "canonical": "drug_a", "records": 0, "class": None},
    {"raw": "Drug B", "kind": "drug", "canonical": "drug_b", "records": 0, "class": None},
    {"raw": "Drug C", "kind": "drug", "canonical": "drug_c", "records": 0, "class": None},
    {"raw": "Drug Class Q", "kind": "class", "canonical": None, "records": 0, "class": "made-up class"},
]}
RAW = {"drug_a": "Drug A", "drug_b": "Drug B", "drug_c": "Drug C"}


def made_up():
    """The made-up release as plain data: reports, records, label table rows, genome table rows, HC20 rows."""
    sp = {"al": ("Synthia alpha", "synthia_alpha"), "be": ("Synthia beta", "synthia_beta"),
          "ga": ("Synthia gamma", "synthia_gamma"), "de": ("Synthia delta", "synthia_delta")}
    datasets = []   # (tag, drug, rows, copy groups, whether the lineage split holds learners)

    # alpha: all three datasets on genomes 1 to 10, alpha's BioProject 1; genome 10's is a placeholder dash.
    datasets.append(("al", "drug_a", plan_rows("al", 1), [[6, 8]], True))      # a copy group: positions 7 and 9
    datasets.append(("al", "drug_b", plan_rows("al", 1), None, True))          # positions 3 and 8 add an I record
    datasets.append(("al", "drug_c", plan_rows("al", 1), None, True))          # the plan
    # beta: drug_a on genomes 1 to 10 (BioProject 1); drug_b and drug_c on genomes 11 to 20 (BioProject 2).
    datasets.append(("be", "drug_a", plan_rows("be", 1), None, True))          # genome 1 blank, genome 5 two values
    datasets.append(("be", "drug_b", plan_rows("be", 11), None, True))         # genome 11 not in the genome table
    datasets.append(("be", "drug_c", plan_rows("be", 11, change={9: {"from_mic": True}}), None, True))
    # gamma: one block of 10 genomes a dataset, BioProjects 1, 2 and 3; genome 1 has a second, blank row.
    datasets.append(("ga", "drug_a", plan_rows("ga", 1), None, True))          # two tuples; case; blank years
    datasets.append(("ga", "drug_b", plan_rows("ga", 11), None, True))         # no record; blank method; blank standard
    datasets.append(("ga", "drug_c", plan_rows("ga", 21), None, True))         # position 1 adds an I record
    # delta: drug_a on genomes 1 to 10 (1 to 6 BioProject 1, 7 to 10 BioProject 2); drug_b on 11 to 22
    # (BioProject 2), 12 rows; drug_c on 23 to 30 without learners, and the second copy lacks genome 30.
    # drug_a's position 10 lacks an MLST type and leaves by the one-outcome rule.
    mlst_a = {1: "ST-1", 2: "ST-1", 3: "ST-3", 4: "ST-4", 5: "ST-5", 6: None, 7: "ST-7", 8: "ST-7", 9: "ST-9",
              10: None}
    datasets.append(("de", "drug_a", plan_rows("de", 1, change={p: {"mlst": m} for p, m in mlst_a.items()}), None, True))
    rows_b = plan_rows("de", 11, mlst="ST-2", change={9: {"mlst": "ST-6"}})
    rows_b += [{"id": gid("de", 21), "resistant": False, "fold": 3, "lineage": 13, "mlst": "ST-6", "from_mic": False},
               {"id": gid("de", 22), "resistant": True, "fold": 4, "lineage": 14, "mlst": "ST-2", "from_mic": False}]
    datasets.append(("de", "drug_b", rows_b, None, True))
    rows_c = [{"id": gid("de", 23 + i), "resistant": i % 2 == 0, "fold": i % 5, "lineage": 20 + i, "mlst": "ST-2",
               "from_mic": False} for i in range(8)]
    datasets.append(("de", "drug_c", rows_c, None, False))

    reports, records, labels = [], [], []
    for tag, drug, rows, copies, has in datasets:
        name, slug = sp[tag]
        dataset = f"{slug}__{drug}"
        reports.append(report(dataset, name, slug, drug, rows, copies, FIVE + ("catalog",) if has else ()))
        for p, r in enumerate(rows, 1):
            call = "Resistant" if r["resistant"] else "Susceptible"
            g = r["id"]
            if r["from_mic"]:
                recs = [record(g, RAW[drug], call="", measurement_value="4", measurement_unit="mg/L")]
                labels.append({"genome_id": g, "species": name, "antibiotic": drug, "reported": None,
                               "conflict": False, "label": "R" if r["resistant"] else "S", "label_from_mic": True})
            else:
                recs = [record(g, RAW[drug], call=call)]
                labels.append({"genome_id": g, "species": name, "antibiotic": drug, "reported": call[0],
                               "conflict": False, "label": call[0], "label_from_mic": False})
            if (tag, drug) == ("al", "drug_b") and p in (3, 8):
                recs.append(record(g, RAW[drug], call="Intermediate"))
            if (tag, drug) == ("ga", "drug_c") and p == 1:
                recs.append(record(g, RAW[drug], call="Intermediate"))
            if (tag, drug) == ("ga", "drug_a"):
                if p == 3:
                    recs.append(record(g, RAW[drug], call=call, method="Disk diffusion"))
                if p == 5:
                    recs = [record(g, RAW[drug], call=call, standard="clsi"),
                            record(g, " " + RAW[drug] + " ", call=call, method=" Broth dilution ", standard="CLSI ")]
                if p in (7, 8):
                    recs = [record(g, RAW[drug], call=call, year="")]
            if (tag, drug) == ("ga", "drug_b"):
                if p == 1:
                    recs = []
                if p == 2:
                    recs.append(record(g, RAW[drug], call=call))
                if p == 3:
                    recs = [record(g, RAW[drug], call=call, method="")]
                if p == 6:
                    recs = [record(g, RAW[drug], call=call, standard=None)]
            records += recs
    # A record of a drug class and one of a genome in no report: both read and passed over.
    records.append(record(gid("al", 1), "Drug Class Q", call="Susceptible"))
    records.append(record(gid("zz", 1), "Drug A", call="Resistant"))
    # A record whose insertion time is a list of two texts, each before the seal.
    records[0]["date_inserted"] = [STAMP, "2026-11-05T19:09:59Z"]

    project = {}
    for n in range(1, 11):
        project[gid("al", n)] = bioproject("al", 1)
        project[gid("be", n)] = bioproject("be", 1)
        project[gid("ga", n)] = bioproject("ga", 1)
        project[gid("de", n)] = bioproject("de", 1 if n <= 6 else 2)
    for n in range(11, 21):
        project[gid("be", n)] = bioproject("be", 2)
        project[gid("ga", n)] = bioproject("ga", 2)
    for n in range(21, 31):
        project[gid("ga", n)] = bioproject("ga", 3)
    for n in range(11, 31):
        project[gid("de", n)] = bioproject("de", 2)
    project[gid("al", 10)] = "-"
    project[gid("be", 1)] = ""
    project[gid("be", 5)] = bioproject("be", 1) + ", " + bioproject("be", 9)
    project[gid("zz", 1)] = bioproject("zz", 1)
    copy = {g: b for g, b in project.items() if g not in (gid("be", 11), gid("de", 30))}
    genomes = [{"genome_id": g, "bioproject_accession": b, "mlst": "made-up", "date_inserted": STAMP}
               for g, b in sorted(copy.items())]
    # A second, blank row of gamma's genome 1, which the reader skips; a row whose insertion time is a list.
    genomes.insert(genomes.index(next(r for r in genomes if r["genome_id"] == gid("ga", 1))) + 1,
                   {"genome_id": gid("ga", 1), "bioproject_accession": " ", "mlst": "made-up",
                    "date_inserted": STAMP})
    next(r for r in genomes if r["genome_id"] == gid("al", 2))["date_inserted"] = [STAMP, "2025-07-01T00:00:00Z"]

    hc20 = {}
    for g in sorted(set(project) | {gid("be", 11)}):
        hc20[g] = "hc-1"
    for n, h in {1: "hc-1", 2: "hc-1", 3: "hc-3", 4: "hc-3", 5: "hc-3", 6: "hc-3", 7: "hc-7", 8: "hc-8",
                 9: "hc-9"}.items():
        hc20[gid("de", n)] = h
    for n in range(11, 31):
        hc20[gid("de", n)] = "hc-2"
    hc20[gid("de", 14)] = ""
    del hc20[gid("de", 10)]
    hc20_rows = [{"genome_id": g, "genome_name": "made-up", "mlst": "made-up", "cgmlst_hc20": h}
                 for g, h in sorted(hc20.items())]
    return {"reports": reports, "records": records, "labels": labels, "genomes": genomes, "hc20": hc20_rows,
            "drug_names": DRUG_NAMES}


def jsonl_bytes(rows):
    return "".join(json.dumps(r, sort_keys=True) + "\n" for r in rows).encode("utf-8")


def gz_bytes(data):
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as f:
        f.write(data)
    return buf.getvalue()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_release(top, release, hc20=True):
    """Write a release's files under <top>; return (argv without --out, pins for the script)."""
    top = Path(top)
    reports = top / "reports"
    reports.mkdir(parents=True, exist_ok=True)
    lines = []
    for doc in release["reports"]:
        data = json.dumps(doc, sort_keys=True).encode("utf-8")
        (reports / f"{doc['id']}.json").write_bytes(data)
        lines.append(f"{sha(data)}  ./datasets/{doc['id']}.json\n")
    sums = "".join(sorted(lines, key=lambda line: line.split("  ")[1])).encode("utf-8")
    (top / "reports.sha256").write_bytes(sums)
    files = {"genome_amr.jsonl": release.get("records_bytes") or jsonl_bytes(release["records"]),
             "drug_names.json": json.dumps(release["drug_names"], sort_keys=True).encode("utf-8"),
             "phenotypes.jsonl": jsonl_bytes(release["labels"]),
             "master-genome-meta.jsonl.gz": release.get("genomes_bytes") or gz_bytes(jsonl_bytes(release["genomes"])),
             "genomes.jsonl": jsonl_bytes(release["hc20"])}
    for name, data in files.items():
        (top / name).write_bytes(data)
    argv = ["--reports", str(reports), "--records", str(top / "genome_amr.jsonl"),
            "--drug-names", str(top / "drug_names.json"), sha(files["drug_names.json"]),
            "--labels", str(top / "phenotypes.jsonl"), sha(files["phenotypes.jsonl"]),
            "--genomes", str(top / "master-genome-meta.jsonl.gz")]
    argv += ["--hc20", str(top / "genomes.jsonl"), sha(files["genomes.jsonl"])] if hc20 else ["--hc20", "-"]
    records = files["genome_amr.jsonl"]
    pins = {"REPORT_SUMS": top / "reports.sha256", "REPORT_SUMS_SHA256": sha(sums),
            "RECORDS_SHA256": sha(records), "RECORDS_LINES": records.count(b"\n"),
            "GENOMES_SHA256": sha(files["master-genome-meta.jsonl.gz"])}
    return argv, pins


def run(module, argv, pins):
    """Run stage1.main with its pins set to the made-up files'; put the pins back after."""
    old = {k: getattr(module, k) for k in pins}
    try:
        for k, v in pins.items():
            setattr(module, k, v)
        return module.main(argv)
    finally:
        for k, v in old.items():
            setattr(module, k, v)


def main(argv):
    if len(argv) not in (2, 3) or (len(argv) == 3 and argv[2] != "--no-hc20"):
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    top, out = Path(argv[0]), Path(argv[1])
    args, pins = write_release(top, made_up(), hc20=len(argv) == 2)
    module = load_script()
    for k, v in pins.items():
        print(f"made-up pin: {k} = {v if k != 'REPORT_SUMS' else Path(v).name}", flush=True)
    return run(module, args + ["--out", str(out)], pins)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
