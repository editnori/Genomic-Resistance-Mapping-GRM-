"""Shared loaders and SYNTHETIC fixtures for the mechanism split (scripts/mechanism, docs/mechanism-split.md).

Every report, table, genome id, BioProject, database row and count here is
made up in the test. No genome, call or label comes from a real dataset, and
the species, dataset and gene names are made up or used only as shapes.
"""

import csv
import gzip
import hashlib
import importlib.util
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
KIT = ROOT / "scripts" / "mechanism"
SPECIES = "synthia_coli"
SPECIES_NAME = "Synthia coli"
DATASET = "synthia_coli__drug_a"
RULE = [["C:QUINOLONE"]]
LEARNERS = ("scm", "cart", "rf", "gbt", "lasso")
TABLE_HEADER = ["Protein id", "Contig id", "Start", "Stop", "Strand", "Element symbol", "Element name", "Scope",
                "Type", "Subtype", "Class", "Subclass", "Method", "Target length", "Reference sequence length",
                "% Coverage of reference", "% Identity to reference", "Alignment length",
                "Closest reference accession", "Closest reference name", "HMM accession", "HMM description"]

if str(KIT) not in sys.path:
    sys.path.insert(0, str(KIT))


def load(name):
    """The script module itself, as the scripts import each other (top-level names)."""
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, KIT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# A table row per kind of genome; all made up.
ROWS = {
    "point": [("gyrS_S83L", "core", "AMR", "POINT", "QUINOLONE", "QUINOLONE", "POINTX")],
    "gene": [("qnrS1", "core", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "EXACTX")],
    "gene_b": [("qnrB4", "core", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "BLASTX")],
    "both": [("gyrS_S83L", "core", "AMR", "POINT", "QUINOLONE", "QUINOLONE", "POINTX"),
             ("qnrS1", "core", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "EXACTX")],
    "neither": [("tetS", "core", "AMR", "AMR", "TETRACYCLINE", "TETRACYCLINE", "EXACTX")],
    "partial": [("qnrB4", "core", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "PARTIALX")],
    "disrupt": [("synR", "core", "AMR", "POINT_DISRUPT", "QUINOLONE", "QUINOLONE", "BLASTX")],
    "plus_only": [("qnrS1", "plus", "AMR", "AMR", "QUINOLONE", "QUINOLONE", "EXACTX")],
    "empty": [],
}
DB_CATALOG = [
    # allele, gene_family, class, subclass
    ("qnrS1", "qnrS", "QUINOLONE", "QUINOLONE"),
    ("qnrB4", "qnrB", "QUINOLONE", "QUINOLONE"),
    ("", "tetS", "TETRACYCLINE", "TETRACYCLINE"),
    ("", "synR", "QUINOLONE", "QUINOLONE"),
    ("blaKPC-2", "blaKPC", "BETA-LACTAM", "CARBAPENEM"),
    ("blaX-1", "blaX", "BETA-LACTAM", "CEPHALOSPORIN"),
]
DB_HIERARCHY = [("blaX-9", "blaX"), ("blaX", "BETA-LACTAM"), ("BETA-LACTAM", "AMR")]


def write_table(path, kind):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(TABLE_HEADER) + "\n")
        for symbol, scope, typ, subtype, cls, subclass, method in ROWS[kind]:
            cells = dict.fromkeys(TABLE_HEADER, "NA")
            cells.update({"Element symbol": symbol, "Element name": symbol, "Scope": scope, "Type": typ,
                          "Subtype": subtype, "Class": cls, "Subclass": subclass, "Method": method})
            f.write("\t".join(cells[h] for h in TABLE_HEADER) + "\n")


def write_db(directory):
    d = pathlib.Path(directory)
    d.mkdir(parents=True, exist_ok=True)
    with open(d / "ReferenceGeneCatalog.txt", "w", encoding="utf-8") as f:
        f.write("allele\tgene_family\twhitelisted_taxa\tproduct_name\tscope\ttype\tsubtype\tclass\tsubclass\n")
        for allele, fam, cls, sub in DB_CATALOG:
            f.write(f"{allele}\t{fam}\t\tmade up\tcore\tAMR\tAMR\t{cls}\t{sub}\n")
    with open(d / "ReferenceGeneHierarchy.txt", "w", encoding="utf-8") as f:
        f.write("node_id\tparent_node_id\tsymbol\n")
        for node, parent in DB_HIERARCHY:
            f.write(f"{node}\t{parent}\t{node}\n")
    return d


def program_calls(kinds):
    """The program's tool_default call per kind under RULE."""
    return [k in ("point", "gene", "gene_b", "both", "partial", "disrupt") for k in kinds]


def confusion(truth, calls):
    c = {"tp": 0, "fn": 0, "tn": 0, "fp": 0}
    for t, p in zip(truth, calls):
        c[("tp" if p else "fn") if t else ("fp" if p else "tn")] += 1
    return c


def report(genomes, dataset=DATASET, scores=None, rule=RULE, drug_class="fluoroquinolone", splits=("random", "lineage")):
    """A SYNTHETIC revision 14 report. `genomes` is a list of dicts with id,
    kind, resistant and optional year, country, lineage; `scores(split,
    learner, i, genome)` gives a learner's score (default: 0.9 for a resistant
    genome, 0.1 otherwise)."""
    scores = scores or (lambda split, learner, i, g: 0.9 if g["resistant"] else 0.1)
    truth = [g["resistant"] for g in genomes]
    catalog = [g.get("catalog", g["resistant"]) for g in genomes]
    out_splits = []
    for split in splits:
        learners = []
        for name in LEARNERS:
            sc = [scores(split, name, i, g) for i, g in enumerate(genomes)]
            learners.append({"learner": name, "scores": sc, "metrics": confusion(truth, [s >= 0.5 for s in sc]),
                             "no_call": {"resistant": 0, "susceptible": 0}})
        sc = [1.0 if c else 0.0 for c in catalog]
        learners.append({"learner": "catalog", "scores": sc, "metrics": confusion(truth, catalog)})
        out_splits.append({"scheme": split, "learners": learners, "folds": 5, "n_groups": 10, "test_fold": []})
    samples = [{"id": g["id"], "resistant": g["resistant"], "catalog_resistant": c, "year": g.get("year"),
                "country": g.get("country"), "lineage": g.get("lineage", 0), "host": None, "mlst": None, "mic": None,
                "from_mic": False, "groups": []} for g, c in zip(genomes, catalog)]
    sp = dataset.split("__")[0]
    return {"id": dataset, "species_slug": sp, "species": SPECIES_NAME if sp == SPECIES else sp,
            "antibiotic": dataset.split("__")[1].replace("_", " "), "drug_class": drug_class,
            "catalog": {"rule": rule, "version": "2026-08-07.1", "groups": [], "metrics": {}},
            "samples": samples, "splits": out_splits}


def lists_row(dataset=DATASET, in_m=False, in_a=True, depends=""):
    sp = dataset.split("__")[0]
    return {"dataset": dataset, "species_slug": sp, "species": SPECIES_NAME if sp == SPECIES else sp,
            "antibiotic": dataset.split("__")[1].replace("_", " "), "in_M": "1" if in_m else "0",
            "in_A": "1" if in_a else "0", "labels_depend_on": depends}


def write_tsv(path, header, rows):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(h, "")) for h in header) + "\n")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))


def write_species(top, datasets, lists=None, adapter_change=None, record_change=None, rules_change=None):
    """Lay out one SYNTHETIC species: reports/, tool/<species>/, lists.tsv,
    drug_rules.tsv, calls.tsv, callers.tsv. `datasets` maps a dataset id to
    (genomes, report kwargs)."""
    amr = load("amr_rule")
    top = pathlib.Path(top)
    (top / "reports").mkdir(parents=True, exist_ok=True)
    calls, record, rules, lrows = [], [], [], []
    written = set()
    for d, (genomes, kw) in datasets.items():
        rep = report(genomes, dataset=d, **kw)
        with open(top / "reports" / f"{d}.json", "w", encoding="utf-8") as f:
            json.dump(rep, f)
        for g in genomes:
            if g["id"] not in written and g["kind"] is not None:
                write_table(top / "tool" / SPECIES / f"{g['id']}.tsv", g["kind"])
                written.add(g["id"])
        prog = []
        for g in genomes:
            if g["kind"] is None:
                prog.append(None)
                continue
            rows = amr.read_table(top / "tool" / SPECIES / f"{g['id']}.tsv")
            prog.append(amr.program_call(rows, rep["catalog"]["rule"]))
        for g, p in zip(genomes, prog):
            call = "no_call" if p is None else "resistant" if p else "susceptible"
            if adapter_change:
                call = adapter_change(d, g, call)
            calls.append({"dataset_id": d, "genome_id": g["id"], "variant": "tool_default", "call": call,
                          "point_screened": "yes", "elements_used": ""})
            calls.append({"dataset_id": d, "genome_id": g["id"], "variant": "genes_core", "call": "susceptible",
                          "point_screened": "yes", "elements_used": ""})
        kept = [(g["resistant"], p) for g, p in zip(genomes, prog) if p is not None]
        c = confusion([k[0] for k in kept], [k[1] for k in kept])
        rec = dict(c, dataset=d, caller="program:tool_default", genomes=len(genomes), scored=len(kept))
        if record_change:
            record_change(rec)
        record.append(rec)
        rule = amr.format_rule(rep["catalog"]["rule"])
        if rules_change:
            rule = rules_change(d, rule)
        rules.append({"dataset": d, "species_slug": rep["species_slug"], "species": rep["species"],
                      "antibiotic": rep["antibiotic"], "drug_class": rep["drug_class"], "rule": rule})
        lrows.append((lists or {}).get(d, lists_row(d)))
    write_tsv(top / "lists.tsv", ["dataset", "species_slug", "species", "antibiotic", "in_M", "in_A",
                                  "labels_depend_on"], lrows)
    write_tsv(top / "drug_rules.tsv", ["dataset", "species_slug", "species", "antibiotic", "drug_class", "rule"],
              rules)
    write_tsv(top / "calls.tsv", ["dataset_id", "genome_id", "variant", "call", "point_screened", "elements_used"],
              calls)
    write_tsv(top / "callers.tsv", ["dataset", "caller", "genomes", "scored", "tp", "fn", "tn", "fp"], record)
    return top


def write_meta(path, bioprojects):
    """A SYNTHETIC BV-BRC genome table: genome id -> BioProject (or None)."""
    with gzip.open(path, "wt", encoding="utf-8") as f:
        for gid, bp in bioprojects.items():
            f.write(json.dumps({"genome_id": gid, "bioproject_accession": bp or ""}) + "\n")
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def genomes(prefix, n, kind, resistant, **extra):
    """n SYNTHETIC genomes with made-up ids <prefix>.<i>."""
    return [dict({"id": f"{prefix}.{i}", "kind": kind, "resistant": resistant}, **extra) for i in range(1, n + 1)]
