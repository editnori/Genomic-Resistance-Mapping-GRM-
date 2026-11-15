"""SYNTHETIC inputs for scripts/other-database/record.py: a small lanes tree, release and data root.

Nothing here is data. Every species gets one synthetic dataset of 60 scored
genomes (30 labelled R, 30 labelled S) and one genome removed at D3; the
pilot species gets a second dataset and four genomes whose labels test the
consolidation rule. S. aureus gets two more datasets on the class line:
edgemycin with exactly 25 R, then five genomes labelled only by the
intermediate wordings, then 30 S; underedgemycin with 24 R. In A. baumannii
the predictions are written and the first label request is made in the
second ADDENDUM2.md was fixed. The figures are chosen so that the cells can
be worked by hand (see test_script.py).

make(root, fault=None) writes the tree under root and returns (lanes,
release, data_root). A fault name breaks one input the way a test of a stop
needs it.
"""

import gzip
import hashlib
import json
import os

SPECIES = [
    ("escherichia_coli", "Escherichia coli"),
    ("salmonella_enterica", "Salmonella enterica"),
    ("klebsiella_pneumoniae", "Klebsiella pneumoniae"),
    ("staphylococcus_aureus", "Staphylococcus aureus"),
    ("campylobacter_jejuni", "Campylobacter jejuni"),
    ("campylobacter_coli", "Campylobacter coli"),
    ("neisseria_gonorrhoeae", "Neisseria gonorrhoeae"),
    ("pseudomonas_aeruginosa", "Pseudomonas aeruginosa"),
    ("acinetobacter_baumannii", "Acinetobacter baumannii"),
]
PILOT = "escherichia_coli"
N = 60  # scored genomes per species: 30 R then 30 S
PLAN_TIME = "2026-01-01T00:00:00Z"
ADD1_TIME = "2026-01-01T01:00:00Z"
ADD2_TIME = "2026-01-01T02:00:00Z"
PRED_TIME = "2026-01-01T03:00:00Z"
LABEL_TIME = "2026-01-01T03:00:05Z"
SAME_SECOND = "acinetobacter_baumannii"  # predictions and first label request in ADDENDUM2.md's own second
# The wordings parse_call reads as intermediate, given in turn to genomes 25 to 29 of edgemycin.
INTERMEDIATE = ("intermediate", "susceptible-dose dependent", "reduced susceptibility")


def write(path, text, mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if path.endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8", newline="") as f:
            f.write(text)
    else:
        with open(path, mode, encoding="utf-8", newline="") as f:
            f.write(text)


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def tsv(cols, rows):
    return "\t".join(cols) + "\n" + "".join("\t".join(str(r[c]) for c in cols) + "\n" for r in rows)


def genomes(k):
    """(assembly, biosample) of the scored genomes of species k, plus the one D3 removes."""
    return [(f"GCA_9{k:02d}{i:04d}.1", f"SAMS{k:02d}{i:04d}") for i in range(N + 1)]


def bioproject(slug, i):
    if slug == PILOT:
        # Two BioProjects tie with 30 genomes each: the lower accession, PRJNA9, is the largest.
        return "PRJNA9" if i % 2 == 0 else "PRJNA10"
    return "PRJNA100" if i < 24 else f"PRJNA{101 + i % 12}"


def scm_calls():
    """SCM on a synthetic dataset: 27 of 30 R called resistant, 28 of 30 S called susceptible."""
    return [1 if i < 27 else 0 for i in range(30)] + [1 if i < 2 else 0 for i in range(30)]


def cart_score(i):
    """CART's p on scored genome i: 0.50 + 0.01 i for the 30 R, 0.01 (i - 29) for the 30 S."""
    return f"{0.5 + 0.01 * i:.2f}" if i < 30 else f"{0.01 * (i - 29):.2f}"


def phenotype(drug, i):
    """The label of scored genome i for a drug: 30 R then 30 S, except on the class line."""
    if drug == "edgemycin":
        return "resistant" if i < 25 else INTERMEDIATE[(i - 25) % len(INTERMEDIATE)] if i < 30 else "susceptible"
    if drug == "underedgemycin":
        return "resistant" if i < 24 else "not defined" if i < 30 else "susceptible"
    return "resistant" if i < 30 else "susceptible"


def chain_log(written, fetched, fault):
    """A run's log: the prediction step, then the label step; the fault log_order writes them the other way."""
    steps = [f"predictions written {written}", f"labels fetched {fetched}"]
    if fault == "log_order":
        steps.reverse()
    return "SYNTHETIC\n" + "\n".join(steps) + "\n"


def lineage(fault):
    """The release's lineage split: 20 R and 20 S samples; SCM recall 16 of 20 and specificity 18 of 20.

    One CART score of an R sample is 0.5 exactly, which the cutoff calls resistant."""
    resistant = [True] * 20 + [False] * 20
    scm = [1.0 if i < 16 else 0.0 for i in range(20)] + [1.0 if i < 2 else 0.0 for i in range(20)]
    cart = [0.8] * 19 + [0.5] + [0.3] * 20
    metrics = {"sensitivity": 0.8 if fault != "lineage_metric" else 0.81, "specificity": 0.9}
    return resistant, [
        {"learner": "scm", "scores": scm, "metrics": metrics},
        {"learner": "cart", "scores": cart, "metrics": {"sensitivity": 1.0, "specificity": 1.0}},
        {"learner": "lasso", "scores": cart, "metrics": {}},
    ]


def make(root, fault=None):
    lanes, release, data = (os.path.join(root, x) for x in ("lanes", "rev-synthetic", "data-root"))
    write(os.path.join(data, "curated", "drug_names.json"), json.dumps({"names": [{"raw": "Synthomycin B", "canonical": "synthomycin", "kind": "drug"}]}))

    # Plan files and their sha256 lines.
    for rel, sums, time in (("outside/PLAN.md", "outside/results/PLAN.sha256", PLAN_TIME), ("outside2/ADDENDUM.md", "outside2/ADDENDUM.sha256", ADD1_TIME),
                            ("outside2/ADDENDUM2.md", "outside2/ADDENDUM2.sha256", ADD2_TIME)):
        path = os.path.join(lanes, rel)
        write(path, f"# SYNTHETIC {os.path.basename(rel)}\n")
        digest = sha(path) if fault != "plan_sha" else "0" * 64
        write(os.path.join(lanes, sums), f"{digest}  {os.path.basename(rel)}\n# fixed {time}, SYNTHETIC\n")

    index = {"evaluation_revision": 14, "datasets": []}
    requests = []
    noacc = []
    for k, (slug, name) in enumerate(SPECIES):
        sdir = os.path.join(lanes, "outside", "results", "pilot") if slug == PILOT else os.path.join(lanes, "outside2", "results", slug)
        gs = genomes(k)
        scored, removed = gs[:N], gs[N]
        drugs = [("synthomycin", "Synthomycin B")]
        if slug == PILOT:
            drugs.append(("amoxicillin/clavulanate", "amoxicillin-clavulanic acid"))
        if slug == "staphylococcus_aureus":
            drugs += [("edgemycin", "edgemycin"), ("underedgemycin", "underedgemycin")]
        # Selected, filter and training files.
        sel = []
        for i, (asm, bs) in enumerate(gs):
            bp = bioproject(slug, i)
            sel.append({"biosample": bs, "assembly": asm, "assembly_keys": "", "runs": "", "bioprojects": bp, "release_drugs": "",
                        "bioproject_shared_with_training": "yes" if bp == "PRJNA100" else "no", "outside_study_overlap": "no", "removed_at": ""})
        write(os.path.join(sdir, "selected.tsv"), tsv(list(sel[0]), sel))
        filt = []
        for i, (asm, _) in enumerate(gs):
            filt.append({"assembly": asm, "nearest_distance": f"{0.0006 + 0.0001 * i:.6f}" if i < N else "0.000400",
                         "removed_at": "" if i < N else ("D3 sequence duplicate" if fault != "d3_count" else "")})
        write(os.path.join(sdir, "filter.tsv"), tsv(["assembly", "nearest_distance", "removed_at"], filt))
        training = [{"genome_id": f"{k}.{j}", "bioprojects": bp} for j, bp in enumerate(["PRJNA100", "PRJNA100", ""])]
        write(os.path.join(sdir, "training-accessions.tsv"), tsv(["genome_id", "bioprojects"], training))
        # Funnel: C1 to P2 label-blind, then the sequence steps.
        c1 = N + 1 + 5 + 20 + 10
        steps = [(f"C1 NCBI AST isolates (BioSamples) of {name}", 0), ("C2 no release drug", 0), ("C3 no GenBank assembly", 10),
                 ("C4 in BV-BRC labelled set", 20), ("D1 training accession", 0), ("D2 training BioSample", 0)]
        if slug != PILOT:
            steps.append(("O1 outside-study candidate", 0))
        steps += [("D4 within-batch duplicate", 0), ("P2 cap", 5 if fault != "funnel" else 6)]
        rows, left = [], c1
        for step, gone in steps:
            left -= gone
            rows.append({"step": step, "removed": gone, "remaining": left})
        if fault == "funnel":
            rows[-1]["remaining"] += 1
        write(os.path.join(sdir, "funnel-label-blind.tsv"), tsv(["step", "removed", "remaining"], rows) + "# SYNTHETIC\n")
        seq = [{"step": "after P2", "removed": 0, "remaining": N + 1}]
        for step, gone in (("F1 FASTA not fetched", 0), ("Q1 QC", 0), ("Q2 species check: not judged", 0), ("Q2 species check: fail", 0), ("D3 sequence duplicate", 1)):
            seq.append({"step": step, "removed": gone, "remaining": seq[-1]["remaining"] - gone})
        write(os.path.join(sdir, "funnel-sequence.tsv"), tsv(["step", "removed", "remaining"], seq))

        # Reports, predictions and labels.
        pins, preds, labels = [], [], []
        for drug, ncbi_name in drugs:
            ds = f"{slug}__{drug.replace('/', '_')}"
            index["datasets"].append({"id": ds, "antibiotic": drug})
            resistant, learners = lineage(fault if (slug, drug) == (PILOT, "synthomycin") else None)
            report = {"id": ds, "samples": [{"resistant": r} for r in resistant], "splits": [{"scheme": "lineage", "learners": learners}]}
            path = os.path.join(release, "datasets", ds + ".json")
            write(path, json.dumps(report))
            digest = sha(path)
            pins.append(f"{digest if fault != 'pin' else '1' * 64}  {ds}.json\n")
            if (slug, drug) == (PILOT, "synthomycin"):
                iv = {"splits": [{"scheme": "lineage", "learners": [{"learner": "scm", "vme": {"lo": 0.05, "hi": 0.4}}]}]}
                write(os.path.join(release, "intervals", ds + ".json"), json.dumps(iv))
            calls = scm_calls()
            for i, (asm, bs) in enumerate(scored):
                preds.append({"dataset": ds, "genome": asm, "learner": "scm", "call": "resistant" if calls[i] else "susceptible", "score": "", "report_sha256": digest})
                # CART: right on every genome, R at p 0.50 to 0.79 and S at p 0.01 to 0.30, so the losses differ.
                preds.append({"dataset": ds, "genome": asm, "learner": "cart", "call": "resistant" if i < 30 else "susceptible",
                              "score": cart_score(i), "report_sha256": digest})
                if slug == PILOT and drug == "synthomycin" and i % 3 == 0:
                    # The lasso on a third of the genomes: 10 R and 10 S, too few for the class line.
                    preds.append({"dataset": ds, "genome": asm, "learner": "lasso", "call": "resistant" if i < 30 else "susceptible",
                                  "score": "0.7" if i < 30 else "0.2", "report_sha256": digest})
                labels.append({"antibiotic": ncbi_name, "biosample_acc": bs.lower() if i == 0 else bs, "phenotype": phenotype(drug, i)})
        if slug == PILOT:
            # Four more scored genomes for the consolidation rule: R and S (conflict), intermediate only, "not defined" only, no row.
            extra = [(f"GCA_9{k:02d}{9000 + j:04d}.1", f"SAMS{k:02d}{9000 + j:04d}") for j in range(4)]
            for j, (asm, bs) in enumerate(extra):
                sel.append({"biosample": bs, "assembly": asm, "assembly_keys": "", "runs": "", "bioprojects": "", "release_drugs": "",
                            "bioproject_shared_with_training": "no", "outside_study_overlap": "no", "removed_at": ""})
                filt.append({"assembly": asm, "nearest_distance": "0.009000", "removed_at": ""})
                for drug, _ in drugs:
                    ds = f"{slug}__{drug.replace('/', '_')}"
                    preds.append({"dataset": ds, "genome": asm, "learner": "scm", "call": "resistant", "score": "", "report_sha256": sha(os.path.join(release, "datasets", ds + ".json"))})
            for phen in ("resistant", "susceptible"):
                labels.append({"antibiotic": "Synthomycin B", "biosample_acc": extra[0][1], "phenotype": phen})
            labels.append({"antibiotic": "Synthomycin B", "biosample_acc": extra[1][1], "phenotype": "intermediate"})
            labels.append({"antibiotic": "Synthomycin B", "biosample_acc": extra[2][1], "phenotype": "not defined"})
            labels.append({"antibiotic": "an antibiotic no release holds", "biosample_acc": extra[3][1], "phenotype": "resistant"})
            write(os.path.join(sdir, "selected.tsv"), tsv(list(sel[0]), sel))
            write(os.path.join(sdir, "filter.tsv"), tsv(["assembly", "nearest_distance", "removed_at"], filt))
            for name_ in ("funnel-label-blind.tsv", "funnel-sequence.tsv"):
                path = os.path.join(sdir, name_)
                with open(path, encoding="utf-8") as f:
                    text = f.read()
                fixed = []
                for line in text.splitlines():
                    parts = line.split("\t")
                    if len(parts) == 3 and parts[2].isdigit():
                        parts[2] = str(int(parts[2]) + 4)
                    fixed.append("\t".join(parts))
                write(path, "\n".join(fixed) + "\n")
        if fault == "stray_genome" and slug == "salmonella_enterica":
            preds.append({"dataset": preds[0]["dataset"], "genome": "GCA_000000000.1", "learner": "scm", "call": "resistant", "score": "", "report_sha256": preds[0]["report_sha256"]})
        write(os.path.join(sdir, "predictions.tsv.gz"), tsv(["dataset", "genome", "learner", "call", "score", "report_sha256"], preds))
        write(os.path.join(sdir, "labels.jsonl.gz"), "".join(json.dumps(x, sort_keys=True) + "\n" for x in labels))
        pins_path = os.path.join(lanes, "outside", "results", "pilot-reports.sha256") if slug == PILOT else os.path.join(sdir, "reports.sha256")
        write(pins_path, "".join(pins))
        noacc.append({"species": name, "labelled_genomes": 100, "neither_biosample_nor_assembly": 25})

        # The order of predictions and label requests.
        scored_bs = [bs for _, bs in scored] + ([f"SAMS{k:02d}{9000 + j:04d}" for j in range(4)] if slug == PILOT else [])
        if slug == PILOT:
            write(os.path.join(lanes, "outside", "results", "logs", "chain.log"), chain_log(PRED_TIME, LABEL_TIME, fault))
            probes = [
                {"utc": "2025-12-31T23:00:00Z", "label": "debug-a", "http": "200", "bytes": "67", "curl_rc": "0", "url": "https://example.invalid/", "error": "POST-test"},
                {"utc": "2025-12-31T23:00:02Z", "label": "debug-b", "http": "200", "bytes": "304", "curl_rc": "0", "url": "https://example.invalid/", "error": "POST-test"},
                {"utc": "2026-01-01T05:00:00Z", "label": "debug-after", "http": "200", "bytes": "9", "curl_rc": "0", "url": "https://example.invalid/", "error": "POST-test"},
                {"utc": LABEL_TIME, "label": "labels 0", "http": "200", "bytes": "100", "curl_rc": "0", "url": "https://example.invalid/", "error": ""},
            ]
            write(os.path.join(lanes, "outside", "results", "probes.tsv"), tsv(list(probes[0]), probes))
        else:
            written, asked = (ADD2_TIME, ADD2_TIME) if slug == SAME_SECOND else (PRED_TIME, LABEL_TIME)
            write(os.path.join(lanes, "outside2", "results", "logs", slug + ".chain.log"), chain_log(written, asked, fault))
            fq = "biosample_acc==[" + ",".join(f'"{b}"' for b in scored_bs) + "]"
            params = json.dumps({"action": "retrieve", "collection": "ast", "fq": fq})
            requests.append([asked, "labels 0", "POST", "https://example.invalid/", "200", "100", params])
            requests.append([LABEL_TIME, "train-iso x 0", "POST", "https://example.invalid/", "200", "100", "{}"])
    write(os.path.join(lanes, "outside2", "results", "requests.tsv"), "".join("\t".join(r) + "\n" for r in requests))
    write(os.path.join(lanes, "outside2", "results", "no-accession-share.tsv"), tsv(["species", "labelled_genomes", "neither_biosample_nor_assembly"], noacc))
    write(os.path.join(release, "index.json"), json.dumps(index))
    if fault == "missing":
        os.remove(os.path.join(lanes, "outside2", "results", "campylobacter_coli", "training-accessions.tsv"))
    return lanes, release, data
