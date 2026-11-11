"""A stand-in for the grm binary in the run script's tests: `fetch fasta`,
`species-check`, `count` and `predict`, each writing the files the real
command writes in the shape the run script reads, and nothing else. It
shares no code with score.py: `predict` applies the final models itself.
Behaviour is set by environment variables:

  FAKE_MISSING       comma-separated ids the database answers with no sequence
  FAKE_TRANSPORT     comma-separated ids whose fetch fails on the connection
  FAKE_SILENT        comma-separated ids fetch neither gets nor records
  FAKE_VERDICTS      JSON {id: verdict} for species-check (default pass)
  FAKE_DISTANCES     JSON {id: distance} for species-check (default 0.01)
  FAKE_UNJUDGED      comma-separated ids species-check leaves out
  FAKE_REFERENCE     a reference id species-check writes for every species
  FAKE_UNCOUNTED     comma-separated ids count fails on
  FAKE_KMERS         a JSON file {id: [k-mer, ...]}: what each genome holds
  FAKE_GENES         a JSON file {id: [gene call, ...]}: count's gene calls
  FAKE_DIGESTS       a JSON file {id: digest}: count's k-mer set digests
  FAKE_CALLS_VERSION the catalog version count records (default the root's)
  FAKE_PREDICT_WRONG comma-separated ids whose model calls predict flips
  FAKE_PREDICT_K     a k predict writes instead of the report's
"""

import hashlib
import json
import math
import os
import sys
from pathlib import Path


def ids_of(root, slug):
    p = root / "species" / slug / "genomes.txt"
    return [x.strip() for x in p.read_text().splitlines() if x.strip()] if p.exists() else []


def listed(name):
    return [x for x in os.environ.get(name, "").split(",") if x]


def json_file(name):
    p = os.environ.get(name)
    return json.loads(Path(p).read_text()) if p else {}


def fetch(root, cmd):
    slug = cmd[cmd.index("--species") + 1]
    missing, transport, silent = listed("FAKE_MISSING"), listed("FAKE_TRANSPORT"), listed("FAKE_SILENT")
    fails = []
    for g in ids_of(root, slug):
        if g in missing:
            fails.append(f"{g}\t{g}: giving up after 5 attempts: empty FASTA response\n")
        elif g in transport:
            fails.append(f"{g}\t{g}: giving up after 5 attempts: error sending request: connection refused\n")
        elif g in silent:
            continue
        else:
            dst = root / "fasta" / f"{g}.fna.gz"
            if not dst.exists():
                dst.write_bytes(b"fasta " + g.encode())
    if fails:
        (root / "species" / slug / "fasta_failures.tsv").write_text("".join(fails))
    print(f"fetched; {len(fails)} failed")


def species_check(root):
    verdicts = json.loads(os.environ.get("FAKE_VERDICTS", "{}"))
    distances = json.loads(os.environ.get("FAKE_DISTANCES", "{}"))
    unjudged = listed("FAKE_UNJUDGED")
    path = root / "species-check" / "species_check.json"
    doc = json.loads(path.read_text())
    rows = []
    for sp in json.loads((root / "curated" / "species.json").read_text()):
        for g in ids_of(root, sp["slug"]):
            if g in unjudged:
                continue
            v = verdicts.get(g, "pass")
            d = distances.get(g, 0.01)
            rows.append({"genome_id": g, "species": sp["species"], "verdict": v, "gc": 0.33, "gc_z": 0.0,
                         "distance": d, "ani": 1.0 - d, "nearest_species": None, "nearest_distance": None,
                         "inside_species": None,
                         "reason": "far from the reference" if v == "fail" else None})
    doc["genomes"] = rows
    if os.environ.get("FAKE_REFERENCE"):
        for sp in doc["species"]:
            sp["reference"] = os.environ["FAKE_REFERENCE"]
    path.write_text(json.dumps(doc))
    print(f"judged {len(rows)}")


def count(root, cmd):
    uncounted = listed("FAKE_UNCOUNTED")
    kmers, genes, digests = json_file("FAKE_KMERS"), json_file("FAKE_GENES"), json_file("FAKE_DIGESTS")
    version = root / "catalog" / "version.txt"
    cv = os.environ.get("FAKE_CALLS_VERSION") or \
        (version.read_text().strip() if version.exists() else "none") + "|group-calls-v2"
    for slug in cmd[cmd.index("--species") + 1].split(","):
        d = root / "species" / slug
        (d / "ksets").mkdir(parents=True, exist_ok=True)
        out = {"genomes": [], "missing_fasta": [], "failed": []}
        for g in ids_of(root, slug):
            if not (root / "fasta" / f"{g}.fna.gz").exists():
                out["missing_fasta"].append(g)
            elif g in uncounted:
                out["failed"].append([g, "truncated gzip"])
            else:
                (d / "ksets" / f"{g}.kset").write_text(json.dumps({"genome": g, "kmers": kmers.get(g, [])}))
                digest = digests[g] if g in digests else "fnv1a64:" + hashlib.sha256(g.encode()).hexdigest()[:16]
                out["genomes"].append({"genome_id": g, "catalog_version": cv, "genes": genes.get(g, []),
                                       "kset_digest": digest})
        (d / "counted.json").write_text(json.dumps(out))
    print("counted")


def apply_model(report, learner, found):
    """(call, score) of one final model, worked out here from the report."""
    m = report["models"]
    if learner == "scm":
        holds = [found[r["feature"]] == r["presence"] for r in m["scm"]["rules"]]
        return ("resistant" if (all(holds) if m["scm"]["kind"] == "conjunction" else any(holds))
                else "susceptible"), None
    if learner == "cart":
        node = m["cart"]["nodes"][0]
        while node.get("split") is not None and node["depth"] < m["cart_depth"]:
            s = node["split"]
            node = m["cart"]["nodes"][s["present"] if found[s["feature"]] else s["absent"]]
        return ("resistant" if node["value"] >= 0.5 else "susceptible"), node["value"]
    lasso = m["lasso"]
    if not lasso["coef"]:
        return None, None
    p = 1.0 / (1.0 + math.exp(-(lasso["intercept"] + sum(w for f, w in lasso["coef"] if found[f]))))
    return ("resistant" if p >= 0.5 else "susceptible"), p


def model_features(report, learner):
    m = report["models"]
    if learner == "scm":
        return [r["feature"] for r in m["scm"]["rules"]]
    if learner == "cart":
        out = []
        for n in m["cart"]["nodes"]:
            if n.get("split") is not None and n["depth"] < m["cart_depth"] and n["split"]["feature"] not in out:
                out.append(n["split"]["feature"])
        return out
    return None if m.get("lasso") is None else [c[0] for c in m["lasso"]["coef"]]


def predict(cmd):
    report_path = Path(cmd[cmd.index("--report") + 1])
    out = Path(cmd[cmd.index("--out") + 1])
    ksets = [Path(cmd[i + 1]) for i, a in enumerate(cmd) if a == "--kset"]
    data = report_path.read_bytes()
    report = json.loads(data)
    feats = {f["feature"]: f for f in report["features"]}
    k = report["counts"]["k"]
    learners = []
    for learner in ("scm", "cart", "lasso"):
        ids = model_features(report, learner)
        if ids is not None and all(f in feats and len(feats[f]["kmer"]) == k for f in ids):
            learners.append(learner)
    wrong = listed("FAKE_PREDICT_WRONG")
    flip = {"resistant": "susceptible", "susceptible": "resistant", None: None}
    rows = []
    for ks in ksets:
        doc = json.loads(ks.read_text())
        g, held = doc["genome"], set(doc["kmers"])
        for learner in learners:
            ids = model_features(report, learner)
            found = {f: feats[f]["kmer"] in held for f in ids}
            call, score = apply_model(report, learner, found)
            if g in wrong:
                call = flip[call]
            rows.append({"contract": "grm-predict-1", "dataset": report["id"],
                         "report_sha256": hashlib.sha256(data).hexdigest(),
                         "input_fingerprint": report["input_fingerprint"],
                         "evaluation_revision": report["evaluation"]["revision"],
                         "k": int(os.environ.get("FAKE_PREDICT_K") or k), "genome": g, "source": "kset",
                         "assembly": None, "learner": learner, "call": call, "score": score,
                         "features": [{"feature": f, "kmer": feats[f]["kmer"], "found": found[f],
                                       "n_kmers": feats[f]["n_kmers"],
                                       "kmers_found": feats[f]["n_kmers"] if found[f] else 0} for f in ids],
                         "partial_features": 0})
    out.write_text("".join(json.dumps(r, sort_keys=True) + "\n" for r in rows))


def main(argv):
    if argv[:1] == ["predict"]:
        predict(argv[1:])
        return 0
    root = Path(argv[argv.index("--data") + 1])
    cmd = argv[argv.index("--data") + 2:]
    if cmd[:2] == ["fetch", "fasta"]:
        fetch(root, cmd)
    elif cmd[:1] == ["species-check"]:
        species_check(root)
    elif cmd[:1] == ["count"]:
        count(root, cmd)
    else:
        print(f"fake grm: unsupported {cmd}", file=sys.stderr)
        return 64
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
