#!/usr/bin/env python3
"""Recompute a GRM report's input_fingerprint from its data root (stdlib only).

usage: python3 -I fingerprint.py DATA_ROOT [--slug mycobacterium_tuberculosis]
           [--check] [--out OUT.json]

`grm learn` stamps every report with `input_fingerprint`
(src/analysis.rs `input_fingerprint`, unchanged since evaluation revision 7):
64-bit FNV-1a (src/fnv.rs) over, in order,

1. `matrix_digest`: the matrix summary (k, genome names, kset digests,
   min_support, pattern and k-mer counts), then every word of patterns.bin,
   support.bin, kmer_offsets.bin and kmers.bin, as stored;
2. the antibiotic (length-prefixed), the dataset's class minimum, its
   breakpoint, the lineage threshold's f32 bits and the catalog calls
   version (none when no catalog loaded);
3. every EvalConfig setting, including the forest and boosting seeds the
   report omits (`#[serde(skip)]`; FOREST_SEED and GBT_SEED below are the
   source defaults `grm learn` runs with);
4. the number of the drug's matrix genomes whose label was inferred from
   MIC (`labels_from_mic_excluded`);
5. per evaluation row in matrix column order (`cohort`: the drug's
   phenotypes on matrix genomes with an unconflicted, reported R or S
   call): genome ID, R flag, MIC, the genome's lineage label, and its
   counted gene calls (family, group, allele, class, subclass, containment).

It reads only the root's files: species/<slug>/matrix/, counted.json,
structure.json, curated/datasets.json, curated/phenotypes.jsonl and
catalog/version.txt. With --check, every saved report under
species/<slug>/results must carry the fingerprint recomputed here, from its
own eval_config; any difference exits 1.
"""

import argparse
import json
import struct
import sys
from pathlib import Path

FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
U64 = (1 << 64) - 1
FOREST_SEED = 17  # src/learn/tree.rs ForestParams::default
GBT_SEED = 23     # src/learn/gbt.rs GbtParams::default
CALLS_SHAPE = "|group-calls-v2"  # src/catalog.rs Catalog::calls_version
CENSOR = {"Exact": 1, "AtMost": 2, "Above": 3}
MATRIX_FILES = ["patterns.bin", "support.bin", "kmer_offsets.bin", "kmers.bin"]


class Fnv:
    """src/fnv.rs: FNV-1a over fixed-width little-endian values."""

    def __init__(self):
        self.h = FNV_OFFSET

    def bytes(self, data):
        h = self.h
        for b in data:
            h = ((h ^ b) * FNV_PRIME) & U64
        self.h = h

    def u64(self, v):
        self.bytes(struct.pack("<Q", v))

    def str(self, s):
        data = s.encode()
        self.u64(len(data))
        self.bytes(data)

    def hex(self):
        return f"{self.h:016x}"


def f32_bits(x):
    return struct.pack("<f", x)


def f64_bits(x):
    return struct.pack("<d", x)


class Inputs:
    """What `grm learn` reads for one species, from DATA_ROOT's files."""

    def __init__(self, data_root: Path, slug: str):
        self.root = Path(data_root)
        sp = self.root / "species" / slug
        self.matrix_dir = sp / "matrix"
        self.summary = json.loads((self.matrix_dir / "matrix.json").read_text())
        self.genomes = self.summary["genomes"]
        self.column = {g: i for i, g in enumerate(self.genomes)}
        counted = json.loads((sp / "counted.json").read_text())
        self.hits = {g["genome_id"]: g.get("genes", []) for g in counted["genomes"]}
        versions = {g["catalog_version"] for g in counted["genomes"]}
        if len(versions) > 1:
            raise ValueError(f"counted.json mixes gene calls from catalogs {sorted(versions)}")
        version = versions.pop() if versions else ""
        self.calls_version = version or None
        if self.calls_version is not None:
            on_disk = (self.root / "catalog" / "version.txt").read_text().strip() + CALLS_SHAPE
            if on_disk != self.calls_version:
                raise ValueError(f"counted.json calls version {self.calls_version} differs from the catalog's {on_disk}")
        structure = json.loads((sp / "structure.json").read_text())
        if structure["genomes"] != self.genomes:
            raise ValueError("structure.json genomes differ from the matrix")
        self.threshold = structure["lineages"]["threshold"]
        self.labels = structure["lineages"]["labels"]
        datasets = json.loads((self.root / "curated" / "datasets.json").read_text())
        self.datasets = {d["id"]: d for d in datasets if d["species_slug"] == slug}
        species = next(iter(self.datasets.values()))["species"] if self.datasets else ""
        self.phenotypes = []
        with (self.root / "curated" / "phenotypes.jsonl").open() as handle:
            for line in handle:
                line = line.rstrip("\n")
                if not line:
                    continue
                if '"genome_id"' in line:
                    p = json.loads(line)
                    if p["species"] == species and p["genome_id"] in self.column:
                        self.phenotypes.append(p)
        self._digest = None

    def matrix_digest(self) -> int:
        """src/analysis.rs `matrix_digest`."""
        if self._digest is None:
            s = self.summary
            h = Fnv()
            h.u64(s["k"])
            h.u64(len(s["genomes"]))
            for g in s["genomes"]:
                h.str(g)
            digests = s.get("kset_digests", [])
            h.u64(len(digests))
            for d in digests:
                h.str(d)
            h.bytes(struct.pack("<I", s["min_support"]))
            for key in ["n_patterns", "kmers_distinct", "kmers_kept", "kmers_rare", "kmers_core"]:
                h.u64(s[key])
            for name in MATRIX_FILES:
                h.bytes((self.matrix_dir / name).read_bytes())
            self._digest = h.h
        return self._digest

    def cohort(self, spec):
        """src/analysis.rs `cohort`: (rows in column order, MIC-inferred count)."""
        drug = [p for p in self.phenotypes if p["antibiotic"] == spec["antibiotic"]]
        rows = [p for p in drug if p.get("reported") in ("R", "S") and not p["conflict"] and not p["label_from_mic"]]
        rows.sort(key=lambda p: self.column[p["genome_id"]])
        return rows, sum(1 for p in drug if p["label_from_mic"])

    def class_counts(self, spec):
        rows, _ = self.cohort(spec)
        r = sum(p["reported"] == "R" for p in rows)
        return r, len(rows) - r

    def fingerprint(self, dataset_id, eval_config) -> str:
        """src/analysis.rs `input_fingerprint` for one curated dataset."""
        spec = self.datasets[dataset_id]
        h = Fnv()
        h.u64(self.matrix_digest())
        h.str(spec["antibiotic"])
        h.u64(spec.get("min_class", 25))
        bp = spec.get("breakpoint")
        if bp is None:
            h.bytes(b"\0")
        else:
            h.bytes(b"\1")
            h.bytes(struct.pack("<i", bp["r_above"]))
            h.bytes(f64_bits(bp["agreement"]))
            h.u64(bp["n"])
        h.bytes(f32_bits(self.threshold))
        if self.calls_version is None:
            h.bytes(b"\0")
        else:
            h.bytes(b"\1")
            h.str(self.calls_version)
        c, f, g = eval_config, eval_config["forest"], eval_config["gbt"]
        for n in [c["folds"], c["inner_folds"], c["scm_max_rules"], f["n_trees"], g["rounds"], g["patience"],
                  c["screened_features"], c["lasso_lambdas"]]:
            h.u64(n)
        for n in [c["cart_max_depth"], f["max_depth"], f["min_leaf"], g["max_depth"]]:
            h.bytes(struct.pack("<I", n))
        for x in [f["subsample"], g["learning_rate"], g["lambda"], g["min_child_weight"], g["colsample"]]:
            h.bytes(f64_bits(x))
        for n in [c["seed"], FOREST_SEED, GBT_SEED]:
            h.u64(n)
        rows, from_mic = self.cohort(spec)
        h.u64(from_mic)
        for p in rows:
            h.str(p["genome_id"])
            h.bytes(b"\1" if p["reported"] == "R" else b"\0")
            mic = p.get("mic")
            if mic is None:
                h.bytes(b"\0")
            else:
                h.bytes(bytes([CENSOR[mic["censor"]]]))
                h.bytes(struct.pack("<i", mic["log2"]))
            h.bytes(struct.pack("<I", self.labels[self.column[p["genome_id"]]]))
            hits = self.hits.get(p["genome_id"], [])
            h.u64(len(hits))
            for hit in hits:
                for key in ["family", "group", "allele", "class", "subclass"]:
                    h.str(hit.get(key, "") if key == "group" else hit[key])
                h.bytes(f32_bits(hit["containment"]))
        return h.hex()


def check_reports(data_root: Path, slug: str):
    inputs = Inputs(data_root, slug)
    rows = []
    for path in sorted((Path(data_root) / "species" / slug / "results").glob("*.json")):
        report = json.loads(path.read_text())
        want = inputs.fingerprint(report["id"], report["eval_config"])
        rows.append({"id": report["id"], "report": report.get("input_fingerprint"), "recomputed": want,
                     "match": report.get("input_fingerprint") == want})
    return {"data_root": str(data_root), "matrix_digest": f"{inputs.matrix_digest():016x}", "calls_version": inputs.calls_version,
            "phenotypes": len(inputs.phenotypes), "reports": rows, "matched": sum(r["match"] for r in rows)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("data_root", type=Path)
    ap.add_argument("--slug", default="mycobacterium_tuberculosis")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    try:
        result = check_reports(args.data_root, args.slug)
    except (OSError, KeyError, ValueError) as err:
        print(f"fingerprint: STOP: {err}", file=sys.stderr)
        return 2
    text = json.dumps(result, indent=1) + "\n"
    if args.out:
        with args.out.open("x") as f:
            f.write(text)
    print(json.dumps({k: v for k, v in result.items() if k != "reports"}, indent=1))
    if args.check and result["matched"] != len(result["reports"]):
        print(f"fingerprint: {len(result['reports']) - result['matched']} reports differ", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
