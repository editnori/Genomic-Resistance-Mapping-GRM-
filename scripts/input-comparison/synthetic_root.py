#!/usr/bin/env python3
"""A synthetic GRM data root for the input comparison's end-to-end test.

Writes what `grm curate` and the counting stages read, with no real data and
no network: BV-BRC-shaped AMR and genome records, one gzipped FASTA per
genome, empty BV-BRC feature caches (so annotation stays offline), an
AMRFinderPlus-shaped catalog, and truth.json naming every planted effect.
Genome records carry a BV-BRC-shaped `mlst` field (the clade's sequence
type), missing on a share of genomes, which `grm curate --selection
group-round-robin-v1` then draws as singleton groups.

The cohort is one Salmonella enterica x meropenem dataset:
- genomes sit in clades: a shared random core, mutated per clade and lightly
  per genome, plus a clade accessory segment; a few genomes are exact copies
  of another genome of their clade (k-mer-identical duplicates);
- every genome carries the gyrA coding sequence H's salmonella_qrdr assay
  probes, as its own contig; a random share, independent of the phenotype,
  carries a QRDR substitution (the null effect C sees);
- resistant genomes carry one allele of a carbapenemase family, blaKPC, as
  their own contig, so no core k-mer spans an insertion site. The family has
  many unrelated alleles, each in so few genomes that its k-mers fall below a
  dataset's per-fit support floor (max(n / 200, 3)): K cannot see the
  resistance, while H's panels and P's group presence call the family (the
  planted effect);
- a small share of labels is flipped.

With --second-drug the root also holds a Salmonella enterica x ciprofloxacin dataset (for `run.py combine`'s
end-to-end test): every genome gets a ciprofloxacin record, resistant when it carries the QRDR substitution, with
the same share of labels flipped from a separate random stream. Without it the root is byte for byte the same.

`--visible V` grades the planted effect for the power curve (main's ruling after the first calibration
failed): each blaKPC carrier stays resistant, but its allele contig is written only with probability V,
drawn from its own random stream, so H and P see a share V of the resistance. At V = 1, the default, the
root is byte-identical to one written without the option.

usage: python3 -I -B synthetic_root.py OUT_DIR [--genomes N] [--clades N] [--resistant F]
           [--alleles N] [--seed N] [--second-drug] [--visible V]
"""
import argparse
import gzip
import json
import os
import random
import sys

SPECIES = "Salmonella enterica"
DRUG = "meropenem"
SECOND_DRUG = "ciprofloxacin"
TAXON = "28901"
# H's salmonella_qrdr wild type (src/probes.rs, SALMONELLA_GYRA): codons 83-87 are probed.
GYRA = ("ATGAGCGACCTTGCGAGAGAAATTACACCGGTCAACATTGAGGAGGAGCTGAAGAGCTCCTATCTGGATTATGCGATGTCGGTCATTGTTGGCCGTGCGC"
        "TGCCGGATGTCCGAGATGGCCTGAAGCCGGTACACCGTCGCGTACTTTACGCCATGAACGTATTGGGCAATGACTGGAACAAAGCCTATAAAAAATCTGCC"
        "CGTGTCGTTGGTGACGTAATCGGTAAATACCATCCCCACGGCGATTCCGCAGTGTATGACACCATCGTTCGTATGGCGCAGCCATTCTCGC")
QRDR_BASE = 3 * (83 - 1) + 1  # codon 83's second base, 0-based
FAMILY = "blaKPC"
CATALOG_VERSION = "2026-01-01.1"
DEFAULTS = {"genomes": 3000, "clades": 40, "resistant": 0.2, "alleles": 300, "core": 20000,
            "accessory": 1500, "allele_length": 900, "clade_rate": 0.01, "genome_rate": 0.0005,
            "duplicates": 0.02, "qrdr": 0.3, "flip": 0.02, "no_mlst": 0.1}


def bases(rng, n):
    return bytearray(rng.choice(b"ACGT") for _ in range(n))


def mutate(rng, seq, rate):
    out = bytearray(seq)
    for i in rng.sample(range(len(out)), round(len(out) * rate)):
        out[i] = rng.choice([b for b in b"ACGT" if b != out[i]])
    return out


def write_gz(path, text):
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as f:
        f.write(text.encode("ascii"))


def generate(out, cfg, seed):
    rng = random.Random(seed)
    for d in ("bvbrc/features", "fasta", "catalog"):
        os.makedirs(os.path.join(out, d), exist_ok=True)
    core = bases(rng, cfg["core"])
    clades = [(mutate(rng, core, cfg["clade_rate"]), bases(rng, cfg["accessory"])) for _ in range(cfg["clades"])]
    alleles = [b"ATG" + bases(rng, cfg["allele_length"] - 3) for _ in range(cfg["alleles"])]
    genomes, amr, meta = [], [], []
    no_mlst = random.Random(f"mlst:{seed}")
    visible = cfg.get("visible", 1.0)
    shown = random.Random(f"visible:{seed}")
    for i in range(cfg["genomes"]):
        gid = f"{TAXON}.{100000 + i}"
        previous = genomes[-1] if genomes else None
        if previous and rng.random() < cfg["duplicates"]:
            g = dict(previous, genome_id=gid, duplicate_of=previous["genome_id"])
        else:
            clade = rng.randrange(cfg["clades"])
            backbone, accessory = clades[clade]
            carrier = rng.random() < cfg["resistant"]
            g = {"genome_id": gid, "clade": clade, "duplicate_of": None,
                 "chromosome": (mutate(rng, backbone, cfg["genome_rate"]) + accessory).decode(),
                 "qrdr_mutant": rng.random() < cfg["qrdr"],
                 "allele": rng.randrange(cfg["alleles"]) if carrier else None,
                 "flipped": rng.random() < cfg["flip"]}
        genomes.append(g)
    if visible < 1:
        # Each carrier's contig is shown with probability `visible`; a copy shows what its original shows.
        for g in genomes:
            hide = (g["allele"] is not None and shown.random() >= visible) if g["duplicate_of"] is None else None
            g["hidden"] = hide if hide is not None else next(o["hidden"] for o in genomes
                                                             if o["genome_id"] == g["duplicate_of"])
    for g in genomes:
        gyra = bytearray(GYRA.encode())
        if g["qrdr_mutant"]:
            gyra[QRDR_BASE] = b"T"[0] if gyra[QRDR_BASE] != b"T"[0] else b"G"[0]
        contigs = [g["chromosome"], gyra.decode()]
        if g["allele"] is not None and not g.get("hidden"):
            contigs.append(alleles[g["allele"]].decode())
        fasta = "".join(f">{g['genome_id']}.con.{j + 1:04d}\n{seq}\n" for j, seq in enumerate(contigs))
        write_gz(os.path.join(out, "fasta", f"{g['genome_id']}.fna.gz"), fasta)
        with open(os.path.join(out, "bvbrc/features", f"{g['genome_id']}.json"), "w", encoding="ascii") as f:
            f.write("[]")
        resistant = (g["allele"] is not None) != g["flipped"]
        amr.append({"genome_id": g["genome_id"], "antibiotic": MISSPELLED if len(amr) % 150 == 7 else DRUG,
                    "resistant_phenotype": "Resistant" if resistant else "Susceptible",
                    "laboratory_typing_method": "Broth dilution"})
        record = {"genome_id": g["genome_id"], "genome_name": f"{SPECIES} synthetic {g['genome_id']}",
                  "species": SPECIES, "taxon_id": TAXON,
                  "genome_length": str(sum(len(c) for c in contigs)), "contigs": str(len(contigs))}
        g["mlst"] = None if no_mlst.random() < cfg["no_mlst"] else f"MLST.Salmonella_enterica.{g['clade'] + 1}"
        if g["mlst"]:
            record["mlst"] = g["mlst"]
        meta.append(record)
        g["resistant"] = resistant
    if cfg.get("second_drug"):
        flips = random.Random(f"second:{seed}")
        for g in genomes:
            g["resistant_second"] = g["qrdr_mutant"] != (flips.random() < cfg["flip"])
            amr.append({"genome_id": g["genome_id"], "antibiotic": cfg["second_drug"],
                        "resistant_phenotype": "Resistant" if g["resistant_second"] else "Susceptible",
                        "laboratory_typing_method": "Broth dilution"})
    with open(os.path.join(out, "bvbrc/genome_amr.jsonl"), "w", encoding="ascii") as f:
        f.writelines(json.dumps(r, sort_keys=True) + "\n" for r in amr)
    with open(os.path.join(out, "bvbrc/genomes.jsonl"), "w", encoding="ascii") as f:
        f.writelines(json.dumps(r, sort_keys=True) + "\n" for r in meta)
    cat = os.path.join(out, "catalog")
    with open(os.path.join(cat, "AMR_CDS.fa"), "w", encoding="ascii") as f:
        for a, seq in enumerate(alleles, 1):
            f.write(f">WP_9{a:05d}.1|NG_9{a:05d}.1|1|1|{FAMILY}-{a}|{FAMILY}|carbapenem-hydrolyzing class A "
                    f"beta-lactamase KPC-{a}\n{seq.decode()}\n")
    with open(os.path.join(cat, "ReferenceGeneCatalog.txt"), "w", encoding="ascii") as f:
        f.write("allele\tgene_family\tproduct_name\ttype\tclass\tsubclass\trefseq_protein_accession\t"
                "genbank_protein_accession\n")
        for a in range(1, len(alleles) + 1):
            f.write(f"{FAMILY}-{a}\t{FAMILY}\tcarbapenem-hydrolyzing class A beta-lactamase KPC-{a}\tAMR\t"
                    f"BETA-LACTAM\tCARBAPENEM\tWP_9{a:05d}.1\t\n")
    with open(os.path.join(cat, "version.txt"), "w", encoding="ascii") as f:
        f.write(CATALOG_VERSION + "\n")
    truth = {
        "format": "grm-input-comparison-synthetic", "version": 1, "seed": seed, "config": cfg,
        "species": SPECIES, "drug": DRUG,
        "planted": {"effect": "resistance is carriage of any blaKPC allele; each allele is too rare for K's "
                              "per-fit support floor, while H and P call the family",
                    "null": "QRDR substitutions (block C) are drawn independently of the phenotype"},
        "genomes": [{k: g[k] for k in ("genome_id", "clade", "duplicate_of", "qrdr_mutant", "allele",
                                       "flipped", "resistant", "mlst", "resistant_second", "hidden") if k in g}
                    for g in genomes],
    }
    if visible < 1:
        truth["planted"]["visible"] = (f"a share {visible} of carriers show their blaKPC contig; the rest are "
                                       "resistant with nothing any block sees")
    curate_inputs(out, genomes)
    truth["curate_inputs"] = {"species_check_fails": [FAILED_GENOME(genomes)],
                              "misspelled": [g["genome_id"] for i, g in enumerate(genomes) if i % 150 == 7],
                              "fold": [MISSPELLED, DRUG]}
    with open(os.path.join(out, "truth.json"), "w", encoding="ascii") as f:
        json.dump(truth, f, sort_keys=True, indent=1)
    return truth


# What revision 12's curate reads beyond BV-BRC (main's ruling on REVIEW-INPUT7 M1): a species check that fails
# one genome the default would select (the lowest-id resistant genome), and a SigDrug-format vocabulary whose
# brand fold turns "synthmerrem", the name every 150th genome's record uses, into meropenem. GRM's own resolver
# can't reach it (its typo match stops at edit distance 2), so only the vocabulary labels those genomes.
MISSPELLED = "synthmerrem"
VOCAB_FILES = {
    "preferred_ingredients.tsv": "value\tclassification\tpreferred_name\n",
    "preferred_base_names.tsv": "value\tclassification\tpreferred_name\n",
    "drug_name_typo_folds.tsv": "observed\tcanonical\treceipt\tdate\tscope_note\n",
    "brand_generic_folds.tsv": "observed\tcanonical\treceipt\tdate\tscope_note\n",
    "ingredient_abbrev_folds.tsv": "observed\tcanonical\treceipt\tdate\tscope_note\n",
    "salt_base_ingredients.tsv": "ingredient_name\tsalt_base_ingredient\tsalt_base_seqn\tsalt_token\tsource\n",
}


def FAILED_GENOME(genomes):
    return min(g["genome_id"] for g in genomes if g["resistant"])


def curate_inputs(out, genomes):
    failed = FAILED_GENOME(genomes)
    os.makedirs(os.path.join(out, "species-check"))
    check = {"k": 21, "sketch_size": 1000,
             "config": {"fail_distance": 0.13, "species_distance": 0.05, "gc_z": 6.0, "gc_points": 1.0,
                        "genomic_groups": []},
             "species": [{"species": SPECIES, "slug": "salmonella_enterica", "reference": genomes[0]["genome_id"],
                          "gc_median": 0.5, "gc_sd": 0.0005, "checked": 2, "failed": 1, "borderline": 0,
                          "unchecked": [], "distance_p50": 0.01, "distance_p99": 0.02, "distance_max_pass": 0.02}],
             "pairs": [],
             "genomes": [{"genome_id": genomes[0]["genome_id"], "species": SPECIES, "verdict": "pass", "gc": 0.5,
                          "gc_z": 0.0, "distance": 0.01, "ani": 0.99, "nearest_species": None,
                          "nearest_distance": None, "inside_species": None, "reason": None},
                         {"genome_id": failed, "species": SPECIES, "verdict": "fail", "gc": 0.5, "gc_z": 0.0,
                          "distance": 0.2, "ani": 0.8, "nearest_species": "Salmonella bongori",
                          "nearest_distance": 0.03, "inside_species": "Salmonella bongori",
                          "reason": "synthetic: nearer another species' reference"}]}
    with open(os.path.join(out, "species-check", "species_check.json"), "w", encoding="ascii") as f:
        json.dump(check, f, sort_keys=True, indent=1)
    os.makedirs(os.path.join(out, "vocabulary"))
    for name, header in VOCAB_FILES.items():
        with open(os.path.join(out, "vocabulary", name), "w", encoding="ascii") as f:
            f.write(header)
            if name == "brand_generic_folds.tsv":
                f.write(f"{MISSPELLED}\t{DRUG}\tsynthetic\t2026-11-06\tsynthetic brand fold\n")


def main(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("out")
    for key in ("genomes", "clades", "alleles"):
        p.add_argument(f"--{key}", type=int, default=DEFAULTS[key])
    p.add_argument("--resistant", type=float, default=DEFAULTS["resistant"])
    p.add_argument("--seed", type=int, default=11)
    p.add_argument("--second-drug", action="store_true", help=f"add a {SECOND_DRUG} dataset (QRDR carriers resistant)")
    p.add_argument("--visible", type=float, default=1.0,
                   help="the share of blaKPC carriers whose allele contig is written (the power curve's size)")
    a = p.parse_args(argv)
    if not 0 < a.visible <= 1:
        sys.exit(f"--visible {a.visible} is not in (0, 1]")
    if os.path.exists(a.out) and os.listdir(a.out):
        sys.exit(f"{a.out} is not empty")
    cfg = dict(DEFAULTS, genomes=a.genomes, clades=a.clades, alleles=a.alleles, resistant=a.resistant)
    if a.second_drug:
        cfg["second_drug"] = SECOND_DRUG
    if a.visible < 1:
        cfg["visible"] = a.visible
    truth = generate(a.out, cfg, a.seed)
    n_r = sum(g["resistant"] for g in truth["genomes"])
    print(f"{len(truth['genomes'])} genomes ({n_r} resistant) in {cfg['clades']} clades written to {a.out}")


if __name__ == "__main__":
    main(sys.argv[1:])
