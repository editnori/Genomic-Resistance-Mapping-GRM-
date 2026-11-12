"""Derive one stage's input manifest from the run's master manifest (lane DATA's, pinned by MASTER_SHA256): every
line the stage reads, in the master's own bytes, so each line's hash is the master's. Taken: the shared inputs
(bvbrc/genome_amr.jsonl, bvbrc/genomes.jsonl, curated/* but grm curate's provenance side file, catalog/*,
species-check/*), and for each species named,
species/<slug>/genomes.txt, then fasta/<id>.fna.gz and bvbrc/features/<id>.json for each id in it (read from
DATA_ROOT, whose genomes.txt must hash as the master says). With CANARY=1 (the name is kept; it now feeds the
per-stage revision-9 diff report, not a gate), CANARY_PREFIX's index.json and the named species' reports under it.
A species without genomes.txt, or a selected genome without its FASTA or annotation line, is
an error (exit 1, message on the last line): the run would fetch it from BV-BRC.

usage: python3 -I stage-manifest.py DATA_ROOT MASTER SPECIES_CSV CANARY_PREFIX CANARY(0|1) > STAGE.sha256"""
import hashlib
import os
import sys

root, master, species, prefix, canary = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5] == "1"
# EXTRA (optional 6th argument): comma list of further genome lists in the master, such as MTB's
# bindings/full-genomes.txt, whose ids' FASTA and annotations the stage also reads
extra = [x for x in (sys.argv[6] if len(sys.argv) > 6 else "").split(",") if x]
# grm curate's provenance side file (GRM docs/provenance-fields.md, about 310 MB on the real root): no stage reads
# it, so a master manifest that lists it never sends it to a box.
SIDE_FILES = ("curated/provenance.jsonl",)
lines = {}
for line in open(master, encoding="utf-8").read().splitlines():
    lines[line[66:]] = line
out, errors = [], []


def take(p):
    if p not in lines:
        errors.append(p)
    else:
        out.append(lines[p])


for p in sorted(lines):
    if p in SIDE_FILES:
        continue
    if p in ("bvbrc/genome_amr.jsonl", "bvbrc/genomes.jsonl") or p.split("/")[0] in ("curated", "catalog", "species-check") \
            or (extra and p.split("/")[0] == "bindings"):
        out.append(lines[p])
for s in [x for x in species.split(",") if x] + ["@" + x for x in extra]:
    g = s[1:] if s.startswith("@") else "species/%s/genomes.txt" % s
    take(g)
    if g not in lines:
        continue
    data = open(os.path.join(root, g), "rb").read()
    if hashlib.sha256(data).hexdigest() != lines[g][:64]:
        errors.append(g + " differs from the master manifest")
        continue
    for i in data.decode().split():
        take("fasta/%s.fna.gz" % i)
        take("bvbrc/features/%s.json" % i)
if canary:
    # the revision-9 reports of this stage's species only (datasets/<slug>__<drug>.json) and the index, for the
    # per-stage diff report; a species with none (new since revision 9) just has none
    root9 = prefix.rstrip("/") + "/"
    want = tuple("%sdatasets/%s__" % (root9, s) for s in species.split(",") if s)
    tree = [p for p in sorted(lines) if p == root9 + "index.json" or (want and p.startswith(want))]
    if root9 + "index.json" not in lines:
        errors.append("no revision-9 index under %s" % prefix)
    out += [lines[p] for p in tree]
if errors:
    print("%d stage inputs missing from the master manifest, first: %s" % (len(errors), errors[0]))
    sys.exit(1)
seen = set()
for line in out:
    if line not in seen:
        seen.add(line)
        print(line)
