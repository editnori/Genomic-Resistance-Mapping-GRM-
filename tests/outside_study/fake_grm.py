"""A stand-in for the grm binary in the run script's tests: `fetch fasta`,
`species-check` and `count`, each writing the files the real command writes
in the shape the run script reads, and nothing else. Behaviour is set by
environment variables:

  FAKE_MISSING   comma-separated ids fetch fails to get
  FAKE_VERDICTS  JSON {id: verdict} for species-check (default pass)
  FAKE_UNCOUNTED comma-separated ids count fails on
"""

import json
import os
import sys
from pathlib import Path


def ids_of(root, slug):
    p = root / "species" / slug / "genomes.txt"
    return [x.strip() for x in p.read_text().splitlines() if x.strip()] if p.exists() else []


def main(argv):
    root = Path(argv[argv.index("--data") + 1])
    cmd = argv[argv.index("--data") + 2:]
    env = lambda k: [x for x in os.environ.get(k, "").split(",") if x]
    if cmd[:2] == ["fetch", "fasta"]:
        slug = cmd[cmd.index("--species") + 1]
        missing = env("FAKE_MISSING")
        fails = []
        for g in ids_of(root, slug):
            if g in missing:
                fails.append(f"{g}\tno contigs returned\n")
                continue
            dst = root / "fasta" / f"{g}.fna.gz"
            if not dst.exists():
                dst.write_bytes(b"fasta " + g.encode())
        if fails:
            (root / "species" / slug / "fasta_failures.tsv").write_text("".join(fails))
        print(f"fetched; {len(fails)} failed")
        return 0
    if cmd[:1] == ["species-check"]:
        verdicts = json.loads(os.environ.get("FAKE_VERDICTS", "{}"))
        path = root / "species-check" / "species_check.json"
        doc = json.loads(path.read_text())
        rows = []
        for sp in json.loads((root / "curated" / "species.json").read_text()):
            for g in ids_of(root, sp["slug"]):
                v = verdicts.get(g, "pass")
                rows.append({"genome_id": g, "species": sp["species"], "verdict": v, "gc": 0.33,
                             "distance": 0.01, "reason": "far from the reference" if v == "fail" else None})
        doc["genomes"] = rows
        path.write_text(json.dumps(doc))
        print(f"judged {len(rows)}")
        return 0
    if cmd[:1] == ["count"]:
        uncounted = env("FAKE_UNCOUNTED")
        version = (root / "catalog" / "version.txt")
        cv = (version.read_text().strip() if version.exists() else "none") + "|group-calls-v2"
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
                    (d / "ksets" / f"{g}.kset").write_bytes(b"kset " + g.encode())
                    out["genomes"].append({"genome_id": g, "catalog_version": cv, "genes": []})
            (d / "counted.json").write_text(json.dumps(out))
        print("counted")
        return 0
    print(f"fake grm: unsupported {cmd}", file=sys.stderr)
    return 64


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
