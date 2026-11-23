"""Checks the fetch of one species' cap-all genomes (docs/cap-arms.md, The fetch).

  python3 -I fetch_check.py <release species json> <cohorts dir> <fasta dir> <master list>
                            <master genome meta .jsonl.gz> <fasta_failures.tsv> <out json> <fetched sha256>

After `grm fetch fasta` over <cohorts dir>/cap-all.txt. A genome failed when
grm fetch fasta wrote it in fasta_failures.tsv or wrote no file for it.

- A release genome (the default's) must be fetched, with the master list's
  sha256 and the release's length, contig count and N50; any failure or
  difference exits 1.
- An added genome in the master genome rows is held to the master list's
  sha256, the count that differs written; any other added genome has its
  sha256 written to <fetched sha256>, which stays on the box.
- A failed added genome leaves every wider cohort list: cap-*.txt are
  written again without it, and cohorts.json gains its count.

Each genome's fetch time is its file's birth to its last write (GNU stat,
%W and %Y). Writes <out json>, counts only: genomes, failures by kind, the
sha256 and statistics checks, bytes, and the seconds per genome (count,
sum, median, 90th centile, largest). Prints it.
"""

import gzip
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def stats(path):
    """(length, contigs, N50) of a gzip FASTA, as the release counts them."""
    lens = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                lens.append(0)
            elif lens:
                lens[-1] += len(line.strip())
    tot, acc, n50 = sum(lens), 0, -1
    for x in sorted(lens, reverse=True):
        acc += x
        if 2 * acc >= tot:
            n50 = x
            break
    return tot, len(lens), n50


def read_lines(path):
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def master_list(path):
    out = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            digest, _, name = line.rstrip("\n").partition("  ")
            out[name] = digest
    return out


def meta_ids(path):
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return {json.loads(line)["genome_id"] for line in f if line.strip()}


def seconds(paths):
    """Per file, its last write less its birth, in seconds; files whose birth the file system does not record
    are left out."""
    out = []
    paths = list(paths)
    for i in range(0, len(paths), 500):
        res = subprocess.run(["stat", "-c", "%.3W %.3Y"] + [str(p) for p in paths[i:i + 500]],
                             check=True, capture_output=True, text=True)
        for line in res.stdout.splitlines():
            born, wrote = (float(x) for x in line.split())
            if born > 0:
                out.append(max(0.0, wrote - born))
    return out


def centile(xs, q):
    xs = sorted(xs)
    return xs[min(len(xs) - 1, int(q * (len(xs) - 1) + 0.5))] if xs else None


def drop_failed(cohorts, failed):
    """Writes each wider cohort list again without the failed genomes; returns how many left each."""
    out = {}
    for p in sorted(Path(cohorts).glob("cap-*.txt")):
        if p.stem == "cap-2000":
            continue
        ids = read_lines(p)
        keep = [g for g in ids if g not in failed]
        out[p.stem] = len(ids) - len(keep)
        p.write_text("".join(g + "\n" for g in keep), encoding="utf-8")
    doc_path = Path(cohorts) / "cohorts.json"
    doc = json.loads(doc_path.read_text(encoding="utf-8"))
    for arm, n in out.items():
        a = doc["arms"][arm]
        a["fetch_failed"] = n
        a["genomes"] -= n
        a["added"] -= n
    doc_path.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    return out


def main(argv):
    if len(argv) != 9:
        raise SystemExit(__doc__)
    release_json, cohorts, fasta, master, meta, failures, out_json, fetched = argv[1:]
    with open(release_json, encoding="utf-8") as f:
        rel = {g["id"]: g for g in json.load(f)["genomes"]}
    ids = read_lines(Path(cohorts) / "cap-all.txt")
    if not set(rel) <= set(ids):
        raise SystemExit("the cap-all list does not hold the release's whole list")
    listed_failed = {line.split("\t")[0] for line in read_lines(failures)} if os.path.isfile(failures) else set()
    path = {g: Path(fasta) / f"{g}.fna.gz" for g in ids}
    failed = {g for g in ids if g in listed_failed or not path[g].is_file()}
    sums = master_list(master)
    in_meta = meta_ids(meta)
    doc = {"genomes": len(ids), "release": len(rel), "added": len(ids) - len(rel),
           "failed_release": sum(g in failed for g in rel), "failed_added": sum(g in failed and g not in rel for g in ids)}
    rel_bytes = rel_stats = 0
    for g in rel:
        if g in failed:
            continue
        rel_bytes += sums.get(f"fasta/{g}.fna.gz") != sha256(path[g])
        rel_stats += stats(path[g]) != (rel[g]["length"], rel[g]["contigs"], rel[g]["n50"])
    doc["release_bytes_differ_from_master"] = rel_bytes
    doc["release_stats_differ_from_release"] = rel_stats
    added = [g for g in ids if g not in rel and g not in failed]
    in_master = [g for g in added if g in in_meta and f"fasta/{g}.fna.gz" in sums]
    doc["added_in_master"] = len(in_master)
    doc["added_bytes_differ_from_master"] = sum(sums[f"fasta/{g}.fna.gz"] != sha256(path[g]) for g in in_master)
    recorded = [g for g in added if g not in set(in_master)]
    with open(fetched, "w", encoding="utf-8") as f:
        for g in recorded:
            f.write(f"{sha256(path[g])}  fasta/{g}.fna.gz\n")
    doc["added_sha256_recorded"] = len(recorded)
    ok = [path[g] for g in ids if g not in failed]
    doc["bytes"] = sum(p.stat().st_size for p in ok)
    secs = seconds(ok)
    doc["seconds_per_genome"] = {"timed": len(secs), "sum": round(sum(secs), 1),
                                 "median": centile(secs, 0.5), "p90": centile(secs, 0.9),
                                 "max": max(secs) if secs else None}
    if doc["failed_added"]:
        doc["left_wider_cohorts"] = drop_failed(cohorts, failed - set(rel))
    Path(out_json).write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(doc, sort_keys=True))
    if doc["failed_release"] or rel_bytes or rel_stats:
        print("a release genome failed or differs from the master list or the release")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
