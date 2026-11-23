"""The BioProjects of a wider arm's added genomes, from BV-BRC's genome table (docs/cap-arms.md, ruling 4).

  python3 -I bioprojects.py <cohorts dir> <master genome meta .jsonl.gz> <answers .jsonl.gz> <counts .json>

Asks BV-BRC's genome table, by genome id, for two fields: genome_id and
bioproject_accession. The ids asked are the cap-all cohort's genomes that
the default's cohort does not hold, the added genomes; nothing else is
asked, no other table or field is read, and no genome is drawn: an answer
for an id not asked stops the fetch. The query is the one `src/bvbrc.rs`
posts to the genome collection (`genome_query`), with its fields cut to
the two:

  POST https://www.bv-brc.org/api/genome/
  Content-Type: application/rqlquery+x-www-form-urlencoded
  in(genome_id,(<ids>))&select(genome_id,bioproject_accession)&limit(<n>)

Writes <answers .jsonl.gz>, one row per answered genome in the master
genome rows' shape, {"genome_id", "bioproject_accession"}, which stays on
the box beside the cohort lists, and <counts .json>, counts only: the
endpoint, the table, the fields, the fetch's UTC date and time, the
batches, the genomes asked, answered, with a BioProject and without, not
answered, and, for the added genomes the master genome rows also name,
how many name the same BioProject and how many another. study.py and
check_folds.py read both files; a genome whose two files name different
BioProjects counts as unknown there. Prints the counts.
"""

import datetime
import gzip
import json
import sys
import time
import urllib.request
from pathlib import Path

API = "https://www.bv-brc.org/api"
TABLE = "genome"
FIELDS = ("genome_id", "bioproject_accession")
BATCH = 500
TRIES = 4


def read_ids(path):
    with open(path, encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def query(ids):
    """(url, body) of one batch: the genome table, by id, the two fields."""
    body = f"in(genome_id,({','.join(ids)}))&select(genome_id,bioproject_accession)&limit({len(ids)})"
    return f"{API}/genome/", body


def check_query(url, body):
    """Raises unless the request asks the genome table for genome_id and bioproject_accession only."""
    if url != f"{API}/{TABLE}/":
        raise ValueError("the query asks another table than genome")
    select = [part for part in body.split("&") if part.startswith("select(")]
    if len(select) != 1 or tuple(select[0][len("select("):-1].split(",")) != FIELDS:
        raise ValueError("the query asks other fields than genome_id and bioproject_accession")


def post(url, body):
    req = urllib.request.Request(url, data=body.encode("utf-8"), method="POST",
                                 headers={"Content-Type": "application/rqlquery+x-www-form-urlencoded",
                                          "Accept": "application/json", "User-Agent": "grm"})
    with urllib.request.urlopen(req, timeout=120) as resp:
        return json.loads(resp.read().decode("utf-8"))


def fetch(ids, send=post, sleep=time.sleep):
    """{genome: BioProject or None} for the answered ids, and the batch count."""
    answers, batches = {}, 0
    asked = set(ids)
    for i in range(0, len(ids), BATCH):
        url, body = query(ids[i:i + BATCH])
        check_query(url, body)
        for attempt in range(TRIES):
            try:
                rows = send(url, body)
                break
            except OSError:
                if attempt == TRIES - 1:
                    raise
                sleep(2 ** attempt)
        batches += 1
        for r in rows:
            if not set(r) <= set(FIELDS):
                raise ValueError("an answer holds a field that was not asked")
            g = r.get("genome_id")
            if g not in asked:
                raise ValueError("an answer names a genome that was not asked")
            if g in answers:
                raise ValueError("a genome is answered twice")
            answers[g] = r.get("bioproject_accession") or None
    return answers, batches


def master(path):
    out = {}
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            if line.strip():
                r = json.loads(line)
                if r.get("bioproject_accession"):
                    out[r["genome_id"]] = r["bioproject_accession"]
    return out


def counts(ids, answers, known, batches, when):
    both = [g for g in ids if answers.get(g) and g in known]
    return {"endpoint": f"{API}/{TABLE}/", "table": TABLE, "fields": list(FIELDS), "fetched_utc": when,
            "batches": batches, "asked": len(ids), "answered": len(answers),
            "with_bioproject": sum(1 for v in answers.values() if v), "without_bioproject": sum(1 for v in answers.values() if not v),
            "not_answered": len(ids) - len(answers),
            "master_same": sum(1 for g in both if known[g] == answers[g]),
            "master_other": sum(1 for g in both if known[g] != answers[g])}


def main(argv, send=post):
    if len(argv) != 5:
        raise SystemExit(__doc__)
    cohorts, meta, out_answers, out_counts = (Path(a) for a in argv[1:])
    base = set(read_ids(cohorts / "cap-2000.txt"))
    ids = sorted(set(read_ids(cohorts / "cap-all.txt")) - base)
    when = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    try:
        answers, batches = fetch(ids, send)
    except ValueError as e:
        raise SystemExit(f"bioprojects: {e}") from None
    with gzip.open(out_answers, "wt", encoding="utf-8") as f:
        for g in sorted(answers):
            if answers[g]:
                f.write(json.dumps({"genome_id": g, "bioproject_accession": answers[g]}) + "\n")
    doc = counts(ids, answers, master(meta), batches, when)
    out_counts.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(doc, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
