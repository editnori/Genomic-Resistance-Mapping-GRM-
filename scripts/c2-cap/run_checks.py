"""The checks of plan block C2's run, before any score is read (docs/cap-arms-result.md, Method).

  python3 -I run_checks.py <runs dir> <release dir> <out txt> [<id list> ...]

<runs dir> holds one folder per species, named by its slug, as each run
came back: run_species.sh's OUT layout less work/ (results/, logs/,
expect.tsv, binary.txt, lists.txt, DONE or FAILED, SHA256SUMS), with
score/ as c2-drive.sh brings it back, score.tar beside the rest and not
unpacked (or, where no score.tar is there, a score/ folder). <release dir>
is the release data-rev14-274-20261112, whose species/<slug>.json lists
each species' 2,000 genomes. It reads no score: score.tar's tables are read
in place, only as bytes, for SHA256SUMS and the genome id count; score.tar
is never read whole, and its score/score.log is never extracted, read,
listed or counted.

Checks, each printed with its counts and none with a genome id:

- species: the runs directory holds the species of lists.tsv's first run,
  each once;
- SHA256SUMS: every line of each species' SHA256SUMS holds, and none
  names score/score.log; c2-drive.sh's back.tar leaves lists.txt out, so a
  copy without it takes run_species.sh's printf of the species' first run
  (lists.tsv's order) and c2-drive.sh's arguments (3 draws, 3 twin slots,
  cap-all, lineage,study, yes), and only when SHA256SUMS' lists.txt line
  holds it; and every table run_species.sh lists (a .json, .tsv or .txt
  file under results/ or score/) is listed;
- DONE: each species has DONE naming its species, binary.txt's grm and
  code sha256 and the sha256 of its SHA256SUMS; a FAILED beside a newer
  DONE is named and set aside by time, and a FAILED with no newer DONE
  fails;
- one binary: every species' binary.txt names the same commit and code
  sha256; else each is named with its species, and the check fails; two
  binary sha256s under one code sha256 are each named with their species,
  and the check does not fail on them (the lead's answer 3);
- expect.py: each expect.tsv is expect.py's from the run's lists.txt,
  roots.json, added-labelled.tsv and the study fold files' table
  (study-files-<arm>.tsv: a dataset with a study plan, and its file
  written or not);
- cohorts: each cohorts.json against the selection rule's counts of
  docs/cap-arms.md (ranked, added at 4,000 and at 4,096), each cohort at
  most 4,096 and each wider one holding the default's 2,000; added
  genomes that failed the fetch (fetch.json) are named;
- genome ids and BioProjects: 0 of the release's genome ids of the 8
  species, or of any <id list> given (one id a line, as the cohort lists
  on the box), and 0 BioProject accessions in every file that will be read
  or published (everything above but score/score.log), with a genome id
  read as digits, a dot and digits with no digit or dot before it
  (scrub.py's token lets no letter or underscore stand there, and misses an
  id joined to a word);
- refused lineage draws: the datasets whose lineage split is unavailable,
  from the fold check's logs, and their refused slots, counted.

Writes the printout to <out txt> and exits 1 when a check fails.
"""

import hashlib
import importlib.util
import json
import posixpath
import re
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCORE_LOG = "score/score.log"
SCORE_TAR = "score.tar"
# docs/cap-arms.md, The arms: ranked by the rule, added at 4,000, added at 4,096.
SELECTION = {
    "escherichia_coli": (8704, 2000, 2096),
    "enterococcus_faecium": (2907, 907, 907),
    "klebsiella_pneumoniae": (6159, 2000, 2096),
    "mycobacterium_tuberculosis": (21295, 2000, 2096),
    "neisseria_gonorrhoeae": (5396, 2000, 2096),
    "salmonella_enterica": (8170, 2000, 2096),
    "staphylococcus_aureus": (5045, 2000, 2096),
    "streptococcus_pneumoniae": (15530, 2000, 2096),
}
DEFAULT_GENOMES, GATE = 2000, 4096
# c2-drive.sh's call of run_species.sh: draws, twin slots, arms, splits, score.
LISTS_ARGS = ("3", "3", "cap-all", "lineage,study", "yes")
BIOPROJECT = re.compile(rb"PRJ[A-Z]{2}[0-9]+")
# A genome id as a token: scrub.py's pattern lets no word character stand before an id, so it misses
# genome_1280.1234 and id1280.1234; here only a digit or a dot may not stand before it.
ID_TOKEN = re.compile(rb"(?<![\d.])\d+\.\d+(?!\d)")
# run_species.sh writes DONE's and FAILED's time with date -u +%FT%TZ; any other form is not a time.
TIME = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ")
REFUSED = re.compile(r"^(\S+): lineage split unavailable: ([0-9]+) default slot plans refused")


def _load(name):
    spec = importlib.util.spec_from_file_location(f"c2_{name}_for_checks", HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


expect, scrub = _load("expect"), _load("scrub")


def read_tsv(path):
    with open(path, encoding="utf-8") as f:
        head = f.readline().rstrip("\n").split("\t")
        return [dict(zip(head, line.rstrip("\n").split("\t"))) for line in f if line.strip()]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def first_run(lrows):
    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]


def published(folder):
    """Every file of a species' copy that will be read or published, but score/score.log and score.tar, whose
    tables score_files reads in place."""
    return sorted(p for p in Path(folder).rglob("*") if p.is_file() and p.relative_to(folder).as_posix() != SCORE_LOG
                  and p.relative_to(folder).as_posix() != SCORE_TAR and "work" not in p.relative_to(folder).parts)


def score_files(folder):
    """{name: bytes} of a copy's score/ tables: read in place from score.tar when the copy holds it, else from a
    score/ folder. score/score.log is never extracted, opened or read."""
    folder, out = Path(folder), {}
    tar = folder / SCORE_TAR
    if tar.is_file():
        with tarfile.open(tar) as archive:
            for member in archive:
                name = posixpath.normpath(member.name)
                if name == SCORE_LOG or not member.isfile() or not name.startswith("score/"):
                    continue
                out[name] = archive.extractfile(member).read()
        return out
    top = folder / "score"
    for p in sorted(top.rglob("*")) if top.is_dir() else []:
        name = p.relative_to(folder).as_posix()
        if p.is_file() and name != SCORE_LOG:
            out[name] = p.read_bytes()
    return out


TABLE_SUFFIXES = (".json", ".tsv", ".txt")


def unlisted(folder, scores):
    """The tables run_species.sh would list (a .json, .tsv or .txt file under results/ or score/) that the copy's
    SHA256SUMS does not."""
    folder = Path(folder)
    listed = {posixpath.normpath(l.split("  ", 1)[-1])
              for l in (folder / "SHA256SUMS").read_text(encoding="utf-8").splitlines() if l.strip()}
    have = {p.relative_to(folder).as_posix() for p in (folder / "results").rglob("*")
            if p.is_file() and p.suffix in TABLE_SUFFIXES} if (folder / "results").is_dir() else set()
    have |= {n for n in scores if n.endswith(TABLE_SUFFIXES)}
    return sorted(have - listed)


def lists_text(folder, slug, lrows):
    """(lists.txt's text, rebuilt): the copy's own, or where the copy has none, run_species.sh's printf of the
    species' first run in lists.tsv's order with c2-drive.sh's arguments."""
    p = Path(folder) / "lists.txt"
    if p.is_file():
        return p.read_text(encoding="utf-8"), False
    datasets = [r["dataset"] for r in first_run(lrows) if r["species_slug"] == slug]
    draws, twin, arms, splits, score = LISTS_ARGS
    return (f"species\t{slug}\ndatasets\t{','.join(datasets)}\ndraws\t{draws}\ntwin\t{twin}\narms\t{arms}\n"
            f"splits\t{splits}\nscore\t{score}\n"), True


def check_sums(folder, scores=None, lists=None):
    """(lines, bad lines, names score.log). A score/ line is held to score_files' bytes, and a rebuilt lists.txt
    (lists_text) to its line, which must be there."""
    scores = score_files(folder) if scores is None else scores
    lines = [l for l in (folder / "SHA256SUMS").read_text(encoding="utf-8").splitlines() if l.strip()]
    bad, log = 0, 0
    rebuilt = lists is not None and lists[1]
    if rebuilt and not any(l.split("  ", 1)[-1] == "lists.txt" for l in lines):
        bad += 1
    for line in lines:
        want, name = line.split("  ", 1)
        name = posixpath.normpath(name)
        if name == SCORE_LOG or name.startswith(("/", "../")) or name == "..":
            # Named, never opened or hashed: score.log stays out of the record, and nothing outside the copy is read.
            log += name == SCORE_LOG
            bad += 1
            continue
        if name == "lists.txt" and rebuilt:
            bad += hashlib.sha256(lists[0].encode()).hexdigest() != want
            continue
        if name.startswith("score/"):
            data = scores.get(name)
            bad += data is None or hashlib.sha256(data).hexdigest() != want
            continue
        p = folder / name
        bad += not p.is_file() or sha256(p) != want
    return len(lines), bad, log


def done_of(folder):
    """{species, time, binary, code, sums} of DONE, or None."""
    p = folder / "DONE"
    if not p.is_file():
        return None
    parts = p.read_text(encoding="utf-8").split()
    if len(parts) < 11 or parts[0] != "DONE" or not TIME.fullmatch(parts[2]):
        return {}
    return {"species": parts[1], "time": parts[2], "binary": parts[4], "code": parts[6], "sums": parts[10]}


def failed_time(folder):
    p = folder / "FAILED"
    if not p.is_file():
        return None
    parts = p.read_text(encoding="utf-8").split()
    return parts[-1] if parts and TIME.fullmatch(parts[-1]) else ""


def binary_of(folder):
    out = {}
    for line in (folder / "binary.txt").read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) >= 2 and parts[0] in ("commit", "grm", "code") and parts[0] not in out:
            out[parts[0]] = parts[1]
    return out


def lists_txt(text):
    out = {}
    for line in text.splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            out[k] = v
    return out


def expected(folder, slug, text):
    """expect.py's rows from the run's own counts-only files and lists.txt's text."""
    lt = lists_txt(text)
    res = folder / "results"
    roots = json.loads((res / "roots.json").read_text(encoding="utf-8"))
    study, written = [], set()
    for arm in roots["fits"]:
        p = res / f"study-files-{arm}.tsv"
        if p.is_file():
            for r in read_tsv(p):
                if r["dataset"] not in study:
                    study.append(r["dataset"])
                if r["written"] == "1":
                    written.add((arm, r["dataset"]))
    rows = expect.rows(slug, expect.comma(lt["datasets"]), sorted(study), int(lt["draws"]), int(lt["twin"]),
                       expect.comma(lt["arms"]), expect.comma(lt["splits"]), roots,
                       lambda a, d: (a, d) in written, expect.none_added(res / "added-labelled.tsv"))
    return [dict(zip(expect.COLUMNS, r)) for r in rows]


def check_cohorts(folder, slug, selection):
    """[problem] of one species' cohorts.json against the selection rule's counts."""
    doc = json.loads((folder / "results" / "cohorts.json").read_text(encoding="utf-8"))
    arms, out = doc.get("arms", {}), []
    ranked, at4000, at4096 = selection[slug]
    if doc.get("slug") != slug:
        out.append("another species' file")
    if doc.get("ranked") != ranked:
        out.append(f"ranked {doc.get('ranked')}, not {ranked}")
    want = {"cap-2000": (DEFAULT_GENOMES, 0), "cap-4000": (DEFAULT_GENOMES + at4000, at4000),
            "cap-all": (DEFAULT_GENOMES + at4096, at4096)}
    for arm, (genomes, added) in want.items():
        a = arms.get(arm)
        if a is None:
            out.append(f"no {arm}")
            continue
        failed = int(a.get("fetch_failed", 0))
        if (a["genomes"] + failed, a["added"] + failed) != (genomes, added):
            out.append(f"{arm} holds {a['genomes']} genomes and adds {a['added']}, not {genomes} and {added}")
        if a["genomes"] > GATE:
            out.append(f"{arm} holds {a['genomes']} genomes, past {GATE:,}")
        if a["holds_default"] != DEFAULT_GENOMES or a["lines"] != a["genomes"]:
            out.append(f"{arm} does not hold the default's {DEFAULT_GENOMES:,} or its list's lines differ")
    return out


def release_ids(release, slugs, lists=()):
    """The release's genome ids of the species, and those of each id list, as bytes."""
    ids = set()
    for slug in slugs:
        doc = json.loads((Path(release) / "species" / f"{slug}.json").read_text(encoding="utf-8"))
        ids.update(g["id"] for g in doc["genomes"])
    for p in lists:
        ids.update(line.strip() for line in Path(p).read_text(encoding="utf-8").splitlines() if line.strip())
    return {i.encode() for i in ids}


def scan(files, ids, texts=()):
    """(genome ids found, BioProject accessions found) in the files and in score.tar's tables' bytes."""
    found_ids = found_bp = 0
    for text in [p.read_bytes() for p in files] + list(texts):
        found_ids += sum(1 for t in ID_TOKEN.findall(text) if t in ids)
        found_bp += len(BIOPROJECT.findall(text))
    return found_ids, found_bp


def run(runs, release, lrows, selection=None, lists=()):
    """(printed lines, all passed)."""
    selection = SELECTION if selection is None else selection
    runs = Path(runs)
    rows = first_run(lrows)
    species = sorted({r["species_slug"] for r in rows})
    out, ok = [], True

    def say(name, passed, text):
        nonlocal ok
        ok = ok and passed
        out.append(f"{'pass' if passed else 'FAIL'} {name}: {text}")

    found = sorted(p.name for p in runs.iterdir() if p.is_dir())
    say("species", found == species, f"{len(found)} copies for the {len(species)} species of the first run"
        + ("" if found == species else f"; missing {sorted(set(species) - set(found))}, more {sorted(set(found) - set(species))}"))
    slugs = [s for s in species if s in found]

    scores = {slug: score_files(runs / slug) for slug in slugs}
    texts = {slug: lists_text(runs / slug, slug, lrows) for slug in slugs}
    total = bad = logs = 0
    for slug in slugs:
        n, b, l = check_sums(runs / slug, scores[slug], texts[slug])
        total, bad, logs = total + n, bad + b, logs + l
    rebuilt = [slug for slug in slugs if texts[slug][1]]
    extra = [f"{slug}/{n}" for slug in slugs for n in unlisted(runs / slug, scores[slug])]
    say("SHA256SUMS", bad == 0 and logs == 0 and not extra, f"{total - bad} of {total} lines hold in {len(slugs)} "
        f"species; score/score.log listed {logs} times; {len(extra)} tables not listed"
        + "".join(f"; not listed: {x}" for x in extra[:5])
        + (f"; lists.txt not in {len(rebuilt)} copies, rebuilt from lists.tsv and c2-drive.sh's arguments and held "
           f"to SHA256SUMS ({', '.join(rebuilt)})" if rebuilt else ""))

    done_ok, set_aside, notes = 0, [], []
    for slug in slugs:
        f = runs / slug
        d, b = done_of(f), binary_of(f)
        t_failed = failed_time(f)
        good = (d is not None and d.get("species") == slug and d.get("binary") == b.get("grm")
                and d.get("code") == b.get("code") and d.get("sums") == sha256(f / "SHA256SUMS"))
        if t_failed is not None:
            if good and t_failed and t_failed < d["time"]:
                set_aside.append(f"{slug} FAILED at {t_failed}, DONE at {d['time']}")
            else:
                good = False
                notes.append(f"{slug} FAILED at {t_failed} with no newer DONE" if t_failed else
                             f"{slug} FAILED with no time in date -u +%FT%TZ form")
        done_ok += good
    say("DONE", done_ok == len(slugs), f"{done_ok} of {len(slugs)} name their species, binary, code and SHA256SUMS"
        + "".join(f"; set aside: {s}" for s in set_aside) + "".join(f"; {s}" for s in notes))

    groups = {}
    for slug in slugs:
        b = binary_of(runs / slug)
        groups.setdefault((b.get("commit"), b.get("code")), []).append(slug)
    grms = {}
    for slug in slugs:
        b = binary_of(runs / slug)
        grms.setdefault((b.get("commit"), b.get("code")), {}).setdefault(b.get("grm"), []).append(slug)
    named = "".join(f"; under code {k}, binary {g}: {len(v)} species ({', '.join(v)})"
                    for (c, k), by in sorted(grms.items(), key=lambda x: str(x[0])) if len(by) > 1
                    for g, v in sorted(by.items(), key=lambda x: str(x[0])))
    whole = all(c and k for c, k in groups)
    say("one binary", len(groups) == 1 and whole,
        "; ".join(f"commit {c} code {k}: {len(v)} species ({', '.join(v)})"
                  for (c, k), v in sorted(groups.items(), key=lambda x: str(x[0]))) + named
        + ("" if whole else "; a binary.txt names no commit or no code sha256"))

    same, rows_n = 0, 0
    for slug in slugs:
        got = read_tsv(runs / slug / "expect.tsv")
        rows_n += len(got)
        try:
            same += got == expected(runs / slug, slug, texts[slug][0])
        except (OSError, KeyError, ValueError):
            pass
    say("expect.py", same == len(slugs), f"{same} of {len(slugs)} expect.tsv files are expect.py's; {rows_n} rows")

    problems = []
    for slug in slugs:
        problems += [f"{slug}: {p}" for p in check_cohorts(runs / slug, slug, selection)]
    failed = {}
    for slug in slugs:
        p = runs / slug / "results" / "fetch.json"
        if p.is_file():
            n = int(json.loads(p.read_text(encoding="utf-8")).get("failed_added", 0) or 0)
            if n:
                failed[slug] = n
    say("cohorts", not problems, f"{len(slugs) - len({p.split(':')[0] for p in problems})} of {len(slugs)} cohorts.json "
        "equal the selection rule's counts" + "".join(f"; {p}" for p in problems)
        + "".join(f"; {s}: {n} added genomes failed the fetch" for s, n in sorted(failed.items())))

    ids = release_ids(release, slugs, lists)
    files = [p for slug in slugs for p in published(runs / slug)]
    tables = [data for slug in slugs for name, data in sorted(scores[slug].items())
              if not (runs / slug / name).is_file()]
    n_ids, n_bp = scan(files, ids, tables)
    say("genome ids and BioProjects", n_ids == 0 and n_bp == 0,
        f"{n_ids} genome ids and {n_bp} BioProject accessions in {len(files) + len(tables)} files, against the "
        f"{len(ids):,} genome ids of the release's {len(slugs)} species lists and {len(lists)} id lists; "
        "score/score.log not read")

    refused = []
    for slug in slugs:
        logs_dir = runs / slug / "logs"
        for p in sorted(logs_dir.glob("folds-*.log")) if logs_dir.is_dir() else []:
            for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
                m = REFUSED.match(line)
                if m:
                    refused.append((m.group(1), int(m.group(2))))
    say("refused lineage draws", True, f"{len(refused)} datasets with the lineage split unavailable, "
        f"{sum(n for _, n in refused)} refused default slots" + "".join(f"; {d}: {n}" for d, n in refused))
    return out, ok


def main(argv):
    if len(argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    lrows = read_tsv(HERE / "lists.tsv")
    lines, ok = run(argv[0], argv[1], lrows, lists=argv[3:])
    text = "\n".join(lines) + "\n"
    Path(argv[2]).write_text(text, encoding="utf-8")
    sys.stdout.write(text)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
