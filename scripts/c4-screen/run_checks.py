"""The checks of plan block C4's run copies, made before any score is read (docs/screen-arms.md).

  python3 -I run_checks.py check <copies dir> <lists.tsv> <release dir> [<checks.txt>]
  python3 -I run_checks.py gather <copies dir> <lists.tsv> <checks.txt> <reading dir>

<copies dir> holds one folder per species, named by its slug, as
run_species.sh left its OUT: results/, check/, score/, logs/, binary.txt,
lists.txt, expect.tsv, SHA256SUMS and DONE (or FAILED). <lists.tsv> is
scripts/c4-screen/lists.tsv: the first run is its rows with first_run 1.
<release dir> is the rebuilt release the run read; its
species/<slug>.json lists the genome ids the count looks for.

`check` prints one line per check, PASS or FAIL with its counts, a SET
ASIDE line for each FAILED file older than its species' DONE, and last
`CHECKS <passed> of <made>`; with <checks.txt> it writes the same lines
there. It exits 1 unless every check passes. It reads no score: of
score/draws.tsv it reads the names and the constant column alone, to count
the constant-model rows against the refusal records, and it never opens
score/score.log, which stays out of the record.

`gather` copies, for each species, score.py's checks.tsv, draws.tsv and
datasets.tsv, the run's logs/fits.tsv (as fits.tsv) and binary.txt into
<reading dir>/score/<slug>/, the layout page_tables.py reads, once
<checks.txt> says that every check passed. It copies nothing else: no
score.log, no result, no log.
"""

import csv
import hashlib
import importlib.util
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


expect = _load("c4_screen_expect_for_checks", HERE / "expect.py")

ARMS = tuple(expect.ARMS)
SEEDS, TWIN_SEEDS = (42, 43, 44), (45, 46, 47)
LEARNERS = 5
TOKEN = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
NAMED = ("binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv", "logs/fits.tsv")
SCORE_LOG = "score/score.log"
REFUSED_WHY = "constant model: refused for no candidate"
# What gather copies from a species' copy into the reading, and where.
GATHER = (("score/checks.tsv", "checks.tsv"), ("score/draws.tsv", "draws.tsv"), ("score/datasets.tsv", "datasets.tsv"),
          ("logs/fits.tsv", "fits.tsv"), ("binary.txt", "binary.txt"))


def load(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def first_run(lrows):
    return [r for r in lrows if r["first_run"] == "1" and r["in_M"] == "0"]


def covered(out):
    """The files run_species.sh's SHA256SUMS covers: every .json, .tsv and .txt under results/, check/ and score/,
    and the five it names."""
    names = [p.relative_to(out).as_posix() for top in ("results", "check", "score") if (out / top).is_dir()
             for p in (out / top).rglob("*") if p.is_file() and p.suffix in (".json", ".tsv", ".txt")]
    return sorted(set(names) | set(NAMED))


def fields(path):
    """binary.txt's commit, grm and code lines."""
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition(" ")
        if key in ("commit", "grm", "code"):
            out.setdefault(key, []).append(value.strip())
    return out


def ended(slug, out):
    """(DONE's fields or None, FAILED's time or None, a FAILED set aside or None, why not ended or None)."""
    done = failed = None
    if (out / "DONE").is_file():
        words = (out / "DONE").read_text(encoding="utf-8").split()
        keys = ["DONE", "species", "time", "binary", "grm", "code", "code_sha256", "results", "n", "SHA256SUMS", "sums"]
        done = dict(zip(keys, words)) if len(words) == len(keys) else {"bad": " ".join(words)}
    if (out / "FAILED").is_file():
        words = (out / "FAILED").read_text(encoding="utf-8").split()
        failed = words[-1] if words else ""
    if done is None:
        return None, failed, None, "no DONE" + (f"; FAILED at {failed}" if failed is not None else "")
    if "bad" in done:
        return None, failed, None, f"DONE does not read as run_species.sh writes it: {done['bad']}"
    if failed is not None:
        if failed < done["time"]:
            return done, failed, f"{slug}/FAILED at {failed}, older than its DONE at {done['time']}", None
        return None, failed, None, f"FAILED at {failed} is not older than DONE at {done['time']}"
    return done, None, None, None


def expect_rows(lrows, slug):
    run = [r for r in first_run(lrows) if r["species_slug"] == slug]
    datasets = [r["dataset"] for r in run]
    study = [r["dataset"] for r in run if r["in_S"] == "1"]
    return [[str(x) for x in row] for row in expect.rows(slug, datasets, study, list(SEEDS), list(TWIN_SEEDS),
                                                         list(ARMS))]


def refused_name(file):
    return file[: -len(".pick.json")] + ".refused.json"


def check(copies, lists_path, release):
    """[(passed, name, detail)] and [set-aside lines] over every species' copy."""
    copies, release = Path(copies), Path(release)
    lrows = load(lists_path)
    run = first_run(lrows)
    species = sorted({r["species_slug"] for r in run})
    lines, aside = [], []

    def ok(name, passed, detail):
        lines.append((bool(passed), name, detail))

    found = sorted(p.name for p in copies.iterdir() if p.is_dir())
    ok("the first run's species, each once", found == species,
       f"{len(found)} folders for {len(species)} species; missing {sorted(set(species) - set(found))}, more "
       f"{sorted(set(found) - set(species))}")
    present = [s for s in species if s in found]

    ends, bad_end = {}, []
    for slug in present:
        done, failed, set_aside, why = ended(slug, copies / slug)
        if why:
            bad_end.append(f"{slug}: {why}")
        else:
            ends[slug] = done
        if set_aside:
            aside.append(set_aside)
    ok("each species ended DONE, any FAILED older than its DONE", not bad_end and len(ends) == len(species),
       f"{len(ends)} of {len(species)}" + ("; " + "; ".join(bad_end) if bad_end else ""))

    bad_sums, listed_n = [], 0
    for slug in present:
        out = copies / slug
        if not (out / "SHA256SUMS").is_file():
            bad_sums.append(f"{slug}: no SHA256SUMS")
            continue
        listed = {}
        for line in (out / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            sha, _, name = line.partition("  ")
            listed[name] = sha
        listed_n += len(listed)
        if sorted(listed) != covered(out):
            bad_sums.append(f"{slug}: SHA256SUMS lists {len(listed)} files, the run's rule {len(covered(out))}")
        wrong = [n for n, sha in listed.items() if not (out / n).is_file() or sha256(out / n) != sha]
        if wrong:
            bad_sums.append(f"{slug}: {len(wrong)} files differ from SHA256SUMS")
    ok("each species' SHA256SUMS holds every file of the run's rule, each at its sha256", not bad_sums,
       f"{listed_n} files" + ("; " + "; ".join(bad_sums) if bad_sums else ""))
    ok("score/score.log stays out of SHA256SUMS", all(
        SCORE_LOG not in (copies / s / "SHA256SUMS").read_text(encoding="utf-8") for s in present
        if (copies / s / "SHA256SUMS").is_file()), "")

    commits, codes, grms, bad_done = set(), set(), set(), []
    for slug in present:
        out = copies / slug
        if not (out / "binary.txt").is_file():
            bad_done.append(f"{slug}: no binary.txt")
            continue
        f = fields(out / "binary.txt")
        if any(len(f.get(k, [])) != 1 for k in ("commit", "grm", "code")):
            bad_done.append(f"{slug}: binary.txt does not name one commit, binary and code")
            continue
        commits.add(f["commit"][0])
        grms.add(f["grm"][0])
        codes.add(f["code"][0])
        done = ends.get(slug)
        if done is None:
            continue
        picks = len(list((out / "results" / "picks").rglob("*.pick.json")))
        if (done["species"], done["grm"], done["code_sha256"], done["n"]) != (slug, f["grm"][0], f["code"][0],
                                                                            str(picks)) \
                or not (out / "SHA256SUMS").is_file() or done["sums"] != sha256(out / "SHA256SUMS"):
            bad_done.append(f"{slug}: DONE does not agree with binary.txt, the results or SHA256SUMS")
    ok("each DONE names its species, its binary.txt's binary and code, its results and its SHA256SUMS", not bad_done,
       "; ".join(bad_done))
    ok("every species ran one commit and one code sha256", len(commits) == 1 and len(codes) == 1,
       f"{len(commits)} commits, {len(codes)} code sha256, {len(grms)} binaries")

    bad_lists = []
    for slug in present:
        out = copies / slug
        want_d = sorted(r["dataset"] for r in run if r["species_slug"] == slug)
        got = {}
        if (out / "lists.txt").is_file():
            for line in (out / "lists.txt").read_text(encoding="utf-8").splitlines():
                key, _, value = line.partition("\t")
                got[key] = value
        if (got.get("species"), sorted(x for x in got.get("datasets", "").split(",") if x), got.get("seeds"),
                got.get("twin"), got.get("arms")) != (slug, want_d, ",".join(map(str, SEEDS)),
                                                      ",".join(map(str, TWIN_SEEDS)), ",".join(ARMS)):
            bad_lists.append(slug)
    ok("each lists.txt names its species' first-run datasets, seeds 42, 43 and 44, twin 45, 46 and 47 and the "
       "three arms", not bad_lists, ", ".join(bad_lists))

    bad_expect, bad_results, rows_n, refusals_n, constant_n, bad_constant = [], [], 0, 0, 0, []
    for slug in present:
        out = copies / slug
        want = expect_rows(lrows, slug)
        rows_n += len(want)
        got = []
        if (out / "expect.tsv").is_file():
            with open(out / "expect.tsv", newline="", encoding="utf-8") as handle:
                got = [r for r in csv.reader(handle, delimiter="\t")]
        if not got or got[0] != list(expect.COLUMNS) or got[1:] != want:
            bad_expect.append(slug)
        picks = out / "results" / "picks"
        cols = list(expect.COLUMNS)
        needed, refused = set(), set()
        for row in want:
            e = dict(zip(cols, row))
            if (picks / e["file"]).is_file():
                needed.add(e["file"])
                if e["inputs"] != "-":
                    needed.add(e["inputs"])
            elif e["arm"] in ARMS and e["split"] == "study" and (picks / refused_name(e["file"])).is_file():
                needed.add(refused_name(e["file"]))
                refused.add((e["dataset"], e["arm"], e["split"], e["seed"]))
            else:
                needed.add(e["file"])
        have = {p.relative_to(picks).as_posix() for p in picks.rglob("*") if p.is_file()} if picks.is_dir() else set()
        if have != needed:
            bad_results.append(f"{slug}: {len(needed - have)} missing, {len(have - needed)} more")
        refusals_n += len(refused)
        marked = {}
        if (out / "score" / "draws.tsv").is_file():
            for r in load(out / "score" / "draws.tsv"):
                if r["constant"] == "1" or r["why"] == REFUSED_WHY:
                    key = (r["dataset"], r["arm"], r["split"], r["seed"])
                    marked[key] = marked.get(key, 0) + (r["constant"] == "1" and r["why"] == REFUSED_WHY)
                    constant_n += 1
        if set(marked) != refused or any(n != LEARNERS for n in marked.values()):
            bad_constant.append(slug)
    ok("each expect.tsv is expect.py's list from lists.tsv", not bad_expect,
       f"{rows_n} results" + (f"; differs in {', '.join(bad_expect)}" if bad_expect else ""))
    ok("each species holds every listed result or refusal record and inputs file, and no other", not bad_results,
       "; ".join(bad_results))
    ok("each refusal record is a study run scored as the constant model on five learners, and each constant-model "
       "row a refusal record's", not bad_constant,
       f"{refusals_n} refusal records, {constant_n} constant-model rows"
       + (f"; differs in {', '.join(bad_constant)}" if bad_constant else ""))

    ids, hits, scanned = set(), 0, 0
    for slug in present:
        path = release / "species" / f"{slug}.json"
        if path.is_file():
            with open(path, encoding="utf-8") as handle:
                ids |= {g["id"] for g in json.load(handle)["genomes"]}
    missing_ids = [s for s in present if not (release / "species" / f"{s}.json").is_file()]
    for slug in present:
        out = copies / slug
        names = set(covered(out)) | {"DONE", "FAILED", "SHA256SUMS"}
        for name in sorted(names):
            if name == SCORE_LOG or not (out / name).is_file():
                continue
            scanned += 1
            text = (out / name).read_text(encoding="utf-8", errors="replace")
            hits += sum(1 for m in TOKEN.finditer(text) if m.group(0) in ids)
    ok("no genome id in any file that will be read or published", hits == 0 and not missing_ids and ids,
       f"{hits} genome ids in {scanned} files, against {len(ids)} ids"
       + (f"; no genome list for {', '.join(missing_ids)}" if missing_ids else ""))
    return lines, aside


def report(lines, aside):
    out = [f"{'PASS' if p else 'FAIL'} {name}" + (f": {detail}" if detail else "") for p, name, detail in lines]
    out += [f"SET ASIDE {a}" for a in aside]
    out.append(f"CHECKS {sum(p for p, _, _ in lines)} of {len(lines)}")
    return out


def all_passed(checks_path):
    lines = Path(checks_path).read_text(encoding="utf-8").splitlines()
    tail = re.fullmatch(r"CHECKS (\d+) of (\d+)", lines[-1]) if lines else None
    return bool(tail) and tail.group(1) == tail.group(2) and not any(x.startswith("FAIL") for x in lines)


def gather(copies, lists_path, checks_path, reading):
    copies, reading = Path(copies), Path(reading)
    if not all_passed(checks_path):
        raise SystemExit(f"{checks_path} does not say that every check passed")
    target = reading / "score"
    if target.exists() and any(target.iterdir()):
        raise SystemExit(f"{target} is not empty")
    species = sorted({r["species_slug"] for r in first_run(load(lists_path))})
    n = 0
    for slug in species:
        (target / slug).mkdir(parents=True, exist_ok=True)
        for src, dst in GATHER:
            (target / slug / dst).write_bytes((copies / slug / src).read_bytes())
            n += 1
    return n


def main(argv):
    if argv and argv[0] == "check" and len(argv) in (4, 5):
        lines, aside = check(argv[1], argv[2], argv[3])
        text = report(lines, aside)
        if len(argv) == 5:
            Path(argv[4]).write_text("\n".join(text) + "\n", encoding="utf-8")
        print("\n".join(text))
        return 0 if all(p for p, _, _ in lines) else 1
    if argv and argv[0] == "gather" and len(argv) == 5:
        n = gather(argv[1], argv[2], argv[3], argv[4])
        print(f"gathered {n} files")
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
