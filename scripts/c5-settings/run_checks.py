"""The checks of plan block C5's run copies, made before any score is read (docs/settings-arms.md).

  python3 -I run_checks.py check <copies dir> <lists.tsv> <release dir> [<checks.txt>]
  python3 -I run_checks.py gather <copies dir> <lists.tsv> <checks.txt> <reading dir>

<copies dir> holds one folder per species, named by its slug, as
run_species.sh left its OUT: results/, check/, score/, logs/, binary.txt,
lists.txt, expect.tsv, SHA256SUMS and DONE (or FAILED). Any other file
beside the folders, such as an archive of them, is not read. <lists.tsv>
is scripts/c5-settings/lists.tsv: the first run is its rows with first_run
1, none in M, and each has a study run. <release dir> is the rebuilt
release the run read; its species/<slug>.json lists the genome ids the
count looks for.

`check` prints one line per check, PASS or FAIL with its counts, a SET
ASIDE line for each FAILED file older than its species' DONE, and last
`CHECKS <passed> of <made>`; with <checks.txt> it writes the same lines
there. It exits 1 unless every check passes. It reads no score: of a
default result it reads the number of its outer folds alone, to know the
choices files the run owes, and of score/draws.tsv the names and the
constant column alone, to count the constant-model rows against the
refusal records. It never opens score/score.log, which stays out of the
record.

`gather` copies, for each species, score.py's checks.tsv, draws.tsv and
datasets.tsv, the run's binary.txt and every choices file into
<reading dir>/score/<slug>/, the layout page_tables.py reads, once
<checks.txt> says that every check passed. It copies nothing else: no
score.log, no other result, no log.
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


expect = _load("c5_settings_expect_for_checks", HERE / "expect.py")

CHOICE_ARMS = tuple(expect.ARMS)
SEEDS, TWIN_SEEDS = (42, 43, 44), (45, 46, 47)
LEARNERS = 5
TOKEN = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
NAMED = ("binary.txt", "expect.tsv", "lists.txt", "logs/timing.tsv", "logs/drive-timing.tsv")
SCORE_LOG = "score/score.log"
STAMP = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ")
REFUSED_WHY = "constant model: refused for no candidate"
DONE_KEYS = ["DONE", "species", "time", "binary", "grm", "code", "code_sha256", "results", "n", "choices", "nc",
             "SHA256SUMS", "sums"]
# What gather copies from a species' copy into the reading, and where.
GATHER = (("score/checks.tsv", "checks.tsv"), ("score/draws.tsv", "draws.tsv"), ("score/datasets.tsv", "datasets.tsv"),
          ("binary.txt", "binary.txt"))
CHOICES = "results/picks/choices"


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
        done = dict(zip(DONE_KEYS, words)) if len(words) == len(DONE_KEYS) else {"bad": " ".join(words)}
    if (out / "FAILED").is_file():
        words = (out / "FAILED").read_text(encoding="utf-8").split()
        failed = words[-1] if words else ""
    if done is None:
        return None, failed, None, "no DONE" + (f"; FAILED at {failed}" if failed is not None else "")
    if "bad" in done:
        return None, failed, None, f"DONE does not read as run_species.sh writes it: {done['bad']}"
    if failed is not None:
        # Both times as run_species.sh writes them (date -u +%FT%TZ); a FAILED with no such time is never set aside.
        if not (STAMP.fullmatch(failed) and STAMP.fullmatch(done["time"])):
            return None, failed, None, f"FAILED names no time to set beside DONE's: '{failed}'"
        if failed < done["time"]:
            return done, failed, f"{slug}/FAILED at {failed}, older than its DONE at {done['time']}", None
        return None, failed, None, f"FAILED at {failed} is not older than DONE at {done['time']}"
    return done, None, None, None


def expect_rows(lrows, slug):
    """expect.py's rows for one species: every first-run dataset is in S and has a study run."""
    run = [r for r in first_run(lrows) if r["species_slug"] == slug]
    datasets = [r["dataset"] for r in run]
    return [[str(x) for x in row] for row in expect.rows(slug, datasets, datasets, list(SEEDS), list(TWIN_SEEDS))]


def fold_count(path):
    """The number of outer folds of a default result: the one figure of it this script reads."""
    with open(path, encoding="utf-8") as handle:
        return len(json.load(handle).get("folds") or [])


def owed(picks, want):
    """The result files the run owes: each default and twin result, and for each choice arm's row one choices file
    per outer fold of the default's result at its dataset, split and seed (five where that result is missing)."""
    cols = list(expect.COLUMNS)
    rows = [dict(zip(cols, r)) for r in want]
    base = {(e["dataset"], e["split"], e["seed"]): e["file"] for e in rows if e["arm"] == "default"}
    needed = set()
    for e in rows:
        if e["arm"] not in CHOICE_ARMS:
            needed.add(e["file"])
            continue
        b = picks / base[(e["dataset"], e["split"], e["seed"])]
        n = fold_count(b) if b.is_file() else expect.FOLDS
        needed |= {e["file"].replace("<k>", str(k)) for k in range(n)}
    return needed


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
        # score/score.log is never opened, not even to hash it when a SHA256SUMS names it (the line above fails then).
        wrong = [n for n, sha in listed.items()
                 if n != SCORE_LOG and (not (out / n).is_file() or sha256(out / n) != sha)]
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
        picks = out / "results" / "picks"
        n_picks = len(list(picks.rglob("*.pick.json")))
        n_choices = len(list((out / CHOICES).rglob("*.choices.json")))
        if (done["species"], done["grm"], done["code_sha256"], done["n"], done["nc"]) \
                != (slug, f["grm"][0], f["code"][0], str(n_picks), str(n_choices)) \
                or not (out / "SHA256SUMS").is_file() or done["sums"] != sha256(out / "SHA256SUMS"):
            bad_done.append(f"{slug}: DONE does not agree with binary.txt, the results, the choices files or "
                            "SHA256SUMS")
    ok("each DONE names its species, its binary.txt's binary and code, its results, its choices files and its "
       "SHA256SUMS", not bad_done, "; ".join(bad_done))
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
        if (set(got) != {"species", "datasets", "seeds", "twin"}
                or (got["species"], sorted(x for x in got["datasets"].split(",") if x), got["seeds"], got["twin"])
                != (slug, want_d, ",".join(map(str, SEEDS)), ",".join(map(str, TWIN_SEEDS)))):
            bad_lists.append(slug)
    ok("each lists.txt names its species' first-run datasets, seeds 42, 43 and 44 and twin 45, 46 and 47",
       not bad_lists, ", ".join(bad_lists))

    bad_expect, bad_results, rows_n, files_n, refusals_n, constant_n, bad_constant = [], [], 0, 0, 0, 0, []
    for slug in present:
        out = copies / slug
        want = expect_rows(lrows, slug)
        rows_n += len(want)
        got = []
        if (out / "expect.tsv").is_file():
            with open(out / "expect.tsv", newline="", encoding="utf-8") as handle:
                got = [r for r in csv.reader(handle, delimiter="\t")]
        if not got or got[0] != list(expect.COLUMNS) or sorted(got[1:]) != sorted(want):
            bad_expect.append(slug)
        picks = out / "results" / "picks"
        needed = owed(picks, want)
        have = {p.relative_to(picks).as_posix() for p in picks.rglob("*") if p.is_file()} if picks.is_dir() else set()
        refused = {n for n in have if n.endswith(".refused.json")}
        files_n += len(needed)
        if have - refused != needed:
            bad_results.append(f"{slug}: {len(needed - have)} missing, {len(have - refused - needed)} more")
        refusals_n += len(refused)
        if (out / "score" / "draws.tsv").is_file():
            n_constant = sum(1 for r in load(out / "score" / "draws.tsv")
                             if r["constant"] == "1" or r["why"] == REFUSED_WHY)
            constant_n += n_constant
            if n_constant != LEARNERS * len(refused):
                bad_constant.append(slug)
        else:
            bad_constant.append(slug)
    ok("each expect.tsv is expect.py's list from lists.tsv", not bad_expect,
       f"{rows_n} results" + (f"; differs in {', '.join(bad_expect)}" if bad_expect else ""))
    ok("each species holds every default and twin result and one choices file per outer fold of each draw and "
       "study run of the default, and no other", not bad_results,
       f"{files_n} files" + ("; " + "; ".join(bad_results) if bad_results else ""))
    ok("the refusal records and the constant-model rows, counted: five constant-model rows per refusal record",
       not bad_constant, f"{refusals_n} refusal records, {constant_n} constant-model rows"
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
        (target / slug / "choices").mkdir(parents=True, exist_ok=True)
        for src, dst in GATHER:
            (target / slug / dst).write_bytes((copies / slug / src).read_bytes())
            n += 1
        for p in sorted((copies / slug / CHOICES).glob("*.choices.json")):
            (target / slug / "choices" / p.name).write_bytes(p.read_bytes())
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
