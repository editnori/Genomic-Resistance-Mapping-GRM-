"""The checks of plan block C7's run, before its count is taken (docs/outside-learners-result.md, Method).

  python3 -I run_checks.py check <back dir> <lists.tsv> <release dir> [--withdrawn ARM REASON TIME] [--ids-also PATH ...]
  python3 -I run_checks.py gather <back dir> <lists.tsv> <release dir> <reading dir> [--withdrawn ARM REASON TIME]
                                  [--ids-also PATH ...]

<back dir> holds one copy per species, as the run's copies came back: each
<species>/ with results/outside/, logs/, lists.txt, libraries.txt,
expect.tsv, export-binary.txt, SHA256SUMS and DONE, and score.tar, which
holds score/ and is never unpacked. A folder <species>.<UTC time> beside
them is a species' first copy, moved aside when the species restarted for
its second pass. When the elastic net is withdrawn after some second passes
ran, a species' copy is the one whose DONE says the elastic net was left
out: <species>/ where no second pass ran, and where one did, its first copy
<species>.<UTC time>/; its second-pass copy, <species>/ and any copy moved
aside after the first, is named, counted and set aside, and never opened.
<lists.tsv> is scripts/c7-outside/lists.tsv. <release
dir> is the public release unpacked, whose species/<species>.json lists
each species' genome ids. --withdrawn names a library arm the lead
withdrew before it ran, with the reason and the UTC time the lead gave.

`check` prints one line per check, PASS or FAIL, with counts only, and
exits 1 when one fails:

  species      the copies are exactly the species of lists.tsv's first run;
               each folder <species>.<UTC time> is a first copy of a species
               of the run, older than its species' DONE: named and set aside;
               with the elastic net withdrawn, a species' first copy whose
               DONE says it was left out is its copy, and each second-pass
               copy beside it is named with its species and counted, never
               opened
  SHA256SUMS   every line of each species' SHA256SUMS holds, score/'s lines
               read from score.tar in place, and DONE names its sha256; the
               files under results/ and score/ are its lines
  DONE         each DONE names its species, the export binary
               export-binary.txt names, its result count and SHA256SUMS;
               one export commit, one grm binary and one code sha256 (the C3
               run's binary.txt, as DONE gives them) in every species; the
               elastic net run in every species, or left out in every
               species when it is withdrawn; a FAILED beside a newer DONE is
               named and set aside, and any other FAILED fails
  expect.py    lists.txt names the species' first-run datasets, the seeds,
               twin seeds, arm seeds and the arms DONE's elasticnet field
               gives, expect.tsv is expect.py's for them (the elastic net's
               rows owed only when it ran), and results/outside/ holds
               exactly the library arms' result files of its rows
  genome ids   no genome id of the species' release lists in any file that
               is read or published: results/, logs/, the copy's top files,
               score.tar's tables, and --ids-also
  features     no unitig or k-mer (31 or more nucleotides) and no column name
               of a booster's matrix (LightGBM's Column_<n>, XGBoost's
               f<n>) in the same files

The last line is `CHECKS <passed> of <made>`. score.tar's score/score.log,
and on a new box its c3/score/, are neither read nor counted nor copied: the
record keeps score.py's three tables alone. back.tar, the archive a copy came
in, is not read. No score is printed.

`gather` refuses a withdrawal whose time is later than its own clock read. It
runs the same checks and, only when every one passes, writes the
reading directory page_tables.py reads: checks.txt (the printout), withdrawn.tsv,
reading.txt (this command's outputs with their sha256), score/<species>/ with
checks.tsv, draws.tsv and datasets.tsv read from score.tar, and run/<species>/
with DONE, libraries.txt, timing.tsv, plans.log and outside/. It prints the
count's command, which writes count/; that command and its output's sha256
go into reading.txt after it.
"""

import csv
import hashlib
import importlib.util
import io
import json
import re
import shutil
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEEDS, TWIN, ARM_SEEDS = "42,43,44", "45,46,47", "42,43,44"
ARMS = ("xgboost", "lightgbm", "elasticnet")
ARMS_BY_PASS = {"run": "xgboost,lightgbm,elasticnet", "leave-out": "xgboost,lightgbm"}
# DONE and FAILED give their time as date -u +%FT%TZ writes it; any other form is no DONE or FAILED line.
DONE_LINE = re.compile(r"DONE (?P<species>\S+) (?P<time>\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ) grm (?P<grm>[0-9a-f]{64}) code (?P<code>[0-9a-f]{64}) "
                       r"export-grm (?P<export>[0-9a-f]{64}) elasticnet (?P<elasticnet>run|leave-out) "
                       r"results (?P<results>\d+) SHA256SUMS (?P<sums>[0-9a-f]{64})")
FAILED_LINE = re.compile(r"FAILED (?P<species>\S+) (?P<how>at|stopped at) (?P<step>\S+) "
                         r"(?P<time>\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ)")
ASIDE = re.compile(r"(?P<slug>[a-z0-9_]+)\.(?P<time>(?P<d>\d{8})T(?P<t>\d{6})Z)")
TIME = re.compile(r"\d{4}-\d\d-\d\dT\d\d:\d\d(:\d\d)?Z")
# A genome id is <taxon>.<n>; a letter or _ before it, as in genome_<id> or id<id>, still leaves it an id.
ID_TOKEN = re.compile(r"(?<![\d.])\d+\.\d+(?!\d)")
# A unitig or k-mer of 31 or more nucleotides, and the column names LightGBM and XGBoost give a matrix's features.
FEATURE_TOKENS = (re.compile(r"[ACGTacgt]{31,}"), re.compile(r"Column_\d+"), re.compile(r"(?<![\w.])f\d+(?![\w.])"))
READ = ("results", "logs")
TOP_FILES = ("lists.txt", "libraries.txt", "expect.tsv", "export-binary.txt", "SHA256SUMS", "DONE")
SCORE_TABLES = ("score/checks.tsv", "score/draws.tsv", "score/datasets.tsv")
# What score.tar may hold besides score.py's tables: its output, and on a new box stage 1's score/. Never read.
NOT_READ = ("score/score.log",)
NOT_READ_PREFIX = ("c3/score/",)


def _expect():
    spec = importlib.util.spec_from_file_location("c7_outside_expect_for_checks", HERE / "expect.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


expect = _expect()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def iso(aside_time):
    """%Y%m%dT%H%M%SZ as %FT%TZ, the time DONE and FAILED give."""
    d, t = aside_time[:8], aside_time[9:15]
    return f"{d[:4]}-{d[4:6]}-{d[6:]}T{t[:2]}:{t[2:4]}:{t[4:]}Z"


def score_files(copy, skipped=None):
    """{name: bytes} of score.py's tables in one copy's score.tar, read in place. score.log and c3/score/ are skipped
    unread, each member's name added to skipped when it is given; any other member stops the check."""
    out = {}
    tar = copy / "score.tar"
    if not tar.is_file():
        return out
    with tarfile.open(tar) as archive:
        for member in archive:
            if not member.isfile():
                continue
            name = member.name[2:] if member.name.startswith("./") else member.name
            if name in NOT_READ or name.startswith(NOT_READ_PREFIX):
                if skipped is not None:
                    skipped.append(name)
                continue
            if not name.startswith("score/") or not name.endswith((".tsv", ".json")):
                raise SystemExit(f"{copy.name}: score.tar holds a file that is no table of score/")
            with archive.extractfile(member) as handle:
                out[name] = handle.read()
    return out


def lists_txt(text):
    return dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)


def expect_rows(rows):
    return [dict(zip(expect.COLUMNS, (str(x) for x in r))) for r in rows]


def load_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def done_of(copy):
    path = copy / "DONE"
    return DONE_LINE.fullmatch(path.read_text(encoding="utf-8").strip()) if path.is_file() else None


class Report:
    def __init__(self):
        self.lines, self.failed, self.checks = [], 0, 0
        self.copy_of = {}

    def check(self, name, good, detail, more=()):
        self.checks += 1
        self.failed += not good
        self.lines.append(f"{'PASS' if good else 'FAIL'} {name}: {detail}")
        self.lines += [f"  {m}" for m in more]


def run(back, l_path, release, withdrawn=None, ids_also=()):
    back, release = Path(back), Path(release)
    withdrawn = withdrawn or {}
    with open(l_path, encoding="utf-8", newline="") as handle:
        lrows = list(csv.DictReader(handle, delimiter="\t"))
    first = [r for r in lrows if r["first_run"] == "1"]
    species = sorted({r["species_slug"] for r in first})
    rep = Report()

    # species
    folders = sorted(p.name for p in back.iterdir() if p.is_dir())
    got = [f for f in folders if ASIDE.fullmatch(f) is None]
    asides = sorted((ASIDE.fullmatch(f) for f in folders if ASIDE.fullmatch(f) is not None), key=lambda m: m["time"])
    copy_of, second = {s: back / s for s in species if (back / s).is_dir()}, {}
    late = "elasticnet" in withdrawn
    if late:
        # Withdrawn after some second passes ran: a species' copy is its first copy whose DONE says the elastic net
        # was left out; <species>/ and every copy moved aside after it are its second pass, set aside unopened.
        for s in species:
            mine = [m for m in asides if m["slug"] == s]
            for i, m in enumerate(mine):
                d = done_of(back / m.group(0))
                if d is not None and d["species"] == s and d["elasticnet"] == "leave-out":
                    copy_of[s] = back / m.group(0)
                    second[s] = [x.group(0) for x in mine[i + 1:]] + ([s] if s in got else [])
                    break
    have = sorted(copy_of)
    dones = {s: done_of(copy_of[s]) for s in have}
    named, bad_aside, set_aside = [], [], {x for v in second.values() for x in v}
    for m in asides:
        done = dones.get(m["slug"])
        if copy_of.get(m["slug"]) == back / m.group(0) or m.group(0) in set_aside:
            continue
        if m["slug"] in species and done is not None and iso(m["time"]) < done["time"]:
            named.append(f"{m.group(0)}: a first copy of {m['slug']}, older than its DONE at {done['time']}: named and "
                         "set aside, never read")
        else:
            bad_aside.append(m.group(0))
    by_time = len(named)
    named += [f"{copy_of[s].name}: the copy of {s}, its DONE saying the elastic net was left out"
              for s in sorted(second)]
    named += [f"{x}: the second pass of {s}, after the elastic net was withdrawn: named, counted and set aside, never "
              "opened" for s in sorted(second) for x in second[s]]
    outside = sorted(set(got) - set(species))
    rep.check("species", have == species and not outside and not bad_aside,
              f"{len(got)} copies, {len(species)} species in lists.tsv's first run; {by_time} first copies set aside "
              f"by time, {len(bad_aside)} beside no newer DONE"
              + (f"; the elastic net withdrawn: {len(second)} species ran a second pass, {len(set_aside)} second-pass "
                 "copies set aside" if late else ""),
              ([f"missing {len(set(species) - set(have))}, more {len(outside)}"] if have != species or outside
               else []) + named + [f"{x}: no newer DONE of a species of the run" for x in bad_aside])
    rep.copy_of = copy_of
    skipped = {s: [] for s in have}
    scores = {s: score_files(copy_of[s], skipped[s]) for s in have}
    stage1 = sum(1 for s in have if any(n.startswith(NOT_READ_PREFIX) for n in skipped[s]))

    # SHA256SUMS
    good, lines, bad_lines, listed_ok = 0, 0, 0, 0
    for s in have:
        copy = copy_of[s]
        sums = copy / "SHA256SUMS"
        if not sums.is_file():
            continue
        ok, listed = True, set()
        for line in sums.read_text(encoding="utf-8").splitlines():
            want, name = line[:64], line[66:]
            listed.add(name)
            lines += 1
            data = scores[s].get(name) if name.startswith("score/") else (
                (copy / name).read_bytes() if (copy / name).is_file() else None)
            if data is None or sha(data) != want:
                bad_lines += 1
                ok = False
        m = dones[s]
        ok &= m is not None and m["sums"] == sha(sums.read_bytes())
        on_disk = {p.relative_to(copy).as_posix() for p in (copy / "results").rglob("*")
                   if p.is_file() and p.suffix in (".json", ".tsv")} if (copy / "results").is_dir() else set()
        on_disk |= set(scores[s])
        in_sums = {n for n in listed if n.split("/")[0] in ("results", "score")}
        same = on_disk == in_sums and set(SCORE_TABLES) <= set(scores[s])
        listed_ok += same
        good += ok and same
    rep.check("SHA256SUMS", good == len(species) and bad_lines == 0,
              f"{good} of {len(species)} species; {lines - bad_lines} of {lines} lines hold; results/ and score.tar's "
              f"tables are the listed files in {listed_ok} of {len(species)}; stage 1's c3/score/ in the score.tar of "
              f"{stage1} species, skipped unread")

    # DONE
    named_ok, by = 0, {"commit": {}, "grm": {}, "code": {}, "export": {}, "elasticnet": {}}
    aside, failed = [], []
    for s in have:
        copy = copy_of[s]
        m = dones[s]
        export = {}
        if (copy / "export-binary.txt").is_file():
            for line in (copy / "export-binary.txt").read_text(encoding="utf-8").splitlines():
                parts = line.split()
                if len(parts) == 2:
                    export[parts[0]] = parts[1]
        results = sum(1 for _ in (copy / "results" / "outside").rglob("*.outside.json")) \
            if (copy / "results" / "outside").is_dir() else -1
        if m and m["species"] == s and m["export"] == export.get("grm") and int(m["results"]) == results:
            named_ok += 1
        if m:
            for k in ("grm", "code", "export", "elasticnet"):
                by[k].setdefault(m[k], []).append(s)
        by["commit"].setdefault(export.get("commit", "none"), []).append(s)
        if (copy / "FAILED").exists():
            f = FAILED_LINE.fullmatch((copy / "FAILED").read_text(encoding="utf-8").strip())
            if f is not None and f["species"] == s and m is not None and f["time"] < m["time"]:
                aside.append(f"{s}: FAILED {f['how']} {f['step']} at {f['time']}, before DONE at {m['time']}: named "
                             "and set aside")
            else:
                failed.append(s)
    want_pass = "leave-out" if "elasticnet" in withdrawn else "run"
    # One export commit, grm binary, code and export binary; the elastic net's pass is held by passes_ok alone.
    one = {k: len(v) == 1 for k, v in by.items() if k != "elasticnet"}
    passes_ok = set(by["elasticnet"]) == {want_pass}
    rep.check("DONE", named_ok == len(species) and not failed and all(one.values()) and passes_ok,
              f"{named_ok} of {len(species)} DONE lines name their species, export binary, result count and "
              f"SHA256SUMS; FAILED not set aside in {len(failed)} species; the elastic net "
              f"{'left out, withdrawn' if want_pass == 'leave-out' else 'run'} in "
              f"{len(by['elasticnet'].get(want_pass, []))} of {len(species)}",
              [f"species per export commit: {'; '.join(f'{c[:7]} {len(v)}' for c, v in sorted(by['commit'].items()))}",
               f"species per grm sha256: {'; '.join(f'{c[:12]} {len(v)}' for c, v in sorted(by['grm'].items()))}",
               f"species per code sha256: {'; '.join(f'{c[:12]} {len(v)}' for c, v in sorted(by['code'].items()))}",
               f"species per export grm sha256: "
               f"{'; '.join(f'{c[:12]} {len(v)}' for c, v in sorted(by['export'].items()))}",
               f"species per elastic net pass: "
               f"{'; '.join(f'{c} {len(v)}' for c, v in sorted(by['elasticnet'].items()))}"]
              + [f"elastic net {p}: {', '.join(v)}" for p, v in sorted(by["elasticnet"].items()) if p != want_pass]
              + aside + [f"FAILED: {s}" for s in failed])

    # expect.py
    ok_species, found, expected, more, datasets, study = 0, 0, 0, 0, 0, 0
    for s in have:
        copy = copy_of[s]
        mine = [r for r in first if r["species_slug"] == s]
        ds = [r["dataset"] for r in mine]
        st = [r["dataset"] for r in mine if r["in_S"] == "1"]
        datasets += len(ds)
        study += len(st)
        m = dones[s]
        try:
            lt = lists_txt((copy / "lists.txt").read_text(encoding="utf-8"))
            arms = tuple(lt.get("arms", "").split(","))
            rows = expect_rows(expect.rows(s, ds, st, [int(x) for x in SEEDS.split(",")],
                                           [int(x) for x in TWIN.split(",")], [int(x) for x in ARM_SEEDS.split(",")],
                                           arms))
            got_rows = load_tsv((copy / "expect.tsv").read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        same_lists = (lt.get("species") == s and sorted(lt.get("datasets", "").split(",")) == sorted(ds)
                      and lt.get("seeds") == SEEDS and lt.get("twin") == TWIN and lt.get("arm_seeds") == ARM_SEEDS
                      and m is not None and lt.get("arms") == ARMS_BY_PASS[m["elasticnet"]])
        same_expect = sorted(map(lambda r: tuple(sorted(r.items())), got_rows)) == \
            sorted(map(lambda r: tuple(sorted(r.items())), rows))
        outside = copy / "results" / "outside"
        owed = {r["file"][len("results/outside/"):] for r in rows if r["file"].startswith("results/outside/")}
        on_disk = {p.relative_to(outside).as_posix() for p in outside.rglob("*") if p.is_file()} \
            if outside.is_dir() else set()
        expected += len(owed)
        found += len(owed & on_disk)
        more += len(on_disk - owed)
        ok_species += same_lists and same_expect and on_disk == owed
    rep.check("expect.py", ok_species == len(species) and found == expected and more == 0,
              f"{ok_species} of {len(species)} species: lists.txt names the first run, the seeds, twin seeds, arm "
              f"seeds and the arms of its pass, expect.tsv is expect.py's, and results/outside/ is exactly its library "
              f"arm rows; {found} of {expected} results found, {more} more",
              [f"datasets {datasets}, study datasets {study}, species {len(have)}"])

    # genome ids and features
    ids, lists_found = set(), 0
    for s in species:
        path = release / "species" / f"{s}.json"
        if path.is_file():
            lists_found += 1
            with open(path, encoding="utf-8") as handle:
                ids |= {g["id"] for g in json.load(handle)["genomes"]}
    texts = []
    for s in have:
        copy = copy_of[s]
        paths = [p for top in READ if (copy / top).is_dir() for p in sorted((copy / top).rglob("*")) if p.is_file()]
        paths += [copy / n for n in TOP_FILES if (copy / n).is_file()]
        texts += [p.read_bytes().decode("utf-8", errors="replace") for p in paths]
        texts += [data.decode("utf-8", errors="replace") for data in scores[s].values()]
    for extra in ids_also:
        extra = Path(extra)
        for p in ([extra] if extra.is_file() else sorted(q for q in extra.rglob("*") if q.is_file())):
            texts.append(p.read_bytes().decode("utf-8", errors="replace"))
    hits = sum(1 for t in texts for m in ID_TOKEN.finditer(t) if m.group(0) in ids)
    rep.check("genome ids", hits == 0 and lists_found == len(species) and len(ids) > 0,
              f"{hits} genome ids in {len(texts)} files, against {len(ids)} release ids of {lists_found} species lists")
    features = [sum(1 for t in texts for _ in pattern.finditer(t)) for pattern in FEATURE_TOKENS]
    rep.check("features", sum(features) == 0,
              f"{features[0]} runs of 31 or more nucleotides, {features[1]} LightGBM and {features[2]} XGBoost column "
              f"names in {len(texts)} files")
    rep.lines.append(f"CHECKS {rep.checks - rep.failed} of {rep.checks}")
    return rep


def later_than_now(time, now=None):
    """True when a withdrawal's UTC time (%Y-%m-%dT%H:%M or with seconds, then Z) is later than the clock."""
    now = now or datetime.now(timezone.utc)
    form = "%Y-%m-%dT%H:%M:%SZ" if len(time) == 20 else "%Y-%m-%dT%H:%MZ"
    return datetime.strptime(time, form).replace(tzinfo=timezone.utc) > now


def gather(back, l_path, release, reading, withdrawn=None, ids_also=()):
    """Write the reading directory from the copies, only when every check passes. Returns the report."""
    back, reading = Path(back), Path(reading)
    if reading.exists() and any(reading.iterdir()):
        raise SystemExit(f"{reading} is not empty")
    now = datetime.now(timezone.utc)
    for arm, (_, time) in sorted((withdrawn or {}).items()):
        if later_than_now(time, now):
            raise SystemExit(f"--withdrawn {arm}: its time {time} is later than this clock's "
                             f"{now.strftime('%Y-%m-%dT%H:%M:%SZ')}; nothing written")
    rep = run(back, l_path, release, withdrawn, ids_also)
    if rep.failed:
        return rep
    with open(l_path, encoding="utf-8", newline="") as handle:
        species = sorted({r["species_slug"] for r in csv.DictReader(handle, delimiter="\t") if r["first_run"] == "1"})
    reading.mkdir(parents=True, exist_ok=True)
    (reading / "checks.txt").write_text("\n".join(rep.lines) + "\n", encoding="utf-8")
    with open(reading / "withdrawn.tsv", "w", encoding="utf-8") as handle:
        handle.write("arm\treason\ttime\n")
        for arm, (reason, time) in sorted((withdrawn or {}).items()):
            handle.write(f"{arm}\t{reason}\t{time}\n")
    for s in species:
        copy = rep.copy_of[s]
        score = reading / "score" / s
        score.mkdir(parents=True)
        tables = score_files(copy)
        for name in SCORE_TABLES:
            (score / Path(name).name).write_bytes(tables[name])
        run_dir = reading / "run" / s
        run_dir.mkdir(parents=True)
        for src, dst in (("DONE", "DONE"), ("libraries.txt", "libraries.txt"), ("logs/timing.tsv", "timing.tsv"),
                         ("logs/plans.log", "plans.log")):
            shutil.copyfile(copy / src, run_dir / dst)
        outside = copy / "results" / "outside"
        if outside.is_dir():
            shutil.copytree(outside, run_dir / "outside")
        else:
            (run_dir / "outside").mkdir()
    files = sorted(p.relative_to(reading).as_posix() for p in reading.rglob("*") if p.is_file())
    with open(reading / "reading.txt", "w", encoding="utf-8") as handle:
        handle.write("run_checks.py gather: " + " ".join(["--withdrawn " + a for a in sorted(withdrawn or {})]).strip()
                     + "\n")
        for name in files:
            handle.write(f"{sha((reading / name).read_bytes())}  {name}\n")
    return rep


def count_command(reading, withdrawn):
    flags = " ".join(f"--withdrawn {a}" for a in sorted(withdrawn))
    return (f"python3 -I scripts/c7-outside/score.py --count {reading}/count {flags} {reading}/score/*/datasets.tsv"
            .replace("  ", " "))


def parse(argv):
    rest, withdrawn, ids_also = [], {}, []
    i = 0
    while i < len(argv):
        if argv[i] == "--withdrawn":
            if i + 4 > len(argv):
                raise SystemExit("--withdrawn takes an arm, a reason and a UTC time")
            arm, reason, time = argv[i + 1: i + 4]
            if arm not in ARMS or not reason.strip() or not TIME.fullmatch(time) or arm in withdrawn:
                raise SystemExit("--withdrawn takes a library arm once, a reason and a UTC time")
            withdrawn[arm] = (reason.strip(), time)
            i += 4
        elif argv[i] == "--ids-also":
            if i + 1 >= len(argv):
                raise SystemExit("--ids-also takes a path")
            ids_also.append(argv[i + 1])
            i += 2
        else:
            rest.append(argv[i])
            i += 1
    return rest, withdrawn, ids_also


def main(argv):
    if not argv or argv[0] not in ("check", "gather"):
        sys.stderr.write(__doc__)
        return 2
    rest, withdrawn, ids_also = parse(argv[1:])
    if argv[0] == "check" and len(rest) == 3:
        rep = run(rest[0], rest[1], rest[2], withdrawn, ids_also)
        print("\n".join(rep.lines))
        return int(rep.failed > 0)
    if argv[0] == "gather" and len(rest) == 4:
        rep = gather(rest[0], rest[1], rest[2], rest[3], withdrawn, ids_also)
        print("\n".join(rep.lines))
        if rep.failed:
            print("gather: nothing written, a check failed")
            return 1
        print(f"gather: wrote {rest[3]}; next: {count_command(rest[3], withdrawn)}")
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
