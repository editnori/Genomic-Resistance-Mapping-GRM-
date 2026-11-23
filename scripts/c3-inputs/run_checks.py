"""The checks of plan block C3's run, before its count is taken (docs/input-arms-result.md, Method).

  python3 -I run_checks.py <back dir> <lists.tsv> <release dir> [--set-aside SPECIES ...] [--ids-also PATH ...]

<back dir> holds one copy per species, as the run's copies came back: each
<species>/ with results/, check/, logs/, expect.tsv, lists.txt, binary.txt,
DONE and SHA256SUMS, and score/ either unpacked or still in score.tar.
<lists.tsv> is scripts/c3-inputs/lists.tsv. <release dir> is the public
release unpacked, whose species/<species>.json lists each species' genome
ids. Prints one line per check, PASS or FAIL, with counts only, and exits 1
when one fails:

  species      the copies are exactly the species of lists.tsv's first run
  SHA256SUMS   every line of each species' SHA256SUMS holds, and DONE names
               its sha256; the files under results/ and check/ are its lines
  DONE         each DONE names its species, binary, code, result count and
               SHA256SUMS, as binary.txt does; species per commit and per
               code; no FAILED but those named with --set-aside, each named,
               set aside and older than its species' DONE
  expect.py    each expect.tsv is expect.py's from lists.tsv's first run,
               lists.txt names the design's seeds, twin seeds and arms, and
               the result files are exactly its rows
  refusals     each refusal record is an arm's study run refused for no
               candidate, and score/draws.tsv marks five constant-model rows,
               one per learner, for each, and no others
  genome ids   no genome id of the species' release lists in any file of the
               copies that is read or published, nor under --ids-also

The score/ files are read from score.tar in place when score/ is not
unpacked: hashed, counted for genome ids, and draws.tsv's `constant` column
counted. No score is printed. score/score.log is never opened, nor taken
out of score.tar: a SHA256SUMS line naming it fails. back.tar, the archive a
copy came in, is neither read nor published and is not checked here.
"""

import csv
import hashlib
import importlib.util
import io
import json
import re
import sys
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SEEDS, TWIN = "42,43,44", "45,46,47"
ARMS = ("program", "k31-program", "k31-counts", "k31-aa", "unitigs")
DEFAULT, TWIN_ARM = "k31", "twin"
REFUSED_FORMAT, NO_CANDIDATE, STUDY_SCHEME, STUDY_SEED = "grm-c3-refused-1", "no candidate", "study-whole", 42
CONSTANT_WHY = "constant model: refused for no candidate"
LEARNERS = 5
# Times as run_species.sh writes them (date -u +%FT%TZ): a FAILED or DONE with no such time is never set aside.
STAMP = r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ"
DONE_LINE = re.compile(rf"DONE (?P<species>\S+) (?P<time>{STAMP}) binary (?P<binary>[0-9a-f]{{64}}) "
                       r"code (?P<code>[0-9a-f]{64}) results (?P<results>\d+) SHA256SUMS (?P<sums>[0-9a-f]{64})")
FAILED_LINE = re.compile(rf"FAILED (?P<species>\S+) (?P<how>at|stopped at) (?P<step>\S+) (?P<time>{STAMP})")
# A genome id is any <digits>.<digits> not inside a longer number: after a letter or an underscore it still counts
# (genome_<id>, id<id>).
ID_TOKEN = re.compile(r"(?<![\d.])\d+\.\d+(?!\d)")
READ = ("results", "check", "logs")
READ_FILES = ("expect.tsv", "lists.txt", "binary.txt", "DONE", "SHA256SUMS")
SCORE_LOG = "score/score.log"


def _expect():
    spec = importlib.util.spec_from_file_location("c3_inputs_expect_for_checks", HERE / "expect.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


expect = _expect()


def load_tsv(text):
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def score_files(copy):
    """{name under score/: bytes} of one copy, score/score.log apart, which is never opened: from score/ when it is
    unpacked, else from score.tar in place."""
    folder = copy / "score"
    if folder.is_dir():
        out = {}
        for p in folder.rglob("*"):
            name = "score/" + p.relative_to(folder).as_posix()
            if p.is_file() and name != SCORE_LOG:
                out[name] = p.read_bytes()
        return out
    out = {}
    tar = copy / "score.tar"
    if not tar.is_file():
        return out
    with tarfile.open(tar) as archive:
        for member in archive:
            if not member.isfile():
                continue
            name = member.name[2:] if member.name.startswith("./") else member.name
            if not name.startswith("score/"):
                raise SystemExit(f"{copy.name}: score.tar holds a file outside score/")
            if name == SCORE_LOG:
                continue
            with archive.extractfile(member) as handle:
                out[name] = handle.read()
    return out


def lists_txt(text):
    return dict(line.split("\t", 1) for line in text.splitlines() if "\t" in line)


def expect_text(rows):
    lines = ["\t".join(expect.COLUMNS)] + ["\t".join(str(x) for x in r) for r in rows]
    return "\n".join(lines) + "\n"


class Report:
    def __init__(self):
        self.lines, self.failed, self.checks = [], 0, 0

    def check(self, name, good, detail, more=()):
        self.checks += 1
        self.failed += not good
        self.lines.append(f"{'PASS' if good else 'FAIL'} {name}: {detail}")
        self.lines += [f"  {m}" for m in more]


def run(back, l_path, release, set_aside=(), ids_also=()):
    back, release = Path(back), Path(release)
    with open(l_path, encoding="utf-8", newline="") as handle:
        lrows = list(csv.DictReader(handle, delimiter="\t"))
    first = [r for r in lrows if r["first_run"] == "1"]
    species = sorted({r["species_slug"] for r in first})
    rep = Report()

    # species
    got = sorted(p.name for p in back.iterdir() if p.is_dir())
    rep.check("species", got == species, f"{len(got)} copies, {len(species)} species in lists.tsv's first run",
              [f"missing {len(set(species) - set(got))}, more {len(set(got) - set(species))}"] if got != species else [])
    have = [s for s in species if (back / s).is_dir()]
    scores = {s: score_files(back / s) for s in have}

    # SHA256SUMS
    good, lines, bad_lines, listed_ok = 0, 0, 0, 0
    for s in have:
        copy = back / s
        sums = copy / "SHA256SUMS"
        if not sums.is_file():
            continue
        ok = True
        listed = set()
        for line in sums.read_text(encoding="utf-8").splitlines():
            want, name = line[:64], line[66:]
            listed.add(name)
            lines += 1
            data = scores[s].get(name) if name.startswith("score/") else (
                (copy / name).read_bytes() if (copy / name).is_file() else None)
            if data is None or sha(data) != want:
                bad_lines += 1
                ok = False
        done = (copy / "DONE").read_text(encoding="utf-8").strip() if (copy / "DONE").is_file() else ""
        m = DONE_LINE.fullmatch(done)
        ok &= m is not None and m["sums"] == sha(sums.read_bytes())
        on_disk = {p.relative_to(copy).as_posix() for top in ("results", "check") if (copy / top).is_dir()
                   for p in (copy / top).rglob("*") if p.is_file() and p.suffix in (".json", ".tsv", ".txt")}
        in_sums = {n for n in listed if n.split("/")[0] in ("results", "check")}
        listed_ok += on_disk == in_sums
        ok &= on_disk == in_sums
        good += ok
    rep.check("SHA256SUMS", good == len(species) and bad_lines == 0,
              f"{good} of {len(species)} species; {lines - bad_lines} of {lines} lines hold; results/ and check/ are "
              f"the listed files in {listed_ok} of {len(species)}")

    # DONE and binary.txt
    named, commits, codes, restarts, aside, failed = 0, {}, {}, 0, [], []
    for s in have:
        copy = back / s
        done = (copy / "DONE").read_text(encoding="utf-8").strip() if (copy / "DONE").is_file() else ""
        m = DONE_LINE.fullmatch(done)
        binary = {}
        if (copy / "binary.txt").is_file():
            for line in (copy / "binary.txt").read_text(encoding="utf-8").splitlines():
                parts = line.split()
                if len(parts) == 2 and parts[0] in ("commit", "grm", "code"):
                    binary[parts[0]] = parts[1]
                restarts += bool(parts) and parts[0] == "restart"
        results = sum(1 for _ in (copy / "results" / "picks").rglob("*.pick.json")) if (copy / "results").is_dir() else -1
        if m and m["species"] == s and (m["binary"], m["code"]) == (binary.get("grm"), binary.get("code")) \
                and int(m["results"]) == results:
            named += 1
        commits.setdefault(binary.get("commit", "none")[:7], []).append(s)
        codes.setdefault(binary.get("code", "none")[:12], []).append(s)
        if (copy / "FAILED").exists():
            f = FAILED_LINE.fullmatch((copy / "FAILED").read_text(encoding="utf-8").strip())
            if s in set_aside and f is not None and f["species"] == s and m is not None and f["time"] < m["time"]:
                aside.append(f"{s}: FAILED {f['how']} {f['step']} at {f['time']}, before DONE at {m['time']}: "
                             "named and set aside")
            else:
                failed.append(s)
    unknown = sorted(set(set_aside) - set(species))
    rep.check("DONE and binary.txt", named == len(species) and not failed and not unknown,
              f"{named} of {len(species)} DONE lines name their species, SHA256SUMS, result count, binary and code as "
              f"binary.txt does; FAILED in {len(failed)} species not set aside",
              [f"species per commit: {'; '.join(f'{c} {len(v)}' for c, v in sorted(commits.items()))}",
               f"species per code sha256: {'; '.join(f'{c} {len(v)}' for c, v in sorted(codes.items()))}",
               f"restart lines in binary.txt: {restarts}"] + aside
              + ([f"--set-aside names no species of the first run: {len(unknown)}"] if unknown else []))

    # expect.py
    ok_species, found, expected, more, datasets, study = 0, 0, 0, 0, 0, 0
    for s in have:
        copy = back / s
        mine = [r for r in first if r["species_slug"] == s]
        ds = [r["dataset"] for r in mine]
        a = [r["dataset"] for r in mine if r["first_run_A"] == "1"]
        st = [r["dataset"] for r in mine if r["in_S"] == "1"]
        datasets += len(ds)
        study += len(st)
        try:
            lt = lists_txt((copy / "lists.txt").read_text(encoding="utf-8"))
            arms = lt.get("arms", "").split(",")
            rows = expect.rows(s, ds, a, st, [int(x) for x in SEEDS.split(",")], [int(x) for x in TWIN.split(",")],
                               arms)
        except (OSError, ValueError):
            continue
        same_lists = (lt.get("species") == s and sorted(lt.get("datasets", "").split(",")) == sorted(ds)
                      and sorted(x for x in lt.get("a", "").split(",") if x) == sorted(a)
                      and lt.get("seeds") == SEEDS and lt.get("twin") == TWIN and sorted(arms) == sorted(ARMS))
        same_expect = (copy / "expect.tsv").is_file() and (copy / "expect.tsv").read_text(encoding="utf-8") == \
            expect_text(rows)
        picks = copy / "results" / "picks"
        listed, missing, both = set(), 0, 0
        for r in rows:
            name, inputs, arm, split = r[7], r[8], r[2], r[3]
            expected += 1
            refused = name[: -len(".pick.json")] + ".refused.json"
            # Only an arm's study run can be refused for no candidate (ruling 16).
            may_refuse = arm not in (DEFAULT, TWIN_ARM) and split == "study"
            here = (picks / name).is_file() or (may_refuse and (picks / refused).is_file())
            both += may_refuse and (picks / name).is_file() and (picks / refused).is_file()
            if inputs != "-":
                here = here and (picks / inputs).is_file()
            found += here
            missing += not here
            listed.add(name)
            if may_refuse:
                listed.add(refused)
        on_disk = {p.relative_to(picks).as_posix() for p in picks.rglob("*")
                   if p.name.endswith((".pick.json", ".refused.json"))} if picks.is_dir() else set()
        extra = len(on_disk - listed) + both
        more += extra
        ok_species += same_lists and same_expect and extra == 0 and missing == 0
    rep.check("expect.py", ok_species == len(species) and found == expected and more == 0,
              f"{ok_species} of {len(species)} species: lists.txt names the first run, the design's seeds, twin seeds "
              f"and arms, expect.tsv is expect.py's, and the result files are exactly its rows; {found} of {expected} "
              f"results found, {more} more",
              [f"datasets {datasets}, study datasets {study}, species {len(have)}"])

    # refusals and constant-model rows
    records, good_records, constant_rows, per_arm = 0, 0, 0, {}
    ok = True
    for s in have:
        picks = back / s / "results" / "picks"
        recs = sorted(picks.rglob("*.refused.json")) if picks.is_dir() else []
        mine = set()
        for p in recs:
            records += 1
            arm, name = p.parent.name, p.name
            dataset = name[: -len(f".{STUDY_SCHEME}-s{STUDY_SEED}.refused.json")]
            with open(p, encoding="utf-8") as handle:
                rec = json.load(handle)
            if (name.endswith(f".{STUDY_SCHEME}-s{STUDY_SEED}.refused.json") and arm in ARMS
                    and (rec.get("format"), rec.get("dataset"), rec.get("arm"), rec.get("scheme"), rec.get("seed"),
                         rec.get("reason")) == (REFUSED_FORMAT, dataset, arm, STUDY_SCHEME, STUDY_SEED, NO_CANDIDATE)):
                good_records += 1
            mine.add((dataset, arm))
            per_arm[arm] = per_arm.get(arm, 0) + 1
        draws = scores[s].get("score/draws.tsv")
        if draws is None:
            ok = False
            continue
        marked = {}
        for r in load_tsv(draws.decode("utf-8")):
            if r["constant"] == "1":
                constant_rows += 1
                key = (r["dataset"], r["arm"])
                marked[key] = marked.get(key, 0) + int(r["split"] == "study" and r["why"] == CONSTANT_WHY)
        ok &= set(marked) == mine and all(v == LEARNERS for v in marked.values())
    rep.check("refusals", ok and good_records == records and constant_rows == LEARNERS * records,
              f"{good_records} of {records} refusal records are an arm's study run refused for no candidate; "
              f"{constant_rows} constant-model rows in draws.tsv, {LEARNERS} per record",
              [f"records per arm: {'; '.join(f'{a} {n}' for a, n in sorted(per_arm.items())) or 'none'}"])

    # genome ids
    ids = set()
    lists_found = 0
    for s in species:
        path = release / "species" / f"{s}.json"
        if path.is_file():
            lists_found += 1
            with open(path, encoding="utf-8") as handle:
                ids |= {g["id"] for g in json.load(handle)["genomes"]}
    hits, files = 0, 0

    def count(data):
        return sum(1 for m in ID_TOKEN.finditer(data.decode("utf-8", errors="replace")) if m.group(0) in ids)

    for s in have:
        copy = back / s
        paths = [p for top in READ if (copy / top).is_dir() for p in (copy / top).rglob("*") if p.is_file()]
        paths += [copy / n for n in READ_FILES if (copy / n).is_file()]
        for p in paths:
            files += 1
            hits += count(p.read_bytes())
        for data in scores[s].values():
            files += 1
            hits += count(data)
    for extra in ids_also:
        extra = Path(extra)
        for p in ([extra] if extra.is_file() else sorted(q for q in extra.rglob("*") if q.is_file())):
            files += 1
            hits += count(p.read_bytes())
    rep.check("genome ids", hits == 0 and lists_found == len(species) and len(ids) > 0,
              f"{hits} genome ids in {files} files, against {len(ids)} release ids of {lists_found} species lists")
    rep.lines.append(f"CHECKS {rep.checks - rep.failed} of {rep.checks}")
    return rep


def main(argv):
    rest, set_aside, ids_also = [], [], []
    i = 0
    while i < len(argv):
        if argv[i] in ("--set-aside", "--ids-also"):
            if i + 1 >= len(argv):
                raise SystemExit(f"{argv[i]} takes a value")
            (set_aside if argv[i] == "--set-aside" else ids_also).append(argv[i + 1])
            i += 2
        else:
            rest.append(argv[i])
            i += 1
    if len(rest) != 3:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    rep = run(rest[0], rest[1], rest[2], set_aside, ids_also)
    print("\n".join(rep.lines))
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
