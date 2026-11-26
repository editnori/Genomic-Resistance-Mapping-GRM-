"""The cost of plan block C3's run on its datasets of M, in the plan's units, before it starts (docs/input-arms.md).

  python3 -I cost.py measured --tables-reused [<times table>]
  python3 -I cost.py measured --program-again [<times table>]
  python3 -I cost.py times <copies dir> <times table>
  python3 -I cost.py model

The plan's units: an xlarge-hour is one hour of a box of 16 CPUs and 32 GB,
16 core-hours, at $0.20 (C3's own budget: 45 xlarge-hours, 720 core-hours,
$9). The entry M run asks the estimate from the first run's measured times.

times writes those times as one table, <times table>: per species of M, each
step of its first-run box's logs/timing.tsv in whole seconds, rounded half
up, with the sha256 of that timing.tsv, which must be the one the copy's
SHA256SUMS gives it (<copies dir>/<species slug>/). The committed table is
reference/input-arms-m/first-run-times.tsv, which measured reads by default.

measured costs each species of M on one box of its own, so a box's wall
seconds are its xlarge-seconds:

- each lineage fit of the first run, scaled by the samples of the species'
  datasets of M over those of its first-run datasets in the arm's set (the
  program arms run on the datasets in A only);
- each study fit, scaled per arm by the samples of its datasets of M in S
  over those of its first-run datasets in S (and in A for the program
  arms). A species with datasets of M in S and no first-run study fit of an
  arm takes that arm's seed-42 lineage fit per sample of its first-run
  datasets, times its samples of M in S, and is named;
- its matrix builds and data steps once each, and the grm build twice (Same
  build's build of the first run's commit and the changed binary's);
- the program on each species with a dataset of M in A: with
  --tables-reused, the step that copies its tables from its first-run box
  (a species whose first run ran the program takes the largest such step of
  the others); with --program-again, the program run again: the measured
  program_env, and the measured program_run scaled by the species'
  prodigal step over that of the species that ran it;
- on the check species, its first-run fits again for Same results.

It names what it leaves out: m_list, Same counts, score.py, second.py, the
pack, the git archive copies, and box start-up and idle time.

model runs on the committed files alone and gives an upper bound: the entry
C3 arms' cost model of the first run (docs/assumption-plan-corrections.md,
C3 arms, Cost), per arm in xlarge-hours, central and upper. Each arm costs
its figure times its share:
the samples of the datasets of M in its set over those of the first run's
datasets in it (scripts/c3-inputs/lists.tsv), and for the three matrix arms,
whose figures hold their builds, the larger of that and the genomes of the
species of M over those of the first run's species
(reference/input-arms/page-builds.tsv), which bounds the model whatever
part of the figure the builds are. The program on every species of A takes
the genomes of the species of M in A over those of every species in A, or
nothing when its tables are reused. The data steps (the run's total less
its arm rows) take 9 of the 23 boxes, one more grm build a box takes the
45.05 s this box measured for a release build of 8bfd4ca at 16 jobs, and
the check species adds its first-run share of every arm for Same results.
Prints each part and the totals.
"""

import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
XLARGE_CPUS = 16
DOLLARS_PER_XLARGE_HOUR = 0.20
# The entry C3 arms, Cost: xlarge-hours of the first run, central and upper.
ARM_COST = {"k31": (12.12, 36.37), "twin": (10.01, 30.02), "k31-counts": (5.92, 17.40), "k31-aa": (6.74, 17.25),
            "unitigs": (7.42, 21.09), "program": (2.01, 6.02), "k31-program": (10.09, 30.28)}
PROGRAM_ON_SPECIES = (24.30, 36.45)
RUN_TOTAL = (82.50, 200.53)
FIRST_RUN_BOXES = 23
MATRIX_ARMS = ("k31-counts", "k31-aa", "unitigs")
PROGRAM_ARMS = ("program", "k31-program")
# A release build of grm at 8bfd4ca, 16 jobs, on the helper box of 26 Nov 2026 (Same build, run by hand).
BUILD_SECONDS = 45.05
CHECK_SPECIES = "neisseria_gonorrhoeae"


class Stop(Exception):
    pass


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def units(xh):
    return f"{xh:.2f} xlarge-hours, {xh * XLARGE_CPUS:.1f} core-hours, ${xh * DOLLARS_PER_XLARGE_HOUR:.2f}"


def in_set(row, arm, first):
    """Whether a lists.tsv row is in an arm's set: of M (first False) or of the first run (first True)."""
    if first:
        return row["first_run_A" if arm in PROGRAM_ARMS else "first_run"] == "1"
    return row["in_M"] == "1" and (arm not in PROGRAM_ARMS or row["in_A"] == "1")


def samples(rows, arm, first, slug=None):
    return sum(int(r["samples"]) for r in rows if in_set(r, arm, first) and (slug is None or r["species_slug"] == slug))


def genomes(builds):
    """Genomes per species, from page-builds.tsv."""
    out = {}
    for r in builds:
        out[r["species_slug"]] = max(out.get(r["species_slug"], 0), int(r["genomes"]))
    return out


# ---------- model ----------

def model(lists, builds, tables_reused=False):
    """[(part, central, upper)] from the committed cost model."""
    g = genomes(builds)
    m_species = sorted({r["species_slug"] for r in lists if r["in_M"] == "1"})
    first_species = sorted({r["species_slug"] for r in lists if r["first_run"] == "1"})
    if len(first_species) != FIRST_RUN_BOXES or set(first_species) != set(g):
        raise Stop(f"lists.tsv's first-run species ({len(first_species)}) are not page-builds.tsv's ({len(g)})")
    build_share = sum(g[s] for s in m_species) / sum(g[s] for s in first_species)
    parts = []
    for arm, (central, upper) in ARM_COST.items():
        list_arm = "k31" if arm == "twin" else arm
        share = samples(lists, list_arm, False) / samples(lists, list_arm, True)
        if arm in MATRIX_ARMS:
            share = max(share, build_share)
        parts.append((f"{arm}: share {share:.4f}", central * share, upper * share))
    a_species = sorted({r["species_slug"] for r in lists if r["in_A"] == "1"})
    m_a = sorted({r["species_slug"] for r in lists if r["in_M"] == "1" and r["in_A"] == "1"})
    p_share = 0.0 if tables_reused else sum(g[s] for s in m_a) / sum(g[s] for s in a_species)
    parts.append((f"the program on {len(m_a)} species of A: share {p_share:.4f}",
                  PROGRAM_ON_SPECIES[0] * p_share, PROGRAM_ON_SPECIES[1] * p_share))
    data = [t - sum(c[i] for c in ARM_COST.values()) - PROGRAM_ON_SPECIES[i] for i, t in enumerate(RUN_TOTAL)]
    boxes = len(m_species)
    parts.append((f"data steps on {boxes} of {FIRST_RUN_BOXES} boxes", data[0] * boxes / FIRST_RUN_BOXES,
                  data[1] * boxes / FIRST_RUN_BOXES))
    parts.append((f"one more grm build on each of {boxes} boxes", boxes * BUILD_SECONDS / 3600,
                  boxes * BUILD_SECONDS / 3600))
    check_c = check_u = 0.0
    for arm, (central, upper) in ARM_COST.items():
        list_arm = "k31" if arm == "twin" else arm
        share = samples(lists, list_arm, True, CHECK_SPECIES) / samples(lists, list_arm, True)
        check_c, check_u = check_c + central * share, check_u + upper * share
    parts.append((f"Same results: the first-run fits of {CHECK_SPECIES} again", check_c, check_u))
    return parts


# ---------- measured ----------

TIMES = REPO / "reference" / "input-arms-m" / "first-run-times.tsv"
TIMES_COLUMNS = ["species_slug", "step", "seconds", "source_sha256"]
STUDY_ARMS = ("k31", "unitigs", "k31-counts", "k31-aa", "program", "k31-program")
LEFT_OUT = ("m_list", "Same counts", "score.py", "second.py", "the pack", "the git archive copies",
            "box start-up and idle time")


def kind(step):
    """(kind, arm) of a timing.tsv step: lineage, study, build (a matrix build), grm, program or data."""
    if step.startswith("fit-") or step.startswith("study-"):
        rest = step.split("-", 1)[1]
        for arm in sorted(list(ARM_COST) + ["k31"], key=len, reverse=True):
            if rest == arm or rest.startswith(arm + "-"):
                return ("lineage" if step.startswith("fit-") else "study"), arm
        raise Stop(f"timing.tsv step {step}: no arm")
    if step in ("inputs_counts",):
        return "build", "k31-counts"
    if step in ("prodigal", "inputs_aa"):
        return "build", "k31-aa"
    if step in ("bcalm", "inputs_unitigs"):
        return "build", "unitigs"
    if step in ("program_env", "program_run", "program_tables"):
        return "program", None
    if step == "build":
        return "grm", None
    return "data", None


def held_timing(copy):
    """(timing.tsv's rows, its sha256), held to the copy's SHA256SUMS."""
    copy = Path(copy)
    sums = {}
    for line in (copy / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        sums[name.lstrip("*")] = digest
    path = copy / "logs" / "timing.tsv"
    if sums.get("logs/timing.tsv") != hashlib.sha256(path.read_bytes()).hexdigest():
        raise Stop(f"{copy.name}: logs/timing.tsv differs from its SHA256SUMS")
    return read_tsv(path), sums["logs/timing.tsv"]


def m_species(lists):
    return sorted({r["species_slug"] for r in lists if r["in_M"] == "1"})


def times_rows(copies, lists):
    """The times table's rows: every step of each species of M, in whole seconds, with its timing.tsv's sha256."""
    out = []
    for slug in m_species(lists):
        rows, digest = held_timing(Path(copies) / slug)
        out += [{"species_slug": slug, "step": r["step"], "seconds": int(float(r["wall_s"]) + 0.5),
                 "source_sha256": digest} for r in rows]
    return out


def write_times(path, rows):
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(TIMES_COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in TIMES_COLUMNS) + "\n")


def in_s(row, arm, first):
    """Whether a lists.tsv row is in an arm's study set: in S and in the arm's set of M or of the first run."""
    return row["in_S"] == "1" and in_set(row, arm, first)


def study_samples(rows, arm, first, slug):
    return sum(int(r["samples"]) for r in rows if in_s(r, arm, first) and r["species_slug"] == slug)


def by_species(times, lists):
    out = {}
    for r in times:
        out.setdefault(r["species_slug"], {})[r["step"]] = int(r["seconds"])
    if set(out) != set(m_species(lists)):
        raise Stop(f"the times table holds {len(out)} species, not the {len(m_species(lists))} species of M")
    if len({(r["species_slug"], r["source_sha256"]) for r in times}) != len(out):
        raise Stop("a species of the times table has more than one source sha256")
    return out


def program_seconds(secs, slug, program_again):
    """(seconds, note) of the program on one species with a dataset of M in A."""
    ran = sorted(s for s, x in secs.items() if "program_run" in x)
    if len(ran) != 1:
        raise Stop(f"{len(ran)} species of the times table ran the program; the estimate reads one")
    ran = ran[0]
    if program_again:
        if slug == ran:
            return secs[ran]["program_env"] + secs[ran]["program_run"], "measured"
        if "prodigal" not in secs[slug]:
            raise Stop(f"{slug}: no prodigal step to scale the program from")
        run = secs[ran]["program_run"] * secs[slug]["prodigal"] / secs[ran]["prodigal"]
        return secs[ran]["program_env"] + run, f"{ran}'s program_run times prodigal {secs[slug]['prodigal']} s over {secs[ran]['prodigal']} s"
    if "program_tables" in secs[slug]:
        return secs[slug]["program_tables"], "measured"
    largest = max(x["program_tables"] for x in secs.values() if "program_tables" in x)
    return largest, f"the largest program_tables step of the others, {largest} s"


def measured(times, lists, program_again=False):
    """([(part, xlarge-hours)] per species of M, [notes]) from the first run's step times, in whole seconds."""
    secs = by_species(times, lists)
    a_species = {r["species_slug"] for r in lists if r["in_M"] == "1" and r["in_A"] == "1"}
    parts, notes = [], []
    for slug in m_species(lists):
        lineage = study = build = data = grm = again = program = 0.0
        study_of = {}
        for step, s in secs[slug].items():
            k, arm = kind(step)
            if k == "lineage":
                list_arm = "k31" if arm == "twin" else arm
                first = samples(lists, list_arm, True, slug)
                m = samples(lists, list_arm, False, slug)
                if m and not first:
                    raise Stop(f"{slug}: {arm} has datasets of M and no first-run dataset to scale from")
                lineage += s * m / first if first else 0.0
                again += s if slug == CHECK_SPECIES else 0.0
            elif k == "study":
                study_of[arm] = study_of.get(arm, 0) + s
                again += s if slug == CHECK_SPECIES else 0.0
            elif k == "build":
                build += s
            elif k == "grm":
                grm += 2 * s
            elif k == "data":
                data += s
        lacking = []
        for arm in STUDY_ARMS:
            m_s, first_s = study_samples(lists, arm, False, slug), study_samples(lists, arm, True, slug)
            if not m_s:
                continue
            if first_s and arm in study_of:
                study += study_of[arm] * m_s / first_s
            else:
                first = samples(lists, arm, True, slug)
                if f"fit-{arm}-s42" not in secs[slug] or not first:
                    raise Stop(f"{slug}: {arm} has datasets of M in S and nothing to scale a study fit from")
                study += secs[slug][f"fit-{arm}-s42"] * m_s / first
                lacking.append(f"{arm} ({m_s} samples)")
        if lacking:
            notes.append(f"{slug}: no first-run study fit for its datasets of M in S, so each study fit is taken as "
                         f"its seed-42 lineage fit per sample: {', '.join(lacking)}")
        if slug in a_species:
            program, how = program_seconds(secs, slug, program_again)
            if how != "measured":
                notes.append(f"{slug}: the program, {how}")
        if not grm:
            raise Stop(f"{slug}: the times table has no build step")
        parts.append((f"{slug}: lineage fits {lineage / 3600:.2f}, study fits {study / 3600:.2f}, "
                      f"builds {build / 3600:.2f}, program {program / 3600:.2f}, data {data / 3600:.2f}, "
                      f"two grm builds {grm / 3600:.2f}, Same results {again / 3600:.2f}",
                      (lineage + study + build + program + data + grm + again) / 3600))
    return parts, notes


def main(argv):
    lists = read_tsv(HERE / "lists.tsv")
    try:
        if argv == ["model"]:
            parts = model(lists, read_tsv(REPO / "reference" / "input-arms" / "page-builds.tsv"))
            for name, c, u in parts:
                print(f"{name}: central {c:.2f}, upper {u:.2f} xlarge-hours")
            c, u = sum(p[1] for p in parts), sum(p[2] for p in parts)
            print(f"total, central: {units(c)}")
            print(f"total, upper: {units(u)}")
            reused = model(lists, read_tsv(REPO / "reference" / "input-arms" / "page-builds.tsv"), tables_reused=True)
            print(f"with the program's tables reused: central {units(sum(p[1] for p in reused))}; "
                  f"upper {units(sum(p[2] for p in reused))}")
            return 0
        if len(argv) in (2, 3) and argv[0] == "measured" and argv[1] in ("--tables-reused", "--program-again"):
            again = argv[1] == "--program-again"
            parts, notes = measured(read_tsv(argv[2] if len(argv) == 3 else TIMES), lists, program_again=again)
            for name, xh in parts:
                print(f"{name}: {xh:.2f} xlarge-hours")
            for note in notes:
                print(note)
            n_a = len({r["species_slug"] for r in lists if r["in_M"] == "1" and r["in_A"] == "1"})
            how = "the program run again" if again else "the program's tables reused"
            print(f"total, {how} on the {n_a} species with a dataset of M in A: {units(sum(p[1] for p in parts))}")
            print("left out: " + ", ".join(LEFT_OUT))
            return 0
        if len(argv) == 3 and argv[0] == "times":
            rows = times_rows(argv[1], lists)
            write_times(argv[2], rows)
            print(f"times: {len(rows)} steps of {len({r['species_slug'] for r in rows})} species, whole seconds")
            return 0
    except Stop as e:
        print(f"cost: stopped: {e}")
        return 1
    print("\n".join(__doc__.strip().splitlines()[2:6]), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
