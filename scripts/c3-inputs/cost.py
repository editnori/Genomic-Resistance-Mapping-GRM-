"""The cost of plan block C3's run on its datasets of M, in the plan's units, before it starts (docs/input-arms.md).

  python3 -I cost.py model
  python3 -I cost.py measured <copies dir>

The plan's units: an xlarge-hour is one hour of a box of 16 CPUs and 32 GB,
16 core-hours, at $0.20 (C3's own budget: 45 xlarge-hours, 720 core-hours,
$9). The entry M run asks the estimate from the first run's measured times.

measured reads them: <copies dir>/<species slug>/logs/timing.tsv of each
first-run species' copy (step, wall seconds, peak KB), held to the sha256
that copy's SHA256SUMS gives it. Each species of M costs what its first-run
box measured, with its fits scaled by the samples of its datasets of M over
those of its first-run datasets in each arm's set (the program arms run on
the datasets in A only); its data steps, builds and program steps once
each; one more grm build for Same build; and on the check species, its
first-run fits again for Same results. One box a species, so a box's wall
seconds are its xlarge-seconds.

model runs on the committed files alone: those times were never copied into
reference/input-arms/, so it rests on the entry C3 arms' cost model of the
first run (docs/assumption-plan-corrections.md, C3 arms, Cost), per arm in
xlarge-hours, central and upper. Each arm costs its figure times its share:
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

def kind(step):
    """(kind, arm) of a timing.tsv step: fit, build, program or data."""
    if step.startswith("fit-") or step.startswith("study-"):
        rest = step.split("-", 1)[1]
        for arm in sorted(list(ARM_COST) + ["k31"], key=len, reverse=True):
            if rest == arm or rest.startswith(arm + "-"):
                return "fit", arm
        raise Stop(f"timing.tsv step {step}: no arm")
    if step in ("inputs_counts",):
        return "build", "k31-counts"
    if step in ("prodigal", "inputs_aa"):
        return "build", "k31-aa"
    if step in ("bcalm", "inputs_unitigs"):
        return "build", "unitigs"
    if step in ("program_env", "program_run", "program_tables"):
        return "program", None
    return "data", None


def held_timing(copy):
    copy = Path(copy)
    sums = {}
    for line in (copy / "SHA256SUMS").read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        sums[name.lstrip("*")] = digest
    path = copy / "logs" / "timing.tsv"
    if sums.get("logs/timing.tsv") != hashlib.sha256(path.read_bytes()).hexdigest():
        raise Stop(f"{copy.name}: logs/timing.tsv differs from its SHA256SUMS")
    return read_tsv(path)


def measured(copies, lists, tables_reused=False):
    """[(part, xlarge-hours)] from the copies' measured times."""
    m_species = sorted({r["species_slug"] for r in lists if r["in_M"] == "1"})
    parts = []
    for slug in m_species:
        rows = held_timing(Path(copies) / slug)
        fit = build = program = data = again = 0.0
        built = 0.0
        for r in rows:
            k, arm = kind(r["step"])
            wall = float(r["wall_s"])
            if k == "fit":
                list_arm = "k31" if arm == "twin" else arm
                first = samples(lists, list_arm, True, slug)
                m = samples(lists, list_arm, False, slug)
                if m and not first:
                    raise Stop(f"{slug}: {arm} has datasets of M and no first-run dataset to scale from")
                fit += wall * m / first if first else 0.0
                again += wall if slug == CHECK_SPECIES else 0.0
            elif k == "build":
                build += wall
            elif k == "program":
                program += 0.0 if tables_reused else wall
            else:
                data += wall
                built += wall if r["step"] == "build" else 0.0
        if not built:
            raise Stop(f"{slug}: timing.tsv has no build step")
        parts.append((f"{slug}: fits {fit / 3600:.2f}, builds {build / 3600:.2f}, program {program / 3600:.2f}, "
                      f"data {data / 3600:.2f}, Same build {built / 3600:.2f}, Same results {again / 3600:.2f}",
                      (fit + build + program + data + built + again) / 3600))
    return parts


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
        if len(argv) == 2 and argv[0] == "measured":
            parts = measured(argv[1], lists)
            for name, xh in parts:
                print(f"{name}: {xh:.2f} xlarge-hours")
            print(f"total: {units(sum(p[1] for p in parts))}")
            return 0
    except Stop as e:
        print(f"cost: stopped: {e}")
        return 1
    print("\n".join(__doc__.strip().splitlines()[2:4]), file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
