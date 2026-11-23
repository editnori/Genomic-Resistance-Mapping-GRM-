"""The dataset lists of plan block C4 (docs/screen-arms.md, Datasets).

  python3 -I lists.py <c3 lists.tsv> <out dir>

<c3 lists.tsv> is L as scripts/c3-inputs/lists.py wrote it
(scripts/c3-inputs/lists.tsv): the release reports with a lineage figure,
with their species and antibiotic in full, their genome count, M, S and the
first run (L less M). Every arm of C4 runs on every dataset of the first
run, so no other set is needed.

Writes <out dir>/lists.tsv (L's rows, the not-M datasets first, with the
default's floor max(3, n/200) at the dataset's genome count, whether the
floor arm moves it, which is when that floor is above 3, n of 800 or more,
and the arms it runs) and <out dir>/species.tsv (per species, in full: L,
L less M, of those in S, of those the floor arm moves, and of those in S),
and prints the totals. A fold file's floor is taken from the genomes it
tests or trains, which can be fewer than the dataset's; the floor arm's
study fits are run on every study dataset all the same.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ARMS = ["k31", "cmh-4000", "cmh-16000", "floor-3", "twin"]
FLOOR = 3
COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
           "in_M", "in_S", "first_run", "default_floor", "floor_moves", "arms"]
SPECIES_COLUMNS = ["species_slug", "species", "L", "L_not_M", "L_not_M_in_S", "L_not_M_floor_moves",
                   "L_not_M_floor_moves_in_S"]


def _c3():
    spec = importlib.util.spec_from_file_location("c3_lists_for_c4", HERE.parent / "c3-inputs" / "lists.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _c3()
read, write = C3.read, C3.write


def default_floor(samples):
    """grm's per-fit floor at a cohort of `samples` genomes (eval::support_floor)."""
    return max(FLOOR, samples // 200)


def lists(l_rows):
    out = []
    for r in l_rows:
        row = {k: r[k] for k in COLUMNS if k in r}
        floor = default_floor(int(r["samples"]))
        row["default_floor"] = floor
        row["floor_moves"] = int(floor != FLOOR)
        row["arms"] = ",".join(ARMS) if r["first_run"] == "1" else "-"
        out.append(row)
    # The first run (not M) first, then M; each part in L's order.
    return [r for r in out if r["first_run"] == "1"] + [r for r in out if r["first_run"] != "1"]


def species_rows(rows):
    out = []
    for sp in sorted({r["species_slug"] for r in rows}):
        rs = [r for r in rows if r["species_slug"] == sp]
        first = [r for r in rs if r["first_run"] == "1"]
        moves = [r for r in first if r["floor_moves"] == 1]
        out.append({
            "species_slug": sp, "species": rs[0]["species"], "L": len(rs), "L_not_M": len(first),
            "L_not_M_in_S": sum(r["in_S"] == "1" for r in first), "L_not_M_floor_moves": len(moves),
            "L_not_M_floor_moves_in_S": sum(r["in_S"] == "1" for r in moves),
        })
    return out


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    rows = lists(read(argv[0]))
    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    write(out / "lists.tsv", rows, COLUMNS)
    sp = species_rows(rows)
    write(out / "species.tsv", sp, SPECIES_COLUMNS)
    tot = {c: sum(r[c] for r in sp) for c in SPECIES_COLUMNS[2:]}
    print(f"L {tot['L']} in {len(sp)} species; L less M {tot['L_not_M']} in "
          f"{sum(r['L_not_M'] > 0 for r in sp)} species, {tot['L_not_M_in_S']} in S; the floor arm moves "
          f"{tot['L_not_M_floor_moves']} in {sum(r['L_not_M_floor_moves'] > 0 for r in sp)} species, "
          f"{tot['L_not_M_floor_moves_in_S']} in S")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
