"""The dataset list of plan block C5 (docs/settings-arms.md, Datasets).

  python3 -I lists.py <c3 lists.tsv> <out dir>

<c3 lists.tsv> is block C3's list (scripts/c3-inputs/lists.tsv): L's rows
with their species and antibiotic in full, M and S. C5 runs on S, the 69
datasets of docs/assumption-plan.md's Sets; its first run is S less M, and
those datasets come first in the list, by dataset id, then the datasets of S
in M, by dataset id.

Writes <out dir>/lists.tsv (S's rows, with first_run) and <out dir>/species.tsv
(per species, in full: S and S less M), and prints the totals. It reads and
writes the lists with block C3's code.
"""

import importlib.util
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
COLUMNS = ["dataset", "species_slug", "species", "antibiotic", "samples", "resistant", "susceptible", "lineages",
           "in_M", "first_run", "pairs_S+R", "labels_depend_on"]
SPECIES_COLUMNS = ["species_slug", "species", "S", "S_not_M"]


def _c3():
    spec = importlib.util.spec_from_file_location("c3_lists_for_c5", HERE.parent / "c3-inputs" / "lists.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C3 = _c3()


def lists(rows):
    """S's rows, S less M first; first_run is 1 for S less M."""
    s = [r for r in rows if r["in_S"] == "1"]
    out = []
    for r in sorted(s, key=lambda r: (r["in_M"] == "1", r["dataset"])):
        row = {k: r[k] for k in COLUMNS if k in r}
        row["first_run"] = int(r["in_M"] != "1")
        out.append(row)
    return out


def species_rows(rows):
    out = []
    for sp in sorted({r["species_slug"] for r in rows}):
        rs = [r for r in rows if r["species_slug"] == sp]
        out.append({"species_slug": sp, "species": rs[0]["species"], "S": len(rs),
                    "S_not_M": sum(r["first_run"] for r in rs)})
    return out


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    rows = lists(C3.read(argv[0]))
    out = Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    C3.write(out / "lists.tsv", rows, COLUMNS)
    sp = species_rows(rows)
    C3.write(out / "species.tsv", sp, SPECIES_COLUMNS)
    print(f"S {len(rows)} in {len(sp)} species; S less M {sum(r['first_run'] for r in rows)} in "
          f"{sum(r['S_not_M'] > 0 for r in sp)} species")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
