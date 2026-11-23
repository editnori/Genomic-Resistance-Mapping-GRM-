"""Fit times of a C4 run, from its inputs files (docs/screen-arms.md).

  python3 -I fit_times.py <picks dir> <out.tsv>

Reads every <arm>/<dataset>.<tag>-s<seed>.inputs.json under <picks dir>
and writes <out.tsv>: per arm, tag, seed, dataset, split, outer fold and
learner, the fit's wall and CPU seconds and the process's peak resident
size after it, in KB. Times and counts only: an inputs file holds no
genome id, feature or score. Prints the row count.
"""

import json
import sys
from pathlib import Path

COLUMNS = ["arm", "rule", "budget", "min_support", "tag", "seed", "dataset", "split", "fold", "learner", "seconds",
           "cpu_seconds", "peak_kb"]


def rows(picks):
    out = []
    for path in sorted(Path(picks).glob("*/*.inputs.json")):
        with open(path) as f:
            doc = json.load(f)
        for t in doc["fit_times"]:
            out.append({"arm": path.parent.name, "rule": doc["rule"], "budget": doc["budget"],
                        "min_support": doc["min_support"], "tag": doc["tag"], "seed": doc["seed"],
                        "dataset": doc["dataset"], "split": t["split"], "fold": t["fold"], "learner": t["learner"],
                        "seconds": t["seconds"], "cpu_seconds": t["cpu_seconds"], "peak_kb": t["peak_kb"]})
    return out


def main(argv):
    if len(argv) != 2:
        print(__doc__.strip().splitlines()[2], file=sys.stderr)
        return 2
    out = rows(argv[0])
    with open(argv[1], "w", encoding="utf-8") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in out:
            f.write("\t".join(str(r[c]) for c in COLUMNS) + "\n")
    print(f"fit times: {len(out)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
