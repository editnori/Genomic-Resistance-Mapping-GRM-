"""Compares a data root's dataset settings with the release's reports, for
the probe of the MIC fold-plan experiment (docs/mic-plan.md).

  python3 -I release_specs.py <datasets.json> <release dir> <datasets> <out dir>

<datasets.json> is the data root's curated/datasets.json, <datasets> a
comma list of dataset ids. The input fingerprint hashes three settings of a
dataset: its antibiotic, its class minimum and its breakpoint. The
breakpoint is fitted at curation on every MIC pair of the species and drug,
so a curated table that passes even one more genome through QC moves its
`n` and `agreement`, and with them the fingerprint, while the dataset's
genomes, calls and MICs stay the same. The probe refits seed 42 on a copy of
the root whose breakpoints are the release's, and the scorer then needs the
release's fingerprint from it.

Writes to <out dir>:
  datasets.json      the root's list, with each listed dataset's breakpoint
                     set to its release report's;
  release-specs.tsv  dataset, breakpoint (equal or differs), and the root's
                     and the release's r_above and n;
  differ.txt         the listed datasets whose breakpoint differs, as a
                     comma list (empty when none).
Exits 1, writing nothing, when a listed dataset is missing from either side
or differs in its antibiotic or class minimum: no probe can stand in then.
"""

import json
import os
import sys

COLUMNS = ["dataset", "breakpoint", "root_r_above", "release_r_above", "root_n", "release_n"]


def field(bp, key):
    return "-" if bp is None else str(bp[key])


def compare(specs, release_dir, listed):
    """(patched list, tsv rows, differing ids); raises ValueError."""
    by_id = {d["id"]: d for d in specs}
    rows, differ, patched = [], [], [dict(d) for d in specs]
    index = {d["id"]: i for i, d in enumerate(patched)}
    for ds in listed:
        if ds not in by_id:
            raise ValueError(f"{ds}: not in the root's datasets.json")
        path = os.path.join(release_dir, "datasets", ds + ".json")
        if not os.path.isfile(path):
            raise ValueError(f"{ds}: the release has no report")
        with open(path) as f:
            report = json.load(f)
        root = by_id[ds]
        if root["antibiotic"] != report["antibiotic"] or root["min_class"] != report["counts"]["min_class"]:
            raise ValueError(f"{ds}: antibiotic or class minimum differs from the release's")
        mine, theirs = root.get("breakpoint"), report["breakpoint"]
        same = mine == theirs
        rows.append([ds, "equal" if same else "differs", field(mine, "r_above"), field(theirs, "r_above"),
                     field(mine, "n"), field(theirs, "n")])
        if not same:
            differ.append(ds)
            patched[index[ds]]["breakpoint"] = theirs
    return patched, rows, differ


def main(argv):
    if len(argv) != 4:
        print(__doc__.strip().splitlines()[3], file=sys.stderr)
        return 3
    specs_path, release_dir, listed, out = argv
    with open(specs_path) as f:
        specs = json.load(f)
    listed = [d for d in listed.split(",") if d and d != "-"]
    try:
        patched, rows, differ = compare(specs, release_dir, listed)
    except ValueError as e:
        print(e, file=sys.stderr)
        return 1
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "datasets.json"), "w") as f:
        json.dump(patched, f)
    with open(os.path.join(out, "release-specs.tsv"), "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    with open(os.path.join(out, "differ.txt"), "w") as f:
        f.write(",".join(differ))
    print(f"datasets whose breakpoint differs from the release's: {len(differ)} of {len(listed)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
