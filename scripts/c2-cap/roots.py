"""The data roots of one species' cap arms (docs/cap-arms.md, Matrices).

  python3 -I roots.py plan <cohorts dir> <arms> <roots.json>
  python3 -I roots.py share <cap-all root> <root> <species slug> <linked list>
  python3 -I roots.py record <cap-all root> <root> <species slug> <linked list>
  python3 -I roots.py counts <roots dir> <species slug> <roots.json>

plan: the roots to build and the root each wider arm of <arms> (a comma
list) fits on. Every run builds cap-2000 and cap-all, whose structure.json
is the species grouping of every arm; cap-4000 gets its own root unless its
cohort equals cap-all's, and then, with cap-all among the arms, it is not
fitted and `same_as` names cap-all (Enterococcus faecium). Writes
<roots.json>: roots, fits ({arm: root}), same_as and each arm's genomes.

share and record: the cap-all root is the species pool of
scripts/input-comparison/run.py (H-1 gap G1). share gives the root a hard
link to each pool kset and a copy of its sidecar (Runner.share_ksets) and
writes the linked genomes to <linked list>, which stays on the box; record,
after grm count, holds every link to the pool's file and every sidecar to
the pool's bytes and writes the root's shared.json (Runner.record_shared).

counts: adds to <roots.json> per root its genomes, ksets, shared ksets,
lineages, duplicate groups of two or more and its matrix's bytes.
"""

import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT, MID, ALL = "cap-2000", "cap-4000", "cap-all"


def _run():
    spec = importlib.util.spec_from_file_location("c2_input_comparison_run",
                                                  HERE.parent / "input-comparison" / "run.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_ids(path):
    return [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines() if line.strip()]


def plan(cohorts, arms):
    """{roots, fits, same_as, genomes} for the comma list `arms`."""
    arms = [a for a in arms.split(",") if a]
    if not arms or any(a not in (MID, ALL) for a in arms) or len(set(arms)) != len(arms):
        raise ValueError(f"the arms are cap-4000, cap-all or both, each once, not {arms}")
    ids = {a: read_ids(Path(cohorts) / f"{a}.txt") for a in (DEFAULT, MID, ALL)}
    roots, fits, same_as = [ALL, DEFAULT], {}, {}
    for a in arms:
        if a == MID and ids[MID] == ids[ALL]:
            if ALL in arms:
                same_as[MID] = ALL
            else:
                fits[MID] = ALL
            continue
        if a == MID:
            roots.append(MID)
        fits[a] = a
    return {"roots": roots, "fits": fits, "same_as": same_as,
            "genomes": {a: len(ids[a]) for a in (DEFAULT, MID, ALL)}}


def pool_runner(pool_root, slug):
    run = _run()
    r = run.Runner.__new__(run.Runner)
    r.slug = slug
    r.pool_ksets = Path(pool_root) / "species" / slug / "ksets"
    return r


def counts(roots_dir, slug, doc):
    built = {}
    for name in doc["roots"]:
        d = Path(roots_dir) / name / "species" / slug
        ksets = list((d / "ksets").glob("*.kset")) if (d / "ksets").is_dir() else []
        s = json.loads((d / "structure.json").read_text(encoding="utf-8"))
        digests = {}
        for g in s["kset_digests"]:
            digests[g] = digests.get(g, 0) + 1
        shared = Path(roots_dir) / name / "shared.json"
        matrix = [p for p in d.rglob("*") if p.is_file() and "ksets" not in p.parts and p.name != "structure.json"]
        built[name] = {"genomes": len(read_ids(d / "genomes.txt")), "ksets": len(ksets),
                       "kset_bytes": sum(p.stat().st_size for p in ksets),
                       "shared": json.loads(shared.read_text(encoding="utf-8"))["shared"] if shared.is_file() else None,
                       "lineages": len(set(s["lineages"]["labels"])),
                       "duplicate_groups": sum(1 for n in digests.values() if n > 1),
                       "other_bytes": sum(p.stat().st_size for p in matrix)}
        if built[name]["genomes"] != len(s["genomes"]):
            raise ValueError(f"root {name}: structure.json holds {len(s['genomes'])} genomes, genomes.txt "
                             f"{built[name]['genomes']}")
    return built


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    what = argv[1]
    if what == "plan" and len(argv) == 5:
        try:
            doc = plan(argv[2], argv[3])
        except ValueError as e:
            raise SystemExit(str(e)) from None
        Path(argv[4]).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
        print(json.dumps(doc))
    elif what == "share" and len(argv) == 6:
        r = pool_runner(argv[2], argv[4])
        linked = r.share_ksets(Path(argv[3]))
        Path(argv[5]).write_text("".join(g + "\n" for g in linked), encoding="utf-8")
        print(f"ksets linked from the pool: {len(linked)}")
    elif what == "record" and len(argv) == 6:
        r = pool_runner(argv[2], argv[4])
        r.record_shared(Path(argv[3]), read_ids(argv[5]))
        print((Path(argv[3]) / "shared.json").read_text(encoding="utf-8"))
    elif what == "counts" and len(argv) == 5:
        doc = json.loads(Path(argv[4]).read_text(encoding="utf-8"))
        try:
            doc["built"] = counts(argv[2], argv[3], doc)
        except ValueError as e:
            raise SystemExit(str(e)) from None
        Path(argv[4]).write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8")
        print(json.dumps(doc["built"]))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv)
