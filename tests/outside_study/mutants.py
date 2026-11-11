"""Breaks the kit on purpose and records that the suite catches each break.

Run from the repository root:

    python3 -I tests/outside_study/mutants.py [--out results.tsv]

Each mutant is a copy of scripts/outside-study, tests/outside_study and the
plan in a temporary directory with one change made; the suite runs there
with `python3 -I -B -m unittest discover -s tests/outside_study -t .`. A
mutant is caught when the suite fails. The copy's FROZEN.sha256 is rewritten
to match the changed file first, so the test that compares the pin file with
the files is never what catches a break. The unchanged copy must pass first.
Not a test module: unittest does not discover it.
"""

import argparse
import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PLAN = "docs/outside-study-test.md"

# (name, file, [(old, new), ...]); each old string must occur exactly once.
MUTANTS = [
    ("TP and FN swapped", "score.py",
     [('c["tp" if call == "resistant" else "fn"] += 1', 'c["fn" if call == "resistant" else "tp"] += 1')]),
    ("the two populations pooled", "score.py",
     [('c["population"] == want and', 'c["population"] in ("primary (other-study)", "same-study") and')]),
    ("an unlabelled (conflict) pair counted as R", "score.py",
     [('        elif p[flag] == "yes":\n            refuse(f"pairs: {key} has call {p[\'call\']!r} yet is flagged {flag}")',
       '        elif p["call"] == "conflict":\n            out[key] = True')]),
    ("the random split's reference read in place of the lineage split's", "score.py",
     [('ref = row["lineage_" + m]', 'ref = row["random_" + m]')]),
    ("interval bounds swapped", "score.py",
     [('out[f"{name}_lo"] = None if ci is None else ci[0]', 'out[f"{name}_lo"] = None if ci is None else ci[1]'),
      ('out[f"{name}_hi"] = None if ci is None else ci[1]', 'out[f"{name}_hi"] = None if ci is None else ci[0]')]),
    ("an excluded genome scored", "score.py",
     [('survivors = [g for g in ids if g not in excluded]', 'survivors = list(ids)')]),
    ("an excluded genome's prediction row accepted", "score.py",
     [('            if g in excluded:\n                refuse(f"{d}: excluded genome',
       '            if False:\n                refuse(f"{d}: excluded genome')]),
    ("the reading's comparison reversed", "score.py",
     [('    if reference > hi:\n        return "below"\n    if reference < lo:\n        return "above"',
       '    if reference > hi:\n        return "above"\n    if reference < lo:\n        return "below"')]),
    ("the row's report hash not checked", "score.py",
     [('    if row.get("report_sha256") != sha:', '    if False:')]),
    ("the report pin not checked", "score.py",
     [('    if sha != want:\n        refuse(f"{dataset}: report sha256', '    if False:\n        refuse(f"{dataset}: report sha256')]),
    ("the frozen tables' pins not checked", "score.py",
     [('            if got != tpins[name]:', '            if False:')]),
    ("the frozen gate skips git", "score.py",
     [('    if git("diff", "--quiet", "HEAD", "--", *tracked).returncode != 0:', '    if False:')]),
    ("25 of a class no longer enough", "score.py",
     [('reportable = n_r >= MIN_CLASS and n_s >= MIN_CLASS', 'reportable = n_r > MIN_CLASS and n_s > MIN_CLASS')]),
    ("leave-one-study-out keeps only the study left out", "score.py",
     [('rest = [g for g in scored if group_of(g, "bioproject") != b]',
       'rest = [g for g in scored if group_of(g, "bioproject") == b]')]),
    ("exclusion order: species check before sequence", "score.py",
     [('        if a["fasta"] != "yes":\n            out[g] = ("no sequence", a["detail"])\n            continue\n',
       ''),
      ('        if v["verdict"] == "fail":',
       '        if a["fasta"] != "yes":\n            out[g] = ("no sequence", a["detail"])\n            continue\n'
       '        if v["verdict"] == "fail":')]),
    ("a call not recomputed from its features", "score.py",
     [('    if row.get("call") != call:', '    if False:')]),
    ("the known-gene table not replayed", "score.py",
     [('        if rows != cat["groups"]:', '        if False:')]),
    ("run.py: a changed earlier file not noticed", "run.py",
     [('                    if self.digest(key) != sha:', '                    if False:')]),
    ("run.py: a cohort genome accepted as a frozen id", "run.py",
     [('            inside = [g for g in ids if g in cohort]', '            inside = []')]),
    ("run.py: a changed recorded species-check verdict not noticed", "run.py",
     [('        if changes:\n            p = cdir', '        if False:\n            p = cdir')]),
    ("run.py: a recorded verdict that moves to or from fail read as a GC move", "run.py",
     [('                and {g["verdict"], n["verdict"]} == {"pass", "borderline"}\n',
       '                and True\n')]),
    ("run.py: a recorded flag that moves with another distance read as a GC move", "run.py",
     [('                and all(k in g and k in n and g[k] == n[k] for k in SAME_DISTANCES)\n',
       '                and True\n')]),
    ("run.py: GC moves of recorded genomes not listed", "run.py",
     [('        S.write_tsv(moves_path, ["genome_id", "species", "recorded", "now"], moves)',
       '        S.write_tsv(moves_path, ["genome_id", "species", "recorded", "now"], [])')]),
]


def copy_tree(dst):
    shutil.copytree(ROOT / "scripts/outside-study", dst / "scripts/outside-study",
                    ignore=shutil.ignore_patterns("__pycache__"))
    shutil.copytree(ROOT / "tests/outside_study", dst / "tests/outside_study",
                    ignore=shutil.ignore_patterns("__pycache__"))
    (dst / "docs").mkdir()
    shutil.copy(ROOT / PLAN, dst / PLAN)


def repin(top):
    """The copy's pin file, rewritten for the files as they now are."""
    pins = top / "scripts/outside-study/FROZEN.sha256"
    names = [line.split(None, 1)[1].strip() for line in pins.read_text().splitlines() if line.strip()]
    pins.write_text("".join(f"{hashlib.sha256((top / n).read_bytes()).hexdigest()}  {n}\n" for n in names))


def run_suite(top):
    p = subprocess.run([sys.executable, "-I", "-B", "-m", "unittest", "discover", "-s", "tests/outside_study",
                        "-t", "."], cwd=top, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "PYTHONHASHSEED": "0"})
    failed = sorted({line.split(" ", 1)[1].split(" ")[0] for line in p.stderr.splitlines()
                     if line.startswith(("FAIL: ", "ERROR: "))})
    return p.returncode, p.stderr.strip().splitlines()[-1], failed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out")
    args = ap.parse_args(argv)
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "unchanged"
        copy_tree(base)
        code, last, _ = run_suite(base)
        if code != 0:
            print(f"the unchanged copy fails: {last}", file=sys.stderr)
            return 1
        print(f"unchanged copy: {last}", file=sys.stderr)
        for i, (name, fname, edits) in enumerate(MUTANTS, start=1):
            top = Path(tmp) / f"m{i}"
            copy_tree(top)
            path = top / "scripts/outside-study" / fname
            text = path.read_text()
            for old, new in edits:
                if text.count(old) != 1:
                    print(f"mutant {i} ({name}): {old!r} occurs {text.count(old)} times", file=sys.stderr)
                    return 1
                text = text.replace(old, new)
            path.write_text(text)
            repin(top)
            code, last, failed = run_suite(top)
            caught = code != 0
            rows.append((i, name, fname, "caught" if caught else "NOT CAUGHT", last, len(failed),
                         ", ".join(failed[:4]) + (" ..." if len(failed) > 4 else "")))
            print(f"{i:2d} {'caught' if caught else 'NOT CAUGHT':10s} {name}: {last}", file=sys.stderr)
            shutil.rmtree(top)
    lines = ["mutant\tchange\tfile\tresult\tsuite_last_line\tfailing_tests\tfirst_failing_tests"]
    lines += ["\t".join(str(x) for x in r) for r in rows]
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    else:
        sys.stdout.write(text)
    return 0 if all(r[3] == "caught" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
