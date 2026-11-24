"""Score every species of plan block C5's run again from its copy, with this checkout's score.py.

  python3 -I rescore.py <copies dir> <lists.tsv> <release dir> <out dir>

For each species of the first run in <lists.tsv> it runs score.py on the
copy's expect.tsv and results/picks/ and the release the run read, writing
checks.tsv, draws.tsv and datasets.tsv into <out dir>/<slug>/. It never
opens, lists or copies anything under a copy's score/: of a finished
species it reads the SHA256SUMS lines of score/checks.tsv, score/draws.tsv
and score/datasets.tsv, and holds the three new tables to them byte for
byte before any table is read. A species with no SHA256SUMS (one that
stopped at its score step) is scored and not held.

Writes <out dir>/rescore.txt: per species a PASS or FAIL line with
score.py's exit status and the tables held, and last `SAME <n> of <held>`.
Exits 1 unless every species was scored and every held table is the same.
"""

import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TABLES = ("checks.tsv", "draws.tsv", "datasets.tsv")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


checks = _load("c5_settings_run_checks_for_rescore", HERE / "run_checks.py")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sums(path):
    """{name: sha256} of a SHA256SUMS file: the copy's list, not a file under score/."""
    out = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        sha, _, name = line.partition("  ")
        out[name] = sha
    return out


def rescore(copies, lists_path, release, out):
    copies, release, out = Path(copies), Path(release), Path(out)
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"{out} is not empty")
    species = sorted({r["species_slug"] for r in checks.first_run(checks.load(lists_path))})
    lines, same, held = [], 0, 0
    for slug in species:
        copy = copies / slug
        dest = out / slug
        dest.mkdir(parents=True)
        run = subprocess.run([sys.executable, "-I", "-B", str(HERE / "score.py"), str(copy / "expect.tsv"),
                              str(copy / "results" / "picks"), str(release), str(dest)],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        # score.py prints its check counts alone; keep them beside its tables.
        (dest / "score-output.txt").write_text(run.stdout, encoding="utf-8")
        if run.returncode != 0 or not all((dest / t).is_file() for t in TABLES):
            lines.append(f"FAIL {slug}: score.py ended {run.returncode}; "
                         f"{sum((dest / t).is_file() for t in TABLES)} of 3 tables written")
            continue
        if not (copy / "SHA256SUMS").is_file():
            lines.append(f"PASS {slug}: score.py ended 0; 3 tables, not held: the copy has no SHA256SUMS")
            continue
        listed = sums(copy / "SHA256SUMS")
        agree = [t for t in TABLES if listed.get(f"score/{t}") == sha256(dest / t)]
        held += len(TABLES)
        same += len(agree)
        lines.append(f"{'PASS' if len(agree) == len(TABLES) else 'FAIL'} {slug}: score.py ended 0; {len(agree)} of "
                     f"3 tables the same as its SHA256SUMS lines"
                     + (f"; differs: {', '.join(t for t in TABLES if t not in agree)}" if len(agree) < 3 else ""))
    lines.append(f"SAME {same} of {held}")
    (out / "rescore.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return lines


def main(argv):
    if len(argv) != 4:
        sys.stderr.write(__doc__)
        return 2
    lines = rescore(*argv)
    print("\n".join(lines))
    return 0 if checks.rescored(argv[3]) else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
