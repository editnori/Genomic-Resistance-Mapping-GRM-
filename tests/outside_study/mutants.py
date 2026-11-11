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
import concurrent.futures
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
    ("the frozen tables not checked against tables.sha256", "score.py",
     [('            if got != tpins[name]:', '            if False:')]),
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
    # The independent review's breaks ("review 1" to "review 16" and "review
    # 1b"). Breaks 7, 12 and 16 are the same faults written against the code
    # as it now is.
    ("review 1: a genome's label taken from the first dataset it is labelled on", "score.py",
     [('            out[key] = p["call"] == "R"',
       '            out[key] = out.setdefault(("*", p["genome_id"]), p["call"] == "R")')]),
    ("review 1b: a genome's first label reused on every dataset", "score.py",
     [('    out, seen = {}, set()\n    for p in pairs:', '    out, seen, first = {}, set(), {}\n    for p in pairs:'),
      ('            out[key] = p["call"] == "R"', '            out[key] = first.setdefault(p["genome_id"], p["call"] == "R")')]),
    ("review 2: the lasso scored with the SCM's calls", "score.py",
     [('            c = confusion((lab[(d, g)], calls[(d, learner, g)]) for g in scored)',
       '            c = confusion((lab[(d, g)], calls[(d, "scm" if learner == "lasso" else learner, g)]) for g in scored)')]),
    ("review 3: known-gene calls read from the first dataset's record on every dataset", "score.py",
     [('            row = kg["calls"].get((d, g))', '            row = kg["calls"].get((datasets[0], g))')]),
    ("review 4: counts only decided on the R genomes labelled before the exclusions", "score.py",
     [('reportable = n_r >= MIN_CLASS and n_s >= MIN_CLASS',
       'reportable = sum(lab[(d, g)] for g in labelled) >= MIN_CLASS and n_s >= MIN_CLASS')]),
    ("review 5: a leave-one-study-out fold needs 25 R only", "score.py",
     [('"reportable": rr >= MIN_CLASS and rs >= MIN_CLASS, **rc}', '"reportable": rr >= MIN_CLASS, **rc}')]),
    ("review 6: a no-call counted as a susceptible call", "score.py",
     [('        if call is None:\n            c["no_call_resistant" if resistant else "no_call_susceptible"] += 1\n        elif resistant:',
       '        if resistant:')]),
    ("review 7: the gate reads the frozen files from the index, not from HEAD", "score.py",
     [('    blob = git("cat-file", "blob", f"HEAD:{rel}")', '    blob = git("cat-file", "blob", f":{rel}")')]),
    ("review 8: a row's k not checked", "score.py",
     [('    if row.get("k") != report["counts"]["k"] or row.get("input_fingerprint")',
       '    if row.get("input_fingerprint")')]),
    ("review 9: run.py predict: the final-model replay on the cohort never refuses", "run.py",
     [('        if problems:\n            p = rdir / "replay-problems.txt"', '        if False:\n            p = rdir / "replay-problems.txt"')]),
    ("review 10: run.py known-gene: the cohort replay never refuses", "run.py",
     [('        if problems:\n            p = kdir / "replay-problems.txt"', '        if False:\n            p = kdir / "replay-problems.txt"')]),
    ("review 11: run.py known-gene: gene calls from another catalog accepted", "run.py",
     [('            if g["catalog_version"] != want_version:', '            if False:')]),
    ("review 12: run.py: the stages after setup skip the frozen gate", "run.py",
     [('        return S.frozen_gate(self.args.expect_commit, self.rehearsal)',
       '        if self.receipt_path("setup").exists():\n'
       '            return json.loads(self.receipt_path("setup").read_text())["kit"]\n'
       '        return S.frozen_gate(self.args.expect_commit, self.rehearsal)')]),
    ("review 13: SCM disjunction evaluated as a conjunction", "score.py",
     [('            resistant = any(holds)', '            resistant = all(holds)')]),
    ("review 14: CART read past cart_depth", "score.py",
     [('            if s is not None and node["depth"] < depth:\n                i = s["present"]',
       '            if s is not None:\n                i = s["present"]')]),
    ("review 15: the reading compares the estimate, not the interval, with the report", "score.py",
     [('                    row[f"reading_{m}"] = reading(ref, lo, hi)',
       '                    row[f"reading_{m}"] = reading(ref, row[m], row[m])')]),
    ("review 16: run.py: a recorded verdict moving between pass and borderline not noticed", "run.py",
     [('            if n is not None and n["verdict"] == g["verdict"]:',
       '            if n is not None and (n["verdict"] == "fail") == (g["verdict"] == "fail"):')]),
    # The frozen gate, one check at a time.
    ("gate 1: git follows replace objects", "score.py",
     [('"-C", str(KIT_ROOT), "--no-replace-objects", *args]', '"-C", str(KIT_ROOT), *args]')]),
    ("gate 2: git inherits the caller's environment", "score.py",
     [('capture_output=True, env=git_env())', 'capture_output=True, env=dict(os.environ, **git_env()))')]),
    ("gate 3: git looked up on the caller's PATH", "score.py",
     [('    for path in GIT_PATHS:\n        if os.path.isfile(path) and os.access(path, os.X_OK):\n            return path',
       '    for path in GIT_PATHS:\n        if os.path.isfile(path) and os.access(path, os.X_OK):\n            return "git"'),
      ('    return {"PATH": "/usr/bin:/bin", "HOME"', '    return {"PATH": os.environ["PATH"], "HOME"')]),
    ("gate 4: FROZEN.sha256 on disk not compared with HEAD's", "score.py",
     [('    if committed_bytes(FROZEN_REL) != frozen_disk:', '    if False:')]),
    ("gate 5: a frozen file's committed bytes not checked against its pin", "score.py",
     [('        got = sha256_bytes(committed_bytes(rel))\n        if got != pins[rel]:',
       '        got = sha256_bytes(committed_bytes(rel))\n        if False:')]),
    ("gate 6: a frozen file on disk not checked against its pin", "score.py",
     [('        got = sha256_bytes(work_tree_bytes(rel))\n        if got != pins[rel]:',
       '        got = sha256_bytes(work_tree_bytes(rel))\n        if False:')]),
    ("gate 7: HEAD not compared with the expected commit", "score.py",
     [('    if expect_commit is not None and commit != expect_commit:', '    if False:')]),
    ("gate 8: a real run goes on without --expect-commit", "score.py",
     [('    if expect_commit is None and not rehearsal:', '    if False:')]),
    ("gate 9: the top of git's work tree not checked", "score.py",
     [('    if Path(top.stdout.decode().strip()).resolve() != KIT_ROOT:', '    if False:')]),
    ("gate 10: symbolic links followed", "score.py",
     [('        if path.is_symlink():', '        if False:')]),
    ("gate 11: the running score.py not checked", "score.py",
     [('    if Path(path).resolve() != KIT_ROOT / rel:', '    if False:')]),
    ("gate 12: run.py: the loaded score.py not checked", "run.py",
     [('    if Path(S.__file__).resolve() != KIT / "score.py":', '    if False:')]),
    ("gate 13: a pin for a file the gate does not know accepted", "score.py",
     [('        if rel not in FROZEN_FILES:\n            refuse', '        if False:\n            refuse')]),
    ("gate 14: score.py takes --tables in a real run", "score.py",
     [('        if args.tables is not None and not args.rehearsal:', '        if False:'),
      ('        if not rehearsal:\n            refuse("--tables is for', '        if False:\n            refuse("--tables is for')]),
    ("gate 15: run.py takes --tables, --ids and --datasets in a real run", "run.py",
     [('        if not args.rehearsal:\n            for opt in', '        if False:\n            for opt in')]),
    # Receipts and stages.
    ("run.py: a stage runs again over its own receipt", "run.py",
     [('        if self.receipt_path(stage).exists():', '        if False:')]),
    ("run.py: a non-empty scratch root accepted", "run.py",
     [('        if self.root.exists() and any(self.root.iterdir()):', '        if False:')]),
    ("run.py: any fetch failure read as an answer of no sequence", "run.py",
     [('            elif g in failed and no_sequence_answer(failed[g]):', '            elif g in failed:')]),
    ("run.py: attempts not logged", "run.py",
     [('        with open(p, "a") as f:\n            f.write(json.dumps(line, sort_keys=True) + "\\n")',
       '        return')]),
    ("score.py: a receipt from another commit accepted", "score.py",
     [('        if r["kit"]["commit"] != gate["commit"]:', '        if False:')]),
    ("score.py: a work file changed after its stage accepted", "score.py",
     [('        if got != want:\n            refuse(f"work/{rel} has sha256', '        if False:\n            refuse(f"work/{rel} has sha256')]),
    ("score.py: a receipt for other arguments accepted", "score.py",
     [('            if r["args"][key] != want:', '            if False:')]),
    ("score.py: predicted genomes not checked against the survivors", "score.py",
     [('    if inputs["predict_survivors"] != survivors:', '    if False:')]),
    # Outputs.
    ("an MLST value with ST '-' read as an ST", "score.py",
     [('in ("", "-")', 'in ("",)')]),
    ("no BioProject recorded counted as a study", "score.py",
     [('            projects = sorted({group_of(g, "bioproject") for g in scored} - {NO_BIOPROJECT})',
       '            projects = sorted({group_of(g, "bioproject") for g in scored})')]),
    ("duplicates of cohort genomes not listed", "score.py",
     [('        by_digest[dg].append((g, "cohort"))', '        pass')]),
    ("paired: learner-only and known-gene-only swapped", "score.py",
     [('                    learner_only += a and not k\n                    known_only += k and not a',
       '                    learner_only += k and not a\n                    known_only += a and not k')]),
    ("the known-gene reference marked as cross-validated", "score.py",
     [('REFERENCE_KIND = {KNOWN: "fixed rule, whole cohort"}', 'REFERENCE_KIND = {}')]),
    ("a rehearsal's outputs without the prefix", "score.py",
     [('    return ("rehearsal-" if res["rehearsal"] else "") + name', '    return name')]),
    ("the reportable leave-one-study-out folds not counted", "score.py",
     [('            folds += rest_r >= MIN_CLASS and rest_s >= MIN_CLASS', '            folds += 1')]),
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
                        "-t", "."], cwd=top, capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"})
    failed = sorted({line.split(" ", 1)[1].split(" ")[0] for line in p.stderr.splitlines()
                     if line.startswith(("FAIL: ", "ERROR: "))})
    return p.returncode, p.stderr.strip().splitlines()[-1], failed


def apply(top, fname, edits):
    path = top / "scripts/outside-study" / fname
    text = path.read_text()
    for old, new in edits:
        if text.count(old) != 1:
            raise SystemExit(f"{fname}: {old!r} occurs {text.count(old)} times")
        text = text.replace(old, new)
    path.write_text(text)


def one(tmp, i, mutant):
    name, fname, edits = mutant
    top = Path(tmp) / f"m{i}"
    copy_tree(top)
    apply(top, fname, edits)
    repin(top)
    code, last, failed = run_suite(top)
    shutil.rmtree(top)
    caught = code != 0
    print(f"{i:2d} {'caught' if caught else 'NOT CAUGHT':10s} {name}: {last}", file=sys.stderr, flush=True)
    return (i, name, fname, "caught" if caught else "NOT CAUGHT", last, len(failed),
            ", ".join(failed[:4]) + (" ..." if len(failed) > 4 else ""))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out")
    ap.add_argument("--jobs", type=int, default=4, help="mutants run at once")
    args = ap.parse_args(argv)
    for name, fname, edits in MUTANTS:  # every edit applies, before any suite runs
        text = (ROOT / "scripts/outside-study" / fname).read_text()
        for old, _ in edits:
            if text.count(old) != 1:
                print(f"{name}: {old!r} occurs {text.count(old)} times", file=sys.stderr)
                return 1
    with tempfile.TemporaryDirectory() as tmp:
        base = Path(tmp) / "unchanged"
        copy_tree(base)
        code, last, _ = run_suite(base)
        if code != 0:
            print(f"the unchanged copy fails: {last}", file=sys.stderr)
            return 1
        print(f"unchanged copy: {last}", file=sys.stderr)
        with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
            rows = list(pool.map(lambda im: one(tmp, *im), enumerate(MUTANTS, start=1)))
    lines = ["mutant\tchange\tfile\tresult\tsuite_last_line\tfailing_tests\tfirst_failing_tests"]
    lines += ["\t".join(str(x) for x in r) for r in rows]
    text = "\n".join(lines) + "\n"
    if args.out:
        Path(args.out).write_text(text)
    else:
        sys.stdout.write(text)
    caught = sum(r[3] == "caught" for r in rows)
    print(f"{caught} of {len(rows)} mutants caught", file=sys.stderr)
    return 0 if caught == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
