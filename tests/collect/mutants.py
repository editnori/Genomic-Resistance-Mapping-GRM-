"""Breaks scripts/collect/collect.py on purpose and records which tests catch each break.

Run from the repository root:

    python3 -I tests/collect/mutants.py [--out results.tsv] [--jobs N]

Each mutant is a copy of scripts/collect, tests/collect and docs/data-collection.md in a temporary directory with
one change made to collect.py; the suite runs there with `python3 -I -B -W error -m unittest discover -s
tests/collect -t .`. A mutant is caught when the suite fails, and the tests that failed are named. The unchanged copy
must pass first. Exit status 1 if a mutant is not caught. Not a test module: unittest does not discover it.
"""

import argparse
import concurrent.futures
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
COPIED = ("scripts/collect", "tests/collect", "docs/data-collection.md")
TOOL = "scripts/collect/collect.py"

# (name, [(old, new), ...]); each old string must occur exactly once in collect.py.
MUTANTS = [
    # Barred words.
    ("a barred word is matched in lower case only",
     [("    folded = unicodedata.normalize(\"NFKC\", path).casefold()\n",
       "    folded = unicodedata.normalize(\"NFKC\", path)\n")]),
    ("a barred word in another width or form is let through",
     [("    folded = unicodedata.normalize(\"NFKC\", path).casefold()\n", "    folded = path.casefold()\n")]),
    ("a joined word matches only as written",
     [('r"[\\W_]*".join(map(re.escape, JOINED[w]))', "re.escape(w)")]),
    ("separators are dropped from the whole path",
     [("    return next((w for w, r in BARRED_RE if r.search(folded)), None)\n",
       "    folded = re.sub(r\"[\\W_]+\", \"\", folded)\n"
       "    return next((w for w, r in BARRED_RE if r.search(folded)), None)\n")]),
    ("one barred word, ks50, is forgotten",
     [('"ks30", "ks50", "vumc"', '"ks30", "vumc"')]),
    ("a barred directory is entered and its files are read",
     [("        if word:\n            self.refused[path] = word\n            return False\n",
       "        if word:\n            self.refused[path] = word\n            return True\n")]),
    ("a barred name under a regenerable directory is passed over",
     [("            if word:\n                self.refused[child] = word\n            elif e.is_dir",
       "            if word:\n                pass\n            elif e.is_dir")]),
    ("a link that leads to a barred place is not resolved",
     [(" or barred_word(os.path.realpath(path))\n", "\n")]),
    ("a path is judged only after its .. is resolved",
     [("    return barred_word(path) or barred_word(os.path.abspath(path)) or",
       "    return barred_word(os.path.abspath(path)) or")]),
    ("a root of pack is not looked at for a barred word",
     [("            raise Stop(3, f\"pack: refused: the root {label}={top} holds the barred word {word!r}\")\n",
       "            pass\n")]),
    ("a barred word in an argument is let through",
     [("            word = barred_on_disk(path)\n            if word:\n                raise Stop(3, f\"{a.command}",
       "            word = None\n            if word:\n                raise Stop(3, f\"{a.command}")]),
    ("refused paths are not counted by top folder",
     [('                + "\\n".join(f"  {top}  {n}" for top, n in sorted(tops.items())))',
       '                + "")')]),
    ("record-refused writes a refused row that holds no barred word",
     [("        if hidden:\n", "        if False:\n")]),
    ("pack does not look for a barred word in --published",
     [("paths=lambda a: [a.list, a.out] + (a.published or []))", "paths=lambda a: [a.list, a.out])")]),
    ("summary does not look for a barred word in --out",
     [(" +\n                   ([a.out] if a.out else []))", ")")]),
    ("list does not look for a barred word in a root's label",
     [("word = barred_word(path) or barred_on_disk(real)", "word = barred_on_disk(real)")]),
    ("summary takes a place name with a barred word",
     [('            raise Stop(3, f"summary: refused: the place name {name} holds the barred word {word!r}")\n',
       "            pass\n")]),
    ("list writes its list under one of its roots",
     [("        if out == real or out.startswith(real.rstrip(os.sep) + os.sep):\n", "        if False:\n")]),
    ("a --strip-prefix may end inside a name",
     [("    if a.strip_prefix and not a.strip_prefix.endswith(\"/\"):\n", "    if False:\n")]),
    ("a table's barred row is listed like any other",
     [("            word = barred_word(given)\n", "            word = None\n")]),
    ("a list may give a barred path a packed class",
     [("        if word and cls != REFUSED:\n", "        if False:\n")]),
    ("check opens a file with a barred name",
     [("                if barred_word(child):\n                    found.refused.append(child)\n                elif",
       "                if False:\n                    found.refused.append(child)\n                elif")]),
    ("a file system encoding that is not UTF-8 is let through",
     [('    if codecs.lookup(encoding).name != "utf-8":\n', "    if False:\n")]),
    # Classes.
    ("the last matching rule gives the class, not the first",
     [("        return next((cls for cls, r in self.files if r.match(path)), None)",
       "        return next((cls for cls, r in reversed(self.files) if r.match(path)), None)")]),
    ("a build directory is walked and listed file by file",
     [("        return next((m for m in REGENERABLE_MARKERS if os.path.lexists(os.path.join(real, m))), None)",
       "        return None")]),
    ("a directory a marker made regenerable is not named",
     [("                    if why in REGENERABLE_MARKERS:\n", "                    if False:\n")]),
    ("a table's directory a marker made regenerable is not named",
     [("        walk.by_marker.append((inside + \"/\", marker, files))\n", "        pass\n")]),
    ("a directory is left out of the archive for being named target",
     [('REGENERABLE_DIRS = ("node_modules", "__pycache__")',
       'REGENERABLE_DIRS = ("node_modules", "__pycache__", "target")')]),
    ("a table's marker file does not mark its directory",
     [("                self.rules.marked[inside] = name", "                pass")]),
    ("a file no rule matches is left out without a word",
     [("    if walk.unclassed:\n", "    if False:\n")]),
    ("a table's file under node_modules takes its class from the rules",
     [("        if any(self.regenerable_dir(\"/\".join(parts[:i])) for i in range(1, len(parts))):",
       "        if False:")]),
    ("pack takes every packed class into one archive",
     [("    rows = [r for r in read_list(a.list, LIST_HEADER) if r[3] == cls]",
       "    rows = [r for r in read_list(a.list, LIST_HEADER) if r[3] in PACKABLE]")]),
    ("a public archive does not ask for the published checksum file",
     [("    if cls == PUBLIC:\n", "    if False:\n")]),
    ("a public file missing from the published checksums is packed anyway",
     [("        missing = [m[0] for m in members if (m[2], base_name(m[0])) not in published]",
       "        missing = []")]),
    ("a public file is matched by its sha256 alone",
     [("        missing = [m[0] for m in members if (m[2], base_name(m[0])) not in published]",
       "        missing = [m[0] for m in members if m[2] not in {d for d, _ in published}]")]),
    ("a published file keeps one name for each sha256",
     [("        published = {(digest, base_name(n)) for p in a.published for digest, n in read_hashes(p)}",
       "        published = set({d: base_name(n) for p in a.published for d, n in read_hashes(p)}.items())")]),
    ("regenerable files can be packed",
     [("    if cls not in PACKABLE:\n        raise Stop(2, f\"pack: {cls} is never packed\")\n", "")]),
    # Sizes and hashes.
    ("only the first block of a file is hashed",
     [("    while block := f.read(BLOCK):\n        h.update(block)\n        n += len(block)\n",
       "    if block := f.read(BLOCK):\n        h.update(block)\n        n += len(block)\n")]),
    ("a part is never closed for being full",
     [("        if used + cost > room:\n", "        if False:\n")]),
    ("a part's room forgets that data may grow when compressed",
     [("    while room > 0 and packed_bound(room, a.compress) > a.max_bytes:\n", "    while False:\n")]),
    ("a file changed since the list is packed as it now is",
     [("                        if reader.h.hexdigest() != digest or f.read(1):\n",
       "                        if False:\n")]),
    ("pack follows a link put in a file's place",
     [("                        fd = os.open(real, os.O_RDONLY | os.O_NOFOLLOW)\n",
       "                        fd = os.open(real, os.O_RDONLY)\n")]),
    ("pack writes over the parts that are there",
     [("        if os.path.lexists(t):\n            raise Stop(2, f\"pack: {t} exists; nothing was written\")\n", "")]),
    ("a stopped pack leaves its half-written parts",
     [("    except BaseException:\n        for t in written:\n", "    except BaseException:\n        for t in []:\n")]),
    # The check.
    ("check does not compare sizes",
     [("        if size != r[1]:\n            wrong_size.append", "        if False:\n            wrong_size.append")]),
    ("check does not compare sha256",
     [("        elif digest != r[2]:\n", "        elif False:\n")]),
    ("check does not look for files the list does not name",
     [("    extra = sorted(set(found.files) - set(listed))", "    extra = []")]),
    ("check does not look for listed files that are missing",
     [("    missing = [r for r in rows if r[0] not in found.files]", "    missing = []")]),
    ("check follows a link in an unpacked archive",
     [("                elif e.is_file(follow_symlinks=False):\n                    paths[child] = e.path\n",
       "                elif e.is_file():\n                    paths[child] = e.path\n")]),
    ("check lets a list name two classes",
     [("    if len(classes) != 1 or classes[0] not in PACKABLE:\n", "    if False:\n")]),
    ("check does not compare a part with its sha256",
     [("        if r[\"sha256\"] is not None and r[\"sha256\"] != present[n][1]:\n", "        if False:\n")]),
    ("check DIR does not compare the folder's name with the list's class",
     [("    elif parts is None and not os.path.basename(os.path.normpath(where)).endswith(\".\" + classes[0]):\n",
       "    elif False:\n")]),
    ("check --parts passes a part file beside it that the parts file does not name",
     [("    if beside:\n", "    if False:\n")]),
    ("check passes a part that is not there",
     [("        found.faults.append(f\"the part {os.path.basename(parts[n][0])} is not beside {base}\")\n",
       "        pass\n")]),
    ("check does not mind a file in another part than listed",
     [("        elif part is not None and part != r[4]:\n", "        elif False:\n")]),
    ("check ignores what lies beside the list and the files",
     [("    if found.stray:\n", "    if False:\n")]),
    ("check takes the last of two LIST.tsv members of a part",
     [('and m.isreg() and out["list"] is None:', "and m.isreg():")]),
    ("check DIR follows a link put in place of LIST.tsv",
     [("if not os.path.isfile(inner) or os.path.islink(inner):", "if not os.path.isfile(inner):")]),
    ("check DIR follows a link put in place of files/",
     [("if os.path.isdir(files) and not os.path.islink(files):", "if os.path.isdir(files):")]),
    ("check DIR passes a directory with no file under files/",
     [("    if found.empty:\n", "    if False:\n")]),
    ("a name that is not UTF-8 found twice is printed as it is",
     [("found.faults.append(f\"{repr(path) if unsafe(path) else path} is in the archive twice\")",
       "found.faults.append(f\"{path} is in the archive twice\")")]),
    ("check reads only the first part's list",
     [("    if len(set(found.lists.values())) > 1:\n", "    if False:\n")]),
    ("check stops reading a part at its tar end mark",
     [("                stream.end(tar.offset)\n", "")]),
    ("check lets bytes follow a part's compressed stream",
     [("        if unz.unused_data or getattr(unz, \"unconsumed_tail\", b\"\") or self.raw.read(1):\n",
       "        if False:\n")]),
    ("check lets bytes that are not zero follow a part's last member",
     [("            if data.count(0) != len(data):\n", "            if False:\n")]),
    ("check lets a part lack its tar end mark",
     [("        if zeros < 2 * tarfile.BLOCKSIZE:\n", "        if False:\n")]),
    # The summary.
    ("copies are counted per path, not per place",
     [("                holders.setdefault(digest, set()).add(place)",
       "                holders.setdefault(digest, set()).add(path)")]),
    ("duplicates are not removed from the byte count",
     [("        return {r[3]: r[2] for r in of if r[3] != \"-\"}",
       "        return {r[:2]: r[2] for r in of if r[3] != \"-\"}")]),
    ("a place known by its hashes is not counted as a copy",
     [("        for d in named & set(holders):\n            holders[d].add(place)\n", "")]),
]


def suite(tree):
    p = subprocess.run([sys.executable, "-I", "-B", "-W", "error", "-m", "unittest", "discover", "-s", "tests/collect",
                        "-t", "."], cwd=tree, capture_output=True, text=True, stdin=subprocess.DEVNULL)
    failed = sorted(set(re.findall(r"^(?:FAIL|ERROR): (\S+) \(tests\.collect\.(\w+)\.", p.stderr, flags=re.M)))
    return p.returncode, [f"{module}.{test}" for test, module in failed], p.stderr.strip().splitlines()[-1:]


def copy(tmp):
    for rel in COPIED:
        src, dst = ROOT / rel, Path(tmp) / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        if src.is_dir():
            shutil.copytree(src, dst, ignore=shutil.ignore_patterns("__pycache__"))
        else:
            shutil.copy2(src, dst)


def one(mutant):
    name, changes = mutant
    with tempfile.TemporaryDirectory(prefix="collect-mutant-") as tmp:
        copy(tmp)
        tool = Path(tmp) / TOOL
        text = tool.read_text(encoding="utf-8")
        for old, new in changes:
            if text.count(old) != 1:
                return name, None, [f"the text to change occurs {text.count(old)} times: {old!r}"]
            text = text.replace(old, new)
        tool.write_text(text, encoding="utf-8")
        status, failed, last = suite(tmp)
        return name, status != 0, failed or last


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    with tempfile.TemporaryDirectory(prefix="collect-mutant-") as tmp:
        copy(tmp)
        status, failed, last = suite(tmp)
        if status != 0:
            sys.exit(f"the unchanged copy does not pass: {failed or last}")
    print(f"the unchanged copy passes: {last[0] if last else ''}")
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, a.jobs)) as pool:
        results = list(pool.map(one, MUTANTS))
    lines = ["mutant\tcaught\tfailing_tests\ttests"]
    for name, caught, tests in results:
        state = "not applied" if caught is None else "yes" if caught else "NO"
        lines.append(f"{name}\t{state}\t{len(tests) if caught else 0}\t{' '.join(tests)}")
    text = "\n".join(lines) + "\n"
    if a.out:
        Path(a.out).write_text(text, encoding="utf-8")
    for name, caught, tests in results:
        state = "not applied" if caught is None else "caught" if caught else "NOT CAUGHT"
        print(f"{state:11} {name}: {len(tests)} tests" + (f", first {tests[0]}" if tests else ""))
    missed = [r for r in results if not r[1]]
    print(f"{len(results) - len(missed)} of {len(results)} mutants caught")
    sys.exit(1 if missed else 0)


if __name__ == "__main__":
    main()
