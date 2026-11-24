"""Private details that no tracked file may hold, and the scan that finds them.

Six kinds:

- K1, a home path under a name that is not a neutral one;
- K2, a Windows user path;
- K3, a link to an agent thread;
- K4, a box id;
- K5, a host name of a rented or private machine, or a public address;
- K6, the name of the project's owner or of one of their accounts, outside a licence or an author line.

The names of K6 are held as sha256 values of lowercased words, so this file does not print them. More names can
come from the environment variable PUBLIC_WORDS_EXTRA, one word a line, which CI fills from a repository secret of
the same name. They are hashed as they are read and no word of them is printed. Without the variable the scan says
"extra words: not checked" and checks the names held here. Each pattern is put together from pieces, so that no
pattern matches its own source and this suite passes its own scan.

    python3 -I tests/public_words/rules.py [ROOT]

prints the count of hits and of files for each kind, whether extra words were checked, then each hit as kind,
path and line number. It never prints the text of a hit, and a barred name in a path reads ***. It exits 1 when
there is a hit, or when a line of the extra words can match nothing.
"""

import collections
import hashlib
import os
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
SUITE = "tests/public_words"
WORKFLOW = ROOT / ".github" / "workflows" / "ci.yml"
CI_DOC = ROOT / "docs" / "ci.md"
REMOVED_PAGE = ("docs", "work-program-" + "20261027.md")
KINDS = ("K1", "K2", "K3", "K4", "K5", "K6")
EXTRA = "PUBLIC_WORDS_EXTRA"
MASK = "***"

S = "/"
NAME = r"([A-Za-z0-9_.-]+)"
HOME = re.compile(S + "home" + S + NAME + "|" + S + "Users" + S + NAME)
NEUTRAL_HOMES = frozenset({"user", "runner", "proj"})
WINDOWS = re.compile(r"[A-Za-z]:\\+" + "Users" + r"\\+" + "|" + S + "mnt" + S + "[a-z]" + S + "Users" + S)
WSL_MOUNT = re.compile(S + "mnt" + S + "[a-z]$")
THREAD = re.compile("t3" + "-" + "[a-z]+" + ":" + S + S)
BOX = re.compile(r"\bbx" + "_" + r"[a-z0-9]{8}\b")
HOST = re.compile(r"\bssh[0-9]*\." + "vast" + r"\.ai\b|\b[a-z0-9-]+\." + "ts" + r"\.net\b|\b(?:DESKTOP|LAPTOP)-[A-Z0-9]{7}\b")
ADDRESS = re.compile(r"(?<![0-9.])([0-9]{1,3})\.([0-9]{1,3})\.([0-9]{1,3})\.([0-9]{1,3})(?![0-9.])")
LOCK_FILES = ("Cargo.lock", "package-lock.json")   # their version numbers can have four parts

# sha256 of each barred name, lowercased.
NAMES = frozenset({
    "d5c8bfa6741c39d916643cdea20893bb1833c7780e13746e1410b67c44ceaaca",
    "4493ad3410c566c8ebb5df5400347d8500532bf1c639134ecf928879e4629850",
})
WORD = re.compile(r"[A-Z]?[a-z]+|[A-Z]+(?![a-z])")
RUN = re.compile(r"[A-Za-z0-9]+")
LICENCE_FILE = re.compile(r"(^|/)(LICENSE|LICENCE|COPYING|NOTICE)(\.[A-Za-z]+)?$")
AUTHOR_LINE = re.compile(r"^\s*(authors?\s*[=:]|\"authors?\"\s*:|Author:|Copyright\b)", re.IGNORECASE)


def sha(word):
    return hashlib.sha256(word.lower().encode("utf-8")).hexdigest()


def lineage_name(match, text):
    """True for a lineage name as TBProfiler writes it (lineage4.1.2.1): the 7 characters before are `lineage`, no
    letter stands right before those 7, the first number is 1 to 9 and each of the other three is under 100."""
    start = match.start()
    if text[max(0, start - 7):start] != "lineage" or (start >= 8 and text[start - 8].isalpha()):
        return False
    first, *rest = match.groups()
    return len(first) == 1 and first != "0" and all(int(g) < 100 for g in rest)


def public_address(match, text):
    """True for an address outside the loopback, private and documentation ranges. An enzyme number is not one,
    nor a lineage name as TBProfiler writes it."""
    parts = tuple(int(g) for g in match.groups())
    if max(parts) > 255 or text[max(0, match.start() - 3):match.start()] == "EC " or lineage_name(match, text):
        return False
    a, b, c, _ = parts
    if a in (0, 10, 127, 255) or (a, b) in ((192, 168), (169, 254)) or (a == 172 and 16 <= b <= 31):
        return False
    return (a, b, c) not in ((192, 0, 2), (198, 51, 100), (203, 0, 113))


def barred_words(text, names=NAMES, seen=None):
    """The words of the text whose sha256 is a barred name: each letter word, split at case changes, and each run
    of letters and digits. The dict `seen` keeps the answer for a word across calls."""
    seen = {} if seen is None else seen
    found = set()
    for word in set(WORD.findall(text)) | set(RUN.findall(text)):
        key = word.lower()
        if key not in seen:
            seen[key] = sha(key) in names
        if seen[key]:
            found.add(word)
    return found


def extra_names(value=None):
    """(sha256 values, lines that can match nothing) for the extra barred names: one word a line in the environment
    variable, or in `value` when given. (None, 0) when there is no word. A word is hashed as it is read. A line
    that is not one run of letters and digits can match no word of a file, so it is counted and not used."""
    value = os.environ.get(EXTRA, "") if value is None else value
    words = [line.strip() for line in value.splitlines() if line.strip()]
    if not words:
        return None, 0
    usable = [word for word in words if RUN.fullmatch(word)]
    return frozenset(sha(word) for word in usable), len(words) - len(usable)


def extra_state(value=None):
    """The line that says whether extra words were checked. It gives counts and no word."""
    names, unusable = extra_names(value)
    if names is None:
        return "extra words: not checked"
    return f"extra words: {len(names)} checked" + (f", not usable: {unusable}" if unusable else "")


def checked_names(value=None):
    """The names held here and the extra ones."""
    return NAMES | (extra_names(value)[0] or frozenset())


def shown(path, names=NAMES, seen=None):
    """The path as it may be printed: a barred name in it reads ***."""
    words = barred_words(path, names, seen)
    if not words:
        return path
    return re.sub("|".join(sorted((re.escape(w) for w in words), key=len, reverse=True)), MASK, path)


def listing(hits, names=NAMES):
    """One line for each hit: kind, path and line number. Never the text, and never a barred name in a path."""
    seen = {}
    return [f"{kind}\t{shown(path, names, seen)}:{number}" for kind, path, number in hits]


def line_of(text, offset):
    return text.count("\n", 0, offset) + 1


def hits_in(path, text, names=NAMES, seen=None):
    """(kind, line number) for each hit in one file's text, given its path in the tree with forward slashes."""
    out = []
    for m in HOME.finditer(text):
        name = m.group(1) or m.group(2)
        if WSL_MOUNT.search(text[max(0, m.start() - 6):m.start()]):
            continue   # the Windows view of a user path: K2 counts it
        if re.search("[A-Za-z]", name) and name not in NEUTRAL_HOMES:
            out.append(("K1", line_of(text, m.start())))
    out += [("K2", line_of(text, m.start())) for m in WINDOWS.finditer(text)]
    out += [("K3", line_of(text, m.start())) for m in THREAD.finditer(text)]
    out += [("K4", line_of(text, m.start())) for m in BOX.finditer(text)]
    out += [("K5", line_of(text, m.start())) for m in HOST.finditer(text)]
    if not path.endswith(LOCK_FILES):
        out += [("K5", line_of(text, m.start())) for m in ADDRESS.finditer(text) if public_address(m, text)]
    words = barred_words(text, names, seen)
    if words and not LICENCE_FILE.search(path):
        pattern = re.compile("|".join(sorted(re.escape(w) for w in words)))
        for number, line in enumerate(text.split("\n"), 1):
            if not AUTHOR_LINE.search(line):
                out += [("K6", number)] * len(pattern.findall(line))
    return sorted(out)


def tracked(root=ROOT):
    """The tracked files of the tree as git lists them. Without git: every file but build and package folders."""
    root = pathlib.Path(root)
    try:
        done = subprocess.run(["git", "-C", str(root), "ls-files", "-z"], capture_output=True, check=True)
        names = [n for n in done.stdout.decode("utf-8").split("\0") if n]
    except (OSError, subprocess.CalledProcessError):
        names = []
    if not names:
        for p in sorted(root.rglob("*")):
            parts = p.relative_to(root).parts
            if p.is_file() and parts[0] not in (".git", "target", "data") \
                    and not {"node_modules", "__pycache__", ".venv"} & set(parts):
                names.append("/".join(parts))
    return names


def findings(root=ROOT, names=None):
    """(kind, path, line number) for each hit in the tree. A path that is itself a hit has line 0. Files that are
    not UTF-8 text are not read; the second value lists them. Without `names`, the names held here and the extra
    ones are checked."""
    root = pathlib.Path(root)
    names = checked_names() if names is None else names
    seen, hits, unread = {}, [], []
    for path in tracked(root):
        hits += [(kind, path, 0) for kind, _ in hits_in("", path, names, seen)]
        file = root / path
        if file.is_symlink() or not file.is_file():
            continue
        try:
            text = file.read_bytes().decode("utf-8")
        except UnicodeDecodeError:
            unread.append(path)
            continue
        hits += [(kind, path, number) for kind, number in hits_in(path, text, names, seen)]
    return hits, unread


def counts(hits):
    """kind -> (hits, files)."""
    n = collections.Counter(kind for kind, _, _ in hits)
    files = collections.defaultdict(set)
    for kind, path, _ in hits:
        files[kind].add(path)
    return {kind: (n[kind], len(files[kind])) for kind in KINDS}


def main(argv):
    root = pathlib.Path(argv[1]) if len(argv) > 1 else ROOT
    names = checked_names()
    hits, unread = findings(root, names)
    for kind, (n, files) in counts(hits).items():
        print(f"{kind}\t{n} hits\t{files} files")
    print(f"files read\t{len(tracked(root)) - len(unread)}\tnot text\t{len(unread)}")
    print(extra_state())
    for line in listing(hits, names):
        print(line)
    return 1 if hits or extra_names()[1] else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
