"""collect.py list|pack|check|summary

One checked collection of a project's files: a list with a size, a sha256 and a class for every file, archives that
hold one class each, and a check of an archive against the list it carries. docs/data-collection.md shows a whole
run. Standard library only; run it as `python3 -I scripts/collect/collect.py ...`.

  list     --rules RULES --out LIST [--root LABEL=PATH]... [--listing FILE [--strip-prefix P]] [--record-refused]
           [--skip-special] [--threads N]
           Walks each PATH (a directory, or one file) and writes LIST, one row per file: path, bytes, sha256, class.
           A path starts with its root's LABEL, so a list reads the same on any machine. --listing takes the rows of
           an existing table (a header naming path, bytes and sha256) in place of a walk: only the table is opened.
  pack     --list LIST --class CLASS --name NAME --out DIR [--root LABEL=PATH]... [--published SHA256FILE]...
           [--max-bytes N] [--compress gz|xz|none]
           Writes the files of one class into NAME.CLASS.partNNN.tar[.gz|.xz], each at most N bytes (1,900,000,000 by
           default), and NAME.CLASS.parts.sha256 beside them. Every part holds the class's whole list as
           NAME.CLASS/LIST.tsv, with the part of each file, and its own files under NAME.CLASS/files/.
  check    DIR | --parts NAME.CLASS.parts.sha256   [--threads N]
           DIR is an unpacked NAME.CLASS folder. Every row of LIST.tsv must be a regular file of that size and
           sha256 under files/, and files/ must hold nothing else, and DIR's name must end in the list's class.
           --parts checks the packed parts where they lie, without unpacking: each part's sha256, then every file
           inside it, read as a stream.
           Each part is read to the end of its gzip or xz stream, so the stream's own check (gzip's CRC and length,
           xz's check) is read too; after the last member only the tar end mark, zero bytes, may follow, and nothing
           may follow the compressed stream.
  summary  --place NAME=LIST... [--hashes NAME=SHA256FILE]... [--out FILE]
           Reads the lists of several places and counts, by sha256, how many places hold each content. A --hashes
           file (sha256sum format) is a place known by its hashes only, such as the checksum file of a release.
           Prints four tables and writes every row with its count to FILE.

Classes. `public-result`, `internal-record` and `public-input` are packed. `regenerable` is listed and never packed:
builds, caches, environments, and an unpacked copy of an archive that is itself listed. Such a directory is one row
with a trailing slash and the total size of its files, and no file in it is opened. It is known by its name
(REGENERABLE_DIRS: node_modules, __pycache__), by a file it holds (REGENERABLE_MARKERS: CACHEDIR.TAG, which cargo,
pytest and other tools write into their build and cache directories, or pyvenv.cfg, a Python environment), or by a
rules pattern that ends in "/". A name such as target or .cargo alone does not make a directory regenerable.
`refused` is any path holding one of the words in BARRED. The path is read after NFKC and casefold, so letter case
and letter width do not matter, and the four joined words held-out, ks30, ks50 and sigdrug match with any run of
characters that are not letters or digits between their halves, or none: held_out, KS-30, sig drug and ks/50 are
refused. The other five match as written; nothing is dropped from the whole path, so tasks/data is not ksd. `list`
names every such path, counts them by top folder, and stops with status 3 before it writes anything; it never opens
one and never enters one. A refused path is in no archive, so look at each by hand. With --record-refused it writes
them as rows with no size and no hash, and goes on, unless a refused path holds no barred word in its list path (a
link or a root that leads to a barred place, or a word that --strip-prefix cut off): a refused row must hold its
word, so that it cannot be told from a row edited by hand to hide a file. `pack` and `check` stop on such
a path too, and so does every command given one as an argument.

RULES gives the other classes: one line per rule, CLASS<TAB>PATTERN, the first match wins, "#" starts a comment.
A pattern is matched against the whole list path with fnmatch, so "*" also crosses "/". A file that no rule matches
stops the walk with its name: end the file with "internal-record<TAB>*" to make that the default.

Two classes never share an archive: `pack` takes one class, and `check` fails a list that names two. A public-result
archive needs a second key. Each of its files must also stand in a --published file (sha256sum format: the checksum
file of a release already published) with its sha256 and its file name together, so a file is not made public by a
pattern alone. The name compared is the base name, the part after the last "/", on both sides; every name a published
file gives for a sha256 counts, and an empty file is matched like any other.

A list holds regular files. A link or a special file stops `list` with its name (--skip-special leaves it out), and
an empty directory is not recorded.

Exit status: 0 done; 1 `check` found a difference; 2 an input the command cannot use, named on stderr; 3 a path with
a barred word.
"""
import argparse
import concurrent.futures
import fnmatch
import gzip
import hashlib
import io
import lzma
import os
import re
import stat
import sys
import tarfile
import unicodedata
import zlib

BARRED = ("sealed", "held-out", "ks30", "ks50", "vumc", "clinical", "private", "ksd", "sigdrug")
# The two halves of each joined word; any run of characters that are not letters or digits may stand between them.
JOINED = {"held-out": ("held", "out"), "ks30": ("ks", "30"), "ks50": ("ks", "50"), "sigdrug": ("sig", "drug")}
BARRED_RE = tuple((w, re.compile(r"[\W_]*".join(map(re.escape, JOINED[w])) if w in JOINED else re.escape(w)))
                  for w in BARRED)
PACKABLE = ("public-result", "internal-record", "public-input")
REGENERABLE, REFUSED = "regenerable", "refused"
CLASSES = PACKABLE + (REGENERABLE, REFUSED)
PUBLIC = "public-result"
REGENERABLE_DIRS = ("node_modules", "__pycache__")
REGENERABLE_MARKERS = ("CACHEDIR.TAG", "pyvenv.cfg")
LIST_HEADER = ("path", "bytes", "sha256", "class")
INNER_HEADER = LIST_HEADER + ("part",)
INNER_LIST = "LIST.tsv"
PARTS_SUFFIX = ".parts.sha256"
PART_BYTES = 1_900_000_000
BLOCK = 1 << 20
TAR_END = 2 * tarfile.BLOCKSIZE + tarfile.RECORDSIZE  # the two closing blocks, and the padding to a whole record
MTIME_MAX = 8 ** 11 - 1  # the largest time a plain tar header holds
EXTENSION = {"gz": ".tar.gz", "xz": ".tar.xz", "none": ".tar"}
STREAM = {"gz": "gzip", "xz": "xz", "none": "tar"}
SHOWN = 20  # names printed per kind of problem
HEX64 = re.compile(r"[0-9a-f]{64}\Z")


class Stop(Exception):
    def __init__(self, status, message):
        super().__init__(message)
        self.status = status


class Damaged(Exception):
    """A part whose stream is not whole."""


def barred_word(path):
    """The first barred word in the path, read after NFKC and casefold, or None."""
    folded = unicodedata.normalize("NFKC", path).casefold()
    return next((w for w, r in BARRED_RE if r.search(folded)), None)


def barred_on_disk(path):
    """barred_word over the path as given, then made absolute, then with its links resolved. The text as given comes
    first, so a path such as sealed/../x is refused before anything on disk is looked at."""
    return barred_word(path) or barred_word(os.path.abspath(path)) or barred_word(os.path.realpath(path))


def unsafe(path):
    """Why a list cannot hold this path, or None. A directory row may end in one slash."""
    try:
        path.encode("utf-8")
    except UnicodeEncodeError:
        return "is not UTF-8"
    if any(ord(c) < 32 or ord(c) == 127 for c in path):
        return "holds a control character"
    parts = (path[:-1] if path.endswith("/") else path).split("/")
    if path.startswith("/") or any(p in ("", ".", "..") for p in parts):
        return "is not a plain relative path"
    return None


def shown(names):
    names = sorted(names)
    more = [f"  ... and {len(names) - SHOWN} more"] if len(names) > SHOWN else []
    return "\n".join([f"  {n!r}" if unsafe(n) else f"  {n}" for n in names[:SHOWN]] + more)


def sha256_stream(f):
    h, n = hashlib.sha256(), 0
    while block := f.read(BLOCK):
        h.update(block)
        n += len(block)
    return h.hexdigest(), n


def sha256_file(path):
    """(sha256, bytes read) of a regular file, opened without following a link."""
    fd = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    with os.fdopen(fd, "rb") as f:
        if not stat.S_ISREG(os.fstat(f.fileno()).st_mode):
            raise OSError(f"{path} is not a regular file")
        return sha256_stream(f)


def in_threads(work, items, threads):
    """work over items, in order. hashlib and zlib release the lock on large blocks, so threads help."""
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, threads)) as pool:
        return list(pool.map(work, items))


# ---- rules ----

class Rules:
    def __init__(self, path):
        self.files, self.dirs, self.marked = [], [], set()
        try:
            with open(path, encoding="utf-8") as f:
                lines = f.read().split("\n")
        except (OSError, UnicodeDecodeError) as e:
            raise Stop(2, f"rules: cannot read {path}: {e}") from None
        for n, line in enumerate(lines, 1):
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.split("\t")
            if len(cells) != 2 or not cells[1]:
                raise Stop(2, f"rules: {path}:{n}: want CLASS<TAB>PATTERN, read {line!r}")
            cls, pattern = cells
            if cls not in PACKABLE + (REGENERABLE,):
                raise Stop(2, f"rules: {path}:{n}: {cls!r} is not a class a rule can give; "
                              f"use one of {', '.join(PACKABLE + (REGENERABLE,))}")
            compiled = re.compile(fnmatch.translate(pattern))
            if pattern.endswith("/"):
                if cls != REGENERABLE:
                    raise Stop(2, f"rules: {path}:{n}: a pattern ending in '/' names a directory, "
                                  f"and only {REGENERABLE} takes one")
                self.dirs.append(compiled)
            else:
                self.files.append((cls, compiled))

    def regenerable_dir(self, path, real=None):
        """path is a directory's list path, without the trailing slash; real is where it lies, if it is on disk."""
        if path.rsplit("/", 1)[-1] in REGENERABLE_DIRS or any(r.match(path + "/") for r in self.dirs):
            return True
        if real is None:
            return path in self.marked
        return any(os.path.lexists(os.path.join(real, m)) for m in REGENERABLE_MARKERS)

    def file_class(self, path):
        return next((cls for cls, r in self.files if r.match(path)), None)

    def listed_class(self, path):
        """The class of a path read from a table: a file under a regenerable directory is regenerable."""
        parts = path.split("/")
        if any(self.regenerable_dir("/".join(parts[:i])) for i in range(1, len(parts))):
            return REGENERABLE
        return self.file_class(path)


# ---- lists ----

def parse_list(text, header, where):
    """The rows of a list as (path, bytes, sha256, class[, part]); bytes is an int, or None on a refused row."""
    lines = text.split("\n")
    if lines[0].split("\t") != list(header) or lines[-1] != "":
        raise Stop(2, f"{where}: not a list: the first line must be {' '.join(header)} and the last line must end")
    rows, seen = [], set()
    for n, line in enumerate(lines[1:-1], 2):
        cells = line.split("\t")
        if len(cells) != len(header):
            raise Stop(2, f"{where}:{n}: {len(cells)} cells, want {len(header)}")
        name, size, digest, cls = cells[:4]
        bad = unsafe(name)
        if bad:
            raise Stop(2, f"{where}:{n}: the path {name!r} {bad}")
        if cls not in CLASSES:
            raise Stop(2, f"{where}:{n}: unknown class {cls!r}")
        word = barred_word(name)
        if word and cls != REFUSED:
            raise Stop(3, f"{where}:{n}: refused: {name} holds the barred word {word!r} and is listed as {cls}")
        if cls == REFUSED and not word:
            raise Stop(2, f"{where}:{n}: {name} is listed as {REFUSED} and holds no barred word")
        if name in seen:
            raise Stop(2, f"{where}:{n}: the path {name} is listed twice")
        seen.add(name)
        if not (size.isascii() and size.isdigit() or size == "-" and cls == REFUSED):
            raise Stop(2, f"{where}:{n}: bytes is {size!r}")
        if not (HEX64.match(digest) or digest == "-" and cls in (REGENERABLE, REFUSED)):
            raise Stop(2, f"{where}:{n}: sha256 is {digest!r}")
        if name.endswith("/") and cls != REGENERABLE:
            raise Stop(2, f"{where}:{n}: {name} is a directory row and its class is {cls}")
        if len(cells) == 5 and not re.fullmatch(r"[0-9]{3}", cells[4]):
            raise Stop(2, f"{where}:{n}: part is {cells[4]!r}")
        rows.append((name, None if size == "-" else int(size), digest, cls) + tuple(cells[4:]))
    return rows


def read_text(path, what):
    try:
        with open(path, encoding="utf-8", newline="") as f:
            return f.read()
    except (OSError, UnicodeDecodeError) as e:
        raise Stop(2, f"cannot read {what} {path}: {e}") from None


def read_list(path, header):
    return parse_list(read_text(path, "the list"), header, path)


def list_text(rows, header):
    out = ["\t".join(header)]
    out += ["\t".join("-" if c is None else str(c) for c in r) for r in sorted(rows)]
    return "\n".join(out) + "\n"


def write_new(path, data):
    """Writes data to path through a temporary name, so a stopped run leaves no half file."""
    tmp = f"{path}.tmp{os.getpid()}"
    try:
        with open(tmp, "wb") as f:
            f.write(data)
        os.replace(tmp, path)
    except OSError as e:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise Stop(2, f"cannot write {path}: {e}") from None


def parse_named(values, what, nested=False):
    """[NAME=VALUE] as a dict in the order given. A name is one plain path component, or a plain path if nested."""
    out = {}
    for v in values or []:
        name, eq, value = v.partition("=")
        if not eq or not value or unsafe(name) or name.endswith("/") or "/" in name and not nested:
            raise Stop(2, f"{what}: want NAME=VALUE with a plain name, read {v!r}")
        if name in out:
            raise Stop(2, f"{what}: {name} is given twice")
        out[name] = value
    return out


def read_hashes(path):
    """[(sha256, name)] of a sha256sum-format file, every line in its order."""
    out = []
    for n, line in enumerate([ln for ln in read_text(path, "the sha256 file").split("\n") if ln], 1):
        if not (HEX64.match(line[:64]) and line[64:66] in ("  ", " *") and len(line) > 66):
            raise Stop(2, f"{path}:{n}: not a sha256sum line: {line[:80]!r}")
        out.append((line[:64], line[66:]))
    return out


def base_name(path):
    return path.rsplit("/", 1)[-1]


# ---- list ----

class Walk:
    def __init__(self, rules):
        self.rules = rules
        self.rows, self.to_hash = [], []
        self.refused, self.special, self.unclassed, self.unnamed = {}, [], [], []

    def meet(self, path, real):
        """True if the entry may be looked at. A barred word in its list path or its path on disk refuses it."""
        if unsafe(path):
            self.unnamed.append(path)
            return False
        word = barred_word(path) or barred_on_disk(real)
        if word:
            self.refused[path] = word
            return False
        return True

    def root(self, label, top):
        if not self.meet(label, top):
            return
        try:
            st = os.lstat(top)
        except OSError as e:
            raise Stop(2, f"list: cannot read the root {label}={top}: {e}") from None
        if stat.S_ISDIR(st.st_mode):
            self.directory(top, label)
        elif stat.S_ISREG(st.st_mode):
            self.file(top, label, st.st_size)
        else:
            self.special.append(label)

    def file(self, real, path, size):
        cls = self.rules.file_class(path)
        if cls is None:
            self.unclassed.append(path)
        elif cls == REGENERABLE:
            self.rows.append((path, size, "-", cls))
        else:
            self.to_hash.append((path, real, cls))

    def directory(self, real, path):
        try:
            with os.scandir(real) as it:
                entries = sorted(it, key=lambda e: e.name)
        except OSError as e:
            raise Stop(2, f"list: cannot read the directory {real}: {e}") from None
        for e in entries:
            child = f"{path}/{e.name}"
            if not self.meet(child, e.path):
                continue
            if e.is_dir(follow_symlinks=False):
                if self.rules.regenerable_dir(child, e.path):
                    self.rows.append((child + "/", self.size_below(e.path, child), "-", REGENERABLE))
                else:
                    self.directory(e.path, child)
            elif e.is_file(follow_symlinks=False):
                self.file(e.path, child, e.stat(follow_symlinks=False).st_size)
            else:
                self.special.append(child)

    def size_below(self, real, path):
        """The bytes of the regular files under a regenerable directory. Names are read; no file is opened."""
        total = 0
        try:
            with os.scandir(real) as it:
                entries = list(it)
        except OSError as e:
            raise Stop(2, f"list: cannot read the directory {real}: {e}") from None
        for e in entries:
            child = f"{path}/{e.name}"
            word = barred_word(child) or barred_word(e.path)
            if word:
                self.refused[child] = word
            elif e.is_dir(follow_symlinks=False):
                total += self.size_below(e.path, child)
            elif e.is_file(follow_symlinks=False):
                total += e.stat(follow_symlinks=False).st_size
        return total

    def table(self, file, strip):
        lines = read_text(file, "the table").split("\n")
        header = lines[0].split("\t")
        if lines[-1] != "" or not {"path", "bytes", "sha256"} <= set(header):
            raise Stop(2, f"list: {file}: the first line must name path, bytes and sha256, and the last line must end")
        col = {k: header.index(k) for k in ("path", "bytes", "sha256")}
        # A table has no directories to look into: a marker file among its rows marks its directory.
        for line in lines[1:-1]:
            given = (line.split("\t") + [""] * len(header))[col["path"]]
            inside, _, name = given[len(strip):].rpartition("/")
            if name in REGENERABLE_MARKERS and inside and given.startswith(strip):
                self.rules.marked.add(inside)
        for n, line in enumerate(lines[1:-1], 2):
            cells = line.split("\t")
            if len(cells) != len(header):
                raise Stop(2, f"list: {file}:{n}: {len(cells)} cells, want {len(header)}")
            given, size, digest = (cells[col[k]] for k in ("path", "bytes", "sha256"))
            if not given.startswith(strip) or not (size.isascii() and size.isdigit()) or not HEX64.match(digest):
                raise Stop(2, f"list: {file}:{n}: want a path under {strip!r}, a size and a sha256, "
                              f"read {given!r} {size!r} {digest!r}")
            path = given[len(strip):]
            if unsafe(path) or path.endswith("/"):
                self.unnamed.append(path)
                continue
            word = barred_word(given)
            if word:
                self.refused[path] = word
                continue
            cls = self.rules.listed_class(path)
            if cls is None:
                self.unclassed.append(path)
            else:
                self.rows.append((path, int(size), digest, cls))


def cmd_list(a):
    rules = Rules(a.rules)
    roots = parse_named(a.root, "--root", nested=True)
    if not roots and not a.listing:
        raise Stop(2, "list: give at least one --root LABEL=PATH or one --listing FILE")
    walk = Walk(rules)
    for label, top in roots.items():
        walk.root(label, top)
    for file in a.listing or []:
        walk.table(file, a.strip_prefix)
    if walk.unnamed:
        raise Stop(2, f"list: {len(walk.unnamed)} names a list cannot hold (a control character, or not UTF-8, "
                      f"or not a plain relative path):\n{shown(walk.unnamed)}")
    if walk.refused:
        lines = "\n".join(f"  {p}  ({w})" for p, w in sorted(walk.refused.items())[:SHOWN])
        more = f"\n  ... and {len(walk.refused) - SHOWN} more" if len(walk.refused) > SHOWN else ""
        tops = {}
        for p in walk.refused:
            tops[p.split("/", 1)[0]] = tops.get(p.split("/", 1)[0], 0) + 1
        text = (f"{len(walk.refused)} paths hold a barred word; none was opened:\n{lines}{more}\n"
                f"list: refused paths by top folder; a refused path is in no archive, so look at it by hand:\n"
                + "\n".join(f"  {top}  {n}" for top, n in sorted(tops.items())))
        if not a.record_refused:
            raise Stop(3, f"list: refused: {text}\nlist: nothing was written")
        # A refused row must hold its word, so that it cannot be told from a row edited by hand to hide a file.
        hidden = sorted(p for p in walk.refused if not barred_word(p))
        if hidden:
            raise Stop(3, f"list: refused: {text}\nlist: {len(hidden)} of them hold no barred word in their list "
                          f"path (their place on disk does, through a link or the root, or --strip-prefix cut it "
                          f"off), and a refused row must hold its word; nothing was written:\n{shown(hidden)}")
        print(f"list: {text}", file=sys.stderr)
        walk.rows += [(p, None, "-", REFUSED) for p in walk.refused]
    if walk.special:
        text = f"{len(walk.special)} links or special files, which a list does not hold:\n{shown(walk.special)}"
        if not a.skip_special:
            raise Stop(2, f"list: {text}\nlist: nothing was written; --skip-special leaves them out")
        print(f"list: left out {text}", file=sys.stderr)
    if walk.unclassed:
        raise Stop(2, f"list: no rule of {a.rules} gives a class to {len(walk.unclassed)} files:\n"
                      f"{shown(walk.unclassed)}\nlist: nothing was written")
    try:
        hashed = in_threads(sha256_file, [real for _, real, _ in walk.to_hash], a.threads)
    except OSError as e:
        raise Stop(2, f"list: cannot read a file: {e}") from None
    walk.rows += [(path, n, digest, cls) for (path, _, cls), (digest, n) in zip(walk.to_hash, hashed)]
    names = set()
    twice = {r[0] for r in walk.rows if r[0] in names or names.add(r[0])}
    if twice:
        raise Stop(2, f"list: {len(twice)} paths are listed by two sources:\n{shown(twice)}")
    write_new(a.out, list_text(walk.rows, LIST_HEADER).encode("utf-8"))
    print("class\tentries\tbytes")
    for cls in CLASSES:
        of = [r for r in walk.rows if r[3] == cls]
        if of:
            print(f"{cls}\t{len(of)}\t{sum(r[1] or 0 for r in of)}")
    print(f"list: wrote {len(walk.rows)} rows to {a.out}")


# ---- pack ----

def tar_info(name, size, mtime=0, mode=0o644):
    info = tarfile.TarInfo(name)
    info.size, info.mode = size, mode
    info.mtime = min(max(int(mtime), 0), MTIME_MAX)  # kept inside the plain header, so its size is known ahead
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    return info


def member_cost(name, size):
    """The tar bytes of one member: its header (with a pax header for a long name or size) and its padded data."""
    header = len(tar_info(name, size).tobuf(tarfile.PAX_FORMAT, "utf-8", "surrogateescape"))
    return header + (size + tarfile.BLOCKSIZE - 1) // tarfile.BLOCKSIZE * tarfile.BLOCKSIZE


def packed_bound(n, compress):
    """The most bytes n tar bytes can come to once compressed: data that does not compress grows a little."""
    if compress == "none":
        return n
    return n + (n >> 12) + (n >> 14) + (n >> 25) + (64 if compress == "gz" else 4096)


class Hashing:
    """Reads a file for tarfile and keeps the sha256 of what was read."""

    def __init__(self, f):
        self.f, self.h = f, hashlib.sha256()

    def read(self, n=-1):
        block = self.f.read(n)
        self.h.update(block)
        return block


def write_part(target, compress, top, inner, members):
    """One archive: the inner list, then each (path, bytes, sha256, real). The source is hashed as it is read."""
    with open(target, "wb") as raw:
        if compress == "gz":
            packed = gzip.GzipFile(filename="", fileobj=raw, mode="wb", compresslevel=6, mtime=0)
        elif compress == "xz":
            packed = lzma.LZMAFile(raw, "wb", preset=6)
        else:
            packed = raw
        try:
            with tarfile.open(fileobj=packed, mode="w", format=tarfile.PAX_FORMAT) as tar:
                tar.copybufsize = BLOCK
                tar.addfile(tar_info(f"{top}/{INNER_LIST}", len(inner)), io.BytesIO(inner))
                for path, size, digest, real in members:
                    try:
                        fd = os.open(real, os.O_RDONLY | os.O_NOFOLLOW)
                    except OSError as e:
                        raise Stop(2, f"pack: cannot open {path}: {e}") from None
                    with os.fdopen(fd, "rb") as f:
                        st = os.fstat(f.fileno())
                        if not stat.S_ISREG(st.st_mode) or st.st_size != size:
                            raise Stop(2, f"pack: {path} changed since it was listed: the list has {size} bytes "
                                          f"and the file is {st.st_size} bytes of mode {stat.filemode(st.st_mode)}")
                        reader = Hashing(f)
                        mode = 0o755 if st.st_mode & 0o111 else 0o644
                        try:
                            tar.addfile(tar_info(f"{top}/files/{path}", size, st.st_mtime, mode), reader)
                        except OSError as e:
                            raise Stop(2, f"pack: {path} changed since it was listed, or the part cannot be "
                                          f"written: {e}") from None
                        if reader.h.hexdigest() != digest or f.read(1):
                            raise Stop(2, f"pack: {path} changed since it was listed: its sha256 is no longer {digest}")
        finally:
            if packed is not raw:
                packed.close()


def cmd_pack(a):
    cls, name = a.cls, a.name
    if cls not in CLASSES:
        raise Stop(2, f"pack: unknown class {cls!r}; the classes are {', '.join(CLASSES)}")
    if cls not in PACKABLE:
        raise Stop(2, f"pack: {cls} is never packed")
    if unsafe(name) or "/" in name:
        raise Stop(2, f"pack: --name must be one plain path component, read {name!r}")
    if barred_word(name):
        raise Stop(3, f"pack: refused: the name {name} holds the barred word {barred_word(name)!r}")
    roots = parse_named(a.root, "--root", nested=True)
    for label, top in roots.items():
        word = barred_word(label) or barred_on_disk(top)
        if word:
            raise Stop(3, f"pack: refused: the root {label}={top} holds the barred word {word!r}")
    rows = [r for r in read_list(a.list, LIST_HEADER) if r[3] == cls]
    if not rows:
        raise Stop(2, f"pack: {a.list} holds no row of class {cls}")
    members = []
    for path, size, digest, _ in sorted(rows):
        label = max((r for r in roots if path == r or path.startswith(r + "/")), key=len, default=None)
        if label is None:
            raise Stop(2, f"pack: no --root LABEL=PATH holds {path}")
        real = os.path.join(roots[label], *path[len(label) + 1:].split("/")) if path != label else roots[label]
        word = barred_on_disk(real)
        if word:
            raise Stop(3, f"pack: refused: {path} is at {real}, which holds the barred word {word!r}")
        members.append((path, size, digest, real))
    if cls == PUBLIC:
        if not a.published:
            raise Stop(2, f"pack: {PUBLIC} needs --published SHA256FILE: a file is public only if a published "
                          f"checksum file names its sha256 and its name")
        published = {(digest, base_name(n)) for p in a.published for digest, n in read_hashes(p)}
        missing = [m[0] for m in members if (m[2], base_name(m[0])) not in published]
        if missing:
            raise Stop(2, f"pack: {len(missing)} files listed as {PUBLIC} stand in no --published file with their "
                          f"sha256 and their name, so they are not packed as public:\n{shown(missing)}")
    top = f"{name}.{cls}"
    # Parts are planned on tar bytes. The list's size does not depend on the plan: a part number has three digits.
    inner_size = len(list_text([m[:3] + (cls, "000") for m in members], INNER_HEADER).encode("utf-8"))
    room = a.max_bytes
    while room > 0 and packed_bound(room, a.compress) > a.max_bytes:
        room -= packed_bound(room, a.compress) - a.max_bytes
    room -= member_cost(f"{top}/{INNER_LIST}", inner_size) + TAR_END
    plan, used = [[]], 0
    for m in members:
        cost = member_cost(f"{top}/files/{m[0]}", m[1])
        if cost > room:
            raise Stop(2, f"pack: {m[0]} needs {cost} bytes with its header and a part has room for {max(room, 0)} "
                          f"(--max-bytes {a.max_bytes}, less the list and the tar end); nothing was written")
        if used + cost > room:
            plan.append([])
            used = 0
        plan[-1].append(m)
        used += cost
    if len(plan) > 999:
        raise Stop(2, f"pack: {len(plan)} parts, and a part number has three digits; raise --max-bytes")
    part_of = {m[0]: f"{i:03d}" for i, part in enumerate(plan, 1) for m in part}
    inner = list_text([m[:3] + (cls, part_of[m[0]]) for m in members], INNER_HEADER).encode("utf-8")
    try:
        os.makedirs(a.out, exist_ok=True)
    except OSError as e:
        raise Stop(2, f"pack: cannot make {a.out}: {e}") from None
    targets = [os.path.join(a.out, f"{top}.part{i:03d}{EXTENSION[a.compress]}") for i in range(1, len(plan) + 1)]
    sums = os.path.join(a.out, top + PARTS_SUFFIX)
    for t in targets + [sums]:
        if os.path.lexists(t):
            raise Stop(2, f"pack: {t} exists; nothing was written")
    written = []
    try:
        for target, part in zip(targets, plan):
            written.append(target)
            try:
                write_part(target, a.compress, top, inner, part)
            except OSError as e:
                raise Stop(2, f"pack: cannot write {target}: {e}") from None
            size = os.path.getsize(target)
            if size > a.max_bytes:
                raise Stop(2, f"pack: {target} came to {size} bytes, over --max-bytes {a.max_bytes}")
        digests = in_threads(sha256_file, targets, 1)
        written.append(sums)
        write_new(sums, "".join(f"{d}  {os.path.basename(t)}\n" for t, (d, _) in zip(targets, digests)).encode())
    except BaseException:
        for t in written:
            if os.path.lexists(t):
                os.unlink(t)
        raise
    print("part\tbytes\tsha256\tfiles")
    for target, part, (digest, size) in zip(targets, plan, digests):
        print(f"{os.path.basename(target)}\t{size}\t{digest}\t{len(part)}")
    print(f"pack: {len(members)} files of class {cls}, {sum(m[1] for m in members)} bytes, in {len(plan)} parts; "
          f"their sha256 are in {sums}")


# ---- check ----

class Found:
    """What a check finds: files is {path: (bytes, sha256, part or None)}; the rest should not be there."""

    def __init__(self):
        self.files, self.lists, self.faults = {}, {}, []
        self.refused, self.odd, self.stray = [], [], []


def find_in_folder(top, threads):
    found = Found()
    inner = os.path.join(top, INNER_LIST)
    if not os.path.isfile(inner) or os.path.islink(inner):
        raise Stop(2, f"check: {inner} is not there: give the unpacked NAME.CLASS folder, or --parts")
    found.lists[top] = read_text(inner, "the list")
    paths = {}

    def walk(real, path):
        with os.scandir(real) as it:
            for e in sorted(it, key=lambda e: e.name):
                child = f"{path}/{e.name}" if path else e.name
                if barred_word(child):
                    found.refused.append(child)
                elif e.is_dir(follow_symlinks=False):
                    walk(e.path, child)
                elif e.is_file(follow_symlinks=False):
                    paths[child] = e.path
                else:
                    found.odd.append(child)

    files = os.path.join(top, "files")
    try:
        if os.path.isdir(files) and not os.path.islink(files):
            walk(files, "")
        found.stray = [e for e in os.listdir(top) if e not in (INNER_LIST, "files")]
        hashed = in_threads(sha256_file, list(paths.values()), threads)
    except OSError as e:
        raise Stop(2, f"check: cannot read {top}: {e}") from None
    found.files = {p: (size, digest, None) for p, (digest, size) in zip(paths, hashed)}
    return found


class TarBytes:
    """The tar bytes of one part, for tarfile's stream mode: its one gzip or xz stream decompressed, or the part as it
    is. The stream is read to its end, where zlib and lzma read its check, and Damaged is raised if it stops short or
    if any byte follows it. The bytes last given out are kept, so what follows tarfile's last member can be read."""

    KEEP = 4 * BLOCK

    def __init__(self, raw, compress):
        self.raw, self.compress = raw, compress
        self.unz = (zlib.decompressobj(16 + zlib.MAX_WBITS) if compress == "gz"
                    else lzma.LZMADecompressor(lzma.FORMAT_XZ) if compress == "xz" else None)
        self.full = False  # zlib gave a whole block last time, so it may hold more output
        self.pending, self.kept, self.kept_at = b"", bytearray(), 0

    def raw_block(self):
        data = self.raw.read(BLOCK)
        if not data:
            raise Damaged("the stream stops before its end")
        return data

    def more(self):
        """The next bytes of tar, b"" at the end."""
        unz = self.unz
        if unz is None:
            return self.raw.read(BLOCK)
        while not unz.eof:
            if self.compress == "gz":
                data = unz.unconsumed_tail or (b"" if self.full else self.raw_block())
                out = unz.decompress(data, BLOCK)
                self.full = len(out) == BLOCK
            else:
                out = unz.decompress(self.raw_block() if unz.needs_input else b"", BLOCK)
            if out:
                return out
        if unz.unused_data or getattr(unz, "unconsumed_tail", b"") or self.raw.read(1):
            raise Damaged("bytes follow the end of the stream")
        return b""

    def read(self, n=BLOCK):
        if n is None or n < 0:
            n = BLOCK
        data = self.pending or self.more()
        data, self.pending = data[:n], data[n:]
        self.kept += data
        if len(self.kept) > 2 * self.KEEP:
            cut = len(self.kept) - self.KEEP
            del self.kept[:cut]
            self.kept_at += cut
        return data

    def end(self, offset):
        """Reads from offset, where tarfile found no next member, to the end: only zero bytes, at least the two
        blocks of the tar end mark."""
        if offset < self.kept_at:
            raise Damaged("its last member could not be found again")
        data, zeros = bytes(self.kept[offset - self.kept_at:]), 0
        while data:
            if data.count(0) != len(data):
                raise Damaged("bytes that are not zero follow its last member")
            zeros += len(data)
            data = self.read(BLOCK)
        if zeros < 2 * tarfile.BLOCKSIZE:
            raise Damaged("the tar end mark is missing or cut short")


def read_part(job):
    """One packed part, read as a stream and never unpacked: its sha256, its list, and each file's size and sha256."""
    file, top, compress = job
    out = {"sha256": None, "list": None, "files": [], "faults": [], "refused": [], "odd": [], "stray": []}
    under = f"{top}/files/"
    try:
        out["sha256"] = sha256_file(file)[0]
        with open(file, "rb") as raw:
            stream = TarBytes(raw, compress)
            with tarfile.open(fileobj=stream, mode="r|", bufsize=BLOCK) as tar:
                for m in tar:
                    if m.name == f"{top}/{INNER_LIST}" and m.isreg() and out["list"] is None:
                        out["list"] = tar.extractfile(m).read().decode("utf-8")
                    elif not m.name.startswith(under):
                        out["stray"].append(m.name)
                    elif barred_word(m.name):
                        out["refused"].append(m.name[len(under):])
                    elif not m.isreg():
                        out["odd"].append(m.name[len(under):])
                    else:
                        digest, size = sha256_stream(tar.extractfile(m))
                        out["files"].append((m.name[len(under):], size, digest))
                stream.end(tar.offset)
    except (Damaged, OSError, EOFError, UnicodeDecodeError, tarfile.TarError, zlib.error, lzma.LZMAError) as e:
        why = str(e) if isinstance(e, Damaged) else f"{type(e).__name__}: {e}"
        out["faults"].append(f"the part {os.path.basename(file)} fails when read to the end of its "
                             f"{STREAM[compress]} stream: {why}")
    return out


def find_in_parts(sums, threads):
    found = Found()
    base = os.path.basename(sums)
    if not base.endswith(PARTS_SUFFIX):
        raise Stop(2, f"check: --parts takes a NAME.CLASS{PARTS_SUFFIX} file, read {base}")
    top = base[:-len(PARTS_SUFFIX)]
    parts = {}
    for digest, name in read_hashes(sums):
        number = re.fullmatch(re.escape(top) + r"\.part([0-9]{3})\.tar(\.gz|\.xz)?", name)
        if not number or number.group(1) in parts:
            raise Stop(2, f"check: {sums}: {name!r} is not one part of {top}, or its number is there twice")
        compress = {".gz": "gz", ".xz": "xz", None: "none"}[number.group(2)]
        parts[number.group(1)] = (os.path.join(os.path.dirname(sums), name), digest, compress)
    if not parts:
        raise Stop(2, f"check: {sums} names no part")
    present = {n: p for n, p in parts.items() if os.path.isfile(p[0]) and not os.path.islink(p[0])}
    for n in sorted(set(parts) - set(present)):
        found.faults.append(f"the part {os.path.basename(parts[n][0])} is not beside {base}")
    results = in_threads(read_part, [(present[n][0], top, present[n][2]) for n in sorted(present)], threads)
    for n, r in zip(sorted(present), results):
        name = os.path.basename(present[n][0])
        found.faults += r["faults"]
        if r["sha256"] is not None and r["sha256"] != present[n][1]:
            found.faults.append(f"the part {name} has the sha256 {r['sha256']}, and {base} has {present[n][1]}")
        if r["list"] is None:
            found.faults.append(f"the part {name} holds no {top}/{INNER_LIST}")
        else:
            found.lists[name] = r["list"]
        for path, size, digest in r["files"]:
            if path in found.files:
                found.faults.append(f"{path} is in the archive twice")
            found.files[path] = (size, digest, n)
        found.refused += r["refused"]
        found.odd += r["odd"]
        found.stray += r["stray"]
    return found, top, sorted(parts)


def cmd_check(a):
    if (a.dir is None) == (a.parts is None):
        raise Stop(2, "check: give one unpacked folder, or --parts NAME.CLASS.parts.sha256")
    if a.parts:
        found, where, parts = find_in_parts(a.parts, a.threads)
    else:
        found, where, parts = find_in_folder(a.dir, a.threads), a.dir, None
    faults = found.faults
    if not found.lists:
        for f in faults:
            print(f"check: {f}")
        print(f"check: FAILED {where}: no list could be read")
        return 1
    if len(set(found.lists.values())) > 1:
        faults.append(f"the parts do not hold the same {INNER_LIST}: {', '.join(sorted(found.lists))}")
    first = sorted(found.lists)[0]
    rows = parse_list(found.lists[first], INNER_HEADER, f"{first}: {INNER_LIST}")
    classes = sorted({r[3] for r in rows})
    if len(classes) != 1 or classes[0] not in PACKABLE:
        faults.append(f"the list names the classes {', '.join(classes) or 'none'}; an archive holds one packed class")
    elif parts is not None and not where.endswith("." + classes[0]):
        faults.append(f"the list is of class {classes[0]} and the archive is named {where}")
    elif parts is None and not os.path.basename(os.path.normpath(where)).endswith("." + classes[0]):
        faults.append(f"the list is of class {classes[0]} and the folder is named "
                      f"{os.path.basename(os.path.normpath(where))}")
    if parts is not None:
        unnamed = sorted({r[4] for r in rows} - set(parts))
        if unnamed:
            faults.append(f"the list has files in the parts {', '.join(unnamed)}, which "
                          f"{os.path.basename(a.parts)} does not name")
    listed = {r[0]: r for r in rows}
    missing = [r for r in rows if r[0] not in found.files]
    extra = sorted(set(found.files) - set(listed))
    wrong_size, wrong_hash, wrong_part = [], [], []
    for path, (size, digest, part) in found.files.items():
        if path not in listed:
            continue
        r = listed[path]
        if size != r[1]:
            wrong_size.append(f"{path} (the list has {r[1]} bytes, the file has {size})")
        elif digest != r[2]:
            wrong_hash.append(path)
        elif part is not None and part != r[4]:
            wrong_part.append(f"{path} (the list has part {r[4]}, the file is in part {part})")
    if found.refused:
        faults.append(f"{len(found.refused)} names hold a barred word; none was opened:\n{shown(found.refused)}")
    if missing:
        by_part = {}
        for r in missing:
            by_part[r[4]] = by_part.get(r[4], 0) + 1
        of = ", ".join(f"{n} of part {p}" for p, n in sorted(by_part.items()))
        faults.append(f"{len(missing)} listed files are missing ({of}):\n{shown(r[0] for r in missing)}")
    if wrong_size:
        faults.append(f"{len(wrong_size)} files have another size:\n{shown(wrong_size)}")
    if wrong_hash:
        faults.append(f"{len(wrong_hash)} files have the listed size and another sha256:\n{shown(wrong_hash)}")
    if wrong_part:
        faults.append(f"{len(wrong_part)} files are in another part than listed:\n{shown(wrong_part)}")
    if extra:
        faults.append(f"{len(extra)} files are not in the list:\n{shown(extra)}")
    if found.odd:
        faults.append(f"{len(found.odd)} links or special files, which an archive does not hold:\n{shown(found.odd)}")
    if found.stray:
        faults.append(f"{len(found.stray)} entries beside {INNER_LIST} and files/:\n{shown(found.stray)}")
    if faults:
        good = len(rows) - len(missing) - len(wrong_size) - len(wrong_hash) - len(wrong_part)
        for f in faults:
            print(f"check: {f}")
        print(f"check: FAILED {where}: {good} of {len(rows)} listed files are as listed")
        return 1
    held = f"{len(parts)} parts with the listed sha256, " if parts is not None else ""
    print(f"check: OK {where}: {held}{len(rows)} files, {sum(r[1] for r in rows)} bytes, class {classes[0]}, "
          f"every size and sha256 as listed, nothing else under files/")
    return 0


# ---- summary ----

def cmd_summary(a):
    places = parse_named(a.place, "--place")
    extra = parse_named(a.hashes, "--hashes")
    if set(places) & set(extra):
        raise Stop(2, f"summary: {sorted(set(places) & set(extra))} is both a --place and a --hashes name")
    rows = []  # (place, path, bytes, sha256, class)
    holders = {}  # sha256 -> the places that hold it
    for place, file in places.items():
        for path, size, digest, cls in read_list(file, LIST_HEADER):
            rows.append((place, path, size, digest, cls))
            if digest != "-":
                holders.setdefault(digest, set()).add(place)
    named_only = {}
    for place, file in extra.items():
        named = {digest for digest, _ in read_hashes(file)}
        named_only[place] = (len(named), sum(d not in holders for d in named))
        for d in named & set(holders):
            holders[d].add(place)
    out = []

    def table(title, header, lines):
        out.append(f"# {title}")
        out.append("\t".join(header))
        out.extend("\t".join(map(str, ln)) for ln in lines)
        out.append("")

    def distinct(of):
        """{sha256: bytes} of the hashed rows."""
        return {r[3]: r[2] for r in of if r[3] != "-"}

    lines = []
    for cls in CLASSES + ("all",):
        of = [r for r in rows if cls in (r[4], "all")]
        if not of:
            continue
        d = distinct(of)
        one = {k: v for k, v in d.items() if len(holders[k]) == 1}
        bare = [r for r in of if r[3] == "-"]
        lines.append((cls, len(of), sum(r[2] or 0 for r in of), len(d), sum(d.values()), len(one), sum(one.values()),
                      len(bare), sum(r[2] or 0 for r in bare)))
    table("By class. A content is one sha256; entries with no sha256 (regenerable directories and files, refused "
          "paths) are counted apart. The row all removes duplicates across classes too.",
          ("class", "entries", "bytes", "contents", "bytes_without_duplicates", "contents_in_one_place",
           "bytes_in_one_place", "entries_not_hashed", "bytes_not_hashed"), lines)
    lines = []
    for cls in CLASSES + ("all",):
        groups = {}
        for k, size in distinct([r for r in rows if cls in (r[4], "all")]).items():
            g = groups.setdefault(tuple(sorted(holders[k])), [0, 0])
            g[0] += 1
            g[1] += size
        lines += [(cls, len(held), "+".join(held), n, size) for held, (n, size) in sorted(groups.items())]
    table("By class and by the places that hold a content.", ("class", "places", "held_by", "contents", "bytes"), lines)
    lines = []
    for place in places:
        for cls in CLASSES:
            of = [r for r in rows if r[0] == place and r[4] == cls]
            if of:
                lines.append((place, cls, len(of), sum(r[2] or 0 for r in of),
                              sum(distinct(of).values()) + sum(r[2] or 0 for r in of if r[3] == "-")))
    table("By place and class.", ("place", "class", "entries", "bytes", "bytes_without_duplicates"), lines)
    table("Places known by their sha256 only.", ("place", "sha256_named", "in_no_list"),
          [(place, n, absent) for place, (n, absent) in named_only.items()])
    print("\n".join(out), end="")
    if a.out:
        body = ["\t".join(("place", "path", "bytes", "sha256", "class", "places", "held_by"))]
        for place, path, size, digest, cls in sorted(rows, key=lambda r: r[:2]):
            held = sorted(holders[digest]) if digest != "-" else []
            body.append("\t".join((place, path, "-" if size is None else str(size), digest, cls,
                                   str(len(held)) if held else "-", "+".join(held) or "-")))
        write_new(a.out, ("\n".join(body) + "\n").encode("utf-8"))
        print(f"summary: wrote {len(rows)} rows to {a.out}")


# ---- command line ----

def parser():
    p = argparse.ArgumentParser(prog="collect.py", description="List, pack and check a collection of files.")
    sub = p.add_subparsers(dest="command", required=True)
    s = sub.add_parser("list")
    s.add_argument("--rules", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--root", action="append", metavar="LABEL=PATH")
    s.add_argument("--listing", action="append", metavar="FILE")
    s.add_argument("--strip-prefix", default="", metavar="P")
    s.add_argument("--record-refused", action="store_true")
    s.add_argument("--skip-special", action="store_true")
    s.add_argument("--threads", type=int, default=4)
    s.set_defaults(run=cmd_list, paths=lambda a: [a.rules, a.out] + (a.listing or []))
    s = sub.add_parser("pack")
    s.add_argument("--list", required=True)
    s.add_argument("--class", dest="cls", required=True)
    s.add_argument("--name", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--root", action="append", metavar="LABEL=PATH")
    s.add_argument("--published", action="append", metavar="SHA256FILE")
    s.add_argument("--max-bytes", type=int, default=PART_BYTES)
    s.add_argument("--compress", choices=tuple(EXTENSION), default="gz")
    s.set_defaults(run=cmd_pack, paths=lambda a: [a.list, a.out] + (a.published or []))
    s = sub.add_parser("check")
    s.add_argument("dir", metavar="DIR", nargs="?")
    s.add_argument("--parts", metavar="PARTS_SHA256")
    s.add_argument("--threads", type=int, default=4)
    s.set_defaults(run=cmd_check, paths=lambda a: [p for p in (a.dir, a.parts) if p])
    s = sub.add_parser("summary")
    s.add_argument("--place", action="append", required=True, metavar="NAME=LIST")
    s.add_argument("--hashes", action="append", metavar="NAME=SHA256FILE")
    s.add_argument("--out")
    s.set_defaults(run=cmd_summary, paths=lambda a: [v.partition("=")[2] for v in a.place + (a.hashes or [])] +
                   ([a.out] if a.out else []))
    return p


def main(argv=None):
    a = parser().parse_args(argv)
    try:
        # A root is refused inside its command, so that `list` can name it with the rest; every other path here.
        for path in a.paths(a):
            word = barred_on_disk(path)
            if word:
                raise Stop(3, f"{a.command}: refused: the argument {path} holds the barred word {word!r}")
        return a.run(a) or 0
    except Stop as e:
        print(e, file=sys.stderr)
        return e.status


if __name__ == "__main__":
    sys.exit(main())
