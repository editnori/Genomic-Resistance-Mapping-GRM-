"""H-1's push-data manifests (G5, G11), each derived on the laptop from a pinned source; every line is
"<sha256>  <relative path>", sorted by path. Nothing is written into the source trees.
  meta MASTER                the master manifest's own lines (verbatim) for bvbrc/genomes.jsonl, bvbrc/genome_amr.jsonl,
                             catalog/* and species-check/* (and a synthetic root's receipts/truth.json): what every root
                             curates from (no features, curated/,
                             rev9/ or bindings/, RUN-PLAN-H1 setup step 5)
  fasta MASTER UNION         the master's line for fasta/<id>.fna.gz of every id in UNION (the box's union.txt); an id
                             without one is an error (the run would refuse at the pool)
  vocab ROOT                 ROOT/vocabulary/*: GRM's SigDrug vocabulary copy (D6) as vocabulary/<name>.tsv lines,
                             hashed here; regular files only, no subfolder or link
  combine OUT MASTER SHA256 PART=DIR...
                             builds OUT (new): parts/<part>/h1work/ hard links of each pulled WORK's records
                             (top-level .json and .txt, arms/, j/, k/, report/; never roots/ or logs/) and a copy of
                             MASTER, which must hash to SHA256; prints the manifest of OUT
Errors: exit 1, the message on the last line of stdout.
usage: python3 -I h1-manifest.py meta|fasta|vocab|combine ..."""
import hashlib
import os
import re
import shutil
import stat
import sys


def die(msg):
    print(msg)
    sys.exit(1)


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def master(path):
    lines = {}
    for n, line in enumerate(open(path, encoding="utf-8").read().splitlines(), 1):
        if not re.fullmatch(r"[0-9a-f]{64}  [^\n]+", line):
            die("master manifest line %d is not '<sha256>  <path>'" % n)
        lines[line[66:]] = line
    return lines


def emit(lines):
    for line in sorted(lines, key=lambda x: x[66:]):
        print(line)


def tree(top):
    """Every regular file under top, relative, refusing links and special files."""
    out = []
    for here, dirs, files in os.walk(top, followlinks=False):
        dirs.sort()
        for n in sorted(dirs + files):
            st = os.lstat(os.path.join(here, n))
            if stat.S_ISLNK(st.st_mode) or not (stat.S_ISREG(st.st_mode) or stat.S_ISDIR(st.st_mode)):
                die("%s is a link or special file" % os.path.join(here, n))
        out += [os.path.relpath(os.path.join(here, n), top) for n in sorted(files)]
    return out


def main(a):
    if a[:1] == ["meta"] and len(a) == 2:
        m = master(a[1])
        # receipts/truth.json: a synthetic root's marker (synthetic_root.py's truth.json, which run.py reads for a
        # spec that says "synthetic": true); h1-box.sh puts it back at RAW/truth.json only for such a spec
        take = [m[p] for p in m if p in ("bvbrc/genomes.jsonl", "bvbrc/genome_amr.jsonl", "receipts/truth.json")
                or p.split("/")[0] in ("catalog", "species-check")]
        for need in ("bvbrc/genomes.jsonl", "bvbrc/genome_amr.jsonl"):
            if need not in m:
                die("the master manifest has no %s" % need)
        if not any(x[66:].startswith("species-check/") for x in take):
            die("the master manifest has no species-check/")
        emit(take)
    elif a[:1] == ["fasta"] and len(a) == 3:
        m = master(a[1])
        ids = [x for x in open(a[2], encoding="utf-8").read().split()]
        if not ids or len(set(ids)) != len(ids):
            die("the union is empty or lists an id twice")
        bad = [i for i in ids if not re.fullmatch(r"[0-9]+\.[0-9]+", i)]
        if bad:
            die("%d union ids are not BV-BRC genome ids, first %r" % (len(bad), bad[0]))
        miss = [i for i in ids if "fasta/%s.fna.gz" % i not in m]
        if miss:
            die("%d of the union's %d genomes have no FASTA line in the master manifest, first %s" % (len(miss), len(ids), miss[0]))
        emit(m["fasta/%s.fna.gz" % i] for i in ids)
    elif a[:1] == ["vocab"] and len(a) == 2:
        v = os.path.join(a[1], "vocabulary")
        if os.path.islink(v) or not os.path.isdir(v):
            die("%s is not a folder" % v)
        files = tree(v)
        if not files or any("/" in f for f in files):
            die("vocabulary/ must hold files only, no subfolder")
        emit("%s  vocabulary/%s" % (sha(os.path.join(v, f)), f) for f in files)
    elif a[:1] == ["combine"] and len(a) >= 5:
        out, mpath, want, parts = a[1], a[2], a[3], a[4:]
        if os.path.exists(out):
            die("%s exists; the drive makes it fresh" % out)
        if sha(mpath) != want:
            die("master manifest %s is not %s" % (mpath, want))
        if not re.fullmatch(r"MASTER-MANIFEST[A-Za-z0-9._-]*\.sha256", os.path.basename(mpath)):
            die("the master manifest's name %s is not MASTER-MANIFEST*.sha256, the runner's class" % os.path.basename(mpath))
        lines = []
        os.makedirs(out)
        shutil.copyfile(mpath, os.path.join(out, os.path.basename(mpath)))
        if sha(os.path.join(out, os.path.basename(mpath))) != want:
            die("the master manifest's copy does not hash to %s" % want)
        lines.append("%s  %s" % (want, os.path.basename(mpath)))
        for spec in parts:
            name, _, d = spec.partition("=")
            if not re.fullmatch(r"[a-z0-9][a-z0-9._-]*", name) or not os.path.isdir(d):
                die("bad part %r (want name=DIR, name lowercase)" % spec)
            for f in tree(d):
                top = f.split("/")[0]
                if "/" in f and top not in ("arms", "j", "k", "report"):
                    continue
                if "/" not in f and not f.endswith((".json", ".txt")):
                    continue
                dst = os.path.join(out, "parts", name, "h1work", f)
                os.makedirs(os.path.dirname(dst), exist_ok=True)
                os.link(os.path.join(d, f), dst)
                lines.append("%s  parts/%s/h1work/%s" % (sha(dst), name, f))
        emit(lines)
    else:
        print(__doc__)
        sys.exit(2)


main(sys.argv[1:])
