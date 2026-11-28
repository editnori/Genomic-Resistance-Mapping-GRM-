#!/usr/bin/env python3
"""Check web/src/flow/flow.json against the files its counts were read from.

    python3 -I scripts/check-flow.py
    python3 -I scripts/check-flow.py --repo . --allow-missing-lead
    python3 -I scripts/check-flow.py --components

Every count is read back from its source at the pinned commit via `git show`
(never from the work tree), and lead-rooted citations are checked against the
lead directory when present. Exits 1 on any difference.
"""
import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FLOW = REPO_ROOT / "web" / "src" / "flow" / "flow.json"
COMPONENTS_DIR = REPO_ROOT / "web" / "src" / "flow"

# Paths this work must never open. A source that matches one is refused.
NEVER = [
    "/run/", "/home/user/scratch", "/latetest/", "/back/", "/score/",
    "/set/", ".tar", ".gz", ".jsonl", "timesplit/REPORT.md", "species-check.md",
    "outside-study", "sampling-audit/README.md", "MANIFEST", ".sha256", "owner2/work",
]

# Words the writing guide says not to use in a string the reader sees.
AVOID = [
    "held", "spot checked", "lane", "census", "receipt", "provenance", "landed",
    "cleanly", "trap", "snag", "breakthrough", "synced", "frame", "keeper",
]

# Shapes that would be a genome id, an accession or a sequence.
ID_SHAPES = [
    (r"\bGC[AF]_\d{6,}", "assembly accession"),
    (r"\b(?:SAM[NED][A-Z]?\d{5,}|[SED]RR\d{5,}|PRJ[NED][A-Z]\d{3,})", "archive accession"),
    (r"\b[A-Z]{1,2}_?\d{6,}\.\d\b", "sequence accession"),
    (r"\b\d{3,7}\.\d{2,6}\b", "genome id shape"),
    (r"[ACGT]{12,}", "sequence"),
]

WORDS = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6,
    "seven": 7, "eight": 8, "nine": 9, "ten": 10, "neither": 0, "no": 0,
}

PLACE = re.compile(r"\{([a-z0-9_]+)\}")

errors = []
skipped_lead = []
_cache = {}


def bad(msg):
    errors.append(msg)


def resolve_root(flow, root, repo_override, lead_override):
    if root == "repo":
        if repo_override:
            return str(Path(repo_override).resolve())
        configured = Path(flow["roots"]["repo"])
        if configured.exists():
            return str(configured)
        return str(REPO_ROOT)
    if root == "lead":
        if lead_override:
            return str(Path(lead_override).resolve())
        return flow["roots"]["lead"]
    return flow["roots"][root]


def read_source(flow, src, repo_override, lead_override, allow_missing_lead):
    """Lines of a source file. 'repo' is read at the pinned commit, never from the work tree."""
    root, rel = src["root"], src["file"]
    key = (root, rel)
    if key in _cache:
        return _cache[key]
    base = resolve_root(flow, root, repo_override, lead_override)
    full = f"{base}/{rel}"
    for n in NEVER:
        if n in full or full.endswith(n):
            raise SystemExit(f"REFUSED: {full} is on the never-open list ({n})")
    if root == "repo":
        out = subprocess.run(
            ["git", "-C", base, "show", f"{flow['pinned']['commit']}:{rel}"],
            capture_output=True, text=True, encoding="utf-8",
        )
        if out.returncode != 0:
            bad(f"{rel}: not in commit {flow['pinned']['commit'][:8]}")
            _cache[key] = None
            return None
        text = out.stdout
    else:
        p = Path(full)
        if not p.is_file():
            if root == "lead" and allow_missing_lead:
                skipped_lead.append(rel)
                _cache[key] = None
                return None
            bad(f"{full}: no such file")
            _cache[key] = None
            return None
        text = p.read_text(encoding="utf-8")
    _cache[key] = text.split("\n")
    return _cache[key]


def numbers_in(text):
    found = set()
    for m in re.finditer(r"\d[\d,]*(?:\.\d+)?", text):
        s = m.group(0).rstrip(",").replace(",", "")
        try:
            found.add(float(s))
        except ValueError:
            pass
    return found


def line_of(lines, n, where):
    if lines is None:
        return None
    if not (1 <= n <= len(lines)):
        bad(f"{where}: line {n} is past the end of the file ({len(lines)} lines)")
        return None
    return lines[n - 1]


def check_src_text(flow, src, where, repo_override, lead_override, allow_missing_lead):
    lines = read_source(flow, src, repo_override, lead_override, allow_missing_lead)
    line = line_of(lines, src["line"], where)
    if line is None:
        return False
    if src["text"] not in line:
        hits = [i + 1 for i, l in enumerate(lines) if src["text"] in l]
        bad(
            f"{where}: {src['file']}:{src['line']} does not hold \"{src['text']}\""
            + (f" (found on line {hits})" if hits else " (found on no line)")
            + f"\n      line reads: {line[:200]}"
        )
        return False
    return True


def check_count(flow, name, c, repo_override, lead_override, allow_missing_lead):
    where = f"count {name}"
    src = c.get("src")
    if not src:
        bad(f"{where}: no source, so it stays out")
        return
    v = c["v"]
    lines = read_source(flow, src, repo_override, lead_override, allow_missing_lead)
    if lines is None:
        return
    if "key" in src:
        try:
            node = json.loads("\n".join(lines))
            for k in src["key"]:
                node = node[k]
        except (KeyError, ValueError) as e:
            bad(f"{where}: key path {src['key']} not in {src['file']} ({e})")
            return
        if node != v:
            bad(f"{where}: flow.json says {v}, {src['file']} says {node}")
        line = line_of(lines, src["line"], where)
        if line is not None and (f"\"{src['key'][-1]}\"" not in line or str(v) not in line):
            bad(
                f"{where}: {src['file']}:{src['line']} is not the line of {src['key'][-1]} = {v}"
                f"\n      line reads: {line[:160]}"
            )
    elif "col" in src:
        head = lines[0].split("\t")
        line = line_of(lines, src["line"], where)
        if line is None:
            return
        row = dict(zip(head, line.split("\t")))
        for k, want in src["where"].items():
            if row.get(k) != want:
                bad(f"{where}: {src['file']}:{src['line']} has {k}={row.get(k)!r}, not {want!r}")
        got = row.get(src["col"])
        if got != str(v):
            bad(f"{where}: flow.json says {v}, {src['file']}:{src['line']} column {src['col']} says {got}")
    elif "distinct" in src:
        head = lines[0].split("\t")
        if src["distinct"] not in head:
            bad(f"{where}: {src['file']} has no column {src['distinct']} (it has {head})")
            return
        i = head.index(src["distinct"])
        vals = {l.split("\t")[i] for l in lines[1:] if l.strip()}
        if len(vals) != v:
            bad(f"{where}: flow.json says {v}, {src['file']} has {len(vals)} distinct {src['distinct']}: {sorted(vals)}")
    elif "text" in src:
        if not check_src_text(flow, src, where, repo_override, lead_override, allow_missing_lead):
            return
        if "regex_count" in src:
            n = sum(1 for l in lines if re.search(src["regex_count"], l))
            if n != v:
                bad(f"{where}: flow.json says {v}, {src['file']} has {n} lines matching {src['regex_count']}")
        elif "word" in src:
            if src["word"] not in src["text"] or WORDS.get(src["word"]) != v:
                bad(f"{where}: the word \"{src['word']}\" does not give {v}")
        elif float(v) not in numbers_in(src["text"]):
            bad(f"{where}: {v} is not a number in the quoted phrase \"{src['text']}\"")
    else:
        bad(f"{where}: source has no key, column or text to check")


def label_strings(node):
    out = [("title", node["title"])]
    for k in ("line", "note", "tag"):
        if node.get(k):
            out.append((k, node[k]))
    for name, text in node.get("fields", []):
        out.append((f"field {name}", name))
        out.append((f"field {name}", text))
    for a in node.get("alts", []):
        out.append(("alternative", a[0]))
    for k in ("from", "cardtitle"):
        if node.get(k):
            out.append((k, node[k]))
    for row in node.get("card", []):
        out += [("card", row[0]), ("card", row[2])]
    for row in node.get("bars", []):
        out.append(("bar", row[0]))
    for row in node.get("stats", []):
        out += [("stat", row[0]), ("stat", row[2])]
    return out


def scan_ids(where, text):
    for pat, what in ID_SHAPES:
        m = re.search(pat, text)
        if m:
            bad(f"{where}: looks like a {what}: \"{m.group(0)}\"")


def check_label(flow, where, text):
    for m in PLACE.finditer(text):
        if m.group(1) not in flow["counts"]:
            bad(f"{where}: {{{m.group(1)}}} is no count in flow.json")
    rest = PLACE.sub("", text)
    for lit in flow["literals"]:
        rest = rest.replace(lit["text"], "")
    if re.search(r"\d", rest):
        bad(f"{where}: a number is typed into the label: \"{text}\"")
    low = re.sub(r"[^a-z ]", " ", text.lower())
    for w in AVOID:
        if re.search(rf"\b{w}\b", low):
            bad(f"{where}: the guide's avoided word \"{w}\" is in \"{text}\"")
    scan_ids(where, text)


def check_components():
    """Ensure React flow components hold no real genome IDs, accessions, or sequences."""
    for p in sorted(COMPONENTS_DIR.glob("*.ts*")):
        if p.name.endswith(".test.ts"):
            continue
        text = p.read_text(encoding="utf-8")
        scan_ids(f"component {p.name}", text)


def main():
    parser = argparse.ArgumentParser(description="Verify web/src/flow/flow.json against pinned source files.")
    parser.add_argument("--repo", default=None, help="Path to git clone of Genomic-Resistance-Mapping-GRM-.")
    parser.add_argument("--lead", default=None, help="Path to lead program directory (optional).")
    parser.add_argument(
        "--allow-missing-lead",
        action="store_true",
        help="Skip lead-rooted citation files that do not exist locally (all 136 repo counts still checked).",
    )
    parser.add_argument("--components", action="store_true", help="Also scan React flow components for ID/sequence shapes.")
    args = parser.parse_args()

    flow = json.loads(FLOW.read_text(encoding="utf-8"))
    counts = flow["counts"]
    repo_dir = resolve_root(flow, "repo", args.repo, args.lead)

    head = subprocess.run(
        ["git", "-C", repo_dir, "rev-parse", "HEAD"], capture_output=True, text=True, encoding="utf-8"
    ).stdout.strip()
    moved = head != flow["pinned"]["commit"]

    for name, c in counts.items():
        check_count(flow, name, c, args.repo, args.lead, args.allow_missing_lead)
        if "of" in c:
            w = counts.get(c["of"])
            if w is None:
                bad(f"count {name}: its whole \"{c['of']}\" is no count")
            elif c["v"] > w["v"]:
                bad(f"count {name}: {c['v']} is more than its whole {c['of']} = {w['v']}")

    for s in flow.get("sums", []):
        total = counts[s["total"]]["v"]
        parts = sum(counts[p]["v"] for p in s["parts"])
        if parts != total:
            bad(f"sum {s['total']}: parts add to {parts}, the total is {total}")

    for lit in flow["literals"]:
        check_src_text(flow, lit["src"], f"literal \"{lit['text']}\"", args.repo, args.lead, args.allow_missing_lead)

    ids = set()
    for n in flow["nodes"]:
        where = f"node {n['id']}"
        if n["id"] in ids:
            bad(f"{where}: id used twice")
        ids.add(n["id"])
        if n["kind"] not in flow["kinds"]:
            bad(f"{where}: kind \"{n['kind']}\" is not one of {flow['kinds']}")
        if n["state"] not in flow["states"]:
            bad(f"{where}: state \"{n['state']}\" is not one of {flow['states']}")
        if not n.get("src"):
            bad(f"{where}: no source")
        for s in n.get("src", []):
            check_src_text(flow, s, where, args.repo, args.lead, args.allow_missing_lead)
        if "count" in n and n["count"] not in counts:
            bad(f"{where}: its count \"{n['count']}\" is no count")
        for part, text in label_strings(n):
            check_label(flow, f"{where} {part}", text)
        quotes = {lit["text"] for lit in flow["literals"]}
        for name, text in n.get("fields", []):
            if name == "Asked" and text.strip('"') not in quotes:
                bad(f"{where}: the Asked line is not a quote with a source: {text}")
        for a in n.get("alts", []):
            if a[1] not in ("no", "pass", "open"):
                bad(f"{where}: alternative state \"{a[1]}\" is not no, pass or open")
        for row in n.get("card", []):
            if row[1] not in ("used", "checked", "kept", "unread"):
                bad(f"{where}: card use \"{row[1]}\" is not known")
        for row in n.get("bars", []) + n.get("stats", []):
            if row[1] not in counts:
                bad(f"{where}: \"{row[1]}\" is no count")
        for t in n.get("tools", []):
            out = subprocess.run(
                ["git", "-C", repo_dir, "cat-file", "-e", f"{flow['pinned']['commit']}:{t}"],
                capture_output=True,
            )
            if out.returncode != 0:
                bad(f"{where}: tool path {t} is not in the commit")

    for n in flow["nodes"]:
        for i in [n.get("test"), n.get("step")] + n.get("decisions", []):
            if i and i not in ids:
                bad(f"node {n['id']}: \"{i}\" is no node")
        if n["kind"] == "assumption" and not n.get("mark"):
            bad(f"node {n['id']}: a choice with no mark")

    for e in flow["edges"]:
        for end in ("from", "to"):
            if e[end] not in ids:
                bad(f"edge {e['from']} -> {e['to']}: \"{e[end]}\" is no node")
        if e.get("label"):
            check_label(flow, f"edge {e['from']} -> {e['to']}", e["label"])

    for vname, view in flow["views"].items():
        for k in ("kicker", "title", "sub"):
            check_label(flow, f"view {vname} {k}", view[k])
        for k, text in view.get("ui", {}).items():
            check_label(flow, f"view {vname} ui {k}", text)
        used = list(view.get("steps", []))
        for row in view.get("rows", []):
            used += [i for i in row if i]
        used += list(view.get("next", {}))
        for step, rows in view.get("trunk", []):
            used += [step] + rows
        for i in list(view.get("steps", [])):
            used += next((n.get("decisions", []) for n in flow["nodes"] if n["id"] == i), [])
        for i in used + view.get("walk", []) + [view["start"]]:
            if i not in ids:
                bad(f"view {vname}: \"{i}\" is no node")
        if "trunk" in view:
            rows = [r for _, rr in view["trunk"] for r in rr]
            got = sorted(r for r in rows if r.startswith("as"))
            want = sorted(n["id"] for n in flow["nodes"] if n["kind"] == "assumption")
            if got != want or len(got) != counts["n_assumptions"]["v"]:
                bad(f"view {vname}: the tree names {len(got)} choices, the plan lists {counts['n_assumptions']['v']}")
            bare = sum(1 for n in flow["nodes"] if n.get("bare"))
            if bare != counts["n_undefended"]["v"]:
                bad(f"view {vname}: {bare} choices are marked with no reason or measure, the plan says {counts['n_undefended']['v']}")
        if "walk" in view and sorted(view["walk"]) != sorted(used):
            bad(f"view {vname}: the walk order and the picture do not name the same steps")
        for i in used:
            if i in ids and not next(n for n in flow["nodes"] if n["id"] == i).get("line"):
                bad(f"view {vname}: step {i} has no line to print")

    if args.components:
        check_components()

    n_files = len({k for k, v in _cache.items() if v is not None})
    if errors:
        print(f"STOP: {len(errors)} difference(s)")
        for e in errors:
            print("  -", e)
        return 1
    note_skip = f" (skipped {len(set(skipped_lead))} non-repo lead file(s): {', '.join(sorted(set(skipped_lead)))})" if skipped_lead else ""
    print(
        f"OK: {len(counts)} counts, {len(flow['literals'])} literals, {len(flow['nodes'])} steps and "
        f"{len(flow['edges'])} links agree with {n_files} source files at commit {flow['pinned']['commit'][:8]}{note_skip}"
    )
    if moved:
        print(f"NOTE: the main repository is now at {head[:8]}; this file is pinned to {flow['pinned']['commit'][:8]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
