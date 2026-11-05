"""Sealed lockbox for GRM multi-task development.

The lockbox is the union of recent genomes and whole held-out lineages,
closed over duplicate groups. define_lockbox returns a manifest holding the
definition, the ids, the definition-time duplicate groups of the whole
input, counts, per-task status, the access log's path and a seal: the sha256
of the manifest's canonical JSON. write_lockbox writes it once, read-only,
and starts the access log with a genesis record bound to the seal.

Guards. require_development claims a genome by lockbox id, held-out lineage,
recency or a lockbox duplicate group, or because an earlier release passed
in releases claimed it by one of those rules, and also every genome
connected to a claimed genome through duplicate groups. The connections are
the union of the definition-time groups of the whole input (stored in the
manifest), the groups of every earlier release passed in releases, and the
groups of the input. Group labels are one namespace across releases, so a
reused label links its members conservatively. Every lockbox id, and every
genome an earlier release claimed, seeds the closure even when it has left
the input. A release the manifest has not seen links nothing unless it is
passed in releases: pass every release used for development since
definition. Run development_genomes on the full release and split its
result, never a pre-filtered subset. A duplicate that upstream clustering
left without a shared duplicate_group cannot be detected here.
assert_untouched guards any later use of ids, under the same genome-id rule.

What the seal and log guarantee. The seal detects an edit to a manifest that
keeps its recorded seal, but anyone can reseal an edited manifest. Comparing
with the seal recorded at seal time refuses a resealed forgery: unseal
requires expected_seal, and verify_seal, load_lockbox, require_development,
development_genomes, assert_untouched and splits.verify_manifest accept it.
The access log is a JSONL hash chain without a key, at an absolute path the
manifest records. It detects an edited, reordered or removed record inside
the chain, a missing log and a log started for another seal. On its own it
cannot detect the tail being cut off or the whole chain being rewritten or
re-created. Those need a head (record count and last hash) recorded outside
the box: write_lockbox and every recorded unseal attempt return the new
head, the caller records it (GRM records it in the work-program ledger), and
unseal requires it as expected_head and refuses unless the log's head is
exactly that. verify_log checks a recorded head on demand. There is no file
lock, so concurrent unseals can break the chain, and write_lockbox is not
atomic.

Not access control. The manifest is readable without unseal: it lists the
lockbox genome ids and publishes the lockbox's per-task class counts. unseal
is a procedural record of approved use, not a lock.

Identifiers (rule 6, v2, shared with Tasks M and N). For `lineage` and `duplicate_group`, a non-null value is first stripped of surrounding whitespace with Python `str.strip()`. The stripped value is refused if it's empty, or if it contains any character outside printable ASCII (U+0021 to U+007E) other than an interior space (U+0020). Two distinct values of one field in one input that are equal under `unicodedata.normalize("NFKC", v).casefold()` are refused as an alias collision, naming both. Values are otherwise kept as given. All refusals are `ValueError` naming the genome or the two values.

Task P additions: the lockbox also runs the alias check across its stored
values (lockbox ids, held-out lineages, stored groups and their members),
the input and any earlier releases together. Genome ids get the same printable-ASCII and alias rules but are
never stripped: surrounding whitespace is refused. Hierarchical lineages are
not prefix-matched, so holding out "4" does not claim "4.2.1".

Genome records hold genome_id, lineage, duplicate_group (or None), an
optional collection_date ("YYYY", "YYYY-MM" or "YYYY-MM-DD") and tasks:
either a list of task names, or a mapping of task name to a 0/1 label
(1 = resistant) or None. A task carries 0/1 labels on every genome or on none.
A partial date counts as its latest possible day, so a genome that might be
recent is treated as recent. Undated genomes are never recent; the manifest
counts them.

Development only: never unseal a real lockbox outside an approved
confirmation.
"""
import calendar
import datetime as dt
import errno
import hashlib
import json
import os
import re
import unicodedata

KIND = "grm.multitask.lockbox"
VERSION = 3
LOG_GENESIS = "0" * 64
RECENT_RULE = "latest possible collection day after recent_after; undated genomes are not recent"
# Placeholders only: set both per protocol for each confirmation.
DEFAULT_MIN_PER_TASK = 20
DEFAULT_MIN_CLASS = 5
_DATE = re.compile(r"(\d{4})(?:-(\d{2})(?:-(\d{2}))?)?")
_HEX64 = re.compile(r"[0-9a-f]{64}")


class LockboxError(Exception):
    """A lockbox genome reached development, or the lockbox was mishandled.

    A refusal that unseal recorded carries the log's new head as log_head; otherwise log_head is None.
    """

    def __init__(self, message, log_head=None):
        super().__init__(message)
        self.log_head = log_head


def canonical_json(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"))


def sha256_json(obj):
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def latest_day(value):
    """The latest calendar day an ISO date or partial date can mean."""
    match = _DATE.fullmatch(value) if isinstance(value, str) else None
    if match is None:
        raise ValueError("collection_date must be YYYY, YYYY-MM or YYYY-MM-DD, got %r" % (value,))
    year, month, day = match.groups()
    year, month = int(year), int(month or 12)
    if day is None:
        day = calendar.monthrange(year, month)[1]
    return dt.date(year, month, int(day))


def full_date(value):
    """A complete YYYY-MM-DD date, as used for cutoffs."""
    if not isinstance(value, str) or len(value) != 10:
        raise ValueError("expected a YYYY-MM-DD date, got %r" % (value,))
    return latest_day(value)


def identifier(value, field, owner=None, strip=True):
    """Rule 6 v2: the value stripped with str.strip(), refused if empty or if it holds anything
    outside printable ASCII other than an interior space.

    With strip=False (genome ids), surrounding whitespace is refused instead of stripped.
    """
    where = "" if owner is None else " in genome %r" % (owner,)
    if not isinstance(value, str):
        raise ValueError("%s must be a string, got %r%s" % (field, value, where))
    stripped = value.strip()
    if not strip and stripped != value:
        raise ValueError("%s %r%s has surrounding whitespace" % (field, value, where))
    if not stripped:
        raise ValueError("%s must be non-empty%s" % (field, where))
    bad = next((c for c in stripped if not ("!" <= c <= "~" or c == " ")), None)
    if bad is not None:
        raise ValueError("%s %r%s holds U+%04X, outside printable ASCII" % (field, value, where, ord(bad)))
    return stripped


def alias_key(value):
    return unicodedata.normalize("NFKC", value).casefold()


def refuse_aliases(field, *values):
    """ValueError naming both values when two distinct values are equal under NFKC casefolding."""
    seen = {}
    for value in sorted(set().union(*values)):
        first = seen.setdefault(alias_key(value), value)
        if first != value:
            raise ValueError("%s alias collision: %r and %r" % (field, first, value))


def _at_least_one(value, name):
    if type(value) is not int or value < 1:
        raise ValueError("%s must be an integer of at least 1, got %r" % (name, value))
    return value


def _label(genome_id, task, label):
    if type(label) not in (int, bool) or label not in (0, 1):
        raise ValueError("%s: task %r label must be 0 or 1, got %r" % (genome_id, task, label))
    return int(label)


def normalize_genomes(genomes):
    """Validated copies sorted by genome_id; tasks become {task: 0, 1 or None}."""
    rows, classed = [], {}
    for genome in genomes:
        gid = identifier(genome.get("genome_id"), "genome_id", strip=False)
        dup = genome.get("duplicate_group")
        date = genome.get("collection_date")
        if date is not None:
            latest_day(date)
        tasks = genome.get("tasks", ())
        if isinstance(tasks, (str, bytes)):
            raise ValueError("%s: tasks must be a list or mapping, not a string" % gid)
        if isinstance(tasks, dict):
            tasks = {task: None if label is None else _label(gid, task, label) for task, label in tasks.items()}
        else:
            tasks = dict.fromkeys(tasks)
        for task, label in tasks.items():
            if not isinstance(task, str) or not task:
                raise ValueError("%s: task names must be non-empty strings" % gid)
            if classed.setdefault(task, label is not None) != (label is not None):
                raise ValueError("task %r mixes 0/1 labels with unlabelled genomes" % task)
        rows.append({"genome_id": gid, "lineage": identifier(genome.get("lineage"), "lineage", gid),
                     "duplicate_group": None if dup is None else identifier(dup, "duplicate_group", gid),
                     "collection_date": date, "tasks": dict(sorted(tasks.items()))})
    rows.sort(key=lambda row: row["genome_id"])
    repeated = sorted({a["genome_id"] for a, b in zip(rows, rows[1:]) if a["genome_id"] == b["genome_id"]})
    if repeated:
        raise ValueError("duplicate genome_id: %s" % ", ".join(repeated[:5]))
    for field in ("genome_id", "lineage", "duplicate_group"):
        refuse_aliases(field, {row[field] for row in rows} - {None})
    return rows


def task_classes(rows, tasks=()):
    """{task: True for 0/1-labelled tasks}; declared tasks without genomes are False."""
    if isinstance(tasks, (str, bytes)):
        raise ValueError("tasks must be a list of task names")
    classed = dict.fromkeys(tasks, False)
    for row in rows:
        for task, label in row["tasks"].items():
            classed[task] = label is not None
    return dict(sorted(classed.items()))


def task_counts(rows, classed):
    """{task: {"n": genomes carrying the task, "classes": {"0": n, "1": n} or None}} for every task in classed.

    For a list-form task, n counts genomes that list it; P holds no labels for it, so classes is None.
    """
    counts = {task: {"n": 0, "classes": {"0": 0, "1": 0} if has else None} for task, has in classed.items()}
    for row in rows:
        for task, label in row["tasks"].items():
            counts[task]["n"] += 1
            if label is not None:
                counts[task]["classes"][str(label)] += 1
    return counts


def evidence_reasons(count, minimum, min_class, where):
    """Why one task's count falls short; empty when the evidence suffices."""
    reasons = []
    if count["n"] < minimum:
        reasons.append("%s has %d labelled genomes, below %d" % (where, count["n"], minimum))
    for name, n in sorted((count["classes"] or {}).items()):
        if n < min_class:
            reasons.append("%s has %d of class %s, below %d" % (where, n, name, min_class))
    return reasons


def _recent(row, cutoff):
    return row["collection_date"] is not None and latest_day(row["collection_date"]) > cutoff


def _groups_of(pairs):
    groups = {}
    for gid, label in pairs:
        if label is not None:
            groups.setdefault(label, []).append(gid)
    return {label: sorted(ids) for label, ids in sorted(groups.items())}


def define_lockbox(genomes, *, recent_after, access_log, holdout_lineages=(), tasks=(),
                   min_per_task=DEFAULT_MIN_PER_TASK, min_class=DEFAULT_MIN_CLASS):
    """The sealed manifest for recent genomes plus whole held-out lineages and their duplicates.

    access_log is the absolute path write_lockbox starts the log at; unseal only ever uses that log.
    """
    cutoff = full_date(recent_after)
    if not isinstance(access_log, str) or not os.path.isabs(access_log) or access_log != access_log.strip():
        raise ValueError("access_log must be an absolute path, got %r" % (access_log,))
    _at_least_one(min_per_task, "min_per_task")
    _at_least_one(min_class, "min_class")
    if isinstance(holdout_lineages, (str, bytes)):
        raise ValueError("holdout_lineages must be a list of lineage names")
    holdout = sorted({identifier(lineage, "holdout lineage") for lineage in holdout_lineages})
    rows = normalize_genomes(genomes)
    refuse_aliases("lineage", holdout, {row["lineage"] for row in rows})
    recent = {row["genome_id"] for row in rows if _recent(row, cutoff)}
    held = {row["genome_id"] for row in rows if row["lineage"] in holdout}
    chosen = recent | held
    groups = {row["duplicate_group"] for row in rows if row["genome_id"] in chosen} - {None}
    ids = chosen | {row["genome_id"] for row in rows if row["duplicate_group"] in groups}
    classed = task_classes(rows, tasks)
    per_task = task_counts([row for row in rows if row["genome_id"] in ids], classed)
    for task, count in per_task.items():
        count["reasons"] = evidence_reasons(count, min_per_task, min_class, "lockbox")
        count["status"] = "insufficient_evidence" if count["reasons"] else "ok"
    body = {
        "kind": KIND,
        "version": VERSION,
        "definition": {"recent_after": cutoff.isoformat(), "holdout_lineages": holdout, "recent_rule": RECENT_RULE,
                       "duplicate_closure": True, "min_per_task": min_per_task, "min_class": min_class},
        "access_log": access_log,
        "genome_ids": sorted(ids),
        "duplicate_groups": sorted(groups),
        "definition_groups": _groups_of((row["genome_id"], row["duplicate_group"]) for row in rows),
        "counts": {"input_genomes": len(rows), "lockbox_genomes": len(ids), "development_genomes": len(rows) - len(ids),
                   "recent": len(recent), "holdout_lineage": len(held), "duplicate_closure": len(ids - chosen),
                   "undated_input": sum(row["collection_date"] is None for row in rows),
                   "holdout_lineages_absent": sorted(set(holdout) - {row["lineage"] for row in rows})},
        "tasks": per_task,
        "input_sha256": sha256_json(rows),
    }
    return dict(body, seal=sha256_json(body))


def verify_seal(manifest, expected_seal=None):
    """The manifest's seal, after checking it matches the contents and, if given, expected_seal."""
    if not isinstance(manifest, dict) or manifest.get("kind") != KIND:
        raise LockboxError("not a GRM lockbox manifest")
    body = {key: value for key, value in manifest.items() if key != "seal"}
    if manifest.get("seal") != sha256_json(body):
        raise LockboxError("lockbox seal does not match its manifest; refusing a tampered lockbox")
    if expected_seal is not None and manifest["seal"] != expected_seal:
        raise LockboxError("lockbox seal %s is not the expected seal %s; refusing a resealed or different lockbox"
                           % (manifest["seal"], expected_seal))
    return manifest["seal"]


def _now(now):
    return (now or dt.datetime.now(dt.timezone.utc)).isoformat(timespec="seconds")


def _write_record(stream, record):
    stream.write(canonical_json(record) + "\n")
    stream.flush()
    os.fsync(stream.fileno())


def write_lockbox(manifest, path, *, actor, now=None, expected_seal=None):
    """Write the sealed manifest once, read-only, and start its access log with a genesis record.

    Refuses, with FileExistsError, if the manifest file or the access log already exists.
    Returns {"seal", "log_head"}; record both outside the box, since the first unseal needs that head.
    """
    seal = verify_seal(manifest, expected_seal)
    if not isinstance(actor, str) or not actor.strip():
        raise ValueError("actor must be a non-blank string")
    log = manifest["access_log"]
    if os.path.lexists(log):
        raise FileExistsError(errno.EEXIST, "access log already exists", log)
    with open(path, "x", encoding="utf-8") as stream:
        stream.write(canonical_json(manifest) + "\n")
    os.chmod(path, 0o444)
    genesis = {"event": "seal", "time": _now(now), "actor": actor, "purpose": "define lockbox", "approval_ref": None,
               "seal": seal, "outcome": "written", "reason": None, "prev_sha256": LOG_GENESIS}
    genesis["record_sha256"] = sha256_json(genesis)
    with open(log, "x", encoding="utf-8") as stream:
        _write_record(stream, genesis)
    return {"seal": seal, "log_head": {"records": 1, "record_sha256": genesis["record_sha256"]}}


def load_lockbox(path, expected_seal=None):
    with open(path, encoding="utf-8") as stream:
        manifest = json.load(stream)
    verify_seal(manifest, expected_seal)
    return manifest


def _releases(releases):
    """Each earlier release as normalized rows; refuses a single release or record passed by mistake."""
    shape = "releases must be a list of releases, each a list of genome records"
    if isinstance(releases, (dict, str, bytes)):
        raise ValueError(shape)
    releases = list(releases)
    if any(isinstance(release, (dict, str, bytes)) for release in releases):
        raise ValueError(shape)
    return [row for release in releases for row in normalize_genomes(release)]


def lockbox_conflicts(genomes, manifest, expected_seal=None, releases=()):
    """Genomes the lockbox claims, with reasons: lockbox_genome, duplicate_group (a lockbox group),
    holdout_lineage, recent, earlier_release (claimed by one of those rules in a release passed in
    releases), or duplicate_closure (connected through duplicate groups to a lockbox genome or
    another claimed genome).

    releases: every earlier release used for development since definition, each a list of genome
    records; their claims seed the closure and their groups join it. Raises ValueError on an alias
    collision between the input, the releases and the lockbox's stored values.
    """
    verify_seal(manifest, expected_seal)
    ids, groups = set(manifest["genome_ids"]), set(manifest["duplicate_groups"])
    holdout = set(manifest["definition"]["holdout_lineages"])
    cutoff = full_date(manifest["definition"]["recent_after"])
    rows = normalize_genomes(genomes)
    history = _releases(releases)
    stored = [(gid, label) for label, members in manifest["definition_groups"].items() for gid in members]
    refuse_aliases("genome_id", ids, {gid for gid, _ in stored}, {row["genome_id"] for row in rows + history})
    refuse_aliases("lineage", holdout, {row["lineage"] for row in rows + history})
    refuse_aliases("duplicate_group", set(manifest["definition_groups"]),
                   {row["duplicate_group"] for row in rows + history} - {None})

    def hits(row):
        rules = (("lockbox_genome", row["genome_id"] in ids), ("duplicate_group", row["duplicate_group"] in groups),
                 ("holdout_lineage", row["lineage"] in holdout), ("recent", _recent(row, cutoff)))
        return [name for name, hit in rules if hit]

    direct = {row["genome_id"]: hits(row) for row in rows}
    earlier = {row["genome_id"] for row in history if hits(row)}
    parent = {}

    def find(node):
        parent.setdefault(node, node)
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    for gid, label in stored + [(row["genome_id"], row["duplicate_group"]) for row in history + rows]:
        if label is not None:
            parent[find(("genome", gid))] = find(("group", label))
    seeds = {find(("genome", gid)) for gid in ids | earlier | {gid for gid, hit in direct.items() if hit}}
    conflicts = []
    for row in rows:
        gid = row["genome_id"]
        reasons = (direct[gid] or (["earlier_release"] if gid in earlier else [])
                   or (["duplicate_closure"] if find(("genome", gid)) in seeds else []))
        if reasons:
            conflicts.append({"genome_id": gid, "reasons": reasons})
    return conflicts


def require_development(genomes, manifest, expected_seal=None, releases=()):
    """The seal, once no genome belongs to the lockbox; otherwise LockboxError."""
    conflicts = lockbox_conflicts(genomes, manifest, expected_seal, releases)
    if conflicts:
        shown = "; ".join("%s (%s)" % (c["genome_id"], ", ".join(c["reasons"])) for c in conflicts[:5])
        raise LockboxError("%d lockbox genomes cannot enter development train or test: %s" % (len(conflicts), shown))
    return manifest["seal"]


def development_genomes(genomes, manifest, expected_seal=None, releases=()):
    """The input records the lockbox does not claim, in input order. Pass the full release."""
    genomes = list(genomes)
    claimed = {c["genome_id"] for c in lockbox_conflicts(genomes, manifest, expected_seal, releases)}
    return [genome for genome in genomes if genome["genome_id"] not in claimed]


def assert_untouched(manifest, used_ids, expected_seal=None):
    """Raise LockboxError if any used id is a lockbox genome.

    Used ids follow the genome-id rule: ValueError for a non-string, surrounding whitespace or a
    character outside printable ASCII, or an alias of a lockbox id or a stored definition-group member.
    """
    verify_seal(manifest, expected_seal)
    if isinstance(used_ids, (str, bytes)):
        raise TypeError("used_ids must be a collection of genome ids, not one string")
    used = {identifier(gid, "genome_id", strip=False) for gid in used_ids}
    ids = set(manifest["genome_ids"])
    refuse_aliases("genome_id", ids, {gid for members in manifest["definition_groups"].values() for gid in members}, used)
    leaked = sorted(ids & used)
    if leaked:
        raise LockboxError("%d lockbox genomes were used: %s" % (len(leaked), ", ".join(leaked[:10])))


def read_access_log(path):
    """The records in order, after checking the chain from a seal genesis record.

    LockboxError if the log is missing, empty, edited inside the chain or not started by a seal record.
    A cut-off tail or a whole rewritten chain reads clean; see verify_log.
    """
    if not isinstance(path, str) or not os.path.isfile(path):
        raise LockboxError("access log %r is missing" % (path,))
    records, previous = [], LOG_GENESIS
    with open(path, encoding="utf-8") as stream:
        for number, line in enumerate(stream, 1):
            try:
                record = json.loads(line)
            except json.JSONDecodeError:
                record = None
            body = {key: value for key, value in record.items() if key != "record_sha256"} if isinstance(record, dict) else None
            if (not line.endswith("\n") or body is None or record.get("prev_sha256") != previous
                    or record.get("record_sha256") != sha256_json(body)):
                raise LockboxError("access log %s line %d breaks the hash chain" % (path, number))
            previous = record["record_sha256"]
            records.append(record)
    if not records or records[0].get("event") != "seal":
        raise LockboxError("access log %s does not start with a seal record" % path)
    return records


def _valid_head(head):
    return (isinstance(head, dict) and set(head) == {"records", "record_sha256"} and type(head["records"]) is int
            and head["records"] >= 1 and isinstance(head["record_sha256"], str) and _HEX64.fullmatch(head["record_sha256"]))


def verify_log(path, expected_head=None, *, seal):
    """The chain head {"records", "record_sha256"} after read_access_log's checks.

    seal: the genesis record must be bound to it. expected_head: a head recorded earlier, outside
    the box; the log must still hold it at the same position, which detects truncation and rewrites
    of everything up to that head.
    """
    if not isinstance(seal, str) or not _HEX64.fullmatch(seal):
        raise ValueError("seal must be the lockbox's sha256 seal, got %r" % (seal,))
    records = read_access_log(path)
    if records[0].get("seal") != seal:
        raise LockboxError("access log %s was started for another seal" % path)
    if expected_head is not None:
        if not _valid_head(expected_head):
            raise ValueError("expected_head must be {\"records\": n >= 1, \"record_sha256\": sha256}")
        count, digest = expected_head["records"], expected_head["record_sha256"]
        if len(records) < count or records[count - 1]["record_sha256"] != digest:
            raise LockboxError("access log %s no longer holds the recorded head (%d records, %.12s): truncated, rewritten "
                               "or re-created" % (path, count, digest))
    return {"records": len(records), "record_sha256": records[-1]["record_sha256"]}


def _append_access(path, record):
    records = read_access_log(path)
    record = dict(record, prev_sha256=records[-1]["record_sha256"])
    record["record_sha256"] = sha256_json(record)
    with open(path, "a", encoding="utf-8") as stream:
        _write_record(stream, record)
    return {"records": len(records) + 1, "record_sha256": record["record_sha256"]}


def unseal(manifest, purpose, approval_ref, *, actor, expected_seal, expected_head, now=None):
    """Release the lockbox ids for one approved purpose.

    Needs a non-blank actor, purpose and approval_ref; a manifest whose seal verifies and equals
    expected_seal, the seal recorded at seal time; and the manifest's own access log, intact, started
    for that seal, and ending exactly at expected_head, the head recorded after the previous unseal
    attempt (or returned by write_lockbox). If the log or its head fails, nothing is recorded and
    unseal refuses. Every other attempt is recorded, granted or refused; a recorded refusal raises
    LockboxError carrying the new head as log_head.

    Returns {"genome_ids", "log_head"}. The caller must record log_head outside the box after every
    recorded attempt (GRM records it in the work-program ledger); the next unseal needs it.

    Development code exercises this only on synthetic fixtures.
    """
    fields = {"actor": actor, "purpose": purpose, "approval_ref": approval_ref, "expected_seal": expected_seal}
    reason = next(("missing " + name for name, value in fields.items()
                   if not isinstance(value, str) or not value.strip()), None)
    if reason is None:
        try:
            verify_seal(manifest, expected_seal)
        except LockboxError as error:
            reason = str(error)
    seal = manifest.get("seal") if isinstance(manifest, dict) else None
    log = manifest.get("access_log") if isinstance(manifest, dict) else None
    try:
        if not isinstance(seal, str) or not _HEX64.fullmatch(seal):
            raise LockboxError("manifest has no seal to find its access log by")
        if not _valid_head(expected_head):
            raise LockboxError("missing or malformed expected_head")
        head = verify_log(log, expected_head, seal=seal)
        if head != expected_head:
            raise LockboxError("access log %s has %d records past the expected head" % (log, head["records"] - expected_head["records"]))
    except LockboxError as error:
        raise LockboxError("unseal refused, not recorded: %s%s" % (reason + "; " if reason else "", error)) from None
    record = {name: value if isinstance(value, str) else None for name, value in fields.items() if name != "expected_seal"}
    record.update(event="unseal", time=_now(now), seal=seal, outcome="refused" if reason else "granted", reason=reason)
    head = _append_access(log, record)
    if reason:
        raise LockboxError("unseal refused: " + reason, log_head=head)
    return {"genome_ids": list(manifest["genome_ids"]), "log_head": head}
