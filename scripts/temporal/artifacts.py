"""Versioned annotation and catalog artifacts, resolved as of a fit cutoff.

A fit at cutoff T may use only the artifact releases that were public at T:
for each kind, the newest release with released_at <= T. resolve(T) refuses
when a kind has no such release or the chosen release has no sha256, and its
manifest is what a fit at T must declare, exactly.

Release rows, as JSON or mappings:
    {"kind": "catalog", "release_id": "who-v2", "released_at": "2026-03-01",
     "sha256": "<64 lowercase hex, or null until known>",
     "vocabulary": ["katG_S315T", ...], "path": "<optional file to verify>"}

Genes or markers absent from T's vocabulary are encoded by encode_genes under
one of two documented policies and are never added to the vocabulary:
"unknown" maps any number of them to the single UNKNOWN token; "drop" omits
them.
"""
from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

KINDS = ("annotation", "catalog")
UNKNOWN = "<unknown>"
UNSEEN_POLICIES = ("unknown", "drop")
_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")
_SHA256 = re.compile(r"[0-9a-f]{64}")
_RELEASE_KEYS = frozenset({"kind", "release_id", "released_at", "sha256", "vocabulary", "path"})
_MANIFEST_KEYS = frozenset({"cutoff", "artifacts"})
_ENTRY_KEYS = frozenset({"release_id", "released_at", "sha256"})


class ArtifactError(ValueError):
    """A registry, resolution or manifest that a fit must not use."""


def parse_date(value, what="date"):
    """Return a datetime.date from a date or a strict YYYY-MM-DD string."""
    if isinstance(value, dt.datetime):
        raise ArtifactError(f"{what}: expected a date, not a datetime: {value!r}")
    if isinstance(value, dt.date):
        return value
    if isinstance(value, str) and _DATE.fullmatch(value):
        try:
            return dt.date.fromisoformat(value)
        except ValueError as error:
            raise ArtifactError(f"{what}: {error}") from None
    raise ArtifactError(f"{what}: expected an ISO date YYYY-MM-DD, got {value!r}")


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _text(value, what):
    if not isinstance(value, str) or not value:
        raise ArtifactError(f"{what}: expected a non-empty string, got {value!r}")
    return value


@dataclass(frozen=True)
class Release:
    """One release of one artifact kind. sha256 may be None until it is known;
    resolve() refuses to hand such a release to a fit."""
    kind: str
    release_id: str
    released_at: dt.date
    sha256: str | None
    vocabulary: frozenset = frozenset()
    path: str | None = None

    def __post_init__(self):
        where = f"release {self.kind!r}:{self.release_id!r}"
        _text(self.kind, "release kind")
        _text(self.release_id, f"{where} release_id")
        object.__setattr__(self, "released_at", parse_date(self.released_at, f"{where} released_at"))
        if self.sha256 == "":
            object.__setattr__(self, "sha256", None)
        if self.sha256 is not None and not (isinstance(self.sha256, str) and _SHA256.fullmatch(self.sha256)):
            raise ArtifactError(f"{where}: sha256 must be 64 lowercase hex characters, got {self.sha256!r}")
        if isinstance(self.vocabulary, str):
            raise ArtifactError(f"{where}: vocabulary must be a collection of strings, not a string")
        vocabulary = frozenset(self.vocabulary)
        for term in vocabulary:
            _text(term, f"{where} vocabulary term")
        if UNKNOWN in vocabulary:
            raise ArtifactError(f"{where}: vocabulary uses the reserved token {UNKNOWN!r}")
        object.__setattr__(self, "vocabulary", vocabulary)
        if self.path is not None:
            object.__setattr__(self, "path", _text(str(self.path), f"{where} path"))

    @classmethod
    def from_mapping(cls, row):
        if not isinstance(row, Mapping):
            raise ArtifactError(f"release row: expected a mapping, got {type(row).__name__}")
        unknown = set(row) - _RELEASE_KEYS
        missing = {"kind", "release_id", "released_at", "sha256"} - set(row)
        if unknown or missing:
            raise ArtifactError(f"release row {row.get('release_id')!r}: "
                                f"unknown keys {sorted(unknown)}, missing keys {sorted(missing)}")
        return cls(row["kind"], row["release_id"], row["released_at"], row["sha256"],
                   row.get("vocabulary", ()), row.get("path"))

    @property
    def ref(self):
        """kind:release_id@released_at, the form canaries report."""
        return f"{self.kind}:{self.release_id}@{self.released_at.isoformat()}"

    def entry(self):
        return {"release_id": self.release_id, "released_at": self.released_at.isoformat(),
                "sha256": self.sha256}


@dataclass(frozen=True)
class Resolution:
    """The releases a fit at cutoff may use: one per kind, sorted by kind."""
    cutoff: dt.date
    releases: tuple

    def release(self, kind):
        for release in self.releases:
            if release.kind == kind:
                return release
        raise ArtifactError(f"no {kind!r} release in this resolution")

    def vocabulary(self, kind=None):
        if kind is not None:
            return self.release(kind).vocabulary
        return frozenset().union(*(release.vocabulary for release in self.releases))

    def manifest(self):
        return {"cutoff": self.cutoff.isoformat(),
                "artifacts": {release.kind: release.entry() for release in self.releases}}


class ArtifactRegistry:
    """Every release of every kind, with the as-of lookups a fit at T may use."""

    def __init__(self, releases=(), kinds=KINDS):
        kinds = tuple(sorted(set(kinds)))
        if not kinds:
            raise ArtifactError("a registry needs at least one artifact kind")
        for kind in kinds:
            _text(kind, "artifact kind")
        self.kinds = kinds
        self._releases = {kind: [] for kind in kinds}
        for release in releases:
            self.add(release)

    @classmethod
    def from_records(cls, rows, kinds=KINDS):
        return cls((Release.from_mapping(row) for row in rows), kinds)

    @classmethod
    def load(cls, path, kinds=KINDS):
        """Load a JSON list of release rows."""
        rows = json.loads(Path(path).read_text())
        if not isinstance(rows, list):
            raise ArtifactError(f"{path}: expected a JSON list of release rows")
        return cls.from_records(rows, kinds)

    def add(self, release):
        if isinstance(release, Mapping):
            release = Release.from_mapping(release)
        if not isinstance(release, Release):
            raise ArtifactError(f"expected a Release or a mapping, got {type(release).__name__}")
        if release.kind not in self._releases:
            raise ArtifactError(f"release {release.release_id!r}: unknown kind {release.kind!r}, "
                                f"expected one of {list(self.kinds)}")
        existing = self._releases[release.kind]
        for other in existing:
            if other.release_id == release.release_id:
                raise ArtifactError(f"duplicate {release.kind} release_id {release.release_id!r}")
            if other.released_at == release.released_at:
                raise ArtifactError(f"{release.kind} releases {other.release_id!r} and {release.release_id!r} "
                                    f"share released_at {release.released_at}; as-of resolution would be ambiguous")
        existing.append(release)
        existing.sort(key=lambda r: r.released_at)
        return release

    def releases(self, kind=None):
        """All releases, oldest first within each kind, kinds in sorted order."""
        if kind is not None:
            if kind not in self._releases:
                raise ArtifactError(f"unknown kind {kind!r}")
            return tuple(self._releases[kind])
        return tuple(release for k in self.kinds for release in self._releases[k])

    def resolve(self, cutoff, verify_files=False):
        """The newest release of each kind with released_at <= cutoff.

        Refuses when a kind has no such release or the chosen release has no
        sha256; it never falls back to an older release. With verify_files,
        each chosen release must name a path whose sha256 matches.
        """
        cutoff = parse_date(cutoff, "cutoff")
        chosen, problems = [], []
        for kind in self.kinds:
            candidates = [r for r in self._releases[kind] if r.released_at <= cutoff]
            if not candidates:
                problems.append(f"{kind}: no release on or before {cutoff}")
                continue
            release = candidates[-1]
            if release.sha256 is None:
                problems.append(f"{release.ref}: sha256 is missing")
            elif verify_files:
                if release.path is None:
                    problems.append(f"{release.ref}: no path to verify")
                elif not Path(release.path).is_file():
                    problems.append(f"{release.ref}: {release.path} is not a file")
                elif sha256_file(release.path) != release.sha256:
                    problems.append(f"{release.ref}: {release.path} does not match its sha256")
            chosen.append(release)
        if problems:
            raise ArtifactError(f"cannot resolve artifacts as of {cutoff}: " + "; ".join(problems))
        return Resolution(cutoff, tuple(chosen))

    def vocabulary(self, cutoff, kind=None):
        """The marker or gene set as of cutoff: one kind's, or the union of all kinds."""
        return self.resolve(cutoff).vocabulary(kind)

    def manifest(self, cutoff):
        return self.resolve(cutoff).manifest()

    def check_declared(self, cutoff, declared):
        """Raise unless declared is exactly the manifest for cutoff."""
        mismatches = manifest_mismatches(self.manifest(cutoff), declared)
        if mismatches:
            raise ArtifactError("declared artifacts do not match the cutoff's resolution: "
                                + "; ".join(f"{where}: {problem}" for where, problem in mismatches))


def manifest_mismatches(expected, declared):
    """Sorted (where, problem) pairs; empty when declared is exactly expected."""
    if not isinstance(declared, Mapping):
        return [("manifest", f"expected a mapping, got {type(declared).__name__}")]
    problems = []
    keys = set(declared)
    if keys != _MANIFEST_KEYS:
        problems.append(("manifest", f"keys {sorted(keys)}, expected {sorted(_MANIFEST_KEYS)}"))
    if declared.get("cutoff") != expected["cutoff"]:
        problems.append(("cutoff", f"declared {declared.get('cutoff')!r}, resolved {expected['cutoff']!r}"))
    entries = declared.get("artifacts")
    if not isinstance(entries, Mapping):
        problems.append(("artifacts", f"expected a mapping, got {type(entries).__name__}"))
        return sorted(problems)
    for kind in sorted(set(entries) - set(expected["artifacts"]), key=str):
        problems.append((str(kind), "declared but not a resolved kind"))
    for kind, want in sorted(expected["artifacts"].items()):
        have = entries.get(kind)
        if have is None:
            problems.append((kind, f"not declared; resolved {want['release_id']}"))
        elif not isinstance(have, Mapping) or set(have) != _ENTRY_KEYS:
            problems.append((kind, f"entry must have exactly keys {sorted(_ENTRY_KEYS)}"))
        else:
            for key in sorted(_ENTRY_KEYS):
                if have[key] != want[key]:
                    problems.append((kind, f"{key} declared {have[key]!r}, resolved {want[key]!r}"))
    return sorted(problems)


def encode_genes(genes, vocabulary, unseen="unknown"):
    """Encode one genome's genes against a fixed vocabulary.

    Returns the sorted, de-duplicated genes found in the vocabulary. Under
    unseen="unknown", any unseen genes add one trailing UNKNOWN token; under
    "drop" they are omitted. The vocabulary is only read, never extended.
    """
    if unseen not in UNSEEN_POLICIES:
        raise ArtifactError(f"unseen policy must be one of {list(UNSEEN_POLICIES)}, got {unseen!r}")
    if isinstance(genes, str):
        raise ArtifactError("genes must be a collection of strings, not a string")
    genes = set(genes)
    known = sorted(gene for gene in genes if gene in vocabulary)
    if unseen == "unknown" and len(known) < len(genes):
        known.append(UNKNOWN)
    return tuple(known)
