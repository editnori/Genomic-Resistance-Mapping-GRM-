"""Fit-isolation canaries for temporal later-label benchmarking.

Each canary checks one way a fit at cutoff T could see data from after T and
returns a CanaryResult: passed, plus the offending genome ids, features or
artifacts. Every canary is meant to run on a clean synthetic control, where it
must pass, and on a planted leak, where it must fail.

This module carries its own loader and validator for the shared record
contract, so it never imports splits.py.

Shared record contract (Task M and Task N document this same text):

- **`GenomeRecord`:** `genome_id` (str), `collection_date` (ISO date or null), `released_at` (ISO date the assembly became public, or null), `features_available_at` (ISO date or null; defaults to `released_at`), `lineage` (str or null), `duplicate_group` (str or null), and `labels` (a list of `LabelEvent`).
- **`LabelEvent`:** `antibiotic`, `value` (R or S), `recorded_at` (ISO date the label became known), `source_release` (str), and `supersedes` (the `recorded_at` of the event it revises, or null). Together these form a bitemporal label history.
- **`LabelAsOf(record, antibiotic, T)`:** the latest event with `recorded_at <= T`, or none.

Binding rulings on the contract (main, 2026-10-30):

1. `recorded_at` null is accepted at load. The event is never knowable, and a record whose label relies on it is excluded with reason `unknown_label_recorded_at`, counted.
2. `labels` is required, and may be an empty list. `source_release` is required and non-empty.
3. Same-day events: agreeing events on one day collapse to one: the first in canonical label order, that is the smallest `source_release`, then an event with null `supersedes` before one with a date, then the earliest `supersedes`. The pick does not depend on input order. Conflicting events on one day make the label ambiguous only within the as-of window that reads them: <= T for the train label, <= `eval_date` for the test label.
4. `supersedes` must name an existing earlier event of the same antibiotic, or the record is refused at load.
5. Availability = max(`released_at`, `features_available_at`) when both are present. Features can't predate public release. Count such records as `inconsistent_dates`. A label recorded before `collection_date` is refused at load.
6. For `lineage` and `duplicate_group`, a non-null value is first stripped of surrounding whitespace with Python `str.strip()`. The stripped value is refused if it's empty, or if it contains any character outside printable ASCII (U+0021 to U+007E) other than an interior space (U+0020). Two distinct values of one field in one input that are equal under `unicodedata.normalize("NFKC", v).casefold()` are refused as an alias collision, naming both. Values are otherwise kept as given. All refusals are `ValueError` naming the genome or the two values.

How this module applies the rulings:

- Rule 1: label_as_of ignores undated events at every date, so an undated
  event never makes a genome or label knowable. The backfill and revision
  canaries name a training label taken from one.
- Rule 2: a record needs `genome_id` and `labels`; an event needs
  `antibiotic`, `value`, `recorded_at` (which may be null) and
  `source_release`. A missing optional key is null; an unknown key is refused.
- Rule 3: records keep every event. label_as_of reads the latest dated day on
  or before its date. Agreeing events there return the first in canonical
  label order: the smallest `source_release`, then `supersedes` null first,
  then the earliest `supersedes` (ruling B, 2026-10-30), so the pick never
  depends on input order. Conflicting events there raise AmbiguousLabel (a
  ValueError). A conflict on an earlier day that a later dated event in the
  window superseded is not ambiguous. Train-side reads use T only, so a
  conflict recorded after T never touches training.
- Rule 4: `supersedes` must name a dated event of the same antibiotic on an
  earlier day than its own (when its own is dated).
- Rule 5: `available_at` is the later of `released_at` and
  `features_available_at` that is present; `inconsistent_dates` flags
  features dated before release.
- Rule 6: GenomeRecord strips each value with str.strip() first, then
  refuses an empty value or any character outside U+0020 to U+007E (after
  stripping, every space left is interior); load_records, which
  every canary uses, checks each field for alias collisions across its
  input. Refusals are RecordError (a ValueError) naming the genome or both
  values.

A genome is knowable at T when `available_at` is on or before T; one with
neither date never is. A training row at T (ruling A, 2026-10-30) is a
genome knowable at T whose label for a fitted antibiotic is knowable and
unambiguous at T: the genomes in Task M's train set at T. A fit at T may use
only training rows: their record fields, their features and their dated label
events recorded on or before T. Genomes that were public at T but had no
usable label at T, which includes the later-label test genomes, must not
influence the fit at all, including unsupervised steps such as vocabulary
building, screening or scaling. This is the fold-local rule GRM's pipeline
uses. Everything else is test-period data: every genome outside the training
rows (never knowable, not yet knowable, ambiguous at T, unlabelled at T, or
later), and on every genome, training rows included, the undated events and
anything recorded after T, including same-day conflicts. The invariance
canary covers it in four ways. training_view_only hands the fit exactly the
training view: the training rows with only their events dated on or before T,
and only their features. A fit whose output on the input it ran on differs
from its output on that input's training view fails, whatever it reads,
including whether any test-period data is present at all, provided the
difference shows in what the fingerprint reads: the allowlisted types and
their value forms, and dtype forms with their alias character (see
invariance_to_test_data; any other type is refused). The canary guards
against honest fits that accidentally read test-period data. It doesn't
defend against a fit built to hide information on purpose in private
(underscore) attributes, object internals, memory layout or padding bytes.
drop_non_training removes every genome outside the training rows, and
drop_post_cutoff_events
removes every undated event and every event recorded after T from every
genome; each leaves the other kind in place, so together they say which kind
a fit reads (value, date, source or supersedes of a later event, or any field
of another genome). The other perturbations change one
kind of test-period data at a time (later and undated labels with their
dates, sources and supersedes, same-day conflicts, genes, lineage, group and
dates, added genomes), which catches a fit whose dependence cancels out when
that data simply disappears. training_view_only also removes every
antibiotic that only test-period records carry, so a fit that names a label
column for each antibiotic in its input is caught. None of them covers a fit
that caches across calls, or a dependence that shows only at other input sizes
than the one the canary ran on (see invariance_to_test_data).
"""
from __future__ import annotations

import dataclasses
import datetime as dt
import hashlib
import importlib.util
import json
import math
import struct
import sys
import unicodedata
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

try:
    import numpy as np
except ImportError:  # the canaries themselves do not need numpy
    np = None

if __package__:
    from . import artifacts
else:
    # Loaded by path: share one module object, so ArtifactError and Release
    # are the same classes for every caller that loads artifacts.py this way.
    artifacts = sys.modules.get("grm_temporal_artifacts")
    if artifacts is None:
        _spec = importlib.util.spec_from_file_location(
            "grm_temporal_artifacts", Path(__file__).with_name("artifacts.py"))
        artifacts = importlib.util.module_from_spec(_spec)
        sys.modules[_spec.name] = artifacts
        _spec.loader.exec_module(artifacts)

ArtifactError = artifacts.ArtifactError
UNKNOWN = artifacts.UNKNOWN
VALUES = ("R", "S")
DAY = dt.timedelta(days=1)
GROUP_FIELDS = ("lineage", "duplicate_group")
_RECORD_KEYS = frozenset({"genome_id", "collection_date", "released_at", "features_available_at",
                          "lineage", "duplicate_group", "labels"})
_EVENT_KEYS = frozenset({"antibiotic", "value", "recorded_at", "source_release", "supersedes"})
_ROW_KEYS = frozenset({"genome_id", "antibiotic", "value", "recorded_at"})


class RecordError(ValueError):
    """A record, label event or training row outside the shared contract."""


def _date(value, what, nullable=True):
    if value is None and nullable:
        return None
    try:
        return artifacts.parse_date(value, what)
    except ArtifactError as error:
        raise RecordError(str(error)) from None


def _text(value, what, nullable=False):
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not value:
        raise RecordError(f"{what}: expected a non-empty string, got {value!r}")
    return value


class AmbiguousLabel(ValueError):
    """The day an as-of read uses holds events with different values (rule 3)."""


def _iso(value):
    return None if value is None else value.isoformat()


def _undated_first(value):
    return (value is not None, value or dt.date.min)


@dataclass(frozen=True)
class LabelEvent:
    antibiotic: str
    value: str
    recorded_at: dt.date | None
    source_release: str
    supersedes: dt.date | None = None

    def __post_init__(self):
        _text(self.antibiotic, "label antibiotic")
        where = f"{self.antibiotic} label"
        if self.value not in VALUES:
            raise RecordError(f"{where}: value must be R or S, got {self.value!r}")
        object.__setattr__(self, "recorded_at", _date(self.recorded_at, f"{where} recorded_at"))
        _text(self.source_release, f"{where} source_release")
        object.__setattr__(self, "supersedes", _date(self.supersedes, f"{where} supersedes"))
        if self.recorded_at is not None and self.supersedes is not None and self.supersedes >= self.recorded_at:
            raise RecordError(f"{where} recorded {self.recorded_at}: supersedes {self.supersedes}, "
                              "which is not earlier")

    def sort_key(self):
        return (self.antibiotic, _undated_first(self.recorded_at), self.value, self.source_release,
                _undated_first(self.supersedes))

    def to_mapping(self):
        return {"antibiotic": self.antibiotic, "value": self.value, "recorded_at": _iso(self.recorded_at),
                "source_release": self.source_release, "supersedes": _iso(self.supersedes)}


@dataclass(frozen=True)
class GenomeRecord:
    genome_id: str
    collection_date: dt.date | None = None
    released_at: dt.date | None = None
    features_available_at: dt.date | None = None
    lineage: str | None = None
    duplicate_group: str | None = None
    labels: tuple = ()

    def __post_init__(self):
        gid = _text(self.genome_id, "genome_id")
        for name in ("collection_date", "released_at", "features_available_at"):
            object.__setattr__(self, name, _date(getattr(self, name), f"{gid} {name}"))
        for name in GROUP_FIELDS:  # rule 6
            value = getattr(self, name)
            if isinstance(value, str):
                value = value.strip()
                object.__setattr__(self, name, value)
            _text(value, f"{gid} {name}", True)
            if value is not None and not all(" " <= c <= "~" for c in value):
                raise RecordError(f"{gid} {name}: {value!r} contains a character outside printable ASCII")
        if not isinstance(self.labels, (list, tuple)):
            raise RecordError(f"{gid}: labels must be a list of label events")
        for event in self.labels:
            if not isinstance(event, LabelEvent):
                raise RecordError(f"{gid}: labels must be LabelEvents, got {type(event).__name__}")
        labels = tuple(sorted(self.labels, key=LabelEvent.sort_key))
        dated = {(e.antibiotic, e.recorded_at) for e in labels if e.recorded_at is not None}
        for event in labels:
            if (event.recorded_at is not None and self.collection_date is not None
                    and event.recorded_at < self.collection_date):
                raise RecordError(f"{gid}: {event.antibiotic} label recorded {event.recorded_at} "
                                  f"before collection {self.collection_date}")
            if event.supersedes is not None and (event.antibiotic, event.supersedes) not in dated:
                raise RecordError(f"{gid}: {event.antibiotic} event recorded {event.recorded_at} supersedes "
                                  f"{event.supersedes}, which has no dated {event.antibiotic} event")
        object.__setattr__(self, "labels", labels)

    @property
    def available_at(self):
        """When this genome's features became usable: the later of released_at
        and features_available_at that is present (rule 5)."""
        return max((d for d in (self.released_at, self.features_available_at) if d is not None), default=None)

    @property
    def inconsistent_dates(self):
        """Features dated before public release (rule 5): a flag, not a refusal."""
        return (self.released_at is not None and self.features_available_at is not None
                and self.features_available_at < self.released_at)

    def to_mapping(self):
        """The contract mapping: ISO dates, features_available_at not defaulted."""
        return {"genome_id": self.genome_id, "collection_date": _iso(self.collection_date),
                "released_at": _iso(self.released_at), "features_available_at": _iso(self.features_available_at),
                "lineage": self.lineage, "duplicate_group": self.duplicate_group,
                "labels": [event.to_mapping() for event in self.labels]}


def load_event(row):
    if isinstance(row, LabelEvent):
        return row
    if not isinstance(row, Mapping):
        raise RecordError(f"label event: expected a mapping, got {type(row).__name__}")
    unknown = set(row) - _EVENT_KEYS
    missing = {"antibiotic", "value", "recorded_at", "source_release"} - set(row)
    if unknown or missing:
        raise RecordError(f"label event: unknown keys {sorted(unknown)}, missing keys {sorted(missing)}")
    return LabelEvent(row["antibiotic"], row["value"], row["recorded_at"], row["source_release"],
                      row.get("supersedes"))


def load_record(row):
    """A GenomeRecord from a contract mapping. genome_id and labels are required;
    absent nullable fields are null; unknown keys are refused."""
    if isinstance(row, GenomeRecord):
        return row
    if not isinstance(row, Mapping):
        raise RecordError(f"genome record: expected a mapping, got {type(row).__name__}")
    unknown = set(row) - _RECORD_KEYS
    missing = {"genome_id", "labels"} - set(row)
    if unknown or missing:
        raise RecordError(f"genome record {row.get('genome_id')!r}: unknown keys {sorted(unknown)}, "
                          f"missing keys {sorted(missing)}")
    labels = row["labels"]
    if not isinstance(labels, (list, tuple)):
        raise RecordError(f"{row['genome_id']!r}: labels must be a list")
    return GenomeRecord(row["genome_id"], row.get("collection_date"), row.get("released_at"),
                        row.get("features_available_at"), row.get("lineage"), row.get("duplicate_group"),
                        tuple(load_event(event) for event in labels))


def alias_key(value):
    """Rule 6: two distinct values with the same key are an alias collision."""
    return unicodedata.normalize("NFKC", value).casefold()


def load_records(rows):
    """Records sorted by genome_id. Duplicate ids and rule-6 alias collisions
    (per field, across the whole input) are refused."""
    records = sorted((load_record(row) for row in rows), key=lambda r: r.genome_id)
    for before, after in zip(records, records[1:]):
        if before.genome_id == after.genome_id:
            raise RecordError(f"duplicate genome_id {after.genome_id!r}")
    for name in GROUP_FIELDS:
        seen = {}
        for record in records:
            value = getattr(record, name)
            if value is None:
                continue
            first = seen.setdefault(alias_key(value), (value, record.genome_id))
            if first[0] != value:
                raise RecordError(f"{name} alias collision: {first[0]!r} ({first[1]}) and {value!r} "
                                  f"({record.genome_id}) are equal under NFKC casefold")
    return tuple(records)


def known_events(record, antibiotic, cutoff):
    """The dated events for antibiotic recorded on or before cutoff, oldest first."""
    cutoff = _date(cutoff, "cutoff", False)
    return [e for e in record.labels
            if e.antibiotic == antibiotic and e.recorded_at is not None and e.recorded_at <= cutoff]


def label_as_of(record, antibiotic, cutoff):
    """LabelAsOf: the latest dated event for antibiotic recorded on or before
    cutoff, or None. Agreeing events on that day return the first in canonical
    label order (smallest source_release, then supersedes null first, then
    earliest); conflicting ones raise AmbiguousLabel."""
    known = known_events(record, antibiotic, cutoff)
    if not known:
        return None
    day = known[-1].recorded_at
    tied = [e for e in known if e.recorded_at == day]
    if len({e.value for e in tied}) > 1:
        raise AmbiguousLabel(f"{record.genome_id}: {antibiotic} has conflicting labels on {day}")
    return min(tied, key=LabelEvent.sort_key)


def knowable(record, cutoff):
    """Whether a fit at cutoff may use this genome at all."""
    cutoff = _date(cutoff, "cutoff", False)
    return record.available_at is not None and record.available_at <= cutoff


def inconsistent_dates(records):
    """Sorted ids of records whose features predate their release (rule 5)."""
    return [r.genome_id for r in load_records(records) if r.inconsistent_dates]


def _antibiotics(records):
    return sorted({e.antibiotic for r in records for e in r.labels}) or ["canary-antibiotic"]


def is_training_row(record, cutoff, antibiotics):
    """Ruling A: knowable at cutoff with an unambiguous label as of cutoff for
    at least one of the antibiotics."""
    if not knowable(record, cutoff):
        return False
    for antibiotic in antibiotics:
        try:
            if label_as_of(record, antibiotic, cutoff) is not None:
                return True
        except AmbiguousLabel:
            pass
    return False


def training_ids(records, cutoff, antibiotics=None):
    """The genome ids a fit at cutoff may use (ruling A). antibiotics=None
    means every antibiotic in the records."""
    records = load_records(records)
    antibiotics = _antibiotics(records) if antibiotics is None else sorted(antibiotics)
    return frozenset(r.genome_id for r in records if is_training_row(r, cutoff, antibiotics))


def training_labels_as_of(records, antibiotic, cutoff):
    """Reference training rows: each knowable genome's label as of cutoff.

    Reads only dated events on or before cutoff. A genome whose label as of
    cutoff is ambiguous, or that has none, gets no row. No embargo is applied.
    """
    rows = []
    for record in load_records(records):
        try:
            event = label_as_of(record, antibiotic, cutoff) if knowable(record, cutoff) else None
        except AmbiguousLabel:
            continue
        if event is not None:
            rows.append({"genome_id": record.genome_id, "antibiotic": antibiotic,
                         "value": event.value, "recorded_at": event.recorded_at.isoformat()})
    return rows


@dataclass(frozen=True)
class CanaryResult:
    canary: str
    passed: bool
    offenders: tuple = ()
    detail: str = ""
    constant_keys: tuple = ()  # the caller's claims (TL-4); only the invariance canary takes any

    def to_dict(self):
        return {"canary": self.canary, "passed": self.passed,
                "offenders": [list(offender) for offender in self.offenders], "detail": self.detail,
                "constant_keys": list(self.constant_keys)}


def _result(canary, offenders, clean, leak, constant_keys=()):
    offenders = tuple(sorted(set(offenders)))
    declared = tuple(sorted(constant_keys))
    claims = f"; constant_keys declared by the caller, each a claim to check: {list(declared)}" if declared else ""
    if offenders:
        return CanaryResult(canary, False, offenders, f"{len(offenders)} offender(s): {leak}{claims}", declared)
    return CanaryResult(canary, True, (), clean + claims, declared)


# --- artifacts -------------------------------------------------------------

def future_marker(registry, cutoff, features):
    """Every feature a fit at cutoff uses must be in the cutoff's vocabulary.

    Offenders are (feature, artifact, reason): the first later release that
    carries the feature, the release that retired it, or "-" if no
    registered release has it. UNKNOWN is the unseen-gene token, not a feature.

    Known limit (ruling C, 2026-10-30): this checks names only. A future
    marker renamed onto a term already in the vocabulary at cutoff is not
    caught here; invariance_to_test_data is the backstop for that.
    """
    cutoff = _date(cutoff, "cutoff", False)
    if isinstance(features, str):
        raise RecordError("features must be a collection of strings, not a string")
    features = set(features)
    try:
        vocabulary = registry.vocabulary(cutoff)
    except ArtifactError as error:
        return _result("future_marker", [("-", "-", str(error))], "", "artifacts do not resolve")
    offenders = []
    for feature in sorted(features - vocabulary - {UNKNOWN}):
        holders = [r for r in registry.releases() if feature in r.vocabulary]
        later = sorted((r for r in holders if r.released_at > cutoff), key=lambda r: (r.released_at, r.kind))
        if later:
            offenders.append((feature, later[0].ref, "only in a release after the cutoff"))
        elif holders:
            last = max(holders, key=lambda r: (r.released_at, r.kind))
            offenders.append((feature, last.ref, "retired before the cutoff"))
        else:
            offenders.append((feature, "-", "in no registered release"))
    return _result("future_marker", offenders,
                   f"all {len(features)} feature(s) are in the vocabulary as of {cutoff}",
                   f"features outside the vocabulary as of {cutoff}")


def artifact_manifest(registry, cutoff, declared):
    """The fit's declared artifacts must be exactly the cutoff's resolution.
    Offenders are (kind or field, problem)."""
    cutoff = _date(cutoff, "cutoff", False)
    try:
        expected = registry.manifest(cutoff)
    except ArtifactError as error:
        return _result("artifact_manifest", [("-", str(error))], "", "artifacts do not resolve")
    return _result("artifact_manifest", artifacts.manifest_mismatches(expected, declared),
                   f"declared artifacts match the resolution as of {cutoff}",
                   f"declared artifacts differ from the resolution as of {cutoff}")


# --- labels ----------------------------------------------------------------

def _training_rows(rows):
    parsed = []
    for row in rows:
        if isinstance(row, Mapping):
            unknown = set(row) - _ROW_KEYS
            if unknown or not {"genome_id", "antibiotic", "value"} <= set(row):
                raise RecordError(f"training row {dict(row)!r}: needs genome_id, antibiotic, value "
                                  "and optionally recorded_at")
            gid, antibiotic, value, recorded_at = (row["genome_id"], row["antibiotic"], row["value"],
                                                   row.get("recorded_at"))
        else:
            gid, antibiotic, value, recorded_at = (tuple(row) + (None,))[:4]
        _text(gid, "training row genome_id")
        _text(antibiotic, f"{gid} training row antibiotic")
        if value not in VALUES:
            raise RecordError(f"{gid} training row: value must be R or S, got {value!r}")
        parsed.append((gid, antibiotic, value, _date(recorded_at, f"{gid} training row recorded_at")))
    return parsed


def _unknowable(record):
    return f"features available {record.available_at or 'never'}, after the cutoff"


def backfill(records, cutoff, training_labels):
    """No training label may be one that was unknown at cutoff.

    A training row whose genome had no dated label for the antibiotic on or
    before cutoff is a leak: typically a genome collected before cutoff whose
    label was recorded after it, or one whose only label is undated (rule 1).
    An ambiguous label as of cutoff did exist; label_revision reports it.
    A row for a genome not knowable at cutoff is a leak whatever its label
    (both label canaries report it). Offenders are (genome_id, antibiotic,
    reason).
    """
    cutoff = _date(cutoff, "cutoff", False)
    index = {r.genome_id: r for r in load_records(records)}
    offenders = []
    for gid, antibiotic, _value, _recorded_at in _training_rows(training_labels):
        record = index.get(gid)
        if record is None:
            offenders.append((gid, antibiotic, "genome not in the records"))
            continue
        if not knowable(record, cutoff):
            offenders.append((gid, antibiotic, _unknowable(record)))
            continue
        if known_events(record, antibiotic, cutoff):
            continue
        events = [e for e in record.labels if e.antibiotic == antibiotic]
        dated = [e for e in events if e.recorded_at is not None]
        collected = record.collection_date.isoformat() if record.collection_date else "unknown"
        if dated:
            offenders.append((gid, antibiotic, f"first recorded {dated[0].recorded_at}, after the cutoff; "
                                               f"collected {collected}"))
        elif events:
            offenders.append((gid, antibiotic, "only undated label events, never knowable"))
        else:
            offenders.append((gid, antibiotic, "no label event at all"))
    return _result("backfill", offenders, f"every training label was recorded on or before {cutoff}",
                   f"training labels not recorded by {cutoff}")


def label_revision(records, cutoff, training_labels):
    """Every training label must equal the label as of cutoff.

    A row that disagrees, or that declares an event recorded after cutoff,
    used a later revision or an undated event. A row whose label as of cutoff
    is ambiguous (rule 3) has no single value to train on, so it is an
    offender too, as is a row for a genome not knowable at cutoff. Rows with
    no label as of cutoff belong to the backfill canary. Offenders are
    (genome_id, antibiotic, reason).
    """
    cutoff = _date(cutoff, "cutoff", False)
    index = {r.genome_id: r for r in load_records(records)}
    offenders = []
    for gid, antibiotic, value, recorded_at in _training_rows(training_labels):
        record = index.get(gid)
        if record is None:
            offenders.append((gid, antibiotic, "genome not in the records"))
            continue
        if not knowable(record, cutoff):
            offenders.append((gid, antibiotic, _unknowable(record)))
            continue
        try:
            known = label_as_of(record, antibiotic, cutoff)
        except AmbiguousLabel:
            day = known_events(record, antibiotic, cutoff)[-1].recorded_at
            offenders.append((gid, antibiotic, f"R and S both recorded {day}, the latest day by the cutoff"))
            continue
        if known is None:
            continue
        as_of = f"as of the cutoff {known.value}, recorded {known.recorded_at}"
        if value != known.value:
            events = [e for e in record.labels if e.antibiotic == antibiotic and e.value == value]
            revisions = [e for e in events if e.recorded_at is not None and e.recorded_at > cutoff]
            if revisions:
                offenders.append((gid, antibiotic, f"{value} from the revision recorded "
                                                   f"{revisions[0].recorded_at}; {as_of}"))
            elif events:
                offenders.append((gid, antibiotic, f"{value} from an undated event, never knowable; {as_of}"))
            else:
                offenders.append((gid, antibiotic, f"{value} is in no later revision; {as_of}"))
        elif recorded_at is not None and recorded_at > cutoff:
            offenders.append((gid, antibiotic, f"declares the event recorded {recorded_at}; {as_of}"))
    return _result("label_revision", offenders, f"every training label equals its label as of {cutoff}",
                   f"training labels revised after {cutoff}")


# --- test-data invariance --------------------------------------------------

FINGERPRINT_TYPES = ("None, bool, int, float, str, bytes, dict, list, tuple, set, frozenset, datetime.date, "
                     "numpy.ndarray (structured too), numpy scalars (a void scalar only if structured)")


def _sorted_canonical(items):
    return sorted(items, key=lambda v: json.dumps(v, sort_keys=True))


def _refuse(value, why=""):
    raise TypeError(f"cannot fingerprint a fitted {type(value).__module__}.{type(value).__qualname__}{why}; "
                    f"a fit may return only {FINGERPRINT_TYPES}, by exact type (MC2-1)")


def _canonical(value):
    """A JSON-ready form of a fitted value (MC2-1, MC3-1, MC4). It fails
    closed: a value is read only if its type is exactly one of
    FINGERPRINT_TYPES, and the form carries that type, so [1], (1,) and
    {0: 1} differ. What goes into the form, and nothing else: None, bool and
    str as they are; an int in hex; a float's hex form, or its bits when not
    finite (NaN sign and payload kept); bytes in hex; a date in ISO form; a
    dict's sorted key-value pairs, a set's or frozenset's sorted elements, and
    a list's or tuple's elements in order, each through this function; an
    ndarray as _canonical_array reads it; a numpy scalar's dtype.str, dtype
    char and bytes (the char tells aliases such as longlong and int64 apart). Anything else, a subclass included
    (namedtuple, defaultdict, OrderedDict, a dict subclass, a dataclass,
    scipy sparse, numpy.matrix, numpy.ma.MaskedArray, pandas), raises
    TypeError."""
    kind = type(value)
    if value is None or kind is bool or kind is str:
        return value
    if kind is int:
        return {"int": hex(value)}
    if kind is float:
        return {"float": value.hex() if math.isfinite(value) else "bits:" + struct.pack(">d", value).hex()}
    if kind is bytes:
        return {"bytes": value.hex()}
    if kind is dt.date:
        return {"date": value.isoformat()}
    if kind is dict:
        return {"dict": _sorted_canonical([_canonical(k), _canonical(v)] for k, v in value.items())}
    if kind is list or kind is tuple:
        return {kind.__name__: [_canonical(v) for v in value]}
    if kind is set or kind is frozenset:
        return {kind.__name__: _sorted_canonical(_canonical(v) for v in value)}
    if np is not None and isinstance(value, np.ndarray):
        return _canonical_array(value)
    if np is not None and kind in _NUMPY_SCALARS:
        _check_dtype(value, value.dtype)
        if value.dtype.names is not None:  # a structured scalar keeps its field names
            return {"structured_scalar": _canonical_array(np.array(value))}
        if kind is np.void:
            _refuse(value, " (an unstructured void scalar)")
        return {"scalar": value.dtype.str, "bytes": value.tobytes().hex()} | {"char": value.dtype.char}
    _refuse(value)


# numpy scalar types by exact type. longdouble and clongdouble are left out:
# their width and padding differ between platforms.
_NUMPY_SCALARS = frozenset() if np is None else frozenset(
    [np.bool_, np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64, np.float16,
     np.float32, np.float64, np.complex64, np.complex128, np.str_, np.bytes_, np.void, np.datetime64,
     np.timedelta64] + [np.dtype(c).type for c in "bhilqBHILQ"])
def _check_dtype(value, dtype):
    """No dtype metadata, on the dtype, on any field's own dtype (a sub-array
    field's included, MC4-3) or on a sub-array's base (MC3-1)."""
    if dtype.metadata is not None:
        _refuse(value, " with dtype metadata")
    for field in (dtype.fields or {}).values():
        _check_dtype(value, field[0])
    if dtype.subdtype is not None:
        _check_dtype(value, dtype.subdtype[0])


def _canonical_array(value):
    """An exact ndarray (MC-1, MC4). What goes into the form, and nothing
    else: its shape as given (a 0-d array is not (1,)); for a plain dtype,
    dtype.str, the dtype char (so "q" and "l" differ) and the sha256 of the
    values' bytes in C order; for an object dtype, each element through
    _canonical; for a structured dtype, each field's name, dtype.str, offset,
    title and dtype char, the itemsize, and each field's values read again
    through this function. Any ndarray subclass (numpy.matrix,
    numpy.ma.MaskedArray, a custom one) and any dtype carrying metadata are
    refused (MC2-1, MC3-1, MC4-1)."""
    if type(value) is not np.ndarray:
        _refuse(value, " (an ndarray subclass)")
    _check_dtype(value, value.dtype)
    if value.dtype.names is not None:
        layout = [[name, value.dtype.fields[name][0].str, value.dtype.fields[name][1],
                   _canonical(value.dtype.fields[name][2] if len(value.dtype.fields[name]) > 2 else None),
                   value.dtype.fields[name][0].char]
                  for name in value.dtype.names]
        return {"structured": layout, "itemsize": value.dtype.itemsize,
                "shape": list(value.shape), "fields": [_canonical_array(value[name]) for name in value.dtype.names]}
    if value.dtype.hasobject:
        return {"ndarray": "object", "shape": list(value.shape), "items": _canonical(value.tolist())}
    data = np.ascontiguousarray(value)
    return {"ndarray": data.dtype.str, "shape": list(value.shape),
            "sha256": hashlib.sha256(data.tobytes()).hexdigest(), "char": value.dtype.char}


def fingerprint(value):
    """sha256 of a fitted value's canonical form (see _canonical): exact
    floats, array bytes in C order, sets and dicts independent of order, the
    exact type of every value, and every dtype's str and char. Types outside
    FINGERPRINT_TYPES raise TypeError."""
    text = json.dumps(_canonical(value), sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()


BARE_FIT_KEY = "<fit>"


def _fingerprints(fitted):
    """Fingerprints by key name. A fitted dict's keys must each be an exact
    str (MC3-1), so no two keys can name the same entry (NLK-2) and no key can
    carry state of its own; BARE_FIT_KEY is reserved for a fit that returns a
    bare value, so the two forms never meet (I-1)."""
    if type(fitted) is not dict:  # a Mapping subclass could hide state: fingerprinted whole, so refused
        return {BARE_FIT_KEY: fingerprint(fitted)}
    out = {}
    for key, value in fitted.items():
        if type(key) is not str:
            raise TypeError(f"fitted key {key!r} is a {type(key).__module__}.{type(key).__qualname__}; a fit's "
                            "keys must each be an exact str (MC3-1)")
        if key == BARE_FIT_KEY:
            raise ValueError(f"fitted key {BARE_FIT_KEY!r} is reserved for a fit that returns a bare value")
        out[key] = fingerprint(value)
    return out


def _changed(before, after):
    return sorted(key for key in set(before) | set(after) if before.get(key) != after.get(key))


def _features(features):
    out = {}
    for gid, genes in features.items():
        _text(gid, "features genome_id")
        if isinstance(genes, str):
            raise RecordError(f"{gid}: features must be a collection of strings, not a string")
        genes = frozenset(genes)
        for gene in genes:
            _text(gene, f"{gid} feature")
        out[gid] = genes
    return dict(sorted(out.items()))


def _training_view(records, features, cutoff, training):
    """Everything a fit at cutoff may see (ruling A): each training row's
    record with its labels recorded on or before cutoff, and its features.
    Perturbations must leave it unchanged."""
    view = {}
    for record in records:
        if record.genome_id in training:
            labels = tuple(e for e in record.labels if e.recorded_at is not None and e.recorded_at <= cutoff)
            view[record.genome_id] = (dataclasses.replace(record, labels=labels),
                                      features.get(record.genome_id))
    return view


def _flip(value):
    return "S" if value == "R" else "R"


def _later(event, cutoff):
    """Recorded after cutoff, or undated and so never knowable."""
    return event.recorded_at is None or event.recorded_at > cutoff


def _next_day(record, events, cutoff):
    """A day after cutoff and after every dated event, and not before
    collection, so an event placed there passes rule 5's load check. A label
    recorded on the last representable date leaves no such day (TL-5)."""
    latest = max([cutoff] + [e.recorded_at for e in events if e.recorded_at is not None])
    if latest == dt.date.max:
        raise RecordError(f"{record.genome_id}: an event or cutoff on {latest} leaves no later day "
                          "for a canary label event")
    return max(latest + DAY, record.collection_date or dt.date.min)


def _every_gene(features):
    return frozenset().union(*features.values()) if features else frozenset()


def _training_metadata(records, training):
    """(lineage, duplicate_group) of each training row, in genome_id order."""
    return [(r.lineage, r.duplicate_group) for r in records if r.genome_id in training]


LATER_SOURCE = "canary-resourced"


def _later_date(value):
    """DATE_SHIFT later, capped at the day before the last date so a canary
    event still fits after it, and never earlier than value."""
    cap = dt.date.max - DAY
    if value >= cap:
        return value
    try:
        return min(value + DATE_SHIFT, cap)
    except OverflowError:
        return cap


def _rewrite_later(event, record, cutoff):
    """NLK-1: a later or undated event gets the other value, a later date (a
    dated event stays after cutoff; an undated one stays undated), source
    LATER_SOURCE, and a different supersedes: null if it had one, otherwise
    the latest dated event of its antibiotic on or before cutoff, if any."""
    if event.supersedes is not None:
        supersedes = None
    else:
        earlier = [e.recorded_at for e in record.labels
                   if e.antibiotic == event.antibiotic and e.recorded_at is not None and e.recorded_at <= cutoff]
        supersedes = max(earlier) if earlier else None
    return dataclasses.replace(event, value=_flip(event.value), source_release=LATER_SOURCE, supersedes=supersedes,
                               recorded_at=None if event.recorded_at is None else _later_date(event.recorded_at))


def perturb_labels(records, features, cutoff, training):
    """Rewrite every label event recorded after cutoff or undated (value,
    date, source and supersedes; see _rewrite_later), then append one more
    revision per genome and antibiotic after the latest dated event (or the
    cutoff): later-recorded labels, later revisions and when and from where
    they were reported. Nothing on or before cutoff can name a later date,
    so moving later events keeps rule 4."""
    antibiotics = _antibiotics(records)
    out = []
    for record in records:
        labels = [_rewrite_later(e, record, cutoff) if _later(e, cutoff) else e for e in record.labels]
        for antibiotic in antibiotics:
            dated = [e for e in labels if e.antibiotic == antibiotic and e.recorded_at is not None]
            latest = max(dated, key=LabelEvent.sort_key) if dated else None
            labels.append(LabelEvent(antibiotic, _flip(latest.value) if latest else "R",
                                     _next_day(record, dated, cutoff), "canary-later-release",
                                     latest.recorded_at if latest else None))
        out.append(dataclasses.replace(record, labels=tuple(labels)))
    return tuple(out), dict(features)


def perturb_label_ties(records, features, cutoff, training):
    """Same-day R/S conflicts recorded after cutoff (rule 3): a conflicting
    twin for every event dated after cutoff, and an R and S pair on a new day
    after the latest dated event, per genome and antibiotic."""
    antibiotics = _antibiotics(records)
    out = []
    for record in records:
        labels = list(record.labels)
        labels += [LabelEvent(e.antibiotic, _flip(e.value), e.recorded_at, "canary-tie")
                   for e in record.labels if e.recorded_at is not None and e.recorded_at > cutoff]
        for antibiotic in antibiotics:
            day = _next_day(record, [e for e in record.labels if e.antibiotic == antibiotic], cutoff)
            labels += [LabelEvent(antibiotic, value, day, f"canary-tie-{value.lower()}") for value in VALUES]
        out.append(dataclasses.replace(record, labels=tuple(labels)))
    return tuple(out), dict(features)


def perturb_undated_labels(records, features, cutoff, training):
    """An undated R and an undated S event per genome and antibiotic (rule 1)."""
    antibiotics = _antibiotics(records)
    out = []
    for record in records:
        labels = list(record.labels) + [LabelEvent(antibiotic, value, None, "canary-undated")
                                        for antibiotic in antibiotics for value in VALUES]
        out.append(dataclasses.replace(record, labels=tuple(labels)))
    return tuple(out), dict(features)


def perturb_features(records, features, cutoff, training):
    """Replace the genes of every genome outside the training rows, knowable
    at cutoff or not: the symmetric difference with every gene seen, plus one
    novel gene."""
    every = _every_gene(features)
    out = dict(features)
    for gid in sorted(set(features) | {r.genome_id for r in records}):
        if gid not in training:
            out[gid] = (every ^ features.get(gid, frozenset())) | {f"canary-novel-gene-{gid}"}
    return tuple(records), out


def perturb_metadata_shared(records, features, cutoff, training):
    """Give every genome outside the training rows the lineage and
    duplicate_group of a training row, cycled in genome_id order (N-3)."""
    shared = _training_metadata(records, training) or [("canary-lineage", "canary-group")]
    out, n = [], 0
    for record in records:
        if record.genome_id not in training:
            lineage, group = shared[n % len(shared)]
            record = dataclasses.replace(record, lineage=lineage, duplicate_group=group)
            n += 1
        out.append(record)
    return tuple(out), dict(features)


def perturb_metadata_fresh(records, features, cutoff, training):
    """Give every genome outside the training rows its own new lineage and
    duplicate_group and a null collection_date (N-3). The values carry a hex
    digest of the genome id, not the id itself, so they always pass rule 6:
    printable ASCII, and distinct ids never alias (AGR-2)."""
    out = [record if record.genome_id in training else dataclasses.replace(
        record, lineage=f"canary-lineage-{_id_digest(record.genome_id)}",
        duplicate_group=f"canary-group-{_id_digest(record.genome_id)}", collection_date=None) for record in records]
    return tuple(out), dict(features)


def _id_digest(genome_id):
    return hashlib.sha256(genome_id.encode("utf-8", "surrogatepass")).hexdigest()[:16]


DATE_SHIFT = dt.timedelta(days=400)


def _shift(value, cutoff):
    """DATE_SHIFT later if after cutoff, DATE_SHIFT earlier if on or before
    it, so the date stays on its side of cutoff (clamped at the date range)."""
    if value is None:
        return None
    step = DATE_SHIFT if value > cutoff else -DATE_SHIFT
    try:
        return value + step
    except OverflowError:
        return dt.date.max if value > cutoff else dt.date.min


def perturb_dates(records, features, cutoff, training):
    """Move released_at and features_available_at of every genome outside the
    training rows by DATE_SHIFT, away from cutoff on the side it is on, and
    null its collection_date (NN-2). Knowability and labels at cutoff do not
    change, so neither does training membership."""
    out = [record if record.genome_id in training else dataclasses.replace(
        record, collection_date=None, released_at=_shift(record.released_at, cutoff),
        features_available_at=_shift(record.features_available_at, cutoff)) for record in records]
    return tuple(out), dict(features)


def add_genomes(records, features, cutoff, training):
    """Add six genomes outside the training rows: one collected and released
    after cutoff; one collected before it but released after it; one released
    by cutoff whose features came only after it; one whose features are dated
    before cutoff but whose release is after it (inconsistent_dates, available
    at the later date under rule 5); one public before cutoff whose only
    label came after it (ruling A: knowable but unlabelled); and one dated a
    day after every date in the records (NN-2). The second, third and fourth
    carry labels recorded on or before cutoff. Each copies the lineage and
    duplicate_group of a training row, cycled (N-3)."""
    shapes = [  # collection_date, released_at, features_available_at, label value, label recorded_at
        (cutoff + DAY, cutoff + DAY, None, "R", cutoff + DAY),
        (cutoff - DAY, cutoff + DAY, None, "S", cutoff),
        (cutoff - DAY, cutoff, cutoff + DAY, "R", cutoff - DAY),
        (cutoff - DAY, cutoff + DAY, cutoff - DAY, "S", cutoff - DAY),
        (cutoff - DAY, cutoff - DAY, None, "R", cutoff + DAY),
    ]
    latest = max([cutoff] + [d for r in records for d in (r.collection_date, r.released_at, r.features_available_at)
                             if d is not None] + [e.recorded_at for r in records for e in r.labels if e.recorded_at])
    last = latest + DAY if latest < dt.date.max else latest
    shapes.append((last, last, None, "R", last))
    taken = {r.genome_id for r in records} | set(features)
    every = _every_gene(features)
    ids = []
    n = 0
    while len(ids) < len(shapes):
        n += 1
        if f"canary-new-{n}" not in taken:
            ids.append(f"canary-new-{n}")
    shared = _training_metadata(records, training) or [("canary-lineage", "canary-group")]
    antibiotics = _antibiotics(records)
    added = []
    out = dict(features)
    for i, (gid, (collected, released, available, value, recorded)) in enumerate(zip(ids, shapes)):
        labels = tuple(LabelEvent(antibiotic, value, recorded, "canary-later-release")
                       for antibiotic in antibiotics)
        lineage, group = shared[i % len(shared)]
        added.append(GenomeRecord(gid, collected, released, available, lineage, group, labels))
        out[gid] = every | {f"canary-novel-gene-{gid}"}
    return tuple(sorted(tuple(records) + tuple(added), key=lambda r: r.genome_id)), out


def perturb_training_labels(records, features, cutoff, training):
    """Positive control (N-2): flip every training row's labels recorded on
    or before cutoff."""
    return tuple(record if record.genome_id not in training else dataclasses.replace(record, labels=tuple(
        e if _later(e, cutoff) else dataclasses.replace(e, value=_flip(e.value)) for e in record.labels))
        for record in records), dict(features)


def perturb_training_features(records, features, cutoff, training):
    """Positive control (N-2): replace every training row's genes like
    perturb_features does for the other genomes."""
    every = _every_gene(features)
    changed = dict(features)
    for gid in sorted(training):
        changed[gid] = (every ^ features.get(gid, frozenset())) | {f"canary-novel-gene-{gid}"}
    return tuple(records), changed


def drop_training_row(records, features, cutoff, training):
    """Positive control (NN-1): remove the last training row by genome_id,
    record and features, so counts and row lists move. With fewer than two
    training rows it changes nothing."""
    if len(training) < 2:
        return tuple(records), dict(features)
    gone = max(training)
    return (tuple(r for r in records if r.genome_id != gone),
            {gid: genes for gid, genes in features.items() if gid != gone})


# Applied one at a time: flipping labels and complementing genes together can
# cancel out exactly, e.g. in a centred linear score X'(y - mean(y)).
TRAINING_CONTROLS = {"training_labels": perturb_training_labels, "training_features": perturb_training_features,
                     "training_rows": drop_training_row}


def drop_post_cutoff_events(records, features, cutoff, training):
    """Remove every label event recorded after cutoff, and every undated event,
    from every genome, training rows included (NLK-1). By ruling A nothing
    recorded after cutoff may influence the fit, so the fit must be the same
    without it: not its value, date, source or supersedes."""
    return (tuple(dataclasses.replace(r, labels=tuple(e for e in r.labels if not _later(e, cutoff)))
                  for r in records), dict(features))


def drop_non_training(records, features, cutoff, training):
    """Remove every genome outside the training rows, record and features
    (TL-1, TL-2). Ruling A taken literally: those genomes must not influence
    the fit at all, so the fit must be the same without them. With
    training_view_only, it is one of two perturbations that reach labels
    recorded on or before cutoff on genomes that are never knowable, not yet
    knowable, or ambiguous at cutoff; the others leave those alone."""
    return (tuple(r for r in records if r.genome_id in training),
            {gid: genes for gid, genes in features.items() if gid in training})


def training_view_only(records, features, cutoff, training):
    """Exactly the training view (LEAK-1): the training rows, each with only
    its dated events recorded on or before cutoff, and only their features.
    Each drop above leaves the other kind of test-period data in place, so a
    fit that asks only whether any test-period data exists gives the same
    answer under both; with all of it gone at once, it must not change."""
    view = _training_view(records, features, cutoff, training)
    return (tuple(record for record, _ in view.values()),
            {gid: genes for gid, (_, genes) in view.items() if genes is not None})


PERTURBATIONS = {"labels": perturb_labels, "label_ties": perturb_label_ties,
                 "undated_labels": perturb_undated_labels, "features": perturb_features,
                 "metadata_shared": perturb_metadata_shared, "metadata_fresh": perturb_metadata_fresh,
                 "dates": perturb_dates, "new_genomes": add_genomes, "drop_non_training": drop_non_training,
                 "drop_post_cutoff_events": drop_post_cutoff_events, "training_view_only": training_view_only}


def invariance_to_test_data(records, features, cutoff, fit, perturbations=None, antibiotics=None,
                            constant_keys=()):
    """Changing test-period records must not change anything fitted at cutoff.

    fit(records, features, cutoff) receives every record (a tuple of
    GenomeRecords sorted by genome_id) and every genome's features
    (genome_id -> frozenset), and must itself keep to the training rows of
    ruling A (see the module docstring). It returns plain data or numpy arrays;
    a dict is compared key by key (e.g. vocabulary, encoders, selected
    features, model), and its keys must each be an exact str, else TypeError
    (MC3-1, which supersedes NLK-2's same-str() refusal); anything else is
    the single key "<fit>", which a dict may not use. A key that only one of
    two runs returns counts as changed. The canary guards against honest fits that accidentally read test-period
    data. It doesn't defend against a fit built to hide information on purpose in
    private (underscore) attributes, object internals, memory layout or padding
    bytes.

    The fingerprint fails closed (MC2-1): every fitted value, and every value
    inside it, must have exactly one of these types: None, bool, int, float,
    str, bytes, dict, list, tuple, set, frozenset, datetime.date, an exact
    numpy.ndarray (structured too), and the numpy scalar types other than
    longdouble and clongdouble (a void scalar only if structured). Anything
    else raises TypeError naming its type, and the canary then fails:
    subclasses of these (namedtuple, defaultdict, OrderedDict, any dict
    subclass, numpy.matrix, numpy.ma.MaskedArray), dataclasses, scipy sparse
    and pandas objects, and any dtype carrying metadata, on an array, a
    scalar, a field, a sub-array field or a sub-array's base. What is
    fingerprinted, and nothing else (see _canonical and _canonical_array):
    each value's exact type; the value forms of the scalars and containers;
    an array's shape and its values' bytes in C order, or each object
    element; every dtype's str and char, so aliases such as "q" and "l"
    differ; and a structured dtype's field names, offsets, titles and
    itemsize, with each field's values. Not read (I-3): whether an array is a
    view (its .base), its strides, C or Fortran order, the writeable and
    aligned flags, a structured dtype's aligned-struct flag and its padding
    bytes. Padding can be filled by accident (np.empty over a reused buffer,
    and np.save keeps it), but it never reaches a value a model reads. The fit
    must be
    a pure function of its three inputs: it is called several times, and a fit
    that caches across calls (by cutoff, path or anything else) cannot be
    certified. In particular it must not cache on anything outside its
    training view: a cache keyed on the training rows alone, whose value was
    computed from other genomes, looks clean to every black-box check here.

    antibiotics names the antibiotics being fitted; a training row needs an
    unambiguous label as of cutoff for one of them. None means every
    antibiotic in the records. training_view_only removes every antibiotic
    that only test-period records carry, and the drops can remove one too, so
    a fit that names its label columns from every record's antibiotics fails
    when that changes its output.

    constant_keys names fitted keys that legitimately do not depend on the
    training rows: hyperparameters, schema versions, or a vocabulary taken
    from the artifact registry at cutoff. Default empty. They are exempt from
    the positive control, and must be identical across every run, controls
    included. Every name in constant_keys is a claim the reviewer must check:
    no canary can tell a truly constant key from one cached on the cutoff.

    1. The fit is repeated on unchanged inputs. If it differs, the canary
       fails with ("refit", key, ...) and attributes nothing to test data.
    2. Positive control, per key: each of TRAINING_CONTROLS (flip the training
       rows' labels; replace their genes; drop one training row) is applied
       alone. Every fitted key outside constant_keys that no control changes
       is an offender ("control", key, ...), because a key that ignores its
       own training data (or replays a cached answer) would pass every
       invariance check vacuously. A constant key that a control changes, or
       that the fit does not return, is an offender too. With no training rows
       the canary fails with ("control", "*", ...) and runs nothing else.
    3. Each perturbation is applied alone, and every fitted key it changes is
       an offender (perturbation, key, ...). A fit that raises an Exception on
       the perturbed data depends on that data, so it is an offender too:
       (perturbation, "*", "fit raised <type>: <message, at most 200
       characters>"). A perturbation that alters the training rows raises
       ValueError: that is a canary bug. Exceptions from the fit on unchanged
       inputs, from the controls, from the guard and from the perturbations
       themselves propagate, so a broken setup is never counted as a leak.

    The built-in PERTURBATIONS, each touching only data outside the training
    rows: labels (later labels and revisions), label_ties (same-day R/S
    conflicts after cutoff), undated_labels, features, metadata_shared and
    metadata_fresh (lineage, duplicate_group and collection_date), dates
    (released_at and features_available_at moved, collection_date nulled),
    new_genomes, drop_non_training (every genome outside the training rows
    removed), drop_post_cutoff_events (every undated event and every event
    recorded after cutoff removed, training rows included), and
    training_view_only (both at once: exactly the training view, so a fit
    whose output differs from its output on the training view fails). A
    perturbation is called as perturbation(records, features, cutoff, training) and returns
    (records, features). A pass covers the input the canary ran on: a fit
    whose dependence on test-period data shows only at other sizes or counts
    (NLK-5) is not caught. With the default perturbations and at least one
    training row, a cutoff of 0001-01-01 raises OverflowError from
    new_genomes, which places a genome the day before it (I-2), so the canary
    never passes there. With no training row there, it fails with ("control",
    "*", ...) and raises nothing; with a perturbation set that leaves out
    new_genomes, it runs, and an honest fit can pass.
    """
    cutoff = _date(cutoff, "cutoff", False)
    records = load_records(records)
    features = _features(features)
    perturbations = PERTURBATIONS if perturbations is None else perturbations
    antibiotics = _antibiotics(records) if antibiotics is None else sorted(antibiotics)
    if isinstance(constant_keys, str):
        raise ValueError("constant_keys must be a collection of key names, not a string")
    constant = frozenset(str(key) for key in constant_keys)
    training = training_ids(records, cutoff, antibiotics)

    def run(run_records, run_features):
        return _fingerprints(fit(tuple(run_records), dict(run_features), cutoff))

    baseline = run(records, features)
    offenders = [("refit", key, "differs on unchanged inputs") for key in _changed(baseline, run(records, features))]
    if not offenders and not training:
        offenders.append(("control", "*", "no training rows at the cutoff; invariance not certified"))
    elif not offenders:
        moved = set()
        for _, control in sorted(TRAINING_CONTROLS.items()):
            moved.update(_changed(baseline, run(*control(records, features, cutoff, training))))
        offenders += [("control", key, "insensitive to its own training data; not certified")
                      for key in set(baseline) - constant - moved]
        offenders += [("control", key, "named in constant_keys but changed by training data")
                      for key in constant & moved]
        offenders += [("control", key, "named in constant_keys but not fitted") for key in constant - set(baseline)]
        view = _training_view(records, features, cutoff, training)
        for name in sorted(perturbations):
            changed_records, changed_features = perturbations[name](records, features, cutoff, training)
            changed_records = load_records(changed_records)
            changed_features = _features(changed_features)
            if (training_ids(changed_records, cutoff, antibiotics) != training
                    or _training_view(changed_records, changed_features, cutoff, training) != view):
                raise ValueError(f"perturbation {name!r} changed the training rows at {cutoff}")
            try:
                fitted = run(changed_records, changed_features)
            except Exception as error:  # NL9: the fit depends on what the perturbation changed
                offenders.append((name, "*", f"fit raised {type(error).__name__}: {str(error)[:200]}"))
                continue
            offenders += [(name, key, "changed by test-period records") for key in _changed(baseline, fitted)]
    return _result("test_data_invariance", offenders,
                   f"{len(set(baseline) - constant)} fitted key(s) respond to training data, "
                   f"{len(constant)} declared constant; all unchanged by "
                   f"{len(perturbations)} test-period perturbation(s)",
                   f"fit at {cutoff} depends on test-period records, ignores its training data or is not "
                   "deterministic", constant)


# --- unseen vocabulary -----------------------------------------------------

def unseen_vocabulary(vocabulary, test_features, encode=None, policy="unknown"):
    """Test genes outside the cutoff's vocabulary are encoded deterministically
    under the documented policy and never join the vocabulary.

    encode(genes, vocabulary) defaults to artifacts.encode_genes under policy.
    It receives a mutable copy of the vocabulary so additions are caught. Each
    genome is encoded from sorted and reversed gene orders and once more. Every
    vocabulary token in the encoding must be a gene the genome carries, so a
    remap of an unseen allele onto a known gene is an offender. Tokens are
    compared as a multiset: each may appear once (the documented policies
    emit each known gene and UNKNOWN at most once), so an unseen allele
    counted as a second copy of a carried gene is an offender too.
    Offenders are (genome_id or "<vocabulary>", gene or "*", reason).
    """
    if policy not in artifacts.UNSEEN_POLICIES:
        raise ArtifactError(f"policy must be one of {list(artifacts.UNSEEN_POLICIES)}, got {policy!r}")
    if encode is None:
        def encode(genes, vocab):
            return artifacts.encode_genes(genes, vocab, policy)
    frozen = frozenset(vocabulary)
    working = set(frozen)
    offenders = []
    for gid, genes in _features(test_features).items():
        ordered = sorted(genes)
        outputs = [tuple(encode(list(ordered), working)), tuple(encode(list(reversed(ordered)), working)),
                   tuple(encode(list(ordered), working))]
        if len(set(outputs)) != 1:
            offenders.append((gid, "*", "encoding depends on gene order or repetition"))
        tokens = outputs[0]
        unseen = genes - frozen
        offenders += [(gid, token, "token outside the vocabulary")
                      for token in set(tokens) - frozen - {UNKNOWN}]
        offenders += [(gid, gene, "known gene not encoded") for gene in (genes & frozen) - set(tokens)]
        offenders += [(gid, token, "known token the genome does not carry")
                      for token in (set(tokens) & frozen) - genes]
        offenders += [(gid, token, "token repeated") for token in set(tokens) if tokens.count(token) > 1]
        if policy == "unknown" and unseen and UNKNOWN not in tokens:
            offenders.append((gid, "*", f"{len(unseen)} unseen gene(s) not mapped to {UNKNOWN}"))
        if UNKNOWN in tokens and (policy == "drop" or not unseen):
            offenders.append((gid, "*", f"{UNKNOWN} emitted under policy {policy!r} "
                                        f"with {len(unseen)} unseen gene(s)"))
    offenders += [("<vocabulary>", gene, "added while encoding test genomes") for gene in working - frozen]
    offenders += [("<vocabulary>", gene, "removed while encoding test genomes") for gene in frozen - working]
    return _result("unseen_vocabulary", offenders,
                   f"unseen genes handled by policy {policy!r}; vocabulary unchanged",
                   "unseen genes leak into encodings or the vocabulary")
