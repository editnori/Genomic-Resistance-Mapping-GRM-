"""Cutoff splits for temporal later-label benchmarking.

Train only on what was knowable at a cutoff T; test on genomes and labels that
appeared later. A split supports later-label prediction claims only.

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

- Rule 1: LabelAsOf ignores undated events at every date. A record whose
  antibiotic has no dated event on or before `eval_date` but has an undated
  one is excluded as `unknown_label_recorded_at` instead of `no_label`. An
  undated event never changes a label read from dated events.
- Rule 3: an as-of read reads the latest dated day in its window. A conflict
  on that day makes the read ambiguous; a conflict on an earlier day that a
  later dated event in the same window superseded does not.
- Rule 5: when only one of `released_at` and `features_available_at` is
  present, it is the availability date. `inconsistent_dates` is a flag, not an
  exclusion; the record is placed using the later date.
- A missing required key or an unknown key raises `ValueError`. Required keys
  are `genome_id` and `labels` for a record, and `antibiotic`, `value`,
  `recorded_at` and `source_release` for an event. A missing optional key is null.
- Rule 6: `str.strip()` also removes surrounding Unicode whitespace such as
  a tab or a no-break space (U+00A0). An interior tab, no-break space,
  zero-width space or any other non-ASCII character is refused. Loaded values
  are ASCII, so the alias check amounts to case-insensitive equality. It
  compares lineages with lineages and duplicate groups with duplicate groups
  across the whole input. Values otherwise match as exact strings, so `4.2`
  and `4.2.1` are different lineages.
- `T`, `horizon` or `embargo` that push a date out of `datetime.date`'s range
  raise `ValueError`.

Python names: `GenomeRecord`, `LabelEvent`, and `label_as_of(record, antibiotic, T)`
for LabelAsOf, which raises `AmbiguousLabel` (a `ValueError`) under rule 3.
Dates are `YYYY-MM-DD` strings or `datetime.date`.

Availability contract of `split_at(records, T, horizon, embargo, eval_date,
antibiotic, policy)`. `horizon` and `embargo` are whole days, and `eval_date`
must be on or after T + horizon. Each record lands in exactly one place,
decided from its own fields, in this order:

1. No dated event for the antibiotic on or before `eval_date`:
   `unknown_label_recorded_at` when an undated event exists, else `no_label`.
2. `unknown_features_available_at`: `features_available_at` and `released_at` are null.
3. A dated event on or before T: `ambiguous_label` when LabelAsOf(T) is
   ambiguous. Otherwise **train** with the LabelAsOf(T) value when features
   became available on or before T - embargo, `embargo` when they became
   available in (T - embargo, T], or `features_after_cutoff` when after T.
   A label known at T is not a later label, so it never enters test.
4. The first dated event is after T: `features_after_eval` when features
   became available after `eval_date`; `unknown_collection_date` or
   `unknown_released_at` when the policy's window date is null;
   `ambiguous_label` when LabelAsOf(`eval_date`) is ambiguous.
5. Window date on or before T: backfill. The genome is flagged `backfill` and,
   under `backfill="exclude"`, excluded with reason `backfill`. Under
   `backfill="test"` it enters **test** when its first dated event is in
   (T, T + horizon], else `backfill_outside_horizon`.
6. Window date in (T, T + horizon]: **test** with the LabelAsOf(`eval_date`)
   value. Later window dates are `after_horizon`.

`SplitPolicy.window` picks the window date: `collection` (the default) uses
`collection_date`, and `release` uses `released_at`. Train never reads the
window date, so a null `collection_date` alone does not remove a train genome.

Train membership and train labels depend only on a record's availability date
and its dated events recorded on or before T. Undated events, revisions and
same-day conflicts recorded after T, `eval_date`, and edits to test records
cannot change train. As T moves later with the other parameters fixed, train
keeps every earlier train genome except one whose LabelAsOf at the later
cutoff is ambiguous.

The strict novel-group arm keeps a test genome only when its lineage and
duplicate group are both known, neither appears in that split's train, and
every train genome has a known lineage and duplicate group. Otherwise the test
genome is excluded from the arm as, in this order, `lineage_seen`,
`duplicate_group_seen`, `group_unknown` (its own lineage or group is null), or
`novel_unverified` (unseen, but some train genome has a null lineage or group
and could share it). `counts.train_group_unknown` counts those train genomes.

What the test set measures, also listed in each manifest's `notes`:

- Revisions recorded after T are never tested. A genome with a label at T is
  train or excluded, so test holds first labels that appeared after T.
- A genome phenotyped on or before T but available after T reaches neither
  arm. It is counted as `features_after_cutoff`.
- Test labels are read as of `eval_date` and train labels as of T, so test
  labels have had longer to be revised.
- Test membership and test labels depend on events recorded in
  (T + horizon, `eval_date`]. A genome whose window date is in the horizon
  enters test on a first label recorded up to `eval_date`, while a backfill
  genome needs its first label in the horizon.
- A test genome with an undated event for the antibiotic enters test on its
  dated label (rule 1), even if the undated event was in fact recorded on or
  before T.

The manifest lists train and test ids with their labels, the novel-group ids,
excluded ids and counts by reason, the `backfill` and `inconsistent_dates`
flags, the parameters, the notes, and `input_sha256`. `canonical_json(obj)` is
`json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
allow_nan=False)`, UTF-8 encoded.
`input_sha256` is the SHA-256 of `canonical_json` of the input records sorted
by `genome_id`, each with all seven keys as loaded (dates ISO, lineage and
duplicate group stripped, `features_available_at` not defaulted) and its labels sorted by
(`antibiotic`, `recorded_at`, `value`, `source_release`, `supersedes`), null
dates sorting first. `sha256` is the SHA-256 of `canonical_json` of every
other manifest key.
"""
import datetime as dt
import hashlib
import json
import re
import unicodedata
from collections.abc import Mapping
from dataclasses import MISSING, dataclass, fields

SCHEMA = "grm.temporal.split.v2"
WINDOWS = ("collection", "release")
BACKFILL_POLICIES = ("exclude", "test")
EXCLUSION_REASONS = (
    "unknown_label_recorded_at",
    "no_label",
    "unknown_features_available_at",
    "ambiguous_label",
    "embargo",
    "features_after_cutoff",
    "features_after_eval",
    "unknown_collection_date",
    "unknown_released_at",
    "backfill",
    "backfill_outside_horizon",
    "after_horizon",
)
NOVEL_GROUP_EXCLUSIONS = (
    "lineage_seen", "duplicate_group_seen", "group_unknown", "novel_unverified")
NOTES = (
    "Revisions recorded after T are never tested: a genome with a label at T is train "
    "or excluded, so test holds first labels that appeared after T.",
    "A genome phenotyped on or before T but available after T reaches neither arm; "
    "see counts.excluded.features_after_cutoff.",
    "Test labels are read as of eval_date and train labels as of T, so test labels "
    "have had longer to be revised.",
    "Test membership and test labels depend on events recorded in (T + horizon, eval_date]: "
    "a genome whose window date is in the horizon enters test on a first label recorded up "
    "to eval_date, while a backfill genome needs its first label in the horizon.",
    "A test genome with an undated event for the antibiotic enters test on its dated label "
    "(rule 1), even if the undated event was in fact recorded on or before T.",
)
ISO_DATE = re.compile(r"\d{4}-\d{2}-\d{2}")


class AmbiguousLabel(ValueError):
    """The day an as-of read uses holds events with different values."""


def _date(value, name):
    if value is None:
        return None
    if isinstance(value, dt.datetime):
        raise ValueError(f"{name} must be a date, not a datetime: {value!r}")
    if isinstance(value, dt.date):
        return value
    if isinstance(value, str) and ISO_DATE.fullmatch(value):
        return dt.date.fromisoformat(value)
    raise ValueError(f"{name} must be an ISO date (YYYY-MM-DD) or null: {value!r}")


def _required_date(value, name):
    if value is None:
        raise ValueError(f"{name} is required")
    return _date(value, name)


def _iso(value):
    return None if value is None else value.isoformat()


def _group_name(genome_id, value, name):
    """Rule 6: strip a non-null lineage or duplicate group; refuse it unless printable ASCII."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"{genome_id}: {name} must be a string or null: {value!r}")
    value = value.strip()
    if not value or any(not " " <= char <= "~" for char in value):
        raise ValueError(f"{genome_id}: {name} must be non-empty printable ASCII after "
                         f"stripping: {value!r}")
    return value


def _refuse_aliases(records):
    """Rule 6: refuse two distinct values of one field that are equal under NFKC casefold."""
    for name in ("lineage", "duplicate_group"):
        seen = {}
        for value in sorted({getattr(record, name) for record in records} - {None}):
            key = unicodedata.normalize("NFKC", value).casefold()
            if key in seen:
                raise ValueError(f"{name} values {seen[key]!r} and {value!r} collide "
                                 f"under NFKC casefold")
            seen[key] = value


def _from_mapping(cls, row):
    if not isinstance(row, Mapping):
        raise ValueError(f"{cls.__name__} must be a mapping: {row!r}")
    names = {field.name for field in fields(cls)}
    required = {field.name for field in fields(cls) if field.default is MISSING}
    unknown, missing = set(row) - names, required - set(row)
    if unknown or missing:
        raise ValueError(f"{cls.__name__} unknown keys {sorted(unknown)}, "
                         f"missing keys {sorted(missing)}")
    return cls(**row)


@dataclass(frozen=True, kw_only=True)
class LabelEvent:
    antibiotic: str
    value: str
    recorded_at: dt.date | None
    source_release: str
    supersedes: dt.date | None = None

    def __post_init__(self):
        if not isinstance(self.antibiotic, str) or not self.antibiotic:
            raise ValueError(f"antibiotic must be a non-empty string: {self.antibiotic!r}")
        if self.value not in ("R", "S"):
            raise ValueError(f"value must be R or S: {self.value!r}")
        if not isinstance(self.source_release, str) or not self.source_release:
            raise ValueError(f"source_release must be a non-empty string: {self.source_release!r}")
        recorded = _date(self.recorded_at, "recorded_at")
        supersedes = _date(self.supersedes, "supersedes")
        if recorded is not None and supersedes is not None and supersedes >= recorded:
            raise ValueError(f"a revision recorded {recorded} cannot supersede {supersedes}")
        object.__setattr__(self, "recorded_at", recorded)
        object.__setattr__(self, "supersedes", supersedes)

    @classmethod
    def from_mapping(cls, row):
        return _from_mapping(cls, row)

    def to_mapping(self):
        return {"antibiotic": self.antibiotic, "value": self.value,
                "recorded_at": _iso(self.recorded_at), "source_release": self.source_release,
                "supersedes": _iso(self.supersedes)}


@dataclass(frozen=True, kw_only=True)
class GenomeRecord:
    genome_id: str
    collection_date: dt.date | None = None
    released_at: dt.date | None = None
    features_available_at: dt.date | None = None
    lineage: str | None = None
    duplicate_group: str | None = None
    labels: tuple

    def __post_init__(self):
        if not isinstance(self.genome_id, str) or not self.genome_id:
            raise ValueError(f"genome_id must be a non-empty string: {self.genome_id!r}")
        for name in ("collection_date", "released_at", "features_available_at"):
            object.__setattr__(self, name, _date(getattr(self, name), name))
        for name in ("lineage", "duplicate_group"):
            object.__setattr__(self, name, _group_name(self.genome_id, getattr(self, name), name))
        if not isinstance(self.labels, (list, tuple)):
            raise ValueError(f"{self.genome_id}: labels must be a list: {self.labels!r}")
        labels = tuple(
            event if isinstance(event, LabelEvent) else LabelEvent.from_mapping(event)
            for event in self.labels
        )
        dated = {(event.antibiotic, event.recorded_at) for event in labels
                 if event.recorded_at is not None}
        for event in labels:
            if (event.recorded_at is not None and self.collection_date is not None
                    and event.recorded_at < self.collection_date):
                raise ValueError(f"{self.genome_id}: {event.antibiotic} label recorded "
                                 f"{event.recorded_at} before collection {self.collection_date}")
            if (event.supersedes is not None
                    and (event.antibiotic, event.supersedes) not in dated):
                raise ValueError(f"{self.genome_id}: {event.antibiotic} event supersedes "
                                 f"{event.supersedes}, which has no {event.antibiotic} event")
        object.__setattr__(self, "labels", labels)

    @property
    def available_at(self):
        """The later of `released_at` and `features_available_at` that is present."""
        present = [day for day in (self.released_at, self.features_available_at)
                   if day is not None]
        return max(present, default=None)

    @property
    def inconsistent_dates(self):
        """Features are dated before public release."""
        return (self.released_at is not None and self.features_available_at is not None
                and self.features_available_at < self.released_at)

    @classmethod
    def from_mapping(cls, row):
        return _from_mapping(cls, row)

    def to_mapping(self):
        labels = sorted((event.to_mapping() for event in self.labels), key=lambda row: (
            row["antibiotic"], row["recorded_at"] or "", row["value"],
            row["source_release"], row["supersedes"] or ""))
        return {"genome_id": self.genome_id, "collection_date": _iso(self.collection_date),
                "released_at": _iso(self.released_at),
                "features_available_at": _iso(self.features_available_at),
                "lineage": self.lineage, "duplicate_group": self.duplicate_group,
                "labels": labels}


@dataclass(frozen=True)
class SplitPolicy:
    window: str = "collection"
    backfill: str = "exclude"

    def __post_init__(self):
        if self.window not in WINDOWS:
            raise ValueError(f"window must be one of {WINDOWS}: {self.window!r}")
        if self.backfill not in BACKFILL_POLICIES:
            raise ValueError(f"backfill must be one of {BACKFILL_POLICIES}: {self.backfill!r}")


def _events(record, antibiotic):
    return [event for event in record.labels if event.antibiotic == antibiotic]


def label_as_of(record, antibiotic, T):
    """LabelAsOf: the latest dated event for `antibiotic` with recorded_at <= T, or None."""
    T = _required_date(T, "T")
    known = [event for event in _events(record, antibiotic)
             if event.recorded_at is not None and event.recorded_at <= T]
    if not known:
        return None
    latest = max(event.recorded_at for event in known)
    tied = [event for event in known if event.recorded_at == latest]
    if len({event.value for event in tied}) > 1:
        raise AmbiguousLabel(f"{record.genome_id}: {antibiotic} has conflicting labels on {latest}")
    return min(tied, key=lambda event: (event.source_release, event.supersedes is not None,
                                        event.supersedes or dt.date.min))


def canonical_json(obj):
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False)


def _sha256(obj):
    return hashlib.sha256(canonical_json(obj).encode("utf-8")).hexdigest()


def manifest_sha256(manifest):
    """SHA-256 of the canonical JSON of every manifest key except `sha256`."""
    return _sha256({key: value for key, value in manifest.items() if key != "sha256"})


def input_sha256(records):
    """SHA-256 of the canonical JSON of the records, as defined in the module docstring."""
    records = _records(records)
    return _sha256([record.to_mapping()
                    for record in sorted(records, key=lambda record: record.genome_id)])


def _records(records):
    out = [record if isinstance(record, GenomeRecord) else GenomeRecord.from_mapping(record)
           for record in records]
    seen, duplicates = set(), set()
    for record in out:
        (duplicates if record.genome_id in seen else seen).add(record.genome_id)
    if duplicates:
        raise ValueError(f"duplicate genome_id: {sorted(duplicates)}")
    _refuse_aliases(out)
    return out


def _days(value, name, minimum):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer number of days >= {minimum}: {value!r}")
    return dt.timedelta(days=value)


def _place(record, T, train_end, horizon_end, eval_date, antibiotic, policy):
    """Return (arm, label value or exclusion reason, backfill flag) for one record."""
    events = _events(record, antibiotic)
    dated = sorted(event.recorded_at for event in events if event.recorded_at is not None)
    if not dated or dated[0] > eval_date:
        if any(event.recorded_at is None for event in events):
            return "excluded", "unknown_label_recorded_at", False
        return "excluded", "no_label", False
    available = record.available_at
    if available is None:
        return "excluded", "unknown_features_available_at", False

    if dated[0] <= T:
        try:
            value = label_as_of(record, antibiotic, T).value
        except AmbiguousLabel:
            return "excluded", "ambiguous_label", False
        if available <= train_end:
            return "train", value, False
        return "excluded", "embargo" if available <= T else "features_after_cutoff", False

    if available > eval_date:
        return "excluded", "features_after_eval", False
    if policy.window == "collection":
        window_date, unknown = record.collection_date, "unknown_collection_date"
    else:
        window_date, unknown = record.released_at, "unknown_released_at"
    if window_date is None:
        return "excluded", unknown, False
    try:
        value = label_as_of(record, antibiotic, eval_date).value
    except AmbiguousLabel:
        return "excluded", "ambiguous_label", False
    if window_date <= T:
        if policy.backfill == "exclude":
            return "excluded", "backfill", True
        if dated[0] <= horizon_end:
            return "test", value, True
        return "excluded", "backfill_outside_horizon", True
    if window_date <= horizon_end:
        return "test", value, False
    return "excluded", "after_horizon", False


def split_at(records, T, horizon, embargo, eval_date, antibiotic, policy=SplitPolicy()):
    """Split records at cutoff T and return the split manifest (see module docstring)."""
    records = _records(records)
    T = _required_date(T, "T")
    eval_date = _required_date(eval_date, "eval_date")
    try:
        horizon_end = T + _days(horizon, "horizon", 1)
        train_end = T - _days(embargo, "embargo", 0)
    except OverflowError as error:
        raise ValueError(f"T {T}, horizon {horizon} or embargo {embargo} leaves the "
                         f"supported date range: {error}") from error
    if eval_date < horizon_end:
        raise ValueError(f"eval_date {eval_date} is before T + horizon {horizon_end}")
    if not isinstance(antibiotic, str) or not antibiotic:
        raise ValueError(f"antibiotic must be a non-empty string: {antibiotic!r}")
    if not isinstance(policy, SplitPolicy):
        raise ValueError(f"policy must be a SplitPolicy: {policy!r}")

    arms = {"train": {}, "test": {}}
    excluded = {reason: [] for reason in EXCLUSION_REASONS}
    backfill = []
    for record in records:
        arm, outcome, is_backfill = _place(
            record, T, train_end, horizon_end, eval_date, antibiotic, policy)
        if arm == "excluded":
            excluded[outcome].append(record.genome_id)
        else:
            arms[arm][record.genome_id] = outcome
        if is_backfill:
            backfill.append(record.genome_id)

    by_id = {record.genome_id: record for record in records}
    train_records = [by_id[genome_id] for genome_id in arms["train"]]
    train_lineages = {record.lineage for record in train_records} - {None}
    train_groups = {record.duplicate_group for record in train_records} - {None}
    train_group_unknown = sum(record.lineage is None or record.duplicate_group is None
                              for record in train_records)
    novel = []
    novel_excluded = {reason: [] for reason in NOVEL_GROUP_EXCLUSIONS}
    for genome_id in arms["test"]:
        record = by_id[genome_id]
        if record.lineage in train_lineages:
            novel_excluded["lineage_seen"].append(genome_id)
        elif record.duplicate_group in train_groups:
            novel_excluded["duplicate_group_seen"].append(genome_id)
        elif record.lineage is None or record.duplicate_group is None:
            novel_excluded["group_unknown"].append(genome_id)
        elif train_group_unknown:
            novel_excluded["novel_unverified"].append(genome_id)
        else:
            novel.append(genome_id)
    inconsistent = sorted(record.genome_id for record in records if record.inconsistent_dates)

    manifest = {
        "schema": SCHEMA,
        "parameters": {
            "antibiotic": antibiotic,
            "T": T.isoformat(),
            "horizon_days": horizon,
            "embargo_days": embargo,
            "eval_date": eval_date.isoformat(),
            "window": policy.window,
            "backfill": policy.backfill,
        },
        "input_sha256": input_sha256(records),
        "notes": list(NOTES),
        "train": arms["train"],
        "test": arms["test"],
        "novel_group_test": sorted(novel),
        "excluded": {reason: sorted(ids) for reason, ids in excluded.items()},
        "novel_group_excluded": {reason: sorted(ids) for reason, ids in novel_excluded.items()},
        "flags": {"backfill": sorted(backfill), "inconsistent_dates": inconsistent},
        "counts": {
            "records": len(records),
            "train": len(arms["train"]),
            "test": len(arms["test"]),
            "novel_group_test": len(novel),
            "backfill": len(backfill),
            "inconsistent_dates": len(inconsistent),
            "train_group_unknown": train_group_unknown,
            "excluded": {reason: len(ids) for reason, ids in excluded.items()},
            "novel_group_excluded": {
                reason: len(ids) for reason, ids in novel_excluded.items()},
        },
    }
    manifest["sha256"] = manifest_sha256(manifest)
    return manifest


def rolling_origin(records, cutoffs, horizon, embargo, eval_date, antibiotic,
                   policy=SplitPolicy()):
    """Return one `split_at` manifest per cutoff; cutoffs must strictly increase."""
    records = _records(records)
    cutoffs = [_required_date(cutoff, "cutoff") for cutoff in cutoffs]
    if not cutoffs:
        raise ValueError("cutoffs is empty")
    if any(later <= earlier for earlier, later in zip(cutoffs, cutoffs[1:])):
        raise ValueError(f"cutoffs must strictly increase: {[c.isoformat() for c in cutoffs]}")
    return [split_at(records, cutoff, horizon, embargo, eval_date, antibiotic, policy)
            for cutoff in cutoffs]
