"""Reference canonical k-mer counter, .kset v1 codec, and exact parity checks.

Matches GRM's Rust semantics in src/dna.rs, src/fnv.rs, src/kset.rs,
src/fasta.rs, src/matrix.rs, and src/input_build.rs using Python's standard
library only.
"""

from collections import Counter
from dataclasses import dataclass
import gzip
import hashlib
import pathlib
import struct
from typing import Dict, Iterable, Iterator, List, Mapping, Optional, Sequence, Set, Tuple, Union


MIN_K = 5
MAX_K = 32
BUCKET_BITS = 10
BUCKETS = 1 << BUCKET_BITS
KSET_MAGIC = b"GRMKSET1"

MASK64 = 0xFFFF_FFFF_FFFF_FFFF
FNV_OFFSET_BASIS = 0xCBF2_9CE4_8422_2325
FNV_PRIME = 0x0000_0100_0000_01B3
SPLITMIX_MUL0 = 0xBF58_476D_1CE4_E5B9
SPLITMIX_MUL1 = 0x94D0_49BB_1331_11EB

_BASES = "ACGT"


def _build_code_table() -> Tuple[int, ...]:
    table = [4] * 256
    for ch, code in ((b"A", 0), (b"a", 0), (b"C", 1), (b"c", 1), (b"G", 2), (b"g", 2), (b"T", 3), (b"t", 3)):
        table[ch[0]] = code
    return tuple(table)


CODE: Tuple[int, ...] = _build_code_table()


def undo_xorshift(y: int, shift: int, bits: int) -> int:
    """Exact inverse of y = x ^ (x >> shift) on a `bits`-wide value (src/dna.rs:152)."""
    x = y
    covered = shift
    while covered < bits:
        x = y ^ (x >> shift)
        covered += shift
    return x


def inverse_mod_2_64(a: int) -> int:
    """Multiplicative inverse of an odd u64 modulo 2^64 (src/dna.rs:164)."""
    if (a & 1) != 1:
        raise ValueError(f"expected odd integer, got {a}")
    a &= MASK64
    inv = a
    for _ in range(6):
        inv = (inv * ((2 - ((a * inv) & MASK64)) & MASK64)) & MASK64
    return inv


class KmerSpace:
    """Exact mirror of src/dna.rs::KmerSpace for k in 5..=32."""

    __slots__ = ("k", "bits", "mask", "shifts", "mul", "mul_inv")

    def __init__(self, k: int) -> None:
        if not (MIN_K <= k <= MAX_K):
            raise ValueError(f"k must be within {MIN_K}..={MAX_K}, got {k}")
        self.k = k
        self.bits = 2 * k
        self.mask = MASK64 if self.bits == 64 else (1 << self.bits) - 1

        def scale(s: int) -> int:
            return max(1, (s * self.bits) // 64)

        self.shifts = (scale(30), scale(27), scale(31))
        m0 = (SPLITMIX_MUL0 & self.mask) | 1
        m1 = (SPLITMIX_MUL1 & self.mask) | 1
        self.mul = (m0, m1)
        self.mul_inv = (inverse_mod_2_64(m0) & self.mask, inverse_mod_2_64(m1) & self.mask)

    def mix(self, x: int) -> int:
        s1, s2, s3 = self.shifts
        mask = self.mask
        m0, m1 = self.mul
        x &= mask
        x ^= x >> s1
        x = (x * m0) & mask
        x ^= x >> s2
        x = (x * m1) & mask
        return (x ^ (x >> s3)) & mask

    def unmix(self, y: int) -> int:
        s1, s2, s3 = self.shifts
        mask = self.mask
        mi0, mi1 = self.mul_inv
        y &= mask
        y = undo_xorshift(y, s3, self.bits)
        y = (y * mi1) & mask
        y = undo_xorshift(y, s2, self.bits)
        y = (y * mi0) & mask
        return undo_xorshift(y, s1, self.bits) & mask

    def iter_canonical(self, seq: Union[bytes, bytearray, memoryview, str]) -> Iterator[int]:
        """Yields every canonical 2-bit packed k-mer in `seq`; non-ACGT bytes reset the window."""
        raw = seq.encode("ascii") if isinstance(seq, str) else seq
        k = self.k
        shift = 2 * (k - 1)
        mask = self.mask
        code_table = CODE
        fwd = 0
        rev = 0
        filled = 0
        for byte in raw:
            c = code_table[byte]
            if c > 3:
                filled = 0
                continue
            fwd = ((fwd << 2) | c) & mask
            rev = (rev >> 2) | ((3 - c) << shift)
            filled += 1
            if filled >= k:
                yield fwd if fwd < rev else rev

    def canonical(self, text: str) -> Optional[int]:
        if len(text) != self.k:
            return None
        out: Optional[int] = None
        for val in self.iter_canonical(text.encode("ascii", errors="replace")):
            out = val
        return out

    def reverse_complement(self, kmer: int) -> int:
        """Scalar reverse complement matching src/dna.rs:133."""
        out = 0
        x = kmer & self.mask
        for _ in range(self.k):
            out = (out << 2) | (3 - (x & 3))
            x >>= 2
        return out

    def reverse_complement_swar(self, kmer: int) -> int:
        """O(1) SWAR bit-parallel reverse complement (Tier-1 optimization reference)."""
        c = (~(kmer & self.mask)) & MASK64
        c = int.from_bytes(c.to_bytes(8, "big"), "little")
        c = ((c & 0x0F0F_0F0F_0F0F_0F0F) << 4) | ((c >> 4) & 0x0F0F_0F0F_0F0F_0F0F)
        c = ((c & 0x3333_3333_3333_3333) << 2) | ((c >> 2) & 0x3333_3333_3333_3333)
        return (c >> (64 - self.bits)) & self.mask

    def decode(self, kmer: int) -> str:
        return "".join(_BASES[(kmer >> (2 * i)) & 3] for i in range(self.k - 1, -1, -1))


class Fnv1a64:
    """Exact mirror of src/fnv.rs::Fnv."""

    __slots__ = ("state",)

    def __init__(self, state: int = FNV_OFFSET_BASIS) -> None:
        self.state = state & MASK64

    def update_bytes(self, data: Union[bytes, bytearray, memoryview]) -> None:
        h = self.state
        for b in data:
            h = ((h ^ b) * FNV_PRIME) & MASK64
        self.state = h

    def update_u64(self, value: int) -> None:
        self.update_bytes(struct.pack("<Q", value & MASK64))

    def update_str(self, text: str) -> None:
        raw = text.encode("utf-8")
        self.update_u64(len(raw))
        self.update_bytes(raw)

    def hex(self) -> str:
        return f"{self.state:016x}"


def kset_digest(sorted_mixed_values: Sequence[int]) -> str:
    """FNV-1a64 digest of sorted distinct mixed k-mer values (src/kset.rs:265)."""
    h = Fnv1a64()
    for v in sorted_mixed_values:
        h.update_u64(v)
    return h.hex()


def sha256_u64_stream(values: Iterable[int]) -> str:
    h = hashlib.sha256()
    buf = bytearray()
    for v in values:
        buf.extend(struct.pack("<Q", v & MASK64))
        if len(buf) >= 65536:
            h.update(buf)
            buf.clear()
    if buf:
        h.update(buf)
    return h.hexdigest()


def sha256_count_spectrum(counts: Mapping[int, int]) -> str:
    """Deterministic SHA-256 over sorted (canonical_u64, count_u64) pairs."""
    h = hashlib.sha256()
    buf = bytearray()
    for kmer in sorted(counts):
        buf.extend(struct.pack("<QQ", kmer & MASK64, counts[kmer] & MASK64))
        if len(buf) >= 65536:
            h.update(buf)
            buf.clear()
    if buf:
        h.update(buf)
    return h.hexdigest()


def _bucket_of(value: int, bits: int) -> int:
    return value >> (bits - BUCKET_BITS)


def _bucket_floor(bucket: int, bits: int) -> int:
    return bucket << (bits - BUCKET_BITS)


def _put_varint(out: bytearray, v: int) -> None:
    while v >= 0x80:
        out.append((v & 0x7F) | 0x80)
        v >>= 7
    out.append(v)


def _get_varint(data: bytes, pos: int) -> Tuple[int, int]:
    v = 0
    shift = 0
    n = len(data)
    while True:
        if pos >= n:
            raise ValueError("truncated varint")
        byte = data[pos]
        pos += 1
        v |= (byte & 0x7F) << shift
        if byte < 0x80:
            return v, pos
        shift += 7
        if shift >= 64:
            raise ValueError("varint overflow")


@dataclass(frozen=True)
class KsetData:
    k: int
    name: str
    values: Tuple[int, ...]
    digest: str


def write_kset(path: pathlib.Path, name: str, space: KmerSpace, sorted_mixed: Sequence[int]) -> int:
    """Writes a byte-compatible GRM .kset v1 file (src/kset.rs:154-246). Returns bytes written."""
    for idx in range(1, len(sorted_mixed)):
        if sorted_mixed[idx - 1] >= sorted_mixed[idx]:
            raise ValueError("kset values must be strictly ascending")
    bits = space.bits
    data = bytearray()
    offsets = [0] * (BUCKETS + 1)
    counts = [0] * BUCKETS
    i = 0
    n = len(sorted_mixed)
    for b in range(BUCKETS):
        offsets[b] = len(data)
        prev = _bucket_floor(b, bits)
        while i < n and _bucket_of(sorted_mixed[i], bits) == b:
            _put_varint(data, sorted_mixed[i] - prev)
            prev = sorted_mixed[i]
            counts[b] += 1
            i += 1
    offsets[BUCKETS] = len(data)
    name_bytes = name.encode("utf-8")
    out = bytearray()
    out.extend(KSET_MAGIC)
    out.extend(struct.pack("<BBHQI", space.k, BUCKET_BITS, 0, n, len(name_bytes)))
    out.extend(name_bytes)
    out.extend(struct.pack(f"<{BUCKETS + 1}Q", *offsets))
    out.extend(struct.pack(f"<{BUCKETS}I", *counts))
    out.extend(data)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(out)
    return len(out)


def read_kset(path: pathlib.Path) -> KsetData:
    """Reads and validates a GRM .kset v1 file (src/kset.rs:285-392)."""
    raw = path.read_bytes()
    if len(raw) < 24:
        raise ValueError(f"{path}: truncated kset header")
    if raw[:8] != KSET_MAGIC:
        raise ValueError(f"{path}: invalid kset magic")
    k, bucket_bits, _reserved, n_kmers, name_len = struct.unpack_from("<BBHQI", raw, 8)
    if bucket_bits != BUCKET_BITS:
        raise ValueError(f"{path}: unsupported bucket width {bucket_bits}")
    if not (MIN_K <= k <= MAX_K):
        raise ValueError(f"{path}: unsupported k {k}")
    table_len = (BUCKETS + 1) * 8 + BUCKETS * 4
    data_start = 24 + name_len + table_len
    if data_start > len(raw):
        raise ValueError(f"{path}: truncated kset table")
    name = raw[24 : 24 + name_len].decode("utf-8")
    off_start = 24 + name_len
    offsets = struct.unpack_from(f"<{BUCKETS + 1}Q", raw, off_start)
    counts = struct.unpack_from(f"<{BUCKETS}I", raw, off_start + (BUCKETS + 1) * 8)
    data_len = len(raw) - data_start
    if offsets[0] != 0 or offsets[BUCKETS] != data_len or sum(counts) != n_kmers:
        raise ValueError(f"{path}: inconsistent bucket table")
    bits = 2 * k
    values: List[int] = []
    for b in range(BUCKETS):
        start, end = offsets[b], offsets[b + 1]
        if end < start or counts[b] > (end - start):
            raise ValueError(f"{path}: inconsistent bucket {b} bounds")
        chunk = raw[data_start + start : data_start + end]
        prev = _bucket_floor(b, bits)
        pos = 0
        for j in range(counts[b]):
            gap, pos = _get_varint(chunk, pos)
            if j > 0 and gap == 0:
                raise ValueError(f"{path}: bucket {b} repeats a k-mer")
            prev += gap
            if prev > MASK64:
                raise ValueError(f"{path}: bucket {b} overflows u64")
            values.append(prev)
        if pos != len(chunk) or (counts[b] > 0 and _bucket_of(prev, bits) != b):
            raise ValueError(f"{path}: bucket {b} does not match count and range")
    return KsetData(k=k, name=name, values=tuple(values), digest=kset_digest(values))


def _is_gzip(path: pathlib.Path) -> bool:
    with path.open("rb") as f:
        return f.read(2) == b"\x1f\x8b"


def iter_sequence_records(path: pathlib.Path) -> Iterator[Tuple[str, bytes]]:
    """Streaming FASTA / FASTQ reader matching src/fasta.rs:28-72 (gzip-aware)."""
    opener = gzip.open if _is_gzip(path) else open
    with opener(path, "rb") as fh:
        name = ""
        seq = bytearray()
        in_record = False
        while True:
            line = fh.readline()
            if not line:
                break
            line = line.rstrip(b"\r\n")
            if not line:
                continue
            first = line[0]
            if first == 62:  # '>'
                if in_record:
                    yield name, bytes(seq)
                header = line[1:].decode("utf-8", errors="replace").split()
                name = header[0] if header else ""
                seq.clear()
                in_record = True
            elif first == 64 and (not in_record or not seq):  # '@'
                header = line[1:].decode("utf-8", errors="replace").split()
                name = header[0] if header else ""
                seq_line = fh.readline().rstrip(b"\r\n")
                yield name, seq_line
                fh.readline()  # '+'
                fh.readline()  # quality
                in_record = False
            elif in_record:
                seq.extend(line)
            else:
                raise ValueError(f"{path}: sequence data before any header")
        if in_record:
            yield name, bytes(seq)


@dataclass(frozen=True)
class CountResult:
    k: int
    min_abundance: int
    total_kmers: int
    distinct_all: int
    distinct_kept: int
    counts: Dict[int, int]
    mixed_sorted: Tuple[int, ...]
    fnv1a64_digest: str
    sha256_canonical: str
    sha256_mixed: str
    sha256_spectrum: str
    histogram: Dict[int, int]


def count_sequences(
    records: Iterable[Tuple[str, bytes]],
    k: int,
    min_abundance: int = 1,
) -> CountResult:
    """Exact canonical k-mer counter across an iterable of (name, sequence_bytes)."""
    if min_abundance < 1:
        raise ValueError(f"min_abundance must be >= 1, got {min_abundance}")
    space = KmerSpace(k)
    raw_counts: Counter[int] = Counter()
    total = 0
    for _name, seq in records:
        for canonical in space.iter_canonical(seq):
            raw_counts[canonical] += 1
            total += 1
    kept = {km: cnt for km, cnt in raw_counts.items() if cnt >= min_abundance}
    hist: Counter[int] = Counter(kept.values())
    mixed_sorted = tuple(sorted(space.mix(km) for km in kept))
    canonical_sorted = sorted(kept)
    return CountResult(
        k=k,
        min_abundance=min_abundance,
        total_kmers=total,
        distinct_all=len(raw_counts),
        distinct_kept=len(kept),
        counts=kept,
        mixed_sorted=mixed_sorted,
        fnv1a64_digest=kset_digest(mixed_sorted),
        sha256_canonical=sha256_u64_stream(canonical_sorted),
        sha256_mixed=sha256_u64_stream(mixed_sorted),
        sha256_spectrum=sha256_count_spectrum(kept),
        histogram=dict(sorted(hist.items())),
    )


def count_files(
    paths: Sequence[pathlib.Path],
    k: int,
    min_abundance: int = 1,
) -> CountResult:
    """Exact canonical k-mer counter across one or more FASTA/FASTQ(.gz) files."""
    def _all_records() -> Iterator[Tuple[str, bytes]]:
        for p in paths:
            yield from iter_sequence_records(p)

    return count_sequences(_all_records(), k=k, min_abundance=min_abundance)


def count_file(path: pathlib.Path, k: int, min_abundance: int = 1) -> CountResult:
    return count_files([path], k=k, min_abundance=min_abundance)


def parse_kmer_text_dump(
    path: pathlib.Path,
    k: int,
    min_abundance: int = 1,
) -> Dict[int, int]:
    """Parses whitespace/TSV '<kmer> <count>' or FASTA '><count>\\n<kmer>' dumps."""
    space = KmerSpace(k)
    counts: Dict[int, int] = {}
    pending_count: Optional[int] = None
    opener = gzip.open if _is_gzip(path) else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for line_no, raw_line in enumerate(fh, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                pending_count = int(line[1:].split()[0])
                continue
            if pending_count is not None:
                km = space.canonical(line)
                if km is None:
                    raise ValueError(f"{path}:{line_no}: invalid {k}-mer {line!r}")
                if pending_count >= min_abundance:
                    counts[km] = counts.get(km, 0) + pending_count
                pending_count = None
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}:{line_no}: expected '<kmer> <count>'")
            km = space.canonical(parts[0])
            if km is None:
                raise ValueError(f"{path}:{line_no}: invalid {k}-mer {parts[0]!r}")
            cnt = int(parts[1])
            if cnt >= min_abundance:
                counts[km] = counts.get(km, 0) + cnt
    return counts


@dataclass(frozen=True)
class SetParityReport:
    ok: bool
    k: int
    expected_distinct: int
    observed_distinct: int
    missing_count: int
    extra_count: int
    expected_fnv1a64: str
    observed_fnv1a64: str
    expected_sha256: str
    observed_sha256: str


def verify_set_parity(
    expected_canonical: Iterable[int],
    observed_canonical: Iterable[int],
    k: int,
) -> SetParityReport:
    """Exact canonical k-mer set parity + FNV-1a64 (kset::digest) + SHA-256 verification."""
    space = KmerSpace(k)
    exp_set = set(expected_canonical)
    obs_set = set(observed_canonical)
    missing = len(exp_set - obs_set)
    extra = len(obs_set - exp_set)
    exp_mixed = sorted(space.mix(v) for v in exp_set)
    obs_mixed = sorted(space.mix(v) for v in obs_set)
    exp_fnv = kset_digest(exp_mixed)
    obs_fnv = kset_digest(obs_mixed)
    exp_sha = sha256_u64_stream(sorted(exp_set))
    obs_sha = sha256_u64_stream(sorted(obs_set))
    ok = (missing == 0) and (extra == 0) and (exp_fnv == obs_fnv) and (exp_sha == obs_sha)
    return SetParityReport(
        ok=ok,
        k=k,
        expected_distinct=len(exp_set),
        observed_distinct=len(obs_set),
        missing_count=missing,
        extra_count=extra,
        expected_fnv1a64=exp_fnv,
        observed_fnv1a64=obs_fnv,
        expected_sha256=exp_sha,
        observed_sha256=obs_sha,
    )


@dataclass(frozen=True)
class AbundanceParityReport:
    ok: bool
    k: int
    min_abundance: int
    max_cap: Optional[int]
    expected_distinct: int
    observed_distinct: int
    missing_count: int
    extra_count: int
    count_mismatch: int
    expected_spectrum_sha256: str
    observed_spectrum_sha256: str


def verify_abundance_parity(
    expected_counts: Mapping[int, int],
    observed_counts: Mapping[int, int],
    k: int,
    min_abundance: int = 1,
    max_cap: Optional[int] = None,
) -> AbundanceParityReport:
    """Exact canonical k-mer abundance parity (optionally capped, e.g. 255 for FastK u8)."""
    def _norm(mapping: Mapping[int, int]) -> Dict[int, int]:
        out: Dict[int, int] = {}
        for km, cnt in mapping.items():
            if cnt >= min_abundance:
                out[km] = min(cnt, max_cap) if max_cap is not None else cnt
        return out

    exp = _norm(expected_counts)
    obs = _norm(observed_counts)
    exp_keys = set(exp)
    obs_keys = set(obs)
    missing = len(exp_keys - obs_keys)
    extra = len(obs_keys - exp_keys)
    shared = exp_keys & obs_keys
    mismatched = sum(1 for km in shared if exp[km] != obs[km])
    exp_sha = sha256_count_spectrum(exp)
    obs_sha = sha256_count_spectrum(obs)
    ok = (missing == 0) and (extra == 0) and (mismatched == 0) and (exp_sha == obs_sha)
    return AbundanceParityReport(
        ok=ok,
        k=k,
        min_abundance=min_abundance,
        max_cap=max_cap,
        expected_distinct=len(exp),
        observed_distinct=len(obs),
        missing_count=missing,
        extra_count=extra,
        count_mismatch=mismatched,
        expected_spectrum_sha256=exp_sha,
        observed_spectrum_sha256=obs_sha,
    )


def compact_unitigs(canonical_kmers: Iterable[int], k: int) -> Tuple[frozenset, ...]:
    """Reference compacted de Bruijn graph matching src/input_build.rs:413-533.

    Returns each unitig as a frozenset of its canonical packed k-mers,
    deterministically ordered by smallest canonical k-mer.
    """
    space = KmerSpace(k)
    kmers = sorted(set(canonical_kmers))
    idx_of = {v: i for i, v in enumerate(kmers)}
    n = len(kmers)
    mask = space.mask
    shift = space.bits - 2

    def canon_of(oriented: int) -> int:
        rc = space.reverse_complement_swar(oriented)
        return oriented if oriented <= rc else rc

    def succ(oriented: int) -> List[Tuple[int, int]]:
        out: List[Tuple[int, int]] = []
        base = (oriented << 2) & mask
        for b in range(4):
            y = base | b
            i = idx_of.get(canon_of(y))
            if i is not None:
                out.append((i, y))
        return out

    def pred(oriented: int) -> List[Tuple[int, int]]:
        out: List[Tuple[int, int]] = []
        base = oriented >> 2
        for b in range(4):
            y = base | (b << shift)
            i = idx_of.get(canon_of(y))
            if i is not None:
                out.append((i, y))
        return out

    parent = list(range(n))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    for i, x in enumerate(kmers):
        s = succ(x)
        if len(s) == 1 and s[0][0] != i:
            j, y = s[0]
            if len(pred(y)) == 1:
                union(i, j)
        p = pred(x)
        if len(p) == 1 and p[0][0] != i:
            j, y = p[0]
            if len(succ(y)) == 1:
                union(i, j)

    groups: Dict[int, Set[int]] = {}
    for i, x in enumerate(kmers):
        groups.setdefault(find(i), set()).add(x)
    return tuple(frozenset(groups[r]) for r in sorted(groups))


def parse_unitigs_fasta(path: pathlib.Path, k: int) -> Tuple[frozenset, ...]:
    """Reads a unitigs FASTA file and returns each unitig as a frozenset of canonical k-mers."""
    space = KmerSpace(k)
    seen: Set[int] = set()
    unitigs: List[frozenset] = []
    for _name, seq in iter_sequence_records(path):
        u_kmers = frozenset(space.iter_canonical(seq))
        if not u_kmers:
            continue
        overlap = seen & u_kmers
        if overlap:
            raise ValueError(f"{path}: {len(overlap)} k-mers appear in multiple unitigs")
        seen.update(u_kmers)
        unitigs.append(u_kmers)
    unitigs.sort(key=lambda s: min(s))
    return tuple(unitigs)


@dataclass(frozen=True)
class UnitigParityReport:
    ok: bool
    kmers_expected: int
    kmers_observed: int
    kmers_both: int
    kmers_expected_only: int
    kmers_observed_only: int
    unitigs_expected: int
    unitigs_observed: int
    unitigs_both: int
    unitigs_expected_only: int
    unitigs_observed_only: int


def verify_unitig_parity(
    expected_unitigs: Sequence[frozenset],
    observed_unitigs: Sequence[frozenset],
) -> UnitigParityReport:
    """Exact set-of-sets unitig comparison matching src/input_build.rs:575-626."""
    exp_kmers: Set[int] = set().union(*expected_unitigs) if expected_unitigs else set()
    obs_kmers: Set[int] = set().union(*observed_unitigs) if observed_unitigs else set()
    exp_set = set(expected_unitigs)
    obs_set = set(observed_unitigs)
    u_both = len(exp_set & obs_set)
    k_both = len(exp_kmers & obs_kmers)
    k_exp_only = len(exp_kmers - obs_kmers)
    k_obs_only = len(obs_kmers - exp_kmers)
    u_exp_only = len(exp_set) - u_both
    u_obs_only = len(obs_set) - u_both
    ok = (k_exp_only == 0) and (k_obs_only == 0) and (u_exp_only == 0) and (u_obs_only == 0)
    return UnitigParityReport(
        ok=ok,
        kmers_expected=len(exp_kmers),
        kmers_observed=len(obs_kmers),
        kmers_both=k_both,
        kmers_expected_only=k_exp_only,
        kmers_observed_only=k_obs_only,
        unitigs_expected=len(exp_set),
        unitigs_observed=len(obs_set),
        unitigs_both=u_both,
        unitigs_expected_only=u_exp_only,
        unitigs_observed_only=u_obs_only,
    )


@dataclass(frozen=True)
class ReferenceMatrix:
    k: int
    n_genomes: int
    min_support: int
    kmers_distinct: int
    kmers_kept: int
    kmers_rare: int
    kmers_core: int
    n_patterns: int
    kmer_to_mask: Dict[int, int]
    pattern_masks: Tuple[int, ...]
    digest_sha256: str


def build_reference_matrix(
    genome_canonical_sets: Sequence[Iterable[int]],
    k: int,
    min_support: int = 2,
) -> ReferenceMatrix:
    """Exact reference presence matrix matching src/matrix.rs::build semantics."""
    space = KmerSpace(k)
    n = len(genome_canonical_sets)
    max_support = max(0, n - min_support)
    mixed_to_mask: Dict[int, int] = {}
    for g_idx, kmers in enumerate(genome_canonical_sets):
        bit = 1 << g_idx
        for km in set(kmers):
            mv = space.mix(km)
            mixed_to_mask[mv] = mixed_to_mask.get(mv, 0) | bit

    distinct = len(mixed_to_mask)
    rare = 0
    core = 0
    kept_map: Dict[int, int] = {}
    pattern_id: Dict[int, int] = {}
    patterns: List[int] = []
    # Match GRM order: ascending mixed k-mer value, patterns numbered by first appearance.
    for mv in sorted(mixed_to_mask):
        mask = mixed_to_mask[mv]
        sup = mask.bit_count()
        if sup < min_support:
            rare += 1
        elif sup > max_support:
            core += 1
        else:
            kept_map[space.unmix(mv)] = mask
            if mask not in pattern_id:
                pattern_id[mask] = len(patterns)
                patterns.append(mask)

    h = hashlib.sha256()
    for mv in sorted(mixed_to_mask):
        mask = mixed_to_mask[mv]
        sup = mask.bit_count()
        if min_support <= sup <= max_support:
            h.update(struct.pack("<Q", mv))
            h.update(mask.to_bytes((n + 7) // 8 or 1, "little"))

    return ReferenceMatrix(
        k=k,
        n_genomes=n,
        min_support=min_support,
        kmers_distinct=distinct,
        kmers_kept=len(kept_map),
        kmers_rare=rare,
        kmers_core=core,
        n_patterns=len(patterns),
        kmer_to_mask=kept_map,
        pattern_masks=tuple(patterns),
        digest_sha256=h.hexdigest(),
    )


@dataclass(frozen=True)
class MatrixParityReport:
    ok: bool
    expected_kept: int
    observed_kept: int
    expected_patterns: int
    observed_patterns: int
    missing_kmers: int
    extra_kmers: int
    mask_mismatch: int


def verify_matrix_parity(
    expected: ReferenceMatrix,
    observed_kmer_to_mask: Mapping[int, int],
) -> MatrixParityReport:
    """Verifies exact presence bitset per kept canonical k-mer and collapsed pattern count."""
    exp = expected.kmer_to_mask
    obs = dict(observed_kmer_to_mask)
    exp_keys = set(exp)
    obs_keys = set(obs)
    missing = len(exp_keys - obs_keys)
    extra = len(obs_keys - exp_keys)
    mismatched = sum(1 for km in (exp_keys & obs_keys) if exp[km] != obs[km])
    obs_patterns = len(set(obs.values()))
    ok = (
        missing == 0
        and extra == 0
        and mismatched == 0
        and obs_patterns == expected.n_patterns
    )
    return MatrixParityReport(
        ok=ok,
        expected_kept=len(exp),
        observed_kept=len(obs),
        expected_patterns=expected.n_patterns,
        observed_patterns=obs_patterns,
        missing_kmers=missing,
        extra_kmers=extra,
        mask_mismatch=mismatched,
    )
