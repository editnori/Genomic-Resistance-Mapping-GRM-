"""Benchmark runners, command builders, resource parsers, and grm-kmer-bench-v1 writer.

Supports workloads T1..T5 across grm, fastk, kmc3, jellyfish, ntcard, kmtricks,
ggcat, cuttlefish2, and bcalm2 using Python's standard library only.
"""

from dataclasses import asdict, dataclass
import gzip
import hashlib
import importlib.util
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import threading
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


def _load_sibling(name: str):
    mod_name = f"_grm_kmer_bench_{name}"
    if mod_name in sys.modules:
        return sys.modules[mod_name]
    path = pathlib.Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


parity = _load_sibling("parity")

BENCH_FORMAT = "grm-kmer-bench-v1"

# Workloads T1 through T5.
WORKLOADS: Dict[str, str] = {
    "T1": "Assembly per-genome canonical k-mer counting (FASTA, min_abundance=1)",
    "T2": "Raw paired-end FASTQ counting and abundance recovery (min_abundance>=2)",
    "T3": "Multi-genome presence matrix construction (min_support=2)",
    "T4": "Compacted de Bruijn graph unitig construction over species union",
    "T5": "Direct-from-FASTQ multi-genome matrix construction (min_abundance>=2)",
}

# Supported tool registry: binary name, version args, workloads, parity modes.
TOOL_REGISTRY: Dict[str, Dict[str, Any]] = {
    "grm": {
        "binary": "grm",
        "version_args": None,
        "workloads": ("T1", "T2", "T3", "T4", "T5"),
        "exact_set_parity": True,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": True,
        "supports_matrix": True,
    },
    "fastk": {
        "binary": "FastK",
        "version_args": ["-h"],
        "workloads": ("T1", "T2"),
        "exact_set_parity": True,
        "exact_abundance_parity": True,
        "abundance_cap": 255,
        "supports_unitigs": False,
        "supports_matrix": False,
    },
    "kmc3": {
        "binary": "kmc",
        "version_args": [],
        "workloads": ("T1", "T2"),
        "exact_set_parity": True,
        "exact_abundance_parity": True,
        "abundance_cap": 65535,
        "supports_unitigs": False,
        "supports_matrix": False,
    },
    "jellyfish": {
        "binary": "jellyfish",
        "version_args": ["--version"],
        "workloads": ("T1", "T2"),
        "exact_set_parity": True,
        "exact_abundance_parity": True,
        "abundance_cap": None,
        "supports_unitigs": False,
        "supports_matrix": False,
    },
    "ntcard": {
        "binary": "ntcard",
        "version_args": ["--version"],
        "workloads": ("T1", "T2"),
        "exact_set_parity": False,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": False,
        "supports_matrix": False,
    },
    "kmtricks": {
        "binary": "kmtricks",
        "version_args": ["--version"],
        "workloads": ("T3", "T5"),
        "exact_set_parity": True,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": False,
        "supports_matrix": True,
    },
    "ggcat": {
        "binary": "ggcat",
        "version_args": ["--version"],
        "workloads": ("T3", "T4", "T5"),
        "exact_set_parity": True,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": True,
        "supports_matrix": True,
    },
    "cuttlefish2": {
        "binary": "cuttlefish",
        "version_args": ["--version"],
        "workloads": ("T4",),
        "exact_set_parity": True,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": True,
        "supports_matrix": False,
    },
    "bcalm2": {
        "binary": "bcalm",
        "version_args": ["-version"],
        "workloads": ("T4",),
        "exact_set_parity": True,
        "exact_abundance_parity": False,
        "abundance_cap": None,
        "supports_unitigs": True,
        "supports_matrix": False,
    },
}

# Genome-id, accession, dash, and private-detail hygiene regexes (aligned with
# scripts/c7-outside/scan.py and tests/public_words/rules.py).
GENOME_ID_TOKEN = re.compile(r"(?<![A-Za-z0-9_.])[A-Za-z0-9_]+(?:\.[A-Za-z0-9_]+)?(?![A-Za-z0-9_.])")
BVBRC_ID_SHAPE = re.compile(r"(?<![\w.:/\-+])[1-9]\d{0,6}\.[1-9]\d{2,6}(?![\w.:/\-%])")
BVBRC_CONTEXT_SHORT = re.compile(
    r"(?i)\b(?:genome(?:_id)?|isolate|sample|strain|taxon|bvbrc|id)[ \t]*[:=#]?[ \t]*([1-9]\d{0,6}\.[1-9]\d{0,6})(?![\w.:/\-%])"
)
ASSEMBLY_ACC_RE = re.compile(r"(?<![A-Za-z0-9_])GC[AF]_\d+(?:\.\d+)?(?![A-Za-z0-9_])")
RUN_ACC_RE = re.compile(r"(?<![A-Za-z0-9_])[SED]RR\d{5,}(?![A-Za-z0-9_])")
BIOSAMPLE_RE = re.compile(r"(?<![A-Za-z0-9_])SAM[NED][A-Z]?\d{5,}(?![A-Za-z0-9_])")
BIOPROJECT_RE = re.compile(r"(?<![A-Za-z0-9_])PRJ[NED][A-Z]\d{4,}(?![A-Za-z0-9_])")

_NON_ID_LABEL_BEFORE_RE = re.compile(
    r"(?i)(?:^|[^A-Za-z0-9_])(?:"
    r"version|ver\.?|release|rel\.?|rev(?:ision)?|build|"
    r"python|rust|rustc|ubuntu|node|gcc|linux|kernel|glibc|"
    r"bcalm2?|fastk|kmc3?|jellyfish|ggcat|cuttlefish2?|kmtricks|ntcard|grm|"
    r"date|time|elapsed|wall(?:_seconds)?|user(?:_seconds)?|sys(?:_seconds)?|seconds|"
    r"rss|max_rss(?:_kb)?|cpu(?:_percent)?|bytes_per_kmer|bpk|b/k-?mer|"
    r"ratio|depth|chrom_median_depth|coverage|read_coverage|error_rate|snp_rate|rate|"
    r"score|accuracy|bal_acc|recall|precision|f1|auc|auroc|loss|log_loss|"
    r"p|pval|p_value|holm|median|mean|std|sd|iqr|min|max|avg|"
    r"cutoff|threshold|delta|diff|margin|speedup|throughput|memory|ram|disk|size|gc|n50|"
    r"plasmid_\d+|row|part|section|table|fig(?:ure)?|chapter|line|step|page|clause|"
    r"metric|val(?:ue)?|num(?:ber)?|total|count|bound|ceiling|floor"
    r")[ \t]*[\"']?[ \t]*[:=]?[ \t]*$"
)
_UNIT_AFTER_RE = re.compile(
    r"^[ \t]*(?:"
    r"ns|us|ms|secs?|seconds?|mins?|minutes?|hrs?|hours?|"
    r"bytes?|B/k-?mer|k-?mers(?:/s)?|"
    r"KiB|MiB|GiB|TiB|KB|MB|GB|TB|"
    r"kbp|Mbp(?:/s)?|Gbp|MB/s|GB/s|KB/s|bp|"
    r"fold|pct|s(?![A-Za-z])|B(?![A-Za-z])|x(?![A-Za-z])"
    r")(?![A-Za-z0-9_])"
)

EM_DASH_ESCAPED_RE = re.compile(r"\\u" + "2014", re.IGNORECASE)
EN_DASH_ESCAPED_RE = re.compile(r"\\u" + "2013", re.IGNORECASE)

_S = "/"
HOME_PATH_RE = re.compile(_S + "home" + _S + r"([A-Za-z0-9_.-]+)|" + _S + "Users" + _S + r"([A-Za-z0-9_.-]+)")
NEUTRAL_HOMES = frozenset({"user", "runner", "proj"})
WINDOWS_USER_RE = re.compile(r"[A-Za-z]:\\+Users\\+|" + _S + "mnt" + _S + r"[a-z]" + _S + "Users" + _S)
THREAD_LINK_RE = re.compile("t3" + "-" + r"[a-z]+:" + _S + _S)
BOX_ID_RE = re.compile(r"\bbx" + "_" + r"[a-z0-9]{8}\b")

_NUMERIC_TSV_COLUMNS = frozenset({
    "k",
    "min_abundance",
    "threads",
    "cohort_genomes",
    "wall_seconds",
    "user_seconds",
    "sys_seconds",
    "max_rss_kb",
    "peak_temp_bytes",
    "final_artifact_bytes",
    "bytes_per_kmer",
})
_PLAIN_NUMBER_RE = re.compile(r"^-?\d+(?:\.\d+)?$")


def _has_bvbrc_id_shape(text: str) -> bool:
    for m in BVBRC_ID_SHAPE.finditer(text):
        line_start = text.rfind("\n", 0, m.start()) + 1
        line_end = text.find("\n", m.end())
        if line_end < 0:
            line_end = len(text)
        before = text[line_start : m.start()]
        after = text[m.end() : line_end]
        if _NON_ID_LABEL_BEFORE_RE.search(before):
            continue
        if _UNIT_AFTER_RE.search(after):
            continue
        return True
    if BVBRC_CONTEXT_SHORT.search(text):
        return True
    return False


def scan_for_leaks(
    text: str,
    forbidden_ids: Optional[Iterable[str]] = None,
    check_bvbrc_shape: bool = True,
) -> List[str]:
    """Returns a list of leak violation tags found in `text` (never prints offending tokens)."""
    violations: List[str] = []
    if "\u2014" in text or EM_DASH_ESCAPED_RE.search(text):
        violations.append("em_dash")
    if "\u2013" in text or EN_DASH_ESCAPED_RE.search(text):
        violations.append("en_dash")
    if WINDOWS_USER_RE.search(text):
        violations.append("windows_user_path")
    if THREAD_LINK_RE.search(text):
        violations.append("thread_link")
    if BOX_ID_RE.search(text):
        violations.append("box_id")
    for m in HOME_PATH_RE.finditer(text):
        name = m.group(1) or m.group(2)
        if name and re.search(r"[A-Za-z]", name) and name not in NEUTRAL_HOMES:
            violations.append("private_home_path")
            break
    if ASSEMBLY_ACC_RE.search(text):
        violations.append("assembly_accession")
    if RUN_ACC_RE.search(text):
        violations.append("run_accession")
    if BIOSAMPLE_RE.search(text):
        violations.append("biosample_accession")
    if BIOPROJECT_RE.search(text):
        violations.append("bioproject_accession")
    if forbidden_ids is not None:
        id_set = {str(x) for x in forbidden_ids if str(x)}
        if id_set and any(m.group(0) in id_set for m in GENOME_ID_TOKEN.finditer(text)):
            violations.append("known_genome_id")
    if check_bvbrc_shape and _has_bvbrc_id_shape(text):
        violations.append("bvbrc_id_shape")
    return violations


def scan_json_for_leaks(
    obj: Any,
    *,
    forbidden_ids: Optional[Iterable[str]] = None,
    check_bvbrc_shape: bool = True,
) -> List[str]:
    """Scans all decoded string keys and values in a JSON-compatible object."""
    id_set = {str(x) for x in forbidden_ids if str(x)} if forbidden_ids is not None else None
    found: List[str] = []
    seen = set()

    def _visit(val: Any) -> None:
        if isinstance(val, str):
            for tag in scan_for_leaks(val, forbidden_ids=id_set, check_bvbrc_shape=check_bvbrc_shape):
                if tag not in seen:
                    seen.add(tag)
                    found.append(tag)
        elif isinstance(val, Mapping):
            for k, v in val.items():
                _visit(str(k))
                _visit(v)
        elif isinstance(val, (list, tuple)):
            for item in val:
                _visit(item)

    _visit(obj)
    return found


def scan_tsv_for_leaks(
    tsv_text: str,
    *,
    forbidden_ids: Optional[Iterable[str]] = None,
    check_bvbrc_shape: bool = True,
) -> List[str]:
    """Scans TSV text without false alarms on plain numeric metric columns."""
    id_set = {str(x) for x in forbidden_ids if str(x)} if forbidden_ids is not None else None
    found: List[str] = []
    seen = set()

    def _add(tags: Iterable[str]) -> None:
        for tag in tags:
            if tag not in seen:
                seen.add(tag)
                found.append(tag)

    # Always scan whole text for dashes, paths, accessions, and explicit forbidden IDs.
    _add(scan_for_leaks(tsv_text, forbidden_ids=id_set, check_bvbrc_shape=False))
    if not check_bvbrc_shape:
        return found

    lines = tsv_text.splitlines()
    if not lines:
        return found
    headers = lines[0].split("\t")
    _add(scan_for_leaks(lines[0], forbidden_ids=id_set, check_bvbrc_shape=True))
    for line in lines[1:]:
        cells = line.split("\t")
        for idx, cell in enumerate(cells):
            col = headers[idx] if idx < len(headers) else ""
            if col in _NUMERIC_TSV_COLUMNS and (not cell or _PLAIN_NUMBER_RE.fullmatch(cell)):
                continue
            _add(scan_for_leaks(cell, forbidden_ids=id_set, check_bvbrc_shape=True))
    return found


@dataclass(frozen=True)
class TimeMetrics:
    wall_seconds: float
    user_seconds: float
    sys_seconds: float
    max_rss_kb: int
    exit_status: int
    cpu_percent: Optional[float] = None
    minor_page_faults: Optional[int] = None
    major_page_faults: Optional[int] = None
    fs_inputs: Optional[int] = None
    fs_outputs: Optional[int] = None


def _parse_elapsed_seconds(raw: str) -> float:
    parts = raw.strip().split(":")
    if len(parts) == 2:
        minutes = int(parts[0])
        seconds = float(parts[1])
        return minutes * 60.0 + seconds
    if len(parts) == 3:
        hours = int(parts[0])
        minutes = int(parts[1])
        seconds = float(parts[2])
        return hours * 3600.0 + minutes * 60.0 + seconds
    raise ValueError(f"unrecognized elapsed time format: {raw!r}")


def parse_gnu_time_v(text: str) -> TimeMetrics:
    """Parses GNU `/usr/bin/time -v` verbose output into structured metrics."""
    user_s: Optional[float] = None
    sys_s: Optional[float] = None
    wall_s: Optional[float] = None
    max_rss_kb: Optional[int] = None
    exit_status: Optional[int] = None
    cpu_pct: Optional[float] = None
    minor_pf: Optional[int] = None
    major_pf: Optional[int] = None
    fs_in: Optional[int] = None
    fs_out: Optional[int] = None

    for raw_line in text.splitlines():
        line = raw_line.strip()
        if line.startswith("User time (seconds):"):
            user_s = float(line.split(":", 1)[1].strip())
        elif line.startswith("System time (seconds):"):
            sys_s = float(line.split(":", 1)[1].strip())
        elif line.startswith("Elapsed (wall clock) time"):
            idx = line.rfind("):")
            if idx < 0:
                raise ValueError(f"malformed elapsed line: {line!r}")
            wall_s = _parse_elapsed_seconds(line[idx + 2 :].strip())
        elif line.startswith("Maximum resident set size (kbytes):"):
            max_rss_kb = int(line.split(":", 1)[1].strip())
        elif line.startswith("Exit status:"):
            exit_status = int(line.split(":", 1)[1].strip())
        elif line.startswith("Percent of CPU this job got:"):
            val = line.split(":", 1)[1].strip().rstrip("%")
            if val and val != "?":
                cpu_pct = float(val)
        elif line.startswith("Minor (reclaiming a frame) page faults:"):
            minor_pf = int(line.split(":", 1)[1].strip())
        elif line.startswith("Major (requiring I/O) page faults:"):
            major_pf = int(line.split(":", 1)[1].strip())
        elif line.startswith("File system inputs:"):
            fs_in = int(line.split(":", 1)[1].strip())
        elif line.startswith("File system outputs:"):
            fs_out = int(line.split(":", 1)[1].strip())

    missing = [
        name
        for name, val in (
            ("user_seconds", user_s),
            ("sys_seconds", sys_s),
            ("wall_seconds", wall_s),
            ("max_rss_kb", max_rss_kb),
        )
        if val is None
    ]
    if missing:
        raise ValueError(f"GNU time -v output missing required fields: {', '.join(missing)}")
    if wall_s is None or user_s is None or sys_s is None or max_rss_kb is None:
        raise AssertionError("unreachable")
    return TimeMetrics(
        wall_seconds=round(wall_s, 6),
        user_seconds=round(user_s, 6),
        sys_seconds=round(sys_s, 6),
        max_rss_kb=max_rss_kb,
        exit_status=0 if exit_status is None else exit_status,
        cpu_percent=cpu_pct,
        minor_page_faults=minor_pf,
        major_page_faults=major_pf,
        fs_inputs=fs_in,
        fs_outputs=fs_out,
    )


def directory_bytes(paths: Sequence[pathlib.Path]) -> int:
    """Computes total regular-file bytes across files and directory trees."""
    total = 0
    seen_inodes = set()
    for root in paths:
        if not root.exists():
            continue
        if root.is_file() and not root.is_symlink():
            try:
                st = root.stat()
                key = (st.st_dev, st.st_ino)
                if key not in seen_inodes:
                    seen_inodes.add(key)
                    total += st.st_size
            except OSError:
                pass
            continue
        for dirpath, _dirnames, filenames in os.walk(root):
            for name in filenames:
                fp = pathlib.Path(dirpath) / name
                if fp.is_symlink():
                    continue
                try:
                    st = fp.stat()
                    key = (st.st_dev, st.st_ino)
                    if key not in seen_inodes:
                        seen_inodes.add(key)
                        total += st.st_size
                except OSError:
                    pass
    return total


@dataclass(frozen=True)
class DiskMetrics:
    peak_temp_bytes: int
    final_artifact_bytes: int
    peak_total_bytes: int
    kmers_counted: int
    bytes_per_kmer: Optional[float]


def compute_disk_metrics(
    peak_temp_bytes: int,
    final_artifact_bytes: int,
    kmers_counted: int,
    peak_total_bytes: Optional[int] = None,
) -> DiskMetrics:
    if peak_temp_bytes < 0 or final_artifact_bytes < 0 or kmers_counted < 0:
        raise ValueError("disk metrics and k-mer counts must be non-negative")
    peak_total = max(peak_temp_bytes + final_artifact_bytes, peak_total_bytes or 0)
    bpk = round(final_artifact_bytes / kmers_counted, 6) if kmers_counted > 0 else None
    return DiskMetrics(
        peak_temp_bytes=peak_temp_bytes,
        final_artifact_bytes=final_artifact_bytes,
        peak_total_bytes=peak_total,
        kmers_counted=kmers_counted,
        bytes_per_kmer=bpk,
    )


class DiskSampler:
    """Background disk-usage sampler for temporary and output directories."""

    def __init__(
        self,
        temp_paths: Sequence[pathlib.Path],
        artifact_paths: Sequence[pathlib.Path],
        poll_seconds: float = 0.05,
    ) -> None:
        self._temp = tuple(temp_paths)
        self._artifacts = tuple(artifact_paths)
        self._poll = poll_seconds
        self._stop = threading.Event()
        self._thread: Optional[threading.Thread] = None
        self.peak_temp_bytes = 0
        self.peak_total_bytes = 0

    def sample_once(self) -> None:
        t_bytes = directory_bytes(self._temp)
        a_bytes = directory_bytes(self._artifacts)
        if t_bytes > self.peak_temp_bytes:
            self.peak_temp_bytes = t_bytes
        if (t_bytes + a_bytes) > self.peak_total_bytes:
            self.peak_total_bytes = t_bytes + a_bytes

    def _run(self) -> None:
        while not self._stop.is_set():
            self.sample_once()
            self._stop.wait(self._poll)
        self.sample_once()

    def start(self) -> None:
        self.sample_once()
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        if self._thread is not None:
            self._thread.join(timeout=5.0)
        self.sample_once()


@dataclass(frozen=True)
class ToolProvenance:
    tool: str
    binary_name: str
    installed: bool
    binary_sha256: Optional[str]
    version_line: Optional[str]
    git_revision: Optional[str] = None


def _sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _probe_git_revision(search_dirs: Sequence[pathlib.Path]) -> Optional[str]:
    for d in search_dirs:
        if not d.exists():
            continue
        try:
            proc = subprocess.run(
                ["git", "-C", str(d), "rev-parse", "--short=12", "HEAD"],
                capture_output=True,
                text=True,
                timeout=2,
                check=False,
            )
            if proc.returncode == 0:
                rev = proc.stdout.strip()
                if re.fullmatch(r"[0-9a-f]{7,40}", rev):
                    return rev
        except (OSError, subprocess.SubprocessError):
            continue
    return None


def capture_tool_provenance(
    tool: str,
    binary_override: Optional[str] = None,
) -> ToolProvenance:
    """Captures binary availability, SHA-256, and version fingerprint without leaking user paths."""
    if tool not in TOOL_REGISTRY:
        raise ValueError(f"unknown tool: {tool!r}")
    spec = TOOL_REGISTRY[tool]
    bin_name = binary_override or spec["binary"]
    resolved = shutil.which(bin_name)
    if resolved is None and binary_override:
        cand = pathlib.Path(binary_override)
        if cand.is_file() and os.access(cand, os.X_OK):
            resolved = str(cand)
    if resolved is None:
        return ToolProvenance(
            tool=tool,
            binary_name=pathlib.Path(bin_name).name,
            installed=False,
            binary_sha256=None,
            version_line=None,
            git_revision=None,
        )
    bin_path = pathlib.Path(resolved).resolve()
    sha = _sha256_file(bin_path) if bin_path.is_file() else None
    git_rev: Optional[str] = None
    version_line: Optional[str] = None

    if tool == "grm":
        # grm does not take --version; fingerprint with binary SHA-256 and git revision.
        repo_root = pathlib.Path(__file__).resolve().parents[2]
        git_rev = _probe_git_revision([bin_path.parent, repo_root])
        parts = [bin_path.name]
        if sha:
            parts.append(f"sha256:{sha[:16]}")
        if git_rev:
            parts.append(f"git:{git_rev}")
        version_line = " ".join(parts)
    else:
        vargs = list(spec["version_args"] or [])
        try:
            proc = subprocess.run(
                [str(bin_path), *vargs],
                capture_output=True,
                text=True,
                timeout=5,
                check=False,
            )
            combined = (proc.stdout + "\n" + proc.stderr).splitlines()
            for ln in combined:
                cleaned = ln.strip()
                if cleaned:
                    if not scan_for_leaks(cleaned, check_bvbrc_shape=False):
                        version_line = cleaned[:160]
                    else:
                        version_line = f"{bin_path.name} (version output redacted)"
                    break
        except (OSError, subprocess.SubprocessError):
            version_line = None

    return ToolProvenance(
        tool=tool,
        binary_name=bin_path.name,
        installed=True,
        binary_sha256=sha,
        version_line=version_line,
        git_revision=git_rev,
    )


@dataclass(frozen=True)
class CommandPlan:
    tool: str
    workload: str
    variant: str
    k: int
    min_abundance: int
    min_support: int
    threads: int
    mem_gb: int
    binary: str
    steps: Tuple[Tuple[str, ...], ...]
    work_dir: str
    grm_data_dir: Optional[str]
    species_slug: str
    artifact_paths: Tuple[str, ...]
    temp_dirs: Tuple[str, ...]
    output_scope: str
    comparable: bool
    not_comparable_reason: Optional[str]
    supports_set_parity: bool
    supports_abundance_parity: bool
    abundance_cap: Optional[int]
    supports_unitig_parity: bool
    supports_matrix_parity: bool
    notes: str


def _resolve_sample_paths(
    cohort_dir: pathlib.Path,
    cohort_manifest: Optional[Mapping[str, Any]],
) -> Tuple[pathlib.Path, pathlib.Path, pathlib.Path]:
    """Resolves genome 0 FASTA and paired FASTQ filenames from manifest or cohort_dir."""
    manifest = cohort_manifest
    if manifest is None:
        mp = cohort_dir / "manifest.json"
        if mp.is_file():
            try:
                manifest = json.loads(mp.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                manifest = None

    if manifest and manifest.get("genomes"):
        g0 = manifest["genomes"][0]
        sample_fasta = (cohort_dir / g0["fasta"]).resolve()
        if "reads_r1" in g0 and "reads_r2" in g0:
            sample_r1 = (cohort_dir / g0["reads_r1"]).resolve()
            sample_r2 = (cohort_dir / g0["reads_r2"]).resolve()
            return sample_fasta, sample_r1, sample_r2
        reads_dir = cohort_dir / "reads"
        gz1 = reads_dir / f"{g0['genome']}_R1.fastq.gz"
        gz2 = reads_dir / f"{g0['genome']}_R2.fastq.gz"
        if gz1.exists() and gz2.exists():
            return sample_fasta, gz1.resolve(), gz2.resolve()
        return (
            sample_fasta,
            (reads_dir / f"{g0['genome']}_R1.fastq").resolve(),
            (reads_dir / f"{g0['genome']}_R2.fastq").resolve(),
        )

    sample_fasta = (cohort_dir / "assemblies" / "syn_g0000.fna").resolve()
    gz1 = cohort_dir / "reads" / "syn_g0000_R1.fastq.gz"
    gz2 = cohort_dir / "reads" / "syn_g0000_R2.fastq.gz"
    if gz1.exists() and gz2.exists():
        return sample_fasta, gz1.resolve(), gz2.resolve()
    return (
        sample_fasta,
        (cohort_dir / "reads" / "syn_g0000_R1.fastq").resolve(),
        (cohort_dir / "reads" / "syn_g0000_R2.fastq").resolve(),
    )


def build_command_plan(
    tool: str,
    workload: str,
    cohort_dir: pathlib.Path,
    run_dir: pathlib.Path,
    *,
    k: int = 31,
    min_abundance: Optional[int] = None,
    min_support: int = 2,
    threads: int = 8,
    mem_gb: int = 16,
    variant: str = "default",
    species_slug: str = "synthetic_sp",
    binary_override: Optional[str] = None,
    cohort_manifest: Optional[Mapping[str, Any]] = None,
) -> CommandPlan:
    """Builds exact argv command steps for one (tool, workload, variant) benchmark case."""
    if tool not in TOOL_REGISTRY:
        raise ValueError(f"unknown tool: {tool!r}")
    if workload not in WORKLOADS:
        raise ValueError(f"unknown workload: {workload!r}")
    meta = TOOL_REGISTRY[tool]
    if workload not in meta["workloads"]:
        raise ValueError(f"tool {tool!r} does not support workload {workload!r}")

    # Resolve paths so dash-leading directory names never read as CLI flags (C10)
    # and commands work from any caller working directory (C4).
    cohort_abs = pathlib.Path(cohort_dir).resolve()
    run_abs = pathlib.Path(run_dir).resolve()

    eff_min = (1 if workload in ("T1", "T3", "T4") else 2) if min_abundance is None else min_abundance
    bin_cmd = binary_override or meta["binary"]
    work = run_abs / f"{workload}_{tool}_{variant}"
    tmp = work / "tmp"
    out = work / "out"
    asm_fof = work / "assemblies.fof"
    reads_kmtricks_fof = work / "reads.kmtricks.fof"
    reads_plain_fof = work / "reads.fof"
    kmc_fof = work / "kmc_fof.txt"

    sample_fasta, sample_r1, sample_r2 = _resolve_sample_paths(cohort_abs, cohort_manifest)

    steps: List[Tuple[str, ...]] = []
    artifacts: List[str] = []
    notes = ""
    grm_data_dir: Optional[str] = None
    output_scope = "cohort"
    comparable = True
    not_comparable_reason: Optional[str] = None

    if tool == "grm":
        # Share one staged GRM data directory for T1, T3 (default/column-blocks), and T4 (C2).
        # Disk-sparse T3 deletes ksets, so it gets its own staged grm_data directory.
        grm_data = (work / "grm_data") if (workload == "T3" and variant == "disk-sparse-256") else (run_abs / "grm_data")
        grm_data_dir = str(grm_data)
        if workload == "T1":
            steps.append((bin_cmd, "--data", str(grm_data), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            artifacts.append(str(grm_data / "species" / species_slug / "ksets"))
            notes = "Counts all cohort assemblies into .kset v1 files with rayon parallelism."
        elif workload == "T2":
            comparable = False
            not_comparable_reason = (
                "stock grm count runs at min_abundance=1 on single files rather than paired reads at min_abundance=2"
            )
            steps.append((bin_cmd, "--data", str(grm_data), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            artifacts.append(str(grm_data / "species" / species_slug / "ksets"))
            notes = not_comparable_reason
        elif workload == "T3":
            cmd = [bin_cmd, "--data", str(grm_data), "matrix", "--species", species_slug, "--min-support", str(min_support)]
            if variant == "column-blocks-64":
                cmd.extend(["--column-blocks", "64"])
            elif variant == "disk-sparse-256":
                cmd.extend(["--disk-sparse", "--block-width", "256", "--allow-missing-catalog"])
            steps.append(tuple(cmd))
            artifacts.append(str(grm_data / "species" / species_slug / "matrix"))
            notes = f"GRM multi-genome pattern matrix build ({variant})."
        elif workload == "T4":
            steps.append((
                bin_cmd, "--data", str(grm_data), "exp-input-build", "unitigs",
                "--species", species_slug, "--graph", "grm", "--out", str(out / "unitigs"),
            ))
            artifacts.append(str(out / "unitigs"))
            notes = "GRM native union + prefix-indexed de Bruijn graph compaction."
        elif workload == "T5":
            comparable = False
            not_comparable_reason = (
                "stock grm count runs at min_abundance=1 on single files rather than paired reads at min_abundance=2"
            )
            steps.append((bin_cmd, "--data", str(grm_data), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            steps.append((bin_cmd, "--data", str(grm_data), "matrix", "--species", species_slug, "--min-support", str(min_support)))
            artifacts.append(str(grm_data / "species" / species_slug / "matrix"))
            notes = not_comparable_reason

    elif tool == "fastk":
        output_scope = "single_genome"
        prefix = out / "fastk_out"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, f"-k{k}", f"-t{eff_min}", f"-T{threads}", f"-M{mem_gb}", f"-P{tmp}", f"-N{prefix}", *inp))
        artifacts.append(f"{prefix}.ktab")
        notes = "FastK super-k-mer table; default .ktab saturates counts at 255."

    elif tool == "kmc3":
        output_scope = "single_genome"
        prefix = out / "kmc_db"
        fmt_flag = "-fm" if workload == "T1" else "-fq"
        inp = str(sample_fasta) if workload == "T1" else f"@{kmc_fof}"
        steps.append((bin_cmd, f"-k{k}", f"-ci{eff_min}", "-cs65535", f"-t{threads}", f"-m{mem_gb}", fmt_flag, inp, str(prefix), str(tmp)))
        artifacts.extend([f"{prefix}.kmc_pre", f"{prefix}.kmc_suf"])
        notes = "KMC3 two-stage disk-binned k-mer counter with u16 saturation cap."

    elif tool == "jellyfish":
        output_scope = "single_genome"
        jf_out = out / "mer_counts.jf"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, "count", "-m", str(k), "-s", "100M", "-t", str(threads), "-C", "-o", str(jf_out), *inp))
        artifacts.append(str(jf_out))
        notes = "Jellyfish 2 lock-free in-memory hash counter (-C canonical)."

    elif tool == "ntcard":
        output_scope = "single_genome"
        prefix = out / "ntcard_hist"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, f"-k{k}", f"-t{threads}", "-p", str(prefix), *inp))
        artifacts.append(f"{prefix}_k{k}.hist")
        notes = "ntCard streaming cardinality/histogram estimator (F0 and f1..f64 only, no exact set dump)."

    elif tool == "kmtricks":
        fof = asm_fof if workload == "T3" else reads_kmtricks_fof
        steps.append((
            bin_cmd, "pipeline", "--file", str(fof), "--run-dir", str(out / "kmtricks_run"),
            "--kmer-size", str(k), "--hard-min", str(eff_min), "--threads", str(threads),
            "--mode", "kmer:bin:hard",
        ))
        artifacts.append(str(out / "kmtricks_run"))
        notes = "kmtricks multi-sample k-mer presence matrix pipeline."

    elif tool == "ggcat":
        fof = asm_fof if workload in ("T3", "T4") else reads_plain_fof
        out_fa = out / ("ggcat.colored.fa" if workload in ("T3", "T5") else "ggcat.unitigs.fa")
        cmd = [bin_cmd, "build", "-k", str(k), "-j", str(threads), "-s", str(eff_min), "-t", str(tmp), "-l", str(fof), "-o", str(out_fa)]
        if workload in ("T3", "T5"):
            cmd.append("-c")
        steps.append(tuple(cmd))
        artifacts.append(str(out_fa))
        notes = "GGCAT parallel compacted/colored de Bruijn graph builder."

    elif tool == "cuttlefish2":
        prefix = out / "cf2"
        steps.append((
            bin_cmd, "build", "-l", str(asm_fof), "-k", str(k), "-t", str(threads),
            "--ref", "-o", str(prefix), "-w", str(tmp),
        ))
        artifacts.append(f"{prefix}.fa")
        notes = "Cuttlefish 2 reference-mode compacted de Bruijn graph builder."

    elif tool == "bcalm2":
        prefix = out / "bcalm_out"
        steps.append((
            bin_cmd, "-in", str(asm_fof), "-kmer-size", str(k),
            "-abundance-min", str(eff_min), "-nb-cores", str(threads), "-out", str(prefix),
        ))
        artifacts.append(f"{prefix}.unitigs.fa")
        notes = "bcalm2 minimizer-partitioned compacted de Bruijn graph builder."

    return CommandPlan(
        tool=tool,
        workload=workload,
        variant=variant,
        k=k,
        min_abundance=eff_min,
        min_support=min_support,
        threads=threads,
        mem_gb=mem_gb,
        binary=bin_cmd,
        steps=tuple(steps),
        work_dir=str(work),
        grm_data_dir=grm_data_dir,
        species_slug=species_slug,
        artifact_paths=tuple(artifacts),
        temp_dirs=(str(tmp),),
        output_scope=output_scope,
        comparable=comparable,
        not_comparable_reason=not_comparable_reason,
        supports_set_parity=bool(meta["exact_set_parity"]),
        supports_abundance_parity=bool(meta["exact_abundance_parity"]) and workload in ("T1", "T2"),
        abundance_cap=meta["abundance_cap"],
        supports_unitig_parity=bool(meta["supports_unitigs"]) and workload == "T4",
        supports_matrix_parity=bool(meta["supports_matrix"]) and workload in ("T3", "T5"),
        notes=notes,
    )


def stage_grm_data_dir(
    cohort_dir: pathlib.Path,
    cohort_manifest: Mapping[str, Any],
    grm_data_dir: pathlib.Path,
    species_slug: str,
) -> pathlib.Path:
    """Stages a synthetic cohort into GRM's exact data layout (<grm_data>/fasta + species/<slug>/genomes.txt)."""
    cohort_abs = pathlib.Path(cohort_dir).resolve()
    data_abs = pathlib.Path(grm_data_dir).resolve()
    fasta_dir = data_abs / "fasta"
    sp_dir = data_abs / "species" / species_slug
    fasta_dir.mkdir(parents=True, exist_ok=True)
    sp_dir.mkdir(parents=True, exist_ok=True)

    genomes = cohort_manifest["genomes"]
    genome_ids = [str(g["genome"]) for g in genomes]
    genomes_txt = sp_dir / "genomes.txt"
    want_txt = "\n".join(genome_ids) + "\n"
    if not genomes_txt.is_file() or genomes_txt.read_text(encoding="utf-8") != want_txt:
        genomes_txt.write_text(want_txt, encoding="utf-8")

    for g in genomes:
        gid = str(g["genome"])
        src = (cohort_abs / g["fasta"]).resolve()
        dst = fasta_dir / f"{gid}.fna.gz"
        if dst.is_file() and dst.stat().st_size > 0:
            continue
        raw = src.read_bytes()
        if raw[:2] == b"\x1f\x8b":
            dst.write_bytes(raw)
        else:
            with dst.open("wb") as raw_fh:
                with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, compresslevel=6, mtime=0) as gz:
                    gz.write(raw)
    return data_abs


def materialize_run_lists(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
) -> Dict[str, pathlib.Path]:
    """Writes resolved absolute paths into run-local FOF list files under plan.work_dir (C4)."""
    cohort_abs = pathlib.Path(cohort_dir).resolve()
    work_dir = pathlib.Path(plan.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    genomes = cohort_manifest["genomes"]
    written: Dict[str, pathlib.Path] = {}

    asm_paths = [str((cohort_abs / g["fasta"]).resolve()) for g in genomes]
    asm_fof = work_dir / "assemblies.fof"
    asm_fof.write_text("\n".join(asm_paths) + "\n", encoding="utf-8")
    written["assemblies.fof"] = asm_fof

    if genomes and "reads_r1" in genomes[0] and "reads_r2" in genomes[0]:
        kmtricks_lines = [
            f"{g['genome']} : {(cohort_abs / g['reads_r1']).resolve()} ; {(cohort_abs / g['reads_r2']).resolve()}"
            for g in genomes
        ]
        kmtricks_fof = work_dir / "reads.kmtricks.fof"
        kmtricks_fof.write_text("\n".join(kmtricks_lines) + "\n", encoding="utf-8")
        written["reads.kmtricks.fof"] = kmtricks_fof

        read_lines: List[str] = []
        for g in genomes:
            read_lines.append(str((cohort_abs / g["reads_r1"]).resolve()))
            read_lines.append(str((cohort_abs / g["reads_r2"]).resolve()))
        reads_fof = work_dir / "reads.fof"
        reads_fof.write_text("\n".join(read_lines) + "\n", encoding="utf-8")
        written["reads.fof"] = reads_fof

        kmc_fof = work_dir / "kmc_fof.txt"
        kmc_fof.write_text(
            f"{(cohort_abs / genomes[0]['reads_r1']).resolve()}\n"
            f"{(cohort_abs / genomes[0]['reads_r2']).resolve()}\n",
            encoding="utf-8",
        )
        written["kmc_fof.txt"] = kmc_fof

    return written


def _mock_timing(tool: str, workload: str, n_kmers: int, threads: int) -> TimeMetrics:
    """Deterministic synthetic /usr/bin/time -v telemetry for CI --mock self-check runs."""
    seed_hash = int(hashlib.sha256(f"{tool}:{workload}:{n_kmers}:{threads}".encode()).hexdigest()[:8], 16)
    base_wall = round(0.05 + (n_kmers / 1_000_000.0) * (1.0 + (seed_hash % 7) * 0.15), 6)
    user_s = round(base_wall * min(threads, 4) * 0.85, 6)
    sys_s = round(base_wall * 0.12, 6)
    rss_kb = 16384 + (seed_hash % 8192) + (n_kmers // 8)
    fake_gnu = (
        f"\tUser time (seconds): {user_s:.2f}\n"
        f"\tSystem time (seconds): {sys_s:.2f}\n"
        f"\tPercent of CPU this job got: {int(round(100 * (user_s + sys_s) / max(base_wall, 1e-6)))}%\n"
        f"\tElapsed (wall clock) time (h:mm:ss or m:ss): 0:{base_wall:05.2f}\n"
        f"\tMaximum resident set size (kbytes): {rss_kb}\n"
        f"\tMinor (reclaiming a frame) page faults: 4096\n"
        f"\tMajor (requiring I/O) page faults: 0\n"
        f"\tFile system inputs: 0\n"
        f"\tFile system outputs: 128\n"
        f"\tExit status: 0\n"
    )
    return parse_gnu_time_v(fake_gnu)


def _display_commands(
    steps: Sequence[Sequence[str]],
    cohort_dir: pathlib.Path,
    run_dir: pathlib.Path,
) -> List[List[str]]:
    """Normalizes absolute host paths to neutral placeholders for deterministic, leak-free reports."""
    c_str = str(cohort_dir)
    r_str = str(run_dir)
    c_res = str(cohort_dir.resolve())
    r_res = str(run_dir.resolve())
    replacements = sorted(
        [(c_str, "<cohort>"), (c_res, "<cohort>"), (r_str, "<work>"), (r_res, "<work>")],
        key=lambda kv: len(kv[0]),
        reverse=True,
    )
    out: List[List[str]] = []
    for step in steps:
        norm_step: List[str] = []
        for idx, token in enumerate(step):
            val = token
            for src, dst in replacements:
                if src and src in val:
                    val = val.replace(src, dst)
            val = val.replace("\\", "/")
            if idx == 0 and ("/" in val) and not val.startswith("<"):
                val = pathlib.Path(val).name
            norm_step.append(val)
        out.append(norm_step)
    return out


def _reference_kmer_count(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
) -> int:
    """Computes total stored reference k-mers across all genomes included in the tool's output (C5)."""
    genomes = cohort_manifest["genomes"]
    manifest_k = int(cohort_manifest["spec"]["k"])
    manifest_support = int(cohort_manifest["spec"].get("min_support", 2))

    if plan.workload == "T1":
        if plan.output_scope == "cohort":
            if plan.k == manifest_k:
                return sum(int(g["assembly_truth"]["distinct_kmers_min1"]) for g in genomes)
            return sum(
                parity.count_file((cohort_dir / g["fasta"]).resolve(), k=plan.k, min_abundance=1).distinct_kept
                for g in genomes
            )
        if plan.k == manifest_k:
            return int(genomes[0]["assembly_truth"]["distinct_kmers_min1"])
        return parity.count_file((cohort_dir / genomes[0]["fasta"]).resolve(), k=plan.k, min_abundance=1).distinct_kept

    if plan.workload == "T2":
        if "reads_truth" not in genomes[0]:
            return 0
        if plan.k == manifest_k and plan.min_abundance == 2:
            return int(genomes[0]["reads_truth"]["distinct_kmers_min2"])
        return parity.count_files(
            [(cohort_dir / genomes[0]["reads_r1"]).resolve(), (cohort_dir / genomes[0]["reads_r2"]).resolve()],
            k=plan.k,
            min_abundance=plan.min_abundance,
        ).distinct_kept

    if plan.workload == "T3":
        if plan.k == manifest_k and plan.min_support == manifest_support:
            return int(cohort_manifest["cohort_truth"]["assembly_matrix"]["kmers_kept"])
        sets = [set(parity.count_file((cohort_dir / g["fasta"]).resolve(), k=plan.k, min_abundance=1).counts) for g in genomes]
        return parity.build_reference_matrix(sets, k=plan.k, min_support=plan.min_support).kmers_kept

    if plan.workload == "T4":
        if plan.k == manifest_k and "unitigs" in cohort_manifest["cohort_truth"]:
            return int(cohort_manifest["cohort_truth"]["unitigs"]["graph_kmers"])
        union_kmers: set = set()
        for g in genomes:
            union_kmers.update(parity.count_file((cohort_dir / g["fasta"]).resolve(), k=plan.k, min_abundance=1).counts)
        return len(union_kmers)

    if plan.k == manifest_k and plan.min_abundance == 2 and plan.min_support == manifest_support:
        return int(cohort_manifest["cohort_truth"].get("reads_matrix_min2", {}).get("kmers_kept", 0))
    sets = [
        set(
            parity.count_files(
                [(cohort_dir / g["reads_r1"]).resolve(), (cohort_dir / g["reads_r2"]).resolve()],
                k=plan.k,
                min_abundance=plan.min_abundance,
            ).counts
        )
        for g in genomes
    ]
    return parity.build_reference_matrix(sets, k=plan.k, min_support=plan.min_support).kmers_kept


def _looks_like_text_kmer_dump(path: pathlib.Path) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    head = path.read_bytes()[:4096]
    if b"\x00" in head:
        return False
    try:
        text = head.decode("utf-8")
    except UnicodeDecodeError:
        return False
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            return line[1:].split()[0].isdigit()
        parts = line.split()
        return len(parts) >= 2 and all(ch in "ACGTacgt" for ch in parts[0]) and parts[1].isdigit()
    return False


def _find_or_dump_external_counts(
    plan: CommandPlan,
    work_dir: pathlib.Path,
) -> Tuple[Optional[pathlib.Path], Optional[str]]:
    """Finds a text/FASTA k-mer dump for external T1/T2 tools, or runs a companion dumper if available."""
    out_dir = work_dir / "out"
    candidates = [
        out_dir / "dump.tsv",
        out_dir / "dump.fa",
        out_dir / "mer_counts.dump",
        out_dir / "mer_counts.tsv",
        out_dir / "kmc_db.dump.tsv",
        out_dir / "fastk_out.dump.tsv",
    ]
    for cand in candidates:
        if cand.is_file() and cand.stat().st_size > 0 and _looks_like_text_kmer_dump(cand):
            return cand, None

    if plan.tool == "jellyfish":
        jf_path = pathlib.Path(plan.artifact_paths[0])
        if _looks_like_text_kmer_dump(jf_path):
            return jf_path, None
        dump_out = work_dir / "verify_jellyfish_dump.tsv"
        try:
            proc = subprocess.run(
                [plan.binary, "dump", "-c", str(jf_path)],
                capture_output=True,
                text=True,
                timeout=60,
                check=False,
            )
            if proc.returncode == 0 and proc.stdout.strip():
                dump_out.write_text(proc.stdout, encoding="utf-8")
                if _looks_like_text_kmer_dump(dump_out):
                    return dump_out, None
        except (OSError, subprocess.SubprocessError):
            pass
        return None, "binary Jellyfish database has no available text dump"

    if plan.tool == "kmc3":
        pre_path = pathlib.Path(plan.artifact_paths[0])
        if _looks_like_text_kmer_dump(pre_path):
            return pre_path, None
        kmc_tools = shutil.which("kmc_tools")
        if kmc_tools:
            prefix = str(out_dir / "kmc_db")
            dump_out = work_dir / "verify_kmc_dump.tsv"
            try:
                proc = subprocess.run(
                    [kmc_tools, "transform", prefix, "dump", str(dump_out)],
                    capture_output=True,
                    text=True,
                    timeout=60,
                    check=False,
                )
                if proc.returncode == 0 and dump_out.is_file() and _looks_like_text_kmer_dump(dump_out):
                    return dump_out, None
            except (OSError, subprocess.SubprocessError):
                pass
        return None, "binary KMC3 database has no available kmc_tools text dump"

    if plan.tool == "fastk":
        ktab_path = pathlib.Path(plan.artifact_paths[0])
        if _looks_like_text_kmer_dump(ktab_path):
            return ktab_path, None
        tabex = shutil.which("Tabex")
        if tabex:
            dump_out = work_dir / "verify_fastk_dump.tsv"
            try:
                proc = subprocess.run(
                    [tabex, str(ktab_path), "LIST"],
                    capture_output=True,
                    text=True,
                    timeout=60,
                    check=False,
                )
                if proc.returncode == 0 and proc.stdout.strip():
                    lines = []
                    for ln in proc.stdout.splitlines():
                        parts = ln.strip().split()
                        if len(parts) >= 3 and parts[1].startswith("="):
                            lines.append(f"{parts[0].rstrip(':')}\t{parts[2]}")
                    if lines:
                        dump_out.write_text("\n".join(lines) + "\n", encoding="utf-8")
                        if _looks_like_text_kmer_dump(dump_out):
                            return dump_out, None
            except (OSError, subprocess.SubprocessError):
                pass
        return None, "binary FastK .ktab has no available text dump"

    return None, "no parser for tool output format"


def _verify_real_output(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
    final_bytes: int,
) -> Tuple[str, Dict[str, Any], Optional[int]]:
    """Parses real tool outputs and runs parity.py checks (C1).

    Returns (status, parity_dict, observed_kmers_for_bpk).
    Never returns 'ok' without actually parsing and verifying non-empty outputs.
    """
    genomes = cohort_manifest["genomes"]
    work_dir = pathlib.Path(plan.work_dir).resolve()
    cohort_abs = pathlib.Path(cohort_dir).resolve()

    if final_bytes <= 0:
        return (
            "failed_empty_output",
            {
                "kind": "missing_or_empty_output",
                "check_type": "real_tool_output",
                "self_check": False,
                "ok": False,
                "checked": True,
                "reason": "required output files missing or empty (0 bytes)",
            },
            None,
        )

    if plan.tool == "grm":
        if plan.grm_data_dir is None:
            raise ValueError("grm plan missing grm_data_dir")
        grm_data = pathlib.Path(plan.grm_data_dir).resolve()
        sp_dir = grm_data / "species" / plan.species_slug

        if plan.workload == "T1":
            kset_dir = sp_dir / "ksets"
            space = parity.KmerSpace(plan.k)
            n_passed = 0
            total_exp = 0
            total_obs = 0
            first_failure: Optional[Dict[str, Any]] = None
            for g in genomes:
                gid = str(g["genome"])
                kset_path = kset_dir / f"{gid}.kset"
                if not kset_path.is_file() or kset_path.stat().st_size == 0:
                    return (
                        "failed_empty_output",
                        {
                            "kind": "set",
                            "check_type": "real_tool_output",
                            "self_check": False,
                            "ok": False,
                            "checked": True,
                            "reason": f"missing or empty .kset for {gid}",
                        },
                        None,
                    )
                try:
                    kset = parity.read_kset(kset_path)
                except ValueError as exc:
                    return (
                        "parity_failed",
                        {
                            "kind": "set",
                            "check_type": "real_tool_output",
                            "self_check": False,
                            "ok": False,
                            "checked": True,
                            "reason": f"corrupt .kset for {gid}: {exc}",
                        },
                        None,
                    )
                if kset.k != plan.k:
                    first_failure = {"genome": gid, "reason": f"kset k={kset.k} != expected k={plan.k}"}
                    break
                ref = parity.count_file((cohort_abs / g["fasta"]).resolve(), k=plan.k, min_abundance=1)
                obs_canonical = [space.unmix(v) for v in kset.values]
                rep = parity.verify_set_parity(ref.counts.keys(), obs_canonical, k=plan.k)
                total_exp += rep.expected_distinct
                total_obs += rep.observed_distinct
                if rep.ok:
                    n_passed += 1
                elif first_failure is None:
                    first_failure = {"genome": gid, "report": asdict(rep)}
            all_ok = (first_failure is None) and (n_passed == len(genomes))
            return (
                "ok" if all_ok else "parity_failed",
                {
                    "kind": "set",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": all_ok,
                    "checked": True,
                    "genomes_checked": len(genomes),
                    "genomes_passed": n_passed,
                    "total_expected_kmers": total_exp,
                    "total_observed_kmers": total_obs,
                    "first_failure": first_failure,
                },
                total_obs if total_obs > 0 else None,
            )

        if plan.workload == "T3":
            matrix_dir = sp_dir / "matrix"
            required_files = ("matrix.json", "patterns.bin", "support.bin", "kmer_offsets.bin", "kmers.bin")
            if not (matrix_dir / "matrix.json").is_file() or (matrix_dir / "matrix.json").stat().st_size == 0:
                return (
                    "failed_empty_output",
                    {
                        "kind": "matrix",
                        "check_type": "real_tool_output",
                        "self_check": False,
                        "ok": False,
                        "checked": True,
                        "reason": "missing or empty matrix.json",
                    },
                    None,
                )
            for fname in required_files:
                if not (matrix_dir / fname).is_file():
                    return (
                        "failed_empty_output",
                        {
                            "kind": "matrix",
                            "check_type": "real_tool_output",
                            "self_check": False,
                            "ok": False,
                            "checked": True,
                            "reason": f"missing matrix file {fname}",
                        },
                        None,
                    )
            try:
                obs_mat = parity.read_grm_matrix(matrix_dir)
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                return (
                    "parity_failed",
                    {
                        "kind": "matrix",
                        "check_type": "real_tool_output",
                        "self_check": False,
                        "ok": False,
                        "checked": True,
                        "reason": f"invalid GRM matrix output: {exc}",
                    },
                    None,
                )
            sets = [
                set(parity.count_file((cohort_abs / g["fasta"]).resolve(), k=plan.k, min_abundance=1).counts)
                for g in genomes
            ]
            ref_mat = parity.build_reference_matrix(sets, k=plan.k, min_support=plan.min_support)
            m_rep = parity.verify_matrix_parity(ref_mat, obs_mat.kmer_to_mask)
            want_names = tuple(str(g["genome"]) for g in genomes)
            meta_ok = (
                obs_mat.k == plan.k
                and obs_mat.genomes == want_names
                and obs_mat.min_support == plan.min_support
                and obs_mat.kmers_distinct == ref_mat.kmers_distinct
                and obs_mat.kmers_kept == ref_mat.kmers_kept
                and obs_mat.kmers_rare == ref_mat.kmers_rare
                and obs_mat.kmers_core == ref_mat.kmers_core
                and obs_mat.n_patterns == ref_mat.n_patterns
            )
            all_ok = m_rep.ok and meta_ok
            return (
                "ok" if all_ok else "parity_failed",
                {
                    "kind": "matrix",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": all_ok,
                    "checked": True,
                    "meta_ok": meta_ok,
                    "matrix": asdict(m_rep),
                },
                obs_mat.kmers_kept if obs_mat.kmers_kept > 0 else None,
            )

        if plan.workload == "T4":
            unitigs_dir = pathlib.Path(plan.artifact_paths[0])
            # Support direct unitigs FASTA if present (e.g. stub or future exporter).
            fasta_candidates = [unitigs_dir / "unitigs.fa", unitigs_dir / "grm.unitigs.fa"]
            found_fa = next((p for p in fasta_candidates if p.is_file() and p.stat().st_size > 0), None)
            sets = [
                set(parity.count_file((cohort_abs / g["fasta"]).resolve(), k=plan.k, min_abundance=1).counts)
                for g in genomes
            ]
            if found_fa is not None:
                union_kmers: set = set().union(*sets) if sets else set()
                ref_u = parity.compact_unitigs(union_kmers, k=plan.k)
                try:
                    obs_u = parity.parse_unitigs_fasta(found_fa, k=plan.k)
                except ValueError as exc:
                    return (
                        "parity_failed",
                        {
                            "kind": "unitigs",
                            "check_type": "real_tool_output",
                            "self_check": False,
                            "ok": False,
                            "checked": True,
                            "reason": f"invalid unitig FASTA: {exc}",
                        },
                        None,
                    )
                u_rep = parity.verify_unitig_parity(ref_u, obs_u)
                return (
                    "ok" if u_rep.ok else "parity_failed",
                    {
                        "kind": "unitigs",
                        "check_type": "real_tool_output",
                        "self_check": False,
                        "ok": u_rep.ok,
                        "checked": True,
                        "unitigs": asdict(u_rep),
                    },
                    u_rep.kmers_observed if u_rep.kmers_observed > 0 else None,
                )

            if not (unitigs_dir / "build.json").is_file() or not (unitigs_dir / "genomes.tsv").is_file():
                return (
                    "failed_empty_output",
                    {
                        "kind": "unitigs",
                        "check_type": "real_tool_output",
                        "self_check": False,
                        "ok": False,
                        "checked": True,
                        "reason": "missing build.json or genomes.tsv in grm unitigs output",
                    },
                    None,
                )
            try:
                g_rep = parity.verify_grm_unitig_build(
                    unitigs_dir,
                    sets,
                    [str(g["genome"]) for g in genomes],
                    k=plan.k,
                )
            except (OSError, ValueError, KeyError, json.JSONDecodeError) as exc:
                return (
                    "parity_failed",
                    {
                        "kind": "unitigs",
                        "check_type": "real_tool_output",
                        "self_check": False,
                        "ok": False,
                        "checked": True,
                        "reason": f"invalid grm unitigs build: {exc}",
                    },
                    None,
                )
            return (
                "ok" if g_rep.ok else "parity_failed",
                {
                    "kind": "unitigs",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": g_rep.ok,
                    "checked": True,
                    "unitigs_build": asdict(g_rep),
                },
                g_rep.graph_kmers_observed if g_rep.graph_kmers_observed > 0 else None,
            )

    # Unitig FASTA tools: bcalm2, cuttlefish2, ggcat (T4).
    if plan.workload == "T4" and plan.tool in ("bcalm2", "cuttlefish2", "ggcat"):
        out_fa = pathlib.Path(plan.artifact_paths[0])
        if not out_fa.is_file() or out_fa.stat().st_size == 0:
            return (
                "failed_empty_output",
                {
                    "kind": "unitigs",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": False,
                    "checked": True,
                    "reason": f"missing or empty unitig FASTA {out_fa.name}",
                },
                None,
            )
        union_kmers = set()
        for g in genomes:
            union_kmers.update(parity.count_file((cohort_abs / g["fasta"]).resolve(), k=plan.k, min_abundance=1).counts)
        ref_u = parity.compact_unitigs(union_kmers, k=plan.k)
        try:
            obs_u = parity.parse_unitigs_fasta(out_fa, k=plan.k)
        except ValueError as exc:
            return (
                "parity_failed",
                {
                    "kind": "unitigs",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": False,
                    "checked": True,
                    "reason": f"invalid unitig FASTA: {exc}",
                },
                None,
            )
        u_rep = parity.verify_unitig_parity(ref_u, obs_u)
        return (
            "ok" if u_rep.ok else "parity_failed",
            {
                "kind": "unitigs",
                "check_type": "real_tool_output",
                "self_check": False,
                "ok": u_rep.ok,
                "checked": True,
                "unitigs": asdict(u_rep),
            },
            u_rep.kmers_observed if u_rep.kmers_observed > 0 else None,
        )

    # External single-genome k-mer counters: fastk, kmc3, jellyfish (T1, T2).
    if plan.tool in ("fastk", "kmc3", "jellyfish") and plan.workload in ("T1", "T2"):
        missing_files = [
            pathlib.Path(a).name for a in plan.artifact_paths if not pathlib.Path(a).is_file() or pathlib.Path(a).stat().st_size == 0
        ]
        if missing_files:
            return (
                "failed_empty_output",
                {
                    "kind": "set_and_abundance",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": False,
                    "checked": True,
                    "reason": f"missing or empty output file(s): {', '.join(missing_files)}",
                },
                None,
            )
        dump_path, unparsed_reason = _find_or_dump_external_counts(plan, work_dir)
        if dump_path is None:
            return (
                "unchecked",
                {
                    "kind": "binary_unparsed",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": None,
                    "checked": False,
                    "reason": unparsed_reason or "binary output has no parser",
                },
                None,
            )
        try:
            obs_counts = parity.parse_kmer_text_dump(dump_path, k=plan.k, min_abundance=plan.min_abundance)
        except ValueError as exc:
            return (
                "parity_failed",
                {
                    "kind": "set_and_abundance",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": False,
                    "checked": True,
                    "reason": f"invalid k-mer dump: {exc}",
                },
                None,
            )
        if plan.workload == "T1":
            ref = parity.count_file((cohort_abs / genomes[0]["fasta"]).resolve(), k=plan.k, min_abundance=plan.min_abundance)
        else:
            ref = parity.count_files(
                [(cohort_abs / genomes[0]["reads_r1"]).resolve(), (cohort_abs / genomes[0]["reads_r2"]).resolve()],
                k=plan.k,
                min_abundance=plan.min_abundance,
            )
        set_rep = parity.verify_set_parity(ref.counts.keys(), obs_counts.keys(), k=plan.k)
        ab_rep = parity.verify_abundance_parity(
            ref.counts, obs_counts, k=plan.k, min_abundance=plan.min_abundance, max_cap=plan.abundance_cap
        )
        ok = set_rep.ok and ab_rep.ok
        return (
            "ok" if ok else "parity_failed",
            {
                "kind": "set_and_abundance",
                "check_type": "real_tool_output",
                "self_check": False,
                "ok": ok,
                "checked": True,
                "set": asdict(set_rep),
                "abundance": asdict(ab_rep),
            },
            len(obs_counts) if obs_counts else None,
        )

    # Cardinality estimator (ntcard).
    if plan.tool == "ntcard":
        hist_path = pathlib.Path(plan.artifact_paths[0])
        if not hist_path.is_file() or hist_path.stat().st_size == 0:
            return (
                "failed_empty_output",
                {
                    "kind": "cardinality_estimator",
                    "check_type": "real_tool_output",
                    "self_check": False,
                    "ok": False,
                    "checked": True,
                    "reason": f"missing or empty histogram {hist_path.name}",
                },
                None,
            )
        return (
            "unchecked",
            {
                "kind": "cardinality_estimator",
                "check_type": "real_tool_output",
                "self_check": False,
                "ok": None,
                "checked": False,
                "reason": "ntCard estimates cardinality/histogram only; exact k-mer parity is not defined",
            },
            None,
        )

    # Remaining non-empty binary/colored formats (kmtricks T3/T5, ggcat T3/T5).
    return (
        "unchecked",
        {
            "kind": "binary_unparsed",
            "check_type": "real_tool_output",
            "self_check": False,
            "ok": None,
            "checked": False,
            "reason": f"no output parity parser implemented for {plan.tool} {plan.workload}",
        },
        None,
    )


def _grm_ksets_ready(grm_data: pathlib.Path, species_slug: str, genomes: Sequence[Mapping[str, Any]], k: int) -> bool:
    counted_path = grm_data / "species" / species_slug / "counted.json"
    if not counted_path.is_file():
        return False
    try:
        counted = json.loads(counted_path.read_text(encoding="utf-8"))
        if int(counted.get("k", -1)) != k:
            return False
    except (OSError, ValueError, json.JSONDecodeError):
        return False
    kset_dir = grm_data / "species" / species_slug / "ksets"
    for g in genomes:
        kp = kset_dir / f"{g['genome']}.kset"
        if not kp.is_file() or kp.stat().st_size == 0:
            return False
    return True


def _prepare_grm_prerequisites(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
) -> Optional[int]:
    """Stages grm_data and runs any un-timed prerequisites (count/matrix) for T3/T4 reuse (C2).

    Returns None on success or non-zero exit code if a prerequisite grm command fails.
    """
    if plan.grm_data_dir is None:
        return None
    grm_data = pathlib.Path(plan.grm_data_dir).resolve()
    stage_grm_data_dir(cohort_dir, cohort_manifest, grm_data, plan.species_slug)
    sp_dir = grm_data / "species" / plan.species_slug
    genomes = cohort_manifest["genomes"]
    work_dir = pathlib.Path(plan.work_dir).resolve()

    if plan.workload == "T1":
        # Clean previous count outputs so T1 benchmarks a genuine fresh count.
        kset_dir = sp_dir / "ksets"
        if kset_dir.exists():
            shutil.rmtree(kset_dir)
        counted = sp_dir / "counted.json"
        if counted.exists():
            counted.unlink()
        return None

    if plan.workload == "T3":
        mat_dir = sp_dir / "matrix"
        if mat_dir.exists():
            shutil.rmtree(mat_dir)
        if plan.variant == "disk-sparse-256":
            kset_dir = sp_dir / "ksets"
            if kset_dir.exists():
                shutil.rmtree(kset_dir)
            counted = sp_dir / "counted.json"
            if counted.exists():
                counted.unlink()
            return None
        if not _grm_ksets_ready(grm_data, plan.species_slug, genomes, plan.k):
            proc = subprocess.run(
                [
                    plan.binary, "--data", str(grm_data), "count",
                    "--species", plan.species_slug, "--k", str(plan.k), "--allow-missing-catalog",
                ],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                return proc.returncode
        return None

    if plan.workload == "T4":
        out_unitigs = pathlib.Path(plan.artifact_paths[0])
        if out_unitigs.exists():
            shutil.rmtree(out_unitigs)
        if not _grm_ksets_ready(grm_data, plan.species_slug, genomes, plan.k):
            proc = subprocess.run(
                [
                    plan.binary, "--data", str(grm_data), "count",
                    "--species", plan.species_slug, "--k", str(plan.k), "--allow-missing-catalog",
                ],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                return proc.returncode
        mat_json = sp_dir / "matrix" / "matrix.json"
        if not mat_json.is_file():
            proc = subprocess.run(
                [
                    plan.binary, "--data", str(grm_data), "matrix",
                    "--species", plan.species_slug, "--min-support", str(plan.min_support),
                ],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            if proc.returncode != 0:
                return proc.returncode
        return None

    return None


def execute_case(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
    *,
    mode: str = "real",
    skip_missing: bool = True,
    inject_fault: Optional[str] = None,
) -> Dict[str, Any]:
    """Runs one CommandPlan in 'dry-run', 'mock' (self-check), or 'real' mode and returns a result record."""
    if mode == "mock_self_check":
        mode = "mock"
    if mode not in ("dry-run", "mock", "real"):
        raise ValueError(f"unsupported mode: {mode!r}")

    cohort_abs = pathlib.Path(cohort_dir).resolve()
    prov = capture_tool_provenance(plan.tool, binary_override=plan.binary)
    genomes = cohort_manifest["genomes"]
    n_genomes = len(genomes)

    ref_kmers = _reference_kmer_count(plan, cohort_manifest, cohort_abs)
    run_parent = pathlib.Path(plan.work_dir).resolve().parent
    record_mode = "dry_run" if mode == "dry-run" else ("mock_self_check" if mode == "mock" else "real")

    base_record: Dict[str, Any] = {
        "workload": plan.workload,
        "tool": plan.tool,
        "variant": plan.variant,
        "k": plan.k,
        "min_abundance": plan.min_abundance,
        "threads": plan.threads,
        "mem_gb": plan.mem_gb,
        "cohort_genomes": n_genomes,
        "mode": record_mode,
        "comparable": plan.comparable,
        "synthetic_placeholder": mode == "mock",
        "provenance": asdict(prov),
        "commands": _display_commands(plan.steps, cohort_abs, run_parent),
        "status": "planned",
        "timing": None,
        "disk": None,
        "parity": None,
        "notes": plan.notes,
    }

    if mode == "dry-run":
        base_record["status"] = "dry_run" if plan.comparable else "not_comparable"
        return base_record

    if not plan.comparable:
        base_record["status"] = "not_comparable"
        base_record["notes"] = plan.not_comparable_reason or plan.notes
        base_record["parity"] = {
            "kind": "not_comparable",
            "check_type": record_mode,
            "self_check": mode == "mock",
            "ok": None,
            "checked": False,
            "reason": plan.not_comparable_reason,
        }
        return base_record

    if mode == "mock":
        timing = _mock_timing(plan.tool, plan.workload, ref_kmers, plan.threads)
        # Self-check only: both sides come from parity.py, not an external tool binary (C3).
        parity_result = _run_mock_self_check(plan, cohort_manifest, cohort_abs, inject_fault=inject_fault)
        artifact_bytes = max(1, ref_kmers * (6 if plan.tool == "grm" else 5))
        temp_bytes = max(1, ref_kmers * 16)
        disk = compute_disk_metrics(temp_bytes, artifact_bytes, ref_kmers)
        base_record["status"] = "ok" if parity_result["ok"] else "parity_failed"
        base_record["timing"] = asdict(timing)
        base_record["disk"] = asdict(disk)
        base_record["parity"] = parity_result
        base_record["notes"] = (
            f"[mock_self_check: synthetic placeholder timings/sizes; self-check only] {plan.notes}"
        )
        return base_record

    # Real execution mode.
    if not prov.installed:
        if skip_missing:
            base_record["status"] = "skipped_missing_binary"
            return base_record
        raise FileNotFoundError(f"required binary for tool {plan.tool!r} not found on PATH")

    time_bin = shutil.which("/usr/bin/time") or shutil.which("time")
    if time_bin is None or not pathlib.Path("/usr/bin/time").exists():
        raise RuntimeError("/usr/bin/time is required for real benchmark execution")

    work_dir = pathlib.Path(plan.work_dir).resolve()
    for d in [work_dir, work_dir / "out", *(pathlib.Path(t).resolve() for t in plan.temp_dirs)]:
        d.mkdir(parents=True, exist_ok=True)

    materialize_run_lists(plan, cohort_manifest, cohort_abs)
    if plan.tool == "grm":
        prereq_rc = _prepare_grm_prerequisites(plan, cohort_manifest, cohort_abs)
        if prereq_rc is not None:
            base_record["timing"] = asdict(
                TimeMetrics(wall_seconds=0.0, user_seconds=0.0, sys_seconds=0.0, max_rss_kb=0, exit_status=prereq_rc)
            )
            base_record["status"] = f"failed_exit_{prereq_rc}"
            return base_record

    sampler = DiskSampler(
        temp_paths=[pathlib.Path(t) for t in plan.temp_dirs],
        artifact_paths=[pathlib.Path(a) for a in plan.artifact_paths],
    )
    sampler.start()
    total_wall = 0.0
    total_user = 0.0
    total_sys = 0.0
    peak_rss = 0
    last_exit = 0
    try:
        for step in plan.steps:
            time_file = work_dir / "step.time"
            if time_file.exists():
                time_file.unlink()
            proc = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(time_file), *step],
                cwd=str(work_dir),
                capture_output=True,
                text=True,
                check=False,
            )
            last_exit = proc.returncode
            if time_file.exists():
                tm = parse_gnu_time_v(time_file.read_text(encoding="utf-8", errors="replace"))
                total_wall += tm.wall_seconds
                total_user += tm.user_seconds
                total_sys += tm.sys_seconds
                peak_rss = max(peak_rss, tm.max_rss_kb)
            if proc.returncode != 0:
                break
    finally:
        sampler.stop()

    final_bytes = directory_bytes([pathlib.Path(a) for a in plan.artifact_paths])
    base_record["timing"] = asdict(
        TimeMetrics(
            wall_seconds=round(total_wall, 6),
            user_seconds=round(total_user, 6),
            sys_seconds=round(total_sys, 6),
            max_rss_kb=peak_rss,
            exit_status=last_exit,
        )
    )

    if last_exit != 0:
        disk = compute_disk_metrics(
            peak_temp_bytes=sampler.peak_temp_bytes,
            final_artifact_bytes=final_bytes,
            kmers_counted=ref_kmers,
            peak_total_bytes=sampler.peak_total_bytes,
        )
        base_record["disk"] = asdict(disk)
        base_record["status"] = f"failed_exit_{last_exit}"
        return base_record

    status, parity_report, observed_kmers = _verify_real_output(plan, cohort_manifest, cohort_abs, final_bytes)
    disk = compute_disk_metrics(
        peak_temp_bytes=sampler.peak_temp_bytes,
        final_artifact_bytes=final_bytes,
        kmers_counted=observed_kmers if observed_kmers is not None else ref_kmers,
        peak_total_bytes=sampler.peak_total_bytes,
    )
    base_record["disk"] = asdict(disk)
    base_record["parity"] = parity_report
    base_record["status"] = status
    return base_record


def _run_mock_self_check(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
    inject_fault: Optional[str] = None,
) -> Dict[str, Any]:
    """Runs a synthetic parity.py self-check for mock mode (NOT external tool verification)."""
    genomes = cohort_manifest["genomes"]
    meta_fields = {
        "check_type": "mock_self_check",
        "self_check": True,
        "note": "synthetic self-check only (both sides computed by parity.py; not tool output)",
    }
    if not plan.supports_set_parity:
        return {"kind": "self_check:cardinality_estimator", "ok": True, "checked": False, **meta_fields}

    if plan.workload == "T1":
        fna = cohort_dir / genomes[0]["fasta"]
        res = parity.count_file(fna, k=plan.k, min_abundance=plan.min_abundance)
        obs = dict(res.counts)
        if inject_fault == "drop_kmer" and obs:
            obs.pop(min(obs))
        elif inject_fault == "alter_count" and obs:
            first = min(obs)
            obs[first] += 1
        set_rep = parity.verify_set_parity(res.counts.keys(), obs.keys(), k=plan.k)
        out: Dict[str, Any] = {"kind": "self_check:set", "ok": set_rep.ok, "set": asdict(set_rep), **meta_fields}
        if plan.supports_abundance_parity:
            ab_rep = parity.verify_abundance_parity(
                res.counts, obs, k=plan.k, min_abundance=plan.min_abundance, max_cap=plan.abundance_cap
            )
            out["kind"] = "self_check:set_and_abundance"
            out["abundance"] = asdict(ab_rep)
            out["ok"] = set_rep.ok and ab_rep.ok
        return out

    if plan.workload == "T2":
        r1 = cohort_dir / genomes[0]["reads_r1"]
        r2 = cohort_dir / genomes[0]["reads_r2"]
        res = parity.count_files([r1, r2], k=plan.k, min_abundance=plan.min_abundance)
        obs = dict(res.counts)
        if inject_fault == "drop_kmer" and obs:
            obs.pop(min(obs))
        elif inject_fault == "alter_count" and obs:
            first = min(obs)
            obs[first] += 1
        set_rep = parity.verify_set_parity(res.counts.keys(), obs.keys(), k=plan.k)
        out = {"kind": "self_check:set", "ok": set_rep.ok, "set": asdict(set_rep), **meta_fields}
        if plan.supports_abundance_parity:
            ab_rep = parity.verify_abundance_parity(
                res.counts, obs, k=plan.k, min_abundance=plan.min_abundance, max_cap=plan.abundance_cap
            )
            out["kind"] = "self_check:set_and_abundance"
            out["abundance"] = asdict(ab_rep)
            out["ok"] = set_rep.ok and ab_rep.ok
        return out

    if plan.workload in ("T3", "T5"):
        if plan.workload == "T3":
            sets = [set(parity.count_file(cohort_dir / g["fasta"], k=plan.k, min_abundance=1).counts) for g in genomes]
        else:
            sets = [
                set(
                    parity.count_files(
                        [cohort_dir / g["reads_r1"], cohort_dir / g["reads_r2"]],
                        k=plan.k,
                        min_abundance=plan.min_abundance,
                    ).counts
                )
                for g in genomes
            ]
        ref_mat = parity.build_reference_matrix(sets, k=plan.k, min_support=plan.min_support)
        obs_map = dict(ref_mat.kmer_to_mask)
        if inject_fault == "flip_bit" and obs_map:
            first = min(obs_map)
            obs_map[first] ^= 1
        m_rep = parity.verify_matrix_parity(ref_mat, obs_map)
        return {"kind": "self_check:matrix", "ok": m_rep.ok, "matrix": asdict(m_rep), **meta_fields}

    # T4 unitigs.
    union_kmers: set = set()
    for g in genomes:
        union_kmers.update(parity.count_file(cohort_dir / g["fasta"], k=plan.k, min_abundance=1).counts)
    ref_u = parity.compact_unitigs(union_kmers, k=plan.k)
    obs_u = list(ref_u)
    if inject_fault == "merge_unitig" and len(obs_u) >= 2:
        obs_u = [obs_u[0] | obs_u[1], *obs_u[2:]]
    u_rep = parity.verify_unitig_parity(ref_u, obs_u)
    return {"kind": "self_check:unitigs", "ok": u_rep.ok, "unitigs": asdict(u_rep), **meta_fields}


def build_bench_report(
    results: Sequence[Mapping[str, Any]],
    cohort_manifest: Mapping[str, Any],
    *,
    strict_hygiene: bool = True,
    forbidden_ids: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Assembles a deterministic, zero-genome-ID grm-kmer-bench-v1 report."""
    ordered = sorted(results, key=lambda r: (r["workload"], r["tool"], r["variant"], r.get("mode", "")))
    spec = cohort_manifest["spec"]
    report: Dict[str, Any] = {
        "format": BENCH_FORMAT,
        "cohort_summary": {
            "cohort_format": cohort_manifest["format"],
            "seed": spec["seed"],
            "n_genomes": spec["n_genomes"],
            "k": spec["k"],
            "min_support": spec["min_support"],
            "core_len": spec["core_len"],
            "n_plasmids": spec["n_plasmids"],
            "read_coverage": spec["read_coverage"] if spec["generate_reads"] else None,
            "error_rate": spec["error_rate"] if spec["generate_reads"] else None,
        },
        "results": ordered,
    }
    if strict_hygiene:
        leaks = scan_json_for_leaks(report, forbidden_ids=forbidden_ids)
        if leaks:
            raise ValueError(f"benchmark report failed hygiene scan: {', '.join(leaks)}")
    return report


TSV_COLUMNS = (
    "workload",
    "tool",
    "variant",
    "mode",
    "status",
    "comparable",
    "k",
    "min_abundance",
    "threads",
    "cohort_genomes",
    "wall_seconds",
    "user_seconds",
    "sys_seconds",
    "max_rss_kb",
    "peak_temp_bytes",
    "final_artifact_bytes",
    "bytes_per_kmer",
    "parity_ok",
    "binary_sha256",
)


def filter_real_comparable_results(results: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    """Returns only real-mode, comparable benchmark rows (excluding mock_self_check and not_comparable)."""
    return [
        r
        for r in results
        if r.get("mode") == "real"
        and r.get("comparable", True) is True
        and r.get("status") != "not_comparable"
    ]


def render_summary_tsv(
    results: Sequence[Mapping[str, Any]],
    *,
    include_mock: bool = True,
    include_not_comparable: bool = True,
) -> str:
    """Renders deterministic TSV rows ordered by (workload, tool, variant, mode)."""
    rows_in = list(results)
    if not include_mock:
        rows_in = [r for r in rows_in if r.get("mode") == "real"]
    if not include_not_comparable:
        rows_in = [r for r in rows_in if r.get("comparable", True) is True and r.get("status") != "not_comparable"]
    ordered = sorted(rows_in, key=lambda r: (r["workload"], r["tool"], r["variant"], r.get("mode", "")))
    lines = ["\t".join(TSV_COLUMNS)]
    for r in ordered:
        t = r.get("timing") or {}
        d = r.get("disk") or {}
        p = r.get("parity") or {}
        prov = r.get("provenance") or {}

        def _fmt(v: Any) -> str:
            if v is None:
                return ""
            if isinstance(v, bool):
                return "true" if v else "false"
            if isinstance(v, float):
                return f"{v:.6f}"
            return str(v)

        row = [
            _fmt(r.get("workload")),
            _fmt(r.get("tool")),
            _fmt(r.get("variant")),
            _fmt(r.get("mode")),
            _fmt(r.get("status")),
            _fmt(r.get("comparable", True)),
            _fmt(r.get("k")),
            _fmt(r.get("min_abundance")),
            _fmt(r.get("threads")),
            _fmt(r.get("cohort_genomes")),
            _fmt(t.get("wall_seconds")),
            _fmt(t.get("user_seconds")),
            _fmt(t.get("sys_seconds")),
            _fmt(t.get("max_rss_kb")),
            _fmt(d.get("peak_temp_bytes")),
            _fmt(d.get("final_artifact_bytes")),
            _fmt(d.get("bytes_per_kmer")),
            _fmt(p.get("ok")),
            _fmt(prov.get("binary_sha256")),
        ]
        lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def render_summary_markdown(
    results: Sequence[Mapping[str, Any]],
    *,
    include_mock: bool = False,
    include_not_comparable: bool = False,
) -> str:
    """Renders a neutral Markdown benchmark summary table.

    By default excludes mock_self_check, dry_run, and not_comparable rows so synthetic
    placeholder metrics never appear in real comparison tables (C3, C12).
    """
    filtered = list(results)
    excluded_mock = sum(1 for r in filtered if r.get("mode") != "real")
    excluded_nc = sum(
        1 for r in filtered if r.get("mode") == "real" and (r.get("comparable", True) is False or r.get("status") == "not_comparable")
    )
    if not include_mock:
        filtered = [r for r in filtered if r.get("mode") == "real"]
    if not include_not_comparable:
        filtered = [r for r in filtered if r.get("comparable", True) is True and r.get("status") != "not_comparable"]

    ordered = sorted(filtered, key=lambda r: (r["workload"], r["tool"], r["variant"], r.get("mode", "")))
    header = "| Workload | Tool | Variant | Mode | Status | Wall (s) | Max RSS (KB) | Artifact (B) | B/k-mer | Parity |"
    sep = "| --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |"
    rows = [header, sep]
    for r in ordered:
        t = r.get("timing") or {}
        d = r.get("disk") or {}
        p = r.get("parity") or {}
        wall = f"{t['wall_seconds']:.3f}" if "wall_seconds" in t and t["wall_seconds"] is not None else "-"
        rss = str(t["max_rss_kb"]) if "max_rss_kb" in t and t["max_rss_kb"] is not None else "-"
        art = str(d["final_artifact_bytes"]) if "final_artifact_bytes" in d and d["final_artifact_bytes"] is not None else "-"
        bpk = f"{d['bytes_per_kmer']:.2f}" if d.get("bytes_per_kmer") is not None else "-"
        if r.get("mode") == "mock_self_check":
            par = "self_check_pass" if p.get("ok") is True else ("self_check_fail" if p.get("ok") is False else "-")
        else:
            par = "pass" if p.get("ok") is True else ("fail" if p.get("ok") is False else "unchecked")
        rows.append(
            f"| {r['workload']} | `{r['tool']}` | `{r['variant']}` | `{r.get('mode', '')}` | {r['status']} | {wall} | {rss} | {art} | {bpk} | {par} |"
        )
    notes: List[str] = []
    if not include_mock and excluded_mock > 0:
        notes.append(f"Excluded {excluded_mock} non-real (mock_self_check/dry_run) row(s) from comparison.")
    if not include_not_comparable and excluded_nc > 0:
        notes.append(f"Excluded {excluded_nc} not_comparable row(s) from comparison.")
    body = "\n".join(rows) + "\n"
    if notes:
        body += "\n" + " ".join(notes) + "\n"
    return body


def write_bench_outputs(
    out_dir: pathlib.Path,
    report: Mapping[str, Any],
    *,
    strict_hygiene: bool = True,
    forbidden_ids: Optional[Iterable[str]] = None,
) -> Tuple[pathlib.Path, pathlib.Path]:
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_text = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n"
    tsv_text = render_summary_tsv(report["results"], include_mock=True, include_not_comparable=True)
    if strict_hygiene:
        json_leaks = scan_json_for_leaks(report, forbidden_ids=forbidden_ids)
        if json_leaks:
            raise ValueError(f"bench.json failed hygiene scan: {', '.join(json_leaks)}")
        tsv_leaks = scan_tsv_for_leaks(tsv_text, forbidden_ids=forbidden_ids)
        if tsv_leaks:
            raise ValueError(f"summary.tsv failed hygiene scan: {', '.join(tsv_leaks)}")
    json_path = out_dir / "bench.json"
    tsv_path = out_dir / "summary.tsv"
    json_path.write_text(json_text, encoding="utf-8")
    tsv_path.write_text(tsv_text, encoding="utf-8")
    return json_path, tsv_path
