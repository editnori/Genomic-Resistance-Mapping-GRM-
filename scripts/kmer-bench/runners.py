"""Benchmark runners, command builders, resource parsers, and grm-kmer-bench-v1 writer.

Supports workloads T1..T5 across grm, fastk, kmc3, jellyfish, ntcard, kmtricks,
ggcat, cuttlefish2, and bcalm2 using Python's standard library only.
"""

from dataclasses import asdict, dataclass
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
import time
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
        "version_args": ["--version"],
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

# Genome-id and private-detail hygiene regexes (aligned with scripts/c7-outside/scan.py
# and tests/public_words/rules.py).
GENOME_ID_TOKEN = re.compile(r"(?<![\w.])\d+\.\d+(?!\d)")
BVBRC_ID_SHAPE = re.compile(r"(?<![\w.])\d{3,7}\.\d{3,6}(?!\d)")
_S = "/"
HOME_PATH_RE = re.compile(_S + "home" + _S + r"([A-Za-z0-9_.-]+)|" + _S + "Users" + _S + r"([A-Za-z0-9_.-]+)")
NEUTRAL_HOMES = frozenset({"user", "runner", "proj"})
WINDOWS_USER_RE = re.compile(r"[A-Za-z]:\\+Users\\+|" + _S + "mnt" + _S + r"[a-z]" + _S + "Users" + _S)
THREAD_LINK_RE = re.compile("t3" + "-" + r"[a-z]+:" + _S + _S)
BOX_ID_RE = re.compile(r"\bbx" + "_" + r"[a-z0-9]{8}\b")


def scan_for_leaks(
    text: str,
    forbidden_ids: Optional[Iterable[str]] = None,
    check_bvbrc_shape: bool = True,
) -> List[str]:
    """Returns a list of leak violation tags found in `text` (never prints offending tokens)."""
    violations: List[str] = []
    if "\u2014" in text:
        violations.append("em_dash")
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
    if forbidden_ids is not None:
        id_set = set(forbidden_ids)
        if any(m.group(0) in id_set for m in GENOME_ID_TOKEN.finditer(text)):
            violations.append("known_genome_id")
    if check_bvbrc_shape and BVBRC_ID_SHAPE.search(text):
        violations.append("bvbrc_id_shape")
    return violations


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


def _sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def capture_tool_provenance(
    tool: str,
    binary_override: Optional[str] = None,
) -> ToolProvenance:
    """Captures binary availability, SHA-256, and version string without leaking user paths."""
    if tool not in TOOL_REGISTRY:
        raise ValueError(f"unknown tool: {tool!r}")
    spec = TOOL_REGISTRY[tool]
    bin_name = binary_override or spec["binary"]
    resolved = shutil.which(bin_name)
    if resolved is None:
        return ToolProvenance(
            tool=tool,
            binary_name=pathlib.Path(bin_name).name,
            installed=False,
            binary_sha256=None,
            version_line=None,
        )
    bin_path = pathlib.Path(resolved).resolve()
    sha = _sha256_file(bin_path) if bin_path.is_file() else None
    version_line: Optional[str] = None
    vargs = list(spec["version_args"])
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
                # Strip any home paths if a binary prints argv[0].
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
    )


@dataclass(frozen=True)
class CommandPlan:
    tool: str
    workload: str
    variant: str
    k: int
    min_abundance: int
    threads: int
    mem_gb: int
    steps: Tuple[Tuple[str, ...], ...]
    work_dir: str
    artifact_paths: Tuple[str, ...]
    temp_dirs: Tuple[str, ...]
    supports_set_parity: bool
    supports_abundance_parity: bool
    abundance_cap: Optional[int]
    supports_unitig_parity: bool
    supports_matrix_parity: bool
    notes: str


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
) -> CommandPlan:
    """Builds exact argv command steps for one (tool, workload, variant) benchmark case."""
    if tool not in TOOL_REGISTRY:
        raise ValueError(f"unknown tool: {tool!r}")
    if workload not in WORKLOADS:
        raise ValueError(f"unknown workload: {workload!r}")
    meta = TOOL_REGISTRY[tool]
    if workload not in meta["workloads"]:
        raise ValueError(f"tool {tool!r} does not support workload {workload!r}")

    eff_min = (1 if workload in ("T1", "T3", "T4") else 2) if min_abundance is None else min_abundance
    bin_cmd = binary_override or meta["binary"]
    work = run_dir / f"{workload}_{tool}_{variant}"
    tmp = work / "tmp"
    out = work / "out"
    asm_fof = cohort_dir / "assemblies.fof"
    reads_fof = cohort_dir / "reads.kmtricks.fof"
    sample_fasta = cohort_dir / "assemblies" / "syn_g0000.fna"
    sample_r1 = cohort_dir / "reads" / "syn_g0000_R1.fastq"
    sample_r2 = cohort_dir / "reads" / "syn_g0000_R2.fastq"

    steps: List[Tuple[str, ...]] = []
    artifacts: List[str] = []
    notes = ""

    if tool == "grm":
        if workload == "T1":
            steps.append((bin_cmd, "--data", str(work / "data"), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            artifacts.append(str(work / "data" / "species" / species_slug / "ksets"))
            notes = "Counts all cohort assemblies into .kset v1 files with rayon parallelism."
        elif workload == "T2":
            steps.append((bin_cmd, "--data", str(work / "data"), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            artifacts.append(str(work / "data" / "species" / species_slug / "ksets"))
            notes = "Baseline GRM count on FASTQ (min_abundance=1 in stock CLI; Tier-2 adds streaming min_abundance>=2 + u16 counts)."
        elif workload == "T3":
            cmd = [bin_cmd, "--data", str(work / "data"), "matrix", "--species", species_slug, "--min-support", str(min_support)]
            if variant == "column-blocks-64":
                cmd.extend(["--column-blocks", "64"])
            elif variant == "disk-sparse-256":
                cmd.extend(["--disk-sparse", "--block-width", "256", "--allow-missing-catalog"])
            steps.append(tuple(cmd))
            artifacts.append(str(work / "data" / "species" / species_slug / "matrix"))
            notes = f"GRM multi-genome pattern matrix build ({variant})."
        elif workload == "T4":
            steps.append((
                bin_cmd, "--data", str(work / "data"), "exp-input-build", "unitigs",
                "--species", species_slug, "--graph", "grm", "--out", str(out / "unitigs"),
            ))
            artifacts.append(str(out / "unitigs"))
            notes = "GRM native union + prefix-indexed de Bruijn graph compaction."
        elif workload == "T5":
            steps.append((bin_cmd, "--data", str(work / "data"), "count", "--species", species_slug, "--k", str(k), "--allow-missing-catalog"))
            steps.append((bin_cmd, "--data", str(work / "data"), "matrix", "--species", species_slug, "--min-support", str(min_support)))
            artifacts.append(str(work / "data" / "species" / species_slug / "matrix"))
            notes = "Two-stage GRM count + matrix baseline on FASTQ inputs."

    elif tool == "fastk":
        prefix = out / "fastk_out"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, f"-k{k}", f"-t{eff_min}", f"-T{threads}", f"-M{mem_gb}", f"-P{tmp}", f"-N{prefix}", *inp))
        artifacts.append(f"{prefix}.ktab")
        notes = "FastK super-k-mer table; default .ktab saturates counts at 255."

    elif tool == "kmc3":
        prefix = out / "kmc_db"
        fmt_flag = "-fm" if workload == "T1" else "-fq"
        inp = str(sample_fasta) if workload == "T1" else f"@{work / 'sample_reads.fof'}"
        steps.append((bin_cmd, f"-k{k}", f"-ci{eff_min}", "-cs65535", f"-t{threads}", f"-m{mem_gb}", fmt_flag, inp, str(prefix), str(tmp)))
        artifacts.extend([f"{prefix}.kmc_pre", f"{prefix}.kmc_suf"])
        notes = "KMC3 two-stage disk-binned k-mer counter with u16 saturation cap."

    elif tool == "jellyfish":
        jf_out = out / "mer_counts.jf"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, "count", "-m", str(k), "-s", "100M", "-t", str(threads), "-C", "-o", str(jf_out), *inp))
        artifacts.append(str(jf_out))
        notes = "Jellyfish 2 lock-free in-memory hash counter (-C canonical)."

    elif tool == "ntcard":
        prefix = out / "ntcard_hist"
        inp = [str(sample_fasta)] if workload == "T1" else [str(sample_r1), str(sample_r2)]
        steps.append((bin_cmd, f"-k{k}", f"-t{threads}", "-p", str(prefix), *inp))
        artifacts.append(f"{prefix}_k{k}.hist")
        notes = "ntCard streaming cardinality/histogram estimator (F0 and f1..f64 only, no exact set dump)."

    elif tool == "kmtricks":
        fof = asm_fof if workload == "T3" else reads_fof
        steps.append((
            bin_cmd, "pipeline", "--file", str(fof), "--run-dir", str(out / "kmtricks_run"),
            "--kmer-size", str(k), "--hard-min", str(eff_min), "--threads", str(threads),
            "--mode", "kmer:bin:hard",
        ))
        artifacts.append(str(out / "kmtricks_run"))
        notes = "kmtricks multi-sample k-mer presence matrix pipeline."

    elif tool == "ggcat":
        fof = asm_fof if workload in ("T3", "T4") else reads_fof
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
        notes = "bcalm2 2.2.3 minimizer-partitioned compacted de Bruijn graph builder."

    return CommandPlan(
        tool=tool,
        workload=workload,
        variant=variant,
        k=k,
        min_abundance=eff_min,
        threads=threads,
        mem_gb=mem_gb,
        steps=tuple(steps),
        work_dir=str(work),
        artifact_paths=tuple(artifacts),
        temp_dirs=(str(tmp),),
        supports_set_parity=bool(meta["exact_set_parity"]),
        supports_abundance_parity=bool(meta["exact_abundance_parity"]) and workload in ("T1", "T2"),
        abundance_cap=meta["abundance_cap"],
        supports_unitig_parity=bool(meta["supports_unitigs"]) and workload == "T4",
        supports_matrix_parity=bool(meta["supports_matrix"]) and workload in ("T3", "T5"),
        notes=notes,
    )


def _mock_timing(tool: str, workload: str, n_kmers: int, threads: int) -> TimeMetrics:
    """Deterministic synthetic /usr/bin/time -v telemetry for CI --mock runs."""
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
        for token in step:
            val = token
            for src, dst in replacements:
                if src and src in val:
                    val = val.replace(src, dst)
            norm_step.append(val.replace("\\", "/"))
        out.append(norm_step)
    return out


def execute_case(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
    *,
    mode: str = "real",
    skip_missing: bool = True,
    inject_fault: Optional[str] = None,
) -> Dict[str, Any]:
    """Runs one CommandPlan in 'dry-run', 'mock', or 'real' mode and returns a result record."""
    if mode not in ("dry-run", "mock", "real"):
        raise ValueError(f"unsupported mode: {mode!r}")

    prov = capture_tool_provenance(plan.tool)
    genomes = cohort_manifest["genomes"]
    n_genomes = len(genomes)

    # Determine reference k-mer count for bytes/k-mer accounting.
    if plan.workload == "T1":
        ref_kmers = int(genomes[0]["assembly_truth"]["distinct_kmers_min1"])
    elif plan.workload == "T2":
        ref_kmers = int(genomes[0]["reads_truth"]["distinct_kmers_min2"]) if "reads_truth" in genomes[0] else 0
    elif plan.workload == "T3":
        ref_kmers = int(cohort_manifest["cohort_truth"]["assembly_matrix"]["kmers_kept"])
    elif plan.workload == "T4":
        ref_kmers = int(cohort_manifest["cohort_truth"].get("unitigs", {}).get("graph_kmers", 0))
    else:
        ref_kmers = int(cohort_manifest["cohort_truth"].get("reads_matrix_min2", {}).get("kmers_kept", 0))

    run_parent = pathlib.Path(plan.work_dir).parent
    base_record: Dict[str, Any] = {
        "workload": plan.workload,
        "tool": plan.tool,
        "variant": plan.variant,
        "k": plan.k,
        "min_abundance": plan.min_abundance,
        "threads": plan.threads,
        "mem_gb": plan.mem_gb,
        "cohort_genomes": n_genomes,
        "mode": mode,
        "provenance": asdict(prov),
        "commands": _display_commands(plan.steps, cohort_dir, run_parent),
        "status": "planned",
        "timing": None,
        "disk": None,
        "parity": None,
        "notes": plan.notes,
    }

    if mode == "dry-run":
        base_record["status"] = "dry_run"
        return base_record

    if mode == "mock":
        timing = _mock_timing(plan.tool, plan.workload, ref_kmers, plan.threads)
        # Perform real reference parity computation on the synthetic cohort in mock mode.
        parity_result = _run_mock_parity(plan, cohort_manifest, cohort_dir, inject_fault=inject_fault)
        artifact_bytes = max(1, ref_kmers * (6 if plan.tool == "grm" else 5))
        temp_bytes = max(1, ref_kmers * 16)
        disk = compute_disk_metrics(temp_bytes, artifact_bytes, ref_kmers)
        base_record["status"] = "ok" if parity_result["ok"] else "parity_failed"
        base_record["timing"] = asdict(timing)
        base_record["disk"] = asdict(disk)
        base_record["parity"] = parity_result
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

    work_dir = pathlib.Path(plan.work_dir)
    for d in [work_dir, work_dir / "out", *(pathlib.Path(t) for t in plan.temp_dirs)]:
        d.mkdir(parents=True, exist_ok=True)
    if plan.tool == "kmc3" and plan.workload == "T2" and "reads_r1" in genomes[0]:
        fof = work_dir / "sample_reads.fof"
        fof.write_text(
            f"{cohort_dir / genomes[0]['reads_r1']}\n{cohort_dir / genomes[0]['reads_r2']}\n",
            encoding="utf-8",
        )

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
            proc = subprocess.run(
                ["/usr/bin/time", "-v", "-o", str(time_file), *step],
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
    disk = compute_disk_metrics(
        peak_temp_bytes=sampler.peak_temp_bytes,
        final_artifact_bytes=final_bytes,
        kmers_counted=ref_kmers,
        peak_total_bytes=sampler.peak_total_bytes,
    )
    base_record["timing"] = asdict(
        TimeMetrics(
            wall_seconds=round(total_wall, 6),
            user_seconds=round(total_user, 6),
            sys_seconds=round(total_sys, 6),
            max_rss_kb=peak_rss,
            exit_status=last_exit,
        )
    )
    base_record["disk"] = asdict(disk)
    base_record["status"] = "ok" if last_exit == 0 else f"failed_exit_{last_exit}"
    return base_record


def _run_mock_parity(
    plan: CommandPlan,
    cohort_manifest: Mapping[str, Any],
    cohort_dir: pathlib.Path,
    inject_fault: Optional[str] = None,
) -> Dict[str, Any]:
    """Executes genuine reference parity verification on synthetic cohort files for mock mode."""
    genomes = cohort_manifest["genomes"]
    if not plan.supports_set_parity:
        return {"kind": "cardinality_estimator", "ok": True, "checked": False}

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
        out: Dict[str, Any] = {"kind": "set", "ok": set_rep.ok, "set": asdict(set_rep)}
        if plan.supports_abundance_parity:
            ab_rep = parity.verify_abundance_parity(
                res.counts, obs, k=plan.k, min_abundance=plan.min_abundance, max_cap=plan.abundance_cap
            )
            out["kind"] = "set_and_abundance"
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
        out = {"kind": "set", "ok": set_rep.ok, "set": asdict(set_rep)}
        if plan.supports_abundance_parity:
            ab_rep = parity.verify_abundance_parity(
                res.counts, obs, k=plan.k, min_abundance=plan.min_abundance, max_cap=plan.abundance_cap
            )
            out["kind"] = "set_and_abundance"
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
        ref_mat = parity.build_reference_matrix(sets, k=plan.k, min_support=2)
        obs_map = dict(ref_mat.kmer_to_mask)
        if inject_fault == "flip_bit" and obs_map:
            first = min(obs_map)
            obs_map[first] ^= 1
        m_rep = parity.verify_matrix_parity(ref_mat, obs_map)
        return {"kind": "matrix", "ok": m_rep.ok, "matrix": asdict(m_rep)}

    # T4 unitigs.
    union_kmers: set = set()
    for g in genomes:
        union_kmers.update(parity.count_file(cohort_dir / g["fasta"], k=plan.k, min_abundance=1).counts)
    ref_u = parity.compact_unitigs(union_kmers, k=plan.k)
    obs_u = list(ref_u)
    if inject_fault == "merge_unitig" and len(obs_u) >= 2:
        obs_u = [obs_u[0] | obs_u[1], *obs_u[2:]]
    u_rep = parity.verify_unitig_parity(ref_u, obs_u)
    return {"kind": "unitigs", "ok": u_rep.ok, "unitigs": asdict(u_rep)}


def build_bench_report(
    results: Sequence[Mapping[str, Any]],
    cohort_manifest: Mapping[str, Any],
    *,
    strict_hygiene: bool = True,
    forbidden_ids: Optional[Iterable[str]] = None,
) -> Dict[str, Any]:
    """Assembles a deterministic, zero-genome-ID grm-kmer-bench-v1 report."""
    ordered = sorted(results, key=lambda r: (r["workload"], r["tool"], r["variant"]))
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
    serialized = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if strict_hygiene:
        leaks = scan_for_leaks(serialized, forbidden_ids=forbidden_ids)
        if leaks:
            raise ValueError(f"benchmark report failed hygiene scan: {', '.join(leaks)}")
    return report


TSV_COLUMNS = (
    "workload",
    "tool",
    "variant",
    "status",
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


def render_summary_tsv(results: Sequence[Mapping[str, Any]]) -> str:
    """Renders deterministic TSV rows ordered by (workload, tool, variant)."""
    ordered = sorted(results, key=lambda r: (r["workload"], r["tool"], r["variant"]))
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
            _fmt(r.get("status")),
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


def render_summary_markdown(results: Sequence[Mapping[str, Any]]) -> str:
    """Renders a neutral Markdown benchmark summary table."""
    ordered = sorted(results, key=lambda r: (r["workload"], r["tool"], r["variant"]))
    header = "| Workload | Tool | Variant | Status | Wall (s) | Max RSS (KB) | Artifact (B) | B/k-mer | Parity |"
    sep = "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | --- |"
    rows = [header, sep]
    for r in ordered:
        t = r.get("timing") or {}
        d = r.get("disk") or {}
        p = r.get("parity") or {}
        wall = f"{t['wall_seconds']:.3f}" if "wall_seconds" in t and t["wall_seconds"] is not None else "-"
        rss = str(t["max_rss_kb"]) if "max_rss_kb" in t and t["max_rss_kb"] is not None else "-"
        art = str(d["final_artifact_bytes"]) if "final_artifact_bytes" in d and d["final_artifact_bytes"] is not None else "-"
        bpk = f"{d['bytes_per_kmer']:.2f}" if d.get("bytes_per_kmer") is not None else "-"
        par = "pass" if p.get("ok") is True else ("fail" if p.get("ok") is False else "-")
        rows.append(
            f"| {r['workload']} | `{r['tool']}` | `{r['variant']}` | {r['status']} | {wall} | {rss} | {art} | {bpk} | {par} |"
        )
    return "\n".join(rows) + "\n"


def write_bench_outputs(
    out_dir: pathlib.Path,
    report: Mapping[str, Any],
    *,
    strict_hygiene: bool = True,
    forbidden_ids: Optional[Iterable[str]] = None,
) -> Tuple[pathlib.Path, pathlib.Path]:
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    tsv_text = render_summary_tsv(report["results"])
    if strict_hygiene:
        for label, body in (("bench.json", json_text), ("summary.tsv", tsv_text)):
            leaks = scan_for_leaks(body, forbidden_ids=forbidden_ids)
            if leaks:
                raise ValueError(f"{label} failed hygiene scan: {', '.join(leaks)}")
    json_path = out_dir / "bench.json"
    tsv_path = out_dir / "summary.tsv"
    json_path.write_text(json_text, encoding="utf-8")
    tsv_path.write_text(tsv_text, encoding="utf-8")
    return json_path, tsv_path
