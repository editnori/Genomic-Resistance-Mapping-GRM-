"""Deterministic synthetic bacterial cohort generator for k-mer benchmarking.

Generates multi-contig bacterial FASTA assemblies (shared core genome with SNPs,
accessory genes, assembly gaps, soft-masked stretches, and collapsed single-copy
plasmids) plus optional paired-end noisy FASTQ reads with configurable coverage,
sequencing error rate, and true multi-copy plasmid ratios, along with exact
ground-truth canonical k-mer counts and digests.

Zero external dependencies; strictly independent of PYTHONHASHSEED.
"""

from dataclasses import asdict, dataclass
import gzip
import importlib.util
import json
import pathlib
import statistics
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple


def _load_parity():
    mod_name = "_grm_kmer_bench_parity"
    if mod_name in sys.modules:
        return sys.modules[mod_name]
    path = pathlib.Path(__file__).resolve().parent / "parity.py"
    spec = importlib.util.spec_from_file_location(mod_name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[mod_name] = mod
    spec.loader.exec_module(mod)
    return mod


parity = _load_parity()

COHORT_FORMAT = "grm-kmer-bench-cohort-v1"
_MASK64 = 0xFFFF_FFFF_FFFF_FFFF
_BASES = "ACGT"
_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


class DeterministicRng:
    """Explicit SplitMix64 PRNG independent of Python hash seed or stdlib random."""

    __slots__ = ("state",)

    def __init__(self, seed: int) -> None:
        self.state = seed & _MASK64

    def next_u64(self) -> int:
        self.state = (self.state + 0x9E37_79B9_7F4A_7C15) & _MASK64
        z = self.state
        z = ((z ^ (z >> 30)) * 0xBF58_476D_1CE4_E5B9) & _MASK64
        z = ((z ^ (z >> 27)) * 0x94D0_49BB_1331_11EB) & _MASK64
        return (z ^ (z >> 31)) & _MASK64

    def randrange(self, stop: int) -> int:
        if stop <= 0:
            raise ValueError(f"stop must be > 0, got {stop}")
        return self.next_u64() % stop

    def uniform(self) -> float:
        return (self.next_u64() >> 11) * (1.0 / (1 << 53))

    def random_dna(self, length: int) -> str:
        out: List[str] = []
        rem = length
        while rem > 0:
            word = self.next_u64()
            take = min(rem, 32)
            for _ in range(take):
                out.append(_BASES[word & 3])
                word >>= 2
            rem -= take
        return "".join(out)


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def mutate_dna(seq: str, snp_rate: float, rng: DeterministicRng) -> str:
    if snp_rate <= 0.0:
        return seq
    chars = list(seq)
    for i, ch in enumerate(chars):
        if ch in _BASES and rng.uniform() < snp_rate:
            alt = ( _BASES.index(ch) + 1 + rng.randrange(3) ) & 3
            chars[i] = _BASES[alt]
    return "".join(chars)


@dataclass(frozen=True)
class CohortSpec:
    seed: int = 42
    n_genomes: int = 8
    k: int = 31
    min_support: int = 2
    core_len: int = 4800
    n_core_contigs: int = 3
    n_accessory_genes: int = 6
    accessory_len: int = 600
    snp_rate: float = 0.005
    n_plasmids: int = 2
    plasmid_len: int = 1000
    plasmid_copy_numbers: Tuple[int, ...] = (4, 8)
    collapse_plasmids_in_assembly: bool = True
    include_n_runs: bool = True
    include_lowercase: bool = True
    generate_reads: bool = True
    read_coverage: float = 12.0
    read_length: int = 120
    insert_size: int = 220
    error_rate: float = 0.005
    gzip_reads: bool = False
    write_truth_tables: bool = True
    compute_unitig_truth: bool = True

    def validate(self) -> None:
        if not (parity.MIN_K <= self.k <= parity.MAX_K):
            raise ValueError(f"k must be in {parity.MIN_K}..={parity.MAX_K}")
        if self.n_genomes < 2:
            raise ValueError("n_genomes must be >= 2")
        if self.core_len < self.n_core_contigs * (self.k + 20):
            raise ValueError("core_len is too small for n_core_contigs and k")
        if self.n_plasmids > 0 and len(self.plasmid_copy_numbers) < self.n_plasmids:
            raise ValueError("plasmid_copy_numbers must provide at least n_plasmids entries")
        if self.generate_reads:
            if self.read_length < self.k:
                raise ValueError("read_length must be >= k")
            if self.insert_size < self.read_length:
                raise ValueError("insert_size must be >= read_length")
            if not (0.0 <= self.error_rate < 0.25):
                raise ValueError("error_rate must be in [0.0, 0.25)")
            if self.read_coverage <= 0.0:
                raise ValueError("read_coverage must be > 0")


def _wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def _write_fasta(path: pathlib.Path, contigs: Sequence[Tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: List[str] = []
    for name, seq in contigs:
        lines.append(f">{name}")
        lines.append(_wrap_fasta(seq))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _simulate_paired_reads(
    templates: Sequence[Tuple[str, str, int]],
    coverage: float,
    read_length: int,
    insert_size: int,
    error_rate: float,
    rng: DeterministicRng,
) -> Tuple[List[Tuple[str, str, str]], List[Tuple[str, str, str]], int]:
    """Simulates deterministic paired-end reads weighted by each template's true copy number."""
    r1: List[Tuple[str, str, str]] = []
    r2: List[Tuple[str, str, str]] = []
    errors = 0
    read_idx = 0

    def _apply_errors(seq: str) -> Tuple[str, str, int]:
        chars = list(seq)
        quals: List[str] = []
        n_err = 0
        for i, ch in enumerate(chars):
            if ch in _BASES and error_rate > 0.0 and rng.uniform() < error_rate:
                alt = (_BASES.index(ch) + 1 + rng.randrange(3)) & 3
                chars[i] = _BASES[alt]
                quals.append("#")  # Phred Q2 at injected error base
                n_err += 1
            else:
                quals.append("I")  # Phred Q40 at clean base
        return "".join(chars), "".join(quals), n_err

    for _tmpl_name, raw_seq, copy_number in templates:
        # Strip N runs for read generation from intact molecules.
        clean = "".join(ch.upper() for ch in raw_seq if ch.upper() in _BASES)
        if len(clean) < read_length:
            continue
        eff_insert = min(insert_size, len(clean))
        n_pairs = max(1, int(round((len(clean) * coverage * copy_number) / (2.0 * read_length))))
        max_start = len(clean) - eff_insert
        for _ in range(n_pairs):
            start = rng.randrange(max_start + 1) if max_start > 0 else 0
            frag = clean[start : start + eff_insert]
            if (rng.next_u64() & 1) == 1:
                frag = revcomp(frag)
            s1, q1, e1 = _apply_errors(frag[:read_length])
            s2, q2, e2 = _apply_errors(revcomp(frag[-read_length:]))
            errors += e1 + e2
            rid = f"syn_read_{read_idx:07d}"
            r1.append((f"{rid}/1", s1, q1))
            r2.append((f"{rid}/2", s2, q2))
            read_idx += 1

    return r1, r2, errors


def _write_fastq(path: pathlib.Path, reads: Sequence[Tuple[str, str, str]], gzipped: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    text = "".join(f"@{hdr}\n{seq}\n+\n{qual}\n" for hdr, seq, qual in reads)
    if gzipped:
        # mtime=0 makes gzip output 100% deterministic across runs.
        with path.open("wb") as raw_fh:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, compresslevel=6, mtime=0) as gz:
                gz.write(text.encode("ascii"))
    else:
        path.write_text(text, encoding="ascii")


def _write_truth_counts(path: pathlib.Path, counts: Dict[int, int], space: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [f"{space.decode(km)}\t{counts[km]}\n" for km in sorted(counts)]
    path.write_text("".join(rows), encoding="ascii")


def generate_cohort(out_dir: pathlib.Path, spec: Optional[CohortSpec] = None) -> Dict[str, Any]:
    """Generates a deterministic synthetic bacterial cohort and ground-truth manifest."""
    spec = spec or CohortSpec()
    spec.validate()
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = DeterministicRng(spec.seed)
    space = parity.KmerSpace(spec.k)

    # Ancestral core contigs, accessory genes, and plasmids.
    base_chunk = spec.core_len // spec.n_core_contigs
    core_ancestors: List[str] = []
    for c in range(spec.n_core_contigs):
        clen = base_chunk + (spec.core_len - base_chunk * spec.n_core_contigs if c == spec.n_core_contigs - 1 else 0)
        core_ancestors.append(rng.random_dna(clen))

    accessory_pool = [rng.random_dna(spec.accessory_len) for _ in range(spec.n_accessory_genes)]
    plasmid_pool = [rng.random_dna(spec.plasmid_len) for _ in range(spec.n_plasmids)]

    # Precompute unique k-mer sets for chromosome vs each plasmid to measure copy-number depth.
    plasmid_kmers = [set(space.iter_canonical(p)) for p in plasmid_pool]

    assemblies_dir = out_dir / "assemblies"
    reads_dir = out_dir / "reads"
    truth_dir = out_dir / "truth"
    assemblies_dir.mkdir(exist_ok=True)
    if spec.generate_reads:
        reads_dir.mkdir(exist_ok=True)
    if spec.write_truth_tables:
        truth_dir.mkdir(exist_ok=True)

    genome_records: List[Dict[str, Any]] = []
    assembly_sets: List[set] = []
    read_sets_min2: List[set] = []
    assembly_files: List[str] = []

    for g in range(spec.n_genomes):
        label = f"syn_g{g:04d}"
        contigs: List[Tuple[str, str]] = []
        templates: List[Tuple[str, str, int]] = []
        chrom_kmers_g: set = set()

        for c_idx, anc in enumerate(core_ancestors):
            mutated = mutate_dna(anc, spec.snp_rate, rng)
            chrom_kmers_g.update(space.iter_canonical(mutated))
            seq_for_fasta = mutated
            if spec.include_n_runs and c_idx == 0 and len(seq_for_fasta) > 200:
                mid = len(seq_for_fasta) // 2
                seq_for_fasta = seq_for_fasta[:mid] + "NNNNN" + seq_for_fasta[mid + 5 :]
            if spec.include_lowercase and c_idx == len(core_ancestors) - 1 and len(seq_for_fasta) > 120:
                seq_for_fasta = seq_for_fasta[:60] + seq_for_fasta[60:110].lower() + seq_for_fasta[110:]
            cname = f"{label}_core_{c_idx}"
            contigs.append((cname, seq_for_fasta))
            templates.append((cname, mutated, 1))

        # Deterministic accessory presence pattern across genomes (variable support).
        acc_present: List[int] = []
        for a_idx, a_seq in enumerate(accessory_pool):
            # Spread across different cohort fractions so patterns vary cleanly.
            if ((g + a_idx) % (a_idx + 2)) == 0 or (a_idx == 0 and g < max(2, spec.n_genomes // 2)):
                acc_present.append(a_idx)
                aseq = mutate_dna(a_seq, spec.snp_rate * 0.5, rng)
                aname = f"{label}_acc_{a_idx}"
                contigs.append((aname, aseq))
                templates.append((aname, aseq, 1))

        # Plasmids: collapsed to 1 copy in assembly FASTA when collapse_plasmids_in_assembly=True,
        # while raw reads sample from true copy number (>1).
        plasmids_present: Dict[str, int] = {}
        for p_idx, p_seq in enumerate(plasmid_pool):
            # Half the cohort carries plasmid 0; every third genome carries plasmid 1, etc.
            carries = (g % (p_idx + 2)) == 0
            if carries:
                true_cn = int(spec.plasmid_copy_numbers[p_idx])
                plasmids_present[f"plasmid_{p_idx}"] = true_cn
                pname = f"{label}_plasmid_{p_idx}"
                if spec.collapse_plasmids_in_assembly:
                    contigs.append((pname, p_seq))
                else:
                    for copy_i in range(true_cn):
                        contigs.append((f"{pname}_copy{copy_i}", p_seq))
                templates.append((pname, p_seq, true_cn))

        fna_path = assemblies_dir / f"{label}.fna"
        _write_fasta(fna_path, contigs)
        assembly_files.append(f"assemblies/{label}.fna")

        asm_min1 = parity.count_file(fna_path, k=spec.k, min_abundance=1)
        asm_min2 = parity.count_file(fna_path, k=spec.k, min_abundance=2)
        assembly_sets.append(set(asm_min1.counts))

        if spec.write_truth_tables:
            _write_truth_counts(truth_dir / f"{label}.assembly.counts.tsv", asm_min1.counts, space)

        rec: Dict[str, Any] = {
            "genome": label,
            "fasta": f"assemblies/{label}.fna",
            "contigs": len(contigs),
            "assembly_bp": sum(len(s) for _, s in contigs),
            "accessory_genes": acc_present,
            "plasmids_true_copy_number": plasmids_present,
            "assembly_truth": {
                "total_kmers": asm_min1.total_kmers,
                "distinct_kmers_min1": asm_min1.distinct_kept,
                "distinct_kmers_min2": asm_min2.distinct_kept,
                "fnv1a64_min1": asm_min1.fnv1a64_digest,
                "fnv1a64_min2": asm_min2.fnv1a64_digest,
                "sha256_canonical_min1": asm_min1.sha256_canonical,
                "sha256_spectrum_min1": asm_min1.sha256_spectrum,
            },
        }

        if spec.generate_reads:
            ext = ".fastq.gz" if spec.gzip_reads else ".fastq"
            r1_path = reads_dir / f"{label}_R1{ext}"
            r2_path = reads_dir / f"{label}_R2{ext}"
            r1_reads, r2_reads, n_errors = _simulate_paired_reads(
                templates=templates,
                coverage=spec.read_coverage,
                read_length=spec.read_length,
                insert_size=spec.insert_size,
                error_rate=spec.error_rate,
                rng=rng,
            )
            _write_fastq(r1_path, r1_reads, gzipped=spec.gzip_reads)
            _write_fastq(r2_path, r2_reads, gzipped=spec.gzip_reads)
            rd_min1 = parity.count_files([r1_path, r2_path], k=spec.k, min_abundance=1)
            rd_min2 = parity.count_files([r1_path, r2_path], k=spec.k, min_abundance=2)
            read_sets_min2.append(set(rd_min2.counts))

            if spec.write_truth_tables:
                _write_truth_counts(truth_dir / f"{label}.reads.counts.tsv", rd_min1.counts, space)

            chrom_depths = [rd_min1.counts[km] for km in chrom_kmers_g if km in rd_min1.counts]
            chrom_med = float(statistics.median(chrom_depths)) if chrom_depths else 0.0
            plasmid_ratios: Dict[str, float] = {}
            for p_idx in range(spec.n_plasmids):
                pkey = f"plasmid_{p_idx}"
                if pkey in plasmids_present and chrom_med > 0.0:
                    p_only = plasmid_kmers[p_idx] - chrom_kmers_g
                    p_depths = [rd_min1.counts[km] for km in p_only if km in rd_min1.counts]
                    if p_depths:
                        plasmid_ratios[pkey] = round(float(statistics.median(p_depths)) / chrom_med, 4)

            rec["reads_r1"] = f"reads/{label}_R1{ext}"
            rec["reads_r2"] = f"reads/{label}_R2{ext}"
            rec["reads_truth"] = {
                "read_pairs": len(r1_reads),
                "injected_base_errors": n_errors,
                "total_kmers": rd_min1.total_kmers,
                "distinct_kmers_min1": rd_min1.distinct_kept,
                "distinct_kmers_min2": rd_min2.distinct_kept,
                "fnv1a64_min1": rd_min1.fnv1a64_digest,
                "fnv1a64_min2": rd_min2.fnv1a64_digest,
                "sha256_canonical_min2": rd_min2.sha256_canonical,
                "sha256_spectrum_min1": rd_min1.sha256_spectrum,
                "sha256_spectrum_min2": rd_min2.sha256_spectrum,
                "chrom_median_depth": round(chrom_med, 4),
                "plasmid_observed_depth_ratio": plasmid_ratios,
            }

        genome_records.append(rec)

    # Write FOF files (relative paths) for external multi-genome tools.
    (out_dir / "assemblies.fof").write_text("\n".join(assembly_files) + "\n", encoding="ascii")
    if spec.generate_reads:
        fof_lines = [f"{r['genome']} : {r['reads_r1']} ; {r['reads_r2']}" for r in genome_records]
        (out_dir / "reads.kmtricks.fof").write_text("\n".join(fof_lines) + "\n", encoding="ascii")

    ref_matrix = parity.build_reference_matrix(assembly_sets, k=spec.k, min_support=spec.min_support)
    cohort_truth: Dict[str, Any] = {
        "assembly_matrix": {
            "min_support": ref_matrix.min_support,
            "kmers_distinct": ref_matrix.kmers_distinct,
            "kmers_kept": ref_matrix.kmers_kept,
            "kmers_rare": ref_matrix.kmers_rare,
            "kmers_core": ref_matrix.kmers_core,
            "n_patterns": ref_matrix.n_patterns,
            "digest_sha256": ref_matrix.digest_sha256,
        },
    }
    if spec.generate_reads:
        rd_matrix = parity.build_reference_matrix(read_sets_min2, k=spec.k, min_support=spec.min_support)
        cohort_truth["reads_matrix_min2"] = {
            "min_support": rd_matrix.min_support,
            "kmers_distinct": rd_matrix.kmers_distinct,
            "kmers_kept": rd_matrix.kmers_kept,
            "kmers_rare": rd_matrix.kmers_rare,
            "kmers_core": rd_matrix.kmers_core,
            "n_patterns": rd_matrix.n_patterns,
            "digest_sha256": rd_matrix.digest_sha256,
        }

    if spec.compute_unitig_truth:
        union_kmers: set = set().union(*assembly_sets) if assembly_sets else set()
        unitigs = parity.compact_unitigs(union_kmers, k=spec.k)
        cohort_truth["unitigs"] = {
            "graph_kmers": len(union_kmers),
            "n_unitigs": len(unitigs),
        }

    manifest: Dict[str, Any] = {
        "format": COHORT_FORMAT,
        "spec": asdict(spec),
        "cohort_truth": cohort_truth,
        "genomes": genome_records,
    }
    manifest_path = out_dir / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest
