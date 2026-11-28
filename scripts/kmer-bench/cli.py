#!/usr/bin/env python3
"""CLI for the GRM k-mer counting and multi-genome matrix benchmark suite.

Subcommands:
  generate       Create a deterministic synthetic bacterial cohort (FASTA + FASTQ + ground truth).
  verify-parity  Check exact canonical k-mer set, abundance, or unitig parity against reference files.
  plan           Emit a deterministic benchmark command plan across workloads T1..T5.
  run            Execute benchmarks (--dry-run, --mock self-check, or real) and write JSON + TSV.
  summarize      Aggregate bench.json files into comparison TSV and Markdown tables.

Pass output directories outside the repository checkout (for example under /tmp/kmer-bench)
so benchmark artifacts never leave untracked files in the git tree.
"""

import argparse
from dataclasses import asdict
import importlib.util
import json
import pathlib
import sys
from typing import Any, Dict, List, Optional, Sequence


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
synthetic = _load_sibling("synthetic")
runners = _load_sibling("runners")


def _parse_csv(raw: str, allowed: Sequence[str], label: str) -> List[str]:
    if raw.strip().lower() == "all":
        return list(allowed)
    items = [x.strip() for x in raw.split(",") if x.strip()]
    bad = [x for x in items if x not in allowed]
    if bad:
        raise ValueError(f"unknown {label}: {', '.join(bad)}")
    return items


def _cmd_generate(args: argparse.Namespace) -> int:
    copy_nums = tuple(int(x) for x in args.plasmid_copies.split(",") if x.strip())
    spec = synthetic.CohortSpec(
        seed=args.seed,
        n_genomes=args.n_genomes,
        k=args.k,
        min_support=args.min_support,
        core_len=args.core_len,
        n_core_contigs=args.n_core_contigs,
        n_accessory_genes=args.n_accessory_genes,
        accessory_len=args.accessory_len,
        snp_rate=args.snp_rate,
        n_plasmids=args.n_plasmids,
        plasmid_len=args.plasmid_len,
        plasmid_copy_numbers=copy_nums,
        collapse_plasmids_in_assembly=not args.keep_plasmid_copies_in_fasta,
        generate_reads=not args.no_reads,
        read_coverage=args.coverage,
        read_length=args.read_length,
        insert_size=args.insert_size,
        error_rate=args.error_rate,
        gzip_reads=args.gzip_reads,
        write_truth_tables=not args.no_truth_tables,
        compute_unitig_truth=not args.no_unitig_truth,
    )
    out_path = pathlib.Path(args.out).resolve()
    manifest = synthetic.generate_cohort(out_path, spec)
    print(
        f"generated {manifest['spec']['n_genomes']} synthetic genomes (k={manifest['spec']['k']}) in {out_path}",
        flush=True,
    )
    return 0


def _cmd_verify_parity(args: argparse.Namespace) -> int:
    ref_paths = [pathlib.Path(p) for p in args.reference]
    k = args.k
    min_ab = args.min_abundance
    result: Dict[str, Any]

    if args.mode == "set":
        ref = parity.count_files(ref_paths, k=k, min_abundance=min_ab)
        if args.candidate_kset:
            kset = parity.read_kset(pathlib.Path(args.candidate_kset))
            if kset.k != k:
                raise ValueError(f"kset k={kset.k} does not match expected k={k}")
            space = parity.KmerSpace(k)
            obs_canonical = [space.unmix(v) for v in kset.values]
        elif args.candidate_dump:
            obs_canonical = list(parity.parse_kmer_text_dump(pathlib.Path(args.candidate_dump), k=k, min_abundance=min_ab).keys())
        else:
            raise ValueError("set parity requires --candidate-kset or --candidate-dump")
        rep = parity.verify_set_parity(ref.counts.keys(), obs_canonical, k=k)
        result = {"mode": "set", "ok": rep.ok, "report": asdict(rep)}

    elif args.mode == "abundance":
        if not args.candidate_dump:
            raise ValueError("abundance parity requires --candidate-dump")
        ref = parity.count_files(ref_paths, k=k, min_abundance=min_ab)
        obs = parity.parse_kmer_text_dump(pathlib.Path(args.candidate_dump), k=k, min_abundance=min_ab)
        rep = parity.verify_abundance_parity(
            ref.counts, obs, k=k, min_abundance=min_ab, max_cap=args.max_cap
        )
        result = {"mode": "abundance", "ok": rep.ok, "report": asdict(rep)}

    elif args.mode == "unitigs":
        if not args.candidate_unitigs:
            raise ValueError("unitigs parity requires --candidate-unitigs")
        union_kmers: set = set()
        for p in ref_paths:
            union_kmers.update(parity.count_file(p, k=k, min_abundance=min_ab).counts)
        ref_u = parity.compact_unitigs(union_kmers, k=k)
        obs_u = parity.parse_unitigs_fasta(pathlib.Path(args.candidate_unitigs), k=k)
        rep = parity.verify_unitig_parity(ref_u, obs_u)
        result = {"mode": "unitigs", "ok": rep.ok, "report": asdict(rep)}

    else:
        raise ValueError(f"unsupported parity mode: {args.mode!r}")

    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        out_path = pathlib.Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0 if result["ok"] else 1


def _build_plans(args: argparse.Namespace, cohort_manifest: Optional[Dict[str, Any]] = None) -> List[Any]:
    workloads = _parse_csv(args.workloads, list(runners.WORKLOADS.keys()), "workloads")
    tools = _parse_csv(args.tools, list(runners.TOOL_REGISTRY.keys()), "tools")
    cohort_dir = pathlib.Path(args.cohort_dir).resolve()
    out_dir = pathlib.Path(args.out_dir).resolve()
    min_support = int(cohort_manifest["spec"].get("min_support", 2)) if cohort_manifest else 2
    plans = []
    for wl in workloads:
        for tool in tools:
            if wl not in runners.TOOL_REGISTRY[tool]["workloads"]:
                continue
            variants = (
                ("default", "column-blocks-64", "disk-sparse-256")
                if (tool == "grm" and wl == "T3" and args.include_grm_variants)
                else ("default",)
            )
            bin_override = args.grm_bin if (tool == "grm" and args.grm_bin) else None
            for variant in variants:
                plans.append(
                    runners.build_command_plan(
                        tool=tool,
                        workload=wl,
                        cohort_dir=cohort_dir,
                        run_dir=out_dir / "work",
                        k=args.k,
                        min_support=min_support,
                        threads=args.threads,
                        mem_gb=args.mem_gb,
                        variant=variant,
                        binary_override=bin_override,
                        cohort_manifest=cohort_manifest,
                    )
                )
    return plans


def _cmd_plan(args: argparse.Namespace) -> int:
    plans = _build_plans(args)
    payload = {
        "format": "grm-kmer-bench-plan-v1",
        "cases": [asdict(p) for p in plans],
    }
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.json_out:
        out_path = pathlib.Path(args.json_out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
    else:
        print(text, end="")
    return 0


def _cmd_run(args: argparse.Namespace) -> int:
    cohort_dir = pathlib.Path(args.cohort_dir).resolve()
    out_dir = pathlib.Path(args.out_dir).resolve()
    manifest_path = cohort_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if args.dry_run and args.mock:
        raise ValueError("cannot combine --dry-run and --mock")
    mode = "dry-run" if args.dry_run else ("mock" if args.mock else "real")
    plans = _build_plans(args, cohort_manifest=manifest)
    results = [
        runners.execute_case(
            plan=p,
            cohort_manifest=manifest,
            cohort_dir=cohort_dir,
            mode=mode,
            skip_missing=not args.fail_missing,
        )
        for p in plans
    ]
    report = runners.build_bench_report(results, manifest, strict_hygiene=not args.no_hygiene_check)
    json_path, tsv_path = runners.write_bench_outputs(
        out_dir, report, strict_hygiene=not args.no_hygiene_check
    )
    print(f"wrote {json_path} and {tsv_path} ({len(results)} cases, mode={mode})", flush=True)
    failed = [r for r in results if r["status"].startswith("failed") or r["status"] == "parity_failed"]
    return 1 if failed else 0


def _cmd_summarize(args: argparse.Namespace) -> int:
    combined: List[Dict[str, Any]] = []
    for p_str in args.inputs:
        data = json.loads(pathlib.Path(p_str).read_text(encoding="utf-8"))
        if data.get("format") != runners.BENCH_FORMAT:
            raise ValueError(f"{p_str}: unexpected format {data.get('format')!r}")
        combined.extend(data["results"])
    tsv_text = runners.render_summary_tsv(
        combined,
        include_mock=args.include_mock,
        include_not_comparable=args.include_not_comparable,
    )
    md_text = runners.render_summary_markdown(
        combined,
        include_mock=args.include_mock,
        include_not_comparable=args.include_not_comparable,
    )
    tsv_leaks = runners.scan_tsv_for_leaks(tsv_text)
    if tsv_leaks:
        raise ValueError(f"summary tsv failed hygiene scan: {', '.join(tsv_leaks)}")
    md_leaks = runners.scan_for_leaks(md_text)
    if md_leaks:
        raise ValueError(f"summary md failed hygiene scan: {', '.join(md_leaks)}")
    if args.tsv_out:
        p = pathlib.Path(args.tsv_out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(tsv_text, encoding="utf-8")
    if args.md_out:
        p = pathlib.Path(args.md_out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(md_text, encoding="utf-8")
    if not args.tsv_out and not args.md_out:
        print(md_text, end="")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="kmer-bench",
        description=(
            "Reproducible k-mer counting and multi-genome matrix benchmark harness. "
            "Use output directories outside the git tree (e.g., /tmp/kmer-bench)."
        ),
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p_gen = sub.add_parser("generate", help="Generate a deterministic synthetic bacterial cohort.")
    p_gen.add_argument(
        "--out",
        required=True,
        help="Output cohort directory outside the repository tree (e.g., /tmp/kmer-bench/cohort).",
    )
    p_gen.add_argument("--seed", type=int, default=42)
    p_gen.add_argument("--n-genomes", type=int, default=8)
    p_gen.add_argument("--k", type=int, default=31)
    p_gen.add_argument("--min-support", type=int, default=2)
    p_gen.add_argument("--core-len", type=int, default=4800)
    p_gen.add_argument("--n-core-contigs", type=int, default=3)
    p_gen.add_argument("--n-accessory-genes", type=int, default=6)
    p_gen.add_argument("--accessory-len", type=int, default=600)
    p_gen.add_argument("--snp-rate", type=float, default=0.005)
    p_gen.add_argument("--n-plasmids", type=int, default=2)
    p_gen.add_argument("--plasmid-len", type=int, default=1000)
    p_gen.add_argument("--plasmid-copies", default="4,8", help="Comma-separated true plasmid copy numbers.")
    p_gen.add_argument("--keep-plasmid-copies-in-fasta", action="store_true")
    p_gen.add_argument("--no-reads", action="store_true")
    p_gen.add_argument("--coverage", type=float, default=12.0)
    p_gen.add_argument("--read-length", type=int, default=120)
    p_gen.add_argument("--insert-size", type=int, default=220)
    p_gen.add_argument("--error-rate", type=float, default=0.005)
    p_gen.add_argument("--gzip-reads", action="store_true")
    p_gen.add_argument("--no-truth-tables", action="store_true")
    p_gen.add_argument("--no-unitig-truth", action="store_true")
    p_gen.set_defaults(func=_cmd_generate)

    p_par = sub.add_parser("verify-parity", help="Verify exact k-mer set, abundance, or unitig parity.")
    p_par.add_argument("--mode", choices=("set", "abundance", "unitigs"), required=True)
    p_par.add_argument("--k", type=int, default=31)
    p_par.add_argument("--min-abundance", type=int, default=1)
    p_par.add_argument("--max-cap", type=int, default=None)
    p_par.add_argument("--reference", nargs="+", required=True, help="Reference FASTA/FASTQ file(s).")
    p_par.add_argument("--candidate-kset", default=None)
    p_par.add_argument("--candidate-dump", default=None)
    p_par.add_argument("--candidate-unitigs", default=None)
    p_par.add_argument("--json-out", default=None)
    p_par.set_defaults(func=_cmd_verify_parity)

    def _add_plan_args(p: argparse.ArgumentParser) -> None:
        p.add_argument(
            "--cohort-dir",
            required=True,
            help="Synthetic cohort directory containing manifest.json (e.g., /tmp/kmer-bench/cohort).",
        )
        p.add_argument(
            "--out-dir",
            required=True,
            help="Benchmark output directory outside the repository tree (e.g., /tmp/kmer-bench/out).",
        )
        p.add_argument("--workloads", default="all", help="Comma list of T1..T5 or 'all'.")
        p.add_argument("--tools", default="all", help="Comma list of tool names or 'all'.")
        p.add_argument("--k", type=int, default=31)
        p.add_argument("--threads", type=int, default=8)
        p.add_argument("--mem-gb", type=int, default=16)
        p.add_argument("--grm-bin", default=None, help="Optional path to grm binary override.")
        p.add_argument("--include-grm-variants", action="store_true")

    p_plan = sub.add_parser("plan", help="Emit benchmark command plan JSON.")
    _add_plan_args(p_plan)
    p_plan.add_argument("--json-out", default=None)
    p_plan.set_defaults(func=_cmd_plan)

    p_run = sub.add_parser("run", help="Run benchmarks and write grm-kmer-bench-v1 JSON + TSV.")
    _add_plan_args(p_run)
    p_run.add_argument("--dry-run", action="store_true")
    p_run.add_argument(
        "--mock",
        action="store_true",
        help="Run parity.py self-check with synthetic placeholder timings/sizes (excluded from real summaries).",
    )
    p_run.add_argument("--fail-missing", action="store_true", help="Fail instead of skipping absent external binaries.")
    p_run.add_argument("--no-hygiene-check", action="store_true")
    p_run.set_defaults(func=_cmd_run)

    p_sum = sub.add_parser(
        "summarize",
        help="Aggregate grm-kmer-bench-v1 JSON files into comparison TSV/Markdown (real comparable rows by default).",
    )
    p_sum.add_argument("--inputs", nargs="+", required=True)
    p_sum.add_argument("--tsv-out", default=None)
    p_sum.add_argument("--md-out", default=None)
    p_sum.add_argument(
        "--include-mock",
        action="store_true",
        help="Include mock_self_check and dry_run rows (excluded from real comparison summaries by default).",
    )
    p_sum.add_argument(
        "--include-not-comparable",
        action="store_true",
        help="Include not_comparable rows (excluded from real comparison summaries by default).",
    )
    p_sum.set_defaults(func=_cmd_summarize)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    sys.exit(main())
