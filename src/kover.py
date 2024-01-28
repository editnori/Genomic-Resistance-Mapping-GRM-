import os
import csv
import pathlib as pl
from typing import Optional
from util import to_linux_path


class Source(str):
    READS = "reads"
    CONTIGS = "contigs"
    K_MER_MATREX = "tsv"


def create_contigs_path_tsv(contigs_path: str, genome_name: str):
    genome_contigs_path = os.path.join(contigs_path, genome_name)

    output_file = genome_contigs_path + "_paths.tsv"

    with open(output_file, "w", newline="", encoding="utf-8") as tsv_file:
        writer = csv.writer(tsv_file, delimiter="\t", lineterminator="\n")

        for file_path in pl.Path(genome_contigs_path).glob("*.fna"):
            writer.writerow([file_path.stem, to_linux_path(file_path)])


def create_command(
    kover_path: str,
    source: Source,
    genomic_data: str,
    output: str,
    phenotype_description: Optional[str] = None,
    phenotype_metadata: Optional[str] = None,
    kmer_size: Optional[int | str] = None,
    singleton_kmers: Optional[bool] = False,
    n_cpu: Optional[int | str] = None,
    compression=None,
    temp_dir: Optional[str] = None,
    x: Optional[bool] = False,
    v: Optional[bool] = False,
):
    command = (
        to_linux_path(kover_path),
        f"dataset create from-{source}",
        f"--genomic-data {to_linux_path(genomic_data)}",
        f"--phenotype-description {to_linux_path(phenotype_description)}"
        if phenotype_description
        else "",
        f"--phenotype-metadata {to_linux_path(phenotype_metadata)}"
        if phenotype_metadata
        else "",
        f"--output {to_linux_path(output)}",
        f"--kmer-size {kmer_size}"
        if source != Source.K_MER_MATREX and kmer_size
        else "",
        "--singleton-kmers" if singleton_kmers else "",
        f"--n-cpu {n_cpu}" if source != Source.K_MER_MATREX and n_cpu else "",
        f"--compression {compression}" if compression else "",
        f"--temp-dir {temp_dir}" if source != Source.K_MER_MATREX and temp_dir else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def split_dataset(
    kover_path: str,
    dataset: str,
    id: str,
    train_size: Optional[str] = None,
    train_ids: Optional[str] = None,
    folds: Optional[int | str] = None,
    random_seed: Optional[int | str] = None,
    x: Optional[bool] = False,
    v: Optional[bool] = False,
):
    command = (
        to_linux_path(kover_path),
        "dataset split",
        f"--dataset {to_linux_path(dataset)}",
        f"--id {id}",
        f"--train-size {train_size}" if train_size else "",
        f"--train-ids {train_ids}" if train_ids else "",
        f"--folds {folds}" if folds else "",
        f"--random-seed {random_seed}" if random_seed else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def info_dataset(
    kover_path: str,
    dataset: str,
    a: Optional[bool] = False,
    genome_type: Optional[bool] = False,
    genome_source: Optional[bool] = False,
    genome_ids: Optional[bool] = False,
    genome_count: Optional[bool] = False,
    kmers: Optional[bool] = False,
    kmer_len: Optional[bool] = False,
    kmer_count: Optional[bool] = False,
    phenotype_description: Optional[bool] = False,
    phenotype_metadata: Optional[bool] = False,
    phenotype_tags: Optional[bool] = False,
    splits: Optional[bool] = False,
    uuid: Optional[bool] = False,
    compression: Optional[bool] = False,
    classification_type: Optional[bool] = False,
):
    command = (
        to_linux_path(kover_path),
        "dataset info",
        f"--dataset {to_linux_path(dataset)}",
        "--all" if a else "",
        "--genome-type" if genome_type else "",
        "--genome-source" if genome_source else "",
        "--genome-ids" if genome_ids else "",
        "--genome-count" if genome_count else "",
        "--kmers" if kmers else "",
        "--kmer-len" if kmer_len else "",
        "--kmer-count" if kmer_count else "",
        "--phenotype-description" if phenotype_description else "",
        "--phenotype-metadata" if phenotype_metadata else "",
        "--phenotype-tags" if phenotype_tags else "",
        "--splits" if splits else "",
        "--uuid" if uuid else "",
        "--compression" if compression else "",
        "--classification-type" if classification_type else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def scm_command(
    kover_path: str,
    dataset: str,
    split: str,
    model_type: Optional[str] = None,
    p: Optional[str] = None,
    kmer_blacklist: Optional[str] = None,
    max_rules: Optional[str | int] = None,
    max_equiv_rules: Optional[str | int] = None,
    hp_choice: Optional[str] = None,
    bound_max_genome_size: Optional[str | int] = None,
    random_seed: Optional[str | int] = None,
    n_cpu: Optional[str | int] = None,
    output_dir: Optional[str] = None,
    x: Optional[bool] = False,
    v: Optional[bool] = False,
):
    command = (
        to_linux_path(kover_path),
        "learn scm",
        f"--dataset {dataset}",
        f"--split {split}",
        # [--model-type {conjunction,disjunction} [{conjunction,disjunction} ...]] if model_type else "",
        # [--p P [P ...]] if p else "",
        f"--kmer-blacklist {kmer_blacklist}" if kmer_blacklist else "",
        f"--max-rules {max_rules}" if max_rules else "",
        f"--max-equiv-rules {max_equiv_rules}" if max_equiv_rules else "",
        # f"--hp-choice {bound,cv,none}" if hp_choice else "",
        f"--bound-max-genome-size {bound_max_genome_size}"
        if bound_max_genome_size
        else "",
        f"--random-seed {random_seed}" if random_seed else "",
        f"--n-cpu {n_cpu}" if n_cpu else "",
        f"--output-dir {output_dir}" if output_dir else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def tree_command(
    kover_path: str,
    dataset: Optional[str] = None,
    split: Optional[str] = None,
    criterion: Optional[str] = None,
    max_depth: Optional[str] = None,
    min_samples_split: Optional[str] = None,
    class_importance: Optional[str] = None,
    kmer_blacklist: Optional[str] = None,
    hp_choice: Optional[str] = None,
    bound_max_genome_size: Optional[str | int] = None,
    n_cpu: Optional[str | int] = None,
    output_dir: Optional[str] = None,
    x: Optional[bool] = False,
    v: Optional[bool] = False,
):
    command = (
        to_linux_path(kover_path),
        "learn tree",
        f"--dataset {dataset}",
        f"--split {split}",
        # f"--criterion {gini,crossentropy} [{gini,crossentropy} ...]" if criterion else "",
        # f"--max-depth MAX_DEPTH [MAX_DEPTH ...]" if max_depth else "",
        # f"--min-samples-split MIN_SAMPLES_SPLIT [MIN_SAMPLES_SPLIT ...]" if min_samples_split else "",
        # f"--class-importance CLASS_IMPORTANCE [CLASS_IMPORTANCE ...]" if class_importance else "",
        f"--kmer-blacklist {kmer_blacklist}" if kmer_blacklist else "",
        # "--hp-choice {bound,cv}" if hp_choice else "",
        f"--bound-max-genome-size {bound_max_genome_size}"
        if bound_max_genome_size
        else "",
        f"--n-cpu {n_cpu}" if n_cpu else "",
        f"--output-dir {output_dir}" if output_dir else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))
