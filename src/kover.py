from dataclasses import dataclass
import os
import csv
import pathlib as pl
from typing import Annotated, Collection, Optional
from util import to_linux_path


class Source(str):
    READS = "reads"
    CONTIGS = "contigs"
    K_MER_MATREX = "tsv"


class ModelType(str):
    CONJUNCTION = "conjunction"
    DISJUNCTION = "disjunction"
    BOTH = "conjunction disjunction"


class HpChoice(str):
    BOUND = "bound"
    CV = "cv"
    NONE = "none"


class Criterion(str):
    GINI = "gini"
    CROSS_ENTROPY = "crossentropy"


@dataclass
class MinLen:
    value: int


DEFAULT = None


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
    phenotype_description: Optional[str] = DEFAULT,
    phenotype_metadata: Optional[str] = DEFAULT,
    kmer_size: Optional[int | str] = 31,
    kmer_min_abundance: Optional[int | str] = 1,
    singleton_kmers: bool = DEFAULT,
    n_cpu: Optional[int | str] = 0,
    compression: Optional[int | str] = 4,
    temp_dir: Optional[str] = DEFAULT,
    x: bool = False,
    v: bool = False,
):
    if temp_dir == "":
        temp_dir = None

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
        f"--kmer-min-abundance {kmer_min_abundance}"
        if source == Source.READS and kmer_min_abundance
        else "",
        "--singleton-kmers" if singleton_kmers else "",
        f"--n-cpu {n_cpu}" if source != Source.K_MER_MATREX and n_cpu else "",
        f"--compression {compression}" if compression else "",
        f"--temp-dir {to_linux_path(temp_dir)}"
        if source != Source.K_MER_MATREX and temp_dir
        else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def split_command(
    kover_path: str,
    dataset: str,
    id: str,
    train_size: Optional[str | float] = 0.5,
    train_ids: Optional[str] = DEFAULT,
    test_ids: Optional[str] = DEFAULT,
    folds: int | str = 0,
    random_seed: Optional[int | str] = DEFAULT,
    x: bool = False,
    v: bool = False,
):
    command = (
        to_linux_path(kover_path),
        "dataset split",
        f"--dataset {to_linux_path(dataset)}",
        f"--id {id}",
        f"--train-size {train_size}" if train_size else "",
        f"--train-ids {to_linux_path(train_ids)}" if train_ids else "",
        f"--test-ids {to_linux_path(test_ids)}" if test_ids else "",
        f"--folds {folds}",
        f"--random-seed {random_seed}" if random_seed else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


def info_command(
    kover_path: str,
    dataset: str,
    a: bool = False,
    genome_type: bool = False,
    genome_source: bool = False,
    genome_ids: bool = False,
    genome_count: bool = False,
    kmers: bool = False,
    kmer_len: bool = False,
    kmer_count: bool = False,
    phenotype_description: bool = False,
    phenotype_metadata: bool = False,
    phenotype_tags: bool = False,
    splits: bool = False,
    uuid: bool = False,
    compression: bool = False,
    classification_type: bool = False,
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


DEFAULT_P = (
    0.1,
    0.178,
    0.316,
    0.562,
    1.0,
    1.778,
    3.162,
    5.623,
    10.0,
    999999.0,
)


def scm_command(
    kover_path: str,
    dataset: str,
    split: str,
    model_type: ModelType,
    p: Annotated[Collection[float | str], MinLen(1)] = DEFAULT_P,
    max_rules: Optional[str | int] = DEFAULT,
    max_equiv_rules: Optional[str | int] = DEFAULT,
    kmer_blacklist: Optional[str] = DEFAULT,
    hp_choice: HpChoice = HpChoice.CV,
    bound_max_genome_size: str | int = DEFAULT,
    random_seed: Optional[str | int] = DEFAULT,
    n_cpu: str | int = DEFAULT,
    output_dir: Optional[str] = DEFAULT,
    x: bool = False,
    v: bool = False,
):
    if not p:
        p = DEFAULT_P
    if max_rules == 0:
        max_rules = None
    if max_equiv_rules == 0:
        max_equiv_rules = None

    command = (
        to_linux_path(kover_path),
        "learn scm",
        f"--dataset {to_linux_path(dataset)}",
        f"--split {split}",
        f"--model-type {model_type}",
        f"--p {' '.join([str(n) for n in p])}",
        f"--kmer-blacklist {to_linux_path(kmer_blacklist)}" if kmer_blacklist else "",
        f"--max-rules {max_rules}" if max_rules else "",
        f"--max-equiv-rules {max_equiv_rules}" if max_equiv_rules else "",
        f"--hp-choice {hp_choice}",
        f"--bound-max-genome-size {bound_max_genome_size}"
        if hp_choice == HpChoice.BOUND and bound_max_genome_size
        else "",
        f"--random-seed {random_seed}" if random_seed else "",
        f"--n-cpu {n_cpu}" if n_cpu else "",
        f"--output-dir {to_linux_path(output_dir)}" if output_dir else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))


DEFAULT_CLASS_IMPORTANCE = (0.25, 0.5, 0.75, 1.0)


def tree_command(
    kover_path: str,
    dataset: str,
    split: str,
    criterion: Criterion = Criterion.GINI,
    max_depth: str | int = 10,
    min_samples_split: str | int = 2,
    class_importance: Collection[str | int] = DEFAULT_CLASS_IMPORTANCE,
    kmer_blacklist: Optional[str] = DEFAULT,
    hp_choice: HpChoice = HpChoice.CV,
    bound_max_genome_size: Optional[str | int] = DEFAULT,
    n_cpu: Optional[str | int] = DEFAULT,
    output_dir: Optional[str] = DEFAULT,
    x: bool = False,
    v: bool = False,
):
    if not class_importance:
        class_importance = DEFAULT_CLASS_IMPORTANCE

    command = (
        to_linux_path(kover_path),
        "learn tree",
        f"--dataset {to_linux_path(dataset)}",
        f"--split {split}",
        f"--criterion {criterion}",
        f"--max-depth {max_depth}",
        f"--min-samples-split {min_samples_split}",
        f"--class-importance {' '.join([str(n) for n in class_importance])}",
        f"--kmer-blacklist {to_linux_path(kmer_blacklist)}" if kmer_blacklist else "",
        f"--hp-choice {hp_choice}",
        f"--bound-max-genome-size {bound_max_genome_size}"
        if bound_max_genome_size
        else "",
        f"--n-cpu {n_cpu}",
        f"--output-dir {to_linux_path(output_dir)}" if output_dir else "",
        "-x" if x else "",
        "-v" if v else "",
    )

    return " ".join(filter(lambda x: x != "", command))
