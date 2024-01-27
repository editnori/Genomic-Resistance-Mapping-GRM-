import subprocess
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
    phenotype_description: Optional[int] = None,
    phenotype_metadata: Optional[int] = None,
    kmer_size: Optional[int | str] = None,
    singleton_kmers: Optional[bool] = None,
    n_cpu: Optional[int] = None,
    compression=None,
    temp_dir: Optional[str] = None,
    x: Optional[bool] = None,
    v: Optional[bool] = None,
):
    command = (
        kover_path,
        f"dataset create from-{source}",
        f"--genomic-data {genomic_data}" if genomic_data else "",
        f"--phenotype-description {phenotype_description}"
        if phenotype_description
        else "",
        f"--phenotype-metadata {phenotype_metadata}" if phenotype_metadata else "",
        f"--output {output}",
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

    return " ".join(command)


class KoverDatasetCreator:
    def split_dataset(
        self,
        dataset,
        id,
        train_size=None,
        folds=2,
        random_seed=None,
        output=None,
        v=False,
        x=True,
    ):
        command = [
            "wsl",
            "/home/mhdeeb/kover/bin/kover",
            "dataset",
            "split",
            "--dataset",
            dataset,
            "--id",
            id,
            "--train-size",
            str(train_size) if train_size else "0.5",
            "--folds",
            folds,
            "--random-seed",
            str(random_seed) if random_seed else "",
        ]

        if x:
            command.append("-x")
        if v:
            command.append("-v")

        return command

    def dataset_info(
        self,
        dataset,
        all=False,
        genome_type=False,
        genome_source=False,
        genome_ids=False,
        genome_count=False,
        kmers=False,
        kmer_len=False,
        kmer_count=False,
        phenotype_description=False,
        phenotype_metadata=False,
        phenotype_tags=False,
        splits=False,
        uuid=False,
        compression=False,
        classification_type=False,
    ):
        command = [
            "wsl",
            "kover",
            "dataset",
            "info",
            "--dataset",
            dataset,
        ]

        if all:
            command.append("--all")
        if genome_type:
            command.append("--genome-type")
        if genome_source:
            command.append("--genome-source")
        if genome_ids:
            command.append("--genome-ids")
        if genome_count:
            command.append("--genome-count")
        if kmers:
            command.append("--kmers")
        if kmer_len:
            command.append("--kmer-len")
        if kmer_count:
            command.append("--kmer-count")
        if phenotype_description:
            command.append("--phenotype-description")
        if phenotype_metadata:
            command.append("--phenotype-metadata")
        if phenotype_tags:
            command.append("--phenotype-tags")
        if splits:
            command.append("--splits")
        if uuid:
            command.append("--uuid")
        if compression:
            command.append("--compression")
        if classification_type:
            command.append("--classification-type")

        subprocess.run(command)

    def kover_learn_scm(
        self,
        dataset,
        splitid,
        model_type,
        hyperparameter,
        max_rules,
        maxequivrules,
        hpchoice,
        bound_max=None,
        random_seed=None,
        output=None,
        x=True,
    ):
        command = [
            "wsl",
            "/home/mhdeeb/kover/bin/kover",
            "learn",
            "scm",
            "--dataset",
            dataset,
            "--split",
            splitid,
            "--model-type",
            str(model_type),
        ]

        # Extend the command list with --p options for each value
        command.extend(
            ["--p"] + [str(value).replace(",", "") for value in hyperparameter]
        )

        command.extend(
            [
                "--max-rules",
                max_rules,
                "--max-equiv-rules",
                maxequivrules,
                "--hp-choice",
                hpchoice,
            ]
        )

        # Include --bound-max-genome-size only if hpchoice is "bound"
        if hpchoice == "bound" and bound_max is not None:
            command.extend(["--bound-max-genome-size", str(bound_max)])

        command.extend(
            [
                "--random-seed",
                str(random_seed) if random_seed else "",
                "--output-dir",
                output,
            ]
        )

        if x:
            command.append("-x")

        command_string = " ".join(map(str, command))

        return command_string

    def kover_learn_cart(
        self,
        dataset,
        splitid,
        model_type,
        hyperparameter,
        max_rules,
        maxequivrules,
        hpchoice,
        bound_max=None,
        output=None,
        x=True,
    ):
        command = [
            "wsl",
            "/home/mhdeeb/kover/bin/kover",
            "learn",
            "tree",
            "--dataset",
            dataset,
            "--split",
            splitid,
            "--criterion",
            str(model_type),
            "--max-depth",
            max_rules,
            "--min-samples-split",
            maxequivrules,
        ]

        # Extend the command list with --p options for each value
        command.extend(
            ["--class-importance"]
            + [str(value).replace(",", "") for value in hyperparameter]
        )

        command.extend(
            [
                "--hp-choice",
                hpchoice,
            ]
        )

        # Include --bound-max-genome-size only if hpchoice is "bound"
        if hpchoice == "bound" and bound_max is not None:
            command.extend(["--bound-max-genome-size", str(bound_max)])

        command.extend(
            [
                "--output-dir",
                output,
            ]
        )

        if x:
            command.append("-x")

        command_string = " ".join(map(str, command))

        return command_string
