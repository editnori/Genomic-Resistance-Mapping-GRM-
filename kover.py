import subprocess
import os
import csv


class KoverDatasetCreator:
    def contigs_parser(self, dataset_path, output_tsv):
        # Implement method logic for parsing contigs and creating a TSV file
        with open(output_tsv, "w", newline="", encoding="utf-8") as tsv_file:
            writer = csv.writer(tsv_file, delimiter="\t")

            for root, dirs, files in os.walk(dataset_path):
                for filename in files:
                    if filename.endswith(".fna"):
                        genome_id = os.path.splitext(filename)[0]
                        fasta_path = os.path.abspath(os.path.join(root, filename))

                        # Replace backslashes with forward slashes
                        fasta_path = fasta_path.replace("\\", "/")

                        # Convert Windows path to WSL-compatible path
                        if fasta_path.startswith("C:"):
                            fasta_path = "/mnt/c" + fasta_path[2:]

                        writer.writerow([genome_id, fasta_path])

    def create_from_reads(
        self,
        genomic_data,
        phenotype_description=None,
        phenotype_metadata=None,
        output=None,
        kmer_size=None,
        kmer_min_abundance=None,
        singleton_kmers=False,
        n_cpu=None,
        compression=None,
        temp_dir=None,
        x=True,
        v=False,
    ):
        command = [
            "wsl",
            "kover",
            "dataset",
            "create",
            "from-reads",
            "--genomic-data",
            genomic_data,
            "--phenotype-description",
            phenotype_description,
            "--phenotype-metadata",
            phenotype_metadata,
            "--output",
            output,
            "--kmer-size",
            int(kmer_size) if kmer_size else 15,
            "--kmer-min-abundance",
            str(kmer_min_abundance) if kmer_min_abundance else "",
            "--singleton-kmers" if singleton_kmers else "",
            "--n-cpu",
            str(n_cpu) if n_cpu else "",
            "--compression",
            compression if compression else "",
            "--temp-dir",
            temp_dir if temp_dir else "",
        ]

        if x:
            command.append("-x")
        if v:
            command.append("-v")

        subprocess.run(command)

    def create_from_contigs(
        self,
        genomic_data,
        phenotype_description,
        phenotype_metadata,
        output=None,
        kmer_size=31,
        singleton_kmers=False,
        n_cpu=None,
        compression=None,
        temp_dir=None,
        x=True,
        v=False,
    ):
        command = [
            "wsl",
            "/home/mhdeeb/kover/bin/kover",
            "dataset",
            "create",
            "from-contigs",
            "--genomic-data",
            genomic_data,
            "--phenotype-description",
            phenotype_description,
            "--phenotype-metadata",
            phenotype_metadata,
            "--output",
            output,
            "--kmer-size",
            str(kmer_size),
        ]

        if x:
            command.append("-x")
        if v:
            command.append("-v")

        return command

    def create_from_tsv(
        self,
        genomic_data,
        phenotype_description=None,
        phenotype_metadata=None,
        output=None,
        compression=None,
        x=True,
        v=False,
    ):
        command = [
            "wsl",
            "/home/mhdeeb/kover/bin/kover",
            "dataset",
            "create",
            "from-tsv",
            "--genomic-data",
            genomic_data,
            "--phenotype-description",
            phenotype_description,
            "--phenotype-metadata",
            phenotype_metadata,
            "--output",
            output,
        ]

        if x:
            command.append("-x")
        if v:
            command.append("-v")

        return command

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
