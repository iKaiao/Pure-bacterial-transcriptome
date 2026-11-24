# Pure-bacterial-transcriptome

[![Rust 1.97.1+](https://img.shields.io/badge/Rust-1.97.1%2B-orange)](Cargo.toml)
[![MIT license](https://img.shields.io/badge/License-MIT-blue)](LICENSE)

Pure-bacterial-transcriptome is a Rust command-line workflow for reference-based differential expression analysis of paired-end bacterial RNA-seq. It takes paired FASTQ reads, a sample manifest, and a reference genome in FASTA format. You can supply a matching GFF3 annotation or generate one with Prokka.

The workflow trims reads with Trimmomatic, maps them with Bowtie2, and counts reads per gene with featureCounts. It uses rsdeseq2 for DESeq2-compatible differential expression and produces raw and normalized count tables, read and alignment quality reports, and PCA, sample correlation, MA, and volcano plots. Differential expression and plotting run directly in Rust.

The project is open source under the [MIT License](LICENSE). For questions and bug reports, use [GitHub Issues](https://github.com/iKaiao/Pure-bacterial-transcriptome/issues).

## Repository files

| File or directory | Purpose |
| --- | --- |
| [src/](src/) | Rust workflow source. The command-line entry point is [main.rs](src/main.rs). |
| [Cargo.toml](Cargo.toml) and [Cargo.lock](Cargo.lock) | Cargo package metadata and locked Rust dependencies. |
| [bacterial_rnaseq.yml](bacterial_rnaseq.yml) | Conda environment for running a Cargo build. |
| [conda/](conda/) | Recipe and compiler configuration for building a local Conda package. |
| [samples.example.tsv](samples.example.tsv) | Template for sample names, conditions, and paired FASTQ paths. |

## Building and installing

Build with [Cargo](INSTALL.md#build-from-source) or [package and install with Conda](INSTALL.md#build-and-install-the-conda-package). The [installation guide](INSTALL.md) lists the runtime tools and platform requirements.

## Usage

- [Prepare the inputs](docs/usage.md#prepare-the-inputs)
- [Run the workflow](docs/usage.md#run-the-workflow)
- [Read the output files](docs/usage.md#output-files)
- [Resume a run](docs/usage.md#resuming-an-interrupted-run)
