# Installation

[Project README](README.md) · [Usage guide](docs/usage.md)

## Build from source

From the repository directory, run:

```bash
conda config --set channel_priority strict
conda env create -f bacterial_rnaseq.yml
conda activate BactRNAseq
```

The environment contains Trimmomatic, Bowtie2, Subread/featureCounts, and a fallback font for the plots. For automatic GFF3 annotation, add Prokka to `PATH` and see [running without a GFF3 annotation](docs/usage.md#without-a-gff3-annotation).

Build the Rust executable with Cargo and Rust 1.97.1 or newer. Run it from the environment above so the analysis tools are on `PATH`.

```bash
cargo build --locked --release
./target/release/pure-bacterial-transcriptome --help
```

On Linux and macOS, use a Conda platform where the analysis dependencies resolve. On Windows, run the full workflow from Ubuntu under WSL 2 with Conda installed inside Ubuntu. Use paths visible to Ubuntu in command arguments and in the sample manifest, for example `/mnt/c/data/reads/A1.R1.fq.gz`.

## Build and install the Conda package

With `conda-build` available, build from the repository root and install the resulting package from the local channel. The recipe invokes Cargo with the versions recorded in `Cargo.lock`.

```bash
CONDA_CHANNEL_PRIORITY=strict conda build --override-channels -c conda-forge -c bioconda \
  --no-anaconda-upload --output-folder conda-bld conda
conda create --override-channels --strict-channel-priority -c "$PWD/conda-bld" \
  -c conda-forge -c bioconda -n BactRNAseqRust pure-bacterial-transcriptome
conda activate BactRNAseqRust
pure-bacterial-transcriptome --version
```

The package installs the executable and the tools needed when providing a GFF3 file. To generate an annotation automatically, make Prokka available on `PATH`. When its Perl dependencies conflict with the analysis environment, install Prokka separately and append its environment's `bin` directory to `PATH` after activating the analysis environment. Rust gives that directory priority only while running Prokka. A Rust toolchain is needed to build from source outside Conda.
