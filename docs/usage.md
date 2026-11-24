# Usage guide

[Project README](../README.md) · [Installation](../INSTALL.md)

## Prepare the inputs

### Reference genome

Use `--genome` to provide the reference genome FASTA. Common extensions such as `.fa`, `.fna`, and `.fasta` are all accepted. This must be the genome sequence to which the RNA-seq reads will be mapped.

### GFF3 annotation (optional)

If a GFF3 annotation is available, provide it with `--gff annotation.gff`. **The contig names in column 1 of the GFF must exactly match the sequence names in the reference FASTA.** Otherwise, featureCounts cannot assign reads to genes.

If no GFF is available, omit `--gff` with Prokka on `PATH`. The workflow copies the reference FASTA into the output directory and annotates that copy. The generated annotation is saved as:

```text
<output_directory>/prokka_annotation/annotation.gff
```

### Adapter FASTA for Trimmomatic

The `--adapter` file is **not a genome or gene annotation**. It is a small FASTA file containing the sequencing-adapter sequences used to construct the library. In short inserts, a read can extend through the cDNA insert and into this adapter sequence. Trimmomatic removes that non-biological sequence before mapping; otherwise, read quality and alignment can be adversely affected.

Choose the adapter according to the **library-preparation kit stated by the sequencing provider or in the library construction report**. It does not depend on the bacterium, the treatment, or whether the reads are agar/glucose samples. For a normal paired-end Illumina TruSeq library, use `TruSeq3-PE.fa`. Do not choose an adapter merely because its filename looks newer.

| Library kit / evidence in the report | Adapter file to use |
| --- | --- |
| Illumina TruSeq paired-end library; no conflicting kit information | `TruSeq3-PE.fa` |
| Explicitly states an older TruSeq 2 paired-end kit | `TruSeq2-PE.fa` |
| Explicitly states Nextera paired-end library preparation | `NexteraPE-PE.fa` |
| A custom kit or the provider supplies adapter sequences | Use the supplied adapter FASTA or request the exact paired-end adapter sequences from the provider |

For paired-end reads, select a file with `PE` in its name, not an `SE` (single-end) file. `TruSeq3-PE.fa` contains the mate-specific adapters needed by Trimmomatic's `ILLUMINACLIP` step.

The executable can locate `TruSeq3-PE.fa` under the active Conda environment's `share` directory. To select a different adapter, supply its path with `--adapter`.

### Sample manifest

Make a copy of [samples.example.tsv](../samples.example.tsv) named `samples.tsv` and edit it for your experiment.

The manifest must be **tab-delimited (TSV)** and contain exactly four columns. Each row is one paired-end library:

```tsv
sample	condition	r1	r2
A1	Agar	reads/A1.R1.fq.gz	reads/A1.R2.fq.gz
A2	Agar	reads/A2.R1.fq.gz	reads/A2.R2.fq.gz
A3	Agar	reads/A3.R1.fq.gz	reads/A3.R2.fq.gz
G1	Glucose	reads/G1.R1.fq.gz	reads/G1.R2.fq.gz
G2	Glucose	reads/G2.R1.fq.gz	reads/G2.R2.fq.gz
G3	Glucose	reads/G3.R1.fq.gz	reads/G3.R2.fq.gz
```

Notes:

- `sample` values must be unique ignoring case and may contain only letters, numbers, `.`, `_`, and `-`.
- For a two-group experiment, `condition` values must match `--condition-a` or `--condition-b` exactly.
- `r1` and `r2` may be absolute paths or paths relative to the working directory.
- Each group needs at least two biological replicates; three or more are recommended.

## Run the workflow

### With an existing GFF3 annotation

This example compares Agar with Glucose. Positive log2 fold changes indicate higher expression in Agar.
Use the executable installed by Conda, or substitute `./target/release/pure-bacterial-transcriptome` for a Cargo build.

```bash
pure-bacterial-transcriptome \
  --genome ./genome.fna \
  --gff ./genome.gff \
  --samples ./samples.tsv \
  --out ./RNAseq_Agar_vs_Glucose \
  --condition-a Agar \
  --condition-b Glucose \
  --strand 0 \
  --threads 20

```

### Without a GFF3 annotation

With Prokka on `PATH`, omit `--gff` to trigger automatic annotation.

```bash
pure-bacterial-transcriptome \
  --genome assembly.fna \
  --samples samples.tsv \
  --out RNAseq_conditionA_vs_conditionB \
  --condition-a conditionA \
  --condition-b conditionB \
  --strand 0 \
  --threads 12
```

For a different library adapter, append this option:

```bash
  --adapter /path/to/adapters.fa
```

## Library strandedness

`--strand` is passed to featureCounts as `-s`:

| Value | Library type |
| --- | --- |
| `0` | Unstranded |
| `1` | Forward-stranded |
| `2` | Reverse-stranded |

Use the library-preparation kit documentation to select the correct value. If it is unknown, run featureCounts with `-s 0`, `-s 1`, and `-s 2` on the same BAM files and generally select the direction with a substantially higher number of assigned reads.

## Output files

The major outputs are written under the directory passed to `--out`:

```text
<out>/
├── clean/                         # Adapter- and quality-trimmed paired/unpaired FASTQ files
├── reference/                     # Reference FASTA used for mapping and Bowtie2 index files
├── bam/                           # Coordinate-sorted and indexed BAM files
├── qc/                            # Rust read-quality, flagstat, and idxstats reports
│   ├── read_qc_report.html
│   ├── read_qc_summary.tsv
│   ├── read_qc_by_cycle.tsv
│   └── read_qc_lengths.tsv
├── counts/
│   ├── raw_counts.sN.txt          # featureCounts raw read-count table
│   └── raw_counts.sN.txt.summary  # Read-assignment summary for each sample
├── deseq2/
│   ├── DESeq2_<A>_vs_<B>.csv      # Differential-expression results
│   ├── raw_counts_matrix.csv      # Raw count matrix
│   ├── library_sizes_raw_counts.csv
│   ├── sample_metadata.csv
│   ├── DESeq2_normalized_counts.csv
│   ├── DESeq2_size_factors.csv
│   ├── PCA.pdf
│   ├── sample_correlation.csv
│   ├── sample_correlation.pdf
│   ├── MA_plot.pdf
│   ├── volcano_DESeq2_<A>_vs_<B>.pdf
│   ├── volcano_DESeq2_<A>_vs_<B>.png
│   ├── volcano_source_data_<A>_vs_<B>.csv
│   └── DESeq2_summary.txt
└── logs/                          # Step-specific logs and the complete pipeline.log
```

The read-quality reports cover trimmed paired FASTQ files. They include read lengths, base composition, Phred quality, Q30 and per-cycle statistics. Use `--no-read-qc` to skip these reports. The older options `--fastqc` and `--no-fastqc` are accepted as aliases for the Rust reports.

The analysis uses **median-of-ratios** size-factor normalization, negative-binomial Wald tests, Cook's distance filtering, and Benjamini–Hochberg multiple-testing adjustment with independent filtering. In a two-condition manifest, a gene must have a raw count of at least 10 in at least as many samples as the smaller group size. With additional conditions, all genes enter the model, matching the previous workflow. Output filenames retain the `DESeq2` prefix for downstream scripts.

- `log2FoldChange > 0` means higher expression in `condition-a`;
- `padj < 0.05` denotes BH-adjusted statistical significance;
- the `direction` column additionally requires `|log2FoldChange| >= 1` to label a gene as upregulated.

## Resuming an interrupted run

Run the same command again with the **same `--out` directory**. The executable uses checkpoints under `.bact-rnaseq/` to reuse completed steps when their inputs and outputs have the same paths, sizes and modification times. Differential-expression result files are regenerated. Use a fresh output directory to preserve a previous analysis.

## Mapping command used by this workflow

The mapping step uses the original Bowtie2 settings, `bowtie2 --very-sensitive -p 8`. Rust streams Bowtie2's SAM output into coordinate-sorted BAM, using temporary runs and a bounded merge for large inputs. Each BAM has a BAI index and alignment QC metrics. Bowtie2's overall alignment rate is logged per sample.

## Troubleshooting

1. **featureCounts reports zero Assigned reads:** first check that GFF and FASTA contig names match exactly; then verify `--strand`.
2. **Trimmomatic cannot find the adapter FASTA:** pass the full path explicitly with `--adapter /full/path/TruSeq3-PE.fa`.
3. **The contrast is reversed:** confirm the order of `--condition-a` and `--condition-b`. Results are always reported as `condition-a / condition-b`, so positive values are higher in the first condition.
4. **No GFF is available:** omit `--gff` to annotate the current reference FASTA with Prokka.
