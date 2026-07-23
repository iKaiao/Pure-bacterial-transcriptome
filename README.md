# Pure-bacterial-transcriptome

A reusable paired-end bacterial RNA-seq differential-expression workflow:

`Trimmomatic → BWA → featureCounts → DESeq2`

The pipeline is intended for bacterial RNA-seq data with a reference genome. Supply a GFF3 annotation when available. If no GFF3 file is supplied, the script automatically annotates the exact reference FASTA used for mapping with Prokka before read counting and differential-expression analysis.

## Repository contents

| File | Purpose |
| --- | --- |
| `run_bacterial_rnaseq_deseq2.sh` | Main workflow script |
| `samples.example.tsv` | Example sample manifest |
| `bacterial_rnaseq.yml` | Conda environment with all required tools |

## 1. Create the Conda environment

From the repository directory, run:

```bash
conda config --set channel_priority strict
conda env create -f bacterial_rnaseq.yml
conda activate BactRNAseq
```

The environment contains Trimmomatic, BWA, SAMtools, Subread/featureCounts, FastQC, Prokka, R, and DESeq2. For later runs, simply activate it with `conda activate BactRNAseq`.

## 2. Prepare the inputs

### Reference genome

Use `--genome` to provide the reference genome FASTA. Common extensions such as `.fa`, `.fna`, and `.fasta` are all accepted. This must be the genome sequence to which the RNA-seq reads will be mapped.

### GFF3 annotation (optional)

If a GFF3 annotation is available, provide it with `--gff annotation.gff`. **The contig names in column 1 of the GFF must exactly match the sequence names in the reference FASTA.** Otherwise, featureCounts cannot assign reads to genes.

If no GFF is available, omit `--gff`. The script copies the reference FASTA into the output directory and annotates that copy with Prokka. The generated annotation is saved as:

```text
<output_directory>/prokka_annotation/annotation.gff
```

### Sample manifest

Copy and edit the template:

```bash
cp samples.example.tsv samples.tsv
```

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

- `sample` values must be unique and may contain only letters, numbers, `.`, `_`, and `-`.
- `condition` values must exactly match `--condition-a` and `--condition-b`.
- `r1` and `r2` may be absolute paths or paths relative to the directory from which the script is launched.
- Each group needs at least two biological replicates; three or more are recommended.

## 3. Run the workflow

### With an existing GFF3 annotation

This example compares Agar with Glucose. Positive log2 fold changes indicate higher expression in Agar.

```bash
bash ./run_bacterial_rnaseq_deseq2.sh \
  --genome ./genome.fna \
  --gff ./genome.gff \
  --samples ./samples.tsv \
  --out ./RNAseq_TS_vs_CK \
  --condition-a TS \
  --condition-b CK \
  --strand 0 \
  --adapter /home/hello2/anaconda3/envs/BactRNAseq/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa \
  --threads 20 \
  --fastqc

```

### Without a GFF3 annotation

Omit `--gff` to trigger automatic Prokka annotation:

```bash
bash run_bacterial_rnaseq_deseq2.sh \
  --genome assembly.fna \
  --samples samples.tsv \
  --out RNAseq_conditionA_vs_conditionB \
  --condition-a conditionA \
  --condition-b conditionB \
  --strand 0 \
  --threads 12
```

The script automatically looks for `TruSeq3-PE.fa` in the active Conda environment. For a different library adapter, provide it explicitly:

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
├── reference/                     # Reference FASTA used for mapping and BWA index files
├── bam/                           # Coordinate-sorted and indexed BAM files
├── qc/                            # flagstat, idxstats, and optional FastQC reports
├── counts/
│   ├── raw_counts.sN.txt          # featureCounts raw read-count table
│   └── raw_counts.sN.txt.summary  # Read-assignment summary for each sample
├── deseq2/
│   ├── DESeq2_<A>_vs_<B>.csv      # Differential-expression results
│   ├── raw_counts_matrix.csv      # Raw count matrix
│   ├── DESeq2_normalized_counts.csv
│   ├── DESeq2_size_factors.csv
│   ├── PCA.pdf
│   ├── sample_correlation.pdf
│   ├── MA_plot.pdf
│   ├── volcano_DESeq2_<A>_vs_<B>.pdf
│   └── DESeq2_summary.txt
└── logs/                          # Step-specific logs and the complete pipeline.log
```

DESeq2 uses its default **median-of-ratios** size-factor normalization, negative-binomial Wald tests, and Benjamini–Hochberg multiple-testing adjustment. Before DESeq2, a gene must have a raw count of at least 10 in at least as many samples as the smaller group size. In the results:

- `log2FoldChange > 0` means higher expression in `condition-a`;
- `padj < 0.05` denotes BH-adjusted statistical significance;
- the `direction` column additionally requires `|log2FoldChange| >= 1` to label a gene as upregulated.

## Resuming an interrupted run

Run the exact same command again with the **same `--out` directory**. The script skips completed trimming outputs, sorted/indexed BAM files, an existing count table, and an existing Prokka annotation. DESeq2 result files are regenerated.

## Troubleshooting

1. **featureCounts reports zero Assigned reads:** first check that GFF and FASTA contig names match exactly; then verify `--strand`.
2. **Trimmomatic cannot find the adapter FASTA:** pass the full path explicitly with `--adapter /full/path/TruSeq3-PE.fa`.
3. **The contrast is reversed:** confirm the order of `--condition-a` and `--condition-b`. Results are always reported as `condition-a / condition-b`, so positive values are higher in the first condition.
4. **No GFF is available:** do not manually rename contigs. Omit `--gff` and let the script annotate the current reference FASTA with Prokka.

