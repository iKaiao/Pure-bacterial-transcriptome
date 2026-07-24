#!/usr/bin/env bash
# Generic paired-end bacterial RNA-seq workflow:
# Trimmomatic -> Bowtie2 -> featureCounts -> DESeq2.
#
# Requirements in the active environment:
#   trimmomatic, bowtie2, samtools, featureCounts, Rscript, and R package DESeq2.
#
# The script is restart-safe: run the same command again with the same --out
# directory to skip completed trimming, mapping and counting steps.

set -euo pipefail
IFS=$'\n\t'

usage() {
    cat <<'USAGE'
Usage:
  bash run_bacterial_rnaseq_deseq2.sh \
    --genome reference.fasta \
    --samples samples.tsv \
    --out analysis_out \
    --condition-a Agar \
    --condition-b Glucose \
    --strand 2 \
    --adapter /path/to/TruSeq3-PE.fa

Required options:
  --genome FILE          Reference genome FASTA.
  --samples FILE         Tab-delimited sample table: sample, condition, r1, r2.
  --out DIR              Output directory (re-use this exact directory to resume).
  --condition-a TEXT     Numerator in DESeq2 contrast A / B; positive log2FC is higher in A.
  --condition-b TEXT     Denominator in DESeq2 contrast A / B.
  --strand 0|1|2         featureCounts strandedness: 0 unstranded, 1 forward, 2 reverse.

Useful optional options:
  --gff FILE             Gene annotation GFF3. If omitted, annotate the reference with Prokka.
  --adapter FILE         Trimmomatic PE adapter FASTA. Default: TruSeq3-PE.fa in the active Conda environment.
  --threads N            Bowtie2 mapping threads (default 8); Trimmomatic and featureCounts use 8, samtools sort uses 4.
  --feature TEXT         GFF feature to count (default: CDS).
  --attribute TEXT       GFF attribute used as gene ID (default: locus_tag).
  --fastqc               Run FastQC and MultiQC on trimmed paired reads (default).
  --no-fastqc            Skip FastQC and MultiQC.
  --help                 Show this help.

samples.tsv example (header required; paths are relative to the directory in which you run):
sample  condition  r1                  r2
A1      Agar       reads/A1.R1.fq.gz   reads/A1.R2.fq.gz
A2      Agar       reads/A2.R1.fq.gz   reads/A2.R2.fq.gz
A3      Agar       reads/A3.R1.fq.gz   reads/A3.R2.fq.gz
G1      Glucose    reads/G1.R1.fq.gz   reads/G1.R2.fq.gz
G2      Glucose    reads/G2.R1.fq.gz   reads/G2.R2.fq.gz
G3      Glucose    reads/G3.R1.fq.gz   reads/G3.R2.fq.gz
USAGE
}

GENOME=""
GFF=""
SAMPLES_FILE=""
OUTDIR=""
CONDITION_A=""
CONDITION_B=""
STRAND=""
ADAPTER=""
THREADS_MAP=8
THREADS_TRIM=8
THREADS_COUNT=8
THREADS_PROKKA=8
THREADS_SORT=4
FEATURE="CDS"
ATTRIBUTE="locus_tag"
RUN_FASTQC=1
JAVA_HEAP="8g"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genome) GENOME="$2"; shift 2 ;;
        --gff) GFF="$2"; shift 2 ;;
        --samples) SAMPLES_FILE="$2"; shift 2 ;;
        --out) OUTDIR="$2"; shift 2 ;;
        --condition-a) CONDITION_A="$2"; shift 2 ;;
        --condition-b) CONDITION_B="$2"; shift 2 ;;
        --strand) STRAND="$2"; shift 2 ;;
        --adapter) ADAPTER="$2"; shift 2 ;;
        --threads)
            THREADS_MAP="$2"
            shift 2
            ;;
        --feature) FEATURE="$2"; shift 2 ;;
        --attribute) ATTRIBUTE="$2"; shift 2 ;;
        --fastqc) RUN_FASTQC=1; shift ;;
        --no-fastqc) RUN_FASTQC=0; shift ;;
        --help|-h) usage; exit 0 ;;
        *) echo "ERROR: Unknown option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

for required in GENOME SAMPLES_FILE OUTDIR CONDITION_A CONDITION_B STRAND; do
    [[ -n "${!required}" ]] || { echo "ERROR: --${required,,} is required." >&2; usage >&2; exit 1; }
done
if [[ -z "$ADAPTER" && -n "${CONDA_PREFIX:-}" && -d "$CONDA_PREFIX/share" ]]; then
    ADAPTER=$(find "$CONDA_PREFIX/share" -type f -path '*/adapters/TruSeq3-PE.fa' -print -quit 2>/dev/null || true)
fi
[[ "$STRAND" =~ ^[012]$ ]] || { echo "ERROR: --strand must be 0, 1 or 2." >&2; exit 1; }
[[ -s "$GENOME" ]] || { echo "ERROR: Genome FASTA not found: $GENOME" >&2; exit 1; }
[[ -z "$GFF" || -s "$GFF" ]] || { echo "ERROR: GFF not found: $GFF" >&2; exit 1; }
[[ -s "$SAMPLES_FILE" ]] || { echo "ERROR: Sample table not found: $SAMPLES_FILE" >&2; exit 1; }
[[ -s "$ADAPTER" ]] || { echo "ERROR: Adapter FASTA not found: $ADAPTER" >&2; exit 1; }

for tool in trimmomatic bowtie2 bowtie2-build bowtie2-inspect samtools featureCounts Rscript; do
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: '$tool' is unavailable. Activate an environment containing all pipeline tools." >&2
        exit 1
    }
done
if [[ -z "$GFF" ]]; then
    command -v prokka >/dev/null 2>&1 || {
        echo "ERROR: No --gff was supplied, but prokka is unavailable. Install Prokka or provide --gff." >&2
        exit 1
    }
fi
if [[ "$RUN_FASTQC" -eq 1 ]]; then
    command -v fastqc >/dev/null 2>&1 || { echo "ERROR: fastqc is unavailable. Use --no-fastqc to skip this step." >&2; exit 1; }
    command -v multiqc >/dev/null 2>&1 || { echo "ERROR: multiqc is unavailable. Use --no-fastqc to skip this step." >&2; exit 1; }
fi

########################################
# Read and validate the sample manifest
########################################
SAMPLE_IDS=()
CONDITIONS=()
R1_FILES=()
R2_FILES=()
declare -A SEEN_SAMPLES=()

while IFS=$'\t' read -r sample condition r1 r2 extra; do
    [[ -z "${sample:-}" || "$sample" =~ ^# ]] && continue
    if [[ "${sample,,}" == "sample" && "${condition,,}" == "condition" ]]; then
        continue
    fi
    [[ -z "${condition:-}" || -z "${r1:-}" || -z "${r2:-}" || -n "${extra:-}" ]] && {
        echo "ERROR: Each sample-table line must have exactly four tab-delimited columns: sample, condition, r1, r2." >&2
        exit 1
    }
    [[ "$sample" =~ ^[A-Za-z0-9._-]+$ ]] || { echo "ERROR: Invalid sample name '$sample'. Use only letters, digits, dot, underscore or hyphen." >&2; exit 1; }
    [[ -z "${SEEN_SAMPLES[$sample]+x}" ]] || { echo "ERROR: Duplicate sample name '$sample'." >&2; exit 1; }
    [[ -s "$r1" && -s "$r2" ]] || { echo "ERROR: FASTQ not found or empty for $sample: $r1 ; $r2" >&2; exit 1; }
    SEEN_SAMPLES[$sample]=1
    SAMPLE_IDS+=("$sample")
    CONDITIONS+=("$condition")
    R1_FILES+=("$r1")
    R2_FILES+=("$r2")
done < "$SAMPLES_FILE"

[[ "${#SAMPLE_IDS[@]}" -ge 4 ]] || { echo "ERROR: At least four samples are required." >&2; exit 1; }
N_A=0
N_B=0
for condition in "${CONDITIONS[@]}"; do
    [[ "$condition" == "$CONDITION_A" ]] && ((N_A += 1))
    [[ "$condition" == "$CONDITION_B" ]] && ((N_B += 1))
done
[[ "$N_A" -ge 2 && "$N_B" -ge 2 ]] || {
    echo "ERROR: DESeq2 needs at least two biological replicates in each contrast group; found $N_A for '$CONDITION_A' and $N_B for '$CONDITION_B'." >&2
    exit 1
}

########################################
# Directories, log and input provenance
########################################
TRIMDIR="$OUTDIR/clean"
REFDIR="$OUTDIR/reference"
BAMDIR="$OUTDIR/bam"
QCDIR="$OUTDIR/qc"
COUNTDIR="$OUTDIR/counts"
DESEQDIR="$OUTDIR/deseq2"
LOGDIR="$OUTDIR/logs"
mkdir -p "$TRIMDIR" "$REFDIR" "$BAMDIR" "$QCDIR" "$COUNTDIR" "$DESEQDIR" "$LOGDIR"
exec > >(tee -a "$LOGDIR/pipeline.log") 2>&1

echo "[$(date '+%F %T')] Start bacterial RNA-seq workflow"
echo "Reference: $GENOME"
if [[ -n "$GFF" ]]; then
    echo "Annotation: $GFF"
else
    echo "Annotation: no GFF supplied; Prokka will be run automatically"
fi
echo "Contrast: $CONDITION_A / $CONDITION_B (positive log2FC = higher in $CONDITION_A)"
echo "featureCounts strandedness: -s $STRAND"
cp "$SAMPLES_FILE" "$OUTDIR/samples.tsv"

########################################
# 1. Prepare reference
########################################
REF="$REFDIR/reference.fasta"
if [[ ! -s "$REF" ]]; then
    cp "$GENOME" "$REF"
fi
BOWTIE_INDEX="$REFDIR/reference"
if ! bowtie2-inspect -s "$BOWTIE_INDEX" >/dev/null 2>&1; then
    echo "[$(date '+%F %T')] Building Bowtie2 index"
    bowtie2-build "$REF" "$BOWTIE_INDEX" > "$LOGDIR/bowtie2-build.log" 2>&1
fi

# If no annotation was supplied, create a Prokka annotation from the exact
# reference FASTA used for mapping. The resulting GFF is then used directly
# by featureCounts.
if [[ -z "$GFF" ]]; then
    PROKKA_DIR="$OUTDIR/prokka_annotation"
    PROKKA_GFF="$PROKKA_DIR/annotation.gff"
    if [[ ! -s "$PROKKA_GFF" ]]; then
        echo "[$(date '+%F %T')] No GFF supplied; annotating reference with Prokka"
        prokka --outdir "$PROKKA_DIR" --prefix annotation --force --cpus "$THREADS_PROKKA" "$REF" \
            > "$LOGDIR/prokka.log" 2>&1
    else
        echo "[$(date '+%F %T')] Existing Prokka annotation found; skipping annotation"
    fi
    GFF="$PROKKA_GFF"
fi
[[ -s "$GFF" ]] || { echo "ERROR: Annotation GFF is missing after setup: $GFF" >&2; exit 1; }

########################################
# 2. Trim adapters and low-quality sequence
########################################
export _JAVA_OPTIONS="-Xmx${JAVA_HEAP}"
CLEAN_R1=()
CLEAN_R2=()
for i in "${!SAMPLE_IDS[@]}"; do
    sample="${SAMPLE_IDS[$i]}"
    r1="${R1_FILES[$i]}"
    r2="${R2_FILES[$i]}"
    p1="$TRIMDIR/${sample}.R1.paired.fq.gz"
    u1="$TRIMDIR/${sample}.R1.unpaired.fq.gz"
    p2="$TRIMDIR/${sample}.R2.paired.fq.gz"
    u2="$TRIMDIR/${sample}.R2.unpaired.fq.gz"
    CLEAN_R1+=("$p1")
    CLEAN_R2+=("$p2")

    if [[ -s "$p1" && -s "$p2" ]] && gzip -t "$p1" "$p2" 2>/dev/null; then
        echo "[$(date '+%F %T')] Trimming already complete for $sample; skipping"
        continue
    fi
    if [[ -e "$p1" || -e "$p2" ]]; then
        echo "[$(date '+%F %T')] Incomplete clean FASTQ detected for $sample; trimming again"
    fi
    echo "[$(date '+%F %T')] Trimming $sample"

    # Write to unique temporary files first. Final clean files are only
    # replaced after the paired outputs pass a gzip integrity check, so an
    # interrupted Trimmomatic run cannot be mistaken for a completed one.
    tmp_tag=".${sample}.trim.$$"
    tp1="$TRIMDIR/${sample}.R1.paired${tmp_tag}.fq.gz"
    tu1="$TRIMDIR/${sample}.R1.unpaired${tmp_tag}.fq.gz"
    tp2="$TRIMDIR/${sample}.R2.paired${tmp_tag}.fq.gz"
    tu2="$TRIMDIR/${sample}.R2.unpaired${tmp_tag}.fq.gz"
    trimmomatic PE -threads "$THREADS_TRIM" -phred33 \
        "$r1" "$r2" "$tp1" "$tu1" "$tp2" "$tu2" \
        "ILLUMINACLIP:${ADAPTER}:2:30:10" \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 \
        > "$LOGDIR/${sample}.trimmomatic.log" 2>&1
    gzip -t "$tp1" "$tp2" 2>/dev/null || {
        echo "ERROR: Trimmomatic created an incomplete paired FASTQ for $sample. See $LOGDIR/${sample}.trimmomatic.log" >&2
        exit 1
    }
    mv -f "$tp1" "$p1"
    mv -f "$tu1" "$u1"
    mv -f "$tp2" "$p2"
    mv -f "$tu2" "$u2"
done
unset _JAVA_OPTIONS

########################################
# 3. FastQC and MultiQC on cleaned paired reads
########################################
if [[ "$RUN_FASTQC" -eq 1 && ! -s "$QCDIR/multiqc_report.html" ]]; then
    echo "[$(date '+%F %T')] FastQC on trimmed paired reads"
    fastqc -t "$THREADS_TRIM" -o "$QCDIR" "${CLEAN_R1[@]}" "${CLEAN_R2[@]}"
    echo "[$(date '+%F %T')] MultiQC summary"
    multiqc "$QCDIR" -o "$QCDIR" > "$LOGDIR/multiqc.log" 2>&1
fi

########################################
# 4. Bowtie2 mapping, sorted BAM and mapping QC
########################################

# The command record uses modern `samtools sort -o`, but some existing
# BactRNAseq installations use the older `<in.bam> <out.prefix>` interface.
# Detect that interface from its own usage string and retain the same Bowtie2
# alignment settings in both cases.
SAMTOOLS_SORT_USAGE="$(samtools sort 2>&1 || true)"
if grep -q '<in.bam> <out.prefix>' <<< "$SAMTOOLS_SORT_USAGE"; then
    SAMTOOLS_LEGACY=1
else
    SAMTOOLS_LEGACY=0
fi

map_and_sort() {
    local sample="$1"
    local r1="$2"
    local r2="$3"
    local bam_out="$4"

    if [[ "$SAMTOOLS_LEGACY" -eq 1 ]]; then
        # Old samtools cannot consume Bowtie2 SAM directly in `sort` and has
        # no `-o` option. Use completed temporary files to avoid stream EOFs.
        local tmp_tag=".${sample}.map.$$"
        local sam_tmp="$BAMDIR/${sample}${tmp_tag}.sam"
        local unsorted_bam="$BAMDIR/${sample}${tmp_tag}.unsorted.bam"

        bowtie2 --very-sensitive -p "$THREADS_MAP" \
            -x "$BOWTIE_INDEX" -1 "$r1" -2 "$r2" \
            > "$sam_tmp" 2> "$LOGDIR/${sample}.bowtie2.log"
        samtools view -bS "$sam_tmp" > "$unsorted_bam"
        samtools sort -@ "$THREADS_SORT" "$unsorted_bam" "${bam_out%.bam}"
        rm -f "$sam_tmp" "$unsorted_bam"
    else
        bowtie2 --very-sensitive -p "$THREADS_MAP" \
            -x "$BOWTIE_INDEX" -1 "$r1" -2 "$r2" \
            2> "$LOGDIR/${sample}.bowtie2.log" \
            | samtools sort -@ "$THREADS_SORT" -o "$bam_out"
    fi
}

check_bam() {
    local bam_in="$1"

    if [[ "$SAMTOOLS_LEGACY" -eq 1 ]]; then
        samtools view -H "$bam_in" >/dev/null
    else
        samtools quickcheck -q "$bam_in"
    fi
}

BAMS=()
for i in "${!SAMPLE_IDS[@]}"; do
    sample="${SAMPLE_IDS[$i]}"
    bam="$BAMDIR/${sample}.sorted.bam"
    BAMS+=("$bam")
    if [[ -s "$bam" && -s "$bam.bai" ]]; then
        echo "[$(date '+%F %T')] Mapping already complete for $sample; skipping"
    else
        echo "[$(date '+%F %T')] Mapping $sample"
        map_and_sort "$sample" "${CLEAN_R1[$i]}" "${CLEAN_R2[$i]}" "$bam"
        samtools index "$bam"
    fi
    check_bam "$bam"
    samtools flagstat "$bam" > "$QCDIR/${sample}.flagstat.txt"
    samtools idxstats "$bam" > "$QCDIR/${sample}.idxstats.txt"
done

grep "overall alignment rate" "$LOGDIR"/*.bowtie2.log || true

########################################
# 5. featureCounts
########################################
COUNTS="$COUNTDIR/raw_counts.s${STRAND}.txt"
if [[ ! -s "$COUNTS" ]]; then
    echo "[$(date '+%F %T')] featureCounts: feature=$FEATURE, attribute=$ATTRIBUTE, strand=$STRAND"
    featureCounts -T "$THREADS_COUNT" \
        -p --countReadPairs -B -C \
        -s "$STRAND" \
        -a "$GFF" -t "$FEATURE" -g "$ATTRIBUTE" \
        -o "$COUNTS" "${BAMS[@]}" \
        > "$LOGDIR/featureCounts.log" 2>&1
else
    echo "[$(date '+%F %T')] Count table already exists; skipping featureCounts"
fi
[[ -s "$COUNTS.summary" ]] || { echo "ERROR: featureCounts summary is missing." >&2; exit 1; }

########################################
# 6. DESeq2 with median-of-ratios normalization and Wald test
########################################
TAG="${CONDITION_A}_vs_${CONDITION_B}"
TAG="${TAG// /_}"
RSCRIPT="$DESEQDIR/run_deseq2_${TAG}.R"
cat > "$RSCRIPT" <<'RSCRIPT_EOF'
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript run_deseq2.R <featureCounts.txt> <samples.tsv> <condition_a> <condition_b> <outdir>")
}
count_file <- args[[1]]
sample_file <- args[[2]]
condition_a <- args[[3]]
condition_b <- args[[4]]
out_dir <- args[[5]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

meta <- read.delim(sample_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
required_meta <- c("sample", "condition", "r1", "r2")
if (!all(required_meta %in% colnames(meta))) {
  stop("Sample manifest must have these exact columns: sample, condition, r1, r2")
}
meta <- meta[, required_meta]
if (anyDuplicated(meta$sample)) stop("Sample names must be unique")
if (!all(c(condition_a, condition_b) %in% meta$condition)) stop("Both contrast conditions must occur in samples.tsv")

fc <- read.delim(count_file, comment.char = "#", check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("Geneid", "Chr", "Start", "End", "Strand", "Length") %in% colnames(fc))) {
  stop("Unexpected featureCounts output layout")
}
sample_columns <- colnames(fc)[7:ncol(fc)]
sample_ids_from_counts <- sub("\\.sorted\\.bam$", "", basename(sample_columns))
if (!identical(sample_ids_from_counts, meta$sample)) {
  stop(paste("Sample order in featureCounts does not equal samples.tsv. Counts:",
             paste(sample_ids_from_counts, collapse = ","), "; manifest:", paste(meta$sample, collapse = ",")))
}

count_matrix <- as.matrix(fc[, sample_columns, drop = FALSE])
storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- fc$Geneid
colnames(count_matrix) <- meta$sample

meta$condition <- factor(meta$condition)
meta$condition <- relevel(meta$condition, ref = condition_b)
rownames(meta) <- meta$sample
write.csv(data.frame(Geneid = rownames(count_matrix), count_matrix, check.names = FALSE),
          file.path(out_dir, "raw_counts_matrix.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(count_matrix), library_size = colSums(count_matrix)),
          file.path(out_dir, "library_sizes_raw_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(meta), condition = meta$condition),
          file.path(out_dir, "sample_metadata.csv"), row.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = meta, design = ~ condition)
target_samples <- meta$condition %in% c(condition_a, condition_b)
n_replicates <- min(table(meta$condition[target_samples]))
dds <- dds[rowSums(counts(dds)[, target_samples, drop = FALSE] >= 10) >= n_replicates, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", condition_a, condition_b), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$Geneid <- rownames(res_df)
res_df$direction <- "Not significant"
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange >= 1] <- paste("Up in", condition_a)
res_df$direction[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange <= -1] <- paste("Up in", condition_b)
res_df <- res_df[, c("Geneid", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "direction")]
res_df <- res_df[order(res_df$padj, na.last = TRUE), ]

tag <- gsub(" ", "_", paste0(condition_a, "_vs_", condition_b))
write.csv(res_df, file.path(out_dir, paste0("DESeq2_", tag, ".csv")), row.names = FALSE)
norm_matrix <- counts(dds, normalized = TRUE)
norm_counts <- as.data.frame(norm_matrix)
norm_counts$Geneid <- rownames(norm_counts)
norm_counts <- norm_counts[, c("Geneid", meta$sample)]
write.csv(norm_counts, file.path(out_dir, "DESeq2_normalized_counts.csv"), row.names = FALSE)
write.csv(data.frame(sample = colnames(dds), size_factor = sizeFactors(dds)),
          file.path(out_dir, "DESeq2_size_factors.csv"), row.names = FALSE)

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
pca <- prcomp(t(mat))
var_pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
condition_levels <- levels(meta$condition)
condition_colours <- setNames(grDevices::hcl.colors(length(condition_levels), "Dark 3"), condition_levels)
point_colours <- condition_colours[as.character(meta$condition)]

pdf(file.path(out_dir, "PCA.pdf"), width = 6, height = 5)
plot(pca$x[, 1], pca$x[, 2], pch = 19, cex = 1.5, col = point_colours,
     xlab = sprintf("PC1 (%.1f%%)", var_pct[1]), ylab = sprintf("PC2 (%.1f%%)", var_pct[2]), main = "VST PCA")
text(pca$x[, 1], pca$x[, 2], labels = colnames(mat), pos = 3, cex = 0.85)
legend("topright", legend = condition_levels, col = condition_colours, pch = 19, bty = "n")
dev.off()

corr <- cor(mat, method = "pearson")
write.csv(corr, file.path(out_dir, "sample_correlation.csv"), row.names = TRUE)
pdf(file.path(out_dir, "sample_correlation.pdf"), width = 6, height = 5)
heatmap(corr, symm = TRUE, scale = "none", col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
        margins = c(8, 8), main = "Sample correlation (VST)")
dev.off()

pdf(file.path(out_dir, "MA_plot.pdf"), width = 6, height = 5)
plotMA(res, ylim = c(-8, 8), alpha = 0.05)
dev.off()

# Volcano plot: five categories, palette, labels and layout matched to Fig. 4.
plot_df <- res_df
plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
plot_df$neg_log10_padj[is.na(plot_df$neg_log10_padj) | !is.finite(plot_df$neg_log10_padj)] <- 0

# "Down" and "Up" denote log2 fold-change direction; the significant variants
# additionally meet padj < 0.05. This preserves the five-category Fig. 4 legend.
plot_df$Differential_Significance <- "Non-significant"
plot_df$Differential_Significance[!is.na(plot_df$log2FoldChange) & plot_df$log2FoldChange <= -1] <- "Down"
plot_df$Differential_Significance[!is.na(plot_df$log2FoldChange) & plot_df$log2FoldChange >= 1] <- "Up"
plot_df$Differential_Significance[!is.na(plot_df$padj) & plot_df$padj < 0.05 & plot_df$log2FoldChange <= -1] <- "Significant Down"
plot_df$Differential_Significance[!is.na(plot_df$padj) & plot_df$padj < 0.05 & plot_df$log2FoldChange >= 1] <- "Significant Up"
plot_df$Differential_Significance <- factor(
  plot_df$Differential_Significance,
  levels = c("Down", "Non-significant", "Significant Down", "Significant Up", "Up")
)

# Keep labels selective: very significant genes, or significant genes whose
# normalized expression is high in both conditions. A trailing J-style locus
# identifier after a pipe/semicolon is removed from the displayed label.
mean_a <- rowMeans(norm_matrix[, as.character(meta$condition) == condition_a, drop = FALSE])
mean_b <- rowMeans(norm_matrix[, as.character(meta$condition) == condition_b, drop = FALSE])
plot_df$mean_expression_a <- mean_a[match(plot_df$Geneid, names(mean_a))]
plot_df$mean_expression_b <- mean_b[match(plot_df$Geneid, names(mean_b))]
plot_df$both_conditions_expression <- pmin(plot_df$mean_expression_a, plot_df$mean_expression_b)
positive_both <- plot_df$both_conditions_expression[is.finite(plot_df$both_conditions_expression) & plot_df$both_conditions_expression > 0]
both_high_cutoff <- if (length(positive_both)) unname(stats::quantile(positive_both, 0.98)) else Inf
very_significant <- !is.na(plot_df$padj) & plot_df$padj < 1e-4 & abs(plot_df$log2FoldChange) >= 1
both_high <- !is.na(plot_df$padj) & plot_df$padj < 0.05 & abs(plot_df$log2FoldChange) >= 1 &
  plot_df$both_conditions_expression >= both_high_cutoff
label_candidates <- which(very_significant | both_high)
label_candidates <- label_candidates[order(plot_df$padj[label_candidates], -abs(plot_df$log2FoldChange[label_candidates]), na.last = TRUE)]
label_candidates <- head(label_candidates, 25L)
plot_df$display_label <- sub("[|;][[:space:]]*J[[:alnum:]_.-]+.*$", "", plot_df$Geneid)
plot_df$display_label <- sub("[[:space:]]+J[[:alnum:]_.-]+$", "", plot_df$display_label)
plot_df$display_label[!nzchar(plot_df$display_label)] <- plot_df$Geneid[!nzchar(plot_df$display_label)]
# A bare locus tag (for example JIJAMHPJ_01268) has no biological meaning in a
# manuscript figure. It remains in the source-data table but is not labelled.
label_candidates <- label_candidates[!grepl("^J[[:alnum:]_.-]+$", plot_df$display_label[label_candidates])]
label_df <- plot_df[label_candidates, , drop = FALSE]

volcano_colours <- c(
  "Down" = "#A6CEE3",
  "Non-significant" = "#999999",
  "Significant Down" = "#377EB8",
  "Significant Up" = "#E41A1C",
  "Up" = "#FB9A99"
)

volcano_plot <- ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj,
                                    colour = Differential_Significance)) +
  geom_point(size = 1.9, alpha = 0.90) +
  scale_colour_manual(values = volcano_colours, drop = FALSE, name = "Differential Significance") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)), breaks = scales::pretty_breaks(n = 4)) +
  labs(
    title = sprintf("Volcano Plot of Gene Expression Differences (%s vs %s)", condition_a, condition_b),
    x = expression(log[2]*"(Fold Change)"),
    y = expression(-log[10]*"(Adjusted P Value)")
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 16)),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 11, colour = "black"),
    panel.grid.major = element_line(colour = "#EBEBEB", linewidth = 0.45),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 11),
    legend.key = element_blank(),
    legend.position = "right"
  )

if (nrow(label_df)) {
  volcano_plot <- volcano_plot +
    geom_text_repel(
      data = label_df,
      aes(x = log2FoldChange, y = neg_log10_padj, label = display_label),
      inherit.aes = FALSE,
      seed = 1,
      size = 3.2,
      colour = "#E41A1C",
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.20,
      min.segment.length = 0,
      segment.colour = "#E41A1C",
      show.legend = FALSE
    )
}

write.csv(
  plot_df[, c("Geneid", "baseMean", "log2FoldChange", "padj", "neg_log10_padj",
              "Differential_Significance", "mean_expression_a", "mean_expression_b",
              "both_conditions_expression", "display_label")],
  file.path(out_dir, paste0("volcano_source_data_", tag, ".csv")),
  row.names = FALSE
)

pdf(file.path(out_dir, paste0("volcano_DESeq2_", tag, ".pdf")), width = 13, height = 9, useDingbats = FALSE, bg = "white")
print(volcano_plot)
dev.off()

summary_lines <- c(
  paste0("Contrast: ", condition_a, " / ", condition_b, " (positive log2FC = higher in ", condition_a, ")"),
  "Normalization: DESeq2 median-of-ratios size factors",
  "Test: negative-binomial Wald test; Benjamini-Hochberg adjusted P values",
  paste("Genes tested after count filter:", nrow(res_df)),
  paste("Up in", condition_a, "(padj < 0.05, log2FC >= 1):", sum(res_df$direction == paste("Up in", condition_a))),
  paste("Up in", condition_b, "(padj < 0.05, log2FC <= -1):", sum(res_df$direction == paste("Up in", condition_b)))
)
writeLines(summary_lines, file.path(out_dir, "DESeq2_summary.txt"))
RSCRIPT_EOF

echo "[$(date '+%F %T')] Running DESeq2"
Rscript "$RSCRIPT" "$COUNTS" "$OUTDIR/samples.tsv" "$CONDITION_A" "$CONDITION_B" "$DESEQDIR" > "$LOGDIR/DESeq2.log" 2>&1

echo "[$(date '+%F %T')] Workflow complete"
echo "DESeq2 results: $DESEQDIR/DESeq2_${TAG}.csv"
echo "Full log: $LOGDIR/pipeline.log"
