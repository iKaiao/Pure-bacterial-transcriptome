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
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cp -- "$SCRIPT_DIR/run_deseq2.R" "$RSCRIPT"

echo "[$(date '+%F %T')] Running DESeq2"
Rscript "$RSCRIPT" "$COUNTS" "$OUTDIR/samples.tsv" "$CONDITION_A" "$CONDITION_B" "$DESEQDIR" > "$LOGDIR/DESeq2.log" 2>&1

echo "[$(date '+%F %T')] Workflow complete"
echo "DESeq2 results: $DESEQDIR/DESeq2_${TAG}.csv"
echo "Full log: $LOGDIR/pipeline.log"
