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

meta <- read.delim(
  sample_file, header = TRUE, check.names = FALSE, colClasses = "character",
  quote = "", na.strings = NULL, fileEncoding = "UTF-8-BOM"
)
required_meta <- c("sample", "condition", "r1", "r2")
if (!all(required_meta %in% colnames(meta))) {
  stop("Sample manifest must have these exact columns: sample, condition, r1, r2")
}
meta <- meta[, required_meta]
if (anyDuplicated(meta$sample)) stop("Sample names must be unique")
if (!all(c(condition_a, condition_b) %in% meta$condition)) stop("Both contrast conditions must occur in samples.tsv")

fc <- read.delim(
  count_file, comment.char = "#", check.names = FALSE,
  stringsAsFactors = FALSE, colClasses = c(Geneid = "character")
)
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

# vst() requires 1,000 genes with mean normalized count greater than five.
vsd <- if (sum(rowMeans(norm_matrix) > 5) >= 1000L) {
  vst(dds, blind = FALSE)
} else {
  varianceStabilizingTransformation(dds, blind = FALSE)
}
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

png(file.path(out_dir, paste0("volcano_DESeq2_", tag, ".png")),
    width = 3900, height = 2700, res = 300, bg = "white")
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
