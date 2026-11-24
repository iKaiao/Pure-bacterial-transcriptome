use crate::{Args, Sample};
use anyhow::{Context, Result, ensure};
use csv::{ReaderBuilder, Writer};
use regex::Regex;
use rsdeseq2::prelude::{
    CooksCutoff, CooksReplacementOptions, CountMatrix, DEFAULT_FAST_VST_NSUB, DeseqBuilder,
    DeseqResultRow, DesignMatrix, FactorLevelContrast, resolve_cooks_cutoff,
};
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::Path;

pub const ALPHA: f64 = 0.05;
const COUNT_COLUMNS: [&str; 6] = ["Geneid", "Chr", "Start", "End", "Strand", "Length"];
const RESULT_COLUMNS: [&str; 8] = [
    "Geneid",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "direction",
];

pub fn run(
    args: &Args,
    counts_file: &Path,
    samples: &[Sample],
    tag: &str,
    out: &Path,
    log: &mut File,
) -> Result<()> {
    let counts = read_counts(counts_file, samples)?;
    let names: Vec<_> = samples.iter().map(|sample| sample.id.as_str()).collect();
    let genes = counts.gene_names().context("Missing gene identifiers")?;
    write_matrix(
        out,
        "raw_counts_matrix.csv",
        genes,
        &names,
        counts.as_slice(),
    )?;

    let indices_a = condition_indices(samples, &args.condition_a);
    let indices_b = condition_indices(samples, &args.condition_b);
    let tested = tested_genes(&counts, &indices_a, &indices_b);
    ensure!(
        !tested.is_empty(),
        "No genes passed the minimum count filter"
    );
    let mut library_sizes = vec![0u64; samples.len()];
    for row in counts.as_slice().chunks_exact(samples.len()) {
        for (total, &value) in library_sizes.iter_mut().zip(row) {
            *total += u64::from(value);
        }
    }
    let mut libraries = Writer::from_path(out.join("library_sizes_raw_counts.csv"))?;
    libraries.write_record(["sample", "library_size"])?;
    let mut metadata = Writer::from_path(out.join("sample_metadata.csv"))?;
    metadata.write_record(["sample", "condition"])?;
    for (sample, total) in samples.iter().zip(library_sizes) {
        libraries.serialize((&sample.id, total))?;
        metadata.write_record([&sample.id, &sample.condition])?;
    }
    libraries.flush()?;
    metadata.flush()?;

    let counts = counts.select_rows(&tested)?;
    let design = factor_design(samples, &args.condition_b)?;
    let sample_levels: Vec<_> = samples
        .iter()
        .map(|sample| sample.condition.clone())
        .collect();
    let contrast = FactorLevelContrast {
        factor: "condition",
        numerator: &args.condition_a,
        denominator: &args.condition_b,
        reference: Some(&args.condition_b),
        sample_levels: &sample_levels,
    };
    let builder = DeseqBuilder::new()
        .automatic_dispersion_means()
        .threads(args.threads as usize)
        .independent_filtering_alpha(ALPHA);
    let cutoff = resolve_cooks_cutoff(
        CooksCutoff::default(),
        design.n_samples(),
        design.n_coefficients(),
    )?
    .context("Cook's distance cutoff could not be calculated")?;
    writeln!(
        log,
        "rsdeseq2 0.5.0, {} genes passed the count filter",
        tested.len()
    )?;
    let fit = builder.fit_wald_glm_mu_factor_level_contrast_with_cooks_replacement(
        &counts,
        &design,
        contrast,
        &CooksReplacementOptions::new(cutoff),
    )?;
    let original = &fit.original_fit;
    let normalized = original.normalized_counts(&counts)?;
    let genes = counts
        .gene_names()
        .context("Missing filtered gene identifiers")?;
    write_matrix(
        out,
        "DESeq2_normalized_counts.csv",
        genes,
        &names,
        normalized.as_slice(),
    )?;
    let mut factors = Writer::from_path(out.join("DESeq2_size_factors.csv"))?;
    factors.write_record(["sample", "size_factor"])?;
    for (name, factor) in names.iter().zip(&original.size_factors) {
        factors.serialize((name, factor))?;
    }
    factors.flush()?;

    let mut rows: Vec<_> = fit.results.rows.iter().enumerate().collect();
    rows.sort_by(|(_, a), (_, b)| missing_last(a.padj).total_cmp(&missing_last(b.padj)));
    let summaries = write_results(
        args,
        tag,
        out,
        &rows,
        normalized.as_slice(),
        &indices_a,
        &indices_b,
    )?;

    let eligible = normalized
        .as_slice()
        .chunks_exact(samples.len())
        .filter(|row| row.iter().sum::<f64>() / samples.len() as f64 > 5.0)
        .count();
    let vst = if eligible >= DEFAULT_FAST_VST_NSUB {
        writeln!(log, "Fast VST ({eligible} eligible genes)")?;
        builder
            .size_factors(original.size_factors.clone())
            .default_fast_vst_glm_mu(&counts, &design)?
            .transformed
    } else {
        writeln!(log, "Full-data VST ({eligible} eligible genes)")?;
        original.variance_stabilizing_transform(&counts)?
    };
    crate::deseq_plots::write(out, tag, args, samples, vst.as_slice(), &summaries)?;
    writeln!(log, "Differential expression and plots complete")?;
    Ok(())
}

fn read_counts(file: &Path, samples: &[Sample]) -> Result<CountMatrix> {
    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .comment(Some(b'#'))
        .quoting(false)
        .from_path(file)?;
    let header = reader.headers()?;
    ensure!(
        header.len() == COUNT_COLUMNS.len() + samples.len()
            && header.iter().take(COUNT_COLUMNS.len()).eq(COUNT_COLUMNS),
        "Unexpected featureCounts output layout in {}",
        file.display()
    );
    for (column, sample) in header.iter().skip(COUNT_COLUMNS.len()).zip(samples) {
        let id = column
            .rsplit(['/', '\\'])
            .next()
            .and_then(|name| name.strip_suffix(".sorted.bam"));
        ensure!(
            id == Some(&sample.id),
            "Sample order in featureCounts does not equal samples.tsv: expected {}, found {column}",
            sample.id
        );
    }
    let mut genes = Vec::new();
    let mut ids = HashSet::new();
    let mut values = Vec::new();
    for record in reader.records() {
        let record = record?;
        let id = &record[0];
        ensure!(
            !id.is_empty() && ids.insert(id.to_owned()),
            "Empty or duplicate gene ID {id:?}"
        );
        genes.push(id.to_owned());
        for value in record.iter().skip(COUNT_COLUMNS.len()) {
            values.push(
                value
                    .parse::<u32>()
                    .with_context(|| format!("Invalid count {value:?} for gene {id:?}"))?,
            );
        }
    }
    Ok(CountMatrix::from_row_major_u32_with_names(
        genes.len(),
        samples.len(),
        values,
        Some(genes),
        Some(samples.iter().map(|sample| sample.id.clone()).collect()),
    )?)
}

fn condition_indices(samples: &[Sample], condition: &str) -> Vec<usize> {
    samples
        .iter()
        .enumerate()
        .filter_map(|(index, sample)| (sample.condition == condition).then_some(index))
        .collect()
}

fn tested_genes(counts: &CountMatrix, a: &[usize], b: &[usize]) -> Vec<usize> {
    // The R condition factor retains unused levels in table() after subsetting.
    // A third condition makes the old prefilter threshold zero.
    if a.len() + b.len() < counts.n_samples() {
        return counts.gene_indices().into_iter().collect();
    }
    let minimum_replicates = a.len().min(b.len());
    counts
        .as_slice()
        .chunks_exact(counts.n_samples())
        .enumerate()
        .filter_map(|(gene, row)| {
            let expressed = a
                .iter()
                .chain(b)
                .filter(|&&column| row[column] >= 10)
                .count();
            (expressed >= minimum_replicates).then_some(gene)
        })
        .collect()
}

fn factor_design(samples: &[Sample], reference: &str) -> Result<DesignMatrix> {
    let mut levels: Vec<_> = samples
        .iter()
        .map(|sample| sample.condition.as_str())
        .filter(|&condition| condition != reference)
        .collect();
    levels.sort_unstable();
    levels.dedup();
    let mut values = Vec::new();
    for sample in samples {
        values.push(1.0);
        values.extend(
            levels
                .iter()
                .map(|&level| f64::from(u8::from(sample.condition == level))),
        );
    }
    let coefficients = std::iter::once("Intercept".to_owned())
        .chain(
            levels
                .iter()
                .map(|level| format!("condition_{level}_vs_{reference}")),
        )
        .collect();
    Ok(DesignMatrix::from_row_major(
        samples.len(),
        1 + levels.len(),
        values,
        Some(coefficients),
    )?)
}

fn write_matrix<T: serde::Serialize>(
    out: &Path,
    name: &str,
    genes: &[String],
    samples: &[&str],
    values: &[T],
) -> Result<()> {
    let mut writer = Writer::from_path(out.join(name))?;
    writer.write_record(std::iter::once("Geneid").chain(samples.iter().copied()))?;
    for (gene, row) in genes.iter().zip(values.chunks_exact(samples.len())) {
        writer.serialize((gene, row))?;
    }
    writer.flush()?;
    Ok(())
}

#[derive(Clone, Copy, PartialEq)]
pub enum Significance {
    Down,
    NonSignificant,
    SignificantDown,
    SignificantUp,
    Up,
}

impl Significance {
    pub fn label(self) -> &'static str {
        match self {
            Self::Down => "Down",
            Self::NonSignificant => "Non-significant",
            Self::SignificantDown => "Significant Down",
            Self::SignificantUp => "Significant Up",
            Self::Up => "Up",
        }
    }
}

pub struct GeneSummary {
    pub base_mean: f64,
    pub lfc: Option<f64>,
    pub padj: Option<f64>,
    pub neg_log_padj: f64,
    pub significance: Significance,
    pub label: String,
    pub show_label: bool,
}

fn write_results(
    args: &Args,
    tag: &str,
    out: &Path,
    rows: &[(usize, &DeseqResultRow)],
    normalized: &[f64],
    indices_a: &[usize],
    indices_b: &[usize],
) -> Result<Vec<GeneSummary>> {
    let mut results = Writer::from_path(out.join(format!("DESeq2_{tag}.csv")))?;
    results.write_record(RESULT_COLUMNS)?;
    let mut volcano = Writer::from_path(out.join(format!("volcano_source_data_{tag}.csv")))?;
    volcano.write_record([
        "Geneid",
        "baseMean",
        "log2FoldChange",
        "padj",
        "neg_log10_padj",
        "Differential_Significance",
        "mean_expression_a",
        "mean_expression_b",
        "both_conditions_expression",
        "display_label",
    ])?;
    let delimited_locus = Regex::new(r"[|;]\s*J[[:alnum:]_.-]+.*$")?;
    let locus_suffix = Regex::new(r"\s+J[[:alnum:]_.-]+$")?;
    let bare_locus = Regex::new(r"^J[[:alnum:]_.-]+$")?;
    let width = normalized.len() / rows.len();
    let mean = |row: &[f64], columns: &[usize]| {
        columns.iter().map(|&column| row[column]).sum::<f64>() / columns.len() as f64
    };
    let up_a = format!("Up in {}", args.condition_a);
    let up_b = format!("Up in {}", args.condition_b);
    let mut counts = (0, 0);
    let mut genes = Vec::with_capacity(rows.len());
    let mut both_means = Vec::with_capacity(rows.len());

    for &(index, row) in rows {
        let gene = row
            .gene
            .as_deref()
            .context("Missing result gene identifier")?;
        let lfc = row.log2_fold_change.filter(|value| value.is_finite());
        let padj = row.padj.filter(|value| value.is_finite());
        let significant = padj.is_some_and(|p| p < ALPHA);
        let (significance, direction) = match (lfc, significant) {
            (Some(v), true) if v >= 1.0 => {
                counts.0 += 1;
                (Significance::SignificantUp, up_a.as_str())
            }
            (Some(v), true) if v <= -1.0 => {
                counts.1 += 1;
                (Significance::SignificantDown, up_b.as_str())
            }
            (Some(v), _) if v >= 1.0 => (Significance::Up, "Not significant"),
            (Some(v), _) if v <= -1.0 => (Significance::Down, "Not significant"),
            _ => (Significance::NonSignificant, "Not significant"),
        };
        results.serialize((
            gene,
            row.base_mean,
            number(lfc),
            number(row.lfc_se),
            number(row.stat),
            number(row.pvalue),
            number(padj),
            direction,
        ))?;
        let expression = &normalized[index * width..(index + 1) * width];
        let mean_a = mean(expression, indices_a);
        let mean_b = mean(expression, indices_b);
        let both = mean_a.min(mean_b);
        both_means.push(both);
        let label = delimited_locus.replace(gene, "");
        let label = locus_suffix.replace(&label, "");
        let label = if label.is_empty() {
            gene.to_owned()
        } else {
            label.into_owned()
        };
        let neg_log_padj = padj.map_or(0.0, |p| -p.max(f64::MIN_POSITIVE).log10());
        volcano.serialize((
            gene,
            row.base_mean,
            number(lfc),
            number(padj),
            neg_log_padj,
            significance.label(),
            mean_a,
            mean_b,
            both,
            &label,
        ))?;
        genes.push(GeneSummary {
            base_mean: row.base_mean,
            lfc,
            padj,
            neg_log_padj,
            significance,
            label,
            show_label: false,
        });
    }
    results.flush()?;
    volcano.flush()?;

    let mut positive: Vec<_> = both_means
        .iter()
        .copied()
        .filter(|v| *v > 0.0 && v.is_finite())
        .collect();
    positive.sort_by(f64::total_cmp);
    let cutoff = if positive.is_empty() {
        f64::INFINITY
    } else {
        let index = (positive.len() - 1) as f64 * 0.98;
        let fraction = index.fract();
        positive[index.floor() as usize] * (1.0 - fraction)
            + positive[index.ceil() as usize] * fraction
    };
    let mut label_indices: Vec<_> = (0..genes.len())
        .filter(|&index| {
            let gene = &genes[index];
            gene.lfc.is_some_and(|v| v.abs() >= 1.0)
                && gene
                    .padj
                    .is_some_and(|p| p < 1e-4 || p < ALPHA && both_means[index] >= cutoff)
        })
        .collect();
    label_indices.sort_by(|&a, &b| {
        missing_last(genes[a].padj)
            .total_cmp(&missing_last(genes[b].padj))
            .then_with(|| {
                genes[b]
                    .lfc
                    .unwrap_or_default()
                    .abs()
                    .total_cmp(&genes[a].lfc.unwrap_or_default().abs())
            })
    });
    for index in label_indices.into_iter().take(25) {
        genes[index].show_label = !bare_locus.is_match(&genes[index].label);
    }

    fs::write(
        out.join("DESeq2_summary.txt"),
        format!(
            "Contrast: {} / {} (positive log2FC = higher in {})\n\
        Normalization: DESeq2 median-of-ratios size factors\n\
        Test: negative-binomial Wald test; Benjamini-Hochberg adjusted P values\n\
        Genes tested after count filter: {}\n\
        {up_a} (padj < 0.05, log2FC >= 1): {}\n\
        {up_b} (padj < 0.05, log2FC <= -1): {}\n",
            args.condition_a,
            args.condition_b,
            args.condition_a,
            rows.len(),
            counts.0,
            counts.1,
        ),
    )?;
    Ok(genes)
}

fn number(value: Option<f64>) -> String {
    value
        .filter(|value| value.is_finite())
        .map_or_else(|| "NA".to_owned(), |value| format!("{value:e}"))
}

fn missing_last(value: Option<f64>) -> f64 {
    value
        .filter(|value| value.is_finite())
        .unwrap_or(f64::INFINITY)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_the_legacy_count_filter_with_and_without_a_third_group() -> Result<()> {
        let two = CountMatrix::from_row_major_u32(2, 4, vec![10, 0, 10, 0, 0, 9, 0, 0])?;
        assert_eq!(tested_genes(&two, &[0, 1], &[2, 3]), vec![0]);
        let three =
            CountMatrix::from_row_major_u32(2, 6, vec![10, 0, 10, 0, 0, 0, 0, 9, 0, 0, 100, 100])?;
        assert_eq!(tested_genes(&three, &[0, 1], &[2, 3]), vec![0, 1]);
        Ok(())
    }
}
