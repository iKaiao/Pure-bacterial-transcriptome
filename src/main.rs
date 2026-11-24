mod bam;
mod process;
mod read_qc;
mod state;
mod workflow;

use anyhow::{Context, Result, ensure};
use clap::{ArgAction, Parser};
use std::collections::HashSet;
use std::fs;
use std::path::{Path, PathBuf, absolute};
use std::process::ExitCode;

#[derive(Parser)]
#[command(version, about, disable_help_subcommand = true)]
pub struct Args {
    /// Reference genome FASTA
    #[arg(long)]
    genome: PathBuf,
    /// Reference GFF3 (omit to annotate with Prokka)
    #[arg(long)]
    gff: Option<PathBuf>,
    /// TSV manifest with sample, condition, r1, and r2 columns
    #[arg(long)]
    samples: PathBuf,
    /// Output directory (reuse it to resume)
    #[arg(long)]
    out: PathBuf,
    /// Numerator in the DESeq2 contrast
    #[arg(long)]
    condition_a: String,
    /// Denominator in the DESeq2 contrast
    #[arg(long)]
    condition_b: String,
    /// featureCounts strandedness: 0 unstranded, 1 forward, 2 reverse
    #[arg(long, value_parser = clap::value_parser!(u8).range(0..=2))]
    strand: u8,
    /// Trimmomatic adapter FASTA (defaults to TruSeq3-PE.fa in Conda)
    #[arg(long)]
    adapter: Option<PathBuf>,
    /// Bowtie2 alignment threads
    #[arg(long, default_value_t = 8, value_parser = clap::value_parser!(u32).range(1..))]
    threads: u32,
    /// GFF feature to count
    #[arg(long, default_value = "CDS")]
    feature: String,
    /// GFF attribute used as the gene identifier
    #[arg(long, default_value = "locus_tag")]
    attribute: String,
    /// Write Rust FASTQ quality reports (enabled by default)
    #[arg(long, alias = "fastqc", action = ArgAction::SetTrue, overrides_with = "no_read_qc")]
    read_qc: bool,
    /// Skip FASTQ quality reports
    #[arg(long, alias = "no-fastqc", action = ArgAction::SetTrue, overrides_with = "read_qc")]
    no_read_qc: bool,
}

pub struct Sample {
    id: String,
    condition: String,
    r1: PathBuf,
    r2: PathBuf,
}

fn input_file(path: &Path, label: &str) -> Result<PathBuf> {
    let path = absolute(path)?;
    let metadata = fs::metadata(&path).with_context(|| format!("{label}: {}", path.display()))?;
    ensure!(
        metadata.is_file() && metadata.len() > 0,
        "{label} is not a nonempty file: {}",
        path.display()
    );
    Ok(path)
}

fn read_samples(args: &Args, manifest: &Path) -> Result<(Vec<Sample>, String)> {
    let tag = format!("{}_vs_{}", args.condition_a, args.condition_b).replace(' ', "_");
    ensure!(
        args.condition_a != args.condition_b,
        "The two contrast conditions must be different"
    );
    for condition in [&args.condition_a, &args.condition_b] {
        ensure!(
            !condition.trim().is_empty(),
            "Contrast conditions cannot be blank"
        );
        ensure!(
            !condition
                .chars()
                .any(|c| c.is_control() || "<>:\"/\\|?*".contains(c)),
            "Condition {condition:?} contains a character that cannot appear in an output filename"
        );
    }
    ensure!(
        !args.feature.is_empty() && !args.attribute.is_empty(),
        "--feature and --attribute cannot be empty"
    );

    let content =
        fs::read_to_string(manifest).with_context(|| format!("Reading {}", manifest.display()))?;
    let mut lines = content
        .trim_start_matches('\u{feff}')
        .lines()
        .enumerate()
        .filter(|(_, line)| !line.trim().is_empty() && !line.starts_with('#'));
    ensure!(
        lines
            .next()
            .is_some_and(|(_, line)| line == "sample\tcondition\tr1\tr2"),
        "Expected the TSV header: sample, condition, r1, r2"
    );

    let mut ids = HashSet::new();
    let mut samples = Vec::new();
    for (index, line) in lines {
        let columns: Vec<_> = line.split('\t').collect();
        ensure!(
            columns.len() == 4 && columns.iter().all(|column| !column.trim().is_empty()),
            "Line {} of samples.tsv must have four nonempty tab-separated fields",
            index + 1
        );
        let id = columns[0];
        ensure!(
            id.bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"._-".contains(&b)),
            "Invalid sample identifier {id:?}; use letters, numbers, dot, underscore or hyphen"
        );
        ensure!(
            ids.insert(id.to_ascii_lowercase()),
            "Duplicate sample identifier {id:?} (case insensitive)"
        );
        samples.push(Sample {
            id: id.to_owned(),
            condition: columns[1].to_owned(),
            r1: input_file(Path::new(columns[2]), &format!("{id} R1 FASTQ"))?,
            r2: input_file(Path::new(columns[3]), &format!("{id} R2 FASTQ"))?,
        });
    }
    for condition in [&args.condition_a, &args.condition_b] {
        let replicates = samples
            .iter()
            .filter(|sample| &sample.condition == condition)
            .count();
        ensure!(
            replicates >= 2,
            "Condition {condition:?} requires at least two samples; found {replicates}"
        );
    }
    Ok((samples, tag))
}

fn main() -> ExitCode {
    match workflow::run(Args::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("ERROR: {error:#}");
            ExitCode::FAILURE
        }
    }
}
