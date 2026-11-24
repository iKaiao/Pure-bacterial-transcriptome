use anyhow::{Context, Result, ensure};
use flate2::read::MultiGzDecoder;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

pub const REPORTS: [&str; 4] = [
    "read_qc_report.html",
    "read_qc_summary.tsv",
    "read_qc_by_cycle.tsv",
    "read_qc_lengths.tsv",
];

#[derive(Default)]
struct Cycle {
    composition: [u64; 5],
    quality_sum: u64,
    q30: u64,
}

impl Cycle {
    fn observe(&mut self, base: u8, quality: u8) {
        let index = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 4,
        };
        self.composition[index] += 1;
        self.quality_sum += u64::from(quality);
        self.q30 += u64::from(quality >= 30);
    }

    fn bases(&self) -> u64 {
        self.composition.iter().sum()
    }
}

#[derive(Default)]
struct ReadQc {
    reads: u64,
    total: Cycle,
    cycles: Vec<Cycle>,
    lengths: BTreeMap<usize, u64>,
}

pub fn report(directory: &Path, paired_reads: &[[PathBuf; 2]]) -> Result<()> {
    let [html, summaries, cycles, lengths] = REPORTS.map(|name| directory.join(name));
    let mut html = BufWriter::new(File::create(html)?);
    let mut summaries = BufWriter::new(File::create(summaries)?);
    let mut cycles = BufWriter::new(File::create(cycles)?);
    let mut lengths = BufWriter::new(File::create(lengths)?);
    writeln!(
        html,
        r#"<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Trimmed read quality</title>
<style>
body{{font:1rem system-ui,sans-serif;color:#192731;max-width:70rem;margin:2rem auto;padding:0 1rem}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:.55rem .7rem;border-bottom:1px solid #cad5da;text-align:right}}
th:first-child,td:first-child{{text-align:left}} thead{{background:#eef5f6}}
a{{color:#00576b}} .table{{overflow-x:auto}}
</style>
<h1>Trimmed read quality</h1>
<p><a href="read_qc_summary.tsv">Summary TSV</a> ·
<a href="read_qc_by_cycle.tsv">Per-cycle TSV</a> ·
<a href="read_qc_lengths.tsv">Read lengths TSV</a></p>
<div class="table"><table><thead><tr>
<th scope="col">Read file</th><th scope="col">Reads</th><th scope="col">Bases</th>
<th scope="col">Mean Phred</th><th scope="col">Q30 %</th><th scope="col">GC %</th>
<th scope="col">Ambiguous %</th></tr></thead><tbody>"#
    )?;
    writeln!(
        summaries,
        "read_file\treads\tbases\tmean_phred\tq30_percent\tgc_percent\tambiguous_percent"
    )?;
    writeln!(
        cycles,
        "read_file\tcycle\tbases\tmean_phred\tq30_percent\tA\tC\tG\tT\tambiguous"
    )?;
    writeln!(lengths, "read_file\tlength\treads")?;

    for pair in paired_reads {
        let mut read_counts = Vec::with_capacity(2);
        for path in pair {
            let name = path
                .file_name()
                .and_then(|name| name.to_str())
                .context("Invalid read QC filename")?;
            let qc = inspect(path).with_context(|| format!("Reading FASTQ {}", path.display()))?;
            read_counts.push(qc.reads);
            let bases = qc.total.bases();
            let [_, cytosine, guanine, _, ambiguous] = qc.total.composition;
            let phred = fraction(qc.total.quality_sum, bases, 1.0);
            let q30 = fraction(qc.total.q30, bases, 100.0);
            let gc = fraction(cytosine + guanine, bases, 100.0);
            let unknown = fraction(ambiguous, bases, 100.0);
            let reads = qc.reads;
            writeln!(
                summaries,
                "{name}\t{reads}\t{bases}\t{phred:.2}\t{q30:.2}\t{gc:.2}\t{unknown:.2}"
            )?;
            writeln!(
                html,
                "<tr><th scope=\"row\">{name}</th><td>{reads}</td><td>{bases}</td><td>{phred:.2}</td><td>{q30:.2}</td><td>{gc:.2}</td><td>{unknown:.2}</td></tr>"
            )?;
            for (index, cycle) in qc.cycles.iter().enumerate() {
                let bases = cycle.bases();
                let phred = fraction(cycle.quality_sum, bases, 1.0);
                let q30 = fraction(cycle.q30, bases, 100.0);
                let [a, c, g, t, ambiguous] = cycle.composition;
                writeln!(
                    cycles,
                    "{name}\t{}\t{bases}\t{phred:.2}\t{q30:.2}\t{a}\t{c}\t{g}\t{t}\t{ambiguous}",
                    index + 1
                )?;
            }
            for (length, count) in qc.lengths {
                writeln!(lengths, "{name}\t{length}\t{count}")?;
            }
        }
        ensure!(
            read_counts[0] == read_counts[1],
            "Trimmed mate counts differ for {}",
            pair[0].display()
        );
    }
    writeln!(html, "</tbody></table></div></html>")?;
    for mut writer in [html, summaries, cycles, lengths] {
        writer.flush()?;
    }
    Ok(())
}

fn inspect(path: &Path) -> Result<ReadQc> {
    let mut reader = BufReader::new(MultiGzDecoder::new(File::open(path)?));
    let mut qc = ReadQc::default();
    let mut header = Vec::new();
    let mut seq = Vec::new();
    let mut sep = Vec::new();
    let mut qual = Vec::new();

    while line(&mut reader, &mut header)? {
        let record = qc.reads + 1;
        ensure!(
            line(&mut reader, &mut seq)?
                && line(&mut reader, &mut sep)?
                && line(&mut reader, &mut qual)?,
            "Truncated FASTQ record {record}"
        );
        ensure!(
            header.starts_with(b"@")
                && sep.starts_with(b"+")
                && seq.len() == qual.len()
                && !seq.is_empty(),
            "Malformed FASTQ record {record}"
        );
        qc.reads += 1;
        *qc.lengths.entry(seq.len()).or_insert(0) += 1;
        if qc.cycles.len() < seq.len() {
            qc.cycles.resize_with(seq.len(), Cycle::default);
        }
        for ((&base, &quality), cycle) in seq.iter().zip(&qual).zip(&mut qc.cycles) {
            ensure!(
                (33..=126).contains(&quality),
                "Invalid Phred+33 byte in FASTQ record {record}"
            );
            cycle.observe(base, quality - 33);
            qc.total.observe(base, quality - 33);
        }
    }
    Ok(qc)
}

fn line(reader: &mut impl BufRead, buffer: &mut Vec<u8>) -> Result<bool> {
    buffer.clear();
    let read = reader.read_until(b'\n', buffer)?;
    if buffer.last() == Some(&b'\n') {
        buffer.pop();
    }
    if buffer.last() == Some(&b'\r') {
        buffer.pop();
    }
    Ok(read != 0)
}

fn fraction(numerator: u64, denominator: u64, scale: f64) -> f64 {
    if denominator == 0 {
        0.0
    } else {
        scale * numerator as f64 / denominator as f64
    }
}
