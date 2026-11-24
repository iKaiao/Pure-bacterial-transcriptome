use crate::process::Runner;
use crate::state::{State, atomic_write, nonempty};
use crate::{Args, Sample, input_file, read_samples};
use crate::{bam, deseq, read_qc};
use anyhow::{Context, Result, bail, ensure};
use flate2::read::MultiGzDecoder;
use std::env;
use std::fs::{self, File};
use std::io::{self, Write};
use std::path::{Path, PathBuf, absolute};
use tempfile::{NamedTempFile, TempDir, tempdir_in};

const FIXED_THREADS: &str = "8";
const INDEX_PARTS: [&str; 6] = ["1", "2", "3", "4", "rev.1", "rev.2"];

pub fn run(args: Args) -> Result<()> {
    let genome = input_file(&args.genome, "Reference FASTA")?;
    let manifest = input_file(&args.samples, "Sample manifest")?;
    let (samples, tag) = read_samples(&args, &manifest)?;
    let adapter = find_adapter(args.adapter.as_deref())?;
    let supplied_gff = args
        .gff
        .as_deref()
        .map(|path| input_file(path, "GFF3"))
        .transpose()?;
    let out = absolute(&args.out)?;
    for directory in [
        "reference",
        "clean",
        "bam",
        "qc",
        "counts",
        "deseq2",
        "logs",
    ] {
        fs::create_dir_all(out.join(directory))?;
    }
    let state = State::open(&out)?;
    let mut tools = vec![
        "trimmomatic",
        "bowtie2",
        "bowtie2-build",
        "bowtie2-inspect",
        "featureCounts",
    ];
    if supplied_gff.is_none() {
        tools.push("prokka");
    }
    let runner = Runner::new(&out, &tools)?;
    let pipeline = Pipeline { out, state, runner };
    pipeline.runner.message(&format!(
        "Reference: {}\nContrast: {} / {}; featureCounts strand {}",
        genome.display(),
        args.condition_a,
        args.condition_b,
        args.strand,
    ))?;
    pipeline.write_manifest(&samples)?;
    let reference = pipeline.prepare_reference(&genome)?;
    let gff = match supplied_gff {
        Some(path) => path,
        None => pipeline.annotate(&reference)?,
    };
    pipeline
        .runner
        .message(&format!("Annotation: {}", gff.display()))?;

    let paired_reads: Vec<_> = samples
        .iter()
        .map(|sample| pipeline.trim(sample, &adapter))
        .collect::<Result<_>>()?;
    if !args.no_read_qc {
        pipeline.read_qc(&paired_reads)?;
    }
    let mut bams = Vec::with_capacity(samples.len());
    for (sample, paired) in samples.iter().zip(&paired_reads) {
        bams.push(pipeline.map(sample, paired, args.threads)?);
    }
    let counts = pipeline.count(&args, &gff, &bams)?;
    pipeline.differential_expression(&args, &counts, &samples, &tag)?;
    pipeline.runner.message(&format!(
        "Workflow complete\nDESeq2 results: {}",
        pipeline
            .out
            .join("deseq2")
            .join(format!("DESeq2_{tag}.csv"))
            .display(),
    ))
}

struct Pipeline {
    out: PathBuf,
    state: State,
    runner: Runner,
}

impl Pipeline {
    fn stage(
        &self,
        name: &str,
        files: &[&Path],
        tools: &[&str],
        options: Vec<String>,
        outputs: &[PathBuf],
        action: impl FnOnce() -> Result<()>,
    ) -> Result<()> {
        let mut inputs = files.to_vec();
        inputs.extend(tools.iter().map(|tool| self.runner.tool(tool)));
        if self.state.fresh(name, &inputs, &options, outputs) {
            return self.runner.message(&format!("{name}: already complete"));
        }
        self.runner.message(name)?;
        self.state.invalidate(name)?;
        action().with_context(|| format!("Step {name}"))?;
        self.state.record(name, &inputs, options, outputs)
    }

    fn write_manifest(&self, samples: &[Sample]) -> Result<()> {
        let path = self.out.join("samples.tsv");
        let mut manifest = b"sample\tcondition\tr1\tr2\n".to_vec();
        for sample in samples {
            writeln!(
                manifest,
                "{}\t{}\t{}\t{}",
                sample.id,
                sample.condition,
                sample.r1.display(),
                sample.r2.display()
            )?;
        }
        if !fs::read(&path).is_ok_and(|existing| existing == manifest) {
            atomic_write(&path, &manifest)?;
        }
        Ok(())
    }

    fn prepare_reference(&self, genome: &Path) -> Result<PathBuf> {
        let reference = self.out.join("reference/reference.fasta");
        self.stage(
            "reference",
            &[genome],
            &[],
            vec![],
            std::slice::from_ref(&reference),
            || {
                let mut copy = NamedTempFile::new_in(self.out.join("reference"))?;
                io::copy(&mut File::open(genome)?, copy.as_file_mut())?;
                copy.persist(&reference)?;
                Ok(())
            },
        )?;
        let prefix = self.out.join("reference/reference");
        let inputs = [&*reference, self.runner.tool("bowtie2-build")];
        let current = complete_index(&prefix);
        if current
            .as_ref()
            .is_some_and(|files| self.state.fresh("bowtie2-index", &inputs, &[], files))
        {
            self.runner.message("bowtie2-index: already complete")?;
            return Ok(reference);
        }
        self.runner.message("bowtie2-index")?;
        self.state.invalidate("bowtie2-index")?;
        let working = tempdir_in(self.out.join("reference"))?;
        let index_prefix = working.path().join("reference");
        self.runner.run(
            self.runner
                .command("bowtie2-build")
                .arg(&reference)
                .arg(&index_prefix),
            "bowtie2-build.log",
            None,
        )?;
        self.runner.run(
            self.runner
                .command("bowtie2-inspect")
                .arg("-s")
                .arg(&index_prefix),
            "bowtie2-inspect.log",
            None,
        )?;
        let generated = complete_index(&index_prefix).context("Bowtie2 index is incomplete")?;
        for file in index_paths(&prefix, "bt2")
            .into_iter()
            .chain(index_paths(&prefix, "bt2l"))
        {
            if file.exists() {
                fs::remove_file(file)?;
            }
        }
        for part in generated {
            fs::rename(
                &part,
                prefix.with_file_name(part.file_name().context("Missing index file name")?),
            )?;
        }
        let index = complete_index(&prefix).context("Bowtie2 index was not published")?;
        self.state
            .record("bowtie2-index", &inputs, vec![], &index)?;
        Ok(reference)
    }

    fn annotate(&self, reference: &Path) -> Result<PathBuf> {
        let directory = self.out.join("prokka_annotation");
        fs::create_dir_all(&directory)?;
        let gff = directory.join("annotation.gff");
        self.stage(
            "prokka",
            &[reference],
            &["prokka"],
            vec![],
            std::slice::from_ref(&gff),
            || {
                let temporary = tempdir_in(&directory)?;
                let result = temporary.path().join("annotation");
                let tool_dir = self
                    .runner
                    .tool("prokka")
                    .parent()
                    .context("Missing Prokka executable directory")?;
                let inherited_path = env::var_os("PATH").unwrap_or_default();
                // Prokka uses env perl and helper programs from its own Conda environment.
                let prokka_path = env::join_paths(
                    std::iter::once(tool_dir.to_owned()).chain(env::split_paths(&inherited_path)),
                )?;
                self.runner.run(
                    self.runner
                        .command("prokka")
                        .env("PATH", prokka_path)
                        .arg("--outdir")
                        .arg(&result)
                        .args(["--prefix", "annotation", "--force", "--cpus", FIXED_THREADS])
                        .arg(reference),
                    "prokka.log",
                    None,
                )?;
                ensure!(
                    nonempty(&result.join("annotation.gff")),
                    "Prokka did not produce annotation.gff"
                );
                for entry in fs::read_dir(result)? {
                    let entry = entry?;
                    fs::rename(entry.path(), directory.join(entry.file_name()))?;
                }
                Ok(())
            },
        )?;
        Ok(gff)
    }

    fn trim(&self, sample: &Sample, adapter: &Path) -> Result<[PathBuf; 2]> {
        let clean = self.out.join("clean");
        let names = ["R1.paired", "R1.unpaired", "R2.paired", "R2.unpaired"]
            .map(|read| format!("{}.{read}.fq.gz", sample.id));
        let targets = names.each_ref().map(|name| clean.join(name));
        self.stage(
            &format!("trim-{}", sample.id),
            &[&sample.r1, &sample.r2, adapter],
            &["trimmomatic"],
            vec![],
            &targets,
            || {
                let temporary = tempdir_in(&clean)?;
                // ILLUMINACLIP uses colons as separators, including on Windows.
                fs::copy(adapter, temporary.path().join("adapter.fa"))?;
                self.runner.run(
                    self.runner
                        .command("trimmomatic")
                        .current_dir(temporary.path())
                        .env("_JAVA_OPTIONS", "-Xmx8g")
                        .args(["PE", "-threads", FIXED_THREADS, "-phred33"])
                        .arg(&sample.r1)
                        .arg(&sample.r2)
                        .args(&names)
                        .args([
                            "ILLUMINACLIP:adapter.fa:2:30:10",
                            "LEADING:3",
                            "TRAILING:3",
                            "SLIDINGWINDOW:4:15",
                            "MINLEN:36",
                        ]),
                    &format!("{}.trimmomatic.log", sample.id),
                    None,
                )?;
                for name in &names {
                    let file = temporary.path().join(name);
                    ensure!(
                        nonempty(&file),
                        "Empty Trimmomatic output: {}",
                        file.display()
                    );
                    io::copy(
                        &mut MultiGzDecoder::new(File::open(&file)?),
                        &mut io::sink(),
                    )
                    .with_context(|| format!("Invalid gzip output {}", file.display()))?;
                }
                publish(&temporary, &targets)
            },
        )?;
        Ok([targets[0].clone(), targets[2].clone()])
    }

    fn read_qc(&self, paired_reads: &[[PathBuf; 2]]) -> Result<()> {
        let qc = self.out.join("qc");
        let inputs: Vec<_> = paired_reads
            .iter()
            .flatten()
            .map(PathBuf::as_path)
            .collect();
        let outputs = read_qc::REPORTS.map(|name| qc.join(name));
        self.stage("read-qc", &inputs, &[], vec![], &outputs, || {
            let temporary = tempdir_in(&qc)?;
            read_qc::report(temporary.path(), paired_reads)?;
            publish(&temporary, &outputs)
        })
    }

    fn map(&self, sample: &Sample, paired: &[PathBuf; 2], threads: u32) -> Result<PathBuf> {
        let directory = self.out.join("bam");
        let name = format!("{}.sorted.bam", sample.id);
        let index_name = format!("{name}.bai");
        let bam = directory.join(&name);
        let bai = directory.join(&index_name);
        let suffix = format!("{}.bowtie2.log", sample.id);
        let prefix = self.out.join("reference/reference");
        let index = complete_index(&prefix).context("Bowtie2 index is incomplete")?;
        let inputs: Vec<_> = paired.iter().chain(&index).map(PathBuf::as_path).collect();
        self.stage(
            &format!("map-{}", sample.id),
            &inputs,
            &["bowtie2"],
            vec![threads.to_string()],
            &[bam.clone(), bai.clone()],
            || {
                let temporary = tempdir_in(&directory)?;
                let mut align = self.runner.command("bowtie2");
                align
                    .args(["--very-sensitive", "-p"])
                    .arg(threads.to_string())
                    .arg("-x")
                    .arg(&prefix)
                    .arg("-1")
                    .arg(&paired[0])
                    .arg("-2")
                    .arg(&paired[1]);
                self.runner.read_stdout(&mut align, &suffix, |stream| {
                    bam::sort_and_index(
                        stream,
                        &temporary.path().join(&name),
                        &temporary.path().join(&index_name),
                    )
                })?;
                publish(&temporary, &[bam.clone(), bai.clone()])
            },
        )?;
        let qc = self.out.join("qc");
        let report_names = ["flagstat", "idxstats"].map(|tool| format!("{}.{tool}.txt", sample.id));
        let reports = report_names.each_ref().map(|name| qc.join(name));
        self.stage(
            &format!("bam-qc-{}", sample.id),
            &[&bam, &bai],
            &[],
            vec![],
            &reports,
            || {
                let temporary = tempdir_in(&qc)?;
                bam::stats(
                    &bam,
                    &temporary.path().join(&report_names[0]),
                    &temporary.path().join(&report_names[1]),
                )?;
                publish(&temporary, &reports)
            },
        )?;
        if let Ok(log) = fs::read_to_string(self.runner.log_path(&suffix)) {
            for line in log
                .lines()
                .filter(|line| line.contains("overall alignment rate"))
            {
                self.runner.message(&format!("{}: {line}", sample.id))?;
            }
        }
        Ok(bam)
    }

    fn count(&self, args: &Args, gff: &Path, bams: &[PathBuf]) -> Result<PathBuf> {
        let directory = self.out.join("counts");
        let filename = format!("raw_counts.s{}.txt", args.strand);
        let count_file = directory.join(&filename);
        let summary = directory.join(format!("{filename}.summary"));
        let inputs: Vec<_> = std::iter::once(gff)
            .chain(bams.iter().map(PathBuf::as_path))
            .collect();
        let options = vec![args.feature.clone(), args.attribute.clone()];
        self.stage(
            &format!("featurecounts-s{}", args.strand),
            &inputs,
            &["featureCounts"],
            options,
            &[count_file.clone(), summary.clone()],
            || {
                let temporary = tempdir_in(&directory)?;
                let relative_bams = bams
                    .iter()
                    .map(|path| path.strip_prefix(&self.out))
                    .collect::<Result<Vec<_>, _>>()?;
                self.runner.run(
                    self.runner
                        .command("featureCounts")
                        .args([
                            "-T",
                            FIXED_THREADS,
                            "-p",
                            "--countReadPairs",
                            "-B",
                            "-C",
                            "-s",
                        ])
                        .arg(args.strand.to_string())
                        .arg("-a")
                        .arg(gff)
                        .arg("-t")
                        .arg(&args.feature)
                        .arg("-g")
                        .arg(&args.attribute)
                        .arg("-o")
                        .arg(temporary.path().join(filename))
                        .args(relative_bams),
                    "featureCounts.log",
                    None,
                )?;
                publish(&temporary, &[count_file.clone(), summary])
            },
        )?;
        Ok(count_file)
    }

    fn differential_expression(
        &self,
        args: &Args,
        counts: &Path,
        samples: &[Sample],
        tag: &str,
    ) -> Result<()> {
        self.runner.message("DESeq2")?;
        let directory = self.out.join("deseq2");
        let temporary = tempdir_in(&directory)?;
        let log_path = self.runner.log_path("DESeq2.log");
        let mut log = File::create(&log_path)?;
        let analysis = deseq::run(args, counts, samples, tag, temporary.path(), &mut log);
        if let Err(error) = &analysis {
            let _ = writeln!(log, "{error:#}");
        }
        analysis.with_context(|| {
            format!("Differential expression failed. See {}", log_path.display())
        })?;
        let mut files = [
            "raw_counts_matrix.csv",
            "library_sizes_raw_counts.csv",
            "sample_metadata.csv",
            "DESeq2_normalized_counts.csv",
            "DESeq2_size_factors.csv",
            "PCA.pdf",
            "sample_correlation.csv",
            "sample_correlation.pdf",
            "MA_plot.pdf",
            "DESeq2_summary.txt",
        ]
        .map(|name| directory.join(name))
        .to_vec();
        files.extend(
            [
                format!("DESeq2_{tag}.csv"),
                format!("volcano_source_data_{tag}.csv"),
                format!("volcano_DESeq2_{tag}.pdf"),
                format!("volcano_DESeq2_{tag}.png"),
            ]
            .map(|name| directory.join(name)),
        );
        publish(&temporary, &files)
    }
}

fn publish(temporary: &TempDir, targets: &[PathBuf]) -> Result<()> {
    let sources = targets
        .iter()
        .map(|target| {
            let filename = target.file_name().context("Missing output filename")?;
            let path = temporary.path().join(filename);
            ensure!(
                nonempty(&path),
                "Missing or empty output: {}",
                target.display()
            );
            Ok(path)
        })
        .collect::<Result<Vec<_>>>()?;
    for (source, target) in sources.iter().zip(targets) {
        fs::rename(source, target).with_context(|| format!("Publishing {}", target.display()))?;
    }
    Ok(())
}

fn index_paths(prefix: &Path, extension: &str) -> Vec<PathBuf> {
    INDEX_PARTS
        .iter()
        .map(|part| prefix.with_extension(format!("{part}.{extension}")))
        .collect()
}

fn complete_index(prefix: &Path) -> Option<Vec<PathBuf>> {
    ["bt2", "bt2l"]
        .into_iter()
        .map(|extension| index_paths(prefix, extension))
        .find(|files| files.iter().all(|path| nonempty(path)))
}

fn find_adapter(explicit: Option<&Path>) -> Result<PathBuf> {
    if let Some(path) = explicit {
        return input_file(path, "Adapter FASTA");
    }
    if let Some(prefix) = env::var_os("CONDA_PREFIX") {
        let share = PathBuf::from(prefix).join("share");
        let direct = share.join("adapters/TruSeq3-PE.fa");
        if nonempty(&direct) {
            return input_file(&direct, "Adapter FASTA");
        }
        if let Ok(entries) = fs::read_dir(share) {
            let mut matches: Vec<_> = entries
                .flatten()
                .map(|entry| entry.path().join("adapters/TruSeq3-PE.fa"))
                .filter(|path| nonempty(path))
                .collect();
            matches.sort();
            if let Some(path) = matches.first() {
                return input_file(path, "Adapter FASTA");
            }
        }
    }
    bail!("TruSeq3-PE.fa was not found in the active Conda environment; pass --adapter FILE")
}
