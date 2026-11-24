use anyhow::{Context, Result, bail};
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{ChildStdout, Command, ExitStatus, Stdio};

pub struct Runner {
    cwd: PathBuf,
    logs: PathBuf,
    pipeline_log: File,
    tools: HashMap<&'static str, PathBuf>,
}

impl Runner {
    pub fn new(cwd: &Path, required: &[&'static str]) -> Result<Self> {
        let mut tools = HashMap::new();
        let mut missing = Vec::new();
        for &program in required {
            match which::which(program) {
                Ok(path) => {
                    tools.insert(program, path);
                }
                Err(_) => missing.push(program),
            }
        }
        if !missing.is_empty() {
            if missing == ["prokka"] {
                bail!(
                    "Prokka is required for automatic annotation. Add it to PATH or provide --gff."
                );
            }
            bail!(
                "Missing tools on PATH: {}. Activate BactRNAseq or install the Conda package.",
                missing.join(", ")
            );
        }
        let logs = cwd.join("logs");
        let pipeline_log = OpenOptions::new()
            .create(true)
            .append(true)
            .open(logs.join("pipeline.log"))?;
        Ok(Self {
            cwd: cwd.to_owned(),
            logs,
            pipeline_log,
            tools,
        })
    }

    pub fn tool(&self, program: &str) -> &Path {
        &self.tools[program]
    }

    pub fn message(&self, message: &str) -> Result<()> {
        println!("{message}");
        writeln!(&self.pipeline_log, "{message}")?;
        Ok(())
    }

    pub fn command(&self, program: &str) -> Command {
        let mut command = Command::new(self.tool(program));
        command.current_dir(&self.cwd);
        command
    }

    pub fn log_path(&self, log: &str) -> PathBuf {
        self.logs.join(log)
    }

    pub fn run(&self, command: &mut Command, log: &str, output: Option<&Path>) -> Result<()> {
        let log_path = self.log_path(log);
        let log_file = File::create(&log_path)?;
        let stdout = match output {
            Some(path) => File::create(path)?,
            None => log_file.try_clone()?,
        };
        let status = command
            .stdin(Stdio::null())
            .stdout(stdout)
            .stderr(log_file)
            .status()
            .with_context(|| format!("Could not start {command:?}"))?;
        check_status(status, command, &log_path)
    }

    pub fn read_stdout(
        &self,
        command: &mut Command,
        log: &str,
        read: impl FnOnce(ChildStdout) -> Result<()>,
    ) -> Result<()> {
        let log_path = self.log_path(log);
        let mut child = command
            .stdin(Stdio::null())
            .stdout(Stdio::piped())
            .stderr(File::create(&log_path)?)
            .spawn()
            .with_context(|| format!("Could not start {command:?}"))?;
        let stream = child.stdout.take().expect("stdout was piped");
        let result = read(stream);
        if result.is_err() {
            let _ = child.kill();
        }
        let status = child.wait()?;
        result.with_context(|| {
            format!(
                "Reading output from {command:?}; see {}",
                log_path.display()
            )
        })?;
        check_status(status, command, &log_path)
    }
}

fn check_status(status: ExitStatus, command: &Command, log: &Path) -> Result<()> {
    if !status.success() {
        bail!("{command:?} failed with {status}. See {}", log.display());
    }
    Ok(())
}
