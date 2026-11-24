use anyhow::{Context, Result, bail, ensure};
use serde::{Deserialize, Serialize};
use std::fs::{self, File, TryLockError};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::time::SystemTime;
use tempfile::NamedTempFile;

#[derive(Deserialize, PartialEq, Serialize)]
struct FileStamp {
    path: PathBuf,
    bytes: u64,
    modified: SystemTime,
}

impl FileStamp {
    fn read(path: &Path) -> Result<Self> {
        let metadata = path
            .metadata()
            .with_context(|| format!("Reading metadata for {}", path.display()))?;
        ensure!(
            metadata.is_file() && metadata.len() > 0,
            "Missing or empty file: {}",
            path.display()
        );
        Ok(Self {
            path: path.to_owned(),
            bytes: metadata.len(),
            modified: metadata.modified()?,
        })
    }
}

#[derive(Deserialize, Serialize)]
struct Checkpoint {
    version: String,
    inputs: Vec<FileStamp>,
    options: Vec<String>,
    outputs: Vec<FileStamp>,
}

pub struct State {
    directory: PathBuf,
    _lock: File,
}

impl State {
    pub fn open(out: &Path) -> Result<Self> {
        let directory = out.join(".bact-rnaseq");
        fs::create_dir_all(&directory)?;
        let lock = File::options()
            .read(true)
            .write(true)
            .create(true)
            .truncate(false)
            .open(directory.join("run.lock"))?;
        match lock.try_lock() {
            Ok(()) => Ok(Self {
                directory,
                _lock: lock,
            }),
            Err(TryLockError::WouldBlock) => {
                bail!("Another Rust workflow is using {}", out.display())
            }
            Err(TryLockError::Error(error)) => {
                Err(error).context("Could not lock the output directory")
            }
        }
    }

    pub fn fresh(
        &self,
        step: &str,
        inputs: &[&Path],
        options: &[String],
        outputs: &[PathBuf],
    ) -> bool {
        let file = self.directory.join(format!("{step}.json"));
        let Ok(previous) = File::open(file).and_then(|file| {
            serde_json::from_reader::<_, Checkpoint>(file).map_err(io::Error::from)
        }) else {
            return false;
        };
        previous.version == env!("CARGO_PKG_VERSION")
            && previous.options == options
            && stamps(inputs.iter().copied()).is_ok_and(|current| current == previous.inputs)
            && stamps(outputs.iter().map(PathBuf::as_path))
                .is_ok_and(|current| current == previous.outputs)
    }

    pub fn invalidate(&self, step: &str) -> Result<()> {
        match fs::remove_file(self.directory.join(format!("{step}.json"))) {
            Err(error) if error.kind() != io::ErrorKind::NotFound => Err(error.into()),
            _ => Ok(()),
        }
    }

    pub fn record(
        &self,
        step: &str,
        inputs: &[&Path],
        options: Vec<String>,
        outputs: &[PathBuf],
    ) -> Result<()> {
        let checkpoint = Checkpoint {
            version: env!("CARGO_PKG_VERSION").to_owned(),
            inputs: stamps(inputs.iter().copied())?,
            options,
            outputs: stamps(outputs.iter().map(PathBuf::as_path))?,
        };
        let mut temporary = NamedTempFile::new_in(&self.directory)?;
        serde_json::to_writer(temporary.as_file_mut(), &checkpoint)?;
        temporary.persist(self.directory.join(format!("{step}.json")))?;
        Ok(())
    }
}

fn stamps<'a>(paths: impl Iterator<Item = &'a Path>) -> Result<Vec<FileStamp>> {
    paths.map(FileStamp::read).collect()
}

pub fn atomic_write(destination: &Path, data: &[u8]) -> Result<()> {
    let parent = destination
        .parent()
        .context("Missing output parent directory")?;
    let mut temporary = NamedTempFile::new_in(parent)?;
    temporary.write_all(data)?;
    temporary.persist(destination)?;
    Ok(())
}

pub fn nonempty(path: &Path) -> bool {
    path.metadata()
        .is_ok_and(|metadata| metadata.is_file() && metadata.len() > 0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checkpoints_require_complete_matching_outputs() -> Result<()> {
        let directory = tempfile::tempdir()?;
        let input = directory.path().join("input");
        let output = directory.path().join("output");
        fs::write(&input, "original")?;
        fs::write(&output, "result")?;
        let state = State::open(directory.path())?;
        let inputs = [input.as_path()];
        let outputs = [output.clone()];
        let options = vec!["strand=2".to_owned()];

        assert!(!state.fresh("counts", &inputs, &options, &outputs));
        state.record("counts", &inputs, options.clone(), &outputs)?;
        assert!(state.fresh("counts", &inputs, &options, &outputs));
        assert!(!state.fresh("counts", &inputs, &["strand=0".to_owned()], &outputs));
        fs::write(&input, "different input")?;
        assert!(!state.fresh("counts", &inputs, &options, &outputs));
        state.record("counts", &inputs, options.clone(), &outputs)?;
        fs::remove_file(output)?;
        assert!(!state.fresh("counts", &inputs, &options, &outputs));
        assert!(State::open(directory.path()).is_err());
        drop(state);
        State::open(directory.path())?;
        Ok(())
    }
}
