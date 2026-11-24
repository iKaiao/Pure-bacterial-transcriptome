//! SAM/BAM processing for the alignments returned by Bowtie2.
use anyhow::{Context, Result, ensure};
use noodles_bam as bam;
use noodles_sam::{self as sam, alignment::RecordBuf, alignment::io::Write as _};
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

const RUN_BYTES: usize = 16 * 1024 * 1024;
const MERGE_FAN_IN: usize = 32;
type Coordinate = (usize, usize, bool);

pub fn sort_and_index(input: impl Read, destination: &Path, index: &Path) -> Result<()> {
    coordinate_sort(input, destination, RUN_BYTES)?;
    let entries = bam::fs::index(destination).context("Indexing the sorted BAM")?;
    bam::bai::fs::write(index, &entries)?;
    Ok(())
}

fn coordinate_sort(input: impl Read, destination: &Path, run_limit: usize) -> Result<()> {
    let mut reader = sam::io::Reader::new(BufReader::new(input));
    let mut header = reader
        .read_header()
        .context("Reading the Bowtie2 SAM header")?;
    ensure!(
        !header.reference_sequences().is_empty(),
        "SAM header contains no reference sequences"
    );
    header
        .header_mut()
        .get_or_insert_with(Default::default)
        .other_fields_mut()
        .insert(
            sam::header::record::value::map::header::tag::SORT_ORDER,
            "coordinate".into(),
        );

    let directory = destination
        .parent()
        .context("Missing BAM output directory")?;
    let mut runs = Vec::new();
    let mut records = Vec::new();
    let mut bytes = 0;
    loop {
        let mut record = RecordBuf::default();
        let size = reader
            .read_record_buf(&header, &mut record)
            .context("Reading Bowtie2 SAM")?;
        if size == 0 {
            break;
        }
        records.push(record);
        bytes += size;
        if bytes >= run_limit {
            let path = directory.join(format!("run-0-{}.bam", runs.len()));
            write_sorted(&header, &mut records, &path)?;
            runs.push(path);
            bytes = 0;
        }
    }
    if runs.is_empty() {
        return write_sorted(&header, &mut records, destination);
    }
    if !records.is_empty() {
        let path = directory.join(format!("run-0-{}.bam", runs.len()));
        write_sorted(&header, &mut records, &path)?;
        runs.push(path);
    }

    let mut pass = 1;
    while runs.len() > MERGE_FAN_IN {
        let mut merged = Vec::new();
        for (group, paths) in runs.chunks(MERGE_FAN_IN).enumerate() {
            let path = directory.join(format!("run-{pass}-{group}.bam"));
            merge(&header, paths, &path)?;
            for old in paths {
                fs::remove_file(old)?;
            }
            merged.push(path);
        }
        runs = merged;
        pass += 1;
    }
    merge(&header, &runs, destination)
}

fn write_sorted(header: &sam::Header, records: &mut Vec<RecordBuf>, path: &Path) -> Result<()> {
    // A stable sort preserves Bowtie2's order for identical coordinates and strands.
    records.sort_by_key(|record| {
        (
            record.reference_sequence_id().unwrap_or(usize::MAX),
            record.alignment_start().map_or(0, usize::from),
            record.flags().is_reverse_complemented(),
        )
    });
    let mut writer = bam::io::Writer::new(File::create(path)?);
    writer.write_header(header)?;
    for record in records.drain(..) {
        writer.write_alignment_record(header, &record)?;
    }
    writer.try_finish()?;
    Ok(())
}

fn coordinate(record: &bam::Record) -> Result<Coordinate> {
    Ok((
        record
            .reference_sequence_id()
            .transpose()?
            .unwrap_or(usize::MAX),
        record.alignment_start().transpose()?.map_or(0, usize::from),
        record.flags().is_reverse_complemented(),
    ))
}

fn merge(header: &sam::Header, sources: &[PathBuf], destination: &Path) -> Result<()> {
    let mut readers = sources
        .iter()
        .map(|source| {
            let mut reader = bam::io::Reader::new(File::open(source)?);
            reader.read_header()?;
            Ok(reader)
        })
        .collect::<Result<Vec<_>>>()?;
    let mut writer = bam::io::Writer::new(File::create(destination)?);
    writer.write_header(header)?;
    // Run order is the secondary key so ties stay stable across merge passes.
    let mut queue = BinaryHeap::new();
    let mut records = vec![bam::Record::default(); readers.len()];
    for (index, (reader, record)) in readers.iter_mut().zip(&mut records).enumerate() {
        if reader.read_record(record)? != 0 {
            queue.push(Reverse((coordinate(record)?, index)));
        }
    }
    while let Some(Reverse((_, index))) = queue.pop() {
        let record = &mut records[index];
        writer.write_record(header, record)?;
        if readers[index].read_record(record)? != 0 {
            queue.push(Reverse((coordinate(record)?, index)));
        }
    }
    writer.try_finish()?;
    Ok(())
}

#[derive(Default)]
struct AlignmentCounts {
    total: u64,
    primary: u64,
    secondary: u64,
    supplementary: u64,
    duplicates: u64,
    primary_duplicates: u64,
    mapped: u64,
    primary_mapped: u64,
    paired: u64,
    read1: u64,
    read2: u64,
    properly_paired: u64,
    both_mapped: u64,
    singletons: u64,
    different_reference: u64,
    different_reference_mapq5: u64,
}

impl AlignmentCounts {
    fn observe(&mut self, record: &bam::Record, reference: Option<usize>) -> Result<()> {
        let flags = record.flags();
        let mapped = !flags.is_unmapped();
        self.total += 1;
        self.mapped += u64::from(mapped);
        self.duplicates += u64::from(flags.is_duplicate());
        if flags.is_secondary() {
            self.secondary += 1;
        } else if flags.is_supplementary() {
            self.supplementary += 1;
        } else {
            self.primary += 1;
            self.primary_mapped += u64::from(mapped);
            self.primary_duplicates += u64::from(flags.is_duplicate());
            if flags.is_segmented() {
                self.paired += 1;
                self.read1 += u64::from(flags.is_first_segment());
                self.read2 += u64::from(flags.is_last_segment());
                self.properly_paired += u64::from(mapped && flags.is_properly_segmented());
                self.singletons += u64::from(mapped && flags.is_mate_unmapped());
                if mapped && !flags.is_mate_unmapped() {
                    self.both_mapped += 1;
                    if reference != record.mate_reference_sequence_id().transpose()? {
                        self.different_reference += 1;
                        let mapq = record.mapping_quality().map_or(255, u8::from);
                        self.different_reference_mapq5 += u64::from(mapq >= 5);
                    }
                }
            }
        }
        Ok(())
    }
}

pub fn stats(input: &Path, flagstat: &Path, idxstats: &Path) -> Result<()> {
    let mut reader = bam::io::Reader::new(File::open(input)?);
    let header = reader.read_header()?;
    let references = header.reference_sequences();
    let mut counts = vec![[0_u64; 2]; references.len()];
    let mut unplaced = 0;
    let mut flags: [AlignmentCounts; 2] = Default::default();
    let mut record = bam::Record::default();
    while reader.read_record(&mut record)? != 0 {
        let ref_id = record.reference_sequence_id().transpose()?;
        flags[usize::from(record.flags().is_qc_fail())].observe(&record, ref_id)?;
        match ref_id {
            Some(id) => {
                let count = counts
                    .get_mut(id)
                    .context("Reference ID is missing from the BAM header")?;
                count[usize::from(record.flags().is_unmapped())] += 1;
            }
            None => unplaced += 1,
        }
    }
    let mut idx = BufWriter::new(File::create(idxstats)?);
    for ((name, reference), [mapped, unmapped]) in references.iter().zip(counts) {
        writeln!(idx, "{name}\t{}\t{mapped}\t{unmapped}", reference.length())?;
    }
    writeln!(idx, "*\t0\t0\t{unplaced}")?;
    idx.flush()?;
    let mut report = BufWriter::new(File::create(flagstat)?);
    write_flags(&mut report, &flags[0], &flags[1])?;
    report.flush()?;
    Ok(())
}

fn write_flags(out: &mut impl Write, pass: &AlignmentCounts, fail: &AlignmentCounts) -> Result<()> {
    for (label, left, right, denominator) in [
        (
            "in total (QC-passed reads + QC-failed reads)",
            pass.total,
            fail.total,
            None,
        ),
        ("primary", pass.primary, fail.primary, None),
        ("secondary", pass.secondary, fail.secondary, None),
        (
            "supplementary",
            pass.supplementary,
            fail.supplementary,
            None,
        ),
        ("duplicates", pass.duplicates, fail.duplicates, None),
        (
            "primary duplicates",
            pass.primary_duplicates,
            fail.primary_duplicates,
            None,
        ),
        (
            "mapped",
            pass.mapped,
            fail.mapped,
            Some((pass.total, fail.total)),
        ),
        (
            "primary mapped",
            pass.primary_mapped,
            fail.primary_mapped,
            Some((pass.primary, fail.primary)),
        ),
        ("paired in sequencing", pass.paired, fail.paired, None),
        ("read1", pass.read1, fail.read1, None),
        ("read2", pass.read2, fail.read2, None),
        (
            "properly paired",
            pass.properly_paired,
            fail.properly_paired,
            Some((pass.paired, fail.paired)),
        ),
        (
            "with itself and mate mapped",
            pass.both_mapped,
            fail.both_mapped,
            None,
        ),
        (
            "singletons",
            pass.singletons,
            fail.singletons,
            Some((pass.paired, fail.paired)),
        ),
        (
            "with mate mapped to a different chr",
            pass.different_reference,
            fail.different_reference,
            None,
        ),
        (
            "with mate mapped to a different chr (mapQ>=5)",
            pass.different_reference_mapq5,
            fail.different_reference_mapq5,
            None,
        ),
    ] {
        write!(out, "{left} + {right} {label}")?;
        if let Some((total_pass, total_fail)) = denominator {
            write!(
                out,
                " ({} : {})",
                percent(left, total_pass),
                percent(right, total_fail)
            )?;
        }
        writeln!(out)?;
    }
    Ok(())
}

fn percent(count: u64, total: u64) -> String {
    if total == 0 {
        "N/A".to_owned()
    } else {
        format!("{:.2}%", f64::from(count as f32 / total as f32) * 100.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as _;

    #[test]
    fn merges_multiple_passes_and_builds_a_bai() -> Result<()> {
        let work = tempfile::tempdir()?;
        let mut sam = String::from("@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:1000\n");
        for position in (1..=40).rev() {
            writeln!(
                sam,
                "read{position}\t0\tref\t{position}\t60\t4M\t*\t0\t0\tACGT\tIIII"
            )?;
        }
        sam.push_str("unplaced\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII\n");
        let bam = work.path().join("sorted.bam");
        coordinate_sort(sam.as_bytes(), &bam, 1)?;
        bam::bai::fs::write(work.path().join("sorted.bam.bai"), &bam::fs::index(&bam)?)?;

        let flagstat = work.path().join("flags");
        let idxstats = work.path().join("index");
        stats(&bam, &flagstat, &idxstats)?;
        assert_eq!(
            fs::read_to_string(idxstats)?,
            "ref\t1000\t40\t0\n*\t0\t0\t1\n"
        );
        assert!(fs::read_to_string(flagstat)?.contains("40 + 0 mapped (97.56% : N/A)"));
        let mut reader = bam::io::Reader::new(File::open(bam)?);
        reader.read_header()?;
        for (index, record) in reader.records().enumerate() {
            let record = record?;
            let expected = (index < 40).then_some(index + 1);
            assert_eq!(
                record.alignment_start().transpose()?.map(usize::from),
                expected
            );
        }
        Ok(())
    }
}
