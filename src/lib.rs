#![deny(missing_docs, clippy::missing_docs_in_private_items)]
//! A python module for quickly and efficiently calculating a dynamic CNV profile from contigs continaed in a  BAM file.

use fnv::FnvHashMap;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use pyo3::prelude::*;
use rayon::prelude::*;
use rust_htslib::bam::{self, FetchDefinition, Read, Record};
use std::{error::Error, path::PathBuf};

/// Number of reads expected in given bin width on the genome
const READS_PER_BIN: usize = 100;
/// Genome bin size for binning read starts
const BIN_SIZE: usize = 1000;
/// Excpected ploidy
const PLOIDY: u16 = 2;

/// Dynamic result type for fancy unwrapping
pub type DynResult<T> = Result<T, Box<dyn Error + 'static>>;
/// Represents a BAM file with its path.
pub struct BamFile {
    /// Path to the BAM file, mainly used to convey intent.
    path: PathBuf,
}

impl BamFile {
    /// Creates a new `BamFile` instance.
    ///
    /// # Arguments
    ///
    /// * `path` - The path to the BAM file.
    ///
    /// # Examples
    ///
    /// ```
    /// use cnv_from_bam::BamFile;
    /// use std::path::PathBuf;
    ///
    /// let path = PathBuf::from("path/to/your.bam");
    /// let bam_file = BamFile::new(path);
    /// ```
    pub fn new(path: PathBuf) -> Self {
        BamFile { path }
    }
}

/// Calculate median value
fn median(numbers: &mut [u16]) -> Option<f64> {
    numbers.par_sort();
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        Some((numbers[mid - 1] as f64 + numbers[mid] as f64) / 2.0)
    } else {
        Some(numbers[mid] as f64)
    }
}

/// Calculate dynamically binned CNV for each contig
pub fn calculate_cnv(
    genome_length: usize,
    number_reads: usize,
    bins: FnvHashMap<String, Vec<u16>>,
) -> (FnvHashMap<String, Vec<f64>>, usize) {
    // bin_size = int(genome_length / (total_map_starts / reads_per_bin))
    let bin_width =
        (genome_length as f64 / (number_reads as f64 / READS_PER_BIN as f64).ceil()) as usize;
    let mut all_values: Vec<u16> = bins
        .values()
        .flat_map(|v| v.iter().copied())
        .collect::<Vec<_>>()
        .chunks(bin_width / BIN_SIZE)
        .map(|chunk| chunk.iter().sum::<u16>())
        .collect::<Vec<_>>();

    let median_value: f64 = median(&mut all_values).unwrap().round();
    println!("Median value: {}", median_value);

    let new_map: FnvHashMap<String, Vec<f64>> = bins
        .into_par_iter()
        .map(|(k, v)| {
            let sums = v
                .chunks(bin_width / BIN_SIZE)
                .map(|chunk| {
                    std::convert::Into::<f64>::into(chunk.iter().sum::<u16>())
                        / (median_value as f64)
                        * (PLOIDY as f64)
                })
                .collect::<Vec<f64>>();
            (k, sums)
        })
        .collect();
    (new_map, bin_width)
}

/// Results struct for python
#[pyclass]
#[derive(Debug)]
pub struct CnvResult {
    /// The CNV per contig
    #[pyo3(get)]
    pub cnv: FnvHashMap<String, Vec<f64>>,
    /// Bin width
    #[pyo3(get)]
    pub bin_width: usize,
}

/// Iterate a bam file
#[pyfunction]
fn iterate_bam_file(bam_file_path: PathBuf, threads: Option<usize>) -> PyResult<CnvResult> {
    let bam_file = BamFile::new(bam_file_path);
    let mut bam_reader = bam::IndexedReader::from_path(&bam_file.path).unwrap_or_else(|_| {
        panic!(
            "Error creating bam reader with file {}",
            bam_file.path.display()
        )
    });
    bam_reader
        .set_threads(threads.unwrap_or(num_cpus::get()))
        .expect("Error setting threads");
    let header = bam_reader.header().to_owned();
    let num_targets = header.target_count();

    // Get total number of reads

    let (genome_length, total_mapped, total_unmapped): (u64, u64, u64) =
        bam_reader.index_stats().unwrap().into_iter().fold(
            (0, 0, 0),
            |(acc2, acc3, acc4), (_, second, third, fourth)| {
                (acc2 + second, acc3 + third, acc4 + fourth)
            },
        );

    let mut results: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
    // Prepare all the CNV bins
    for target in 0..num_targets {
        let target_name = std::str::from_utf8(header.tid2name(target)).unwrap();
        let target_length: usize = header.target_len(target).unwrap().try_into().unwrap();
        results.insert(
            target_name.to_string(),
            vec![0; (target_length / BIN_SIZE) + 1],
        );
        println!("Added Target name: {}", target_name);
        println!("Target length: {}", target_length);
    }
    println!("Genome length: {}", genome_length);
    // number of reads in the bam file
    let mut number_reads: usize = 0;
    bam_reader.fetch(FetchDefinition::All).unwrap();
    let mut record = Record::new();
    let bar = ProgressBar::new(total_mapped + total_unmapped);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    while let Some(result) = bam_reader.read(&mut record) {
        match result {
            Ok(_) => {
                if !record.is_unmapped() {
                    let target_name =
                        std::str::from_utf8(header.tid2name(record.tid() as u32)).unwrap();

                    results.get_mut(target_name).unwrap()[record.pos() as usize / BIN_SIZE] += 1;
                }
                number_reads += 1;
            }
            Err(_) => panic!("BAM parsing failed..."),
        }
        bar.inc(1);
    }
    println!("Number of reads: {}", number_reads);
    let (cnv_profile, bin_width) = calculate_cnv(genome_length as usize, number_reads, results);
    let result = CnvResult {
        cnv: cnv_profile,
        bin_width,
    };
    Ok(result)
}

/// A Python module implemented in Rust.
#[pymodule]
fn cnv_from_bam(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(iterate_bam_file, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iterate_bam_file() {
        // Call the function with the temporary BAM file
        let bam_file_path = PathBuf::from("test/test.bam");
        let result = iterate_bam_file(bam_file_path.clone(), None).expect("Function failed");
        println!("{:#?}", result);
        // Assert that the function returns a CnvBins instance
        // assert_eq!(result.bins, vec![0, 0]);
    }
}
