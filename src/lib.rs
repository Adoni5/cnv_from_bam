#![deny(missing_docs, clippy::missing_docs_in_private_items)]
//! # CNV From BAM
//!
//! `cnv_from_bam` is a Rust library designed to efficiently calculate dynamic Copy Number Variation (CNV)
//! profiles from sequence alignments contained in BAM files. It is implemented as a Python module using PyO3,
//! making it easy to integrate into Python-based genomic analysis workflows.
//!
//! The library processes genomic data represented in BAM format to identify and quantify variations in copy
//! number across different regions of the genome. This functionality is especially relevant in fields like
//! cancer research, where CNVs are a key area of study.
//!
//! ## Features
//!
//! - **Efficient Processing**: Leverages the performance capabilities of Rust for handling large BAM files.
//! - **Python Integration**: Facilitates integration with Python, enabling use in Python-based genomic workflows.
//! - **Multithreading**: Improves performance through multithreaded processing, particularly advantageous for large datasets.
//! - **Dynamic Binning**: Adapts binning dynamically based on total read counts and genome length, suitable for a variety of genomic sizes and complexities.
//! - **CNV Computation**: Computes CNV values for each bin in each contig, aiding in the analysis of genomic copy number variations.
//!
//! ## Usage
//!
//! The library provides several functions, including:
//!
//! - `iterate_bam_file`: Iterates over a BAM file to filter alignments and compute CNV profiles with dynamic bin sizing.
//! - `calculate_cnv`: Calculates CNVs for each contig based on genome length, read counts, and read bins.
//! - `median`: Utility function to determine the median of a `u16` slice, used in CNV calculations.
//!
//! ### Example
//!
//! Using `iterate_bam_file` to process a BAM file:
//!
//! ```rust
//! use cnv_from_bam::iterate_bam_file;
//! use std::path::PathBuf;
//!
//! let bam_path = PathBuf::from("path/to/bam/file.bam");
//! let result = iterate_bam_file(bam_path, Some(4), Some(60));
//! // Handle result...
//! ```
//!
//! ## Installation
//!
//! Include this library as a dependency in your Rust project's `Cargo.toml` and build it to use in your Python scripts.
//!
//! ## Testing
//!
//! The module includes a suite of unit tests to ensure functionality reliability and correctness.
//!
//! ---
//!

use fnv::FnvHashMap;
use indicatif::{ProgressBar, ProgressStyle};
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

/// Calculates the median of a slice of `u16` numbers.
///
/// The function first sorts the given slice in place and then computes the median.
/// If the slice length is even, the median is the average of the two middle values.
/// If the slice length is odd, the median is the middle value.
///
/// # Arguments
///
/// * `numbers` - A mutable reference to a slice of `u16` numbers.
///
/// # Returns
///
/// An `Option<f64>` containing the median. If the slice is empty, it returns `None`.
///
/// # Examples
///
/// ```rust,ignore
/// # fn main() {
/// // Even number of elements
/// let mut even_numbers = [2, 3, 1, 4];
/// assert_eq!(median::median(&mut even_numbers), Some(2.5));
///
/// // Odd number of elements
/// let mut odd_numbers = [1, 2, 3, 4, 5];
/// assert_eq!(median::median(&mut odd_numbers), Some(3.0));
///
/// // Empty slice
/// let mut no_numbers: [u16; 0] = [];
/// assert_eq!(median::median(&mut no_numbers), None);
/// # }
/// ```
///
fn median(numbers: &mut [u16]) -> Option<f64> {
    if numbers.is_empty() {
        return None;
    }
    numbers.par_sort();
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        Some((numbers[mid - 1] as f64 + numbers[mid] as f64) / 2.0)
    } else {
        Some(numbers[mid] as f64)
    }
}

/// Calculates Copy Number Variations (CNVs) for each contig in a genome.
///
/// This function dynamically bins the genome based on the total number of reads and
/// the length of the genome. It then calculates the CNV value for each bin in each contig.
///
/// # Arguments
///
/// * `genome_length` - A `usize` representing the total length of the genome.
/// * `number_reads` - A `usize` representing the total number of reads.
/// * `bins` - A `FnvHashMap<String, Vec<u16>>` where keys are contig names and values
///            are vectors containing read counts for each bin in that contig.
///
/// # Returns
///
/// Returns a tuple containing:
/// - A `FnvHashMap<String, Vec<f64>>` where each key is a contig name and the value
///   is a vector of CNV values for each bin in that contig.
/// - A `usize` representing the width of each bin.
///
/// # Examples
///
/// ```rust
/// use fnv::FnvHashMap;
/// use cnv_from_bam::calculate_cnv;  // Replace with the actual crate name
///
/// // Example data
/// let genome_length = 3000;
/// let number_reads = 1500;
/// let mut bins = FnvHashMap::default();
/// bins.insert("contig1".to_string(), vec![100, 200, 150, 250]);
/// bins.insert("contig2".to_string(), vec![300, 400, 350, 450]);
///
/// let (cnv_map, bin_width) = calculate_cnv(genome_length, number_reads, bins);
///
/// // Use cnv_map and bin_width as needed
/// // e.g., print the CNV values for contig1
/// println!("Contig1 CNVs: {:?}", cnv_map.get("contig1"));
/// ```
///
pub fn calculate_cnv(
    genome_length: usize,
    number_reads: usize,
    bins: FnvHashMap<String, Vec<u16>>,
) -> (FnvHashMap<String, Vec<f64>>, usize) {
    // bin_size = int(genome_length / (total_map_starts / reads_per_bin))
    let bin_width =
        (genome_length as f64 / (number_reads as f64 / READS_PER_BIN as f64).ceil()) as usize;
    let chunk_size = bin_width / BIN_SIZE;
    let chunk_size = if chunk_size == 0 { 1 } else { chunk_size };
    let mut all_values: Vec<u16> = bins
        .values()
        .flat_map(|v| v.iter().copied())
        .collect::<Vec<_>>()
        .chunks(chunk_size)
        .map(|chunk| chunk.iter().sum::<u16>())
        .collect::<Vec<_>>();

    let median_value: f64 = median(&mut all_values).unwrap().round();
    println!("Median value: {}", median_value);

    let new_map: FnvHashMap<String, Vec<f64>> = bins
        .into_par_iter()
        .map(|(k, v)| {
            let sums = v
                .chunks(chunk_size)
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
    /// Genome length
    #[pyo3(get)]
    pub genome_length: usize,
}

/// Iterates over a BAM file to compute Copy Number Variation (CNV) profiles for each contig.
///
/// This function reads through a BAM file, filters the alignments based on mapping quality,
/// and calculates CNV profiles using dynamically sized bins. It uses multithreading for improved performance.
///
/// # Arguments
///
/// * `bam_file_path` - A `PathBuf` specifying the path to the BAM file.
/// * `threads` - An `Option<usize>` specifying the number of threads to use for reading the BAM file.
///              If `None`, the number of CPUs is used.
/// * `mapq_filter` - An `Option<u8>` specifying the minimum mapping quality to consider.
///                   Alignments with a mapping quality below this value are ignored. Default is 60.
///
/// # Returns
///
/// A `PyResult<CnvResult>` containing the CNV profiles and additional information about the analysis.
///
/// # Examples
///
/// ```rust,ignore
/// use std::path::PathBuf;
/// use cnv_from_bam::iterate_bam_file;  // Replace with the actual crate name
///
/// // Example usage of iterate_bam_file
/// let bam_path = PathBuf::from("path/to/bam/file.bam");
/// let result = iterate_bam_file(bam_path, Some(4), Some(60));
/// match result {
///     Ok(cnv_result) => {
///         // Process the CNV results
///         println!("Bin width: {}", cnv_result.bin_width);
///     },
///     Err(e) => println!("Error: {:?}", e),
/// }
/// ```
///
/// Note: This function is designed to be used in a Python context through PyO3 bindings.
#[pyfunction]
fn iterate_bam_file(
    bam_file_path: PathBuf,
    threads: Option<usize>,
    mapq_filter: Option<u8>,
) -> PyResult<CnvResult> {
    let mapq_filter = mapq_filter.unwrap_or(60);
    let bam_file = BamFile::new(bam_file_path);
    // Open bam file
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

    // Get total number of reads, mapped and unmapped
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
    // number of reads in the bam file that match our criteria
    let mut valid_number_reads: usize = 0;
    bam_reader.fetch(FetchDefinition::All).unwrap();
    // Setup progress bar
    let bar = ProgressBar::new(total_mapped + total_unmapped);
    bar.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
        )
        .unwrap()
        .progress_chars("##-"),
    );
    // Allocate record
    let mut record = Record::new();
    while let Some(result) = bam_reader.read(&mut record) {
        match result {
            Ok(_) => {
                if !record.is_unmapped()
                    & (record.mapq() >= mapq_filter)
                    & !record.is_secondary()
                    & !record.is_supplementary()
                {
                    let target_name =
                        std::str::from_utf8(header.tid2name(record.tid() as u32)).unwrap();

                    results.get_mut(target_name).unwrap()[record.pos() as usize / BIN_SIZE] += 1;
                    valid_number_reads += 1;
                }
            }
            Err(_) => panic!("BAM parsing failed..."),
        }
        bar.inc(1);
    }
    println!(
        "Number of reads matching criteria mapq >{} and aligned: {}",
        mapq_filter, valid_number_reads
    );
    let (cnv_profile, bin_width) =
        calculate_cnv(genome_length as usize, valid_number_reads, results);
    let result = CnvResult {
        cnv: cnv_profile,
        bin_width,
        genome_length: genome_length as usize,
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
        let result =
            iterate_bam_file(bam_file_path.clone(), None, Some(60)).expect("Function failed");
        println!("{:#?}", result);
        // Assert that the function returns a CnvBins instance
        // assert_eq!(result.bins, vec![0, 0]);
    }

    #[test]
    fn test_median_even_number_of_elements() {
        let mut numbers = [2, 3, 1, 4];
        let result = median(&mut numbers).unwrap();
        assert_eq!(result, 2.5);
    }

    #[test]
    fn test_median_odd_number_of_elements() {
        let mut numbers = [1, 2, 3, 4, 5];
        let result = median(&mut numbers).unwrap();
        assert_eq!(result, 3.0);
    }

    #[test]
    fn test_median_empty_slice() {
        let mut numbers: [u16; 0] = [];
        assert!(median(&mut numbers).is_none());
    }

    #[test]
    fn test_median_with_repeated_elements() {
        let mut numbers = [1, 1, 1, 1];
        let result = median(&mut numbers).unwrap();
        assert_eq!(result, 1.0);
    }

    #[test]
    fn test_median_with_single_element() {
        let mut numbers = [42];
        let result = median(&mut numbers).unwrap();
        assert_eq!(result, 42.0);
    }
}
