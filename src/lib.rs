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
//! ```rust,ignore
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
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};

use noodles::bam::bai;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rayon::prelude::*;
use std::collections::HashSet;
use std::{error::Error, fs, fs::File, io, num::NonZeroUsize, path::PathBuf, thread};

use noodles::{
    bam, csi,
    sam::{self, record::MappingQuality},
};
use noodles_bgzf as bgzf;

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

/// Iterates over a BAM file, filters reads based on the mapping quality (`mapq_filter`),
/// and populates the provided frequency map with read counts for each genomic bin.
///
/// # Arguments
///
/// * `bam_file_path` - Path to the BAM file.
/// * `mapq_filter` - Minimum mapping quality to consider for read filtering.
/// * `genome_length` - A mutable reference to a usize to store the total length of the genome.
/// * `valid_number_reads` - A mutable reference to a usize to store the count of valid reads.
/// * `frequencies` - A mutable reference to a hash map to store read counts for each genomic bin.
///
/// # Returns
///
/// A result indicating successful execution or an error.
///
/// # Example
///
/// ```rust,ignore
/// use std::path::PathBuf;
/// use fnv::FnvHashMap;
///
/// let bam_file_path = PathBuf::from("path/to/sample.bam");
/// let mapq_filter = 30;
/// let mut genome_length = 0;
/// let mut valid_number_reads = 0;
/// let mut frequencies: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
///
/// let result = _iterate_bam_file(bam_file_path, mapq_filter, &mut genome_length, &mut valid_number_reads, &mut frequencies);
/// assert!(result.is_ok());
/// ```
///
fn _iterate_bam_file(
    bam_file_path: PathBuf,
    mapq_filter: u8,
    genome_length: &mut usize,
    valid_number_reads: &mut usize,
    frequencies: &mut FnvHashMap<String, Vec<u16>>,
    contigs: Option<&mut HashSet<String>>,
) -> Result<(), PyErr> {
    let bam_file = BamFile::new(bam_file_path);
    // Open bam file

    let file = File::open(&bam_file.path)?;

    let worker_count = thread::available_parallelism().unwrap_or(NonZeroUsize::MIN);
    let decoder = bgzf::MultithreadedReader::with_worker_count(worker_count, file);

    let mut bam_reader = bam::Reader::from(decoder);
    let header = bam_reader.read_header()?;
    let reference_sequences = header.reference_sequences();

    if contigs.as_ref().is_some_and(|c| c.is_empty()) {
        let contigs = contigs.unwrap();
        for (name, reference_sequence) in reference_sequences {
            contigs.insert(name.to_string());
            let len = (usize::from(reference_sequence.length()) / BIN_SIZE) + 1;
            let bins = vec![0; len];
            frequencies.insert(name.to_string(), bins);
        }
        // Calculate total genome length
        let _genome_length: usize = reference_sequences
            .iter()
            .map(|(_rs_name, rs)| usize::from(rs.length()))
            .sum();
        // Append total genome length
        *genome_length += _genome_length;
    } else if contigs.as_ref().is_some() {
        let contigs = contigs.unwrap();
        let contig_names = contigs.iter().cloned().collect::<HashSet<_>>();
        if contig_names != *contigs {
            return Err(PyRuntimeError::new_err(format!(
                "Reference contigs for bam file {:?} do not match previous bam file",
                bam_file.path
            )));
        }
    }
    // horrendous repetitive, dirty match to try a .csi index, then a .bai index, then no index
    let bar = match File::open(format!("{}.csi", bam_file.path.to_str().unwrap()))
        .map(csi::Reader::new)
    {
        Ok(mut index) => {
            let index = index.read_index().unwrap();
            let (total_mapped, total_unmapped): (u64, u64) = index
                .reference_sequences()
                .iter()
                .map(|x| x.metadata().unwrap())
                .fold(
                    (0, 0),
                    |(mapped_count_acc, unmapped_count_acc), metadata| {
                        (
                            mapped_count_acc + metadata.mapped_record_count(),
                            unmapped_count_acc + metadata.unmapped_record_count(),
                        )
                    },
                );
            let bar = ProgressBar::with_draw_target(
                Some(total_mapped + total_unmapped),
                ProgressDrawTarget::stdout(),
            )
            .with_message("BAM Records");
            bar.set_style(
                ProgressStyle::with_template(
                    "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
                )
                .unwrap()
                .progress_chars("##-"),
            );
            bar
        }
        Err(_) => {
            let bar = match bai::read(format!("{}.bai", bam_file.path.to_str().unwrap())) {
                Ok(index) => {
                    let (total_mapped, total_unmapped): (u64, u64) = index
                        .reference_sequences()
                        .iter()
                        .map(|x| x.metadata().unwrap())
                        .fold(
                            (0, 0),
                            |(mapped_count_acc, unmapped_count_acc), metadata| {
                                (
                                    mapped_count_acc + metadata.mapped_record_count(),
                                    unmapped_count_acc + metadata.unmapped_record_count(),
                                )
                            },
                        );
                    let bar = ProgressBar::with_draw_target(
                        Some(total_mapped + total_unmapped),
                        ProgressDrawTarget::stdout(),
                    )
                    .with_message("BAM Records");
                    bar.set_style(
                        ProgressStyle::with_template(
                            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}",
                        )
                        .unwrap()
                        .progress_chars("##-"),
                    );
                    bar
                }
                Err(_) => {
                    println!("No index found for bam file: {:?}", bam_file.path);
                    let bar = ProgressBar::with_draw_target(None, ProgressDrawTarget::stdout())
                        .with_message("BAM Records");
                    bar.set_style(
                        ProgressStyle::with_template(
                            "[{elapsed_precise}] {spinner} {pos:>7} {msg}",
                        )
                        .unwrap()
                        // For more spinners check out the cli-spinners project:
                        // https://github.com/sindresorhus/cli-spinners/blob/master/spinners.json
                        .tick_strings(&["⣾", "⣽", "⣻", "⢿", "⡿", "⣟", "⣯", "⣷"]),
                    );
                    bar
                }
            };
            bar
        }
    };

    // Get total number of reads, mapped and unmapped

    // number of reads in the bam file that match our criteria
    let mut _valid_number_reads: usize = 0;
    // Setup progress bar

    // Allocate record
    let mut record = bam::lazy::Record::default();

    while bam_reader.read_lazy_record(&mut record)? != 0 {
        if filter(&record, mapq_filter) {
            let reference_sequence_name =
                reference_sequence_name(&record, &header)?.expect("missing reference sequence ID");
            let start = record.alignment_start()?.expect("missing alignment start");

            let bins = frequencies
                .get_mut(reference_sequence_name)
                .expect("missing frequencies");
            let bin = (usize::from(start) - 1) / BIN_SIZE;
            bins[bin] += 1;
            _valid_number_reads += 1;
        }
        bar.inc(1)
    }
    bar.finish();
    println!(
        "Number of reads matching criteria mapq >{} and aligned: {} for bam file: {:?}",
        mapq_filter, _valid_number_reads, bam_file.path
    );
    // Accumulate the valid numebr of reads
    *valid_number_reads += _valid_number_reads;
    Ok(())
}

/// Iterates over a BAM file to compute Copy Number Variation (CNV) profiles for each contig.
///
/// This function reads through a BAM file, filters the alignments based on mapping quality,
/// and calculates CNV profiles using dynamically sized bins. It uses multithreading for improved performance.
///
/// # Arguments
///
/// * `bam_file_path` - A `PathBuf` specifying the path to the BAM file.
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
/// let result = iterate_bam_file(bam_path, Some(60));
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
    _threads: Option<u8>,
    mapq_filter: Option<u8>,
) -> PyResult<CnvResult> {
    let mapq_filter = mapq_filter.unwrap_or(60);
    let mut frequencies: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
    let genome_length = &mut 0_usize;
    let valid_number_reads = &mut 0_usize;
    let contigs: &mut HashSet<String> = &mut HashSet::new();
    if bam_file_path.is_dir() {
        println!("Process directory: {:?}", &bam_file_path);
        // Iterate over all files in the directory
        for entry in fs::read_dir(&bam_file_path)? {
            let entry = entry?;
            let path = entry.path();

            // Check if the entry is a file and ends with .bam
            if path.is_file() && path.extension().and_then(|s| s.to_str()) == Some("bam") {
                // Process the .bam file
                println!("Processing BAM file: {:?}", path);
                // Insert your BAM file processing logic here
                _iterate_bam_file(
                    path,
                    mapq_filter,
                    genome_length,
                    valid_number_reads,
                    &mut frequencies,
                    Some(contigs),
                )
                .unwrap();
            }
        }
    } else if bam_file_path.is_file()
        && bam_file_path.extension().and_then(|s| s.to_str()).unwrap() == "bam"
    {
        // The path is a single .bam file
        println!("Processing single BAM file: {:?}", &bam_file_path);
        // Insert your BAM file processing logic here
        _iterate_bam_file(
            bam_file_path,
            mapq_filter,
            genome_length,
            valid_number_reads,
            &mut frequencies,
            Some(contigs),
        )?;
    } else {
        // The path is neither a directory nor a .bam file
        println!("The path is neither a directory nor a .bam file.");
    }

    let (cnv_profile, bin_width) = calculate_cnv(*genome_length, *valid_number_reads, frequencies);
    let result = CnvResult {
        cnv: cnv_profile,
        bin_width,
        genome_length: *genome_length,
    };
    Ok(result)
}

/// Filters a BAM record based on mapping quality and flags.
fn filter(record: &bam::lazy::Record, mapping_quality: u8) -> bool {
    let min_mapping_quality: MappingQuality = match MappingQuality::new(mapping_quality) {
        Some(mapq) => mapq,
        None => unreachable!(),
    };

    let flags = record.flags();

    !flags.is_unmapped()
        && record
            .mapping_quality()
            .map(|mapq| mapq >= min_mapping_quality)
            .unwrap_or(true)
        && !flags.is_secondary()
        && !flags.is_supplementary()
}

/// Get the ref seq name
fn reference_sequence_name<'a>(
    record: &bam::lazy::Record,
    header: &'a sam::Header,
) -> io::Result<Option<&'a str>> {
    record
        .reference_sequence_id()?
        .map(|id| {
            header
                .reference_sequences()
                .get_index(id)
                .map(|(name, _)| name.as_ref())
                .ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidData, "invalid reference sequence ID")
                })
        })
        .transpose()
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
