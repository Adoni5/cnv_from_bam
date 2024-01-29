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

use log::LevelFilter;
use log::{error, info, warn};
use log::{Level, Metadata, Record};
use natord::compare;
use noodles::bam::bai;
use noodles::{
    bam, csi,
    sam::{self, record::MappingQuality},
};
use noodles_bgzf as bgzf;
use once_cell::sync::Lazy;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use rayon::prelude::*;
use std::collections::HashSet;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use std::{env, error::Error, fs, fs::File, io, num::NonZeroUsize, path::PathBuf, thread};

/// Number of reads expected in given bin width on the genome
const READS_PER_BIN: usize = 100;
/// Genome bin size for binning read starts
const BIN_SIZE: usize = 1000;
/// Excpected ploidy
const PLOIDY: u16 = 2;

/// Global static variable to indicate if the handler is set
static EXIT: Lazy<Arc<AtomicBool>> = Lazy::new(|| Arc::new(AtomicBool::new(false)));
/// Rust based logger
static LOGGER: SimpleLogger = SimpleLogger;
/// Dynamic result type for fancy unwrapping
pub type DynResult<T> = Result<T, Box<dyn Error + 'static>>;
/// Represents a BAM file with its path.
pub struct BamFile {
    /// Path to the BAM file, mainly used to convey intent.
    path: PathBuf,
}

/// Setup the log for the library
struct SimpleLogger;

impl log::Log for SimpleLogger {
    fn enabled(&self, metadata: &Metadata) -> bool {
        metadata.level() <= Level::Info
    }

    fn log(&self, record: &Record) {
        if self.enabled(record.metadata()) {
            println!(
                "{} ({}) - {} - {}",
                record.level(),
                record.target(),
                record.line().unwrap_or(0),
                record.args()
            );
        }
    }

    fn flush(&self) {}
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

/// Calculates the variance of a dataset.
///
/// The variance is computed as the average of the squared differences from the Mean.
/// Returns `None` if the dataset is empty.
///
/// # Arguments
///
/// * `data` - A slice of usize values representing the dataset.
///
/// # Returns
///
/// An `Option<f64>` representing the variance of the dataset.
///
/// # Examples
///
/// Basic usage:
///
/// ```rust
/// # use cnv_from_bam::calculate_variance;
/// let data = vec![1, 2, 3, 4, 5];
/// let variance = calculate_variance(&data).unwrap();
/// assert_eq!(variance, 2.0);
/// ```
///
/// When the dataset is empty:
///
/// ```
/// # use cnv_from_bam::calculate_variance;
/// let empty: Vec<usize> = vec![];
/// assert!(calculate_variance(&empty).is_none());
/// ```
pub fn calculate_variance<'a>(data: impl Iterator<Item = &'a f64>) -> Option<f64> {
    let (count, sum, sum_sq) = data.fold((0, 0.0, 0.0), |(count, sum, sum_sq), &value| {
        (count + 1, sum + value, sum_sq + value * value)
    });

    if count == 0 {
        None
    } else {
        let mean = sum / count as f64;
        let variance = (sum_sq / count as f64) - (mean * mean);
        Some(variance)
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
    // info!("{bins.keys():#?}");
    let median_value: f64 = median(&mut all_values).unwrap().round();
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
    pub cnv: PyObject,
    /// Bin width
    #[pyo3(get)]
    pub bin_width: usize,
    /// Genome length
    #[pyo3(get)]
    pub genome_length: usize,
    /// Variance of the whole genome
    #[pyo3(get)]
    pub variance: f64,
}

/// Get the total number of sequences from a BAM index, returning `Some(count)` if the operation
/// is successful, or `None` if there is an issue opening any of the sequences.
///
/// # Examples
///
/// ```rust,ignore
/// use noodles::csi;
///
/// // Create a mock BAM index for testing
/// let index = csi::Index::default();
///
/// // Get the total sequences
/// let result = total_sequences(index);
///
/// // Since the mock index has no actual data, the result should be `Some(0)`
/// assert_eq!(result, Some(0));
/// ```
///
/// ```rust,ignore
/// use noodles::csi;
///
/// // Create a mock BAM index with some metadata for testing
/// let mut index = csi::Index::default();
/// index.reference_sequences_mut().push(csi::ReferenceSequence::new(0, 10, 5, 5, 0));
///
/// // Get the total sequences
/// let result = total_sequences(index);
///
/// // The total sequences in the mock index is 5, so the result should be `Some(5)`
/// assert_eq!(result, Some(5));
/// ```
fn total_sequences(index: csi::Index) -> Option<u64> {
    let (mut total_mapped, total_unmapped): (u64, u64) = index
        .reference_sequences()
        .iter()
        .map(|x| x.metadata())
        .fold(
            (0, 0),
            |(mapped_count_acc, unmapped_count_acc), metadata_option| {
                metadata_option
                    .map(|metadata| {
                        (
                            mapped_count_acc + metadata.mapped_record_count(),
                            unmapped_count_acc + metadata.unmapped_record_count(),
                        )
                    })
                    .unwrap_or((mapped_count_acc, unmapped_count_acc))
            },
        );
    total_mapped += total_unmapped;
    if total_mapped == 0 {
        None
    } else {
        Some(total_mapped)
    }
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
    exclude_supplementary: Option<bool>,
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
    // Get total number of reads, mapped and unmapped
    // horrendous repetitive, dirty match to try a .csi index, then a .bai index, then no index
    let bar = match File::open(format!("{}.csi", bam_file.path.to_str().unwrap()))
        .map(csi::Reader::new)
    {
        Ok(mut index) => {
            let index = index.read_index().unwrap();

            match total_sequences(index) {
                Some(total_sequences) => {
                    let bar = ProgressBar::with_draw_target(
                        Some(total_sequences),
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
                None => {
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
            }
        }
        Err(_) => {
            match bai::read(format!("{}.bai", bam_file.path.to_str().unwrap())) {
                Ok(index) => {
                    match total_sequences(index) {
                        Some(total_sequences) => {
                            let bar = ProgressBar::with_draw_target(
                                Some(total_sequences),
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
                        None => {
                            let bar =
                                ProgressBar::with_draw_target(None, ProgressDrawTarget::stdout())
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
                    }
                }
                Err(_) => {
                    warn!("No index found for bam file: {:?}", bam_file.path);
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
            }
        }
    };
    if env::var("CI").is_ok() {
        bar.set_draw_target(ProgressDrawTarget::hidden())
    }

    // number of reads in the bam file that match our criteria
    let mut _valid_number_reads: usize = 0;

    // Allocate record
    let mut record = bam::lazy::Record::default();

    while (bam_reader.read_lazy_record(&mut record)? != 0) & !EXIT.load(Ordering::SeqCst) {
        if filter(&record, mapq_filter, exclude_supplementary.unwrap_or(false)) {
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
    if EXIT.load(Ordering::SeqCst) {
        return Err(PyErr::new::<pyo3::exceptions::PySystemExit, _>(
            "SIGINT or SIGTERM intercepted, Terminating function",
        ));
    }
    bar.finish();
    info!(
        "Number of reads matching criteria mapq >{} and aligned: {} for bam file: {:?}",
        mapq_filter, _valid_number_reads, bam_file.path
    );
    // Accumulate the valid numebr of reads
    *valid_number_reads += _valid_number_reads;
    Ok(())
}

/// Merge two vectors, summing the values at each index.
///
/// The function takes two slices of u16 values (`vec1` and `vec2`) and returns a new vector
/// where the values at each index are the sum of the corresponding values from both input vectors.
/// If one vector is shorter than the other, the missing elements are assumed to be 0.
///
/// # Arguments
///
/// * `vec1` - A slice of u16 values representing the first vector.
/// * `vec2` - A slice of u16 values representing the second vector.
///
/// # Returns
///
/// A new vector containing the sum of values at each index from both input vectors.
///
/// # Examples
///
/// ```rust,ignore
/// let vec1 = vec![1, 2, 3];
/// let vec2 = vec![4, 5, 6];
///
/// let result = merge_vecs(&vec1, &vec2);
/// assert_eq!(result, vec![5, 7, 9]);
/// ```
///
/// ```
/// let vec1 = vec![1, 2, 3];
/// let vec2 = vec![4, 5, 6, 7];
///
/// let result = merge_vecs(&vec1, &vec2);
/// assert_eq!(result, vec![5, 7, 9, 7]);
/// ```
fn merge_vecs(vec1: &[u16], vec2: &[u16]) -> Vec<u16> {
    // Determine the length of the merged vector
    let merged_len = vec1.len().max(vec2.len());

    // Use zip to iterate over both vectors simultaneously, summing the values at each index
    let merged_vec: Vec<u16> = vec1
        .iter()
        .chain(std::iter::repeat(&0).take(merged_len - vec1.len()))
        .zip(vec2.iter())
        .map(|(a, b)| a + b)
        .collect();

    merged_vec
}

/// Merge a pydict with a hashmap, and then reset the pydict to have the updated values
fn merge_and_reset(
    py_dict: &PyDict,
    rust_map: FnvHashMap<String, Vec<u16>>,
) -> PyResult<FnvHashMap<String, Vec<u16>>> {
    // Acquire the GIL
    Python::with_gil(|py| -> PyResult<FnvHashMap<String, Vec<u16>>> {
        // Convert PyDict to HashMap
        let py_dict_map: FnvHashMap<String, Vec<u16>> = py_dict
            .iter()
            .map(|(key, value)| {
                let key_str = key.to_string();
                let value_list: Vec<u16> = value.extract().unwrap_or_default();
                (key_str, value_list)
            })
            .collect();

        // Merge values from rust_map into py_dict_map
        let merged_map: FnvHashMap<String, Vec<u16>> = rust_map
            .into_iter()
            .map(|(rust_key, rust_value)| {
                let merged_values: Vec<u16> =
                    merge_vecs(py_dict_map.get(&rust_key).unwrap_or(&vec![]), &rust_value);
                (rust_key, merged_values)
            })
            .collect();

        // Reset values of the Python dictionary to the updated HashMap values
        py_dict.clear();
        for (key, values) in merged_map.clone() {
            let py_key = key.to_object(py);
            let py_values = values.to_object(py);
            py_dict.set_item(py_key, py_values)?;
        }
        Ok(merged_map)
    })
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
/// * `_threads` - An `Option<u8>` specifying the number of threads to use. Default is the available parallelism - a sensible number.
/// * `log_level` - An `Option<u8>` specifying the log level. Default is 20 (INFO).
/// * `exclude_supplementary` - An `Option<bool>` specifying whether to filter out supplementary alignments. Default is true.
/// * `copy_numbers` - A python dictionary (Optional) - this will be mutated in place with the intermediate binned
///                    Read starts per contig, so we can iteratively update them.
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
#[pyo3(signature = (bam_file_path, _threads=None, mapq_filter=60, log_level=None, exclude_supplementary=true, copy_numbers=None), text_signature = None,)]
fn iterate_bam_file(
    bam_file_path: PathBuf,
    _threads: Option<u8>,
    mapq_filter: Option<u8>,
    log_level: Option<u8>,
    exclude_supplementary: Option<bool>,
    copy_numbers: Option<&PyDict>,
) -> PyResult<CnvResult> {
    let mapq_filter = mapq_filter.unwrap_or(60);
    let mut frequencies: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
    let genome_length = &mut 0_usize;
    let valid_number_reads = &mut 0_usize;
    let contigs: &mut HashSet<String> = &mut HashSet::new();
    let level_filter = match log_level {
        Some(level) => match level {
            0 => LevelFilter::Off,
            40 => LevelFilter::Error,
            30 => LevelFilter::Warn,
            20 => LevelFilter::Info,
            10 => LevelFilter::Debug,
            5 => LevelFilter::Trace,
            _ => LevelFilter::Info,
        },
        None => LevelFilter::Info,
    };
    if let Ok(()) = log::set_logger(&LOGGER) {
        log::set_max_level(level_filter)
    }

    if bam_file_path.is_dir() {
        if copy_numbers.is_some() {
            error!("Please refrain from passing `copy_numbers` if `bam_path` is A DIRECTORY.");
            return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                "Please refrain from passing `copy_numbers` if `bam_path` is A DIRECTORY.",
            ));
        }
        info!("Processing directory: {:?}", &bam_file_path);
        // Iterate over all files in the directory
        for entry in fs::read_dir(&bam_file_path)? {
            let entry = entry?;
            let path = entry.path();

            // Check if the entry is a file and ends with .bam
            if path.is_file() && path.extension().and_then(|s| s.to_str()) == Some("bam") {
                // Process the .bam file
                info!("Processing BAM file: {:?}", path);
                // Insert your BAM file processing logic here
                _iterate_bam_file(
                    path,
                    mapq_filter,
                    genome_length,
                    valid_number_reads,
                    &mut frequencies,
                    Some(contigs),
                    exclude_supplementary,
                )
                .unwrap();
            }
        }
    } else if bam_file_path.is_file()
        && bam_file_path.extension().and_then(|s| s.to_str()).unwrap() == "bam"
    {
        // The path is a single .bam file
        info!("Processing single BAM file: {:?}", &bam_file_path);
        // Insert your BAM file processing logic here
        _iterate_bam_file(
            bam_file_path,
            mapq_filter,
            genome_length,
            valid_number_reads,
            &mut frequencies,
            Some(contigs),
            exclude_supplementary,
        )?;
    } else {
        // The path is neither a directory nor a .bam file
        error!("The path is neither a directory nor a .bam file.");
        return Err(PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
            "Please refrain from passing `copy_numbers` if `bam_path` is A DIRECTORY.",
        ));
    }
    // If we have been given a data structure to use for the frequencies, update it here, and use the
    // combined values to calculate CNV
    let frequencies = if let Some(copy_numbers) = copy_numbers {
        let frequencies = merge_and_reset(copy_numbers, frequencies).unwrap();
        // # Update the number of reads that we have seen, to include the new and old counts
        *valid_number_reads = frequencies
            .values()
            .map(|x| x.par_iter().map(|x| *x as usize).sum::<usize>())
            .sum();
        frequencies
    } else {
        frequencies
    };

    // println!("{:#?}", frequencies["NC_000001.11"].iter().sum::<u16>());
    let (cnv_profile, bin_width) = calculate_cnv(*genome_length, *valid_number_reads, frequencies);
    let variance = cnv_profile.values().flatten();
    let variance = calculate_variance(variance).unwrap_or(0.0);
    // let variance = cnv_profile.values().flatten();
    // let variance = calculate_variance(variance).unwrap_or(0.0);
    let result = Python::with_gil(|py| -> PyResult<CnvResult> {
        let mut sorted_keys: Vec<_> = cnv_profile.keys().collect();
        sorted_keys.sort_by(|a, b| compare(a, b));
        let py_dict = PyDict::new(py);
        for key in sorted_keys {
            let value = cnv_profile.get(key).unwrap();
            py_dict.set_item(key, value)?;
        }

        let result = CnvResult {
            cnv: py_dict.into(),
            bin_width,
            genome_length: *genome_length,
            variance,
        };
        Ok(result)
    });
    result
}

/// Filters a BAM record based on mapping quality and flags.
fn filter(record: &bam::lazy::Record, mapping_quality: u8, exclude_supplementary: bool) -> bool {
    let min_mapping_quality: MappingQuality = match MappingQuality::new(mapping_quality) {
        Some(mapq) => mapq,
        None => unreachable!(),
    };

    let flags = record.flags();

    let filter = !flags.is_unmapped()
        && record
            .mapping_quality()
            .map(|mapq| mapq >= min_mapping_quality)
            .unwrap_or(true)
        && !flags.is_secondary();
    if exclude_supplementary {
        filter && !flags.is_supplementary()
    } else {
        filter
    }
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
    let r = EXIT.clone();
    ctrlc::set_handler(move || {
        r.store(true, Ordering::SeqCst);
        println!("{}[2J", 27 as char);
    })
    .expect("Error setting Ctrl-C handler");
    m.add_function(wrap_pyfunction!(iterate_bam_file, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_iterate_bam_file() {
        // Call the function with the temporary BAM file
        let bam_file_path = PathBuf::from("test/NA12878_4000_test.bam");
        let genome_length = &mut 0_usize;
        let valid_number_reads = &mut 0_usize;
        let mut frequencies: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
        let contigs: &mut HashSet<String> = &mut HashSet::new();
        _iterate_bam_file(
            bam_file_path,
            60,
            genome_length,
            valid_number_reads,
            &mut frequencies,
            Some(contigs),
            Some(false),
        )
        .unwrap();
        assert_eq!(valid_number_reads, &mut 3702_usize);
        // info!("{:#?}", result);
        // Assert that the function returns a CnvBins instance
        // assert_eq!(result.bins, vec![0, 0]);
    }

    #[test]
    fn test_iterate_bam_file_no_supp() {
        // Call the function with the temporary BAM file
        let bam_file_path = PathBuf::from("test/NA12878_4000_test.bam");
        let genome_length = &mut 0_usize;
        let valid_number_reads = &mut 0_usize;
        let mut frequencies: FnvHashMap<String, Vec<u16>> = FnvHashMap::default();
        let contigs: &mut HashSet<String> = &mut HashSet::new();
        _iterate_bam_file(
            bam_file_path,
            60,
            genome_length,
            valid_number_reads,
            &mut frequencies,
            Some(contigs),
            Some(true),
        )
        .unwrap();
        assert_eq!(valid_number_reads, &mut 3561_usize);
        // info!("{:#?}", result);
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

    #[test]
    fn test_merge_vecs() {
        // Test with equal-length vectors
        let vec1 = vec![1, 2, 3];
        let vec2 = vec![4, 5, 6];
        let result = merge_vecs(&vec1, &vec2);
        assert_eq!(result, vec![5, 7, 9]);

        // Test with the first vector shorter than the second one
        let vec1 = vec![1, 2, 3];
        let vec2 = vec![4, 5, 6, 7];
        let result = merge_vecs(&vec1, &vec2);
        assert_eq!(result, vec![5, 7, 9, 7]);

        // Test with an empty first vector
        let vec1: Vec<u16> = Vec::new();
        let vec2 = vec![4, 5, 6];
        let result = merge_vecs(&vec1, &vec2);
        assert_eq!(result, vec![4, 5, 6]);
    }

    #[test]
    fn test_total_sequences_empty_index() {
        // Create a mock BAM index for testing
        let index = csi::Index::default();

        // Get the total sequences
        let result = total_sequences(index);

        // Since the mock index has no actual data, the result should be `Some(0)`
        assert_eq!(result, None);
    }

    #[test]
    fn test_total_sequences_non_empty_index() {
        use noodles::csi::{
            self as csi,
            index::{reference_sequence, ReferenceSequence},
        };
        // Create a mock BAM index with some metadata for testing
        let index = csi::Index::builder();
        let metadata = reference_sequence::Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            55,
        );
        let reference_sequences = vec![ReferenceSequence::new(
            Default::default(),
            Default::default(),
            Some(metadata),
        )];
        let index = index.set_reference_sequences(reference_sequences);
        let index = index.build();
        // Get the total sequences
        let result = total_sequences(index);

        // The total sequences in the mock index is 5, so the result should be `Some(5)`
        assert_eq!(result, Some(110));
    }

    #[test]
    fn test_total_sequences_non_empty_index_no_metadata() {
        use noodles::csi::{self as csi, index::ReferenceSequence};
        // Create a mock BAM index with some metadata for testing
        let index = csi::Index::builder();
        let reference_sequences = vec![ReferenceSequence::new(
            Default::default(),
            Default::default(),
            None,
        )];
        let index = index.set_reference_sequences(reference_sequences);
        let index = index.build();
        // Get the total sequences
        let result = total_sequences(index);

        // The total sequences in the mock index is 5, so the result should be `Some(5)`
        assert_eq!(result, None);
    }

    #[test]
    fn test_total_sequences_mixed_metadata() {
        use noodles::csi::{
            self as csi,
            index::{reference_sequence, ReferenceSequence},
        };
        // Create a mock BAM index with some metadata for testing
        let index = csi::Index::builder();
        let metadata = reference_sequence::Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            55,
            55,
        );
        let metadata2 = reference_sequence::Metadata::new(
            bgzf::VirtualPosition::from(610),
            bgzf::VirtualPosition::from(1597),
            7,
            14,
        );
        let reference_sequences = vec![
            ReferenceSequence::new(Default::default(), Default::default(), Some(metadata)),
            ReferenceSequence::new(Default::default(), Default::default(), None),
            ReferenceSequence::new(Default::default(), Default::default(), Some(metadata2)),
        ];
        let index = index.set_reference_sequences(reference_sequences);
        let index = index.build();
        // Get the total sequences
        let result = total_sequences(index);

        // The total sequences in the mock index is 5, so the result should be `Some(5)`
        assert_eq!(result, Some(131));
    }
}
