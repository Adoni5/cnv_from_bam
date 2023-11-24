#![deny(missing_docs, clippy::missing_docs_in_private_items)]
//! A python module for quickly and efficiently calculating a dynamic CNV profile from contigs continaed in a  BAM file.

use fnv::FnvHashMap;
use pyo3::prelude::*;
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
fn median(numbers: &mut [u16]) -> Option<f32> {
    numbers.sort_unstable();
    let mid = numbers.len() / 2;
    if numbers.len() % 2 == 0 {
        Some((numbers[mid - 1] as f32 + numbers[mid] as f32) / 2.0)
    } else {
        Some(numbers[mid] as f32)
    }
}

/// Calculate dynamically binned CNV for each contig
pub fn calculate_cnv(
    genome_length: usize,
    number_reads: usize,
    bins: FnvHashMap<String, Vec<u16>>,
) -> FnvHashMap<String, Vec<u16>> {
    // bin_size = int(genome_length / (total_map_starts / reads_per_bin))
    let bin_width =
        (genome_length as f64 / (number_reads as f64 / READS_PER_BIN as f64).ceil()) as usize;
    let mut all_values: Vec<u16> = bins.values().flat_map(|v| v.iter().copied()).collect();

    let median_value: usize = median(&mut all_values).unwrap().round() as usize;

    let new_map: FnvHashMap<String, Vec<u16>> = bins
        .into_iter()
        .map(|(k, v)| {
            let sums = v
                .chunks(bin_width)
                .map(|chunk| {
                    chunk
                        .iter()
                        .sum::<u16>()
                        .checked_div(median_value as u16)
                        .unwrap_or(0_u16)
                        .checked_mul(PLOIDY)
                        .unwrap_or(0)
                })
                .collect::<Vec<u16>>();
            (k, sums)
        })
        .collect();
    new_map
}

/// Results struct for python
#[pyclass]
#[derive(Debug)]
pub struct CnvResult {
    /// The CNV per contig
    #[pyo3(get)]
    pub cnv: FnvHashMap<String, Vec<u16>>,
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
    // calculate whole genome length
    let genome_length: u64 = (0..num_targets)
        .map(|tid| header.target_len(tid).unwrap())
        .sum();

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
    while let Some(result) = bam_reader.read(&mut record) {
        match result {
            Ok(_) => {
                println!("Read sequence: {:?}", record.pos());
                let target_name =
                    std::str::from_utf8(header.tid2name(record.tid() as u32)).unwrap();
                results.get_mut(target_name).unwrap()[record.pos() as usize / BIN_SIZE] += 1;
                number_reads += 1;
            }
            Err(_) => panic!("BAM parsing failed..."),
        }
    }
    println!("Number of reads: {}", number_reads);
    let cnv_profile = calculate_cnv(genome_length as usize, number_reads, results);
    let result = CnvResult {
        cnv: cnv_profile,
        bin_width: BIN_SIZE,
    };
    // let counts = bam_reader
    //     .rc_records()
    //     .map(|x| x.expect("Failure parsing Bam file")) // Converts Iterator<Result<Record, Error>> to Iterator<Record>, ignoring errors // Keep only first occurrences
    //     .inspect(|x| println!("{:#?}", x.flags()))
    //     .map(|record| record.pos() as usize / BIN_SIZE)
    //     .inspect(|x| println!("{:#?}", x)) // Map each record to its bin
    //     .fold(HashMap::new(), |mut acc, bin| {
    //         // Count occurrences in each bin
    //         *acc.entry(bin).or_insert(0) += 1;
    //         acc
    //     });
    // println!("{:#?}", counts);
    // println!("{:#?}", cnv_profile);
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
