# CNV From BAM

`cnv_from_bam` is a Rust library developed to efficiently calculate dynamic Copy Number Variation (CNV) profiles from sequence alignments contained in BAM files. It seamlessly integrates with Python using PyO3, making it an excellent choice for bioinformatics workflows involving genomic data analysis.

## Features

- **Efficient Processing**: Optimized for handling large genomic datasets in BAM format.
- **Python Integration**: Built with PyO3 for easy integration into Python-based genomic analysis workflows.
- **Multithreading Support**: Utilizes Rust's powerful concurrency model for improved performance.
- **Dynamic Binning**: Bins the genome dynamically based on total read counts and genome length.
- **CNV Calculation**: Accurately calculates CNV values for each bin across different contigs.
- **Directory Support**: Supports processing of multiple BAM files in a directory. (Requires alignment to the same reference in all BAM files)

## Installation

To use `cnv_from_bam` in your Rust project, add the following to your `Cargo.toml` file:

```toml
[dependencies]
cnv_from_bam = "0.1.0"  # Replace with the latest version
```

## Calculation explainer

`cnv_from_bam` calculates copy number from binned mapping starts.
Mapping starts from across the genome are binned into 1000 base wide windows.

The bin width of genomic bases for calculating the Copy Number for is calculated as follows -

```rust
let bin_width = genome_length as f64 / (number_mapped_reads as f64 / 100_f64).ceil()
let bin_width =
    (bin_width + (cnv_params.bin_size / 2)) / cnv_params.bin_size * cnv_params.bin_size;
```
Which effectively bins the reads so that we have roughly 100 reads in each bin. This assumes random sampling of reads across the genome.
The bin width is then binned to the nearest bin size for the mapped read starts, which by default is 1000 bases.

We then sum the number of binned reads in this "super" bin, so if the binwidth is 10000, we sum the values of the 10 1000 base mapped start counts in this larger bin.
We calculate the median for this new set of values.
We then finally calculate the CNV across this set of bins as so:

```rust
let new_map: FnvHashMap<String, Vec<f64>> = bins
    .into_par_iter()
    .map(|(k, v)| {
        let sums = v
            .chunks(chunk_size)
            .map(|chunk| {
                let x = std::convert::Into::<f64>::into(chunk.iter().sum::<u16>())
                    / median_value as f64
                    * (cnv_params.ploidy as f64);
            })
            .collect::<Vec<f64>>();
        (k, sums)
    })
    .collect();
```

We sum the binned read counts (`bins`) int. the chunks,
Chunk size is the `bin_width` from above divided by 1000, so it is scaled to the binned read counts from before.
We then divide this by the median value of read counts in the bins across the genome, and multiply that by the expected ploidy.
This gives us the Copy Number for each bin on the chromosome/contig.


## Usage

Here's a quick example of how to use the `iterate_bam_file` function:

```rust
use cnv_from_bam::iterate_bam_file;
use std::path::PathBuf;

let bam_path = PathBuf::from("path/to/bam/file.bam");
// Iterate over the BAM file and calculate CNV values for each bin. Number of threads is set to 4 and mapping quality filter is set to 60.
// If number of threads is not specified, it defaults to the number of logical cores on the machine.
let result = iterate_bam_file(bam_path, Some(4), Some(60), None, None);
// Process the result...
```

The results in this case are returned as a CnvResult, which has the following structure:

/// Results struct for python
```
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
```

Where `result.cnv` is a Python dict `PyObject` containing the Copy Number for each bin of `bin_width` bases for each contig in the reference genome, `result.bin_width` is the width of the bins in bases, `result.genome_length` is the total length of the genome and `result.variance` is a measure of the variance across the whole genome.

Variance is calculated as the average of the squared differences from the Mean.

> [!NOTE]
> **Note**: Only the main primary mapping alignment start is binned, Supplementary and Secondary alignments are ignored. Supplementary alignments can be included by setting `exclude_supplementary`

**Directory analysis**
To analyse a directory of BAM files, use the `iterate_bam_dir` function:

```rust
use cnv_from_bam::iterate_bam_dir;
use std::path::PathBuf;
let bam_path = PathBuf::from("path/to/bam_directory/");
// Iterate over the BAM files in the directory and calculate CNV values for the whole. Number of threads is set to 4 and mapping quality filter is set to 60.
// If number of threads is not specified, it defaults to the number of logical cores on the machine.
let result = iterate_bam_file(bam_path, Some(4), Some(60));
```

This again returns a CnvResult, but this time the CNV values are summed across all BAM files in the directory. The bin width and genome length are calculated based on the first BAM file in the directory.

> [!NOTE]
> **Note**: All BAM files in the directory must be aligned to the same reference genome.

## Python Integration

`cnv_from_bam` can be used in Python using the PyO3 bindings. To install the Python bindings, run:

`pip install cnv_from_bam`

The same `iterate_bam_file`  is availabl
e in python, accepting a path to a BAM file or a directory of BAM files, the number of threads (set to `None` to use the optimal number of threads for the machine), and the mapping quality filter.

Example simple plot in python, you will need `matplotlib` an `numpy` installed (`pip install matplotlib numpy`)
```python
from matplotlib import pyplot as plt
import matplotlib as mpl
from pathlib import Path
from cnv_from_bam import iterate_bam_file
import numpy as np
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 3))
total = 0
bam_path = Path("path/to/bam/file.bam");
# Iterate over the BAM file and calculate CNV values for each bin. Number of threads is set to 4 and mapping quality filter is set to 60.
# If number of threads is not specified, it defaults to the optimal number of threads for the machine.
result = iterate_bam_file(bam_path, _threads=4, mapq_filter=60);
for contig, cnv in result.cnv.items():
    ax.scatter(x=np.arange(len(cnv)) + total, y=cnv, s =0.1)
    total += len(cnv)

ax.set_ylim((0,8))
ax.set_xlim((0, total))
fig.savefig("Example_cnv_plot.png")
```
Should look something like this. Obviously the cnv data is just a dictionary of lists, so you can do whatever you want with it vis a vis matplotlib, seaborn, etc.
![example cnv plot](https://github.com/Adoni5/cnv_from_bam/blob/10a2b00a8832b46cacbff0e2f775a4f440844da0/example_cnv.png?raw=true)

### Python keyword arguments
It is possible to fix some of the dynamically calculated values, for comparison purposes or whatever strikes your fancy.

This parameter accepts a dictionary of keyword arguments that are passed through to CNV calculation functions. These are primarily used to customize the CNV profile calculation according to specific needs. Here are the keys that can be included in `py_kwargs`:

- `bin_size`: An override for the size of the bins which the BAM mapping starts are binned into. The default for this value is 1000.
- `bin_width`: The width of the bins used to segment the genome for CNV calculations. If not given, this is dynamically calculated as described in [Calculation explainer](#calculation-explainer).
- `genome_length`: The total length of the genome, which is used to calculate the number of bins across the genome. If not listed, this is calculated from the BAM file header by summing the lenght of each contig.
- `ploidy`: The expected number of copies of each chromosome in a normal, non-variant condition. This is used in normalization and analysis of CNV data. The default for this is 2.
- `target_reads`: The target number of reads per bin, which can be used to adjust the resolution of the CNV analysis based on sequencing depth. This is 100 by default.

Each key is optional and has default values or a value is calculated if not specified. When these arguments are provided, they override the default or calculated parameters used in CNV calculations.


### Iterative use
It is possible to iteratively add bam files to a continuing count. By passing a dictionary to iterate_bam_file, the intermediate mapping start counts are kept in this dictionary.

This is limited to parsing files one at a time, rather than by directory.


```python
from cnv_from_bam import iterate_bam_file
update = {}
bam_path = Path("path/to/bam/file.bam");
result = iterate_bam_file(bam_path, _threads=4, mapq_filter=60, copy_numbers=update);

bam_path_2 = Path("path/to/bam/file2.bam");
result = iterate_bam_file(bam_path_2, _threads=4, mapq_filter=60, copy_numbers=update);
# Result now contains the copy number as inferred by both BAMS
```

## Terminal/stdout Output
This is new in version >= 0.3. If you just want raw stdout from rust and no faffing with loggers, use v0.2.
### Progress Bar
By default, a progress bar is displayed, showing the progress of the iteration of each BAM file. To disable the progress bar, set the `CI` environment variable to `1` in your python script:

```python
import os
os.environ["CI"] = "1"
```

### Logging
We use the `log` crate for logging. By default, the log level is set to `INFO`, which means that the program will output the progress of the iteration of each BAM file. To disable all but warning and error logging, set the log level to `WARN` on the `iterate_bam_file` function:


```python


import logging
from cnv_from_bam import iterate_bam_file
iterate_bam_file(bam_path, _threads=4, mapq_filter=60, log_level=int(logging.WARN))

```

`getLevelName` is a function from the `logging` module that converts the log level to the integer value of the level. These values are

| Level | Value |
|-------|-------|
| CRITICAL | 50 |
| ERROR | 40 |
| WARNING | 30 |
| INFO | 20 |
| DEBUG | 10 |
| NOTSET | 0 |

> [!NOTE]
> In v0.3 a regression was introduced, whereby keeping the GIL for logging meant that BAM reading was suddenly single threaded again. Whilst it was possible to fix this and keep `PyO3-log`, I decided to go for truly maximum speed instead. The only drawback to the removal of `PyO3-log` in (v0.4+) is that log messages will not be handled by python loggers, so they won't be written out by a file handler, for example.


## Documentation

To generate the documentation, run:

```bash
cargo doc --open
```

## Contributing

Contributions to `cnv_from_bam` are welcome!

We use pre-commit hooks (particularly `cargo-fmt` and `ruff`) to ensure that code is formatted correctly and passes all tests before being committed. To install the pre-commit hooks, run:

```bash
git clone https://github.com/Adoni5/cnv_from_bam.git
cd cnv_from_bam
pip install -e .[dev]
pre-commit install -t pre-commit -t post-checkout -t post-merge
pre-commit run --all-files
```
## Changelog
### [v0.4.3](https://github.com/Adoni5/cnv_from_bam/releases/tag/v0.4.3)
### Iterative use
* It is possible to iteratively add bam files to a continuing count. By passing a dictionary to iterate_bam_file, the intermediate mapping start counts are kept in this dictionary.
This is limited to parsing files one at a time, rather than by directory.
See example above under [iterative use.](#iterative-use)

* Catches bug where metadata (mapped and unmapped count) for a reference sequence in a BAI or CSI file would return None, and crash the calculation. As this
was used to calculate the Progress bar total, just skips the offending reference sequence, returning - for both counts. May mean the progress bar can have a lower total than number of reads, but won't matter to final numbers.

* Adds Cargo tests to CI


### v0.4.2
* Returns the contig names naturally sorted, rather than in random order!! For example `chr1, chr2, chr3...chr22,chrM,chrX,chrY`!
Huge, will prevent some people getting repeatedly confused about expected CNV vs. Visualised and wasting an hour debugging a non existing issue.
* Returns variance across the whole genome in the CNV result struct.

### v0.4.1
* Add `exclude_supplementary` parameter to `iterate_bam_file`, to exclude supplementary alignments (default True)

### v0.4.0
* Remove `PyO3-log` for maximum speed. This means that log messages will not be handled by python loggers. Can set log level on call to `iterate_bam_file`

### [v0.3.0](https://github.com/Adoni5/cnv_from_bam/releases/tag/v0.0.3)
* Introduce `PyO3-log` for logging. This means that log messages can be handled by python loggers, so they can be written out by a file handler, for example.
* **HAS A LARGE PERFORMANCE ISSUE**
* Can disable progress bar display by setting `CI` environment variable to `1` in python script.

### [v0.2.0](https://github.com/Adoni5/cnv_from_bam/releases/tag/v0.0.2)
* Purely rust based BAM parsing, using noodles.
* Uses a much more sensible number for threading if not provided.
* Allows iteration of BAMS in a directory

### [v0.1.0](https://github.com/Adoni5/cnv_from_bam/releases/tag/v0.0.1)
* Initial release
* Uses `rust-bio/rust-htslib` for BAM parsing. Has to bind C code, is a faff.


## License

This project is licensed under the [Mozilla Public License 2.0](LICENSE).
