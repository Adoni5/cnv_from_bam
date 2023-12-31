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
// Iterate over the BAM files in teh directory and calculate CNV values for the whole. Number of threads is set to 4 and mapping quality filter is set to 60.
// If number of threads is not specified, it defaults to the number of logical cores on the machine.
let result = iterate_bam_file(bam_path, Some(4), Some(60));
```

This again returns a CnvResult, but this time the CNV values are summed across all BAM files in the directory. The bin width and genome length are calculated based on the first BAM file in the directory.

> [!NOTE]
> **Note**: All BAM files in the directory must be aligned to the same reference genome.

## Python Integration

`cnv_from_bam` can be used in Python using the PyO3 bindings. To install the Python bindings, run:

`pip install cnv_from_bam`

The same `iterate_bam_file`  is available in python, accepting a path to a BAM file or a directory of BAM files, the number of threads (set to `None` to use the optimal number of threads for the machine), and the mapping quality filter.
```python

Example simple plot in python
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
```
Should look something like this. Obviously the cnv data is just a dictionary of lists, so you can do whatever you want with it vis a vis matplotlib, seaborn, etc.
![example cnv plot](https://github.com/Adoni5/cnv_from_bam/blob/10a2b00a8832b46cacbff0e2f775a4f440844da0/example_cnv.png?raw=true)

## Output
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
