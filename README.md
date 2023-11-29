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
let result = iterate_bam_file(bam_path, Some(4), Some(60));
// Process the result...
```

The results in this case are returned as a CnvResult, which has the following structure:

```rust
pub struct CnvResult {
    pub cnv: FnvHashMap<String, Vec<f64>>,
    pub bin_width: usize,
    pub genome_length: usize,
}
```

Where `result.cnv` is a hash map containing the Copy Number for each bin of `bin_width` bases for each contig in the reference genome, `result.bin_width` is the width of the bins in bases, and `result.genome_length` is the total length of the genome.

> [!NOTE]
> **Note**: Only the main primary mapping alignment start is binned, Supplementary and Secondary alignments are ignored.

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
We use PyO3-log for rust/python interop logging. By default, the log level is set to `INFO`.

> [!WARNING]
> It is required to set up a logger before a call to `iterate_bam_file` is made. If no logger is set up, the program will not output anything. To set up a logger, run the following code before calling `iterate_bam_file`:

```python
import logging
import sys
FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("cnv_from_bam")
logger.handlers.clear()
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
```

It is possible to hide the log messages by setting the log level to `WARN`:

```python
import logging
import sys
FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger("cnv_from_bam")
logger.handlers.clear()
logger.setLevel(logging.WARN)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
```

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



## License

This project is licensed under the [Mozilla Public License 2.0](LICENSE).
