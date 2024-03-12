import pytest
from cnv_from_bam import iterate_bam_file
from pathlib import Path
from pprint import pprint


@pytest.fixture
def bam_file():
    return "/home/adoni5/Documents/Bioinformatics/test_data/HG00729.sorted.bam"


@pytest.fixture
def small_bam_file():
    return "/home/adoni5//Documents/Bioinformatics/test_data/NA12878_4000_test.bam"


def iterate_bam(bam):
    print(bam)
    result = iterate_bam_file(Path(bam), log_level=30)
    print(len(result.cnv))


def test_iterate_bam_fixed_bin_width(small_bam_file):
    result = iterate_bam_file(
        Path(small_bam_file), mapq_filter=60, log_level=20, ploidy=2, bin_width=10000
    )
    # Chr1 is 248Mb, which should be 2490 10kb bins
    pprint(result.cnv["NC_000001.11"][:10])
    assert len(result.cnv["NC_000001.11"]) == 24896, "Unexpected width of bins"


def test_iterate_bam_update_things_fixed_bin_size(small_bam_file):
    update = {}
    test_sum = [(0, 0)]
    for i in range(10):
        result = iterate_bam_file(
            Path(small_bam_file),
            mapq_filter=60,
            log_level=20,
            copy_numbers=update,
            bin_size=10000,
        )
        test_sum.append((sum(update["NC_000001.11"]), result.cnv["NC_000001.11"]))
    pprint(test_sum)


def test_iterate_bam_update_things(small_bam_file):
    update = {}
    test_sum = [(0, 0)]
    for i in range(10):
        result = iterate_bam_file(
            Path(small_bam_file), mapq_filter=60, log_level=20, copy_numbers=update
        )
        test_sum.append((sum(update["NC_000001.11"]), result.cnv["NC_000001.11"]))
    pprint(test_sum)


def test_iterate_bam_file(benchmark, bam_file):
    benchmark.pedantic(iterate_bam, args=(bam_file,), iterations=2, rounds=2)
