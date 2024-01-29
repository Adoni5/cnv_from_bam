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
