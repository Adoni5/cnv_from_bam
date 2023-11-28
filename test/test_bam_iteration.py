import pytest
from cnv_from_bam import iterate_bam_file
from pathlib import Path


@pytest.fixture
def bam_file():
    return "/home/adoni5/Documents/Bioinformatics/test_data/HG00729.sorted.bam"


def iterate_bam(bam):
    print(bam)
    result = iterate_bam_file(Path(bam))
    print(len(result.cnv))


def test_iterate_bam_file(benchmark, bam_file):
    benchmark.pedantic(iterate_bam, args=(bam_file,), iterations=2, rounds=2)
