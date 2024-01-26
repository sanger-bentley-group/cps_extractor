import pytest
import os

from lib.gap_filler import GapFiller


@pytest.fixture
def gap_filler():
    gap_filler = GapFiller(
        "tests/test_data/test.log",
        "tests/test_data/test.gff3",
        "tests/test_data/ref.fa",
        "read_1.fastq.gz",
        "read_2.fastq.gz",
        "tests/test_data/cps_seq.fa",
    )
    return gap_filler


@pytest.fixture
def hits_list(gap_filler):
    hits_list = gap_filler.read_hits_list()
    return hits_list


@pytest.fixture
def cps_cds_regions(gap_filler):
    cps_cds_regions = gap_filler.get_cps_cds_regions()
    return cps_cds_regions


def test_sort_hits_list(gap_filler):
    unsorted_list = [
        {
            "hit_start": 10563,
            "hit_end": 7590,
            "hit_frame": -1,
            "seq_length": 2973,
            "query_id": ".15682_4_63.185",
            "hit_def": "12F",
        },
        {
            "hit_start": 5094,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 5093,
            "query_id": ".15682_4_63.10",
            "hit_def": "12F",
        },
    ]
    sorted_list = gap_filler.sort_hits_list(unsorted_list)
    assert sorted_list == unsorted_list[::-1]


def test_read_hits_list(gap_filler):
    hits_list = gap_filler.read_hits_list()
    assert hits_list == [
        {
            "hit_start": 5094,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 5093,
            "query_id": ".15682_4_63.10",
            "hit_def": "12F",
        },
        {
            "hit_start": 10563,
            "hit_end": 7590,
            "hit_frame": -1,
            "seq_length": 2973,
            "query_id": ".15682_4_63.185",
            "hit_def": "12F",
        },
        {
            "hit_start": 19081,
            "hit_end": 12214,
            "hit_frame": -1,
            "seq_length": 6867,
            "query_id": ".15682_4_63.76",
            "hit_def": "12F",
        },
    ]


def test_get_cps_cds_regions(gap_filler):
    cds_regions = gap_filler.get_cps_cds_regions()
    assert cds_regions == [
        {293: 466},
        {459: 782},
        {1560: 3005},
        {3007: 3738},
        {4449: 5138},
        {5153: 6520},
        {6525: 7004},
        {7005: 7484},
        {7497: 8561},
        {8561: 10021},
        {10034: 10762},
        {10746: 11909},
        {11978: 13216},
    ]


def test_get_gaps(gap_filler):
    hits_list = [
        {
            "hit_start": 1,
            "hit_end": 5094,
            "hit_frame": 1,
        },
        {
            "hit_start": 7590,
            "hit_end": 10563,
            "hit_frame": 1,
        },
    ]
    gaps = gap_filler.get_gaps(100, hits_list)
    assert gaps == [{5094: 7590}]


def test_do_dicts_overlap(gap_filler):
    dict1 = {0: 100}
    dict2 = {10: 20}
    assert gap_filler.do_dicts_overlap(dict2, dict1)


def test_do_dicts_not_overlap(gap_filler):
    dict1 = {0: 100}
    dict2 = {10: 20}
    assert not gap_filler.do_dicts_overlap(dict1, dict2)


def test_check_gaps(gap_filler, cps_cds_regions):
    gaps_list = [{4500: 5000}]
    gaps_checked = gap_filler.check_gaps(cps_cds_regions, gaps_list)
    assert gaps_checked == gaps_list


def test_check_gaps_multi(gap_filler, cps_cds_regions, hits_list):
    gaps = gap_filler.get_gaps(100, hits_list)
    gaps_checked = gap_filler.check_gaps(cps_cds_regions, gaps)
    assert gaps_checked == gaps


def test_check_gaps_no_gaps(gap_filler, cps_cds_regions):
    gap = [{100000: 110000}]
    gaps_checked = gap_filler.check_gaps(cps_cds_regions, gap)
    assert gaps_checked == list()


def test_get_sequence(gap_filler, tmp_path):
    # write tmp fasta file
    temp_fasta = f"{tmp_path}/tmp.fasta"
    with open(temp_fasta, "w") as fasta:
        fasta.write(">seq\n")
        fasta.write("AAACCCTTTGGG")

    # get sequence
    sequence = gap_filler.get_sequence(temp_fasta)
    assert sequence == "AAACCCTTTGGG"


def test_subset_reference(gap_filler):
    gaps_dict = {1: 5}
    ref_seq = "AAACCCTTTGGG"
    subset_ref = gap_filler.subset_reference(ref_seq, gaps_dict)
    assert subset_ref == "AAC"


def test_write_subset_file(gap_filler):
    subset_seq = "AAC"
    gap_filler.write_subset_file(subset_seq)
    with open("subset_ref.fa", "r") as f:
        data = f.readlines()
        assert len(data) == 2
        assert data[0] == ">subset_seq\n"
        assert data[1] == "AAC"
    os.remove("subset_ref.fa")


def test_get_gap_length(gap_filler):
    gaps = {1: 5}
    length = gap_filler.get_gap_length(gaps)
    assert length == 3


def test_filter_consensus_seq(gap_filler):
    gaps = [{1: 5}]
    seq = "AAACCCTTTGGG"
    filtered_seq = gap_filler.filter_consensus_seq(seq, gaps, 0)
    assert filtered_seq == "AAA"


def test_fill_sequence_gap(gap_filler):
    cps_seq = "AAACCCTTTGGG"
    gap_fill_seq = "GGG"
    gaps = [{3: 7}]
    full_seq = gap_filler.fill_sequence_gap(cps_seq, gap_fill_seq, gaps, 0, 0)
    assert full_seq == "AAAGGGCCCTTTGGG"


### missing unit tests for functions which call external programs (bwa, samtools)
### these functions probably will change when the gap filling code is integrated into the main pipeline
### if they don't, write some tests!
