import pytest
import filecmp
import os

from lib.annotation import Annotation


@pytest.fixture
def annotator():
    annotator = Annotation("tests/test_data/test_cps.fa", "tests/test_data")
    return annotator


def test_check_sequence_completeness_no_mutations(annotator):
    sequence = "AAAAAAAAA"
    completeness = annotator.check_sequence_completeness(sequence, "sequence_id")
    assert completeness == (True, "")


def test_check_sequence_completeness_wrong_no_bases(annotator):
    sequence = "AAAAAAAAAA"
    completeness = annotator.check_sequence_completeness(sequence, "sequence_id")
    assert completeness == (False, "")


def test_check_sequence_completeness_stop_codon(annotator):
    sequence = "AAAAAATAA"
    completeness = annotator.check_sequence_completeness(sequence, "sequence_id")
    assert completeness == (True, "")


def test_check_sequence_completeness_mutation(annotator):
    sequence = "AAATAATAA"
    completeness = annotator.check_sequence_completeness(sequence, "sequence_id")
    assert completeness == (False, "stop codon TAA at position 4-6 in sequence_id")


def test_get_cds_annotations(annotator, tmp_path):
    annotator.get_cds_annotations("test.gff3", "out.gff3")
    test_gff = "tests/test_data/out.gff3"
    acc_gff = "tests/test_data/cds.gff3"
    assert filecmp.cmp(test_gff, acc_gff, shallow=False)
    # remove tmp file after comparison
    os.remove(test_gff)


def test_get_cds_fna(annotator):
    # need bedtools in path for this test
    annotator.get_cds_fna("cds.gff3", "test.fna", "out.fna")
    acc_cds_fna = "tests/test_data/test_cds.fna"
    test_cds_fna = "tests/test_data/out.fna"
    assert filecmp.cmp(test_cds_fna, acc_cds_fna, shallow=False)
    # remove tmp file after comparison
    os.remove(test_cds_fna)


def test_find_mutations_no_mutations(annotator):
    mutation_list = annotator.find_mutations("test_cds.fna")
    assert mutation_list == list()


def test_find_mutations(annotator):
    mutation_list = annotator.find_mutations("test_cds_mutation.fna")
    assert mutation_list == [{"contig_1:292-466(+)": "7-9"}]


def test_write_disruptive_mutations_file(annotator, tmp_path):
    mutation_list = annotator.find_mutations("test_cds_mutation.fna")
    mutation_file = f"{tmp_path}/mutations.csv"
    annotator.write_disruptive_mutations_file(mutation_file, mutation_list)
    with open(mutation_file, "r") as mutations:
        lines = mutations.readlines()
        assert len(lines) == 2
        assert lines[0] == "position,contig\n"
        assert lines[1] == "7-9,contig_1:292-466(+)\n"
