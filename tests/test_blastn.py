import pytest
import logging

from lib.blastn import Blast


@pytest.fixture
def blast():
    # this blast fixture has a truncated sequence in it for readability
    blast_tester = Blast("tests/test_data/blast_1_hit.xml", 2500)
    return blast_tester


@pytest.fixture
def blast_real_results():
    real_blast = Blast("tests/test_data/real_blast_results.xml", 2500)
    result = real_blast.parse_blast_results()
    result_no_sequence = real_blast.parse_blast_results_dev(result)
    return result_no_sequence


@pytest.fixture
def blast_2_hits():
    # this blast fixture has a truncated sequence in it for readability
    blast_tester = Blast("tests/test_data/blast_2_hits.xml", 2500)
    return blast_tester


@pytest.fixture
def blast_result_1_hit_reverse(blast):
    result = blast.parse_blast_results()
    return result


@pytest.fixture
def overlapping_blast_results():
    overlapping_results = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "query_id": ".21127_1_98.1",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "ATATAGTAA",
            "reference_length": 13844,
            "e_value": 1.02314e-45,
        },
        {
            "hit_start": 1992,
            "hit_end": 2000,
            "hit_frame": 1,
            "seq_length": 9,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CAATAATGTCACG",
            "reference_length": 13844,
            "e_value": 0.0,
        },
    ]
    return overlapping_results


def test_parse_blast_results_1_hit(blast):
    blast_results = blast.parse_blast_results()
    assert blast_results == [
        {
            "hit_start": 13844,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 13844,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CAATAATGTCACG",
            "reference_length": 13844,
            "e_value": 0.0,
        }
    ]


@pytest.fixture
def blast_results_no_overlap(blast_2_hits):
    blast_results = blast_2_hits.parse_blast_results()
    return blast_results


def test_parse_blast_results_2_hits(blast_2_hits):
    blast_results = blast_2_hits.parse_blast_results()
    assert blast_results == [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": -1,
            "seq_length": 9892,
            "query_id": ".21127_1_98.1",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "ATATAGTAA",
            "reference_length": 13844,
            "e_value": 1.02314e-45,
        },
        {
            "hit_start": 13844,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 13844,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CAATAATGTCACG",
            "reference_length": 13844,
            "e_value": 0.0,
        },
    ]


def test_do_dicts_overlap(overlapping_blast_results, blast):
    overlap = blast.do_dicts_overlap(
        overlapping_blast_results[0], overlapping_blast_results[1]
    )
    assert overlap


def test_do_dicts_not_overlap(overlapping_blast_results, blast):
    overlap = blast.do_dicts_overlap(
        overlapping_blast_results[1], overlapping_blast_results[0]
    )
    assert not overlap


def test_get_best_hit_wrong_serotype(blast, overlapping_blast_results, caplog):
    serotype = "non_existent_serotype"
    with pytest.raises(SystemExit) as exit_info:
        with caplog.at_level(logging.ERROR):
            blast.get_best_hit(overlapping_blast_results, serotype)
    assert (
        "No results found for non_existent_serotype, please check the blast results file"
        in caplog.text
    )
    assert exit_info.value.code == 1


def test_get_best_hit_empty_blast_result(blast, caplog):
    serotype = "non_existent_serotype"
    empty_blast_results = list()
    with pytest.raises(SystemExit) as exit_info:
        with caplog.at_level(logging.ERROR):
            blast.get_best_hit(empty_blast_results, serotype)
    assert (
        "The blast results file is empty, please check your blast database and input sequence"
        in caplog.text
    )
    assert exit_info.value.code == 1


def test_get_best_hit_known_sero(blast, blast_real_results):
    serotype = "16F"
    result = blast.get_best_hit(blast_real_results, serotype)
    assert result == serotype


def test_get_best_hit_unknown_sero(blast, blast_real_results):
    actual_serotype = "19F"
    result = blast.get_best_hit(blast_real_results, None)
    assert result == actual_serotype


def test_compare_blast_dicts_no_overlap_same_sero(blast):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 1,
            "hit_end": 10,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
    ]
    final_result = blast.compare_blast_dicts(result, serotype)
    assert final_result == result[::-1]


def test_compare_blast_dicts_reverse_hits(blast):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": -1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 1,
            "hit_end": 10,
            "hit_frame": -1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
    ]
    final_result = blast.compare_blast_dicts(result, serotype)
    assert final_result == result[::-1]


def test_compare_blast_dicts_overlap_same_sero(blast):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 10000,
            "hit_end": 10010,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
    ]
    final_result = blast.compare_blast_dicts(result, serotype)
    assert len(final_result) == 1
    assert final_result[0] == result[0]


def test_compare_blast_dicts_bad_e_value(blast, caplog):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-5,
        },
        {
            "hit_start": 1,
            "hit_end": 10,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-5,
        },
    ]
    with pytest.raises(SystemExit) as exit_info:
        with caplog.at_level(logging.ERROR):
            blast.compare_blast_dicts(result, serotype)
    assert (
        "No results found for 19F, please check the blast results file" in caplog.text
    )
    assert exit_info.value.code == 1


def test_compare_blast_dicts_different_sero(blast):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 10000,
            "hit_end": 10010,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "01",
            "e_value": 1.02314e-55,
        },
    ]
    final_result = blast.compare_blast_dicts(result, serotype)
    assert len(final_result) == 1
    assert final_result[0] == result[0]


def test_compare_blast_dicts_multi_overlap(blast):
    serotype = "19F"
    result = [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 10000,
            "hit_end": 10010,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 100000,
            "hit_end": 100010,
            "hit_frame": 1,
            "seq_length": 9,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 99995,
            "hit_end": 100100,
            "hit_frame": 1,
            "seq_length": 100,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
    ]
    final_result = blast.compare_blast_dicts(result, serotype)
    assert final_result == [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
        {
            "hit_start": 99995,
            "hit_end": 100100,
            "hit_frame": 1,
            "seq_length": 100,
            "hit_def": "19F",
            "e_value": 1.02314e-55,
        },
    ]


def test_check_partial_overlap(blast, overlapping_blast_results):
    overlap = blast.check_partial_overlap(
        overlapping_blast_results[0], overlapping_blast_results[1]
    )
    assert overlap


def test_check_partial_overlap_no_overlap(blast, blast_results_no_overlap):
    overlap = blast.check_partial_overlap(
        blast_results_no_overlap[0], blast_results_no_overlap[1]
    )
    assert not overlap


def test_reverse_complement(blast):
    seq = "ATTTGGGCC"
    rev_comp = blast.reverse_complement(seq)
    assert rev_comp == "GGCCCAAAT"


def test_reverse_complement_with_gaps(blast):
    seq = "ATT-GGG-C"
    rev_comp = blast.reverse_complement(seq)
    assert rev_comp == "G-CCC-AAT"


def test_reverse_complement_n_and_gaps(blast):
    seq = "ANT-GNG-C"
    rev_comp = blast.reverse_complement(seq)
    assert rev_comp == "G-CNC-ANT"


def test_reverse_complement_hits(blast, blast_result_1_hit_reverse):
    result = blast.reverse_complement_hits(blast_result_1_hit_reverse)
    assert result == [
        {
            "hit_start": 13844,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 13844,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CGTGACATTATTG",
            "reference_length": 13844,
            "e_value": 0.0,
        }
    ]


def test_reverse_complement_hits_fwd(blast, blast_result_1_hit_reverse):
    # change hit to forward, check that the sequence is not reverse complemented
    blast_result_fwd = blast_result_1_hit_reverse.copy()
    blast_result_fwd[0]["hit_frame"] = 1
    result = blast.reverse_complement_hits(blast_result_fwd)
    assert result == blast_result_fwd


def test_reverse_complement_multi_hits(blast):
    multi_hits = [
        {
            "hit_start": 13844,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 13844,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CGTGACATTATTG",
            "reference_length": 13844,
            "e_value": 0.0,
        },
        {
            "hit_start": 13844,
            "hit_end": 1,
            "hit_frame": -1,
            "seq_length": 13844,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "seq": "CGCATTATTG",
            "reference_length": 13844,
            "e_value": 0.0,
        },
    ]

    result = blast.reverse_complement_hits(multi_hits)
    assert result[0]["seq"] == "CAATAATGTCACG"
    assert result[1]["seq"] == "CAATAATGCG"


def test_curate_sequence_empty_results(blast, caplog):
    empty_results = list()
    with pytest.raises(SystemExit) as exit_info:
        with caplog.at_level(logging.ERROR):
            blast.curate_sequence(empty_results)
    assert "No blast hits found, please check the blast XML file" in caplog.text
    assert exit_info.value.code == 1


def test_curate_sequence_1_hit(blast, blast_result_1_hit_reverse):
    sequence = blast.curate_sequence(blast_result_1_hit_reverse)
    assert sequence == blast_result_1_hit_reverse[0]["seq"]


def test_curate_sequence_2_hit_no_overlap(blast):
    non_overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
        {
            "hit_start": 20,
            "hit_end": 33,
            "seq_length": 13,
            "seq": "CAATAATGTCACG",
        },
    ]

    sequence = blast.curate_sequence(non_overlapping_blast_results)
    assert sequence == "ATATAGTAACAATAATGTCACG"


def test_curate_sequence_2_hit_overlap(blast):
    overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
        {
            "hit_start": 8,
            "hit_end": 21,
            "seq_length": 13,
            "seq": "CAATAATGTCACG",
        },
    ]

    sequence = blast.curate_sequence(overlapping_blast_results)
    assert sequence == "ATATAGTCAATAATGTCACG"


def test_curate_sequence_2_hit_overlap_reverse(blast):
    overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 21,
            "seq_length": 13,
            "seq": "CAATAATGTCACG",
        },
        {
            "hit_start": 20,
            "hit_end": 29,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
    ]

    sequence = blast.curate_sequence(overlapping_blast_results)
    assert sequence == "CAATAATGTCACGATAGTAA"


def test_curate_sequence_3_hit_overlap(blast):
    overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
        {
            "hit_start": 8,
            "hit_end": 21,
            "seq_length": 13,
            "seq": "CAATAATGTCACG",
        },
        {
            "hit_start": 21,
            "hit_end": 23,
            "seq_length": 3,
            "seq": "CAA",
        },
    ]

    sequence = blast.curate_sequence(overlapping_blast_results)
    assert sequence == "ATATAGTCAATAATGTCACGAA"


def test_curate_sequence_3_hit_overlap_reverse(blast):
    overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
        {
            "hit_start": 11,
            "hit_end": 14,
            "seq_length": 3,
            "seq": "CAA",
        },
        {
            "hit_start": 12,
            "hit_end": 26,
            "seq_length": 13,
            "seq": "TTTTAATGTCACG",
        },
    ]

    sequence = blast.curate_sequence(overlapping_blast_results)
    assert sequence == "ATATAGTAACTTTTAATGTCACG"


def test_curate_sequence_3_hit_2_overlap(blast):
    overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATAGTAA",
        },
        {
            "hit_start": 11,
            "hit_end": 14,
            "seq_length": 3,
            "seq": "CAA",
        },
        {
            "hit_start": 15,
            "hit_end": 29,
            "seq_length": 13,
            "seq": "TTTTAATGTCACG",
        },
    ]

    sequence = blast.curate_sequence(overlapping_blast_results)
    assert sequence == "ATATAGTAACAATTTTAATGTCACG"


def test_curate_sequence_remove_n_gaps(blast):
    non_overlapping_blast_results = [
        {
            "hit_start": 1,
            "hit_end": 10,
            "seq_length": 9,
            "seq": "ATATANTAA",
        },
        {
            "hit_start": 20,
            "hit_end": 33,
            "seq_length": 13,
            "seq": "CAA-AATGT-ACG",
        },
    ]

    sequence = blast.curate_sequence(non_overlapping_blast_results)
    assert sequence == "ATATATAACAAAATGTACG"


def test_write_fasta(blast, tmp_path):
    output_file = f"{tmp_path}/test_cps.fa"
    sequence = "AAATTTCCCGGG"
    blast.write_fasta(sequence, output_file)

    with open(output_file, "r") as output:
        lines = output.readlines()
        assert len(lines) == 2
        assert lines[0] == ">test_cps\n"
        assert lines[1] == sequence


def test_parse_blast_results_dev(blast, overlapping_blast_results):
    result = blast.parse_blast_results_dev(overlapping_blast_results)
    assert result == [
        {
            "hit_start": 1449,
            "hit_end": 11340,
            "hit_frame": 1,
            "seq_length": 9892,
            "query_id": ".21127_1_98.1",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "reference_length": 13844,
            "e_value": 1.02314e-45,
        },
        {
            "hit_start": 1992,
            "hit_end": 2000,
            "hit_frame": 1,
            "seq_length": 9,
            "query_id": ".21127_1_98.2",
            "hit_def": "ENA|CR931644|CR931644.1 Streptococcus pneumoniae strain 573/62 (serotype 8).",
            "reference_length": 13844,
            "e_value": 0.0,
        },
    ]
