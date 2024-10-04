import pytest

from lib.argparser import (
    AnnotationParser,
    BlastParser,
    GapFillerParser,
    GeneOrderParser,
)


def test_annotation_args_valid():
    vargs = ["-c", "cps_file.txt", "-b", "bakta_input_folder", "-m", "8000"]
    args = AnnotationParser.parse_args(vargs)

    assert args.cps_sequence == "cps_file.txt"
    assert args.bakta_input == "bakta_input_folder"
    assert args.min_length == 8000


def test_annotation_args_missing_required_argument():
    vargs = ["-c", "cps_file.txt"]

    with pytest.raises(SystemExit):
        AnnotationParser.parse_args(vargs)


def test_annotation_args_invalid_args():
    vargs = ["-c", "cps_file.txt", "-b", "bakta_input_folder", "--custom", "value"]

    with pytest.raises(SystemExit):
        AnnotationParser.parse_args(vargs)


def test_blast_args_valid():
    vargs = [
        "-b",
        "blast_results.txt",
        "-o",
        "output_file.txt",
        "-l",
        "2000",
        "-s",
        "01",
    ]
    args = BlastParser.parse_args(vargs)

    assert args.blast_results_file == "blast_results.txt"
    assert args.output == "output_file.txt"
    assert args.hit_length == 2000
    assert args.serotype == "01"


def test_blast_args_missing_required_argument():
    vargs = ["-b", "blast_results.txt"]

    # Use pytest.raises to catch the argparse error when a required argument is missing
    with pytest.raises(SystemExit):
        BlastParser.parse_args(vargs)


def test_blast_args_default_values():
    vargs = ["-b", "blast_results.txt", "-o", "output_file.txt"]
    args = BlastParser.parse_args(vargs)

    assert args.hit_length == 2500
    assert args.serotype is None


def test_blast_args_invalid_hit_length():
    vargs = ["-b", "blast_results.txt", "-o", "output_file.txt", "-l", "invalid_value"]

    with pytest.raises(SystemExit):
        BlastParser.parse_args(vargs)


def test_blast_args_no_args():
    with pytest.raises(SystemExit):
        BlastParser.parse_args()


def test_blast_args_invalid_args():
    vargs = ["-x", "invalid_argument"]

    with pytest.raises(SystemExit):
        BlastParser.parse_args(vargs)


def test_gap_filler_args_valid():
    vargs = [
        "-l",
        "cps_extractor.log",
        "-a",
        "reference_annotation.gff",
        "-r1",
        "read1.fastq",
        "-r2",
        "read2.fastq",
        "-r",
        "reference.fasta",
        "-i",
        "input_cps_sequence.txt",
        "-g",
        "150",
        "-m",
        "75",
    ]
    args = GapFillerParser.parse_args(vargs)

    assert args.log_file == "cps_extractor.log"
    assert args.annotation == "reference_annotation.gff"
    assert args.read_1 == "read1.fastq"
    assert args.read_2 == "read2.fastq"
    assert args.reference == "reference.fasta"
    assert args.input_sequence == "input_cps_sequence.txt"
    assert args.gap_length == 150
    assert args.minimum_reads == 75


def test_gap_filler_args_missing_required_argument():
    vargs = [
        "-l",
        "cps_extractor.log",
        "-a",
        "reference_annotation.gff",
        "-r1",
        "read1.fastq",
        "-r2",
        "read2.fastq",
        "-r",
        "reference.fasta",
        "-g",
        "150",
        "-m",
        "75",
    ]

    # Use pytest.raises to catch the argparse error when a required argument is missing
    with pytest.raises(SystemExit):
        GapFillerParser.parse_args(vargs)


def test_gap_filler_args_default_values():
    vargs = [
        "-l",
        "cps_extractor.log",
        "-a",
        "reference_annotation.gff",
        "-r1",
        "read1.fastq",
        "-r2",
        "read2.fastq",
        "-r",
        "reference.fasta",
        "-i",
        "input_cps_sequence.txt",
    ]
    args = GapFillerParser.parse_args(vargs)

    assert args.gap_length == 100
    assert args.minimum_reads == 500


def test_gap_filler_args_invalid_gap_length():
    vargs = [
        "-l",
        "cps_extractor.log",
        "-a",
        "reference_annotation.gff",
        "-r1",
        "read1.fastq",
        "-r2",
        "read2.fastq",
        "-r",
        "reference.fasta",
        "-i",
        "input_cps_sequence.txt",
        "-g",
        "invalid_value",
    ]

    # Use pytest.raises to catch the argparse error when an invalid value is provided for gap_length
    with pytest.raises(SystemExit):
        GapFillerParser.parse_args(vargs)


def test_gap_filler_args_invalid_minimum_reads():
    vargs = [
        "-l",
        "cps_extractor.log",
        "-a",
        "reference_annotation.gff",
        "-r1",
        "read1.fastq",
        "-r2",
        "read2.fastq",
        "-r",
        "reference.fasta",
        "-i",
        "input_cps_sequence.txt",
        "-m",
        "invalid_value",
    ]

    # Use pytest.raises to catch the argparse error when an invalid value is provided for minimum_reads
    with pytest.raises(SystemExit):
        GapFillerParser.parse_args(vargs)


def test_gap_filler_args_no_args():
    # Test with no arguments
    with pytest.raises(SystemExit):
        GapFillerParser.parse_args()


def test_gap_filler_args_invalid_args():
    vargs = ["-x", "invalid_argument"]

    # Use pytest.raises to catch the argparse error when an invalid argument is provided
    with pytest.raises(SystemExit):
        GapFillerParser.parse_args(vargs)


def test_gene_order_args_no_args():
    # Test with no arguments
    with pytest.raises(SystemExit):
        GeneOrderParser.parse_args()


def test_gene_order_args_invalid_args():
    vargs = ["-x", "invalid_argument"]

    # Use pytest.raises to catch the argparse error when an invalid argument is provided
    with pytest.raises(SystemExit):
        GeneOrderParser.parse_args(vargs)


def test_gene_order_args_valid():
    vargs = ["-i", "input.gff", "-r", "reference.gff", "-o", "output_folder"]
    args = GeneOrderParser.parse_args(vargs)

    assert args.input_annotation == "input.gff"
    assert args.reference_annotation == "reference.gff"
    assert args.output == "output_folder"
