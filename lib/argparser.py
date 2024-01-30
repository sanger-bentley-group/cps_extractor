#!/usr/bin/env python3
import argparse


class BlastParser:
    @classmethod
    def parse_args(cls, vargs=None):
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        required = parser.add_argument_group("required")
        required.add_argument(
            "-b",
            "--blast-results-file",
            required=True,
            help="Path to blast results file",
        )
        required.add_argument(
            "-o",
            "--output",
            required=True,
            help="Name of output cps sequence file",
        )
        optional = parser.add_argument_group("optional")
        optional.add_argument(
            "-l",
            "--hit-length",
            required=False,
            help="Length of blast hits to add to the cps sequence",
            default=2500,
            type=int,
        )

        optional.add_argument(
            "-s",
            "--serotype",
            required=False,
            help="Provide the serotype if it is known",
            default=None,
        )

        args = parser.parse_args(vargs)

        return args


class AnnotationParser:
    @classmethod
    def parse_args(cls, vargs=None):
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        required = parser.add_argument_group("required")
        required.add_argument(
            "-c",
            "--cps-sequence",
            required=True,
            help="Path to cps sequence file",
        )
        required.add_argument(
            "-b",
            "--bakta-input",
            required=True,
            help="Path to bakta input folder",
        )

        args = parser.parse_args(vargs)

        return args


class GapFillerParser:
    @classmethod
    def parse_args(cls, vargs=None):
        parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
        required = parser.add_argument_group("required")
        required.add_argument(
            "-l",
            "--log-file",
            required=True,
            help="Path to cps extractor log file",
        )
        required.add_argument(
            "-a",
            "--annotation",
            required=True,
            help="Path to reference annotation file",
        )
        required.add_argument(
            "-r1",
            "--read-1",
            required=True,
            help="Path to first read file",
        )
        required.add_argument(
            "-r2",
            "--read-2",
            required=True,
            help="Path to second read file",
        )
        required.add_argument(
            "-r",
            "--reference",
            required=True,
            help="Path to reference fasta file",
        )
        required.add_argument(
            "-i",
            "--input-sequence",
            required=True,
            help="Path to input cps sequence file",
        )
        optional = parser.add_argument_group("optional")
        optional.add_argument(
            "-g",
            "--gap-length",
            required=False,
            help="Minimum length of gaps that are to be filled",
            default=100,
            type=int,
        )

        optional.add_argument(
            "-m",
            "--minimum-reads",
            required=False,
            help="Minimum number of mapped reads required to fill gap",
            default=500,
            type=int,
        )

        args = parser.parse_args(vargs)

        return args
