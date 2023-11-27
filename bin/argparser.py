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