#!/usr/bin/env python3
from blastn import Blast
from argparser import BlastParser


def main(args):
    blast = Blast(args.blast_results_file, args.hit_length)

    blast_results = blast.parse_blast_results()

    final_results = blast.compare_blast_dicts(blast_results)

    sorted_results = blast.sort_and_reverse_complement_hits(final_results)

    sequence = blast.curate_sequence(sorted_results)

    blast.write_fasta(sequence, args.output)

    blast.parse_blast_results_dev()


if __name__ == "__main__":
    args = BlastParser.parse_args(vargs=None)
    main(args)
