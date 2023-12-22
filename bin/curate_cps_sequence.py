#!/usr/bin/env python3
from lib.blastn import Blast
from lib.argparser import BlastParser


def main(args):
    blast = Blast(args.blast_results_file, args.hit_length)

    blast_results = blast.parse_blast_results()

    final_results = blast.compare_blast_dicts(blast_results, args.serotype)

    sorted_results = blast.reverse_complement_hits(final_results)

    sequence = blast.curate_sequence(sorted_results)

    blast.write_fasta(sequence, args.output)

    blast.parse_blast_results_dev(blast_results)


if __name__ == "__main__":
    args = BlastParser.parse_args(vargs=None)
    main(args)
