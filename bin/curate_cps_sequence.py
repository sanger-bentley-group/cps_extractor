#!/usr/bin/env python3
from lib.blastn import Blast
from lib.argparser import BlastParser


def main(args):
    blast = Blast(args.blast_results_file, args.genome, args.hit_length)

    blast_results = blast.parse_blast_results(args.blast_results_file)

    dexb_results = blast.parse_blast_results(args.dexb)

    alia_results = blast.parse_blast_results(args.alia)

    final_results = blast.compare_blast_dicts(blast_results, args.serotype)

    final_dexb_results = blast.compare_blast_dicts(dexb_results, args.serotype)

    final_alia_results = blast.compare_blast_dicts(alia_results, args.serotype)

    sorted_results = blast.reverse_complement_hits(final_results)

    sorted_dexb_results = blast.reverse_complement_hits(final_dexb_results)

    sorted_alia_results = blast.reverse_complement_hits(final_alia_results)

    sequence = blast.curate_sequence(
        sorted_results, sorted_dexb_results, sorted_alia_results
    )

    blast.write_fasta(sequence, args.output)

    blast.parse_blast_results_dev(blast_results)


if __name__ == "__main__":
    args = BlastParser.parse_args(vargs=None)
    main(args)
