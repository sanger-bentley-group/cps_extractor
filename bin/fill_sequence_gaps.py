#!/usr/bin/env python3

import logging

from lib.argparser import GapFillerParser
from lib.gap_filler import GapFiller


def check_gaps(args, gap_filler) -> list:
    hits_list = gap_filler.read_hits_list()
    cps_cds_regions = gap_filler.get_cps_cds_regions()
    gaps = gap_filler.get_gaps(args.gap_length, hits_list)
    gaps_to_fill = gap_filler.check_gaps(cps_cds_regions, gaps)
    return gaps_to_fill


def fill_gaps(args, gap_filler, gaps_to_fill):
    gap_length = 0
    gap_added = 0
    print(gaps_to_fill)
    for i in range(0, len(gaps_to_fill)):
        gap_length = gap_filler.get_gap_length(gaps_to_fill[i])
        reference = gap_filler.get_sequence(gap_filler.reference)
        ref_to_map = gap_filler.subset_reference(reference, gaps_to_fill[i])
        gap_filler.write_subset_file(ref_to_map)
        sam_file = gap_filler.map_to_subset_ref("subset_ref.fa")
        bam_file = gap_filler.filter_mapping(sam_file)
        # check there is a useful number of reads from mapping before proceeding to assembly
        print(f"read count: {gap_filler.count_reads_bam(bam_file)}")
        if gap_filler.count_reads_bam(bam_file) < args.minimum_reads:
            logging.info(f"Too few reads for {gaps_to_fill[i]}")
            continue

        # filtered_reads = gap_filler.bam_to_fastq(bam_file)
        consensus_sequence = gap_filler.samtools_consensus(bam_file)
        consensus_seq = gap_filler.get_sequence(consensus_sequence)

        # basic check to see if there is a useful number of bases in the consensus to fill a gap
        if len(str(consensus_seq)) < (gap_length / 2):
            logging.info(
                f"very low number of bases in consensus sequence for {gaps_to_fill[i]}: {len(str(consensus_seq))}"
            )
            continue

        gap_filling_seq = gap_filler.filter_consensus_seq(
            consensus_seq, gaps_to_fill, i
        )
        if i == 0:
            cps_seq = gap_filler.get_sequence(gap_filler.cps_sequence)
        else:
            cps_seq = gap_filler.get_sequence("gap_filled_seq.fa")

        gap_fill_seq = gap_filler.fill_sequence_gap(
            cps_seq, gap_filling_seq, gaps_to_fill, i, gap_added
        )

        gap_added += len(str(gap_filling_seq))
        print(f"gap added: {gap_added}")
        with open("gap_filled_seq.fa", "w") as f:
            f.write(">gap_filled_seq\n")
            f.write(gap_fill_seq)


def main(args):
    gap_filler = GapFiller(
        args.log_file,
        args.annotation,
        args.reference,
        args.read_1,
        args.read_2,
        args.input_sequence,
    )

    gaps_to_fill = check_gaps(args, gap_filler)

    if len(gaps_to_fill) > 0:
        fill_gaps(args, gap_filler, gaps_to_fill)


if __name__ == "__main__":
    args = GapFillerParser.parse_args(vargs=None)
    main(args)
