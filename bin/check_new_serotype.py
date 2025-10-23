#!/usr/bin/env python3
from lib.check_new_serotype import NewSerotype
from lib.argparser import NewSerotypeParser


def main(args):
    new_serotype = NewSerotype(
        args.serotype, args.known_disruptions, args.disrupted_genes
    )

    known_disruptions_df = new_serotype.read_csv(new_serotype.known_disruptions)

    disrupted_gene_df = new_serotype.read_csv(new_serotype.disrupted_genes)

    filtered_df = new_serotype.remove_known_disruptions(
        disrupted_gene_df, known_disruptions_df
    )

    unique_disruptions = new_serotype.find_unique_samples(filtered_df)

    for sample in unique_disruptions:
        print(sample)


if __name__ == "__main__":
    args = NewSerotypeParser.parse_args(vargs=None)
    main(args)
