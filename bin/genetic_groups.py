#!/usr/bin/env python3

from lib.argparser import GeneticGroupParser
from lib.genetic_variants import GeneticVariants


def main(args):
    genetic_variants = GeneticVariants()
    genetic_variants.assign_groups(args.groups)


if __name__ == "__main__":
    args = GeneticGroupParser.parse_args(vargs=None)
    main(args)
