#!/usr/bin/env python3

from lib.argparser import GeneticVariantParser
from lib.genetic_variants import GeneticVariants


def main(args):
    genetic_variants = GeneticVariants()
    features = genetic_variants.get_features(args.annotation)

    for i in features:
        print(i["name"])


if __name__ == "__main__":
    args = GeneticVariantParser.parse_args(vargs=None)
    main(args)
