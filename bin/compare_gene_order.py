#!/usr/bin/env python3

from lib.argparser import GeneOrderParser
from lib.check_gene_order import CheckGeneOrder
import pandas as pd


def main(args):
    gene_checker = CheckGeneOrder(args.input_annotation, args.reference_annotation)
    features = gene_checker.get_features(args.input_annotation, False)
    features_ref = gene_checker.get_features(args.reference_annotation, True)

    combined_data = gene_checker.combine_ref_sample_data(features, features_ref)

    combined_df = pd.DataFrame(combined_data)

    # fill in any blanks
    combined_df.fillna("NA", inplace=True)

    final_df = gene_checker.gene_comparison(combined_df)

    final_df_reordered = final_df[
        [
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
            "gene_in_ref_and_sample",
        ]
    ]

    final_df_reordered.to_csv(args.output, index=False)


if __name__ == "__main__":
    args = GeneOrderParser.parse_args(vargs=None)
    main(args)
