#!/usr/bin/env python3

from BCBio import GFF
import pandas as pd


class CheckGeneOrder:
    def __init__(self, sample_annotation: str, reference_annotation: str):
        self.sample_annotation = sample_annotation
        self.reference_annotation = reference_annotation

    def sort_for_uniref(self, item: str) -> tuple:
        if item.startswith("UniRef:UniRef100"):
            return (0, item)
        elif item.startswith("UniRef:UniRef90"):
            return (0, item)
        elif item.startswith("UniRef:UniRef50"):
            return (0, item)
        else:
            return (1, item)

    def get_features(self, annotation_file: str, ref: bool) -> dict:
        # get key features out of annotation
        gene_info = list()

        with open(annotation_file) as annotation:
            for rec in GFF.parse(annotation):
                for feature in rec.features:
                    if feature.type != "region":
                        if "gene" in feature.qualifiers.keys():
                            gene_name = " ".join(feature.qualifiers["gene"])
                        else:
                            gene_name = "unidentified"
                        try:
                            db_references = sorted(
                                feature.qualifiers["Dbxref"],
                                key=self.sort_for_uniref,
                            )
                            refseq_id = db_references[0].split("UniRef:")[-1]
                        except KeyError:
                            refseq_id = "NA"

                        product = " ".join(feature.qualifiers["product"])
                        if ref:
                            gene_info.append(
                                {
                                    "name_ref": gene_name,
                                    "uniref_id_ref": refseq_id,
                                    "product_ref": product,
                                    "location_ref": str(feature.location),
                                }
                            )
                        else:
                            gene_info.append(
                                {
                                    "name": gene_name,
                                    "uniref_id": refseq_id,
                                    "product": product,
                                    "location": str(feature.location),
                                }
                            )
        return gene_info

    def combine_ref_sample_data(self, sample_data: list, ref_data: list) -> list:
        if len(sample_data) >= len(ref_data):
            for i in range(0, len(sample_data)):
                if i < len(ref_data):
                    sample_data[i].update(ref_data[i])
            return sample_data
        else:
            for i in range(0, len(ref_data)):
                if i < len(sample_data):
                    ref_data[i].update(sample_data[i])
            return ref_data

    def check_gene_in_ref_and_sample(
        self, row: pd.Series, combined_df: pd.DataFrame
    ) -> str:
        # if the genes are not in the same order, is the gene somewhere else in the reference? (id)
        if row["uniref_id"] != row["uniref_id_ref"] and row["uniref_id"] != "NA":
            contains_ref = any(
                combined_df[
                    combined_df["uniref_id_ref"].str.contains(row["uniref_id"])
                ]["uniref_id_ref"]
            )
            if contains_ref:
                return "Y"
            # does the name appear elsewhere?
            elif row["name"] != "NA" and row["name"] != "unidentified":
                contains_name = any(
                    combined_df[combined_df["name_ref"].str.contains(row["name"])][
                        "name_ref"
                    ]
                )
                if contains_name:
                    return "Y"
                else:
                    return "N"
            # does the EXACT product appear elsewhere, assume same gene in such case - usually if name is "unidentified"
            elif row["product"] != "NA" and row["product"] != "unidentified":
                contains_product = any(
                    combined_df[
                        combined_df["product_ref"].str.contains(row["product"])
                    ]["product_ref"]
                )
                if contains_product:
                    return "Y"
                else:
                    return "N"
            else:
                return "N"
        else:
            if row["uniref_id"] != "NA":
                return "Y"
            else:
                return "N"

    def gene_comparison(self, combined_df: pd.DataFrame) -> pd.DataFrame:
        # direct gene content comparison
        comparison_df = combined_df.copy()
        comparison_df["equal_gene"] = "N"
        # are the names and IDs identical?
        comparison_df.loc[
            (comparison_df["uniref_id"] == comparison_df["uniref_id_ref"])
            & (comparison_df["uniref_id"] != "NA"),
            "equal_gene",
        ] = "Y"
        # are the products identical?
        comparison_df.loc[
            (comparison_df["product"] == comparison_df["product_ref"]),
            "equal_gene",
        ] = "Y"
        # if the genes are not the same, is the gene found elsewhere in the reference?
        comparison_df["gene_in_ref_and_sample"] = comparison_df.apply(
            self.check_gene_in_ref_and_sample, args=(comparison_df,), axis=1
        )
        # sanity check to ensure if the genes are equal it is always found in the ref
        comparison_df.loc[
            comparison_df["equal_gene"] == "Y",
            "gene_in_ref_and_sample",
        ] = "Y"
        return comparison_df
