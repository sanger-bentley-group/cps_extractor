from BCBio import GFF
import pandas as pd


class GeneticVariants:
    def get_features(self, annotation_file) -> dict:
        # get key features out of annotation
        gene_info = list()

        with open(annotation_file) as annotation:
            for rec in GFF.parse(annotation):
                for feature in rec.features:
                    if feature.type != "region":
                        if "gene" in feature.qualifiers.keys():
                            gene_name = " ".join(feature.qualifiers["gene"])
                        else:
                            gene_name = "tnp"

                        product = " ".join(feature.qualifiers["product"])
                        if (
                            "transposase" in product.lower()
                            or "lipo" in product.lower()
                            or "intron" in product.lower()
                            or "chlorohydrolase" in product.lower()
                            or "family element" in product.lower()
                            or "is66" in product.lower()
                            or "srn266" in product.lower()
                        ):
                            gene_name = "tnp"
                        if "glf" in product.lower():
                            gene_name = "glf"
                        if "hypothetical" in product.lower():
                            gene_name = "tnp"
                        if gene_name.lower() == "oppa":
                            gene_name = "aliA"
                        if (
                            gene_name.lower().startswith("wzx")
                            and gene_name[3:].isdigit()
                        ):
                            gene_name = "wzx"
                        if (
                            gene_name.lower().startswith("wzy")
                            and gene_name[3:].isdigit()
                        ):
                            gene_name = "wzy"

                        gene_info.append(
                            {
                                "name": gene_name,
                                "product": product,
                            }
                        )
        return gene_info

    def assign_groups(self, groups_file):
        data = pd.read_csv(groups_file)
        # assign groups for unique genetic variants based on md5sum of gene list
        codes, uniques = pd.factorize(data["md5"])
        data["genetic_group"] = ["group" + str(i + 1) for i in codes]
        data = data.drop(columns=["md5"])
        data.to_csv("genetic_groups.csv", index=False)
        return data
