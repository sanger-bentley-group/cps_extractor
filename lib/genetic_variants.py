from BCBio import GFF


class GeneticVariants:
    def __init__(self, annotation_file: str):
        self.annotation_file = annotation_file

    def get_features(self) -> dict:
        # get key features out of annotation
        gene_info = list()

        with open(self.annotation_file) as annotation:
            for rec in GFF.parse(annotation):
                for feature in rec.features:
                    if feature.type != "region":
                        if "gene" in feature.qualifiers.keys():
                            gene_name = " ".join(feature.qualifiers["gene"])
                        else:
                            gene_name = "unidentified"

                        product = " ".join(feature.qualifiers["product"])
                        if (
                            "transposase" in product.lower()
                            or "lipo" in product.lower()
                            or "intron" in product.lower()
                            or "chlorohydrolase" in product.lower()
                            or "family element" in product.lower()
                            or "is66" in product.lower()
                        ):
                            gene_name = "tnp"
                        if "glf" in product.lower():
                            gene_name = "glf"
                        if "hypothetical" in product.lower():
                            gene_name = "hypothetical protein"
                        if gene_name.lower() == "oppa":
                            gene_name = "aliA"

                        gene_info.append(
                            {
                                "name": gene_name,
                                "product": product,
                            }
                        )
        return gene_info
