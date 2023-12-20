#!/usr/bin/env python3

from lib.annotation import Annotation
from lib.argparser import AnnotationParser

def main(args):
    Annotator = Annotation(args.cps_sequence, args.bakta_input)

    sample_name = args.cps_sequence.split(".fa")[0]
    cds_gff = Annotator.get_cds_annotations(
        f"{sample_name}.gff3", f"{sample_name}_cds.gff3"
    )

    cds_fna = Annotator.get_cds_fna(
        cds_gff, f"{sample_name}.fna", f"{sample_name}_cds.fna"
    )

    mutations = Annotator.find_mutations(cds_fna)
    Annotator.write_disruptive_mutations_file(
        f"{sample_name}_mutations.csv", mutations
    )

if __name__ == "__main__":
    args = AnnotationParser.parse_args(vargs=None)
    main(args)