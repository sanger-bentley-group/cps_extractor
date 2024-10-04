#!/usr/bin/env python3

import logging

from lib.annotation import Annotation
from lib.argparser import AnnotationParser

logging.basicConfig(
    filename=f"cps_extractor.log",
    encoding="utf-8",
    level=logging.ERROR,
)


def main(args):
    Annotator = Annotation(args.cps_sequence, args.bakta_input)

    cps_seq_length = Annotator.check_seq_length()
    # basic sanity check to make sure there is a genuine blast hit to the reference DB
    if cps_seq_length < args.min_length:
        logging.error(
            f"The CPS sequence length is unusually low ({cps_seq_length} bases), please check the blast results file you may have a non capsulated sample or a pneumo 'like' sample"
        )
        raise SystemExit(1)

    sample_name = args.cps_sequence.split(".fa")[0]
    cds_gff = Annotator.get_cds_annotations(
        f"{sample_name}.gff3", f"{sample_name}_cds.gff3"
    )

    cds_fna = Annotator.get_cds_fna(
        cds_gff, f"{sample_name}.fna", f"{sample_name}_cds.fna"
    )

    mutations = Annotator.find_mutations(cds_fna)
    Annotator.write_disruptive_mutations_file(f"{sample_name}_mutations.csv", mutations)


if __name__ == "__main__":
    args = AnnotationParser.parse_args(vargs=None)
    main(args)
