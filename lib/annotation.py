#!/usr/bin/env python3

import logging
import csv
import os
import pybedtools
from Bio import SeqIO


class Annotation:
    def __init__(self, cps_fasta: str, bakta_outdir: str):
        self.cps_fasta = cps_fasta
        self.bakta_outdir = bakta_outdir
        self.reference_basename = self.cps_fasta.split(".fa")[0]

    def check_seq_length(self) -> int:
        # basic check to see if the CPS sequence length is at least 10000 bases
        sequence = SeqIO.parse(self.cps_fasta, "fasta")
        for record in sequence:
            return len(record.seq)

    def check_sequence_completeness(
        self, sequence: str, sequence_id: str
    ) -> (bool, str):
        mutation = str()
        stop_codons = {"TAA", "TGA", "TAG"}
        complete = True
        if len(sequence) % 3 != 0:
            complete = False
        i = 0
        # -6 so you don't check the last codon (which should be a stop codon)
        while i <= len(sequence) - 6 and complete:
            codon = sequence[i : i + 3]
            if codon in stop_codons:
                mutation = (
                    f"stop codon {codon} at position {i + 1}-{i + 3} in {sequence_id}"
                )
                logging.info(f"{mutation}")
                complete = False
            i += 3
        return complete, mutation

    def get_cds_annotations(self, input_gff: str, output_gff: str) -> str:
        with open(f"{self.bakta_outdir}/{input_gff}", "r", newline="") as gff_file:
            gff_reader = csv.reader(gff_file, delimiter="\t")
            with open(f"{self.bakta_outdir}/{output_gff}", "w") as cds_gff:
                for row in gff_reader:
                    if len(row) >= 3 and row[2] == "CDS":
                        cds_gff.write("\t".join(row))
                        cds_gff.write("\n")
        return output_gff

    def get_cds_fna(self, input_gff: str, input_fna: str, output_fna: str) -> str:
        wd = os.getcwd()
        cds = pybedtools.BedTool(f"{self.bakta_outdir}/{input_gff}")
        # hack because pybedtools uses its temp dir as pwd so can't find input files
        fasta = pybedtools.example_filename(f"{wd}/{self.bakta_outdir}/{input_fna}")
        cds = cds.sequence(fi=fasta, s=True)
        with open(cds.seqfn, "r") as cds_seqs:
            with open(f"{self.bakta_outdir}/{output_fna}", "w") as cds_seqs_fa:
                cds_seqs_fa.write(cds_seqs.read())
        return output_fna

    def find_mutations(self, cds_fna: str) -> list:
        disruptive_mutations_list = list()
        with open(f"{self.bakta_outdir}/{cds_fna}", "r") as cds_fna:
            fasta_sequences = SeqIO.parse(cds_fna, "fasta")
            for fasta in fasta_sequences:
                name = fasta.id
                sequence = str(fasta.seq)
                disruptive_mutations = self.check_sequence_completeness(sequence, name)
                if disruptive_mutations[1]:
                    disruptive_mutations = disruptive_mutations[1].split(" ")
                    disruptive_mutations_dict = {
                        disruptive_mutations[-1]: disruptive_mutations[5]
                    }
                    disruptive_mutations_list.append(disruptive_mutations_dict)
        return disruptive_mutations_list

    def write_disruptive_mutations_file(self, mutations_file: str, mutations: list):
        with open(mutations_file, "w") as mut:
            mut.write("position,contig\n")
            for i in range(0, len(mutations)):
                mutation_data = list(mutations[i].items())
                mut.write(f"{mutation_data[0][1]},{mutation_data[0][0]}\n")
