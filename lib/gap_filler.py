#!/usr/bin/env python3

import ast
import csv
import subprocess
from Bio import SeqIO


class GapFiller:
    def __init__(
        self,
        logging_file: str,
        annotation: str,
        reference: str,
        read_1: str,
        read_2: str,
        cps_sequence: str,
    ):
        self.logging_file = logging_file
        self.annotation = annotation
        self.reference = reference
        self.read_1 = read_1
        self.read_2 = read_2
        self.cps_sequence = cps_sequence

    def read_hits_list(self) -> list:
        with open(self.logging_file, "r") as f:
            content = f.read()
            content = content.replace("INFO:root:", "")
            hits_list = ast.literal_eval(content)
            # remove sequence from dict as it no longer necessary
            for d in hits_list:
                del d["seq"]
            return hits_list

    def get_cps_cds_regions(self) -> list:
        # extract cps CDS regions
        coding_regions = list()
        # also remove transposable elements for gap filling purposes
        with open(self.annotation, "r", newline="") as gff_file:
            gff_reader = csv.reader(gff_file, delimiter="\t")
            for row in gff_reader:
                if (
                    len(row) >= 3
                    and row[2] == "CDS"
                    and "transposase" not in row[-1].lower()
                ):
                    coding_dict = {int(row[3]): int(row[4])}
                    coding_regions.append(coding_dict)
        return coding_regions

    def get_gaps(self, min_gap_length: int, hits_list: list) -> list:
        # identify gaps from the blast hits results >= n bases
        gaps = list()
        for i in range(0, (len(hits_list) - 1)):
            if (
                hits_list[i + 1]["hit_start"] - hits_list[i]["hit_end"]
                >= min_gap_length
            ):
                gaps_dict = {hits_list[i]["hit_end"]: hits_list[i + 1]["hit_start"]}
                gaps.append(gaps_dict)
        return gaps

    def do_dicts_overlap(self, dict1: dict, dict2: dict) -> bool:
        # check if a dict is contained within a larger dict
        overlap = False
        for k, v in dict1.items():
            for i, j in dict2.items():
                if k > i and v < j:
                    overlap = True
        return overlap

    def check_gaps(self, cps_cds_regions: list, gaps_list: list) -> list:
        # check if a coding region in the reference annotation is entirely contained within a gap in the blast hits\
        # or if a gap is entirely contained within a coding region, if it is return it
        gap_checking_regions = list()
        for i in range(0, len(gaps_list)):
            for j in range(0, len(cps_cds_regions)):
                # is a gap entirely contained within coding region of reference
                if self.do_dicts_overlap(gaps_list[i], cps_cds_regions[j]):
                    if gaps_list[i] not in gap_checking_regions:
                        gap_checking_regions.append(gaps_list[i])
                        break
                # does a gap overlap an entire coding region of the reference sequence
                elif self.do_dicts_overlap(cps_cds_regions[j], gaps_list[i]):
                    if gaps_list[i] not in gap_checking_regions:
                        gap_checking_regions.append(gaps_list[i])
                        break
                else:
                    for k, v in gaps_list[i].items():
                        if j < (len(cps_cds_regions) - 1):
                            # does a gap partially overlap two coding regions in the reference
                            if (
                                k > list(cps_cds_regions[j].items())[0][0]
                                and v > list(cps_cds_regions[j + 1].items())[0][0]
                                and v < list(cps_cds_regions[j + 1].items())[0][1]
                            ):
                                if gaps_list[i] not in gap_checking_regions:
                                    gap_checking_regions.append(gaps_list[i])
                                    break
        return gap_checking_regions

    def get_sequence(self, input_file: str) -> str:
        ref_seq = str()
        fasta = SeqIO.parse(input_file, "fasta")
        for record in fasta:
            ref_seq += record.seq
        return ref_seq

    def subset_reference(self, ref_seq: str, gaps_dict: dict) -> str:
        for k, v in gaps_dict.items():
            ref_seq = ref_seq[k : v - 1]
        return ref_seq

    def write_subset_file(self, subset_seq: str):
        with open("subset_ref.fa", "w") as f:
            f.write(">subset_seq\n")
            f.write(str(subset_seq))

    def map_to_subset_ref(self, subset_ref_file: str) -> str:
        sample_id = self.read_1.split("_1.fastq.gz")[0]
        sam_file = f"{sample_id}.sam"
        index_cmd = f"bwa index {subset_ref_file}"
        subprocess.check_output(index_cmd, shell=True)
        map_cmd = f"bwa mem {subset_ref_file} {self.read_1} {self.read_2} > {sam_file}"
        subprocess.check_output(map_cmd, shell=True)
        return sam_file

    def filter_mapping(self, sam_file: str) -> str:
        sample_id = sam_file.split(".sam")[0]
        bam_file = f"{sample_id}_filtered.bam"
        filter_cmd = f"samtools view -h -f 0x2 -q 60 -F 4 {sam_file} | samtools sort -o {bam_file}"
        subprocess.check_output(filter_cmd, shell=True)
        return bam_file

    def count_reads_bam(self, bam_file: str) -> int:
        read_count_cmd = f"samtools view -c {bam_file}"
        read_count = subprocess.check_output(read_count_cmd, shell=True)
        return int(read_count.strip())

    def bcftools_mpileup(self, bam_file: str) -> str:
        sample_id = bam_file.split("_filtered.bam")[0]
        bcftools_mpileup_cmd = f"bcftools mpileup -Ou -f subset_ref.fa {bam_file} | bcftools call --ploidy 1 -Ou -mv | bcftools norm -f subset_ref.fa -Oz -o {sample_id}.vcf.gz"
        subprocess.check_output(bcftools_mpileup_cmd, shell=True)

        return f"{sample_id}.vcf.gz"

    def bcftools_filter(self, vcf_file: str) -> str:
        sample_id = vcf_file.split(".vcf.gz")[0]
        bcftools_filter_cmd = f"bcftools view -i 'QUAL >= 50 || DP > 100 || MQBZ > -3 || RPBZ > -3 || RPBZ < 3 || SCBZ < 3' {vcf_file} | bgzip > {sample_id}_filtered.vcf.gz"
        subprocess.check_output(bcftools_filter_cmd, shell=True)

        return f"{sample_id}_filtered.vcf.gz"

    def index_vcf(self, vcf_file: str):
        index_cmd = f"tabix {vcf_file}"
        subprocess.check_output(index_cmd, shell=True)

    def bcftools_consensus(self, vcf_file: str) -> str:
        sample_id = vcf_file.split("_filtered.vcf.gz")[0]
        bcftools_consensus_cmd = (
            f"bcftools consensus -f subset_ref.fa {vcf_file} > {sample_id}_gap_fill.fa"
        )
        subprocess.check_output(bcftools_consensus_cmd, shell=True)

        return f"{sample_id}_gap_fill.fa"

    def get_gap_length(self, gaps: dict) -> int:
        # -1 to not include the bases at either end of the blast hit region
        for k, v in gaps.items():
            gap_length = v - k - 1
        return gap_length

    def filter_consensus_seq(self, seq: str, gaps: list, iteration: int) -> str:
        for k, v in gaps[iteration].items():
            gap_length = v - k
        if len(seq) > gap_length:
            seq = seq[0 : gap_length - 1]
        return seq

    def fill_sequence_gap(
        self,
        seq_to_fill: str,
        gap_seq: str,
        gap_data: list,
        iteration: int,
        seq_added: int,
    ) -> str:
        filled_seq = str()
        for k, v in gap_data[iteration].items():
            # remove any trailing sequence if the sequence is longer than the gap
            ref_seq_start = seq_to_fill[0 : (k + seq_added)]
            ref_seq_end = seq_to_fill[(k + seq_added) :]
            filled_seq = f"{ref_seq_start}{gap_seq}{ref_seq_end}"
        return filled_seq
