#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
import difflib
from Bio import SeqIO
from Bio.Seq import Seq

logging.basicConfig(
    filename=f"cps_extractor.log",
    encoding="utf-8",
    level=logging.INFO,
)


class Blast:
    def __init__(self, blast_results_file: str, assembly_file: str, hit_length: int):
        self.blast_results_file = blast_results_file
        self.assembly_file = assembly_file
        self.hit_length = hit_length

    def parse_blast_results(self, results_file: str) -> list:
        blast_results = list()

        tree = ET.parse(results_file)
        root = tree.getroot()

        if results_file != self.blast_results_file:
            self.hit_length = 700

        for query in root.findall(".//Iteration"):
            query_id = query.find(".//Iteration_query-def").text

            for hit in query.findall(".//Hit"):
                hit_length = int(hit.find(".//Hit_len").text)
                hit_def = hit.find(".//Hit_def").text
                if int(hit_length) >= self.hit_length:
                    for hsp in hit.findall(".//Hsp"):
                        q_sequence = hsp.find(".//Hsp_qseq").text
                        e_value = float(hsp.find(".//Hsp_evalue").text)
                        hit_start = int(hsp.find(".//Hsp_hit-from").text)
                        hit_end = int(hsp.find(".//Hsp_hit-to").text)
                        aln_len = hsp.find(".//Hsp_align-len").text
                        if hit_start < hit_end:
                            seq_length = int(hit_end) - int(hit_start) + 1
                        else:
                            seq_length = int(hit_start) - int(hit_end) + 1
                        hit_frame = hsp.find(".//Hsp_hit-frame").text
                        if int(aln_len) >= self.hit_length:
                            blast_result = {
                                "hit_start": hit_start,
                                "hit_end": hit_end,
                                "hit_frame": int(hit_frame),
                                "seq_length": seq_length,
                                "query_id": query_id,
                                "hit_def": hit_def,
                                "seq": q_sequence,
                                "reference_length": hit_length,
                                "e_value": e_value,
                            }

                            blast_results.append(blast_result)
        return blast_results

    def do_dicts_overlap(self, dict1: dict, dict2: dict) -> bool:
        # check if one blast sequence is entirely contained in a larger hit
        overlap = False

        if int(dict1["seq_length"]) > int(dict2["seq_length"]):
            if int(dict1["hit_start"]) <= int(dict1["hit_end"]) and int(
                dict2["hit_start"]
            ) <= int(dict2["hit_end"]):
                if int(dict2["hit_start"]) >= int(dict1["hit_start"]) and int(
                    dict2["hit_end"]
                ) <= int(dict1["hit_end"]):
                    overlap = True
        return overlap

    def get_best_hit(self, final_blast_results: list, serotype: str) -> str:
        max_seq_ratio = 0
        max_seq_length = 0
        best_serotype = str()
        best_serotype_results = list()

        # basic check to see if there are any blast hits before running the code
        if len(final_blast_results) == 0:
            logging.error(
                "The blast results contain no hits, please check your blast database and input sequence. Your input sequence may be a non encapsulated strain."
            )
            raise SystemExit(1)

        # if the serotype is known, only select results for that serotype
        if serotype is not None:
            for result in final_blast_results:
                if serotype.lower() in result["hit_def"].lower() and result[
                    "e_value"
                ] < float(10**-50):
                    best_serotype_results.append(result)
            if len(best_serotype_results) == 0:
                logging.error(
                    f"No results found for {serotype}, please check the blast results file"
                )
                raise SystemExit(1)

        # get reference with best hit length to total length ratio, 'best guess' to determine serotype
        if serotype is None:
            for result in final_blast_results:
                seq_length = int(result["seq_length"])
                ref_length = int(result["reference_length"])
                seq_ratio = seq_length / ref_length
                if seq_ratio > max_seq_ratio:
                    max_seq_ratio = seq_ratio
                    best_serotype = result["hit_def"].strip().split()[0]

            for result in final_blast_results:
                if best_serotype in result["hit_def"] and result["e_value"] < float(
                    10**-1
                ):
                    best_serotype_results.append(result)

        for result in best_serotype_results:
            seq_length = int(result["seq_length"])
            if seq_length > max_seq_length:
                max_seq_length = seq_length
                max_hit_def = result["hit_def"]
        return max_hit_def

    def compare_blast_dicts(self, blast_results: list, serotype: str) -> list:
        overlaps = list()

        # get the longest hit, remove hits from other references in the blastdb
        largest_hit_ref = self.get_best_hit(blast_results, serotype)
        blast_results = [
            item for item in blast_results if item["hit_def"] == largest_hit_ref
        ]
        # sort hits based on start position, and reverse positions of reverse hits for sorting
        for i in range(0, len(blast_results)):
            if blast_results[i]["hit_frame"] == -1:
                end = blast_results[i]["hit_start"]
                start = blast_results[i]["hit_end"]
                blast_results[i]["hit_start"] = start
                blast_results[i]["hit_end"] = end
        blast_results = sorted(blast_results, key=lambda x: x["hit_start"])

        # check for any sequences that are contained entirely by larger hits and remove them
        for i in range(len(blast_results)):
            for j in range(i + 1, len(blast_results)):
                dict1 = blast_results[i]
                dict2 = blast_results[j]
                if self.do_dicts_overlap(dict1, dict2):
                    if int(dict1["seq_length"]) > int(dict2["seq_length"]):
                        overlaps.append(dict2)
                    else:
                        overlaps.append(dict1)

        final_results = [item for item in blast_results if item not in overlaps]
        return final_results

    def check_partial_overlap(self, dict1: dict, dict2: dict) -> bool:
        # check if one sequence partially overlaps another sequence
        overlap = False
        if int(dict1["hit_end"]) >= int(dict2["hit_start"]):
            overlap = True
        return overlap

    def reverse_complement(self, sequence: str) -> str:
        # rev comp, keep gaps as they are
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-", "N": "N"}
        reversed_sequence = sequence[::-1]
        reverse_complement_sequence = "".join(
            [complement[base] for base in reversed_sequence]
        )
        return reverse_complement_sequence

    def reverse_complement_hits(self, blast_hits_dict: list) -> list:
        # reverse complement hits if needed
        for i in range(0, len(blast_hits_dict)):
            if blast_hits_dict[i]["hit_frame"] == -1:
                rc = self.reverse_complement(blast_hits_dict[i]["seq"])
                blast_hits_dict[i]["seq"] = rc

        return blast_hits_dict

    def join_overlap_sequences(self, s1: str, s2: str) -> str:
        # blast hit positions are not always accurate, check if the sequences overlap with difflib
        s1_end = str(s1[(len(s1) - 200) :])
        s2_start = str(s2[:200])

        # find the best sequence match between sequences
        s = difflib.SequenceMatcher(None, s1_end, s2_start, autojunk=False)
        pos_a, pos_b, size = s.find_longest_match(0, len(s1_end), 0, len(s2_start))
        best_overlap = s1_end[pos_a : pos_a + size]

        # if the best match is at the start of the second sequence, remove it so not to duplicate the sequence
        if best_overlap == s2[: len(best_overlap)] and len(best_overlap) >= 5:
            s2 = s2[len(best_overlap) :]

        # join the sequences
        joined_sequence = s1 + s2
        return joined_sequence

    def are_hits_same_contig(self, blast_results: list) -> bool:
        results_length = len(blast_results)
        if results_length == 2:
            if blast_results[0]["query_id"] == blast_results[1]["query_id"]:
                return True
        if results_length == 3:
            if (
                blast_results[0]["query_id"] == blast_results[1]["query_id"]
                and blast_results[1]["query_id"] == blast_results[2]["query_id"]
            ):
                return True
        return False

    def join_sequence(
        self, start_seq: str, end_seq: str, query_id: str, hit_orientation: int
    ) -> str:
        start_seq = start_seq.replace("-", "")
        end_seq = end_seq.replace("-", "")
        # join up  sequences using the assembly (gene insertions come up as multiple blast hits to same contig or dexB-aliA joins)
        for record in SeqIO.parse(self.assembly_file, "fasta"):
            if hit_orientation == 1:
                seq = str(record.seq)
                start_index = seq.find(start_seq[:100])
                end_index = seq.find(
                    end_seq[(len(end_seq) - 100) :], start_index + len(start_seq)
                )
                end_index += 100
                if record.name == query_id.split(" ")[0]:
                    break
            elif hit_orientation == -1:
                seq = str(record.seq)
                start_seq_rc = str(Seq(start_seq).reverse_complement())
                end_seq_rc = str(Seq(end_seq).reverse_complement())
                start_index = seq.find(end_seq_rc[:100])
                end_index = seq.find(
                    start_seq_rc[(len(start_seq_rc) - 100) :],
                    start_index + len(end_seq_rc),
                )
                end_index += 100
                if record.name == query_id.split(" ")[0]:
                    break
            else:
                return str()

        if start_index != -1 and end_index != -1:
            if hit_orientation == 1:
                extracted_region = seq[start_index:end_index]
                # create an empty file to show gap is filled
                with open("gap_filled", "w") as fp:
                    pass
                return extracted_region
            else:
                extracted_region = seq[start_index:end_index]
                # create an empty file to show gap is filled
                with open("gap_filled", "w") as fp:
                    pass
                return str(Seq(extracted_region).reverse_complement())

        return str()

    def curate_sequence(
        self, sorted_data: list, dexb_data: list, alia_data: list
    ) -> str:
        # curate blast sequence from final blast results
        # if there are more than 3 separate blast hits, the CPS sequence will not be constructed due to data quality issues in the assembly
        logging.info(sorted_data)
        seq = str()
        seq_1 = str()
        alia_seq = str()
        dexb_seq = str()

        if len(sorted_data) > 3:
            logging.error(
                f"There are a large number of blast hits for the CPS region ({len(sorted_data)} hits),\
            please check the quality of your input data"
            )
            raise SystemExit(1)
        if len(sorted_data) == 0:
            logging.error(
                "No blast hits were found for the CPS region, please check the blast results file for more information.\
                 You may have a non encapsulated strain of S.pneumoniae"
            )
            raise SystemExit(1)
        # join sequence between dexB and aliA if they are on the same contig
        elif (
            len(alia_data) >= 1
            and len(dexb_data) >= 1
            and alia_data[0]["query_id"] == dexb_data[0]["query_id"]
        ):

            seq = self.join_sequence(
                dexb_data[0]["seq"],
                alia_data[0]["seq"],
                dexb_data[0]["query_id"],
                dexb_data[0]["hit_frame"],
            )

        # join cps to aliA if in the same contig
        elif (
            len(alia_data) >= 1
            and len(sorted_data) == 1
            and alia_data[0]["query_id"] == sorted_data[0]["query_id"]
        ):

            seq = self.join_sequence(
                sorted_data[0]["seq"],
                alia_data[0]["seq"],
                alia_data[0]["query_id"],
                alia_data[0]["hit_frame"],
            )
            logging.info(
                "Warning: aliA and dexB are fragmented across multiple contigs, please note that this may be a partial CPS sequence"
            )

        # join cps to dexB if in the same contig
        elif (
            len(dexb_data) >= 1
            and len(sorted_data) == 1
            and dexb_data[0]["query_id"] == sorted_data[0]["query_id"]
        ):

            seq = self.join_sequence(
                dexb_data[0]["seq"],
                sorted_data[0]["seq"],
                dexb_data[0]["query_id"],
                dexb_data[0]["hit_frame"],
            )
            logging.info(
                "Warning: aliA and dexB are fragmented across multiple contigs, please note that this may be a partial CPS sequence"
            )

        elif len(sorted_data) == 1:
            seq = sorted_data[0]["seq"]
        elif len(sorted_data) == 2:
            if self.are_hits_same_contig(sorted_data):
                seq = self.join_sequence(
                    sorted_data[0]["seq"],
                    sorted_data[1]["seq"],
                    sorted_data[0]["query_id"],
                    sorted_data[0]["hit_frame"],
                )
                if (
                    len(dexb_data) >= 1
                    and dexb_data[0]["query_id"] == sorted_data[0]["query_id"]
                ):
                    seq = self.join_sequence(
                        dexb_data[0]["seq"],
                        seq,
                        dexb_data[0]["query_id"],
                        dexb_data[0]["hit_frame"],
                    )
                if (
                    len(alia_data) >= 1
                    and alia_data[0]["query_id"] == sorted_data[0]["query_id"]
                ):
                    seq = self.join_sequence(
                        seq,
                        alia_data[0]["seq"],
                        alia_data[0]["query_id"],
                        alia_data[0]["hit_frame"],
                    )
            if not seq:
                if (
                    len(dexb_data) >= 1
                    and dexb_data[0]["query_id"] == sorted_data[0]["query_id"]
                ):
                    dexb_seq = self.join_sequence(
                        dexb_data[0]["seq"],
                        sorted_data[0]["seq"],
                        dexb_data[0]["query_id"],
                        dexb_data[0]["hit_frame"],
                    )
                if (
                    len(alia_data) >= 1
                    and alia_data[0]["query_id"] == sorted_data[1]["query_id"]
                ):
                    alia_seq = self.join_sequence(
                        sorted_data[1]["seq"],
                        alia_data[0]["seq"],
                        alia_data[0]["query_id"],
                        alia_data[0]["hit_frame"],
                    )
                if dexb_seq and alia_seq:
                    seq = self.join_overlap_sequences(dexb_seq, alia_seq)
                elif dexb_seq and not alia_seq:
                    seq = self.join_overlap_sequences(dexb_seq, sorted_data[1]["seq"])
                elif not dexb_seq and alia_seq:
                    seq = self.join_overlap_sequences(sorted_data[0]["seq"], alia_seq)
                else:
                    seq = self.join_overlap_sequences(
                        sorted_data[0]["seq"], sorted_data[1]["seq"]
                    )
            if not dexb_seq or not alia_seq:
                logging.info(
                    "Warning: aliA and dexB are fragmented across multiple contigs, please note that this may be a partial CPS sequence"
                )
            logging.info(
                "Warning: The CPS sequence for this sample is fragmented across 2 blast hits - there may be a data quality issue"
            )
        elif len(sorted_data) == 3:
            if not self.are_hits_same_contig(
                sorted_data[:2]
            ) or not self.are_hits_same_contig(sorted_data[1:3]):
                logging.error(
                    "There are 3 blast hits split across multiple contigs, there is a data quality issue."
                )
                raise SystemExit(1)
            seq_1 = self.join_sequence(
                sorted_data[0]["seq"],
                sorted_data[1]["seq"],
                sorted_data[0]["query_id"],
                sorted_data[0]["hit_frame"],
            )

            seq = self.join_sequence(
                seq_1,
                sorted_data[2]["seq"],
                sorted_data[2]["query_id"],
                sorted_data[2]["hit_frame"],
            )
            logging.info(
                "Warning: The CPS sequence for this sample is fragmented across 3 blast hits - there may be a data quality issue"
            )
        else:
            logging.error(
                f"There are a large number of blast hits for the CPS region ({len(sorted_data)} hits),\
            please check the quality of your input data"
            )
            raise SystemExit(1)

        # remove gaps and Ns
        seq = seq.replace("-", "")
        seq = seq.replace("N", "")
        return seq

    def write_fasta(self, sequence: str, output_file: str):
        fasta_output = output_file.split(".")[0]
        if "/" in fasta_output:
            fasta_output = fasta_output.split("/")[-1]
        with open(output_file, "w") as fasta:
            fasta.write(f">{fasta_output}\n")
            fasta.write(sequence)

    def parse_blast_results_dev(self, blast_results) -> list:
        # function to pass the blast results to a log file without the sequences for readablility
        for result in blast_results:
            del result["seq"]
        logging.info(blast_results)
        return blast_results
