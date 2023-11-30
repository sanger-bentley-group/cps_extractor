#!/usr/bin/env python3
import xml.etree.ElementTree as ET
import logging
import datetime

logging.basicConfig(
    filename=f"cps_extractor_{datetime.datetime.now().strftime('%Y-%m-%d %H:%M').replace(' ','-')}.log",
    encoding="utf-8",
    level=logging.INFO,
)


class Blast:
    def __init__(self, blast_results_file: str, hit_length: int):
        self.blast_results_file = blast_results_file
        self.hit_length = hit_length

    def parse_blast_results(self) -> list:
        blast_results = list()

        tree = ET.parse(self.blast_results_file)
        root = tree.getroot()

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
            if int(dict1["hit_start"]) < int(dict1["hit_end"]) and int(
                dict2["hit_start"]
            ) < int(dict2["hit_end"]):
                if int(dict2["hit_start"]) > int(dict1["hit_start"]) and int(
                    dict2["hit_end"]
                ) < int(dict1["hit_end"]):
                    overlap = True
        return overlap

    def get_best_hit(self, final_blast_results: list, serotype: str) -> str:
        max_seq_ratio = 0
        max_seq_length = 0
        best_serotype = str()
        best_serotype_results = list()

        # if the serotype is known, only select results for that serotype
        if serotype is not None:
            for result in final_blast_results:
                if serotype in result["hit_def"] and result["e_value"] < float(10 ** -50):
                    best_serotype_results.append(result)
            if len(best_serotype_results) == 0:
                logging.error(f"No results found for {serotype}, please check the blast results file")
                raise SystemExit(1)

        # get reference with best hit length to total length ratio, 'best guess' to determine serotype
        if serotype is None:
            for result in final_blast_results:
                seq_length = int(result["seq_length"])
                ref_length = int(result["reference_length"])
                seq_ratio = seq_length / ref_length
                print(f"{result['hit_def']} {seq_ratio}")
                if seq_ratio > max_seq_ratio:
                    max_seq_ratio = seq_ratio
                    best_serotype = result["hit_def"].strip().split()[0]
                    print(best_serotype)

            for result in final_blast_results:
                if best_serotype in result["hit_def"] and result["e_value"] < float(10**-50):
                    best_serotype_results.append(result)

        for result in best_serotype_results:
            seq_length = int(result["seq_length"])
            print(f"{seq_length} {result['hit_def']}")
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
        complement = {"A": "T", "T": "A", "C": "G", "G": "C", "-": "-"}
        reversed_sequence = sequence[::-1]
        reverse_complement_sequence = "".join(
            [complement[base] for base in reversed_sequence]
        )
        return reverse_complement_sequence

    def sort_and_reverse_complement_hits(self, blast_hits_dict: list) -> list:
        # sort hits based on start position, remove gaps and rev comp if needed
        for i in range(0, len(blast_hits_dict)):
            # rev comp if needed
            if blast_hits_dict[i]["hit_frame"] == -1:
                end = blast_hits_dict[i]["hit_start"]
                start = blast_hits_dict[i]["hit_end"]
                rc = self.reverse_complement(blast_hits_dict[i]["seq"])
                blast_hits_dict[i]["hit_start"] = start
                blast_hits_dict[i]["hit_end"] = end
                blast_hits_dict[i]["seq"] = rc

        # sort the data list on start pos
        sorted_data = sorted(blast_hits_dict, key=lambda x: x["hit_start"])
        return sorted_data

    def curate_sequence(self, sorted_data: list) -> str:
        # curate blast sequence from final blast results
        # if sequences don't overlap, join them together
        # if sequences do overlap, prioritise the match from the sequence with the higher length in overlapping regions
        logging.info(sorted_data)
        seq = str()
        if len(sorted_data) == 0:
            logging.error("No blast hits found, please check the blast XML file")
            raise SystemExit(1)
        elif len(sorted_data) == 1:
            seq = sorted_data[0]["seq"]
        else:
            for i in range(0, (len(sorted_data) - 1)):
                if i == 0:
                    if self.check_partial_overlap(sorted_data[i], sorted_data[i + 1]):
                        if int(sorted_data[i]["seq_length"]) > int(
                            sorted_data[i + 1]["seq_length"]
                        ):
                            seq += sorted_data[i]["seq"]
                            overlap_index = (
                                1
                                + int(sorted_data[i]["hit_end"])
                                - int(sorted_data[i + 1]["hit_start"])
                            )
                            seq += sorted_data[i + 1]["seq"][overlap_index::]
                        else:
                            seq += sorted_data[i]["seq"][
                                0 : int((sorted_data[i + 1]["hit_start"]) - 1)
                            ]
                            seq += sorted_data[i + 1]["seq"]
                    else:
                        seq += sorted_data[i]["seq"]
                        seq += sorted_data[i + 1]["seq"]
                else:
                    if self.check_partial_overlap(sorted_data[i], sorted_data[i + 1]):
                        if int(sorted_data[i]["seq_length"]) > int(
                            sorted_data[i + 1]["seq_length"]
                        ):
                            overlap_index = (
                                1
                                + int(sorted_data[i]["hit_end"])
                                - int(sorted_data[i + 1]["hit_start"])
                            )
                            seq += sorted_data[i + 1]["seq"][overlap_index::]
                        else:
                            overlap_index = (
                                1
                                + int(sorted_data[i]["hit_end"])
                                - int(sorted_data[i + 1]["hit_start"])
                            )
                            seq = seq[:-overlap_index]
                            seq += sorted_data[i + 1]["seq"]

                    else:
                        seq += sorted_data[i + 1]["seq"]

        # remove gaps
        seq = seq.replace("-", "")
        return seq

    def write_fasta(self, sequence: str, output_file: str):
        fasta_output = output_file.split(".")[0]
        if "/" in fasta_output:
            fasta_output = fasta_output.split("/")[-1]
        with open(output_file, "w") as fasta:
            fasta.write(f">{fasta_output}_cps\n")
            fasta.write(sequence)

    def parse_blast_results_dev(self) -> list:
        blast_results = list()

        tree = ET.parse(self.blast_results_file)
        root = tree.getroot()

        for query in root.findall(".//Iteration"):
            query_id = query.find(".//Iteration_query-def").text

            for hit in query.findall(".//Hit"):
                hit_length = int(hit.find(".//Hit_len").text)
                hit_def = hit.find(".//Hit_def").text
                if int(hit_length) >= self.hit_length:
                    for hsp in hit.findall(".//Hsp"):
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
                                "reference_length": hit_length,
                                "e_value": e_value,
                            }

                            blast_results.append(blast_result)
        logging.info(blast_results)
        return blast_results
