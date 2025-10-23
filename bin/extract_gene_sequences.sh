#!/bin/bash

# extract the nucleotide sequence for CDS regions in the reference annotations
# create ARIBA metadata files for ARIBA database construction per serotype

sero=$1
ref_db=$2

bedtools_get_cds() {
  # extract CDS regions using bedtools and avoid genes such as transposases, so the ARIBA reports are more readable
  bedtools getfasta -fi ${ref_db}/fasta/${sero}.fasta \
  -bed <(sed '/FASTA/q' ${ref_db}/annotation/${sero}.gff \
   | sed "s/contig_1/${sero}/g" | grep CDS | grep -v -i intron \
   | grep -v -i lipoprotein | grep -v -i chlorohydrolase \
   | grep -v -i tnp | grep -v -i "family element" | grep -v -i transposase \
   | grep -v -i aliB | grep -v -i aliA | grep -v -i amiA | grep -v -i dexB \
   | grep -v -i lipoprotein | grep -v -i "family orf" \
   | grep -v -i "helix-turn-helix" | grep -v -i "hypothetical protein" \
   | grep -v -i  "aerotactic transducer" | grep -v -i "Mobile genetic element" \
   | grep -v -i "pseudo" | grep -w -v -f ${ref_db}/ariba_genes.txt) > ${sero}_cds.fna 
}


rename_genes() {
  # rename the genes to the actual gene names for use in the ARIBA report
  while read line; do
    if [[ $line == ">"* ]]; then
        start=$(echo ${line} | awk -F ":" '{ print $NF }' | awk -F "-" '{ print $1 }')
        # add 1 so bedtools and gff positions align
        start=$((start+1))
        anno=$(grep -w ${start} ${ref_db}/annotation/${sero}.gff)
        if [[ $anno == *"gene="* ]]; then
            anno_id=$(echo ${anno} | awk -F "gene=" '{ print $NF }')
        else
            anno_id=$(echo ${anno} | awk -F "ID=" '{ print $NF }' | awk -F ";" '{ print $1 }')
        fi
        echo ">${anno_id}" > ${sero}_${anno_id}.fa
        sed -n "/^${line}/,/^>/p" ${sero}_cds.fna | grep -v ">" >> ${sero}_${anno_id}.fa
    fi
  done < ${sero}_cds.fna
  cat ${sero}_*.fa > ${sero}_all.fasta
}

create_ariba_meta() {
  # create metadata file for ariba
  while read line
  do
    if [[ $line == ">"* ]]; then
      gene=$(echo ${line} | sed 's|>||g')
      echo -e "${gene}\t1\t0\t.\t.\t." >> ${sero}_ariba_meta.tsv
    fi
  done < ${sero}_all.fasta
}

bedtools_get_cds

rename_genes

create_ariba_meta