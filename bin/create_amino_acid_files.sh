#!/bin/bash

# script to create files containing amino acid sequences for each gene in the CPS region from gffread output file

proteins_file=$1
sample=$2

# get unique ids
ids=$(grep ">" ${proteins_file} | awk '{ print $1 }')

# get all genes in the sample
get_genes() {
  for id in ${ids}
  do
    identifier=$(grep -w -m 1 ${id} ${proteins_file})
    if [[ "${identifier}" == *";gene="* ]]
    then
      genes+=$(echo ${identifier} | awk -F ";gene=" '{ print $NF }' | sed 's|$| |g')
    fi
  done
  echo ${genes}
}

# get duplicate genes in the sample
get_duplicate_genes() {
for gene in ${gene_list}
  do 
    occurences=$(echo ${gene_list} | grep -w -o ${gene} | wc -l)
    if [ "${occurences}" -gt 1 ]
    then
      match=$(echo ${duplicate_list} | grep -w -o ${gene} | wc -l)
      if [ "$match" -eq 0 ]
      then
        gene=$(echo ${gene} | sed 's|$| |g')
        duplicate_list+=${gene}
      fi
    fi
  done
  echo ${duplicate_list}
}

make_protein_files() {
  for id in ${ids}
  do
    aa_seq=$(sed -n -e "/${id}/,/>/ p" ${proteins_file} | grep -v ">" | tr -d '\n')
    identifier=$(grep -w -m 1 ${id} ${proteins_file})
    # if a gene name is in the annotation, use it to name the file
    # otherwise name using the unique ID which can be looked up in the annotation
    if [[ "${identifier}" == *";gene="* ]]
    then
      gene=$(echo ${identifier} | awk -F ";gene=" '{ print $NF }')
      match=$(echo ${duplicate_gene_list} | grep -w -o ${gene} | wc -l)
      # if the gene only occurs once, write the file with just the gene name
      # otherwise, write the gene name and id in the filename for duplicate genes
      if [ "$match" -eq 0 ]
      then
        echo ${identifier} > ${sample}-${gene}_protein.fa
        echo ${aa_seq} >> ${sample}-${gene}_protein.fa
      else
        id_str=$(echo ${id} | sed 's|>||g')
        echo ${identifier} > ${sample}-${gene}_${id_str}_protein.fa
        echo ${aa_seq} >> ${sample}-${gene}_${id_str}_protein.fa
      fi
    else
      id_str=$(echo ${id} | sed 's|>||g')
      echo ${identifier} > ${sample}-${id_str}_protein.fa
      echo ${aa_seq} >> ${sample}-${id_str}_protein.fa
    fi
  done
}

gene_list=$(get_genes)

duplicate_gene_list=$(get_duplicate_genes)

make_protein_files
