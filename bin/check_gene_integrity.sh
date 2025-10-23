#!/bin/bash

ariba_genes_file=$1
ariba_report_file=$2

assembled_genes=$(zcat ${ariba_genes_file} | grep ">")

echo "gene_name,gene_integrity"

for assembled_gene in ${assembled_genes}
do
  id=$(echo ${assembled_gene} | awk -F "." '{ print $1"."$2"."$3 }' | sed 's|>||g')
  gene=$(grep -m1 ${id} ${ariba_report_file} | awk -F "\t" '{ print $1 }')
  if [[ "${assembled_gene}" == *"HAS_STOP"* ]]
  then
    echo "${gene},disrupted"
  else
    echo "${gene},intact"
  fi  
done
