# rename panaroo .aln files with an ID that can be looked up in the annotation if they have no gene name

for f in ${ALIGNMENT_FOLDER}/group_*
do
  marker=$(grep -w ${SAMPLE}_cps ${f} | awk -F ";" '{ print $NF }')
  # account for edge case where panaroo has found extra genes that are in the reference only
  # in this case the annotation ID from the reference is in the name rather than the sample
  if [[ ${marker} == *"refound"* ]]
  then
    ref_marker=$(head -1 ${f} | awk -F ";" '{ print $NF }')
    id=$(grep -w ${ref_marker} ${GENE_DATA_FILE} | awk -F "," '{ print $4 }')
    mv ${f} ${ALIGNMENT_FOLDER}/${id}_refound.aln.fas
  else
    id=$(grep -w ${marker} ${GENE_DATA_FILE} | awk -F "," '{ print $4 }')
    mv ${f} ${ALIGNMENT_FOLDER}/${id}.aln.fas
  fi
done