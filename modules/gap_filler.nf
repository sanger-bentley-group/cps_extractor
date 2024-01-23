process GAP_FILLER {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*_cps.fa"

    label 'gap_filler'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(cps_sequence), path(reads), val(serotype)
    path log_file
    path reference_database

    output:
    tuple val(sample_id), path(cps_sequence), emit: gap_filled_ch

    script:
    read1="${reads[0]}"
    read2="${reads[1]}"
    """
    reference=\$(grep -i -v error $log_file | head -1 | awk -F ":" '{ print \$9 }' | awk -F "," '{ print \$1 }' | sed "s|'||g" | sed 's| ||g')
    grep -i -v error $log_file | head -1 > blast_results.log
    fill_sequence_gaps.py -l blast_results.log -a ${reference_database}/annotation/\${reference}.gff \
    -r ${reference_database}/fasta/\${reference}.fasta -r1 $read1 -r2 $read2 -i $cps_sequence

    if [ -e "gap_filled_seq.fa" ]
    then
        tail -n +2 gap_filled_seq.fa > tmp.fa
        echo -e ">${sample_id}_cps" > ${cps_sequence}
        cat tmp.fa >> ${cps_sequence}
    fi
    """
}