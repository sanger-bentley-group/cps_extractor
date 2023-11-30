// Curate CPS sequence
process CURATE_CPS_SEQUENCE {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*{_cps.fa,log}"

    label 'cps_extractor_python'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(assembly), path(blast_results)

    output:
    tuple val(sample_id), path(cps_sequence), emit: cps_sequence_ch
    path("*.log"), emit: log_ch

    script:
    cps_sequence="${sample_id}_cps.fa"
    """
    if [ "${params.serotype}" == "" ]
    then
        curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa
    else
        curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa -s ${params.serotype}
    fi
    """
}