// Check CPS sequence for disruptive mutations
process CHECK_CPS_SEQUENCE {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.{csv,gff3}"

    label 'cps_extractor_python'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bakta_results), path(cps_sequence)

    output:
    tuple val(sample_id), path(annotation_file), path(mutation_file), emit: results_ch

    script:
    annotation_file="${sample_id}_cps.gff3"
    mutation_file="${sample_id}_cps_mutations.csv"
    """
    check_cps_sequence.py -c ${cps_sequence} -b ${bakta_results}
    mv ${bakta_results}/${sample_id}_cps.gff3 .
    """
}