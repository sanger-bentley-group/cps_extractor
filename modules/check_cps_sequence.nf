// Check CPS sequence for disruptive mutations
process CHECK_CPS_SEQUENCE {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.csv"
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "**_cps.gff3", saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }

    label 'cps_extractor_python_container'
    label 'farm_low_fallible'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bakta_results), path(cps_sequence), val(reference)
    val(results_dir)

    output:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference), path(cps_sequence), emit: results_ch

    script:
    annotation_file="${bakta_results}/${sample_id}_cps.gff3"
    mutation_file="${sample_id}_cps_mutations.csv"
    """
    # catch the sanity check for sequence length and append the error message to the log file
    check_cps_sequence.py -c ${cps_sequence} -b ${bakta_results} || { cat cps_extractor.log >> ${results_dir}/${sample_id}/cps_extractor.log; exit 1; }
    """
}