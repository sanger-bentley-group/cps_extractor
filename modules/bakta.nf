// Annotate cps sequence using bakta
process BAKTA {

    label 'bakta_container'
    label 'farm_high_mem'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(cps_sequence), val(reference)
    path(prodigal_training_file)
    path(bakta_db)

    output:
    tuple val(sample_id), path(bakta_results), path(cps_sequence), val(reference), emit: bakta_results_ch

    script:
    bakta_results="${sample_id}_bakta"
    """
    bakta --db ${bakta_db} -t ${params.bakta_threads} -o ${sample_id}_bakta --prodigal-tf ${prodigal_training_file} --skip-plot ${cps_sequence}
    """
}