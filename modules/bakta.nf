// Annotate cps sequence using bakta
process BAKTA {

    label 'bakta_container'
    label 'farm_high_mem'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(cps_sequence), val(reference)
    path(prodigal_training_file)
    path(bakta_db)
    path(reference_database)

    output:
    tuple val(sample_id), path(bakta_results), path(cps_sequence), path(annotation_file), path(gb_file), val(reference), emit: bakta_results_ch

    script:
    bakta_results="${sample_id}_bakta"
    annotation_file="${bakta_results}/${sample_id}_cps.gff3"
    gb_file="${bakta_results}/${sample_id}_cps.gbff"
    """
    bakta --db ${bakta_db} -t ${params.bakta_threads} -o ${sample_id}_bakta --prodigal-tf ${prodigal_training_file} --proteins ${reference_database}/proteins/${reference}_proteins.txt --skip-plot ${cps_sequence}
    """
}