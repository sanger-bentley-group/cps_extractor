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
    tuple val(sample_id), path(bakta_results), path(cps), path(annotation_file), path(gb_file), val(reference), emit: bakta_results_ch

    script:
    cps="${sample_id}_cps.fa"
    bakta_results="${sample_id}_bakta"
    annotation_file="${bakta_results}/${sample_id}_cps.gff3"
    gb_file="${bakta_results}/${sample_id}_cps.gbff"
    """
    # rename cps with copy to avoid caching issues
    cp ${cps_sequence} ${sample_id}_cps.fa
    bakta --db ${bakta_db} -t "`nproc`" -o ${sample_id}_bakta --prodigal-tf ${prodigal_training_file} --proteins ${reference_database}/proteins/${reference}_proteins.txt --skip-plot ${sample_id}_cps.fa
    """
}