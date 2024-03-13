// Use panaroo to assess gene content difference vs reference
process PANAROO_REF_COMPARISON {

    label 'farm_scratchless'
    label 'panaroo_container'
    label 'farm_low_fallible'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference)
    path reference_database

    output:
    tuple val(sample_id), path(panaroo_results), val(reference), path(annotation_file), emit: panaroo_results_ch

    script:
    panaroo_results="${sample_id}_panaroo_results"
    """
    panaroo -i ${annotation_file} ${reference_database}/annotation/${reference}.gff -o ${sample_id}_panaroo_results --clean-mode strict -a core
    """
}

// Check the gene order of the isolate vs the reference
process CHECK_GENE_ORDER {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.csv"
    label 'farm_scratchless'
    label 'check_gene_content_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(panaroo_results), val(reference), path(annotation_file)
    path reference_database

    output:
    tuple val(sample_id), path(gene_comparison), emit: gene_comparison_ch
    script:
    gene_comparison="${sample_id}_gene_comparison.csv"
    """
    compare_gene_order.py -i ${annotation_file} -r ${reference_database}/annotation/${reference}.gff -o ${sample_id}_gene_comparison.csv
    """
}