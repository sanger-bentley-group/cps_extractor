// Use panaroo to assess gene content difference vs reference
process PANAROO_REF_COMPARISON {

    label 'panaroo_container'
    label 'farm_low_fallible'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference)
    path reference_database

    output:
    tuple val(sample_id), path(panaroo_results), val(reference), path(annotation_file), emit: panaroo_results_ch
    tuple val(sample_id), path(gene_alignment_results), emit: gene_alignment_results

    script:
    panaroo_results="${sample_id}_panaroo_results"
    gene_alignment_results="${sample_id}_panaroo_results/aligned_gene_sequences/*.fas"
    """
    panaroo -i ${annotation_file} ${reference_database}/annotation/${reference}.gff -o ${sample_id}_panaroo_results --clean-mode strict -a core
    """
}

// Check the gene order of the isolate vs the reference
process CHECK_GENE_ORDER {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.csv"
    label 'check_gene_content_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), path(mutation_file), val(reference)
    path reference_database

    output:
    tuple val(sample_id), path(gene_comparison), emit: gene_comparison_ch
    script:
    gene_comparison="${sample_id}_gene_comparison.csv"
    """
    compare_gene_order.py -i ${annotation_file} -r ${reference_database}/annotation/${reference}.gff -o ${sample_id}_gene_comparison.csv
    """
}

// Check the gene order of the isolate vs the reference
process SNP_DISTS {
    publishDir "${params.output}/${sample_id}/snp_dists", mode: 'copy', overwrite: true, pattern: "*.csv"
    label 'snp_dists_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(gene_alignment_file)

    output:
    path("*.csv"), emit: snp_dists_ch
    script:
    """
    gene_name=\$(echo ${gene_alignment_file} | awk -F ".aln.fas" '{ print \$1 }')
    snp-dists -c ${gene_alignment_file} > \${gene_name}_snp_dists.csv
    """
}