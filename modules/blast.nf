// Return blast results
process BLASTN {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*_blast_results.xml"

    label 'blast_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(assembly), val(sero), path(reads)
    val(blast_db)

    output:
    tuple val(sample_id), path(assembly), path(blast_results), val(sero), path(reads), emit: blast_results_ch

    script:
    blast_results="${sample_id}_blast_results.xml"
    """
    blast_db_full=\$(echo $blast_db | awk -F "." '{ print \$1 }')
    blastn -query ${assembly} -db \${blast_db_full} -out ${sample_id}_blast_results.xml -outfmt 5
    """
}