// Return blast results
process BLASTN {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*_blast_results.xml"

    label 'blast_container'
    label 'farm_low'
    label 'farm_scratchless'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(assembly), val(sero), path(reads)
    val(blast_db)
    val(cps_reference_db)

    output:
    tuple val(sample_id), path(assembly), path(blast_results), path(dexb_results), path(alia_results), val(sero), path(reads), emit: blast_results_ch

    script:
    blast_results="${sample_id}_blast_results.xml"
    dexb_results="${sample_id}_dexb_blast_results.xml"
    alia_results="${sample_id}_aliA_blast_results.xml"
    """
    blast_db_full=\$(echo $blast_db | awk -F "." '{ print \$1 }')
    blastn -query ${assembly} -db \${blast_db_full} -out ${sample_id}_blast_results.xml -outfmt 5
    blastn -query ${assembly} -db ${cps_reference_db}/dexbdb -out ${sample_id}_dexb_blast_results.xml -outfmt 5
    blastn -query ${assembly} -db ${cps_reference_db}/aliadb -out ${sample_id}_aliA_blast_results.xml -outfmt 5
    """
}