// Use panaroo to assess gene content difference vs reference
process PANAROO_REF_COMPARISON {

    label 'panaroo_container'
    label 'farm_low_fallible'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), val(reference)
    path reference_database

    output:
    tuple val(sample_id), path(panaroo_results), val(reference), path(annotation_file), emit: panaroo_results_ch

    script:
    panaroo_results="${sample_id}_panaroo_results"
    """
    panaroo -i ${annotation_file} ${reference_database}/annotation/${reference}.gff -o ${sample_id}_panaroo_results --clean-mode strict -a core
    """
}

process PANAROO_ALL {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "panaroo_pangenome_results"

    label 'panaroo_container'
    label 'farm_high'

    input:
    path annotations
    path reference_database
    val serotype

    output:
    path(panaroo_results)

    script:
    panaroo_results="panaroo_pangenome_results"
    """
    cp ${reference_database}/annotation/${serotype}.gff ${serotype}.gff3
    panaroo -i *.gff3 -o panaroo_pangenome_results --clean-mode strict -a pan --threads 32
    """
}

// Check the gene order of the isolate vs the reference
process CHECK_GENE_ORDER {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.tsv"
    label 'check_gene_content_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(annotation_file), val(reference)
    path reference_database

    output:
    tuple val(sample_id), path(gene_order), emit: gene_order_ch
    script:
    gene_order="${sample_id}_gene_order.tsv"
    """
    variants.py -a ${annotation_file} > ${sample_id}_order.txt
    variants.py -a ${reference_database}/annotation/${reference}.gff > ${reference}_order.txt
    echo -e "${sample_id}_genes\t${reference}_reference_genes" > ${sample_id}_gene_order.tsv
    paste ${sample_id}_order.txt ${reference}_order.txt >> ${sample_id}_gene_order.tsv
    """
}

// Create a plot of gene alignments using clinker
process CLINKER {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.html"
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "**_cps.gff3", saveAs: { fn -> fn.substring(fn.lastIndexOf('/')+1) }
    label 'clinker_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(bakta_results), path(cps_sequence), path(annotation_file), path(gb_file), val(reference)
    path reference_database

    output:
    path("*.html")
    path(annotation_file)

    script:
    """
    clinker ${reference_database}/genbank/${reference}.gb ${gb_file} -p ${sample_id}_plot.html -gf ${reference_database}/clinker_descriptions/${reference}_gene_info.csv
    """
}

// Create a plot of gene alignments using clinker
process FIND_POTENTIAL_NEW_SEROTYPES {
    label 'check_gene_content_container'
    label 'farm_low'

    input:
    path(disrupted_genes_file)
    path(reference_database)

    output:
    path("potential_new_serotypes.txt"), emit: new_sero_ch

    script:
    """
    check_new_serotype.py -s ${params.serotype} -k ${reference_database}/disrupted_genes.csv -d ${disrupted_genes_file} > potential_new_serotypes.txt
    """
}

// Create a plot of gene alignments using clinker
process CLINKER_NEW_SEROTYPES {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "*.html"
    label 'clinker_container'
    label 'farm_high_slow'

    input:
    path(gb_files)
    path(reference_database)
    val(reference)
    path(disrupted_sample_list)

    output:
    path("*.html"), optional: true
    script:
    """
    if [[ -s "${disrupted_sample_list}" ]]
    then
        sample_list=\$(cat potential_new_serotypes.txt | sed 's|\$|_cps.gbff|g')
        echo \${sample_list}
        clinker ${reference_database}/genbank/${reference}.gb \${sample_list} -p potential_new_serotypes_plot.html -gf ${reference_database}/clinker_descriptions/${reference}_gene_info.csv
    fi
    """
}