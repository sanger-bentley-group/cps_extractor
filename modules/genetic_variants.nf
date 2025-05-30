process GET_GENETIC_VARIANTS {
    label 'check_gene_content_container'
    label 'farm_low'

    input:
    path(annotations)
    val(serotype)
    path(cps_reference_database)

    output:
    path("genetic_variants.txt"), emit: genetic_variants_ch

    script:
    """
    grep -rli "dexb" *_cps.gff3 | while read -r file; do
    if grep -qi "alia" "\$file" || grep -qi "oppa" "\$file"; then
    echo "\$file" >> full_cps.txt
    fi
    done
    
    while read file
    do
      sample=\$(echo \$file | awk -F "_cps.gff3" '{ print \$1 }')
      variants.py -a \$file > \${sample}_genes.txt
    done < full_cps.txt
    md5sum *_genes.txt | sort -k1,1 -u | awk '{ print \$NF }' | sed 's|_genes.txt||g' > genetic_variants.txt
    """
}

// Create a plot of gene alignments using clinker
process CLINKER_GENETIC_VARIANTS {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "*.html"
    label 'clinker_container'
    label 'farm_high_slow'

    input:
    path(gb_files)
    path(reference_database)
    val(reference)
    path(genetic_variant_list)

    output:
    path("*.html"), optional: true
    script:
    """
    if [[ -s "${genetic_variant_list}" ]]
    then
        sample_list=\$(cat genetic_variants.txt | sed 's|\$|_cps.gbff|g')
        echo \${sample_list}
        clinker ${reference_database}/genbank/${reference}.gb \${sample_list} -p genetic_variants_plot.html -gf ${reference_database}/clinker_descriptions/${reference}_gene_info.csv
    fi
    """
}