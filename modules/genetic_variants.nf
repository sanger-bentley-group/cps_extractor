process GET_GENETIC_VARIANTS {
    label 'check_gene_content_container'
    label 'farm_low'

    input:
    path(annotations)
    val(serotype)
    path(cps_reference_database)

    output:
    path("unique_genetic_variants.txt"), emit: genetic_variants_ch, optional: true
    path("*_genes.txt"), emit: genes_ch, optional: true

    script:
    """
    grep -rli "dexb" *_cps.gff3 | while read -r file; do
    if grep -qi "alia" "\$file" || grep -qi "oppa" "\$file"; then
    echo "\$file" >> full_cps.txt
    fi
    done
    if [[ -s "full_cps.txt" ]]
    then
        while read file
        do
        sample=\$(echo \$file | awk -F "_cps.gff3" '{ print \$1 }')
        variants.py -a \$file > \${sample}_genes.txt
        done < full_cps.txt
        md5sum *_genes.txt | sort -k1,1 -u | awk '{ print \$NF }' | sed 's|_genes.txt||g' > unique_genetic_variants.txt
    fi
    """
}

process ASSIGN_VARIANT_GROUPS {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "*s.csv"
    label 'check_gene_content_container'
    label 'farm_low'

    input:
    path(gene_lists)

    output:
    tuple path("genetic_groups.csv"), path("genetic_clusters.csv"), emit: genetic_groups_ch

    script:
    """
    echo md5,sample > md5.csv
    md5sum *_genes.txt | sed 's|  *|,|g' | sed 's|_genes.txt||g' >> md5.csv
    genetic_groups.py -g md5.csv
    tail -n +2 genetic_groups.csv | sort -t , -k2 -u | while read line
    do
        sample=\$(echo \${line} | awk -F "," '{ print \$1 }')
        group=\$(echo \${line} | awk -F "," '{ print \$2 }')
        cat \${sample}_genes.txt | tr '\n' '-' | rev | sed 's|-||' | rev > \${group}_gene_cluster.txt
    done
    
    for f in *_gene_cluster.txt
    do
      group=\$(echo \$f | awk -F "_gene_cluster.txt" '{ print \$1 }')
      gene_cluster=\$(cat \$f)
      echo -e "\${group},\${gene_cluster}" >> genetic_clusters.csv
    done
      
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
        sample_list=\$(cat ${genetic_variant_list} | sed 's|\$|_cps.gbff|g')
        echo \${sample_list}
        clinker ${reference_database}/genbank/${reference}.gb \${sample_list} -p genetic_variants_plot.html -gf ${reference_database}/clinker_descriptions/${reference}_gene_info.csv -cm ${reference_database}/clinker_descriptions/colours.csv
    fi
    """
}
