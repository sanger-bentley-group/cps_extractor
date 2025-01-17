// Run ARIBA to check for mutations in key cps genes
process ARIBA {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: '*/report.tsv', saveAs: { filename -> "ariba_report.tsv" }

    label 'ariba_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads), val(sero)
    val(ariba_db)

    output:
    tuple val(sample_id), path(ariba_report), path(ariba_genes_file), emit: ariba_results_ch

    script:
    read1="${reads[0]}"
    read2="${reads[1]}"
    ariba_report="${sample_id}_ariba_results/report.tsv"
    ariba_genes_file="${sample_id}_ariba_results/assembled_genes.fa.gz"
    """
    # strip genetic variant back down to serotype, probably make into python script so it is tidier at some point..
    serotype=\$(echo ${sero} | awk -F "/" '{ print \$1 }' | awk -F "-" '{ print \$1 }' | awk -F "(" '{ print \$NF }' | sed 's|)||g' | awk -F "_" '{ print \$1 }' | awk '{if (\$0 ~ /[a-zA-Z][0-9]\$/) sub(/[0-9]\$/, ""); print}' | sed 's|19AF|19F|g')
    ariba_sero=\$(grep -m1 \${serotype} ${projectDir}/cps_reference_database/serotypes.txt)
    if [[ \${ariba_sero} == *"("* ]]; then
        ariba_sero=\$(echo \${ariba_sero} | awk -F "(" '{ print \$NF }' | sed 's|)||g')
    fi
    ariba_db_full="${ariba_db}/\${ariba_sero}_ariba_db"
    ariba run \${ariba_db_full} ${read1} ${read2} ${sample_id}_ariba_results --threads 4
    """
}

process FIND_KEY_MUTATIONS {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "key_ariba_mutations.tsv"

    label 'bash_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(ariba_report), path(ariba_genes_file)

    output:
    path(ariba_key_mutations), emit: ariba_mutations_ch
    path(key_mutations_per_sample), emit: key_mutations_sample_ch, optional: true

    script:
    ariba_key_mutations="key_ariba_mutations.tsv"
    key_mutations_per_sample="${sample_id}_key_ariba_mutations.tsv"
    """
    head -1 ${ariba_report} > key_ariba_mutations.tsv
    # catch exit 1 - grep exits with code 1 if there is no match
    grep -i -e "fshift" -e "trunc" -e "ins" -e "del" -e "indels" ${ariba_report} >> key_ariba_mutations.tsv || [[ \$? == 1 ]]
    if [[ \$(wc -l <key_ariba_mutations.tsv) -gt 1 ]]
    then
      tail +2 key_ariba_mutations.tsv | sed "s|^|${sample_id}\t|g" > ${sample_id}_key_ariba_mutations.tsv
    fi
    """
}

process COLLECT_KEY_MUTATIONS {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "key_ariba_mutations_all_samples.tsv"

    label 'bash_container'
    label 'farm_low'

    input:
    path(mutations)

    output:
    path(ariba_key_mutations), emit: ariba_mutations_ch

    script:
    ariba_key_mutations="key_ariba_mutations_all_samples.tsv"
    """
    echo -e "sample\tariba_ref_name\tref_name\tgene\tvar_only\tflag\
    \treads\tcluster\tref_len\tref_base_assembled\tpc_ident\tctg\
    \tctg_len\tctg_cov\tknown_var\tvar_type\tvar_seq_type\
    \tknown_var_change\thas_known_var\tref_ctg_change\t\
    ref_ctg_effect\tref_start\tref_end\tref_nt\tctg_start\
    \tctg_end\tctg_nt\tsmtls_total_depth\tsmtls_nts\
    \tsmtls_nts_depth\tvar_description\tfree_text" > key_ariba_mutations_all_samples.tsv
    cat *_key_ariba_mutations.tsv >> key_ariba_mutations_all_samples.tsv
    """
}

process CHECK_GENE_INTEGRITY {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "gene_integrity_all_samples.csv"

    label 'bash_container'
    label 'farm_low'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(ariba_report), path(ariba_genes_file)

    output:
    path(gene_integrity), emit: gene_integrity_ch
    path(disrupted_genes), emit: disrupted_genes_ch, optional: true

    script:
    gene_integrity="gene_integrity.csv"
    disrupted_genes="${sample_id}_disrupted_genes.csv"
    """
    check_gene_integrity.sh ${ariba_genes_file} ${ariba_report} > gene_integrity.csv
    if grep -q ",disrupted" gene_integrity.csv
    then
      grep ",disrupted" gene_integrity.csv | awk -v sample="${sample_id}" '{print sample"," \$0}' >> ${sample_id}_disrupted_genes.csv
    fi
    """
}

process COLLATE_DISRUPTED_GENES {
    publishDir "${params.output}", mode: 'copy', overwrite: true, pattern: "disrupted_genes_all_samples.csv"

    label 'bash_container'
    label 'farm_low'

    input:
    path(disrupted_genes)

    output:
    path(disrupted_genes_file), emit: disrupted_genes_ch

    script:
    disrupted_genes_file="disrupted_genes_all_samples.csv"
    """
    echo "sample,gene,gene_integrity" > disrupted_genes_all_samples.csv
    cat *_disrupted_genes.csv >> disrupted_genes_all_samples.csv
    """
}

