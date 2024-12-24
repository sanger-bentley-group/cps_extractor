// Run ARIBA to check for mutations in key cps genes
process ARIBA {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, saveAs: { filename -> "ariba_report.tsv" }

    label 'ariba_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads), val(sero)
    val(ariba_db)

    output:
    path(ariba_report), emit: ariba_results_ch

    script:
    read1="${reads[0]}"
    read2="${reads[1]}"
    ariba_report="${sample_id}_ariba_results/report.tsv"
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