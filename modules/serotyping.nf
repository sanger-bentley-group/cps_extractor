// Serotype sample using seroBA
process SEROBA {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*_serotype_report.csv"
    label 'seroba_container'
    label 'farm_mid'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path(reads), env(SEROTYPE), emit: reads_sero_ch
    path(serotype_report)

    script:
    serotype_report='seroba_serotype_report.csv'
    read1="${reads[0]}"
    read2="${reads[1]}"
    """
    SEROBA_DB="/seroba/database"
    READ1="$read1"
    READ2="$read2"
    SAMPLE_ID="$sample_id"
    SEROTYPE_REPORT="$serotype_report"

    source get_serotype.sh
    SEROTYPE=\$(tail -1 $serotype_report)
    """
}