// Curate CPS sequence
process CURATE_CPS_SEQUENCE {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*{_cps.fa,log}"

    label 'cps_extractor_python'
    label 'farm_low'
    // ignore samples which don't have a CPS sequence so the pipeline doesn't terminate
    errorStrategy 'ignore'
    tag "$sample_id"

    input:
    tuple val(sample_id), path(assembly), path(blast_results), val(serotype)

    output:
    tuple val(sample_id), path(cps_sequence), emit: cps_sequence_ch
    path("*.log"), emit: log_ch

    script:
    cps_sequence="${sample_id}_cps.fa"
    """
    set +eu
    sero=\$(echo "${serotype}" | cut -c1-1)
    # take the first serotype if seroBA gives multiple e.g. 6A/6B, use deepseq variant to extract multiple serotypes
    sero_final=\$(echo "${serotype}" | sed 's|[^a-zA-Z0-9].*||')
    # check sero is an integer, seroBA can output some possible/maybe types
    if [ "\${sero}" -eq "\${sero}" ] 2>/dev/null
    then
        # if no cps sequence can be found for the serotype specified, run without specifying the serotype
        if ! curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa -s \${sero_final}
        then
          echo "No cps sequence found for \${sero_final}, running without specifying the serotype, please check the log file for more information."
          curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa
        fi
    else
        curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa
    fi
    """
}