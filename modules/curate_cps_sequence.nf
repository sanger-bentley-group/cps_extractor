// Curate CPS sequence
process CURATE_CPS_SEQUENCE {
    publishDir "${params.output}/${sample_id}", mode: 'copy', overwrite: true, pattern: "*.log"

    label 'cps_extractor_python_container'
    label 'farm_low_fallible'
    label 'farm_scratchless'

    cache 'lenient'

    errorStrategy 'ignore'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(assembly), path(blast_results), path(dexb_results), path(alia_results), val(serotype), path(reads)
    val (results_dir)

    output:
    tuple val(sample_id), path(cps_sequence), path(reads), val(serotype), env(filled), emit: cps_sequence_ch
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
        if ! curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa -s \${sero_final} -g ${assembly} -a ${alia_results} -d ${dexb_results}
        then
          echo "No cps sequence found for \${sero_final}, running without specifying the serotype, please check the log file for more information."
          curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa -g ${assembly} -a ${alia_results} -d ${dexb_results} || { cp cps_extractor.log ${results_dir}/${sample_id}; exit 1; }
        fi
    else
        curate_cps_sequence.py -b ${blast_results} -o ${sample_id}_cps.fa -g ${assembly} -a ${alia_results} -d ${dexb_results} || { cp cps_extractor.log ${results_dir}/${sample_id}; exit 1; }
    fi

    # if a gap is in the same contig, fill it using position indexing rather than consensus method
    if [[ -f "gap_filled" ]]
    then
        filled=true
    else
        filled=false
    fi
    echo array
    """
}