// Check CPS sequence for disruptive mutations
process CHECK_CPS_SEQUENCE {
    
    label 'bash_container'
    label 'farm_low_fallible'

    errorStrategy 'ignore'

    cache 'lenient'

    tag "$sample_id"

    input:
    tuple val(sample_id), path(cps_sequence), val(reference)
    path(results_dir)

    output:
    tuple val(sample_id), path(cps_sequence), val(reference), emit: results_ch

    script:
    """
    # check length of cps sequence is not unusually low before carrying on with the main pipeline
    base_count=\$(grep -v ">" ${cps_sequence} | grep -o -E "A|C|G|T" | wc -l)
    echo \${base_count}
    if [ "\${base_count}" -lt "${params.minimum_cps_length}" ]
    then
      echo -e "The CPS sequence length is unusually low (\${base_count} bases), please check the blast results file you may have a non capsulated sample or a pneumo 'like' sample" >> ${results_dir}/${sample_id}/cps_extractor.log
      exit 1
    fi
    """
}