// Download standard reference database
process DOWNLOAD_REFERENCE_DATABASE {
    publishDir "$projectDir", mode: 'copy', overwrite: true, pattern: "cps_reference_database"
    label 'git_container'
    label 'farm_low'


    output:
    path "cps_reference_database", optional: true
    val true, emit: ready_ch

    script:
    """
    if [ ! -d "${projectDir}/cps_reference_database/annotation" ] || \
       [ ! -d "${projectDir}/cps_reference_database/clinker_descriptions" ] || \
       [ ! -d "${projectDir}/cps_reference_database/fasta" ] || \
       [ ! -d "${projectDir}/cps_reference_database/genbank" ] || \
       [ ! -d "${projectDir}/cps_reference_database/proteins" ] || \
       [ ! -f "${projectDir}/cps_reference_database/all.trn" ] || \
       [ ! compgen -G "${projectDir}/cps_reference_database/cps_blastdb*" > /dev/null ]; then
         git clone https://github.com/GlobalPneumoSeq/cps_reference_database.git
    fi    
    """
}

// Download standard reference database
process DOWNLOAD_BAKTA_DATABASE {
    publishDir "${projectDir}/cps_reference_database", mode: 'copy', overwrite: true, pattern: "bakta_db"
    label 'bakta_container'
    label 'farm_low'

    output:
    path "bakta_db", optional: true

    script:
    """
    if [ ! -d "$projectDir/cps_reference_database/bakta_db/db" ]; then
        bakta_db download -o bakta_db
    fi   
    """
}


// Download standard reference database
process EXTRACT_CDS_FNA {
    label 'bedtools_container'
    label 'farm_low'

    input:
    val download_ready

    output:
    path("*_meta.tsv"), emit: meta_ch
    path("*.fasta"), emit: fasta_ch

    script:
    """
    while read line
    do
      if [[ \$line == *"("* ]]; then
        sero=\$(echo \${line} | awk -F "(" '{ print \$NF }' | sed 's|)||g')
      else
        sero=\${line}
      fi
      extract_gene_sequences.sh \${sero} $projectDir/cps_reference_database
    done < $projectDir/cps_reference_database/serotypes.txt
    """
}

// Download standard reference database
process ARIBA_PREPARE_REF {
    publishDir "${projectDir}/cps_reference_database/ariba_databases", mode: 'copy', overwrite: true, pattern: "*_ariba_db"
    label 'ariba_container'
    label 'farm_low'

    input:
    tuple (val(sero), path(fasta), path(metadata))

    output:
    path(ariba_db), emit: db_ch

    script:
    ariba_db="${sero}_ariba_db"
    """
    ariba prepareref -f ${fasta} -m ${metadata} ${sero}_ariba_db
    """
}