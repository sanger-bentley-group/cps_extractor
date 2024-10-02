// Download standard reference database
process DOWNLOAD_REFERENCE_DATABASE {
    publishDir "$projectDir", mode: 'copy', overwrite: true, pattern: "cps_reference_database"
    label 'git_container'
    label 'farm_low'


    output:
    path "cps_reference_database", optional: true

    script:
    """
    if [ ! -d "${projectDir}/cps_reference_database/annotation" ] || \
       [ ! -d "${projectDir}/cps_reference_database/clinker_descriptions" ] || \
       [ ! -d "${projectDir}/cps_reference_database/fasta" ] || \
       [ ! -d "${projectDir}/cps_reference_database/genbank" ] || \
       [ ! -d "${projectDir}/cps_reference_database/proteins" ] || \
       [ ! -f "${projectDir}/cps_reference_database/all.trn" ] || \
       [ ! compgen -G "${projectDir}/cps_reference_database/cps_blastdb*" > /dev/null ]; then
         git clone https://github.com/Oliver-Lorenz-dev/cps_reference_database.git
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