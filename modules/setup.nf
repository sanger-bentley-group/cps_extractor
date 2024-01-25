// Download standard reference database
process DOWNLOAD_REFERENCE_DATABASE {
    publishDir "$projectDir", mode: 'copy', overwrite: true, pattern: "cps_reference_database"
    label 'git_container'
    label 'farm_scratchless'
    label 'farm_low'


    output:
    path "cps_reference_database"

    script:
    """
    git clone https://github.com/Oliver-Lorenz-dev/cps_reference_database.git
    rm -rf cps_reference_database/.git
    """
}

// Download standard reference database
process DOWNLOAD_BAKTA_DATABASE {
    publishDir "${projectDir}/cps_reference_database", mode: 'copy', overwrite: true, pattern: "bakta_db"
    label 'bakta_container'
    label 'farm_scratchless'
    label 'farm_low'

    output:
    path "bakta_db"

    script:
    """
    bakta_db download
    mv db bakta_db
    """
}