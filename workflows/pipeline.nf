include { BAKTA } from "$projectDir/modules/bakta"
include { BLASTN } from "$projectDir/modules/blast"
include { CHECK_CPS_SEQUENCE } from "$projectDir/modules/check_cps_sequence"
include { CURATE_CPS_SEQUENCE } from "$projectDir/modules/curate_cps_sequence"

// Main pipeline workflow
workflow PIPELINE {

    main:
        blast_db = file( "${params.blastdb}.n*" )
        blast_db_ch = Channel.fromPath( blast_db )

        prodigal_training_file = Channel.fromPath( params.prodigal_training_file )

        bakta_db = Channel.fromPath( params.bakta_db )

        input_ch = Channel.fromPath( "$params.input/*.{fa,fasta}", checkIfExists: true )
                   .map { tuple( it.name.split('.fa')[0], it ) }

        BLASTN( input_ch, blast_db_ch )

        CURATE_CPS_SEQUENCE( BLASTN.out.blast_results_ch )

        BAKTA( CURATE_CPS_SEQUENCE.out.cps_sequence_ch, prodigal_training_file.first(), bakta_db.first() )

        CHECK_CPS_SEQUENCE( BAKTA.out.bakta_results_ch )
}