include { ASSEMBLY_SHOVILL } from "$projectDir/modules/assembly"
include { BAKTA } from "$projectDir/modules/bakta"
include { BLASTN } from "$projectDir/modules/blast"
include { CHECK_CPS_SEQUENCE } from "$projectDir/modules/check_cps_sequence"
include { CURATE_CPS_SEQUENCE } from "$projectDir/modules/curate_cps_sequence"
include { GAP_FILLER } from "$projectDir/modules/gap_filler"
include { SEROBA } from "$projectDir/modules/serotyping"

// Main pipeline workflow
workflow PIPELINE {

    main:
        blast_db = file( "${params.blastdb}.n*" )
        blast_db_ch = Channel.fromPath( blast_db )

        prodigal_training_file = Channel.fromPath( params.prodigal_training_file )

        bakta_db = Channel.fromPath( params.bakta_db )

        reads_ch = Channel.fromFilePairs( "$params.input/*_{,R}{1,2}{,_001}.{fq.gz,fastq.gz}", checkIfExists: true )

        reference_db_ch = Channel.fromPath( params.reference_database )

        if ( !params.serotype ) {
            SEROBA( reads_ch )
            reads_sero_ch = SEROBA.out.reads_sero_ch
        } else {
            reads_sero_ch = reads_ch.map { it -> it + "$params.serotype" }
        }

        assembly_ch = ASSEMBLY_SHOVILL( reads_sero_ch, params.min_contig_length )

        BLASTN( assembly_ch, blast_db_ch.first() )

        CURATE_CPS_SEQUENCE( BLASTN.out.blast_results_ch )

        GAP_FILLER( CURATE_CPS_SEQUENCE.out.cps_sequence_ch, CURATE_CPS_SEQUENCE.out.log_ch, reference_db_ch.first() )

        BAKTA( GAP_FILLER.out.gap_filled_ch, prodigal_training_file.first(), bakta_db.first() )

        CHECK_CPS_SEQUENCE( BAKTA.out.bakta_results_ch )
}