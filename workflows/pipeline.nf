include { ASSEMBLY_SHOVILL } from "$projectDir/modules/assembly"
include { BAKTA } from "$projectDir/modules/bakta"
include { BLASTN } from "$projectDir/modules/blast"
include { CHECK_CPS_SEQUENCE } from "$projectDir/modules/check_cps_sequence"
include { CURATE_CPS_SEQUENCE } from "$projectDir/modules/curate_cps_sequence"
include { GAP_FILLER } from "$projectDir/modules/gap_filler"
include { CHECK_GENE_ORDER; CLINKER; CLINKER_ALL; PANAROO_ALL; PANAROO_REF_COMPARISON; RENAME_PANAROO_ALIGNMENTS; SNP_DISTS } from "$projectDir/modules/gene_comparison"
include { CONCAT_PROTEIN_SEQUENCES; CREATE_PROTEIN_FILES; EXTRACT_PROTEIN_SEQUENCES } from "$projectDir/modules/proteins"
include { SEROBA } from "$projectDir/modules/serotyping"


// Main pipeline workflow
workflow PIPELINE {

    main:
        blast_db = file( "${params.blastdb}.n*" )
        blast_db_ch = Channel.fromPath( blast_db )

        results_dir_ch = Channel.fromPath( params.output )

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

        CURATE_CPS_SEQUENCE( BLASTN.out.blast_results_ch, results_dir_ch.first() )

        GAP_FILLER( CURATE_CPS_SEQUENCE.out.cps_sequence_ch, CURATE_CPS_SEQUENCE.out.log_ch, reference_db_ch.first() )

        BAKTA( GAP_FILLER.out.gap_filled_ch, prodigal_training_file.first(), bakta_db.first(), reference_db_ch.first() )

        CHECK_CPS_SEQUENCE( BAKTA.out.bakta_results_ch, results_dir_ch.first() )

        EXTRACT_PROTEIN_SEQUENCES( CHECK_CPS_SEQUENCE.out.results_ch )

        CREATE_PROTEIN_FILES( EXTRACT_PROTEIN_SEQUENCES.out.proteins_ch )

        CHECK_GENE_ORDER( EXTRACT_PROTEIN_SEQUENCES.out.results_ch, reference_db_ch.first() )

        CLINKER( CHECK_CPS_SEQUENCE.out.results_ch, reference_db_ch.first() )

        PANAROO_REF_COMPARISON( EXTRACT_PROTEIN_SEQUENCES.out.results_ch, reference_db_ch.first() )

        RENAME_PANAROO_ALIGNMENTS( PANAROO_REF_COMPARISON.out.panaroo_results_ch )

        SNP_DISTS( RENAME_PANAROO_ALIGNMENTS.out.gene_alignment_results.transpose() )

        if ( params.serotype ) {
            // Collect all annotations and run panaroo on them
            annotation_ch = CHECK_CPS_SEQUENCE.out.results_ch.map { it -> it[1] }.collect()
            PANAROO_ALL( annotation_ch, reference_db_ch, params.serotype )

            // Collect all genbank files and generate plot using clinker
            genbank_ch = CHECK_CPS_SEQUENCE.out.results_ch.map { it -> it[2] }.collect()
            CLINKER_ALL( genbank_ch, reference_db_ch, params.serotype )

            // Collect all protein sequences and concatenate them
            proteins_ch = CREATE_PROTEIN_FILES.out.proteins_ch
                                                              .collect()
                                                              .flatten()
                                                              .map { file -> tuple(file.baseName.split('-')[-1], file) }
                                                              .groupTuple()
                                                              
            CONCAT_PROTEIN_SEQUENCES( proteins_ch )
        }
}