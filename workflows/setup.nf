include { ARIBA_PREPARE_REF; DOWNLOAD_REFERENCE_DATABASE; DOWNLOAD_BAKTA_DATABASE; EXTRACT_CDS_FNA } from "$projectDir/modules/setup"
include { GET_DOCKER_COMPOSE; PULL_IMAGES} from "$projectDir/modules/docker"

// Setup workflow
workflow SETUP {

    main:
        DOWNLOAD_REFERENCE_DATABASE()

        DOWNLOAD_BAKTA_DATABASE()

        EXTRACT_CDS_FNA(DOWNLOAD_REFERENCE_DATABASE.out.ready_ch)

        ariba_meta_ch = EXTRACT_CDS_FNA.out.meta_ch.flatten()
                                           .map{ it -> tuple(it.baseName.split('_')[0], it) }

        ariba_ref_ch = EXTRACT_CDS_FNA.out.fasta_ch.flatten()
                                          .map{ it -> tuple(it.baseName.split('_')[0], it) }

        ariba_ch = ariba_ref_ch.join(ariba_meta_ch)

        ARIBA_PREPARE_REF(ariba_ch)

        // Pull all Docker images used in the workflow if using Docker
        if (workflow.containerEngine === 'docker') {
            GET_DOCKER_COMPOSE(
                Channel.fromList(workflow.container.collect { it.value })
                    .unique()
                    .collectFile(name: 'containersList.txt', newLine: true)
            )

            PULL_IMAGES(GET_DOCKER_COMPOSE.out.compose)
        }
}