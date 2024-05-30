include { DOWNLOAD_REFERENCE_DATABASE; DOWNLOAD_BAKTA_DATABASE } from "$projectDir/modules/setup"
include { GET_DOCKER_COMPOSE; PULL_IMAGES} from "$projectDir/modules/docker"
// Setup workflow
workflow SETUP {

    main:
        DOWNLOAD_REFERENCE_DATABASE()

        DOWNLOAD_BAKTA_DATABASE()

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