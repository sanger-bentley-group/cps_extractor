include { DOWNLOAD_REFERENCE_DATABASE; DOWNLOAD_BAKTA_DATABASE } from "$projectDir/modules/setup"

// Main pipeline workflow
workflow SETUP {

    main:
        DOWNLOAD_REFERENCE_DATABASE()

        DOWNLOAD_BAKTA_DATABASE()
}