include { IMAGES; TOOLS; COMBINE_INFO; PARSE; PRINT; SAVE; PYTHON_VERSION; BWA_VERSION; SAMTOOLS_VERSION; BLAST_VERSION; BEDTOOLS_VERSION; BAKTA_VERSION } from "$projectDir/modules/info"

// Alternative workflow that prints versions of pipeline and tools
workflow PRINT_VERSION {
    take:
        pipeline_version

    main:
        GET_VERSION(
            pipeline_version
        ) \
        | PARSE \
        | PRINT
}

// Sub-workflow of PIPELINE workflow the save versions of pipeline and tools, and QC parameters to info.txt at output dir
workflow SAVE_INFO {
    take:
        pipeline_version

    main:
        GET_VERSION(
            pipeline_version
        ) \
       | PARSE \
       | SAVE
}

// Sub-workflow for generating a json that contains versions of pipeline and tools
workflow GET_VERSION {
    take:
        pipeline_version

    main:
        IMAGES(Channel.fromPath("${workflow.configFiles[0]}"))

        nextflow_version = "$nextflow.version"

        BAKTA_VERSION()
        BEDTOOLS_VERSION()
        BLAST_VERSION()
        BWA_VERSION()
        PYTHON_VERSION()
        SAMTOOLS_VERSION()

        TOOLS(
            BAKTA_VERSION.out,
            BEDTOOLS_VERSION.out,
            BLAST_VERSION.out,
            BWA_VERSION.out,
            PYTHON_VERSION.out,
            SAMTOOLS_VERSION.out
        )

        COMBINE_INFO(
            pipeline_version,
            nextflow_version,
            IMAGES.out.json,
            TOOLS.out.json
        )

    emit:
        COMBINE_INFO.out.json
}
