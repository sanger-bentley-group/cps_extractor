// Start message
void startMessage(String pipelineVersion) {
    log.info( 
        $/
        |
        |╔══════════════════════════════════════════════════════════════════════════════════════════╗
        |║                                                                                          ║░
        |║  ____  _            _ _                                                                  ║░
        |║ |  _ \(_)_ __   ___| (_)_ __   ___                                                       ║░
        |║ | |_) | | '_ \ / _ | | | '_ \ / _ \                                                      ║░
        |║ |  __/| | |_) |  __| | | | | |  __/                                                      ║░
        |║ |_|   |_| .__/ \___|_|_|_| |_|\___|                                                      ║░
        |${String.format('║  v %-57s |_|                                                            ║░', pipelineVersion)}
        |╚══════════════════════════════════════════════════════════════════════════════════════════╝░
        |  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
       /$.stripMargin()
    )
}

// Help message
void helpMessage() {
    log.info(
        '''
        |This is a Nextflow Pipeline for extracting CPS sequences from S pneumoniae genome assemblies
        |
        |Usage:
        |nextflow run . [option] [value]
        |
        |--bakta_db [PATH]               Path to bakta database. Default: /data/pam/software/bakta/v5
        |--blastdb [PATH]                Path to blast database. Default: cps_blastdb
        |--input [PATH]                  Path to the input directory that contains reads to be processed. Default: ./input
        |--output [PATH]                 Path to the output directory that save the results. Default: output
        |--serotype [STR]                Serotype (if known). Default: None
        |--prodigal_training_file [PATH] Path to prodigal training file used in annotation. Default: all.trn
        |--version                       Alternative workflow for getting versions of pipeline, container images, tools and databases
        |--help                          Print this help message
        |
        |For more information, please refer to README.md
        '''.stripMargin()
    )
}

// Workflow selection message
void workflowSelectMessage(String selectedWorkflow) {
    String message
    File inputDir = new File(params.input)
    File outputDir = new File(params.output)

    switch (selectedWorkflow) {
        case 'pipeline':
            message = """
            |The main pipeline workflow has been selected.
            |
            |Input Directory: ${inputDir.canonicalPath}
            |Output Directory: ${outputDir.canonicalPath}
            """.stripMargin()
            break
        case 'version':
            message = '''
            |The alternative workflow for getting versions of pipeline, tools and databases has been selected.
            '''.stripMargin()
            break
    }

    Date date = new Date()
    String dateStr = date.format('yyyy-MM-dd')
    String timeStr = date.format('HH:mm:ss z')

    log.info(
        """
        |${message}
        |The workflow started at ${dateStr} ${timeStr}.
        |
        |Current Progress:
        """.stripMargin()
    )
}

// End message
void endMessage(String selectedWorkflow) {
    String successMessage
    String failMessage
    File outputDir = new File(params.output)

    switch (selectedWorkflow) {
        case 'pipeline':
            successMessage = """
                |The pipeline has been completed successfully.
                |The output is located at ${outputDir.canonicalPath}.
                """.stripMargin()
            failMessage = '''
                |The pipeline has failed.
                |If you think it is caused by a bug, contact Oliver Lorenz (ol6@sanger.ac.uk)\".
                '''.stripMargin()
            break
        case 'version':
            successMessage = '''
                |All the version information is printed above.
                '''.stripMargin()
            failMessage = '''
                |Failed to get version information on pipeline, tools or databases.
                |If you think it is caused by a bug, contact Oliver Lorenz (ol6@sanger.ac.uk)\"
                '''.stripMargin()
            break
    }

    Date date = new Date()
    String dateStr = date.format('yyyy-MM-dd')
    String timeStr = date.format('HH:mm:ss z')

    log.info(
        """
        |The workflow ended at ${dateStr} ${timeStr}. The duration was ${workflow.duration}
        |${workflow.success ? successMessage : failMessage}
        """.stripMargin()
    )
}
