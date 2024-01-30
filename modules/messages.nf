// Start message
void startMessage(String pipelineVersion) {
    log.info( 
        $/
        |╔═════════════════════════════════════════════════════════════════════════════╗
        |║   ____ ____  ____    _______  _______ ____      _    ____ _____ ___  ____   ║
        |║  / ___|  _ \/ ___|  | ____\ \/ /_   _|  _ \    / \  / ___|_   _/ _ \|  _ \  ║
        |║ | |   | |_) \___ \  |  _|  \  /  | | | |_) |  / _ \| |     | || | | | |_) | ║
        |║ | |___|  __/ ___) | | |___ /  \  | | |  _ <  / ___ \ |___  | || |_| |  _ <  ║
        |║  \____|_|   |____/  |_____/_/\_\ |_| |_| \_\/_/   \_\____| |_| \___/|_| \_\ ║
        |${String.format('║  v %-57s                ║', pipelineVersion)}
        |╚═════════════════════════════════════════════════════════════════════════════╝
       /$.stripMargin()
    )
}

// Help message
void helpMessage() {
    log.info(
        '''
        |This is a Nextflow Pipeline for extracting CPS sequences from S pneumoniae reads
        |
        |Usage:
        |./run_cps_extractor [option] [value]
        |
        |--input [PATH]                  Path to the input directory that contains reads to be processed. Default: ./input
        |--output [PATH]                 Path to the output directory that save the results. Default: output
        |--serotype [STR]                Serotype (if known). Default: None
        |--setup                         Alternative workflow for setting up the required databases.
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
        case 'setup':
            message = '''
            |The alternative workflow for setting up the database has been selected.
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
        case 'setup':
            successMessage = '''
                |The databases have been correctly setup.
                '''.stripMargin()
            failMessage = '''
                |Failed to setup the database.
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
