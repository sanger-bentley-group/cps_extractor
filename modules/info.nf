// Import for PARSE process
import groovy.json.JsonSlurper

// Extract containers information from nextflow.config and save into a JSON file
process IMAGES {
    label 'bash_container'
    label 'farm_low'

    input:
    path nextflowConfig

    output:
    path(json), emit: json

    script:
    json='images.json'
    """
    NEXTFLOW_CONFIG="$nextflowConfig"
    JSON_FILE="$json"

    source save_images_info.sh
    """
}

// Save received tools versions into a JSON file
process TOOLS {
    label 'bash_container'
    label 'farm_low'

    input:
    val bakta_version
    val bcftools_version
    val bedtools_version
    val blast_version
    val bwa_version
    val panaroo_version
    val python_version
    val samtools_version
    val seroba_version
    val shovill_version
    val snpdists_version

    output:
    path(json), emit: json

    script:
    json='tools.json'
    """
    BAKTA_VERSION="$bakta_version"
    BCFTOOLS_VERSION="$bcftools_version"
    BEDTOOLS_VERSION="$bedtools_version"
    BWA_VERSION="$bwa_version"
    SAMTOOLS_VERSION="$samtools_version"
    BLAST_VERSION="$blast_version"
    PANAROO_VERSION="$panaroo_version"
    PYTHON_VERSION="$python_version"
    SEROBA_VERSION="$seroba_version"
    SHOVILL_VERSION="$shovill_version"
    SNPDISTS_VERSION="$snpdists_version"
    JSON_FILE="$json"
                
    source save_tools_info.sh
    """
}

// Combine pipeline version, Nextflow version, databases information, container images, tools version JSON files into the a single JSON file
process COMBINE_INFO {
    label 'bash_container'
    label 'farm_low'

    input:
    val pipeline_version
    val nextflow_version
    path images
    path tools

    output:
    path(json), emit: json

    script:
    json='result.json'
    """
    PIPELINE_VERSION="$pipeline_version"
    NEXTFLOW_VERSION="$nextflow_version"
    IMAGES="$images"
    TOOLS="$tools"
    JSON_FILE="$json"

    source save_combined_info.sh
    """
}

// Parse information from JSON into human-readable tables
process PARSE {
    label 'farm_local'

    input:
    val json_file

    output:
    val coreText
    val toolText
    val imageText

    exec:
    def jsonSlurper = new JsonSlurper()

    def json = jsonSlurper.parse(new File("${json_file}"))

    def textRow = { leftSpace, rightSpace, leftContent, rightContent ->
        String.format("║ %-${leftSpace}s │ %-${rightSpace}s ║", leftContent, rightContent)
    }

    def coreTextRow = { leftContent, rightContent ->
        textRow(25, 67, leftContent, rightContent)
    }

    coreText = """\
        |┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈ Core Software Versions ┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
        |╔═══════════════════════════╤═════════════════════════════════════════════════════════════════════╗
        |${coreTextRow('Software', 'Version')}
        |╠═══════════════════════════╪═════════════════════════════════════════════════════════════════════╣
        |${coreTextRow('CPS Extractor Pipeline', json.pipeline.version)}
        |${coreTextRow('Nextflow', json.nextflow.version)}
        |╚═══════════════════════════╧═════════════════════════════════════════════════════════════════════╝
        |""".stripMargin()

    def getVersion = { tool ->
        if (json[tool] && json[tool]['version']) {
            return json[tool]['version']
        }

        return 'no version information'
    }

    def toolTextRow = { leftContent, rightContent ->
        textRow(30, 62, leftContent, getVersion(rightContent))
    }

    toolText = """\
        |┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈ Tool Versions ┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
        |╔════════════════════════════════╤════════════════════════════════════════════════════════════════╗
        |${textRow(30, 62, 'Tool', 'Version')}
        |╠════════════════════════════════╪════════════════════════════════════════════════════════════════╣
        |${toolTextRow('Bakta', 'bakta')}
        |${toolTextRow('Bcftools', 'bcftools')}
        |${toolTextRow('Bedtools', 'bedtools')}
        |${toolTextRow('Blast', 'blast')}
        |${toolTextRow('BWA', 'bwa')}
        |${toolTextRow('Panaroo', 'panaroo')}
        |${toolTextRow('python', 'python')}
        |${toolTextRow('SAMtools', 'samtools')}
        |${toolTextRow('SeroBA', 'seroba')}
        |${toolTextRow('Shovill', 'shovill')}
        |${toolTextRow('SNP_Dists', 'snpdists')}
        |╚════════════════════════════════╧════════════════════════════════════════════════════════════════╝
        |""".stripMargin()

    def getImage = { tool ->
        if (json[tool] && json[tool]['container']) {
            return json[tool]['container']
        }

        return 'no image information'
    }

    def imageTextRow = { leftContent, rightContent ->
        textRow(30, 62, leftContent, getImage(rightContent))
    }

    imageText = """\
        |┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈ Container Images ┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
        |╔════════════════════════════════╤════════════════════════════════════════════════════════════════╗
        |${textRow(30, 62, 'Environment For', 'Image')}
        |╠════════════════════════════════╪════════════════════════════════════════════════════════════════╣
        |${imageTextRow('Bakta', 'bakta')}
        |${imageTextRow('Bash', 'bash')}
        |${imageTextRow('Blast', 'blast')}
        |${imageTextRow('Check Gene Content', 'check_gene_content')}
        |${imageTextRow('CPS extractor python', 'cps_extractor_python')}
        |${imageTextRow('Gap Filler', 'gap_filler')}
        |${imageTextRow('Panaroo', 'panaroo')}
        |${imageTextRow('SeroBA', 'seroba')}
        |${imageTextRow('Shovill', 'shovill')}
        |${imageTextRow('SNP_Dists', 'snp_dists')}
        |╚════════════════════════════════╧════════════════════════════════════════════════════════════════╝
        |""".stripMargin()
}

// Print parsed version information
process PRINT {
    label 'farm_local'

    input:
    val coreText
    val toolText
    val imageText

    exec:
    log.info(
        """
        |╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍
        |╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍ Version Information ╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍
        |╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍╍
        |
        |${coreText}
        |${toolText}
        |${imageText}
        |""".stripMargin()
    )
}

// Save core software, I/O, assembler, QC parameters, databases, tools, container engine and images information to info.txt at output dir
process SAVE {
    label 'farm_local'
    
    publishDir "${params.output}", mode: "copy"

    input:
    val coreText
    val toolText
    val imageText

    output:
    path "info.txt", emit: info

    exec:
    File inputDir = new File(params.input)
    File outputDir = new File(params.output)

    def textRow = { leftSpace, rightSpace, leftContent, rightContent ->
        String.format("║ %-${leftSpace}s │ %-${rightSpace}s ║", leftContent, rightContent)
    }

    def ioTextRow = { leftContent, rightContent ->
        textRow(8, 84, leftContent, rightContent)
    }

    String ioText = """\
    |┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈ Input and Output ┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
    |╔══════════╤══════════════════════════════════════════════════════════════════════════════════════╗
    |${ioTextRow('Type', 'Path')}
    |╠══════════╪══════════════════════════════════════════════════════════════════════════════════════╣
    |${ioTextRow('Input', inputDir.canonicalPath)}
    |${ioTextRow('Output', outputDir.canonicalPath)}
    |╚══════════╧══════════════════════════════════════════════════════════════════════════════════════╝
    |""".stripMargin()

    def containerEngineTextRow = { leftContent, rightContent ->
        textRow(25, 67, leftContent, rightContent)
    }

    String containerEngineText = """\
    |┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈ Container Engine Options ┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
    |╔═══════════════════════════╤═════════════════════════════════════════════════════════════════════╗
    |${containerEngineTextRow('Option', 'Value')}
    |╠═══════════════════════════╪═════════════════════════════════════════════════════════════════════╣
    |${containerEngineTextRow('Container Engine', workflow.containerEngine.capitalize())}
    |╚═══════════════════════════╧═════════════════════════════════════════════════════════════════════╝
    |""".stripMargin()

    File output = new File("${task.workDir}/info.txt")
    output.write(
        """\
        |${coreText}
        |${ioText}
        |${toolText}
        |${containerEngineText}
        |${imageText}
        |""".stripMargin()
    )
}

// Below processes get tool versions within container images by running their containers

process PYTHON_VERSION {
    label 'cps_extractor_python_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    $/
    VERSION=$(python3 --version | sed -r "s/.*\s(.+)/\1/")
    /$
}

process PANAROO_VERSION {
    label 'panaroo_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(panaroo --version | awk '{ print $NF }')
    '''
}

process SNP_DISTS_VERSION {
    label 'snp_dists_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(snp-dists -v | awk '{ print $NF }')
    '''
}

process BWA_VERSION {
    label 'gap_filler_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    $/
    VERSION=$(bwa 2>&1 | grep Version | sed -r "s/.*:\s(.+)/\1/")
    /$
}

process SAMTOOLS_VERSION {
    label 'gap_filler_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    $/
    VERSION=$(samtools 2>&1 | grep Version | sed -r "s/.*:\s(.+)\s\(.+/\1/")
    /$
}

process BCFTOOLS_VERSION {
    label 'gap_filler_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    $/
    VERSION=$(bcftools 2>&1 | grep Version | sed -r "s/.*:\s(.+)\s\(.+/\1/")
    /$
}

process BLAST_VERSION {
    label 'blast_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(blastn -version | head -1 | awk '{ print $NF }' | sed 's|+||g')
    '''
}

process BEDTOOLS_VERSION {
    label 'cps_extractor_python_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(bedtools --version | awk '{ print $NF }' | sed 's|v||g')
    '''
}

process BAKTA_VERSION {
    label 'bakta_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(bakta --version | awk '{ print $NF }')
    '''
}

process SHOVILL_VERSION {
    label 'shovill_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(shovill --version | awk '{ print $NF }')
    '''
}

process SEROBA_VERSION {
    label 'seroba_container'
    label 'farm_low'

    output:
    env VERSION

    shell:
    '''
    VERSION=$(seroba version)
    '''
}