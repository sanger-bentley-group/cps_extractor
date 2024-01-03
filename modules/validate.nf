import groovy.json.JsonSlurper

// Map of valid parameters and their value types
validParams = [
    help: 'boolean',
    version: 'boolean',
    bakta_db: 'path_exist',
    bakta_threads: 'int',
    blastdb: 'path_blast_db',
    input: 'path_exist',
    output: 'path',
    lite: 'boolean',
    singularity_cachedir: 'path',
    maxretries: 'int',
    prodigal_training_file: 'path_exist',
    min_contig_length: 'int',
    serotype: 'string'
]

// Validate whether all provided parameters are valid
void validate(Map params) {
    // Ensure only one or none of the alternative workflows is selected
    if ([params.help, params.version].count { it } > 1) {
        log.error('''
            |More than one alternative workflow is selected, please only select one of them
            |(Only one of --version, --help should be used at one time)
            '''.stripMargin())
        System.exit(1)
}

    // Skip validation when help option is used
    if (params.help) {
        return
    }

    // Skip validation when version option is used
    if (params.version) {
        return
    }

    // Add params.singularity_cachedir when workflow.containerEngine == 'singularity'
    if (workflow.containerEngine == 'singularity') {
        validParams.put("singularity_cachedir", "path")
    }

    // Add params.maxretries when workflow.profile contains 'lsf' 
    if (workflow.profile.split(',').contains('lsf')) {
        validParams.put("maxretries", "int")
    }

    // For version, skip all file paths related checks
    skippedParams = []
    if (params.version) {
        validParams.each {
            key, value ->
            if (['path', 'path_exist', 'path_fasta'].contains(value)) {
                skippedParams.add(key)
            }
        }
    }
    skippedParams.each { key -> validParams[key] = 'skip' }

    // To save invalid parameters in this list
    invalidParams = []
    // To save invalid parameter values as "parameter : [value, issue]" in this map
    invalidValues = [:]

    params.each {
        key, value ->

        // If parameter is invalid, add it to invalidParams list and skip the following checks
        if (!validParams.keySet().contains(key)) {
            invalidParams.add(key)
            return
        }

        // Based on the value type of the parameter, perform the appropriate check
        switch (validParams[key]) {
            case 'skip':
                break

            case 'boolean':
                if (value !instanceof Boolean) {
                    invalidValues[key] = [value, 'boolean value']
                }
                break

            case 'int':
                if (value !instanceof Integer) {
                    invalidValues[key] = [value, 'integer value']
                }
                break

            case 'string':
                if (value !instanceof String) {
                    invalidValues[key] = [value, 'string value']
                }
                break

            case 'int_float':
                if (value !instanceof Integer && value !instanceof BigDecimal && value !instanceof Double) {
                    invalidValues[key] = [value, 'integer or float value']
                }
                break

            case 'publish_mode':
                if (!['link', 'symlink', 'copy'].contains(value)) {
                    invalidValues[key] = [value, 'Nextflow publish mode']
                }
                break

            case 'path_exist':
                File dir = new File(value)
                if (!dir.exists()) {
                    invalidValues[key] = [value, 'existing directory']
                }
                break

            case 'path':
                File dir = new File(value)
                if (!(dir.exists() || dir.mkdirs())) {
                    invalidValues[key] = [value, 'directory path (invalid path or insufficient permissions)']
                }
                break

            case 'path_blast_db':
                // check blast db exists with glob matching

                // split path to get base directory
                if (value.contains("/")) {
                def parts = value.split('/')
                def baseDir = parts[0..-2]

                // get db name
                def db_name = parts[-1]

                // check if blast db exists
                def basePath = baseDir.join('/')
                File dir = new File(basePath)
                def files = dir.list()
                def exists = files.any { it.contains(db_name) }
                if (!exists) {
                    invalidValues[key] = [value, 'Blast db does not exist']
                }
                }
                else {
                    File dir = new File(".")
                    def files = dir.list()
                    def exists = files.any { it.contains(value) }
                    if (!exists) {
                        invalidValues[key] = [value, 'Blast db does not exist']
                    }
                }
                break

            case 'path_fasta':
                File fasta = new File(value)
                if (!fasta.exists()) {
                    invalidValues[key] = [value, 'path to a fasta file (file does not exist)']
                } else if (!(value ==~ /.+\.(fa|fasta)$/)) {
                    invalidValues[key] = [value, 'path to a fasta file (file does not have an filename extension of .fasta or .fa)']
                }
                break
            
            case 'path_tsv':
                File tsv = new File(value)
                if (!tsv.exists()) {
                    invalidValues[key] = [value, 'path to a TSV file (file does not exist)']
                } else if (!(value ==~ /.+\.tsv$/)) {
                    invalidValues[key] = [value, 'path to a TSV file (file does not have an filename extension of .tsv)']
                }
                break

            case 'url_targz':
                if (!(value ==~ /^(https?:\/\/)?(?:www\.)?[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}\b(?:[-a-zA-Z0-9()@:%_\+.~#?&\/=]*)\.(tar\.gz|tgz)$/)) {
                    invalidValues[key] = [value, 'URL that points a .tar.gz file (valid URL ending with .tar.gz or .tgz)']
                }
                break

            case 'url_csv':
                if (!(value ==~ /^(https?:\/\/)?(?:www\.)?[-a-zA-Z0-9@:%._\+~#=]{1,256}\.[a-zA-Z0-9()]{1,6}\b(?:[-a-zA-Z0-9()@:%_\+.~#?&\/=]*)\.csv$/)) {
                    invalidValues[key] = [value, 'URL that points a .csv file (valid URL ending with .csv)']
                }
                break

            // Should only reach this statement if a new value type is added to validParams without adding its case above
            default:
                log.error("""
                    |Unknown value type \"${validParams[key]}\"
                    |Please submit an issue at \"https://github.com/HarryHung/gps-unified-pipeline/issues\"}
                    """.stripMargin())
                System.exit(1)
        }
    }

    // If invalidParams list or invalidValues map is not empty, log error messages and terminate the pipeline
    if (invalidParams || invalidValues) {
        log.error('The pipeline will now be terminated due to the following critical error(s):')

        if (invalidParams) {
            log.error("The following invalid option(s) were provided: --${invalidParams.join(', --')}.")
        }

        if (invalidValues) {
            invalidValues.each {
                key, values ->
                log.error("The provided value \"${values[0]}\" for option --${key} is not a valid ${values[1]}.")
            }
        }

        System.exit(1)
    }
}
