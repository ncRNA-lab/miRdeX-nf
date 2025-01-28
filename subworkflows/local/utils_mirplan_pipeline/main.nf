//
// Subworkflow with functionality specific to the mirplan-nf pipeline. This
// subworkflow is based on the utils_nfcore_rnaseq_pipeline subworkflow from
// the nf-core/rnaseq pipeline (https://github.com/nf-core/rnaseq), with minor
// modifications to adapt it for the miRPlan-nf pipeline.
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    
    take:
    version             // boolean: Display version and exit
    validate_params     // boolean: Boolean whether to validate parameters against the schema at runtime
    nextflow_cli_args   //   array: List of positional nextflow CLI args
    outdir              //  string: The output directory where the results will be saved

    main:

    // Empty channel for versions
    ch_versions = Channel.empty()

    // Print version and exit if required and dump pipeline parameters to JSON file
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda']).size() >= 1
    )

    // Validate parameters and generate parameter summary to stdout
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    // Custom validation for pipeline parameters
    validateInputParameters()

    emit:
    version = ch_versions

}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Check and validate pipeline parameters
//

def validateInputParameters() {

    // Create an empty list for the parameters that should be ignored in the
    // conditionals.
    def ignoreParams = []

    // Check the correct usage of the annotation parameters
    if (params.annotation_unitas) {
        // If annotation_unitas is true, the rest of the annotation parameters must be false.
        if (params.annotation_mirbase || params.annotation_srnaanno || params.annotation_pmiren) {
            log.error("If annotation_unitas is true, all other annotation " +
                "parameters (--annotation_mirbase, --annotation_srnaanno, " +
                "--annotation_pmiren) must be false.")
        }
    }
    
    // The parameters --from_counts and --only_preprocessing cannot be used
    // together
    if (params.from_counts && params.only_preprocessing) {
        // Add the ignored params to the list
        ignoreParams.addAll(['from_counts', 'only_preprocessing'])
        
        // Show the error
        log.error("The parameters 'from_counts' and 'only_preprocessing' are " +
            "mutually exclusive and cannot be used together. Please select " +
            "only one of them.")
    }

    // If counts_project_matrix is true, show the corresponding warning.
    if (params.counts_project_matrix){
        createProjectTableWarn()
    }

    // from_counts parameter warning
    if (params.from_counts) {

        // Get the CLI params
        def cliParams = workflow.commandLine.findAll(/--\S+/).collect { it.replaceFirst(/^--/, '') }

        // Parameters not compatible with '--from_counts'
        def previousStepsParams = [
            'extra_fastqc_args',
            'extra_multiqc_args',
            'trimming_adapters',
            'trimming_min_len',
            'trimming_max_len',
            'custom_fastp_args',
            'validation_depth',
            'filtering_db_mismatches',
            'filtering_genome_mismatches',
            'skip_qc_trim',
            'skip_fastqc',
            'skip_multiqc',
            'skip_filt_db',
            'skip_filt_genome'
        ]

        // Check if any of the parameters associated with previous steps of the
        // pipeline have been provided
        def filteredParams = cliParams.findAll { it in previousStepsParams }
        
        // Add the ignored params to the list
        ignoreParams.addAll(filteredParams)

        // If that condition is TRUE, show the corresponding warning...
        if (filteredParams && !ignoreParams.contains('from_counts')) {
            fromCountsWarn(filteredParams)
        }
    }

    // from_counts parameter warning
    if (params.only_preprocessing) {
        
        // Get the CLI params
        def cliParams = workflow.commandLine.findAll(/--\S+/).collect { it.replaceFirst(/^--/, '') }

        // Parameters not compatible with '--only_preprocessing'
        def nextStepsParams = [
            'skip_quantification',
            'skip_dea',
            'skip_annotation',
            'annotation_mismatches',
            'annotation_mirbase',
            'annotation_srnaanno',
            'annotation_pmiren',
            'annotation_unitas',
            'min_counts',
            'min_samples',
            'ea_p_value',
            'dea_alpha',
            'counts_project_matrix'
        ]

        // Check if any of the parameters associated with previous steps of the
        // pipeline have been provided
        def filteredParams = cliParams.findAll { it in nextStepsParams }

        // Add the ignored params to the list
        ignoreParams.addAll(filteredParams)
        
        // If that condition is TRUE, show the corresponding warning...
        if (filteredParams && !ignoreParams.contains('only_preprocessing')) {
            onlyPreprocessingWarn(filteredParams)
        }
    }

    // skip_dea parameter warning.
    if (params.skip_dea) {

        // Get the CLI params
        def cliParams = workflow.commandLine.findAll(/--\S+/).collect { it.replaceFirst(/^--/, '') }
        
        // Parameters not compatible with '--skip_dea'
        def nextStepsParams = [
            'annotation_mismatches',
            'annotation_mirbase',
            'annotation_srnaanno',
            'annotation_pmiren',
            'annotation_unitas'
        ]

        // Check if any of the parameters associated with previous steps of the
        // pipeline have been provided
        def filteredParams = cliParams.findAll { it in nextStepsParams }

        // Add the ignored params to the list
        ignoreParams.addAll(filteredParams)
        
        // If that condition is TRUE, show the corresponding warning...
        if (filteredParams && !ignoreParams.contains('skip_dea')) {
            skipDeaWarn(filteredParams)
        }
    }
}


//
// Print a warning if using '--counts_project_matrix'
//

def createProjectTableWarn() {
    log.warn """ '--counts_project_matrix' has been provided. Generating count tables
        at the project level may significantly slow down the pipeline if
        the number of samples in the project(s) is too high.
    """.stripIndent(true)
}


//
// Print a warning if using '--from_counts'
//

def fromCountsWarn(providedParams) {
    log.warn """ '--from_counts' has been provided. The pipeline will start at the
        quantification section. Any provided parameter related to the
        preprocessing section will be ignored:
        ${providedParams.collect { "--$it" }.join(', ')}
    """.stripIndent(true)
}


//
// Print a warning if using '--only_preprocessing'
//

def onlyPreprocessingWarn(providedParams) {
    log.warn """ '--only_preprocessing' has been provided. Only the preprocessing
        section of the pipeline will be executed. Any provided parameter
        related to other sections of the pipeline will be ignored:
        ${providedParams.collect { "--$it" }.join(', ')}
    """.stripIndent(true)
}


//
// Print a warning if using '--skip_dea'
//

def skipDeaWarn (providedParams){
    log.warn """ '--skip_dea' parameter has been provided. The pipeline steps after
        differential expression analysis (DEA) will neither be executed. Any
        provided parameter related to these steps will be ignored:
        ${providedParams.collect { "--$it" }.join(', ')}
    """.stripIndent(true)
}


//
// Verifies whether the input samples have an associated genome or if one can
// be assigned. The input is the content of the samplesheet, and the output is
// the same, but with the genome properly linked.
//

def validateAndAssignGenome (item) {
    
    // Get some fields from the input map
    def genome_path = item[4]
    def species_name = item[0]
    def predetermined_genome

    // If no genome is provided, look it up in the configuration file
    if (genome_path == []) {
        predetermined_genome = params.genomes.get(species_name, null)?.fasta
    } else {
        // Use the provided genome path
        predetermined_genome = genome_path
    }

    // If no valid genome found, throw an error
    if (predetermined_genome == null && !params.skip_filt_genome) {
        log.error("There is no genome associated with the following file:\n${item[2]}")
    }

    // Return (species, project, metadata_path, genome, file, groups)
    return [item[0], item[1], item[3], predetermined_genome, item[2], item[5]]
}


//
// Throw an error if the '--from_counts' parameter is provided but the input
// files are not in TSV format, or if TSV files are provided but the
// '--from_counts' parameter is missing.
//

def notTsvFilesError(files) {

    if (files.isEmpty() && params.from_counts) {
        // Throw an exception with the list of invalid files
        log.error("The --from_counts parameter must " +
            "be used only when counts matrix files in TSV " +
            "format are provided in the samplesheet.\n")
    } else if (!files.isEmpty() && !params.from_counts) {
        log.error("Counts matrices in TSV format " +
            "have been provided in the samplesheet, but the " +
            "--from_counts parameter has not been specified. Please " +
            "include this parameter to ensure correct processing.\n")
    }
}


//
// Throw an error if the 'from_counts' parameter is provided and any input file
// does not have an associated group. Additionally, print a warning if the input
// files have an associated group but the '--from_counts' parameter is not
// provided
//

def validateGroupInputUsage (item) {

    // Collect valid groups (non-empty elements at index 5)
    def groupsProvided = item.findAll { it != [] }
    
    // Check if all elements have valid groups
    if(groupsProvided.size() != item.size()) {
        if (params.from_counts) {
            log.error("Not all input files have an  " +
                "associated group, even though the --from_counts " +
                "parameter has been specified. Make sure to use the " +
                "--from_counts parameter and provide a group only when " +
                "the input file is a counts matrix.\n")
        }
    }

    if (groupsProvided.size() > 0 && !params.from_counts) {
        // Print a warning if '--from_counts' is false but groups were provided
        log.warn "Groups were detected in the input, but will " +
            "be ignored because '--from_counts' was not specified. Make " +
            "sure to use the --from_counts parameter and provide a group " + 
            "only when the input file is a counts matrix.\n"
    }
}


//
// This function is used to validate Accession list files.
//

def validateAccessionList(file_path) {

    // List of SRA patterns
    def sra_patterns = ['^SRR', '^ERR', '^DRR']
    def is_accession_list = true

    // Input file
    def input_file =  file(file_path)

    try {
        // Read lines
        input_file.eachLine { line ->

            // Check if the Runs match the SRA patterns
            def matches = sra_patterns.any { pattern ->
                line.trim().matches(pattern + '.*')
            }
            if (!matches) {
                is_accession_list = false
                return
            }
        }
        
        return is_accession_list

    } catch (Exception e) {
        return false
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// SEGUIR AQUI CON EL WORKFLOW PIPELINE_COMPLETION. CON EL QUE HAREMOS COSAS
// TRAS HABERSE EJECUTADO EL PIPELINE PRINCIPAL

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

// // Print help message, supply typical command line usage for the pipeline
// if (params.help) {
//    log.info paramsHelp("nextflow run mirplan.nf --input samplesheet.csv")
//    exit 0
// }

// // Validate input parameters
// validateParameters()

// // Print summary of supplied parameters
// log.info paramsSummaryLog(workflow)