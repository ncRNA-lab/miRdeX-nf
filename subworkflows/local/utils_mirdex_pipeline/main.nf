//
// Subworkflow with functionality specific to the miRdeX-nf pipeline. This
// subworkflow is based on the utils_nfcore_rnaseq_pipeline subworkflow from
// the nf-core/rnaseq pipeline (https://github.com/nf-core/rnaseq), with minor
// modifications to adapt it for the miRdeX-nf pipeline.
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { softwareVersionsToYAML    } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../nf-core/utils_nextflow_pipeline'    

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {
    
    take:
    version             // boolean: Display version and exit
    validate_params     // boolean: Boolean whether to validate parameters against the schema at runtime
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
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow PIPELINE_COMPLETION {

    take:
    summary
    versions
    email               //  string: email address
    email_on_fail       //  string: email address sent on pipeline failure
    plaintext_email     // boolean: Send plain-text email instead of HTML
    outdir              //    path: Path to output directory where results will be published

    main:

    // Get the parameters from the pipeline
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    
    // Completion email
    workflow.onComplete {
        if (email || email_on_fail) {
            sendCompletionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir
            )
        }
    }

    workflow.onError {
        log.error "Pipeline failed."
    }
    
    // Save a file with the software versions
    softwareVersionsToYAML(versions)
        .collectFile(storeDir: "${params.outdir}/02-Results", name: "mirdex_versions_${workflow.nextflow.timestamp.replaceAll(/[\\s:]/, '_')}.yml", sort: true, newLine: true)
    
    // Sort the summary chanel
    summary
        .toSortedList { a, b -> 
            a.species <=> b.species ?: a.group <=> b.group ?: a.comparison_id <=> b.comparison_id ?: a.sample <=> b.sample
        }
        .flatMap()
        .set{ch_sorted_summary}

    // Save a file with the summary information
    summaryToTsv(ch_sorted_summary)
        .collectFile(storeDir: "${params.outdir}/02-Results", name: "mirdex_summary_${workflow.nextflow.timestamp.replaceAll(/[\\s:]/, '_')}.tsv", newLine: true, sort: false)
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

    // Validate annotation parameters
    if (!params.annotation_unitas) {

        // Validate 'databases' parameter
        if (params.databases) {
            def validNames = ['mirbase', 'pmiren', 'srnaanno'] as Set

            // Comprobar espacios
            if (params.databases.contains(' ')) {
                throw new IllegalArgumentException("The 'databases' parameter must not contain spaces. Use comma-separated names without spaces.")
            }

            def dbList = params.databases
                .toLowerCase()
                .split(',')
                .collect { it.trim() }

            // Comprobar que no hay duplicados
            if (dbList.size() != dbList.toSet().size()) {
                throw new IllegalArgumentException("Duplicate entries found in 'databases' parameter: ${params.databases}")
            }

            // Comprobar que todos los nombres sean válidos
            def unknown = dbList.findAll { !validNames.contains(it) }
            if (!unknown.isEmpty()) {
                throw new IllegalArgumentException("Invalid database names found in 'databases' parameter: ${unknown.join(', ')}. Allowed values are: miRBase, PmiREN, sRNAanno.")
            }
        }
        
        // List of numeric params
        def numericParams = [
            'log2fc_threshold'            : params.log2fc_threshold,
            'substitutions'               : params.substitutions,
            'three_add'                   : params.three_add,
            'five_add'                    : params.five_add,
            'ends_modification'           : params.ends_modification,
            'min_counts_filt'             : params.min_counts_filt,
            'min_rpm_filt'                : params.min_rpm_filt,
            'min_relative_abundance_filt' : params.min_relative_abundance_filt

        ]

        // Check if the params provided are positive numbers
        numericParams.each { name, value ->
            if (value != null && value < 0) {
                throw new IllegalArgumentException("The parameter '${name}' must be a non-negative number. Current value: ${value}")
            }
        }

        // Create a list with the three filter parameters
        def active_filters = [
            params.min_counts_filt,
            params.min_rpm_filt,
            params.min_relative_abundance_filt
        ]
        // Count how many filters are active (value > 0)
        .count { it > 0 }

        // If more than one filter is active, throw an error and stop the pipeline
        if (active_filters > 1) {
            throw new IllegalArgumentException("Only one abundance filter can be used at a time: " +
                "'min_counts_filt', 'min_rpm_filt' or 'min_relative_abundance_filt'. ")
        }
    }
    
    // The parameters --from_counts and --only_preprocessing cannot be used
    // together
    if (params.from_counts && params.only_preprocessing) {
        // Add the ignored params to the list
        ignoreParams.addAll(['from_counts', 'only_preprocessing'])
        
        // Show the error
        throw new IllegalArgumentException("The parameters 'from_counts' and 'only_preprocessing' are " +
            "mutually exclusive and cannot be used together. Please select " +
            "only one of them.")
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
            'trimming_custom_args',
            'validation_depth',
            'filtering_db_mismatches',
            'filtering_genome_mismatches',
            'skip_qc_trim',
            'skip_fastqc',
            'skip_multiqc',
            'filt_db',
            'filt_genome'
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
            'skip_annotation',
            'databases',
            'mirna_classes',
            'substitutions',
            'five_add',
            'three_add',
            'ends_modification',
            'annotation_unitas',
            'bowtie_nti_ext_args',
            'counts',
            'rpm',
            'min_counts',
            'min_samples',
            'ea_p_value',
            'dea_alpha',
            'counts_not_memory',
            'global_matrix',
            'global_fields'
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

    // skip_annotation parameter warning.
    if (params.skip_annotation) {

        // Get the CLI params
        def cliParams = workflow.commandLine.findAll(/--\S+/).collect { it.replaceFirst(/^--/, '') }

        // Parameters not compatible with '--skip_annotation'
        def nextStepsParams = [
            'databases',
            'mirna_classes',
            'substitutions',
            'five_add',
            'three_add',
            'ends_modification',
            'annotation_unitas',
            'bowtie_nti_ext_args',
            'counts',
            'rpm',
            'counts_not_memory',
            'global_matrix'
        ]
        
        // Check if any of the parameters associated with previous steps of the
        // pipeline have been provided
        def filteredParams = cliParams.findAll { it in nextStepsParams }

        // Add the ignored params to the list
        ignoreParams.addAll(filteredParams)
        
        // If that condition is TRUE, show the corresponding warning...
        if (filteredParams && !ignoreParams.contains('skip_annotation')) {
            skipAnnotationWarn(filteredParams)
        }
    }
}

//
// It organises the input channel to the workflow.
//

def organiseInputChannel(input_channel, key) {

    // Get the species and project from metadata
    ch_sp_pj = input_channel
        .map { item -> item[2] } // metadata file
        .unique()
        .splitCsv(header: true, sep: '\t')
        .map { meta -> [meta[key], meta['Species'], meta['Project']] }
        .unique()

    // Add this information to the original channel
    return input_channel
        .combine(ch_sp_pj, by:0)
        .map{ id, file, meta, genome, group, sp, project ->
            [[id: id, species: sp, project: project, metadata: meta, genome: genome, single_end: true, group_id:group], file]
        }
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
// Print a warning if using '--skip_annotation'
//

def skipAnnotationWarn(providedParams) {
    log.warn """'--skip_annotation' has been provided. The annotation step and all downstream 
        processes that depend on it will be skipped. The following parameters will be ignored:
        ${providedParams.collect { "--$it" }.join(', ')}
    """.stripIndent(true)
}


//
// Verifies whether the input samples have an associated genome or if one can
// be assigned. The input is the content of the samplesheet, and the output is
// the same, but with the genome properly linked.
//

def validateAndAssignGenome (species, genome) {
    
    // Get some fields from the input map
    def predetermined_genome

    // If no genome is provided, look it up in the configuration file
    if (genome == []) {
        predetermined_genome = params.genomes.get(species, null)?.fasta
    } else {
        // Use the provided genome path
        predetermined_genome = genome
    }

    // If no valid genome found, throw an error
    if (predetermined_genome == null && (params.filt_genome || params.mirna_classes != 'ref_miRNA')) {
        throw new IllegalArgumentException("There is no genome associated with the following species:\n${species}")
    }

    // Return (species, project, metadata_path, genome, file, groups)
    return predetermined_genome
}


//
// Throw an error if the '--from_counts' parameter is provided but the input
// files are not in TSV format, or if TSV files are provided but the
// '--from_counts' parameter is missing.
//

def notTsvFilesError(files) {

    if (files.isEmpty() && params.from_counts) {
        // Throw an exception with the list of invalid files
        throw new IllegalArgumentException("The --from_counts parameter must " +
            "be used only when counts matrix files in TSV " +
            "format are provided in the samplesheet.\n")
    } else if (!files.isEmpty() && !params.from_counts) {
        throw new IllegalArgumentException("Counts matrices in TSV format " +
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
            throw new IllegalArgumentException("Not all input files have an  " +
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

    } catch (Exception _e) {
        return false
    }

}

//
// Separate the paths to the mature and precursor miRNA databases
// into distinct elements within the same channel.
//

def formatDbChannel(ch, id_name) {
    return ch.flatMap { list ->
        def (mature, precursor) = list
        return [
            [ [id: "${id_name}_mature"], mature ],
            [ [id: "${id_name}_precursor"], precursor ]
        ]
    }
}

//
// Prepare the summary channel so that its content is written to a TSV file.
//

def summaryToTsv(ch) {
    return ch
        .collect()
        .map { rows -> 
            def keys = rows[0].keySet()
            def header = keys.join('\t')
            def lines = rows.collect { row -> 
                keys.collect { k -> row[k] }.join('\t')
            }
            return ([header] + lines)
        }
        .flatten()
}

//
// Construct and send completion email (based on CompletionEmail function from
// Nf-core)
//

def sendCompletionEmail(summary_params, email, email_on_fail, plaintext_email, outdir) {

    // Construir asunto del correo
    def subject = "[${workflow.manifest.name}] Successful: ${workflow.runName}"
    if (!workflow.success) {
        subject = "[${workflow.manifest.name}] FAILED: ${workflow.runName}"
    }

    // Procesar resumen en un solo mapa
    def summary = [:]
    summary_params
        .keySet()
        .sort()
        .each { group ->
            summary << summary_params[group]
        }

    // Campos adicionales del workflow
    def misc_fields = [:]
    misc_fields['Date Started']              = workflow.start
    misc_fields['Date Completed']            = workflow.complete
    misc_fields['Pipeline script file path'] = workflow.scriptFile
    misc_fields['Pipeline script hash ID']   = workflow.scriptId
    if (workflow.repository) {
        misc_fields['Pipeline repository Git URL']    = workflow.repository
    }
    if (workflow.commitId) {
        misc_fields['Pipeline repository Git Commit'] = workflow.commitId
    }
    if (workflow.revision) {
        misc_fields['Pipeline Git branch/tag']        = workflow.revision
    }
    misc_fields['Nextflow Version']           = workflow.nextflow.version
    misc_fields['Nextflow Build']             = workflow.nextflow.build
    misc_fields['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Variables para la plantilla
    def email_fields = [
        version      : getWorkflowVersion(),
        runName      : workflow.runName,
        success      : workflow.success,
        dateComplete : workflow.complete,
        duration     : workflow.duration,
        exitStatus   : workflow.exitStatus,
        errorMessage : workflow.errorMessage ?: 'None',
        errorReport  : workflow.errorReport ?: 'None',
        commandLine  : workflow.commandLine,
        projectDir   : workflow.projectDir,
        summary      : summary << misc_fields
    ]

    // Determinar si enviar correo y a quién
    def email_address = email
    if (!email && email_on_fail && !workflow.success) {
        email_address = email_on_fail
    }

    // No hacer nada si no hay destinatario
    if (!email_address) return

    // Render the TXT template
    def tf           = new File("${workflow.projectDir}/assets/email_template.txt")
    def templateText = tf.text
    def engine       = new groovy.text.SimpleTemplateEngine()
    def email_txt    = engine.createTemplate(templateText).make(email_fields).toString()

    // Render the HTML template
    def hf            = new File("${workflow.projectDir}/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def rendered      = html_template.toString()

    // Check if an email address was provided
    if (email_address) {

        // Check if the mail must be sent as plaintext
        if (plaintext_email) {
            rendered = email_txt
        }

        // Send the email
        try {
            sendMail(
                to: email_address,
                subject: subject,
                body: rendered
            )
        }
        catch (Exception all) {
            log.error("Failed to send email: ${all.message}")
            throw new RuntimeException("No se pudo enviar el correo: ${all.message}", all)
        }
        
        // Save a copy of the email to the output directory
        def output_tf = new File(workflow.launchDir.toString(), ".pipeline_report_${workflow.nextflow.timestamp}.txt")
        output_tf.withWriter { w -> w << rendered }
        nextflow.extension.FilesEx.copyTo(output_tf.toPath(), "${outdir}/pipeline_info/pipeline_report_${workflow.nextflow.timestamp}.txt")
        output_tf.delete()
    }
}

/**
* Filters a channel of DEA results based on a Mann–Whitney–Wilcoxon (MWW)
* p-value threshold.
*
* This function receives a channel of tuples [meta, dea_file, ea_file], where:
*   - `meta` contains sample metadata, including a comparison `id`
*   - `dea_file` is not used in the filtering but kept for output
*   - `ea_file` is a tab-delimited file containing statistical results
*
* It performs the following steps:
*   - Parses `ea_file` to extract the p-value in column 9 (index 8)
*   - Retains only rows with MWW p-value less than or equal to the threshold
*   - Matches valid comparisons back to the original input using a shared ID
*   - Returns a channel of filtered tuples: [meta, dea_file, ea_file]
*
* The `group_id` used for matching is extracted from `meta.id`, removing any
* trailing "_number" suffix (e.g., "comparison_1" becomes "comparison").
*
* @param ch_input   Channel with tuples: [meta, dea_file, ea_file]
* @param threshold  Maximum MWW p-value to consider a comparison valid
* @return           Filtered channel with [meta, dea_file, ea_file] tuples
*/

def filterByMwwPvalue(ch_input, threshold) {

    // Filter comparisons using the threshold
    def valid_mww_comparisons = ch_input
        .map{ _meta, _dea_file, ea_file -> ea_file }
        .splitCsv(sep: '\t', skip: 1)
        .map { item ->
            def mww_p_value = item[8].toDouble()
            if (mww_p_value <= threshold) {
                return [item[1], item[8]]
            }
        }

    // Remove those comparisons that do not meet the threshold
    def ch_valid_comparisons = ch_input
        .map{ meta, dea_file, ea_file ->
            def group_id = meta.id.replaceAll(/_\d+$/, '')
            [group_id, meta, dea_file, ea_file]
        }
        .join(valid_mww_comparisons)
        .map{item -> [item[1], item[2], item[3]]}
    
    return ch_valid_comparisons
}


/**
* Writes a sample sheet (CSV format) listing metadata and file paths for a
* set of processed samples.
*
* This function takes a Nextflow channel containing tuples of [meta, file],
* where:
*   - `meta` is a map with sample metadata (species, project, metadata path,
*     genome, group ID)
*   - `file` is the result file generated by a process (usually in work/)
*
* It reconstructs the final published path of each file using `publishDir`,
* assuming a structure:
*   ${publishDir}/${species}/${project}/${filename}
*
* Then it:
*   - Groups rows by [species, project]
*   - Writes them into a CSV file with header:
*     Species, Project, File, Metadata, Genome, Group
*   - Merges with existing file (if any) and removes duplicates
*
* @param ch_to_write   Channel with tuples: [meta, file]
* @param publishDir    Directory where files were published (not from work/)
* @param output_path   Path to the final output samplesheet CSV file
*/

def writeSampleSheet(ch_to_write, publishDir, output_path) {

    ch_to_write
        .map { meta, file ->
            def filename = file.getName()
            def published_path = "${publishDir}/${meta.species.replaceAll(' ', '_')}/${meta.project}/${filename}"
            [
                meta.species.toString(),
                meta.project.toString(),
                published_path.toString(),
                meta.metadata.toString(),
                meta.genome.toString(),
                meta.group_id.toString()
            ]
        }
        .groupTuple(by: [0,1])
        .subscribe { item ->
            def species = item[0]
            def project = item[1]
            def records = item.drop(2)

            def filePath = new File(output_path).parent
            new File(filePath).mkdirs()

            def header = "Species,Project,File,Metadata,Genome,Group\n"
            def fileContent = []

            def numRows = records[0].size()
            fileContent = (0..<numRows).collect { i ->
                ([species, project] + records.collect { it[i] }).join(',')
            }

            // Añadir contenido existente si el archivo ya existe
            if (new File(output_path).exists()) {
                def existingContent = new File(output_path).text.readLines().drop(1)
                fileContent = existingContent + fileContent
            }

            fileContent = new LinkedHashSet(fileContent).toList()
            def fileContentString = fileContent.join('\n')

            new File(output_path).withWriter('UTF-8') { writer ->
                writer << header
                writer << fileContentString + '\n'
            }
        }
}