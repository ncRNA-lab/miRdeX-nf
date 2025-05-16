//
// Subworkflow with functionality specific to the mirnote-nf pipeline. This
// subworkflow is based on the utils_nfcore_rnaseq_pipeline subworkflow from
// the nf-core/rnaseq pipeline (https://github.com/nf-core/rnaseq), with minor
// modifications to adapt it for the miRNote-nf pipeline.
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'
include { getWorkflowVersion        } from '../../nf-core/utils_nfcore_pipeline'


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
    
    // List of numeric params
    def numericParams = [
        'substitutions'      : params.substitutions,
        'three_add'          : params.three_add,
        'five_add'           : params.five_add,
        'ends_modification'  : params.ends_modification,
        'rpm'                : params.rpm,
        'raw_counts'         : params.counts
    ]

    // Check if the params provided are positive numbers
    numericParams.each { name, value ->
        if (value != null && value < 0) {
            log.error("The parameter '${name}' must be a non-negative number. Current value: ${value}")
        }
    }

    // Check if both counts and rpm params have been used (not  allowed)
    if (params.counts > 0 && params.rpm > 0) {
        // Show the error
        log.error("The parameters 'counts' and 'rpm' are " +
            "mutually exclusive and cannot be used together. Please select " +
            "only one of them.")
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
        def output_tf = new File(workflow.launchDir.toString(), ".pipeline_report.txt")
        output_tf.withWriter { w -> w << rendered }
        nextflow.extension.FilesEx.copyTo(output_tf.toPath(), "${outdir}/pipeline_info/pipeline_report.txt")
        output_tf.delete()
    }
}