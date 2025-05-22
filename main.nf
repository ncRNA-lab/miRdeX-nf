#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    antoglz/mirplan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/antoglz/miRPlan-nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MIRPLAN                   } from './workflows/mirplan'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_mirplan_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_mirplan_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAIN_MIRPLAN {

    main:
    
    // Create an empty channel for versions 
    ch_versions = Channel.empty()
    
    // Run the workflow
    MIRPLAN (params.input, ch_versions)

    // Get the versions channel from the main workflow
    ch_versions = ch_versions.mix(MIRPLAN.out.versions)

    emit:
    summary  = MIRPLAN.out.summary
    versions = MIRPLAN.out.versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
        
        //
        // SUBWORKFLOW: Run initialisation tasks
        //
        PIPELINE_INITIALISATION(
            params.version,
            params.validate_params,
            args,
            params.outdir
        )

        //
        // SUBWORKFLOW: Run the main workflow
        //
        MAIN_MIRPLAN ()
    
        //
        // SUBWORKFLOW: Run completion tasks
        //
        PIPELINE_COMPLETION(
            MAIN_MIRPLAN.out.summary,
            MAIN_MIRPLAN.out.versions,
            params.email,
            params.email_on_fail,
            params.plaintext_email,
            params.outdir
        )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
