#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    antoglz/miRdeX-nf
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/antoglz/miRdeX-nf
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MIRDEX                    } from './workflows/mirdex'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_mirdex_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_mirdex_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAIN_MIRDEX {

    main:
    
    // Create an empty channel for versions 
    ch_versions = Channel.empty()
    
    // Run the workflow
    MIRDEX (params.input, ch_versions)

    // Get the versions channel from the main workflow
    ch_versions = ch_versions.mix(MIRDEX.out.versions)

    emit:
    summary  = MIRDEX.out.summary
    versions = MIRDEX.out.versions
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
            params.outdir
        )

        //
        // SUBWORKFLOW: Run the main workflow
        //
        MAIN_MIRDEX ()
    
        //
        // SUBWORKFLOW: Run completion tasks
        //
        PIPELINE_COMPLETION(
            MAIN_MIRDEX.out.summary,
            MAIN_MIRDEX.out.versions,
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
