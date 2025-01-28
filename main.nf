#!/usr/bin/env nextflow

/*
========================================================================================
    miRPlan Nextflow Workflow
========================================================================================
    Github   :
    Contact  :
----------------------------------------------------------------------------------------
*/

// println """\
//          M I R P L A N - N F   P I P E L I N E
//          ===================================
//          genome       : ${params.genome}
//          reads        : ${params.reads}
//          outdir       : ${params.outdir}
//          """
//          .stripIndent()


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MIRPLAN                   } from './workflows/mirplan'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_mirplan_pipeline'


//
// WORKFLOW: Run main nf-mirplan pipeline
//
workflow MAIN_MIRPLAN {

    main:
    
    // Create an empty channel for versions 
    ch_versions = Channel.empty()

    //
    // WORKFLOW: Run mirplan workflow
    //

    // Check if the input file exists
    ch_samplesheet = Channel.value(file(params.input, checkIfExists: true))
    
    // Run the workflow
    MIRPLAN ()

    // // Get the versions channel from the main workflow
    // ch_versions = ch_versions.mix(MIRPLAN.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
//
workflow {

    main:
        println(params.version)
        
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
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
