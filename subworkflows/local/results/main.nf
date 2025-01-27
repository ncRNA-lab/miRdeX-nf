#!/usr/bin/env nextflow

/*
========================================================================================
    RESULTS_GENERATION Sub-Workflow
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { PREABS_MATRIX } from "../../../modules/local/matrix"

/*
========================================================================================
    Workflow RESULTS_GENERATION
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

workflow RESULTS_GENERATION {
    take:
        dea_sum
        annot_table

    main:

        PREABS_MATRIX(
            dea_sum,
            annot_table
        )

    emit:
        preabs = PREABS_MATRIX.out.preabs
}
