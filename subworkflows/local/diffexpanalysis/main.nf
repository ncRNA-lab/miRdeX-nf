#!/usr/bin/env nextflow

/*
========================================================================================
    DIFFEXPANALYSIS Sub-Workflow
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { DEA; DEA_SUM } from "../../../modules/local/diffexpanalysis"

/*
========================================================================================
    Workflow DIFFEXPANALYSIS
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

workflow DIFFEXPANALYSIS {
    take:
        ch_input
        alpha
        min_counts
        min_samples

    main:

        // Perform exploratory and differential expression analyses.
        DEA(ch_input, alpha, min_counts, min_samples)

        // Create summary files for EA and DEA
        //DEA_SUM(DEA.out.deasum.collect(), DEA.out.easum.collect())

    emit:
        raw = DEA.out.raw
        sig = DEA.out.sig
        easum = DEA.out.easum
        deasum = DEA.out.deasum

}
