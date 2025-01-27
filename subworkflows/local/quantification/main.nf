#!/usr/bin/env nextflow

/*
========================================================================================
    QUANTIFICATION Sub-Workflow
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { COUNTS; COUNTS_MATRIX; RPM } from "../../../modules/local/counts"

/*
========================================================================================
    Workflow QUANTIFICATION
========================================================================================
*/

workflow QUANTIFICATION {

    take:
        libraries

    main:

        // Calculate the absolute counts
        COUNTS(libraries)

        // Change the meta.id from file to project.
        COUNTS.out.abs
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = meta.project
                tuple(meta.species, meta.project, updatedMeta, file)
            }
            .groupTuple(by: [0,1])
            .map{ item ->
                [item[2][0], item[3]]
            }
            .set{ch_abs_counts}

        // ch_abs_counts.view()
        // Calculate the RPM
        //if (params.norm_rpm) {
        //    RPM(COUNTS.out.abs)
        //}

        // Create the count matrix
        COUNTS_MATRIX(ch_abs_counts)

    emit:
        subproject = COUNTS_MATRIX.out.matrix
}
