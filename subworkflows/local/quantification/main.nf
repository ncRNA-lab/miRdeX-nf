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

include { COUNTS        } from "../../../modules/local/counts"
include { RPM           } from "../../../modules/local/rpm"
include { COUNTS_MATRIX } from "../../../modules/local/counts_matrix"

/*
========================================================================================
    Workflow QUANTIFICATION
========================================================================================
*/

workflow QUANTIFICATION {

    take:
        libraries
        type            // value: 'raw' or 'rpm'

    main:

        // Calculate the raw counts
        COUNTS(libraries)

        // Change the meta.id from file to project.
        ch_counts = COUNTS.out.raw
    
        // Calculate the RPM
        if (type == 'rpm') {
           RPM(COUNTS.out.raw)
           ch_counts = RPM.out.rpm
        }
 
        // Change the meta.id from file to project.
        ch_counts_by_project = ch_counts
            .map { meta, file -> 
                def updatedMeta = meta + [ id: meta.project ]
                return [updatedMeta.id, updatedMeta, file]
            }
            .groupTuple(by: [0,1])
            .map { it -> [it[1], it[2]] }

        ch_counts_by_project.view()

        // Create the count matrix
        COUNTS_MATRIX(ch_counts_by_project, type)

        // Create a new ID and set the group_id
        COUNTS_MATRIX.out.matrix
            .map{ meta, file ->

                // Get the filename
                //def fileName = file.toString().split('/').last().replace('.counts.tsv', '')
                def fileName = file.toString().split('/').last().replaceAll(/\.raw\.tsv|\.rpm\.tsv$/, '')

                // Create the id and group_id using the filename
                def group_id = fileName.split('_')[1]

                // Assign this values to the corresponding fields of the map
                meta.id = fileName
                meta.group_id = group_id
                return [meta, file]
            }
            .set{ ch_counts_matrix }

    emit:
        group_matrix = ch_counts_matrix
}
