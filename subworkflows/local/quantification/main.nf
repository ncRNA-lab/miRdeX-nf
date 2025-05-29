/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { COUNTS        } from "../../../modules/local/counts"
include { RPM           } from "../../../modules/local/rpm"
include { COUNTS_MATRIX } from "../../../modules/local/counts_matrix"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW QUANTIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUANTIFICATION {

    take:
        ch_input
        type                    // value: 'raw' or 'rpm'
        counts_project_matrix   // value: true or false

    main:

        // Create empty channel for versions
        ch_versions             = Channel.empty()

        // Calculate the raw counts
        COUNTS(ch_input)

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
                def updatedMeta = meta.clone()
                updatedMeta.id = updatedMeta.project
                return [updatedMeta.id, updatedMeta, file]
            }
            .groupTuple(by:[0,1])
            .map { it -> [it[1], it[2], it[1].metadata, it[1].valid_groups] }
            
        // Create the count matrix
        COUNTS_MATRIX(ch_counts_by_project, type, counts_project_matrix)

        // Add the software version
        ch_versions = ch_versions.mix(COUNTS_MATRIX.out.versions)

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
        group_matrix    = ch_counts_matrix
        versions        = ch_versions      // channel: [ path(versions.yml) ]
}
