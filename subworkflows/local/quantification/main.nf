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
        not_in_memory           // value: true or false

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

        // create a channel at the group level
        ch_counts
            .flatMap { meta, file ->
                meta.groups.collect { group -> 
                    def new_meta = meta.clone()
                    new_meta.id = group
                    new_meta.group_id = group.split('_')[-1]
                    new_meta.remove('groups')
                    [new_meta, file]
                }
            }
            .groupTuple(by:0)
            .map{ meta, files -> [meta, files, meta.metadata]}
            .set{ ch_counts_by_group }
            
        // Create the count matrix
        COUNTS_MATRIX(ch_counts_by_group, type, not_in_memory)

        // Add the software version
        ch_versions = ch_versions.mix(COUNTS_MATRIX.out.versions)

    emit:
        group_matrix    = COUNTS_MATRIX.out.matrix
        versions        = ch_versions      // channel: [ path(versions.yml) ]
}
