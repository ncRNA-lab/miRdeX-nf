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
include { COUNTS_MATRIX as COUNTS_MATRIX_RAW } from "../../../modules/local/counts_matrix"
include { COUNTS_MATRIX as COUNTS_MATRIX_RPM } from "../../../modules/local/counts_matrix"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Build the group-level channel from an input channel [meta, file]
def groupSamplesByGroup(ch_counts) {

    // 1) Build a per-group metadata lookup channel: [group, meta_map]
    def ch_metadata = ch_counts
        .flatMap { meta, _file ->
            meta.groups.collect { group ->
                def group_id = group.split('_')[-1]
                [ group,
                  [ id: group,
                    project: meta.project,
                    species: meta.species,
                    metadata: meta.metadata,
                    group_id: group_id ] ]
            }
        }
        .unique()

    // 2) Group files by group and attach metadata → [meta_map, files, metadata_path]
    def ch_counts_by_group = ch_counts
        .flatMap { meta, file ->
            meta.groups.collect { group -> [group, file] }
        }
        .groupTuple(by: 0)
        .combine(ch_metadata, by: 0)
        .map { _group, files, meta_map ->
            [ meta_map, files, meta_map.metadata ]
        }

    return ch_counts_by_group
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW QUANTIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUANTIFICATION {

    take:
        ch_input                // channel: [ meta, path(fastq.gz) ]
        calculate_rpm           // value: true or false
        not_in_memory           // value: true or false


    main:

        // Create required empty channels
        ch_versions             = Channel.empty()
        ch_rpm_matrix           = Channel.empty()

        // Calculate the raw counts
        COUNTS(ch_input)
    
        // Calculate the RPM (Optional)
        if (calculate_rpm) {
            // Calculate sequences RPM per sample
            RPM(COUNTS.out.raw)
            
            // Group samples by analysis group
            ch_rpm_by_group = groupSamplesByGroup(RPM.out.rpm)

            // Create the count matrix
            COUNTS_MATRIX_RPM(ch_rpm_by_group, 'rpm', not_in_memory)

            // Save the results into the output channel
            ch_rpm_matrix = COUNTS_MATRIX_RPM.out.matrix
        }

        // Group samples by analysis group
        ch_raw_by_group = groupSamplesByGroup(COUNTS.out.raw)

        // Create the count matrix
        COUNTS_MATRIX_RAW(ch_raw_by_group, 'raw', not_in_memory)

        // Add the software version
        ch_versions = ch_versions.mix(COUNTS_MATRIX_RAW.out.versions)

    emit:
        raw_matrix      = COUNTS_MATRIX_RAW.out.matrix  // channel: [ meta, path(tsv) ]
        rpm_matrix      = ch_rpm_matrix                 // channel: [ meta, path(tsv) ]
        versions        = ch_versions                   // channel: [ path(versions.yml) ]
}
