/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { LIBRARIES_VALIDATION } from "../../../modules/local/libraries_validation"
include { COUNTS_VALIDATION    } from "../../../modules/local/counts_validation"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW VALIDATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VALIDATION {
    take:
        ch_input                    // channel: [[id:(counts_id/project_id), metadata:path(metadata)], counts_file/[lib_file1, lib_file2...]]
        type                        // value: 'libraries' or 'counts'
        replicates_threshold        // integer: > 2
        depth_threshold             // integer: > 0 (It is not used when type = 'counts')

    main:

        // Create empty channel for versions
        ch_versions             = Channel.empty()

        // Branch the workflow based on the value of "type"
        if (type == 'libraries') {

            // Extract the metadata from the map object.
            ch_input
                .map{ meta, file ->
                    return [meta, file, meta.metadata]
                }
                .set { ch_input_validation_lib }

            // Validate the library depth and the number of project replicates
            LIBRARIES_VALIDATION(ch_input_validation_lib, depth_threshold, replicates_threshold)

            // Create the libraries summary channel
            LIBRARIES_VALIDATION.out.sumlibraries
                .splitCsv(sep: '\t' )
                .map { item ->
                    def run = item[0].replaceFirst(/(\.fastp)?\.fastq(\.gz)?$/, '')
                    return [run, item[1], item[2], item [3]]
                }
                .set { ch_libraries_summary }

            // Create the projects summary channel
            LIBRARIES_VALIDATION.out.sumprojects
                .splitCsv(sep: '\t' )
                .map{ item ->
                    // Create the subproject name
                    def group_name = "${item[0]}_${item[1]}"

                    // Create a map with all the project information
                    [
                        project: item[0],
                        group: group_name,
                        group_id:item[1],
                        num_valid_samples: item[2],
                        num_notvalid_samples: item[3], 
                        validity: item[4]
                    ]
                }
                .tap { ch_projects }
                // Create a list with valid subprojects and another with non-valid
                // groups
                .map{ item -> [item.project, item]}
                .groupTuple(by:0)
                .map { project, groups ->
                    def valid_groups = groups.findAll { it.validity == 'valid' }*.group_id
                    def notvalid_groups = groups.findAll { it.validity == 'not-valid' }*.group_id
                    [project, [valid_groups: valid_groups, notvalid_groups: notvalid_groups]]
                }
                .set{ ch_projects_to_lib }

            LIBRARIES_VALIDATION.out.sumprojects.view()
            ch_libraries_summary.view()

            // Combine the valid and invalid libraries into the same channel.
            LIBRARIES_VALIDATION.out.valid
                .concat(LIBRARIES_VALIDATION.out.notvalid)
                .flatMap { meta, files ->
                    def fileList = files instanceof List ? files : [files]
                    fileList.collect { [meta, it] }
                }
                .map { meta, file ->
                    def updatedMeta = meta.clone()
                    def idParent = meta.id
                    updatedMeta.id = file.getName().replaceFirst(/\.(valid|notvalid)\.(fastq(\.gz)?|fq(\.gz)?)$/, '')
                    return [updatedMeta.id, updatedMeta, file, idParent]
                }
                // Add the summary_libraries_ch information to the metadata of the main channel.
                .combine(ch_libraries_summary, by: 0)
                .map { _id, meta, file, idparent, depth, valid1, valid2 ->
                    [idparent, meta + [depth: "${depth}", depth_validity: "${valid1}", replicates_validity: "${valid2}"], file]
                }
                // Add the lists of valid and non-valid subprojects to the meta.
                .combine(ch_projects_to_lib, by: 0)
                .map { _id, meta, file, lists ->
                    [meta + [ 
                        valid_groups: lists.valid_groups,
                        notvalid_groups: lists.notvalid_groups
                    ], file]
                }
                .set { ch_files }

            // Save the software version
            ch_versions = ch_versions.mix(LIBRARIES_VALIDATION.out.versions)

        } else if (type == 'counts'){

            // Extract the metadata and the group_id from the map object.
            ch_input
                .map{ meta, file ->
                    return [meta, file, meta.metadata, meta.group_id]
                }
                .set { ch_input_validation_counts }
            
            // Validate the counts matrices
            COUNTS_VALIDATION(ch_input_validation_counts, replicates_threshold)
            
            // Create a channel with the project info
            COUNTS_VALIDATION.out.summary
                .splitCsv(sep:'\t')
                .map{ item ->
                    // Create a map with all the project information
                    [
                        project: item[1],
                        group: item[0],
                        group_id: item[2],
                        num_valid_samples: item[5],
                        num_notvalid_samples: item[6], 
                        validity: item[3],
                        cause: item[4]
                    ]
                }
                .set{ch_projects}
            
            // Create an intermediate channel with the project info and the group as id
            ch_projects
                .map{item -> [item.group, item]}
                .set { ch_projects_to_file }

            // Add the project info to the output files channel
            COUNTS_VALIDATION.out.valid
                .concat(COUNTS_VALIDATION.out.notvalid)
                .map { meta, file ->
                    def updatedMeta = meta.clone()
                    updatedMeta.id = file.getName().replaceFirst(/\.(valid|notvalid)\.tsv$/, '')
                    return [updatedMeta.id, updatedMeta, file]
                }
                .combine(ch_projects_to_file, by: 0)
                .map{ _id, meta, file, meta_group ->
                    [ meta +
                        [
                        depth: "NULL",
                        depth_validity: "NULL",
                        replicates_validity: meta_group.validity,
                        group:meta_group.group,
                        group_id: meta_group.group_id,
                        valid_groups: "NULL",
                        notvalid_groups: "NULL",
                        num_valid_samples: meta_group.num_valid_samples,
                        num_notvalid_samples:meta_group.num_notvalid_samples
                        ],
                        file]
                }
                .set{ ch_files }
        }

    emit:
        files = ch_files            // channel: [ val(meta), path(fastq) ]
        projects = ch_projects      // channel: [ val(summary) ]
        versions = ch_versions      // channel: [ path(versions.yml) ]
}
