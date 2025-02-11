#!/usr/bin/env nextflow

/*
========================================================================================
    QUALITYCONTROL Sub-Workflow
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { FASTQC } from "../../../modules/nf-core/fastqc"
include { MULTIQC } from "../../../modules/local/multiqc"

/*
========================================================================================
    Workflow QUALITYCONTROL
========================================================================================
*/

workflow QUALITY_CONTROL {
    take:
        ch_input            // channel: [[species:val(species), project:val(project)], path(file)]
        val_skip_multiqc    // boolean: true/false
        val_type            // value: Raw/Trimmed

    main:

        // Versions channels
        ch_fastqc_versions  = Channel.empty()
        ch_versions         = Channel.empty()

        // MULTIQC channels
        multiqc_report      = Channel.empty()

        // Add the type to the meta.id
        ch_input
            .map { meta, file ->
                def updatedMeta = meta.clone()
                updatedMeta.id = "${meta.id}_${val_type}"
                return [updatedMeta + [type: "${val_type}"], file]
            }
            .set {ch_input_adapted}

        // Execute FASTQC
        FASTQC(ch_input_adapted)

        // Save results in output channels
        ch_fastqc_versions = FASTQC.out.versions
        ch_versions        = ch_versions.mix(FASTQC.out.versions.first())

        // If val_skip_multiqc != true...
        if (!val_skip_multiqc) {
            
            // Update meta.id to be the project identifier
            FASTQC.out.zip 
                .map { meta, file ->
                    def updatedMeta = meta.clone()
                    updatedMeta.id = meta.project
                    return [updatedMeta, file]
            }.set {ch_fastqc_meta_updated}

            // Prepare input channel for MULTIQC
            // [[Glycine max, PRJNA720229], [/paht/file1.zip, /paht/file1.zip...] ]
            ch_fastqc_meta_updated
                .groupTuple(by:0, sort:true)
                .set {ch_fastqc_zip_by_project}

            // Execute MULTIQC
            MULTIQC(ch_fastqc_zip_by_project)

            // Save results in output channels
            multiqc_report = MULTIQC.out.report
            ch_versions    = ch_versions.mix(MULTIQC.out.versions.first())

        }

    emit:
        fastqc_zip = FASTQC.out.zip                                     // channel: [[species:val(species), project:val(project)], path(.zip file)]
        fastqc_html = FASTQC.out.html                                   // channel: [[species:val(species), project:val(project)], path(.html file)]
        multiqc_report = multiqc_report.ifEmpty(null)                   // channel: [ multiqc_report.html ]
        versions = ch_versions.ifEmpty(null)                            // channel: [ versions.yml ]
}
