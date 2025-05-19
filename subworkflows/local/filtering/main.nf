#!/usr/bin/env nextflow

/*
========================================================================================
    FILTERING Sub-Workflow
========================================================================================
*/

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { BOWTIE_BUILD  } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_DB } from '../../../modules/local/bowtie/align'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_GENOME } from '../../../modules/local/bowtie/align'

/*
========================================================================================
    Workflow FILTERING
========================================================================================
*/

workflow FILTERING {
    take:
        input                   // channel: [[id:val(id), project:val(projec), species:val(species), genome:val(genome)], file]
        ref                     // value: 'database' or 'genome'
        database                // path: reference fasta file (optional)
        ch_versions             // channel: [ path(versions.yml) ]

    main:
    
        // Create a channel for the outputs file
        ch_bowtie_bam       = Channel.empty()
        ch_bowtie_aligned   = Channel.empty()
        ch_bowtie_unaligned = Channel.empty()

        // Branch the workflow based on the value of "ref"
        if (ref == 'database') {

            ch_input_build = Channel.value([[id:'Filtering database'], database])
            
            // Index the unique reference file (database)
            BOWTIE_BUILD(ch_input_build)

            // ADd the software version
            ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)

            // Add the ref field to the meta
            input
                .map{ meta, file ->
                    def meta_with_ref = meta + [ref: ref]
                    return [meta_with_ref, file]
                }
                .set{ch_input_ref}
            
            // Run Bowtie for each entry in the "input" channel
            BOWTIE_ALIGN_DB(ch_input_ref, BOWTIE_BUILD.out.index, true, true)

            // Add the software version
            ch_versions = ch_versions.mix(BOWTIE_ALIGN_DB.out.versions)

            // Get the alignment results (summary)
            BOWTIE_ALIGN_DB.out.log
                .map { meta, file ->
                    // Leer el contenido del archivo
                    def content = file.text

                    // Extraer las líneas relevantes del archivo
                    def lines = content.readLines()
                    def total = lines[0].find(/: (\d+)/) { match -> match[1] } // Total reads
                    def aligned = lines[1].find(/: (\d+ \(\d+\.\d+%\))/) { match -> match[1] } // Reads with alignment
                    def failed = lines[2].find(/: (\d+ \(\d+\.\d+%\))/) { match -> match[1] } // Reads that failed

                    // Crear un Map con los campos requeridos
                    [meta.id, meta + [
                        filtering_db_total: total,
                        filtering_db_only_align: aligned,
                        filtering_db_failed: failed
                    ], file]
                }
                .set {alingment_summary}
            
            // Save the output into the output channels
            ch_bowtie_bam       = BOWTIE_ALIGN_DB.out.bam
            ch_bowtie_aligned   = BOWTIE_ALIGN_DB.out.aligned
            ch_bowtie_unaligned = BOWTIE_ALIGN_DB.out.unaligned

        } else if (ref == 'genome') {

            // Modify the channel to later combine it with the indexed genome.
            input
                .map{ meta, file ->
                    return [[id: meta.species], meta, file]
                }
                .set{ch_files_to_align}

            // Create a channel for the genomes [[id: Arabidopsis thaliana, etc], genome_file]
            input
                .map{ meta, _file ->
                    return [[id: meta.species], meta.genome]
                }
                .unique()
                .set{ch_genomes}

            // Index the reference genomes
            BOWTIE_BUILD(ch_genomes)

            // Add the software version
            ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)

            // Combine the genome and file channels 
            BOWTIE_BUILD.out.index
                .combine(ch_files_to_align, by:0)
                .multiMap{ it ->
                    // Create two channels: the query and the reference
                    query: [it[2] + [ref: "genome"], it[3]]
                    reference: [it[2] + [ref: "genome"], it [1]]
                }
                .set{ch_input_alignment}

            // Run Bowtie for each entry in the "input" channel
            BOWTIE_ALIGN_GENOME(ch_input_alignment.query, ch_input_alignment.reference, true, true)

            // Add the software version
            ch_versions = ch_versions.mix(BOWTIE_ALIGN_GENOME.out.versions)

            // Get the alignment results (summary)
            BOWTIE_ALIGN_GENOME.out.log
                .map { meta, file ->

                    // Read the file
                    def content = file.text

                    // Relevant rows of th file
                    def lines = content.readLines()
                    def total = lines[0].find(/: (\d+)/) { match -> match[1] } // Total reads
                    def aligned = lines[1].find(/: (\d+ \(\d+\.\d+%\))/) { match -> match[1] } // Reads with alignment
                    def failed = lines[2].find(/: (\d+ \(\d+\.\d+%\))/) { match -> match[1] } // Reads that failed

                    // Add the results to the meta
                    [meta.id, meta + [
                        filtering_genome_total: total,
                        filtering_genome_only_align: aligned,
                        filtering_genome_failed: failed
                    ], file]
                }
                .set {alingment_summary}

            // Save the output into the output channels
            ch_bowtie_bam       = BOWTIE_ALIGN_GENOME.out.bam
            ch_bowtie_aligned   = BOWTIE_ALIGN_GENOME.out.aligned
            ch_bowtie_unaligned = BOWTIE_ALIGN_GENOME.out.unaligned

        }
        
        // Add the summary results to the meta section of BAM channel
        ch_bowtie_bam
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_bam_out }

        // Add the summary results to the meta section of aligned channel
        ch_bowtie_aligned
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_aligned_out }
        
        // Add the summary results to the meta section of unaligned channel
        ch_bowtie_unaligned
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_unaligned_out }

    emit:
        bam         = ch_bam_out
        aligned     = ch_aligned_out
        unaligned   = ch_unaligned_out
        versions    = ch_versions      // channel: [ path(versions.yml) ]
}
