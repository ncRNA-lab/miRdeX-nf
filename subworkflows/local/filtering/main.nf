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

include { BOWTIE_BUILD } from '../../../modules/nf-core/bowtie/build'
include { BOWTIE_ALIGN       } from '../../../modules/local/bowtie/align'
/*
========================================================================================
    Workflow FILTERING
========================================================================================
*/

/// ATENCIÓN. ESTE SUBWORKFLOW ESTÁ SIENDO EDITADO. TERMINAR.
// Specify DSL2
// nextflow.enable.dsl=2

workflow FILTERING {
    take:
        input                   // channel: [[species:val(species), species_id:val(species_id), project:val(project)], metadata:val(metadata)], genome:val(genome)], file]
        mismatches              // value: number of mismatches
        type                    // value: 'database' or 'genome'
        database                // path: reference fasta file (optional)

    main:

        // Branch the workflow based on the value of "type"
        if (type == 'database') {

            ch_input_build = Channel.value([[id:'Filtering database'], database])
            
            // Index the unique reference file (database)
            BOWTIE_BUILD(ch_input_build)

            // Add the type field to the meta
            input
                .map{ meta, file ->
                    def meta_with_type = meta + [type: type]
                    return [meta_with_type, file]
                }
                .set{ch_input_type}

            // Run Bowtie for each entry in the "input" channel
            BOWTIE_ALIGN(ch_input_type, BOWTIE_BUILD.out.index, true, true)

            // Get the alignment results (summary)
            BOWTIE_ALIGN.out.log
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

        } else if (type == 'genome') {

            // Modify the channel to later combine it with the indexed genome.
            input
                .map{ meta, file ->
                    return [[id: meta.species_id], meta, file]
                }
                .set{ch_files_to_align}

            // Create a channel for the genomes [[id: Arabidopsis thaliana, etc], genome_file]
            input
                .map{ meta, file ->
                    return [[id: meta.species_id], meta.genome]
                }
                .unique()
                .set{ch_genomes}

            // Index the reference genomes
            BOWTIE_BUILD(ch_genomes)

            // Combine the genome and file channels 
            BOWTIE_BUILD.out.index
                .combine(ch_files_to_align, by:0)
                .multiMap{ it ->
                    // Create two channels: the query and the reference
                    query: [it[2] + [type: "genome"], it[3]]
                    reference: [it[2] + [type: "genome"], it [1]]
                }
                .set{ch_input_alignment}

            // Run Bowtie for each entry in the "input" channel
            BOWTIE_ALIGN(ch_input_alignment.query, ch_input_alignment.reference, true, true)

            // Get the alignment results (summary)
            BOWTIE_ALIGN.out.log
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

        }
        
        // Add the summary results to the meta section of BAM channel
        BOWTIE_ALIGN.out.bam
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_bam_out }

        // Add the summary results to the meta section of aligned channel
        BOWTIE_ALIGN.out.aligned
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_aligned_out }
        
        // Add the summary results to the meta section of unaligned channel
        BOWTIE_ALIGN.out.unaligned
            .map{ meta, file -> [meta.id, meta, file] }
            .combine(alingment_summary, by:0)
            .map{item -> [item[3], item[2]]}
            .set{ ch_unaligned_out }

    emit:
        bam = ch_bam_out
        aligned = ch_aligned_out
        unaligned = ch_unaligned_out

}
