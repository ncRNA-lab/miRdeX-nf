/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { TSV_TO_FASTA                                    } from '../../../modules/local/tsv_to_fasta'
include { BOWTIE_ALIGN as BOWTIE_ALIGN_NONTEMPLATED       } from '../../../modules/local/bowtie/align'
include { FILTER_TSV_BY_FASTA                             } from '../../../modules/local/filter_tsv_by_fasta'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SEQKIT_RMDUP } from '../../../modules/nf-core/seqkit/rmdup/main'
include { BOWTIE_BUILD } from '../../../modules/nf-core/bowtie/build'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW FILTER_NONTEMPLATED_ISOMIRS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FILTER_NONTEMPLATED_ISOMIRS {
    take:
        ch_input           // channel: [[id:val(id)], path(blast_tsv), path(genome)]

    main:

        // Create empty channel for versions
        ch_versions  = Channel.empty()
        
        /*
        ========================================================================
            1. Create the query fasta file
        ========================================================================
        */

        // Prepare the input channel
        ch_input
            .map{ meta, file, genome ->
                def id = genome.toString().tokenize('/')[-1].replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
                [meta + [genome_id_filt: id], file]}
            .set{ ch_input_file }

        // Create a FASTA file with the non-templated sequences
        TSV_TO_FASTA(ch_input_file, 1, 14, false)

        // Remove duplicates
        SEQKIT_RMDUP(TSV_TO_FASTA.out.fasta)

        // Save the software version
        ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

        /*
        ========================================================================
            2. Index the reference genome
        ========================================================================
        */

        // Get the genomes to be indexed
        ch_input
            .map { _meta, _file, genome ->
                def id = genome.toString().tokenize('/')[-1].replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
                return [[id: id], genome]
            }
            .unique()
            .set{ ch_genomes }

        // Index the genomes
        BOWTIE_BUILD(ch_genomes)

        // Save the software version
        ch_versions = ch_versions.mix(BOWTIE_BUILD.out.versions)

        /*
        ========================================================================
            3. Align non-templated sequences with the reference genome
        ========================================================================
        */

        // Add the indexed genome to the main channel
        BOWTIE_BUILD.out.index
           .map{ meta, file -> [meta.id, file]}
           .set{ ch_genomes_indexed }

        // Prepare the input channel for alignment with Bowtie
        SEQKIT_RMDUP.out.fastx
            .map { meta, file -> [meta.genome_id_filt, meta, file] }
            .combine(ch_genomes_indexed, by:0)
            .map{ _id, meta, file, idx -> [meta, file, idx]}
            .set{ ch_nontemplated_bowtie_input }
        
        // Create the query channel
        ch_nontemplated_bowtie_input
            .map { item ->
                return [item[0], item[1]]
            }
            .set{ ch_nontemplated_bowtie_input_query }

        // Create the reference channel
        ch_nontemplated_bowtie_input
            .map { item ->
                return [item[0], item[2]]
            }
            .set{ ch_nontemplated_bowtie_input_ref }

        // Run bowtie
        BOWTIE_ALIGN_NONTEMPLATED(
            ch_nontemplated_bowtie_input_query,
            ch_nontemplated_bowtie_input_ref,
            true,
            true
        )
        
        // Save the software version
        ch_versions = ch_versions.mix(BOWTIE_ALIGN_NONTEMPLATED.out.versions)

        /*
        ========================================================================
            4. Remove non-templated isomiRs that align to other regions of
               the genome from the original file
        ========================================================================
        */

        // Prepare the input channel to filter non-templated sequences
        BOWTIE_ALIGN_NONTEMPLATED.out.unaligned
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_nontemplated_unalig_fasta }

        ch_input
            .map{ meta, file, _genome -> [meta.id, meta, file]}
            .combine(ch_nontemplated_unalig_fasta, by:0)
            .map{ _id, meta1, tsv, _meta2, fasta -> [meta1, tsv, fasta]}
            .set{ ch_to_filter_non_templated }

        // Filter non-templated sequences
        FILTER_TSV_BY_FASTA(ch_to_filter_non_templated, 14, true)

    emit:
        iso      = FILTER_TSV_BY_FASTA.out.filt_tsv      // channel: [[id:val(id)], path(tsv_file)]
        versions = ch_versions                           // channel: [ path(versions.yml) ]
}
