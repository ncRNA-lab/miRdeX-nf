/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { MERGE_AND_FILTER_MATURE_PRECURSOR_BLAST } from '../../../modules/local/merge_and_filter_mature_precursor_blast'
include { ISOMIRS_PRECURSOR_CLASSIFICATION        } from '../../../modules/local/isomirs_precursor_classification'
include { TSV_TO_FASTA                            } from '../../../modules/local/tsv_to_fasta'
include { BOWTIE_ALIGN                            } from '../../../modules/local/bowtie/align'
include { FILTER_TSV_BY_FASTA                     } from '../../../modules/local/filter_tsv_by_fasta'
include { CONCAT_TSV                              } from '../../../modules/local/concat_tsv'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//

include { COMPLETE_BLASTN as COMPLETE_BLASTN_MATURE    } from '../../../subworkflows/local/complete_blastn'
include { COMPLETE_BLASTN as COMPLETE_BLASTN_PRECURSOR } from '../../../subworkflows/local/complete_blastn'
include { FILTER_NONTEMPLATED_ISOMIRS                  } from '../../../subworkflows/local/filter_nontemplated_isomirs'

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
    WORKFLOW ISOMIRS_IDENTIFICATION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow ISOMIRS_IDENTIFICATION {
    take:
        ch_input           // channel: [[id:val(id)], path(query_fasta), path(mature_fasta), path(precursor_fasta), path(genome)]
        ch_versions        // channel: [ path(versions.yml) ]

    main:

        /*
        ========================================================================
            Align sequences against the databases using blastn
        ========================================================================
        */

        //
        // 1. Alignment with the mature miRNA database.
        //

        // Prepare the input channel for the blastn subworkflow (mature)
        ch_input
            .map { meta, q_file, m_file, _p_file, genome ->
                def id = genome.toString().tokenize('/')[-1].replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
                return [meta + [genome_id: id], q_file, m_file]
            }
            .set { ch_mature_to_blast }

        // Run Blast using mature miRNAs as reference
        COMPLETE_BLASTN_MATURE(ch_mature_to_blast, ch_versions)
    
        // Save the software version
        ch_versions = ch_versions.mix(COMPLETE_BLASTN_MATURE.out.versions)

        // Prepare the channel for merging
        COMPLETE_BLASTN_MATURE.out.blast
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_mature_to_merge }

        //
        // 2. Alignment with the miRNA precursor database.
        //

        // Prepare the input channel for the blastn subworkflow (precursor)
        ch_input
            .map { meta, q_file, _m_file, p_file, genome -> 
                def id = genome.toString().tokenize('/')[-1].replaceAll(/\.(fa|fasta|fna)(\.gz)?$/, '')
                return [meta + [genome_id: id], q_file, p_file]
            }
            .set { ch_precursor_to_blast }

        // Run Blast using miRNA precursors as reference
        COMPLETE_BLASTN_PRECURSOR(ch_precursor_to_blast, ch_versions)

        // Prepare the channel for the MERGE_BLAST process
        COMPLETE_BLASTN_PRECURSOR.out.blast
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_mature_to_merge, by: 0)
            .map{ _id, _meta_pre, blast_pre, meta_mat, blast_mat -> [meta_mat, blast_mat, blast_pre]}
            .set{ ch_to_merge_blast_results }

        // Get a final BLASTN results table
        MERGE_AND_FILTER_MATURE_PRECURSOR_BLAST(ch_to_merge_blast_results)

        /*
        ========================================================================
            Perform a global classification of the sequences
                - Canonical miRNAs
                - Non-canonical templated miRNAs
                - Non-canonical non-templated miRNAs
        ========================================================================
        */

        //
        // 1. Classificate the sequences
        //

        ISOMIRS_PRECURSOR_CLASSIFICATION(MERGE_AND_FILTER_MATURE_PRECURSOR_BLAST.out.mpblast)
    
        //
        // 2. Remove non-templated sequences that align with the genome.
        //

        // Prepare the input channel for FILTER_NONTEMPLATED_ISOMIRS process
        ch_input
            .map{meta, _file, _mat, _pre, genome -> [meta.id, genome]}
            .set{ ch_genome }

        ISOMIRS_PRECURSOR_CLASSIFICATION.out.nontemplated
            .map{ meta, file  -> [meta.id, meta, file]}
            .combine(ch_genome, by:0)
            .map{ _id, meta, file, genome -> [meta, file, genome]}
            .set{ ch_nontemplated_and_genome }

        // Remove non-templated isomiRs that align to other regions of the genome.
        FILTER_NONTEMPLATED_ISOMIRS(ch_nontemplated_and_genome, ch_versions)

        // Save the software version
        ch_versions = ch_versions.mix(FILTER_NONTEMPLATED_ISOMIRS.out.versions)

        //
        // 3. Concat the canonical, templated, and non-templated miRNA files
        //

        // Prepare the non-templated isomirs channel
        FILTER_NONTEMPLATED_ISOMIRS.out.iso
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_non_templated }

        // Prepare the templated isomirs channel
        ISOMIRS_PRECURSOR_CLASSIFICATION.out.templated
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_templated }

        // Prepare the canonical miRNAs channel and combine it with the two previous ones.
        ISOMIRS_PRECURSOR_CLASSIFICATION.out.canon
            .map{ meta, file -> [meta.id, meta, file]}
            .join(ch_templated, by:0, remainder:true)
            .map{ item ->
                if (item[3] == null){ 
                    return [item[0], item[1], item[2], null]
                } else if (item[1] == null){
                    return [item[0], item[2], null, item[3]]
                } else {
                    return [item[0], item[1], item[2], item[4]]
                }
            }
            .join(ch_non_templated, by:0, remainder:true)
            .map{ item ->
                if (item[1] == null){
                    return [item[2] + [id: "${item[2].id}.iso", prev_id: item[2].id], [item[3]]]
                } else if (item[4] == null){ 
                    return [item[1] + [id: "${item[1].id}.iso", prev_id: item[1].id], [item[2], item[3]]]
                } else {
                    return [item[1] + [id: "${item[1].id}.iso", prev_id: item[1].id], [item[2], item[3], item[5]]]
                }
            }
            .set { ch_files_to_concat }

        // Concatenate canonical, templated and non-templated sequences files
        CONCAT_TSV(ch_files_to_concat, false)

        // Set the original id
        CONCAT_TSV.out.concat
            .map { meta, file ->
                def updatedMeta = meta + [id: meta.prev_id]
                updatedMeta.remove('prev_id') 
                return [updatedMeta, file]
            }
            .set { ch_output }

    emit:
        iso      = ch_output        // channel: [[id:val(id)], path(tsv_file)]
        versions = ch_versions      // channel: [ path(versions.yml) ]
}
