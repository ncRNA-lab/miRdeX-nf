/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { ADD_SEQUENCES_BLAST } from '../../../modules/local/add_sequences_blast'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { BLAST_MAKEBLASTDB } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN      } from '../../../modules/nf-core/blast/blastn/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    WORKFLOW COMPLETE_BLASTN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPLETE_BLASTN {
    take:
        ch_input           // channel: [[id:val(id)], path(query_fasta), path(reference_fasta)]
        ch_versions        // channel: [ path(versions.yml) ]

    main:
            
        /*
        ========================================================================
            Index reference database
        ========================================================================
        */

        // Create a channel for the reference file (not indexed)
        ch_input
            .map { meta, q_file, r_file ->
                def id = r_file.simpleName
                return tuple(id, meta + [ref: r_file, ref_id: id], q_file, r_file)
            }
            .tap{ ch_input_to_merge }
            .map{ item -> [[id:item[0]], item[3]]}
            .unique()
            .set { ch_db_to_index }
        
        // Index the reference file
        BLAST_MAKEBLASTDB(ch_db_to_index)

        // Save the software version
        ch_versions = ch_versions.mix(BLAST_MAKEBLASTDB.out.versions)

        // Adapt the channel 
        BLAST_MAKEBLASTDB.out.db
            .map{meta, file -> [meta.id, meta, file] }
            .set{ch_db_indexed}

        // Add the indexed database to the original channel
        ch_input_to_merge
            .combine(ch_db_indexed, by:0)
            .map{ _id, meta1, file, _ref, _meta2, ref_idx -> [meta1, file, ref_idx] }
            .set{ ch_to_blast }

        // Split the channel in two channels (query)
        ch_to_blast
            .map{ meta, file, _ref_idx -> [meta.id, meta, file]}
            .tap{ ch_to_blast_query_merge }
            .map{ _id, meta, file -> [meta, file]}  
            .set{ ch_to_blast_query }

        // Split the channel in two channels (reference)
        ch_to_blast
            .map{ meta, _file, ref_idx -> [meta, ref_idx]}
            .set{ ch_to_blast_ref }
        
        /*
        ========================================================================
            Run BLASTN
        ========================================================================
        */
        
        // Run blastn
        BLAST_BLASTN(ch_to_blast_query, ch_to_blast_ref)

        // Save the software version
        ch_versions = ch_versions.mix(BLAST_BLASTN.out.versions)

        // Prepare the channel for adding sequences to the BLAST results
        BLAST_BLASTN.out.txt
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_to_blast_query_merge, by:0)
            .map{_id, meta1, file1, _meta2, file2 ->
                [meta1 + [id: "${meta1.id}_${meta1.ref_id}", prev_id: meta1.id], file1, file2, meta1.ref]
            }
            .set{ ch_to_add_sequences }

        /*
        ========================================================================
            Add query and reference sequences to BLASTN results
        ========================================================================
        */

        // Add the query and subject sequences to the blast results table
        ADD_SEQUENCES_BLAST(ch_to_add_sequences)

        // Save the software version
        ch_versions = ch_versions.mix(ADD_SEQUENCES_BLAST.out.versions)

        // Prepare the output channel
        ADD_SEQUENCES_BLAST.out.bseqs
            .map{meta, file -> 
                def updatedMeta = meta + [id: meta.prev_id]
                updatedMeta.remove('prev_id')  // Remove prev_id field
                updatedMeta.remove('ref_id')   // Remove ref_id field
                updatedMeta.remove('ref')      // Remove ref field
                [updatedMeta, file]
            }
            .set{ ch_output }

    emit:
        blast    = ch_output   // channel: [[id:val(id)], path(txt_file)]
        versions = ch_versions // channel: [ path(versions.yml) ]
}
