/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//

include { COUNTS } from '../../../modules/local/counts'
include { RPM    } from '../../../modules/local/rpm'
include { ADD_COUNTS_TO_ISOMIRS_DF as ADD_RAW_COUNTS_TO_ISOMIRS_DF  } from '../../../modules/local/add_counts_to_isomirs_df'
include { ADD_COUNTS_TO_ISOMIRS_DF as ADD_RPM_TO_ISOMIRS_DF         } from '../../../modules/local/add_counts_to_isomirs_df'
include { ISOMIRS_MIRNA_CLASSIFICATION                              } from '../../../modules/local/isomirs_mirna_classification'

//
// SUBWORKFLOW: Loaded from subworkflows/local/
//

include { summaryToTsv } from '../utils_mirplan_pipeline'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { PREPARE_MIRNA_DATABASES } from '../prepare_mirna_databases'
include { ISOMIRS_IDENTIFICATION  } from '../isomirs_identification'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { SEQKIT_FQ2FA } from '../../../modules/nf-core/seqkit/fq2fa/main'
include { SEQKIT_RMDUP } from '../../../modules/nf-core/seqkit/rmdup/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_MATURE    } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_MAKEBLASTDB as BLAST_MAKEBLASTDB_PRECURSOR } from '../../../modules/nf-core/blast/makeblastdb/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MIRNOTE {

    take:
    ch_input                // channel: [[id:val(id), species:val(species), genome:val(genome)], fasta/fastq]
    databases               // string: "mirbase,srnaanno,pmiren"
    substitutions           // integer: 1
    five_add                // integer: 0
    three_add               // integer: 3
    ends_modification       // integer: 4

    main:
    
    // Create an empty channel for summary info
    ch_pipeline_summary  = Channel.empty()
    ch_versions  = Channel.empty()
    
    /*
    ========================================================================================
        1. Prepare miRNA databases
    ========================================================================================
    */

    // Get the input species names
    ch_input
        .map{it[0].species}
        .unique()
        .collect()
        .set{ch_species_names}

    // Prepare the databases for the annotation
    PREPARE_MIRNA_DATABASES(ch_species_names, databases)

    // Use the species_id to combine the databases and input channels
    PREPARE_MIRNA_DATABASES.out.dbs
        .map{ meta -> [meta.species, meta]}
        .set{ ch_species_db }

    // Add info to summary channel
    ch_input
        .map{ meta, file -> [meta.species, meta, file]}
        .combine(ch_species_db)
        .map { item ->
            // Create required variables
            def meta = item[1]
            def updatedMeta = meta.clone()
            def species_db_info = item[4]

            // Remove undesired fields for summary channel
            updatedMeta.remove('genome')
            updatedMeta.remove('single_end') 

            // Add fields based on species_db_info
            def additionalFields = [:]
            if (species_db_info) {
                additionalFields['species_db'] = 'PASS'
                additionalFields['database'] = species_db_info.database
            } else {
                additionalFields['species_db'] = 'REJECT'
                additionalFields['database'] = 'NA'
            }

            // Combine fields
            [updatedMeta.id, updatedMeta + additionalFields]
        }
        .set{ ch_pipeline_summary }

    // Combine databases and input channels
    ch_input
        .map{ meta, file -> [meta.species, meta, file]}
        .combine(ch_species_db, by:0)
        .map{_species, meta1, file, meta2 -> [meta1 + meta2, file]}
        .set{ ch_input }


    /*
    ========================================================================================
        2. Quantify sequences
    ========================================================================================
    */

    // Check whether a raw counts or RPM threshold is to be used.
    if (params.counts > 0 || params.rpm > 0){
        
        // Calculate the raw counts
        COUNTS(ch_input)

        // Save the raw counts inot a channel
        COUNTS.out.raw.set{ ch_raw_counts }

        // Check whether a RPM threshold is to be used.
        if (params.rpm  > 0){

            // Calculate RPM
            RPM(ch_raw_counts)

            // Save the raw counts inot a channel
            RPM.out.rpm.set{ ch_rpm }
        }
    }

    /*
    ========================================================================================
        3. Prepare input files
    ========================================================================================
    */

    // Identify FASTA and FASTQ files in the input channel
    sub_fastq_ch = ch_input.filter{ _meta, file -> file.name ==~ /.*\.(fastq|fq)(\.gz)?$/ }
    sub_fasta_ch = ch_input.filter{ _meta, file -> file.name ==~ /.*\.(fasta|fa)(\.gz)?$/ }

    // Convert FASTQ files to FASTA
    SEQKIT_FQ2FA(sub_fastq_ch)

    // Save the software version
    ch_versions = ch_versions.mix(SEQKIT_FQ2FA.out.versions)

    // Merge the contents of the two previous channels into the same channel
    sub_fasta_ch
        .mix(SEQKIT_FQ2FA.out.fasta)
        .map{ meta, file -> [ meta + [id: "${meta.id}.rmdup", run: meta.id], file]}
        .set{ch_input}

    // Remove duplicates in the input files
    SEQKIT_RMDUP(ch_input)

    // Save the software version
    ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions)

    // Save the results in the input channel
    SEQKIT_RMDUP.out.fastx
        .map{ meta, file -> [meta + [id: meta.run], file]}
        .set{ch_input}

    /*
    ========================================================================================
        4. Identify isomiRs
    ========================================================================================
    */

    // Preapare the input channel for ISOMIRS_IDENTIFICATION worflow
    ch_input
        .map{ meta, file -> [meta, file, meta.mature, meta.precursor, meta.genome]}
        .set{ ch_to_identify_isomirs }
    
    // Identify isomiRs
    ISOMIRS_IDENTIFICATION(ch_to_identify_isomirs, ch_versions)

    // Save the software version
    ch_versions = ch_versions.mix(ISOMIRS_IDENTIFICATION.out.versions)

    // Get the numnber of sequences identify as potential isomiRs (summary)
    ISOMIRS_IDENTIFICATION.out.iso
        .map { meta, file ->
            def iso = file.readLines()
                        .collect { it.split('\t', -1)[13] }
                        .findAll { it }
                        .toSet()
            [meta.id, iso.size()]
        }
        .join(ch_pipeline_summary, remainder:true)
        .map{ item ->
            // Required items
            def id = item[0]
            def num_iso = item[1]
            def meta = item[2]

            // Fields to add
            def additionalFields = num_iso ? [ num_pot_isomirs: num_iso ] : [ num_pot_isomirs: 'NA']
            
            // Return
            [id, meta + additionalFields]
        }
        .set{ ch_pipeline_summary }

    /*
    ========================================================================================
        5. Add raw counts and RPM to the Blastn isomiRs dataframe
    ========================================================================================
    */

    // Input channel for ISOMIRS_MIRNA_CLASSIFICATION module
    ch_isomirs = ISOMIRS_IDENTIFICATION.out.iso

    // Add counts to isomiRs dataframe
    if (params.counts  > 0 || params.rpm  > 0) {

        // Prepare the channel for merging
        ch_isomirs
            .map{ meta, file -> [meta.id, meta, file]}
            .set{ ch_isomirs }

        // Prepare the channel for merging
        ch_raw_counts
            .map{ meta, file -> [meta.id, meta, file]}
            .combine(ch_isomirs, by:0)
            .map{ _id, _meta_c, file_c, meta_i, file_i ->
                [meta_i + [id: "${meta_i.id}.rawc", prev_id:meta_i.id], file_i, file_c]
            }
            .set { ch_isomirs_add_raw_counts }
        
        // Add raw counts to isomiRs dataframe
        ADD_RAW_COUNTS_TO_ISOMIRS_DF(ch_isomirs_add_raw_counts, params.counts)


        // Add RPM to isomiRs dataframe
        if (params.rpm  > 0){

            // Set the original id and prepare the channel for merging
            ADD_RAW_COUNTS_TO_ISOMIRS_DF.out.isocounts
                .map { meta, file ->
                    def updatedMeta = meta + [id: meta.prev_id]
                    updatedMeta.remove('prev_id') 
                    return [updatedMeta.id, updatedMeta, file]
                }
                .set { ch_isomirs }

            // Prepare the channel for merging
            ch_rpm
                .map{ meta, file -> [meta.id, meta, file]}
                .combine(ch_isomirs, by:0)
                .map{ _id, _meta_c, file_c, meta_i, file_i ->
                    [meta_i + [id: "${meta_i.id}.rpmc", prev_id:meta_i.id], file_i, file_c]
                }
                .set { ch_isomirs_add_rpm }
    
            // Add RPM to isomiRs dataframe
            ADD_RPM_TO_ISOMIRS_DF(ch_isomirs_add_rpm, params.rpm)
            
            // Set the original id and prepare the channel for merging
            ADD_RPM_TO_ISOMIRS_DF.out.isocounts
                .map { meta, file ->
                    def updatedMeta = meta + [id: meta.prev_id]
                    updatedMeta.remove('prev_id') 
                    return [updatedMeta, file]
                }
                .set { ch_isomirs }    

        } else {

            // Set the original id and prepare the channel for isomir classifcation
            ADD_RAW_COUNTS_TO_ISOMIRS_DF.out.isocounts
                .map { meta, file ->
                    def updatedMeta = meta + [id: meta.prev_id]
                    updatedMeta.remove('prev_id') 
                    return [updatedMeta, file]
                }
                .set { ch_isomirs }
        }
    }

    /*
    ========================================================================================
        6. Classify isomiRs
    ========================================================================================
    */
    
    // Classify isomiRs
    ISOMIRS_MIRNA_CLASSIFICATION(
        ch_isomirs,
        substitutions,
        five_add,
        three_add,
        ends_modification
    )

    // Save the software version
    ch_versions = ch_versions.mix(ISOMIRS_MIRNA_CLASSIFICATION.out.versions)

    ISOMIRS_MIRNA_CLASSIFICATION.out.sum
        .splitCsv( header: true, sep: '\t' )
        .map{ meta, info -> [meta.id, info]}
        .join(ch_pipeline_summary, remainder:true)
        .map { item ->
            def iso_info = item[1]
            def pip_summary = item[2]

            // Fields to add
            def additionalFields = iso_info ? [
                num_isomirs: iso_info.num_isomirs,
                ref_miRNA: iso_info.ref_miRNA,
                iso_5p: iso_info.iso_5p,
                iso_3p: iso_info.iso_3p,
                iso_add3p: iso_info.iso_add3p,
                iso_add5p: iso_info.iso_add5p,
                iso_snv_seed: iso_info.iso_snv_seed,
                iso_snv_central_offset: iso_info.iso_snv_central_offset,
                iso_snv_central: iso_info.iso_snv_central,
                iso_snv_central_supp: iso_info.iso_snv_central_supp,
                iso_snv: iso_info.iso_snv,
                mixed: iso_info.mixed,
                mixed_shift: iso_info.mixed_shift,
                undefined: iso_info.undefined
            ] : [
                num_isomirs: 'NA',
                ref_miRNA: 'NA',
                iso_5p: 'NA',
                iso_3p: 'NA',
                iso_add3p: 'NA',
                iso_add5p: 'NA',
                iso_snv_seed: 'NA',
                iso_snv_central_offset: 'NA',
                iso_snv_central: 'NA',
                iso_snv_central_supp: 'NA',
                iso_snv: 'NA',
                mixed: 'NA',
                mixed_shift: 'NA',
                undefined: 'NA'
            ]
            pip_summary + additionalFields
        }
        .set{ ch_pipeline_summary }
        
    emit:
    annotation     = ISOMIRS_MIRNA_CLASSIFICATION.out.gff3   // channel: [meta, file.gff3]
    summary        = ch_pipeline_summary                     // channel: [id:, num_isomirs:, etc]
    versions       = ch_versions                             // channel: [ path(versions.yml) ]

}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/