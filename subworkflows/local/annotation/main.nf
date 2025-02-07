#!/usr/bin/env nextflow

/*
========================================================================================
    ANNOTATION Sub-Workflow
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { ID_RESOLUTION     } from "../idresolution"
include { MIRNA_ANNOTATION  } from "../../../modules/local/annotation"
include { GROUP_BY_FAMILY   } from "../../../modules/local/grouping"

/*
========================================================================================
    Workflow ANNOTATION
========================================================================================
*/


workflow ANNOTATION {
    take:
        dea_files
        mirbase_annot
        mirbase_taxon
        srnaanno_annot
        pmiren_annot
        mismatches
        ea_summary_file
        mww_pvalue_thrshld
        min_num_db

    main:

        // Get the input species names
        dea_files
            .map{it[0].species}
            .unique()
            .collect()
            .set{ch_sp_names}

        // Prepare identifiers for the species
        ID_RESOLUTION(
            ch_sp_names,
            mirbase_annot,
            mirbase_taxon,
            srnaanno_annot,
            pmiren_annot
        )
        
        dea_files
            .map{ item -> [item[0].species, item[0], item[1]]}
            .combine(ID_RESOLUTION.out.species_ids, by: 0)
            .map{_, meta, file, species_id ->
                [meta + [species_id: species_id], file]
            }
            .set{dea_files}

        // Identify which differentially expressed sRNA sequences are miRNAs.
        MIRNA_ANNOTATION(
            dea_files,
            ID_RESOLUTION.out.mirbase_mature,
            ID_RESOLUTION.out.pmiren_mature,
            ID_RESOLUTION.out.srnaanno_mature,
            mismatches,
            ea_summary_file,
            mww_pvalue_thrshld,
            min_num_db
        )
        
        // Combine both dea_files and miRNA_annot_filt channel
        // Prepare ch_miRNA_annot_filt channel
        MIRNA_ANNOTATION.out.annotfilt
            .map{meta, file -> return[meta.id, meta, file]}
            .set{ch_mirna_annot_filt}
        // Prepate dea_files channel
        dea_files
            .map{meta, file -> return[meta.id, meta, file]}
            .set{dea_files_ch}
        // Combine both channels
        dea_files_ch.combine(ch_mirna_annot_filt, by:0)
            .map{item -> return[item[1], item[2], item[4]]}
            .set{ch_group_miRNAs_input}

        ch_group_miRNAs_input.view()
        
        // Group miRNAs into families
        GROUP_BY_FAMILY(ch_group_miRNAs_input)

    emit:
        annot = MIRNA_ANNOTATION.out.annot
        annotfilt = MIRNA_ANNOTATION.out.annotfilt
        annot_sum = MIRNA_ANNOTATION.out.sum
        annot_sumlen = MIRNA_ANNOTATION.out.sumlen
        fam_annot = GROUP_BY_FAMILY.out.fam_annot
        fam_boxplot = GROUP_BY_FAMILY.out.fam_boxplot
        fam_sum = GROUP_BY_FAMILY.out.fam_sum
}
