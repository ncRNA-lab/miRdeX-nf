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

include { MIRNA_ANNOTATION; MIRNA_ANNOTATION_SUM    } from "../../../modules/local/annotation"
include { GROUP_BY_FAMILY; GROUP_BY_FAMILY_SUM      } from "../../../modules/local/grouping"

/*
========================================================================================
    Workflow ANNOTATION
========================================================================================
*/


workflow ANNOTATION {
    take:
        dea_files
        mirbase
        pmiren
        srnaanno
        mismatches
        ea_summary_file
        mww_pvalue_thrshld

    main:

        // Identify which differentially expressed sRNA sequences are miRNAs.
        MIRNA_ANNOTATION(
            dea_files,
            mirbase,
            pmiren,
            srnaanno,
            mismatches,
            ea_summary_file,
            mww_pvalue_thrshld
        )
        
        // Combine both dea_files and miRNA_annot_filt channel
        // Prepate miRNA_annot_filt channel
        MIRNA_ANNOTATION.out.annotfiltmat
            .map{meta, file -> return[meta.id, meta, file]}
            .set{mirna_annot_filt_ch}
        // Prepate dea_files channel
        dea_files
            .map{meta, file -> return[meta.id, meta, file]}
            .set{dea_files_ch}

        // Combine both channels
        dea_files_ch.combine(mirna_annot_filt_ch, by:0).
            map{item -> return[item[1], item[2], item[4]]}
            .set{group_miRNAs_input_ch}

    //     // Generate an annotation summary file
    //     MIRNA_ANNOTATION_SUM(
    //         MIRNA_ANNOTATION.out.sum_mat.collect(),
    //         MIRNA_ANNOTATION.out.sum_mat_len.collect(),
    //         MIRNA_ANNOTATION.out.sum_hair.collect(),
    //         MIRNA_ANNOTATION.out.sum_hair_len.collect()
    //     )

        // Group miRNAs into families
        GROUP_BY_FAMILY(group_miRNAs_input_ch)

    //     // Group miRNAs into families (summary)
    //     GROUP_BY_FAMILY_SUM(GROUP_BY_FAMILY.out.fam_sum.collect())

    emit:
        annotmat = MIRNA_ANNOTATION.out.annotmat
        annotfiltmat = MIRNA_ANNOTATION.out.annotfiltmat
        annothair = MIRNA_ANNOTATION.out.annothair
        annotfilthair = MIRNA_ANNOTATION.out.annotfilthair
        fam_annot = GROUP_BY_FAMILY.out.fam_annot
        fam_boxplot = GROUP_BY_FAMILY.out.fam_boxplot
        fam_sum = GROUP_BY_FAMILY.out.fam_sum
}
