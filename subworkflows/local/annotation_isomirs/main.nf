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

include { ID_RESOLUTION     } from "../idresolution_isomirs"
include { MIRNA_ANNOTATION  } from "../../../modules/local/annotation_isomirs"
include { GROUP_MIRNAS_BY_FAMILY   } from "../../../modules/local/group_mirnas_by_family"

/*
========================================================================================
    Workflow ANNOTATION
========================================================================================
*/


workflow ANNOTATION {
    take:
        dea_files
        databases
        mirbase_taxon

    main:

        // Get the input species names
        dea_files
            .map{it[0].species}
            .unique()
            .collect()
            .set{ch_sp_names}

        def db_priority = databases.split(',').collect { it.trim().toLowerCase() }

        // Create a boolean variable for each db
        def use_mirbase  = 'mirbase'  in db_priority
        def use_srnaanno = 'srnaanno' in db_priority
        def use_pmiren   = 'pmiren'   in db_priority

         // Prepare identifiers for the species
        ID_RESOLUTION(
            ch_sp_names,
            use_mirbase,
            mirbase_taxon,
            use_srnaanno,
            use_pmiren
        )

        // println(db_priority)

        // dea_ea_files.view()
        ID_RESOLUTION.out.species_ids
            .map{ item -> [item.species, item]}
            .set{ch_species_ids}

        // Procesamos los archivos para agregar la ruta de la base de datos en el campo `mature_db`
        dea_files
            .map{item -> [item[0].species, item[0], item[1]]}
            .combine(ch_species_ids, by:0)
            .map{ species, meta, dea_file, db_ids ->
                // Iterate through db_priority list items
                db_priority.find { db ->
                    // If the identifier is != NUll, we asign the database
                    if (db_ids[db] != "NULL") {
                        if (db == "mirbase") {
                            meta.annot_mature_db = ID_RESOLUTION.out.mirbase_mature.val
                            meta.annot_hairpin_db = ID_RESOLUTION.out.mirbase_hairpin.val
                            meta.species_id = db_ids.final_id
                            return true
                        }
                        if (db == "srnaanno") {
                            meta.annot_mature_db = ID_RESOLUTION.out.srnaanno_mature.val
                            meta.annot_hairpin_db = ID_RESOLUTION.out.srnaanno_hairpin.val
                            meta.species_id = db_ids.final_id
                            return true
                        }
                        if (db == "pmiren") {
                            meta.annot_mature_db = ID_RESOLUTION.out.pmiren_mature.val
                            meta.annot_hairpin_db = ID_RESOLUTION.out.pmiren_hairpin.val
                            meta.species_id = db_ids.final_id
                            return true
                        }
                        meta.annot_mature_db = null
                        meta.annot_hairpin_db = null
                        meta.species_id = null
                    }
                }
                return [meta, dea_file]
            }
            .set { dea_files }
            
        // Identify which differentially expressed sRNA sequences are miRNAs
        MIRNA_ANNOTATION(dea_files)
        MIRNA_ANNOTATION.out.isomirs.view()

    //     // Identify which differentially expressed sRNA sequences are miRNAs.
    //     MIRNA_ANNOTATION(
    //         dea_ea_files,
    //         ID_RESOLUTION.out.mirbase_mature,
    //         ID_RESOLUTION.out.pmiren_mature,
    //         ID_RESOLUTION.out.srnaanno_mature,
    //         mismatches,
    //         mww_pvalue_thrshld,
    //         min_num_db
    //     )
        
    //     // Prepare ch_miRNA_annot_filt channel
    //     MIRNA_ANNOTATION.out.annotfilt
    //         .map{meta, file -> return[meta.id, meta, file]}
    //         .set{ch_mirna_annot_filt}

    //     // Prepate dea_files channel
    //     dea_ea_files
    //         .map{meta, file, ea_file -> return[meta.id, meta, file, ea_file]}
    //         .set{ dea_ea_files }

    //     // Combine both channels
    //     dea_ea_files
    //         .combine(ch_mirna_annot_filt, by:0)
    //         .map{item -> return[item[1], item[2], item[5]]}
    //         .set{ch_group_miRNAs_input}

    //     // Group miRNAs into families
    //     GROUP_MIRNAS_BY_FAMILY(ch_group_miRNAs_input)

    // emit:
    //     annot = MIRNA_ANNOTATION.out.annot
    //     annotfilt = MIRNA_ANNOTATION.out.annotfilt
    //     annot_sum = MIRNA_ANNOTATION.out.sum
    //     annot_sumlen = MIRNA_ANNOTATION.out.sumlen
    //     fam_annot = GROUP_MIRNAS_BY_FAMILY.out.fam_annot
    //     fam_boxplot = GROUP_MIRNAS_BY_FAMILY.out.fam_boxplot
    //     fam_sum = GROUP_MIRNAS_BY_FAMILY.out.fam_sum
}
