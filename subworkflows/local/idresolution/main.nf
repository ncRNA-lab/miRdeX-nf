#!/usr/bin/env nextflow

/*
========================================================================================
    ID_RESOLUTION Sub-Workflow
========================================================================================
*/

// Specify DSL2
nextflow.enable.dsl=2

/*
========================================================================================
    Include Modules
========================================================================================
*/

include { DOWNLOADDB } from '../../../modules/local/downloaddb'
include { DB_ID_RESOLVER } from "../../../modules/local/idresolver"
include { CHECK_SPECIES_IDS } from "../../../modules/local/idchecking"

/*
========================================================================================
    Workflow ID_RESOLUTION
========================================================================================
*/

workflow ID_RESOLUTION {
    take:
        input_species_names // value: species names

    main:
        // Check if the databases have been prepared in a previous run...
        if (file(params.databases_local.annotation.mirbase.all.mature).exists()) {
            
            // Use the previously created databases
            mirbase_mature = Channel.fromPath(params.databases_local.annotation.mirbase.all.mature).collect()
            mirbase_hairpin = Channel.fromPath(params.databases_local.annotation.mirbase.all.hairpin).collect()
            srnaanno_mature = Channel.fromPath(params.databases_local.annotation.srnaanno.all.mature).collect()
            srnaanno_hairpin = Channel.fromPath(params.databases_local.annotation.srnaanno.all.hairpin).collect()
            pmiren_mature = Channel.fromPath(params.databases_local.annotation.pmiren.all.mature).collect()
            pmiren_hairpin = Channel.fromPath(params.databases_local.annotation.pmiren.all.hairpin).collect()


            // Get the previously created species ids
            ch_old_species_ids = Channel.fromPath(params.databases_local.annotation.species_ids).collect()

            // Check the ids of the input species
            CHECK_SPECIES_IDS(input_species_names, ch_old_species_ids)

            // Use an output channel with the final ids file
            CHECK_SPECIES_IDS.out.species_ids
                .splitCsv(header: true)
                .map { row -> [row.Species_name, row.Final_id] }
                .set { species_ids }

        } else {

            // Download the databases
            DOWNLOADDB()

            // Mature miRNAs channels
            ch_mirbase_all_mature = DOWNLOADDB.out.mirbase_mature.collect()
            ch_pmiren_sp_mature= DOWNLOADDB.out.pmiren_mature.collect()

            // Hairpin channels
            ch_mirbase_all_hairpin = DOWNLOADDB.out.mirbase_hairpin.collect()
            ch_pmiren_sp_hairpin = DOWNLOADDB.out.pmiren_hairpin.collect()

            // Mature miRNAs + hairpin channels
            ch_srnaanno_sp_mature_hairpin = DOWNLOADDB.out.srnaanno_all.collect()

            // Other channels (CAMBIAR. PONER LA RUTA EN EL FICHERO CONFIG PRINCIPAL)
            ch_viridiplantae_mirbase_ids = Channel.fromPath(params.databases.annotation.mirbase.all.plant_id_file).collect()

            // Execute PREPARE_DATABASES process
            DB_ID_RESOLVER(
                input_species_names,
                ch_mirbase_all_mature,
                ch_pmiren_sp_mature,
                ch_srnaanno_sp_mature_hairpin,
                ch_mirbase_all_hairpin,
                ch_pmiren_sp_hairpin,
                ch_viridiplantae_mirbase_ids
            )

            // Prepare species_ids channel
            DB_ID_RESOLVER.out.species_ids
                .splitCsv(header: true)
                .map { row -> [row.Species_name, row.Final_id] }
                .set { ch_species_ids }
            
            // Prepare output channels
            mirbase_mature = DB_ID_RESOLVER.out.mirbase_mature
            mirbase_hairpin = DB_ID_RESOLVER.out.mirbase_hairpin
            srnaanno_mature = DB_ID_RESOLVER.out.srnaanno_mature
            srnaanno_hairpin = DB_ID_RESOLVER.out.srnaanno_hairpin
            pmiren_mature = DB_ID_RESOLVER.out.pmiren_mature
            pmiren_hairpin = DB_ID_RESOLVER.out.pmiren_hairpin
            species_ids = ch_species_ids

        }

    emit:
        mirbase_mature
        mirbase_hairpin
        srnaanno_mature
        srnaanno_hairpin
        pmiren_mature
        pmiren_hairpin
        species_ids
}

