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

/*
========================================================================================
    Workflow ID_RESOLUTION
========================================================================================
*/

workflow ID_RESOLUTION {
    take:
        input_species_names     // value : species names
        mirbase                 // boolean: true/false
        mirbase_taxon           // value: species taxon
        srnaanno                // boolean: true/false
        pmiren                  // boolean: true/false

    main:

        // Download the databases
        DOWNLOADDB(mirbase, srnaanno, pmiren)

        // Resolve conflicts in identifiers between databases
        DB_ID_RESOLVER(
            input_species_names,
            DOWNLOADDB.out.mirbase_mature,
            DOWNLOADDB.out.mirbase_hairpin,
            DOWNLOADDB.out.mirbase_species,
            mirbase_taxon,
            DOWNLOADDB.out.srnaanno_all,
            DOWNLOADDB.out.pmiren_mature,
            DOWNLOADDB.out.pmiren_hairpin
        )

        // Prepare species_ids channel
        DB_ID_RESOLVER.out.species_ids
            .splitCsv(header: true)
            .set { ch_species_ids }
    emit:
        mirbase_mature =  DB_ID_RESOLVER.out.mirbase_mature
        mirbase_hairpin = DB_ID_RESOLVER.out.mirbase_hairpin
        srnaanno_mature = DB_ID_RESOLVER.out.srnaanno_mature
        srnaanno_hairpin = DB_ID_RESOLVER.out.srnaanno_hairpin
        pmiren_mature = DB_ID_RESOLVER.out.pmiren_mature
        pmiren_hairpin = DB_ID_RESOLVER.out.pmiren_hairpin
        species_ids = ch_species_ids
}

