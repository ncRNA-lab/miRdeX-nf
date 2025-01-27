#!/usr/bin/env nextflow

/*
========================================================================================
    DATABASES module
========================================================================================
*/
// Specify DSL2
nextflow.enable.dsl=2

process DB_ID_RESOLVER {

    debug true

    input:
        val input_species_names
        path mirbase
        path pmiren
        path srnaanno
        path mirbase_hairpin
        path pmiren_hairpin
        path v_ids_mirbase

    output:
        path "01-Mod_databases/miRBase/mirbase_mature.fa", emit: mirbase_mature
        path "01-Mod_databases/miRBase/mirbase_hairpin.fa", emit: mirbase_hairpin
        path "01-Mod_databases/sRNAanno/srnaanno_mature.fa", emit: srnaanno_mature
        path "01-Mod_databases/sRNAanno/srnaanno_hairpin.fa", emit: srnaanno_hairpin
        path "01-Mod_databases/PmiREN/pmiren_mature.fa", emit: pmiren_mature
        path "01-Mod_databases/PmiREN/pmiren_hairpin.fa", emit: pmiren_hairpin
        path "species_ids_db.csv", emit: species_ids

    script:
    """
    ## Create a string to include all files in the 'srnaanno' argument
    srnaanno_files=''
    first='true'
    for file in ${srnaanno}
    do
        if [ "\$first" = "true" ]
        then
            srnaanno_files+="\$file"
            first='false'
        else
            srnaanno_files+=" --srnaanno \$file"
        fi
    done

    ## Do the same for pmiren files
    pmiren_files=''
    first='true'
    for file in ${pmiren}
    do
        if [ "\$first" = "true" ]
        then
            pmiren_files+="\$file"
            first='false'
        else
            pmiren_files+=" --pmiren \$file"
        fi
    done

    ## Do the same for pmiren hairpin files
    pmiren_hairpin_files=''
    first='true'
    for file in ${pmiren_hairpin}
    do
        if [ "\$first" = "true" ]
        then
            pmiren_hairpin_files+="\$file"
            first='false'
        else
            pmiren_hairpin_files+=" --phairpin \$file"
        fi
    done

    # Execute the programm
    01-Database_id_resolver.sh \
        --inputsps "${input_species_names}" \
        --mirbase ${mirbase} \
        --srnaanno \$srnaanno_files \
        --pmiren \$pmiren_files \
        --mhairpin ${mirbase_hairpin} \
        --phairpin \$pmiren_hairpin_files \
        --mirbaseplants ${v_ids_mirbase}

    """
}