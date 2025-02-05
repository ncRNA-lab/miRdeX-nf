#!/usr/bin/env nextflow

/*
========================================================================================
    DATABASES module
========================================================================================
*/
// Specify DSL2
nextflow.enable.dsl=2

process DB_ID_RESOLVER {

    input:
        val input_species_names
        path mirbase
        path mirbase_hairpin
        path mirbase_species
        val mirbase_taxon
        path srnaanno
        path pmiren
        path pmiren_hairpin


    output:
        path "species_ids_db.csv"                                , emit: species_ids
        path "01-Mod_databases/miRBase/*mirbase_mature.fa"       , emit: mirbase_mature
        path "01-Mod_databases/miRBase/*mirbase_hairpin.fa"      , emit: mirbase_hairpin
        path "01-Mod_databases/sRNAanno/*srnaanno_mature.fa"     , emit: srnaanno_mature
        path "01-Mod_databases/sRNAanno/*srnaanno_hairpin.fa"    , emit: srnaanno_hairpin
        path "01-Mod_databases/PmiREN/*pmiren_mature.fa"         , emit: pmiren_mature
        path "01-Mod_databases/PmiREN/*pmiren_hairpin.fa"        , emit: pmiren_hairpin

    script:
    def mirbase_validity = mirbase.name != 'EMPTY_mirbase.mature.fa' ? true : false
    def pmiren_validity = pmiren[0].name != 'EMPTY_pmiren.mature.fa' ? true : false
    def srnaanno_validity = srnaanno[0].name != 'EMPTY_srnaanno.miRNA.gff3' ? true : false

    // All databases will be used in the annotation process
    if (mirbase_validity && pmiren_validity && srnaanno_validity) {
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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

        # Obtain the species from miRBase of the requested taxon.
        grep '${mirbase_taxon}' mirbase_species.txt | awk -F'\t' '{print \$2","\$4}' | sort -t',' -k1,1 > mirbase_species_taxon.txt

        # Execute the programm
        01-Database_id_resolver.sh \
            --inputsps "${input_species_names}" \
            --mirbase ${mirbase} \
            --srnaanno \$srnaanno_files \
            --pmiren \$pmiren_files \
            --mhairpin ${mirbase_hairpin} \
            --phairpin \$pmiren_hairpin_files \
            --mirbaseplants mirbase_species_taxon.txt

        """
    // Only miRBase and sRNAanno will be used for the annotation process.
    } else if (mirbase_validity && srnaanno_validity && !pmiren_validity){
        """

        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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

        # Obtain the species from miRBase of the requested taxon.
        grep '${mirbase_taxon}' mirbase_species.txt | awk -F'\t' '{print \$2","\$4}' | sort -t',' -k1,1 > mirbase_species_taxon.txt

        # Execute the programm
        01-Database_id_resolver.sh \
            --inputsps "${input_species_names}" \
            --mirbase ${mirbase} \
            --srnaanno \$srnaanno_files \
            --mhairpin ${mirbase_hairpin} \
            --mirbaseplants mirbase_species_taxon.txt
        
        # Create an empty PmiREN database
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa

        """
    // Only miRBase and PmiREN will be used for the annotation process.
    } else if (mirbase_validity && pmiren_validity && !srnaanno_validity){
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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

        # Obtain the species from miRBase of the requested taxon.
        grep '${mirbase_taxon}' mirbase_species.txt | awk -F'\t' '{print \$2","\$4}' | sort -t',' -k1,1 > mirbase_species_taxon.txt

        # Execute the programm
        01-Database_id_resolver.sh \
            --inputsps "${input_species_names}" \
            --mirbase ${mirbase} \
            --pmiren \$pmiren_files \
            --mhairpin ${mirbase_hairpin} \
            --phairpin \$pmiren_hairpin_files \
            --mirbaseplants mirbase_species_taxon.txt

        # Create empty databases
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa

        """
    // Only PmiREN and sRNAanno will be used for the annotation process.
    } else if (srnaanno_validity && pmiren_validity && !mirbase_validity){
        
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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
            --srnaanno \$srnaanno_files \
            --pmiren \$pmiren_files \
            --phairpin \$pmiren_hairpin_files

        # Create empty databases
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa

        """
    // Only miRBase will be used for the annotation process.
    } else if (mirbase_validity && !srnaanno_validity && !pmiren_validity){
        
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

        # Obtain the species from miRBase of the requested taxon.
        grep '${mirbase_taxon}' mirbase_species.txt | awk -F'\t' '{print \$2","\$4}' | sort -t',' -k1,1 > mirbase_species_taxon.txt

        # Execute the programm
        01-Database_id_resolver.sh \
            --inputsps "${input_species_names}" \
            --mirbase ${mirbase} \
            --mhairpin ${mirbase_hairpin} \
            --mirbaseplants mirbase_species_taxon.txt
        
        # Create empty databases
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa

        """
    // Only PmiREN will be used for the annotation process.
    } else if (pmiren_validity && !mirbase_validity && !srnaanno_validity){
        
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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
            --pmiren \$pmiren_files \
            --phairpin \$pmiren_hairpin_files
        
        # Create empty databases
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa


        """
    } else if (srnaanno_validity && !pmiren_validity && !mirbase_validity){
        
        """
        # Create the output directories
        mkdir -p 01-Mod_databases/miRBase/
        mkdir -p 01-Mod_databases/sRNAanno/
        mkdir -p 01-Mod_databases/PmiREN/

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

        # Execute the programm
        01-Database_id_resolver.sh \
            --inputsps "${input_species_names}" \
            --srnaanno \$srnaanno_files
        
        # Create empty databases
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa

        """
    }
}