#!/usr/bin/env nextflow

/*
========================================================================================
    DATABASES module
========================================================================================
*/
// Specify DSL2
nextflow.enable.dsl=2

process DB_ID_RESOLVER {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4e/4e63ee3931d2e9db0b876dcbb3afaa456e57251a665edec3d44811990f06fa55/data' :
        'community.wave.seqera.io/library/seqkit:2.9.0--e0e29e1f5c28842a' }"

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

    stub:
    def mirbase_validity = mirbase.name != 'EMPTY_mirbase.mature.fa' ? true : false
    def pmiren_validity = pmiren[0].name != 'EMPTY_pmiren.mature.fa' ? true : false
    def srnaanno_validity = srnaanno[0].name != 'EMPTY_srnaanno.miRNA.gff3' ? true : false

    """
    # Create the output directories
    mkdir -p 01-Mod_databases/miRBase/
    mkdir -p 01-Mod_databases/sRNAanno/
    mkdir -p 01-Mod_databases/PmiREN/

    # Define booleans for the validity of each database
    mirbase_validity=${mirbase_validity}
    pmiren_validity=${pmiren_validity}
    srnaanno_validity=${srnaanno_validity}

    # Define the databases and their specific species
    declare -A species_db=(
        ["miRBase"]="ath bna bol"    # Species for miRBase
        ["sRNAanno"]="ath bna bol"   # Species for sRNAanno
        ["PmiREN"]="ath bna bol"     # Species for PmiREN
    )

    # Function to generate files for a given database
    generate_files() {
        # Arguments
        local db=\$1

        # Define dummy sequences
        mature_seq="ATGCATGCATGCATGCATGC"
        hairpin_seq="TGCAATGCAATGCAAAGCTAAAGGTTTGCATAGCTG"
        
        # Define the databases and their specific species
        declare -A species_db=(
            ["miRBase"]="ath bna bol"    # Species for miRBase
            ["sRNAanno"]="ath tae bol"   # Species for sRNAanno
            ["PmiREN"]="gma bna ttu"     # Species for PmiREN
        )

        # Create necessary directories
        mkdir -p "01-Mod_databases/\$db"

        # Get the species for the current database
        species_list=(\${species_db[\$db]})

        # Generate sequences for each species in the database
        for species_code in "\${species_list[@]}"; do
            echo \$species_code
            # Create the sequences for mature.fa
            echo -e ">\$species_code-miR1\n\$mature_seq" >> "01-Mod_databases/\$db/\${db,,}_mature.fa"
            echo -e ">\$species_code-miR2\n\$mature_seq" >> "01-Mod_databases/\$db/\${db,,}_mature.fa"

            # Create the sequences for hairpin.fa
            echo -e ">\$species_code-miR1\n\$hairpin_seq" >> "01-Mod_databases/\$db/\${db,,}_hairpin.fa"
            echo -e ">\$species_code-miR2\n\$hairpin_seq" >> "01-Mod_databases/\$db/\${db,,}_hairpin.fa"
        done
    }

    # Only generate files for the databases that are valid (true)
    if [ "\$mirbase_validity" = true ]; then
        # Pass the species for miRBase as arguments to generate_files
        generate_files "miRBase" \${species_db["miRBase"]}
    else
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_mature.fa
        touch 01-Mod_databases/miRBase/EMPTY_mirbase_hairpin.fa
    fi

    if [ "\$srnaanno_validity" = true ]; then
        # Pass the species for sRNAanno as arguments to generate_files
        generate_files "sRNAanno" \${species_db["sRNAanno"]}
    else
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_mature.fa
        touch 01-Mod_databases/sRNAanno/EMPTY_srnaanno_hairpin.fa
    fi

    if [ "\$pmiren_validity" = true ]; then
        # Pass the species for PmiREN as arguments to generate_files
        generate_files "PmiREN" \${species_db["PmiREN"]}
    else
        # Create an empty PmiREN database
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_mature.fa
        touch 01-Mod_databases/PmiREN/EMPTY_pmiren_hairpin.fa
    fi

    # Create species_ids_db.csv file
    echo "species,mirbase,pmiren,srnaanno,final_id" > "species_ids_db.csv"
    echo "Arabidopsis thaliana,ath,NULL,ath,ath" >> "species_ids_db.csv"
    echo "Glycine max,NULL,gma,NULL,gma" >> "species_ids_db.csv"
    echo "Brassica napus,bna,bna,NULL,bna" >> "species_ids_db.csv"
    echo "Brassica olaracea,bol,NULL,bol,bol" >> "species_ids_db.csv"
    echo "Triticum aestivum,NULL,NULL,tae,tae" >> "species_ids_db.csv"
    echo "Triticum turgidum,NULL,ttu,NULL,ttu" >> "species_ids_db.csv"

    """

}