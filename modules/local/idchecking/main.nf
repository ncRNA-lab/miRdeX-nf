#!/usr/bin/env nextflow

/*
========================================================================================
    DATABASES module
========================================================================================
*/
// Specify DSL2
nextflow.enable.dsl=2

process CHECK_SPECIES_IDS {

    cache 'lenient'

    input:
        val input_species_names
        path ids_table

    output:
        path "species_ids_updated.csv", emit: species_ids

    script:
    """
    change_the_id_using_name(){

        # Arguments
        local species_name="\${1}"
        local species_ids_file="\${2}"
        local id_col="\${3}"

        # Get the two words from the species name
        genus=\$(echo "\$species_name" | awk '{print \$1}')
        specific_epithet=\$(echo "\$species_name" | awk '{print \$2}')

        # Generate the initial code
        new_id=\$(echo "\${genus:0:1}\${specific_epithet:0:2}" | tr '[:upper:]' '[:lower:]')

        # Auto-increment variables
        increment=0
        increment2=0
        sl2_beginning=1

        # Check if the code already exists in the id_table.csv
        while awk -F',' -v col="\$id_col" -v code="\$new_id" '{ if (\$col == code) print }' "\$species_ids_file"  | grep -q .
        do  
            # Check if the variable 'increment' is equal to the length of the specific epithet
            if [ \$increment -eq \$(( \${#specific_epithet} - \$sl2_beginning )) ]
            then
                # Update variables
                let "increment2++"
                increment=0
                sl2_beginning=0
            fi

            # Modify the code if it already exists
            specific_epithet_letter1="\${specific_epithet:0+increment2:1}"
            specific_epithet_letter2="\${specific_epithet:sl2_beginning+increment:1}"

            # Generate the new code
            new_id=\$(echo "\${genus:0:1}\${specific_epithet_letter1}\${specific_epithet_letter2}" | tr '[:upper:]' '[:lower:]')

            # Increment the variable
            let "increment++"
        done

        echo \$new_id
    }

    # Create an output species_ids file
    cp ${ids_table} species_ids_updated.csv

    # Create an array using the input: ['Glycine max', 'Homo sapiens']
    input_species_pre=\$(echo ${input_species_names} | tr -d '[]')
    IFS=',' read -ra input_species_array <<< "\$input_species_pre"
    for species_name in "\${input_species_array[@]}"
    do
        # Remove spaces at the beginning and end of the name
        stripped_sp_name=\$(echo "\$species_name" | sed 's/^ *//; s/ *\$//')
        
        # If the species is not in the databases...
        if ! awk -F',' -v name="\$stripped_sp_name" '\$1 == name' species_ids_updated.csv | grep -q .
        then
            # Get the two words from the species name
            genus=\$(echo "\$stripped_sp_name" | awk '{print \$1}')
            specific_epithet=\$(echo "\$stripped_sp_name" | awk '{print \$2}')

            # Generate an indentifier
            col_ids=5
            new_id=\$(echo "\${genus:0:1}\${specific_epithet:0:2}" | tr '[:upper:]' '[:lower:]')
            species_id=\$(change_the_id_using_name "\$stripped_sp_name" species_ids_updated.csv \$col_ids)

            # Add the species
            echo "\$stripped_sp_name,NULL,NULL,NULL,\$species_id" >> species_ids_updated.csv
        fi
    done

    """
}
