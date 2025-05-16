process MERGE_AND_FILTER_MATURE_PRECURSOR_BLAST {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(mature), path(precursor)

    output:
    tuple val(meta), path("*mpblast.tsv"), emit: mpblast

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/bin/bash
    
    ### FUNCTIONS
    expand_and_clean_tsv (){
        
        # Arguments
        local input_file="\${1}"
        local ouput_file="\${2}"

        # Expand and clean tsv table
        awk -F '\t' '{
            split(\$2, ids, "|");
            for (i in ids) {
                id = ids[i];
                gsub(/-?(Known|5p|3p|mature|star)-?/, "-", id);
                gsub(/^-+|-+\$/, "", id);
                gsub(/--+/, "-", id);
                print \$1, id, \$0;
            }
        }' OFS='\t' "\$input_file" > "\$ouput_file"

    }

    create_key_column (){
            
        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        awk -F'\t' '{
            key = tolower(\$1 "|" \$2);
            print key, \$0;
        }' OFS='\t' "\$input_file" > "\$output_file"

    }

    select_most_likely_references(){

        # Arguments
        local input_file="\${1}"
        local output_file="\${2}"

        awk -F'\\t' '
        {
            key = \$14;
            val = \$12;

            if (!(key in min)) {
                min[key] = val;
                lines[key] = \$0;
            } else if (val < min[key]) {
                min[key] = val;
                lines[key] = \$0;
            } else if (val == min[key]) {
                lines[key] = lines[key] "\\n" \$0;
            }
        }
        END {
            for (k in lines) {
                print lines[k];
            }
        }
        '  "\$input_file" > "\$output_file"
    }

    ### MAIN
    main () {
        # Create temporary directory
        mkdir -p tmp

        # Expand and clean TSV tables (Prepare mature miRNA identifiers to find their precursors)
        expand_and_clean_tsv "${mature}" tmp/mature_expanded_cleaned.tsv
        expand_and_clean_tsv "${precursor}" tmp/precursor_expanded_cleaned.tsv

        # Create the key column for the files to be merged
        create_key_column tmp/mature_expanded_cleaned.tsv tmp/mature_with_key.tsv
        create_key_column tmp/precursor_expanded_cleaned.tsv tmp/precursor_with_key.tsv

        # Sort the TSV files to be merged
        sort -k1,1 tmp/mature_with_key.tsv > tmp/mature_sorted.tsv
        sort -k1,1 tmp/precursor_with_key.tsv > tmp/precursor_sorted.tsv

        # Inner Join
        join -t  \$'\t' -1 1 -2 1 tmp/mature_sorted.tsv tmp/precursor_sorted.tsv > tmp/joined_table.tsv

        # Remove undesired columns
        cut --complement -f1,4,18,20,32 tmp/joined_table.tsv > tmp/isomirs_all.tsv
        
        # Select most likely references for each sequence
        select_most_likely_references tmp/isomirs_all.tsv "${prefix}.mpblast.tsv"

        # Add empty suffix if the file is empty
        if [ ! -s "${prefix}.mpblast.tsv" ]; then
            mv "${prefix}.mpblast.tsv" "${prefix}.EMPTY_mpblast.tsv"
        fi
    }
    main "\$@"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output file
    touch "${prefix}.mpblast.tsv"
    """
}