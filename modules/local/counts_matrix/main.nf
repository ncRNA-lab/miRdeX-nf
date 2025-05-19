process COUNTS_MATRIX {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70135ad1633874b556e06df70a86b8529ba96e2f43e7dcd1284f01489e3e4141/data' :
       'community.wave.seqera.io/library/python_pip_numpy_pandas:5731791ee246815e' }"

    input:
    tuple val(meta), path(counts)
    val type
    val counts_project_matrix

    output:
    tuple val(meta), path("${meta.id}*.{raw,rpm}.tsv"), emit: matrix
    path  "versions.yml"                              , emit: versions

    script:
    def create_project_matrix = counts_project_matrix ? "--project-table" : ""
    """
    if [ "${type}" == "raw" ]; then
    04-Create_counts_matrix.py \
        --project ${meta.id} \
        --counts-tsv ${counts} \
        --metadata ${meta.metadata} \
        --valid-groups ${meta.valid_groups} \
        ${create_project_matrix}
    else
    04-Create_counts_matrix.py \
        --project ${meta.id} \
        --counts-tsv ${counts} \
        --metadata ${meta.metadata} \
        --valid-groups ${meta.valid_groups} \
        --rpm \
        ${create_project_matrix}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        numpy: \$(python3 -c "import numpy as np; print(np.__version__)")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """
    stub:
    """
    # Extract the filenames without the extension
    columns="seq"
    for file in ${counts}; do
        name=\$(echo "\$file" | cut -d'.' -f1)  # Remove the file extension
        columns="\$columns	\$name"  # Add filenames as columns
    done

    # Generate the output filename dynamically
    output_file="${meta.id}_1.${type}.tsv"
    echo -e "\$columns" > "\$output_file"

    # Function to generate random DNA sequences
    generate_dna_sequence() {
        local length=\$((RANDOM % 6 + 20))  # Random length between 20 and 25
        local sequence=""
        for i in \$(seq 1 \$length); do
            sequence="\$sequence\$(echo 'ATCG' | fold -w1 | shuf -n1)"  # Choose random nucleotide
        done
        echo "\$sequence"
    }

    # Generate fake count data for 10 rows
    for i in {1..10}; do
        seq_id=\$(generate_dna_sequence)  # Generate a random DNA sequence
        row="\$seq_id"  # Add the sequence in the first column
        for file in ${counts}; do
            row="\$row	\$((RANDOM % 100))"  # Generate random count values
        done
        echo -e "\$row" >> "\$output_file"
    done
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        numpy: \$(python3 -c "import numpy as np; print(np.__version__)")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """

}