process BUILD_MIRNA_EVENT_MATRIX {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/09/09d962c742a1be43eff7a362658618cec47d333f0f09b34f493fc960b05a27df/data' :
        'community.wave.seqera.io/library/python_pip_natsort_pandas:4aa565244c22e36b' }"

    input:
    tuple path(annot_file_list), path(metadata_file_list), val(samples_list)
    val fields

    output:
    path "presence_absence_table.tsv"   , emit: preabs
    path "shrunken_log2fc_table.tsv"    , emit: log2fc
    path "ids.tsv"                      , emit: ids
    path "comparisons_ids.tsv"          , emit: compids
    path  "versions.yml"                , emit: versions

    script:
    def files = annot_file_list.flatten().join(' ')
    def meta_files = metadata_file_list.flatten().join(' ')
    def samples = samples_list.flatten().join(' ')
    """
    # Input
    arr_files=(${files})
    arr_meta_files=(${meta_files})
    arr_samples=(${samples})

    # Prepare arguments string
    args=""
    for i in "\${!arr_files[@]}"; do
        args_element="--annotation \${arr_files[\$i]} --metadata \${arr_meta_files[\$i]} --samples \${arr_samples[\$i]}"
        args+="\$args_element "
    done

    # Run the python script
    07-miRNA_event_matrix_builder.py \\
        \$args \\
        --fields ${fields}

    # Versions
    echo "${task.process}:" > versions.yml
    echo "    python: \$(python3 -c 'import platform; print(platform.python_version())')" >> versions.yml
    echo "    pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')" >> versions.yml
    echo "    natsort: \$(python3 -c 'import natsort; print(natsort.__version__)')" >> versions.yml
    """

    stub:
    """
    # Create output files
    touch presence_absence_table.tsv
    touch shrunken_log2fc_table.tsv
    touch ids.tsv
    touch comparisons_ids.tsv

    # Create versions file
    echo "${task.process}:" > versions.yml
    echo "    python: \$(python3 -c 'import platform; print(platform.python_version())')" >> versions.yml
    echo "    pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')" >> versions.yml
    echo "    natsort: \$(python3 -c 'import natsort; print(natsort.__version__)')" >> versions.yml
    """
}