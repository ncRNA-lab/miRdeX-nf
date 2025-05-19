process BUILD_MIRNA_EVENT_MATRIX {

    input:
    tuple path(annot_file_list), path(metadata_file_list), val(samples_list)
    val fields

    output:
    path "presence_absence_table.tsv"   , emit: preabs
    path "shrunken_log2fc_table.tsv"    , emit: log2fc
    path "id_correspondence.tsv"        , emit: ids
    path  "versions.yml"                , emit: versions

    script:
    def annotList = annot_file_list.toList()
    def metaList = metadata_file_list.toList()
    def samplesList = samples_list.toList()
    def args = (0..<annotList.size()).collect { i ->
        "--annotation \"${annotList[i]}\" --metadata \"${metaList[i]}\" --samples \"${samplesList[i]}\""
    }.join(' ')
    """
    08-miRNA_event_matrix_builder.py \\
        ${args} \\
        --fields ${fields}

    # Create versions file
    echo "${task.process}:" > versions.yml
    echo "    python: \$(python3 -c 'import platform; print(platform.python_version())')" >> versions.yml
    echo "    pandas: \$(python3 -c 'import pandas as pd; print(pd.__version__)')" >> versions.yml
    echo "    natsort: \$(python3 -c 'import natsort; print(natsort.__version__)')" >> versions.yml
    """
}