process CONCAT_TSV {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(tsv_files)
    val header

    output:
    tuple val(meta), path("*.tsv"), emit: concat

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create a Bash array from the list of input files
    files=( ${tsv_files} )

    if [[ "${header}" == "true" ]]; then
        {
            head -n 1 "\${files[0]}"
            for f in "\${files[@]:1}"; do
                tail -n +2 "\$f"
            done
        } > "${prefix}.tsv"
    else
        cat "\${files[@]}" > "${prefix}.tsv"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output file
    touch "${prefix}.tsv"
    """
}
