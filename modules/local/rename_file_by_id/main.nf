process RENAME_FILE_BY_ID {
    
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("${meta.id}.*"), emit: renamed

    script:
    """
    filename=\$(basename "${file}")
    extension=\${filename##*.}
    cp "${file}" "${meta.id}.\${extension}"
    """
    stub:
    """
    filename=\$(basename "${file}")
    extension=\${filename##*.}
    cp "${file}" "${meta.id}.\${extension}"
    """


}