process DEA_TO_FASTA {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(dea)

    output:
    tuple val(meta), path("*.fasta"), emit: fasta
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'NR>1 {printf(">sequence%d\\n%s\\n", NR-1, \$1)}' ${dea} > ${prefix}.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fasta
    """
}
