process ADD_COUNTS_TO_ISOMIRS_DF {
    
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"

    input:
    tuple val(meta), path(blast_tsv), path(counts_tsv)

    output:
    tuple val(meta), path('*.tsv'), emit: isocounts

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    awk 'FNR==NR { map[\$1]=\$2; next } { seq_value = map[\$14]; print \$0 "\t" seq_value }' "${counts_tsv}" "${blast_tsv}" > "${prefix}.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tsv"
    """


}


