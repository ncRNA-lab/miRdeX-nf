process TSV_TO_FASTA {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"

    input:
    tuple val(meta), path(tsv)
    val col1
    val col2
    val header

    output:
    tuple val(meta), path("*.fa"), emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def use_index_as_id = (col1 == null || col1 < 1)
    def awk_cmd = header ? "NR>1" : "1"

    def awk_script = use_index_as_id
        ? """
        awk -F '\\t' '${awk_cmd} {printf(">sequence%d\\n%s\\n", NR - (${header ? 1 : 0}), \$${col2})}' ${tsv} > ${prefix}.fa
        """
        : """
        awk -F '\\t' '${awk_cmd} {print ">"\$${col1}"\\n"\$${col2}}' ${tsv} > ${prefix}.fa
        """
    """
    ${awk_script}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}

