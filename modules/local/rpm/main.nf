process RPM {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.rpm.tsv"), emit: rpm
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create the count table using bash. It is faster
    echo -e "seq\trpm" > ${prefix}.rpm.tsv

    # Perform the sum of all absolute counts (skipping the header)
    total=\$(awk 'BEGIN {FS="\\t"; total=0; getline} {total+=\$2} END {print total}' "${counts}")

    # Calculate the rpm (skipping the header)
    awk -v total="\${total}" 'BEGIN {FS="\\t"; OFS="\\t"; getline} {rpm = (\$2 / total) * 1000000; print \$1, rpm}' "${counts}" >> ${prefix}.rpm.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output file
    touch ${prefix}.rpm.tsv
    """
}
