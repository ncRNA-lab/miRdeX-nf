process COUNTS {

    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9bfad58c74343625d23685a5ea7006c3c154eec7ad85584b8474d7bd8ec956c/data' :
        'community.wave.seqera.io/library/gzip:1.14--19aaa2c84c85ddbc' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.raw.tsv"), emit: raw
    
    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_d = "$file".toString().replaceAll('.gz$', '')
    """
    # Decompress file 1
    if [[ ${file} == *.gz ]]; then
        gzip -d -f ${file}
    fi

    # Create the count table using bash. It is faster
    echo -e "seq\traw" > ${prefix}.raw.tsv

    # Calculate absolute counts
    workdir=\$(pwd)
    mkdir -p "\$workdir/tmp"
    awk 'NR%4==2' "${file_d}" | sort -T "\$workdir/tmp"| uniq -c | sort -nr -T "\$workdir/tmp" | awk 'BEGIN{FS=" "; OFS="\\t"} {print \$2, \$1}' >> ${prefix}.raw.tsv
    rm -r "\$workdir/tmp"
    """
    
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.tsv
    """
}
