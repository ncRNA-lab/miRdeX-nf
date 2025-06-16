process SPLIT_DB_BY_SPECIES {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), val(species_id), path(database)

    output:
    tuple val(meta), path("*.{${species_id},EMPTY}.{fa,fasta}"), emit: species_db

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # get the file extension
    filename=\$(basename "${database}")
    extension="\${filename##*.}"

    # Filter the database by species
    awk -v id="${species_id}" '/^>/{p=(\$0 ~ "^>" id)} p' ${database} > ${prefix}.${species_id}.\$extension

    # Check if the output file is empty
    if [ ! -s "${prefix}.${species_id}.fa" ]; then
        mv ${prefix}.${species_id}.fa ${prefix}.EMPTY.\$extension
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # get the file extension
    filename=\$(basename "${database}")
    extension="\${filename##*.}"

    # Create output file
    touch "${prefix}.${species_id}.\$extension"
    """
}