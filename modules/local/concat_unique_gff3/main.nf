process CONCAT_UNIQUE_GFF3 {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(gff_files)

    output:
    tuple val(meta), path("*.concat.gff3"), emit: gff3

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filesStr = gff_files.join(' ')
    """
    # Extract header lines from the first file, filtering needed lines
    header=\$(grep '^##' "${gff_files[0]}" | grep -E 'mirGFF3|source-ontology|TOOLS')

    # Append COLDATA line with meta.id
    header="\$header\n## COLDATA: ${meta.id}"

    # Concatenate all files skipping first 4 lines, remove duplicates
    awk 'FNR > 4' ${filesStr} | awk '!seen[\$0]++' > "${prefix}.concat.nodup.gff3"

    # Remove rejected sequences
    grep -v 'Filter=REJECT' "${prefix}.concat.nodup.gff3" > "${prefix}.valid.gff3"

    # Write header and deduplicated content to final file
    {
    echo "\$header"
    cat "${prefix}.valid.gff3"
    } > "${prefix}.concat.gff3"
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output file
    touch "${prefix}.concat.gff3"
    """
}
