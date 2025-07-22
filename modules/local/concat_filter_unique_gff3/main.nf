process CONCAT_FILTER_UNIQUE_GFF3 {

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

    # Concatenate all files skipping first 4 lines, normalize fields, and remove duplicates
    awk -F"\t" '
    FNR > 4 {

        # Ignore REJECT sequences
        if (\$9 ~ /Filter=REJECT/) next

        # Reemplazar valores por NA
        gsub(/Expression=[^;]+/, "Expression=NA", \$9)
        gsub(/Norm=[^;]+/, "Norm=NA", \$9)
        gsub(/Filter=[^;]+/, "Filter=NA", \$9)

        # Asegurar que los atributos estén separados por "; "
        gsub(/; */, "; ", \$9)
        sub(/; \$/, "", \$9)

        # Reconstruir línea con tabuladores reales
        line = \$1 FS \$2 FS \$3 FS \$4 FS \$5 FS \$6 FS \$7 FS \$8 FS \$9
        if (!seen[line]++) print \$1, \$2, \$3, \$4, \$5,\$6, \$7, \$8, \$9
    }' OFS="\t" ${filesStr} > "${prefix}.valid.gff3"

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
