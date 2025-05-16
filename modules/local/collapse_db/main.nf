process COLLAPSE_DB {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9bfad58c74343625d23685a5ea7006c3c154eec7ad85584b8474d7bd8ec956c/data' :
        'community.wave.seqera.io/library/gzip:1.14--19aaa2c84c85ddbc' }"
        
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*collapsed.{fa,fasta}"), emit: coldb

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    #!/bin/bash

    # Decompress input file
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    # Get the file extension
    filename=\$(basename "${fasta_name}")
    extension="\${filename##*.}"

    # Use awk to process the fasta file
    awk '
    BEGIN {
        FS=" ";  # Field separator for headers (assuming spaces between header elements)
    }

    {
        if (\$0 ~ /^>/) {
            header = \$0;
            # Extract the first part of the header (before any space, or use the entire header if no space exists)
            if (index(header, " ") > 0) {
                header_part = substr(header, 2, index(header, " ") - 2);  # Extract first part before the space
            } else {
                header_part = substr(header, 2);  # If no space, take the entire header after ">"
            }
        } else {
            seq = \$0;
            # If the sequence is already in the array, append the new header to the list
            if (seq in sequences) {
                sequences[seq] = sequences[seq] "|" header_part;  # Add header part, separated by "|"
            } else {
                sequences[seq] = header_part;  # Save first part of header
            }
        }
    }

    # At the end, print the result
    END {
        for (seq in sequences) {
            # The header will be the concatenation of all headers, separated by "|"
            print ">" sequences[seq];
            print seq;
        }
    }
    ' "${fasta_name}" > "${prefix}.collapsed.\$extension"

    # Check if the output file is empty
    if [ ! -s "${prefix}.collapsed.\$extension" ]; then
        mv "${prefix}.collapsed.\$extension" "${prefix}.EMPTY_collapsed.\$extension"
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    # Get the file extension
    filename=\$(basename "${fasta_name}")
    extension="\${filename##*.}"

    # Create output file
    touch "${prefix}.collapsed.\$extension"
    """
}