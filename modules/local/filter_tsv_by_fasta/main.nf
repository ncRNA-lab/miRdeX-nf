process FILTER_TSV_BY_FASTA {

    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f9/f9bfad58c74343625d23685a5ea7006c3c154eec7ad85584b8474d7bd8ec956c/data' :
        'community.wave.seqera.io/library/gzip:1.14--19aaa2c84c85ddbc' }"
        
    input:
    tuple val(meta), path(tsv), path(fasta)
    val col
    val keep_matches 

    output:
    tuple val(meta), path("*.filtered.tsv")    , emit: filt_tsv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Check if the FASTA file is compressed
    if [[ "${fasta}" == *.gz ]]; then
        gzip -dc "${fasta}" | grep -v '^>' > seq_ids.txt
    else
        grep -v '^>' "${fasta}" > seq_ids.txt
    fi

    awk -F '\t' -v col=${col} -v keep=${keep_matches} '
    BEGIN {
        # Load all sequence IDs into a map
        while ((getline id < "seq_ids.txt") > 0) {
            gsub(/\r/, "", id)      # Remove carriage returns if any
            gsub(/^[ \t]+|[ \t]+\$/, "", id)  # Trim whitespace
            seq[id] = 1
        }
    }
    {
        value = \$col
        if ((keep == "true" && (value in seq)) ||
            (keep == "false" && !(value in seq))) {
            print \$0
        }
    }' "${tsv}" > "${prefix}.filtered.tsv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.filtered.tsv"
    """
}
