process ISOMIRS_PRECURSOR_CLASSIFICATION {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(mpblast_tsv)
    
    output:
    tuple val(meta), path("*.{EMPTY_canonical,canonical}.tsv")          , emit: canon
    tuple val(meta), path("*.{EMPTY_templated,templated}.tsv")          , emit: templated
    tuple val(meta), path("*.{EMPTY_non_templated,non_templated}.tsv")  , emit: nontemplated

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    
    # Create temporary directory
    workdir=\$(pwd)
    mkdir -p "\$workdir/tmp"

    # Save the canonical miRNAs in a table
    awk -F'\t' '\$14 == \$15' "${mpblast_tsv}" > "${prefix}.canonical.tsv"

    # If there are canonical sequences...
    if [ -s "${prefix}.canonical.tsv" ]; then
        # Save only the sequences in a new TXT file
        cut -f14 "${prefix}.canonical.tsv" | sort -T "\$workdir/tmp" | uniq > \$workdir/tmp/canonical_seqs.txt
        # Remove any row from the original table that contains a canonical sequence.
        awk -F'\t' 'NR==FNR {canon[\$1]; next} !(\$14 in canon)' \$workdir/tmp/canonical_seqs.txt "${mpblast_tsv}"  > "${prefix}.non_canonical.tsv"
    else
        mv "${prefix}.canonical.tsv" "${prefix}.EMPTY_canonical.tsv"
        cat "${mpblast_tsv}" > "${prefix}.non_canonical.tsv"
    fi

    ## 2. Get non-canonical templated sequences
    ################################################################################
    # Save the templated miRNAs in a table
    awk -F'\t' 'index(\$28, \$14) > 0' "${prefix}.non_canonical.tsv" > "${prefix}.templated.tsv"

    # If there are non-canonical templated sequences...
    if [ -s "${prefix}.templated.tsv" ]; then
        # Save only the sequences in a new TXT file
        cut -f14 "${prefix}.templated.tsv" | sort -T "\$workdir/tmp" | uniq > \$workdir/tmp/templated_seqs.txt
        # Remove any row from the original table that contains a templated sequence.
        awk -F'\t' 'NR==FNR {seq[\$1]; next} !(\$14 in seq)' \$workdir/tmp/templated_seqs.txt "${prefix}.non_canonical.tsv" > "${prefix}.non_templated.tsv"
    else
        mv "${prefix}.templated.tsv" "${prefix}.EMPTY_templated.tsv"
        cat "${prefix}.non_canonical.tsv" > "${prefix}.non_templated.tsv"
    fi

    # If there are not non-canonical non-templated sequences...
    if [ ! -s "${prefix}.non_templated.tsv" ]; then
        mv "${prefix}.non_templated.tsv" "${prefix}.EMPTY_non_templated.tsv"
    fi

    rm -r "\$workdir/tmp"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output files
    touch "${prefix}.canonical.tsv"
    touch "${prefix}.templated.tsv"
    touch "${prefix}.non_templated.tsv"
    """
}