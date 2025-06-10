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
    val split_by_column

    output:
    tuple val(meta), path("*.fa"), emit: fasta

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def awk_cmd = header ? "NR>1" : "1"
    def col2_idx = col2
    def col1_idx = col1 ?: 0
    def header_offset = header ? 1 : 0

    if (split_by_column) {
        """
        col2=${col2}
        header=${header_offset}

        ncols=\$(head -n 1 ${tsv} | awk -F '\\t' '{print NF}')

        if [ \$header -eq 1 ]; then
            colnames=(\$(head -n 1 ${tsv} | tr '\\t' '\\n'))
        else
            # Colnames (col1, col2, col3,...)
            colnames=()
            for ((j=1; j<=\$ncols; j++)); do
                colnames+=(col\$j)
            done
        fi

        for ((i=1; i<=ncols; i++)); do
            if [ \$i -ne \$col2 ]; then
                sample_name=\${colnames[\$((i-1))]}
                out_file=\$([ \$header -eq 1 ] && echo "\$sample_name.fa" || echo "${prefix}_\${sample_name}.fa")

                awk -F '\\t' -v col2=\$col2 -v colidx=\$i -v offset=\$header -v sname="\$sample_name" '
                    NR > offset {
                        if (\$colidx > 0)
                            printf(">%s_sequence%d\\n%s\\n", sname, NR - offset, \$col2)
                    }
                ' ${tsv} > \$out_file
            fi
        done
        """
    } else {
        def use_index_as_id = (col1 == null || col1 < 1)
        if (use_index_as_id) {
            """
            awk -F '\\t' -v offset=${header_offset} -v col2=${col2_idx} '
                NR > offset {
                    printf(">sequence%d\\n%s\\n", NR - offset, \$col2)
                }
            ' ${tsv} > ${prefix}.fa
            """
        } else {
            """
            awk -F '\\t' -v offset=${header_offset} -v col1=${col1_idx} -v col2=${col2_idx} '
                NR > offset {
                    print ">"\$col1"\\n"\$col2
                }
            ' ${tsv} > ${prefix}.fa
            """
        }
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa
    """
}

