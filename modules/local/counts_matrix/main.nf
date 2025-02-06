process COUNTS_MATRIX {

    tag "$meta.id"
    conda "${modulesDir}/environment.yml"
    debug true

    input:
    tuple val(meta), path(counts)
    val type

    output:
    tuple val(meta), path("${meta.id}_*.{raw,rpm}.tsv"), emit: matrix

    script:
    // RPM matrix
    if (type == 'rpm') {
        """
        04-Create_counts_matrix.py \
            --project ${meta.id} \
            --counts-tsv ${counts} \
            --metadata ${meta.metadata} \
            --valid-groups ${meta.valid_groups} \
            --rpm
        """
    // Raw counts matrix
    } else {
        """
        04-Create_counts_matrix.py \
            --project ${meta.id} \
            --counts-tsv ${counts} \
            --metadata ${meta.metadata} \
            --valid-groups ${meta.valid_groups}
        """
    }
}