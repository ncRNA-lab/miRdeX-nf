process COUNTS_VALIDATION {

    tag "$meta.species-$meta.project"
    //conda "${modulesDir}/environment.yml"

    input:
    tuple val(meta), path(counts_files)
    val rep_threshold

    output:
    tuple val(meta), path('*.valid.tsv')    , emit: valid, optional: true
    tuple val(meta), path('*.notvalid.tsv') , emit: notvalid, optional: true
    path '*.sum.tsv'                        , emit: summary

    script:
    """
    03-Validate_counts_matrix.py \
        -i ${meta.id} \
        -g ${meta.group_id} \
        -c ${counts_files} \
        -m ${meta.metadata} \
        -r ${rep_threshold} \
    """

    stub:
    """
    touch ${meta.project}.counts_validity.tsv
    """
}


