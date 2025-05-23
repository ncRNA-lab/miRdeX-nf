process COUNTS_VALIDATION {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70135ad1633874b556e06df70a86b8529ba96e2f43e7dcd1284f01489e3e4141/data' :
       'community.wave.seqera.io/library/python_pip_numpy_pandas:5731791ee246815e' }"

    input:
    tuple val(meta), path(counts_files)
    val rep_threshold

    output:
    tuple val(meta), path('*.valid.tsv')    , emit: valid, optional: true
    tuple val(meta), path('*.notvalid.tsv') , emit: notvalid, optional: true
    path '*.sum.tsv'                        , emit: summary

    script:
    """
    02-Validate_counts_matrix.py \
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


