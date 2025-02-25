process COUNTS_MATRIX {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70135ad1633874b556e06df70a86b8529ba96e2f43e7dcd1284f01489e3e4141/data' :
       'community.wave.seqera.io/library/python_pip_numpy_pandas:5731791ee246815e' }"

    input:
    tuple val(meta), path(counts)
    val type

    output:
    tuple val(meta), path("${meta.id}_*.{raw,rpm}.tsv"), emit: matrix

    script:
    """
    if [ "${type}" == "raw" ]; then
    04-Create_counts_matrix.py \
        --project ${meta.id} \
        --counts-tsv ${counts} \
        --metadata ${meta.metadata} \
        --valid-groups ${meta.valid_groups}
    else
    04-Create_counts_matrix.py \
        --project ${meta.id} \
        --counts-tsv ${counts} \
        --metadata ${meta.metadata} \
        --valid-groups ${meta.valid_groups} \
        --rpm
    fi
    """
}

// """
// if [ "${type}" == "raw" ]; then
// 04-Create_counts_matrix.py \
//     --project ${meta.id} \
//     --counts-tsv ${counts} \
//     --metadata ${meta.metadata} \
//     --valid-groups ${meta.valid_groups}
// else
// 04-Create_counts_matrix.py \
//     --project ${meta.id} \
//     --counts-tsv ${counts} \
//     --metadata ${meta.metadata} \
//     --valid-groups ${meta.valid_groups} \
//     --rpm
// fi
// """