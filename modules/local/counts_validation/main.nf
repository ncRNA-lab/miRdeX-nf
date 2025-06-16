process COUNTS_VALIDATION {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70135ad1633874b556e06df70a86b8529ba96e2f43e7dcd1284f01489e3e4141/data' :
       'community.wave.seqera.io/library/python_pip_numpy_pandas:5731791ee246815e' }"

    input:
    tuple val(meta), path(counts), path(metadata), val(group_id)
    val rep_threshold

    output:
    tuple val(meta), path('*.valid.tsv')    , emit: valid, optional: true
    tuple val(meta), path('*.notvalid.tsv') , emit: notvalid, optional: true
    path '*.sum.tsv'                        , emit: summary
    path  "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    02-Validate_counts_matrix.py \
        -i ${prefix} \
        -g ${group_id} \
        -c ${counts} \
        -m ${metadata} \
        -r ${rep_threshold} \

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        numpy: \$(python3 -c "import numpy as np; print(np.__version__)")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_${group_id}.valid.tsv
    touch ${prefix}_${group_id}.notvalid.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        numpy: \$(python3 -c "import numpy as np; print(np.__version__)")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
    END_VERSIONS
    """
}


