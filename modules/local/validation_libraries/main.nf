process LIBRARIES_VALIDATION {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/70/70135ad1633874b556e06df70a86b8529ba96e2f43e7dcd1284f01489e3e4141/data' :
        'community.wave.seqera.io/library/python_pip_numpy_pandas:5731791ee246815e' }"

    input:
    tuple val(meta), path(fastq_files)
    val depth_threshold
    val rep_threshold

    output:
    tuple val(meta), path('*.valid.fastq.gz')       , emit: valid, optional: true
    tuple val(meta), path('*.notvalid.fastq.gz')    , emit: notvalid, optional: true
    path '*.sum_projects.tsv'                       , emit: sumprojects
    path '*.sum_libraries.tsv'                      , emit: sumlibraries

    script:
    """
    03-Validate_libraries.py \
        -i ${fastq_files} \
        -j ${meta.id} \
        -m ${meta.metadata} \
        -d ${depth_threshold} \
        -r ${rep_threshold} \
        -p ${task.cpus}
    """

    stub:
    """
    # Iterate through the input files
    for file in ${fastq_files}; do
        # Get the input file name
        filename=\$(basename \$file)
        name=\${filename%%.*}

        # Create a valid empty file
        touch \$name.valid.fastq.gz

        # Create the libraries summary file
        echo -e "\$name.fastq.gz\t\$name\tvalid\tvalid" >> ${meta.id}.sum_libraries.tsv
    done

    #!/bin/bash

    # Create the project summary files
    echo -e "${meta.id}\t1\t50\t10\tvalid" > ${meta.id}.sum_projects.tsv
    """
}

