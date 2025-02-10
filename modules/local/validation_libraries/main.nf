process LIBRARIES_VALIDATION {

    tag "$meta.species-$meta.project"

    cpus 2
    //conda "${modulesDir}/environment.yml"

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
    touch ${meta.id}_1.filt.fastq.gz
    touch ${meta.id}_1.sum_projects.tsv'
    touch ${meta.id}_1.sum_libraries.tsv'
    """
}

