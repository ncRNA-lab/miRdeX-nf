#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to execute 03-Filter_by_depth_rep.py
process LIBRARIES_VALIDATION {

    tag "$meta.species-$meta.project"

    cpus 2
    //conda "${modulesDir}/environment.yml"

    input:
    tuple val(meta), path(fastq_files)
    val depth_threshold
    val rep_threshold

    output:
    tuple val(meta), path('*.valid.fastq.gz')        , emit: valid, optional: true
    tuple val(meta), path('*.notvalid.fastq.gz')   , emit: notvalid, optional: true
    path '*.sum_projects.tsv'                       , emit: sumprojects
    path '*.sum_libraries.tsv'                      , emit: sumlibraries

    script:
    """
    03-Validate_libraries.py \
        -i ${fastq_files} \
        -j ${meta.project} \
        -s ${meta.species_id} \
        -m ${meta.metadata} \
        -d ${depth_threshold} \
        -r ${rep_threshold} \
        -p ${task.cpus}
    """

    stub:
    """
    touch ${meta.project}_1.filt.fastq.gz
    touch ${meta.project}_1.sum_projects.tsv'
    touch ${meta.project}_1.sum_libraries.tsv'
    """
}

// Define the process to summarize the results obtained by 03-Filter_by_depth_rep.py
process VALIDATION_SUM {

    input:
    path sum_projects
    path sum_libraries

    output:
    path 'sum_projects.tsv', emit: projects
    path 'sum_libraries.tsv', emit: libraries
    
    shell:
    """
    # Concatenate files with the results of each specie
    awk '{ print \$1\"\t\"\$2\"\t\"\$3\"\t\"\$4\"\t\"\$5\"\t\"\$6}' !{sum_projects} > sum_projects_temp.tsv
    sort -t'\t' -k1,1 -k2  -n sum_projects_temp.tsv > sum_projects.tsv
    sed -i '1s/^/Specie\tProject\tGroup\tFiltered\tExcluded\tValidity\\n/' sum_projects.tsv

    # Concatenate summary files of filtered and excluded libraries
    [ `ls -1 !{sum_libraries} 2>/dev/null | wc -l` -gt 0 ] && awk '{print}' !{sum_libraries} > sum_libraries_temp.tsv
    sort -t'\t' -k1 -n sum_libraries_temp.tsv > sum_libraries.tsv
    sed -i '1s/^/File\tSequencing_depth\tDepth_sequencing_filter\tRep_threshold_filter\\n/' sum_libraries.tsv

    # Remove specie files
    rm -f *_sum_projects.tsv
    rm -f *_sum_libraries.tsv
    rm -f *_temp.tsv
    """
    stub:
    """
    touch sum_projects.tsv
    touch sum_libraries.tsv
    """
}


