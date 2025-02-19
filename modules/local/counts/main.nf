#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to calculate the absolute counts.
process COUNTS {

    tag "$meta.id"
    
    conda "${modulesDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11' :
        'quay.io/biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.raw.tsv"), emit: raw
    
    script:
    def file_d = "$file".toString().replaceAll('.gz$', '')
    """
    # Decompress file 1
    if [[ ${file} == *.gz ]]; then
        gzip -d -f ${file}
    fi

    # Create the count table using bash. It is faster
    echo -e "seq\traw" > ${meta.id}.raw.tsv

    # Calculate absolute counts
    awk 'NR%4==2' "${file_d}" | sort | uniq -c | sort -nr | awk 'BEGIN{FS=" "; OFS="\\t"} {print \$2, \$1}' >> ${meta.id}.raw.tsv
    """
}
