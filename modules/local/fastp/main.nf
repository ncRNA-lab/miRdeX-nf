#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to execute fastp
process FASTP {

    tag "$meta.species-$meta.project"

    //conda "${modulesDir}/environment.yml"

    input:
    tuple val(meta), path(fastqgz)
    path adapters
    val min_len
    val max_len

    output:
    tuple val(meta), path('*tr.fastq.gz'), emit: trimmed
    
    script:
    def args = task.ext.args ?: "--length_required $min_len --length_limit $max_len"
    """
    # Get the file name and SRR
    file=\$(basename "${fastqgz}")
    srr=\$(echo "\$file" | awk -F '[_.]' '{print \$1}')
    
    fastp --adapter_fasta ${adapters} \
            -i ${fastqgz} \
            -o ./\$srr"_tr.fastq.gz" \
            $args
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS

    """

    stub:
    """
    touch ${meta.project}_A_tr.fastq.gz
    touch ${meta.project}_B_tr.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
    END_VERSIONS
    """
}
