#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to execute fastqc
process FASTQC {

    tag "$meta.species-$meta.project"

    //conda "${modulesDir}/environment.yml"

    input:
    tuple val(meta), path(fastq)
    val type

    output:
    tuple val(meta), path('*.zip') , emit : zip
    tuple val(meta), path('*.html'), emit : html
    path  "versions.yml"           , emit : versions

    when:
    meta.size() > 0 

    script:
    def args = task.ext.args ?: ''
    """
    # Execute fastqc
    fastqc \\
        $args \\
        ${fastq}

    # Add the suffix 'raw' or 'trimmed' to the output files
    for file in *.{html,zip}
    do
        # Get the file extension
        extension="\${file##*.}"
        # Rename the file
        mv "\$file" "${type}_\${file%.*}.\$extension"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | grep -e 'FastQC' | sed 's/^.*FastQC //' )
    END_VERSIONS

    """
    
    stub:
    """
    touch ${meta.project}_1.zip
    touch ${meta.project}_2.zip
    touch ${meta.project}_1.html
    touch ${meta.project}_2.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | grep -e 'FastQC' | sed 's/^.*FastQC //' )
    END_VERSIONS
    """
}