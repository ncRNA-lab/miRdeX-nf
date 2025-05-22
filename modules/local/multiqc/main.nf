process MULTIQC {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.22.3--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.22.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files)

    output:
    path "*.html"              , emit: report
    path "*_data"              , optional:true, emit: data
    path "*_plots"             , optional:true, emit: plots
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ? "--filename ${meta.id}.html" : ''
    """
    multiqc \\
        --force \\
        $args \\
        $prefix \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}

