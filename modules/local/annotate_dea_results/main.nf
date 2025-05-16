process ANNOTATE_DEA_RESULTS {

    // Process tag
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://antoglz/mirnas_analysis/annotate_dea_results:latest' :
        'docker.io/antoglz/annotate_dea_results:latest' }"

    input:
    tuple val(meta), path(dea_file), path(annot_file)
    val classes

    output:
    tuple val(meta), path("*.all.tsv")      , emit: all
    tuple val(meta), path("*.unique.tsv")   , emit: unique
    tuple val(meta), path("*.boxplot.png")  , emit: boxplot
    path "*.summary.tsv"                    , emit: fam_sum

    script:
    """
    Annotate_dea_results.r \
        --id ${meta.id} \
        --dea ${dea_file} \
        --annotation ${annot_file} \
        --classes ${classes}
    """
}
