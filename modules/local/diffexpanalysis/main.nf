#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DIFFEXPANALYSIS {

    // Process tag
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://antoglz/mirnas_analysis/diffexp:latest' :
        'docker.io/antoglz/diffexp:latest' }"

    input:
    tuple val(meta), path(matrix)
    val alpha
    val min_counts
    val min_samples

    output:
    tuple val(meta), path("*_raw.tsv")  , emit: raw
    tuple val(meta), path("*_sig.tsv")  , emit: sig
    tuple val(meta), path("${meta.id}.ea_summary.tsv"), emit: easum
    tuple val(meta), path("${meta.id}.dea_summary.tsv")   , emit: deasum

    script:
    """
    05-Diff_exp_analysis.r \
        --id ${meta.id} \
        --group_id ${meta.group_id} \
        --counts ${matrix} \
        --metadata ${meta.metadata} \
        --alpha ${alpha} \
        --min_counts ${min_counts} \
        --min_samples ${min_samples}
    """

    stub:
    """
    # CAMBIAR AL TERMINAR LA PRUEBA DE RESUME
    touch ${meta.project}_3.${type}.tsv
    touch ${meta.project}_3.${type}.tsv
    touch ${meta.project}_3.${type}.tsv
    """
}