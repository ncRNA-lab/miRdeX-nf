#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DIFFEXPANALYSIS {

    // Process tag
    tag "$meta.id"

    input:
    tuple val(meta), path(matrix)
    val alpha
    val min_counts
    val min_samples

    output:
    tuple val(meta), path("*_raw.tsv")  , emit: raw
    tuple val(meta), path("*_sig.tsv")  , emit: sig
    path "${meta.id}.ea_summary.tsv"    , emit: easum
    path "${meta.id}.dea_summary.tsv"   , emit: deasum

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
}