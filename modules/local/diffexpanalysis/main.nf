#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process DEA {

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
        --counts ${matrix} \
        --metadata ${meta.metadata} \
        --alpha ${alpha} \
        --min_counts ${min_counts} \
        --min_samples ${min_samples}
    """
}

process DEA_SUM {

    input:
    path sum_dea
    path sum_ea

    output:
    path 'ea_summary.tsv', emit: easum
    path 'dea_summary.tsv', emit: deasum

    script:
    """
    ################################################################################
    #               2. CREATE EXPLORATORY ANALYSIS SUMMARY FILE                    #
    ################################################################################

    # If there is more than one file...
    if [[ \$(ls ${sum_ea} | wc -l) > 1 ]]
    then
        # Concatenate exporatory analysis summary files
        cat \$(ls ${sum_ea} | head -n1) > ea_summary.tsv && tail -n +2 -q ${sum_ea} >> ea_summary.tsv
        # Delete individual species files
        rm -f ${sum_ea}
    # If there is only one file...
    else
        # Rename it
        mv ${sum_ea} ea_summary.tsv
    fi


    ################################################################################
    #           3. CREATE DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY FILE            #
    ################################################################################

    # If there is more than one file...
    if [[ \$(ls ${sum_dea} | wc -l) > 1 ]]
    then
        # Concatenate differential expression analysis summary files
        cat \$(ls ${sum_dea} | head -n1) > dea_summary.tsv && tail -n +2 -q \$(ls ${sum_dea} | tail -n+2) >> dea_summary.tsv
        # Delete individual species files
        rm -f ${sum_dea}
    # If there is only one file...
    else
        # Rename it
        mv ${sum_dea} dea_summary.tsv
    fi
    """
}