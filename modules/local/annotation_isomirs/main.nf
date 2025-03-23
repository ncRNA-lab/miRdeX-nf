#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process MIRNA_ANNOTATION {

    // Process tag
    tag "$meta.id"

    debug true

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/56/5679100c87339bf761da90857a16e6ced1865542dac922cad1a460396f3f0bd1/data' :
        'community.wave.seqera.io/library/blast_python_pip_biopython_pandas:89b2eff8957e8c23' }"

    input:
    tuple val(meta), path(dea_file)

    output:
    tuple val(meta), path("*.dea_isomirs.tsv")      , emit: isomirs, optional: true
    tuple val(meta), path("*.blast_isomirs.tsv")    , emit: blast, optional: true
    tuple val(meta), path("*.summary_isomirs.tsv")  , emit: summary, optional: true

    script:
    """
    # Run the required script
    06-IsomiRs_annotation.py \
        --id ${meta.id} \
        --input ${dea_file} \
        --species-id ${meta.species_id} \
        --mature-db ${meta.annot_mature_db} \
        --hairpin-db ${meta.annot_hairpin_db}
    """

    stub:
    """
    # Create the main output files
    touch ${meta.id}.annot_all.tsv
    touch ${meta.id}.annot_filt.tsv

    # Create the summary files
    file_name=\$(basename "${mirbase}")
    reference_name=\$(echo "\$file_name" | grep -oE "mature|hairpin")
    touch ${meta.id}.\$reference_name"_summary.tsv"
    touch ${meta.id}.\$reference_name"_summary_len.tsv"
    """
}