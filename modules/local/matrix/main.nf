#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Build the miRNA-vs-Stress matrix (presence-absence)
process PREABS_MATRIX {

    input:
    path sum_dea_file
    path fam_annot_files

    output:
    path "*_table.tsv", emit: preabs

    script:
    """
    08-Build_miRNA_vs_stress_tables.py \
        --path-dea ${sum_dea_file} \
        --annot-files ${fam_annot_files}
    """

}
