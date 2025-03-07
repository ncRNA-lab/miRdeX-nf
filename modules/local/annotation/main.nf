#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

process MIRNA_ANNOTATION {

    // Process tag
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/85/85c083c5e9b2c93a87c941de745f9c7a0497a842219ad328dc4ec11afd5753cb/data' :
        'community.wave.seqera.io/library/bowtie_bc:1b7543aadb5dcbcd' }"

    input:
    tuple val(meta), path(dea_files), path(ea_file)
    path mirbase
    path pmiren
    path srnaanno
    val mismatches
    val mww_pvalue
    val min_num_db

    output:
    tuple val(meta), path("*.annot_all.tsv"), emit: annot, optional: true
    tuple val(meta), path("*.annot_filt.tsv"), emit: annotfilt, optional: true
    path "*_summary.tsv", emit: sum
    path "*_summary_len.tsv", emit: sumlen, optional: true

    script:
    def mirbase_in = mirbase.name != 'EMPTY_mirbase_mature.fa' ? "--mirbase ${mirbase}" : ""
    def srnaanno_in = srnaanno.name != 'EMPTY_srnaanno_mature.fa' ? "--srnaanno ${srnaanno}" : ""
    def pmiren_in = pmiren.name != 'EMPTY_pmiren_mature.fa' ? "--pmiren ${pmiren}" : ""
    //PROBAR ESTO PARA PERMITIR HAIRPIN ADEMAS DE MATURE
    //def mirbase_in = mirbase.name.matches('EMPTY_.*\\.fa') ? "" : "--mirbase ${mirbase}"
    //def srnaanno_in = srnaanno.name.matches('EMPTY_.*\\.fa') ? "" : "--srnaanno ${srnaanno}"
    //def pmiren_in = pmiren.name.matches('EMPTY_.*\\.fa') ? "" : "--pmiren ${pmiren}"
    """
    ## Create a string to include all files in the 'input' argument
    input_files=''
    first='true'
    for file in ${dea_files}
    do
        if [ "\$first" = "true" ]
        then
            input_files+="\$file"
            first='false'
        else
            input_files+=" --input \$file"
        fi
    done

    ## Execute the programm
    06-miRNAs_annotation.sh \
        --input \$input_files \
        --id ${meta.id} \
        --species ${meta.species_id} \
        $mirbase_in \
        $pmiren_in \
        $srnaanno_in \
        --mismatches ${mismatches} \
        --threads ${task.cpus} \
        --ea-table ${ea_file} \
        --mww-pvalue ${mww_pvalue} \
        --min-num-db ${min_num_db} 
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

process MIRNA_ANNOTATION_SUM {

    debug true

    input:
    path sum_mat_files
    path sum_mat_len_files
    path sum_hair_files
    path sum_hair_len_files

    output:
    path "summary_annot.tsv"
    path "summary_len.tsv"

    script:
    """
    # Concatenate summary files
    cat ${sum_mat_files} > sum_mat.tsv
    cat ${sum_mat_len_files} > sum_mat_len.tsv
    cat ${sum_hair_files} > sum_hair.tsv
    cat ${sum_hair_len_files} > sum_hair_len.tsv

    # Sort files to be joined
    LANG=en_EN sort -k 2 -t \$'\\t' sum_mat.tsv -o mature_summary_sort.tsv
    LANG=en_EN sort -k 2 -t \$'\\t' sum_hair.tsv -o hairpin_summary_sort.tsv

    # Join summary files
    LANG=en_EN join -1 2 -2 2 -t \$'\\t' \
            -o 1.1,1.2,1.3,1.4,2.3,2.4  \
            mature_summary_sort.tsv \
            hairpin_summary_sort.tsv  > summary_tmp.tsv
    LANG=en_EN sort -k 1 -t \$'\\t' summary_tmp.tsv -o summary_annot.tsv
    sed -i '1iSpecies\tStress_event\tAnnot_miRNAs\tAnnot_miRNAs_filtered\tAnnot_precursor\tAnnot_precursor_filtered' summary_annot.tsv

    # Sort files to be joined
    LANG=en_EN sort -k 2 -t \$'\\t' sum_mat_len.tsv -o mature_summary_len_sort.tsv
    LANG=en_EN sort -k 2 -t \$'\\t' sum_hair_len.tsv -o hairpin_summary_len_sort.tsv

    # Join summary files
    LANG=en_EN join -1 2 -2 2 -t \$'\\t' \
            -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.3,2.4,2.5,2.6,2.7,2.8,2.9,2.10,2.11,2.12,2.13,2.14 \
            mature_summary_len_sort.tsv \
            hairpin_summary_len_sort.tsv  > summary_len_tmp.tsv
    LANG=en_EN sort -k 1 -t \$'\\t' summary_len_tmp.tsv -o summary_len.tsv
    sed -i '1iSpecies\tStress_event\t20_miRNAs\t21_miRNAs\t22_miRNAs\t23_miRNAs\t24_miRNAs\t25_miRNAs\t20_miRNAs_filt\t21_miRNAs_filt\t22_miRNAs_filt\t23_miRNAs_filt\t24_miRNAs_filt\t25_miRNAs_filt\t20_precursor\t21_precursor\t22_precursor\t23_precursor\t24_precursor\t25_precursor\t20_precursor_filt\t21_precursor_filt\t22_precursor_filt\t23_precursor_filt\t24_precursor_filt\t25_precursor_filt' summary_len.tsv
    """
}