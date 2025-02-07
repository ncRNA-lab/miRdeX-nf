process GROUP_BY_FAMILY {

    input:
    tuple val(meta), path(dea_file), path(annot_file)

    output:
    tuple val(meta), path("*fam_annot.tsv"), emit: fam_annot
    tuple val(meta), path("*boxplot.png"), emit: fam_boxplot
    path "*summary.tsv", emit: fam_sum

    script:
    """
    07-Group_miRNAs_by_family.r \
        --id ${meta.id} \
        --dea ${dea_file} \
        --annotation ${annot_file}
    """
}