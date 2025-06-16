process ANNOTATE_DEA_RESULTS {

    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://antoglz/mirnas_analysis/annotate_dea_results:latest' :
        'docker.io/antoglz/annotate_dea_results:latest' }"

    input:
    tuple val(meta), path(dea_file), path(annot_file)
    val classes

    output:
    tuple val(meta), path("*.all.tsv")      , emit: all
    tuple val(meta), path("*.unique.tsv")   , emit: unique
    tuple val(meta), path("*.boxplot*")     , emit: boxplot
    path "*.summary.tsv"                    , emit: fam_sum
    path  "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    06-Annotate_dea_results.r \
        --id ${prefix} \
        --dea ${dea_file} \
        --annotation ${annot_file} \
        --classes ${classes}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        argparse: \$(Rscript -e "suppressMessages(library(argparse)); cat(as.character(packageVersion('argparse')))") 
        ape: \$(Rscript -e "suppressMessages(library(ape)); cat(as.character(packageVersion('ape')))") 
        tidyverse: \$(Rscript -e "suppressMessages(library(tidyverse)); cat(as.character(packageVersion('tidyverse')))") 
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output files 
    touch ${prefix}.all.tsv
    touch ${prefix}.unique.tsv
    touch ${prefix}.boxplot_fam.png
    touch ${prefix}.boxplot_name.png
    touch ${prefix}.summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        argparse: \$(Rscript -e "suppressMessages(library(argparse)); cat(as.character(packageVersion('argparse')))") 
        ape: \$(Rscript -e "suppressMessages(library(ape)); cat(as.character(packageVersion('ape')))") 
        tidyverse: \$(Rscript -e "suppressMessages(library(tidyverse)); cat(as.character(packageVersion('tidyverse')))") 
    END_VERSIONS
    """
}
