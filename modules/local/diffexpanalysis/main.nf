process DIFFEXPANALYSIS {

    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://antoglz/mirnas_analysis/diffexp:latest' :
        'docker.io/antoglz/diffexp:latest' }"

    input:
    tuple val(meta), path(matrix), path(metadata), val(group_id)
    val alpha
    val min_counts
    val min_samples
    val lfc_threshold

    output:
    tuple val(meta), path("DEA/*_raw*.tsv")                             , emit: raw
    tuple val(meta), path("DEA/*_sig*.tsv")                             , emit: sig
    tuple val(meta), path("DEA/*.volcano.png")                          , emit: volcano
    tuple val(meta), path("Exploratory_analysis/01-PCA")                , emit: pca
    tuple val(meta), path("Exploratory_analysis/02-Mean_vs_variance")   , emit: var
    tuple val(meta), path("Exploratory_analysis/*.ea_summary.tsv")      , emit: easum
    tuple val(meta), path("DEA/*.dea_summary.tsv")                      , emit: deasum
    path  "versions.yml"                                                , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def lfc_th_arg = (lfc_threshold != 0) ? "--lfc_threshold ${lfc_threshold}" : ""
    """
    04-Diff_exp_analysis.r \
        --id ${prefix} \
        --group_id ${group_id} \
        --counts ${matrix} \
        --metadata ${metadata} \
        --alpha ${alpha} \
        --min_counts ${min_counts} \
        --min_samples ${min_samples} \
        ${lfc_th_arg}

    # Check if the output files are empty
    for file in DEA/*.dea_*.tsv; do
        [[ "\$file" == *EMPTY* ]] && continue
        (( \$(wc -l < "\$file") <= 1 )) && mv "\$file" "\${file%.tsv}_EMPTY.tsv"
    done
    
    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        argparse: \$(Rscript -e "suppressMessages(library(argparse)); cat(as.character(packageVersion('argparse')))") 
        dendextend: \$(Rscript -e "suppressMessages(library(dendextend)); cat(as.character(packageVersion('dendextend')))") 
        DESeq2: \$(Rscript -e "suppressMessages(library(DESeq2)); cat(as.character(packageVersion('DESeq2')))") 
        ff: \$(Rscript -e "suppressMessages(library(ff)); cat(as.character(packageVersion('ff')))") 
        htmltools: \$(Rscript -e "suppressMessages(library(htmltools)); cat(as.character(packageVersion('htmltools')))") 
        plotly: \$(Rscript -e "suppressMessages(library(plotly)); cat(as.character(packageVersion('plotly')))") 
        RColorBrewer: \$(Rscript -e "suppressMessages(library(RColorBrewer)); cat(as.character(packageVersion('RColorBrewer')))") 
        readr: \$(Rscript -e "suppressMessages(library(readr)); cat(as.character(packageVersion('readr')))") 
        SARTools: \$(Rscript -e "suppressMessages(library(SARTools)); cat(as.character(packageVersion('SARTools')))") 
        stringr: \$(Rscript -e "suppressMessages(library(stringr)); cat(as.character(packageVersion('stringr')))") 
        tibble: \$(Rscript -e "suppressMessages(library(tibble)); cat(as.character(packageVersion('tibble')))") 
        tidyverse: \$(Rscript -e "suppressMessages(library(tidyverse)); cat(as.character(packageVersion('tidyverse')))") 
        ggthemes: \$(Rscript -e "suppressMessages(library(ggthemes)); cat(as.character(packageVersion('ggthemes')))") 
        DEGreport: \$(Rscript -e "suppressMessages(library(DEGreport)); cat(as.character(packageVersion('DEGreport')))") 
        paletteer: \$(Rscript -e "suppressMessages(library(paletteer)); cat(as.character(packageVersion('paletteer')))") 
        gtools: \$(Rscript -e "suppressMessages(library(gtools)); cat(as.character(packageVersion('gtools')))") 
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Create output directories
    mkdir -p DEA
    mkdir -p Exploratory_analysis/01-PCA
    mkdir -p Exploratory_analysis/02-Mean_vs_variance

    # Create raw and sig tables
    touch DEA/${prefix}_1.dea_raw.tsv
    touch DEA/${prefix}_1.dea_sig.tsv
    touch DEA/${prefix}_1.volcano.png

    # Create EA summary file
    ea_summary=Exploratory_analysis/${meta.id}.ea_summary.tsv
    echo -e "Group_id\tGroup\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tP-value(MWW)" > \$ea_summary
    echo -e "1\t${meta.id}\t37.93\t22.45\t15.52\t13.8\t10.3\t0\t0.00759240759240759" >> \$ea_summary
    
    # Get the sample names from the matrix file (column names excluding 'seq')
    sample_names=\$(head -n 1 ${matrix} | tr '\\t' '\\n' | grep -v '^seq\$' | tr '\\n' ',' | sed 's/,\$//')

    # Create DEA summary file
    dea_summary=DEA/${meta.id}.dea_summary.tsv
    echo -e "Group\tTest\tPadj<alpha\tTotal\tCoefficient\tContrast\tContrast_coefficient\tSamples" > \$dea_summary
    echo -e "${meta.id}_1\tWald\t4503\t91672\tTime_24h_vs_0h\tNo contrast\tNo contrast\t\$sample_names" >> \$dea_summary

    # Create versions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        argparse: \$(Rscript -e "suppressMessages(library(argparse)); cat(as.character(packageVersion('argparse')))") 
        dendextend: \$(Rscript -e "suppressMessages(library(dendextend)); cat(as.character(packageVersion('dendextend')))") 
        DESeq2: \$(Rscript -e "suppressMessages(library(DESeq2)); cat(as.character(packageVersion('DESeq2')))") 
        ff: \$(Rscript -e "suppressMessages(library(ff)); cat(as.character(packageVersion('ff')))") 
        htmltools: \$(Rscript -e "suppressMessages(library(htmltools)); cat(as.character(packageVersion('htmltools')))") 
        plotly: \$(Rscript -e "suppressMessages(library(plotly)); cat(as.character(packageVersion('plotly')))") 
        RColorBrewer: \$(Rscript -e "suppressMessages(library(RColorBrewer)); cat(as.character(packageVersion('RColorBrewer')))") 
        readr: \$(Rscript -e "suppressMessages(library(readr)); cat(as.character(packageVersion('readr')))") 
        SARTools: \$(Rscript -e "suppressMessages(library(SARTools)); cat(as.character(packageVersion('SARTools')))") 
        stringr: \$(Rscript -e "suppressMessages(library(stringr)); cat(as.character(packageVersion('stringr')))") 
        tibble: \$(Rscript -e "suppressMessages(library(tibble)); cat(as.character(packageVersion('tibble')))") 
        tidyverse: \$(Rscript -e "suppressMessages(library(tidyverse)); cat(as.character(packageVersion('tidyverse')))") 
        ggthemes: \$(Rscript -e "suppressMessages(library(ggthemes)); cat(as.character(packageVersion('ggthemes')))") 
        DEGreport: \$(Rscript -e "suppressMessages(library(DEGreport)); cat(as.character(packageVersion('DEGreport')))") 
        paletteer: \$(Rscript -e "suppressMessages(library(paletteer)); cat(as.character(packageVersion('paletteer')))") 
        gtools: \$(Rscript -e "suppressMessages(library(gtools)); cat(as.character(packageVersion('gtools')))") 
    END_VERSIONS
    """


}