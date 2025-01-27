#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to calculate the absolute counts.
process COUNTS {

    tag "$speciesID-$projectID"

    input:
    tuple val(speciesID), val(projectID), path(fastq_file)

    output:
    tuple val(speciesID), val(projectID), path("*_abs.tsv"), emit: abs
    
    script:
    """
    # Get the file name
    file_name=\$(basename "${fastq_file}" | cut -d'_' -f1)

    # Create the count table using bash. It is faster
    echo -e "seq\tcounts" > \${file_name}_abs.tsv

    # Calculate absolute counts
    awk 'NR%4==2' "${fastq_file}" | sort | uniq -c | sort -nr | awk 'BEGIN{FS=" "; OFS="\\t"} {print \$2, \$1}' >> \${file_name}_abs.tsv
    """
}

// Define the process to calculate the RPM
process RPM {

    tag "$speciesID-$projectID"

    input:
    tuple val(speciesID), val(projectID), path(abs_counts)

    output:
    tuple val(speciesID), val(projectID), path("*_rpm.tsv"), emit: rpm
    
    script:
    """
    # Get the file name
    file_name=\$(basename "${abs_counts}" | sed 's/_abs\\?\\.[^.]*\$//')

    # Create the count table using bash. It is faster
    echo -e "seq\tRPM" > \${file_name}_rpm.tsv

    # Perform the sum of all absolute counts (skipping the header)
    total=\$(awk 'BEGIN {FS="\\t"; total=0; getline} {total+=\$2} END {print total}' "${abs_counts}")

    # Calculate the rpm (skipping the header)
    awk -v total="\${total}" 'BEGIN {FS="\\t"; OFS="\\t"; getline} {rpm = (\$2 / total) * 1000000; print \$1, rpm}' "${abs_counts}" >> \${file_name}_rpm.tsv
    """
}

// Define the process to calculate the RPM
process COUNTS_MATRIX {

    tag "$speciesID-$projectID"
    conda "${modulesDir}/environment.yml"

    input:
    tuple val(speciesID), val(projectID), path(counts_files)
    path metadata_dir
    val filt_counts
    val filt_samples
    val rpm
    val project_table_bool
    val project_avg_bool
    val subproject_avg_bool

    output:
    tuple val(speciesID), val(projectID), path("${projectID}_counts.tsv"), emit: pabs, optional: true
    tuple val(speciesID), val(projectID), path("${projectID}_RPM.tsv"), emit: prpm, optional: true
    tuple val(speciesID), val(projectID), path("${projectID}_mean.tsv"), emit: pmean, optional: true
    tuple val(speciesID), val(projectID), path("${projectID}_*_counts.tsv"), emit: spabs, optional: true
    tuple val(speciesID), val(projectID), path("${projectID}_*_RPM.tsv"), emit: sprpm, optional: true
    tuple val(speciesID), val(projectID), path("${projectID}_*_mean.tsv"), emit: spmean, optional: true

    script:
    """
    # Create an empty string variable
    options=""

    # Check rpm option
    if [[ ${rpm} == "true" ]]
    then
        options+=" --rpm"
    fi
    # Check project-table option
    if [[ ${project_table_bool} == "true" ]]
    then
        options+=" --project-table"
    fi
    # Check avg-project option
    if [[ ${project_avg_bool} == "true" ]]
    then
        options+=" --avg-project"
    fi
    # Check avg-subproject option
    if [[ ${subproject_avg_bool} == "true" ]]
    then
        options+=" --avg-subproject"
    fi

    04-Create_counts_matrix.py \
        --species ${speciesID} \
        --project ${projectID} \
        --counts-tsv ${counts_files} \
        --metadata-dir ${metadata_dir} \
        --filt-counts ${filt_counts} \
        --filt-samples ${filt_samples} \
        \$options &
    """
}