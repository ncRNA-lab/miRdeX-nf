#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to calculate the absolute counts.
process COUNTS {

    tag "$meta.id"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("*.abs.tsv"), emit: abs
    
    script:
    def file_d = "$file".toString().replaceAll('.gz$', '')
    """
    # Decompress file 1
    if [[ ${file} == *.gz ]]; then
        gzip -d -f ${file}
    fi

    # Create the count table using bash. It is faster
    echo -e "seq\tcounts" > ${meta.id}.abs.tsv

    # Calculate absolute counts
    awk 'NR%4==2' "${file_d}" | sort | uniq -c | sort -nr | awk 'BEGIN{FS=" "; OFS="\\t"} {print \$2, \$1}' >> ${meta.id}.abs.tsv
    """
}

// Define the process to calculate the RPM
process RPM {

    tag "$meta.id"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.rpm.tsv"), emit: rpm
    
    script:
    """
    # Create the count table using bash. It is faster
    echo -e "seq\tRPM" > ${meta.id}_rpm.tsv

    # Perform the sum of all absolute counts (skipping the header)
    total=\$(awk 'BEGIN {FS="\\t"; total=0; getline} {total+=\$2} END {print total}' "${counts}")

    # Calculate the rpm (skipping the header)
    awk -v total="\${total}" 'BEGIN {FS="\\t"; OFS="\\t"; getline} {rpm = (\$2 / total) * 1000000; print \$1, rpm}' "${counts}" >> ${meta.id}.rpm.tsv
    """
}

// Define the process to calculate the RPM
process COUNTS_MATRIX {

    tag "$meta.id"
    conda "${modulesDir}/environment.yml"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${meta.id}_*.counts.tsv"), emit: matrix

    script:
    """
    04-Create_counts_matrix.py \
        --species ${meta.species_id} \
        --project ${meta.id} \
        --counts-tsv ${counts} \
        --metadata ${meta.metadata} \
        --valid-groups ${meta.valid_groups}
    """
}