process RPM {

    tag "$meta.id"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*.rpm.tsv"), emit: rpm
    
    script:
    """
    # Create the count table using bash. It is faster
    echo -e "seq\trpm" > ${meta.id}.rpm.tsv

    # Perform the sum of all absolute counts (skipping the header)
    total=\$(awk 'BEGIN {FS="\\t"; total=0; getline} {total+=\$2} END {print total}' "${counts}")

    # Calculate the rpm (skipping the header)
    awk -v total="\${total}" 'BEGIN {FS="\\t"; OFS="\\t"; getline} {rpm = (\$2 / total) * 1000000; print \$1, rpm}' "${counts}" >> ${meta.id}.rpm.tsv
    """
}
