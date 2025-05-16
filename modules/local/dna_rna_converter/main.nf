process DNA_RNA_CONVERTER {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd8170c44903910daa9e40d71124e4ccfb6e070bb6dc4c28e6928af9c6e3be2b/data' :
        'community.wave.seqera.io/library/bash:5.2.21--5bc877f5b6cf0654' }"
        
    input:
    tuple val(meta), path(file)
    val rna_to_dna

    output:
    tuple val(meta), path("*.{dna,rna}.{fq,fa,fastq,fasta}"), emit: file_conv

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # get the file extension
    filename=\$(basename "${file}")
    extension="\${filename##*.}"

    # Perform conversion based on format and rna_to_dna flag
    if [[ ${rna_to_dna} == true ]]; then
        if [[ "\$extension" == "fq" || "\$extension" == "fastq" ]]; then
            awk 'NR % 4 == 2 { gsub(/[Uu]/, "T") } { print }' "${file}" > "${prefix}.dna.\$extension"
        elif [[ "\$extension" == "fa" || "\$extension" == "fasta" ]]; then
            awk '/^>/ {print; next} {gsub(/[Uu]/, "T"); print}' "${file}" > "${prefix}.dna.\$extension"
        else
            echo "Error: Unsupported file format '\${extension}'. Must be .fa, .fasta, .fq, or .fastq."
            exit 1
        fi
    else
        if [[ "\$extension" == "fq" || "\$extension" == "fastq" ]]; then
            awk 'NR % 4 == 2 { gsub(/[Tt]/, "U") } { print }' "${file}" > "${prefix}.rna.\$extension"
        elif [[ "\$extension" == "fa" || "\$extension" == "fasta" ]]; then
            awk '/^>/ {print; next} {gsub(/[Tt]/, "U"); print}' "${file}" > "${prefix}.rna.\$extension"
        else
            echo "Error: Unsupported file format '\${extension}'. Must be .fa, .fasta, .fq, or .fastq."
            exit 1
        fi
    fi
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # get the file extension
    filename=\$(basename "${file}")
    extension="\${filename##*.}"

    if [[ ${rna_to_dna} == true ]]; then
        if [[ "\$extension" == "fq" || "\$extension" == "fastq" || "\$extension" == "fa" || "\$extension" == "fasta" ]]; then
            touch "${prefix}.dna.\$extension"
        else
            echo "Error: Unsupported file format '\${extension}'. Must be .fa, .fasta, .fq, or .fastq."
            exit 1
        fi
    else
        if [[ "\$extension" == "fq" || "\$extension" == "fastq" || "\$extension" == "fa" || "\$extension" == "fasta" ]]; then
            touch "${prefix}.rna.\$extension"
        else
            echo "Error: Unsupported file format '\${extension}'. Must be .fa, .fasta, .fq, or .fastq."
            exit 1
        fi
    fi

    """
}