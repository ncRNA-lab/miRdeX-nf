process BOWTIE_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8c0819a9b1f520c49c933e667ae50de2a0730ece4c8b9efe79ac5e403963a9f/data' :
        'community​.wave​.seqera​.io/library/bowtie_samtools:e1a14e1ce4e0170d' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val (save_aligned) // quitar parentesis
    val (save_unaligned) // quitar parentesis

    output:
    tuple val(meta), path('*.bam')     , emit: bam
    tuple val(meta), path('*.out')     , emit: log
    tuple val(meta), path("*.${meta.type}.aligned.fastq.gz") , emit: aligned, optional : true
    tuple val(meta), path("*.${meta.type}.unaligned.fastq.gz") , emit: unaligned, optional : true
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligned = save_aligned ? "--al ${prefix}.${meta.type}.aligned.fastq" : ''
    def unaligned = save_unaligned ? "--un ${prefix}.${meta.type}.unaligned.fastq" : ''
    def endedness = meta.single_end ? "$reads".toString().replaceAll('.gz$', '') : "-1 ${reads[0].toString().replaceAll('.gz$', '')} -2  ${reads[1].toString().replaceAll('.gz$', '')}"
    """
    # Decompress query files
    if [[ ${meta.single_end} == true ]]; then
        gzip -d -f ${reads}
    else
        # Decompress file 1
        if [[ ${reads[1]} == *.gz ]]; then
            gzip -d -f ${reads[1]}
        fi

        # Decompress file 2
        if [[ ${reads[2]} == *.gz ]]; then
            gzip -d -f ${reads[2]}
        fi
    fi

    INDEX=\$(find -L ./ -name "*.3.ebwt" | sed 's/\\.3.ebwt\$//')
    bowtie \\
        --threads $task.cpus \\
        --sam \\
        -x \$INDEX \\
        -q \\
        $aligned \\
        $unaligned \\
        $args \\
        $endedness \\
        2> >(tee ${prefix}.out >&2) \\
        | samtools view $args2 -@ $task.cpus -bS -o ${prefix}.bam -

    # Compress aligned sequences file
    if [ -f ${prefix}.${meta.type}.aligned.fastq ]; then
        gzip ${prefix}.${meta.type}.aligned.fastq
    fi
    if [ -f ${prefix}.${meta.type}.aligned_1.fastq ]; then
        gzip ${prefix}.${meta.type}.aligned_1.fastq
        gzip ${prefix}.${meta.type}.aligned_2.fastq
    fi

    # Compress unaligned sequences file
    if [ -f ${prefix}.${meta.type}.unaligned.fastq ]; then
        gzip ${prefix}.${meta.type}.unaligned.fastq
    fi
    if [ -f ${prefix}.${meta.type}.unaligned_1.fastq ]; then
        gzip ${prefix}.${meta.type}.unaligned_1.fastq
        gzip ${prefix}.${meta.type}.unaligned_2.fastq
    fi



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def aligned = save_aligned ?
                    meta.single_end ? "echo '' | gzip > ${prefix}.aligned.fastq.gz" :
                        "echo '' | gzip > ${prefix}.aligned_1.fastq.gz; echo '' | gzip > ${prefix}.aligned_2.fastq.gz"
                    : ''
    def unaligned = save_unaligned ?
                    meta.single_end ? "echo '' | gzip > ${prefix}.unaligned.fastq.gz" :
                        "echo '' | gzip > ${prefix}.unaligned_1.fastq.gz; echo '' | gzip > ${prefix}.unaligned_2.fastq.gz"
                    : ''
    """
    touch ${prefix}.bam
    touch ${prefix}.out
    $aligned
    $unaligned

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """


}
