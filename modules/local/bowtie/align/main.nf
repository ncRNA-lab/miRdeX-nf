process BOWTIE_ALIGN {
    
    tag "$meta.id"
    
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c8/c8c0819a9b1f520c49c933e667ae50de2a0730ece4c8b9efe79ac5e403963a9f/data' :
        'community.wave.seqera.io/library/bowtie_samtools:e1a14e1ce4e0170d' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)
    val save_aligned
    val save_unaligned

    output:
    tuple val(meta), path('*.bam') , emit: bam
    tuple val(meta), path('*.out') , emit: log
    tuple val(meta), path("*.aligned*.{fa,fasta,fq,fastq}.gz")   , emit: aligned, optional : true
    tuple val(meta), path("*.unaligned*.{fa,fasta,fq,fastq}.gz") , emit: unaligned, optional : true
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ref = task.ext.ref ?: ''
    def ext = meta.single_end ? reads.name.replaceAll(/\.gz$/, '').tokenize('.')[-1] : reads[0].name.replaceAll(/\.gz$/, '').tokenize('.')[-1]
    def aligned = save_aligned ? "--al ${prefix}.${ref}.aligned.${ext}" : ''
    def unaligned = save_unaligned ? "--un ${prefix}.${ref}.unaligned.${ext}" : ''
    def endedness = meta.single_end ? "$reads".toString().replaceAll('.gz$', '') : "-1 ${reads[0].toString().replaceAll('.gz$', '')} -2  ${reads[1].toString().replaceAll('.gz$', '')}"
    """
    # Decompress query files
    if [[ ${meta.single_end} == true ]]; then
        # Decompress file if it is required
        if [[ ${reads} == *.gz ]]; then
            gzip -d -f ${reads}
        fi
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
    if [ -f ${prefix}.${ref}.aligned.${ext} ]; then
        gzip ${prefix}.${ref}.aligned.${ext}
    fi
    if [ -f ${prefix}.${ref}.aligned_1.${ext} ]; then
        gzip ${prefix}.${ref}.aligned_1.${ext}
        gzip ${prefix}.${ref}.aligned_2.${ext}
    fi

    # Compress unaligned sequences file
    if [ -f ${prefix}.${ref}.unaligned.${ext} ]; then
        gzip ${prefix}.${ref}.unaligned.${ext}
    fi
    if [ -f ${prefix}.${ref}.unaligned_1.${ext} ]; then
        gzip ${prefix}.${ref}.unaligned_1.${ext}
        gzip ${prefix}.${ref}.unaligned_2.${ext}
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version) | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ext = meta.single_end ? reads.name.replaceAll(/\.gz$/, '').tokenize('.')[-1] : reads[0].name.replaceAll(/\.gz$/, '').tokenize('.')[-1]
    def aligned = save_aligned ?
                    meta.single_end ? "touch ${prefix}.aligned.${ext}.gz" :
                        "touch ${prefix}.aligned_1.${ext}.gz; touch ${prefix}.aligned_2.${ext}.gz"
                    : ''
    def unaligned = save_unaligned ?
                    meta.single_end ? "touch ${prefix}.unaligned.${ext}.gz" :
                        "touch ${prefix}.unaligned_1.${ext}.gz; touch ${prefix}.unaligned_2.${ext}.gz"
                    : ''
    """
    touch ${prefix}.bam
    ${aligned}
    ${unaligned}

    # Create the log file
    echo "# reads processed: 921451" > ${prefix}.out
    echo "# reads with at least one alignment: 640255 (69.48%)" >> ${prefix}.out
    echo "# reads that failed to align: 281196 (30.52%)" >> ${prefix}.out
    echo "# reads with alignments suppressed due to -m: 2459 (0.27%)" >> ${prefix}.out
    echo "# Reported 637796 alignments" >> ${prefix}.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bowtie: \$(echo \$(bowtie --version 2>&1) | sed 's/^.*bowtie-align-s version //; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """


}