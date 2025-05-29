process ADD_SEQUENCES_BLAST {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f8/f8dc20302635f868211541a2a2537104da9c75ce7efd85e81cdb93947652f265/data' :
        'community.wave.seqera.io/library/seqkit_gzip:17e1ce4d127a4305' }"
        
    input:
    tuple val(meta), path(blast_tsv), path(query_fasta), path(subject_fasta)

    output:
    tuple val(meta), path("*.blast.tsv") , emit: bseqs
    path "versions.yml"                  , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def q_is_compressed = query_fasta.getExtension() == "gz" ? true : false
    def q_fasta_name = q_is_compressed ? query_fasta.getBaseName() : query_fasta
    def s_is_compressed = subject_fasta.getExtension() == "gz" ? true : false
    def s_fasta_name = s_is_compressed ? subject_fasta.getBaseName() : subject_fasta
    """
    #!/bin/bash
    if [ "${q_is_compressed}" == "true" ]; then
        gzip -c -d ${query_fasta} > ${q_fasta_name}
    fi
    if [ "${s_is_compressed}" == "true" ]; then
        gzip -c -d ${subject_fasta} > ${s_fasta_name}
    fi

    # Get the sequence ids
    cut -f1 "${blast_tsv}" | sort | uniq > query_ids.txt
    cut -f2 "${blast_tsv}" | sort | uniq > subject_ids.txt

    # Get the sequences
    seqkit grep -f query_ids.txt "${query_fasta}" | seqkit fx2tab | awk -F'\t' '{split(\$1, a, " "); print a[1] "\t" \$2}' > query_seqs.tsv
    seqkit grep -f subject_ids.txt "${subject_fasta}"  | seqkit fx2tab | awk -F'\t' '{split(\$1, a, " "); print a[1] "\t" \$2}' > subject_seqs.tsv

    # Merge the data
    awk 'FNR==NR {seqs[\$1]=\$2; next} {print \$0 "\t" seqs[\$1]}' query_seqs.tsv "${blast_tsv}" > tmp_with_q.tsv
    awk 'FNR==NR {seqs[\$1]=\$2; next} {print \$0 "\t" seqs[\$2]}' subject_seqs.tsv tmp_with_q.tsv > "${prefix}.blast.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/^.*seqkit //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.blast.tsv"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$(seqkit version 2>&1 | sed 's/^.*seqkit //; s/ .*\$//')
    END_VERSIONS
    """
}