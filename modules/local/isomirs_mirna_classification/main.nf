process ISOMIRS_MIRNA_CLASSIFICATION {

    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/ee1cb5ff22056795ddb2bb8ddc2f755a32f7c463c12665f5f246919ce595307e/data' :
        'community.wave.seqera.io/library/python_pip_biopython_pandas:5b18ba7531011a4b' }"
        
    input:
    tuple val(meta), path(blast_tsv)
    val substitutions
    val five_add
    val three_add
    val ends_modification

    output:
    tuple val(meta), path("*gff3")         , emit: gff3
    tuple val(meta), path("*.summary.tsv") , emit: sum
    path "versions.yml"                    , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Execute isomir classification script
    05-IsomiRs_classification.py \
        --id ${prefix} \
        --database ${meta.database} \
        --input ${blast_tsv} \
        --substitutions ${substitutions} \
        --five_add ${five_add} \
        --three_add ${three_add} \
        --ends_modification ${ends_modification}
    
    # Create vesions file
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gff3

    echo -e "num_isomirs\\tref_miRNA\\tiso_5p\\tiso_3p\\tiso_add3p\\tiso_add5p\\tiso_snv_seed" \\
    "\\tiso_snv_central_offset\\tiso_snv_central\\tiso_snv_central_supp\\tiso_snv\\tmixed" \\
    "\\tmixed_shift\\tundefined\\n2119\\t207\\t245\\t346\\t60\\t0\\t1\\t0\\t0\\t0\\t0\\t274\\t966" \\
    "\\t20" > ${prefix}.summary.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c "import platform; print(platform.python_version())")
        pandas: \$(python3 -c "import pandas as pd; print(pd.__version__)")
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

}