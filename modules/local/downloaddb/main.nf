#!/usr/bin/env nextflow

process DOWNLOADDB {

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data' :
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"
        
    input:
    val mirbase
    val srnaanno
    val pmiren

    output:
    path "mirbase/*mature.fa" , emit: mirbase_mature
    path "mirbase/*hairpin.fa" , emit: mirbase_hairpin
    path "mirbase/*mirbase_species.txt" , emit: mirbase_species
    path "srnaanno/*.miRNA.gff3" , emit: srnaanno_all
    path "pmiren/*/*mature.fa"  , emit: pmiren_mature
    path "pmiren/*/*hairpin.fa" , emit: pmiren_hairpin


    script:
    def script_lines = []

    // Download miRBase
    if (mirbase) {
        script_lines.add("""
        # Download miRBase FASTA files
        wget --no-check-certificate -P mirbase https://mirbase.org/download/mature.fa
        wget --no-check-certificate -P mirbase https://mirbase.org/download/hairpin.fa
        wget --no-check-certificate -P mirbase https://mirbase.org/download/CURRENT/database_files/mirna_species.txt
        awk '{gsub("<p>", ""); gsub("</p>", ""); gsub("<br>", "\\n"); print}' mirbase/mirna_species.txt > mirbase/mirbase_species.txt
        """)
    } else {
        script_lines.add("""
        mkdir -p mirbase
        touch mirbase/EMPTY_mirbase.mature.fa
        touch mirbase/EMPTY_mirbase.hairpin.fa
        touch mirbase/EMPTY_mirbase_species.txt
        """)
    }

    // Download sRNAanno
    if (srnaanno) {
        script_lines.add("""
        # Download sRNAanno GFF3 files
        wget --no-check-certificate http://121.37.229.61:84/sRNAannoDATA/miRNAGFF3/miRNA.gff3.tar.gz
        mkdir -p srnaanno && tar -xzvf miRNA.gff3.tar.gz -C srnaanno
        rm miRNA.gff3.tar.gz
        """)
    } else {
        script_lines.add("""
        mkdir -p srnaanno
        touch srnaanno/EMPTY_srnaanno.miRNA.gff3
        """)
    }

    // Download PmiREN
    if (pmiren) {
        script_lines.add("""
        # Download PmiREN FASTA files
        wget --no-check-certificate -r --no-parent -nH --cut-dirs=1 -P pmiren --accept "*.fa" https://www.pmiren.com/ftp-download
        """)
    } else {
        script_lines.add("""
        mkdir -p pmiren/empty
        touch pmiren/empty/EMPTY_pmiren.mature.fa
        touch pmiren/empty/EMPTY_pmiren.hairpin.fa
        """)
    }

    // Execute only the required lines
    script_lines.join("\n")

    stub:
    """
    # Create directories for mirbase, srnaanno, and pmiren
    mkdir -p mirbase srnaanno pmiren/Arabidopsis_thaliana_Ath

    # Simulate mirbase files (empty files if mirbase is not selected)
    if [ "${mirbase}" == "true" ]; then
        # Create db files
        touch mirbase/mature.fa
        touch mirbase/hairpin.fa
        touch mirbase/mirbase_species.txt

    else
        touch mirbase/EMPTY_mirbase.mature.fa
        touch mirbase/EMPTY_mirbase.hairpin.fa
        touch mirbase/EMPTY_mirbase_species.txt
    fi

    # Simulate srnaanno files (empty files if srnaanno is not selected)
    if [ "${srnaanno}" == "true" ]; then
        # Create db files
        touch srnaanno/Arabidopsis_thaliana.miRNA.gff3
        touch srnaanno/Brassica_napus.miRNA.gff3
    else
        touch srnaanno/EMPTY_srnaanno.miRNA.gff3
    fi

    # Simulate pmiren files (empty files if pmiren is not selected)
    if [ "${pmiren}" == "true" ]; then
        # Create db files
        touch pmiren/Arabidopsis_thaliana_Ath/Arabidopsis_thaliana_mature.fa
        touch pmiren/Arabidopsis_thaliana_Ath/Arabidopsis_thaliana_hairpin.fa
    else
        echo "Simulating empty files for pmiren..."
        touch pmiren/empty/EMPTY_pmiren.mature.fa
        touch pmiren/empty/EMPTY_pmiren.hairpin.fa
    fi
    """
}
