#!/usr/bin/env nextflow

process DOWNLOADDB {

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

}
