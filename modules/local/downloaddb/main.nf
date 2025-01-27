#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

// Define the process to execute Prefetch.sh.
process DOWNLOADDB {

    output:
    path "pmiren/*/*mature.fa"  , emit: pmiren_mature
    path "pmiren/*/*hairpin.fa" , emit: pmiren_hairpin
    path "mirbase/mature.fa" , emit: mirbase_mature
    path "mirbase/hairpin.fa" , emit: mirbase_hairpin
    path "srnaanno/*.miRNA.gff3" , emit: srnaanno_all

    script:
    """
    # Download miRBase FASTA files
    wget --no-check-certificate -P mirbase https://mirbase.org/download/mature.fa
    wget --no-check-certificate -P mirbase https://mirbase.org/download/hairpin.fa

    # Download PmiREN FASTA files
    wget --no-check-certificate -r --no-parent -nH --cut-dirs=1 -P pmiren --accept "*.fa" https://www.pmiren.com/ftp-download

    # Download sRNAanno GFF3 files
    wget --no-check-certificate http://121.37.229.61:84/sRNAannoDATA/miRNAGFF3/miRNA.gff3.tar.gz
    mkdir -p srnaanno && tar -xzvf miRNA.gff3.tar.gz -C srnaanno
    rm miRNA.gff3.tar.gz
    """
}

//    # Download miRBase FASTA files
//    wget --no-check-certificate -P mirbase https://mirbase.org/download/mature.fa
//    wget --no-check-certificate -P mirbase https://mirbase.org/download/hairpin.fa

//# Download PmiREN FASTA files
//    # wget --no-check-certificate -r --no-parent -nH --cut-dirs=1 -P pmiren --accept "*.fa" https://www.pmiren.com/ftp-download