#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   isomiR_Classification.py
#
#   This program identifies potential isomiRs by performing a BLASTN search 
#   between a query file (containing sequences of interest) and a reference 
#   file (containing canonical miRNA sequences). After obtaining the BLASTN 
#   results, the program classifies the identified isomiRs, annotating them 
#   based on the type of modification (e.g., sequence shifts, truncations, 
#   extensions, or substitutions) relative to the canonical miRNA.
#
#   The following modifications can occur in isomiRs relative to the canonical
#   miRNA:
#
#   1. 5' Overhang: Sequence present at the 5' end of the isomiR that is absent 
#      in the canonical miRNA.
#
#   2. 5' Deletion: Sequence present in the canonical miRNA but absent at the
#      5' end of the isomiR.
#
#   3. 5' Variation: Sequence variation between the 5' ends of the isomiR and 
#      the canonical miRNA, including mismatches or insertions.
#
#   4. 3' Overhang: Sequence present at the 3' end of the isomiR that is absent 
#      in the canonical miRNA.
#
#   5. 3' Deletion: Sequence present in the canonical miRNA but absent at the 3' 
#      end of the isomiR.
#
#   6. 3' Variation: Sequence variation between the 3' ends of the isomiR and 
#      the canonical miRNA, including mismatches or insertions.
#
#   7. Polymorphic Changes: Sequence differences between the isomiR and the
#      canonical miRNA that result in substitutions or point mutations.
#
#   The following steps are performed:
#
#   1. **Index the reference file** using BLAST's `makeblastdb`.
#
#   The reference file (containing canonical miRNA sequences) is indexed using 
#   the `makeblastdb` tool, which prepares the reference file for fast sequence 
#   comparison using BLASTN.
#
#   2. Run BLASTN to identify potential isomiRs.
#
#   The query file containing the sequences of interest is aligned against the 
#   indexed reference file using BLASTN. The result of this alignment is used 
#   to identify potential isomiRs, with each query sequence being compared 
#   to a reference miRNA sequence.
#
#   3. Classify the identified isomiRs.
#
#   Once potential isomiRs are identified, the program classifies each sequence 
#   based on the type of modification it has undergone in comparison to the 
#   canonical miRNA sequence. This classification provides important information 
#   about the nature of the sequence modifications, such as truncation, 
#   extension, or sequence shift.
#
#   4. Annotate the isomiRs with classification results.
#
#   The classification results, which describe the modifications each isomiR 
#   has undergone, are combined with the BLASTN output. This annotated table 
#   provides the final information about each isomiR, including the alignment 
#   data from BLAST and the classification of the isomiRs based on their modifications.
#
#   5. Return the annotated isomiRs.
#
#   The final output is a table containing the annotated isomiRs, which includes 
#   both the BLASTN results (e.g., alignment identity, mismatches, gaps, etc.) 
#   and the classification of the isomiRs based on their modifications relative 
#   to the canonical miRNA.
#
#   Authors: Antonio González Sánchez
#   Date: March 17, 2025
#   Version: 1.0
#
#******************************************************************************

## IMPORT MODULES
import argparse
import subprocess
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import pandas as pd
import os
import itertools
import sys
import re

## CLASSES
class IsomiRs:
    """
    A class to analyze and classify isomiR sequences, which are sequence variants of miRNAs.
    
    This class provides methods to identify sequence modifications such as 5' and 3' overhangs,
    deletions, variations, and polymorphic changes. It also assigns a classification identifier
    based on these modifications.

    Attributes:
        seq (str): The sequence of the isomiR.
        canonical_name (str): The canonical miRNA name(s) associated with this isomiR.
        canonical_seq (str): The reference sequence of the miRNA.
        hairpin_list (list): List of hairpin precursors associated with this isomiR.
        aligner (PairwiseAligner): The alignment tool used for sequence comparison.
        one_alignment: The best alignment found between the isomiR and the canonical sequence.
        seq_start (int): Start position of the aligned region in the isomiR sequence.
        can_start (int): Start position of the aligned region in the canonical sequence.

    Methods:
        swalignment():
            Performs sequence alignment between the isomiR and the canonical sequence.
        five_prime_overhang():
            Identifies a 5' overhang in the isomiR sequence, if present.
        five_prime_deletion():
            Identifies a 5' deletion in the isomiR sequence, if present.
        five_prime_variation():
            Identifies a 5' variation between the isomiR and the canonical sequence.
        _remove_5prime_overhang():
            Removes the 5' overhang from the isomiR sequence for further analysis.
        _get_aligned_regions():
            Extracts the aligned regions from the sequence alignment.
        _get_clean_aligned_regions():
            Retrieves the aligned sequences without gaps.
        three_prime_overhang():
            Identifies a 3' overhang in the isomiR sequence, if present.
        three_prime_deletion():
            Identifies a 3' deletion in the isomiR sequence, if present.
        three_prime_variation():
            Identifies a 3' variation between the isomiR and the canonical sequence.
        polymorphic_changes():
            Detects polymorphic nucleotide changes between the isomiR and the canonical sequence.
        is_template():
            Checks if the isomiR sequence is present in any of the provided hairpin sequences.
        _get_isomir_classification_suffix():
            Generates an isomiR classification identifier based on detected modifications.
        isomir_name():
            Generates a standardized isomiR name based on associated miRNA variants.
        isomir_classification_results():
            Returns a dictionary containing classification details of the isomiR.
    """

    def __init__(self, seq, canonical_name, canonical_seq, hairpin_list):
        """
        Initializes the IsomiRs class for the analysis of IsomiR sequences.
        
        :param seq: Sequence to analyze.
        :param canonical_seq: Canonical reference sequence.
        :param hairpin_list: List of hairpin structures.
        """
        self.seq = seq
        self.canonical_name = canonical_name
        self.canonical_seq = canonical_seq
        self.hairpin_list = hairpin_list

        # Initialize alignment-related attributes
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'local'
        self.aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
        self.aligner.open_gap_score = -8
        self.aligner.extend_gap_score = -8
        self.one_alignment = None
        self.seq_start = None
        self.can_start = None
        self.seq_wo_five_prime_overhang = None
        self.can_wo_five_prime_overhang = None
        self.seq_aligned = None
        self.can_aligned = None
        self.seq_aligned_clean = None
        self.can_aligned_clean = None

        # Execute the alignment
        self.swalignment()

    def swalignment(self):
        """Performs the alignment between the sequences."""
        alignments = self.aligner.align(self.seq, self.canonical_seq)
        self.one_alignment = max(alignments, key=lambda aln: aln.score)
        self.seq_start = self.one_alignment.aligned[0][0][0]
        self.can_start = self.one_alignment.aligned[1][0][0]

    def five_prime_overhang(self):
        """Gets the 5' overhang sequence if it exists."""
        # Check if a 5' overhang exists.
        if self.seq_start > 0 and self.can_start == 0:
            five_overhang_seq = self.seq[:self.seq_start]
            return five_overhang_seq
        return None

    def five_prime_deletion(self):
        """Gets the 5' deleted sequence if it exists."""
        # Check if a 5' deletion exists.
        if self.seq_start == 0 and self.can_start > 0:
            five_deletion_seq = self.canonical_seq[:self.can_start]
            return five_deletion_seq
        return None

    def five_prime_variation(self):
        """Gets the 5' variation between the sequence and the canonical sequence if it exists."""
        # Check if a 5' variation exists.
        if self.seq_start > 0 and self.can_start > 0:
            five_variation_seq = self.seq[:self.seq_start]
            five_variation_can = self.canonical_seq[:self.can_start]
            return (five_variation_seq, five_variation_can)
        return None

    def three_prime_overhang(self):
        """Gets the 3' overhang sequence if it exists."""
        # Execute the requiered methods
        seq_wo_five_prime_overhang, can_wo_five_prime_overhang = self._remove_5prime_overhang()
        seq_aligned_clean, _ = self._get_clean_aligned_regions()
        # Check if a 3' overghang exists.
        if len(seq_wo_five_prime_overhang) > len(seq_aligned_clean) and len(seq_aligned_clean) == len(can_wo_five_prime_overhang):
            three_overhang_seq = seq_wo_five_prime_overhang[len(seq_aligned_clean):]
            return three_overhang_seq
        return None

    def three_prime_deletion(self):
        """Gets the 3' deleted sequence if it exists."""
        # Execute the required methods
        seq_wo_five_prime_overhang, can_wo_five_prime_overhang = self._remove_5prime_overhang()
        seq_aligned_clean, can_aligned_clean = self._get_clean_aligned_regions()
        # Check if a 3' deletion exists.
        if len(seq_wo_five_prime_overhang) == len(seq_aligned_clean) and len(can_aligned_clean) < len(can_wo_five_prime_overhang):
            three_deletion_seq = can_wo_five_prime_overhang[len(can_aligned_clean):]
            return three_deletion_seq
        return None

    def three_prime_variation(self):
        """Gets the 3' variation between the sequence and the canonical sequence if it exists."""
        # Exectue the required methods
        seq_wo_five_prime_overhang, can_wo_five_prime_overhang = self._remove_5prime_overhang()
        seq_aligned_clean, can_aligned_clean = self._get_clean_aligned_regions()
        # Check if a 3' variation exists.
        if len(seq_wo_five_prime_overhang) > len(seq_aligned_clean) and len(can_aligned_clean) < len(can_wo_five_prime_overhang):
            three_variation_seq = seq_wo_five_prime_overhang[len(seq_aligned_clean):]
            three_variation_can = can_wo_five_prime_overhang[len(can_aligned_clean):]
            three_variation = (three_variation_seq, three_variation_can)
            return three_variation
        return None
    
    def polymorphic_changes(self):
        """Gets the polymorphic changes between the sequence and the canonical sequence."""
        # Get 5' modifications and ensure they are not None
        five_overhang_seq = '' if self.five_prime_overhang() is None else self.five_prime_overhang()
        five_deletion_seq = '' if self.five_prime_deletion() is None else self.five_prime_deletion()
        # Ensure aligned sequences are obtained
        seq_aligned, can_aligned = self._get_aligned_regions()
        polymorphic_changes = []
        # Iterate through the aligned sequences
        for i in range(len(seq_aligned)):
            seq_nuc = seq_aligned[i]
            can_nuc = can_aligned[i]
            # If there's a polymorphic difference
            if seq_nuc != can_nuc:
                real_query_nuc_position = i + len(five_overhang_seq)
                real_ref_nuc_position = i + len(five_deletion_seq)
                polymorphic_changes.append((seq_nuc, real_query_nuc_position + 1,
                                            can_nuc, real_ref_nuc_position + 1))        
        # Check if a change in the sequence exits
        if len(polymorphic_changes) > 0:
            return polymorphic_changes
        return None
    
    def is_template(self):
        """Checks if the sequence is present in any hairpin from the hairpin list."""
        return any(self.seq in hairpin for hairpin in self.hairpin_list)
    
    def isomir_name(self):
        """
        Generates a unique isomiR name based on the canonical miRNA variants it is associated with.

        This method processes the canonical miRNA name, which may contain multiple variants separated by "|",
        and extracts relevant information to create a unified isomiR name. Since different miRNA variants from
        distinct precursors can share the same sequence, this method ensures a comprehensive representation.

        Example of a canonical miRNA name:
            "ath-miR156a-5p|ath-miR156b-5p|ath-miR156c-5p|ath-miR156d-5p|ath-miR156e|ath-miR156f-5p"
        
        This method processes such names and generates a standardized identifier for the isomiR.

        Example output for the above canonical miRNA name:
            "miR156abcdef-5p-<classification_suffix>"
        
        Where:
        - "miR156abcdef-5p" represents the unified name including all miRNA variants.
        - "<classification_suffix>" is the identifier generated by `_get_isomir_classification_suffix()`,
          which represents the sequence modifications. If no modifications are present, the suffix will be "0:0:0:0:0:0:0".

        Returns:
            str: A unique isomiR name including all associated miRNA variants and the classification suffix.
        """
        
        # Get the classification suffix
        class_suffix = self._get_isomir_classification_suffix()

        # Split the string by the delimiter "|"
        miRNA_vars = self.canonical_name.split('|')
        miRNAs_vars_dic = {}
        # Preprocess to extract the relevant information
        for miRNA in miRNA_vars:

            # Get the miRNA parts (ath-miR398a-5p -> ['ath', 'miR398a', '5p'])
            parts = miRNA.split('-')

            # Get the simple miRNA name (['ath', 'miR398a', '5p'] -> miR398a)
            miRNA_simple = next((item for item in parts if 'mir' in item.lower()), None)
            
            # If no valid miRNA name is found, skip to the next one
            if not miRNA_simple:
                continue
            
            # Get the miRNA name (miR398a -> miR398)
            miRNA_name = re.search(r'(miR\d+(\.\d+)?)', miRNA_simple, re.IGNORECASE).group(1)
            
            # Extract the variant (letter) (miR398a -> a)
            var_match = re.search(r'\d([a-zA-Z])', miRNA_simple)
            var = var_match.group(1) if var_match else ''
            
            # Extract the strand (5p or 3p)
            #strand = next((part for part in parts if part.endswith("p")), '')
            strand = next((part for part in parts if re.search(r"(5p|3p)$", part)), '')

            # Store in the dictionary
            if miRNA_name not in miRNAs_vars_dic:
                miRNAs_vars_dic[miRNA_name] = {'vars': [], 'strands': []}
            
            miRNAs_vars_dic[miRNA_name]['vars'].append(var)
            miRNAs_vars_dic[miRNA_name]['strands'].append(strand)

        # Create miRNA names for isomiR classification
        miRNA_class_names_list = []

        for miRNA, elements in miRNAs_vars_dic.items():
            # Concatenate all variants
            vars = "".join(elements['vars'])

            # Get the strand, if exists
            strand = next((item for item in elements['strands'] if item), None)

            # Format the name
            if strand:
                miRNA_class_names_list.append(f'{miRNA}{vars}-{strand}')
            else:
                miRNA_class_names_list.append(f'{miRNA}{vars}')
        
        # Build the isomiR name
        isomir_name = f'{"|".join(miRNA_class_names_list)}-{class_suffix}'

        return isomir_name

    def isomir_classification_results(self):
        """
        Returns a dictionary containing the classification results of the isomiR based on its modifications.

        This method compiles the relevant sequence modifications and classification details of the isomiR,
        including variations at the 5' and 3' ends, polymorphic changes, and whether the isomiR follows a template.
        
        The output includes:
        - The original sequence.
        - The computed classification name of the isomiR.
        - Various modifications related to the 5' and 3' ends, including overhangs, deletions, and nucleotide variations.
        - Polymorphic changes within the sequence.
        - Canonical miRNA information.
        - A list of associated hairpin precursor sequences.
        
        The format for nucleotide variations follows "X>Y", where:
        - `X` is the original nucleotide.
        - `Y` is the modified nucleotide.
        
        Example of a polymorphic change format:
        - "G5>A" indicates a guanine at position 5 being changed to adenine.

        Returns:
            dict: A dictionary containing the classification results of the isomiR.
        """

        # Execute variations related methods
        five_prime_variation = self.five_prime_variation()
        three_prime_variation = self.three_prime_variation()
        polymorphic_changes = self.polymorphic_changes()

        # Prepare the output for these methods
        if five_prime_variation is not None:
            five_prime_var_out = f'{five_prime_variation[1]}>{five_prime_variation[0]}'
        else:
            five_prime_var_out = five_prime_variation
        if three_prime_variation is not None:
            three_prime_var_out = f'{three_prime_variation[1]}>{three_prime_variation[0]}'
        else:
            three_prime_var_out = three_prime_variation
        if polymorphic_changes is not None:
            polymorphic_ch_out_list = []
            for change in polymorphic_changes:
                polymorphic_ch_out_list.append(f'{change[3]}{change[2]}>{change[0]}')
            polymorphic_ch_out = "|".join(polymorphic_ch_out_list)
        else:
            polymorphic_ch_out = polymorphic_changes

        return {
            "sequence": self.seq,
            "classification": self.isomir_name(),
            "five_prime_overhang": self.five_prime_overhang(),
            "five_prime_deletion": self.five_prime_deletion(),
            "five_prime_variation": five_prime_var_out,
            "three_prime_overhang": self.three_prime_overhang(),
            "three_prime_deletion": self.three_prime_deletion(),
            "three_prime_variation": three_prime_var_out,
            "polymorphic_changes": polymorphic_ch_out,
            "is_template": self.is_template(),
            "canonical_name": self.canonical_name,
            "canonical_sequence": self.canonical_seq,
            "hairpin_list": "|".join(self.hairpin_list)
        }
    def _remove_5prime_overhang(self):
        """Removes the 5' overhang from the original sequences."""
        seq_wo_five_prime_overhang = self.seq[self.seq_start:]
        can_wo_five_prime_overhang = self.canonical_seq[self.can_start:]
        return (seq_wo_five_prime_overhang, can_wo_five_prime_overhang)
    
    def _get_aligned_regions(self):
        """Gets the aligned regions between the query and reference sequences."""
        seq_aligned = self.one_alignment[0]
        can_aligned = self.one_alignment[1]
        return (seq_aligned, can_aligned)

    def _get_clean_aligned_regions(self):
        """Gets the aligned regions without dashes."""
        seq_aligned_clean = self.one_alignment[0].replace('-', '')
        can_aligned_clean = self.one_alignment[1].replace('-', '')
        return (seq_aligned_clean, can_aligned_clean)
    
    def _get_isomir_classification_suffix(self):
        """
        Generates an isomiR classification identifier based on its modifications.

        This method gathers results from other methods that indicate the different types of modifications
        present in the isomiR compared to the canonical reference sequence.

        It returns a string representing the modifications following the pattern:
        "0:0:0:0:0:0:0" where each position represents a specific modification:
        
        1. 5' overhang modifications
        2. 5' deletions
        3. 5' variations
        4. 3' overhang modifications
        5. 3' deletions
        6. 3' variations
        7. Polymorphic changes within the sequence
        
        If all positions contain "0", the sequence is identical to the canonical one.
        Otherwise, the specific modification is displayed. For 5' and 3' overhang modifications,
        the suffix "-T" is added if the change is template-based with respect to the precursor,
        and "-NT" if it is not template-based.

        Returns:
            str: An isomiR classification identifier in the format "X:X:X:X:X:X:X".
        """

        # Execute the required methods
        five_overhang = self.five_prime_overhang()
        five_deletion = self.five_prime_deletion()
        five_variation = self.five_prime_variation()
        three_overhang = self.three_prime_overhang()
        three_deletion = self.three_prime_deletion() 
        three_variation = self.three_prime_variation()
        polymorphic = self.polymorphic_changes()
        template = self.is_template()

        # Create an isomiR classification name
        classification_parts = []
        
        # 5' overhang
        if five_overhang is not None:
            if template:
                classification_parts.append(f'{five_overhang}-T')
            else:
                classification_parts.append(f'{five_overhang}-NT')
        else:
            classification_parts.append('0')

        # 5' deletion
        if five_deletion is not None:
            classification_parts.append(f'{five_deletion.lower()}')
        else:
            classification_parts.append('0')
        # 5' variation
        if five_variation is not None:
            classification_parts.append(f'{five_variation[1]}1{five_variation[0]}')
        else:
            classification_parts.append('0')
        # 3' overhang
        if three_overhang is not None:
            if template:
                classification_parts.append(f'{three_overhang}-T')
            else:
                classification_parts.append(f'{three_overhang}-NT')
        else:
            classification_parts.append('0')
        # 3' deletion
        if three_deletion is not None:
            classification_parts.append(f'{three_deletion.lower()}')
        else:
            classification_parts.append('0')
        # 3' variation
        if three_variation is not None:
            classification_parts.append(f'{three_variation[1]}{len(self.canonical_seq)}{three_variation[0]}')
        else:
            classification_parts.append('0')

        # Polymorphic changes
        if polymorphic is not None:
            # Iterate through the different changes
            for change in polymorphic:
                classification_parts.append(f'{change[2]}{change[3]}{change[0]}')
        else:
            classification_parts.append('0')
        
        # Check if the query sequence is canonical
        if len(classification_parts) == 0:
            classification_id = 'canonical'
        else:
            classification_id = ":".join(classification_parts)
        
        return classification_id


## FUNCTIONS
def extract_species_sequences(species_id, path_file_in, path_file_out, rna_to_dna=False):
    """
    This function extracts sequences from a FASTA file that match a 
    specific species based on the species ID. It then writes these 
    sequences to a new FASTA file and returns the number of sequences 
    extracted. If the argument `rna_to_dna` is set to True, the function 
    will replace uracil ('U') with thymine ('T') to convert RNA sequences 
    to DNA sequences.

    Parameters
    ----------
    path_file_in : str
        The path to the input FASTA file containing sequences.
    
    species_id : str
        The ID of the species to search for in the FASTA file headers.
    
    path_file_out : str
        The path to the output file where the matching sequences will 
        be written.

    rna_to_dna : bool, optional
        If True, converts RNA sequences (where 'U' is present) to DNA sequences 
        by replacing 'U' with 'T'. Default is False, meaning no conversion is made.

    Returns
    -------
    int
        The number of sequences that were extracted and written to the 
        output file. If no sequences are found, returns 0.
    
    Notes
    -----
    The function scans the headers of the input FASTA file for the 
    specified species ID. Only sequences whose headers contain the 
    species ID will be written to the output file.
    
    If `rna_to_dna` is True, the function will convert RNA sequences to 
    DNA by replacing 'U' with 'T'.
    """
    
    # Convert the identifier to lowercase
    species_id_lower = species_id.lower()

    # Initialize lists to hold sequences
    species_sequences = []

    # Iterate through fasta file
    for header, seq in parse_fasta_file(path_file_in):
        
        # Select only the important part of the header.
        clean_header = header.split(" ")[0]
        
        # Check if the header contains the species ID
        if species_id_lower in clean_header.lower():
            if rna_to_dna:
                # Convert RNA to DNA by replacing 'U' with 'T'
                species_sequences.append((clean_header, seq.replace('U', 'T')))
            else:
                species_sequences.append((clean_header, seq))
    
    # Get the number of selected sequences
    number_of_sequences = len(species_sequences)

    # Write the species sequences to the output file
    if number_of_sequences > 0:
        with open(path_file_out, "w") as species_file:
            for header, seq in species_sequences:
                species_file.write(f">{header}\n{seq}\n")

    return number_of_sequences


def get_diffexp_sequences(dea_file: str) -> pd.DataFrame:
    """
    This function extracts sequences from a differential expression analysis 
    (DEA) results file, typically generated by DESeq2 in R. It reads the file, 
    retrieves the sequences, and assigns invented names to them. A new DataFrame 
    is returned containing the invented names and their corresponding sequences.

    Parameters
    ----------
    dea_file : str
        Path to the DEA results file (usually a .tsv file) generated by DESeq2.
        The file should have sequences in the first column.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing two columns: 'Name' with invented names for the 
        sequences (e.g., Sequence1, Sequence2), and 'Sequence' with the actual 
        sequences extracted from the input file.
    """

    # Create a pandas dataframe
    df = pd.read_csv(dea_file, sep='\t')

    # Get the sequences from the table
    sequences = df.iloc[:, 0]

    # Create a list with invented names for the sequences
    names = [f"Sequence{i+1}" for i in range(len(sequences))]

    # Create a dataframe with invented names and the sequences
    sequences_df = pd.DataFrame({
        'Name': names,
        'Sequence': sequences
    })
    return sequences_df


def sequences_to_fasta(df: pd.DataFrame, fasta_file: str) -> None:
    """
    This function converts a DataFrame containing sequence names and sequences 
    into a FASTA file format. The DataFrame is iterated row by row, and for each 
    row, the sequence name is written as a header in the FASTA format followed 
    by the corresponding sequence.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame with two columns: 'Name' containing the sequence names and 
        'Sequence' containing the corresponding sequences.

    fasta_file : str
        Path to the output FASTA file where the sequences will be written.

    Returns
    -------
    None
    """

    # Open the output file
    with open(fasta_file, 'w') as file:
        # Iterate through the df rows
        for row in df.itertuples(index=False):
            # Write the header
            file.write(f">{row.Name}\n")
            # Write the sequence
            file.write(f"{row.Sequence}\n")  


def makeblastdb(input_file: str, db_output_dir: str, db_name: str) -> None:
    """
    This function creates a nucleotide BLAST database using the `makeblastdb` 
    command-line tool. It generates the database from an input file (such as a 
    FASTA file) and stores it in the specified output directory with the provided 
    database name.

    Parameters
    ----------
    input_file : str
        Path to the input file (e.g., a FASTA file) that will be used to create 
        the BLAST database.

    db_output_dir : str
        Directory where the BLAST database files will be saved.

    db_name : str
        Name to be assigned to the BLAST database. The function will append the 
        necessary file extensions for the database format.

    Returns
    -------
    None
        This function does not return any value. It executes the `makeblastdb` 
        command to create the database in the specified directory.

    Raises
    ------
    subprocess.CalledProcessError
        If the `makeblastdb` command fails, an error message will be printed.
    """

    # Create the makeblastdb command
    command = (
        f"makeblastdb -in {input_file} -dbtype nucl "
        f"-out {db_output_dir}/{db_name}"
    )
    
    # Ejecutar el comando con subprocess
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        print(f"Error creating the database with makeblastdb: {e}")
        sys.exit(1)


def isomir_blastn(query_file: str, reference_file: str, output_file: str, 
                  task: str = "blastn-short", outfmt: int = 6, word_size: int = 13, 
                  gapopen: int = 5, gapextend: int = 2, strand: str = "plus",
                  evalue: float = 0.01) -> pd.DataFrame:
    """
    This function performs a BLASTN search to find potential isomiRs by aligning 
    a query file (containing the sequences of interest) against a reference file. 
    The reference file is indexed using `makeblastdb`, and the resulting BLASTN 
    alignment is saved to an output file in a specified format. The function 
    processes the BLASTN output and returns the results as a DataFrame.

    Parameters
    ----------
    query_file : str
        Path to the query file containing sequences to be aligned (e.g., a FASTA file).

    reference_file : str
        Path to the reference file against which the query file will be aligned (e.g., a FASTA file).

    output_file : str
        Path to the output file where the BLASTN results will be saved.

    task : str, optional, default "blastn-short"
        The BLASTN task to be executed. Default is "blastn-short" for isomiR detection.

    outfmt : int, optional, default 6
        Output format for BLASTN results. Format 6 provides tabular output with standard fields.

    word_size : int, optional, default 13
        Word size for the BLASTN search. Default is 13, which is typically used for isomiR detection.

    gapopen : int, optional, default 5
        Cost to open a gap in the alignment. Default is 5.

    gapextend : int, optional, default 2
        Cost to extend a gap in the alignment. Default is 2.

    strand : str, optional, default "plus"
        The strand to search against. Options are "plus" or "both". Default is "plus".

    evalue : float, optional, default 0.05
        Expectation value (E-value) threshold for reporting hits. Default is 0.05.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the BLASTN results with the following columns:
        'Query', 'Subject', 'Identity', 'Length', 'Mismatches', 'Gap_Openings', 
        'Q_Start', 'Q_End', 'S_Start', 'S_End', 'E_value', 'Bit_Score'.

    Raises
    ------
    subprocess.CalledProcessError
        If the BLASTN command fails, an error message will be printed.
    """

    
    # Obtain the reference directory and filename
    ref_dir = os.path.dirname(reference_file)
    ref_basename = os.path.splitext(os.path.basename(reference_file))[0]
    
    # Index the reference
    makeblastdb(reference_file, f'{ref_dir}/{ref_basename}_idx', ref_basename)

    # Create the blastn command
    command = (
        f"blastn -task {task} -query {query_file} -db {f'{ref_dir}/{ref_basename}_idx/{ref_basename}'} "
        f"-outfmt {outfmt} -out {output_file} -word_size {word_size} "
        f"-gapopen {gapopen} -gapextend {gapextend} -strand {strand} "
        f"-evalue {evalue}"
    )
    
    # Run the command
    try:
        subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL)
        
    except subprocess.CalledProcessError as e:
        print(f"Error executing blastn: {e}")
        sys.exit(1)
    else:
        # Read the output file
        df = pd.read_csv(output_file, sep='\t', header=None)

        # Add colnames
        df.columns = ['Query', 'Subject', 'Identity', 'Length', 'Mismatches', 'Gap_Openings',
                    'Q_Start', 'Q_End', 'S_Start', 'S_End', 'E_value', 'Bit_Score']

    return df


def get_hairpin_sequences_from_fasta(seq_name: str, fasta: str) -> str:
    """
    This function retrieves the hairpin sequences from a FASTA file based on the 
    provided sequence name (typically the mature miRNA name). It uses an `awk` 
    command to format the FASTA file, then filters it using `grep` to find the 
    relevant sequence and returns it.

    Parameters
    ----------
    seq_name : str
        The name of the mature miRNA (e.g., 'miRNA1') whose corresponding hairpin 
        sequence is to be retrieved from the FASTA file.

    fasta : str
        Path to the input FASTA file containing the sequences.

    Returns
    -------
    str
        The hairpin sequence corresponding to the provided miRNA name, or `None` 
        if the sequence is not found.

    Raises
    ------
    subprocess.CalledProcessError
        If the command to extract the sequence fails, an error message will be printed.
    """
    try:
        # Build the command
        command = f"awk '/^>/ {{printf(\"%s%s\t\",(N>0?\"\\n\":\"\"),$0);N++;next;}} {{printf(\"%s\",$0);}} END {{printf(\"\\n\");}}' {fasta} | grep -P \"^>{seq_name}\\s*\" | cut -f2 "
        
        # Run the grep command
        result = subprocess.run(command, shell=True, text=True, capture_output=True)

        # If the sequence was found...
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return None

    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")
        return None


def add_hairpin_sequence(df: pd.DataFrame, mature_name_col: str, hairpin_file: str) -> pd.DataFrame:
    """
    This function adds the hairpin precursor sequences to a DataFrame based on 
    the provided names of mature miRNAs. The function first retrieves unique 
    mature miRNA names from the specified column, then queries the corresponding 
    hairpin sequences from a FASTA file and appends these sequences to the 
    DataFrame in a new column.

    Parameters
    ----------
    df : pd.DataFrame
        The input DataFrame containing miRNA information. The DataFrame must 
        include a column with names of mature miRNAs.

    mature_name_col : str
        The name of the column in the DataFrame that contains the names of 
        mature miRNAs (e.g., 'miRNA_name').

    hairpin_file : str
        Path to the FASTA file containing the precursor (hairpin) sequences. 
        Each sequence should correspond to a mature miRNA name.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with an additional column 'Hairpin_sequence', 
        which contains the corresponding hairpin sequences for each mature miRNA. 
        If no hairpin sequence is found, the cell will contain `None`.

    Notes
    -----
    The function assumes that the mature miRNA names in the input DataFrame 
    follow a specific format where the precursor name can be derived by 
    splitting the mature miRNA name (e.g., 'miR-123a' becomes 'MIR-123').
    """

    # Get a list of uniue mature miRNA names
    uniq_miRNAs_mat = df[mature_name_col].unique().tolist()

    # Iterate through mature miRNA names
    hairpin_seq_list = []
    for mat in uniq_miRNAs_mat:

        # Create the precursor seq name
        hairpin_seq_name = "-".join([mat.split('-')[0], mat.split('-')[1].split('.')[0]]).replace('miR', 'MIR')

        # Get the hairpin sequence
        hairpin_seq = get_hairpin_sequences_from_fasta(hairpin_seq_name, hairpin_file)

        # Add the hairpin sequence to the list
        hairpin_seq_list.append(hairpin_seq)

    # Iterate through miRNA-names and hairpin-sequences lists
    for mat, hairpin_seq in zip(uniq_miRNAs_mat, hairpin_seq_list):
        # Find the index of those rows with the miRNA name in the column 'mature_name_col'
        idx = df[df[mature_name_col] == mat].index

        # Add the hairpin sequences to the dataframe
        if not idx.empty:
            df.loc[idx, 'Hairpin_sequence'] = hairpin_seq
    return df


def parse_fasta_file(fasta_file: str):
    """
    This function parse a fasta file returning a generator object containing
    the headers and the sequences. If this object is iterated outside the
    function, we can access a new sequence and its header at each iteration
    of the loop.

    Parameters
    ----------
    fasta_file : str
        Absolute path of the fasta file
    """
    try:
        with open(fasta_file) as fasta:
            iterator = (x[1] for x in itertools.groupby (fasta, lambda line: line[0] == '>'))
            for header in iterator:
                # Drop the ">"
                headerStr = header.__next__()[1:].strip()
                # join all sequence lines to one
                seq = ''.join(s.strip() for s in iterator.__next__())
                yield(headerStr,seq)

    except StopIteration:
        print('Error. Check that the files entered are fasta')
        sys.exit(1)


def collapse_fasta_sequences(fasta, output_file=''):
    """
    This function processes a FASTA file and creates a new FASTA file where
    identical sequences are collapsed together. For each unique sequence, all
    headers corresponding to that sequence are combined into a single header line,
    separated by a pipe (`|`). The new sequences are then written to the specified
    output file if the `output_file` parameter is provided. If no `output_file` is
    provided, no file will be written.

    Parameters
    ----------
    fasta : str
        Absolute path of the input FASTA file that contains sequences and headers.
    output_file : str, optional
        Absolute path of the output FASTA file where the collapsed sequences will
        be written. If not provided (defaults to an empty string), no output file
        will be generated.

    Returns
    -------
    dict
        A dictionary where the keys are sequences (as strings) and the values
        are lists of headers that correspond to each sequence. Each list contains
        the headers of all entries that have the same sequence.
    """
    
    # Seq-[Headers] dictionary
    sequence_to_headers = {}
    
    # Iterate through fasta file
    for header, seq in parse_fasta_file(fasta):

        # If the sequence has been previosuly added...
        if seq in sequence_to_headers:
            sequence_to_headers[seq].append(header)
        # Add the sequences into the dictionary
        else:
            sequence_to_headers[seq] = [header]
    
    # Write the new fasta file and modify the header_to_sequence dictionary
    header_to_sequence = {}
    with open(output_file, 'w') as out:
        for seq, headers in sequence_to_headers.items():
            out.write(f">{'|'.join(headers)}\n{seq}\n")
            header_to_sequence['|'.join(headers)] = seq

    return header_to_sequence


def fasta_to_dic(fasta):
    """
    This function processes a FASTA file and converts it into a dictionary
    where each header is associated with its corresponding sequence. The
    headers in the file become the keys in the dictionary, and the sequences
    are stored as the corresponding values.

    Parameters
    ----------
    fasta : str
        Absolute path of the input FASTA file that contains sequences and headers.

    Returns
    -------
    dict
        A dictionary where the keys are the headers (as strings) and the values
        are the sequences (as strings) corresponding to each header in the FASTA file.
    """
    
    # Header-seq dictionary
    sequences_dic = {}

    # Iterate through fasta file
    for header, seq in parse_fasta_file(fasta):
        # Add the sequences into the dictionary
        sequences_dic[header] = seq.replace('U', 'T')

    return sequences_dic


def get_miRNA_base(miRNA: str) -> str:
    """
    This function extracts the base miRNA name from a given mature miRNA name. 
    It uses a regular expression to find and return the miRNA base name, 
    which typically follows the pattern 'miR' followed by digits (and optionally a version number).

    Parameters
    ----------
    miRNA : str
        The name of the mature miRNA, which may include a version number (e.g., 'miR-21', 'miR-21.1').

    Returns
    -------
    str
        The base name of the miRNA, extracted from the input string (e.g., 'miR21' from 'miR-21' or 'miR21.1').

    Raises
    ------
    AttributeError
        If the input string does not match the expected miRNA name pattern, 
        the function will raise an `AttributeError` when trying to access the group.
    """

    return re.search(r'(miR\d+(\.\d+)?)', miRNA, re.IGNORECASE).group(1)


def best_alignment_within_hairpin(df: pd.DataFrame) -> pd.DataFrame:
    """
    This function resolves cases where a potential isomiR aligns to 
    two or more canonical mature miRNAs that belong to the same 
    precursor. It compares the alignments based on e-value and selects 
    the one with the best match (lowest e-value). Sequences that do 
    not match this scenario are not modified.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the results of a BLAST search, where 
        each row represents an alignment between a potential isomiR and 
        a canonical mature miRNA. The DataFrame must include the following 
        columns:
        - 'Query': the name of the isomiR (potential sequence).
        - 'Subject': the name of the canonical mature miRNA(s), which may 
          include variant names.
        - 'E_value': the e-value for the alignment, used to determine the 
          quality of the match.

    Returns
    -------
    pd.DataFrame
        The input DataFrame with only the best matching alignment 
        (lowest e-value) retained for those isomiRs that align with 
        multiple mature miRNAs from the same precursor. For other 
        sequences, no changes are made. Rows with less optimal alignments 
        (based on e-value) are removed for those isomiRs that meet this 
        condition, ensuring that only the most accurate match for each 
        isomiR is kept.

    Notes
    -----
    This function is useful in cases where a potential isomiR aligns 
    with multiple mature miRNAs from the same precursor, such as when 
    two mature miRNAs are close to each other in the precursor sequence 
    and an isomiR aligns to both. For those isomiRs, the function ensures 
    that only the best matching mature miRNA (based on the lowest e-value) 
    is retained, while others are removed. Sequences that do not align to 
    multiple miRNAs from the same precursor are not affected.
    """

    #  Identify sequences with more than one match
    duplicated_queries = df["Query"].value_counts()
    duplicated_queries = duplicated_queries[duplicated_queries > 1].index 

    # Iterate through sequences with more than one match
    for sequence in duplicated_queries:
        # Get the rows asscoiated to 'sequence'
        df_sequence = df[df["Query"] == sequence]

        # Obtain the miRNA names from 'Subject' columnas.
        subject_miRNAs = df_sequence["Subject"].tolist()
        
        # Iterate through subject names to find variants
        miRNA_variants_list = [] 
        for miRNA_set in subject_miRNAs:
            # Find if exists a miRNA variant
            miRNA_list = miRNA_set.split('|')
            miRNA_variants = [(miRNA, get_miRNA_base(miRNA)) for miRNA in miRNA_list if re.match(r'\S+\.\d+', miRNA)]
            miRNA_variants_list.extend(miRNA_variants)

        # Variants groups
        groups = {}
        
        # Group tuples by miRNA base name
        for complete_miRNA, miRNA in miRNA_variants_list:
            # Check that it is a variant within the same hairpin
            if re.search(r'\.(\d+)$', miRNA):
                # Get the base
                base = '.'.join(miRNA.split('.')[:-1])
                # Save results into the dictionary
                if base not in groups:
                    groups[base] = [] 
                groups[base].append(complete_miRNA)

        # Check which variant has the lowest e-value and remove the rest
        for miRNAs_names in groups.values():
            filtered_df = df[(df['Query'] == sequence) & 
                        (df['Subject'].str.contains("|".join(miRNAs_names)))]
            # Select the row with the lowest evalue
            min_e_value_row = filtered_df.loc[filtered_df['E_value'].idxmin()]
            
            # Remove the unselected rows
            df = df.drop(filtered_df[filtered_df['E_value'] != min_e_value_row['E_value']].index)
    return df


def classify_mono_or_multiple(value):
    """
    This function classifies a given value as either 'Mono' or 'Multiple'. 

    If the value is 'None', it returns None. If the value consists of a 
    single unique character, it returns that character (e.g., 'A', 'T', 'G', etc.). 
    Otherwise, it classifies the value as 'Multiple' when there are multiple 
    distinct characters in the input.

    Parameters
    ----------
    value : str
        A string representing a nucleotide or sequence of nucleotides. It can be:
        - A string containing a single character (e.g., 'A', 'T', 'G', 'C').
        - A string containing multiple distinct characters (e.g., 'AT', 'AG', etc.).
        - The string "None" to represent a null value.

    Returns
    -------
    str or None
        - If the input value is "None", returns None.
        - If the input value consists of a single unique character, returns that character.
        - If the input value consists of multiple distinct characters, returns 'Multiple'.
    """
    if value == "None":
        return None
    elif len(set(value)) == 1:
        return value 
    else:
        return "Multiple"


def isomirs_summary_table(df):
    """
    This function generates a summary table that provides information about the 
    distribution of various isomiRs, their variants, polymorphic changes, and 
    mixed isomiRs based on a given DataFrame. The summary includes counts for 
    annotated sequences, isomiRs, specific variant types (such as 5' and 3' overhangs, 
    deletions, and variations), polymorphic changes, and mixed isomiRs.

    The function processes the DataFrame to:
    - Count the number of annotated sequences and isomiRs.
    - Classify isomiRs into various categories based on the presence of 5' and 3' 
      overhangs, deletions, variations, and polymorphic changes.
    - Identify isomiRs with polymorphic changes and classify the types of 
      polymorphic changes (e.g., A>T, G>A).
    - Count mixed isomiRs that have more than one variant.
    - Count canonical miRNAs where no variants are present.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing the isomiR data. The DataFrame must include the following columns:
        - 'seq': the sequence of the annotated isomiRs (optional).
        - 'classification': the classification of the isomiRs (e.g., isomiR type).
        - Columns representing variants: 'five_prime_overhang', 'five_prime_deletion', 
          'five_prime_variation', 'three_prime_overhang', 'three_prime_deletion', 
          'three_prime_variation', 'polymorphic_changes'.

    Returns
    -------
    pd.DataFrame
        A DataFrame containing the summary of isomiRs
    
    Notes
    -----
    The output DataFrame contains the following information:
        - 'Num_annotated_seq': Number of unique annotated sequences.
        - 'Num_isomirs': Number of isomiRs (no unique sequences).
        - 'Num_canonical': Number of canonical miRNAs.
        - 'Num_five_prime_overhang', 'Num_five_prime_deletion', 'Num_five_prime_variation', 
          'Num_three_prime_overhang', 'Num_three_prime_deletion', 'Num_three_prime_variation': 
          Counts of isomiRs with each specific variant type.
        - 'Num_polymorphic_changes': Number of isomiRs with polymorphic changes.
        - Polymorphic change columns (e.g., 'polymorphic_A>T', 'polymorphic_G>A') for 
          counting the specific polymorphic changes.
        - 'Num_mixed': Number of mixed isomiRs with more than one variant.
        - This DataFrame also contains several columns that count the number
          of specific changes for each type of variant (such as addition,
          deletion, or change of specific nucleotides).
    """

    # Create a summary dictionary
    summary = {}
    
    # Get the number of annotated sequences
    if 'seq' in df.columns:
        summary['Num_annotated_seq'] = df['seq'].nunique()
    else:
        summary['Num_annotated_seq'] = None
    
    # Get the number of isomirs
    summary['Num_isomirs'] = df['classification'].nunique()

    
    ## 2. 5' AND 3' ISOMIRS
    ###########################################################################

    # Columns representing each variant (except polymorphic changes)
    variant_cols = [
        'five_prime_overhang', 'five_prime_deletion', 'five_prime_variation',
        'three_prime_overhang', 'three_prime_deletion', 'three_prime_variation',
        'polymorphic_changes'
    ]
    
    # Count de number of sequences for each variant
    for col in variant_cols:

        # Skip 'polymorphic_changes'
        if col != 'polymorphic_changes':
            # Select those isomirs with only 'col' changes
            variant_mask = (df[col].ne("None")) & (df[[c for c in variant_cols if c != col]].eq("None").all(axis=1))
            filt_df = df[variant_mask]
            summary[f'Num_{col}'] = filt_df.shape[0]
            
            # Classify in mono or multiple
            values = filt_df.loc[filt_df[col].ne("None"), col].apply(classify_mono_or_multiple)
            summary[f'{col}_mono_A'] = (values == 'A').sum()
            summary[f'{col}_mono_T'] = (values == 'T').sum()
            summary[f'{col}_mono_C'] = (values == 'C').sum()
            summary[f'{col}_mono_G'] = (values == 'G').sum()
            summary[f'{col}_multiple'] = (values == 'Multiple').sum()
    
    ## 2. POLYMORPHIC CHANGES
    ###########################################################################

    # Create the required columns for polymorphic changes summary
    nucleotides = ['A', 'C', 'T', 'G']
    change_columns = [f'polymorphic_{x}>{y}' for x, y in itertools.product(nucleotides, repeat=2) if x != y]
    for col in change_columns:
        summary[col] = 0
    summary['polymorphic_others'] = 0

    # Select only those isomirs with polymorphic changes
    variant_mask = (df['polymorphic_changes'].ne("None")) & (df[variant_cols].eq("None").all(axis=1))
    filt_df = df[variant_mask]
    summary[f'Num_polymorphic_changes'] = filt_df.shape[0]

    # Iterate through those isomirs with polymorphic changes
    for _, row in filt_df.iterrows():
        # Get the change
        change = row['polymorphic_changes']
        # Check if the element is different to None
        if ">" in change:
            # Get the original and new nucleotides
            parts = change.split('>')
            # Remove position from original nucleotides string
            original_nuc = re.sub(r'\d+', '', parts[0])
            # Verify if the change implies more than one nucleotide
            if len(original_nuc) == 1 and len(parts[1]) == 1:
                change_pair = f'polymorphic_{original_nuc}>{parts[1]}'
                if change_pair in change_columns:
                    summary[change_pair] += 1
            else:
                summary['polymorphic_others'] += 1
    
    ## 3. MIXED ISOMIRS
    ###########################################################################

    # Count the number of mixed isomirs
    summary['Num_mixed'] = df[df[variant_cols].ne("None").sum(axis=1) > 1].shape[0]

    ## 4. CANONICAL MIRNAS
    ###########################################################################

    # Count the number of canonical miRNAs
    summary['Num_canonical'] = df[df[variant_cols].eq("None").all(axis=1)].shape[0]

    ## 5. SORT THE FINAL DATAFRAME
    ###########################################################################

    # Create the output dataframe
    summary_df = pd.DataFrame([summary])

    # First columns
    base_columns = [
        'Num_annotated_seq', 'Num_isomirs', 'Num_canonical', 
        'Num_five_prime_overhang', 'Num_five_prime_deletion', 'Num_five_prime_variation',
        'Num_three_prime_overhang', 'Num_three_prime_deletion', 'Num_three_prime_variation',
        'Num_polymorphic_changes', 'Num_mixed'
    ]

    # Get the dataframe columns
    all_columns = summary_df.columns.tolist()

    # Create the new order for the columns of the final dataframe
    ordered_columns = [col for col in base_columns if col in all_columns]
    remaining_columns = [col for col in all_columns if col not in base_columns]
    final_column_order = ordered_columns + remaining_columns

    # Change the order of the columns in the dataframe
    summary_df = summary_df[final_column_order]
    
    return summary_df


def get_mirna_family(annotation: str) -> str:
    """
    This function extracts the miRNA family name from a given miRNA annotation string.
    The annotation string is expected to follow a naming convention where the family 
    name is composed of a precursor identifier (e.g., "miR") followed by a number 
    and, optionally, a letter (e.g., "156a" or "156").

    Parameters
    ----------
    annotation : str
        A string representing the miRNA annotation. Examples: ath-miR1a-5p,
        ath-miR-1a-5p, miR1a-3p, miR-1a-3p, cel-let7, etc.

    Returns
    -------
    str
        The miRNA family name derived from the annotation Examples: miR1, let7,
        etc.

    Notes
    -----
    This function is useful for parsing miRNA annotations where the family 
    name needs to be extracted for further analysis, such as grouping or 
    classification of miRNAs.
    """

    # Divide the miRNA name in elements
    annot_elements = annotation.split('-')

    # Find the first element in the list tha contains a number
    for i, element in enumerate(annot_elements):
        # Find a number in annot_elements
        match_num = re.search(r'\d', element)

        # If a number is found...
        if match_num:
            # Get the position of the first number found
            number_pos = match_num.start()
            
            # The number is in the first position (e.g 156a)
            if number_pos == 0:
                # Found the first letter in the string
                match_letter = re.search(r'[a-zA-Z]', element)

                # If a letter is found...
                if match_letter:
                    letter_pos = match_letter.start()
                    # Get a string with only the numeric part
                    numeric_part = element[:letter_pos]
                else:
                    numeric_part = element

                # Build the miRNA family name
                family_id = f'{annot_elements[i-1]}{numeric_part}'
            
            # The number is NOT in the first position (e.g miR156a)
            else:
                # Get the numeric raw part (e.g. 156ab. clean= 156)
                numeric_raw_part = element[number_pos:]

                # Found the first letter in the string
                match_letter = re.search(r'[a-zA-Z]', numeric_raw_part)
                
                # If a letter is found...
                if match_letter:
                    # Get the position of the letter
                    letter_pos = match_letter.start()
                    # Build the miRNA family name
                    family_id = element[:number_pos + letter_pos]
                else:
                    family_id = element
            break

    return family_id

## MAIN
def main():
    '''
    Main program
    '''

    # Define the argument parser
    parser = argparse.ArgumentParser(prog='isomiRs_identification', 
                                    description='''This program performs the
                                    identification and classification of
                                    isomiRs using sequences from FASTA files
                                    or differential expression analysis
                                    results.''',  
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Argument for specifying the output file name (id)
    parser.add_argument('-id', '--id', type=str, required=True, 
                        help='File ID.')

    # Argument for specifying the species ID (species_id)
    parser.add_argument('-s', '--species-id', type=str, required=True, 
                        help='Species ID.')

    # Argument for specifying the input file (either a FASTA or a TSV file)
    parser.add_argument('-i', '--input', type=str, required=True, 
                        help='Input file containing either a FASTA file or a'
                        'TSV file with differential expression results.')

    # Argument for specifying the mature miRNA sequences FASTA file
    parser.add_argument('-m', '--mature-db', type=str, required=True, 
                        help='FASTA file containing mature miRNA sequences'
                        'from a reference database (e.g., miRBase, PmiREN, sRNAanno).')

    # Argument for specifying the hairpin miRNA sequences FASTA file
    parser.add_argument('-p', '--hairpin-db', type=str, required=True, 
                        help='FASTA file containing hairpin miRNA sequences'
                        'from a reference database (e.g., miRBase, PmiREN, sRNAanno).')

    # Parse the arguments
    args = parser.parse_args()

    ###########################################################################
    #                        1. CHECK ARGUMENTS                               #
    ###########################################################################

    try:
        # Store the arguments in variables
        output_id = args.id
        species_id = args.species_id
        input_file = args.input
        mature_fasta = args.mature_db
        hairpin_fasta = args.hairpin_db

    except Exception as e:
        print(f'ERROR: {e}')
        print('Please check the provided parameters or make sure no parameter is missing.')
        parser.print_help()
        sys.exit(1)

    # Create temporary directory
    os.makedirs('./tmp', exist_ok=True)

    ###########################################################################
    ##  1. PREPARE THE SEQUENCES OF THE DATABASE                             ##
    ###########################################################################

    # DB temporary files
    tmp_mature_db = f'./tmp/{species_id}_db_mature.fasta'
    tmp_mature_db_collapse = f'./tmp/{species_id}_db_mature_collapse.fasta'
    tmp_hairpin_db = f'./tmp/{species_id}_db_hairpin.fasta'

    # Other temporary files
    tmp_query_fasta = './tmp/input_tmp.fasta'
    tmp_blast_results = './tmp/blast_results_tmp.fasta'

    # Create a fasta file only with the "species" sequences (mature)
    number_of_sp_seqs = extract_species_sequences(species_id, mature_fasta, tmp_mature_db, True)

    # Check if the species exists in the database
    if number_of_sp_seqs > 0:
        
        # Create a fasta file only with the "species" sequences (hairpin)
        number_of_sp_seqs = extract_species_sequences(species_id, hairpin_fasta, tmp_hairpin_db, True)

        # Collapse canonical miRNAs fasta
        collapsed_canonical_fasta_dic = collapse_fasta_sequences(tmp_mature_db, tmp_mature_db_collapse)

        # Save the hairpin sequences into a dictionary
        hairpin_dic = fasta_to_dic(tmp_hairpin_db)

        #######################################################################
        ##  2. RUN BLASTN TO IDENTIFY ISOMIRS                                ##
        #######################################################################

        # Read the DEA table and obtain their sequences
        df = get_diffexp_sequences(input_file)

        # Create a fasta file with these sequences
        sequences_to_fasta(df, tmp_query_fasta)

        # Execute the blastn command
        blatn_results_df = isomir_blastn(tmp_query_fasta, tmp_mature_db_collapse, tmp_blast_results)

        #######################################################################
        ##  3. ADD THE CANONICAL AND HAIRPIN SEQUENCES TO THE BLASTN RESULTS ##
        #######################################################################

        # 1. Filtrar las secuencias duplicadas en la columna "Query"
        duplicated_queries = blatn_results_df["Query"].value_counts()
        duplicated_queries = duplicated_queries[duplicated_queries > 1].index
        
        # 2. Crear una nueva columna que elimine la versión del miRNA (e.g., ath-miR161.1 -> ath-miR161)
        blatn_results_df["Subject_Base"] = blatn_results_df["Subject"].apply(lambda x: re.sub(r'\.\d+', '', x))  

        # Select the best alignment if the sequence aligns with different variants
        # that belong to the same precursor.
        blatn_results_df_clean = best_alignment_within_hairpin(blatn_results_df)

        # Get the sequences of each match
        query_sequences = []
        subject_sequences = []
        hairpin_sequences = []
        for row in blatn_results_df_clean.itertuples(index=False):
            # Get the sequences names (query and subject)
            query = row.Query
            subject = row.Subject

            # Add the sequences to the lists
            query_sequences.append(df.loc[df['Name'] == query, 'Sequence'].values[0])
            subject_sequences.append(collapsed_canonical_fasta_dic[subject])

            #hairpin_sequences.append("|".join([hairpin_dic["-".join(name.split("-")[0:2])] for name in subject.split('|')]))
            hairpin_sequences.append("|".join([
                hairpin_dic["-".join(re.split(r'[.-]', name)[0:2])] for name in subject.split('|')
            ]))

        # Add the new columns 
        blatn_results_df_clean['Query_sequence'] = query_sequences
        blatn_results_df_clean['Subject_sequence'] = subject_sequences
        blatn_results_df_clean['Hairpin_sequences'] = hairpin_sequences

        # Filtrar las filas donde Query_sequence es igual a Subject_sequence y obtener la lista
        query_canonical_sequences = blatn_results_df_clean.loc[blatn_results_df_clean['Query_sequence'] == blatn_results_df_clean['Subject_sequence'], 'Query_sequence'].tolist()

        # Remove the rows whose sequences have already been previously annotated as canonical.
        isomirs_sequences_df = blatn_results_df_clean[
            ~((blatn_results_df_clean['Query_sequence'].isin(query_canonical_sequences)) & 
            (blatn_results_df_clean['Query_sequence'] != blatn_results_df_clean['Subject_sequence']))
        ]
        
        #######################################################################
        ##  3. ISOMIR CLASIFICATION                                          ##
        #######################################################################

        ## 3.1 Classificate the isomiRs
        ########################################################################

        data = []
        # isomiR classification
        for row in isomirs_sequences_df.itertuples(index=False):

            # Get the required sequences
            query_seq = row.Query_sequence
            canonical_name = row.Subject
            canonical_Seq = row.Subject_sequence
            hairpin_seq = row.Hairpin_sequences.split("|")

            # Classification
            isomirs = IsomiRs(query_seq, canonical_name, canonical_Seq, hairpin_seq)
            classification_res = isomirs.isomir_classification_results()
            classification_res_clean = {k: v if v is not None else "None" for k, v in classification_res.items()}
            data.append(classification_res_clean)

        # Save the results table
        isomir_classification_df = pd.DataFrame(data)


        ## 3.2 Get the miRNA family name for each isomiR
        ########################################################################

        # Create a list with the names of the canonical miRNAs
        canonical_list_names = isomir_classification_df['canonical_name'].tolist()
        
        # Iterate through canonical_list_names ['gma-miR167e|gma-miR167f',...]
        miRNA_families = []
        for canonical_names in canonical_list_names:
            # Split canonical names (e.g. 'gma-miR167e|gma-miR167f' -> ['gma-miR167e', 'gma-miR167f'])
            canonical_names_list = canonical_names.split("|")
            # Get unique sorted family IDs
            unique_family_ids = sorted({get_mirna_family(canonical_name) for canonical_name in canonical_names_list})
            # Construct final family name based on the length of unique family IDs
            final_family_name = "/".join(unique_family_ids) if len(unique_family_ids) > 1 else unique_family_ids[0]
            miRNA_families.append(final_family_name)

        # Add the miRNA family ids to the isomir classification dataframe
        isomir_classification_df['miRNA_family'] = miRNA_families

        ## 3.3 Combine isomiR classification results with other dataframes
        ########################################################################

        # Add the blastn results to the output dataframe
        isomir_class_with_blast_df = pd.merge(isomir_classification_df, blatn_results_df_clean, 
                  left_on=['sequence', 'canonical_sequence'], 
                  right_on=['Query_sequence', 'Subject_sequence'], 
                  how='left')

        # Remove undesired sequences
        isomir_class_with_blast_df = isomir_class_with_blast_df.drop(
            columns=['Query', 'Subject', 'Subject_Base', 'Query_sequence',
                     'Subject_sequence', 'Hairpin_sequences']
            )

        # Create the output files        
        isomir_class_with_blast_df.to_csv(f'{output_id}.blast_isomirs.tsv', sep="\t", index=False)

        # Create a dataframe with the DEA results
        dea_df = pd.read_csv(input_file, sep='\t')

        # Realizar la unión conservando todas las coincidencias de df2
        dea_isomirs_annotated_df= dea_df.merge(isomir_classification_df, left_on='seq', right_on='sequence', how='right')

        # Remove undesired sequences
        dea_isomirs_annotated_df = dea_isomirs_annotated_df.drop(columns=['sequence'])

        # Mostrar resultado
        dea_isomirs_annotated_df.to_csv(f'{output_id}.dea_isomirs.tsv', sep="\t", index=False)


        #######################################################################
        ##  4. ISOMIR CLASSIFICATION SUMMARY TABLE                           ##
        #######################################################################

        # Create the summary table 
        summary_df = isomirs_summary_table(dea_isomirs_annotated_df)

        # Save the summary table
        summary_df.to_csv(f'{output_id}.summary_isomirs.tsv', sep="\t", index=False)
    
    # Remove temporary directory
    #os.system('rm -r tmp')

if __name__ == "__main__":
    main()
