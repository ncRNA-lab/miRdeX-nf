#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   isomiR_Classification.py
#
#   This program identifies and classifies potential isomiRs based on alignment
#   results obtained from a previous BLASTN search.
#
#   The script uses a BLASTN result file containing alignment results between
#   query sequences and a miRNA database, including both mature miRNAs and
#   miRNA precursors. It then classifies the sequences that aligned to a
#   canonical miRNA according to the types of modifications they present.
#
#   The classification follows the mirGFF3 standard described at:
#   https://github.com/miRTop/mirGFF3/blob/master/definition.md
#
#   Supported isomiR modification types include:
#
#   - iso_5p:+/-N           Shift in the 5' end of the sequence.
#   - iso_3p:+/-N           Shift in the 3' end of the sequence.
#   - iso_add5p:N           Non-templated nucleotide additions at the 5' end.
#   - iso_add3p:N           Non-templated nucleotide additions at the 3' end.
#   - iso_snv_seed          Substitutions in positions 2–7 (seed region).
#   - iso_snv_central       Substitutions in positions 9–12.
#   - iso_snv_central_offset Substitutions in position 8.
#   - iso_snv_central_supp  Substitutions in positions 13–17.
#   - iso_snv               Substitutions outside the regions above.
#
#   The output is a DataFrame containing classification results for each
#   aligned sequence, which can later be used for downstream analysis or
#   exported to formats like GFF3.
#
#   Authors: Antonio González Sánchez
#   Date: March 17, 2025
#   Version: 1.0
#
#******************************************************************************

## IMPORT MODULES

import argparse
from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
from datetime import datetime
from mintplates import convert
import pandas as pd
import re
import sys
import os

## CLASSES

class IsomiRs:
    """
    A class for the classification of IsomiR sequences in relation to a
    canonical miRNA reference sequence.

    This class provides methods to classify a sequence as either an IsomiR 
    or a reference miRNA based on the presence of 5' and 3' overhangs, 
    deletions, substitutions, and single nucleotide variations (SNVs).
    It performs alignment with the canonical miRNA sequence and identifies 
    variations that distinguish the IsomiR from the reference miRNA.

    Attributes
    ----------
    seq : str
        The sequence to analyze (IsomiR).
    
    mirna_seq : str
        The canonical reference miRNA sequence.
    
    precursor_seq : str
        The precursor sequence (from which the miRNA is derived).
    
    alignment_done : bool
        A flag indicating whether the alignment between the IsomiR sequence 
        and the canonical miRNA has been completed.
    
    one_alignment_seq : tuple or None
        Stores the result of the alignment between the IsomiR and canonical 
        miRNA reference sequence. It contains two elements: the query 
        (IsomiR) sequence and the reference sequence.
    
    seq_mirna_start : int or None
        The start position of the IsomiR sequence in the reference miRNA 
        alignment.
    
    seq_mirna_ref_start : int or None
        The start position of the reference miRNA in the alignment.
    
    seq_pre_start : int or None
        The start position of the IsomiR sequence in the precursor alignment.
    
    seq_pre_ref_start : int or None
        The start position of the precursor reference miRNA in the alignment.

    Methods
    -------
    __init__(seq, mirna_seq, precursor_seq)
        Initializes the IsomiRs class for analyzing IsomiR sequences by 
        providing the sequence to be analyzed, the reference miRNA sequence, 
        and the precursor sequence.
    
    get_classification()
        Classifies the given sequence as either an IsomiR or a reference 
        miRNA based on the identified variations.

    get_variant()
        Identifies and returns the variations present in the IsomiR sequence 
        (such as overhangs, deletions, and substitutions).
    
    five_prime_overhang()
        Returns the 5' overhang sequence if it exists.
    
    five_prime_deletion()
        Returns the 5' deleted sequence if it exists.
    
    three_prime_overhang()
        Returns the 3' overhang sequence if it exists.
    
    three_prime_deletion()
        Returns the 3' deleted sequence if it exists.
    
    substitutions()
        Returns a list of polymorphic changes (substitutions) between the 
        IsomiR sequence and the canonical miRNA sequence.
    
    is_template()
        Checks if the IsomiR sequence is present in any hairpin structure 
        from the precursor sequence.

    """

    def __init__(self, seq, mirna_seq, precursor_seq):
        """
        Initializes the IsomiRs class for the analysis and classification 
        of IsomiR sequences in relation to a canonical miRNA reference 
        sequence.

        :param seq: str
            The sequence to classify (IsomiR). This is the miRNA sequence 
            that will be compared to the canonical miRNA for variation 
            analysis.
            
        :param mirna_seq: str
            The reference miRNA sequence. This is used for comparison to detect
            variations in the provided IsomiR sequence.
            
        :param precursor_seq: str
            The precursor sequence from which the canonical miRNA and 
            the IsomiR are derived. It helps to understand the context 
            in which the miRNA sequence is embedded.
        
        Initializes the following attributes:
        - seq: The IsomiR sequence.
        - mirna_seq: The canonical reference miRNA sequence.
        - precursor_seq: The precursor sequence.
        - one_alignment_seq: Holds the result of the alignment between
            the IsomiR and the canonical miRNA.
        - seq_mirna_start: The start position of the IsomiR sequence in 
            the reference miRNA alignment.
        - seq_mirna_ref_start: The start position of the reference miRNA
            in the alignment.
        - seq_pre_start: The start position of the IsomiR sequence in the
            precursor alignment.
        - seq_pre_ref_start: The start position of the precursor reference
            miRNA in the alignment.
        """
        # Input attributes
        self.seq = seq
        self.mirna_seq = mirna_seq
        self.precursor_seq = precursor_seq
        
        # Alignment attribute
        self.one_alignment_seq = None
        self.seq_mirna_start = None
        self.seq_mirna_ref_start = None
        self.seq_pre_start = None
        self.seq_pre_ref_start = None


    def get_isomir_id(self):
        """
        Create a unique identifier for the isomiR.

        This method executes the "convert" function from the Python script
        mintplate.py for the construction of the identifier. This script has
        been obtained from:
        https://github.com/miRTop/mirtop/blob/dev/mirtop/mirna/mintplates.py
        Additionally, this code is inspired by MINTplate:
        Inspired by MINTplate: https://cm.jefferson.edu/MINTbase
        https://github.com/TJU-CMC-Org/MINTmap/tree/master/MINTplates

        Args:
            *seq(str)*: nucleotides sequences.

        Returns:
            *idName(str)*: unique identifier for the sequence.
        """
        try:
            identifier = convert(self.seq, True, 'iso')
        except KeyError as e:
            print(f'ERROR: {e}')
            raise
        return identifier

    def align_to_miRNA(self):
        """
        This method performs a local sequence alignment of the object's
        sequence (`seq`) with a reference miRNA sequence (`mirna_seq`) 
        using the Smith-Waterman algorithm. The alignment result is stored, 
        and specific start positions of the aligned sequences are extracted 
        for further use.

        Parameters
        ----------
        None

        Returns
        -------
        alignment : object
            The result of the Smith-Waterman alignment, which includes the
            aligned sequences and their respective positions.

        Attributes
        ----------
        one_alignment_seq : object
            Stores the result of the alignment.
        
        seq_mirna_start : int
            The start position of the sequence in the reference miRNA 
            alignment.
        
        seq_mirna_ref_start : int
            The start position of the reference miRNA in the alignment.
        """
        alignment = self._swalignment(self.seq, self.mirna_seq)
        self.one_alignment_seq = alignment
        self.seq_mirna_start = alignment.aligned[0][0][0]
        self.seq_mirna_ref_start = alignment.aligned[1][0][0]
        return alignment

    def five_prime_overhang(self):
        """
        This method retrieves the 5' overhang sequence from the object's
        sequence (`seq`) if it exists, based on the alignment with the 
        reference miRNA sequence.

        Parameters
        ----------
        None

        Returns
        -------
        five_overhang_seq : str or None
            The 5' overhang sequence if it exists, otherwise None.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.
        """
        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")
        
        # Check if a 5' variation exists
        five_var = self._five_prime_variation()

        # Check if a 5' overhang exists.
        if self.seq_mirna_start > 0 and self.seq_mirna_ref_start == 0:
            five_overhang_seq = self.seq[:self.seq_mirna_start]
            return five_overhang_seq
        elif five_var:
            if len(five_var[0]) > len(five_var[1]):
                five_overhang_seq = five_var[0][:-len(five_var[1])]
                return five_overhang_seq
        return None

    def five_prime_deletion(self):
        """
        This method retrieves the 5' deleted sequence from the reference 
        miRNA sequence (`mirna_seq`) if it exists, based on the alignment 
        with the object's sequence (`seq`). 

        Parameters
        ----------
        None

        Returns
        -------
        five_deletion_seq : str or None
            The 5' deleted sequence if it exists, otherwise None.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.
        """

        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")

        # Check if a 5' variation exists
        five_var = self._five_prime_variation()
        
        # Check if a 5' deletion exists.
        if self.seq_mirna_start == 0 and self.seq_mirna_ref_start > 0:
            five_deletion_seq = self.mirna_seq[:self.seq_mirna_ref_start]
            return five_deletion_seq
        elif five_var:
            if len(five_var[0]) < len(five_var[1]):
                five_overhang_seq = five_var[1][:-len(five_var[0])]
                return five_overhang_seq
        return None

    def three_prime_overhang(self):
        """
        This method retrieves the 3' overhang sequence from the object's
        sequence (`seq`) if it exists, based on the alignment with the 
        reference miRNA sequence (`mirna_seq`).

        Parameters
        ----------
        None

        Returns
        -------
        three_overhang_seq : str or None
            The 3' overhang sequence if it exists, otherwise None.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.
        """

        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")
        
        # Check if a 3' variation exists
        three_var = self._three_prime_variation()

        # Execute the requiered methods
        seq_wo_five_prime_overhang, can_wo_five_prime_overhang = self._remove_5prime_overhang()
        seq_aligned_clean, _ = self._get_clean_aligned_regions()
        # Check if a 3' overghang exists.
        if len(seq_wo_five_prime_overhang) > len(seq_aligned_clean) and len(seq_aligned_clean) == len(can_wo_five_prime_overhang):
            three_overhang_seq = seq_wo_five_prime_overhang[len(seq_aligned_clean):]
            return three_overhang_seq
        elif three_var:
            if len(three_var[0]) > len(three_var[1]):
                three_overhang_seq = three_var[0][:-len(three_var[1])]
                return three_overhang_seq
        return None

    def three_prime_deletion(self):
        """
        This method retrieves the 3' deleted sequence from the reference 
        miRNA sequence (`mirna_seq`) if it exists, based on the alignment 
        with the object's sequence (`seq`).

        Parameters
        ----------
        None

        Returns
        -------
        three_deletion_seq : str or None
            The 3' deleted sequence if it exists, otherwise None.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.
        """

        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")

        # Check if a 3' variation exists
        three_var = self._three_prime_variation()

        # Execute the required methods
        seq_wo_five_prime_overhang, can_wo_five_prime_overhang = self._remove_5prime_overhang()
        seq_aligned_clean, can_aligned_clean = self._get_clean_aligned_regions()
        # Check if a 3' deletion exists.
        if len(seq_wo_five_prime_overhang) == len(seq_aligned_clean) and len(can_aligned_clean) < len(can_wo_five_prime_overhang):
            three_deletion_seq = can_wo_five_prime_overhang[len(can_aligned_clean):]
            return three_deletion_seq
        elif three_var:
            if len(three_var[0]) < len(three_var[1]):
                three_overhang_seq = three_var[1][:-len(three_var[0])]
                return three_overhang_seq
        return None

    def substitutions(self):
        """
        This method retrieves the polymorphic changes (substitutions) 
        between the object's sequence (`seq`) and the reference miRNA 
        sequence (`mirna_seq`). It compares the aligned sequences and 
        identifies the positions where nucleotide differences occur.

        Parameters
        ----------
        None

        Returns
        -------
        substitutions : list of tuples
            A list of tuples representing the polymorphic changes. Each 
            tuple contains:
                - seq_nuc: The nucleotide in the query sequence.
                - query_position: The 1-based position in the query sequence.
                - can_nuc: The nucleotide in the canonical reference sequence.
                - ref_position: The 1-based position in the reference sequence.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.
        """

        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")

        # Empty variables
        five_seq = ''
        three_seq = ''
        five_can = ''
        three_can = ''

        # Ensure aligned sequences are obtained
        seq_aligned, can_aligned = self._get_aligned_regions()
        
        # Get the alignment start position for each sequence
        qseq_start = self.seq_mirna_start
        cseq_start = self.seq_mirna_ref_start

        # Get the 5' and 3' variation
        five_var = '' if self._five_prime_variation() is None else self._five_prime_variation()
        three_var = '' if self._three_prime_variation() is None else self._three_prime_variation()

        # This is done because if there is variation, there may be trimming or
        # tailing..
        if five_var:
            # Get the changed nucleotides in 5' (only those nucleotides changed)
            min_5 = len(min(five_var, key=len))
            five_seq = five_var[0][:min_5]
            five_can = five_var[1][:min_5]

            # Modify the start position
            qseq_start -= min_5
            cseq_start -= min_5

        # The same for 3' modifications...
        if three_var:
            min_3 = len(min(three_var, key=len))
            three_seq = three_var[0][:min_3]
            three_can = three_var[1][:min_3]

        # Recreate the alignment
        query_alig_seq = f'{five_seq}{seq_aligned}{three_seq}'
        can_alig_seq = f'{five_can}{can_aligned}{three_can}'

        # Get the substitutions
        substitutions = []
        for i in range(len(can_alig_seq)):

            # Get the nucleotide in "i" position
            seq_nuc = query_alig_seq[i]
            can_nuc = can_alig_seq[i]

            # If there's a polymorphic difference
            if seq_nuc != can_nuc:
                real_query_nuc_position = i + qseq_start
                real_ref_nuc_position = i + cseq_start
                substitutions.append((seq_nuc, int(real_query_nuc_position + 1),
                                            can_nuc, int(real_ref_nuc_position + 1)))
        return substitutions
    
    def is_template(self):
        """
        This method checks if the object's sequence (`seq`) is present 
        in the precursor sequence.

        Parameters
        ----------
        None

        Returns
        -------
        bool
            True if the sequence is found in the precursor sequence, 
            otherwise False.
        """
        return self.seq in self.precursor_seq
    
    def get_variant(self):
        """
        This method retrieves the different types of variants present in 
        the alignment between the object's sequence (`seq`) and the 
        reference miRNA sequence (`mirna_seq`). It identifies 5' and 3' 
        overhangs, deletions, and single nucleotide variants (SNVs), 
        categorizing them based on their templated or non-templated nature. 
        The variants are labeled according to the specific type of 
        modification and the position in the sequence.

        Variant types:
        --------------
        - iso_5p:+/-N : The 5' end shift. '+' indicates a shift to the right 
        (towards the 3' direction), and '-' indicates a shift to the left 
        (towards the 5' direction). N represents the number of nucleotides 
        difference.
        - iso_3p:+/-N : Same as iso_5p but for the 3' end.
        - iso_add5p:N : Non-template nucleotides added at the 5' end (N is the 
        number of nucleotides).
        - iso_add3p:N : Non-template nucleotides added at the 3' end (N is the 
        number of nucleotides).
        - iso_snv_seed : Single nucleotide variants (SNVs) affecting nucleotides 
        between positions 2-7.
        - iso_snv_central_offset : SNVs affecting nucleotide position 8.
        - iso_snv_central : SNVs affecting nucleotides between positions 9-12.
        - iso_snv_central_supp : SNVs affecting nucleotides between positions 
        13-17.
        - iso_snv : SNVs affecting nucleotides outside the above-mentioned 
        ranges.

        Change descriptions:
        --------------------
        - Additions are represented in capital letters (e.g., 'GTC').
        - Deletions are represented in lowercase letters (e.g., 'tt').
        Example: Changes iso_3p:TT,iso_add3p:GTC would indicate a 
        variant of iso_3p:+2 and iso_add3p:3.

        Parameters
        ----------
        None

        Returns
        -------
        variant_list : str
            A comma-separated string representing the types of variants 
            (e.g., iso_5p:+X,iso_add3p:Y) identified between the sequence 
            and the reference miRNA.

        variant_changes_list : str
            A comma-separated string describing the specific nucleotide 
            changes or additions for each variant (e.g., iso_5p:X,iso_add3p:Y).

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.

        Attributes
        ----------
        one_alignment_seq : object
            Stores the result of the alignment between the sequence and 
            the reference miRNA.

        seq_pre_ref_start : int
            The start position of the precursor reference sequence for 
            overhang calculation.
        """


        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")
        
        # Align the canonical miRNA to the precursor
        self._align_miRNA_to_precursor()

        # Get the modifications
        five_over = self.five_prime_overhang()
        five_del = self.five_prime_deletion()
        three_over = self.three_prime_overhang()
        three_del = self.three_prime_deletion()
        polymorphic = self.substitutions()

        # Variants lists
        variant_list = []
        variant_changes_list = []

        # 5' Overhang
        ########################################################################
        if five_over:
            # Get the overhang sequence in the precursor
            overhan_pre_start_5 = self.seq_pre_ref_start - len(five_over)
            overhan_pre_end_5 = self.seq_pre_ref_start
            pre_over_5 = self.precursor_seq[overhan_pre_start_5:overhan_pre_end_5]
            # Check if the overhang sequence is templated or non-templated
            if five_over == pre_over_5:
                # Templated
                variant_list.append(f'iso_5p:+{len(five_over)}')
                variant_changes_list.append(f'iso_5p:{five_over}')
            else:
                # Non-templated
                variant_list.append(f'iso_add5p:{len(five_over)}')
                variant_changes_list.append(f'iso_add5p:{five_over}')

        # 5' deletion
        ########################################################################
        if five_del:
            variant_list.append(f'iso_5p:-{len(five_del)}')
            variant_changes_list.append(f'iso_5p:{five_del.lower()}')

        # 3' Overhang
        ########################################################################
        if three_over:
            # Get the overhang sequence in the precursor
            overhang_pre_start_3 = self.seq_pre_ref_start + len(self.mirna_seq)
            overhang_pre_end_3 = self.seq_pre_ref_start + len(self.mirna_seq) + len(three_over)
            pre_over_3 = self.precursor_seq[overhang_pre_start_3:overhang_pre_end_3]
            # Check if the overhang sequence is templated or non-templated
            if three_over == pre_over_3:
                # Templated
                variant_list.append(f'iso_3p:+{len(three_over)}')
                variant_changes_list.append(f'iso_3p:{three_over}')
            else:
                # Non-templated
                variant_list.append(f'iso_add3p:{len(three_over)}')
                variant_changes_list.append(f'iso_add3p:{three_over}')

        # 3' deletion
        ########################################################################
        if three_del:
            variant_list.append(f'iso_3p:-{len(three_del)}')
            variant_changes_list.append(f'iso_3p:{three_del.lower()}')
        
        # SNV Variants
        ########################################################################
        if polymorphic:
            # Iterate through polymorphic changes:
            for change in polymorphic:
                # Get the position changed
                canonical_pos = change[3]
                # Within the seed
                if canonical_pos > 1 and canonical_pos < 8:
                    variant_list.append(f'iso_snv_seed')
                    variant_changes_list.append(f'iso_snv_seed:{change[2]}{change[1]}{change[0]}')
                # Central offset
                elif canonical_pos == 8:
                    variant_list.append(f'iso_snv_central_offset')
                    variant_changes_list.append(f'iso_snv_central_offset:{change[2]}{change[1]}{change[0]}')
                # Central
                elif canonical_pos > 8 and canonical_pos <= 13:
                    variant_list.append(f'iso_snv_central')
                    variant_changes_list.append(f'iso_snv_central:{change[2]}{change[1]}{change[0]}')
                # Central Supp
                elif canonical_pos > 13 and canonical_pos <= 17:
                    variant_list.append(f'iso_snv_central_supp')
                    variant_changes_list.append(f'iso_snv_central_supp:{change[2]}{change[1]}{change[0]}')
                else:
                    variant_list.append(f'iso_snv')
                    variant_changes_list.append(f'iso_snv:{change[2]}{change[1]}{change[0]}')
        
        # Remove duplicates from variant_list
        variant_list = list(set(variant_list))

        # Sort both lists
        variant_list.sort()
        variant_changes_list.sort()
        
        return (",".join(variant_list), ",".join(variant_changes_list))


    def get_classification(self):
        """
        This method classifies the sequence as either an isomiR or a 
        reference miRNA (ref_miRNA). It gathers information about the 
        variations present in the sequence, including 5' and 3' overhangs, 
        deletions, and single nucleotide variants (SNVs). It returns a 
        dictionary containing the classification details along with the 
        identified variants and changes.

        Classification types:
        ---------------------
        - 'isomiR' : If any variations (e.g., overhangs, deletions, SNVs) 
        are detected in the sequence.
        - 'ref_miRNA' : If the sequence matches the reference miRNA without 
        significant variation.

        Parameters
        ----------
        None

        Returns
        -------
        dict
            A dictionary containing classification details, including:
            - "id": A unique identifier for the isomiR (generated by 
            `get_isomir_id()`).
            - "sequence": The input sequence.
            - "variant": A comma-separated string representing the identified 
            variants.
            - "changes": A comma-separated string describing the specific 
            nucleotide changes.
            - "type": Either 'isomiR' or 'ref_miRNA', indicating the type of 
            the sequence.
            - "five_prime_overhang": The 5' overhang sequence, if present.
            - "five_prime_deletion": The 5' deletion sequence, if present.
            - "three_prime_overhang": The 3' overhang sequence, if present.
            - "three_prime_deletion": The 3' deletion sequence, if present.
            - "substitutions": A string representing polymorphic substitutions, 
            if present.
            - "is_template": Boolean indicating whether the sequence is part of 
            the precursor sequence.
            - "ref_miRNA": The reference miRNA sequence.
            - "precursor": The precursor sequence.

        Raises
        ------
        RuntimeError
            If the alignment has not been executed yet, a RuntimeError 
            is raised, instructing the user to call `align_to_miRNA()` first.

        Attributes
        ----------
        one_alignment_seq : object
            Stores the result of the alignment between the sequence and 
            the reference miRNA.

        seq : str
            The sequence to be classified.

        mirna_seq : str
            The reference miRNA sequence.

        precursor_seq : str
            The precursor sequence for the alignment.
        """


        # Chek if the alignment has been previously executed 
        if self.one_alignment_seq is None:
            raise RuntimeError("Alignment not yet run. Call `align_to_miRNA()` first.")

        # Execute variations related methods
        five_over = self.five_prime_overhang()
        five_del = self.five_prime_deletion()
        three_over = self.three_prime_overhang()
        three_del = self.three_prime_deletion()
        polymorphic = self.substitutions()
        variant = self.get_variant()

        # Create the variant and changes strings
        variant_str = None if variant[0] == '' else variant[0]
        changes_str = None if variant[1] == '' else variant[1]

        # Create a string for each polymorphic change ('T', 3, 'A', 3) -> A3T
        polymorphic_str = [f'{change[2]}{change[3]}{change[0]}' for change in polymorphic]
        polymorphic_str = None if len(polymorphic_str) == 0 else polymorphic_str
        polymorphic_str = ",".join(polymorphic_str) if polymorphic_str != None else None

        # Check if the sequence is isomiR
        is_isomir = any(var is not None for var in [five_over, five_del, three_over, three_del, polymorphic_str])
        if is_isomir:
            type = 'isomiR'
        else:
            type = 'ref_miRNA'

        return {
            "id": self.get_isomir_id(),
            "sequence": self.seq,
            "variant": variant_str,
            "changes": changes_str,
            "type": type,
            "five_prime_overhang": five_over,
            "five_prime_deletion": five_del,
            "three_prime_overhang": three_over,
            "three_prime_deletion": three_del,
            "substitutions": polymorphic_str,
            "is_template": self.is_template(),
            "ref_miRNA": self.mirna_seq,
            "precursor": self.precursor_seq
        }
    
    def _create_aligner(self):
        """
        Creates and configures a local pairwise aligner using the 
        'NUC.4.4' substitution matrix, and specific gap penalties for 
        alignment.

        This method is used to create an aligner object that will be used 
        for local sequence alignments. The aligner is configured to use 
        the specified substitution matrix and gap penalties.

        Parameters
        ----------
        None

        Returns
        -------
        aligner : PairwiseAligner
            A PairwiseAligner object configured with local alignment mode, 
            the 'NUC.4.4' substitution matrix, and gap penalties.
        """
        aligner = PairwiseAligner()
        aligner.mode = 'local'
        aligner.substitution_matrix = substitution_matrices.load("NUC.4.4")
        aligner.open_gap_score = -8
        aligner.extend_gap_score = -8
        return aligner
    
    def _swalignment(self, seq1, seq2):
        """
        Aligns two sequences using the Smith-Waterman algorithm for local 
        alignment and returns the alignment with the highest score.

        Parameters
        ----------
        seq1 : str
            The first sequence to align.
        
        seq2 : str
            The second sequence to align.

        Returns
        -------
        alignment : Alignment
            The alignment object that has the highest score from the 
            Smith-Waterman local sequence alignment.
        """
        aligner = self._create_aligner()
        alignments = aligner.align(seq1, seq2)
        return max(alignments, key=lambda aln: aln.score)

    def _align_miRNA_to_precursor(self):
        """
        Aligns the reference miRNA sequence to the precursor sequence using the 
        Smith-Waterman algorithm and stores the alignment start positions.

        Parameters
        ----------
        None

        Returns
        -------
        alignment : Alignment
            The alignment object resulting from the Smith-Waterman alignment 
            of the miRNA sequence to the precursor sequence.
        """
        alignment = self._swalignment(self.mirna_seq, self.precursor_seq)
        self.seq_pre_start = alignment.aligned[0][0][0]
        self.seq_pre_ref_start = alignment.aligned[1][0][0]
        return alignment

    def _five_prime_variation(self):
        """
        Gets the 5' variation between the sequence and the reference miRNA 
        sequence if it exists.

        Parameters
        ----------
        None

        Returns
        -------
        tuple or None
            A tuple containing two strings: the variation at the 5' end 
            of the sequence and the reference miRNA sequence. Returns None if 
            no variation exists.
        """
        # Check if a 5' variation exists.
        if self.seq_mirna_start > 0 and self.seq_mirna_ref_start > 0:
            five_variation_seq = self.seq[:self.seq_mirna_start]
            five_variation_can = self.mirna_seq[:self.seq_mirna_ref_start]
            return (five_variation_seq, five_variation_can)
        return None

    def _three_prime_variation(self):
        """
        Gets the 3' variation between the sequence and the reference miRNA 
        sequence if it exists.

        Parameters
        ----------
        None

        Returns
        -------
        tuple or None
            A tuple containing two strings: the variation at the 3' end 
            of the sequence and the reference miRNA sequence. Returns None if 
            no variation exists.
        """
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
    
    def _remove_5prime_overhang(self):
        """
        This method removes the 5' overhang from both the sequence and 
        the canonical miRNA by slicing the sequences starting from their 
        respective alignment start positions.

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            A tuple containing two strings: the sequence and the reference 
            miRNA sequence with the 5' overhang removed.
        """
        seq_wo_five_prime_overhang = self.seq[self.seq_mirna_start:]
        can_wo_five_prime_overhang = self.mirna_seq[self.seq_mirna_ref_start:]
        return (seq_wo_five_prime_overhang, can_wo_five_prime_overhang)
    
    def _get_aligned_regions(self):
        """
        This method retrieves the aligned regions from the stored alignment 
        result, where the first element is the aligned query sequence and 
        the second element is the aligned reference miRNA sequence.

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            A tuple containing two strings: the aligned query sequence and 
            the aligned reference miRNA sequence.
        """
        seq_aligned = self.one_alignment_seq[0]
        can_aligned = self.one_alignment_seq[1]
        return (seq_aligned, can_aligned)

    def _get_clean_aligned_regions(self):
        """
        This method retrieves the aligned regions from the stored alignment 
        result and removes any gaps (dashes) from both the query and 
        reference sequences.

        Parameters
        ----------
        None

        Returns
        -------
        tuple
            A tuple containing two strings: the cleaned aligned query sequence 
            and the cleaned aligned reference miRNA sequence with gaps removed.
        """
        seq_aligned_clean = self.one_alignment_seq[0].replace('-', '')
        can_aligned_clean = self.one_alignment_seq[1].replace('-', '')
        return (seq_aligned_clean, can_aligned_clean)

## FUNCTIONS

def read_blastn_df(file: str) -> pd.DataFrame:
    """
    Reads a BLASTN result file containing alignment results between query
    sequences and a miRNA database, including both mature miRNAs and miRNA
    precursors. The function returns a pandas DataFrame with the appropriate
    column names and data.

    The input file is the result of performing BLASTN alignments where
    query sequences are compared against a reference database of mature miRNAs
    and precursor miRNAs. The file contains two sets of results:
    1. Alignments to mature miRNAs.
    2. Alignments to precursor miRNAs.
    
    This function reads the file, checks the number of columns, assigns the
    correct column  names based on the number of columns in the file, and
    returns a pandas DataFrame.

    Parameters
    ----------
    file : str
        Path to the BLASTN result file (in tab-separated format) to be read.
        The file should  contain alignment results for both mature miRNAs and
        miRNA precursors.

    Returns
    -------
    pandas.DataFrame
        A DataFrame containing the BLASTN results with the appropriate column
        names.
    """

    # Read the input file
    df = pd.read_csv(file, sep='\t', header=None)
    
    # Check the number of columns of the input dataframe
    num_col = df.shape[1]

    # List with the dataframe colnames
    colnames = ['Query_id', 'miRNA_id', 'miRNA_group_id', 'Identity_m',
                'Length_m', 'Mismatches_m', 'Gap_Openings_m', 'Q_Start_m',
                'Q_End_m', 'S_Start_m', 'S_End_m', 'E_value_m', 'Bit_Score_m',
                'Query_seq', 'miRNA_seq', 'Precursor_id', 'Precursor_group_id',
                'Identity_p', 'Length_p', 'Mismatches_p', 'Gap_Openings_p',
                'Q_Start_p', 'Q_End_p', 'S_Start_p', 'S_End_p', 'E_value_p',
                'Bit_Score_p', 'Precursor_seq']
    
    # The df has raw counts columns
    if num_col == 30:
        colnames.extend(['Raw_counts', 'Raw_counts_filt'])
    # The df has RPM columns
    if num_col == 32:
        colnames.extend(['Raw_counts', 'Raw_counts_filt', 'RPM', 'RPM_filt'])
    # Add colnames
    df.columns = colnames

    return df


def filter_isomirs(classification, substitutions_inside=1, nt_5add=0, nt_3add=3, diff_ends=4):
    """
    Filters isomiRs based on allowed sequence variations with respect to a 
    reference miRNA.

    This function evaluates whether an isomiR meets a set of criteria to be
    considered  valid based on its sequence variations compared to a reference
    miRNA. It checks for substitutions within the sequence, non-templated
    additions at the 5' and 3' ends, and the number of modifications at the
    terminal positions. The filtering is intended to remove isomiRs with
    excessive or undesired variation.

    Parameters
    ----------
    classification : dict
        A dictionary representing the classification of a single isomiR.
        It must contain at least the following keys:
        - 'substitutions': a comma-separated string of nucleotide substitutions 
          (e.g., "A2G,C12U") or 'None'.
        - 'ref_miRNA': the reference mature miRNA sequence.
        - 'changes': a comma-separated list of change annotations 
          (e.g., "iso_add3p:U,iso_snv:5:G>A").

    substitutions_inside : int, optional
        Maximum number of allowed internal substitutions (not located at the
        sequence ends).
        Default is 1.

    nt_5add : int, optional
        Maximum number of allowed non-templated nucleotides added at the 5' end.
        Default is 0.

    nt_3add : int, optional
        Maximum number of allowed non-templated nucleotides added at the 3' end.
        Default is 3.

    diff_ends : int, optional
        Maximum number of total modifications (substitutions) allowed across
        both ends (5' + 3').
        Default is 4.

    Returns
    -------
    bool
        True if the isomiR passes all the filtering criteria; False otherwise.
    """

    # Get the required information
    substitutions_str = None if classification['substitutions'] == 'None' else classification['substitutions']
    sequence_length = len(classification['ref_miRNA'])
    changes = classification['changes']

    ## 1. 5' or 3' ends modifications (not addition or deletion)
    ############################################################################

    # Booleans
    subs_inside_validity = True
    nt_5add_validity = True
    nt_3add_validity = True
    diff_ends_validity = True
    isomir_valid = True


    # Check if there are substitutions within the sequence
    if substitutions_str:

        # Check whether there are insertions or deletions within the sequence (not allowed)
        if "-" not in substitutions_str:

            # Convert the substitutions into a list of numerical positions
            positions = [int(sub[1:-1]) for sub in substitutions_str.split(',')]
            
            # Sort the positions to check if they are consecutive
            positions.sort()

            # Check if there is a modification in any of the ends
            mod_5 = True if 1 in positions else False
            mod_3 = True if sequence_length in positions else False

            # Counter
            subs_5 = 1
            subs_3 = 1

            # If there is a modification at 5'
            ends_modifications = []
            if mod_5:
                for i in range(1, len(positions)):
                    if positions[i] == positions[i-1] + 1:
                        subs_5 += 1
                    else:
                        break
                    
            # If there is a modification at 3'
            if mod_3:
                for i in range(len(positions)-2, -1, -1):
                    if positions[i] == positions[i+1] - 1:
                        subs_3 += 1
                    else:
                        break
                    
            # Remove from positions list those that have been classified as end modifications
            positions_inside = [pos for pos in positions if pos not in ends_modifications]

            # If there is more than one modification inside the sequence...
            if len(positions_inside) > substitutions_inside:
                subs_inside_validity = False

            ## 2. Polymorphic changes within the sequence
            ############################################################################
            # Check the diff_ends condition
            if (subs_5 + subs_3) > diff_ends:
                diff_ends_validity = False
        else:
            isomir_valid = False

    ## 3. Non-templated additions
    ############################################################################

    if changes:

        # Split the string
        changes_list = changes.split(',')
        
        # Iterate through the changes
        for change in changes_list:

            # Get the non-templated additions
            nt_5_add = change if 'iso_add5p' in change else None
            nt_3_add = change if 'iso_add3p' in change else None
            
            # 5' non-templated addition
            if nt_5_add:
                num_nuc_5_nt = len(nt_5_add.split(":")[1])
                if num_nuc_5_nt > nt_5add:
                    nt_5add_validity = False
            # 3' non-templated addition
            if nt_3_add:
                num_nuc_3_nt = len(nt_3_add.split(":")[1])
                if num_nuc_3_nt > nt_3add:
                    nt_3add_validity = False

    ## 3. Check if the sequence meets all selection criteria.
    ############################################################################

    if not all([subs_inside_validity, nt_5add_validity, nt_3add_validity, diff_ends_validity]):
        isomir_valid = False
        
    return isomir_valid


def isomir_group_name(name):
    """
    Generates a unique isomiR name based on the canonical miRNA variants it is
    associated with. This method processes the canonical miRNA name, which may
    contain multiple variants separated by "|", and extracts relevant
    information to create a unified isomiR name. Since different miRNA variants
    from distinct precursors can share the same sequence, this method ensures a
    comprehensive representation.

    Example of a canonical miRNA name:
        "ath-miR156a-5p|ath-miR156b-5p|ath-miR156c-5p|ath-miR156d-5p"
    
    This method processes such names and generates a standardized identifier
    for the isomiR.

    Example output for the above canonical miRNA name:
        "miR156abcd-5p"
    
    Where:
    - "miR156abcd-5p" represents the unified name including all miRNA variants.

    Returns:
        str: A unique isomiR name including all associated miRNA variants.
    """

    # Split the string by the delimiter "|"
    miRNA_vars = name.split('|')
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
    isomir_name = f'{"|".join(miRNA_class_names_list)}'

    return isomir_name


def get_mirna_family(annotation: str) -> str:
    """
    This function extracts the miRNA family name from a given miRNA annotation
    string. The annotation string is expected to follow a naming convention
    where the family name is composed of a precursor identifier (e.g., "miR")
    followed by a number and, optionally, a letter (e.g., "156a" or "156").

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


def simplify_variant(variant_str: str) -> set:
    """
    This fucntion simplifies the variant string by extracting only the variant
    type, removing any suffixes after the colon.

    Parameters
    ----------
    variant_str : str or None
        A string representing variant annotations, separated by commas.
        Each annotation follows the format: <variant_type>:<details>,
        e.g., "iso_3p:-1,iso_5p:+2".

    Returns
    -------
    set
        A set of simplified variant types (e.g., {'iso_3p', 'iso_5p'}).
        Returns an empty set if the input is NaN.
    """
    if pd.isna(variant_str):
        return set()
    variant_items = variant_str.split(',')
    simplified = set(item.split(':')[0] for item in variant_items)
    return simplified


def assign_sequence_class(df:pd.DataFrame) -> pd.DataFrame:
    """
    This function assigns a classification to each unique sequence based on
    their type and variant information from the input DataFrame. Sequences with
    'RPM_filt' value other than 'PASS' will be marked as None.

    Classification Rules:
    - If all entries of a sequence are 'ref_miRNA': assign 'ref_miRNA'.
    - If all entries are 'isomiR':
        - If all have the same single variant type: assign that type.
        - If all have consistent multiple variant types: assign 'Mixed'.
        - If there is inconsistency across variant types: assign 'Undefined'.
    - If mixed types ('ref_miRNA' and 'isomiR') are found: assign 'ref_miRNA'.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing at least the columns 'sequence', 'type',
        'variant', and 'RPM_filt'.

    Returns
    -------
    pandas.DataFrame
        A copy of the input DataFrame with an added column 'Class',
        representing the classification for each sequence.
    """
    # Check if the dataframe must be filtered before the class asignment
    if 'RPM_filt' in df.columns:
        if 'Raw_counts_filt' in df.columns:
            df_pass = df[df['Raw_counts_filt'] == 'PASS']
        else:
            df_pass = df[df['RPM_filt'] == 'PASS']
    else:
        df_pass = df

    # Prepare a dict to store results
    assignment = {}

    # Group by sequence
    for seq, group in df_pass.groupby('sequence'):
        
        # Get the sequence type
        types = set(group['type'])

        # If the sequence is a canonical miRNA...
        if types == {'ref_miRNA'}:
            assignment[seq] = 'ref_miRNA'
        # If the sequence is an isomiR...
        elif types == {'isomiR'}:
            
            # Get the variants
            simplified_variants = group['variant'].apply(simplify_variant)
            
            # Convert list of sets to a list of sorted tuples for consistent comparison
            sorted_variants = [tuple(sorted(variants)) for variants in simplified_variants]

            # Check if all the types of variants are the same that the first one
            if all(variants == sorted_variants[0] for variants in sorted_variants):
                # Get the variants present in the sequence
                unique_variants = set(sorted_variants[0])
                # Only one variant type
                if len(unique_variants) == 1:
                    assignment[seq] = sorted_variants[0][0]
                # Multiple variant types but consistent across all
                else:
                    # Check for the special case: only iso_3p and iso_5p
                    if unique_variants == {"iso_3p", "iso_5p"}:
                        assignment[seq] = "mixed_shift"
                    else:
                        assignment[seq] = "mixed"
            else:
                assignment[seq] = 'undefined'
        # If it is both canonical and isomiR (it must not)...
        else:
            assignment[seq] = 'ref_miRNA'

    # Return the final df
    df = df.copy()
    df['Class'] = df['sequence'].map(assignment)
    return df

    
def build_attributes(row, hit_counts):
    """
    Constructs the 'attributes' field for a GFF3 entry specific to isomiR
    annotations.

    This function gathers relevant information from a row of the isomiR
    classification DataFrame, including read sequence, isomiR ID, reference
    miRNA data, variant details, expression values, and filtering status.

    Parameters
    ----------
    row : pandas.Series
        A single row from the isomiR classification DataFrame.

    hit_counts : dict
        A dictionary mapping each isomiR ID (`id`) to the number of unique 
        `miRNA_group_id` values it is associated with.

    Returns
    -------
    str
        A semicolon-separated string formatted for the GFF3 'attributes' field,
        representing isomiR annotation metadata.
    """

    # Default variables
    rpm='None'
    counts='None'
    filt='None'

    # Check if Raw_counts column exists in df
    if 'Raw_counts' in row:
        counts=row['Raw_counts']
        if row['Raw_counts_filt'] == 'REJECT':
            filt=f"{row['Raw_counts_filt']}:lowRawCounts"
        else:
            filt = row['Raw_counts_filt']

    # Check if RPM column exists in df
    if 'RPM' in row:
        rpm=row['RPM']
        if row['RPM_filt'] == 'REJECT':
            filt = f"{row['RPM_filt']}:lowRPM"
        else:
            filt = row['RPM_filt']
        
    # Create the attributes section
    parts = [
        f"Read={row['sequence']}",
        f"UID={row['id']}",
        f"Name={row['ref_miRNA_name']}",
        f"Parent={row['ref_precursor_name']}",
        f"Variant={row['variant']}",
        f"Changes={row['changes']}",
        f"Cigar=None",
        f"Hits={hit_counts.get(row['id'], 0)}",
        f"Expression={counts}",
        f"Norm={rpm}",
        f"Filter={filt}",
        f"miRNA_fam={row['miRNA_family']}",
        f"Class={row['Class']}",
        f"miRNA_seq={row['ref_miRNA']}",
        f"Parent_seq={row['precursor']}",
    ]
    return "; ".join(parts)


def create_isomirs_gff3(df, sample, database, path_out):
    """
    Generates a GFF3 file containing isomiR annotations based on a classification
    DataFrame.

    This function transforms a DataFrame of isomiR classification results into
    GFF3 format. It constructs a mirGFF3-compliant file following the specifications
    defined at: https://github.com/miRTop/mirGFF3/blob/master/definition.md.

    Additionally, several custom fields are added to the attributes column, which
    are not part of the original mirGFF3 specification:

    - Norm:         Normalized expression values (e.g., RPM), complementing raw
                    counts in 'Expression'.
    - miRNA_seq:    The mature miRNA sequence to which the query aligns.
    - Parent_seq:   The precursor miRNA sequence associated with the canonical miRNA.
    - miRNA_fam:    The reference miRNA family name.
    - Class:        The classification of the sequence based on its variants

    These enhancements support downstream analysis and improve visualization and interpretation
    of isomiR data.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing annotated isomiR information, with required fields such as
        'ref_miRNA_name', 'type', 'S_Start_p', 'S_End_p', 'Bit_Score_p', etc.
    sample : str
        Identifier of the sample being processed, included in the GFF3 header metadata.
    database : str
        Source database of the annotations (e.g., 'miRBase').
    path_out : str
        Output file path for saving the generated GFF3 file.

    Returns
    -------
    None
    """

    # Calculate the number of unique miRNA_group_id per id
    hit_counts = df.groupby('id')['miRNA_group_id'].nunique()

    # Build GFF3 DataFrame
    gff3_df = pd.DataFrame({
        'seqid': df['ref_miRNA_name'],
        'source': database,
        'type': df['type'],
        'start': df['S_Start_p'],
        'end': df['S_End_p'],
        'score': df['Bit_Score_p'],
        'strand': '+',
        'phase': '.',
        'attributes': df.apply(lambda row: build_attributes(row, hit_counts), axis=1)
    })

    # Write gff3 file
    date = datetime.now().strftime("(%d-%m-%Y)")
    with open(path_out, "w") as f:
        f.write("## mirGFF3. VERSION 1.2 modified\n")
        f.write(f"## source-ontology: {database} {date}\n")
        f.write("## TOOLS: Blastn\n")
        f.write(f"## COLDATA: {sample}\n")
        for _, row in gff3_df.iterrows():
            f.write("\t".join(map(str, row)) + "\n")


def create_summary_df(df:pd.DataFrame) -> pd.DataFrame:
    """
    This function generates a one-row summary DataFrame showing the count of
    each isomiR class from the input DataFrame, optionally filtering based on
    quality flags.

    This function counts the number of occurrences of each classification type
    found in the 'Class' column of the input DataFrame. If the DataFrame contains
    either 'RPM_filt' or 'Raw_counts_filt', only rows marked as 'PASS' will be 
    considered for the summary.

    Expected values in the 'Class' column include:
    - ref_miRNA
    - iso_5p, iso_3p, iso_add3p, iso_add5p
    - iso_snv_seed, iso_snv_central_offset, iso_snv_central, iso_snv_central_supp,
      iso_snv
    - Mixed
    - Mixed_shift
    - Undefined

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing at least the column 'Class', and optionally 
        'RPM_filt' or 'Raw_counts_filt' to determine which rows to include 
        based on filtering.

    Returns
    -------
    pandas.DataFrame
        A one-row DataFrame with each column representing a class type and 
        its corresponding count.
    """

    # Check if the dataframe must be filtered before the class asignment
    if 'Raw_counts_filt' in df.columns:
        if 'RPM_filt' in df.columns:
            df_pass = df[df['RPM_filt'] == 'PASS']
        else:
            df_pass = df[df['Raw_counts_filt'] == 'PASS']
    else:
        df_pass = df

    # Remove duplicate rows based on the 'sequence' column
    df_pass = df_pass.drop_duplicates(subset=['sequence'])

    # List of possible classes
    classes = [
        "ref_miRNA", "iso_5p", "iso_3p", "iso_add3p", "iso_add5p", "iso_snv_seed",
        "iso_snv_central_offset", "iso_snv_central", "iso_snv_central_supp",
        "iso_snv", "mixed", "mixed_shift", "undefined"
    ]

    # Count the number of appearances of each class
    class_counts = df_pass['Class'].value_counts()

    # Create the summary dataframe
    summary_data = {
        class_name: class_counts.get(class_name, 0) for class_name in classes
    }
    # Get the total isomirs value
    summary_data['num_isomirs'] = sum(summary_data.values())
    summary_df = pd.DataFrame([summary_data])

    # Sort the columns
    cols = ['num_isomirs'] + [c for c in summary_df.columns if c != 'num_isomirs']
    summary_df = summary_df[cols]

    return summary_df



## MAIN
def main():
    # '''
    # Main program
    # '''

    # Define the argument parser
    parser = argparse.ArgumentParser(prog='isomiRs_identification', 
                                    description='''This program performs the
                                    identification and classification of
                                    isomiRs using sequences from FASTA files
                                    or differential expression analysis
                                    results.''',  
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Argument for specifying the input file ID
    parser.add_argument('-i', '--id', type=str, required=True, 
                        help='File ID.')
    # Argument for specifying the name of the database used
    parser.add_argument('-d', '--database', type=str, required=True, 
                        help='Name of the database used for the annotation of' \
                        'the sequences.')
    # Argument for specifying the input file.
    parser.add_argument('-f', '--input', type=str, required=True, 
                        help='BLASTN result file in TSV format containing alignment' \
                        'results between query sequences and a miRNA database, including' \
                        'both mature miRNAs and miRNA precursors.')
    # Argument for specifying the number of allowed substitutions
    parser.add_argument('-s', '--substitutions', type=int, required=True, 
                        help='Maximum number of allowed internal substitutions'
                        '(not located at the sequence ends.')
    # Argument for specifying the number of allowed substitutions
    parser.add_argument('-b', '--five_add', type=int, required=True, 
                        help='Maximum number of allowed non-templated nucleotides' \
                        'added at the 5 end')
    # Argument for specifying the number of allowed substitutions
    parser.add_argument('-e', '--three_add', type=int, required=True, 
                        help='Maximum number of allowed non-templated nucleotides' \
                        'added at the 3 end.')
    parser.add_argument('-x', '--ends_modification', type=int, required=True, 
                        help='Maximum number of total modifications (substitutions)' \
                        'allowed across both ends (5 + 3).')
    # Parse the arguments
    args = parser.parse_args()

    ## 1. Check arguments
    #######################################################################

    try:
        id = args.id
        database = args.database
        input_file = args.input
        substitutions = args.substitutions
        nt_5add = args.five_add
        nt_3add = args.three_add
        diff_ends = args.ends_modification

    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()    

    ## 2. IsomiRs classification
    #######################################################################
    
    # Read the input file
    blastn_df = read_blastn_df(input_file)

    data = []
    # isomiR classification
    for row in blastn_df.itertuples(index=False):

        # Get the required sequences
        query_seq = row.Query_seq
        mirna_seq = row.miRNA_seq
        precursor_seq = row.Precursor_seq

        # Classification
        isomirs = IsomiRs(query_seq, mirna_seq, precursor_seq)
        isomirs.align_to_miRNA()
        classification_res = isomirs.get_classification()
        classification_res_clean = {k: v if v is not None else "None" for k, v in classification_res.items()}

        # Check if the isomiR es valid:
        is_valid = filter_isomirs(classification_res_clean, substitutions, nt_5add, nt_3add, diff_ends)

        # Save the results of valid isomiRs
        if is_valid:
            # Add the name of the canonical and precursor sequences
            classification_res_clean['ref_miRNA_name'] =  row.miRNA_id
            classification_res_clean['ref_miRNA_group'] = isomir_group_name(row.miRNA_group_id)
            classification_res_clean['ref_precursor_name'] = row.Precursor_group_id
            # Append the results to the classification_res_clean list
            data.append(classification_res_clean)

    # Save the results table
    isomir_classification_df = pd.DataFrame(data)

    ## 3 Get the miRNA family name for each isomiR
    ########################################################################

    # Create a list with the names of the canonical miRNAs
    canonical_list_names = isomir_classification_df['ref_miRNA_group'].tolist()
    
    # Iterate through canonical_list_names ['gma-miR167e|gma-miR167f',...]
    miRNA_families = []
    for canonical_names in canonical_list_names:
        # Split canonical names (e.g. 'gma-miR167e|gma-miR167f' -> ['gma-miR167e', 'gma-miR167f'])
        canonical_names_list = canonical_names.split("|")
        # Get unique sorted family IDsref_miRNA
        unique_family_ids = sorted({get_mirna_family(canonical_name) for canonical_name in canonical_names_list})
        # Construct final family name based on the length of unique family IDs
        final_family_name = "/".join(unique_family_ids) if len(unique_family_ids) > 1 else unique_family_ids[0]
        miRNA_families.append(final_family_name)

    # Add the miRNA family ids to the isomir classification dataframe
    isomir_classification_df['miRNA_family'] = miRNA_families

    ## 4. Combine isomiR classification results with other dataframes
    ########################################################################

    # Add the blastn results to the output dataframe
    isomir_class_with_blast_df = pd.merge(isomir_classification_df, blastn_df, 
                left_on=['sequence', 'ref_miRNA_name'], 
                right_on=['Query_seq', 'miRNA_id'], 
                how='left')

    # Remove undesired sequences
    isomir_class_with_blast_df = isomir_class_with_blast_df.drop(
        columns=['Query_seq', 'miRNA_seq', 'Precursor_id',
                'Precursor_seq', 'miRNA_id']
        )
    other_columns_to_drop = isomir_class_with_blast_df.filter(regex=r'^(Identity_|Length_|Mismatches_|Gap_Openings_|Q_Start_|Q_End_).{1}$').columns
    isomir_class_with_blast_df.drop(columns=other_columns_to_drop, inplace=True)

    ## 5. Create the output GFF3 file
    ########################################################################

    # Get the class of the sequences
    isomir_complete_df = assign_sequence_class(isomir_class_with_blast_df)
    
    # Create a gff3 file
    create_isomirs_gff3(isomir_complete_df, id, database, f'{id}.gff3')

    ## 6. Create a summary file
    ########################################################################
    summary_df = create_summary_df(isomir_complete_df)
    summary_df.to_csv(f'{id}.summary.tsv', sep='\t', index=False)

    # Remove __pycache__ directory
    os.system(f"rm -rf ../../../bin/__pycache__")

if __name__ == "__main__":
    main()