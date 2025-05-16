#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   08-miRNA_event_matrix_builder.py
#
#   This script generates two summary matrices (presence/absence and Shrunken
#   log2FoldChange) that describe the differential expression of miRNA families
#   across multiple stress conditions or experimental events, based on annotated
#   DESeq2 outputs and sample metadata files.
#
#   The process is structured as follows:
#   
#   1. Assign a unique identifier to each stress event.
#
#   For each annotation file provided (corresponding to a stress event), the script
#   reads its associated metadata and the list of samples used for differential
#   expression analysis. It then creates a unique hierarchical identifier based on
#   metadata fields such as species, treatment, and user-defined descriptors. These
#   identifiers are saved in a correspondence file (`id_correspondence.tsv`).
#
#   2. Extract the most representative miRNA family per stress event.
#
#   From each annotation file, the script selects the most representative member of
#   each miRNA family (based on maximum baseMean value) and retrieves its Shrunken
#   log2FoldChange. This results in a dictionary linking each original file/event to
#   a list of differentially expressed miRNA families and their log2fc values.
#
#   3. Rename the stress events using the unique identifiers.
#
#   Using the previously generated ID mapping, each original file identifier is
#   replaced with its corresponding hierarchical ID to ensure consistency across
#   the final matrices.
#
#   4. Generate two summary matrices.
#
#   Two final tables are constructed:
#   - `presence_absence_table.tsv`: Binary matrix indicating the presence (1) or
#     absence (0) of each miRNA family across events.
#   - `shrunken_log2fc_table.tsv`: Matrix containing the Shrunken log2FoldChange
#     of each miRNA family in each event. If a family is not present, the value
#     is marked as 'NA'.
#
#
#   Authors: Antonio Gonzalez Sanchez
#   Date: 16/05/2025
#   Version: 1.0
#
#******************************************************************************


## IMPORT MODULES

import argparse
import csv
from natsort import natsorted
import pandas as pd
import sys
import warnings
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


## FUNCTIONS

def get_events_ids(samples_list:str, metadata_list:int, samples_list_str:str, fields:list, save_id_file:bool = True):
    '''
    This function generates a unique numeric identifier for each specific stress
    event (i.e., each comparison performed in the differential expression
    analysis). The identifier is composed of integers, where each number
    represents a unique value found in one of the metadata fields specified in
    the `fields` argument, following the same order.

    For example, if the fields are ['Species', 'Treatment', 'Tissue'], an
    identifier such as "1.2.1.1" could mean:
    - Species = Arabidopsis thaliana → 1
    - Treatment = drought → 2
    - Tissue = seedling → 1
    - ".1" is added at the end to distinguish between comparisons with identical metadata.

    Additionally, the function can generate a TSV file (`id_correspondence.tsv`)
    that maps each metadata field-value combination to its assigned numeric ID.

    Parameters
    ----------
    samples_list : list of str
        List of annotation file paths, each representing a specific
        stress comparison in the differential expression analysis.

    metadata_list : list of str
        List of metadata file paths (TSV format), one per annotation file.
        These contain experimental information such as sample group, condition,
        species, etc.

    samples_list_str : list of str
        List of comma-separated sample names used for the comparison, each
        containing the samples identifiers used in the corresponding annotation
        file.

    fields : list of str
        List of metadata column names to be used when building the identifier.
        Each unique value in these columns is assigned a number.

    save_id_file : bool, optional
        If True (default), a TSV file will be saved documenting the correspondence
        between field values and their assigned numeric IDs.

    Returns
    -------
    dict
        Dictionary mapping each annotation file name (e.g., PRJNA277424_3_0) to
        its unique numeric identifier (e.g., '1.2.1.1'), based on the metadata
        content of the comparison.
    '''

    # Empty variables
    ids_dic = {}
    field_id_dic = {}
    ids_done = []

    # Create a dictionary to map field names to unique IDs
    for anno, meta, samp in zip(samples_list, metadata_list, samples_list_str):

        # Create a list with the samples of the current annotation file
        samp_list = samp.split(',')

        # Get the sample group from the annotation file name
        samples_group = int(anno.split('.')[0].split('_')[1])

        # Read the metadata file
        df_metadata = pd.read_csv(meta, sep='\t')

        # Select those rows belonging to the group of samples
        group_rows = df_metadata[
            (df_metadata['Group'] == samples_group) &
            (df_metadata['Run'].isin(samp_list))
        ]
        # Select those fields that are in the metadata file
        all_fields = [col for col in fields
                    if col.lower() in [c.lower() for c in group_rows.columns]]
        id_parts = []
        for col in all_fields:
            # Get the real column name (respecting case)
            real_col = next(c for c in group_rows.columns if c.lower() == col.lower())

            # Get unique values from that column (ignoring NaN)
            unique_vals = group_rows[real_col].dropna().unique()
            
            # If there is only one unique value, assign an ID
            if len(unique_vals) == 1:

                # Check if this id has been already assigned
                dic_id = f'{real_col}_{unique_vals[0]}'
                if dic_id not in field_id_dic:
                    # If the dictionary is empty, assign the ID
                    if not field_id_dic:
                        # Assign the ID
                        field_id_dic[dic_id] = 1
                        # Append the new id to the list of id parts
                        id_parts.append(1)
                    else:
                        # Get the previously assigned ids for the levels of this field
                        ids_asigned_for_this_field = [field_id_dic[key] for key in field_id_dic if key.startswith(real_col)]

                        # Check if other levels have been saved in the dictionary
                        if ids_asigned_for_this_field:
                            # Get the maximum id assigned for this field
                            last_num_id = max(ids_asigned_for_this_field)
                            # Assign the ID
                            field_id_dic[dic_id] = last_num_id + 1
                            # Append the new id to the list of id parts
                            id_parts.append(last_num_id + 1)
                        else:
                            # Assign the ID
                            field_id_dic[dic_id] = 1
                            # Append the new id to the list of id parts
                            id_parts.append(1)
                else:
                    # Append the new id to the list of id parts
                    id_parts.append(field_id_dic[dic_id])
            else:
                raise ValueError(f"Expected only one unique value for column '{real_col}', found: {unique_vals}")
            
        # Create the identifier
        id = '.'.join(map(str, id_parts))

        # Create the identifier for the comparison
        i = 1
        while True:
            new_id = f'{id}.{i}'
            if new_id not in ids_done:
                break
            i += 1

        # Save the identifier in the dictionary
        ids_dic['_'.join(anno.split('.')[0].split('_'))] = new_id
        # Save the identifier in the "done" list
        ids_done.append(new_id)
        
    # Convert to list of rows
    rows = []
    for key, value in field_id_dic.items():
        type_, name = key.split('_', 1)
        rows.append((type_, name, value))

    # Sort by Type and then ID
    rows.sort(key=lambda x: (x[0], x[2]))

    if save_id_file:
        # Write to TSV file
        with open('id_correspondence.tsv', 'w', newline='') as tsvfile:
            writer = csv.writer(tsvfile, delimiter='\t')
            writer.writerow(['Type', 'Name', 'ID'])
            writer.writerows(rows)
            return ids_dic
    
    return ids_dic


def get_miRNA_families(annot_files_list: list):
    '''
    This function processes a list of annotation files containing differentially
    expressed miRNAs for various stress events or comparisons. Each file is
    expected to contain miRNA family annotations along with expression
    statistics (e.g., baseMean, log2FoldChange, etc.).

    For each file, the function selects a single representative miRNA family
    per family group — specifically, the one with the highest `baseMean` value.
    It then extracts the `Shrunkenlog2FoldChange` associated with that miRNA
    family, assuming it reflects the expression change for the most
    representative member.

    The function returns two objects:
    - A dictionary mapping the experiment identifier (derived from the file name)
      to a list of tuples containing miRNA family names and their corresponding
      log2FoldChange values.
    - A sorted list of all unique miRNA families found across all input files.

    Parameters
    ----------
    annot_files_list : list of str
        List of file paths to the annotation tables (TSV format), each
        representing the results of differential expression analysis for a given
        stress event or comparison.

    Returns
    -------
    dict
        Dictionary where keys are experiment identifiers (extracted from the file
        name before the first dot), and values are lists of tuples. Each tuple
        contains:
            - miRNA family name (str)
            - Shrunken log2FoldChange (float) from the most representative member.
    
    list
        Sorted list of all unique miRNA family names found in the input files.
    '''

    # Iterate files
    exp_miRNAs = {}
    for file in annot_files_list:

        # Get the experiment name
        file_id = file.split('.')[0]

        # Read the file
        df = pd.read_csv(file, sep='\t')
        
        # Group by 'general_annot' and find the index of the row with the maximum value in 'baseMean'
        indices = df.groupby('miRNA_fam')['baseMean'].idxmax()

        # Complete the dictionary with the corresponding values
        miRNAs = []
        for index in indices:
            # Retrieve the row in the DataFrame at the specified index
            row = df.loc[index]
            # Extract the values in the 'general_annot' and 'log2FoldChange' columns for the current row
            general_annot = row['miRNA_fam']
            log2FoldChange = row['Shrunkenlog2FoldChange']
            # Append a tuple containing 'general_annot' and 'log2FoldChange' to the dictonary
            if file_id not in exp_miRNAs:
                # If it doesn't exist, initialize it as an empty list.
                exp_miRNAs[file_id] = []
            # Add elements
            exp_miRNAs[file_id].append((general_annot, log2FoldChange))
        
        # Get unique miRNAs
        miRNAs += set(df['miRNA_fam'].unique().tolist())

    # Obtain unique miRNA families
    miRNAs = sorted(list(set(miRNAs)))

    return exp_miRNAs, miRNAs


def create_miRNA_tables(dic_miRNAs_event: dict, miRNAs_list: list):
    '''
    This function generates two summary tables based on a dictionary that maps
    stress event identifiers to lists of differentially expressed miRNA families
    and their corresponding Shrunken log2FoldChange values.

    The two output tables are:
    1. A **presence-absence table**, indicating whether a given miRNA family is
       differentially expressed in each stress event (1 = present, 0 = absent).
    2. A **log2FoldChange table**, storing the Shrunken log2FoldChange value
       for each miRNA family in each stress event. If a miRNA family is not
       present in a specific event, 'NA' is used.

    These tables are useful for summarizing expression patterns across
    comparisons, and are saved as TSV files.

    Parameters
    ----------
    dic_miRNAs_event : dict
        Dictionary mapping stress event identifiers (str) to lists of tuples.
        Each tuple contains:
            - miRNA family name (str)
            - Shrunken log2FoldChange value (float)

    miRNAs_list : list of str
        List of all miRNA families to include as rows in the output tables. This
        ensures a uniform structure across all columns (stress events), even if
        a family is not present in some of them.

    Returns
    -------
    None
        The function writes two `.tsv` files to disk:
        - 'presence_absence_table.tsv'
        - 'shrunken_log2fc_table.tsv'
    '''

    # Sort miRNAs and ids alphabetically
    ids_list = natsorted(list(dic_miRNAs_event.keys()))
    miRNAs_list = natsorted(miRNAs_list)

    # Create final tables
    pre_abs_table = pd.DataFrame()
    log2fc_table = pd.DataFrame()

    # Iterate events ids
    for id in ids_list:
        
        # List of values associated with a specific stress event for each of the tables.
        pre_abs_values_list = []
        log2fc_values_list = []

        # Iterate miRNAs
        for miRNA in miRNAs_list:  
            
            ## 1. Obtain the values for the presence-absence table.
            # Search for the miRNA in the dictionary.
            found = any(miRNA == tupla[0] for tupla in dic_miRNAs_event[id])

            # Select a value (1 or 0) depending on whether the miRNA is represented or not.
            if found:
                pre_abs_values_list.append(1)
            else:
                pre_abs_values_list.append(0)
            
            ## 2. Obtain the values for the log2fc table.
            # Find the position in the list where the tuple associated with the desired miRNA is located
            list_index = [i for i, tuple_miR in enumerate(dic_miRNAs_event[id]) if tuple_miR[0] == miRNA]

            # If this tuple exists...
            if list_index:
                # Get log2fc value from tuple
                log2fc = dic_miRNAs_event[id][list_index[0]][1]
                log2fc_values_list.append(log2fc)
            else:
                log2fc_values_list.append('NA')    

        # Add event column to final table
        pre_abs_table[id] = pre_abs_values_list
        log2fc_table[id] = log2fc_values_list
    
    # Add miRNA families as row names (pre-abs table)
    pre_abs_table.index = miRNAs_list
    pre_abs_table.index.name = 'miRNA_fam'

    # Add miRNA families as row names (log2fc table)
    log2fc_table.index = miRNAs_list
    log2fc_table.index.name = 'miRNA_fam'

    # Save table in csv file
    pre_abs_table.to_csv('presence_absence_table.tsv', sep='\t') 
    log2fc_table.to_csv('shrunken_log2fc_table.tsv', sep="\t")


## MAIN PROGRAM

def main():
    '''
    Main program
    '''
    
    # Arguments
    parser = argparse.ArgumentParser(prog='08-Build_miRNA_vs_stress_tables.py', 
                                     description='''
                                        This function creates two types of tables: a presence-absence table and a \
                                        Shrunken Log2FC table. The presence-absence table represents whether a miRNA \
                                        is present (1) or absent (0) in a specific stress event. In contrast, the \
                                        Shrunken Log2FC table includes the Shrunken Log2FC values of those miRNAs \
                                        that are represented in the respective stress event. In both tables the \
                                        events are placed in the columns using identifiers created by the program, \
                                        where the species, the stress and the event in question are considered. \
                                        Meanwhile, the rows show the miRNA families annotated in all the experiments \
                                        considered. This table indicates that at least one member of the miRNA \
                                        family has been differentially expressed in a given stress event, \
                                        represented by a 1 (or its Shrunken log2FC), while the absence of these \
                                        miRNAs is indicated by a 0 (or NA). This program also generates a table \
                                        in which the identifiers created from the program are related to the \
                                        stress events they are linked to, among other things. ''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--annotation',
        action='append',
        help='Path to an annotation file (in TSV format). Must be provided' \
            'alongside a corresponding --metadata file.'
    )
    parser.add_argument(
        '--metadata',
        action='append',
        help='Path to the metadata file corresponding to the given' \
            'annotation file.'
    )
    parser.add_argument(
        '--samples',
        action='append',
        help='String representing the list of samples used in the' \
        ' contrast performed in the differential expression analysis' \
        '(The sample identifiers are separated by commas).'
    )
    parser.add_argument(
        '--fields',
        type=str,
        required=True,
        help=(
            'Comma-separated list of metadata fields to be considered for generating unique '
            'identifiers for each stress event (e.g., "Tissue,Genotype"). These fields '
            'will be combined with mandatory fields like Species and Treatment.'
        )
    )
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    # Check the arguments and store them in a variable
    try:
        files_list = args.annotation
        meta_list = args.metadata
        samples_list = args.samples
        fields = args.fields.split(',')
    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()

    # These two fields must be considered to create the identifiers
    mandatory_fields = ['Species', 'Treatment']

    # Combine mandatory fields with user-defined fields, ensuring no duplicates
    combined_fields = []
    seen = set()
    for field in mandatory_fields + fields:
        field_low = field.lower()
        if field_low not in seen:
            combined_fields.append(field_low)
            seen.add(field_low)
    
    # Create an identifier for each stress event (Comparison)
    print('Creanting stress event ids...')
    ids_dic = get_events_ids(files_list, meta_list, samples_list, combined_fields)
    print('Done!')

    # Get the list of miRNAs differentially expressed in each stress event
    print('Relating stress event identifiers to miRNA families...')
    files_miRNA_list, miRNA_list = get_miRNA_families(files_list)
    print('Done!')

    # Create the event_id-miRNAs dictionary (e.g. {'1.1.1.1': [('miR156', 3.23), ('miR472', -1.43)]})
    eventid_miRNAs_dic = {}
    for key, new_key in ids_dic.items():
        if key in files_miRNA_list:
            eventid_miRNAs_dic[new_key] = files_miRNA_list[key]

    # Create the both presence-absence and Log2fc matrices
    print('Creating miRNA-events matrices...')
    create_miRNA_tables(eventid_miRNAs_dic, miRNA_list)
    print('Done!')

## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
