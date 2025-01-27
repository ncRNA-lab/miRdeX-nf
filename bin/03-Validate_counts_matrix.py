#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   03-Validate_counts_matrix.py
#
#   This program validates a given counts matrix by checking its structural 
#   integrity and ensuring it meets the necessary criteria for differential 
#   expression analysis. The script verifies that:
#   1. The first column is labeled 'seq' and contains nucleotide sequences.
#   2. All remaining columns contain integer values corresponding to absolute 
#      counts for each sample.
#   3. The sample columns match the samples in the provided metadata.
#   4. The samples in the counts matrix form at least two groups with the 
#      required minimum number of replicates for differential expression analysis.
#
#   The script generates a TSV file indicating the validity of each group within 
#   the project based on the samples present in the counts matrix.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 5/11/2024
#   Version: 1.0
#
#******************************************************************************

### IMPORT MODULES

import argparse
import numpy as np
import pandas as pd
import subprocess
import sys
import re


### FUNCTIONS

# Función para comprobar si un valor en texto puede convertirse en un entero
def string_is_valid_integer(value):
    '''
    This function checks whether the elements of the list are valid integers 
    represented as strings.

    Parameters
    ----------
    value : list
        A string that is expected to represent an integer.

    Returns
    -------
    bool
        Returns True if the value can be successfully converted to integers,
        otherwise False.
    '''

    try:
        # Intentamos convertir la cadena a un entero
        int(value)
        return True
    except ValueError:
        # Si no se puede convertir, no es un entero válido
        return False


def check_counts_matrix_structure(path_counts_matrix):
    '''
    This function checks the structure of a counts matrix in TSV format to ensure
    it meets the required conditions: the first column should contain nucleotide
    sequences, and all other columns should contain integer values. 

    It also extracts the sample names from the header of the matrix file.

    Parameters
    ----------
    path_counts_matrix : str
        Absolute path of the TSV file containing the counts matrix.

    Returns
    -------
    validity : bool
        Returns True if the counts matrix structure is valid (first column is 'seq',
        first row contains nucleotide sequences, and all other columns contain integers), 
        otherwise False.
    
    sample_names : list
        A list containing the sample names from the header of the counts matrix. If the 
        matrix is invalid, the list will be empty.
    '''
    # Default variables
    validity = True
    valid_nucleotides = {"A", "T", "C", "G"}

    # Read the two first rows of the counts matrix file
    with open(path_counts_matrix, mode='r') as file:
        
        # Get the header of the file to check the columns validity
        header = file.readline().strip().split('\t')

        # Get the first row of the file (excluding the header) to check if it is an absolute count matrix.
        first_row = file.readline().strip().split('\t')

    # If the first column is not 'seq', validity = false
    if header[0] != 'seq':
        validity = False
    
    # Check if the first column contains nucleotide sequences
    valid_sequence = all(base in valid_nucleotides for base in first_row[0])

    # Check if the values of the other columns are integers.
    valid_counts_values = all(string_is_valid_integer(x) for x in first_row[1:])

    # If either of these two variables is false, set validity to false.
    if not valid_sequence or not valid_counts_values:
        validity = False
        
    # Get the sample names list (e.g. ['SRR14182750', 'SRR14182749', ...])
    sample_names = header[1:] if validity else []

    return validity, sample_names

def get_project_metadata(path):
    '''
    This function extracts the sample metadata from the metadata table of the
    project to which it belongs. 

    Parameters
    ----------
    path : str
        Absolute path of the txt file containing the project metadata.

    Returns
    -------
    df : pd.DataFrame
        Metadata dataframe of the sample of interest.
    '''

    try:
        # Read metadata csv file
        metadata_table = np.genfromtxt(path, delimiter='\t', dtype=None, encoding=None,  names=True)

        # Create dataframe
        df = pd.DataFrame(metadata_table)

    except Exception as e:
        print('Unable to read metadata table.')
        print('Exception: ', e)
        sys.exit()
    else:
        return df

    
def get_valid_libraries_by_rep (libraries, metadata, threshold):
    """
    This function receives a list of Runs and checks if the sample groups to
    which they belong have the required minimum number of biological
    replicates using the project metadata they belong to. The function
    returns two lists: one with the list of valid samples and another
    with the valid groups of the project.

    Parameters
    ----------
    libraries : list
        A list of Runs to check for valid samples.
    metadata : pd.DataFrame
        A DataFrame containing the project metadata, including sample group
            information.
    threshold : int
        The minimum number of biological replicates required for a sample
            group to be considered valid.

    Returns
    -------
    tuple
        A tuple containing two lists:
            - valid_samples: A list of valid samples that meet the biological replicate threshold.
            - valid_groups: A list of valid groups that meet the biological replicate threshold.
    """

    # Get the group of samples
    grouped = metadata.groupby('Group')

    filtered_libraries = []
    filtered_groups = []
    for (_, group_df) in grouped:

        # Get the number of the group
        group = group_df['Group'].unique()[0]

        # Check if a specific contrast is specified 
        contrast_is_specified = (group_df['Contrast'] != 'CT.0').any()

        # Column info dictionary
        col_info_dic = {}

        # Get the designRef information of the group
        design_ref = group_df['DesignRef'].unique()[0]

        # Get the design information of each relevant column
        design_col_info = design_ref.split(":")

        # Iterate through the relevant columns information
        for col_info in design_col_info:
            
            # Separate the column name from the reference value
            splitted_col_info = re.split(r'[()]', col_info)

            # Get the colname and the reference value
            colname = splitted_col_info[0]
            ref_value = splitted_col_info[1]

            # Save the information in the dictionary
            col_info_dic[colname] = ref_value

        # Get the relevant columns
        factors = list(col_info_dic.keys())

        # Create a list with the relevant columns, and the run
        factors_and_run = factors.copy()
        factors_and_run.append('Run')

        # Create the new dataframe
        relevant_df = group_df[factors_and_run]


        ## 1. More than one factor or with a custom contrast.
        ####################################################################
        
        # If there is more than one factor, all sample groups must be valid !!
        if len(factors) > 1 or contrast_is_specified:
            
            # Create a column by combining the levels of the different factors.
            relevant_df_multiple_factor = relevant_df.copy()
            relevant_df_multiple_factor.loc[:, 'Sample_groups'] = relevant_df_multiple_factor[factors].astype(str).agg('_'.join, axis=1)

            # Filter the DataFrame to keep only the samples in 'libraries'
            filtered_df = relevant_df_multiple_factor[relevant_df_multiple_factor['Run'].isin(libraries)]

            # Group by the factor column and count the number of samples in each time group
            count_by_factor = filtered_df.groupby('Sample_groups').size().reindex(relevant_df_multiple_factor['Sample_groups'].unique(), fill_value=0).reset_index(name='Sample_Count')

            # Check if there is any level that does not meet the minimum number
            # of replicates.
            has_less_than_two = (count_by_factor['Sample_Count'] < 2).any()
            if not has_less_than_two:

                # Get the factor levels with at least rep_threshold replicates (dataframe)
                valid_levels_df = count_by_factor[count_by_factor['Sample_Count'] >= threshold]

                # Get the valid levels
                valid_levels_values = valid_levels_df['Sample_groups'].unique()

                # Get the samples from the valid levels of the factor
                samples_valid_by_rep = list(filtered_df[filtered_df['Sample_groups'].isin(valid_levels_values)]['Run'])
                
                # Add the samples to the filtered_samples list
                filtered_libraries.extend(samples_valid_by_rep)
                filtered_groups.append(group)
            
        
        ## 2. Only 1 factor 
        ####################################################################

        # If there is only one factor, there must be at least repthreshold valid sample
        # groups (one of them being the baseline level)!!
        else:

            # Filter the DataFrame to keep only the samples in 'libraries'
            filtered_df = relevant_df[relevant_df['Run'].isin(libraries)]

            # Group by the factor column and count the number of samples in each time group
            count_by_factor = filtered_df.groupby(factors[0]).size().reindex(relevant_df[factors[0]].unique(), fill_value=0).reset_index(name='Sample_Count')

            # Get the baseline level of the factor
            factor_reference = col_info_dic[factors[0]]
                
            # Get the factor levels with at least rep_threshold replicates (dataframe)
            valid_levels_df = count_by_factor[count_by_factor['Sample_Count'] >= threshold]

            # If there are valid levels
            if not valid_levels_df.empty:

                # Get the valid levels
                valid_levels_values = valid_levels_df[factors[0]].unique()

                # The baseline level must be in the list of valid levels and must not be the only one.
                if len(valid_levels_values) > 1 and factor_reference in valid_levels_values:

                    # Get the samples from the valid levels of the factor
                    samples_valid_by_rep = list(filtered_df[filtered_df[factors[0]].isin(valid_levels_values)]['Run'])

                    # Add the samples to the filtered_samples list
                    filtered_libraries.extend(samples_valid_by_rep)
                    filtered_groups.append(group)

    return sorted(list(set(filtered_libraries))), filtered_groups


## MAIN PROGRAM

def main():
    '''
    Main program
    '''

    parser = argparse.ArgumentParser(prog='03-Validate_counts_matrix.py', 
                                     description='''
                                     This program validates a given counts matrix by checking its structural
                                     integrity and ensuring it meets the necessary criteria for differential
                                     expression analysis. The script verifies that:
                                        1. The first column is labeled 'seq' and contains nucleotide sequences.
                                        2. All remaining columns contain integer values corresponding to absolute
                                           counts for each sample.
                                        3. The sample columns match the samples in the provided metadata.
                                        4. The samples in the counts matrix form at least two groups with the
                                           required minimum number of replicates for differential expression analysis.
                                    The script generates a TSV file indicating the validity of each group within
                                    the project based on the samples present in the counts matrix.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-i', '--id', type=str, nargs=1,
                        help='Project')
    parser.add_argument('-g', '--group', type=int, nargs=1,
                        help='Number of the group')
    parser.add_argument('-c', '--counts_matrix', type=str, nargs=1,  
                        help='Counts matrix file path')
    parser.add_argument('-m', '--metadata', type=str, nargs=1,
                        help='Path to the project metadata file.')
    parser.add_argument('-r', '--rep-threshold', type=int, nargs=1,
                        help='Number of replicates threshold. This argument \
                            must be a positive numerical value.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    
    ## 1. CHECK ARGUMENTS
    #######################################################################

    try:
        id = args.id[0]
        group = args.group[0]
        counts_file = args.counts_matrix[0]
        path_metadata = args.metadata[0]
        rep_threshold = args.rep_threshold[0]

    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()

    # Create group id
    group_id = f'{id}_{str(group)}'

    # Create output summary paths
    results_s_path = f'{group_id}.sum.tsv'

    ## 2. CHECK COUNTS MATRIX STRUCTURE
    ###################################################################

    validity, sample_names = check_counts_matrix_structure(counts_file)

    # If the matrix is structurally valid...
    if validity:
    
        ## 4. CHECK IF THE MATRIX MEETS THE MINIMUM NUMBER OF REPLICATES
        ###################################################################

        # Get project metadata table
        project_metadata = get_project_metadata(path_metadata)

        # Get the metadata of the group
        group_metadata = project_metadata[project_metadata['Group'] == group]

        # Check if the matrix is valid (replicates)
        _, valid_groups = get_valid_libraries_by_rep(sample_names, group_metadata, rep_threshold)

        # Get all the runs of the group
        all_runs_group = project_metadata.loc[project_metadata['Group'] == group, 'Run'].tolist()

        # If the group is valid...
        if (len(valid_groups) > 0):
            
            # Number of valid and not-valid samples
            num_valid_samples = len(sample_names)
            num_notvalid_samples = len(all_runs_group) - len(sample_names)

            # Change the name of the input file
            subprocess.run(f'mv {counts_file} {group_id}.valid.tsv', shell=True)

            # Write results in summary file
            with open(results_s_path, 'a') as results:
                results.write(f'{group_id}\t{id}\t{str(group)}\tvalid\tNA\t{num_valid_samples}\t{num_notvalid_samples}\n')

        # The group is not valid...
        else:

            # Change the name of the input file
            subprocess.run(f'mv {counts_file} {group_id}.notvalid.tsv', shell=True)

            # Write results in summary file
            with open(results_s_path, 'a') as results:
                results.write(f'{group_id}\t{id}\t{str(group)}\tnot-valid\treplicates\t0\t{len(all_runs_group)}\n')
    
    # If the matrix is not structurally valid
    else:
        # Change the name of the input file
        subprocess.run(f'mv {counts_file} {group_id}.valid.tsv', shell=True)

        # Write results in summary file
        with open(results_s_path, 'a') as results:
            print(f'{group_id}\t{id}\t{str(group)}\tnot-valid\tstructure\tNA\tNA\n')
   
## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
