#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   Filter_by_depth_rep.py
#
#   This program filters the libraries of a given project by selecting those
#   that meet a minimum threshold for sequencing depth and number of replicates
#   specified by the user. The files selected as valid will have the suffix
#   ".valid.fastq.gz". This function generates a file that indicates which
#   libraries are valid and which are not, along with the reasons; and
#   another file that specifies which groups of the project are valid and
#   which are not.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 5/11/2024
#   Version: 3.0 
#
#
#******************************************************************************


### IMPORT MODULES

import argparse
import numpy as np
import multiprocessing
import os
import pandas as pd
import subprocess
import sys
import re


### FUNCTIONS

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


def get_libraries_depth_processing(libraries_list, metadata,
                                 depth_threshold, num_process):
    '''
    This function is used to parallelize the GetLibrariesDept function, which
    obtains the sequencing depth of each library within a group of libraries.

    Parameters
    ----------
    libraries_list : list
        List with the names of compressed FASTQ files (.gz) from a group
        of libraries.
    path_project_metadata : str
        Project's metadata table to which the libraries to be
        used belong.
    depth_threshold: int
        Sequencing depth threshold.
    num_process : int
        Number of processes.

    '''

    # Create shared dictionary
    manager = multiprocessing.Manager()
    depth_info_lib = manager.dict() 


    # Calculate process distribution
    number_files = len(libraries_list)
    # If the number of libraries is <= to the max number of processes...
    if number_files <= num_process:
        # Number of processes = number of libraries
        num_process = number_files
        # And, therefore, 1 library for each process
        distribution = 1
    # If there are more than 30 libraries, num_process will always be the
    # maximum (30 in this case).
    else:
        # The libraries in the list are distributed among 30 processes
        distribution = number_files // num_process

    # "pi" would be the starting point of the list section assigned to each
    # process and "pf" the end point.
    processes = [] 
    pi = 0
    for i in range(num_process):
        
        if i == num_process - 1:
            pf = number_files
        else:
            pf = pi + distribution
        
        processes.append(multiprocessing.Process(target = get_libraries_depth,
                                                 args=(libraries_list[pi:pf],
                                                       metadata,
                                                       depth_threshold,
                                                       depth_info_lib) )) 
        processes[i].start()
        print(f'Process {i} launched.')
        
        pi = pf
        
    for p in processes:
        p.join()
    
    return depth_info_lib


def get_libraries_depth(libraries_list, metadata, depth_threshold, depth_info_lib):
    """
    This function obtains the sequencing depth of each library within a group of
    libraries (understood as the total number of sequences in the library). Once
    obtained, it enters the value into a dictionary shared between processes
    along with a string that indicates whether the library is valid (valid) or
    not (not-valid) based on a depth threshold.

    Parameters
    ----------
    libraries_list : list
        List with the names of compressed FASTQ files (.gz) from a group
        of libraries.
    metadata : pd.DataFrame
        Project's metadata table to which the libraries to be
        used belong.
    depth_threshold: int
        Sequencing depth threshold.
    depth_info_lib : multiprocessing.Manager()
        Shared dictionary among processes.
    """

    # Iterate libraries
    for path_library in libraries_list:

        # Get library name
        fastq_gz_name = os.path.basename(path_library)

        # Check if the library is in the metadata
        run_name = fastq_gz_name.split('.')[0]
        lib_in_metadata = metadata.loc[metadata['Run'].str.contains(run_name, case=False)].any().any()
        
        # The library is in the metadata
        if lib_in_metadata:
    
            # Calculate sequencing depth
            cmd = f"zcat {path_library} | wc -l"
            depth = int(subprocess.check_output(cmd, shell=True).decode('utf-8').strip()) // 4

            # Filter by sequencing depth
            if depth >= depth_threshold:
                depth_info_lib[fastq_gz_name] = [depth, 'valid']
            else:
                depth_info_lib[fastq_gz_name] = [depth, 'not-valid']
        else:
            # Write not-in-metadata) if the sample is not found in the metadata
            depth_info_lib[fastq_gz_name] = ['NA', 'not-in-metadata']
    
    return depth_info_lib

    
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

        # If there is only one factor, there must be at least two valid sample
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

    parser = argparse.ArgumentParser(prog='03-Filter_by_depth_rep.py', 
                                     description='''This program filters the
                                     libraries of a given project by selecting
                                     those that meet a minimum threshold for
                                     sequencing depth and number of replicates
                                     specified by the user. The files selected 
                                     s valid will have the suffix
                                     ".valid.fastq.gz". This function generates
                                     a file that indicates which libraries
                                     are valid and which are not, along with
                                     the reasons; and another file that
                                     specifies which groups of the project
                                     are valid and which are not.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('-i', '--project-files', type=str, nargs='+',  
                        help='List of project libraries \
                            (SRR10747100_tr.fastq.gz SRR10747101.fastq.gz...).')
    parser.add_argument('-j', '--project', type=str, nargs=1,
                        help='Project to which the libraries belong.')
    parser.add_argument('-s', '--species', type=str, nargs=1,
                        help='Species to which the libraries belong.')
    parser.add_argument('-m', '--metadata', type=str, nargs=1,
                        help='Path to the project metadata file.')
    parser.add_argument('-d', '--depth-threshold', type=int, nargs=1, 
                        help='Sequencing depth threshold. This argument must \
                            be a positive numerical value.')
    parser.add_argument('-r', '--rep-threshold', type=int, nargs=1,
                        help='Number of replicates threshold. This argument \
                            must be a positive numerical value.')
    parser.add_argument('-p', '--processes', type=int, nargs=1, default=1,
                        help='Number of processes (default: %(default)s)')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    
    ## 1. CHECK ARGUMENTS
    #######################################################################

    try:
        project_files = args.project_files
        project = args.project[0]
        species = args.species[0]
        path_metadata = args.metadata[0]
        depth_threshold = args.depth_threshold[0]
        rep_threshold = args.rep_threshold[0]
        num_process = args.processes[0]

    except:

        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()

    # Create output summary paths
    results_s_path = f'{species}_{project}.sum_projects.tsv'
    summary_s_path = f'{species}_{project}.sum_libraries.tsv'


    ## 2. CHECK AND FILTER BY SEQUENCING DEPTH
    ###################################################################

    # Get project metadata table
    project_metadata = get_project_metadata(path_metadata)
    
    # Check the sequencing depth of the libraries
    depth_info_lib = get_libraries_depth_processing(project_files, project_metadata, depth_threshold, num_process)

    ## 3. SELECT VALID LIBRARIES ACCORDING TO THEIR SEQUENCING DEPTH
    ###################################################################

    # Lists with all libraries and with valid libraries according to their depth.
    total_libraries = []
    lib_valid_names = []

    # Map the file name to the run name in a dictionary
    run_file_dic = {}

    # Iterate through the names of the runs
    for file in depth_info_lib.keys():

        # Remove the suffix/extension
        run = file.split('.')[0]

        # Update the total libraries list
        total_libraries.append(run)

        # Update the valid libraries list
        if depth_info_lib[file][1] == "valid":
            lib_valid_names.append(run)

        # Update the dictionary
        run_file_dic[run] = file

    
    ## 4. CHECK AND FILTER BY NUMBER OF REPLICATES
    ###################################################################

    # Get a list of valid libraries and groups
    filtered_libraries, filtered_groups = get_valid_libraries_by_rep(lib_valid_names, project_metadata,rep_threshold)
    print(filtered_groups)

    # Get a list of discarded libraries
    discarded_libraries = list(set(total_libraries) - set(filtered_libraries))

    ## 4. SAVE THE SELECTED FILES AND CREATE THE SUMMARY FILES
    ###################################################################

    # If valid libraries exist...
    if len(filtered_libraries) > 0:

        # Save filtered libraries
        for lib_name in filtered_libraries:

            # Create original fastq.gz name
            original_fastq_lib = run_file_dic[lib_name]

            # Add the .filt suffix
            os.system(f'mv {original_fastq_lib} {lib_name}.valid.fastq.gz')

            # Write the path of the selected files.
            with open(summary_s_path, 'a') as filtered:

                # Get depth value from dictionary of selected files
                filtered.write(f'{original_fastq_lib}\t{depth_info_lib[original_fastq_lib][0]}\t{depth_info_lib[original_fastq_lib][1]}\tvalid\n')

    # If discarded libraries exist...
    if len(discarded_libraries) > 0:

        # Save discarded libraries
        for lib_name in discarded_libraries:

            # Create original fastq.gz name
            original_fastq_lib = run_file_dic[lib_name]

            # Remove discarded libraries
            #os.system(f'rm {original_fastq_lib}')
            os.system(f'mv {original_fastq_lib} {lib_name}.notvalid.fastq.gz')

            # Write the path of the discarded files.
            with open(summary_s_path, 'a') as excluded:

                # Get depth value from dictionary of selected files
                excluded.write(f'{original_fastq_lib}\t{depth_info_lib[original_fastq_lib][0]}\t{depth_info_lib[original_fastq_lib][1]}\tnot-valid\n')

    # Save global results
    total_groups = project_metadata['Group'].unique().tolist()
    for group in total_groups:

        # All the runs of the group
        all_runs_group = project_metadata.loc[project_metadata['Group'] == group, 'Run'].tolist()

        # Valid samples from the subproject
        subproject_valid_samples = set(all_runs_group) & set(filtered_libraries)

        # Number of valid and not-valid samples
        num_valid_samples = len(subproject_valid_samples)
        num_notvalid_samples = len(all_runs_group) - len(subproject_valid_samples)

        # Determine if the group is valid or not
        validity = 'valid' if group in filtered_groups else 'not-valid'

        # Write the results to the file
        with open(results_s_path, 'a') as results:
            results.write(f'{species}\t{project}\t{group}\t{num_valid_samples}\t{num_notvalid_samples}\t{validity}\n')
     
    print(f'{project} ({species}) done!')


## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
