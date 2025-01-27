#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   08-Build_presence_absence_table.py
#
#   This function creates two types of tables: a presence-absence table
#   and a Shrunken Log2FC table. The presence-absence table represents
#   whether a miRNA is present (1) or absent (0) in a specific stress
#   event. In contrast, the Shrunken Log2FC table includes the Shrunken
#   Log2FC values of those miRNAs that are represented in the respective
#   stress event. In both tables the events are placed in the columns
#   using identifiers created by the program, where the species, the
#   stress and the event in question are considered. Meanwhile, the
#   rows show the miRNA families annotated in all the experiments
#   considered. This table indicates that at least one member of the
#   miRNA family has been differentially expressed in a given stress
#   event, represented by a 1 (or its Shrunken log2FC), while the
#   absence of these miRNAs is indicated by a 0 (or NA). This program
#   also generates a table in which the identifiers created from the
#   program are related to the stress events they are linked to, among
#   other things. 
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 10/26/2023
#   Version: 2.0
#
#******************************************************************************


## IMPORT MODULES

import argparse
from natsort import natsorted
import pandas as pd
import sys
import warnings
import re
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)


## FUNCTIONS

def get_events_ids(path_in: str, path_out: str):
    '''
    This function creates an identifier for each stress event. These
    identifiers consist of 3 numbers separated by ".". The first one refers to
    the species, the second one to the type of stress and the third one to the
    stress event. Finally, this function returns a dictionary whose key is the
    name of the file associated to that stress event (e.g. PRJNA277424_1_2)
    and as value a list containing the new generated identifier, the species,
    the stress and the stress event itself.

    Parameters
    ----------
    path_in : str
        Absolute path of the differential expression analysis summary file,
        which relates the file name to the stress event. 
    path_out : str
        Absolute path of the output file of the table with the new identifiers,
        files name, species, stresses and stress events.

    Returns
    -------
    dict
        Dictionary with the new identifiers, the files name, species, stresses
        and stress events.
    '''

    # Read DEA summary table
    comp_table = pd.read_csv(path_in, sep='\t', usecols=[0,1,2,3])
    comp_table = comp_table.reset_index()  # make sure indexes pair with number of rows

    # Get species-id dictionary
    species_list = comp_table['Species'].unique().tolist()
    dic_species = { sp : i + 1 for i, sp in enumerate(species_list) }

    # Get stress-id dictionary
    comp_array = comp_table['Comparison']
    stress_list  = sorted(set([c.split("_")[2] for c in comp_array]))
    dic_stress = { sts : j + 1 for j, sts in enumerate(stress_list) }

    # Get id-event dictionary
    dic_events = {}
    dic_count = {}
    for _, row in comp_table.iterrows():

        # Get species and stress
        sp = row['Species']
        stress = row['Comparison'].split("_")[2]
        
        # Build species-stress id
        id = f'{str(dic_species[sp])}.{str(dic_stress[stress])}'
        
        # Build event id
        if id not in dic_count:
            dic_count[id] = 1
            event_id = f'{id}.{str(dic_count[id])}'
        else:
            dic_count[id] += 1
            event_id = f'{id}.{str(dic_count[id])}'
        
        # Get experiment file name
        exp_filename = row["Event"]

        # Save results in dictionary
        dic_events[exp_filename] = (event_id, sp, stress, row['Comparison'])

    with open(path_out, 'a') as file_out:
        # Write header
        file_out.write('File\tID\tSpecies\tStress\tComparison\n')
        # Write body
        for event_id in dic_events:
            file_out.write(event_id + '\t' + dic_events[event_id][0]+ '\t' +
                        dic_events[event_id][1]+ '\t' +
                        dic_events[event_id][2] + '\t' +
                        dic_events[event_id][3] + '\n')
    
    return dic_events


def get_miRNA_families(annot_files_list: list, dic_id_event: dict):
    '''
    This function receives as parameters the absolute path of the differentially
    expressed miRNAs annotation directory and a dictionary in which the name of
    the file associated with a given stress event (e.g. PRJNA277424_1_2) is
    related to its new identifier (1.4.7), the species, the stress and the
    stress event itself. From these two variables, this function creates a 
    dictionary with the miRNA families (values) that are differentially
    expressed in each stress event (key) and a list of tuples with the miRNA
    family represented in that stress event and its corresponding Shrunken
    log2FoldChange value (selected from the most representative member of the
    family). This information is useful when the desired table includes
    log2fc values instead of presence-absence data.

    Parameters
    ----------
    path_in_annot : str
        Absolute path of the differentially expressed miRNAs annotation directory.
    dic_id_event : dict
        Dictionary in which the name of the file associated with a given stress
        event (key) is related to its new identifier (1.4.7), the species, the
        stress and the stress event (Values stored in a list in the same order).

    Returns
    -------
    dict
        Dictionary whose keys are the identifiers of the stress events in
        question, and the value associated with each of these keys is a list of 
        tuples. These tuples contain, in the first element, one of the miRNAs
        represented in that stress event, and in the second element, the
        corresponding Shrunken log2FoldChange value (selected from the most
        representative member of the family).
    list
        List with all the annotated miRNA families.
    '''

    # Iterate species
    miRNAs = []
    exp_miRNAs = {}
    id_pattern = re.compile(r"^(PRJNA\d+(?:_\d+){0,2})")

    for file in annot_files_list:

        # Get the experiment name
        match = id_pattern.match(file)
        id = match.group(1)

        print(id)

        # Read the file
        df = pd.read_csv(file, sep='\t')

        print(df)

        # Group by 'general_annot' and find the index of the row with the maximum value in 'baseMean'
        indices = df.groupby('general_annot')['baseMean'].idxmax()

        print(indices)

        # Complete the dictionary with the corresponding values
        for index in indices:
            # Retrieve the row in the DataFrame at the specified index
            row = df.loc[index]
            # Extract the values in the 'general_annot' and 'log2FoldChange' columns for the current row
            general_annot = row['general_annot']
            log2FoldChange = row['Shrunkenlog2FoldChange']
            # Append a tuple containing 'general_annot' and 'log2FoldChange' to the dictonary
            if id not in exp_miRNAs:
                # If it doesn't exist, initialize it as an empty list.
                exp_miRNAs[id] = []
            # Add elements
            exp_miRNAs[id].append((general_annot, log2FoldChange))
        
        # Get unique miRNAs
        miRNAs += set(df['general_annot'].unique().tolist())

    # Obtain unique miRNA families
    miRNAs = sorted(list(set(miRNAs)))

    return exp_miRNAs, miRNAs


def create_miRNA_tables(dic_miRNAs_event: dict, miRNAs_list: list):
    '''
    This function creates two types of tables: a presence-absence table and a
    Shrunken Log2FC table. The presence-absence table represents whether a miRNA
    is present (1) or absent (0) in a specific stress event. In contrast, the
    Shrunken Log2FC table includes the Shrunken Log2FC values of those miRNAs
    that are represented in the respective stress event. In both tables the
    events are placed in the columns using identifiers created by the program,
    where the species, the stress and the event in question are considered.
    Meanwhile, the rows show the miRNA families annotated in all the experiments
    considered. This table indicates that at least one member of the miRNA
    family has been differentially expressed in a given stress event,
    represented by a 1 (or its Shrunken log2FC), while the absence of these
    miRNAs is indicated by a 0 (or NA).

    Parameters
    ----------
    dic_miRNAs_event : dict
        Dictionary with the miRNA families (values) that are differentially
        expressed in each stress event (key)
    miRNAs_list : list
        List with all the annotated miRNA families.
    path_dir_out : str
        Path of the directory where the presence-absence tables and Shrunken
        log2fc will be stored.
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
    
    parser.add_argument('-d', '--path-dea', type=str, nargs=1,  
                        help='Absolute path of the differential expression \
                        analysis summary file, which relates the file name \
                        to the stress event. ')
    parser.add_argument('-a', '--annot-files', type=str, nargs='+', 
                        help='List of files with the annotation and \
                        classification into families of differentially \
                        expressed miRNAs.')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    # Check the arguments and store them in a variable
    try:
        path_in_sum_dea = args.path_dea[0]
        annot_files_list = args.annot_files

    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()
    
    # Get ids-event dictionary
    print('Creanting stress event ids...')
    dic_id_event = get_events_ids(path_in_sum_dea, 'ids_table.tsv')
    print('Done!')
    
    # Get miRNAs and the events in which they are differentially expressed
    print('Relating stress event identifiers to miRNA families...')
    dic_miRNAs_event, miRNAs_list = get_miRNA_families(annot_files_list, dic_id_event)
    print('Done!')

    # Create Presence-Absence table
    print('Creating presence-absence table...')
    create_miRNA_tables(dic_miRNAs_event, miRNAs_list)
    print('Done!')


## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
