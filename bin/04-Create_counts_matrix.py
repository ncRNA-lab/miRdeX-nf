#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#******************************************************************************
#  
#   sRNA_Counts.py
#
#   This program generates the absolute counts and Reads Per Million (RPM)
#   tables for a specific project using the trimmed and filtered libraries,
#   also calculating the averages of both types of counts in the different
#   replicates of each condition for the sequences analysed. To do this,
#   a series of procedures are carried out:
#   
#   1. Filter the libraries by RNAcentral (Optional)
#
#   Align the sequences of each library with the rRNA, tRNA, snRNA and
#   snoRNA sequences contained in the RNAcentral database for identification
#   and elimination.
#
#   2. Create the absolute counts and RPM tables for each library.
#
#   Read filtered libraries and create absolute counts and RPM tables
#   in CSV format. These tables have two columns: seq and counts (absolute
#   counts)/RPM (Reads per Million).
#   
#               RPM = absolute count * 1000000 / size library.
#  
#   3. Create count tables for each project and subproject.
#
#   Once the absolute counts and RPM have been calculated for the sequences
#   of each library, a table is created for each type of count by joining
#   the results of the libraries belonging to the same project. Additionally,
#   count tables are also generated for each of the subprojects within the
#   project in question. These tables consist of the libraries from one of
#   the control groups and the libraries from the treatment groups associated
#   with that control group. The subproject tables are filtered to remove
#   sequences with a low count number. Sequences that do not have more
#   than 5 counts in at least 5 libraries are discarded. If the project has
#   fewer than 5 libraries, this condition is applied using the maximum
#   number of libraries in the project.
#
#
#   Authors: Antonio Gonzalez Sanchez, Pascual Villalba Bermell (pvbermell)
#   Date: 02/12/2024
#   Version: 3.0
#
#******************************************************************************


## IMPORT MODULES

import argparse
import csv
import itertools
import numpy as np
from numpy import genfromtxt
import pandas as pd
import re
import sqlite3
from sqlite3 import Error
import sys
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


## FUNCTIONS

### 1. SQLITE CONNECT, INSERT AND JOIN FUNCTIONS
def connect_to_database (database_name: str):
    '''
    This function establishes a connection to an SQlite database, creating it
    if it does not exist.

    Parameters
    ----------
    database_name : str
        Absolute database path

    Returns
    -------
    Connection
        Conexion
    Cursor
        Cursor
    '''
    try:
        ## Connect to SQlite
        sqliteConnection = sqlite3.connect(database_name)
        cursor = sqliteConnection.cursor()
        
    except Error as error:
        print('ERROR. Unable to connect to the database')
        print(f'ERROR:\n{error}\n')
        sys.exit()

    else:
        return sqliteConnection, cursor

def insert_to_database (database_name: str, table_name: str, data_path: str,
                        type_data: str='raw', sep: str='\t', chunksize=1000000) -> None:
    '''
    This function inserts data into a specified SQLite database table.

    Parameters
    ----------
    database_name : str
        SQLite database name (filename)
    table_name : str
        Name of the table in which the data will be inserted
    data_path : str
        Absolute path of the file that contains the data to insert
    type_data : str
        Type of data to be inserted. In this program you can insert a library
        of sequences in .tsv format (type_data="sequences"), a table of
        absolute counts (type_data="raw") or a table of rpm
        ("type_data="rpm"). By default "sequences.
    '''
    # Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)

    # Select type of data
    if type_data == 'raw':
        query_section = '(seq TEXT, raw INT);'
    elif type_data == 'rpm':
        query_section = '(seq TEXT, rpm REAL);'
    
    # Read data to insert
    data_to_insert = pd.read_csv(data_path, sep=sep, chunksize=chunksize, low_memory=False)

    try:
        # Pragma adjust
        cursor.execute('PRAGMA synchronous = OFF')
        cursor.execute('PRAGMA journal_mode = OFF')

        # Create table
        cursor.execute(f'CREATE TABLE IF NOT EXISTS {table_name}{query_section}')

        # Insert chunks in table
        for chunk in data_to_insert:
            # Insert chunk data into table
            chunk.to_sql(name=table_name, con=sqliteConnection, if_exists='append', index=False)       
        print(f'Data inserted in {table_name} correctly!')
        sys.stdout.flush()

        # Create table index to increase query speed
        print(f'Creating index of {table_name}')
        sys.stdout.flush()
        cursor.execute(f'CREATE INDEX idx_{table_name} ON {table_name} (seq);')
        print('DONE!\n')
        sys.stdout.flush()
        exit = False

    except Error as error:
        print('ERROR: Something went wrong during the data insertion.')
        print(f'ERROR:\n{error}\n')
        exit = True
    
    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()


def merge_counts_tables (database_name: str, data_in: list, type_data: str,
                         type_tables: str, mode: str='outer', final_table = 'project_table') -> None:
    '''
    This function generates and executes the necessary SQLite queries to join in
    the same table the different replicates of the same condition or multiple
    tables of different conditions that contain the data of all their replicates.
    The first of the two processes can be performed using the INNER JOIN (PCA)
    or FULL OUTER JOIN (Differential Expression Analysis) methods. However, the
    second one is done by default using the FULL OUTER JOIN method.

    Parameters
    ----------
    database_name : str
        Absolute path to the database where the tables to be joined are
        located.
    data_in : list
        List with the name of replicates of each condition (t_1_r_1, t_1_r_2,
        t_2_r_1, t_2_r_2)
    type_data : str
        Type of data. Type_data can be "raw" or "rpm".
    type_tables:
        Type of tables to join. type_tables can be "replicates" or "conditions".
    mode : str
        Method by which replicates of the same condition will join if
        type_tables is "replicates". If mode = "inner" the INNER JOIN method
        will be used, while if mode = "outer" the FULL OUTER JOIN method will
        be used. This parameter defaults to "inner".
    final_table : str
        Name of the final table. This argument should only be used when
        type_tables = 'conditions'. By default it is 'project_table'.
    
    Returns
    -------
    None
    '''


    # Variables
    query_list = []
    sample_dict = {}
    if mode == 'outer':
        mode = ' FULL OUTER JOIN '
    elif mode == 'inner':
        mode = ' INNER JOIN '

    # Create dictionary with conditions names as keys and replicates list of
    # the condition as values.
    for sample in data_in:
        # Create condition's name
        sample_elements = sample.split("_")
        del sample_elements[3]
        new_sample_name = '_'.join(sample_elements)
        # Count the number of replicates of each condition
        if new_sample_name in sample_dict:
            sample_dict[new_sample_name].append(sample)
        else:
            sample_dict[new_sample_name] = [sample]   

        
    ### 1. MERGE REPLICATE TABLES
    ###########################################################################
    if type_tables == 'replicates':

        for key in sample_dict:

            # Samples list of the same condition
            sample_list = sample_dict[key]
            sample_list.sort()

            # New table = condition name
            new_table = key

            # Create query
            for j in range(len(sample_list)):
                if j == 0:
                    # First iteration
                    create = f'CREATE TABLE {new_table} AS '
                    select_seq = 'SELECT COALESCE(t1.seq' 
                    select_counts = f't1.{type_data} AS {sample_list[j]}'
                    on_section = '(t1.seq'
                    from_q = f' FROM {sample_list[j]} t1'
                    join = str()
                else:
                    # Complete query sections with the info from the rest of the replicates
                    select_counts = f'{select_counts}, t{str(j + 1)}.{type_data} AS {sample_list[j]}'
                    select_seq = f'{select_seq}, t{str(j + 1)}.seq'
                    join = f'{join}{mode}{sample_list[j]} t{str(j + 1)} ON {on_section}) = t{str(j + 1)}.seq'
                    # If second iteration
                    if j == 1:
                        on_section = f'COALESCE(t1.seq, t{str(j + 1)}.seq'
                    else:
                        on_section = f'{on_section}, t{str(j + 1)}.seq'
                   
        
            # Save query in query_list
            query = f'{create}{select_seq}) AS seq, {select_counts}{from_q}{join}'
            query_list.append(query.strip() + ";")

            # Remove tables used for joining
            for old_table in sample_list:
                drop_query = f'DROP TABLE {old_table};'
                query_list.append(drop_query)

            # Create Index of the table
            query_list.append(f'CREATE INDEX idx_{new_table} ON {new_table} (seq);')

    ### 2. MERGE CONDITION TABLES
    ###########################################################################      
    else:
        # First step in the query build
        first = True
        first_on = True

        # Indicates the current condition
        count_key = 1

        # Iterate conditions-replicates dictionary (t_1_r: t_1_r_1, t_1_r_2)
        for key in sample_dict:

            # Samples list of the same condition
            sample_list = sample_dict[key]
            sample_list.sort()

            # Create query
            for i in range(len(sample_list)):

                # Write only when a new dictionary key is used and it is not the first one.
                if i == 0 and not first:
                    select_counts = f'{select_counts}, t{str(count_key + 1)}.{sample_list[i]} AS {sample_list[i]}'
                    select_seq = f'{select_seq}, t{str(count_key + 1)}.seq'
                    # First appearance of FULL OUTER JOIN in the query
                    if first_on:
                        join = f'{join} FULL OUTER JOIN {key} t{str(count_key + 1)} ON {on_section}) = t{str(count_key + 1)}.seq'
                        first_on = False
                    # Others
                    else:
                        join = f'{join} FULL OUTER JOIN {key} t{str(count_key + 1)} ON COALESCE{on_section}) = t{str(count_key + 1)}.seq'
                    on_section = f'{on_section}, t{str(count_key + 1)}.seq'
                    count_key += 1
                    
                # If it is the first step in the query build
                elif first:
                    create = f'CREATE TABLE {final_table} AS '
                    select_seq = 'SELECT COALESCE(t1.seq' 
                    select_counts = f't1.{sample_list[i]} AS {sample_list[i]}'
                    on_section = '(t1.seq'
                    from_q = f' FROM {key} t1 '
                    join = str()
                    first = False
                    
                # All other cases
                else:
                    select_counts = f'{select_counts}, t{str(count_key)}.{sample_list[i]} AS {sample_list[i]}'
        
        # Save query in query_list
        query = f'{create}{select_seq}) AS seq, {select_counts}{from_q}{join}'
        query_list.append(query.strip() + ";")

         # Create Index of the table
        query_list.append(f'CREATE INDEX idx_{final_table} ON {final_table} (seq);')


    ### 3. EXECUTE SQLITE QUERIES
    ###########################################################################

    # Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)
    
    try:
        for query in query_list:
            cursor.execute(query)
        exit = False

    except Error as error:
        print("ERROR: Something went wrong during the join of the tables.")
        print(f'ERROR:\n{error}\n')
        exit = True
    
    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()

        # If it fails, exit the program
        if exit:
            sys.exit()


def rep_counts_avg(database_name: str, table: str, new_table: str) -> str:
    """
    This function calculates the average counts (avg) of each sequence 
    belonging to the same condition (e.g. Control or Treated), that is,
    the average of the replicates (e.g. control1, control2). To do this, it
    accesses a table contained in a specified database, calculates the average
    of the counts of each sequence, and generates a table in the same database
    with the result (table_average).

    Parameters
    ----------
    database_name : str
        Absolute database path
    table : str
        Name of the input table.

    Returns
    -------
    new_table : str
        Name of the table containing the average of the counts of each sequence
        in the different conditions.
    """
        
    ## Connect to database
    sqliteConnection, cursor = connect_to_database(database_name)

    ## Get the counts table
    try:
        # Execute query
        table_sql = cursor.execute(f'SELECT * FROM {table}')

    except Error as error:
        # Error
        print(f'ERROR:\n{error}\n')

        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()

        # If it fails, exit the program
        sys.exit()
    
    ## Get condition replicates columns
    dic_colnames = {}
    seq_index = []
    for i in range(len(table_sql.description)):
        # Get column name
        col = table_sql.description[i][0]
        # Discard sequence columns
        if col[:3] == "seq":
            seq_index.append(i)
            continue
        # Remove replicate from the name
        col_elements = col.split("_")
        del col_elements[3]
        new_col = "_".join(col_elements)

        # Store in dictionary the columns of each condition (column = replicate)
        if new_col not in dic_colnames:
            dic_colnames[new_col] = [col]
        else:
            dic_colnames[new_col].append(col)

    ## Build Query
    final_columns_query = str()
    first_con = True

    # Iterate condition dictionary
    for con in dic_colnames:
        count_reps = int()
        new_colum_query = str()

        # Iterate replicates of each condition
        for j in range(len(dic_colnames[con])):
            count_reps += 1
            if j == 0:
                new_colum_query += f'CAST(({dic_colnames[con][j]}'
            else:
                new_colum_query += f' +  {dic_colnames[con][j]}'

        # new_column_query example: (control1 + control2) / 2 AS control
        new_colum_query += f') AS REAL ) / {str(count_reps)} AS {con}'

        if first_con:
            final_columns_query += new_colum_query
            first_con = False
        else:
            final_columns_query += ', ' + new_colum_query

    # Create final query
    query = f'CREATE TABLE IF NOT EXISTS {new_table} AS SELECT seq, {final_columns_query} FROM {table}'

    ## Execute final query
    try:
        # Execute query
        cursor.execute(query)
        exit = False

    except Error as error:
        print(f'ERROR:\n{error}\n')
        exit = True

    finally:
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()

    return new_table


### 2. FORMAT CONVERSION FUNCTIONS

def write_from_sql(database: str, table: str, path_out: str, columns: list,
                   shortened_sample_dic: dict, filter: bool=False,
                   filter_num_counts: int= 5, filter_num_samples: int=5) -> None:
    """
    This function writes counts tables from a specific Sqlite database in .csv
    file.

    Parameters
    ----------
    database : str
        Absolute database path
    table : str
        Database table to write in .csv file.
    path_out : str
        Absolute path of the output file.
    columns : list
        List of table columns to be written to the new .csv file.
    shortened_sample_dic:
        Dictionary containing the original sample names and their associated
        shortened names.
    """

    ## 1. Connect to database
    sqliteConnection, cursor = connect_to_database(database)

    ## 2. Create query
    # Create a string with the columns to write from the table
    columns_str = 'seq'
    for col in columns:
        columns_str += f', ifnull({col}, 0) AS {col}'

    # Complete query with columns_str
    query = f'SELECT {columns_str} FROM {table}'

    ## 3. Execute query
    try:
        cursor.execute(query)
        exit = False

    except Error as error:
        print('ERROR. Something went wrong with the creation of the new .csv file')
        print(f'ERROR:\n{error}\n')
        exit = True  

    else:
        ## 4. Write table in .tsv file
        print('Exporting data into TSV...')
        sys.stdout.flush()

        with open(path_out, 'w') as csv_file:
            csv_writer = csv.writer(csv_file, delimiter='\t', lineterminator='\n')

            # Get original columns names
            original_columns = ['seq']
            for col in columns:
                red_name = col # p.e. t_1_r_3
                sample_name = shortened_sample_dic[red_name] # p.e. control_salt_1_5d
                original_columns.append(sample_name)

            # Write the header
            csv_writer.writerow(original_columns)
            
            # Sequences will be filtered to remove those with low count numbers
            if filter:

                # Iterate SQLite cursor object          
                for i, line in enumerate(cursor):

                    # Select sequences that have more than filter_num_counts
                    # counts in at least filter_num_samples samples
                    count = len([j for j in line[1:] if j > filter_num_counts])
                    if count >= filter_num_samples:
                        csv_writer.writerow(line)

                    # If there are less than filter_num_samples samples, check that the condition is
                    # fulfilled for all of them.
                    elif count < filter_num_samples and count == len(line) - 1:
                        csv_writer.writerow(line)
            else:
                # Write columns names and data
                csv_writer.writerows(cursor)
    
    finally:    
        # Commit work and close connection
        sqliteConnection.commit()
        sqliteConnection.close()
        
        # If it fails, exit the program
        if exit:
            sys.exit()
    

### 5. METADATA FUNCTIONS
                    
def get_project_metadata (path: str) -> pd.DataFrame:
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

        # Remove spaces (strip)
        df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    except Exception as e:
        print('Unable to read metadata table.')
        print(f'Exception: {e}')
        sys.exit()
    else:
        return df


### 6. OTHER FUNCTIONS

def parse_fasta_file(fasta_file: str):
    '''
    This function parse a fasta file returning a generator object containing
    the headers and the sequences. If this object is iterated outside the
    function, we can access a new sequence and its header at each iteration
    of the loop.

    Parameters
    ----------
    fasta_file : str
        Absolute path of the fasta file
    '''
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
        sys.exit()


## MAIN PROGRAM

def main():
    '''
    Main program
    '''
    parser = argparse.ArgumentParser(prog='sRNA_counts', 
                                     description='''This program filters sRNA \
                                        sequences using RNAcentral to discard \
                                        rRNA, tRNA, snRNA and snoRNA sequences \
                                        present in the libraries of a specific \
                                        project, generates tables of absolute \
                                        counts and reads per million for each \
                                        of them, and then joins these tables \
                                        to form a single table for the project.''',  
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-p', '--project', type=str, nargs=1, 
                        help='Project id.')
    parser.add_argument('-c', '--counts-tsv', type=str, nargs='+',  
                        help='List of TSV files with absolute counts or  Reads per million (RPM) of \
                            a set of sequences. These tables should consist of \
                            two columns, seq and counts for absolute counts and \
                            seq and RPM for RPM.')
    parser.add_argument('-m', '--metadata', type=str, nargs=1,  
                        help='Absolute path of metadata samples table.')
    parser.add_argument('-g', '--valid-groups', type=str, nargs='+',  
                        help='List of valid groups of the project.')
    parser.add_argument('-r', '--rpm', action='store_true',
                        help='This option should be specified when the input tables used are RPM-based and not absolute counts')
    parser.add_argument('-j', '--project-table', action='store_true',
                        help='This option should be added if one wishes to \
                            generate a table of counts or RPM for the entire \
                            project. By default, this option is disabled.')
    parser.add_argument('-a', '--avg-project', action='store_true',  
                        help='This option should be added if one wishes to \
                            generate tables with the averages of the replicates\
                            for each conditionat at the project level. By default, this option is disabled.')
    parser.add_argument('-v', '--avg-subproject', action='store_true',  
                        help='This option should be added if one wishes to \
                            generate tables with the averages of the replicates\
                            for each condition at the subproject level. By default, this option is disabled.')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    args = parser.parse_args()
    
    
    ###########################################################################
    #                        1. CHECK ARGUMENTS                               #
    ###########################################################################

    try:
        project = args.project[0]
        counts_tsv_files = args.counts_tsv
        metadata = args.metadata[0]
        valid_groups = args.valid_groups
        rpm = args.rpm
        project_table = args.project_table
        avg_table_project = args.avg_project
        avg_table_subproject = args.avg_subproject

    except:
        print('ERROR: You have inserted a wrong parameter or you are missing a parameter.')
        parser.print_help()
        sys.exit()    

    ############################################################################
    #   1. BUILD SHORTENED NAMES AND LOAD COUNTS TABLE INTO SQLITE DATABASE    #
    ############################################################################

    # To verify the type of the input data.
    if rpm:
        data_type = 'rpm'
    else:
        data_type = 'raw'

    print(data_type)

    # Variables
    shortened_sample_dic = {}           # t_1_r_3 = SRRXXXXX1
    shortened_sample_list = []          # [t_1_r_1, t_1_r_2, ...]
    shortened_condition_dic = {}        # t_1_r = SRRXXXXX1/SRRXXXXX2/SRRXXXXX3
    samples_groups = {}                 # t_1_r = ["SRRXXXXX1", "SRRXXXXX2", "SRRXXXXX3"]
    mode = 'outer'                      # Join mode


    # 1.1 Find out which sample group each sample belongs to
    #####################################################################

    # Get the metadata dataframe
    metadata_df = get_project_metadata(metadata)

    # Group by the relevant columns (remembering that pandas starts at 0)
    # The columns we're interested in are the ones with indices: 4 (Treatment), 6 (Level), 7 (Time),
    # 8 (Cultivar), 9 (Tissue), 10 (Stage), 11 (genotype)

    grouped = metadata_df.groupby(['Group','Treatment', 'Level', 'Time',
                                   'Cultivar', 'Tissue', 'Stage', 'Genotype'])

    # Iterate over the groups and add each group to the dictionary with its respective name and 'Run' list
    for counter, (_, group_df) in enumerate(grouped, start=1):

        # Create the name of the sample group
        group_name = f"t_{counter}_r"
        
        # Get the runs that belong to the group
        run_list = group_df['Run'].tolist()
        replicate = group_df['Replicate'].tolist()

        # Get one of the samples of the group
        run_sample = run_list[0]

        # Check if the run was not previously saved
        if run_sample not in shortened_sample_dic.values():

            # Add the RUNs (list) to the dictionary
            samples_groups[group_name] = run_list

            # Add the RUNs (string) to the dictionary (for AVG tables)
            shortened_condition_dic[group_name] = "/".join(run_list)

            # Link the run name with a simplified name of the sample
            for i in range(len(run_list)):
                shortened_sample_dic[f'{group_name}_{replicate[i]}'] = run_list[i]
    
    ## 1.2 Insert the counts tables into the SQLite database
    #######################################################################

    # Sort the files list
    counts_tsv_files.sort()

    # Input samples list
    run_input_list = []

    # Iterate a number of times equal to len(list_abs) and len(list_rpm)
    for counts_file in counts_tsv_files:
        
        # Get the run name from the file name (delete the suffix)
        run_name = re.split('.raw|.rpm', counts_file)[0]

        # Add the run name to the run names list
        run_input_list.append(run_name)

        # Get the shortened sample name of the RUN
        shortened_sample_name = next((key for key, value in shortened_sample_dic.items() if value == run_name), None)

        # Save shortened_sample_name into a list
        shortened_sample_list.append(shortened_sample_name)

        # Create the name of the SQLite database
        db_name = f'db_{mode}_{data_type}.db' # e.g. db_outer.db

        # Insert absolute counts of each sample in db (1 sample = 1 table in db)
        insert_to_database(db_name, shortened_sample_name, counts_file, data_type, sep='\t')


    ############################################################################
    #                      2. JOIN THE COUNTS TABLES                           # 
    ############################################################################

    ## 2.1. Join the different replicates of each condition in a single table 
    ###########################################################################

    merge_counts_tables(db_name,                        # SQLite database
                        shortened_sample_list,          # Shortened names of all the project samples
                        data_type,                      # Two options: counts or RPM
                        'replicates',                   # Two options: replicates or conditions
                        mode)                           # Two options: outer (default) or inner

    ## 2.2. Create SUBPROJECTS groups (always)
    ###########################################################################

    # Remove undesired characters from valid groups list
    valid_groups = [int(re.sub(r"[,\[\]]", "", elem)) for elem in valid_groups]

    # This dictionary associates the list of samples with the corresponding group number
    project_groups_run_dic = metadata_df.groupby('Group')['Run'].apply(list).to_dict()

    # Remove those groups that are not valid from project_groups_run_dic 
    project_groups_run_dic = {k: v for k, v in project_groups_run_dic.items() if k in valid_groups}

    # This dictionary associates the list of samples (shortened names) with the
    # corresponding group number
    shortened_names_groups_dic = {}

    # Iterate through the groups information
    for group_num, group_runs_list in project_groups_run_dic.items():

        # Get a list with the valid Runs of the group
        group_runs_list_valid = list(set(group_runs_list) & set(run_input_list))

        # Get a list with the shortened names of the valid Runs
        group_runs_list_valid_shortened = [key for key, value in shortened_sample_dic.items() if value in group_runs_list_valid]

        # Match the group number with the list of abbreviated names of its samples
        shortened_names_groups_dic[group_num] = group_runs_list_valid_shortened


    # Check if it is desired to generate a count table for the entire project.
    if project_table:

        ## 2.3.A Create count table for the PROJECT
        ########################################################################
        
        # Final table name
        final_table = f'project_table_{data_type}'

        # Join
        print(f'Creating {project} table...')
        merge_counts_tables(db_name,
                            shortened_sample_list,
                            data_type,
                            'conditions',
                            mode='outer',
                            final_table=final_table)
        print('Done!\n')
        
        ## 2.4.A Write PROJECT table into a TSV file
        ########################################################################

        # Output paths
        table_file = f'{project}.{data_type}.tsv'

        # Sort sample list
        shortened_sample_list.sort()

        # Write SQL table into a TSV file.
        print(f'Writing {project} table to a TSV file...')
        write_from_sql(db_name,                                 # SQLite database
                       final_table,                             # Final table (Project)
                       table_file,                              # Final table output file
                       shortened_sample_list,                   # Shortened names of all the project samples
                       shortened_sample_dic)                    # Dictionary t_1_r_3 = SRRXXXXXX (key = value)
        print('Done!\n')
        
        # Check if the avg_table option has been enabled
        if avg_table_project:

            ## 2.5.B Calculate replicates average
            ################################################################

            # Out paths
            project_table_mean_file = f'{project}.{data_type}_mean.tsv'

            # Calculate average
            print(f'Creating AVG table from {project} table ...')
            rep_counts_avg(db_name,
                           final_table,
                           project + '_avg')
            print('Done!\n')
            
            # List shortened condition names
            list_condition = list(shortened_condition_dic.keys())

            # Write SQL table into a TSV file.
            print(f'Writing {project} AVG table to a TSV file...')
            write_from_sql(db_name,
                            project + '_avg',
                            project_table_mean_file,
                            list_condition,
                            shortened_condition_dic)
            print('Done!\n')
                
        ## 2.5.A Write SUBPROJECT tables into different TSV files
        ########################################################################

        # Iterate through the groups with valid samples
        for group_num, shortened_group in shortened_names_groups_dic.items():

            # Create the subproject name
            subproject_name = f'{project}_{group_num}'

            # Write subproject table (Uses the PROJECT table to generate the SUBPROJECT tables).
            print(f'Creating {subproject_name} table...')
            write_from_sql(db_name,                                                     # SQLite database
                        final_table,                                                    # Final table (Project)
                        f'{subproject_name}.{data_type}.tsv',                           # Final table output file
                        shortened_group,                                                # Group of shortened samples names
                        shortened_sample_dic)                                           # Dictionary t_1_r_3 = SRRXXXXXX (key = value)
            print('Done!\n')

            # Check if the avg_table option has been enabled
            if avg_table_subproject:

                ## 2.5.B Calculate replicates average
                ################################################################

                # Out paths
                subproject_table_mean_file = f'{subproject_name}.{data_type}_mean.tsv'

                # It only executes if the table for the entire project has not
                # been previously created and if it is the first iteration of
                # the loop (since the table for the project is used to generate
                # the tables for the subprojects)
                if not avg_table_project and i == 0:

                    # Calculate average
                    print(f'Creating AVG table from {project} table (Necessary to create the AVG tables for the subprojects) ...')
                    rep_counts_avg(db_name,
                                   final_table,
                                   project + '_avg')
                    print('Done!\n')
                
                # Create shortened condition names using the shortened sample names
                shortened_group_conditions_dup = [element[:-2] for element in shortened_group]

                # Get unique names
                shortened_group_conditions_list = list(set(shortened_group_conditions_dup))

                # Write SQL table into a TSV file.
                print(f'Writing {subproject_name} AVG table to a TSV file...')
                write_from_sql(db_name,
                               project + '_avg',
                               subproject_table_mean_file,
                               shortened_group_conditions_list,
                               shortened_condition_dic)
                print('Done!\n')

    # Alternatively, generate count tables at the subproject level only.
    else:

        # Iterate through the groups with valid samples
        for group_num, shortened_group in shortened_names_groups_dic.items():

            ## 2.3.B Create count table for the SUBPROJECT
            ####################################################################

            # Create the subproject name
            subproject_name = f'{project}_{group_num}'

            # Join all the subproject count in a single table
            print(f'Creating {subproject_name} table...')
            merge_counts_tables(db_name,
                                shortened_group,
                                data_type,
                                'conditions',
                                mode='outer',
                                final_table=subproject_name)
            print('Done!\n')
            
            ## 2.4.B Write table into a TSV file
            ####################################################################
            
            # Sort sample list
            shortened_group.sort()

            # Write subproject table (Uses the SUBPROJECT table)
            print(f'Writing {subproject_name} table to a TSV file...')
            write_from_sql(db_name,                                                     # SQLite database
                        subproject_name,                                                # Final table (Subproject)
                        f'{subproject_name}.{data_type}.tsv',                           # Final table output file
                        shortened_group,                                                # Group of shortened samples names
                        shortened_sample_dic)                                           # Dictionary t_1_r_3 = SRRXXXXXX (key = value)
            print('Done!\n')

            # Check if the avg_table option has been enabled
            if avg_table_subproject:

                ## 3.5.B Calculate replicates average
                ################################################################

                # Out paths
                subproject_table_mean_file = f'{subproject_name}.{data_type}_mean.tsv'

                # Calculate average
                print(f'Creating AVG table from {subproject_name} table ...')
                rep_counts_avg(db_name,
                            subproject_name,
                            subproject_name + '_avg')
                print('Done!\n')
                
                # Create shortened condition names using the shortened sample names
                shortened_group_conditions_dup = [element[:-2] for element in shortened_group]

                # Get unique names
                shortened_group_conditions_list = list(set(shortened_group_conditions_dup))

                # Write SQL table into a TSV file.
                print(f'Writing {subproject_name} AVG table to a TSV file...')
                write_from_sql(db_name,
                               subproject_name + '_avg',
                               subproject_table_mean_file,
                               shortened_group_conditions_list,
                               shortened_condition_dic)
                print('Done!\n')

## CALL THE MAIN PROGRAM

if __name__ == '__main__':
    '''
    Call to the main program
    '''
    main()
