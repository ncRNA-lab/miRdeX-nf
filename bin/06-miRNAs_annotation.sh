#!/bin/bash

#******************************************************************************
#  
#   miRNAs_annotation.sh
#
#   This program annotates significant differentially expressed miRNAs
#   using the miRBase, PmiREN, and sRNAanno databases without considering a
#   threshold value of log2FC. To do this, it generates a fasta file with the
#   sequences and aligns them with the miRBase, PmiREN, and sRNAanno databases
#   to identify the mature miRNAs and also their precursors. Once annotated,
#   the results of all three alignments are merged into a single five-column
#   table: sequence, Annotation with miRBase using only sequences of specific
#   species, Annotation with miRBase using the others of the species, Annotation
#   with PmiREN using only sequences of specific species, Annotation with
#   PmiREN, and Annotation with sRNAanno using only sequences of specific
#   species, Annotation with sRNAanno using the others of the species. If any of
#   the sequences has not obtained results in the alignment with one of the
#   three databases, it will be represented as NULL. Additionally, the program
#   generates files where sequences are filtered based on whether they have
#   been annotated in at least two of the three databases. This program also
#   generates a summary file with the number of sequences annotated for each of
#   the three cases: miRNAs, and precursors; and summary file with the number
#   of sequences of each length within a range of 20-25nt.
#
#   Author: Antonio Gonzalez Sanchez
#   Date: 22/12/2023
#   Version: 3.0
#
#
#******************************************************************************



### FUNCTIONS

#################################################
#
#   This function contains the usage message
#
#################################################

usage() {
    
    help='''
    usage: miRNAs_annotation.sh [options] ...
    options:
        -h                      help
        -i/--input              Path to the directory where the species
                                directories are located. Tables with the
                                results of the differential expression
                                analysis can be found in each of these
                                directories.
        -s/--species            Species identifier.
        -p/--project            Project identifier.
        -o/--output             Path to the directory where the miRNA annotation
                                results will be stored.
        -v/--mismatches         Number of mismatches allowed per Bowtie.
        -m/--mirbase            Path to the directory where the downloaded
                                miRNAs and precursors (hairpin) files from
                                miRBase are located.
        -e/--pmiren             Path to the directory where the downloaded
                                miRNAs and precursors (hairpin) files from
                                Plant miRNA ENcyclopedia (PmiREN) are located.
        -a/--srnaanno           Path to the directory where the downloaded
                                miRNAs and precursors (hairpin) files from
                                SRNAanno are located.
        -d/--species-ids        Path to the file (table) that establishes the
                                correspondence between the three-letter code
                                and four-letter code used in this pipeline to
                                name the species.
        -r/--threads            Number of threads.
    '''
    echo "$help"
    exit 0
}


#####################################################
#
#   This function manages the input of arguments
#   to the program.
#
#####################################################

arguments_management() {

    ## 1. GET ARGUMENTS
    ############################################################################

    # Set a default value for optional arguments
    projects_valid_list=''
    ea_table=''
    mww_pvalue=-1
    threads=1
    min_num_db=2
    path_mirbase=''
    path_PmiREN=''
    path_sRNAanno=''

    # Inicialize the files variable.
    files=''

    # Read the options
    TEMP=$(getopt -o h::i:d:s:m:e:a:v:n:t:l:a:w: --long help::,input:,id:,species:,mirbase:,pmiren:,srnaanno:,mismatches:,min-num-db:,threads:,projects-list:,ea-table:,mww-pvalue: -- "$@")

    # Check if the arguments are valid
    VALID_ARGUMENTS=$?
    if [ "$VALID_ARGUMENTS" != 0 ];
    then
        usage
    fi

    # Convert options to script arguments
    eval set -- "$TEMP"

    # Extract options and their arguments into variables.
    while true
    do
        case "$1" in
            -h|--help)
                usage ;;
            -i|--input)
                files+=" $2"; shift 2 ;;
            -d|--id)
                id="$2"; shift 2 ;;
            -s|--species)
                species="$2"; shift 2 ;;
            -m|--mirbase)
                # Check if srnaanno miRBase is provided
                if [ -n "$2" ] ; then
                    path_mirbase="$2"
                fi
                shift 2 ;;
            -e|--pmiren)
                # Check if srnaanno PmiREN is provided
                if [ -n "$2" ] ; then
                    path_PmiREN="$2"
                fi
                shift 2 ;;
            -a|--srnaanno)
                # Check if srnaanno database is provided
                if [ -n "$2" ] ; then
                    path_sRNAanno="$2"
                fi
                shift 2 ;;
            -v|--mismatches)
                mismatches="$2"; shift 2 ;;
            -n|--min-num-db)
                # Check if a value for -n is provided
                if [ -n "$2" ] ; then
                    min_num_db="$2"
                fi
                shift 2 ;;
            -t|--threads)
                threads="$2";
                # Check if -d is a positive number
                case "$threads" in
                    # Argument is not a positive number
                    ''|*[!0-9]*)
                        echo "Error: -t (--threads) must be a positive number"
                        usage ;;
                esac
                shift 2;;
            -l|--projects-list)
                # Check if a value for -l is provided
                if [ -n "$2" ] ; then
                    projects_valid_list="$2"
                fi
                shift 2 ;;
            -a|--ea-table)
                # Check if a value for -l is provided
                if [ -n "$2" ] ; then
                    ea_table="$2"
                fi
                shift 2 ;;
            -w|--mww-pvalue)
                # Check if a value for -p is provided
                if [ -n "$2" ]; then
                    mww_pvalue="$2"
                    # Check if -p is a positive number (including decimals)
                    if ! [[ "$mww_pvalue" =~ ^[0-9]+(\.[0-9]+)?$ ]]; then
                        echo "Error: -p (--mww-pvalue) must be a positive number"
                        usage
                    fi
                fi
                shift 2 ;;
            # -- meands the end of the arguments; drop this, and break out the while loop
            --) shift ; break ;;
            # If invalid options were passed...
            *) echo "Unexpected option: $1 - this should not happen."
                usage ;;
        esac
    done


    ## 2. CHECK THE STRATEGY THAT WILL BE USED
    ############################################################################
    # Way 1 -> Use all the subprojects (Default)
    # Way 2 -> Use the subprojects from the list
    # Way 3 -> Use the subprojects that obtained a p-value (Exploratory Analysis)
    #          lower than 'mww_pvalue'

    # Default way
    way=1

    # Is mww_pvalue less than 0? (yes = 1, no = 0)
    result_pvalue=$(echo "$mww_pvalue < 0" | bc -l)

    # Check that the variable project_list is empty
    [[ -n $projects_valid_list ]] && way=2

    # Check that variable ea_table or variable mww_pvalue is not empty.
    [[ -n "$ea_table" || "$result_pvalue" -eq 0 ]] && way=3

    # Check that variables from way 2 and 3 have not been provided.
    [[ -n "$projects_valid_list" && ( "$result_pvalue" -eq 0 || -n "$ea_table" ) ]] &&
        { echo 'The optional arguments -l (--projects-list) and -p (--mww-pvalue) / -a (--ea-table) are mutually exclusive. Choose one option or none, but not both'; exit; }

    # Check if only one of the way 3 variables has been provided
    [[( -z "$ea_table" || "$result_pvalue" -eq 1 ) && "$way" -eq 3 ]] &&
        { echo "You cannot provide only one of the arguments ea-table and mww-pvalue. You must provide both or none"; exit; }
}


#####################################################
#
#   This function separates a database (.fasta)
#   into two different files, one with all the
#   sequences belonging to a specific species and
#   the other with the sequences of the others of
#   the species in the database. To do this, it
#   receives the four-letter "id" of the species
#   in question, obtains the three-letter code of
#   that species establishing a correspondence
#   between both codes thanks to the file
#   'ids_table_file' and, finally, performs the
#   separation of the database file.
#
#   Arguments:
#       species id
#       Ids table file
#       Input file path
#       Output directory path
#
#####################################################

create_temporary_databases(){

    # Arguments
    local database="${1}"
    local species_id="${2}"
    local path_file_in="${3}"
    local path_dir_out="${4}"

    # Other variables
    file_out_species=$path_dir_out/$database"_species.fasta"
    file_out_others=$path_dir_out/$database"_others.fasta"
    path_tmp=$path_dir_out/tmp_databases

    # Create temporary directory
    mkdir -p $path_tmp
    
    # Get the required ids
    id_up=`echo ${species_id:0:1} | tr '[a-z]' '[A-Z]'`${species_id:1}

    # Replace uracils with thymines
    cat $path_file_in | cut -d " " -f 1 | awk 'BEGIN{RS=">"} {gsub("U", "T", $2); print ">"$1"\n"$2}' > $path_tmp/temp_file.fasta
    tail -n +3 $path_tmp/temp_file.fasta > $path_tmp/temp2_file.fasta

    # Separate database by species and others.
    cat $path_tmp/temp2_file.fasta | cut -d " " -f 1 | awk -v id="$species_id" \
                                              -v id_up="$id_up" \
                                              -v file_sp="$file_out_species" \
                                              -v file_o="$path_tmp/temp_others_file.fasta" \
                                              'BEGIN{RS=">"} {if ($0 ~ id || $0 ~ id_up) {print ">"$1"\n"$2 > file_sp } else {print ">"$1"\n"$2 > file_o}}'

    # Remove the first two faulty lines from the 'others' file.
    tail -n +3 $path_tmp/temp_others_file.fasta > $file_out_others

    # Remove temporary files
    rm -r $path_tmp
    
}


#####################################################
#
#   This function creates a fasta file from
#   the sequences located in the first column
#   of a results table obtained from a
#   differential expression analysis (DEA).
#   The name assigned to each sequence is
#   formed by the string 'sequence' followed
#   by a number that increments with each
#   sequence.
#
#   Arguments:
#       DEA results file path.
#       Output directory path.
#
#####################################################

table_to_fasta_and_tsv () {

    # Arguments
    local path_file="${1}"
    local path_dir_out="${2}"
    local path_file_out_id="${3}"

    # Sequence count
    seqnum=1

    # Read file
    while IFS=$'\t' read -r seq other
    do
        # Discard header
        if [[ $seq != "seq" ]]
        then
            # Create fasta file with sequences
            echo ">sequence"$seqnum >> $path_dir_out/$path_file_out_id.fasta
            echo $seq >> $path_dir_out/$path_file_out_id.fasta

            # Create tsv file with sequences
            printf "sequence"$seqnum"\t"$seq"\n" >> $path_dir_out/$path_file_out_id.tsv
            let "seqnum++" 
        fi
    done < $path_file

}


#####################################################
#
#   This function joins two TSV tables by
#   INNER JOIN or FULL OUTER JOIN, generating
#   an output file with the same format. It
#   receives the absolute paths of both TSV
#   files, the key column of each table to
#   perform the join, the columns of each
#   table to be saved in the output file,
#   the absolute path of the output file and
#   an argument specifying whether the
#   procedure will be done by INNER JOIN or
#   FULL OUTER JOIN.
#
#   Arguments:
#       Absolute path of TSV file 1
#       Absolute path of TSV file 1
#       Key column of TSV file 1
#       Key column of TSV file 2
#       Output file path
#       String specifying whether the join
#           will be made by INNER JOIN
#           (FALSE) or FULL OUTER JOIN (TRUE)
#           (Optional).
#
#####################################################

merge_tsv_files () {

     # Arguments
    local path_file_1="${1}"
    local path_file_2="${2}"
    local key_file_1="${3}"
    local key_file_2="${4}"
    local output_cols="${5}"
    local output_file="${6}"
    local outer_join="${7:-FALSE}"

    # Files name
    filename1_path=${path_file_1%.tsv}
    filename2_path=${path_file_2%.tsv}

    # Sort files to be joined
    LANG=C.UTF-8 sort -k $key_file_1 -t$'\t' $path_file_1 -o $filename1_path"_sort.tsv"
    LANG=C.UTF-8 sort -k $key_file_2 -t$'\t' $path_file_2 -o $filename2_path"_sort.tsv"
    
    # Inner join
    if [ $outer_join == "FALSE" ]
    then
        echo "LANG=C.UTF-8 join -1 $key_file_1 -2 $key_file_2 -t\$'\t' -o $output_cols -e \"NULL\" $filename1_path\"_sort.tsv\" $filename2_path\"_sort.tsv\" > $output_file"

        LANG=C.UTF-8 join -1 $key_file_1 -2 $key_file_2 -t$'\t' \
            -o $output_cols -e "NULL" \
            $filename1_path"_sort.tsv" \
            $filename2_path"_sort.tsv"  > $output_file
        

    # Full outer join
    elif [ $outer_join == "TRUE" ]
    then
        echo "LANG=C.UTF-8 join -1 $key_file_1 -2 $key_file_2 -t\$'\t' -o $output_cols -e \"NULL\" -a 1 -a 2 $filename1_path\"_sort.tsv\" $filename2_path\"_sort.tsv\" > $output_file"

        LANG=C.UTF-8 join -1 $key_file_1 -2 $key_file_2 -t$'\t' \
            -o $output_cols -e "NULL" \
            -a 1 -a 2 \
            $filename1_path"_sort.tsv" \
            $filename2_path"_sort.tsv"  > $output_file

    else
        echo "Error: Parameter incorrectly entered in merge_tsv_files function!"
        echo "Please specify whether the union of TSV files will be done by INNER JOIN (FALSE) or FULL OUTER JOIN (TRUE)."
    fi

}


sort_table_by_colnames() {
    # Arguments:
    # $1 -> path to the file (input file)
    # $2 -> output file name (without extension)
    # $3 -> Desired column order (array)

    # Read the file paths and desired order
    local path_in="${1}"
    local -n order="${2}"
    local path_out="${3}"

    # Read the header to get the index of each column
    read -r header < "$path_in"
    IFS=$'\t' read -r -a columns <<< "$header"

    # Create an array with the indices of the columns in the desired order
    indices=()
    for col in "${order[@]}"; do
        for i in "${!columns[@]}"; do
            if [[ "${columns[i]}" == "$col" ]]; then
                indices+=($i)
                break
            fi
        done
    done

    # Write the reordered header to the output file
    reordered_header=""
    len=${#indices[@]}
    for i in "${!indices[@]}"; do
        if [[ $i -eq $((len - 1)) ]]; then
            reordered_header+="${columns[${indices[i]}]}"
        else
            reordered_header+="${columns[${indices[i]}]}\t"
        fi
    done
    echo -e "$reordered_header" > "$path_out"

    # Reorder columns for each data line and write to the output file
    first_line=true
    while read -r line; do
        # Skip the header
        if $first_line; then
            first_line=false
            continue 
        fi

        IFS=$'\t' read -r -a data <<< "$line"
        reordered=""
        for i in "${!indices[@]}"; do
            if [[ $i -gt 0 ]]; then
                reordered+=$'\t'  # Add tab before all columns except the first
            fi
            reordered+="${data[${indices[i]}]}"
        done
        # Write the reordered line to the output file
        echo -e "$reordered" >> "$path_out"
    done < "$path_in"

}


#####################################################
#
#   This function creates an annotation table
#   from 6 SAM files obtained from the
#   alignment of a fasta file with the miRBase,
#   PmiREN and sRNAanno databases (sequences
#   belonging to the species and sequences
#   belonging to other species separately).
#   This table consists of 8 columns: seq, miRBase_
#   species, miRBase_others, PmiREN_species,
#   PmiREN_others, sRNAanno_species, sRNAanno_others
#   and length. The first column contains
#   the nucleotide sequence, the next 6 columns
#   represent the sequence name according to
#   miRBase (species and other species), PmiREN
#   (species and other species), and sRNAanno
#   (species and other species), respectively;
#   and the last column shows the length of the
#   nucleotide sequence in question.
#
#   Arguments:
#       Array with SAM files
#       Two-column TSV file with nucleotide
#           sequences of significantly
#           differentially expressed sequences
#           (column 2) and a name used for
#           sorting them (column 1) (e.g.,
#           sequence1 TTGAAAGTGACTACATCGGGG).
#       Output file name
#       Output directory path.
#
#####################################################

get_miRNAs_annotation () {

    # Arguments
    local -n array="${1}"
    local path_tsv_seq="${2}"
    local out_name="${3}"
    local path_out="${4}"

    # Create an array to know which files must be merged
    declare -A files_to_merge

    echo "IMPORTANTEEE"
    # Iterate through the array elements
    for db in "${!array[@]}"; do
        
        # Get the SAM file from the array
        sam_file=${array[$db]}

        # Get annotation tables
        [[ -f $sam_file ]] &&
            cat $sam_file | cut -f1,3 | awk -F "\t" '{if (substr($1,1,1) != "@"){print}}' > $path_out/$out_name"_"$db"_tmp.tsv" ||
            touch $path_out/$out_name"_"$db"_tmp.tsv"

        # Join nucleotide sequence with identifier
        merge_tsv_files $path_out/$out_name"_"$db"_tmp.tsv" $path_tsv_seq 1 1 1.2,2.2 $path_out/$out_name"_"$db"_annot.tsv"

        # Save the output file into the files_to_merge array
        files_to_merge["${db}"]=$path_out/$out_name"_"$db"_annot.tsv"

    done

    # Required vars
    output_columns="0,2.1"
    counter=1
    key_file1=2
    key_file_2=2
    final_colnames='1i\seq'

    # Iterate through the array elements (sorted!)
    for filename in $(echo "${!files_to_merge[@]}" | tr ' ' '\n' | sort -t_ -k1,1 -k2,2r); do

        # Get the file path
        current_file=${files_to_merge[$filename]}
        
        # Firs iteration
        if [[ $counter -eq 1 ]]; then
            merged_file=$current_file
            ((counter++))
            final_colnames+="\t${filename}"
            continue
        fi

        # Get the output columns depending on the merge step
        if [[ $counter -eq 3 ]]; then
            # Second merge
            output_columns="0,1.2,1.3,2.1"
        else
            next_column_to_add=$(echo $output_columns | awk -F, '{print $(NF-1)}' | awk '{if ($1 == 0) print 1.1; else print $1+0.1}' | bc)
            output_columns=$(echo $output_columns | sed 's/,[^,]*$//' | awk -v nc="$next_column_to_add" '{print $0 "," nc ",2.1"}')
        fi
        
        # Create the name of the new tmp file
        next_temp_file="$path_out/$out_name"_tmp_$counter.tsv

        # Merge tsv files
        merge_tsv_files \
            $merged_file \
            $current_file \
            $key_file1 \
            $key_file_2 \
            $output_columns \
            $next_temp_file \
            TRUE

        # Update some variables
        ((counter++))
        key_file1=1
        key_file_2=2
        merged_file=$next_temp_file
        final_colnames+="\t${filename}"
    done

    # # Sort the TSV file
    sort -k 1b,1 -t$'\t' $merged_file -o $path_out/$out_name"_annot_tmp.tsv"
    # #sed -i $final_colnames $path_out/$out_name"_annot.tsv"

    # If miRBase has not been provided...
    if [[ ! -v array["mirbase_species"] && ! -v array["mirbase_others"] ]]; then
        # Add the missing column names.
        final_colnames+="\tmirbase_species\tmirbase_others"
        # Add two NULL columns at the end of the table.
        awk 'BEGIN {OFS="\t"} {print $0, "NULL", "NULL"}' $path_out/$out_name"_annot_tmp.tsv" > tmp.csv && mv tmp.csv $path_out/$out_name"_annot_tmp.tsv"
    fi

    # If sRNAanno has not been provided...
    if [[ ! -v array["srnaanno_species"] && ! -v array["srnaanno_others"] ]]; then
        # Add the missing column names.
        final_colnames+="\tsrnaanno_species\tsrnaanno_others"
        # Add two NULL columns at the end of the table.
        awk 'BEGIN {OFS="\t"} {print $0, "NULL", "NULL"}' $path_out/$out_name"_annot_tmp.tsv" > tmp.csv && mv tmp.csv $path_out/$out_name"_annot_tmp.tsv"
    fi

    # If PmiREN has not been provided...
    if [[ ! -v array["pmiren_species"] && ! -v array["pmiren_others"] ]]; then
        # Add the missing column names.
        final_colnames+="\tpmiren_species\tpmiren_others"
        # Add two NULL columns at the end of the table.
        awk 'BEGIN {OFS="\t"} {print $0, "NULL", "NULL"}' $path_out/$out_name"_annot_tmp.tsv" > tmp.csv && mv tmp.csv $path_out/$out_name"_annot_tmp.tsv"
    fi
    
    # Sort the table
    sort -t$'\t' -k1,1 $path_out/$out_name"_annot_tmp.tsv" > tmp.csv && mv tmp.csv $path_out/$out_name"_annot_tmp.tsv"

    # Add the header to the annotation file
    sed -i $final_colnames $path_out/$out_name"_annot_tmp.tsv"

    ## Order the columns: miRBase, sRNAanno, PmiREN
    # Desired column order (by name)
    order=("seq" "mirbase_species" "mirbase_others" "srnaanno_species" "srnaanno_others" "pmiren_species" "pmiren_others")
    sort_table_by_colnames $path_out/$out_name"_annot_tmp.tsv" order $path_out/$out_name"_annot.tsv"
}


#####################################################
#
#   This function filters annotated sequences,
#   selecting only those that have been annotated
#   in at least two out of the three databases:
#   miRBase, PmiREN, and sRNAanno. It processes
#   the CSV file corresponding to the annotation
#   table containing the following fields: seq,
#   miRBase_species, miRBase_others, PmiREN_species,
#   PmiREN_others, sRNAanno_species, sRNAanno_others,
#   and length. The function checks each line in
#   the input file to determine if there are
#   annotations in at least two of the three
#   databases. Specifically, it checks if at
#   least one field in both miRBase (miRBase_species
#   and miRBase_others) and PmiREN (PmiREN_species
#   and PmiREN_others) or sRNAanno (sRNAanno_species
#   and sRNAanno_others) is not equal to "NULL".
#   If this condition is met, the function writes
#   the line to a new file (path_file_out). This
#   function effectively filters and saves lines
#   from the input CSV file based on the updated
#   criteria.
#
#   Arguments:
#       Input Annotation File Path.
#       Output Path for Filtered Annotation File.
#
#####################################################

filter_annotated_sequences () {

    # Arguments
    local path_file_in="${1}"
    local min_num_db="${2}"
    local path_file_out="${3}"

    # Read the first line to determine column positions
    IFS=$'\t' read -r header < "$path_file_in"

    # Convert the header into an array
    read -r -a columns <<< "$header"

    # Identify the indices of the relevant columns (excluding seq and length)
    declare -a db_columns=()
    for ((i = 1; i < ${#columns[@]} - 1; i++)); do
        db_columns+=("$i")
    done

    # Calculate the number of databases (each database has two columns: species & others)
    num_databases=$(( ${#db_columns[@]} / 2 ))

    # Print the header to the output file
    echo -e "$header" > "$path_file_out"

    # Read the original table and iterate through each line (skipping header)
    tail -n +2 "$path_file_in" | while IFS=$'\t' read -r -a line; do  
        # Number of databases where the sequence has been annotated
        points=0
        
        # Check each database (each has two associated columns)
        for ((i = 0; i < ${#db_columns[@]}; i+=2)); do
            species_col=${db_columns[$i]}
            others_col=${db_columns[$i+1]}

            # Check if the sequence is annotated at least once in this database
            [[ ( "${line[$species_col]}" != "NULL" || "${line[$others_col]}" != "NULL" ) ]] && ((points++))
        done
        
        # Print the line that meets the criteria and save it to the output file (using tab separator)
        if [ $points -ge $min_num_db ]; then
            printf "%s\n" "$(IFS=$'\t'; echo "${line[*]}")" >> "$path_file_out"
        fi

    done
}


#####################################################
#
#   This function takes the absolute path of the
#   summary table obtained from the Diff_exp_
#   analysis.r program and filters projects based
#   on a p-value threshold (MWW). It returns a
#   list of valid projects.
#
#   Arguments:
#       Input ea_summary.tsv file path.
#       MWW p-value threshold
#
#####################################################

filter_subprojects_by_mww_pvalue (){
    
    # Arguments
    local path_file_in="${1}"
    local mww_pvalue="${2}"

    # Filter sequences by pvalue (MWW)
    filtered_projects=$(tail -n +2 $path_file_in | awk -F '\t' -v var=$mww_pvalue '{if ($NF <= var) print $2}')
    
    # Return
    echo $filtered_projects

}


### MAIN
main () {

    # Get arguments
    arguments_management "$@"

    # Create an array with the provided databases
    declare -A databases=(
        ["mirbase"]="$path_mirbase"
        ["pmiren"]="$path_PmiREN"
        ["srnaanno"]="$path_sRNAanno"
    )

    # Create directory for temporary files
    mkdir -p ./tmp

    # Way 3: Get the valid subprojects from the provided list
    [ "$way" -eq 3 ] && valid_subprojects=$(filter_subprojects_by_mww_pvalue "$ea_table" "$mww_pvalue")
    
    # Create output directory
    file_name=$(basename "$path_mirbase") # hairpin.fa or mature.fa
    reference_name=$(echo "$file_name" | grep -oE "mature|hairpin")
    mkdir -p tmp/$reference_name

    ## Create temporary databases
    # Iterate through the array items (mirbase, pmiren, srnaanno)
    for db in "${!databases[@]}"; do

        # Check if the database has been provided
        if [ -n "${databases[$db]}" ]; then

            # Create the temporary database only if the key exists.
            create_temporary_databases "$db" "$species" "${databases[$db]}" "tmp/$reference_name"

            # Create Index directories
            mkdir -p tmp/$reference_name/"$db"_species_idx
            mkdir -p tmp/$reference_name/"$db"_others_idx

            # Index databases
            bowtie-build --threads $threads tmp/$reference_name/"${db}"_species.fasta tmp/$reference_name/"${db}"_species_idx/"${db}"_species > /dev/null 2>&1
            bowtie-build --threads $threads tmp/$reference_name/"$db"_others.fasta tmp/$reference_name/"$db"_others_idx/"$db"_others > /dev/null 2>&1
        fi
    done

    # Create an array to store the alignment results (SAM files).
    declare -A sam_array

    # Define the database types
    types=("species" "others")

    # Iterate through project files
    for file in $files
    do
        # Get data suffix (raw or sig)
        suffix=$(echo $file | awk -F '[_.]' '{print $(NF-1)}')

        # Way 2: Filtering using subprojects list
        if [ "$way" -eq 2 ]
        then
            # Check if subproject is valid
            check_subproject_valid=$( grep $file $projects_valid_list )

            # If it is not valid, pass to the next subproject
            [[ -z $check_subproject_valid ]] && valid="false" || valid="true"

            # Save results in summary file and pass to the next subproject.
            if [[ $valid == "false" ]]
            then
                echo -e "$species\t$id\tNA\tNA" >> tmp/$reference_name/$id"."$reference_name"_summary.tsv"
                continue
            fi

        # Way 3: Filtering using MWW p-value
        elif [ "$way" -eq 3 ]
        then
            # Delete suffix
            file_suffix_removed=$(echo "$file" | sed "s/\.dea_$suffix\.tsv$//" | cut -d'_' -f1-2)

            # If it is not valid, pass to the next subproject
            [[ $valid_subprojects != *$file_suffix_removed* ]] && valid="false" || valid="true"
                
            # Save results in summary file and pass to the next subproject.
            if [[ $valid == "false" ]]
            then
                echo -e "$species\t$id\tNA\tNA" >> tmp/$reference_name/$id"."$reference_name"_summary.tsv"
                continue
            fi
        fi

        printf "\n########################### File: $file ($species) ###########################\n\n"

        # Remove header and create temporary file
        tail -n +2 $file > tmp/$reference_name/$id"_temp_file.tsv"
        
        # Check if the file does not contain differentially expressed sRNAs (file empty)
        if [ -s "tmp/${reference_name}/${id}_temp_file.tsv" ]
        then

            # Create Fasta
            table_to_fasta_and_tsv $file tmp/$reference_name $id

            # Check if fasta file exist
            if test -f "tmp/${reference_name}/${id}.fasta"
            then

                # Iterate through databases names
                for db in "${!databases[@]}"; do

                    for type in "${types[@]}"; do

                        # Check if the database has been provided
                        if [ -n "${databases[$db]}" ]; then

                            # Create the required files names
                            index=tmp/$reference_name/"${db}_${type}_idx/${db}_${type}"
                            fasta_file=tmp/$reference_name/$id".fasta"
                            sam_file=tmp/$reference_name/$id"_${db}_${type}.sam"

                            # Alignment with miRBase (species databasw)
                            printf "Bowtie alignment with ${db} (${type} database)...\n"
                            bowtie -x $index --best -v $mismatches -k1 \
                                --no-unal -p $threads -f $fasta_file \
                                -S $sam_file 2>&1

                            # Check if an alignment error occurred due to the non-existence of the database.
                            [ $? -ne 0 ] && echo "Species $species not found in the ${db} database!" || printf "Done!\n"

                            # Save the sam files paths into the array
                            sam_array["${db}_${type}"]=$sam_file

                        fi
                    done
                done


                # Create annotation table from sam files
                printf "\nCreating annotation table...\n"
                get_miRNAs_annotation \
                    "sam_array" \
                    tmp/$reference_name/$id".tsv" \
                    $id \
                    tmp/$reference_name
                printf "Done!\n"
                
                # Calculate the sequences length and add a length column to annotated file
                awk 'BEGIN{ FS=OFS="\t" } {if (NR==1) {print $0, "length"} else {print $0, length($1)} }' tmp/$reference_name/$id"_annot.tsv" > $id".annot_all.tsv"
                #rm -r tmp/$reference_name/*_annot.tsv

                printf "\nCreating filtered annotation tables....\n"
                filter_annotated_sequences \
                    $id".annot_all.tsv" \
                    $min_num_db \
                    $id".annot_filt.tsv"
                printf "Done!\n"

                # Count the number of annotated sequences and how many of them have been selected.
                num_miRNAs=$(tail -n +2 $id".annot_all.tsv" | wc -l)
                num_miRNAs_filtered=$(tail -n +2 $id".annot_filt.tsv" | wc -l)

                ## Get the number of annotated sequences for each length (20-25)
                awk -F '\t' -v species="$species" -v project="$id" 'NR>1 {count[$8]++} END {printf "%s\t%s\t", species, project; for (i=20; i<=25; i++) printf "%s%s", count[i] ? count[i] : 0, (i<25) ? "\t" : ""; printf "\n"}' $id".annot_all.tsv" >> tmp/$reference_name/summary_len_"$reference_name"_tmp.tsv
                awk -F '\t' -v species="$species" -v project="$id" 'NR>1 {count[$8]++} END {printf "%s\t%s\t", species, project; for (i=20; i<=25; i++) printf "%s%s", count[i] ? count[i] : 0, (i<25) ? "\t" : ""; printf "\n"}' $id".annot_filt.tsv" >> tmp/$reference_name/summary_len_"$reference_name"_filt_tmp.tsv
                merge_tsv_files tmp/$reference_name/summary_len_"$reference_name"_tmp.tsv tmp/$reference_name/summary_len_"$reference_name"_filt_tmp.tsv 2 2 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.3,2.4,2.5,2.6,2.7,2.8 $id"."$reference_name"_summary_len.tsv" 
                sed -i '1i\Species_id\tGroup_comparison\t20nt_all\t21nt_all\t22nt_all\t23nt_all\t24nt_all\t25nt_all\t20nt_filt\t21nt_filt\t22nt_filt\t23nt_filt\t24nt_filt\t25nt_filt' $id"."$reference_name"_summary_len.tsv"
                #rm ./tmp/*tmp.tsv

                # Save it in summary file
                echo -e "$species\t$id\t$num_miRNAs\t$num_miRNAs_filtered" >> $id"."$reference_name"_summary.tsv"
                sed -i '1i\Species_id\tGroup_comparison\tNum_seq_all\tNum_seq_filt' $id"."$reference_name"_summary.tsv"
                
                
            else
                printf "NOTE: $id"".dea_"$suffix".fasta does not exist\n"
                echo -e "$species\t$id\t0\t0" >> $id"."$reference_name"_summary.tsv"
            fi

        # The file has no differentially expressed sRNAs. 
        else
            printf "NOTE: the file has no differentially expressed sRNAs (File is empty)\n"
            echo -e "$species\t$id\t0\t0" >> $id"."$reference_name"_summary.tsv"
        fi
    done

    #rm -r tmp

}
main "$@"
