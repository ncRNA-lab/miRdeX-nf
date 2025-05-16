#!/usr/bin/env Rscript

################################################################################
##                                                                            
##  Annotate_dea_results.r
##
##  This program annotates the output tables from DESeq2, linking each sequence 
##  identified as a miRNA to its corresponding annotation based on a GFF3 file. 
##  It generates two output tables: one containing all possible annotations 
##  for each sequence (accounting for cases where a single sequence could be 
##  an isomiR of multiple reference miRNAs), and another table containing 
##  unique miRNA sequences with their associated isomiR type.
##
##  The script also checks whether all the sequences belonging to the same 
##  miRNA family exhibit the same differential expression trend. It produces 
##  a summary table with this information, and a boxplot to visually explore 
##  the expression patterns across families.
##                                                                            
##                                                                            
##  Author: Antonio Gonzalez Sanchez                                         
##  Date: 05/014/2025
##  Version: 1.0
##                                                                            
################################################################################


################################################################################
################################## MODULES #####################################
################################################################################

suppressMessages(library(argparse))
suppressMessages(library(ape))
suppressMessages(library(tidyverse))

################################################################################
################################## FUNCTIONS ###################################
################################################################################


#' Get the command line arguments
#' 
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values
#' 

getArguments <- function(){
  
  # create parser object
  parser <- ArgumentParser(prog='Group_miRNAs_by_family.r',
                           description= '
   This program groups the miRNAs identified in families to study whether all
   the miRNAs belonging to the same family follow the same trend in terms of
   differential expression using boxplots. In addition, it generates a summary
   file that lists those miRNA families whose members do not show the same
   trend.',
                           formatter_class= 'argparse.RawTextHelpFormatter')
  
  required = parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option
  required$add_argument('-i', '--id',
                        nargs = 1,
                        type = 'character',
                        help = 'Job identifier',
                        required = TRUE)
  required$add_argument('-d', '--dea',
                        nargs = "*",
                        type = 'character',
                        help = 'Result files from the differential expression analysis',
                        required = TRUE)
  required$add_argument('-a', '--annotation',
                        nargs = "*",
                        type = 'character',
                        help = 'miRNA annotation gff3 files',
                        required = TRUE)
  parser$add_argument('-c', '--classes',
                        nargs = 1,
                        type = 'character',
                        help = 'Comma-separated list of miRNA classes to include (e.g., ref_miRNA,iso_5p)',
                        required = FALSE)
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('id', 'dea', 'annotation')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  # If the --classes argument is specified, convert it to a vector
  if (!is.null(args$classes)) {
    args$selected_classes <- unlist(strsplit(args$classes, ","))
  } else {
    args$selected_classes <- NULL
  }
  
  return(args)
}


#' Read and preprocess a miRNA GFF3 file
#'
#' This function reads a GFF3 file containing miRNA annotations and processes
#' the `attributes` field by converting it into separate columns. The resulting
#' data frame  contains the original GFF3 fields along with the extracted
#' attributes for easier access.
#'
#' @param gff3_path Path to the GFF3 file to be read.
#'
#' @return A data frame containing the GFF3 data, with each attribute moved to
#'         its own column.
#'
#' @examples
#' gff3_df <- read_miRNA_gff3("path/to/miRNA_annotation.gff3")

read_miRNA_gff3 <- function(gff3_path) {
  
  # Read the miRNA gff3 file
  gff3 <- read.gff(gff3_path, GFF3 = TRUE)
  
  # Move attributes to columns
  gff3_df <- gff3 %>%
    select(-phase) %>% 
    separate_rows(attributes, sep = ";\\s*") %>%
    separate(attributes, into = c("attr_name", "attr_value"), sep = "=", fill = "right") %>%
    pivot_wider(names_from = attr_name, values_from = attr_value)
  
  return(gff3_df)
}


#' Create size-adjusted boxplot
#' 
#' This function creates a boxplot by separating it into different panels
#' if the number of labels on the x-axis exceeds a certain threshold, thus
#' avoiding an overlapping of the labels or a bad display of the plot.
#' 
#' @param data A dataframe
#' @param x Column to be set on the x-axis
#' @param y Column to be set on the y-axis
#' @param max_labels X-axis label threshold
#' @param x_lab X label
#' @param y_lab Y label
#' @param z Column used to set the colors of the points.
#' @param legend_title Title of the legend. This argument will only be used if
#'                     the 'z' argument has also been provided.
#' @return  Ggplot2 boxplot
#' @examples 
#' createBoxplot(df, 'Column1', 'Column2', 40,'miRNAs','Slog2FC')

createBoxplot <- function(data, x, y, max_labels, x_lab, y_lab, z, legend_title) {
  
  # Get number of labels
  num_labels <- length(unique(data[[x]]))
  
  # Split data into subsets for plotting
  if (num_labels > max_labels) {
    n_splits <- ceiling(num_labels / max_labels)
    split_labels <- cut(as.numeric(factor(data[[x]], levels = unique(data[[x]]))), breaks = n_splits, labels = FALSE)
    data$split <- split_labels
  } else {
    data$split <- 1
  }
  
  # Create plot
  p <- ggplot(data, aes(x=!!sym(x), !!sym(y))) +
    geom_boxplot(color='#333333', fill='#FFFAFA') +
    scale_y_continuous(breaks = round(seq(min(data[[y]]), max(data[[y]]), by = 2))) +
    labs(x = x_lab, y = y_lab) +
    theme_linedraw() +
    theme(axis.text.x = element_text(angle=90),
          axis.title.y = element_text(margin = margin(r = 10)),
          strip.text = element_blank()) +
    geom_hline(yintercept=0, color = '#333333') +
    facet_wrap(~ split, ncol = 1, scales = "free_x")
  
  # Add optional options
  if (!missing(z)) {
    # Add colors to points
    p <- p +
      geom_point(aes(colour = as.factor(.data[[z]])), size=2, show.legend = TRUE) +
      scale_color_viridis_d(option = "D")
    # Add legend title
    if (!missing(legend_title)){
      p <- p + labs(colour = legend_title)
    }
  }
  
  
  
  # Return plot
  return(p)
}


################################################################################
##################################### MAIN #####################################
################################################################################

# Get program arguments
args <- getArguments()

# Save the rest of the arguments in variables
id <- args$id
dea_path <- args$dea
gff3_path <- args$annotation
classes <- args$selected_classes

### 2. CREATE THE ANNOTATED DEA RESULTS DATAFRAMES
################################################################################

# Read the input files
dea_df <- read.table(dea_path, header = TRUE)
gff3_df <- read_miRNA_gff3(gff3_path)

# Remove those sequences whose class is not in classes
if (!is.null(classes)) {
  classes_v <- unlist(strsplit(classes, ","))
  gff3_df <- gff3_df[gff3_df$Class %in% classes_v, ]
}

# Check if all the duplications has the same isomiR class
gff3_check_df <- gff3_df %>%
  group_by(Read) %>%
  mutate(
    Class_check = if (n_distinct(Class, na.rm = TRUE) == 1) first(Class) else "Undefined"
  ) %>%
  ungroup()

# Merge DEA results dataframe with gff3 dataframe (Keep duplicated sequences (IsomiRs from multiple reference miRNAs)
dea_annotated_all_isomirs <- merge(dea_df, gff3_check_df, by.x = 'seq', by.y = 'Read', all.y = TRUE)

# Save table
file_name_table_all = paste0(id, '.all.tsv')
write.table(dea_annotated_all_isomirs, file = file_name_table_all, sep = "\t", quote = FALSE, row.names = FALSE)

# Remove duplicated sequences 
gff3_uniq <- gff3_check_df %>%
  select(Read, UID, miRNA_fam, Class_check) %>%
  distinct()

# Merge DEA results dataframe with gff3 dataframe
dea_annotated_uniq <- merge(dea_df, gff3_uniq, by.x = 'seq', by.y = 'Read')

# Save table
file_name_table_unique = paste0(id, '.unique.tsv')
write.table(dea_annotated_uniq, file = file_name_table_unique, sep = "\t", quote = FALSE, row.names = FALSE)

### 3. CREATE BOXPLOT
################################################################################

# Create output file name
file_name_plot = paste0(id, '.boxplot.png')

# Create the boxplot
p <- createBoxplot(dea_annotated_uniq, 'miRNA_fam', 'Shrunkenlog2FoldChange', 40, '', 'Log2FC', z='Class_check', legend_title = 'isomiR class')

# Save plot in output directory
ggsave(file_name_plot, p)


### 4. TABLE WITH DIFFUSE-TREND MIRNAS FAMILY
################################################################################

# Get miRNA family names
miRNAs_v <- unique(dea_annotated_uniq$miRNA_fam)

# Select miRNAs families with diffuse trend
miRNAs_var <- c()
for (miRNA in miRNAs_v){
  
  # Get the miRNA Shrunkenlog2FoldChange vector
  lfc_v <- na.omit(dea_annotated_uniq[dea_annotated_uniq$miRNA_fam == miRNA,]$Shrunkenlog2FoldChange)
  
  # If there are positives and negatives
  if (any(lfc_v > 0) && any(lfc_v < 0)) {
    miRNAs_var <- c(miRNAs_var, miRNA)
  }
}

# If there are no miRNAs families with diffuse trend...
if (is.null(miRNAs_var)){
  miRNAs_var <- "NULL"
  num_miRNAs_d <- 0
} else {
  num_miRNAs_d <- length(miRNAs_var)
}

# Save it in summary file
miRNAs_txt <- paste(miRNAs_var, collapse = '/')
text <- paste(id, length(miRNAs_v), num_miRNAs_d, miRNAs_txt, sep='\t')
cat(text, file=paste0(id,'.summary.tsv'), append=TRUE, sep='\n')


