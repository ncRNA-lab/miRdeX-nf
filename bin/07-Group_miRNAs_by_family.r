#!/usr/bin/env Rscript

################################################################################
##                                                                            
##  Group_miRNAs_by_family.r
##
##  This program groups the miRNAs identified in families to study whether     
##  all the miRNAs belonging to the same family follow the same trend in
##  terms of differential expression using boxplots. This program also
##  generates a new table where the results of the differential expression
##  analysis are put together with the annotation results, including the
##  family to which each miRNA belongs.  It also generates a summary file
##  that lists those miRNA families whose members do not show the same trend.                                                
##                                                                            
##                                                                            
##  Author: Antonio Gonzalez Sanchez                                         
##  Date: 12/20/2023
##  Version: 2.0
##                                                                            
################################################################################


################################### MODULES ####################################

suppressMessages(library(argparse))
suppressMessages(library(tidyverse))
#suppressMessages(library(dplyr)) EStos tres paquetes están en tidyverse
#suppressMessages(library(stringr))
#suppressMessages(library(ggplot2))


################################## FUNCTIONS ###################################

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
                        help = 'miRNA annotation files',
                        required = TRUE)
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('id', 'dea', 'annotation')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  return(args)
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
    p <- p + geom_point(aes(colour = as.factor(.data[[z]])), size=0.5, show.legend = TRUE) +
      scale_color_manual(values = c("#CB0000", "blue"))
    # Add legend title
    if (!missing(legend_title)){
      p <- p + labs(colour = legend_title)
    }
  }

  
  # Return plot
  return(p)
}


##################################### MAIN #####################################

# Get program arguments
args <- getArguments()

# Save the rest of the arguments in variables
id <- args$id
files_dea_list <- args$dea
files_annot_list <- args$annotation

# Iterate files
for (i in 1:length(files_annot_list)){
  
  ### 1. CREATE FAMILY MIRNAS NAMES
  ##############################################################################

  # Get files
  file_annot <- files_annot_list[i]
  file_dea <- files_dea_list[i]
  
  # Read files
  annot_table <- read.csv2(file_annot, sep = '\t')
  dea_table <- read.csv2(file_dea, sep = '\t')
  
  # If dataframe is not empty
  if (nrow(annot_table) != 0){
    # Reorder columns and replace NUll with NA
    annot_table <- annot_table[, c(1,2,4,6,3,5,7,8)] %>% replace(.=='NULL', NA)
    
    # Iterate rows
    fam_id_v = c()
    for (row in 1:nrow(annot_table)){
      # Iterate columns
      for (col in 2:ncol(annot_table)){
        # Check if the sequence is annotated
        if (!is.na(annot_table[row,col])){
          annotation <- annot_table[row,col]
          break
        }
      }
      # Remove species id and '-'
      annot_elements <- strsplit(annotation, '-', fixed =TRUE)
      if (length(annot_elements[[1]]) > 3){
        annot_wsp <- paste0(annot_elements[[1]][c(-1,-4)], collapse ='')
      } else {
        # Detects the presence of a number in the string
        check <- suppressWarnings(min(which(!is.na(as.numeric(str_split(annot_elements[[1]][2], '', simplify = TRUE))))))
        if (check == "Inf"){
          # The identifier number is in element 3.
          annot_wsp <- paste0(annot_elements[[1]][-1], collapse ='')
        } else {
          # The identifier number is in element 2, delete element 3.
          annot_wsp <- paste0(annot_elements[[1]][c(-1,-3)], collapse ='')
        }
      }
      
      # Detects the position of the first numeric value within the string
      pos_num <- suppressWarnings(min(which(!is.na(as.numeric(str_split(annot_wsp, '', simplify = TRUE))))))
      
      if ( pos_num == 'In'){
        fam_id_v = c(fam_id_v, annot_wsp)
      } else {
        
        # Create a new string from this position
        subs <- str_sub(annot_wsp, pos_num, -1)
        
        # Detects the position of the first non-numeric value within the new string
        pos_let <- suppressWarnings(min(which(is.na(as.numeric(str_split(subs, "", simplify = TRUE))))))
        
        # Create the miRNA family id
        if ( pos_let == 'Inf'){
          fam_id_v = c(fam_id_v, annot_wsp)
        } else {
          fam_id <- str_sub(annot_wsp, 0, pos_let + pos_num - 2)
          fam_id_v = c(fam_id_v, fam_id)
        }
      }
    }
    
    # Add family names
    annot_table <- cbind(annot_table, general_annot = fam_id_v)
    
    ### 2. ADD ANNOTATION TO THE DEA-TABLE
    ######################################################################
    
    # Create final table
    final_table = merge(dea_table, annot_table, by = 'seq')
    final_table$Shrunkenlog2FoldChange <- as.numeric(as.character(final_table$Shrunkenlog2FoldChange))
    
    # Save table
    file_name_table = paste0(id, '.fam_annot.tsv')
    write.table(final_table, file = file_name_table, sep = "\t", quote = FALSE, row.names = FALSE)
	    
    ### 3. CREATE THE PLOT
    ######################################################################
    
    # Create boxplot
    final_table$length_condition <- ifelse(final_table$length <= 22, "<= 22", "> 22")
    p <- createBoxplot(final_table, 'general_annot', 'Shrunkenlog2FoldChange', 40,'','Slog2FC', z='length_condition', legend_title = 'Sequence length')

    # Create output file name
    file_name_plot = paste0(id, '.fam_annot_boxplot.png')
  
    # Save plot in output directory
    ggsave(file_name_plot, p)

  
    ### 4. TABLE WITH DIFFUSE-TREND MIRNAS FAMILY
    ######################################################################
  
    # Get miRNA family names
    miRNAs_v <- unique(final_table$general_annot)
  
    # Select miRNAs families with diffuse trend
    miRNAs_var <- c()
    for (miRNA in miRNAs_v){
    
      # Get the miRNA Shrunkenlog2FoldChange vector
      lfc_v <- na.omit(final_table[final_table$general_annot == miRNA,]$Shrunkenlog2FoldChange)
    
      # If there are positives and negatives
      if (any(lfc_v > 0) && any(lfc_v < 0)) {
        miRNAs_var <- c(miRNAs_var, miRNA)
      }
    }
  
    # If there are no miRNAs families with diffuse trend...
    if (length(miRNAs_var) == 0){
      miRNAs_var <- 'NULL'
      num_miRNAs_d <- 0
    } else {
      num_miRNAs_d <- length(miRNAs_var)
    }
  
    # Save it in summary file
    miRNAs_txt <- paste(miRNAs_var, collapse = '/')
    text <- paste(id, length(miRNAs_v), num_miRNAs_d, miRNAs_txt, sep='\t')
    cat(text, file=paste0(id,'.annot_fam_summary.tsv'), append=TRUE, sep='\n')
  }
}

