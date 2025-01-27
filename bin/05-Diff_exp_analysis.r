#!/usr/bin/env Rscript

################################################################################
################################## MODULES #####################################
################################################################################

suppressMessages(library(argparse))
suppressMessages(library(dendextend))
suppressMessages(library(DESeq2))
suppressMessages(library(ff))
suppressMessages(library(htmltools))
suppressMessages(library(plotly))
suppressMessages(library(RColorBrewer))
suppressMessages(library(readr))
suppressMessages(library(SARTools))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(tidyverse))
suppressMessages(library(ggthemes))
suppressMessages(library(DEGreport))
suppressMessages(library(paletteer))
suppressMessages(library(gtools))

# Instlacion DEGreport para la imagen
#install.packages("https://cran.r-project.org/src/contrib/Archive/lasso2/lasso2_1.2-22.tar.gz", type = "source", repos = NULL)
#BiocManager::install("DEGreport")


################################################################################
################################# FUNCTIONS ####################################
################################################################################


#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = '05-Diff_exp_analysis.r',
                           description = '
    This program takes the tables of absolute counts from each project and
    1. Exploratory analysis
    
    This program takes the tables of absolute counts from a project and
    performs a Principal Component Analysis (PCA) for each of the stress
    events considered in that project. From the results of this analysis,
    it takes the coordinates generated for each sample from the values of the
    first three principal components and calculates the Euclidean distances
    between samples of the same condition or INTRA-group (e.g. treated1 - 
    treated2) and the distances between samples of different conditions or
    INTER-group (e.g. treated1 - control1). Then, it performs a Mann-whitney-
    wilcoxon test to check if there are differences between the INTRA and
    INTER-group distances previously calculated. Additionally, creates a
    dendrogram using the Simple Error Rate Estimation (SERE) obtained from
    comparing the different samples of the absolute counts table and a plot
    of mean vs variance comparison.
     
    2. Differential expression analysis
     
    Then, the program performs a differential expression analysis using
    DESeq2. The absolute counts tables contain a group of control samples and
    different treatment samples to which they are related. The differential
    expression analysis is performed considering all possible combinations of
    control vs treated (c_vs_t1, c_vs_t2, etc), so the program returns a result
    table for each of them. The results table contains all the information
    provided by the results() function of DESeq2 together with the
    log2FoldChange and lfcSE from lfcShrink. In addition to the raw data
    obtained in the analysis, this script also provides tables with those
    sequences with an adjusted p-value lower than 0.05.',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--id',
                        type = 'character',
                        help = 'Job identifier.',
                        required = TRUE)
  required$add_argument('-c', '--counts',
                        type = 'character',
                        help = 'Absolute counts matrix file.',
                        required = TRUE)
  required$add_argument('-m', '--metadata',
                        type = 'character',
                        help = 'Metadata directory path.',
                        required = TRUE)
  parser$add_argument('-a', '--alpha',
                      type = 'double',
                      default = 0.05,
                      help = 'Alpha significance level  (default: 0.05)')
  parser$add_argument('-f', '--min_counts',
                      type = "integer",
                      help = 'Minimum count threshold per sample (default: 5)',
                      default = 5)
  parser$add_argument('-s', '--min_samples',
                      type = "integer",
                      help = 'Minimum number of samples that must meet the threshold (default: 5)',
                      default = 5)
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('id', 'counts', 'metadata', 'alpha', 'min_counts', 'min_samples')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  # Check if the input directory exists
  if (!file.exists(args$counts)) {
    stop('Error. The input counts matrix does not exist.')
  }
  
  return(args)
}


#' Filter and Prepare DESeq2 Dataset
#'
#' This function takes as input a count matrix file path, a metadata file path,
#' and two filtering thresholds: the minimum count value (`min_counts`) and the
#' minimum number of samples (`min_samples`) that must meet the count threshold.
#' The function first filters the count matrix to remove sequences (rows) with
#' low counts. Then, it creates a DESeqDataSet object using the filtered count
#' matrix and metadata.
#'
#' @param counts_path A character string representing the file path to the
#'                    count matrix (TSV format).
#' @param metadata_path A character string representing the file path to the
#'                      metadata file associated with the samples (TSV format).
#' @param min_counts An integer specifying the minimum count threshold for a
#'                    sequence to not be removed.
#' @param min_samples An integer specifying the minimum number of samples
#'                    that must meet the `min_counts` threshold for a sequence
#'                    to not be removed.
#' 
#' @return A DESeqDataSet object
#' 
#' @examples
#' # Assuming 'counts.tsv' is your count matrix file and 'metadata.tsv' is
#' # your metadata file with minimum count of 5 and minimum samples of 3.
#' dds <- create_DeseqDataSet("counts.tsv", "metadata.tsv", min_counts = 5, min_samples = 3)
#' 
 
create_DeseqDataSet <- function(id, file, metadata, min_counts=5, min_samples=5){
  
  ###################### Prepare the counts matrix #############################
  
  # Read the absolute counts file using fread
  counts_tb <- read_tsv(file, col_names = TRUE)
  
  # The df is not empty
  if (nrow(counts_tb) > 0) {
    
    # Move th seq column to the rownames
    counts_df <- column_to_rownames(counts_tb, "seq")
    
    # Filter the count matrix by low counts
    if (ncol(counts_df) < min_samples) {
      # If there are fewer samples than min_samples, require all samples to meet the threshold
      counts_df_filt <- counts_df[rowSums(counts_df >= min_counts) == ncol(counts_df), ]
    } else {
      # Otherwise, apply the standard filter
      counts_df_filt <- counts_df[rowSums(counts_df >= min_counts) >= min_samples, ]
    }
  
    ############# Create a default ColData using the metadata file ###############
  
    # Read the metadata file
    metadata_df <- read.csv(metadata, sep = "\t", header = TRUE)
    
    # Select only the runs present in the counts matrix
    metadata_df <- metadata_df[metadata_df$Run %in% colnames(counts_df), ]
    metadata_df <- metadata_df[match(colnames(counts_df), metadata_df$Run), ]
    
    # Apply trimws to each element of the dataframe (strip)
    metadata_df <- as.data.frame(lapply(metadata_df, function(col) {
      if (is.character(col)) { return(trimws(col)) } else { return(col) }
    }), stringsAsFactors = FALSE)
    
    # Get the group of samples
    group_of_samples <- str_split(id, "_")[[1]][2]
    
    # Get the metadata of the group of samples
    group_metadata <- metadata_df[metadata_df$Group == group_of_samples,]
    print(group_of_samples)
    print(group_metadata)
    
    # Get the design formula
    design_formula <- unique(group_metadata$Design)
    
    # Create a default ColdData dataframe
    coldata <- data.frame(
      Group = group_metadata$Group,
      species = group_metadata$Species,
      Project = group_metadata$Project,
      Run = group_metadata$Run,
      Replicate = group_metadata$Replicate, 
      Test = group_metadata$Test,
      Design = group_metadata$Design,
      row.names = 'Run')
    
    # Add the reduced formula if the test selected is LRT
    if (unique(group_metadata$Test) == "LRT") {
      coldata$DesignRed <- group_metadata$DesignRed
    }
    
    # Add the contrast if it is specified
    if (unique(group_metadata$Contrast) != "CT.0") {
      coldata$Contrast<- group_metadata$Contrast
    }
    
    ##################### Add the factors to the colData #########################
  
    # Obtain the factors involved in the design.
    design_factors <- strsplit(unique(group_metadata$DesignRef), ":")[[1]]
    
    # Iterate through factors
    factors_reference_list <- c()
    for (factor_info in design_factors) {
      # Get the factor name and its reference
      factor_name <- strsplit(factor_info, "[(|)]")[[1]][1]
      factor_ref <- strsplit(factor_info, "[(|)]")[[1]][2]
      
      # Add both elements to the list
      factors_reference_list[[factor_name]] <- factor_ref
      
      # Add the factor values in the coldata table (as factor)
      coldata[[factor_name]] <- as.factor(group_metadata[[factor_name]])
    }
    
    ######################### Create the DESeqDataSet ############################
  
    print(colnames(as.matrix(counts_df_filt)))
    print(coldata)
    
    # Create DESeqDataSet
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = as.matrix(counts_df_filt),
                                                   colData = coldata,
                                                   design = as.formula(design_formula)))
    
    
    ################## Set the reference level for each factor ###################
  
    # Iterate through factors
    for (factor_name in names(factors_reference_list)) {
      # Set the reference values as reference for the rest of the factors
      dds[[factor_name]] <- relevel(dds[[factor_name]], ref = factors_reference_list[[factor_name]])
    }
  
  # The counts matrix is
  } else{
    dds <- -1
  }
  
  # Return the DeseqDataSet
  return(dds)
}

#' Row variance
#' 
#' This function calculates the variance of each row in a matrix or data frame
#' "a" using the R function "apply".
#'
#' @param a A DataFrame or Matrix object
#' @return A vector with the variances of each row of the DataFrame or Matrix.
#' @examples
#' row_var(df)

row_var <- function(a) {
  apply(a, 1, var)
}

#' Euclidean distance
#' 
#' This function calculates the Euclidean distance between two points using
#' the formula: \code{sqrt(sum((a - b) ^ 2))}
#'
#' @param a A numeric vector
#' @param b A numeric vector
#' @return Euclidean distance between a and b
#' @examples
#' euclidean_dist(c(3, 15, 4), c(23,11,8))

euclidean_dist <- function(a, b) {
  return(sqrt(sum((a - b) ^ 2)))
}


#' Creates a dendrogram using the Simple Error Rate Estimation (SERE) obtained
#' from comparing the different samples of a DESeqDataSet
#' 
#' This function recieves a DESeq Dataset and uses the counts matrix to calculate
#' the Simple Error Ratio Estimate (SERE) for each sample comparison within the
#' dataset, a statistic that can determine whether two RNA-seq libraries are
#' faithful replicates or globally different (Schulze, Kanwar & Gölzenleuchter,
#' 2012). With the SERE obtained from these comparisons, it creates a distance
#' matrix and generates a dendrogram, which it saves in the output directory
#' together with the node points associated with the dendrogram (points at which
#' a branch is bifurcates) and the SERE values obtained (matrix).
#'
#' @param dds DESeqDataSet
#' @param path_dir_out Output directory path
#' @return  A Matrix with the SERE values of each comparison
#' @examples 
#' csv_to_deseq_dataset("/home/minimind/Desktop/Results/arth/PRJNA277424_1.tsv")
#' @references
#' Schulze, S. K., Kanwar, R., Gölzenleuchter, M., Therneau, T. M., & Beutler,
#' A. S. (2012). SERE: single-parameter quality control and sample comparison
#' for RNA-Seq. BMC genomics, 13, 524. https://doi.org/10.1186/1471-2164-13-524

sere_dendrogram <- function(dds, path_dir_out){
  
  ###################### PREPARE THE SAMPLE NAMES ##############################
  # Reduce the sample name if any of its fields exceeds 40 characters.
  # Get the counts matrix
  df <- counts(dds)
  
  # Get the colData
  coldata_df <- as.data.frame(colData(dds))
  
  # Get the design
  design_elements <- str_split(as.character(design(dds))[2], " \\+ ")[[1]]
  
  # Remove the interaction if it exists
  design_elements <- design_elements[!grepl(":", design_elements)]
  
  # Add a new column joining the factors of the design formula (Sample name. Factors + replicate)
  coldata_df$Sample_name <- apply(coldata_df[c(design_elements, "Replicate")], 1, paste, collapse = "_")
  
  # Add a new column joining the factors of the design formula (Condition name. Factors)
  coldata_df$Condition_name <- apply(coldata_df[design_elements], 1, paste, collapse = "_")
  
  ################### BUILD THE MATRIX WITH SERE VALUES ########################
  
  # Create an empty matrix
  sere_dist <- matrix(NA, nrow = ncol(df), ncol = ncol(df))
  colnames(sere_dist) <- coldata_df$Sample_name
  rownames(sere_dist) <- coldata_df$Sample_name
  
  # Calculate SERE coefficent for each sample comparison
  for (i in 1:(ncol(df) - 1)){
    for (j in (i + 1):ncol(df)){
      two_samples_counts <- cbind(df[, i], df[, j])
      sere_coefficient <- SERE(two_samples_counts)
      sere_dist[i, j] <- sere_coefficient
      sere_dist[j, i] <- sere_coefficient
    }
  }
  
  # Save sere_dist object
  write.table(sere_dist, file = paste0(path_dir_out,'/01-SERE_matrix.tsv'), sep = "\t", quote = FALSE)
  
  ############################ BUILD THE DENDROGRAM ############################

  ## 3.1. Build the dendrogram
  dend <- as.dist(sere_dist) %>% hclust %>% as.dendrogram
  
  # Save node points
  write.table(hclust(as.dist(sere_dist))$height,
              file = paste0(path_dir_out,'/02-SERE_dendrogram_node_points.tsv'),
              quote = FALSE,
              row.names = FALSE,
              col.names = "Node points", sep = ",")
  
  ## 3.2. Let's add some color
  # Get unique names
  unique_names <- unique(coldata_df$Condition_name)
  
  # Check if i only have two levels
  if (length(unique_names) == 2) {
    # Use the 'Set1' palette for two levels
    color_mapping <- setNames(c("#1B9E77", "#D95F02"), unique_names)
  } else {
    # Use the 'Dark2' palette for more than two levels
    color_mapping <- setNames(brewer.pal(length(unique_names), 'Dark2'), unique_names)
  }
  
  # Create a new vector of colors based on the original vector
  assigned_colors <- color_mapping[coldata_df$Condition_name]
  
  # But sort them based on their order in dend:
  colors_to_use <- assigned_colors[order.dendrogram(dend)]
  
  ## 3.3. Add labels and colors
  labels(dend) <- coldata_df$Sample_name[order.dendrogram(dend)]
  labels_colors(dend) <- colors_to_use
  dend <- dend %>% set("labels_cex", 0.5)
  
  ## 4.4. Create and save the plot
  pdf(paste(path_dir_out, '/03-Dendrogram.pdf', sep = ''))
  right_m <- round(max(nchar(coldata_df$Sample_name))/4)
  par(mar = c(4, 4, 2, right_m))
  plot(dend, horiz= TRUE)
  abline(v = 1, col = "#666666", lty = 2)
  dev.off()
  
  return(sere_dist)
}


#' Exploratory Analysis
#' 
#' This function receives a DESeqDataSet object and performs an individual
#' exploratory analysis for each of the stress events it contains. In this
#' analysis, a variance stabilisation transformation (VST) is performed on the
#' absolute counts data to subsequently perform a principal component analysis
#' (PCA). To discriminate whether the samples are grouped by conditions
#' (control, treated...), the Mann-Whitney-Wilcoxon test (MWW) is performed to
#' check i there are significant differences between intra-group (replicates of
#' the same condition e.g. control1 vs control2) and inter-group (samples that
#' are associated to different conditions e.g. control vs treated) distances.
#' The function returns a vector whose first 6 elements indicate the variance
#' explained by each of the first 6 principal components and the last element
#' represents the p-value obtained in the MWW test.
#'
#' @param dds A DESeqDataSet object
#' @param path_out Path of the output directory to save the generated plots
#'                 (PCA and barplot of the variance explained by each component).
#' @param file_name Name for output files
#' @return  A numeric vector with the variance explained by the first 6
#'          principal components and de p-value obtained in the MWW test.
#' @examples
#' euclidean_dist(dds,/home/minimind/Desktop/Results, PRJNA277424_1)
#'

exploratory_analysis <- function(dds, file_name) {
  
  # Create output directories
  path_dir_pca_out <- "01-PCA"
  path_dir_sere_out <-"02-SERE_dendrogram"
  path_dir_mvv_out <- "03-Mean_vs_variance"
  dir.create(path_dir_pca_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_dir_sere_out, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_dir_mvv_out, recursive = TRUE, showWarnings = FALSE)
  
  # Extract the colData of the DESeqDataSet object
  colData_df <- as.data.frame(colData(dds))
  rownames(colData_df) <- rownames(colData(dds))
  
  # Get the Design formula
  dds_design <- as.character(design(dds))[2]
  
  # Get the test
  test <- unique(colData_df$Test)
  
  # Check the type of analysis and if there are interactions.
  design_elements <- str_split(dds_design, " \\+ ")[[1]]
  with_interaction <- any(grep(":", design_elements) > 0)
  
  # Boolean variables with this information
  multifactor <- length(design_elements) > 1
  interaction <- ""
  
  # Check the type of analysis and if there are interactions.
  if (with_interaction) {
    interaction <- "with_interaction"
  }
  
  ## Get the comparison information
  # One factor
  if (!multifactor){
    comp <- paste0("One_Factor_", test, "_", paste0(rownames(colData_df), collapse = "_"))
    
    # Multi-factor
  } else {
    comp <- paste0("Multiple_Factor_", interaction, "_", test, "_", paste0(rownames(colData_df), collapse = "_"))
  }
  
  ############################ VST NORMALIZATION ###############################

  # Absolute counts normalization for mean vs variance plot (blind = FALSE)
  deseqds <- suppressMessages(DESeq2::estimateSizeFactors(dds))
  assay(deseqds, 'counts.norm.VST.false') <- as.data.frame(assay(varianceStabilizingTransformation(deseqds, blind = FALSE)))
  
  ### 1.2 MEAN VS VARIANCE PLOT
  png(file = paste(path_dir_mvv_out, '/', file_name, '_meanvsvar.ea.png', sep = ''),
      width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  par(mfrow = c(1, 2))
  plot(log10(rowMeans(assay(deseqds, 'counts')) + 1),
       log10(row_var(assay(deseqds, 'counts')) + 1),
       xlab = expression('Log'[10] ~ 'Mean count'),
       ylab = expression('Log'[10] ~ 'Variance'),
       main = 'Counts')
  plot(rowMeans(assay(deseqds, 'counts.norm.VST.false')),
       row_var(assay(deseqds, 'counts.norm.VST.false')),
       xlab = 'Mean count',
       ylab = 'Variance',
       main = 'VST')
  dev.off()
  
  ############################# SERE DENDROGRAM ################################
  
  # Create a SERE dendrogram
  sere_dendrogram(dds, path_dir_sere_out)
  
  ###################### PRINCIPAL COMPONENT ANALYSIS (PCA) ####################

  # Recalculate vst, this time with blind = TRUE for a fully unsupervised calculation
  assay(deseqds, 'counts.norm.VST.true') <- as.data.frame(assay(varianceStabilizingTransformation(deseqds, blind = TRUE)))
  
  # Perform the PCA
  pca_res <- prcomp(x=t(assay(deseqds,'counts.norm.VST.true')), rank. = 6)
  pca_df <- as.data.frame(pca_res$x)
  
  # Get the design
  design_elements <- str_split(as.character(design(deseqds))[2], " \\+ ")[[1]]
  
  # Remove iteraction elements
  design_elements <- design_elements[!grepl(":", design_elements)]
  
  # Create the conditions column for coloring the plot
  condition_df <- as.data.frame(apply(colData_df[design_elements], 1, paste, collapse = "_"))
  colnames(condition_df) <- c("Condition_name")
  
  # Bind condition column to PCA counts dataframe (to specify sample color)
  pca_df <- merge(pca_df, condition_df, by=0, all=TRUE)
  rownames(pca_df) <- pca_df$Row.names
  pca_df$Row.names<- NULL
  
  # Create interactive plot
  plot <- plot_ly(pca_df,
                  x = pca_df[, 1],
                  y = pca_df[, 2],
                  z = pca_df[, 3],
                  type = 'scatter3d',
                  mode = 'markers',
                  text = ~rownames(pca_df),
                  color = ~pca_df$Condition_name,
                  colors = 'Paired') %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')),
           showlegend = TRUE,
           legend = list(font = list(size = 20)))
  
  # Create html file with interactive plot
  htmlwidgets::saveWidget(widget = plot,
                          file = paste(path_dir_pca_out, '/', file_name, '.ea.html', sep=''),
                          selfcontained = FALSE)
  
  ############# EUCLIDEAN AND INTRA-/INTER-GROUP DISTANCES ##################### 
  
  # Create Null matrix
  pca_comp <- pca_res$x
  distance_matrix <-  matrix(data = rep(0, nrow(pca_comp) * nrow(pca_comp)),
                             nrow = nrow(pca_comp),
                             ncol = nrow(pca_comp)
  )
  rownames(distance_matrix) <- rownames(pca_comp)
  colnames(distance_matrix) <- rownames(pca_comp)
  
  # INTRA- and INTER-group distance vectors
  intra_group <- c()
  inter_group <- c()
  
  # Complete matrix and vectors with the corresponding distances
  for (i in 1:nrow(pca_comp)) {
    for (j in 1:nrow(pca_comp)) {
      # Save distances in matrix
      distance_matrix[i,j] <- euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]);
      # Prevent the presence of repeated distances in the vectors.
      if (j > i) {
        # INTRA-group distances
        if (pca_df$Condition_name[i] == pca_df$Condition_name[j]) {
          intra_group <- c(intra_group, euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }
        # INTER-group distances
        else{
          inter_group <- c(inter_group, euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }
      }
    }
  }
  
  ######################## MANN-WHITNEY-WILCOXON TEST ########################## 

  # Perform a Mann-Whitney-Wilcoxon test
  mww_res <- wilcox.test(inter_group, intra_group, paired = FALSE)
  p_value <- mww_res[3][[1]]
  
  ############# PROPORTION OF VARIANCE EXPLAINED BY EACH COMPONENT #############

  # Create plot
  png(file = paste(path_dir_pca_out, '/', file_name, '_variance.ea.png', sep = ''),
      width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  data <- head(round(pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2) * 100, 2), n = 6)
  plot <- barplot(data, ylim = c(0, max(data) * 1.2),
                  las = 2,
                  names.arg = colnames(pca_res$x),
                  ylab = '% variance explained')
  text(x = plot, y = data, labels = data, pos = 3)
  dev.off()
  
  # Complete vector until it has 6 elements.
  # If there are no more principal components, add 0
  n_ceros <- 6 - length(data)
  ceros <- rep(0.00, n_ceros)
  data_six <- c(data, ceros)
  
  ## RESULTS
  # Return subfile name, % variance explained and MWW p-value (List)
  return(c(comp, unique(colData_df$Group), file_name, data_six, p_value, paste(rownames(colData_df), collapse = ',')))
  
}


#' Volcano plot
#' 
#' This function receives a dataframe containing the results of a DESeq2
#' differential expression analysis, which includes a column with Shrunken
#' Log2FC values, and generates a volcano plot. It also receives the alpha
#' value used in the differential expression analysis to specify which
#' sequences are differentially expressed and the path to the output file
#' where the generated plot will be saved.
#'
#' @param deseq_results  A dataframe containing the results of a DESeq2 (and
#'                       Shrunken Log2FC values).
#' @param alpha Alpha used in the differential expression analysis.
#' @param output_file Path for the output file (PNG format).
#' @return Volcano plot in PNG format.
#' @examples
#' volcano_plot(deseq_results, 0.05, /home/user/results/volcano.png)
#'

volcano_plot <- function(deseq_results, alpha, output_file) {
  
  # Create the output directory if it does not exist
  output_dir <- dirname(output_file)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Add a column with the type of sequence (up-/down-regulated or Not differentially expressed)
  deseq_results$diffexpressed <- "Not differentially expressed"
  deseq_results$diffexpressed[ deseq_results$log2FoldChange > 0 & deseq_results$padj < alpha ] <- "Up-regulated"
  deseq_results$diffexpressed[ deseq_results$log2FoldChange < 0 & deseq_results$padj < alpha ] <- "Down-regulated"
  
  # Count the number of up-regulated and down-regulated sequences
  counts <- deseq_results %>%
    as_tibble() %>%
    group_by(diffexpressed) %>%
    summarise(count = n()) %>%
    ungroup()
  
  # Create the labels for the legend
  labels <- counts %>%
    mutate(label = paste0(count, " ", diffexpressed)) %>%
    pull(label)
  names(labels) <- counts$diffexpressed
  
  # Get the maximum value in Shrunkenlog2FoldChange column
  max_abs_value <- max(abs(deseq_results$Shrunkenlog2FoldChange))
  
  # Create the volcano plot
  volcano <- ggplot(data = as.data.frame(deseq_results), aes(x = Shrunkenlog2FoldChange, y = -log10(padj), col = diffexpressed)) +
    geom_point() +
    geom_vline(xintercept = 0, col = "grey") +
    geom_hline(yintercept = -log10(alpha), col = "grey") +
    scale_color_manual(values = c("Not differentially expressed" = "snow3", "Up-regulated" = "#FF6F61", "Down-regulated" = "#6EC5E9"),
                       name = "Differential Expression",
                       labels = labels) +
    labs(x = "Shrunken Log2FoldChange",
         y = "-Log10 (Adjusted p-value)") +
    scale_x_continuous(limits = c(-max_abs_value - 0.2, max_abs_value + 0.2)) +
    guides(color = guide_legend(title = NULL)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme_bw()
  
  # Save the volcano plot
  ggsave(output_file, plot = volcano, width = 8, height = 6, dpi = 500)
  
}


#' Split Contrast
#'
#' This function extracts the elements corresponding to a specified contrast
#' from a DESeqDataSet object. It takes as input a DESeqDataSet, a string
#' representing a numerical contrast, and a vector of factors involved in
#' the contrast. The function identifies the elements of the contrast and
#' returns a vector with those elements.
#'
#' @param dds A DESeqDataSet object.
#' @param contrast A character string specifying the contrast to be
#'                 split (e.g., (50h&12h - 0h) - ((50h - 0h) + (12h - 0h))).
#' @param factors A character vector of the factors involved in the
#'                contrast (e.g., c("Time", etc)).
#' 
#' @return A vector containing the elements of the specified contrast.
#' 
#' @examples
#' contrast_elements <- split_contrast(dds,
#'   "(50h&12h - 0h) - ((50h - 0h) + (12h - 0h))",
#'   "Time"
#'   )
#' 
split_contrast <- function(dds, contrast, factors){
  
  # Remove the spaces from the contrast
  contrast_for_no_spaces <- gsub(" ", "", contrast)
  
  # Iterate through the factor columns
  factor_levels <- c()
  for (factor in factors) {
    factor_levels <- c(factor_levels, as.character(unique(colData(dds)[[factor]])))
  }
  
  # Iterate through the factor levels
  final_positions <- c()
  for (element_to_search in factor_levels){
    
    # For the pattern to be valid, it must be surrounded by _, parentheses or operators.
    pattern <- paste0("([_\\(\\+\\-])", element_to_search, "([_\\)\\+\\-])")
    
    # Obtain the match positions
    matches <- gregexpr(pattern, contrast_for_no_spaces)
    
    # If match
    if (unlist(matches)[1] != -1) {
      
      # Add a position (the specified position is that of the parenthesis or the operator).
      positions <- unlist(matches) + 1
      # Create an empty vector to store new positions
      new_positions <- c()
      
      # If the nchar of the level is 1
      if (nchar(element_to_search) == 1) {
        
        new_positions <- positions
        
      } else {
        
        # Loop through the original positions
        for (pos in positions) {
          # Obtain the positions in the formula that are part of the element being
          # searched for.
          new_numbers <- seq(pos + 1, pos + nchar(element_to_search) -1 )
          new_positions <- c(new_positions, new_numbers)  # Append to the new positions vector
        }
        new_positions <- c(positions,new_positions)
      }
      
      # Combine the original and new positions
      final_positions <- sort(c(final_positions, new_positions))
    }
  }
  
  # Split the formula
  splitted_formula <- c()
  level_name <- ""
  for (pos in 1:nchar(contrast_for_no_spaces)){
    
    # Save the parenthesis and operators in the splited formula vector
    if (!pos %in% final_positions){
      
      if (nchar(level_name) != 0){
        splitted_formula <- c(splitted_formula, level_name)
        level_name <- ""
      }
      
      # Add the operator/parenthesis
      splitted_formula <- c(splitted_formula, substr(contrast_for_no_spaces, pos, pos))
      
    } else {
      #splitted_formula <- c(splitted_formula, substr(contrast_for_no_spaces, pos, 1))
      level_name <- paste0(level_name, substr(contrast_for_no_spaces, pos, pos))
    }
  }
  
  # Initialize a vector to store the result
  result <- c()
  
  # If there are multiple factors
  if (any(grepl("_", splitted_formula))){
    
    # If there are multiple factors (the names of the conditions will appear as
    # the combination of the two factors, e.g., condition_12h).
    
    # Iterate over the original vector
    i <- 1
    while (i <= length(splitted_formula)) {
      # Check if the current element is "_"
      if (splitted_formula[i] == "_") {
        # Join the elements before and after the "_"
        result[length(result)] <- paste0(result[length(result)], "_", splitted_formula[i + 1])
        i <- i + 2  # Skip the next element since it has been joined
      } else {
        # Add the element to the result
        result <- c(result, splitted_formula[i])
        i <- i + 1
      }
    }
  } else {
    result <- splitted_formula
  }

  return(result)
}


#' Custom Contrast Differential Expression Analysis
#'
#' This function extracts the results of a differential expression analysis
#' performed with a custom contrast in DESeq2. The DESeq2 analysis must be
#' performed outside this function (e.g., using `DESeq()`), and this function
#' simply retrieves the results for a specified contrast, applying the given
#' significance threshold (alpha).
#'
#' @param dds A DESeqDataSet object containing the count data and associated
#'            metadata.
#' @param contrast A character string specifying the custom contrast to be
#'                 analyzed (e.g., "(50h&12h - 0h) - ((50h - 0h) + (12h - 0h))").
#' @param alpha A numeric value specifying the significance threshold
#'              for the results (commonly 0.05).
#' @param test A character string specifying the type of test to
#'             perform (e.g., "Wald" or "LRT").
#' 
#' @return A list containing the results of the differential expression analysis.
#' 
#' @examples
#' results_list <- custom_contrast_DEA(dds,
#'   "(50h&12h - 0h) - ((50h - 0h) + (12h - 0h))",
#'   alpha = 0.05,
#'   test = "Wald"
#' )
#' 

custom_contrast_DEA <- function(dds, contrast, alpha, test_deseq) {
  
  # Split the design formula
  design_formula_elements <- str_split(as.character(design(dds))[2], " \\+ ")[[1]]
  
  # Remove interaction terms
  factors_v <- design_formula_elements[!grepl(":", design_formula_elements)]
  
  # Remove the value associated with the intercept of the equation if it exists (0 or 1)
  if (factors_v[1] %in% c(0, 1)) factors_v <- factors_v[-1]
  
  # Create the model matrix
  mod_mat <- model.matrix(design(dds), colData(dds))
  
  # Iterate through the rows
  coldata <- as.data.frame(colData(dds))
  factor_combination <- c()
  for (i in 1:nrow(coldata)){
    
    # Select the factor columns
    row <- coldata[i, factors_v]
    
    # Create a name combining the information of all the factors
    row_string <- paste(unlist(row), collapse = "_")
    
    # Save it
    if (!(row_string %in% factor_combination)) {
      factor_combination <- c(factor_combination, row_string)
    }
  }
  
  # Calculate the coefficients
  coefficients_list <- list()
  for (conditon in factor_combination){
    # Create the string to define the coefficient vectors for each condition
    splitted_condition <- str_split(conditon,"_")[[1]]
    
    # Iterate through condition levels
    for (i in 1:length(splitted_condition)){
      # Get the level and the factor
      level <- splitted_condition[i]
      factor <- factors_v[i]
      
      # Add the condition to the string
      if (i == 1){
        string_to_select_samples <- paste0('dds$', factor, ' == "', level, '"')
      } else {
        string_to_select_samples <- paste0(string_to_select_samples, ' & dds$', factor, ' == "', level, '"')
      }
    }
    
    # Save the condition coefficients in condition list
    coefficients_list[[conditon]] <- colMeans(mod_mat[eval(parse(text=string_to_select_samples)), ])
  }
  
  # Split the contrast
  splitted_contrast <- split_contrast(dds, contrast, factors_v)
  
  # Create the contrast using the previously calculated coefficients
  # Iterate through the elements of the contrast
  final_contrast <- ""
  for (element in splitted_contrast){
    
    # Split the element (vector with their elements)
    splitted_element <- str_split(element, "_")[[1]]
    
    # Iterate through the coefficients
    for (coefficient in names(coefficients_list)){
      
      # Split the coefficient name
      splitted_coefficient <- str_split(coefficient, "_")[[1]]
      
      # Check if coefficient and element are the same thing
      if (setequal(splitted_element, splitted_coefficient)) {
        coefficient_name <- paste0('coefficients_list$`', coefficient, '`')
        condition <- TRUE
        break
      } else {
        condition <- FALSE
      }
    }
    # Check if element is a condition
    if (condition){
      final_contrast <- paste0(final_contrast, coefficient_name)
    } else {
      final_contrast <- paste0(final_contrast, element)
    }
  }
  
  # Execute the contrast
  res_deseq <- results(dds, contrast = eval(parse(text=final_contrast)), alpha = alpha, test = test_deseq)
  
  # Create the output list
  output_list <- list(
    DESeq_results = res_deseq,
    Coefficients = coefficients_list,
    Contrast = final_contrast
  )
  
  # Return the results
  return(output_list)
  
}

#' Get DESeq results
#'
#' This function extracts the results of a differential expression analysis
#' that has already been performed using DESeq2. The function retrieves the
#' results based on a specified contrast or coefficient, applying a significance
#' threshold (`alpha`) and the test type (Wald or LRT). The function returns a
#' dataframe, generates files with the results, and returns a dataframe with a
#' summary of them.
#'
#' @param dds A DESeqDataSet object containing the count data and associated
#'            metadata.
#' @param alpha A numeric value specifying the significance threshold for the
#'              results (commonly 0.05).
#' @param test A character string specifying the type of test to perform
#'             ("Wald" or "LRT").
#' @param summary_df A dataframe that accumulates the summary of differential
#'                   expression analysis results. It will be updated with
#'                   the new results.
#' @param id A string used to name the output files generated by the function.
#' @param contrast_dres A character string specifying a
#'                      custom contrast to be used for the analysis
#'                      (e.g., "(50h&12h - 0h) - ((50h - 0h) + (12h - 0h))").
#'                      Optional.
#'                      If provided, the function uses this custom contrast for the analysis.
#' @param coefficient_dres A character string or integer specifying the
#'                      coefficient to be used in the results argument of
#'                      DESeq2 (i.e., for `name` in `results()`).
#' 
#' @return A summary dataframe of the differential expression analysis results.
#' 
#' @examples
#' results_df <- get_DESeq_results(dds, alpha = 0.05, test = "Wald", 
#'                                 summary_df = summary_df, id = "experiment1",
#'                                 contrast_dres = "(50h&12h - 0h) - ((50h - 0h) + (12h - 0h))"
#'                                 )
#'

get_DESeq_results <- function(dds, alpha, test, summary_df, output_file_id,
                              contrast_dres=NULL, coefficient_dres=NULL) {
  
  # If the selected test is LRT...
  if (toupper(test) == "LRT"){
    
    # Get the results of the LRT test. The Log2FC has no value when performing a LRT!!!!
    deseq_results <- results(dds, alpha = alpha)
    
    # Significant sequences
    significant_table_lrt <- deseq_results[which(deseq_results$padj < alpha), ]
    
    # Check if any of the time columns are in the design.
    if ("Time" %in% colnames(colData(dds))) {
      time_column <- "Time"
    } else if ("Stage" %in% colnames(colData(dds))) {
      time_column <- "Stage"
    } else {
      time_column <- NA
    }

    # Rownames to column
    raw_table_lrt <- rownames_to_column(as.data.frame(deseq_results), var = "seq")
    significant_table_lrt <- rownames_to_column(as.data.frame(significant_table_lrt), var = "seq")
    
    # Save the results from LRT (raw and sig)
    write.table(as.data.frame(raw_table_lrt), 
                file = paste(output_file_id, '.general_dea_raw.tsv', sep = ''), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE)
    write.table(as.data.frame(significant_table_lrt), 
                file = paste(output_file_id, '.general_dea_sig.tsv', sep = ''), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE)
    
    
    # Add data to summary dataframe
    summary_df <- rbind(summary_df, c(output_file_id, "LRT_general", nrow(significant_table_lrt), nrow(deseq_results), "No coefficient", "No contrast", "No contrast", paste(rownames(colData(dds)), collapse = ',')))
    
    # Run the clustering of differentially expressed sequences (if time or Stage
    # are in the design)
    if (!is.na(time_column)){
      
      # Obtain the column that will be used to separate the samples in the clustering
      # The first column of the design formula
      sep_samples_col <- str_split(as.character(design(dds)), " \\+ ")[[2]][1]
      
      # Execute the clustering
      sRNA_cluster_profile(dds, deseq_results, alpha, time_column, sep_samples_col)
    }
    
  # If it is other type of analysis. Use the Wald test.
  } else {
    
    # If there is custom numerical contrast...
    if (!is.null(contrast_dres)) {
      
      # Execute the custom contrast
      custom_contrast_res <- custom_contrast_DEA(dds, contrast, alpha, "Wald")
      
      # Get the DESeq results
      deseq_results <- custom_contrast_res$DESeq_results
      
      # Execute the lfcShrink function
      resLFC <- suppressMessages(lfcShrink(dds, contrast = custom_contrast_res$Contrast, res = custom_contrast_res$DESeq_results, type = 'ashr'))
      
      # String for the summary
      custom_contrast_final  <- custom_contrast_res$Contrast
      coefficient_str <- "No coefficient"
      
      # If there is a coefficient (from result)...
    } else {
      
      # Get results. It is necessary to specify the use of the Wald test
      # (necessary to see individual comparisons when the LRT test has been
      # previously used).
      deseq_results <- results(dds, name = coefficient_dres, alpha = alpha, test = "Wald")
      
      # Execute the lfcShrink function
      resLFC <- suppressMessages(lfcShrink(dds, coef = coefficient_dres, res = deseq_results))
      
      # Strings for the sum dataframe
      coefficient_str <- coefficient_dres
      contrast <- "No contrast"
      custom_contrast_final <- "No contrast"
    }
    
    # Add Shrunken LFC to results
    final_res <- na.omit(cbind(deseq_results,
                               Shrunkenlog2FoldChange = resLFC$log2FoldChange,
                               ShrunkenlfcSE = resLFC$lfcSE))
    
    # Create a volcano plot
    volcano_plot(final_res, alpha, paste0(output_file_id, ".volcano.png"))
    
    # Extract significant differentially expressed miRNAs
    final_res_sig <- final_res %>%
      data.frame() %>%
      rownames_to_column(var = 'seq') %>%
      as_tibble() %>%
      dplyr::filter(padj <= alpha)
    
    ### RESULTS TABLES
    # Save raw and sig results
    final_res <- rownames_to_column(as.data.frame(final_res), var = 'seq')
    write.table(final_res, 
                file = paste(output_file_id, '.dea_raw.tsv', sep = ''), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE)
    write.table(as.data.frame(final_res_sig), 
                file = paste(output_file_id, '.dea_sig.tsv', sep = ''), 
                sep = "\t", 
                quote = FALSE, 
                row.names = FALSE)
    
    ### SUMMARY TABLE ###
    # Get experiment
    file_elements <- strsplit(output_file_id, split = '_')
    exp <- paste(file_elements[[1]][2], file_elements[[1]][3], sep = '.')
    
    # Significant miRNAs (p.adj < 0.05)
    final_res_sig <- as.data.frame(final_res_sig)
    sig <- nrow(final_res_sig)
    
    # Add data to dataframe
    summary_df <- rbind(summary_df, c(output_file_id, "Wald", sig, nrow(final_res), coefficient_str, contrast, custom_contrast_final,  paste(rownames(colData(dds)), collapse = ',')))
  }
  
  return(summary_df)
}


#' sRNA cluster profile
#' The function selects the sequences with an adjusted p-value lower than the
#' alpha value provided as an argument to the function and generates clusters of
#' these sequences based on their expression profile using the degPatterns
#' function from the DEGreport package. Subsequently, the expression profiles
#' of clusters with more than 10 sequences are graphically represented. The
#' function also generates a TXT file for each cluster with the associated
#' sequences.
#'
#' @param dds DeseqDataSet object
#' @param deseq_results Results from running the DESeq function using the
#'                      Likelihood Ratio Test (LRT).
#' @param alpha Adjusted p-value (padj) threshold.
#' @param time_column Character column of the DeseqDataSet object's ColData
#'                    that will be used as a variable that changes (normally a
#'                    time variable).
#' @param condition_column Character column of the DeseqDataSet object's ColData
#'                         that will be used to separate samples (normally
#'                         control/treated or control/mutant).
#' @param output_dir Path for the output directory.
#' @return No value is returned.
#' @examples
#' sRNA_cluster_profile(dds, deseq_results, 0.05, "Stage", "Condition")
#'

sRNA_cluster_profile <- function(dds, deseq_results, alpha, time_column, condition_column) {
  
  # Create a tibble for LRT results
  res_LRT_tb <- deseq_results %>%
    data.frame() %>%
    rownames_to_column(var="sRNA") %>% 
    as_tibble()
  
  # Subset to return sRNA with padj < alpha
  sigLRT_sRNAs <- res_LRT_tb %>% 
    filter(padj < alpha)
  
  # Normalize the count matrix using the rlog function.
  rld_mat <- rlog(dds)
  
  # Filter the counts matrix selecting the significant sRNAs
  rld_mat_filt <- rld_mat[sigLRT_sRNAs$sRNA, ]
  
  # Create clusters of sRNA based on their expression profile. 
  clusters <- suppressMessages(degPatterns(assay(rld_mat_filt),
                          colData(dds),
                          minc = 10,
                          pattern = NULL,
                          time = time_column,
                          col = condition_column,
                          plot = FALSE))
  
  ################### Create the title of the clusters #########################
  
  # Get the summary table
  summary_table <- clusters$summarise
  
  # Get the clusters
  clusters_v <- unique(summary_table$cluster)
  
  # Iterate through the clusters
  plots_titles <- character(length(clusters_v))
  clusters_names <- character(length(clusters_v))
  for (i in 1:length(clusters_v)) {
    # Create the new cluster name
    cluster_name_new <- paste0("Cluster ", i)
    
    # Get the number of sequences of the cluster
    num_seqs <- unique(summary_table[summary_table$cluster == clusters_v[i],]$n_genes)
    
    # Store the title in the vector
    plots_titles[i] <- paste0(cluster_name_new, ". Nº seqs: ", num_seqs)
    clusters_names[i] <- cluster_name_new
    
    # Assign the name in the same iteration
    names(plots_titles)[i] <- clusters_v[i]
    names(clusters_names)[i] <- clusters_v[i]
  }
  
  ################## Specify the elements of the legend ########################
  
  # Get the colors
  number_of_levels <- length(unique(colData(dds)[[condition_column]]))
  colors <- as.character(paletteer_d("ggthemes::Tableau_10")[1:number_of_levels])
  
  # Create the default plot
  default_plot <- degPlotCluster(clusters[["normalized"]], time_column, condition_column)
  
  # Custom plot (Standard error)
  plot_se <- ggplot(clusters[["normalized"]],
                    aes(!!sym(time_column),
                        value,
                        color = !!sym(condition_column),
                        fill = !!sym(condition_column),
                        group = !!sym(condition_column))) +
    theme(legend.title=element_blank()) +
    geom_vline(xintercept = as.factor(18),
               linetype = 3) +
    stat_summary(fun = mean,
                 geom = "line") + 
    stat_summary(fun = mean,
                 geom = "point") +
    stat_summary(fun.data = mean_se,
                 geom = "errorbar",
                 width = 0.2,
                 alpha = 0.8) +
    facet_wrap(~ cluster,
               labeller = labeller(cluster = plots_titles)) +
    ylab("Z-score") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_few() +
    guides(color = guide_legend(title = NULL),
           fill = guide_legend(title = NULL)) +
    theme(
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20))
    )
  
  # Custom plot (Standard deviation)
  plot_sdl <- ggplot(clusters[["normalized"]],
                     aes(!!sym(time_column),
                         value,
                         color = !!sym(condition_column),
                         fill = !!sym(condition_column),
                         group = !!sym(condition_column))) +
    theme(legend.title=element_blank()) +
    geom_vline(xintercept = as.factor(18),
               linetype = 3) +
    stat_summary(fun = mean,
                 geom = "line") + 
    stat_summary(fun = mean,
                 geom = "point") +
    stat_summary(fun.data = mean_sdl,
                 geom = "errorbar",
                 width = 0.2,
                 alpha = 0.8) +
    facet_wrap(~ cluster,
               labeller = labeller(cluster = plots_titles)) +
    ylab("Z-score") +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    theme_few() +
    guides(color = guide_legend(title = NULL),
           fill = guide_legend(title = NULL)) +
    theme(
      axis.title.x = element_text(margin = margin(t = 20)),
      axis.title.y = element_text(margin = margin(r = 20))
    )
  
  ########################## Save the output plots #############################
  
  # Save the plots
  ggsave(paste0(path_plots, "/default.cluster.png"), plot = default_plot, width = 8, height = 6, dpi = 300)
  ggsave(paste0(path_plots, "/standard_error.cluster.png"), plot = plot_se, width = 8, height = 6, dpi = 300)
  ggsave(paste0(path_plots, "/standard_deviation.cluster.png"), plot = plot_sdl, width = 8, height = 6, dpi = 300)

  ########################## Save the output tables ############################

  # Get the sequences belonging to each cluster
  for (cluster_id in names(clusters_names)) {
    
    # Get the sequences of the cluster
    sequences_of_cluster_df <- clusters$df[clusters$df$cluster == cluster_id,]
    
    # Save the sequences in a TXT file
    writeLines(sequences_of_cluster_df$genes, paste0(gsub(" ", "_", clusters_names[cluster_id]), "_sequences.cluster.txt"))
  }
}


################################################################################
################################### MAIN #######################################
################################################################################

# Get programm arguments
args <- get_arguments()

# Save the rest of the arguments in variables
subproject <- args$id
file <- args$counts
metadata <- args$metadata
alpha <- args$alpha
min_counts <- args$min_counts
min_samples <- args$min_samples

######################### Create the DeseqDataSet ##############################

# Filter the counts matrix by low counts and create the DeseqDataSet
dds <- create_DeseqDataSet(subproject, file, metadata, min_counts, min_samples)

########################## Exploratory analysis ################################

# Perform an exploratory analysis.
ea_results <- exploratory_analysis(dds, subproject)

# Get the test
test <- unique(colData(dds)$Test)

cat("\n######################## ", subproject, " ########################\n\n")
cat("- Test: ", test, "\n")

###################### Differential expression analysis ########################
sum <- data.frame()

# If the selected test is LRT...
if (toupper(test) == "LRT"){
  
  # Get the reduced formula from the colData
  reduced_formula <- unique(colData(dds)$DesignRed)
  
  # Execute the DEA (LRT)
  dds <- suppressMessages(DESeq(dds, test = "LRT", reduced = as.formula(reduced_formula)))
  
  # Get the DESeq results
  sum <- get_DESeq_results(dds, alpha, test, sum, paste0(subproject,"_0"))
  
  # Print some information
  cat("- Full model: ", unique(colData(dds)$Design), "\n")
  cat("- Reduced model: ", reduced_formula,  "\n")
  
# If it is other type of analysis. Use the Wald test.
} else {
  
  # Execute the DEA (Wald test)
  dds <- suppressMessages(DESeq(dds))
  
  # Print some information
  cat("- Design: ", unique(colData(dds)$Design), "\n")
}

# Only for designs with more than one comparison
subfile_num <- 1

# If there are custom numerical contrasts...
if ("Contrast" %in% colnames(colData(dds))){
  
  # Print some information
  cat("- Contrasts:\n")
  
  # Obtain the set of custom contrasts
  group_of_contrast_str <- unique(colData(dds)$Contrast)
  
  # Get a vector with the contrasts
  group_of_contrast_v <- str_split(group_of_contrast_str, ":")[[1]]
  
  # Iterate through the contrasts
  for (contrast in group_of_contrast_v) {
    
    # If there is more than one comparison
    if (length(group_of_contrast_v) > 1) {
      output_file_id <- paste0(subproject,"_", subfile_num)
      
      # If the formula uses only one factor ...
    } else {
      output_file_id <- paste0(subproject,"_0")
    }
    
    # Print some information
    cat("\t ", output_file_id, " -> ", contrast, "\n")
    
    # Get the DESeq results
    sum <- get_DESeq_results(dds, alpha, "Wald", sum, output_file_id, contrast_dres=contrast)
    
    # Increment the subfile_num variable
    subfile_num <- subfile_num + 1
    
  }

# If no contrast is specified, use the coefficients returned by resultsNames()...
} else {
  
  # Print some information
  cat("- Coefficients:\n")
  
  # Iterate through comparisons
  results_names <- resultsNames(dds)
  for (comp in results_names){
    if (comp != 'Intercept') {
      
      # If there is more than one comparison
      if (length(results_names) > 2) {
        output_file_id <- paste0(subproject,"_", subfile_num)
        
        # If the formula uses only one factor ...
      } else {
        output_file_id <- paste0(subproject,"_0")
      }
      
      # Print some information
      cat("\t ", output_file_id, " -> ", comp, "\n")
      
      # Get the DESeq results
      sum <- get_DESeq_results(dds, alpha, "Wald", sum, output_file_id, coefficient_dres=comp)
      
      # Increment the subfile_num variable
      subfile_num <- subfile_num + 1
      
    }
  }
}

# end files loop
cat('\n')

# Save Differential expression analysis summary file
colnames(sum) <- c('Group', 'Test', 'Padj<0.05', 'Total', 'Coefficient', 'Contrast', 'Contrast_coefficient', 'Samples')
write.table(sum,
            file=paste0(subproject, '.dea_summary.tsv', sep=""),
            quote=FALSE,
            sep='\t',
            row.names = FALSE)


# Save the proportion of variance explained by each component and the p-value obtained in the MWW test in a .csv file
# Save the exploratory analysis results
ea_df <- t(as.data.frame(ea_results[2:10]))
colnames(ea_df) <- c('Group_id', 'Group', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'P-value(MWW)')
write.table(ea_df,
            file = paste0(subproject, '.ea_summary.tsv', sep=""),
            quote=FALSE,
            sep='\t',
            row.names = FALSE)


